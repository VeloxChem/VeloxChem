//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include <hip/hip_runtime.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DftFunc.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "GtoValuesGgaGPU.hpp"
#include "GtoValuesLdaGPU.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MpiFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"
#include "XCIntegratorGPU.hpp"

#define TILE_DIM 16

#define hipSafe(e)                                                                                                       \
    {                                                                                                                     \
        hipError_t err = (e);                                                                                            \
        if (err != hipSuccess)                                                                                           \
        {                                                                                                                 \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                      \
        }                                                                                                                 \
    }

namespace gpu {  // gpu namespace

__global__ void
getSubDensityMatrix(double* d_den_mat, const double* d_den_mat_full, const uint32_t naos, const uint32_t* d_ao_inds, const uint32_t aocount)
{
    const uint32_t row = blockDim.y * blockIdx.y + threadIdx.y;
    const uint32_t col = blockDim.x * blockIdx.x + threadIdx.x;

    if ((row < aocount) && (col < aocount))
    {
        const auto row_orig = d_ao_inds[row];
        const auto col_orig = d_ao_inds[col];

        d_den_mat[row * aocount + col] = d_den_mat_full[row_orig * naos + col_orig];
    }
}

__global__ void
zeroKohnShamMatrix(double* d_mat_Vxc_full, const uint32_t naos)
{
    const uint32_t row = blockDim.y * blockIdx.y + threadIdx.y;
    const uint32_t col = blockDim.x * blockIdx.x + threadIdx.x;

    if ((row < naos) && (col < naos))
    {
        d_mat_Vxc_full[row * naos + col] = 0.0;
    }
}

__global__ void
cudaDensityOnGrids(double* d_rho, const double* d_mat_F, const double* d_gto_values, const uint32_t aocount, const uint32_t npoints)
{
    __shared__ double As[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM][TILE_DIM + 1];

    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    double rho_a = 0.0;

    for (uint32_t k = 0; k < (aocount + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((k * TILE_DIM + threadIdx.y < aocount) && (g < npoints))
        {
            As[threadIdx.y][threadIdx.x] = d_mat_F[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs[threadIdx.y][threadIdx.x] = d_gto_values[(k * TILE_DIM + threadIdx.y) * npoints + g];
        }
        else
        {
            As[threadIdx.y][threadIdx.x] = 0.0;
            Bs[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

        if (threadIdx.y == 0)
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                rho_a += As[m][threadIdx.x] * Bs[m][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (g < npoints))
    {
        d_rho[2 * g + 0] = rho_a;
        d_rho[2 * g + 1] = rho_a;
    }
}

__global__ void
cudaDensitySigmaOnGrids(double*        d_rho,
                        double*        d_rhograd,
                        double*        d_sigma,
                        const double*  d_mat_F,
                        const double*  d_gto_values,
                        const double*  d_gto_values_x,
                        const double*  d_gto_values_y,
                        const double*  d_gto_values_z,
                        const uint32_t aocount,
                        const uint32_t npoints)
{
    __shared__ double As[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs_x[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs_y[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs_z[TILE_DIM][TILE_DIM + 1];

    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    double rho_a = 0.0;
    double rho_a_x = 0.0;
    double rho_a_y = 0.0;
    double rho_a_z = 0.0;

    for (uint32_t k = 0; k < (aocount + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((k * TILE_DIM + threadIdx.y < aocount) && (g < npoints))
        {
            As[threadIdx.y][threadIdx.x] = d_mat_F[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs[threadIdx.y][threadIdx.x] = d_gto_values[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_x[threadIdx.y][threadIdx.x] = d_gto_values_x[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_y[threadIdx.y][threadIdx.x] = d_gto_values_y[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_z[threadIdx.y][threadIdx.x] = d_gto_values_z[(k * TILE_DIM + threadIdx.y) * npoints + g];
        }
        else
        {
            As[threadIdx.y][threadIdx.x] = 0.0;
            Bs[threadIdx.y][threadIdx.x] = 0.0;
            Bs_x[threadIdx.y][threadIdx.x] = 0.0;
            Bs_y[threadIdx.y][threadIdx.x] = 0.0;
            Bs_z[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

        if (threadIdx.y == 0)
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                rho_a += As[m][threadIdx.x] * Bs[m][threadIdx.x];
                rho_a_x += 2.0 * (As[m][threadIdx.x] * Bs_x[m][threadIdx.x]);
                rho_a_y += 2.0 * (As[m][threadIdx.x] * Bs_y[m][threadIdx.x]);
                rho_a_z += 2.0 * (As[m][threadIdx.x] * Bs_z[m][threadIdx.x]);
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (g < npoints))
    {
        d_rho[2 * g + 0] = rho_a;
        d_rho[2 * g + 1] = rho_a;

        d_rhograd[6 * g + 0] = rho_a_x;
        d_rhograd[6 * g + 1] = rho_a_y;
        d_rhograd[6 * g + 2] = rho_a_z;
        d_rhograd[6 * g + 3] = rho_a_x;
        d_rhograd[6 * g + 4] = rho_a_y;
        d_rhograd[6 * g + 5] = rho_a_z;

        d_sigma[3 * g + 0] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
        d_sigma[3 * g + 1] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
        d_sigma[3 * g + 2] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
    }
}

__global__ void
cudaGetLdaMatrixG(double*        d_mat_G,
                  const double*  d_grid_w,
                  const uint32_t grid_offset,
                  const uint32_t npoints,
                  const double*  d_gto_values,
                  const uint32_t aocount,
                  const double*  d_vrho)
{
    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((i < aocount) && (g < npoints))
    {
        d_mat_G[g + i * npoints] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0] * d_gto_values[g + i * npoints];
    }
}

__global__ void
cudaGetGgaMatrixG(double*        d_mat_G,
                  const double*  d_grid_w,
                  const uint32_t grid_offset,
                  const uint32_t npoints,
                  const double*  d_gto_values,
                  const double*  d_gto_values_x,
                  const double*  d_gto_values_y,
                  const double*  d_gto_values_z,
                  const uint32_t aocount,
                  const double*  d_rhograd,
                  const double*  d_vrho,
                  const double*  d_vsigma)
{
    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    // Note: Assuming symmetric density and KohnSham matrices. Only works for ground state.

    if ((i < aocount) && (g < npoints))
    {
        const double vx = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 0] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 3];
        const double vy = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 1] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 4];
        const double vz = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 2] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 5];

        d_mat_G[g + i * npoints] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0] * d_gto_values[g + i * npoints] +
                                   2.0 * (d_grid_w[g + grid_offset] * (vx * d_gto_values_x[g + i * npoints] + vy * d_gto_values_y[g + i * npoints] +
                                                                       vz * d_gto_values_z[g + i * npoints]));
    }
}

__global__ void
matmulAtB(double* C, double* A, double* B, uint32_t aocount, uint32_t npoints)
{
    __shared__ double As[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM][TILE_DIM + 1];

    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    double val = 0.0;

    for (uint32_t k = 0; k < (aocount + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((i < aocount) && (k * TILE_DIM + threadIdx.x < aocount))
        {
            // note: transpose of A
            As[threadIdx.x][threadIdx.y] = A[(k * TILE_DIM + threadIdx.x) * aocount + i];
        }
        else
        {
            As[threadIdx.x][threadIdx.y] = 0.0;
        }

        if ((k * TILE_DIM + threadIdx.y < aocount) && (g < npoints))
        {
            Bs[threadIdx.y][threadIdx.x] = B[(k * TILE_DIM + threadIdx.y) * npoints + g];
        }
        else
        {
            Bs[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

        for (uint32_t m = 0; m < TILE_DIM; m++)
        {
            // note: transpose of A
            val += As[m][threadIdx.y] * Bs[m][threadIdx.x];
        }

        __syncthreads();
    }

    if ((i < aocount) && (g < npoints))
    {
        C[i * npoints + g] = val;
    }
}

__global__ void
matmulABtKohnShamMatrix(double*         d_mat_Vxc_full,
                        const uint32_t  naos,
                        const double*   A,
                        const double*   B,
                        const uint32_t* d_ao_inds,
                        const uint32_t  aocount,
                        const uint32_t  npoints)
{
    __shared__ double As[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM][TILE_DIM + 1];

    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t j = blockIdx.x * blockDim.x + threadIdx.x;

    double val = 0.0;

    for (uint32_t k = 0; k < (npoints + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((i < aocount) && (k * TILE_DIM + threadIdx.x < npoints))
        {
            As[threadIdx.y][threadIdx.x] = A[i * npoints + (k * TILE_DIM + threadIdx.x)];
        }
        else
        {
            As[threadIdx.y][threadIdx.x] = 0.0;
        }

        if ((k * TILE_DIM + threadIdx.y < npoints) && (j < aocount))
        {
            // note: transpose of B
            Bs[threadIdx.x][threadIdx.y] = B[j * npoints + (k * TILE_DIM + threadIdx.y)];
        }
        else
        {
            Bs[threadIdx.x][threadIdx.y] = 0.0;
        }

        __syncthreads();

        for (uint32_t m = 0; m < TILE_DIM; m++)
        {
            // note: transpose of B
            val += As[threadIdx.y][m] * Bs[threadIdx.x][m];
        }

        __syncthreads();
    }

    if ((i < aocount) && (j < aocount))
    {
        d_mat_Vxc_full[d_ao_inds[i] * naos + d_ao_inds[j]] += val;
    }
}

static auto
getGtoInfo(const CGtoBlock gto_block, const std::vector<int64_t>& gtos_mask) -> std::vector<double>
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // set up data on host and device

    std::vector<double> gto_info(5 * nrows * npgtos);

    auto gto_info_ptr = gto_info.data();

    for (int64_t i = 0, irow = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            const auto r_x = gto_coords[i][0];
            const auto r_y = gto_coords[i][1];
            const auto r_z = gto_coords[i][2];

            for (int64_t j = 0; j < npgtos; j++)
            {
                const auto fexp  = gto_exps[j * ncgtos + i];
                const auto fnorm = gto_norms[j * ncgtos + i];

                gto_info_ptr[irow + j * nrows + npgtos * nrows * 0] = fexp;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 1] = fnorm;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 2] = r_x;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 3] = r_y;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 4] = r_z;
            }

            irow++;
        }
    }

    return gto_info;
}

static auto
getGtoValuesForLda(double*                     d_gto_values,
                   const int64_t               row_offset,
                   double*                     d_gto_info,
                   const CGtoBlock&            gto_block,
                   const double*               d_grid_x,
                   const double*               d_grid_y,
                   const double*               d_grid_z,
                   const int64_t               grid_offset,
                   const int64_t               n_grid_points,
                   const std::vector<int64_t>& gtos_mask) -> void
{
    // number of useful CGTOs

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    // number of primitives per CGTO

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // number of grid points

    const auto ncols = n_grid_points;

    // prepare GTO information

    auto gto_info = gpu::getGtoInfo(gto_block, gtos_mask);

    hipSafe(hipMemcpy(d_gto_info, gto_info.data(), gto_info.size() * sizeof(double), hipMemcpyHostToDevice));

    // evaluate GTO values on grid points

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((ncols + threads_per_block.x - 1) / threads_per_block.x, (nrows + threads_per_block.y - 1) / threads_per_block.y);

    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        hipLaunchKernelGGL(gpu::gtoValuesLdaRecS, num_blocks, threads_per_block, 0, 0, d_gto_values,
                                                                 static_cast<uint32_t>(row_offset),
                                                                 d_gto_info,
                                                                 d_grid_x,
                                                                 d_grid_y,
                                                                 d_grid_z,
                                                                 static_cast<uint32_t>(grid_offset),
                                                                 static_cast<uint32_t>(nrows),
                                                                 static_cast<uint32_t>(npgtos),
                                                                 static_cast<uint32_t>(ncols));
    }
    else if (gto_ang == 1)
    {
        hipLaunchKernelGGL(gpu::gtoValuesLdaRecP, num_blocks, threads_per_block, 0, 0, d_gto_values,
                                                                 static_cast<uint32_t>(row_offset),
                                                                 d_gto_info,
                                                                 d_grid_x,
                                                                 d_grid_y,
                                                                 d_grid_z,
                                                                 static_cast<uint32_t>(grid_offset),
                                                                 static_cast<uint32_t>(nrows),
                                                                 static_cast<uint32_t>(npgtos),
                                                                 static_cast<uint32_t>(ncols));
    }
    else if (gto_ang == 2)
    {
        const double f2_3 = 2.0 * std::sqrt(3.0);

        hipLaunchKernelGGL(gpu::gtoValuesLdaRecD, num_blocks, threads_per_block, 0, 0, d_gto_values,
                                                                 static_cast<uint32_t>(row_offset),
                                                                 f2_3,
                                                                 d_gto_info,
                                                                 d_grid_x,
                                                                 d_grid_y,
                                                                 d_grid_z,
                                                                 static_cast<uint32_t>(grid_offset),
                                                                 static_cast<uint32_t>(nrows),
                                                                 static_cast<uint32_t>(npgtos),
                                                                 static_cast<uint32_t>(ncols));
    }
    else
    {
        std::string err_ang("gpu::getGtoValuesForLda: Only implemented for s, p and d-orbitals");

        errors::assertMsgCritical(false, err_ang);
    }
}

static auto
getGtoValuesForGga(double*                     d_gto_values_0,
                   double*                     d_gto_values_x,
                   double*                     d_gto_values_y,
                   double*                     d_gto_values_z,
                   const int64_t               row_offset,
                   double*                     d_gto_info,
                   const CGtoBlock&            gto_block,
                   const double*               d_grid_x,
                   const double*               d_grid_y,
                   const double*               d_grid_z,
                   const int64_t               grid_offset,
                   const int64_t               n_grid_points,
                   const std::vector<int64_t>& gtos_mask) -> void
{
    // number of useful CGTOs

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    // number of primitives per CGTO

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // number of grid points

    const auto ncols = n_grid_points;

    // prepare GTO information

    auto gto_info = gpu::getGtoInfo(gto_block, gtos_mask);

    hipSafe(hipMemcpy(d_gto_info, gto_info.data(), gto_info.size() * sizeof(double), hipMemcpyHostToDevice));

    // evaluate GTO values on grid points

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((ncols + threads_per_block.x - 1) / threads_per_block.x, (nrows + threads_per_block.y - 1) / threads_per_block.y);

    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        hipLaunchKernelGGL(gpu::gtoValuesGgaRecS, num_blocks, threads_per_block, 0, 0, d_gto_values_0,
                                                                 d_gto_values_x,
                                                                 d_gto_values_y,
                                                                 d_gto_values_z,
                                                                 static_cast<uint32_t>(row_offset),
                                                                 d_gto_info,
                                                                 d_grid_x,
                                                                 d_grid_y,
                                                                 d_grid_z,
                                                                 static_cast<uint32_t>(grid_offset),
                                                                 static_cast<uint32_t>(nrows),
                                                                 static_cast<uint32_t>(npgtos),
                                                                 static_cast<uint32_t>(ncols));
    }
    else if (gto_ang == 1)
    {
        hipLaunchKernelGGL(gpu::gtoValuesGgaRecP, num_blocks, threads_per_block, 0, 0, d_gto_values_0,
                                                                 d_gto_values_x,
                                                                 d_gto_values_y,
                                                                 d_gto_values_z,
                                                                 static_cast<uint32_t>(row_offset),
                                                                 d_gto_info,
                                                                 d_grid_x,
                                                                 d_grid_y,
                                                                 d_grid_z,
                                                                 static_cast<uint32_t>(grid_offset),
                                                                 static_cast<uint32_t>(nrows),
                                                                 static_cast<uint32_t>(npgtos),
                                                                 static_cast<uint32_t>(ncols));
    }
    else if (gto_ang == 2)
    {
        const double f2_3 = 2.0 * std::sqrt(3.0);

        hipLaunchKernelGGL(gpu::gtoValuesGgaRecD, num_blocks, threads_per_block, 0, 0, d_gto_values_0,
                                                                 d_gto_values_x,
                                                                 d_gto_values_y,
                                                                 d_gto_values_z,
                                                                 static_cast<uint32_t>(row_offset),
                                                                 f2_3,
                                                                 d_gto_info,
                                                                 d_grid_x,
                                                                 d_grid_y,
                                                                 d_grid_z,
                                                                 static_cast<uint32_t>(grid_offset),
                                                                 static_cast<uint32_t>(nrows),
                                                                 static_cast<uint32_t>(npgtos),
                                                                 static_cast<uint32_t>(ncols));
    }
    else
    {
        std::string err_ang("gpu::getGtoValuesForGga: Only implemented for s, p and d-orbitals");

        errors::assertMsgCritical(false, err_ang);
    }
}

auto
computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> CDenseMatrix
{
    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    double* d_gto_info;

    hipSafe(hipMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double* d_gaos;

    hipSafe(hipMalloc(&d_gaos, naos * max_npoints_per_box * sizeof(double)));

    // coordinates of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z;

    hipSafe(hipMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));

    hipSafe(hipMemcpy(d_grid_x, xcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_y, ycoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_z, zcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening

        std::vector<std::vector<int64_t>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int64_t> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 0, 1.0e-12, boxdim);

            auto cgto_mask   = std::get<0>(pre_scr_info);
            auto pre_ao_inds = std::get<1>(pre_scr_info);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int64_t>(aoinds.size());

        // GTO values on grid points

        CDenseMatrix mat_chi(aocount, npoints);

        const auto grid_x_ptr = xcoords + gridblockpos;
        const auto grid_y_ptr = ycoords + gridblockpos;
        const auto grid_z_ptr = zcoords + gridblockpos;

        std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);
        std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);
        std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);

        // go through GTO blocks

        int64_t row_offset = 0;

        for (size_t i_block = 0; i_block < gto_blocks.size(); i_block++)
        {
            const auto& gto_block = gto_blocks[i_block];

            const auto& cgto_mask = cgto_mask_blocks[i_block];

            const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

            gpu::getGtoValuesForLda(d_gaos, row_offset, d_gto_info, gto_block, d_grid_x, d_grid_y, d_grid_z, gridblockpos, npoints, cgto_mask);

            row_offset += static_cast<int64_t>(pre_ao_inds.size());
        }

        hipSafe(hipMemcpy(mat_chi.values(), d_gaos, aocount * npoints * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t nu = 0; nu < aocount; nu++)
        {
            std::memcpy(allgtovalues.row(aoinds[nu]) + gridblockpos, mat_chi.row(nu), npoints * sizeof(double));
        }
    }

    hipSafe(hipFree(d_gto_info));
    hipSafe(hipFree(d_gaos));
    hipSafe(hipFree(d_grid_x));
    hipSafe(hipFree(d_grid_y));
    hipSafe(hipFree(d_grid_z));

    return allgtovalues;
}

auto
computeGtoValuesAndDerivativesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid)
    -> std::vector<CDenseMatrix>
{
    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    double* d_gto_info;

    hipSafe(hipMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // GTO values on grid points

    CDenseMatrix allgtovalues_0(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_x(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_y(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_z(naos, molecularGrid.getNumberOfGridPoints());

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double *d_gaos, *d_gaox, *d_gaoy, *d_gaoz;

    hipSafe(hipMalloc(&d_gaos, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_gaox, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_gaoy, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_gaoz, naos * max_npoints_per_box * sizeof(double)));

    // coordinates of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z;

    hipSafe(hipMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));

    hipSafe(hipMemcpy(d_grid_x, xcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_y, ycoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_z, zcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening

        std::vector<std::vector<int64_t>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int64_t> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 0, 1.0e-12, boxdim);

            auto cgto_mask   = std::get<0>(pre_scr_info);
            auto pre_ao_inds = std::get<1>(pre_scr_info);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int64_t>(aoinds.size());

        // GTO values on grid points

        CDenseMatrix mat_chi_0(aocount, npoints);
        CDenseMatrix mat_chi_x(aocount, npoints);
        CDenseMatrix mat_chi_y(aocount, npoints);
        CDenseMatrix mat_chi_z(aocount, npoints);

        const auto grid_x_ptr = xcoords + gridblockpos;
        const auto grid_y_ptr = ycoords + gridblockpos;
        const auto grid_z_ptr = zcoords + gridblockpos;

        std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);
        std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);
        std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);

        // go through GTO blocks

        int64_t row_offset = 0;

        for (size_t i_block = 0; i_block < gto_blocks.size(); i_block++)
        {
            const auto& gto_block = gto_blocks[i_block];

            const auto& cgto_mask = cgto_mask_blocks[i_block];

            const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

            gpu::getGtoValuesForGga(
                d_gaos, d_gaox, d_gaoy, d_gaoz, row_offset, d_gto_info, gto_block, d_grid_x, d_grid_y, d_grid_z, gridblockpos, npoints, cgto_mask);

            row_offset += static_cast<int64_t>(pre_ao_inds.size());
        }

        hipSafe(hipMemcpy(mat_chi_0.values(), d_gaos, aocount * npoints * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_chi_x.values(), d_gaox, aocount * npoints * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_chi_y.values(), d_gaoy, aocount * npoints * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_chi_z.values(), d_gaoz, aocount * npoints * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t nu = 0; nu < aocount; nu++)
        {
            std::memcpy(allgtovalues_0.row(aoinds[nu]) + gridblockpos, mat_chi_0.row(nu), npoints * sizeof(double));
            std::memcpy(allgtovalues_x.row(aoinds[nu]) + gridblockpos, mat_chi_x.row(nu), npoints * sizeof(double));
            std::memcpy(allgtovalues_y.row(aoinds[nu]) + gridblockpos, mat_chi_y.row(nu), npoints * sizeof(double));
            std::memcpy(allgtovalues_z.row(aoinds[nu]) + gridblockpos, mat_chi_z.row(nu), npoints * sizeof(double));
        }
    }

    hipSafe(hipFree(d_gto_info));
    hipSafe(hipFree(d_gaos));
    hipSafe(hipFree(d_grid_x));
    hipSafe(hipFree(d_grid_y));
    hipSafe(hipFree(d_grid_z));

    return std::vector<CDenseMatrix>({allgtovalues_0, allgtovalues_x, allgtovalues_y, allgtovalues_z});
}

static auto
integrateVxcFockForLDA(const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CAODensityMatrix& densityMatrix,
                       const CMolecularGrid&   molecularGrid,
                       const CXCFunctional&    xcFunctional,
                       const std::string&      flag) -> CAOKohnShamMatrix
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    if (rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (rank % 8 == 7) hipSafe(hipSetDevice(1));

    CMultiTimer timer;

    timer.start("Total timing");

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("gpu::integrateVxcFockForLDA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    double* d_gto_info;

    hipSafe(hipMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    uint32_t* d_ao_inds;

    hipSafe(hipMalloc(&d_ao_inds, naos * sizeof(uint32_t)));

    // Kohn-Sham matrix

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double *d_mat_Vxc_full, *d_den_mat_full, *d_den_mat, *d_gto_values, *d_mat_F;

    hipSafe(hipMalloc(&d_den_mat, naos * naos * sizeof(double)));
    hipSafe(hipMalloc(&d_gto_values, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_F, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_den_mat_full, naos * naos * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Vxc_full, naos * naos * sizeof(double)));

    hipSafe(hipMemcpy(d_den_mat_full, densityMatrix.alphaDensity(0), naos * naos * sizeof(double), hipMemcpyHostToDevice));

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    hipLaunchKernelGGL(gpu::zeroKohnShamMatrix, num_blocks, threads_per_block, 0, 0, d_mat_Vxc_full, static_cast<uint32_t>(naos));

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    double *d_rho, *d_exc, *d_vrho;

    hipSafe(hipMalloc(&d_rho, dim->rho * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_exc, dim->zk * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double)));

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    hipSafe(hipMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_w, n_total_grid_points * sizeof(double)));

    hipSafe(hipMemcpy(d_grid_x, xcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_y, ycoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_z, zcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_w, weights, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        std::vector<std::vector<int64_t>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int64_t> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 0, 1.0e-12, boxdim);

            auto cgto_mask   = std::get<0>(pre_scr_info);
            auto pre_ao_inds = std::get<1>(pre_scr_info);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int64_t>(aoinds.size());

        if (aocount == 0) continue;

        // GTO values on grid points

        int64_t row_offset = 0;

        for (size_t i_block = 0; i_block < gto_blocks.size(); i_block++)
        {
            const auto& gto_block = gto_blocks[i_block];

            const auto& cgto_mask = cgto_mask_blocks[i_block];

            const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

            gpu::getGtoValuesForLda(d_gto_values, row_offset, d_gto_info, gto_block, d_grid_x, d_grid_y, d_grid_z, gridblockpos, npoints, cgto_mask);

            row_offset += static_cast<int64_t>(pre_ao_inds.size());
        }

        // generate sub density matrix and density grid

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        hipSafe(hipMemcpy(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), hipMemcpyHostToDevice));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::getSubDensityMatrix, num_blocks, threads_per_block, 0, 0, 
            d_den_mat, d_den_mat_full, static_cast<uint32_t>(naos), d_ao_inds, static_cast<uint32_t>(aocount));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        // Note: we should use A * B here but since A is symmetric we can also use A^T * B
        hipLaunchKernelGGL(gpu::matmulAtB, num_blocks, threads_per_block, 0, 0, 
            d_mat_F, d_den_mat, d_gto_values, static_cast<uint32_t>(aocount), static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(gpu::cudaDensityOnGrids, num_blocks, threads_per_block, 0, 0, 
            d_rho, d_mat_F, d_gto_values, static_cast<uint32_t>(aocount), static_cast<uint32_t>(npoints));

        hipSafe(hipMemcpy(rho, d_rho, 2 * npoints * sizeof(double), hipMemcpyDeviceToHost));

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        hipSafe(hipMemcpy(d_exc, exc, dim->zk * npoints * sizeof(double), hipMemcpyHostToDevice));
        hipSafe(hipMemcpy(d_vrho, vrho, dim->vrho * npoints * sizeof(double), hipMemcpyHostToDevice));

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // reuse d_den_mat and d_mat_F as working space
        auto d_mat_G = d_mat_F;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::cudaGetLdaMatrixG, num_blocks, threads_per_block, 0, 0, d_mat_G,
                                                                  d_grid_w,
                                                                  static_cast<uint32_t>(gridblockpos),
                                                                  static_cast<uint32_t>(npoints),
                                                                  d_gto_values,
                                                                  static_cast<uint32_t>(aocount),
                                                                  d_vrho);

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::matmulABtKohnShamMatrix, num_blocks, threads_per_block, 0, 0, d_mat_Vxc_full,
                                                                        static_cast<uint32_t>(naos),
                                                                        d_gto_values,
                                                                        d_mat_G,
                                                                        d_ao_inds,
                                                                        static_cast<uint32_t>(aocount),
                                                                        static_cast<uint32_t>(npoints));

        // compute partial contribution to XC energy

        for (int64_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += weights[g + gridblockpos] * rho_total;

            xcene += weights[g + gridblockpos] * exc[g] * rho_total;
        }
    }

    hipSafe(hipMemcpy(mat_Vxc.getPointerToAlphaValues(), d_mat_Vxc_full, naos * naos * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_gto_info));
    hipSafe(hipFree(d_ao_inds));

    hipSafe(hipFree(d_den_mat));
    hipSafe(hipFree(d_den_mat_full));
    hipSafe(hipFree(d_gto_values));
    hipSafe(hipFree(d_mat_F));

    hipSafe(hipFree(d_rho));
    hipSafe(hipFree(d_exc));
    hipSafe(hipFree(d_vrho));

    hipSafe(hipFree(d_grid_x));
    hipSafe(hipFree(d_grid_y));
    hipSafe(hipFree(d_grid_z));
    hipSafe(hipFree(d_grid_w));

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

    std::cout << "\nTiming of GPU integrator\n";
    std::cout << "------------------------\n";
    std::cout << timer.getSummary() << std::endl;

    return mat_Vxc;
}

static auto
integrateVxcFockForGGA(const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CAODensityMatrix& densityMatrix,
                       const CMolecularGrid&   molecularGrid,
                       const CXCFunctional&    xcFunctional,
                       const std::string&      flag) -> CAOKohnShamMatrix
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    if (rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (rank % 8 == 7) hipSafe(hipSetDevice(1));

    CMultiTimer timer;

    timer.start("Total timing");

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("gpu::integrateVxcFockForLDA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    double* d_gto_info;

    hipSafe(hipMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    uint32_t* d_ao_inds;

    hipSafe(hipMalloc(&d_ao_inds, naos * sizeof(uint32_t)));

    // Kohn-Sham matrix

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double *d_mat_Vxc_full, *d_den_mat_full, *d_den_mat, *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z, *d_mat_F;

    hipSafe(hipMalloc(&d_den_mat, naos * naos * sizeof(double)));
    hipSafe(hipMalloc(&d_gto_values, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_F, naos * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_den_mat_full, naos * naos * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Vxc_full, naos * naos * sizeof(double)));

    hipSafe(hipMemcpy(d_den_mat_full, densityMatrix.alphaDensity(0), naos * naos * sizeof(double), hipMemcpyHostToDevice));

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    hipLaunchKernelGGL(gpu::zeroKohnShamMatrix, num_blocks, threads_per_block, 0, 0, d_mat_Vxc_full, static_cast<uint32_t>(naos));

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto exc    = exc_data.data();
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    double *d_rho, *d_rhograd, *d_sigma, *d_exc, *d_vrho, *d_vsigma;

    hipSafe(hipMalloc(&d_rho, dim->rho * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_rhograd, dim->rho * 3 * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_sigma, dim->sigma * max_npoints_per_box * sizeof(double)));

    hipSafe(hipMalloc(&d_exc, dim->zk * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double)));
    hipSafe(hipMalloc(&d_vsigma, dim->vsigma * max_npoints_per_box * sizeof(double)));

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    hipSafe(hipMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_w, n_total_grid_points * sizeof(double)));

    hipSafe(hipMemcpy(d_grid_x, xcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_y, ycoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_z, zcoords, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_w, weights, n_total_grid_points * sizeof(double), hipMemcpyHostToDevice));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        std::vector<std::vector<int64_t>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int64_t> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 1st order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 1, 1.0e-12, boxdim);

            auto cgto_mask   = std::get<0>(pre_scr_info);
            auto pre_ao_inds = std::get<1>(pre_scr_info);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int64_t>(aoinds.size());

        if (aocount == 0) continue;

        // GTO values on grid points

        int64_t row_offset = 0;

        for (size_t i_block = 0; i_block < gto_blocks.size(); i_block++)
        {
            const auto& gto_block = gto_blocks[i_block];

            const auto& cgto_mask = cgto_mask_blocks[i_block];

            const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

            gpu::getGtoValuesForGga(d_gto_values,
                                    d_gto_values_x,
                                    d_gto_values_y,
                                    d_gto_values_z,
                                    row_offset,
                                    d_gto_info,
                                    gto_block,
                                    d_grid_x,
                                    d_grid_y,
                                    d_grid_z,
                                    gridblockpos,
                                    npoints,
                                    cgto_mask);

            row_offset += static_cast<int64_t>(pre_ao_inds.size());
        }

        // generate sub density matrix and density grid

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        hipSafe(hipMemcpy(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), hipMemcpyHostToDevice));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::getSubDensityMatrix, num_blocks, threads_per_block, 0, 0, 
            d_den_mat, d_den_mat_full, static_cast<uint32_t>(naos), d_ao_inds, static_cast<uint32_t>(aocount));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        // Note: we should use A * B here but since A is symmetric we can also use A^T * B
        hipLaunchKernelGGL(gpu::matmulAtB, num_blocks, threads_per_block, 0, 0, 
            d_mat_F, d_den_mat, d_gto_values, static_cast<uint32_t>(aocount), static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(gpu::cudaDensitySigmaOnGrids, num_blocks, threads_per_block, 0, 0, d_rho,
                                                                           d_rhograd,
                                                                           d_sigma,
                                                                           d_mat_F,
                                                                           d_gto_values,
                                                                           d_gto_values_x,
                                                                           d_gto_values_y,
                                                                           d_gto_values_z,
                                                                           static_cast<uint32_t>(aocount),
                                                                           static_cast<uint32_t>(npoints));

        hipSafe(hipMemcpy(rho, d_rho, dim->rho * npoints * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(rhograd, d_rhograd, dim->rho * 3 * npoints * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(sigma, d_sigma, dim->sigma * npoints * sizeof(double), hipMemcpyDeviceToHost));

        xcFunctional.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        hipSafe(hipMemcpy(d_exc, exc, dim->zk * npoints * sizeof(double), hipMemcpyHostToDevice));
        hipSafe(hipMemcpy(d_vrho, vrho, dim->vrho * npoints * sizeof(double), hipMemcpyHostToDevice));
        hipSafe(hipMemcpy(d_vsigma, vsigma, dim->vsigma * npoints * sizeof(double), hipMemcpyHostToDevice));

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // reuse d_den_mat and d_mat_F as working space
        auto d_mat_G = d_mat_F;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::cudaGetGgaMatrixG, num_blocks, threads_per_block, 0, 0, d_mat_G,
                                                                  d_grid_w,
                                                                  static_cast<uint32_t>(gridblockpos),
                                                                  static_cast<uint32_t>(npoints),
                                                                  d_gto_values,
                                                                  d_gto_values_x,
                                                                  d_gto_values_y,
                                                                  d_gto_values_z,
                                                                  static_cast<uint32_t>(aocount),
                                                                  d_rhograd,
                                                                  d_vrho,
                                                                  d_vsigma);

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::matmulABtKohnShamMatrix, num_blocks, threads_per_block, 0, 0, d_mat_Vxc_full,
                                                                        static_cast<uint32_t>(naos),
                                                                        d_gto_values,
                                                                        d_mat_G,
                                                                        d_ao_inds,
                                                                        static_cast<uint32_t>(aocount),
                                                                        static_cast<uint32_t>(npoints));

        // compute partial contribution to XC energy

        for (int64_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += weights[g + gridblockpos] * rho_total;

            xcene += weights[g + gridblockpos] * exc[g] * rho_total;
        }
    }

    hipSafe(hipMemcpy(mat_Vxc.getPointerToAlphaValues(), d_mat_Vxc_full, naos * naos * sizeof(double), hipMemcpyDeviceToHost));

    // Note: symmetrize mat_Vxc

    mat_Vxc.symmetrizeAndScale(0.5);

    hipSafe(hipFree(d_gto_info));
    hipSafe(hipFree(d_ao_inds));

    hipSafe(hipFree(d_den_mat));
    hipSafe(hipFree(d_den_mat_full));
    hipSafe(hipFree(d_gto_values));
    hipSafe(hipFree(d_mat_F));

    hipSafe(hipFree(d_rho));
    hipSafe(hipFree(d_exc));
    hipSafe(hipFree(d_vrho));

    hipSafe(hipFree(d_grid_x));
    hipSafe(hipFree(d_grid_y));
    hipSafe(hipFree(d_grid_z));
    hipSafe(hipFree(d_grid_w));

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

    std::cout << "\nTiming of GPU integrator\n";
    std::cout << "------------------------\n";
    std::cout << timer.getSummary() << std::endl;

    return mat_Vxc;
}

auto
integrateVxcFock(const CMolecule&        molecule,
                 const CMolecularBasis&  basis,
                 const CAODensityMatrix& densityMatrix,
                 const CMolecularGrid&   molecularGrid,
                 const std::string&      xcFuncLabel) -> CAOKohnShamMatrix
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    auto flag = densityMatrix.isClosedShell() ? std::string("CLOSEDSHELL") : std::string("OPENSHELL");

    std::string erropenshell("gpu::integrateVxcFock: Only implemented for closed-shell");

    errors::assertMsgCritical(densityMatrix.isClosedShell(), erropenshell);

    if (xcfuntype == xcfun::lda) return gpu::integrateVxcFockForLDA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    if (xcfuntype == xcfun::gga) return gpu::integrateVxcFockForGGA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    /*
    if (xcfuntype == xcfun::mgga) return gpu::integrateVxcFockForMGGA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    std::string errxcfuntype("gpu::integrateVxcFock: Only implemented for LDA/GGA/meta-GGA");

    errors::assertMsgCritical(false, errxcfuntype);
    */

    return CAOKohnShamMatrix();
}

}  // namespace gpu
