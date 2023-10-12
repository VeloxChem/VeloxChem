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

#include <cublas_v2.h>
#include <cuda_runtime.h>

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
#include "GtoValues.hpp"
#include "GtoValuesRecF.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"
#include "XCIntegratorGPU.hpp"

#define cudaSafe(e)                                                                                                       \
    {                                                                                                                     \
        cudaError_t err = (e);                                                                                            \
        if (err != cudaSuccess)                                                                                           \
        {                                                                                                                 \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                      \
        }                                                                                                                 \
    }

#define cublasSafe(e)                                                                            \
    {                                                                                            \
        cublasStatus_t err = (e);                                                                \
        if (err != CUBLAS_STATUS_SUCCESS)                                                        \
        {                                                                                        \
            std::cerr << "cuBLAS error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                             \
        }                                                                                        \
    }

namespace gpu {  // gpu namespace

__global__ void
cudaLdaValuesDirectRecS(double*        gto_values,
                        const uint32_t row_offset,
                        const double*  gto_info,
                        const double*  grid_x,
                        const double*  grid_y,
                        const double*  grid_z,
                        const uint32_t grid_offset,
                        const uint32_t nrows,
                        const uint32_t npgtos,
                        const uint32_t ncols)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t g = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        double s0 = 0.0;

        const auto g_x = grid_x[g + grid_offset];
        const auto g_y = grid_y[g + grid_offset];
        const auto g_z = grid_z[g + grid_offset];

        for (uint32_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_info[i + j * nrows + npgtos * nrows * 0];
            const auto fnorm = gto_info[i + j * nrows + npgtos * nrows * 1];
            const auto r_x   = gto_info[i + j * nrows + npgtos * nrows * 2];
            const auto r_y   = gto_info[i + j * nrows + npgtos * nrows * 3];
            const auto r_z   = gto_info[i + j * nrows + npgtos * nrows * 4];

            const auto gr_x = g_x - r_x;
            const auto gr_y = g_y - r_y;
            const auto gr_z = g_z - r_z;

            s0 += fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
        }

        gto_values[g + (i + row_offset) * ncols] = s0;
    }
}

__global__ void
cudaLdaValuesDirectRecP(double*        gto_values_p3,
                        const uint32_t row_offset,
                        const double*  gto_info,
                        const double*  grid_x,
                        const double*  grid_y,
                        const double*  grid_z,
                        const uint32_t grid_offset,
                        const uint32_t nrows,
                        const uint32_t npgtos,
                        const uint32_t ncols)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t g = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        double px = 0.0;
        double py = 0.0;
        double pz = 0.0;

        const auto g_x = grid_x[g + grid_offset];
        const auto g_y = grid_y[g + grid_offset];
        const auto g_z = grid_z[g + grid_offset];

        for (uint32_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_info[i + j * nrows + npgtos * nrows * 0];
            const auto fnorm = gto_info[i + j * nrows + npgtos * nrows * 1];
            const auto r_x   = gto_info[i + j * nrows + npgtos * nrows * 2];
            const auto r_y   = gto_info[i + j * nrows + npgtos * nrows * 3];
            const auto r_z   = gto_info[i + j * nrows + npgtos * nrows * 4];

            const auto gr_x = g_x - r_x;
            const auto gr_y = g_y - r_y;
            const auto gr_z = g_z - r_z;

            const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            px += gr_x * fss;
            py += gr_y * fss;
            pz += gr_z * fss;
        }

        // p-1: py
        // p_0: pz
        // p+1: px

        gto_values_p3[g + (i + row_offset) * ncols + nrows * ncols * 0] = py;
        gto_values_p3[g + (i + row_offset) * ncols + nrows * ncols * 1] = pz;
        gto_values_p3[g + (i + row_offset) * ncols + nrows * ncols * 2] = px;
    }
}

__global__ void
cudaLdaValuesDirectRecD(double*        gto_values_d5,
                        const uint32_t row_offset,
                        const double   f2_3,
                        const double*  gto_info,
                        const double*  grid_x,
                        const double*  grid_y,
                        const double*  grid_z,
                        const uint32_t grid_offset,
                        const uint32_t nrows,
                        const uint32_t npgtos,
                        const uint32_t ncols)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t g = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        double dxx = 0.0;
        double dxy = 0.0;
        double dxz = 0.0;
        double dyy = 0.0;
        double dyz = 0.0;
        double dzz = 0.0;

        const auto g_x = grid_x[g + grid_offset];
        const auto g_y = grid_y[g + grid_offset];
        const auto g_z = grid_z[g + grid_offset];

        for (uint32_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_info[i + j * nrows + npgtos * nrows * 0];
            const auto fnorm = gto_info[i + j * nrows + npgtos * nrows * 1];
            const auto r_x   = gto_info[i + j * nrows + npgtos * nrows * 2];
            const auto r_y   = gto_info[i + j * nrows + npgtos * nrows * 3];
            const auto r_z   = gto_info[i + j * nrows + npgtos * nrows * 4];

            const auto gr_x = g_x - r_x;
            const auto gr_y = g_y - r_y;
            const auto gr_z = g_z - r_z;

            const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            dxx += gr_x * gr_x * fss;
            dxy += gr_x * gr_y * fss;
            dxz += gr_x * gr_z * fss;
            dyy += gr_y * gr_y * fss;
            dyz += gr_y * gr_z * fss;
            dzz += gr_z * gr_z * fss;
        }

        // d-2: dxy * f2_3
        // d-1: dyz * f2_3
        // d_0: dzz * 2.0 - dxx - dyy
        // d+1: dxz * f2_3
        // d+2: (dxx - dyy) * 0.5 * f2_3

        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 0] = dxy * f2_3;
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 1] = dyz * f2_3;
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 2] = (dzz * 2.0 - dxx - dyy);
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 3] = dxz * f2_3;
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 4] = (dxx - dyy) * 0.5 * f2_3;
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

                gto_info_ptr[irow + j * ncgtos + npgtos * ncgtos * 0] = fexp;
                gto_info_ptr[irow + j * ncgtos + npgtos * ncgtos * 1] = fnorm;
                gto_info_ptr[irow + j * ncgtos + npgtos * ncgtos * 2] = r_x;
                gto_info_ptr[irow + j * ncgtos + npgtos * ncgtos * 3] = r_y;
                gto_info_ptr[irow + j * ncgtos + npgtos * ncgtos * 4] = r_z;
            }

            irow++;
        }
    }

    return gto_info;
}

static auto
getGtoValuesForLdaDirect(double*                     d_gto_values,
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

    cudaSafe(cudaMemcpy(d_gto_info, gto_info.data(), 5 * nrows * npgtos * sizeof(double), cudaMemcpyHostToDevice));

    // evaluate GTO values on grid points

    dim3 threads_per_block(8, 32);

    dim3 num_blocks((nrows + threads_per_block.x - 1) / threads_per_block.x, (ncols + threads_per_block.y - 1) / threads_per_block.y);

    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        gpu::cudaLdaValuesDirectRecS<<<num_blocks, threads_per_block>>>(d_gto_values,
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
        gpu::cudaLdaValuesDirectRecP<<<num_blocks, threads_per_block>>>(d_gto_values,
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

        gpu::cudaLdaValuesDirectRecD<<<num_blocks, threads_per_block>>>(d_gto_values,
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
        std::string err_ang("gpu::getGtoValuesForLdaDirect: Only implemented for s, p and d-orbitals");

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

    cudaSafe(cudaMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double* d_gaos;

    cudaSafe(cudaMalloc(&d_gaos, naos * max_npoints_per_box * sizeof(double)));

    // coordinates of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z;

    cudaSafe(cudaMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));

    cudaSafe(cudaMemcpy(d_grid_x, xcoords, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_y, ycoords, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_z, zcoords, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));

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

            gpu::getGtoValuesForLdaDirect(d_gaos, row_offset, d_gto_info, gto_block, d_grid_x, d_grid_y, d_grid_z, gridblockpos, npoints, cgto_mask);

            row_offset += static_cast<int64_t>(pre_ao_inds.size());
        }

        cudaSafe(cudaMemcpy(mat_chi.values(), d_gaos, aocount * npoints * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t nu = 0; nu < aocount; nu++)
        {
            std::memcpy(allgtovalues.row(aoinds[nu]) + gridblockpos, mat_chi.row(nu), npoints * sizeof(double));
        }
    }

    cudaSafe(cudaFree(d_gto_info));
    cudaSafe(cudaFree(d_gaos));
    cudaSafe(cudaFree(d_grid_x));
    cudaSafe(cudaFree(d_grid_y));
    cudaSafe(cudaFree(d_grid_z));

    return allgtovalues;
}

static auto
generateDensityForLDA(double*             rho,
                      const CDenseMatrix& densityMatrix,
                      const CDenseMatrix& gtoValues,
                      const double*       d_gto_values,
                      const int64_t       aocount,
                      const int64_t       npoints,
                      CMultiTimer&        timer) -> void
{
    timer.start("Density grid matmul");

    // density matrix: nao x nao
    auto narow = static_cast<uint32_t>(aocount);
    auto nacol = static_cast<uint32_t>(aocount);

    // GTO values: nao x npoints
    // auto nbrow = static_cast<uint32_t>(aocount);
    auto nbcol = static_cast<uint32_t>(npoints);

    // prepare matrices

    double *d_den_mat, *d_mat_F;

    cudaSafe(cudaMalloc(&d_den_mat, narow * nacol * sizeof(double)));
    cudaSafe(cudaMalloc(&d_mat_F, narow * nbcol * sizeof(double)));

    cudaSafe(cudaMemcpy(d_den_mat, densityMatrix.values(), narow * nacol * sizeof(double), cudaMemcpyHostToDevice));

    // use cublas to get multAB(densityMatrix, gtoValues)

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto m = narow, k = nacol, n = nbcol;

    // we want row-major C = A * B but cublas is column-major.
    // so we do C^T = B^T * A^T instead.
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, k, &alpha, d_gto_values, n, d_den_mat, k, &beta, d_mat_F, n));

    cublasDestroy(handle);

    CDenseMatrix mat_F(narow, nbcol);

    cudaSafe(cudaMemcpy(mat_F.values(), d_mat_F, narow * nbcol * sizeof(double), cudaMemcpyDeviceToHost));

    cudaSafe(cudaFree(d_den_mat));
    cudaSafe(cudaFree(d_mat_F));

    timer.stop("Density grid matmul");

    timer.start("Density grid rho");

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    for (int64_t g = 0; g < npoints; g++)
    {
        // rho_alpha
        rho[2 * g + 0] = 0.0;
    }

    for (int64_t nu = 0; nu < aocount; nu++)
    {
        auto nu_offset = nu * npoints;

        for (int64_t g = 0; g < npoints; g++)
        {
            // rho_alpha
            rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
        }
    }

    for (int64_t g = 0; g < npoints; g++)
    {
        // rho_beta
        rho[2 * g + 1] = rho[2 * g + 0];
    }
}

static auto
integratePartialVxcFockForLDA(const double* weights, const CDenseMatrix& gtoValues, const double* vrho, CMultiTimer timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    timer.start("Vxc matrix G");

    auto chi_val = gtoValues.values();

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

    {
        for (int64_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            for (int64_t g = 0; g < npoints; g++)
            {
                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = denblas::multABt(gtoValues, mat_G);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

static auto
integrateVxcFockForLDA(const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CAODensityMatrix& densityMatrix,
                       const CMolecularGrid&   molecularGrid,
                       const CXCFunctional&    xcFunctional,
                       const std::string&      flag) -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

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

    cudaSafe(cudaMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // Kohn-Sham matrix

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double* d_gaos;

    cudaSafe(cudaMalloc(&d_gaos, naos * max_npoints_per_box * sizeof(double)));

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    cudaSafe(cudaMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_w, n_total_grid_points * sizeof(double)));

    cudaSafe(cudaMemcpy(d_grid_x, xcoords, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_y, ycoords, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_z, zcoords, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_w, weights, n_total_grid_points * sizeof(double), cudaMemcpyHostToDevice));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.stop("Preparation");

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        timer.start("GTO pre-screening");

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

        timer.stop("GTO pre-screening");

        if (aocount == 0) continue;

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

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

            gpu::getGtoValuesForLdaDirect(d_gaos, row_offset, d_gto_info, gto_block, d_grid_x, d_grid_y, d_grid_z, gridblockpos, npoints, cgto_mask);

            row_offset += static_cast<int64_t>(pre_ao_inds.size());
        }

        cudaSafe(cudaMemcpy(mat_chi.values(), d_gaos, aocount * npoints * sizeof(double), cudaMemcpyDeviceToHost));

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(densityMatrix, 0, std::string("ALPHA"), aoinds);

            timer.stop("Density matrix slicing");

            gpu::generateDensityForLDA(rho, sub_dens_mat, mat_chi, d_gaos, aocount, npoints, timer);
        }
        else
        {
            // TODO: openshell
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc = gpu::integratePartialVxcFockForLDA(local_weights, mat_chi, vrho, timer);

            timer.start("Vxc matrix dist.");

            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            // TODO: openshell
        }

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int64_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
    }

    cudaSafe(cudaFree(d_gto_info));
    cudaSafe(cudaFree(d_gaos));
    cudaSafe(cudaFree(d_grid_x));
    cudaSafe(cudaFree(d_grid_y));
    cudaSafe(cudaFree(d_grid_z));
    cudaSafe(cudaFree(d_grid_w));

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

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

    /*
    if (xcfuntype == xcfun::gga) return gpu::integrateVxcFockForGGA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    if (xcfuntype == xcfun::mgga) return gpu::integrateVxcFockForMGGA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    std::string errxcfuntype("gpu::integrateVxcFock: Only implemented for LDA/GGA/meta-GGA");

    errors::assertMsgCritical(false, errxcfuntype);
    */

    return CAOKohnShamMatrix();
}

}  // namespace gpu
