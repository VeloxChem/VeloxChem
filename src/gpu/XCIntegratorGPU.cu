//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "GpuRuntime.hpp"


#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "AOIndices.hpp"
#include "ChunkedMemcpyGPU.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "GpuDevices.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "GtoValuesLdaGPU.hpp"
#include "GtoValuesGgaGPU.hpp"
#include "GtoValuesMggaGPU.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"
#include "XCIntegratorGPU.hpp"

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
getSymmetricSubDensityMatrix(double* d_den_mat, const double* d_den_mat_full, const uint32_t naos, const uint32_t* d_ao_inds, const uint32_t aocount)
{
    const uint32_t row = blockDim.y * blockIdx.y + threadIdx.y;
    const uint32_t col = blockDim.x * blockIdx.x + threadIdx.x;

    if ((row < aocount) && (col < aocount))
    {
        const auto row_orig = d_ao_inds[row];
        const auto col_orig = d_ao_inds[col];

        d_den_mat[row * aocount + col] = 0.5 * (d_den_mat_full[row_orig * naos + col_orig] + d_den_mat_full[col_orig * naos + row_orig]);
    }
}

__global__ void
zeroMatrix(double* d_matrix, const uint32_t nrows, const uint32_t ncols)
{
    const uint32_t row = blockDim.y * blockIdx.y + threadIdx.y;
    const uint32_t col = blockDim.x * blockIdx.x + threadIdx.x;

    if ((row < nrows) && (col < ncols))
    {
        d_matrix[row * ncols + col] = 0.0;
    }
}

__global__ void
getDensityOnGrids(double* d_rho, const double* d_mat_F, const double* d_gto_values, const uint32_t aocount, const uint32_t npoints)
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
getDensityOnGridsOpenShell(double*        d_rho,
                           const double*  d_mat_F_a,
                           const double*  d_mat_F_b,
                           const double*  d_gto_values,
                           const uint32_t aocount,
                           const uint32_t npoints)
{
    __shared__ double As_a[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_b[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM][TILE_DIM + 1];

    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    double rho_a = 0.0;
    double rho_b = 0.0;

    for (uint32_t k = 0; k < (aocount + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((k * TILE_DIM + threadIdx.y < aocount) && (g < npoints))
        {
            As_a[threadIdx.y][threadIdx.x] = d_mat_F_a[(k * TILE_DIM + threadIdx.y) * npoints + g];
            As_b[threadIdx.y][threadIdx.x] = d_mat_F_b[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs[threadIdx.y][threadIdx.x] = d_gto_values[(k * TILE_DIM + threadIdx.y) * npoints + g];
        }
        else
        {
            As_a[threadIdx.y][threadIdx.x] = 0.0;
            As_b[threadIdx.y][threadIdx.x] = 0.0;
            Bs[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

        if (threadIdx.y == 0)
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                rho_a += As_a[m][threadIdx.x] * Bs[m][threadIdx.x];
                rho_b += As_b[m][threadIdx.x] * Bs[m][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (g < npoints))
    {
        d_rho[2 * g + 0] = rho_a;
        d_rho[2 * g + 1] = rho_b;
    }
}

__global__ void
getDensitySigmaOnGrids(double*        d_rho,
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

    double rho_a   = 0.0;
    double rho_a_x = 0.0;
    double rho_a_y = 0.0;
    double rho_a_z = 0.0;

    for (uint32_t k = 0; k < (aocount + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((k * TILE_DIM + threadIdx.y < aocount) && (g < npoints))
        {
            As[threadIdx.y][threadIdx.x]   = d_mat_F[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs[threadIdx.y][threadIdx.x]   = d_gto_values[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_x[threadIdx.y][threadIdx.x] = d_gto_values_x[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_y[threadIdx.y][threadIdx.x] = d_gto_values_y[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_z[threadIdx.y][threadIdx.x] = d_gto_values_z[(k * TILE_DIM + threadIdx.y) * npoints + g];
        }
        else
        {
            As[threadIdx.y][threadIdx.x]   = 0.0;
            Bs[threadIdx.y][threadIdx.x]   = 0.0;
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

        if (d_sigma != nullptr)
        {
            d_sigma[3 * g + 0] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
            d_sigma[3 * g + 1] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
            d_sigma[3 * g + 2] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
        }
    }
}

__global__ void
getDensitySigmaOnGridsOpenShell(double*        d_rho,
                                double*        d_rhograd,
                                double*        d_sigma,
                                const double*  d_mat_F_a,
                                const double*  d_mat_F_b,
                                const double*  d_gto_values,
                                const double*  d_gto_values_x,
                                const double*  d_gto_values_y,
                                const double*  d_gto_values_z,
                                const uint32_t aocount,
                                const uint32_t npoints)
{
    __shared__ double As_a[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_b[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs_x[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs_y[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs_z[TILE_DIM][TILE_DIM + 1];

    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    double rho_a   = 0.0;
    double rho_a_x = 0.0;
    double rho_a_y = 0.0;
    double rho_a_z = 0.0;

    double rho_b   = 0.0;
    double rho_b_x = 0.0;
    double rho_b_y = 0.0;
    double rho_b_z = 0.0;

    for (uint32_t k = 0; k < (aocount + TILE_DIM - 1) / TILE_DIM; k++)
    {
        if ((k * TILE_DIM + threadIdx.y < aocount) && (g < npoints))
        {
            As_a[threadIdx.y][threadIdx.x] = d_mat_F_a[(k * TILE_DIM + threadIdx.y) * npoints + g];
            As_b[threadIdx.y][threadIdx.x] = d_mat_F_b[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs[threadIdx.y][threadIdx.x]   = d_gto_values[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_x[threadIdx.y][threadIdx.x] = d_gto_values_x[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_y[threadIdx.y][threadIdx.x] = d_gto_values_y[(k * TILE_DIM + threadIdx.y) * npoints + g];
            Bs_z[threadIdx.y][threadIdx.x] = d_gto_values_z[(k * TILE_DIM + threadIdx.y) * npoints + g];
        }
        else
        {
            As_a[threadIdx.y][threadIdx.x] = 0.0;
            As_b[threadIdx.y][threadIdx.x] = 0.0;
            Bs[threadIdx.y][threadIdx.x]   = 0.0;
            Bs_x[threadIdx.y][threadIdx.x] = 0.0;
            Bs_y[threadIdx.y][threadIdx.x] = 0.0;
            Bs_z[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

        if (threadIdx.y == 0)
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                rho_a += As_a[m][threadIdx.x] * Bs[m][threadIdx.x];
                rho_b += As_b[m][threadIdx.x] * Bs[m][threadIdx.x];

                rho_a_x += 2.0 * (As_a[m][threadIdx.x] * Bs_x[m][threadIdx.x]);
                rho_a_y += 2.0 * (As_a[m][threadIdx.x] * Bs_y[m][threadIdx.x]);
                rho_a_z += 2.0 * (As_a[m][threadIdx.x] * Bs_z[m][threadIdx.x]);
                rho_b_x += 2.0 * (As_b[m][threadIdx.x] * Bs_x[m][threadIdx.x]);
                rho_b_y += 2.0 * (As_b[m][threadIdx.x] * Bs_y[m][threadIdx.x]);
                rho_b_z += 2.0 * (As_b[m][threadIdx.x] * Bs_z[m][threadIdx.x]);
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (g < npoints))
    {
        d_rho[2 * g + 0] = rho_a;
        d_rho[2 * g + 1] = rho_b;

        d_rhograd[6 * g + 0] = rho_a_x;
        d_rhograd[6 * g + 1] = rho_a_y;
        d_rhograd[6 * g + 2] = rho_a_z;
        d_rhograd[6 * g + 3] = rho_b_x;
        d_rhograd[6 * g + 4] = rho_b_y;
        d_rhograd[6 * g + 5] = rho_b_z;

        if (d_sigma != nullptr)
        {
            d_sigma[3 * g + 0] = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z;
            d_sigma[3 * g + 1] = rho_a_x * rho_b_x + rho_a_y * rho_b_y + rho_a_z * rho_b_z;
            d_sigma[3 * g + 2] = rho_b_x * rho_b_x + rho_b_y * rho_b_y + rho_b_z * rho_b_z;
        }
    }
}

__global__ void
getLdaDensityGradientOnGrids(double* d_dengrad_x,
                             double* d_dengrad_y,
                             double* d_dengrad_z,
                             const double* d_mat_F,
                             const double* d_gto_values_x,
                             const double* d_gto_values_y,
                             const double* d_gto_values_z,
                             const uint32_t aocount,
                             const uint32_t* d_ao_to_atom_ids,
                             const uint32_t natoms,
                             const uint32_t npoints)
{
    const uint32_t a = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((a < natoms) && (g < npoints))
    {
        double dgx = 0.0, dgy = 0.0, dgz = 0.0;

        for (uint32_t i = 0; i < aocount; i++)
        {
            if (a == d_ao_to_atom_ids[i])
            {
                dgx -= 2.0 * d_mat_F[i * npoints + g] * d_gto_values_x[i * npoints + g];
                dgy -= 2.0 * d_mat_F[i * npoints + g] * d_gto_values_y[i * npoints + g];
                dgz -= 2.0 * d_mat_F[i * npoints + g] * d_gto_values_z[i * npoints + g];
            }
        }

        d_dengrad_x[a * npoints + g] = dgx;
        d_dengrad_y[a * npoints + g] = dgy;
        d_dengrad_z[a * npoints + g] = dgz;
    }
}

__global__ void
getLdaDensityGradientOnGridsOpenShell(double* d_dengrad_x_a,
                                      double* d_dengrad_x_b,
                                      double* d_dengrad_y_a,
                                      double* d_dengrad_y_b,
                                      double* d_dengrad_z_a,
                                      double* d_dengrad_z_b,
                                      const double* d_mat_F_a,
                                      const double* d_mat_F_b,
                                      const double* d_gto_values_x,
                                      const double* d_gto_values_y,
                                      const double* d_gto_values_z,
                                      const uint32_t aocount,
                                      const uint32_t* d_ao_to_atom_ids,
                                      const uint32_t natoms,
                                      const uint32_t npoints)
{
    const uint32_t a = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((a < natoms) && (g < npoints))
    {
        double dgx_a = 0.0, dgy_a = 0.0, dgz_a = 0.0;
        double dgx_b = 0.0, dgy_b = 0.0, dgz_b = 0.0;

        for (uint32_t i = 0; i < aocount; i++)
        {
            if (a == d_ao_to_atom_ids[i])
            {
                dgx_a -= 2.0 * d_mat_F_a[i * npoints + g] * d_gto_values_x[i * npoints + g];
                dgy_a -= 2.0 * d_mat_F_a[i * npoints + g] * d_gto_values_y[i * npoints + g];
                dgz_a -= 2.0 * d_mat_F_a[i * npoints + g] * d_gto_values_z[i * npoints + g];

                dgx_b -= 2.0 * d_mat_F_b[i * npoints + g] * d_gto_values_x[i * npoints + g];
                dgy_b -= 2.0 * d_mat_F_b[i * npoints + g] * d_gto_values_y[i * npoints + g];
                dgz_b -= 2.0 * d_mat_F_b[i * npoints + g] * d_gto_values_z[i * npoints + g];
            }
        }

        d_dengrad_x_a[a * npoints + g] = dgx_a;
        d_dengrad_y_a[a * npoints + g] = dgy_a;
        d_dengrad_z_a[a * npoints + g] = dgz_a;

        d_dengrad_x_b[a * npoints + g] = dgx_b;
        d_dengrad_y_b[a * npoints + g] = dgy_b;
        d_dengrad_z_b[a * npoints + g] = dgz_b;
    }
}

__global__ void
getGgaDensityGradientOnGrids(double* d_dengrad_x,
                             double* d_dengrad_y,
                             double* d_dengrad_z,
                             double* d_dengrad_xx,
                             double* d_dengrad_xy,
                             double* d_dengrad_xz,
                             double* d_dengrad_yx,
                             double* d_dengrad_yy,
                             double* d_dengrad_yz,
                             double* d_dengrad_zx,
                             double* d_dengrad_zy,
                             double* d_dengrad_zz,
                             const double* d_mat_F,
                             const double* d_mat_F_x,
                             const double* d_mat_F_y,
                             const double* d_mat_F_z,
                             const double* d_gto_values_x,
                             const double* d_gto_values_y,
                             const double* d_gto_values_z,
                             const double* d_gto_values_xx,
                             const double* d_gto_values_xy,
                             const double* d_gto_values_xz,
                             const double* d_gto_values_yy,
                             const double* d_gto_values_yz,
                             const double* d_gto_values_zz,
                             const uint32_t aocount,
                             const uint32_t* d_ao_to_atom_ids,
                             const uint32_t natoms,
                             const uint32_t npoints)
{
    const uint32_t a = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((a < natoms) && (g < npoints))
    {
        double dgx = 0.0, dgy = 0.0, dgz = 0.0;

        double dgxx = 0.0, dgxy = 0.0, dgxz = 0.0;
        double dgyx = 0.0, dgyy = 0.0, dgyz = 0.0;
        double dgzx = 0.0, dgzy = 0.0, dgzz = 0.0;

        for (uint32_t i = 0; i < aocount; i++)
        {
            if (a == d_ao_to_atom_ids[i])
            {
                const uint32_t ig = i * npoints + g;

                dgx -= 2.0 * d_mat_F[ig] * d_gto_values_x[ig];
                dgy -= 2.0 * d_mat_F[ig] * d_gto_values_y[ig];
                dgz -= 2.0 * d_mat_F[ig] * d_gto_values_z[ig];

                dgxx -= 2.0 * (d_mat_F_x[ig] * d_gto_values_x[ig] + d_mat_F[ig] * d_gto_values_xx[ig]);
                dgxy -= 2.0 * (d_mat_F_x[ig] * d_gto_values_y[ig] + d_mat_F[ig] * d_gto_values_xy[ig]);
                dgxz -= 2.0 * (d_mat_F_x[ig] * d_gto_values_z[ig] + d_mat_F[ig] * d_gto_values_xz[ig]);

                dgyx -= 2.0 * (d_mat_F_y[ig] * d_gto_values_x[ig] + d_mat_F[ig] * d_gto_values_xy[ig]);
                dgyy -= 2.0 * (d_mat_F_y[ig] * d_gto_values_y[ig] + d_mat_F[ig] * d_gto_values_yy[ig]);
                dgyz -= 2.0 * (d_mat_F_y[ig] * d_gto_values_z[ig] + d_mat_F[ig] * d_gto_values_yz[ig]);

                dgzx -= 2.0 * (d_mat_F_z[ig] * d_gto_values_x[ig] + d_mat_F[ig] * d_gto_values_xz[ig]);
                dgzy -= 2.0 * (d_mat_F_z[ig] * d_gto_values_y[ig] + d_mat_F[ig] * d_gto_values_yz[ig]);
                dgzz -= 2.0 * (d_mat_F_z[ig] * d_gto_values_z[ig] + d_mat_F[ig] * d_gto_values_zz[ig]);
            }
        }

        d_dengrad_x[a * npoints + g] = dgx;
        d_dengrad_y[a * npoints + g] = dgy;
        d_dengrad_z[a * npoints + g] = dgz;

        d_dengrad_xx[a * npoints + g] = dgxx;
        d_dengrad_xy[a * npoints + g] = dgxy;
        d_dengrad_xz[a * npoints + g] = dgxz;

        d_dengrad_yx[a * npoints + g] = dgyx;
        d_dengrad_yy[a * npoints + g] = dgyy;
        d_dengrad_yz[a * npoints + g] = dgyz;

        d_dengrad_zx[a * npoints + g] = dgzx;
        d_dengrad_zy[a * npoints + g] = dgzy;
        d_dengrad_zz[a * npoints + g] = dgzz;
    }
}

__global__ void
getGgaDensityGradientOnGridsOpenShell(double* d_dengrad_x_a,
                                      double* d_dengrad_x_b,
                                      double* d_dengrad_y_a,
                                      double* d_dengrad_y_b,
                                      double* d_dengrad_z_a,
                                      double* d_dengrad_z_b,
                                      double* d_dengrad_xx_a,
                                      double* d_dengrad_xx_b,
                                      double* d_dengrad_xy_a,
                                      double* d_dengrad_xy_b,
                                      double* d_dengrad_xz_a,
                                      double* d_dengrad_xz_b,
                                      double* d_dengrad_yx_a,
                                      double* d_dengrad_yx_b,
                                      double* d_dengrad_yy_a,
                                      double* d_dengrad_yy_b,
                                      double* d_dengrad_yz_a,
                                      double* d_dengrad_yz_b,
                                      double* d_dengrad_zx_a,
                                      double* d_dengrad_zx_b,
                                      double* d_dengrad_zy_a,
                                      double* d_dengrad_zy_b,
                                      double* d_dengrad_zz_a,
                                      double* d_dengrad_zz_b,
                                      const double* d_mat_F_a,
                                      const double* d_mat_F_b,
                                      const double* d_mat_F_x_a,
                                      const double* d_mat_F_x_b,
                                      const double* d_mat_F_y_a,
                                      const double* d_mat_F_y_b,
                                      const double* d_mat_F_z_a,
                                      const double* d_mat_F_z_b,
                                      const double* d_gto_values_x,
                                      const double* d_gto_values_y,
                                      const double* d_gto_values_z,
                                      const double* d_gto_values_xx,
                                      const double* d_gto_values_xy,
                                      const double* d_gto_values_xz,
                                      const double* d_gto_values_yy,
                                      const double* d_gto_values_yz,
                                      const double* d_gto_values_zz,
                                      const uint32_t aocount,
                                      const uint32_t* d_ao_to_atom_ids,
                                      const uint32_t natoms,
                                      const uint32_t npoints)
{
    const uint32_t a = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((a < natoms) && (g < npoints))
    {
        double dgx_a = 0.0, dgy_a = 0.0, dgz_a = 0.0;
        double dgx_b = 0.0, dgy_b = 0.0, dgz_b = 0.0;

        double dgxx_a = 0.0, dgxy_a = 0.0, dgxz_a = 0.0;
        double dgyx_a = 0.0, dgyy_a = 0.0, dgyz_a = 0.0;
        double dgzx_a = 0.0, dgzy_a = 0.0, dgzz_a = 0.0;

        double dgxx_b = 0.0, dgxy_b = 0.0, dgxz_b = 0.0;
        double dgyx_b = 0.0, dgyy_b = 0.0, dgyz_b = 0.0;
        double dgzx_b = 0.0, dgzy_b = 0.0, dgzz_b = 0.0;

        for (uint32_t i = 0; i < aocount; i++)
        {
            if (a == d_ao_to_atom_ids[i])
            {
                const uint32_t ig = i * npoints + g;

                dgx_a -= 2.0 * d_mat_F_a[ig] * d_gto_values_x[ig];
                dgy_a -= 2.0 * d_mat_F_a[ig] * d_gto_values_y[ig];
                dgz_a -= 2.0 * d_mat_F_a[ig] * d_gto_values_z[ig];

                dgx_b -= 2.0 * d_mat_F_b[ig] * d_gto_values_x[ig];
                dgy_b -= 2.0 * d_mat_F_b[ig] * d_gto_values_y[ig];
                dgz_b -= 2.0 * d_mat_F_b[ig] * d_gto_values_z[ig];

                dgxx_a -= 2.0 * (d_mat_F_x_a[ig] * d_gto_values_x[ig] + d_mat_F_a[ig] * d_gto_values_xx[ig]);
                dgxy_a -= 2.0 * (d_mat_F_x_a[ig] * d_gto_values_y[ig] + d_mat_F_a[ig] * d_gto_values_xy[ig]);
                dgxz_a -= 2.0 * (d_mat_F_x_a[ig] * d_gto_values_z[ig] + d_mat_F_a[ig] * d_gto_values_xz[ig]);

                dgxx_b -= 2.0 * (d_mat_F_x_b[ig] * d_gto_values_x[ig] + d_mat_F_b[ig] * d_gto_values_xx[ig]);
                dgxy_b -= 2.0 * (d_mat_F_x_b[ig] * d_gto_values_y[ig] + d_mat_F_b[ig] * d_gto_values_xy[ig]);
                dgxz_b -= 2.0 * (d_mat_F_x_b[ig] * d_gto_values_z[ig] + d_mat_F_b[ig] * d_gto_values_xz[ig]);

                dgyx_a -= 2.0 * (d_mat_F_y_a[ig] * d_gto_values_x[ig] + d_mat_F_a[ig] * d_gto_values_xy[ig]);
                dgyy_a -= 2.0 * (d_mat_F_y_a[ig] * d_gto_values_y[ig] + d_mat_F_a[ig] * d_gto_values_yy[ig]);
                dgyz_a -= 2.0 * (d_mat_F_y_a[ig] * d_gto_values_z[ig] + d_mat_F_a[ig] * d_gto_values_yz[ig]);

                dgyx_b -= 2.0 * (d_mat_F_y_b[ig] * d_gto_values_x[ig] + d_mat_F_b[ig] * d_gto_values_xy[ig]);
                dgyy_b -= 2.0 * (d_mat_F_y_b[ig] * d_gto_values_y[ig] + d_mat_F_b[ig] * d_gto_values_yy[ig]);
                dgyz_b -= 2.0 * (d_mat_F_y_b[ig] * d_gto_values_z[ig] + d_mat_F_b[ig] * d_gto_values_yz[ig]);

                dgzx_a -= 2.0 * (d_mat_F_z_a[ig] * d_gto_values_x[ig] + d_mat_F_a[ig] * d_gto_values_xz[ig]);
                dgzy_a -= 2.0 * (d_mat_F_z_a[ig] * d_gto_values_y[ig] + d_mat_F_a[ig] * d_gto_values_yz[ig]);
                dgzz_a -= 2.0 * (d_mat_F_z_a[ig] * d_gto_values_z[ig] + d_mat_F_a[ig] * d_gto_values_zz[ig]);

                dgzx_b -= 2.0 * (d_mat_F_z_b[ig] * d_gto_values_x[ig] + d_mat_F_b[ig] * d_gto_values_xz[ig]);
                dgzy_b -= 2.0 * (d_mat_F_z_b[ig] * d_gto_values_y[ig] + d_mat_F_b[ig] * d_gto_values_yz[ig]);
                dgzz_b -= 2.0 * (d_mat_F_z_b[ig] * d_gto_values_z[ig] + d_mat_F_b[ig] * d_gto_values_zz[ig]);
            }
        }

        d_dengrad_x_a[a * npoints + g] = dgx_a;
        d_dengrad_y_a[a * npoints + g] = dgy_a;
        d_dengrad_z_a[a * npoints + g] = dgz_a;

        d_dengrad_x_b[a * npoints + g] = dgx_b;
        d_dengrad_y_b[a * npoints + g] = dgy_b;
        d_dengrad_z_b[a * npoints + g] = dgz_b;

        d_dengrad_xx_a[a * npoints + g] = dgxx_a;
        d_dengrad_xy_a[a * npoints + g] = dgxy_a;
        d_dengrad_xz_a[a * npoints + g] = dgxz_a;

        d_dengrad_xx_b[a * npoints + g] = dgxx_b;
        d_dengrad_xy_b[a * npoints + g] = dgxy_b;
        d_dengrad_xz_b[a * npoints + g] = dgxz_b;

        d_dengrad_yx_a[a * npoints + g] = dgyx_a;
        d_dengrad_yy_a[a * npoints + g] = dgyy_a;
        d_dengrad_yz_a[a * npoints + g] = dgyz_a;

        d_dengrad_yx_b[a * npoints + g] = dgyx_b;
        d_dengrad_yy_b[a * npoints + g] = dgyy_b;
        d_dengrad_yz_b[a * npoints + g] = dgyz_b;

        d_dengrad_zx_a[a * npoints + g] = dgzx_a;
        d_dengrad_zy_a[a * npoints + g] = dgzy_a;
        d_dengrad_zz_a[a * npoints + g] = dgzz_a;

        d_dengrad_zx_b[a * npoints + g] = dgzx_b;
        d_dengrad_zy_b[a * npoints + g] = dgzy_b;
        d_dengrad_zz_b[a * npoints + g] = dgzz_b;
    }
}

__global__ void
getLdaVxcMatrixG(double*        d_mat_G,
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
getLdaVxcMatrixGOpenShell(double*        d_mat_G_a,
                          double*        d_mat_G_b,
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
        d_mat_G_a[g + i * npoints] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0] * d_gto_values[g + i * npoints];
        d_mat_G_b[g + i * npoints] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 1] * d_gto_values[g + i * npoints];
    }
}

__global__ void
getLdaFxcMatrixG(double*        d_mat_G,
                 const double*  d_grid_w,
                 const uint32_t grid_offset,
                 const uint32_t npoints,
                 const double*  d_gto_values,
                 const uint32_t aocount,
                 const double*  d_rhow,
                 const uint32_t dim_rhow,
                 const double*  d_v2rho2,
                 const uint32_t dim_v2rho2)
{
    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((i < aocount) && (g < npoints))
    {
        auto rhow_a = d_rhow[dim_rhow * g + 0];
        auto rhow_b = d_rhow[dim_rhow * g + 1];

        // functional derivatives

        // second-order

        auto v2rho2_aa = d_v2rho2[dim_v2rho2 * g + 0];
        auto v2rho2_ab = d_v2rho2[dim_v2rho2 * g + 1];

        d_mat_G[g + i * npoints] = d_grid_w[g + grid_offset] * (v2rho2_aa * rhow_a + v2rho2_ab * rhow_b) * d_gto_values[g + i * npoints];
    }
}

__global__ void
getLdaVxcGradient(double*        d_mol_grad,
                  const double*  d_grid_w,
                  const uint32_t grid_offset,
                  const uint32_t npoints,
                  const double*  d_dengrad_x,
                  const double*  d_dengrad_y,
                  const double*  d_dengrad_z,
                  const uint32_t natoms,
                  const double*  d_vrho)
{
    __shared__ double As_x[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_y[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_z[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM + 1];

    const uint32_t a = blockDim.x * blockIdx.x + threadIdx.x;

    double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

    for (uint32_t k = 0; k < (npoints + TILE_DIM - 1) / TILE_DIM; k++)
    {
        const uint32_t g = k * TILE_DIM + threadIdx.y;

        if ((g < npoints) && (a < natoms))
        {
            As_x[threadIdx.y][threadIdx.x] = d_dengrad_x[a * npoints + g];
            As_y[threadIdx.y][threadIdx.x] = d_dengrad_y[a * npoints + g];
            As_z[threadIdx.y][threadIdx.x] = d_dengrad_z[a * npoints + g];

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0];
            }
        }
        else
        {
            As_x[threadIdx.y][threadIdx.x] = 0.0;
            As_y[threadIdx.y][threadIdx.x] = 0.0;
            As_z[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x[m][threadIdx.x] * Bs[m];
                gatmy += As_y[m][threadIdx.x] * Bs[m];
                gatmz += As_z[m][threadIdx.x] * Bs[m];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (a < natoms))
    {
        d_mol_grad[a * 3 + 0] += 2.0 * gatmx;
        d_mol_grad[a * 3 + 1] += 2.0 * gatmy;
        d_mol_grad[a * 3 + 2] += 2.0 * gatmz;
    }
}

__global__ void
getLdaVxcGradientOpenShell(double*        d_mol_grad,
                           const double*  d_grid_w,
                           const uint32_t grid_offset,
                           const uint32_t npoints,
                           const double*  d_dengrad_x_a,
                           const double*  d_dengrad_x_b,
                           const double*  d_dengrad_y_a,
                           const double*  d_dengrad_y_b,
                           const double*  d_dengrad_z_a,
                           const double*  d_dengrad_z_b,
                           const uint32_t natoms,
                           const double*  d_vrho)
{
    __shared__ double As_x_a[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_y_a[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_z_a[TILE_DIM][TILE_DIM + 1];

    __shared__ double As_x_b[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_y_b[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_z_b[TILE_DIM][TILE_DIM + 1];

    __shared__ double Bs_a[TILE_DIM + 1];
    __shared__ double Bs_b[TILE_DIM + 1];

    const uint32_t a = blockDim.x * blockIdx.x + threadIdx.x;

    double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

    for (uint32_t k = 0; k < (npoints + TILE_DIM - 1) / TILE_DIM; k++)
    {
        const uint32_t g = k * TILE_DIM + threadIdx.y;

        if ((g < npoints) && (a < natoms))
        {
            As_x_a[threadIdx.y][threadIdx.x] = d_dengrad_x_a[a * npoints + g];
            As_y_a[threadIdx.y][threadIdx.x] = d_dengrad_y_a[a * npoints + g];
            As_z_a[threadIdx.y][threadIdx.x] = d_dengrad_z_a[a * npoints + g];

            As_x_b[threadIdx.y][threadIdx.x] = d_dengrad_x_b[a * npoints + g];
            As_y_b[threadIdx.y][threadIdx.x] = d_dengrad_y_b[a * npoints + g];
            As_z_b[threadIdx.y][threadIdx.x] = d_dengrad_z_b[a * npoints + g];

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0];
                Bs_b[threadIdx.y] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 1];
            }
        }
        else
        {
            As_x_a[threadIdx.y][threadIdx.x] = 0.0;
            As_y_a[threadIdx.y][threadIdx.x] = 0.0;
            As_z_a[threadIdx.y][threadIdx.x] = 0.0;

            As_x_b[threadIdx.y][threadIdx.x] = 0.0;
            As_y_b[threadIdx.y][threadIdx.x] = 0.0;
            As_z_b[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = 0.0;
                Bs_b[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x_a[m][threadIdx.x] * Bs_a[m];
                gatmy += As_y_a[m][threadIdx.x] * Bs_a[m];
                gatmz += As_z_a[m][threadIdx.x] * Bs_a[m];

                gatmx += As_x_b[m][threadIdx.x] * Bs_b[m];
                gatmy += As_y_b[m][threadIdx.x] * Bs_b[m];
                gatmz += As_z_b[m][threadIdx.x] * Bs_b[m];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (a < natoms))
    {
        d_mol_grad[a * 3 + 0] += gatmx;
        d_mol_grad[a * 3 + 1] += gatmy;
        d_mol_grad[a * 3 + 2] += gatmz;
    }
}

__global__ void
getGgaVxcMatrixG(double*        d_mat_G,
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
getGgaVxcMatrixGOpenShell(double*        d_mat_G_a,
                          double*        d_mat_G_b,
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
        const double vxa = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 0] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 3];
        const double vya = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 1] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 4];
        const double vza = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 2] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 5];

        const double vxb = 2.0 * d_vsigma[3 * g + 2] * d_rhograd[6 * g + 3] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 0];
        const double vyb = 2.0 * d_vsigma[3 * g + 2] * d_rhograd[6 * g + 4] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 1];
        const double vzb = 2.0 * d_vsigma[3 * g + 2] * d_rhograd[6 * g + 5] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 2];

        // adding G_a_val and 2 * G_a_gga_val together (only works for ground state Vxc)
        d_mat_G_a[g + i * npoints] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0] * d_gto_values[g + i * npoints] +

                                     2.0 * (d_grid_w[g + grid_offset] * (vxa * d_gto_values_x[g + i * npoints] +
                                                                         vya * d_gto_values_y[g + i * npoints] +
                                                                         vza * d_gto_values_z[g + i * npoints]));

        // adding G_b_val and 2 * G_b_gga_val together (only works for ground state Vxc)
        d_mat_G_b[g + i * npoints] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 1] * d_gto_values[g + i * npoints] +

                                     2.0 * (d_grid_w[g + grid_offset] * (vxb * d_gto_values_x[g + i * npoints] +
                                                                         vyb * d_gto_values_y[g + i * npoints] +
                                                                         vzb * d_gto_values_z[g + i * npoints]));
    }
}

__global__ void
getGgaFxcMatrixG(double*        d_mat_G,
                 double*        d_mat_G_gga,
                 const double*  d_grid_w,
                 const uint32_t grid_offset,
                 const uint32_t npoints,
                 const double*  d_gto_values,
                 const double*  d_gto_values_x,
                 const double*  d_gto_values_y,
                 const double*  d_gto_values_z,
                 const uint32_t aocount,
                 const double*  d_rhow,
                 const double*  d_rhograd,
                 const double*  d_rhowgrad,
                 const uint32_t dim_rho,
                 const double*  d_vsigma,
                 const uint32_t dim_vsigma,
                 const double*  d_v2rho2,
                 const uint32_t dim_v2rho2,
                 const double*  d_v2rhosigma,
                 const uint32_t dim_v2rhosigma,
                 const double*  d_v2sigma2,
                 const uint32_t dim_v2sigma2)
{
    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t g = blockIdx.x * blockDim.x + threadIdx.x;

    if ((i < aocount) && (g < npoints))
    {
        const auto w = d_grid_w[g + grid_offset];

        const auto rhow_a = d_rhow[2 * g + 0];
        const auto rhow_b = d_rhow[2 * g + 1];

        const auto rhowgrad_ax = d_rhowgrad[6 * g + 0];
        const auto rhowgrad_ay = d_rhowgrad[6 * g + 1];
        const auto rhowgrad_az = d_rhowgrad[6 * g + 2];
        const auto rhowgrad_bx = d_rhowgrad[6 * g + 3];
        const auto rhowgrad_by = d_rhowgrad[6 * g + 4];
        const auto rhowgrad_bz = d_rhowgrad[6 * g + 5];

        const auto rhograd_ax = d_rhograd[6 * g + 0];
        const auto rhograd_ay = d_rhograd[6 * g + 1];
        const auto rhograd_az = d_rhograd[6 * g + 2];
        const auto rhograd_bx = d_rhograd[6 * g + 3];
        const auto rhograd_by = d_rhograd[6 * g + 4];
        const auto rhograd_bz = d_rhograd[6 * g + 5];

        const auto grhow_grho_aa = 2.0 * (rhowgrad_ax * rhograd_ax + rhowgrad_ay * rhograd_ay + rhowgrad_az * rhograd_az);

        const auto grhow_grho_bb = 2.0 * (rhowgrad_bx * rhograd_bx + rhowgrad_by * rhograd_by + rhowgrad_bz * rhograd_bz);

        const auto grhow_grho_ab = (rhowgrad_ax * rhograd_bx + rhowgrad_ay * rhograd_by + rhowgrad_az * rhograd_bz +
                                    rhowgrad_bx * rhograd_ax + rhowgrad_by * rhograd_ay + rhowgrad_bz * rhograd_az);

        // functional derivatives

        // first-order

        const auto vsigma_a = d_vsigma[3 * g + 0];
        const auto vsigma_c = d_vsigma[3 * g + 1];

        // second-order

        const auto v2rho2_aa = d_v2rho2[3 * g + 0];
        const auto v2rho2_ab = d_v2rho2[3 * g + 1];

        const auto v2rhosigma_aa = d_v2rhosigma[6 * g + 0];
        const auto v2rhosigma_ac = d_v2rhosigma[6 * g + 1];
        const auto v2rhosigma_ab = d_v2rhosigma[6 * g + 2];
        const auto v2rhosigma_ba = d_v2rhosigma[6 * g + 3];
        const auto v2rhosigma_bc = d_v2rhosigma[6 * g + 4];

        const auto v2sigma2_aa = d_v2sigma2[6 * g + 0];
        const auto v2sigma2_ac = d_v2sigma2[6 * g + 1];
        const auto v2sigma2_ab = d_v2sigma2[6 * g + 2];
        const auto v2sigma2_cc = d_v2sigma2[6 * g + 3];
        const auto v2sigma2_cb = d_v2sigma2[6 * g + 4];

        // scalar contribution

        const auto f_0 = v2rho2_aa * rhow_a + v2rho2_ab * rhow_b + v2rhosigma_aa * grhow_grho_aa + v2rhosigma_ac * grhow_grho_ab + v2rhosigma_ab * grhow_grho_bb;

        d_mat_G[g + i * npoints] = w * f_0 * d_gto_values[g + i * npoints];

        // vector contribution

        const auto f_aa = v2rhosigma_aa * rhow_a + v2rhosigma_ba * rhow_b + v2sigma2_aa * grhow_grho_aa + v2sigma2_ac * grhow_grho_ab + v2sigma2_ab * grhow_grho_bb;
        const auto f_ab = v2rhosigma_ac * rhow_a + v2rhosigma_bc * rhow_b + v2sigma2_ac * grhow_grho_aa + v2sigma2_cc * grhow_grho_ab + v2sigma2_cb * grhow_grho_bb;

        double xcomp = 2.0 * f_aa * rhograd_ax + f_ab * rhograd_bx;
        double ycomp = 2.0 * f_aa * rhograd_ay + f_ab * rhograd_by;
        double zcomp = 2.0 * f_aa * rhograd_az + f_ab * rhograd_bz;

        xcomp += 2.0 * vsigma_a * rhowgrad_ax + vsigma_c * rhowgrad_bx;
        ycomp += 2.0 * vsigma_a * rhowgrad_ay + vsigma_c * rhowgrad_by;
        zcomp += 2.0 * vsigma_a * rhowgrad_az + vsigma_c * rhowgrad_bz;

        d_mat_G_gga[g + i * npoints] = w * (xcomp * d_gto_values_x[g + i * npoints] +
                                            ycomp * d_gto_values_y[g + i * npoints] +
                                            zcomp * d_gto_values_z[g + i * npoints]);
    }
}

__global__ void
getGgaVxcGradient(double*        d_mol_grad,
                  const double*  d_grid_w,
                  const uint32_t grid_offset,
                  const uint32_t npoints,
                  const double*  d_dengrad_x,
                  const double*  d_dengrad_y,
                  const double*  d_dengrad_z,
                  const double*  d_dengrad_xx,
                  const double*  d_dengrad_xy,
                  const double*  d_dengrad_xz,
                  const double*  d_dengrad_yx,
                  const double*  d_dengrad_yy,
                  const double*  d_dengrad_yz,
                  const double*  d_dengrad_zx,
                  const double*  d_dengrad_zy,
                  const double*  d_dengrad_zz,
                  const uint32_t natoms,
                  const double*  d_vrho,
                  const double*  d_vsigma,
                  const double*  d_rhograd)
{
    __shared__ double As_x[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_y[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_z[TILE_DIM][TILE_DIM + 1];
    __shared__ double Bs[TILE_DIM + 1];

    const uint32_t a = blockDim.x * blockIdx.x + threadIdx.x;

    double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

    for (uint32_t k = 0; k < (npoints + TILE_DIM - 1) / TILE_DIM; k++)
    {
        const uint32_t g = k * TILE_DIM + threadIdx.y;

        // vrho contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x[threadIdx.y][threadIdx.x] = d_dengrad_x[a * npoints + g];
            As_y[threadIdx.y][threadIdx.x] = d_dengrad_y[a * npoints + g];
            As_z[threadIdx.y][threadIdx.x] = d_dengrad_z[a * npoints + g];

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0];
            }
        }
        else
        {
            As_x[threadIdx.y][threadIdx.x] = 0.0;
            As_y[threadIdx.y][threadIdx.x] = 0.0;
            As_z[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x[m][threadIdx.x] * Bs[m];
                gatmy += As_y[m][threadIdx.x] * Bs[m];
                gatmz += As_z[m][threadIdx.x] * Bs[m];
            }
        }

        __syncthreads();

        // vx contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x[threadIdx.y][threadIdx.x] = d_dengrad_xx[a * npoints + g];
            As_y[threadIdx.y][threadIdx.x] = d_dengrad_xy[a * npoints + g];
            As_z[threadIdx.y][threadIdx.x] = d_dengrad_xz[a * npoints + g];

            if (threadIdx.x == 0)
            {
                const double vx = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 0] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 3];

                Bs[threadIdx.y] = d_grid_w[g + grid_offset] * vx;
            }
        }
        else
        {
            As_x[threadIdx.y][threadIdx.x] = 0.0;
            As_y[threadIdx.y][threadIdx.x] = 0.0;
            As_z[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if (threadIdx.y == 0)
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x[m][threadIdx.x] * Bs[m];
                gatmy += As_y[m][threadIdx.x] * Bs[m];
                gatmz += As_z[m][threadIdx.x] * Bs[m];
            }
        }

        __syncthreads();

        // vy contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x[threadIdx.y][threadIdx.x] = d_dengrad_yx[a * npoints + g];
            As_y[threadIdx.y][threadIdx.x] = d_dengrad_yy[a * npoints + g];
            As_z[threadIdx.y][threadIdx.x] = d_dengrad_yz[a * npoints + g];

            if (threadIdx.x == 0)
            {
                const double vy = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 1] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 4];

                Bs[threadIdx.y] = d_grid_w[g + grid_offset] * vy;
            }
        }
        else
        {
            As_x[threadIdx.y][threadIdx.x] = 0.0;
            As_y[threadIdx.y][threadIdx.x] = 0.0;
            As_z[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x[m][threadIdx.x] * Bs[m];
                gatmy += As_y[m][threadIdx.x] * Bs[m];
                gatmz += As_z[m][threadIdx.x] * Bs[m];
            }
        }

        __syncthreads();

        // vz contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x[threadIdx.y][threadIdx.x] = d_dengrad_zx[a * npoints + g];
            As_y[threadIdx.y][threadIdx.x] = d_dengrad_zy[a * npoints + g];
            As_z[threadIdx.y][threadIdx.x] = d_dengrad_zz[a * npoints + g];

            if (threadIdx.x == 0)
            {
                const double vz = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 2] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 5];

                Bs[threadIdx.y] = d_grid_w[g + grid_offset] * vz;
            }
        }
        else
        {
            As_x[threadIdx.y][threadIdx.x] = 0.0;
            As_y[threadIdx.y][threadIdx.x] = 0.0;
            As_z[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x[m][threadIdx.x] * Bs[m];
                gatmy += As_y[m][threadIdx.x] * Bs[m];
                gatmz += As_z[m][threadIdx.x] * Bs[m];
            }
        }

        __syncthreads();
    }

    // accumulate gradient

    if ((threadIdx.y == 0) && (a < natoms))
    {
        d_mol_grad[a * 3 + 0] += 2.0 * gatmx;
        d_mol_grad[a * 3 + 1] += 2.0 * gatmy;
        d_mol_grad[a * 3 + 2] += 2.0 * gatmz;
    }
}

__global__ void
getGgaVxcGradientOpenShell(double*        d_mol_grad,
                           const double*  d_grid_w,
                           const uint32_t grid_offset,
                           const uint32_t npoints,
                           const double*  d_dengrad_x_a,
                           const double*  d_dengrad_x_b,
                           const double*  d_dengrad_y_a,
                           const double*  d_dengrad_y_b,
                           const double*  d_dengrad_z_a,
                           const double*  d_dengrad_z_b,
                           const double*  d_dengrad_xx_a,
                           const double*  d_dengrad_xx_b,
                           const double*  d_dengrad_xy_a,
                           const double*  d_dengrad_xy_b,
                           const double*  d_dengrad_xz_a,
                           const double*  d_dengrad_xz_b,
                           const double*  d_dengrad_yx_a,
                           const double*  d_dengrad_yx_b,
                           const double*  d_dengrad_yy_a,
                           const double*  d_dengrad_yy_b,
                           const double*  d_dengrad_yz_a,
                           const double*  d_dengrad_yz_b,
                           const double*  d_dengrad_zx_a,
                           const double*  d_dengrad_zx_b,
                           const double*  d_dengrad_zy_a,
                           const double*  d_dengrad_zy_b,
                           const double*  d_dengrad_zz_a,
                           const double*  d_dengrad_zz_b,
                           const uint32_t natoms,
                           const double*  d_vrho,
                           const double*  d_vsigma,
                           const double*  d_rhograd)
{
    __shared__ double As_x_a[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_y_a[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_z_a[TILE_DIM][TILE_DIM + 1];

    __shared__ double As_x_b[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_y_b[TILE_DIM][TILE_DIM + 1];
    __shared__ double As_z_b[TILE_DIM][TILE_DIM + 1];

    __shared__ double Bs_a[TILE_DIM + 1];
    __shared__ double Bs_b[TILE_DIM + 1];

    const uint32_t a = blockDim.x * blockIdx.x + threadIdx.x;

    double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

    for (uint32_t k = 0; k < (npoints + TILE_DIM - 1) / TILE_DIM; k++)
    {
        const uint32_t g = k * TILE_DIM + threadIdx.y;

        // vrho contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x_a[threadIdx.y][threadIdx.x] = d_dengrad_x_a[a * npoints + g];
            As_y_a[threadIdx.y][threadIdx.x] = d_dengrad_y_a[a * npoints + g];
            As_z_a[threadIdx.y][threadIdx.x] = d_dengrad_z_a[a * npoints + g];

            As_x_b[threadIdx.y][threadIdx.x] = d_dengrad_x_b[a * npoints + g];
            As_y_b[threadIdx.y][threadIdx.x] = d_dengrad_y_b[a * npoints + g];
            As_z_b[threadIdx.y][threadIdx.x] = d_dengrad_z_b[a * npoints + g];

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 0];
                Bs_b[threadIdx.y] = d_grid_w[g + grid_offset] * d_vrho[2 * g + 1];
            }
        }
        else
        {
            As_x_a[threadIdx.y][threadIdx.x] = 0.0;
            As_y_a[threadIdx.y][threadIdx.x] = 0.0;
            As_z_a[threadIdx.y][threadIdx.x] = 0.0;

            As_x_b[threadIdx.y][threadIdx.x] = 0.0;
            As_y_b[threadIdx.y][threadIdx.x] = 0.0;
            As_z_b[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = 0.0;
                Bs_b[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x_a[m][threadIdx.x] * Bs_a[m];
                gatmy += As_y_a[m][threadIdx.x] * Bs_a[m];
                gatmz += As_z_a[m][threadIdx.x] * Bs_a[m];

                gatmx += As_x_b[m][threadIdx.x] * Bs_b[m];
                gatmy += As_y_b[m][threadIdx.x] * Bs_b[m];
                gatmz += As_z_b[m][threadIdx.x] * Bs_b[m];
            }
        }

        __syncthreads();

        // vx contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x_a[threadIdx.y][threadIdx.x] = d_dengrad_xx_a[a * npoints + g];
            As_y_a[threadIdx.y][threadIdx.x] = d_dengrad_xy_a[a * npoints + g];
            As_z_a[threadIdx.y][threadIdx.x] = d_dengrad_xz_a[a * npoints + g];

            As_x_b[threadIdx.y][threadIdx.x] = d_dengrad_xx_b[a * npoints + g];
            As_y_b[threadIdx.y][threadIdx.x] = d_dengrad_xy_b[a * npoints + g];
            As_z_b[threadIdx.y][threadIdx.x] = d_dengrad_xz_b[a * npoints + g];

            if (threadIdx.x == 0)
            {
                const double vxa = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 0] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 3];
                const double vxb = 2.0 * d_vsigma[3 * g + 2] * d_rhograd[6 * g + 3] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 0];

                Bs_a[threadIdx.y] = d_grid_w[g + grid_offset] * vxa;
                Bs_b[threadIdx.y] = d_grid_w[g + grid_offset] * vxb;
            }
        }
        else
        {
            As_x_a[threadIdx.y][threadIdx.x] = 0.0;
            As_y_a[threadIdx.y][threadIdx.x] = 0.0;
            As_z_a[threadIdx.y][threadIdx.x] = 0.0;

            As_x_b[threadIdx.y][threadIdx.x] = 0.0;
            As_y_b[threadIdx.y][threadIdx.x] = 0.0;
            As_z_b[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = 0.0;
                Bs_b[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if (threadIdx.y == 0)
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x_a[m][threadIdx.x] * Bs_a[m];
                gatmy += As_y_a[m][threadIdx.x] * Bs_a[m];
                gatmz += As_z_a[m][threadIdx.x] * Bs_a[m];

                gatmx += As_x_b[m][threadIdx.x] * Bs_b[m];
                gatmy += As_y_b[m][threadIdx.x] * Bs_b[m];
                gatmz += As_z_b[m][threadIdx.x] * Bs_b[m];
            }
        }

        __syncthreads();

        // vy contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x_a[threadIdx.y][threadIdx.x] = d_dengrad_yx_a[a * npoints + g];
            As_y_a[threadIdx.y][threadIdx.x] = d_dengrad_yy_a[a * npoints + g];
            As_z_a[threadIdx.y][threadIdx.x] = d_dengrad_yz_a[a * npoints + g];

            As_x_b[threadIdx.y][threadIdx.x] = d_dengrad_yx_b[a * npoints + g];
            As_y_b[threadIdx.y][threadIdx.x] = d_dengrad_yy_b[a * npoints + g];
            As_z_b[threadIdx.y][threadIdx.x] = d_dengrad_yz_b[a * npoints + g];

            if (threadIdx.x == 0)
            {
                const double vya = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 1] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 4];
                const double vyb = 2.0 * d_vsigma[3 * g + 2] * d_rhograd[6 * g + 4] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 1];

                Bs_a[threadIdx.y] = d_grid_w[g + grid_offset] * vya;
                Bs_b[threadIdx.y] = d_grid_w[g + grid_offset] * vyb;
            }
        }
        else
        {
            As_x_a[threadIdx.y][threadIdx.x] = 0.0;
            As_y_a[threadIdx.y][threadIdx.x] = 0.0;
            As_z_a[threadIdx.y][threadIdx.x] = 0.0;

            As_x_b[threadIdx.y][threadIdx.x] = 0.0;
            As_y_b[threadIdx.y][threadIdx.x] = 0.0;
            As_z_b[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = 0.0;
                Bs_b[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x_a[m][threadIdx.x] * Bs_a[m];
                gatmy += As_y_a[m][threadIdx.x] * Bs_a[m];
                gatmz += As_z_a[m][threadIdx.x] * Bs_a[m];

                gatmx += As_x_b[m][threadIdx.x] * Bs_b[m];
                gatmy += As_y_b[m][threadIdx.x] * Bs_b[m];
                gatmz += As_z_b[m][threadIdx.x] * Bs_b[m];
            }
        }

        __syncthreads();

        // vz contrib.

        if ((g < npoints) && (a < natoms))
        {
            As_x_a[threadIdx.y][threadIdx.x] = d_dengrad_zx_a[a * npoints + g];
            As_y_a[threadIdx.y][threadIdx.x] = d_dengrad_zy_a[a * npoints + g];
            As_z_a[threadIdx.y][threadIdx.x] = d_dengrad_zz_a[a * npoints + g];

            As_x_b[threadIdx.y][threadIdx.x] = d_dengrad_zx_b[a * npoints + g];
            As_y_b[threadIdx.y][threadIdx.x] = d_dengrad_zy_b[a * npoints + g];
            As_z_b[threadIdx.y][threadIdx.x] = d_dengrad_zz_b[a * npoints + g];

            if (threadIdx.x == 0)
            {
                const double vza = 2.0 * d_vsigma[3 * g + 0] * d_rhograd[6 * g + 2] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 5];
                const double vzb = 2.0 * d_vsigma[3 * g + 2] * d_rhograd[6 * g + 5] + d_vsigma[3 * g + 1] * d_rhograd[6 * g + 2];

                Bs_a[threadIdx.y] = d_grid_w[g + grid_offset] * vza;
                Bs_b[threadIdx.y] = d_grid_w[g + grid_offset] * vzb;
            }
        }
        else
        {
            As_x_a[threadIdx.y][threadIdx.x] = 0.0;
            As_y_a[threadIdx.y][threadIdx.x] = 0.0;
            As_z_a[threadIdx.y][threadIdx.x] = 0.0;

            As_x_b[threadIdx.y][threadIdx.x] = 0.0;
            As_y_b[threadIdx.y][threadIdx.x] = 0.0;
            As_z_b[threadIdx.y][threadIdx.x] = 0.0;

            if (threadIdx.x == 0)
            {
                Bs_a[threadIdx.y] = 0.0;
                Bs_b[threadIdx.y] = 0.0;
            }
        }

        __syncthreads();

        if ((threadIdx.y == 0) && (g < npoints))
        {
            for (uint32_t m = 0; m < TILE_DIM; m++)
            {
                gatmx += As_x_a[m][threadIdx.x] * Bs_a[m];
                gatmy += As_y_a[m][threadIdx.x] * Bs_a[m];
                gatmz += As_z_a[m][threadIdx.x] * Bs_a[m];

                gatmx += As_x_b[m][threadIdx.x] * Bs_b[m];
                gatmy += As_y_b[m][threadIdx.x] * Bs_b[m];
                gatmz += As_z_b[m][threadIdx.x] * Bs_b[m];
            }
        }

        __syncthreads();
    }

    // accumulate gradient

    if ((threadIdx.y == 0) && (a < natoms))
    {
        d_mol_grad[a * 3 + 0] += gatmx;
        d_mol_grad[a * 3 + 1] += gatmy;
        d_mol_grad[a * 3 + 2] += gatmz;
    }
}

__global__ void
matmulAB(double* C, double* A, double* B, uint32_t aocount, uint32_t npoints)
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
            As[threadIdx.y][threadIdx.x] = A[i * aocount + (k * TILE_DIM + threadIdx.x)];
        }
        else
        {
            As[threadIdx.y][threadIdx.x] = 0.0;
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
            val += As[threadIdx.y][m] * Bs[m][threadIdx.x];
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
            // Note: transpose of B
            Bs[threadIdx.x][threadIdx.y] = B[j * npoints + (k * TILE_DIM + threadIdx.y)];
        }
        else
        {
            Bs[threadIdx.x][threadIdx.y] = 0.0;
        }

        __syncthreads();

        for (uint32_t m = 0; m < TILE_DIM; m++)
        {
            // Note: transpose of B
            val += As[threadIdx.y][m] * Bs[threadIdx.x][m];
        }

        __syncthreads();
    }

    if ((i < aocount) && (j < aocount))
    {
        d_mat_Vxc_full[d_ao_inds[i] * naos + d_ao_inds[j]] += val;
    }
}

__global__ void
matmulABt(double* C, double* A, double* B, uint32_t aocount, uint32_t npoints)
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
            // Note: transpose of B
            Bs[threadIdx.x][threadIdx.y] = B[j * npoints + (k * TILE_DIM + threadIdx.y)];
        }
        else
        {
            Bs[threadIdx.x][threadIdx.y] = 0.0;
        }

        __syncthreads();

        for (uint32_t m = 0; m < TILE_DIM; m++)
        {
            // Note: transpose of B
            val += As[threadIdx.y][m] * Bs[threadIdx.x][m];
        }

        __syncthreads();
    }

    if ((i < aocount) && (j < aocount))
    {
        C[i * aocount + j] = val;
    }
}

__global__ void
distributeGgaSubmMatrix(double*         d_mat_full,
                        const uint32_t  naos,
                        const double*   d_mat_Fxc,
                        const double*   d_mat_Fxc_gga,
                        const uint32_t* d_ao_inds,
                        const uint32_t  aocount)
{
    const uint32_t i = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t j = blockIdx.x * blockDim.x + threadIdx.x;

    if ((i < aocount) && (j < aocount))
    {
        d_mat_full[d_ao_inds[i] * naos + d_ao_inds[j]] += (d_mat_Fxc[i * aocount + j] +
                                                           d_mat_Fxc_gga[i * aocount + j] +
                                                           d_mat_Fxc_gga[j * aocount + i]);
    }
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

    auto gto_info = gtoinfo::getGtoInfo(gto_block, gtos_mask);

    gpu::chunkedMemcpyHostToDevice<double>(d_gto_info, gto_info.data(), gto_info.size());

    // evaluate GTO values on grid points

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((ncols + threads_per_block.x - 1) / threads_per_block.x, (nrows + threads_per_block.y - 1) / threads_per_block.y);

    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        gpu::gtoValuesLdaRecS<<<num_blocks, threads_per_block>>>(d_gto_values,
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
        gpu::gtoValuesLdaRecP<<<num_blocks, threads_per_block>>>(d_gto_values,
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

        gpu::gtoValuesLdaRecD<<<num_blocks, threads_per_block>>>(d_gto_values,
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

    auto gto_info = gtoinfo::getGtoInfo(gto_block, gtos_mask);

    gpu::chunkedMemcpyHostToDevice<double>(d_gto_info, gto_info.data(), gto_info.size());

    // evaluate GTO values on grid points

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((ncols + threads_per_block.x - 1) / threads_per_block.x, (nrows + threads_per_block.y - 1) / threads_per_block.y);

    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        gpu::gtoValuesGgaRecS<<<num_blocks, threads_per_block>>>(d_gto_values_0,
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
        gpu::gtoValuesGgaRecP<<<num_blocks, threads_per_block>>>(d_gto_values_0,
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

        gpu::gtoValuesGgaRecD<<<num_blocks, threads_per_block>>>(d_gto_values_0,
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

static auto
getGtoValuesForMgga(double*                     d_gto_values_0,
                    double*                     d_gto_values_x,
                    double*                     d_gto_values_y,
                    double*                     d_gto_values_z,
                    double*                     d_gto_values_xx,
                    double*                     d_gto_values_xy,
                    double*                     d_gto_values_xz,
                    double*                     d_gto_values_yy,
                    double*                     d_gto_values_yz,
                    double*                     d_gto_values_zz,
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

    auto gto_info = gtoinfo::getGtoInfo(gto_block, gtos_mask);

    gpu::chunkedMemcpyHostToDevice<double>(d_gto_info, gto_info.data(), gto_info.size());

    // evaluate GTO values on grid points

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((ncols + threads_per_block.x - 1) / threads_per_block.x, (nrows + threads_per_block.y - 1) / threads_per_block.y);

    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        gpu::gtoValuesMggaRecS<<<num_blocks, threads_per_block>>>(
                                                                 d_gto_values_0,
                                                                 d_gto_values_x,
                                                                 d_gto_values_y,
                                                                 d_gto_values_z,
                                                                 d_gto_values_xx,
                                                                 d_gto_values_xy,
                                                                 d_gto_values_xz,
                                                                 d_gto_values_yy,
                                                                 d_gto_values_yz,
                                                                 d_gto_values_zz,
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
        gpu::gtoValuesMggaRecP<<<num_blocks, threads_per_block>>>(
                                                                 d_gto_values_0,
                                                                 d_gto_values_x,
                                                                 d_gto_values_y,
                                                                 d_gto_values_z,
                                                                 d_gto_values_xx,
                                                                 d_gto_values_xy,
                                                                 d_gto_values_xz,
                                                                 d_gto_values_yy,
                                                                 d_gto_values_yz,
                                                                 d_gto_values_zz,
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

        gpu::gtoValuesMggaRecD<<<num_blocks, threads_per_block>>>(
                                                                 d_gto_values_0,
                                                                 d_gto_values_x,
                                                                 d_gto_values_y,
                                                                 d_gto_values_z,
                                                                 d_gto_values_xx,
                                                                 d_gto_values_xy,
                                                                 d_gto_values_xz,
                                                                 d_gto_values_yy,
                                                                 d_gto_values_yz,
                                                                 d_gto_values_zz,
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
        std::string err_ang("gpu::getGtoValuesForMgga: Only implemented for s, p and d-orbitals");

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

    gpuSafe(gpuMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double* d_gaos;

    gpuSafe(gpuMalloc(&d_gaos, naos * max_npoints_per_box * sizeof(double)));

    // coordinates of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z;

    gpuSafe(gpuMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    gpuSafe(gpuMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    gpuSafe(gpuMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_grid_x, xcoords, n_total_grid_points);
    gpu::chunkedMemcpyHostToDevice<double>(d_grid_y, ycoords, n_total_grid_points);
    gpu::chunkedMemcpyHostToDevice<double>(d_grid_z, zcoords, n_total_grid_points);

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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_chi.values(), d_gaos, aocount * npoints);

        for (int64_t nu = 0; nu < aocount; nu++)
        {
            std::memcpy(allgtovalues.row(aoinds[nu]) + gridblockpos, mat_chi.row(nu), npoints * sizeof(double));
        }
    }

    gpuSafe(gpuFree(d_gto_info));
    gpuSafe(gpuFree(d_gaos));
    gpuSafe(gpuFree(d_grid_x));
    gpuSafe(gpuFree(d_grid_y));
    gpuSafe(gpuFree(d_grid_z));

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

    gpuSafe(gpuMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // GTO values on grid points

    CDenseMatrix allgtovalues_0(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_x(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_y(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_z(naos, molecularGrid.getNumberOfGridPoints());

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    double *d_gaos, *d_gaox, *d_gaoy, *d_gaoz;

    gpuSafe(gpuMalloc(&d_gaos, naos * max_npoints_per_box * sizeof(double)));
    gpuSafe(gpuMalloc(&d_gaox, naos * max_npoints_per_box * sizeof(double)));
    gpuSafe(gpuMalloc(&d_gaoy, naos * max_npoints_per_box * sizeof(double)));
    gpuSafe(gpuMalloc(&d_gaoz, naos * max_npoints_per_box * sizeof(double)));

    // coordinates of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z;

    gpuSafe(gpuMalloc(&d_grid_x, n_total_grid_points * sizeof(double)));
    gpuSafe(gpuMalloc(&d_grid_y, n_total_grid_points * sizeof(double)));
    gpuSafe(gpuMalloc(&d_grid_z, n_total_grid_points * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_grid_x, xcoords, n_total_grid_points);
    gpu::chunkedMemcpyHostToDevice<double>(d_grid_y, ycoords, n_total_grid_points);
    gpu::chunkedMemcpyHostToDevice<double>(d_grid_z, zcoords, n_total_grid_points);

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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_chi_0.values(), d_gaos, aocount * npoints);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_chi_x.values(), d_gaox, aocount * npoints);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_chi_y.values(), d_gaoy, aocount * npoints);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_chi_z.values(), d_gaoz, aocount * npoints);

        for (int64_t nu = 0; nu < aocount; nu++)
        {
            std::memcpy(allgtovalues_0.row(aoinds[nu]) + gridblockpos, mat_chi_0.row(nu), npoints * sizeof(double));
            std::memcpy(allgtovalues_x.row(aoinds[nu]) + gridblockpos, mat_chi_x.row(nu), npoints * sizeof(double));
            std::memcpy(allgtovalues_y.row(aoinds[nu]) + gridblockpos, mat_chi_y.row(nu), npoints * sizeof(double));
            std::memcpy(allgtovalues_z.row(aoinds[nu]) + gridblockpos, mat_chi_z.row(nu), npoints * sizeof(double));
        }
    }

    gpuSafe(gpuFree(d_gto_info));
    gpuSafe(gpuFree(d_gaos));
    gpuSafe(gpuFree(d_gaox));
    gpuSafe(gpuFree(d_gaoy));
    gpuSafe(gpuFree(d_gaoz));
    gpuSafe(gpuFree(d_grid_x));
    gpuSafe(gpuFree(d_grid_y));
    gpuSafe(gpuFree(d_grid_z));

    return std::vector<CDenseMatrix>({allgtovalues_0, allgtovalues_x, allgtovalues_y, allgtovalues_z});
}

static auto
integrateVxcFockForLdaClosedShell(const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const CAODensityMatrix& densityMatrix,
                                  const CMolecularGrid&   molecularGrid,
                                  const CXCFunctional&    xcFunctional,
                                  const int64_t           numGpusPerNode,
                                  const int64_t           rank,
                                  const int64_t           nnodes) -> CAOKohnShamMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    errors::assertMsgCritical(
        (naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)),
        std::string(__func__) + std::string(": Inconsistent number of AOs"));

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // Kohn-Sham matrix

    auto is_closed_shell = densityMatrix.isClosedShell();

    CAOKohnShamMatrix mat_Vxc_sum(naos, naos, is_closed_shell);

    mat_Vxc_sum.zero();

    std::vector<CAOKohnShamMatrix> mat_Vxc_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Vxc_omp[gpu_id] = CAOKohnShamMatrix(naos, naos, is_closed_shell);
    }

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    // Kohn-Sham matrix

    mat_Vxc_omp[gpu_id].zero();

    // GTOs on grid points

    double *d_mat_Vxc_full, *d_den_mat_full, *d_den_mat, *d_gto_values, *d_mat_F;

    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Vxc_full, naos * naos * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_den_mat_full, densityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ldafunc = xcfun_copy.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    double *d_rho, *d_exc, *d_vrho;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_exc, dim->zk * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Vxc_full, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // LDA density for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensityOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_mat_F,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // functional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        gpuSafe(gpuMemcpyAsync(d_exc, exc, dim->zk * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // reusing d_mat_F
        auto d_mat_G = d_mat_F;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getLdaVxcMatrixG<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_G,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           d_vrho);

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Vxc_full,
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

    gpuSafe(gpuMemcpyAsync(mat_Vxc_omp[gpu_id].alphaValues(), d_mat_Vxc_full, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Vxc_full, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_exc, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());

    mat_Vxc_omp[gpu_id].setNumberOfElectrons(nele);
    mat_Vxc_omp[gpu_id].setExchangeCorrelationEnergy(xcene);
    }

    auto p_mat_Vxc = mat_Vxc_sum.alphaValues();

    double nele_sum = 0.0, xcene_sum = 0.0;

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_v = mat_Vxc_omp[gpu_id].alphaValues();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_Vxc[ind] += p_mat_v[ind];
        }

        nele_sum += mat_Vxc_omp[gpu_id].getNumberOfElectrons();

        xcene_sum += mat_Vxc_omp[gpu_id].getExchangeCorrelationEnergy();
    }

    mat_Vxc_sum.setNumberOfElectrons(nele_sum);

    mat_Vxc_sum.setExchangeCorrelationEnergy(xcene_sum);

    return mat_Vxc_sum;
}

static auto
integrateVxcFockForLdaOpenShell(const CMolecule&        molecule,
                                const CMolecularBasis&  basis,
                                const CAODensityMatrix& densityMatrix,
                                const CMolecularGrid&   molecularGrid,
                                const CXCFunctional&    xcFunctional,
                                const int64_t           numGpusPerNode,
                                const int64_t           rank,
                                const int64_t           nnodes) -> CAOKohnShamMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    errors::assertMsgCritical(
        (naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)),
        std::string(__func__) + std::string(": Inconsistent number of AOs"));

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // Kohn-Sham matrix

    auto is_closed_shell = densityMatrix.isClosedShell();

    CAOKohnShamMatrix mat_Vxc_sum(naos, naos, is_closed_shell);

    mat_Vxc_sum.zero();

    std::vector<CAOKohnShamMatrix> mat_Vxc_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Vxc_omp[gpu_id] = CAOKohnShamMatrix(naos, naos, is_closed_shell);
    }

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    // Kohn-Sham matrix

    mat_Vxc_omp[gpu_id].zero();

    // GTOs on grid points

    double *d_mat_Vxc_full, *d_den_mat_full, *d_den_mat, *d_gto_values, *d_mat_F;

    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat_full, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Vxc_full, naos * naos * 2 * sizeof(double), stream));

    double* d_den_mat_a = d_den_mat;
    double* d_den_mat_b = d_den_mat_a + naos * naos;

    double* d_den_mat_full_a = d_den_mat_full;
    double* d_den_mat_full_b = d_den_mat_full_a + naos * naos;

    double* d_mat_F_a = d_mat_F;
    double* d_mat_F_b = d_mat_F_a + naos * max_npoints_per_box;

    double* d_mat_Vxc_full_a = d_mat_Vxc_full;
    double* d_mat_Vxc_full_b = d_mat_Vxc_full_a + naos * naos;

    gpuSafe(gpuMemcpyAsync(d_den_mat_full_a, densityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_den_mat_full_b, densityMatrix.betaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ldafunc = xcfun_copy.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    double *d_rho, *d_exc, *d_vrho;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_exc, dim->zk * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Vxc_full, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_a,
                           d_den_mat_full_a,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_b,
                           d_den_mat_full_b,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // LDA density for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_a,
                           d_den_mat_a,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_b,
                           d_den_mat_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensityOnGridsOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_mat_F_a,
                           d_mat_F_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // functional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        gpuSafe(gpuMemcpyAsync(d_exc, exc, dim->zk * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // reusing d_mat_F
        auto d_mat_G_a = d_mat_F_a;
        auto d_mat_G_b = d_mat_F_b;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getLdaVxcMatrixGOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_G_a,
                           d_mat_G_b,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           d_vrho);

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Vxc_full_a,
                           static_cast<uint32_t>(naos),
                           d_gto_values,
                           d_mat_G_a,
                           d_ao_inds,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Vxc_full_b,
                           static_cast<uint32_t>(naos),
                           d_gto_values,
                           d_mat_G_b,
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

    gpuSafe(gpuMemcpyAsync(mat_Vxc_omp[gpu_id].alphaValues(), d_mat_Vxc_full_a, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));
    gpuSafe(gpuMemcpyAsync(mat_Vxc_omp[gpu_id].betaValues(), d_mat_Vxc_full_b, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Vxc_full, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_exc, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());

    mat_Vxc_omp[gpu_id].setNumberOfElectrons(nele);
    mat_Vxc_omp[gpu_id].setExchangeCorrelationEnergy(xcene);
    }

    auto p_mat_Vxc_a = mat_Vxc_sum.alphaValues();
    auto p_mat_Vxc_b = mat_Vxc_sum.betaValues();

    double nele_sum = 0.0, xcene_sum = 0.0;

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_v_a = mat_Vxc_omp[gpu_id].alphaValues();
        auto p_mat_v_b = mat_Vxc_omp[gpu_id].betaValues();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_Vxc_a[ind] += p_mat_v_a[ind];
            p_mat_Vxc_b[ind] += p_mat_v_b[ind];
        }

        nele_sum += mat_Vxc_omp[gpu_id].getNumberOfElectrons();

        xcene_sum += mat_Vxc_omp[gpu_id].getExchangeCorrelationEnergy();
    }

    mat_Vxc_sum.setNumberOfElectrons(nele_sum);

    mat_Vxc_sum.setExchangeCorrelationEnergy(xcene_sum);

    return mat_Vxc_sum;
}

static auto
integrateVxcFockForGgaClosedShell(const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const CAODensityMatrix& densityMatrix,
                                  const CMolecularGrid&   molecularGrid,
                                  const CXCFunctional&    xcFunctional,
                                  const int64_t           numGpusPerNode,
                                  const int64_t           rank,
                                  const int64_t           nnodes) -> CAOKohnShamMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    errors::assertMsgCritical(
        (naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)),
        std::string(__func__) + std::string(": Inconsistent number of AOs"));

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // Kohn-Sham matrix

    auto is_closed_shell = densityMatrix.isClosedShell();

    CAOKohnShamMatrix mat_Vxc_sum(naos, naos, is_closed_shell);

    mat_Vxc_sum.zero();

    std::vector<CAOKohnShamMatrix> mat_Vxc_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Vxc_omp[gpu_id] = CAOKohnShamMatrix(naos, naos, is_closed_shell);
    }

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    // Kohn-Sham matrix

    mat_Vxc_omp[gpu_id].zero();

    // GTOs on grid points

    double *d_mat_Vxc_full, *d_mat_Vxc, *d_den_mat_full, *d_den_mat, *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z, *d_mat_F;

    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Vxc_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Vxc, naos * naos * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_den_mat_full, densityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ggafunc = xcfun_copy.getFunctionalPointerToGgaComponent();
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

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhograd, dim->rho * 3 * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_sigma, dim->sigma * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_exc, dim->zk * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vsigma, dim->vsigma * max_npoints_per_box * sizeof(double), stream));

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Vxc_full, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // GGA density for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensitySigmaOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_rhograd,
                           d_sigma,
                           d_mat_F,
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // funtional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(rhograd, d_rhograd, dim->rho * 3 * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(sigma, d_sigma, dim->sigma * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        gpuSafe(gpuMemcpyAsync(d_exc, exc, dim->zk * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vsigma, vsigma, dim->vsigma * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // reusing d_mat_F
        auto d_mat_G = d_mat_F;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getGgaVxcMatrixG<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_G,
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

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Vxc_full,
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

    gpuSafe(gpuMemcpyAsync(mat_Vxc_omp[gpu_id].alphaValues(), d_mat_Vxc_full, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    // Note: symmetrize mat_Vxc

    mat_Vxc_omp[gpu_id].symmetrizeAndScale(0.5);

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Vxc_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Vxc, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_rhograd, stream));
    gpuSafe(gpuFreeAsync(d_sigma, stream));
    gpuSafe(gpuFreeAsync(d_exc, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));
    gpuSafe(gpuFreeAsync(d_vsigma, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());

    mat_Vxc_omp[gpu_id].setNumberOfElectrons(nele);

    mat_Vxc_omp[gpu_id].setExchangeCorrelationEnergy(xcene);
    }

    auto p_mat_Vxc = mat_Vxc_sum.alphaValues();

    double nele_sum = 0.0, xcene_sum = 0.0;

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_v = mat_Vxc_omp[gpu_id].alphaValues();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_Vxc[ind] += p_mat_v[ind];
        }

        nele_sum += mat_Vxc_omp[gpu_id].getNumberOfElectrons();

        xcene_sum += mat_Vxc_omp[gpu_id].getExchangeCorrelationEnergy();
    }

    mat_Vxc_sum.setNumberOfElectrons(nele_sum);

    mat_Vxc_sum.setExchangeCorrelationEnergy(xcene_sum);

    return mat_Vxc_sum;
}

static auto
integrateVxcFockForGgaOpenShell(const CMolecule&        molecule,
                                const CMolecularBasis&  basis,
                                const CAODensityMatrix& densityMatrix,
                                const CMolecularGrid&   molecularGrid,
                                const CXCFunctional&    xcFunctional,
                                const int64_t           numGpusPerNode,
                                const int64_t           rank,
                                const int64_t           nnodes) -> CAOKohnShamMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    errors::assertMsgCritical(
        (naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)),
        std::string(__func__) + std::string(": Inconsistent number of AOs"));

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // Kohn-Sham matrix

    auto is_closed_shell = densityMatrix.isClosedShell();

    CAOKohnShamMatrix mat_Vxc_sum(naos, naos, is_closed_shell);

    mat_Vxc_sum.zero();

    std::vector<CAOKohnShamMatrix> mat_Vxc_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Vxc_omp[gpu_id] = CAOKohnShamMatrix(naos, naos, is_closed_shell);
    }

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    // Kohn-Sham matrix

    mat_Vxc_omp[gpu_id].zero();

    // GTOs on grid points

    double *d_mat_Vxc_full, *d_mat_Vxc, *d_den_mat_full, *d_den_mat, *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z, *d_mat_F;

    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Vxc, naos * naos * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat_full, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Vxc_full, naos * naos * 2 * sizeof(double), stream));

    double* d_den_mat_a = d_den_mat;
    double* d_den_mat_b = d_den_mat_a + naos * naos;

    double* d_den_mat_full_a = d_den_mat_full;
    double* d_den_mat_full_b = d_den_mat_full_a + naos * naos;

    double* d_mat_F_a = d_mat_F;
    double* d_mat_F_b = d_mat_F_a + naos * max_npoints_per_box;

    double* d_mat_Vxc_full_a = d_mat_Vxc_full;
    double* d_mat_Vxc_full_b = d_mat_Vxc_full_a + naos * naos;

    gpuSafe(gpuMemcpyAsync(d_den_mat_full_a, densityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_den_mat_full_b, densityMatrix.betaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ggafunc = xcfun_copy.getFunctionalPointerToGgaComponent();
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

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhograd, dim->rho * 3 * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_sigma, dim->sigma * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_exc, dim->zk * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vsigma, dim->vsigma * max_npoints_per_box * sizeof(double), stream));

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Vxc_full_a, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));
    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Vxc_full_b, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_a,
                           d_den_mat_full_a,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_b,
                           d_den_mat_full_b,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // GGA density for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_a,
                           d_den_mat_a,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_b,
                           d_den_mat_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensitySigmaOnGridsOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_rhograd,
                           d_sigma,
                           d_mat_F_a,
                           d_mat_F_b,
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // funtional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(rhograd, d_rhograd, dim->rho * 3 * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(sigma, d_sigma, dim->sigma * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        gpuSafe(gpuMemcpyAsync(d_exc, exc, dim->zk * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vsigma, vsigma, dim->vsigma * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // reusing d_mat_F
        auto d_mat_G_a = d_mat_F_a;
        auto d_mat_G_b = d_mat_F_b;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getGgaVxcMatrixGOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_G_a,
                           d_mat_G_b,
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

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Vxc_full_a,
                           static_cast<uint32_t>(naos),
                           d_gto_values,
                           d_mat_G_a,
                           d_ao_inds,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Vxc_full_b,
                           static_cast<uint32_t>(naos),
                           d_gto_values,
                           d_mat_G_b,
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

    gpuSafe(gpuMemcpyAsync(mat_Vxc_omp[gpu_id].alphaValues(), d_mat_Vxc_full_a, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));
    gpuSafe(gpuMemcpyAsync(mat_Vxc_omp[gpu_id].betaValues(),  d_mat_Vxc_full_b, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    // Note: symmetrize mat_Vxc

    mat_Vxc_omp[gpu_id].symmetrizeAndScale(0.5);

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Vxc_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Vxc, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_rhograd, stream));
    gpuSafe(gpuFreeAsync(d_sigma, stream));
    gpuSafe(gpuFreeAsync(d_exc, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));
    gpuSafe(gpuFreeAsync(d_vsigma, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());

    mat_Vxc_omp[gpu_id].setNumberOfElectrons(nele);
    mat_Vxc_omp[gpu_id].setExchangeCorrelationEnergy(xcene);
    }

    auto p_mat_Vxc_a = mat_Vxc_sum.alphaValues();
    auto p_mat_Vxc_b = mat_Vxc_sum.betaValues();

    double nele_sum = 0.0, xcene_sum = 0.0;

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_v_a = mat_Vxc_omp[gpu_id].alphaValues();
        auto p_mat_v_b = mat_Vxc_omp[gpu_id].betaValues();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_Vxc_a[ind] += p_mat_v_a[ind];
            p_mat_Vxc_b[ind] += p_mat_v_b[ind];
        }

        nele_sum += mat_Vxc_omp[gpu_id].getNumberOfElectrons();

        xcene_sum += mat_Vxc_omp[gpu_id].getExchangeCorrelationEnergy();
    }

    mat_Vxc_sum.setNumberOfElectrons(nele_sum);

    mat_Vxc_sum.setExchangeCorrelationEnergy(xcene_sum);

    return mat_Vxc_sum;
}

auto
integrateVxcFock(const CMolecule&        molecule,
                 const CMolecularBasis&  basis,
                 const CAODensityMatrix& densityMatrix,
                 const CMolecularGrid&   molecularGrid,
                 const std::string&      xcFuncLabel,
                 const int64_t           numGpusPerNode,
                 const int64_t           rank,
                 const int64_t           nnodes) -> CAOKohnShamMatrix
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    auto flag = densityMatrix.isClosedShell() ? std::string("CLOSEDSHELL") : std::string("OPENSHELL");

    if (xcfuntype == xcfun::lda)
    {
        if (flag == std::string("CLOSEDSHELL"))
        {
            return gpu::integrateVxcFockForLdaClosedShell(molecule, basis, densityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
        else
        {
            return gpu::integrateVxcFockForLdaOpenShell(molecule, basis, densityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
    }
    else if (xcfuntype == xcfun::gga)
    {
        if (flag == std::string("CLOSEDSHELL"))
        {
            return gpu::integrateVxcFockForGgaClosedShell(molecule, basis, densityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
        else
        {
            return gpu::integrateVxcFockForGgaOpenShell(molecule, basis, densityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
    }
    else
    {
        std::string errxcfuntype("gpu::integrateVxcFock: Only implemented for LDA/GGA");

        errors::assertMsgCritical(false, errxcfuntype);

        return CAOKohnShamMatrix();
    }
}

static auto
integrateFxcFockForLDA(CDenseMatrix&           aoFockMatrix,
                       const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CAODensityMatrix& rwDensityMatrix,
                       const CAODensityMatrix& gsDensityMatrix,
                       const CMolecularGrid&   molecularGrid,
                       const CXCFunctional&    xcFunctional,
                       const int64_t           numGpusPerNode,
                       const int64_t           rank,
                       const int64_t           nnodes) -> void
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("gpu::integrateFxcFockForLDA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == gsDensityMatrix.getNumberOfRows(0)) && (naos == gsDensityMatrix.getNumberOfColumns(0)), errnaos);
    errors::assertMsgCritical((naos == rwDensityMatrix.getNumberOfRows(0)) && (naos == rwDensityMatrix.getNumberOfColumns(0)), errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // Fxc contribution to Fock matrix

    std::vector<CDenseMatrix> mat_Fxc_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Fxc_omp[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    // Fxc matrix

    mat_Fxc_omp[gpu_id].zero();

    // GTOs on grid points

    double *d_mat_Fxc_full, *d_gs_den_mat_full, *d_rw_den_mat_full, *d_den_mat, *d_gto_values, *d_mat_F;

    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gs_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rw_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Fxc_full, naos * naos * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full, gsDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full, rwDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ldafunc = xcfun_copy.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhow_data(dim->rho * max_npoints_per_box);

    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);

    auto rho = rho_data.data();
    auto rhow = rhow_data.data();

    auto v2rho2 = v2rho2_data.data();

    double *d_rho, *d_rhow, *d_v2rho2;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhow, dim->rho * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_v2rho2, dim->v2rho2 * max_npoints_per_box * sizeof(double), stream));

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Fxc_full, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_gs_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // LDA density for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensityOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_mat_F,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // functional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_fxc_for_lda(npoints, rho, v2rho2);

        gpuSafe(gpuMemcpyAsync(d_v2rho2, v2rho2, dim->v2rho2 * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // compute partial contribution to Fxc matrix and distribute partial

        // sub density matrix for rho_w
        // reusing d_den_mat

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_rw_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // LDA density for rho_w
        // reusing d_mat_F

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensityOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rhow,
                           d_mat_F,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // integrate LDA Fxc contribution to Fock

        // reusing d_mat_F
        auto d_mat_G = d_mat_F;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getLdaFxcMatrixG<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_G,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           d_rhow,
                           static_cast<uint32_t>(dim->rho),
                           d_v2rho2,
                           static_cast<uint32_t>(dim->v2rho2));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulABtKohnShamMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Fxc_full,
                           static_cast<uint32_t>(naos),
                           d_gto_values,
                           d_mat_G,
                           d_ao_inds,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));
    }

    gpuSafe(gpuMemcpyAsync(mat_Fxc_omp[gpu_id].values(), d_mat_Fxc_full, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_gs_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_rw_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Fxc_full, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_rhow, stream));
    gpuSafe(gpuFreeAsync(d_v2rho2, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());
    }

    auto p_mat_Fxc = aoFockMatrix.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_f = mat_Fxc_omp[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_Fxc[ind] += p_mat_f[ind];
        }
    }
}

static auto
integrateFxcFockForGGA(CDenseMatrix&           aoFockMatrix,
                       const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CAODensityMatrix& rwDensityMatrix,
                       const CAODensityMatrix& gsDensityMatrix,
                       const CMolecularGrid&   molecularGrid,
                       const CXCFunctional&    xcFunctional,
                       const int64_t           numGpusPerNode,
                       const int64_t           rank,
                       const int64_t           nnodes) -> void
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("gpu::integrateFxcFockForGGA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == gsDensityMatrix.getNumberOfRows(0)) && (naos == gsDensityMatrix.getNumberOfColumns(0)), errnaos);
    errors::assertMsgCritical((naos == rwDensityMatrix.getNumberOfRows(0)) && (naos == rwDensityMatrix.getNumberOfColumns(0)), errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // Fxc contribution to Fock matrix

    std::vector<CDenseMatrix> mat_Fxc_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Fxc_omp[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    // Fxc matrix

    mat_Fxc_omp[gpu_id].zero();

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ggafunc = xcfun_copy.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> rhow_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhowgrad_data(dim->rho * 3 * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto rhow     = rhow_data.data();
    auto rhowgrad = rhowgrad_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2     = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2   = v2sigma2_data.data();

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    // data on device

    double* d_gto_info       ;
    double* d_den_mat        ;
    double* d_gto_values     ;
    double* d_gto_values_x   ;
    double* d_gto_values_y   ;
    double* d_gto_values_z   ;
    double* d_mat_F          ;
    double* d_mat_G_gga      ;
    double* d_mat_Fxc        ;
    double* d_mat_Fxc_gga    ;
    double* d_gs_den_mat_full;
    double* d_rw_den_mat_full;
    double* d_mat_Fxc_full   ;
    double* d_rho            ;
    double* d_rhograd        ;
    double* d_sigma          ;
    double* d_rhow           ;
    double* d_rhowgrad       ;
    double* d_vsigma         ;
    double* d_v2rho2         ;
    double* d_v2rhosigma     ;
    double* d_v2sigma2       ;
    double* d_grid_x         ;
    double* d_grid_y         ;
    double* d_grid_z         ;
    double* d_grid_w         ;

    gpuSafe(gpuMallocAsync(&d_gto_info,          5 * max_ncgtos * max_npgtos           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat,           naos * naos                           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values,        naos * max_npoints_per_box            * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x,      naos * max_npoints_per_box            * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y,      naos * max_npoints_per_box            * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z,      naos * max_npoints_per_box            * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F,             naos * max_npoints_per_box            * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_G_gga,         naos * max_npoints_per_box            * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Fxc,           naos * naos                           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Fxc_gga,       naos * naos                           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gs_den_mat_full,   naos * naos                           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rw_den_mat_full,   naos * naos                           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_Fxc_full,      naos * naos                           * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rho,               dim->rho * max_npoints_per_box        * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhograd,           dim->rho * 3 * max_npoints_per_box    * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_sigma,             dim->sigma * max_npoints_per_box      * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhow,              dim->rho * max_npoints_per_box        * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhowgrad,          dim->rho * 3 * max_npoints_per_box    * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vsigma,            dim->vsigma * max_npoints_per_box     * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_v2rho2,            dim->v2rho2 * max_npoints_per_box     * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_v2rhosigma,        dim->v2rhosigma * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_v2sigma2,          dim->v2sigma2 * max_npoints_per_box   * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_x,            n_total_grid_points                   * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y,            n_total_grid_points                   * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z,            n_total_grid_points                   * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w,            n_total_grid_points                   * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full, gsDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full, rwDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((naos + threads_per_block.x - 1) / threads_per_block.x, (naos + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mat_Fxc_full, static_cast<uint32_t>(naos), static_cast<uint32_t>(naos));

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        // Note: symmetrize density matrix for GGA 
        gpu::getSymmetricSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_gs_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // GGA density for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensitySigmaOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_rhograd,
                           d_sigma,
                           d_mat_F,
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // functional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(rhograd, d_rhograd, dim->rho * 3 * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(sigma, d_sigma, dim->sigma * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcfun_copy.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        gpuSafe(gpuMemcpyAsync(d_vsigma, vsigma, dim->vsigma * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        gpuSafe(gpuMemcpyAsync(d_v2rho2, v2rho2, dim->v2rho2 * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_v2rhosigma, v2rhosigma, dim->v2rhosigma * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_v2sigma2, v2sigma2, dim->v2sigma2 * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // compute partial contribution to Fxc matrix and distribute partial

        // sub density matrix for rho_w
        // reusing d_den_mat

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        // Note: symmetrize density matrix for GGA
        gpu::getSymmetricSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_rw_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // GGA density for rho_w
        // reusing d_mat_F

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensitySigmaOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rhow,
                           d_rhowgrad,
                           nullptr,
                           d_mat_F,
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // integrate GGA Fxc contribution to Fock

        // reusing d_mat_F
        auto d_mat_G = d_mat_F;

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getGgaFxcMatrixG<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_G,
                           d_mat_G_gga,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           d_rhow,
                           d_rhograd,
                           d_rhowgrad,
                           static_cast<uint32_t>(dim->rho),
                           d_vsigma,
                           static_cast<uint32_t>(dim->vsigma),
                           d_v2rho2,
                           static_cast<uint32_t>(dim->v2rho2),
                           d_v2rhosigma,
                           static_cast<uint32_t>(dim->v2rhosigma),
                           d_v2sigma2,
                           static_cast<uint32_t>(dim->v2sigma2));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulABt<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Fxc,
                           d_gto_values,
                           d_mat_G,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulABt<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Fxc_gga,
                           d_gto_values,
                           d_mat_G_gga,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::distributeGgaSubmMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_Fxc_full,
                           static_cast<uint32_t>(naos),
                           d_mat_Fxc,
                           d_mat_Fxc_gga,
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));
    }

    gpuSafe(gpuMemcpyAsync(mat_Fxc_omp[gpu_id].values(), d_mat_Fxc_full, naos * naos * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info       , stream));
    gpuSafe(gpuFreeAsync(d_den_mat        , stream));
    gpuSafe(gpuFreeAsync(d_gto_values     , stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x   , stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y   , stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z   , stream));
    gpuSafe(gpuFreeAsync(d_mat_F          , stream));
    gpuSafe(gpuFreeAsync(d_mat_G_gga      , stream));
    gpuSafe(gpuFreeAsync(d_mat_Fxc        , stream));
    gpuSafe(gpuFreeAsync(d_mat_Fxc_gga    , stream));
    gpuSafe(gpuFreeAsync(d_gs_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_rw_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_mat_Fxc_full   , stream));
    gpuSafe(gpuFreeAsync(d_rho            , stream));
    gpuSafe(gpuFreeAsync(d_rhograd        , stream));
    gpuSafe(gpuFreeAsync(d_sigma          , stream));
    gpuSafe(gpuFreeAsync(d_rhow           , stream));
    gpuSafe(gpuFreeAsync(d_rhowgrad       , stream));
    gpuSafe(gpuFreeAsync(d_vsigma         , stream));
    gpuSafe(gpuFreeAsync(d_v2rho2         , stream));
    gpuSafe(gpuFreeAsync(d_v2rhosigma     , stream));
    gpuSafe(gpuFreeAsync(d_v2sigma2       , stream));
    gpuSafe(gpuFreeAsync(d_grid_x         , stream));
    gpuSafe(gpuFreeAsync(d_grid_y         , stream));
    gpuSafe(gpuFreeAsync(d_grid_z         , stream));
    gpuSafe(gpuFreeAsync(d_grid_w         , stream));

    gpuSafe(gpuFreeAsync(d_ao_inds, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());
    }

    auto p_mat_Fxc = aoFockMatrix.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_f = mat_Fxc_omp[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_Fxc[ind] += p_mat_f[ind];
        }
    }
}

auto
integrateFxcFock(CDenseMatrix&           aoFockMatrix,
                 const CMolecule&        molecule,
                 const CMolecularBasis&  basis,
                 const CAODensityMatrix& rwDensityMatrix,
                 const CAODensityMatrix& gsDensityMatrix,
                 const CMolecularGrid&   molecularGrid,
                 const std::string&      xcFuncLabel,
                 const int64_t           numGpusPerNode,
                 const int64_t           rank,
                 const int64_t           nnodes) -> void
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    std::string erropenshell("gpu::integrateFxcFock: Only implemented for closed-shell");

    const auto rw_closed_shell = rwDensityMatrix.isClosedShell();
    const auto gs_closed_shell = gsDensityMatrix.isClosedShell();

    errors::assertMsgCritical((rw_closed_shell && gs_closed_shell), erropenshell);

    if (xcfuntype == xcfun::lda)
    {
        integrateFxcFockForLDA(aoFockMatrix, molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
    }
    else if (xcfuntype == xcfun::gga)
    {
        integrateFxcFockForGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
    }
    /*
    else if (xcfuntype == xcfun::mgga)
    {
        integrateFxcFockForMGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode);
    }
    */
    else
    {
        std::string errxcfuntype("gpu::integrateFxcFock: Only implemented for LDA/GGA");

        errors::assertMsgCritical(false, errxcfuntype);
    }
}

static auto
integrateVxcGradientForLdaClosedShell(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& rwDensityMatrix,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecularGrid&   molecularGrid,
                                      const CXCFunctional&    xcFunctional,
                                      const int64_t           numGpusPerNode,
                                      const int64_t           rank,
                                      const int64_t           nnodes) -> CDenseMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    const auto errnaos = std::string(__func__) + std::string(": Inconsistent number of AOs");

    errors::assertMsgCritical((naos == rwDensityMatrix.getNumberOfRows(0)) &&
                              (naos == rwDensityMatrix.getNumberOfColumns(0)) &&
                              (naos == gsDensityMatrix.getNumberOfRows(0)) &&
                              (naos == gsDensityMatrix.getNumberOfColumns(0)),
                              errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    // AO-to-atom mapping

    std::vector<int64_t> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix molgrad_omp(num_gpus_per_node, natoms * 3);

    molgrad_omp.zero();

    // grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    uint32_t* d_ao_to_atom_ids;

    gpuSafe(gpuMallocAsync(&d_ao_to_atom_ids, naos * sizeof(uint32_t), stream));

    // GTOs on grid points

    double *d_mol_grad, *d_rw_den_mat_full, *d_gs_den_mat_full, *d_den_mat, *d_mat_F;
    double *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z;
    double *d_dengrad_x, *d_dengrad_y, *d_dengrad_z;

    gpuSafe(gpuMallocAsync(&d_mol_grad, natoms * 3 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rw_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gs_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_x, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_y, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_z, natoms * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full, rwDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full, gsDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ldafunc = xcfun_copy.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    double *d_rho, *d_exc, *d_vrho;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_exc, dim->zk * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(1, TILE_DIM);

    dim3 num_blocks(1, (natoms * 3 + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mol_grad, static_cast<uint32_t>(natoms * 3), 1);

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices and AO-to-Atom indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        std::vector<uint32_t> ao_to_atom_ids_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);

            ao_to_atom_ids_int32[ind] = static_cast<uint32_t>(ao_to_atom_ids[aoinds[ind]]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        gpuSafe(gpuMemcpyAsync(d_ao_to_atom_ids, ao_to_atom_ids_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_gs_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // LDA density for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensityOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_mat_F,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // density gradient on grid points

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_rw_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (natoms + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getLdaDensityGradientOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_dengrad_x,
                           d_dengrad_y,
                           d_dengrad_z,
                           d_mat_F,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           d_ao_to_atom_ids,
                           static_cast<uint32_t>(natoms),
                           static_cast<uint32_t>(npoints));

        // functional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuDeviceSynchronize());

        xcfun_copy.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        gpuSafe(gpuMemcpyAsync(d_exc, exc, dim->zk * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // accumulate partial contribution to Vxc gradient

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((natoms + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getLdaVxcGradient<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mol_grad,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_dengrad_x,
                           d_dengrad_y,
                           d_dengrad_z,
                           static_cast<uint32_t>(natoms),
                           d_vrho);
    }

    // copy final gradient back
    gpuSafe(gpuMemcpyAsync(molgrad_omp.row(gpu_id), d_mol_grad, natoms * 3 * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));
    gpuSafe(gpuFreeAsync(d_ao_to_atom_ids, stream));

    gpuSafe(gpuFreeAsync(d_mol_grad, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_rw_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_gs_den_mat_full, stream));

    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_x, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_y, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_z, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_exc, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());
    }

    CDenseMatrix molgrad(natoms, 3);

    molgrad.zero();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        for (int64_t a = 0; a < natoms; a++)
        {
            for (int64_t d = 0; d < 3; d++)
            {
                molgrad.row(a)[d] += molgrad_omp.row(gpu_id)[a * 3 + d];
            }
        }
    }

    return molgrad;
}

static auto
integrateVxcGradientForLdaOpenShell(const CMolecule&        molecule,
                                    const CMolecularBasis&  basis,
                                    const CAODensityMatrix& rwDensityMatrix,
                                    const CAODensityMatrix& gsDensityMatrix,
                                    const CMolecularGrid&   molecularGrid,
                                    const CXCFunctional&    xcFunctional,
                                    const int64_t           numGpusPerNode,
                                    const int64_t           rank,
                                    const int64_t           nnodes) -> CDenseMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    const auto errnaos = std::string(__func__) + std::string(": Inconsistent number of AOs");

    errors::assertMsgCritical((naos == rwDensityMatrix.getNumberOfRows(0)) &&
                              (naos == rwDensityMatrix.getNumberOfColumns(0)) &&
                              (naos == gsDensityMatrix.getNumberOfRows(0)) &&
                              (naos == gsDensityMatrix.getNumberOfColumns(0)),
                              errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    // AO-to-atom mapping

    std::vector<int64_t> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix molgrad_omp(num_gpus_per_node, natoms * 3);

    molgrad_omp.zero();

    // grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    uint32_t* d_ao_to_atom_ids;

    gpuSafe(gpuMallocAsync(&d_ao_to_atom_ids, naos * sizeof(uint32_t), stream));

    // GTOs on grid points

    double *d_mol_grad, *d_rw_den_mat_full, *d_gs_den_mat_full, *d_den_mat, *d_mat_F;
    double *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z;
    double *d_dengrad_x, *d_dengrad_y, *d_dengrad_z;

    gpuSafe(gpuMallocAsync(&d_mol_grad, natoms * 3 * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_rw_den_mat_full, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gs_den_mat_full, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * 2 * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * 2 * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_x, natoms * max_npoints_per_box * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_y, natoms * max_npoints_per_box * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_z, natoms * max_npoints_per_box * 2 * sizeof(double), stream));

    double* d_rw_den_mat_full_a = d_rw_den_mat_full;
    double* d_rw_den_mat_full_b = d_rw_den_mat_full_a + naos * naos;

    double* d_gs_den_mat_full_a = d_gs_den_mat_full;
    double* d_gs_den_mat_full_b = d_gs_den_mat_full_a + naos * naos;

    double* d_den_mat_a = d_den_mat;
    double* d_den_mat_b = d_den_mat_a + naos * naos;

    double* d_mat_F_a = d_mat_F;
    double* d_mat_F_b = d_mat_F_a + naos * max_npoints_per_box;

    double* d_dengrad_x_a = d_dengrad_x;
    double* d_dengrad_x_b = d_dengrad_x_a + natoms * max_npoints_per_box;
    double* d_dengrad_y_a = d_dengrad_y;
    double* d_dengrad_y_b = d_dengrad_y_a + natoms * max_npoints_per_box;
    double* d_dengrad_z_a = d_dengrad_z;
    double* d_dengrad_z_b = d_dengrad_z_a + natoms * max_npoints_per_box;

    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full_a, rwDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full_b, rwDensityMatrix.betaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full_a, gsDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full_b, gsDensityMatrix.betaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ldafunc = xcfun_copy.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    double *d_rho, *d_exc, *d_vrho;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_exc, dim->zk * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(1, TILE_DIM);

    dim3 num_blocks(1, (natoms * 3 + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mol_grad, static_cast<uint32_t>(natoms * 3), 1);

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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

        // AO indices and AO-to-Atom indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        std::vector<uint32_t> ao_to_atom_ids_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);

            ao_to_atom_ids_int32[ind] = static_cast<uint32_t>(ao_to_atom_ids[aoinds[ind]]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        gpuSafe(gpuMemcpyAsync(d_ao_to_atom_ids, ao_to_atom_ids_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_a,
                           d_gs_den_mat_full_a,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_b,
                           d_gs_den_mat_full_b,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // LDA density for ground-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_a,
                           d_den_mat_a,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_b,
                           d_den_mat_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensityOnGridsOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_mat_F_a,
                           d_mat_F_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // density gradient on grid points

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_a,
                           d_rw_den_mat_full_a,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_b,
                           d_rw_den_mat_full_b,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_a,
                           d_den_mat_a,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_b,
                           d_den_mat_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (natoms + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getLdaDensityGradientOnGridsOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_dengrad_x_a,
                           d_dengrad_x_b,
                           d_dengrad_y_a,
                           d_dengrad_y_b,
                           d_dengrad_z_a,
                           d_dengrad_z_b,
                           d_mat_F_a,
                           d_mat_F_b,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           d_ao_to_atom_ids,
                           static_cast<uint32_t>(natoms),
                           static_cast<uint32_t>(npoints));

        // functional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        gpuSafe(gpuMemcpyAsync(d_exc, exc, dim->zk * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // accumulate partial contribution to Vxc gradient

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((natoms + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getLdaVxcGradientOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mol_grad,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_dengrad_x_a,
                           d_dengrad_x_b,
                           d_dengrad_y_a,
                           d_dengrad_y_b,
                           d_dengrad_z_a,
                           d_dengrad_z_b,
                           static_cast<uint32_t>(natoms),
                           d_vrho);
    }

    // copy final gradient back
    gpuSafe(gpuMemcpyAsync(molgrad_omp.row(gpu_id), d_mol_grad, natoms * 3 * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));
    gpuSafe(gpuFreeAsync(d_ao_to_atom_ids, stream));

    gpuSafe(gpuFreeAsync(d_mol_grad, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_rw_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_gs_den_mat_full, stream));

    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_x, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_y, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_z, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_exc, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());
    }

    CDenseMatrix molgrad(natoms, 3);

    molgrad.zero();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        for (int64_t a = 0; a < natoms; a++)
        {
            for (int64_t d = 0; d < 3; d++)
            {
                molgrad.row(a)[d] += molgrad_omp.row(gpu_id)[a * 3 + d];
            }
        }
    }

    return molgrad;
}

static auto
integrateVxcGradientForGgaClosedShell(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& rwDensityMatrix,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecularGrid&   molecularGrid,
                                      const CXCFunctional&    xcFunctional,
                                      const int64_t           numGpusPerNode,
                                      const int64_t           rank,
                                      const int64_t           nnodes) -> CDenseMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    const auto errnaos = std::string(__func__) + std::string(": Inconsistent number of AOs");

    errors::assertMsgCritical((naos == rwDensityMatrix.getNumberOfRows(0)) &&
                              (naos == rwDensityMatrix.getNumberOfColumns(0)) &&
                              (naos == gsDensityMatrix.getNumberOfRows(0)) &&
                              (naos == gsDensityMatrix.getNumberOfColumns(0)),
                              errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    // AO-to-atom mapping

    std::vector<int64_t> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix molgrad_omp(num_gpus_per_node, natoms * 3);

    molgrad_omp.zero();

    // grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    uint32_t* d_ao_to_atom_ids;

    gpuSafe(gpuMallocAsync(&d_ao_to_atom_ids, naos * sizeof(uint32_t), stream));

    // GTOs on grid points

    double *d_mol_grad, *d_rw_den_mat_full, *d_gs_den_mat_full, *d_den_mat;
    double *d_mat_F, *d_mat_F_x, *d_mat_F_y, *d_mat_F_z;
    double *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z;
    double *d_gto_values_xx, *d_gto_values_xy, *d_gto_values_xz, *d_gto_values_yy, *d_gto_values_yz, *d_gto_values_zz;

    double *d_dengrad_x, *d_dengrad_y, *d_dengrad_z;
    double *d_dengrad_xx, *d_dengrad_xy, *d_dengrad_xz;
    double *d_dengrad_yx, *d_dengrad_yy, *d_dengrad_yz;
    double *d_dengrad_zx, *d_dengrad_zy, *d_dengrad_zz;

    gpuSafe(gpuMallocAsync(&d_mol_grad, natoms * 3 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rw_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gs_den_mat_full, naos * naos * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F_z, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_gto_values_xx, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_xy, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_xz, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_yy, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_yz, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_zz, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_x, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_y, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_z, natoms * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_xx, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_xy, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_xz, natoms * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_yx, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_yy, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_yz, natoms * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_zx, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_zy, natoms * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_zz, natoms * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full, rwDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full, gsDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ggafunc = xcfun_copy.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    double *d_rho, *d_rhograd, *d_sigma, *d_vrho, *d_vsigma;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhograd, dim->rho * 3 * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_sigma, dim->sigma * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vsigma, dim->vsigma * max_npoints_per_box * sizeof(double), stream));

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(1, TILE_DIM);

    dim3 num_blocks(1, (natoms * 3 + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mol_grad, static_cast<uint32_t>(natoms * 3), 1);

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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
            // 2nd order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 2, 1.0e-12, boxdim);

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

            gpu::getGtoValuesForMgga(d_gto_values,
                                     d_gto_values_x,
                                     d_gto_values_y,
                                     d_gto_values_z,
                                     d_gto_values_xx,
                                     d_gto_values_xy,
                                     d_gto_values_xz,
                                     d_gto_values_yy,
                                     d_gto_values_yz,
                                     d_gto_values_zz,
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        std::vector<uint32_t> ao_to_atom_ids_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);

            ao_to_atom_ids_int32[ind] = static_cast<uint32_t>(ao_to_atom_ids[aoinds[ind]]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        gpuSafe(gpuMemcpyAsync(d_ao_to_atom_ids, ao_to_atom_ids_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_gs_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // GGA density for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensitySigmaOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_rhograd,
                           d_sigma,
                           d_mat_F,
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // density gradient on grid points

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat,
                           d_rw_den_mat_full,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F,
                           d_den_mat,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_x,
                           d_den_mat,
                           d_gto_values_x,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_y,
                           d_den_mat,
                           d_gto_values_y,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_z,
                           d_den_mat,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (natoms + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getGgaDensityGradientOnGrids<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_dengrad_x,
                           d_dengrad_y,
                           d_dengrad_z,
                           d_dengrad_xx,
                           d_dengrad_xy,
                           d_dengrad_xz,
                           d_dengrad_yx,
                           d_dengrad_yy,
                           d_dengrad_yz,
                           d_dengrad_zx,
                           d_dengrad_zy,
                           d_dengrad_zz,
                           d_mat_F,
                           d_mat_F_x,
                           d_mat_F_y,
                           d_mat_F_z,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           d_gto_values_xx,
                           d_gto_values_xy,
                           d_gto_values_xz,
                           d_gto_values_yy,
                           d_gto_values_yz,
                           d_gto_values_zz,
                           static_cast<uint32_t>(aocount),
                           d_ao_to_atom_ids,
                           static_cast<uint32_t>(natoms),
                           static_cast<uint32_t>(npoints));

        // funtional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(rhograd, d_rhograd, dim->rho * 3 * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(sigma, d_sigma, dim->sigma * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vsigma, vsigma, dim->vsigma * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // accumulate partial contribution to Vxc gradient

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((natoms + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getGgaVxcGradient<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mol_grad,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_dengrad_x,
                           d_dengrad_y,
                           d_dengrad_z,
                           d_dengrad_xx,
                           d_dengrad_xy,
                           d_dengrad_xz,
                           d_dengrad_yx,
                           d_dengrad_yy,
                           d_dengrad_yz,
                           d_dengrad_zx,
                           d_dengrad_zy,
                           d_dengrad_zz,
                           static_cast<uint32_t>(natoms),
                           d_vrho,
                           d_vsigma,
                           d_rhograd);
    }

    // copy final gradient back
    gpuSafe(gpuMemcpyAsync(molgrad_omp.row(gpu_id), d_mol_grad, natoms * 3 * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));
    gpuSafe(gpuFreeAsync(d_ao_to_atom_ids, stream));

    gpuSafe(gpuFreeAsync(d_mol_grad, stream));

    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z, stream));

    gpuSafe(gpuFreeAsync(d_gto_values_xx, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_xy, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_xz, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_yy, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_yz, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_zz, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_x, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_y, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_z, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_xx, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_xy, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_xz, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_yx, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_yy, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_yz, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_zx, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_zy, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_zz, stream));

    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_mat_F_x, stream));
    gpuSafe(gpuFreeAsync(d_mat_F_y, stream));
    gpuSafe(gpuFreeAsync(d_mat_F_z, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_rw_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_gs_den_mat_full, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_rhograd, stream));
    gpuSafe(gpuFreeAsync(d_sigma, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));
    gpuSafe(gpuFreeAsync(d_vsigma, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());
    }

    CDenseMatrix molgrad(natoms, 3);

    molgrad.zero();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        for (int64_t a = 0; a < natoms; a++)
        {
            for (int64_t d = 0; d < 3; d++)
            {
                molgrad.row(a)[d] += molgrad_omp.row(gpu_id)[a * 3 + d];
            }
        }
    }

    return molgrad;
}

static auto
integrateVxcGradientForGgaOpenShell(const CMolecule&        molecule,
                                    const CMolecularBasis&  basis,
                                    const CAODensityMatrix& rwDensityMatrix,
                                    const CAODensityMatrix& gsDensityMatrix,
                                    const CMolecularGrid&   molecularGrid,
                                    const CXCFunctional&    xcFunctional,
                                    const int64_t           numGpusPerNode,
                                    const int64_t           rank,
                                    const int64_t           nnodes) -> CDenseMatrix
{
    const auto num_gpus_per_node = numGpusPerNode;

    CGpuDevices gpu_devices;
    const auto ndevices = gpu_devices.getNumberOfDevices();

    const auto nthreads = omp_get_max_threads();

    errors::assertMsgCritical(
        static_cast<int64_t>(ndevices) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of devices does not match numGpusPerNode."));

    errors::assertMsgCritical(
        static_cast<int64_t>(nthreads) == num_gpus_per_node,
        std::string(__func__) + std::string(": Number of OMP threads does not match numGpusPerNode."));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    const auto errnaos = std::string(__func__) + std::string(": Inconsistent number of AOs");

    errors::assertMsgCritical((naos == rwDensityMatrix.getNumberOfRows(0)) &&
                              (naos == rwDensityMatrix.getNumberOfColumns(0)) &&
                              (naos == gsDensityMatrix.getNumberOfRows(0)) &&
                              (naos == gsDensityMatrix.getNumberOfColumns(0)),
                              errnaos);

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    // AO-to-atom mapping

    std::vector<int64_t> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix molgrad_omp(num_gpus_per_node, natoms * 3);

    molgrad_omp.zero();

    // grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

#pragma omp parallel
    {
    const auto thread_id = omp_get_thread_num();

    const auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    double* d_gto_info;

    gpuSafe(gpuMallocAsync(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double), stream));

    uint32_t* d_ao_inds;

    gpuSafe(gpuMallocAsync(&d_ao_inds, naos * sizeof(uint32_t), stream));

    uint32_t* d_ao_to_atom_ids;

    gpuSafe(gpuMallocAsync(&d_ao_to_atom_ids, naos * sizeof(uint32_t), stream));

    // GTOs on grid points

    double *d_mol_grad, *d_rw_den_mat_full, *d_gs_den_mat_full, *d_den_mat;
    double *d_mat_F, *d_mat_F_x, *d_mat_F_y, *d_mat_F_z;
    double *d_gto_values, *d_gto_values_x, *d_gto_values_y, *d_gto_values_z;
    double *d_gto_values_xx, *d_gto_values_xy, *d_gto_values_xz, *d_gto_values_yy, *d_gto_values_yz, *d_gto_values_zz;

    double *d_dengrad_x, *d_dengrad_y, *d_dengrad_z;
    double *d_dengrad_xx, *d_dengrad_xy, *d_dengrad_xz;
    double *d_dengrad_yx, *d_dengrad_yy, *d_dengrad_yz;
    double *d_dengrad_zx, *d_dengrad_zy, *d_dengrad_zz;

    gpuSafe(gpuMallocAsync(&d_mol_grad, natoms * 3 * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_rw_den_mat_full, naos * naos * 2 * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gs_den_mat_full, naos * naos * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_den_mat, naos * naos * 2  * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_mat_F, naos * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F_x, naos * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F_y, naos * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_mat_F_z, naos * max_npoints_per_box * 2  * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_gto_values, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_x, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_y, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_z, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_gto_values_xx, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_xy, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_xz, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_yy, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_yz, naos * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_gto_values_zz, naos * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_x, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_y, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_z, natoms * max_npoints_per_box * 2  * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_xx, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_xy, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_xz, natoms * max_npoints_per_box * 2  * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_yx, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_yy, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_yz, natoms * max_npoints_per_box * 2  * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_dengrad_zx, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_zy, natoms * max_npoints_per_box * 2  * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_dengrad_zz, natoms * max_npoints_per_box * 2  * sizeof(double), stream));

    double* d_rw_den_mat_full_a = d_rw_den_mat_full;
    double* d_rw_den_mat_full_b = d_rw_den_mat_full_a + naos * naos;

    double* d_gs_den_mat_full_a = d_gs_den_mat_full;
    double* d_gs_den_mat_full_b = d_gs_den_mat_full_a + naos * naos;

    double* d_den_mat_a = d_den_mat;
    double* d_den_mat_b = d_den_mat_a + naos * naos;

    double* d_mat_F_a = d_mat_F;
    double* d_mat_F_b = d_mat_F_a + naos * max_npoints_per_box;

    double* d_mat_F_x_a = d_mat_F_x;
    double* d_mat_F_x_b = d_mat_F_x_a + naos * max_npoints_per_box;
    double* d_mat_F_y_a = d_mat_F_y;
    double* d_mat_F_y_b = d_mat_F_y_a + naos * max_npoints_per_box;
    double* d_mat_F_z_a = d_mat_F_z;
    double* d_mat_F_z_b = d_mat_F_z_a + naos * max_npoints_per_box;

    double* d_dengrad_x_a = d_dengrad_x;
    double* d_dengrad_x_b = d_dengrad_x_a + natoms * max_npoints_per_box;
    double* d_dengrad_y_a = d_dengrad_y;
    double* d_dengrad_y_b = d_dengrad_y_a + natoms * max_npoints_per_box;
    double* d_dengrad_z_a = d_dengrad_z;
    double* d_dengrad_z_b = d_dengrad_z_a + natoms * max_npoints_per_box;

    double* d_dengrad_xx_a = d_dengrad_xx;
    double* d_dengrad_xx_b = d_dengrad_xx_a + natoms * max_npoints_per_box;
    double* d_dengrad_xy_a = d_dengrad_xy;
    double* d_dengrad_xy_b = d_dengrad_xy_a + natoms * max_npoints_per_box;
    double* d_dengrad_xz_a = d_dengrad_xz;
    double* d_dengrad_xz_b = d_dengrad_xz_a + natoms * max_npoints_per_box;

    double* d_dengrad_yx_a = d_dengrad_yx;
    double* d_dengrad_yx_b = d_dengrad_yx_a + natoms * max_npoints_per_box;
    double* d_dengrad_yy_a = d_dengrad_yy;
    double* d_dengrad_yy_b = d_dengrad_yy_a + natoms * max_npoints_per_box;
    double* d_dengrad_yz_a = d_dengrad_yz;
    double* d_dengrad_yz_b = d_dengrad_yz_a + natoms * max_npoints_per_box;

    double* d_dengrad_zx_a = d_dengrad_zx;
    double* d_dengrad_zx_b = d_dengrad_zx_a + natoms * max_npoints_per_box;
    double* d_dengrad_zy_a = d_dengrad_zy;
    double* d_dengrad_zy_b = d_dengrad_zy_a + natoms * max_npoints_per_box;
    double* d_dengrad_zz_a = d_dengrad_zz;
    double* d_dengrad_zz_b = d_dengrad_zz_a + natoms * max_npoints_per_box;

    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full_a, rwDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_rw_den_mat_full_b, rwDensityMatrix.betaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full_a, gsDensityMatrix.alphaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_gs_den_mat_full_b, gsDensityMatrix.betaDensity(0), naos * naos * sizeof(double), gpuMemcpyHostToDevice, stream));

    // density and functional derivatives

    auto xcfun_copy = vxcfuncs::getExchangeCorrelationFunctional(xcFunctional.getFunctionalLabel());

    auto       ggafunc = xcfun_copy.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    double *d_rho, *d_rhograd, *d_sigma, *d_vrho, *d_vsigma;

    gpuSafe(gpuMallocAsync(&d_rho, dim->rho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_rhograd, dim->rho * 3 * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_sigma, dim->sigma * max_npoints_per_box * sizeof(double), stream));

    gpuSafe(gpuMallocAsync(&d_vrho, dim->vrho * max_npoints_per_box * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_vsigma, dim->vsigma * max_npoints_per_box * sizeof(double), stream));

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    auto n_total_grid_points = molecularGrid.getNumberOfGridPoints();

    double *d_grid_x, *d_grid_y, *d_grid_z, *d_grid_w;

    gpuSafe(gpuMallocAsync(&d_grid_x, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_y, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_z, n_total_grid_points * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_grid_w, n_total_grid_points * sizeof(double), stream));

    gpuSafe(gpuMemcpyAsync(d_grid_x, xcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_y, ycoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_z, zcoords, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyAsync(d_grid_w, weights, n_total_grid_points * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    dim3 threads_per_block(1, TILE_DIM);

    dim3 num_blocks(1, (natoms * 3 + threads_per_block.y - 1) / threads_per_block.y);

    gpu::zeroMatrix<<<num_blocks, threads_per_block, 0, stream>>>(d_mol_grad, static_cast<uint32_t>(natoms * 3), 1);

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (size_t box_id = gpu_id; box_id < counts.size(); box_id += num_gpus_per_node)
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
            // 2nd order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 2, 1.0e-12, boxdim);

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

            gpu::getGtoValuesForMgga(d_gto_values,
                                     d_gto_values_x,
                                     d_gto_values_y,
                                     d_gto_values_z,
                                     d_gto_values_xx,
                                     d_gto_values_xy,
                                     d_gto_values_xz,
                                     d_gto_values_yy,
                                     d_gto_values_yz,
                                     d_gto_values_zz,
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

        // AO indices

        std::vector<uint32_t> ao_inds_int32(aocount);

        std::vector<uint32_t> ao_to_atom_ids_int32(aocount);

        for (int64_t ind = 0; ind < aocount; ind++)
        {
            ao_inds_int32[ind] = static_cast<uint32_t>(aoinds[ind]);

            ao_to_atom_ids_int32[ind] = static_cast<uint32_t>(ao_to_atom_ids[aoinds[ind]]);
        }

        gpuSafe(gpuMemcpyAsync(d_ao_inds, ao_inds_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        gpuSafe(gpuMemcpyAsync(d_ao_to_atom_ids, ao_to_atom_ids_int32.data(), aocount * sizeof(uint32_t), gpuMemcpyHostToDevice, stream));

        // sub density matrix for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_a,
                           d_gs_den_mat_full_a,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_b,
                           d_gs_den_mat_full_b,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        // GGA density for groud-state rho

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_a,
                           d_den_mat_a,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_b,
                           d_den_mat_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        // Note: one block in y dimension
        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getDensitySigmaOnGridsOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_rho,
                           d_rhograd,
                           d_sigma,
                           d_mat_F_a,
                           d_mat_F_b,
                           d_gto_values,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        // density gradient on grid points

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((aocount + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_a,
                           d_rw_den_mat_full_a,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        gpu::getSubDensityMatrix<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_den_mat_b,
                           d_rw_den_mat_full_b,
                           static_cast<uint32_t>(naos),
                           d_ao_inds,
                           static_cast<uint32_t>(aocount));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (aocount + threads_per_block.y - 1) / threads_per_block.y);

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_a,
                           d_den_mat_a,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_b,
                           d_den_mat_b,
                           d_gto_values,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_x_a,
                           d_den_mat_a,
                           d_gto_values_x,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_x_b,
                           d_den_mat_b,
                           d_gto_values_x,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_y_a,
                           d_den_mat_a,
                           d_gto_values_y,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_y_b,
                           d_den_mat_b,
                           d_gto_values_y,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_z_a,
                           d_den_mat_a,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        gpu::matmulAB<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mat_F_z_b,
                           d_den_mat_b,
                           d_gto_values_z,
                           static_cast<uint32_t>(aocount),
                           static_cast<uint32_t>(npoints));

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((npoints + threads_per_block.x - 1) / threads_per_block.x, (natoms + threads_per_block.y - 1) / threads_per_block.y);

        gpu::getGgaDensityGradientOnGridsOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_dengrad_x_a,
                           d_dengrad_x_b,
                           d_dengrad_y_a,
                           d_dengrad_y_b,
                           d_dengrad_z_a,
                           d_dengrad_z_b,
                           d_dengrad_xx_a,
                           d_dengrad_xx_b,
                           d_dengrad_xy_a,
                           d_dengrad_xy_b,
                           d_dengrad_xz_a,
                           d_dengrad_xz_b,
                           d_dengrad_yx_a,
                           d_dengrad_yx_b,
                           d_dengrad_yy_a,
                           d_dengrad_yy_b,
                           d_dengrad_yz_a,
                           d_dengrad_yz_b,
                           d_dengrad_zx_a,
                           d_dengrad_zx_b,
                           d_dengrad_zy_a,
                           d_dengrad_zy_b,
                           d_dengrad_zz_a,
                           d_dengrad_zz_b,
                           d_mat_F_a,
                           d_mat_F_b,
                           d_mat_F_x_a,
                           d_mat_F_x_b,
                           d_mat_F_y_a,
                           d_mat_F_y_b,
                           d_mat_F_z_a,
                           d_mat_F_z_b,
                           d_gto_values_x,
                           d_gto_values_y,
                           d_gto_values_z,
                           d_gto_values_xx,
                           d_gto_values_xy,
                           d_gto_values_xz,
                           d_gto_values_yy,
                           d_gto_values_yz,
                           d_gto_values_zz,
                           static_cast<uint32_t>(aocount),
                           d_ao_to_atom_ids,
                           static_cast<uint32_t>(natoms),
                           static_cast<uint32_t>(npoints));

        // funtional evaluation

        gpuSafe(gpuMemcpyAsync(rho, d_rho, dim->rho * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(rhograd, d_rhograd, dim->rho * 3 * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));
        gpuSafe(gpuMemcpyAsync(sigma, d_sigma, dim->sigma * npoints * sizeof(double), gpuMemcpyDeviceToHost, stream));

        gpuSafe(gpuStreamSynchronize(stream));

        xcfun_copy.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        gpuSafe(gpuMemcpyAsync(d_vrho, vrho, dim->vrho * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));
        gpuSafe(gpuMemcpyAsync(d_vsigma, vsigma, dim->vsigma * npoints * sizeof(double), gpuMemcpyHostToDevice, stream));

        // accumulate partial contribution to Vxc gradient

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((natoms + threads_per_block.x - 1) / threads_per_block.x, 1);

        gpu::getGgaVxcGradientOpenShell<<<num_blocks, threads_per_block, 0, stream>>>(
                           d_mol_grad,
                           d_grid_w,
                           static_cast<uint32_t>(gridblockpos),
                           static_cast<uint32_t>(npoints),
                           d_dengrad_x_a,
                           d_dengrad_x_b,
                           d_dengrad_y_a,
                           d_dengrad_y_b,
                           d_dengrad_z_a,
                           d_dengrad_z_b,
                           d_dengrad_xx_a,
                           d_dengrad_xx_b,
                           d_dengrad_xy_a,
                           d_dengrad_xy_b,
                           d_dengrad_xz_a,
                           d_dengrad_xz_b,
                           d_dengrad_yx_a,
                           d_dengrad_yx_b,
                           d_dengrad_yy_a,
                           d_dengrad_yy_b,
                           d_dengrad_yz_a,
                           d_dengrad_yz_b,
                           d_dengrad_zx_a,
                           d_dengrad_zx_b,
                           d_dengrad_zy_a,
                           d_dengrad_zy_b,
                           d_dengrad_zz_a,
                           d_dengrad_zz_b,
                           static_cast<uint32_t>(natoms),
                           d_vrho,
                           d_vsigma,
                           d_rhograd);
    }

    // copy final gradient back
    gpuSafe(gpuMemcpyAsync(molgrad_omp.row(gpu_id), d_mol_grad, natoms * 3 * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

#pragma omp barrier

    gpuSafe(gpuFreeAsync(d_gto_info, stream));
    gpuSafe(gpuFreeAsync(d_ao_inds, stream));
    gpuSafe(gpuFreeAsync(d_ao_to_atom_ids, stream));

    gpuSafe(gpuFreeAsync(d_mol_grad, stream));

    gpuSafe(gpuFreeAsync(d_gto_values, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_x, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_y, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_z, stream));

    gpuSafe(gpuFreeAsync(d_gto_values_xx, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_xy, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_xz, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_yy, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_yz, stream));
    gpuSafe(gpuFreeAsync(d_gto_values_zz, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_x, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_y, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_z, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_xx, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_xy, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_xz, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_yx, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_yy, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_yz, stream));

    gpuSafe(gpuFreeAsync(d_dengrad_zx, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_zy, stream));
    gpuSafe(gpuFreeAsync(d_dengrad_zz, stream));

    gpuSafe(gpuFreeAsync(d_mat_F, stream));
    gpuSafe(gpuFreeAsync(d_mat_F_x, stream));
    gpuSafe(gpuFreeAsync(d_mat_F_y, stream));
    gpuSafe(gpuFreeAsync(d_mat_F_z, stream));

    gpuSafe(gpuFreeAsync(d_den_mat, stream));
    gpuSafe(gpuFreeAsync(d_rw_den_mat_full, stream));
    gpuSafe(gpuFreeAsync(d_gs_den_mat_full, stream));

    gpuSafe(gpuFreeAsync(d_rho, stream));
    gpuSafe(gpuFreeAsync(d_rhograd, stream));
    gpuSafe(gpuFreeAsync(d_sigma, stream));
    gpuSafe(gpuFreeAsync(d_vrho, stream));
    gpuSafe(gpuFreeAsync(d_vsigma, stream));

    gpuSafe(gpuFreeAsync(d_grid_x, stream));
    gpuSafe(gpuFreeAsync(d_grid_y, stream));
    gpuSafe(gpuFreeAsync(d_grid_z, stream));
    gpuSafe(gpuFreeAsync(d_grid_w, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(gpuDeviceSynchronize());
    }

    CDenseMatrix molgrad(natoms, 3);

    molgrad.zero();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        for (int64_t a = 0; a < natoms; a++)
        {
            for (int64_t d = 0; d < 3; d++)
            {
                molgrad.row(a)[d] += molgrad_omp.row(gpu_id)[a * 3 + d];
            }
        }
    }

    return molgrad;
}

auto
integrateVxcGradient(const CMolecule&        molecule,
                     const CMolecularBasis&  basis,
                     const CAODensityMatrix& rwDensityMatrix,
                     const CAODensityMatrix& gsDensityMatrix,
                     const CMolecularGrid&   molecularGrid,
                     const std::string&      xcFuncLabel,
                     const int64_t           numGpusPerNode,
                     const int64_t           rank,
                     const int64_t           nnodes) -> CDenseMatrix
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    auto flag = (rwDensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell()) ? std::string("CLOSEDSHELL") : std::string("OPENSHELL");

    if (xcfuntype == xcfun::lda)
    {
        if (flag == std::string("CLOSEDSHELL"))
        {
            return gpu::integrateVxcGradientForLdaClosedShell(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
        else
        {
            return gpu::integrateVxcGradientForLdaOpenShell(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
    }
    else if (xcfuntype == xcfun::gga)
    {
        if (flag == std::string("CLOSEDSHELL"))
        {
            return gpu::integrateVxcGradientForGgaClosedShell(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
        else
        {
            return gpu::integrateVxcGradientForGgaOpenShell(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc, numGpusPerNode, rank, nnodes);
        }
    }
    else
    {
        std::string errxcfuntype("gpu::integrateVxcGradient: Only implemented for LDA/GGA");

        errors::assertMsgCritical(false, errxcfuntype);

        return CDenseMatrix();
    }
}

}  // namespace gpu
