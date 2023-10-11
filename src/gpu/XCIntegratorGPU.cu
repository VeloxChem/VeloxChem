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

#include <omp.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "DftFunc.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "GtoValuesRecF.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "Prescreener.hpp"
#include "XCIntegratorGPU.hpp"

#define cudaSafe(err)                                                                                                     \
    {                                                                                                                     \
        if (err != cudaSuccess)                                                                                           \
        {                                                                                                                 \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                      \
        }                                                                                                                 \
    }

namespace gpu {  // gpu namespace

__global__ void
cudaLdaValuesRecS(double*        gto_values,
                  const double*  gto_info,
                  const double*  grid_xyz,
                  const uint32_t nrows,
                  const uint32_t npgtos,
                  const uint32_t ncols)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t g = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        const auto g_x = grid_xyz[g + ncols * 0];
        const auto g_y = grid_xyz[g + ncols * 1];
        const auto g_z = grid_xyz[g + ncols * 2];

        gto_values[g + i * ncols] = 0.0;

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

            gto_values[g + i * ncols] += fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
        }
    }
}

__global__ void
cudaLdaValuesRecP(double*        gto_values_p3,
                  const double*  gto_info,
                  const double*  grid_xyz,
                  const uint32_t nrows,
                  const uint32_t npgtos,
                  const uint32_t ncols)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t g = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        const auto g_x = grid_xyz[g + ncols * 0];
        const auto g_y = grid_xyz[g + ncols * 1];
        const auto g_z = grid_xyz[g + ncols * 2];

        gto_values_p3[g + i * ncols + nrows * ncols * 0] = 0.0;
        gto_values_p3[g + i * ncols + nrows * ncols * 1] = 0.0;
        gto_values_p3[g + i * ncols + nrows * ncols * 2] = 0.0;

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

            gto_values_p3[g + i * ncols + nrows * ncols * 0] += gr_x * fss;
            gto_values_p3[g + i * ncols + nrows * ncols * 1] += gr_y * fss;
            gto_values_p3[g + i * ncols + nrows * ncols * 2] += gr_z * fss;
        }
    }
}

__global__ void
cudaLdaValuesRecD(double*        gto_values_d6,
                  const double*  gto_info,
                  const double*  grid_xyz,
                  const uint32_t nrows,
                  const uint32_t npgtos,
                  const uint32_t ncols)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t g = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        const auto g_x = grid_xyz[g];
        const auto g_y = grid_xyz[g + ncols];
        const auto g_z = grid_xyz[g + ncols * 2];

        gto_values_d6[g + i * ncols + nrows * ncols * 0] = 0.0;
        gto_values_d6[g + i * ncols + nrows * ncols * 1] = 0.0;
        gto_values_d6[g + i * ncols + nrows * ncols * 2] = 0.0;
        gto_values_d6[g + i * ncols + nrows * ncols * 3] = 0.0;
        gto_values_d6[g + i * ncols + nrows * ncols * 4] = 0.0;
        gto_values_d6[g + i * ncols + nrows * ncols * 5] = 0.0;

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

            gto_values_d6[g + i * ncols + nrows * ncols * 0] += gr_x * gr_x * fss;
            gto_values_d6[g + i * ncols + nrows * ncols * 1] += gr_x * gr_y * fss;
            gto_values_d6[g + i * ncols + nrows * ncols * 2] += gr_x * gr_z * fss;
            gto_values_d6[g + i * ncols + nrows * ncols * 3] += gr_y * gr_y * fss;
            gto_values_d6[g + i * ncols + nrows * ncols * 4] += gr_y * gr_z * fss;
            gto_values_d6[g + i * ncols + nrows * ncols * 5] += gr_z * gr_z * fss;
        }
    }
}

auto
getLdaValuesRecS(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int64_t>(grid_coords_x.size());

    auto gto_values = matfunc::makeMatrix("LDA", nrows, ncols);

    auto submat = gto_values.getSubMatrix({0, 0});

    submat->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // set up data on host and device

    double *h_gto_values, *h_gto_info, *h_grid_xyz;

    cudaSafe(cudaMallocHost(&h_gto_values, nrows * ncols * sizeof(double)));
    cudaSafe(cudaMallocHost(&h_gto_info, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMallocHost(&h_grid_xyz, 3 * ncols * sizeof(double)));

    double *d_gto_values, *d_gto_info, *d_grid_xyz;

    cudaSafe(cudaMalloc(&d_gto_values, nrows * ncols * sizeof(double)));
    cudaSafe(cudaMalloc(&d_gto_info, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_xyz, 3 * ncols * sizeof(double)));

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

                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 0] = fexp;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 1] = fnorm;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 2] = r_x;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 3] = r_y;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 4] = r_z;
            }

            irow++;
        }
    }

    for (int64_t k = 0; k < ncols; k++)
    {
        h_grid_xyz[k + ncols * 0] = g_x[k];
        h_grid_xyz[k + ncols * 1] = g_y[k];
        h_grid_xyz[k + ncols * 2] = g_z[k];
    }

    cudaSafe(cudaMemcpy(d_gto_info, h_gto_info, 5 * nrows * npgtos * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_xyz, h_grid_xyz, 3 * ncols * sizeof(double), cudaMemcpyHostToDevice));

    dim3 threads_per_block(8, 32);
    dim3 num_blocks((nrows + threads_per_block.x - 1) / threads_per_block.x, (ncols + threads_per_block.y - 1) / threads_per_block.y);

    gpu::cudaLdaValuesRecS<<<num_blocks, threads_per_block>>>(
        d_gto_values, d_gto_info, d_grid_xyz, static_cast<uint32_t>(nrows), static_cast<uint32_t>(npgtos), static_cast<uint32_t>(ncols));

    cudaSafe(cudaMemcpy(h_gto_values, d_gto_values, nrows * ncols * sizeof(double), cudaMemcpyDeviceToHost));

    for (int64_t irow = 0; irow < nrows; irow++)
    {
        for (int64_t k = 0; k < ncols; k++)
        {
            submat->at(irow, k, false) += h_gto_values[k + irow * ncols];
        }
    }

    cudaSafe(cudaFreeHost(h_gto_values));
    cudaSafe(cudaFreeHost(h_gto_info));
    cudaSafe(cudaFreeHost(h_grid_xyz));

    cudaSafe(cudaFree(d_gto_values));
    cudaSafe(cudaFree(d_gto_info));
    cudaSafe(cudaFree(d_grid_xyz));

    return gto_values;
}

auto
getLdaValuesRecP(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int64_t>(grid_coords_x.size());

    auto gto_values = matfunc::makeMatrix("LDA", 3 * nrows, ncols);

    auto submat = gto_values.getSubMatrix({0, 0});

    submat->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // set up data on host and device

    double *h_gto_values, *h_gto_info, *h_grid_xyz;

    cudaSafe(cudaMallocHost(&h_gto_values, 3 * nrows * ncols * sizeof(double)));
    cudaSafe(cudaMallocHost(&h_gto_info, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMallocHost(&h_grid_xyz, 3 * ncols * sizeof(double)));

    double *d_gto_values, *d_gto_info, *d_grid_xyz;

    cudaSafe(cudaMalloc(&d_gto_values, 3 * nrows * ncols * sizeof(double)));
    cudaSafe(cudaMalloc(&d_gto_info, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_xyz, 3 * ncols * sizeof(double)));

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

                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 0] = fexp;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 1] = fnorm;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 2] = r_x;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 3] = r_y;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 4] = r_z;
            }

            irow++;
        }
    }

    for (int64_t k = 0; k < ncols; k++)
    {
        h_grid_xyz[k + ncols * 0] = g_x[k];
        h_grid_xyz[k + ncols * 1] = g_y[k];
        h_grid_xyz[k + ncols * 2] = g_z[k];
    }

    cudaSafe(cudaMemcpy(d_gto_info, h_gto_info, 5 * nrows * npgtos * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_xyz, h_grid_xyz, 3 * ncols * sizeof(double), cudaMemcpyHostToDevice));

    dim3 threads_per_block(8, 32);
    dim3 num_blocks((nrows + threads_per_block.x - 1) / threads_per_block.x, (ncols + threads_per_block.y - 1) / threads_per_block.y);

    gpu::cudaLdaValuesRecP<<<num_blocks, threads_per_block>>>(
        d_gto_values, d_gto_info, d_grid_xyz, static_cast<uint32_t>(nrows), static_cast<uint32_t>(npgtos), static_cast<uint32_t>(ncols));

    cudaSafe(cudaMemcpy(h_gto_values, d_gto_values, 3 * nrows * ncols * sizeof(double), cudaMemcpyDeviceToHost));

    for (int64_t irow = 0; irow < nrows; irow++)
    {
        for (int64_t k = 0; k < ncols; k++)
        {
            // buffer_x, 2 * nrows + irow
            // buffer_y,             irow
            // buffer_z,     nrows + irow

            submat->at(irow + nrows * 2, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 0];
            submat->at(irow + nrows * 0, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 1];
            submat->at(irow + nrows * 1, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 2];
        }
    }

    cudaSafe(cudaFreeHost(h_gto_values));
    cudaSafe(cudaFreeHost(h_gto_info));
    cudaSafe(cudaFreeHost(h_grid_xyz));

    cudaSafe(cudaFree(d_gto_values));
    cudaSafe(cudaFree(d_gto_info));
    cudaSafe(cudaFree(d_grid_xyz));

    return gto_values;
}

auto
getLdaValuesRecD(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int64_t>(grid_coords_x.size());

    auto gto_values = matfunc::makeMatrix("LDA", 5 * nrows, ncols);

    auto submat = gto_values.getSubMatrix({0, 0});

    submat->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // set up data on host and device

    double *h_gto_values, *h_gto_info, *h_grid_xyz;

    cudaSafe(cudaMallocHost(&h_gto_values, 6 * nrows * ncols * sizeof(double)));
    cudaSafe(cudaMallocHost(&h_gto_info, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMallocHost(&h_grid_xyz, 3 * ncols * sizeof(double)));

    double *d_gto_values, *d_gto_info, *d_grid_xyz;

    cudaSafe(cudaMalloc(&d_gto_values, 6 * nrows * ncols * sizeof(double)));
    cudaSafe(cudaMalloc(&d_gto_info, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMalloc(&d_grid_xyz, 3 * ncols * sizeof(double)));

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

                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 0] = fexp;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 1] = fnorm;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 2] = r_x;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 3] = r_y;
                h_gto_info[irow + j * ncgtos + npgtos * ncgtos * 4] = r_z;
            }

            irow++;
        }
    }

    for (int64_t k = 0; k < ncols; k++)
    {
        h_grid_xyz[k + ncols * 0] = g_x[k];
        h_grid_xyz[k + ncols * 1] = g_y[k];
        h_grid_xyz[k + ncols * 2] = g_z[k];
    }

    cudaSafe(cudaMemcpy(d_gto_info, h_gto_info, 5 * nrows * npgtos * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_grid_xyz, h_grid_xyz, 3 * ncols * sizeof(double), cudaMemcpyHostToDevice));

    dim3 threads_per_block(8, 32);
    dim3 num_blocks((nrows + threads_per_block.x - 1) / threads_per_block.x, (ncols + threads_per_block.y - 1) / threads_per_block.y);

    gpu::cudaLdaValuesRecD<<<num_blocks, threads_per_block>>>(
        d_gto_values, d_gto_info, d_grid_xyz, static_cast<uint32_t>(nrows), static_cast<uint32_t>(npgtos), static_cast<uint32_t>(ncols));

    cudaSafe(cudaMemcpy(h_gto_values, d_gto_values, 6 * nrows * ncols * sizeof(double), cudaMemcpyDeviceToHost));

    for (int64_t irow = 0; irow < nrows; irow++)
    {
        for (int64_t k = 0; k < ncols; k++)
        {
            // buffer_xx, -1.0,        2 * nrows + irow
            // buffer_xx, 0.5 * f2_3,  4 * nrows + irow
            // buffer_xy, f2_3,                    irow
            // buffer_xz, f2_3,        3 * nrows + irow
            // buffer_yy, -1.0,        2 * nrows + irow
            // buffer_yy, -0.5 * f2_3, 4 * nrows + irow
            // buffer_yz, f2_3,            nrows + irow
            // buffer_zz, 2.0,         2 * nrows + irow

            submat->at(irow + nrows * 2, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 0] * (-1.0);
            submat->at(irow + nrows * 4, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 0] * 0.5 * f2_3;
            submat->at(irow + nrows * 0, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 1] * f2_3;
            submat->at(irow + nrows * 3, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 2] * f2_3;
            submat->at(irow + nrows * 2, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 3] * (-1.0);
            submat->at(irow + nrows * 4, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 3] * (-0.5) * f2_3;
            submat->at(irow + nrows * 1, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 4] * f2_3;
            submat->at(irow + nrows * 2, k, false) += h_gto_values[k + irow * ncols + nrows * ncols * 5] * 2.0;
        }
    }

    cudaSafe(cudaFreeHost(h_gto_values));
    cudaSafe(cudaFreeHost(h_gto_info));
    cudaSafe(cudaFreeHost(h_grid_xyz));

    cudaSafe(cudaFree(d_gto_values));
    cudaSafe(cudaFree(d_gto_info));
    cudaSafe(cudaFree(d_grid_xyz));

    return gto_values;
}

auto
getGtoValuesForLda(const CGtoBlock&            gto_block,
                   const std::vector<double>&  grid_coords_x,
                   const std::vector<double>&  grid_coords_y,
                   const std::vector<double>&  grid_coords_z,
                   const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        return gpu::getLdaValuesRecS(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gpu::getLdaValuesRecP(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gpu::getLdaValuesRecD(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 3)
    {
        return gtoval::getLdaValuesRecF(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
}

auto
computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> CDenseMatrix
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

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

        // compute GTO values on grid points

        const auto grid_x_ptr = xcoords + gridblockpos;
        const auto grid_y_ptr = ycoords + gridblockpos;
        const auto grid_z_ptr = zcoords + gridblockpos;

        std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);
        std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);
        std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);

        // go through GTO blocks

        for (const auto& gto_block : gto_blocks)
        {
            // prescreen GTO block

            // 0th order GTO derivative
            auto pre_scr_info = prescr::preScreenGtoBlock(gto_block, 0, 1.0e-12, boxdim);

            auto cgto_mask   = std::get<0>(pre_scr_info);
            auto pre_ao_inds = std::get<1>(pre_scr_info);

            // GTO values on grid points

            auto cmat = gpu::getGtoValuesForLda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

            auto submat_ptr = cmat.getSubMatrix({0, 0});

            auto subgaos_ptr = submat_ptr->getData();

            for (int64_t nu = 0; nu < static_cast<int64_t>(pre_ao_inds.size()); nu++)
            {
                std::memcpy(allgtovalues.row(pre_ao_inds[nu]) + gridblockpos, subgaos_ptr + nu * npoints, npoints * sizeof(double));
            }
        }
    }

    return allgtovalues;
}

}  // namespace gpu
