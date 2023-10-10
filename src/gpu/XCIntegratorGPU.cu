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
#include "GtoValuesRecD.hpp"
#include "GtoValuesRecF.hpp"
#include "GtoValuesRecP.hpp"
#include "GtoValuesRecS.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "Prescreener.hpp"
#include "XCIntegratorGPU.hpp"

// clang-format off
#define cudaSafe(err) { gpu::cudaCheck(err); }
// clang-format on

namespace gpu {  // gpu namespace

inline void
cudaCheck(cudaError_t err)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA error " << __FILE__ << ":" << __LINE__ << " : " << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

__global__ void
comp_gto_on_grids(double*        gto_values,
                  const double*  gto_info,
                  const double*  grid_xyz,
                  const uint32_t nrows,
                  const uint32_t npgtos,
                  const uint32_t ncols)
{
    const uint32_t g = blockDim.x * blockIdx.x + threadIdx.x;

    if (g < ncols)
    {
        auto g_x = grid_xyz[g];
        auto g_y = grid_xyz[g + ncols];
        auto g_z = grid_xyz[g + ncols * 2];

        for (uint32_t i = 0; i < nrows * npgtos; i++)
        {
            auto fexp  = gto_info[i];
            auto fnorm = gto_info[i + nrows * npgtos];
            auto r_x   = gto_info[i + nrows * npgtos * 2];
            auto r_y   = gto_info[i + nrows * npgtos * 3];
            auto r_z   = gto_info[i + nrows * npgtos * 4];

            auto gr_x = g_x - r_x;
            auto gr_y = g_y - r_y;
            auto gr_z = g_z - r_z;

            // TODO use +=

            gto_values[i * ncols + g] = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
        }
    }
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

        {
            // auto thread_id = omp_get_thread_num();

            // auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            // auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            int64_t grid_batch_size   = npoints;
            int64_t grid_batch_offset = 0;

            const auto grid_x_ptr = xcoords + gridblockpos + grid_batch_offset;
            const auto grid_y_ptr = ycoords + gridblockpos + grid_batch_offset;
            const auto grid_z_ptr = zcoords + gridblockpos + grid_batch_offset;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + grid_batch_size);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + grid_batch_size);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + grid_batch_size);

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
                    std::memcpy(allgtovalues.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                }
            }
        }
    }

    return allgtovalues;
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
        return gtoval::getLdaValuesRecP(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::getLdaValuesRecD(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 3)
    {
        return gtoval::getLdaValuesRecF(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
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

    double *gto_values_h, *gto_info_h, *grid_xyz_h;

    cudaSafe(cudaMallocHost(&gto_values_h, nrows * npgtos * ncols * sizeof(double)));
    cudaSafe(cudaMallocHost(&gto_info_h, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMallocHost(&grid_xyz_h, 3 * ncols * sizeof(double)));

    double *gto_values_d, *gto_info_d, *grid_xyz_d;

    cudaSafe(cudaMalloc(&gto_values_d, nrows * npgtos * ncols * sizeof(double)));
    cudaSafe(cudaMalloc(&gto_info_d, 5 * nrows * npgtos * sizeof(double)));
    cudaSafe(cudaMalloc(&grid_xyz_d, 3 * ncols * sizeof(double)));

    for (int64_t i = 0, idx = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            const auto r_x = gto_coords[i][0];
            const auto r_y = gto_coords[i][1];
            const auto r_z = gto_coords[i][2];

            for (int64_t j = 0; j < npgtos; j++, idx++)
            {
                const auto fexp  = gto_exps[j * ncgtos + i];
                const auto fnorm = gto_norms[j * ncgtos + i];

                gto_info_h[idx]                      = fexp;
                gto_info_h[idx + nrows * npgtos]     = fnorm;
                gto_info_h[idx + nrows * npgtos * 2] = r_x;
                gto_info_h[idx + nrows * npgtos * 3] = r_y;
                gto_info_h[idx + nrows * npgtos * 4] = r_z;
            }
        }
    }

    for (int64_t k = 0; k < ncols; k++)
    {
        grid_xyz_h[k]             = g_x[k];
        grid_xyz_h[k + ncols]     = g_y[k];
        grid_xyz_h[k + ncols * 2] = g_z[k];
    }

    cudaSafe(cudaMemcpy(gto_info_d, gto_info_h, 5 * nrows * npgtos * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(grid_xyz_d, grid_xyz_h, 3 * ncols * sizeof(double), cudaMemcpyHostToDevice));

    dim3 threads_per_block(256);
    dim3 nblocks((ncols + 255) / threads_per_block.x);

    gpu::comp_gto_on_grids<<<nblocks, threads_per_block>>>(
        gto_values_d, gto_info_d, grid_xyz_d, static_cast<uint32_t>(nrows), static_cast<uint32_t>(npgtos), static_cast<uint32_t>(ncols));

    cudaSafe(cudaMemcpy(gto_values_h, gto_values_d, nrows * npgtos * ncols * sizeof(double), cudaMemcpyDeviceToHost));

    for (int64_t i = 0, irow = 0, idx = 0; i < nrows; i++)
    {
        if (gtos_mask[i] == 1)
        {
            for (int64_t j = 0; j < npgtos; j++, idx++)
            {
                for (int64_t k = 0; k < ncols; k++)
                {
                    submat->at(irow, k, false) += gto_values_h[idx * ncols + k];
                }
            }

            irow++;
        }
    }

    cudaSafe(cudaFreeHost(gto_values_h));
    cudaSafe(cudaFreeHost(gto_info_h));
    cudaSafe(cudaFreeHost(grid_xyz_h));

    cudaSafe(cudaFree(gto_values_d));
    cudaSafe(cudaFree(gto_info_d));
    cudaSafe(cudaFree(grid_xyz_d));

    return gto_values;
}

}  // namespace gpu
