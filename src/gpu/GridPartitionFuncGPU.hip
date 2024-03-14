//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#include <omp.h>

#include <iostream>

#include "GridPartitionFuncGPU.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "MathFunc.hpp"
#include "MpiFunc.hpp"

namespace gpu {  // gpu namespace

__device__ double
distance(const double aCoordX, const double aCoordY, const double aCoordZ, const double bCoordX, const double bCoordY, const double bCoordZ)
{
    const double rx = aCoordX - bCoordX;
    const double ry = aCoordY - bCoordY;
    const double rz = aCoordZ - bCoordZ;

    return sqrt(rx * rx + ry * ry + rz * rz);
}

__device__ double
zeta(const double eRadius)
{
    double val = 0.0;

    // SSF parameter a = 0.64

    if (eRadius <= -0.64)
    {
        val = -1.0;
    }
    else if (eRadius >= 0.64)
    {
        val = 1.0;
    }
    else
    {
        const double mab = 1.5625 * eRadius;
        const double mab2 = mab * mab;
        val = 0.0625 * mab * (35.0 + mab2 * (-35.0 + mab2 * (21.0 - 5.0 * mab2)));
    }

    return val;
}

__global__ void
ssf(const double*   gridx,
    const double*   gridy,
    const double*   gridz,
    const uint32_t* atomIdsOfGridPoints,
    const uint32_t  nGridPoints,
    const double*   atomCoordinatesX,
    const double*   atomCoordinatesY,
    const double*   atomCoordinatesZ,
    const uint32_t  nAtoms,
    double*         pweights)
{
    __shared__ double pweights_local[TILE_DIM][TILE_DIM + 1];
    __shared__ double pweights_sum[TILE_DIM][TILE_DIM + 1];

    const uint32_t i = blockDim.y * blockIdx.y + threadIdx.y;  // grid

    if (threadIdx.x == 0)
    {
        pweights_local[threadIdx.y][TILE_DIM] = 0.0;
        pweights_sum[threadIdx.y][TILE_DIM] = 0.0;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (nAtoms + TILE_DIM - 1) / TILE_DIM; m++)
    {
        pweights_local[threadIdx.y][threadIdx.x] = 0.0;
        pweights_sum[threadIdx.y][threadIdx.x] = 0.0;

        __syncthreads();

        const uint32_t j = m * TILE_DIM + threadIdx.x;  // atom

        if ((i < nGridPoints) && (j < nAtoms))
        {
            const uint32_t idAtomic = atomIdsOfGridPoints[i];

            // grid coordinates

            const double rgx = gridx[i];
            const double rgy = gridy[i];
            const double rgz = gridz[i];

            const double rax = atomCoordinatesX[j];
            const double ray = atomCoordinatesY[j];
            const double raz = atomCoordinatesZ[j];

            const double rag = gpu::distance(rax, ray, raz, rgx, rgy, rgz);

            double pweight_j = 1.0;

            for (uint32_t k = j + 1; k < j + nAtoms; k++)
            {
                // Note: atomCoordinatesX/Y/Z have the size of nAtoms*2
                const double rbx = atomCoordinatesX[k];
                const double rby = atomCoordinatesY[k];
                const double rbz = atomCoordinatesZ[k];

                const double rbg = gpu::distance(rbx, rby, rbz, rgx, rgy, rgz);

                const double rab = gpu::distance(rax, ray, raz, rbx, rby, rbz);

                const double mab = (rag - rbg) / rab;

                pweight_j *= 0.5 * (1.0 - gpu::zeta(mab));
            }

            pweights_local[threadIdx.y][threadIdx.x] += static_cast<double>(j == idAtomic) * pweight_j;
            pweights_sum[threadIdx.y][threadIdx.x] += pweight_j;
        }

        __syncthreads();

        if (threadIdx.x == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                pweights_local[threadIdx.y][TILE_DIM] += pweights_local[threadIdx.y][n];
                pweights_sum[threadIdx.y][TILE_DIM]   += pweights_sum[threadIdx.y][n];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.x == 0) && (i < nGridPoints))
    {
        pweights[i] = pweights_local[threadIdx.y][TILE_DIM] / pweights_sum[threadIdx.y][TILE_DIM];
    }
}

auto
applyGridPartitionFunc(CDenseMatrix*                rawGridPoints,
                       const std::vector<uint32_t>& atomIdsOfGridPoints,
                       const std::vector<double>&   atomMinDistances,
                       const int64_t                nGridPoints,
                       const TPoint3D*              atomCoordinates,
                       const int64_t                nAtoms,
                       const int64_t                numGpusPerNode) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = numGpusPerNode;
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    // TODO: make number of GPUs per compute node (not per MPI process) adjustable
    hipSafe(hipSetDevice(gpu_rank % 8));

    // grid and atom data

    auto grid_batch_size = mathfunc::batch_size(nGridPoints, gpu_id, num_gpus_per_node);
    auto grid_batch_offset = mathfunc::batch_offset(nGridPoints, gpu_id, num_gpus_per_node);

    auto gridx = rawGridPoints->row(0) + grid_batch_offset;
    auto gridy = rawGridPoints->row(1) + grid_batch_offset;
    auto gridz = rawGridPoints->row(2) + grid_batch_offset;
    auto gridw = rawGridPoints->row(3) + grid_batch_offset;

    auto atom_ids = atomIdsOfGridPoints.data() + grid_batch_offset;
    auto atom_min_dist = atomMinDistances.data() + grid_batch_offset;

    std::vector<double> atom_coords(3 * (nAtoms * 2));

    for (int64_t j = 0; j < nAtoms; j++)
    {
        atom_coords[0 * nAtoms * 2 + j] = atomCoordinates[j][0];
        atom_coords[1 * nAtoms * 2 + j] = atomCoordinates[j][1];
        atom_coords[2 * nAtoms * 2 + j] = atomCoordinates[j][2];

        atom_coords[0 * nAtoms * 2 + j + nAtoms] = atomCoordinates[j][0];
        atom_coords[1 * nAtoms * 2 + j + nAtoms] = atomCoordinates[j][1];
        atom_coords[2 * nAtoms * 2 + j + nAtoms] = atomCoordinates[j][2];
    }

    std::vector<double> partial_weights(grid_batch_size);

    // grid and atom data on device

    double *d_grid_x, *d_grid_y, *d_grid_z;

    hipSafe(hipMalloc(&d_grid_x, grid_batch_size * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_y, grid_batch_size * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_z, grid_batch_size * sizeof(double)));

    hipSafe(hipMemcpy(d_grid_x, gridx, grid_batch_size * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_y, gridy, grid_batch_size * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_z, gridz, grid_batch_size * sizeof(double), hipMemcpyHostToDevice));

    uint32_t *d_atom_ids_of_points;

    hipSafe(hipMalloc(&d_atom_ids_of_points, grid_batch_size * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_atom_ids_of_points, atom_ids, grid_batch_size * sizeof(uint32_t), hipMemcpyHostToDevice));

    double *d_partial_weights;

    hipSafe(hipMalloc(&d_partial_weights, grid_batch_size * sizeof(double)));

    double *d_atom_x, *d_atom_y, *d_atom_z;

    hipSafe(hipMalloc(&d_atom_x, nAtoms * 2 * sizeof(double)));
    hipSafe(hipMalloc(&d_atom_y, nAtoms * 2 * sizeof(double)));
    hipSafe(hipMalloc(&d_atom_z, nAtoms * 2 * sizeof(double)));

    hipSafe(hipMemcpy(d_atom_x, atom_coords.data() + 0 * nAtoms * 2, nAtoms * 2 * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_atom_y, atom_coords.data() + 1 * nAtoms * 2, nAtoms * 2 * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_atom_z, atom_coords.data() + 2 * nAtoms * 2, nAtoms * 2 * sizeof(double), hipMemcpyHostToDevice));

    // update weights using GPU

    hipSafe(hipDeviceSynchronize());

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks(1, (grid_batch_size + threads_per_block.x - 1) / threads_per_block.x);

    hipLaunchKernelGGL(gpu::ssf, num_blocks, threads_per_block, 0, 0,
                       d_grid_x, d_grid_y, d_grid_z, d_atom_ids_of_points, static_cast<uint32_t>(grid_batch_size),
                       d_atom_x, d_atom_y, d_atom_z, static_cast<uint32_t>(nAtoms),
                       d_partial_weights);

    hipSafe(hipDeviceSynchronize());

    hipSafe(hipMemcpy(partial_weights.data(), d_partial_weights, grid_batch_size * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_grid_x));
    hipSafe(hipFree(d_grid_y));
    hipSafe(hipFree(d_grid_z));

    hipSafe(hipFree(d_atom_ids_of_points));

    hipSafe(hipFree(d_partial_weights));

    hipSafe(hipFree(d_atom_x));
    hipSafe(hipFree(d_atom_y));
    hipSafe(hipFree(d_atom_z));

    for (int64_t i = 0; i < grid_batch_size; i++)
    {
        auto idAtomic = atom_ids[i];

        auto minDistance = atom_min_dist[i];

        double rig = mathfunc::distance(atomCoordinates[idAtomic][0],
                                        atomCoordinates[idAtomic][1],
                                        atomCoordinates[idAtomic][2],
                                        gridx[i],
                                        gridy[i],
                                        gridz[i]);

        if (rig < 0.18 * minDistance) continue;

        gridw[i] *= partial_weights[i];
    }

    }}
}

}  // namespace gpu
