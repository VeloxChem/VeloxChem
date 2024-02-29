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
#include "MathFunc.hpp"
#include "MpiFunc.hpp"

#define TILE_DIM 16

#define hipSafe(e)                                                                                                      \
    {                                                                                                                   \
        hipError_t err = (e);                                                                                           \
        if (err != hipSuccess)                                                                                          \
        {                                                                                                               \
            std::cerr << "HIP error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                    \
        }                                                                                                               \
    }

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
    // SSF parameter a = 0.64

    if (eRadius <= -0.64) return -1.0;
    if (eRadius >= 0.64) return 1.0;

    // middle interval

    const double mab = 1.5625 * eRadius;
    const double mab2 = mab * mab;
    return 0.0625 * mab * (35.0 + mab2 * (-35.0 + mab2 * (21.0 - 5.0 * mab2)));
}

__global__ void
ssf(const double*  gridx,
    const double*  gridy,
    const double*  gridz,
    const uint32_t nGridPoints,
    const double*  atomCoordinatesX,
    const double*  atomCoordinatesY,
    const double*  atomCoordinatesZ,
    const uint32_t nAtoms,
    double*        pweights)
{
    const uint32_t i = blockDim.y * blockIdx.y + threadIdx.y;  // grid
    const uint32_t j = blockDim.x * blockIdx.x + threadIdx.x;  // atom

    if ((i < nGridPoints) && (j < nAtoms))
    {
        // grid coordinates

        double rgx = gridx[i];
        double rgy = gridy[i];
        double rgz = gridz[i];

        double rax = atomCoordinatesX[j];
        double ray = atomCoordinatesY[j];
        double raz = atomCoordinatesZ[j];

        double rag = gpu::distance(rax, ray, raz, rgx, rgy, rgz);

        pweights[i * nAtoms + j] = 1.0;

        for (uint32_t k = 0; k < nAtoms; k++)
        {
            if (k != j)
            {
                double rbx = atomCoordinatesX[k];
                double rby = atomCoordinatesY[k];
                double rbz = atomCoordinatesZ[k];

                double rbg = gpu::distance(rbx, rby, rbz, rgx, rgy, rgz);

                double rab = gpu::distance(rax, ray, raz, rbx, rby, rbz);

                double mab = (rag - rbg) / rab;

                pweights[i * nAtoms + j] *= 0.5 * (1.0 - gpu::zeta(mab));
            }
        }
    }
}

auto
applyGridPartitionFunc(CDenseMatrix*   rawGridPoints,
                       const double    minDistance,
                       const int64_t   gridOffset,
                       const int64_t   nGridPoints,
                       const TPoint3D* atomCoordinates,
                       const int64_t   nAtoms,
                       const int64_t   idAtomic,
                       const int64_t   numGpusPerNode) -> void
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

    auto gridx = rawGridPoints->row(0) + gridOffset + grid_batch_offset;
    auto gridy = rawGridPoints->row(1) + gridOffset + grid_batch_offset;
    auto gridz = rawGridPoints->row(2) + gridOffset + grid_batch_offset;
    auto gridw = rawGridPoints->row(3) + gridOffset + grid_batch_offset;

    std::vector<double> atom_coords(3 * nAtoms);

    for (int64_t j = 0; j < nAtoms; j++)
    {
        atom_coords[0 * nAtoms + j] = atomCoordinates[j][0];
        atom_coords[1 * nAtoms + j] = atomCoordinates[j][1];
        atom_coords[2 * nAtoms + j] = atomCoordinates[j][2];
    }

    std::vector<double> partial_weights(grid_batch_size * nAtoms);

    // grid and atom data on device

    double *d_grid_x, *d_grid_y, *d_grid_z;

    hipSafe(hipMalloc(&d_grid_x, grid_batch_size * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_y, grid_batch_size * sizeof(double)));
    hipSafe(hipMalloc(&d_grid_z, grid_batch_size * sizeof(double)));

    hipSafe(hipMemcpy(d_grid_x, gridx, grid_batch_size * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_y, gridy, grid_batch_size * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_grid_z, gridz, grid_batch_size * sizeof(double), hipMemcpyHostToDevice));

    double *d_partial_weights;

    hipSafe(hipMalloc(&d_partial_weights, grid_batch_size * nAtoms * sizeof(double)));

    double *d_atom_x, *d_atom_y, *d_atom_z;

    hipSafe(hipMalloc(&d_atom_x, nAtoms * sizeof(double)));
    hipSafe(hipMalloc(&d_atom_y, nAtoms * sizeof(double)));
    hipSafe(hipMalloc(&d_atom_z, nAtoms * sizeof(double)));

    hipSafe(hipMemcpy(d_atom_x, atom_coords.data() + 0 * nAtoms, nAtoms * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_atom_y, atom_coords.data() + 1 * nAtoms, nAtoms * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_atom_z, atom_coords.data() + 2 * nAtoms, nAtoms * sizeof(double), hipMemcpyHostToDevice));

    // update weights using GPU

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((nAtoms + threads_per_block.x - 1) / threads_per_block.x,
                    (grid_batch_size + threads_per_block.y - 1) / threads_per_block.y);

    hipLaunchKernelGGL(gpu::ssf, num_blocks, threads_per_block, 0, 0,
                       d_grid_x, d_grid_y, d_grid_z, static_cast<uint32_t>(grid_batch_size),
                       d_atom_x, d_atom_y, d_atom_z, static_cast<uint32_t>(nAtoms),
                       d_partial_weights);

    hipSafe(hipMemcpy(partial_weights.data(), d_partial_weights, grid_batch_size * nAtoms * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipDeviceSynchronize());

    hipSafe(hipFree(d_grid_x));
    hipSafe(hipFree(d_grid_y));
    hipSafe(hipFree(d_grid_z));

    hipSafe(hipFree(d_partial_weights));

    hipSafe(hipFree(d_atom_x));
    hipSafe(hipFree(d_atom_y));
    hipSafe(hipFree(d_atom_z));

    for (int64_t i = 0; i < grid_batch_size; i++)
    {
        double rig = mathfunc::distance(atomCoordinates[idAtomic][0],
                                        atomCoordinates[idAtomic][1],
                                        atomCoordinates[idAtomic][2],
                                        gridx[i],
                                        gridy[i],
                                        gridz[i]);

        if (rig < 0.18 * minDistance) continue;

        double wsum = 0.0;

        for (int64_t j = 0; j < nAtoms; j++)
        {
            wsum += partial_weights[i * nAtoms + j];
        }

        gridw[i] *= partial_weights[i * nAtoms + idAtomic] / wsum;
    }

    }}
}

}  // namespace gpu
