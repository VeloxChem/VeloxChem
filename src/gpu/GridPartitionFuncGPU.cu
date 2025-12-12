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

#include <iostream>

#include "ChunkedMemcpyGPU.hpp"
#include "ErrorHandler.hpp"
#include "GridPartitionFuncGPU.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "GpuDevices.hpp"
#include "MathFunc.hpp"

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

__global__ static void
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
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x; // grid

    if (i < nGridPoints)
    {
        const uint32_t idAtomic = atomIdsOfGridPoints[i];

        // grid coordinates

        const double rgx = gridx[i];
        const double rgy = gridy[i];
        const double rgz = gridz[i];

        double pweights_local = 0.0;
        double pweights_sum = 0.0;

        for (uint32_t j = 0; j < nAtoms; j++)
        {
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

            pweights_local += static_cast<double>(j == idAtomic) * pweight_j;
            pweights_sum += pweight_j;
        }

        pweights[i] = pweights_local / pweights_sum;
    }
}

auto
applyGridPartitionFunc(CDenseMatrix*                rawGridPoints,
                       const std::vector<uint32_t>& atomIdsOfGridPoints,
                       const std::vector<double>&   atomMinDistances,
                       const int64_t                nGridPoints,
                       const TPoint3D*              atomCoordinates,
                       const int64_t                nAtoms,
                       const int64_t                numGpusPerNode,
                       const int64_t                rank,
                       const int64_t                nnodes) -> void
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

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

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

    std::vector<double> partial_weights(grid_batch_size, 1.0);

    // grid and atom data on device

    double* d_data_double;
    uint32_t* d_data_uint32;

    auto data_double_count = grid_batch_size * 4 + nAtoms * 2 * 3;
    auto data_uint32_count = grid_batch_size;

    gpuSafe(gpuMalloc(&d_data_double, data_double_count * sizeof(double)));
    gpuSafe(gpuMalloc(&d_data_uint32, data_uint32_count * sizeof(uint32_t)));

    double* d_grid_x = d_data_double;
    double* d_grid_y = d_grid_x + grid_batch_size;
    double* d_grid_z = d_grid_y + grid_batch_size;
    double* d_partial_weights = d_grid_z + grid_batch_size;
    double* d_atom_x = d_partial_weights + grid_batch_size;
    double* d_atom_y = d_atom_x + nAtoms * 2;
    double* d_atom_z = d_atom_y + nAtoms * 2;

    uint32_t* d_atom_ids_of_points = d_data_uint32;

    gpu::chunkedMemcpyHostToDevice<double>(d_grid_x, gridx, grid_batch_size);
    gpu::chunkedMemcpyHostToDevice<double>(d_grid_y, gridy, grid_batch_size);
    gpu::chunkedMemcpyHostToDevice<double>(d_grid_z, gridz, grid_batch_size);

    gpu::chunkedMemcpyHostToDevice<double>(d_atom_x, atom_coords.data() + 0 * nAtoms * 2, nAtoms * 2);
    gpu::chunkedMemcpyHostToDevice<double>(d_atom_y, atom_coords.data() + 1 * nAtoms * 2, nAtoms * 2);
    gpu::chunkedMemcpyHostToDevice<double>(d_atom_z, atom_coords.data() + 2 * nAtoms * 2, nAtoms * 2);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_atom_ids_of_points, atom_ids, grid_batch_size);

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // update weights using GPU

    dim3 threads_per_block(TILE_DIM * TILE_DIM);

    dim3 num_blocks((grid_batch_size + threads_per_block.x - 1) / threads_per_block.x);

    gpu::ssf<<<num_blocks, threads_per_block>>>(
                       d_grid_x, d_grid_y, d_grid_z, d_atom_ids_of_points, static_cast<uint32_t>(grid_batch_size),
                       d_atom_x, d_atom_y, d_atom_z, static_cast<uint32_t>(nAtoms),
                       d_partial_weights);

    gpuSafe(gpuDeviceSynchronize());

    gpu::chunkedMemcpyDeviceToHost<double>(partial_weights.data(), d_partial_weights, grid_batch_size);

#pragma omp barrier

    gpuSafe(gpuFree(d_data_double));
    gpuSafe(gpuFree(d_data_uint32));

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

    }
}

}  // namespace gpu
