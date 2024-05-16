//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
#include <cusolverDn.h>

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ScreeningData.hpp"
#include "BoysFuncTable.hpp"
#include "DenseLinearAlgebra.hpp"
#include "FockDriverGPU.hpp"
#include "EriCoulomb.hpp"
#include "EriExchange.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuDevices.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MpiFunc.hpp"
#include "MultiTimer.hpp"
#include "OneElectronIntegrals.hpp"
#include "StringFormat.hpp"

namespace gpu {  // gpu namespace

__global__ void
zeroData(double* d_data, const uint32_t n)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < n) d_data[i] = 0.0;
}

auto
computeQMatrixOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> CDenseMatrix
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    const auto all_prim_count = s_prim_count + p_prim_count * 3 + d_prim_count *6;

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> Q_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        Q_omp[gpu_id] = CDenseMatrix(all_prim_count, all_prim_count);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    // TODO: double check by using different number of gpus per MPI
    cudaSafe(cudaSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    cudaSafe(cudaMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), cudaMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    cudaSafe(cudaMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), cudaMemcpyHostToDevice));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d[cart_ind] = sph_ind_coef;
            }
        }
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    double*   d_s_prim_info;

    cudaSafe(cudaMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    cudaSafe(cudaMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    cudaSafe(cudaMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // GTO block pairs

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(gpu_id);
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(gpu_id);
    const auto& sd_first_inds_local = screening.get_sd_first_inds_local(gpu_id);
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(gpu_id);
    const auto& pd_first_inds_local = screening.get_pd_first_inds_local(gpu_id);
    const auto& dd_first_inds_local = screening.get_dd_first_inds_local(gpu_id);

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(gpu_id);
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(gpu_id);
    const auto& sd_second_inds_local = screening.get_sd_second_inds_local(gpu_id);
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(gpu_id);
    const auto& pd_second_inds_local = screening.get_pd_second_inds_local(gpu_id);
    const auto& dd_second_inds_local = screening.get_dd_second_inds_local(gpu_id);

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto sd_prim_pair_count_local = static_cast<int64_t>(sd_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());
    const auto pd_prim_pair_count_local = static_cast<int64_t>(pd_first_inds_local.size());
    const auto dd_prim_pair_count_local = static_cast<int64_t>(dd_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local, sd_prim_pair_count_local,
                                                     pp_prim_pair_count_local, pd_prim_pair_count_local, dd_prim_pair_count_local});

    std::vector<double> h_mat_Q(max_prim_pair_count_local);

    // Q on device

    double *d_mat_Q;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    cudaSafe(cudaMalloc(&d_mat_Q, max_prim_pair_count_local * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaDeviceSynchronize());

    Q_omp[gpu_id].zero();

    auto& mat_Q_omp = Q_omp[gpu_id];

    // compute Q

    // Q: SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSS<<<num_blocks,threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(h_mat_Q.data(), d_mat_Q, ss_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            mat_Q_omp.row(i)[j] = h_mat_Q[ij];

            if (i != j) mat_Q_omp.row(j)[i] = h_mat_Q[ij];
        }
    }

    // Q: SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSP<<<num_blocks,threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(h_mat_Q.data(), d_mat_Q, sp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            mat_Q_omp.row(i)[s_prim_count + j] = h_mat_Q[ij];
            mat_Q_omp.row(s_prim_count + j)[i] = h_mat_Q[ij];
        }
    }

    // Q: SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSD<<<num_blocks,threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(h_mat_Q.data(), d_mat_Q, sd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sd_prim_pair_count_local; ij++)
        {
            const auto i = sd_first_inds_local[ij];
            const auto j = sd_second_inds_local[ij];

            mat_Q_omp.row(i)[s_prim_count + p_prim_count * 3 + j] = h_mat_Q[ij];
            mat_Q_omp.row(s_prim_count + p_prim_count * 3 + j)[i] = h_mat_Q[ij];
        }
    }

    // Q: PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixPP<<<num_blocks,threads_per_block>>>(
                           d_mat_Q,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(h_mat_Q.data(), d_mat_Q, pp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            mat_Q_omp.row(s_prim_count + i)[s_prim_count + j] = h_mat_Q[ij];

            if (i != j) mat_Q_omp.row(s_prim_count + j)[s_prim_count + i] = h_mat_Q[ij];
        }
    }

    // Q: PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixPD<<<num_blocks,threads_per_block>>>(
                           d_mat_Q,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(h_mat_Q.data(), d_mat_Q, pd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pd_prim_pair_count_local; ij++)
        {
            const auto i = pd_first_inds_local[ij];
            const auto j = pd_second_inds_local[ij];

            mat_Q_omp.row(s_prim_count + i)[s_prim_count + p_prim_count * 3 + j] = h_mat_Q[ij];
            mat_Q_omp.row(s_prim_count + p_prim_count * 3 + j)[s_prim_count + i] = h_mat_Q[ij];
        }
    }

    // Q: DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixDD<<<num_blocks,threads_per_block>>>(
                           d_mat_Q,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(h_mat_Q.data(), d_mat_Q, dd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            mat_Q_omp.row(s_prim_count + p_prim_count * 3 + i)[s_prim_count + p_prim_count * 3 + j] = h_mat_Q[ij];

            if (i != j) mat_Q_omp.row(s_prim_count + p_prim_count * 3 + j)[s_prim_count + p_prim_count * 3 + i] = h_mat_Q[ij];
        }
    }

    cudaSafe(cudaDeviceSynchronize());

    cudaSafe(cudaFree(d_boys_func_table));
    cudaSafe(cudaFree(d_boys_func_ft));

    cudaSafe(cudaFree(d_s_prim_info));
    cudaSafe(cudaFree(d_p_prim_info));
    cudaSafe(cudaFree(d_d_prim_info));

    cudaSafe(cudaFree(d_mat_Q));

    cudaSafe(cudaFree(d_ss_first_inds_local));
    cudaSafe(cudaFree(d_ss_second_inds_local));
    cudaSafe(cudaFree(d_sp_first_inds_local));
    cudaSafe(cudaFree(d_sp_second_inds_local));
    cudaSafe(cudaFree(d_sd_first_inds_local));
    cudaSafe(cudaFree(d_sd_second_inds_local));
    cudaSafe(cudaFree(d_pp_first_inds_local));
    cudaSafe(cudaFree(d_pp_second_inds_local));
    cudaSafe(cudaFree(d_pd_first_inds_local));
    cudaSafe(cudaFree(d_pd_second_inds_local));
    cudaSafe(cudaFree(d_dd_first_inds_local));
    cudaSafe(cudaFree(d_dd_second_inds_local));

    }}

    CDenseMatrix Q_matrix_sum(all_prim_count, all_prim_count);

    Q_matrix_sum.zero();

    auto p_mat_Q_sum = Q_matrix_sum.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_Q_omp = Q_omp[gpu_id].values();

        for (int64_t ind = 0; ind < all_prim_count * all_prim_count; ind++)
        {
            p_mat_Q_sum[ind] += p_mat_Q_omp[ind];
        }
    }

    return Q_matrix_sum;
}

auto
computeOverlapAndKineticEnergyIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> S_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> T_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        S_matrices[gpu_id] = CDenseMatrix(naos, naos);
        T_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    cudaSafe(cudaSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d[cart_ind] = sph_ind_coef;
            }
        }
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    double*   d_s_prim_info;

    cudaSafe(cudaMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    cudaSafe(cudaMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    cudaSafe(cudaMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // GTO block pairs

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(gpu_id);
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(gpu_id);
    const auto& sd_first_inds_local = screening.get_sd_first_inds_local(gpu_id);
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(gpu_id);
    const auto& pd_first_inds_local = screening.get_pd_first_inds_local(gpu_id);
    const auto& dd_first_inds_local = screening.get_dd_first_inds_local(gpu_id);

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(gpu_id);
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(gpu_id);
    const auto& sd_second_inds_local = screening.get_sd_second_inds_local(gpu_id);
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(gpu_id);
    const auto& pd_second_inds_local = screening.get_pd_second_inds_local(gpu_id);
    const auto& dd_second_inds_local = screening.get_dd_second_inds_local(gpu_id);

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto sd_prim_pair_count_local = static_cast<int64_t>(sd_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());
    const auto pd_prim_pair_count_local = static_cast<int64_t>(pd_first_inds_local.size());
    const auto dd_prim_pair_count_local = static_cast<int64_t>(dd_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local, sd_prim_pair_count_local,
                                                     pp_prim_pair_count_local, pd_prim_pair_count_local, dd_prim_pair_count_local});

    std::vector<double> mat_S(max_prim_pair_count_local);
    std::vector<double> mat_T(max_prim_pair_count_local);

    // S and T on device

    double *d_mat_S, *d_mat_T;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    cudaSafe(cudaMalloc(&d_mat_S, max_prim_pair_count_local * sizeof(double)));
    cudaSafe(cudaMalloc(&d_mat_T, max_prim_pair_count_local * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    S_matrices[gpu_id].zero();
    T_matrices[gpu_id].zero();

    auto& mat_overlap = S_matrices[gpu_id];
    auto& mat_kinetic_energy = T_matrices[gpu_id];

    cudaSafe(cudaDeviceSynchronize());

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergySS<<<num_blocks,threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local));

        cudaSafe(cudaMemcpy(mat_S.data(), d_mat_S, ss_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));
        cudaSafe(cudaMemcpy(mat_T.data(), d_mat_T, ss_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_overlap.row(i_cgto)[j_cgto] += mat_S[ij];

            if (i != j) mat_overlap.row(j_cgto)[i_cgto] += mat_S[ij];

            mat_kinetic_energy.row(i_cgto)[j_cgto] += mat_T[ij];

            if (i != j) mat_kinetic_energy.row(j_cgto)[i_cgto] += mat_T[ij];
        }
    }

    // SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergySP<<<num_blocks,threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local));

        cudaSafe(cudaMemcpy(mat_S.data(), d_mat_S, sp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));
        cudaSafe(cudaMemcpy(mat_T.data(), d_mat_T, sp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_overlap.row(i_cgto)[j_cgto_sph] += mat_S[ij] * j_coef_sph;
                mat_overlap.row(j_cgto_sph)[i_cgto] += mat_S[ij] * j_coef_sph;

                mat_kinetic_energy.row(i_cgto)[j_cgto_sph] += mat_T[ij] * j_coef_sph;
                mat_kinetic_energy.row(j_cgto_sph)[i_cgto] += mat_T[ij] * j_coef_sph;
            }
        }
    }

    // SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergySD<<<num_blocks,threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local));

        cudaSafe(cudaMemcpy(mat_S.data(), d_mat_S, sd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));
        cudaSafe(cudaMemcpy(mat_T.data(), d_mat_T, sd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sd_prim_pair_count_local; ij++)
        {
            const auto i = sd_first_inds_local[ij];
            const auto j = sd_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_overlap.row(i_cgto)[j_cgto_sph] += mat_S[ij] * j_coef_sph;
                mat_overlap.row(j_cgto_sph)[i_cgto] += mat_S[ij] * j_coef_sph;

                mat_kinetic_energy.row(i_cgto)[j_cgto_sph] += mat_T[ij] * j_coef_sph;
                mat_kinetic_energy.row(j_cgto_sph)[i_cgto] += mat_T[ij] * j_coef_sph;
            }
        }
    }

    // PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergyPP<<<num_blocks,threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local));

        cudaSafe(cudaMemcpy(mat_S.data(), d_mat_S, pp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));
        cudaSafe(cudaMemcpy(mat_T.data(), d_mat_T, pp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(i_cgto_sph)[j_cgto_sph] += mat_S[ij] * coef_sph;

                    if (i != j) mat_overlap.row(j_cgto_sph)[i_cgto_sph] += mat_S[ij] * coef_sph;

                    mat_kinetic_energy.row(i_cgto_sph)[j_cgto_sph] += mat_T[ij] * coef_sph;

                    if (i != j) mat_kinetic_energy.row(j_cgto_sph)[i_cgto_sph] += mat_T[ij] * coef_sph;
                }
            }
        }
    }

    // PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergyPD<<<num_blocks,threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local));

        cudaSafe(cudaMemcpy(mat_S.data(), d_mat_S, pd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));
        cudaSafe(cudaMemcpy(mat_T.data(), d_mat_T, pd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pd_prim_pair_count_local; ij++)
        {
            const auto i = pd_first_inds_local[ij];
            const auto j = pd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(i_cgto_sph)[j_cgto_sph] += mat_S[ij] * coef_sph;
                    mat_overlap.row(j_cgto_sph)[i_cgto_sph] += mat_S[ij] * coef_sph;

                    mat_kinetic_energy.row(i_cgto_sph)[j_cgto_sph] += mat_T[ij] * coef_sph;
                    mat_kinetic_energy.row(j_cgto_sph)[i_cgto_sph] += mat_T[ij] * coef_sph;
                }
            }
        }
    }

    // DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergyDD<<<num_blocks,threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local));

        cudaSafe(cudaMemcpy(mat_S.data(), d_mat_S, dd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));
        cudaSafe(cudaMemcpy(mat_T.data(), d_mat_T, dd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(i_cgto_sph)[j_cgto_sph] += mat_S[ij] * coef_sph;

                    if (i != j) mat_overlap.row(j_cgto_sph)[i_cgto_sph] += mat_S[ij] * coef_sph;

                    mat_kinetic_energy.row(i_cgto_sph)[j_cgto_sph] += mat_T[ij] * coef_sph;

                    if (i != j) mat_kinetic_energy.row(j_cgto_sph)[i_cgto_sph] += mat_T[ij] * coef_sph;
                }
            }
        }
    }

    cudaSafe(cudaDeviceSynchronize());

    cudaSafe(cudaFree(d_s_prim_info));
    cudaSafe(cudaFree(d_p_prim_info));
    cudaSafe(cudaFree(d_d_prim_info));

    cudaSafe(cudaFree(d_mat_S));
    cudaSafe(cudaFree(d_mat_T));

    cudaSafe(cudaFree(d_ss_first_inds_local));
    cudaSafe(cudaFree(d_ss_second_inds_local));
    cudaSafe(cudaFree(d_sp_first_inds_local));
    cudaSafe(cudaFree(d_sp_second_inds_local));
    cudaSafe(cudaFree(d_sd_first_inds_local));
    cudaSafe(cudaFree(d_sd_second_inds_local));
    cudaSafe(cudaFree(d_pp_first_inds_local));
    cudaSafe(cudaFree(d_pp_second_inds_local));
    cudaSafe(cudaFree(d_pd_first_inds_local));
    cudaSafe(cudaFree(d_pd_second_inds_local));
    cudaSafe(cudaFree(d_dd_first_inds_local));
    cudaSafe(cudaFree(d_dd_second_inds_local));

    }}

    std::vector<CDenseMatrix> ST_matrices(2);

    ST_matrices[0] = CDenseMatrix(naos, naos);
    ST_matrices[1] = CDenseMatrix(naos, naos);

    ST_matrices[0].zero();
    ST_matrices[1].zero();

    auto p_mat_S = ST_matrices[0].values();
    auto p_mat_T = ST_matrices[1].values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_overlap = S_matrices[gpu_id].values();
        auto p_mat_kinetic_energy = T_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_S[ind] += p_mat_overlap[ind];
            p_mat_T[ind] += p_mat_kinetic_energy[ind];
        }
    }

    return ST_matrices;
}

auto
computeNuclearPotentialIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> CDenseMatrix
{
    const auto mol_charges = molecule.getCharges();
    const auto mol_coords = molecule.getCoordinates(std::string("BOHR"));
    const auto natoms = molecule.getNumberOfAtoms();

    std::vector<double> points_info(natoms * 4);

    for (int64_t a = 0; a < natoms; a++)
    {
        points_info[a + natoms * 0] = mol_coords[a][0];
        points_info[a + natoms * 1] = mol_coords[a][1];
        points_info[a + natoms * 2] = mol_coords[a][2];
        points_info[a + natoms * 3] = mol_charges[a];
    }

    return gpu::computePointChargesIntegralsOnGPU(molecule, basis, screening, points_info.data(), natoms);
}

auto
computePointChargesIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening, const double* points_info_ptr, const int64_t npoints) -> CDenseMatrix
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> V_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        V_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    cudaSafe(cudaSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    cudaSafe(cudaMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), cudaMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    cudaSafe(cudaMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), cudaMemcpyHostToDevice));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d[cart_ind] = sph_ind_coef;
            }
        }
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    double*   d_s_prim_info;

    cudaSafe(cudaMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    cudaSafe(cudaMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    cudaSafe(cudaMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));

    // GTO block pairs

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(gpu_id);
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(gpu_id);
    const auto& sd_first_inds_local = screening.get_sd_first_inds_local(gpu_id);
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(gpu_id);
    const auto& pd_first_inds_local = screening.get_pd_first_inds_local(gpu_id);
    const auto& dd_first_inds_local = screening.get_dd_first_inds_local(gpu_id);

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(gpu_id);
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(gpu_id);
    const auto& sd_second_inds_local = screening.get_sd_second_inds_local(gpu_id);
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(gpu_id);
    const auto& pd_second_inds_local = screening.get_pd_second_inds_local(gpu_id);
    const auto& dd_second_inds_local = screening.get_dd_second_inds_local(gpu_id);

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto sd_prim_pair_count_local = static_cast<int64_t>(sd_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());
    const auto pd_prim_pair_count_local = static_cast<int64_t>(pd_first_inds_local.size());
    const auto dd_prim_pair_count_local = static_cast<int64_t>(dd_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local, sd_prim_pair_count_local,
                                                     pp_prim_pair_count_local, pd_prim_pair_count_local, dd_prim_pair_count_local});

    std::vector<double> mat_V(max_prim_pair_count_local);

    // V on device

    double *d_mat_V;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    cudaSafe(cudaMalloc(&d_mat_V, max_prim_pair_count_local * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    V_matrices[gpu_id].zero();

    auto& mat_nuclear_potential = V_matrices[gpu_id];

    double *d_points_info;

    cudaSafe(cudaMalloc(&d_points_info, npoints * 4 * sizeof(double)));

    cudaSafe(cudaMemcpy(d_points_info, points_info_ptr, npoints * 4 * sizeof(double), cudaMemcpyHostToDevice));

    cudaSafe(cudaDeviceSynchronize());

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialSS<<<num_blocks,threads_per_block>>>(
                           d_mat_V,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local),
                           d_points_info,
                           static_cast<uint32_t>(npoints),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(mat_V.data(), d_mat_V, ss_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_nuclear_potential.row(i_cgto)[j_cgto] += mat_V[ij];

            if (i != j) mat_nuclear_potential.row(j_cgto)[i_cgto] += mat_V[ij];
        }
    }

    // SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialSP<<<num_blocks,threads_per_block>>>(
                           d_mat_V,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local),
                           d_points_info,
                           static_cast<uint32_t>(npoints),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(mat_V.data(), d_mat_V, sp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_nuclear_potential.row(i_cgto)[j_cgto_sph] += mat_V[ij] * j_coef_sph;
                mat_nuclear_potential.row(j_cgto_sph)[i_cgto] += mat_V[ij] * j_coef_sph;
            }
        }
    }

    // SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialSD<<<num_blocks,threads_per_block>>>(
                           d_mat_V,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local),
                           d_points_info,
                           static_cast<uint32_t>(npoints),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(mat_V.data(), d_mat_V, sd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sd_prim_pair_count_local; ij++)
        {
            const auto i = sd_first_inds_local[ij];
            const auto j = sd_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_nuclear_potential.row(i_cgto)[j_cgto_sph] += mat_V[ij] * j_coef_sph;
                mat_nuclear_potential.row(j_cgto_sph)[i_cgto] += mat_V[ij] * j_coef_sph;
            }
        }
    }

    // PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialPP<<<num_blocks,threads_per_block>>>(
                           d_mat_V,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local),
                           d_points_info,
                           static_cast<uint32_t>(npoints),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(mat_V.data(), d_mat_V, pp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_nuclear_potential.row(i_cgto_sph)[j_cgto_sph] += mat_V[ij] * coef_sph;

                    if (i != j) mat_nuclear_potential.row(j_cgto_sph)[i_cgto_sph] += mat_V[ij] * coef_sph;
                }
            }
        }
    }

    // PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialPD<<<num_blocks,threads_per_block>>>(
                           d_mat_V,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local),
                           d_points_info,
                           static_cast<uint32_t>(npoints),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(mat_V.data(), d_mat_V, pd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pd_prim_pair_count_local; ij++)
        {
            const auto i = pd_first_inds_local[ij];
            const auto j = pd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_nuclear_potential.row(i_cgto_sph)[j_cgto_sph] += mat_V[ij] * coef_sph;
                    mat_nuclear_potential.row(j_cgto_sph)[i_cgto_sph] += mat_V[ij] * coef_sph;
                }
            }
        }
    }

    // DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialDD<<<num_blocks,threads_per_block>>>(
                           d_mat_V,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local),
                           d_points_info,
                           static_cast<uint32_t>(npoints),
                           d_boys_func_table,
                           d_boys_func_ft);

        cudaSafe(cudaMemcpy(mat_V.data(), d_mat_V, dd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_nuclear_potential.row(i_cgto_sph)[j_cgto_sph] += mat_V[ij] * coef_sph;

                    if (i != j) mat_nuclear_potential.row(j_cgto_sph)[i_cgto_sph] += mat_V[ij] * coef_sph;
                }
            }
        }
    }

    cudaSafe(cudaDeviceSynchronize());

    cudaSafe(cudaFree(d_boys_func_table));
    cudaSafe(cudaFree(d_boys_func_ft));

    cudaSafe(cudaFree(d_s_prim_info));
    cudaSafe(cudaFree(d_p_prim_info));
    cudaSafe(cudaFree(d_d_prim_info));

    cudaSafe(cudaFree(d_mat_V));

    cudaSafe(cudaFree(d_ss_first_inds_local));
    cudaSafe(cudaFree(d_ss_second_inds_local));
    cudaSafe(cudaFree(d_sp_first_inds_local));
    cudaSafe(cudaFree(d_sp_second_inds_local));
    cudaSafe(cudaFree(d_sd_first_inds_local));
    cudaSafe(cudaFree(d_sd_second_inds_local));
    cudaSafe(cudaFree(d_pp_first_inds_local));
    cudaSafe(cudaFree(d_pp_second_inds_local));
    cudaSafe(cudaFree(d_pd_first_inds_local));
    cudaSafe(cudaFree(d_pd_second_inds_local));
    cudaSafe(cudaFree(d_dd_first_inds_local));
    cudaSafe(cudaFree(d_dd_second_inds_local));

    cudaSafe(cudaFree(d_points_info));

    }}

    CDenseMatrix V_matrix(naos, naos);

    V_matrix.zero();

    auto p_mat_V = V_matrix.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_nuclear_potential = V_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_V[ind] += p_mat_nuclear_potential[ind];
        }
    }

    return V_matrix;
}

auto
transformDensity(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> CDenseMatrix
{
    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    int64_t s_prim_count = 0, s_ao_count = 0;
    int64_t p_prim_count = 0, p_ao_count = 0;
    int64_t d_prim_count = 0, d_ao_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0)
        {
            s_prim_count += npgtos * ncgtos;
            s_ao_count += ncgtos;
        }
        else if (gto_ang == 1)
        {
            p_prim_count += npgtos * ncgtos;
            p_ao_count += ncgtos;
        }
        else if (gto_ang == 2)
        {
            d_prim_count += npgtos * ncgtos;
            d_ao_count += ncgtos;
        }
    }

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> sph_cart_map;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        auto m = gto_block.getSphericalToCartesianMapping();

        for (const auto& [sph_ind, cart_ind_coef] : m)
        {
            sph_cart_map[sph_ind] = cart_ind_coef;
        }
    }

    const auto cart_naos = s_ao_count + p_ao_count * 3 + d_ao_count * 6; 

    CDenseMatrix cart_dens_mat(cart_naos, cart_naos);
    cart_dens_mat.zero();

    auto cart_dens_ptr = cart_dens_mat.values();
    auto sph_dens_ptr = densityMatrix.alphaDensity(0);

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_map_i, cart_sph_map_j;

    for (int64_t i_cgto = 0; i_cgto < naos; i_cgto++)
    {
        for (int64_t j_cgto = 0; j_cgto < naos; j_cgto++)
        {
            for (const auto& i_cgto_cart_ind_coef : sph_cart_map[i_cgto])
            {
                auto i_cgto_cart = i_cgto_cart_ind_coef.first;
                auto i_coef_cart = i_cgto_cart_ind_coef.second;

                for (const auto& j_cgto_cart_ind_coef : sph_cart_map[j_cgto])
                {
                    auto j_cgto_cart = j_cgto_cart_ind_coef.first;
                    auto j_coef_cart = j_cgto_cart_ind_coef.second;

                    cart_dens_ptr[i_cgto_cart * cart_naos + j_cgto_cart] += 
                        sph_dens_ptr[i_cgto * naos + j_cgto] * i_coef_cart * j_coef_cart;
                }
            }
        }
    }

    return cart_dens_mat;
}

auto
computeFockOnGPU(const              CMolecule& molecule,
                 const              CMolecularBasis& basis,
                 const              CAODensityMatrix& densityMatrix,
                 const double       prefac_coulomb,
                 const double       frac_exact_exchange,
                 const double       omega,
                 const std::string& flag_K,
                 const double       eri_threshold,
                 const double       prelink_threshold,
                 CScreeningData&    screening) -> CDenseMatrix
{
    // TODO sanity check for flag_K: SYMM or ANTISYMM

    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    auto gpu_rank = rank * num_gpus_per_node;

    cudaSafe(cudaSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    CMultiTimer timer;

    timer.start("Prep. blocks");

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    int64_t s_prim_count = 0, s_ao_count = 0;
    int64_t p_prim_count = 0, p_ao_count = 0;
    int64_t d_prim_count = 0, d_ao_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0)
        {
            s_prim_count += npgtos * ncgtos;
            s_ao_count += ncgtos;
        }
        else if (gto_ang == 1)
        {
            p_prim_count += npgtos * ncgtos;
            p_ao_count += ncgtos;
        }
        else if (gto_ang == 2)
        {
            d_prim_count += npgtos * ncgtos;
            d_ao_count += ncgtos;
        }
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d[cart_ind] = sph_ind_coef;
            }
        }
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    timer.stop("Prep. blocks");

    timer.start("Prep. SphToCart");

    // spherical to Cartesian index mapping

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> sph_cart_map;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        auto m = gto_block.getSphericalToCartesianMapping();

        for (const auto& [sph_ind, cart_ind_coef] : m)
        {
            sph_cart_map[sph_ind] = cart_ind_coef;
        }
    }

    const auto cart_naos = s_ao_count + p_ao_count * 3 + d_ao_count * 6; 

    CDenseMatrix cart_dens_mat(cart_naos, cart_naos);
    cart_dens_mat.zero();

    auto cart_dens_ptr = cart_dens_mat.values();
    auto sph_dens_ptr = densityMatrix.alphaDensity(0);

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_map_i, cart_sph_map_j;

    for (int64_t i_cgto = 0; i_cgto < naos; i_cgto++)
    {
        for (int64_t j_cgto = 0; j_cgto < naos; j_cgto++)
        {
            for (const auto& i_cgto_cart_ind_coef : sph_cart_map[i_cgto])
            {
                auto i_cgto_cart = i_cgto_cart_ind_coef.first;
                auto i_coef_cart = i_cgto_cart_ind_coef.second;

                for (const auto& j_cgto_cart_ind_coef : sph_cart_map[j_cgto])
                {
                    auto j_cgto_cart = j_cgto_cart_ind_coef.first;
                    auto j_coef_cart = j_cgto_cart_ind_coef.second;

                    cart_dens_ptr[i_cgto_cart * cart_naos + j_cgto_cart] += 
                        sph_dens_ptr[i_cgto * naos + j_cgto] * i_coef_cart * j_coef_cart;
                }
            }
        }
    }

    timer.stop("Prep. SphToCart");

    timer.start("Prep. sortQD");

    screening.sortQD(s_prim_count, p_prim_count, d_prim_count,
                     s_prim_aoinds, p_prim_aoinds, d_prim_aoinds,
                     s_prim_info, p_prim_info, d_prim_info,
                     cart_naos, cart_dens_ptr, eri_threshold);

    timer.stop("Prep. sortQD");

    CTimer prelink_timer;

    prelink_timer.start();

    timer.start("Prep. Q_prime");

    // preLinK
    // J. Chem. Phys. 138, 134114 (2013)

    // TODO distribute computation of Q_prime

    double *d_matrix_A, *d_matrix_B, *d_matrix_C;

    auto mat_full = screening.get_mat_Q_full(s_prim_count, p_prim_count, d_prim_count);

    cudaSafe(cudaMalloc(&d_matrix_A, mat_full.getNumberOfElements() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_matrix_B, mat_full.getNumberOfElements() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_matrix_C, mat_full.getNumberOfElements() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_matrix_A, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), cudaMemcpyHostToDevice));

    mat_full = screening.get_mat_D_abs_full(s_prim_count, p_prim_count, d_prim_count, s_prim_aoinds, p_prim_aoinds, d_prim_aoinds, cart_naos, cart_dens_ptr);

    cudaSafe(cudaMemcpy(d_matrix_B, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), cudaMemcpyHostToDevice));

    const auto all_prim_count = mat_full.getNumberOfRows();

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto n = static_cast<int32_t>(all_prim_count);

    // compute A^T * (B^T * A^T) since cublas is column-major
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha, d_matrix_B, n, d_matrix_A, n, &beta, d_matrix_C, n));

    cudaSafe(cudaDeviceSynchronize());

    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha, d_matrix_A, n, d_matrix_C, n, &beta, d_matrix_B, n));

    cudaSafe(cudaMemcpy(mat_full.values(), d_matrix_B, mat_full.getNumberOfElements() * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_matrix_A));
    cudaSafe(cudaFree(d_matrix_B));
    cudaSafe(cudaFree(d_matrix_C));

    timer.stop("Prep. Q_prime");

    prelink_timer.stop();

    auto prelink_elapsed_time = prelink_timer.getElapsedTime();

    screening.setPreLinkTime(prelink_elapsed_time);

    timer.start("Prep. preLinK");

    screening.form_pair_inds_for_K(s_prim_count, p_prim_count, d_prim_count, mat_full, prelink_threshold);

    mat_full = CDenseMatrix();

    // max densities are needed by exchange Fock
    screening.findMaxDensities(s_prim_count, p_prim_count, d_prim_count,
                               s_prim_aoinds, p_prim_aoinds, d_prim_aoinds,
                               cart_naos, cart_dens_ptr);

    timer.stop("Prep. preLinK");

    std::string errnaos("gpu::computeFockGPU: Inconsistent number of AOs");
    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    std::vector<CDenseMatrix> mat_Fock_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Fock_omp[gpu_id] = CDenseMatrix(naos, naos);
    }

    screening.initTimers(num_gpus_per_node);

    // std::cout << "\nTiming of prep. on rank " << rank << "\n";
    // std::cout << "-------------------------\n";
    // std::cout << timer.getSummary() << std::endl;

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    cudaSafe(cudaSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Boys func. prep.");

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    cudaSafe(cudaMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), cudaMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    cudaSafe(cudaMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), cudaMemcpyHostToDevice));

    timer.stop("Boys func. prep.");

    timer.start("GTO block prep.");

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    double*   d_s_prim_info;
    uint32_t* d_s_prim_aoinds;

    cudaSafe(cudaMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_s_prim_aoinds, s_prim_aoinds.size() * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;
    uint32_t* d_p_prim_aoinds;

    cudaSafe(cudaMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_p_prim_aoinds, p_prim_aoinds.size() * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    double*   d_d_prim_info;
    uint32_t* d_d_prim_aoinds;

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    cudaSafe(cudaMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_d_prim_aoinds, d_prim_aoinds.size() * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_d_prim_aoinds, d_prim_aoinds.data(), d_prim_aoinds.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    timer.stop("GTO block prep.");

    // GTO block pairs

    timer.start("Coulomb prep.");

    const auto ss_mat_Q_orig = screening.getQMatrixSS();
    const auto sp_mat_Q_orig = screening.getQMatrixSP();
    const auto sd_mat_Q_orig = screening.getQMatrixSD();
    const auto pp_mat_Q_orig = screening.getQMatrixPP();
    const auto pd_mat_Q_orig = screening.getQMatrixPD();
    const auto dd_mat_Q_orig = screening.getQMatrixDD();

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(gpu_id);
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(gpu_id);
    const auto& sd_first_inds_local = screening.get_sd_first_inds_local(gpu_id);
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(gpu_id);
    const auto& pd_first_inds_local = screening.get_pd_first_inds_local(gpu_id);
    const auto& dd_first_inds_local = screening.get_dd_first_inds_local(gpu_id);

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(gpu_id);
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(gpu_id);
    const auto& sd_second_inds_local = screening.get_sd_second_inds_local(gpu_id);
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(gpu_id);
    const auto& pd_second_inds_local = screening.get_pd_second_inds_local(gpu_id);
    const auto& dd_second_inds_local = screening.get_dd_second_inds_local(gpu_id);

    const auto& ss_mat_Q_local = screening.get_ss_mat_Q_local(gpu_id);
    const auto& sp_mat_Q_local = screening.get_sp_mat_Q_local(gpu_id);
    const auto& sd_mat_Q_local = screening.get_sd_mat_Q_local(gpu_id);
    const auto& pp_mat_Q_local = screening.get_pp_mat_Q_local(gpu_id);
    const auto& pd_mat_Q_local = screening.get_pd_mat_Q_local(gpu_id);
    const auto& dd_mat_Q_local = screening.get_dd_mat_Q_local(gpu_id);

    const auto& ss_pair_data_local = screening.get_ss_pair_data_local(gpu_id);
    const auto& sp_pair_data_local = screening.get_sp_pair_data_local(gpu_id);
    const auto& sd_pair_data_local = screening.get_sd_pair_data_local(gpu_id);
    const auto& pp_pair_data_local = screening.get_pp_pair_data_local(gpu_id);
    const auto& pd_pair_data_local = screening.get_pd_pair_data_local(gpu_id);
    const auto& dd_pair_data_local = screening.get_dd_pair_data_local(gpu_id);

    const auto& ss_first_inds = screening.get_ss_first_inds();
    const auto& sp_first_inds = screening.get_sp_first_inds();
    const auto& sd_first_inds = screening.get_sd_first_inds();
    const auto& pp_first_inds = screening.get_pp_first_inds();
    const auto& pd_first_inds = screening.get_pd_first_inds();
    const auto& dd_first_inds = screening.get_dd_first_inds();

    const auto& ss_second_inds = screening.get_ss_second_inds();
    const auto& sp_second_inds = screening.get_sp_second_inds();
    const auto& sd_second_inds = screening.get_sd_second_inds();
    const auto& pp_second_inds = screening.get_pp_second_inds();
    const auto& pd_second_inds = screening.get_pd_second_inds();
    const auto& dd_second_inds = screening.get_dd_second_inds();

    const auto& ss_mat_Q = screening.get_ss_mat_Q(); 
    const auto& sp_mat_Q = screening.get_sp_mat_Q(); 
    const auto& sd_mat_Q = screening.get_sd_mat_Q(); 
    const auto& pp_mat_Q = screening.get_pp_mat_Q(); 
    const auto& pd_mat_Q = screening.get_pd_mat_Q(); 
    const auto& dd_mat_Q = screening.get_dd_mat_Q(); 

    const auto& ss_mat_D = screening.get_ss_mat_D(); 
    const auto& sp_mat_D = screening.get_sp_mat_D(); 
    const auto& sd_mat_D = screening.get_sd_mat_D(); 
    const auto& pp_mat_D = screening.get_pp_mat_D(); 
    const auto& pd_mat_D = screening.get_pd_mat_D(); 
    const auto& dd_mat_D = screening.get_dd_mat_D(); 

    const auto& ss_pair_data = screening.get_ss_pair_data(); 
    const auto& sp_pair_data = screening.get_sp_pair_data(); 
    const auto& sd_pair_data = screening.get_sd_pair_data(); 
    const auto& pp_pair_data = screening.get_pp_pair_data(); 
    const auto& pd_pair_data = screening.get_pd_pair_data(); 
    const auto& dd_pair_data = screening.get_dd_pair_data(); 

    const auto ss_prim_pair_count = static_cast<int64_t>(ss_first_inds.size());
    const auto sp_prim_pair_count = static_cast<int64_t>(sp_first_inds.size());
    const auto sd_prim_pair_count = static_cast<int64_t>(sd_first_inds.size());
    const auto pp_prim_pair_count = static_cast<int64_t>(pp_first_inds.size());
    const auto pd_prim_pair_count = static_cast<int64_t>(pd_first_inds.size());
    const auto dd_prim_pair_count = static_cast<int64_t>(dd_first_inds.size());

    const auto max_prim_pair_count = std::max({ss_prim_pair_count, sp_prim_pair_count, sd_prim_pair_count,
                                               pp_prim_pair_count, pd_prim_pair_count, dd_prim_pair_count});

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto sd_prim_pair_count_local = static_cast<int64_t>(sd_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());
    const auto pd_prim_pair_count_local = static_cast<int64_t>(pd_first_inds_local.size());
    const auto dd_prim_pair_count_local = static_cast<int64_t>(dd_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local,
                                                     sd_prim_pair_count_local, pp_prim_pair_count_local,
                                                     pd_prim_pair_count_local, dd_prim_pair_count_local});

    std::vector<double> mat_J(max_prim_pair_count_local);

    // sorted Q, D, and indices on device

    double *d_mat_J, *d_mat_D;
    double *d_ss_mat_Q, *d_sp_mat_Q, *d_sd_mat_Q, *d_pp_mat_Q, *d_pd_mat_Q, *d_dd_mat_Q;

    uint32_t *d_ss_first_inds, *d_ss_second_inds;
    uint32_t *d_sp_first_inds, *d_sp_second_inds;
    uint32_t *d_sd_first_inds, *d_sd_second_inds;
    uint32_t *d_pp_first_inds, *d_pp_second_inds;
    uint32_t *d_pd_first_inds, *d_pd_second_inds;
    uint32_t *d_dd_first_inds, *d_dd_second_inds;

    double *d_ss_pair_data, *d_sp_pair_data, *d_sd_pair_data,
           *d_pp_pair_data, *d_pd_pair_data, *d_dd_pair_data;

    double *d_ss_mat_Q_local, *d_sp_mat_Q_local, *d_sd_mat_Q_local, *d_pp_mat_Q_local, *d_pd_mat_Q_local, *d_dd_mat_Q_local;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    double *d_ss_pair_data_local, *d_sp_pair_data_local, *d_sd_pair_data_local,
           *d_pp_pair_data_local, *d_pd_pair_data_local, *d_dd_pair_data_local;

    cudaSafe(cudaMalloc(&d_mat_D, max_prim_pair_count * sizeof(double)));
    cudaSafe(cudaMalloc(&d_mat_J, max_prim_pair_count_local * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_mat_Q, ss_prim_pair_count * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sp_mat_Q, sp_prim_pair_count * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sd_mat_Q, sd_prim_pair_count * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pp_mat_Q, pp_prim_pair_count * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pd_mat_Q, pd_prim_pair_count * sizeof(double)));
    cudaSafe(cudaMalloc(&d_dd_mat_Q, dd_prim_pair_count * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_first_inds, ss_prim_pair_count * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_ss_second_inds, ss_prim_pair_count * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sp_first_inds, sp_prim_pair_count * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sp_second_inds, sp_prim_pair_count * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sd_first_inds, sd_prim_pair_count * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sd_second_inds, sd_prim_pair_count * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pp_first_inds, pp_prim_pair_count * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pp_second_inds, pp_prim_pair_count * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pd_first_inds, pd_prim_pair_count * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pd_second_inds, pd_prim_pair_count * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_dd_first_inds, dd_prim_pair_count * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_dd_second_inds, dd_prim_pair_count * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_ss_pair_data, ss_pair_data.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sp_pair_data, sp_pair_data.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sd_pair_data, sd_pair_data.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pp_pair_data, pp_pair_data.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pd_pair_data, pd_pair_data.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_dd_pair_data, dd_pair_data.size() * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_mat_Q_local, ss_prim_pair_count_local * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sp_mat_Q_local, sp_prim_pair_count_local * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sd_mat_Q_local, sd_prim_pair_count_local * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pp_mat_Q_local, pp_prim_pair_count_local * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pd_mat_Q_local, pd_prim_pair_count_local * sizeof(double)));
    cudaSafe(cudaMalloc(&d_dd_mat_Q_local, dd_prim_pair_count_local * sizeof(double)));

    cudaSafe(cudaMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_ss_pair_data_local, ss_pair_data_local.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sp_pair_data_local, sp_pair_data_local.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_sd_pair_data_local, sd_pair_data_local.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pp_pair_data_local, pp_pair_data_local.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pd_pair_data_local, pd_pair_data_local.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_dd_pair_data_local, dd_pair_data_local.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_ss_mat_Q, ss_mat_Q.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_mat_Q, sp_mat_Q.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_mat_Q, sd_mat_Q.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_mat_Q, pp_mat_Q.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_mat_Q, pd_mat_Q.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_mat_Q, dd_mat_Q.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_ss_first_inds, ss_first_inds.data(), ss_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_ss_second_inds, ss_second_inds.data(), ss_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sp_first_inds, sp_first_inds.data(), sp_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_second_inds, sp_second_inds.data(), sp_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sd_first_inds, sd_first_inds.data(), sd_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_second_inds, sd_second_inds.data(), sd_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pp_first_inds, pp_first_inds.data(), pp_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_second_inds, pp_second_inds.data(), pp_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pd_first_inds, pd_first_inds.data(), pd_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_second_inds, pd_second_inds.data(), pd_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_dd_first_inds, dd_first_inds.data(), dd_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_second_inds, dd_second_inds.data(), dd_prim_pair_count * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_ss_pair_data, ss_pair_data.data(), ss_pair_data.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_pair_data, sp_pair_data.data(), sp_pair_data.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_pair_data, sd_pair_data.data(), sd_pair_data.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_pair_data, pp_pair_data.data(), pp_pair_data.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_pair_data, pd_pair_data.data(), pd_pair_data.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_pair_data, dd_pair_data.data(), dd_pair_data.size() * sizeof(double), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_ss_mat_Q_local, ss_mat_Q_local.data(), ss_prim_pair_count_local * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_mat_Q_local, sp_mat_Q_local.data(), sp_prim_pair_count_local * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_mat_Q_local, sd_mat_Q_local.data(), sd_prim_pair_count_local * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_mat_Q_local, pp_mat_Q_local.data(), pp_prim_pair_count_local * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_mat_Q_local, pd_mat_Q_local.data(), pd_prim_pair_count_local * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_mat_Q_local, dd_mat_Q_local.data(), dd_prim_pair_count_local * sizeof(double), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_ss_pair_data_local, ss_pair_data_local.data(), ss_pair_data_local.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sp_pair_data_local, sp_pair_data_local.data(), sp_pair_data_local.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_sd_pair_data_local, sd_pair_data_local.data(), sd_pair_data_local.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pp_pair_data_local, pp_pair_data_local.data(), pp_pair_data_local.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pd_pair_data_local, pd_pair_data_local.data(), pd_pair_data_local.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_dd_pair_data_local, dd_pair_data_local.data(), dd_pair_data_local.size() * sizeof(double), cudaMemcpyHostToDevice));

    timer.stop("Coulomb prep.");

    mat_Fock_omp[gpu_id].zero();

    cudaSafe(cudaDeviceSynchronize());

    CTimer coulomb_timer;

    coulomb_timer.start();

    timer.start("J computation");

    // compute J

    if (std::fabs(prefac_coulomb) > 1.0e-13)
    {

    // J: S-S block

    if (ss_prim_pair_count_local > 0)
    {
        timer.start("  J block SS");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_J,static_cast<uint32_t>(ss_prim_pair_count_local));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SS|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            timer.start("    J block SSSS");

            gpu::computeCoulombFockSSSS<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_mat_D,
                               d_ss_mat_Q_local,
                               d_ss_mat_Q,
                               d_ss_first_inds_local,
                               d_ss_second_inds_local,
                               d_ss_pair_data_local,
                               static_cast<uint32_t>(ss_prim_pair_count_local),
                               d_ss_first_inds,
                               d_ss_second_inds,
                               d_ss_pair_data,
                               static_cast<uint32_t>(ss_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block SSSS");
        }

        // J: (SS|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSSSP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_ss_mat_Q_local,
                               d_sp_mat_Q,
                               d_ss_first_inds_local,
                               d_ss_second_inds_local,
                               d_ss_pair_data_local,
                               static_cast<uint32_t>(ss_prim_pair_count_local),
                               d_sp_first_inds,
                               d_sp_second_inds,
                               d_sp_pair_data,
                               static_cast<uint32_t>(sp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SS|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSSSD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_ss_mat_Q_local,
                               d_sd_mat_Q,
                               d_ss_first_inds_local,
                               d_ss_second_inds_local,
                               d_ss_pair_data_local,
                               static_cast<uint32_t>(ss_prim_pair_count_local),
                               d_sd_first_inds,
                               d_sd_second_inds,
                               d_sd_pair_data,
                               static_cast<uint32_t>(sd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SS|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSSPP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_ss_mat_Q_local,
                               d_pp_mat_Q,
                               d_ss_first_inds_local,
                               d_ss_second_inds_local,
                               d_ss_pair_data_local,
                               static_cast<uint32_t>(ss_prim_pair_count_local),
                               d_pp_first_inds,
                               d_pp_second_inds,
                               d_pp_pair_data,
                               static_cast<uint32_t>(pp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SS|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSSPD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_ss_mat_Q_local,
                               d_pd_mat_Q,
                               d_ss_first_inds_local,
                               d_ss_second_inds_local,
                               d_ss_pair_data_local,
                               static_cast<uint32_t>(ss_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SS|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSSDD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_ss_mat_Q_local,
                               d_dd_mat_Q,
                               d_ss_first_inds_local,
                               d_ss_second_inds_local,
                               d_ss_pair_data_local,
                               static_cast<uint32_t>(ss_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        cudaSafe(cudaMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * prefac_coulomb;

            if (i != j) mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * prefac_coulomb;
        }

        timer.stop("  J block SS");
    }

    // J: S-P block

    if (sp_prim_pair_count_local > 0)
    {
        timer.start("  J block SP");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_J,static_cast<uint32_t>(sp_prim_pair_count_local));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSPSS<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_sp_mat_Q_local,
                               d_ss_mat_Q,
                               d_sp_first_inds_local,
                               d_sp_second_inds_local,
                               d_sp_pair_data_local,
                               static_cast<uint32_t>(sp_prim_pair_count_local),
                               d_ss_first_inds,
                               d_ss_second_inds,
                               d_ss_pair_data,
                               static_cast<uint32_t>(ss_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSPSP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_sp_mat_Q_local,
                               d_sp_mat_Q,
                               d_sp_first_inds_local,
                               d_sp_second_inds_local,
                               d_sp_pair_data_local,
                               static_cast<uint32_t>(sp_prim_pair_count_local),
                               d_sp_first_inds,
                               d_sp_second_inds,
                               d_sp_pair_data,
                               static_cast<uint32_t>(sp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSPSD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sp_mat_Q_local,
                               d_sd_mat_Q,
                               d_sp_first_inds_local,
                               d_sp_second_inds_local,
                               d_sp_pair_data_local,
                               static_cast<uint32_t>(sp_prim_pair_count_local),
                               d_sd_first_inds,
                               d_sd_second_inds,
                               d_sd_pair_data,
                               static_cast<uint32_t>(sd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSPPP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_sp_mat_Q_local,
                               d_pp_mat_Q,
                               d_sp_first_inds_local,
                               d_sp_second_inds_local,
                               d_sp_pair_data_local,
                               static_cast<uint32_t>(sp_prim_pair_count_local),
                               d_pp_first_inds,
                               d_pp_second_inds,
                               d_pp_pair_data,
                               static_cast<uint32_t>(pp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSPPD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sp_mat_Q_local,
                               d_pd_mat_Q,
                               d_sp_first_inds_local,
                               d_sp_second_inds_local,
                               d_sp_pair_data_local,
                               static_cast<uint32_t>(sp_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSPDD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sp_mat_Q_local,
                               d_dd_mat_Q,
                               d_sp_first_inds_local,
                               d_sp_second_inds_local,
                               d_sp_pair_data_local,
                               static_cast<uint32_t>(sp_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        cudaSafe(cudaMemcpy(mat_J.data(), d_mat_J, sp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto_sph] += mat_J[ij] * j_coef_sph * prefac_coulomb;
                mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto] += mat_J[ij] * j_coef_sph * prefac_coulomb;
            }
        }

        timer.stop("  J block SP");
    }

    // J: P-P block

    if (pp_prim_pair_count_local > 0)
    {
        timer.start("  J block PP");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_J,static_cast<uint32_t>(pp_prim_pair_count_local));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPPSS<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_pp_mat_Q_local,
                               d_ss_mat_Q,
                               d_pp_first_inds_local,
                               d_pp_second_inds_local,
                               d_pp_pair_data_local,
                               static_cast<uint32_t>(pp_prim_pair_count_local),
                               d_ss_first_inds,
                               d_ss_second_inds,
                               d_ss_pair_data,
                               static_cast<uint32_t>(ss_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPPSP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_pp_mat_Q_local,
                               d_sp_mat_Q,
                               d_pp_first_inds_local,
                               d_pp_second_inds_local,
                               d_pp_pair_data_local,
                               static_cast<uint32_t>(pp_prim_pair_count_local),
                               d_sp_first_inds,
                               d_sp_second_inds,
                               d_sp_pair_data,
                               static_cast<uint32_t>(sp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPPSD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pp_mat_Q_local,
                               d_sd_mat_Q,
                               d_pp_first_inds_local,
                               d_pp_second_inds_local,
                               d_pp_pair_data_local,
                               static_cast<uint32_t>(pp_prim_pair_count_local),
                               d_sd_first_inds,
                               d_sd_second_inds,
                               d_sd_pair_data,
                               static_cast<uint32_t>(sd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPPPP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D,
                               d_pp_mat_Q_local,
                               d_pp_mat_Q,
                               d_pp_first_inds_local,
                               d_pp_second_inds_local,
                               d_pp_pair_data_local,
                               static_cast<uint32_t>(pp_prim_pair_count_local),
                               d_pp_first_inds,
                               d_pp_second_inds,
                               d_pp_pair_data,
                               static_cast<uint32_t>(pp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPPPD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pp_mat_Q_local,
                               d_pd_mat_Q,
                               d_pp_first_inds_local,
                               d_pp_second_inds_local,
                               d_pp_pair_data_local,
                               static_cast<uint32_t>(pp_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block PPDD");

            cudaSafe(cudaMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPPDD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pp_mat_Q_local,
                               d_dd_mat_Q,
                               d_pp_first_inds_local,
                               d_pp_second_inds_local,
                               d_pp_pair_data_local,
                               static_cast<uint32_t>(pp_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block PPDD");
        }

        cudaSafe(cudaMemcpy(mat_J.data(), d_mat_J, pp_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[j_cgto_sph] += mat_J[ij] * coef_sph * prefac_coulomb;

                    if (i != j) mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto_sph] += mat_J[ij] * coef_sph * prefac_coulomb;
                }
            }
        }

        timer.stop("  J block PP");
    }

    // J: S-D block

    if (sd_prim_pair_count_local > 0)
    {
        timer.start("  J block SD");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_J,static_cast<uint32_t>(sd_prim_pair_count_local));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSDSS<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sd_mat_Q_local,
                               d_ss_mat_Q,
                               d_sd_first_inds_local,
                               d_sd_second_inds_local,
                               d_sd_pair_data_local,
                               static_cast<uint32_t>(sd_prim_pair_count_local),
                               d_ss_first_inds,
                               d_ss_second_inds,
                               d_ss_pair_data,
                               static_cast<uint32_t>(ss_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSDSP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sd_mat_Q_local,
                               d_sp_mat_Q,
                               d_sd_first_inds_local,
                               d_sd_second_inds_local,
                               d_sd_pair_data_local,
                               static_cast<uint32_t>(sd_prim_pair_count_local),
                               d_sp_first_inds,
                               d_sp_second_inds,
                               d_sp_pair_data,
                               static_cast<uint32_t>(sp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSDSD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sd_mat_Q_local,
                               d_sd_mat_Q,
                               d_sd_first_inds_local,
                               d_sd_second_inds_local,
                               d_sd_pair_data_local,
                               static_cast<uint32_t>(sd_prim_pair_count_local),
                               d_sd_first_inds,
                               d_sd_second_inds,
                               d_sd_pair_data,
                               static_cast<uint32_t>(sd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSDPP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sd_mat_Q_local,
                               d_pp_mat_Q,
                               d_sd_first_inds_local,
                               d_sd_second_inds_local,
                               d_sd_pair_data_local,
                               static_cast<uint32_t>(sd_prim_pair_count_local),
                               d_pp_first_inds,
                               d_pp_second_inds,
                               d_pp_pair_data,
                               static_cast<uint32_t>(pp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSDPD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sd_mat_Q_local,
                               d_pd_mat_Q,
                               d_sd_first_inds_local,
                               d_sd_second_inds_local,
                               d_sd_pair_data_local,
                               static_cast<uint32_t>(sd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (SD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block SDDD");

            cudaSafe(cudaMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockSDDD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_sd_mat_Q_local,
                               d_dd_mat_Q,
                               d_sd_first_inds_local,
                               d_sd_second_inds_local,
                               d_sd_pair_data_local,
                               static_cast<uint32_t>(sd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block SDDD");
        }

        cudaSafe(cudaMemcpy(mat_J.data(), d_mat_J, sd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sd_prim_pair_count_local; ij++)
        {
            const auto i = sd_first_inds_local[ij];
            const auto j = sd_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto_sph] += mat_J[ij] * j_coef_sph * prefac_coulomb;
                mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto] += mat_J[ij] * j_coef_sph * prefac_coulomb;
            }
        }

        timer.stop("  J block SD");
    }

    // J: P-D block

    if (pd_prim_pair_count_local > 0)
    {
        timer.start("  J block PD");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_J,static_cast<uint32_t>(pd_prim_pair_count_local));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPDSS<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_ss_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_ss_first_inds,
                               d_ss_second_inds,
                               d_ss_pair_data,
                               static_cast<uint32_t>(ss_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            cudaSafe(cudaMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPDSP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_sp_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_sp_first_inds,
                               d_sp_second_inds,
                               d_sp_pair_data,
                               static_cast<uint32_t>(sp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());
        }

        // J: (PD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            timer.start("    J block PDSD");

            cudaSafe(cudaMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPDSD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_sd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_sd_first_inds,
                               d_sd_second_inds,
                               d_sd_pair_data,
                               static_cast<uint32_t>(sd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block PDSD");
        }

        // J: (PD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            timer.start("    J block PDPP");

            cudaSafe(cudaMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPDPP<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_pp_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_pp_first_inds,
                               d_pp_second_inds,
                               d_pp_pair_data,
                               static_cast<uint32_t>(pp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block PDPP");
        }

        // J: (PD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            timer.start("    J block PDPD");

            cudaSafe(cudaMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPDPD<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_pd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block PDPD");
        }

        // J: (PD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block PDDD");

            cudaSafe(cudaMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            gpu::computeCoulombFockPDDD0<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_dd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockPDDD1<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_dd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockPDDD2<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_dd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockPDDD3<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_dd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockPDDD4<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_dd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockPDDD5<<<num_blocks,threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_pd_mat_Q_local,
                               d_dd_mat_Q,
                               d_pd_first_inds_local,
                               d_pd_second_inds_local,
                               d_pd_pair_data_local,
                               static_cast<uint32_t>(pd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block PDDD");
        }

        cudaSafe(cudaMemcpy(mat_J.data(), d_mat_J, pd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pd_prim_pair_count_local; ij++)
        {
            const auto i = pd_first_inds_local[ij];
            const auto j = pd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[j_cgto_sph] += mat_J[ij] * coef_sph * prefac_coulomb;
                    mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto_sph] += mat_J[ij] * coef_sph * prefac_coulomb;
                }
            }
        }

        timer.stop("  J block PD");
    }

    // J: D-D block

    if (dd_prim_pair_count_local > 0)
    {
        timer.start("  J block DD");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_J,static_cast<uint32_t>(dd_prim_pair_count_local));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (DD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            timer.start("    J block DDSS");

            cudaSafe(cudaMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDSS<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_ss_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_ss_first_inds,
                               d_ss_second_inds,
                               d_ss_pair_data,
                               static_cast<uint32_t>(ss_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block DDSS");
        }

        // J: (DD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            timer.start("    J block DDSP");

            cudaSafe(cudaMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDSP<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_sp_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_sp_first_inds,
                               d_sp_second_inds,
                               d_sp_pair_data,
                               static_cast<uint32_t>(sp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block DDSP");
        }

        // J: (DD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            timer.start("    J block DDSD");

            cudaSafe(cudaMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDSD<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_sd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_sd_first_inds,
                               d_sd_second_inds,
                               d_sd_pair_data,
                               static_cast<uint32_t>(sd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block DDSD");
        }

        // J: (DD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            timer.start("    J block DDPP");

            cudaSafe(cudaMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDPP<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pp_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pp_first_inds,
                               d_pp_second_inds,
                               d_pp_pair_data,
                               static_cast<uint32_t>(pp_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block DDPP");
        }

        // J: (DD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            timer.start("    J block DDPD");

            cudaSafe(cudaMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDPD0<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD1<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD2<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD3<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD4<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD5<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD6<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD7<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD8<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDPD9<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_pd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_pd_first_inds,
                               d_pd_second_inds,
                               d_pd_pair_data,
                               static_cast<uint32_t>(pd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block DDPD");
        }

        // J: (DD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block DDDD");

            cudaSafe(cudaMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), cudaMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDDD0<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD1<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD2<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD3<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD4<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD5<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD6<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD7<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD8<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD9<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD10<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD11<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD12<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD13<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD14<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD15<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD16<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD17<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD18<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD19<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD20<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD21<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD22<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            gpu::computeCoulombFockDDDD23<<<dd_num_blocks,dd_threads_per_block>>>(
                               d_mat_J,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D,
                               d_dd_mat_Q_local,
                               d_dd_mat_Q,
                               d_dd_first_inds_local,
                               d_dd_second_inds_local,
                               d_dd_pair_data_local,
                               static_cast<uint32_t>(dd_prim_pair_count_local),
                               d_dd_first_inds,
                               d_dd_second_inds,
                               d_dd_pair_data,
                               static_cast<uint32_t>(dd_prim_pair_count),
                               d_boys_func_table,
                               d_boys_func_ft,
                               eri_threshold);

            cudaSafe(cudaDeviceSynchronize());

            timer.stop("    J block DDDD");
        }

        cudaSafe(cudaMemcpy(mat_J.data(), d_mat_J, dd_prim_pair_count_local * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[j_cgto_sph] += mat_J[ij] * coef_sph * prefac_coulomb;

                    if (i != j) mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto_sph] += mat_J[ij] * coef_sph * prefac_coulomb;
                }
            }
        }

        timer.stop("  J block DD");
    }

    }  // end of compute J

    timer.stop("J computation");

    cudaSafe(cudaDeviceSynchronize());

    coulomb_timer.stop();

    auto coulomb_elapsed_time = coulomb_timer.getElapsedTime();

    screening.setCoulombTime(gpu_id, coulomb_elapsed_time);

    cudaSafe(cudaFree(d_mat_D));
    cudaSafe(cudaFree(d_mat_J));

    cudaSafe(cudaFree(d_ss_mat_Q));
    cudaSafe(cudaFree(d_sp_mat_Q));
    cudaSafe(cudaFree(d_sd_mat_Q));
    cudaSafe(cudaFree(d_pp_mat_Q));
    cudaSafe(cudaFree(d_pd_mat_Q));
    cudaSafe(cudaFree(d_dd_mat_Q));

    cudaSafe(cudaFree(d_ss_first_inds));
    cudaSafe(cudaFree(d_ss_second_inds));
    cudaSafe(cudaFree(d_sp_first_inds));
    cudaSafe(cudaFree(d_sp_second_inds));
    cudaSafe(cudaFree(d_sd_first_inds));
    cudaSafe(cudaFree(d_sd_second_inds));
    cudaSafe(cudaFree(d_pp_first_inds));
    cudaSafe(cudaFree(d_pp_second_inds));
    cudaSafe(cudaFree(d_pd_first_inds));
    cudaSafe(cudaFree(d_pd_second_inds));
    cudaSafe(cudaFree(d_dd_first_inds));
    cudaSafe(cudaFree(d_dd_second_inds));

    cudaSafe(cudaFree(d_ss_pair_data));
    cudaSafe(cudaFree(d_sp_pair_data));
    cudaSafe(cudaFree(d_sd_pair_data));
    cudaSafe(cudaFree(d_pp_pair_data));
    cudaSafe(cudaFree(d_pd_pair_data));
    cudaSafe(cudaFree(d_dd_pair_data));

    cudaSafe(cudaFree(d_ss_mat_Q_local));
    cudaSafe(cudaFree(d_sp_mat_Q_local));
    cudaSafe(cudaFree(d_sd_mat_Q_local));
    cudaSafe(cudaFree(d_pp_mat_Q_local));
    cudaSafe(cudaFree(d_pd_mat_Q_local));
    cudaSafe(cudaFree(d_dd_mat_Q_local));

    cudaSafe(cudaFree(d_ss_first_inds_local));
    cudaSafe(cudaFree(d_ss_second_inds_local));
    cudaSafe(cudaFree(d_sp_first_inds_local));
    cudaSafe(cudaFree(d_sp_second_inds_local));
    cudaSafe(cudaFree(d_sd_first_inds_local));
    cudaSafe(cudaFree(d_sd_second_inds_local));
    cudaSafe(cudaFree(d_pp_first_inds_local));
    cudaSafe(cudaFree(d_pp_second_inds_local));
    cudaSafe(cudaFree(d_pd_first_inds_local));
    cudaSafe(cudaFree(d_pd_second_inds_local));
    cudaSafe(cudaFree(d_dd_first_inds_local));
    cudaSafe(cudaFree(d_dd_second_inds_local));

    cudaSafe(cudaFree(d_ss_pair_data_local));
    cudaSafe(cudaFree(d_sp_pair_data_local));
    cudaSafe(cudaFree(d_sd_pair_data_local));
    cudaSafe(cudaFree(d_pp_pair_data_local));
    cudaSafe(cudaFree(d_pd_pair_data_local));
    cudaSafe(cudaFree(d_dd_pair_data_local));

    timer.start("Exchange prep.");

    // K preparation

    const auto& Q_K_ss = screening.get_Q_K_ss();
    const auto& Q_K_sp = screening.get_Q_K_sp();
    const auto& Q_K_ps = screening.get_Q_K_ps();
    const auto& Q_K_sd = screening.get_Q_K_sd();
    const auto& Q_K_ds = screening.get_Q_K_ds();
    const auto& Q_K_pp = screening.get_Q_K_pp();
    const auto& Q_K_pd = screening.get_Q_K_pd();
    const auto& Q_K_dp = screening.get_Q_K_dp();
    const auto& Q_K_dd = screening.get_Q_K_dd();

    const auto& D_inds_K_ss = screening.get_D_inds_K_ss();
    const auto& D_inds_K_sp = screening.get_D_inds_K_sp();
    const auto& D_inds_K_ps = screening.get_D_inds_K_ps();
    const auto& D_inds_K_sd = screening.get_D_inds_K_sd();
    const auto& D_inds_K_ds = screening.get_D_inds_K_ds();
    const auto& D_inds_K_pp = screening.get_D_inds_K_pp();
    const auto& D_inds_K_pd = screening.get_D_inds_K_pd();
    const auto& D_inds_K_dp = screening.get_D_inds_K_dp();
    const auto& D_inds_K_dd = screening.get_D_inds_K_dd();

    const auto& pair_displs_K_ss = screening.get_pair_displs_K_ss();
    const auto& pair_displs_K_sp = screening.get_pair_displs_K_sp();
    const auto& pair_displs_K_ps = screening.get_pair_displs_K_ps();
    const auto& pair_displs_K_sd = screening.get_pair_displs_K_sd();
    const auto& pair_displs_K_ds = screening.get_pair_displs_K_ds();
    const auto& pair_displs_K_pp = screening.get_pair_displs_K_pp();
    const auto& pair_displs_K_pd = screening.get_pair_displs_K_pd();
    const auto& pair_displs_K_dp = screening.get_pair_displs_K_dp();
    const auto& pair_displs_K_dd = screening.get_pair_displs_K_dd();

    const auto& pair_counts_K_ss = screening.get_pair_counts_K_ss();
    const auto& pair_counts_K_sp = screening.get_pair_counts_K_sp();
    const auto& pair_counts_K_ps = screening.get_pair_counts_K_ps();
    const auto& pair_counts_K_sd = screening.get_pair_counts_K_sd();
    const auto& pair_counts_K_ds = screening.get_pair_counts_K_ds();
    const auto& pair_counts_K_pp = screening.get_pair_counts_K_pp();
    const auto& pair_counts_K_pd = screening.get_pair_counts_K_pd();
    const auto& pair_counts_K_dp = screening.get_pair_counts_K_dp();
    const auto& pair_counts_K_dd = screening.get_pair_counts_K_dd();

    const auto& pair_data_K_ss = screening.get_pair_data_K_ss();
    const auto& pair_data_K_sp = screening.get_pair_data_K_sp();
    const auto& pair_data_K_ps = screening.get_pair_data_K_ps();
    const auto& pair_data_K_sd = screening.get_pair_data_K_sd();
    const auto& pair_data_K_ds = screening.get_pair_data_K_ds();
    const auto& pair_data_K_pp = screening.get_pair_data_K_pp();
    const auto& pair_data_K_pd = screening.get_pair_data_K_pd();
    const auto& pair_data_K_dp = screening.get_pair_data_K_dp();
    const auto& pair_data_K_dd = screening.get_pair_data_K_dd();

    const auto ss_max_D = screening.get_ss_max_D(); 
    const auto sp_max_D = screening.get_sp_max_D(); 
    const auto sd_max_D = screening.get_sd_max_D(); 
    const auto pp_max_D = screening.get_pp_max_D(); 
    const auto pd_max_D = screening.get_pd_max_D(); 
    const auto dd_max_D = screening.get_dd_max_D(); 

    // Note: assuming symmetric D
    // TODO: double check if this needs to be changed as it requires
    // the density to be either symmetric or antisymmetric
    const auto ps_max_D = screening.get_sp_max_D(); 
    const auto ds_max_D = screening.get_sd_max_D(); 
    const auto dp_max_D = screening.get_pd_max_D(); 

    const auto& pair_inds_i_for_K_ss = screening.get_local_pair_inds_i_for_K_ss(gpu_id);
    const auto& pair_inds_k_for_K_ss = screening.get_local_pair_inds_k_for_K_ss(gpu_id);

    const auto& pair_inds_i_for_K_sp = screening.get_local_pair_inds_i_for_K_sp(gpu_id);
    const auto& pair_inds_k_for_K_sp = screening.get_local_pair_inds_k_for_K_sp(gpu_id);

    const auto& pair_inds_i_for_K_sd = screening.get_local_pair_inds_i_for_K_sd(gpu_id);
    const auto& pair_inds_k_for_K_sd = screening.get_local_pair_inds_k_for_K_sd(gpu_id);

    const auto& pair_inds_i_for_K_pp = screening.get_local_pair_inds_i_for_K_pp(gpu_id);
    const auto& pair_inds_k_for_K_pp = screening.get_local_pair_inds_k_for_K_pp(gpu_id);

    const auto& pair_inds_i_for_K_pd = screening.get_local_pair_inds_i_for_K_pd(gpu_id);
    const auto& pair_inds_k_for_K_pd = screening.get_local_pair_inds_k_for_K_pd(gpu_id);

    const auto& pair_inds_i_for_K_dd = screening.get_local_pair_inds_i_for_K_dd(gpu_id);
    const auto& pair_inds_k_for_K_dd = screening.get_local_pair_inds_k_for_K_dd(gpu_id);

    auto pair_inds_count_for_K_ss = static_cast<uint32_t>(pair_inds_i_for_K_ss.size());
    auto pair_inds_count_for_K_sp = static_cast<uint32_t>(pair_inds_i_for_K_sp.size());
    auto pair_inds_count_for_K_sd = static_cast<uint32_t>(pair_inds_i_for_K_sd.size());
    auto pair_inds_count_for_K_pp = static_cast<uint32_t>(pair_inds_i_for_K_pp.size());
    auto pair_inds_count_for_K_pd = static_cast<uint32_t>(pair_inds_i_for_K_pd.size());
    auto pair_inds_count_for_K_dd = static_cast<uint32_t>(pair_inds_i_for_K_dd.size());

    double*   d_mat_K;
    uint32_t *d_pair_inds_i_for_K_ss, *d_pair_inds_k_for_K_ss;
    uint32_t *d_pair_inds_i_for_K_sp, *d_pair_inds_k_for_K_sp;
    uint32_t *d_pair_inds_i_for_K_sd, *d_pair_inds_k_for_K_sd;
    uint32_t *d_pair_inds_i_for_K_pp, *d_pair_inds_k_for_K_pp;
    uint32_t *d_pair_inds_i_for_K_pd, *d_pair_inds_k_for_K_pd;
    uint32_t *d_pair_inds_i_for_K_dd, *d_pair_inds_k_for_K_dd;

    const auto max_pair_inds_count = std::max({pair_inds_count_for_K_ss, pair_inds_count_for_K_sp, pair_inds_count_for_K_pp,
                                               pair_inds_count_for_K_sd, pair_inds_count_for_K_pd, pair_inds_count_for_K_dd});

    std::vector<double> mat_K(max_pair_inds_count);

    cudaSafe(cudaMalloc(&d_mat_K, max_pair_inds_count * sizeof(double)));

    cudaSafe(cudaMalloc(&d_pair_inds_i_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_inds_k_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_inds_i_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_inds_k_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_inds_i_for_K_sd, pair_inds_count_for_K_sd * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_inds_k_for_K_sd, pair_inds_count_for_K_sd * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_inds_i_for_K_pp, pair_inds_count_for_K_pp * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_inds_k_for_K_pp, pair_inds_count_for_K_pp * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_inds_i_for_K_pd, pair_inds_count_for_K_pd * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_inds_k_for_K_pd, pair_inds_count_for_K_pd * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_inds_i_for_K_dd, pair_inds_count_for_K_dd * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_inds_k_for_K_dd, pair_inds_count_for_K_dd * sizeof(uint32_t)));

    cudaSafe(cudaMemcpy(d_pair_inds_i_for_K_ss, pair_inds_i_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_inds_k_for_K_ss, pair_inds_k_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_inds_i_for_K_sp, pair_inds_i_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_inds_k_for_K_sp, pair_inds_k_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_inds_i_for_K_sd, pair_inds_i_for_K_sd.data(), pair_inds_count_for_K_sd * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_inds_k_for_K_sd, pair_inds_k_for_K_sd.data(), pair_inds_count_for_K_sd * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_inds_i_for_K_pp, pair_inds_i_for_K_pp.data(), pair_inds_count_for_K_pp * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_inds_k_for_K_pp, pair_inds_k_for_K_pp.data(), pair_inds_count_for_K_pp * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_inds_i_for_K_pd, pair_inds_i_for_K_pd.data(), pair_inds_count_for_K_pd * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_inds_k_for_K_pd, pair_inds_k_for_K_pd.data(), pair_inds_count_for_K_pd * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_inds_i_for_K_dd, pair_inds_i_for_K_dd.data(), pair_inds_count_for_K_dd * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_inds_k_for_K_dd, pair_inds_k_for_K_dd.data(), pair_inds_count_for_K_dd * sizeof(uint32_t), cudaMemcpyHostToDevice));

    double* d_mat_D_full_AO;

    double   *d_Q_K_ss, *d_Q_K_sp, *d_Q_K_ps,
             *d_Q_K_sd, *d_Q_K_ds, *d_Q_K_pp,
             *d_Q_K_pd, *d_Q_K_dp, *d_Q_K_dd;

    uint32_t *d_D_inds_K_ss, *d_D_inds_K_sp, *d_D_inds_K_ps,
             *d_D_inds_K_sd, *d_D_inds_K_ds, *d_D_inds_K_pp,
             *d_D_inds_K_pd, *d_D_inds_K_dp, *d_D_inds_K_dd;

    uint32_t *d_pair_displs_K_ss, *d_pair_displs_K_sp, *d_pair_displs_K_ps,
             *d_pair_displs_K_sd, *d_pair_displs_K_ds, *d_pair_displs_K_pp,
             *d_pair_displs_K_pd, *d_pair_displs_K_dp, *d_pair_displs_K_dd;

    uint32_t *d_pair_counts_K_ss, *d_pair_counts_K_sp, *d_pair_counts_K_ps,
             *d_pair_counts_K_sd, *d_pair_counts_K_ds, *d_pair_counts_K_pp,
             *d_pair_counts_K_pd, *d_pair_counts_K_dp, *d_pair_counts_K_dd;

    double   *d_pair_data_K_ss, *d_pair_data_K_sp, *d_pair_data_K_ps,
             *d_pair_data_K_sd, *d_pair_data_K_ds, *d_pair_data_K_pp,
             *d_pair_data_K_pd, *d_pair_data_K_dp, *d_pair_data_K_dd;

    cudaSafe(cudaMalloc(&d_mat_D_full_AO, cart_naos * cart_naos * sizeof(double)));

    cudaSafe(cudaMalloc(&d_Q_K_ss, Q_K_ss.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_sp, Q_K_sp.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_ps, Q_K_ps.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_sd, Q_K_sd.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_ds, Q_K_ds.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_pp, Q_K_pp.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_pd, Q_K_pd.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_dp, Q_K_dp.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Q_K_dd, Q_K_dd.size() * sizeof(double)));

    cudaSafe(cudaMalloc(&d_D_inds_K_ss, D_inds_K_ss.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_sp, D_inds_K_sp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_ps, D_inds_K_ps.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_sd, D_inds_K_sd.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_ds, D_inds_K_ds.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_pp, D_inds_K_pp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_pd, D_inds_K_pd.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_dp, D_inds_K_dp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_D_inds_K_dd, D_inds_K_dd.size() * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_displs_K_ss, pair_displs_K_ss.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_sp, pair_displs_K_sp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_ps, pair_displs_K_ps.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_sd, pair_displs_K_sd.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_ds, pair_displs_K_ds.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_pp, pair_displs_K_pp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_pd, pair_displs_K_pd.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_dp, pair_displs_K_dp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_displs_K_dd, pair_displs_K_dd.size() * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_counts_K_ss, pair_counts_K_ss.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_sp, pair_counts_K_sp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_ps, pair_counts_K_ps.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_sd, pair_counts_K_sd.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_ds, pair_counts_K_ds.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_pp, pair_counts_K_pp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_pd, pair_counts_K_pd.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_dp, pair_counts_K_dp.size() * sizeof(uint32_t)));
    cudaSafe(cudaMalloc(&d_pair_counts_K_dd, pair_counts_K_dd.size() * sizeof(uint32_t)));

    cudaSafe(cudaMalloc(&d_pair_data_K_ss, pair_data_K_ss.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_sp, pair_data_K_sp.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_ps, pair_data_K_ps.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_sd, pair_data_K_sd.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_ds, pair_data_K_ds.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_pp, pair_data_K_pp.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_pd, pair_data_K_pd.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_dp, pair_data_K_dp.size() * sizeof(double)));
    cudaSafe(cudaMalloc(&d_pair_data_K_dd, pair_data_K_dd.size() * sizeof(double)));

    cudaSafe(cudaMemcpy(d_mat_D_full_AO, cart_dens_ptr, cart_naos * cart_naos * sizeof(double), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_Q_K_ss, Q_K_ss.data(), Q_K_ss.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_sp, Q_K_sp.data(), Q_K_sp.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_ps, Q_K_ps.data(), Q_K_ps.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_sd, Q_K_sd.data(), Q_K_sd.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_ds, Q_K_ds.data(), Q_K_ds.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_pp, Q_K_pp.data(), Q_K_pp.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_pd, Q_K_pd.data(), Q_K_pd.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_dp, Q_K_dp.data(), Q_K_dp.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_Q_K_dd, Q_K_dd.data(), Q_K_dd.size() * sizeof(double), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_D_inds_K_ss, D_inds_K_ss.data(), D_inds_K_ss.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_sp, D_inds_K_sp.data(), D_inds_K_sp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_ps, D_inds_K_ps.data(), D_inds_K_ps.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_sd, D_inds_K_sd.data(), D_inds_K_sd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_ds, D_inds_K_ds.data(), D_inds_K_ds.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_pp, D_inds_K_pp.data(), D_inds_K_pp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_pd, D_inds_K_pd.data(), D_inds_K_pd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_dp, D_inds_K_dp.data(), D_inds_K_dp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_D_inds_K_dd, D_inds_K_dd.data(), D_inds_K_dd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_displs_K_ss, pair_displs_K_ss.data(), pair_displs_K_ss.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_sp, pair_displs_K_sp.data(), pair_displs_K_sp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_ps, pair_displs_K_ps.data(), pair_displs_K_ps.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_sd, pair_displs_K_sd.data(), pair_displs_K_sd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_ds, pair_displs_K_ds.data(), pair_displs_K_ds.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_pp, pair_displs_K_pp.data(), pair_displs_K_pp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_pd, pair_displs_K_pd.data(), pair_displs_K_pd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_dp, pair_displs_K_dp.data(), pair_displs_K_dp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_displs_K_dd, pair_displs_K_dd.data(), pair_displs_K_dd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_counts_K_ss, pair_counts_K_ss.data(), pair_counts_K_ss.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_sp, pair_counts_K_sp.data(), pair_counts_K_sp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_ps, pair_counts_K_ps.data(), pair_counts_K_ps.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_sd, pair_counts_K_sd.data(), pair_counts_K_sd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_ds, pair_counts_K_ds.data(), pair_counts_K_ds.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_pp, pair_counts_K_pp.data(), pair_counts_K_pp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_pd, pair_counts_K_pd.data(), pair_counts_K_pd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_dp, pair_counts_K_dp.data(), pair_counts_K_dp.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_counts_K_dd, pair_counts_K_dd.data(), pair_counts_K_dd.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    cudaSafe(cudaMemcpy(d_pair_data_K_ss, pair_data_K_ss.data(), pair_data_K_ss.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_sp, pair_data_K_sp.data(), pair_data_K_sp.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_ps, pair_data_K_ps.data(), pair_data_K_ps.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_sd, pair_data_K_sd.data(), pair_data_K_sd.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_ds, pair_data_K_ds.data(), pair_data_K_ds.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_pp, pair_data_K_pp.data(), pair_data_K_pp.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_pd, pair_data_K_pd.data(), pair_data_K_pd.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_dp, pair_data_K_dp.data(), pair_data_K_dp.size() * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_pair_data_K_dd, pair_data_K_dd.data(), pair_data_K_dd.size() * sizeof(double), cudaMemcpyHostToDevice));

    timer.stop("Exchange prep.");

    cudaSafe(cudaDeviceSynchronize());

    CTimer exchange_timer;

    exchange_timer.start();

    timer.start("K computation");

    // compute K

    if (std::fabs(frac_exact_exchange) > 1.0e-13)
    {

    // K: S-S block

    if (pair_inds_count_for_K_ss > 0)
    {
        timer.start("  K block SS");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_ss + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_K,static_cast<uint32_t>(pair_inds_count_for_K_ss));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_ss, 1);

        // K: (SS|SS)
        //     *  *

        gpu::computeExchangeFockSSSS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           ss_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_D_inds_K_ss,
                           d_pair_displs_K_ss,
                           d_pair_counts_K_ss,
                           d_pair_data_K_ss,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SS|SP)
        //     *  *

        gpu::computeExchangeFockSSSP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           sp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_sp,
                           d_D_inds_K_ss,
                           d_D_inds_K_sp,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_sp,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_sp,
                           d_pair_data_K_ss,
                           d_pair_data_K_sp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|SS)
        //     *  *

        gpu::computeExchangeFockSPSS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           ps_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_ss,
                           d_D_inds_K_sp,
                           d_D_inds_K_ss,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_ss,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_ss,
                           d_pair_data_K_sp,
                           d_pair_data_K_ss,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|SP)
        //     *  *

        gpu::computeExchangeFockSPSP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           pp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_D_inds_K_sp,
                           d_pair_displs_K_sp,
                           d_pair_counts_K_sp,
                           d_pair_data_K_sp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SS|SD)
        //     *  *

        gpu::computeExchangeFockSSSD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_sd,
                           d_D_inds_K_ss,
                           d_D_inds_K_sd,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_sd,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_sd,
                           d_pair_data_K_ss,
                           d_pair_data_K_sd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SD|SS)
        //     *  *

        gpu::computeExchangeFockSDSS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ds_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_ss,
                           d_D_inds_K_sd,
                           d_D_inds_K_ss,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_ss,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_ss,
                           d_pair_data_K_sd,
                           d_pair_data_K_ss,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|SD)
        //     *  *

        gpu::computeExchangeFockSPSD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_sd,
                           d_D_inds_K_sp,
                           d_D_inds_K_sd,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_sd,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_sd,
                           d_pair_data_K_sp,
                           d_pair_data_K_sd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SD|SP)
        //     *  *

        gpu::computeExchangeFockSDSP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_sp,
                           d_D_inds_K_sd,
                           d_D_inds_K_sp,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_sp,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_sp,
                           d_pair_data_K_sd,
                           d_pair_data_K_sp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SD|SD)
        //     *  *

        gpu::computeExchangeFockSDSD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_ss,
                           d_pair_inds_k_for_K_ss,
                           static_cast<uint32_t>(pair_inds_count_for_K_ss),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_D_inds_K_sd,
                           d_pair_displs_K_sd,
                           d_pair_counts_K_sd,
                           d_pair_data_K_sd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0) * frac_exact_exchange;

            const auto symm_pref = (fstr::upcase(flag_K) == "SYMM") ? 1.0 : -1.0;

            if (i != k) mat_Fock_omp[gpu_id].row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0) * frac_exact_exchange * symm_pref;
        }

        timer.stop("  K block SS");
    }

    // K: S-P block

    if (pair_inds_count_for_K_sp > 0)
    {
        timer.start("  K block SP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sp + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_K,static_cast<uint32_t>(pair_inds_count_for_K_sp));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_sp, 1);

        // K: (SS|PS)
        //     *  *

        gpu::computeExchangeFockSSPS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           ss_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_ps,
                           d_D_inds_K_ss,
                           d_D_inds_K_ps,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_ps,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_ps,
                           d_pair_data_K_ss,
                           d_pair_data_K_ps,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SS|PP)
        //     *  *

        gpu::computeExchangeFockSSPP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           sp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_pp,
                           d_D_inds_K_ss,
                           d_D_inds_K_pp,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_pp,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_pp,
                           d_pair_data_K_ss,
                           d_pair_data_K_pp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|PS)
        //     *  *

        gpu::computeExchangeFockSPPS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           ps_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_ps,
                           d_D_inds_K_sp,
                           d_D_inds_K_ps,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_ps,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_ps,
                           d_pair_data_K_sp,
                           d_pair_data_K_ps,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|PP)
        //     *  *

        gpu::computeExchangeFockSPPP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           pp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_pp,
                           d_D_inds_K_sp,
                           d_D_inds_K_pp,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_pp,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_pp,
                           d_pair_data_K_sp,
                           d_pair_data_K_pp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SS|PD)
        //     *  *

        gpu::computeExchangeFockSSPD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_pd,
                           d_D_inds_K_ss,
                           d_D_inds_K_pd,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_pd,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_pd,
                           d_pair_data_K_ss,
                           d_pair_data_K_pd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SD|PS)
        //     *  *

        gpu::computeExchangeFockSDPS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ds_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_ps,
                           d_D_inds_K_sd,
                           d_D_inds_K_ps,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_ps,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_ps,
                           d_pair_data_K_sd,
                           d_pair_data_K_ps,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|PD)
        //     *  *

        gpu::computeExchangeFockSPPD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_pd,
                           d_D_inds_K_sp,
                           d_D_inds_K_pd,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_pd,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_pd,
                           d_pair_data_K_sp,
                           d_pair_data_K_pd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SD|PP)
        //     *  *

        gpu::computeExchangeFockSDPP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_pp,
                           d_D_inds_K_sd,
                           d_D_inds_K_pp,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_pp,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_pp,
                           d_pair_data_K_sd,
                           d_pair_data_K_pp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SD|PD)
        //     *  *

        gpu::computeExchangeFockSDPD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sp,
                           d_pair_inds_k_for_K_sp,
                           static_cast<uint32_t>(pair_inds_count_for_K_sp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_pd,
                           d_D_inds_K_sd,
                           d_D_inds_K_pd,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_pd,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_pd,
                           d_pair_data_K_sd,
                           d_pair_data_K_pd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sp; ik++)
        {
            const auto i = pair_inds_i_for_K_sp[ik];
            const auto k = pair_inds_k_for_K_sp[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            // Cartesian to spherical
            for (const auto& k_cgto_sph_ind_coef : cart_sph_p[k_cgto])
            {
                auto k_cgto_sph = k_cgto_sph_ind_coef.first;
                auto k_coef_sph = k_cgto_sph_ind_coef.second;

                mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto_sph] += mat_K[ik] * k_coef_sph * (-1.0) * frac_exact_exchange;

                const auto symm_pref = (fstr::upcase(flag_K) == "SYMM") ? 1.0 : -1.0;

                mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto] += mat_K[ik] * k_coef_sph * (-1.0) * frac_exact_exchange * symm_pref;
            }
        }

        timer.stop("  K block SP");
    }

    // K: P-P block

    if (pair_inds_count_for_K_pp > 0)
    {
        timer.start("  K block PP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pp + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_K,static_cast<uint32_t>(pair_inds_count_for_K_pp));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_pp, 1);

        // K: (PS|PS)
        //     *  *

        timer.start("    K block PSPS");

        gpu::computeExchangeFockPSPS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           ss_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ps,
                           d_D_inds_K_ps,
                           d_pair_displs_K_ps,
                           d_pair_counts_K_ps,
                           d_pair_data_K_ps,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PSPS");

        // K: (PS|PP)
        //     *  *

        timer.start("    K block PSPP");

        gpu::computeExchangeFockPSPP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           sp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ps,
                           d_Q_K_pp,
                           d_D_inds_K_ps,
                           d_D_inds_K_pp,
                           d_pair_displs_K_ps,
                           d_pair_displs_K_pp,
                           d_pair_counts_K_ps,
                           d_pair_counts_K_pp,
                           d_pair_data_K_ps,
                           d_pair_data_K_pp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PSPP");

        // K: (PP|PS)
        //     *  *

        timer.start("    K block PPPS");

        gpu::computeExchangeFockPPPS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           ps_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pp,
                           d_Q_K_ps,
                           d_D_inds_K_pp,
                           d_D_inds_K_ps,
                           d_pair_displs_K_pp,
                           d_pair_displs_K_ps,
                           d_pair_counts_K_pp,
                           d_pair_counts_K_ps,
                           d_pair_data_K_pp,
                           d_pair_data_K_ps,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PPPS");

        // K: (PP|PP)
        //     *  *

        timer.start("    K block PPPP");

        gpu::computeExchangeFockPPPP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           pp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pp,
                           d_D_inds_K_pp,
                           d_pair_displs_K_pp,
                           d_pair_counts_K_pp,
                           d_pair_data_K_pp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PPPP");

        // K: (PS|PD)
        //     *  *

        timer.start("    K block PSPD");

        gpu::computeExchangeFockPSPD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ps,
                           d_Q_K_pd,
                           d_D_inds_K_ps,
                           d_D_inds_K_pd,
                           d_pair_displs_K_ps,
                           d_pair_displs_K_pd,
                           d_pair_counts_K_ps,
                           d_pair_counts_K_pd,
                           d_pair_data_K_ps,
                           d_pair_data_K_pd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PSPD");

        // K: (PD|PS)
        //     *  *

        timer.start("    K block PDPS");

        gpu::computeExchangeFockPDPS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ds_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_ps,
                           d_D_inds_K_pd,
                           d_D_inds_K_ps,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_ps,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_ps,
                           d_pair_data_K_pd,
                           d_pair_data_K_ps,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PDPS");

        // K: (PP|PD)
        //     *  *

        timer.start("    K block PPPD");

        gpu::computeExchangeFockPPPD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pp,
                           d_Q_K_pd,
                           d_D_inds_K_pp,
                           d_D_inds_K_pd,
                           d_pair_displs_K_pp,
                           d_pair_displs_K_pd,
                           d_pair_counts_K_pp,
                           d_pair_counts_K_pd,
                           d_pair_data_K_pp,
                           d_pair_data_K_pd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PPPD");

        // K: (PD|PP)
        //     *  *

        timer.start("    K block PDPP");

        gpu::computeExchangeFockPDPP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_pp,
                           d_D_inds_K_pd,
                           d_D_inds_K_pp,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_pp,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_pp,
                           d_pair_data_K_pd,
                           d_pair_data_K_pp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PDPP");

        // K: (PD|PD)
        //     *  *

        timer.start("    K block PDPD");

        gpu::computeExchangeFockPDPD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pp,
                           d_pair_inds_k_for_K_pp,
                           static_cast<uint32_t>(pair_inds_count_for_K_pp),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_D_inds_K_pd,
                           d_pair_displs_K_pd,
                           d_pair_counts_K_pd,
                           d_pair_data_K_pd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PDPD");

        cudaSafe(cudaMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pp; ik++)
        {
            const auto i = pair_inds_i_for_K_pp[ik];
            const auto k = pair_inds_k_for_K_pp[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& k_cgto_sph_ind_coef : cart_sph_p[k_cgto])
                {
                    auto k_cgto_sph = k_cgto_sph_ind_coef.first;
                    auto k_coef_sph = k_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * k_coef_sph;

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[k_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange;

                    const auto symm_pref = (fstr::upcase(flag_K) == "SYMM") ? 1.0 : -1.0;

                    if (i != k) mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange * symm_pref;
                }
            }
        }

        timer.stop("  K block PP");
    }

    // K: S-D block

    if (pair_inds_count_for_K_sd > 0)
    {
        timer.start("  K block SD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sd + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_K,static_cast<uint32_t>(pair_inds_count_for_K_sd));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_sd, 1);

        // K: (SS|DS)
        //     *  *

        gpu::computeExchangeFockSSDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ss_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_ds,
                           d_D_inds_K_ss,
                           d_D_inds_K_ds,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_ds,
                           d_pair_data_K_ss,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SS|DP)
        //     *  *

        gpu::computeExchangeFockSSDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_dp,
                           d_D_inds_K_ss,
                           d_D_inds_K_dp,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_dp,
                           d_pair_data_K_ss,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|DS)
        //     *  *

        gpu::computeExchangeFockSPDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ps_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_ds,
                           d_D_inds_K_sp,
                           d_D_inds_K_ds,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_ds,
                           d_pair_data_K_sp,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (SP|DP)
        //     *  *

        timer.start("    K block SPDP");

        gpu::computeExchangeFockSPDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_dp,
                           d_D_inds_K_sp,
                           d_D_inds_K_dp,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_dp,
                           d_pair_data_K_sp,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block SPDP");

        // K: (SS|DD)
        //     *  *

        timer.start("    K block SSDD");

        gpu::computeExchangeFockSSDD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ss,
                           d_Q_K_dd,
                           d_D_inds_K_ss,
                           d_D_inds_K_dd,
                           d_pair_displs_K_ss,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_ss,
                           d_pair_counts_K_dd,
                           d_pair_data_K_ss,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block SSDD");

        // K: (SD|DS)
        //     *  *

        timer.start("    K block SDDS");

        gpu::computeExchangeFockSDDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ds_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_ds,
                           d_D_inds_K_sd,
                           d_D_inds_K_ds,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_ds,
                           d_pair_data_K_sd,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block SDDS");

        // K: (SP|DD)
        //     *  *

        timer.start("    K block SPDD");

        gpu::computeExchangeFockSPDD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sp,
                           d_Q_K_dd,
                           d_D_inds_K_sp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_sp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_sp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_sp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block SPDD");

        // K: (SD|DP)
        //     *  *

        timer.start("    K block SDDP");

        gpu::computeExchangeFockSDDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_dp,
                           d_D_inds_K_sd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_sd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block SDDP");

        // K: (SD|DD)
        //     *  *

        timer.start("    K block SDDD");

        gpu::computeExchangeFockSDDD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_sd,
                           d_pair_inds_k_for_K_sd,
                           static_cast<uint32_t>(pair_inds_count_for_K_sd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_sd,
                           d_Q_K_dd,
                           d_D_inds_K_sd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_sd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_sd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_sd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block SDDD");

        cudaSafe(cudaMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sd * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sd; ik++)
        {
            const auto i = pair_inds_i_for_K_sd[ik];
            const auto k = pair_inds_k_for_K_sd[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = d_prim_aoinds[(k / 6) + d_prim_count * (k % 6)];

            // Cartesian to spherical
            for (const auto& k_cgto_sph_ind_coef : cart_sph_d[k_cgto])
            {
                auto k_cgto_sph = k_cgto_sph_ind_coef.first;
                auto k_coef_sph = k_cgto_sph_ind_coef.second;

                mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto_sph] += mat_K[ik] * k_coef_sph * (-1.0) * frac_exact_exchange;

                const auto symm_pref = (fstr::upcase(flag_K) == "SYMM") ? 1.0 : -1.0;

                mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto] += mat_K[ik] * k_coef_sph * (-1.0) * frac_exact_exchange * symm_pref;
            }
        }

        timer.stop("  K block SD");
    }

    // K: P-D block

    if (pair_inds_count_for_K_pd > 0)
    {
        timer.start("  K block PD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pd + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_K,static_cast<uint32_t>(pair_inds_count_for_K_pd));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_pd, 1);

        // K: (PS|DS)
        //     *  *

        timer.start("    K block PSDS");

        gpu::computeExchangeFockPSDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ss_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ps,
                           d_Q_K_ds,
                           d_D_inds_K_ps,
                           d_D_inds_K_ds,
                           d_pair_displs_K_ps,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_ps,
                           d_pair_counts_K_ds,
                           d_pair_data_K_ps,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PSDS");

        // K: (PS|DP)
        //     *  *

        timer.start("    K block PSDP");

        gpu::computeExchangeFockPSDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ps,
                           d_Q_K_dp,
                           d_D_inds_K_ps,
                           d_D_inds_K_dp,
                           d_pair_displs_K_ps,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_ps,
                           d_pair_counts_K_dp,
                           d_pair_data_K_ps,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PSDP");

        // K: (PP|DS)
        //     *  *

        timer.start("    K block PPDS");

        gpu::computeExchangeFockPPDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ps_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pp,
                           d_Q_K_ds,
                           d_D_inds_K_pp,
                           d_D_inds_K_ds,
                           d_pair_displs_K_pp,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_pp,
                           d_pair_counts_K_ds,
                           d_pair_data_K_pp,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PPDS");

        // K: (PS|DD)
        //     *  *

        timer.start("    K block PSDD");

        gpu::computeExchangeFockPSDD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ps,
                           d_Q_K_dd,
                           d_D_inds_K_ps,
                           d_D_inds_K_dd,
                           d_pair_displs_K_ps,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_ps,
                           d_pair_counts_K_dd,
                           d_pair_data_K_ps,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PSDD");

        // K: (PD|DS)
        //     *  *

        timer.start("    K block PDDS");

        gpu::computeExchangeFockPDDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ds_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_ds,
                           d_D_inds_K_pd,
                           d_D_inds_K_ds,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_ds,
                           d_pair_data_K_pd,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PDDS");

        // K: (PP|DP)
        //     *  *

        timer.start("    K block PPDP");

        gpu::computeExchangeFockPPDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pp,
                           d_Q_K_dp,
                           d_D_inds_K_pp,
                           d_D_inds_K_dp,
                           d_pair_displs_K_pp,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_pp,
                           d_pair_counts_K_dp,
                           d_pair_data_K_pp,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PPDP");

        // K: (PP|DD)
        //     *  *

        timer.start("    K block PPDD");

        gpu::computeExchangeFockPPDD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pp,
                           d_Q_K_dd,
                           d_D_inds_K_pp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PPDD");

        // K: (PD|DP)
        //     *  *

        timer.start("    K block PDDP");

        gpu::computeExchangeFockPDDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dp,
                           d_D_inds_K_pd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_pd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PDDP");

        // K: (PD|DD)
        //     *  *

        timer.start("    K block PDDD");

        gpu::computeExchangeFockPDDD0<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockPDDD1<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockPDDD2<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockPDDD3<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockPDDD4<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockPDDD5<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockPDDD6<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_pd,
                           d_pair_inds_k_for_K_pd,
                           static_cast<uint32_t>(pair_inds_count_for_K_pd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_pd,
                           d_Q_K_dd,
                           d_D_inds_K_pd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_pd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_pd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_pd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block PDDD");

        cudaSafe(cudaMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pd * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pd; ik++)
        {
            const auto i = pair_inds_i_for_K_pd[ik];
            const auto k = pair_inds_k_for_K_pd[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = d_prim_aoinds[(k / 6) + d_prim_count * (k % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& k_cgto_sph_ind_coef : cart_sph_d[k_cgto])
                {
                    auto k_cgto_sph = k_cgto_sph_ind_coef.first;
                    auto k_coef_sph = k_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * k_coef_sph;

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[k_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange;

                    const auto symm_pref = (fstr::upcase(flag_K) == "SYMM") ? 1.0 : -1.0;

                    mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange * symm_pref;
                }
            }
        }

        timer.stop("  K block PD");
    }

    // K: D-D block

    if (pair_inds_count_for_K_dd > 0)
    {
        timer.start("  K block DD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_dd + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks,threads_per_block>>>(d_mat_K,static_cast<uint32_t>(pair_inds_count_for_K_dd));

        cudaSafe(cudaDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_dd, 1);

        // K: (DS|DS)
        //     *  *

        gpu::computeExchangeFockDSDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ss_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ds,
                           d_D_inds_K_ds,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_ds,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (DS|DP)
        //     *  *

        gpu::computeExchangeFockDSDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ds,
                           d_Q_K_dp,
                           d_D_inds_K_ds,
                           d_D_inds_K_dp,
                           d_pair_displs_K_ds,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_ds,
                           d_pair_counts_K_dp,
                           d_pair_data_K_ds,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (DP|DS)
        //     *  *

        gpu::computeExchangeFockDPDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ps_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_ds,
                           d_D_inds_K_dp,
                           d_D_inds_K_ds,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_ds,
                           d_pair_data_K_dp,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        // K: (DS|DD)
        //     *  *

        timer.start("    K block DSDD");

        gpu::computeExchangeFockDSDD<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           sd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_ds,
                           d_Q_K_dd,
                           d_D_inds_K_ds,
                           d_D_inds_K_dd,
                           d_pair_displs_K_ds,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_ds,
                           d_pair_counts_K_dd,
                           d_pair_data_K_ds,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block DSDD");

        // K: (DD|DS)
        //     *  *

        timer.start("    K block DDDS");

        gpu::computeExchangeFockDDDS<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_s_prim_info,
                           d_s_prim_aoinds,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           ds_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_ds,
                           d_D_inds_K_dd,
                           d_D_inds_K_ds,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_ds,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_ds,
                           d_pair_data_K_dd,
                           d_pair_data_K_ds,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block DDDS");

        // K: (DP|DP)
        //     *  *

        timer.start("    K block DPDP");

        gpu::computeExchangeFockDPDP<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block DPDP");

        // K: (DP|DD)
        //     *  *

        timer.start("    K block DPDD");

        gpu::computeExchangeFockDPDD0<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDPDD1<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDPDD2<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDPDD3<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDPDD4<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDPDD5<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDPDD6<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           pd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dp,
                           d_Q_K_dd,
                           d_D_inds_K_dp,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dp,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block DPDD");

        // K: (DD|DP)
        //     *  *

        timer.start("    K block DDDP");

        gpu::computeExchangeFockDDDP0<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDP1<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDP2<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDP3<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDP4<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDP5<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDP6<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_p_prim_info,
                           d_p_prim_aoinds,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dp_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_Q_K_dp,
                           d_D_inds_K_dd,
                           d_D_inds_K_dp,
                           d_pair_displs_K_dd,
                           d_pair_displs_K_dp,
                           d_pair_counts_K_dd,
                           d_pair_counts_K_dp,
                           d_pair_data_K_dd,
                           d_pair_data_K_dp,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block DDDP");

        // K: (DD|DD)
        //     *  *

        timer.start("    K block DDDD");

        gpu::computeExchangeFockDDDD0<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD1<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD2<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD3<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD4<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD5<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD6<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD7<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD8<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD9<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD10<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD11<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD12<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD13<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        gpu::computeExchangeFockDDDD14<<<num_blocks,threads_per_block>>>(
                           d_mat_K,
                           d_pair_inds_i_for_K_dd,
                           d_pair_inds_k_for_K_dd,
                           static_cast<uint32_t>(pair_inds_count_for_K_dd),
                           d_d_prim_info,
                           d_d_prim_aoinds,
                           static_cast<uint32_t>(d_prim_count),
                           dd_max_D,
                           d_mat_D_full_AO,
                           static_cast<uint32_t>(cart_naos),
                           d_Q_K_dd,
                           d_D_inds_K_dd,
                           d_pair_displs_K_dd,
                           d_pair_counts_K_dd,
                           d_pair_data_K_dd,
                           d_boys_func_table,
                           d_boys_func_ft,
                           omega,
                           eri_threshold);

        cudaSafe(cudaDeviceSynchronize());

        timer.stop("    K block DDDD");

        cudaSafe(cudaMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_dd * sizeof(double), cudaMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_dd; ik++)
        {
            const auto i = pair_inds_i_for_K_dd[ik];
            const auto k = pair_inds_k_for_K_dd[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
            const auto k_cgto = d_prim_aoinds[(k / 6) + d_prim_count * (k % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& k_cgto_sph_ind_coef : cart_sph_d[k_cgto])
                {
                    auto k_cgto_sph = k_cgto_sph_ind_coef.first;
                    auto k_coef_sph = k_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * k_coef_sph;

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[k_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange;

                    const auto symm_pref = (fstr::upcase(flag_K) == "SYMM") ? 1.0 : -1.0;

                    if (i != k) mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange * symm_pref;
                }
            }
        }

        timer.stop("  K block DD");
    }

    }  // end of compute K

    timer.stop("K computation");

    cudaSafe(cudaDeviceSynchronize());

    exchange_timer.stop();

    auto exchange_elapsed_time = exchange_timer.getElapsedTime();

    screening.setExchangeTime(gpu_id, exchange_elapsed_time);

    cudaSafe(cudaFree(d_boys_func_table));
    cudaSafe(cudaFree(d_boys_func_ft));

    cudaSafe(cudaFree(d_s_prim_info));
    cudaSafe(cudaFree(d_s_prim_aoinds));

    cudaSafe(cudaFree(d_p_prim_info));
    cudaSafe(cudaFree(d_p_prim_aoinds));

    cudaSafe(cudaFree(d_d_prim_info));
    cudaSafe(cudaFree(d_d_prim_aoinds));

    cudaSafe(cudaFree(d_mat_K));

    cudaSafe(cudaFree(d_pair_inds_i_for_K_ss));
    cudaSafe(cudaFree(d_pair_inds_k_for_K_ss));
    cudaSafe(cudaFree(d_pair_inds_i_for_K_sp));
    cudaSafe(cudaFree(d_pair_inds_k_for_K_sp));
    cudaSafe(cudaFree(d_pair_inds_i_for_K_sd));
    cudaSafe(cudaFree(d_pair_inds_k_for_K_sd));
    cudaSafe(cudaFree(d_pair_inds_i_for_K_pp));
    cudaSafe(cudaFree(d_pair_inds_k_for_K_pp));
    cudaSafe(cudaFree(d_pair_inds_i_for_K_pd));
    cudaSafe(cudaFree(d_pair_inds_k_for_K_pd));
    cudaSafe(cudaFree(d_pair_inds_i_for_K_dd));
    cudaSafe(cudaFree(d_pair_inds_k_for_K_dd));

    cudaSafe(cudaFree(d_mat_D_full_AO));

    cudaSafe(cudaFree(d_Q_K_ss));
    cudaSafe(cudaFree(d_Q_K_sp));
    cudaSafe(cudaFree(d_Q_K_ps));
    cudaSafe(cudaFree(d_Q_K_sd));
    cudaSafe(cudaFree(d_Q_K_ds));
    cudaSafe(cudaFree(d_Q_K_pp));
    cudaSafe(cudaFree(d_Q_K_pd));
    cudaSafe(cudaFree(d_Q_K_dp));
    cudaSafe(cudaFree(d_Q_K_dd));

    cudaSafe(cudaFree(d_D_inds_K_ss));
    cudaSafe(cudaFree(d_D_inds_K_sp));
    cudaSafe(cudaFree(d_D_inds_K_ps));
    cudaSafe(cudaFree(d_D_inds_K_sd));
    cudaSafe(cudaFree(d_D_inds_K_ds));
    cudaSafe(cudaFree(d_D_inds_K_pp));
    cudaSafe(cudaFree(d_D_inds_K_pd));
    cudaSafe(cudaFree(d_D_inds_K_dp));
    cudaSafe(cudaFree(d_D_inds_K_dd));

    cudaSafe(cudaFree(d_pair_displs_K_ss));
    cudaSafe(cudaFree(d_pair_displs_K_sp));
    cudaSafe(cudaFree(d_pair_displs_K_ps));
    cudaSafe(cudaFree(d_pair_displs_K_sd));
    cudaSafe(cudaFree(d_pair_displs_K_ds));
    cudaSafe(cudaFree(d_pair_displs_K_pp));
    cudaSafe(cudaFree(d_pair_displs_K_pd));
    cudaSafe(cudaFree(d_pair_displs_K_dp));
    cudaSafe(cudaFree(d_pair_displs_K_dd));

    cudaSafe(cudaFree(d_pair_counts_K_ss));
    cudaSafe(cudaFree(d_pair_counts_K_sp));
    cudaSafe(cudaFree(d_pair_counts_K_ps));
    cudaSafe(cudaFree(d_pair_counts_K_sd));
    cudaSafe(cudaFree(d_pair_counts_K_ds));
    cudaSafe(cudaFree(d_pair_counts_K_pp));
    cudaSafe(cudaFree(d_pair_counts_K_pd));
    cudaSafe(cudaFree(d_pair_counts_K_dp));
    cudaSafe(cudaFree(d_pair_counts_K_dd));

    cudaSafe(cudaFree(d_pair_data_K_ss));
    cudaSafe(cudaFree(d_pair_data_K_sp));
    cudaSafe(cudaFree(d_pair_data_K_ps));
    cudaSafe(cudaFree(d_pair_data_K_sd));
    cudaSafe(cudaFree(d_pair_data_K_ds));
    cudaSafe(cudaFree(d_pair_data_K_pp));
    cudaSafe(cudaFree(d_pair_data_K_pd));
    cudaSafe(cudaFree(d_pair_data_K_dp));
    cudaSafe(cudaFree(d_pair_data_K_dd));

    timer.stop("Total timing");

    // std::cout << "\nTiming of ERIs on GPU " << gpu_rank << "\n";
    // std::cout << "-----------------------\n";
    // std::cout << timer.getSummary() << std::endl;

    }}

    CDenseMatrix mat_Fock_sum (naos, naos);

    mat_Fock_sum.zero();

    auto p_mat_F = mat_Fock_sum.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_fock = mat_Fock_omp[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_F[ind] += p_mat_fock[ind];
        }
    }

    return mat_Fock_sum;
}

auto
computeDotProduct(const double* A, const double* B, const int64_t size) -> double
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeDotProduct: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_A, *d_B;

    cudaSafe(cudaMalloc(&d_A, n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_B, n * sizeof(double)));

    cudaSafe(cudaMemcpy(d_A, A, n * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_B, B, n * sizeof(double), cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double dot_product;
    cublasSafe(cublasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_B));

    return dot_product;
}

auto
computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& data_pointers, const int64_t size) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeWeightedSum: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_X, *d_Y;

    cudaSafe(cudaMalloc(&d_X, n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Y, n * sizeof(double)));

    cudaSafe(cudaMemcpy(d_Y, weighted_data, n * sizeof(double), cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        cudaSafe(cudaMemcpy(d_X, data_pointers[i], n * sizeof(double), cudaMemcpyHostToDevice));

        cublasSafe(cublasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
    }

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaMemcpy(weighted_data, d_Y, n * sizeof(double), cudaMemcpyDeviceToHost));

    cudaSafe(cudaFree(d_X));
    cudaSafe(cudaFree(d_Y));
}

auto
computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                   const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeErrorVector: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_A, *d_B, *d_C;

    cudaSafe(cudaMalloc(&d_A, nao * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_B, nao * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_C, nao * nao * sizeof(double)));

    cudaSafe(cudaMemcpy(d_A, F, nao * nao * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_B, D, nao * nao * sizeof(double), cudaMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

    // Note: we compute C^T = B^T * A^T since cublas is column-major

    // D^T F^T (=> FD)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));

    // S^T(FD)^T (=> FDS)
    cudaSafe(cudaMemcpy(d_A, S, nao * nao * sizeof(double), cudaMemcpyHostToDevice));

    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));

    cudaSafe(cudaDeviceSynchronize());

    // FDS - (FDS)^T
    beta = -1.0;

    cublasSafe(cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_T, nao, nao, &alpha, d_B, nao, &beta, d_B, nao, d_C, nao));

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    cudaSafe(cudaMemcpy(d_X, X, nmo * nao * sizeof(double), cudaMemcpyHostToDevice));

    auto op_X  = (trans_X == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? CUBLAS_OP_T : CUBLAS_OP_N;

    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    beta = 0.0;

    // TODO: double check basis set with linear dependency

    // let E == (FDS - SDF)
    // X^T E^T (=> EX)
    cublasSafe(cublasDgemm(handle, op_X, CUBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_C, nao, &beta, d_Y, nmo));

    cudaSafe(cudaDeviceSynchronize());

    // (EX)^T X (=> X^T(E)X)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_C, nmo));

    cudaSafe(cudaMemcpy(errvec, d_C, nmo * nmo * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_B));
    cudaSafe(cudaFree(d_C));
}

auto
transformMatrix(double* transformed_F, const double* X, const double* F,
                const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::transformMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_F, *d_X, *d_Y;

    cudaSafe(cudaMalloc(&d_F, nao * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_X, nmo * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Y, nmo * nao * sizeof(double)));

    cudaSafe(cudaMemcpy(d_F, F, nao * nao * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_X, X, nmo * nao * sizeof(double), cudaMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

    auto op_X = (trans_X == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    auto op_XT = (trans_X == std::string("N")) ? CUBLAS_OP_T : CUBLAS_OP_N;
    auto lda_XT = (trans_X == std::string("N")) ? nao : nmo;

    // Note: we compute C^T = B^T * A^T since cublas is column-major

    // X^T F^T (=> FX)
    cublasSafe(cublasDgemm(handle, op_X, CUBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_F, nao, &beta, d_Y, nmo));

    cudaSafe(cudaDeviceSynchronize());

    // (FX)^T X (=> X^T(F)X)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_F, nmo));

    cudaSafe(cudaMemcpy(transformed_F, d_F, nmo * nmo * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_F));
    cudaSafe(cudaFree(d_X));
    cudaSafe(cudaFree(d_Y));
}

auto
computeMatrixMultiplication(double* C, const double* A, const double* B, const std::string& trans_A, const std::string& trans_B,
                            const int64_t m_inp, const int64_t k_inp, const int64_t n_inp) -> void
{
    // TODO allow computeMatrixMultiplication on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeMatrixMultiplication: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    // TODO: matmul on multiple GPUs

    cudaSafe(cudaSetDevice(0));

    auto m = static_cast<int32_t>(m_inp);
    auto k = static_cast<int32_t>(k_inp);
    auto n = static_cast<int32_t>(n_inp);

    double *d_A, *d_B, *d_C;

    cudaSafe(cudaMalloc(&d_A, m * k * sizeof(double)));
    cudaSafe(cudaMalloc(&d_B, k * n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_C, m * n * sizeof(double)));

    cudaSafe(cudaMemcpy(d_A, A, m * k * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_B, B, k * n * sizeof(double), cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto op_A = (trans_A == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;

    // Note: we compute C^T = B^T * A^T since cublas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

    cublasSafe(cublasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));

    cudaSafe(cudaMemcpy(C, d_C, m * n * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_B));
    cudaSafe(cudaFree(d_C));
}

auto
diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void
{
    cudaSafe(cudaSetDevice(0));

    auto n = static_cast<int32_t>(nrows_A);
    int32_t lwork;

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    cudaSafe(cudaMalloc(&d_A, n * n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_D, n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_info, sizeof(int32_t)));

    cudaSafe(cudaMemcpy(d_A, A, n * n * sizeof(double), cudaMemcpyHostToDevice));

    cusolverDnHandle_t handle;
    cusolverSafe(cusolverDnCreate(&handle));

    cusolverSafe(cusolverDnDsyevd_bufferSize(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    cudaSafe(cudaMalloc(&d_work, lwork));

    cusolverSafe(cusolverDnDsyevd(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    cusolverSafe(cusolverDnDestroy(handle));

    cudaSafe(cudaMemcpy(A, d_A, n * n * sizeof(double), cudaMemcpyDeviceToHost));
    cudaSafe(cudaMemcpy(D, d_D, n * sizeof(double), cudaMemcpyDeviceToHost));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_D));
    cudaSafe(cudaFree(d_info));
    cudaSafe(cudaFree(d_work));
}

}  // namespace gpu
