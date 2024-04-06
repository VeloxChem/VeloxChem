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
#include <hipblas.h>
//#include <hipsolver.h>
#include <magma_v2.h>

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

    // std::stringstream ss;
    // ss << "rank " << rank << ",  thread " << thread_id << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    // TODO: double check by using different number of gpus per MPI
    hipSafe(hipSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    hipSafe(hipMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    hipSafe(hipMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));

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

    hipSafe(hipMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    hipSafe(hipMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    hipSafe(hipMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));

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

    hipSafe(hipMalloc(&d_mat_Q, max_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipDeviceSynchronize());

    Q_omp[gpu_id].zero();

    auto& mat_Q_omp = Q_omp[gpu_id];

    // compute Q

    // Q: SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeQMatrixSS, num_blocks, threads_per_block, 0, 0,
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        hipSafe(hipMemcpy(h_mat_Q.data(), d_mat_Q, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

        hipLaunchKernelGGL(gpu::computeQMatrixSP, num_blocks, threads_per_block, 0, 0,
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

        hipSafe(hipMemcpy(h_mat_Q.data(), d_mat_Q, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

        hipLaunchKernelGGL(gpu::computeQMatrixSD, num_blocks, threads_per_block, 0, 0,
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

        hipSafe(hipMemcpy(h_mat_Q.data(), d_mat_Q, sd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

        hipLaunchKernelGGL(gpu::computeQMatrixPP, num_blocks, threads_per_block, 0, 0,
                           d_mat_Q,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        hipSafe(hipMemcpy(h_mat_Q.data(), d_mat_Q, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

        hipLaunchKernelGGL(gpu::computeQMatrixPD, num_blocks, threads_per_block, 0, 0,
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

        hipSafe(hipMemcpy(h_mat_Q.data(), d_mat_Q, pd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

        hipLaunchKernelGGL(gpu::computeQMatrixDD, num_blocks, threads_per_block, 0, 0,
                           d_mat_Q,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        hipSafe(hipMemcpy(h_mat_Q.data(), d_mat_Q, dd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            mat_Q_omp.row(s_prim_count + p_prim_count * 3 + i)[s_prim_count + p_prim_count * 3 + j] = h_mat_Q[ij];

            if (i != j) mat_Q_omp.row(s_prim_count + p_prim_count * 3 + j)[s_prim_count + p_prim_count * 3 + i] = h_mat_Q[ij];
        }
    }

    hipSafe(hipDeviceSynchronize());

    hipSafe(hipFree(d_boys_func_table));
    hipSafe(hipFree(d_boys_func_ft));

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_p_prim_info));
    hipSafe(hipFree(d_d_prim_info));

    hipSafe(hipFree(d_mat_Q));

    hipSafe(hipFree(d_ss_first_inds_local));
    hipSafe(hipFree(d_ss_second_inds_local));
    hipSafe(hipFree(d_sp_first_inds_local));
    hipSafe(hipFree(d_sp_second_inds_local));
    hipSafe(hipFree(d_sd_first_inds_local));
    hipSafe(hipFree(d_sd_second_inds_local));
    hipSafe(hipFree(d_pp_first_inds_local));
    hipSafe(hipFree(d_pp_second_inds_local));
    hipSafe(hipFree(d_pd_first_inds_local));
    hipSafe(hipFree(d_pd_second_inds_local));
    hipSafe(hipFree(d_dd_first_inds_local));
    hipSafe(hipFree(d_dd_second_inds_local));

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
computeOneElectronIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>
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
    std::vector<CDenseMatrix> V_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        S_matrices[gpu_id] = CDenseMatrix(naos, naos);
        T_matrices[gpu_id] = CDenseMatrix(naos, naos);
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

    // std::stringstream ss;
    // ss << "rank " << rank << ",  thread " << thread_id << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    hipSafe(hipSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    hipSafe(hipMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    hipSafe(hipMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));

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

    hipSafe(hipMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    hipSafe(hipMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    hipSafe(hipMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));

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
    std::vector<double> mat_V(max_prim_pair_count_local);

    // S on device

    double *d_mat_S, *d_mat_T, *d_mat_V;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    hipSafe(hipMalloc(&d_mat_S, max_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_T, max_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_V, max_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    S_matrices[gpu_id].zero();
    T_matrices[gpu_id].zero();
    V_matrices[gpu_id].zero();

    auto& mat_overlap = S_matrices[gpu_id];
    auto& mat_kinetic_energy = T_matrices[gpu_id];
    auto& mat_nuclear_potential = V_matrices[gpu_id];

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

    double *d_points_info;

    hipSafe(hipMalloc(&d_points_info, natoms * 4 * sizeof(double)));

    hipSafe(hipMemcpy(d_points_info, points_info.data(), natoms * 4 * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipDeviceSynchronize());

    // compute S

    // S: SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergySS, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialSS, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

            mat_nuclear_potential.row(i_cgto)[j_cgto] += mat_V[ij];

            if (i != j) mat_nuclear_potential.row(j_cgto)[i_cgto] += mat_V[ij];
        }
    }

    // S: SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergySP, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialSP, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                mat_nuclear_potential.row(i_cgto)[j_cgto_sph] += mat_V[ij] * j_coef_sph;
                mat_nuclear_potential.row(j_cgto_sph)[i_cgto] += mat_V[ij] * j_coef_sph;
            }
        }
    }

    // S: SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergySD, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_d_prim_info,
                                                     static_cast<uint32_t>(d_prim_count),
                                                     d_sd_first_inds_local,
                                                     d_sd_second_inds_local,
                                                     static_cast<uint32_t>(sd_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, sd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, sd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialSD, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_d_prim_info,
                                                     static_cast<uint32_t>(d_prim_count),
                                                     d_sd_first_inds_local,
                                                     d_sd_second_inds_local,
                                                     static_cast<uint32_t>(sd_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, sd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                mat_nuclear_potential.row(i_cgto)[j_cgto_sph] += mat_V[ij] * j_coef_sph;
                mat_nuclear_potential.row(j_cgto_sph)[i_cgto] += mat_V[ij] * j_coef_sph;
            }
        }
    }

    // S: PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergyPP, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialPP, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                    mat_nuclear_potential.row(i_cgto_sph)[j_cgto_sph] += mat_V[ij] * coef_sph;

                    if (i != j) mat_nuclear_potential.row(j_cgto_sph)[i_cgto_sph] += mat_V[ij] * coef_sph;
                }
            }
        }
    }

    // S: PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergyPD, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_d_prim_info,
                                                     static_cast<uint32_t>(d_prim_count),
                                                     d_pd_first_inds_local,
                                                     d_pd_second_inds_local,
                                                     static_cast<uint32_t>(pd_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, pd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, pd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialPD, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_d_prim_info,
                                                     static_cast<uint32_t>(d_prim_count),
                                                     d_pd_first_inds_local,
                                                     d_pd_second_inds_local,
                                                     static_cast<uint32_t>(pd_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, pd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                    mat_nuclear_potential.row(i_cgto_sph)[j_cgto_sph] += mat_V[ij] * coef_sph;
                    mat_nuclear_potential.row(j_cgto_sph)[i_cgto_sph] += mat_V[ij] * coef_sph;
                }
            }
        }
    }

    // S: DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergyDD, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_d_prim_info,
                                                     static_cast<uint32_t>(d_prim_count),
                                                     d_dd_first_inds_local,
                                                     d_dd_second_inds_local,
                                                     static_cast<uint32_t>(dd_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, dd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, dd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialDD, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_d_prim_info,
                                                     static_cast<uint32_t>(d_prim_count),
                                                     d_dd_first_inds_local,
                                                     d_dd_second_inds_local,
                                                     static_cast<uint32_t>(dd_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, dd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                    mat_nuclear_potential.row(i_cgto_sph)[j_cgto_sph] += mat_V[ij] * coef_sph;

                    if (i != j) mat_nuclear_potential.row(j_cgto_sph)[i_cgto_sph] += mat_V[ij] * coef_sph;
                }
            }
        }
    }

    hipSafe(hipDeviceSynchronize());

    hipSafe(hipFree(d_boys_func_table));
    hipSafe(hipFree(d_boys_func_ft));

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_p_prim_info));
    hipSafe(hipFree(d_d_prim_info));

    hipSafe(hipFree(d_mat_S));
    hipSafe(hipFree(d_mat_T));
    hipSafe(hipFree(d_mat_V));

    hipSafe(hipFree(d_ss_first_inds_local));
    hipSafe(hipFree(d_ss_second_inds_local));
    hipSafe(hipFree(d_sp_first_inds_local));
    hipSafe(hipFree(d_sp_second_inds_local));
    hipSafe(hipFree(d_sd_first_inds_local));
    hipSafe(hipFree(d_sd_second_inds_local));
    hipSafe(hipFree(d_pp_first_inds_local));
    hipSafe(hipFree(d_pp_second_inds_local));
    hipSafe(hipFree(d_pd_first_inds_local));
    hipSafe(hipFree(d_pd_second_inds_local));
    hipSafe(hipFree(d_dd_first_inds_local));
    hipSafe(hipFree(d_dd_second_inds_local));

    hipSafe(hipFree(d_points_info));

    }}

    std::vector<CDenseMatrix> STV_matrices(3);

    STV_matrices[0] = CDenseMatrix(naos, naos);
    STV_matrices[1] = CDenseMatrix(naos, naos);
    STV_matrices[2] = CDenseMatrix(naos, naos);

    STV_matrices[0].zero();
    STV_matrices[1].zero();
    STV_matrices[2].zero();

    auto p_mat_S = STV_matrices[0].values();
    auto p_mat_T = STV_matrices[1].values();
    auto p_mat_V = STV_matrices[2].values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_overlap = S_matrices[gpu_id].values();
        auto p_mat_kinetic_energy = T_matrices[gpu_id].values();
        auto p_mat_nuclear_potential = V_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_S[ind] += p_mat_overlap[ind];
            p_mat_T[ind] += p_mat_kinetic_energy[ind];
            p_mat_V[ind] += p_mat_nuclear_potential[ind];
        }
    }

    return STV_matrices;
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
computeFockOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix, const double frac_exact_exchange, const double eri_threshold, const double prelink_threshold, CScreeningData& screening) -> CDenseMatrix
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    auto gpu_rank = rank * num_gpus_per_node;

    // std::stringstream ss;
    // ss << "\n(Prep.) rank " << rank << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    hipSafe(hipSetDevice(gpu_rank % total_num_gpus_per_compute_node));

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

    timer.start("Prep. Q_prime");

    // preLinK
    // J. Chem. Phys. 138, 134114 (2013)

    // TODO distribute computation of Q_prime

    double *d_matrix_A, *d_matrix_B, *d_matrix_C;

    auto mat_full = screening.get_mat_Q_full(s_prim_count, p_prim_count, d_prim_count);

    hipSafe(hipMalloc(&d_matrix_A, mat_full.getNumberOfElements() * sizeof(double)));
    hipSafe(hipMalloc(&d_matrix_B, mat_full.getNumberOfElements() * sizeof(double)));
    hipSafe(hipMalloc(&d_matrix_C, mat_full.getNumberOfElements() * sizeof(double)));

    hipSafe(hipMemcpy(d_matrix_A, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), hipMemcpyHostToDevice));

    mat_full = screening.get_mat_D_abs_full(s_prim_count, p_prim_count, d_prim_count, s_prim_aoinds, p_prim_aoinds, d_prim_aoinds, cart_naos, cart_dens_ptr);

    hipSafe(hipMemcpy(d_matrix_B, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), hipMemcpyHostToDevice));

    const auto all_prim_count = mat_full.getNumberOfRows();

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto n = static_cast<int32_t>(all_prim_count);

    // compute A^T * (B^T * A^T) since hipblas is column-major
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, n, n, n, &alpha, d_matrix_B, n, d_matrix_A, n, &beta, d_matrix_C, n));

    hipSafe(hipDeviceSynchronize());

    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, n, n, n, &alpha, d_matrix_A, n, d_matrix_C, n, &beta, d_matrix_B, n));

    hipSafe(hipMemcpy(mat_full.values(), d_matrix_B, mat_full.getNumberOfElements() * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_matrix_A));
    hipSafe(hipFree(d_matrix_B));
    hipSafe(hipFree(d_matrix_C));

    timer.stop("Prep. Q_prime");

    timer.start("Prep. preLinK");

    screening.form_pair_inds_for_K(s_prim_count, p_prim_count, d_prim_count, mat_full, prelink_threshold);

    mat_full = CDenseMatrix();

    timer.stop("Prep. preLinK");

    std::string errnaos("gpu::computeFockGPU: Inconsistent number of AOs");
    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    std::vector<CDenseMatrix> mat_Fock_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Fock_omp[gpu_id] = CDenseMatrix(naos, naos);
    }

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

    // std::stringstream ss;
    // ss << "rank " << rank << ",  thread " << thread_id << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    hipSafe(hipSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Boys func. prep.");

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    hipSafe(hipMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    hipSafe(hipMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));

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

    hipSafe(hipMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_s_prim_aoinds, s_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;
    uint32_t* d_p_prim_aoinds;

    hipSafe(hipMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_p_prim_aoinds, p_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    double*   d_d_prim_info;
    uint32_t* d_d_prim_aoinds;

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    hipSafe(hipMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_d_prim_aoinds, d_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_d_prim_aoinds, d_prim_aoinds.data(), d_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

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

    hipSafe(hipMalloc(&d_mat_D, max_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_J, max_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_mat_Q, ss_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_sp_mat_Q, sp_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_sd_mat_Q, sd_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_pp_mat_Q, pp_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_pd_mat_Q, pd_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_dd_mat_Q, dd_prim_pair_count * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds, ss_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds, ss_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds, sp_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds, sp_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sd_first_inds, sd_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sd_second_inds, sd_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds, pp_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds, pp_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pd_first_inds, pd_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pd_second_inds, pd_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_dd_first_inds, dd_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_dd_second_inds, dd_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_ss_pair_data, ss_pair_data.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_sp_pair_data, sp_pair_data.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_sd_pair_data, sd_pair_data.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pp_pair_data, pp_pair_data.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pd_pair_data, pd_pair_data.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_dd_pair_data, dd_pair_data.size() * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_mat_Q_local, ss_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_sp_mat_Q_local, sp_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_sd_mat_Q_local, sd_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_pp_mat_Q_local, pp_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_pd_mat_Q_local, pd_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_dd_mat_Q_local, dd_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_ss_pair_data_local, ss_pair_data_local.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_sp_pair_data_local, sp_pair_data_local.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_sd_pair_data_local, sd_pair_data_local.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pp_pair_data_local, pp_pair_data_local.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pd_pair_data_local, pd_pair_data_local.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_dd_pair_data_local, dd_pair_data_local.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_ss_mat_Q, ss_mat_Q.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_mat_Q, sp_mat_Q.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_mat_Q, sd_mat_Q.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_mat_Q, pp_mat_Q.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_mat_Q, pd_mat_Q.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_mat_Q, dd_mat_Q.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_first_inds, ss_first_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds, ss_second_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds, sp_first_inds.data(), sp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds, sp_second_inds.data(), sp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sd_first_inds, sd_first_inds.data(), sd_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_second_inds, sd_second_inds.data(), sd_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds, pp_first_inds.data(), pp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds, pp_second_inds.data(), pp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pd_first_inds, pd_first_inds.data(), pd_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_second_inds, pd_second_inds.data(), pd_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_dd_first_inds, dd_first_inds.data(), dd_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_second_inds, dd_second_inds.data(), dd_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_pair_data, ss_pair_data.data(), ss_pair_data.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_pair_data, sp_pair_data.data(), sp_pair_data.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_pair_data, sd_pair_data.data(), sd_pair_data.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_pair_data, pp_pair_data.data(), pp_pair_data.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_pair_data, pd_pair_data.data(), pd_pair_data.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_pair_data, dd_pair_data.data(), dd_pair_data.size() * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_mat_Q_local, ss_mat_Q_local.data(), ss_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_mat_Q_local, sp_mat_Q_local.data(), sp_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_mat_Q_local, sd_mat_Q_local.data(), sd_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_mat_Q_local, pp_mat_Q_local.data(), pp_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_mat_Q_local, pd_mat_Q_local.data(), pd_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_mat_Q_local, dd_mat_Q_local.data(), dd_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_pair_data_local, ss_pair_data_local.data(), ss_pair_data_local.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_pair_data_local, sp_pair_data_local.data(), sp_pair_data_local.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sd_pair_data_local, sd_pair_data_local.data(), sd_pair_data_local.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_pair_data_local, pp_pair_data_local.data(), pp_pair_data_local.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pd_pair_data_local, pd_pair_data_local.data(), pd_pair_data_local.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_dd_pair_data_local, dd_pair_data_local.data(), dd_pair_data_local.size() * sizeof(double), hipMemcpyHostToDevice));

    timer.stop("Coulomb prep.");

    mat_Fock_omp[gpu_id].zero();

    hipSafe(hipDeviceSynchronize());

    timer.start("J computation");

    // compute J

    // J: S-S block

    if (ss_prim_pair_count_local > 0)
    {
        timer.start("  J block SS");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(ss_prim_pair_count_local));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SS|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            timer.start("    J block SSSS");

            hipLaunchKernelGGL(gpu::computeCoulombFockSSSS, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block SSSS");
        }

        // J: (SS|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSSP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SS|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSSD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SS|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSPP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SS|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSPD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SS|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSDD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
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

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(sp_prim_pair_count_local));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPSS, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPSP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPSD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPPP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPPD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPDD0, num_blocks, threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSPDD1, num_blocks, threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSPDD2, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto_sph] += mat_J[ij] * j_coef_sph * 2.0;
                mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto] += mat_J[ij] * j_coef_sph * 2.0;
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

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(pp_prim_pair_count_local));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPSS, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPSP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPSD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPPP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPPD0, num_blocks, threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPPD1, num_blocks, threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPPD2, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block PPDD");

            hipSafe(hipMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 pp_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 pp_num_blocks ((pp_prim_pair_count_local + pp_threads_per_block.x - 1) / pp_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD0, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD1, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD2, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD3, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD4, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD5, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD6, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD7, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPPDD8, pp_num_blocks, pp_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block PPDD");
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[j_cgto_sph] += mat_J[ij] * coef_sph * 2.0;

                    if (i != j) mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto_sph] += mat_J[ij] * coef_sph * 2.0;
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

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(sd_prim_pair_count_local));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSDSS, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSDSP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSDSD, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSDPP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSDPD0, num_blocks, threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDPD1, num_blocks, threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDPD2, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (SD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block SDDD");

            hipSafe(hipMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 sd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 sd_num_blocks ((sd_prim_pair_count_local + sd_threads_per_block.x - 1) / sd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockSDDD0, sd_num_blocks, sd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDDD1, sd_num_blocks, sd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDDD2, sd_num_blocks, sd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDDD3, sd_num_blocks, sd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDDD4, sd_num_blocks, sd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockSDDD5, sd_num_blocks, sd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block SDDD");
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, sd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto_sph] += mat_J[ij] * j_coef_sph * 2.0;
                mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto] += mat_J[ij] * j_coef_sph * 2.0;
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

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(pd_prim_pair_count_local));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPDSS, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPDSP, num_blocks, threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());
        }

        // J: (PD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            timer.start("    J block PDSD");

            hipSafe(hipMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 pd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 pd_num_blocks ((pd_prim_pair_count_local + pd_threads_per_block.x - 1) / pd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockPDSD0, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDSD1, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block PDSD");
        }

        // J: (PD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            timer.start("    J block PDPP");

            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 pd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 pd_num_blocks ((pd_prim_pair_count_local + pd_threads_per_block.x - 1) / pd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPP0, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPP1, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block PDPP");
        }

        // J: (PD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            timer.start("    J block PDPD");

            hipSafe(hipMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 pd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 pd_num_blocks ((pd_prim_pair_count_local + pd_threads_per_block.x - 1) / pd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD0, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD1, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD2, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD3, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD4, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD5, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDPD6, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block PDPD");
        }

        // J: (PD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block PDDD");

            hipSafe(hipMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 pd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 pd_num_blocks ((pd_prim_pair_count_local + pd_threads_per_block.x - 1) / pd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD0, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD1, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD2, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD3, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD4, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD5, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD6, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD7, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD8, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD9, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD10, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD11, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD12, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD13, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD14, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD15, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD16, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD17, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD18, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockPDDD19, pd_num_blocks, pd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block PDDD");
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, pd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[j_cgto_sph] += mat_J[ij] * coef_sph * 2.0;
                    mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto_sph] += mat_J[ij] * coef_sph * 2.0;
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

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(dd_prim_pair_count_local));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (DD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            timer.start("    J block DDSS");

            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSS, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block DDSS");
        }

        // J: (DD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            timer.start("    J block DDSP");

            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSP0, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSP1, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSP2, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block DDSP");
        }

        // J: (DD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            timer.start("    J block DDSD");

            hipSafe(hipMemcpy(d_mat_D, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSD0, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSD1, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSD2, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSD3, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSD4, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDSD5, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block DDSD");
        }

        // J: (DD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            timer.start("    J block DDPP");

            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPP0, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPP1, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPP2, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPP3, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPP4, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPP5, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block DDPP");
        }

        // J: (DD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            timer.start("    J block DDPD");

            hipSafe(hipMemcpy(d_mat_D, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD0, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD1, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD2, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD3, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD4, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD5, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD6, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD7, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD8, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD9, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD10, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD11, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD12, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD13, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD14, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD15, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD16, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD17, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDPD18, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block DDPD");
        }

        // J: (DD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            timer.start("    J block DDDD");

            hipSafe(hipMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD0, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD1, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD2, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD3, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD4, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD5, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD6, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD7, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD8, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD9, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD10, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD11, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD12, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD13, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD14, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD15, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD16, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD17, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD18, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD19, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD20, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD21, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD22, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD23, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD24, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD25, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD26, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD27, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD28, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD29, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD30, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD31, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD32, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD33, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD34, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD35, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD36, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD37, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD38, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD39, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD40, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD41, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipLaunchKernelGGL(gpu::computeCoulombFockDDDD42, dd_num_blocks, dd_threads_per_block, 0, 0,
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

            hipSafe(hipDeviceSynchronize());

            timer.stop("    J block DDDD");
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, dd_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

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

                    mat_Fock_omp[gpu_id].row(i_cgto_sph)[j_cgto_sph] += mat_J[ij] * coef_sph * 2.0;

                    if (i != j) mat_Fock_omp[gpu_id].row(j_cgto_sph)[i_cgto_sph] += mat_J[ij] * coef_sph * 2.0;
                }
            }
        }

        timer.stop("  J block DD");
    }

    timer.stop("J computation");

    hipSafe(hipDeviceSynchronize());

    hipSafe(hipFree(d_mat_D));
    hipSafe(hipFree(d_mat_J));

    hipSafe(hipFree(d_ss_mat_Q));
    hipSafe(hipFree(d_sp_mat_Q));
    hipSafe(hipFree(d_sd_mat_Q));
    hipSafe(hipFree(d_pp_mat_Q));
    hipSafe(hipFree(d_pd_mat_Q));
    hipSafe(hipFree(d_dd_mat_Q));

    hipSafe(hipFree(d_ss_first_inds));
    hipSafe(hipFree(d_ss_second_inds));
    hipSafe(hipFree(d_sp_first_inds));
    hipSafe(hipFree(d_sp_second_inds));
    hipSafe(hipFree(d_sd_first_inds));
    hipSafe(hipFree(d_sd_second_inds));
    hipSafe(hipFree(d_pp_first_inds));
    hipSafe(hipFree(d_pp_second_inds));
    hipSafe(hipFree(d_pd_first_inds));
    hipSafe(hipFree(d_pd_second_inds));
    hipSafe(hipFree(d_dd_first_inds));
    hipSafe(hipFree(d_dd_second_inds));

    hipSafe(hipFree(d_ss_pair_data));
    hipSafe(hipFree(d_sp_pair_data));
    hipSafe(hipFree(d_sd_pair_data));
    hipSafe(hipFree(d_pp_pair_data));
    hipSafe(hipFree(d_pd_pair_data));
    hipSafe(hipFree(d_dd_pair_data));

    hipSafe(hipFree(d_ss_mat_Q_local));
    hipSafe(hipFree(d_sp_mat_Q_local));
    hipSafe(hipFree(d_sd_mat_Q_local));
    hipSafe(hipFree(d_pp_mat_Q_local));
    hipSafe(hipFree(d_pd_mat_Q_local));
    hipSafe(hipFree(d_dd_mat_Q_local));

    hipSafe(hipFree(d_ss_first_inds_local));
    hipSafe(hipFree(d_ss_second_inds_local));
    hipSafe(hipFree(d_sp_first_inds_local));
    hipSafe(hipFree(d_sp_second_inds_local));
    hipSafe(hipFree(d_sd_first_inds_local));
    hipSafe(hipFree(d_sd_second_inds_local));
    hipSafe(hipFree(d_pp_first_inds_local));
    hipSafe(hipFree(d_pp_second_inds_local));
    hipSafe(hipFree(d_pd_first_inds_local));
    hipSafe(hipFree(d_pd_second_inds_local));
    hipSafe(hipFree(d_dd_first_inds_local));
    hipSafe(hipFree(d_dd_second_inds_local));

    hipSafe(hipFree(d_ss_pair_data_local));
    hipSafe(hipFree(d_sp_pair_data_local));
    hipSafe(hipFree(d_sd_pair_data_local));
    hipSafe(hipFree(d_pp_pair_data_local));
    hipSafe(hipFree(d_pd_pair_data_local));
    hipSafe(hipFree(d_dd_pair_data_local));

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

    hipSafe(hipMalloc(&d_mat_K, max_pair_inds_count * sizeof(double)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_sd, pair_inds_count_for_K_sd * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_sd, pair_inds_count_for_K_sd * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_pp, pair_inds_count_for_K_pp * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_pp, pair_inds_count_for_K_pp * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_pd, pair_inds_count_for_K_pd * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_pd, pair_inds_count_for_K_pd * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_dd, pair_inds_count_for_K_dd * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_dd, pair_inds_count_for_K_dd * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_ss, pair_inds_i_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_ss, pair_inds_k_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_sp, pair_inds_i_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_sp, pair_inds_k_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_sd, pair_inds_i_for_K_sd.data(), pair_inds_count_for_K_sd * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_sd, pair_inds_k_for_K_sd.data(), pair_inds_count_for_K_sd * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_pp, pair_inds_i_for_K_pp.data(), pair_inds_count_for_K_pp * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_pp, pair_inds_k_for_K_pp.data(), pair_inds_count_for_K_pp * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_pd, pair_inds_i_for_K_pd.data(), pair_inds_count_for_K_pd * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_pd, pair_inds_k_for_K_pd.data(), pair_inds_count_for_K_pd * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_dd, pair_inds_i_for_K_dd.data(), pair_inds_count_for_K_dd * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_dd, pair_inds_k_for_K_dd.data(), pair_inds_count_for_K_dd * sizeof(uint32_t), hipMemcpyHostToDevice));

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

    hipSafe(hipMalloc(&d_mat_D_full_AO, cart_naos * cart_naos * sizeof(double)));

    hipSafe(hipMalloc(&d_Q_K_ss, Q_K_ss.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_sp, Q_K_sp.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_ps, Q_K_ps.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_sd, Q_K_sd.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_ds, Q_K_ds.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_pp, Q_K_pp.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_pd, Q_K_pd.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_dp, Q_K_dp.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_Q_K_dd, Q_K_dd.size() * sizeof(double)));

    hipSafe(hipMalloc(&d_D_inds_K_ss, D_inds_K_ss.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_sp, D_inds_K_sp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_ps, D_inds_K_ps.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_sd, D_inds_K_sd.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_ds, D_inds_K_ds.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_pp, D_inds_K_pp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_pd, D_inds_K_pd.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_dp, D_inds_K_dp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_D_inds_K_dd, D_inds_K_dd.size() * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_displs_K_ss, pair_displs_K_ss.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_sp, pair_displs_K_sp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_ps, pair_displs_K_ps.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_sd, pair_displs_K_sd.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_ds, pair_displs_K_ds.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_pp, pair_displs_K_pp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_pd, pair_displs_K_pd.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_dp, pair_displs_K_dp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_displs_K_dd, pair_displs_K_dd.size() * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_counts_K_ss, pair_counts_K_ss.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_sp, pair_counts_K_sp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_ps, pair_counts_K_ps.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_sd, pair_counts_K_sd.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_ds, pair_counts_K_ds.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_pp, pair_counts_K_pp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_pd, pair_counts_K_pd.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_dp, pair_counts_K_dp.size() * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_counts_K_dd, pair_counts_K_dd.size() * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_data_K_ss, pair_data_K_ss.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_sp, pair_data_K_sp.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_ps, pair_data_K_ps.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_sd, pair_data_K_sd.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_ds, pair_data_K_ds.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_pp, pair_data_K_pp.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_pd, pair_data_K_pd.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_dp, pair_data_K_dp.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_pair_data_K_dd, pair_data_K_dd.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_mat_D_full_AO, cart_dens_ptr, cart_naos * cart_naos * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_Q_K_ss, Q_K_ss.data(), Q_K_ss.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_sp, Q_K_sp.data(), Q_K_sp.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_ps, Q_K_ps.data(), Q_K_ps.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_sd, Q_K_sd.data(), Q_K_sd.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_ds, Q_K_ds.data(), Q_K_ds.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_pp, Q_K_pp.data(), Q_K_pp.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_pd, Q_K_pd.data(), Q_K_pd.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_dp, Q_K_dp.data(), Q_K_dp.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_Q_K_dd, Q_K_dd.data(), Q_K_dd.size() * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_D_inds_K_ss, D_inds_K_ss.data(), D_inds_K_ss.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_sp, D_inds_K_sp.data(), D_inds_K_sp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_ps, D_inds_K_ps.data(), D_inds_K_ps.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_sd, D_inds_K_sd.data(), D_inds_K_sd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_ds, D_inds_K_ds.data(), D_inds_K_ds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_pp, D_inds_K_pp.data(), D_inds_K_pp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_pd, D_inds_K_pd.data(), D_inds_K_pd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_dp, D_inds_K_dp.data(), D_inds_K_dp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_D_inds_K_dd, D_inds_K_dd.data(), D_inds_K_dd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_displs_K_ss, pair_displs_K_ss.data(), pair_displs_K_ss.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_sp, pair_displs_K_sp.data(), pair_displs_K_sp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_ps, pair_displs_K_ps.data(), pair_displs_K_ps.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_sd, pair_displs_K_sd.data(), pair_displs_K_sd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_ds, pair_displs_K_ds.data(), pair_displs_K_ds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_pp, pair_displs_K_pp.data(), pair_displs_K_pp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_pd, pair_displs_K_pd.data(), pair_displs_K_pd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_dp, pair_displs_K_dp.data(), pair_displs_K_dp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_displs_K_dd, pair_displs_K_dd.data(), pair_displs_K_dd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_counts_K_ss, pair_counts_K_ss.data(), pair_counts_K_ss.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_sp, pair_counts_K_sp.data(), pair_counts_K_sp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_ps, pair_counts_K_ps.data(), pair_counts_K_ps.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_sd, pair_counts_K_sd.data(), pair_counts_K_sd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_ds, pair_counts_K_ds.data(), pair_counts_K_ds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_pp, pair_counts_K_pp.data(), pair_counts_K_pp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_pd, pair_counts_K_pd.data(), pair_counts_K_pd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_dp, pair_counts_K_dp.data(), pair_counts_K_dp.size() * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_counts_K_dd, pair_counts_K_dd.data(), pair_counts_K_dd.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_data_K_ss, pair_data_K_ss.data(), pair_data_K_ss.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_sp, pair_data_K_sp.data(), pair_data_K_sp.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_ps, pair_data_K_ps.data(), pair_data_K_ps.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_sd, pair_data_K_sd.data(), pair_data_K_sd.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_ds, pair_data_K_ds.data(), pair_data_K_ds.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_pp, pair_data_K_pp.data(), pair_data_K_pp.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_pd, pair_data_K_pd.data(), pair_data_K_pd.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_dp, pair_data_K_dp.data(), pair_data_K_dp.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_data_K_dd, pair_data_K_dd.data(), pair_data_K_dd.size() * sizeof(double), hipMemcpyHostToDevice));

    timer.stop("Exchange prep.");

    hipSafe(hipDeviceSynchronize());

    timer.start("K computation");

    // compute K

    // K: S-S block

    if ((pair_inds_count_for_K_ss > 0) && (std::fabs(frac_exact_exchange) > 1.0e-13))
    {
        timer.start("  K block SS");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_ss + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_ss));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_ss, 1);

        // K: (SS|SS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSSS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|SP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSSP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|SS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPSS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|SP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPSP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|SD)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSSD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SD|SS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSDSS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|SD)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPSD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SD|SP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSDSP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SD|SD)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSDSD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0) * frac_exact_exchange;

            if (i != k) mat_Fock_omp[gpu_id].row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0) * frac_exact_exchange;
        }

        timer.stop("  K block SS");
    }

    // K: S-P block

    if ((pair_inds_count_for_K_sp > 0) && (std::fabs(frac_exact_exchange) > 1.0e-13))
    {
        timer.start("  K block SP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sp + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_sp));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_sp, 1);

        // K: (SS|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSPS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSPP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPPS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPPP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|PD)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSPD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SD|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSDPS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|PD)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPPD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SD|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSDPP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SD|PD)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSDPD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDPD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp * sizeof(double), hipMemcpyDeviceToHost));

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
                mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto] += mat_K[ik] * k_coef_sph * (-1.0) * frac_exact_exchange;
            }
        }

        timer.stop("  K block SP");
    }

    // K: P-P block

    if ((pair_inds_count_for_K_pp > 0) && (std::fabs(frac_exact_exchange) > 1.0e-13))
    {
        timer.start("  K block PP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pp + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_pp));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_pp, 1);

        // K: (PS|PS)
        //     *  *

        timer.start("    K block PSPS");

        hipLaunchKernelGGL(gpu::computeExchangeFockPSPS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PSPS");

        // K: (PS|PP)
        //     *  *

        timer.start("    K block PSPP");

        hipLaunchKernelGGL(gpu::computeExchangeFockPSPP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PSPP");

        // K: (PP|PS)
        //     *  *

        timer.start("    K block PPPS");

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PPPS");

        // K: (PP|PP)
        //     *  *

        timer.start("    K block PPPP");

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PPPP");

        // K: (PS|PD)
        //     *  *

        timer.start("    K block PSPD");

        hipLaunchKernelGGL(gpu::computeExchangeFockPSPD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PSPD");

        // K: (PD|PS)
        //     *  *

        timer.start("    K block PDPS");

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PDPS");

        // K: (PP|PD)
        //     *  *

        timer.start("    K block PPPD");

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PPPD");

        // K: (PD|PP)
        //     *  *

        timer.start("    K block PDPP");

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PDPP");

        // K: (PD|PD)
        //     *  *

        timer.start("    K block PDPD");

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDPD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PDPD");

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp * sizeof(double), hipMemcpyDeviceToHost));

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

                    if (i != k) mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange;
                }
            }
        }

        timer.stop("  K block PP");
    }

    // K: S-D block

    if ((pair_inds_count_for_K_sd > 0) && (std::fabs(frac_exact_exchange) > 1.0e-13))
    {
        timer.start("  K block SD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sd + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_sd));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_sd, 1);

        // K: (SS|DS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSDS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|DP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSDP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|DS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPDS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|DP)
        //     *  *

        timer.start("    K block SPDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockSPDP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block SPDP");

        // K: (SS|DD)
        //     *  *

        timer.start("    K block SSDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockSSDD, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block SSDD");

        // K: (SD|DS)
        //     *  *

        timer.start("    K block SDDS");

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block SDDS");

        // K: (SP|DD)
        //     *  *

        timer.start("    K block SPDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockSPDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSPDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSPDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block SPDD");

        // K: (SD|DP)
        //     *  *

        timer.start("    K block SDDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block SDDP");

        // K: (SD|DD)
        //     *  *

        timer.start("    K block SDDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockSDDD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block SDDD");

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sd * sizeof(double), hipMemcpyDeviceToHost));

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
                mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto] += mat_K[ik] * k_coef_sph * (-1.0) * frac_exact_exchange;
            }
        }

        timer.stop("  K block SD");
    }

    // K: P-D block

    if ((pair_inds_count_for_K_pd > 0) && (std::fabs(frac_exact_exchange) > 1.0e-13))
    {
        timer.start("  K block PD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pd + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_pd));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_pd, 1);

        // K: (PS|DS)
        //     *  *

        timer.start("    K block PSDS");

        hipLaunchKernelGGL(gpu::computeExchangeFockPSDS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PSDS");

        // K: (PS|DP)
        //     *  *

        timer.start("    K block PSDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockPSDP, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PSDP");

        // K: (PP|DS)
        //     *  *

        timer.start("    K block PPDS");

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PPDS");

        // K: (PS|DD)
        //     *  *

        timer.start("    K block PSDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockPSDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPSDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPSDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PSDD");

        // K: (PD|DS)
        //     *  *

        timer.start("    K block PDDS");

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDS0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDS1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDS2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PDDS");

        // K: (PP|DP)
        //     *  *

        timer.start("    K block PPDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PPDP");

        // K: (PP|DD)
        //     *  *

        timer.start("    K block PPDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPPDD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PPDD");

        // K: (PD|DP)
        //     *  *

        timer.start("    K block PDDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDP3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDP4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PDDP");

        // K: (PD|DD)
        //     *  *

        timer.start("    K block PDDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD6, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD7, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD8, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD9, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD10, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD11, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD12, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD13, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD14, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD15, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD16, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockPDDD17, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block PDDD");

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pd * sizeof(double), hipMemcpyDeviceToHost));

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
                    mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange;
                }
            }
        }

        timer.stop("  K block PD");
    }

    // K: D-D block

    if ((pair_inds_count_for_K_dd > 0) && (std::fabs(frac_exact_exchange) > 1.0e-13))
    {
        timer.start("  K block DD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_dd + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_dd));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_dd, 1);

        // K: (DS|DS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDS, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (DS|DP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (DP|DS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDS0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDS1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDS2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        // K: (DS|DD)
        //     *  *

        timer.start("    K block DSDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDSDD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block DSDD");

        // K: (DD|DS)
        //     *  *

        timer.start("    K block DDDS");

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDS0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDS1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDS2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDS3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDS4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDS5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block DDDS");

        // K: (DP|DP)
        //     *  *

        timer.start("    K block DPDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDP3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDP4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDP5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block DPDP");

        // K: (DP|DD)
        //     *  *

        timer.start("    K block DPDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD6, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD7, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD8, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD9, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD10, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD11, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD12, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD13, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD14, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD15, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDPDD16, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block DPDD");

        // K: (DD|DP)
        //     *  *

        timer.start("    K block DDDP");

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP6, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP7, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP8, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP9, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP10, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP11, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP12, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP13, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP14, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP15, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP16, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP17, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP18, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDP19, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block DDDP");

        // K: (DD|DD)
        //     *  *

        timer.start("    K block DDDD");

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD0, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD1, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD2, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD3, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD4, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD5, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD6, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD7, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD8, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD9, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD10, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD11, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD12, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD13, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD14, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD15, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD16, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD17, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD18, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD19, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD20, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD21, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD22, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD23, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD24, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD25, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD26, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD27, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD28, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD29, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD30, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD31, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD32, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD33, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD34, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD35, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD36, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD37, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD38, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD39, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD40, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipLaunchKernelGGL(gpu::computeExchangeFockDDDD41, num_blocks, threads_per_block, 0, 0,
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
                           eri_threshold);

        hipSafe(hipDeviceSynchronize());

        timer.stop("    K block DDDD");

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_dd * sizeof(double), hipMemcpyDeviceToHost));

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

                    if (i != k) mat_Fock_omp[gpu_id].row(k_cgto_sph)[i_cgto_sph] += mat_K[ik] * coef_sph * (-1.0) * frac_exact_exchange;
                }
            }
        }

        timer.stop("  K block DD");
    }

    timer.stop("K computation");

    hipSafe(hipDeviceSynchronize());

    hipSafe(hipFree(d_boys_func_table));
    hipSafe(hipFree(d_boys_func_ft));

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_s_prim_aoinds));

    hipSafe(hipFree(d_p_prim_info));
    hipSafe(hipFree(d_p_prim_aoinds));

    hipSafe(hipFree(d_d_prim_info));
    hipSafe(hipFree(d_d_prim_aoinds));

    hipSafe(hipFree(d_mat_K));

    hipSafe(hipFree(d_pair_inds_i_for_K_ss));
    hipSafe(hipFree(d_pair_inds_k_for_K_ss));
    hipSafe(hipFree(d_pair_inds_i_for_K_sp));
    hipSafe(hipFree(d_pair_inds_k_for_K_sp));
    hipSafe(hipFree(d_pair_inds_i_for_K_sd));
    hipSafe(hipFree(d_pair_inds_k_for_K_sd));
    hipSafe(hipFree(d_pair_inds_i_for_K_pp));
    hipSafe(hipFree(d_pair_inds_k_for_K_pp));
    hipSafe(hipFree(d_pair_inds_i_for_K_pd));
    hipSafe(hipFree(d_pair_inds_k_for_K_pd));
    hipSafe(hipFree(d_pair_inds_i_for_K_dd));
    hipSafe(hipFree(d_pair_inds_k_for_K_dd));

    hipSafe(hipFree(d_mat_D_full_AO));

    hipSafe(hipFree(d_Q_K_ss));
    hipSafe(hipFree(d_Q_K_sp));
    hipSafe(hipFree(d_Q_K_ps));
    hipSafe(hipFree(d_Q_K_sd));
    hipSafe(hipFree(d_Q_K_ds));
    hipSafe(hipFree(d_Q_K_pp));
    hipSafe(hipFree(d_Q_K_pd));
    hipSafe(hipFree(d_Q_K_dp));
    hipSafe(hipFree(d_Q_K_dd));

    hipSafe(hipFree(d_D_inds_K_ss));
    hipSafe(hipFree(d_D_inds_K_sp));
    hipSafe(hipFree(d_D_inds_K_ps));
    hipSafe(hipFree(d_D_inds_K_sd));
    hipSafe(hipFree(d_D_inds_K_ds));
    hipSafe(hipFree(d_D_inds_K_pp));
    hipSafe(hipFree(d_D_inds_K_pd));
    hipSafe(hipFree(d_D_inds_K_dp));
    hipSafe(hipFree(d_D_inds_K_dd));

    hipSafe(hipFree(d_pair_displs_K_ss));
    hipSafe(hipFree(d_pair_displs_K_sp));
    hipSafe(hipFree(d_pair_displs_K_ps));
    hipSafe(hipFree(d_pair_displs_K_sd));
    hipSafe(hipFree(d_pair_displs_K_ds));
    hipSafe(hipFree(d_pair_displs_K_pp));
    hipSafe(hipFree(d_pair_displs_K_pd));
    hipSafe(hipFree(d_pair_displs_K_dp));
    hipSafe(hipFree(d_pair_displs_K_dd));

    hipSafe(hipFree(d_pair_counts_K_ss));
    hipSafe(hipFree(d_pair_counts_K_sp));
    hipSafe(hipFree(d_pair_counts_K_ps));
    hipSafe(hipFree(d_pair_counts_K_sd));
    hipSafe(hipFree(d_pair_counts_K_ds));
    hipSafe(hipFree(d_pair_counts_K_pp));
    hipSafe(hipFree(d_pair_counts_K_pd));
    hipSafe(hipFree(d_pair_counts_K_dp));
    hipSafe(hipFree(d_pair_counts_K_dd));

    hipSafe(hipFree(d_pair_data_K_ss));
    hipSafe(hipFree(d_pair_data_K_sp));
    hipSafe(hipFree(d_pair_data_K_ps));
    hipSafe(hipFree(d_pair_data_K_sd));
    hipSafe(hipFree(d_pair_data_K_ds));
    hipSafe(hipFree(d_pair_data_K_pp));
    hipSafe(hipFree(d_pair_data_K_pd));
    hipSafe(hipFree(d_pair_data_K_dp));
    hipSafe(hipFree(d_pair_data_K_dd));

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

    hipSafe(hipSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_A, *d_B;

    hipSafe(hipMalloc(&d_A, n * sizeof(double)));
    hipSafe(hipMalloc(&d_B, n * sizeof(double)));

    hipSafe(hipMemcpy(d_A, A, n * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, B, n * sizeof(double), hipMemcpyHostToDevice));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double dot_product;
    hipblasSafe(hipblasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));

    return dot_product;
}

auto
computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& data_pointers, const int64_t size) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeWeightedSum: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_X, *d_Y;

    hipSafe(hipMalloc(&d_X, n * sizeof(double)));
    hipSafe(hipMalloc(&d_Y, n * sizeof(double)));

    hipSafe(hipMemcpy(d_Y, weighted_data, n * sizeof(double), hipMemcpyHostToDevice));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        hipSafe(hipMemcpy(d_X, data_pointers[i], n * sizeof(double), hipMemcpyHostToDevice));

        hipblasSafe(hipblasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
    }

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipMemcpy(weighted_data, d_Y, n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_X));
    hipSafe(hipFree(d_Y));
}

auto
computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                   const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeErrorVector: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_A, *d_B, *d_C;

    hipSafe(hipMalloc(&d_A, nao * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_B, nao * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_C, nao * nao * sizeof(double)));

    hipSafe(hipMemcpy(d_A, F, nao * nao * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, D, nao * nao * sizeof(double), hipMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major

    // D^T F^T (=> FD)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));

    // S^T(FD)^T (=> FDS)
    hipSafe(hipMemcpy(d_A, S, nao * nao * sizeof(double), hipMemcpyHostToDevice));

    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));

    hipSafe(hipDeviceSynchronize());

    // FDS - (FDS)^T
    beta = -1.0;

    hipblasSafe(hipblasDgeam(handle, HIPBLAS_OP_N, HIPBLAS_OP_T, nao, nao, &alpha, d_B, nao, &beta, d_B, nao, d_C, nao));

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    hipSafe(hipMemcpy(d_X, X, nmo * nao * sizeof(double), hipMemcpyHostToDevice));

    auto op_X  = (trans_X == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? HIPBLAS_OP_T : HIPBLAS_OP_N;

    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    beta = 0.0;

    // let E == (FDS - SDF)
    // X^T E^T (=> EX)
    hipblasSafe(hipblasDgemm(handle, op_X, HIPBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_C, nao, &beta, d_Y, nmo));

    hipSafe(hipDeviceSynchronize());

    // (EX)^T X (=> X^T(E)X)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_C, nmo));

    hipSafe(hipMemcpy(errvec, d_C, nmo * nmo * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));
    hipSafe(hipFree(d_C));
}

auto
transformMatrix(double* transformed_F, const double* X, const double* F,
                const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::transformMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_F, *d_X, *d_Y;

    hipSafe(hipMalloc(&d_F, nao * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_X, nmo * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_Y, nmo * nao * sizeof(double)));

    hipSafe(hipMemcpy(d_F, F, nao * nao * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_X, X, nmo * nao * sizeof(double), hipMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

    auto op_X = (trans_X == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    auto op_XT = (trans_X == std::string("N")) ? HIPBLAS_OP_T : HIPBLAS_OP_N;
    auto lda_XT = (trans_X == std::string("N")) ? nao : nmo;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major

    // X^T F^T (=> FX)
    hipblasSafe(hipblasDgemm(handle, op_X, HIPBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_F, nao, &beta, d_Y, nmo));

    hipSafe(hipDeviceSynchronize());

    // (FX)^T X (=> X^T(F)X)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_F, nmo));

    hipSafe(hipMemcpy(transformed_F, d_F, nmo * nmo * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_F));
    hipSafe(hipFree(d_X));
    hipSafe(hipFree(d_Y));
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

    hipSafe(hipSetDevice(0));

    auto m = static_cast<int32_t>(m_inp);
    auto k = static_cast<int32_t>(k_inp);
    auto n = static_cast<int32_t>(n_inp);

    double *d_A, *d_B, *d_C;

    hipSafe(hipMalloc(&d_A, m * k * sizeof(double)));
    hipSafe(hipMalloc(&d_B, k * n * sizeof(double)));
    hipSafe(hipMalloc(&d_C, m * n * sizeof(double)));

    hipSafe(hipMemcpy(d_A, A, m * k * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, B, k * n * sizeof(double), hipMemcpyHostToDevice));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto op_A = (trans_A == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

    hipblasSafe(hipblasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));

    hipSafe(hipMemcpy(C, d_C, m * n * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));
    hipSafe(hipFree(d_C));
}

auto
diagonalizeMatrix(double* A, double* W, const int64_t nrows_A) -> void
{
    // TODO allow diagonalizeMatrix on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::diagonalizeMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    magmaSafe(magma_init());

    magma_setdevice(0);

    magma_int_t n = static_cast<magma_int_t>(nrows_A);
    //magma_int_t ldda = magma_roundup(n, 32);

    double *d_A;
    //hipSafe(hipMalloc(&d_A, n * ldda * sizeof(double)));
    //hipblasSafe(hipblasSetMatrix(n, n, sizeof(double), A, n, d_A, ldda));
    hipSafe(hipMalloc(&d_A, n * n * sizeof(double)));
    hipSafe(hipMemcpy(d_A, A, n * n * sizeof(double), hipMemcpyHostToDevice));

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n * n * sizeof(double)));
    hipSafe(hipHostMalloc(&work, lwork * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, liwork * sizeof(magma_int_t)));

    magma_int_t info;
    //magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, ldda, W, wA, n, work, lwork, iwork, liwork, &info);
    magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, n, W, wA, n, work, lwork, iwork, liwork, &info);

    if (info != 0)
    {
        std::stringstream ss;
        ss << "gpu::diagonalizeMatrix: (magma error) " << magma_strerror(info);
        errors::assertMsgCritical(false, ss.str());
    }

    //hipblasSafe(hipblasGetMatrix(n, n, sizeof(double), d_A, ldda, A, n));
    hipSafe(hipMemcpy(A, d_A, n * n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));

    hipSafe(hipHostFree(wA));
    hipSafe(hipHostFree(work));
    hipSafe(hipHostFree(iwork));

    magmaSafe(magma_finalize());
}

auto
diagonalizeMatrixMultiGPU(double* A, double* W, const int64_t nrows_A, const int64_t num_gpus_per_node) -> void
{
    // TODO allow diagonalizeMatrix on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::diagonalizeMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    magmaSafe(magma_init());

    magma_int_t n = static_cast<magma_int_t>(nrows_A);

    magma_int_t ngpu = static_cast<magma_int_t>(num_gpus_per_node);

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n * n * sizeof(double)));
    hipSafe(hipHostMalloc(&work, lwork * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, liwork * sizeof(magma_int_t)));

    magma_int_t info;
    magma_dsyevd_m(ngpu, MagmaVec, MagmaUpper, n, A, n, W, work, lwork, iwork, liwork, &info);

    if (info != 0)
    {
        std::stringstream ss;
        ss << "gpu::diagonalizeMatrixMultiGPU: (magma error) " << magma_strerror(info);
        errors::assertMsgCritical(false, ss.str());
    }

    hipSafe(hipHostFree(wA));
    hipSafe(hipHostFree(work));
    hipSafe(hipHostFree(iwork));

    magmaSafe(magma_finalize());
}

/*
auto
diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void
{
    hipSafe(hipSetDevice(0));

    auto n = static_cast<int32_t>(nrows_A);
    int32_t lwork;

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    hipSafe(hipMalloc(&d_A, n * n * sizeof(double)));
    hipSafe(hipMalloc(&d_D, n * sizeof(double)));
    hipSafe(hipMalloc(&d_info, sizeof(int32_t)));

    hipSafe(hipMemcpy(d_A, A, n * n * sizeof(double), hipMemcpyHostToDevice));

    hipsolverHandle_t handle;
    hipsolverSafe(hipsolverCreate(&handle));

    hipsolverSafe(hipsolverDsyevd_bufferSize(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    hipSafe(hipMalloc(&d_work, lwork));

    hipsolverSafe(hipsolverDsyevd(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    hipsolverSafe(hipsolverDestroy(handle));

    // TODO check h_info
    int32_t h_info;

    hipSafe(hipMemcpy(&h_info, d_info, sizeof(int32_t), hipMemcpyDeviceToHost));
    hipSafe(hipMemcpy(A, d_A, n * n * sizeof(double), hipMemcpyDeviceToHost));
    hipSafe(hipMemcpy(D, d_D, n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_D));
    hipSafe(hipFree(d_info));
    hipSafe(hipFree(d_work));
}
*/

}  // namespace gpu
