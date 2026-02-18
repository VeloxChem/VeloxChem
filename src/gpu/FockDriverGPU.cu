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
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ScreeningData.hpp"
#include "BoysFuncTable.hpp"
#include "ChunkedMemcpyGPU.hpp"
#include "FockDriverGPU.hpp"
#include "ElectricDipoleIntegrals.hpp"
#include "LinearMomentumIntegrals.hpp"
#include "AngularMomentumIntegrals.hpp"
#include "EriCoulomb.hpp"
#include "EriExchange.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "GpuDevices.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "OneElectronIntegrals.hpp"
#include "StringFormat.hpp"
#include "PrecisionCut.hpp"

namespace gpu {  // gpu namespace

static inline std::vector<float>
to_float_vec(const std::vector<double>& src)
{
    std::vector<float> dst(src.size());
    std::transform(src.begin(), src.end(), dst.begin(),
                   [](double x) { return static_cast<float>(x); });
    return dst;
}

__global__ void
zeroData(double* d_data, const uint32_t n)
{
    const uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < n) d_data[i] = 0.0;
}

auto
computeQMatrixOnGPU(const CMolecule& molecule,
                    const CMolecularBasis& basis,
                    const CScreeningData& screening) -> CDenseMatrix
{
    CGpuDevices gpu_devices;

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

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> Q_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        Q_omp[gpu_id] = CDenseMatrix(all_prim_count, all_prim_count);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();
    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    const auto data_boys_func_size = boys_func_table.size() + boys_func_ft.size();

    double* d_data_boys_func;
    gpuSafe(gpuMalloc(&d_data_boys_func, data_boys_func_size * sizeof(double)));

    double* d_boys_func_table = d_data_boys_func;
    double* d_boys_func_ft = d_boys_func_table + boys_func_table.size();

    gpu::chunkedMemcpyHostToDevice<double>(d_boys_func_table, boys_func_table.data(), boys_func_table.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size());

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

    // S, P, D gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<double>   d_prim_info(5 * d_prim_count);

    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);
    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);
    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    const auto data_spd_prim_info_size = s_prim_info.size() + p_prim_info.size() + d_prim_info.size();

    double* d_data_spd_prim_info;
    gpuSafe(gpuMalloc(&d_data_spd_prim_info, data_spd_prim_info_size * sizeof(double)));

    double* d_s_prim_info = d_data_spd_prim_info;
    double* d_p_prim_info = d_s_prim_info + s_prim_info.size();
    double* d_d_prim_info = d_p_prim_info + p_prim_info.size();

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

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

    gpuSafe(gpuMalloc(&d_mat_Q, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    gpuSafe(gpuDeviceSynchronize());

    Q_omp[gpu_id].zero();

    auto& mat_Q_omp = Q_omp[gpu_id];

    // compute Q

    // Q: SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSS<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpu::chunkedMemcpyDeviceToHost<double>(h_mat_Q.data(), d_mat_Q, ss_prim_pair_count_local);

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

        gpu::computeQMatrixSP<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(h_mat_Q.data(), d_mat_Q, sp_prim_pair_count_local);

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

        gpu::computeQMatrixSD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(h_mat_Q.data(), d_mat_Q, sd_prim_pair_count_local);

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

        gpu::computeQMatrixPP<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpu::chunkedMemcpyDeviceToHost<double>(h_mat_Q.data(), d_mat_Q, pp_prim_pair_count_local);

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

        gpu::computeQMatrixPD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(h_mat_Q.data(), d_mat_Q, pd_prim_pair_count_local);

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

        gpu::computeQMatrixDD<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpu::chunkedMemcpyDeviceToHost<double>(h_mat_Q.data(), d_mat_Q, dd_prim_pair_count_local);

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            mat_Q_omp.row(s_prim_count + p_prim_count * 3 + i)[s_prim_count + p_prim_count * 3 + j] = h_mat_Q[ij];

            if (i != j) mat_Q_omp.row(s_prim_count + p_prim_count * 3 + j)[s_prim_count + p_prim_count * 3 + i] = h_mat_Q[ij];
        }
    }

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_data_boys_func));
    gpuSafe(gpuFree(d_data_spd_prim_info));

    gpuSafe(gpuFree(d_mat_Q));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));
    }

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
computeMixedBasisOverlapIntegralsOnGPU(const CMolecule&       molecule,
                                       const CMolecularBasis& basis_1,
                                       const CMolecularBasis& basis_2,
                                       const int64_t          num_gpus_per_node,
                                       const int64_t          rank,
                                       const int64_t          nnodes) -> CDenseMatrix
{
    const auto gto_blocks_1 = gtofunc::makeGtoBlocks(basis_1, molecule);
    const auto naos_1 = gtofunc::getNumberOfAtomicOrbitals(gto_blocks_1);

    const auto gto_blocks_2 = gtofunc::makeGtoBlocks(basis_2, molecule);
    const auto naos_2 = gtofunc::getNumberOfAtomicOrbitals(gto_blocks_2);

    std::string err_f_orb("F-orbital is not yet supported");

    for (const auto& gto_block : gto_blocks_1)
    {
        const auto gto_ang = gto_block.getAngularMomentum();
        errors::assertMsgCritical(gto_ang <= 2, err_f_orb);
    }

    for (const auto& gto_block : gto_blocks_2)
    {
        const auto gto_ang = gto_block.getAngularMomentum();
        errors::assertMsgCritical(gto_ang <= 2, err_f_orb);
    }

    std::vector<CDenseMatrix> S_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        S_matrices[gpu_id] = CDenseMatrix(naos_1, naos_2);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_id));

    // GTOs blocks and number of AOs

    const auto gto_blocks_1 = gtofunc::makeGtoBlocks(basis_1, molecule);
    const auto naos_1 = gtofunc::getNumberOfAtomicOrbitals(gto_blocks_1);

    const auto gto_blocks_2 = gtofunc::makeGtoBlocks(basis_2, molecule);
    const auto naos_2 = gtofunc::getNumberOfAtomicOrbitals(gto_blocks_2);

    // gto blocks

    int64_t s_prim_count_1 = 0;
    int64_t p_prim_count_1 = 0;
    int64_t d_prim_count_1 = 0;

    for (const auto& gto_block : gto_blocks_1)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count_1 += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count_1 += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count_1 += npgtos * ncgtos;
    }

    int64_t s_prim_count_2 = 0;
    int64_t p_prim_count_2 = 0;
    int64_t d_prim_count_2 = 0;

    for (const auto& gto_block : gto_blocks_2)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count_2 += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count_2 += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count_2 += npgtos * ncgtos;
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p_1;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d_1;

    for (const auto& gto_block : gto_blocks_1)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p_1[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d_1[cart_ind] = sph_ind_coef;
            }
        }
    }

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p_2;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d_2;

    for (const auto& gto_block : gto_blocks_2)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p_2[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d_2[cart_ind] = sph_ind_coef;
            }
        }
    }

    // S gto block

    std::vector<double>   s_prim_info_1(5 * s_prim_count_1);
    std::vector<double>   s_prim_info_2(5 * s_prim_count_2);

    std::vector<uint32_t> s_prim_aoinds_1(1 * s_prim_count_1);
    std::vector<uint32_t> s_prim_aoinds_2(1 * s_prim_count_2);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info_1.data(), s_prim_aoinds_1.data(), s_prim_count_1, gto_blocks_1);
    gtoinfo::updatePrimitiveInfoForS(s_prim_info_2.data(), s_prim_aoinds_2.data(), s_prim_count_2, gto_blocks_2);

    // P gto block

    std::vector<double>   p_prim_info_1(5 * p_prim_count_1);
    std::vector<double>   p_prim_info_2(5 * p_prim_count_2);

    std::vector<uint32_t> p_prim_aoinds_1(3 * p_prim_count_1);
    std::vector<uint32_t> p_prim_aoinds_2(3 * p_prim_count_2);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info_1.data(), p_prim_aoinds_1.data(), p_prim_count_1, gto_blocks_1);
    gtoinfo::updatePrimitiveInfoForP(p_prim_info_2.data(), p_prim_aoinds_2.data(), p_prim_count_2, gto_blocks_2);

    // D gto block

    std::vector<double>   d_prim_info_1(5 * d_prim_count_1);
    std::vector<double>   d_prim_info_2(5 * d_prim_count_2);

    std::vector<uint32_t> d_prim_aoinds_1(6 * d_prim_count_1);
    std::vector<uint32_t> d_prim_aoinds_2(6 * d_prim_count_2);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info_1.data(), d_prim_aoinds_1.data(), d_prim_count_1, gto_blocks_1);
    gtoinfo::updatePrimitiveInfoForD(d_prim_info_2.data(), d_prim_aoinds_2.data(), d_prim_count_2, gto_blocks_2);

    double*   d_s_prim_info_1;
    double*   d_p_prim_info_1;
    double*   d_d_prim_info_1;

    gpuSafe(gpuMalloc(&d_s_prim_info_1, s_prim_info_1.size() * sizeof(double)));
    gpuSafe(gpuMalloc(&d_p_prim_info_1, p_prim_info_1.size() * sizeof(double)));
    gpuSafe(gpuMalloc(&d_d_prim_info_1, d_prim_info_1.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info_1, s_prim_info_1.data(), s_prim_info_1.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info_1, p_prim_info_1.data(), p_prim_info_1.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info_1, d_prim_info_1.data(), d_prim_info_1.size());

    double*   d_s_prim_info_2;
    double*   d_p_prim_info_2;
    double*   d_d_prim_info_2;

    gpuSafe(gpuMalloc(&d_s_prim_info_2, s_prim_info_2.size() * sizeof(double)));
    gpuSafe(gpuMalloc(&d_p_prim_info_2, p_prim_info_2.size() * sizeof(double)));
    gpuSafe(gpuMalloc(&d_d_prim_info_2, d_prim_info_2.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info_2, s_prim_info_2.data(), s_prim_info_2.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info_2, p_prim_info_2.data(), p_prim_info_2.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info_2, d_prim_info_2.data(), d_prim_info_2.size());

    // GTO block pairs

    std::vector<uint32_t> ss_first_inds_local;
    std::vector<uint32_t> sp_first_inds_local;
    std::vector<uint32_t> ps_first_inds_local;
    std::vector<uint32_t> sd_first_inds_local;
    std::vector<uint32_t> ds_first_inds_local;
    std::vector<uint32_t> pp_first_inds_local;
    std::vector<uint32_t> pd_first_inds_local;
    std::vector<uint32_t> dp_first_inds_local;
    std::vector<uint32_t> dd_first_inds_local;

    std::vector<uint32_t> ss_second_inds_local;
    std::vector<uint32_t> sp_second_inds_local;
    std::vector<uint32_t> ps_second_inds_local;
    std::vector<uint32_t> sd_second_inds_local;
    std::vector<uint32_t> ds_second_inds_local;
    std::vector<uint32_t> pp_second_inds_local;
    std::vector<uint32_t> pd_second_inds_local;
    std::vector<uint32_t> dp_second_inds_local;
    std::vector<uint32_t> dd_second_inds_local;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = gpu_rank; i < s_prim_count_1; i+=gpu_count)
    {
        // S-S gto block pair

        for (int64_t j = 0; j < s_prim_count_2; j++)
        {
            ss_first_inds_local.push_back(i);
            ss_second_inds_local.push_back(j);
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count_2; j++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                sp_first_inds_local.push_back(i);
                sp_second_inds_local.push_back(j * 3 + s);
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count_2; j++)
        {
            for (int64_t s = 0; s < 6; s++)
            {
                sd_first_inds_local.push_back(i);
                sd_second_inds_local.push_back(j * 6 + s);
            }
        }
    }

    // P-P gto block pair and P-D gto block pair

    for (int64_t i = gpu_rank; i < p_prim_count_1; i+=gpu_count)
    {
        // P-S gto block pair

        for (int64_t j = 0; j < s_prim_count_2; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                ps_first_inds_local.push_back(i * 3 + i_cart);
                ps_second_inds_local.push_back(j);
            }
        }

        // P-P gto block pair

        for (int64_t j = 0; j < p_prim_count_2; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 3; j_cart++)
                {
                    pp_first_inds_local.push_back(i * 3 + i_cart);
                    pp_second_inds_local.push_back(j * 3 + j_cart);
                }
            }
        }

        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count_2; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    pd_first_inds_local.push_back(i * 3 + i_cart);
                    pd_second_inds_local.push_back(j * 6 + j_cart);
                }
            }
        }
    }

    // D-D gto block pair

    for (int64_t i = gpu_rank; i < d_prim_count_1; i+=gpu_count)
    {
        // D-S gto block pair

        for (int64_t j = 0; j < s_prim_count_2; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                ds_first_inds_local.push_back(i * 6 + i_cart);
                ds_second_inds_local.push_back(j);
            }
        }

        // D-P gto block pair

        for (int64_t j = 0; j < p_prim_count_2; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 3; j_cart++)
                {
                    dp_first_inds_local.push_back(i * 6 + i_cart);
                    dp_second_inds_local.push_back(j * 3 + j_cart);
                }
            }
        }

        // D-D gto block pair

        for (int64_t j = 0; j < d_prim_count_2; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    dd_first_inds_local.push_back(i * 6 + i_cart);
                    dd_second_inds_local.push_back(j * 6 + j_cart);
                }
            }
        }
    }

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto ps_prim_pair_count_local = static_cast<int64_t>(ps_first_inds_local.size());
    const auto sd_prim_pair_count_local = static_cast<int64_t>(sd_first_inds_local.size());
    const auto ds_prim_pair_count_local = static_cast<int64_t>(ds_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());
    const auto pd_prim_pair_count_local = static_cast<int64_t>(pd_first_inds_local.size());
    const auto dp_prim_pair_count_local = static_cast<int64_t>(dp_first_inds_local.size());
    const auto dd_prim_pair_count_local = static_cast<int64_t>(dd_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local, ps_prim_pair_count_local,
                                                     sd_prim_pair_count_local, ds_prim_pair_count_local, pp_prim_pair_count_local,
                                                     pd_prim_pair_count_local, dp_prim_pair_count_local, dd_prim_pair_count_local});

    std::vector<double> mat_S(max_prim_pair_count_local);

    // S and T on device

    double *d_mat_S;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_ps_first_inds_local, *d_ps_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_ds_first_inds_local, *d_ds_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dp_first_inds_local, *d_dp_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    gpuSafe(gpuMalloc(&d_mat_S, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_ps_first_inds_local, ps_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ps_second_inds_local, ps_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_ds_first_inds_local, ds_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ds_second_inds_local, ds_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dp_first_inds_local, dp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dp_second_inds_local, dp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ps_first_inds_local, ps_first_inds_local.data(), ps_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ps_second_inds_local, ps_second_inds_local.data(), ps_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ds_first_inds_local, ds_first_inds_local.data(), ds_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ds_second_inds_local, ds_second_inds_local.data(), ds_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dp_first_inds_local, dp_first_inds_local.data(), dp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dp_second_inds_local, dp_second_inds_local.data(), dp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    S_matrices[gpu_id].zero();

    auto& mat_overlap = S_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeMixedBasisOverlapSS<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_s_prim_info_1,
                           static_cast<uint32_t>(s_prim_count_1),
                           d_s_prim_info_2,
                           static_cast<uint32_t>(s_prim_count_2),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, ss_prim_pair_count_local);

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds_1[i];
            const auto j_cgto = s_prim_aoinds_2[j];

            mat_overlap.row(i_cgto)[j_cgto] += mat_S[ij];
        }
    }

    // SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeMixedBasisOverlapSP<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_s_prim_info_1,
                           static_cast<uint32_t>(s_prim_count_1),
                           d_p_prim_info_2,
                           static_cast<uint32_t>(p_prim_count_2),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, sp_prim_pair_count_local);

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds_1[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds_2[(j / 3) + p_prim_count_2 * (j % 3)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_p_2[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_overlap.row(i_cgto)[j_cgto_sph] += mat_S[ij] * j_coef_sph;
            }
        }
    }

    // PS

    if (ps_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ps_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        // Note: use SP kernel to compute PS
        gpu::computeMixedBasisOverlapSP<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_s_prim_info_2,
                           static_cast<uint32_t>(s_prim_count_2),
                           d_p_prim_info_1,
                           static_cast<uint32_t>(p_prim_count_1),
                           d_ps_second_inds_local,
                           d_ps_first_inds_local,
                           static_cast<uint32_t>(ps_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, ps_prim_pair_count_local);

        for (int64_t ij = 0; ij < ps_prim_pair_count_local; ij++)
        {
            const auto i = ps_second_inds_local[ij];
            const auto j = ps_first_inds_local[ij];

            const auto i_cgto = s_prim_aoinds_2[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds_1[(j / 3) + p_prim_count_1 * (j % 3)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_p_1[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_overlap.row(j_cgto_sph)[i_cgto] += mat_S[ij] * j_coef_sph;
            }
        }
    }

    // SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeMixedBasisOverlapSD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_s_prim_info_1,
                           static_cast<uint32_t>(s_prim_count_1),
                           d_d_prim_info_2,
                           static_cast<uint32_t>(d_prim_count_2),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, sd_prim_pair_count_local);

        for (int64_t ij = 0; ij < sd_prim_pair_count_local; ij++)
        {
            const auto i = sd_first_inds_local[ij];
            const auto j = sd_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds_1[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = d_prim_aoinds_2[(j / 6) + d_prim_count_2 * (j % 6)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_d_2[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_overlap.row(i_cgto)[j_cgto_sph] += mat_S[ij] * j_coef_sph;
            }
        }
    }

    // DS

    if (ds_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ds_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        // Note: use SD kernel to compute DS
        gpu::computeMixedBasisOverlapSD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_s_prim_info_2,
                           static_cast<uint32_t>(s_prim_count_2),
                           d_d_prim_info_1,
                           static_cast<uint32_t>(d_prim_count_1),
                           d_ds_second_inds_local,
                           d_ds_first_inds_local,
                           static_cast<uint32_t>(ds_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, ds_prim_pair_count_local);

        for (int64_t ij = 0; ij < ds_prim_pair_count_local; ij++)
        {
            const auto i = ds_second_inds_local[ij];
            const auto j = ds_first_inds_local[ij];

            const auto i_cgto = s_prim_aoinds_2[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = d_prim_aoinds_1[(j / 6) + d_prim_count_1 * (j % 6)];

            // Cartesian to spherical
            for (const auto& j_cgto_sph_ind_coef : cart_sph_d_1[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                mat_overlap.row(j_cgto_sph)[i_cgto] += mat_S[ij] * j_coef_sph;
            }
        }
    }

    // PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeMixedBasisOverlapPP<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_p_prim_info_1,
                           static_cast<uint32_t>(p_prim_count_1),
                           d_p_prim_info_2,
                           static_cast<uint32_t>(p_prim_count_2),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, pp_prim_pair_count_local);

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds_1[(i / 3) + p_prim_count_1 * (i % 3)];
            const auto j_cgto = p_prim_aoinds_2[(j / 3) + p_prim_count_2 * (j % 3)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p_1[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_p_2[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(i_cgto_sph)[j_cgto_sph] += mat_S[ij] * coef_sph;
                }
            }
        }
    }

    // PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeMixedBasisOverlapPD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_p_prim_info_1,
                           static_cast<uint32_t>(p_prim_count_1),
                           d_d_prim_info_2,
                           static_cast<uint32_t>(d_prim_count_2),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, pd_prim_pair_count_local);

        for (int64_t ij = 0; ij < pd_prim_pair_count_local; ij++)
        {
            const auto i = pd_first_inds_local[ij];
            const auto j = pd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds_1[(i / 3) + p_prim_count_1 * (i % 3)];
            const auto j_cgto = d_prim_aoinds_2[(j / 6) + d_prim_count_2 * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p_1[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d_2[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(i_cgto_sph)[j_cgto_sph] += mat_S[ij] * coef_sph;
                }
            }
        }
    }

    // DP

    if (dp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        // Note: use PD kernel to compute DP
        gpu::computeMixedBasisOverlapPD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_p_prim_info_2,
                           static_cast<uint32_t>(p_prim_count_2),
                           d_d_prim_info_1,
                           static_cast<uint32_t>(d_prim_count_1),
                           d_dp_second_inds_local,
                           d_dp_first_inds_local,
                           static_cast<uint32_t>(dp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, dp_prim_pair_count_local);

        for (int64_t ij = 0; ij < dp_prim_pair_count_local; ij++)
        {
            const auto i = dp_second_inds_local[ij];
            const auto j = dp_first_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds_2[(i / 3) + p_prim_count_2 * (i % 3)];
            const auto j_cgto = d_prim_aoinds_1[(j / 6) + d_prim_count_1 * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_p_2[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d_1[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(j_cgto_sph)[i_cgto_sph] += mat_S[ij] * coef_sph;
                }
            }
        }
    }

    // DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeMixedBasisOverlapDD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_d_prim_info_1,
                           static_cast<uint32_t>(d_prim_count_1),
                           d_d_prim_info_2,
                           static_cast<uint32_t>(d_prim_count_2),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, dd_prim_pair_count_local);

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = d_prim_aoinds_1[(i / 6) + d_prim_count_1 * (i % 6)];
            const auto j_cgto = d_prim_aoinds_2[(j / 6) + d_prim_count_2 * (j % 6)];

            // Cartesian to spherical
            for (const auto& i_cgto_sph_ind_coef : cart_sph_d_1[i_cgto])
            {
                auto i_cgto_sph = i_cgto_sph_ind_coef.first;
                auto i_coef_sph = i_cgto_sph_ind_coef.second;

                for (const auto& j_cgto_sph_ind_coef : cart_sph_d_2[j_cgto])
                {
                    auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                    auto j_coef_sph = j_cgto_sph_ind_coef.second;

                    auto coef_sph = i_coef_sph * j_coef_sph;

                    mat_overlap.row(i_cgto_sph)[j_cgto_sph] += mat_S[ij] * coef_sph;
                }
            }
        }
    }

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_s_prim_info_1));
    gpuSafe(gpuFree(d_p_prim_info_1));
    gpuSafe(gpuFree(d_d_prim_info_1));

    gpuSafe(gpuFree(d_s_prim_info_2));
    gpuSafe(gpuFree(d_p_prim_info_2));
    gpuSafe(gpuFree(d_d_prim_info_2));

    gpuSafe(gpuFree(d_mat_S));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_ps_first_inds_local));
    gpuSafe(gpuFree(d_ps_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_ds_first_inds_local));
    gpuSafe(gpuFree(d_ds_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dp_first_inds_local));
    gpuSafe(gpuFree(d_dp_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));
    }

    CDenseMatrix S(naos_1, naos_2);

    S.zero();

    auto p_mat_S = S.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_overlap = S_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos_1 * naos_2; ind++)
        {
            p_mat_S[ind] += p_mat_overlap[ind];
        }
    }

    return S;
}

auto
computeOverlapAndKineticEnergyIntegralsOnGPU(const CMolecule& molecule,
                                             const CMolecularBasis& basis,
                                             const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

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

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

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

    gpuSafe(gpuMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

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

    gpuSafe(gpuMalloc(&d_mat_S, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_T, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    S_matrices[gpu_id].zero();
    T_matrices[gpu_id].zero();

    auto& mat_overlap = S_matrices[gpu_id];
    auto& mat_kinetic_energy = T_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeOverlapAndKineticEnergySS<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_T.data(), d_mat_T, ss_prim_pair_count_local);

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

        gpu::computeOverlapAndKineticEnergySP<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_T.data(), d_mat_T, sp_prim_pair_count_local);

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

        gpu::computeOverlapAndKineticEnergySD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_T.data(), d_mat_T, sd_prim_pair_count_local);

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

        gpu::computeOverlapAndKineticEnergyPP<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_T.data(), d_mat_T, pp_prim_pair_count_local);

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

        gpu::computeOverlapAndKineticEnergyPD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_T.data(), d_mat_T, pd_prim_pair_count_local);

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

        gpu::computeOverlapAndKineticEnergyDD<<<num_blocks, threads_per_block>>>(
                           d_mat_S,
                           d_mat_T,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_S.data(), d_mat_S, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_T.data(), d_mat_T, dd_prim_pair_count_local);

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

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_s_prim_info));
    gpuSafe(gpuFree(d_p_prim_info));
    gpuSafe(gpuFree(d_d_prim_info));

    gpuSafe(gpuFree(d_mat_S));
    gpuSafe(gpuFree(d_mat_T));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));
    }

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
computeNuclearPotentialIntegralsOnGPU(const CMolecule& molecule,
                                      const CMolecularBasis& basis,
                                      const CScreeningData& screening) -> CDenseMatrix
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
computePointChargesIntegralsOnGPU(const CMolecule& molecule,
                                  const CMolecularBasis& basis,
                                  const CScreeningData& screening,
                                  const double* points_info_ptr,
                                  const int64_t npoints) -> CDenseMatrix
{
    CGpuDevices gpu_devices;

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> V_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        V_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    gpuSafe(gpuMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_boys_func_table, boys_func_table.data(), boys_func_table.size());

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    gpuSafe(gpuMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size());

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

    gpuSafe(gpuMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

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

    gpuSafe(gpuMalloc(&d_mat_V, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    V_matrices[gpu_id].zero();

    auto& mat_nuclear_potential = V_matrices[gpu_id];

    double *d_points_info;

    gpuSafe(gpuMalloc(&d_points_info, npoints * 4 * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_points_info, points_info_ptr, npoints * 4);

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeNuclearPotentialSS<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_V.data(), d_mat_V, ss_prim_pair_count_local);

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

        gpu::computeNuclearPotentialSP<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_V.data(), d_mat_V, sp_prim_pair_count_local);

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

        gpu::computeNuclearPotentialSD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_V.data(), d_mat_V, sd_prim_pair_count_local);

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

        gpu::computeNuclearPotentialPP<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_V.data(), d_mat_V, pp_prim_pair_count_local);

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

        gpu::computeNuclearPotentialPD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_V.data(), d_mat_V, pd_prim_pair_count_local);

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

        gpu::computeNuclearPotentialDD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_V.data(), d_mat_V, dd_prim_pair_count_local);

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

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_boys_func_table));
    gpuSafe(gpuFree(d_boys_func_ft));

    gpuSafe(gpuFree(d_s_prim_info));
    gpuSafe(gpuFree(d_p_prim_info));
    gpuSafe(gpuFree(d_d_prim_info));

    gpuSafe(gpuFree(d_mat_V));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));

    gpuSafe(gpuFree(d_points_info));
    }

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
computeElectricDipoleIntegralsOnGPU(const CMolecule& molecule,
                                    const CMolecularBasis& basis,
                                    const std::vector<double>& origin,
                                    const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> MX_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> MY_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> MZ_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        MX_matrices[gpu_id] = CDenseMatrix(naos, naos);
        MY_matrices[gpu_id] = CDenseMatrix(naos, naos);
        MZ_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

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

    gpuSafe(gpuMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

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

    std::vector<double> mat_MX(max_prim_pair_count_local);
    std::vector<double> mat_MY(max_prim_pair_count_local);
    std::vector<double> mat_MZ(max_prim_pair_count_local);

    // MX/MY/MZ on device

    double *d_mat_MX, *d_mat_MY, *d_mat_MZ;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    gpuSafe(gpuMalloc(&d_mat_MX, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_MY, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_MZ, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    MX_matrices[gpu_id].zero();
    MY_matrices[gpu_id].zero();
    MZ_matrices[gpu_id].zero();

    auto& mat_mu_x = MX_matrices[gpu_id];
    auto& mat_mu_y = MY_matrices[gpu_id];
    auto& mat_mu_z = MZ_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeElectricDipoleSS<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, ss_prim_pair_count_local);

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_mu_x.row(i_cgto)[j_cgto] += mat_MX[ij];

            if (i != j) mat_mu_x.row(j_cgto)[i_cgto] += mat_MX[ij];

            mat_mu_y.row(i_cgto)[j_cgto] += mat_MY[ij];

            if (i != j) mat_mu_y.row(j_cgto)[i_cgto] += mat_MY[ij];

            mat_mu_z.row(i_cgto)[j_cgto] += mat_MZ[ij];

            if (i != j) mat_mu_z.row(j_cgto)[i_cgto] += mat_MZ[ij];
        }
    }

    // SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeElectricDipoleSP<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, sp_prim_pair_count_local);

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

                mat_mu_x.row(i_cgto)[j_cgto_sph] += mat_MX[ij] * j_coef_sph;
                mat_mu_x.row(j_cgto_sph)[i_cgto] += mat_MX[ij] * j_coef_sph;

                mat_mu_y.row(i_cgto)[j_cgto_sph] += mat_MY[ij] * j_coef_sph;
                mat_mu_y.row(j_cgto_sph)[i_cgto] += mat_MY[ij] * j_coef_sph;

                mat_mu_z.row(i_cgto)[j_cgto_sph] += mat_MZ[ij] * j_coef_sph;
                mat_mu_z.row(j_cgto_sph)[i_cgto] += mat_MZ[ij] * j_coef_sph;
            }
        }
    }

    // SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeElectricDipoleSD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, sd_prim_pair_count_local);

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

                mat_mu_x.row(i_cgto)[j_cgto_sph] += mat_MX[ij] * j_coef_sph;
                mat_mu_x.row(j_cgto_sph)[i_cgto] += mat_MX[ij] * j_coef_sph;

                mat_mu_y.row(i_cgto)[j_cgto_sph] += mat_MY[ij] * j_coef_sph;
                mat_mu_y.row(j_cgto_sph)[i_cgto] += mat_MY[ij] * j_coef_sph;

                mat_mu_z.row(i_cgto)[j_cgto_sph] += mat_MZ[ij] * j_coef_sph;
                mat_mu_z.row(j_cgto_sph)[i_cgto] += mat_MZ[ij] * j_coef_sph;
            }
        }
    }

    // PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeElectricDipolePP<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, pp_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;

                    if (i != j) mat_mu_x.row(j_cgto_sph)[i_cgto_sph] += mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;

                    if (i != j) mat_mu_y.row(j_cgto_sph)[i_cgto_sph] += mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;

                    if (i != j) mat_mu_z.row(j_cgto_sph)[i_cgto_sph] += mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    // PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeElectricDipolePD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, pd_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;
                    mat_mu_x.row(j_cgto_sph)[i_cgto_sph] += mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;
                    mat_mu_y.row(j_cgto_sph)[i_cgto_sph] += mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;
                    mat_mu_z.row(j_cgto_sph)[i_cgto_sph] += mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    // DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeElectricDipoleDD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, dd_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;

                    if (i != j) mat_mu_x.row(j_cgto_sph)[i_cgto_sph] += mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;

                    if (i != j) mat_mu_y.row(j_cgto_sph)[i_cgto_sph] += mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;

                    if (i != j) mat_mu_z.row(j_cgto_sph)[i_cgto_sph] += mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_s_prim_info));
    gpuSafe(gpuFree(d_p_prim_info));
    gpuSafe(gpuFree(d_d_prim_info));

    gpuSafe(gpuFree(d_mat_MX));
    gpuSafe(gpuFree(d_mat_MY));
    gpuSafe(gpuFree(d_mat_MZ));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));
    }

    std::vector<CDenseMatrix> edip_matrices(3);

    edip_matrices[0] = CDenseMatrix(naos, naos);
    edip_matrices[1] = CDenseMatrix(naos, naos);
    edip_matrices[2] = CDenseMatrix(naos, naos);

    edip_matrices[0].zero();
    edip_matrices[1].zero();
    edip_matrices[2].zero();

    auto p_mat_MX = edip_matrices[0].values();
    auto p_mat_MY = edip_matrices[1].values();
    auto p_mat_MZ = edip_matrices[2].values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_mu_x = MX_matrices[gpu_id].values();
        auto p_mat_mu_y = MY_matrices[gpu_id].values();
        auto p_mat_mu_z = MZ_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_MX[ind] += p_mat_mu_x[ind];
            p_mat_MY[ind] += p_mat_mu_y[ind];
            p_mat_MZ[ind] += p_mat_mu_z[ind];
        }
    }

    return edip_matrices;
}

auto
computeLinearMomentumIntegralsOnGPU(const CMolecule& molecule,
                                    const CMolecularBasis& basis,
                                    const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> MX_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> MY_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> MZ_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        MX_matrices[gpu_id] = CDenseMatrix(naos, naos);
        MY_matrices[gpu_id] = CDenseMatrix(naos, naos);
        MZ_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

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

    gpuSafe(gpuMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

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

    std::vector<double> mat_MX(max_prim_pair_count_local);
    std::vector<double> mat_MY(max_prim_pair_count_local);
    std::vector<double> mat_MZ(max_prim_pair_count_local);

    // MX/MY/MZ on device

    double *d_mat_MX, *d_mat_MY, *d_mat_MZ;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    gpuSafe(gpuMalloc(&d_mat_MX, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_MY, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_MZ, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    MX_matrices[gpu_id].zero();
    MY_matrices[gpu_id].zero();
    MZ_matrices[gpu_id].zero();

    auto& mat_mu_x = MX_matrices[gpu_id];
    auto& mat_mu_y = MY_matrices[gpu_id];
    auto& mat_mu_z = MZ_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeLinearMomentumSS<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, ss_prim_pair_count_local);

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_mu_x.row(i_cgto)[j_cgto] += mat_MX[ij];

            if (i != j) mat_mu_x.row(j_cgto)[i_cgto] -= mat_MX[ij];

            mat_mu_y.row(i_cgto)[j_cgto] += mat_MY[ij];

            if (i != j) mat_mu_y.row(j_cgto)[i_cgto] -= mat_MY[ij];

            mat_mu_z.row(i_cgto)[j_cgto] += mat_MZ[ij];

            if (i != j) mat_mu_z.row(j_cgto)[i_cgto] -= mat_MZ[ij];
        }
    }

    // SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeLinearMomentumSP<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, sp_prim_pair_count_local);

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

                mat_mu_x.row(i_cgto)[j_cgto_sph] += mat_MX[ij] * j_coef_sph;
                mat_mu_x.row(j_cgto_sph)[i_cgto] -= mat_MX[ij] * j_coef_sph;

                mat_mu_y.row(i_cgto)[j_cgto_sph] += mat_MY[ij] * j_coef_sph;
                mat_mu_y.row(j_cgto_sph)[i_cgto] -= mat_MY[ij] * j_coef_sph;

                mat_mu_z.row(i_cgto)[j_cgto_sph] += mat_MZ[ij] * j_coef_sph;
                mat_mu_z.row(j_cgto_sph)[i_cgto] -= mat_MZ[ij] * j_coef_sph;
            }
        }
    }

    // SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeLinearMomentumSD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, sd_prim_pair_count_local);

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

                mat_mu_x.row(i_cgto)[j_cgto_sph] += mat_MX[ij] * j_coef_sph;
                mat_mu_x.row(j_cgto_sph)[i_cgto] -= mat_MX[ij] * j_coef_sph;

                mat_mu_y.row(i_cgto)[j_cgto_sph] += mat_MY[ij] * j_coef_sph;
                mat_mu_y.row(j_cgto_sph)[i_cgto] -= mat_MY[ij] * j_coef_sph;

                mat_mu_z.row(i_cgto)[j_cgto_sph] += mat_MZ[ij] * j_coef_sph;
                mat_mu_z.row(j_cgto_sph)[i_cgto] -= mat_MZ[ij] * j_coef_sph;
            }
        }
    }

    // PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeLinearMomentumPP<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, pp_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;

                    if (i != j) mat_mu_x.row(j_cgto_sph)[i_cgto_sph] -= mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;

                    if (i != j) mat_mu_y.row(j_cgto_sph)[i_cgto_sph] -= mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;

                    if (i != j) mat_mu_z.row(j_cgto_sph)[i_cgto_sph] -= mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    // PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeLinearMomentumPD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, pd_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;
                    mat_mu_x.row(j_cgto_sph)[i_cgto_sph] -= mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;
                    mat_mu_y.row(j_cgto_sph)[i_cgto_sph] -= mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;
                    mat_mu_z.row(j_cgto_sph)[i_cgto_sph] -= mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    // DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeLinearMomentumDD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, dd_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;

                    if (i != j) mat_mu_x.row(j_cgto_sph)[i_cgto_sph] -= mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;

                    if (i != j) mat_mu_y.row(j_cgto_sph)[i_cgto_sph] -= mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;

                    if (i != j) mat_mu_z.row(j_cgto_sph)[i_cgto_sph] -= mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_s_prim_info));
    gpuSafe(gpuFree(d_p_prim_info));
    gpuSafe(gpuFree(d_d_prim_info));

    gpuSafe(gpuFree(d_mat_MX));
    gpuSafe(gpuFree(d_mat_MY));
    gpuSafe(gpuFree(d_mat_MZ));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));
    }

    std::vector<CDenseMatrix> lmom_matrices(3);

    lmom_matrices[0] = CDenseMatrix(naos, naos);
    lmom_matrices[1] = CDenseMatrix(naos, naos);
    lmom_matrices[2] = CDenseMatrix(naos, naos);

    lmom_matrices[0].zero();
    lmom_matrices[1].zero();
    lmom_matrices[2].zero();

    auto p_mat_MX = lmom_matrices[0].values();
    auto p_mat_MY = lmom_matrices[1].values();
    auto p_mat_MZ = lmom_matrices[2].values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_mu_x = MX_matrices[gpu_id].values();
        auto p_mat_mu_y = MY_matrices[gpu_id].values();
        auto p_mat_mu_z = MZ_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_MX[ind] += p_mat_mu_x[ind];
            p_mat_MY[ind] += p_mat_mu_y[ind];
            p_mat_MZ[ind] += p_mat_mu_z[ind];
        }
    }

    return lmom_matrices;
}

auto
computeAngularMomentumIntegralsOnGPU(const CMolecule& molecule,
                                     const CMolecularBasis& basis,
                                     const std::vector<double>& origin,
                                     const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> MX_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> MY_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> MZ_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        MX_matrices[gpu_id] = CDenseMatrix(naos, naos);
        MY_matrices[gpu_id] = CDenseMatrix(naos, naos);
        MZ_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

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

    gpuSafe(gpuMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

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

    std::vector<double> mat_MX(max_prim_pair_count_local);
    std::vector<double> mat_MY(max_prim_pair_count_local);
    std::vector<double> mat_MZ(max_prim_pair_count_local);

    // MX/MY/MZ on device

    double *d_mat_MX, *d_mat_MY, *d_mat_MZ;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_sd_first_inds_local, *d_sd_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;
    uint32_t *d_pd_first_inds_local, *d_pd_second_inds_local;
    uint32_t *d_dd_first_inds_local, *d_dd_second_inds_local;

    gpuSafe(gpuMalloc(&d_mat_MX, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_MY, max_prim_pair_count_local * sizeof(double)));
    gpuSafe(gpuMalloc(&d_mat_MZ, max_prim_pair_count_local * sizeof(double)));

    gpuSafe(gpuMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_sd_first_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_sd_second_inds_local, sd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_pd_first_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_pd_second_inds_local, pd_prim_pair_count_local * sizeof(uint32_t)));

    gpuSafe(gpuMalloc(&d_dd_first_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));
    gpuSafe(gpuMalloc(&d_dd_second_inds_local, dd_prim_pair_count_local * sizeof(uint32_t)));

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local);

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local);
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local);

    MX_matrices[gpu_id].zero();
    MY_matrices[gpu_id].zero();
    MZ_matrices[gpu_id].zero();

    auto& mat_mu_x = MX_matrices[gpu_id];
    auto& mat_mu_y = MY_matrices[gpu_id];
    auto& mat_mu_z = MZ_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    // SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeAngularMomentumSS<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds_local,
                           d_ss_second_inds_local,
                           static_cast<uint32_t>(ss_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, ss_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, ss_prim_pair_count_local);

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_mu_x.row(i_cgto)[j_cgto] += mat_MX[ij];

            if (i != j) mat_mu_x.row(j_cgto)[i_cgto] -= mat_MX[ij];

            mat_mu_y.row(i_cgto)[j_cgto] += mat_MY[ij];

            if (i != j) mat_mu_y.row(j_cgto)[i_cgto] -= mat_MY[ij];

            mat_mu_z.row(i_cgto)[j_cgto] += mat_MZ[ij];

            if (i != j) mat_mu_z.row(j_cgto)[i_cgto] -= mat_MZ[ij];
        }
    }

    // SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeAngularMomentumSP<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds_local,
                           d_sp_second_inds_local,
                           static_cast<uint32_t>(sp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, sp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, sp_prim_pair_count_local);

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

                mat_mu_x.row(i_cgto)[j_cgto_sph] += mat_MX[ij] * j_coef_sph;
                mat_mu_x.row(j_cgto_sph)[i_cgto] -= mat_MX[ij] * j_coef_sph;

                mat_mu_y.row(i_cgto)[j_cgto_sph] += mat_MY[ij] * j_coef_sph;
                mat_mu_y.row(j_cgto_sph)[i_cgto] -= mat_MY[ij] * j_coef_sph;

                mat_mu_z.row(i_cgto)[j_cgto_sph] += mat_MZ[ij] * j_coef_sph;
                mat_mu_z.row(j_cgto_sph)[i_cgto] -= mat_MZ[ij] * j_coef_sph;
            }
        }
    }

    // SD

    if (sd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeAngularMomentumSD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds_local,
                           d_sd_second_inds_local,
                           static_cast<uint32_t>(sd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, sd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, sd_prim_pair_count_local);

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

                mat_mu_x.row(i_cgto)[j_cgto_sph] += mat_MX[ij] * j_coef_sph;
                mat_mu_x.row(j_cgto_sph)[i_cgto] -= mat_MX[ij] * j_coef_sph;

                mat_mu_y.row(i_cgto)[j_cgto_sph] += mat_MY[ij] * j_coef_sph;
                mat_mu_y.row(j_cgto_sph)[i_cgto] -= mat_MY[ij] * j_coef_sph;

                mat_mu_z.row(i_cgto)[j_cgto_sph] += mat_MZ[ij] * j_coef_sph;
                mat_mu_z.row(j_cgto_sph)[i_cgto] -= mat_MZ[ij] * j_coef_sph;
            }
        }
    }

    // PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeAngularMomentumPP<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds_local,
                           d_pp_second_inds_local,
                           static_cast<uint32_t>(pp_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, pp_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, pp_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;

                    if (i != j) mat_mu_x.row(j_cgto_sph)[i_cgto_sph] -= mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;

                    if (i != j) mat_mu_y.row(j_cgto_sph)[i_cgto_sph] -= mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;

                    if (i != j) mat_mu_z.row(j_cgto_sph)[i_cgto_sph] -= mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    // PD

    if (pd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeAngularMomentumPD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds_local,
                           d_pd_second_inds_local,
                           static_cast<uint32_t>(pd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, pd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, pd_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;
                    mat_mu_x.row(j_cgto_sph)[i_cgto_sph] -= mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;
                    mat_mu_y.row(j_cgto_sph)[i_cgto_sph] -= mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;
                    mat_mu_z.row(j_cgto_sph)[i_cgto_sph] -= mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    // DD

    if (dd_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeAngularMomentumDD<<<num_blocks, threads_per_block>>>(
                           d_mat_MX,
                           d_mat_MY,
                           d_mat_MZ,
                           origin[0],
                           origin[1],
                           origin[2],
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds_local,
                           d_dd_second_inds_local,
                           static_cast<uint32_t>(dd_prim_pair_count_local));

        gpu::chunkedMemcpyDeviceToHost<double>(mat_MX.data(), d_mat_MX, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MY.data(), d_mat_MY, dd_prim_pair_count_local);
        gpu::chunkedMemcpyDeviceToHost<double>(mat_MZ.data(), d_mat_MZ, dd_prim_pair_count_local);

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

                    mat_mu_x.row(i_cgto_sph)[j_cgto_sph] += mat_MX[ij] * coef_sph;

                    if (i != j) mat_mu_x.row(j_cgto_sph)[i_cgto_sph] -= mat_MX[ij] * coef_sph;

                    mat_mu_y.row(i_cgto_sph)[j_cgto_sph] += mat_MY[ij] * coef_sph;

                    if (i != j) mat_mu_y.row(j_cgto_sph)[i_cgto_sph] -= mat_MY[ij] * coef_sph;

                    mat_mu_z.row(i_cgto_sph)[j_cgto_sph] += mat_MZ[ij] * coef_sph;

                    if (i != j) mat_mu_z.row(j_cgto_sph)[i_cgto_sph] -= mat_MZ[ij] * coef_sph;
                }
            }
        }
    }

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_s_prim_info));
    gpuSafe(gpuFree(d_p_prim_info));
    gpuSafe(gpuFree(d_d_prim_info));

    gpuSafe(gpuFree(d_mat_MX));
    gpuSafe(gpuFree(d_mat_MY));
    gpuSafe(gpuFree(d_mat_MZ));

    gpuSafe(gpuFree(d_ss_first_inds_local));
    gpuSafe(gpuFree(d_ss_second_inds_local));
    gpuSafe(gpuFree(d_sp_first_inds_local));
    gpuSafe(gpuFree(d_sp_second_inds_local));
    gpuSafe(gpuFree(d_sd_first_inds_local));
    gpuSafe(gpuFree(d_sd_second_inds_local));
    gpuSafe(gpuFree(d_pp_first_inds_local));
    gpuSafe(gpuFree(d_pp_second_inds_local));
    gpuSafe(gpuFree(d_pd_first_inds_local));
    gpuSafe(gpuFree(d_pd_second_inds_local));
    gpuSafe(gpuFree(d_dd_first_inds_local));
    gpuSafe(gpuFree(d_dd_second_inds_local));
    }

    std::vector<CDenseMatrix> amom_matrices(3);

    amom_matrices[0] = CDenseMatrix(naos, naos);
    amom_matrices[1] = CDenseMatrix(naos, naos);
    amom_matrices[2] = CDenseMatrix(naos, naos);

    amom_matrices[0].zero();
    amom_matrices[1].zero();
    amom_matrices[2].zero();

    auto p_mat_MX = amom_matrices[0].values();
    auto p_mat_MY = amom_matrices[1].values();
    auto p_mat_MZ = amom_matrices[2].values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_mu_x = MX_matrices[gpu_id].values();
        auto p_mat_mu_y = MY_matrices[gpu_id].values();
        auto p_mat_mu_z = MZ_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_MX[ind] += p_mat_mu_x[ind];
            p_mat_MY[ind] += p_mat_mu_y[ind];
            p_mat_MZ[ind] += p_mat_mu_z[ind];
        }
    }

    return amom_matrices;
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
                 const std::vector<double>& frac_exact_exchange_values,
                 const std::vector<double>& omega_values,
                 const std::string& flag_K,
                 const double       eri_threshold,
                 const double       prelink_threshold,
                 const int32_t*     Q_prime_row_ptr,
                 const int32_t*     Q_prime_col_ptr,
                 const int32_t      Q_prime_ind_count,
                 CScreeningData&    screening) -> CDenseMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    // TODO sanity check for flag_K: SYMM or ANTISYMM

    CGpuDevices gpu_devices;

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();

    gpuSafe(gpuSetDevice(0));

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

    // preLinK
    // J. Chem. Phys. 138, 134114 (2013)

    timer.start("Prep. preLinK 1");

    screening.form_pair_inds_for_K(s_prim_count, p_prim_count, d_prim_count,
                                   Q_prime_row_ptr, Q_prime_col_ptr, Q_prime_ind_count);

    timer.stop("Prep. preLinK 1");

    timer.start("Prep. preLinK 2");

    // max densities are needed by exchange Fock
    screening.findMaxDensities(s_prim_count, p_prim_count, d_prim_count,
                               s_prim_aoinds, p_prim_aoinds, d_prim_aoinds,
                               cart_naos, cart_dens_ptr);

    timer.stop("Prep. preLinK 2");

    timer.start("Prep. Fockmat");

    std::string errnaos("gpu::computeFockGPU: Inconsistent number of AOs");
    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    std::vector<CDenseMatrix> mat_Fock_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Fock_omp[gpu_id] = CDenseMatrix(naos, naos);
    }

    screening.initGpuTimers(num_gpus_per_node);

    timer.stop("Prep. Fockmat");

    timer.start("Compute Fockmat");

    std::vector<CMultiTimer> omptimers(nthreads);

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    auto gpu_id = thread_id;

    gpuSafe(gpuSetDevice(gpu_id));

    omptimers[thread_id].start("Boys func. prep.");

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();
    const auto boys_func_ft    = boysfunc::getBoysFuncFactors();

    const auto data_boys_func_size = boys_func_table.size() + boys_func_ft.size();

    double* d_data_boys_func;
    gpuSafe(gpuMalloc(&d_data_boys_func, data_boys_func_size * sizeof(double)));

    double* d_boys_func_table = d_data_boys_func;
    double* d_boys_func_ft = d_boys_func_table + boys_func_table.size();

    gpu::chunkedMemcpyHostToDevice<double>(d_boys_func_table, boys_func_table.data(), boys_func_table.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size());

    // computeFockOnGPU boys
    // FP32 Boys function

    std::vector<float> boys_func_table_f;
    std::vector<float> boys_func_ft_f;

    boys_func_table_f.resize(boys_func_table.size());
    boys_func_ft_f.resize(boys_func_ft.size());

    for (size_t i = 0; i < boys_func_table.size(); ++i)
        boys_func_table_f[i] = static_cast<float>(boys_func_table[i]);

    for (size_t i = 0; i < boys_func_ft.size(); ++i)
        boys_func_ft_f[i] = static_cast<float>(boys_func_ft[i]);

    float* d_data_boys_func_f;
    gpuSafe(gpuMalloc(&d_data_boys_func_f,
        (boys_func_table_f.size() + boys_func_ft_f.size()) * sizeof(float)));

    float* d_boys_func_table_f = d_data_boys_func_f;
    float* d_boys_func_ft_f    = d_boys_func_table_f + boys_func_table_f.size();

    gpu::chunkedMemcpyHostToDevice<float>(
    d_boys_func_table_f,
    boys_func_table_f.data(),
    boys_func_table_f.size());

    gpu::chunkedMemcpyHostToDevice<float>(
        d_boys_func_ft_f,
        boys_func_ft_f.data(),
        boys_func_ft_f.size());

    omptimers[thread_id].stop("Boys func. prep.");

    omptimers[thread_id].start("GTO block prep.");

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

    // S, P and D gto blocks

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<double>   d_prim_info(5 * d_prim_count);

    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);
    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);
    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    const auto data_spd_prim_info_size = s_prim_info.size() + p_prim_info.size() + d_prim_info.size();
    const auto data_spd_prim_aoinds_size = s_prim_aoinds.size() + p_prim_aoinds.size() + d_prim_aoinds.size();

    double* d_data_spd_prim_info;
    uint32_t* d_data_spd_prim_aoinds;

    gpuSafe(gpuMalloc(&d_data_spd_prim_info, data_spd_prim_info_size * sizeof(double)));
    gpuSafe(gpuMalloc(&d_data_spd_prim_aoinds, data_spd_prim_aoinds_size * sizeof(uint32_t)));

    double*   d_s_prim_info = d_data_spd_prim_info;
    double*   d_p_prim_info = d_s_prim_info + s_prim_info.size();
    double*   d_d_prim_info = d_p_prim_info + p_prim_info.size();

    uint32_t* d_s_prim_aoinds = d_data_spd_prim_aoinds;
    uint32_t* d_p_prim_aoinds = d_s_prim_aoinds + s_prim_aoinds.size();
    uint32_t* d_d_prim_aoinds = d_p_prim_aoinds + p_prim_aoinds.size();

    gpu::chunkedMemcpyHostToDevice<double>(d_s_prim_info, s_prim_info.data(), s_prim_info.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_p_prim_info, p_prim_info.data(), p_prim_info.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_d_prim_info, d_prim_info.data(), d_prim_info.size());

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_d_prim_aoinds, d_prim_aoinds.data(), d_prim_aoinds.size());

    omptimers[thread_id].stop("GTO block prep.");

#pragma omp barrier

    // GTO block pairs

    omptimers[thread_id].start("J prep.");

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

    size_t data_J_double_size = 0;

    data_J_double_size += static_cast<size_t>(max_prim_pair_count);         // d_mat_D
    data_J_double_size += static_cast<size_t>(max_prim_pair_count_local);   // d_mat_J
    data_J_double_size += static_cast<size_t>(ss_prim_pair_count);          // d_ss_mat_Q
    data_J_double_size += static_cast<size_t>(sp_prim_pair_count);          // d_sp_mat_Q
    data_J_double_size += static_cast<size_t>(sd_prim_pair_count);          // d_sd_mat_Q
    data_J_double_size += static_cast<size_t>(pp_prim_pair_count);          // d_pp_mat_Q
    data_J_double_size += static_cast<size_t>(pd_prim_pair_count);          // d_pd_mat_Q
    data_J_double_size += static_cast<size_t>(dd_prim_pair_count);          // d_dd_mat_Q
    data_J_double_size += ss_pair_data.size();                              // d_ss_pair_data
    data_J_double_size += sp_pair_data.size();                              // d_sp_pair_data
    data_J_double_size += sd_pair_data.size();                              // d_sd_pair_data
    data_J_double_size += pp_pair_data.size();                              // d_pp_pair_data
    data_J_double_size += pd_pair_data.size();                              // d_pd_pair_data
    data_J_double_size += dd_pair_data.size();                              // d_dd_pair_data
    data_J_double_size += static_cast<size_t>(ss_prim_pair_count_local);    // d_ss_mat_Q_local
    data_J_double_size += static_cast<size_t>(sp_prim_pair_count_local);    // d_sp_mat_Q_local
    data_J_double_size += static_cast<size_t>(sd_prim_pair_count_local);    // d_sd_mat_Q_local
    data_J_double_size += static_cast<size_t>(pp_prim_pair_count_local);    // d_pp_mat_Q_local
    data_J_double_size += static_cast<size_t>(pd_prim_pair_count_local);    // d_pd_mat_Q_local
    data_J_double_size += static_cast<size_t>(dd_prim_pair_count_local);    // d_dd_mat_Q_local
    data_J_double_size += ss_pair_data_local.size();                        // d_ss_pair_data_local
    data_J_double_size += sp_pair_data_local.size();                        // d_sp_pair_data_local
    data_J_double_size += sd_pair_data_local.size();                        // d_sd_pair_data_local
    data_J_double_size += pp_pair_data_local.size();                        // d_pp_pair_data_local
    data_J_double_size += pd_pair_data_local.size();                        // d_pd_pair_data_local
    data_J_double_size += dd_pair_data_local.size();                        // d_dd_pair_data_local

    size_t data_J_uint32_size = 0;

    data_J_uint32_size += static_cast<size_t>(ss_prim_pair_count);          // d_ss_first_inds
    data_J_uint32_size += static_cast<size_t>(ss_prim_pair_count);          // d_ss_second_inds
    data_J_uint32_size += static_cast<size_t>(sp_prim_pair_count);          // d_sp_first_inds 
    data_J_uint32_size += static_cast<size_t>(sp_prim_pair_count);          // d_sp_second_inds
    data_J_uint32_size += static_cast<size_t>(sd_prim_pair_count);          // d_sd_first_inds 
    data_J_uint32_size += static_cast<size_t>(sd_prim_pair_count);          // d_sd_second_inds
    data_J_uint32_size += static_cast<size_t>(pp_prim_pair_count);          // d_pp_first_inds 
    data_J_uint32_size += static_cast<size_t>(pp_prim_pair_count);          // d_pp_second_inds
    data_J_uint32_size += static_cast<size_t>(pd_prim_pair_count);          // d_pd_first_inds 
    data_J_uint32_size += static_cast<size_t>(pd_prim_pair_count);          // d_pd_second_inds
    data_J_uint32_size += static_cast<size_t>(dd_prim_pair_count);          // d_dd_first_inds 
    data_J_uint32_size += static_cast<size_t>(dd_prim_pair_count);          // d_dd_second_inds
    data_J_uint32_size += static_cast<size_t>(ss_prim_pair_count_local);    // d_ss_first_inds_local 
    data_J_uint32_size += static_cast<size_t>(ss_prim_pair_count_local);    // d_ss_second_inds_local
    data_J_uint32_size += static_cast<size_t>(sp_prim_pair_count_local);    // d_sp_first_inds_local 
    data_J_uint32_size += static_cast<size_t>(sp_prim_pair_count_local);    // d_sp_second_inds_local
    data_J_uint32_size += static_cast<size_t>(sd_prim_pair_count_local);    // d_sd_first_inds_local 
    data_J_uint32_size += static_cast<size_t>(sd_prim_pair_count_local);    // d_sd_second_inds_local
    data_J_uint32_size += static_cast<size_t>(pp_prim_pair_count_local);    // d_pp_first_inds_local 
    data_J_uint32_size += static_cast<size_t>(pp_prim_pair_count_local);    // d_pp_second_inds_local
    data_J_uint32_size += static_cast<size_t>(pd_prim_pair_count_local);    // d_pd_first_inds_local 
    data_J_uint32_size += static_cast<size_t>(pd_prim_pair_count_local);    // d_pd_second_inds_local
    data_J_uint32_size += static_cast<size_t>(dd_prim_pair_count_local);    // d_dd_first_inds_local 
    data_J_uint32_size += static_cast<size_t>(dd_prim_pair_count_local);    // d_dd_second_inds_local

    // J data on device

#pragma omp barrier

    screening.resize_devptr_double(gpu_id, data_J_double_size);
    screening.resize_devptr_uint32(gpu_id, data_J_uint32_size);

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    auto d_data_J_double = screening.get_devptr_double(gpu_id);
    auto d_data_J_uint32 = screening.get_devptr_uint32(gpu_id);

    double* d_mat_D              = d_data_J_double;
    double* d_mat_J              = d_mat_D + max_prim_pair_count;
    double* d_ss_mat_Q           = d_mat_J + max_prim_pair_count_local;
    double* d_sp_mat_Q           = d_ss_mat_Q + ss_prim_pair_count;
    double* d_sd_mat_Q           = d_sp_mat_Q + sp_prim_pair_count;
    double* d_pp_mat_Q           = d_sd_mat_Q + sd_prim_pair_count;
    double* d_pd_mat_Q           = d_pp_mat_Q + pp_prim_pair_count;
    double* d_dd_mat_Q           = d_pd_mat_Q + pd_prim_pair_count;
    double* d_ss_pair_data       = d_dd_mat_Q + dd_prim_pair_count;
    double* d_sp_pair_data       = d_ss_pair_data + ss_pair_data.size();
    double* d_sd_pair_data       = d_sp_pair_data + sp_pair_data.size();
    double* d_pp_pair_data       = d_sd_pair_data + sd_pair_data.size();
    double* d_pd_pair_data       = d_pp_pair_data + pp_pair_data.size();
    double* d_dd_pair_data       = d_pd_pair_data + pd_pair_data.size();
    double* d_ss_mat_Q_local     = d_dd_pair_data + dd_pair_data.size();
    double* d_sp_mat_Q_local     = d_ss_mat_Q_local + ss_prim_pair_count_local;
    double* d_sd_mat_Q_local     = d_sp_mat_Q_local + sp_prim_pair_count_local;
    double* d_pp_mat_Q_local     = d_sd_mat_Q_local + sd_prim_pair_count_local;
    double* d_pd_mat_Q_local     = d_pp_mat_Q_local + pp_prim_pair_count_local;
    double* d_dd_mat_Q_local     = d_pd_mat_Q_local + pd_prim_pair_count_local;
    double* d_ss_pair_data_local = d_dd_mat_Q_local + dd_prim_pair_count_local;
    double* d_sp_pair_data_local = d_ss_pair_data_local + ss_pair_data_local.size();
    double* d_sd_pair_data_local = d_sp_pair_data_local + sp_pair_data_local.size();
    double* d_pp_pair_data_local = d_sd_pair_data_local + sd_pair_data_local.size();
    double* d_pd_pair_data_local = d_pp_pair_data_local + pp_pair_data_local.size();
    double* d_dd_pair_data_local = d_pd_pair_data_local + pd_pair_data_local.size();

    uint32_t* d_ss_first_inds        = d_data_J_uint32;
    uint32_t* d_ss_second_inds       = d_ss_first_inds  + ss_prim_pair_count;
    uint32_t* d_sp_first_inds        = d_ss_second_inds + ss_prim_pair_count;
    uint32_t* d_sp_second_inds       = d_sp_first_inds  + sp_prim_pair_count;
    uint32_t* d_sd_first_inds        = d_sp_second_inds + sp_prim_pair_count;
    uint32_t* d_sd_second_inds       = d_sd_first_inds  + sd_prim_pair_count;
    uint32_t* d_pp_first_inds        = d_sd_second_inds + sd_prim_pair_count;
    uint32_t* d_pp_second_inds       = d_pp_first_inds  + pp_prim_pair_count;
    uint32_t* d_pd_first_inds        = d_pp_second_inds + pp_prim_pair_count;
    uint32_t* d_pd_second_inds       = d_pd_first_inds  + pd_prim_pair_count;
    uint32_t* d_dd_first_inds        = d_pd_second_inds + pd_prim_pair_count;
    uint32_t* d_dd_second_inds       = d_dd_first_inds  + dd_prim_pair_count;
    uint32_t* d_ss_first_inds_local  = d_dd_second_inds + dd_prim_pair_count;
    uint32_t* d_ss_second_inds_local = d_ss_first_inds_local  + ss_prim_pair_count_local;
    uint32_t* d_sp_first_inds_local  = d_ss_second_inds_local + ss_prim_pair_count_local;
    uint32_t* d_sp_second_inds_local = d_sp_first_inds_local  + sp_prim_pair_count_local;
    uint32_t* d_sd_first_inds_local  = d_sp_second_inds_local + sp_prim_pair_count_local;
    uint32_t* d_sd_second_inds_local = d_sd_first_inds_local  + sd_prim_pair_count_local;
    uint32_t* d_pp_first_inds_local  = d_sd_second_inds_local + sd_prim_pair_count_local;
    uint32_t* d_pp_second_inds_local = d_pp_first_inds_local  + pp_prim_pair_count_local;
    uint32_t* d_pd_first_inds_local  = d_pp_second_inds_local + pp_prim_pair_count_local;
    uint32_t* d_pd_second_inds_local = d_pd_first_inds_local  + pd_prim_pair_count_local;
    uint32_t* d_dd_first_inds_local  = d_pd_second_inds_local + pd_prim_pair_count_local;
    uint32_t* d_dd_second_inds_local = d_dd_first_inds_local  + dd_prim_pair_count_local;

    gpu::chunkedMemcpyHostToDevice<double>(d_ss_mat_Q, ss_mat_Q.data(), ss_mat_Q.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sp_mat_Q, sp_mat_Q.data(), sp_mat_Q.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sd_mat_Q, sd_mat_Q.data(), sd_mat_Q.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pp_mat_Q, pp_mat_Q.data(), pp_mat_Q.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pd_mat_Q, pd_mat_Q.data(), pd_mat_Q.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_dd_mat_Q, dd_mat_Q.data(), dd_mat_Q.size());

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds,  ss_first_inds.data(),  ss_first_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds, ss_second_inds.data(), ss_second_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds,  sp_first_inds.data(),  sp_first_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds, sp_second_inds.data(), sp_second_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds,  sd_first_inds.data(),  sd_first_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds, sd_second_inds.data(), sd_second_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds,  pp_first_inds.data(),  pp_first_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds, pp_second_inds.data(), pp_second_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds,  pd_first_inds.data(),  pd_first_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds, pd_second_inds.data(), pd_second_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds,  dd_first_inds.data(),  dd_first_inds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds, dd_second_inds.data(), dd_second_inds.size());

    gpu::chunkedMemcpyHostToDevice<double>(d_ss_pair_data, ss_pair_data.data(), ss_pair_data.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sp_pair_data, sp_pair_data.data(), sp_pair_data.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sd_pair_data, sd_pair_data.data(), sd_pair_data.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pp_pair_data, pp_pair_data.data(), pp_pair_data.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pd_pair_data, pd_pair_data.data(), pd_pair_data.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_dd_pair_data, dd_pair_data.data(), dd_pair_data.size());

    gpu::chunkedMemcpyHostToDevice<double>(d_ss_mat_Q_local, ss_mat_Q_local.data(), ss_mat_Q_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sp_mat_Q_local, sp_mat_Q_local.data(), sp_mat_Q_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sd_mat_Q_local, sd_mat_Q_local.data(), sd_mat_Q_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pp_mat_Q_local, pp_mat_Q_local.data(), pp_mat_Q_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pd_mat_Q_local, pd_mat_Q_local.data(), pd_mat_Q_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_dd_mat_Q_local, dd_mat_Q_local.data(), dd_mat_Q_local.size());

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_first_inds_local,  ss_first_inds_local.data(),  ss_first_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_ss_second_inds_local, ss_second_inds_local.data(), ss_second_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_first_inds_local,  sp_first_inds_local.data(),  sp_first_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sp_second_inds_local, sp_second_inds_local.data(), sp_second_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_first_inds_local,  sd_first_inds_local.data(),  sd_first_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_sd_second_inds_local, sd_second_inds_local.data(), sd_second_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_first_inds_local,  pp_first_inds_local.data(),  pp_first_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pp_second_inds_local, pp_second_inds_local.data(), pp_second_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_first_inds_local,  pd_first_inds_local.data(),  pd_first_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pd_second_inds_local, pd_second_inds_local.data(), pd_second_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_first_inds_local,  dd_first_inds_local.data(),  dd_first_inds_local.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_dd_second_inds_local, dd_second_inds_local.data(), dd_second_inds_local.size());

    gpu::chunkedMemcpyHostToDevice<double>(d_ss_pair_data_local, ss_pair_data_local.data(), ss_pair_data_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sp_pair_data_local, sp_pair_data_local.data(), sp_pair_data_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_sd_pair_data_local, sd_pair_data_local.data(), sd_pair_data_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pp_pair_data_local, pp_pair_data_local.data(), pp_pair_data_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pd_pair_data_local, pd_pair_data_local.data(), pd_pair_data_local.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_dd_pair_data_local, dd_pair_data_local.data(), dd_pair_data_local.size());

    mat_Fock_omp[gpu_id].zero();

    gpuSafe(gpuDeviceSynchronize());

    omptimers[thread_id].stop("J prep.");

#pragma omp barrier

    omptimers[thread_id].start("J compute");

    // compute J

    if (std::fabs(prefac_coulomb) > 1.0e-13)
    {

    // J: S-S block

    if (ss_prim_pair_count_local > 0)
    {
        //omptimers[thread_id].start("  J block SS");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J, static_cast<uint32_t>(ss_prim_pair_count_local));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SS|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, ss_mat_D.data(), ss_prim_pair_count);

            //omptimers[thread_id].start("    J block SSSS");

            gpu::computeCoulombFockSSSS<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block SSSS");
        }

        // J: (SS|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sp_mat_D.data(), sp_prim_pair_count);

            gpu::computeCoulombFockSSSP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SS|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sd_mat_D.data(), sd_prim_pair_count);

            gpu::computeCoulombFockSSSD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SS|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pp_mat_D.data(), pp_prim_pair_count);

            gpu::computeCoulombFockSSPP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SS|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pd_mat_D.data(), pd_prim_pair_count);

            gpu::computeCoulombFockSSPD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SS|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, dd_mat_D.data(), dd_prim_pair_count);

            gpu::computeCoulombFockSSDD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        gpu::chunkedMemcpyDeviceToHost<double>(mat_J.data(), d_mat_J, ss_prim_pair_count_local);

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * prefac_coulomb;

            if (i != j) mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * prefac_coulomb;
        }

        //omptimers[thread_id].stop("  J block SS");
    }

    // J: S-P block

    if (sp_prim_pair_count_local > 0)
    {
        //omptimers[thread_id].start("  J block SP");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J, static_cast<uint32_t>(sp_prim_pair_count_local));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, ss_mat_D.data(), ss_prim_pair_count);

            gpu::computeCoulombFockSPSS<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sp_mat_D.data(), sp_prim_pair_count);

            gpu::computeCoulombFockSPSP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sd_mat_D.data(), sd_prim_pair_count);

            gpu::computeCoulombFockSPSD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pp_mat_D.data(), pp_prim_pair_count);

            gpu::computeCoulombFockSPPP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pd_mat_D.data(), pd_prim_pair_count);

            gpu::computeCoulombFockSPPD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, dd_mat_D.data(), dd_prim_pair_count);

            gpu::computeCoulombFockSPDD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        gpu::chunkedMemcpyDeviceToHost<double>(mat_J.data(), d_mat_J, sp_prim_pair_count_local);

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

        //omptimers[thread_id].stop("  J block SP");
    }

    // J: P-P block

    if (pp_prim_pair_count_local > 0)
    {
        //omptimers[thread_id].start("  J block PP");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J, static_cast<uint32_t>(pp_prim_pair_count_local));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, ss_mat_D.data(), ss_prim_pair_count);

            gpu::computeCoulombFockPPSS<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sp_mat_D.data(), sp_prim_pair_count);

            gpu::computeCoulombFockPPSP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sd_mat_D.data(), sd_prim_pair_count);

            gpu::computeCoulombFockPPSD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {

            double *d_mat_J2, *d_mat_J2_2kernels, *d_mat_J2_alld, *d_mat_J2_inside, *d_mat_J2_ref;
            gpuSafe(gpuMalloc(&d_mat_J2, sizeof(double) * pp_prim_pair_count_local));
            gpuSafe(gpuMalloc(&d_mat_J2_2kernels, sizeof(double) * pp_prim_pair_count_local));
            gpuSafe(gpuMalloc(&d_mat_J2_alld, sizeof(double) * pp_prim_pair_count_local));
            gpuSafe(gpuMalloc(&d_mat_J2_ref, sizeof(double) * pp_prim_pair_count_local));
            gpuSafe(gpuMalloc(&d_mat_J2_inside, sizeof(double) * pp_prim_pair_count_local));
            

            double *d_mat_D2;
            gpuMalloc(&d_mat_D2, sizeof(double)*pp_prim_pair_count);
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D2, pp_mat_D.data(), pp_prim_pair_count);

            gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J2, static_cast<uint32_t>(pp_prim_pair_count_local));
            gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J2_2kernels, static_cast<uint32_t>(pp_prim_pair_count_local));
            gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J2_alld, static_cast<uint32_t>(pp_prim_pair_count_local));
            gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J2_ref, static_cast<uint32_t>(pp_prim_pair_count_local));
            gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J2_inside, static_cast<uint32_t>(pp_prim_pair_count_local));


            const double tau_precision = 1e-6;
            int32_t* d_cut_ij_tile = nullptr;

            // 1) CPU 端算 cut_ij_tile（per block/ij-tile）
            auto cut_ij_tile_h = build_cut_ij_tile(
                pp_mat_Q_local,
                pp_mat_Q,
                pp_mat_D,
                (uint32_t)pp_prim_pair_count_local,
                (uint32_t)pp_prim_pair_count,
                TILE_DIM,
                tau_precision);

            // 2) 分配 + 拷到 GPU
            const uint32_t nij_tiles = (pp_prim_pair_count_local + TILE_DIM - 1) / TILE_DIM;

            gpuSafe(gpuMalloc((void**)&d_cut_ij_tile, nij_tiles * sizeof(int32_t)));

            gpu::chunkedMemcpyHostToDevice<int32_t>(d_cut_ij_tile, cut_ij_tile_h.data(), nij_tiles);

            // ---- prepare float versions on host ----
            std::vector<float> p_prim_info_f       = to_float_vec(p_prim_info); 
            std::vector<float> pp_mat_D_f          = to_float_vec(pp_mat_D);
            std::vector<float> pp_mat_Q_f          = to_float_vec(pp_mat_Q);
            std::vector<float> pp_mat_Q_local_f    = to_float_vec(pp_mat_Q_local);
            std::vector<float> pp_pair_data_f      = to_float_vec(pp_pair_data);
            std::vector<float> pp_pair_data_local_f= to_float_vec(pp_pair_data_local);

            float* d_p_prim_info_f        = nullptr;
            float* d_pp_mat_D_f           = nullptr;
            float* d_pp_mat_Q_f           = nullptr;
            float* d_pp_mat_Q_local_f     = nullptr;
            float* d_pp_pair_data_f       = nullptr;
            float* d_pp_pair_data_local_f = nullptr;

            gpuSafe(gpuMalloc(&d_p_prim_info_f,
                    sizeof(float) * p_prim_info_f.size()));
            gpuSafe(gpuMalloc(&d_pp_mat_D_f,
                    sizeof(float) * pp_mat_D_f.size()));
            gpuSafe(gpuMalloc(&d_pp_mat_Q_f,
                    sizeof(float) * pp_mat_Q_f.size()));
            gpuSafe(gpuMalloc(&d_pp_mat_Q_local_f,
                    sizeof(float) * pp_mat_Q_local_f.size()));
            gpuSafe(gpuMalloc(&d_pp_pair_data_f,
                    sizeof(float) * pp_pair_data_f.size()));
            gpuSafe(gpuMalloc(&d_pp_pair_data_local_f,
                    sizeof(float) * pp_pair_data_local_f.size()));

            gpu::chunkedMemcpyHostToDevice<float>(
                d_p_prim_info_f,
                p_prim_info_f.data(),
                p_prim_info_f.size());

            gpu::chunkedMemcpyHostToDevice<float>(
                d_pp_mat_D_f,
                pp_mat_D_f.data(),
                pp_mat_D_f.size());

            gpu::chunkedMemcpyHostToDevice<float>(
                d_pp_mat_Q_f,
                pp_mat_Q_f.data(),
                pp_mat_Q_f.size());

            gpu::chunkedMemcpyHostToDevice<float>(
                d_pp_mat_Q_local_f,
                pp_mat_Q_local_f.data(),
                pp_mat_Q_local_f.size());

            gpu::chunkedMemcpyHostToDevice<float>(
                d_pp_pair_data_f,
                pp_pair_data_f.data(),
                pp_pair_data_f.size());

            gpu::chunkedMemcpyHostToDevice<float>(
                d_pp_pair_data_local_f,
                pp_pair_data_local_f.data(),
                pp_pair_data_local_f.size());
            
            gpu::computeCoulombFockPPPP_MP_alld<<<num_blocks, threads_per_block>>>(
                               d_mat_J2_alld,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D2,
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
                               eri_threshold,
                               d_cut_ij_tile,
                               tau_precision);

            gpu::computeCoulombFockPPPP_ori<<<num_blocks, threads_per_block>>>(
                               d_mat_J2_ref,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D2,
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

            gpu::computeCoulombFockPPPP_MP_inside_cast<<<num_blocks, threads_per_block>>>(
                               d_mat_J2_inside,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D2,
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
                               d_boys_func_table_f,
                               d_boys_func_ft_f,
                               eri_threshold,
                               d_cut_ij_tile,
                               tau_precision);

            gpu::computeCoulombFockPPPP_MP<<<num_blocks, threads_per_block>>>(
                               d_mat_J2,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D2,
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
                               eri_threshold,
                               d_cut_ij_tile,
                               tau_precision,
                               d_p_prim_info_f,
                               d_pp_mat_D_f,
                               d_pp_mat_Q_local_f,
                               d_pp_mat_Q_f,
                               d_pp_pair_data_local_f,
                               d_pp_pair_data_f,
                               d_boys_func_table_f,
                               d_boys_func_ft_f);
            
            gpu::computeCoulombFockPPPP_FP64<<<num_blocks, threads_per_block>>>(
                               d_mat_J2_2kernels,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D2,
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
                               eri_threshold,
                               d_cut_ij_tile,
                               tau_precision,
                               d_p_prim_info_f,
                               d_pp_mat_D_f,
                               d_pp_mat_Q_local_f,
                               d_pp_mat_Q_f,
                               d_pp_pair_data_local_f,
                               d_pp_pair_data_f,
                               d_boys_func_table_f,
                               d_boys_func_ft_f);
            
            gpu::computeCoulombFockPPPP_FP32<<<num_blocks, threads_per_block>>>(
                               d_mat_J2_2kernels,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D2,
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
                               eri_threshold,
                               d_cut_ij_tile,
                               tau_precision,
                               d_p_prim_info_f,
                               d_pp_mat_D_f,
                               d_pp_mat_Q_local_f,
                               d_pp_mat_Q_f,
                               d_pp_pair_data_local_f,
                               d_pp_pair_data_f,
                               d_boys_func_table_f,
                               d_boys_func_ft_f);

            gpuSafe(gpuDeviceSynchronize());

            std::vector<double> h_mat_J2(pp_prim_pair_count_local, 0.0);
            std::vector<double> h_mat_J2_2kernels(pp_prim_pair_count_local, 0.0);
            std::vector<double> h_mat_J2_alld(pp_prim_pair_count_local, 0.0);
            std::vector<double> h_mat_J2_ref(pp_prim_pair_count_local, 0.0);
            std::vector<double> h_mat_J2_inside(pp_prim_pair_count_local, 0.0);

            gpu::chunkedMemcpyDeviceToHost<double>(h_mat_J2.data(), d_mat_J2, pp_prim_pair_count_local);
            gpu::chunkedMemcpyDeviceToHost<double>(h_mat_J2_2kernels.data(), d_mat_J2_2kernels, pp_prim_pair_count_local);
            gpu::chunkedMemcpyDeviceToHost<double>(h_mat_J2_alld.data(), d_mat_J2_alld, pp_prim_pair_count_local);
            gpu::chunkedMemcpyDeviceToHost<double>(h_mat_J2_ref.data(), d_mat_J2_ref, pp_prim_pair_count_local);
            gpu::chunkedMemcpyDeviceToHost<double>(h_mat_J2_inside.data(), d_mat_J2_inside, pp_prim_pair_count_local);

            gpuSafe(gpuFree(d_mat_J2));
            gpuSafe(gpuFree(d_mat_J2_2kernels));
            gpuSafe(gpuFree(d_mat_J2_alld));
            gpuSafe(gpuFree(d_mat_J2_ref));
            gpuSafe(gpuFree(d_mat_J2_inside));
            gpuSafe(gpuFree(d_cut_ij_tile));
            gpuSafe(gpuFree(d_mat_D2));

            gpuSafe(gpuFree(d_p_prim_info_f));
            gpuSafe(gpuFree(d_pp_mat_D_f));
            gpuSafe(gpuFree(d_pp_mat_Q_f));
            gpuSafe(gpuFree(d_pp_mat_Q_local_f));
            gpuSafe(gpuFree(d_pp_pair_data_f));
            gpuSafe(gpuFree(d_pp_pair_data_local_f));

            auto check_J_against_ref =
            [&](const char* tag,
                const std::vector<double> &J,
                const std::vector<double> &J_ref,
                uint32_t n)
            {
                double max_abs_err = 0.0;
                double rms_err = 0.0;
                double max_ref = 0.0;
                int    idx_max = -1;
                uint32_t nan_count = 0;
                uint32_t nan_idx = 0;

                for (uint32_t i = 0; i < n; ++i)
                {
                    const double ref  = J_ref[i];
                    const double val  = J[i];

                    if (!std::isfinite(ref) || !std::isfinite(val))
                    {
                        if (nan_count == 0) nan_idx = i;
                        nan_count++;
                        continue;
                    }

                    const double diff = std::abs(val - ref);

                    if (diff > max_abs_err)
                    {
                        max_abs_err = diff;
                        idx_max = i;
                    }

                    rms_err += diff * diff;
                    max_ref = std::max(max_ref, std::abs(ref));
                }

                rms_err = std::sqrt(rms_err / n);
                const double rel_err = (max_ref > 0.0) ? (max_abs_err / max_ref) : 0.0;

                std::cout << std::scientific;
                std::cout << "=== " << tag << " ===\n";
                std::cout << "  max |ΔJ|   = " << max_abs_err << "\n";
                std::cout << "  rms |ΔJ|   = " << rms_err << "\n";
                std::cout << "  max |Jref| = " << max_ref << "\n";
                std::cout << "  rel error  = " << rel_err << "\n";

                if (idx_max >= 0)
                {
                    std::cout << "  worst idx  = " << idx_max
                            << "  J = " << J[idx_max]
                            << "  J_ref = " << J_ref[idx_max]
                            << "\n";
                }

                if (nan_count > 0)
                {
                    std::cout << "  NaN entries = " << nan_count
                            << " (first idx = " << nan_idx << ")\n";
                }

                std::cout << "============================\n";
            };

            check_J_against_ref("PPPP_MP ALL-D kernel check (J2_alld vs ref)",
                                h_mat_J2_alld,
                                h_mat_J2_ref,
                                (uint32_t)pp_prim_pair_count_local);
            check_J_against_ref("Two seperate kernels check with original (J2_2kernels vs ref)",
                                h_mat_J2_2kernels,
                                h_mat_J2_ref,
                                (uint32_t)pp_prim_pair_count_local);       
            check_J_against_ref("Two seperate kernels check with one MP kernel (PPPP_MP)  (J2_2kernels vs J2)",
                                h_mat_J2_2kernels,
                                h_mat_J2,
                                (uint32_t)pp_prim_pair_count_local);            
            check_J_against_ref("PPPP_MP kernel check (J2 vs ref)",
                    h_mat_J2,
                    h_mat_J2_ref,
                    (uint32_t)pp_prim_pair_count_local);
            check_J_against_ref("PPPP_MP kernel check (J2_inside vs ref)",
                    h_mat_J2_inside,
                    h_mat_J2_ref,
                    (uint32_t)pp_prim_pair_count_local);
            check_J_against_ref("PPPP_MP kernel check (J2_inside vs J2)",
                    h_mat_J2_inside,
                    h_mat_J2,
                    (uint32_t)pp_prim_pair_count_local);

            // cut status
            {
                const int32_t m_tiles = (int32_t)((pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM);

                int32_t cut_min = INT_MAX;
                int32_t cut_max = 0;
                long long cut_sum = 0;

                // 可选：看有多少 tile cut=0 / cut=m_tiles
                int32_t n_cut0 = 0;
                int32_t n_cutfull = 0;

                for (int32_t c : cut_ij_tile_h) {
                    cut_min = std::min(cut_min, c);
                    cut_max = std::max(cut_max, c);
                    cut_sum += c;
                    if (c == 0) n_cut0++;
                    if (c >= m_tiles) n_cutfull++;
                }

                const double cut_avg = (double)cut_sum / (double)cut_ij_tile_h.size();
                const double fp32_frac = (m_tiles > 0) ? (double)(m_tiles - cut_avg) / (double)m_tiles : 0.0;

                printf("=== PPPP MP cut stats (host) ===\n");
                printf("  ij_tiles         = %zu\n", cut_ij_tile_h.size());
                printf("  kl_tiles (m_tiles)= %d\n", m_tiles);
                printf("  cut min / max    = %d / %d\n", cut_min, cut_max);
                printf("  cut avg          = %.2f\n", cut_avg);
                printf("  FP32 fraction    = %.1f %%\n", fp32_frac * 100.0);
                printf("  cut==0 tiles     = %d\n", n_cut0);
                printf("  cut>=m_tiles     = %d\n", n_cutfull);
                printf("===============================\n");
            }

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pp_mat_D.data(), pp_prim_pair_count);

            gpu::computeCoulombFockPPPP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pd_mat_D.data(), pd_prim_pair_count);

            gpu::computeCoulombFockPPPD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block PPDD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, dd_mat_D.data(), dd_prim_pair_count);

            gpu::computeCoulombFockPPDD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block PPDD");
        }

        gpu::chunkedMemcpyDeviceToHost<double>(mat_J.data(), d_mat_J, pp_prim_pair_count_local);

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

        //omptimers[thread_id].stop("  J block PP");
    }

    // J: S-D block

    if (sd_prim_pair_count_local > 0)
    {
        //omptimers[thread_id].start("  J block SD");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J, static_cast<uint32_t>(sd_prim_pair_count_local));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, ss_mat_D.data(), ss_prim_pair_count);

            gpu::computeCoulombFockSDSS<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sp_mat_D.data(), sp_prim_pair_count);

            gpu::computeCoulombFockSDSP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sd_mat_D.data(), sd_prim_pair_count);

            gpu::computeCoulombFockSDSD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pp_mat_D.data(), pp_prim_pair_count);

            gpu::computeCoulombFockSDPP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pd_mat_D.data(), pd_prim_pair_count);

            gpu::computeCoulombFockSDPD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (SD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block SDDD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, dd_mat_D.data(), dd_prim_pair_count);

            gpu::computeCoulombFockSDDD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block SDDD");
        }

        gpu::chunkedMemcpyDeviceToHost<double>(mat_J.data(), d_mat_J, sd_prim_pair_count_local);

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

        //omptimers[thread_id].stop("  J block SD");
    }

    // J: P-D block

    if (pd_prim_pair_count_local > 0)
    {
        //omptimers[thread_id].start("  J block PD");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J, static_cast<uint32_t>(pd_prim_pair_count_local));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, ss_mat_D.data(), ss_prim_pair_count);

            gpu::computeCoulombFockPDSS<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sp_mat_D.data(), sp_prim_pair_count);

            gpu::computeCoulombFockPDSP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());
        }

        // J: (PD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block PDSD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sd_mat_D.data(), sd_prim_pair_count);

            gpu::computeCoulombFockPDSD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block PDSD");
        }

        // J: (PD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block PDPP");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pp_mat_D.data(), pp_prim_pair_count);

            gpu::computeCoulombFockPDPP<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block PDPP");
        }

        // J: (PD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block PDPD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pd_mat_D.data(), pd_prim_pair_count);

            gpu::computeCoulombFockPDPD<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block PDPD");
        }

        // J: (PD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block PDDD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, dd_mat_D.data(), dd_prim_pair_count);

            gpu::computeCoulombFockPDDD0<<<num_blocks, threads_per_block>>>(
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

            gpu::computeCoulombFockPDDD1<<<num_blocks, threads_per_block>>>(
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

            gpu::computeCoulombFockPDDD2<<<num_blocks, threads_per_block>>>(
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

            gpu::computeCoulombFockPDDD3<<<num_blocks, threads_per_block>>>(
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

            gpu::computeCoulombFockPDDD4<<<num_blocks, threads_per_block>>>(
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

            gpu::computeCoulombFockPDDD5<<<num_blocks, threads_per_block>>>(
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

            gpu::computeCoulombFockPDDD6<<<num_blocks, threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block PDDD");
        }

        gpu::chunkedMemcpyDeviceToHost<double>(mat_J.data(), d_mat_J, pd_prim_pair_count_local);

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

        //omptimers[thread_id].stop("  J block PD");
    }

    // J: D-D block

    if (dd_prim_pair_count_local > 0)
    {
        //omptimers[thread_id].start("  J block DD");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_J, static_cast<uint32_t>(dd_prim_pair_count_local));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (DD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block DDSS");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, ss_mat_D.data(), ss_prim_pair_count);

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDSS<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block DDSS");
        }

        // J: (DD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block DDSP");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sp_mat_D.data(), sp_prim_pair_count);

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDSP<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block DDSP");
        }

        // J: (DD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block DDSD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, sd_mat_D.data(), sd_prim_pair_count);

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDSD<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block DDSD");
        }

        // J: (DD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block DDPP");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pp_mat_D.data(), pp_prim_pair_count);

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDPP<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block DDPP");
        }

        // J: (DD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block DDPD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, pd_mat_D.data(), pd_prim_pair_count);

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDPD0<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD1<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD2<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD3<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD4<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD5<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD6<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD7<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD8<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDPD9<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block DDPD");
        }

        // J: (DD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            //omptimers[thread_id].start("    J block DDDD");

            gpu::chunkedMemcpyHostToDevice<double>(d_mat_D, dd_mat_D.data(), dd_prim_pair_count);

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDDD0<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD1<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD2<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD3<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD4<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD5<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD6<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD7<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD8<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD9<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD10<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD11<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD12<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD13<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD14<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD15<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD16<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD17<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD18<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD19<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD20<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD21<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD22<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD23<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD24<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD25<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD26<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD27<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD28<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpu::computeCoulombFockDDDD29<<<dd_num_blocks, dd_threads_per_block>>>(
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

            gpuSafe(gpuDeviceSynchronize());

            //omptimers[thread_id].stop("    J block DDDD");
        }

        gpu::chunkedMemcpyDeviceToHost<double>(mat_J.data(), d_mat_J, dd_prim_pair_count_local);

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

        //omptimers[thread_id].stop("  J block DD");
    }

    }  // end of compute J

    gpuSafe(gpuDeviceSynchronize());

    omptimers[thread_id].stop("J compute");

#pragma omp barrier

    omptimers[thread_id].start("K prep.");

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

    const auto max_pair_inds_count = std::max({pair_inds_count_for_K_ss, pair_inds_count_for_K_sp, pair_inds_count_for_K_pp,
                                               pair_inds_count_for_K_sd, pair_inds_count_for_K_pd, pair_inds_count_for_K_dd});

    std::vector<double> mat_K(max_pair_inds_count);

    size_t data_K_double_size = 0;

    data_K_double_size += static_cast<size_t>(max_pair_inds_count);         // d_mat_K
    data_K_double_size += static_cast<size_t>(cart_naos * cart_naos);       // d_mat_D_full_AO
    data_K_double_size += Q_K_ss.size();                                    // d_Q_K_ss
    data_K_double_size += Q_K_sp.size();                                    // d_Q_K_sp
    data_K_double_size += Q_K_ps.size();                                    // d_Q_K_ps
    data_K_double_size += Q_K_sd.size();                                    // d_Q_K_sd
    data_K_double_size += Q_K_ds.size();                                    // d_Q_K_ds
    data_K_double_size += Q_K_pp.size();                                    // d_Q_K_pp
    data_K_double_size += Q_K_pd.size();                                    // d_Q_K_pd
    data_K_double_size += Q_K_dp.size();                                    // d_Q_K_dp
    data_K_double_size += Q_K_dd.size();                                    // d_Q_K_dd
    data_K_double_size += pair_data_K_ss.size();                            // d_pair_data_K_ss
    data_K_double_size += pair_data_K_sp.size();                            // d_pair_data_K_sp
    data_K_double_size += pair_data_K_ps.size();                            // d_pair_data_K_ps
    data_K_double_size += pair_data_K_sd.size();                            // d_pair_data_K_sd
    data_K_double_size += pair_data_K_ds.size();                            // d_pair_data_K_ds
    data_K_double_size += pair_data_K_pp.size();                            // d_pair_data_K_pp
    data_K_double_size += pair_data_K_pd.size();                            // d_pair_data_K_pd
    data_K_double_size += pair_data_K_dp.size();                            // d_pair_data_K_dp
    data_K_double_size += pair_data_K_dd.size();                            // d_pair_data_K_dd

    size_t data_K_uint32_size = 0;

    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_ss);    // d_pair_inds_i_for_K_ss
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_ss);    // d_pair_inds_k_for_K_ss
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_sp);    // d_pair_inds_i_for_K_sp
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_sp);    // d_pair_inds_k_for_K_sp
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_sd);    // d_pair_inds_i_for_K_sd
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_sd);    // d_pair_inds_k_for_K_sd
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_pp);    // d_pair_inds_i_for_K_pp
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_pp);    // d_pair_inds_k_for_K_pp
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_pd);    // d_pair_inds_i_for_K_pd
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_pd);    // d_pair_inds_k_for_K_pd
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_dd);    // d_pair_inds_i_for_K_dd
    data_K_uint32_size += static_cast<size_t>(pair_inds_count_for_K_dd);    // d_pair_inds_k_for_K_dd
    data_K_uint32_size += D_inds_K_ss.size();                               // d_D_inds_K_ss
    data_K_uint32_size += D_inds_K_sp.size();                               // d_D_inds_K_sp
    data_K_uint32_size += D_inds_K_ps.size();                               // d_D_inds_K_ps
    data_K_uint32_size += D_inds_K_sd.size();                               // d_D_inds_K_sd
    data_K_uint32_size += D_inds_K_ds.size();                               // d_D_inds_K_ds
    data_K_uint32_size += D_inds_K_pp.size();                               // d_D_inds_K_pp
    data_K_uint32_size += D_inds_K_pd.size();                               // d_D_inds_K_pd
    data_K_uint32_size += D_inds_K_dp.size();                               // d_D_inds_K_dp
    data_K_uint32_size += D_inds_K_dd.size();                               // d_D_inds_K_dd
    data_K_uint32_size += pair_displs_K_ss.size();                          // d_pair_displs_K_ss
    data_K_uint32_size += pair_displs_K_sp.size();                          // d_pair_displs_K_sp
    data_K_uint32_size += pair_displs_K_ps.size();                          // d_pair_displs_K_ps
    data_K_uint32_size += pair_displs_K_sd.size();                          // d_pair_displs_K_sd
    data_K_uint32_size += pair_displs_K_ds.size();                          // d_pair_displs_K_ds
    data_K_uint32_size += pair_displs_K_pp.size();                          // d_pair_displs_K_pp
    data_K_uint32_size += pair_displs_K_pd.size();                          // d_pair_displs_K_pd
    data_K_uint32_size += pair_displs_K_dp.size();                          // d_pair_displs_K_dp
    data_K_uint32_size += pair_displs_K_dd.size();                          // d_pair_displs_K_dd
    data_K_uint32_size += pair_counts_K_ss.size();                          // d_pair_counts_K_ss
    data_K_uint32_size += pair_counts_K_sp.size();                          // d_pair_counts_K_sp
    data_K_uint32_size += pair_counts_K_ps.size();                          // d_pair_counts_K_ps
    data_K_uint32_size += pair_counts_K_sd.size();                          // d_pair_counts_K_sd
    data_K_uint32_size += pair_counts_K_ds.size();                          // d_pair_counts_K_ds
    data_K_uint32_size += pair_counts_K_pp.size();                          // d_pair_counts_K_pp
    data_K_uint32_size += pair_counts_K_pd.size();                          // d_pair_counts_K_pd
    data_K_uint32_size += pair_counts_K_dp.size();                          // d_pair_counts_K_dp
    data_K_uint32_size += pair_counts_K_dd.size();                          // d_pair_counts_K_dd

    // K data on device

#pragma omp barrier

    screening.resize_devptr_double(gpu_id, data_K_double_size);
    screening.resize_devptr_uint32(gpu_id, data_K_uint32_size);

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    auto d_data_K_double = screening.get_devptr_double(gpu_id);
    auto d_data_K_uint32 = screening.get_devptr_uint32(gpu_id);

    double* d_mat_K          = d_data_K_double;
    double* d_mat_D_full_AO  = d_mat_K          + max_pair_inds_count;
    double* d_Q_K_ss         = d_mat_D_full_AO  + cart_naos * cart_naos;
    double* d_Q_K_sp         = d_Q_K_ss         + Q_K_ss.size();
    double* d_Q_K_ps         = d_Q_K_sp         + Q_K_sp.size();
    double* d_Q_K_sd         = d_Q_K_ps         + Q_K_ps.size();
    double* d_Q_K_ds         = d_Q_K_sd         + Q_K_sd.size();
    double* d_Q_K_pp         = d_Q_K_ds         + Q_K_ds.size();
    double* d_Q_K_pd         = d_Q_K_pp         + Q_K_pp.size();
    double* d_Q_K_dp         = d_Q_K_pd         + Q_K_pd.size();
    double* d_Q_K_dd         = d_Q_K_dp         + Q_K_dp.size();
    double* d_pair_data_K_ss = d_Q_K_dd         + Q_K_dd.size();
    double* d_pair_data_K_sp = d_pair_data_K_ss + pair_data_K_ss.size();
    double* d_pair_data_K_ps = d_pair_data_K_sp + pair_data_K_sp.size();
    double* d_pair_data_K_sd = d_pair_data_K_ps + pair_data_K_ps.size();
    double* d_pair_data_K_ds = d_pair_data_K_sd + pair_data_K_sd.size();
    double* d_pair_data_K_pp = d_pair_data_K_ds + pair_data_K_ds.size();
    double* d_pair_data_K_pd = d_pair_data_K_pp + pair_data_K_pp.size();
    double* d_pair_data_K_dp = d_pair_data_K_pd + pair_data_K_pd.size();
    double* d_pair_data_K_dd = d_pair_data_K_dp + pair_data_K_dp.size();

    uint32_t* d_pair_inds_i_for_K_ss = d_data_K_uint32;
    uint32_t* d_pair_inds_k_for_K_ss = d_pair_inds_i_for_K_ss + pair_inds_count_for_K_ss;
    uint32_t* d_pair_inds_i_for_K_sp = d_pair_inds_k_for_K_ss + pair_inds_count_for_K_ss;
    uint32_t* d_pair_inds_k_for_K_sp = d_pair_inds_i_for_K_sp + pair_inds_count_for_K_sp;
    uint32_t* d_pair_inds_i_for_K_sd = d_pair_inds_k_for_K_sp + pair_inds_count_for_K_sp;
    uint32_t* d_pair_inds_k_for_K_sd = d_pair_inds_i_for_K_sd + pair_inds_count_for_K_sd;
    uint32_t* d_pair_inds_i_for_K_pp = d_pair_inds_k_for_K_sd + pair_inds_count_for_K_sd;
    uint32_t* d_pair_inds_k_for_K_pp = d_pair_inds_i_for_K_pp + pair_inds_count_for_K_pp;
    uint32_t* d_pair_inds_i_for_K_pd = d_pair_inds_k_for_K_pp + pair_inds_count_for_K_pp;
    uint32_t* d_pair_inds_k_for_K_pd = d_pair_inds_i_for_K_pd + pair_inds_count_for_K_pd;
    uint32_t* d_pair_inds_i_for_K_dd = d_pair_inds_k_for_K_pd + pair_inds_count_for_K_pd;
    uint32_t* d_pair_inds_k_for_K_dd = d_pair_inds_i_for_K_dd + pair_inds_count_for_K_dd;
    uint32_t* d_D_inds_K_ss          = d_pair_inds_k_for_K_dd + pair_inds_count_for_K_dd;
    uint32_t* d_D_inds_K_sp          = d_D_inds_K_ss          + D_inds_K_ss.size();
    uint32_t* d_D_inds_K_ps          = d_D_inds_K_sp          + D_inds_K_sp.size();
    uint32_t* d_D_inds_K_sd          = d_D_inds_K_ps          + D_inds_K_ps.size();
    uint32_t* d_D_inds_K_ds          = d_D_inds_K_sd          + D_inds_K_sd.size();
    uint32_t* d_D_inds_K_pp          = d_D_inds_K_ds          + D_inds_K_ds.size();
    uint32_t* d_D_inds_K_pd          = d_D_inds_K_pp          + D_inds_K_pp.size();
    uint32_t* d_D_inds_K_dp          = d_D_inds_K_pd          + D_inds_K_pd.size();
    uint32_t* d_D_inds_K_dd          = d_D_inds_K_dp          + D_inds_K_dp.size();
    uint32_t* d_pair_counts_K_ss     = d_D_inds_K_dd          + D_inds_K_dd.size();
    uint32_t* d_pair_counts_K_sp     = d_pair_counts_K_ss     + pair_counts_K_ss.size();
    uint32_t* d_pair_counts_K_ps     = d_pair_counts_K_sp     + pair_counts_K_sp.size();
    uint32_t* d_pair_counts_K_sd     = d_pair_counts_K_ps     + pair_counts_K_ps.size();
    uint32_t* d_pair_counts_K_ds     = d_pair_counts_K_sd     + pair_counts_K_sd.size();
    uint32_t* d_pair_counts_K_pp     = d_pair_counts_K_ds     + pair_counts_K_ds.size();
    uint32_t* d_pair_counts_K_pd     = d_pair_counts_K_pp     + pair_counts_K_pp.size();
    uint32_t* d_pair_counts_K_dp     = d_pair_counts_K_pd     + pair_counts_K_pd.size();
    uint32_t* d_pair_counts_K_dd     = d_pair_counts_K_dp     + pair_counts_K_dp.size();
    uint32_t* d_pair_displs_K_ss     = d_pair_counts_K_dd     + pair_counts_K_dd.size();
    uint32_t* d_pair_displs_K_sp     = d_pair_displs_K_ss     + pair_displs_K_ss.size();
    uint32_t* d_pair_displs_K_ps     = d_pair_displs_K_sp     + pair_displs_K_sp.size();
    uint32_t* d_pair_displs_K_sd     = d_pair_displs_K_ps     + pair_displs_K_ps.size();
    uint32_t* d_pair_displs_K_ds     = d_pair_displs_K_sd     + pair_displs_K_sd.size();
    uint32_t* d_pair_displs_K_pp     = d_pair_displs_K_ds     + pair_displs_K_ds.size();
    uint32_t* d_pair_displs_K_pd     = d_pair_displs_K_pp     + pair_displs_K_pp.size();
    uint32_t* d_pair_displs_K_dp     = d_pair_displs_K_pd     + pair_displs_K_pd.size();
    uint32_t* d_pair_displs_K_dd     = d_pair_displs_K_dp     + pair_displs_K_dp.size();

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_i_for_K_ss, pair_inds_i_for_K_ss.data(), pair_inds_i_for_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_k_for_K_ss, pair_inds_k_for_K_ss.data(), pair_inds_k_for_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_i_for_K_sp, pair_inds_i_for_K_sp.data(), pair_inds_i_for_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_k_for_K_sp, pair_inds_k_for_K_sp.data(), pair_inds_k_for_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_i_for_K_sd, pair_inds_i_for_K_sd.data(), pair_inds_i_for_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_k_for_K_sd, pair_inds_k_for_K_sd.data(), pair_inds_k_for_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_i_for_K_pp, pair_inds_i_for_K_pp.data(), pair_inds_i_for_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_k_for_K_pp, pair_inds_k_for_K_pp.data(), pair_inds_k_for_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_i_for_K_pd, pair_inds_i_for_K_pd.data(), pair_inds_i_for_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_k_for_K_pd, pair_inds_k_for_K_pd.data(), pair_inds_k_for_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_i_for_K_dd, pair_inds_i_for_K_dd.data(), pair_inds_i_for_K_dd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_inds_k_for_K_dd, pair_inds_k_for_K_dd.data(), pair_inds_k_for_K_dd.size());

    gpu::chunkedMemcpyHostToDevice<double>(d_mat_D_full_AO, cart_dens_ptr, cart_naos * cart_naos);

    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_ss, Q_K_ss.data(), Q_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_sp, Q_K_sp.data(), Q_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_ps, Q_K_ps.data(), Q_K_ps.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_sd, Q_K_sd.data(), Q_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_ds, Q_K_ds.data(), Q_K_ds.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_pp, Q_K_pp.data(), Q_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_pd, Q_K_pd.data(), Q_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_dp, Q_K_dp.data(), Q_K_dp.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_Q_K_dd, Q_K_dd.data(), Q_K_dd.size());

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_ss, D_inds_K_ss.data(), D_inds_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_sp, D_inds_K_sp.data(), D_inds_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_ps, D_inds_K_ps.data(), D_inds_K_ps.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_sd, D_inds_K_sd.data(), D_inds_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_ds, D_inds_K_ds.data(), D_inds_K_ds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_pp, D_inds_K_pp.data(), D_inds_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_pd, D_inds_K_pd.data(), D_inds_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_dp, D_inds_K_dp.data(), D_inds_K_dp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_D_inds_K_dd, D_inds_K_dd.data(), D_inds_K_dd.size());

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_ss, pair_displs_K_ss.data(), pair_displs_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_sp, pair_displs_K_sp.data(), pair_displs_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_ps, pair_displs_K_ps.data(), pair_displs_K_ps.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_sd, pair_displs_K_sd.data(), pair_displs_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_ds, pair_displs_K_ds.data(), pair_displs_K_ds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_pp, pair_displs_K_pp.data(), pair_displs_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_pd, pair_displs_K_pd.data(), pair_displs_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_dp, pair_displs_K_dp.data(), pair_displs_K_dp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_displs_K_dd, pair_displs_K_dd.data(), pair_displs_K_dd.size());

    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_ss, pair_counts_K_ss.data(), pair_counts_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_sp, pair_counts_K_sp.data(), pair_counts_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_ps, pair_counts_K_ps.data(), pair_counts_K_ps.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_sd, pair_counts_K_sd.data(), pair_counts_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_ds, pair_counts_K_ds.data(), pair_counts_K_ds.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_pp, pair_counts_K_pp.data(), pair_counts_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_pd, pair_counts_K_pd.data(), pair_counts_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_dp, pair_counts_K_dp.data(), pair_counts_K_dp.size());
    gpu::chunkedMemcpyHostToDevice<uint32_t>(d_pair_counts_K_dd, pair_counts_K_dd.data(), pair_counts_K_dd.size());

    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_ss, pair_data_K_ss.data(), pair_data_K_ss.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_sp, pair_data_K_sp.data(), pair_data_K_sp.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_ps, pair_data_K_ps.data(), pair_data_K_ps.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_sd, pair_data_K_sd.data(), pair_data_K_sd.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_ds, pair_data_K_ds.data(), pair_data_K_ds.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_pp, pair_data_K_pp.data(), pair_data_K_pp.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_pd, pair_data_K_pd.data(), pair_data_K_pd.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_dp, pair_data_K_dp.data(), pair_data_K_dp.size());
    gpu::chunkedMemcpyHostToDevice<double>(d_pair_data_K_dd, pair_data_K_dd.data(), pair_data_K_dd.size());

    gpuSafe(gpuDeviceSynchronize());

    omptimers[thread_id].stop("K prep.");

#pragma omp barrier

    omptimers[thread_id].start("K compute");

    // compute K

    for (int64_t exch_idx = 0; exch_idx < static_cast<int64_t>(frac_exact_exchange_values.size()); exch_idx++)
    {
    const auto frac_exact_exchange = frac_exact_exchange_values[exch_idx];
    const auto omega = omega_values[exch_idx];

    if (std::fabs(frac_exact_exchange) > 1.0e-13)
    {

    // K: S-S block

    if (pair_inds_count_for_K_ss > 0)
    {
        //omptimers[thread_id].start("  K block SS");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_ss + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_ss));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_ss, 1);

        // K: (SS|SS)
        //     *  *

        gpu::computeExchangeFockSSSS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SS|SP)
        //     *  *

        gpu::computeExchangeFockSSSP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|SS)
        //     *  *

        gpu::computeExchangeFockSPSS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|SP)
        //     *  *

        gpu::computeExchangeFockSPSP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SS|SD)
        //     *  *

        gpu::computeExchangeFockSSSD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SD|SS)
        //     *  *

        gpu::computeExchangeFockSDSS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|SD)
        //     *  *

        gpu::computeExchangeFockSPSD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SD|SP)
        //     *  *

        gpu::computeExchangeFockSDSP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SD|SD)
        //     *  *

        gpu::computeExchangeFockSDSD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss);

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

        //omptimers[thread_id].stop("  K block SS");
    }

    // K: S-P block

    if (pair_inds_count_for_K_sp > 0)
    {
        //omptimers[thread_id].start("  K block SP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sp + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_sp));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_sp, 1);

        // K: (SS|PS)
        //     *  *

        gpu::computeExchangeFockSSPS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SS|PP)
        //     *  *

        gpu::computeExchangeFockSSPP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|PS)
        //     *  *

        gpu::computeExchangeFockSPPS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|PP)
        //     *  *

        gpu::computeExchangeFockSPPP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SS|PD)
        //     *  *

        gpu::computeExchangeFockSSPD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SD|PS)
        //     *  *

        gpu::computeExchangeFockSDPS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|PD)
        //     *  *

        gpu::computeExchangeFockSPPD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SD|PP)
        //     *  *

        gpu::computeExchangeFockSDPP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SD|PD)
        //     *  *

        gpu::computeExchangeFockSDPD<<<num_blocks, threads_per_block>>>(
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

        gpu::chunkedMemcpyDeviceToHost<double>(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp);

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

        //omptimers[thread_id].stop("  K block SP");
    }

    // K: P-P block

    if (pair_inds_count_for_K_pp > 0)
    {
        //omptimers[thread_id].start("  K block PP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pp + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_pp));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_pp, 1);

        // K: (PS|PS)
        //     *  *

        //omptimers[thread_id].start("    K block PSPS");

        gpu::computeExchangeFockPSPS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PSPS");

        // K: (PS|PP)
        //     *  *

        //omptimers[thread_id].start("    K block PSPP");

        gpu::computeExchangeFockPSPP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PSPP");

        // K: (PP|PS)
        //     *  *

        //omptimers[thread_id].start("    K block PPPS");

        gpu::computeExchangeFockPPPS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PPPS");

        // K: (PP|PP)
        //     *  *

        //omptimers[thread_id].start("    K block PPPP");

        gpu::computeExchangeFockPPPP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PPPP");

        // K: (PS|PD)
        //     *  *

        //omptimers[thread_id].start("    K block PSPD");

        gpu::computeExchangeFockPSPD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PSPD");

        // K: (PD|PS)
        //     *  *

        //omptimers[thread_id].start("    K block PDPS");

        gpu::computeExchangeFockPDPS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PDPS");

        // K: (PP|PD)
        //     *  *

        //omptimers[thread_id].start("    K block PPPD");

        gpu::computeExchangeFockPPPD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PPPD");

        // K: (PD|PP)
        //     *  *

        //omptimers[thread_id].start("    K block PDPP");

        gpu::computeExchangeFockPDPP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PDPP");

        // K: (PD|PD)
        //     *  *

        //omptimers[thread_id].start("    K block PDPD");

        gpu::computeExchangeFockPDPD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PDPD");

        gpu::chunkedMemcpyDeviceToHost<double>(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp);

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

        //omptimers[thread_id].stop("  K block PP");
    }

    // K: S-D block

    if (pair_inds_count_for_K_sd > 0)
    {
        //omptimers[thread_id].start("  K block SD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sd + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_sd));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_sd, 1);

        // K: (SS|DS)
        //     *  *

        gpu::computeExchangeFockSSDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SS|DP)
        //     *  *

        gpu::computeExchangeFockSSDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|DS)
        //     *  *

        gpu::computeExchangeFockSPDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (SP|DP)
        //     *  *

        //omptimers[thread_id].start("    K block SPDP");

        gpu::computeExchangeFockSPDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block SPDP");

        // K: (SS|DD)
        //     *  *

        //omptimers[thread_id].start("    K block SSDD");

        gpu::computeExchangeFockSSDD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block SSDD");

        // K: (SD|DS)
        //     *  *

        //omptimers[thread_id].start("    K block SDDS");

        gpu::computeExchangeFockSDDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block SDDS");

        // K: (SP|DD)
        //     *  *

        //omptimers[thread_id].start("    K block SPDD");

        gpu::computeExchangeFockSPDD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block SPDD");

        // K: (SD|DP)
        //     *  *

        //omptimers[thread_id].start("    K block SDDP");

        gpu::computeExchangeFockSDDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block SDDP");

        // K: (SD|DD)
        //     *  *

        //omptimers[thread_id].start("    K block SDDD");

        gpu::computeExchangeFockSDDD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block SDDD");

        gpu::chunkedMemcpyDeviceToHost<double>(mat_K.data(), d_mat_K, pair_inds_count_for_K_sd);

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

        //omptimers[thread_id].stop("  K block SD");
    }

    // K: P-D block

    if (pair_inds_count_for_K_pd > 0)
    {
        //omptimers[thread_id].start("  K block PD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pd + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_pd));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_pd, 1);

        // K: (PS|DS)
        //     *  *

        //omptimers[thread_id].start("    K block PSDS");

        gpu::computeExchangeFockPSDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PSDS");

        // K: (PS|DP)
        //     *  *

        //omptimers[thread_id].start("    K block PSDP");

        gpu::computeExchangeFockPSDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PSDP");

        // K: (PP|DS)
        //     *  *

        //omptimers[thread_id].start("    K block PPDS");

        gpu::computeExchangeFockPPDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PPDS");

        // K: (PS|DD)
        //     *  *

        //omptimers[thread_id].start("    K block PSDD");

        gpu::computeExchangeFockPSDD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PSDD");

        // K: (PD|DS)
        //     *  *

        //omptimers[thread_id].start("    K block PDDS");

        gpu::computeExchangeFockPDDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PDDS");

        // K: (PP|DP)
        //     *  *

        //omptimers[thread_id].start("    K block PPDP");

        gpu::computeExchangeFockPPDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PPDP");

        // K: (PP|DD)
        //     *  *

        //omptimers[thread_id].start("    K block PPDD");

        gpu::computeExchangeFockPPDD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PPDD");

        // K: (PD|DP)
        //     *  *

        //omptimers[thread_id].start("    K block PDDP");

        gpu::computeExchangeFockPDDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PDDP");

        // K: (PD|DD)
        //     *  *

        //omptimers[thread_id].start("    K block PDDD");

        gpu::computeExchangeFockPDDD0<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD1<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD2<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD3<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD4<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD5<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD6<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockPDDD7<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block PDDD");

        gpu::chunkedMemcpyDeviceToHost<double>(mat_K.data(), d_mat_K, pair_inds_count_for_K_pd);

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

        //omptimers[thread_id].stop("  K block PD");
    }

    // K: D-D block

    if (pair_inds_count_for_K_dd > 0)
    {
        //omptimers[thread_id].start("  K block DD");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_dd + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block>>>(d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_dd));

        gpuSafe(gpuDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        num_blocks = dim3(pair_inds_count_for_K_dd, 1);

        // K: (DS|DS)
        //     *  *

        gpu::computeExchangeFockDSDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (DS|DP)
        //     *  *

        gpu::computeExchangeFockDSDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (DP|DS)
        //     *  *

        gpu::computeExchangeFockDPDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        // K: (DS|DD)
        //     *  *

        //omptimers[thread_id].start("    K block DSDD");

        gpu::computeExchangeFockDSDD<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block DSDD");

        // K: (DD|DS)
        //     *  *

        //omptimers[thread_id].start("    K block DDDS");

        gpu::computeExchangeFockDDDS<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block DDDS");

        // K: (DP|DP)
        //     *  *

        //omptimers[thread_id].start("    K block DPDP");

        gpu::computeExchangeFockDPDP<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block DPDP");

        // K: (DP|DD)
        //     *  *

        //omptimers[thread_id].start("    K block DPDD");

        gpu::computeExchangeFockDPDD0<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDPDD1<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDPDD2<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDPDD3<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDPDD4<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDPDD5<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDPDD6<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block DPDD");

        // K: (DD|DP)
        //     *  *

        //omptimers[thread_id].start("    K block DDDP");

        gpu::computeExchangeFockDDDP0<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDP1<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDP2<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDP3<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDP4<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDP5<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDP6<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block DDDP");

        // K: (DD|DD)
        //     *  *

        //omptimers[thread_id].start("    K block DDDD");

        gpu::computeExchangeFockDDDD0<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD1<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD2<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD3<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD4<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD5<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD6<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD7<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD8<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD9<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD10<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD11<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD12<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD13<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD14<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD15<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD16<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD17<<<num_blocks, threads_per_block>>>(
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

        gpu::computeExchangeFockDDDD18<<<num_blocks, threads_per_block>>>(
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

        gpuSafe(gpuDeviceSynchronize());

        //omptimers[thread_id].stop("    K block DDDD");

        gpu::chunkedMemcpyDeviceToHost<double>(mat_K.data(), d_mat_K, pair_inds_count_for_K_dd);

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

        //omptimers[thread_id].stop("  K block DD");
    }

    }
    }  // end of compute K

    gpuSafe(gpuDeviceSynchronize());

#pragma omp barrier

    gpuSafe(gpuFree(d_data_boys_func));

    gpuSafe(gpuFree(d_data_spd_prim_info));
    gpuSafe(gpuFree(d_data_spd_prim_aoinds));

    // TODO: keep devptr
    screening.clear_devptr_double(gpu_id);
    screening.clear_devptr_uint32(gpu_id);

#pragma omp barrier

    omptimers[thread_id].stop("K compute");
    }

    timer.stop("Compute Fockmat");

    timer.start("Accumulate Fockmat");

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

    timer.stop("Accumulate Fockmat");

    timer.stop("Total timing");

    auto timer_summary = timer.getSummary();
    screening.setTimerSummary(timer_summary);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        auto gpu_timer_summary = omptimers[thread_id].getSummary();
        screening.setGpuTimerSummary(thread_id, gpu_timer_summary);
    }

    return mat_Fock_sum;
}

}  // namespace gpu
