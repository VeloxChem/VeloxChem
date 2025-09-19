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

#if defined(USE_CUDA)

#include <cublas_v2.h>
#include <cusolverDn.h>

#elif defined(USE_HIP)

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
//#include <hipsolver/hipsolver.h>
#include <magma_v2.h>

#endif

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

#include "Allocator.hpp"
#include "ScreeningData.hpp"
#include "BoysFuncTable.hpp"
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
#include "MatrixFunc.hpp"
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

struct EventWrapper
{
    gpuEvent_t event;

    inline void destroyEvent()
    {
        gpuSafe(gpuEventDestroy(event));
    }

    inline void createEvent()
    {
        gpuSafe(gpuEventCreateWithFlags(&event, gpuEventDisableTiming));
    }

    inline void waitForCompletion()
    {
        gpuSafe(gpuEventSynchronize(event));
    }

    inline void markEvent()
    {
        gpuSafe(gpuEventRecord(event));
    }

    inline void markStreamEvent(gpuStream_t stream)
    {
        gpuSafe(gpuEventRecordStream(event, stream));
    }

    inline void markStreamEventAndWait(gpuStream_t stream)
    {
        markStreamEvent(stream);
        waitForCompletion();

    }

};

struct StreamWrapper
{
    gpuStream_t stream;

    inline void destroyStream()
    {
        gpuSafe(gpuStreamDestroy(stream));
    }

    inline void createStream()
    {
        gpuSafe(gpuStreamCreate(&stream));
    }

    inline void createHighPriorityStream()
    {
        int priority;
        gpuSafe(gpuDeviceGetStreamPriorityRange(nullptr, &priority));
        gpuSafe(gpuStreamCreateWithPriority(&stream, gpuStreamDefault, priority));
    }
};

auto
computeQMatrixOnGPU(const CMolecule& molecule,
                    const CMolecularBasis& basis,
                    const CScreeningData& screening,
                    const int64_t rank,
                    const int64_t nnodes) -> CDenseMatrix
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

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    // TODO: double check by using different number of gpus per MPI
    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    gpuSafe(gpuMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), gpuMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    gpuSafe(gpuMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

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

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count_local; ij++)
        {
            const auto i = dd_first_inds_local[ij];
            const auto j = dd_second_inds_local[ij];

            mat_Q_omp.row(s_prim_count + p_prim_count * 3 + i)[s_prim_count + p_prim_count * 3 + j] = h_mat_Q[ij];

            if (i != j) mat_Q_omp.row(s_prim_count + p_prim_count * 3 + j)[s_prim_count + p_prim_count * 3 + i] = h_mat_Q[ij];
        }
    }

    gpuSafe(gpuDeviceSynchronize());

    gpuSafe(gpuFree(d_boys_func_table));
    gpuSafe(gpuFree(d_boys_func_ft));

    gpuSafe(gpuFree(d_s_prim_info));
    gpuSafe(gpuFree(d_p_prim_info));
    gpuSafe(gpuFree(d_d_prim_info));

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
computeOverlapAndKineticEnergyIntegralsOnGPU(const CMolecule& molecule,
                                             const CMolecularBasis& basis,
                                             const CScreeningData& screening,
                                             const int64_t rank,
                                             const int64_t nnodes) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

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

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

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

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    S_matrices[gpu_id].zero();
    T_matrices[gpu_id].zero();

    auto& mat_overlap = S_matrices[gpu_id];
    auto& mat_kinetic_energy = T_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

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

        gpuSafe(gpuMemcpy(mat_S.data(), d_mat_S, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_T.data(), d_mat_T, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_S.data(), d_mat_S, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_T.data(), d_mat_T, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_S.data(), d_mat_S, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_T.data(), d_mat_T, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_S.data(), d_mat_S, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_T.data(), d_mat_T, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_S.data(), d_mat_S, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_T.data(), d_mat_T, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_S.data(), d_mat_S, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_T.data(), d_mat_T, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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
                                      const CScreeningData& screening,
                                      const int64_t rank,
                                      const int64_t nnodes) -> CDenseMatrix
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

    return gpu::computePointChargesIntegralsOnGPU(molecule, basis, screening, points_info.data(), natoms, rank, nnodes);
}

auto
computePointChargesIntegralsOnGPU(const CMolecule& molecule,
                                  const CMolecularBasis& basis,
                                  const CScreeningData& screening,
                                  const double* points_info_ptr,
                                  const int64_t npoints,
                                  const int64_t rank,
                                  const int64_t nnodes) -> CDenseMatrix
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

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

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    // setting up the two device streams for Computation and Copying
    std::vector<StreamWrapper> streamList;
    streamList.reserve(2);
    // first stream, high priority, kernels
    streamList.emplace_back(StreamWrapper{});
    streamList.back().createHighPriorityStream();
    // second stream, default priority, data copying
    streamList.emplace_back(StreamWrapper{});
    streamList.back().createStream();

    // vector of hipEvent handles (one for each computation block)
    // once those get marked as completed, the code below can start
    // moving data back to the host to perform the post processing
    std::vector<EventWrapper> eventList;
    EventWrapper copyEvent;
    eventList.reserve(6);

    for (uint32_t i = 0; i < 6; i++)
    {
        eventList.emplace_back(EventWrapper{});
        eventList.back().createEvent();
    }
    copyEvent.createEvent();

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    gpuSafe(gpuMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), gpuMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    gpuSafe(gpuMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    V_matrices[gpu_id].zero();

    auto& mat_nuclear_potential = V_matrices[gpu_id];

    double *d_points_info;

    gpuSafe(gpuMalloc(&d_points_info, npoints * 4 * sizeof(double)));

    gpuSafe(gpuMemcpy(d_points_info, points_info_ptr, npoints * 4 * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuDeviceSynchronize());

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

        gpuSafe(gpuMemcpy(mat_V.data(), d_mat_V, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_V.data(), d_mat_V, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_V.data(), d_mat_V, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_V.data(), d_mat_V, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_V.data(), d_mat_V, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_V.data(), d_mat_V, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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
                                    const CScreeningData& screening,
                                    const int64_t rank,
                                    const int64_t nnodes) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

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

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

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

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    MX_matrices[gpu_id].zero();
    MY_matrices[gpu_id].zero();
    MZ_matrices[gpu_id].zero();

    auto& mat_mu_x = MX_matrices[gpu_id];
    auto& mat_mu_y = MY_matrices[gpu_id];
    auto& mat_mu_z = MZ_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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
                                    const CScreeningData& screening,
                                    const int64_t rank,
                                    const int64_t nnodes) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

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

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

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

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    MX_matrices[gpu_id].zero();
    MY_matrices[gpu_id].zero();
    MZ_matrices[gpu_id].zero();

    auto& mat_mu_x = MX_matrices[gpu_id];
    auto& mat_mu_y = MY_matrices[gpu_id];
    auto& mat_mu_z = MZ_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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
                                    const CScreeningData& screening,
                                    const int64_t rank,
                                    const int64_t nnodes) -> std::vector<CDenseMatrix>
{
    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

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

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

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

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;

    gpuSafe(gpuMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    // D gto block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    double*   d_d_prim_info;

    gpuSafe(gpuMalloc(&d_d_prim_info, d_prim_info.size() * sizeof(double)));

    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    gpuSafe(gpuMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds_local, sd_first_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds_local, pd_first_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds_local, dd_first_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_prim_pair_count_local * sizeof(uint32_t), gpuMemcpyHostToDevice));

    MX_matrices[gpu_id].zero();
    MY_matrices[gpu_id].zero();
    MZ_matrices[gpu_id].zero();

    auto& mat_mu_x = MX_matrices[gpu_id];
    auto& mat_mu_y = MY_matrices[gpu_id];
    auto& mat_mu_z = MZ_matrices[gpu_id];

    gpuSafe(gpuDeviceSynchronize());

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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

        gpuSafe(gpuMemcpy(mat_MX.data(), d_mat_MX, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MY.data(), d_mat_MY, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));
        gpuSafe(gpuMemcpy(mat_MZ.data(), d_mat_MZ, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost));

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
                 CScreeningData&    screening,
                 const int64_t      rank,
                 const int64_t      nnodes) -> CDenseMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Set device");

    // TODO sanity check for flag_K: SYMM or ANTISYMM

    CGpuDevices gpu_devices;

    const auto total_num_gpus_per_compute_node = gpu_devices.getNumberOfDevices();

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    // auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    auto gpu_rank = rank * num_gpus_per_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    timer.stop("Set device");

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

    CTimer prelink_timer;

    prelink_timer.start();

    // preLinK
    // J. Chem. Phys. 138, 134114 (2013)

    // TODO distribute computation of Q_prime

    auto mat_full = screening.get_mat_Q_full(s_prim_count, p_prim_count, d_prim_count);

    double *d_data_matrices_ABC;
    gpuSafe(gpuMalloc(&d_data_matrices_ABC, 3 * mat_full.getNumberOfElements() * sizeof(double)));

    double *d_matrix_A = d_data_matrices_ABC;
    double *d_matrix_B = d_matrix_A + mat_full.getNumberOfElements();
    double *d_matrix_C = d_matrix_B + mat_full.getNumberOfElements();

    gpuSafe(gpuMemcpy(d_matrix_A, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), gpuMemcpyHostToDevice));

    mat_full = screening.get_mat_D_abs_full(s_prim_count, p_prim_count, d_prim_count, s_prim_aoinds, p_prim_aoinds, d_prim_aoinds, cart_naos, cart_dens_ptr);

    gpuSafe(gpuMemcpy(d_matrix_B, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), gpuMemcpyHostToDevice));

    const auto all_prim_count = mat_full.getNumberOfRows();

#if defined(USE_CUDA)

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

#elif defined(USE_HIP)

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

#endif

    gpuSafe(gpuFree(d_data_matrices_ABC));

    timer.stop("Prep. Q_prime");

    timer.start("Prep. preLinK 1");

    prelink_timer.stop();

    auto prelink_elapsed_time = prelink_timer.getElapsedTime();

    screening.setPreLinkTime(prelink_elapsed_time);

    screening.form_pair_inds_for_K(s_prim_count, p_prim_count, d_prim_count, mat_full, prelink_threshold);

    timer.stop("Prep. preLinK 1");

    timer.start("Prep. preLinK 2");

    mat_full = CDenseMatrix();

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

    screening.initTimers(num_gpus_per_node);

    timer.stop("Prep. Fockmat");

    std::vector<CMultiTimer> omptimers(nthreads);

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id < num_gpus_per_node)
    {
    auto gpu_id = thread_id;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    // auto gpu_count = nnodes * num_gpus_per_node;

    omptimers[thread_id].start("Set device");

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    omptimers[thread_id].stop("Set device");


    // setting up the two device streams for Coulomb and Exchange
    std::vector<StreamWrapper> streamList;
    streamList.reserve(3);
    // first stream, high priority, Coulomb
    streamList.emplace_back(StreamWrapper{});
    streamList.back().createHighPriorityStream();
    // second stream, high priority, Exchange
    streamList.emplace_back(StreamWrapper{});
    streamList.back().createHighPriorityStream();
    // third stream, default priority, data copying
    streamList.emplace_back(StreamWrapper{});
    streamList.back().createStream();

    // vector of hipEvent handles (one for each computation block)
    // once those get marked as completed, the code below can start
    // moving data back to the host to perform the post processing
    std::vector<EventWrapper> eventListCoulomb;
    std::vector<EventWrapper> eventListExchange;
    eventListCoulomb.reserve(6);
    eventListExchange.reserve(6);
    for (uint32_t i = 0; i < 6; i++)
    {
        eventListCoulomb.emplace_back(EventWrapper{});
        eventListCoulomb.back().createEvent();
        eventListExchange.emplace_back(EventWrapper{});
        eventListExchange.back().createEvent();
    }
    // extra event for handling copying back the data
    EventWrapper copyEvent;
    copyEvent.createEvent();

    omptimers[thread_id].start("Boys func. prep.");

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();
    const auto boys_func_ft    = boysfunc::getBoysFuncFactors();

    double *d_data_boys_func;
    gpuSafe(gpuMalloc(&d_data_boys_func, (boys_func_table.size() + boys_func_table.size()) * sizeof(double)));

    double* d_boys_func_table = d_data_boys_func;
    double* d_boys_func_ft = d_boys_func_table + boys_func_table.size();

    gpuSafe(gpuMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), gpuMemcpyHostToDevice));

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

    double*   d_data_spd_prim_info;
    gpuSafe(gpuMalloc(&d_data_spd_prim_info, (s_prim_info.size() + p_prim_info.size() + d_prim_info.size()) * sizeof(double)));

    uint32_t* d_data_spd_prim_aoinds;
    gpuSafe(gpuMalloc(&d_data_spd_prim_aoinds, (s_prim_aoinds.size() + p_prim_aoinds.size() + d_prim_aoinds.size())* sizeof(uint32_t)));

    double*   d_s_prim_info = d_data_spd_prim_info;
    double*   d_p_prim_info = d_s_prim_info + s_prim_info.size();
    double*   d_d_prim_info = d_p_prim_info + p_prim_info.size();

    uint32_t* d_s_prim_aoinds = d_data_spd_prim_aoinds;
    uint32_t* d_p_prim_aoinds = d_s_prim_aoinds + s_prim_aoinds.size();
    uint32_t* d_d_prim_aoinds = d_p_prim_aoinds + p_prim_aoinds.size();

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_d_prim_aoinds, d_prim_aoinds.data(), d_prim_aoinds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    omptimers[thread_id].stop("GTO block prep.");

    // GTO block pairs

    omptimers[thread_id].start("Coulomb prep.");

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

    VeloxHostVector<double> mat_J(max_prim_pair_count_local);

    // sorted Q, D, and indices on device

    const uint32_t numCalculationBlocksCoulomb = 6;

    double *d_data_mat_D_J;
    gpuSafe(gpuMalloc(&d_data_mat_D_J, (max_prim_pair_count + max_prim_pair_count_local) * sizeof(double) * numCalculationBlocksCoulomb));

    double *d_mat_D = d_data_mat_D_J;
    double *d_mat_J = d_mat_D + (max_prim_pair_count * numCalculationBlocksCoulomb);

    double *d_data_mat_Q;
    gpuSafe(gpuMalloc(&d_data_mat_Q, (ss_prim_pair_count +
                                      sp_prim_pair_count +
                                      sd_prim_pair_count +
                                      pp_prim_pair_count +
                                      pd_prim_pair_count +
                                      dd_prim_pair_count) * sizeof(double)));

    double *d_ss_mat_Q = d_data_mat_Q;
    double *d_sp_mat_Q = d_ss_mat_Q + ss_prim_pair_count;
    double *d_sd_mat_Q = d_sp_mat_Q + sp_prim_pair_count;
    double *d_pp_mat_Q = d_sd_mat_Q + sd_prim_pair_count;
    double *d_pd_mat_Q = d_pp_mat_Q + pp_prim_pair_count;
    double *d_dd_mat_Q = d_pd_mat_Q + pd_prim_pair_count;

    uint32_t *d_data_first_second_inds;
    gpuSafe(gpuMalloc(&d_data_first_second_inds, (ss_prim_pair_count +
                                                  ss_prim_pair_count +
                                                  sp_prim_pair_count +
                                                  sp_prim_pair_count +
                                                  sd_prim_pair_count +
                                                  sd_prim_pair_count +
                                                  pp_prim_pair_count +
                                                  pp_prim_pair_count +
                                                  pd_prim_pair_count +
                                                  pd_prim_pair_count +
                                                  dd_prim_pair_count +
                                                  dd_prim_pair_count) * sizeof(uint32_t)));

    uint32_t *d_ss_first_inds  = d_data_first_second_inds;
    uint32_t *d_ss_second_inds = d_ss_first_inds  + ss_prim_pair_count;
    uint32_t *d_sp_first_inds  = d_ss_second_inds + ss_prim_pair_count;
    uint32_t *d_sp_second_inds = d_sp_first_inds  + sp_prim_pair_count;
    uint32_t *d_sd_first_inds  = d_sp_second_inds + sp_prim_pair_count;
    uint32_t *d_sd_second_inds = d_sd_first_inds  + sd_prim_pair_count;
    uint32_t *d_pp_first_inds  = d_sd_second_inds + sd_prim_pair_count;
    uint32_t *d_pp_second_inds = d_pp_first_inds  + pp_prim_pair_count;
    uint32_t *d_pd_first_inds  = d_pp_second_inds + pp_prim_pair_count;
    uint32_t *d_pd_second_inds = d_pd_first_inds  + pd_prim_pair_count;
    uint32_t *d_dd_first_inds  = d_pd_second_inds + pd_prim_pair_count;
    uint32_t *d_dd_second_inds = d_dd_first_inds  + dd_prim_pair_count;

    double *d_data_pair_data;
    gpuSafe(gpuMalloc(&d_data_pair_data, (ss_pair_data.size() +
                                          sp_pair_data.size() +
                                          sd_pair_data.size() +
                                          pp_pair_data.size() +
                                          pd_pair_data.size() +
                                          dd_pair_data.size()) * sizeof(double)));

    double *d_ss_pair_data = d_data_pair_data;
    double *d_sp_pair_data = d_ss_pair_data + ss_pair_data.size();
    double *d_sd_pair_data = d_sp_pair_data + sp_pair_data.size();
    double *d_pp_pair_data = d_sd_pair_data + sd_pair_data.size();
    double *d_pd_pair_data = d_pp_pair_data + pp_pair_data.size();
    double *d_dd_pair_data = d_pd_pair_data + pd_pair_data.size();

    double *d_data_mat_Q_local;
    gpuSafe(gpuMalloc(&d_data_mat_Q_local, (ss_prim_pair_count_local +
                                            sp_prim_pair_count_local +
                                            sd_prim_pair_count_local +
                                            pp_prim_pair_count_local +
                                            pd_prim_pair_count_local +
                                            dd_prim_pair_count_local) * sizeof(double)));

    double *d_ss_mat_Q_local = d_data_mat_Q_local;
    double *d_sp_mat_Q_local = d_ss_mat_Q_local + ss_prim_pair_count_local;
    double *d_sd_mat_Q_local = d_sp_mat_Q_local + sp_prim_pair_count_local;
    double *d_pp_mat_Q_local = d_sd_mat_Q_local + sd_prim_pair_count_local;
    double *d_pd_mat_Q_local = d_pp_mat_Q_local + pp_prim_pair_count_local;
    double *d_dd_mat_Q_local = d_pd_mat_Q_local + pd_prim_pair_count_local;

    uint32_t *d_data_first_second_inds_local;
    gpuSafe(gpuMalloc(&d_data_first_second_inds_local, (ss_prim_pair_count_local +
                                                        ss_prim_pair_count_local +
                                                        sp_prim_pair_count_local +
                                                        sp_prim_pair_count_local +
                                                        sd_prim_pair_count_local +
                                                        sd_prim_pair_count_local +
                                                        pp_prim_pair_count_local +
                                                        pp_prim_pair_count_local +
                                                        pd_prim_pair_count_local +
                                                        pd_prim_pair_count_local +
                                                        dd_prim_pair_count_local +
                                                        dd_prim_pair_count_local) * sizeof(uint32_t)));

    uint32_t *d_ss_first_inds_local  = d_data_first_second_inds_local;
    uint32_t *d_ss_second_inds_local = d_ss_first_inds_local  + ss_prim_pair_count_local;
    uint32_t *d_sp_first_inds_local  = d_ss_second_inds_local + ss_prim_pair_count_local;
    uint32_t *d_sp_second_inds_local = d_sp_first_inds_local  + sp_prim_pair_count_local;
    uint32_t *d_sd_first_inds_local  = d_sp_second_inds_local + sp_prim_pair_count_local;
    uint32_t *d_sd_second_inds_local = d_sd_first_inds_local  + sd_prim_pair_count_local;
    uint32_t *d_pp_first_inds_local  = d_sd_second_inds_local + sd_prim_pair_count_local;
    uint32_t *d_pp_second_inds_local = d_pp_first_inds_local  + pp_prim_pair_count_local;
    uint32_t *d_pd_first_inds_local  = d_pp_second_inds_local + pp_prim_pair_count_local;
    uint32_t *d_pd_second_inds_local = d_pd_first_inds_local  + pd_prim_pair_count_local;
    uint32_t *d_dd_first_inds_local  = d_pd_second_inds_local + pd_prim_pair_count_local;
    uint32_t *d_dd_second_inds_local = d_dd_first_inds_local  + dd_prim_pair_count_local;

    double *d_data_pair_data_local;
    gpuSafe(gpuMalloc(&d_data_pair_data_local, (ss_pair_data_local.size() +
                                                sp_pair_data_local.size() +
                                                sd_pair_data_local.size() +
                                                pp_pair_data_local.size() +
                                                pd_pair_data_local.size() +
                                                dd_pair_data_local.size()) * sizeof(double)));

    double *d_ss_pair_data_local = d_data_pair_data_local;
    double *d_sp_pair_data_local = d_ss_pair_data_local + ss_pair_data_local.size();
    double *d_sd_pair_data_local = d_sp_pair_data_local + sp_pair_data_local.size();
    double *d_pp_pair_data_local = d_sd_pair_data_local + sd_pair_data_local.size();
    double *d_pd_pair_data_local = d_pp_pair_data_local + pp_pair_data_local.size();
    double *d_dd_pair_data_local = d_pd_pair_data_local + pd_pair_data_local.size();

    {
        const uint32_t offset = max_prim_pair_count * 0;
        gpuSafe(gpuMemcpy(d_mat_D + offset, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));
    }
    {
        const uint32_t offset = max_prim_pair_count * 1;
        gpuSafe(gpuMemcpy(d_mat_D + offset, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));
    }
    {
        const uint32_t offset = max_prim_pair_count * 2;
        gpuSafe(gpuMemcpy(d_mat_D + offset, sd_mat_D.data(), sd_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));
    }
    {
        const uint32_t offset = max_prim_pair_count * 3;
        gpuSafe(gpuMemcpy(d_mat_D + offset, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));
    }
    {
        const uint32_t offset = max_prim_pair_count * 4;
        gpuSafe(gpuMemcpy(d_mat_D + offset, pd_mat_D.data(), pd_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));
    }
    {
        const uint32_t offset = max_prim_pair_count * 5;
        gpuSafe(gpuMemcpy(d_mat_D + offset, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));
    }

    gpuSafe(gpuMemcpy(d_ss_mat_Q, ss_mat_Q.data(), ss_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_mat_Q, sp_mat_Q.data(), sp_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_mat_Q, sd_mat_Q.data(), sd_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_mat_Q, pp_mat_Q.data(), pp_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_mat_Q, pd_mat_Q.data(), pd_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_mat_Q, dd_mat_Q.data(), dd_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_ss_first_inds,  ss_first_inds.data(),  ss_first_inds.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds, ss_second_inds.data(), ss_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_first_inds,  sp_first_inds.data(),  sp_first_inds.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds, sp_second_inds.data(), sp_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_first_inds,  sd_first_inds.data(),  sd_first_inds.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds, sd_second_inds.data(), sd_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_first_inds,  pp_first_inds.data(),  pp_first_inds.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds, pp_second_inds.data(), pp_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_first_inds,  pd_first_inds.data(),  pd_first_inds.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds, pd_second_inds.data(), pd_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_first_inds,  dd_first_inds.data(),  dd_first_inds.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds, dd_second_inds.data(), dd_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_ss_pair_data, ss_pair_data.data(), ss_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_pair_data, sp_pair_data.data(), sp_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_pair_data, sd_pair_data.data(), sd_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_pair_data, pp_pair_data.data(), pp_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_pair_data, pd_pair_data.data(), pd_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_pair_data, dd_pair_data.data(), dd_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_ss_mat_Q_local, ss_mat_Q_local.data(), ss_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_mat_Q_local, sp_mat_Q_local.data(), sp_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_mat_Q_local, sd_mat_Q_local.data(), sd_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_mat_Q_local, pp_mat_Q_local.data(), pp_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_mat_Q_local, pd_mat_Q_local.data(), pd_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_mat_Q_local, dd_mat_Q_local.data(), dd_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_ss_first_inds_local,  ss_first_inds_local.data(),  ss_first_inds_local.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_second_inds_local.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_first_inds_local,  sp_first_inds_local.data(),  sp_first_inds_local.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_second_inds_local.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_first_inds_local,  sd_first_inds_local.data(),  sd_first_inds_local.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds_local, sd_second_inds_local.data(), sd_second_inds_local.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_first_inds_local,  pp_first_inds_local.data(),  pp_first_inds_local.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_second_inds_local.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_first_inds_local,  pd_first_inds_local.data(),  pd_first_inds_local.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds_local, pd_second_inds_local.data(), pd_second_inds_local.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_first_inds_local,  dd_first_inds_local.data(),  dd_first_inds_local.size()  * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds_local, dd_second_inds_local.data(), dd_second_inds_local.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_ss_pair_data_local, ss_pair_data_local.data(), ss_pair_data_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_pair_data_local, sp_pair_data_local.data(), sp_pair_data_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_pair_data_local, sd_pair_data_local.data(), sd_pair_data_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_pair_data_local, pp_pair_data_local.data(), pp_pair_data_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_pair_data_local, pd_pair_data_local.data(), pd_pair_data_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_pair_data_local, dd_pair_data_local.data(), dd_pair_data_local.size() * sizeof(double), gpuMemcpyHostToDevice));

    omptimers[thread_id].stop("Coulomb prep.");

    omptimers[thread_id].start("Exchange prep.");

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

    VeloxHostVector<double> mat_K(max_pair_inds_count);

    const uint32_t numCalculationBlocksExchange = 6;

    double*   d_mat_K;
    gpuSafe(gpuMalloc(&d_mat_K, max_pair_inds_count * sizeof(double) * numCalculationBlocksExchange));

    uint32_t *d_data_pair_inds_for_K;
    gpuSafe(gpuMalloc(&d_data_pair_inds_for_K, (pair_inds_count_for_K_ss +
                                                pair_inds_count_for_K_ss +
                                                pair_inds_count_for_K_sp +
                                                pair_inds_count_for_K_sp +
                                                pair_inds_count_for_K_sd +
                                                pair_inds_count_for_K_sd +
                                                pair_inds_count_for_K_pp +
                                                pair_inds_count_for_K_pp +
                                                pair_inds_count_for_K_pd +
                                                pair_inds_count_for_K_pd +
                                                pair_inds_count_for_K_dd +
                                                pair_inds_count_for_K_dd) * sizeof(uint32_t)));

    uint32_t *d_pair_inds_i_for_K_ss = d_data_pair_inds_for_K;
    uint32_t *d_pair_inds_k_for_K_ss = d_pair_inds_i_for_K_ss + pair_inds_count_for_K_ss;
    uint32_t *d_pair_inds_i_for_K_sp = d_pair_inds_k_for_K_ss + pair_inds_count_for_K_ss;
    uint32_t *d_pair_inds_k_for_K_sp = d_pair_inds_i_for_K_sp + pair_inds_count_for_K_sp;
    uint32_t *d_pair_inds_i_for_K_sd = d_pair_inds_k_for_K_sp + pair_inds_count_for_K_sp;
    uint32_t *d_pair_inds_k_for_K_sd = d_pair_inds_i_for_K_sd + pair_inds_count_for_K_sd;
    uint32_t *d_pair_inds_i_for_K_pp = d_pair_inds_k_for_K_sd + pair_inds_count_for_K_sd;
    uint32_t *d_pair_inds_k_for_K_pp = d_pair_inds_i_for_K_pp + pair_inds_count_for_K_pp;
    uint32_t *d_pair_inds_i_for_K_pd = d_pair_inds_k_for_K_pp + pair_inds_count_for_K_pp;
    uint32_t *d_pair_inds_k_for_K_pd = d_pair_inds_i_for_K_pd + pair_inds_count_for_K_pd;
    uint32_t *d_pair_inds_i_for_K_dd = d_pair_inds_k_for_K_pd + pair_inds_count_for_K_pd;
    uint32_t *d_pair_inds_k_for_K_dd = d_pair_inds_i_for_K_dd + pair_inds_count_for_K_dd;

    gpuSafe(gpuMemcpy(d_pair_inds_i_for_K_ss, pair_inds_i_for_K_ss.data(), pair_inds_i_for_K_ss.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_k_for_K_ss, pair_inds_k_for_K_ss.data(), pair_inds_k_for_K_ss.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_i_for_K_sp, pair_inds_i_for_K_sp.data(), pair_inds_i_for_K_sp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_k_for_K_sp, pair_inds_k_for_K_sp.data(), pair_inds_k_for_K_sp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_i_for_K_sd, pair_inds_i_for_K_sd.data(), pair_inds_i_for_K_sd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_k_for_K_sd, pair_inds_k_for_K_sd.data(), pair_inds_k_for_K_sd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_i_for_K_pp, pair_inds_i_for_K_pp.data(), pair_inds_i_for_K_pp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_k_for_K_pp, pair_inds_k_for_K_pp.data(), pair_inds_k_for_K_pp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_i_for_K_pd, pair_inds_i_for_K_pd.data(), pair_inds_i_for_K_pd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_k_for_K_pd, pair_inds_k_for_K_pd.data(), pair_inds_k_for_K_pd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_i_for_K_dd, pair_inds_i_for_K_dd.data(), pair_inds_i_for_K_dd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_inds_k_for_K_dd, pair_inds_k_for_K_dd.data(), pair_inds_k_for_K_dd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    double* d_mat_D_full_AO;
    gpuSafe(gpuMalloc(&d_mat_D_full_AO, cart_naos * cart_naos * sizeof(double)));

    double *d_data_Q_K;
    gpuSafe(gpuMalloc(&d_data_Q_K, (Q_K_ss.size() +
                                    Q_K_sp.size() +
                                    Q_K_ps.size() +
                                    Q_K_sd.size() +
                                    Q_K_ds.size() +
                                    Q_K_pp.size() +
                                    Q_K_pd.size() +
                                    Q_K_dp.size() +
                                    Q_K_dd.size()) * sizeof(double)));

    double *d_Q_K_ss = d_data_Q_K;
    double *d_Q_K_sp = d_Q_K_ss + Q_K_ss.size();
    double *d_Q_K_ps = d_Q_K_sp + Q_K_sp.size();
    double *d_Q_K_sd = d_Q_K_ps + Q_K_ps.size();
    double *d_Q_K_ds = d_Q_K_sd + Q_K_sd.size();
    double *d_Q_K_pp = d_Q_K_ds + Q_K_ds.size();
    double *d_Q_K_pd = d_Q_K_pp + Q_K_pp.size();
    double *d_Q_K_dp = d_Q_K_pd + Q_K_pd.size();
    double *d_Q_K_dd = d_Q_K_dp + Q_K_dp.size();

    uint32_t *d_data_D_inds_K;
    gpuSafe(gpuMalloc(&d_data_D_inds_K, (D_inds_K_ss.size() +
                                         D_inds_K_sp.size() +
                                         D_inds_K_ps.size() +
                                         D_inds_K_sd.size() +
                                         D_inds_K_ds.size() +
                                         D_inds_K_pp.size() +
                                         D_inds_K_pd.size() +
                                         D_inds_K_dp.size() +
                                         D_inds_K_dd.size()) * sizeof(uint32_t)));

    uint32_t *d_D_inds_K_ss = d_data_D_inds_K;
    uint32_t *d_D_inds_K_sp = d_D_inds_K_ss + D_inds_K_ss.size();
    uint32_t *d_D_inds_K_ps = d_D_inds_K_sp + D_inds_K_sp.size();
    uint32_t *d_D_inds_K_sd = d_D_inds_K_ps + D_inds_K_ps.size();
    uint32_t *d_D_inds_K_ds = d_D_inds_K_sd + D_inds_K_sd.size();
    uint32_t *d_D_inds_K_pp = d_D_inds_K_ds + D_inds_K_ds.size();
    uint32_t *d_D_inds_K_pd = d_D_inds_K_pp + D_inds_K_pp.size();
    uint32_t *d_D_inds_K_dp = d_D_inds_K_pd + D_inds_K_pd.size();
    uint32_t *d_D_inds_K_dd = d_D_inds_K_dp + D_inds_K_dp.size();

    uint32_t *d_data_pair_counts_displs_K;
    gpuSafe(gpuMalloc(&d_data_pair_counts_displs_K, (pair_counts_K_ss.size() +
                                                     pair_counts_K_sp.size() +
                                                     pair_counts_K_ps.size() +
                                                     pair_counts_K_sd.size() +
                                                     pair_counts_K_ds.size() +
                                                     pair_counts_K_pp.size() +
                                                     pair_counts_K_pd.size() +
                                                     pair_counts_K_dp.size() +
                                                     pair_counts_K_dd.size() +
                                                     pair_displs_K_ss.size() +
                                                     pair_displs_K_sp.size() +
                                                     pair_displs_K_ps.size() +
                                                     pair_displs_K_sd.size() +
                                                     pair_displs_K_ds.size() +
                                                     pair_displs_K_pp.size() +
                                                     pair_displs_K_pd.size() +
                                                     pair_displs_K_dp.size() +
                                                     pair_displs_K_dd.size()) * sizeof(uint32_t)));

    uint32_t *d_pair_counts_K_ss = d_data_pair_counts_displs_K;
    uint32_t *d_pair_counts_K_sp = d_pair_counts_K_ss + pair_counts_K_ss.size();
    uint32_t *d_pair_counts_K_ps = d_pair_counts_K_sp + pair_counts_K_sp.size();
    uint32_t *d_pair_counts_K_sd = d_pair_counts_K_ps + pair_counts_K_ps.size();
    uint32_t *d_pair_counts_K_ds = d_pair_counts_K_sd + pair_counts_K_sd.size();
    uint32_t *d_pair_counts_K_pp = d_pair_counts_K_ds + pair_counts_K_ds.size();
    uint32_t *d_pair_counts_K_pd = d_pair_counts_K_pp + pair_counts_K_pp.size();
    uint32_t *d_pair_counts_K_dp = d_pair_counts_K_pd + pair_counts_K_pd.size();
    uint32_t *d_pair_counts_K_dd = d_pair_counts_K_dp + pair_counts_K_dp.size();
    uint32_t *d_pair_displs_K_ss = d_pair_counts_K_dd + pair_counts_K_dd.size();
    uint32_t *d_pair_displs_K_sp = d_pair_displs_K_ss + pair_displs_K_ss.size();
    uint32_t *d_pair_displs_K_ps = d_pair_displs_K_sp + pair_displs_K_sp.size();
    uint32_t *d_pair_displs_K_sd = d_pair_displs_K_ps + pair_displs_K_ps.size();
    uint32_t *d_pair_displs_K_ds = d_pair_displs_K_sd + pair_displs_K_sd.size();
    uint32_t *d_pair_displs_K_pp = d_pair_displs_K_ds + pair_displs_K_ds.size();
    uint32_t *d_pair_displs_K_pd = d_pair_displs_K_pp + pair_displs_K_pp.size();
    uint32_t *d_pair_displs_K_dp = d_pair_displs_K_pd + pair_displs_K_pd.size();
    uint32_t *d_pair_displs_K_dd = d_pair_displs_K_dp + pair_displs_K_dp.size();

    double *d_data_pair_data_K;
    gpuSafe(gpuMalloc(&d_data_pair_data_K, (pair_data_K_ss.size() +
                                            pair_data_K_sp.size() +
                                            pair_data_K_ps.size() +
                                            pair_data_K_sd.size() +
                                            pair_data_K_ds.size() +
                                            pair_data_K_pp.size() +
                                            pair_data_K_pd.size() +
                                            pair_data_K_dp.size() +
                                            pair_data_K_dd.size()) * sizeof(double)));

    double *d_pair_data_K_ss = d_data_pair_data_K;
    double *d_pair_data_K_sp = d_pair_data_K_ss + pair_data_K_ss.size();
    double *d_pair_data_K_ps = d_pair_data_K_sp + pair_data_K_sp.size();
    double *d_pair_data_K_sd = d_pair_data_K_ps + pair_data_K_ps.size();
    double *d_pair_data_K_ds = d_pair_data_K_sd + pair_data_K_sd.size();
    double *d_pair_data_K_pp = d_pair_data_K_ds + pair_data_K_ds.size();
    double *d_pair_data_K_pd = d_pair_data_K_pp + pair_data_K_pp.size();
    double *d_pair_data_K_dp = d_pair_data_K_pd + pair_data_K_pd.size();
    double *d_pair_data_K_dd = d_pair_data_K_dp + pair_data_K_dp.size();

    gpuSafe(gpuMemcpy(d_mat_D_full_AO, cart_dens_ptr, cart_naos * cart_naos * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_Q_K_ss, Q_K_ss.data(), Q_K_ss.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_sp, Q_K_sp.data(), Q_K_sp.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_ps, Q_K_ps.data(), Q_K_ps.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_sd, Q_K_sd.data(), Q_K_sd.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_ds, Q_K_ds.data(), Q_K_ds.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_pp, Q_K_pp.data(), Q_K_pp.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_pd, Q_K_pd.data(), Q_K_pd.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_dp, Q_K_dp.data(), Q_K_dp.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_Q_K_dd, Q_K_dd.data(), Q_K_dd.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_D_inds_K_ss, D_inds_K_ss.data(), D_inds_K_ss.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_sp, D_inds_K_sp.data(), D_inds_K_sp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_ps, D_inds_K_ps.data(), D_inds_K_ps.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_sd, D_inds_K_sd.data(), D_inds_K_sd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_ds, D_inds_K_ds.data(), D_inds_K_ds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_pp, D_inds_K_pp.data(), D_inds_K_pp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_pd, D_inds_K_pd.data(), D_inds_K_pd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_dp, D_inds_K_dp.data(), D_inds_K_dp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_D_inds_K_dd, D_inds_K_dd.data(), D_inds_K_dd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pair_displs_K_ss, pair_displs_K_ss.data(), pair_displs_K_ss.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_sp, pair_displs_K_sp.data(), pair_displs_K_sp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_ps, pair_displs_K_ps.data(), pair_displs_K_ps.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_sd, pair_displs_K_sd.data(), pair_displs_K_sd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_ds, pair_displs_K_ds.data(), pair_displs_K_ds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_pp, pair_displs_K_pp.data(), pair_displs_K_pp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_pd, pair_displs_K_pd.data(), pair_displs_K_pd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_dp, pair_displs_K_dp.data(), pair_displs_K_dp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_displs_K_dd, pair_displs_K_dd.data(), pair_displs_K_dd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pair_counts_K_ss, pair_counts_K_ss.data(), pair_counts_K_ss.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_sp, pair_counts_K_sp.data(), pair_counts_K_sp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_ps, pair_counts_K_ps.data(), pair_counts_K_ps.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_sd, pair_counts_K_sd.data(), pair_counts_K_sd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_ds, pair_counts_K_ds.data(), pair_counts_K_ds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_pp, pair_counts_K_pp.data(), pair_counts_K_pp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_pd, pair_counts_K_pd.data(), pair_counts_K_pd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_dp, pair_counts_K_dp.data(), pair_counts_K_dp.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_counts_K_dd, pair_counts_K_dd.data(), pair_counts_K_dd.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pair_data_K_ss, pair_data_K_ss.data(), pair_data_K_ss.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_sp, pair_data_K_sp.data(), pair_data_K_sp.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_ps, pair_data_K_ps.data(), pair_data_K_ps.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_sd, pair_data_K_sd.data(), pair_data_K_sd.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_ds, pair_data_K_ds.data(), pair_data_K_ds.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_pp, pair_data_K_pp.data(), pair_data_K_pp.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_pd, pair_data_K_pd.data(), pair_data_K_pd.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_dp, pair_data_K_dp.data(), pair_data_K_dp.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pair_data_K_dd, pair_data_K_dd.data(), pair_data_K_dd.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuDeviceSynchronize());

    omptimers[thread_id].stop("Exchange prep.");

    mat_Fock_omp[gpu_id].zero();

    omptimers[thread_id].start("J computation");


    CTimer coulomb_timer;

    coulomb_timer.start();

    // compute J

    if (std::fabs(prefac_coulomb) > 1.0e-13)
    {

    {
        omptimers[thread_id].start("  J block SS");
        const dim3 threads_per_block(TILE_DIM * TILE_DIM);

        const dim3 num_blocks(((max_prim_pair_count_local * numCalculationBlocksCoulomb) + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(d_mat_J, static_cast<uint32_t>(max_prim_pair_count_local * numCalculationBlocksCoulomb));
    }

    // J: S-S block

    if (ss_prim_pair_count_local > 0)
    {
        timer.start("  J block SS");

        const uint32_t offsetJ = 0 * max_prim_pair_count_local;

        // set up thread blocks for J

        const dim3 threads_per_block = dim3(TILE_DIM, TILE_DIM);

        const dim3 num_blocks = dim3((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SS|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            const uint32_t offsetD = 0 * max_prim_pair_count;

            //omptimers[thread_id].start("    J block SSSS");

            gpu::computeCoulombFockSSSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block SSSS");
        }

        // J: (SS|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 1 * max_prim_pair_count;

            gpu::computeCoulombFockSSSP<<<num_blocks, threads_per_block,0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SS|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 2 * max_prim_pair_count;

            gpu::computeCoulombFockSSSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SS|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 3 * max_prim_pair_count;

            gpu::computeCoulombFockSSPP<<<num_blocks, threads_per_block,0 , streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SS|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 4 * max_prim_pair_count;

            gpu::computeCoulombFockSSPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SS|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 5 * max_prim_pair_count;

            gpu::computeCoulombFockSSDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        eventListCoulomb[0].markStreamEvent(streamList[0].stream);

        omptimers[thread_id].stop("  J block SS");
    }

    // J: S-P block

    if (sp_prim_pair_count_local > 0)
    {
        omptimers[thread_id].start("  J block SP");

        const uint32_t offsetJ = 1 * max_prim_pair_count_local;

        // set up thread blocks for J

        const dim3 threads_per_block = dim3(TILE_DIM, TILE_DIM);

        const dim3 num_blocks = dim3((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            const uint32_t offsetD = 0 * max_prim_pair_count;
            gpu::computeCoulombFockSPSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 1 * max_prim_pair_count;
            gpu::computeCoulombFockSPSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 2 * max_prim_pair_count;
            gpu::computeCoulombFockSPSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 3 * max_prim_pair_count;
            gpu::computeCoulombFockSPPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 4 * max_prim_pair_count;
            gpu::computeCoulombFockSPPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 5 * max_prim_pair_count;
            gpu::computeCoulombFockSPDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        eventListCoulomb[1].markStreamEvent(streamList[0].stream);

        omptimers[thread_id].stop("  J block SP");
    }

    // J: P-P block

    if (pp_prim_pair_count_local > 0)
    {
        omptimers[thread_id].start("  J block PP");

        const uint32_t offsetJ = 2 * max_prim_pair_count_local;

        // set up thread blocks for J

        const dim3 threads_per_block = dim3(TILE_DIM, TILE_DIM);

        const dim3 num_blocks = dim3((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            const uint32_t offsetD = 0 * max_prim_pair_count;
            gpu::computeCoulombFockPPSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 1 * max_prim_pair_count;
            gpu::computeCoulombFockPPSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PP|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 2 * max_prim_pair_count;
            gpu::computeCoulombFockPPSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 3 * max_prim_pair_count;
            gpu::computeCoulombFockPPPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PP|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 4 * max_prim_pair_count;
        gpu::computeCoulombFockPPPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PP|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 5 * max_prim_pair_count;
            omptimers[thread_id].start("    J block PPDD");

            gpu::computeCoulombFockPPDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block PPDD");
        }

        eventListCoulomb[2].markStreamEvent(streamList[0].stream);

        omptimers[thread_id].stop("  J block PP");
    }

    // J: S-D block

    if (sd_prim_pair_count_local > 0)
    {
        omptimers[thread_id].start("  J block SD");

        const uint32_t offsetJ = 3 * max_prim_pair_count_local;

        // set up thread blocks for J

        const dim3 threads_per_block = dim3(TILE_DIM, TILE_DIM);

        const dim3 num_blocks = dim3((sd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            const uint32_t offsetD = 0 * max_prim_pair_count;
            gpu::computeCoulombFockSDSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 1 * max_prim_pair_count;
            gpu::computeCoulombFockSDSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 2 * max_prim_pair_count;
            gpu::computeCoulombFockSDSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 3 * max_prim_pair_count;
            gpu::computeCoulombFockSDPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 4 * max_prim_pair_count;
            gpu::computeCoulombFockSDPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (SD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 5 * max_prim_pair_count;
            omptimers[thread_id].start("    J block SDDD");

            gpu::computeCoulombFockSDDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block SDDD");
        }

        eventListCoulomb[3].markStreamEvent(streamList[0].stream);

        omptimers[thread_id].stop("  J block SD");
    }

    // J: P-D block

    if (pd_prim_pair_count_local > 0)
    {
        omptimers[thread_id].start("  J block PD");

        const uint32_t offsetJ = 4 * max_prim_pair_count_local;

        // set up thread blocks for J

        const dim3 threads_per_block = dim3(TILE_DIM, TILE_DIM);

        const dim3 num_blocks = dim3((pd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            const uint32_t offsetD = 0 * max_prim_pair_count;
            gpu::computeCoulombFockPDSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 1 * max_prim_pair_count;
            gpu::computeCoulombFockPDSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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
        }

        // J: (PD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 2 * max_prim_pair_count;
            omptimers[thread_id].start("    J block PDSD");

            gpu::computeCoulombFockPDSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block PDSD");
        }

        // J: (PD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 3 * max_prim_pair_count;
            omptimers[thread_id].start("    J block PDPP");

            gpu::computeCoulombFockPDPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block PDPP");
        }

        // J: (PD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 4 * max_prim_pair_count;
            omptimers[thread_id].start("    J block PDPD");

            gpu::computeCoulombFockPDPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block PDPD");
        }

        // J: (PD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 5 * max_prim_pair_count;
            omptimers[thread_id].start("    J block PDDD");

            gpu::computeCoulombFockPDDD0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockPDDD1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockPDDD2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockPDDD3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockPDDD4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockPDDD5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockPDDD6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block PDDD");
        }

        eventListCoulomb[4].markStreamEvent(streamList[0].stream);

        omptimers[thread_id].stop("  J block PD");
    }

    // J: D-D block

    if (dd_prim_pair_count_local > 0)
    {
        omptimers[thread_id].start("  J block DD");

        const uint32_t offsetJ = 5 * max_prim_pair_count_local;

        // set up thread blocks for J

        const dim3 threads_per_block = dim3(TILE_DIM_SMALL, TILE_DIM_LARGE);

        const dim3 num_blocks = dim3((dd_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (DD|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            const uint32_t offsetD = 0 * max_prim_pair_count;
            omptimers[thread_id].start("    J block DDSS");

            gpu::computeCoulombFockDDSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block DDSS");
        }

        // J: (DD|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 1 * max_prim_pair_count;
            omptimers[thread_id].start("    J block DDSP");

            gpu::computeCoulombFockDDSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block DDSP");
        }

        // J: (DD|SD)
        //     **

        if (sd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 2 * max_prim_pair_count;
            omptimers[thread_id].start("    J block DDSD");

            gpu::computeCoulombFockDDSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_s_prim_info,
                               static_cast<uint32_t>(s_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block DDSD");
        }

        // J: (DD|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            const uint32_t offsetD = 3 * max_prim_pair_count;
            omptimers[thread_id].start("    J block DDPP");

            gpu::computeCoulombFockDDPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block DDPP");
        }

        // J: (DD|PD)
        //     **

        if (pd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 4 * max_prim_pair_count;
            omptimers[thread_id].start("    J block DDPD");

            gpu::computeCoulombFockDDPD0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J +offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD7<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD8<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDPD9<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_p_prim_info,
                               static_cast<uint32_t>(p_prim_count),
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block DDPD");
        }

        // J: (DD|DD)
        //     **

        if (dd_prim_pair_count > 0)
        {
            const uint32_t offsetD = 5 * max_prim_pair_count;
            omptimers[thread_id].start("    J block DDDD");

            gpuSafe(gpuMemcpy(d_mat_D, dd_mat_D.data(), dd_prim_pair_count * sizeof(double), gpuMemcpyHostToDevice));

            dim3 dd_threads_per_block (TILE_DIM_SMALL, TILE_DIM_LARGE);

            dim3 dd_num_blocks ((dd_prim_pair_count_local + dd_threads_per_block.x - 1) / dd_threads_per_block.x, 1);

            gpu::computeCoulombFockDDDD0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            gpu::computeCoulombFockDDDD1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD7<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD8<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD9<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD10<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD11<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD12<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD13<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD14<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD15<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD16<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD17<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD18<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD19<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD20<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD21<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD22<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD23<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD24<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD25<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD26<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD27<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD28<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

        gpu::computeCoulombFockDDDD29<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                               d_mat_J + offsetJ,
                               d_d_prim_info,
                               static_cast<uint32_t>(d_prim_count),
                               d_mat_D + offsetD,
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

            omptimers[thread_id].stop("    J block DDDD");
        }
        eventListCoulomb[5].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  J block DD");
    }

    }  // end of compute J

    coulomb_timer.stop();

    auto coulomb_elapsed_time = coulomb_timer.getElapsedTime();

    screening.setCoulombTime(gpu_id, coulomb_elapsed_time);

    omptimers[thread_id].stop("J computation");

    CTimer exchange_timer;

    exchange_timer.start();

    // compute K

    if (std::fabs(frac_exact_exchange) > 1.0e-13)
    {

    // Zero all blocks
    {
        const dim3 threads_per_block(TILE_DIM * TILE_DIM);

        const dim3 num_blocks(((max_pair_inds_count * numCalculationBlocksExchange) + threads_per_block.x - 1) / threads_per_block.x);

        gpu::zeroData<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(d_mat_K, static_cast<uint32_t>(max_pair_inds_count) * numCalculationBlocksExchange);
    }

    // K: S-S block

    if (pair_inds_count_for_K_ss > 0)
    {
        omptimers[thread_id].start("  K block SS");

        const uint32_t offset = 0 * max_pair_inds_count;

        // set up thread blocks for K

        dim3 threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        dim3 num_blocks = dim3(pair_inds_count_for_K_ss, 1);

        // K: (SS|SS)
        //     *  *

        gpu::computeExchangeFockSSSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SS|SP)
        //     *  *

        gpu::computeExchangeFockSSSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|SS)
        //     *  *

        gpu::computeExchangeFockSPSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|SP)
        //     *  *

        gpu::computeExchangeFockSPSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SS|SD)
        //     *  *

        gpu::computeExchangeFockSSSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SD|SS)
        //     *  *

        gpu::computeExchangeFockSDSS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|SD)
        //     *  *

        gpu::computeExchangeFockSPSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SD|SP)
        //     *  *

        gpu::computeExchangeFockSDSP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SD|SD)
        //     *  *

        gpu::computeExchangeFockSDSD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        eventListExchange[0].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  K block SS");

    }

    if (pair_inds_count_for_K_sp > 0)
    {
        omptimers[thread_id].start("  K block SP");

        const uint32_t offset = 1 * max_pair_inds_count;

        // set up thread blocks for K

        dim3 threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        dim3 num_blocks = dim3(pair_inds_count_for_K_sp, 1);

        // K: (SS|PS)
        //     *  *

        gpu::computeExchangeFockSSPS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SS|PP)
        //     *  *

        gpu::computeExchangeFockSSPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|PS)
        //     *  *

        gpu::computeExchangeFockSPPS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|PP)
        //     *  *

        gpu::computeExchangeFockSPPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SS|PD)
        //     *  *

        gpu::computeExchangeFockSSPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SD|PS)
        //     *  *

        gpu::computeExchangeFockSDPS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|PD)
        //     *  *

        gpu::computeExchangeFockSPPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SD|PP)
        //     *  *

        gpu::computeExchangeFockSDPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SD|PD)
        //     *  *

        gpu::computeExchangeFockSDPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        eventListExchange[1].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  K block SP");
    }

    if (pair_inds_count_for_K_pp > 0)
    {
        omptimers[thread_id].start("  K block PP");
        const uint32_t offset = 2 * max_pair_inds_count;

        // set up thread blocks for K

        dim3 threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        dim3 num_blocks = dim3(pair_inds_count_for_K_pp, 1);

        // K: (PS|PS)
        //     *  *

        omptimers[thread_id].start("    K block PSPS");

        gpu::computeExchangeFockPSPS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PSPS");

        // K: (PS|PP)
        //     *  *

        omptimers[thread_id].start("    K block PSPP");

        gpu::computeExchangeFockPSPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PSPP");

        // K: (PP|PS)
        //     *  *

        omptimers[thread_id].start("    K block PPPS");

        gpu::computeExchangeFockPPPS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PPPS");

        // K: (PP|PP)
        //     *  *

        omptimers[thread_id].start("    K block PPPP");

        gpu::computeExchangeFockPPPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PPPP");

        // K: (PS|PD)
        //     *  *

        omptimers[thread_id].start("    K block PSPD");

        gpu::computeExchangeFockPSPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PSPD");

        // K: (PD|PS)
        //     *  *

        omptimers[thread_id].start("    K block PDPS");

        gpu::computeExchangeFockPDPS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PDPS");

        // K: (PP|PD)
        //     *  *

        omptimers[thread_id].start("    K block PPPD");

        gpu::computeExchangeFockPPPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PPPD");

        // K: (PD|PP)
        //     *  *

        omptimers[thread_id].start("    K block PDPP");

        gpu::computeExchangeFockPDPP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PDPP");

        // K: (PD|PD)
        //     *  *

        omptimers[thread_id].start("    K block PDPD");

        gpu::computeExchangeFockPDPD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PDPD");
        eventListExchange[2].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  K block PP");
    }

    if (pair_inds_count_for_K_sd > 0)
    {
        omptimers[thread_id].start("  K block SD");

        // set up thread blocks for K

        dim3 threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        dim3 num_blocks = dim3(pair_inds_count_for_K_sd, 1);

        const uint32_t offset = 3 * max_pair_inds_count;

        // K: (SS|DS)
        //     *  *

        gpu::computeExchangeFockSSDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SS|DP)
        //     *  *

        gpu::computeExchangeFockSSDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|DS)
        //     *  *

        gpu::computeExchangeFockSPDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (SP|DP)
        //     *  *

        omptimers[thread_id].start("    K block SPDP");

        gpu::computeExchangeFockSPDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block SPDP");

        // K: (SS|DD)
        //     *  *

        omptimers[thread_id].start("    K block SSDD");

        gpu::computeExchangeFockSSDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block SSDD");

        // K: (SD|DS)
        //     *  *

        omptimers[thread_id].start("    K block SDDS");

        gpu::computeExchangeFockSDDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block SDDS");

        // K: (SP|DD)
        //     *  *

        omptimers[thread_id].start("    K block SPDD");

        gpu::computeExchangeFockSPDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block SPDD");

        // K: (SD|DP)
        //     *  *

        omptimers[thread_id].start("    K block SDDP");

        gpu::computeExchangeFockSDDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block SDDP");

        // K: (SD|DD)
        //     *  *

        omptimers[thread_id].start("    K block SDDD");

        gpu::computeExchangeFockSDDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block SDDD");
        eventListExchange[3].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  K block SD");
    }

    if (pair_inds_count_for_K_pd > 0)
    {
        omptimers[thread_id].start("  K block PD");

        // set up thread blocks for K

        dim3 threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        dim3 num_blocks = dim3(pair_inds_count_for_K_pd, 1);

        const uint32_t offset = 4 * max_pair_inds_count;

        // K: (PS|DS)
        //     *  *

        omptimers[thread_id].start("    K block PSDS");

        gpu::computeExchangeFockPSDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PSDS");

        // K: (PS|DP)
        //     *  *

        omptimers[thread_id].start("    K block PSDP");

        gpu::computeExchangeFockPSDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PSDP");

        // K: (PP|DS)
        //     *  *

        omptimers[thread_id].start("    K block PPDS");

        gpu::computeExchangeFockPPDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PPDS");

        // K: (PS|DD)
        //     *  *

        omptimers[thread_id].start("    K block PSDD");

        gpu::computeExchangeFockPSDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PSDD");

        // K: (PD|DS)
        //     *  *

        omptimers[thread_id].start("    K block PDDS");

        gpu::computeExchangeFockPDDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PDDS");

        // K: (PP|DP)
        //     *  *

        omptimers[thread_id].start("    K block PPDP");

        gpu::computeExchangeFockPPDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PPDP");

        // K: (PP|DD)
        //     *  *

        omptimers[thread_id].start("    K block PPDD");

        gpu::computeExchangeFockPPDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PPDD");

        // K: (PD|DP)
        //     *  *

        omptimers[thread_id].start("    K block PDDP");

        gpu::computeExchangeFockPDDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PDDP");

        // K: (PD|DD)
        //     *  *

        omptimers[thread_id].start("    K block PDDD");

        gpu::computeExchangeFockPDDD0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockPDDD7<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block PDDD");
        eventListExchange[4].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  K block PD");
    }

    if (pair_inds_count_for_K_dd > 0)
    {
        omptimers[thread_id].start("  K block DD");

        const uint32_t offset = 5 * max_pair_inds_count;

        // set up thread blocks for K

        dim3 threads_per_block = dim3(TILE_DIM_X_K, TILE_DIM_Y_K);

        dim3 num_blocks = dim3(pair_inds_count_for_K_dd, 1);

        // K: (DS|DS)
        //     *  *

        gpu::computeExchangeFockDSDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (DS|DP)
        //     *  *

        gpu::computeExchangeFockDSDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (DP|DS)
        //     *  *

        gpu::computeExchangeFockDPDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        // K: (DS|DD)
        //     *  *

        omptimers[thread_id].start("    K block DSDD");

        gpu::computeExchangeFockDSDD<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block DSDD");

        // K: (DD|DS)
        //     *  *

        omptimers[thread_id].start("    K block DDDS");

        gpu::computeExchangeFockDDDS<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block DDDS");

        // K: (DP|DP)
        //     *  *

        omptimers[thread_id].start("    K block DPDP");

        gpu::computeExchangeFockDPDP<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block DPDP");

        // K: (DP|DD)
        //     *  *

        omptimers[thread_id].start("    K block DPDD");

        gpu::computeExchangeFockDPDD0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDPDD1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDPDD2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDPDD3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDPDD4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDPDD5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDPDD6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block DPDD");

        // K: (DD|DP)
        //     *  *

        omptimers[thread_id].start("    K block DDDP");

        gpu::computeExchangeFockDDDP0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDP1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDP2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDP3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDP4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDP5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDP6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block DDDP");

        // K: (DD|DD)
        //     *  *

        omptimers[thread_id].start("    K block DDDD");

        gpu::computeExchangeFockDDDD0<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD1<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD2<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD3<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD4<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD5<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD6<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD7<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD8<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD9<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD10<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD11<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD12<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD13<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD14<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD15<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD16<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD17<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        gpu::computeExchangeFockDDDD18<<<num_blocks, threads_per_block, 0, streamList[0].stream>>>(
                           d_mat_K + offset,
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

        omptimers[thread_id].stop("    K block DDDD");
        eventListExchange[5].markStreamEvent(streamList[0].stream);
        omptimers[thread_id].stop("  K block DD");
    }

    } // end of K kernel launches

    if (std::fabs(prefac_coulomb) > 1.0e-13)
    {

    if (ss_prim_pair_count_local > 0)
    {
        eventListCoulomb[0].waitForCompletion();
        const uint32_t offset = 0 * max_prim_pair_count_local;
        gpuSafe(gpuMemcpyAsync(mat_J.data(), d_mat_J + offset, ss_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * prefac_coulomb;

            if (i != j) mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * prefac_coulomb;
        }
    }

    if (sp_prim_pair_count_local > 0)
    {
        eventListCoulomb[1].waitForCompletion();
        const uint32_t offset = 1 * max_prim_pair_count_local;
        gpuSafe(gpuMemcpyAsync(mat_J.data(), d_mat_J + offset, sp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pp_prim_pair_count_local > 0)
    {
        eventListCoulomb[2].waitForCompletion();
        const uint32_t offset = 2 * max_prim_pair_count_local;
        gpuSafe(gpuMemcpyAsync(mat_J.data(), d_mat_J + offset, pp_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (sd_prim_pair_count_local > 0)
    {
        eventListCoulomb[3].waitForCompletion();
        const uint32_t offset = 3 * max_prim_pair_count_local;
        gpuSafe(gpuMemcpyAsync(mat_J.data(), d_mat_J + offset, sd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pd_prim_pair_count_local > 0)
    {
        eventListCoulomb[4].waitForCompletion();
        const uint32_t offset = 4 * max_prim_pair_count_local;
        gpuSafe(gpuMemcpyAsync(mat_J.data(), d_mat_J + offset, pd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (dd_prim_pair_count_local > 0)
    {
        eventListCoulomb[5].waitForCompletion();
        const uint32_t offset = 5 * max_prim_pair_count_local;
        gpuSafe(gpuMemcpyAsync(mat_J.data(), d_mat_J + offset, dd_prim_pair_count_local * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    }  // end of compute J post processing

    if (std::fabs(frac_exact_exchange) > 1.0e-13)
    {

    if (pair_inds_count_for_K_ss > 0)
    {
        eventListExchange[0].waitForCompletion();
        const uint32_t offset = 0 * max_pair_inds_count;

        gpuSafe(gpuMemcpyAsync(mat_K.data(), d_mat_K + offset, pair_inds_count_for_K_ss * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pair_inds_count_for_K_sp > 0)
    {
        eventListExchange[1].waitForCompletion();
        const uint32_t offset = 1 * max_pair_inds_count;

        gpuSafe(gpuMemcpyAsync(mat_K.data(), d_mat_K + offset, pair_inds_count_for_K_sp * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pair_inds_count_for_K_pp > 0)
    {
        eventListExchange[2].waitForCompletion();
        const uint32_t offset = 2 * max_pair_inds_count;

        gpuSafe(gpuMemcpyAsync(mat_K.data(), d_mat_K + offset, pair_inds_count_for_K_pp * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pair_inds_count_for_K_sd > 0)
    {
        eventListExchange[3].waitForCompletion();
        const uint32_t offset = 3 * max_pair_inds_count;

        gpuSafe(gpuMemcpyAsync(mat_K.data(), d_mat_K + offset, pair_inds_count_for_K_sd * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pair_inds_count_for_K_pd > 0)
    {
        eventListExchange[4].waitForCompletion();
        const uint32_t offset = 4 * max_pair_inds_count;

        gpuSafe(gpuMemcpyAsync(mat_K.data(), d_mat_K + offset, pair_inds_count_for_K_pd * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    if (pair_inds_count_for_K_dd > 0)
    {
        eventListExchange[5].waitForCompletion();
        const uint32_t offset = 5 * max_pair_inds_count;

        gpuSafe(gpuMemcpyAsync(mat_K.data(), d_mat_K + offset, pair_inds_count_for_K_dd * sizeof(double), gpuMemcpyDeviceToHost, streamList[2].stream));
        copyEvent.markStreamEventAndWait(streamList[2].stream);

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
    }

    }  // end of compute K post processing

    for (uint32_t i = 0; i < 6; i++)
    {
        eventListExchange[i].destroyEvent();
        eventListCoulomb[i].destroyEvent();
    }
    copyEvent.destroyEvent();

    gpuSafe(gpuDeviceSynchronize());

    exchange_timer.stop();

    auto exchange_elapsed_time = exchange_timer.getElapsedTime();

    screening.setExchangeTime(gpu_id, exchange_elapsed_time);

    streamList[0].destroyStream();
    streamList[1].destroyStream();
    streamList[2].destroyStream();
    omptimers[thread_id].stop("K computation");

    omptimers[thread_id].start("J finalize");

    gpuSafe(gpuFree(d_data_mat_D_J));
    gpuSafe(gpuFree(d_data_mat_Q));
    gpuSafe(gpuFree(d_data_first_second_inds));
    gpuSafe(gpuFree(d_data_pair_data));
    gpuSafe(gpuFree(d_data_mat_Q_local));
    gpuSafe(gpuFree(d_data_first_second_inds_local));
    gpuSafe(gpuFree(d_data_pair_data_local));

    omptimers[thread_id].stop("J finalize");

    omptimers[thread_id].start("K finalize");

    gpuSafe(gpuFree(d_data_boys_func));

    gpuSafe(gpuFree(d_data_spd_prim_info));
    gpuSafe(gpuFree(d_data_spd_prim_aoinds));

    gpuSafe(gpuFree(d_mat_K));
    gpuSafe(gpuFree(d_data_pair_inds_for_K));
    gpuSafe(gpuFree(d_mat_D_full_AO));
    gpuSafe(gpuFree(d_data_Q_K));
    gpuSafe(gpuFree(d_data_D_inds_K));
    gpuSafe(gpuFree(d_data_pair_counts_displs_K));
    gpuSafe(gpuFree(d_data_pair_data_K));

    omptimers[thread_id].stop("K finalize");
    }
    }

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

    // std::cout << "Timing of Fock build" << std::endl;
    // std::cout << "--------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;
    // std::cout << "GPU timing" << std::endl;
    // std::cout << "----------" << std::endl;
    // for (int thread_id = 0; thread_id < nthreads; thread_id++)
    // {
    //     std::cout << "GPU " << thread_id << std::endl;
    //     std::cout << omptimers[thread_id].getSummary() << std::endl;
    // }

    return mat_Fock_sum;
}

}  // namespace gpu
