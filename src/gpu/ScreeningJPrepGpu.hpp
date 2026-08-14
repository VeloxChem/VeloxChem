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

#ifndef ScreeningJPrepGpu_hpp
#define ScreeningJPrepGpu_hpp

#include <algorithm>
#include <cstdint>
#include <vector>

#include "GpuRuntime.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "ScreeningData.hpp"

namespace gpujprep {

struct JPrepDeviceData
{
    int64_t ss_prim_pair_count{0};
    int64_t sp_prim_pair_count{0};
    int64_t sd_prim_pair_count{0};
    int64_t pp_prim_pair_count{0};
    int64_t pd_prim_pair_count{0};
    int64_t dd_prim_pair_count{0};

    int64_t ss_prim_pair_count_local{0};
    int64_t sp_prim_pair_count_local{0};
    int64_t sd_prim_pair_count_local{0};
    int64_t pp_prim_pair_count_local{0};
    int64_t pd_prim_pair_count_local{0};
    int64_t dd_prim_pair_count_local{0};

    int64_t max_prim_pair_count{0};
    int64_t max_prim_pair_count_local{0};

    double* d_data_mat_D_J{nullptr};
    double* d_mat_D{nullptr};
    double* d_mat_J{nullptr};

    double* d_data_mat_Q{nullptr};
    double* d_ss_mat_Q{nullptr};
    double* d_sp_mat_Q{nullptr};
    double* d_sd_mat_Q{nullptr};
    double* d_pp_mat_Q{nullptr};
    double* d_pd_mat_Q{nullptr};
    double* d_dd_mat_Q{nullptr};

    uint32_t* d_data_first_second_inds{nullptr};
    uint32_t* d_ss_first_inds{nullptr};
    uint32_t* d_ss_second_inds{nullptr};
    uint32_t* d_sp_first_inds{nullptr};
    uint32_t* d_sp_second_inds{nullptr};
    uint32_t* d_sd_first_inds{nullptr};
    uint32_t* d_sd_second_inds{nullptr};
    uint32_t* d_pp_first_inds{nullptr};
    uint32_t* d_pp_second_inds{nullptr};
    uint32_t* d_pd_first_inds{nullptr};
    uint32_t* d_pd_second_inds{nullptr};
    uint32_t* d_dd_first_inds{nullptr};
    uint32_t* d_dd_second_inds{nullptr};

    double* d_data_pair_data{nullptr};
    double* d_ss_pair_data{nullptr};
    double* d_sp_pair_data{nullptr};
    double* d_sd_pair_data{nullptr};
    double* d_pp_pair_data{nullptr};
    double* d_pd_pair_data{nullptr};
    double* d_dd_pair_data{nullptr};

    double* d_data_mat_Q_local{nullptr};
    double* d_ss_mat_Q_local{nullptr};
    double* d_sp_mat_Q_local{nullptr};
    double* d_sd_mat_Q_local{nullptr};
    double* d_pp_mat_Q_local{nullptr};
    double* d_pd_mat_Q_local{nullptr};
    double* d_dd_mat_Q_local{nullptr};

    uint32_t* d_data_first_second_inds_local{nullptr};
    uint32_t* d_ss_first_inds_local{nullptr};
    uint32_t* d_ss_second_inds_local{nullptr};
    uint32_t* d_sp_first_inds_local{nullptr};
    uint32_t* d_sp_second_inds_local{nullptr};
    uint32_t* d_sd_first_inds_local{nullptr};
    uint32_t* d_sd_second_inds_local{nullptr};
    uint32_t* d_pp_first_inds_local{nullptr};
    uint32_t* d_pp_second_inds_local{nullptr};
    uint32_t* d_pd_first_inds_local{nullptr};
    uint32_t* d_pd_second_inds_local{nullptr};
    uint32_t* d_dd_first_inds_local{nullptr};
    uint32_t* d_dd_second_inds_local{nullptr};

    double* d_data_pair_data_local{nullptr};
    double* d_ss_pair_data_local{nullptr};
    double* d_sp_pair_data_local{nullptr};
    double* d_sd_pair_data_local{nullptr};
    double* d_pp_pair_data_local{nullptr};
    double* d_pd_pair_data_local{nullptr};
    double* d_dd_pair_data_local{nullptr};
};

inline auto uploadJPrepDeviceData(const CScreeningData& screening,
                                  const int64_t         gpu_id,
                                  const int             rank,
                                  const int64_t         total_num_gpus_per_compute_node) -> JPrepDeviceData
{
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

    const auto& ss_pair_data = screening.get_ss_pair_data();
    const auto& sp_pair_data = screening.get_sp_pair_data();
    const auto& sd_pair_data = screening.get_sd_pair_data();
    const auto& pp_pair_data = screening.get_pp_pair_data();
    const auto& pd_pair_data = screening.get_pd_pair_data();
    const auto& dd_pair_data = screening.get_dd_pair_data();

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

    JPrepDeviceData data;

    data.ss_prim_pair_count = static_cast<int64_t>(ss_first_inds.size());
    data.sp_prim_pair_count = static_cast<int64_t>(sp_first_inds.size());
    data.sd_prim_pair_count = static_cast<int64_t>(sd_first_inds.size());
    data.pp_prim_pair_count = static_cast<int64_t>(pp_first_inds.size());
    data.pd_prim_pair_count = static_cast<int64_t>(pd_first_inds.size());
    data.dd_prim_pair_count = static_cast<int64_t>(dd_first_inds.size());

    data.ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    data.sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    data.sd_prim_pair_count_local = static_cast<int64_t>(sd_first_inds_local.size());
    data.pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());
    data.pd_prim_pair_count_local = static_cast<int64_t>(pd_first_inds_local.size());
    data.dd_prim_pair_count_local = static_cast<int64_t>(dd_first_inds_local.size());

    data.max_prim_pair_count = std::max({data.ss_prim_pair_count,
                                         data.sp_prim_pair_count,
                                         data.sd_prim_pair_count,
                                         data.pp_prim_pair_count,
                                         data.pd_prim_pair_count,
                                         data.dd_prim_pair_count});

    data.max_prim_pair_count_local = std::max({data.ss_prim_pair_count_local,
                                               data.sp_prim_pair_count_local,
                                               data.sd_prim_pair_count_local,
                                               data.pp_prim_pair_count_local,
                                               data.pd_prim_pair_count_local,
                                               data.dd_prim_pair_count_local});

    const auto gpu_rank = gpu_id + rank * total_num_gpus_per_compute_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    gpuSafe(gpuMalloc(&data.d_data_mat_D_J, (data.max_prim_pair_count + data.max_prim_pair_count_local) * sizeof(double)));

    data.d_mat_D = data.d_data_mat_D_J;
    data.d_mat_J = data.d_mat_D + data.max_prim_pair_count;

    gpuSafe(gpuMalloc(&data.d_data_mat_Q,
                      (data.ss_prim_pair_count + data.sp_prim_pair_count + data.sd_prim_pair_count + data.pp_prim_pair_count +
                       data.pd_prim_pair_count + data.dd_prim_pair_count) *
                          sizeof(double)));

    data.d_ss_mat_Q = data.d_data_mat_Q;
    data.d_sp_mat_Q = data.d_ss_mat_Q + data.ss_prim_pair_count;
    data.d_sd_mat_Q = data.d_sp_mat_Q + data.sp_prim_pair_count;
    data.d_pp_mat_Q = data.d_sd_mat_Q + data.sd_prim_pair_count;
    data.d_pd_mat_Q = data.d_pp_mat_Q + data.pp_prim_pair_count;
    data.d_dd_mat_Q = data.d_pd_mat_Q + data.pd_prim_pair_count;

    gpuSafe(gpuMalloc(&data.d_data_first_second_inds,
                      (data.ss_prim_pair_count + data.ss_prim_pair_count + data.sp_prim_pair_count + data.sp_prim_pair_count +
                       data.sd_prim_pair_count + data.sd_prim_pair_count + data.pp_prim_pair_count + data.pp_prim_pair_count +
                       data.pd_prim_pair_count + data.pd_prim_pair_count + data.dd_prim_pair_count + data.dd_prim_pair_count) *
                          sizeof(uint32_t)));

    data.d_ss_first_inds  = data.d_data_first_second_inds;
    data.d_ss_second_inds = data.d_ss_first_inds + data.ss_prim_pair_count;
    data.d_sp_first_inds  = data.d_ss_second_inds + data.ss_prim_pair_count;
    data.d_sp_second_inds = data.d_sp_first_inds + data.sp_prim_pair_count;
    data.d_sd_first_inds  = data.d_sp_second_inds + data.sp_prim_pair_count;
    data.d_sd_second_inds = data.d_sd_first_inds + data.sd_prim_pair_count;
    data.d_pp_first_inds  = data.d_sd_second_inds + data.sd_prim_pair_count;
    data.d_pp_second_inds = data.d_pp_first_inds + data.pp_prim_pair_count;
    data.d_pd_first_inds  = data.d_pp_second_inds + data.pp_prim_pair_count;
    data.d_pd_second_inds = data.d_pd_first_inds + data.pd_prim_pair_count;
    data.d_dd_first_inds  = data.d_pd_second_inds + data.pd_prim_pair_count;
    data.d_dd_second_inds = data.d_dd_first_inds + data.dd_prim_pair_count;

    gpuSafe(gpuMalloc(&data.d_data_pair_data,
                      (ss_pair_data.size() + sp_pair_data.size() + sd_pair_data.size() + pp_pair_data.size() + pd_pair_data.size() +
                       dd_pair_data.size()) *
                          sizeof(double)));

    data.d_ss_pair_data = data.d_data_pair_data;
    data.d_sp_pair_data = data.d_ss_pair_data + ss_pair_data.size();
    data.d_sd_pair_data = data.d_sp_pair_data + sp_pair_data.size();
    data.d_pp_pair_data = data.d_sd_pair_data + sd_pair_data.size();
    data.d_pd_pair_data = data.d_pp_pair_data + pp_pair_data.size();
    data.d_dd_pair_data = data.d_pd_pair_data + pd_pair_data.size();

    gpuSafe(gpuMalloc(&data.d_data_mat_Q_local,
                      (data.ss_prim_pair_count_local + data.sp_prim_pair_count_local + data.sd_prim_pair_count_local +
                       data.pp_prim_pair_count_local + data.pd_prim_pair_count_local + data.dd_prim_pair_count_local) *
                          sizeof(double)));

    data.d_ss_mat_Q_local = data.d_data_mat_Q_local;
    data.d_sp_mat_Q_local = data.d_ss_mat_Q_local + data.ss_prim_pair_count_local;
    data.d_sd_mat_Q_local = data.d_sp_mat_Q_local + data.sp_prim_pair_count_local;
    data.d_pp_mat_Q_local = data.d_sd_mat_Q_local + data.sd_prim_pair_count_local;
    data.d_pd_mat_Q_local = data.d_pp_mat_Q_local + data.pp_prim_pair_count_local;
    data.d_dd_mat_Q_local = data.d_pd_mat_Q_local + data.pd_prim_pair_count_local;

    gpuSafe(gpuMalloc(&data.d_data_first_second_inds_local,
                      (data.ss_prim_pair_count_local + data.ss_prim_pair_count_local + data.sp_prim_pair_count_local +
                       data.sp_prim_pair_count_local + data.sd_prim_pair_count_local + data.sd_prim_pair_count_local +
                       data.pp_prim_pair_count_local + data.pp_prim_pair_count_local + data.pd_prim_pair_count_local +
                       data.pd_prim_pair_count_local + data.dd_prim_pair_count_local + data.dd_prim_pair_count_local) *
                          sizeof(uint32_t)));

    data.d_ss_first_inds_local  = data.d_data_first_second_inds_local;
    data.d_ss_second_inds_local = data.d_ss_first_inds_local + data.ss_prim_pair_count_local;
    data.d_sp_first_inds_local  = data.d_ss_second_inds_local + data.ss_prim_pair_count_local;
    data.d_sp_second_inds_local = data.d_sp_first_inds_local + data.sp_prim_pair_count_local;
    data.d_sd_first_inds_local  = data.d_sp_second_inds_local + data.sp_prim_pair_count_local;
    data.d_sd_second_inds_local = data.d_sd_first_inds_local + data.sd_prim_pair_count_local;
    data.d_pp_first_inds_local  = data.d_sd_second_inds_local + data.sd_prim_pair_count_local;
    data.d_pp_second_inds_local = data.d_pp_first_inds_local + data.pp_prim_pair_count_local;
    data.d_pd_first_inds_local  = data.d_pp_second_inds_local + data.pp_prim_pair_count_local;
    data.d_pd_second_inds_local = data.d_pd_first_inds_local + data.pd_prim_pair_count_local;
    data.d_dd_first_inds_local  = data.d_pd_second_inds_local + data.pd_prim_pair_count_local;
    data.d_dd_second_inds_local = data.d_dd_first_inds_local + data.dd_prim_pair_count_local;

    gpuSafe(gpuMalloc(&data.d_data_pair_data_local,
                      (ss_pair_data_local.size() + sp_pair_data_local.size() + sd_pair_data_local.size() + pp_pair_data_local.size() +
                       pd_pair_data_local.size() + dd_pair_data_local.size()) *
                          sizeof(double)));

    data.d_ss_pair_data_local = data.d_data_pair_data_local;
    data.d_sp_pair_data_local = data.d_ss_pair_data_local + ss_pair_data_local.size();
    data.d_sd_pair_data_local = data.d_sp_pair_data_local + sp_pair_data_local.size();
    data.d_pp_pair_data_local = data.d_sd_pair_data_local + sd_pair_data_local.size();
    data.d_pd_pair_data_local = data.d_pp_pair_data_local + pp_pair_data_local.size();
    data.d_dd_pair_data_local = data.d_pd_pair_data_local + pd_pair_data_local.size();

    if (!ss_mat_Q.empty())
    {
        gpuSafe(gpuMemcpy(data.d_ss_mat_Q, ss_mat_Q.data(), ss_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_mat_Q, sp_mat_Q.data(), sp_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_mat_Q, sd_mat_Q.data(), sd_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_mat_Q, pp_mat_Q.data(), pp_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_mat_Q, pd_mat_Q.data(), pd_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_mat_Q, dd_mat_Q.data(), dd_mat_Q.size() * sizeof(double), gpuMemcpyHostToDevice));
    }

    if (!ss_first_inds.empty())
    {
        gpuSafe(gpuMemcpy(data.d_ss_first_inds, ss_first_inds.data(), ss_first_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_ss_second_inds, ss_second_inds.data(), ss_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_first_inds, sp_first_inds.data(), sp_first_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_second_inds, sp_second_inds.data(), sp_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_first_inds, sd_first_inds.data(), sd_first_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_second_inds, sd_second_inds.data(), sd_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_first_inds, pp_first_inds.data(), pp_first_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_second_inds, pp_second_inds.data(), pp_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_first_inds, pd_first_inds.data(), pd_first_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_second_inds, pd_second_inds.data(), pd_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_first_inds, dd_first_inds.data(), dd_first_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_second_inds, dd_second_inds.data(), dd_second_inds.size() * sizeof(uint32_t), gpuMemcpyHostToDevice));
    }

    if (!ss_pair_data.empty())
    {
        gpuSafe(gpuMemcpy(data.d_ss_pair_data, ss_pair_data.data(), ss_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_pair_data, sp_pair_data.data(), sp_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_pair_data, sd_pair_data.data(), sd_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_pair_data, pp_pair_data.data(), pp_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_pair_data, pd_pair_data.data(), pd_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_pair_data, dd_pair_data.data(), dd_pair_data.size() * sizeof(double), gpuMemcpyHostToDevice));
    }

    if (!ss_mat_Q_local.empty())
    {
        gpuSafe(gpuMemcpy(data.d_ss_mat_Q_local, ss_mat_Q_local.data(), ss_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_mat_Q_local, sp_mat_Q_local.data(), sp_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_mat_Q_local, sd_mat_Q_local.data(), sd_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_mat_Q_local, pp_mat_Q_local.data(), pp_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_mat_Q_local, pd_mat_Q_local.data(), pd_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_mat_Q_local, dd_mat_Q_local.data(), dd_mat_Q_local.size() * sizeof(double), gpuMemcpyHostToDevice));
    }

    if (!ss_first_inds_local.empty())
    {
        gpuSafe(gpuMemcpy(data.d_ss_first_inds_local, ss_first_inds_local.data(), ss_first_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_ss_second_inds_local, ss_second_inds_local.data(), ss_second_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_first_inds_local, sp_first_inds_local.data(), sp_first_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_second_inds_local, sp_second_inds_local.data(), sp_second_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_first_inds_local, sd_first_inds_local.data(), sd_first_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_second_inds_local, sd_second_inds_local.data(), sd_second_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_first_inds_local, pp_first_inds_local.data(), pp_first_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_second_inds_local, pp_second_inds_local.data(), pp_second_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_first_inds_local, pd_first_inds_local.data(), pd_first_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_second_inds_local, pd_second_inds_local.data(), pd_second_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_first_inds_local, dd_first_inds_local.data(), dd_first_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_second_inds_local, dd_second_inds_local.data(), dd_second_inds_local.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
    }

    if (!ss_pair_data_local.empty())
    {
        gpuSafe(gpuMemcpy(data.d_ss_pair_data_local, ss_pair_data_local.data(), ss_pair_data_local.size() * sizeof(double),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sp_pair_data_local, sp_pair_data_local.data(), sp_pair_data_local.size() * sizeof(double),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_sd_pair_data_local, sd_pair_data_local.data(), sd_pair_data_local.size() * sizeof(double),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pp_pair_data_local, pp_pair_data_local.data(), pp_pair_data_local.size() * sizeof(double),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_pd_pair_data_local, pd_pair_data_local.data(), pd_pair_data_local.size() * sizeof(double),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(data.d_dd_pair_data_local, dd_pair_data_local.data(), dd_pair_data_local.size() * sizeof(double),
                          gpuMemcpyHostToDevice));
    }

    return data;
}

inline auto uploadJPrepForAllDevices(const CScreeningData& screening,
                                     const int64_t         num_gpus_per_node,
                                     const int             rank,
                                     const int64_t         total_num_gpus_per_compute_node) -> std::vector<JPrepDeviceData>
{
    std::vector<JPrepDeviceData> device_data(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        device_data[gpu_id] = uploadJPrepDeviceData(screening, gpu_id, rank, total_num_gpus_per_compute_node);
    }

    return device_data;
}

inline auto freeJPrepDeviceData(const int64_t gpu_id,
                                const int     rank,
                                const int64_t total_num_gpus_per_compute_node,
                                JPrepDeviceData& data) -> void
{
    if (data.d_data_mat_D_J == nullptr) return;

    const auto gpu_rank = gpu_id + rank * total_num_gpus_per_compute_node;

    gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));

    gpuSafe(gpuFree(data.d_data_mat_D_J));
    gpuSafe(gpuFree(data.d_data_mat_Q));
    gpuSafe(gpuFree(data.d_data_first_second_inds));
    gpuSafe(gpuFree(data.d_data_pair_data));
    gpuSafe(gpuFree(data.d_data_mat_Q_local));
    gpuSafe(gpuFree(data.d_data_first_second_inds_local));
    gpuSafe(gpuFree(data.d_data_pair_data_local));

    data = JPrepDeviceData{};
}

inline auto freeJPrepForAllDevices(const int64_t                   num_gpus_per_node,
                                   const int                       rank,
                                   const int64_t                   total_num_gpus_per_compute_node,
                                   std::vector<JPrepDeviceData>& device_data) -> void
{
    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        freeJPrepDeviceData(gpu_id, rank, total_num_gpus_per_compute_node, device_data[gpu_id]);
    }
}

}  // namespace gpujprep

#endif /* ScreeningJPrepGpu_hpp */
