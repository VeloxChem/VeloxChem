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

#ifndef GpuScreeningDeviceStorage_hpp
#define GpuScreeningDeviceStorage_hpp

#include <cstdint>
#include <vector>

#include "BoysFuncTableGpu.hpp"
#include "GpuDeviceBufferStorage.hpp"
#include "ScreeningJPrepGpu.hpp"

class CScreeningData;

struct PrimInfoDeviceViews
{
    double*   d_data_spd_prim_info{nullptr};
    uint32_t* d_data_spd_prim_aoinds{nullptr};

    double*   d_s_prim_info{nullptr};
    double*   d_p_prim_info{nullptr};
    double*   d_d_prim_info{nullptr};

    uint32_t* d_s_prim_aoinds{nullptr};
    uint32_t* d_p_prim_aoinds{nullptr};
    uint32_t* d_d_prim_aoinds{nullptr};
};

/**
 Persistent per-GPU device buffers for geometry-dependent screening data.
 Buffers are reused across Fock builds and structure optimization steps; they
 grow only when the required size increases.
 */
class CGpuScreeningDeviceStorage
{
    struct GpuDeviceSlot
    {
        int64_t gpu_id{0};

        int device_id{0};

        CGpuDeviceBuffer boys_buffer;

        boysfunc::DeviceTables boys_tables;

        gpujprep::JPrepDeviceBufferPools j_prep_pools;

        gpujprep::JPrepDeviceData j_prep;

        CGpuDeviceBuffer prim_info_buffer;

        CGpuDeviceBuffer prim_aoinds_buffer;

        PrimInfoDeviceViews prim_views;
    };

    int     _rank{0};

    int64_t _num_gpus_per_node{0};

    int64_t _total_num_gpus_per_compute_node{0};

    bool _configured{false};

    std::vector<GpuDeviceSlot> _slots;

    auto _slot(const int64_t gpu_id) -> GpuDeviceSlot&;

    auto _slot(const int64_t gpu_id) const -> const GpuDeviceSlot&;

    auto _deviceId(const int64_t gpu_id) const -> int;

   public:
    CGpuScreeningDeviceStorage() = default;

    CGpuScreeningDeviceStorage(const CGpuScreeningDeviceStorage&)            = delete;
    CGpuScreeningDeviceStorage& operator=(const CGpuScreeningDeviceStorage&) = delete;

    ~CGpuScreeningDeviceStorage();

    auto configure(const int rank, const int64_t num_gpus_per_node, const int64_t total_num_gpus_per_compute_node) -> void;

    auto ensureBoysTables() -> void;

    auto ensureJPrep(const CScreeningData& screening) -> void;

    auto ensurePrimInfo(const int64_t                              gpu_id,
                        const std::vector<double>&                 s_prim_info,
                        const std::vector<double>&                 p_prim_info,
                        const std::vector<double>&                 d_prim_info,
                        const std::vector<uint32_t>&               s_prim_aoinds,
                        const std::vector<uint32_t>&               p_prim_aoinds,
                        const std::vector<uint32_t>&               d_prim_aoinds) -> PrimInfoDeviceViews&;

    auto getBoysTables(const int64_t gpu_id) const -> const boysfunc::DeviceTables&;

    auto getJPrep(const int64_t gpu_id) -> gpujprep::JPrepDeviceData&;

    auto getJPrep(const int64_t gpu_id) const -> const gpujprep::JPrepDeviceData&;

    auto getPrimViews(const int64_t gpu_id) const -> const PrimInfoDeviceViews&;

    auto releaseAll() -> void;
};

#endif /* GpuScreeningDeviceStorage_hpp */
