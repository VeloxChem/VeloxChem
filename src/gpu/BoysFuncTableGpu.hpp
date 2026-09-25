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

#ifndef BoysFuncTableGpu_hpp
#define BoysFuncTableGpu_hpp

#include <cstdint>
#include <vector>

#include "BoysFuncTable.hpp"
#include "GpuDeviceBufferStorage.hpp"
#include "GpuRuntime.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"

namespace boysfunc {

struct DeviceTables
{
    double* table{nullptr};
    double* ft{nullptr};
    double* data{nullptr};
};

inline auto ensureBoysFuncTables(CGpuDeviceBuffer& buffer,
                                 const int         device_id,
                                 DeviceTables&     tables) -> void
{
    const auto& boys_func_table = getFullBoysFuncTable();
    const auto& boys_func_ft    = getBoysFuncFactors();
    const auto  table_size      = static_cast<int64_t>(boys_func_table.size());
    const auto  ft_size         = static_cast<int64_t>(boys_func_ft.size());
    const auto  nbytes          = static_cast<size_t>(table_size + ft_size) * sizeof(double);

    if (nbytes == 0)
    {
        return;
    }

    if (buffer.ptr() != nullptr && static_cast<int>(buffer.deviceId()) == device_id && buffer.capacity() >= nbytes)
    {
        tables.data  = static_cast<double*>(buffer.ptr());
        tables.table = tables.data;
        tables.ft    = tables.data + table_size;
        return;
    }

    tables.data  = static_cast<double*>(buffer.ensure(device_id, nbytes));
    tables.table = tables.data;
    tables.ft    = tables.data + table_size;

    gpuSafe(gpuSetDevice(device_id));
    gpuSafe(gpuMemcpy(tables.table, boys_func_table.data(), table_size * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(tables.ft, boys_func_ft.data(), ft_size * sizeof(double), gpuMemcpyHostToDevice));
}

inline auto uploadFullBoysFuncTables(const int64_t num_gpus_per_node,
                                     const int rank,
                                     const int64_t total_num_gpus_per_compute_node) -> std::vector<DeviceTables>
{
    std::vector<DeviceTables> device_tables(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        const auto gpu_rank = gpu_id + rank * num_gpus_per_node;

        CGpuDeviceBuffer buffer;

        ensureBoysFuncTables(buffer, static_cast<int>(gpu_rank % total_num_gpus_per_compute_node), device_tables[gpu_id]);
    }

    return device_tables;
}

inline auto freeFullBoysFuncTables(const int64_t num_gpus_per_node,
                                   const int rank,
                                   const int64_t total_num_gpus_per_compute_node,
                                   std::vector<DeviceTables>& device_tables) -> void
{
    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        if (device_tables[gpu_id].data == nullptr) continue;

        const auto gpu_rank = gpu_id + rank * num_gpus_per_node;

        gpuSafe(gpuSetDevice(gpu_rank % total_num_gpus_per_compute_node));
        gpuSafe(gpuFree(device_tables[gpu_id].data));

        device_tables[gpu_id] = DeviceTables{};
    }
}

}  // namespace boysfunc

#endif /* BoysFuncTableGpu_hpp */
