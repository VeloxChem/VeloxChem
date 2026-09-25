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

#include "GpuScreeningDeviceStorage.hpp"

#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "ScreeningData.hpp"

CGpuScreeningDeviceStorage::~CGpuScreeningDeviceStorage() { releaseAll(); }

auto CGpuScreeningDeviceStorage::_deviceId(const int64_t gpu_id) const -> int
{
    return static_cast<int>((gpu_id + _rank * _num_gpus_per_node) % _total_num_gpus_per_compute_node);
}

auto CGpuScreeningDeviceStorage::_slot(const int64_t gpu_id) -> GpuDeviceSlot&
{
    return _slots.at(static_cast<size_t>(gpu_id));
}

auto CGpuScreeningDeviceStorage::_slot(const int64_t gpu_id) const -> const GpuDeviceSlot&
{
    return _slots.at(static_cast<size_t>(gpu_id));
}

auto CGpuScreeningDeviceStorage::configure(const int     rank,
                                           const int64_t num_gpus_per_node,
                                           const int64_t total_num_gpus_per_compute_node) -> void
{
    if (_configured && _rank == rank && _num_gpus_per_node == num_gpus_per_node &&
        _total_num_gpus_per_compute_node == total_num_gpus_per_compute_node &&
        static_cast<int64_t>(_slots.size()) == num_gpus_per_node)
    {
        return;
    }

    releaseAll();

    _rank                            = rank;
    _num_gpus_per_node               = num_gpus_per_node;
    _total_num_gpus_per_compute_node = total_num_gpus_per_compute_node;
    _slots                           = std::vector<GpuDeviceSlot>(static_cast<size_t>(num_gpus_per_node));
    _configured                      = true;

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        _slots[static_cast<size_t>(gpu_id)].gpu_id    = gpu_id;
        _slots[static_cast<size_t>(gpu_id)].device_id = _deviceId(gpu_id);
    }
}

auto CGpuScreeningDeviceStorage::ensureBoysTables() -> void
{
    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        auto& slot = _slot(gpu_id);

        boysfunc::ensureBoysFuncTables(slot.boys_buffer, slot.device_id, slot.boys_tables);
    }
}

auto CGpuScreeningDeviceStorage::ensureJPrep(const CScreeningData& screening) -> void
{
    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        auto& slot = _slot(gpu_id);

        slot.j_prep = gpujprep::ensureJPrepDeviceData(slot.j_prep_pools, screening, gpu_id, _rank, _total_num_gpus_per_compute_node);
    }
}

auto CGpuScreeningDeviceStorage::ensurePrimInfo(const int64_t                gpu_id,
                                                const std::vector<double>&   s_prim_info,
                                                const std::vector<double>&   p_prim_info,
                                                const std::vector<double>&   d_prim_info,
                                                const std::vector<uint32_t>& s_prim_aoinds,
                                                const std::vector<uint32_t>& p_prim_aoinds,
                                                const std::vector<uint32_t>& d_prim_aoinds) -> PrimInfoDeviceViews&
{
    auto& slot = _slot(gpu_id);

    const auto prim_info_bytes   = static_cast<size_t>(s_prim_info.size() + p_prim_info.size() + d_prim_info.size()) * sizeof(double);
    const auto prim_aoinds_bytes = static_cast<size_t>(s_prim_aoinds.size() + p_prim_aoinds.size() + d_prim_aoinds.size()) * sizeof(uint32_t);

    slot.prim_views.d_data_spd_prim_info = static_cast<double*>(slot.prim_info_buffer.ensure(slot.device_id, prim_info_bytes));
    slot.prim_views.d_data_spd_prim_aoinds =
        static_cast<uint32_t*>(slot.prim_aoinds_buffer.ensure(slot.device_id, prim_aoinds_bytes));

    slot.prim_views.d_s_prim_info = slot.prim_views.d_data_spd_prim_info;
    slot.prim_views.d_p_prim_info = slot.prim_views.d_s_prim_info + s_prim_info.size();
    slot.prim_views.d_d_prim_info = slot.prim_views.d_p_prim_info + p_prim_info.size();

    slot.prim_views.d_s_prim_aoinds = slot.prim_views.d_data_spd_prim_aoinds;
    slot.prim_views.d_p_prim_aoinds = slot.prim_views.d_s_prim_aoinds + s_prim_aoinds.size();
    slot.prim_views.d_d_prim_aoinds = slot.prim_views.d_p_prim_aoinds + p_prim_aoinds.size();

    gpuSafe(gpuSetDevice(slot.device_id));

    if (!s_prim_info.empty())
    {
        gpuSafe(gpuMemcpy(slot.prim_views.d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(slot.prim_views.d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(slot.prim_views.d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
    }

    if (!s_prim_aoinds.empty())
    {
        gpuSafe(gpuMemcpy(slot.prim_views.d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(slot.prim_views.d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
        gpuSafe(gpuMemcpy(slot.prim_views.d_d_prim_aoinds, d_prim_aoinds.data(), d_prim_aoinds.size() * sizeof(uint32_t),
                          gpuMemcpyHostToDevice));
    }

    return slot.prim_views;
}

auto CGpuScreeningDeviceStorage::getBoysTables(const int64_t gpu_id) const -> const boysfunc::DeviceTables&
{
    return _slot(gpu_id).boys_tables;
}

auto CGpuScreeningDeviceStorage::getJPrep(const int64_t gpu_id) -> gpujprep::JPrepDeviceData&
{
    return _slot(gpu_id).j_prep;
}

auto CGpuScreeningDeviceStorage::getJPrep(const int64_t gpu_id) const -> const gpujprep::JPrepDeviceData&
{
    return _slot(gpu_id).j_prep;
}

auto CGpuScreeningDeviceStorage::getPrimViews(const int64_t gpu_id) const -> const PrimInfoDeviceViews&
{
    return _slot(gpu_id).prim_views;
}

auto CGpuScreeningDeviceStorage::releaseAll() -> void
{
    for (auto& slot : _slots)
    {
        slot.boys_buffer.release();
        slot.boys_tables = boysfunc::DeviceTables{};

        slot.j_prep_pools.mat_D_J.release();
        slot.j_prep_pools.mat_Q.release();
        slot.j_prep_pools.first_second_inds.release();
        slot.j_prep_pools.pair_data.release();
        slot.j_prep_pools.mat_Q_local.release();
        slot.j_prep_pools.first_second_inds_local.release();
        slot.j_prep_pools.pair_data_local.release();
        slot.j_prep = gpujprep::JPrepDeviceData{};

        slot.prim_info_buffer.release();
        slot.prim_aoinds_buffer.release();
        slot.prim_views = PrimInfoDeviceViews{};
    }

    _slots.clear();
    _configured = false;
}

CScreeningData::~CScreeningData()
{
    releaseDeviceStorage();
}

auto CScreeningData::getOrCreateDeviceStorage() const -> CGpuScreeningDeviceStorage&
{
    if (_device_storage == nullptr)
    {
        _device_storage = new CGpuScreeningDeviceStorage();
    }

    return *_device_storage;
}

auto CScreeningData::detachDeviceStorage() -> CGpuScreeningDeviceStorage*
{
    auto* storage   = _device_storage;
    _device_storage = nullptr;

    return storage;
}

auto CScreeningData::attachDeviceStorage(CGpuScreeningDeviceStorage* storage) -> void
{
    _device_storage = storage;
}

auto CScreeningData::releaseDeviceStorage() -> void
{
    if (_device_storage != nullptr)
    {
        _device_storage->releaseAll();
        delete _device_storage;
        _device_storage = nullptr;
    }
}

auto destroyGpuScreeningDeviceStorage(void* storage) -> void
{
    if (storage == nullptr)
    {
        return;
    }

    auto* device_storage = static_cast<CGpuScreeningDeviceStorage*>(storage);

    device_storage->releaseAll();
    delete device_storage;
}
