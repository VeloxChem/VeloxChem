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

#ifndef GpuDeviceBufferStorage_hpp
#define GpuDeviceBufferStorage_hpp

#include <cstddef>
#include <cstdint>

#include "GpuRuntime.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"

/**
 Grow-only device buffer. Reuses the same allocation when the requested size fits
 in the current capacity. Reallocates only when a larger size is required or when
 the target device changes.
 */
class CGpuDeviceBuffer
{
    void*   _ptr{nullptr};
    size_t  _capacity_bytes{0};
    int     _device_id{-1};

   public:
    CGpuDeviceBuffer() = default;

    CGpuDeviceBuffer(const CGpuDeviceBuffer&)            = delete;
    CGpuDeviceBuffer& operator=(const CGpuDeviceBuffer&) = delete;

    CGpuDeviceBuffer(CGpuDeviceBuffer&& other) noexcept
        : _ptr(other._ptr),
          _capacity_bytes(other._capacity_bytes),
          _device_id(other._device_id)
    {
        other._ptr             = nullptr;
        other._capacity_bytes  = 0;
        other._device_id       = -1;
    }

    auto operator=(CGpuDeviceBuffer&& other) noexcept -> CGpuDeviceBuffer&
    {
        if (this != &other)
        {
            release();

            _ptr            = other._ptr;
            _capacity_bytes = other._capacity_bytes;
            _device_id      = other._device_id;

            other._ptr             = nullptr;
            other._capacity_bytes  = 0;
            other._device_id       = -1;
        }

        return *this;
    }

    ~CGpuDeviceBuffer() { release(); }

    auto ptr() const -> void* { return _ptr; }

    auto capacity() const -> size_t { return _capacity_bytes; }

    auto deviceId() const -> int { return _device_id; }

    auto empty() const -> bool { return _ptr == nullptr; }

    auto ensure(const int device_id, const size_t nbytes) -> void*
    {
        if (nbytes == 0)
        {
            return _ptr;
        }

        if (_ptr != nullptr && device_id == _device_id && nbytes <= _capacity_bytes)
        {
            return _ptr;
        }

        release();

        gpuSafe(gpuSetDevice(device_id));
        gpuSafe(gpuMalloc(&_ptr, nbytes));

        _capacity_bytes = nbytes;
        _device_id      = device_id;

        return _ptr;
    }

    auto upload(const int device_id, const void* host, const size_t nbytes) -> void*
    {
        auto* device_ptr = ensure(device_id, nbytes);

        if (host != nullptr && nbytes > 0)
        {
            gpuSafe(gpuSetDevice(device_id));
            gpuSafe(gpuMemcpy(device_ptr, host, nbytes, gpuMemcpyHostToDevice));
        }

        return device_ptr;
    }

    auto release() -> void
    {
        if (_ptr == nullptr)
        {
            return;
        }

        if (_device_id >= 0)
        {
            gpuSafe(gpuSetDevice(_device_id));
        }

        gpuSafe(gpuFree(_ptr));

        _ptr            = nullptr;
        _capacity_bytes = 0;
        _device_id      = -1;
    }

    auto releaseOwnership() -> void*
    {
        void* ptr = _ptr;

        _ptr            = nullptr;
        _capacity_bytes = 0;
        _device_id      = -1;

        return ptr;
    }
};

#endif /* GpuDeviceBufferStorage_hpp */
