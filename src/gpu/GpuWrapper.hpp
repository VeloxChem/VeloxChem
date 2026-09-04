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

#ifndef GpuWrapper_hpp
#define GpuWrapper_hpp

#if defined(USE_CUDA)

    #define gpuSafe(e)                          cudaSafe(e)
    #define gpublasSafe(e)                      cublasSafe(e)

    #define gpuDeviceProp                       cudaDeviceProp
    #define gpuGetDeviceCount(ptr)              cudaGetDeviceCount(ptr)
    #define gpuGetDeviceProperties(ptr, idx)    cudaGetDeviceProperties(ptr, idx)
    #define gpuSetDevice(idx)                   cudaSetDevice(idx)
    #define gpuMalloc(ptr, size)                cudaMalloc(ptr, size)
    #define gpuMallocAsync(ptr, size, s)        cudaMallocAsync(ptr, size, s)
    #define gpuFree(ptr)                        cudaFree(ptr)
    #define gpuFreeAsync(ptr, s)                cudaFreeAsync(ptr, s)
    #define gpuDeviceSynchronize()              cudaDeviceSynchronize()
    #define gpuMemGetInfo(p_free, p_total)      cudaMemGetInfo(p_free, p_total)
    #define gpuMemcpy(dst, src, size, kind)     cudaMemcpy(dst, src, size, kind)
    #define gpuMemcpyHostToDevice               cudaMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost               cudaMemcpyDeviceToHost
    #define gpuError_t                          cudaError_t
    #define gpuSuccess                          cudaSuccess
    #define gpuMemcpyKind                       cudaMemcpyKind
    // portable pinning: the per-thread staging buffer must stay page-locked for
    // whichever CUDA context issues the transfer, not only the allocating one
    #define gpuHostMalloc(ptr, size)            cudaHostAlloc(ptr, size, cudaHostAllocPortable)
    #define gpuHostFree(ptr)                    cudaFreeHost(ptr)
    #define gpuMemcpyAsync(dst, src, size, kind, s)        cudaMemcpyAsync(dst, src, size, kind, s)
    #define gpuMemcpyStaged(dst, src, size, kind, s)       ::gpu::stagedMemcpy(dst, src, size, kind, s)

    #define gpuStream_t                         cudaStream_t
    #define gpuStreamCreate(ptr)                cudaStreamCreate(ptr)
    #define gpuStreamSynchronize(s)             cudaStreamSynchronize(s)
    #define gpuStreamDestroy(s)                 cudaStreamDestroy(s)

    #define gpublasHandle_t                     cublasHandle_t
    #define gpublasCreate(ptr)                  cublasCreate(ptr)
    #define gpublasSetStream(h, s)              cublasSetStream(h, s)
    #define gpublasDestroy(h)                   cublasDestroy(h)

#elif defined(USE_HIP)

    #define gpuSafe(e)                          hipSafe(e)
    #define gpublasSafe(e)                      hipblasSafe(e)

    #define gpuDeviceProp                       hipDeviceProp_t
    #define gpuGetDeviceCount(ptr)              hipGetDeviceCount(ptr)
    #define gpuGetDeviceProperties(ptr, idx)    hipGetDeviceProperties(ptr, idx)
    #define gpuSetDevice(idx)                   hipSetDevice(idx)
    #define gpuMalloc(ptr, size)                hipMalloc(ptr, size)
    #define gpuMallocAsync(ptr, size, s)        hipMallocAsync(ptr, size, s)
    #define gpuFree(ptr)                        hipFree(ptr)
    #define gpuFreeAsync(ptr, s)                hipFreeAsync(ptr, s)
    #define gpuDeviceSynchronize()              hipDeviceSynchronize()
    #define gpuMemGetInfo(p_free, p_total)      hipMemGetInfo(p_free, p_total)
    #define gpuMemcpy(dst, src, size, kind)     hipMemcpy(dst, src, size, kind)
    #define gpuMemcpyHostToDevice               hipMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost               hipMemcpyDeviceToHost
    #define gpuError_t                          hipError_t
    #define gpuSuccess                          hipSuccess
    #define gpuMemcpyKind                       hipMemcpyKind
    // portable pinning: the per-thread staging buffer must stay page-locked for
    // whichever HIP device issues the transfer, not only the allocating one
    #define gpuHostMalloc(ptr, size)            hipHostMalloc(ptr, size, hipHostMallocPortable)
    #define gpuHostFree(ptr)                    hipHostFree(ptr)
    #define gpuMemcpyAsync(dst, src, size, kind, s)        hipMemcpyAsync(dst, src, size, kind, s)
    #define gpuMemcpyStaged(dst, src, size, kind, s)       ::gpu::stagedMemcpy(dst, src, size, kind, s)

    #define gpuStream_t                         hipStream_t
    #define gpuStreamCreate(ptr)                hipStreamCreate(ptr)
    #define gpuStreamSynchronize(s)             hipStreamSynchronize(s)
    #define gpuStreamDestroy(s)                 hipStreamDestroy(s)

    #define gpublasHandle_t                     hipblasHandle_t
    #define gpublasCreate(ptr)                  hipblasCreate(ptr)
    #define gpublasSetStream(h, s)              hipblasSetStream(h, s)
    #define gpublasDestroy(h)                   hipblasDestroy(h)

#else

  #error "Please define either USE_CUDA or USE_HIP"

#endif

// Pinned-staging implementation of gpuMemcpyStaged.  It is kept in this header
// next to the gpu* backend macros it is built on.

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <string>

#include "ErrorHandler.hpp"
#include "GpuRuntime.hpp"
#include "GpuSafeChecks.hpp"

namespace gpu {  // gpu namespace

// gpuMemcpyStaged host<->device copies are staged in chunks through a pinned
// host buffer: pageable host memory is unreliable for asynchronous device
// transfers on multi-NUMA systems.  Each chunk is completed before the buffer
// is reused, and the whole copy is complete when the call returns.

// maximum number of bytes staged through pinned host memory per chunk
constexpr size_t kStagedMemcpyChunkBytes = 8 * 1024 * 1024;  // 8 MB

// Copy sessions: every driver session (an omp parallel GPU region, or a serial
// GPU driver function) opens with preparePinnedMemcpyBuffer after gpuSetDevice
// and closes with releasePinnedMemcpyBuffer after its last copy; all chunked
// copies of the session share the calling thread's pinned buffer.  A copy
// outside a session is a programming error (abort).
inline char*&
stagedMemcpyPinnedBuffer()
{
    static thread_local char* h_pinned = nullptr;

    return h_pinned;
}

// Open the copy session of the calling thread, allocating its pinned staging
// buffer on first use.  Call after gpuSetDevice at the start of the driver
// session; allocating on the thread that drives the device gives best-effort
// NUMA locality.  Idempotent; every session must end with
// releasePinnedMemcpyBuffer.
inline gpuError_t
preparePinnedMemcpyBuffer()
{
    if (stagedMemcpyPinnedBuffer() == nullptr)
    {
        gpuSafe(gpuHostMalloc(reinterpret_cast<void**>(&stagedMemcpyPinnedBuffer()), kStagedMemcpyChunkBytes));
    }

    return gpuSuccess;
}

// Close the copy session of the calling thread: free its pinned staging
// buffer after the session's last copy (null-safe and idempotent).
inline gpuError_t
releasePinnedMemcpyBuffer()
{
    if (stagedMemcpyPinnedBuffer() != nullptr)
    {
        gpuSafe(gpuHostFree(stagedMemcpyPinnedBuffer()));

        stagedMemcpyPinnedBuffer() = nullptr;
    }

    return gpuSuccess;
}

inline gpuError_t
stagedMemcpyHostToDevice(void* d_ptr, const void* h_ptr, const size_t byte_count, gpuStream_t stream)
{
    if (byte_count == 0) return gpuSuccess;

    char* h_pinned = stagedMemcpyPinnedBuffer();

    if (h_pinned == nullptr)
    {
        // copies require an open session, see preparePinnedMemcpyBuffer
        errors::assertMsgCritical(false,
            std::string("gpu::") + std::string(__func__) +
                ": no pinned host staging buffer on this thread; call gpu::preparePinnedMemcpyBuffer() "
                "before the first gpuMemcpyStaged of the driver session");
    }

    auto*       d_bytes = static_cast<char*>(d_ptr);
    const auto* h_bytes = static_cast<const char*>(h_ptr);

    size_t offset = 0;

    while (offset < byte_count)
    {
        const size_t copy_bytes = std::min(kStagedMemcpyChunkBytes, byte_count - offset);

        if (offset > 0)
        {
            // wait for the previous staged transfer before refilling the buffer
            gpuSafe(gpuStreamSynchronize(stream));
        }

        std::memcpy(h_pinned, h_bytes + offset, copy_bytes);

        gpuSafe(gpuMemcpyAsync(d_bytes + offset, h_pinned, copy_bytes, gpuMemcpyHostToDevice, stream));

        offset += copy_bytes;
    }

    // the last staged transfer must be complete before the buffer is reused
    gpuSafe(gpuStreamSynchronize(stream));

    return gpuSuccess;
}

inline gpuError_t
stagedMemcpyDeviceToHost(void* h_ptr, const void* d_ptr, const size_t byte_count, gpuStream_t stream)
{
    if (byte_count == 0) return gpuSuccess;

    char* h_pinned = stagedMemcpyPinnedBuffer();

    if (h_pinned == nullptr)
    {
        // copies require an open session, see preparePinnedMemcpyBuffer
        errors::assertMsgCritical(false,
            std::string("gpu::") + std::string(__func__) +
                ": no pinned host staging buffer on this thread; call gpu::preparePinnedMemcpyBuffer() "
                "before the first gpuMemcpyStaged of the driver session");
    }

    auto*       h_bytes = static_cast<char*>(h_ptr);
    const auto* d_bytes = static_cast<const char*>(d_ptr);

    size_t offset = 0;

    while (offset < byte_count)
    {
        const size_t copy_bytes = std::min(kStagedMemcpyChunkBytes, byte_count - offset);

        gpuSafe(gpuMemcpyAsync(h_pinned, d_bytes + offset, copy_bytes, gpuMemcpyDeviceToHost, stream));

        // the staged chunk must be fully arrived before it is copied out or refilled
        gpuSafe(gpuStreamSynchronize(stream));

        std::memcpy(h_bytes + offset, h_pinned, copy_bytes);

        offset += copy_bytes;
    }

    return gpuSuccess;
}

// Memcpy dispatcher: gpuMemcpyStaged handles host<->device transfers, staged in
// chunks through the pinned-staging copies above.
inline gpuError_t
stagedMemcpy(void* dst, const void* src, const size_t byte_count, const gpuMemcpyKind kind, gpuStream_t stream)
{
    if (kind == gpuMemcpyHostToDevice)
    {
        return stagedMemcpyHostToDevice(dst, src, byte_count, stream);
    }
    else if (kind == gpuMemcpyDeviceToHost)
    {
        return stagedMemcpyDeviceToHost(dst, src, byte_count, stream);
    }

    errors::assertMsgCritical(false,
        std::string("gpu::") + std::string(__func__) +
            ": unsupported transfer kind; gpuMemcpyStaged handles host<->device copies only");

    return gpuSuccess;  // never reached: assertMsgCritical aborts
}

}  // namespace gpu

#endif /* GpuWrapper_hpp */
