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

    #define gpuDeviceProp                       cudaDeviceProp
    #define gpuGetDeviceCount(ptr)              cudaGetDeviceCount(ptr)
    #define gpuGetDeviceProperties(ptr, idx)    cudaGetDeviceProperties(ptr, idx)
    #define gpuSetDevice(idx)                   cudaSetDevice(idx)
    #define gpuMalloc(ptr, size)                cudaMalloc(ptr, size)
    #define gpuFree(ptr)                        cudaFree(ptr)
    #define gpuDeviceSynchronize()              cudaDeviceSynchronize()
    #define gpuMemcpy(dst, src, size, kind)     cudaMemcpy(dst, src, size, kind)
    #define gpuMemcpyHostToDevice               cudaMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost               cudaMemcpyDeviceToHost
    #define gpuHostMallocDefault                cudaHostMallocDefault
    #define gpuHostMalloc(ptr, size, f)         cudaHostMalloc(ptr, size, f)
    #define gpuHostFree(ptr)                    cudaHostFree(ptr)
    #define gpuMemcpyAsync(dst, src, size, kind, stream) \
            cudaMemcpyAsync(dst, src, size, kind, stream)
    #define gpuEvent_t                          cudaEvent_t
    #define gpuEventDestroy(e)                  cudaEventDestroy(e)
    #define gpuEventCreateWithFlags             cudaEventCreateWithFlags
    #define gpuEventDisableTiming               cudaEventDisableTiming
    #define gpuEventSynchronize(e)              cudaEventSynchronize(e)
    #define gpuEventRecord(e)                   cudaEventRecord(e)
    #define gpuEventRecordStream(e, s)          cudaEventRecord(e, s)
    #define gpuStream_t                         cudaStream_t
    #define gpuStreamDestroy(s)                 cudaStreamDestroy(s)
    #define gpuStreamCreate(s)                  cudaStreamCreate(s)
    #define gpuDeviceGetStreamPriorityRange(ptr, p) \
            cudaDeviceGetStreamPriorityRange(ptr, p)
    #define gpuStreamDefault                    cudaStreamDefault
    #define gpuStreamCreateWithPriority(s, f, p) \
            cudaStreamCreateWithPriority(s, f, p)

#elif defined(USE_HIP)

    #define gpuSafe(e)                          hipSafe(e)

    #define gpuDeviceProp                       hipDeviceProp_t
    #define gpuGetDeviceCount(ptr)              hipGetDeviceCount(ptr)
    #define gpuGetDeviceProperties(ptr, idx)    hipGetDeviceProperties(ptr, idx)
    #define gpuSetDevice(idx)                   hipSetDevice(idx)
    #define gpuMalloc(ptr, size)                hipMalloc(ptr, size)
    #define gpuFree(ptr)                        hipFree(ptr)
    #define gpuDeviceSynchronize()              hipDeviceSynchronize()
    #define gpuMemcpy(dst, src, size, kind)     hipMemcpy(dst, src, size, kind)
    #define gpuMemcpyHostToDevice               hipMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost               hipMemcpyDeviceToHost
    #define gpuHostMallocDefault                hipHostMallocDefault
    #define gpuHostMalloc(ptr, size, f)         hipHostMalloc(ptr, size, f)
    #define gpuHostFree(ptr)                    hipHostFree(ptr)
    #define gpuMemcpyAsync(dst, src, size, kind, stream) \
            hipMemcpyAsync(dst, src, size, kind, stream)
    #define gpuEvent_t                          hipEvent_t
    #define gpuEventDestroy(e)                  hipEventDestroy(e)
    #define gpuEventCreateWithFlags             hipEventCreateWithFlags
    #define gpuEventDisableTiming               hipEventDisableTiming
    #define gpuEventSynchronize(e)              hipEventSynchronize(e)
    #define gpuEventRecord(e)                   hipEventRecord(e)
    #define gpuEventRecordStream(e, s)          hipEventRecord(e, s)
    #define gpuStream_t                         hipStream_t
    #define gpuStreamDestroy(s)                 hipStreamDestroy(s)
    #define gpuStreamCreate(s)                  hipStreamCreate(s)
    #define gpuDeviceGetStreamPriorityRange(ptr, p) \
            hipDeviceGetStreamPriorityRange(ptr, p)
    #define gpuStreamDefault                    hipStreamDefault
    #define gpuStreamCreateWithPriority(s, f, p) \
            hipStreamCreateWithPriority(s, f, p)
#else

  #error "Please define either USE_CUDA or USE_HIP"

#endif

#endif /* GpuWrapper_hpp */
