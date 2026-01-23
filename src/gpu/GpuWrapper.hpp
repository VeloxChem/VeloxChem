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

// TODO: Use hipFreeAsync for gpuFreeAsync with newer rocm
//       At the moment we use hipFree due to rocm 6.3

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
    #define gpuMemcpyAsync(dst, src, size, kind, s)     cudaMemcpyAsync(dst, src, size, kind, s)
    #define gpuMemcpyHostToDevice               cudaMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost               cudaMemcpyDeviceToHost

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
    #define gpuFreeAsync(ptr, s)                hipFree(ptr)
    #define gpuDeviceSynchronize()              hipDeviceSynchronize()
    #define gpuMemGetInfo(p_free, p_total)      hipMemGetInfo(p_free, p_total)
    #define gpuMemcpy(dst, src, size, kind)     hipMemcpy(dst, src, size, kind)
    #define gpuMemcpyAsync(dst, src, size, kind, s)     hipMemcpyAsync(dst, src, size, kind, s)
    #define gpuMemcpyHostToDevice               hipMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost               hipMemcpyDeviceToHost

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

#endif /* GpuWrapper_hpp */
