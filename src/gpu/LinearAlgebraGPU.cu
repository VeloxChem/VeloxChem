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
#include <hipsolver/hipsolver.h>

#endif

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "LinearAlgebraGPU.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "GpuDevices.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "StringFormat.hpp"

namespace gpu {  // gpu namespace

auto
getAvailableGpuMemory() -> double
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    size_t mem_free_bytes = 0;
    size_t mem_total_bytes = 0;

    gpuSafe(gpuMemGetInfo(&mem_free_bytes, &mem_total_bytes));

    return static_cast<double>(mem_free_bytes);
}

auto
computeDotProduct(const double* A, const double* B, const int64_t size_int64) -> double
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuSafe(preparePinnedMemcpyBuffer());  // per-thread pinned staging buffer for chunked copies

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto size = static_cast<size_t>(size_int64);

    double* d_A;
    double* d_B;

    gpuSafe(gpuMallocAsync(&d_A, size * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_B, size * sizeof(double), stream));

    gpuSafe(gpuMemcpyStaged(d_A, A, size * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyStaged(d_B, B, size * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpublasHandle_t handle;
    gpublasSafe(gpublasCreate(&handle));
    gpublasSafe(gpublasSetStream(handle, stream));

    auto n = static_cast<int32_t>(size);

    double dot_product;

#if defined(USE_CUDA)
    cublasSafe(cublasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));
#endif

    gpuSafe(gpuStreamSynchronize(stream));

    gpuSafe(gpuFreeAsync(d_A, stream));
    gpuSafe(gpuFreeAsync(d_B, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpublasSafe(gpublasDestroy(handle));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(releasePinnedMemcpyBuffer());  // session end: free per-thread pinned staging buffer
    gpuSafe(gpuDeviceSynchronize());

    return dot_product;
}

auto
computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& data_pointers, const int64_t size_int64) -> void
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuSafe(preparePinnedMemcpyBuffer());  // per-thread pinned staging buffer for chunked copies

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto size = static_cast<size_t>(size_int64);

    double* d_X;
    double* d_Y;

    gpuSafe(gpuMallocAsync(&d_X, size * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_Y, size * sizeof(double), stream));

    gpuSafe(gpuMemcpyStaged(d_Y, weighted_data, size * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpublasHandle_t handle;
    gpublasSafe(gpublasCreate(&handle));
    gpublasSafe(gpublasSetStream(handle, stream));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        gpuSafe(gpuMemcpyStaged(d_X, data_pointers[i], size * sizeof(double), gpuMemcpyHostToDevice, stream));

        auto n = static_cast<int32_t>(size);

#if defined(USE_CUDA)
        cublasSafe(cublasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
#elif defined(USE_HIP)
        hipblasSafe(hipblasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
#endif
    }

    gpuSafe(gpuMemcpyStaged(weighted_data, d_Y, size * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    gpuSafe(gpuFreeAsync(d_X, stream));
    gpuSafe(gpuFreeAsync(d_Y, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpublasSafe(gpublasDestroy(handle));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(releasePinnedMemcpyBuffer());  // session end: free per-thread pinned staging buffer
    gpuSafe(gpuDeviceSynchronize());
}

auto
computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                   const int64_t nmo_int64, const int64_t nao_int64, const std::string& trans_X) -> void
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuSafe(preparePinnedMemcpyBuffer());  // per-thread pinned staging buffer for chunked copies

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto nmo_size = static_cast<size_t>(nmo_int64);
    auto nao_size = static_cast<size_t>(nao_int64);

    double* d_A;
    double* d_B;
    double* d_C;

    gpuSafe(gpuMallocAsync(&d_A, nao_size * nao_size * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_B, nao_size * nao_size * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_C, nao_size * nao_size * sizeof(double), stream));

    gpuSafe(gpuMemcpyStaged(d_A, F, nao_size * nao_size * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyStaged(d_B, D, nao_size * nao_size * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpublasHandle_t handle;
    gpublasSafe(gpublasCreate(&handle));
    gpublasSafe(gpublasSetStream(handle, stream));

    double alpha = 1.0, beta = 0.0;

    auto nmo = static_cast<int32_t>(nmo_int64);
    auto nao = static_cast<int32_t>(nao_int64);

    // Note: we compute C^T = B^T * A^T since cublas/hipblas is column-major

    // D^T F^T (=> FD)
#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));
#endif

    // S^T(FD)^T (=> FDS)
    gpuSafe(gpuMemcpyStaged(d_A, S, nao_size * nao_size * sizeof(double), gpuMemcpyHostToDevice, stream));

#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));
#endif

    // FDS - (FDS)^T
    double geam_beta = -1.0;

#if defined(USE_CUDA)
    cublasSafe(cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_T, nao, nao, &alpha, d_B, nao, &geam_beta, d_B, nao, d_C, nao));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgeam(handle, HIPBLAS_OP_N, HIPBLAS_OP_T, nao, nao, &alpha, d_B, nao, &geam_beta, d_B, nao, d_C, nao));
#endif

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    gpuSafe(gpuMemcpyStaged(d_X, X, nmo_size * nao_size * sizeof(double), gpuMemcpyHostToDevice, stream));

#if defined(USE_CUDA)
    auto op_X  = (trans_X == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? CUBLAS_OP_T : CUBLAS_OP_N;
#elif defined(USE_HIP)
    auto op_X  = (trans_X == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? HIPBLAS_OP_T : HIPBLAS_OP_N;
#endif

    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    // TODO: double check basis set with linear dependency

    // let E == (FDS - SDF)
    // X^T E^T (=> EX)
#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, op_X, CUBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_C, nao, &beta, d_Y, nmo));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, op_X, HIPBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_C, nao, &beta, d_Y, nmo));
#endif

    // (EX)^T X (=> X^T(E)X)
#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_C, nmo));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_C, nmo));
#endif

    gpuSafe(gpuMemcpyStaged(errvec, d_C, nmo_size * nmo_size * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    gpuSafe(gpuFreeAsync(d_A, stream));
    gpuSafe(gpuFreeAsync(d_B, stream));
    gpuSafe(gpuFreeAsync(d_C, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpublasSafe(gpublasDestroy(handle));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(releasePinnedMemcpyBuffer());  // session end: free per-thread pinned staging buffer
    gpuSafe(gpuDeviceSynchronize());
}

auto
transformMatrix(double* transformed_F, const double* X, const double* F,
                const int64_t nmo_int64, const int64_t nao_int64, const std::string& trans_X) -> void
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuSafe(preparePinnedMemcpyBuffer());  // per-thread pinned staging buffer for chunked copies

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto nmo_size = static_cast<size_t>(nmo_int64);
    auto nao_size = static_cast<size_t>(nao_int64);

    double* d_F;
    double* d_X;
    double* d_Y;

    gpuSafe(gpuMallocAsync(&d_F, (nao_size * nao_size) * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_X, (nmo_size * nao_size) * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_Y, (nmo_size * nao_size) * sizeof(double), stream));

    gpuSafe(gpuMemcpyStaged(d_F, F, nao_size * nao_size * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyStaged(d_X, X, nmo_size * nao_size * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpublasHandle_t handle;
    gpublasSafe(gpublasCreate(&handle));
    gpublasSafe(gpublasSetStream(handle, stream));

    double alpha = 1.0, beta = 0.0;

    auto nmo = static_cast<int32_t>(nmo_int64);
    auto nao = static_cast<int32_t>(nao_int64);

#if defined(USE_CUDA)
    auto op_X = (trans_X == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? CUBLAS_OP_T : CUBLAS_OP_N;
#elif defined(USE_HIP)
    auto op_X = (trans_X == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? HIPBLAS_OP_T : HIPBLAS_OP_N;
#endif

    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;
    auto lda_XT = (trans_X == std::string("N")) ? nao : nmo;

    // Note: we compute C^T = B^T * A^T since cublas/hipblas is column-major

    // X^T F^T (=> FX)
#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, op_X, CUBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_F, nao, &beta, d_Y, nmo));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, op_X, HIPBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_F, nao, &beta, d_Y, nmo));
#endif

    // (FX)^T X (=> X^T(F)X)
#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_F, nmo));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_F, nmo));
#endif

    gpuSafe(gpuMemcpyStaged(transformed_F, d_F, nmo_size * nmo_size * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    gpuSafe(gpuFreeAsync(d_F, stream));
    gpuSafe(gpuFreeAsync(d_X, stream));
    gpuSafe(gpuFreeAsync(d_Y, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpublasSafe(gpublasDestroy(handle));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(releasePinnedMemcpyBuffer());  // session end: free per-thread pinned staging buffer
    gpuSafe(gpuDeviceSynchronize());
}

auto
computeMatrixMultiplication(double* C, const double* A, const double* B, const std::string& trans_A, const std::string& trans_B,
                            const int64_t m_int64, const int64_t k_int64, const int64_t n_int64) -> void
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuSafe(preparePinnedMemcpyBuffer());  // per-thread pinned staging buffer for chunked copies

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto m_size = static_cast<size_t>(m_int64);
    auto k_size = static_cast<size_t>(k_int64);
    auto n_size = static_cast<size_t>(n_int64);

    double* d_A;
    double* d_B;
    double* d_C;

    gpuSafe(gpuMallocAsync(&d_A, (m_size * k_size) * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_B, (k_size * n_size) * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_C, (m_size * n_size) * sizeof(double), stream));

    gpuSafe(gpuMemcpyStaged(d_A, A, m_size * k_size * sizeof(double), gpuMemcpyHostToDevice, stream));
    gpuSafe(gpuMemcpyStaged(d_B, B, k_size * n_size * sizeof(double), gpuMemcpyHostToDevice, stream));

    gpublasHandle_t handle;
    gpublasSafe(gpublasCreate(&handle));
    gpublasSafe(gpublasSetStream(handle, stream));

    double alpha = 1.0, beta = 0.0;

    auto m = static_cast<int32_t>(m_int64);
    auto k = static_cast<int32_t>(k_int64);
    auto n = static_cast<int32_t>(n_int64);

#if defined(USE_CUDA)
    auto op_A = (trans_A == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
#elif defined(USE_HIP)
    auto op_A = (trans_A == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
#endif

    // Note: we compute C^T = B^T * A^T since cublas/hipblas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

#if defined(USE_CUDA)
    cublasSafe(cublasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));
#elif defined(USE_HIP)
    hipblasSafe(hipblasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));
#endif

    gpuSafe(gpuMemcpyStaged(C, d_C, m_size * n_size * sizeof(double), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    gpuSafe(gpuFreeAsync(d_A, stream));
    gpuSafe(gpuFreeAsync(d_B, stream));
    gpuSafe(gpuFreeAsync(d_C, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpublasSafe(gpublasDestroy(handle));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(releasePinnedMemcpyBuffer());  // session end: free per-thread pinned staging buffer
    gpuSafe(gpuDeviceSynchronize());
}

auto
diagonalizeMatrix(double* A, double* D, const int64_t n_int64) -> void
{
    gpuSafe(gpuSetDevice(0));
    gpuSafe(gpuDeviceSynchronize());  // early context initialization after setdevice

    gpuSafe(preparePinnedMemcpyBuffer());  // per-thread pinned staging buffer for chunked copies

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto n_size = static_cast<size_t>(n_int64);

    gpuStream_t stream;
    gpuSafe(gpuStreamCreate(&stream));

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    gpuSafe(gpuMallocAsync(&d_A, n_size * n_size * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_D, n_size * sizeof(double), stream));
    gpuSafe(gpuMallocAsync(&d_info, sizeof(int32_t), stream));

    gpuSafe(gpuMemcpyStaged(d_A, A, n_size * n_size * sizeof(double), gpuMemcpyHostToDevice, stream));

    auto n = static_cast<int32_t>(n_int64);

#if defined(USE_CUDA)

    int32_t lwork, info;

    cusolverDnHandle_t handle;
    cusolverSafe(cusolverDnCreate(&handle));
    cusolverDnSetStream(handle, stream);

    cusolverSafe(cusolverDnDsyevd_bufferSize(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    gpuSafe(gpuMallocAsync(&d_work, static_cast<size_t>(lwork) * sizeof(double), stream));

    cusolverSafe(cusolverDnDsyevd(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    gpuSafe(gpuMemcpyStaged(A, d_A, n_size * n_size * sizeof(double), gpuMemcpyDeviceToHost, stream));
    gpuSafe(gpuMemcpyStaged(D, d_D, n_size * sizeof(double), gpuMemcpyDeviceToHost, stream));
    gpuSafe(gpuMemcpyStaged(&info, d_info, 1 * sizeof(int32_t), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    // TODO: check info

    // TODO: gpu wrapper for cusolver
    cusolverSafe(cusolverDnDestroy(handle));

#elif defined(USE_HIP)

    int32_t lwork, info;

    hipsolverHandle_t handle;
    hipsolverSafe(hipsolverCreate(&handle));
    hipsolverSetStream(handle, stream);

    hipsolverSafe(hipsolverDsyevd_bufferSize(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    //gpuSafe(gpuMallocAsync(&d_work, static_cast<size_t>(lwork), stream));
    gpuSafe(gpuMallocAsync(&d_work, static_cast<size_t>(lwork) * sizeof(double), stream));

    hipsolverSafe(hipsolverDsyevd(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    gpuSafe(gpuMemcpyStaged(A, d_A, n_size * n_size * sizeof(double), gpuMemcpyDeviceToHost, stream));
    gpuSafe(gpuMemcpyStaged(D, d_D, n_size * sizeof(double), gpuMemcpyDeviceToHost, stream));
    gpuSafe(gpuMemcpyStaged(&info, d_info, 1 * sizeof(int32_t), gpuMemcpyDeviceToHost, stream));

    gpuSafe(gpuStreamSynchronize(stream));

    // TODO check h_info

    // TODO: gpu wrapper for hipsolver
    hipsolverSafe(hipsolverDestroy(handle));

#endif

    gpuSafe(gpuFreeAsync(d_A, stream));
    gpuSafe(gpuFreeAsync(d_D, stream));
    gpuSafe(gpuFreeAsync(d_info, stream));
    gpuSafe(gpuFreeAsync(d_work, stream));

    gpuSafe(gpuStreamSynchronize(stream));
    gpuSafe(gpuStreamDestroy(stream));
    gpuSafe(releasePinnedMemcpyBuffer());  // session end: free per-thread pinned staging buffer
    gpuSafe(gpuDeviceSynchronize());
}

}  // namespace gpu
