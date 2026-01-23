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

#include "ChunkedMemcpyGPU.hpp"
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

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto size = static_cast<size_t>(size_int64);

    double* d_data;

    gpuSafe(gpuMalloc(&d_data, size * 2 * sizeof(double)));

    double* d_A = d_data;
    double* d_B = d_A + size;

    gpu::chunkedMemcpyHostToDevice<double>(d_A, A, size);
    gpu::chunkedMemcpyHostToDevice<double>(d_B, B, size);

#if defined(USE_CUDA)

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    auto n = static_cast<int32_t>(size);

    double dot_product;
    cublasSafe(cublasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));

    cublasSafe(cublasDestroy(handle));

#elif defined(USE_HIP)

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    auto n = static_cast<int32_t>(size);

    double dot_product;
    hipblasSafe(hipblasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));

    hipblasSafe(hipblasDestroy(handle));

#endif

    gpuSafe(gpuFree(d_data));

    return dot_product;
}

auto
computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& data_pointers, const int64_t size_int64) -> void
{
    gpuSafe(gpuSetDevice(0));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto size = static_cast<size_t>(size_int64);

    double* d_data;

    gpuSafe(gpuMalloc(&d_data, size * 2 * sizeof(double)));

    double* d_X = d_data;
    double* d_Y = d_X + size;

    gpu::chunkedMemcpyHostToDevice<double>(d_Y, weighted_data, size);

#if defined(USE_CUDA)

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        gpu::chunkedMemcpyHostToDevice<double>(d_X, data_pointers[i], size);

        auto n = static_cast<int32_t>(size);

        cublasSafe(cublasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
    }

    cublasSafe(cublasDestroy(handle));

#elif defined(USE_HIP)

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        gpu::chunkedMemcpyHostToDevice<double>(d_X, data_pointers[i], size);

        auto n = static_cast<int32_t>(size);

        hipblasSafe(hipblasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
    }

    hipblasSafe(hipblasDestroy(handle));

#endif

    gpu::chunkedMemcpyDeviceToHost<double>(weighted_data, d_Y, size);

    gpuSafe(gpuFree(d_data));
}

auto
computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                   const int64_t nmo_int64, const int64_t nao_int64, const std::string& trans_X) -> void
{
    gpuSafe(gpuSetDevice(0));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto nmo_size = static_cast<size_t>(nmo_int64);
    auto nao_size = static_cast<size_t>(nao_int64);

    double* d_data;

    gpuSafe(gpuMalloc(&d_data, nao_size * nao_size * 3 * sizeof(double)));

    double* d_A = d_data;
    double* d_B = d_A + nao_size * nao_size;
    double* d_C = d_B + nao_size * nao_size;

    gpu::chunkedMemcpyHostToDevice<double>(d_A, F, nao_size * nao_size);
    gpu::chunkedMemcpyHostToDevice<double>(d_B, D, nao_size * nao_size);

#if defined(USE_CUDA)

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto nmo = static_cast<int32_t>(nmo_int64);
    auto nao = static_cast<int32_t>(nao_int64);

    // Note: we compute C^T = B^T * A^T since cublas is column-major

    // D^T F^T (=> FD)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));

    // S^T(FD)^T (=> FDS)
    gpu::chunkedMemcpyHostToDevice<double>(d_A, S, nao_size * nao_size);

    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));

    cudaSafe(cudaDeviceSynchronize());

    // FDS - (FDS)^T
    beta = -1.0;

    cublasSafe(cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_T, nao, nao, &alpha, d_B, nao, &beta, d_B, nao, d_C, nao));

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    gpu::chunkedMemcpyHostToDevice<double>(d_X, X, nmo_size * nao_size);

    auto op_X  = (trans_X == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? CUBLAS_OP_T : CUBLAS_OP_N;

    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    beta = 0.0;

    // TODO: double check basis set with linear dependency

    // let E == (FDS - SDF)
    // X^T E^T (=> EX)
    cublasSafe(cublasDgemm(handle, op_X, CUBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_C, nao, &beta, d_Y, nmo));

    cudaSafe(cudaDeviceSynchronize());

    // (EX)^T X (=> X^T(E)X)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_C, nmo));

    gpu::chunkedMemcpyDeviceToHost<double>(errvec, d_C, nmo_size * nmo_size);

    cublasSafe(cublasDestroy(handle));

#elif defined(USE_HIP)

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto nmo = static_cast<int32_t>(nmo_int64);
    auto nao = static_cast<int32_t>(nao_int64);

    // Note: we compute C^T = B^T * A^T since hipblas is column-major

    // D^T F^T (=> FD)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));

    // S^T(FD)^T (=> FDS)
    gpu::chunkedMemcpyHostToDevice<double>(d_A, S, nao_size * nao_size);

    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));

    hipSafe(hipDeviceSynchronize());

    // FDS - (FDS)^T
    beta = -1.0;

    hipblasSafe(hipblasDgeam(handle, HIPBLAS_OP_N, HIPBLAS_OP_T, nao, nao, &alpha, d_B, nao, &beta, d_B, nao, d_C, nao));

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    gpu::chunkedMemcpyHostToDevice<double>(d_X, X, nmo_size * nao_size);

    auto op_X  = (trans_X == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_XT = (trans_X == std::string("N")) ? HIPBLAS_OP_T : HIPBLAS_OP_N;

    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    beta = 0.0;

    // TODO: double check basis set with linear dependency

    // let E == (FDS - SDF)
    // X^T E^T (=> EX)
    hipblasSafe(hipblasDgemm(handle, op_X, HIPBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_C, nao, &beta, d_Y, nmo));

    hipSafe(hipDeviceSynchronize());

    // (EX)^T X (=> X^T(E)X)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_C, nmo));

    gpu::chunkedMemcpyDeviceToHost<double>(errvec, d_C, nmo_size * nmo_size);

    hipblasSafe(hipblasDestroy(handle));

#endif

    gpuSafe(gpuFree(d_data));
}

auto
transformMatrix(double* transformed_F, const double* X, const double* F,
                const int64_t nmo_int64, const int64_t nao_int64, const std::string& trans_X) -> void
{
    gpuSafe(gpuSetDevice(0));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto nmo_size = static_cast<size_t>(nmo_int64);
    auto nao_size = static_cast<size_t>(nao_int64);

    double* d_data;

    gpuSafe(gpuMalloc(&d_data, (nao_size * nao_size + nmo_size * nao_size + nmo_size * nao_size) * sizeof(double)));

    double* d_F = d_data;
    double* d_X = d_F + nao_size * nao_size;
    double* d_Y = d_X + nmo_size * nao_size;

    gpu::chunkedMemcpyHostToDevice<double>(d_F, F, nao_size * nao_size);
    gpu::chunkedMemcpyHostToDevice<double>(d_X, X, nmo_size * nao_size);

#if defined(USE_CUDA)

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto nmo = static_cast<int32_t>(nmo_int64);
    auto nao = static_cast<int32_t>(nao_int64);

    auto op_X = (trans_X == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    auto op_XT = (trans_X == std::string("N")) ? CUBLAS_OP_T : CUBLAS_OP_N;
    auto lda_XT = (trans_X == std::string("N")) ? nao : nmo;

    // Note: we compute C^T = B^T * A^T since cublas is column-major

    // X^T F^T (=> FX)
    cublasSafe(cublasDgemm(handle, op_X, CUBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_F, nao, &beta, d_Y, nmo));

    cudaSafe(cudaDeviceSynchronize());

    // (FX)^T X (=> X^T(F)X)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_F, nmo));

    gpu::chunkedMemcpyDeviceToHost<double>(transformed_F, d_F, nmo_size * nmo_size);

    cublasSafe(cublasDestroy(handle));

#elif defined(USE_HIP)

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto nmo = static_cast<int32_t>(nmo_int64);
    auto nao = static_cast<int32_t>(nao_int64);

    auto op_X = (trans_X == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto lda_X = (trans_X == std::string("N")) ? nmo : nao;

    auto op_XT = (trans_X == std::string("N")) ? HIPBLAS_OP_T : HIPBLAS_OP_N;
    auto lda_XT = (trans_X == std::string("N")) ? nao : nmo;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major

    // X^T F^T (=> FX)
    hipblasSafe(hipblasDgemm(handle, op_X, HIPBLAS_OP_N, nmo, nao, nao, &alpha, d_X, lda_X, d_F, nao, &beta, d_Y, nmo));

    hipSafe(hipDeviceSynchronize());

    // (FX)^T X (=> X^T(F)X)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, op_XT, nmo, nmo, nao, &alpha, d_Y, nmo, d_X, lda_X, &beta, d_F, nmo));

    gpu::chunkedMemcpyDeviceToHost<double>(transformed_F, d_F, nmo_size * nmo_size);

    hipblasSafe(hipblasDestroy(handle));

#endif

    gpuSafe(gpuFree(d_data));
}

auto
computeMatrixMultiplication(double* C, const double* A, const double* B, const std::string& trans_A, const std::string& trans_B,
                            const int64_t m_int64, const int64_t k_int64, const int64_t n_int64) -> void
{
    gpuSafe(gpuSetDevice(0));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto m_size = static_cast<size_t>(m_int64);
    auto k_size = static_cast<size_t>(k_int64);
    auto n_size = static_cast<size_t>(n_int64);

    double* d_data;

    gpuSafe(gpuMalloc(&d_data, (m_size * k_size + k_size * n_size + m_size * n_size) * sizeof(double)));

    double* d_A = d_data;
    double* d_B = d_A + m_size * k_size;
    double* d_C = d_B + k_size * n_size;

    gpu::chunkedMemcpyHostToDevice<double>(d_A, A, m_size * k_size);
    gpu::chunkedMemcpyHostToDevice<double>(d_B, B, k_size * n_size);

#if defined(USE_CUDA)

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto m = static_cast<int32_t>(m_int64);
    auto k = static_cast<int32_t>(k_int64);
    auto n = static_cast<int32_t>(n_int64);

    auto op_A = (trans_A == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;

    // Note: we compute C^T = B^T * A^T since cublas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

    cublasSafe(cublasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));

    gpu::chunkedMemcpyDeviceToHost<double>(C, d_C, m_size * n_size);

    cublasSafe(cublasDestroy(handle));

#elif defined(USE_HIP)

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto m = static_cast<int32_t>(m_int64);
    auto k = static_cast<int32_t>(k_int64);
    auto n = static_cast<int32_t>(n_int64);

    auto op_A = (trans_A == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

    hipblasSafe(hipblasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));

    gpu::chunkedMemcpyDeviceToHost<double>(C, d_C, m_size * n_size);

    hipblasSafe(hipblasDestroy(handle));

#endif

    gpuSafe(gpuFree(d_data));
}

auto
diagonalizeMatrix(double* A, double* D, const int64_t n_int64) -> void
{
    gpuSafe(gpuSetDevice(0));

    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto n_size = static_cast<size_t>(n_int64);

#if defined(USE_CUDA)

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    cudaSafe(cudaMalloc(&d_A, n_size * n_size * sizeof(double)));
    cudaSafe(cudaMalloc(&d_D, n_size * sizeof(double)));
    cudaSafe(cudaMalloc(&d_info, sizeof(int32_t)));

    gpu::chunkedMemcpyHostToDevice<double>(d_A, A, n_size * n_size);

    auto n = static_cast<int32_t>(n_int64);
    int32_t lwork, info;

    cusolverDnHandle_t handle;
    cusolverSafe(cusolverDnCreate(&handle));

    cusolverSafe(cusolverDnDsyevd_bufferSize(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    cudaSafe(cudaMalloc(&d_work, static_cast<size_t>(lwork) * sizeof(double)));

    cusolverSafe(cusolverDnDsyevd(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    cusolverSafe(cusolverDnDestroy(handle));

    gpu::chunkedMemcpyDeviceToHost<double>(A, d_A, n_size * n_size);
    gpu::chunkedMemcpyDeviceToHost<double>(D, d_D, n_size);
    gpu::chunkedMemcpyDeviceToHost<int32_t>(&info, d_info, 1);

    // TODO: check info

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_D));
    cudaSafe(cudaFree(d_info));
    cudaSafe(cudaFree(d_work));

#elif defined(USE_HIP)

    magmaSafe(magma_init());

    magma_setdevice(0);

    magma_int_t n = static_cast<magma_int_t>(n_int64);
    //magma_int_t ldda = magma_roundup(n, 32);

    double *d_A;
    //hipSafe(hipMalloc(&d_A, n * ldda * sizeof(double)));
    //hipblasSafe(hipblasSetMatrix(n, n, sizeof(double), A, n, d_A, ldda));
    hipSafe(hipMalloc(&d_A, n_size * n_size * sizeof(double)));
    gpu::chunkedMemcpyHostToDevice<double>(d_A, A, n_size * n_size);

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n_size * n_size * sizeof(double)));
    hipSafe(hipHostMalloc(&work, static_cast<size_t>(lwork) * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, static_cast<size_t>(liwork) * sizeof(magma_int_t)));

    magma_int_t info;
    //magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, ldda, D, wA, n, work, lwork, iwork, liwork, &info);
    magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, n, D, wA, n, work, lwork, iwork, liwork, &info);

    if (info != 0)
    {
        std::stringstream ss;
        ss << "gpu::diagonalizeMatrix: (magma error) " << magma_strerror(info);
        errors::assertMsgCritical(false, ss.str());
    }

    //hipblasSafe(hipblasGetMatrix(n, n, sizeof(double), d_A, ldda, A, n));
    gpu::chunkedMemcpyDeviceToHost<double>(A, d_A, n_size * n_size);

    hipSafe(hipFree(d_A));

    hipSafe(hipHostFree(wA));
    hipSafe(hipHostFree(work));
    hipSafe(hipHostFree(iwork));

    magmaSafe(magma_finalize());

#endif
}

#if defined(USE_HIP)
auto
diagonalizeMatrixMultiGPU(double* A, double* D, const int64_t n_int64, const int64_t num_gpus_per_node) -> void
{
    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string(__func__) + std::string(": should not be called in omp parallel region"));

    auto n_size = static_cast<size_t>(n_int64);

    magmaSafe(magma_init());

    magma_int_t n = static_cast<magma_int_t>(n_int64);

    magma_int_t ngpu = static_cast<magma_int_t>(num_gpus_per_node);

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n_size * n_size * sizeof(double)));
    hipSafe(hipHostMalloc(&work, static_cast<size_t>(lwork) * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, static_cast<size_t>(liwork) * sizeof(magma_int_t)));

    magma_int_t info;
    magma_dsyevd_m(ngpu, MagmaVec, MagmaUpper, n, A, n, D, work, lwork, iwork, liwork, &info);

    if (info != 0)
    {
        std::stringstream ss;
        ss << "gpu::diagonalizeMatrixMultiGPU: (magma error) " << magma_strerror(info);
        errors::assertMsgCritical(false, ss.str());
    }

    hipSafe(hipHostFree(wA));
    hipSafe(hipHostFree(work));
    hipSafe(hipHostFree(iwork));

    magmaSafe(magma_finalize());
}
#endif

}  // namespace gpu
