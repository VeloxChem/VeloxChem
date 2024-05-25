//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include <cublas_v2.h>
#include <cusolverDn.h>

#include "ErrorHandler.hpp"
#include "GpuSafeChecks.hpp"
#include "LinearAlgebraGPU.hpp"
#include "MpiFunc.hpp"

namespace gpu {  // gpu namespace

auto
computeDotProduct(const double* A, const double* B, const int64_t size) -> double
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeDotProduct: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_A, *d_B;

    cudaSafe(cudaMalloc(&d_A, n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_B, n * sizeof(double)));

    cudaSafe(cudaMemcpy(d_A, A, n * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_B, B, n * sizeof(double), cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double dot_product;
    cublasSafe(cublasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_B));

    return dot_product;
}

auto
computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& data_pointers, const int64_t size) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeWeightedSum: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_X, *d_Y;

    cudaSafe(cudaMalloc(&d_X, n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Y, n * sizeof(double)));

    cudaSafe(cudaMemcpy(d_Y, weighted_data, n * sizeof(double), cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        cudaSafe(cudaMemcpy(d_X, data_pointers[i], n * sizeof(double), cudaMemcpyHostToDevice));

        cublasSafe(cublasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
    }

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaMemcpy(weighted_data, d_Y, n * sizeof(double), cudaMemcpyDeviceToHost));

    cudaSafe(cudaFree(d_X));
    cudaSafe(cudaFree(d_Y));
}

auto
computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                   const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeErrorVector: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_A, *d_B, *d_C;

    cudaSafe(cudaMalloc(&d_A, nao * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_B, nao * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_C, nao * nao * sizeof(double)));

    cudaSafe(cudaMemcpy(d_A, F, nao * nao * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_B, D, nao * nao * sizeof(double), cudaMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

    // Note: we compute C^T = B^T * A^T since cublas is column-major

    // D^T F^T (=> FD)
    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));

    // S^T(FD)^T (=> FDS)
    cudaSafe(cudaMemcpy(d_A, S, nao * nao * sizeof(double), cudaMemcpyHostToDevice));

    cublasSafe(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));

    cudaSafe(cudaDeviceSynchronize());

    // FDS - (FDS)^T
    beta = -1.0;

    cublasSafe(cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_T, nao, nao, &alpha, d_B, nao, &beta, d_B, nao, d_C, nao));

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    cudaSafe(cudaMemcpy(d_X, X, nmo * nao * sizeof(double), cudaMemcpyHostToDevice));

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

    cudaSafe(cudaMemcpy(errvec, d_C, nmo * nmo * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_B));
    cudaSafe(cudaFree(d_C));
}

auto
transformMatrix(double* transformed_F, const double* X, const double* F,
                const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::transformMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    cudaSafe(cudaSetDevice(0));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_F, *d_X, *d_Y;

    cudaSafe(cudaMalloc(&d_F, nao * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_X, nmo * nao * sizeof(double)));
    cudaSafe(cudaMalloc(&d_Y, nmo * nao * sizeof(double)));

    cudaSafe(cudaMemcpy(d_F, F, nao * nao * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_X, X, nmo * nao * sizeof(double), cudaMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

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

    cudaSafe(cudaMemcpy(transformed_F, d_F, nmo * nmo * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_F));
    cudaSafe(cudaFree(d_X));
    cudaSafe(cudaFree(d_Y));
}

auto
computeMatrixMultiplication(double* C, const double* A, const double* B, const std::string& trans_A, const std::string& trans_B,
                            const int64_t m_inp, const int64_t k_inp, const int64_t n_inp) -> void
{
    // TODO allow computeMatrixMultiplication on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeMatrixMultiplication: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    // TODO: matmul on multiple GPUs

    cudaSafe(cudaSetDevice(0));

    auto m = static_cast<int32_t>(m_inp);
    auto k = static_cast<int32_t>(k_inp);
    auto n = static_cast<int32_t>(n_inp);

    double *d_A, *d_B, *d_C;

    cudaSafe(cudaMalloc(&d_A, m * k * sizeof(double)));
    cudaSafe(cudaMalloc(&d_B, k * n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_C, m * n * sizeof(double)));

    cudaSafe(cudaMemcpy(d_A, A, m * k * sizeof(double), cudaMemcpyHostToDevice));
    cudaSafe(cudaMemcpy(d_B, B, k * n * sizeof(double), cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    cublasSafe(cublasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto op_A = (trans_A == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? CUBLAS_OP_N : CUBLAS_OP_T;

    // Note: we compute C^T = B^T * A^T since cublas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

    cublasSafe(cublasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));

    cudaSafe(cudaMemcpy(C, d_C, m * n * sizeof(double), cudaMemcpyDeviceToHost));

    cublasSafe(cublasDestroy(handle));

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_B));
    cudaSafe(cudaFree(d_C));
}

auto
diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void
{
    cudaSafe(cudaSetDevice(0));

    auto n = static_cast<int32_t>(nrows_A);
    int32_t lwork, info;

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    cudaSafe(cudaMalloc(&d_A, n * n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_D, n * sizeof(double)));
    cudaSafe(cudaMalloc(&d_info, sizeof(int32_t)));

    cudaSafe(cudaMemcpy(d_A, A, n * n * sizeof(double), cudaMemcpyHostToDevice));

    cusolverDnHandle_t handle;
    cusolverSafe(cusolverDnCreate(&handle));

    cusolverSafe(cusolverDnDsyevd_bufferSize(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    cudaSafe(cudaMalloc(&d_work, lwork * sizeof(double)));

    cusolverSafe(cusolverDnDsyevd(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    cusolverSafe(cusolverDnDestroy(handle));

    cudaSafe(cudaMemcpy(A, d_A, n * n * sizeof(double), cudaMemcpyDeviceToHost));
    cudaSafe(cudaMemcpy(D, d_D, n * sizeof(double), cudaMemcpyDeviceToHost));
    cudaSafe(cudaMemcpy(&info, d_info, sizeof(int32_t), cudaMemcpyDeviceToHost));

    // TODO: check info

    cudaSafe(cudaFree(d_A));
    cudaSafe(cudaFree(d_D));
    cudaSafe(cudaFree(d_info));
    cudaSafe(cudaFree(d_work));
}

}  // namespace gpu
