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

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
//#include <hipsolver/hipsolver.h>
#include <magma_v2.h>

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

    hipSafe(hipSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_A, *d_B;

    hipSafe(hipMalloc(&d_A, n * sizeof(double)));
    hipSafe(hipMalloc(&d_B, n * sizeof(double)));

    hipSafe(hipMemcpy(d_A, A, n * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, B, n * sizeof(double), hipMemcpyHostToDevice));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double dot_product;
    hipblasSafe(hipblasDdot(handle, n, d_A, 1, d_B, 1, &dot_product));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));

    return dot_product;
}

auto
computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& data_pointers, const int64_t size) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeWeightedSum: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    auto n = static_cast<int32_t>(size);

    double *d_X, *d_Y;

    hipSafe(hipMalloc(&d_X, n * sizeof(double)));
    hipSafe(hipMalloc(&d_Y, n * sizeof(double)));

    hipSafe(hipMemcpy(d_Y, weighted_data, n * sizeof(double), hipMemcpyHostToDevice));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    for (size_t i = 0; i < data_pointers.size(); i++)
    {
        double alpha = weights[i];

        hipSafe(hipMemcpy(d_X, data_pointers[i], n * sizeof(double), hipMemcpyHostToDevice));

        hipblasSafe(hipblasDaxpy(handle, n, &alpha, d_X, 1, d_Y, 1));
    }

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipMemcpy(weighted_data, d_Y, n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_X));
    hipSafe(hipFree(d_Y));
}

auto
computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                   const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeErrorVector: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_A, *d_B, *d_C;

    hipSafe(hipMalloc(&d_A, nao * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_B, nao * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_C, nao * nao * sizeof(double)));

    hipSafe(hipMemcpy(d_A, F, nao * nao * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, D, nao * nao * sizeof(double), hipMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major

    // D^T F^T (=> FD)
    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_B, nao, d_A, nao, &beta, d_C, nao));

    // S^T(FD)^T (=> FDS)
    hipSafe(hipMemcpy(d_A, S, nao * nao * sizeof(double), hipMemcpyHostToDevice));

    hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, nao, nao, nao, &alpha, d_A, nao, d_C, nao, &beta, d_B, nao));

    hipSafe(hipDeviceSynchronize());

    // FDS - (FDS)^T
    beta = -1.0;

    hipblasSafe(hipblasDgeam(handle, HIPBLAS_OP_N, HIPBLAS_OP_T, nao, nao, &alpha, d_B, nao, &beta, d_B, nao, d_C, nao));

    // X^T (FDS - (FDS)^T) X
    double* d_X = d_A;  // note: nao >= nmo
    double* d_Y = d_B;  // note: nao >= nmo

    hipSafe(hipMemcpy(d_X, X, nmo * nao * sizeof(double), hipMemcpyHostToDevice));

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

    hipSafe(hipMemcpy(errvec, d_C, nmo * nmo * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));
    hipSafe(hipFree(d_C));
}

auto
transformMatrix(double* transformed_F, const double* X, const double* F,
                const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::transformMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    auto nmo = static_cast<int32_t>(nmo_inp);
    auto nao = static_cast<int32_t>(nao_inp);

    double *d_F, *d_X, *d_Y;

    hipSafe(hipMalloc(&d_F, nao * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_X, nmo * nao * sizeof(double)));
    hipSafe(hipMalloc(&d_Y, nmo * nao * sizeof(double)));

    hipSafe(hipMemcpy(d_F, F, nao * nao * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_X, X, nmo * nao * sizeof(double), hipMemcpyHostToDevice));

    double alpha = 1.0, beta = 0.0;

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

    hipSafe(hipMemcpy(transformed_F, d_F, nmo * nmo * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_F));
    hipSafe(hipFree(d_X));
    hipSafe(hipFree(d_Y));
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

    hipSafe(hipSetDevice(0));

    auto m = static_cast<int32_t>(m_inp);
    auto k = static_cast<int32_t>(k_inp);
    auto n = static_cast<int32_t>(n_inp);

    double *d_A, *d_B, *d_C;

    hipSafe(hipMalloc(&d_A, m * k * sizeof(double)));
    hipSafe(hipMalloc(&d_B, k * n * sizeof(double)));
    hipSafe(hipMalloc(&d_C, m * n * sizeof(double)));

    hipSafe(hipMemcpy(d_A, A, m * k * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, B, k * n * sizeof(double), hipMemcpyHostToDevice));

    hipblasHandle_t handle;
    hipblasSafe(hipblasCreate(&handle));

    double alpha = 1.0, beta = 0.0;

    auto op_A = (trans_A == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
    auto op_B = (trans_B == std::string("N")) ? HIPBLAS_OP_N : HIPBLAS_OP_T;

    // Note: we compute C^T = B^T * A^T since hipblas is column-major
    // Also need to adjust lda accordingly

    auto lda_A = (trans_A == std::string("N")) ? k : m;
    auto lda_B = (trans_B == std::string("N")) ? n : k;

    hipblasSafe(hipblasDgemm(handle, op_B, op_A, n, m, k, &alpha, d_B, lda_B, d_A, lda_A, &beta, d_C, n));

    hipSafe(hipMemcpy(C, d_C, m * n * sizeof(double), hipMemcpyDeviceToHost));

    hipblasSafe(hipblasDestroy(handle));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));
    hipSafe(hipFree(d_C));
}

auto
diagonalizeMatrix(double* A, double* W, const int64_t nrows_A) -> void
{
    // TODO allow diagonalizeMatrix on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::diagonalizeMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(0));

    magmaSafe(magma_init());

    magma_setdevice(0);

    magma_int_t n = static_cast<magma_int_t>(nrows_A);
    //magma_int_t ldda = magma_roundup(n, 32);

    double *d_A;
    //hipSafe(hipMalloc(&d_A, n * ldda * sizeof(double)));
    //hipblasSafe(hipblasSetMatrix(n, n, sizeof(double), A, n, d_A, ldda));
    hipSafe(hipMalloc(&d_A, n * n * sizeof(double)));
    hipSafe(hipMemcpy(d_A, A, n * n * sizeof(double), hipMemcpyHostToDevice));

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n * n * sizeof(double)));
    hipSafe(hipHostMalloc(&work, lwork * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, liwork * sizeof(magma_int_t)));

    magma_int_t info;
    //magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, ldda, W, wA, n, work, lwork, iwork, liwork, &info);
    magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, n, W, wA, n, work, lwork, iwork, liwork, &info);

    if (info != 0)
    {
        std::stringstream ss;
        ss << "gpu::diagonalizeMatrix: (magma error) " << magma_strerror(info);
        errors::assertMsgCritical(false, ss.str());
    }

    //hipblasSafe(hipblasGetMatrix(n, n, sizeof(double), d_A, ldda, A, n));
    hipSafe(hipMemcpy(A, d_A, n * n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));

    hipSafe(hipHostFree(wA));
    hipSafe(hipHostFree(work));
    hipSafe(hipHostFree(iwork));

    magmaSafe(magma_finalize());
}

auto
diagonalizeMatrixMultiGPU(double* A, double* W, const int64_t nrows_A, const int64_t num_gpus_per_node) -> void
{
    // TODO allow diagonalizeMatrix on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::diagonalizeMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    magmaSafe(magma_init());

    magma_int_t n = static_cast<magma_int_t>(nrows_A);

    magma_int_t ngpu = static_cast<magma_int_t>(num_gpus_per_node);

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n * n * sizeof(double)));
    hipSafe(hipHostMalloc(&work, lwork * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, liwork * sizeof(magma_int_t)));

    magma_int_t info;
    magma_dsyevd_m(ngpu, MagmaVec, MagmaUpper, n, A, n, W, work, lwork, iwork, liwork, &info);

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

/*
auto
diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void
{
    hipSafe(hipSetDevice(0));

    auto n = static_cast<int32_t>(nrows_A);
    int32_t lwork;

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    hipSafe(hipMalloc(&d_A, n * n * sizeof(double)));
    hipSafe(hipMalloc(&d_D, n * sizeof(double)));
    hipSafe(hipMalloc(&d_info, sizeof(int32_t)));

    hipSafe(hipMemcpy(d_A, A, n * n * sizeof(double), hipMemcpyHostToDevice));

    hipsolverHandle_t handle;
    hipsolverSafe(hipsolverCreate(&handle));

    hipsolverSafe(hipsolverDsyevd_bufferSize(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    hipSafe(hipMalloc(&d_work, lwork));

    hipsolverSafe(hipsolverDsyevd(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    hipsolverSafe(hipsolverDestroy(handle));

    // TODO check h_info
    int32_t h_info;

    hipSafe(hipMemcpy(&h_info, d_info, sizeof(int32_t), hipMemcpyDeviceToHost));
    hipSafe(hipMemcpy(A, d_A, n * n * sizeof(double), hipMemcpyDeviceToHost));
    hipSafe(hipMemcpy(D, d_D, n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_D));
    hipSafe(hipFree(d_info));
    hipSafe(hipFree(d_work));
}
*/

}  // namespace gpu
