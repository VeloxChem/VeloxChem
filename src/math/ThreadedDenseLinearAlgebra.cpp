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

#include "ThreadedDenseLinearAlgebra.hpp"

#include <omp.h>

#include "ErrorHandler.hpp"
#include "SerialDenseLinearAlgebra.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#endif

namespace tdenblas {  // tdenblas namespace

auto
threadedMultAtB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void
{
#ifdef VLX_USE_MATHLIB

    // use the math library only outside OpenMP parallel regions, and fall
    // back to the serial implementation otherwise

    if (!omp_in_parallel())
    {
        auto narow = matrixA.number_of_rows();
        auto nacol = matrixA.number_of_columns();

        auto nbrow = matrixB.number_of_rows();
        auto nbcol = matrixB.number_of_columns();

        auto ncrow = matrixC.number_of_rows();
        auto nccol = matrixC.number_of_columns();

        errors::assertMsgCritical(narow == nbrow, "tdenblas::threadedMultAtB: Inconsistent sizes in matrix multiplication");

        errors::assertMsgCritical(ncrow == nacol, "tdenblas::threadedMultAtB: Inconsistent sizes in matrix multiplication");

        errors::assertMsgCritical(nccol == nbcol, "tdenblas::threadedMultAtB: Inconsistent sizes in matrix multiplication");

        if ((narow == 0) || (nacol == 0) || (nbcol == 0)) return;

        auto A = matrixA.data();
        auto B = matrixB.data();
        auto C = matrixC.data();

        // C^T = B^T * A, over the columns

        const lapack_int_t mdim = static_cast<lapack_int_t>(nbcol);
        const lapack_int_t ndim = static_cast<lapack_int_t>(nacol);
        const lapack_int_t kdim = static_cast<lapack_int_t>(narow);

        const double alpha = 1.0, beta = 1.0;

        dgemm_("N", "T", &mdim, &ndim, &kdim, &alpha, B, &mdim, A, &ndim, &beta, C, &mdim);

        return;
    }

#endif

    sdenblas::serialMultAtB(matrixC, matrixA, matrixB);
}

auto
threadedMultAB(const size_t  nrows,
               const size_t  ncols,
               const size_t  nsums,
               const double  alpha,
               const double *matrixA,
               const size_t  lda,
               const double *matrixB,
               const size_t  ldb,
               const double  beta,
               double       *matrixC,
               const size_t  ldc) -> void
{
#ifdef VLX_USE_MATHLIB

    // use the math library only outside OpenMP parallel regions, and fall back
    // to the serial implementation otherwise

    if (!omp_in_parallel())
    {
        // NOTE: the library is column major and the column major matrix of a row
        // major array is its transpose, so the product of the row major arrays is
        // the product of the two in the other order with the rows and the columns
        // swapped.

        const char trans = 'N';

        auto m_arg = static_cast<lapack_int_t>(ncols);

        auto n_arg = static_cast<lapack_int_t>(nrows);

        auto k_arg = static_cast<lapack_int_t>(nsums);

        auto lda_arg = static_cast<lapack_int_t>(lda);

        auto ldb_arg = static_cast<lapack_int_t>(ldb);

        auto ldc_arg = static_cast<lapack_int_t>(ldc);

        dgemm_(&trans, &trans, &m_arg, &n_arg, &k_arg, &alpha, matrixB, &ldb_arg, matrixA, &lda_arg, &beta, matrixC,
               &ldc_arg);

        return;
    }

#endif

    sdenblas::serialMultAB(nrows, ncols, nsums, alpha, matrixA, lda, matrixB, ldb, beta, matrixC, ldc);
}

auto
threadedMultABt(const size_t  nrows,
                const size_t  ncols,
                const size_t  nsums,
                const double  alpha,
                const double *matrixA,
                const size_t  lda,
                const double *matrixB,
                const size_t  ldb,
                const double  beta,
                double       *matrixC,
                const size_t  ldc) -> void
{
#ifdef VLX_USE_MATHLIB

    // use the math library only outside OpenMP parallel regions, and fall back
    // to the serial implementation otherwise

    if (!omp_in_parallel())
    {
        // NOTE: the library is column major and the column major matrix of a row
        // major array is its transpose, so the product of the row major arrays is
        // the product of the two in the other order with the rows and the columns
        // swapped.

        const char trans_t = 'T';

        const char trans_n = 'N';

        auto m_arg = static_cast<lapack_int_t>(ncols);

        auto n_arg = static_cast<lapack_int_t>(nrows);

        auto k_arg = static_cast<lapack_int_t>(nsums);

        auto lda_arg = static_cast<lapack_int_t>(lda);

        auto ldb_arg = static_cast<lapack_int_t>(ldb);

        auto ldc_arg = static_cast<lapack_int_t>(ldc);

        dgemm_(&trans_t, &trans_n, &m_arg, &n_arg, &k_arg, &alpha, matrixB, &ldb_arg, matrixA, &lda_arg, &beta,
               matrixC, &ldc_arg);

        return;
    }

#endif

    sdenblas::serialMultABt(nrows, ncols, nsums, alpha, matrixA, lda, matrixB, ldb, beta, matrixC, ldc);
}

auto
threadedRankUpdate(const size_t  n,
                   const size_t  k,
                   const double  alpha,
                   const double *matrixA,
                   const size_t  lda,
                   double       *matrixC,
                   const size_t  ldc) -> void
{
#ifdef VLX_USE_MATHLIB

    // use the math library only outside OpenMP parallel regions, and fall back
    // to the serial implementation otherwise

    if (!omp_in_parallel())
    {
        // NOTE: the library is column major and the column major matrix of a row
        // major array is its transpose, so the upper triangle of the library is
        // the lower triangle of the array. The update of the transposed stored
        // matrix is therefore the update of the row major array.

        const char uplo = 'U';

        const char trans = 'T';

        auto n_arg = static_cast<lapack_int_t>(n);

        auto k_arg = static_cast<lapack_int_t>(k);

        auto lda_arg = static_cast<lapack_int_t>(lda);

        auto ldc_arg = static_cast<lapack_int_t>(ldc);

        const double one = 1.0;

        dsyrk_(&uplo, &trans, &n_arg, &k_arg, &alpha, matrixA, &lda_arg, &one, matrixC, &ldc_arg);

        return;
    }

#endif

    sdenblas::serialRankUpdate(n, k, alpha, matrixA, lda, matrixC, ldc);
}

}  // namespace tdenblas
