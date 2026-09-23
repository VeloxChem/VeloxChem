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

}  // namespace tdenblas
