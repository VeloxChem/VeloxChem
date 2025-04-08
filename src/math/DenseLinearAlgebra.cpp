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

#include "DenseLinearAlgebra.hpp"

#ifdef ENABLE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#include "ErrorHandler.hpp"

namespace denblas {  // denblas namespace

auto
multAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbrow, "denblas::multAB: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbcol);

    if ((narow == 0) || (nbcol == 0)) return mat;

    // compute matrix-matrix multiplication

    auto narow_int32 = static_cast<int32_t>(narow);
    auto nacol_int32 = static_cast<int32_t>(nacol);
    auto nbcol_int32 = static_cast<int32_t>(nbcol);

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                narow_int32,
                nbcol_int32,
                nacol_int32,
                1.0,
                matrixA.values(),
                nacol_int32,
                matrixB.values(),
                nbcol_int32,
                0.0,
                mat.values(),
                nbcol_int32);

    return mat;
}

auto
multABt(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbcol, "denblas::multABt: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbrow);

    if ((narow == 0) || (nbrow == 0)) return mat;

    // compute matrix-matrix multiplcation

    auto narow_int32 = static_cast<int32_t>(narow);
    auto nbrow_int32 = static_cast<int32_t>(nbrow);
    auto nacol_int32 = static_cast<int32_t>(nacol);
    auto nbcol_int32 = static_cast<int32_t>(nbcol);

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasTrans,
                narow_int32,
                nbrow_int32,
                nacol_int32,
                1.0,
                matrixA.values(),
                nacol_int32,
                matrixB.values(),
                nbcol_int32,
                0.0,
                mat.values(),
                nbrow_int32);

    return mat;
}

auto
multAtB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(narow == nbrow, "denblas::multAtB: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(nacol, nbcol);

    if ((nacol == 0) || (nbcol == 0)) return mat;

    // compute matrix-matrix multiplication

    auto narow_int32 = static_cast<int32_t>(narow);
    auto nacol_int32 = static_cast<int32_t>(nacol);
    auto nbcol_int32 = static_cast<int32_t>(nbcol);

    cblas_dgemm(CblasRowMajor,
                CblasTrans,
                CblasNoTrans,
                nacol_int32,
                nbcol_int32,
                narow_int32,
                1.0,
                matrixA.values(),
                nacol_int32,
                matrixB.values(),
                nbcol_int32,
                0.0,
                mat.values(),
                nbcol_int32);

    return mat;
}

auto
subAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    errors::assertMsgCritical(matrixA.getNumberOfElements() == matrixB.getNumberOfElements(),
                              "denblas::subAB: Inconsistent sizes in matrix subtraction");

    // copy matrix

    CDenseMatrix mat = matrixA;

    // substract matrix

    auto nelems_int32 = static_cast<int32_t>(mat.getNumberOfElements());

    cblas_daxpy(nelems_int32, -1.0, matrixB.values(), 1, mat.values(), 1);

    return mat;
}

auto
addAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor) -> CDenseMatrix
{
    errors::assertMsgCritical(matrixA.getNumberOfElements() == matrixB.getNumberOfElements(),
                              "denblas::addAB: Inconsistent sizes in matrix addition");

    // copy matrix

    CDenseMatrix mat = matrixA;

    // add scaled matrix

    auto nelems_int32 = static_cast<int32_t>(mat.getNumberOfElements());

    cblas_daxpy(nelems_int32, factor, matrixB.values(), 1, mat.values(), 1);

    return mat;
}

auto
dot(const std::vector<double>& vectorA, const std::vector<double>& vectorB) -> double
{
    errors::assertMsgCritical(vectorA.size() == vectorB.size(), "denblas::dot: Inconsistent sizes in dot product of vectors");

    return cblas_ddot(vectorA.size(), vectorA.data(), 1, vectorB.data(), 1);
}

auto
trace(const CDenseMatrix& matrix) -> double
{
    errors::assertMsgCritical(matrix.getNumberOfColumns() == matrix.getNumberOfRows(), "denblas::trace: Expecting a square matrix");

    auto pvals = matrix.values();

    auto mdim = matrix.getNumberOfRows();

    double fsum = 0.0;

    for (int i = 0; i < mdim; i++)
    {
        fsum += pvals[i * mdim + i];
    }

    return fsum;
}

}  // namespace denblas
