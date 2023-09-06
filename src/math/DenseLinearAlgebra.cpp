//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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
multDiagByA(const std::vector<double>& diagonal, const CDenseMatrix& matrix) -> CDenseMatrix
{
    // set up dimensions of matrix

    auto nrow = matrix.getNumberOfRows();
    auto ncol = matrix.getNumberOfColumns();

    // allocate results matrix

    CDenseMatrix mat(nrow, ncol);

    // set up pointers to matrices

    auto mval = mat.values();
    auto sval = matrix.values();

    auto dval = diagonal.data();

    // compute matrix multiplication

    for (int64_t i = 0; i < nrow; i++)
    {
        // set up local pointers to rows

        auto cmval = mval + i * ncol;
        auto csval = sval + i * ncol;

        // fetch value of diagonal

        auto f = dval[i];

        for (int64_t j = 0; j < ncol; j++)
        {
            cmval[j] = f * csval[j];
        }
    }

    return mat;
}

auto
multDiagByAt(const std::vector<double>& diagonal, const CDenseMatrix& matrix) -> CDenseMatrix
{
    // set up dimensions of matrix

    auto nrow = matrix.getNumberOfRows();
    auto ncol = matrix.getNumberOfColumns();

    // allocate results matrix

    CDenseMatrix mat(ncol, nrow);

    // set up pointers to matrices

    auto mval = mat.values();
    auto sval = matrix.values();

    auto dval = diagonal.data();

    // compute matrix multiplication

    for (int64_t i = 0; i < ncol; i++)
    {
        auto ioff = i * nrow;

        for (int64_t j = 0; j < nrow; j++)
        {
            mval[ioff + j] = dval[i] * sval[j * ncol + i];
        }
    }

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

    for (int64_t i = 0; i < mdim; i++)
    {
        fsum += pvals[i * mdim + i];
    }

    return fsum;
}

}  // namespace denblas
