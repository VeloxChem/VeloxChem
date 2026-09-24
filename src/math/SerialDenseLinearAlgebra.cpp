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

// NOTE: the serial routines must not start threads of their own, as their callers
// parallelize above them. Eigen parallelizes with OpenMP when this file is
// compiled with it, so its parallelizer is disabled here.

#define EIGEN_DONT_PARALLELIZE

#include "SerialDenseLinearAlgebra.hpp"

#include "Eigen/Dense"

#include "ErrorHandler.hpp"

namespace sdenblas {  // sdenblas namespace

auto
serialMultAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbrow, "sdenblas::serialMultAB: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbcol);

    if ((narow == 0) || (nbcol == 0)) return mat;

    // compute matrix-matrix multiplication

    auto A = matrixA.values();
    auto B = matrixB.values();

    mat.zero();
    auto C = mat.values();

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, narow, nbcol);

    ematC.noalias() = ematA * ematB;

    return mat;
}

auto
serialMultABt(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbcol, "sdenblas::serialMultABt: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbrow);

    if ((narow == 0) || (nbrow == 0) || (nacol == 0)) return mat;

    // compute matrix-matrix multiplcation

    auto A = matrixA.values();
    auto B = matrixB.values();

    mat.zero();
    auto C = mat.values();

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, narow, nbrow);

    ematC.noalias() = ematA * ematB.transpose();

    return mat;
}

auto
serialMultAtB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(narow == nbrow, "sdenblas::serialMultAtB: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(nacol, nbcol);

    if ((narow == 0) || (nbcol == 0)) return mat;

    // compute matrix-matrix multiplication

    auto A = matrixA.values();
    auto B = matrixB.values();

    mat.zero();
    auto C = mat.values();

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, nacol, nbcol);

    ematC.noalias() = ematA.transpose() * ematB;

    return mat;
}

auto
serialAddAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor) -> CDenseMatrix
{
    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical((narow == nbrow) && (nacol == nbcol),
                              "sdenblas::serialAddAB: Inconsistent sizes in matrix addition");

    auto A = matrixA.values();
    auto B = matrixB.values();

    CDenseMatrix mat(narow, nacol);
    mat.zero();
    auto C = mat.values();

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, narow, nacol);

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, narow, nacol);

    ematC = ematA + factor * ematB;

    return mat;
}

auto
serialInPlaceAddAB(CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor) -> void
{
    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical((narow == nbrow) && (nacol == nbcol),
                              "sdenblas::serialInPlaceAddAB: Inconsistent sizes in matrix addition");

    auto A = matrixA.values();
    auto B = matrixB.values();

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, narow, nacol);

    ematA += factor * ematB;
}

auto
serialSolve(const CDenseMatrix& mat, const std::vector<double>& vec) -> std::vector<double>
{
    auto narow = mat.getNumberOfRows();
    auto nacol = mat.getNumberOfColumns();

    auto nbsize = static_cast<int>(vec.size());

    errors::assertMsgCritical((narow == nacol) && (narow == nbsize),
                              "sdenblas::serialSolve: Inconsistent sizes in matrix and vector");

    std::vector<double> sol(nbsize);

    auto A = mat.values();
    auto B = vec.data();
    auto C = sol.data();

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned> evecB(B, nbsize);
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned> evecC(C, nbsize);

    evecC = ematA.colPivHouseholderQr().solve(evecB);

    return sol;
}

auto
serialMultAB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void
{
    // set up dimensions of matrix A

    auto narow = matrixA.number_of_rows();
    
    auto nacol = matrixA.number_of_columns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.number_of_rows();
    
    auto nbcol = matrixB.number_of_columns();
    
    errors::assertMsgCritical(nacol == nbrow, "sdenblas::serialMultAB: Inconsistent sizes in matrix multiplication");
    
    // set up dimensions of matrix C

    auto ncrow = matrixC.number_of_rows();
    
    auto nccol = matrixC.number_of_columns();
    
    errors::assertMsgCritical(ncrow == narow, "sdenblas::serialMultAB: Inconsistent sizes in matrix multiplication");
    
    errors::assertMsgCritical(nccol == nbcol, "sdenblas::serialMultAB: Inconsistent sizes in matrix multiplication");
    
    // set up pointers to data

    auto A = matrixA.data();
    
    auto B = matrixB.data();
    
    auto C = matrixC.data();

    // compute matrix-matrix multiplication
    
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, ncrow, nccol);
    
    ematC.noalias() += ematA * ematB;
}

auto
serialMultAtB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void
{
    // set up dimensions of matrix A

    auto narow = matrixA.number_of_rows();
    
    auto nacol = matrixA.number_of_columns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.number_of_rows();
    
    auto nbcol = matrixB.number_of_columns();
    
    errors::assertMsgCritical(narow == nbrow, "sdenblas::serialMultAtB: Inconsistent sizes in matrix multiplication");
    
    // set up dimensions of matrix C

    auto ncrow = matrixC.number_of_rows();
    
    auto nccol = matrixC.number_of_columns();
    
    errors::assertMsgCritical(ncrow == nacol, "sdenblas::serialMultAtB: Inconsistent sizes in matrix multiplication");
    
    errors::assertMsgCritical(nccol == nbcol, "sdenblas::serialMultAtB: Inconsistent sizes in matrix multiplication");
    
    // set up pointers to data

    auto A = matrixA.data();
    
    auto B = matrixB.data();
    
    auto C = matrixC.data();

    // compute matrix-matrix multiplication

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, ncrow, nccol);
    
    ematC.noalias() += ematA.transpose() * ematB;
}

auto
serialMultAB(const size_t  nrows,
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
    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using RowMajorStride = Eigen::Stride<Eigen::Dynamic, 1>;

    const auto rows = static_cast<Eigen::Index>(nrows);

    const auto cols = static_cast<Eigen::Index>(ncols);

    const auto sums = static_cast<Eigen::Index>(nsums);

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> ematA(matrixA, rows, sums,
                                                              RowMajorStride(static_cast<Eigen::Index>(lda), 1));

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> ematB(matrixB, sums, cols,
                                                              RowMajorStride(static_cast<Eigen::Index>(ldb), 1));

    Eigen::Map<RowMajorMatrix, 0, RowMajorStride> ematC(matrixC, rows, cols,
                                                        RowMajorStride(static_cast<Eigen::Index>(ldc), 1));

    if (beta == 0.0)
    {
        ematC.noalias() = alpha * ematA * ematB;
    }
    else
    {
        if (beta != 1.0) ematC *= beta;

        ematC.noalias() += alpha * ematA * ematB;
    }
}

auto
serialMultABt(const size_t  nrows,
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
    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using RowMajorStride = Eigen::Stride<Eigen::Dynamic, 1>;

    const auto rows = static_cast<Eigen::Index>(nrows);

    const auto cols = static_cast<Eigen::Index>(ncols);

    const auto sums = static_cast<Eigen::Index>(nsums);

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> ematA(matrixA, rows, sums,
                                                              RowMajorStride(static_cast<Eigen::Index>(lda), 1));

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> ematB(matrixB, cols, sums,
                                                              RowMajorStride(static_cast<Eigen::Index>(ldb), 1));

    Eigen::Map<RowMajorMatrix, 0, RowMajorStride> ematC(matrixC, rows, cols,
                                                        RowMajorStride(static_cast<Eigen::Index>(ldc), 1));

    if (beta == 0.0)
    {
        ematC.noalias() = alpha * (ematA * ematB.transpose());
    }
    else
    {
        if (beta != 1.0) ematC *= beta;

        ematC.noalias() += alpha * (ematA * ematB.transpose());
    }
}

auto
serialRankUpdate(const size_t  n,
                 const size_t  k,
                 const double  alpha,
                 const double *matrixA,
                 const size_t  lda,
                 double       *matrixC,
                 const size_t  ldc) -> void
{
    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using RowMajorStride = Eigen::Stride<Eigen::Dynamic, 1>;

    const auto rows = static_cast<Eigen::Index>(n);

    const auto sums = static_cast<Eigen::Index>(k);

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> ematA(matrixA, rows, sums,
                                                              RowMajorStride(static_cast<Eigen::Index>(lda), 1));

    Eigen::Map<RowMajorMatrix, 0, RowMajorStride> ematC(matrixC, rows, rows,
                                                        RowMajorStride(static_cast<Eigen::Index>(ldc), 1));

    ematC.template selfadjointView<Eigen::Lower>().rankUpdate(ematA, alpha);
}

auto
serialSolveTriangular(const size_t  nrows,
                      const size_t  ncols,
                      const double *factor,
                      const size_t  ldf,
                      double       *values,
                      const size_t  ldv,
                      const bool    transposed) -> void
{
    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using RowMajorStride = Eigen::Stride<Eigen::Dynamic, 1>;

    const auto rows = static_cast<Eigen::Index>(nrows);

    const auto cols = static_cast<Eigen::Index>(ncols);

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> ematL(factor, rows, rows,
                                                             RowMajorStride(static_cast<Eigen::Index>(ldf), 1));

    Eigen::Map<RowMajorMatrix, 0, RowMajorStride> ematB(values, rows, cols,
                                                        RowMajorStride(static_cast<Eigen::Index>(ldv), 1));

    if (transposed)
    {
        ematL.transpose().template triangularView<Eigen::Upper>().solveInPlace(ematB);
    }
    else
    {
        ematL.template triangularView<Eigen::Lower>().solveInPlace(ematB);
    }
}


}  // namespace sdenblas
