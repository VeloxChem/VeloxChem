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

#include "SerialDenseLinearAlgebra.hpp"

#include "Eigen/Dense"

#include "ErrorHandler.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#endif

/// @note The matrices here are stored by rows, and the math library works by
/// columns, so a matrix is seen by it as its own transpose. A row major
/// C = op(A) * op(B) is therefore asked for as the column major
/// C^T = op(B)^T * op(A)^T, with the operands handed over in the other order.
///
/// @note Which of the equivalent spellings is used is not a matter of taste.
/// Measured on one core with the shapes the exchange-correlation quadrature
/// produces, a contracted basis of 248 functions on a box of 225 points, the
/// library reaches 440 Gflop/s against 58 for the expression templates, and on
/// fourteen cores 112 against 36. But the plain A * B of a square matrix with a
/// wide one, asked for in that order, collapses to 25 -- below the expression
/// templates -- while the identical product with the square operand marked as
/// transposed holds 92. The square operand here is a density matrix, which is
/// symmetric, so marking it transposed costs nothing and buys the difference.

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

    if ((narow == 0) || (nbrow == 0)) return mat;

    // compute matrix-matrix multiplcation

    auto A = matrixA.values();
    auto B = matrixB.values();

    mat.zero();
    auto C = mat.values();

#ifdef VLX_USE_MATHLIB

    // C^T = B * A^T, over the columns

    const lapack_int_t mdim = static_cast<lapack_int_t>(nbrow);
    const lapack_int_t ndim = static_cast<lapack_int_t>(narow);
    const lapack_int_t kdim = static_cast<lapack_int_t>(nacol);

    const double alpha = 1.0, beta = 0.0;

    dgemm_("T", "N", &mdim, &ndim, &kdim, &alpha, B, &kdim, A, &kdim, &beta, C, &mdim);

#else

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, narow, nbrow);

    ematC.noalias() = ematA * ematB.transpose();

#endif

    return mat;
}

auto
serialMultSymAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();
    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();
    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(narow == nacol, "sdenblas::serialMultSymAB: Matrix A is not square");

    errors::assertMsgCritical(nacol == nbrow, "sdenblas::serialMultSymAB: Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbcol);

    if ((narow == 0) || (nbcol == 0)) return mat;

    // compute matrix-matrix multiplication

    auto A = matrixA.values();
    auto B = matrixB.values();

    mat.zero();
    auto C = mat.values();

#ifdef VLX_USE_MATHLIB

    // C^T = B^T * A, over the columns, with A marked as transposed: it is
    // symmetric, so this is the same matrix, and the shape the library is fast at

    const lapack_int_t mdim = static_cast<lapack_int_t>(nbcol);
    const lapack_int_t ndim = static_cast<lapack_int_t>(narow);
    const lapack_int_t kdim = static_cast<lapack_int_t>(nacol);

    const double alpha = 1.0, beta = 0.0;

    dgemm_("N", "T", &mdim, &ndim, &kdim, &alpha, B, &mdim, A, &ndim, &beta, C, &mdim);

#else

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematA(A, narow, nacol);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematB(B, nbrow, nbcol);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned> ematC(C, narow, nbcol);

    ematC.noalias() = ematA * ematB;

#endif

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


}  // namespace sdenblas
