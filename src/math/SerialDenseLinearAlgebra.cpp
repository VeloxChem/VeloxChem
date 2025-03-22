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

    if ((narow == 0) || (nbrow == 0)) return mat;

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

}  // namespace sdenblas
