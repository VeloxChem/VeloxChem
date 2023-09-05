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

#ifndef DenseMatrix_hpp
#define DenseMatrix_hpp

#include <cstdint>
#include <vector>

/**
 Class CDenseMatrix stores dense matrix in coordinate format (zero-based
 indexing scheme) and provides set of methods for manipulating dense matrix
 data.

 @author Z. Rinkevicius
 */
class CDenseMatrix
{
    /**
     The number of rows.
     */
    int32_t _nRows;

    /**
     The number of columns.
     */
    int32_t _nColumns;

    /**
     The matrix element values.
     */
    std::vector<double> _values;

   public:
    /**
     Creates an empty dense matrix object.
     */
    CDenseMatrix();

    /**
     Creates a dense matrix object.

     @param values the vector of matrix elements.
     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     */
    CDenseMatrix(const std::vector<double>& values, const int64_t nRows, const int64_t nColumns);

    /**
     Creates a dense matrix object.

     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     */
    CDenseMatrix(const int64_t nRows, const int64_t nColumns);

    /**
     Creates a dense matrix object by copying other dense matrix object.

     @param source the dense matrix object.
     */
    CDenseMatrix(const CDenseMatrix& source);

    /**
     Creates a dense matrix object by moving other dense matrix object.

     @param source the dense matrix object.
     */
    CDenseMatrix(CDenseMatrix&& source) noexcept;

    /**
     Destroys a dense matrix object.
     */
    ~CDenseMatrix();

    /**
     Assigns a dense matrix object by copying other dense matrix object.

     @param source the dense matrix object.
     */
    CDenseMatrix& operator=(const CDenseMatrix& source);

    /**
     Assigns a dense matrix object by moving other dense matrix object.

     @param source the dense matrix object.
     */
    CDenseMatrix& operator=(CDenseMatrix&& source) noexcept;

    /**
     Sets all values in dense matrix to zero.
     */
    auto zero() -> void;

    /**
     Creates transpose dense matrix.

     @return the transpose dense matrix.
     */
    auto transpose() const -> CDenseMatrix;

    /**
     Symmetrizes elements of square matrix: a_ij = a_ji = (a_ij + a_ji).
     */
    auto symmetrize() -> void;

    /**
     Symmetrizes elements of square matrix: a_ij = a_ji = factor * (a_ij + a_ji).

     @param factor the factor.
     */
    auto symmetrizeAndScale(const double factor) -> void;

    /**
     Gets number of rows in dense matrix.

     @return the number of rows.
     */
    auto getNumberOfRows() const -> int32_t;

    /**
     Gets number of columns in dense matrix.

     @return the number of columns.
     */
    auto getNumberOfColumns() const -> int32_t;

    /**
     Gets number of elements in dense matrix.

     @return the number of elements.
     */
    auto getNumberOfElements() const -> int32_t;

    /**
     Gets constant pointer to first element of dense matrix.

     @return the constant pointer to first element of dense matrix.
     */
    auto values() const -> const double*;

    /**
     Gets pointer to first element of dense matrix.

     @return the pointer to first element of dense matrix.
     */
    auto values() -> double*;

    /**
     Gets constant pointer to first element of specific row in dense matrix.

     @return the constant pointer to first element of specific row.
     */
    auto row(const int32_t iRow) const -> const double*;

    /**
     Gets pointer to first element of specific row in dense matrix.

     @return the pointer to first element of specific row.
     */
    auto row(const int32_t iRow) -> double*;
};

#endif /* DenseMatrix_hpp */
