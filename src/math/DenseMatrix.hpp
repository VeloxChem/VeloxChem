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

#ifndef DenseMatrix_hpp
#define DenseMatrix_hpp

#include <cstdint>
#include <vector>

#include "SubMatrix.hpp"

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
    int64_t _nRows;

    /**
     The number of columns.
     */
    int64_t _nColumns;

    /**
     The matrix element values.
     */
    CSubMatrix _values;

   public:
    /**
     Creates an empty dense matrix object.
     */
    CDenseMatrix();

    /**
     Creates a dense matrix object.

     @param values the SubMatrix object to copy data from.
     */
    CDenseMatrix(const CSubMatrix& submat);

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
    auto operator=(const CDenseMatrix& source) -> CDenseMatrix&;

    /**
     Assigns a dense matrix object by moving other dense matrix object.

     @param source the dense matrix object.
     */
    auto operator=(CDenseMatrix&& source) noexcept -> CDenseMatrix&;

    /**
     Compares dense matrix object with other dense matrix object.

     @param other the dense matrix object.
     @return true if dense matrix objects are equal, false otherwise.
     */
    auto operator==(const CDenseMatrix& other) const -> bool;

    /**
     Compares dense matrix object with other dense matrix object.

     @param other the dense matrix object.
     @return true if dense matrix objects are not equal, false otherwise.
     */
    auto operator!=(const CDenseMatrix& other) const -> bool;

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
    auto getNumberOfRows() const -> int64_t;

    /**
     Gets number of columns in dense matrix.

     @return the number of columns.
     */
    auto getNumberOfColumns() const -> int64_t;

    /**
     Gets number of elements in dense matrix.

     @return the number of elements.
     */
    auto getNumberOfElements() const -> int64_t;

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
    auto row(const int64_t iRow) const -> const double*;

    /**
     Gets pointer to first element of specific row in dense matrix.

     @return the pointer to first element of specific row.
     */
    auto row(const int64_t iRow) -> double*;

    /**
     Creates dense matrix object by slicing part of this dense matrix
     object.

     @param iPosition the position of first column.
     @param nElements the number of columns to be sliced.
     @return the dense matrix object.
     */
    auto slice(const int64_t iPosition, const int64_t nElements) const -> CDenseMatrix;
};

#endif /* DenseMatrix_hpp */
