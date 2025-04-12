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
