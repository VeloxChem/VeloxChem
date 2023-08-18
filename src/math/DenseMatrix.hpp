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

#include "SubMatrix.hpp"
#include "T4Index.hpp"

/**
 Class CDenseMatrix stores dense matrix in coordinate format (zero-based
 indexing scheme) and provides set of methods for manipulating dense matrix
 data.

 @author Z. Rinkevicius
 */
class CDenseMatrix
{
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

     @param values the vector of matrix elements.
     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     */
    CDenseMatrix(const std::vector<double>& values, const int32_t nRows, const int32_t nColumns);

    /**
     Creates a dense matrix object.

     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     */
    CDenseMatrix(const int32_t nRows, const int32_t nColumns);

    /**
     Creates an dense matrix object.

     @param other the dense matrix to copy from.
     */
    CDenseMatrix(const CDenseMatrix& other);

    /**
     Destroys a dense matrix object.
     */
    ~CDenseMatrix();

    /**
     Sets all values in dense matrix to zero.
     */
    void zero();

    /**
     Creates transpose dense matrix.

     @return the transpose dense matrix.
     */
    CDenseMatrix transpose() const;

    /**
     Symmetrizes elements of square matrix: a_ij = a_ji = (a_ij + a_ji).
     */
    void symmetrize();

    /**
     Symmetrizes elements of square matrix: a_ij = a_ji = factor * (a_ij + a_ji).

     @param factor the factor.
     */
    void symmetrizeAndScale(const double factor);

    /**
     Gets number of rows in dense matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in dense matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in dense matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of dense matrix.

     @return the constant pointer to first element of dense matrix.
     */
    const double* values() const;

    /**
     Gets pointer to first element of dense matrix.

     @return the pointer to first element of dense matrix.
     */
    double* values();

    /**
     Gets constant pointer to first element of specific row in dense matrix.

     @return the constant pointer to first element of specific row.
     */
    const double* row(const int32_t iRow) const;

    /**
     Gets pointer to first element of specific row in dense matrix.

     @return the pointer to first element of specific row.
     */
    double* row(const int32_t iRow);
};

#endif /* DenseMatrix_hpp */
