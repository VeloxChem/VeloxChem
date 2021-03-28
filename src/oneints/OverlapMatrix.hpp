//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#ifndef OverlapMatrix_hpp
#define OverlapMatrix_hpp

#include <cstdint>
#include <string>

#include "DenseMatrix.hpp"

/**
 Class COverlapMatrix stores general one electron overlap matrix and provides
 set of methods for handling of overlap matrix data.

 @author Z. Rinkevicius
 */
class COverlapMatrix
{
    /**
     The generic dense overlap matrix (rectangular or square).
     */
    CDenseMatrix _matrix;

   public:
    /**
     Creates an empty overlap matrix object.
     */
    COverlapMatrix();

    /**
     Creates a overlap matrix object.

     @param matrix the dense matrix with overlap integrals.
     */
    COverlapMatrix(const CDenseMatrix& matrix);

    /**
     Creates a overlap matrix object by copying other overlap matrix object.

     @param source the overlap matrix object.
     */
    COverlapMatrix(const COverlapMatrix& source);

    /**
     Creates a overlap matrix object by moving other overlap matrix object.

     @param source the overlap matrix object.
     */
    COverlapMatrix(COverlapMatrix&& source) noexcept;

    /**
     Destroys a overlap matrix object.
     */
    ~COverlapMatrix();

    /**
     Assigns a overlap matrix object by copying other overlap matrix object.

     @param source the overlap matrix object.
     */
    COverlapMatrix& operator=(const COverlapMatrix& source);

    /**
     Assigns a overlap matrix object by moving other overlap matrix object.

     @param source the overlap matrix object.
     */
    COverlapMatrix& operator=(COverlapMatrix&& source) noexcept;

    /**
     Compares overlap matrix object with other overlap matrix object.

     @param other the overlap matrix object.
     @return true if overlap matrix objects are equal, false otherwise.
     */
    bool operator==(const COverlapMatrix& other) const;

    /**
     Compares overlap matrix object with other overlap matrix object.

     @param other the overlap matrix object.
     @return true if overlap matrix objects are not equal, false otherwise.
     */
    bool operator!=(const COverlapMatrix& other) const;

    /**
     Gets string representation of overlap matrix.

     @return a string for printing the overlap matrix.
     */
    std::string getString() const;

    /**
     Gets number of rows in overlap matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in overlap matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in overlap matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of overlap matrix.

     @return the constant pointer to first element of overlap matrix.
     */
    const double* values() const;

    /**
     Gets Löwdin symmetric orthogonalization matrix if overalp matrix
     eigenvalues is linearly independent, or canonical orthogonalization
     matrix if overlap matrix eigenvalues is linearly dependent.

     @param threshold the linear dependence threshold.
     @return the orthogonalization matrix.
     */
    CDenseMatrix getOrthogonalizationMatrix(const double threshold) const;

    /**
     Converts overlap matrix object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the overlap matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const COverlapMatrix& source);
};

#endif /* OverlapMatrix_hpp */
