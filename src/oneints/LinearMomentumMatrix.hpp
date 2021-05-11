//
//                           VELOXCHEM 1.0-RC2
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

#ifndef LinearMomentumMatrix_hpp
#define LinearMomentumMatrix_hpp

#include <array>
#include <cstdint>
#include <string>

#include "CartesianComponents.hpp"
#include "DenseMatrix.hpp"

/**
 Class CLinearMomentumMatrix stores linear momentum matrix and provides
 set of methods for handling of linear momentum matrix data.

 @author Z. Rinkevicius
 */
class CLinearMomentumMatrix
{
    /**
     The generic dense x component of linear momentum matrix (rectangular or
     square).
     */
    CDenseMatrix _xMatrix;

    /**
     The generic dense Y component of linear momentum matrix (rectangular or
     square).
     */
    CDenseMatrix _yMatrix;

    /**
     The generic dense Z component of linear momentum matrix (rectangular or
     square).
     */
    CDenseMatrix _zMatrix;

   public:
    /**
     Creates an empty linear momentum matrix object.
     */
    CLinearMomentumMatrix();

    /**
     Creates a linear momentum matrix object.

     @param xMatrix the dense matrix with X component of linear momentum integrals.
     @param yMatrix the dense matrix with Y component of linear momentum integrals.
     @param zMatrix the dense matrix with Z component of linear momentum integrals.
     */
    CLinearMomentumMatrix(const CDenseMatrix& xMatrix, const CDenseMatrix& yMatrix, const CDenseMatrix& zMatrix);

    /**
     Creates a linear momentum matrix object.

     @param matrices cartesian components of linear momentum integrals.
     */
    CLinearMomentumMatrix(const std::array<CDenseMatrix, 3>& matrices);

    /**
     Creates a linear momentum matrix object by copying other linear momentum
     matrix object.

     @param source the linear momentum matrix object.
     */
    CLinearMomentumMatrix(const CLinearMomentumMatrix& source);

    /**
     Creates a linear momentum  matrix object by moving other linear momentum
     matrix object.

     @param source the linear momentum matrix object.
     */
    CLinearMomentumMatrix(CLinearMomentumMatrix&& source) noexcept;

    /**
     Destroys a linear momentum matrix object.
     */
    ~CLinearMomentumMatrix();

    /**
     Assigns a linear momentum  matrix object by copying other linear momentum
     matrix object.

     @param source the linear momentum matrix object.
     */
    CLinearMomentumMatrix& operator=(const CLinearMomentumMatrix& source);

    /**
     Assigns a linear momentum matrix object by moving other linear momentum
     matrix object.

     @param source the linear momentum matrix object.
     */
    CLinearMomentumMatrix& operator=(CLinearMomentumMatrix&& source) noexcept;

    /**
     Compares linear momentum matrix object with other linear momentum matrix
     object.

     @param other the linear momentum matrix object.
     @return true if linear momentum matrix objects are equal, false otherwise.
     */
    bool operator==(const CLinearMomentumMatrix& other) const;

    /**
     Compares linear momentum matrix object with other linear momentum matrix
     object.

     @param other the linear momentum matrix object.
     @return true if linear momentum matrix objects are not equal, false
     otherwise.
     */
    bool operator!=(const CLinearMomentumMatrix& other) const;

    /**
     Gets string representation of linear momentum matrix.

     @return a string for printing the linear momentum matrix.
     */
    std::string getString() const;

    /**
     Gets string representation of X component of linear momentum matrix.

     @return a string for printing the X component of linear momentum matrix.
     */
    std::string getStringForComponentX() const;

    /**
     Gets string representation of Y component of linear momentum matrix.

     @return a string for printing the Y component of linear momentum matrix.
     */
    std::string getStringForComponentY() const;

    /**
     Gets string representation of Z component of linear momentum matrix.

     @return a string for printing the Z component of linear momentum matrix.
     */
    std::string getStringForComponentZ() const;

    /**
     Gets number of rows in linear momentum matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in linear momentum matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in linear momentum matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of X component of linear momentum
     matrix.

     @return the constant pointer to first element of X component of electric
     dipole matrix.
     */
    const double* xvalues() const;

    /**
     Gets constant pointer to first element of Y component of linear momentum
     matrix.

     @return the constant pointer to first element of Y component of electric
     dipole matrix.
     */
    const double* yvalues() const;

    /**
     Gets constant pointer to first element of Z component of linear momentum
     matrix.

     @return the constant pointer to first element of Z component of electric
     dipole matrix.
     */
    const double* zvalues() const;

    /**
     Gets constant pointer to first element of requested component of linear momentum
     matrix.

     @param cart requested Cartesian component of the linear momentum integrals matrix

     @return the constant pointer to first element of requested component of linear momentum
     matrix.
     */
    const double* values(cartesians cart) const;

    /**
     Gets constant pointer to first element of requested component of linear momentum
     matrix.

     @param cart requested Cartesian component of the linear momentum integrals matrix

     @return the constant pointer to first element of requested component of linear momentum
     matrix.
     */
    const double* values(int32_t cart) const;

    /**
     Converts linear momentum matrix object to text output.
     */
    std::string repr() const;
};

/**
 Converts linear momentum matrix object to text output and insert it into
 output text stream.

 @param output the output text stream.
 @param source the linear momentum matrix object.
 */
std::ostream& operator<<(std::ostream& output, const CLinearMomentumMatrix& source);

#endif /* LinearMomentumMatrix_hpp */
