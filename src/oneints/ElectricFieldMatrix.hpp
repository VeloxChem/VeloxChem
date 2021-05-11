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

#ifndef ElectricFieldMatrix_hpp
#define ElectricFieldMatrix_hpp

#include <array>
#include <cstdint>
#include <string>

#include "CartesianComponents.hpp"
#include "DenseMatrix.hpp"

/**
 Class CElectricFieldMatrix stores electric field matrix and provides
 set of methods for handling of electric field matrix data.

 @author Z. Rinkevicius
 */
class CElectricFieldMatrix
{
    /**
     The generic dense x component of electric field matrix (rectangular or
     square).
     */
    CDenseMatrix _xMatrix;

    /**
     The generic dense Y component of electric field matrix (rectangular or
     square).
     */
    CDenseMatrix _yMatrix;

    /**
     The generic dense Z component of electric field matrix (rectangular or
     square).
     */
    CDenseMatrix _zMatrix;

   public:
    /**
     Creates an empty electric field  matrix object.
     */
    CElectricFieldMatrix();

    /**
     Creates a electric field matrix object.

     @param xMatrix the dense matrix with X component of electric field
     integrals.
     @param yMatrix the dense matrix with Y component of electric field
     integrals.
     @param zMatrix the dense matrix with Z component of electric field
     integrals.
     */
    CElectricFieldMatrix(const CDenseMatrix& xMatrix, const CDenseMatrix& yMatrix, const CDenseMatrix& zMatrix);

    /**
     Creates a electric field matrix object.

     @param matrices array of dense matrices with components of electric field
     integrals.
     */
    CElectricFieldMatrix(const std::array<CDenseMatrix, 3>& matrix);

    /**
     Creates a electric field matrix object by copying other electric field
     matrix object.

     @param source the electric field matrix object.
     */
    CElectricFieldMatrix(const CElectricFieldMatrix& source);

    /**
     Creates a electric field matrix object by moving other electric field
     matrix object.

     @param source the electric field matrix object.
     */
    CElectricFieldMatrix(CElectricFieldMatrix&& source) noexcept;

    /**
     Destroys a electric field matrix object.
     */
    ~CElectricFieldMatrix();

    /**
     Assigns a electric field matrix object by copying other electric field
     matrix object.

     @param source the electric field matrix object.
     */
    CElectricFieldMatrix& operator=(const CElectricFieldMatrix& source);

    /**
     Assigns a electric field matrix object by moving other electric field
     matrix object.

     @param source the electric field matrix object.
     */
    CElectricFieldMatrix& operator=(CElectricFieldMatrix&& source) noexcept;

    /**
     Compares electric field matrix object with other electric field matrix
     object.

     @param other the electric field matrix object.
     @return true if electric field matrix objects are equal, false otherwise.
     */
    bool operator==(const CElectricFieldMatrix& other) const;

    /**
     Compares electric field matrix object with other electric field matrix
     object.

     @param other the electric field matrix object.
     @return true if electric field matrix objects are not equal, false
     otherwise.
     */
    bool operator!=(const CElectricFieldMatrix& other) const;

    /**
     Gets string representation of X component of electric field matrix.

     @return a string for printing the X component of electric field matrix.
     */
    std::string getStringForComponentX() const;

    /**
     Gets string representation of Y component of electric field matrix.

     @return a string for printing the Y component of electric field matrix.
     */
    std::string getStringForComponentY() const;

    /**
     Gets string representation of Z component of electric field matrix.

     @return a string for printing the Z component of electric field matrix.
     */
    std::string getStringForComponentZ() const;

    /**
     Gets string representation of electric field matrix.

     @return a string for printing the components of electric field matrix.
     */
    std::string getString() const;

    /**
     Gets number of rows in electric field matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in electric field matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in electric field matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of X component of electric field
     matrix.

     @return the constant pointer to first element of X component of electric
     field matrix.
     */
    const double* xvalues() const;

    /**
     Gets constant pointer to first element of Y component of electric field
     matrix.

     @return the constant pointer to first element of Y component of electric
     field matrix.
     */
    const double* yvalues() const;

    /**
     Gets constant pointer to first element of Z component of electric field
     matrix.

     @return the constant pointer to first element of Z component of electric
     field matrix.
     */
    const double* zvalues() const;

    /**
     Gets constant pointer to first element of requested component of electric field
     matrix.

     @param cart requested Cartesian component of the electric field integrals matrix

     @return the constant pointer to first element of requested component of electric field
     matrix.
     */
    const double* values(cartesians cart) const;

    /**
     Gets constant pointer to first element of requested component of electric field
     matrix.

     @param cart requested Cartesian component of the electric field integrals matrix

     @return the constant pointer to first element of requested component of electric field
     matrix.
     */
    const double* values(int32_t cart) const;

    /**
     Converts electric field matrix object to text output.
     */
    std::string repr() const;
};

/**
 Converts electric field matrix object to text output and insert it into
 output text stream.

 @param output the output text stream.
 @param source the electric field matrix object.
 */
std::ostream& operator<<(std::ostream& output, const CElectricFieldMatrix& source);

#endif /* ElectricFieldMatrix_hpp */
