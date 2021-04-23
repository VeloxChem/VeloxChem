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

#ifndef ElectricFieldGradientMatrix_hpp
#define ElectricFieldGradientMatrix_hpp

#include <cstdint>
#include <string>

#include "DenseMatrix.hpp"

/**
 Class CElectricFieldMatrix stores electric field gradient matrix and provides
 set of methods for handling of electric field gradient matrix data.

 @author Z. Rinkevicius
 */
class CElectricFieldGradientMatrix
{
    /**
     The generic dense xx component of electric field gradient matrix (rectangular or
     square).
     */
    CDenseMatrix _xxMatrix;

    /**
     The generic dense xy component of electric field gradient matrix (rectangular or
     square).
     */
    CDenseMatrix _xyMatrix;

    /**
     The generic dense xz component of electric field gradient matrix (rectangular or
     square).
     */
    CDenseMatrix _xzMatrix;

    /**
     The generic dense yy component of electric field gradient matrix (rectangular or
     square).
     */
    CDenseMatrix _yyMatrix;

    /**
     The generic dense yz component of electric field gradient matrix (rectangular or
     square).
     */
    CDenseMatrix _yzMatrix;

    /**
     The generic dense zz component of electric field gradient matrix (rectangular or
     square).
     */
    CDenseMatrix _zzMatrix;

   public:
    /**
     Creates an empty electric field gradient matrix object.
     */
    CElectricFieldGradientMatrix();

    /**
     Creates a electric field geadient matrix object.

     @param xxMatrix the dense matrix with XX component of electric field
     integrals.
     @param xyMatrix the dense matrix with XY component of electric field
     integrals.
     @param xzMatrix the dense matrix with XZ component of electric field
     integrals.
     @param yyMatrix the dense matrix with YY component of electric field
     integrals.
     @param yzMatrix the dense matrix with YZ component of electric field
     integrals.
     @param zzMatrix the dense matrix with ZZ component of electric field
     integrals.
     */
    CElectricFieldGradientMatrix(const CDenseMatrix& xxMatrix,
                                 const CDenseMatrix& xyMatrix,
                                 const CDenseMatrix& xzMatrix,
                                 const CDenseMatrix& yyMatrix,
                                 const CDenseMatrix& yzMatrix,
                                 const CDenseMatrix& zzMatrix);

    /**
     Creates a electric field gradient matrix object by copying other electric field
     gradient matrix object.

     @param source the electric field gradient matrix object.
     */
    CElectricFieldGradientMatrix(const CElectricFieldGradientMatrix& source);

    /**
     Creates a electric field gradient matrix object by moving other electric field
     gradient matrix object.

     @param source the electric field gradient matrix object.
     */
    CElectricFieldGradientMatrix(CElectricFieldGradientMatrix&& source) noexcept;

    /**
     Destroys a electric field gradient matrix object.
     */
    ~CElectricFieldGradientMatrix();

    /**
     Assigns a electric field gradient matrix object by copying other electric field
     gradient matrix object.

     @param source the electric field gradient matrix object.
     */
    CElectricFieldGradientMatrix& operator=(const CElectricFieldGradientMatrix& source);

    /**
     Assigns a electric field gradient matrix object by moving other electric field
     gradient matrix object.

     @param source the electric field gradient matrix object.
     */
    CElectricFieldGradientMatrix& operator=(CElectricFieldGradientMatrix&& source) noexcept;

    /**
     Compares electric field gradient matrix object with other electric field gradient matrix
     object.

     @param other the electric field gradient matrix object.
     @return true if electric field gradient matrix objects are equal, false otherwise.
     */
    bool operator==(const CElectricFieldGradientMatrix& other) const;

    /**
     Compares electric field gradient matrix object with other electric field gradient matrix
     object.

     @param other the electric field gradient matrix object.
     @return true if electric field gradient matrix objects are not equal, false
     otherwise.
     */
    bool operator!=(const CElectricFieldGradientMatrix& other) const;

    /**
     Gets string representation of XX component of electric field gradient matrix.

     @return a string for printing the XX component of electric field gradient matrix.
     */
    std::string getStringForComponentXX() const;

    /**
     Gets string representation of XY component of electric field gradient matrix.

     @return a string for printing the XY component of electric field gradient matrix.
     */
    std::string getStringForComponentXY() const;

    /**
     Gets string representation of XZ component of electric field gradient matrix.

     @return a string for printing the XZ component of electric field gradient matrix.
     */
    std::string getStringForComponentXZ() const;

    /**
     Gets string representation of YY component of electric field gradient matrix.

     @return a string for printing the YY component of electric field gradient matrix.
     */
    std::string getStringForComponentYY() const;

    /**
     Gets string representation of YZ component of electric field gradient matrix.

     @return a string for printing the YZ component of electric field gradient matrix.
     */
    std::string getStringForComponentYZ() const;

    /**
     Gets string representation of ZZ component of electric field gradient matrix.

     @return a string for printing the ZZ component of electric field gradient matrix.
     */
    std::string getStringForComponentZZ() const;

    /**
     Gets number of rows in electric field gradient matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in electric field gradient matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in electric field gradient matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of XX component of electric field gradient
     matrix.

     @return the constant pointer to first element of XX component of electric
     field gradient matrix.
     */
    const double* xxvalues() const;

    /**
     Gets constant pointer to first element of XY component of electric field gradient
     matrix.

     @return the constant pointer to first element of XY component of electric
     field gradient matrix.
     */
    const double* xyvalues() const;

    /**
     Gets constant pointer to first element of XZ component of electric field gradient
     matrix.

     @return the constant pointer to first element of XZ component of electric
     field gradient matrix.
     */
    const double* xzvalues() const;

    /**
     Gets constant pointer to first element of YY component of electric field gradient
     matrix.

     @return the constant pointer to first element of YY component of electric
     field gradient matrix.
     */
    const double* yyvalues() const;

    /**
     Gets constant pointer to first element of YZ component of electric field gradient
     matrix.

     @return the constant pointer to first element of YZ component of electric
     field gradient matrix.
     */
    const double* yzvalues() const;

    /**
     Gets constant pointer to first element of ZZ component of electric field gradient
     matrix.

     @return the constant pointer to first element of ZZ component of electric
     field gradient matrix.
     */
    const double* zzvalues() const;

    /**
     Converts electric field gradient matrix object to text output and insert it into
     output text stream.

     @param output the output text stream.
     @param source the electric field gradient matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CElectricFieldGradientMatrix& source);
};

#endif /* ElectricFieldGradientMatrix_hpp */
