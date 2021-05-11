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

#ifndef KineticEnergyMatrix_hpp
#define KineticEnergyMatrix_hpp

#include <cstdint>
#include <string>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"

/**
 Class CKineticEnergyMatrix stores general kinetic energy matrix and provides
 set of methods for handling of kinetic energy matrix data.

 @author Z. Rinkevicius
 */
class CKineticEnergyMatrix
{
    /**
     The generic dense kinetic energy matrix (rectangular or square).
     */
    CDenseMatrix _matrix;

   public:
    /**
     Creates an empty kinetic energy matrix object.
     */
    CKineticEnergyMatrix();

    /**
     Creates a kinetic energy matrix object.

     @param matrix the dense matrix with kinetic energy integrals.
     */
    CKineticEnergyMatrix(const CDenseMatrix& matrix);

    /**
     Creates a kinetic energy matrix object by copying other kinetic energy
     matrix object.

     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix(const CKineticEnergyMatrix& source);

    /**
     Creates a kinetic energy matrix object by moving other kinetic energy
     matrix object.

     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix(CKineticEnergyMatrix&& source) noexcept;

    /**
     Destroys a kinetic energy matrix object.
     */
    ~CKineticEnergyMatrix();

    /**
     Assigns a kinetic energy matrix object by copying other kinetic energy
     matrix object.

     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix& operator=(const CKineticEnergyMatrix& source);

    /**
     Assigns a kinetic energy matrix object by moving other kinetic energy
     matrix object.

     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix& operator=(CKineticEnergyMatrix&& source) noexcept;

    /**
     Compares kinetic energy matrix object with other kinetic energy matrix
     object.

     @param other the kinetic energy matrix object.
     @return true if kinetic energy matrix objects are equal, false otherwise.
     */
    bool operator==(const CKineticEnergyMatrix& other) const;

    /**
     Compares kinetic energy matrix object with other kinetic energy matrix
     object.

     @param other the kinetic energy matrix object.
     @return true if kinetic energy matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CKineticEnergyMatrix& other) const;

    /**
     Gets string representation of kinetic energy matrix.

     @return a string for printing the kinetic energy matrix.
     */
    std::string getString() const;

    /**
     Gets number of rows in kinetic energy matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in kinetic energy matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in kinetic energy matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of kinetic energy matrix.

     @return the constant pointer to first element of kinetic energy matrix.
     */
    const double* values() const;

    /**
     Computes kinetic energy for specific AO density matrix.

     @param aoDensityMatrix the AO density matrix object.
     @param iDensityMatrix the index of AO density matrix in AO density matrix object.
     @return the kinetic energy.
     */
    double getKineticEnergy(const CAODensityMatrix& aoDensityMatrix, const int32_t iDensityMatrix) const;

    /**
     Converts kinetic energy matrix object to text output and insert it into
     output text stream.

     @param output the output text stream.
     @param source the kinetic energy matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CKineticEnergyMatrix& source);
};

#endif /* KineticEnergyMatrix_hpp */
