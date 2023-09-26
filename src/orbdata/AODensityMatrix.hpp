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

#ifndef AODensityMatrix_hpp
#define AODensityMatrix_hpp

#include <cstdint>
#include <vector>

#include "DenseMatrix.hpp"

/**
 Enumerate class denmat:

 Defines supported density matrix types:
 denmat::rest   - the restricted density matrix
 denmat::unrest - the unrestricted density matrix
 */
enum class denmat
{
    rest,
    unrest
};

/**
 Class CAODensityMatrix stores set of AO density matrices and provides
 set of methods for handling of AO density matrices data.

 @author Z. Rinkevicius
 */
class CAODensityMatrix
{
    /**
     The set of density matrices.
     */
    std::vector<CDenseMatrix> _denMatrices;

    /**
     The type of density matrices.
     */
    denmat _denType;

   public:
    /**
     Creates an empty AO density matrix object.
     */
    CAODensityMatrix();

    /**
     Creates a AO density matrix object.

     @param denMatrices the set of density matrices.
     @param denType the type (restricted, unrestricted, etc) of density matrices.
     */
    CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices, const denmat denType);

    /**
     Creates a AO density matrix object by copying other AO density matrix
     object.

     @param source the AO density matrix object.
     */
    CAODensityMatrix(const CAODensityMatrix& source);

    /**
     Creates a AO density matrix object by moving other AO density matrix
     object.

     @param source the AO density matrix object.
     */
    CAODensityMatrix(CAODensityMatrix&& source) noexcept;

    /**
     Destroys a AO density matrix object.
     */
    ~CAODensityMatrix();

    /**
     Assigns a AO density matrix object by copying other AO density matrix
     object.

     @param source the AO density matrix object.
     */
    CAODensityMatrix& operator=(const CAODensityMatrix& source);

    /**
     Assigns a AO density matrix object by moving other AO density matrix
     object.

     @param source the AO density matrix object.
     */
    CAODensityMatrix& operator=(CAODensityMatrix&& source) noexcept;

    /**
     Compares AO density matrix object with other AO density matrix object.

     @param other the AO density matrix object.
     @return true if AO density matrix objects are equal, false otherwise.
     */
    bool operator==(const CAODensityMatrix& other) const;

    /**
     Compares AO density matrix object with other AO density matrix object.

     @param other the AO density matrix object.
     @return true if AO density matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CAODensityMatrix& other) const;

    /**
     Gets type of density matrix.

     @return the type of density matrix.
     */
    denmat getDensityType() const;

    /**
     Checks if AO density matrix is of closed-shell type.

     @return true if AO density matrix is of closed-shell type.
     */
    bool isClosedShell() const;

    /**
     Gets number of density matrices.

     @return the number of density matrices.
     */
    int64_t getNumberOfDensityMatrices() const;

    /**
     Gets number of rows in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of rows.
     */
    int64_t getNumberOfRows(const int64_t iDensityMatrix) const;

    /**
     Gets number of columns in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of columns.
     */
    int64_t getNumberOfColumns(const int64_t iDensityMatrix) const;

    /**
     Gets number of elements in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of elements.
     */
    int64_t getNumberOfElements(const int64_t iDensityMatrix) const;

    /**
     Gets constant pointer to first element of spin-alpha density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of spin-alpha density matrix.
     */
    const double* alphaDensity(const int64_t iDensityMatrix) const;

    /**
     Gets constant pointer to first element of spin-beta density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of spin-beta density matrix.
     */
    const double* betaDensity(const int64_t iDensityMatrix) const;

    /**
     Gets constant reference to spin-alpha density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant reference to density matrix.
     */
    const CDenseMatrix& getReferenceToAlphaDensity(const int64_t iDensityMatrix) const;

    /**
     Gets constant reference to spin-beta density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant reference to density matrix.
     */
    const CDenseMatrix& getReferenceToBetaDensity(const int64_t iDensityMatrix) const;
};

#endif /* AODensityMatrix_hpp */
