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
     @param denType the type (restricted or unrestricted) of density matrices.
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
    int getNumberOfDensityMatrices() const;

    /**
     Gets number of rows in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of rows.
     */
    int getNumberOfRows(const int iDensityMatrix) const;

    /**
     Gets number of columns in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of columns.
     */
    int getNumberOfColumns(const int iDensityMatrix) const;

    /**
     Gets number of elements in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of elements.
     */
    int getNumberOfElements(const int iDensityMatrix) const;

    /**
     Gets constant pointer to first element of spin-alpha density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of spin-alpha density matrix.
     */
    const double* alphaDensity(const int iDensityMatrix) const;

    /**
     Gets constant pointer to first element of spin-beta density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of spin-beta density matrix.
     */
    const double* betaDensity(const int iDensityMatrix) const;

    /**
     Gets constant reference to density matrix as dense matrix object.

     @param iDensityMatrix the index of density matrix.
     @return the constant reference to density matrix.
     */
    const CDenseMatrix& getReferenceToDensity(const int iDensityMatrix) const;
};

#endif /* AODensityMatrix_hpp */
