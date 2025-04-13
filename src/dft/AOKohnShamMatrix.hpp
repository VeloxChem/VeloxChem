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

#ifndef AOKohnShamMatrix_hpp
#define AOKohnShamMatrix_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"

/**
 Class CAOKohnShamMatrix stores set of AO Kohn-Sham matrices and provides set of methods
 for handling of AO Kohn-Sham matrices data.
 */
class CAOKohnShamMatrix
{
    /**
     The set of exchange-correlation matrices.
     */
    std::vector<CDenseMatrix> _xcMatrices;

    /**
     The flag indicating form of Kohn-Sham matrix: restricted or unrestricted.
     */
    bool _xcRestricted;

    /**
     The number of electrons obtained during integration of exchange-correlation matrices.
     */
    double _xcElectrons;

    /**
     The exchange-correlation energy associated with Kohn-Sham matrix.
     */
    double _xcEnergy;

   public:
    /**
     Creates an empty AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix();

    /**
    Creates a AO Kohn-Sham matrix object.

    @param nRows the number of rows in exchange-correlation matrices.
    @param nColumns the numner of columns in exchange-correlation matrices.
    @param xcRestricted the flag indicating restricted on unrestricted form of Kohn-Sham matrix.
    */
    CAOKohnShamMatrix(const int64_t nRows, const int64_t nColumns, const bool xcRestricted);

    /**
     Creates a AO Kohn-Sham matrix object by copying other AO Kohn-Sham matrix object.

     @param source the AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix(const CAOKohnShamMatrix& source);

    /**
     Creates a AO Kohn-Sham matrix object by moving other AO Kohn-Sham matrix object.

     @param source the AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix(CAOKohnShamMatrix&& source) noexcept;

    /**
     Destroys a AO Kohn-Sham matrix object.
     */
    ~CAOKohnShamMatrix();

    /**
     Assigns a AO Kohn-Sham matrix object by copying other AO Kohn-Sham matrix object.

     @param source the AO Kohn-Sham matrix object.
     */
    auto operator=(const CAOKohnShamMatrix& source) -> CAOKohnShamMatrix&;

    /**
     Assigns a AO Kohn-Sham matrix object by moving other AO Kohn-Sham matrix object.

     @param source the AO Kohn-Sham matrix object.
     */
    auto operator=(CAOKohnShamMatrix&& source) noexcept -> CAOKohnShamMatrix&;

    /**
     Compares AO Kohn-Sham matrix object with other AO Kohn-Sham matrix object.

     @param other the AO Kohn-Sham matrix object.
     @return true if AO Kohn-Sham matrix objects are equal, false otherwise.
     */
    auto operator==(const CAOKohnShamMatrix& other) const -> bool;

    /**
     Compares AO Kohn-Sham matrix object with other AO Kohn-Sham matrix object.

     @param other the AO Kohn-Sham matrix object.
     @return true if AO Kohn-Sham matrix objects are not equal, false otherwise.
     */
    auto operator!=(const CAOKohnShamMatrix& other) const -> bool;

    /**
     Resets all elements of AO Kohn-Sham matrix to zero.
     */
    auto zero() -> void;

    /**
     Symmetrizes and scales AO Kohn-Sham matrix.
     */
    auto symmetrizeAndScale(const double factor) -> void;

    /**
     Sets number of electron obtained by integrating Kohn-Sham matrix.

     @param xcElectrons the number of electrons.
     */
    auto setNumberOfElectrons(const double xcElectrons) -> void;

    /**
     Set exchange-correlation energy associated with Kohn-Sham matrix.

     @param xcEnergy the exchange-correlation energy.
     */
    auto setExchangeCorrelationEnergy(const double xcEnergy) -> void;

    /**
     Checks if the Fock matrices are restricted.

     @return true if the Fock matrices are restricted.
     */
    auto isRestricted() const -> bool;

    /**
     Gets number of electrons obtained by integrating Kohn-Sham matrix.

     @return the number of electrons.
     */
    auto getNumberOfElectrons() const -> double;

    /**
     Gets exchange-correlation energy associated with Kohn-Sham matrix.

     @return the exchange-correlation energy.
     */
    auto getExchangeCorrelationEnergy() const -> double;

    /**
     Gets number of rows in Kohn-Sham matrix.

     @return the number of rows.
     */
    auto getNumberOfRows() const -> int64_t;

    /**
     Gets number of columns in Kohn-Sham matrix.

     @return the number of columns.
     */
    auto getNumberOfColumns() const -> int64_t;

    /**
     Gets number of elements in Kohn-Sham matrix.

     @return the number of elements.
     */
    auto getNumberOfElements() const -> int64_t;

    /**
     Gets constant reference to alpha-spin Kohn-Sham matrix.

     @return the constant reference to alpha-spin Kohn-Sham matrix.
     */
    auto getReferenceToAlphaKohnSham() const -> const CDenseMatrix&;

    /**
     Gets constant reference to beta-spin Kohn-Sham matrix.

     @return the constant reference to beta-spin Kohn-Sham matrix.
     */
    auto getReferenceToBetaKohnSham() const -> const CDenseMatrix&;

    /**
     Gets const pointer to alpha-spin Kohn-Sham matrix data.

     @return the const pointer to alpha-spin Kohn-Sham matrix data.
     */
    auto getPointerToAlphaValues() const -> const double*;

    /**
     Gets pointer to alpha-spin Kohn-Sham matrix data.

     @return the pointer to alpha-spin Kohn-Sham matrix data.
     */
    auto getPointerToAlphaValues() -> double*;

    /**
     Gets const pointer to beta-spin Kohn-Sham matrix data.

     @return the const pointer to beta-spin Kohn-Sham matrix data.
     */
    auto getPointerToBetaValues() const -> const double*;

    /**
     Gets pointer to beta-spin Kohn-Sham matrix data.

     @return the pointer to beta-spin Kohn-Sham matrix data.
     */
    auto getPointerToBetaValues() -> double*;
};

#endif /* AOKohnShamMatrix_hpp */
