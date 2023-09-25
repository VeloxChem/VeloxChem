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

#ifndef AOKohnShamMatrix_hpp
#define AOKohnShamMatrix_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"

/**
 Class CAOKohnShamMatrix stores set of AO Kohn-Sham matrices and provides set of methods
 for handling of AO Kohn-Sham matrices data.

 @author Z. Rinkevicius
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
     Gets pointer to alpha-spin Kohn-Sham matrix data.

     @return the pointer to alpha-spin Kohn-Sham matrix data.
     */
    auto getPointerToAlphaValues() -> double*;

    /**
     Gets pointer to beta-spin Kohn-Sham matrix data.

     @return the pointer to beta-spin Kohn-Sham matrix data.
     */
    auto getPointerToBetaValues() -> double*;
};

#endif /* AOKohnShamMatrix_hpp */
