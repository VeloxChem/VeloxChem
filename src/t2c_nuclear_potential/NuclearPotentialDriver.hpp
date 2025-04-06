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

#ifndef NuclearPotentialDriver_hpp
#define NuclearPotentialDriver_hpp

#include <vector>

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CNuclearPotentialDriver provides methods for computing two-center nuclear potential integrals.
class CNuclearPotentialDriver
{
   public:
    /// @brief Creates an nuclear potential integrals driver.
    CNuclearPotentialDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The nuclear potential integrals driver to be copied.
    CNuclearPotentialDriver(const CNuclearPotentialDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The nuclear potential integrals driver  to be moved.
    CNuclearPotentialDriver(CNuclearPotentialDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CNuclearPotentialDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The nuclear potential integrals driver to be copy assigned.
    /// @return The assigned nuclear potential integrals driver.
    auto operator=(const CNuclearPotentialDriver &other) -> CNuclearPotentialDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The nuclear potential integrals driver to be move assigned.
    /// @return The assigned nuclear potential integrals driver .
    auto operator=(CNuclearPotentialDriver &&other) noexcept -> CNuclearPotentialDriver & = delete;

    /// @brief The equality operator.
    /// @param other The nuclear potential integrals driver  to be compared.
    /// @return True if nuclear potential integrals drivers  are equal, False otherwise.
    auto operator==(const CNuclearPotentialDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The nuclear potential integrals driver to be compared.
    /// @return True if nuclear potential integrals drivers  are not equal, False otherwise.
    auto operator!=(const CNuclearPotentialDriver &other) const -> bool = delete;

    /// @brief Computes nuclear potential matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule) const -> CMatrix;

    /// @brief Computes nuclear potential matrix for given set of external charges,  molecule and molecular basis.
    /// @param charges The vector of external charges.
    /// @param coordinates The vector of coordinates of external charges.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The nuclear potential matrix.
    auto compute(const std::vector<double>         &charges,
                 const std::vector<TPoint<double>> &coordinates,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule) const -> CMatrix;
};

#endif /* NuclearPotentialDriver_hpp */
