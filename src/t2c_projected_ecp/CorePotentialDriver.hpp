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

#ifndef CorePotentialDriver_hpp
#define CorePotentialDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "AtomCorePotential.hpp"

/// @brief Class CCorePotentialDriver provides methods for computing two-center projected ECP integrals.
class CCorePotentialDriver
{
   public:
    /// @brief Creates a ECP integrals driver.
    CCorePotentialDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The ECP integrals driver to be copied.
    CCorePotentialDriver(const CCorePotentialDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The ECP integrals driver  to be moved.
    CCorePotentialDriver(CCorePotentialDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CCorePotentialDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The ECP integrals driver to be copy assigned.
    /// @return The assigned ECP integrals driver.
    auto operator=(const CCorePotentialDriver &other) -> CCorePotentialDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The ECP integrals driver to be move assigned.
    /// @return The assigned ECP integrals driver .
    auto operator=(CCorePotentialDriver &&other) noexcept -> CCorePotentialDriver & = delete;

    /// @brief The equality operator.
    /// @param other The ECP integrals driver  to be compared.
    /// @return True if ECP integrals drivers  are equal, False otherwise.
    auto operator==(const CCorePotentialDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The ECP integrals driver to be compared.
    /// @return True if ECP integrals drivers  are not equal, False otherwise.
    auto operator!=(const CCorePotentialDriver &other) const -> bool = delete;

    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential) const -> CMatrix;
};

#endif /* CorePotentialDriver_hpp */
