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

#ifndef MolecularCorePotential_hpp
#define MolecularCorePotential_hpp

#include <vector>

#include "AtomCorePotential.hpp"

/// @brief Class CMolecularCorePotential stores data about molecular core potential and provides
/// set of methods for handling of molecular core potential data.
class CMolecularCorePotential
{
   public:
    /// @brief The default constructor.
    CMolecularCorePotential();

    /// @brief The constructor with vector of unique atom core potentials, vector of
    /// atom core potential indices, and atom indices.
    /// @param core_potentials The vector of unique atom core potentials.
    /// @param indices The vector of atom core potential indices.
    /// @param atom_indices The vector of atom indices.
    CMolecularCorePotential(const std::vector<CAtomCorePotential> &core_potentials, const std::vector<int> &indices, const std::vector<int> &atom_indices);

    /// @brief The default copy constructor.
    /// @param other The molecular core potential to be copied.
    CMolecularCorePotential(const CMolecularCorePotential &other);

    /// @brief The default move constructor.
    /// @param other The molecular core potential to be moved.
    CMolecularCorePotential(CMolecularCorePotential &&other) noexcept;

    /// @brief The default destructor.
    ~CMolecularCorePotential() = default;

    /// @brief The default copy assignment operator.
    /// @param other The molecular core potential to be copy assigned.
    /// @return The assigned molecular core potential.
    auto operator=(const CMolecularCorePotential &other) -> CMolecularCorePotential &;

    /// @brief The default move assignment operator.
    /// @param other The molecular core potential to be move assigned.
    /// @return The assigned molecular core potential.
    auto operator=(CMolecularCorePotential &&other) noexcept -> CMolecularCorePotential &;

    /// @brief The equality operator.
    /// @param other The molecular core potential to be compared.
    /// @return True if molecular core potentials are equal, False otherwise.
    auto operator==(const CMolecularCorePotential &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The molecular core potential to be compared.
    /// @return True if molecular core potentials are not equal, False otherwise.
    auto operator!=(const CMolecularCorePotential &other) const -> bool;

    /// @brief Adds atom core potential to molecular core potential.
    /// @param core_potential The atom core potential to be added.
    /// @param iatom The index of atom associated with core potential.
    auto add(const CAtomCorePotential &core_potential, const int iatom) -> void;
    
    /// @brief Gets vector of core potentials.
    /// @return The vector of core potentials.
    auto core_potentials() const -> std::vector<CAtomCorePotential>;
    
    /// @brief Gets vector of core potential indices.
    /// @return The vector of core potential indices.
    auto indices() const -> std::vector<int>;
    
    /// @brief Gets vector of atomic indices.
    /// @return The vector of atomic indices.
    auto atomic_indices() const -> std::vector<int>;
    
    /// @brief Gets specific core potential.
    /// @param index The index of core potential.
    /// @return The specific core potential.
    auto core_potential(const int index) const -> CAtomCorePotential; 

   private:
    /// @brief The vector of atom core potentials.
    std::vector<CAtomCorePotential> _core_potentials;

    /// @brief The vector of atom core potential indices.
    std::vector<int> _indices;

    /// @brief The vector of atom indices.
    std::vector<int> _atom_indices;
};

#endif /* MolecularCorePotential_hpp */
