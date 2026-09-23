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


#ifndef PolarizableEmbedding_hpp
#define PolarizableEmbedding_hpp

#include <array>
#include <utility>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"

#include "EmbeddingRegion.hpp"

/// @brief Class CPolarizableEmbedding holds the environment a solute is embedded
/// in, as two regions: one whose molecules polarize and one whose molecules do
/// not.
///
/// @note The solute is not here. This is what surrounds it, and it is held apart
/// from the molecule the wave function is of.
///
/// @note The two regions differ in one thing: the nonpolarizable one refuses a
/// force field which carries a polarizability, so a force field put in the wrong
/// region is refused rather than quietly stripped of the part of it which would
/// have mattered. A species which is to sit in the nonpolarizable region is given
/// a force field without polarizabilities, and that is a statement about the
/// calculation rather than something to be inferred from where it was put.
class CPolarizableEmbedding
{
   public:
    /// @brief The default constructor, which makes an empty environment.
    CPolarizableEmbedding();

    /// @brief Gets the polarizable region.
    auto get_polarizable_region() const -> const CEmbeddingRegion &;

    /// @brief Gets the polarizable region, to add to.
    auto polarizable_region() -> CEmbeddingRegion &;

    /// @brief Gets the nonpolarizable region.
    auto get_nonpolarizable_region() const -> const CEmbeddingRegion &;

    /// @brief Gets the nonpolarizable region, to add to.
    auto nonpolarizable_region() -> CEmbeddingRegion &;

    /// @brief Gets the number of molecules of both regions.
    auto number_of_molecules() const -> int;

    /// @brief Gets the number of sites of both regions.
    auto number_of_sites() const -> int;

    /// @brief Gets the number of sites which carry a polarizability.
    /// @return The number of them, which only the polarizable region contributes
    /// to.
    auto number_of_polarizable_sites() const -> int;

    /// @brief Checks whether the environment polarizes at all.
    /// @return True if it does.
    /// @note An environment which does not is an ordinary one of fixed multipoles,
    /// and the induced dipoles need not be solved for.
    auto is_polarizable() const -> bool;

    /// @brief Gets the permanent multipoles of an order of the whole
    /// environment, put onto the molecules which carry them.
    /// @param order The order: zero for the charges, one for the dipoles, two
    /// for the quadrupoles.
    /// @return The multipoles of both regions, the polarizable one first.
    /// @note Both regions, because the permanent multipoles of the whole
    /// environment act on the wave function. What tells the two regions apart
    /// is the induction which is added on top of these, not these.
    auto permanent_multipoles(const int order) const -> const TPermanentMultipoles &;

    /// @brief Computes what the permanent charges of the environment do to the
    /// nuclei of the quantum region.
    /// @param molecule The molecule of the quantum region.
    /// @param basis The molecular basis of it, which says what is left of each
    /// nucleus once a core potential describes the rest.
    /// @return The energy.
    /// @note The other half of the permanent electrostatics. The Fock matrix
    /// carries what the environment does to the electrons and this is what it
    /// does to the nuclei; a total energy without it is wrong by the whole of
    /// this term.
    /// @note The effective charges and not the bare ones, which for an atom
    /// whose core is a potential are not the same number.
    auto permanent_nuclear_energy(const CMolecule &molecule, const CMolecularBasis &basis) const -> double;

   private:
    /// @brief The multipoles of each order of both regions together.
    mutable std::array<TPermanentMultipoles, 3> _multipoles;

    /// @brief The versions of the two regions each order was last built from.
    mutable std::array<std::pair<size_t, size_t>, 3> _multipoles_version;

    /// @brief Whether each order's multipoles have ever been built.
    mutable std::array<bool, 3> _has_multipoles = {false, false, false};

    /// @brief The region whose molecules polarize.
    CEmbeddingRegion _polarizable;

    /// @brief The region whose molecules do not polarize.
    CEmbeddingRegion _nonpolarizable;
};

#endif /* PolarizableEmbedding_hpp */
