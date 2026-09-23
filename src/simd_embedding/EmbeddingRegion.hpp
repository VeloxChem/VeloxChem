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


#ifndef EmbeddingRegion_hpp
#define EmbeddingRegion_hpp

#include <string>
#include <vector>

#include "Molecule.hpp"
#include "PolarizableForceField.hpp"

/// @brief Class CEmbeddingRegion holds one region of a polarizable environment:
/// the molecules which make it up and the force fields they are described by.
///
/// @note The molecules are instances and the force fields are kinds. Each
/// molecule carries an identifier which is the index of its force field, so a box
/// of five hundred waters is five hundred geometries, five hundred identifiers
/// which are all the same, and one force field. The force fields are unique: one
/// name stands for one set of parameters, and a second force field of a name
/// already held is either the same one, and reuses it, or is different, and is
/// refused.
///
/// @note A molecule and its identifier are only ever added together, so the two
/// cannot fall out of step. A force field is checked against the molecule it is
/// given to: it carries one site for each atom, and it is the order of the atoms
/// which says which site belongs to which. A count which does not match is the
/// only sign of a force field applied to the wrong kind of molecule that can be
/// had from the parameters alone, so it is taken here.
class CEmbeddingRegion
{
   public:
    /// @brief The default constructor, which makes a region whose molecules may
    /// be polarizable.
    CEmbeddingRegion() = default;

    /// @brief The constructor.
    /// @param allows_polarizabilities Whether a force field of this region may
    /// carry a polarizability.
    CEmbeddingRegion(const bool allows_polarizabilities);

    /// @brief Checks whether a force field of this region may carry a
    /// polarizability.
    auto allows_polarizabilities() const -> bool;

    /// @brief Adds a force field, or finds the one already held.
    /// @param force_field The force field.
    /// @return The index of the force field, which is the identifier a molecule
    /// of this kind is added with.
    /// @note A force field of a name already held is not added again. If it says
    /// the same thing as the one held, the index of that one is answered; if it
    /// says something else, it is refused, because the two could not then be told
    /// apart by name.
    auto add_force_field(const CPolarizableForceField &force_field) -> int;

    /// @brief Adds a molecule of a force field already held.
    /// @param molecule The molecule.
    /// @param identifier The index of its force field.
    auto add_molecule(const CMolecule &molecule, const int identifier) -> void;

    /// @brief Adds a molecule and the force field which describes it.
    /// @param molecule The molecule.
    /// @param force_field The force field, which is added if it is a new one and
    /// found if it is not.
    auto add_molecule(const CMolecule &molecule, const CPolarizableForceField &force_field) -> void;

    /// @brief Gets the index of the force field of this name.
    /// @param name The name.
    /// @return The index, or minus one if the region holds no force field of that
    /// name.
    auto index_of_force_field(const std::string &name) const -> int;

    /// @brief Gets the number of molecules.
    auto number_of_molecules() const -> int;

    /// @brief Gets the number of force fields, which is the number of kinds of
    /// molecule the region holds.
    auto number_of_force_fields() const -> int;

    /// @brief Gets a molecule.
    /// @param index The index of the molecule.
    auto get_molecule(const int index) const -> const CMolecule &;

    /// @brief Gets the identifier of a molecule, which is the index of its force
    /// field.
    /// @param index The index of the molecule.
    auto get_identifier(const int index) const -> int;

    /// @brief Gets a force field.
    /// @param index The index of the force field.
    auto get_force_field(const int index) const -> const CPolarizableForceField &;

    /// @brief Gets the force field of a molecule.
    /// @param index The index of the molecule.
    auto force_field_of(const int index) const -> const CPolarizableForceField &;

    /// @brief Gets the molecules.
    auto get_molecules() const -> const std::vector<CMolecule> &;

    /// @brief Gets the identifiers, one for each molecule.
    auto get_identifiers() const -> const std::vector<int> &;

    /// @brief Gets the force fields.
    auto get_force_fields() const -> const std::vector<CPolarizableForceField> &;

    /// @brief Gets the number of atoms of all the molecules, which is the number
    /// of sites the region carries.
    auto number_of_sites() const -> int;

    /// @brief Gets the number of sites which carry a polarizability.
    /// @return The number of them.
    /// @note This is what the induced dipoles are solved for, and it is not the
    /// number of sites: a force field may leave some of its atoms unpolarizable.
    auto number_of_polarizable_sites() const -> int;

    /// @brief Checks whether any molecule of the region is polarizable.
    auto is_polarizable() const -> bool;

   private:
    /// @brief Whether a force field of this region may carry a polarizability.
    bool _allows_polarizabilities = true;

    /// @brief The molecules of the region.
    std::vector<CMolecule> _molecules;

    /// @brief The identifiers, one for each molecule, which are the indices of
    /// their force fields.
    std::vector<int> _identifiers;

    /// @brief The force fields, one for each kind of molecule the region holds.
    std::vector<CPolarizableForceField> _force_fields;
};

#endif /* EmbeddingRegion_hpp */
