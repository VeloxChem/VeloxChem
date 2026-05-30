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

#ifndef newints_MolecularBasisOutline_hpp
#define newints_MolecularBasisOutline_hpp

#include <cstddef>
#include <vector>

class CMolecularBasis;

namespace newints {

/// @brief Describes the contracted-GTO indexing of a molecular basis.
///
/// Each contracted basis function (shell) is assigned a global index equal to
/// the starting offset of its angular-expanded block: the index advances by
/// (2 l + 1) per shell, in atom-major order (within an atom, the basis's
/// natural angular-momentum order). These indices are the keys used by the
/// sparse matrices of the newints subsystem; a block at index i with angular
/// momentum l_a spans atomic-orbital rows [i, i + 2 l_a + 1).
///
/// For water / STO-3G (O:[2s,1p], H:[1s], H:[1s]) the shell indices are
/// 0, 1, 2, 5, 6 and the total angular-expanded dimension is 7.
class MolecularBasisOutline
{
   public:
    /// @brief The default constructor, creating an empty outline.
    MolecularBasisOutline() = default;

    /// @brief The constructor building the indexing map from a molecular basis.
    /// @param basis The molecular basis.
    explicit MolecularBasisOutline(const CMolecularBasis &basis);

    /// @brief The equality operator.
    /// @param other The outline to compare with.
    /// @return True if both outlines have identical indexing maps.
    auto operator==(const MolecularBasisOutline &other) const -> bool;

    /// @brief Gets the number of atoms.
    /// @return The number of atoms.
    auto number_of_atoms() const -> std::size_t;

    /// @brief Gets the number of contracted basis functions (shells).
    /// @return The number of contracted basis functions.
    auto number_of_basis_functions() const -> std::size_t;

    /// @brief Gets the total angular-expanded dimension, i.e. the number of
    /// atomic orbitals.
    /// @return The number of atomic orbitals.
    auto number_of_atomic_orbitals() const -> std::size_t;

    /// @brief Gets the contracted-GTO indices of one atom.
    /// @param atom The atom index.
    /// @return The vector of contracted-GTO (shell) indices on the atom.
    auto basis_function_indices(const int atom) const -> std::vector<int>;

    /// @brief Gets the contracted-GTO index of each shell (atom-major order).
    /// @return The vector of shell indices.
    auto indices() const -> std::vector<int>;

    /// @brief Gets the angular momentum of each shell (atom-major order).
    /// @return The vector of angular momenta.
    auto angular_momenta() const -> std::vector<int>;

    /// @brief Gets the owning atom of each shell (atom-major order).
    /// @return The vector of atom indices.
    auto atom_indices() const -> std::vector<int>;

   private:
    /// @brief The contracted-GTO index (angular-expanded offset) of each shell.
    std::vector<int> _indices;

    /// @brief The angular momentum of each shell.
    std::vector<int> _angular_momenta;

    /// @brief The owning atom of each shell.
    std::vector<int> _atom_indices;

    /// @brief The per-atom slice boundaries into the shell arrays; atom a owns
    /// shell positions [_atom_shell_offsets[a], _atom_shell_offsets[a + 1]).
    std::vector<int> _atom_shell_offsets;
};

}  // namespace newints

#endif /* newints_MolecularBasisOutline_hpp */
