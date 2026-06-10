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

#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <vector>

#include "AtomBasis.hpp"
#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"
#include "ScreenedBasisFunctionPair.hpp"

namespace gtofunc {  // gtofunc namespace

/// @brief Creates vector of basis functions blocks.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @return The vector of basis functions blocks.
auto make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule) -> std::vector<CGtoBlock>;

/// @brief Creates vector of basis functions blocks for selected atoms.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param atoms The vector of atoms to select.
/// @return The vector of basis functions blocks.
auto make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int> &atoms) -> std::vector<CGtoBlock>;

/**
 Gets number of atomic orbitals from vector of contracted GTOs blocks.
 @param gto_blocks the vector of contracted GTOs blocks.
 @return the number of atomic orbitals.
 */
auto getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int;

/// @brief Gets Cartesian coordinates (in au) of atoms associated with a
/// specific atom basis in molecule.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param atom_basis The atom basis to select atoms for; matched by basis name
/// and chemical element identifier.
/// @return The vector of atom Cartesian coordinates.
auto atom_coordinates(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomBasis &atom_basis) -> std::vector<TPoint<double>>;

/// @brief Gets indices of atoms associated with a specific atom basis in
/// molecular basis. The atoms are listed in the same order as the coordinates
/// returned by atom_coordinates for the same atom basis.
/// @param basis The molecular basis.
/// @param atom_basis The atom basis to select atoms for; matched by basis name
/// and chemical element identifier.
/// @return The vector of atom indices.
auto atom_indices(const CMolecularBasis &basis, const CAtomBasis &atom_basis) -> std::vector<int>;

/// @brief Creates vector of screened basis function pairs for a molecular basis
/// and molecule. The unique symmetric pairs of atom bases are enumerated; for
/// each atom basis pair the atom coordinates and indices are retrieved and a
/// screened basis function pair is built for every basis function pair
/// combination. When the two atom bases are identical a triangular loop over
/// basis functions is used (with the symmetric, single-function constructor on
/// the diagonal), otherwise a full loop is used.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param keep_pair The predicate selecting important atom pairs; called as
/// keep_pair(bra_center, ket_center) and returning true to keep the pair.
/// @return The vector of screened basis function pairs.
template <typename Predicate>
auto
make_screened_basis_function_pairs(const CMolecularBasis &basis, const CMolecule &molecule, Predicate &&keep_pair)
    -> std::vector<CScreenedBasisFunctionPair>
{
    std::vector<CScreenedBasisFunctionPair> screened_pairs;

    for (const auto &[bra_atom_basis, ket_atom_basis] : basis.unique_basis_pairs())
    {
        const auto bra_coords = atom_coordinates(basis, molecule, bra_atom_basis);

        const auto bra_atoms = atom_indices(basis, bra_atom_basis);

        const auto bra_functions = bra_atom_basis.basis_functions();

        const auto same_atom_basis =
            (bra_atom_basis.get_name() == ket_atom_basis.get_name()) && (bra_atom_basis.get_identifier() == ket_atom_basis.get_identifier());

        if (same_atom_basis)
        {
            const auto nfuncs = bra_functions.size();

            for (size_t i = 0; i < nfuncs; i++)
            {
                for (size_t j = i; j < nfuncs; j++)
                {
                    if (i == j)
                    {
                        // identical function on the identical atom set: triangular atom pairs

                        screened_pairs.push_back(
                            CScreenedBasisFunctionPair(bra_functions[i], static_cast<int>(i), bra_coords, bra_atoms, keep_pair));
                    }
                    else
                    {
                        // distinct functions on the same atom set: full atom product, shared coordinates

                        screened_pairs.push_back(CScreenedBasisFunctionPair(
                            bra_functions[i], static_cast<int>(i), bra_coords, bra_atoms, bra_functions[j], static_cast<int>(j), bra_coords, bra_atoms, keep_pair));
                    }
                }
            }
        }
        else
        {
            const auto ket_coords = atom_coordinates(basis, molecule, ket_atom_basis);

            const auto ket_atoms = atom_indices(basis, ket_atom_basis);

            const auto ket_functions = ket_atom_basis.basis_functions();

            for (size_t i = 0; i < bra_functions.size(); i++)
            {
                for (size_t j = 0; j < ket_functions.size(); j++)
                {
                    screened_pairs.push_back(CScreenedBasisFunctionPair(
                        bra_functions[i], static_cast<int>(i), bra_coords, bra_atoms, ket_functions[j], static_cast<int>(j), ket_coords, ket_atoms, keep_pair));
                }
            }
        }
    }

    return screened_pairs;
}

}  // namespace gtofunc

#endif /* GtoFunc_hpp */
