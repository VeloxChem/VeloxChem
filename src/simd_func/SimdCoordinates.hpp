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


#ifndef SimdCoordinates_hpp
#define SimdCoordinates_hpp

#include "AtomBasisPairGroup.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "AtomBasisTripleSparsity.hpp"
#include "Molecule.hpp"
#include "Point.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the coordinates of the atom pairs of a sparsity pattern.
/// @param sparsity The sparsity pattern of the atom pairs.
/// @param molecule The molecule to take the atomic coordinates from.
/// @return The matrix of coordinates with seven rows and as many columns as there
/// are atom pairs, holding the coordinates of the atoms on bra side in rows zero
/// to two, of the atoms on ket side in rows three to five, and the squared
/// distance of the atom pair in row six.
/// @note The coordinates are given in atomic units. The padding of the rows is
/// left undefined, as the atom pairs beyond the number of atom pairs are not
/// part of the sparsity pattern.
auto make_coordinates(const CAtomBasisPairSparsity &sparsity, const CMolecule &molecule) -> CSimdMatrix;

/// @brief Creates the coordinates of the atom pairs of an atom basis pair group.
/// @param group The atom basis pair group to create the coordinates of.
/// @param molecule The molecule to take the atomic coordinates from.
/// @return The matrix of seven rows and as many columns as the atom pairs of the
/// group, with the coordinates of the atoms on bra side in rows zero to two, of
/// the atoms on ket side in rows three to five, and the squared distances of the
/// atom pairs in row six.
/// @note This is the overload for the integrals which are not screened, whose
/// blocks are atom basis pair groups and carry no sparsity pattern.
auto make_coordinates(const CAtomBasisPairGroup &group, const CMolecule &molecule) -> CSimdMatrix;

/// @brief Creates the coordinates of the atom pairs on a and b sides of a
/// sparsity pattern.
/// @param sparsity The sparsity pattern of the atom pairs.
/// @param molecule The molecule to take the atomic coordinates from.
/// @return The matrix of coordinates with seven rows and as many columns as there
/// are atom pairs, holding the coordinates of the atoms on a side in rows zero to
/// two, of the atoms on b side in rows three to five, and the squared distance of
/// the atom pair in row six.
/// @note The atoms on c side are not part of the matrix, as they form a separate
/// and shorter dimension than the atom pairs on a and b sides.
auto make_coordinates(const CAtomBasisTripleSparsity &sparsity, const CMolecule &molecule) -> CSimdMatrix;

/// @brief Creates the coordinates of the atoms on b side of a set of atom pairs
/// and of one atom on c side.
/// @param coordinates The coordinates of the atom pairs on a and b sides, of
/// seven rows as the factories above lay them out.
/// @param center The coordinates of the atom on c side, in atomic units.
/// @return The matrix of coordinates with seven rows and as many columns as the
/// atom pairs, holding the coordinates of the atoms on b side in rows zero to
/// two, of the atom on c side in rows three to five, and the squared distance
/// between them in row six.
/// @note The layout is the one the atom pairs are given, so the solid harmonics
/// of the vectors from the atoms on b side to the atom on c side are formed by
/// the same recursions, with the same convention of the first center minus the
/// second, and no factory of the harmonics needs an overload.
/// @note The coordinates of the atom on c side are the same in every column, as
/// one atom is paired with every atom pair. They are broadcast rather than
/// looked up, so the rows stay contiguous and the recursions load them with
/// aligned SIMD instructions like any other row.
auto make_coordinates(const CSimdMatrix &coordinates, const TPoint<double> &center) -> CSimdMatrix;

}  // namespace simdfunc

#endif /* SimdCoordinates_hpp */
