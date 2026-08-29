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

#include "AtomBasisPairSparsity.hpp"
#include "Molecule.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the coordinates of the atom pairs of a sparsity pattern.
/// @param sparsity The sparsity pattern of the atom pairs.
/// @param molecule The molecule to take the atomic coordinates from.
/// @return The matrix of coordinates with six rows and as many columns as there
/// are atom pairs, holding the coordinates of the atoms on bra side in rows zero
/// to two and of the atoms on ket side in rows three to five.
/// @note The coordinates are given in atomic units. The padding of the rows is
/// left undefined, as the atom pairs beyond the number of atom pairs are not
/// part of the sparsity pattern.
auto make_coordinates(const CAtomBasisPairSparsity &sparsity, const CMolecule &molecule) -> CSimdMatrix;

}  // namespace simdfunc

#endif /* SimdCoordinates_hpp */
