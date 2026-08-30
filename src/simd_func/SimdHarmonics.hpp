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



#ifndef SimdHarmonics_hpp
#define SimdHarmonics_hpp

#include <vector>

#include "SimdMatrix.hpp"
#include "SimdHarmonicsRecP.hpp"
#include "SimdHarmonicsRecD.hpp"
#include "SimdHarmonicsRecF.hpp"
#include "SimdHarmonicsRecG.hpp"
#include "SimdHarmonicsRecH.hpp"
#include "SimdHarmonicsRecI.hpp"
#include "SimdHarmonicsRecK.hpp"
#include "SimdHarmonicsRecL.hpp"
#include "SimdHarmonicsRecM.hpp"
#include "SimdHarmonicsRecN.hpp"
#include "SimdHarmonicsRecO.hpp"
#include "SimdHarmonicsRecQ.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the real solid harmonics of angular momentum one of the vectors
/// between the atoms of the atom pairs.
/// @param coordinates The coordinates of the atom pairs, as seven rows holding
/// the coordinates of the atoms on bra side in rows zero to two and of the atoms
/// on ket side in rows three to five.
/// @return The matrix of three rows and as many columns as the coordinates, with
/// the row of index m + 1 holding the solid harmonic of order m.
/// @note The vector between the atoms is taken from the atom on ket side to the
/// atom on bra side, as the squared distance of the coordinates is, so the sign
/// of the solid harmonics of odd angular momenta follows that convention.
/// @note The solid harmonics are ordered and normalized as the spherical
/// components of a basis function of the same angular momentum, i.e. the row of
/// index zero holds y, the row of index one holds z and the row of index two
/// holds x. The padding of the rows is left undefined, as the atom pairs beyond
/// the number of columns are not part of the sparsity pattern.
/// @brief Creates the real solid harmonics of all angular momenta up to the given
/// one of the vectors between the atoms of the atom pairs.
/// @param coordinates The coordinates of the atom pairs, as seven rows.
/// @param lmax The highest angular momentum to create the harmonics of.
/// @return The vector of matrices, with the element of index l - 1 holding the
/// harmonics of angular momentum l as 2 l + 1 rows. The vector is empty for a
/// highest angular momentum below one, as the harmonic of angular momentum zero
/// is one for every atom pair and is not stored.
/// @note The harmonics of an angular momentum are formed from those of the two
/// angular momenta below it, so all of them are created on the way up and kept.
auto make_solid_harmonics(const CSimdMatrix &coordinates, const int lmax) -> std::vector<CSimdMatrix>;

}  // namespace simdfunc

#endif /* SimdHarmonics_hpp */
