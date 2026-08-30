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



#ifndef SimdHarmonicsRecK_hpp
#define SimdHarmonicsRecK_hpp

#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the real solid harmonics of angular momentum seven of the
/// vectors between the atoms of the atom pairs.
/// @param i_harmonics The solid harmonics of angular momentum six.
/// @param h_harmonics The solid harmonics of angular momentum five.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 15 rows and as many columns as the coordinates, with
/// the row of index m + 7 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_k_solid_harmonics(const CSimdMatrix &i_harmonics, const CSimdMatrix &h_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

}  // namespace simdfunc

#endif /* SimdHarmonicsRecK_hpp */
