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

#ifndef newints_CoulombDiagonal_hpp
#define newints_CoulombDiagonal_hpp

class CBasisFunction;

namespace newints {

/// @brief Same-center (concentric, R_AB = 0) two-center Coulomb value for two shells
/// of equal angular momentum l.
///
/// On a common center the Coulomb block is diagonal in (l, m) and independent of m.
/// For concentric solid-harmonic Gaussians the momentum-space Coulomb integral gives,
/// in VeloxChem's unnormalized (Racah, int Y^2 dOmega = 4 pi / (2l+1)) convention, the
/// per-primitive-pair value
///   2 pi^{5/2} (2l-1)!! / ((2l+1) 2^l) * 1 / (alpha beta p^l sqrt(p)),  p = alpha + beta,
/// summed over the contraction. (Reduces to 2 pi^{5/2} / (alpha beta sqrt(p)) for l = 0,
/// matching the (s|s) seed.) The two shells must share angular momentum (the caller
/// guarantees l_a == l_b).
///
/// This single value is the (m, m) diagonal entry of the block, used both by the
/// ElectronRepulsionDriver same-atom blocks and by the per-shell Cauchy-Schwarz factor
/// Q_P = sqrt((P|P)) in CauchySchwarzData.
auto coulomb_concentric_value(const CBasisFunction &bra, const CBasisFunction &ket) -> double;

}  // namespace newints

#endif /* newints_CoulombDiagonal_hpp */
