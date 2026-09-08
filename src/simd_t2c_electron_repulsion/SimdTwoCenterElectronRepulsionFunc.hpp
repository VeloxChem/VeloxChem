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



#ifndef SimdTwoCenterElectronRepulsionFunc_hpp
#define SimdTwoCenterElectronRepulsionFunc_hpp

#include <cstddef>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simderi {  // simderi namespace

/// @brief Computes the two-center electron repulsion integral of two basis
/// functions centered on the same atom.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @return The two-center electron repulsion integral of the basis functions.
/// @note The Coulomb operator is spherically symmetric about the atom the basis
/// functions are centered on, so the integral is diagonal in the angular components
/// and does not depend on the component or on the position of the atom, as the
/// overlap and the kinetic energy are. Basis functions of different angular momenta
/// have no integral and give zero.
/// @note The integral is the value the integrals of two atoms approach as the atoms
/// meet, which the kernels of the atom pairs do not carry: they hold the solid
/// harmonic of the vector between the atoms, which vanishes there. It is a closed
/// formula in the exponents alone,
///
///   (l|J|l')_{m,m'} = 2 pi^(5/2) (2l - 1)!! / (alpha beta (2l + 1) (2p)^l sqrt(p))
///
/// for equal angular momenta and equal components, and zero otherwise.
auto one_center_electron_repulsion(const CBasisFunction &bra, const CBasisFunction &ket) -> double;

/// @brief Computes the two-center electron repulsion integrals of a combination of
/// basis functions over the atom pairs of a block.
/// @param values The values of the combination of basis functions, as one row of
/// nvalues columns for each of the (2 l_bra + 1) (2 l_ket + 1) spherical components,
/// with the components of the bra side running slowest.
/// @param nvalues The number of atom pairs of the block.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, as ten rows ordered by
/// ascending interatomic distance, holding the vector between the atoms in rows six
/// to eight and its squared length in row nine.
/// @note There is no screening threshold here, unlike the overlap and the kinetic
/// energy. The Coulomb operator decays as the inverse of the interatomic distance,
/// so neither an atom pair nor a pair of primitives falls below a threshold at any
/// distance a molecule reaches, and every combination of basis functions is computed
/// over every atom pair of every block.
/// @note The kernels of the atom pairs are not written yet, so every combination
/// stops with an error rather than returning values it did not compute. The integral
/// of two basis functions on the same atom is unaffected and is computed by
/// one_center_electron_repulsion above.
auto compute_electron_repulsion(double               *values,
                                const size_t          nvalues,
                                const CBasisFunction &bra,
                                const CBasisFunction &ket,
                                const CSimdMatrix    &coordinates) -> void;

}  // namespace simderi

#endif /* SimdTwoCenterElectronRepulsionFunc_hpp */
