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
#include <vector>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

/// @brief Computes the two-center electron repulsion integral of two basis
/// functions centered on the same atom.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @return The two-center electron repulsion integral of the basis functions.
/// @note The Coulomb operator is spherically symmetric about the atom the basis
/// functions are centered on, so the integral is diagonal in the angular
/// components and does not depend on the component or on the position of the
/// atom, as the overlap and the kinetic energy are. Basis functions of different
/// angular momenta have no integral and give zero.
/// @note The kernels are not implemented yet and this returns zero.
auto one_center_electron_repulsion(const CBasisFunction &bra, const CBasisFunction &ket) -> double;

/// @brief Computes the two-center electron repulsion integrals of a combination
/// of basis functions over the atom pairs of a block.
/// @param values The values of the combination of basis functions to compute, as
/// one row of nvalues columns per pair of angular components, with the component
/// on bra side as the slowest running index.
/// @param nvalues The number of atom pairs of the block.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param harmonics The solid harmonics of the vectors between the atoms.
/// @param coordinates The coordinates of the atom pairs.
/// @note There is no screening threshold, unlike the overlap and the kinetic
/// energy. The Coulomb operator decays as the inverse of the interatomic
/// distance, so neither an atom pair nor a pair of primitives falls below a
/// threshold at any distance a molecule reaches, and every combination of basis
/// functions is computed over every atom pair of the block.
/// @note The kernels are not implemented yet and this sets the values to zero.
auto compute_electron_repulsion(double                         *values,
                                const size_t                    nvalues,
                                const CBasisFunction           &bra,
                                const CBasisFunction           &ket,
                                const std::vector<CSimdMatrix> &harmonics,
                                const CSimdMatrix              &coordinates) -> void;

}  // namespace simdt2ceri

#endif /* SimdTwoCenterElectronRepulsionFunc_hpp */
