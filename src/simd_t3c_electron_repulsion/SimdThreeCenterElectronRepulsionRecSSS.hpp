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



#ifndef SimdThreeCenterElectronRepulsionRecSSS_hpp
#define SimdThreeCenterElectronRepulsionRecSSS_hpp

#include <cstddef>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

/// @brief Computes the three-center electron repulsion integrals of three basis
/// functions of zero angular momentum, over the atom pairs of a block and for one
/// atom on c side.
/// @param values The values of the combination of basis functions, whose slice of
/// the atom on c side this kernel writes.
/// @param npairs The number of surviving atom pairs of the combination.
/// @param natoms The number of atoms on c side of the block.
/// @param iatom The index of the atom on c side among the atoms of the block.
/// @param a_function The basis function on a side.
/// @param b_function The basis function on b side.
/// @param c_function The basis function on c side.
/// @param ab_coordinates The coordinates of the atom pairs, whose rows zero to
/// two carry the atoms on a side and whose row six carries their squared
/// distances.
/// @param bc_coordinates The coordinates of the atoms on b side in rows zero to
/// two and of the atom on c side in rows three to five.
/// @param threshold The screening threshold of the integrals.
/// @note The three basis functions carry no angular components, so the solid
/// harmonics of the block do not enter and the combination holds one value for
/// each atom pair and atom on c side.
auto compute_sss_electron_repulsion(double               *values,
                                    const size_t          npairs,
                                    const size_t          natoms,
                                    const size_t          iatom,
                                    const CBasisFunction &a_function,
                                    const CBasisFunction &b_function,
                                    const CBasisFunction &c_function,
                                    const CSimdMatrix    &ab_coordinates,
                                    const CSimdMatrix    &bc_coordinates,
                                    const double          threshold) -> void;

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionRecSSS_hpp */
