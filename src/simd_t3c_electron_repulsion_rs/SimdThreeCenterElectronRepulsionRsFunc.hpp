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




#ifndef SimdThreeCenterElectronRepulsionRsFunc_hpp
#define SimdThreeCenterElectronRepulsionRsFunc_hpp

#include <cstddef>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

/// @brief Computes the Coulomb and the range separated Coulomb three-center
/// integrals of a combination of basis functions over the atom pairs and the atoms
/// on c side of a block.
/// @param values The values of the combination. The block of the Coulomb integrals
/// comes first and the block of the range separated ones follows it, each of
/// (2 la + 1) (2 lb + 1) (2 lc + 1) rows of natoms times npairs columns.
/// @param npairs The number of atom pairs surviving the screening of the combination.
/// @param natoms The number of atoms on c side of the block.
/// @param a_function The basis function on a side.
/// @param b_function The basis function on b side.
/// @param c_function The basis function on c side.
/// @param coordinates The coordinates of the atom pairs.
/// @param c_coordinates The coordinates of the atoms on c side.
/// @param buffer The buffer of the combination of basis functions.
/// @param omega The range separation parameter.
/// @param threshold The screening threshold.
/// @note The two operators are computed together and not one after the other. They
/// share the pattern, the coordinates, the recurrence, the primitives and the
/// transformation, and differ only in the Boys function each of them carries.
/// @note The screening is the Coulomb bound and is applied to both. The attenuated
/// operator is bounded by the plain one everywhere, erf(omega r) / r being at most
/// 1 / r, so a pair the Coulomb bound keeps is a pair the attenuated integrals
/// cannot need more of, and a pair it drops carries nothing in either.
auto compute_rs_electron_repulsion(double               *values,
                                   const size_t          npairs,
                                   const size_t          natoms,
                                   const CBasisFunction &a_function,
                                   const CBasisFunction &b_function,
                                   const CBasisFunction &c_function,
                                   const CSimdMatrix    &coordinates,
                                   const CSimdMatrix    &c_coordinates,
                                   CSimdMatrix          &buffer,
                                   const double          omega,
                                   const double          threshold) -> void;

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionRsFunc_hpp */
