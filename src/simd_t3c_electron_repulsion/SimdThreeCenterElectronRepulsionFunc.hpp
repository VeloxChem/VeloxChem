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



#ifndef SimdThreeCenterElectronRepulsionFunc_hpp
#define SimdThreeCenterElectronRepulsionFunc_hpp

#include <cstddef>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

/// @brief Computes the three-center electron repulsion integrals of a combination of
/// basis functions over the atom pairs and the atoms on c side of a block.
/// @param values The values of the combination of basis functions to compute.
/// @param npairs The number of atom pairs surviving the screening of the combination.
/// @param natoms The number of atoms on c side of the block.
/// @param a_function The basis function on a side.
/// @param b_function The basis function on b side.
/// @param c_function The basis function on c side.
/// @param coordinates The coordinates of the atom pairs on a and b sides, as
/// SimdCoordinates lays them out, ten rows ordered by ascending interatomic distance.
/// @param c_coordinates The coordinates of the atoms on c side, as three rows of as
/// many columns as there are atoms.
/// @note The atom pairs run fastest and are contiguous, the atom on c side runs next
/// and the triple of angular components runs slowest, with the component on a side
/// the slowest of the three and the one on c side the fastest. The value of the atom
/// pair k, the atom ic on c side and the components ma, mb and mc is therefore at
/// (((ma * ncomps_b + mb) * ncomps_c + mc) * natoms + ic) * npairs + k.
/// @note The atoms on c side are a separate and shorter dimension than the atom
/// pairs, so the values of one triple of angular components hold the atom pairs
/// contiguously and repeat them for every atom on c side.
/// @param threshold The screening threshold of the integrals.
/// @note The threshold screens the triples of primitives here, as the blocks were
/// screened with it on the atom pairs. The atom on c side does not screen, as the
/// Coulomb operator decays as the inverse of the distance to it, so the surviving
/// triples are the same for every atom on that side and are counted once.
/// @note Only the combination of three basis functions of zero angular momentum is
/// implemented; the remaining kernels are generated elsewhere and are not in the tree
/// yet, so they stop with an error rather than returning values they did not
/// compute.
auto compute_electron_repulsion(double               *values,
                                const size_t          npairs,
                                const size_t          natoms,
                                const CBasisFunction &a_function,
                                const CBasisFunction &b_function,
                                const CBasisFunction &c_function,
                                const CSimdMatrix    &coordinates,
                                const CSimdMatrix    &c_coordinates,
                                CSimdMatrix          &buffer,
                                const double          threshold) -> void;

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionFunc_hpp */
