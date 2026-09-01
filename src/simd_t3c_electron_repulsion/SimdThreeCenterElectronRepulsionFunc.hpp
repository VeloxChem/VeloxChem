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
#include <vector>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

/// @brief Computes the three-center electron repulsion integrals of a
/// combination of basis functions over the atom pairs of a block, for one atom
/// on c side.
/// @param values The values of the combination of basis functions to compute.
/// The atom pairs run fastest and are contiguous, the atom on c side runs next,
/// and the triple of angular components runs slowest, with the component on a
/// side the slowest of the three and the one on c side the fastest. The value of
/// the atom pair k and the components ma, mb and mc is therefore at
/// (((ma * ncomps_b + mb) * ncomps_c + mc) * natoms + iatom) * npairs + k, so
/// the atom pairs of one triple of components are contiguous and the triples are
/// natoms * npairs apart.
/// @param npairs The number of surviving atom pairs of the combination.
/// @param natoms The number of atoms on c side of the block.
/// @param iatom The index of the atom on c side among the atoms of the block.
/// @param a_function The basis function on a side.
/// @param b_function The basis function on b side.
/// @param c_function The basis function on c side.
/// @param ab_harmonics The solid harmonics of the vectors between the atoms of
/// the atom pairs, reaching the sum of the angular momenta on the a and b sides.
/// @param bc_harmonics The solid harmonics of the vectors from the atoms on b
/// side to the atom on c side, reaching the sum of the angular momenta on the b
/// and c sides.
/// @param ab_coordinates The coordinates of the atom pairs on a and b sides, as
/// SimdCoordinates lays them out.
/// @param bc_coordinates The vectors from the atoms on b side to the atom on c
/// side and their squared lengths, as SimdCoordinates lays them out.
/// @note One atom on c side is computed at a time, so the geometry of the three
/// centers is fixed for the whole call and the primitives can be contracted
/// before the angular components are formed.
/// @note The surviving atom pairs of a combination are a leading subrange of the
/// atom pairs of the block, which the ordering by interatomic distance and the
/// bisection of the screening give, so the coordinates and the harmonics of the
/// block are read over their first npairs columns and no gathering is needed.
/// @note This is the stub of the skeleton. It writes a value which encodes the
/// position of every element it is responsible for, so that the layout of the
/// values blocks and the loops of the driver are checked before the kernels
/// exist. It is not the integral.
auto compute_electron_repulsion(double                         *values,
                                const size_t                    npairs,
                                const size_t                    natoms,
                                const size_t                    iatom,
                                const CBasisFunction           &a_function,
                                const CBasisFunction           &b_function,
                                const CBasisFunction           &c_function,
                                const std::vector<CSimdMatrix> &ab_harmonics,
                                const std::vector<CSimdMatrix> &bc_harmonics,
                                const CSimdMatrix              &ab_coordinates,
                                const CSimdMatrix              &bc_coordinates) -> void;

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionFunc_hpp */
