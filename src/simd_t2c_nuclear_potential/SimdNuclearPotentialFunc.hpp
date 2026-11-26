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


#ifndef SimdNuclearPotentialFunc_hpp
#define SimdNuclearPotentialFunc_hpp

#include <cstddef>
#include <vector>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdnpot {  // simdnpot namespace

/// @brief Computes the nuclear potential integrals of a combination of basis
/// functions by dispatching to the kernel of their angular momenta.
/// @param values The values of the combination of basis functions in the values
/// block of the sparsity pattern.
/// @param nvalues The number of values to compute, i.e. the number of atom pairs
/// surviving the screening of the combination of basis functions.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, as ten rows ordered by
/// ascending interatomic distance, the last of which holds the squared distance of
/// the atom pair.
/// @param charges The magnitude of each point charge.
/// @param points The position of each point charge, as a flat array of three values
/// per charge.
/// @param buffer The buffer of the block, which spans every combination the block
/// carries and is sized from the highest angular momenta it holds.
/// @param threshold The screening threshold of the integrals.
/// @note The values of the combination of basis functions are stored as one row of
/// nvalues columns for each of the (2 l_bra + 1) (2 l_ket + 1) spherical components,
/// with the components of the bra side running slowest.
/// @note The kernels reach angular momentum six. A combination above it stops with
/// an error rather than returning the zeros of an unwritten kernel, which a caller
/// could not tell from integrals which are genuinely zero.
/// @note Unlike the overlap and the kinetic energy, there is no closed form for two
/// basis functions on the same atom: the operator is centred on the charges and not
/// on the atom, so the diagonal blocks are computed by these kernels as well.
auto compute_nuclear_potential(double                    *values,
                               const size_t               nvalues,
                               const CBasisFunction      &bra,
                               const CBasisFunction      &ket,
                               const CSimdMatrix         &coordinates,
                               const std::vector<double> &charges,
                               const std::vector<double> &points,
                               CSimdMatrix               &buffer,
                               const double               threshold) -> void;

}  // namespace simdnpot

#endif /* SimdNuclearPotentialFunc_hpp */
