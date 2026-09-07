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



#ifndef SimdOverlapRecSLH_hpp
#define SimdOverlapRecSLH_hpp

#include <cstddef>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdovl {  // simdovl namespace

/// @brief Computes the overlap integrals of a combination of one basis function
/// of zero angular momentum and one of angular momentum five, in either order.
/// @param values The values of the combination of basis functions in the values
/// block of the sparsity pattern.
/// @param nvalues The number of values to compute, i.e. the number of atom pairs
/// surviving the screening of the combination of basis functions.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, as ten rows ordered by
/// ascending interatomic distance, holding the vector between the atoms in rows
/// six to eight and its squared length in row nine.
/// @param threshold The screening threshold of the integrals.
/// @note The angular half is the same for both orders, as the harmonic is the same
/// polynomial of the vector between the atoms either way. The orders differ only in
/// the prefactor, which is selected once for the whole combination.
/// @note One term survives the integration over the Gaussian product center, so the
/// buffer holds a single accumulator and the integrals of the angular components are
/// formed straight into the values.
auto compute_slh_overlap(double               *values,
                         const size_t          nvalues,
                         const CBasisFunction &bra,
                         const CBasisFunction &ket,
                         const CSimdMatrix    &coordinates,
                         const double          threshold) -> void;

}  // namespace simdovl

#endif /* SimdOverlapRecSLH_hpp */
