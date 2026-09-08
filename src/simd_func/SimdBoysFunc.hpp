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



#ifndef SimdBoysFunc_hpp
#define SimdBoysFunc_hpp

#include <cstddef>
#include <initializer_list>

#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Computes the values of Boys function of every order up to the requested
/// one, for one pair of primitives.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to. The values follow
/// it, so the routine writes order + 2 rows in all.
/// @param order The highest order of Boys function to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fj The factor every value is scaled by.
/// @param mu The factor the squared distance is scaled by to give the argument.
/// @note The row at target holds the argument, mu times the squared distance, and
/// the rows from target + 1 to target + order + 1 hold the values of order zero to
/// order, each scaled by fj. The argument is kept rather than recomputed, as the
/// recursions below read it once per order.
/// @note The values are computed from the Taylor expansion on a grid of arguments
/// below the end of the grid and from the asymptotic expansion above it, with the
/// remaining orders following from the downward and the upward recursion
/// respectively. The relative accuracy is better than 1.0e-14 for every order up to
/// the highest supported one and for every argument.
auto compute_full_boys_function(CSimdMatrix       &buffer,
                                const CSimdMatrix &coordinates,
                                const size_t       target,
                                const size_t       order,
                                const size_t       ncols,
                                const double       fj,
                                const double       mu) -> void;

/// @brief Computes the values of Boys function of the requested orders alone, for
/// one pair of primitives.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to. The values follow
/// it, one row per requested order and in the order they are requested, so the
/// routine writes one row more than there are orders.
/// @param orders The orders of Boys function to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fj The factor every value is scaled by.
/// @param mu The factor the squared distance is scaled by to give the argument.
/// @note A recursion carries the values from the order the expansion gives to the
/// order which is wanted, so the orders below the highest requested one are computed
/// whether they are asked for or not. What this form saves is the rows of the
/// buffer, not the arithmetic: a kernel which needs two orders out of twelve holds
/// three rows rather than fourteen.
auto compute_boys_function(CSimdMatrix                        &buffer,
                           const CSimdMatrix                  &coordinates,
                           const size_t                        target,
                           const std::initializer_list<size_t> orders,
                           const size_t                        ncols,
                           const double                        fj,
                           const double                        mu) -> void;

}  // namespace simdfunc

#endif /* SimdBoysFunc_hpp */
