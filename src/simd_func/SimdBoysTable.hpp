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


#ifndef SimdBoysTable_hpp
#define SimdBoysTable_hpp

namespace simdfunc {  // simdfunc namespace

/// @brief Gets the highest order of Boys function which can be computed.
/// @return The highest order of Boys function.
constexpr auto
max_boys_order() -> int
{
    return 32;
}

/// @brief Gets the number of grid points of the Taylor expansion of Boys function.
/// @return The number of grid points.
/// @note The grid spans the arguments from zero to thirty in steps of one tenth.
/// The upward recursion used above the grid loses accuracy for the highest orders
/// where the values of Boys function approach the exponential term subtracted in
/// it, which happens below an argument of twenty three for order thirty two, so
/// the grid extends well past it.
constexpr auto
boys_grid_points() -> int
{
    return 351;
}

/// @brief Gets the number of Taylor expansion coefficients of Boys function.
/// @return The number of coefficients.
constexpr auto
boys_coefficients() -> int
{
    return 8;
}

/// @brief Gets the Taylor expansion coefficients of Boys function on the grid of
/// arguments.
/// @return The pointer to the coefficients, indexed as
/// table[(order * boys_grid_points() + point) * boys_coefficients() + index].
auto boys_table() -> const double *;

/// @brief Gets the factors of the downward recursion of Boys function, i.e. the
/// reciprocals of the odd numbers.
/// @return The pointer to the factors, indexed by the order of Boys function.
auto boys_factors() -> const double *;

}  // namespace simdfunc

#endif /* SimdBoysTable_hpp */
