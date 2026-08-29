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


#ifndef ScreeningFunc_hpp
#define ScreeningFunc_hpp

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <vector>

#include "BasisFunction.hpp"
#include "MathConst.hpp"

namespace screenfunc {  // screenfunc namespace

/// @brief Counts leading elements of range satisfying predicate.
/// @param values The range of values to count leading elements in.
/// @param predicate The predicate to partition range with.
/// @return The number of leading elements satisfying predicate.
/// @note The range must be partitioned with respect to predicate, i.e. all
/// elements satisfying predicate must precede those which do not. Bisection is
/// used to locate the partition point, so the predicate is evaluated a
/// logarithmic number of times. A range which is not partitioned yields a wrong
/// count instead of an error, as detecting it requires the linear scan avoided
/// here.
template <typename T, typename P>
inline auto
number_of_leading(const std::vector<T> &values, P predicate) -> size_t
{
    return static_cast<size_t>(std::ranges::distance(values.begin(), std::ranges::partition_point(values, predicate)));
}

/// @brief Computes upper bound of two-center overlap integral between basis
/// functions on bra and ket sides.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param r The distance between the atoms carrying the basis functions.
/// @return The upper bound of two-center overlap integral.
/// @note The bound rises with distance beyond unit distance before the Gaussian
/// decay takes over, so it is not monotone. It is nonetheless partitioned by any
/// threshold below its value at unit distance, which is the regime of screening
/// thresholds, and may thus be bisected with number_of_leading.
inline auto
two_center_overlap_bound(const CBasisFunction &bra, const CBasisFunction &ket, const double r) -> double
{
    const auto fexp = bra.smallest_exponent() + ket.smallest_exponent();

    const auto fmu = bra.smallest_exponent() * ket.smallest_exponent() / fexp;

    const auto fpi = mathconst::pi_value() / fexp;

    auto fact = bra.sum_of_absolute_norms() * ket.sum_of_absolute_norms();

    // NOTE: the distance prefactor is clamped to one below unit distance, where
    // it would lower the bound instead of raising it. Repeated multiplication is
    // used, as the total angular momentum is small.

    if (r > 1.0)
    {
        const auto lsum = bra.get_angular_momentum() + ket.get_angular_momentum();

        for (int i = 0; i < lsum; i++)
        {
            fact *= r;
        }
    }

    return fact * fpi * std::sqrt(fpi) * std::exp(-fmu * r * r);
}

}  // namespace screenfunc

#endif /* ScreeningFunc_hpp */
