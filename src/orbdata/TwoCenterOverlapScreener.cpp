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

#include "TwoCenterOverlapScreener.hpp"

#include <algorithm>
#include <cmath>

#include "MathConst.hpp"

CTwoCenterOverlapScreener::CTwoCenterOverlapScreener(const CBasisFunction &bra_function, const CBasisFunction &ket_function, const double threshold)

    : _factor(0.0)

    , _exponent(0.0)

    , _total_angular_momentum(bra_function.get_angular_momentum() + ket_function.get_angular_momentum())

    , _threshold(threshold)
{
    // smallest exponent and its coefficient are the last primitive on each side
    // (basis functions are ordered from largest to smallest exponent)

    const auto a_min = bra_function.get_exponents().back();

    const auto b_min = ket_function.get_exponents().back();

    const auto c_a = bra_function.get_normalization_factors().back();

    const auto c_b = ket_function.get_normalization_factors().back();

    const auto fpi = mathconst::pi_value();

    const auto sab = a_min + b_min;

    _factor = (fpi / sab) * std::sqrt(fpi / sab) * std::abs(c_a) * std::abs(c_b);

    _exponent = a_min * b_min / sab;
}

auto
CTwoCenterOverlapScreener::operator()(const TPoint<double> &bra_center, const TPoint<double> &ket_center) const -> bool
{
    // squared distance between the two centers

    const auto r_a = bra_center.coordinates();

    const auto r_b = ket_center.coordinates();

    const auto rx = r_a[0] - r_b[0];

    const auto ry = r_a[1] - r_b[1];

    const auto rz = r_a[2] - r_b[2];

    const auto r2 = rx * rx + ry * ry + rz * rz;

    // angular factor R^{la+lb} by repeated multiplication, floored at 1.0

    const auto r = std::sqrt(r2);

    auto rpow = 1.0;

    for (int k = 0; k < _total_angular_momentum; k++)
    {
        rpow *= r;
    }

    const auto angular = std::max(1.0, rpow);

    // largest overlap contribution estimate

    const auto estimate = _factor * angular * std::exp(-_exponent * r2);

    return estimate >= _threshold;
}

auto
CTwoCenterOverlapScreener::threshold() const -> double
{
    return _threshold;
}
