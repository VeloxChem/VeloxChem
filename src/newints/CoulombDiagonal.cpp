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

#include "CoulombDiagonal.hpp"

#include <cmath>
#include <cstddef>

#include "BasisFunction.hpp"
#include "MathConst.hpp"

namespace newints {

auto
coulomb_concentric_value(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    const auto l = bra.get_angular_momentum();  // == ket angular momentum (caller guarantees l_a == l_b)

    const auto &exps_a = bra.exponents();

    const auto &coefs_a = bra.normalization_factors();

    const auto &exps_b = ket.exponents();

    const auto &coefs_b = ket.normalization_factors();

    // (2l - 1)!! with (-1)!! = 1
    auto dfact = 1.0;

    for (int k = 2 * l - 1; k > 0; k -= 2) dfact *= static_cast<double>(k);

    const auto two_l = static_cast<double>(1 << l);

    const auto pi = mathconst::pi_value();

    const auto two_pi52 = 2.0 * pi * pi * std::sqrt(pi);

    // 2 pi^{5/2} (2l-1)!! / ((2l+1) 2^l)
    const auto prefac = two_pi52 * dfact / (static_cast<double>(2 * l + 1) * two_l);

    auto vab = 0.0;

    for (std::size_t i = 0; i < exps_a.size(); i++)
    {
        for (std::size_t j = 0; j < exps_b.size(); j++)
        {
            const auto alpha = exps_a[i];

            const auto beta = exps_b[j];

            const auto p = alpha + beta;

            // p^l by repeated multiplication (l = 0..6), then 1 / (alpha beta p^l sqrt(p))
            auto pl = 1.0;

            for (int k = 0; k < l; k++) pl *= p;

            const auto denom = alpha * beta * pl * std::sqrt(p);

            vab += coefs_a[i] * coefs_b[j] * prefac / denom;
        }
    }

    return vab;
}

}  // namespace newints
