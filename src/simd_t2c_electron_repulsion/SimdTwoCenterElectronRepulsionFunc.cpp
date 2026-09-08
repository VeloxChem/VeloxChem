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



#include "SimdTwoCenterElectronRepulsionFunc.hpp"

#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"

namespace simderi {  // simderi namespace

auto
one_center_electron_repulsion(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if (lbra != lket) return 0.0;

    // NOTE: the atoms meet, so the argument of the Boys function is zero and the
    // solid harmonic of the vector between them is one for the angular momentum
    // zero and vanishes above it. What is left is a closed formula in the
    // exponents alone, diagonal in the angular components and independent of
    // which component it is.

    constexpr auto fpi = mathconst::pi_value();

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    // NOTE: the double factorial of twice the angular momentum less one over
    // twice it plus one is the factor the angular momentum contributes, and is
    // one for the angular momentum zero.

    const auto fang = mathfunc::double_factorial(2 * lbra - 1) / static_cast<double>(2 * lbra + 1);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    auto fsum = 0.0;

    for (size_t i = 0; i < a_exps.size(); i++)
    {
        for (size_t j = 0; j < b_exps.size(); j++)
        {
            const auto fexp = a_exps[i] + b_exps[j];

            auto fden = a_exps[i] * b_exps[j] * std::sqrt(fexp);

            // NOTE: twice the sum of the exponents is raised to the angular
            // momentum by repeated multiplication, as the angular momentum is
            // small.

            for (int l = 0; l < lbra; l++)
            {
                fden *= 2.0 * fexp;
            }

            fsum += a_norms[i] * b_norms[j] / fden;
        }
    }

    return fcoul * fang * fsum;
}

auto
compute_electron_repulsion(double               *values,
                           const size_t          nvalues,
                           const CBasisFunction &bra,
                           const CBasisFunction &ket,
                           const CSimdMatrix    &coordinates) -> void
{
    // NOTE: the kernels of the atom pairs are generated elsewhere and are not in the
    // tree yet, so no combination of basis functions is computed here. The
    // combination stops rather than leaving the values of the matrix unwritten,
    // which is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(
        false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_electron_repulsion: Integrals are not implemented"));
}

}  // namespace simderi
