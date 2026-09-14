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


#include "SimdNuclearPotentialFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"

#include "SimdNuclearPotentialRecSS.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_nuclear_potential(double                    *values,
                          const size_t               nvalues,
                          const CBasisFunction      &bra,
                          const CBasisFunction      &ket,
                          const CSimdMatrix         &coordinates,
                          const std::vector<double> &charges,
                          const std::vector<double> &points,
                          CSimdMatrix               &buffer,
                          const double               threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);

        return;
    }

    // NOTE: falling through to nothing would leave the values as the caller found
    // them and say nothing, and a caller cannot tell integrals which were not
    // computed from integrals which are zero. The kernels are added one pair of
    // angular momenta at a time and this list grows with them.

    errors::assertMsgCritical(false,
                              std::string("compute_nuclear_potential: No kernel for the combination of angular momenta ") +
                                  std::to_string(lbra) + std::string(" and ") + std::to_string(lket) +
                                  std::string("; only S S is written"));
}

}  // namespace simdnpot
