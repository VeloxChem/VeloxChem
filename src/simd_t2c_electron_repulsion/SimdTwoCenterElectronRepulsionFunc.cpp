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

#include <algorithm>
#include <string>

#include "ErrorHandler.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
one_center_electron_repulsion(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    if (bra.get_angular_momentum() != ket.get_angular_momentum()) return 0.0;

    // NOTE: the kernels are not implemented yet. The dispatch over the angular
    // momenta of the basis functions replaces this.

    return 0.0;
}

auto
compute_electron_repulsion(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false,
            std::string("SimdTwoCenterElectronRepulsionFunc.compute_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the geometry of the integrals reaches the sum of the angular momenta
    // of the basis functions, so the harmonics of the block must reach it too.

    if (const auto lsum = lbra + lket; (lsum > 0) && (harmonics.size() < static_cast<size_t>(lsum)))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_electron_repulsion: Harmonics do not reach the angular momenta"));
    }

    // NOTE: the kernels are not implemented yet. The dispatch over the angular
    // momenta of the basis functions replaces this, and the values are set to
    // zero until then so that the matrix the driver returns is defined.

    const auto ncomps = static_cast<size_t>((2 * lbra + 1) * (2 * lket + 1));

    std::fill(values, values + ncomps * nvalues, 0.0);
}

}  // namespace simdt2ceri
