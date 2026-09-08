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



#include "SimdThreeCenterElectronRepulsionFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_electron_repulsion(double               *values,
                           const size_t          npairs,
                           const size_t          natoms,
                           const CBasisFunction &a_function,
                           const CBasisFunction &b_function,
                           const CBasisFunction &c_function,
                           const CSimdMatrix    &coordinates,
                           const CSimdMatrix    &c_coordinates,
                           const double          threshold) -> void
{
    const auto la = a_function.get_angular_momentum();

    const auto lb = b_function.get_angular_momentum();

    const auto lc = c_function.get_angular_momentum();

    if ((la == 0) && (lb == 0) && (lc == 0))
    {
        compute_sss_three_center_electron_repulsion(
            values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, threshold);

        return;
    }

    // NOTE: the remaining kernels are generated elsewhere and are not in the tree
    // yet. The combination stops rather than leaving the values of the tensor
    // unwritten, which is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(
        false, std::string("SimdThreeCenterElectronRepulsionFunc.compute_electron_repulsion: Integrals are not implemented"));
}

}  // namespace simdt3ceri
