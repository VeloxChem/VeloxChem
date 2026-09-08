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

#include <array>

#include "TensorComponents.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_electron_repulsion(double               *values,
                           const size_t          npairs,
                           const size_t          natoms,
                           const CBasisFunction &a_function,
                           const CBasisFunction &b_function,
                           const CBasisFunction &c_function,
                           const CSimdMatrix    &coordinates,
                           const CSimdMatrix    &c_coordinates) -> void
{
    // NOTE: the kernels of the atom pairs are generated elsewhere and are not in the
    // tree yet. Rather than stopping, this writes a value which encodes the position
    // of every element it is responsible for, so that the layout of the values and
    // the loops of the driver can be checked before the kernels exist. A caller which
    // reads these as integrals gets numbers which are obviously not integrals.

    const auto ncomps_a =
        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{a_function.get_angular_momentum()}));

    const auto ncomps_b =
        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{b_function.get_angular_momentum()}));

    const auto ncomps_c =
        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{c_function.get_angular_momentum()}));

    for (size_t ma = 0; ma < ncomps_a; ma++)
    {
        for (size_t mb = 0; mb < ncomps_b; mb++)
        {
            for (size_t mc = 0; mc < ncomps_c; mc++)
            {
                const auto icomp = (ma * ncomps_b + mb) * ncomps_c + mc;

                for (size_t ic = 0; ic < natoms; ic++)
                {
                    auto *row = values + (icomp * natoms + ic) * npairs;

                    for (size_t k = 0; k < npairs; k++)
                    {
                        // NOTE: the components, the atom on c side and the atom pair
                        // are packed into one number, so that a value read back names
                        // the place it was written to.

                        row[k] = static_cast<double>(icomp) * 1.0e6 + static_cast<double>(ic) * 1.0e3 +
                                 static_cast<double>(k);
                    }
                }
            }
        }
    }
}

}  // namespace simdt3ceri
