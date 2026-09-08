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
#include <string>

#include "ErrorHandler.hpp"
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
    // tree yet. Rather than stopping, this fills every element the block hands it
    // with a value which encodes that element's position, so that the loops of the
    // driver and the offsets of the tensor can be checked before the kernels exist.
    // A caller which reads these as integrals gets numbers which are obviously not
    // integrals.

    // NOTE: this walks the values as one range and decomposes the index only because
    // the encoding has to name a position. It is not the shape of a kernel: a kernel
    // loops over the pairs of primitives, builds a buffer, contracts it and writes
    // the values once through a transform. Nothing here anticipates that, and the
    // only contract it keeps is to fill the elements it was handed.

    const auto ncomps_a =
        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{a_function.get_angular_momentum()}));

    const auto ncomps_b =
        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{b_function.get_angular_momentum()}));

    const auto ncomps_c =
        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{c_function.get_angular_momentum()}));

    // NOTE: the three indices are packed into one number by decimal place, so a
    // block wider than a thousand atom pairs or carrying more than a thousand atoms
    // on c side would let them collide and the check would stop distinguishing
    // positions. It stops here rather than reporting agreement it did not establish.

    errors::assertMsgCritical(
        (npairs < 1000) && (natoms < 1000),
        std::string("SimdThreeCenterElectronRepulsionFunc.compute_electron_repulsion: The stub cannot encode a block this large"));

    const auto nvalues = ncomps_a * ncomps_b * ncomps_c * natoms * npairs;

    for (size_t n = 0; n < nvalues; n++)
    {
        const auto k = n % npairs;

        const auto ic = (n / npairs) % natoms;

        const auto icomp = n / (npairs * natoms);

        values[n] = static_cast<double>(icomp) * 1.0e6 + static_cast<double>(ic) * 1.0e3 + static_cast<double>(k);
    }
}

}  // namespace simdt3ceri
