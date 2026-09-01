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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_electron_repulsion(double                         *values,
                           const size_t                    npairs,
                           const size_t                    natoms,
                           const size_t                    iatom,
                           const CBasisFunction           &a_function,
                           const CBasisFunction           &b_function,
                           const CBasisFunction           &c_function,
                           const std::vector<CSimdMatrix> &ab_harmonics,
                           const std::vector<CSimdMatrix> &bc_harmonics,
                           const CSimdMatrix              &ab_coordinates,
                           const CSimdMatrix              &bc_coordinates) -> void
{
    // NOTE: this is the stub of the skeleton and not the integral. It writes the
    // position of every element within the values of the combination, so that
    // the sizes and the offsets the sparsity pattern lays out are checked before
    // the kernels exist: the values of a combination must read 0, 1, 2 and so on
    // over its whole extent once every atom on c side has been visited, which no
    // two combinations can do if their storage overlaps or is sized wrongly.

    const auto ncomps_a = static_cast<size_t>(2 * a_function.get_angular_momentum() + 1);

    const auto ncomps_b = static_cast<size_t>(2 * b_function.get_angular_momentum() + 1);

    const auto ncomps_c = static_cast<size_t>(2 * c_function.get_angular_momentum() + 1);

    for (size_t ma = 0; ma < ncomps_a; ma++)
    {
        for (size_t mb = 0; mb < ncomps_b; mb++)
        {
            for (size_t mc = 0; mc < ncomps_c; mc++)
            {
                const auto offset = (((ma * ncomps_b + mb) * ncomps_c + mc) * natoms + iatom) * npairs;

                auto *prim = values + offset;

                for (size_t k = 0; k < npairs; k++)
                {
                    prim[k] = static_cast<double>(offset + k);
                }
            }
        }
    }
}

}  // namespace simdt3ceri
