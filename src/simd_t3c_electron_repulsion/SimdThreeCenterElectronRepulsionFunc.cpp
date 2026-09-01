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
                           const CBasisFunction           &a_function,
                           const CBasisFunction           &b_function,
                           const CBasisFunction           &c_function,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const CSimdMatrix              &c_coordinates) -> void
{
    // NOTE: this is the stub of the skeleton and not the integral. It writes the
    // position of every element within the values of the combination, so that
    // the sizes and the offsets the sparsity pattern lays out are checked before
    // the kernels exist: the values of a combination must read 0, 1, 2 and so on
    // over its whole extent, which no two combinations can do if their storage
    // overlaps or is sized wrongly.

    const auto ncomps = static_cast<size_t>(2 * a_function.get_angular_momentum() + 1) *
                        static_cast<size_t>(2 * b_function.get_angular_momentum() + 1) *
                        static_cast<size_t>(2 * c_function.get_angular_momentum() + 1);

    const auto nvalues = ncomps * natoms * npairs;

    for (size_t i = 0; i < nvalues; i++)
    {
        values[i] = static_cast<double>(i);
    }
}

}  // namespace simdt3ceri
