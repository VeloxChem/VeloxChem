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



#include "SimdHarmonicsRecP.hpp"

#include <cmath>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_p_solid_harmonics(const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_p_solid_harmonics: Coordinates must have seven rows"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(3, ncols);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);
    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    auto *s_m1 = matrix.data(0);
    auto *s_z0 = matrix.data(1);
    auto *s_p1 = matrix.data(2);

    // NOTE: the real solid harmonics of angular momentum one are the components
    // of the vector between the atoms, taken in the order of the spherical
    // components, and need no normalization factor.

    // NOTE: the rows of the coordinates and of the harmonics start at a cache
    // line boundary, so the loop is vectorized with aligned loads and stores.

#pragma omp simd aligned(s_m1, s_z0, s_p1, a_x, a_y, a_z, b_x, b_y, b_z : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        s_m1[i] = a_y[i] - b_y[i];

        s_z0[i] = a_z[i] - b_z[i];

        s_p1[i] = a_x[i] - b_x[i];
    }

    return matrix;
}

}  // namespace simdfunc
