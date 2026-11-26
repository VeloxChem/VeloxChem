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


#ifndef SimdTwoCenterElectronRepulsionBufferRows_hpp
#define SimdTwoCenterElectronRepulsionBufferRows_hpp

#include <array>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

/// @brief Gets the number of rows of the buffer a combination of basis functions
/// needs.
/// @param bra_angular_momentum The angular momentum of basis function on bra side.
/// @param ket_angular_momentum The angular momentum of basis function on ket side.
/// @return The number of rows.
/// @note The buffer belongs to the block and not to the combination, so a caller
/// sizes it once from the highest angular momenta the block carries. The numbers
/// below do not decrease with either angular momentum, so the largest combination of
/// a block is the one of its highest momenta and no other combination needs more.
/// @note This table is written with the kernels and describes them. Regenerating the
/// kernels without regenerating it leaves a buffer which is too small, which the
/// assertions of simdfunc::prepare_buffer report rather than let pass.
inline auto
number_of_buffer_rows(const int bra_angular_momentum, const int ket_angular_momentum) -> size_t
{
    constexpr std::array<std::array<size_t, 9>, 9> rows{{
        {       2,      11,      22,      42,      76,     129,     208,     320,     474},
        {      11,      39,      85,     153,     254,     397,     592,     850,    1183},
        {      22,      85,     207,     418,     775,    1328,    2151,    3314,    4911},
        {      42,     159,     395,     812,    1525,    2632,    4285,    6622,    9835},
        {      76,     269,     836,    1741,    3246,    5556,    8940,   13673,   20094},
        {     129,     424,    1393,    2912,    5436,    9298,   14942,   22813,   33467},
        {     208,     634,    2330,    4881,    9074,   15450,   24691,   37506,   54745},
        {     320,     910,    3530,    7395,   13728,   23322,   37181,   56338,   82037},
        {     474,    1264,    5277,   11050,   20447,   34618,   54980,   83013,  120464}
    }};

    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum < 9) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum < 9),
                              std::string("SimdTwoCenterElectronRepulsionBufferRows.number_of_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(bra_angular_momentum)][static_cast<size_t>(ket_angular_momentum)];
}

}  // namespace simdt2ceri

#endif /* SimdTwoCenterElectronRepulsionBufferRows_hpp */
