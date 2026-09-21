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


#ifndef SimdTwoCenterElectronRepulsionGeom10RsBufferRows_hpp
#define SimdTwoCenterElectronRepulsionGeom10RsBufferRows_hpp

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

// NOTE: the name is qualified by which kernels the table describes. Every set in
// this namespace declares its own table as an inline function, so a shared name is
// one definition of several different things: the program is ill formed, the linker
// keeps whichever it saw first, and the tables differ both in their numbers and in
// their shape. The generator still emits the bare name; until it does not, a
// regenerated set has to be renamed again.
inline auto
number_of_geom_10_rs_buffer_rows(const int bra_angular_momentum, const int ket_angular_momentum) -> size_t
{
    constexpr std::array<std::array<size_t, 9>, 9> rows{{
        {      25,      72,     170,     318,     538,     848,    1268,    1820,    2528},
        {      65,     215,     525,    1051,    1903,    3191,    5059,    7657,   11169},
        {     133,     444,    1104,    2192,    3954,    6590,   10408,   15688,   22818},
        {     239,     796,    2106,    4236,    7656,   12776,   20134,   30280,   43892},
        {     393,    1289,    3473,    7009,   12693,   21191,   33391,   50183,   72679},
        {     609,    1949,    5655,   11471,   20767,   34625,   54409,   81537,  117709},
        {     901,    2804,    8418,   17116,   31014,   51698,   81176,  121514,  175200},
        {    1287,    3884,   12330,   25136,   45532,   75818,  118828,  177522,  255394},
        {    1785,    5221,   17197,   35105,   63597,  105829,  165689,  247213,  355169}
    }};

    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum < 9) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum < 9),
                              std::string("SimdTwoCenterElectronRepulsionGeom10RsBufferRows.number_of_geom_10_rs_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(bra_angular_momentum)][static_cast<size_t>(ket_angular_momentum)];
}

}  // namespace simdt2ceri

#endif /* SimdTwoCenterElectronRepulsionGeom10RsBufferRows_hpp */
