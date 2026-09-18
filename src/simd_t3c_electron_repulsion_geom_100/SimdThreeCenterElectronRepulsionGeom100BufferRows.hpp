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


#ifndef SimdThreeCenterElectronRepulsionGeom100BufferRows_hpp
#define SimdThreeCenterElectronRepulsionGeom100BufferRows_hpp

#include <array>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

/// @brief Gets the number of rows of the buffer a combination of basis functions
/// needs.
/// @param a_angular_momentum The angular momentum of basis function on a side.
/// @param b_angular_momentum The angular momentum of basis function on b side.
/// @param c_angular_momentum The angular momentum of basis function on c side.
/// @return The number of rows.
/// @note The buffer belongs to the block and not to the combination, so a caller
/// sizes it once from the highest angular momenta the block carries. The numbers
/// below do not decrease with either angular momentum, so the largest combination of
/// a block is the one of its highest momenta and no other combination needs more.
/// @note This table is written with the kernels and describes them. Regenerating the
/// kernels without regenerating it leaves a buffer which is too small, which the
/// assertions of simdfunc::prepare_buffer report rather than let pass.
/// @note A combination with momentum on both functions of the bra has no kernel
/// written for it yet. Its entry is what the largest combination below it needs,
/// so that an arena sized from it holds every combination of the block that does
/// have one.

// NOTE: the name is qualified by which kernels the table describes. Every set in
// this namespace declares its own table as an inline function, so a shared name is
// one definition of several different things: the program is ill formed, the linker
// keeps whichever it saw first, and the tables differ both in their numbers and in
// their shape. A gradient set given the plain set's table would size its buffer
// several times too small and overrun it. The generator still emits the bare name;
// until it does not, a regenerated set has to be renamed again.
inline auto
number_of_geom_100_buffer_rows(const int a_angular_momentum, const int b_angular_momentum,
                      const int c_angular_momentum) -> size_t
{
    constexpr std::array<std::array<std::array<size_t, 7>, 7>, 7> rows{{
        {{
            {      18,      39,      95,     199,     354,     605,     957},
            {      46,     121,     399,     801,    1493,    2515,    4035},
            {     107,     311,    1030,    2114,    3928,    6654,   10650},
            {     218,     656,    2174,    4502,    8354,   14162,   22614},
            {     400,    1222,    4049,    8417,   15594,   26407,   42059},
            {     400,    1222,    4049,    8417,   15594,   26407,   42059},
            {     400,    1222,    4049,    8417,   15594,   26407,   42059}
        }},
        {{
            {      40,     106,     300,     588,    1049,    1738,    2730},
            {     165,     471,    1069,    1899,    3121,    4850,    7241},
            {     425,    1226,    2640,    4492,    7097,   10665,   15476},
            {     966,    2800,    5767,    9475,   14484,   21144,   29917},
            {    2001,    5821,   11559,   18459,   27445,   39063,   54027},
            {    2001,    5821,   11559,   18459,   27445,   39063,   54027},
            {    2001,    5821,   11559,   18459,   27445,   39063,   54027}
        }},
        {{
            {      78,     207,     643,    1269,    2266,    3749,    5873},
            {     305,     863,    2019,    3589,    5891,    9135,   13601},
            {     744,    2131,    4636,    7858,   12360,   18492,   26716},
            {    1620,    4675,    9633,   15729,   23890,   34662,   48759},
            {    3240,    9403,   18594,   29484,   43516,   61500,   84486},
            {    3240,    9403,   18594,   29484,   43516,   61500,   84486},
            {    3240,    9403,   18594,   29484,   43516,   61500,   84486}
        }},
        {{
            {     137,     350,    1187,    2357,    4217,    6977,   10917},
            {     498,    1384,    3364,    6010,    9894,   15366,   22888},
            {    1166,    3304,    7321,   12425,   19552,   29248,   42227},
            {    2464,    7066,   14672,   23926,   36280,   52544,   73768},
            {    4812,   13915,   27571,   43599,   64156,   90397,  123807},
            {    4812,   13915,   27571,   43599,   64156,   90397,  123807},
            {    4812,   13915,   27571,   43599,   64156,   90397,  123807}
        }},
        {{
            {     223,     541,    1994,    3980,    7134,   11806,   18458},
            {     751,    2041,    5176,    9310,   15397,   23983,   35782},
            {    1699,    4753,   10777,   18361,   28975,   43429,   62773},
            {    3507,    9982,   20976,   34254,   51991,   75342,  105792},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922}
        }},
        {{
            {     223,     541,    1994,    3980,    7134,   11806,   18458},
            {     751,    2041,    5176,    9310,   15397,   23983,   35782},
            {    1699,    4753,   10777,   18361,   28975,   43429,   62773},
            {    3507,    9982,   20976,   34254,   51991,   75342,  105792},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922}
        }},
        {{
            {     223,     541,    1994,    3980,    7134,   11806,   18458},
            {     751,    2041,    5176,    9310,   15397,   23983,   35782},
            {    1699,    4753,   10777,   18361,   28975,   43429,   62773},
            {    3507,    9982,   20976,   34254,   51991,   75342,  105792},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922},
            {    6727,   19367,   38592,   61012,   89737,  126362,  172922}
        }}
    }};

    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 7) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum < 7) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum < 7),
                              std::string("SimdThreeCenterElectronRepulsionGeom100BufferRows.number_of_geom_100_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(a_angular_momentum)][static_cast<size_t>(b_angular_momentum)][static_cast<size_t>(c_angular_momentum)];
}

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionGeom100BufferRows_hpp */
