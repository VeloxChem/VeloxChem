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


#ifndef SimdThreeCenterElectronRepulsionGeom010BufferRows_hpp
#define SimdThreeCenterElectronRepulsionGeom010BufferRows_hpp

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
number_of_geom_010_buffer_rows(const int a_angular_momentum, const int b_angular_momentum,
                      const int c_angular_momentum) -> size_t
{
    constexpr std::array<std::array<std::array<size_t, 7>, 7>, 7> rows{{
        {{
            {      18,      39,      95,     175,     354,     581,     957},
            {      40,     106,     300,     570,    1049,    1720,    2730},
            {      78,     207,     643,    1251,    2266,    3731,    5873},
            {     137,     350,    1187,    2339,    4217,    6959,   10917},
            {     223,     541,    1994,    3962,    7134,   11788,   18458},
            {     223,     541,    1994,    3962,    7134,   11788,   18458},
            {     223,     541,    1994,    3962,    7134,   11788,   18458}
        }},
        {{
            {      46,     121,     399,     801,    1493,    2515,    4035},
            {     165,     471,    1069,    1881,    3121,    4832,    7241},
            {     302,     854,    2004,    3550,    5864,    9084,   13562},
            {     489,    1357,    3319,    5929,    9813,   15249,   22771},
            {     733,    1987,    5086,    9166,   15235,   23767,   35548},
            {     733,    1987,    5086,    9166,   15235,   23767,   35548},
            {     733,    1987,    5086,    9166,   15235,   23767,   35548}
        }},
        {{
            {     116,     314,    1036,    2114,    3934,    6654,   10656},
            {     428,    1235,    2655,    4495,    7124,   10680,   15515},
            {     744,    2131,    4636,    7840,   12360,   18474,   26716},
            {    1158,    3280,    7281,   12351,   19480,   29142,   42123},
            {    1678,    4690,   10672,   18196,   28786,   43180,   62500},
            {    1678,    4690,   10672,   18196,   28786,   43180,   62500},
            {    1678,    4690,   10672,   18196,   28786,   43180,   62500}
        }},
        {{
            {     230,     659,    2180,    4502,    8360,   14162,   22620},
            {     975,    2827,    5812,    9520,   14565,   21225,   30034},
            {    1628,    4699,    9673,   15767,   23962,   34732,   48863},
            {    2464,    7066,   14672,   23908,   36280,   52526,   73768},
            {    3492,    9937,   20901,   34131,   51856,   75159,  105597},
            {    3492,    9937,   20901,   34131,   51856,   75159,  105597},
            {    3492,    9937,   20901,   34131,   51856,   75159,  105597}
        }},
        {{
            {     415,    1225,    4055,    8417,   15600,   26407,   42065},
            {    2019,    5875,   11649,   18567,   27607,   39243,   54261},
            {    3261,    9466,   18699,   29613,   43705,   61713,   84759},
            {    4827,   13960,   27646,   43686,   64291,   90544,  124002},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922}
        }},
        {{
            {     415,    1225,    4055,    8417,   15600,   26407,   42065},
            {    2019,    5875,   11649,   18567,   27607,   39243,   54261},
            {    3261,    9466,   18699,   29613,   43705,   61713,   84759},
            {    4827,   13960,   27646,   43686,   64291,   90544,  124002},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922}
        }},
        {{
            {     415,    1225,    4055,    8417,   15600,   26407,   42065},
            {    2019,    5875,   11649,   18567,   27607,   39243,   54261},
            {    3261,    9466,   18699,   29613,   43705,   61713,   84759},
            {    4827,   13960,   27646,   43686,   64291,   90544,  124002},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922},
            {    6727,   19367,   38592,   60994,   89737,  126344,  172922}
        }}
    }};

    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 7) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum < 7) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum < 7),
                              std::string("SimdThreeCenterElectronRepulsionGeom010BufferRows.number_of_geom_010_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(a_angular_momentum)][static_cast<size_t>(b_angular_momentum)][static_cast<size_t>(c_angular_momentum)];
}

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionGeom010BufferRows_hpp */
