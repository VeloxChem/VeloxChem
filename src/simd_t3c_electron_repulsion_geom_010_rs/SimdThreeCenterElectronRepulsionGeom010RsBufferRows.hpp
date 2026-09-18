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


#ifndef SimdThreeCenterElectronRepulsionGeom010RsBufferRows_hpp
#define SimdThreeCenterElectronRepulsionGeom010RsBufferRows_hpp

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
// their shape. The generator still emits the bare name; until it does not, a
// regenerated set has to be renamed again.
inline auto
number_of_geom_010_rs_buffer_rows(const int a_angular_momentum, const int b_angular_momentum,
                      const int c_angular_momentum) -> size_t
{
    constexpr std::array<std::array<std::array<size_t, 7>, 7>, 7> rows{{
        {{
            {      30,      72,     184,     344,     702,    1156,    1908},
            {      74,     197,     579,    1113,    2065,    3401,    5415},
            {     150,     390,    1250,    2454,    4472,    7390,   11662},
            {     268,     664,    2318,    4602,    8338,   13802,   21698},
            {     440,    1031,    3907,    7813,   14127,   23405,   36715},
            {     440,    1031,    3907,    7813,   14127,   23405,   36715},
            {     440,    1031,    3907,    7813,   14127,   23405,   36715}
        }},
        {{
            {      83,     224,     774,    1572,    2950,    4988,    8022},
            {     315,     909,    2087,    3693,    6155,    9559,   14359},
            {     583,    1657,    3927,    6989,   11587,   17997,   26923},
            {     951,    2645,    6527,   11705,   19431,   30261,   45263},
            {    1433,    3887,   10031,   18137,   30221,   47231,   70739},
            {    1433,    3887,   10031,   18137,   30221,   47231,   70739},
            {    1433,    3887,   10031,   18137,   30221,   47231,   70739}
        }},
        {{
            {     223,     601,    2033,    4177,    7805,   13233,   21225},
            {     832,    2410,    5214,    8858,   14080,   21156,   30790},
            {    1452,    4166,    9116,   15464,   24444,   36612,   53036},
            {    2268,    6428,   14346,   24402,   38576,   57816,   83694},
            {    3296,    9212,   21068,   36008,   57080,   85760,  124292},
            {    3296,    9212,   21068,   36008,   57080,   85760,  124292},
            {    3296,    9212,   21068,   36008,   57080,   85760,  124292}
        }},
        {{
            {     451,    1279,    4301,    8925,   16621,   28205,   45101},
            {    1914,    5558,   11468,   18824,   28854,   42114,   59672},
            {    3200,    9242,   19090,   31178,   47468,   68908,   97070},
            {    4852,   13916,   28988,   47320,   71924,  104276,  146620},
            {    6888,   19598,   41346,   67626,  102896,  149322,  210018},
            {    6888,   19598,   41346,   67626,  102896,  149322,  210018},
            {    6888,   19598,   41346,   67626,  102896,  149322,  210018}
        }},
        {{
            {     821,    2396,    8026,   16720,   31056,   52640,   83926},
            {    3987,   11609,   23067,   36813,   54803,   77985,  107931},
            {    6441,   18701,   37017,   58695,   86729,  122595,  168537},
            {    9543,   27599,   54761,   86631,  127631,  179927,  246633},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083}
        }},
        {{
            {     821,    2396,    8026,   16720,   31056,   52640,   83926},
            {    3987,   11609,   23067,   36813,   54803,   77985,  107931},
            {    6441,   18701,   37017,   58695,   86729,  122595,  168537},
            {    9543,   27599,   54761,   86631,  127631,  179927,  246633},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083}
        }},
        {{
            {     821,    2396,    8026,   16720,   31056,   52640,   83926},
            {    3987,   11609,   23067,   36813,   54803,   77985,  107931},
            {    6441,   18701,   37017,   58695,   86729,  122595,  168537},
            {    9543,   27599,   54761,   86631,  127631,  179927,  246633},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083},
            {   13313,   38323,   76503,  121037,  178253,  251197,  344083}
        }}
    }};

    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 7) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum < 7) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum < 7),
                              std::string("SimdThreeCenterElectronRepulsionGeom010RsBufferRows.number_of_geom_010_rs_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(a_angular_momentum)][static_cast<size_t>(b_angular_momentum)][static_cast<size_t>(c_angular_momentum)];
}

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionGeom010RsBufferRows_hpp */
