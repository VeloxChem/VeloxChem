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


#ifndef SimdThreeCenterElectronRepulsionGeom100RsBufferRows_hpp
#define SimdThreeCenterElectronRepulsionGeom100RsBufferRows_hpp

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
number_of_geom_100_rs_buffer_rows(const int a_angular_momentum, const int b_angular_momentum,
                      const int c_angular_momentum) -> size_t
{
    constexpr std::array<std::array<std::array<size_t, 7>, 7>, 7> rows{{
        {{
            {      30,      72,     184,     392,     702,    1204,    1908},
            {      83,     224,     774,    1572,    2950,    4988,    8022},
            {     205,     595,    2021,    4177,    7793,   13233,   21213},
            {     427,    1273,    4289,    8925,   16609,   28205,   45089},
            {     791,    2390,    8014,   16720,   31044,   52640,   83914},
            {     791,    2390,    8014,   16720,   31044,   52640,   83914},
            {     791,    2390,    8014,   16720,   31044,   52640,   83914}
        }},
        {{
            {      74,     197,     579,    1149,    2065,    3437,    5415},
            {     315,     909,    2087,    3729,    6155,    9595,   14359},
            {     829,    2401,    5199,    8873,   14053,   21159,   30751},
            {    1905,    5531,   11423,   18797,   28773,   42051,   59555},
            {    3969,   11555,   22977,   36723,   54641,   77823,  107697},
            {    3969,   11555,   22977,   36723,   54641,   77823,  107697},
            {    3969,   11555,   22977,   36723,   54641,   77823,  107697}
        }},
        {{
            {     150,     390,    1250,    2490,    4472,    7426,   11662},
            {     586,    1666,    3942,    7046,   11614,   18066,   26962},
            {    1452,    4166,    9116,   15500,   24444,   36648,   53036},
            {    3192,    9218,   19050,   31158,   47396,   68856,   96966},
            {    6420,   18638,   36912,   58584,   86540,  122400,  168264},
            {    6420,   18638,   36912,   58584,   86540,  122400,  168264},
            {    6420,   18638,   36912,   58584,   86540,  122400,  168264}
        }},
        {{
            {     268,     664,    2318,    4638,    8338,   13838,   21698},
            {     960,    2672,    6572,   11804,   19512,   30396,   45380},
            {    2276,    6452,   14386,   24494,   38648,   57940,   83798},
            {    4852,   13916,   28988,   47356,   71924,  104312,  146620},
            {    9528,   27554,   54686,   86562,  127496,  179798,  246438},
            {    9528,   27554,   54686,   86562,  127496,  179798,  246438},
            {    9528,   27554,   54686,   86562,  127496,  179798,  246438}
        }},
        {{
            {     440,    1031,    3907,    7849,   14127,   23441,   36715},
            {    1451,    3941,   10121,   18299,   30383,   47465,   70973},
            {    3317,    9275,   21173,   36191,   57269,   86027,  124565},
            {    6903,   19643,   41421,   67767,  103031,  149523,  210213},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083}
        }},
        {{
            {     440,    1031,    3907,    7849,   14127,   23441,   36715},
            {    1451,    3941,   10121,   18299,   30383,   47465,   70973},
            {    3317,    9275,   21173,   36191,   57269,   86027,  124565},
            {    6903,   19643,   41421,   67767,  103031,  149523,  210213},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083}
        }},
        {{
            {     440,    1031,    3907,    7849,   14127,   23441,   36715},
            {    1451,    3941,   10121,   18299,   30383,   47465,   70973},
            {    3317,    9275,   21173,   36191,   57269,   86027,  124565},
            {    6903,   19643,   41421,   67767,  103031,  149523,  210213},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083},
            {   13313,   38323,   76503,  121073,  178253,  251233,  344083}
        }}
    }};

    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 7) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum < 7) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum < 7),
                              std::string("SimdThreeCenterElectronRepulsionGeom100RsBufferRows.number_of_geom_100_rs_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(a_angular_momentum)][static_cast<size_t>(b_angular_momentum)][static_cast<size_t>(c_angular_momentum)];
}

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionGeom100RsBufferRows_hpp */
