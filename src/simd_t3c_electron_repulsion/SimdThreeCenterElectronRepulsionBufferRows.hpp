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


#ifndef SimdThreeCenterElectronRepulsionBufferRows_hpp
#define SimdThreeCenterElectronRepulsionBufferRows_hpp

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
inline auto
number_of_buffer_rows(const int a_angular_momentum, const int b_angular_momentum,
                      const int c_angular_momentum) -> size_t
{
    constexpr std::array<std::array<std::array<size_t, 9>, 9>, 9> rows{{
        {{
            {       5,      11,      22,      42,      76,     129,     208,     320,     474},
            {      15,      39,      92,     166,     336,     551,     912,    1376,    2052},
            {      28,      76,     243,     471,     896,    1501,    2433,    3695,    5454},
            {      52,     132,     507,    1019,    1912,    3229,    5197,    7899,   11602},
            {      92,     212,     942,    1924,    3587,    6069,    9722,   14754,   21587},
            {     154,     322,    1610,    3314,    6153,   10405,   16604,   25140,   36659},
            {     245,     469,    2583,    5337,    9877,   16677,   26523,   40057,   58233},
            {     245,     469,    2583,    5337,    9877,   16677,   26523,   40057,   58233},
            {     245,     469,    2583,    5337,    9877,   16677,   26523,   40057,   58233}
        }},
        {{
            {      15,      39,      92,     190,     336,     575,     912,    1400,    2052},
            {      58,     157,     366,     666,    1130,    1813,    2790,    4136,    5946},
            {     107,     285,     738,    1352,    2353,    3784,    5872,    8700,   12535},
            {     178,     455,    1307,    2451,    4286,    6950,   10795,   16029,   23074},
            {     277,     673,    2135,    4073,    7161,   11677,   18155,   26985,   38813},
            {     411,     946,    3294,    6366,   11245,   18405,   28632,   42568,   61167},
            {     588,    1282,    4866,    9498,   16840,   27630,   42990,   63898,   91716},
            {     588,    1282,    4866,    9498,   16840,   27630,   42990,   63898,   91716},
            {     588,    1282,    4866,    9498,   16840,   27630,   42990,   63898,   91716}
        }},
        {{
            {      28,      76,     243,     489,     896,    1519,    2433,    3713,    5454},
            {     110,     294,     753,    1391,    2380,    3835,    5911,    8763,   12586},
            {     253,     698,    1688,    3012,    4997,    7853,   11860,   17298,   24517},
            {     402,    1078,    2770,    4978,    8346,   13152,   19930,   29070,   41218},
            {     597,    1549,    4239,    7713,   13009,   20601,   31275,   45673,   64749},
            {     846,    2119,    6177,   11367,   19288,   30678,   46659,   68209,   96690},
            {    1158,    2797,    8676,   16128,   27520,   43935,   66930,   97918,  138786},
            {    1158,    2797,    8676,   16128,   27520,   43935,   66930,   97918,  138786},
            {    1158,    2797,    8676,   16128,   27520,   43935,   66930,   97918,  138786}
        }},
        {{
            {      52,     132,     507,    1037,    1912,    3247,    5197,    7917,   11602},
            {     187,     482,    1352,    2532,    4367,    7067,   10912,   16182,   23227},
            {     410,    1102,    2810,    5052,    8418,   13258,   20034,   29208,   41354},
            {     823,    2257,    5379,    9343,   15103,   23205,   34363,   49291,   68871},
            {    1185,    3181,    7887,   13785,   22429,   34557,   51291,   73609,  102873},
            {    1631,    4279,   11062,   19502,   31903,   49348,   73394,  105454,  147415},
            {    2171,    5561,   15006,   26684,   43897,   68168,  101604,  146168,  204407},
            {    2171,    5561,   15006,   26684,   43897,   68168,  101604,  146168,  204407},
            {    2171,    5561,   15006,   26684,   43897,   68168,  101604,  146168,  204407}
        }},
        {{
            {      92,     212,     942,    1942,    3587,    6087,    9722,   14772,   21587},
            {     295,     727,    2225,    4217,    7323,   11893,   18389,   27273,   39119},
            {     618,    1612,    4344,    7878,   13198,   20850,   31548,   46006,   65106},
            {    1200,    3226,    7962,   13908,   22564,   34740,   51486,   73852,  103128},
            {    2196,    6019,   13902,   23520,   37078,   55731,   80964,  114262,  157440},
            {    2962,    7997,   18982,   32284,   51142,   77079,  112202,  158474,  218442},
            {    3869,   10281,   25124,   43014,   68463,  103544,  151046,  213614,  294609},
            {    3869,   10281,   25124,   43014,   68463,  103544,  151046,  213614,  294609},
            {    3869,   10281,   25124,   43014,   68463,  103544,  151046,  213614,  294609}
        }},
        {{
            {     154,     322,    1610,    3332,    6153,   10423,   16604,   25158,   36659},
            {     441,    1036,    3444,    6594,   11515,   18753,   29022,   43036,   61677},
            {     885,    2236,    6372,   11658,   19639,   31125,   47166,   68812,   97353},
            {    1668,    4390,   11247,   19779,   32236,   49773,   73875,  106027,  148044},
            {    2986,    8069,   19102,   32470,   51358,   77361,  112514,  158852,  218850},
            {    5102,   14043,   31310,   51726,   79686,  117335,  167390,  232568,  316158},
            {    6573,   17899,   40698,   67532,  104447,  154192,  220388,  306512,  416913},
            {    6573,   17899,   40698,   67532,  104447,  154192,  220388,  306512,  416913},
            {    6573,   17899,   40698,   67532,  104447,  154192,  220388,  306512,  416913}
        }},
        {{
            {     245,     469,    2583,    5355,    9877,   16695,   26523,   40075,   58233},
            {     633,    1417,    5091,    9831,   17245,   28143,   43575,   64591,   92481},
            {    1220,    2983,    8986,   16580,   28078,   44635,   67736,   98866,  139840},
            {    2237,    5759,   15336,   27164,   44491,   68912,  102462,  147176,  205529},
            {    3926,   10452,   25409,   43431,   68976,  104189,  151787,  214487,  295578},
            {    6608,   18004,   40873,   67795,  104762,  154595,  220843,  307055,  417508},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552}
        }},
        {{
            {     245,     469,    2583,    5355,    9877,   16695,   26523,   40075,   58233},
            {     633,    1417,    5091,    9831,   17245,   28143,   43575,   64591,   92481},
            {    1220,    2983,    8986,   16580,   28078,   44635,   67736,   98866,  139840},
            {    2237,    5759,   15336,   27164,   44491,   68912,  102462,  147176,  205529},
            {    3926,   10452,   25409,   43431,   68976,  104189,  151787,  214487,  295578},
            {    6608,   18004,   40873,   67795,  104762,  154595,  220843,  307055,  417508},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552}
        }},
        {{
            {     245,     469,    2583,    5355,    9877,   16695,   26523,   40075,   58233},
            {     633,    1417,    5091,    9831,   17245,   28143,   43575,   64591,   92481},
            {    1220,    2983,    8986,   16580,   28078,   44635,   67736,   98866,  139840},
            {    2237,    5759,   15336,   27164,   44491,   68912,  102462,  147176,  205529},
            {    3926,   10452,   25409,   43431,   68976,  104189,  151787,  214487,  295578},
            {    6608,   18004,   40873,   67795,  104762,  154595,  220843,  307055,  417508},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552},
            {   10696,   29631,   63861,  103327,  155932,  225316,  316029,  432621,  580552}
        }}
    }};

    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 9) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum < 9) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum < 9),
                              std::string("SimdThreeCenterElectronRepulsionBufferRows.number_of_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(a_angular_momentum)][static_cast<size_t>(b_angular_momentum)][static_cast<size_t>(c_angular_momentum)];
}

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionBufferRows_hpp */
