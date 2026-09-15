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
            {       6,      12,      23,      43,      77,     130,     209,     321,     475},
            {      16,      40,      93,     167,     337,     552,     913,    1377,    2053},
            {      29,      77,     244,     472,     897,    1502,    2434,    3696,    5455},
            {      53,     133,     508,    1020,    1913,    3230,    5198,    7900,   11603},
            {      93,     213,     943,    1925,    3588,    6070,    9723,   14755,   21588},
            {     155,     323,    1611,    3315,    6154,   10406,   16605,   25141,   36660},
            {     246,     470,    2584,    5338,    9878,   16678,   26524,   40058,   58234},
            {     246,     470,    2584,    5338,    9878,   16678,   26524,   40058,   58234},
            {     246,     470,    2584,    5338,    9878,   16678,   26524,   40058,   58234}
        }},
        {{
            {      16,      40,      93,     191,     337,     576,     913,    1401,    2053},
            {      59,     158,     367,     667,    1131,    1814,    2791,    4137,    5947},
            {     108,     286,     739,    1353,    2354,    3785,    5873,    8701,   12536},
            {     179,     456,    1308,    2452,    4287,    6951,   10796,   16030,   23075},
            {     278,     674,    2136,    4074,    7162,   11678,   18156,   26986,   38814},
            {     412,     947,    3295,    6367,   11246,   18406,   28633,   42569,   61168},
            {     589,    1283,    4867,    9499,   16841,   27631,   42991,   63899,   91717},
            {     589,    1283,    4867,    9499,   16841,   27631,   42991,   63899,   91717},
            {     589,    1283,    4867,    9499,   16841,   27631,   42991,   63899,   91717}
        }},
        {{
            {      29,      77,     244,     490,     897,    1520,    2434,    3714,    5455},
            {     111,     295,     754,    1392,    2381,    3836,    5912,    8764,   12587},
            {     254,     699,    1689,    3013,    4998,    7854,   11861,   17299,   24518},
            {     403,    1079,    2771,    4979,    8347,   13153,   19931,   29071,   41219},
            {     598,    1550,    4240,    7714,   13010,   20602,   31276,   45674,   64750},
            {     847,    2120,    6178,   11368,   19289,   30679,   46660,   68210,   96691},
            {    1159,    2798,    8677,   16129,   27521,   43936,   66931,   97919,  138787},
            {    1159,    2798,    8677,   16129,   27521,   43936,   66931,   97919,  138787},
            {    1159,    2798,    8677,   16129,   27521,   43936,   66931,   97919,  138787}
        }},
        {{
            {      53,     133,     508,    1038,    1913,    3248,    5198,    7918,   11603},
            {     188,     483,    1353,    2533,    4368,    7068,   10913,   16183,   23228},
            {     411,    1103,    2811,    5053,    8419,   13259,   20035,   29209,   41355},
            {     824,    2258,    5380,    9344,   15104,   23206,   34364,   49292,   68872},
            {    1186,    3182,    7888,   13786,   22430,   34558,   51292,   73610,  102874},
            {    1632,    4280,   11063,   19503,   31904,   49349,   73395,  105455,  147416},
            {    2172,    5562,   15007,   26685,   43898,   68169,  101605,  146169,  204408},
            {    2172,    5562,   15007,   26685,   43898,   68169,  101605,  146169,  204408},
            {    2172,    5562,   15007,   26685,   43898,   68169,  101605,  146169,  204408}
        }},
        {{
            {      93,     213,     943,    1943,    3588,    6088,    9723,   14773,   21588},
            {     296,     728,    2226,    4218,    7324,   11894,   18390,   27274,   39120},
            {     619,    1613,    4345,    7879,   13199,   20851,   31549,   46007,   65107},
            {    1201,    3227,    7963,   13909,   22565,   34741,   51487,   73853,  103129},
            {    2197,    6020,   13903,   23521,   37079,   55732,   80965,  114263,  157441},
            {    2963,    7998,   18983,   32285,   51143,   77080,  112203,  158475,  218443},
            {    3870,   10282,   25125,   43015,   68464,  103545,  151047,  213615,  294610},
            {    3870,   10282,   25125,   43015,   68464,  103545,  151047,  213615,  294610},
            {    3870,   10282,   25125,   43015,   68464,  103545,  151047,  213615,  294610}
        }},
        {{
            {     155,     323,    1611,    3333,    6154,   10424,   16605,   25159,   36660},
            {     442,    1037,    3445,    6595,   11516,   18754,   29023,   43037,   61678},
            {     886,    2237,    6373,   11659,   19640,   31126,   47167,   68813,   97354},
            {    1669,    4391,   11248,   19780,   32237,   49774,   73876,  106028,  148045},
            {    2987,    8070,   19103,   32471,   51359,   77362,  112515,  158853,  218851},
            {    5103,   14044,   31311,   51727,   79687,  117336,  167391,  232569,  316159},
            {    6574,   17900,   40699,   67533,  104448,  154193,  220389,  306513,  416914},
            {    6574,   17900,   40699,   67533,  104448,  154193,  220389,  306513,  416914},
            {    6574,   17900,   40699,   67533,  104448,  154193,  220389,  306513,  416914}
        }},
        {{
            {     246,     470,    2584,    5356,    9878,   16696,   26524,   40076,   58234},
            {     634,    1418,    5092,    9832,   17246,   28144,   43576,   64592,   92482},
            {    1221,    2984,    8987,   16581,   28079,   44636,   67737,   98867,  139841},
            {    2238,    5760,   15337,   27165,   44492,   68913,  102463,  147177,  205530},
            {    3927,   10453,   25410,   43432,   68977,  104190,  151788,  214488,  295579},
            {    6609,   18005,   40874,   67796,  104763,  154596,  220844,  307056,  417509},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553}
        }},
        {{
            {     246,     470,    2584,    5356,    9878,   16696,   26524,   40076,   58234},
            {     634,    1418,    5092,    9832,   17246,   28144,   43576,   64592,   92482},
            {    1221,    2984,    8987,   16581,   28079,   44636,   67737,   98867,  139841},
            {    2238,    5760,   15337,   27165,   44492,   68913,  102463,  147177,  205530},
            {    3927,   10453,   25410,   43432,   68977,  104190,  151788,  214488,  295579},
            {    6609,   18005,   40874,   67796,  104763,  154596,  220844,  307056,  417509},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553}
        }},
        {{
            {     246,     470,    2584,    5356,    9878,   16696,   26524,   40076,   58234},
            {     634,    1418,    5092,    9832,   17246,   28144,   43576,   64592,   92482},
            {    1221,    2984,    8987,   16581,   28079,   44636,   67737,   98867,  139841},
            {    2238,    5760,   15337,   27165,   44492,   68913,  102463,  147177,  205530},
            {    3927,   10453,   25410,   43432,   68977,  104190,  151788,  214488,  295579},
            {    6609,   18005,   40874,   67796,  104763,  154596,  220844,  307056,  417509},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553},
            {   10697,   29632,   63862,  103328,  155933,  225317,  316030,  432622,  580553}
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
