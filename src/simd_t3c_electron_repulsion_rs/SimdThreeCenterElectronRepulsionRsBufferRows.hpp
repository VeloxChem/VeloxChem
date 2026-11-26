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


#ifndef SimdThreeCenterElectronRepulsionRsBufferRows_hpp
#define SimdThreeCenterElectronRepulsionRsBufferRows_hpp

#include <array>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

// NOTE: the name differs from the table of the unattenuated kernels on purpose, as
// it does in the two-center range separated set. Both tables live in this namespace
// and both are inline, so a shared name is one definition of two different things:
// a translation unit which includes both does not compile, and one which includes
// either takes whichever the linker kept. The range separated buffer is the larger
// of the two, so taking the wrong table is a buffer overrun and not a wrong answer.
// This wants fixing in the generator; it has now been renamed by hand twice.

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
number_of_rs_buffer_rows(const int a_angular_momentum, const int b_angular_momentum,
                      const int c_angular_momentum) -> size_t
{
    constexpr std::array<std::array<std::array<size_t, 9>, 9>, 9> rows{{
        {{
            {       7,      19,      41,      81,     149,     255,     413,     637,     945},
            {      24,      63,     163,     305,     639,    1063,    1779,    2701,    4047},
            {      50,     128,     450,     894,    1732,    2930,    4782,    7294,   10800},
            {      98,     228,     958,    1962,    3728,    6342,   10258,   15642,   23028},
            {     178,     373,    1803,    3737,    7033,   11967,   19243,   29277,   42913},
            {     302,     575,    3109,    6475,   12111,   20573,   32929,   49959,   72955},
            {     484,     848,    5020,   10472,   19496,   33040,   52676,   79688,  115984},
            {     484,     848,    5020,   10472,   19496,   33040,   52676,   79688,  115984},
            {     484,     848,    5020,   10472,   19496,   33040,   52676,   79688,  115984}
        }},
        {{
            {      24,      63,     163,     353,     639,    1111,    1779,    2749,    4047},
            {     101,     281,     681,    1263,    2173,    3521,    5457,    8131,   11733},
            {     193,     519,    1395,    2593,    4565,    7397,   11543,   17169,   24809},
            {     329,     841,    2503,    4749,    8377,   13663,   21311,   31737,   45785},
            {     521,    1259,    4129,    7951,   14073,   23051,   35953,   53559,   77161},
            {     783,    1787,    6417,   12495,   22187,   36441,   56829,   84635,  121767},
            {    1131,    2441,    9531,   18717,   33323,   54825,   85467,  127205,  182763},
            {    1131,    2441,    9531,   18717,   33323,   54825,   85467,  127205,  182763},
            {    1131,    2441,    9531,   18717,   33323,   54825,   85467,  127205,  182763}
        }},
        {{
            {      50,     128,     450,     930,    1732,    2966,    4782,    7330,   10800},
            {     196,     528,    1410,    2650,    4592,    7466,   11582,   17250,   24860},
            {     470,    1300,    3220,    5808,    9718,   15370,   23324,   34140,   48518},
            {     756,    2024,    5324,    9656,   16308,   25836,   39308,   57504,   81716},
            {    1134,    2930,    8202,   15042,   25526,   40602,   61842,   90530,  128574},
            {    1620,    4034,   12018,   22266,   37976,   60624,   92454,  135422,  192252},
            {    2232,    5354,   16956,   31704,   54332,   87006,  132840,  194660,  276240},
            {    2232,    5354,   16956,   31704,   54332,   87006,  132840,  194660,  276240},
            {    2232,    5354,   16956,   31704,   54332,   87006,  132840,  194660,  276240}
        }},
        {{
            {      98,     228,     958,    1998,    3728,    6378,   10258,   15678,   23028},
            {     338,     868,    2548,    4848,    8458,   13798,   21428,   31908,   45938},
            {     764,    2048,    5364,    9748,   16380,   25960,   39412,   57660,   81852},
            {    1570,    4298,   10402,   18190,   29570,   45634,   67810,   97526,  136546},
            {    2274,    6086,   15318,   26934,   44042,   68118,  101406,  145862,  204210},
            {    3146,    8222,   21568,   38228,   62810,   97480,  145352,  209252,  292954},
            {    4206,   10726,   29356,   52452,   86618,  134900,  201512,  290380,  406598},
            {    4206,   10726,   29356,   52452,   86618,  134900,  201512,  290380,  406598},
            {    4206,   10726,   29356,   52452,   86618,  134900,  201512,  290380,  406598}
        }},
        {{
            {     178,     373,    1803,    3773,    7033,   12003,   19243,   29313,   42913},
            {     539,    1313,    4219,    8113,   14235,   23285,   36187,   53865,   77467},
            {    1155,    2993,    8307,   15225,   25715,   40869,   62115,   90881,  128931},
            {    2289,    6131,   15393,   27075,   44177,   68319,  101601,  146123,  204465},
            {    4251,   11627,   27123,   46089,   72935,  109971,  160167,  226493,  312579},
            {    5753,   15493,   37133,   63407,  100793,  152337,  222253,  314467,  434073},
            {    7537,   19971,   49267,   84657,  135165,  204937,  299551,  424297,  585897},
            {    7537,   19971,   49267,   84657,  135165,  204937,  299551,  424297,  585897},
            {    7537,   19971,   49267,   84657,  135165,  204937,  299551,  424297,  585897}
        }},
        {{
            {     302,     575,    3109,    6511,   12111,   20609,   32929,   49995,   72955},
            {     813,    1877,    6567,   12741,   22457,   36807,   57219,   85121,  122277},
            {    1659,    4151,   12213,   22575,   38327,   61089,   92961,  136043,  192915},
            {    3183,    8333,   21753,   38523,   63143,   97923,  145833,  209843,  293583},
            {    5777,   15565,   37253,   63611,  101009,  152637,  222565,  314863,  434481},
            {    9967,   27387,   61459,  101829,  157287,  232123,  331771,  461665,  628383},
            {   12867,   34973,   80025,  133147,  206431,  305375,  437221,  608923,  829179},
            {   12867,   34973,   80025,  133147,  206431,  305375,  437221,  608923,  829179},
            {   12867,   34973,   80025,  133147,  206431,  305375,  437221,  608923,  829179}
        }},
        {{
            {     484,     848,    5020,   10508,   19496,   33076,   52676,   79724,  115984},
            {    1176,    2576,    9756,   19068,   33728,   55356,   86052,  127916,  183528},
            {    2294,    5540,   17266,   32174,   54890,   87724,  133646,  195626,  277294},
            {    4272,   10924,   29686,   52950,   87212,  135662,  202370,  291406,  407720},
            {    7594,   20142,   49552,   85092,  135678,  205600,  300292,  425188,  586866},
            {   12902,   35078,   80200,  133428,  206746,  305796,  437676,  609484,  829774},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910}
        }},
        {{
            {     484,     848,    5020,   10508,   19496,   33076,   52676,   79724,  115984},
            {    1176,    2576,    9756,   19068,   33728,   55356,   86052,  127916,  183528},
            {    2294,    5540,   17266,   32174,   54890,   87724,  133646,  195626,  277294},
            {    4272,   10924,   29686,   52950,   87212,  135662,  202370,  291406,  407720},
            {    7594,   20142,   49552,   85092,  135678,  205600,  300292,  425188,  586866},
            {   12902,   35078,   80200,  133428,  206746,  305796,  437676,  609484,  829774},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910}
        }},
        {{
            {     484,     848,    5020,   10508,   19496,   33076,   52676,   79724,  115984},
            {    1176,    2576,    9756,   19068,   33728,   55356,   86052,  127916,  183528},
            {    2294,    5540,   17266,   32174,   54890,   87724,  133646,  195626,  277294},
            {    4272,   10924,   29686,   52950,   87212,  135662,  202370,  291406,  407720},
            {    7594,   20142,   49552,   85092,  135678,  205600,  300292,  425188,  586866},
            {   12902,   35078,   80200,  133428,  206746,  305796,  437676,  609484,  829774},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910},
            {   21022,   58164,  125896,  204100,  308582,  446622,  627320,  859776, 1154910}
        }}
    }};

    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 9) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum < 9) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum < 9),
                              std::string("SimdThreeCenterElectronRepulsionRsBufferRows.number_of_rs_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(a_angular_momentum)][static_cast<size_t>(b_angular_momentum)][static_cast<size_t>(c_angular_momentum)];
}

}  // namespace simdt3ceri

#endif /* SimdThreeCenterElectronRepulsionRsBufferRows_hpp */
