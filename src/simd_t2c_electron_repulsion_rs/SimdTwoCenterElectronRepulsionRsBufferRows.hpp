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


#ifndef SimdTwoCenterElectronRepulsionRsBufferRows_hpp
#define SimdTwoCenterElectronRepulsionRsBufferRows_hpp

#include <array>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

// NOTE: the name differs from the table of the unattenuated kernels on purpose.
// Both tables live in this namespace and both are inline, so a shared name is one
// definition of two different things: a translation unit which includes both does
// not compile, and one which includes either takes whichever the linker kept. The
// range separated buffer is about twice the unattenuated one, so taking the wrong
// table is a buffer overrun and not a wrong answer.

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
number_of_rs_buffer_rows(const int bra_angular_momentum, const int ket_angular_momentum) -> size_t
{
    constexpr std::array<std::array<size_t, 9>, 9> rows{{
        {       4,      19,      41,      81,     149,     255,     413,     637,     945},
        {      19,      63,     149,     279,     475,     755,    1139,    1649,    2309},
        {      41,     146,     378,     788,    1490,    2584,    4218,    6532,    9714},
        {      81,     282,     734,    1548,    2954,    5148,    8434,   13088,   19494},
        {     149,     487,    1591,    3371,    6351,   10941,   17679,   27115,   39927},
        {     255,     779,    2675,    5671,   10677,   18359,   29605,   45305,   66571},
        {     413,    1178,    4514,    9560,   17890,   30586,   49012,   74586,  109008},
        {     637,    1706,    6874,   14532,   27126,   46242,   73888,  112130,  163456},
        {     945,    2387,   10323,   21779,   40483,   68735,  109369,  165345,  240157}
    }};

    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum < 9) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum < 9),
                              std::string("SimdTwoCenterElectronRepulsionRsBufferRows.number_of_rs_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(bra_angular_momentum)][static_cast<size_t>(ket_angular_momentum)];
}

}  // namespace simdt2ceri

#endif /* SimdTwoCenterElectronRepulsionRsBufferRows_hpp */
