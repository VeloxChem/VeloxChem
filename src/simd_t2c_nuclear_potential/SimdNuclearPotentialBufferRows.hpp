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


#ifndef SimdNuclearPotentialBufferRows_hpp
#define SimdNuclearPotentialBufferRows_hpp

#include <array>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"

namespace simdnpot {  // simdnpot namespace

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
    constexpr std::array<std::array<size_t, 7>, 7> rows{{
        {       6,      16,      29,      53,      93,     155,     246},
        {      16,      50,      92,     154,     242,     363,     525},
        {      29,      95,     223,     357,     534,     762,    1050},
        {      53,     163,     365,     750,    1086,    1502,    2008},
        {      93,     260,     555,    1101,    2052,    2778,    3640},
        {     155,     393,     801,    1539,    2802,    4852,    6266},
        {     246,     570,    1112,    2074,    3697,    6301,   10298}
    }};

    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum < 7) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum < 7),
                              std::string("SimdNuclearPotentialBufferRows.number_of_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(bra_angular_momentum)][static_cast<size_t>(ket_angular_momentum)];
}

}  // namespace simdnpot

#endif /* SimdNuclearPotentialBufferRows_hpp */
