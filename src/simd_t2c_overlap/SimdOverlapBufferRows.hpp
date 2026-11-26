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


#ifndef SimdOverlapBufferRows_hpp
#define SimdOverlapBufferRows_hpp

#include <array>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"

namespace simdovl {  // simdovl namespace

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
        {       0,      10,      19,      33,      53,      80,     115},
        {      10,      40,      72,     114,     167,     232,     310},
        {      19,      75,     183,     282,     403,     547,     715},
        {      33,     123,     290,     619,     871,    1167,    1508},
        {      53,     185,     424,     886,    1717,    2278,    2920},
        {      80,     262,     586,    1204,    2302,    4132,    5260},
        {     115,     355,     777,    1574,    2977,    5295,    8928}
    }};

    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum < 7) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum < 7),
                              std::string("SimdOverlapBufferRows.number_of_buffer_rows: Angular momentum is out of range"));

    return rows[static_cast<size_t>(bra_angular_momentum)][static_cast<size_t>(ket_angular_momentum)];
}

}  // namespace simdovl

#endif /* SimdOverlapBufferRows_hpp */
