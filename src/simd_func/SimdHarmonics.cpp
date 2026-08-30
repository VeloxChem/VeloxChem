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



#include "SimdHarmonics.hpp"

#include <cstddef>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_solid_harmonics(const CSimdMatrix &coordinates, const int lmax) -> std::vector<CSimdMatrix>
{
    errors::assertMsgCritical(lmax <= 12,
                              std::string("SimdHarmonics.make_solid_harmonics: Solid harmonics above angular momentum twelve are not implemented"));

    std::vector<CSimdMatrix> harmonics;

    if (lmax < 1) return harmonics;

    // NOTE: the storage is reserved before the harmonics are created, so that the
    // matrices of the lower angular momenta are not moved while they are read by
    // the recursion.

    harmonics.reserve(static_cast<size_t>(lmax));

    harmonics.push_back(make_p_solid_harmonics(coordinates));

    if (lmax > 1) harmonics.push_back(make_d_solid_harmonics(harmonics[0], coordinates));

    if (lmax > 2) harmonics.push_back(make_f_solid_harmonics(harmonics[1], harmonics[0], coordinates));

    if (lmax > 3) harmonics.push_back(make_g_solid_harmonics(harmonics[2], harmonics[1], harmonics[0], coordinates));

    if (lmax > 4) harmonics.push_back(make_h_solid_harmonics(harmonics[3], harmonics[2], harmonics[0], coordinates));

    if (lmax > 5) harmonics.push_back(make_i_solid_harmonics(harmonics[4], harmonics[3], harmonics[0], coordinates));

    if (lmax > 6) harmonics.push_back(make_k_solid_harmonics(harmonics[5], harmonics[4], harmonics[0], coordinates));

    if (lmax > 7) harmonics.push_back(make_l_solid_harmonics(harmonics[6], harmonics[5], harmonics[0], coordinates));

    if (lmax > 8) harmonics.push_back(make_m_solid_harmonics(harmonics[7], harmonics[6], harmonics[0], coordinates));

    if (lmax > 9) harmonics.push_back(make_n_solid_harmonics(harmonics[8], harmonics[7], harmonics[0], coordinates));

    if (lmax > 10) harmonics.push_back(make_o_solid_harmonics(harmonics[9], harmonics[8], harmonics[0], coordinates));

    if (lmax > 11) harmonics.push_back(make_q_solid_harmonics(harmonics[10], harmonics[9], harmonics[0], coordinates));

    return harmonics;
}

}  // namespace simdfunc
