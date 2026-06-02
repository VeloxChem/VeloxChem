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

#include "CauchySchwarzPairData.hpp"

#include <algorithm>
#include <utility>

namespace newints {

namespace {

/// @brief Lexicographic order on (i, j); the canonical layout of the pair list.
inline auto
pair_less(const ShellPair &a, const ShellPair &b) -> bool
{
    return (a.i != b.i) ? (a.i < b.i) : (a.j < b.j);
}

}  // namespace

CauchySchwarzPairData::CauchySchwarzPairData(std::size_t nshells, double threshold, std::vector<ShellPair> pairs)
    : nshells_(nshells)
    , threshold_(threshold)
    , pairs_(std::move(pairs))
{
    // keep the list sorted by (i, j) so q(i, j) can binary-search it
    std::sort(pairs_.begin(), pairs_.end(), pair_less);
}

auto
CauchySchwarzPairData::q(std::size_t i, std::size_t j) const -> double
{
    if (i < j) std::swap(i, j);  // canonical i >= j

    const ShellPair key{static_cast<int>(i), static_cast<int>(j), 0.0};

    const auto it = std::lower_bound(pairs_.begin(), pairs_.end(), key, pair_less);

    if (it != pairs_.end() && it->i == key.i && it->j == key.j) return it->q;

    return 0.0;
}

}  // namespace newints
