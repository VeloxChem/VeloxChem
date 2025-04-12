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

#ifndef CantorFunc_hpp
#define CantorFunc_hpp

#include <cmath>

#include "T2Pair.hpp"

namespace mathfunc {  // mathfunc namespace

/**
 Computes Cantor pairing index for two integer non-negative numbers.

 @return the Cantor index.
 */
inline auto
getCantorIndex(const T2Pair& pair) -> int64_t
{
    const auto x = pair.first;

    const auto y = pair.second;

    return (x + y + 1) * (x + y) / 2 + y;
}

/**
 Reduces Cantor pairing index to Cantor pair i.e. two integer non-negative numbers.

 @return the Cantor pair.
 */
inline auto
getCantorPair(const int64_t index) -> T2Pair
{
    const auto w = static_cast<int64_t>(std::floor(0.5 * (std::sqrt(8 * index + 1) - 1)));

    const auto y = index - (w * w + w) / 2;

    return {w - y, y};
}

}  // namespace mathfunc

#endif /* CantorFunc_hpp */
