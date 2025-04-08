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

#ifndef TensorComponents_hpp
#define TensorComponents_hpp

#include <array>
#include <cstddef>
#include <numeric>

namespace tensor {  // tensor

/// @brief Gets compound number of Cartesian components of canonical tensors array.
/// @tparam N The size of canonical tensors array.
/// @param orders The array of cononical tensor orders.
/// @return The number of Cartesian components.
template <size_t N>
auto
number_of_cartesian_components(const std::array<int, N> &orders) -> int
{
    return std::accumulate(orders.begin(), orders.end(), 1, [](const int prod, const int order) { return prod * (order + 1) * (order + 2) / 2; });
}

/// @brief Gets compound number of spherical components of canonical tensors array.
/// @tparam N The size of canonical tensors array.
/// @param orders The array of cononical tensor orders.
/// @return The number of spherical components.
template <size_t N>
auto
number_of_spherical_components(const std::array<int, N> &orders) -> int
{
    return std::accumulate(orders.begin(), orders.end(), 1, [](const int prod, const int order) { return prod * (2 * order + 1); });
}

}  // namespace tensor

#endif /* TensorComponents_hpp */
