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

#ifndef CustomViews_hpp
#define CustomViews_hpp

#include <ranges>
#include <utility>

#include "CustomConstrains.hpp"

namespace views {  // views

/// @brief Creates indices (i,j) view of rectangular matrix.
/// @tparam T The integral indexing type.
/// @param nrows The number of rows in matrix.
/// @param ncolumns The number of columns in matrix.
/// @return The indices (i,j) view of rectangular matrix.
template <Integral T>
inline auto
rectangular(const T nrows, const T ncolumns) -> decltype(auto)
{
    auto indices = [=](const T index) { return std::make_pair(index / ncolumns, index % ncolumns); };

    return std::views::iota(T{0}, nrows * ncolumns) | std::views::transform(indices);
}

/// @brief Creates triangular indices (i,j) view of square matrix.
/// @tparam T The integral indexing type.
/// @param nrows The number of rows in matrix.
/// @return The triangular indices (i,j) view of square matrix.
template <Integral T>
inline auto
triangular(const T nrows) -> decltype(auto)
{
    auto rec_indices = [](const auto &index) { return index.first <= index.second; };

    return views::rectangular(nrows, nrows) | std::views::filter(rec_indices);
}

/// @brief Creates upper triangular indices (i,j) view of square matrix.
/// @tparam T The integral indexing type.
/// @param nrows The number of rows in matrix.
/// @return The upper triangular indices (i,j) view of square matrix.
template <Integral T>
inline auto
upper_triangular(const T nrows) -> decltype(auto)
{
    auto rec_up_indices = [](const auto &index) { return index.first < index.second; };

    return views::rectangular(nrows, nrows) | std::views::filter(rec_up_indices);
}

}  // namespace views

#endif /* CustomViews_hpp */
