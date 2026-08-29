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


#ifndef SimdDimensions_hpp
#define SimdDimensions_hpp

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <string>
#include <vector>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the number of atom pairs surviving screening for each pair of
/// primitives of the basis functions on bra and ket sides.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param npairs The number of atom pairs to screen.
/// @param coordinates The coordinates of the atom pairs, as six rows of npairs
/// columns, ordered by ascending interatomic distance.
/// @param bound The integral bound, evaluated as bound(bra, iprim, ket, jprim,
/// distance).
/// @param threshold The screening threshold.
/// @return The vector of the number of surviving atom pairs of each pair of
/// primitives, with the primitives on bra side as the slowest running index.
/// @note The primitives of a basis function are kept sorted by descending
/// exponent, so the tightest primitives come first and the returned numbers do
/// not decrease along either index.
template <typename B>
inline auto
make_column_dimensions(const CBasisFunction &bra,
                       const CBasisFunction &ket,
                       const size_t          npairs,
                       const CSimdMatrix    &coordinates,
                       const B              &bound,
                       const double          threshold) -> std::vector<size_t>
{
    errors::assertMsgCritical(coordinates.number_of_rows() == 6,
                              std::string("SimdDimensions.make_column_dimensions: Coordinates must have six rows"));

    errors::assertMsgCritical(npairs <= coordinates.number_of_columns(),
                              std::string("SimdDimensions.make_column_dimensions: Number of atom pairs exceeds coordinates"));

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);
    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    // NOTE: the interatomic distance is recomputed from the coordinates of an
    // atom pair, so that the bisection below evaluates it a logarithmic number
    // of times rather than storing it for every atom pair.

    const auto distance = [&](const size_t i) {
        const auto rx = a_x[i] - b_x[i];

        const auto ry = a_y[i] - b_y[i];

        const auto rz = a_z[i] - b_z[i];

        return std::sqrt(rx * rx + ry * ry + rz * rz);
    };

    const auto nprim_a = bra.exponents().size();

    const auto nprim_b = ket.exponents().size();

    std::vector<size_t> dimensions;

    dimensions.reserve(nprim_a * nprim_b);

    std::ranges::for_each(std::views::iota(size_t{0}, nprim_a), [&](const auto i) {
        std::ranges::for_each(std::views::iota(size_t{0}, nprim_b), [&](const auto j) {
            dimensions.push_back(screenfunc::number_of_leading(
                npairs, [&](const size_t k) { return bound(bra, i, ket, j, distance(k)) > threshold; }));
        });
    });

    return dimensions;
}

/// @brief Creates the number of atom pairs surviving screening for each triple of
/// primitives of the basis functions on a, b and c sides.
/// @param a The basis function on a side.
/// @param b The basis function on b side.
/// @param c The basis function on c side.
/// @param npairs The number of atom pairs to screen.
/// @param coordinates The coordinates of the atom pairs on a and b sides, as six
/// rows of npairs columns, ordered by ascending interatomic distance.
/// @param bound The integral bound, evaluated as bound(a, iprim, b, jprim, c,
/// kprim, distance).
/// @param threshold The screening threshold.
/// @return The vector of the number of surviving atom pairs of each triple of
/// primitives, with the primitives on a side as the slowest and on c side as the
/// fastest running index.
/// @note Only the distance between the atoms on a and b sides enters, as the
/// dependence on the position of the atom on c side is neglected. The primitives
/// of a basis function are kept sorted by descending exponent, so the returned
/// numbers do not decrease along any of the three indices.
template <typename B>
inline auto
make_column_dimensions(const CBasisFunction &a,
                       const CBasisFunction &b,
                       const CBasisFunction &c,
                       const size_t          npairs,
                       const CSimdMatrix    &coordinates,
                       const B              &bound,
                       const double          threshold) -> std::vector<size_t>
{
    errors::assertMsgCritical(coordinates.number_of_rows() == 6,
                              std::string("SimdDimensions.make_column_dimensions: Coordinates must have six rows"));

    errors::assertMsgCritical(npairs <= coordinates.number_of_columns(),
                              std::string("SimdDimensions.make_column_dimensions: Number of atom pairs exceeds coordinates"));

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);
    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    const auto distance = [&](const size_t i) {
        const auto rx = a_x[i] - b_x[i];

        const auto ry = a_y[i] - b_y[i];

        const auto rz = a_z[i] - b_z[i];

        return std::sqrt(rx * rx + ry * ry + rz * rz);
    };

    const auto nprim_a = a.exponents().size();

    const auto nprim_b = b.exponents().size();

    const auto nprim_c = c.exponents().size();

    std::vector<size_t> dimensions;

    dimensions.reserve(nprim_a * nprim_b * nprim_c);

    std::ranges::for_each(std::views::iota(size_t{0}, nprim_a), [&](const auto i) {
        std::ranges::for_each(std::views::iota(size_t{0}, nprim_b), [&](const auto j) {
            std::ranges::for_each(std::views::iota(size_t{0}, nprim_c), [&](const auto k) {
                dimensions.push_back(screenfunc::number_of_leading(
                    npairs, [&](const size_t l) { return bound(a, i, b, j, c, k, distance(l)) > threshold; }));
            });
        });
    });

    return dimensions;
}

}  // namespace simdfunc

#endif /* SimdDimensions_hpp */
