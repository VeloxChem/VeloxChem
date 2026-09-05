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



#include "SimdOverlapRecHH.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hh_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHH.compute_hh_overlap: Basis functions must be of angular momenta five and five"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHH.compute_hh_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto nprims = bra.exponents().size() * ket.exponents().size();

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 6);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 121 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);

    // NOTE: the components of the vector between the atoms and its squared length
    // are carried by the coordinates, so the angular half below reads rows which
    // are already in place.

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ab_2 = coordinates.data(9);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each term over the pairs of primitives

    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {
        const auto ncols = pair.ncols;

        const auto fexp = pair.aexp + pair.bexp;

        const auto fmu = pair.aexp * pair.bexp / fexp;

        const auto fovl = fpi / fexp;

        const auto fbase = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);

        // NOTE: the Gaussian product center is displaced from the atom on bra side
        // by fal times the vector between the atoms and from the atom on ket side by
        // fbe times it, and fh is the second moment the integration over that center
        // leaves behind.

        const auto fal = -pair.bexp / fexp;

        const auto fbe = pair.aexp / fexp;

        const auto fh = 0.5 / fexp;

        const auto f_0 = fbase * fal * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fal * fbe * fbe * fbe * fh * fh;

        const auto f_3 = fbase * fal * fal * fbe * fbe * fh * fh * fh;

        const auto f_4 = fbase * fal * fbe * fh * fh * fh * fh;

        const auto f_5 = fbase * fh * fh * fh * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
            pe_3[k] += f_3 * fss;
            pe_4[k] += f_4 * fss;
            pe_5[k] += f_5 * fss;
        }
    });

    // NOTE: the rows of the values are not aligned, as they start at the offset of
    // this combination of basis functions in the values block, so they are kept out
    // of the aligned clauses below.

    auto *pc_0 = values + 0 * nvalues;
    auto *pc_1 = values + 1 * nvalues;
    auto *pc_2 = values + 2 * nvalues;
    auto *pc_3 = values + 3 * nvalues;
    auto *pc_4 = values + 4 * nvalues;
    auto *pc_5 = values + 5 * nvalues;
    auto *pc_6 = values + 6 * nvalues;
    auto *pc_7 = values + 7 * nvalues;
    auto *pc_8 = values + 8 * nvalues;
    auto *pc_9 = values + 9 * nvalues;
    auto *pc_10 = values + 10 * nvalues;
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 18 * nvalues;
    auto *pc_18 = values + 19 * nvalues;
    auto *pc_19 = values + 20 * nvalues;
    auto *pc_20 = values + 21 * nvalues;
    auto *pc_21 = values + 24 * nvalues;
    auto *pc_22 = values + 25 * nvalues;
    auto *pc_23 = values + 26 * nvalues;
    auto *pc_24 = values + 27 * nvalues;
    auto *pc_25 = values + 28 * nvalues;
    auto *pc_26 = values + 29 * nvalues;
    auto *pc_27 = values + 30 * nvalues;
    auto *pc_28 = values + 31 * nvalues;
    auto *pc_29 = values + 32 * nvalues;
    auto *pc_30 = values + 36 * nvalues;
    auto *pc_31 = values + 37 * nvalues;
    auto *pc_32 = values + 38 * nvalues;
    auto *pc_33 = values + 39 * nvalues;
    auto *pc_34 = values + 40 * nvalues;
    auto *pc_35 = values + 41 * nvalues;
    auto *pc_36 = values + 42 * nvalues;
    auto *pc_37 = values + 43 * nvalues;
    auto *pc_38 = values + 48 * nvalues;
    auto *pc_39 = values + 49 * nvalues;
    auto *pc_40 = values + 50 * nvalues;
    auto *pc_41 = values + 51 * nvalues;
    auto *pc_42 = values + 52 * nvalues;
    auto *pc_43 = values + 53 * nvalues;
    auto *pc_44 = values + 54 * nvalues;
    auto *pc_45 = values + 60 * nvalues;
    auto *pc_46 = values + 61 * nvalues;
    auto *pc_47 = values + 62 * nvalues;
    auto *pc_48 = values + 63 * nvalues;
    auto *pc_49 = values + 64 * nvalues;
    auto *pc_50 = values + 65 * nvalues;
    auto *pc_51 = values + 72 * nvalues;
    auto *pc_52 = values + 73 * nvalues;
    auto *pc_53 = values + 74 * nvalues;
    auto *pc_54 = values + 75 * nvalues;
    auto *pc_55 = values + 76 * nvalues;
    auto *pc_56 = values + 84 * nvalues;
    auto *pc_57 = values + 85 * nvalues;
    auto *pc_58 = values + 86 * nvalues;
    auto *pc_59 = values + 87 * nvalues;
    auto *pc_60 = values + 96 * nvalues;
    auto *pc_61 = values + 97 * nvalues;
    auto *pc_62 = values + 98 * nvalues;
    auto *pc_63 = values + 108 * nvalues;
    auto *pc_64 = values + 109 * nvalues;
    auto *pc_65 = values + 120 * nvalues;

    // NOTE: the components are formed in 17 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_0[k] = e_0 * (12.3046875 * x * x * x * x * x * x * x * x * y * y - 49.21875 * x * x * x * x * x * x * y * y * y * y + 54.140625 * x * x * x * x * y * y * y * y * y * y - 9.84375 * x * x * y * y * y * y * y * y * y * y + 0.4921875 * y * y * y * y * y * y * y * y * y * y) + e_1 * (12.3046875 * x * x * x * x * x * x * x * x + 49.21875 * x * x * x * x * x * x * y * y + 73.828125 * x * x * x * x * y * y * y * y + 49.21875 * x * x * y * y * y * y * y * y + 12.3046875 * y * y * y * y * y * y * y * y) + e_2 * (196.875 * x * x * x * x * x * x + 590.625 * x * x * x * x * y * y + 590.625 * x * x * y * y * y * y + 196.875 * y * y * y * y * y * y) + e_3 * (1181.25 * x * x * x * x + 2362.5 * x * x * y * y + 1181.25 * y * y * y * y) + e_4 * (2362.5 * x * x + 2362.5 * y * y) + e_5 * (945.0);

        pc_1[k] = e_0 * (std::sqrt(968.994140625) * x * x * x * x * x * x * x * y * y * z - std::sqrt(8720.947265625) * x * x * x * x * x * y * y * y * y * z + std::sqrt(4689.931640625) * x * x * x * y * y * y * y * y * y * z - std::sqrt(38.759765625) * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(968.994140625) * x * x * x * x * x * x * x * z + std::sqrt(8720.947265625) * x * x * x * x * x * y * y * z + std::sqrt(8720.947265625) * x * x * x * y * y * y * y * z + std::sqrt(968.994140625) * x * y * y * y * y * y * y * z) + e_2 * (std::sqrt(139535.15625) * x * x * x * x * x * z + std::sqrt(558140.625) * x * x * x * y * y * z + std::sqrt(139535.15625) * x * y * y * y * y * z) + e_3 * (std::sqrt(2232562.5) * x * x * x * z + std::sqrt(2232562.5) * x * y * y * z) + e_4 * (std::sqrt(2232562.5) * x * z);

        pc_2[k] = e_0 * (-std::sqrt(30.28106689453125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(53.8330078125) * x * x * x * x * x * x * y * y * y * y + std::sqrt(1937.98828125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(65.137939453125) * x * x * x * x * y * y * y * y * y * y - std::sqrt(10551.26953125) * x * x * x * x * y * y * y * y * z * z - std::sqrt(19.3798828125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(1455.64453125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(0.13458251953125) * y * y * y * y * y * y * y * y * y * y - std::sqrt(8.61328125) * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(30.28106689453125) * x * x * x * x * x * x * x * x - std::sqrt(4360.4736328125) * x * x * x * x * x * x * y * y + std::sqrt(1937.98828125) * x * x * x * x * x * x * z * z + std::sqrt(16486.358642578125) * x * x * x * x * y * y * y * y + std::sqrt(1937.98828125) * x * x * x * x * y * y * z * z - std::sqrt(1345.8251953125) * x * x * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * y * y * y * y * z * z + std::sqrt(84.11407470703125) * y * y * y * y * y * y * y * y - std::sqrt(1937.98828125) * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(7751.953125) * x * x * x * x * x * x - std::sqrt(7751.953125) * x * x * x * x * y * y + std::sqrt(124031.25) * x * x * x * x * z * z + std::sqrt(7751.953125) * x * x * y * y * y * y + std::sqrt(7751.953125) * y * y * y * y * y * y - std::sqrt(124031.25) * y * y * y * y * z * z) + e_3 * (-std::sqrt(124031.25) * x * x * x * x + std::sqrt(496125.0) * x * x * z * z + std::sqrt(124031.25) * y * y * y * y - std::sqrt(496125.0) * y * y * z * z) + e_4 * (-std::sqrt(124031.25) * x * x + std::sqrt(124031.25) * y * y);

        pc_3[k] = e_0 * (-std::sqrt(322.998046875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(322.998046875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1291.9921875) * x * x * x * x * x * y * y * z * z * z + std::sqrt(1046.513671875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(5167.96875) * x * x * x * y * y * y * y * z * z * z - std::sqrt(12.919921875) * x * y * y * y * y * y * y * y * y * z + std::sqrt(51.6796875) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(322.998046875) * x * x * x * x * x * x * x * z - std::sqrt(26162.841796875) * x * x * x * x * x * y * y * z + std::sqrt(1291.9921875) * x * x * x * x * x * z * z * z + std::sqrt(201873.779296875) * x * x * x * y * y * y * y * z - std::sqrt(5167.96875) * x * x * x * y * y * z * z * z + std::sqrt(322.998046875) * x * y * y * y * y * y * y * z - std::sqrt(11627.9296875) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(46511.71875) * x * x * x * x * x * z + std::sqrt(186046.875) * x * x * x * y * y * z + std::sqrt(20671.875) * x * x * x * z * z * z + std::sqrt(418605.46875) * x * y * y * y * y * z - std::sqrt(186046.875) * x * y * y * z * z * z) + e_3 * (-std::sqrt(186046.875) * x * x * x * z + std::sqrt(1674421.875) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_4[k] = e_0 * (std::sqrt(2.8839111328125) * x * x * x * x * x * x * x * x * y * y - std::sqrt(415.283203125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(22.60986328125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(415.283203125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(184.5703125) * x * x * x * x * y * y * z * z * z * z - std::sqrt(7.3828125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(1345.517578125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(738.28125) * x * x * y * y * y * y * z * z * z * z + std::sqrt(0.1153564453125) * y * y * y * y * y * y * y * y * y * y - std::sqrt(16.611328125) * y * y * y * y * y * y * y * y * z * z + std::sqrt(7.3828125) * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(2.8839111328125) * x * x * x * x * x * x * x * x + std::sqrt(738.28125) * x * x * x * x * x * x * y * y - std::sqrt(415.283203125) * x * x * x * x * x * x * z * z - std::sqrt(2595.52001953125) * x * x * x * x * y * y * y * y - std::sqrt(10382.080078125) * x * x * x * x * y * y * z * z + std::sqrt(184.5703125) * x * x * x * x * z * z * z * z - std::sqrt(4614.2578125) * x * x * y * y * y * y * y * y + std::sqrt(259552.001953125) * x * x * y * y * y * y * z * z - std::sqrt(6644.53125) * x * x * y * y * z * z * z * z + std::sqrt(72.0977783203125) * y * y * y * y * y * y * y * y - std::sqrt(3737.548828125) * y * y * y * y * y * y * z * z + std::sqrt(184.5703125) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(738.28125) * x * x * x * x * x * x - std::sqrt(26578.125) * x * x * x * x * z * z - std::sqrt(166113.28125) * x * x * y * y * y * y + std::sqrt(956812.5) * x * x * y * y * z * z + std::sqrt(2953.125) * y * y * y * y * y * y - std::sqrt(26578.125) * y * y * y * y * z * z) + e_3 * (std::sqrt(6644.53125) * x * x * x * x - std::sqrt(239203.125) * x * x * y * y + std::sqrt(6644.53125) * y * y * y * y);

        pc_5[k] = e_0 * (std::sqrt(43.2586669921875) * x * x * x * x * x * x * x * x * y * z - std::sqrt(307.6171875) * x * x * x * x * x * x * y * z * z * z - std::sqrt(339.14794921875) * x * x * x * x * y * y * y * y * y * z + std::sqrt(307.6171875) * x * x * x * x * y * y * y * z * z * z + std::sqrt(12.3046875) * x * x * x * x * y * z * z * z * z * z - std::sqrt(110.7421875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(996.6796875) * x * x * y * y * y * y * y * z * z * z - std::sqrt(49.21875) * x * x * y * y * y * z * z * z * z * z + std::sqrt(1.7303466796875) * y * y * y * y * y * y * y * y * y * z - std::sqrt(12.3046875) * y * y * y * y * y * y * y * z * z * z + std::sqrt(0.4921875) * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(17303.466796875) * x * x * x * x * x * x * y * z - std::sqrt(17303.466796875) * x * x * x * x * y * y * y * z - std::sqrt(30761.71875) * x * x * x * x * y * z * z * z - std::sqrt(56063.232421875) * x * x * y * y * y * y * y * z + std::sqrt(123046.875) * x * x * y * y * y * z * z * z + std::sqrt(692.138671875) * y * y * y * y * y * y * y * z - std::sqrt(1230.46875) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(276855.46875) * x * x * x * x * y * z - std::sqrt(1107421.875) * x * x * y * y * y * z + std::sqrt(11074.21875) * y * y * y * y * y * z);

        pc_6[k] = e_0 * (std::sqrt(2.8839111328125) * x * x * x * x * x * x * x * x * x * y - std::sqrt(415.283203125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(22.60986328125) * x * x * x * x * x * y * y * y * y * y + std::sqrt(415.283203125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(184.5703125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(7.3828125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(1345.517578125) * x * x * x * y * y * y * y * y * z * z - std::sqrt(738.28125) * x * x * x * y * y * y * z * z * z * z + std::sqrt(0.1153564453125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(16.611328125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(7.3828125) * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(1661.1328125) * x * x * x * x * x * x * x * y - std::sqrt(738.28125) * x * x * x * x * x * y * y * y - std::sqrt(81395.5078125) * x * x * x * x * x * y * z * z - std::sqrt(4614.2578125) * x * x * x * y * y * y * y * y + std::sqrt(166113.28125) * x * x * x * y * y * y * z * z + std::sqrt(2953.125) * x * x * x * y * z * z * z * z + std::sqrt(1661.1328125) * x * y * y * y * y * y * z * z - std::sqrt(2953.125) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(59800.78125) * x * x * x * x * x * y - std::sqrt(73828.125) * x * x * x * y * y * y - std::sqrt(425250.0) * x * x * x * y * z * z - std::sqrt(6644.53125) * x * y * y * y * y * y + std::sqrt(425250.0) * x * y * y * y * z * z) + e_3 * (std::sqrt(106312.5) * x * x * x * y - std::sqrt(106312.5) * x * y * y * y);

        pc_7[k] = e_0 * (-std::sqrt(80.74951171875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(322.998046875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(322.998046875) * x * x * x * x * x * x * y * z * z * z + std::sqrt(51.6796875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(2906.982421875) * x * x * x * x * y * y * y * z * z * z - std::sqrt(322.998046875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1563.310546875) * x * x * y * y * y * y * y * z * z * z + std::sqrt(3.22998046875) * y * y * y * y * y * y * y * y * y * z - std::sqrt(12.919921875) * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(20671.875) * x * x * x * x * x * x * y * z + std::sqrt(32299.8046875) * x * x * x * x * y * y * y * z + std::sqrt(11627.9296875) * x * x * x * x * y * z * z * z - std::sqrt(46511.71875) * x * x * y * y * y * y * y * z + std::sqrt(5167.96875) * x * x * y * y * y * z * z * z + std::sqrt(1291.9921875) * y * y * y * y * y * y * y * z - std::sqrt(1291.9921875) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(418605.46875) * x * x * x * x * y * z - std::sqrt(186046.875) * x * x * y * y * y * z + std::sqrt(186046.875) * x * x * y * z * z * z + std::sqrt(46511.71875) * y * y * y * y * y * z - std::sqrt(20671.875) * y * y * y * z * z * z) + e_3 * (-std::sqrt(1674421.875) * x * x * y * z + std::sqrt(186046.875) * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_8[k] = e_0 * (-std::sqrt(3.36456298828125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(53.8330078125) * x * x * x * x * x * x * x * y * y * y + std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(4.844970703125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(5383.30078125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(105.5126953125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(8277.36328125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(1.21124267578125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(77.51953125) * x * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(861.328125) * x * x * x * x * x * x * x * y + std::sqrt(3445.3125) * x * x * x * x * x * y * y * y + std::sqrt(7751.953125) * x * x * x * x * x * y * z * z - std::sqrt(21533.203125) * x * x * x * y * y * y * y * y + std::sqrt(31007.8125) * x * x * x * y * y * y * z * z + std::sqrt(7751.953125) * x * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(31007.8125) * x * x * x * x * x * y - std::sqrt(124031.25) * x * x * x * y * y * y + std::sqrt(496125.0) * x * x * x * y * z * z - std::sqrt(31007.8125) * x * y * y * y * y * y + std::sqrt(496125.0) * x * y * y * y * z * z) + e_3 * (-std::sqrt(496125.0) * x * x * x * y - std::sqrt(496125.0) * x * y * y * y + std::sqrt(1984500.0) * x * y * z * z) + e_4 * (-std::sqrt(496125.0) * x * y);

        pc_9[k] = e_0 * (std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(3875.9765625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(10552.34619140625) * x * x * x * x * y * y * y * y * y * z - std::sqrt(620.15625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2.4224853515625) * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(968.994140625) * x * x * x * x * x * x * y * z + std::sqrt(8720.947265625) * x * x * x * x * y * y * y * z + std::sqrt(8720.947265625) * x * x * y * y * y * y * y * z + std::sqrt(968.994140625) * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(139535.15625) * x * x * x * x * y * z + std::sqrt(558140.625) * x * x * y * y * y * z + std::sqrt(139535.15625) * y * y * y * y * y * z) + e_3 * (std::sqrt(2232562.5) * x * x * y * z + std::sqrt(2232562.5) * y * y * y * z) + e_4 * (std::sqrt(2232562.5) * y * z);

        pc_10[k] = e_0 * (2.4609375 * x * x * x * x * x * x * x * x * x * y - 29.53125 * x * x * x * x * x * x * x * y * y * y + 62.015625 * x * x * x * x * x * y * y * y * y * y - 29.53125 * x * x * x * y * y * y * y * y * y * y + 2.4609375 * x * y * y * y * y * y * y * y * y * y);

        pc_11[k] = e_0 * (78.75 * x * x * x * x * x * x * y * y * z * z - 157.5 * x * x * x * x * y * y * y * y * z * z + 78.75 * x * x * y * y * y * y * y * y * z * z) + e_1 * (78.75 * x * x * x * x * x * x * y * y + 78.75 * x * x * x * x * x * x * z * z - 157.5 * x * x * x * x * y * y * y * y + 236.25 * x * x * x * x * y * y * z * z + 78.75 * x * x * y * y * y * y * y * y + 236.25 * x * x * y * y * y * y * z * z + 78.75 * y * y * y * y * y * y * z * z) + e_2 * (78.75 * x * x * x * x * x * x + 236.25 * x * x * x * x * y * y + 708.75 * x * x * x * x * z * z + 236.25 * x * x * y * y * y * y + 1417.5 * x * x * y * y * z * z + 78.75 * y * y * y * y * y * y + 708.75 * y * y * y * y * z * z) + e_3 * (708.75 * x * x * x * x + 1417.5 * x * x * y * y + 1890.0 * x * x * z * z + 708.75 * y * y * y * y + 1890.0 * y * y * z * z) + e_4 * (1890.0 * x * x + 1890.0 * y * y + 945.0 * z * z) + e_5 * (945.0);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_12[k] = e_0 * (-std::sqrt(193.798828125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(21.533203125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(12403.125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(193.798828125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(22050.0) * x * x * x * y * y * y * y * z * z * z - std::sqrt(21.533203125) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1378.125) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(193.798828125) * x * x * x * x * x * x * x * z + std::sqrt(4844.970703125) * x * x * x * x * x * y * y * z + std::sqrt(12403.125) * x * x * x * x * x * z * z * z - std::sqrt(36197.314453125) * x * x * x * y * y * y * y * z + std::sqrt(49612.5) * x * x * x * y * y * z * z * z + std::sqrt(538.330078125) * x * y * y * y * y * y * y * z + std::sqrt(12403.125) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * x * z + std::sqrt(12403.125) * x * x * x * y * y * z + std::sqrt(446512.5) * x * x * x * z * z * z + std::sqrt(3100.78125) * x * y * y * y * y * z + std::sqrt(446512.5) * x * y * y * z * z * z) + e_3 * (std::sqrt(793800.0) * x * x * x * z + std::sqrt(793800.0) * x * y * y * z + std::sqrt(793800.0) * x * z * z * z) + e_4 * (std::sqrt(2431012.5) * x * z);

        pc_13[k] = e_0 * (-std::sqrt(2067.1875) * x * x * x * x * x * x * y * y * z * z + std::sqrt(8268.75) * x * x * x * x * y * y * z * z * z * z + std::sqrt(2067.1875) * x * x * y * y * y * y * y * y * z * z - std::sqrt(8268.75) * x * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(2067.1875) * x * x * x * x * x * x * y * y - std::sqrt(2067.1875) * x * x * x * x * x * x * z * z - std::sqrt(18604.6875) * x * x * x * x * y * y * z * z + std::sqrt(8268.75) * x * x * x * x * z * z * z * z + std::sqrt(2067.1875) * x * x * y * y * y * y * y * y + std::sqrt(18604.6875) * x * x * y * y * y * y * z * z + std::sqrt(2067.1875) * y * y * y * y * y * y * z * z - std::sqrt(8268.75) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(2067.1875) * x * x * x * x * x * x - std::sqrt(167442.1875) * x * x * x * x * y * y - std::sqrt(18604.6875) * x * x * x * x * z * z + std::sqrt(167442.1875) * x * x * y * y * y * y + std::sqrt(74418.75) * x * x * z * z * z * z + std::sqrt(2067.1875) * y * y * y * y * y * y + std::sqrt(18604.6875) * y * y * y * y * z * z - std::sqrt(74418.75) * y * y * z * z * z * z) + e_3 * (-std::sqrt(167442.1875) * x * x * x * x + std::sqrt(74418.75) * x * x * z * z + std::sqrt(167442.1875) * y * y * y * y - std::sqrt(74418.75) * y * y * z * z) + e_4 * (-std::sqrt(297675.0) * x * x + std::sqrt(297675.0) * y * y);

        pc_14[k] = e_0 * (std::sqrt(18.45703125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(18.45703125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(2657.8125) * x * x * x * x * x * y * y * z * z * z - std::sqrt(18.45703125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(1181.25) * x * x * x * y * y * z * z * z * z * z - std::sqrt(18.45703125) * x * y * y * y * y * y * y * y * y * z + std::sqrt(2657.8125) * x * y * y * y * y * y * y * z * z * z - std::sqrt(1181.25) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(18.45703125) * x * x * x * x * x * x * x * z - std::sqrt(1495.01953125) * x * x * x * x * x * y * y * z - std::sqrt(2657.8125) * x * x * x * x * x * z * z * z - std::sqrt(461.42578125) * x * x * x * y * y * y * y * z - std::sqrt(29531.25) * x * x * x * y * y * z * z * z + std::sqrt(1181.25) * x * x * x * z * z * z * z * z + std::sqrt(461.42578125) * x * y * y * y * y * y * y * z + std::sqrt(184570.3125) * x * y * y * y * y * z * z * z - std::sqrt(10631.25) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2657.8125) * x * x * x * x * x * z - std::sqrt(265781.25) * x * x * x * y * y * z - std::sqrt(29531.25) * x * x * x * z * z * z + std::sqrt(598007.8125) * x * y * y * y * y * z + std::sqrt(265781.25) * x * y * y * z * z * z) + e_3 * (-std::sqrt(265781.25) * x * x * x * z + std::sqrt(2392031.25) * x * y * y * z);

        pc_15[k] = e_0 * (std::sqrt(276.85546875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(276.85546875) * x * x * x * x * x * y * y * y * z * z - std::sqrt(1968.75) * x * x * x * x * x * y * z * z * z * z - std::sqrt(276.85546875) * x * x * x * y * y * y * y * y * z * z + std::sqrt(78.75) * x * x * x * y * z * z * z * z * z * z - std::sqrt(276.85546875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(1968.75) * x * y * y * y * y * y * z * z * z * z - std::sqrt(78.75) * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(276.85546875) * x * x * x * x * x * x * x * y + std::sqrt(276.85546875) * x * x * x * x * x * y * y * y + std::sqrt(17718.75) * x * x * x * x * x * y * z * z - std::sqrt(276.85546875) * x * x * x * y * y * y * y * y - std::sqrt(96468.75) * x * x * x * y * z * z * z * z - std::sqrt(276.85546875) * x * y * y * y * y * y * y * y - std::sqrt(17718.75) * x * y * y * y * y * y * z * z + std::sqrt(96468.75) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(70875.0) * x * x * x * x * x * y - std::sqrt(70875.0) * x * x * x * y * z * z - std::sqrt(70875.0) * x * y * y * y * y * y + std::sqrt(70875.0) * x * y * y * y * z * z) + e_3 * (std::sqrt(637875.0) * x * x * x * y - std::sqrt(637875.0) * x * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_16[k] = e_0 * (std::sqrt(18.45703125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(18.45703125) * x * x * x * x * x * x * y * y * y * z - std::sqrt(2657.8125) * x * x * x * x * x * x * y * z * z * z - std::sqrt(18.45703125) * x * x * x * x * y * y * y * y * y * z + std::sqrt(1181.25) * x * x * x * x * y * z * z * z * z * z - std::sqrt(18.45703125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2657.8125) * x * x * y * y * y * y * y * z * z * z - std::sqrt(1181.25) * x * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(461.42578125) * x * x * x * x * x * x * y * z + std::sqrt(461.42578125) * x * x * x * x * y * y * y * z - std::sqrt(184570.3125) * x * x * x * x * y * z * z * z + std::sqrt(1495.01953125) * x * x * y * y * y * y * y * z + std::sqrt(29531.25) * x * x * y * y * y * z * z * z + std::sqrt(10631.25) * x * x * y * z * z * z * z * z - std::sqrt(18.45703125) * y * y * y * y * y * y * y * z + std::sqrt(2657.8125) * y * y * y * y * y * z * z * z - std::sqrt(1181.25) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(598007.8125) * x * x * x * x * y * z + std::sqrt(265781.25) * x * x * y * y * y * z - std::sqrt(265781.25) * x * x * y * z * z * z + std::sqrt(2657.8125) * y * y * y * y * y * z + std::sqrt(29531.25) * y * y * y * z * z * z) + e_3 * (-std::sqrt(2392031.25) * x * x * y * z + std::sqrt(265781.25) * y * y * y * z);

        pc_17[k] = e_0 * (-std::sqrt(516.796875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(516.796875) * x * x * x * x * x * y * y * y * z * z + std::sqrt(2067.1875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(516.796875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(8268.75) * x * x * x * y * y * y * z * z * z * z - std::sqrt(516.796875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(2067.1875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(516.796875) * x * x * x * x * x * x * x * y + std::sqrt(516.796875) * x * x * x * x * x * y * y * y - std::sqrt(18604.6875) * x * x * x * x * x * y * z * z + std::sqrt(516.796875) * x * x * x * y * y * y * y * y - std::sqrt(8268.75) * x * x * x * y * y * y * z * z + std::sqrt(33075.0) * x * x * x * y * z * z * z * z - std::sqrt(516.796875) * x * y * y * y * y * y * y * y - std::sqrt(18604.6875) * x * y * y * y * y * y * z * z + std::sqrt(33075.0) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(74418.75) * x * x * x * x * x * y + std::sqrt(33075.0) * x * x * x * y * y * y - std::sqrt(74418.75) * x * x * x * y * z * z - std::sqrt(74418.75) * x * y * y * y * y * y - std::sqrt(74418.75) * x * y * y * y * z * z + std::sqrt(297675.0) * x * y * z * z * z * z) + e_3 * (-std::sqrt(669768.75) * x * x * x * y - std::sqrt(669768.75) * x * y * y * y + std::sqrt(297675.0) * x * y * z * z) + e_4 * (-std::sqrt(1190700.0) * x * y);

        pc_18[k] = e_0 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(193.798828125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1378.125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(21.533203125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(22050.0) * x * x * x * x * y * y * y * z * z * z - std::sqrt(193.798828125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(12403.125) * x * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(538.330078125) * x * x * x * x * x * x * y * z - std::sqrt(36197.314453125) * x * x * x * x * y * y * y * z + std::sqrt(12403.125) * x * x * x * x * y * z * z * z + std::sqrt(4844.970703125) * x * x * y * y * y * y * y * z + std::sqrt(49612.5) * x * x * y * y * y * z * z * z - std::sqrt(193.798828125) * y * y * y * y * y * y * y * z + std::sqrt(12403.125) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * y * z + std::sqrt(12403.125) * x * x * y * y * y * z + std::sqrt(446512.5) * x * x * y * z * z * z + std::sqrt(3100.78125) * y * y * y * y * y * z + std::sqrt(446512.5) * y * y * y * z * z * z) + e_3 * (std::sqrt(793800.0) * x * x * y * z + std::sqrt(793800.0) * y * y * y * z + std::sqrt(793800.0) * y * z * z * z) + e_4 * (std::sqrt(2431012.5) * y * z);

        pc_19[k] = e_0 * (19.6875 * x * x * x * x * x * x * x * y * z * z - 137.8125 * x * x * x * x * x * y * y * y * z * z + 137.8125 * x * x * x * y * y * y * y * y * z * z - 19.6875 * x * y * y * y * y * y * y * y * z * z) + e_1 * (19.6875 * x * x * x * x * x * x * x * y - 137.8125 * x * x * x * x * x * y * y * y + 137.8125 * x * x * x * y * y * y * y * y - 19.6875 * x * y * y * y * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_20[k] = e_0 * (std::sqrt(38.759765625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(4689.931640625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(8720.947265625) * x * x * x * x * y * y * y * y * y * z - std::sqrt(968.994140625) * x * x * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(968.994140625) * x * x * x * x * x * x * y * z - std::sqrt(8720.947265625) * x * x * x * x * y * y * y * z - std::sqrt(8720.947265625) * x * x * y * y * y * y * y * z - std::sqrt(968.994140625) * y * y * y * y * y * y * y * z) + e_2 * (-std::sqrt(139535.15625) * x * x * x * x * y * z - std::sqrt(558140.625) * x * x * y * y * y * z - std::sqrt(139535.15625) * y * y * y * y * y * z) + e_3 * (-std::sqrt(2232562.5) * x * x * y * z - std::sqrt(2232562.5) * y * y * y * z) + e_4 * (-std::sqrt(2232562.5) * y * z);

        pc_21[k] = e_0 * (2.4609375 * x * x * x * x * x * x * x * x * y * y + 3.28125 * x * x * x * x * x * x * y * y * y * y - 39.375 * x * x * x * x * x * x * y * y * z * z - 0.546875 * x * x * x * x * y * y * y * y * y * y - 13.125 * x * x * x * x * y * y * y * y * z * z + 157.5 * x * x * x * x * y * y * z * z * z * z - 1.09375 * x * x * y * y * y * y * y * y * y * y + 21.875 * x * x * y * y * y * y * y * y * z * z - 105.0 * x * x * y * y * y * y * z * z * z * z + 0.2734375 * y * y * y * y * y * y * y * y * y * y - 4.375 * y * y * y * y * y * y * y * y * z * z + 17.5 * y * y * y * y * y * y * z * z * z * z) + e_1 * (2.4609375 * x * x * x * x * x * x * x * x + 49.21875 * x * x * x * x * x * x * y * y - 39.375 * x * x * x * x * x * x * z * z + 27.890625 * x * x * x * x * y * y * y * y + 275.625 * x * x * x * x * y * y * z * z + 157.5 * x * x * x * x * z * z * z * z - 12.03125 * x * x * y * y * y * y * y * y - 380.625 * x * x * y * y * y * y * z * z + 315.0 * x * x * y * y * z * z * z * z + 6.8359375 * y * y * y * y * y * y * y * y + 4.375 * y * y * y * y * y * y * z * z + 157.5 * y * y * y * y * z * z * z * z) + e_2 * (39.375 * x * x * x * x * x * x + 590.625 * x * x * x * x * y * y + 315.0 * x * x * x * x * z * z - 196.875 * x * x * y * y * y * y + 630.0 * x * x * y * y * z * z + 630.0 * x * x * z * z * z * z + 91.875 * y * y * y * y * y * y + 315.0 * y * y * y * y * z * z + 630.0 * y * y * z * z * z * z) + e_3 * (498.75 * x * x * x * x + 997.5 * x * x * y * y + 2100.0 * x * x * z * z + 498.75 * y * y * y * y + 2100.0 * y * y * z * z + 420.0 * z * z * z * z) + e_4 * (1522.5 * x * x + 1522.5 * y * y + 1680.0 * z * z) + e_5 * (945.0);

        pc_22[k] = e_0 * (std::sqrt(64.599609375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(179.443359375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(6459.9609375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(7.177734375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(2871.09375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(16537.5) * x * x * x * y * y * z * z * z * z * z - std::sqrt(7.177734375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(717.7734375) * x * y * y * y * y * y * y * z * z * z - std::sqrt(1837.5) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(64.599609375) * x * x * x * x * x * x * x * z + std::sqrt(64.599609375) * x * x * x * x * x * y * y * z - std::sqrt(6459.9609375) * x * x * x * x * x * z * z * z + std::sqrt(179.443359375) * x * x * x * y * y * y * y * z + std::sqrt(25839.84375) * x * x * x * y * y * z * z * z + std::sqrt(16537.5) * x * x * x * z * z * z * z * z + std::sqrt(179.443359375) * x * y * y * y * y * y * y * z - std::sqrt(35170.8984375) * x * y * y * y * y * z * z * z + std::sqrt(16537.5) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(1033.59375) * x * x * x * x * x * z + std::sqrt(103359.375) * x * x * x * y * y * z + std::sqrt(103359.375) * x * x * x * z * z * z - std::sqrt(25839.84375) * x * y * y * y * y * z + std::sqrt(103359.375) * x * y * y * z * z * z + std::sqrt(66150.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(103359.375) * x * x * x * z + std::sqrt(103359.375) * x * y * y * z + std::sqrt(1653750.0) * x * z * z * z) + e_4 * (std::sqrt(1653750.0) * x * z);

        pc_23[k] = e_0 * (-std::sqrt(0.5767822265625) * x * x * x * x * x * x * x * x * y * y - std::sqrt(4.1015625) * x * x * x * x * x * x * y * y * y * y + std::sqrt(230.712890625) * x * x * x * x * x * x * y * y * z * z - std::sqrt(2.30712890625) * x * x * x * x * y * y * y * y * y * y + std::sqrt(640.869140625) * x * x * x * x * y * y * y * y * z * z - std::sqrt(6238.4765625) * x * x * x * x * y * y * z * z * z * z + std::sqrt(25.634765625) * x * x * y * y * y * y * y * y * z * z - std::sqrt(2772.65625) * x * x * y * y * y * y * z * z * z * z + std::sqrt(2362.5) * x * x * y * y * z * z * z * z * z * z + std::sqrt(0.0640869140625) * y * y * y * y * y * y * y * y * y * y - std::sqrt(25.634765625) * y * y * y * y * y * y * y * y * z * z + std::sqrt(693.1640625) * y * y * y * y * y * y * z * z * z * z - std::sqrt(262.5) * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(0.5767822265625) * x * x * x * x * x * x * x * x - std::sqrt(332.2265625) * x * x * x * x * x * x * y * y + std::sqrt(230.712890625) * x * x * x * x * x * x * z * z - std::sqrt(775.45166015625) * x * x * x * x * y * y * y * y - std::sqrt(8868.603515625) * x * x * x * x * y * y * z * z - std::sqrt(6238.4765625) * x * x * x * x * z * z * z * z - std::sqrt(16.40625) * x * x * y * y * y * y * y * y - std::sqrt(7761.181640625) * x * x * y * y * y * y * z * z - std::sqrt(3691.40625) * x * x * y * y * z * z * z * z + std::sqrt(2362.5) * x * x * z * z * z * z * z * z + std::sqrt(40.0543212890625) * y * y * y * y * y * y * y * y + std::sqrt(452.197265625) * y * y * y * y * y * y * z * z + std::sqrt(9847.8515625) * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(147.65625) * x * x * x * x * x * x - std::sqrt(71465.625) * x * x * x * x * y * y - std::sqrt(28940.625) * x * x * x * x * z * z - std::sqrt(24953.90625) * x * x * y * y * y * y - std::sqrt(531562.5) * x * x * y * y * z * z + std::sqrt(9450.0) * x * x * z * z * z * z + std::sqrt(9450.0) * y * y * y * y * y * y + std::sqrt(170690.625) * y * y * y * y * z * z - std::sqrt(9450.0) * y * y * z * z * z * z) + e_3 * (-std::sqrt(42672.65625) * x * x * x * x - std::sqrt(1196015.625) * x * x * y * y - std::sqrt(151200.0) * x * x * z * z + std::sqrt(326172.65625) * y * y * y * y + std::sqrt(151200.0) * y * y * z * z) + e_4 * (-std::sqrt(463050.0) * x * x + std::sqrt(463050.0) * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_24[k] = e_0 * (-std::sqrt(8.6517333984375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(61.5234375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(984.375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(34.60693359375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(2734.375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(4136.8359375) * x * x * x * x * y * z * z * z * z * z + std::sqrt(109.375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(1838.59375) * x * x * y * y * y * z * z * z * z * z + std::sqrt(157.5) * x * x * y * z * z * z * z * z * z * z + std::sqrt(0.9613037109375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(109.375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(459.6484375) * y * y * y * y * y * z * z * z * z * z - std::sqrt(17.5) * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(138.427734375) * x * x * x * x * x * x * y * z - std::sqrt(384.521484375) * x * x * x * x * y * y * y * z - std::sqrt(246.09375) * x * x * x * x * y * z * z * z - std::sqrt(15.380859375) * x * x * y * y * y * y * y * z - std::sqrt(109.375) * x * x * y * y * y * z * z * z - std::sqrt(63000.0) * x * x * y * z * z * z * z * z + std::sqrt(15.380859375) * y * y * y * y * y * y * y * z + std::sqrt(27.34375) * y * y * y * y * y * z * z * z + std::sqrt(7000.0) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(19933.59375) * x * x * x * x * y * z - std::sqrt(8859.375) * x * x * y * y * y * z - std::sqrt(1736437.5) * x * x * y * z * z * z + std::sqrt(2214.84375) * y * y * y * y * y * z + std::sqrt(192937.5) * y * y * y * z * z * z) + e_3 * (-std::sqrt(2870437.5) * x * x * y * z + std::sqrt(318937.5) * y * y * y * z);

        pc_25[k] = e_0 * (-std::sqrt(0.5767822265625) * x * x * x * x * x * x * x * x * x * y - std::sqrt(4.1015625) * x * x * x * x * x * x * x * y * y * y + std::sqrt(230.712890625) * x * x * x * x * x * x * x * y * z * z - std::sqrt(2.30712890625) * x * x * x * x * x * y * y * y * y * y + std::sqrt(640.869140625) * x * x * x * x * x * y * y * y * z * z - std::sqrt(6238.4765625) * x * x * x * x * x * y * z * z * z * z + std::sqrt(25.634765625) * x * x * x * y * y * y * y * y * z * z - std::sqrt(2772.65625) * x * x * x * y * y * y * z * z * z * z + std::sqrt(2362.5) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.0640869140625) * x * y * y * y * y * y * y * y * y * y - std::sqrt(25.634765625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(693.1640625) * x * y * y * y * y * y * z * z * z * z - std::sqrt(262.5) * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(332.2265625) * x * x * x * x * x * x * x * y - std::sqrt(1050.0) * x * x * x * x * x * y * y * y - std::sqrt(6238.4765625) * x * x * x * x * x * y * z * z - std::sqrt(102.5390625) * x * x * x * y * y * y * y * y - std::sqrt(147.65625) * x * x * x * y * y * y * z * z - std::sqrt(47840.625) * x * x * x * y * z * z * z * z + std::sqrt(16.40625) * x * y * y * y * y * y * y * y + std::sqrt(4466.6015625) * x * y * y * y * y * y * z * z - std::sqrt(18965.625) * x * y * y * y * z * z * z * z + std::sqrt(9450.0) * x * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(78110.15625) * x * x * x * x * x * y - std::sqrt(47840.625) * x * x * x * y * y * y - std::sqrt(1143450.0) * x * x * x * y * z * z + std::sqrt(3691.40625) * x * y * y * y * y * y - std::sqrt(9450.0) * x * y * y * y * z * z + std::sqrt(37800.0) * x * y * z * z * z * z) + e_3 * (-std::sqrt(2270362.5) * x * x * x * y - std::sqrt(2362.5) * x * y * y * y - std::sqrt(604800.0) * x * y * z * z) + e_4 * (-std::sqrt(1852200.0) * x * y);

        pc_26[k] = e_0 * (std::sqrt(16.14990234375) * x * x * x * x * x * x * x * x * y * z + std::sqrt(7.177734375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(1614.990234375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(28.7109375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(179.443359375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(4134.375) * x * x * x * x * y * z * z * z * z * z - std::sqrt(7.177734375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1614.990234375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(7350.0) * x * x * y * y * y * z * z * z * z * z + std::sqrt(1.79443359375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(179.443359375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(459.375) * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(717.7734375) * x * x * x * x * y * y * y * z + std::sqrt(6459.9609375) * x * x * x * x * y * z * z * z + std::sqrt(1033.59375) * x * x * y * y * y * y * y * z - std::sqrt(140683.59375) * x * x * y * y * y * z * z * z + std::sqrt(16537.5) * x * x * y * z * z * z * z * z + std::sqrt(28.7109375) * y * y * y * y * y * y * y * z - std::sqrt(717.7734375) * y * y * y * y * y * z * z * z + std::sqrt(16537.5) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(25839.84375) * x * x * x * x * y * z - std::sqrt(103359.375) * x * x * y * y * y * z + std::sqrt(103359.375) * x * x * y * z * z * z + std::sqrt(1033.59375) * y * y * y * y * y * z + std::sqrt(103359.375) * y * y * y * z * z * z + std::sqrt(66150.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(103359.375) * x * x * y * z + std::sqrt(103359.375) * y * y * y * z + std::sqrt(1653750.0) * y * z * z * z) + e_4 * (std::sqrt(1653750.0) * y * z);

        pc_27[k] = e_0 * (0.8203125 * x * x * x * x * x * x * x * x * x * y - 1.09375 * x * x * x * x * x * x * x * y * y * y - 13.125 * x * x * x * x * x * x * x * y * z * z - 3.828125 * x * x * x * x * x * y * y * y * y * y + 30.625 * x * x * x * x * x * y * y * y * z * z + 52.5 * x * x * x * x * x * y * z * z * z * z - 1.09375 * x * x * x * y * y * y * y * y * y * y + 30.625 * x * x * x * y * y * y * y * y * z * z - 175.0 * x * x * x * y * y * y * z * z * z * z + 0.8203125 * x * y * y * y * y * y * y * y * y * y - 13.125 * x * y * y * y * y * y * y * y * z * z + 52.5 * x * y * y * y * y * y * z * z * z * z) + e_1 * (13.125 * x * x * x * x * x * x * x * y - 30.625 * x * x * x * x * x * y * y * y + 131.25 * x * x * x * x * x * y * z * z - 30.625 * x * x * x * y * y * y * y * y - 437.5 * x * x * x * y * y * y * z * z + 13.125 * x * y * y * y * y * y * y * y + 131.25 * x * y * y * y * y * y * z * z) + e_2 * (157.5 * x * x * x * x * x * y - 525.0 * x * x * x * y * y * y + 157.5 * x * y * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_28[k] = e_0 * (-std::sqrt(12.1124267578125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(344.53125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(775.1953125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(134.58251953125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(31093.9453125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(86.1328125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(6976.7578125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(1.3458251953125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(86.1328125) * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(1744.189453125) * x * x * x * x * x * x * y * z - std::sqrt(18109.423828125) * x * x * x * x * y * y * y * z - std::sqrt(12403.125) * x * x * x * x * y * z * z * z + std::sqrt(15697.705078125) * x * x * y * y * y * y * y * z - std::sqrt(49612.5) * x * x * y * y * y * z * z * z + std::sqrt(21.533203125) * y * y * y * y * y * y * y * z - std::sqrt(12403.125) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(3100.78125) * x * x * x * x * y * z - std::sqrt(12403.125) * x * x * y * y * y * z - std::sqrt(446512.5) * x * x * y * z * z * z - std::sqrt(3100.78125) * y * y * y * y * y * z - std::sqrt(446512.5) * y * y * y * z * z * z) + e_3 * (-std::sqrt(793800.0) * x * x * y * z - std::sqrt(793800.0) * y * y * y * z - std::sqrt(793800.0) * y * z * z * z) + e_4 * (-std::sqrt(2431012.5) * y * z);

        pc_29[k] = e_0 * (-std::sqrt(1.21124267578125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(105.5126953125) * x * x * x * x * x * x * x * y * y * y + std::sqrt(77.51953125) * x * x * x * x * x * x * x * y * z * z + std::sqrt(4.844970703125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(8277.36328125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(53.8330078125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(5383.30078125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(3.36456298828125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(215.33203125) * x * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(21533.203125) * x * x * x * x * x * y * y * y - std::sqrt(7751.953125) * x * x * x * x * x * y * z * z - std::sqrt(3445.3125) * x * x * x * y * y * y * y * y - std::sqrt(31007.8125) * x * x * x * y * y * y * z * z + std::sqrt(861.328125) * x * y * y * y * y * y * y * y - std::sqrt(7751.953125) * x * y * y * y * y * y * z * z) + e_2 * (std::sqrt(31007.8125) * x * x * x * x * x * y + std::sqrt(124031.25) * x * x * x * y * y * y - std::sqrt(496125.0) * x * x * x * y * z * z + std::sqrt(31007.8125) * x * y * y * y * y * y - std::sqrt(496125.0) * x * y * y * y * z * z) + e_3 * (std::sqrt(496125.0) * x * x * x * y + std::sqrt(496125.0) * x * y * y * y - std::sqrt(1984500.0) * x * y * z * z) + e_4 * (std::sqrt(496125.0) * x * y);

        pc_30[k] = e_0 * (26.25 * x * x * x * x * x * x * y * y * z * z + 52.5 * x * x * x * x * y * y * y * y * z * z - 105.0 * x * x * x * x * y * y * z * z * z * z + 26.25 * x * x * y * y * y * y * y * y * z * z - 105.0 * x * x * y * y * y * y * z * z * z * z + 105.0 * x * x * y * y * z * z * z * z * z * z) + e_1 * (26.25 * x * x * x * x * x * x * y * y + 26.25 * x * x * x * x * x * x * z * z + 52.5 * x * x * x * x * y * y * y * y + 78.75 * x * x * x * x * y * y * z * z - 105.0 * x * x * x * x * z * z * z * z + 26.25 * x * x * y * y * y * y * y * y + 78.75 * x * x * y * y * y * y * z * z + 315.0 * x * x * y * y * z * z * z * z + 105.0 * x * x * z * z * z * z * z * z + 26.25 * y * y * y * y * y * y * z * z - 105.0 * y * y * y * y * z * z * z * z + 105.0 * y * y * z * z * z * z * z * z) + e_2 * (26.25 * x * x * x * x * x * x + 393.75 * x * x * x * x * y * y - 78.75 * x * x * x * x * z * z + 393.75 * x * x * y * y * y * y + 1417.5 * x * x * y * y * z * z + 630.0 * x * x * z * z * z * z + 26.25 * y * y * y * y * y * y - 78.75 * y * y * y * y * z * z + 630.0 * y * y * z * z * z * z + 105.0 * z * z * z * z * z * z) + e_3 * (236.25 * x * x * x * x + 2047.5 * x * x * y * y + 1575.0 * x * x * z * z + 236.25 * y * y * y * y + 1575.0 * y * y * z * z + 945.0 * z * z * z * z) + e_4 * (1260.0 * x * x + 1260.0 * y * y + 2205.0 * z * z) + e_5 * (945.0);

        pc_31[k] = e_0 * (-std::sqrt(6.15234375) * x * x * x * x * x * x * x * y * y * z - std::sqrt(55.37109375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1205.859375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(55.37109375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(4823.4375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(6300.0) * x * x * x * y * y * z * z * z * z * z - std::sqrt(6.15234375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1205.859375) * x * y * y * y * y * y * y * z * z * z - std::sqrt(6300.0) * x * y * y * y * y * z * z * z * z * z + std::sqrt(1575.0) * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(6.15234375) * x * x * x * x * x * x * x * z + std::sqrt(55.37109375) * x * x * x * x * x * y * y * z + std::sqrt(1205.859375) * x * x * x * x * x * z * z * z + std::sqrt(498.33984375) * x * x * x * y * y * y * y * z - std::sqrt(885.9375) * x * x * x * y * y * z * z * z - std::sqrt(6300.0) * x * x * x * z * z * z * z * z + std::sqrt(153.80859375) * x * y * y * y * y * y * y * z - std::sqrt(4158.984375) * x * y * y * y * y * z * z * z + std::sqrt(14175.0) * x * y * y * z * z * z * z * z + std::sqrt(1575.0) * x * z * z * z * z * z * z * z) + e_2 * (std::sqrt(885.9375) * x * x * x * x * x * z + std::sqrt(3543.75) * x * x * x * y * y * z - std::sqrt(56700.0) * x * x * x * z * z * z + std::sqrt(885.9375) * x * y * y * y * y * z + std::sqrt(127575.0) * x * y * y * z * z * z + std::sqrt(127575.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(14175.0) * x * x * x * z + std::sqrt(226800.0) * x * y * y * z + std::sqrt(1148175.0) * x * z * z * z) + e_4 * (std::sqrt(694575.0) * x * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_32[k] = e_0 * (-std::sqrt(92.28515625) * x * x * x * x * x * x * x * y * z * z - std::sqrt(830.56640625) * x * x * x * x * x * y * y * y * z * z + std::sqrt(2009.765625) * x * x * x * x * x * y * z * z * z * z - std::sqrt(830.56640625) * x * x * x * y * y * y * y * y * z * z + std::sqrt(8039.0625) * x * x * x * y * y * y * z * z * z * z - std::sqrt(3176.25) * x * x * x * y * z * z * z * z * z * z - std::sqrt(92.28515625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(2009.765625) * x * y * y * y * y * y * z * z * z * z - std::sqrt(3176.25) * x * y * y * y * z * z * z * z * z * z + std::sqrt(105.0) * x * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(92.28515625) * x * x * x * x * x * x * x * y - std::sqrt(830.56640625) * x * x * x * x * x * y * y * y - std::sqrt(369.140625) * x * x * x * x * x * y * z * z - std::sqrt(830.56640625) * x * x * x * y * y * y * y * y - std::sqrt(1476.5625) * x * x * x * y * y * y * z * z - std::sqrt(16406.25) * x * x * x * y * z * z * z * z - std::sqrt(92.28515625) * x * y * y * y * y * y * y * y - std::sqrt(369.140625) * x * y * y * y * y * y * z * z - std::sqrt(16406.25) * x * y * y * y * z * z * z * z - std::sqrt(2625.0) * x * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(23625.0) * x * x * x * x * x * y - std::sqrt(94500.0) * x * x * x * y * y * y - std::sqrt(289406.25) * x * x * x * y * z * z - std::sqrt(23625.0) * x * y * y * y * y * y - std::sqrt(289406.25) * x * y * y * y * z * z - std::sqrt(590625.0) * x * y * z * z * z * z) + e_3 * (-std::sqrt(998156.25) * x * x * x * y - std::sqrt(998156.25) * x * y * y * y - std::sqrt(6827625.0) * x * y * z * z) + e_4 * (-std::sqrt(4630500.0) * x * y);

        pc_33[k] = e_0 * (-std::sqrt(6.15234375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(55.37109375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1205.859375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(55.37109375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(4823.4375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(6300.0) * x * x * x * x * y * z * z * z * z * z - std::sqrt(6.15234375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1205.859375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(6300.0) * x * x * y * y * y * z * z * z * z * z + std::sqrt(1575.0) * x * x * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(153.80859375) * x * x * x * x * x * x * y * z + std::sqrt(498.33984375) * x * x * x * x * y * y * y * z - std::sqrt(4158.984375) * x * x * x * x * y * z * z * z + std::sqrt(55.37109375) * x * x * y * y * y * y * y * z - std::sqrt(885.9375) * x * x * y * y * y * z * z * z + std::sqrt(14175.0) * x * x * y * z * z * z * z * z - std::sqrt(6.15234375) * y * y * y * y * y * y * y * z + std::sqrt(1205.859375) * y * y * y * y * y * z * z * z - std::sqrt(6300.0) * y * y * y * z * z * z * z * z + std::sqrt(1575.0) * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(885.9375) * x * x * x * x * y * z + std::sqrt(3543.75) * x * x * y * y * y * z + std::sqrt(127575.0) * x * x * y * z * z * z + std::sqrt(885.9375) * y * y * y * y * y * z - std::sqrt(56700.0) * y * y * y * z * z * z + std::sqrt(127575.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(226800.0) * x * x * y * z - std::sqrt(14175.0) * y * y * y * z + std::sqrt(1148175.0) * y * z * z * z) + e_4 * (std::sqrt(694575.0) * y * z);

        pc_34[k] = e_0 * (13.125 * x * x * x * x * x * x * x * y * z * z + 13.125 * x * x * x * x * x * y * y * y * z * z - 52.5 * x * x * x * x * x * y * z * z * z * z - 13.125 * x * x * x * y * y * y * y * y * z * z + 52.5 * x * x * x * y * z * z * z * z * z * z - 13.125 * x * y * y * y * y * y * y * y * z * z + 52.5 * x * y * y * y * y * y * z * z * z * z - 52.5 * x * y * y * y * z * z * z * z * z * z) + e_1 * (13.125 * x * x * x * x * x * x * x * y + 13.125 * x * x * x * x * x * y * y * y - 13.125 * x * x * x * y * y * y * y * y + 262.5 * x * x * x * y * z * z * z * z - 13.125 * x * y * y * y * y * y * y * y - 262.5 * x * y * y * y * z * z * z * z) + e_2 * (157.5 * x * x * x * x * x * y + 787.5 * x * x * x * y * z * z - 157.5 * x * y * y * y * y * y - 787.5 * x * y * y * y * z * z) + e_3 * (787.5 * x * x * x * y - 787.5 * x * y * y * y);

        pc_35[k] = e_0 * (std::sqrt(7.177734375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(7.177734375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(717.7734375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(179.443359375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(2871.09375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(1837.5) * x * x * x * x * y * z * z * z * z * z - std::sqrt(64.599609375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(6459.9609375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(16537.5) * x * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(179.443359375) * x * x * x * x * x * x * y * z - std::sqrt(179.443359375) * x * x * x * x * y * y * y * z + std::sqrt(35170.8984375) * x * x * x * x * y * z * z * z - std::sqrt(64.599609375) * x * x * y * y * y * y * y * z - std::sqrt(25839.84375) * x * x * y * y * y * z * z * z - std::sqrt(16537.5) * x * x * y * z * z * z * z * z - std::sqrt(64.599609375) * y * y * y * y * y * y * y * z + std::sqrt(6459.9609375) * y * y * y * y * y * z * z * z - std::sqrt(16537.5) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(25839.84375) * x * x * x * x * y * z - std::sqrt(103359.375) * x * x * y * y * y * z - std::sqrt(103359.375) * x * x * y * z * z * z + std::sqrt(1033.59375) * y * y * y * y * y * z - std::sqrt(103359.375) * y * y * y * z * z * z - std::sqrt(66150.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(103359.375) * x * x * y * z - std::sqrt(103359.375) * y * y * y * z - std::sqrt(1653750.0) * y * z * z * z) + e_4 * (-std::sqrt(1653750.0) * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_36[k] = e_0 * (-std::sqrt(129.19921875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(3229.98046875) * x * x * x * x * x * y * y * y * z * z + std::sqrt(516.796875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(3229.98046875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(18604.6875) * x * x * x * y * y * y * z * z * z * z - std::sqrt(129.19921875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(516.796875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(129.19921875) * x * x * x * x * x * x * x * y + std::sqrt(3229.98046875) * x * x * x * x * x * y * y * y + std::sqrt(4651.171875) * x * x * x * x * x * y * z * z + std::sqrt(3229.98046875) * x * x * x * y * y * y * y * y + std::sqrt(101292.1875) * x * x * x * y * y * y * z * z - std::sqrt(33075.0) * x * x * x * y * z * z * z * z - std::sqrt(129.19921875) * x * y * y * y * y * y * y * y + std::sqrt(4651.171875) * x * y * y * y * y * y * z * z - std::sqrt(33075.0) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(529200.0) * x * x * x * y * y * y + std::sqrt(74418.75) * x * x * x * y * z * z + std::sqrt(74418.75) * x * y * y * y * z * z - std::sqrt(297675.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(669768.75) * x * x * x * y + std::sqrt(669768.75) * x * y * y * y - std::sqrt(297675.0) * x * y * z * z) + e_4 * (std::sqrt(1190700.0) * x * y);

        pc_37[k] = e_0 * (-std::sqrt(12.919921875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(1046.513671875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(51.6796875) * x * x * x * x * x * x * y * z * z * z + std::sqrt(322.998046875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(5167.96875) * x * x * x * x * y * y * y * z * z * z - std::sqrt(322.998046875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1291.9921875) * x * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(322.998046875) * x * x * x * x * x * x * y * z + std::sqrt(201873.779296875) * x * x * x * x * y * y * y * z - std::sqrt(11627.9296875) * x * x * x * x * y * z * z * z - std::sqrt(26162.841796875) * x * x * y * y * y * y * y * z - std::sqrt(5167.96875) * x * x * y * y * y * z * z * z - std::sqrt(322.998046875) * y * y * y * y * y * y * y * z + std::sqrt(1291.9921875) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(418605.46875) * x * x * x * x * y * z + std::sqrt(186046.875) * x * x * y * y * y * z - std::sqrt(186046.875) * x * x * y * z * z * z - std::sqrt(46511.71875) * y * y * y * y * y * z + std::sqrt(20671.875) * y * y * y * z * z * z) + e_3 * (std::sqrt(1674421.875) * x * x * y * z - std::sqrt(186046.875) * y * y * y * z);

        pc_38[k] = e_0 * (0.234375 * x * x * x * x * x * x * x * x * y * y + 0.9375 * x * x * x * x * x * x * y * y * y * y - 5.625 * x * x * x * x * x * x * y * y * z * z + 1.40625 * x * x * x * x * y * y * y * y * y * y - 16.875 * x * x * x * x * y * y * y * y * z * z + 37.5 * x * x * x * x * y * y * z * z * z * z + 0.9375 * x * x * y * y * y * y * y * y * y * y - 16.875 * x * x * y * y * y * y * y * y * z * z + 75.0 * x * x * y * y * y * y * z * z * z * z - 45.0 * x * x * y * y * z * z * z * z * z * z + 0.234375 * y * y * y * y * y * y * y * y * y * y - 5.625 * y * y * y * y * y * y * y * y * z * z + 37.5 * y * y * y * y * y * y * z * z * z * z - 45.0 * y * y * y * y * z * z * z * z * z * z + 15.0 * y * y * z * z * z * z * z * z * z * z) + e_1 * (0.234375 * x * x * x * x * x * x * x * x + 6.5625 * x * x * x * x * x * x * y * y - 5.625 * x * x * x * x * x * x * z * z + 18.28125 * x * x * x * x * y * y * y * y + 39.375 * x * x * x * x * y * y * z * z + 37.5 * x * x * x * x * z * z * z * z + 17.8125 * x * x * y * y * y * y * y * y + 95.625 * x * x * y * y * y * y * z * z - 45.0 * x * x * z * z * z * z * z * z + 5.859375 * y * y * y * y * y * y * y * y + 50.625 * y * y * y * y * y * y * z * z - 37.5 * y * y * y * y * z * z * z * z + 105.0 * y * y * z * z * z * z * z * z + 15.0 * z * z * z * z * z * z * z * z) + e_2 * (3.75 * x * x * x * x * x * x + 123.75 * x * x * x * x * y * y + 90.0 * x * x * x * x * z * z + 236.25 * x * x * y * y * y * y + 405.0 * x * x * y * y * z * z - 225.0 * x * x * z * z * z * z + 116.25 * y * y * y * y * y * y + 315.0 * y * y * y * y * z * z + 675.0 * y * y * z * z * z * z + 240.0 * z * z * z * z * z * z) + e_3 * (90.0 * x * x * x * x + 855.0 * x * x * y * y - 135.0 * x * x * z * z + 765.0 * y * y * y * y + 2115.0 * y * y * z * z + 1350.0 * z * z * z * z) + e_4 * (315.0 * x * x + 1890.0 * y * y + 2520.0 * z * z) + e_5 * (945.0);

        pc_39[k] = e_0 * (std::sqrt(0.823974609375) * x * x * x * x * x * x * x * x * y * z + std::sqrt(13.18359375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(177.24609375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(29.6630859375) * x * x * x * x * y * y * y * y * y * z - std::sqrt(1595.21484375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(1353.75) * x * x * x * x * y * z * z * z * z * z + std::sqrt(13.18359375) * x * x * y * y * y * y * y * y * y * z - std::sqrt(1595.21484375) * x * x * y * y * y * y * y * z * z * z + std::sqrt(5415.0) * x * x * y * y * y * z * z * z * z * z - std::sqrt(633.75) * x * x * y * z * z * z * z * z * z * z + std::sqrt(0.823974609375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(177.24609375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(1353.75) * y * y * y * y * y * z * z * z * z * z - std::sqrt(633.75) * y * y * y * z * z * z * z * z * z * z + std::sqrt(15.0) * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(13.18359375) * x * x * x * x * x * x * y * z - std::sqrt(118.65234375) * x * x * x * x * y * y * y * z + std::sqrt(2343.75) * x * x * x * x * y * z * z * z - std::sqrt(118.65234375) * x * x * y * y * y * y * y * z + std::sqrt(9375.0) * x * x * y * y * y * z * z * z - std::sqrt(7593.75) * x * x * y * z * z * z * z * z - std::sqrt(13.18359375) * y * y * y * y * y * y * y * z + std::sqrt(2343.75) * y * y * y * y * y * z * z * z - std::sqrt(7593.75) * y * y * y * z * z * z * z * z + std::sqrt(1500.0) * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(843.75) * x * x * x * x * y * z + std::sqrt(3375.0) * x * x * y * y * y * z - std::sqrt(21093.75) * x * x * y * z * z * z + std::sqrt(843.75) * y * y * y * y * y * z - std::sqrt(21093.75) * y * y * y * z * z * z + std::sqrt(54000.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(843.75) * x * x * y * z - std::sqrt(843.75) * y * y * y * z + std::sqrt(337500.0) * y * z * z * z) + e_4 * (std::sqrt(165375.0) * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_40[k] = e_0 * (0.234375 * x * x * x * x * x * x * x * x * x * y + 0.9375 * x * x * x * x * x * x * x * y * y * y - 5.625 * x * x * x * x * x * x * x * y * z * z + 1.40625 * x * x * x * x * x * y * y * y * y * y - 16.875 * x * x * x * x * x * y * y * y * z * z + 37.5 * x * x * x * x * x * y * z * z * z * z + 0.9375 * x * x * x * y * y * y * y * y * y * y - 16.875 * x * x * x * y * y * y * y * y * z * z + 75.0 * x * x * x * y * y * y * z * z * z * z - 45.0 * x * x * x * y * z * z * z * z * z * z + 0.234375 * x * y * y * y * y * y * y * y * y * y - 5.625 * x * y * y * y * y * y * y * y * z * z + 37.5 * x * y * y * y * y * y * z * z * z * z - 45.0 * x * y * y * y * z * z * z * z * z * z + 15.0 * x * y * z * z * z * z * z * z * z * z) + e_1 * (5.625 * x * x * x * x * x * x * x * y + 16.875 * x * x * x * x * x * y * y * y + 56.25 * x * x * x * x * x * y * z * z + 16.875 * x * x * x * y * y * y * y * y + 112.5 * x * x * x * y * y * y * z * z - 75.0 * x * x * x * y * z * z * z * z + 5.625 * x * y * y * y * y * y * y * y + 56.25 * x * y * y * y * y * y * z * z - 75.0 * x * y * y * y * z * z * z * z + 150.0 * x * y * z * z * z * z * z * z) + e_2 * (112.5 * x * x * x * x * x * y + 225.0 * x * x * x * y * y * y + 225.0 * x * x * x * y * z * z + 112.5 * x * y * y * y * y * y + 225.0 * x * y * y * y * z * z + 900.0 * x * y * z * z * z * z) + e_3 * (675.0 * x * x * x * y + 675.0 * x * y * y * y + 2250.0 * x * y * z * z) + e_4 * (1575.0 * x * y);

        pc_41[k] = e_0 * (-std::sqrt(1.5380859375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(6.15234375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(301.46484375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(301.46484375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(1575.0) * x * x * x * x * y * z * z * z * z * z + std::sqrt(6.15234375) * x * x * y * y * y * y * y * y * y * z - std::sqrt(301.46484375) * x * x * y * y * y * y * y * z * z * z + std::sqrt(393.75) * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.5380859375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(301.46484375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(1575.0) * y * y * y * y * y * z * z * z * z * z - std::sqrt(393.75) * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(98.4375) * x * x * x * x * x * x * y * z + std::sqrt(221.484375) * x * x * x * x * y * y * y * z - std::sqrt(7112.109375) * x * x * x * x * y * z * z * z - std::sqrt(4823.4375) * x * x * y * y * y * z * z * z + std::sqrt(31893.75) * x * x * y * z * z * z * z * z - std::sqrt(24.609375) * y * y * y * y * y * y * y * z + std::sqrt(221.484375) * y * y * y * y * y * z * z * z - std::sqrt(393.75) * y * y * y * z * z * z * z * z - std::sqrt(1575.0) * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(885.9375) * x * x * x * x * y * z - std::sqrt(3543.75) * x * x * y * y * y * z + std::sqrt(287043.75) * x * x * y * z * z * z - std::sqrt(885.9375) * y * y * y * y * y * z - std::sqrt(3543.75) * y * y * y * z * z * z - std::sqrt(127575.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(173643.75) * x * x * y * z - std::sqrt(31893.75) * y * y * y * z - std::sqrt(1148175.0) * y * z * z * z) + e_4 * (-std::sqrt(694575.0) * y * z);

        pc_42[k] = e_0 * (-std::sqrt(0.0640869140625) * x * x * x * x * x * x * x * x * x * y + std::sqrt(25.634765625) * x * x * x * x * x * x * x * y * z * z + std::sqrt(2.30712890625) * x * x * x * x * x * y * y * y * y * y - std::sqrt(25.634765625) * x * x * x * x * x * y * y * y * z * z - std::sqrt(693.1640625) * x * x * x * x * x * y * z * z * z * z + std::sqrt(4.1015625) * x * x * x * y * y * y * y * y * y * y - std::sqrt(640.869140625) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2772.65625) * x * x * x * y * y * y * z * z * z * z + std::sqrt(262.5) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.5767822265625) * x * y * y * y * y * y * y * y * y * y - std::sqrt(230.712890625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(6238.4765625) * x * y * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(16.40625) * x * x * x * x * x * x * x * y + std::sqrt(102.5390625) * x * x * x * x * x * y * y * y - std::sqrt(4466.6015625) * x * x * x * x * x * y * z * z + std::sqrt(1050.0) * x * x * x * y * y * y * y * y + std::sqrt(147.65625) * x * x * x * y * y * y * z * z + std::sqrt(18965.625) * x * x * x * y * z * z * z * z + std::sqrt(332.2265625) * x * y * y * y * y * y * y * y + std::sqrt(6238.4765625) * x * y * y * y * y * y * z * z + std::sqrt(47840.625) * x * y * y * y * z * z * z * z - std::sqrt(9450.0) * x * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(3691.40625) * x * x * x * x * x * y + std::sqrt(47840.625) * x * x * x * y * y * y + std::sqrt(9450.0) * x * x * x * y * z * z + std::sqrt(78110.15625) * x * y * y * y * y * y + std::sqrt(1143450.0) * x * y * y * y * z * z - std::sqrt(37800.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(2362.5) * x * x * x * y + std::sqrt(2270362.5) * x * y * y * y + std::sqrt(604800.0) * x * y * z * z) + e_4 * (std::sqrt(1852200.0) * x * y);

        pc_43[k] = e_0 * (std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(18.45703125) * x * x * x * x * x * x * y * y * y * z - std::sqrt(166.11328125) * x * x * x * x * x * x * y * z * z * z - std::sqrt(115.3564453125) * x * x * x * x * y * y * y * y * y * z + std::sqrt(4152.83203125) * x * x * x * x * y * y * y * z * z * z + std::sqrt(73.828125) * x * x * x * x * y * z * z * z * z * z - std::sqrt(18.45703125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(4152.83203125) * x * x * y * y * y * y * y * z * z * z - std::sqrt(2657.8125) * x * x * y * y * y * z * z * z * z * z + std::sqrt(1.153564453125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(166.11328125) * y * y * y * y * y * y * y * z * z * z + std::sqrt(73.828125) * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(461.42578125) * x * x * x * x * x * x * y * z + std::sqrt(461.42578125) * x * x * x * x * y * y * y * z + std::sqrt(7382.8125) * x * x * x * x * y * z * z * z + std::sqrt(1495.01953125) * x * x * y * y * y * y * y * z + std::sqrt(265781.25) * x * x * y * y * y * z * z * z - std::sqrt(10631.25) * x * x * y * z * z * z * z * z - std::sqrt(18.45703125) * y * y * y * y * y * y * y * z - std::sqrt(14470.3125) * y * y * y * y * y * z * z * z + std::sqrt(1181.25) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(1063125.0) * x * x * y * y * y * z + std::sqrt(265781.25) * x * x * y * z * z * z - std::sqrt(42525.0) * y * y * y * y * y * z - std::sqrt(29531.25) * y * y * y * z * z * z) + e_3 * (std::sqrt(2392031.25) * x * x * y * z - std::sqrt(265781.25) * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_44[k] = e_0 * (std::sqrt(0.1153564453125) * x * x * x * x * x * x * x * x * x * y - std::sqrt(7.3828125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(16.611328125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(22.60986328125) * x * x * x * x * x * y * y * y * y * y + std::sqrt(1345.517578125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(7.3828125) * x * x * x * x * x * y * z * z * z * z + std::sqrt(415.283203125) * x * x * x * y * y * y * y * y * z * z - std::sqrt(738.28125) * x * x * x * y * y * y * z * z * z * z + std::sqrt(2.8839111328125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(415.283203125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(184.5703125) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(4614.2578125) * x * x * x * x * x * y * y * y + std::sqrt(1661.1328125) * x * x * x * x * x * y * z * z - std::sqrt(738.28125) * x * x * x * y * y * y * y * y + std::sqrt(166113.28125) * x * x * x * y * y * y * z * z - std::sqrt(2953.125) * x * x * x * y * z * z * z * z + std::sqrt(1661.1328125) * x * y * y * y * y * y * y * y - std::sqrt(81395.5078125) * x * y * y * y * y * y * z * z + std::sqrt(2953.125) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(6644.53125) * x * x * x * x * x * y - std::sqrt(73828.125) * x * x * x * y * y * y + std::sqrt(425250.0) * x * x * x * y * z * z + std::sqrt(59800.78125) * x * y * y * y * y * y - std::sqrt(425250.0) * x * y * y * y * z * z) + e_3 * (-std::sqrt(106312.5) * x * x * x * y + std::sqrt(106312.5) * x * y * y * y);

        pc_45[k] = e_0 * (3.515625 * x * x * x * x * x * x * x * x * z * z + 14.0625 * x * x * x * x * x * x * y * y * z * z - 18.75 * x * x * x * x * x * x * z * z * z * z + 21.09375 * x * x * x * x * y * y * y * y * z * z - 56.25 * x * x * x * x * y * y * z * z * z * z + 28.75 * x * x * x * x * z * z * z * z * z * z + 14.0625 * x * x * y * y * y * y * y * y * z * z - 56.25 * x * x * y * y * y * y * z * z * z * z + 57.5 * x * x * y * y * z * z * z * z * z * z - 10.0 * x * x * z * z * z * z * z * z * z * z + 3.515625 * y * y * y * y * y * y * y * y * z * z - 18.75 * y * y * y * y * y * y * z * z * z * z + 28.75 * y * y * y * y * z * z * z * z * z * z - 10.0 * y * y * z * z * z * z * z * z * z * z + z * z * z * z * z * z * z * z * z * z) + e_1 * (3.515625 * x * x * x * x * x * x * x * x + 14.0625 * x * x * x * x * x * x * y * y + 21.09375 * x * x * x * x * y * y * y * y + 93.75 * x * x * x * x * z * z * z * z + 14.0625 * x * x * y * y * y * y * y * y + 187.5 * x * x * y * y * z * z * z * z - 50.0 * x * x * z * z * z * z * z * z + 3.515625 * y * y * y * y * y * y * y * y + 93.75 * y * y * y * y * z * z * z * z - 50.0 * y * y * z * z * z * z * z * z + 25.0 * z * z * z * z * z * z * z * z) + e_2 * (56.25 * x * x * x * x * x * x + 168.75 * x * x * x * x * y * y + 281.25 * x * x * x * x * z * z + 168.75 * x * x * y * y * y * y + 562.5 * x * x * y * y * z * z + 56.25 * y * y * y * y * y * y + 281.25 * y * y * y * y * z * z + 300.0 * z * z * z * z * z * z) + e_3 * (431.25 * x * x * x * x + 862.5 * x * x * y * y + 750.0 * x * x * z * z + 431.25 * y * y * y * y + 750.0 * y * y * z * z + 1500.0 * z * z * z * z) + e_4 * (1050.0 * x * x + 1050.0 * y * y + 2625.0 * z * z) + e_5 * (945.0);

        pc_46[k] = e_0 * (std::sqrt(0.823974609375) * x * x * x * x * x * x * x * x * x * z + std::sqrt(13.18359375) * x * x * x * x * x * x * x * y * y * z - std::sqrt(177.24609375) * x * x * x * x * x * x * x * z * z * z + std::sqrt(29.6630859375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(1595.21484375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(1353.75) * x * x * x * x * x * z * z * z * z * z + std::sqrt(13.18359375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(1595.21484375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(5415.0) * x * x * x * y * y * z * z * z * z * z - std::sqrt(633.75) * x * x * x * z * z * z * z * z * z * z + std::sqrt(0.823974609375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(177.24609375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(1353.75) * x * y * y * y * y * z * z * z * z * z - std::sqrt(633.75) * x * y * y * z * z * z * z * z * z * z + std::sqrt(15.0) * x * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(13.18359375) * x * x * x * x * x * x * x * z - std::sqrt(118.65234375) * x * x * x * x * x * y * y * z + std::sqrt(2343.75) * x * x * x * x * x * z * z * z - std::sqrt(118.65234375) * x * x * x * y * y * y * y * z + std::sqrt(9375.0) * x * x * x * y * y * z * z * z - std::sqrt(7593.75) * x * x * x * z * z * z * z * z - std::sqrt(13.18359375) * x * y * y * y * y * y * y * z + std::sqrt(2343.75) * x * y * y * y * y * z * z * z - std::sqrt(7593.75) * x * y * y * z * z * z * z * z + std::sqrt(1500.0) * x * z * z * z * z * z * z * z) + e_2 * (std::sqrt(843.75) * x * x * x * x * x * z + std::sqrt(3375.0) * x * x * x * y * y * z - std::sqrt(21093.75) * x * x * x * z * z * z + std::sqrt(843.75) * x * y * y * y * y * z - std::sqrt(21093.75) * x * y * y * z * z * z + std::sqrt(54000.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(843.75) * x * x * x * z - std::sqrt(843.75) * x * y * y * z + std::sqrt(337500.0) * x * z * z * z) + e_4 * (std::sqrt(165375.0) * x * z);

        pc_47[k] = e_0 * (-std::sqrt(23.0712890625) * x * x * x * x * x * x * x * x * z * z - std::sqrt(92.28515625) * x * x * x * x * x * x * y * y * z * z + std::sqrt(502.44140625) * x * x * x * x * x * x * z * z * z * z + std::sqrt(502.44140625) * x * x * x * x * y * y * z * z * z * z - std::sqrt(794.0625) * x * x * x * x * z * z * z * z * z * z + std::sqrt(92.28515625) * x * x * y * y * y * y * y * y * z * z - std::sqrt(502.44140625) * x * x * y * y * y * y * z * z * z * z + std::sqrt(26.25) * x * x * z * z * z * z * z * z * z * z + std::sqrt(23.0712890625) * y * y * y * y * y * y * y * y * z * z - std::sqrt(502.44140625) * y * y * y * y * y * y * z * z * z * z + std::sqrt(794.0625) * y * y * y * y * z * z * z * z * z * z - std::sqrt(26.25) * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(23.0712890625) * x * x * x * x * x * x * x * x - std::sqrt(92.28515625) * x * x * x * x * x * x * y * y - std::sqrt(92.28515625) * x * x * x * x * x * x * z * z - std::sqrt(92.28515625) * x * x * x * x * y * y * z * z - std::sqrt(4101.5625) * x * x * x * x * z * z * z * z + std::sqrt(92.28515625) * x * x * y * y * y * y * y * y + std::sqrt(92.28515625) * x * x * y * y * y * y * z * z - std::sqrt(656.25) * x * x * z * z * z * z * z * z + std::sqrt(23.0712890625) * y * y * y * y * y * y * y * y + std::sqrt(92.28515625) * y * y * y * y * y * y * z * z + std::sqrt(4101.5625) * y * y * y * y * z * z * z * z + std::sqrt(656.25) * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(5906.25) * x * x * x * x * x * x - std::sqrt(5906.25) * x * x * x * x * y * y - std::sqrt(72351.5625) * x * x * x * x * z * z + std::sqrt(5906.25) * x * x * y * y * y * y - std::sqrt(147656.25) * x * x * z * z * z * z + std::sqrt(5906.25) * y * y * y * y * y * y + std::sqrt(72351.5625) * y * y * y * y * z * z + std::sqrt(147656.25) * y * y * z * z * z * z) + e_3 * (-std::sqrt(249539.0625) * x * x * x * x - std::sqrt(1706906.25) * x * x * z * z + std::sqrt(249539.0625) * y * y * y * y + std::sqrt(1706906.25) * y * y * z * z) + e_4 * (-std::sqrt(1157625.0) * x * x + std::sqrt(1157625.0) * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_48[k] = e_0 * (-std::sqrt(0.9613037109375) * x * x * x * x * x * x * x * x * x * z + std::sqrt(109.375) * x * x * x * x * x * x * x * z * z * z + std::sqrt(34.60693359375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(109.375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(459.6484375) * x * x * x * x * x * z * z * z * z * z + std::sqrt(61.5234375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(2734.375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(1838.59375) * x * x * x * y * y * z * z * z * z * z + std::sqrt(17.5) * x * x * x * z * z * z * z * z * z * z + std::sqrt(8.6517333984375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(984.375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(4136.8359375) * x * y * y * y * y * z * z * z * z * z - std::sqrt(157.5) * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(15.380859375) * x * x * x * x * x * x * x * z + std::sqrt(15.380859375) * x * x * x * x * x * y * y * z - std::sqrt(27.34375) * x * x * x * x * x * z * z * z + std::sqrt(384.521484375) * x * x * x * y * y * y * y * z + std::sqrt(109.375) * x * x * x * y * y * z * z * z - std::sqrt(7000.0) * x * x * x * z * z * z * z * z + std::sqrt(138.427734375) * x * y * y * y * y * y * y * z + std::sqrt(246.09375) * x * y * y * y * y * z * z * z + std::sqrt(63000.0) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x * z + std::sqrt(8859.375) * x * x * x * y * y * z - std::sqrt(192937.5) * x * x * x * z * z * z + std::sqrt(19933.59375) * x * y * y * y * y * z + std::sqrt(1736437.5) * x * y * y * z * z * z) + e_3 * (-std::sqrt(318937.5) * x * x * x * z + std::sqrt(2870437.5) * x * y * y * z);

        pc_49[k] = e_0 * (std::sqrt(17.303466796875) * x * x * x * x * x * x * x * x * z * z - std::sqrt(276.85546875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(123.046875) * x * x * x * x * x * x * z * z * z * z - std::sqrt(1730.3466796875) * x * x * x * x * y * y * y * y * z * z + std::sqrt(3076.171875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(4.921875) * x * x * x * x * z * z * z * z * z * z - std::sqrt(276.85546875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(3076.171875) * x * x * y * y * y * y * z * z * z * z - std::sqrt(177.1875) * x * x * y * y * z * z * z * z * z * z + std::sqrt(17.303466796875) * y * y * y * y * y * y * y * y * z * z - std::sqrt(123.046875) * y * y * y * y * y * y * z * z * z * z + std::sqrt(4.921875) * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(17.303466796875) * x * x * x * x * x * x * x * x - std::sqrt(276.85546875) * x * x * x * x * x * x * y * y + std::sqrt(1107.421875) * x * x * x * x * x * x * z * z - std::sqrt(1730.3466796875) * x * x * x * x * y * y * y * y - std::sqrt(27685.546875) * x * x * x * x * y * y * z * z - std::sqrt(6029.296875) * x * x * x * x * z * z * z * z - std::sqrt(276.85546875) * x * x * y * y * y * y * y * y - std::sqrt(27685.546875) * x * x * y * y * y * y * z * z + std::sqrt(217054.6875) * x * x * y * y * z * z * z * z + std::sqrt(17.303466796875) * y * y * y * y * y * y * y * y + std::sqrt(1107.421875) * y * y * y * y * y * y * z * z - std::sqrt(6029.296875) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(4429.6875) * x * x * x * x * x * x - std::sqrt(110742.1875) * x * x * x * x * y * y - std::sqrt(4429.6875) * x * x * x * x * z * z - std::sqrt(110742.1875) * x * x * y * y * y * y + std::sqrt(159468.75) * x * x * y * y * z * z + std::sqrt(4429.6875) * y * y * y * y * y * y - std::sqrt(4429.6875) * y * y * y * y * z * z) + e_3 * (std::sqrt(39867.1875) * x * x * x * x - std::sqrt(1435218.75) * x * x * y * y + std::sqrt(39867.1875) * y * y * y * y);

        pc_50[k] = e_0 * (std::sqrt(1.7303466796875) * x * x * x * x * x * x * x * x * x * z - std::sqrt(110.7421875) * x * x * x * x * x * x * x * y * y * z - std::sqrt(12.3046875) * x * x * x * x * x * x * x * z * z * z - std::sqrt(339.14794921875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(996.6796875) * x * x * x * x * x * y * y * z * z * z + std::sqrt(0.4921875) * x * x * x * x * x * z * z * z * z * z + std::sqrt(307.6171875) * x * x * x * y * y * y * y * z * z * z - std::sqrt(49.21875) * x * x * x * y * y * z * z * z * z * z + std::sqrt(43.2586669921875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(307.6171875) * x * y * y * y * y * y * y * z * z * z + std::sqrt(12.3046875) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(692.138671875) * x * x * x * x * x * x * x * z - std::sqrt(56063.232421875) * x * x * x * x * x * y * y * z - std::sqrt(1230.46875) * x * x * x * x * x * z * z * z - std::sqrt(17303.466796875) * x * x * x * y * y * y * y * z + std::sqrt(123046.875) * x * x * x * y * y * z * z * z + std::sqrt(17303.466796875) * x * y * y * y * y * y * y * z - std::sqrt(30761.71875) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(11074.21875) * x * x * x * x * x * z - std::sqrt(1107421.875) * x * x * x * y * y * z + std::sqrt(276855.46875) * x * y * y * y * y * z);

        pc_51[k] = e_0 * (0.234375 * x * x * x * x * x * x * x * x * x * x + 0.9375 * x * x * x * x * x * x * x * x * y * y - 5.625 * x * x * x * x * x * x * x * x * z * z + 1.40625 * x * x * x * x * x * x * y * y * y * y - 16.875 * x * x * x * x * x * x * y * y * z * z + 37.5 * x * x * x * x * x * x * z * z * z * z + 0.9375 * x * x * x * x * y * y * y * y * y * y - 16.875 * x * x * x * x * y * y * y * y * z * z + 75.0 * x * x * x * x * y * y * z * z * z * z - 45.0 * x * x * x * x * z * z * z * z * z * z + 0.234375 * x * x * y * y * y * y * y * y * y * y - 5.625 * x * x * y * y * y * y * y * y * z * z + 37.5 * x * x * y * y * y * y * z * z * z * z - 45.0 * x * x * y * y * z * z * z * z * z * z + 15.0 * x * x * z * z * z * z * z * z * z * z) + e_1 * (5.859375 * x * x * x * x * x * x * x * x + 17.8125 * x * x * x * x * x * x * y * y + 50.625 * x * x * x * x * x * x * z * z + 18.28125 * x * x * x * x * y * y * y * y + 95.625 * x * x * x * x * y * y * z * z - 37.5 * x * x * x * x * z * z * z * z + 6.5625 * x * x * y * y * y * y * y * y + 39.375 * x * x * y * y * y * y * z * z + 105.0 * x * x * z * z * z * z * z * z + 0.234375 * y * y * y * y * y * y * y * y - 5.625 * y * y * y * y * y * y * z * z + 37.5 * y * y * y * y * z * z * z * z - 45.0 * y * y * z * z * z * z * z * z + 15.0 * z * z * z * z * z * z * z * z) + e_2 * (116.25 * x * x * x * x * x * x + 236.25 * x * x * x * x * y * y + 315.0 * x * x * x * x * z * z + 123.75 * x * x * y * y * y * y + 405.0 * x * x * y * y * z * z + 675.0 * x * x * z * z * z * z + 3.75 * y * y * y * y * y * y + 90.0 * y * y * y * y * z * z - 225.0 * y * y * z * z * z * z + 240.0 * z * z * z * z * z * z) + e_3 * (765.0 * x * x * x * x + 855.0 * x * x * y * y + 2115.0 * x * x * z * z + 90.0 * y * y * y * y - 135.0 * y * y * z * z + 1350.0 * z * z * z * z) + e_4 * (1890.0 * x * x + 315.0 * y * y + 2520.0 * z * z) + e_5 * (945.0);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_52[k] = e_0 * (-std::sqrt(1.5380859375) * x * x * x * x * x * x * x * x * x * z - std::sqrt(6.15234375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(301.46484375) * x * x * x * x * x * x * x * z * z * z + std::sqrt(301.46484375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(1575.0) * x * x * x * x * x * z * z * z * z * z + std::sqrt(6.15234375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(301.46484375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(393.75) * x * x * x * z * z * z * z * z * z * z + std::sqrt(1.5380859375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(301.46484375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(1575.0) * x * y * y * y * y * z * z * z * z * z - std::sqrt(393.75) * x * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(24.609375) * x * x * x * x * x * x * x * z - std::sqrt(221.484375) * x * x * x * x * x * z * z * z - std::sqrt(221.484375) * x * x * x * y * y * y * y * z + std::sqrt(4823.4375) * x * x * x * y * y * z * z * z + std::sqrt(393.75) * x * x * x * z * z * z * z * z - std::sqrt(98.4375) * x * y * y * y * y * y * y * z + std::sqrt(7112.109375) * x * y * y * y * y * z * z * z - std::sqrt(31893.75) * x * y * y * z * z * z * z * z + std::sqrt(1575.0) * x * z * z * z * z * z * z * z) + e_2 * (std::sqrt(885.9375) * x * x * x * x * x * z + std::sqrt(3543.75) * x * x * x * y * y * z + std::sqrt(3543.75) * x * x * x * z * z * z + std::sqrt(885.9375) * x * y * y * y * y * z - std::sqrt(287043.75) * x * y * y * z * z * z + std::sqrt(127575.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(31893.75) * x * x * x * z - std::sqrt(173643.75) * x * y * y * z + std::sqrt(1148175.0) * x * z * z * z) + e_4 * (std::sqrt(694575.0) * x * z);

        pc_53[k] = e_0 * (-std::sqrt(0.0640869140625) * x * x * x * x * x * x * x * x * x * x + std::sqrt(25.634765625) * x * x * x * x * x * x * x * x * z * z + std::sqrt(2.30712890625) * x * x * x * x * x * x * y * y * y * y - std::sqrt(25.634765625) * x * x * x * x * x * x * y * y * z * z - std::sqrt(693.1640625) * x * x * x * x * x * x * z * z * z * z + std::sqrt(4.1015625) * x * x * x * x * y * y * y * y * y * y - std::sqrt(640.869140625) * x * x * x * x * y * y * y * y * z * z + std::sqrt(2772.65625) * x * x * x * x * y * y * z * z * z * z + std::sqrt(262.5) * x * x * x * x * z * z * z * z * z * z + std::sqrt(0.5767822265625) * x * x * y * y * y * y * y * y * y * y - std::sqrt(230.712890625) * x * x * y * y * y * y * y * y * z * z + std::sqrt(6238.4765625) * x * x * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * x * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(40.0543212890625) * x * x * x * x * x * x * x * x + std::sqrt(16.40625) * x * x * x * x * x * x * y * y - std::sqrt(452.197265625) * x * x * x * x * x * x * z * z + std::sqrt(775.45166015625) * x * x * x * x * y * y * y * y + std::sqrt(7761.181640625) * x * x * x * x * y * y * z * z - std::sqrt(9847.8515625) * x * x * x * x * z * z * z * z + std::sqrt(332.2265625) * x * x * y * y * y * y * y * y + std::sqrt(8868.603515625) * x * x * y * y * y * y * z * z + std::sqrt(3691.40625) * x * x * y * y * z * z * z * z + std::sqrt(2362.5) * x * x * z * z * z * z * z * z + std::sqrt(0.5767822265625) * y * y * y * y * y * y * y * y - std::sqrt(230.712890625) * y * y * y * y * y * y * z * z + std::sqrt(6238.4765625) * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(9450.0) * x * x * x * x * x * x + std::sqrt(24953.90625) * x * x * x * x * y * y - std::sqrt(170690.625) * x * x * x * x * z * z + std::sqrt(71465.625) * x * x * y * y * y * y + std::sqrt(531562.5) * x * x * y * y * z * z + std::sqrt(9450.0) * x * x * z * z * z * z + std::sqrt(147.65625) * y * y * y * y * y * y + std::sqrt(28940.625) * y * y * y * y * z * z - std::sqrt(9450.0) * y * y * z * z * z * z) + e_3 * (-std::sqrt(326172.65625) * x * x * x * x + std::sqrt(1196015.625) * x * x * y * y - std::sqrt(151200.0) * x * x * z * z + std::sqrt(42672.65625) * y * y * y * y + std::sqrt(151200.0) * y * y * z * z) + e_4 * (-std::sqrt(463050.0) * x * x + std::sqrt(463050.0) * y * y);

        pc_54[k] = e_0 * (std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * x * z - std::sqrt(18.45703125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(166.11328125) * x * x * x * x * x * x * x * z * z * z - std::sqrt(115.3564453125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(4152.83203125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(73.828125) * x * x * x * x * x * z * z * z * z * z - std::sqrt(18.45703125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(4152.83203125) * x * x * x * y * y * y * y * z * z * z - std::sqrt(2657.8125) * x * x * x * y * y * z * z * z * z * z + std::sqrt(1.153564453125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(166.11328125) * x * y * y * y * y * y * y * z * z * z + std::sqrt(73.828125) * x * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(18.45703125) * x * x * x * x * x * x * x * z + std::sqrt(1495.01953125) * x * x * x * x * x * y * y * z - std::sqrt(14470.3125) * x * x * x * x * x * z * z * z + std::sqrt(461.42578125) * x * x * x * y * y * y * y * z + std::sqrt(265781.25) * x * x * x * y * y * z * z * z + std::sqrt(1181.25) * x * x * x * z * z * z * z * z - std::sqrt(461.42578125) * x * y * y * y * y * y * y * z + std::sqrt(7382.8125) * x * y * y * y * y * z * z * z - std::sqrt(10631.25) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(42525.0) * x * x * x * x * x * z + std::sqrt(1063125.0) * x * x * x * y * y * z - std::sqrt(29531.25) * x * x * x * z * z * z + std::sqrt(265781.25) * x * y * y * z * z * z) + e_3 * (-std::sqrt(265781.25) * x * x * x * z + std::sqrt(2392031.25) * x * y * y * z);

        pc_55[k] = e_0 * (std::sqrt(0.1153564453125) * x * x * x * x * x * x * x * x * x * x - std::sqrt(7.3828125) * x * x * x * x * x * x * x * x * y * y - std::sqrt(16.611328125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(22.60986328125) * x * x * x * x * x * x * y * y * y * y + std::sqrt(1345.517578125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(7.3828125) * x * x * x * x * x * x * z * z * z * z + std::sqrt(415.283203125) * x * x * x * x * y * y * y * y * z * z - std::sqrt(738.28125) * x * x * x * x * y * y * z * z * z * z + std::sqrt(2.8839111328125) * x * x * y * y * y * y * y * y * y * y - std::sqrt(415.283203125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(184.5703125) * x * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(72.0977783203125) * x * x * x * x * x * x * x * x - std::sqrt(4614.2578125) * x * x * x * x * x * x * y * y - std::sqrt(3737.548828125) * x * x * x * x * x * x * z * z - std::sqrt(2595.52001953125) * x * x * x * x * y * y * y * y + std::sqrt(259552.001953125) * x * x * x * x * y * y * z * z + std::sqrt(184.5703125) * x * x * x * x * z * z * z * z + std::sqrt(738.28125) * x * x * y * y * y * y * y * y - std::sqrt(10382.080078125) * x * x * y * y * y * y * z * z - std::sqrt(6644.53125) * x * x * y * y * z * z * z * z + std::sqrt(2.8839111328125) * y * y * y * y * y * y * y * y - std::sqrt(415.283203125) * y * y * y * y * y * y * z * z + std::sqrt(184.5703125) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(2953.125) * x * x * x * x * x * x - std::sqrt(166113.28125) * x * x * x * x * y * y - std::sqrt(26578.125) * x * x * x * x * z * z + std::sqrt(956812.5) * x * x * y * y * z * z + std::sqrt(738.28125) * y * y * y * y * y * y - std::sqrt(26578.125) * y * y * y * y * z * z) + e_3 * (std::sqrt(6644.53125) * x * x * x * x - std::sqrt(239203.125) * x * x * y * y + std::sqrt(6644.53125) * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_56[k] = e_0 * (6.5625 * x * x * x * x * x * x * x * x * z * z - 26.25 * x * x * x * x * x * x * z * z * z * z - 13.125 * x * x * x * x * y * y * y * y * z * z + 26.25 * x * x * x * x * y * y * z * z * z * z + 26.25 * x * x * x * x * z * z * z * z * z * z + 26.25 * x * x * y * y * y * y * z * z * z * z - 52.5 * x * x * y * y * z * z * z * z * z * z + 6.5625 * y * y * y * y * y * y * y * y * z * z - 26.25 * y * y * y * y * y * y * z * z * z * z + 26.25 * y * y * y * y * z * z * z * z * z * z) + e_1 * (6.5625 * x * x * x * x * x * x * x * x + 26.25 * x * x * x * x * x * x * z * z - 13.125 * x * x * x * x * y * y * y * y + 78.75 * x * x * x * x * y * y * z * z + 26.25 * x * x * x * x * z * z * z * z + 78.75 * x * x * y * y * y * y * z * z - 472.5 * x * x * y * y * z * z * z * z + 105.0 * x * x * z * z * z * z * z * z + 6.5625 * y * y * y * y * y * y * y * y + 26.25 * y * y * y * y * y * y * z * z + 26.25 * y * y * y * y * z * z * z * z + 105.0 * y * y * z * z * z * z * z * z) + e_2 * (105.0 * x * x * x * x * x * x + 315.0 * x * x * x * x * z * z - 945.0 * x * x * y * y * z * z + 630.0 * x * x * z * z * z * z + 105.0 * y * y * y * y * y * y + 315.0 * y * y * y * y * z * z + 630.0 * y * y * z * z * z * z + 105.0 * z * z * z * z * z * z) + e_3 * (630.0 * x * x * x * x - 315.0 * x * x * y * y + 1575.0 * x * x * z * z + 630.0 * y * y * y * y + 1575.0 * y * y * z * z + 945.0 * z * z * z * z) + e_4 * (1260.0 * x * x + 1260.0 * y * y + 2205.0 * z * z) + e_5 * (945.0);

        pc_57[k] = e_0 * (std::sqrt(1.79443359375) * x * x * x * x * x * x * x * x * x * z - std::sqrt(7.177734375) * x * x * x * x * x * x * x * y * y * z - std::sqrt(179.443359375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(28.7109375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1614.990234375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(459.375) * x * x * x * x * x * z * z * z * z * z + std::sqrt(7.177734375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(179.443359375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(7350.0) * x * x * x * y * y * z * z * z * z * z + std::sqrt(16.14990234375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(1614.990234375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(4134.375) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(28.7109375) * x * x * x * x * x * x * x * z + std::sqrt(1033.59375) * x * x * x * x * x * y * y * z - std::sqrt(717.7734375) * x * x * x * x * x * z * z * z + std::sqrt(717.7734375) * x * x * x * y * y * y * y * z - std::sqrt(140683.59375) * x * x * x * y * y * z * z * z + std::sqrt(16537.5) * x * x * x * z * z * z * z * z + std::sqrt(6459.9609375) * x * y * y * y * y * z * z * z + std::sqrt(16537.5) * x * y * y * z * z * z * z * z) + e_2 * (std::sqrt(1033.59375) * x * x * x * x * x * z - std::sqrt(103359.375) * x * x * x * y * y * z + std::sqrt(103359.375) * x * x * x * z * z * z + std::sqrt(25839.84375) * x * y * y * y * y * z + std::sqrt(103359.375) * x * y * y * z * z * z + std::sqrt(66150.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(103359.375) * x * x * x * z + std::sqrt(103359.375) * x * y * y * z + std::sqrt(1653750.0) * x * z * z * z) + e_4 * (std::sqrt(1653750.0) * x * z);

        pc_58[k] = e_0 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * z * z + std::sqrt(1162.79296875) * x * x * x * x * x * x * y * y * z * z + std::sqrt(129.19921875) * x * x * x * x * x * x * z * z * z * z - std::sqrt(6330.76171875) * x * x * x * x * y * y * z * z * z * z - std::sqrt(1162.79296875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(6330.76171875) * x * x * y * y * y * y * z * z * z * z + std::sqrt(32.2998046875) * y * y * y * y * y * y * y * y * z * z - std::sqrt(129.19921875) * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x + std::sqrt(1162.79296875) * x * x * x * x * x * x * y * y - std::sqrt(3229.98046875) * x * x * x * x * x * x * z * z + std::sqrt(1162.79296875) * x * x * x * x * y * y * z * z + std::sqrt(8268.75) * x * x * x * x * z * z * z * z - std::sqrt(1162.79296875) * x * x * y * y * y * y * y * y - std::sqrt(1162.79296875) * x * x * y * y * y * y * z * z + std::sqrt(32.2998046875) * y * y * y * y * y * y * y * y + std::sqrt(3229.98046875) * y * y * y * y * y * y * z * z - std::sqrt(8268.75) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(8268.75) * x * x * x * x * x * x + std::sqrt(74418.75) * x * x * x * x * y * y - std::sqrt(18604.6875) * x * x * x * x * z * z - std::sqrt(74418.75) * x * x * y * y * y * y + std::sqrt(74418.75) * x * x * z * z * z * z + std::sqrt(8268.75) * y * y * y * y * y * y + std::sqrt(18604.6875) * y * y * y * y * z * z - std::sqrt(74418.75) * y * y * z * z * z * z) + e_3 * (-std::sqrt(167442.1875) * x * x * x * x + std::sqrt(74418.75) * x * x * z * z + std::sqrt(167442.1875) * y * y * y * y - std::sqrt(74418.75) * y * y * z * z) + e_4 * (-std::sqrt(297675.0) * x * x + std::sqrt(297675.0) * y * y);

        pc_59[k] = e_0 * (-std::sqrt(3.22998046875) * x * x * x * x * x * x * x * x * x * z + std::sqrt(322.998046875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(12.919921875) * x * x * x * x * x * x * x * z * z * z - std::sqrt(51.6796875) * x * x * x * x * x * y * y * y * y * z - std::sqrt(1563.310546875) * x * x * x * x * x * y * y * z * z * z - std::sqrt(322.998046875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(2906.982421875) * x * x * x * y * y * y * y * z * z * z + std::sqrt(80.74951171875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(322.998046875) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1291.9921875) * x * x * x * x * x * x * x * z + std::sqrt(46511.71875) * x * x * x * x * x * y * y * z + std::sqrt(1291.9921875) * x * x * x * x * x * z * z * z - std::sqrt(32299.8046875) * x * x * x * y * y * y * y * z - std::sqrt(5167.96875) * x * x * x * y * y * z * z * z + std::sqrt(20671.875) * x * y * y * y * y * y * y * z - std::sqrt(11627.9296875) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(46511.71875) * x * x * x * x * x * z + std::sqrt(186046.875) * x * x * x * y * y * z + std::sqrt(20671.875) * x * x * x * z * z * z + std::sqrt(418605.46875) * x * y * y * y * y * z - std::sqrt(186046.875) * x * y * y * z * z * z) + e_3 * (-std::sqrt(186046.875) * x * x * x * z + std::sqrt(1674421.875) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_60[k] = e_0 * (0.2734375 * x * x * x * x * x * x * x * x * x * x - 1.09375 * x * x * x * x * x * x * x * x * y * y - 4.375 * x * x * x * x * x * x * x * x * z * z - 0.546875 * x * x * x * x * x * x * y * y * y * y + 21.875 * x * x * x * x * x * x * y * y * z * z + 17.5 * x * x * x * x * x * x * z * z * z * z + 3.28125 * x * x * x * x * y * y * y * y * y * y - 13.125 * x * x * x * x * y * y * y * y * z * z - 105.0 * x * x * x * x * y * y * z * z * z * z + 2.4609375 * x * x * y * y * y * y * y * y * y * y - 39.375 * x * x * y * y * y * y * y * y * z * z + 157.5 * x * x * y * y * y * y * z * z * z * z) + e_1 * (6.8359375 * x * x * x * x * x * x * x * x - 12.03125 * x * x * x * x * x * x * y * y + 4.375 * x * x * x * x * x * x * z * z + 27.890625 * x * x * x * x * y * y * y * y - 380.625 * x * x * x * x * y * y * z * z + 157.5 * x * x * x * x * z * z * z * z + 49.21875 * x * x * y * y * y * y * y * y + 275.625 * x * x * y * y * y * y * z * z + 315.0 * x * x * y * y * z * z * z * z + 2.4609375 * y * y * y * y * y * y * y * y - 39.375 * y * y * y * y * y * y * z * z + 157.5 * y * y * y * y * z * z * z * z) + e_2 * (91.875 * x * x * x * x * x * x - 196.875 * x * x * x * x * y * y + 315.0 * x * x * x * x * z * z + 590.625 * x * x * y * y * y * y + 630.0 * x * x * y * y * z * z + 630.0 * x * x * z * z * z * z + 39.375 * y * y * y * y * y * y + 315.0 * y * y * y * y * z * z + 630.0 * y * y * z * z * z * z) + e_3 * (498.75 * x * x * x * x + 997.5 * x * x * y * y + 2100.0 * x * x * z * z + 498.75 * y * y * y * y + 2100.0 * y * y * z * z + 420.0 * z * z * z * z) + e_4 * (1522.5 * x * x + 1522.5 * y * y + 1680.0 * z * z) + e_5 * (945.0);

        pc_61[k] = e_0 * (-std::sqrt(1.3458251953125) * x * x * x * x * x * x * x * x * x * z + std::sqrt(86.1328125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(86.1328125) * x * x * x * x * x * x * x * z * z * z - std::sqrt(134.58251953125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(6976.7578125) * x * x * x * x * x * y * y * z * z * z - std::sqrt(344.53125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(31093.9453125) * x * x * x * y * y * y * y * z * z * z + std::sqrt(12.1124267578125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(775.1953125) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * x * z - std::sqrt(15697.705078125) * x * x * x * x * x * y * y * z + std::sqrt(12403.125) * x * x * x * x * x * z * z * z + std::sqrt(18109.423828125) * x * x * x * y * y * y * y * z + std::sqrt(49612.5) * x * x * x * y * y * z * z * z - std::sqrt(1744.189453125) * x * y * y * y * y * y * y * z + std::sqrt(12403.125) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * x * z + std::sqrt(12403.125) * x * x * x * y * y * z + std::sqrt(446512.5) * x * x * x * z * z * z + std::sqrt(3100.78125) * x * y * y * y * y * z + std::sqrt(446512.5) * x * y * y * z * z * z) + e_3 * (std::sqrt(793800.0) * x * x * x * z + std::sqrt(793800.0) * x * y * y * z + std::sqrt(793800.0) * x * z * z * z) + e_4 * (std::sqrt(2431012.5) * x * z);

        pc_62[k] = e_0 * (-std::sqrt(0.13458251953125) * x * x * x * x * x * x * x * x * x * x + std::sqrt(19.3798828125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(8.61328125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(65.137939453125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(1455.64453125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(53.8330078125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(10551.26953125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(30.28106689453125) * x * x * y * y * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(84.11407470703125) * x * x * x * x * x * x * x * x + std::sqrt(1345.8251953125) * x * x * x * x * x * x * y * y + std::sqrt(1937.98828125) * x * x * x * x * x * x * z * z - std::sqrt(16486.358642578125) * x * x * x * x * y * y * y * y + std::sqrt(1937.98828125) * x * x * x * x * y * y * z * z + std::sqrt(4360.4736328125) * x * x * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * y * y * y * y * z * z + std::sqrt(30.28106689453125) * y * y * y * y * y * y * y * y - std::sqrt(1937.98828125) * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(7751.953125) * x * x * x * x * x * x - std::sqrt(7751.953125) * x * x * x * x * y * y + std::sqrt(124031.25) * x * x * x * x * z * z + std::sqrt(7751.953125) * x * x * y * y * y * y + std::sqrt(7751.953125) * y * y * y * y * y * y - std::sqrt(124031.25) * y * y * y * y * z * z) + e_3 * (-std::sqrt(124031.25) * x * x * x * x + std::sqrt(496125.0) * x * x * z * z + std::sqrt(124031.25) * y * y * y * y - std::sqrt(496125.0) * y * y * z * z) + e_4 * (-std::sqrt(124031.25) * x * x + std::sqrt(124031.25) * y * y);

        pc_63[k] = e_0 * (4.921875 * x * x * x * x * x * x * x * x * z * z - 59.0625 * x * x * x * x * x * x * y * y * z * z + 187.03125 * x * x * x * x * y * y * y * y * z * z - 59.0625 * x * x * y * y * y * y * y * y * z * z + 4.921875 * y * y * y * y * y * y * y * y * z * z) + e_1 * (4.921875 * x * x * x * x * x * x * x * x - 59.0625 * x * x * x * x * x * x * y * y + 78.75 * x * x * x * x * x * x * z * z + 187.03125 * x * x * x * x * y * y * y * y + 236.25 * x * x * x * x * y * y * z * z - 59.0625 * x * x * y * y * y * y * y * y + 236.25 * x * x * y * y * y * y * z * z + 4.921875 * y * y * y * y * y * y * y * y + 78.75 * y * y * y * y * y * y * z * z) + e_2 * (78.75 * x * x * x * x * x * x + 236.25 * x * x * x * x * y * y + 708.75 * x * x * x * x * z * z + 236.25 * x * x * y * y * y * y + 1417.5 * x * x * y * y * z * z + 78.75 * y * y * y * y * y * y + 708.75 * y * y * y * y * z * z) + e_3 * (708.75 * x * x * x * x + 1417.5 * x * x * y * y + 1890.0 * x * x * z * z + 708.75 * y * y * y * y + 1890.0 * y * y * z * z) + e_4 * (1890.0 * x * x + 1890.0 * y * y + 945.0 * z * z) + e_5 * (945.0);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_64[k] = e_0 * (std::sqrt(2.4224853515625) * x * x * x * x * x * x * x * x * x * z - std::sqrt(620.15625) * x * x * x * x * x * x * x * y * y * z + std::sqrt(10552.34619140625) * x * x * x * x * x * y * y * y * y * z - std::sqrt(3875.9765625) * x * x * x * y * y * y * y * y * y * z + std::sqrt(60.5621337890625) * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(968.994140625) * x * x * x * x * x * x * x * z + std::sqrt(8720.947265625) * x * x * x * x * x * y * y * z + std::sqrt(8720.947265625) * x * x * x * y * y * y * y * z + std::sqrt(968.994140625) * x * y * y * y * y * y * y * z) + e_2 * (std::sqrt(139535.15625) * x * x * x * x * x * z + std::sqrt(558140.625) * x * x * x * y * y * z + std::sqrt(139535.15625) * x * y * y * y * y * z) + e_3 * (std::sqrt(2232562.5) * x * x * x * z + std::sqrt(2232562.5) * x * y * y * z) + e_4 * (std::sqrt(2232562.5) * x * z);

        pc_65[k] = e_0 * (0.4921875 * x * x * x * x * x * x * x * x * x * x - 9.84375 * x * x * x * x * x * x * x * x * y * y + 54.140625 * x * x * x * x * x * x * y * y * y * y - 49.21875 * x * x * x * x * y * y * y * y * y * y + 12.3046875 * x * x * y * y * y * y * y * y * y * y) + e_1 * (12.3046875 * x * x * x * x * x * x * x * x + 49.21875 * x * x * x * x * x * x * y * y + 73.828125 * x * x * x * x * y * y * y * y + 49.21875 * x * x * y * y * y * y * y * y + 12.3046875 * y * y * y * y * y * y * y * y) + e_2 * (196.875 * x * x * x * x * x * x + 590.625 * x * x * x * x * y * y + 590.625 * x * x * y * y * y * y + 196.875 * y * y * y * y * y * y) + e_3 * (1181.25 * x * x * x * x + 2362.5 * x * x * y * y + 1181.25 * y * y * y * y) + e_4 * (2362.5 * x * x + 2362.5 * y * y) + e_5 * (945.0);
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[121] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 2, 13, 24, 25, 26, 27, 28, 29, 30, 31, 32, 3, 14, 25, 36, 37, 38, 39, 40, 41, 42, 43, 4, 15, 26, 37, 48, 49, 50, 51, 52, 53, 54, 5, 16, 27, 38, 49, 60, 61, 62, 63, 64, 65, 6, 17, 28, 39, 50, 61, 72, 73, 74, 75, 76, 7, 18, 29, 40, 51, 62, 73, 84, 85, 86, 87, 8, 19, 30, 41, 52, 63, 74, 85, 96, 97, 98, 9, 20, 31, 42, 53, 64, 75, 86, 97, 108, 109, 10, 21, 32, 43, 54, 65, 76, 87, 98, 109, 120};

    for (size_t m = 0; m < 121; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
