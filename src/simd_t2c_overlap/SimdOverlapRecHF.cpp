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



#include "SimdOverlapRecHF.hpp"

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
compute_hf_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHF.compute_hf_overlap: Basis functions must be of angular momenta five and three"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHF.compute_hf_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 4);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 77 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fal * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fal * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fal * fbe * fh * fh;

        const auto f_3 = fbase * fal * fal * fh * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
            pe_3[k] += f_3 * fss;
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
    auto *pc_11 = values + 11 * nvalues;
    auto *pc_12 = values + 12 * nvalues;
    auto *pc_13 = values + 13 * nvalues;
    auto *pc_14 = values + 14 * nvalues;
    auto *pc_15 = values + 15 * nvalues;
    auto *pc_16 = values + 16 * nvalues;
    auto *pc_17 = values + 17 * nvalues;
    auto *pc_18 = values + 18 * nvalues;
    auto *pc_19 = values + 19 * nvalues;
    auto *pc_20 = values + 20 * nvalues;
    auto *pc_21 = values + 21 * nvalues;
    auto *pc_22 = values + 22 * nvalues;
    auto *pc_23 = values + 23 * nvalues;
    auto *pc_24 = values + 24 * nvalues;
    auto *pc_25 = values + 25 * nvalues;
    auto *pc_26 = values + 26 * nvalues;
    auto *pc_27 = values + 27 * nvalues;
    auto *pc_28 = values + 28 * nvalues;
    auto *pc_29 = values + 29 * nvalues;
    auto *pc_30 = values + 30 * nvalues;
    auto *pc_31 = values + 31 * nvalues;
    auto *pc_32 = values + 32 * nvalues;
    auto *pc_33 = values + 33 * nvalues;
    auto *pc_34 = values + 34 * nvalues;
    auto *pc_35 = values + 35 * nvalues;
    auto *pc_36 = values + 36 * nvalues;
    auto *pc_37 = values + 37 * nvalues;
    auto *pc_38 = values + 38 * nvalues;
    auto *pc_39 = values + 39 * nvalues;
    auto *pc_40 = values + 40 * nvalues;
    auto *pc_41 = values + 41 * nvalues;
    auto *pc_42 = values + 42 * nvalues;
    auto *pc_43 = values + 43 * nvalues;
    auto *pc_44 = values + 44 * nvalues;
    auto *pc_45 = values + 45 * nvalues;
    auto *pc_46 = values + 46 * nvalues;
    auto *pc_47 = values + 47 * nvalues;
    auto *pc_48 = values + 48 * nvalues;
    auto *pc_49 = values + 49 * nvalues;
    auto *pc_50 = values + 50 * nvalues;
    auto *pc_51 = values + 51 * nvalues;
    auto *pc_52 = values + 52 * nvalues;
    auto *pc_53 = values + 53 * nvalues;
    auto *pc_54 = values + 54 * nvalues;
    auto *pc_55 = values + 55 * nvalues;
    auto *pc_56 = values + 56 * nvalues;
    auto *pc_57 = values + 57 * nvalues;
    auto *pc_58 = values + 58 * nvalues;
    auto *pc_59 = values + 59 * nvalues;
    auto *pc_60 = values + 60 * nvalues;
    auto *pc_61 = values + 61 * nvalues;
    auto *pc_62 = values + 62 * nvalues;
    auto *pc_63 = values + 63 * nvalues;
    auto *pc_64 = values + 64 * nvalues;
    auto *pc_65 = values + 65 * nvalues;
    auto *pc_66 = values + 66 * nvalues;
    auto *pc_67 = values + 67 * nvalues;
    auto *pc_68 = values + 68 * nvalues;
    auto *pc_69 = values + 69 * nvalues;
    auto *pc_70 = values + 70 * nvalues;
    auto *pc_71 = values + 71 * nvalues;
    auto *pc_72 = values + 72 * nvalues;
    auto *pc_73 = values + 73 * nvalues;
    auto *pc_74 = values + 74 * nvalues;
    auto *pc_75 = values + 75 * nvalues;
    auto *pc_76 = values + 76 * nvalues;

    // NOTE: the components are formed in 20 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

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

        pc_0[k] = e_0 * (std::sqrt(69.2138671875) * x * x * x * x * x * x * y * y - std::sqrt(376.8310546875) * x * x * x * x * y * y * y * y + std::sqrt(51.9873046875) * x * x * y * y * y * y * y * y - std::sqrt(0.3076171875) * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(69.2138671875) * x * x * x * x * x * x + std::sqrt(69.2138671875) * x * x * x * x * y * y - std::sqrt(69.2138671875) * x * x * y * y * y * y - std::sqrt(69.2138671875) * y * y * y * y * y * y) + e_2 * (std::sqrt(4429.6875) * x * x * x * x - std::sqrt(4429.6875) * y * y * y * y) + e_3 * (std::sqrt(17718.75) * x * x - std::sqrt(17718.75) * y * y);

        pc_1[k] = e_0 * (std::sqrt(184.5703125) * x * x * x * x * x * y * y * z - std::sqrt(738.28125) * x * x * x * y * y * y * y * z + std::sqrt(7.3828125) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(184.5703125) * x * x * x * x * x * z - std::sqrt(738.28125) * x * x * x * y * y * z - std::sqrt(1661.1328125) * x * y * y * y * y * z) + e_2 * (std::sqrt(2953.125) * x * x * x * z - std::sqrt(26578.125) * x * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x * y * y + std::sqrt(4.6142578125) * x * x * x * x * y * y * y * y + std::sqrt(73.828125) * x * x * x * x * y * y * z * z + std::sqrt(14.9501953125) * x * x * y * y * y * y * y * y - std::sqrt(295.3125) * x * x * y * y * y * y * z * z - std::sqrt(0.1845703125) * y * y * y * y * y * y * y * y + std::sqrt(2.953125) * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x - std::sqrt(115.3564453125) * x * x * x * x * y * y + std::sqrt(73.828125) * x * x * x * x * z * z + std::sqrt(2883.9111328125) * x * x * y * y * y * y - std::sqrt(2657.8125) * x * x * y * y * z * z - std::sqrt(41.5283203125) * y * y * y * y * y * y + std::sqrt(73.828125) * y * y * y * y * z * z) + e_2 * (-std::sqrt(295.3125) * x * x * x * x + std::sqrt(10631.25) * x * x * y * y - std::sqrt(295.3125) * y * y * y * y);

        pc_3[k] = e_0 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * y * z + std::sqrt(27.685546875) * x * x * x * x * y * y * y * z + std::sqrt(12.3046875) * x * x * x * x * y * z * z * z + std::sqrt(89.701171875) * x * x * y * y * y * y * y * z - std::sqrt(49.21875) * x * x * y * y * y * z * z * z - std::sqrt(1.107421875) * y * y * y * y * y * y * y * z + std::sqrt(0.4921875) * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(2768.5546875) * x * x * x * x * y * z + std::sqrt(11074.21875) * x * x * y * y * y * z - std::sqrt(110.7421875) * y * y * y * y * y * z);
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

        pc_4[k] = e_0 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x * x * y + std::sqrt(4.6142578125) * x * x * x * x * x * y * y * y + std::sqrt(73.828125) * x * x * x * x * x * y * z * z + std::sqrt(14.9501953125) * x * x * x * y * y * y * y * y - std::sqrt(295.3125) * x * x * x * y * y * y * z * z - std::sqrt(0.1845703125) * x * y * y * y * y * y * y * y + std::sqrt(2.953125) * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(904.39453125) * x * x * x * x * x * y + std::sqrt(1845.703125) * x * x * x * y * y * y + std::sqrt(1181.25) * x * x * x * y * z * z + std::sqrt(18.45703125) * x * y * y * y * y * y - std::sqrt(1181.25) * x * y * y * y * z * z) + e_2 * (-std::sqrt(4725.0) * x * x * x * y + std::sqrt(4725.0) * x * y * y * y);

        pc_5[k] = e_0 * (std::sqrt(46.142578125) * x * x * x * x * x * x * y * z - std::sqrt(415.283203125) * x * x * x * x * y * y * y * z + std::sqrt(223.330078125) * x * x * y * y * y * y * y * z - std::sqrt(1.845703125) * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1661.1328125) * x * x * x * x * y * z + std::sqrt(738.28125) * x * x * y * y * y * z - std::sqrt(184.5703125) * y * y * y * y * y * z) + e_2 * (std::sqrt(26578.125) * x * x * y * z - std::sqrt(2953.125) * y * y * y * z);

        pc_6[k] = e_0 * (std::sqrt(7.6904296875) * x * x * x * x * x * x * x * y - std::sqrt(192.2607421875) * x * x * x * x * x * y * y * y + std::sqrt(295.6201171875) * x * x * x * y * y * y * y * y - std::sqrt(2.7685546875) * x * y * y * y * y * y * y * y) + e_1 * (std::sqrt(276.85546875) * x * x * x * x * x * y + std::sqrt(1107.421875) * x * x * x * y * y * y + std::sqrt(276.85546875) * x * y * y * y * y * y) + e_2 * (std::sqrt(17718.75) * x * x * x * y + std::sqrt(17718.75) * x * y * y * y) + e_3 * (std::sqrt(70875.0) * x * y);

        pc_7[k] = e_0 * (std::sqrt(442.96875) * x * x * x * x * x * y * y * z - std::sqrt(787.5) * x * x * x * y * y * y * y * z + std::sqrt(49.21875) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(442.96875) * x * x * x * x * x * z + std::sqrt(1771.875) * x * x * x * y * y * z + std::sqrt(442.96875) * x * y * y * y * y * z) + e_2 * (std::sqrt(15946.875) * x * x * x * z + std::sqrt(15946.875) * x * y * y * z) + e_3 * (std::sqrt(28350.0) * x * z);
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

        pc_8[k] = e_0 * (std::sqrt(1181.25) * x * x * x * x * y * y * z * z - std::sqrt(1181.25) * x * x * y * y * y * y * z * z) + e_1 * (std::sqrt(1181.25) * x * x * x * x * y * y + std::sqrt(1181.25) * x * x * x * x * z * z - std::sqrt(1181.25) * x * x * y * y * y * y - std::sqrt(1181.25) * y * y * y * y * z * z) + e_2 * (std::sqrt(1181.25) * x * x * x * x + std::sqrt(10631.25) * x * x * z * z - std::sqrt(1181.25) * y * y * y * y - std::sqrt(10631.25) * y * y * z * z) + e_3 * (std::sqrt(10631.25) * x * x - std::sqrt(10631.25) * y * y);

        pc_9[k] = e_0 * (-std::sqrt(29.53125) * x * x * x * x * x * y * y * z + std::sqrt(472.5) * x * x * x * y * y * z * z * z + std::sqrt(29.53125) * x * y * y * y * y * y * y * z - std::sqrt(472.5) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(29.53125) * x * x * x * x * x * z + std::sqrt(118.125) * x * x * x * y * y * z + std::sqrt(472.5) * x * x * x * z * z * z + std::sqrt(265.78125) * x * y * y * y * y * z - std::sqrt(4252.5) * x * y * y * z * z * z) + e_2 * (std::sqrt(118.125) * x * x * x * z - std::sqrt(1063.125) * x * y * y * z);

        pc_10[k] = e_0 * (-std::sqrt(177.1875) * x * x * x * x * x * y * z * z + std::sqrt(78.75) * x * x * x * y * z * z * z * z + std::sqrt(177.1875) * x * y * y * y * y * y * z * z - std::sqrt(78.75) * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(177.1875) * x * x * x * x * x * y - std::sqrt(6378.75) * x * x * x * y * z * z + std::sqrt(177.1875) * x * y * y * y * y * y + std::sqrt(6378.75) * x * y * y * y * z * z) + e_2 * (-std::sqrt(11340.0) * x * x * x * y + std::sqrt(11340.0) * x * y * y * y);

        pc_11[k] = e_0 * (-std::sqrt(29.53125) * x * x * x * x * x * x * y * z + std::sqrt(472.5) * x * x * x * x * y * z * z * z + std::sqrt(29.53125) * x * x * y * y * y * y * y * z - std::sqrt(472.5) * x * x * y * y * y * z * z * z) + e_1 * (-std::sqrt(265.78125) * x * x * x * x * y * z - std::sqrt(118.125) * x * x * y * y * y * z + std::sqrt(4252.5) * x * x * y * z * z * z + std::sqrt(29.53125) * y * y * y * y * y * z - std::sqrt(472.5) * y * y * y * z * z * z) + e_2 * (std::sqrt(1063.125) * x * x * y * z - std::sqrt(118.125) * y * y * y * z);
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

        pc_12[k] = e_0 * (std::sqrt(295.3125) * x * x * x * x * x * y * z * z - std::sqrt(1181.25) * x * x * x * y * y * y * z * z + std::sqrt(295.3125) * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(295.3125) * x * x * x * x * x * y - std::sqrt(1181.25) * x * x * x * y * y * y + std::sqrt(4725.0) * x * x * x * y * z * z + std::sqrt(295.3125) * x * y * y * y * y * y + std::sqrt(4725.0) * x * y * y * y * z * z) + e_2 * (std::sqrt(4725.0) * x * x * x * y + std::sqrt(4725.0) * x * y * y * y + std::sqrt(42525.0) * x * y * z * z) + e_3 * (std::sqrt(42525.0) * x * y);

        pc_13[k] = e_0 * (std::sqrt(49.21875) * x * x * x * x * x * x * y * z - std::sqrt(787.5) * x * x * x * x * y * y * y * z + std::sqrt(442.96875) * x * x * y * y * y * y * y * z) + e_1 * (std::sqrt(442.96875) * x * x * x * x * y * z + std::sqrt(1771.875) * x * x * y * y * y * z + std::sqrt(442.96875) * y * y * y * y * y * z) + e_2 * (std::sqrt(15946.875) * x * x * y * z + std::sqrt(15946.875) * y * y * y * z) + e_3 * (std::sqrt(28350.0) * y * z);

        pc_14[k] = e_0 * (-std::sqrt(13.8427734375) * x * x * x * x * x * x * y * y - std::sqrt(1.5380859375) * x * x * x * x * y * y * y * y + std::sqrt(885.9375) * x * x * x * x * y * y * z * z + std::sqrt(4.2724609375) * x * x * y * y * y * y * y * y - std::sqrt(393.75) * x * x * y * y * y * y * z * z - std::sqrt(0.1708984375) * y * y * y * y * y * y * y * y + std::sqrt(10.9375) * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(13.8427734375) * x * x * x * x * x * x - std::sqrt(1121.2646484375) * x * x * x * x * y * y + std::sqrt(885.9375) * x * x * x * x * z * z + std::sqrt(13.8427734375) * x * x * y * y * y * y + std::sqrt(3543.75) * x * x * y * y * z * z - std::sqrt(38.4521484375) * y * y * y * y * y * y + std::sqrt(885.9375) * y * y * y * y * z * z) + e_2 * (-std::sqrt(885.9375) * x * x * x * x - std::sqrt(3543.75) * x * x * y * y + std::sqrt(14175.0) * x * x * z * z - std::sqrt(885.9375) * y * y * y * y + std::sqrt(14175.0) * y * y * z * z) + e_3 * (-std::sqrt(1575.0) * x * x - std::sqrt(1575.0) * y * y + std::sqrt(6300.0) * z * z);

        pc_15[k] = e_0 * (-std::sqrt(36.9140625) * x * x * x * x * x * y * y * z - std::sqrt(16.40625) * x * x * x * y * y * y * y * z + std::sqrt(2362.5) * x * x * x * y * y * z * z * z + std::sqrt(4.1015625) * x * y * y * y * y * y * y * z - std::sqrt(262.5) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(36.9140625) * x * x * x * x * x * z + std::sqrt(3691.40625) * x * x * x * y * y * z + std::sqrt(2362.5) * x * x * x * z * z * z - std::sqrt(922.8515625) * x * y * y * y * y * z + std::sqrt(2362.5) * x * y * y * z * z * z) + e_2 * (std::sqrt(5315.625) * x * x * x * z + std::sqrt(5315.625) * x * y * y * z + std::sqrt(9450.0) * x * z * z * z) + e_3 * (std::sqrt(37800.0) * x * z);
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

        pc_16[k] = e_0 * (std::sqrt(0.9228515625) * x * x * x * x * x * x * y * y + std::sqrt(2.5634765625) * x * x * x * x * y * y * y * y - std::sqrt(132.890625) * x * x * x * x * y * y * z * z + std::sqrt(0.1025390625) * x * x * y * y * y * y * y * y - std::sqrt(59.0625) * x * x * y * y * y * y * z * z + std::sqrt(945.0) * x * x * y * y * z * z * z * z - std::sqrt(0.1025390625) * y * y * y * y * y * y * y * y + std::sqrt(14.765625) * y * y * y * y * y * y * z * z - std::sqrt(105.0) * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(0.9228515625) * x * x * x * x * x * x + std::sqrt(155.9619140625) * x * x * x * x * y * y - std::sqrt(132.890625) * x * x * x * x * z * z + std::sqrt(45.2197265625) * x * x * y * y * y * y + std::sqrt(4784.0625) * x * x * y * y * z * z + std::sqrt(945.0) * x * x * z * z * z * z - std::sqrt(23.0712890625) * y * y * y * y * y * y - std::sqrt(132.890625) * y * y * y * y * z * z - std::sqrt(945.0) * y * y * z * z * z * z) + e_2 * (std::sqrt(59.0625) * x * x * x * x + std::sqrt(8505.0) * x * x * y * y + std::sqrt(8505.0) * x * x * z * z - std::sqrt(1476.5625) * y * y * y * y - std::sqrt(8505.0) * y * y * z * z) + e_3 * (std::sqrt(5906.25) * x * x - std::sqrt(5906.25) * y * y);

        pc_17[k] = e_0 * (std::sqrt(5.537109375) * x * x * x * x * x * x * y * z + std::sqrt(15.380859375) * x * x * x * x * y * y * y * z - std::sqrt(415.8984375) * x * x * x * x * y * z * z * z + std::sqrt(0.615234375) * x * x * y * y * y * y * y * z - std::sqrt(184.84375) * x * x * y * y * y * z * z * z + std::sqrt(157.5) * x * x * y * z * z * z * z * z - std::sqrt(0.615234375) * y * y * y * y * y * y * y * z + std::sqrt(46.2109375) * y * y * y * y * y * z * z * z - std::sqrt(17.5) * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(199.3359375) * x * x * x * x * y * z - std::sqrt(88.59375) * x * x * y * y * y * z - std::sqrt(1417.5) * x * x * y * z * z * z + std::sqrt(22.1484375) * y * y * y * y * y * z + std::sqrt(157.5) * y * y * y * z * z * z) + e_2 * (-std::sqrt(12757.5) * x * x * y * z + std::sqrt(1417.5) * y * y * y * z);

        pc_18[k] = e_0 * (std::sqrt(0.9228515625) * x * x * x * x * x * x * x * y + std::sqrt(2.5634765625) * x * x * x * x * x * y * y * y - std::sqrt(132.890625) * x * x * x * x * x * y * z * z + std::sqrt(0.1025390625) * x * x * x * y * y * y * y * y - std::sqrt(59.0625) * x * x * x * y * y * y * z * z + std::sqrt(945.0) * x * x * x * y * z * z * z * z - std::sqrt(0.1025390625) * x * y * y * y * y * y * y * y + std::sqrt(14.765625) * x * y * y * y * y * y * z * z - std::sqrt(105.0) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(180.87890625) * x * x * x * x * x * y + std::sqrt(132.890625) * x * x * x * y * y * y + std::sqrt(2126.25) * x * x * x * y * z * z - std::sqrt(3.69140625) * x * y * y * y * y * y - std::sqrt(2126.25) * x * y * y * y * z * z + std::sqrt(3780.0) * x * y * z * z * z * z) + e_2 * (std::sqrt(11576.25) * x * x * x * y - std::sqrt(236.25) * x * y * y * y + std::sqrt(34020.0) * x * y * z * z) + e_3 * (std::sqrt(23625.0) * x * y);

        pc_19[k] = e_0 * (-std::sqrt(9.228515625) * x * x * x * x * x * x * y * z + std::sqrt(1.025390625) * x * x * x * x * y * y * y * z + std::sqrt(590.625) * x * x * x * x * y * z * z * z + std::sqrt(9.228515625) * x * x * y * y * y * y * y * z - std::sqrt(1050.0) * x * x * y * y * y * z * z * z - std::sqrt(1.025390625) * y * y * y * y * y * y * y * z + std::sqrt(65.625) * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(922.8515625) * x * x * x * x * y * z - std::sqrt(3691.40625) * x * x * y * y * y * z + std::sqrt(2362.5) * x * x * y * z * z * z + std::sqrt(36.9140625) * y * y * y * y * y * z + std::sqrt(2362.5) * y * y * y * z * z * z) + e_2 * (std::sqrt(5315.625) * x * x * y * z + std::sqrt(5315.625) * y * y * y * z + std::sqrt(9450.0) * y * z * z * z) + e_3 * (std::sqrt(37800.0) * y * z);
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

        pc_20[k] = e_0 * (-std::sqrt(1.5380859375) * x * x * x * x * x * x * x * y + std::sqrt(8.3740234375) * x * x * x * x * x * y * y * y + std::sqrt(98.4375) * x * x * x * x * x * y * z * z + std::sqrt(8.3740234375) * x * x * x * y * y * y * y * y - std::sqrt(1093.75) * x * x * x * y * y * y * z * z - std::sqrt(1.5380859375) * x * y * y * y * y * y * y * y + std::sqrt(98.4375) * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(55.37109375) * x * x * x * x * x * y + std::sqrt(615.234375) * x * x * x * y * y * y - std::sqrt(55.37109375) * x * y * y * y * y * y);

        pc_21[k] = e_0 * (-std::sqrt(147.65625) * x * x * x * x * x * y * y * z - std::sqrt(65.625) * x * x * x * y * y * y * y * z + std::sqrt(590.625) * x * x * x * y * y * z * z * z + std::sqrt(16.40625) * x * y * y * y * y * y * y * z - std::sqrt(65.625) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(147.65625) * x * x * x * x * x * z - std::sqrt(9450.0) * x * x * x * y * y * z + std::sqrt(590.625) * x * x * x * z * z * z + std::sqrt(147.65625) * x * y * y * y * y * z + std::sqrt(590.625) * x * y * y * z * z * z) + e_2 * (-std::sqrt(5315.625) * x * x * x * z - std::sqrt(5315.625) * x * y * y * z + std::sqrt(2362.5) * x * z * z * z) + e_3 * (-std::sqrt(2362.5) * x * z);

        pc_22[k] = e_0 * (-std::sqrt(393.75) * x * x * x * x * y * y * z * z - std::sqrt(393.75) * x * x * y * y * y * y * z * z + std::sqrt(1575.0) * x * x * y * y * z * z * z * z) + e_1 * (-std::sqrt(393.75) * x * x * x * x * y * y - std::sqrt(393.75) * x * x * x * x * z * z - std::sqrt(393.75) * x * x * y * y * y * y + std::sqrt(1575.0) * x * x * z * z * z * z - std::sqrt(393.75) * y * y * y * y * z * z + std::sqrt(1575.0) * y * y * z * z * z * z) + e_2 * (-std::sqrt(393.75) * x * x * x * x - std::sqrt(14175.0) * x * x * y * y + std::sqrt(3543.75) * x * x * z * z - std::sqrt(393.75) * y * y * y * y + std::sqrt(3543.75) * y * y * z * z + std::sqrt(1575.0) * z * z * z * z) + e_3 * (-std::sqrt(3543.75) * x * x - std::sqrt(3543.75) * y * y + std::sqrt(14175.0) * z * z);

        pc_23[k] = e_0 * (std::sqrt(9.84375) * x * x * x * x * x * y * y * z + std::sqrt(39.375) * x * x * x * y * y * y * y * z - std::sqrt(354.375) * x * x * x * y * y * z * z * z + std::sqrt(9.84375) * x * y * y * y * y * y * y * z - std::sqrt(354.375) * x * y * y * y * y * z * z * z + std::sqrt(630.0) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(9.84375) * x * x * x * x * x * z + std::sqrt(157.5) * x * x * x * y * y * z - std::sqrt(354.375) * x * x * x * z * z * z + std::sqrt(88.59375) * x * y * y * y * y * z + std::sqrt(6654.375) * x * y * y * z * z * z + std::sqrt(630.0) * x * z * z * z * z * z) + e_2 * (-std::sqrt(39.375) * x * x * x * z + std::sqrt(28704.375) * x * y * y * z + std::sqrt(19057.5) * x * z * z * z) + e_3 * (std::sqrt(35437.5) * x * z);
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

        pc_24[k] = e_0 * (std::sqrt(59.0625) * x * x * x * x * x * y * z * z + std::sqrt(236.25) * x * x * x * y * y * y * z * z - std::sqrt(420.0) * x * x * x * y * z * z * z * z + std::sqrt(59.0625) * x * y * y * y * y * y * z * z - std::sqrt(420.0) * x * y * y * y * z * z * z * z + std::sqrt(105.0) * x * y * z * z * z * z * z * z) + e_1 * (std::sqrt(59.0625) * x * x * x * x * x * y + std::sqrt(236.25) * x * x * x * y * y * y + std::sqrt(59.0625) * x * y * y * y * y * y + std::sqrt(945.0) * x * y * z * z * z * z) + e_2 * (std::sqrt(3780.0) * x * x * x * y + std::sqrt(3780.0) * x * y * y * y + std::sqrt(8505.0) * x * y * z * z) + e_3 * (std::sqrt(23625.0) * x * y);

        pc_25[k] = e_0 * (std::sqrt(9.84375) * x * x * x * x * x * x * y * z + std::sqrt(39.375) * x * x * x * x * y * y * y * z - std::sqrt(354.375) * x * x * x * x * y * z * z * z + std::sqrt(9.84375) * x * x * y * y * y * y * y * z - std::sqrt(354.375) * x * x * y * y * y * z * z * z + std::sqrt(630.0) * x * x * y * z * z * z * z * z) + e_1 * (std::sqrt(88.59375) * x * x * x * x * y * z + std::sqrt(157.5) * x * x * y * y * y * z + std::sqrt(6654.375) * x * x * y * z * z * z + std::sqrt(9.84375) * y * y * y * y * y * z - std::sqrt(354.375) * y * y * y * z * z * z + std::sqrt(630.0) * y * z * z * z * z * z) + e_2 * (std::sqrt(28704.375) * x * x * y * z - std::sqrt(39.375) * y * y * y * z + std::sqrt(19057.5) * y * z * z * z) + e_3 * (std::sqrt(35437.5) * y * z);

        pc_26[k] = e_0 * (-std::sqrt(98.4375) * x * x * x * x * x * y * z * z + std::sqrt(393.75) * x * x * x * y * z * z * z * z + std::sqrt(98.4375) * x * y * y * y * y * y * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(98.4375) * x * x * x * x * x * y + std::sqrt(393.75) * x * x * x * y * z * z + std::sqrt(98.4375) * x * y * y * y * y * y - std::sqrt(393.75) * x * y * y * y * z * z) + e_2 * (-std::sqrt(1575.0) * x * x * x * y + std::sqrt(1575.0) * x * y * y * y);

        pc_27[k] = e_0 * (-std::sqrt(16.40625) * x * x * x * x * x * x * y * z + std::sqrt(65.625) * x * x * x * x * y * y * y * z + std::sqrt(65.625) * x * x * x * x * y * z * z * z + std::sqrt(147.65625) * x * x * y * y * y * y * y * z - std::sqrt(590.625) * x * x * y * y * y * z * z * z) + e_1 * (-std::sqrt(147.65625) * x * x * x * x * y * z + std::sqrt(9450.0) * x * x * y * y * y * z - std::sqrt(590.625) * x * x * y * z * z * z + std::sqrt(147.65625) * y * y * y * y * y * z - std::sqrt(590.625) * y * y * y * z * z * z) + e_2 * (std::sqrt(5315.625) * x * x * y * z + std::sqrt(5315.625) * y * y * y * z - std::sqrt(2362.5) * y * z * z * z) + e_3 * (std::sqrt(2362.5) * y * z);
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

        pc_28[k] = e_0 * (std::sqrt(1.318359375) * x * x * x * x * x * x * y * y + std::sqrt(3.662109375) * x * x * x * x * y * y * y * y - std::sqrt(189.84375) * x * x * x * x * y * y * z * z + std::sqrt(0.146484375) * x * x * y * y * y * y * y * y - std::sqrt(84.375) * x * x * y * y * y * y * z * z + std::sqrt(84.375) * x * x * y * y * z * z * z * z - std::sqrt(0.146484375) * y * y * y * y * y * y * y * y + std::sqrt(21.09375) * y * y * y * y * y * y * z * z - std::sqrt(9.375) * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(1.318359375) * x * x * x * x * x * x + std::sqrt(222.802734375) * x * x * x * x * y * y - std::sqrt(189.84375) * x * x * x * x * z * z + std::sqrt(64.599609375) * x * x * y * y * y * y - std::sqrt(6834.375) * x * x * y * y * z * z + std::sqrt(84.375) * x * x * z * z * z * z - std::sqrt(32.958984375) * y * y * y * y * y * y + std::sqrt(1708.59375) * y * y * y * y * z * z - std::sqrt(84.375) * y * y * z * z * z * z) + e_2 * (std::sqrt(84.375) * x * x * x * x + std::sqrt(759.375) * x * x * y * y - std::sqrt(3037.5) * x * x * z * z - std::sqrt(337.5) * y * y * y * y + std::sqrt(3037.5) * y * y * z * z) + e_3 * (std::sqrt(84.375) * x * x - std::sqrt(84.375) * y * y);

        pc_29[k] = e_0 * (1.875 * x * x * x * x * x * y * y * z + 3.75 * x * x * x * y * y * y * y * z - 22.5 * x * x * x * y * y * z * z * z + 1.875 * x * y * y * y * y * y * y * z - 22.5 * x * y * y * y * y * z * z * z + 15.0 * x * y * y * z * z * z * z * z) + e_1 * (1.875 * x * x * x * x * x * z - 26.25 * x * x * x * y * y * z - 22.5 * x * x * x * z * z * z - 28.125 * x * y * y * y * y * z - 52.5 * x * y * y * z * z * z + 15.0 * x * z * z * z * z * z) + e_2 * (-37.5 * x * x * x * z - 202.5 * x * y * y * z + 15.0 * x * z * z * z) + e_3 * (-90.0 * x * z);

        pc_30[k] = e_0 * (-std::sqrt(0.087890625) * x * x * x * x * x * x * y * y - std::sqrt(0.791015625) * x * x * x * x * y * y * y * y + std::sqrt(22.5) * x * x * x * x * y * y * z * z - std::sqrt(0.791015625) * x * x * y * y * y * y * y * y + std::sqrt(90.0) * x * x * y * y * y * y * z * z - std::sqrt(275.625) * x * x * y * y * z * z * z * z - std::sqrt(0.087890625) * y * y * y * y * y * y * y * y + std::sqrt(22.5) * y * y * y * y * y * y * z * z - std::sqrt(275.625) * y * y * y * y * z * z * z * z + std::sqrt(90.0) * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(0.087890625) * x * x * x * x * x * x - std::sqrt(25.400390625) * x * x * x * x * y * y + std::sqrt(22.5) * x * x * x * x * z * z - std::sqrt(84.462890625) * x * x * y * y * y * y - std::sqrt(202.5) * x * x * y * y * z * z - std::sqrt(275.625) * x * x * z * z * z * z - std::sqrt(19.775390625) * y * y * y * y * y * y - std::sqrt(360.0) * y * y * y * y * z * z + std::sqrt(680.625) * y * y * z * z * z * z + std::sqrt(90.0) * z * z * z * z * z * z) + e_2 * (-std::sqrt(5.625) * x * x * x * x - std::sqrt(2480.625) * x * x * y * y - std::sqrt(1822.5) * x * x * z * z - std::sqrt(2250.0) * y * y * y * y + std::sqrt(202.5) * y * y * z * z + std::sqrt(5760.0) * z * z * z * z) + e_3 * (-std::sqrt(1265.625) * x * x - std::sqrt(11390.625) * y * y + std::sqrt(20250.0) * z * z);

        pc_31[k] = e_0 * (-std::sqrt(0.52734375) * x * x * x * x * x * x * y * z - std::sqrt(4.74609375) * x * x * x * x * y * y * y * z + std::sqrt(84.609375) * x * x * x * x * y * z * z * z - std::sqrt(4.74609375) * x * x * y * y * y * y * y * z + std::sqrt(338.4375) * x * x * y * y * y * z * z * z - std::sqrt(135.0) * x * x * y * z * z * z * z * z - std::sqrt(0.52734375) * y * y * y * y * y * y * y * z + std::sqrt(84.609375) * y * y * y * y * y * z * z * z - std::sqrt(135.0) * y * y * y * z * z * z * z * z + std::sqrt(15.0) * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(103.359375) * x * x * x * x * y * z + std::sqrt(413.4375) * x * x * y * y * y * z - std::sqrt(33.75) * x * x * y * z * z * z + std::sqrt(103.359375) * y * y * y * y * y * z - std::sqrt(33.75) * y * y * y * z * z * z + std::sqrt(1215.0) * y * z * z * z * z * z) + e_2 * (std::sqrt(2733.75) * x * x * y * z + std::sqrt(2733.75) * y * y * y * z + std::sqrt(26460.0) * y * z * z * z) + e_3 * (std::sqrt(54000.0) * y * z);
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

        pc_32[k] = e_0 * (-std::sqrt(0.087890625) * x * x * x * x * x * x * x * y - std::sqrt(0.791015625) * x * x * x * x * x * y * y * y + std::sqrt(22.5) * x * x * x * x * x * y * z * z - std::sqrt(0.791015625) * x * x * x * y * y * y * y * y + std::sqrt(90.0) * x * x * x * y * y * y * z * z - std::sqrt(275.625) * x * x * x * y * z * z * z * z - std::sqrt(0.087890625) * x * y * y * y * y * y * y * y + std::sqrt(22.5) * x * y * y * y * y * y * z * z - std::sqrt(275.625) * x * y * y * y * z * z * z * z + std::sqrt(90.0) * x * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(17.2265625) * x * x * x * x * x * y - std::sqrt(68.90625) * x * x * x * y * y * y - std::sqrt(562.5) * x * x * x * y * z * z - std::sqrt(17.2265625) * x * y * y * y * y * y - std::sqrt(562.5) * x * y * y * y * z * z + std::sqrt(1822.5) * x * y * z * z * z * z) + e_2 * (-std::sqrt(2030.625) * x * x * x * y - std::sqrt(2030.625) * x * y * y * y + std::sqrt(3240.0) * x * y * z * z) + e_3 * (-std::sqrt(5062.5) * x * y);

        pc_33[k] = e_0 * (0.9375 * x * x * x * x * x * x * y * z + 0.9375 * x * x * x * x * y * y * y * z - 11.25 * x * x * x * x * y * z * z * z - 0.9375 * x * x * y * y * y * y * y * z + 7.5 * x * x * y * z * z * z * z * z - 0.9375 * y * y * y * y * y * y * y * z + 11.25 * y * y * y * y * y * z * z * z - 7.5 * y * y * y * z * z * z * z * z) + e_1 * (-16.875 * x * x * x * x * y * z - 3.75 * x * x * y * y * y * z + 7.5 * x * x * y * z * z * z + 13.125 * y * y * y * y * y * z + 37.5 * y * y * y * z * z * z - 15.0 * y * z * z * z * z * z) + e_2 * (-45.0 * x * x * y * z + 120.0 * y * y * y * z - 15.0 * y * z * z * z) + e_3 * (90.0 * y * z);

        pc_34[k] = e_0 * (std::sqrt(0.146484375) * x * x * x * x * x * x * x * y - std::sqrt(0.146484375) * x * x * x * x * x * y * y * y - std::sqrt(21.09375) * x * x * x * x * x * y * z * z - std::sqrt(3.662109375) * x * x * x * y * y * y * y * y + std::sqrt(84.375) * x * x * x * y * y * y * z * z + std::sqrt(9.375) * x * x * x * y * z * z * z * z - std::sqrt(1.318359375) * x * y * y * y * y * y * y * y + std::sqrt(189.84375) * x * y * y * y * y * y * z * z - std::sqrt(84.375) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(5.2734375) * x * x * x * x * x * y - std::sqrt(189.84375) * x * x * x * y * y * y - std::sqrt(258.3984375) * x * y * y * y * y * y + std::sqrt(12150.0) * x * y * y * y * z * z - std::sqrt(337.5) * x * y * z * z * z * z) + e_2 * (-std::sqrt(84.375) * x * x * x * y - std::sqrt(2109.375) * x * y * y * y + std::sqrt(12150.0) * x * y * z * z) + e_3 * (-std::sqrt(337.5) * x * y);

        pc_35[k] = e_0 * (std::sqrt(19.775390625) * x * x * x * x * x * x * y * z + std::sqrt(54.931640625) * x * x * x * x * y * y * y * z - std::sqrt(140.625) * x * x * x * x * y * z * z * z + std::sqrt(2.197265625) * x * x * y * y * y * y * y * z - std::sqrt(62.5) * x * x * y * y * y * z * z * z + std::sqrt(5.625) * x * x * y * z * z * z * z * z - std::sqrt(2.197265625) * y * y * y * y * y * y * y * z + std::sqrt(15.625) * y * y * y * y * y * z * z * z - std::sqrt(0.625) * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(2847.65625) * x * x * x * x * y * z + std::sqrt(1265.625) * x * x * y * y * y * z - std::sqrt(5062.5) * x * x * y * z * z * z - std::sqrt(316.40625) * y * y * y * y * y * z + std::sqrt(562.5) * y * y * y * z * z * z) + e_2 * (std::sqrt(11390.625) * x * x * y * z - std::sqrt(1265.625) * y * y * y * z);
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

        pc_36[k] = e_0 * (std::sqrt(52.734375) * x * x * x * x * x * y * z * z + std::sqrt(210.9375) * x * x * x * y * y * y * z * z - std::sqrt(375.0) * x * x * x * y * z * z * z * z + std::sqrt(52.734375) * x * y * y * y * y * y * z * z - std::sqrt(375.0) * x * y * y * y * z * z * z * z + std::sqrt(15.0) * x * y * z * z * z * z * z * z) + e_1 * (std::sqrt(52.734375) * x * x * x * x * x * y + std::sqrt(210.9375) * x * x * x * y * y * y + std::sqrt(52.734375) * x * y * y * y * y * y - std::sqrt(3375.0) * x * y * z * z * z * z) + e_2 * (std::sqrt(3375.0) * x * x * x * y + std::sqrt(3375.0) * x * y * y * y - std::sqrt(30375.0) * x * y * z * z) + e_3 * (std::sqrt(3375.0) * x * y);

        pc_37[k] = e_0 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * y * z - std::sqrt(11.865234375) * x * x * x * x * y * y * y * z + std::sqrt(58.59375) * x * x * x * x * y * z * z * z - std::sqrt(11.865234375) * x * x * y * y * y * y * y * z + std::sqrt(234.375) * x * x * y * y * y * z * z * z - std::sqrt(165.375) * x * x * y * z * z * z * z * z - std::sqrt(1.318359375) * y * y * y * y * y * y * y * z + std::sqrt(58.59375) * y * y * y * y * y * z * z * z - std::sqrt(165.375) * y * y * y * z * z * z * z * z + std::sqrt(6.0) * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(21.09375) * x * x * x * x * y * z - std::sqrt(84.375) * x * x * y * y * y * z - std::sqrt(1350.0) * x * x * y * z * z * z - std::sqrt(21.09375) * y * y * y * y * y * z - std::sqrt(1350.0) * y * y * y * z * z * z) + e_2 * (-std::sqrt(6834.375) * x * x * y * z - std::sqrt(6834.375) * y * y * y * z - std::sqrt(5400.0) * y * z * z * z) + e_3 * (-std::sqrt(33750.0) * y * z);

        pc_38[k] = e_0 * (-2.8125 * x * x * x * x * x * x * z * z - 8.4375 * x * x * x * x * y * y * z * z + 9.375 * x * x * x * x * z * z * z * z - 8.4375 * x * x * y * y * y * y * z * z + 18.75 * x * x * y * y * z * z * z * z - 6.5 * x * x * z * z * z * z * z * z - 2.8125 * y * y * y * y * y * y * z * z + 9.375 * y * y * y * y * z * z * z * z - 6.5 * y * y * z * z * z * z * z * z + z * z * z * z * z * z * z * z) + e_1 * (-2.8125 * x * x * x * x * x * x - 8.4375 * x * x * x * x * y * y + 5.625 * x * x * x * x * z * z - 8.4375 * x * x * y * y * y * y + 11.25 * x * x * y * y * z * z - 22.5 * x * x * z * z * z * z - 2.8125 * y * y * y * y * y * y + 5.625 * y * y * y * y * z * z - 22.5 * y * y * z * z * z * z + 15.0 * z * z * z * z * z * z) + e_2 * (-22.5 * x * x * x * x - 45.0 * x * x * y * y - 45.0 * x * x * z * z - 22.5 * y * y * y * y - 45.0 * y * y * z * z + 90.0 * z * z * z * z) + e_3 * (-75.0 * x * x - 75.0 * y * y + 150.0 * z * z);

        pc_39[k] = e_0 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * x * z - std::sqrt(11.865234375) * x * x * x * x * x * y * y * z + std::sqrt(58.59375) * x * x * x * x * x * z * z * z - std::sqrt(11.865234375) * x * x * x * y * y * y * y * z + std::sqrt(234.375) * x * x * x * y * y * z * z * z - std::sqrt(165.375) * x * x * x * z * z * z * z * z - std::sqrt(1.318359375) * x * y * y * y * y * y * y * z + std::sqrt(58.59375) * x * y * y * y * y * z * z * z - std::sqrt(165.375) * x * y * y * z * z * z * z * z + std::sqrt(6.0) * x * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(21.09375) * x * x * x * x * x * z - std::sqrt(84.375) * x * x * x * y * y * z - std::sqrt(1350.0) * x * x * x * z * z * z - std::sqrt(21.09375) * x * y * y * y * y * z - std::sqrt(1350.0) * x * y * y * z * z * z) + e_2 * (-std::sqrt(6834.375) * x * x * x * z - std::sqrt(6834.375) * x * y * y * z - std::sqrt(5400.0) * x * z * z * z) + e_3 * (-std::sqrt(33750.0) * x * z);
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

        pc_40[k] = e_0 * (std::sqrt(13.18359375) * x * x * x * x * x * x * z * z + std::sqrt(13.18359375) * x * x * x * x * y * y * z * z - std::sqrt(93.75) * x * x * x * x * z * z * z * z - std::sqrt(13.18359375) * x * x * y * y * y * y * z * z + std::sqrt(3.75) * x * x * z * z * z * z * z * z - std::sqrt(13.18359375) * y * y * y * y * y * y * z * z + std::sqrt(93.75) * y * y * y * y * z * z * z * z - std::sqrt(3.75) * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(13.18359375) * x * x * x * x * x * x + std::sqrt(13.18359375) * x * x * x * x * y * y - std::sqrt(13.18359375) * x * x * y * y * y * y - std::sqrt(843.75) * x * x * z * z * z * z - std::sqrt(13.18359375) * y * y * y * y * y * y + std::sqrt(843.75) * y * y * z * z * z * z) + e_2 * (std::sqrt(843.75) * x * x * x * x - std::sqrt(7593.75) * x * x * z * z - std::sqrt(843.75) * y * y * y * y + std::sqrt(7593.75) * y * y * z * z) + e_3 * (std::sqrt(843.75) * x * x - std::sqrt(843.75) * y * y);

        pc_41[k] = e_0 * (std::sqrt(2.197265625) * x * x * x * x * x * x * x * z - std::sqrt(2.197265625) * x * x * x * x * x * y * y * z - std::sqrt(15.625) * x * x * x * x * x * z * z * z - std::sqrt(54.931640625) * x * x * x * y * y * y * y * z + std::sqrt(62.5) * x * x * x * y * y * z * z * z + std::sqrt(0.625) * x * x * x * z * z * z * z * z - std::sqrt(19.775390625) * x * y * y * y * y * y * y * z + std::sqrt(140.625) * x * y * y * y * y * z * z * z - std::sqrt(5.625) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(316.40625) * x * x * x * x * x * z - std::sqrt(1265.625) * x * x * x * y * y * z - std::sqrt(562.5) * x * x * x * z * z * z - std::sqrt(2847.65625) * x * y * y * y * y * z + std::sqrt(5062.5) * x * y * y * z * z * z) + e_2 * (std::sqrt(1265.625) * x * x * x * z - std::sqrt(11390.625) * x * y * y * z);

        pc_42[k] = e_0 * (std::sqrt(1.318359375) * x * x * x * x * x * x * x * y + std::sqrt(3.662109375) * x * x * x * x * x * y * y * y - std::sqrt(189.84375) * x * x * x * x * x * y * z * z + std::sqrt(0.146484375) * x * x * x * y * y * y * y * y - std::sqrt(84.375) * x * x * x * y * y * y * z * z + std::sqrt(84.375) * x * x * x * y * z * z * z * z - std::sqrt(0.146484375) * x * y * y * y * y * y * y * y + std::sqrt(21.09375) * x * y * y * y * y * y * z * z - std::sqrt(9.375) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(258.3984375) * x * x * x * x * x * y + std::sqrt(189.84375) * x * x * x * y * y * y - std::sqrt(12150.0) * x * x * x * y * z * z - std::sqrt(5.2734375) * x * y * y * y * y * y + std::sqrt(337.5) * x * y * z * z * z * z) + e_2 * (std::sqrt(2109.375) * x * x * x * y + std::sqrt(84.375) * x * y * y * y - std::sqrt(12150.0) * x * y * z * z) + e_3 * (std::sqrt(337.5) * x * y);

        pc_43[k] = e_0 * (1.875 * x * x * x * x * x * x * y * z + 3.75 * x * x * x * x * y * y * y * z - 22.5 * x * x * x * x * y * z * z * z + 1.875 * x * x * y * y * y * y * y * z - 22.5 * x * x * y * y * y * z * z * z + 15.0 * x * x * y * z * z * z * z * z) + e_1 * (-28.125 * x * x * x * x * y * z - 26.25 * x * x * y * y * y * z - 52.5 * x * x * y * z * z * z + 1.875 * y * y * y * y * y * z - 22.5 * y * y * y * z * z * z + 15.0 * y * z * z * z * z * z) + e_2 * (-202.5 * x * x * y * z - 37.5 * y * y * y * z + 15.0 * y * z * z * z) + e_3 * (-90.0 * y * z);
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

        pc_44[k] = e_0 * (-std::sqrt(0.087890625) * x * x * x * x * x * x * x * y - std::sqrt(0.791015625) * x * x * x * x * x * y * y * y + std::sqrt(22.5) * x * x * x * x * x * y * z * z - std::sqrt(0.791015625) * x * x * x * y * y * y * y * y + std::sqrt(90.0) * x * x * x * y * y * y * z * z - std::sqrt(275.625) * x * x * x * y * z * z * z * z - std::sqrt(0.087890625) * x * y * y * y * y * y * y * y + std::sqrt(22.5) * x * y * y * y * y * y * z * z - std::sqrt(275.625) * x * y * y * y * z * z * z * z + std::sqrt(90.0) * x * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(17.2265625) * x * x * x * x * x * y - std::sqrt(68.90625) * x * x * x * y * y * y - std::sqrt(562.5) * x * x * x * y * z * z - std::sqrt(17.2265625) * x * y * y * y * y * y - std::sqrt(562.5) * x * y * y * y * z * z + std::sqrt(1822.5) * x * y * z * z * z * z) + e_2 * (-std::sqrt(2030.625) * x * x * x * y - std::sqrt(2030.625) * x * y * y * y + std::sqrt(3240.0) * x * y * z * z) + e_3 * (-std::sqrt(5062.5) * x * y);

        pc_45[k] = e_0 * (-std::sqrt(0.52734375) * x * x * x * x * x * x * x * z - std::sqrt(4.74609375) * x * x * x * x * x * y * y * z + std::sqrt(84.609375) * x * x * x * x * x * z * z * z - std::sqrt(4.74609375) * x * x * x * y * y * y * y * z + std::sqrt(338.4375) * x * x * x * y * y * z * z * z - std::sqrt(135.0) * x * x * x * z * z * z * z * z - std::sqrt(0.52734375) * x * y * y * y * y * y * y * z + std::sqrt(84.609375) * x * y * y * y * y * z * z * z - std::sqrt(135.0) * x * y * y * z * z * z * z * z + std::sqrt(15.0) * x * z * z * z * z * z * z * z) + e_1 * (std::sqrt(103.359375) * x * x * x * x * x * z + std::sqrt(413.4375) * x * x * x * y * y * z - std::sqrt(33.75) * x * x * x * z * z * z + std::sqrt(103.359375) * x * y * y * y * y * z - std::sqrt(33.75) * x * y * y * z * z * z + std::sqrt(1215.0) * x * z * z * z * z * z) + e_2 * (std::sqrt(2733.75) * x * x * x * z + std::sqrt(2733.75) * x * y * y * z + std::sqrt(26460.0) * x * z * z * z) + e_3 * (std::sqrt(54000.0) * x * z);

        pc_46[k] = e_0 * (-std::sqrt(0.087890625) * x * x * x * x * x * x * x * x - std::sqrt(0.791015625) * x * x * x * x * x * x * y * y + std::sqrt(22.5) * x * x * x * x * x * x * z * z - std::sqrt(0.791015625) * x * x * x * x * y * y * y * y + std::sqrt(90.0) * x * x * x * x * y * y * z * z - std::sqrt(275.625) * x * x * x * x * z * z * z * z - std::sqrt(0.087890625) * x * x * y * y * y * y * y * y + std::sqrt(22.5) * x * x * y * y * y * y * z * z - std::sqrt(275.625) * x * x * y * y * z * z * z * z + std::sqrt(90.0) * x * x * z * z * z * z * z * z) + e_1 * (-std::sqrt(19.775390625) * x * x * x * x * x * x - std::sqrt(84.462890625) * x * x * x * x * y * y - std::sqrt(360.0) * x * x * x * x * z * z - std::sqrt(25.400390625) * x * x * y * y * y * y - std::sqrt(202.5) * x * x * y * y * z * z + std::sqrt(680.625) * x * x * z * z * z * z - std::sqrt(0.087890625) * y * y * y * y * y * y + std::sqrt(22.5) * y * y * y * y * z * z - std::sqrt(275.625) * y * y * z * z * z * z + std::sqrt(90.0) * z * z * z * z * z * z) + e_2 * (-std::sqrt(2250.0) * x * x * x * x - std::sqrt(2480.625) * x * x * y * y + std::sqrt(202.5) * x * x * z * z - std::sqrt(5.625) * y * y * y * y - std::sqrt(1822.5) * y * y * z * z + std::sqrt(5760.0) * z * z * z * z) + e_3 * (-std::sqrt(11390.625) * x * x - std::sqrt(1265.625) * y * y + std::sqrt(20250.0) * z * z);

        pc_47[k] = e_0 * (0.9375 * x * x * x * x * x * x * x * z + 0.9375 * x * x * x * x * x * y * y * z - 11.25 * x * x * x * x * x * z * z * z - 0.9375 * x * x * x * y * y * y * y * z + 7.5 * x * x * x * z * z * z * z * z - 0.9375 * x * y * y * y * y * y * y * z + 11.25 * x * y * y * y * y * z * z * z - 7.5 * x * y * y * z * z * z * z * z) + e_1 * (-13.125 * x * x * x * x * x * z + 3.75 * x * x * x * y * y * z - 37.5 * x * x * x * z * z * z + 16.875 * x * y * y * y * y * z - 7.5 * x * y * y * z * z * z + 15.0 * x * z * z * z * z * z) + e_2 * (-120.0 * x * x * x * z + 45.0 * x * y * y * z + 15.0 * x * z * z * z) + e_3 * (-90.0 * x * z);
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

        pc_48[k] = e_0 * (std::sqrt(0.146484375) * x * x * x * x * x * x * x * x - std::sqrt(0.146484375) * x * x * x * x * x * x * y * y - std::sqrt(21.09375) * x * x * x * x * x * x * z * z - std::sqrt(3.662109375) * x * x * x * x * y * y * y * y + std::sqrt(84.375) * x * x * x * x * y * y * z * z + std::sqrt(9.375) * x * x * x * x * z * z * z * z - std::sqrt(1.318359375) * x * x * y * y * y * y * y * y + std::sqrt(189.84375) * x * x * y * y * y * y * z * z - std::sqrt(84.375) * x * x * y * y * z * z * z * z) + e_1 * (std::sqrt(32.958984375) * x * x * x * x * x * x - std::sqrt(64.599609375) * x * x * x * x * y * y - std::sqrt(1708.59375) * x * x * x * x * z * z - std::sqrt(222.802734375) * x * x * y * y * y * y + std::sqrt(6834.375) * x * x * y * y * z * z + std::sqrt(84.375) * x * x * z * z * z * z - std::sqrt(1.318359375) * y * y * y * y * y * y + std::sqrt(189.84375) * y * y * y * y * z * z - std::sqrt(84.375) * y * y * z * z * z * z) + e_2 * (std::sqrt(337.5) * x * x * x * x - std::sqrt(759.375) * x * x * y * y - std::sqrt(3037.5) * x * x * z * z - std::sqrt(84.375) * y * y * y * y + std::sqrt(3037.5) * y * y * z * z) + e_3 * (std::sqrt(84.375) * x * x - std::sqrt(84.375) * y * y);

        pc_49[k] = e_0 * (-std::sqrt(36.9140625) * x * x * x * x * x * x * y * z + std::sqrt(4.1015625) * x * x * x * x * y * y * y * z + std::sqrt(147.65625) * x * x * x * x * y * z * z * z + std::sqrt(36.9140625) * x * x * y * y * y * y * y * z - std::sqrt(262.5) * x * x * y * y * y * z * z * z - std::sqrt(4.1015625) * y * y * y * y * y * y * y * z + std::sqrt(16.40625) * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(2362.5) * x * x * x * x * y * z + std::sqrt(590.625) * x * x * y * y * y * z + std::sqrt(590.625) * x * x * y * z * z * z - std::sqrt(590.625) * y * y * y * y * y * z + std::sqrt(590.625) * y * y * y * z * z * z) + e_2 * (-std::sqrt(5315.625) * x * x * y * z - std::sqrt(5315.625) * y * y * y * z + std::sqrt(2362.5) * y * z * z * z) + e_3 * (-std::sqrt(2362.5) * y * z);

        pc_50[k] = e_0 * (-std::sqrt(98.4375) * x * x * x * x * x * y * z * z + std::sqrt(393.75) * x * x * x * y * z * z * z * z + std::sqrt(98.4375) * x * y * y * y * y * y * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(98.4375) * x * x * x * x * x * y + std::sqrt(393.75) * x * x * x * y * z * z + std::sqrt(98.4375) * x * y * y * y * y * y - std::sqrt(393.75) * x * y * y * y * z * z) + e_2 * (-std::sqrt(1575.0) * x * x * x * y + std::sqrt(1575.0) * x * y * y * y);

        pc_51[k] = e_0 * (std::sqrt(2.4609375) * x * x * x * x * x * x * y * z + std::sqrt(2.4609375) * x * x * x * x * y * y * y * z - std::sqrt(88.59375) * x * x * x * x * y * z * z * z - std::sqrt(2.4609375) * x * x * y * y * y * y * y * z + std::sqrt(157.5) * x * x * y * z * z * z * z * z - std::sqrt(2.4609375) * y * y * y * y * y * y * y * z + std::sqrt(88.59375) * y * y * y * y * y * z * z * z - std::sqrt(157.5) * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(39.375) * x * x * y * y * y * z + std::sqrt(4764.375) * x * x * y * z * z * z - std::sqrt(39.375) * y * y * y * y * y * z - std::sqrt(984.375) * y * y * y * z * z * z - std::sqrt(630.0) * y * z * z * z * z * z) + e_2 * (std::sqrt(8859.375) * x * x * y * z - std::sqrt(6654.375) * y * y * y * z - std::sqrt(19057.5) * y * z * z * z) + e_3 * (-std::sqrt(35437.5) * y * z);
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

        pc_52[k] = e_0 * (std::sqrt(14.765625) * x * x * x * x * x * x * z * z + std::sqrt(14.765625) * x * x * x * x * y * y * z * z - std::sqrt(105.0) * x * x * x * x * z * z * z * z - std::sqrt(14.765625) * x * x * y * y * y * y * z * z + std::sqrt(26.25) * x * x * z * z * z * z * z * z - std::sqrt(14.765625) * y * y * y * y * y * y * z * z + std::sqrt(105.0) * y * y * y * y * z * z * z * z - std::sqrt(26.25) * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(14.765625) * x * x * x * x * x * x + std::sqrt(14.765625) * x * x * x * x * y * y - std::sqrt(14.765625) * x * x * y * y * y * y + std::sqrt(236.25) * x * x * z * z * z * z - std::sqrt(14.765625) * y * y * y * y * y * y - std::sqrt(236.25) * y * y * z * z * z * z) + e_2 * (std::sqrt(945.0) * x * x * x * x + std::sqrt(2126.25) * x * x * z * z - std::sqrt(945.0) * y * y * y * y - std::sqrt(2126.25) * y * y * z * z) + e_3 * (std::sqrt(5906.25) * x * x - std::sqrt(5906.25) * y * y);

        pc_53[k] = e_0 * (std::sqrt(2.4609375) * x * x * x * x * x * x * x * z + std::sqrt(2.4609375) * x * x * x * x * x * y * y * z - std::sqrt(88.59375) * x * x * x * x * x * z * z * z - std::sqrt(2.4609375) * x * x * x * y * y * y * y * z + std::sqrt(157.5) * x * x * x * z * z * z * z * z - std::sqrt(2.4609375) * x * y * y * y * y * y * y * z + std::sqrt(88.59375) * x * y * y * y * y * z * z * z - std::sqrt(157.5) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(39.375) * x * x * x * x * x * z + std::sqrt(39.375) * x * x * x * y * y * z + std::sqrt(984.375) * x * x * x * z * z * z - std::sqrt(4764.375) * x * y * y * z * z * z + std::sqrt(630.0) * x * z * z * z * z * z) + e_2 * (std::sqrt(6654.375) * x * x * x * z - std::sqrt(8859.375) * x * y * y * z + std::sqrt(19057.5) * x * z * z * z) + e_3 * (std::sqrt(35437.5) * x * z);

        pc_54[k] = e_0 * (-std::sqrt(24.609375) * x * x * x * x * x * x * z * z + std::sqrt(24.609375) * x * x * x * x * y * y * z * z + std::sqrt(98.4375) * x * x * x * x * z * z * z * z + std::sqrt(24.609375) * x * x * y * y * y * y * z * z - std::sqrt(393.75) * x * x * y * y * z * z * z * z - std::sqrt(24.609375) * y * y * y * y * y * y * z * z + std::sqrt(98.4375) * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(24.609375) * x * x * x * x * x * x + std::sqrt(24.609375) * x * x * x * x * y * y - std::sqrt(98.4375) * x * x * x * x * z * z + std::sqrt(24.609375) * x * x * y * y * y * y - std::sqrt(3543.75) * x * x * y * y * z * z + std::sqrt(1575.0) * x * x * z * z * z * z - std::sqrt(24.609375) * y * y * y * y * y * y - std::sqrt(98.4375) * y * y * y * y * z * z + std::sqrt(1575.0) * y * y * z * z * z * z) + e_2 * (-std::sqrt(1575.0) * x * x * x * x + std::sqrt(3543.75) * x * x * z * z - std::sqrt(1575.0) * y * y * y * y + std::sqrt(3543.75) * y * y * z * z + std::sqrt(1575.0) * z * z * z * z) + e_3 * (-std::sqrt(3543.75) * x * x - std::sqrt(3543.75) * y * y + std::sqrt(14175.0) * z * z);

        pc_55[k] = e_0 * (-std::sqrt(4.1015625) * x * x * x * x * x * x * x * z + std::sqrt(36.9140625) * x * x * x * x * x * y * y * z + std::sqrt(16.40625) * x * x * x * x * x * z * z * z + std::sqrt(4.1015625) * x * x * x * y * y * y * y * z - std::sqrt(262.5) * x * x * x * y * y * z * z * z - std::sqrt(36.9140625) * x * y * y * y * y * y * y * z + std::sqrt(147.65625) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(590.625) * x * x * x * x * x * z + std::sqrt(590.625) * x * x * x * y * y * z + std::sqrt(590.625) * x * x * x * z * z * z - std::sqrt(2362.5) * x * y * y * y * y * z + std::sqrt(590.625) * x * y * y * z * z * z) + e_2 * (-std::sqrt(5315.625) * x * x * x * z - std::sqrt(5315.625) * x * y * y * z + std::sqrt(2362.5) * x * z * z * z) + e_3 * (-std::sqrt(2362.5) * x * z);
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

        pc_56[k] = e_0 * (-std::sqrt(1.5380859375) * x * x * x * x * x * x * x * y + std::sqrt(8.3740234375) * x * x * x * x * x * y * y * y + std::sqrt(98.4375) * x * x * x * x * x * y * z * z + std::sqrt(8.3740234375) * x * x * x * y * y * y * y * y - std::sqrt(1093.75) * x * x * x * y * y * y * z * z - std::sqrt(1.5380859375) * x * y * y * y * y * y * y * y + std::sqrt(98.4375) * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(55.37109375) * x * x * x * x * x * y + std::sqrt(615.234375) * x * x * x * y * y * y - std::sqrt(55.37109375) * x * y * y * y * y * y);

        pc_57[k] = e_0 * (-std::sqrt(4.1015625) * x * x * x * x * x * x * y * z + std::sqrt(16.40625) * x * x * x * x * y * y * y * z + std::sqrt(262.5) * x * x * x * x * y * z * z * z + std::sqrt(36.9140625) * x * x * y * y * y * y * y * z - std::sqrt(2362.5) * x * x * y * y * y * z * z * z) + e_1 * (std::sqrt(922.8515625) * x * x * x * x * y * z - std::sqrt(3691.40625) * x * x * y * y * y * z - std::sqrt(2362.5) * x * x * y * z * z * z + std::sqrt(36.9140625) * y * y * y * y * y * z - std::sqrt(2362.5) * y * y * y * z * z * z) + e_2 * (-std::sqrt(5315.625) * x * x * y * z - std::sqrt(5315.625) * y * y * y * z - std::sqrt(9450.0) * y * z * z * z) + e_3 * (-std::sqrt(37800.0) * y * z);

        pc_58[k] = e_0 * (std::sqrt(0.1025390625) * x * x * x * x * x * x * x * y - std::sqrt(0.1025390625) * x * x * x * x * x * y * y * y - std::sqrt(14.765625) * x * x * x * x * x * y * z * z - std::sqrt(2.5634765625) * x * x * x * y * y * y * y * y + std::sqrt(59.0625) * x * x * x * y * y * y * z * z + std::sqrt(105.0) * x * x * x * y * z * z * z * z - std::sqrt(0.9228515625) * x * y * y * y * y * y * y * y + std::sqrt(132.890625) * x * y * y * y * y * y * z * z - std::sqrt(945.0) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(3.69140625) * x * x * x * x * x * y - std::sqrt(132.890625) * x * x * x * y * y * y + std::sqrt(2126.25) * x * x * x * y * z * z - std::sqrt(180.87890625) * x * y * y * y * y * y - std::sqrt(2126.25) * x * y * y * y * z * z - std::sqrt(3780.0) * x * y * z * z * z * z) + e_2 * (std::sqrt(236.25) * x * x * x * y - std::sqrt(11576.25) * x * y * y * y - std::sqrt(34020.0) * x * y * z * z) + e_3 * (-std::sqrt(23625.0) * x * y);

        pc_59[k] = e_0 * (std::sqrt(0.615234375) * x * x * x * x * x * x * x * z - std::sqrt(0.615234375) * x * x * x * x * x * y * y * z - std::sqrt(46.2109375) * x * x * x * x * x * z * z * z - std::sqrt(15.380859375) * x * x * x * y * y * y * y * z + std::sqrt(184.84375) * x * x * x * y * y * z * z * z + std::sqrt(17.5) * x * x * x * z * z * z * z * z - std::sqrt(5.537109375) * x * y * y * y * y * y * y * z + std::sqrt(415.8984375) * x * y * y * y * y * z * z * z - std::sqrt(157.5) * x * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(22.1484375) * x * x * x * x * x * z + std::sqrt(88.59375) * x * x * x * y * y * z - std::sqrt(157.5) * x * x * x * z * z * z + std::sqrt(199.3359375) * x * y * y * y * y * z + std::sqrt(1417.5) * x * y * y * z * z * z) + e_2 * (-std::sqrt(1417.5) * x * x * x * z + std::sqrt(12757.5) * x * y * y * z);
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

        pc_60[k] = e_0 * (std::sqrt(0.1025390625) * x * x * x * x * x * x * x * x - std::sqrt(0.1025390625) * x * x * x * x * x * x * y * y - std::sqrt(14.765625) * x * x * x * x * x * x * z * z - std::sqrt(2.5634765625) * x * x * x * x * y * y * y * y + std::sqrt(59.0625) * x * x * x * x * y * y * z * z + std::sqrt(105.0) * x * x * x * x * z * z * z * z - std::sqrt(0.9228515625) * x * x * y * y * y * y * y * y + std::sqrt(132.890625) * x * x * y * y * y * y * z * z - std::sqrt(945.0) * x * x * y * y * z * z * z * z) + e_1 * (std::sqrt(23.0712890625) * x * x * x * x * x * x - std::sqrt(45.2197265625) * x * x * x * x * y * y + std::sqrt(132.890625) * x * x * x * x * z * z - std::sqrt(155.9619140625) * x * x * y * y * y * y - std::sqrt(4784.0625) * x * x * y * y * z * z + std::sqrt(945.0) * x * x * z * z * z * z - std::sqrt(0.9228515625) * y * y * y * y * y * y + std::sqrt(132.890625) * y * y * y * y * z * z - std::sqrt(945.0) * y * y * z * z * z * z) + e_2 * (std::sqrt(1476.5625) * x * x * x * x - std::sqrt(8505.0) * x * x * y * y + std::sqrt(8505.0) * x * x * z * z - std::sqrt(59.0625) * y * y * y * y - std::sqrt(8505.0) * y * y * z * z) + e_3 * (std::sqrt(5906.25) * x * x - std::sqrt(5906.25) * y * y);

        pc_61[k] = e_0 * (-std::sqrt(1.025390625) * x * x * x * x * x * x * x * z + std::sqrt(9.228515625) * x * x * x * x * x * y * y * z + std::sqrt(65.625) * x * x * x * x * x * z * z * z + std::sqrt(1.025390625) * x * x * x * y * y * y * y * z - std::sqrt(1050.0) * x * x * x * y * y * z * z * z - std::sqrt(9.228515625) * x * y * y * y * y * y * y * z + std::sqrt(590.625) * x * y * y * y * y * z * z * z) + e_1 * (std::sqrt(36.9140625) * x * x * x * x * x * z - std::sqrt(3691.40625) * x * x * x * y * y * z + std::sqrt(2362.5) * x * x * x * z * z * z + std::sqrt(922.8515625) * x * y * y * y * y * z + std::sqrt(2362.5) * x * y * y * z * z * z) + e_2 * (std::sqrt(5315.625) * x * x * x * z + std::sqrt(5315.625) * x * y * y * z + std::sqrt(9450.0) * x * z * z * z) + e_3 * (std::sqrt(37800.0) * x * z);

        pc_62[k] = e_0 * (-std::sqrt(0.1708984375) * x * x * x * x * x * x * x * x + std::sqrt(4.2724609375) * x * x * x * x * x * x * y * y + std::sqrt(10.9375) * x * x * x * x * x * x * z * z - std::sqrt(1.5380859375) * x * x * x * x * y * y * y * y - std::sqrt(393.75) * x * x * x * x * y * y * z * z - std::sqrt(13.8427734375) * x * x * y * y * y * y * y * y + std::sqrt(885.9375) * x * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(38.4521484375) * x * x * x * x * x * x + std::sqrt(13.8427734375) * x * x * x * x * y * y + std::sqrt(885.9375) * x * x * x * x * z * z - std::sqrt(1121.2646484375) * x * x * y * y * y * y + std::sqrt(3543.75) * x * x * y * y * z * z - std::sqrt(13.8427734375) * y * y * y * y * y * y + std::sqrt(885.9375) * y * y * y * y * z * z) + e_2 * (-std::sqrt(885.9375) * x * x * x * x - std::sqrt(3543.75) * x * x * y * y + std::sqrt(14175.0) * x * x * z * z - std::sqrt(885.9375) * y * y * y * y + std::sqrt(14175.0) * y * y * z * z) + e_3 * (-std::sqrt(1575.0) * x * x - std::sqrt(1575.0) * y * y + std::sqrt(6300.0) * z * z);

        pc_63[k] = e_0 * (std::sqrt(27.685546875) * x * x * x * x * x * x * y * z - std::sqrt(1110.498046875) * x * x * x * x * y * y * y * z + std::sqrt(249.169921875) * x * x * y * y * y * y * y * z - std::sqrt(3.076171875) * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(442.96875) * x * x * x * x * y * z - std::sqrt(1771.875) * x * x * y * y * y * z - std::sqrt(442.96875) * y * y * y * y * y * z) + e_2 * (-std::sqrt(15946.875) * x * x * y * z - std::sqrt(15946.875) * y * y * y * z) + e_3 * (-std::sqrt(28350.0) * y * z);
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

        pc_64[k] = e_0 * (std::sqrt(73.828125) * x * x * x * x * x * y * z * z - std::sqrt(2657.8125) * x * x * x * y * y * y * z * z + std::sqrt(73.828125) * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(73.828125) * x * x * x * x * x * y - std::sqrt(2657.8125) * x * x * x * y * y * y - std::sqrt(4725.0) * x * x * x * y * z * z + std::sqrt(73.828125) * x * y * y * y * y * y - std::sqrt(4725.0) * x * y * y * y * z * z) + e_2 * (-std::sqrt(4725.0) * x * x * x * y - std::sqrt(4725.0) * x * y * y * y - std::sqrt(42525.0) * x * y * z * z) + e_3 * (-std::sqrt(42525.0) * x * y);

        pc_65[k] = e_0 * (-std::sqrt(1.845703125) * x * x * x * x * x * x * y * z + std::sqrt(46.142578125) * x * x * x * x * y * y * y * z + std::sqrt(29.53125) * x * x * x * x * y * z * z * z + std::sqrt(46.142578125) * x * x * y * y * y * y * y * z - std::sqrt(1063.125) * x * x * y * y * y * z * z * z - std::sqrt(1.845703125) * y * y * y * y * y * y * y * z + std::sqrt(29.53125) * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(265.78125) * x * x * x * x * y * z + std::sqrt(118.125) * x * x * y * y * y * z - std::sqrt(4252.5) * x * x * y * z * z * z - std::sqrt(29.53125) * y * y * y * y * y * z + std::sqrt(472.5) * y * y * y * z * z * z) + e_2 * (-std::sqrt(1063.125) * x * x * y * z + std::sqrt(118.125) * y * y * y * z);

        pc_66[k] = e_0 * (-std::sqrt(11.07421875) * x * x * x * x * x * x * z * z + std::sqrt(276.85546875) * x * x * x * x * y * y * z * z + std::sqrt(4.921875) * x * x * x * x * z * z * z * z + std::sqrt(276.85546875) * x * x * y * y * y * y * z * z - std::sqrt(177.1875) * x * x * y * y * z * z * z * z - std::sqrt(11.07421875) * y * y * y * y * y * y * z * z + std::sqrt(4.921875) * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(11.07421875) * x * x * x * x * x * x + std::sqrt(276.85546875) * x * x * x * x * y * y - std::sqrt(398.671875) * x * x * x * x * z * z + std::sqrt(276.85546875) * x * x * y * y * y * y + std::sqrt(14352.1875) * x * x * y * y * z * z - std::sqrt(11.07421875) * y * y * y * y * y * y - std::sqrt(398.671875) * y * y * y * y * z * z) + e_2 * (-std::sqrt(708.75) * x * x * x * x + std::sqrt(25515.0) * x * x * y * y - std::sqrt(708.75) * y * y * y * y);

        pc_67[k] = e_0 * (-std::sqrt(1.845703125) * x * x * x * x * x * x * x * z + std::sqrt(46.142578125) * x * x * x * x * x * y * y * z + std::sqrt(29.53125) * x * x * x * x * x * z * z * z + std::sqrt(46.142578125) * x * x * x * y * y * y * y * z - std::sqrt(1063.125) * x * x * x * y * y * z * z * z - std::sqrt(1.845703125) * x * y * y * y * y * y * y * z + std::sqrt(29.53125) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(29.53125) * x * x * x * x * x * z + std::sqrt(118.125) * x * x * x * y * y * z + std::sqrt(472.5) * x * x * x * z * z * z + std::sqrt(265.78125) * x * y * y * y * y * z - std::sqrt(4252.5) * x * y * y * z * z * z) + e_2 * (std::sqrt(118.125) * x * x * x * z - std::sqrt(1063.125) * x * y * y * z);
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

        pc_68[k] = e_0 * (std::sqrt(18.45703125) * x * x * x * x * x * x * z * z - std::sqrt(904.39453125) * x * x * x * x * y * y * z * z + std::sqrt(904.39453125) * x * x * y * y * y * y * z * z - std::sqrt(18.45703125) * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(18.45703125) * x * x * x * x * x * x - std::sqrt(904.39453125) * x * x * x * x * y * y + std::sqrt(1181.25) * x * x * x * x * z * z + std::sqrt(904.39453125) * x * x * y * y * y * y - std::sqrt(18.45703125) * y * y * y * y * y * y - std::sqrt(1181.25) * y * y * y * y * z * z) + e_2 * (std::sqrt(1181.25) * x * x * x * x + std::sqrt(10631.25) * x * x * z * z - std::sqrt(1181.25) * y * y * y * y - std::sqrt(10631.25) * y * y * z * z) + e_3 * (std::sqrt(10631.25) * x * x - std::sqrt(10631.25) * y * y);

        pc_69[k] = e_0 * (std::sqrt(3.076171875) * x * x * x * x * x * x * x * z - std::sqrt(249.169921875) * x * x * x * x * x * y * y * z + std::sqrt(1110.498046875) * x * x * x * y * y * y * y * z - std::sqrt(27.685546875) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(442.96875) * x * x * x * x * x * z + std::sqrt(1771.875) * x * x * x * y * y * z + std::sqrt(442.96875) * x * y * y * y * y * z) + e_2 * (std::sqrt(15946.875) * x * x * x * z + std::sqrt(15946.875) * x * y * y * z) + e_3 * (std::sqrt(28350.0) * x * z);

        pc_70[k] = e_0 * (std::sqrt(2.7685546875) * x * x * x * x * x * x * x * y - std::sqrt(295.6201171875) * x * x * x * x * x * y * y * y + std::sqrt(192.2607421875) * x * x * x * y * y * y * y * y - std::sqrt(7.6904296875) * x * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(276.85546875) * x * x * x * x * x * y - std::sqrt(1107.421875) * x * x * x * y * y * y - std::sqrt(276.85546875) * x * y * y * y * y * y) + e_2 * (-std::sqrt(17718.75) * x * x * x * y - std::sqrt(17718.75) * x * y * y * y) + e_3 * (-std::sqrt(70875.0) * x * y);

        pc_71[k] = e_0 * (std::sqrt(7.3828125) * x * x * x * x * x * x * y * z - std::sqrt(738.28125) * x * x * x * x * y * y * y * z + std::sqrt(184.5703125) * x * x * y * y * y * y * y * z) + e_1 * (-std::sqrt(1661.1328125) * x * x * x * x * y * z - std::sqrt(738.28125) * x * x * y * y * y * z + std::sqrt(184.5703125) * y * y * y * y * y * z) + e_2 * (-std::sqrt(26578.125) * x * x * y * z + std::sqrt(2953.125) * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_72[k] = e_0 * (-std::sqrt(0.1845703125) * x * x * x * x * x * x * x * y + std::sqrt(14.9501953125) * x * x * x * x * x * y * y * y + std::sqrt(2.953125) * x * x * x * x * x * y * z * z + std::sqrt(4.6142578125) * x * x * x * y * y * y * y * y - std::sqrt(295.3125) * x * x * x * y * y * y * z * z - std::sqrt(4.6142578125) * x * y * y * y * y * y * y * y + std::sqrt(73.828125) * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(18.45703125) * x * x * x * x * x * y + std::sqrt(1845.703125) * x * x * x * y * y * y - std::sqrt(1181.25) * x * x * x * y * z * z - std::sqrt(904.39453125) * x * y * y * y * y * y + std::sqrt(1181.25) * x * y * y * y * z * z) + e_2 * (std::sqrt(4725.0) * x * x * x * y - std::sqrt(4725.0) * x * y * y * y);

        pc_73[k] = e_0 * (-std::sqrt(1.107421875) * x * x * x * x * x * x * x * z + std::sqrt(89.701171875) * x * x * x * x * x * y * y * z + std::sqrt(0.4921875) * x * x * x * x * x * z * z * z + std::sqrt(27.685546875) * x * x * x * y * y * y * y * z - std::sqrt(49.21875) * x * x * x * y * y * z * z * z - std::sqrt(27.685546875) * x * y * y * y * y * y * y * z + std::sqrt(12.3046875) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(110.7421875) * x * x * x * x * x * z + std::sqrt(11074.21875) * x * x * x * y * y * z - std::sqrt(2768.5546875) * x * y * y * y * y * z);

        pc_74[k] = e_0 * (-std::sqrt(0.1845703125) * x * x * x * x * x * x * x * x + std::sqrt(14.9501953125) * x * x * x * x * x * x * y * y + std::sqrt(2.953125) * x * x * x * x * x * x * z * z + std::sqrt(4.6142578125) * x * x * x * x * y * y * y * y - std::sqrt(295.3125) * x * x * x * x * y * y * z * z - std::sqrt(4.6142578125) * x * x * y * y * y * y * y * y + std::sqrt(73.828125) * x * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(41.5283203125) * x * x * x * x * x * x + std::sqrt(2883.9111328125) * x * x * x * x * y * y + std::sqrt(73.828125) * x * x * x * x * z * z - std::sqrt(115.3564453125) * x * x * y * y * y * y - std::sqrt(2657.8125) * x * x * y * y * z * z - std::sqrt(4.6142578125) * y * y * y * y * y * y + std::sqrt(73.828125) * y * y * y * y * z * z) + e_2 * (-std::sqrt(295.3125) * x * x * x * x + std::sqrt(10631.25) * x * x * y * y - std::sqrt(295.3125) * y * y * y * y);

        pc_75[k] = e_0 * (std::sqrt(1.845703125) * x * x * x * x * x * x * x * z - std::sqrt(223.330078125) * x * x * x * x * x * y * y * z + std::sqrt(415.283203125) * x * x * x * y * y * y * y * z - std::sqrt(46.142578125) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(184.5703125) * x * x * x * x * x * z - std::sqrt(738.28125) * x * x * x * y * y * z - std::sqrt(1661.1328125) * x * y * y * y * y * z) + e_2 * (std::sqrt(2953.125) * x * x * x * z - std::sqrt(26578.125) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_76[k] = e_0 * (std::sqrt(0.3076171875) * x * x * x * x * x * x * x * x - std::sqrt(51.9873046875) * x * x * x * x * x * x * y * y + std::sqrt(376.8310546875) * x * x * x * x * y * y * y * y - std::sqrt(69.2138671875) * x * x * y * y * y * y * y * y) + e_1 * (std::sqrt(69.2138671875) * x * x * x * x * x * x + std::sqrt(69.2138671875) * x * x * x * x * y * y - std::sqrt(69.2138671875) * x * x * y * y * y * y - std::sqrt(69.2138671875) * y * y * y * y * y * y) + e_2 * (std::sqrt(4429.6875) * x * x * x * x - std::sqrt(4429.6875) * y * y * y * y) + e_3 * (std::sqrt(17718.75) * x * x - std::sqrt(17718.75) * y * y);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 77; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
