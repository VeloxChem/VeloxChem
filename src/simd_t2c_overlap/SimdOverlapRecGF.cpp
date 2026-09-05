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



#include "SimdOverlapRecGF.hpp"

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
compute_gf_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGF.compute_gf_overlap: Basis functions must be of angular momenta four and three"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGF.compute_gf_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 63 * nvalues, 0.0);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fbe * fh * fh;

        const auto f_3 = fbase * fal * fh * fh * fh;

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

    // NOTE: the components are formed in 16 loops, as the vectorizer runs out
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

        pc_0[k] = e_0 * (std::sqrt(49.21875) * x * x * x * x * x * y * y - std::sqrt(87.5) * x * x * x * y * y * y * y + std::sqrt(5.46875) * x * y * y * y * y * y * y) + e_1 * (std::sqrt(49.21875) * x * x * x * x * x + std::sqrt(196.875) * x * x * x * y * y + std::sqrt(49.21875) * x * y * y * y * y) + e_2 * (std::sqrt(1771.875) * x * x * x + std::sqrt(1771.875) * x * y * y) + e_3 * (std::sqrt(3150.0) * x);

        pc_1[k] = e_0 * (std::sqrt(131.25) * x * x * x * x * y * y * z - std::sqrt(131.25) * x * x * y * y * y * y * z) + e_1 * (std::sqrt(131.25) * x * x * x * x * z - std::sqrt(131.25) * y * y * y * y * z) + e_2 * (std::sqrt(1181.25) * x * x * z - std::sqrt(1181.25) * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(3.28125) * x * x * x * x * x * y * y + std::sqrt(52.5) * x * x * x * y * y * z * z + std::sqrt(3.28125) * x * y * y * y * y * y * y - std::sqrt(52.5) * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(3.28125) * x * x * x * x * x - std::sqrt(118.125) * x * x * x * y * y + std::sqrt(52.5) * x * x * x * z * z + std::sqrt(397.03125) * x * y * y * y * y - std::sqrt(472.5) * x * y * y * z * z) + e_2 * (-std::sqrt(118.125) * x * x * x + std::sqrt(1063.125) * x * y * y);

        pc_3[k] = e_0 * (-std::sqrt(19.6875) * x * x * x * x * x * y * z + std::sqrt(8.75) * x * x * x * y * z * z * z + std::sqrt(19.6875) * x * y * y * y * y * y * z - std::sqrt(8.75) * x * y * y * y * z * z * z) + e_1 * (-std::sqrt(1260.0) * x * x * x * y * z + std::sqrt(1260.0) * x * y * y * y * z);
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

        pc_4[k] = e_0 * (-std::sqrt(3.28125) * x * x * x * x * x * x * y + std::sqrt(52.5) * x * x * x * x * y * z * z + std::sqrt(3.28125) * x * x * y * y * y * y * y - std::sqrt(52.5) * x * x * y * y * y * z * z) + e_1 * (-std::sqrt(397.03125) * x * x * x * x * y + std::sqrt(118.125) * x * x * y * y * y + std::sqrt(472.5) * x * x * y * z * z + std::sqrt(3.28125) * y * y * y * y * y - std::sqrt(52.5) * y * y * y * z * z) + e_2 * (-std::sqrt(1063.125) * x * x * y + std::sqrt(118.125) * y * y * y);

        pc_5[k] = e_0 * (std::sqrt(32.8125) * x * x * x * x * x * y * z - std::sqrt(131.25) * x * x * x * y * y * y * z + std::sqrt(32.8125) * x * y * y * y * y * y * z) + e_1 * (std::sqrt(525.0) * x * x * x * y * z + std::sqrt(525.0) * x * y * y * y * z) + e_2 * (std::sqrt(4725.0) * x * y * z);

        pc_6[k] = e_0 * (std::sqrt(5.46875) * x * x * x * x * x * x * y - std::sqrt(87.5) * x * x * x * x * y * y * y + std::sqrt(49.21875) * x * x * y * y * y * y * y) + e_1 * (std::sqrt(49.21875) * x * x * x * x * y + std::sqrt(196.875) * x * x * y * y * y + std::sqrt(49.21875) * y * y * y * y * y) + e_2 * (std::sqrt(1771.875) * x * x * y + std::sqrt(1771.875) * y * y * y) + e_3 * (std::sqrt(3150.0) * y);

        pc_7[k] = e_0 * (std::sqrt(221.484375) * x * x * x * x * y * y * z - std::sqrt(98.4375) * x * x * y * y * y * y * z + std::sqrt(2.734375) * y * y * y * y * y * y * z) + e_1 * (std::sqrt(221.484375) * x * x * x * x * z + std::sqrt(885.9375) * x * x * y * y * z + std::sqrt(221.484375) * y * y * y * y * z) + e_2 * (std::sqrt(3543.75) * x * x * z + std::sqrt(3543.75) * y * y * z) + e_3 * (std::sqrt(1575.0) * z);
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

        pc_8[k] = e_0 * (std::sqrt(590.625) * x * x * x * y * y * z * z - std::sqrt(65.625) * x * y * y * y * y * z * z) + e_1 * (std::sqrt(590.625) * x * x * x * y * y + std::sqrt(590.625) * x * x * x * z * z - std::sqrt(65.625) * x * y * y * y * y + std::sqrt(590.625) * x * y * y * z * z) + e_2 * (std::sqrt(590.625) * x * x * x + std::sqrt(590.625) * x * y * y + std::sqrt(2362.5) * x * z * z) + e_3 * (std::sqrt(2362.5) * x);

        pc_9[k] = e_0 * (-std::sqrt(14.765625) * x * x * x * x * y * y * z - std::sqrt(6.5625) * x * x * y * y * y * y * z + std::sqrt(236.25) * x * x * y * y * z * z * z + std::sqrt(1.640625) * y * y * y * y * y * y * z - std::sqrt(26.25) * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(14.765625) * x * x * x * x * z + std::sqrt(59.0625) * x * x * y * y * z + std::sqrt(236.25) * x * x * z * z * z + std::sqrt(1.640625) * y * y * y * y * z - std::sqrt(236.25) * y * y * z * z * z) + e_2 * (std::sqrt(236.25) * x * x * z - std::sqrt(236.25) * y * y * z);

        pc_10[k] = e_0 * (-std::sqrt(88.59375) * x * x * x * x * y * z * z - std::sqrt(39.375) * x * x * y * y * y * z * z + std::sqrt(39.375) * x * x * y * z * z * z * z + std::sqrt(9.84375) * y * y * y * y * y * z * z - std::sqrt(4.375) * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(88.59375) * x * x * x * x * y - std::sqrt(39.375) * x * x * y * y * y - std::sqrt(1417.5) * x * x * y * z * z + std::sqrt(9.84375) * y * y * y * y * y + std::sqrt(157.5) * y * y * y * z * z) + e_2 * (-std::sqrt(3189.375) * x * x * y + std::sqrt(354.375) * y * y * y);

        pc_11[k] = e_0 * (-std::sqrt(14.765625) * x * x * x * x * x * y * z - std::sqrt(6.5625) * x * x * x * y * y * y * z + std::sqrt(236.25) * x * x * x * y * z * z * z + std::sqrt(1.640625) * x * y * y * y * y * y * z - std::sqrt(26.25) * x * y * y * y * z * z * z) + e_1 * (-std::sqrt(105.0) * x * y * y * y * z + std::sqrt(945.0) * x * y * z * z * z) + e_2 * (std::sqrt(945.0) * x * y * z);
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

        pc_12[k] = e_0 * (std::sqrt(147.65625) * x * x * x * x * y * z * z - std::sqrt(262.5) * x * x * y * y * y * z * z + std::sqrt(16.40625) * y * y * y * y * y * z * z) + e_1 * (std::sqrt(147.65625) * x * x * x * x * y - std::sqrt(262.5) * x * x * y * y * y + std::sqrt(590.625) * x * x * y * z * z + std::sqrt(16.40625) * y * y * y * y * y + std::sqrt(590.625) * y * y * y * z * z) + e_2 * (std::sqrt(590.625) * x * x * y + std::sqrt(590.625) * y * y * y + std::sqrt(2362.5) * y * z * z) + e_3 * (std::sqrt(2362.5) * y);

        pc_13[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * y * z - std::sqrt(273.4375) * x * x * x * y * y * y * z + std::sqrt(24.609375) * x * y * y * y * y * y * z);

        pc_14[k] = e_0 * (-std::sqrt(7.03125) * x * x * x * x * x * y * y - std::sqrt(3.125) * x * x * x * y * y * y * y + std::sqrt(253.125) * x * x * x * y * y * z * z + std::sqrt(0.78125) * x * y * y * y * y * y * y - std::sqrt(28.125) * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(7.03125) * x * x * x * x * x - std::sqrt(450.0) * x * x * x * y * y + std::sqrt(253.125) * x * x * x * z * z + std::sqrt(7.03125) * x * y * y * y * y + std::sqrt(253.125) * x * y * y * z * z) + e_2 * (-std::sqrt(253.125) * x * x * x - std::sqrt(253.125) * x * y * y + std::sqrt(1012.5) * x * z * z) + e_3 * (-std::sqrt(112.5) * x);

        pc_15[k] = e_0 * (-std::sqrt(18.75) * x * x * x * x * y * y * z - std::sqrt(18.75) * x * x * y * y * y * y * z + std::sqrt(675.0) * x * x * y * y * z * z * z) + e_1 * (-std::sqrt(18.75) * x * x * x * x * z + std::sqrt(675.0) * x * x * y * y * z + std::sqrt(675.0) * x * x * z * z * z - std::sqrt(18.75) * y * y * y * y * z + std::sqrt(675.0) * y * y * z * z * z) + e_2 * (std::sqrt(1518.75) * x * x * z + std::sqrt(1518.75) * y * y * z + std::sqrt(675.0) * z * z * z) + e_3 * (std::sqrt(2700.0) * z);
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

        pc_16[k] = e_0 * (std::sqrt(0.46875) * x * x * x * x * x * y * y + std::sqrt(1.875) * x * x * x * y * y * y * y - std::sqrt(46.875) * x * x * x * y * y * z * z + std::sqrt(0.46875) * x * y * y * y * y * y * y - std::sqrt(46.875) * x * y * y * y * y * z * z + std::sqrt(270.0) * x * y * y * z * z * z * z) + e_1 * (std::sqrt(0.46875) * x * x * x * x * x + std::sqrt(67.5) * x * x * x * y * y - std::sqrt(46.875) * x * x * x * z * z + std::sqrt(56.71875) * x * y * y * y * y + std::sqrt(1366.875) * x * y * y * z * z + std::sqrt(270.0) * x * z * z * z * z) + e_2 * (std::sqrt(16.875) * x * x * x + std::sqrt(2851.875) * x * y * y + std::sqrt(3307.5) * x * z * z) + e_3 * (std::sqrt(1687.5) * x);

        pc_17[k] = e_0 * (std::sqrt(2.8125) * x * x * x * x * x * y * z + std::sqrt(11.25) * x * x * x * y * y * y * z - std::sqrt(125.0) * x * x * x * y * z * z * z + std::sqrt(2.8125) * x * y * y * y * y * y * z - std::sqrt(125.0) * x * y * y * y * z * z * z + std::sqrt(45.0) * x * y * z * z * z * z * z) + e_1 * (-std::sqrt(45.0) * x * x * x * y * z - std::sqrt(45.0) * x * y * y * y * z) + e_2 * (-std::sqrt(405.0) * x * y * z);

        pc_18[k] = e_0 * (std::sqrt(0.46875) * x * x * x * x * x * x * y + std::sqrt(1.875) * x * x * x * x * y * y * y - std::sqrt(46.875) * x * x * x * x * y * z * z + std::sqrt(0.46875) * x * x * y * y * y * y * y - std::sqrt(46.875) * x * x * y * y * y * z * z + std::sqrt(270.0) * x * x * y * z * z * z * z) + e_1 * (std::sqrt(56.71875) * x * x * x * x * y + std::sqrt(67.5) * x * x * y * y * y + std::sqrt(1366.875) * x * x * y * z * z + std::sqrt(0.46875) * y * y * y * y * y - std::sqrt(46.875) * y * y * y * z * z + std::sqrt(270.0) * y * z * z * z * z) + e_2 * (std::sqrt(2851.875) * x * x * y + std::sqrt(16.875) * y * y * y + std::sqrt(3307.5) * y * z * z) + e_3 * (std::sqrt(1687.5) * y);

        pc_19[k] = e_0 * (-std::sqrt(4.6875) * x * x * x * x * x * y * z + std::sqrt(168.75) * x * x * x * y * z * z * z + std::sqrt(4.6875) * x * y * y * y * y * y * z - std::sqrt(168.75) * x * y * y * y * z * z * z) + e_1 * (std::sqrt(300.0) * x * x * x * y * z - std::sqrt(300.0) * x * y * y * y * z);
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

        pc_20[k] = e_0 * (-std::sqrt(0.78125) * x * x * x * x * x * x * y + std::sqrt(3.125) * x * x * x * x * y * y * y + std::sqrt(28.125) * x * x * x * x * y * z * z + std::sqrt(7.03125) * x * x * y * y * y * y * y - std::sqrt(253.125) * x * x * y * y * y * z * z) + e_1 * (-std::sqrt(7.03125) * x * x * x * x * y + std::sqrt(450.0) * x * x * y * y * y - std::sqrt(253.125) * x * x * y * z * z + std::sqrt(7.03125) * y * y * y * y * y - std::sqrt(253.125) * y * y * y * z * z) + e_2 * (std::sqrt(253.125) * x * x * y + std::sqrt(253.125) * y * y * y - std::sqrt(1012.5) * y * z * z) + e_3 * (std::sqrt(112.5) * y);

        pc_21[k] = e_0 * (-5.625 * x * x * x * x * y * y * z - 3.75 * x * x * y * y * y * y * z + 7.5 * x * x * y * y * z * z * z + 1.875 * y * y * y * y * y * y * z - 2.5 * y * y * y * y * z * z * z) + e_1 * (-5.625 * x * x * x * x * z - 33.75 * x * x * y * y * z + 7.5 * x * x * z * z * z + 16.875 * y * y * y * y * z - 7.5 * y * y * z * z * z) + e_2 * (-22.5 * x * x * z + 22.5 * y * y * z);

        pc_22[k] = e_0 * (-std::sqrt(84.375) * x * x * x * y * y * z * z - std::sqrt(84.375) * x * y * y * y * y * z * z + std::sqrt(150.0) * x * y * y * z * z * z * z) + e_1 * (-std::sqrt(84.375) * x * x * x * y * y - std::sqrt(84.375) * x * x * x * z * z - std::sqrt(84.375) * x * y * y * y * y - std::sqrt(84.375) * x * y * y * z * z + std::sqrt(150.0) * x * z * z * z * z) + e_2 * (-std::sqrt(84.375) * x * x * x - std::sqrt(2109.375) * x * y * y + std::sqrt(337.5) * x * z * z) + e_3 * (-std::sqrt(337.5) * x);

        pc_23[k] = e_0 * (std::sqrt(2.109375) * x * x * x * x * y * y * z + std::sqrt(8.4375) * x * x * y * y * y * y * z - std::sqrt(60.0) * x * x * y * y * z * z * z + std::sqrt(2.109375) * y * y * y * y * y * y * z - std::sqrt(60.0) * y * y * y * y * z * z * z + std::sqrt(60.0) * y * y * z * z * z * z * z) + e_1 * (std::sqrt(2.109375) * x * x * x * x * z + std::sqrt(8.4375) * x * x * y * y * z - std::sqrt(60.0) * x * x * z * z * z + std::sqrt(2.109375) * y * y * y * y * z + std::sqrt(540.0) * y * y * z * z * z + std::sqrt(60.0) * z * z * z * z * z) + e_2 * (-std::sqrt(33.75) * x * x * z + std::sqrt(1653.75) * y * y * z + std::sqrt(2160.0) * z * z * z) + e_3 * (std::sqrt(3375.0) * z);
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

        pc_24[k] = e_0 * (std::sqrt(12.65625) * x * x * x * x * y * z * z + std::sqrt(50.625) * x * x * y * y * y * z * z - std::sqrt(50.625) * x * x * y * z * z * z * z + std::sqrt(12.65625) * y * y * y * y * y * z * z - std::sqrt(50.625) * y * y * y * z * z * z * z + std::sqrt(10.0) * y * z * z * z * z * z * z) + e_1 * (std::sqrt(12.65625) * x * x * x * x * y + std::sqrt(50.625) * x * x * y * y * y + std::sqrt(12.65625) * y * y * y * y * y + std::sqrt(360.0) * y * z * z * z * z) + e_2 * (std::sqrt(455.625) * x * x * y + std::sqrt(455.625) * y * y * y + std::sqrt(3240.0) * y * z * z) + e_3 * (std::sqrt(2250.0) * y);

        pc_25[k] = e_0 * (std::sqrt(2.109375) * x * x * x * x * x * y * z + std::sqrt(8.4375) * x * x * x * y * y * y * z - std::sqrt(60.0) * x * x * x * y * z * z * z + std::sqrt(2.109375) * x * y * y * y * y * y * z - std::sqrt(60.0) * x * y * y * y * z * z * z + std::sqrt(60.0) * x * y * z * z * z * z * z) + e_1 * (std::sqrt(960.0) * x * y * z * z * z) + e_2 * (std::sqrt(2160.0) * x * y * z);

        pc_26[k] = e_0 * (-std::sqrt(21.09375) * x * x * x * x * y * z * z + std::sqrt(37.5) * x * x * y * z * z * z * z + std::sqrt(21.09375) * y * y * y * y * y * z * z - std::sqrt(37.5) * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(21.09375) * x * x * x * x * y + std::sqrt(84.375) * x * x * y * z * z + std::sqrt(21.09375) * y * y * y * y * y + std::sqrt(84.375) * y * y * y * z * z - std::sqrt(150.0) * y * z * z * z * z) + e_2 * (-std::sqrt(84.375) * x * x * y + std::sqrt(759.375) * y * y * y - std::sqrt(337.5) * y * z * z) + e_3 * (std::sqrt(337.5) * y);

        pc_27[k] = e_0 * (-1.875 * x * x * x * x * x * y * z + 3.75 * x * x * x * y * y * y * z + 2.5 * x * x * x * y * z * z * z + 5.625 * x * y * y * y * y * y * z - 7.5 * x * y * y * y * z * z * z) + e_1 * (45.0 * x * y * y * y * z - 15.0 * x * y * z * z * z) + e_2 * (45.0 * x * y * z);
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

        pc_28[k] = e_0 * (std::sqrt(0.791015625) * x * x * x * x * x * x * y + std::sqrt(2.197265625) * x * x * x * x * y * y * y - std::sqrt(50.625) * x * x * x * x * y * z * z + std::sqrt(0.087890625) * x * x * y * y * y * y * y - std::sqrt(22.5) * x * x * y * y * y * z * z + std::sqrt(5.625) * x * x * y * z * z * z * z - std::sqrt(0.087890625) * y * y * y * y * y * y * y + std::sqrt(5.625) * y * y * y * y * y * z * z - std::sqrt(0.625) * y * y * y * z * z * z * z) + e_1 * (std::sqrt(113.90625) * x * x * x * x * y + std::sqrt(50.625) * x * x * y * y * y - std::sqrt(1822.5) * x * x * y * z * z - std::sqrt(12.65625) * y * y * y * y * y + std::sqrt(202.5) * y * y * y * z * z) + e_2 * (std::sqrt(455.625) * x * x * y - std::sqrt(50.625) * y * y * y);

        pc_29[k] = e_0 * (std::sqrt(2.109375) * x * x * x * x * x * y * z + std::sqrt(8.4375) * x * x * x * y * y * y * z - std::sqrt(135.0) * x * x * x * y * z * z * z + std::sqrt(2.109375) * x * y * y * y * y * y * z - std::sqrt(135.0) * x * y * y * y * z * z * z + std::sqrt(15.0) * x * y * z * z * z * z * z) + e_1 * (-std::sqrt(135.0) * x * x * x * y * z - std::sqrt(135.0) * x * y * y * y * z - std::sqrt(960.0) * x * y * z * z * z) + e_2 * (-std::sqrt(6615.0) * x * y * z);

        pc_30[k] = e_0 * (-std::sqrt(0.052734375) * x * x * x * x * x * x * y - std::sqrt(0.474609375) * x * x * x * x * y * y * y + std::sqrt(7.59375) * x * x * x * x * y * z * z - std::sqrt(0.474609375) * x * x * y * y * y * y * y + std::sqrt(30.375) * x * x * y * y * y * z * z - std::sqrt(63.375) * x * x * y * z * z * z * z - std::sqrt(0.052734375) * y * y * y * y * y * y * y + std::sqrt(7.59375) * y * y * y * y * y * z * z - std::sqrt(63.375) * y * y * y * z * z * z * z + std::sqrt(6.0) * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(7.59375) * x * x * x * x * y - std::sqrt(30.375) * x * x * y * y * y - std::sqrt(216.0) * x * x * y * z * z - std::sqrt(7.59375) * y * y * y * y * y - std::sqrt(216.0) * y * y * y * z * z + std::sqrt(24.0) * y * z * z * z * z) + e_2 * (-std::sqrt(570.375) * x * x * y - std::sqrt(570.375) * y * y * y - std::sqrt(216.0) * y * z * z) + e_3 * (-std::sqrt(1350.0) * y);

        pc_31[k] = e_0 * (-0.5625 * x * x * x * x * x * x * z - 1.6875 * x * x * x * x * y * y * z + 4.875 * x * x * x * x * z * z * z - 1.6875 * x * x * y * y * y * y * z + 9.75 * x * x * y * y * z * z * z - 4.5 * x * x * z * z * z * z * z - 0.5625 * y * y * y * y * y * y * z + 4.875 * y * y * y * y * z * z * z - 4.5 * y * y * z * z * z * z * z + z * z * z * z * z * z * z) + e_1 * (4.5 * x * x * x * x * z + 9.0 * x * x * y * y * z - 6.0 * x * x * z * z * z + 4.5 * y * y * y * y * z - 6.0 * y * y * z * z * z + 12.0 * z * z * z * z * z) + e_2 * (9.0 * x * x * z + 9.0 * y * y * z + 54.0 * z * z * z) + e_3 * (60.0 * z);
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

        pc_32[k] = e_0 * (-std::sqrt(0.052734375) * x * x * x * x * x * x * x - std::sqrt(0.474609375) * x * x * x * x * x * y * y + std::sqrt(7.59375) * x * x * x * x * x * z * z - std::sqrt(0.474609375) * x * x * x * y * y * y * y + std::sqrt(30.375) * x * x * x * y * y * z * z - std::sqrt(63.375) * x * x * x * z * z * z * z - std::sqrt(0.052734375) * x * y * y * y * y * y * y + std::sqrt(7.59375) * x * y * y * y * y * z * z - std::sqrt(63.375) * x * y * y * z * z * z * z + std::sqrt(6.0) * x * z * z * z * z * z * z) + e_1 * (-std::sqrt(7.59375) * x * x * x * x * x - std::sqrt(30.375) * x * x * x * y * y - std::sqrt(216.0) * x * x * x * z * z - std::sqrt(7.59375) * x * y * y * y * y - std::sqrt(216.0) * x * y * y * z * z + std::sqrt(24.0) * x * z * z * z * z) + e_2 * (-std::sqrt(570.375) * x * x * x - std::sqrt(570.375) * x * y * y - std::sqrt(216.0) * x * z * z) + e_3 * (-std::sqrt(1350.0) * x);

        pc_33[k] = e_0 * (std::sqrt(0.52734375) * x * x * x * x * x * x * z + std::sqrt(0.52734375) * x * x * x * x * y * y * z - std::sqrt(33.75) * x * x * x * x * z * z * z - std::sqrt(0.52734375) * x * x * y * y * y * y * z + std::sqrt(3.75) * x * x * z * z * z * z * z - std::sqrt(0.52734375) * y * y * y * y * y * y * z + std::sqrt(33.75) * y * y * y * y * z * z * z - std::sqrt(3.75) * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(33.75) * x * x * x * x * z - std::sqrt(240.0) * x * x * z * z * z + std::sqrt(33.75) * y * y * y * y * z + std::sqrt(240.0) * y * y * z * z * z) + e_2 * (-std::sqrt(1653.75) * x * x * z + std::sqrt(1653.75) * y * y * z);

        pc_34[k] = e_0 * (std::sqrt(0.087890625) * x * x * x * x * x * x * x - std::sqrt(0.087890625) * x * x * x * x * x * y * y - std::sqrt(5.625) * x * x * x * x * x * z * z - std::sqrt(2.197265625) * x * x * x * y * y * y * y + std::sqrt(22.5) * x * x * x * y * y * z * z + std::sqrt(0.625) * x * x * x * z * z * z * z - std::sqrt(0.791015625) * x * y * y * y * y * y * y + std::sqrt(50.625) * x * y * y * y * y * z * z - std::sqrt(5.625) * x * y * y * z * z * z * z) + e_1 * (std::sqrt(12.65625) * x * x * x * x * x - std::sqrt(50.625) * x * x * x * y * y - std::sqrt(202.5) * x * x * x * z * z - std::sqrt(113.90625) * x * y * y * y * y + std::sqrt(1822.5) * x * y * y * z * z) + e_2 * (std::sqrt(50.625) * x * x * x - std::sqrt(455.625) * x * y * y);

        pc_35[k] = e_0 * (-5.625 * x * x * x * x * x * y * z - 3.75 * x * x * x * y * y * y * z + 7.5 * x * x * x * y * z * z * z + 1.875 * x * y * y * y * y * y * z - 2.5 * x * y * y * y * z * z * z) + e_1 * (-45.0 * x * x * x * y * z + 15.0 * x * y * z * z * z) + e_2 * (-45.0 * x * y * z);
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

        pc_36[k] = e_0 * (-std::sqrt(84.375) * x * x * x * x * y * z * z - std::sqrt(84.375) * x * x * y * y * y * z * z + std::sqrt(150.0) * x * x * y * z * z * z * z) + e_1 * (-std::sqrt(84.375) * x * x * x * x * y - std::sqrt(84.375) * x * x * y * y * y - std::sqrt(84.375) * x * x * y * z * z - std::sqrt(84.375) * y * y * y * z * z + std::sqrt(150.0) * y * z * z * z * z) + e_2 * (-std::sqrt(2109.375) * x * x * y - std::sqrt(84.375) * y * y * y + std::sqrt(337.5) * y * z * z) + e_3 * (-std::sqrt(337.5) * y);

        pc_37[k] = e_0 * (std::sqrt(2.109375) * x * x * x * x * x * y * z + std::sqrt(8.4375) * x * x * x * y * y * y * z - std::sqrt(60.0) * x * x * x * y * z * z * z + std::sqrt(2.109375) * x * y * y * y * y * y * z - std::sqrt(60.0) * x * y * y * y * z * z * z + std::sqrt(60.0) * x * y * z * z * z * z * z) + e_1 * (std::sqrt(960.0) * x * y * z * z * z) + e_2 * (std::sqrt(2160.0) * x * y * z);

        pc_38[k] = e_0 * (std::sqrt(12.65625) * x * x * x * x * x * z * z + std::sqrt(50.625) * x * x * x * y * y * z * z - std::sqrt(50.625) * x * x * x * z * z * z * z + std::sqrt(12.65625) * x * y * y * y * y * z * z - std::sqrt(50.625) * x * y * y * z * z * z * z + std::sqrt(10.0) * x * z * z * z * z * z * z) + e_1 * (std::sqrt(12.65625) * x * x * x * x * x + std::sqrt(50.625) * x * x * x * y * y + std::sqrt(12.65625) * x * y * y * y * y + std::sqrt(360.0) * x * z * z * z * z) + e_2 * (std::sqrt(455.625) * x * x * x + std::sqrt(455.625) * x * y * y + std::sqrt(3240.0) * x * z * z) + e_3 * (std::sqrt(2250.0) * x);

        pc_39[k] = e_0 * (std::sqrt(2.109375) * x * x * x * x * x * x * z + std::sqrt(8.4375) * x * x * x * x * y * y * z - std::sqrt(60.0) * x * x * x * x * z * z * z + std::sqrt(2.109375) * x * x * y * y * y * y * z - std::sqrt(60.0) * x * x * y * y * z * z * z + std::sqrt(60.0) * x * x * z * z * z * z * z) + e_1 * (std::sqrt(2.109375) * x * x * x * x * z + std::sqrt(8.4375) * x * x * y * y * z + std::sqrt(540.0) * x * x * z * z * z + std::sqrt(2.109375) * y * y * y * y * z - std::sqrt(60.0) * y * y * z * z * z + std::sqrt(60.0) * z * z * z * z * z) + e_2 * (std::sqrt(1653.75) * x * x * z - std::sqrt(33.75) * y * y * z + std::sqrt(2160.0) * z * z * z) + e_3 * (std::sqrt(3375.0) * z);
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

        pc_40[k] = e_0 * (-std::sqrt(21.09375) * x * x * x * x * x * z * z + std::sqrt(37.5) * x * x * x * z * z * z * z + std::sqrt(21.09375) * x * y * y * y * y * z * z - std::sqrt(37.5) * x * y * y * z * z * z * z) + e_1 * (-std::sqrt(21.09375) * x * x * x * x * x - std::sqrt(84.375) * x * x * x * z * z + std::sqrt(21.09375) * x * y * y * y * y - std::sqrt(84.375) * x * y * y * z * z + std::sqrt(150.0) * x * z * z * z * z) + e_2 * (-std::sqrt(759.375) * x * x * x + std::sqrt(84.375) * x * y * y + std::sqrt(337.5) * x * z * z) + e_3 * (-std::sqrt(337.5) * x);

        pc_41[k] = e_0 * (-1.875 * x * x * x * x * x * x * z + 3.75 * x * x * x * x * y * y * z + 2.5 * x * x * x * x * z * z * z + 5.625 * x * x * y * y * y * y * z - 7.5 * x * x * y * y * z * z * z) + e_1 * (-16.875 * x * x * x * x * z + 33.75 * x * x * y * y * z + 7.5 * x * x * z * z * z + 5.625 * y * y * y * y * z - 7.5 * y * y * z * z * z) + e_2 * (-22.5 * x * x * z + 22.5 * y * y * z);

        pc_42[k] = e_0 * (-std::sqrt(1.7578125) * x * x * x * x * x * x * y + std::sqrt(0.1953125) * x * x * x * x * y * y * y + std::sqrt(63.28125) * x * x * x * x * y * z * z + std::sqrt(1.7578125) * x * x * y * y * y * y * y - std::sqrt(112.5) * x * x * y * y * y * z * z - std::sqrt(0.1953125) * y * y * y * y * y * y * y + std::sqrt(7.03125) * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(112.5) * x * x * x * x * y + std::sqrt(28.125) * x * x * y * y * y + std::sqrt(253.125) * x * x * y * z * z - std::sqrt(28.125) * y * y * y * y * y + std::sqrt(253.125) * y * y * y * z * z) + e_2 * (-std::sqrt(253.125) * x * x * y - std::sqrt(253.125) * y * y * y + std::sqrt(1012.5) * y * z * z) + e_3 * (-std::sqrt(112.5) * y);

        pc_43[k] = e_0 * (-std::sqrt(4.6875) * x * x * x * x * x * y * z + std::sqrt(168.75) * x * x * x * y * z * z * z + std::sqrt(4.6875) * x * y * y * y * y * y * z - std::sqrt(168.75) * x * y * y * y * z * z * z) + e_1 * (std::sqrt(300.0) * x * x * x * y * z - std::sqrt(300.0) * x * y * y * y * z);
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

        pc_44[k] = e_0 * (std::sqrt(0.1171875) * x * x * x * x * x * x * y + std::sqrt(0.1171875) * x * x * x * x * y * y * y - std::sqrt(11.71875) * x * x * x * x * y * z * z - std::sqrt(0.1171875) * x * x * y * y * y * y * y + std::sqrt(67.5) * x * x * y * z * z * z * z - std::sqrt(0.1171875) * y * y * y * y * y * y * y + std::sqrt(11.71875) * y * y * y * y * y * z * z - std::sqrt(67.5) * y * y * y * z * z * z * z) + e_1 * (std::sqrt(7.5) * x * x * x * x * y - std::sqrt(1.875) * x * x * y * y * y + std::sqrt(826.875) * x * x * y * z * z - std::sqrt(16.875) * y * y * y * y * y - std::sqrt(226.875) * y * y * y * z * z - std::sqrt(270.0) * y * z * z * z * z) + e_2 * (std::sqrt(421.875) * x * x * y - std::sqrt(826.875) * y * y * y - std::sqrt(3307.5) * y * z * z) + e_3 * (-std::sqrt(1687.5) * y);

        pc_45[k] = e_0 * (std::sqrt(0.703125) * x * x * x * x * x * x * z + std::sqrt(0.703125) * x * x * x * x * y * y * z - std::sqrt(31.25) * x * x * x * x * z * z * z - std::sqrt(0.703125) * x * x * y * y * y * y * z + std::sqrt(11.25) * x * x * z * z * z * z * z - std::sqrt(0.703125) * y * y * y * y * y * y * z + std::sqrt(31.25) * y * y * y * y * z * z * z - std::sqrt(11.25) * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(11.25) * x * x * x * x * z + std::sqrt(11.25) * y * y * y * y * z) + e_2 * (-std::sqrt(101.25) * x * x * z + std::sqrt(101.25) * y * y * z);

        pc_46[k] = e_0 * (std::sqrt(0.1171875) * x * x * x * x * x * x * x + std::sqrt(0.1171875) * x * x * x * x * x * y * y - std::sqrt(11.71875) * x * x * x * x * x * z * z - std::sqrt(0.1171875) * x * x * x * y * y * y * y + std::sqrt(67.5) * x * x * x * z * z * z * z - std::sqrt(0.1171875) * x * y * y * y * y * y * y + std::sqrt(11.71875) * x * y * y * y * y * z * z - std::sqrt(67.5) * x * y * y * z * z * z * z) + e_1 * (std::sqrt(16.875) * x * x * x * x * x + std::sqrt(1.875) * x * x * x * y * y + std::sqrt(226.875) * x * x * x * z * z - std::sqrt(7.5) * x * y * y * y * y - std::sqrt(826.875) * x * y * y * z * z + std::sqrt(270.0) * x * z * z * z * z) + e_2 * (std::sqrt(826.875) * x * x * x - std::sqrt(421.875) * x * y * y + std::sqrt(3307.5) * x * z * z) + e_3 * (std::sqrt(1687.5) * x);

        pc_47[k] = e_0 * (-std::sqrt(1.171875) * x * x * x * x * x * x * z + std::sqrt(1.171875) * x * x * x * x * y * y * z + std::sqrt(42.1875) * x * x * x * x * z * z * z + std::sqrt(1.171875) * x * x * y * y * y * y * z - std::sqrt(168.75) * x * x * y * y * z * z * z - std::sqrt(1.171875) * y * y * y * y * y * y * z + std::sqrt(42.1875) * y * y * y * y * z * z * z) + e_1 * (std::sqrt(18.75) * x * x * x * x * z - std::sqrt(675.0) * x * x * y * y * z + std::sqrt(675.0) * x * x * z * z * z + std::sqrt(18.75) * y * y * y * y * z + std::sqrt(675.0) * y * y * z * z * z) + e_2 * (std::sqrt(1518.75) * x * x * z + std::sqrt(1518.75) * y * y * z + std::sqrt(675.0) * z * z * z) + e_3 * (std::sqrt(2700.0) * z);
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

        pc_48[k] = e_0 * (-std::sqrt(0.1953125) * x * x * x * x * x * x * x + std::sqrt(1.7578125) * x * x * x * x * x * y * y + std::sqrt(7.03125) * x * x * x * x * x * z * z + std::sqrt(0.1953125) * x * x * x * y * y * y * y - std::sqrt(112.5) * x * x * x * y * y * z * z - std::sqrt(1.7578125) * x * y * y * y * y * y * y + std::sqrt(63.28125) * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(28.125) * x * x * x * x * x + std::sqrt(28.125) * x * x * x * y * y + std::sqrt(253.125) * x * x * x * z * z - std::sqrt(112.5) * x * y * y * y * y + std::sqrt(253.125) * x * y * y * z * z) + e_2 * (-std::sqrt(253.125) * x * x * x - std::sqrt(253.125) * x * y * y + std::sqrt(1012.5) * x * z * z) + e_3 * (-std::sqrt(112.5) * x);

        pc_49[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * y * z - std::sqrt(273.4375) * x * x * x * y * y * y * z + std::sqrt(24.609375) * x * y * y * y * y * y * z);

        pc_50[k] = e_0 * (std::sqrt(65.625) * x * x * x * x * y * z * z - std::sqrt(590.625) * x * x * y * y * y * z * z) + e_1 * (std::sqrt(65.625) * x * x * x * x * y - std::sqrt(590.625) * x * x * y * y * y - std::sqrt(590.625) * x * x * y * z * z - std::sqrt(590.625) * y * y * y * z * z) + e_2 * (-std::sqrt(590.625) * x * x * y - std::sqrt(590.625) * y * y * y - std::sqrt(2362.5) * y * z * z) + e_3 * (-std::sqrt(2362.5) * y);

        pc_51[k] = e_0 * (-std::sqrt(1.640625) * x * x * x * x * x * y * z + std::sqrt(6.5625) * x * x * x * y * y * y * z + std::sqrt(26.25) * x * x * x * y * z * z * z + std::sqrt(14.765625) * x * y * y * y * y * y * z - std::sqrt(236.25) * x * y * y * y * z * z * z) + e_1 * (std::sqrt(105.0) * x * x * x * y * z - std::sqrt(945.0) * x * y * z * z * z) + e_2 * (-std::sqrt(945.0) * x * y * z);
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

        pc_52[k] = e_0 * (-std::sqrt(9.84375) * x * x * x * x * x * z * z + std::sqrt(39.375) * x * x * x * y * y * z * z + std::sqrt(4.375) * x * x * x * z * z * z * z + std::sqrt(88.59375) * x * y * y * y * y * z * z - std::sqrt(39.375) * x * y * y * z * z * z * z) + e_1 * (-std::sqrt(9.84375) * x * x * x * x * x + std::sqrt(39.375) * x * x * x * y * y - std::sqrt(157.5) * x * x * x * z * z + std::sqrt(88.59375) * x * y * y * y * y + std::sqrt(1417.5) * x * y * y * z * z) + e_2 * (-std::sqrt(354.375) * x * x * x + std::sqrt(3189.375) * x * y * y);

        pc_53[k] = e_0 * (-std::sqrt(1.640625) * x * x * x * x * x * x * z + std::sqrt(6.5625) * x * x * x * x * y * y * z + std::sqrt(26.25) * x * x * x * x * z * z * z + std::sqrt(14.765625) * x * x * y * y * y * y * z - std::sqrt(236.25) * x * x * y * y * z * z * z) + e_1 * (-std::sqrt(1.640625) * x * x * x * x * z - std::sqrt(59.0625) * x * x * y * y * z + std::sqrt(236.25) * x * x * z * z * z + std::sqrt(14.765625) * y * y * y * y * z - std::sqrt(236.25) * y * y * z * z * z) + e_2 * (std::sqrt(236.25) * x * x * z - std::sqrt(236.25) * y * y * z);

        pc_54[k] = e_0 * (std::sqrt(16.40625) * x * x * x * x * x * z * z - std::sqrt(262.5) * x * x * x * y * y * z * z + std::sqrt(147.65625) * x * y * y * y * y * z * z) + e_1 * (std::sqrt(16.40625) * x * x * x * x * x - std::sqrt(262.5) * x * x * x * y * y + std::sqrt(590.625) * x * x * x * z * z + std::sqrt(147.65625) * x * y * y * y * y + std::sqrt(590.625) * x * y * y * z * z) + e_2 * (std::sqrt(590.625) * x * x * x + std::sqrt(590.625) * x * y * y + std::sqrt(2362.5) * x * z * z) + e_3 * (std::sqrt(2362.5) * x);

        pc_55[k] = e_0 * (std::sqrt(2.734375) * x * x * x * x * x * x * z - std::sqrt(98.4375) * x * x * x * x * y * y * z + std::sqrt(221.484375) * x * x * y * y * y * y * z) + e_1 * (std::sqrt(221.484375) * x * x * x * x * z + std::sqrt(885.9375) * x * x * y * y * z + std::sqrt(221.484375) * y * y * y * y * z) + e_2 * (std::sqrt(3543.75) * x * x * z + std::sqrt(3543.75) * y * y * z) + e_3 * (std::sqrt(1575.0) * z);
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

        pc_56[k] = e_0 * (std::sqrt(3.076171875) * x * x * x * x * x * x * y - std::sqrt(123.388671875) * x * x * x * x * y * y * y + std::sqrt(27.685546875) * x * x * y * y * y * y * y - std::sqrt(0.341796875) * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(49.21875) * x * x * x * x * y - std::sqrt(196.875) * x * x * y * y * y - std::sqrt(49.21875) * y * y * y * y * y) + e_2 * (-std::sqrt(1771.875) * x * x * y - std::sqrt(1771.875) * y * y * y) + e_3 * (-std::sqrt(3150.0) * y);

        pc_57[k] = e_0 * (std::sqrt(8.203125) * x * x * x * x * x * y * z - std::sqrt(295.3125) * x * x * x * y * y * y * z + std::sqrt(8.203125) * x * y * y * y * y * y * z) + e_1 * (-std::sqrt(525.0) * x * x * x * y * z - std::sqrt(525.0) * x * y * y * y * z) + e_2 * (-std::sqrt(4725.0) * x * y * z);

        pc_58[k] = e_0 * (-std::sqrt(0.205078125) * x * x * x * x * x * x * y + std::sqrt(5.126953125) * x * x * x * x * y * y * y + std::sqrt(3.28125) * x * x * x * x * y * z * z + std::sqrt(5.126953125) * x * x * y * y * y * y * y - std::sqrt(118.125) * x * x * y * y * y * z * z - std::sqrt(0.205078125) * y * y * y * y * y * y * y + std::sqrt(3.28125) * y * y * y * y * y * z * z) + e_1 * (std::sqrt(3.28125) * x * x * x * x * y + std::sqrt(643.125) * x * x * y * y * y - std::sqrt(472.5) * x * x * y * z * z - std::sqrt(29.53125) * y * y * y * y * y + std::sqrt(52.5) * y * y * y * z * z) + e_2 * (std::sqrt(1063.125) * x * x * y - std::sqrt(118.125) * y * y * y);

        pc_59[k] = e_0 * (-std::sqrt(1.23046875) * x * x * x * x * x * x * z + std::sqrt(30.76171875) * x * x * x * x * y * y * z + std::sqrt(0.546875) * x * x * x * x * z * z * z + std::sqrt(30.76171875) * x * x * y * y * y * y * z - std::sqrt(19.6875) * x * x * y * y * z * z * z - std::sqrt(1.23046875) * y * y * y * y * y * y * z + std::sqrt(0.546875) * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(78.75) * x * x * x * x * z + std::sqrt(2835.0) * x * x * y * y * z - std::sqrt(78.75) * y * y * y * y * z);
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

        pc_60[k] = e_0 * (-std::sqrt(0.205078125) * x * x * x * x * x * x * x + std::sqrt(5.126953125) * x * x * x * x * x * y * y + std::sqrt(3.28125) * x * x * x * x * x * z * z + std::sqrt(5.126953125) * x * x * x * y * y * y * y - std::sqrt(118.125) * x * x * x * y * y * z * z - std::sqrt(0.205078125) * x * y * y * y * y * y * y + std::sqrt(3.28125) * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(29.53125) * x * x * x * x * x + std::sqrt(643.125) * x * x * x * y * y + std::sqrt(52.5) * x * x * x * z * z + std::sqrt(3.28125) * x * y * y * y * y - std::sqrt(472.5) * x * y * y * z * z) + e_2 * (-std::sqrt(118.125) * x * x * x + std::sqrt(1063.125) * x * y * y);

        pc_61[k] = e_0 * (std::sqrt(2.05078125) * x * x * x * x * x * x * z - std::sqrt(100.48828125) * x * x * x * x * y * y * z + std::sqrt(100.48828125) * x * x * y * y * y * y * z - std::sqrt(2.05078125) * y * y * y * y * y * y * z) + e_1 * (std::sqrt(131.25) * x * x * x * x * z - std::sqrt(131.25) * y * y * y * y * z) + e_2 * (std::sqrt(1181.25) * x * x * z - std::sqrt(1181.25) * y * y * z);

        pc_62[k] = e_0 * (std::sqrt(0.341796875) * x * x * x * x * x * x * x - std::sqrt(27.685546875) * x * x * x * x * x * y * y + std::sqrt(123.388671875) * x * x * x * y * y * y * y - std::sqrt(3.076171875) * x * y * y * y * y * y * y) + e_1 * (std::sqrt(49.21875) * x * x * x * x * x + std::sqrt(196.875) * x * x * x * y * y + std::sqrt(49.21875) * x * y * y * y * y) + e_2 * (std::sqrt(1771.875) * x * x * x + std::sqrt(1771.875) * x * y * y) + e_3 * (std::sqrt(3150.0) * x);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 63; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
