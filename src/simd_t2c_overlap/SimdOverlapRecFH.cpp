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



#include "SimdOverlapRecFH.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_fh_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFH.compute_fh_overlap: Basis functions must be of angular momenta three and five"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFH.compute_fh_overlap: Harmonics must reach angular momentum eight"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFH.compute_fh_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    const auto nmax = dimensions.back();

    if (nmax == 0)
    {
        std::fill(values, values + 77 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = CSimdMatrix(4, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);

    const auto *ab_2 = coordinates.data(6);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each term over the pairs of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto bexp = b_exps[j];

            const auto fexp = aexp + bexp;

            const auto fmu = aexp * bexp / fexp;

            const auto fovl = fpi / fexp;

            const auto fbase = anorm * b_norms[j] * fovl * std::sqrt(fovl);

            const auto f_0 = fbase * aexp * aexp * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * aexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * aexp * aexp * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * aexp * aexp / fexp / fexp / fexp / fexp / fexp;

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
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

    const auto *ph2_m2 = harmonics[1].data(0);
    const auto *ph2_m1 = harmonics[1].data(1);
    const auto *ph2_0 = harmonics[1].data(2);
    const auto *ph2_p1 = harmonics[1].data(3);
    const auto *ph2_p2 = harmonics[1].data(4);
    const auto *ph4_m4 = harmonics[3].data(0);
    const auto *ph4_m3 = harmonics[3].data(1);
    const auto *ph4_m2 = harmonics[3].data(2);
    const auto *ph4_m1 = harmonics[3].data(3);
    const auto *ph4_0 = harmonics[3].data(4);
    const auto *ph4_p1 = harmonics[3].data(5);
    const auto *ph4_p2 = harmonics[3].data(6);
    const auto *ph4_p3 = harmonics[3].data(7);
    const auto *ph4_p4 = harmonics[3].data(8);
    const auto *ph6_m6 = harmonics[5].data(0);
    const auto *ph6_m5 = harmonics[5].data(1);
    const auto *ph6_m4 = harmonics[5].data(2);
    const auto *ph6_m3 = harmonics[5].data(3);
    const auto *ph6_m2 = harmonics[5].data(4);
    const auto *ph6_m1 = harmonics[5].data(5);
    const auto *ph6_0 = harmonics[5].data(6);
    const auto *ph6_p1 = harmonics[5].data(7);
    const auto *ph6_p2 = harmonics[5].data(8);
    const auto *ph6_p3 = harmonics[5].data(9);
    const auto *ph6_p4 = harmonics[5].data(10);
    const auto *ph6_p5 = harmonics[5].data(11);
    const auto *ph6_p6 = harmonics[5].data(12);
    const auto *ph8_m8 = harmonics[7].data(0);
    const auto *ph8_m7 = harmonics[7].data(1);
    const auto *ph8_m6 = harmonics[7].data(2);
    const auto *ph8_m5 = harmonics[7].data(3);
    const auto *ph8_m4 = harmonics[7].data(4);
    const auto *ph8_m3 = harmonics[7].data(5);
    const auto *ph8_m2 = harmonics[7].data(6);
    const auto *ph8_m1 = harmonics[7].data(7);
    const auto *ph8_0 = harmonics[7].data(8);
    const auto *ph8_p1 = harmonics[7].data(9);
    const auto *ph8_p2 = harmonics[7].data(10);
    const auto *ph8_p3 = harmonics[7].data(11);
    const auto *ph8_p4 = harmonics[7].data(12);
    const auto *ph8_p5 = harmonics[7].data(13);
    const auto *ph8_p6 = harmonics[7].data(14);
    const auto *ph8_p7 = harmonics[7].data(15);
    const auto *ph8_p8 = harmonics[7].data(16);

    // NOTE: the rows of the values are not aligned, as they start at the offset
    // of this combination of basis functions in the values block, so they are kept
    // out of the aligned clauses below.

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_2025_112 = std::sqrt(2025.0 / 112.0);
    const auto fs_30375_112 = std::sqrt(30375.0 / 112.0);
    const auto fs_75_968 = std::sqrt(75.0 / 968.0);
    const auto fs_2025_847 = std::sqrt(2025.0 / 847.0);
    const auto fs_375_28 = std::sqrt(375.0 / 28.0);
    const auto fs_1_40898 = std::sqrt(1.0 / 40898.0);
    const auto fs_28_143 = std::sqrt(28.0 / 143.0);
    const auto fs_1_726 = std::sqrt(1.0 / 726.0);
    const auto fs_2025_143143 = std::sqrt(2025.0 / 143143.0);
    const auto fs_125_2541 = std::sqrt(125.0 / 2541.0);
    const auto fs_23625_64 = std::sqrt(369.140625);
    const auto fs_3645_112 = std::sqrt(3645.0 / 112.0);
    const auto fs_6075_56 = std::sqrt(6075.0 / 56.0);
    const auto fs_75_242 = std::sqrt(75.0 / 242.0);
    const auto fs_3645_847 = std::sqrt(3645.0 / 847.0);
    const auto fs_75_14 = std::sqrt(75.0 / 14.0);
    const auto fs_7_40898 = std::sqrt(7.0 / 40898.0);
    const auto fs_35_286 = std::sqrt(35.0 / 286.0);
    const auto fs_2_363 = std::sqrt(2.0 / 363.0);
    const auto fs_3645_143143 = std::sqrt(3645.0 / 143143.0);
    const auto fs_50_2541 = std::sqrt(50.0 / 2541.0);
    const auto fs_4725_32 = std::sqrt(147.65625);
    const auto fs_2025_28 = std::sqrt(2025.0 / 28.0);
    const auto fs_175_121 = std::sqrt(175.0 / 121.0);
    const auto fs_75_88 = std::sqrt(75.0 / 88.0);
    const auto fs_8100_847 = std::sqrt(8100.0 / 847.0);
    const auto fs_25_7 = std::sqrt(25.0 / 7.0);
    const auto fs_28_20449 = std::sqrt(28.0 / 20449.0);
    const auto fs_21_286 = std::sqrt(21.0 / 286.0);
    const auto fs_28_1089 = std::sqrt(28.0 / 1089.0);
    const auto fs_1_66 = std::sqrt(1.0 / 66.0);
    const auto fs_8100_143143 = std::sqrt(8100.0 / 143143.0);
    const auto fs_100_7623 = std::sqrt(100.0 / 7623.0);
    const auto fs_1575_16 = std::sqrt(98.4375);
    const auto fs_3375_112 = std::sqrt(3375.0 / 112.0);
    const auto fs_2025_224 = std::sqrt(2025.0 / 224.0);
    const auto fs_1225_968 = std::sqrt(1225.0 / 968.0);
    const auto fs_75_44 = std::sqrt(75.0 / 44.0);
    const auto fs_3375_847 = std::sqrt(3375.0 / 847.0);
    const auto fs_25_56 = std::sqrt(25.0 / 56.0);
    const auto fs_42_20449 = std::sqrt(42.0 / 20449.0);
    const auto fs_6_143 = std::sqrt(6.0 / 143.0);
    const auto fs_49_2178 = std::sqrt(49.0 / 2178.0);
    const auto fs_1_33 = std::sqrt(1.0 / 33.0);
    const auto fs_3375_143143 = std::sqrt(3375.0 / 143143.0);
    const auto fs_25_15246 = std::sqrt(25.0 / 15246.0);
    const auto fs_1575_128 = std::sqrt(12.3046875);
    const auto fs_30375_1568 = std::sqrt(30375.0 / 1568.0);
    const auto fs_135_56 = std::sqrt(135.0 / 56.0);
    const auto fs_2025_1568 = std::sqrt(2025.0 / 1568.0);
    const auto fs_875_484 = std::sqrt(875.0 / 484.0);
    const auto fs_525_242 = std::sqrt(525.0 / 242.0);
    const auto fs_30375_11858 = std::sqrt(30375.0 / 11858.0);
    const auto fs_270_847 = std::sqrt(270.0 / 847.0);
    const auto fs_25_392 = std::sqrt(25.0 / 392.0);
    const auto fs_105_20449 = std::sqrt(105.0 / 20449.0);
    const auto fs_42_1859 = std::sqrt(42.0 / 1859.0);
    const auto fs_35_1089 = std::sqrt(35.0 / 1089.0);
    const auto fs_14_363 = std::sqrt(14.0 / 363.0);
    const auto fs_30375_2004002 = std::sqrt(30375.0 / 2004002.0);
    const auto fs_270_143143 = std::sqrt(270.0 / 143143.0);
    const auto fs_25_106722 = std::sqrt(25.0 / 106722.0);
    const auto fs_225_128 = std::sqrt(1.7578125);
    const auto fs_525_121 = std::sqrt(525.0 / 121.0);
    const auto fs_28_363 = std::sqrt(28.0 / 363.0);
    const auto fs_675_16 = std::sqrt(42.1875);
    const auto f_15_22 = 15.0 / 22.0;
    const auto fs_675_121 = std::sqrt(675.0 / 121.0);
    const auto fs_1_3718 = std::sqrt(1.0 / 3718.0);
    const auto f_1_11 = 1.0 / 11.0;
    const auto fs_675_20449 = std::sqrt(675.0 / 20449.0);
    const auto fs_135_7 = std::sqrt(135.0 / 7.0);
    const auto fs_18225_112 = std::sqrt(18225.0 / 112.0);
    const auto fs_1125_968 = std::sqrt(1125.0 / 968.0);
    const auto fs_225_88 = std::sqrt(225.0 / 88.0);
    const auto fs_2160_847 = std::sqrt(2160.0 / 847.0);
    const auto fs_225_28 = std::sqrt(225.0 / 28.0);
    const auto fs_30_20449 = std::sqrt(30.0 / 20449.0);
    const auto fs_14_143 = std::sqrt(14.0 / 143.0);
    const auto fs_5_242 = std::sqrt(5.0 / 242.0);
    const auto fs_1_22 = std::sqrt(1.0 / 22.0);
    const auto fs_2160_143143 = std::sqrt(2160.0 / 143143.0);
    const auto fs_25_847 = std::sqrt(25.0 / 847.0);
    const auto fs_14175_64 = std::sqrt(221.484375);
    const auto fs_135_112 = std::sqrt(135.0 / 112.0);
    const auto fs_2025_14 = std::sqrt(2025.0 / 14.0);
    const auto fs_200_121 = std::sqrt(200.0 / 121.0);
    const auto fs_135_847 = std::sqrt(135.0 / 847.0);
    const auto fs_50_7 = std::sqrt(50.0 / 7.0);
    const auto fs_189_40898 = std::sqrt(189.0 / 40898.0);
    const auto fs_27_286 = std::sqrt(27.0 / 286.0);
    const auto fs_32_1089 = std::sqrt(32.0 / 1089.0);
    const auto fs_135_143143 = std::sqrt(135.0 / 143143.0);
    const auto fs_200_7623 = std::sqrt(200.0 / 7623.0);
    const auto fs_1575_8 = std::sqrt(196.875);
    const auto fs_45_4 = std::sqrt(11.25);
    const auto fs_1575_484 = std::sqrt(1575.0 / 484.0);
    const auto fs_900_847 = std::sqrt(900.0 / 847.0);
    const auto fs_180_121 = std::sqrt(180.0 / 121.0);
    const auto fs_448_20449 = std::sqrt(448.0 / 20449.0);
    const auto fs_144_1859 = std::sqrt(144.0 / 1859.0);
    const auto fs_7_121 = std::sqrt(7.0 / 121.0);
    const auto fs_900_143143 = std::sqrt(900.0 / 143143.0);
    const auto fs_180_20449 = std::sqrt(180.0 / 20449.0);
    const auto fs_28125_1568 = std::sqrt(28125.0 / 1568.0);
    const auto fs_5445_224 = std::sqrt(5445.0 / 224.0);
    const auto fs_6075_196 = std::sqrt(6075.0 / 196.0);
    const auto fs_525_484 = std::sqrt(525.0 / 484.0);
    const auto fs_28125_11858 = std::sqrt(28125.0 / 11858.0);
    const auto fs_45_14 = std::sqrt(45.0 / 14.0);
    const auto fs_75_49 = std::sqrt(75.0 / 49.0);
    const auto f_21_143 = 21.0 / 143.0;
    const auto fs_105_1859 = std::sqrt(105.0 / 1859.0);
    const auto fs_7_363 = std::sqrt(7.0 / 363.0);
    const auto fs_28125_2004002 = std::sqrt(28125.0 / 2004002.0);
    const auto fs_45_2366 = std::sqrt(45.0 / 2366.0);
    const auto fs_100_17787 = std::sqrt(100.0 / 17787.0);
    const auto fs_2700_49 = std::sqrt(2700.0 / 49.0);
    const auto fs_10125_784 = std::sqrt(10125.0 / 784.0);
    const auto fs_175_242 = std::sqrt(175.0 / 242.0);
    const auto fs_43200_5929 = std::sqrt(43200.0 / 5929.0);
    const auto fs_125_196 = std::sqrt(125.0 / 196.0);
    const auto fs_1512_20449 = std::sqrt(1512.0 / 20449.0);
    const auto fs_14_1089 = std::sqrt(14.0 / 1089.0);
    const auto fs_43200_1002001 = std::sqrt(43200.0 / 1002001.0);
    const auto fs_125_53361 = std::sqrt(125.0 / 53361.0);
    const auto fs_1125_64 = std::sqrt(17.578125);
    const auto fs_135_4 = std::sqrt(33.75);
    const auto fs_675_484 = std::sqrt(675.0 / 484.0);
    const auto fs_540_121 = std::sqrt(540.0 / 121.0);
    const auto fs_3_1859 = std::sqrt(3.0 / 1859.0);
    const auto fs_7_286 = std::sqrt(7.0 / 286.0);
    const auto fs_3_121 = std::sqrt(3.0 / 121.0);
    const auto fs_540_20449 = std::sqrt(540.0 / 20449.0);
    const auto fs_27_16 = std::sqrt(1.6875);
    const auto f_15_11 = 15.0 / 11.0;
    const auto fs_27_121 = std::sqrt(27.0 / 121.0);
    const auto fs_25_3718 = std::sqrt(25.0 / 3718.0);
    const auto fs_15_286 = std::sqrt(15.0 / 286.0);
    const auto f_2_11 = 2.0 / 11.0;
    const auto fs_27_20449 = std::sqrt(27.0 / 20449.0);
    const auto fs_2187_112 = std::sqrt(2187.0 / 112.0);
    const auto fs_27 = std::sqrt(27.0);
    const auto fs_10125_112 = std::sqrt(10125.0 / 112.0);
    const auto fs_375_484 = std::sqrt(375.0 / 484.0);
    const auto fs_2187_847 = std::sqrt(2187.0 / 847.0);
    const auto fs_432_121 = std::sqrt(432.0 / 121.0);
    const auto fs_125_28 = std::sqrt(125.0 / 28.0);
    const auto fs_675_40898 = std::sqrt(675.0 / 40898.0);
    const auto fs_135_1859 = std::sqrt(135.0 / 1859.0);
    const auto fs_5_363 = std::sqrt(5.0 / 363.0);
    const auto fs_2187_143143 = std::sqrt(2187.0 / 143143.0);
    const auto fs_432_20449 = std::sqrt(432.0 / 20449.0);
    const auto fs_125_7623 = std::sqrt(125.0 / 7623.0);
    const auto fs_7875_64 = std::sqrt(123.046875);
    const auto fs_144_7 = std::sqrt(144.0 / 7.0);
    const auto f_21_4 = 5.25;
    const auto fs_30375_224 = std::sqrt(30375.0 / 224.0);
    const auto fs_375_968 = std::sqrt(375.0 / 968.0);
    const auto fs_2304_847 = std::sqrt(2304.0 / 847.0);
    const auto f_21_11 = 21.0 / 11.0;
    const auto fs_375_56 = std::sqrt(375.0 / 56.0);
    const auto fs_630_20449 = std::sqrt(630.0 / 20449.0);
    const auto fs_150_1859 = std::sqrt(150.0 / 1859.0);
    const auto fs_5_726 = std::sqrt(5.0 / 726.0);
    const auto fs_2304_143143 = std::sqrt(2304.0 / 143143.0);
    const auto fs_125_5082 = std::sqrt(125.0 / 5082.0);
    const auto fs_23625_128 = std::sqrt(184.5703125);
    const auto fs_5445_392 = std::sqrt(5445.0 / 392.0);
    const auto fs_16641_1568 = std::sqrt(16641.0 / 1568.0);
    const auto fs_91125_392 = std::sqrt(91125.0 / 392.0);
    const auto fs_90_49 = std::sqrt(90.0 / 49.0);
    const auto fs_16641_11858 = std::sqrt(16641.0 / 11858.0);
    const auto fs_1125_98 = std::sqrt(1125.0 / 98.0);
    const auto fs_375_392 = std::sqrt(375.0 / 392.0);
    const auto fs_1960_20449 = std::sqrt(1960.0 / 20449.0);
    const auto fs_1575_20449 = std::sqrt(1575.0 / 20449.0);
    const auto fs_90_8281 = std::sqrt(90.0 / 8281.0);
    const auto fs_16641_2004002 = std::sqrt(16641.0 / 2004002.0);
    const auto fs_250_5929 = std::sqrt(250.0 / 5929.0);
    const auto fs_125_35574 = std::sqrt(125.0 / 35574.0);
    const auto fs_10125_32 = std::sqrt(316.40625);
    const auto fs_3375_128 = std::sqrt(26.3671875);
    const auto fs_135_784 = std::sqrt(135.0 / 784.0);
    const auto fs_50625_392 = std::sqrt(50625.0 / 392.0);
    const auto fs_135_5929 = std::sqrt(135.0 / 5929.0);
    const auto fs_625_98 = std::sqrt(625.0 / 98.0);
    const auto fs_2646_20449 = std::sqrt(2646.0 / 20449.0);
    const auto fs_135_1002001 = std::sqrt(135.0 / 1002001.0);
    const auto fs_1250_53361 = std::sqrt(1250.0 / 53361.0);
    const auto fs_5625_32 = std::sqrt(175.78125);
    const auto fs_225_44 = std::sqrt(225.0 / 44.0);
    const auto fs_2_143 = std::sqrt(2.0 / 143.0);
    const auto fs_1_11 = std::sqrt(1.0 / 11.0);
    const auto f_9 = 9.0;
    const auto fs_1125_484 = std::sqrt(1125.0 / 484.0);
    const auto f_36_11 = 36.0 / 11.0;
    const auto fs_80_1859 = std::sqrt(80.0 / 1859.0);
    const auto fs_5_121 = std::sqrt(5.0 / 121.0);
    const auto f_36_143 = 36.0 / 143.0;
    const auto f_9_2 = 4.5;
    const auto fs_75_484 = std::sqrt(75.0 / 484.0);
    const auto f_18_11 = 18.0 / 11.0;
    const auto fs_1_363 = std::sqrt(1.0 / 363.0);
    const auto f_18_143 = 18.0 / 143.0;
    const auto fs_27_28 = std::sqrt(27.0 / 28.0);
    const auto fs_50_121 = std::sqrt(50.0 / 121.0);
    const auto fs_108_847 = std::sqrt(108.0 / 847.0);
    const auto fs_2400_20449 = std::sqrt(2400.0 / 20449.0);
    const auto fs_8_1089 = std::sqrt(8.0 / 1089.0);
    const auto fs_108_143143 = std::sqrt(108.0 / 143143.0);
    const auto fs_9747_392 = std::sqrt(9747.0 / 392.0);
    const auto fs_10125_49 = std::sqrt(10125.0 / 49.0);
    const auto fs_19494_5929 = std::sqrt(19494.0 / 5929.0);
    const auto fs_500_49 = std::sqrt(500.0 / 49.0);
    const auto fs_2940_20449 = std::sqrt(2940.0 / 20449.0);
    const auto fs_19494_1002001 = std::sqrt(19494.0 / 1002001.0);
    const auto fs_2000_53361 = std::sqrt(2000.0 / 53361.0);
    const auto fs_1125_4 = std::sqrt(281.25);
    const auto f_45_7 = 45.0 / 7.0;
    const auto f_225_14 = 225.0 / 14.0;
    const auto f_35_22 = 35.0 / 22.0;
    const auto f_180_77 = 180.0 / 77.0;
    const auto f_25_7 = 25.0 / 7.0;
    const auto f_56_143 = 56.0 / 143.0;
    const auto f_7_33 = 7.0 / 33.0;
    const auto f_180_1001 = 180.0 / 1001.0;
    const auto f_50_231 = 50.0 / 231.0;
    const auto f_75_4 = 18.75;

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph8_p7, ph8_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * (fs_2025_112 * h4_p2 - fs_30375_112 * r_2 * h2_p2) + e_1 * (fs_75_968 * h6_p2 - fs_2025_847 * r_2 * h4_p2 + fs_375_28 * r_4 * h2_p2) + e_2 * (fs_1_40898 * h8_p2 + fs_28_143 * h8_p8 - fs_1_726 * r_2 * h6_p2 + fs_2025_143143 * r_4 * h4_p2 - fs_125_2541 * r_6 * h2_p2) + fs_23625_64 * e_3 * h2_p2;

        pc_1[k] = e_0 * (fs_3645_112 * h4_p1 - fs_6075_56 * r_2 * h2_p1) + e_1 * (fs_75_242 * h6_p1 - fs_3645_847 * r_2 * h4_p1 + fs_75_14 * r_4 * h2_p1) + e_2 * (fs_7_40898 * h8_p1 + fs_35_286 * h8_p7 - fs_2_363 * r_2 * h6_p1 + fs_3645_143143 * r_4 * h4_p1 - fs_50_2541 * r_6 * h2_p1) + fs_4725_32 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_2[k] = e_0 * (fs_2025_28 * h4_0 - fs_2025_28 * r_2 * h2_0) + e_1 * (fs_175_121 * h6_0 + fs_75_88 * h6_p6 - fs_8100_847 * r_2 * h4_0 + fs_25_7 * r_4 * h2_0) + e_2 * (fs_28_20449 * h8_0 + fs_21_286 * h8_p6 - fs_28_1089 * r_2 * h6_0 - fs_1_66 * r_2 * h6_p6 + fs_8100_143143 * r_4 * h4_0 - fs_100_7623 * r_6 * h2_0) + fs_1575_16 * e_3 * h2_0;

        pc_3[k] = e_0 * (-fs_3375_112 * h4_p1 + fs_2025_224 * r_2 * h2_p1) + e_1 * (-fs_1225_968 * h6_p1 + fs_75_44 * h6_p5 + fs_3375_847 * r_2 * h4_p1 - fs_25_56 * r_4 * h2_p1) + e_2 * (-fs_42_20449 * h8_p1 + fs_6_143 * h8_p5 + fs_49_2178 * r_2 * h6_p1 - fs_1_33 * r_2 * h6_p5 - fs_3375_143143 * r_4 * h4_p1 + fs_25_15246 * r_6 * h2_p1) - fs_1575_128 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_m3, ph4_p2, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_4[k] = e_0 * (fs_30375_1568 * h4_p2 + fs_135_56 * h4_p4 - fs_2025_1568 * r_2 * h2_p2) + e_1 * (fs_875_484 * h6_p2 + fs_525_242 * h6_p4 - fs_30375_11858 * r_2 * h4_p2 - fs_270_847 * r_2 * h4_p4 + fs_25_392 * r_4 * h2_p2) + e_2 * (fs_105_20449 * h8_p2 + fs_42_1859 * h8_p4 - fs_35_1089 * r_2 * h6_p2 - fs_14_363 * r_2 * h6_p4 + fs_30375_2004002 * r_4 * h4_p2 + fs_270_143143 * r_4 * h4_p4 - fs_25_106722 * r_6 * h2_p2) + fs_225_128 * e_3 * h2_p2;

        pc_5[k] = -fs_2025_112 * e_0 * h4_m3 + e_1 * (-fs_525_121 * h6_m3 + fs_2025_847 * r_2 * h4_m3) + e_2 * (-fs_42_1859 * h8_m3 + fs_28_363 * r_2 * h6_m3 - fs_2025_143143 * r_4 * h4_m3);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_6[k] = e_0 * (-fs_135_56 * h4_m4 + fs_30375_1568 * h4_m2 - fs_2025_1568 * r_2 * h2_m2) + e_1 * (-fs_525_242 * h6_m4 + fs_875_484 * h6_m2 + fs_270_847 * r_2 * h4_m4 - fs_30375_11858 * r_2 * h4_m2 + fs_25_392 * r_4 * h2_m2) + e_2 * (-fs_42_1859 * h8_m4 + fs_105_20449 * h8_m2 + fs_14_363 * r_2 * h6_m4 - fs_35_1089 * r_2 * h6_m2 - fs_270_143143 * r_4 * h4_m4 + fs_30375_2004002 * r_4 * h4_m2 - fs_25_106722 * r_6 * h2_m2) + fs_225_128 * e_3 * h2_m2;

        pc_7[k] = e_0 * (-fs_3375_112 * h4_m1 + fs_2025_224 * r_2 * h2_m1) + e_1 * (-fs_75_44 * h6_m5 - fs_1225_968 * h6_m1 + fs_3375_847 * r_2 * h4_m1 - fs_25_56 * r_4 * h2_m1) + e_2 * (-fs_6_143 * h8_m5 - fs_42_20449 * h8_m1 + fs_1_33 * r_2 * h6_m5 + fs_49_2178 * r_2 * h6_m1 - fs_3375_143143 * r_4 * h4_m1 + fs_25_15246 * r_6 * h2_m1) - fs_1575_128 * e_3 * h2_m1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_8[k] = -fs_75_88 * e_1 * h6_m6 + e_2 * (-fs_21_286 * h8_m6 + fs_1_66 * r_2 * h6_m6);

        pc_9[k] = e_0 * (-fs_3645_112 * h4_m1 + fs_6075_56 * r_2 * h2_m1) + e_1 * (-fs_75_242 * h6_m1 + fs_3645_847 * r_2 * h4_m1 - fs_75_14 * r_4 * h2_m1) + e_2 * (-fs_35_286 * h8_m7 - fs_7_40898 * h8_m1 + fs_2_363 * r_2 * h6_m1 - fs_3645_143143 * r_4 * h4_m1 + fs_50_2541 * r_6 * h2_m1) - fs_4725_32 * e_3 * h2_m1;

        pc_10[k] = e_0 * (-fs_2025_112 * h4_m2 + fs_30375_112 * r_2 * h2_m2) + e_1 * (-fs_75_968 * h6_m2 + fs_2025_847 * r_2 * h4_m2 - fs_375_28 * r_4 * h2_m2) + e_2 * (-fs_28_143 * h8_m8 - fs_1_40898 * h8_m2 + fs_1_726 * r_2 * h6_m2 - fs_2025_143143 * r_4 * h4_m2 + fs_125_2541 * r_6 * h2_m2) - fs_23625_64 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p6, ph8_p2, ph8_p3, ph8_p6, ph8_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_11[k] = -fs_675_16 * e_0 * h4_p3 + e_1 * (-f_15_22 * h6_p3 + fs_675_121 * r_2 * h4_p3) + e_2 * (-fs_1_3718 * h8_p3 + fs_21_286 * h8_p7 + f_1_11 * r_2 * h6_p3 - fs_675_20449 * r_4 * h4_p3);

        pc_12[k] = e_0 * (-fs_135_7 * h4_p2 - fs_18225_112 * r_2 * h2_p2) + e_1 * (-fs_1125_968 * h6_p2 - fs_225_88 * h6_p6 + fs_2160_847 * r_2 * h4_p2 + fs_225_28 * r_4 * h2_p2) + e_2 * (-fs_30_20449 * h8_p2 + fs_14_143 * h8_p6 + fs_5_242 * r_2 * h6_p2 + fs_1_22 * r_2 * h6_p6 - fs_2160_143143 * r_4 * h4_p2 - fs_25_847 * r_6 * h2_p2) + fs_14175_64 * e_3 * h2_p2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph4_p4, ph6_0, ph6_p1, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p4, ph8_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_13[k] = e_0 * (-fs_135_112 * h4_p1 - fs_2025_14 * r_2 * h2_p1) + e_1 * (-fs_200_121 * h6_p1 - fs_75_44 * h6_p5 + fs_135_847 * r_2 * h4_p1 + fs_50_7 * r_4 * h2_p1) + e_2 * (-fs_189_40898 * h8_p1 + fs_27_286 * h8_p5 + fs_32_1089 * r_2 * h6_p1 + fs_1_33 * r_2 * h6_p5 - fs_135_143143 * r_4 * h4_p1 - fs_200_7623 * r_6 * h2_p1) + fs_1575_8 * e_3 * h2_p1;

        pc_14[k] = e_0 * (fs_225_28 * h4_0 - fs_45_4 * h4_p4 - fs_18225_112 * r_2 * h2_0) + e_1 * (-fs_1575_484 * h6_0 - f_15_22 * h6_p4 - fs_900_847 * r_2 * h4_0 + fs_180_121 * r_2 * h4_p4 + fs_225_28 * r_4 * h2_0) + e_2 * (-fs_448_20449 * h8_0 + fs_144_1859 * h8_p4 + fs_7_121 * r_2 * h6_0 + f_1_11 * r_2 * h6_p4 + fs_900_143143 * r_4 * h4_0 - fs_180_20449 * r_4 * h4_p4 - fs_25_847 * r_6 * h2_0) + fs_14175_64 * e_3 * h2_0;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph4_p3, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph8_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_15[k] = e_0 * (-fs_28125_1568 * h4_p1 - fs_5445_224 * h4_p3 + fs_6075_196 * r_2 * h2_p1) + e_1 * (fs_525_484 * h6_p1 + fs_28125_11858 * r_2 * h4_p1 + fs_45_14 * r_2 * h4_p3 - fs_75_49 * r_4 * h2_p1) + e_2 * (f_21_143 * h8_p1 + fs_105_1859 * h8_p3 - fs_7_363 * r_2 * h6_p1 - fs_28125_2004002 * r_4 * h4_p1 - fs_45_2366 * r_4 * h4_p3 + fs_100_17787 * r_6 * h2_p1) - fs_675_16 * e_3 * h2_p1;

        pc_16[k] = e_0 * (fs_2700_49 * h4_m2 - fs_10125_784 * r_2 * h2_m2) + e_1 * (-fs_175_242 * h6_m2 - fs_43200_5929 * r_2 * h4_m2 + fs_125_196 * r_4 * h2_m2) + e_2 * (-fs_1512_20449 * h8_m2 + fs_14_1089 * r_2 * h6_m2 + fs_43200_1002001 * r_4 * h4_m2 - fs_125_53361 * r_6 * h2_m2) + fs_1125_64 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_17[k] = e_0 * (fs_5445_224 * h4_m3 - fs_28125_1568 * h4_m1 + fs_6075_196 * r_2 * h2_m1) + e_1 * (fs_525_484 * h6_m1 - fs_45_14 * r_2 * h4_m3 + fs_28125_11858 * r_2 * h4_m1 - fs_75_49 * r_4 * h2_m1) + e_2 * (-fs_105_1859 * h8_m3 + f_21_143 * h8_m1 - fs_7_363 * r_2 * h6_m1 + fs_45_2366 * r_4 * h4_m3 - fs_28125_2004002 * r_4 * h4_m1 + fs_100_17787 * r_6 * h2_m1) - fs_675_16 * e_3 * h2_m1;

        pc_18[k] = fs_45_4 * e_0 * h4_m4 + e_1 * (f_15_22 * h6_m4 - fs_180_121 * r_2 * h4_m4) + e_2 * (-fs_144_1859 * h8_m4 - f_1_11 * r_2 * h6_m4 + fs_180_20449 * r_4 * h4_m4);

        pc_19[k] = e_0 * (fs_135_112 * h4_m1 + fs_2025_14 * r_2 * h2_m1) + e_1 * (fs_75_44 * h6_m5 + fs_200_121 * h6_m1 - fs_135_847 * r_2 * h4_m1 - fs_50_7 * r_4 * h2_m1) + e_2 * (-fs_27_286 * h8_m5 + fs_189_40898 * h8_m1 - fs_1_33 * r_2 * h6_m5 - fs_32_1089 * r_2 * h6_m1 + fs_135_143143 * r_4 * h4_m1 + fs_200_7623 * r_6 * h2_m1) - fs_1575_8 * e_3 * h2_m1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m3, ph4_m2, ph6_m6, ph6_m3, ph6_m2, ph8_m7, ph8_m6, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_20[k] = e_0 * (fs_135_7 * h4_m2 + fs_18225_112 * r_2 * h2_m2) + e_1 * (fs_225_88 * h6_m6 + fs_1125_968 * h6_m2 - fs_2160_847 * r_2 * h4_m2 - fs_225_28 * r_4 * h2_m2) + e_2 * (-fs_14_143 * h8_m6 + fs_30_20449 * h8_m2 - fs_1_22 * r_2 * h6_m6 - fs_5_242 * r_2 * h6_m2 + fs_2160_143143 * r_4 * h4_m2 + fs_25_847 * r_6 * h2_m2) - fs_14175_64 * e_3 * h2_m2;

        pc_21[k] = fs_675_16 * e_0 * h4_m3 + e_1 * (f_15_22 * h6_m3 - fs_675_121 * r_2 * h4_m3) + e_2 * (-fs_21_286 * h8_m7 + fs_1_3718 * h8_m3 - f_1_11 * r_2 * h6_m3 + fs_675_20449 * r_4 * h4_m3);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_22[k] = fs_135_4 * e_0 * h4_p4 + e_1 * (fs_675_484 * h6_p4 + fs_225_88 * h6_p6 - fs_540_121 * r_2 * h4_p4) + e_2 * (fs_3_1859 * h8_p4 + fs_7_286 * h8_p6 - fs_3_121 * r_2 * h6_p4 - fs_1_22 * r_2 * h6_p6 + fs_540_20449 * r_4 * h4_p4);

        pc_23[k] = -fs_27_16 * e_0 * h4_p3 + e_1 * (f_15_11 * h6_p3 + fs_27_121 * r_2 * h4_p3) + e_2 * (fs_25_3718 * h8_p3 + fs_15_286 * h8_p5 - f_2_11 * r_2 * h6_p3 - fs_27_20449 * r_4 * h4_p3);

        pc_24[k] = e_0 * (-fs_2187_112 * h4_p2 + fs_27 * h4_p4 - fs_10125_112 * r_2 * h2_p2) + e_1 * (fs_1225_968 * h6_p2 - fs_375_484 * h6_p4 + fs_2187_847 * r_2 * h4_p2 - fs_432_121 * r_2 * h4_p4 + fs_125_28 * r_4 * h2_p2) + e_2 * (fs_675_40898 * h8_p2 + fs_135_1859 * h8_p4 - fs_49_2178 * r_2 * h6_p2 + fs_5_363 * r_2 * h6_p4 - fs_2187_143143 * r_4 * h4_p2 + fs_432_20449 * r_4 * h4_p4 - fs_125_7623 * r_6 * h2_p2) + fs_7875_64 * e_3 * h2_p2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_25[k] = e_0 * (-fs_144_7 * h4_p1 + f_21_4 * h4_p3 - fs_30375_224 * r_2 * h2_p1) + e_1 * (fs_375_968 * h6_p1 - fs_675_484 * h6_p3 + fs_2304_847 * r_2 * h4_p1 - f_21_11 * r_2 * h4_p3 + fs_375_56 * r_4 * h2_p1) + e_2 * (fs_630_20449 * h8_p1 + fs_150_1859 * h8_p3 - fs_5_726 * r_2 * h6_p1 + fs_3_121 * r_2 * h6_p3 - fs_2304_143143 * r_4 * h4_p1 + f_21_143 * r_4 * h4_p3 - fs_125_5082 * r_6 * h2_p1) + fs_23625_128 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph2_0, ph2_p2, ph4_m1, ph4_0, ph4_p2, ph6_m1, ph6_p2, ph8_m1, ph8_0, ph8_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_26[k] = e_0 * (-fs_5445_392 * h4_0 + fs_16641_1568 * h4_p2 - fs_91125_392 * r_2 * h2_0 - fs_30375_1568 * r_2 * h2_p2) + e_1 * (-fs_525_484 * h6_p2 + fs_90_49 * r_2 * h4_0 - fs_16641_11858 * r_2 * h4_p2 + fs_1125_98 * r_4 * h2_0 + fs_375_392 * r_4 * h2_p2) + e_2 * (fs_1960_20449 * h8_0 + fs_1575_20449 * h8_p2 + fs_7_363 * r_2 * h6_p2 - fs_90_8281 * r_4 * h4_0 + fs_16641_2004002 * r_4 * h4_p2 - fs_250_5929 * r_6 * h2_0 - fs_125_35574 * r_6 * h2_p2) + e_3 * (fs_10125_32 * h2_0 + fs_3375_128 * h2_p2);

        pc_27[k] = e_0 * (-fs_135_784 * h4_m1 + fs_50625_392 * r_2 * h2_m1) + e_1 * (fs_175_242 * h6_m1 + fs_135_5929 * r_2 * h4_m1 - fs_625_98 * r_4 * h2_m1) + e_2 * (-fs_2646_20449 * h8_m1 - fs_14_1089 * r_2 * h6_m1 - fs_135_1002001 * r_4 * h4_m1 + fs_1250_53361 * r_6 * h2_m1) - fs_5625_32 * e_3 * h2_m1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_28[k] = e_0 * (-fs_16641_1568 * h4_m2 + fs_30375_1568 * r_2 * h2_m2) + e_1 * (fs_525_484 * h6_m2 + fs_16641_11858 * r_2 * h4_m2 - fs_375_392 * r_4 * h2_m2) + e_2 * (-fs_1575_20449 * h8_m2 - fs_7_363 * r_2 * h6_m2 - fs_16641_2004002 * r_4 * h4_m2 + fs_125_35574 * r_6 * h2_m2) - fs_3375_128 * e_3 * h2_m2;

        pc_29[k] = e_0 * (-f_21_4 * h4_m3 + fs_144_7 * h4_m1 + fs_30375_224 * r_2 * h2_m1) + e_1 * (fs_675_484 * h6_m3 - fs_375_968 * h6_m1 + f_21_11 * r_2 * h4_m3 - fs_2304_847 * r_2 * h4_m1 - fs_375_56 * r_4 * h2_m1) + e_2 * (-fs_150_1859 * h8_m3 - fs_630_20449 * h8_m1 - fs_3_121 * r_2 * h6_m3 + fs_5_726 * r_2 * h6_m1 - f_21_143 * r_4 * h4_m3 + fs_2304_143143 * r_4 * h4_m1 + fs_125_5082 * r_6 * h2_m1) - fs_23625_128 * e_3 * h2_m1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m4, ph6_m3, ph6_m2, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_30[k] = e_0 * (-fs_27 * h4_m4 + fs_2187_112 * h4_m2 + fs_10125_112 * r_2 * h2_m2) + e_1 * (fs_375_484 * h6_m4 - fs_1225_968 * h6_m2 + fs_432_121 * r_2 * h4_m4 - fs_2187_847 * r_2 * h4_m2 - fs_125_28 * r_4 * h2_m2) + e_2 * (-fs_135_1859 * h8_m4 - fs_675_40898 * h8_m2 - fs_5_363 * r_2 * h6_m4 + fs_49_2178 * r_2 * h6_m2 - fs_432_20449 * r_4 * h4_m4 + fs_2187_143143 * r_4 * h4_m2 + fs_125_7623 * r_6 * h2_m2) - fs_7875_64 * e_3 * h2_m2;

        pc_31[k] = fs_27_16 * e_0 * h4_m3 + e_1 * (-f_15_11 * h6_m3 - fs_27_121 * r_2 * h4_m3) + e_2 * (-fs_15_286 * h8_m5 - fs_25_3718 * h8_m3 + f_2_11 * r_2 * h6_m3 + fs_27_20449 * r_4 * h4_m3);

        pc_32[k] = -fs_135_4 * e_0 * h4_m4 + e_1 * (-fs_225_88 * h6_m6 - fs_675_484 * h6_m4 + fs_540_121 * r_2 * h4_m4) + e_2 * (-fs_7_286 * h8_m6 - fs_3_1859 * h8_m4 + fs_1_22 * r_2 * h6_m6 + fs_3_121 * r_2 * h6_m4 - fs_540_20449 * r_4 * h4_m4);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_33[k] = -fs_225_44 * e_1 * h6_m5 + e_2 * (-fs_2_143 * h8_m5 + fs_1_11 * r_2 * h6_m5);

        pc_34[k] = f_9 * e_0 * h4_m4 + e_1 * (-fs_1125_484 * h6_m4 - f_36_11 * r_2 * h4_m4) + e_2 * (-fs_80_1859 * h8_m4 + fs_5_121 * r_2 * h6_m4 + f_36_143 * r_4 * h4_m4);

        pc_35[k] = f_9_2 * e_0 * h4_m3 + e_1 * (-fs_75_484 * h6_m3 - f_18_11 * r_2 * h4_m3) + e_2 * (-fs_150_1859 * h8_m3 + fs_1_363 * r_2 * h6_m3 + f_18_143 * r_4 * h4_m3);

        pc_36[k] = e_0 * (-fs_27_28 * h4_m2 - fs_10125_112 * r_2 * h2_m2) + e_1 * (fs_50_121 * h6_m2 + fs_108_847 * r_2 * h4_m2 + fs_125_28 * r_4 * h2_m2) + e_2 * (-fs_2400_20449 * h8_m2 - fs_8_1089 * r_2 * h6_m2 - fs_108_143143 * r_4 * h4_m2 - fs_125_7623 * r_6 * h2_m2) + fs_7875_64 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph2_0, ph2_p1, ph4_m1, ph4_0, ph4_p1, ph6_m1, ph6_0, ph6_p1, ph8_m1, ph8_0, ph8_p1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];

        pc_37[k] = e_0 * (-fs_9747_392 * h4_m1 - fs_10125_49 * r_2 * h2_m1) + e_1 * (fs_875_484 * h6_m1 + fs_19494_5929 * r_2 * h4_m1 + fs_500_49 * r_4 * h2_m1) + e_2 * (-fs_2940_20449 * h8_m1 - fs_35_1089 * r_2 * h6_m1 - fs_19494_1002001 * r_4 * h4_m1 - fs_2000_53361 * r_6 * h2_m1) + fs_1125_4 * e_3 * h2_m1;

        pc_38[k] = e_0 * (-f_45_7 * h4_0 - f_225_14 * r_2 * h2_0) + e_1 * (f_35_22 * h6_0 + f_180_77 * r_2 * h4_0 + f_25_7 * r_4 * h2_0) + e_2 * (-f_56_143 * h8_0 - f_7_33 * r_2 * h6_0 - f_180_1001 * r_4 * h4_0 - f_50_231 * r_6 * h2_0) + f_75_4 * e_3 * h2_0;

        pc_39[k] = e_0 * (-fs_9747_392 * h4_p1 - fs_10125_49 * r_2 * h2_p1) + e_1 * (fs_875_484 * h6_p1 + fs_19494_5929 * r_2 * h4_p1 + fs_500_49 * r_4 * h2_p1) + e_2 * (-fs_2940_20449 * h8_p1 - fs_35_1089 * r_2 * h6_p1 - fs_19494_1002001 * r_4 * h4_p1 - fs_2000_53361 * r_6 * h2_p1) + fs_1125_4 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_40[k] = e_0 * (-fs_27_28 * h4_p2 - fs_10125_112 * r_2 * h2_p2) + e_1 * (fs_50_121 * h6_p2 + fs_108_847 * r_2 * h4_p2 + fs_125_28 * r_4 * h2_p2) + e_2 * (-fs_2400_20449 * h8_p2 - fs_8_1089 * r_2 * h6_p2 - fs_108_143143 * r_4 * h4_p2 - fs_125_7623 * r_6 * h2_p2) + fs_7875_64 * e_3 * h2_p2;

        pc_41[k] = f_9_2 * e_0 * h4_p3 + e_1 * (-fs_75_484 * h6_p3 - f_18_11 * r_2 * h4_p3) + e_2 * (-fs_150_1859 * h8_p3 + fs_1_363 * r_2 * h6_p3 + f_18_143 * r_4 * h4_p3);

        pc_42[k] = f_9 * e_0 * h4_p4 + e_1 * (-fs_1125_484 * h6_p4 - f_36_11 * r_2 * h4_p4) + e_2 * (-fs_80_1859 * h8_p4 + fs_5_121 * r_2 * h6_p4 + f_36_143 * r_4 * h4_p4);

        pc_43[k] = -fs_225_44 * e_1 * h6_p5 + e_2 * (-fs_2_143 * h8_p5 + fs_1_11 * r_2 * h6_p5);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m4, ph6_m3, ph6_m2, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_44[k] = fs_135_4 * e_0 * h4_m4 + e_1 * (-fs_225_88 * h6_m6 + fs_675_484 * h6_m4 - fs_540_121 * r_2 * h4_m4) + e_2 * (-fs_7_286 * h8_m6 + fs_3_1859 * h8_m4 + fs_1_22 * r_2 * h6_m6 - fs_3_121 * r_2 * h6_m4 + fs_540_20449 * r_4 * h4_m4);

        pc_45[k] = -fs_27_16 * e_0 * h4_m3 + e_1 * (f_15_11 * h6_m3 + fs_27_121 * r_2 * h4_m3) + e_2 * (-fs_15_286 * h8_m5 + fs_25_3718 * h8_m3 - f_2_11 * r_2 * h6_m3 - fs_27_20449 * r_4 * h4_m3);

        pc_46[k] = e_0 * (-fs_27 * h4_m4 - fs_2187_112 * h4_m2 - fs_10125_112 * r_2 * h2_m2) + e_1 * (fs_375_484 * h6_m4 + fs_1225_968 * h6_m2 + fs_432_121 * r_2 * h4_m4 + fs_2187_847 * r_2 * h4_m2 + fs_125_28 * r_4 * h2_m2) + e_2 * (-fs_135_1859 * h8_m4 + fs_675_40898 * h8_m2 - fs_5_363 * r_2 * h6_m4 - fs_49_2178 * r_2 * h6_m2 - fs_432_20449 * r_4 * h4_m4 - fs_2187_143143 * r_4 * h4_m2 - fs_125_7623 * r_6 * h2_m2) + fs_7875_64 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_47[k] = e_0 * (-f_21_4 * h4_m3 - fs_144_7 * h4_m1 - fs_30375_224 * r_2 * h2_m1) + e_1 * (fs_675_484 * h6_m3 + fs_375_968 * h6_m1 + f_21_11 * r_2 * h4_m3 + fs_2304_847 * r_2 * h4_m1 + fs_375_56 * r_4 * h2_m1) + e_2 * (-fs_150_1859 * h8_m3 + fs_630_20449 * h8_m1 - fs_3_121 * r_2 * h6_m3 - fs_5_726 * r_2 * h6_m1 - f_21_143 * r_4 * h4_m3 - fs_2304_143143 * r_4 * h4_m1 - fs_125_5082 * r_6 * h2_m1) + fs_23625_128 * e_3 * h2_m1;

        pc_48[k] = e_0 * (-fs_16641_1568 * h4_m2 + fs_30375_1568 * r_2 * h2_m2) + e_1 * (fs_525_484 * h6_m2 + fs_16641_11858 * r_2 * h4_m2 - fs_375_392 * r_4 * h2_m2) + e_2 * (-fs_1575_20449 * h8_m2 - fs_7_363 * r_2 * h6_m2 - fs_16641_2004002 * r_4 * h4_m2 + fs_125_35574 * r_6 * h2_m2) - fs_3375_128 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_0, ph8_p1, ph8_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_49[k] = e_0 * (-fs_135_784 * h4_p1 + fs_50625_392 * r_2 * h2_p1) + e_1 * (fs_175_242 * h6_p1 + fs_135_5929 * r_2 * h4_p1 - fs_625_98 * r_4 * h2_p1) + e_2 * (-fs_2646_20449 * h8_p1 - fs_14_1089 * r_2 * h6_p1 - fs_135_1002001 * r_4 * h4_p1 + fs_1250_53361 * r_6 * h2_p1) - fs_5625_32 * e_3 * h2_p1;

        pc_50[k] = e_0 * (-fs_5445_392 * h4_0 - fs_16641_1568 * h4_p2 - fs_91125_392 * r_2 * h2_0 + fs_30375_1568 * r_2 * h2_p2) + e_1 * (fs_525_484 * h6_p2 + fs_90_49 * r_2 * h4_0 + fs_16641_11858 * r_2 * h4_p2 + fs_1125_98 * r_4 * h2_0 - fs_375_392 * r_4 * h2_p2) + e_2 * (fs_1960_20449 * h8_0 - fs_1575_20449 * h8_p2 - fs_7_363 * r_2 * h6_p2 - fs_90_8281 * r_4 * h4_0 - fs_16641_2004002 * r_4 * h4_p2 - fs_250_5929 * r_6 * h2_0 + fs_125_35574 * r_6 * h2_p2) + e_3 * (fs_10125_32 * h2_0 - fs_3375_128 * h2_p2);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_51[k] = e_0 * (-fs_144_7 * h4_p1 - f_21_4 * h4_p3 - fs_30375_224 * r_2 * h2_p1) + e_1 * (fs_375_968 * h6_p1 + fs_675_484 * h6_p3 + fs_2304_847 * r_2 * h4_p1 + f_21_11 * r_2 * h4_p3 + fs_375_56 * r_4 * h2_p1) + e_2 * (fs_630_20449 * h8_p1 - fs_150_1859 * h8_p3 - fs_5_726 * r_2 * h6_p1 - fs_3_121 * r_2 * h6_p3 - fs_2304_143143 * r_4 * h4_p1 - f_21_143 * r_4 * h4_p3 - fs_125_5082 * r_6 * h2_p1) + fs_23625_128 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_52[k] = e_0 * (-fs_2187_112 * h4_p2 - fs_27 * h4_p4 - fs_10125_112 * r_2 * h2_p2) + e_1 * (fs_1225_968 * h6_p2 + fs_375_484 * h6_p4 + fs_2187_847 * r_2 * h4_p2 + fs_432_121 * r_2 * h4_p4 + fs_125_28 * r_4 * h2_p2) + e_2 * (fs_675_40898 * h8_p2 - fs_135_1859 * h8_p4 - fs_49_2178 * r_2 * h6_p2 - fs_5_363 * r_2 * h6_p4 - fs_2187_143143 * r_4 * h4_p2 - fs_432_20449 * r_4 * h4_p4 - fs_125_7623 * r_6 * h2_p2) + fs_7875_64 * e_3 * h2_p2;

        pc_53[k] = -fs_27_16 * e_0 * h4_p3 + e_1 * (f_15_11 * h6_p3 + fs_27_121 * r_2 * h4_p3) + e_2 * (fs_25_3718 * h8_p3 - fs_15_286 * h8_p5 - f_2_11 * r_2 * h6_p3 - fs_27_20449 * r_4 * h4_p3);

        pc_54[k] = fs_135_4 * e_0 * h4_p4 + e_1 * (fs_675_484 * h6_p4 - fs_225_88 * h6_p6 - fs_540_121 * r_2 * h4_p4) + e_2 * (fs_3_1859 * h8_p4 - fs_7_286 * h8_p6 - fs_3_121 * r_2 * h6_p4 + fs_1_22 * r_2 * h6_p6 + fs_540_20449 * r_4 * h4_p4);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m3, ph4_m2, ph6_m6, ph6_m3, ph6_m2, ph8_m7, ph8_m6, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_55[k] = -fs_675_16 * e_0 * h4_m3 + e_1 * (-f_15_22 * h6_m3 + fs_675_121 * r_2 * h4_m3) + e_2 * (-fs_21_286 * h8_m7 - fs_1_3718 * h8_m3 + f_1_11 * r_2 * h6_m3 - fs_675_20449 * r_4 * h4_m3);

        pc_56[k] = e_0 * (-fs_135_7 * h4_m2 - fs_18225_112 * r_2 * h2_m2) + e_1 * (fs_225_88 * h6_m6 - fs_1125_968 * h6_m2 + fs_2160_847 * r_2 * h4_m2 + fs_225_28 * r_4 * h2_m2) + e_2 * (-fs_14_143 * h8_m6 - fs_30_20449 * h8_m2 - fs_1_22 * r_2 * h6_m6 + fs_5_242 * r_2 * h6_m2 - fs_2160_143143 * r_4 * h4_m2 - fs_25_847 * r_6 * h2_m2) + fs_14175_64 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_57[k] = e_0 * (-fs_135_112 * h4_m1 - fs_2025_14 * r_2 * h2_m1) + e_1 * (fs_75_44 * h6_m5 - fs_200_121 * h6_m1 + fs_135_847 * r_2 * h4_m1 + fs_50_7 * r_4 * h2_m1) + e_2 * (-fs_27_286 * h8_m5 - fs_189_40898 * h8_m1 - fs_1_33 * r_2 * h6_m5 + fs_32_1089 * r_2 * h6_m1 - fs_135_143143 * r_4 * h4_m1 - fs_200_7623 * r_6 * h2_m1) + fs_1575_8 * e_3 * h2_m1;

        pc_58[k] = fs_45_4 * e_0 * h4_m4 + e_1 * (f_15_22 * h6_m4 - fs_180_121 * r_2 * h4_m4) + e_2 * (-fs_144_1859 * h8_m4 - f_1_11 * r_2 * h6_m4 + fs_180_20449 * r_4 * h4_m4);

        pc_59[k] = e_0 * (fs_5445_224 * h4_m3 + fs_28125_1568 * h4_m1 - fs_6075_196 * r_2 * h2_m1) + e_1 * (-fs_525_484 * h6_m1 - fs_45_14 * r_2 * h4_m3 - fs_28125_11858 * r_2 * h4_m1 + fs_75_49 * r_4 * h2_m1) + e_2 * (-fs_105_1859 * h8_m3 - f_21_143 * h8_m1 + fs_7_363 * r_2 * h6_m1 + fs_45_2366 * r_4 * h4_m3 + fs_28125_2004002 * r_4 * h4_m1 - fs_100_17787 * r_6 * h2_m1) + fs_675_16 * e_3 * h2_m1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph4_p3, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph8_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_60[k] = e_0 * (fs_2700_49 * h4_p2 - fs_10125_784 * r_2 * h2_p2) + e_1 * (-fs_175_242 * h6_p2 - fs_43200_5929 * r_2 * h4_p2 + fs_125_196 * r_4 * h2_p2) + e_2 * (-fs_1512_20449 * h8_p2 + fs_14_1089 * r_2 * h6_p2 + fs_43200_1002001 * r_4 * h4_p2 - fs_125_53361 * r_6 * h2_p2) + fs_1125_64 * e_3 * h2_p2;

        pc_61[k] = e_0 * (-fs_28125_1568 * h4_p1 + fs_5445_224 * h4_p3 + fs_6075_196 * r_2 * h2_p1) + e_1 * (fs_525_484 * h6_p1 + fs_28125_11858 * r_2 * h4_p1 - fs_45_14 * r_2 * h4_p3 - fs_75_49 * r_4 * h2_p1) + e_2 * (f_21_143 * h8_p1 - fs_105_1859 * h8_p3 - fs_7_363 * r_2 * h6_p1 - fs_28125_2004002 * r_4 * h4_p1 + fs_45_2366 * r_4 * h4_p3 + fs_100_17787 * r_6 * h2_p1) - fs_675_16 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph4_p4, ph6_0, ph6_p1, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p4, ph8_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_62[k] = e_0 * (fs_225_28 * h4_0 + fs_45_4 * h4_p4 - fs_18225_112 * r_2 * h2_0) + e_1 * (-fs_1575_484 * h6_0 + f_15_22 * h6_p4 - fs_900_847 * r_2 * h4_0 - fs_180_121 * r_2 * h4_p4 + fs_225_28 * r_4 * h2_0) + e_2 * (-fs_448_20449 * h8_0 - fs_144_1859 * h8_p4 + fs_7_121 * r_2 * h6_0 - f_1_11 * r_2 * h6_p4 + fs_900_143143 * r_4 * h4_0 + fs_180_20449 * r_4 * h4_p4 - fs_25_847 * r_6 * h2_0) + fs_14175_64 * e_3 * h2_0;

        pc_63[k] = e_0 * (-fs_135_112 * h4_p1 - fs_2025_14 * r_2 * h2_p1) + e_1 * (-fs_200_121 * h6_p1 + fs_75_44 * h6_p5 + fs_135_847 * r_2 * h4_p1 + fs_50_7 * r_4 * h2_p1) + e_2 * (-fs_189_40898 * h8_p1 - fs_27_286 * h8_p5 + fs_32_1089 * r_2 * h6_p1 - fs_1_33 * r_2 * h6_p5 - fs_135_143143 * r_4 * h4_p1 - fs_200_7623 * r_6 * h2_p1) + fs_1575_8 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p6, ph8_p2, ph8_p3, ph8_p6, ph8_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_64[k] = e_0 * (-fs_135_7 * h4_p2 - fs_18225_112 * r_2 * h2_p2) + e_1 * (-fs_1125_968 * h6_p2 + fs_225_88 * h6_p6 + fs_2160_847 * r_2 * h4_p2 + fs_225_28 * r_4 * h2_p2) + e_2 * (-fs_30_20449 * h8_p2 - fs_14_143 * h8_p6 + fs_5_242 * r_2 * h6_p2 - fs_1_22 * r_2 * h6_p6 - fs_2160_143143 * r_4 * h4_p2 - fs_25_847 * r_6 * h2_p2) + fs_14175_64 * e_3 * h2_p2;

        pc_65[k] = -fs_675_16 * e_0 * h4_p3 + e_1 * (-f_15_22 * h6_p3 + fs_675_121 * r_2 * h4_p3) + e_2 * (-fs_1_3718 * h8_p3 - fs_21_286 * h8_p7 + f_1_11 * r_2 * h6_p3 - fs_675_20449 * r_4 * h4_p3);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_66[k] = e_0 * (fs_2025_112 * h4_m2 - fs_30375_112 * r_2 * h2_m2) + e_1 * (fs_75_968 * h6_m2 - fs_2025_847 * r_2 * h4_m2 + fs_375_28 * r_4 * h2_m2) + e_2 * (-fs_28_143 * h8_m8 + fs_1_40898 * h8_m2 - fs_1_726 * r_2 * h6_m2 + fs_2025_143143 * r_4 * h4_m2 - fs_125_2541 * r_6 * h2_m2) + fs_23625_64 * e_3 * h2_m2;

        pc_67[k] = e_0 * (fs_3645_112 * h4_m1 - fs_6075_56 * r_2 * h2_m1) + e_1 * (fs_75_242 * h6_m1 - fs_3645_847 * r_2 * h4_m1 + fs_75_14 * r_4 * h2_m1) + e_2 * (-fs_35_286 * h8_m7 + fs_7_40898 * h8_m1 - fs_2_363 * r_2 * h6_m1 + fs_3645_143143 * r_4 * h4_m1 - fs_50_2541 * r_6 * h2_m1) + fs_4725_32 * e_3 * h2_m1;

        pc_68[k] = -fs_75_88 * e_1 * h6_m6 + e_2 * (-fs_21_286 * h8_m6 + fs_1_66 * r_2 * h6_m6);
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_69[k] = e_0 * (fs_3375_112 * h4_m1 - fs_2025_224 * r_2 * h2_m1) + e_1 * (-fs_75_44 * h6_m5 + fs_1225_968 * h6_m1 - fs_3375_847 * r_2 * h4_m1 + fs_25_56 * r_4 * h2_m1) + e_2 * (-fs_6_143 * h8_m5 + fs_42_20449 * h8_m1 + fs_1_33 * r_2 * h6_m5 - fs_49_2178 * r_2 * h6_m1 + fs_3375_143143 * r_4 * h4_m1 - fs_25_15246 * r_6 * h2_m1) + fs_1575_128 * e_3 * h2_m1;

        pc_70[k] = e_0 * (-fs_135_56 * h4_m4 - fs_30375_1568 * h4_m2 + fs_2025_1568 * r_2 * h2_m2) + e_1 * (-fs_525_242 * h6_m4 - fs_875_484 * h6_m2 + fs_270_847 * r_2 * h4_m4 + fs_30375_11858 * r_2 * h4_m2 - fs_25_392 * r_4 * h2_m2) + e_2 * (-fs_42_1859 * h8_m4 - fs_105_20449 * h8_m2 + fs_14_363 * r_2 * h6_m4 + fs_35_1089 * r_2 * h6_m2 - fs_270_143143 * r_4 * h4_m4 - fs_30375_2004002 * r_4 * h4_m2 + fs_25_106722 * r_6 * h2_m2) - fs_225_128 * e_3 * h2_m2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph8_p2, ph8_p3, ph8_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_71[k] = -fs_2025_112 * e_0 * h4_p3 + e_1 * (-fs_525_121 * h6_p3 + fs_2025_847 * r_2 * h4_p3) + e_2 * (-fs_42_1859 * h8_p3 + fs_28_363 * r_2 * h6_p3 - fs_2025_143143 * r_4 * h4_p3);

        pc_72[k] = e_0 * (fs_30375_1568 * h4_p2 - fs_135_56 * h4_p4 - fs_2025_1568 * r_2 * h2_p2) + e_1 * (fs_875_484 * h6_p2 - fs_525_242 * h6_p4 - fs_30375_11858 * r_2 * h4_p2 + fs_270_847 * r_2 * h4_p4 + fs_25_392 * r_4 * h2_p2) + e_2 * (fs_105_20449 * h8_p2 - fs_42_1859 * h8_p4 - fs_35_1089 * r_2 * h6_p2 + fs_14_363 * r_2 * h6_p4 + fs_30375_2004002 * r_4 * h4_p2 - fs_270_143143 * r_4 * h4_p4 - fs_25_106722 * r_6 * h2_p2) + fs_225_128 * e_3 * h2_p2;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ph8_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_73[k] = e_0 * (-fs_3375_112 * h4_p1 + fs_2025_224 * r_2 * h2_p1) + e_1 * (-fs_1225_968 * h6_p1 - fs_75_44 * h6_p5 + fs_3375_847 * r_2 * h4_p1 - fs_25_56 * r_4 * h2_p1) + e_2 * (-fs_42_20449 * h8_p1 - fs_6_143 * h8_p5 + fs_49_2178 * r_2 * h6_p1 + fs_1_33 * r_2 * h6_p5 - fs_3375_143143 * r_4 * h4_p1 + fs_25_15246 * r_6 * h2_p1) - fs_1575_128 * e_3 * h2_p1;

        pc_74[k] = e_0 * (fs_2025_28 * h4_0 - fs_2025_28 * r_2 * h2_0) + e_1 * (fs_175_121 * h6_0 - fs_75_88 * h6_p6 - fs_8100_847 * r_2 * h4_0 + fs_25_7 * r_4 * h2_0) + e_2 * (fs_28_20449 * h8_0 - fs_21_286 * h8_p6 - fs_28_1089 * r_2 * h6_0 + fs_1_66 * r_2 * h6_p6 + fs_8100_143143 * r_4 * h4_0 - fs_100_7623 * r_6 * h2_0) + fs_1575_16 * e_3 * h2_0;

        pc_75[k] = e_0 * (fs_3645_112 * h4_p1 - fs_6075_56 * r_2 * h2_p1) + e_1 * (fs_75_242 * h6_p1 - fs_3645_847 * r_2 * h4_p1 + fs_75_14 * r_4 * h2_p1) + e_2 * (fs_7_40898 * h8_p1 - fs_35_286 * h8_p7 - fs_2_363 * r_2 * h6_p1 + fs_3645_143143 * r_4 * h4_p1 - fs_50_2541 * r_6 * h2_p1) + fs_4725_32 * e_3 * h2_p1;
    }

    // NOTE: the rows are formed in 33 loops, as the vectorizer runs out of
    // registers with all 77 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];

        pc_76[k] = e_0 * (fs_2025_112 * h4_p2 - fs_30375_112 * r_2 * h2_p2) + e_1 * (fs_75_968 * h6_p2 - fs_2025_847 * r_2 * h4_p2 + fs_375_28 * r_4 * h2_p2) + e_2 * (fs_1_40898 * h8_p2 - fs_28_143 * h8_p8 - fs_1_726 * r_2 * h6_p2 + fs_2025_143143 * r_4 * h4_p2 - fs_125_2541 * r_6 * h2_p2) + fs_23625_64 * e_3 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[77] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76};

    for (size_t m = 0; m < 77; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
