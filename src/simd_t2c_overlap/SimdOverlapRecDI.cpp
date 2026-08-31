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



#include "SimdOverlapRecDI.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_di_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDI.compute_di_overlap: Basis functions must be of angular momenta two and six"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDI.compute_di_overlap: Harmonics must reach angular momentum eight"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDI.compute_di_overlap: Number of values exceeds number of atom pairs"));
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

    // NOTE: the buffer spans the atom pairs reached by the pair of primitives
    // reaching furthest, which is searched for rather than assumed. The
    // primitives are sorted by descending exponent, but the bound of a pair of
    // primitives carries their prefactor as well as their decay, so a tighter
    // pair with a larger prefactor reaches further than a more diffuse pair with
    // a smaller one, and the last pair is not always the furthest reaching.

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 65 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = CSimdMatrix(3, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);

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

            const auto f_0 = fbase * aexp * aexp * aexp * aexp * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * aexp * aexp * aexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * aexp * aexp * aexp * aexp / fexp / fexp / fexp / fexp / fexp / fexp;

            // NOTE: the exponential depends on the pair of primitives alone, so it is
            // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fss = std::exp(-fmu * ab_2[k]);

                pe_0[k] += f_0 * fss;
                pe_1[k] += f_1 * fss;
                pe_2[k] += f_2 * fss;
            }
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_9_22 = std::sqrt(9.0 / 22.0);
    const auto fs_405_22 = std::sqrt(405.0 / 22.0);
    const auto fs_1_8450 = std::sqrt(1.0 / 8450.0);
    const auto fs_14_65 = std::sqrt(14.0 / 65.0);
    const auto fs_2_275 = std::sqrt(2.0 / 275.0);
    const auto fs_405_3718 = std::sqrt(405.0 / 3718.0);
    const auto fs_4455_32 = std::sqrt(139.21875);
    const auto fs_45_44 = std::sqrt(45.0 / 44.0);
    const auto fs_135_11 = std::sqrt(135.0 / 11.0);
    const auto fs_1_1690 = std::sqrt(1.0 / 1690.0);
    const auto fs_21_130 = std::sqrt(21.0 / 130.0);
    const auto fs_1_55 = std::sqrt(1.0 / 55.0);
    const auto fs_135_1859 = std::sqrt(135.0 / 1859.0);
    const auto fs_1485_16 = std::sqrt(92.8125);
    const auto fs_405_242 = std::sqrt(405.0 / 242.0);
    const auto fs_945_121 = std::sqrt(945.0 / 121.0);
    const auto fs_3_1690 = std::sqrt(3.0 / 1690.0);
    const auto fs_77_650 = std::sqrt(77.0 / 650.0);
    const auto fs_18_605 = std::sqrt(18.0 / 605.0);
    const auto fs_945_20449 = std::sqrt(945.0 / 20449.0);
    const auto fs_945_16 = std::sqrt(59.0625);
    const auto fs_270_121 = std::sqrt(270.0 / 121.0);
    const auto fs_567_121 = std::sqrt(567.0 / 121.0);
    const auto fs_7_1690 = std::sqrt(7.0 / 1690.0);
    const auto fs_11_130 = std::sqrt(11.0 / 130.0);
    const auto fs_24_605 = std::sqrt(24.0 / 605.0);
    const auto fs_567_20449 = std::sqrt(567.0 / 20449.0);
    const auto fs_567_16 = std::sqrt(35.4375);
    const auto fs_630_121 = std::sqrt(630.0 / 121.0);
    const auto fs_9_242 = std::sqrt(9.0 / 242.0);
    const auto fs_14_845 = std::sqrt(14.0 / 845.0);
    const auto fs_99_1690 = std::sqrt(99.0 / 1690.0);
    const auto fs_56_605 = std::sqrt(56.0 / 605.0);
    const auto fs_630_20449 = std::sqrt(630.0 / 20449.0);
    const auto fs_9_40898 = std::sqrt(9.0 / 40898.0);
    const auto fs_315_8 = std::sqrt(39.375);
    const auto fs_9_32 = std::sqrt(0.28125);
    const auto fs_1323_484 = std::sqrt(1323.0 / 484.0);
    const auto fs_315_242 = std::sqrt(315.0 / 242.0);
    const auto fs_45_242 = std::sqrt(45.0 / 242.0);
    const auto fs_63_4225 = std::sqrt(63.0 / 4225.0);
    const auto fs_33_845 = std::sqrt(33.0 / 845.0);
    const auto fs_147_3025 = std::sqrt(147.0 / 3025.0);
    const auto fs_315_40898 = std::sqrt(315.0 / 40898.0);
    const auto fs_45_40898 = std::sqrt(45.0 / 40898.0);
    const auto fs_315_32 = std::sqrt(9.84375);
    const auto fs_45_32 = std::sqrt(1.40625);
    const auto fs_135_121 = std::sqrt(135.0 / 121.0);
    const auto fs_42_845 = std::sqrt(42.0 / 845.0);
    const auto fs_135_20449 = std::sqrt(135.0 / 20449.0);
    const auto fs_135_16 = std::sqrt(8.4375);
    const auto f_3_2 = 1.5;
    const auto fs_1_650 = std::sqrt(1.0 / 650.0);
    const auto fs_7_130 = std::sqrt(7.0 / 130.0);
    const auto f_1_5 = 1.0 / 5.0;
    const auto fs_243_88 = std::sqrt(243.0 / 88.0);
    const auto fs_135_22 = std::sqrt(135.0 / 22.0);
    const auto fs_24_4225 = std::sqrt(24.0 / 4225.0);
    const auto fs_28_325 = std::sqrt(28.0 / 325.0);
    const auto fs_27_550 = std::sqrt(27.0 / 550.0);
    const auto fs_135_3718 = std::sqrt(135.0 / 3718.0);
    const auto fs_1485_32 = std::sqrt(46.40625);
    const auto fs_2205_968 = std::sqrt(2205.0 / 968.0);
    const auto fs_1080_121 = std::sqrt(1080.0 / 121.0);
    const auto fs_11_845 = std::sqrt(11.0 / 845.0);
    const auto fs_33_325 = std::sqrt(33.0 / 325.0);
    const auto fs_49_1210 = std::sqrt(49.0 / 1210.0);
    const auto fs_1080_20449 = std::sqrt(1080.0 / 20449.0);
    const auto fs_135_2 = std::sqrt(67.5);
    const auto fs_675_484 = std::sqrt(675.0 / 484.0);
    const auto fs_1134_121 = std::sqrt(1134.0 / 121.0);
    const auto fs_81_242 = std::sqrt(81.0 / 242.0);
    const auto f_2_13 = 2.0 / 13.0;
    const auto fs_88_845 = std::sqrt(88.0 / 845.0);
    const auto fs_3_121 = std::sqrt(3.0 / 121.0);
    const auto fs_1134_20449 = std::sqrt(1134.0 / 20449.0);
    const auto fs_81_40898 = std::sqrt(81.0 / 40898.0);
    const auto fs_567_8 = std::sqrt(70.875);
    const auto fs_81_32 = std::sqrt(2.53125);
    const auto fs_135_242 = std::sqrt(135.0 / 242.0);
    const auto fs_1008_121 = std::sqrt(1008.0 / 121.0);
    const auto f_12_11 = 12.0 / 11.0;
    const auto fs_63_1690 = std::sqrt(63.0 / 1690.0);
    const auto fs_33_338 = std::sqrt(33.0 / 338.0);
    const auto fs_6_605 = std::sqrt(6.0 / 605.0);
    const auto fs_1008_20449 = std::sqrt(1008.0 / 20449.0);
    const auto f_12_143 = 12.0 / 143.0;
    const auto fs_63 = std::sqrt(63.0);
    const auto f_3 = 3.0;
    const auto fs_63_484 = std::sqrt(63.0 / 484.0);
    const auto fs_1575_121 = std::sqrt(1575.0 / 121.0);
    const auto fs_315_121 = std::sqrt(315.0 / 121.0);
    const auto fs_448_4225 = std::sqrt(448.0 / 4225.0);
    const auto fs_72_845 = std::sqrt(72.0 / 845.0);
    const auto fs_7_3025 = std::sqrt(7.0 / 3025.0);
    const auto fs_1575_20449 = std::sqrt(1575.0 / 20449.0);
    const auto fs_315_20449 = std::sqrt(315.0 / 20449.0);
    const auto fs_1575_16 = std::sqrt(98.4375);
    const auto fs_315_16 = std::sqrt(19.6875);
    const auto fs_588_4225 = std::sqrt(588.0 / 4225.0);
    const auto fs_7_325 = std::sqrt(7.0 / 325.0);
    const auto f_2_5 = 2.0 / 5.0;
    const auto fs_18_325 = std::sqrt(18.0 / 325.0);
    const auto f_3_11 = 3.0 / 11.0;
    const auto fs_405_121 = std::sqrt(405.0 / 121.0);
    const auto fs_396_4225 = std::sqrt(396.0 / 4225.0);
    const auto f_2_55 = 2.0 / 55.0;
    const auto fs_405_20449 = std::sqrt(405.0 / 20449.0);
    const auto fs_405_16 = std::sqrt(25.3125);
    const auto f_15_22 = 15.0 / 22.0;
    const auto fs_972_121 = std::sqrt(972.0 / 121.0);
    const auto fs_22_169 = std::sqrt(22.0 / 169.0);
    const auto f_1_11 = 1.0 / 11.0;
    const auto fs_972_20449 = std::sqrt(972.0 / 20449.0);
    const auto fs_243_4 = std::sqrt(60.75);
    const auto f_15_11 = 15.0 / 11.0;
    const auto fs_1512_121 = std::sqrt(1512.0 / 121.0);
    const auto fs_27_169 = std::sqrt(27.0 / 169.0);
    const auto f_2_11 = 2.0 / 11.0;
    const auto fs_1512_20449 = std::sqrt(1512.0 / 20449.0);
    const auto fs_189_2 = std::sqrt(94.5);
    const auto f_39_22 = 39.0 / 22.0;
    const auto fs_1890_121 = std::sqrt(1890.0 / 121.0);
    const auto fs_756_4225 = std::sqrt(756.0 / 4225.0);
    const auto f_13_55 = 13.0 / 55.0;
    const auto fs_1890_20449 = std::sqrt(1890.0 / 20449.0);
    const auto fs_945_8 = std::sqrt(118.125);
    const auto f_21_11 = 21.0 / 11.0;
    const auto f_45_11 = 45.0 / 11.0;
    const auto f_28_65 = 28.0 / 65.0;
    const auto f_14_55 = 14.0 / 55.0;
    const auto f_45_143 = 45.0 / 143.0;
    const auto f_45_4 = 11.25;

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ph8_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * (fs_9_22 * h6_p4 - fs_405_22 * r_2 * h4_p4) + e_1 * (fs_1_8450 * h8_p4 - fs_14_65 * h8_p8 - fs_2_275 * r_2 * h6_p4 + fs_405_3718 * r_4 * h4_p4) + fs_4455_32 * e_2 * h4_p4;

        pc_1[k] = e_0 * (fs_45_44 * h6_p3 - fs_135_11 * r_2 * h4_p3) + e_1 * (fs_1_1690 * h8_p3 - fs_21_130 * h8_p7 - fs_1_55 * r_2 * h6_p3 + fs_135_1859 * r_4 * h4_p3) + fs_1485_16 * e_2 * h4_p3;

        pc_2[k] = e_0 * (fs_405_242 * h6_p2 - fs_9_22 * h6_p6 - fs_945_121 * r_2 * h4_p2) + e_1 * (fs_3_1690 * h8_p2 - fs_77_650 * h8_p6 - fs_18_605 * r_2 * h6_p2 + fs_2_275 * r_2 * h6_p6 + fs_945_20449 * r_4 * h4_p2) + fs_945_16 * e_2 * h4_p2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p3, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p3, ph8_p4, ph8_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_3[k] = e_0 * (fs_270_121 * h6_p1 - fs_45_44 * h6_p5 - fs_567_121 * r_2 * h4_p1) + e_1 * (fs_7_1690 * h8_p1 - fs_11_130 * h8_p5 - fs_24_605 * r_2 * h6_p1 + fs_1_55 * r_2 * h6_p5 + fs_567_20449 * r_4 * h4_p1) + fs_567_16 * e_2 * h4_p1;

        pc_4[k] = e_0 * (fs_630_121 * h6_0 - fs_405_242 * h6_p4 - fs_630_121 * r_2 * h4_0 + fs_9_242 * r_2 * h4_p4) + e_1 * (fs_14_845 * h8_0 - fs_99_1690 * h8_p4 - fs_56_605 * r_2 * h6_0 + fs_18_605 * r_2 * h6_p4 + fs_630_20449 * r_4 * h4_0 - fs_9_40898 * r_4 * h4_p4) + e_2 * (fs_315_8 * h4_0 - fs_9_32 * h4_p4);

        pc_5[k] = e_0 * (-fs_1323_484 * h6_p1 - fs_270_121 * h6_p3 + fs_315_242 * r_2 * h4_p1 + fs_45_242 * r_2 * h4_p3) + e_1 * (-fs_63_4225 * h8_p1 - fs_33_845 * h8_p3 + fs_147_3025 * r_2 * h6_p1 + fs_24_605 * r_2 * h6_p3 - fs_315_40898 * r_4 * h4_p1 - fs_45_40898 * r_4 * h4_p3) + e_2 * (-fs_315_32 * h4_p1 - fs_45_32 * h4_p3);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_6[k] = e_0 * (fs_630_121 * h6_m2 - fs_135_121 * r_2 * h4_m2) + e_1 * (fs_42_845 * h8_m2 - fs_56_605 * r_2 * h6_m2 + fs_135_20449 * r_4 * h4_m2) + fs_135_16 * e_2 * h4_m2;

        pc_7[k] = e_0 * (fs_270_121 * h6_m3 - fs_1323_484 * h6_m1 - fs_45_242 * r_2 * h4_m3 + fs_315_242 * r_2 * h4_m1) + e_1 * (fs_33_845 * h8_m3 - fs_63_4225 * h8_m1 - fs_24_605 * r_2 * h6_m3 + fs_147_3025 * r_2 * h6_m1 + fs_45_40898 * r_4 * h4_m3 - fs_315_40898 * r_4 * h4_m1) + e_2 * (fs_45_32 * h4_m3 - fs_315_32 * h4_m1);

        pc_8[k] = e_0 * (fs_405_242 * h6_m4 - fs_9_242 * r_2 * h4_m4) + e_1 * (fs_99_1690 * h8_m4 - fs_18_605 * r_2 * h6_m4 + fs_9_40898 * r_4 * h4_m4) + fs_9_32 * e_2 * h4_m4;

        pc_9[k] = e_0 * (fs_45_44 * h6_m5 - fs_270_121 * h6_m1 + fs_567_121 * r_2 * h4_m1) + e_1 * (fs_11_130 * h8_m5 - fs_7_1690 * h8_m1 - fs_1_55 * r_2 * h6_m5 + fs_24_605 * r_2 * h6_m1 - fs_567_20449 * r_4 * h4_m1) - fs_567_16 * e_2 * h4_m1;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m4, ph6_m3, ph6_m2, ph8_m8, ph8_m7, ph8_m6, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_10[k] = e_0 * (fs_9_22 * h6_m6 - fs_405_242 * h6_m2 + fs_945_121 * r_2 * h4_m2) + e_1 * (fs_77_650 * h8_m6 - fs_3_1690 * h8_m2 - fs_2_275 * r_2 * h6_m6 + fs_18_605 * r_2 * h6_m2 - fs_945_20449 * r_4 * h4_m2) - fs_945_16 * e_2 * h4_m2;

        pc_11[k] = e_0 * (-fs_45_44 * h6_m3 + fs_135_11 * r_2 * h4_m3) + e_1 * (fs_21_130 * h8_m7 - fs_1_1690 * h8_m3 + fs_1_55 * r_2 * h6_m3 - fs_135_1859 * r_4 * h4_m3) - fs_1485_16 * e_2 * h4_m3;

        pc_12[k] = e_0 * (-fs_9_22 * h6_m4 + fs_405_22 * r_2 * h4_m4) + e_1 * (fs_14_65 * h8_m8 - fs_1_8450 * h8_m4 + fs_2_275 * r_2 * h6_m4 - fs_405_3718 * r_4 * h4_m4) - fs_4455_32 * e_2 * h4_m4;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ph8_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_13[k] = -f_3_2 * e_0 * h6_p5 + e_1 * (-fs_1_650 * h8_p5 - fs_7_130 * h8_p7 + f_1_5 * r_2 * h6_p5);

        pc_14[k] = e_0 * (-fs_243_88 * h6_p4 + f_3_2 * h6_p6 - fs_135_22 * r_2 * h4_p4) + e_1 * (-fs_24_4225 * h8_p4 - fs_28_325 * h8_p6 + fs_27_550 * r_2 * h6_p4 - f_1_5 * r_2 * h6_p6 + fs_135_3718 * r_4 * h4_p4) + fs_1485_32 * e_2 * h4_p4;

        pc_15[k] = e_0 * (-fs_2205_968 * h6_p3 + fs_243_88 * h6_p5 - fs_1080_121 * r_2 * h4_p3) + e_1 * (-fs_11_845 * h8_p3 - fs_33_325 * h8_p5 + fs_49_1210 * r_2 * h6_p3 - fs_27_550 * r_2 * h6_p5 + fs_1080_20449 * r_4 * h4_p3) + fs_135_2 * e_2 * h4_p3;

        pc_16[k] = e_0 * (-fs_675_484 * h6_p2 + fs_2205_968 * h6_p4 - fs_1134_121 * r_2 * h4_p2 - fs_81_242 * r_2 * h4_p4) + e_1 * (-f_2_13 * h8_p2 - fs_88_845 * h8_p4 + fs_3_121 * r_2 * h6_p2 - fs_49_1210 * r_2 * h6_p4 + fs_1134_20449 * r_4 * h4_p2 + fs_81_40898 * r_4 * h4_p4) + e_2 * (fs_567_8 * h4_p2 + fs_81_32 * h4_p4);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph8_0, ph8_p1, ph8_p2, ph8_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_17[k] = e_0 * (-fs_135_242 * h6_p1 + fs_675_484 * h6_p3 - fs_1008_121 * r_2 * h4_p1 - f_12_11 * r_2 * h4_p3) + e_1 * (-fs_63_1690 * h8_p1 - fs_33_338 * h8_p3 + fs_6_605 * r_2 * h6_p1 - fs_3_121 * r_2 * h6_p3 + fs_1008_20449 * r_4 * h4_p1 + f_12_143 * r_4 * h4_p3) + e_2 * (fs_63 * h4_p1 + f_3 * h4_p3);

        pc_18[k] = e_0 * (-fs_63_484 * h6_0 + fs_135_242 * h6_p2 - fs_1575_121 * r_2 * h4_0 - fs_315_121 * r_2 * h4_p2) + e_1 * (-fs_448_4225 * h8_0 - fs_72_845 * h8_p2 + fs_7_3025 * r_2 * h6_0 - fs_6_605 * r_2 * h6_p2 + fs_1575_20449 * r_4 * h4_0 + fs_315_20449 * r_4 * h4_p2) + e_2 * (fs_1575_16 * h4_0 + fs_315_16 * h4_p2);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ph8_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_19[k] = e_0 * (-fs_63_484 * h6_m1 + fs_1080_121 * r_2 * h4_m1) + e_1 * (fs_588_4225 * h8_m1 + fs_7_3025 * r_2 * h6_m1 - fs_1080_20449 * r_4 * h4_m1) - fs_135_2 * e_2 * h4_m1;

        pc_20[k] = e_0 * (-fs_135_242 * h6_m2 + fs_315_121 * r_2 * h4_m2) + e_1 * (fs_72_845 * h8_m2 + fs_6_605 * r_2 * h6_m2 - fs_315_20449 * r_4 * h4_m2) - fs_315_16 * e_2 * h4_m2;

        pc_21[k] = e_0 * (-fs_675_484 * h6_m3 + fs_135_242 * h6_m1 + f_12_11 * r_2 * h4_m3 + fs_1008_121 * r_2 * h4_m1) + e_1 * (fs_33_338 * h8_m3 + fs_63_1690 * h8_m1 + fs_3_121 * r_2 * h6_m3 - fs_6_605 * r_2 * h6_m1 - f_12_143 * r_4 * h4_m3 - fs_1008_20449 * r_4 * h4_m1) + e_2 * (-f_3 * h4_m3 - fs_63 * h4_m1);

        pc_22[k] = e_0 * (-fs_2205_968 * h6_m4 + fs_675_484 * h6_m2 + fs_81_242 * r_2 * h4_m4 + fs_1134_121 * r_2 * h4_m2) + e_1 * (fs_88_845 * h8_m4 + f_2_13 * h8_m2 + fs_49_1210 * r_2 * h6_m4 - fs_3_121 * r_2 * h6_m2 - fs_81_40898 * r_4 * h4_m4 - fs_1134_20449 * r_4 * h4_m2) + e_2 * (-fs_81_32 * h4_m4 - fs_567_8 * h4_m2);

        pc_23[k] = e_0 * (-fs_243_88 * h6_m5 + fs_2205_968 * h6_m3 + fs_1080_121 * r_2 * h4_m3) + e_1 * (fs_33_325 * h8_m5 + fs_11_845 * h8_m3 + fs_27_550 * r_2 * h6_m5 - fs_49_1210 * r_2 * h6_m3 - fs_1080_20449 * r_4 * h4_m3) - fs_135_2 * e_2 * h4_m3;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m7, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_24[k] = e_0 * (-f_3_2 * h6_m6 + fs_243_88 * h6_m4 + fs_135_22 * r_2 * h4_m4) + e_1 * (fs_28_325 * h8_m6 + fs_24_4225 * h8_m4 + f_1_5 * r_2 * h6_m6 - fs_27_550 * r_2 * h6_m4 - fs_135_3718 * r_4 * h4_m4) - fs_1485_32 * e_2 * h4_m4;

        pc_25[k] = f_3_2 * e_0 * h6_m5 + e_1 * (fs_7_130 * h8_m7 + fs_1_650 * h8_m5 - f_1_5 * r_2 * h6_m5);

        pc_26[k] = f_3 * e_0 * h6_m6 + e_1 * (fs_7_325 * h8_m6 - f_2_5 * r_2 * h6_m6);

        pc_27[k] = f_3_2 * e_0 * h6_m5 + e_1 * (fs_18_325 * h8_m5 - f_1_5 * r_2 * h6_m5);

        pc_28[k] = e_0 * (f_3_11 * h6_m4 - fs_405_121 * r_2 * h4_m4) + e_1 * (fs_396_4225 * h8_m4 - f_2_55 * r_2 * h6_m4 + fs_405_20449 * r_4 * h4_m4) + fs_405_16 * e_2 * h4_m4;

        pc_29[k] = e_0 * (-f_15_22 * h6_m3 - fs_972_121 * r_2 * h4_m3) + e_1 * (fs_22_169 * h8_m3 + f_1_11 * r_2 * h6_m3 + fs_972_20449 * r_4 * h4_m3) + fs_243_4 * e_2 * h4_m3;

        pc_30[k] = e_0 * (-f_15_11 * h6_m2 - fs_1512_121 * r_2 * h4_m2) + e_1 * (fs_27_169 * h8_m2 + f_2_11 * r_2 * h6_m2 + fs_1512_20449 * r_4 * h4_m2) + fs_189_2 * e_2 * h4_m2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m1, ph4_0, ph4_p1, ph4_p2, ph6_m1, ph6_0, ph6_p1, ph6_p2, ph8_m1, ph8_0, ph8_p1, ph8_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_31[k] = e_0 * (-f_39_22 * h6_m1 - fs_1890_121 * r_2 * h4_m1) + e_1 * (fs_756_4225 * h8_m1 + f_13_55 * r_2 * h6_m1 + fs_1890_20449 * r_4 * h4_m1) + fs_945_8 * e_2 * h4_m1;

        pc_32[k] = e_0 * (-f_21_11 * h6_0 - f_45_11 * r_2 * h4_0) + e_1 * (f_28_65 * h8_0 + f_14_55 * r_2 * h6_0 + f_45_143 * r_4 * h4_0) + f_45_4 * e_2 * h4_0;

        pc_33[k] = e_0 * (-f_39_22 * h6_p1 - fs_1890_121 * r_2 * h4_p1) + e_1 * (fs_756_4225 * h8_p1 + f_13_55 * r_2 * h6_p1 + fs_1890_20449 * r_4 * h4_p1) + fs_945_8 * e_2 * h4_p1;

        pc_34[k] = e_0 * (-f_15_11 * h6_p2 - fs_1512_121 * r_2 * h4_p2) + e_1 * (fs_27_169 * h8_p2 + f_2_11 * r_2 * h6_p2 + fs_1512_20449 * r_4 * h4_p2) + fs_189_2 * e_2 * h4_p2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p3, ph4_p4, ph6_m5, ph6_p3, ph6_p4, ph6_p5, ph6_p6, ph8_m7, ph8_m5, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_35[k] = e_0 * (-f_15_22 * h6_p3 - fs_972_121 * r_2 * h4_p3) + e_1 * (fs_22_169 * h8_p3 + f_1_11 * r_2 * h6_p3 + fs_972_20449 * r_4 * h4_p3) + fs_243_4 * e_2 * h4_p3;

        pc_36[k] = e_0 * (f_3_11 * h6_p4 - fs_405_121 * r_2 * h4_p4) + e_1 * (fs_396_4225 * h8_p4 - f_2_55 * r_2 * h6_p4 + fs_405_20449 * r_4 * h4_p4) + fs_405_16 * e_2 * h4_p4;

        pc_37[k] = f_3_2 * e_0 * h6_p5 + e_1 * (fs_18_325 * h8_p5 - f_1_5 * r_2 * h6_p5);

        pc_38[k] = f_3 * e_0 * h6_p6 + e_1 * (fs_7_325 * h8_p6 - f_2_5 * r_2 * h6_p6);

        pc_39[k] = -f_3_2 * e_0 * h6_m5 + e_1 * (fs_7_130 * h8_m7 - fs_1_650 * h8_m5 + f_1_5 * r_2 * h6_m5);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_40[k] = e_0 * (-f_3_2 * h6_m6 - fs_243_88 * h6_m4 - fs_135_22 * r_2 * h4_m4) + e_1 * (fs_28_325 * h8_m6 - fs_24_4225 * h8_m4 + f_1_5 * r_2 * h6_m6 + fs_27_550 * r_2 * h6_m4 + fs_135_3718 * r_4 * h4_m4) + fs_1485_32 * e_2 * h4_m4;

        pc_41[k] = e_0 * (-fs_243_88 * h6_m5 - fs_2205_968 * h6_m3 - fs_1080_121 * r_2 * h4_m3) + e_1 * (fs_33_325 * h8_m5 - fs_11_845 * h8_m3 + fs_27_550 * r_2 * h6_m5 + fs_49_1210 * r_2 * h6_m3 + fs_1080_20449 * r_4 * h4_m3) + fs_135_2 * e_2 * h4_m3;

        pc_42[k] = e_0 * (-fs_2205_968 * h6_m4 - fs_675_484 * h6_m2 + fs_81_242 * r_2 * h4_m4 - fs_1134_121 * r_2 * h4_m2) + e_1 * (fs_88_845 * h8_m4 - f_2_13 * h8_m2 + fs_49_1210 * r_2 * h6_m4 + fs_3_121 * r_2 * h6_m2 - fs_81_40898 * r_4 * h4_m4 + fs_1134_20449 * r_4 * h4_m2) + e_2 * (-fs_81_32 * h4_m4 + fs_567_8 * h4_m2);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m3, ph4_m2, ph4_m1, ph4_p1, ph6_m3, ph6_m2, ph6_m1, ph6_p1, ph8_m3, ph8_m2, ph8_m1, ph8_p1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_p1 = ph8_p1[k];

        pc_43[k] = e_0 * (-fs_675_484 * h6_m3 - fs_135_242 * h6_m1 + f_12_11 * r_2 * h4_m3 - fs_1008_121 * r_2 * h4_m1) + e_1 * (fs_33_338 * h8_m3 - fs_63_1690 * h8_m1 + fs_3_121 * r_2 * h6_m3 + fs_6_605 * r_2 * h6_m1 - f_12_143 * r_4 * h4_m3 + fs_1008_20449 * r_4 * h4_m1) + e_2 * (-f_3 * h4_m3 + fs_63 * h4_m1);

        pc_44[k] = e_0 * (-fs_135_242 * h6_m2 + fs_315_121 * r_2 * h4_m2) + e_1 * (fs_72_845 * h8_m2 + fs_6_605 * r_2 * h6_m2 - fs_315_20449 * r_4 * h4_m2) - fs_315_16 * e_2 * h4_m2;

        pc_45[k] = e_0 * (-fs_63_484 * h6_p1 + fs_1080_121 * r_2 * h4_p1) + e_1 * (fs_588_4225 * h8_p1 + fs_7_3025 * r_2 * h6_p1 - fs_1080_20449 * r_4 * h4_p1) - fs_135_2 * e_2 * h4_p1;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph8_0, ph8_p1, ph8_p2, ph8_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_46[k] = e_0 * (-fs_63_484 * h6_0 - fs_135_242 * h6_p2 - fs_1575_121 * r_2 * h4_0 + fs_315_121 * r_2 * h4_p2) + e_1 * (-fs_448_4225 * h8_0 + fs_72_845 * h8_p2 + fs_7_3025 * r_2 * h6_0 + fs_6_605 * r_2 * h6_p2 + fs_1575_20449 * r_4 * h4_0 - fs_315_20449 * r_4 * h4_p2) + e_2 * (fs_1575_16 * h4_0 - fs_315_16 * h4_p2);

        pc_47[k] = e_0 * (-fs_135_242 * h6_p1 - fs_675_484 * h6_p3 - fs_1008_121 * r_2 * h4_p1 + f_12_11 * r_2 * h4_p3) + e_1 * (-fs_63_1690 * h8_p1 + fs_33_338 * h8_p3 + fs_6_605 * r_2 * h6_p1 + fs_3_121 * r_2 * h6_p3 + fs_1008_20449 * r_4 * h4_p1 - f_12_143 * r_4 * h4_p3) + e_2 * (fs_63 * h4_p1 - f_3 * h4_p3);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ph8_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_48[k] = e_0 * (-fs_675_484 * h6_p2 - fs_2205_968 * h6_p4 - fs_1134_121 * r_2 * h4_p2 + fs_81_242 * r_2 * h4_p4) + e_1 * (-f_2_13 * h8_p2 + fs_88_845 * h8_p4 + fs_3_121 * r_2 * h6_p2 + fs_49_1210 * r_2 * h6_p4 + fs_1134_20449 * r_4 * h4_p2 - fs_81_40898 * r_4 * h4_p4) + e_2 * (fs_567_8 * h4_p2 - fs_81_32 * h4_p4);

        pc_49[k] = e_0 * (-fs_2205_968 * h6_p3 - fs_243_88 * h6_p5 - fs_1080_121 * r_2 * h4_p3) + e_1 * (-fs_11_845 * h8_p3 + fs_33_325 * h8_p5 + fs_49_1210 * r_2 * h6_p3 + fs_27_550 * r_2 * h6_p5 + fs_1080_20449 * r_4 * h4_p3) + fs_135_2 * e_2 * h4_p3;

        pc_50[k] = e_0 * (-fs_243_88 * h6_p4 - f_3_2 * h6_p6 - fs_135_22 * r_2 * h4_p4) + e_1 * (-fs_24_4225 * h8_p4 + fs_28_325 * h8_p6 + fs_27_550 * r_2 * h6_p4 + f_1_5 * r_2 * h6_p6 + fs_135_3718 * r_4 * h4_p4) + fs_1485_32 * e_2 * h4_p4;

        pc_51[k] = -f_3_2 * e_0 * h6_p5 + e_1 * (-fs_1_650 * h8_p5 + fs_7_130 * h8_p7 + f_1_5 * r_2 * h6_p5);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m4, ph6_m3, ph6_m2, ph8_m8, ph8_m7, ph8_m6, ph8_m4, ph8_m3, ph8_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_52[k] = e_0 * (fs_9_22 * h6_m4 - fs_405_22 * r_2 * h4_m4) + e_1 * (fs_14_65 * h8_m8 + fs_1_8450 * h8_m4 - fs_2_275 * r_2 * h6_m4 + fs_405_3718 * r_4 * h4_m4) + fs_4455_32 * e_2 * h4_m4;

        pc_53[k] = e_0 * (fs_45_44 * h6_m3 - fs_135_11 * r_2 * h4_m3) + e_1 * (fs_21_130 * h8_m7 + fs_1_1690 * h8_m3 - fs_1_55 * r_2 * h6_m3 + fs_135_1859 * r_4 * h4_m3) + fs_1485_16 * e_2 * h4_m3;

        pc_54[k] = e_0 * (fs_9_22 * h6_m6 + fs_405_242 * h6_m2 - fs_945_121 * r_2 * h4_m2) + e_1 * (fs_77_650 * h8_m6 + fs_3_1690 * h8_m2 - fs_2_275 * r_2 * h6_m6 - fs_18_605 * r_2 * h6_m2 + fs_945_20449 * r_4 * h4_m2) + fs_945_16 * e_2 * h4_m2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m1, ph4_p2, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph6_p2, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ph8_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_55[k] = e_0 * (fs_45_44 * h6_m5 + fs_270_121 * h6_m1 - fs_567_121 * r_2 * h4_m1) + e_1 * (fs_11_130 * h8_m5 + fs_7_1690 * h8_m1 - fs_1_55 * r_2 * h6_m5 - fs_24_605 * r_2 * h6_m1 + fs_567_20449 * r_4 * h4_m1) + fs_567_16 * e_2 * h4_m1;

        pc_56[k] = e_0 * (fs_405_242 * h6_m4 - fs_9_242 * r_2 * h4_m4) + e_1 * (fs_99_1690 * h8_m4 - fs_18_605 * r_2 * h6_m4 + fs_9_40898 * r_4 * h4_m4) + fs_9_32 * e_2 * h4_m4;

        pc_57[k] = e_0 * (fs_270_121 * h6_m3 + fs_1323_484 * h6_m1 - fs_45_242 * r_2 * h4_m3 - fs_315_242 * r_2 * h4_m1) + e_1 * (fs_33_845 * h8_m3 + fs_63_4225 * h8_m1 - fs_24_605 * r_2 * h6_m3 - fs_147_3025 * r_2 * h6_m1 + fs_45_40898 * r_4 * h4_m3 + fs_315_40898 * r_4 * h4_m1) + e_2 * (fs_45_32 * h4_m3 + fs_315_32 * h4_m1);

        pc_58[k] = e_0 * (fs_630_121 * h6_p2 - fs_135_121 * r_2 * h4_p2) + e_1 * (fs_42_845 * h8_p2 - fs_56_605 * r_2 * h6_p2 + fs_135_20449 * r_4 * h4_p2) + fs_135_16 * e_2 * h4_p2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p3, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p3, ph8_p4, ph8_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_59[k] = e_0 * (-fs_1323_484 * h6_p1 + fs_270_121 * h6_p3 + fs_315_242 * r_2 * h4_p1 - fs_45_242 * r_2 * h4_p3) + e_1 * (-fs_63_4225 * h8_p1 + fs_33_845 * h8_p3 + fs_147_3025 * r_2 * h6_p1 - fs_24_605 * r_2 * h6_p3 - fs_315_40898 * r_4 * h4_p1 + fs_45_40898 * r_4 * h4_p3) + e_2 * (-fs_315_32 * h4_p1 + fs_45_32 * h4_p3);

        pc_60[k] = e_0 * (fs_630_121 * h6_0 + fs_405_242 * h6_p4 - fs_630_121 * r_2 * h4_0 - fs_9_242 * r_2 * h4_p4) + e_1 * (fs_14_845 * h8_0 + fs_99_1690 * h8_p4 - fs_56_605 * r_2 * h6_0 - fs_18_605 * r_2 * h6_p4 + fs_630_20449 * r_4 * h4_0 + fs_9_40898 * r_4 * h4_p4) + e_2 * (fs_315_8 * h4_0 + fs_9_32 * h4_p4);

        pc_61[k] = e_0 * (fs_270_121 * h6_p1 + fs_45_44 * h6_p5 - fs_567_121 * r_2 * h4_p1) + e_1 * (fs_7_1690 * h8_p1 + fs_11_130 * h8_p5 - fs_24_605 * r_2 * h6_p1 - fs_1_55 * r_2 * h6_p5 + fs_567_20449 * r_4 * h4_p1) + fs_567_16 * e_2 * h4_p1;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 65 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ph8_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_62[k] = e_0 * (fs_405_242 * h6_p2 + fs_9_22 * h6_p6 - fs_945_121 * r_2 * h4_p2) + e_1 * (fs_3_1690 * h8_p2 + fs_77_650 * h8_p6 - fs_18_605 * r_2 * h6_p2 - fs_2_275 * r_2 * h6_p6 + fs_945_20449 * r_4 * h4_p2) + fs_945_16 * e_2 * h4_p2;

        pc_63[k] = e_0 * (fs_45_44 * h6_p3 - fs_135_11 * r_2 * h4_p3) + e_1 * (fs_1_1690 * h8_p3 + fs_21_130 * h8_p7 - fs_1_55 * r_2 * h6_p3 + fs_135_1859 * r_4 * h4_p3) + fs_1485_16 * e_2 * h4_p3;

        pc_64[k] = e_0 * (fs_9_22 * h6_p4 - fs_405_22 * r_2 * h4_p4) + e_1 * (fs_1_8450 * h8_p4 + fs_14_65 * h8_p8 - fs_2_275 * r_2 * h6_p4 + fs_405_3718 * r_4 * h4_p4) + fs_4455_32 * e_2 * h4_p4;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[65] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64};

    for (size_t m = 0; m < 65; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
