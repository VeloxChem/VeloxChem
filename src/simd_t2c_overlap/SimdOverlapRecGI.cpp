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



#include "SimdOverlapRecGI.hpp"

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
compute_gi_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGI.compute_gi_overlap: Basis functions must be of angular momenta four and six"));
    }

    if (harmonics.size() < 10)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGI.compute_gi_overlap: Harmonics must reach angular momentum ten"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGI.compute_gi_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 117 * nvalues, 0.0);

        return;
    }

    // NOTE: the first five rows accumulate the contracted prefactors of the terms,
    // and the remaining 117 rows hold the integrals of the combinations of angular
    // components.

    auto buffer = CSimdMatrix(122, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);

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

            const auto f_0 = fbase * aexp * aexp * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * aexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * aexp * aexp * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * aexp * aexp * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_4 = fbase * aexp * aexp / fexp / fexp / fexp / fexp / fexp / fexp;

            // NOTE: the exponential depends on the pair of primitives alone, so it is
            // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fss = std::exp(-fmu * ab_2[k]);

                pe_0[k] += f_0 * fss;
                pe_1[k] += f_1 * fss;
                pe_2[k] += f_2 * fss;
                pe_3[k] += f_3 * fss;
                pe_4[k] += f_4 * fss;
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
    const auto *ph10_m10 = harmonics[9].data(0);
    const auto *ph10_m9 = harmonics[9].data(1);
    const auto *ph10_m8 = harmonics[9].data(2);
    const auto *ph10_m7 = harmonics[9].data(3);
    const auto *ph10_m6 = harmonics[9].data(4);
    const auto *ph10_m5 = harmonics[9].data(5);
    const auto *ph10_m4 = harmonics[9].data(6);
    const auto *ph10_m3 = harmonics[9].data(7);
    const auto *ph10_m2 = harmonics[9].data(8);
    const auto *ph10_m1 = harmonics[9].data(9);
    const auto *ph10_0 = harmonics[9].data(10);
    const auto *ph10_p1 = harmonics[9].data(11);
    const auto *ph10_p2 = harmonics[9].data(12);
    const auto *ph10_p3 = harmonics[9].data(13);
    const auto *ph10_p4 = harmonics[9].data(14);
    const auto *ph10_p5 = harmonics[9].data(15);
    const auto *ph10_p6 = harmonics[9].data(16);
    const auto *ph10_p7 = harmonics[9].data(17);
    const auto *ph10_p8 = harmonics[9].data(18);
    const auto *ph10_p9 = harmonics[9].data(19);
    const auto *ph10_p10 = harmonics[9].data(20);

    auto *pc_0 = buffer.data(5);
    auto *pc_1 = buffer.data(6);
    auto *pc_2 = buffer.data(7);
    auto *pc_3 = buffer.data(8);
    auto *pc_4 = buffer.data(9);
    auto *pc_5 = buffer.data(10);
    auto *pc_6 = buffer.data(11);
    auto *pc_7 = buffer.data(12);
    auto *pc_8 = buffer.data(13);
    auto *pc_9 = buffer.data(14);
    auto *pc_10 = buffer.data(15);
    auto *pc_11 = buffer.data(16);
    auto *pc_12 = buffer.data(17);
    auto *pc_13 = buffer.data(18);
    auto *pc_14 = buffer.data(19);
    auto *pc_15 = buffer.data(20);
    auto *pc_16 = buffer.data(21);
    auto *pc_17 = buffer.data(22);
    auto *pc_18 = buffer.data(23);
    auto *pc_19 = buffer.data(24);
    auto *pc_20 = buffer.data(25);
    auto *pc_21 = buffer.data(26);
    auto *pc_22 = buffer.data(27);
    auto *pc_23 = buffer.data(28);
    auto *pc_24 = buffer.data(29);
    auto *pc_25 = buffer.data(30);
    auto *pc_26 = buffer.data(31);
    auto *pc_27 = buffer.data(32);
    auto *pc_28 = buffer.data(33);
    auto *pc_29 = buffer.data(34);
    auto *pc_30 = buffer.data(35);
    auto *pc_31 = buffer.data(36);
    auto *pc_32 = buffer.data(37);
    auto *pc_33 = buffer.data(38);
    auto *pc_34 = buffer.data(39);
    auto *pc_35 = buffer.data(40);
    auto *pc_36 = buffer.data(41);
    auto *pc_37 = buffer.data(42);
    auto *pc_38 = buffer.data(43);
    auto *pc_39 = buffer.data(44);
    auto *pc_40 = buffer.data(45);
    auto *pc_41 = buffer.data(46);
    auto *pc_42 = buffer.data(47);
    auto *pc_43 = buffer.data(48);
    auto *pc_44 = buffer.data(49);
    auto *pc_45 = buffer.data(50);
    auto *pc_46 = buffer.data(51);
    auto *pc_47 = buffer.data(52);
    auto *pc_48 = buffer.data(53);
    auto *pc_49 = buffer.data(54);
    auto *pc_50 = buffer.data(55);
    auto *pc_51 = buffer.data(56);
    auto *pc_52 = buffer.data(57);
    auto *pc_53 = buffer.data(58);
    auto *pc_54 = buffer.data(59);
    auto *pc_55 = buffer.data(60);
    auto *pc_56 = buffer.data(61);
    auto *pc_57 = buffer.data(62);
    auto *pc_58 = buffer.data(63);
    auto *pc_59 = buffer.data(64);
    auto *pc_60 = buffer.data(65);
    auto *pc_61 = buffer.data(66);
    auto *pc_62 = buffer.data(67);
    auto *pc_63 = buffer.data(68);
    auto *pc_64 = buffer.data(69);
    auto *pc_65 = buffer.data(70);
    auto *pc_66 = buffer.data(71);
    auto *pc_67 = buffer.data(72);
    auto *pc_68 = buffer.data(73);
    auto *pc_69 = buffer.data(74);
    auto *pc_70 = buffer.data(75);
    auto *pc_71 = buffer.data(76);
    auto *pc_72 = buffer.data(77);
    auto *pc_73 = buffer.data(78);
    auto *pc_74 = buffer.data(79);
    auto *pc_75 = buffer.data(80);
    auto *pc_76 = buffer.data(81);
    auto *pc_77 = buffer.data(82);
    auto *pc_78 = buffer.data(83);
    auto *pc_79 = buffer.data(84);
    auto *pc_80 = buffer.data(85);
    auto *pc_81 = buffer.data(86);
    auto *pc_82 = buffer.data(87);
    auto *pc_83 = buffer.data(88);
    auto *pc_84 = buffer.data(89);
    auto *pc_85 = buffer.data(90);
    auto *pc_86 = buffer.data(91);
    auto *pc_87 = buffer.data(92);
    auto *pc_88 = buffer.data(93);
    auto *pc_89 = buffer.data(94);
    auto *pc_90 = buffer.data(95);
    auto *pc_91 = buffer.data(96);
    auto *pc_92 = buffer.data(97);
    auto *pc_93 = buffer.data(98);
    auto *pc_94 = buffer.data(99);
    auto *pc_95 = buffer.data(100);
    auto *pc_96 = buffer.data(101);
    auto *pc_97 = buffer.data(102);
    auto *pc_98 = buffer.data(103);
    auto *pc_99 = buffer.data(104);
    auto *pc_100 = buffer.data(105);
    auto *pc_101 = buffer.data(106);
    auto *pc_102 = buffer.data(107);
    auto *pc_103 = buffer.data(108);
    auto *pc_104 = buffer.data(109);
    auto *pc_105 = buffer.data(110);
    auto *pc_106 = buffer.data(111);
    auto *pc_107 = buffer.data(112);
    auto *pc_108 = buffer.data(113);
    auto *pc_109 = buffer.data(114);
    auto *pc_110 = buffer.data(115);
    auto *pc_111 = buffer.data(116);
    auto *pc_112 = buffer.data(117);
    auto *pc_113 = buffer.data(118);
    auto *pc_114 = buffer.data(119);
    auto *pc_115 = buffer.data(120);
    auto *pc_116 = buffer.data(121);

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_7425_8 = std::sqrt(928.125);
    const auto fs_111375_8 = std::sqrt(13921.875);
    const auto fs_1575_176 = std::sqrt(1575.0 / 176.0);
    const auto fs_6075_22 = std::sqrt(6075.0 / 22.0);
    const auto fs_12375_8 = std::sqrt(1546.875);
    const auto fs_21_1859 = std::sqrt(21.0 / 1859.0);
    const auto fs_7_11 = std::sqrt(7.0 / 11.0);
    const auto fs_12150_1859 = std::sqrt(12150.0 / 1859.0);
    const auto fs_250_11 = std::sqrt(250.0 / 11.0);
    const auto fs_49_35263202 = std::sqrt(49.0 / 35263202.0);
    const auto fs_735_4199 = std::sqrt(735.0 / 4199.0);
    const auto fs_84_671099 = std::sqrt(84.0 / 671099.0);
    const auto fs_7_3179 = std::sqrt(7.0 / 3179.0);
    const auto fs_24_1859 = std::sqrt(24.0 / 1859.0);
    const auto fs_125_3718 = std::sqrt(125.0 / 3718.0);
    const auto fs_5457375_512 = std::sqrt(10658.935546875);
    const auto fs_22275_16 = std::sqrt(1392.1875);
    const auto fs_37125_8 = std::sqrt(4640.625);
    const auto fs_2625_88 = std::sqrt(2625.0 / 88.0);
    const auto fs_18225_44 = std::sqrt(18225.0 / 44.0);
    const auto fs_4125_8 = std::sqrt(515.625);
    const auto fs_245_3718 = std::sqrt(245.0 / 3718.0);
    const auto fs_70_33 = std::sqrt(70.0 / 33.0);
    const auto fs_18225_1859 = std::sqrt(18225.0 / 1859.0);
    const auto fs_250_33 = std::sqrt(250.0 / 33.0);
    const auto fs_441_35263202 = std::sqrt(441.0 / 35263202.0);
    const auto fs_441_4199 = std::sqrt(441.0 / 4199.0);
    const auto fs_490_671099 = std::sqrt(490.0 / 671099.0);
    const auto fs_70_9537 = std::sqrt(70.0 / 9537.0);
    const auto fs_36_1859 = std::sqrt(36.0 / 1859.0);
    const auto fs_125_11154 = std::sqrt(125.0 / 11154.0);
    const auto fs_1819125_512 = std::sqrt(3552.978515625);
    const auto fs_10125_4 = std::sqrt(2531.25);
    const auto fs_55125_484 = std::sqrt(55125.0 / 484.0);
    const auto fs_91125_121 = std::sqrt(91125.0 / 121.0);
    const auto fs_1125_4 = std::sqrt(281.25);
    const auto fs_8820_20449 = std::sqrt(8820.0 / 20449.0);
    const auto fs_196_143 = std::sqrt(196.0 / 143.0);
    const auto fs_980_121 = std::sqrt(980.0 / 121.0);
    const auto fs_364500_20449 = std::sqrt(364500.0 / 20449.0);
    const auto fs_500_121 = std::sqrt(500.0 / 121.0);
    const auto fs_2205_17631601 = std::sqrt(2205.0 / 17631601.0);
    const auto fs_4851_79781 = std::sqrt(4851.0 / 79781.0);
    const auto fs_35280_7382089 = std::sqrt(35280.0 / 7382089.0);
    const auto fs_784_51623 = std::sqrt(784.0 / 51623.0);
    const auto fs_980_34969 = std::sqrt(980.0 / 34969.0);
    const auto fs_720_20449 = std::sqrt(720.0 / 20449.0);
    const auto fs_125_20449 = std::sqrt(125.0 / 20449.0);
    const auto fs_496125_256 = std::sqrt(1937.98828125);
    const auto fs_3375_4 = std::sqrt(843.75);
    const auto fs_2025_8 = std::sqrt(253.125);
    const auto fs_77175_968 = std::sqrt(77175.0 / 968.0);
    const auto fs_30375_121 = std::sqrt(30375.0 / 121.0);
    const auto fs_225_8 = std::sqrt(28.125);
    const auto fs_10584_20449 = std::sqrt(10584.0 / 20449.0);
    const auto fs_735_286 = std::sqrt(735.0 / 286.0);
    const auto fs_686_121 = std::sqrt(686.0 / 121.0);
    const auto fs_121500_20449 = std::sqrt(121500.0 / 20449.0);
    const auto fs_50_121 = std::sqrt(50.0 / 121.0);
    const auto fs_8085_35263202 = std::sqrt(8085.0 / 35263202.0);
    const auto fs_2695_79781 = std::sqrt(2695.0 / 79781.0);
    const auto fs_42336_7382089 = std::sqrt(42336.0 / 7382089.0);
    const auto fs_1470_51623 = std::sqrt(1470.0 / 51623.0);
    const auto fs_686_34969 = std::sqrt(686.0 / 34969.0);
    const auto fs_240_20449 = std::sqrt(240.0 / 20449.0);
    const auto fs_25_40898 = std::sqrt(25.0 / 40898.0);
    const auto fs_99225_512 = std::sqrt(193.798828125);
    const auto fs_3375_8 = std::sqrt(421.875);
    const auto fs_42875_484 = std::sqrt(42875.0 / 484.0);
    const auto fs_30375_242 = std::sqrt(30375.0 / 242.0);
    const auto fs_25_8 = std::sqrt(3.125);
    const auto fs_20580_20449 = std::sqrt(20580.0 / 20449.0);
    const auto fs_441_143 = std::sqrt(441.0 / 143.0);
    const auto fs_6860_1089 = std::sqrt(6860.0 / 1089.0);
    const auto fs_60750_20449 = std::sqrt(60750.0 / 20449.0);
    const auto fs_50_1089 = std::sqrt(50.0 / 1089.0);
    const auto fs_24255_35263202 = std::sqrt(24255.0 / 35263202.0);
    const auto fs_24255_1356277 = std::sqrt(24255.0 / 1356277.0);
    const auto fs_82320_7382089 = std::sqrt(82320.0 / 7382089.0);
    const auto fs_1764_51623 = std::sqrt(1764.0 / 51623.0);
    const auto fs_6860_314721 = std::sqrt(6860.0 / 314721.0);
    const auto fs_120_20449 = std::sqrt(120.0 / 20449.0);
    const auto fs_25_368082 = std::sqrt(25.0 / 368082.0);
    const auto fs_11025_512 = std::sqrt(21.533203125);
    const auto fs_4725_32 = std::sqrt(147.65625);
    const auto fs_42525_968 = std::sqrt(42525.0 / 968.0);
    const auto fs_3087_1859 = std::sqrt(3087.0 / 1859.0);
    const auto fs_420_143 = std::sqrt(420.0 / 143.0);
    const auto fs_42525_40898 = std::sqrt(42525.0 / 40898.0);
    const auto fs_4851_2712554 = std::sqrt(4851.0 / 2712554.0);
    const auto fs_24255_2712554 = std::sqrt(24255.0 / 2712554.0);
    const auto fs_12348_671099 = std::sqrt(12348.0 / 671099.0);
    const auto fs_1680_51623 = std::sqrt(1680.0 / 51623.0);
    const auto fs_42_20449 = std::sqrt(42.0 / 20449.0);
    const auto f_15_2 = 7.5;
    const auto f_45_11 = 45.0 / 11.0;
    const auto fs_8820_1859 = std::sqrt(8820.0 / 1859.0);
    const auto f_90_143 = 90.0 / 143.0;
    const auto fs_11319_1356277 = std::sqrt(11319.0 / 1356277.0);
    const auto fs_35280_671099 = std::sqrt(35280.0 / 671099.0);
    const auto f_4_143 = 4.0 / 143.0;
    const auto fs_51975_32 = std::sqrt(1624.21875);
    const auto fs_14175_352 = std::sqrt(14175.0 / 352.0);
    const auto fs_42525_88 = std::sqrt(42525.0 / 88.0);
    const auto fs_63_676 = std::sqrt(63.0 / 676.0);
    const auto fs_63_22 = std::sqrt(63.0 / 22.0);
    const auto fs_42525_3718 = std::sqrt(42525.0 / 3718.0);
    const auto fs_49_2712554 = std::sqrt(49.0 / 2712554.0);
    const auto fs_294_4199 = std::sqrt(294.0 / 4199.0);
    const auto fs_63_61009 = std::sqrt(63.0 / 61009.0);
    const auto fs_63_6358 = std::sqrt(63.0 / 6358.0);
    const auto fs_42_1859 = std::sqrt(42.0 / 1859.0);
    const auto fs_22275_64 = std::sqrt(348.046875);
    const auto fs_37125_4 = std::sqrt(9281.25);
    const auto fs_25725_352 = std::sqrt(25725.0 / 352.0);
    const auto fs_18225_176 = std::sqrt(18225.0 / 176.0);
    const auto fs_4125_4 = std::sqrt(1031.25);
    const auto fs_5887_14872 = std::sqrt(5887.0 / 14872.0);
    const auto fs_49_13 = std::sqrt(49.0 / 13.0);
    const auto fs_343_66 = std::sqrt(343.0 / 66.0);
    const auto fs_18225_7436 = std::sqrt(18225.0 / 7436.0);
    const auto fs_500_33 = std::sqrt(500.0 / 33.0);
    const auto fs_2352_17631601 = std::sqrt(2352.0 / 17631601.0);
    const auto fs_7056_79781 = std::sqrt(7056.0 / 79781.0);
    const auto fs_5887_1342198 = std::sqrt(5887.0 / 1342198.0);
    const auto fs_196_4693 = std::sqrt(196.0 / 4693.0);
    const auto fs_343_19074 = std::sqrt(343.0 / 19074.0);
    const auto fs_9_1859 = std::sqrt(9.0 / 1859.0);
    const auto fs_125_5577 = std::sqrt(125.0 / 5577.0);
    const auto fs_1819125_256 = std::sqrt(7105.95703125);
    const auto f_45_8 = 5.625;
    const auto fs_6750 = std::sqrt(6750.0);
    const auto fs_65625_968 = std::sqrt(65625.0 / 968.0);
    const auto f_135_44 = 135.0 / 44.0;
    const auto fs_750 = std::sqrt(750.0);
    const auto fs_153125_163592 = std::sqrt(153125.0 / 163592.0);
    const auto fs_2401_1144 = std::sqrt(2401.0 / 1144.0);
    const auto fs_1750_363 = std::sqrt(1750.0 / 363.0);
    const auto f_135_286 = 135.0 / 286.0;
    const auto fs_4000_363 = std::sqrt(4000.0 / 363.0);
    const auto fs_9702_17631601 = std::sqrt(9702.0 / 17631601.0);
    const auto fs_6468_79781 = std::sqrt(6468.0 / 79781.0);
    const auto fs_153125_14764178 = std::sqrt(153125.0 / 14764178.0);
    const auto fs_2401_103246 = std::sqrt(2401.0 / 103246.0);
    const auto fs_1750_104907 = std::sqrt(1750.0 / 104907.0);
    const auto f_3_143 = 3.0 / 143.0;
    const auto fs_1000_61347 = std::sqrt(1000.0 / 61347.0);
    const auto fs_165375_32 = std::sqrt(5167.96875);
    const auto fs_16875_16 = std::sqrt(1054.6875);
    const auto fs_6075 = std::sqrt(6075.0);
    const auto fs_33075_484 = std::sqrt(33075.0 / 484.0);
    const auto fs_151875_484 = std::sqrt(151875.0 / 484.0);
    const auto fs_675 = std::sqrt(675.0);
    const auto fs_64827_20449 = std::sqrt(64827.0 / 20449.0);
    const auto fs_441_1144 = std::sqrt(441.0 / 1144.0);
    const auto fs_588_121 = std::sqrt(588.0 / 121.0);
    const auto fs_151875_20449 = std::sqrt(151875.0 / 20449.0);
    const auto fs_1200_121 = std::sqrt(1200.0 / 121.0);
    const auto fs_58800_17631601 = std::sqrt(58800.0 / 17631601.0);
    const auto fs_86240_1356277 = std::sqrt(86240.0 / 1356277.0);
    const auto fs_259308_7382089 = std::sqrt(259308.0 / 7382089.0);
    const auto fs_441_103246 = std::sqrt(441.0 / 103246.0);
    const auto fs_588_34969 = std::sqrt(588.0 / 34969.0);
    const auto fs_300_20449 = std::sqrt(300.0 / 20449.0);
    const auto fs_297675_64 = std::sqrt(4651.171875);
    const auto fs_30375_32 = std::sqrt(949.21875);
    const auto f_30 = 30.0;
    const auto fs_8575_1936 = std::sqrt(8575.0 / 1936.0);
    const auto fs_273375_968 = std::sqrt(273375.0 / 968.0);
    const auto f_10 = 10.0;
    const auto fs_42483_20449 = std::sqrt(42483.0 / 20449.0);
    const auto fs_21_572 = std::sqrt(21.0 / 572.0);
    const auto fs_343_1089 = std::sqrt(343.0 / 1089.0);
    const auto fs_273375_40898 = std::sqrt(273375.0 / 40898.0);
    const auto f_40_33 = 40.0 / 33.0;
    const auto fs_72765_17631601 = std::sqrt(72765.0 / 17631601.0);
    const auto fs_121275_2712554 = std::sqrt(121275.0 / 2712554.0);
    const auto fs_169932_7382089 = std::sqrt(169932.0 / 7382089.0);
    const auto fs_21_51623 = std::sqrt(21.0 / 51623.0);
    const auto fs_343_314721 = std::sqrt(343.0 / 314721.0);
    const auto fs_270_20449 = std::sqrt(270.0 / 20449.0);
    const auto f_20_429 = 20.0 / 429.0;
    const auto f_105_4 = 26.25;
    const auto fs_114075_128 = std::sqrt(891.2109375);
    const auto fs_1125_8 = std::sqrt(140.625);
    const auto fs_1026675_3872 = std::sqrt(1026675.0 / 3872.0);
    const auto fs_125_8 = std::sqrt(15.625);
    const auto fs_1029_484 = std::sqrt(1029.0 / 484.0);
    const auto fs_2625_3718 = std::sqrt(2625.0 / 3718.0);
    const auto fs_6075_968 = std::sqrt(6075.0 / 968.0);
    const auto fs_250_1089 = std::sqrt(250.0 / 1089.0);
    const auto fs_155232_17631601 = std::sqrt(155232.0 / 17631601.0);
    const auto fs_38808_1356277 = std::sqrt(38808.0 / 1356277.0);
    const auto fs_1029_43681 = std::sqrt(1029.0 / 43681.0);
    const auto fs_5250_671099 = std::sqrt(5250.0 / 671099.0);
    const auto fs_3_242 = std::sqrt(3.0 / 242.0);
    const auto fs_125_368082 = std::sqrt(125.0 / 368082.0);
    const auto fs_55125_512 = std::sqrt(107.666015625);
    const auto f_255_8 = 31.875;
    const auto f_765_44 = 765.0 / 44.0;
    const auto fs_11907_3718 = std::sqrt(11907.0 / 3718.0);
    const auto f_765_286 = 765.0 / 286.0;
    const auto fs_45276_1356277 = std::sqrt(45276.0 / 1356277.0);
    const auto fs_23814_671099 = std::sqrt(23814.0 / 671099.0);
    const auto f_17_143 = 17.0 / 143.0;
    const auto fs_30375_352 = std::sqrt(30375.0 / 352.0);
    const auto fs_135_338 = std::sqrt(135.0 / 338.0);
    const auto fs_42_13 = std::sqrt(42.0 / 13.0);
    const auto fs_135_22 = std::sqrt(135.0 / 22.0);
    const auto fs_343_2712554 = std::sqrt(343.0 / 2712554.0);
    const auto fs_2058_79781 = std::sqrt(2058.0 / 79781.0);
    const auto fs_270_61009 = std::sqrt(270.0 / 61009.0);
    const auto fs_168_4693 = std::sqrt(168.0 / 4693.0);
    const auto fs_135_6358 = std::sqrt(135.0 / 6358.0);
    const auto fs_675_11 = std::sqrt(675.0 / 11.0);
    const auto fs_1587_1352 = std::sqrt(1587.0 / 1352.0);
    const auto fs_7_104 = std::sqrt(7.0 / 104.0);
    const auto fs_48_11 = std::sqrt(48.0 / 11.0);
    const auto fs_1029_1356277 = std::sqrt(1029.0 / 1356277.0);
    const auto fs_4116_79781 = std::sqrt(4116.0 / 79781.0);
    const auto fs_1587_122018 = std::sqrt(1587.0 / 122018.0);
    const auto fs_7_9386 = std::sqrt(7.0 / 9386.0);
    const auto fs_48_3179 = std::sqrt(48.0 / 3179.0);
    const auto fs_14175_16 = std::sqrt(885.9375);
    const auto fs_23625_4 = std::sqrt(5906.25);
    const auto fs_39675_3872 = std::sqrt(39675.0 / 3872.0);
    const auto fs_127575_484 = std::sqrt(127575.0 / 484.0);
    const auto fs_2625_4 = std::sqrt(656.25);
    const auto fs_38642_20449 = std::sqrt(38642.0 / 20449.0);
    const auto fs_210_143 = std::sqrt(210.0 / 143.0);
    const auto fs_529_726 = std::sqrt(529.0 / 726.0);
    const auto fs_127575_20449 = std::sqrt(127575.0 / 20449.0);
    const auto fs_3500_363 = std::sqrt(3500.0 / 363.0);
    const auto fs_45276_17631601 = std::sqrt(45276.0 / 17631601.0);
    const auto fs_90552_1356277 = std::sqrt(90552.0 / 1356277.0);
    const auto fs_154568_7382089 = std::sqrt(154568.0 / 7382089.0);
    const auto fs_840_51623 = std::sqrt(840.0 / 51623.0);
    const auto fs_529_209814 = std::sqrt(529.0 / 209814.0);
    const auto fs_252_20449 = std::sqrt(252.0 / 20449.0);
    const auto fs_875_61347 = std::sqrt(875.0 / 61347.0);
    const auto fs_1157625_256 = std::sqrt(4521.97265625);
    const auto fs_23625_64 = std::sqrt(369.140625);
    const auto fs_14175_2 = std::sqrt(7087.5);
    const auto fs_450_121 = std::sqrt(450.0 / 121.0);
    const auto fs_212625_1936 = std::sqrt(212625.0 / 1936.0);
    const auto fs_1575_2 = std::sqrt(787.5);
    const auto fs_338709_163592 = std::sqrt(338709.0 / 163592.0);
    const auto fs_2187_1144 = std::sqrt(2187.0 / 1144.0);
    const auto fs_32_121 = std::sqrt(32.0 / 121.0);
    const auto fs_212625_81796 = std::sqrt(212625.0 / 81796.0);
    const auto fs_1400_121 = std::sqrt(1400.0 / 121.0);
    const auto fs_113190_17631601 = std::sqrt(113190.0 / 17631601.0);
    const auto fs_94325_1356277 = std::sqrt(94325.0 / 1356277.0);
    const auto fs_338709_14764178 = std::sqrt(338709.0 / 14764178.0);
    const auto fs_2187_103246 = std::sqrt(2187.0 / 103246.0);
    const auto fs_32_34969 = std::sqrt(32.0 / 34969.0);
    const auto fs_105_20449 = std::sqrt(105.0 / 20449.0);
    const auto fs_350_20449 = std::sqrt(350.0 / 20449.0);
    const auto fs_694575_128 = std::sqrt(5426.3671875);
    const auto fs_9450 = std::sqrt(9450.0);
    const auto fs_525_8 = std::sqrt(65.625);
    const auto fs_1050 = std::sqrt(1050.0);
    const auto fs_378_121 = std::sqrt(378.0 / 121.0);
    const auto fs_4107_3718 = std::sqrt(4107.0 / 3718.0);
    const auto fs_14_3 = std::sqrt(14.0 / 3.0);
    const auto fs_5600_363 = std::sqrt(5600.0 / 363.0);
    const auto fs_463050_17631601 = std::sqrt(463050.0 / 17631601.0);
    const auto fs_169785_2712554 = std::sqrt(169785.0 / 2712554.0);
    const auto fs_1512_43681 = std::sqrt(1512.0 / 43681.0);
    const auto fs_8214_671099 = std::sqrt(8214.0 / 671099.0);
    const auto fs_14_867 = std::sqrt(14.0 / 867.0);
    const auto fs_1400_61347 = std::sqrt(1400.0 / 61347.0);
    const auto fs_231525_32 = std::sqrt(7235.15625);
    const auto fs_42525_128 = std::sqrt(332.2265625);
    const auto fs_7875_4 = std::sqrt(1968.75);
    const auto fs_6125_121 = std::sqrt(6125.0 / 121.0);
    const auto fs_382725_3872 = std::sqrt(382725.0 / 3872.0);
    const auto fs_875_4 = std::sqrt(218.75);
    const auto fs_55545_81796 = std::sqrt(55545.0 / 81796.0);
    const auto fs_9_44 = std::sqrt(9.0 / 44.0);
    const auto fs_3920_1089 = std::sqrt(3920.0 / 1089.0);
    const auto fs_382725_163592 = std::sqrt(382725.0 / 163592.0);
    const auto fs_3500_1089 = std::sqrt(3500.0 / 1089.0);
    const auto fs_407484_17631601 = std::sqrt(407484.0 / 17631601.0);
    const auto fs_67914_1356277 = std::sqrt(67914.0 / 1356277.0);
    const auto fs_55545_7382089 = std::sqrt(55545.0 / 7382089.0);
    const auto fs_9_3971 = std::sqrt(9.0 / 3971.0);
    const auto fs_3920_314721 = std::sqrt(3920.0 / 314721.0);
    const auto fs_189_40898 = std::sqrt(189.0 / 40898.0);
    const auto fs_875_184041 = std::sqrt(875.0 / 184041.0);
    const auto fs_385875_256 = std::sqrt(1507.32421875);
    const auto f_165_4 = 41.25;
    const auto f_45_2 = 22.5;
    const auto fs_375_4 = std::sqrt(93.75);
    const auto fs_2016_20449 = std::sqrt(2016.0 / 20449.0);
    const auto f_45_13 = 45.0 / 13.0;
    const auto fs_500_363 = std::sqrt(500.0 / 363.0);
    const auto fs_1267728_17631601 = std::sqrt(1267728.0 / 17631601.0);
    const auto fs_8064_7382089 = std::sqrt(8064.0 / 7382089.0);
    const auto f_2_13 = 2.0 / 13.0;
    const auto fs_125_61347 = std::sqrt(125.0 / 61347.0);
    const auto fs_165375_256 = std::sqrt(645.99609375);
    const auto fs_3375_32 = std::sqrt(105.46875);
    const auto fs_15_13 = std::sqrt(15.0 / 13.0);
    const auto fs_189_52 = std::sqrt(189.0 / 52.0);
    const auto fs_15_2 = std::sqrt(7.5);
    const auto fs_1715_2712554 = std::sqrt(1715.0 / 2712554.0);
    const auto fs_686_79781 = std::sqrt(686.0 / 79781.0);
    const auto fs_60_4693 = std::sqrt(60.0 / 4693.0);
    const auto fs_189_4693 = std::sqrt(189.0 / 4693.0);
    const auto fs_15_578 = std::sqrt(15.0 / 578.0);
    const auto fs_1125_176 = std::sqrt(1125.0 / 176.0);
    const auto fs_1445_676 = std::sqrt(1445.0 / 676.0);
    const auto fs_105_104 = std::sqrt(105.0 / 104.0);
    const auto fs_5_11 = std::sqrt(5.0 / 11.0);
    const auto fs_4116_1356277 = std::sqrt(4116.0 / 1356277.0);
    const auto fs_32928_1356277 = std::sqrt(32928.0 / 1356277.0);
    const auto fs_1445_61009 = std::sqrt(1445.0 / 61009.0);
    const auto fs_105_9386 = std::sqrt(105.0 / 9386.0);
    const auto fs_5_3179 = std::sqrt(5.0 / 3179.0);
    const auto fs_33075_1936 = std::sqrt(33075.0 / 1936.0);
    const auto fs_30603_14872 = std::sqrt(30603.0 / 14872.0);
    const auto fs_5_1144 = std::sqrt(5.0 / 1144.0);
    const auto fs_147_121 = std::sqrt(147.0 / 121.0);
    const auto fs_56595_1356277 = std::sqrt(56595.0 / 1356277.0);
    const auto fs_30603_1342198 = std::sqrt(30603.0 / 1342198.0);
    const auto fs_5_103246 = std::sqrt(5.0 / 103246.0);
    const auto fs_147_34969 = std::sqrt(147.0 / 34969.0);
    const auto fs_14175_4 = std::sqrt(3543.75);
    const auto fs_190125_3872 = std::sqrt(190125.0 / 3872.0);
    const auto fs_1575_4 = std::sqrt(393.75);
    const auto fs_184815_163592 = std::sqrt(184815.0 / 163592.0);
    const auto fs_5547_7436 = std::sqrt(5547.0 / 7436.0);
    const auto fs_845_242 = std::sqrt(845.0 / 242.0);
    const auto fs_700_121 = std::sqrt(700.0 / 121.0);
    const auto fs_301840_17631601 = std::sqrt(301840.0 / 17631601.0);
    const auto fs_75460_1356277 = std::sqrt(75460.0 / 1356277.0);
    const auto fs_184815_14764178 = std::sqrt(184815.0 / 14764178.0);
    const auto fs_5547_671099 = std::sqrt(5547.0 / 671099.0);
    const auto fs_845_69938 = std::sqrt(845.0 / 69938.0);
    const auto fs_175_20449 = std::sqrt(175.0 / 20449.0);
    const auto fs_694575_256 = std::sqrt(2713.18359375);
    const auto fs_23625_32 = std::sqrt(738.28125);
    const auto fs_6300 = std::sqrt(6300.0);
    const auto f_265_44 = 265.0 / 44.0;
    const auto fs_212625_968 = std::sqrt(212625.0 / 968.0);
    const auto fs_700 = std::sqrt(700.0);
    const auto fs_17661_81796 = std::sqrt(17661.0 / 81796.0);
    const auto fs_2880_1859 = std::sqrt(2880.0 / 1859.0);
    const auto f_53_33 = 53.0 / 33.0;
    const auto fs_212625_40898 = std::sqrt(212625.0 / 40898.0);
    const auto fs_11200_1089 = std::sqrt(11200.0 / 1089.0);
    const auto fs_509355_17631601 = std::sqrt(509355.0 / 17631601.0);
    const auto fs_17661_7382089 = std::sqrt(17661.0 / 7382089.0);
    const auto fs_11520_671099 = std::sqrt(11520.0 / 671099.0);
    const auto f_53_561 = 53.0 / 561.0;
    const auto fs_210_20449 = std::sqrt(210.0 / 20449.0);
    const auto fs_2800_184041 = std::sqrt(2800.0 / 184041.0);
    const auto fs_77175_16 = std::sqrt(4823.4375);
    const auto fs_23625_2 = std::sqrt(11812.5);
    const auto fs_7875_8 = std::sqrt(984.375);
    const auto fs_2625_242 = std::sqrt(2625.0 / 242.0);
    const auto fs_2625_2 = std::sqrt(1312.5);
    const auto fs_875_8 = std::sqrt(109.375);
    const auto fs_1890_20449 = std::sqrt(1890.0 / 20449.0);
    const auto fs_118803_81796 = std::sqrt(118803.0 / 81796.0);
    const auto fs_280_363 = std::sqrt(280.0 / 363.0);
    const auto fs_7000_363 = std::sqrt(7000.0 / 363.0);
    const auto fs_1750_1089 = std::sqrt(1750.0 / 1089.0);
    const auto fs_1481760_17631601 = std::sqrt(1481760.0 / 17631601.0);
    const auto fs_1086624_17631601 = std::sqrt(1086624.0 / 17631601.0);
    const auto fs_7560_7382089 = std::sqrt(7560.0 / 7382089.0);
    const auto fs_118803_7382089 = std::sqrt(118803.0 / 7382089.0);
    const auto fs_280_104907 = std::sqrt(280.0 / 104907.0);
    const auto fs_1750_61347 = std::sqrt(1750.0 / 61347.0);
    const auto fs_875_368082 = std::sqrt(875.0 / 368082.0);
    const auto fs_1157625_128 = std::sqrt(9043.9453125);
    const auto fs_385875_512 = std::sqrt(753.662109375);
    const auto f_15_8 = 1.875;
    const auto f_45_44 = 45.0 / 44.0;
    const auto fs_55125_40898 = std::sqrt(55125.0 / 40898.0);
    const auto f_45_286 = 45.0 / 286.0;
    const auto fs_1901592_17631601 = std::sqrt(1901592.0 / 17631601.0);
    const auto fs_110250_7382089 = std::sqrt(110250.0 / 7382089.0);
    const auto f_1_143 = 1.0 / 143.0;
    const auto f_45_4 = 11.25;
    const auto fs_63_13 = std::sqrt(63.0 / 13.0);
    const auto f_3 = 3.0;
    const auto fs_6860_1356277 = std::sqrt(6860.0 / 1356277.0);
    const auto fs_252_4693 = std::sqrt(252.0 / 4693.0);
    const auto f_3_17 = 3.0 / 17.0;
    const auto fs_121_26 = std::sqrt(121.0 / 26.0);
    const auto f_2 = 2.0;
    const auto fs_25725_1356277 = std::sqrt(25725.0 / 1356277.0);
    const auto fs_242_4693 = std::sqrt(242.0 / 4693.0);
    const auto f_2_17 = 2.0 / 17.0;
    const auto f_120_11 = 120.0 / 11.0;
    const auto fs_3364_1859 = std::sqrt(3364.0 / 1859.0);
    const auto f_32_11 = 32.0 / 11.0;
    const auto fs_13456_671099 = std::sqrt(13456.0 / 671099.0);
    const auto f_32_187 = 32.0 / 187.0;
    const auto f_135_22 = 135.0 / 22.0;
    const auto fs_225_3718 = std::sqrt(225.0 / 3718.0);
    const auto f_18_11 = 18.0 / 11.0;
    const auto fs_450_671099 = std::sqrt(450.0 / 671099.0);
    const auto f_18_187 = 18.0 / 187.0;
    const auto fs_7875_2 = std::sqrt(3937.5);
    const auto f_5_4 = 1.25;
    const auto fs_875_2 = std::sqrt(437.5);
    const auto fs_75_121 = std::sqrt(75.0 / 121.0);
    const auto f_1_3 = 1.0 / 3.0;
    const auto fs_7000_1089 = std::sqrt(7000.0 / 1089.0);
    const auto fs_1697850_17631601 = std::sqrt(1697850.0 / 17631601.0);
    const auto fs_300_43681 = std::sqrt(300.0 / 43681.0);
    const auto f_1_51 = 1.0 / 51.0;
    const auto fs_1750_184041 = std::sqrt(1750.0 / 184041.0);
    const auto fs_385875_128 = std::sqrt(3014.6484375);
    const auto fs_39375_4 = std::sqrt(9843.75);
    const auto f_80_11 = 80.0 / 11.0;
    const auto fs_4375_4 = std::sqrt(1093.75);
    const auto fs_46389_20449 = std::sqrt(46389.0 / 20449.0);
    const auto f_64_33 = 64.0 / 33.0;
    const auto fs_17500_1089 = std::sqrt(17500.0 / 1089.0);
    const auto fs_2037420_17631601 = std::sqrt(2037420.0 / 17631601.0);
    const auto fs_185556_7382089 = std::sqrt(185556.0 / 7382089.0);
    const auto f_64_561 = 64.0 / 561.0;
    const auto fs_4375_184041 = std::sqrt(4375.0 / 184041.0);
    const auto fs_1929375_256 = std::sqrt(7536.62109375);
    const auto f_75_2 = 37.5;
    const auto f_225_2 = 112.5;
    const auto f_105_11 = 105.0 / 11.0;
    const auto f_225_11 = 225.0 / 11.0;
    const auto f_252_143 = 252.0 / 143.0;
    const auto f_28_11 = 28.0 / 11.0;
    const auto f_450_143 = 450.0 / 143.0;
    const auto f_50_11 = 50.0 / 11.0;
    const auto f_1470_4199 = 1470.0 / 4199.0;
    const auto f_504_2717 = 504.0 / 2717.0;
    const auto f_28_187 = 28.0 / 187.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_25_143 = 25.0 / 143.0;
    const auto f_1575_16 = 98.4375;

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph10_p1, ph10_p2, ph10_p9, ph10_p10, ab_2, pc_0, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h10_p10 = ph10_p10[k];

        pc_0[k] = e_0 * (fs_7425_8 * h4_p2 - fs_111375_8 * r_2 * h2_p2) + e_1 * (fs_1575_176 * h6_p2 - fs_6075_22 * r_2 * h4_p2 + fs_12375_8 * r_4 * h2_p2) + e_2 * (fs_21_1859 * h8_p2 - fs_7_11 * r_2 * h6_p2 + fs_12150_1859 * r_4 * h4_p2 - fs_250_11 * r_6 * h2_p2) + e_3 * (fs_49_35263202 * h10_p2 - fs_735_4199 * h10_p10 - fs_84_671099 * r_2 * h8_p2 + fs_7_3179 * r_4 * h6_p2 - fs_24_1859 * r_6 * h4_p2 + fs_125_3718 * r_8 * h2_p2) + fs_5457375_512 * e_4 * h2_p2;

        pc_1[k] = e_0 * (fs_22275_16 * h4_p1 - fs_37125_8 * r_2 * h2_p1) + e_1 * (fs_2625_88 * h6_p1 - fs_18225_44 * r_2 * h4_p1 + fs_4125_8 * r_4 * h2_p1) + e_2 * (fs_245_3718 * h8_p1 - fs_70_33 * r_2 * h6_p1 + fs_18225_1859 * r_4 * h4_p1 - fs_250_33 * r_6 * h2_p1) + e_3 * (fs_441_35263202 * h10_p1 - fs_441_4199 * h10_p9 - fs_490_671099 * r_2 * h8_p1 + fs_70_9537 * r_4 * h6_p1 - fs_36_1859 * r_6 * h4_p1 + fs_125_11154 * r_8 * h2_p1) + fs_1819125_512 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_2[k] = e_0 * (fs_10125_4 * h4_0 - fs_10125_4 * r_2 * h2_0) + e_1 * (fs_55125_484 * h6_0 - fs_91125_121 * r_2 * h4_0 + fs_1125_4 * r_4 * h2_0) + e_2 * (fs_8820_20449 * h8_0 - fs_196_143 * h8_p8 - fs_980_121 * r_2 * h6_0 + fs_364500_20449 * r_4 * h4_0 - fs_500_121 * r_6 * h2_0) + e_3 * (fs_2205_17631601 * h10_0 - fs_4851_79781 * h10_p8 - fs_35280_7382089 * r_2 * h8_0 + fs_784_51623 * r_2 * h8_p8 + fs_980_34969 * r_4 * h6_0 - fs_720_20449 * r_6 * h4_0 + fs_125_20449 * r_8 * h2_0) + fs_496125_256 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_3[k] = e_0 * (-fs_3375_4 * h4_p1 + fs_2025_8 * r_2 * h2_p1) + e_1 * (-fs_77175_968 * h6_p1 + fs_30375_121 * r_2 * h4_p1 - fs_225_8 * r_4 * h2_p1) + e_2 * (-fs_10584_20449 * h8_p1 - fs_735_286 * h8_p7 + fs_686_121 * r_2 * h6_p1 - fs_121500_20449 * r_4 * h4_p1 + fs_50_121 * r_6 * h2_p1) + e_3 * (-fs_8085_35263202 * h10_p1 - fs_2695_79781 * h10_p7 + fs_42336_7382089 * r_2 * h8_p1 + fs_1470_51623 * r_2 * h8_p7 - fs_686_34969 * r_4 * h6_p1 + fs_240_20449 * r_6 * h4_p1 - fs_25_40898 * r_8 * h2_p1) - fs_99225_512 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_4[k] = e_0 * (fs_3375_8 * h4_p2 - fs_225_8 * r_2 * h2_p2) + e_1 * (fs_42875_484 * h6_p2 - fs_1575_176 * h6_p6 - fs_30375_242 * r_2 * h4_p2 + fs_25_8 * r_4 * h2_p2) + e_2 * (fs_20580_20449 * h8_p2 - fs_441_143 * h8_p6 - fs_6860_1089 * r_2 * h6_p2 + fs_7_11 * r_2 * h6_p6 + fs_60750_20449 * r_4 * h4_p2 - fs_50_1089 * r_6 * h2_p2) + e_3 * (fs_24255_35263202 * h10_p2 - fs_24255_1356277 * h10_p6 - fs_82320_7382089 * r_2 * h8_p2 + fs_1764_51623 * r_2 * h8_p6 + fs_6860_314721 * r_4 * h6_p2 - fs_7_3179 * r_4 * h6_p6 - fs_120_20449 * r_6 * h4_p2 + fs_25_368082 * r_8 * h2_p2) + fs_11025_512 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_p3, ph6_m4, ph6_p3, ph6_p5, ph8_m4, ph8_p3, ph8_p5, ph10_m4, ph10_p3, ph10_p5, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_5[k] = -fs_4725_32 * e_0 * h4_p3 + e_1 * (-fs_77175_968 * h6_p3 - fs_2625_88 * h6_p5 + fs_42525_968 * r_2 * h4_p3) + e_2 * (-fs_3087_1859 * h8_p3 - fs_420_143 * h8_p5 + fs_686_121 * r_2 * h6_p3 + fs_70_33 * r_2 * h6_p5 - fs_42525_40898 * r_4 * h4_p3) + e_3 * (-fs_4851_2712554 * h10_p3 - fs_24255_2712554 * h10_p5 + fs_12348_671099 * r_2 * h8_p3 + fs_1680_51623 * r_2 * h8_p5 - fs_686_34969 * r_4 * h6_p3 - fs_70_9537 * r_4 * h6_p5 + fs_42_20449 * r_6 * h4_p3);

        pc_6[k] = f_15_2 * e_0 * h4_m4 + e_1 * (fs_55125_484 * h6_m4 - f_45_11 * r_2 * h4_m4) + e_2 * (fs_8820_1859 * h8_m4 - fs_980_121 * r_2 * h6_m4 + f_90_143 * r_4 * h4_m4) + e_3 * (fs_11319_1356277 * h10_m4 - fs_35280_671099 * r_2 * h8_m4 + fs_980_34969 * r_4 * h6_m4 - f_4_143 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_7[k] = -fs_4725_32 * e_0 * h4_m3 + e_1 * (fs_2625_88 * h6_m5 - fs_77175_968 * h6_m3 + fs_42525_968 * r_2 * h4_m3) + e_2 * (fs_420_143 * h8_m5 - fs_3087_1859 * h8_m3 - fs_70_33 * r_2 * h6_m5 + fs_686_121 * r_2 * h6_m3 - fs_42525_40898 * r_4 * h4_m3) + e_3 * (fs_24255_2712554 * h10_m5 - fs_4851_2712554 * h10_m3 - fs_1680_51623 * r_2 * h8_m5 + fs_12348_671099 * r_2 * h8_m3 + fs_70_9537 * r_4 * h6_m5 - fs_686_34969 * r_4 * h6_m3 + fs_42_20449 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_8[k] = e_0 * (fs_3375_8 * h4_m2 - fs_225_8 * r_2 * h2_m2) + e_1 * (fs_1575_176 * h6_m6 + fs_42875_484 * h6_m2 - fs_30375_242 * r_2 * h4_m2 + fs_25_8 * r_4 * h2_m2) + e_2 * (fs_441_143 * h8_m6 + fs_20580_20449 * h8_m2 - fs_7_11 * r_2 * h6_m6 - fs_6860_1089 * r_2 * h6_m2 + fs_60750_20449 * r_4 * h4_m2 - fs_50_1089 * r_6 * h2_m2) + e_3 * (fs_24255_1356277 * h10_m6 + fs_24255_35263202 * h10_m2 - fs_1764_51623 * r_2 * h8_m6 - fs_82320_7382089 * r_2 * h8_m2 + fs_7_3179 * r_4 * h6_m6 + fs_6860_314721 * r_4 * h6_m2 - fs_120_20449 * r_6 * h4_m2 + fs_25_368082 * r_8 * h2_m2) + fs_11025_512 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m9, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_9, pc_10, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_9[k] = e_0 * (-fs_3375_4 * h4_m1 + fs_2025_8 * r_2 * h2_m1) + e_1 * (-fs_77175_968 * h6_m1 + fs_30375_121 * r_2 * h4_m1 - fs_225_8 * r_4 * h2_m1) + e_2 * (fs_735_286 * h8_m7 - fs_10584_20449 * h8_m1 + fs_686_121 * r_2 * h6_m1 - fs_121500_20449 * r_4 * h4_m1 + fs_50_121 * r_6 * h2_m1) + e_3 * (fs_2695_79781 * h10_m7 - fs_8085_35263202 * h10_m1 - fs_1470_51623 * r_2 * h8_m7 + fs_42336_7382089 * r_2 * h8_m1 - fs_686_34969 * r_4 * h6_m1 + fs_240_20449 * r_6 * h4_m1 - fs_25_40898 * r_8 * h2_m1) - fs_99225_512 * e_4 * h2_m1;

        pc_10[k] = fs_196_143 * e_2 * h8_m8 + e_3 * (fs_4851_79781 * h10_m8 - fs_784_51623 * r_2 * h8_m8);

        pc_11[k] = e_0 * (-fs_22275_16 * h4_m1 + fs_37125_8 * r_2 * h2_m1) + e_1 * (-fs_2625_88 * h6_m1 + fs_18225_44 * r_2 * h4_m1 - fs_4125_8 * r_4 * h2_m1) + e_2 * (-fs_245_3718 * h8_m1 + fs_70_33 * r_2 * h6_m1 - fs_18225_1859 * r_4 * h4_m1 + fs_250_33 * r_6 * h2_m1) + e_3 * (fs_441_4199 * h10_m9 - fs_441_35263202 * h10_m1 + fs_490_671099 * r_2 * h8_m1 - fs_70_9537 * r_4 * h6_m1 + fs_36_1859 * r_6 * h4_m1 - fs_125_11154 * r_8 * h2_m1) - fs_1819125_512 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph4_p3, ph6_m2, ph6_p3, ph8_m2, ph8_p3, ph10_m10, ph10_m2, ph10_p3, ph10_p9, ab_2, pc_12, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];

        pc_12[k] = e_0 * (-fs_7425_8 * h4_m2 + fs_111375_8 * r_2 * h2_m2) + e_1 * (-fs_1575_176 * h6_m2 + fs_6075_22 * r_2 * h4_m2 - fs_12375_8 * r_4 * h2_m2) + e_2 * (-fs_21_1859 * h8_m2 + fs_7_11 * r_2 * h6_m2 - fs_12150_1859 * r_4 * h4_m2 + fs_250_11 * r_6 * h2_m2) + e_3 * (fs_735_4199 * h10_m10 - fs_49_35263202 * h10_m2 + fs_84_671099 * r_2 * h8_m2 - fs_7_3179 * r_4 * h6_m2 + fs_24_1859 * r_6 * h4_m2 - fs_125_3718 * r_8 * h2_m2) - fs_5457375_512 * e_4 * h2_m2;

        pc_13[k] = -fs_51975_32 * e_0 * h4_p3 + e_1 * (-fs_14175_352 * h6_p3 + fs_42525_88 * r_2 * h4_p3) + e_2 * (-fs_63_676 * h8_p3 + fs_63_22 * r_2 * h6_p3 - fs_42525_3718 * r_4 * h4_p3) + e_3 * (-fs_49_2712554 * h10_p3 - fs_294_4199 * h10_p9 + fs_63_61009 * r_2 * h8_p3 - fs_63_6358 * r_4 * h6_p3 + fs_42_1859 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_14[k] = e_0 * (-fs_22275_64 * h4_p2 - fs_37125_4 * r_2 * h2_p2) + e_1 * (-fs_25725_352 * h6_p2 + fs_18225_176 * r_2 * h4_p2 + fs_4125_4 * r_4 * h2_p2) + e_2 * (-fs_5887_14872 * h8_p2 + fs_49_13 * h8_p8 + fs_343_66 * r_2 * h6_p2 - fs_18225_7436 * r_4 * h4_p2 - fs_500_33 * r_6 * h2_p2) + e_3 * (-fs_2352_17631601 * h10_p2 - fs_7056_79781 * h10_p8 + fs_5887_1342198 * r_2 * h8_p2 - fs_196_4693 * r_2 * h8_p8 - fs_343_19074 * r_4 * h6_p2 + fs_9_1859 * r_6 * h4_p2 + fs_125_5577 * r_8 * h2_p2) + fs_1819125_256 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_15[k] = e_0 * (f_45_8 * h4_p1 - fs_6750 * r_2 * h2_p1) + e_1 * (-fs_65625_968 * h6_p1 - f_135_44 * r_2 * h4_p1 + fs_750 * r_4 * h2_p1) + e_2 * (-fs_153125_163592 * h8_p1 + fs_2401_1144 * h8_p7 + fs_1750_363 * r_2 * h6_p1 + f_135_286 * r_4 * h4_p1 - fs_4000_363 * r_6 * h2_p1) + e_3 * (-fs_9702_17631601 * h10_p1 - fs_6468_79781 * h10_p7 + fs_153125_14764178 * r_2 * h8_p1 - fs_2401_103246 * r_2 * h8_p7 - fs_1750_104907 * r_4 * h6_p1 - f_3_143 * r_6 * h4_p1 + fs_1000_61347 * r_8 * h2_p1) + fs_165375_32 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_16[k] = e_0 * (fs_16875_16 * h4_0 - fs_6075 * r_2 * h2_0) + e_1 * (-fs_33075_484 * h6_0 + fs_14175_352 * h6_p6 - fs_151875_484 * r_2 * h4_0 + fs_675 * r_4 * h2_0) + e_2 * (-fs_64827_20449 * h8_0 + fs_441_1144 * h8_p6 + fs_588_121 * r_2 * h6_0 - fs_63_22 * r_2 * h6_p6 + fs_151875_20449 * r_4 * h4_0 - fs_1200_121 * r_6 * h2_0) + e_3 * (-fs_58800_17631601 * h10_0 - fs_86240_1356277 * h10_p6 + fs_259308_7382089 * r_2 * h8_0 - fs_441_103246 * r_2 * h8_p6 - fs_588_34969 * r_4 * h6_0 + fs_63_6358 * r_4 * h6_p6 - fs_300_20449 * r_6 * h4_0 + fs_300_20449 * r_8 * h2_0) + fs_297675_64 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_17[k] = e_0 * (-fs_30375_32 * h4_p1 + f_30 * r_2 * h2_p1) + e_1 * (fs_8575_1936 * h6_p1 + fs_25725_352 * h6_p5 + fs_273375_968 * r_2 * h4_p1 - f_10 * r_4 * h2_p1) + e_2 * (fs_42483_20449 * h8_p1 - fs_21_572 * h8_p5 - fs_343_1089 * r_2 * h6_p1 - fs_343_66 * r_2 * h6_p5 - fs_273375_40898 * r_4 * h4_p1 + f_40_33 * r_6 * h2_p1) + e_3 * (fs_72765_17631601 * h10_p1 - fs_121275_2712554 * h10_p5 - fs_169932_7382089 * r_2 * h8_p1 + fs_21_51623 * r_2 * h8_p5 + fs_343_314721 * r_4 * h6_p1 + fs_343_19074 * r_4 * h6_p5 + fs_270_20449 * r_6 * h4_p1 - f_20_429 * r_8 * h2_p1) - f_105_4 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_18[k] = e_0 * (fs_114075_128 * h4_p2 + fs_4725_32 * h4_p4 - fs_1125_8 * r_2 * h2_p2) + e_1 * (fs_8575_1936 * h6_p2 + fs_65625_968 * h6_p4 - fs_1026675_3872 * r_2 * h4_p2 - fs_42525_968 * r_2 * h4_p4 + fs_125_8 * r_4 * h2_p2) + e_2 * (-fs_1029_484 * h8_p2 - fs_2625_3718 * h8_p4 - fs_343_1089 * r_2 * h6_p2 - fs_1750_363 * r_2 * h6_p4 + fs_6075_968 * r_4 * h4_p2 + fs_42525_40898 * r_4 * h4_p4 - fs_250_1089 * r_6 * h2_p2) + e_3 * (-fs_155232_17631601 * h10_p2 - fs_38808_1356277 * h10_p4 + fs_1029_43681 * r_2 * h8_p2 + fs_5250_671099 * r_2 * h8_p4 + fs_343_314721 * r_4 * h6_p2 + fs_1750_104907 * r_4 * h6_p4 - fs_3_242 * r_6 * h4_p2 - fs_42_20449 * r_6 * h4_p4 + fs_125_368082 * r_8 * h2_p2) + fs_55125_512 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ab_2, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];

        pc_19[k] = -f_255_8 * e_0 * h4_m3 + e_1 * (-fs_33075_484 * h6_m3 + f_765_44 * r_2 * h4_m3) + e_2 * (fs_11907_3718 * h8_m3 + fs_588_121 * r_2 * h6_m3 - f_765_286 * r_4 * h4_m3) + e_3 * (fs_45276_1356277 * h10_m3 - fs_23814_671099 * r_2 * h8_m3 - fs_588_34969 * r_4 * h6_m3 + f_17_143 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_20[k] = e_0 * (-fs_4725_32 * h4_m4 + fs_114075_128 * h4_m2 - fs_1125_8 * r_2 * h2_m2) + e_1 * (-fs_65625_968 * h6_m4 + fs_8575_1936 * h6_m2 + fs_42525_968 * r_2 * h4_m4 - fs_1026675_3872 * r_2 * h4_m2 + fs_125_8 * r_4 * h2_m2) + e_2 * (fs_2625_3718 * h8_m4 - fs_1029_484 * h8_m2 + fs_1750_363 * r_2 * h6_m4 - fs_343_1089 * r_2 * h6_m2 - fs_42525_40898 * r_4 * h4_m4 + fs_6075_968 * r_4 * h4_m2 - fs_250_1089 * r_6 * h2_m2) + e_3 * (fs_38808_1356277 * h10_m4 - fs_155232_17631601 * h10_m2 - fs_5250_671099 * r_2 * h8_m4 + fs_1029_43681 * r_2 * h8_m2 - fs_1750_104907 * r_4 * h6_m4 + fs_343_314721 * r_4 * h6_m2 + fs_42_20449 * r_6 * h4_m4 - fs_3_242 * r_6 * h4_m2 + fs_125_368082 * r_8 * h2_m2) + fs_55125_512 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ph10_m6, ph10_m5, ph10_m1, ab_2, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_21[k] = e_0 * (-fs_30375_32 * h4_m1 + f_30 * r_2 * h2_m1) + e_1 * (-fs_25725_352 * h6_m5 + fs_8575_1936 * h6_m1 + fs_273375_968 * r_2 * h4_m1 - f_10 * r_4 * h2_m1) + e_2 * (fs_21_572 * h8_m5 + fs_42483_20449 * h8_m1 + fs_343_66 * r_2 * h6_m5 - fs_343_1089 * r_2 * h6_m1 - fs_273375_40898 * r_4 * h4_m1 + f_40_33 * r_6 * h2_m1) + e_3 * (fs_121275_2712554 * h10_m5 + fs_72765_17631601 * h10_m1 - fs_21_51623 * r_2 * h8_m5 - fs_169932_7382089 * r_2 * h8_m1 - fs_343_19074 * r_4 * h6_m5 + fs_343_314721 * r_4 * h6_m1 + fs_270_20449 * r_6 * h4_m1 - f_20_429 * r_8 * h2_m1) - f_105_4 * e_4 * h2_m1;

        pc_22[k] = -fs_14175_352 * e_1 * h6_m6 + e_2 * (-fs_441_1144 * h8_m6 + fs_63_22 * r_2 * h6_m6) + e_3 * (fs_86240_1356277 * h10_m6 + fs_441_103246 * r_2 * h8_m6 - fs_63_6358 * r_4 * h6_m6);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_23[k] = e_0 * (-f_45_8 * h4_m1 + fs_6750 * r_2 * h2_m1) + e_1 * (fs_65625_968 * h6_m1 + f_135_44 * r_2 * h4_m1 - fs_750 * r_4 * h2_m1) + e_2 * (-fs_2401_1144 * h8_m7 + fs_153125_163592 * h8_m1 - fs_1750_363 * r_2 * h6_m1 - f_135_286 * r_4 * h4_m1 + fs_4000_363 * r_6 * h2_m1) + e_3 * (fs_6468_79781 * h10_m7 + fs_9702_17631601 * h10_m1 + fs_2401_103246 * r_2 * h8_m7 - fs_153125_14764178 * r_2 * h8_m1 + fs_1750_104907 * r_4 * h6_m1 + f_3_143 * r_6 * h4_m1 - fs_1000_61347 * r_8 * h2_m1) - fs_165375_32 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m3, ph4_m2, ph6_m3, ph6_m2, ph8_m8, ph8_m3, ph8_m2, ph10_m9, ph10_m8, ph10_m3, ph10_m2, ab_2, pc_24, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m2 = ph10_m2[k];

        pc_24[k] = e_0 * (fs_22275_64 * h4_m2 + fs_37125_4 * r_2 * h2_m2) + e_1 * (fs_25725_352 * h6_m2 - fs_18225_176 * r_2 * h4_m2 - fs_4125_4 * r_4 * h2_m2) + e_2 * (-fs_49_13 * h8_m8 + fs_5887_14872 * h8_m2 - fs_343_66 * r_2 * h6_m2 + fs_18225_7436 * r_4 * h4_m2 + fs_500_33 * r_6 * h2_m2) + e_3 * (fs_7056_79781 * h10_m8 + fs_2352_17631601 * h10_m2 + fs_196_4693 * r_2 * h8_m8 - fs_5887_1342198 * r_2 * h8_m2 + fs_343_19074 * r_4 * h6_m2 - fs_9_1859 * r_6 * h4_m2 - fs_125_5577 * r_8 * h2_m2) - fs_1819125_256 * e_4 * h2_m2;

        pc_25[k] = fs_51975_32 * e_0 * h4_m3 + e_1 * (fs_14175_352 * h6_m3 - fs_42525_88 * r_2 * h4_m3) + e_2 * (fs_63_676 * h8_m3 - fs_63_22 * r_2 * h6_m3 + fs_42525_3718 * r_4 * h4_m3) + e_3 * (fs_294_4199 * h10_m9 + fs_49_2712554 * h10_m3 - fs_63_61009 * r_2 * h8_m3 + fs_63_6358 * r_4 * h6_m3 - fs_42_1859 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ph8_p7, ph8_p8, ph10_p3, ph10_p4, ph10_p7, ph10_p8, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h10_p8 = ph10_p8[k];

        pc_26[k] = fs_7425_8 * e_0 * h4_p4 + e_1 * (fs_30375_352 * h6_p4 - fs_6075_22 * r_2 * h4_p4) + e_2 * (fs_135_338 * h8_p4 - fs_42_13 * h8_p8 - fs_135_22 * r_2 * h6_p4 + fs_12150_1859 * r_4 * h4_p4) + e_3 * (fs_343_2712554 * h10_p4 - fs_2058_79781 * h10_p8 - fs_270_61009 * r_2 * h8_p4 + fs_168_4693 * r_2 * h8_p8 + fs_135_6358 * r_4 * h6_p4 - fs_24_1859 * r_6 * h4_p4);

        pc_27[k] = -fs_22275_64 * e_0 * h4_p3 + e_1 * (fs_675_11 * h6_p3 + fs_18225_176 * r_2 * h4_p3) + e_2 * (fs_1587_1352 * h8_p3 + fs_7_104 * h8_p7 - fs_48_11 * r_2 * h6_p3 - fs_18225_7436 * r_4 * h4_p3) + e_3 * (fs_1029_1356277 * h10_p3 - fs_4116_79781 * h10_p7 - fs_1587_122018 * r_2 * h8_p3 - fs_7_9386 * r_2 * h8_p7 + fs_48_3179 * r_4 * h6_p3 + fs_9_1859 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_28[k] = e_0 * (-fs_14175_16 * h4_p2 - fs_23625_4 * r_2 * h2_p2) + e_1 * (fs_39675_3872 * h6_p2 - fs_30375_352 * h6_p6 + fs_127575_484 * r_2 * h4_p2 + fs_2625_4 * r_4 * h2_p2) + e_2 * (fs_38642_20449 * h8_p2 + fs_210_143 * h8_p6 - fs_529_726 * r_2 * h6_p2 + fs_135_22 * r_2 * h6_p6 - fs_127575_20449 * r_4 * h4_p2 - fs_3500_363 * r_6 * h2_p2) + e_3 * (fs_45276_17631601 * h10_p2 - fs_90552_1356277 * h10_p6 - fs_154568_7382089 * r_2 * h8_p2 - fs_840_51623 * r_2 * h8_p6 + fs_529_209814 * r_4 * h6_p2 - fs_135_6358 * r_4 * h6_p6 + fs_252_20449 * r_6 * h4_p2 + fs_875_61347 * r_8 * h2_p2) + fs_1157625_256 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_29[k] = e_0 * (-fs_23625_64 * h4_p1 - fs_14175_2 * r_2 * h2_p1) + e_1 * (-fs_450_121 * h6_p1 - fs_675_11 * h6_p5 + fs_212625_1936 * r_2 * h4_p1 + fs_1575_2 * r_4 * h2_p1) + e_2 * (fs_338709_163592 * h8_p1 + fs_2187_1144 * h8_p5 + fs_32_121 * r_2 * h6_p1 + fs_48_11 * r_2 * h6_p5 - fs_212625_81796 * r_4 * h4_p1 - fs_1400_121 * r_6 * h2_p1) + e_3 * (fs_113190_17631601 * h10_p1 - fs_94325_1356277 * h10_p5 - fs_338709_14764178 * r_2 * h8_p1 - fs_2187_103246 * r_2 * h8_p5 - fs_32_34969 * r_4 * h6_p1 - fs_48_3179 * r_4 * h6_p5 + fs_105_20449 * r_6 * h4_p1 + fs_350_20449 * r_8 * h2_p1) + fs_694575_128 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_30[k] = e_0 * (-fs_3375_8 * h4_p4 - fs_9450 * r_2 * h2_0) + e_1 * (-fs_525_8 * h6_0 - fs_39675_3872 * h6_p4 + fs_30375_242 * r_2 * h4_p4 + fs_1050 * r_4 * h2_0) + e_2 * (fs_378_121 * h8_0 + fs_4107_3718 * h8_p4 + fs_14_3 * r_2 * h6_0 + fs_529_726 * r_2 * h6_p4 - fs_60750_20449 * r_4 * h4_p4 - fs_5600_363 * r_6 * h2_0) + e_3 * (fs_463050_17631601 * h10_0 - fs_169785_2712554 * h10_p4 - fs_1512_43681 * r_2 * h8_0 - fs_8214_671099 * r_2 * h8_p4 - fs_14_867 * r_4 * h6_0 - fs_529_209814 * r_4 * h6_p4 + fs_120_20449 * r_6 * h4_p4 + fs_1400_61347 * r_8 * h2_0) + fs_231525_32 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_31[k] = e_0 * (-fs_42525_128 * h4_p1 - fs_114075_128 * h4_p3 + fs_7875_4 * r_2 * h2_p1) + e_1 * (fs_6125_121 * h6_p1 + fs_450_121 * h6_p3 + fs_382725_3872 * r_2 * h4_p1 + fs_1026675_3872 * r_2 * h4_p3 - fs_875_4 * r_4 * h2_p1) + e_2 * (-fs_55545_81796 * h8_p1 + fs_9_44 * h8_p3 - fs_3920_1089 * r_2 * h6_p1 - fs_32_121 * r_2 * h6_p3 - fs_382725_163592 * r_4 * h4_p1 - fs_6075_968 * r_4 * h4_p3 + fs_3500_1089 * r_6 * h2_p1) + e_3 * (-fs_407484_17631601 * h10_p1 - fs_67914_1356277 * h10_p3 + fs_55545_7382089 * r_2 * h8_p1 - fs_9_3971 * r_2 * h8_p3 + fs_3920_314721 * r_4 * h6_p1 + fs_32_34969 * r_4 * h6_p3 + fs_189_40898 * r_6 * h4_p1 + fs_3_242 * r_6 * h4_p3 - fs_875_184041 * r_8 * h2_p1) - fs_385875_256 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ab_2, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];

        pc_32[k] = e_0 * (f_165_4 * h4_m2 - fs_3375_4 * r_2 * h2_m2) + e_1 * (-fs_525_8 * h6_m2 - f_45_2 * r_2 * h4_m2 + fs_375_4 * r_4 * h2_m2) + e_2 * (fs_2016_20449 * h8_m2 + fs_14_3 * r_2 * h6_m2 + f_45_13 * r_4 * h4_m2 - fs_500_363 * r_6 * h2_m2) + e_3 * (fs_1267728_17631601 * h10_m2 - fs_8064_7382089 * r_2 * h8_m2 - fs_14_867 * r_4 * h6_m2 - f_2_13 * r_6 * h4_m2 + fs_125_61347 * r_8 * h2_m2) + fs_165375_256 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_33[k] = e_0 * (fs_114075_128 * h4_m3 - fs_42525_128 * h4_m1 + fs_7875_4 * r_2 * h2_m1) + e_1 * (-fs_450_121 * h6_m3 + fs_6125_121 * h6_m1 - fs_1026675_3872 * r_2 * h4_m3 + fs_382725_3872 * r_2 * h4_m1 - fs_875_4 * r_4 * h2_m1) + e_2 * (-fs_9_44 * h8_m3 - fs_55545_81796 * h8_m1 + fs_32_121 * r_2 * h6_m3 - fs_3920_1089 * r_2 * h6_m1 + fs_6075_968 * r_4 * h4_m3 - fs_382725_163592 * r_4 * h4_m1 + fs_3500_1089 * r_6 * h2_m1) + e_3 * (fs_67914_1356277 * h10_m3 - fs_407484_17631601 * h10_m1 + fs_9_3971 * r_2 * h8_m3 + fs_55545_7382089 * r_2 * h8_m1 - fs_32_34969 * r_4 * h6_m3 + fs_3920_314721 * r_4 * h6_m1 - fs_3_242 * r_6 * h4_m3 + fs_189_40898 * r_6 * h4_m1 - fs_875_184041 * r_8 * h2_m1) - fs_385875_256 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m1, ph10_m5, ph10_m4, ph10_m1, ab_2, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m1 = ph10_m1[k];

        pc_34[k] = fs_3375_8 * e_0 * h4_m4 + e_1 * (fs_39675_3872 * h6_m4 - fs_30375_242 * r_2 * h4_m4) + e_2 * (-fs_4107_3718 * h8_m4 - fs_529_726 * r_2 * h6_m4 + fs_60750_20449 * r_4 * h4_m4) + e_3 * (fs_169785_2712554 * h10_m4 + fs_8214_671099 * r_2 * h8_m4 + fs_529_209814 * r_4 * h6_m4 - fs_120_20449 * r_6 * h4_m4);

        pc_35[k] = e_0 * (fs_23625_64 * h4_m1 + fs_14175_2 * r_2 * h2_m1) + e_1 * (fs_675_11 * h6_m5 + fs_450_121 * h6_m1 - fs_212625_1936 * r_2 * h4_m1 - fs_1575_2 * r_4 * h2_m1) + e_2 * (-fs_2187_1144 * h8_m5 - fs_338709_163592 * h8_m1 - fs_48_11 * r_2 * h6_m5 - fs_32_121 * r_2 * h6_m1 + fs_212625_81796 * r_4 * h4_m1 + fs_1400_121 * r_6 * h2_m1) + e_3 * (fs_94325_1356277 * h10_m5 - fs_113190_17631601 * h10_m1 + fs_2187_103246 * r_2 * h8_m5 + fs_338709_14764178 * r_2 * h8_m1 + fs_48_3179 * r_4 * h6_m5 + fs_32_34969 * r_4 * h6_m1 - fs_105_20449 * r_6 * h4_m1 - fs_350_20449 * r_8 * h2_m1) - fs_694575_128 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_36[k] = e_0 * (fs_14175_16 * h4_m2 + fs_23625_4 * r_2 * h2_m2) + e_1 * (fs_30375_352 * h6_m6 - fs_39675_3872 * h6_m2 - fs_127575_484 * r_2 * h4_m2 - fs_2625_4 * r_4 * h2_m2) + e_2 * (-fs_210_143 * h8_m6 - fs_38642_20449 * h8_m2 - fs_135_22 * r_2 * h6_m6 + fs_529_726 * r_2 * h6_m2 + fs_127575_20449 * r_4 * h4_m2 + fs_3500_363 * r_6 * h2_m2) + e_3 * (fs_90552_1356277 * h10_m6 - fs_45276_17631601 * h10_m2 + fs_840_51623 * r_2 * h8_m6 + fs_154568_7382089 * r_2 * h8_m2 + fs_135_6358 * r_4 * h6_m6 - fs_529_209814 * r_4 * h6_m2 - fs_252_20449 * r_6 * h4_m2 - fs_875_61347 * r_8 * h2_m2) - fs_1157625_256 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph6_m4, ph6_m3, ph8_m8, ph8_m7, ph8_m4, ph8_m3, ph10_m8, ph10_m7, ph10_m4, ph10_m3, ab_2, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m3 = ph10_m3[k];

        pc_37[k] = fs_22275_64 * e_0 * h4_m3 + e_1 * (-fs_675_11 * h6_m3 - fs_18225_176 * r_2 * h4_m3) + e_2 * (-fs_7_104 * h8_m7 - fs_1587_1352 * h8_m3 + fs_48_11 * r_2 * h6_m3 + fs_18225_7436 * r_4 * h4_m3) + e_3 * (fs_4116_79781 * h10_m7 - fs_1029_1356277 * h10_m3 + fs_7_9386 * r_2 * h8_m7 + fs_1587_122018 * r_2 * h8_m3 - fs_48_3179 * r_4 * h6_m3 - fs_9_1859 * r_6 * h4_m3);

        pc_38[k] = -fs_7425_8 * e_0 * h4_m4 + e_1 * (-fs_30375_352 * h6_m4 + fs_6075_22 * r_2 * h4_m4) + e_2 * (fs_42_13 * h8_m8 - fs_135_338 * h8_m4 + fs_135_22 * r_2 * h6_m4 - fs_12150_1859 * r_4 * h4_m4) + e_3 * (fs_2058_79781 * h10_m8 - fs_343_2712554 * h10_m4 - fs_168_4693 * r_2 * h8_m8 + fs_270_61009 * r_2 * h8_m4 - fs_135_6358 * r_4 * h6_m4 + fs_24_1859 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p4, ph6_p4, ph6_p5, ph6_p6, ph8_p4, ph8_p5, ph8_p6, ph8_p7, ph10_p4, ph10_p5, ph10_p6, ph10_p7, ab_2, pc_39, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h10_p7 = ph10_p7[k];

        pc_39[k] = -fs_3375_32 * e_1 * h6_p5 + e_2 * (-fs_15_13 * h8_p5 - fs_189_52 * h8_p7 + fs_15_2 * r_2 * h6_p5) + e_3 * (-fs_1715_2712554 * h10_p5 - fs_686_79781 * h10_p7 + fs_60_4693 * r_2 * h8_p5 + fs_189_4693 * r_2 * h8_p7 - fs_15_578 * r_4 * h6_p5);

        pc_40[k] = fs_22275_16 * e_0 * h4_p4 + e_1 * (-fs_1125_176 * h6_p4 + fs_3375_32 * h6_p6 - fs_18225_44 * r_2 * h4_p4) + e_2 * (-fs_1445_676 * h8_p4 - fs_105_104 * h8_p6 + fs_5_11 * r_2 * h6_p4 - fs_15_2 * r_2 * h6_p6 + fs_18225_1859 * r_4 * h4_p4) + e_3 * (-fs_4116_1356277 * h10_p4 - fs_32928_1356277 * h10_p6 + fs_1445_61009 * r_2 * h8_p4 + fs_105_9386 * r_2 * h8_p6 - fs_5_3179 * r_4 * h6_p4 + fs_15_578 * r_4 * h6_p6 - fs_36_1859 * r_6 * h4_p4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_41[k] = f_45_8 * e_0 * h4_p3 + e_1 * (fs_33075_1936 * h6_p3 + fs_1125_176 * h6_p5 - f_135_44 * r_2 * h4_p3) + e_2 * (-fs_30603_14872 * h8_p3 + fs_5_1144 * h8_p5 - fs_147_121 * r_2 * h6_p3 - fs_5_11 * r_2 * h6_p5 + f_135_286 * r_4 * h4_p3) + e_3 * (-fs_11319_1356277 * h10_p3 - fs_56595_1356277 * h10_p5 + fs_30603_1342198 * r_2 * h8_p3 - fs_5_103246 * r_2 * h8_p5 + fs_147_34969 * r_4 * h6_p3 + fs_5_3179 * r_4 * h6_p5 - f_3_143 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_42[k] = e_0 * (-fs_23625_64 * h4_p2 + fs_3375_4 * h4_p4 - fs_14175_4 * r_2 * h2_p2) + e_1 * (fs_190125_3872 * h6_p2 - fs_33075_1936 * h6_p4 + fs_212625_1936 * r_2 * h4_p2 - fs_30375_121 * r_2 * h4_p4 + fs_1575_4 * r_4 * h2_p2) + e_2 * (-fs_184815_163592 * h8_p2 + fs_5547_7436 * h8_p4 - fs_845_242 * r_2 * h6_p2 + fs_147_121 * r_2 * h6_p4 - fs_212625_81796 * r_4 * h4_p2 + fs_121500_20449 * r_4 * h4_p4 - fs_700_121 * r_6 * h2_p2) + e_3 * (-fs_301840_17631601 * h10_p2 - fs_75460_1356277 * h10_p4 + fs_184815_14764178 * r_2 * h8_p2 - fs_5547_671099 * r_2 * h8_p4 + fs_845_69938 * r_4 * h6_p2 - fs_147_34969 * r_4 * h6_p4 + fs_105_20449 * r_6 * h4_p2 - fs_240_20449 * r_6 * h4_p4 + fs_175_20449 * r_8 * h2_p2) + fs_694575_256 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_43[k] = e_0 * (-fs_23625_32 * h4_p1 + fs_30375_32 * h4_p3 - fs_6300 * r_2 * h2_p1) + e_1 * (f_265_44 * h6_p1 - fs_190125_3872 * h6_p3 + fs_212625_968 * r_2 * h4_p1 - fs_273375_968 * r_2 * h4_p3 + fs_700 * r_4 * h2_p1) + e_2 * (-fs_17661_81796 * h8_p1 + fs_2880_1859 * h8_p3 - f_53_33 * r_2 * h6_p1 + fs_845_242 * r_2 * h6_p3 - fs_212625_40898 * r_4 * h4_p1 + fs_273375_40898 * r_4 * h4_p3 - fs_11200_1089 * r_6 * h2_p1) + e_3 * (-fs_509355_17631601 * h10_p1 - fs_169785_2712554 * h10_p3 + fs_17661_7382089 * r_2 * h8_p1 - fs_11520_671099 * r_2 * h8_p3 + f_53_561 * r_4 * h6_p1 - fs_845_69938 * r_4 * h6_p3 + fs_210_20449 * r_6 * h4_p1 - fs_270_20449 * r_6 * h4_p3 + fs_2800_184041 * r_8 * h2_p1) + fs_77175_16 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p2 = ph10_p2[k];

        pc_44[k] = e_0 * (-fs_23625_32 * h4_0 + fs_42525_128 * h4_p2 - fs_23625_2 * r_2 * h2_0 - fs_7875_8 * r_2 * h2_p2) + e_1 * (fs_2625_242 * h6_0 - f_265_44 * h6_p2 + fs_212625_968 * r_2 * h4_0 - fs_382725_3872 * r_2 * h4_p2 + fs_2625_2 * r_4 * h2_0 + fs_875_8 * r_4 * h2_p2) + e_2 * (fs_1890_20449 * h8_0 + fs_118803_81796 * h8_p2 - fs_280_363 * r_2 * h6_0 + f_53_33 * r_2 * h6_p2 - fs_212625_40898 * r_4 * h4_0 + fs_382725_163592 * r_4 * h4_p2 - fs_7000_363 * r_6 * h2_0 - fs_1750_1089 * r_6 * h2_p2) + e_3 * (-fs_1481760_17631601 * h10_0 - fs_1086624_17631601 * h10_p2 - fs_7560_7382089 * r_2 * h8_0 - fs_118803_7382089 * r_2 * h8_p2 + fs_280_104907 * r_4 * h6_0 - f_53_561 * r_4 * h6_p2 + fs_210_20449 * r_6 * h4_0 - fs_189_40898 * r_6 * h4_p2 + fs_1750_61347 * r_8 * h2_0 + fs_875_368082 * r_8 * h2_p2) + e_4 * (fs_1157625_128 * h2_0 + fs_385875_512 * h2_p2);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_m1 = ph10_m1[k];

        pc_45[k] = e_0 * (f_15_8 * h4_m1 + fs_6750 * r_2 * h2_m1) + e_1 * (fs_2625_242 * h6_m1 - f_45_44 * r_2 * h4_m1 - fs_750 * r_4 * h2_m1) + e_2 * (-fs_55125_40898 * h8_m1 - fs_280_363 * r_2 * h6_m1 + f_45_286 * r_4 * h4_m1 + fs_4000_363 * r_6 * h2_m1) + e_3 * (fs_1901592_17631601 * h10_m1 + fs_110250_7382089 * r_2 * h8_m1 + fs_280_104907 * r_4 * h6_m1 - f_1_143 * r_6 * h4_m1 - fs_1000_61347 * r_8 * h2_m1) - fs_165375_32 * e_4 * h2_m1;

        pc_46[k] = e_0 * (-fs_42525_128 * h4_m2 + fs_7875_8 * r_2 * h2_m2) + e_1 * (f_265_44 * h6_m2 + fs_382725_3872 * r_2 * h4_m2 - fs_875_8 * r_4 * h2_m2) + e_2 * (-fs_118803_81796 * h8_m2 - f_53_33 * r_2 * h6_m2 - fs_382725_163592 * r_4 * h4_m2 + fs_1750_1089 * r_6 * h2_m2) + e_3 * (fs_1086624_17631601 * h10_m2 + fs_118803_7382089 * r_2 * h8_m2 + f_53_561 * r_4 * h6_m2 + fs_189_40898 * r_6 * h4_m2 - fs_875_368082 * r_8 * h2_m2) - fs_385875_512 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_47[k] = e_0 * (-fs_30375_32 * h4_m3 + fs_23625_32 * h4_m1 + fs_6300 * r_2 * h2_m1) + e_1 * (fs_190125_3872 * h6_m3 - f_265_44 * h6_m1 + fs_273375_968 * r_2 * h4_m3 - fs_212625_968 * r_2 * h4_m1 - fs_700 * r_4 * h2_m1) + e_2 * (-fs_2880_1859 * h8_m3 + fs_17661_81796 * h8_m1 - fs_845_242 * r_2 * h6_m3 + f_53_33 * r_2 * h6_m1 - fs_273375_40898 * r_4 * h4_m3 + fs_212625_40898 * r_4 * h4_m1 + fs_11200_1089 * r_6 * h2_m1) + e_3 * (fs_169785_2712554 * h10_m3 + fs_509355_17631601 * h10_m1 + fs_11520_671099 * r_2 * h8_m3 - fs_17661_7382089 * r_2 * h8_m1 + fs_845_69938 * r_4 * h6_m3 - f_53_561 * r_4 * h6_m1 + fs_270_20449 * r_6 * h4_m3 - fs_210_20449 * r_6 * h4_m1 - fs_2800_184041 * r_8 * h2_m1) - fs_77175_16 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_48[k] = e_0 * (-fs_3375_4 * h4_m4 + fs_23625_64 * h4_m2 + fs_14175_4 * r_2 * h2_m2) + e_1 * (fs_33075_1936 * h6_m4 - fs_190125_3872 * h6_m2 + fs_30375_121 * r_2 * h4_m4 - fs_212625_1936 * r_2 * h4_m2 - fs_1575_4 * r_4 * h2_m2) + e_2 * (-fs_5547_7436 * h8_m4 + fs_184815_163592 * h8_m2 - fs_147_121 * r_2 * h6_m4 + fs_845_242 * r_2 * h6_m2 - fs_121500_20449 * r_4 * h4_m4 + fs_212625_81796 * r_4 * h4_m2 + fs_700_121 * r_6 * h2_m2) + e_3 * (fs_75460_1356277 * h10_m4 + fs_301840_17631601 * h10_m2 + fs_5547_671099 * r_2 * h8_m4 - fs_184815_14764178 * r_2 * h8_m2 + fs_147_34969 * r_4 * h6_m4 - fs_845_69938 * r_4 * h6_m2 + fs_240_20449 * r_6 * h4_m4 - fs_105_20449 * r_6 * h4_m2 - fs_175_20449 * r_8 * h2_m2) - fs_694575_256 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_49[k] = -f_45_8 * e_0 * h4_m3 + e_1 * (-fs_1125_176 * h6_m5 - fs_33075_1936 * h6_m3 + f_135_44 * r_2 * h4_m3) + e_2 * (-fs_5_1144 * h8_m5 + fs_30603_14872 * h8_m3 + fs_5_11 * r_2 * h6_m5 + fs_147_121 * r_2 * h6_m3 - f_135_286 * r_4 * h4_m3) + e_3 * (fs_56595_1356277 * h10_m5 + fs_11319_1356277 * h10_m3 + fs_5_103246 * r_2 * h8_m5 - fs_30603_1342198 * r_2 * h8_m3 - fs_5_3179 * r_4 * h6_m5 - fs_147_34969 * r_4 * h6_m3 + f_3_143 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph6_m6, ph6_m5, ph6_m4, ph8_m7, ph8_m6, ph8_m5, ph8_m4, ph10_m7, ph10_m6, ph10_m5, ph10_m4, ab_2, pc_50, pc_51, pc_52, pc_53, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];

        pc_50[k] = -fs_22275_16 * e_0 * h4_m4 + e_1 * (-fs_3375_32 * h6_m6 + fs_1125_176 * h6_m4 + fs_18225_44 * r_2 * h4_m4) + e_2 * (fs_105_104 * h8_m6 + fs_1445_676 * h8_m4 + fs_15_2 * r_2 * h6_m6 - fs_5_11 * r_2 * h6_m4 - fs_18225_1859 * r_4 * h4_m4) + e_3 * (fs_32928_1356277 * h10_m6 + fs_4116_1356277 * h10_m4 - fs_105_9386 * r_2 * h8_m6 - fs_1445_61009 * r_2 * h8_m4 - fs_15_578 * r_4 * h6_m6 + fs_5_3179 * r_4 * h6_m4 + fs_36_1859 * r_6 * h4_m4);

        pc_51[k] = fs_3375_32 * e_1 * h6_m5 + e_2 * (fs_189_52 * h8_m7 + fs_15_13 * h8_m5 - fs_15_2 * r_2 * h6_m5) + e_3 * (fs_686_79781 * h10_m7 + fs_1715_2712554 * h10_m5 - fs_189_4693 * r_2 * h8_m7 - fs_60_4693 * r_2 * h8_m5 + fs_15_578 * r_4 * h6_m5);

        pc_52[k] = f_45_4 * e_1 * h6_m6 + e_2 * (fs_63_13 * h8_m6 - f_3 * r_2 * h6_m6) + e_3 * (fs_6860_1356277 * h10_m6 - fs_252_4693 * r_2 * h8_m6 + f_3_17 * r_4 * h6_m6);

        pc_53[k] = -f_15_2 * e_1 * h6_m5 + e_2 * (fs_121_26 * h8_m5 + f_2 * r_2 * h6_m5) + e_3 * (fs_25725_1356277 * h10_m5 - fs_242_4693 * r_2 * h8_m5 - f_2_17 * r_4 * h6_m5);

        pc_54[k] = fs_10125_4 * e_0 * h4_m4 + e_1 * (-f_120_11 * h6_m4 - fs_91125_121 * r_2 * h4_m4) + e_2 * (fs_3364_1859 * h8_m4 + f_32_11 * r_2 * h6_m4 + fs_364500_20449 * r_4 * h4_m4) + e_3 * (fs_56595_1356277 * h10_m4 - fs_13456_671099 * r_2 * h8_m4 - f_32_187 * r_4 * h6_m4 - fs_720_20449 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m3, ph6_m3, ph6_m2, ph8_m3, ph8_m2, ph10_m3, ph10_m2, ab_2, pc_55, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m2 = ph10_m2[k];

        pc_55[k] = fs_16875_16 * e_0 * h4_m3 + e_1 * (-f_135_22 * h6_m3 - fs_151875_484 * r_2 * h4_m3) + e_2 * (fs_225_3718 * h8_m3 + f_18_11 * r_2 * h6_m3 + fs_151875_20449 * r_4 * h4_m3) + e_3 * (fs_94325_1356277 * h10_m3 - fs_450_671099 * r_2 * h8_m3 - f_18_187 * r_4 * h6_m3 - fs_300_20449 * r_6 * h4_m3);

        pc_56[k] = -fs_7875_2 * e_0 * r_2 * h2_m2 + e_1 * (f_5_4 * h6_m2 + fs_875_2 * r_4 * h2_m2) + e_2 * (-fs_75_121 * h8_m2 - f_1_3 * r_2 * h6_m2 - fs_7000_1089 * r_6 * h2_m2) + e_3 * (fs_1697850_17631601 * h10_m2 + fs_300_43681 * r_2 * h8_m2 + f_1_51 * r_4 * h6_m2 + fs_1750_184041 * r_8 * h2_m2) + fs_385875_128 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph2_0, ph4_m1, ph4_0, ph6_m1, ph6_0, ph8_m1, ph8_0, ph10_m1, ph10_0, ab_2, pc_57, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h10_0 = ph10_0[k];

        pc_57[k] = e_0 * (-fs_23625_32 * h4_m1 - fs_39375_4 * r_2 * h2_m1) + e_1 * (f_80_11 * h6_m1 + fs_212625_968 * r_2 * h4_m1 + fs_4375_4 * r_4 * h2_m1) + e_2 * (-fs_46389_20449 * h8_m1 - f_64_33 * r_2 * h6_m1 - fs_212625_40898 * r_4 * h4_m1 - fs_17500_1089 * r_6 * h2_m1) + e_3 * (fs_2037420_17631601 * h10_m1 + fs_185556_7382089 * r_2 * h8_m1 + f_64_561 * r_4 * h6_m1 + fs_210_20449 * r_6 * h4_m1 + fs_4375_184041 * r_8 * h2_m1) + fs_1929375_256 * e_4 * h2_m1;

        pc_58[k] = e_0 * (-f_75_2 * h4_0 - f_225_2 * r_2 * h2_0) + e_1 * (f_105_11 * h6_0 + f_225_11 * r_2 * h4_0 + f_75_2 * r_4 * h2_0) + e_2 * (-f_252_143 * h8_0 - f_28_11 * r_2 * h6_0 - f_450_143 * r_4 * h4_0 - f_50_11 * r_6 * h2_0) + e_3 * (f_1470_4199 * h10_0 + f_504_2717 * r_2 * h8_0 + f_28_187 * r_4 * h6_0 + f_20_143 * r_6 * h4_0 + f_25_143 * r_8 * h2_0) + f_1575_16 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph10_p1, ph10_p2, ab_2, pc_59, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p2 = ph10_p2[k];

        pc_59[k] = e_0 * (-fs_23625_32 * h4_p1 - fs_39375_4 * r_2 * h2_p1) + e_1 * (f_80_11 * h6_p1 + fs_212625_968 * r_2 * h4_p1 + fs_4375_4 * r_4 * h2_p1) + e_2 * (-fs_46389_20449 * h8_p1 - f_64_33 * r_2 * h6_p1 - fs_212625_40898 * r_4 * h4_p1 - fs_17500_1089 * r_6 * h2_p1) + e_3 * (fs_2037420_17631601 * h10_p1 + fs_185556_7382089 * r_2 * h8_p1 + f_64_561 * r_4 * h6_p1 + fs_210_20449 * r_6 * h4_p1 + fs_4375_184041 * r_8 * h2_p1) + fs_1929375_256 * e_4 * h2_p1;

        pc_60[k] = -fs_7875_2 * e_0 * r_2 * h2_p2 + e_1 * (f_5_4 * h6_p2 + fs_875_2 * r_4 * h2_p2) + e_2 * (-fs_75_121 * h8_p2 - f_1_3 * r_2 * h6_p2 - fs_7000_1089 * r_6 * h2_p2) + e_3 * (fs_1697850_17631601 * h10_p2 + fs_300_43681 * r_2 * h8_p2 + f_1_51 * r_4 * h6_p2 + fs_1750_184041 * r_8 * h2_p2) + fs_385875_128 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph6_p5, ph8_p3, ph8_p4, ph8_p5, ph10_p3, ph10_p4, ph10_p5, ab_2, pc_61, pc_62, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];

        pc_61[k] = fs_16875_16 * e_0 * h4_p3 + e_1 * (-f_135_22 * h6_p3 - fs_151875_484 * r_2 * h4_p3) + e_2 * (fs_225_3718 * h8_p3 + f_18_11 * r_2 * h6_p3 + fs_151875_20449 * r_4 * h4_p3) + e_3 * (fs_94325_1356277 * h10_p3 - fs_450_671099 * r_2 * h8_p3 - f_18_187 * r_4 * h6_p3 - fs_300_20449 * r_6 * h4_p3);

        pc_62[k] = fs_10125_4 * e_0 * h4_p4 + e_1 * (-f_120_11 * h6_p4 - fs_91125_121 * r_2 * h4_p4) + e_2 * (fs_3364_1859 * h8_p4 + f_32_11 * r_2 * h6_p4 + fs_364500_20449 * r_4 * h4_p4) + e_3 * (fs_56595_1356277 * h10_p4 - fs_13456_671099 * r_2 * h8_p4 - f_32_187 * r_4 * h6_p4 - fs_720_20449 * r_6 * h4_p4);

        pc_63[k] = -f_15_2 * e_1 * h6_p5 + e_2 * (fs_121_26 * h8_p5 + f_2 * r_2 * h6_p5) + e_3 * (fs_25725_1356277 * h10_p5 - fs_242_4693 * r_2 * h8_p5 - f_2_17 * r_4 * h6_p5);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph6_m5, ph6_p6, ph8_m7, ph8_m5, ph8_p6, ph10_m7, ph10_m5, ph10_p6, ab_2, pc_64, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h6_m5 = ph6_m5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_p6 = ph10_p6[k];

        pc_64[k] = f_45_4 * e_1 * h6_p6 + e_2 * (fs_63_13 * h8_p6 - f_3 * r_2 * h6_p6) + e_3 * (fs_6860_1356277 * h10_p6 - fs_252_4693 * r_2 * h8_p6 + f_3_17 * r_4 * h6_p6);

        pc_65[k] = -fs_3375_32 * e_1 * h6_m5 + e_2 * (fs_189_52 * h8_m7 - fs_15_13 * h8_m5 + fs_15_2 * r_2 * h6_m5) + e_3 * (fs_686_79781 * h10_m7 - fs_1715_2712554 * h10_m5 - fs_189_4693 * r_2 * h8_m7 + fs_60_4693 * r_2 * h8_m5 - fs_15_578 * r_4 * h6_m5);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ab_2, pc_66 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];

        pc_66[k] = fs_22275_16 * e_0 * h4_m4 + e_1 * (-fs_3375_32 * h6_m6 - fs_1125_176 * h6_m4 - fs_18225_44 * r_2 * h4_m4) + e_2 * (fs_105_104 * h8_m6 - fs_1445_676 * h8_m4 + fs_15_2 * r_2 * h6_m6 + fs_5_11 * r_2 * h6_m4 + fs_18225_1859 * r_4 * h4_m4) + e_3 * (fs_32928_1356277 * h10_m6 - fs_4116_1356277 * h10_m4 - fs_105_9386 * r_2 * h8_m6 + fs_1445_61009 * r_2 * h8_m4 - fs_15_578 * r_4 * h6_m6 - fs_5_3179 * r_4 * h6_m4 - fs_36_1859 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_67[k] = f_45_8 * e_0 * h4_m3 + e_1 * (-fs_1125_176 * h6_m5 + fs_33075_1936 * h6_m3 - f_135_44 * r_2 * h4_m3) + e_2 * (-fs_5_1144 * h8_m5 - fs_30603_14872 * h8_m3 + fs_5_11 * r_2 * h6_m5 - fs_147_121 * r_2 * h6_m3 + f_135_286 * r_4 * h4_m3) + e_3 * (fs_56595_1356277 * h10_m5 - fs_11319_1356277 * h10_m3 + fs_5_103246 * r_2 * h8_m5 + fs_30603_1342198 * r_2 * h8_m3 - fs_5_3179 * r_4 * h6_m5 + fs_147_34969 * r_4 * h6_m3 - f_3_143 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_68[k] = e_0 * (-fs_3375_4 * h4_m4 - fs_23625_64 * h4_m2 - fs_14175_4 * r_2 * h2_m2) + e_1 * (fs_33075_1936 * h6_m4 + fs_190125_3872 * h6_m2 + fs_30375_121 * r_2 * h4_m4 + fs_212625_1936 * r_2 * h4_m2 + fs_1575_4 * r_4 * h2_m2) + e_2 * (-fs_5547_7436 * h8_m4 - fs_184815_163592 * h8_m2 - fs_147_121 * r_2 * h6_m4 - fs_845_242 * r_2 * h6_m2 - fs_121500_20449 * r_4 * h4_m4 - fs_212625_81796 * r_4 * h4_m2 - fs_700_121 * r_6 * h2_m2) + e_3 * (fs_75460_1356277 * h10_m4 - fs_301840_17631601 * h10_m2 + fs_5547_671099 * r_2 * h8_m4 + fs_184815_14764178 * r_2 * h8_m2 + fs_147_34969 * r_4 * h6_m4 + fs_845_69938 * r_4 * h6_m2 + fs_240_20449 * r_6 * h4_m4 + fs_105_20449 * r_6 * h4_m2 + fs_175_20449 * r_8 * h2_m2) + fs_694575_256 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_69[k] = e_0 * (-fs_30375_32 * h4_m3 - fs_23625_32 * h4_m1 - fs_6300 * r_2 * h2_m1) + e_1 * (fs_190125_3872 * h6_m3 + f_265_44 * h6_m1 + fs_273375_968 * r_2 * h4_m3 + fs_212625_968 * r_2 * h4_m1 + fs_700 * r_4 * h2_m1) + e_2 * (-fs_2880_1859 * h8_m3 - fs_17661_81796 * h8_m1 - fs_845_242 * r_2 * h6_m3 - f_53_33 * r_2 * h6_m1 - fs_273375_40898 * r_4 * h4_m3 - fs_212625_40898 * r_4 * h4_m1 - fs_11200_1089 * r_6 * h2_m1) + e_3 * (fs_169785_2712554 * h10_m3 - fs_509355_17631601 * h10_m1 + fs_11520_671099 * r_2 * h8_m3 + fs_17661_7382089 * r_2 * h8_m1 + fs_845_69938 * r_4 * h6_m3 + f_53_561 * r_4 * h6_m1 + fs_270_20449 * r_6 * h4_m3 + fs_210_20449 * r_6 * h4_m1 + fs_2800_184041 * r_8 * h2_m1) + fs_77175_16 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph10_m2, ph10_p1, ab_2, pc_70, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_p1 = ph10_p1[k];

        pc_70[k] = e_0 * (-fs_42525_128 * h4_m2 + fs_7875_8 * r_2 * h2_m2) + e_1 * (f_265_44 * h6_m2 + fs_382725_3872 * r_2 * h4_m2 - fs_875_8 * r_4 * h2_m2) + e_2 * (-fs_118803_81796 * h8_m2 - f_53_33 * r_2 * h6_m2 - fs_382725_163592 * r_4 * h4_m2 + fs_1750_1089 * r_6 * h2_m2) + e_3 * (fs_1086624_17631601 * h10_m2 + fs_118803_7382089 * r_2 * h8_m2 + f_53_561 * r_4 * h6_m2 + fs_189_40898 * r_6 * h4_m2 - fs_875_368082 * r_8 * h2_m2) - fs_385875_512 * e_4 * h2_m2;

        pc_71[k] = e_0 * (f_15_8 * h4_p1 + fs_6750 * r_2 * h2_p1) + e_1 * (fs_2625_242 * h6_p1 - f_45_44 * r_2 * h4_p1 - fs_750 * r_4 * h2_p1) + e_2 * (-fs_55125_40898 * h8_p1 - fs_280_363 * r_2 * h6_p1 + f_45_286 * r_4 * h4_p1 + fs_4000_363 * r_6 * h2_p1) + e_3 * (fs_1901592_17631601 * h10_p1 + fs_110250_7382089 * r_2 * h8_p1 + fs_280_104907 * r_4 * h6_p1 - f_1_143 * r_6 * h4_p1 - fs_1000_61347 * r_8 * h2_p1) - fs_165375_32 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p2 = ph10_p2[k];

        pc_72[k] = e_0 * (-fs_23625_32 * h4_0 - fs_42525_128 * h4_p2 - fs_23625_2 * r_2 * h2_0 + fs_7875_8 * r_2 * h2_p2) + e_1 * (fs_2625_242 * h6_0 + f_265_44 * h6_p2 + fs_212625_968 * r_2 * h4_0 + fs_382725_3872 * r_2 * h4_p2 + fs_2625_2 * r_4 * h2_0 - fs_875_8 * r_4 * h2_p2) + e_2 * (fs_1890_20449 * h8_0 - fs_118803_81796 * h8_p2 - fs_280_363 * r_2 * h6_0 - f_53_33 * r_2 * h6_p2 - fs_212625_40898 * r_4 * h4_0 - fs_382725_163592 * r_4 * h4_p2 - fs_7000_363 * r_6 * h2_0 + fs_1750_1089 * r_6 * h2_p2) + e_3 * (-fs_1481760_17631601 * h10_0 + fs_1086624_17631601 * h10_p2 - fs_7560_7382089 * r_2 * h8_0 + fs_118803_7382089 * r_2 * h8_p2 + fs_280_104907 * r_4 * h6_0 + f_53_561 * r_4 * h6_p2 + fs_210_20449 * r_6 * h4_0 + fs_189_40898 * r_6 * h4_p2 + fs_1750_61347 * r_8 * h2_0 - fs_875_368082 * r_8 * h2_p2) + e_4 * (fs_1157625_128 * h2_0 - fs_385875_512 * h2_p2);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_73[k] = e_0 * (-fs_23625_32 * h4_p1 - fs_30375_32 * h4_p3 - fs_6300 * r_2 * h2_p1) + e_1 * (f_265_44 * h6_p1 + fs_190125_3872 * h6_p3 + fs_212625_968 * r_2 * h4_p1 + fs_273375_968 * r_2 * h4_p3 + fs_700 * r_4 * h2_p1) + e_2 * (-fs_17661_81796 * h8_p1 - fs_2880_1859 * h8_p3 - f_53_33 * r_2 * h6_p1 - fs_845_242 * r_2 * h6_p3 - fs_212625_40898 * r_4 * h4_p1 - fs_273375_40898 * r_4 * h4_p3 - fs_11200_1089 * r_6 * h2_p1) + e_3 * (-fs_509355_17631601 * h10_p1 + fs_169785_2712554 * h10_p3 + fs_17661_7382089 * r_2 * h8_p1 + fs_11520_671099 * r_2 * h8_p3 + f_53_561 * r_4 * h6_p1 + fs_845_69938 * r_4 * h6_p3 + fs_210_20449 * r_6 * h4_p1 + fs_270_20449 * r_6 * h4_p3 + fs_2800_184041 * r_8 * h2_p1) + fs_77175_16 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_74[k] = e_0 * (-fs_23625_64 * h4_p2 - fs_3375_4 * h4_p4 - fs_14175_4 * r_2 * h2_p2) + e_1 * (fs_190125_3872 * h6_p2 + fs_33075_1936 * h6_p4 + fs_212625_1936 * r_2 * h4_p2 + fs_30375_121 * r_2 * h4_p4 + fs_1575_4 * r_4 * h2_p2) + e_2 * (-fs_184815_163592 * h8_p2 - fs_5547_7436 * h8_p4 - fs_845_242 * r_2 * h6_p2 - fs_147_121 * r_2 * h6_p4 - fs_212625_81796 * r_4 * h4_p2 - fs_121500_20449 * r_4 * h4_p4 - fs_700_121 * r_6 * h2_p2) + e_3 * (-fs_301840_17631601 * h10_p2 + fs_75460_1356277 * h10_p4 + fs_184815_14764178 * r_2 * h8_p2 + fs_5547_671099 * r_2 * h8_p4 + fs_845_69938 * r_4 * h6_p2 + fs_147_34969 * r_4 * h6_p4 + fs_105_20449 * r_6 * h4_p2 + fs_240_20449 * r_6 * h4_p4 + fs_175_20449 * r_8 * h2_p2) + fs_694575_256 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_75[k] = f_45_8 * e_0 * h4_p3 + e_1 * (fs_33075_1936 * h6_p3 - fs_1125_176 * h6_p5 - f_135_44 * r_2 * h4_p3) + e_2 * (-fs_30603_14872 * h8_p3 - fs_5_1144 * h8_p5 - fs_147_121 * r_2 * h6_p3 + fs_5_11 * r_2 * h6_p5 + f_135_286 * r_4 * h4_p3) + e_3 * (-fs_11319_1356277 * h10_p3 + fs_56595_1356277 * h10_p5 + fs_30603_1342198 * r_2 * h8_p3 + fs_5_103246 * r_2 * h8_p5 + fs_147_34969 * r_4 * h6_p3 - fs_5_3179 * r_4 * h6_p5 - f_3_143 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p4, ph6_p4, ph6_p5, ph6_p6, ph8_p4, ph8_p5, ph8_p6, ph8_p7, ph10_p4, ph10_p5, ph10_p6, ph10_p7, ab_2, pc_76, pc_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h10_p7 = ph10_p7[k];

        pc_76[k] = fs_22275_16 * e_0 * h4_p4 + e_1 * (-fs_1125_176 * h6_p4 - fs_3375_32 * h6_p6 - fs_18225_44 * r_2 * h4_p4) + e_2 * (-fs_1445_676 * h8_p4 + fs_105_104 * h8_p6 + fs_5_11 * r_2 * h6_p4 + fs_15_2 * r_2 * h6_p6 + fs_18225_1859 * r_4 * h4_p4) + e_3 * (-fs_4116_1356277 * h10_p4 + fs_32928_1356277 * h10_p6 + fs_1445_61009 * r_2 * h8_p4 - fs_105_9386 * r_2 * h8_p6 - fs_5_3179 * r_4 * h6_p4 - fs_15_578 * r_4 * h6_p6 - fs_36_1859 * r_6 * h4_p4);

        pc_77[k] = -fs_3375_32 * e_1 * h6_p5 + e_2 * (-fs_15_13 * h8_p5 + fs_189_52 * h8_p7 + fs_15_2 * r_2 * h6_p5) + e_3 * (-fs_1715_2712554 * h10_p5 + fs_686_79781 * h10_p7 + fs_60_4693 * r_2 * h8_p5 - fs_189_4693 * r_2 * h8_p7 - fs_15_578 * r_4 * h6_p5);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph6_m4, ph6_m3, ph8_m8, ph8_m7, ph8_m4, ph8_m3, ph10_m8, ph10_m7, ph10_m4, ph10_m3, ab_2, pc_78, pc_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m3 = ph10_m3[k];

        pc_78[k] = fs_7425_8 * e_0 * h4_m4 + e_1 * (fs_30375_352 * h6_m4 - fs_6075_22 * r_2 * h4_m4) + e_2 * (fs_42_13 * h8_m8 + fs_135_338 * h8_m4 - fs_135_22 * r_2 * h6_m4 + fs_12150_1859 * r_4 * h4_m4) + e_3 * (fs_2058_79781 * h10_m8 + fs_343_2712554 * h10_m4 - fs_168_4693 * r_2 * h8_m8 - fs_270_61009 * r_2 * h8_m4 + fs_135_6358 * r_4 * h6_m4 - fs_24_1859 * r_6 * h4_m4);

        pc_79[k] = -fs_22275_64 * e_0 * h4_m3 + e_1 * (fs_675_11 * h6_m3 + fs_18225_176 * r_2 * h4_m3) + e_2 * (-fs_7_104 * h8_m7 + fs_1587_1352 * h8_m3 - fs_48_11 * r_2 * h6_m3 - fs_18225_7436 * r_4 * h4_m3) + e_3 * (fs_4116_79781 * h10_m7 + fs_1029_1356277 * h10_m3 + fs_7_9386 * r_2 * h8_m7 - fs_1587_122018 * r_2 * h8_m3 + fs_48_3179 * r_4 * h6_m3 + fs_9_1859 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_80[k] = e_0 * (-fs_14175_16 * h4_m2 - fs_23625_4 * r_2 * h2_m2) + e_1 * (fs_30375_352 * h6_m6 + fs_39675_3872 * h6_m2 + fs_127575_484 * r_2 * h4_m2 + fs_2625_4 * r_4 * h2_m2) + e_2 * (-fs_210_143 * h8_m6 + fs_38642_20449 * h8_m2 - fs_135_22 * r_2 * h6_m6 - fs_529_726 * r_2 * h6_m2 - fs_127575_20449 * r_4 * h4_m2 - fs_3500_363 * r_6 * h2_m2) + e_3 * (fs_90552_1356277 * h10_m6 + fs_45276_17631601 * h10_m2 + fs_840_51623 * r_2 * h8_m6 - fs_154568_7382089 * r_2 * h8_m2 + fs_135_6358 * r_4 * h6_m6 + fs_529_209814 * r_4 * h6_m2 + fs_252_20449 * r_6 * h4_m2 + fs_875_61347 * r_8 * h2_m2) + fs_1157625_256 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m1, ph10_m5, ph10_m4, ph10_m1, ab_2, pc_81, pc_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m1 = ph10_m1[k];

        pc_81[k] = e_0 * (-fs_23625_64 * h4_m1 - fs_14175_2 * r_2 * h2_m1) + e_1 * (fs_675_11 * h6_m5 - fs_450_121 * h6_m1 + fs_212625_1936 * r_2 * h4_m1 + fs_1575_2 * r_4 * h2_m1) + e_2 * (-fs_2187_1144 * h8_m5 + fs_338709_163592 * h8_m1 - fs_48_11 * r_2 * h6_m5 + fs_32_121 * r_2 * h6_m1 - fs_212625_81796 * r_4 * h4_m1 - fs_1400_121 * r_6 * h2_m1) + e_3 * (fs_94325_1356277 * h10_m5 + fs_113190_17631601 * h10_m1 + fs_2187_103246 * r_2 * h8_m5 - fs_338709_14764178 * r_2 * h8_m1 + fs_48_3179 * r_4 * h6_m5 - fs_32_34969 * r_4 * h6_m1 + fs_105_20449 * r_6 * h4_m1 + fs_350_20449 * r_8 * h2_m1) + fs_694575_128 * e_4 * h2_m1;

        pc_82[k] = fs_3375_8 * e_0 * h4_m4 + e_1 * (fs_39675_3872 * h6_m4 - fs_30375_242 * r_2 * h4_m4) + e_2 * (-fs_4107_3718 * h8_m4 - fs_529_726 * r_2 * h6_m4 + fs_60750_20449 * r_4 * h4_m4) + e_3 * (fs_169785_2712554 * h10_m4 + fs_8214_671099 * r_2 * h8_m4 + fs_529_209814 * r_4 * h6_m4 - fs_120_20449 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_83[k] = e_0 * (fs_114075_128 * h4_m3 + fs_42525_128 * h4_m1 - fs_7875_4 * r_2 * h2_m1) + e_1 * (-fs_450_121 * h6_m3 - fs_6125_121 * h6_m1 - fs_1026675_3872 * r_2 * h4_m3 - fs_382725_3872 * r_2 * h4_m1 + fs_875_4 * r_4 * h2_m1) + e_2 * (-fs_9_44 * h8_m3 + fs_55545_81796 * h8_m1 + fs_32_121 * r_2 * h6_m3 + fs_3920_1089 * r_2 * h6_m1 + fs_6075_968 * r_4 * h4_m3 + fs_382725_163592 * r_4 * h4_m1 - fs_3500_1089 * r_6 * h2_m1) + e_3 * (fs_67914_1356277 * h10_m3 + fs_407484_17631601 * h10_m1 + fs_9_3971 * r_2 * h8_m3 - fs_55545_7382089 * r_2 * h8_m1 - fs_32_34969 * r_4 * h6_m3 - fs_3920_314721 * r_4 * h6_m1 - fs_3_242 * r_6 * h4_m3 - fs_189_40898 * r_6 * h4_m1 + fs_875_184041 * r_8 * h2_m1) + fs_385875_256 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ab_2, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];

        pc_84[k] = e_0 * (f_165_4 * h4_p2 - fs_3375_4 * r_2 * h2_p2) + e_1 * (-fs_525_8 * h6_p2 - f_45_2 * r_2 * h4_p2 + fs_375_4 * r_4 * h2_p2) + e_2 * (fs_2016_20449 * h8_p2 + fs_14_3 * r_2 * h6_p2 + f_45_13 * r_4 * h4_p2 - fs_500_363 * r_6 * h2_p2) + e_3 * (fs_1267728_17631601 * h10_p2 - fs_8064_7382089 * r_2 * h8_p2 - fs_14_867 * r_4 * h6_p2 - f_2_13 * r_6 * h4_p2 + fs_125_61347 * r_8 * h2_p2) + fs_165375_256 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_85[k] = e_0 * (-fs_42525_128 * h4_p1 + fs_114075_128 * h4_p3 + fs_7875_4 * r_2 * h2_p1) + e_1 * (fs_6125_121 * h6_p1 - fs_450_121 * h6_p3 + fs_382725_3872 * r_2 * h4_p1 - fs_1026675_3872 * r_2 * h4_p3 - fs_875_4 * r_4 * h2_p1) + e_2 * (-fs_55545_81796 * h8_p1 - fs_9_44 * h8_p3 - fs_3920_1089 * r_2 * h6_p1 + fs_32_121 * r_2 * h6_p3 - fs_382725_163592 * r_4 * h4_p1 + fs_6075_968 * r_4 * h4_p3 + fs_3500_1089 * r_6 * h2_p1) + e_3 * (-fs_407484_17631601 * h10_p1 + fs_67914_1356277 * h10_p3 + fs_55545_7382089 * r_2 * h8_p1 + fs_9_3971 * r_2 * h8_p3 + fs_3920_314721 * r_4 * h6_p1 - fs_32_34969 * r_4 * h6_p3 + fs_189_40898 * r_6 * h4_p1 - fs_3_242 * r_6 * h4_p3 - fs_875_184041 * r_8 * h2_p1) - fs_385875_256 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_86[k] = e_0 * (fs_3375_8 * h4_p4 - fs_9450 * r_2 * h2_0) + e_1 * (-fs_525_8 * h6_0 + fs_39675_3872 * h6_p4 - fs_30375_242 * r_2 * h4_p4 + fs_1050 * r_4 * h2_0) + e_2 * (fs_378_121 * h8_0 - fs_4107_3718 * h8_p4 + fs_14_3 * r_2 * h6_0 - fs_529_726 * r_2 * h6_p4 + fs_60750_20449 * r_4 * h4_p4 - fs_5600_363 * r_6 * h2_0) + e_3 * (fs_463050_17631601 * h10_0 + fs_169785_2712554 * h10_p4 - fs_1512_43681 * r_2 * h8_0 + fs_8214_671099 * r_2 * h8_p4 - fs_14_867 * r_4 * h6_0 + fs_529_209814 * r_4 * h6_p4 - fs_120_20449 * r_6 * h4_p4 + fs_1400_61347 * r_8 * h2_0) + fs_231525_32 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_87[k] = e_0 * (-fs_23625_64 * h4_p1 - fs_14175_2 * r_2 * h2_p1) + e_1 * (-fs_450_121 * h6_p1 + fs_675_11 * h6_p5 + fs_212625_1936 * r_2 * h4_p1 + fs_1575_2 * r_4 * h2_p1) + e_2 * (fs_338709_163592 * h8_p1 - fs_2187_1144 * h8_p5 + fs_32_121 * r_2 * h6_p1 - fs_48_11 * r_2 * h6_p5 - fs_212625_81796 * r_4 * h4_p1 - fs_1400_121 * r_6 * h2_p1) + e_3 * (fs_113190_17631601 * h10_p1 + fs_94325_1356277 * h10_p5 - fs_338709_14764178 * r_2 * h8_p1 + fs_2187_103246 * r_2 * h8_p5 - fs_32_34969 * r_4 * h6_p1 + fs_48_3179 * r_4 * h6_p5 + fs_105_20449 * r_6 * h4_p1 + fs_350_20449 * r_8 * h2_p1) + fs_694575_128 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_88[k] = e_0 * (-fs_14175_16 * h4_p2 - fs_23625_4 * r_2 * h2_p2) + e_1 * (fs_39675_3872 * h6_p2 + fs_30375_352 * h6_p6 + fs_127575_484 * r_2 * h4_p2 + fs_2625_4 * r_4 * h2_p2) + e_2 * (fs_38642_20449 * h8_p2 - fs_210_143 * h8_p6 - fs_529_726 * r_2 * h6_p2 - fs_135_22 * r_2 * h6_p6 - fs_127575_20449 * r_4 * h4_p2 - fs_3500_363 * r_6 * h2_p2) + e_3 * (fs_45276_17631601 * h10_p2 + fs_90552_1356277 * h10_p6 - fs_154568_7382089 * r_2 * h8_p2 + fs_840_51623 * r_2 * h8_p6 + fs_529_209814 * r_4 * h6_p2 + fs_135_6358 * r_4 * h6_p6 + fs_252_20449 * r_6 * h4_p2 + fs_875_61347 * r_8 * h2_p2) + fs_1157625_256 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ph8_p7, ph8_p8, ph10_p3, ph10_p4, ph10_p7, ph10_p8, ab_2, pc_89, pc_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h10_p8 = ph10_p8[k];

        pc_89[k] = -fs_22275_64 * e_0 * h4_p3 + e_1 * (fs_675_11 * h6_p3 + fs_18225_176 * r_2 * h4_p3) + e_2 * (fs_1587_1352 * h8_p3 - fs_7_104 * h8_p7 - fs_48_11 * r_2 * h6_p3 - fs_18225_7436 * r_4 * h4_p3) + e_3 * (fs_1029_1356277 * h10_p3 + fs_4116_79781 * h10_p7 - fs_1587_122018 * r_2 * h8_p3 + fs_7_9386 * r_2 * h8_p7 + fs_48_3179 * r_4 * h6_p3 + fs_9_1859 * r_6 * h4_p3);

        pc_90[k] = fs_7425_8 * e_0 * h4_p4 + e_1 * (fs_30375_352 * h6_p4 - fs_6075_22 * r_2 * h4_p4) + e_2 * (fs_135_338 * h8_p4 + fs_42_13 * h8_p8 - fs_135_22 * r_2 * h6_p4 + fs_12150_1859 * r_4 * h4_p4) + e_3 * (fs_343_2712554 * h10_p4 + fs_2058_79781 * h10_p8 - fs_270_61009 * r_2 * h8_p4 - fs_168_4693 * r_2 * h8_p8 + fs_135_6358 * r_4 * h6_p4 - fs_24_1859 * r_6 * h4_p4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m3, ph4_m2, ph6_m3, ph6_m2, ph8_m8, ph8_m3, ph8_m2, ph10_m9, ph10_m8, ph10_m3, ph10_m2, ab_2, pc_91, pc_92 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m2 = ph10_m2[k];

        pc_91[k] = -fs_51975_32 * e_0 * h4_m3 + e_1 * (-fs_14175_352 * h6_m3 + fs_42525_88 * r_2 * h4_m3) + e_2 * (-fs_63_676 * h8_m3 + fs_63_22 * r_2 * h6_m3 - fs_42525_3718 * r_4 * h4_m3) + e_3 * (fs_294_4199 * h10_m9 - fs_49_2712554 * h10_m3 + fs_63_61009 * r_2 * h8_m3 - fs_63_6358 * r_4 * h6_m3 + fs_42_1859 * r_6 * h4_m3);

        pc_92[k] = e_0 * (-fs_22275_64 * h4_m2 - fs_37125_4 * r_2 * h2_m2) + e_1 * (-fs_25725_352 * h6_m2 + fs_18225_176 * r_2 * h4_m2 + fs_4125_4 * r_4 * h2_m2) + e_2 * (-fs_49_13 * h8_m8 - fs_5887_14872 * h8_m2 + fs_343_66 * r_2 * h6_m2 - fs_18225_7436 * r_4 * h4_m2 - fs_500_33 * r_6 * h2_m2) + e_3 * (fs_7056_79781 * h10_m8 - fs_2352_17631601 * h10_m2 + fs_196_4693 * r_2 * h8_m8 + fs_5887_1342198 * r_2 * h8_m2 - fs_343_19074 * r_4 * h6_m2 + fs_9_1859 * r_6 * h4_m2 + fs_125_5577 * r_8 * h2_m2) + fs_1819125_256 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m1, ph8_m7, ph8_m6, ph8_m1, ph10_m7, ph10_m6, ph10_m1, ab_2, pc_93, pc_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m1 = ph10_m1[k];

        pc_93[k] = e_0 * (f_45_8 * h4_m1 - fs_6750 * r_2 * h2_m1) + e_1 * (-fs_65625_968 * h6_m1 - f_135_44 * r_2 * h4_m1 + fs_750 * r_4 * h2_m1) + e_2 * (-fs_2401_1144 * h8_m7 - fs_153125_163592 * h8_m1 + fs_1750_363 * r_2 * h6_m1 + f_135_286 * r_4 * h4_m1 - fs_4000_363 * r_6 * h2_m1) + e_3 * (fs_6468_79781 * h10_m7 - fs_9702_17631601 * h10_m1 + fs_2401_103246 * r_2 * h8_m7 + fs_153125_14764178 * r_2 * h8_m1 - fs_1750_104907 * r_4 * h6_m1 - f_3_143 * r_6 * h4_m1 + fs_1000_61347 * r_8 * h2_m1) + fs_165375_32 * e_4 * h2_m1;

        pc_94[k] = -fs_14175_352 * e_1 * h6_m6 + e_2 * (-fs_441_1144 * h8_m6 + fs_63_22 * r_2 * h6_m6) + e_3 * (fs_86240_1356277 * h10_m6 + fs_441_103246 * r_2 * h8_m6 - fs_63_6358 * r_4 * h6_m6);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ab_2, pc_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_95[k] = e_0 * (fs_30375_32 * h4_m1 - f_30 * r_2 * h2_m1) + e_1 * (-fs_25725_352 * h6_m5 - fs_8575_1936 * h6_m1 - fs_273375_968 * r_2 * h4_m1 + f_10 * r_4 * h2_m1) + e_2 * (fs_21_572 * h8_m5 - fs_42483_20449 * h8_m1 + fs_343_66 * r_2 * h6_m5 + fs_343_1089 * r_2 * h6_m1 + fs_273375_40898 * r_4 * h4_m1 - f_40_33 * r_6 * h2_m1) + e_3 * (fs_121275_2712554 * h10_m5 - fs_72765_17631601 * h10_m1 - fs_21_51623 * r_2 * h8_m5 + fs_169932_7382089 * r_2 * h8_m1 - fs_343_19074 * r_4 * h6_m5 - fs_343_314721 * r_4 * h6_m1 - fs_270_20449 * r_6 * h4_m1 + f_20_429 * r_8 * h2_m1) + f_105_4 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_96[k] = e_0 * (-fs_4725_32 * h4_m4 - fs_114075_128 * h4_m2 + fs_1125_8 * r_2 * h2_m2) + e_1 * (-fs_65625_968 * h6_m4 - fs_8575_1936 * h6_m2 + fs_42525_968 * r_2 * h4_m4 + fs_1026675_3872 * r_2 * h4_m2 - fs_125_8 * r_4 * h2_m2) + e_2 * (fs_2625_3718 * h8_m4 + fs_1029_484 * h8_m2 + fs_1750_363 * r_2 * h6_m4 + fs_343_1089 * r_2 * h6_m2 - fs_42525_40898 * r_4 * h4_m4 - fs_6075_968 * r_4 * h4_m2 + fs_250_1089 * r_6 * h2_m2) + e_3 * (fs_38808_1356277 * h10_m4 + fs_155232_17631601 * h10_m2 - fs_5250_671099 * r_2 * h8_m4 - fs_1029_43681 * r_2 * h8_m2 - fs_1750_104907 * r_4 * h6_m4 - fs_343_314721 * r_4 * h6_m2 + fs_42_20449 * r_6 * h4_m4 + fs_3_242 * r_6 * h4_m2 - fs_125_368082 * r_8 * h2_m2) - fs_55125_512 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ab_2, pc_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];

        pc_97[k] = -f_255_8 * e_0 * h4_p3 + e_1 * (-fs_33075_484 * h6_p3 + f_765_44 * r_2 * h4_p3) + e_2 * (fs_11907_3718 * h8_p3 + fs_588_121 * r_2 * h6_p3 - f_765_286 * r_4 * h4_p3) + e_3 * (fs_45276_1356277 * h10_p3 - fs_23814_671099 * r_2 * h8_p3 - fs_588_34969 * r_4 * h6_p3 + f_17_143 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_98[k] = e_0 * (fs_114075_128 * h4_p2 - fs_4725_32 * h4_p4 - fs_1125_8 * r_2 * h2_p2) + e_1 * (fs_8575_1936 * h6_p2 - fs_65625_968 * h6_p4 - fs_1026675_3872 * r_2 * h4_p2 + fs_42525_968 * r_2 * h4_p4 + fs_125_8 * r_4 * h2_p2) + e_2 * (-fs_1029_484 * h8_p2 + fs_2625_3718 * h8_p4 - fs_343_1089 * r_2 * h6_p2 + fs_1750_363 * r_2 * h6_p4 + fs_6075_968 * r_4 * h4_p2 - fs_42525_40898 * r_4 * h4_p4 - fs_250_1089 * r_6 * h2_p2) + e_3 * (-fs_155232_17631601 * h10_p2 + fs_38808_1356277 * h10_p4 + fs_1029_43681 * r_2 * h8_p2 - fs_5250_671099 * r_2 * h8_p4 + fs_343_314721 * r_4 * h6_p2 - fs_1750_104907 * r_4 * h6_p4 - fs_3_242 * r_6 * h4_p2 + fs_42_20449 * r_6 * h4_p4 + fs_125_368082 * r_8 * h2_p2) + fs_55125_512 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_99[k] = e_0 * (-fs_30375_32 * h4_p1 + f_30 * r_2 * h2_p1) + e_1 * (fs_8575_1936 * h6_p1 - fs_25725_352 * h6_p5 + fs_273375_968 * r_2 * h4_p1 - f_10 * r_4 * h2_p1) + e_2 * (fs_42483_20449 * h8_p1 + fs_21_572 * h8_p5 - fs_343_1089 * r_2 * h6_p1 + fs_343_66 * r_2 * h6_p5 - fs_273375_40898 * r_4 * h4_p1 + f_40_33 * r_6 * h2_p1) + e_3 * (fs_72765_17631601 * h10_p1 + fs_121275_2712554 * h10_p5 - fs_169932_7382089 * r_2 * h8_p1 - fs_21_51623 * r_2 * h8_p5 + fs_343_314721 * r_4 * h6_p1 - fs_343_19074 * r_4 * h6_p5 + fs_270_20449 * r_6 * h4_p1 - f_20_429 * r_8 * h2_p1) - f_105_4 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_100[k] = e_0 * (fs_16875_16 * h4_0 - fs_6075 * r_2 * h2_0) + e_1 * (-fs_33075_484 * h6_0 - fs_14175_352 * h6_p6 - fs_151875_484 * r_2 * h4_0 + fs_675 * r_4 * h2_0) + e_2 * (-fs_64827_20449 * h8_0 - fs_441_1144 * h8_p6 + fs_588_121 * r_2 * h6_0 + fs_63_22 * r_2 * h6_p6 + fs_151875_20449 * r_4 * h4_0 - fs_1200_121 * r_6 * h2_0) + e_3 * (-fs_58800_17631601 * h10_0 + fs_86240_1356277 * h10_p6 + fs_259308_7382089 * r_2 * h8_0 + fs_441_103246 * r_2 * h8_p6 - fs_588_34969 * r_4 * h6_0 - fs_63_6358 * r_4 * h6_p6 - fs_300_20449 * r_6 * h4_0 + fs_300_20449 * r_8 * h2_0) + fs_297675_64 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_101[k] = e_0 * (f_45_8 * h4_p1 - fs_6750 * r_2 * h2_p1) + e_1 * (-fs_65625_968 * h6_p1 - f_135_44 * r_2 * h4_p1 + fs_750 * r_4 * h2_p1) + e_2 * (-fs_153125_163592 * h8_p1 - fs_2401_1144 * h8_p7 + fs_1750_363 * r_2 * h6_p1 + f_135_286 * r_4 * h4_p1 - fs_4000_363 * r_6 * h2_p1) + e_3 * (-fs_9702_17631601 * h10_p1 + fs_6468_79781 * h10_p7 + fs_153125_14764178 * r_2 * h8_p1 + fs_2401_103246 * r_2 * h8_p7 - fs_1750_104907 * r_4 * h6_p1 - f_3_143 * r_6 * h4_p1 + fs_1000_61347 * r_8 * h2_p1) + fs_165375_32 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph8_p8, ph10_p2, ph10_p3, ph10_p8, ph10_p9, ab_2, pc_102, pc_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h10_p9 = ph10_p9[k];

        pc_102[k] = e_0 * (-fs_22275_64 * h4_p2 - fs_37125_4 * r_2 * h2_p2) + e_1 * (-fs_25725_352 * h6_p2 + fs_18225_176 * r_2 * h4_p2 + fs_4125_4 * r_4 * h2_p2) + e_2 * (-fs_5887_14872 * h8_p2 - fs_49_13 * h8_p8 + fs_343_66 * r_2 * h6_p2 - fs_18225_7436 * r_4 * h4_p2 - fs_500_33 * r_6 * h2_p2) + e_3 * (-fs_2352_17631601 * h10_p2 + fs_7056_79781 * h10_p8 + fs_5887_1342198 * r_2 * h8_p2 + fs_196_4693 * r_2 * h8_p8 - fs_343_19074 * r_4 * h6_p2 + fs_9_1859 * r_6 * h4_p2 + fs_125_5577 * r_8 * h2_p2) + fs_1819125_256 * e_4 * h2_p2;

        pc_103[k] = -fs_51975_32 * e_0 * h4_p3 + e_1 * (-fs_14175_352 * h6_p3 + fs_42525_88 * r_2 * h4_p3) + e_2 * (-fs_63_676 * h8_p3 + fs_63_22 * r_2 * h6_p3 - fs_42525_3718 * r_4 * h4_p3) + e_3 * (-fs_49_2712554 * h10_p3 + fs_294_4199 * h10_p9 + fs_63_61009 * r_2 * h8_p3 - fs_63_6358 * r_4 * h6_p3 + fs_42_1859 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m10, ph10_m9, ph10_m2, ph10_m1, ab_2, pc_104, pc_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_m1 = ph10_m1[k];

        pc_104[k] = e_0 * (fs_7425_8 * h4_m2 - fs_111375_8 * r_2 * h2_m2) + e_1 * (fs_1575_176 * h6_m2 - fs_6075_22 * r_2 * h4_m2 + fs_12375_8 * r_4 * h2_m2) + e_2 * (fs_21_1859 * h8_m2 - fs_7_11 * r_2 * h6_m2 + fs_12150_1859 * r_4 * h4_m2 - fs_250_11 * r_6 * h2_m2) + e_3 * (fs_735_4199 * h10_m10 + fs_49_35263202 * h10_m2 - fs_84_671099 * r_2 * h8_m2 + fs_7_3179 * r_4 * h6_m2 - fs_24_1859 * r_6 * h4_m2 + fs_125_3718 * r_8 * h2_m2) + fs_5457375_512 * e_4 * h2_m2;

        pc_105[k] = e_0 * (fs_22275_16 * h4_m1 - fs_37125_8 * r_2 * h2_m1) + e_1 * (fs_2625_88 * h6_m1 - fs_18225_44 * r_2 * h4_m1 + fs_4125_8 * r_4 * h2_m1) + e_2 * (fs_245_3718 * h8_m1 - fs_70_33 * r_2 * h6_m1 + fs_18225_1859 * r_4 * h4_m1 - fs_250_33 * r_6 * h2_m1) + e_3 * (fs_441_4199 * h10_m9 + fs_441_35263202 * h10_m1 - fs_490_671099 * r_2 * h8_m1 + fs_70_9537 * r_4 * h6_m1 - fs_36_1859 * r_6 * h4_m1 + fs_125_11154 * r_8 * h2_m1) + fs_1819125_512 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_106, pc_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_106[k] = fs_196_143 * e_2 * h8_m8 + e_3 * (fs_4851_79781 * h10_m8 - fs_784_51623 * r_2 * h8_m8);

        pc_107[k] = e_0 * (fs_3375_4 * h4_m1 - fs_2025_8 * r_2 * h2_m1) + e_1 * (fs_77175_968 * h6_m1 - fs_30375_121 * r_2 * h4_m1 + fs_225_8 * r_4 * h2_m1) + e_2 * (fs_735_286 * h8_m7 + fs_10584_20449 * h8_m1 - fs_686_121 * r_2 * h6_m1 + fs_121500_20449 * r_4 * h4_m1 - fs_50_121 * r_6 * h2_m1) + e_3 * (fs_2695_79781 * h10_m7 + fs_8085_35263202 * h10_m1 - fs_1470_51623 * r_2 * h8_m7 - fs_42336_7382089 * r_2 * h8_m1 + fs_686_34969 * r_4 * h6_m1 - fs_240_20449 * r_6 * h4_m1 + fs_25_40898 * r_8 * h2_m1) + fs_99225_512 * e_4 * h2_m1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_108[k] = e_0 * (-fs_3375_8 * h4_m2 + fs_225_8 * r_2 * h2_m2) + e_1 * (fs_1575_176 * h6_m6 - fs_42875_484 * h6_m2 + fs_30375_242 * r_2 * h4_m2 - fs_25_8 * r_4 * h2_m2) + e_2 * (fs_441_143 * h8_m6 - fs_20580_20449 * h8_m2 - fs_7_11 * r_2 * h6_m6 + fs_6860_1089 * r_2 * h6_m2 - fs_60750_20449 * r_4 * h4_m2 + fs_50_1089 * r_6 * h2_m2) + e_3 * (fs_24255_1356277 * h10_m6 - fs_24255_35263202 * h10_m2 - fs_1764_51623 * r_2 * h8_m6 + fs_82320_7382089 * r_2 * h8_m2 + fs_7_3179 * r_4 * h6_m6 - fs_6860_314721 * r_4 * h6_m2 + fs_120_20449 * r_6 * h4_m2 - fs_25_368082 * r_8 * h2_m2) - fs_11025_512 * e_4 * h2_m2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_p4, ph6_m5, ph6_m3, ph6_p4, ph8_m5, ph8_m3, ph8_p4, ph10_m5, ph10_m3, ph10_p4, ab_2, pc_109, pc_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_p4 = ph10_p4[k];

        pc_109[k] = fs_4725_32 * e_0 * h4_m3 + e_1 * (fs_2625_88 * h6_m5 + fs_77175_968 * h6_m3 - fs_42525_968 * r_2 * h4_m3) + e_2 * (fs_420_143 * h8_m5 + fs_3087_1859 * h8_m3 - fs_70_33 * r_2 * h6_m5 - fs_686_121 * r_2 * h6_m3 + fs_42525_40898 * r_4 * h4_m3) + e_3 * (fs_24255_2712554 * h10_m5 + fs_4851_2712554 * h10_m3 - fs_1680_51623 * r_2 * h8_m5 - fs_12348_671099 * r_2 * h8_m3 + fs_70_9537 * r_4 * h6_m5 + fs_686_34969 * r_4 * h6_m3 - fs_42_20449 * r_6 * h4_m3);

        pc_110[k] = f_15_2 * e_0 * h4_p4 + e_1 * (fs_55125_484 * h6_p4 - f_45_11 * r_2 * h4_p4) + e_2 * (fs_8820_1859 * h8_p4 - fs_980_121 * r_2 * h6_p4 + f_90_143 * r_4 * h4_p4) + e_3 * (fs_11319_1356277 * h10_p4 - fs_35280_671099 * r_2 * h8_p4 + fs_980_34969 * r_4 * h6_p4 - f_4_143 * r_6 * h4_p4);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_111[k] = -fs_4725_32 * e_0 * h4_p3 + e_1 * (-fs_77175_968 * h6_p3 + fs_2625_88 * h6_p5 + fs_42525_968 * r_2 * h4_p3) + e_2 * (-fs_3087_1859 * h8_p3 + fs_420_143 * h8_p5 + fs_686_121 * r_2 * h6_p3 - fs_70_33 * r_2 * h6_p5 - fs_42525_40898 * r_4 * h4_p3) + e_3 * (-fs_4851_2712554 * h10_p3 + fs_24255_2712554 * h10_p5 + fs_12348_671099 * r_2 * h8_p3 - fs_1680_51623 * r_2 * h8_p5 - fs_686_34969 * r_4 * h6_p3 + fs_70_9537 * r_4 * h6_p5 + fs_42_20449 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_112[k] = e_0 * (fs_3375_8 * h4_p2 - fs_225_8 * r_2 * h2_p2) + e_1 * (fs_42875_484 * h6_p2 + fs_1575_176 * h6_p6 - fs_30375_242 * r_2 * h4_p2 + fs_25_8 * r_4 * h2_p2) + e_2 * (fs_20580_20449 * h8_p2 + fs_441_143 * h8_p6 - fs_6860_1089 * r_2 * h6_p2 - fs_7_11 * r_2 * h6_p6 + fs_60750_20449 * r_4 * h4_p2 - fs_50_1089 * r_6 * h2_p2) + e_3 * (fs_24255_35263202 * h10_p2 + fs_24255_1356277 * h10_p6 - fs_82320_7382089 * r_2 * h8_p2 - fs_1764_51623 * r_2 * h8_p6 + fs_6860_314721 * r_4 * h6_p2 + fs_7_3179 * r_4 * h6_p6 - fs_120_20449 * r_6 * h4_p2 + fs_25_368082 * r_8 * h2_p2) + fs_11025_512 * e_4 * h2_p2;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_113[k] = e_0 * (-fs_3375_4 * h4_p1 + fs_2025_8 * r_2 * h2_p1) + e_1 * (-fs_77175_968 * h6_p1 + fs_30375_121 * r_2 * h4_p1 - fs_225_8 * r_4 * h2_p1) + e_2 * (-fs_10584_20449 * h8_p1 + fs_735_286 * h8_p7 + fs_686_121 * r_2 * h6_p1 - fs_121500_20449 * r_4 * h4_p1 + fs_50_121 * r_6 * h2_p1) + e_3 * (-fs_8085_35263202 * h10_p1 + fs_2695_79781 * h10_p7 + fs_42336_7382089 * r_2 * h8_p1 - fs_1470_51623 * r_2 * h8_p7 - fs_686_34969 * r_4 * h6_p1 + fs_240_20449 * r_6 * h4_p1 - fs_25_40898 * r_8 * h2_p1) - fs_99225_512 * e_4 * h2_p1;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_114 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_114[k] = e_0 * (fs_10125_4 * h4_0 - fs_10125_4 * r_2 * h2_0) + e_1 * (fs_55125_484 * h6_0 - fs_91125_121 * r_2 * h4_0 + fs_1125_4 * r_4 * h2_0) + e_2 * (fs_8820_20449 * h8_0 + fs_196_143 * h8_p8 - fs_980_121 * r_2 * h6_0 + fs_364500_20449 * r_4 * h4_0 - fs_500_121 * r_6 * h2_0) + e_3 * (fs_2205_17631601 * h10_0 + fs_4851_79781 * h10_p8 - fs_35280_7382089 * r_2 * h8_0 - fs_784_51623 * r_2 * h8_p8 + fs_980_34969 * r_4 * h6_0 - fs_720_20449 * r_6 * h4_0 + fs_125_20449 * r_8 * h2_0) + fs_496125_256 * e_4 * h2_0;
    }

    // NOTE: the rows are formed in 83 loops, as the vectorizer runs out of
    // registers with all 117 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph10_p1, ph10_p2, ph10_p9, ph10_p10, ab_2, pc_115, pc_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h10_p10 = ph10_p10[k];

        pc_115[k] = e_0 * (fs_22275_16 * h4_p1 - fs_37125_8 * r_2 * h2_p1) + e_1 * (fs_2625_88 * h6_p1 - fs_18225_44 * r_2 * h4_p1 + fs_4125_8 * r_4 * h2_p1) + e_2 * (fs_245_3718 * h8_p1 - fs_70_33 * r_2 * h6_p1 + fs_18225_1859 * r_4 * h4_p1 - fs_250_33 * r_6 * h2_p1) + e_3 * (fs_441_35263202 * h10_p1 + fs_441_4199 * h10_p9 - fs_490_671099 * r_2 * h8_p1 + fs_70_9537 * r_4 * h6_p1 - fs_36_1859 * r_6 * h4_p1 + fs_125_11154 * r_8 * h2_p1) + fs_1819125_512 * e_4 * h2_p1;

        pc_116[k] = e_0 * (fs_7425_8 * h4_p2 - fs_111375_8 * r_2 * h2_p2) + e_1 * (fs_1575_176 * h6_p2 - fs_6075_22 * r_2 * h4_p2 + fs_12375_8 * r_4 * h2_p2) + e_2 * (fs_21_1859 * h8_p2 - fs_7_11 * r_2 * h6_p2 + fs_12150_1859 * r_4 * h4_p2 - fs_250_11 * r_6 * h2_p2) + e_3 * (fs_49_35263202 * h10_p2 + fs_735_4199 * h10_p10 - fs_84_671099 * r_2 * h8_p2 + fs_7_3179 * r_4 * h6_p2 - fs_24_1859 * r_6 * h4_p2 + fs_125_3718 * r_8 * h2_p2) + fs_5457375_512 * e_4 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};

    for (size_t m = 0; m < 117; m++)
    {
        const auto *pc = buffer.data(5 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
