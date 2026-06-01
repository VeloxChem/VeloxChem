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

#include "OverlapABRecII.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_i_i(
    const CBasisFunction &bra,
    const CBasisFunction &ket,
    const TPoint<double> &bra_center,
    const TPoint<double> &ket_center,
    double *buffer
) -> void
{
    // ---- Phase 1: geometry ----
    const auto a_xyz = bra_center.coordinates();
    const auto b_xyz = ket_center.coordinates();
    const auto AB_x  = a_xyz[0] - b_xyz[0];
    const auto AB_y  = a_xyz[1] - b_xyz[1];
    const auto AB_z  = a_xyz[2] - b_xyz[2];
    const auto AB_xx = AB_x * AB_x;
    const auto AB_yy = AB_y * AB_y;
    const auto AB_zz = AB_z * AB_z;
    const auto AB_xy = AB_x * AB_y;
    const auto AB_xz = AB_x * AB_z;
    const auto AB_yz = AB_y * AB_z;
    const auto R2    = AB_xx + AB_yy + AB_zz;

    const auto sqrt2 = std::sqrt(2.0);
    const auto sqrt3 = std::sqrt(3.0);
    const auto sqrt5 = std::sqrt(5.0);
    const auto sqrt6 = std::sqrt(6.0);
    const auto sqrt7 = std::sqrt(7.0);
    const auto sqrt10 = std::sqrt(10.0);
    const auto sqrt11 = std::sqrt(11.0);
    const auto sqrt13 = std::sqrt(13.0);
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt22 = std::sqrt(22.0);
    const auto sqrt26 = std::sqrt(26.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt39 = std::sqrt(39.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt55 = std::sqrt(55.0);
    const auto sqrt65 = std::sqrt(65.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt77 = std::sqrt(77.0);
    const auto sqrt78 = std::sqrt(78.0);
    const auto sqrt91 = std::sqrt(91.0);
    const auto sqrt130 = std::sqrt(130.0);
    const auto sqrt143 = std::sqrt(143.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt182 = std::sqrt(182.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt221 = std::sqrt(221.0);
    const auto sqrt273 = std::sqrt(273.0);
    const auto sqrt286 = std::sqrt(286.0);
    const auto sqrt330 = std::sqrt(330.0);
    const auto sqrt385 = std::sqrt(385.0);
    const auto sqrt429 = std::sqrt(429.0);
    const auto sqrt442 = std::sqrt(442.0);
    const auto sqrt455 = std::sqrt(455.0);
    const auto sqrt462 = std::sqrt(462.0);
    const auto sqrt663 = std::sqrt(663.0);
    const auto sqrt715 = std::sqrt(715.0);
    const auto sqrt770 = std::sqrt(770.0);
    const auto sqrt858 = std::sqrt(858.0);
    const auto sqrt1001 = std::sqrt(1001.0);
    const auto sqrt1155 = std::sqrt(1155.0);
    const auto sqrt1430 = std::sqrt(1430.0);
    const auto sqrt2145 = std::sqrt(2145.0);
    const auto sqrt2310 = std::sqrt(2310.0);
    const auto sqrt2431 = std::sqrt(2431.0);
    const auto sqrt3003 = std::sqrt(3003.0);
    const auto sqrt3315 = std::sqrt(3315.0);
    const auto sqrt4290 = std::sqrt(4290.0);
    const auto sqrt8398 = std::sqrt(8398.0);
    const auto sqrt12155 = std::sqrt(12155.0);
    const auto sqrt12597 = std::sqrt(12597.0);
    const auto sqrt15015 = std::sqrt(15015.0);
    const auto sqrt20995 = std::sqrt(20995.0);
    const auto sqrt92378 = std::sqrt(92378.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^6 · β^6 · p^{-12} · (s|s)
    //   V[1] ↔ α^5 · β^5 · p^{-11} · (s|s)
    //   V[2] ↔ α^4 · β^4 · p^{-10} · (s|s)
    //   V[3] ↔ α^3 · β^3 · p^{-9} · (s|s)
    //   V[4] ↔ α^2 · β^2 · p^{-8} · (s|s)
    //   V[5] ↔ α · β · p^{-7} · (s|s)
    //   V[6] ↔ p^{-6} · (s|s)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 7> V = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < exps_a.size(); ++i)
    {
        const auto alpha = exps_a[i];
        const auto ca    = coefs_a[i];
        const auto alpha2 = alpha * alpha;
        const auto alpha3 = alpha2 * alpha;
        const auto alpha4 = alpha3 * alpha;
        const auto alpha5 = alpha4 * alpha;
        const auto alpha6 = alpha5 * alpha;
        for (std::size_t j = 0; j < exps_b.size(); ++j)
        {
            const auto beta = exps_b[j];
            const auto cb   = coefs_b[j];
            const auto beta2 = beta * beta;
            const auto beta3 = beta2 * beta;
            const auto beta4 = beta3 * beta;
            const auto beta5 = beta4 * beta;
            const auto beta6 = beta5 * beta;

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto pinv7 = pinv6 * pinv;
            const auto pinv8 = pinv7 * pinv;
            const auto pinv9 = pinv8 * pinv;
            const auto pinv10 = pinv9 * pinv;
            const auto pinv11 = pinv10 * pinv;
            const auto pinv12 = pinv11 * pinv;

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha6 * beta6 * pinv12;
            V[1] += cab_ss * alpha5 * beta5 * pinv11;
            V[2] += cab_ss * alpha4 * beta4 * pinv10;
            V[3] += cab_ss * alpha3 * beta3 * pinv9;
            V[4] += cab_ss * alpha2 * beta2 * pinv8;
            V[5] += cab_ss * alpha * beta * pinv7;
            V[6] += cab_ss * pinv6;
        }
    }

    // ---- Phase 3: fused M·V → 13 × 13 spherical block ----
    const auto Y2_n2 = harm::Y_ll_2_m_n2(AB_x, AB_y, AB_z);
    const auto Y2_n1 = harm::Y_ll_2_m_n1(AB_x, AB_y, AB_z);
    const auto Y2_p0 = harm::Y_ll_2_m_p0(AB_x, AB_y, AB_z);
    const auto Y2_p1 = harm::Y_ll_2_m_p1(AB_x, AB_y, AB_z);
    const auto Y2_p2 = harm::Y_ll_2_m_p2(AB_x, AB_y, AB_z);
    const auto Y4_n4 = harm::Y_ll_4_m_n4(AB_x, AB_y, AB_z);
    const auto Y4_n3 = harm::Y_ll_4_m_n3(AB_x, AB_y, AB_z);
    const auto Y4_n2 = harm::Y_ll_4_m_n2(AB_x, AB_y, AB_z);
    const auto Y4_n1 = harm::Y_ll_4_m_n1(AB_x, AB_y, AB_z);
    const auto Y4_p0 = harm::Y_ll_4_m_p0(AB_x, AB_y, AB_z);
    const auto Y4_p1 = harm::Y_ll_4_m_p1(AB_x, AB_y, AB_z);
    const auto Y4_p2 = harm::Y_ll_4_m_p2(AB_x, AB_y, AB_z);
    const auto Y4_p3 = harm::Y_ll_4_m_p3(AB_x, AB_y, AB_z);
    const auto Y4_p4 = harm::Y_ll_4_m_p4(AB_x, AB_y, AB_z);
    const auto Y6_n6 = harm::Y_ll_6_m_n6(AB_x, AB_y, AB_z);
    const auto Y6_n5 = harm::Y_ll_6_m_n5(AB_x, AB_y, AB_z);
    const auto Y6_n4 = harm::Y_ll_6_m_n4(AB_x, AB_y, AB_z);
    const auto Y6_n3 = harm::Y_ll_6_m_n3(AB_x, AB_y, AB_z);
    const auto Y6_n2 = harm::Y_ll_6_m_n2(AB_x, AB_y, AB_z);
    const auto Y6_n1 = harm::Y_ll_6_m_n1(AB_x, AB_y, AB_z);
    const auto Y6_p0 = harm::Y_ll_6_m_p0(AB_x, AB_y, AB_z);
    const auto Y6_p1 = harm::Y_ll_6_m_p1(AB_x, AB_y, AB_z);
    const auto Y6_p2 = harm::Y_ll_6_m_p2(AB_x, AB_y, AB_z);
    const auto Y6_p3 = harm::Y_ll_6_m_p3(AB_x, AB_y, AB_z);
    const auto Y6_p4 = harm::Y_ll_6_m_p4(AB_x, AB_y, AB_z);
    const auto Y6_p5 = harm::Y_ll_6_m_p5(AB_x, AB_y, AB_z);
    const auto Y6_p6 = harm::Y_ll_6_m_p6(AB_x, AB_y, AB_z);
    const auto Y8_n8 = harm::Y_ll_8_m_n8(AB_x, AB_y, AB_z);
    const auto Y8_n7 = harm::Y_ll_8_m_n7(AB_x, AB_y, AB_z);
    const auto Y8_n6 = harm::Y_ll_8_m_n6(AB_x, AB_y, AB_z);
    const auto Y8_n5 = harm::Y_ll_8_m_n5(AB_x, AB_y, AB_z);
    const auto Y8_n4 = harm::Y_ll_8_m_n4(AB_x, AB_y, AB_z);
    const auto Y8_n3 = harm::Y_ll_8_m_n3(AB_x, AB_y, AB_z);
    const auto Y8_n2 = harm::Y_ll_8_m_n2(AB_x, AB_y, AB_z);
    const auto Y8_n1 = harm::Y_ll_8_m_n1(AB_x, AB_y, AB_z);
    const auto Y8_p0 = harm::Y_ll_8_m_p0(AB_x, AB_y, AB_z);
    const auto Y8_p1 = harm::Y_ll_8_m_p1(AB_x, AB_y, AB_z);
    const auto Y8_p2 = harm::Y_ll_8_m_p2(AB_x, AB_y, AB_z);
    const auto Y8_p3 = harm::Y_ll_8_m_p3(AB_x, AB_y, AB_z);
    const auto Y8_p4 = harm::Y_ll_8_m_p4(AB_x, AB_y, AB_z);
    const auto Y8_p5 = harm::Y_ll_8_m_p5(AB_x, AB_y, AB_z);
    const auto Y8_p6 = harm::Y_ll_8_m_p6(AB_x, AB_y, AB_z);
    const auto Y8_p7 = harm::Y_ll_8_m_p7(AB_x, AB_y, AB_z);
    const auto Y8_p8 = harm::Y_ll_8_m_p8(AB_x, AB_y, AB_z);
    const auto Y10_n10 = harm::Y_ll_10_m_n10(AB_x, AB_y, AB_z);
    const auto Y10_n9 = harm::Y_ll_10_m_n9(AB_x, AB_y, AB_z);
    const auto Y10_n8 = harm::Y_ll_10_m_n8(AB_x, AB_y, AB_z);
    const auto Y10_n7 = harm::Y_ll_10_m_n7(AB_x, AB_y, AB_z);
    const auto Y10_n6 = harm::Y_ll_10_m_n6(AB_x, AB_y, AB_z);
    const auto Y10_n5 = harm::Y_ll_10_m_n5(AB_x, AB_y, AB_z);
    const auto Y10_n4 = harm::Y_ll_10_m_n4(AB_x, AB_y, AB_z);
    const auto Y10_n3 = harm::Y_ll_10_m_n3(AB_x, AB_y, AB_z);
    const auto Y10_n2 = harm::Y_ll_10_m_n2(AB_x, AB_y, AB_z);
    const auto Y10_n1 = harm::Y_ll_10_m_n1(AB_x, AB_y, AB_z);
    const auto Y10_p0 = harm::Y_ll_10_m_p0(AB_x, AB_y, AB_z);
    const auto Y10_p1 = harm::Y_ll_10_m_p1(AB_x, AB_y, AB_z);
    const auto Y10_p2 = harm::Y_ll_10_m_p2(AB_x, AB_y, AB_z);
    const auto Y10_p3 = harm::Y_ll_10_m_p3(AB_x, AB_y, AB_z);
    const auto Y10_p4 = harm::Y_ll_10_m_p4(AB_x, AB_y, AB_z);
    const auto Y10_p5 = harm::Y_ll_10_m_p5(AB_x, AB_y, AB_z);
    const auto Y10_p6 = harm::Y_ll_10_m_p6(AB_x, AB_y, AB_z);
    const auto Y10_p7 = harm::Y_ll_10_m_p7(AB_x, AB_y, AB_z);
    const auto Y10_p8 = harm::Y_ll_10_m_p8(AB_x, AB_y, AB_z);
    const auto Y10_p9 = harm::Y_ll_10_m_p9(AB_x, AB_y, AB_z);
    const auto Y10_p10 = harm::Y_ll_10_m_p10(AB_x, AB_y, AB_z);
    const auto R4 = R2 * R2;
    const auto R6 = R4 * R2;
    const auto R8 = R6 * R2;
    const auto R10 = R8 * R2;
    auto *d = buffer;
    d[84] = Y6_p0 * Y6_p0 * V[0] + ((-7938.0 / 4199.0) * Y10_p0 + (-7350.0 / 2717.0) * Y8_p0 * R2 + (-600.0 / 187.0) * Y6_p0 * R4 + (-504.0 / 143.0) * Y4_p0 * R6 + (-525.0 / 143.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((3675.0 / 286.0) * Y8_p0 + (300.0 / 11.0) * Y6_p0 * R2 + (5670.0 / 143.0) * Y4_p0 * R4 + (525.0 / 11.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + ((-750.0 / 11.0) * Y6_p0 + (-1890.0 / 11.0) * Y4_p0 * R2 - 262.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (236.25 * Y4_p0 + 590.625 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (-413.4375 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[85] = d[97] = Y6_p0 * Y6_p1 * V[0] + ((-189.0 / 4199.0) * sqrt1155 * Y10_p1 + (-1050.0 / 2717.0) * sqrt21 * Y8_p1 * R2 + (-300.0 / 187.0) * Y6_p1 * R4 + (-12.0 / 143.0) * sqrt210 * Y4_p1 * R6 + (-75.0 / 286.0) * sqrt7 * Y2_p1 * R8) * V[1] + ((525.0 / 286.0) * sqrt21 * Y8_p1 + (150.0 / 11.0) * Y6_p1 * R2 + (135.0 / 143.0) * sqrt210 * Y4_p1 * R4 + (75.0 / 22.0) * sqrt7 * Y2_p1 * R6) * V[2] + ((-375.0 / 11.0) * Y6_p1 + (-45.0 / 11.0) * sqrt210 * Y4_p1 * R2 - 18.75 * sqrt7 * Y2_p1 * R4) * V[3] + (5.625 * sqrt210 * Y4_p1 + 42.1875 * sqrt7 * Y2_p1 * R2) * V[4] - 29.53125 * sqrt7 * Y2_p1 * V[5];
    d[86] = d[110] = Y6_p0 * Y6_p2 * V[0] + ((-189.0 / 4199.0) * sqrt154 * Y10_p2 + (735.0 / 2717.0) * sqrt3 * Y8_p2 * R2 + (30.0 / 17.0) * Y6_p2 * R4 + (6.0 / 13.0) * sqrt42 * Y4_p2 * R6 + (75.0 / 143.0) * sqrt70 * Y2_p2 * R8) * V[1] + ((-735.0 / 572.0) * sqrt3 * Y8_p2 - 15.0 * Y6_p2 * R2 + (-135.0 / 26.0) * sqrt42 * Y4_p2 * R4 + (-75.0 / 11.0) * sqrt70 * Y2_p2 * R6) * V[2] + (37.5 * Y6_p2 + 22.5 * sqrt42 * Y4_p2 * R2 + 37.5 * sqrt70 * Y2_p2 * R4) * V[3] + (-30.9375 * sqrt42 * Y4_p2 - 84.375 * sqrt70 * Y2_p2 * R2) * V[4] + 59.0625 * sqrt70 * Y2_p2 * V[5];
    d[87] = d[123] = Y6_p0 * Y6_p3 * V[0] + ((189.0 / 8398.0) * sqrt1001 * Y10_p3 + (1470.0 / 2717.0) * sqrt22 * Y8_p3 * R2 + (645.0 / 187.0) * Y6_p3 * R4 + (252.0 / 143.0) * sqrt3 * Y4_p3 * R6) * V[1] + ((-735.0 / 286.0) * sqrt22 * Y8_p3 + (-645.0 / 22.0) * Y6_p3 * R2 + (-2835.0 / 143.0) * sqrt3 * Y4_p3 * R4) * V[2] + ((3225.0 / 44.0) * Y6_p3 + (945.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[3] - 118.125 * sqrt3 * Y4_p3 * V[4];
    d[88] = d[136] = Y6_p0 * Y6_p4 * V[0] + ((63.0 / 4199.0) * sqrt15015 * Y10_p4 + (2205.0 / 2717.0) * sqrt11 * Y8_p4 * R2 + (120.0 / 187.0) * Y6_p4 * R4 + (-252.0 / 143.0) * sqrt5 * Y4_p4 * R6) * V[1] + ((-2205.0 / 572.0) * sqrt11 * Y8_p4 + (-60.0 / 11.0) * Y6_p4 * R2 + (2835.0 / 143.0) * sqrt5 * Y4_p4 * R4) * V[2] + ((150.0 / 11.0) * Y6_p4 + (-945.0 / 11.0) * sqrt5 * Y4_p4 * R2) * V[3] + 118.125 * sqrt5 * Y4_p4 * V[4];
    d[89] = d[149] = Y6_p0 * Y6_p5 * V[0] + ((63.0 / 442.0) * sqrt273 * Y10_p5 + (-75.0 / 17.0) * Y6_p5 * R4) * V[1] + 37.5 * Y6_p5 * R2 * V[2] - 93.75 * Y6_p5 * V[3];
    d[90] = d[162] = Y6_p0 * Y6_p6 * V[0] + ((378.0 / 4199.0) * sqrt455 * Y10_p6 + (-105.0 / 247.0) * sqrt91 * Y8_p6 * R2 + (30.0 / 17.0) * Y6_p6 * R4) * V[1] + ((105.0 / 52.0) * sqrt91 * Y8_p6 - 15.0 * Y6_p6 * R2) * V[2] + 37.5 * Y6_p6 * V[3];
    d[83] = d[71] = Y6_p0 * Y6_n1 * V[0] + ((-189.0 / 4199.0) * sqrt1155 * Y10_n1 + (-1050.0 / 2717.0) * sqrt21 * Y8_n1 * R2 + (-300.0 / 187.0) * Y6_n1 * R4 + (-12.0 / 143.0) * sqrt210 * Y4_n1 * R6 + (-75.0 / 286.0) * sqrt7 * Y2_n1 * R8) * V[1] + ((525.0 / 286.0) * sqrt21 * Y8_n1 + (150.0 / 11.0) * Y6_n1 * R2 + (135.0 / 143.0) * sqrt210 * Y4_n1 * R4 + (75.0 / 22.0) * sqrt7 * Y2_n1 * R6) * V[2] + ((-375.0 / 11.0) * Y6_n1 + (-45.0 / 11.0) * sqrt210 * Y4_n1 * R2 - 18.75 * sqrt7 * Y2_n1 * R4) * V[3] + (5.625 * sqrt210 * Y4_n1 + 42.1875 * sqrt7 * Y2_n1 * R2) * V[4] - 29.53125 * sqrt7 * Y2_n1 * V[5];
    d[82] = d[58] = Y6_p0 * Y6_n2 * V[0] + ((-189.0 / 4199.0) * sqrt154 * Y10_n2 + (735.0 / 2717.0) * sqrt3 * Y8_n2 * R2 + (30.0 / 17.0) * Y6_n2 * R4 + (6.0 / 13.0) * sqrt42 * Y4_n2 * R6 + (75.0 / 143.0) * sqrt70 * Y2_n2 * R8) * V[1] + ((-735.0 / 572.0) * sqrt3 * Y8_n2 - 15.0 * Y6_n2 * R2 + (-135.0 / 26.0) * sqrt42 * Y4_n2 * R4 + (-75.0 / 11.0) * sqrt70 * Y2_n2 * R6) * V[2] + (37.5 * Y6_n2 + 22.5 * sqrt42 * Y4_n2 * R2 + 37.5 * sqrt70 * Y2_n2 * R4) * V[3] + (-30.9375 * sqrt42 * Y4_n2 - 84.375 * sqrt70 * Y2_n2 * R2) * V[4] + 59.0625 * sqrt70 * Y2_n2 * V[5];
    d[81] = d[45] = Y6_p0 * Y6_n3 * V[0] + ((189.0 / 8398.0) * sqrt1001 * Y10_n3 + (1470.0 / 2717.0) * sqrt22 * Y8_n3 * R2 + (645.0 / 187.0) * Y6_n3 * R4 + (252.0 / 143.0) * sqrt3 * Y4_n3 * R6) * V[1] + ((-735.0 / 286.0) * sqrt22 * Y8_n3 + (-645.0 / 22.0) * Y6_n3 * R2 + (-2835.0 / 143.0) * sqrt3 * Y4_n3 * R4) * V[2] + ((3225.0 / 44.0) * Y6_n3 + (945.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[3] - 118.125 * sqrt3 * Y4_n3 * V[4];
    d[80] = d[32] = Y6_p0 * Y6_n4 * V[0] + ((63.0 / 4199.0) * sqrt15015 * Y10_n4 + (2205.0 / 2717.0) * sqrt11 * Y8_n4 * R2 + (120.0 / 187.0) * Y6_n4 * R4 + (-252.0 / 143.0) * sqrt5 * Y4_n4 * R6) * V[1] + ((-2205.0 / 572.0) * sqrt11 * Y8_n4 + (-60.0 / 11.0) * Y6_n4 * R2 + (2835.0 / 143.0) * sqrt5 * Y4_n4 * R4) * V[2] + ((150.0 / 11.0) * Y6_n4 + (-945.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[3] + 118.125 * sqrt5 * Y4_n4 * V[4];
    d[79] = d[19] = Y6_p0 * Y6_n5 * V[0] + ((63.0 / 442.0) * sqrt273 * Y10_n5 + (-75.0 / 17.0) * Y6_n5 * R4) * V[1] + 37.5 * Y6_n5 * R2 * V[2] - 93.75 * Y6_n5 * V[3];
    d[78] = d[6] = Y6_p0 * Y6_n6 * V[0] + ((378.0 / 4199.0) * sqrt455 * Y10_n6 + (-105.0 / 247.0) * sqrt91 * Y8_n6 * R2 + (30.0 / 17.0) * Y6_n6 * R4) * V[1] + ((105.0 / 52.0) * sqrt91 * Y8_n6 - 15.0 * Y6_n6 * R2) * V[2] + 37.5 * Y6_n6 * V[3];
    d[98] = Y6_p1 * Y6_p1 * V[0] + ((189.0 / 323.0) * Y10_p0 + (-441.0 / 4199.0) * sqrt165 * Y10_p2 + (-1050.0 / 2717.0) * Y8_p0 * R2 + (-630.0 / 2717.0) * sqrt70 * Y8_p2 * R2 + (-300.0 / 187.0) * Y6_p0 * R4 + (-30.0 / 187.0) * sqrt210 * Y6_p2 * R4 + (-384.0 / 143.0) * Y4_p0 * R6 + (-168.0 / 143.0) * sqrt5 * Y4_p2 * R6 + (-75.0 / 22.0) * Y2_p0 * R8 + (-525.0 / 286.0) * sqrt3 * Y2_p2 * R8 - 3.0 * R10) * V[1] + ((525.0 / 286.0) * Y8_p0 + (315.0 / 286.0) * sqrt70 * Y8_p2 + (150.0 / 11.0) * Y6_p0 * R2 + (15.0 / 11.0) * sqrt210 * Y6_p2 * R2 + (4320.0 / 143.0) * Y4_p0 * R4 + (1890.0 / 143.0) * sqrt5 * Y4_p2 * R4 + (975.0 / 22.0) * Y2_p0 * R6 + (525.0 / 22.0) * sqrt3 * Y2_p2 * R6 + 41.25 * R8) * V[2] + ((-375.0 / 11.0) * Y6_p0 + (-75.0 / 22.0) * sqrt210 * Y6_p2 + (-1440.0 / 11.0) * Y4_p0 * R2 + (-630.0 / 11.0) * sqrt5 * Y4_p2 * R2 - 243.75 * Y2_p0 * R4 - 131.25 * sqrt3 * Y2_p2 * R4 - 247.5 * R6) * V[3] + (180.0 * Y4_p0 + 78.75 * sqrt5 * Y4_p2 + 548.4375 * Y2_p0 * R2 + 295.3125 * sqrt3 * Y2_p2 * R2 + 649.6875 * R4) * V[4] + (-383.90625 * Y2_p0 - 206.71875 * sqrt3 * Y2_p2 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[99] = d[111] = Y6_p1 * Y6_p2 * V[0] + ((-378.0 / 4199.0) * sqrt22 * Y10_p1 + (-441.0 / 8398.0) * sqrt429 * Y10_p3 + (-2835.0 / 5434.0) * sqrt10 * Y8_p1 * R2 + (-315.0 / 5434.0) * sqrt462 * Y8_p3 * R2 + (-30.0 / 187.0) * sqrt210 * Y6_p1 * R4 + (-45.0 / 187.0) * sqrt21 * Y6_p3 * R4 + (-318.0 / 143.0) * Y4_p1 * R6 + (-42.0 / 143.0) * sqrt7 * Y4_p3 * R6 + (-75.0 / 286.0) * sqrt30 * Y2_p1 * R8) * V[1] + ((2835.0 / 1144.0) * sqrt10 * Y8_p1 + (315.0 / 1144.0) * sqrt462 * Y8_p3 + (15.0 / 11.0) * sqrt210 * Y6_p1 * R2 + (45.0 / 22.0) * sqrt21 * Y6_p3 * R2 + (7155.0 / 286.0) * Y4_p1 * R4 + (945.0 / 286.0) * sqrt7 * Y4_p3 * R4 + (75.0 / 22.0) * sqrt30 * Y2_p1 * R6) * V[2] + ((-75.0 / 22.0) * sqrt210 * Y6_p1 + (-225.0 / 44.0) * sqrt21 * Y6_p3 + (-2385.0 / 22.0) * Y4_p1 * R2 + (-315.0 / 22.0) * sqrt7 * Y4_p3 * R2 - 18.75 * sqrt30 * Y2_p1 * R4) * V[3] + (149.0625 * Y4_p1 + 19.6875 * sqrt7 * Y4_p3 + 42.1875 * sqrt30 * Y2_p1 * R2) * V[4] - 29.53125 * sqrt30 * Y2_p1 * V[5];
    d[100] = d[124] = Y6_p1 * Y6_p3 * V[0] + ((-63.0 / 442.0) * sqrt66 * Y10_p2 + (-189.0 / 16796.0) * sqrt858 * Y10_p4 + (-105.0 / 143.0) * sqrt7 * Y8_p2 * R2 + (105.0 / 5434.0) * sqrt770 * Y8_p4 * R2 + (-45.0 / 187.0) * sqrt21 * Y6_p2 * R4 + (75.0 / 374.0) * sqrt70 * Y6_p4 * R4 + (72.0 / 143.0) * sqrt2 * Y4_p2 * R6 + (126.0 / 143.0) * sqrt14 * Y4_p4 * R6 + (75.0 / 143.0) * sqrt30 * Y2_p2 * R8) * V[1] + ((1995.0 / 572.0) * sqrt7 * Y8_p2 + (-105.0 / 1144.0) * sqrt770 * Y8_p4 + (45.0 / 22.0) * sqrt21 * Y6_p2 * R2 + (-75.0 / 44.0) * sqrt70 * Y6_p4 * R2 + (-810.0 / 143.0) * sqrt2 * Y4_p2 * R4 + (-2835.0 / 286.0) * sqrt14 * Y4_p4 * R4 + (-75.0 / 11.0) * sqrt30 * Y2_p2 * R6) * V[2] + ((-225.0 / 44.0) * sqrt21 * Y6_p2 + (375.0 / 88.0) * sqrt70 * Y6_p4 + (270.0 / 11.0) * sqrt2 * Y4_p2 * R2 + (945.0 / 22.0) * sqrt14 * Y4_p4 * R2 + 37.5 * sqrt30 * Y2_p2 * R4) * V[3] + (-33.75 * sqrt2 * Y4_p2 - 59.0625 * sqrt14 * Y4_p4 - 84.375 * sqrt30 * Y2_p2 * R2) * V[4] + 59.0625 * sqrt30 * Y2_p2 * V[5];
    d[101] = d[137] = Y6_p1 * Y6_p4 * V[0] + ((-693.0 / 16796.0) * sqrt1430 * Y10_p3 + (693.0 / 16796.0) * sqrt286 * Y10_p5 + (-105.0 / 2717.0) * sqrt385 * Y8_p3 * R2 + (105.0 / 2717.0) * sqrt3003 * Y8_p5 * R2 + (75.0 / 374.0) * sqrt70 * Y6_p3 * R4 + (45.0 / 374.0) * sqrt462 * Y6_p5 * R4 + (30.0 / 143.0) * sqrt210 * Y4_p3 * R6) * V[1] + ((105.0 / 572.0) * sqrt385 * Y8_p3 + (-105.0 / 572.0) * sqrt3003 * Y8_p5 + (-75.0 / 44.0) * sqrt70 * Y6_p3 * R2 + (-45.0 / 44.0) * sqrt462 * Y6_p5 * R2 + (-675.0 / 286.0) * sqrt210 * Y4_p3 * R4) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_p3 + (225.0 / 88.0) * sqrt462 * Y6_p5 + (225.0 / 22.0) * sqrt210 * Y4_p3 * R2) * V[3] - 14.0625 * sqrt210 * Y4_p3 * V[4];
    d[102] = d[150] = Y6_p1 * Y6_p5 * V[0] + ((-2205.0 / 16796.0) * sqrt130 * Y10_p4 + (63.0 / 323.0) * sqrt65 * Y10_p6 + (105.0 / 494.0) * sqrt42 * Y8_p4 * R2 + (105.0 / 247.0) * sqrt13 * Y8_p6 * R2 + (45.0 / 374.0) * sqrt462 * Y6_p4 * R4 + (-15.0 / 17.0) * sqrt7 * Y6_p6 * R4 + (-6.0 / 143.0) * sqrt2310 * Y4_p4 * R6) * V[1] + ((-105.0 / 104.0) * sqrt42 * Y8_p4 + (-105.0 / 52.0) * sqrt13 * Y8_p6 + (-45.0 / 44.0) * sqrt462 * Y6_p4 * R2 + 7.5 * sqrt7 * Y6_p6 * R2 + (135.0 / 286.0) * sqrt2310 * Y4_p4 * R4) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_p4 - 18.75 * sqrt7 * Y6_p6 + (-45.0 / 22.0) * sqrt2310 * Y4_p4 * R2) * V[3] + 2.8125 * sqrt2310 * Y4_p4 * V[4];
    d[103] = d[163] = Y6_p1 * Y6_p6 * V[0] + ((-1323.0 / 8398.0) * sqrt39 * Y10_p5 + (126.0 / 4199.0) * sqrt3315 * Y10_p7 + (105.0 / 494.0) * sqrt182 * Y8_p5 * R2 + (-105.0 / 494.0) * sqrt130 * Y8_p7 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_p5 * R4) * V[1] + ((-105.0 / 104.0) * sqrt182 * Y8_p5 + (105.0 / 104.0) * sqrt130 * Y8_p7 + 7.5 * sqrt7 * Y6_p5 * R2) * V[2] - 18.75 * sqrt7 * Y6_p5 * V[3];
    d[96] = d[72] = Y6_p1 * Y6_n1 * V[0] + ((-441.0 / 4199.0) * sqrt165 * Y10_n2 + (-630.0 / 2717.0) * sqrt70 * Y8_n2 * R2 + (-30.0 / 187.0) * sqrt210 * Y6_n2 * R4 + (-168.0 / 143.0) * sqrt5 * Y4_n2 * R6 + (-525.0 / 286.0) * sqrt3 * Y2_n2 * R8) * V[1] + ((315.0 / 286.0) * sqrt70 * Y8_n2 + (15.0 / 11.0) * sqrt210 * Y6_n2 * R2 + (1890.0 / 143.0) * sqrt5 * Y4_n2 * R4 + (525.0 / 22.0) * sqrt3 * Y2_n2 * R6) * V[2] + ((-75.0 / 22.0) * sqrt210 * Y6_n2 + (-630.0 / 11.0) * sqrt5 * Y4_n2 * R2 - 131.25 * sqrt3 * Y2_n2 * R4) * V[3] + (78.75 * sqrt5 * Y4_n2 + 295.3125 * sqrt3 * Y2_n2 * R2) * V[4] - 206.71875 * sqrt3 * Y2_n2 * V[5];
    d[95] = d[59] = Y6_p1 * Y6_n2 * V[0] + ((-441.0 / 8398.0) * sqrt429 * Y10_n3 + (-378.0 / 4199.0) * sqrt22 * Y10_n1 + (-315.0 / 5434.0) * sqrt462 * Y8_n3 * R2 + (-2835.0 / 5434.0) * sqrt10 * Y8_n1 * R2 + (-45.0 / 187.0) * sqrt21 * Y6_n3 * R4 + (-30.0 / 187.0) * sqrt210 * Y6_n1 * R4 + (-42.0 / 143.0) * sqrt7 * Y4_n3 * R6 + (-318.0 / 143.0) * Y4_n1 * R6 + (-75.0 / 286.0) * sqrt30 * Y2_n1 * R8) * V[1] + ((315.0 / 1144.0) * sqrt462 * Y8_n3 + (2835.0 / 1144.0) * sqrt10 * Y8_n1 + (45.0 / 22.0) * sqrt21 * Y6_n3 * R2 + (15.0 / 11.0) * sqrt210 * Y6_n1 * R2 + (945.0 / 286.0) * sqrt7 * Y4_n3 * R4 + (7155.0 / 286.0) * Y4_n1 * R4 + (75.0 / 22.0) * sqrt30 * Y2_n1 * R6) * V[2] + ((-225.0 / 44.0) * sqrt21 * Y6_n3 + (-75.0 / 22.0) * sqrt210 * Y6_n1 + (-315.0 / 22.0) * sqrt7 * Y4_n3 * R2 + (-2385.0 / 22.0) * Y4_n1 * R2 - 18.75 * sqrt30 * Y2_n1 * R4) * V[3] + (19.6875 * sqrt7 * Y4_n3 + 149.0625 * Y4_n1 + 42.1875 * sqrt30 * Y2_n1 * R2) * V[4] - 29.53125 * sqrt30 * Y2_n1 * V[5];
    d[94] = d[46] = Y6_p1 * Y6_n3 * V[0] + ((-189.0 / 16796.0) * sqrt858 * Y10_n4 + (-63.0 / 442.0) * sqrt66 * Y10_n2 + (105.0 / 5434.0) * sqrt770 * Y8_n4 * R2 + (-105.0 / 143.0) * sqrt7 * Y8_n2 * R2 + (75.0 / 374.0) * sqrt70 * Y6_n4 * R4 + (-45.0 / 187.0) * sqrt21 * Y6_n2 * R4 + (126.0 / 143.0) * sqrt14 * Y4_n4 * R6 + (72.0 / 143.0) * sqrt2 * Y4_n2 * R6 + (75.0 / 143.0) * sqrt30 * Y2_n2 * R8) * V[1] + ((-105.0 / 1144.0) * sqrt770 * Y8_n4 + (1995.0 / 572.0) * sqrt7 * Y8_n2 + (-75.0 / 44.0) * sqrt70 * Y6_n4 * R2 + (45.0 / 22.0) * sqrt21 * Y6_n2 * R2 + (-2835.0 / 286.0) * sqrt14 * Y4_n4 * R4 + (-810.0 / 143.0) * sqrt2 * Y4_n2 * R4 + (-75.0 / 11.0) * sqrt30 * Y2_n2 * R6) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_n4 + (-225.0 / 44.0) * sqrt21 * Y6_n2 + (945.0 / 22.0) * sqrt14 * Y4_n4 * R2 + (270.0 / 11.0) * sqrt2 * Y4_n2 * R2 + 37.5 * sqrt30 * Y2_n2 * R4) * V[3] + (-59.0625 * sqrt14 * Y4_n4 - 33.75 * sqrt2 * Y4_n2 - 84.375 * sqrt30 * Y2_n2 * R2) * V[4] + 59.0625 * sqrt30 * Y2_n2 * V[5];
    d[93] = d[33] = Y6_p1 * Y6_n4 * V[0] + ((693.0 / 16796.0) * sqrt286 * Y10_n5 + (-693.0 / 16796.0) * sqrt1430 * Y10_n3 + (105.0 / 2717.0) * sqrt3003 * Y8_n5 * R2 + (-105.0 / 2717.0) * sqrt385 * Y8_n3 * R2 + (45.0 / 374.0) * sqrt462 * Y6_n5 * R4 + (75.0 / 374.0) * sqrt70 * Y6_n3 * R4 + (30.0 / 143.0) * sqrt210 * Y4_n3 * R6) * V[1] + ((-105.0 / 572.0) * sqrt3003 * Y8_n5 + (105.0 / 572.0) * sqrt385 * Y8_n3 + (-45.0 / 44.0) * sqrt462 * Y6_n5 * R2 + (-75.0 / 44.0) * sqrt70 * Y6_n3 * R2 + (-675.0 / 286.0) * sqrt210 * Y4_n3 * R4) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_n5 + (375.0 / 88.0) * sqrt70 * Y6_n3 + (225.0 / 22.0) * sqrt210 * Y4_n3 * R2) * V[3] - 14.0625 * sqrt210 * Y4_n3 * V[4];
    d[92] = d[20] = Y6_p1 * Y6_n5 * V[0] + ((63.0 / 323.0) * sqrt65 * Y10_n6 + (-2205.0 / 16796.0) * sqrt130 * Y10_n4 + (105.0 / 247.0) * sqrt13 * Y8_n6 * R2 + (105.0 / 494.0) * sqrt42 * Y8_n4 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_n6 * R4 + (45.0 / 374.0) * sqrt462 * Y6_n4 * R4 + (-6.0 / 143.0) * sqrt2310 * Y4_n4 * R6) * V[1] + ((-105.0 / 52.0) * sqrt13 * Y8_n6 + (-105.0 / 104.0) * sqrt42 * Y8_n4 + 7.5 * sqrt7 * Y6_n6 * R2 + (-45.0 / 44.0) * sqrt462 * Y6_n4 * R2 + (135.0 / 286.0) * sqrt2310 * Y4_n4 * R4) * V[2] + (-18.75 * sqrt7 * Y6_n6 + (225.0 / 88.0) * sqrt462 * Y6_n4 + (-45.0 / 22.0) * sqrt2310 * Y4_n4 * R2) * V[3] + 2.8125 * sqrt2310 * Y4_n4 * V[4];
    d[91] = d[7] = Y6_p1 * Y6_n6 * V[0] + ((126.0 / 4199.0) * sqrt3315 * Y10_n7 + (-1323.0 / 8398.0) * sqrt39 * Y10_n5 + (-105.0 / 494.0) * sqrt130 * Y8_n7 * R2 + (105.0 / 494.0) * sqrt182 * Y8_n5 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_n5 * R4) * V[1] + ((105.0 / 104.0) * sqrt130 * Y8_n7 + (-105.0 / 104.0) * sqrt182 * Y8_n5 + 7.5 * sqrt7 * Y6_n5 * R2) * V[2] - 18.75 * sqrt7 * Y6_n5 * V[3];
    d[112] = Y6_p2 * Y6_p2 * V[0] + ((6615.0 / 4199.0) * Y10_p0 + (-126.0 / 4199.0) * sqrt2145 * Y10_p4 + (7455.0 / 2717.0) * Y8_p0 * R2 + (-630.0 / 2717.0) * sqrt77 * Y8_p4 * R2 + (30.0 / 17.0) * Y6_p0 * R4 + (-180.0 / 187.0) * sqrt7 * Y6_p4 * R4 + (-6.0 / 13.0) * Y4_p0 * R6 + (-84.0 / 143.0) * sqrt35 * Y4_p4 * R6 + (-375.0 / 143.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((-7455.0 / 572.0) * Y8_p0 + (315.0 / 286.0) * sqrt77 * Y8_p4 - 15.0 * Y6_p0 * R2 + (90.0 / 11.0) * sqrt7 * Y6_p4 * R2 + (135.0 / 26.0) * Y4_p0 * R4 + (945.0 / 143.0) * sqrt35 * Y4_p4 * R4 + (375.0 / 11.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + (37.5 * Y6_p0 + (-225.0 / 11.0) * sqrt7 * Y6_p4 - 22.5 * Y4_p0 * R2 + (-315.0 / 11.0) * sqrt35 * Y4_p4 * R2 - 187.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (30.9375 * Y4_p0 + 39.375 * sqrt35 * Y4_p4 + 421.875 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (-295.3125 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[113] = d[125] = Y6_p2 * Y6_p3 * V[0] + ((1701.0 / 8398.0) * sqrt55 * Y10_p1 + (-315.0 / 8398.0) * sqrt858 * Y10_p5 + (210.0 / 209.0) * Y8_p1 * R2 + (-105.0 / 2717.0) * sqrt1001 * Y8_p5 * R2 + (-45.0 / 187.0) * sqrt21 * Y6_p1 * R4 + (-15.0 / 187.0) * sqrt154 * Y6_p5 * R4 + (-9.0 / 11.0) * sqrt10 * Y4_p1 * R6 + (-375.0 / 286.0) * sqrt3 * Y2_p1 * R8) * V[1] + ((-105.0 / 22.0) * Y8_p1 + (105.0 / 572.0) * sqrt1001 * Y8_p5 + (45.0 / 22.0) * sqrt21 * Y6_p1 * R2 + (15.0 / 22.0) * sqrt154 * Y6_p5 * R2 + (405.0 / 44.0) * sqrt10 * Y4_p1 * R4 + (375.0 / 22.0) * sqrt3 * Y2_p1 * R6) * V[2] + ((-225.0 / 44.0) * sqrt21 * Y6_p1 + (-75.0 / 44.0) * sqrt154 * Y6_p5 + (-1755.0 / 44.0) * sqrt10 * Y4_p1 * R2 - 93.75 * sqrt3 * Y2_p1 * R4) * V[3] + (54.84375 * sqrt10 * Y4_p1 + 210.9375 * sqrt3 * Y2_p1 * R2) * V[4] - 147.65625 * sqrt3 * Y2_p1 * V[5];
    d[114] = d[138] = Y6_p2 * Y6_p4 * V[0] + ((2709.0 / 8398.0) * sqrt22 * Y10_p2 + (-63.0 / 4199.0) * sqrt143 * Y10_p6 + (-420.0 / 2717.0) * sqrt21 * Y8_p2 * R2 + (105.0 / 2717.0) * sqrt715 * Y8_p6 * R2 + (-180.0 / 187.0) * sqrt7 * Y6_p2 * R4 + (30.0 / 187.0) * sqrt385 * Y6_p6 * R4 + (-69.0 / 143.0) * sqrt6 * Y4_p2 * R6 + (225.0 / 286.0) * sqrt10 * Y2_p2 * R8) * V[1] + ((105.0 / 143.0) * sqrt21 * Y8_p2 + (-105.0 / 572.0) * sqrt715 * Y8_p6 + (90.0 / 11.0) * sqrt7 * Y6_p2 * R2 + (-15.0 / 11.0) * sqrt385 * Y6_p6 * R2 + (3105.0 / 572.0) * sqrt6 * Y4_p2 * R4 + (-225.0 / 22.0) * sqrt10 * Y2_p2 * R6) * V[2] + ((-225.0 / 11.0) * sqrt7 * Y6_p2 + (75.0 / 22.0) * sqrt385 * Y6_p6 + (-1035.0 / 44.0) * sqrt6 * Y4_p2 * R2 + 56.25 * sqrt10 * Y2_p2 * R4) * V[3] + (32.34375 * sqrt6 * Y4_p2 - 126.5625 * sqrt10 * Y2_p2 * R2) * V[4] + 88.59375 * sqrt10 * Y2_p2 * V[5];
    d[115] = d[151] = Y6_p2 * Y6_p5 * V[0] + ((1953.0 / 8398.0) * sqrt26 * Y10_p3 + (441.0 / 8398.0) * sqrt442 * Y10_p7 + (-210.0 / 247.0) * sqrt7 * Y8_p3 * R2 + (105.0 / 247.0) * sqrt39 * Y8_p7 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_p3 * R4 + (21.0 / 143.0) * sqrt462 * Y4_p3 * R6) * V[1] + ((105.0 / 26.0) * sqrt7 * Y8_p3 + (-105.0 / 52.0) * sqrt39 * Y8_p7 + (15.0 / 22.0) * sqrt154 * Y6_p3 * R2 + (-945.0 / 572.0) * sqrt462 * Y4_p3 * R4) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_p3 + (315.0 / 44.0) * sqrt462 * Y4_p3 * R2) * V[3] - 9.84375 * sqrt462 * Y4_p3 * V[4];
    d[116] = d[164] = Y6_p2 * Y6_p6 * V[0] + ((441.0 / 4199.0) * sqrt39 * Y10_p4 + (567.0 / 4199.0) * sqrt221 * Y10_p8 + (-105.0 / 247.0) * sqrt35 * Y8_p4 * R2 + (-105.0 / 247.0) * sqrt13 * Y8_p8 * R2 + (30.0 / 187.0) * sqrt385 * Y6_p4 * R4 + (-18.0 / 143.0) * sqrt77 * Y4_p4 * R6) * V[1] + ((105.0 / 52.0) * sqrt35 * Y8_p4 + (105.0 / 52.0) * sqrt13 * Y8_p8 + (-15.0 / 11.0) * sqrt385 * Y6_p4 * R2 + (405.0 / 286.0) * sqrt77 * Y4_p4 * R4) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_p4 + (-135.0 / 22.0) * sqrt77 * Y4_p4 * R2) * V[3] + 8.4375 * sqrt77 * Y4_p4 * V[4];
    d[109] = d[73] = Y6_p2 * Y6_n1 * V[0] + ((-441.0 / 8398.0) * sqrt429 * Y10_n3 + (378.0 / 4199.0) * sqrt22 * Y10_n1 + (-315.0 / 5434.0) * sqrt462 * Y8_n3 * R2 + (2835.0 / 5434.0) * sqrt10 * Y8_n1 * R2 + (-45.0 / 187.0) * sqrt21 * Y6_n3 * R4 + (30.0 / 187.0) * sqrt210 * Y6_n1 * R4 + (-42.0 / 143.0) * sqrt7 * Y4_n3 * R6 + (318.0 / 143.0) * Y4_n1 * R6 + (75.0 / 286.0) * sqrt30 * Y2_n1 * R8) * V[1] + ((315.0 / 1144.0) * sqrt462 * Y8_n3 + (-2835.0 / 1144.0) * sqrt10 * Y8_n1 + (45.0 / 22.0) * sqrt21 * Y6_n3 * R2 + (-15.0 / 11.0) * sqrt210 * Y6_n1 * R2 + (945.0 / 286.0) * sqrt7 * Y4_n3 * R4 + (-7155.0 / 286.0) * Y4_n1 * R4 + (-75.0 / 22.0) * sqrt30 * Y2_n1 * R6) * V[2] + ((-225.0 / 44.0) * sqrt21 * Y6_n3 + (75.0 / 22.0) * sqrt210 * Y6_n1 + (-315.0 / 22.0) * sqrt7 * Y4_n3 * R2 + (2385.0 / 22.0) * Y4_n1 * R2 + 18.75 * sqrt30 * Y2_n1 * R4) * V[3] + (19.6875 * sqrt7 * Y4_n3 - 149.0625 * Y4_n1 - 42.1875 * sqrt30 * Y2_n1 * R2) * V[4] + 29.53125 * sqrt30 * Y2_n1 * V[5];
    d[108] = d[60] = Y6_p2 * Y6_n2 * V[0] + ((-126.0 / 4199.0) * sqrt2145 * Y10_n4 + (-630.0 / 2717.0) * sqrt77 * Y8_n4 * R2 + (-180.0 / 187.0) * sqrt7 * Y6_n4 * R4 + (-84.0 / 143.0) * sqrt35 * Y4_n4 * R6) * V[1] + ((315.0 / 286.0) * sqrt77 * Y8_n4 + (90.0 / 11.0) * sqrt7 * Y6_n4 * R2 + (945.0 / 143.0) * sqrt35 * Y4_n4 * R4) * V[2] + ((-225.0 / 11.0) * sqrt7 * Y6_n4 + (-315.0 / 11.0) * sqrt35 * Y4_n4 * R2) * V[3] + 39.375 * sqrt35 * Y4_n4 * V[4];
    d[107] = d[47] = Y6_p2 * Y6_n3 * V[0] + ((-315.0 / 8398.0) * sqrt858 * Y10_n5 + (1701.0 / 8398.0) * sqrt55 * Y10_n1 + (-105.0 / 2717.0) * sqrt1001 * Y8_n5 * R2 + (210.0 / 209.0) * Y8_n1 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_n5 * R4 + (-45.0 / 187.0) * sqrt21 * Y6_n1 * R4 + (-9.0 / 11.0) * sqrt10 * Y4_n1 * R6 + (-375.0 / 286.0) * sqrt3 * Y2_n1 * R8) * V[1] + ((105.0 / 572.0) * sqrt1001 * Y8_n5 + (-105.0 / 22.0) * Y8_n1 + (15.0 / 22.0) * sqrt154 * Y6_n5 * R2 + (45.0 / 22.0) * sqrt21 * Y6_n1 * R2 + (405.0 / 44.0) * sqrt10 * Y4_n1 * R4 + (375.0 / 22.0) * sqrt3 * Y2_n1 * R6) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_n5 + (-225.0 / 44.0) * sqrt21 * Y6_n1 + (-1755.0 / 44.0) * sqrt10 * Y4_n1 * R2 - 93.75 * sqrt3 * Y2_n1 * R4) * V[3] + (54.84375 * sqrt10 * Y4_n1 + 210.9375 * sqrt3 * Y2_n1 * R2) * V[4] - 147.65625 * sqrt3 * Y2_n1 * V[5];
    d[106] = d[34] = Y6_p2 * Y6_n4 * V[0] + ((-63.0 / 4199.0) * sqrt143 * Y10_n6 + (2709.0 / 8398.0) * sqrt22 * Y10_n2 + (105.0 / 2717.0) * sqrt715 * Y8_n6 * R2 + (-420.0 / 2717.0) * sqrt21 * Y8_n2 * R2 + (30.0 / 187.0) * sqrt385 * Y6_n6 * R4 + (-180.0 / 187.0) * sqrt7 * Y6_n2 * R4 + (-69.0 / 143.0) * sqrt6 * Y4_n2 * R6 + (225.0 / 286.0) * sqrt10 * Y2_n2 * R8) * V[1] + ((-105.0 / 572.0) * sqrt715 * Y8_n6 + (105.0 / 143.0) * sqrt21 * Y8_n2 + (-15.0 / 11.0) * sqrt385 * Y6_n6 * R2 + (90.0 / 11.0) * sqrt7 * Y6_n2 * R2 + (3105.0 / 572.0) * sqrt6 * Y4_n2 * R4 + (-225.0 / 22.0) * sqrt10 * Y2_n2 * R6) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_n6 + (-225.0 / 11.0) * sqrt7 * Y6_n2 + (-1035.0 / 44.0) * sqrt6 * Y4_n2 * R2 + 56.25 * sqrt10 * Y2_n2 * R4) * V[3] + (32.34375 * sqrt6 * Y4_n2 - 126.5625 * sqrt10 * Y2_n2 * R2) * V[4] + 88.59375 * sqrt10 * Y2_n2 * V[5];
    d[105] = d[21] = Y6_p2 * Y6_n5 * V[0] + ((441.0 / 8398.0) * sqrt442 * Y10_n7 + (1953.0 / 8398.0) * sqrt26 * Y10_n3 + (105.0 / 247.0) * sqrt39 * Y8_n7 * R2 + (-210.0 / 247.0) * sqrt7 * Y8_n3 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_n3 * R4 + (21.0 / 143.0) * sqrt462 * Y4_n3 * R6) * V[1] + ((-105.0 / 52.0) * sqrt39 * Y8_n7 + (105.0 / 26.0) * sqrt7 * Y8_n3 + (15.0 / 22.0) * sqrt154 * Y6_n3 * R2 + (-945.0 / 572.0) * sqrt462 * Y4_n3 * R4) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_n3 + (315.0 / 44.0) * sqrt462 * Y4_n3 * R2) * V[3] - 9.84375 * sqrt462 * Y4_n3 * V[4];
    d[104] = d[8] = Y6_p2 * Y6_n6 * V[0] + ((567.0 / 4199.0) * sqrt221 * Y10_n8 + (441.0 / 4199.0) * sqrt39 * Y10_n4 + (-105.0 / 247.0) * sqrt13 * Y8_n8 * R2 + (-105.0 / 247.0) * sqrt35 * Y8_n4 * R2 + (30.0 / 187.0) * sqrt385 * Y6_n4 * R4 + (-18.0 / 143.0) * sqrt77 * Y4_n4 * R6) * V[1] + ((105.0 / 52.0) * sqrt13 * Y8_n8 + (105.0 / 52.0) * sqrt35 * Y8_n4 + (-15.0 / 11.0) * sqrt385 * Y6_n4 * R2 + (405.0 / 286.0) * sqrt77 * Y4_n4 * R4) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_n4 + (-135.0 / 22.0) * sqrt77 * Y4_n4 * R2) * V[3] + 8.4375 * sqrt77 * Y4_n4 * V[4];
    d[126] = Y6_p3 * Y6_p3 * V[0] + ((-945.0 / 442.0) * Y10_p0 + (-189.0 / 8398.0) * sqrt4290 * Y10_p6 + (105.0 / 143.0) * Y8_p0 * R2 + (-210.0 / 2717.0) * sqrt858 * Y8_p6 * R2 + (645.0 / 187.0) * Y6_p0 * R4 + (-30.0 / 187.0) * sqrt462 * Y6_p6 * R4 + (324.0 / 143.0) * Y4_p0 * R6 + (-375.0 / 286.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((-1995.0 / 572.0) * Y8_p0 + (105.0 / 286.0) * sqrt858 * Y8_p6 + (-645.0 / 22.0) * Y6_p0 * R2 + (15.0 / 11.0) * sqrt462 * Y6_p6 * R2 + (-3645.0 / 143.0) * Y4_p0 * R4 + (375.0 / 22.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + ((3225.0 / 44.0) * Y6_p0 + (-75.0 / 22.0) * sqrt462 * Y6_p6 + (1215.0 / 11.0) * Y4_p0 * R2 - 93.75 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (-151.875 * Y4_p0 + 210.9375 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (-147.65625 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[127] = d[139] = Y6_p3 * Y6_p4 * V[0] + ((-2583.0 / 16796.0) * sqrt66 * Y10_p1 + (-189.0 / 8398.0) * sqrt2431 * Y10_p7 + (945.0 / 2717.0) * sqrt30 * Y8_p1 * R2 + (-105.0 / 2717.0) * sqrt858 * Y8_p7 * R2 + (75.0 / 374.0) * sqrt70 * Y6_p1 * R4 + (-126.0 / 143.0) * sqrt3 * Y4_p1 * R6 + (-525.0 / 572.0) * sqrt10 * Y2_p1 * R8) * V[1] + ((-945.0 / 572.0) * sqrt30 * Y8_p1 + (105.0 / 572.0) * sqrt858 * Y8_p7 + (-75.0 / 44.0) * sqrt70 * Y6_p1 * R2 + (2835.0 / 286.0) * sqrt3 * Y4_p1 * R4 + (525.0 / 44.0) * sqrt10 * Y2_p1 * R6) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_p1 + (-945.0 / 22.0) * sqrt3 * Y4_p1 * R2 - 65.625 * sqrt10 * Y2_p1 * R4) * V[3] + (59.0625 * sqrt3 * Y4_p1 + 147.65625 * sqrt10 * Y2_p1 * R2) * V[4] - 103.359375 * sqrt10 * Y2_p1 * V[5];
    d[128] = d[152] = Y6_p3 * Y6_p5 * V[0] + ((-6993.0 / 8398.0) * Y10_p2 + (63.0 / 8398.0) * sqrt663 * Y10_p8 + (315.0 / 2717.0) * sqrt462 * Y8_p2 * R2 + (105.0 / 247.0) * sqrt39 * Y8_p8 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_p2 * R4 + (-72.0 / 143.0) * sqrt33 * Y4_p2 * R6 + (75.0 / 286.0) * sqrt55 * Y2_p2 * R8) * V[1] + ((-315.0 / 572.0) * sqrt462 * Y8_p2 + (-105.0 / 52.0) * sqrt39 * Y8_p8 + (15.0 / 22.0) * sqrt154 * Y6_p2 * R2 + (810.0 / 143.0) * sqrt33 * Y4_p2 * R4 + (-75.0 / 22.0) * sqrt55 * Y2_p2 * R6) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_p2 + (-270.0 / 11.0) * sqrt33 * Y4_p2 * R2 + 18.75 * sqrt55 * Y2_p2 * R4) * V[3] + (33.75 * sqrt33 * Y4_p2 - 42.1875 * sqrt55 * Y2_p2 * R2) * V[4] + 29.53125 * sqrt55 * Y2_p2 * V[5];
    d[129] = d[165] = Y6_p3 * Y6_p6 * V[0] + ((-189.0 / 4199.0) * sqrt78 * Y10_p3 + (189.0 / 8398.0) * sqrt8398 * Y10_p9 + (105.0 / 247.0) * sqrt21 * Y8_p3 * R2 + (-30.0 / 187.0) * sqrt462 * Y6_p3 * R4 + (27.0 / 143.0) * sqrt154 * Y4_p3 * R6) * V[1] + ((-105.0 / 52.0) * sqrt21 * Y8_p3 + (15.0 / 11.0) * sqrt462 * Y6_p3 * R2 + (-1215.0 / 572.0) * sqrt154 * Y4_p3 * R4) * V[2] + ((-75.0 / 22.0) * sqrt462 * Y6_p3 + (405.0 / 44.0) * sqrt154 * Y4_p3 * R2) * V[3] - 12.65625 * sqrt154 * Y4_p3 * V[4];
    d[122] = d[74] = Y6_p3 * Y6_n1 * V[0] + ((-189.0 / 16796.0) * sqrt858 * Y10_n4 + (63.0 / 442.0) * sqrt66 * Y10_n2 + (105.0 / 5434.0) * sqrt770 * Y8_n4 * R2 + (105.0 / 143.0) * sqrt7 * Y8_n2 * R2 + (75.0 / 374.0) * sqrt70 * Y6_n4 * R4 + (45.0 / 187.0) * sqrt21 * Y6_n2 * R4 + (126.0 / 143.0) * sqrt14 * Y4_n4 * R6 + (-72.0 / 143.0) * sqrt2 * Y4_n2 * R6 + (-75.0 / 143.0) * sqrt30 * Y2_n2 * R8) * V[1] + ((-105.0 / 1144.0) * sqrt770 * Y8_n4 + (-1995.0 / 572.0) * sqrt7 * Y8_n2 + (-75.0 / 44.0) * sqrt70 * Y6_n4 * R2 + (-45.0 / 22.0) * sqrt21 * Y6_n2 * R2 + (-2835.0 / 286.0) * sqrt14 * Y4_n4 * R4 + (810.0 / 143.0) * sqrt2 * Y4_n2 * R4 + (75.0 / 11.0) * sqrt30 * Y2_n2 * R6) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_n4 + (225.0 / 44.0) * sqrt21 * Y6_n2 + (945.0 / 22.0) * sqrt14 * Y4_n4 * R2 + (-270.0 / 11.0) * sqrt2 * Y4_n2 * R2 - 37.5 * sqrt30 * Y2_n2 * R4) * V[3] + (-59.0625 * sqrt14 * Y4_n4 + 33.75 * sqrt2 * Y4_n2 + 84.375 * sqrt30 * Y2_n2 * R2) * V[4] - 59.0625 * sqrt30 * Y2_n2 * V[5];
    d[121] = d[61] = Y6_p3 * Y6_n2 * V[0] + ((-315.0 / 8398.0) * sqrt858 * Y10_n5 + (-1701.0 / 8398.0) * sqrt55 * Y10_n1 + (-105.0 / 2717.0) * sqrt1001 * Y8_n5 * R2 + (-210.0 / 209.0) * Y8_n1 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_n5 * R4 + (45.0 / 187.0) * sqrt21 * Y6_n1 * R4 + (9.0 / 11.0) * sqrt10 * Y4_n1 * R6 + (375.0 / 286.0) * sqrt3 * Y2_n1 * R8) * V[1] + ((105.0 / 572.0) * sqrt1001 * Y8_n5 + (105.0 / 22.0) * Y8_n1 + (15.0 / 22.0) * sqrt154 * Y6_n5 * R2 + (-45.0 / 22.0) * sqrt21 * Y6_n1 * R2 + (-405.0 / 44.0) * sqrt10 * Y4_n1 * R4 + (-375.0 / 22.0) * sqrt3 * Y2_n1 * R6) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_n5 + (225.0 / 44.0) * sqrt21 * Y6_n1 + (1755.0 / 44.0) * sqrt10 * Y4_n1 * R2 + 93.75 * sqrt3 * Y2_n1 * R4) * V[3] + (-54.84375 * sqrt10 * Y4_n1 - 210.9375 * sqrt3 * Y2_n1 * R2) * V[4] + 147.65625 * sqrt3 * Y2_n1 * V[5];
    d[120] = d[48] = Y6_p3 * Y6_n3 * V[0] + ((-189.0 / 8398.0) * sqrt4290 * Y10_n6 + (-210.0 / 2717.0) * sqrt858 * Y8_n6 * R2 + (-30.0 / 187.0) * sqrt462 * Y6_n6 * R4) * V[1] + ((105.0 / 286.0) * sqrt858 * Y8_n6 + (15.0 / 11.0) * sqrt462 * Y6_n6 * R2) * V[2] + (-75.0 / 22.0) * sqrt462 * Y6_n6 * V[3];
    d[119] = d[35] = Y6_p3 * Y6_n4 * V[0] + ((-189.0 / 8398.0) * sqrt2431 * Y10_n7 + (-2583.0 / 16796.0) * sqrt66 * Y10_n1 + (-105.0 / 2717.0) * sqrt858 * Y8_n7 * R2 + (945.0 / 2717.0) * sqrt30 * Y8_n1 * R2 + (75.0 / 374.0) * sqrt70 * Y6_n1 * R4 + (-126.0 / 143.0) * sqrt3 * Y4_n1 * R6 + (-525.0 / 572.0) * sqrt10 * Y2_n1 * R8) * V[1] + ((105.0 / 572.0) * sqrt858 * Y8_n7 + (-945.0 / 572.0) * sqrt30 * Y8_n1 + (-75.0 / 44.0) * sqrt70 * Y6_n1 * R2 + (2835.0 / 286.0) * sqrt3 * Y4_n1 * R4 + (525.0 / 44.0) * sqrt10 * Y2_n1 * R6) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_n1 + (-945.0 / 22.0) * sqrt3 * Y4_n1 * R2 - 65.625 * sqrt10 * Y2_n1 * R4) * V[3] + (59.0625 * sqrt3 * Y4_n1 + 147.65625 * sqrt10 * Y2_n1 * R2) * V[4] - 103.359375 * sqrt10 * Y2_n1 * V[5];
    d[118] = d[22] = Y6_p3 * Y6_n5 * V[0] + ((63.0 / 8398.0) * sqrt663 * Y10_n8 + (-6993.0 / 8398.0) * Y10_n2 + (105.0 / 247.0) * sqrt39 * Y8_n8 * R2 + (315.0 / 2717.0) * sqrt462 * Y8_n2 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_n2 * R4 + (-72.0 / 143.0) * sqrt33 * Y4_n2 * R6 + (75.0 / 286.0) * sqrt55 * Y2_n2 * R8) * V[1] + ((-105.0 / 52.0) * sqrt39 * Y8_n8 + (-315.0 / 572.0) * sqrt462 * Y8_n2 + (15.0 / 22.0) * sqrt154 * Y6_n2 * R2 + (810.0 / 143.0) * sqrt33 * Y4_n2 * R4 + (-75.0 / 22.0) * sqrt55 * Y2_n2 * R6) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_n2 + (-270.0 / 11.0) * sqrt33 * Y4_n2 * R2 + 18.75 * sqrt55 * Y2_n2 * R4) * V[3] + (33.75 * sqrt33 * Y4_n2 - 42.1875 * sqrt55 * Y2_n2 * R2) * V[4] + 29.53125 * sqrt55 * Y2_n2 * V[5];
    d[117] = d[9] = Y6_p3 * Y6_n6 * V[0] + ((189.0 / 8398.0) * sqrt8398 * Y10_n9 + (-189.0 / 4199.0) * sqrt78 * Y10_n3 + (105.0 / 247.0) * sqrt21 * Y8_n3 * R2 + (-30.0 / 187.0) * sqrt462 * Y6_n3 * R4 + (27.0 / 143.0) * sqrt154 * Y4_n3 * R6) * V[1] + ((-105.0 / 52.0) * sqrt21 * Y8_n3 + (15.0 / 11.0) * sqrt462 * Y6_n3 * R2 + (-1215.0 / 572.0) * sqrt154 * Y4_n3 * R4) * V[2] + ((-75.0 / 22.0) * sqrt462 * Y6_n3 + (405.0 / 44.0) * sqrt154 * Y4_n3 * R2) * V[3] - 12.65625 * sqrt154 * Y4_n3 * V[4];
    d[140] = Y6_p4 * Y6_p4 * V[0] + ((5229.0 / 4199.0) * Y10_p0 + (-63.0 / 4199.0) * sqrt12155 * Y10_p8 + (-9345.0 / 2717.0) * Y8_p0 * R2 + (-315.0 / 2717.0) * sqrt715 * Y8_p8 * R2 + (120.0 / 187.0) * Y6_p0 * R4 + (576.0 / 143.0) * Y4_p0 * R6 + (75.0 / 143.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((9345.0 / 572.0) * Y8_p0 + (315.0 / 572.0) * sqrt715 * Y8_p8 + (-60.0 / 11.0) * Y6_p0 * R2 + (-6480.0 / 143.0) * Y4_p0 * R4 + (-75.0 / 11.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + ((150.0 / 11.0) * Y6_p0 + (2160.0 / 11.0) * Y4_p0 * R2 + 37.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (-270.0 * Y4_p0 - 84.375 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (59.0625 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[141] = d[153] = Y6_p4 * Y6_p5 * V[0] + ((2709.0 / 16796.0) * sqrt10 * Y10_p1 + (-63.0 / 8398.0) * sqrt20995 * Y10_p9 + (-1260.0 / 2717.0) * sqrt22 * Y8_p1 * R2 + (45.0 / 374.0) * sqrt462 * Y6_p1 * R4 + (18.0 / 143.0) * sqrt55 * Y4_p1 * R6 + (-225.0 / 572.0) * sqrt66 * Y2_p1 * R8) * V[1] + ((315.0 / 143.0) * sqrt22 * Y8_p1 + (-45.0 / 44.0) * sqrt462 * Y6_p1 * R2 + (-405.0 / 286.0) * sqrt55 * Y4_p1 * R4 + (225.0 / 44.0) * sqrt66 * Y2_p1 * R6) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_p1 + (135.0 / 22.0) * sqrt55 * Y4_p1 * R2 - 28.125 * sqrt66 * Y2_p1 * R4) * V[3] + (-8.4375 * sqrt55 * Y4_p1 + 63.28125 * sqrt66 * Y2_p1 * R2) * V[4] - 44.296875 * sqrt66 * Y2_p1 * V[5];
    d[142] = d[166] = Y6_p4 * Y6_p6 * V[0] + ((567.0 / 8398.0) * sqrt10 * Y10_p2 + (63.0 / 4199.0) * sqrt12597 * Y10_p10 + (-105.0 / 2717.0) * sqrt1155 * Y8_p2 * R2 + (30.0 / 187.0) * sqrt385 * Y6_p2 * R4 + (-27.0 / 143.0) * sqrt330 * Y4_p2 * R6 + (75.0 / 286.0) * sqrt22 * Y2_p2 * R8) * V[1] + ((105.0 / 572.0) * sqrt1155 * Y8_p2 + (-15.0 / 11.0) * sqrt385 * Y6_p2 * R2 + (1215.0 / 572.0) * sqrt330 * Y4_p2 * R4 + (-75.0 / 22.0) * sqrt22 * Y2_p2 * R6) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_p2 + (-405.0 / 44.0) * sqrt330 * Y4_p2 * R2 + 18.75 * sqrt22 * Y2_p2 * R4) * V[3] + (12.65625 * sqrt330 * Y4_p2 - 42.1875 * sqrt22 * Y2_p2 * R2) * V[4] + 29.53125 * sqrt22 * Y2_p2 * V[5];
    d[135] = d[75] = Y6_p4 * Y6_n1 * V[0] + ((693.0 / 16796.0) * sqrt286 * Y10_n5 + (693.0 / 16796.0) * sqrt1430 * Y10_n3 + (105.0 / 2717.0) * sqrt3003 * Y8_n5 * R2 + (105.0 / 2717.0) * sqrt385 * Y8_n3 * R2 + (45.0 / 374.0) * sqrt462 * Y6_n5 * R4 + (-75.0 / 374.0) * sqrt70 * Y6_n3 * R4 + (-30.0 / 143.0) * sqrt210 * Y4_n3 * R6) * V[1] + ((-105.0 / 572.0) * sqrt3003 * Y8_n5 + (-105.0 / 572.0) * sqrt385 * Y8_n3 + (-45.0 / 44.0) * sqrt462 * Y6_n5 * R2 + (75.0 / 44.0) * sqrt70 * Y6_n3 * R2 + (675.0 / 286.0) * sqrt210 * Y4_n3 * R4) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_n5 + (-375.0 / 88.0) * sqrt70 * Y6_n3 + (-225.0 / 22.0) * sqrt210 * Y4_n3 * R2) * V[3] + 14.0625 * sqrt210 * Y4_n3 * V[4];
    d[134] = d[62] = Y6_p4 * Y6_n2 * V[0] + ((-63.0 / 4199.0) * sqrt143 * Y10_n6 + (-2709.0 / 8398.0) * sqrt22 * Y10_n2 + (105.0 / 2717.0) * sqrt715 * Y8_n6 * R2 + (420.0 / 2717.0) * sqrt21 * Y8_n2 * R2 + (30.0 / 187.0) * sqrt385 * Y6_n6 * R4 + (180.0 / 187.0) * sqrt7 * Y6_n2 * R4 + (69.0 / 143.0) * sqrt6 * Y4_n2 * R6 + (-225.0 / 286.0) * sqrt10 * Y2_n2 * R8) * V[1] + ((-105.0 / 572.0) * sqrt715 * Y8_n6 + (-105.0 / 143.0) * sqrt21 * Y8_n2 + (-15.0 / 11.0) * sqrt385 * Y6_n6 * R2 + (-90.0 / 11.0) * sqrt7 * Y6_n2 * R2 + (-3105.0 / 572.0) * sqrt6 * Y4_n2 * R4 + (225.0 / 22.0) * sqrt10 * Y2_n2 * R6) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_n6 + (225.0 / 11.0) * sqrt7 * Y6_n2 + (1035.0 / 44.0) * sqrt6 * Y4_n2 * R2 - 56.25 * sqrt10 * Y2_n2 * R4) * V[3] + (-32.34375 * sqrt6 * Y4_n2 + 126.5625 * sqrt10 * Y2_n2 * R2) * V[4] - 88.59375 * sqrt10 * Y2_n2 * V[5];
    d[133] = d[49] = Y6_p4 * Y6_n3 * V[0] + ((-189.0 / 8398.0) * sqrt2431 * Y10_n7 + (2583.0 / 16796.0) * sqrt66 * Y10_n1 + (-105.0 / 2717.0) * sqrt858 * Y8_n7 * R2 + (-945.0 / 2717.0) * sqrt30 * Y8_n1 * R2 + (-75.0 / 374.0) * sqrt70 * Y6_n1 * R4 + (126.0 / 143.0) * sqrt3 * Y4_n1 * R6 + (525.0 / 572.0) * sqrt10 * Y2_n1 * R8) * V[1] + ((105.0 / 572.0) * sqrt858 * Y8_n7 + (945.0 / 572.0) * sqrt30 * Y8_n1 + (75.0 / 44.0) * sqrt70 * Y6_n1 * R2 + (-2835.0 / 286.0) * sqrt3 * Y4_n1 * R4 + (-525.0 / 44.0) * sqrt10 * Y2_n1 * R6) * V[2] + ((-375.0 / 88.0) * sqrt70 * Y6_n1 + (945.0 / 22.0) * sqrt3 * Y4_n1 * R2 + 65.625 * sqrt10 * Y2_n1 * R4) * V[3] + (-59.0625 * sqrt3 * Y4_n1 - 147.65625 * sqrt10 * Y2_n1 * R2) * V[4] + 103.359375 * sqrt10 * Y2_n1 * V[5];
    d[132] = d[36] = Y6_p4 * Y6_n4 * V[0] + ((-63.0 / 4199.0) * sqrt12155 * Y10_n8 + (-315.0 / 2717.0) * sqrt715 * Y8_n8 * R2) * V[1] + (315.0 / 572.0) * sqrt715 * Y8_n8 * V[2];
    d[131] = d[23] = Y6_p4 * Y6_n5 * V[0] + ((-63.0 / 8398.0) * sqrt20995 * Y10_n9 + (2709.0 / 16796.0) * sqrt10 * Y10_n1 + (-1260.0 / 2717.0) * sqrt22 * Y8_n1 * R2 + (45.0 / 374.0) * sqrt462 * Y6_n1 * R4 + (18.0 / 143.0) * sqrt55 * Y4_n1 * R6 + (-225.0 / 572.0) * sqrt66 * Y2_n1 * R8) * V[1] + ((315.0 / 143.0) * sqrt22 * Y8_n1 + (-45.0 / 44.0) * sqrt462 * Y6_n1 * R2 + (-405.0 / 286.0) * sqrt55 * Y4_n1 * R4 + (225.0 / 44.0) * sqrt66 * Y2_n1 * R6) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_n1 + (135.0 / 22.0) * sqrt55 * Y4_n1 * R2 - 28.125 * sqrt66 * Y2_n1 * R4) * V[3] + (-8.4375 * sqrt55 * Y4_n1 + 63.28125 * sqrt66 * Y2_n1 * R2) * V[4] - 44.296875 * sqrt66 * Y2_n1 * V[5];
    d[130] = d[10] = Y6_p4 * Y6_n6 * V[0] + ((63.0 / 4199.0) * sqrt12597 * Y10_n10 + (567.0 / 8398.0) * sqrt10 * Y10_n2 + (-105.0 / 2717.0) * sqrt1155 * Y8_n2 * R2 + (30.0 / 187.0) * sqrt385 * Y6_n2 * R4 + (-27.0 / 143.0) * sqrt330 * Y4_n2 * R6 + (75.0 / 286.0) * sqrt22 * Y2_n2 * R8) * V[1] + ((105.0 / 572.0) * sqrt1155 * Y8_n2 + (-15.0 / 11.0) * sqrt385 * Y6_n2 * R2 + (1215.0 / 572.0) * sqrt330 * Y4_n2 * R4 + (-75.0 / 22.0) * sqrt22 * Y2_n2 * R6) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_n2 + (-405.0 / 44.0) * sqrt330 * Y4_n2 * R2 + 18.75 * sqrt22 * Y2_n2 * R4) * V[3] + (12.65625 * sqrt330 * Y4_n2 - 42.1875 * sqrt22 * Y2_n2 * R2) * V[4] + 29.53125 * sqrt22 * Y2_n2 * V[5];
    d[154] = Y6_p5 * Y6_p5 * V[0] + ((-3087.0 / 8398.0) * Y10_p0 + (-63.0 / 8398.0) * sqrt92378 * Y10_p10 + (525.0 / 247.0) * Y8_p0 * R2 + (-75.0 / 17.0) * Y6_p0 * R4 + (36.0 / 13.0) * Y4_p0 * R6 + (75.0 / 26.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((-525.0 / 52.0) * Y8_p0 + 37.5 * Y6_p0 * R2 + (-405.0 / 13.0) * Y4_p0 * R4 - 37.5 * Y2_p0 * R6 + 41.25 * R8) * V[2] + (-93.75 * Y6_p0 + 135.0 * Y4_p0 * R2 + 206.25 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (-185.625 * Y4_p0 - 464.0625 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (324.84375 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[155] = d[167] = Y6_p5 * Y6_p6 * V[0] + ((-63.0 / 8398.0) * sqrt165 * Y10_p1 + (105.0 / 247.0) * sqrt3 * Y8_p1 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_p1 * R4 + (9.0 / 13.0) * sqrt30 * Y4_p1 * R6 + (-75.0 / 26.0) * Y2_p1 * R8) * V[1] + ((-105.0 / 52.0) * sqrt3 * Y8_p1 + 7.5 * sqrt7 * Y6_p1 * R2 + (-405.0 / 52.0) * sqrt30 * Y4_p1 * R4 + 37.5 * Y2_p1 * R6) * V[2] + (-18.75 * sqrt7 * Y6_p1 + 33.75 * sqrt30 * Y4_p1 * R2 - 206.25 * Y2_p1 * R4) * V[3] + (-46.40625 * sqrt30 * Y4_p1 + 464.0625 * Y2_p1 * R2) * V[4] - 324.84375 * Y2_p1 * V[5];
    d[148] = d[76] = Y6_p5 * Y6_n1 * V[0] + ((63.0 / 323.0) * sqrt65 * Y10_n6 + (2205.0 / 16796.0) * sqrt130 * Y10_n4 + (105.0 / 247.0) * sqrt13 * Y8_n6 * R2 + (-105.0 / 494.0) * sqrt42 * Y8_n4 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_n6 * R4 + (-45.0 / 374.0) * sqrt462 * Y6_n4 * R4 + (6.0 / 143.0) * sqrt2310 * Y4_n4 * R6) * V[1] + ((-105.0 / 52.0) * sqrt13 * Y8_n6 + (105.0 / 104.0) * sqrt42 * Y8_n4 + 7.5 * sqrt7 * Y6_n6 * R2 + (45.0 / 44.0) * sqrt462 * Y6_n4 * R2 + (-135.0 / 286.0) * sqrt2310 * Y4_n4 * R4) * V[2] + (-18.75 * sqrt7 * Y6_n6 + (-225.0 / 88.0) * sqrt462 * Y6_n4 + (45.0 / 22.0) * sqrt2310 * Y4_n4 * R2) * V[3] - 2.8125 * sqrt2310 * Y4_n4 * V[4];
    d[147] = d[63] = Y6_p5 * Y6_n2 * V[0] + ((441.0 / 8398.0) * sqrt442 * Y10_n7 + (-1953.0 / 8398.0) * sqrt26 * Y10_n3 + (105.0 / 247.0) * sqrt39 * Y8_n7 * R2 + (210.0 / 247.0) * sqrt7 * Y8_n3 * R2 + (15.0 / 187.0) * sqrt154 * Y6_n3 * R4 + (-21.0 / 143.0) * sqrt462 * Y4_n3 * R6) * V[1] + ((-105.0 / 52.0) * sqrt39 * Y8_n7 + (-105.0 / 26.0) * sqrt7 * Y8_n3 + (-15.0 / 22.0) * sqrt154 * Y6_n3 * R2 + (945.0 / 572.0) * sqrt462 * Y4_n3 * R4) * V[2] + ((75.0 / 44.0) * sqrt154 * Y6_n3 + (-315.0 / 44.0) * sqrt462 * Y4_n3 * R2) * V[3] + 9.84375 * sqrt462 * Y4_n3 * V[4];
    d[146] = d[50] = Y6_p5 * Y6_n3 * V[0] + ((63.0 / 8398.0) * sqrt663 * Y10_n8 + (6993.0 / 8398.0) * Y10_n2 + (105.0 / 247.0) * sqrt39 * Y8_n8 * R2 + (-315.0 / 2717.0) * sqrt462 * Y8_n2 * R2 + (15.0 / 187.0) * sqrt154 * Y6_n2 * R4 + (72.0 / 143.0) * sqrt33 * Y4_n2 * R6 + (-75.0 / 286.0) * sqrt55 * Y2_n2 * R8) * V[1] + ((-105.0 / 52.0) * sqrt39 * Y8_n8 + (315.0 / 572.0) * sqrt462 * Y8_n2 + (-15.0 / 22.0) * sqrt154 * Y6_n2 * R2 + (-810.0 / 143.0) * sqrt33 * Y4_n2 * R4 + (75.0 / 22.0) * sqrt55 * Y2_n2 * R6) * V[2] + ((75.0 / 44.0) * sqrt154 * Y6_n2 + (270.0 / 11.0) * sqrt33 * Y4_n2 * R2 - 18.75 * sqrt55 * Y2_n2 * R4) * V[3] + (-33.75 * sqrt33 * Y4_n2 + 42.1875 * sqrt55 * Y2_n2 * R2) * V[4] - 29.53125 * sqrt55 * Y2_n2 * V[5];
    d[145] = d[37] = Y6_p5 * Y6_n4 * V[0] + ((-63.0 / 8398.0) * sqrt20995 * Y10_n9 + (-2709.0 / 16796.0) * sqrt10 * Y10_n1 + (1260.0 / 2717.0) * sqrt22 * Y8_n1 * R2 + (-45.0 / 374.0) * sqrt462 * Y6_n1 * R4 + (-18.0 / 143.0) * sqrt55 * Y4_n1 * R6 + (225.0 / 572.0) * sqrt66 * Y2_n1 * R8) * V[1] + ((-315.0 / 143.0) * sqrt22 * Y8_n1 + (45.0 / 44.0) * sqrt462 * Y6_n1 * R2 + (405.0 / 286.0) * sqrt55 * Y4_n1 * R4 + (-225.0 / 44.0) * sqrt66 * Y2_n1 * R6) * V[2] + ((-225.0 / 88.0) * sqrt462 * Y6_n1 + (-135.0 / 22.0) * sqrt55 * Y4_n1 * R2 + 28.125 * sqrt66 * Y2_n1 * R4) * V[3] + (8.4375 * sqrt55 * Y4_n1 - 63.28125 * sqrt66 * Y2_n1 * R2) * V[4] + 44.296875 * sqrt66 * Y2_n1 * V[5];
    d[144] = d[24] = Y6_p5 * Y6_n5 * V[0] + (-63.0 / 8398.0) * sqrt92378 * Y10_n10 * V[1];
    d[143] = d[11] = Y6_p5 * Y6_n6 * V[0] + ((-63.0 / 8398.0) * sqrt165 * Y10_n1 + (105.0 / 247.0) * sqrt3 * Y8_n1 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_n1 * R4 + (9.0 / 13.0) * sqrt30 * Y4_n1 * R6 + (-75.0 / 26.0) * Y2_n1 * R8) * V[1] + ((-105.0 / 52.0) * sqrt3 * Y8_n1 + 7.5 * sqrt7 * Y6_n1 * R2 + (-405.0 / 52.0) * sqrt30 * Y4_n1 * R4 + 37.5 * Y2_n1 * R6) * V[2] + (-18.75 * sqrt7 * Y6_n1 + 33.75 * sqrt30 * Y4_n1 * R2 - 206.25 * Y2_n1 * R4) * V[3] + (-46.40625 * sqrt30 * Y4_n1 + 464.0625 * Y2_n1 * R2) * V[4] - 324.84375 * Y2_n1 * V[5];
    d[168] = Y6_p6 * Y6_p6 * V[0] + ((189.0 / 4199.0) * Y10_p0 + (-105.0 / 247.0) * Y8_p0 * R2 + (30.0 / 17.0) * Y6_p0 * R4 + (-54.0 / 13.0) * Y4_p0 * R6 + (75.0 / 13.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((105.0 / 52.0) * Y8_p0 - 15.0 * Y6_p0 * R2 + (1215.0 / 26.0) * Y4_p0 * R4 - 75.0 * Y2_p0 * R6 + 41.25 * R8) * V[2] + (37.5 * Y6_p0 - 202.5 * Y4_p0 * R2 + 412.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (278.4375 * Y4_p0 - 928.125 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (649.6875 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[161] = d[77] = Y6_p6 * Y6_n1 * V[0] + ((126.0 / 4199.0) * sqrt3315 * Y10_n7 + (1323.0 / 8398.0) * sqrt39 * Y10_n5 + (-105.0 / 494.0) * sqrt130 * Y8_n7 * R2 + (-105.0 / 494.0) * sqrt182 * Y8_n5 * R2 + (15.0 / 17.0) * sqrt7 * Y6_n5 * R4) * V[1] + ((105.0 / 104.0) * sqrt130 * Y8_n7 + (105.0 / 104.0) * sqrt182 * Y8_n5 - 7.5 * sqrt7 * Y6_n5 * R2) * V[2] + 18.75 * sqrt7 * Y6_n5 * V[3];
    d[160] = d[64] = Y6_p6 * Y6_n2 * V[0] + ((567.0 / 4199.0) * sqrt221 * Y10_n8 + (-441.0 / 4199.0) * sqrt39 * Y10_n4 + (-105.0 / 247.0) * sqrt13 * Y8_n8 * R2 + (105.0 / 247.0) * sqrt35 * Y8_n4 * R2 + (-30.0 / 187.0) * sqrt385 * Y6_n4 * R4 + (18.0 / 143.0) * sqrt77 * Y4_n4 * R6) * V[1] + ((105.0 / 52.0) * sqrt13 * Y8_n8 + (-105.0 / 52.0) * sqrt35 * Y8_n4 + (15.0 / 11.0) * sqrt385 * Y6_n4 * R2 + (-405.0 / 286.0) * sqrt77 * Y4_n4 * R4) * V[2] + ((-75.0 / 22.0) * sqrt385 * Y6_n4 + (135.0 / 22.0) * sqrt77 * Y4_n4 * R2) * V[3] - 8.4375 * sqrt77 * Y4_n4 * V[4];
    d[159] = d[51] = Y6_p6 * Y6_n3 * V[0] + ((189.0 / 8398.0) * sqrt8398 * Y10_n9 + (189.0 / 4199.0) * sqrt78 * Y10_n3 + (-105.0 / 247.0) * sqrt21 * Y8_n3 * R2 + (30.0 / 187.0) * sqrt462 * Y6_n3 * R4 + (-27.0 / 143.0) * sqrt154 * Y4_n3 * R6) * V[1] + ((105.0 / 52.0) * sqrt21 * Y8_n3 + (-15.0 / 11.0) * sqrt462 * Y6_n3 * R2 + (1215.0 / 572.0) * sqrt154 * Y4_n3 * R4) * V[2] + ((75.0 / 22.0) * sqrt462 * Y6_n3 + (-405.0 / 44.0) * sqrt154 * Y4_n3 * R2) * V[3] + 12.65625 * sqrt154 * Y4_n3 * V[4];
    d[158] = d[38] = Y6_p6 * Y6_n4 * V[0] + ((63.0 / 4199.0) * sqrt12597 * Y10_n10 + (-567.0 / 8398.0) * sqrt10 * Y10_n2 + (105.0 / 2717.0) * sqrt1155 * Y8_n2 * R2 + (-30.0 / 187.0) * sqrt385 * Y6_n2 * R4 + (27.0 / 143.0) * sqrt330 * Y4_n2 * R6 + (-75.0 / 286.0) * sqrt22 * Y2_n2 * R8) * V[1] + ((-105.0 / 572.0) * sqrt1155 * Y8_n2 + (15.0 / 11.0) * sqrt385 * Y6_n2 * R2 + (-1215.0 / 572.0) * sqrt330 * Y4_n2 * R4 + (75.0 / 22.0) * sqrt22 * Y2_n2 * R6) * V[2] + ((-75.0 / 22.0) * sqrt385 * Y6_n2 + (405.0 / 44.0) * sqrt330 * Y4_n2 * R2 - 18.75 * sqrt22 * Y2_n2 * R4) * V[3] + (-12.65625 * sqrt330 * Y4_n2 + 42.1875 * sqrt22 * Y2_n2 * R2) * V[4] - 29.53125 * sqrt22 * Y2_n2 * V[5];
    d[157] = d[25] = Y6_p6 * Y6_n5 * V[0] + ((63.0 / 8398.0) * sqrt165 * Y10_n1 + (-105.0 / 247.0) * sqrt3 * Y8_n1 * R2 + (15.0 / 17.0) * sqrt7 * Y6_n1 * R4 + (-9.0 / 13.0) * sqrt30 * Y4_n1 * R6 + (75.0 / 26.0) * Y2_n1 * R8) * V[1] + ((105.0 / 52.0) * sqrt3 * Y8_n1 - 7.5 * sqrt7 * Y6_n1 * R2 + (405.0 / 52.0) * sqrt30 * Y4_n1 * R4 - 37.5 * Y2_n1 * R6) * V[2] + (18.75 * sqrt7 * Y6_n1 - 33.75 * sqrt30 * Y4_n1 * R2 + 206.25 * Y2_n1 * R4) * V[3] + (46.40625 * sqrt30 * Y4_n1 - 464.0625 * Y2_n1 * R2) * V[4] + 324.84375 * Y2_n1 * V[5];
    d[156] = d[12] = Y6_p6 * Y6_n6 * V[0];
    d[70] = Y6_n1 * Y6_n1 * V[0] + ((189.0 / 323.0) * Y10_p0 + (441.0 / 4199.0) * sqrt165 * Y10_p2 + (-1050.0 / 2717.0) * Y8_p0 * R2 + (630.0 / 2717.0) * sqrt70 * Y8_p2 * R2 + (-300.0 / 187.0) * Y6_p0 * R4 + (30.0 / 187.0) * sqrt210 * Y6_p2 * R4 + (-384.0 / 143.0) * Y4_p0 * R6 + (168.0 / 143.0) * sqrt5 * Y4_p2 * R6 + (-75.0 / 22.0) * Y2_p0 * R8 + (525.0 / 286.0) * sqrt3 * Y2_p2 * R8 - 3.0 * R10) * V[1] + ((525.0 / 286.0) * Y8_p0 + (-315.0 / 286.0) * sqrt70 * Y8_p2 + (150.0 / 11.0) * Y6_p0 * R2 + (-15.0 / 11.0) * sqrt210 * Y6_p2 * R2 + (4320.0 / 143.0) * Y4_p0 * R4 + (-1890.0 / 143.0) * sqrt5 * Y4_p2 * R4 + (975.0 / 22.0) * Y2_p0 * R6 + (-525.0 / 22.0) * sqrt3 * Y2_p2 * R6 + 41.25 * R8) * V[2] + ((-375.0 / 11.0) * Y6_p0 + (75.0 / 22.0) * sqrt210 * Y6_p2 + (-1440.0 / 11.0) * Y4_p0 * R2 + (630.0 / 11.0) * sqrt5 * Y4_p2 * R2 - 243.75 * Y2_p0 * R4 + 131.25 * sqrt3 * Y2_p2 * R4 - 247.5 * R6) * V[3] + (180.0 * Y4_p0 - 78.75 * sqrt5 * Y4_p2 + 548.4375 * Y2_p0 * R2 - 295.3125 * sqrt3 * Y2_p2 * R2 + 649.6875 * R4) * V[4] + (-383.90625 * Y2_p0 + 206.71875 * sqrt3 * Y2_p2 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[69] = d[57] = Y6_n1 * Y6_n2 * V[0] + ((-378.0 / 4199.0) * sqrt22 * Y10_p1 + (441.0 / 8398.0) * sqrt429 * Y10_p3 + (-2835.0 / 5434.0) * sqrt10 * Y8_p1 * R2 + (315.0 / 5434.0) * sqrt462 * Y8_p3 * R2 + (-30.0 / 187.0) * sqrt210 * Y6_p1 * R4 + (45.0 / 187.0) * sqrt21 * Y6_p3 * R4 + (-318.0 / 143.0) * Y4_p1 * R6 + (42.0 / 143.0) * sqrt7 * Y4_p3 * R6 + (-75.0 / 286.0) * sqrt30 * Y2_p1 * R8) * V[1] + ((2835.0 / 1144.0) * sqrt10 * Y8_p1 + (-315.0 / 1144.0) * sqrt462 * Y8_p3 + (15.0 / 11.0) * sqrt210 * Y6_p1 * R2 + (-45.0 / 22.0) * sqrt21 * Y6_p3 * R2 + (7155.0 / 286.0) * Y4_p1 * R4 + (-945.0 / 286.0) * sqrt7 * Y4_p3 * R4 + (75.0 / 22.0) * sqrt30 * Y2_p1 * R6) * V[2] + ((-75.0 / 22.0) * sqrt210 * Y6_p1 + (225.0 / 44.0) * sqrt21 * Y6_p3 + (-2385.0 / 22.0) * Y4_p1 * R2 + (315.0 / 22.0) * sqrt7 * Y4_p3 * R2 - 18.75 * sqrt30 * Y2_p1 * R4) * V[3] + (149.0625 * Y4_p1 - 19.6875 * sqrt7 * Y4_p3 + 42.1875 * sqrt30 * Y2_p1 * R2) * V[4] - 29.53125 * sqrt30 * Y2_p1 * V[5];
    d[68] = d[44] = Y6_n1 * Y6_n3 * V[0] + ((-63.0 / 442.0) * sqrt66 * Y10_p2 + (189.0 / 16796.0) * sqrt858 * Y10_p4 + (-105.0 / 143.0) * sqrt7 * Y8_p2 * R2 + (-105.0 / 5434.0) * sqrt770 * Y8_p4 * R2 + (-45.0 / 187.0) * sqrt21 * Y6_p2 * R4 + (-75.0 / 374.0) * sqrt70 * Y6_p4 * R4 + (72.0 / 143.0) * sqrt2 * Y4_p2 * R6 + (-126.0 / 143.0) * sqrt14 * Y4_p4 * R6 + (75.0 / 143.0) * sqrt30 * Y2_p2 * R8) * V[1] + ((1995.0 / 572.0) * sqrt7 * Y8_p2 + (105.0 / 1144.0) * sqrt770 * Y8_p4 + (45.0 / 22.0) * sqrt21 * Y6_p2 * R2 + (75.0 / 44.0) * sqrt70 * Y6_p4 * R2 + (-810.0 / 143.0) * sqrt2 * Y4_p2 * R4 + (2835.0 / 286.0) * sqrt14 * Y4_p4 * R4 + (-75.0 / 11.0) * sqrt30 * Y2_p2 * R6) * V[2] + ((-225.0 / 44.0) * sqrt21 * Y6_p2 + (-375.0 / 88.0) * sqrt70 * Y6_p4 + (270.0 / 11.0) * sqrt2 * Y4_p2 * R2 + (-945.0 / 22.0) * sqrt14 * Y4_p4 * R2 + 37.5 * sqrt30 * Y2_p2 * R4) * V[3] + (-33.75 * sqrt2 * Y4_p2 + 59.0625 * sqrt14 * Y4_p4 - 84.375 * sqrt30 * Y2_p2 * R2) * V[4] + 59.0625 * sqrt30 * Y2_p2 * V[5];
    d[67] = d[31] = Y6_n1 * Y6_n4 * V[0] + ((-693.0 / 16796.0) * sqrt1430 * Y10_p3 + (-693.0 / 16796.0) * sqrt286 * Y10_p5 + (-105.0 / 2717.0) * sqrt385 * Y8_p3 * R2 + (-105.0 / 2717.0) * sqrt3003 * Y8_p5 * R2 + (75.0 / 374.0) * sqrt70 * Y6_p3 * R4 + (-45.0 / 374.0) * sqrt462 * Y6_p5 * R4 + (30.0 / 143.0) * sqrt210 * Y4_p3 * R6) * V[1] + ((105.0 / 572.0) * sqrt385 * Y8_p3 + (105.0 / 572.0) * sqrt3003 * Y8_p5 + (-75.0 / 44.0) * sqrt70 * Y6_p3 * R2 + (45.0 / 44.0) * sqrt462 * Y6_p5 * R2 + (-675.0 / 286.0) * sqrt210 * Y4_p3 * R4) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_p3 + (-225.0 / 88.0) * sqrt462 * Y6_p5 + (225.0 / 22.0) * sqrt210 * Y4_p3 * R2) * V[3] - 14.0625 * sqrt210 * Y4_p3 * V[4];
    d[66] = d[18] = Y6_n1 * Y6_n5 * V[0] + ((-2205.0 / 16796.0) * sqrt130 * Y10_p4 + (-63.0 / 323.0) * sqrt65 * Y10_p6 + (105.0 / 494.0) * sqrt42 * Y8_p4 * R2 + (-105.0 / 247.0) * sqrt13 * Y8_p6 * R2 + (45.0 / 374.0) * sqrt462 * Y6_p4 * R4 + (15.0 / 17.0) * sqrt7 * Y6_p6 * R4 + (-6.0 / 143.0) * sqrt2310 * Y4_p4 * R6) * V[1] + ((-105.0 / 104.0) * sqrt42 * Y8_p4 + (105.0 / 52.0) * sqrt13 * Y8_p6 + (-45.0 / 44.0) * sqrt462 * Y6_p4 * R2 - 7.5 * sqrt7 * Y6_p6 * R2 + (135.0 / 286.0) * sqrt2310 * Y4_p4 * R4) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_p4 + 18.75 * sqrt7 * Y6_p6 + (-45.0 / 22.0) * sqrt2310 * Y4_p4 * R2) * V[3] + 2.8125 * sqrt2310 * Y4_p4 * V[4];
    d[65] = d[5] = Y6_n1 * Y6_n6 * V[0] + ((-1323.0 / 8398.0) * sqrt39 * Y10_p5 + (-126.0 / 4199.0) * sqrt3315 * Y10_p7 + (105.0 / 494.0) * sqrt182 * Y8_p5 * R2 + (105.0 / 494.0) * sqrt130 * Y8_p7 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_p5 * R4) * V[1] + ((-105.0 / 104.0) * sqrt182 * Y8_p5 + (-105.0 / 104.0) * sqrt130 * Y8_p7 + 7.5 * sqrt7 * Y6_p5 * R2) * V[2] - 18.75 * sqrt7 * Y6_p5 * V[3];
    d[56] = Y6_n2 * Y6_n2 * V[0] + ((6615.0 / 4199.0) * Y10_p0 + (126.0 / 4199.0) * sqrt2145 * Y10_p4 + (7455.0 / 2717.0) * Y8_p0 * R2 + (630.0 / 2717.0) * sqrt77 * Y8_p4 * R2 + (30.0 / 17.0) * Y6_p0 * R4 + (180.0 / 187.0) * sqrt7 * Y6_p4 * R4 + (-6.0 / 13.0) * Y4_p0 * R6 + (84.0 / 143.0) * sqrt35 * Y4_p4 * R6 + (-375.0 / 143.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((-7455.0 / 572.0) * Y8_p0 + (-315.0 / 286.0) * sqrt77 * Y8_p4 - 15.0 * Y6_p0 * R2 + (-90.0 / 11.0) * sqrt7 * Y6_p4 * R2 + (135.0 / 26.0) * Y4_p0 * R4 + (-945.0 / 143.0) * sqrt35 * Y4_p4 * R4 + (375.0 / 11.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + (37.5 * Y6_p0 + (225.0 / 11.0) * sqrt7 * Y6_p4 - 22.5 * Y4_p0 * R2 + (315.0 / 11.0) * sqrt35 * Y4_p4 * R2 - 187.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (30.9375 * Y4_p0 - 39.375 * sqrt35 * Y4_p4 + 421.875 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (-295.3125 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[55] = d[43] = Y6_n2 * Y6_n3 * V[0] + ((1701.0 / 8398.0) * sqrt55 * Y10_p1 + (315.0 / 8398.0) * sqrt858 * Y10_p5 + (210.0 / 209.0) * Y8_p1 * R2 + (105.0 / 2717.0) * sqrt1001 * Y8_p5 * R2 + (-45.0 / 187.0) * sqrt21 * Y6_p1 * R4 + (15.0 / 187.0) * sqrt154 * Y6_p5 * R4 + (-9.0 / 11.0) * sqrt10 * Y4_p1 * R6 + (-375.0 / 286.0) * sqrt3 * Y2_p1 * R8) * V[1] + ((-105.0 / 22.0) * Y8_p1 + (-105.0 / 572.0) * sqrt1001 * Y8_p5 + (45.0 / 22.0) * sqrt21 * Y6_p1 * R2 + (-15.0 / 22.0) * sqrt154 * Y6_p5 * R2 + (405.0 / 44.0) * sqrt10 * Y4_p1 * R4 + (375.0 / 22.0) * sqrt3 * Y2_p1 * R6) * V[2] + ((-225.0 / 44.0) * sqrt21 * Y6_p1 + (75.0 / 44.0) * sqrt154 * Y6_p5 + (-1755.0 / 44.0) * sqrt10 * Y4_p1 * R2 - 93.75 * sqrt3 * Y2_p1 * R4) * V[3] + (54.84375 * sqrt10 * Y4_p1 + 210.9375 * sqrt3 * Y2_p1 * R2) * V[4] - 147.65625 * sqrt3 * Y2_p1 * V[5];
    d[54] = d[30] = Y6_n2 * Y6_n4 * V[0] + ((2709.0 / 8398.0) * sqrt22 * Y10_p2 + (63.0 / 4199.0) * sqrt143 * Y10_p6 + (-420.0 / 2717.0) * sqrt21 * Y8_p2 * R2 + (-105.0 / 2717.0) * sqrt715 * Y8_p6 * R2 + (-180.0 / 187.0) * sqrt7 * Y6_p2 * R4 + (-30.0 / 187.0) * sqrt385 * Y6_p6 * R4 + (-69.0 / 143.0) * sqrt6 * Y4_p2 * R6 + (225.0 / 286.0) * sqrt10 * Y2_p2 * R8) * V[1] + ((105.0 / 143.0) * sqrt21 * Y8_p2 + (105.0 / 572.0) * sqrt715 * Y8_p6 + (90.0 / 11.0) * sqrt7 * Y6_p2 * R2 + (15.0 / 11.0) * sqrt385 * Y6_p6 * R2 + (3105.0 / 572.0) * sqrt6 * Y4_p2 * R4 + (-225.0 / 22.0) * sqrt10 * Y2_p2 * R6) * V[2] + ((-225.0 / 11.0) * sqrt7 * Y6_p2 + (-75.0 / 22.0) * sqrt385 * Y6_p6 + (-1035.0 / 44.0) * sqrt6 * Y4_p2 * R2 + 56.25 * sqrt10 * Y2_p2 * R4) * V[3] + (32.34375 * sqrt6 * Y4_p2 - 126.5625 * sqrt10 * Y2_p2 * R2) * V[4] + 88.59375 * sqrt10 * Y2_p2 * V[5];
    d[53] = d[17] = Y6_n2 * Y6_n5 * V[0] + ((1953.0 / 8398.0) * sqrt26 * Y10_p3 + (-441.0 / 8398.0) * sqrt442 * Y10_p7 + (-210.0 / 247.0) * sqrt7 * Y8_p3 * R2 + (-105.0 / 247.0) * sqrt39 * Y8_p7 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_p3 * R4 + (21.0 / 143.0) * sqrt462 * Y4_p3 * R6) * V[1] + ((105.0 / 26.0) * sqrt7 * Y8_p3 + (105.0 / 52.0) * sqrt39 * Y8_p7 + (15.0 / 22.0) * sqrt154 * Y6_p3 * R2 + (-945.0 / 572.0) * sqrt462 * Y4_p3 * R4) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_p3 + (315.0 / 44.0) * sqrt462 * Y4_p3 * R2) * V[3] - 9.84375 * sqrt462 * Y4_p3 * V[4];
    d[52] = d[4] = Y6_n2 * Y6_n6 * V[0] + ((441.0 / 4199.0) * sqrt39 * Y10_p4 + (-567.0 / 4199.0) * sqrt221 * Y10_p8 + (-105.0 / 247.0) * sqrt35 * Y8_p4 * R2 + (105.0 / 247.0) * sqrt13 * Y8_p8 * R2 + (30.0 / 187.0) * sqrt385 * Y6_p4 * R4 + (-18.0 / 143.0) * sqrt77 * Y4_p4 * R6) * V[1] + ((105.0 / 52.0) * sqrt35 * Y8_p4 + (-105.0 / 52.0) * sqrt13 * Y8_p8 + (-15.0 / 11.0) * sqrt385 * Y6_p4 * R2 + (405.0 / 286.0) * sqrt77 * Y4_p4 * R4) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_p4 + (-135.0 / 22.0) * sqrt77 * Y4_p4 * R2) * V[3] + 8.4375 * sqrt77 * Y4_p4 * V[4];
    d[42] = Y6_n3 * Y6_n3 * V[0] + ((-945.0 / 442.0) * Y10_p0 + (189.0 / 8398.0) * sqrt4290 * Y10_p6 + (105.0 / 143.0) * Y8_p0 * R2 + (210.0 / 2717.0) * sqrt858 * Y8_p6 * R2 + (645.0 / 187.0) * Y6_p0 * R4 + (30.0 / 187.0) * sqrt462 * Y6_p6 * R4 + (324.0 / 143.0) * Y4_p0 * R6 + (-375.0 / 286.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((-1995.0 / 572.0) * Y8_p0 + (-105.0 / 286.0) * sqrt858 * Y8_p6 + (-645.0 / 22.0) * Y6_p0 * R2 + (-15.0 / 11.0) * sqrt462 * Y6_p6 * R2 + (-3645.0 / 143.0) * Y4_p0 * R4 + (375.0 / 22.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + ((3225.0 / 44.0) * Y6_p0 + (75.0 / 22.0) * sqrt462 * Y6_p6 + (1215.0 / 11.0) * Y4_p0 * R2 - 93.75 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (-151.875 * Y4_p0 + 210.9375 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (-147.65625 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[41] = d[29] = Y6_n3 * Y6_n4 * V[0] + ((-2583.0 / 16796.0) * sqrt66 * Y10_p1 + (189.0 / 8398.0) * sqrt2431 * Y10_p7 + (945.0 / 2717.0) * sqrt30 * Y8_p1 * R2 + (105.0 / 2717.0) * sqrt858 * Y8_p7 * R2 + (75.0 / 374.0) * sqrt70 * Y6_p1 * R4 + (-126.0 / 143.0) * sqrt3 * Y4_p1 * R6 + (-525.0 / 572.0) * sqrt10 * Y2_p1 * R8) * V[1] + ((-945.0 / 572.0) * sqrt30 * Y8_p1 + (-105.0 / 572.0) * sqrt858 * Y8_p7 + (-75.0 / 44.0) * sqrt70 * Y6_p1 * R2 + (2835.0 / 286.0) * sqrt3 * Y4_p1 * R4 + (525.0 / 44.0) * sqrt10 * Y2_p1 * R6) * V[2] + ((375.0 / 88.0) * sqrt70 * Y6_p1 + (-945.0 / 22.0) * sqrt3 * Y4_p1 * R2 - 65.625 * sqrt10 * Y2_p1 * R4) * V[3] + (59.0625 * sqrt3 * Y4_p1 + 147.65625 * sqrt10 * Y2_p1 * R2) * V[4] - 103.359375 * sqrt10 * Y2_p1 * V[5];
    d[40] = d[16] = Y6_n3 * Y6_n5 * V[0] + ((-6993.0 / 8398.0) * Y10_p2 + (-63.0 / 8398.0) * sqrt663 * Y10_p8 + (315.0 / 2717.0) * sqrt462 * Y8_p2 * R2 + (-105.0 / 247.0) * sqrt39 * Y8_p8 * R2 + (-15.0 / 187.0) * sqrt154 * Y6_p2 * R4 + (-72.0 / 143.0) * sqrt33 * Y4_p2 * R6 + (75.0 / 286.0) * sqrt55 * Y2_p2 * R8) * V[1] + ((-315.0 / 572.0) * sqrt462 * Y8_p2 + (105.0 / 52.0) * sqrt39 * Y8_p8 + (15.0 / 22.0) * sqrt154 * Y6_p2 * R2 + (810.0 / 143.0) * sqrt33 * Y4_p2 * R4 + (-75.0 / 22.0) * sqrt55 * Y2_p2 * R6) * V[2] + ((-75.0 / 44.0) * sqrt154 * Y6_p2 + (-270.0 / 11.0) * sqrt33 * Y4_p2 * R2 + 18.75 * sqrt55 * Y2_p2 * R4) * V[3] + (33.75 * sqrt33 * Y4_p2 - 42.1875 * sqrt55 * Y2_p2 * R2) * V[4] + 29.53125 * sqrt55 * Y2_p2 * V[5];
    d[39] = d[3] = Y6_n3 * Y6_n6 * V[0] + ((-189.0 / 4199.0) * sqrt78 * Y10_p3 + (-189.0 / 8398.0) * sqrt8398 * Y10_p9 + (105.0 / 247.0) * sqrt21 * Y8_p3 * R2 + (-30.0 / 187.0) * sqrt462 * Y6_p3 * R4 + (27.0 / 143.0) * sqrt154 * Y4_p3 * R6) * V[1] + ((-105.0 / 52.0) * sqrt21 * Y8_p3 + (15.0 / 11.0) * sqrt462 * Y6_p3 * R2 + (-1215.0 / 572.0) * sqrt154 * Y4_p3 * R4) * V[2] + ((-75.0 / 22.0) * sqrt462 * Y6_p3 + (405.0 / 44.0) * sqrt154 * Y4_p3 * R2) * V[3] - 12.65625 * sqrt154 * Y4_p3 * V[4];
    d[28] = Y6_n4 * Y6_n4 * V[0] + ((5229.0 / 4199.0) * Y10_p0 + (63.0 / 4199.0) * sqrt12155 * Y10_p8 + (-9345.0 / 2717.0) * Y8_p0 * R2 + (315.0 / 2717.0) * sqrt715 * Y8_p8 * R2 + (120.0 / 187.0) * Y6_p0 * R4 + (576.0 / 143.0) * Y4_p0 * R6 + (75.0 / 143.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((9345.0 / 572.0) * Y8_p0 + (-315.0 / 572.0) * sqrt715 * Y8_p8 + (-60.0 / 11.0) * Y6_p0 * R2 + (-6480.0 / 143.0) * Y4_p0 * R4 + (-75.0 / 11.0) * Y2_p0 * R6 + 41.25 * R8) * V[2] + ((150.0 / 11.0) * Y6_p0 + (2160.0 / 11.0) * Y4_p0 * R2 + 37.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (-270.0 * Y4_p0 - 84.375 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (59.0625 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[27] = d[15] = Y6_n4 * Y6_n5 * V[0] + ((2709.0 / 16796.0) * sqrt10 * Y10_p1 + (63.0 / 8398.0) * sqrt20995 * Y10_p9 + (-1260.0 / 2717.0) * sqrt22 * Y8_p1 * R2 + (45.0 / 374.0) * sqrt462 * Y6_p1 * R4 + (18.0 / 143.0) * sqrt55 * Y4_p1 * R6 + (-225.0 / 572.0) * sqrt66 * Y2_p1 * R8) * V[1] + ((315.0 / 143.0) * sqrt22 * Y8_p1 + (-45.0 / 44.0) * sqrt462 * Y6_p1 * R2 + (-405.0 / 286.0) * sqrt55 * Y4_p1 * R4 + (225.0 / 44.0) * sqrt66 * Y2_p1 * R6) * V[2] + ((225.0 / 88.0) * sqrt462 * Y6_p1 + (135.0 / 22.0) * sqrt55 * Y4_p1 * R2 - 28.125 * sqrt66 * Y2_p1 * R4) * V[3] + (-8.4375 * sqrt55 * Y4_p1 + 63.28125 * sqrt66 * Y2_p1 * R2) * V[4] - 44.296875 * sqrt66 * Y2_p1 * V[5];
    d[26] = d[2] = Y6_n4 * Y6_n6 * V[0] + ((567.0 / 8398.0) * sqrt10 * Y10_p2 + (-63.0 / 4199.0) * sqrt12597 * Y10_p10 + (-105.0 / 2717.0) * sqrt1155 * Y8_p2 * R2 + (30.0 / 187.0) * sqrt385 * Y6_p2 * R4 + (-27.0 / 143.0) * sqrt330 * Y4_p2 * R6 + (75.0 / 286.0) * sqrt22 * Y2_p2 * R8) * V[1] + ((105.0 / 572.0) * sqrt1155 * Y8_p2 + (-15.0 / 11.0) * sqrt385 * Y6_p2 * R2 + (1215.0 / 572.0) * sqrt330 * Y4_p2 * R4 + (-75.0 / 22.0) * sqrt22 * Y2_p2 * R6) * V[2] + ((75.0 / 22.0) * sqrt385 * Y6_p2 + (-405.0 / 44.0) * sqrt330 * Y4_p2 * R2 + 18.75 * sqrt22 * Y2_p2 * R4) * V[3] + (12.65625 * sqrt330 * Y4_p2 - 42.1875 * sqrt22 * Y2_p2 * R2) * V[4] + 29.53125 * sqrt22 * Y2_p2 * V[5];
    d[14] = Y6_n5 * Y6_n5 * V[0] + ((-3087.0 / 8398.0) * Y10_p0 + (63.0 / 8398.0) * sqrt92378 * Y10_p10 + (525.0 / 247.0) * Y8_p0 * R2 + (-75.0 / 17.0) * Y6_p0 * R4 + (36.0 / 13.0) * Y4_p0 * R6 + (75.0 / 26.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((-525.0 / 52.0) * Y8_p0 + 37.5 * Y6_p0 * R2 + (-405.0 / 13.0) * Y4_p0 * R4 - 37.5 * Y2_p0 * R6 + 41.25 * R8) * V[2] + (-93.75 * Y6_p0 + 135.0 * Y4_p0 * R2 + 206.25 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (-185.625 * Y4_p0 - 464.0625 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (324.84375 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
    d[13] = d[1] = Y6_n5 * Y6_n6 * V[0] + ((-63.0 / 8398.0) * sqrt165 * Y10_p1 + (105.0 / 247.0) * sqrt3 * Y8_p1 * R2 + (-15.0 / 17.0) * sqrt7 * Y6_p1 * R4 + (9.0 / 13.0) * sqrt30 * Y4_p1 * R6 + (-75.0 / 26.0) * Y2_p1 * R8) * V[1] + ((-105.0 / 52.0) * sqrt3 * Y8_p1 + 7.5 * sqrt7 * Y6_p1 * R2 + (-405.0 / 52.0) * sqrt30 * Y4_p1 * R4 + 37.5 * Y2_p1 * R6) * V[2] + (-18.75 * sqrt7 * Y6_p1 + 33.75 * sqrt30 * Y4_p1 * R2 - 206.25 * Y2_p1 * R4) * V[3] + (-46.40625 * sqrt30 * Y4_p1 + 464.0625 * Y2_p1 * R2) * V[4] - 324.84375 * Y2_p1 * V[5];
    d[0] = Y6_n6 * Y6_n6 * V[0] + ((189.0 / 4199.0) * Y10_p0 + (-105.0 / 247.0) * Y8_p0 * R2 + (30.0 / 17.0) * Y6_p0 * R4 + (-54.0 / 13.0) * Y4_p0 * R6 + (75.0 / 13.0) * Y2_p0 * R8 - 3.0 * R10) * V[1] + ((105.0 / 52.0) * Y8_p0 - 15.0 * Y6_p0 * R2 + (1215.0 / 26.0) * Y4_p0 * R4 - 75.0 * Y2_p0 * R6 + 41.25 * R8) * V[2] + (37.5 * Y6_p0 - 202.5 * Y4_p0 * R2 + 412.5 * Y2_p0 * R4 - 247.5 * R6) * V[3] + (278.4375 * Y4_p0 - 928.125 * Y2_p0 * R2 + 649.6875 * R4) * V[4] + (649.6875 * Y2_p0 - 649.6875 * R2) * V[5] + 162.421875 * V[6];
}

}  // namespace ovlab