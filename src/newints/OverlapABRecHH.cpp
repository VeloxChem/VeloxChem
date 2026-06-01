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

#include "OverlapABRecHH.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_h_h(
    const CBasisFunction &bra,
    const CBasisFunction &ket,
    const TPoint<double> &bra_center,
    const TPoint<double> &ket_center
) -> newints::Block
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
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt55 = std::sqrt(55.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt77 = std::sqrt(77.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt143 = std::sqrt(143.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt231 = std::sqrt(231.0);
    const auto sqrt286 = std::sqrt(286.0);
    const auto sqrt330 = std::sqrt(330.0);
    const auto sqrt429 = std::sqrt(429.0);
    const auto sqrt462 = std::sqrt(462.0);
    const auto sqrt715 = std::sqrt(715.0);
    const auto sqrt858 = std::sqrt(858.0);
    const auto sqrt1430 = std::sqrt(1430.0);
    const auto sqrt5005 = std::sqrt(5005.0);
    const auto sqrt6006 = std::sqrt(6006.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^5 · β^5 · p^{-10} · (s|s)
    //   V[1] ↔ α^4 · β^4 · p^{-9} · (s|s)
    //   V[2] ↔ α^3 · β^3 · p^{-8} · (s|s)
    //   V[3] ↔ α^2 · β^2 · p^{-7} · (s|s)
    //   V[4] ↔ α · β · p^{-6} · (s|s)
    //   V[5] ↔ p^{-5} · (s|s)
    const auto exps_a  = bra.get_exponents();
    const auto coefs_a = bra.get_normalization_factors();
    const auto exps_b  = ket.get_exponents();
    const auto coefs_b = ket.get_normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 6> V = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < exps_a.size(); ++i)
    {
        const auto alpha = exps_a[i];
        const auto ca    = coefs_a[i];
        const auto alpha2 = alpha * alpha;
        const auto alpha3 = alpha2 * alpha;
        const auto alpha4 = alpha3 * alpha;
        const auto alpha5 = alpha4 * alpha;
        for (std::size_t j = 0; j < exps_b.size(); ++j)
        {
            const auto beta = exps_b[j];
            const auto cb   = coefs_b[j];
            const auto beta2 = beta * beta;
            const auto beta3 = beta2 * beta;
            const auto beta4 = beta3 * beta;
            const auto beta5 = beta4 * beta;

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

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha5 * beta5 * pinv10;
            V[1] += cab_ss * alpha4 * beta4 * pinv9;
            V[2] += cab_ss * alpha3 * beta3 * pinv8;
            V[3] += cab_ss * alpha2 * beta2 * pinv7;
            V[4] += cab_ss * alpha * beta * pinv6;
            V[5] += cab_ss * pinv5;
        }
    }

    // ---- Phase 3: fused M·V → 11 × 11 spherical block ----
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
    const auto Y5_n5 = harm::Y_ll_5_m_n5(AB_x, AB_y, AB_z);
    const auto Y5_n4 = harm::Y_ll_5_m_n4(AB_x, AB_y, AB_z);
    const auto Y5_n3 = harm::Y_ll_5_m_n3(AB_x, AB_y, AB_z);
    const auto Y5_n2 = harm::Y_ll_5_m_n2(AB_x, AB_y, AB_z);
    const auto Y5_n1 = harm::Y_ll_5_m_n1(AB_x, AB_y, AB_z);
    const auto Y5_p0 = harm::Y_ll_5_m_p0(AB_x, AB_y, AB_z);
    const auto Y5_p1 = harm::Y_ll_5_m_p1(AB_x, AB_y, AB_z);
    const auto Y5_p2 = harm::Y_ll_5_m_p2(AB_x, AB_y, AB_z);
    const auto Y5_p3 = harm::Y_ll_5_m_p3(AB_x, AB_y, AB_z);
    const auto Y5_p4 = harm::Y_ll_5_m_p4(AB_x, AB_y, AB_z);
    const auto Y5_p5 = harm::Y_ll_5_m_p5(AB_x, AB_y, AB_z);
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
    const auto R4 = R2 * R2;
    const auto R6 = R4 * R2;
    const auto R8 = R6 * R2;
    newints::Block out{11, 11, std::vector<double>(121, 0.0)};
    auto *d = out.data.data();
    d[60] = -Y5_p0 * Y5_p0 * V[0] + ((245.0 / 143.0) * Y8_p0 + (80.0 / 33.0) * Y6_p0 * R2 + (405.0 / 143.0) * Y4_p0 * R4 + (100.0 / 33.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((-100.0 / 11.0) * Y6_p0 + (-405.0 / 22.0) * Y4_p0 * R2 - 25.0 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (33.75 * Y4_p0 + 75.0 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (-65.625 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[61] = d[71] = -Y5_p0 * Y5_p1 * V[0] + ((49.0 / 143.0) * sqrt15 * Y8_p1 + (8.0 / 33.0) * sqrt35 * Y6_p1 * R2 + (135.0 / 286.0) * sqrt6 * Y4_p1 * R4 + (10.0 / 33.0) * sqrt5 * Y2_p1 * R6) * V[1] + ((-10.0 / 11.0) * sqrt35 * Y6_p1 + (-135.0 / 44.0) * sqrt6 * Y4_p1 * R2 - 2.5 * sqrt5 * Y2_p1 * R4) * V[2] + (5.625 * sqrt6 * Y4_p1 + 7.5 * sqrt5 * Y2_p1 * R2) * V[3] - 6.5625 * sqrt5 * Y2_p1 * V[4];
    d[62] = d[82] = -Y5_p0 * Y5_p2 * V[0] + ((35.0 / 286.0) * sqrt6 * Y8_p2 + (-20.0 / 33.0) * sqrt2 * Y6_p2 * R2 + (-135.0 / 286.0) * sqrt21 * Y4_p2 * R4 + (-20.0 / 33.0) * sqrt35 * Y2_p2 * R6) * V[1] + ((25.0 / 11.0) * sqrt2 * Y6_p2 + (135.0 / 44.0) * sqrt21 * Y4_p2 * R2 + 5.0 * sqrt35 * Y2_p2 * R4) * V[2] + (-5.625 * sqrt21 * Y4_p2 - 15.0 * sqrt35 * Y2_p2 * R2) * V[3] + 13.125 * sqrt35 * Y2_p2 * V[4];
    d[63] = d[93] = -Y5_p0 * Y5_p3 * V[0] + ((-35.0 / 286.0) * sqrt66 * Y8_p3 + (-50.0 / 33.0) * sqrt3 * Y6_p3 * R2 + (-405.0 / 143.0) * Y4_p3 * R4) * V[1] + ((125.0 / 22.0) * sqrt3 * Y6_p3 + (405.0 / 22.0) * Y4_p3 * R2) * V[2] - 33.75 * Y4_p3 * V[3];
    d[64] = d[104] = -Y5_p0 * Y5_p4 * V[0] + ((-7.0 / 26.0) * sqrt55 * Y8_p4 + (-8.0 / 11.0) * sqrt5 * Y6_p4 * R2 + (405.0 / 143.0) * Y4_p4 * R4) * V[1] + ((30.0 / 11.0) * sqrt5 * Y6_p4 + (-405.0 / 22.0) * Y4_p4 * R2) * V[2] + 33.75 * Y4_p4 * V[3];
    d[65] = d[115] = -Y5_p0 * Y5_p5 * V[0] + ((-35.0 / 286.0) * sqrt286 * Y8_p5 + (10.0 / 11.0) * sqrt11 * Y6_p5 * R2) * V[1] + (-75.0 / 22.0) * sqrt11 * Y6_p5 * V[2];
    d[59] = d[49] = -Y5_p0 * Y5_n1 * V[0] + ((49.0 / 143.0) * sqrt15 * Y8_n1 + (8.0 / 33.0) * sqrt35 * Y6_n1 * R2 + (135.0 / 286.0) * sqrt6 * Y4_n1 * R4 + (10.0 / 33.0) * sqrt5 * Y2_n1 * R6) * V[1] + ((-10.0 / 11.0) * sqrt35 * Y6_n1 + (-135.0 / 44.0) * sqrt6 * Y4_n1 * R2 - 2.5 * sqrt5 * Y2_n1 * R4) * V[2] + (5.625 * sqrt6 * Y4_n1 + 7.5 * sqrt5 * Y2_n1 * R2) * V[3] - 6.5625 * sqrt5 * Y2_n1 * V[4];
    d[58] = d[38] = -Y5_p0 * Y5_n2 * V[0] + ((35.0 / 286.0) * sqrt6 * Y8_n2 + (-20.0 / 33.0) * sqrt2 * Y6_n2 * R2 + (-135.0 / 286.0) * sqrt21 * Y4_n2 * R4 + (-20.0 / 33.0) * sqrt35 * Y2_n2 * R6) * V[1] + ((25.0 / 11.0) * sqrt2 * Y6_n2 + (135.0 / 44.0) * sqrt21 * Y4_n2 * R2 + 5.0 * sqrt35 * Y2_n2 * R4) * V[2] + (-5.625 * sqrt21 * Y4_n2 - 15.0 * sqrt35 * Y2_n2 * R2) * V[3] + 13.125 * sqrt35 * Y2_n2 * V[4];
    d[57] = d[27] = -Y5_p0 * Y5_n3 * V[0] + ((-35.0 / 286.0) * sqrt66 * Y8_n3 + (-50.0 / 33.0) * sqrt3 * Y6_n3 * R2 + (-405.0 / 143.0) * Y4_n3 * R4) * V[1] + ((125.0 / 22.0) * sqrt3 * Y6_n3 + (405.0 / 22.0) * Y4_n3 * R2) * V[2] - 33.75 * Y4_n3 * V[3];
    d[56] = d[16] = -Y5_p0 * Y5_n4 * V[0] + ((-7.0 / 26.0) * sqrt55 * Y8_n4 + (-8.0 / 11.0) * sqrt5 * Y6_n4 * R2 + (405.0 / 143.0) * Y4_n4 * R4) * V[1] + ((30.0 / 11.0) * sqrt5 * Y6_n4 + (-405.0 / 22.0) * Y4_n4 * R2) * V[2] + 33.75 * Y4_n4 * V[3];
    d[55] = d[5] = -Y5_p0 * Y5_n5 * V[0] + ((-35.0 / 286.0) * sqrt286 * Y8_n5 + (10.0 / 11.0) * sqrt11 * Y6_n5 * R2) * V[1] + (-75.0 / 22.0) * sqrt11 * Y6_n5 * V[2];
    d[72] = -Y5_p1 * Y5_p1 * V[0] + ((-49.0 / 143.0) * Y8_p0 + (21.0 / 143.0) * sqrt70 * Y8_p2 + (8.0 / 11.0) * Y6_p0 * R2 + (4.0 / 33.0) * sqrt210 * Y6_p2 * R2 + (270.0 / 143.0) * Y4_p0 * R4 + (135.0 / 143.0) * sqrt5 * Y4_p2 * R4 + (30.0 / 11.0) * Y2_p0 * R6 + (50.0 / 33.0) * sqrt3 * Y2_p2 * R6 + 2.5 * R8) * V[1] + ((-30.0 / 11.0) * Y6_p0 + (-5.0 / 11.0) * sqrt210 * Y6_p2 + (-135.0 / 11.0) * Y4_p0 * R2 + (-135.0 / 22.0) * sqrt5 * Y4_p2 * R2 - 22.5 * Y2_p0 * R4 - 12.5 * sqrt3 * Y2_p2 * R4 - 22.5 * R6) * V[2] + (22.5 * Y4_p0 + 11.25 * sqrt5 * Y4_p2 + 67.5 * Y2_p0 * R2 + 37.5 * sqrt3 * Y2_p2 * R2 + 78.75 * R4) * V[3] + (-59.0625 * Y2_p0 - 32.8125 * sqrt3 * Y2_p2 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[73] = d[83] = -Y5_p1 * Y5_p2 * V[0] + ((63.0 / 286.0) * sqrt7 * Y8_p1 + (21.0 / 286.0) * sqrt165 * Y8_p3 + (32.0 / 33.0) * sqrt3 * Y6_p1 * R2 + (2.0 / 11.0) * sqrt30 * Y6_p3 * R2 + (135.0 / 572.0) * sqrt70 * Y4_p1 * R4 + (135.0 / 572.0) * sqrt10 * Y4_p3 * R4 + (10.0 / 33.0) * sqrt21 * Y2_p1 * R6) * V[1] + ((-40.0 / 11.0) * sqrt3 * Y6_p1 + (-15.0 / 22.0) * sqrt30 * Y6_p3 + (-135.0 / 88.0) * sqrt70 * Y4_p1 * R2 + (-135.0 / 88.0) * sqrt10 * Y4_p3 * R2 - 2.5 * sqrt21 * Y2_p1 * R4) * V[2] + (2.8125 * sqrt70 * Y4_p1 + 2.8125 * sqrt10 * Y4_p3 + 7.5 * sqrt21 * Y2_p1 * R2) * V[3] - 6.5625 * sqrt21 * Y2_p1 * V[4];
    d[74] = d[94] = -Y5_p1 * Y5_p3 * V[0] + ((7.0 / 22.0) * sqrt15 * Y8_p2 + (7.0 / 572.0) * sqrt66 * Y8_p4 + (2.0 / 3.0) * sqrt5 * Y6_p2 * R2 + (-13.0 / 33.0) * sqrt6 * Y6_p4 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_p4 * R4 + (-20.0 / 33.0) * sqrt14 * Y2_p2 * R6) * V[1] + (-2.5 * sqrt5 * Y6_p2 + (65.0 / 44.0) * sqrt6 * Y6_p4 + (135.0 / 44.0) * sqrt30 * Y4_p4 * R2 + 5.0 * sqrt14 * Y2_p2 * R4) * V[2] + (-5.625 * sqrt30 * Y4_p4 - 15.0 * sqrt14 * Y2_p2 * R2) * V[3] + 13.125 * sqrt14 * Y2_p2 * V[4];
    d[75] = d[105] = -Y5_p1 * Y5_p4 * V[0] + ((28.0 / 143.0) * sqrt55 * Y8_p3 + (-7.0 / 143.0) * sqrt429 * Y8_p5 + (-1.0 / 11.0) * sqrt10 * Y6_p3 * R2 + (-3.0 / 11.0) * sqrt66 * Y6_p5 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_p3 * R4) * V[1] + ((15.0 / 44.0) * sqrt10 * Y6_p3 + (45.0 / 44.0) * sqrt66 * Y6_p5 + (135.0 / 44.0) * sqrt30 * Y4_p3 * R2) * V[2] - 5.625 * sqrt30 * Y4_p3 * V[3];
    d[76] = d[116] = -Y5_p1 * Y5_p5 * V[0] + ((35.0 / 572.0) * sqrt330 * Y8_p4 + (-7.0 / 286.0) * sqrt5005 * Y8_p6 + (-5.0 / 11.0) * sqrt30 * Y6_p4 * R2 + (2.0 / 11.0) * sqrt55 * Y6_p6 * R2 + (135.0 / 286.0) * sqrt6 * Y4_p4 * R4) * V[1] + ((75.0 / 44.0) * sqrt30 * Y6_p4 + (-15.0 / 22.0) * sqrt55 * Y6_p6 + (-135.0 / 44.0) * sqrt6 * Y4_p4 * R2) * V[2] + 5.625 * sqrt6 * Y4_p4 * V[3];
    d[70] = d[50] = -Y5_p1 * Y5_n1 * V[0] + ((21.0 / 143.0) * sqrt70 * Y8_n2 + (4.0 / 33.0) * sqrt210 * Y6_n2 * R2 + (135.0 / 143.0) * sqrt5 * Y4_n2 * R4 + (50.0 / 33.0) * sqrt3 * Y2_n2 * R6) * V[1] + ((-5.0 / 11.0) * sqrt210 * Y6_n2 + (-135.0 / 22.0) * sqrt5 * Y4_n2 * R2 - 12.5 * sqrt3 * Y2_n2 * R4) * V[2] + (11.25 * sqrt5 * Y4_n2 + 37.5 * sqrt3 * Y2_n2 * R2) * V[3] - 32.8125 * sqrt3 * Y2_n2 * V[4];
    d[69] = d[39] = -Y5_p1 * Y5_n2 * V[0] + ((21.0 / 286.0) * sqrt165 * Y8_n3 + (63.0 / 286.0) * sqrt7 * Y8_n1 + (2.0 / 11.0) * sqrt30 * Y6_n3 * R2 + (32.0 / 33.0) * sqrt3 * Y6_n1 * R2 + (135.0 / 572.0) * sqrt10 * Y4_n3 * R4 + (135.0 / 572.0) * sqrt70 * Y4_n1 * R4 + (10.0 / 33.0) * sqrt21 * Y2_n1 * R6) * V[1] + ((-15.0 / 22.0) * sqrt30 * Y6_n3 + (-40.0 / 11.0) * sqrt3 * Y6_n1 + (-135.0 / 88.0) * sqrt10 * Y4_n3 * R2 + (-135.0 / 88.0) * sqrt70 * Y4_n1 * R2 - 2.5 * sqrt21 * Y2_n1 * R4) * V[2] + (2.8125 * sqrt10 * Y4_n3 + 2.8125 * sqrt70 * Y4_n1 + 7.5 * sqrt21 * Y2_n1 * R2) * V[3] - 6.5625 * sqrt21 * Y2_n1 * V[4];
    d[68] = d[28] = -Y5_p1 * Y5_n3 * V[0] + ((7.0 / 572.0) * sqrt66 * Y8_n4 + (7.0 / 22.0) * sqrt15 * Y8_n2 + (-13.0 / 33.0) * sqrt6 * Y6_n4 * R2 + (2.0 / 3.0) * sqrt5 * Y6_n2 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_n4 * R4 + (-20.0 / 33.0) * sqrt14 * Y2_n2 * R6) * V[1] + ((65.0 / 44.0) * sqrt6 * Y6_n4 - 2.5 * sqrt5 * Y6_n2 + (135.0 / 44.0) * sqrt30 * Y4_n4 * R2 + 5.0 * sqrt14 * Y2_n2 * R4) * V[2] + (-5.625 * sqrt30 * Y4_n4 - 15.0 * sqrt14 * Y2_n2 * R2) * V[3] + 13.125 * sqrt14 * Y2_n2 * V[4];
    d[67] = d[17] = -Y5_p1 * Y5_n4 * V[0] + ((-7.0 / 143.0) * sqrt429 * Y8_n5 + (28.0 / 143.0) * sqrt55 * Y8_n3 + (-3.0 / 11.0) * sqrt66 * Y6_n5 * R2 + (-1.0 / 11.0) * sqrt10 * Y6_n3 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_n3 * R4) * V[1] + ((45.0 / 44.0) * sqrt66 * Y6_n5 + (15.0 / 44.0) * sqrt10 * Y6_n3 + (135.0 / 44.0) * sqrt30 * Y4_n3 * R2) * V[2] - 5.625 * sqrt30 * Y4_n3 * V[3];
    d[66] = d[6] = -Y5_p1 * Y5_n5 * V[0] + ((-7.0 / 286.0) * sqrt5005 * Y8_n6 + (35.0 / 572.0) * sqrt330 * Y8_n4 + (2.0 / 11.0) * sqrt55 * Y6_n6 * R2 + (-5.0 / 11.0) * sqrt30 * Y6_n4 * R2 + (135.0 / 286.0) * sqrt6 * Y4_n4 * R4) * V[1] + ((-15.0 / 22.0) * sqrt55 * Y6_n6 + (75.0 / 44.0) * sqrt30 * Y6_n4 + (-135.0 / 44.0) * sqrt6 * Y4_n4 * R2) * V[2] + 5.625 * sqrt6 * Y4_n4 * V[3];
    d[84] = -Y5_p2 * Y5_p2 * V[0] + ((-238.0 / 143.0) * Y8_p0 + (21.0 / 143.0) * sqrt77 * Y8_p4 + (-24.0 / 11.0) * Y6_p0 * R2 + (8.0 / 11.0) * sqrt7 * Y6_p4 * R2 + (-135.0 / 286.0) * Y4_p0 * R4 + (135.0 / 286.0) * sqrt35 * Y4_p4 * R4 + (20.0 / 11.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((90.0 / 11.0) * Y6_p0 + (-30.0 / 11.0) * sqrt7 * Y6_p4 + (135.0 / 44.0) * Y4_p0 * R2 + (-135.0 / 44.0) * sqrt35 * Y4_p4 * R2 - 15.0 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (-5.625 * Y4_p0 + 5.625 * sqrt35 * Y4_p4 + 45.0 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (-39.375 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[85] = d[95] = -Y5_p2 * Y5_p3 * V[0] + ((-329.0 / 572.0) * sqrt6 * Y8_p1 + (7.0 / 572.0) * sqrt6006 * Y8_p5 + (-2.0 / 33.0) * sqrt14 * Y6_p1 * R2 + (2.0 / 33.0) * sqrt231 * Y6_p5 * R2 + (135.0 / 286.0) * sqrt15 * Y4_p1 * R4 + (50.0 / 33.0) * sqrt2 * Y2_p1 * R6) * V[1] + ((5.0 / 22.0) * sqrt14 * Y6_p1 + (-5.0 / 22.0) * sqrt231 * Y6_p5 + (-135.0 / 44.0) * sqrt15 * Y4_p1 * R2 - 12.5 * sqrt2 * Y2_p1 * R4) * V[2] + (5.625 * sqrt15 * Y4_p1 + 37.5 * sqrt2 * Y2_p1 * R2) * V[3] - 32.8125 * sqrt2 * Y2_p1 * V[4];
    d[86] = d[106] = -Y5_p2 * Y5_p4 * V[0] + ((-49.0 / 572.0) * sqrt210 * Y8_p2 + (-7.0 / 572.0) * sqrt286 * Y8_p6 + (2.0 / 11.0) * sqrt70 * Y6_p2 * R2 + (-2.0 / 11.0) * sqrt154 * Y6_p6 * R2 + (135.0 / 286.0) * sqrt15 * Y4_p2 * R4 + (-20.0 / 11.0) * Y2_p2 * R6) * V[1] + ((-15.0 / 22.0) * sqrt70 * Y6_p2 + (15.0 / 22.0) * sqrt154 * Y6_p6 + (-135.0 / 44.0) * sqrt15 * Y4_p2 * R2 + 15.0 * Y2_p2 * R4) * V[2] + (5.625 * sqrt15 * Y4_p2 - 45.0 * Y2_p2 * R2) * V[3] + 39.375 * Y2_p2 * V[4];
    d[87] = d[117] = -Y5_p2 * Y5_p5 * V[0] + ((-35.0 / 572.0) * sqrt154 * Y8_p3 + (-35.0 / 572.0) * sqrt858 * Y8_p7 + (10.0 / 11.0) * sqrt7 * Y6_p3 * R2 + (-135.0 / 286.0) * sqrt21 * Y4_p3 * R4) * V[1] + ((-75.0 / 22.0) * sqrt7 * Y6_p3 + (135.0 / 44.0) * sqrt21 * Y4_p3 * R2) * V[2] - 5.625 * sqrt21 * Y4_p3 * V[3];
    d[81] = d[51] = -Y5_p2 * Y5_n1 * V[0] + ((21.0 / 286.0) * sqrt165 * Y8_n3 + (-63.0 / 286.0) * sqrt7 * Y8_n1 + (2.0 / 11.0) * sqrt30 * Y6_n3 * R2 + (-32.0 / 33.0) * sqrt3 * Y6_n1 * R2 + (135.0 / 572.0) * sqrt10 * Y4_n3 * R4 + (-135.0 / 572.0) * sqrt70 * Y4_n1 * R4 + (-10.0 / 33.0) * sqrt21 * Y2_n1 * R6) * V[1] + ((-15.0 / 22.0) * sqrt30 * Y6_n3 + (40.0 / 11.0) * sqrt3 * Y6_n1 + (-135.0 / 88.0) * sqrt10 * Y4_n3 * R2 + (135.0 / 88.0) * sqrt70 * Y4_n1 * R2 + 2.5 * sqrt21 * Y2_n1 * R4) * V[2] + (2.8125 * sqrt10 * Y4_n3 - 2.8125 * sqrt70 * Y4_n1 - 7.5 * sqrt21 * Y2_n1 * R2) * V[3] + 6.5625 * sqrt21 * Y2_n1 * V[4];
    d[80] = d[40] = -Y5_p2 * Y5_n2 * V[0] + ((21.0 / 143.0) * sqrt77 * Y8_n4 + (8.0 / 11.0) * sqrt7 * Y6_n4 * R2 + (135.0 / 286.0) * sqrt35 * Y4_n4 * R4) * V[1] + ((-30.0 / 11.0) * sqrt7 * Y6_n4 + (-135.0 / 44.0) * sqrt35 * Y4_n4 * R2) * V[2] + 5.625 * sqrt35 * Y4_n4 * V[3];
    d[79] = d[29] = -Y5_p2 * Y5_n3 * V[0] + ((7.0 / 572.0) * sqrt6006 * Y8_n5 + (-329.0 / 572.0) * sqrt6 * Y8_n1 + (2.0 / 33.0) * sqrt231 * Y6_n5 * R2 + (-2.0 / 33.0) * sqrt14 * Y6_n1 * R2 + (135.0 / 286.0) * sqrt15 * Y4_n1 * R4 + (50.0 / 33.0) * sqrt2 * Y2_n1 * R6) * V[1] + ((-5.0 / 22.0) * sqrt231 * Y6_n5 + (5.0 / 22.0) * sqrt14 * Y6_n1 + (-135.0 / 44.0) * sqrt15 * Y4_n1 * R2 - 12.5 * sqrt2 * Y2_n1 * R4) * V[2] + (5.625 * sqrt15 * Y4_n1 + 37.5 * sqrt2 * Y2_n1 * R2) * V[3] - 32.8125 * sqrt2 * Y2_n1 * V[4];
    d[78] = d[18] = -Y5_p2 * Y5_n4 * V[0] + ((-7.0 / 572.0) * sqrt286 * Y8_n6 + (-49.0 / 572.0) * sqrt210 * Y8_n2 + (-2.0 / 11.0) * sqrt154 * Y6_n6 * R2 + (2.0 / 11.0) * sqrt70 * Y6_n2 * R2 + (135.0 / 286.0) * sqrt15 * Y4_n2 * R4 + (-20.0 / 11.0) * Y2_n2 * R6) * V[1] + ((15.0 / 22.0) * sqrt154 * Y6_n6 + (-15.0 / 22.0) * sqrt70 * Y6_n2 + (-135.0 / 44.0) * sqrt15 * Y4_n2 * R2 + 15.0 * Y2_n2 * R4) * V[2] + (5.625 * sqrt15 * Y4_n2 - 45.0 * Y2_n2 * R2) * V[3] + 39.375 * Y2_n2 * V[4];
    d[77] = d[7] = -Y5_p2 * Y5_n5 * V[0] + ((-35.0 / 572.0) * sqrt858 * Y8_n7 + (-35.0 / 572.0) * sqrt154 * Y8_n3 + (10.0 / 11.0) * sqrt7 * Y6_n3 * R2 + (-135.0 / 286.0) * sqrt21 * Y4_n3 * R4) * V[1] + ((-75.0 / 22.0) * sqrt7 * Y6_n3 + (135.0 / 44.0) * sqrt21 * Y4_n3 * R2) * V[2] - 5.625 * sqrt21 * Y4_n3 * V[3];
    d[96] = -Y5_p3 * Y5_p3 * V[0] + ((511.0 / 286.0) * Y8_p0 + (7.0 / 143.0) * sqrt858 * Y8_p6 + (-58.0 / 33.0) * Y6_p0 * R2 + (4.0 / 33.0) * sqrt462 * Y6_p6 * R2 + (-405.0 / 143.0) * Y4_p0 * R4 + (10.0 / 33.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((145.0 / 22.0) * Y6_p0 + (-5.0 / 11.0) * sqrt462 * Y6_p6 + (405.0 / 22.0) * Y4_p0 * R2 - 2.5 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (-33.75 * Y4_p0 + 7.5 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (-6.5625 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[97] = d[107] = -Y5_p3 * Y5_p4 * V[0] + ((7.0 / 11.0) * sqrt2 * Y8_p1 + (7.0 / 286.0) * sqrt1430 * Y8_p7 + (-1.0 / 3.0) * sqrt42 * Y6_p1 * R2 + (35.0 / 33.0) * sqrt6 * Y2_p1 * R6) * V[1] + (1.25 * sqrt42 * Y6_p1 - 8.75 * sqrt6 * Y2_p1 * R4) * V[2] + 26.25 * sqrt6 * Y2_p1 * R2 * V[3] - 22.96875 * sqrt6 * Y2_p1 * V[4];
    d[98] = d[118] = -Y5_p3 * Y5_p5 * V[0] + ((35.0 / 286.0) * sqrt14 * Y8_p2 + (-35.0 / 286.0) * sqrt143 * Y8_p8 + (-10.0 / 33.0) * sqrt42 * Y6_p2 * R2 + (405.0 / 143.0) * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt15 * Y2_p2 * R6) * V[1] + ((25.0 / 22.0) * sqrt42 * Y6_p2 + (-405.0 / 22.0) * Y4_p2 * R2 + 2.5 * sqrt15 * Y2_p2 * R4) * V[2] + (33.75 * Y4_p2 - 7.5 * sqrt15 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt15 * Y2_p2 * V[4];
    d[92] = d[52] = -Y5_p3 * Y5_n1 * V[0] + ((7.0 / 572.0) * sqrt66 * Y8_n4 + (-7.0 / 22.0) * sqrt15 * Y8_n2 + (-13.0 / 33.0) * sqrt6 * Y6_n4 * R2 + (-2.0 / 3.0) * sqrt5 * Y6_n2 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_n4 * R4 + (20.0 / 33.0) * sqrt14 * Y2_n2 * R6) * V[1] + ((65.0 / 44.0) * sqrt6 * Y6_n4 + 2.5 * sqrt5 * Y6_n2 + (135.0 / 44.0) * sqrt30 * Y4_n4 * R2 - 5.0 * sqrt14 * Y2_n2 * R4) * V[2] + (-5.625 * sqrt30 * Y4_n4 + 15.0 * sqrt14 * Y2_n2 * R2) * V[3] - 13.125 * sqrt14 * Y2_n2 * V[4];
    d[91] = d[41] = -Y5_p3 * Y5_n2 * V[0] + ((7.0 / 572.0) * sqrt6006 * Y8_n5 + (329.0 / 572.0) * sqrt6 * Y8_n1 + (2.0 / 33.0) * sqrt231 * Y6_n5 * R2 + (2.0 / 33.0) * sqrt14 * Y6_n1 * R2 + (-135.0 / 286.0) * sqrt15 * Y4_n1 * R4 + (-50.0 / 33.0) * sqrt2 * Y2_n1 * R6) * V[1] + ((-5.0 / 22.0) * sqrt231 * Y6_n5 + (-5.0 / 22.0) * sqrt14 * Y6_n1 + (135.0 / 44.0) * sqrt15 * Y4_n1 * R2 + 12.5 * sqrt2 * Y2_n1 * R4) * V[2] + (-5.625 * sqrt15 * Y4_n1 - 37.5 * sqrt2 * Y2_n1 * R2) * V[3] + 32.8125 * sqrt2 * Y2_n1 * V[4];
    d[90] = d[30] = -Y5_p3 * Y5_n3 * V[0] + ((7.0 / 143.0) * sqrt858 * Y8_n6 + (4.0 / 33.0) * sqrt462 * Y6_n6 * R2) * V[1] + (-5.0 / 11.0) * sqrt462 * Y6_n6 * V[2];
    d[89] = d[19] = -Y5_p3 * Y5_n4 * V[0] + ((7.0 / 286.0) * sqrt1430 * Y8_n7 + (7.0 / 11.0) * sqrt2 * Y8_n1 + (-1.0 / 3.0) * sqrt42 * Y6_n1 * R2 + (35.0 / 33.0) * sqrt6 * Y2_n1 * R6) * V[1] + (1.25 * sqrt42 * Y6_n1 - 8.75 * sqrt6 * Y2_n1 * R4) * V[2] + 26.25 * sqrt6 * Y2_n1 * R2 * V[3] - 22.96875 * sqrt6 * Y2_n1 * V[4];
    d[88] = d[8] = -Y5_p3 * Y5_n5 * V[0] + ((-35.0 / 286.0) * sqrt143 * Y8_n8 + (35.0 / 286.0) * sqrt14 * Y8_n2 + (-10.0 / 33.0) * sqrt42 * Y6_n2 * R2 + (405.0 / 143.0) * Y4_n2 * R4 + (-10.0 / 33.0) * sqrt15 * Y2_n2 * R6) * V[1] + ((25.0 / 22.0) * sqrt42 * Y6_n2 + (-405.0 / 22.0) * Y4_n2 * R2 + 2.5 * sqrt15 * Y2_n2 * R4) * V[2] + (33.75 * Y4_n2 - 7.5 * sqrt15 * Y2_n2 * R2) * V[3] + 6.5625 * sqrt15 * Y2_n2 * V[4];
    d[108] = -Y5_p4 * Y5_p4 * V[0] + ((-217.0 / 286.0) * Y8_p0 + (21.0 / 286.0) * sqrt715 * Y8_p8 + (32.0 / 11.0) * Y6_p0 * R2 + (-405.0 / 143.0) * Y4_p0 * R4 + (-20.0 / 11.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((-120.0 / 11.0) * Y6_p0 + (405.0 / 22.0) * Y4_p0 * R2 + 15.0 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (-33.75 * Y4_p0 - 45.0 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (39.375 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[109] = d[119] = -Y5_p4 * Y5_p5 * V[0] + ((-21.0 / 286.0) * sqrt10 * Y8_p1 + (1.0 / 11.0) * sqrt210 * Y6_p1 * R2 + (-405.0 / 143.0) * Y4_p1 * R4 + (5.0 / 11.0) * sqrt30 * Y2_p1 * R6) * V[1] + ((-15.0 / 44.0) * sqrt210 * Y6_p1 + (405.0 / 22.0) * Y4_p1 * R2 - 3.75 * sqrt30 * Y2_p1 * R4) * V[2] + (-33.75 * Y4_p1 + 11.25 * sqrt30 * Y2_p1 * R2) * V[3] - 9.84375 * sqrt30 * Y2_p1 * V[4];
    d[103] = d[53] = -Y5_p4 * Y5_n1 * V[0] + ((-7.0 / 143.0) * sqrt429 * Y8_n5 + (-28.0 / 143.0) * sqrt55 * Y8_n3 + (-3.0 / 11.0) * sqrt66 * Y6_n5 * R2 + (1.0 / 11.0) * sqrt10 * Y6_n3 * R2 + (135.0 / 286.0) * sqrt30 * Y4_n3 * R4) * V[1] + ((45.0 / 44.0) * sqrt66 * Y6_n5 + (-15.0 / 44.0) * sqrt10 * Y6_n3 + (-135.0 / 44.0) * sqrt30 * Y4_n3 * R2) * V[2] + 5.625 * sqrt30 * Y4_n3 * V[3];
    d[102] = d[42] = -Y5_p4 * Y5_n2 * V[0] + ((-7.0 / 572.0) * sqrt286 * Y8_n6 + (49.0 / 572.0) * sqrt210 * Y8_n2 + (-2.0 / 11.0) * sqrt154 * Y6_n6 * R2 + (-2.0 / 11.0) * sqrt70 * Y6_n2 * R2 + (-135.0 / 286.0) * sqrt15 * Y4_n2 * R4 + (20.0 / 11.0) * Y2_n2 * R6) * V[1] + ((15.0 / 22.0) * sqrt154 * Y6_n6 + (15.0 / 22.0) * sqrt70 * Y6_n2 + (135.0 / 44.0) * sqrt15 * Y4_n2 * R2 - 15.0 * Y2_n2 * R4) * V[2] + (-5.625 * sqrt15 * Y4_n2 + 45.0 * Y2_n2 * R2) * V[3] - 39.375 * Y2_n2 * V[4];
    d[101] = d[31] = -Y5_p4 * Y5_n3 * V[0] + ((7.0 / 286.0) * sqrt1430 * Y8_n7 + (-7.0 / 11.0) * sqrt2 * Y8_n1 + (1.0 / 3.0) * sqrt42 * Y6_n1 * R2 + (-35.0 / 33.0) * sqrt6 * Y2_n1 * R6) * V[1] + (-1.25 * sqrt42 * Y6_n1 + 8.75 * sqrt6 * Y2_n1 * R4) * V[2] - 26.25 * sqrt6 * Y2_n1 * R2 * V[3] + 22.96875 * sqrt6 * Y2_n1 * V[4];
    d[100] = d[20] = -Y5_p4 * Y5_n4 * V[0] + (21.0 / 286.0) * sqrt715 * Y8_n8 * V[1];
    d[99] = d[9] = -Y5_p4 * Y5_n5 * V[0] + ((-21.0 / 286.0) * sqrt10 * Y8_n1 + (1.0 / 11.0) * sqrt210 * Y6_n1 * R2 + (-405.0 / 143.0) * Y4_n1 * R4 + (5.0 / 11.0) * sqrt30 * Y2_n1 * R6) * V[1] + ((-15.0 / 44.0) * sqrt210 * Y6_n1 + (405.0 / 22.0) * Y4_n1 * R2 - 3.75 * sqrt30 * Y2_n1 * R4) * V[2] + (-33.75 * Y4_n1 + 11.25 * sqrt30 * Y2_n1 * R2) * V[3] - 9.84375 * sqrt30 * Y2_n1 * V[4];
    d[120] = -Y5_p5 * Y5_p5 * V[0] + ((35.0 / 286.0) * Y8_p0 + (-10.0 / 11.0) * Y6_p0 * R2 + (405.0 / 143.0) * Y4_p0 * R4 + (-50.0 / 11.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((75.0 / 22.0) * Y6_p0 + (-405.0 / 22.0) * Y4_p0 * R2 + 37.5 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (33.75 * Y4_p0 - 112.5 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (98.4375 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[114] = d[54] = -Y5_p5 * Y5_n1 * V[0] + ((-7.0 / 286.0) * sqrt5005 * Y8_n6 + (-35.0 / 572.0) * sqrt330 * Y8_n4 + (2.0 / 11.0) * sqrt55 * Y6_n6 * R2 + (5.0 / 11.0) * sqrt30 * Y6_n4 * R2 + (-135.0 / 286.0) * sqrt6 * Y4_n4 * R4) * V[1] + ((-15.0 / 22.0) * sqrt55 * Y6_n6 + (-75.0 / 44.0) * sqrt30 * Y6_n4 + (135.0 / 44.0) * sqrt6 * Y4_n4 * R2) * V[2] - 5.625 * sqrt6 * Y4_n4 * V[3];
    d[113] = d[43] = -Y5_p5 * Y5_n2 * V[0] + ((-35.0 / 572.0) * sqrt858 * Y8_n7 + (35.0 / 572.0) * sqrt154 * Y8_n3 + (-10.0 / 11.0) * sqrt7 * Y6_n3 * R2 + (135.0 / 286.0) * sqrt21 * Y4_n3 * R4) * V[1] + ((75.0 / 22.0) * sqrt7 * Y6_n3 + (-135.0 / 44.0) * sqrt21 * Y4_n3 * R2) * V[2] + 5.625 * sqrt21 * Y4_n3 * V[3];
    d[112] = d[32] = -Y5_p5 * Y5_n3 * V[0] + ((-35.0 / 286.0) * sqrt143 * Y8_n8 + (-35.0 / 286.0) * sqrt14 * Y8_n2 + (10.0 / 33.0) * sqrt42 * Y6_n2 * R2 + (-405.0 / 143.0) * Y4_n2 * R4 + (10.0 / 33.0) * sqrt15 * Y2_n2 * R6) * V[1] + ((-25.0 / 22.0) * sqrt42 * Y6_n2 + (405.0 / 22.0) * Y4_n2 * R2 - 2.5 * sqrt15 * Y2_n2 * R4) * V[2] + (-33.75 * Y4_n2 + 7.5 * sqrt15 * Y2_n2 * R2) * V[3] - 6.5625 * sqrt15 * Y2_n2 * V[4];
    d[111] = d[21] = -Y5_p5 * Y5_n4 * V[0] + ((21.0 / 286.0) * sqrt10 * Y8_n1 + (-1.0 / 11.0) * sqrt210 * Y6_n1 * R2 + (405.0 / 143.0) * Y4_n1 * R4 + (-5.0 / 11.0) * sqrt30 * Y2_n1 * R6) * V[1] + ((15.0 / 44.0) * sqrt210 * Y6_n1 + (-405.0 / 22.0) * Y4_n1 * R2 + 3.75 * sqrt30 * Y2_n1 * R4) * V[2] + (33.75 * Y4_n1 - 11.25 * sqrt30 * Y2_n1 * R2) * V[3] + 9.84375 * sqrt30 * Y2_n1 * V[4];
    d[110] = d[10] = -Y5_p5 * Y5_n5 * V[0];
    d[48] = -Y5_n1 * Y5_n1 * V[0] + ((-49.0 / 143.0) * Y8_p0 + (-21.0 / 143.0) * sqrt70 * Y8_p2 + (8.0 / 11.0) * Y6_p0 * R2 + (-4.0 / 33.0) * sqrt210 * Y6_p2 * R2 + (270.0 / 143.0) * Y4_p0 * R4 + (-135.0 / 143.0) * sqrt5 * Y4_p2 * R4 + (30.0 / 11.0) * Y2_p0 * R6 + (-50.0 / 33.0) * sqrt3 * Y2_p2 * R6 + 2.5 * R8) * V[1] + ((-30.0 / 11.0) * Y6_p0 + (5.0 / 11.0) * sqrt210 * Y6_p2 + (-135.0 / 11.0) * Y4_p0 * R2 + (135.0 / 22.0) * sqrt5 * Y4_p2 * R2 - 22.5 * Y2_p0 * R4 + 12.5 * sqrt3 * Y2_p2 * R4 - 22.5 * R6) * V[2] + (22.5 * Y4_p0 - 11.25 * sqrt5 * Y4_p2 + 67.5 * Y2_p0 * R2 - 37.5 * sqrt3 * Y2_p2 * R2 + 78.75 * R4) * V[3] + (-59.0625 * Y2_p0 + 32.8125 * sqrt3 * Y2_p2 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[47] = d[37] = -Y5_n1 * Y5_n2 * V[0] + ((63.0 / 286.0) * sqrt7 * Y8_p1 + (-21.0 / 286.0) * sqrt165 * Y8_p3 + (32.0 / 33.0) * sqrt3 * Y6_p1 * R2 + (-2.0 / 11.0) * sqrt30 * Y6_p3 * R2 + (135.0 / 572.0) * sqrt70 * Y4_p1 * R4 + (-135.0 / 572.0) * sqrt10 * Y4_p3 * R4 + (10.0 / 33.0) * sqrt21 * Y2_p1 * R6) * V[1] + ((-40.0 / 11.0) * sqrt3 * Y6_p1 + (15.0 / 22.0) * sqrt30 * Y6_p3 + (-135.0 / 88.0) * sqrt70 * Y4_p1 * R2 + (135.0 / 88.0) * sqrt10 * Y4_p3 * R2 - 2.5 * sqrt21 * Y2_p1 * R4) * V[2] + (2.8125 * sqrt70 * Y4_p1 - 2.8125 * sqrt10 * Y4_p3 + 7.5 * sqrt21 * Y2_p1 * R2) * V[3] - 6.5625 * sqrt21 * Y2_p1 * V[4];
    d[46] = d[26] = -Y5_n1 * Y5_n3 * V[0] + ((7.0 / 22.0) * sqrt15 * Y8_p2 + (-7.0 / 572.0) * sqrt66 * Y8_p4 + (2.0 / 3.0) * sqrt5 * Y6_p2 * R2 + (13.0 / 33.0) * sqrt6 * Y6_p4 * R2 + (135.0 / 286.0) * sqrt30 * Y4_p4 * R4 + (-20.0 / 33.0) * sqrt14 * Y2_p2 * R6) * V[1] + (-2.5 * sqrt5 * Y6_p2 + (-65.0 / 44.0) * sqrt6 * Y6_p4 + (-135.0 / 44.0) * sqrt30 * Y4_p4 * R2 + 5.0 * sqrt14 * Y2_p2 * R4) * V[2] + (5.625 * sqrt30 * Y4_p4 - 15.0 * sqrt14 * Y2_p2 * R2) * V[3] + 13.125 * sqrt14 * Y2_p2 * V[4];
    d[45] = d[15] = -Y5_n1 * Y5_n4 * V[0] + ((28.0 / 143.0) * sqrt55 * Y8_p3 + (7.0 / 143.0) * sqrt429 * Y8_p5 + (-1.0 / 11.0) * sqrt10 * Y6_p3 * R2 + (3.0 / 11.0) * sqrt66 * Y6_p5 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_p3 * R4) * V[1] + ((15.0 / 44.0) * sqrt10 * Y6_p3 + (-45.0 / 44.0) * sqrt66 * Y6_p5 + (135.0 / 44.0) * sqrt30 * Y4_p3 * R2) * V[2] - 5.625 * sqrt30 * Y4_p3 * V[3];
    d[44] = d[4] = -Y5_n1 * Y5_n5 * V[0] + ((35.0 / 572.0) * sqrt330 * Y8_p4 + (7.0 / 286.0) * sqrt5005 * Y8_p6 + (-5.0 / 11.0) * sqrt30 * Y6_p4 * R2 + (-2.0 / 11.0) * sqrt55 * Y6_p6 * R2 + (135.0 / 286.0) * sqrt6 * Y4_p4 * R4) * V[1] + ((75.0 / 44.0) * sqrt30 * Y6_p4 + (15.0 / 22.0) * sqrt55 * Y6_p6 + (-135.0 / 44.0) * sqrt6 * Y4_p4 * R2) * V[2] + 5.625 * sqrt6 * Y4_p4 * V[3];
    d[36] = -Y5_n2 * Y5_n2 * V[0] + ((-238.0 / 143.0) * Y8_p0 + (-21.0 / 143.0) * sqrt77 * Y8_p4 + (-24.0 / 11.0) * Y6_p0 * R2 + (-8.0 / 11.0) * sqrt7 * Y6_p4 * R2 + (-135.0 / 286.0) * Y4_p0 * R4 + (-135.0 / 286.0) * sqrt35 * Y4_p4 * R4 + (20.0 / 11.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((90.0 / 11.0) * Y6_p0 + (30.0 / 11.0) * sqrt7 * Y6_p4 + (135.0 / 44.0) * Y4_p0 * R2 + (135.0 / 44.0) * sqrt35 * Y4_p4 * R2 - 15.0 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (-5.625 * Y4_p0 - 5.625 * sqrt35 * Y4_p4 + 45.0 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (-39.375 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[35] = d[25] = -Y5_n2 * Y5_n3 * V[0] + ((-329.0 / 572.0) * sqrt6 * Y8_p1 + (-7.0 / 572.0) * sqrt6006 * Y8_p5 + (-2.0 / 33.0) * sqrt14 * Y6_p1 * R2 + (-2.0 / 33.0) * sqrt231 * Y6_p5 * R2 + (135.0 / 286.0) * sqrt15 * Y4_p1 * R4 + (50.0 / 33.0) * sqrt2 * Y2_p1 * R6) * V[1] + ((5.0 / 22.0) * sqrt14 * Y6_p1 + (5.0 / 22.0) * sqrt231 * Y6_p5 + (-135.0 / 44.0) * sqrt15 * Y4_p1 * R2 - 12.5 * sqrt2 * Y2_p1 * R4) * V[2] + (5.625 * sqrt15 * Y4_p1 + 37.5 * sqrt2 * Y2_p1 * R2) * V[3] - 32.8125 * sqrt2 * Y2_p1 * V[4];
    d[34] = d[14] = -Y5_n2 * Y5_n4 * V[0] + ((-49.0 / 572.0) * sqrt210 * Y8_p2 + (7.0 / 572.0) * sqrt286 * Y8_p6 + (2.0 / 11.0) * sqrt70 * Y6_p2 * R2 + (2.0 / 11.0) * sqrt154 * Y6_p6 * R2 + (135.0 / 286.0) * sqrt15 * Y4_p2 * R4 + (-20.0 / 11.0) * Y2_p2 * R6) * V[1] + ((-15.0 / 22.0) * sqrt70 * Y6_p2 + (-15.0 / 22.0) * sqrt154 * Y6_p6 + (-135.0 / 44.0) * sqrt15 * Y4_p2 * R2 + 15.0 * Y2_p2 * R4) * V[2] + (5.625 * sqrt15 * Y4_p2 - 45.0 * Y2_p2 * R2) * V[3] + 39.375 * Y2_p2 * V[4];
    d[33] = d[3] = -Y5_n2 * Y5_n5 * V[0] + ((-35.0 / 572.0) * sqrt154 * Y8_p3 + (35.0 / 572.0) * sqrt858 * Y8_p7 + (10.0 / 11.0) * sqrt7 * Y6_p3 * R2 + (-135.0 / 286.0) * sqrt21 * Y4_p3 * R4) * V[1] + ((-75.0 / 22.0) * sqrt7 * Y6_p3 + (135.0 / 44.0) * sqrt21 * Y4_p3 * R2) * V[2] - 5.625 * sqrt21 * Y4_p3 * V[3];
    d[24] = -Y5_n3 * Y5_n3 * V[0] + ((511.0 / 286.0) * Y8_p0 + (-7.0 / 143.0) * sqrt858 * Y8_p6 + (-58.0 / 33.0) * Y6_p0 * R2 + (-4.0 / 33.0) * sqrt462 * Y6_p6 * R2 + (-405.0 / 143.0) * Y4_p0 * R4 + (10.0 / 33.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((145.0 / 22.0) * Y6_p0 + (5.0 / 11.0) * sqrt462 * Y6_p6 + (405.0 / 22.0) * Y4_p0 * R2 - 2.5 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (-33.75 * Y4_p0 + 7.5 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (-6.5625 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[23] = d[13] = -Y5_n3 * Y5_n4 * V[0] + ((7.0 / 11.0) * sqrt2 * Y8_p1 + (-7.0 / 286.0) * sqrt1430 * Y8_p7 + (-1.0 / 3.0) * sqrt42 * Y6_p1 * R2 + (35.0 / 33.0) * sqrt6 * Y2_p1 * R6) * V[1] + (1.25 * sqrt42 * Y6_p1 - 8.75 * sqrt6 * Y2_p1 * R4) * V[2] + 26.25 * sqrt6 * Y2_p1 * R2 * V[3] - 22.96875 * sqrt6 * Y2_p1 * V[4];
    d[22] = d[2] = -Y5_n3 * Y5_n5 * V[0] + ((35.0 / 286.0) * sqrt14 * Y8_p2 + (35.0 / 286.0) * sqrt143 * Y8_p8 + (-10.0 / 33.0) * sqrt42 * Y6_p2 * R2 + (405.0 / 143.0) * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt15 * Y2_p2 * R6) * V[1] + ((25.0 / 22.0) * sqrt42 * Y6_p2 + (-405.0 / 22.0) * Y4_p2 * R2 + 2.5 * sqrt15 * Y2_p2 * R4) * V[2] + (33.75 * Y4_p2 - 7.5 * sqrt15 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt15 * Y2_p2 * V[4];
    d[12] = -Y5_n4 * Y5_n4 * V[0] + ((-217.0 / 286.0) * Y8_p0 + (-21.0 / 286.0) * sqrt715 * Y8_p8 + (32.0 / 11.0) * Y6_p0 * R2 + (-405.0 / 143.0) * Y4_p0 * R4 + (-20.0 / 11.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((-120.0 / 11.0) * Y6_p0 + (405.0 / 22.0) * Y4_p0 * R2 + 15.0 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (-33.75 * Y4_p0 - 45.0 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (39.375 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];
    d[11] = d[1] = -Y5_n4 * Y5_n5 * V[0] + ((-21.0 / 286.0) * sqrt10 * Y8_p1 + (1.0 / 11.0) * sqrt210 * Y6_p1 * R2 + (-405.0 / 143.0) * Y4_p1 * R4 + (5.0 / 11.0) * sqrt30 * Y2_p1 * R6) * V[1] + ((-15.0 / 44.0) * sqrt210 * Y6_p1 + (405.0 / 22.0) * Y4_p1 * R2 - 3.75 * sqrt30 * Y2_p1 * R4) * V[2] + (-33.75 * Y4_p1 + 11.25 * sqrt30 * Y2_p1 * R2) * V[3] - 9.84375 * sqrt30 * Y2_p1 * V[4];
    d[0] = -Y5_n5 * Y5_n5 * V[0] + ((35.0 / 286.0) * Y8_p0 + (-10.0 / 11.0) * Y6_p0 * R2 + (405.0 / 143.0) * Y4_p0 * R4 + (-50.0 / 11.0) * Y2_p0 * R6 + 2.5 * R8) * V[1] + ((75.0 / 22.0) * Y6_p0 + (-405.0 / 22.0) * Y4_p0 * R2 + 37.5 * Y2_p0 * R4 - 22.5 * R6) * V[2] + (33.75 * Y4_p0 - 112.5 * Y2_p0 * R2 + 78.75 * R4) * V[3] + (98.4375 * Y2_p0 - 98.4375 * R2) * V[4] + 29.53125 * V[5];

    return out;
}

}  // namespace ovlab