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

#include "OverlapABRecHG.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_h_g(
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
    const auto sqrt22 = std::sqrt(22.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt858 = std::sqrt(858.0);
    const auto sqrt2002 = std::sqrt(2002.0);
    const auto sqrt4290 = std::sqrt(4290.0);
    const auto sqrt6006 = std::sqrt(6006.0);
    const auto sqrt30030 = std::sqrt(30030.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^4 · β^5 · p^{-9} · (s|s)
    //   V[1] ↔ α^3 · β^4 · p^{-8} · (s|s)
    //   V[2] ↔ α^2 · β^3 · p^{-7} · (s|s)
    //   V[3] ↔ α · β^2 · p^{-6} · (s|s)
    //   V[4] ↔ β · p^{-5} · (s|s)
    const auto exps_a  = bra.get_exponents();
    const auto coefs_a = bra.get_normalization_factors();
    const auto exps_b  = ket.get_exponents();
    const auto coefs_b = ket.get_normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 5> V = {0.0, 0.0, 0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < exps_a.size(); ++i)
    {
        const auto alpha = exps_a[i];
        const auto ca    = coefs_a[i];
        const auto alpha2 = alpha * alpha;
        const auto alpha3 = alpha2 * alpha;
        const auto alpha4 = alpha3 * alpha;
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

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha4 * beta5 * pinv9;
            V[1] += cab_ss * alpha3 * beta4 * pinv8;
            V[2] += cab_ss * alpha2 * beta3 * pinv7;
            V[3] += cab_ss * alpha * beta2 * pinv6;
            V[4] += cab_ss * beta * pinv5;
        }
    }

    // ---- Phase 3: fused M·V → 11 × 9 spherical block ----
    const auto Y1_n1 = harm::Y_ll_1_m_n1(AB_x, AB_y, AB_z);
    const auto Y1_p0 = harm::Y_ll_1_m_p0(AB_x, AB_y, AB_z);
    const auto Y1_p1 = harm::Y_ll_1_m_p1(AB_x, AB_y, AB_z);
    const auto Y3_n3 = harm::Y_ll_3_m_n3(AB_x, AB_y, AB_z);
    const auto Y3_n2 = harm::Y_ll_3_m_n2(AB_x, AB_y, AB_z);
    const auto Y3_n1 = harm::Y_ll_3_m_n1(AB_x, AB_y, AB_z);
    const auto Y3_p0 = harm::Y_ll_3_m_p0(AB_x, AB_y, AB_z);
    const auto Y3_p1 = harm::Y_ll_3_m_p1(AB_x, AB_y, AB_z);
    const auto Y3_p2 = harm::Y_ll_3_m_p2(AB_x, AB_y, AB_z);
    const auto Y3_p3 = harm::Y_ll_3_m_p3(AB_x, AB_y, AB_z);
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
    const auto Y7_n7 = harm::Y_ll_7_m_n7(AB_x, AB_y, AB_z);
    const auto Y7_n6 = harm::Y_ll_7_m_n6(AB_x, AB_y, AB_z);
    const auto Y7_n5 = harm::Y_ll_7_m_n5(AB_x, AB_y, AB_z);
    const auto Y7_n4 = harm::Y_ll_7_m_n4(AB_x, AB_y, AB_z);
    const auto Y7_n3 = harm::Y_ll_7_m_n3(AB_x, AB_y, AB_z);
    const auto Y7_n2 = harm::Y_ll_7_m_n2(AB_x, AB_y, AB_z);
    const auto Y7_n1 = harm::Y_ll_7_m_n1(AB_x, AB_y, AB_z);
    const auto Y7_p0 = harm::Y_ll_7_m_p0(AB_x, AB_y, AB_z);
    const auto Y7_p1 = harm::Y_ll_7_m_p1(AB_x, AB_y, AB_z);
    const auto Y7_p2 = harm::Y_ll_7_m_p2(AB_x, AB_y, AB_z);
    const auto Y7_p3 = harm::Y_ll_7_m_p3(AB_x, AB_y, AB_z);
    const auto Y7_p4 = harm::Y_ll_7_m_p4(AB_x, AB_y, AB_z);
    const auto Y7_p5 = harm::Y_ll_7_m_p5(AB_x, AB_y, AB_z);
    const auto Y7_p6 = harm::Y_ll_7_m_p6(AB_x, AB_y, AB_z);
    const auto Y7_p7 = harm::Y_ll_7_m_p7(AB_x, AB_y, AB_z);
    const auto R4 = R2 * R2;
    const auto R6 = R4 * R2;
    newints::Block out{11, 9, std::vector<double>(99, 0.0)};
    auto *d = out.data.data();
    d[49] = -Y5_p0 * Y4_p0 * V[0] + ((700.0 / 429.0) * Y7_p0 + (30.0 / 13.0) * Y5_p0 * R2 + (30.0 / 11.0) * Y3_p0 * R4 + (10.0 / 3.0) * Y1_p0 * R6) * V[1] + (-7.5 * Y5_p0 - 15.0 * Y3_p0 * R2 - 22.5 * Y1_p0 * R4) * V[2] + (22.5 * Y3_p0 + 52.5 * Y1_p0 * R2) * V[3] - 32.8125 * Y1_p0 * V[4];
    d[50] = -Y5_p0 * Y4_p1 * V[0] + ((115.0 / 858.0) * sqrt70 * Y7_p1 + (5.0 / 13.0) * sqrt6 * Y5_p1 * R2 + (1.0 / 22.0) * sqrt15 * Y3_p1 * R4 + (-2.0 / 3.0) * sqrt10 * Y1_p1 * R6) * V[1] + (-1.25 * sqrt6 * Y5_p1 - 0.25 * sqrt15 * Y3_p1 * R2 + 4.5 * sqrt10 * Y1_p1 * R4) * V[2] + (0.375 * sqrt15 * Y3_p1 - 10.5 * sqrt10 * Y1_p1 * R2) * V[3] + 6.5625 * sqrt10 * Y1_p1 * V[4];
    d[51] = -Y5_p0 * Y4_p2 * V[0] + ((-5.0 / 429.0) * sqrt210 * Y7_p2 + (-5.0 / 13.0) * sqrt21 * Y5_p2 * R2 + (-20.0 / 11.0) * sqrt3 * Y3_p2 * R4) * V[1] + (1.25 * sqrt21 * Y5_p2 + 10.0 * sqrt3 * Y3_p2 * R2) * V[2] - 15.0 * sqrt3 * Y3_p2 * V[3];
    d[52] = -Y5_p0 * Y4_p3 * V[0] + ((-245.0 / 858.0) * sqrt30 * Y7_p3 + (-30.0 / 13.0) * Y5_p3 * R2 + (15.0 / 22.0) * sqrt7 * Y3_p3 * R4) * V[1] + (7.5 * Y5_p3 - 3.75 * sqrt7 * Y3_p3 * R2) * V[2] + 5.625 * sqrt7 * Y3_p3 * V[3];
    d[53] = -Y5_p0 * Y4_p4 * V[0] + ((-70.0 / 429.0) * sqrt165 * Y7_p4 + (30.0 / 13.0) * Y5_p4 * R2) * V[1] - 7.5 * Y5_p4 * V[2];
    d[48] = -Y5_p0 * Y4_n1 * V[0] + ((115.0 / 858.0) * sqrt70 * Y7_n1 + (5.0 / 13.0) * sqrt6 * Y5_n1 * R2 + (1.0 / 22.0) * sqrt15 * Y3_n1 * R4 + (-2.0 / 3.0) * sqrt10 * Y1_n1 * R6) * V[1] + (-1.25 * sqrt6 * Y5_n1 - 0.25 * sqrt15 * Y3_n1 * R2 + 4.5 * sqrt10 * Y1_n1 * R4) * V[2] + (0.375 * sqrt15 * Y3_n1 - 10.5 * sqrt10 * Y1_n1 * R2) * V[3] + 6.5625 * sqrt10 * Y1_n1 * V[4];
    d[47] = -Y5_p0 * Y4_n2 * V[0] + ((-5.0 / 429.0) * sqrt210 * Y7_n2 + (-5.0 / 13.0) * sqrt21 * Y5_n2 * R2 + (-20.0 / 11.0) * sqrt3 * Y3_n2 * R4) * V[1] + (1.25 * sqrt21 * Y5_n2 + 10.0 * sqrt3 * Y3_n2 * R2) * V[2] - 15.0 * sqrt3 * Y3_n2 * V[3];
    d[46] = -Y5_p0 * Y4_n3 * V[0] + ((-245.0 / 858.0) * sqrt30 * Y7_n3 + (-30.0 / 13.0) * Y5_n3 * R2 + (15.0 / 22.0) * sqrt7 * Y3_n3 * R4) * V[1] + (7.5 * Y5_n3 - 3.75 * sqrt7 * Y3_n3 * R2) * V[2] + 5.625 * sqrt7 * Y3_n3 * V[3];
    d[45] = -Y5_p0 * Y4_n4 * V[0] + ((-70.0 / 429.0) * sqrt165 * Y7_n4 + (30.0 / 13.0) * Y5_n4 * R2) * V[1] - 7.5 * Y5_n4 * V[2];
    d[58] = -Y5_p1 * Y4_p0 * V[0] + ((5.0 / 39.0) * sqrt105 * Y7_p1 + (20.0 / 13.0) * Y5_p1 * R2 + 0.5 * sqrt10 * Y3_p1 * R4 + (2.0 / 3.0) * sqrt15 * Y1_p1 * R6) * V[1] + (-5.0 * Y5_p1 - 2.75 * sqrt10 * Y3_p1 * R2 - 4.5 * sqrt15 * Y1_p1 * R4) * V[2] + (4.125 * sqrt10 * Y3_p1 + 10.5 * sqrt15 * Y1_p1 * R2) * V[3] - 6.5625 * sqrt15 * Y1_p1 * V[4];
    d[59] = -Y5_p1 * Y4_p1 * V[0] + ((-35.0 / 429.0) * sqrt6 * Y7_p0 + (125.0 / 286.0) * sqrt7 * Y7_p2 + (5.0 / 13.0) * sqrt6 * Y5_p0 * R2 + (5.0 / 26.0) * sqrt70 * Y5_p2 * R2 + (19.0 / 22.0) * sqrt6 * Y3_p0 * R4 + (25.0 / 44.0) * sqrt10 * Y3_p2 * R4 + (4.0 / 3.0) * sqrt6 * Y1_p0 * R6) * V[1] + (-1.25 * sqrt6 * Y5_p0 - 0.625 * sqrt70 * Y5_p2 - 4.75 * sqrt6 * Y3_p0 * R2 - 3.125 * sqrt10 * Y3_p2 * R2 - 9.0 * sqrt6 * Y1_p0 * R4) * V[2] + (7.125 * sqrt6 * Y3_p0 + 4.6875 * sqrt10 * Y3_p2 + 21.0 * sqrt6 * Y1_p0 * R2) * V[3] - 13.125 * sqrt6 * Y1_p0 * V[4];
    d[60] = -Y5_p1 * Y4_p2 * V[0] + ((145.0 / 858.0) * sqrt21 * Y7_p1 + (5.0 / 22.0) * sqrt7 * Y7_p3 + (10.0 / 13.0) * sqrt5 * Y5_p1 * R2 + (43.0 / 44.0) * sqrt2 * Y3_p1 * R4 + (-15.0 / 44.0) * sqrt30 * Y3_p3 * R4 + (-2.0 / 3.0) * sqrt3 * Y1_p1 * R6) * V[1] + (-2.5 * sqrt5 * Y5_p1 - 5.375 * sqrt2 * Y3_p1 * R2 + 1.875 * sqrt30 * Y3_p3 * R2 + 4.5 * sqrt3 * Y1_p1 * R4) * V[2] + (8.0625 * sqrt2 * Y3_p1 - 2.8125 * sqrt30 * Y3_p3 - 10.5 * sqrt3 * Y1_p1 * R2) * V[3] + 6.5625 * sqrt3 * Y1_p1 * V[4];
    d[61] = -Y5_p1 * Y4_p3 * V[0] + ((35.0 / 26.0) * Y7_p2 + (-35.0 / 286.0) * sqrt22 * Y7_p4 + (5.0 / 26.0) * sqrt10 * Y5_p2 * R2 + (-5.0 / 13.0) * sqrt30 * Y5_p4 * R2 - 0.25 * sqrt70 * Y3_p2 * R4) * V[1] + (-0.625 * sqrt10 * Y5_p2 + 1.25 * sqrt30 * Y5_p4 + 1.375 * sqrt70 * Y3_p2 * R2) * V[2] - 2.0625 * sqrt70 * Y3_p2 * V[3];
    d[62] = -Y5_p1 * Y4_p4 * V[0] + ((175.0 / 143.0) * Y7_p3 + (-70.0 / 143.0) * sqrt11 * Y7_p5 + (-5.0 / 13.0) * sqrt30 * Y5_p3 * R2 + (5.0 / 13.0) * sqrt6 * Y5_p5 * R2 + (1.0 / 22.0) * sqrt210 * Y3_p3 * R4) * V[1] + (1.25 * sqrt30 * Y5_p3 - 1.25 * sqrt6 * Y5_p5 - 0.25 * sqrt210 * Y3_p3 * R2) * V[2] + 0.375 * sqrt210 * Y3_p3 * V[3];
    d[57] = -Y5_p1 * Y4_n1 * V[0] + ((125.0 / 286.0) * sqrt7 * Y7_n2 + (5.0 / 26.0) * sqrt70 * Y5_n2 * R2 + (25.0 / 44.0) * sqrt10 * Y3_n2 * R4) * V[1] + (-0.625 * sqrt70 * Y5_n2 - 3.125 * sqrt10 * Y3_n2 * R2) * V[2] + 4.6875 * sqrt10 * Y3_n2 * V[3];
    d[56] = -Y5_p1 * Y4_n2 * V[0] + ((5.0 / 22.0) * sqrt7 * Y7_n3 + (145.0 / 858.0) * sqrt21 * Y7_n1 + (10.0 / 13.0) * sqrt5 * Y5_n1 * R2 + (-15.0 / 44.0) * sqrt30 * Y3_n3 * R4 + (43.0 / 44.0) * sqrt2 * Y3_n1 * R4 + (-2.0 / 3.0) * sqrt3 * Y1_n1 * R6) * V[1] + (-2.5 * sqrt5 * Y5_n1 + 1.875 * sqrt30 * Y3_n3 * R2 - 5.375 * sqrt2 * Y3_n1 * R2 + 4.5 * sqrt3 * Y1_n1 * R4) * V[2] + (-2.8125 * sqrt30 * Y3_n3 + 8.0625 * sqrt2 * Y3_n1 - 10.5 * sqrt3 * Y1_n1 * R2) * V[3] + 6.5625 * sqrt3 * Y1_n1 * V[4];
    d[55] = -Y5_p1 * Y4_n3 * V[0] + ((-35.0 / 286.0) * sqrt22 * Y7_n4 + (35.0 / 26.0) * Y7_n2 + (-5.0 / 13.0) * sqrt30 * Y5_n4 * R2 + (5.0 / 26.0) * sqrt10 * Y5_n2 * R2 - 0.25 * sqrt70 * Y3_n2 * R4) * V[1] + (1.25 * sqrt30 * Y5_n4 - 0.625 * sqrt10 * Y5_n2 + 1.375 * sqrt70 * Y3_n2 * R2) * V[2] - 2.0625 * sqrt70 * Y3_n2 * V[3];
    d[54] = -Y5_p1 * Y4_n4 * V[0] + ((-70.0 / 143.0) * sqrt11 * Y7_n5 + (175.0 / 143.0) * Y7_n3 + (5.0 / 13.0) * sqrt6 * Y5_n5 * R2 + (-5.0 / 13.0) * sqrt30 * Y5_n3 * R2 + (1.0 / 22.0) * sqrt210 * Y3_n3 * R4) * V[1] + (-1.25 * sqrt6 * Y5_n5 + 1.25 * sqrt30 * Y5_n3 - 0.25 * sqrt210 * Y3_n3 * R2) * V[2] + 0.375 * sqrt210 * Y3_n3 * V[3];
    d[67] = -Y5_p2 * Y4_p0 * V[0] + ((20.0 / 143.0) * sqrt10 * Y7_p2 + (-5.0 / 13.0) * Y5_p2 * R2 + (-5.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[1] + (1.25 * Y5_p2 + 2.5 * sqrt7 * Y3_p2 * R2) * V[2] - 3.75 * sqrt7 * Y3_p2 * V[3];
    d[68] = -Y5_p2 * Y4_p1 * V[0] + ((205.0 / 858.0) * sqrt6 * Y7_p1 + (215.0 / 286.0) * sqrt2 * Y7_p3 + (5.0 / 26.0) * sqrt70 * Y5_p1 * R2 + (5.0 / 13.0) * sqrt15 * Y5_p3 * R2 + (8.0 / 11.0) * sqrt7 * Y3_p1 * R4 + (5.0 / 22.0) * sqrt105 * Y3_p3 * R4 + (1.0 / 3.0) * sqrt42 * Y1_p1 * R6) * V[1] + (-0.625 * sqrt70 * Y5_p1 - 1.25 * sqrt15 * Y5_p3 - 4.0 * sqrt7 * Y3_p1 * R2 - 1.25 * sqrt105 * Y3_p3 * R2 - 2.25 * sqrt42 * Y1_p1 * R4) * V[2] + (6.0 * sqrt7 * Y3_p1 + 1.875 * sqrt105 * Y3_p3 + 5.25 * sqrt42 * Y1_p1 * R2) * V[3] - 3.28125 * sqrt42 * Y1_p1 * V[4];
    d[69] = -Y5_p2 * Y4_p2 * V[0] + ((-160.0 / 429.0) * sqrt21 * Y7_p0 + (50.0 / 143.0) * sqrt11 * Y7_p4 + (-5.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (5.0 / 13.0) * sqrt15 * Y5_p4 * R2 + (1.0 / 11.0) * sqrt21 * Y3_p0 * R4 + (2.0 / 3.0) * sqrt21 * Y1_p0 * R6) * V[1] + (1.25 * sqrt21 * Y5_p0 - 1.25 * sqrt15 * Y5_p4 - 0.5 * sqrt21 * Y3_p0 * R2 - 4.5 * sqrt21 * Y1_p0 * R4) * V[2] + (0.75 * sqrt21 * Y3_p0 + 10.5 * sqrt21 * Y1_p0 * R2) * V[3] - 6.5625 * sqrt21 * Y1_p0 * V[4];
    d[70] = -Y5_p2 * Y4_p3 * V[0] + ((-175.0 / 858.0) * sqrt42 * Y7_p1 + (5.0 / 286.0) * sqrt154 * Y7_p5 + (5.0 / 26.0) * sqrt10 * Y5_p1 * R2 + (-5.0 / 13.0) * sqrt21 * Y5_p5 * R2 + (49.0 / 22.0) * Y3_p1 * R4 + (-1.0 / 3.0) * sqrt6 * Y1_p1 * R6) * V[1] + (-0.625 * sqrt10 * Y5_p1 + 1.25 * sqrt21 * Y5_p5 - 12.25 * Y3_p1 * R2 + 2.25 * sqrt6 * Y1_p1 * R4) * V[2] + (18.375 * Y3_p1 - 5.25 * sqrt6 * Y1_p1 * R2) * V[3] + 3.28125 * sqrt6 * Y1_p1 * V[4];
    d[71] = -Y5_p2 * Y4_p4 * V[0] + ((-35.0 / 143.0) * sqrt14 * Y7_p2 + (-5.0 / 143.0) * sqrt2002 * Y7_p6 + (5.0 / 13.0) * sqrt35 * Y5_p2 * R2 + (-7.0 / 11.0) * sqrt5 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt35 * Y5_p2 + 3.5 * sqrt5 * Y3_p2 * R2) * V[2] - 5.25 * sqrt5 * Y3_p2 * V[3];
    d[66] = -Y5_p2 * Y4_n1 * V[0] + ((215.0 / 286.0) * sqrt2 * Y7_n3 + (-205.0 / 858.0) * sqrt6 * Y7_n1 + (5.0 / 13.0) * sqrt15 * Y5_n3 * R2 + (-5.0 / 26.0) * sqrt70 * Y5_n1 * R2 + (5.0 / 22.0) * sqrt105 * Y3_n3 * R4 + (-8.0 / 11.0) * sqrt7 * Y3_n1 * R4 + (-1.0 / 3.0) * sqrt42 * Y1_n1 * R6) * V[1] + (-1.25 * sqrt15 * Y5_n3 + 0.625 * sqrt70 * Y5_n1 - 1.25 * sqrt105 * Y3_n3 * R2 + 4.0 * sqrt7 * Y3_n1 * R2 + 2.25 * sqrt42 * Y1_n1 * R4) * V[2] + (1.875 * sqrt105 * Y3_n3 - 6.0 * sqrt7 * Y3_n1 - 5.25 * sqrt42 * Y1_n1 * R2) * V[3] + 3.28125 * sqrt42 * Y1_n1 * V[4];
    d[65] = -Y5_p2 * Y4_n2 * V[0] + ((50.0 / 143.0) * sqrt11 * Y7_n4 + (5.0 / 13.0) * sqrt15 * Y5_n4 * R2) * V[1] - 1.25 * sqrt15 * Y5_n4 * V[2];
    d[64] = -Y5_p2 * Y4_n3 * V[0] + ((5.0 / 286.0) * sqrt154 * Y7_n5 + (-175.0 / 858.0) * sqrt42 * Y7_n1 + (-5.0 / 13.0) * sqrt21 * Y5_n5 * R2 + (5.0 / 26.0) * sqrt10 * Y5_n1 * R2 + (49.0 / 22.0) * Y3_n1 * R4 + (-1.0 / 3.0) * sqrt6 * Y1_n1 * R6) * V[1] + (1.25 * sqrt21 * Y5_n5 - 0.625 * sqrt10 * Y5_n1 - 12.25 * Y3_n1 * R2 + 2.25 * sqrt6 * Y1_n1 * R4) * V[2] + (18.375 * Y3_n1 - 5.25 * sqrt6 * Y1_n1 * R2) * V[3] + 3.28125 * sqrt6 * Y1_n1 * V[4];
    d[63] = -Y5_p2 * Y4_n4 * V[0] + ((-5.0 / 143.0) * sqrt2002 * Y7_n6 + (-35.0 / 143.0) * sqrt14 * Y7_n2 + (5.0 / 13.0) * sqrt35 * Y5_n2 * R2 + (-7.0 / 11.0) * sqrt5 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt35 * Y5_n2 + 3.5 * sqrt5 * Y3_n2 * R2) * V[2] - 5.25 * sqrt5 * Y3_n2 * V[3];
    d[76] = -Y5_p3 * Y4_p0 * V[0] + ((-115.0 / 858.0) * sqrt30 * Y7_p3 + (-30.0 / 13.0) * Y5_p3 * R2 + (-15.0 / 11.0) * sqrt7 * Y3_p3 * R4) * V[1] + (7.5 * Y5_p3 + 7.5 * sqrt7 * Y3_p3 * R2) * V[2] - 11.25 * sqrt7 * Y3_p3 * V[3];
    d[77] = -Y5_p3 * Y4_p1 * V[0] + ((815.0 / 1716.0) * sqrt6 * Y7_p2 + (5.0 / 66.0) * sqrt33 * Y7_p4 + (5.0 / 13.0) * sqrt15 * Y5_p2 * R2 + (1.0 / 22.0) * sqrt105 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt15 * Y5_p2 - 0.25 * sqrt105 * Y3_p2 * R2) * V[2] + 0.375 * sqrt105 * Y3_p2 * V[3];
    d[78] = -Y5_p3 * Y4_p2 * V[0] + ((-125.0 / 132.0) * sqrt2 * Y7_p1 + (265.0 / 1716.0) * sqrt66 * Y7_p5 + (30.0 / 13.0) * Y5_p5 * R2 + (9.0 / 22.0) * sqrt21 * Y3_p1 * R4 + (2.0 / 3.0) * sqrt14 * Y1_p1 * R6) * V[1] + (-7.5 * Y5_p5 - 2.25 * sqrt21 * Y3_p1 * R2 - 4.5 * sqrt14 * Y1_p1 * R4) * V[2] + (3.375 * sqrt21 * Y3_p1 + 10.5 * sqrt14 * Y1_p1 * R2) * V[3] - 6.5625 * sqrt14 * Y1_p1 * V[4];
    d[79] = -Y5_p3 * Y4_p3 * V[0] + ((665.0 / 429.0) * Y7_p0 + (25.0 / 1716.0) * sqrt6006 * Y7_p6 + (-30.0 / 13.0) * Y5_p0 * R2 + (-21.0 / 11.0) * Y3_p0 * R4 + (8.0 / 3.0) * Y1_p0 * R6) * V[1] + (7.5 * Y5_p0 + 10.5 * Y3_p0 * R2 - 18.0 * Y1_p0 * R4) * V[2] + (-15.75 * Y3_p0 + 42.0 * Y1_p0 * R2) * V[3] - 26.25 * Y1_p0 * V[4];
    d[80] = -Y5_p3 * Y4_p4 * V[0] + ((70.0 / 429.0) * sqrt14 * Y7_p1 + (-35.0 / 858.0) * sqrt858 * Y7_p7 + (-5.0 / 13.0) * sqrt30 * Y5_p1 * R2 + (14.0 / 11.0) * sqrt3 * Y3_p1 * R4 + (-1.0 / 3.0) * sqrt2 * Y1_p1 * R6) * V[1] + (1.25 * sqrt30 * Y5_p1 - 7.0 * sqrt3 * Y3_p1 * R2 + 2.25 * sqrt2 * Y1_p1 * R4) * V[2] + (10.5 * sqrt3 * Y3_p1 - 5.25 * sqrt2 * Y1_p1 * R2) * V[3] + 3.28125 * sqrt2 * Y1_p1 * V[4];
    d[75] = -Y5_p3 * Y4_n1 * V[0] + ((5.0 / 66.0) * sqrt33 * Y7_n4 + (-815.0 / 1716.0) * sqrt6 * Y7_n2 + (-5.0 / 13.0) * sqrt15 * Y5_n2 * R2 + (-1.0 / 22.0) * sqrt105 * Y3_n2 * R4) * V[1] + (1.25 * sqrt15 * Y5_n2 + 0.25 * sqrt105 * Y3_n2 * R2) * V[2] - 0.375 * sqrt105 * Y3_n2 * V[3];
    d[74] = -Y5_p3 * Y4_n2 * V[0] + ((265.0 / 1716.0) * sqrt66 * Y7_n5 + (125.0 / 132.0) * sqrt2 * Y7_n1 + (30.0 / 13.0) * Y5_n5 * R2 + (-9.0 / 22.0) * sqrt21 * Y3_n1 * R4 + (-2.0 / 3.0) * sqrt14 * Y1_n1 * R6) * V[1] + (-7.5 * Y5_n5 + 2.25 * sqrt21 * Y3_n1 * R2 + 4.5 * sqrt14 * Y1_n1 * R4) * V[2] + (-3.375 * sqrt21 * Y3_n1 - 10.5 * sqrt14 * Y1_n1 * R2) * V[3] + 6.5625 * sqrt14 * Y1_n1 * V[4];
    d[73] = -Y5_p3 * Y4_n3 * V[0] + (25.0 / 1716.0) * sqrt6006 * Y7_n6 * V[1];
    d[72] = -Y5_p3 * Y4_n4 * V[0] + ((-35.0 / 858.0) * sqrt858 * Y7_n7 + (70.0 / 429.0) * sqrt14 * Y7_n1 + (-5.0 / 13.0) * sqrt30 * Y5_n1 * R2 + (14.0 / 11.0) * sqrt3 * Y3_n1 * R4 + (-1.0 / 3.0) * sqrt2 * Y1_n1 * R6) * V[1] + (1.25 * sqrt30 * Y5_n1 - 7.0 * sqrt3 * Y3_n1 * R2 + 2.25 * sqrt2 * Y1_n1 * R4) * V[2] + (10.5 * sqrt3 * Y3_n1 - 5.25 * sqrt2 * Y1_n1 * R2) * V[3] + 3.28125 * sqrt2 * Y1_n1 * V[4];
    d[85] = -Y5_p4 * Y4_p0 * V[0] + ((-20.0 / 143.0) * sqrt165 * Y7_p4 + (-30.0 / 13.0) * Y5_p4 * R2) * V[1] + 7.5 * Y5_p4 * V[2];
    d[86] = -Y5_p4 * Y4_p1 * V[0] + ((25.0 / 44.0) * sqrt6 * Y7_p3 + (-45.0 / 572.0) * sqrt66 * Y7_p5 + (-30.0 / 13.0) * Y5_p5 * R2 + (-9.0 / 22.0) * sqrt35 * Y3_p3 * R4) * V[1] + (7.5 * Y5_p5 + 2.25 * sqrt35 * Y3_p3 * R2) * V[2] - 3.375 * sqrt35 * Y3_p3 * V[3];
    d[87] = -Y5_p4 * Y4_p2 * V[0] + ((-135.0 / 286.0) * sqrt6 * Y7_p2 + (5.0 / 286.0) * sqrt858 * Y7_p6 + (5.0 / 13.0) * sqrt15 * Y5_p2 * R2 + (2.0 / 11.0) * sqrt105 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt15 * Y5_p2 - sqrt105 * Y3_p2 * R2) * V[2] + 1.5 * sqrt105 * Y3_p2 * V[3];
    d[88] = -Y5_p4 * Y4_p3 * V[0] + ((115.0 / 572.0) * sqrt14 * Y7_p1 + (35.0 / 572.0) * sqrt858 * Y7_p7 + (-5.0 / 13.0) * sqrt30 * Y5_p1 * R2 + (7.0 / 22.0) * sqrt3 * Y3_p1 * R4 + 2.0 * sqrt2 * Y1_p1 * R6) * V[1] + (1.25 * sqrt30 * Y5_p1 - 1.75 * sqrt3 * Y3_p1 * R2 - 13.5 * sqrt2 * Y1_p1 * R4) * V[2] + (2.625 * sqrt3 * Y3_p1 + 31.5 * sqrt2 * Y1_p1 * R2) * V[3] - 19.6875 * sqrt2 * Y1_p1 * V[4];
    d[89] = -Y5_p4 * Y4_p4 * V[0] + ((-70.0 / 143.0) * Y7_p0 + (30.0 / 13.0) * Y5_p0 * R2 + (-42.0 / 11.0) * Y3_p0 * R4 + 2.0 * Y1_p0 * R6) * V[1] + (-7.5 * Y5_p0 + 21.0 * Y3_p0 * R2 - 13.5 * Y1_p0 * R4) * V[2] + (-31.5 * Y3_p0 + 31.5 * Y1_p0 * R2) * V[3] - 19.6875 * Y1_p0 * V[4];
    d[84] = -Y5_p4 * Y4_n1 * V[0] + ((-45.0 / 572.0) * sqrt66 * Y7_n5 + (-25.0 / 44.0) * sqrt6 * Y7_n3 + (-30.0 / 13.0) * Y5_n5 * R2 + (9.0 / 22.0) * sqrt35 * Y3_n3 * R4) * V[1] + (7.5 * Y5_n5 - 2.25 * sqrt35 * Y3_n3 * R2) * V[2] + 3.375 * sqrt35 * Y3_n3 * V[3];
    d[83] = -Y5_p4 * Y4_n2 * V[0] + ((5.0 / 286.0) * sqrt858 * Y7_n6 + (135.0 / 286.0) * sqrt6 * Y7_n2 + (-5.0 / 13.0) * sqrt15 * Y5_n2 * R2 + (-2.0 / 11.0) * sqrt105 * Y3_n2 * R4) * V[1] + (1.25 * sqrt15 * Y5_n2 + sqrt105 * Y3_n2 * R2) * V[2] - 1.5 * sqrt105 * Y3_n2 * V[3];
    d[82] = -Y5_p4 * Y4_n3 * V[0] + ((35.0 / 572.0) * sqrt858 * Y7_n7 + (-115.0 / 572.0) * sqrt14 * Y7_n1 + (5.0 / 13.0) * sqrt30 * Y5_n1 * R2 + (-7.0 / 22.0) * sqrt3 * Y3_n1 * R4 - 2.0 * sqrt2 * Y1_n1 * R6) * V[1] + (-1.25 * sqrt30 * Y5_n1 + 1.75 * sqrt3 * Y3_n1 * R2 + 13.5 * sqrt2 * Y1_n1 * R4) * V[2] + (-2.625 * sqrt3 * Y3_n1 - 31.5 * sqrt2 * Y1_n1 * R2) * V[3] + 19.6875 * sqrt2 * Y1_n1 * V[4];
    d[81] = -Y5_p4 * Y4_n4 * V[0];
    d[94] = -Y5_p5 * Y4_p0 * V[0] + ((-75.0 / 286.0) * sqrt66 * Y7_p5 + (30.0 / 13.0) * Y5_p5 * R2) * V[1] - 7.5 * Y5_p5 * V[2];
    d[95] = -Y5_p5 * Y4_p1 * V[0] + ((25.0 / 286.0) * sqrt165 * Y7_p4 + (-15.0 / 572.0) * sqrt4290 * Y7_p6 + (-30.0 / 13.0) * Y5_p4 * R2) * V[1] + 7.5 * Y5_p4 * V[2];
    d[96] = -Y5_p5 * Y4_p2 * V[0] + ((-75.0 / 572.0) * sqrt30 * Y7_p3 + (-5.0 / 572.0) * sqrt30030 * Y7_p7 + (30.0 / 13.0) * Y5_p3 * R2 + (-15.0 / 22.0) * sqrt7 * Y3_p3 * R4) * V[1] + (-7.5 * Y5_p3 + 3.75 * sqrt7 * Y3_p3 * R2) * V[2] - 5.625 * sqrt7 * Y3_p3 * V[3];
    d[97] = -Y5_p5 * Y4_p3 * V[0] + ((15.0 / 572.0) * sqrt210 * Y7_p2 + (-5.0 / 13.0) * sqrt21 * Y5_p2 * R2 + (35.0 / 22.0) * sqrt3 * Y3_p2 * R4) * V[1] + (1.25 * sqrt21 * Y5_p2 - 8.75 * sqrt3 * Y3_p2 * R2) * V[2] + 13.125 * sqrt3 * Y3_p2 * V[3];
    d[98] = -Y5_p5 * Y4_p4 * V[0] + ((-5.0 / 286.0) * sqrt70 * Y7_p1 + (5.0 / 13.0) * sqrt6 * Y5_p1 * R2 + (-7.0 / 11.0) * sqrt15 * Y3_p1 * R4 + sqrt10 * Y1_p1 * R6) * V[1] + (-1.25 * sqrt6 * Y5_p1 + 3.5 * sqrt15 * Y3_p1 * R2 - 6.75 * sqrt10 * Y1_p1 * R4) * V[2] + (-5.25 * sqrt15 * Y3_p1 + 15.75 * sqrt10 * Y1_p1 * R2) * V[3] - 9.84375 * sqrt10 * Y1_p1 * V[4];
    d[93] = -Y5_p5 * Y4_n1 * V[0] + ((-15.0 / 572.0) * sqrt4290 * Y7_n6 + (-25.0 / 286.0) * sqrt165 * Y7_n4 + (30.0 / 13.0) * Y5_n4 * R2) * V[1] - 7.5 * Y5_n4 * V[2];
    d[92] = -Y5_p5 * Y4_n2 * V[0] + ((-5.0 / 572.0) * sqrt30030 * Y7_n7 + (75.0 / 572.0) * sqrt30 * Y7_n3 + (-30.0 / 13.0) * Y5_n3 * R2 + (15.0 / 22.0) * sqrt7 * Y3_n3 * R4) * V[1] + (7.5 * Y5_n3 - 3.75 * sqrt7 * Y3_n3 * R2) * V[2] + 5.625 * sqrt7 * Y3_n3 * V[3];
    d[91] = -Y5_p5 * Y4_n3 * V[0] + ((-15.0 / 572.0) * sqrt210 * Y7_n2 + (5.0 / 13.0) * sqrt21 * Y5_n2 * R2 + (-35.0 / 22.0) * sqrt3 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt21 * Y5_n2 + 8.75 * sqrt3 * Y3_n2 * R2) * V[2] - 13.125 * sqrt3 * Y3_n2 * V[3];
    d[90] = -Y5_p5 * Y4_n4 * V[0] + ((5.0 / 286.0) * sqrt70 * Y7_n1 + (-5.0 / 13.0) * sqrt6 * Y5_n1 * R2 + (7.0 / 11.0) * sqrt15 * Y3_n1 * R4 - sqrt10 * Y1_n1 * R6) * V[1] + (1.25 * sqrt6 * Y5_n1 - 3.5 * sqrt15 * Y3_n1 * R2 + 6.75 * sqrt10 * Y1_n1 * R4) * V[2] + (5.25 * sqrt15 * Y3_n1 - 15.75 * sqrt10 * Y1_n1 * R2) * V[3] + 9.84375 * sqrt10 * Y1_n1 * V[4];
    d[40] = -Y5_n1 * Y4_p0 * V[0] + ((5.0 / 39.0) * sqrt105 * Y7_n1 + (20.0 / 13.0) * Y5_n1 * R2 + 0.5 * sqrt10 * Y3_n1 * R4 + (2.0 / 3.0) * sqrt15 * Y1_n1 * R6) * V[1] + (-5.0 * Y5_n1 - 2.75 * sqrt10 * Y3_n1 * R2 - 4.5 * sqrt15 * Y1_n1 * R4) * V[2] + (4.125 * sqrt10 * Y3_n1 + 10.5 * sqrt15 * Y1_n1 * R2) * V[3] - 6.5625 * sqrt15 * Y1_n1 * V[4];
    d[41] = -Y5_n1 * Y4_p1 * V[0] + ((125.0 / 286.0) * sqrt7 * Y7_n2 + (5.0 / 26.0) * sqrt70 * Y5_n2 * R2 + (25.0 / 44.0) * sqrt10 * Y3_n2 * R4) * V[1] + (-0.625 * sqrt70 * Y5_n2 - 3.125 * sqrt10 * Y3_n2 * R2) * V[2] + 4.6875 * sqrt10 * Y3_n2 * V[3];
    d[42] = -Y5_n1 * Y4_p2 * V[0] + ((5.0 / 22.0) * sqrt7 * Y7_n3 + (-145.0 / 858.0) * sqrt21 * Y7_n1 + (-10.0 / 13.0) * sqrt5 * Y5_n1 * R2 + (-15.0 / 44.0) * sqrt30 * Y3_n3 * R4 + (-43.0 / 44.0) * sqrt2 * Y3_n1 * R4 + (2.0 / 3.0) * sqrt3 * Y1_n1 * R6) * V[1] + (2.5 * sqrt5 * Y5_n1 + 1.875 * sqrt30 * Y3_n3 * R2 + 5.375 * sqrt2 * Y3_n1 * R2 - 4.5 * sqrt3 * Y1_n1 * R4) * V[2] + (-2.8125 * sqrt30 * Y3_n3 - 8.0625 * sqrt2 * Y3_n1 + 10.5 * sqrt3 * Y1_n1 * R2) * V[3] - 6.5625 * sqrt3 * Y1_n1 * V[4];
    d[43] = -Y5_n1 * Y4_p3 * V[0] + ((-35.0 / 286.0) * sqrt22 * Y7_n4 + (-35.0 / 26.0) * Y7_n2 + (-5.0 / 13.0) * sqrt30 * Y5_n4 * R2 + (-5.0 / 26.0) * sqrt10 * Y5_n2 * R2 + 0.25 * sqrt70 * Y3_n2 * R4) * V[1] + (1.25 * sqrt30 * Y5_n4 + 0.625 * sqrt10 * Y5_n2 - 1.375 * sqrt70 * Y3_n2 * R2) * V[2] + 2.0625 * sqrt70 * Y3_n2 * V[3];
    d[44] = -Y5_n1 * Y4_p4 * V[0] + ((-70.0 / 143.0) * sqrt11 * Y7_n5 + (-175.0 / 143.0) * Y7_n3 + (5.0 / 13.0) * sqrt6 * Y5_n5 * R2 + (5.0 / 13.0) * sqrt30 * Y5_n3 * R2 + (-1.0 / 22.0) * sqrt210 * Y3_n3 * R4) * V[1] + (-1.25 * sqrt6 * Y5_n5 - 1.25 * sqrt30 * Y5_n3 + 0.25 * sqrt210 * Y3_n3 * R2) * V[2] - 0.375 * sqrt210 * Y3_n3 * V[3];
    d[39] = -Y5_n1 * Y4_n1 * V[0] + ((-35.0 / 429.0) * sqrt6 * Y7_p0 + (-125.0 / 286.0) * sqrt7 * Y7_p2 + (5.0 / 13.0) * sqrt6 * Y5_p0 * R2 + (-5.0 / 26.0) * sqrt70 * Y5_p2 * R2 + (19.0 / 22.0) * sqrt6 * Y3_p0 * R4 + (-25.0 / 44.0) * sqrt10 * Y3_p2 * R4 + (4.0 / 3.0) * sqrt6 * Y1_p0 * R6) * V[1] + (-1.25 * sqrt6 * Y5_p0 + 0.625 * sqrt70 * Y5_p2 - 4.75 * sqrt6 * Y3_p0 * R2 + 3.125 * sqrt10 * Y3_p2 * R2 - 9.0 * sqrt6 * Y1_p0 * R4) * V[2] + (7.125 * sqrt6 * Y3_p0 - 4.6875 * sqrt10 * Y3_p2 + 21.0 * sqrt6 * Y1_p0 * R2) * V[3] - 13.125 * sqrt6 * Y1_p0 * V[4];
    d[38] = -Y5_n1 * Y4_n2 * V[0] + ((145.0 / 858.0) * sqrt21 * Y7_p1 + (-5.0 / 22.0) * sqrt7 * Y7_p3 + (10.0 / 13.0) * sqrt5 * Y5_p1 * R2 + (43.0 / 44.0) * sqrt2 * Y3_p1 * R4 + (15.0 / 44.0) * sqrt30 * Y3_p3 * R4 + (-2.0 / 3.0) * sqrt3 * Y1_p1 * R6) * V[1] + (-2.5 * sqrt5 * Y5_p1 - 5.375 * sqrt2 * Y3_p1 * R2 - 1.875 * sqrt30 * Y3_p3 * R2 + 4.5 * sqrt3 * Y1_p1 * R4) * V[2] + (8.0625 * sqrt2 * Y3_p1 + 2.8125 * sqrt30 * Y3_p3 - 10.5 * sqrt3 * Y1_p1 * R2) * V[3] + 6.5625 * sqrt3 * Y1_p1 * V[4];
    d[37] = -Y5_n1 * Y4_n3 * V[0] + ((35.0 / 26.0) * Y7_p2 + (35.0 / 286.0) * sqrt22 * Y7_p4 + (5.0 / 26.0) * sqrt10 * Y5_p2 * R2 + (5.0 / 13.0) * sqrt30 * Y5_p4 * R2 - 0.25 * sqrt70 * Y3_p2 * R4) * V[1] + (-0.625 * sqrt10 * Y5_p2 - 1.25 * sqrt30 * Y5_p4 + 1.375 * sqrt70 * Y3_p2 * R2) * V[2] - 2.0625 * sqrt70 * Y3_p2 * V[3];
    d[36] = -Y5_n1 * Y4_n4 * V[0] + ((175.0 / 143.0) * Y7_p3 + (70.0 / 143.0) * sqrt11 * Y7_p5 + (-5.0 / 13.0) * sqrt30 * Y5_p3 * R2 + (-5.0 / 13.0) * sqrt6 * Y5_p5 * R2 + (1.0 / 22.0) * sqrt210 * Y3_p3 * R4) * V[1] + (1.25 * sqrt30 * Y5_p3 + 1.25 * sqrt6 * Y5_p5 - 0.25 * sqrt210 * Y3_p3 * R2) * V[2] + 0.375 * sqrt210 * Y3_p3 * V[3];
    d[31] = -Y5_n2 * Y4_p0 * V[0] + ((20.0 / 143.0) * sqrt10 * Y7_n2 + (-5.0 / 13.0) * Y5_n2 * R2 + (-5.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[1] + (1.25 * Y5_n2 + 2.5 * sqrt7 * Y3_n2 * R2) * V[2] - 3.75 * sqrt7 * Y3_n2 * V[3];
    d[32] = -Y5_n2 * Y4_p1 * V[0] + ((215.0 / 286.0) * sqrt2 * Y7_n3 + (205.0 / 858.0) * sqrt6 * Y7_n1 + (5.0 / 13.0) * sqrt15 * Y5_n3 * R2 + (5.0 / 26.0) * sqrt70 * Y5_n1 * R2 + (5.0 / 22.0) * sqrt105 * Y3_n3 * R4 + (8.0 / 11.0) * sqrt7 * Y3_n1 * R4 + (1.0 / 3.0) * sqrt42 * Y1_n1 * R6) * V[1] + (-1.25 * sqrt15 * Y5_n3 - 0.625 * sqrt70 * Y5_n1 - 1.25 * sqrt105 * Y3_n3 * R2 - 4.0 * sqrt7 * Y3_n1 * R2 - 2.25 * sqrt42 * Y1_n1 * R4) * V[2] + (1.875 * sqrt105 * Y3_n3 + 6.0 * sqrt7 * Y3_n1 + 5.25 * sqrt42 * Y1_n1 * R2) * V[3] - 3.28125 * sqrt42 * Y1_n1 * V[4];
    d[33] = -Y5_n2 * Y4_p2 * V[0] + ((50.0 / 143.0) * sqrt11 * Y7_n4 + (5.0 / 13.0) * sqrt15 * Y5_n4 * R2) * V[1] - 1.25 * sqrt15 * Y5_n4 * V[2];
    d[34] = -Y5_n2 * Y4_p3 * V[0] + ((5.0 / 286.0) * sqrt154 * Y7_n5 + (175.0 / 858.0) * sqrt42 * Y7_n1 + (-5.0 / 13.0) * sqrt21 * Y5_n5 * R2 + (-5.0 / 26.0) * sqrt10 * Y5_n1 * R2 + (-49.0 / 22.0) * Y3_n1 * R4 + (1.0 / 3.0) * sqrt6 * Y1_n1 * R6) * V[1] + (1.25 * sqrt21 * Y5_n5 + 0.625 * sqrt10 * Y5_n1 + 12.25 * Y3_n1 * R2 - 2.25 * sqrt6 * Y1_n1 * R4) * V[2] + (-18.375 * Y3_n1 + 5.25 * sqrt6 * Y1_n1 * R2) * V[3] - 3.28125 * sqrt6 * Y1_n1 * V[4];
    d[35] = -Y5_n2 * Y4_p4 * V[0] + ((-5.0 / 143.0) * sqrt2002 * Y7_n6 + (35.0 / 143.0) * sqrt14 * Y7_n2 + (-5.0 / 13.0) * sqrt35 * Y5_n2 * R2 + (7.0 / 11.0) * sqrt5 * Y3_n2 * R4) * V[1] + (1.25 * sqrt35 * Y5_n2 - 3.5 * sqrt5 * Y3_n2 * R2) * V[2] + 5.25 * sqrt5 * Y3_n2 * V[3];
    d[30] = -Y5_n2 * Y4_n1 * V[0] + ((205.0 / 858.0) * sqrt6 * Y7_p1 + (-215.0 / 286.0) * sqrt2 * Y7_p3 + (5.0 / 26.0) * sqrt70 * Y5_p1 * R2 + (-5.0 / 13.0) * sqrt15 * Y5_p3 * R2 + (8.0 / 11.0) * sqrt7 * Y3_p1 * R4 + (-5.0 / 22.0) * sqrt105 * Y3_p3 * R4 + (1.0 / 3.0) * sqrt42 * Y1_p1 * R6) * V[1] + (-0.625 * sqrt70 * Y5_p1 + 1.25 * sqrt15 * Y5_p3 - 4.0 * sqrt7 * Y3_p1 * R2 + 1.25 * sqrt105 * Y3_p3 * R2 - 2.25 * sqrt42 * Y1_p1 * R4) * V[2] + (6.0 * sqrt7 * Y3_p1 - 1.875 * sqrt105 * Y3_p3 + 5.25 * sqrt42 * Y1_p1 * R2) * V[3] - 3.28125 * sqrt42 * Y1_p1 * V[4];
    d[29] = -Y5_n2 * Y4_n2 * V[0] + ((-160.0 / 429.0) * sqrt21 * Y7_p0 + (-50.0 / 143.0) * sqrt11 * Y7_p4 + (-5.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (-5.0 / 13.0) * sqrt15 * Y5_p4 * R2 + (1.0 / 11.0) * sqrt21 * Y3_p0 * R4 + (2.0 / 3.0) * sqrt21 * Y1_p0 * R6) * V[1] + (1.25 * sqrt21 * Y5_p0 + 1.25 * sqrt15 * Y5_p4 - 0.5 * sqrt21 * Y3_p0 * R2 - 4.5 * sqrt21 * Y1_p0 * R4) * V[2] + (0.75 * sqrt21 * Y3_p0 + 10.5 * sqrt21 * Y1_p0 * R2) * V[3] - 6.5625 * sqrt21 * Y1_p0 * V[4];
    d[28] = -Y5_n2 * Y4_n3 * V[0] + ((-175.0 / 858.0) * sqrt42 * Y7_p1 + (-5.0 / 286.0) * sqrt154 * Y7_p5 + (5.0 / 26.0) * sqrt10 * Y5_p1 * R2 + (5.0 / 13.0) * sqrt21 * Y5_p5 * R2 + (49.0 / 22.0) * Y3_p1 * R4 + (-1.0 / 3.0) * sqrt6 * Y1_p1 * R6) * V[1] + (-0.625 * sqrt10 * Y5_p1 - 1.25 * sqrt21 * Y5_p5 - 12.25 * Y3_p1 * R2 + 2.25 * sqrt6 * Y1_p1 * R4) * V[2] + (18.375 * Y3_p1 - 5.25 * sqrt6 * Y1_p1 * R2) * V[3] + 3.28125 * sqrt6 * Y1_p1 * V[4];
    d[27] = -Y5_n2 * Y4_n4 * V[0] + ((-35.0 / 143.0) * sqrt14 * Y7_p2 + (5.0 / 143.0) * sqrt2002 * Y7_p6 + (5.0 / 13.0) * sqrt35 * Y5_p2 * R2 + (-7.0 / 11.0) * sqrt5 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt35 * Y5_p2 + 3.5 * sqrt5 * Y3_p2 * R2) * V[2] - 5.25 * sqrt5 * Y3_p2 * V[3];
    d[22] = -Y5_n3 * Y4_p0 * V[0] + ((-115.0 / 858.0) * sqrt30 * Y7_n3 + (-30.0 / 13.0) * Y5_n3 * R2 + (-15.0 / 11.0) * sqrt7 * Y3_n3 * R4) * V[1] + (7.5 * Y5_n3 + 7.5 * sqrt7 * Y3_n3 * R2) * V[2] - 11.25 * sqrt7 * Y3_n3 * V[3];
    d[23] = -Y5_n3 * Y4_p1 * V[0] + ((5.0 / 66.0) * sqrt33 * Y7_n4 + (815.0 / 1716.0) * sqrt6 * Y7_n2 + (5.0 / 13.0) * sqrt15 * Y5_n2 * R2 + (1.0 / 22.0) * sqrt105 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt15 * Y5_n2 - 0.25 * sqrt105 * Y3_n2 * R2) * V[2] + 0.375 * sqrt105 * Y3_n2 * V[3];
    d[24] = -Y5_n3 * Y4_p2 * V[0] + ((265.0 / 1716.0) * sqrt66 * Y7_n5 + (-125.0 / 132.0) * sqrt2 * Y7_n1 + (30.0 / 13.0) * Y5_n5 * R2 + (9.0 / 22.0) * sqrt21 * Y3_n1 * R4 + (2.0 / 3.0) * sqrt14 * Y1_n1 * R6) * V[1] + (-7.5 * Y5_n5 - 2.25 * sqrt21 * Y3_n1 * R2 - 4.5 * sqrt14 * Y1_n1 * R4) * V[2] + (3.375 * sqrt21 * Y3_n1 + 10.5 * sqrt14 * Y1_n1 * R2) * V[3] - 6.5625 * sqrt14 * Y1_n1 * V[4];
    d[25] = -Y5_n3 * Y4_p3 * V[0] + (25.0 / 1716.0) * sqrt6006 * Y7_n6 * V[1];
    d[26] = -Y5_n3 * Y4_p4 * V[0] + ((-35.0 / 858.0) * sqrt858 * Y7_n7 + (-70.0 / 429.0) * sqrt14 * Y7_n1 + (5.0 / 13.0) * sqrt30 * Y5_n1 * R2 + (-14.0 / 11.0) * sqrt3 * Y3_n1 * R4 + (1.0 / 3.0) * sqrt2 * Y1_n1 * R6) * V[1] + (-1.25 * sqrt30 * Y5_n1 + 7.0 * sqrt3 * Y3_n1 * R2 - 2.25 * sqrt2 * Y1_n1 * R4) * V[2] + (-10.5 * sqrt3 * Y3_n1 + 5.25 * sqrt2 * Y1_n1 * R2) * V[3] - 3.28125 * sqrt2 * Y1_n1 * V[4];
    d[21] = -Y5_n3 * Y4_n1 * V[0] + ((815.0 / 1716.0) * sqrt6 * Y7_p2 + (-5.0 / 66.0) * sqrt33 * Y7_p4 + (5.0 / 13.0) * sqrt15 * Y5_p2 * R2 + (1.0 / 22.0) * sqrt105 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt15 * Y5_p2 - 0.25 * sqrt105 * Y3_p2 * R2) * V[2] + 0.375 * sqrt105 * Y3_p2 * V[3];
    d[20] = -Y5_n3 * Y4_n2 * V[0] + ((-125.0 / 132.0) * sqrt2 * Y7_p1 + (-265.0 / 1716.0) * sqrt66 * Y7_p5 + (-30.0 / 13.0) * Y5_p5 * R2 + (9.0 / 22.0) * sqrt21 * Y3_p1 * R4 + (2.0 / 3.0) * sqrt14 * Y1_p1 * R6) * V[1] + (7.5 * Y5_p5 - 2.25 * sqrt21 * Y3_p1 * R2 - 4.5 * sqrt14 * Y1_p1 * R4) * V[2] + (3.375 * sqrt21 * Y3_p1 + 10.5 * sqrt14 * Y1_p1 * R2) * V[3] - 6.5625 * sqrt14 * Y1_p1 * V[4];
    d[19] = -Y5_n3 * Y4_n3 * V[0] + ((665.0 / 429.0) * Y7_p0 + (-25.0 / 1716.0) * sqrt6006 * Y7_p6 + (-30.0 / 13.0) * Y5_p0 * R2 + (-21.0 / 11.0) * Y3_p0 * R4 + (8.0 / 3.0) * Y1_p0 * R6) * V[1] + (7.5 * Y5_p0 + 10.5 * Y3_p0 * R2 - 18.0 * Y1_p0 * R4) * V[2] + (-15.75 * Y3_p0 + 42.0 * Y1_p0 * R2) * V[3] - 26.25 * Y1_p0 * V[4];
    d[18] = -Y5_n3 * Y4_n4 * V[0] + ((70.0 / 429.0) * sqrt14 * Y7_p1 + (35.0 / 858.0) * sqrt858 * Y7_p7 + (-5.0 / 13.0) * sqrt30 * Y5_p1 * R2 + (14.0 / 11.0) * sqrt3 * Y3_p1 * R4 + (-1.0 / 3.0) * sqrt2 * Y1_p1 * R6) * V[1] + (1.25 * sqrt30 * Y5_p1 - 7.0 * sqrt3 * Y3_p1 * R2 + 2.25 * sqrt2 * Y1_p1 * R4) * V[2] + (10.5 * sqrt3 * Y3_p1 - 5.25 * sqrt2 * Y1_p1 * R2) * V[3] + 3.28125 * sqrt2 * Y1_p1 * V[4];
    d[13] = -Y5_n4 * Y4_p0 * V[0] + ((-20.0 / 143.0) * sqrt165 * Y7_n4 + (-30.0 / 13.0) * Y5_n4 * R2) * V[1] + 7.5 * Y5_n4 * V[2];
    d[14] = -Y5_n4 * Y4_p1 * V[0] + ((-45.0 / 572.0) * sqrt66 * Y7_n5 + (25.0 / 44.0) * sqrt6 * Y7_n3 + (-30.0 / 13.0) * Y5_n5 * R2 + (-9.0 / 22.0) * sqrt35 * Y3_n3 * R4) * V[1] + (7.5 * Y5_n5 + 2.25 * sqrt35 * Y3_n3 * R2) * V[2] - 3.375 * sqrt35 * Y3_n3 * V[3];
    d[15] = -Y5_n4 * Y4_p2 * V[0] + ((5.0 / 286.0) * sqrt858 * Y7_n6 + (-135.0 / 286.0) * sqrt6 * Y7_n2 + (5.0 / 13.0) * sqrt15 * Y5_n2 * R2 + (2.0 / 11.0) * sqrt105 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt15 * Y5_n2 - sqrt105 * Y3_n2 * R2) * V[2] + 1.5 * sqrt105 * Y3_n2 * V[3];
    d[16] = -Y5_n4 * Y4_p3 * V[0] + ((35.0 / 572.0) * sqrt858 * Y7_n7 + (115.0 / 572.0) * sqrt14 * Y7_n1 + (-5.0 / 13.0) * sqrt30 * Y5_n1 * R2 + (7.0 / 22.0) * sqrt3 * Y3_n1 * R4 + 2.0 * sqrt2 * Y1_n1 * R6) * V[1] + (1.25 * sqrt30 * Y5_n1 - 1.75 * sqrt3 * Y3_n1 * R2 - 13.5 * sqrt2 * Y1_n1 * R4) * V[2] + (2.625 * sqrt3 * Y3_n1 + 31.5 * sqrt2 * Y1_n1 * R2) * V[3] - 19.6875 * sqrt2 * Y1_n1 * V[4];
    d[17] = -Y5_n4 * Y4_p4 * V[0];
    d[12] = -Y5_n4 * Y4_n1 * V[0] + ((25.0 / 44.0) * sqrt6 * Y7_p3 + (45.0 / 572.0) * sqrt66 * Y7_p5 + (30.0 / 13.0) * Y5_p5 * R2 + (-9.0 / 22.0) * sqrt35 * Y3_p3 * R4) * V[1] + (-7.5 * Y5_p5 + 2.25 * sqrt35 * Y3_p3 * R2) * V[2] - 3.375 * sqrt35 * Y3_p3 * V[3];
    d[11] = -Y5_n4 * Y4_n2 * V[0] + ((-135.0 / 286.0) * sqrt6 * Y7_p2 + (-5.0 / 286.0) * sqrt858 * Y7_p6 + (5.0 / 13.0) * sqrt15 * Y5_p2 * R2 + (2.0 / 11.0) * sqrt105 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt15 * Y5_p2 - sqrt105 * Y3_p2 * R2) * V[2] + 1.5 * sqrt105 * Y3_p2 * V[3];
    d[10] = -Y5_n4 * Y4_n3 * V[0] + ((115.0 / 572.0) * sqrt14 * Y7_p1 + (-35.0 / 572.0) * sqrt858 * Y7_p7 + (-5.0 / 13.0) * sqrt30 * Y5_p1 * R2 + (7.0 / 22.0) * sqrt3 * Y3_p1 * R4 + 2.0 * sqrt2 * Y1_p1 * R6) * V[1] + (1.25 * sqrt30 * Y5_p1 - 1.75 * sqrt3 * Y3_p1 * R2 - 13.5 * sqrt2 * Y1_p1 * R4) * V[2] + (2.625 * sqrt3 * Y3_p1 + 31.5 * sqrt2 * Y1_p1 * R2) * V[3] - 19.6875 * sqrt2 * Y1_p1 * V[4];
    d[9] = -Y5_n4 * Y4_n4 * V[0] + ((-70.0 / 143.0) * Y7_p0 + (30.0 / 13.0) * Y5_p0 * R2 + (-42.0 / 11.0) * Y3_p0 * R4 + 2.0 * Y1_p0 * R6) * V[1] + (-7.5 * Y5_p0 + 21.0 * Y3_p0 * R2 - 13.5 * Y1_p0 * R4) * V[2] + (-31.5 * Y3_p0 + 31.5 * Y1_p0 * R2) * V[3] - 19.6875 * Y1_p0 * V[4];
    d[4] = -Y5_n5 * Y4_p0 * V[0] + ((-75.0 / 286.0) * sqrt66 * Y7_n5 + (30.0 / 13.0) * Y5_n5 * R2) * V[1] - 7.5 * Y5_n5 * V[2];
    d[5] = -Y5_n5 * Y4_p1 * V[0] + ((-15.0 / 572.0) * sqrt4290 * Y7_n6 + (25.0 / 286.0) * sqrt165 * Y7_n4 + (-30.0 / 13.0) * Y5_n4 * R2) * V[1] + 7.5 * Y5_n4 * V[2];
    d[6] = -Y5_n5 * Y4_p2 * V[0] + ((-5.0 / 572.0) * sqrt30030 * Y7_n7 + (-75.0 / 572.0) * sqrt30 * Y7_n3 + (30.0 / 13.0) * Y5_n3 * R2 + (-15.0 / 22.0) * sqrt7 * Y3_n3 * R4) * V[1] + (-7.5 * Y5_n3 + 3.75 * sqrt7 * Y3_n3 * R2) * V[2] - 5.625 * sqrt7 * Y3_n3 * V[3];
    d[7] = -Y5_n5 * Y4_p3 * V[0] + ((15.0 / 572.0) * sqrt210 * Y7_n2 + (-5.0 / 13.0) * sqrt21 * Y5_n2 * R2 + (35.0 / 22.0) * sqrt3 * Y3_n2 * R4) * V[1] + (1.25 * sqrt21 * Y5_n2 - 8.75 * sqrt3 * Y3_n2 * R2) * V[2] + 13.125 * sqrt3 * Y3_n2 * V[3];
    d[8] = -Y5_n5 * Y4_p4 * V[0] + ((-5.0 / 286.0) * sqrt70 * Y7_n1 + (5.0 / 13.0) * sqrt6 * Y5_n1 * R2 + (-7.0 / 11.0) * sqrt15 * Y3_n1 * R4 + sqrt10 * Y1_n1 * R6) * V[1] + (-1.25 * sqrt6 * Y5_n1 + 3.5 * sqrt15 * Y3_n1 * R2 - 6.75 * sqrt10 * Y1_n1 * R4) * V[2] + (-5.25 * sqrt15 * Y3_n1 + 15.75 * sqrt10 * Y1_n1 * R2) * V[3] - 9.84375 * sqrt10 * Y1_n1 * V[4];
    d[3] = -Y5_n5 * Y4_n1 * V[0] + ((25.0 / 286.0) * sqrt165 * Y7_p4 + (15.0 / 572.0) * sqrt4290 * Y7_p6 + (-30.0 / 13.0) * Y5_p4 * R2) * V[1] + 7.5 * Y5_p4 * V[2];
    d[2] = -Y5_n5 * Y4_n2 * V[0] + ((-75.0 / 572.0) * sqrt30 * Y7_p3 + (5.0 / 572.0) * sqrt30030 * Y7_p7 + (30.0 / 13.0) * Y5_p3 * R2 + (-15.0 / 22.0) * sqrt7 * Y3_p3 * R4) * V[1] + (-7.5 * Y5_p3 + 3.75 * sqrt7 * Y3_p3 * R2) * V[2] - 5.625 * sqrt7 * Y3_p3 * V[3];
    d[1] = -Y5_n5 * Y4_n3 * V[0] + ((15.0 / 572.0) * sqrt210 * Y7_p2 + (-5.0 / 13.0) * sqrt21 * Y5_p2 * R2 + (35.0 / 22.0) * sqrt3 * Y3_p2 * R4) * V[1] + (1.25 * sqrt21 * Y5_p2 - 8.75 * sqrt3 * Y3_p2 * R2) * V[2] + 13.125 * sqrt3 * Y3_p2 * V[3];
    d[0] = -Y5_n5 * Y4_n4 * V[0] + ((-5.0 / 286.0) * sqrt70 * Y7_p1 + (5.0 / 13.0) * sqrt6 * Y5_p1 * R2 + (-7.0 / 11.0) * sqrt15 * Y3_p1 * R4 + sqrt10 * Y1_p1 * R6) * V[1] + (-1.25 * sqrt6 * Y5_p1 + 3.5 * sqrt15 * Y3_p1 * R2 - 6.75 * sqrt10 * Y1_p1 * R4) * V[2] + (-5.25 * sqrt15 * Y3_p1 + 15.75 * sqrt10 * Y1_p1 * R2) * V[3] - 9.84375 * sqrt10 * Y1_p1 * V[4];

    return out;
}

}  // namespace ovlab