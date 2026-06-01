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

#include "OverlapABRecFI.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_f_i(
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
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt22 = std::sqrt(22.0);
    const auto sqrt26 = std::sqrt(26.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt55 = std::sqrt(55.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt77 = std::sqrt(77.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt231 = std::sqrt(231.0);
    const auto sqrt273 = std::sqrt(273.0);
    const auto sqrt286 = std::sqrt(286.0);
    const auto sqrt330 = std::sqrt(330.0);
    const auto sqrt385 = std::sqrt(385.0);
    const auto sqrt462 = std::sqrt(462.0);
    const auto sqrt770 = std::sqrt(770.0);
    const auto sqrt910 = std::sqrt(910.0);
    const auto sqrt1155 = std::sqrt(1155.0);
    const auto sqrt1430 = std::sqrt(1430.0);
    const auto sqrt2310 = std::sqrt(2310.0);
    const auto sqrt30030 = std::sqrt(30030.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^6 · β^3 · p^{-9} · (s|s)
    //   V[1] ↔ α^5 · β^2 · p^{-8} · (s|s)
    //   V[2] ↔ α^4 · β · p^{-7} · (s|s)
    //   V[3] ↔ α^3 · p^{-6} · (s|s)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 4> V = {0.0, 0.0, 0.0, 0.0};

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

            V[0] += cab_ss * alpha6 * beta3 * pinv9;
            V[1] += cab_ss * alpha5 * beta2 * pinv8;
            V[2] += cab_ss * alpha4 * beta * pinv7;
            V[3] += cab_ss * alpha3 * pinv6;
        }
    }

    // ---- Phase 3: fused M·V → 7 × 13 spherical block ----
    const auto Y3_n3 = harm::Y_ll_3_m_n3(AB_x, AB_y, AB_z);
    const auto Y3_n2 = harm::Y_ll_3_m_n2(AB_x, AB_y, AB_z);
    const auto Y3_n1 = harm::Y_ll_3_m_n1(AB_x, AB_y, AB_z);
    const auto Y3_p0 = harm::Y_ll_3_m_p0(AB_x, AB_y, AB_z);
    const auto Y3_p1 = harm::Y_ll_3_m_p1(AB_x, AB_y, AB_z);
    const auto Y3_p2 = harm::Y_ll_3_m_p2(AB_x, AB_y, AB_z);
    const auto Y3_p3 = harm::Y_ll_3_m_p3(AB_x, AB_y, AB_z);
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
    auto *d = buffer;
    d[45] = -Y3_p0 * Y6_p0 * V[0] + ((252.0 / 143.0) * Y7_p0 + (35.0 / 13.0) * Y5_p0 * R2 + (50.0 / 11.0) * Y3_p0 * R4) * V[1] + (-8.75 * Y5_p0 - 25.0 * Y3_p0 * R2) * V[2] + 37.5 * Y3_p0 * V[3];
    d[46] = -Y3_p0 * Y6_p1 * V[0] + ((129.0 / 143.0) * sqrt3 * Y7_p1 + (5.0 / 13.0) * sqrt35 * Y5_p1 * R2 + (25.0 / 22.0) * sqrt14 * Y3_p1 * R4) * V[1] + (-1.25 * sqrt35 * Y5_p1 - 6.25 * sqrt14 * Y3_p1 * R2) * V[2] + 9.375 * sqrt14 * Y3_p1 * V[3];
    d[47] = -Y3_p0 * Y6_p2 * V[0] + ((63.0 / 143.0) * sqrt5 * Y7_p2 + (10.0 / 13.0) * sqrt2 * Y5_p2 * R2 + (10.0 / 11.0) * sqrt14 * Y3_p2 * R4) * V[1] + (-2.5 * sqrt2 * Y5_p2 - 5.0 * sqrt14 * Y3_p2 * R2) * V[2] + 7.5 * sqrt14 * Y3_p2 * V[3];
    d[48] = -Y3_p0 * Y6_p3 * V[0] + ((9.0 / 286.0) * sqrt10 * Y7_p3 + (-5.0 / 13.0) * sqrt3 * Y5_p3 * R2 + (5.0 / 11.0) * sqrt21 * Y3_p3 * R4) * V[1] + (1.25 * sqrt3 * Y5_p3 - 2.5 * sqrt21 * Y3_p3 * R2) * V[2] + 3.75 * sqrt21 * Y3_p3 * V[3];
    d[49] = -Y3_p0 * Y6_p4 * V[0] + ((-24.0 / 143.0) * sqrt33 * Y7_p4 + (-15.0 / 13.0) * sqrt5 * Y5_p4 * R2) * V[1] + 3.75 * sqrt5 * Y5_p4 * V[2];
    d[50] = -Y3_p0 * Y6_p5 * V[0] + ((-21.0 / 26.0) * sqrt6 * Y7_p5 + (-15.0 / 13.0) * sqrt11 * Y5_p5 * R2) * V[1] + 3.75 * sqrt11 * Y5_p5 * V[2];
    d[51] = -Y3_p0 * Y6_p6 * V[0] + (-9.0 / 13.0) * sqrt13 * Y7_p6 * V[1];
    d[44] = -Y3_p0 * Y6_n1 * V[0] + ((129.0 / 143.0) * sqrt3 * Y7_n1 + (5.0 / 13.0) * sqrt35 * Y5_n1 * R2 + (25.0 / 22.0) * sqrt14 * Y3_n1 * R4) * V[1] + (-1.25 * sqrt35 * Y5_n1 - 6.25 * sqrt14 * Y3_n1 * R2) * V[2] + 9.375 * sqrt14 * Y3_n1 * V[3];
    d[43] = -Y3_p0 * Y6_n2 * V[0] + ((63.0 / 143.0) * sqrt5 * Y7_n2 + (10.0 / 13.0) * sqrt2 * Y5_n2 * R2 + (10.0 / 11.0) * sqrt14 * Y3_n2 * R4) * V[1] + (-2.5 * sqrt2 * Y5_n2 - 5.0 * sqrt14 * Y3_n2 * R2) * V[2] + 7.5 * sqrt14 * Y3_n2 * V[3];
    d[42] = -Y3_p0 * Y6_n3 * V[0] + ((9.0 / 286.0) * sqrt10 * Y7_n3 + (-5.0 / 13.0) * sqrt3 * Y5_n3 * R2 + (5.0 / 11.0) * sqrt21 * Y3_n3 * R4) * V[1] + (1.25 * sqrt3 * Y5_n3 - 2.5 * sqrt21 * Y3_n3 * R2) * V[2] + 3.75 * sqrt21 * Y3_n3 * V[3];
    d[41] = -Y3_p0 * Y6_n4 * V[0] + ((-24.0 / 143.0) * sqrt33 * Y7_n4 + (-15.0 / 13.0) * sqrt5 * Y5_n4 * R2) * V[1] + 3.75 * sqrt5 * Y5_n4 * V[2];
    d[40] = -Y3_p0 * Y6_n5 * V[0] + ((-21.0 / 26.0) * sqrt6 * Y7_n5 + (-15.0 / 13.0) * sqrt11 * Y5_n5 * R2) * V[1] + 3.75 * sqrt11 * Y5_n5 * V[2];
    d[39] = -Y3_p0 * Y6_n6 * V[0] + (-9.0 / 13.0) * sqrt13 * Y7_n6 * V[1];
    d[58] = -Y3_p1 * Y6_p0 * V[0] + ((3.0 / 22.0) * sqrt42 * Y7_p1 + (-75.0 / 22.0) * Y3_p1 * R4) * V[1] + 18.75 * Y3_p1 * R2 * V[2] - 28.125 * Y3_p1 * V[3];
    d[59] = -Y3_p1 * Y6_p1 * V[0] + ((-3.0 / 143.0) * sqrt14 * Y7_p0 + (177.0 / 286.0) * sqrt3 * Y7_p2 + (5.0 / 13.0) * sqrt14 * Y5_p0 * R2 + (5.0 / 26.0) * sqrt30 * Y5_p2 * R2 + (25.0 / 22.0) * sqrt14 * Y3_p0 * R4 + (-5.0 / 44.0) * sqrt210 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt14 * Y5_p0 - 0.625 * sqrt30 * Y5_p2 - 6.25 * sqrt14 * Y3_p0 * R2 + 0.625 * sqrt210 * Y3_p2 * R2) * V[2] + (9.375 * sqrt14 * Y3_p0 - 0.9375 * sqrt210 * Y3_p2) * V[3];
    d[60] = -Y3_p1 * Y6_p2 * V[0] + ((69.0 / 286.0) * sqrt5 * Y7_p1 + (48.0 / 143.0) * sqrt15 * Y7_p3 + (5.0 / 13.0) * sqrt21 * Y5_p1 * R2 + (35.0 / 26.0) * sqrt2 * Y5_p3 * R2 + (5.0 / 22.0) * sqrt210 * Y3_p1 * R4 + (-5.0 / 22.0) * sqrt14 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt21 * Y5_p1 - 4.375 * sqrt2 * Y5_p3 - 1.25 * sqrt210 * Y3_p1 * R2 + 1.25 * sqrt14 * Y3_p3 * R2) * V[2] + (1.875 * sqrt210 * Y3_p1 - 1.875 * sqrt14 * Y3_p3) * V[3];
    d[61] = -Y3_p1 * Y6_p3 * V[0] + ((111.0 / 572.0) * sqrt30 * Y7_p2 + (27.0 / 286.0) * sqrt165 * Y7_p4 + (15.0 / 13.0) * sqrt3 * Y5_p2 * R2 + (30.0 / 13.0) * Y5_p4 * R2 + (15.0 / 22.0) * sqrt21 * Y3_p2 * R4) * V[1] + (-3.75 * sqrt3 * Y5_p2 - 7.5 * Y5_p4 - 3.75 * sqrt21 * Y3_p2 * R2) * V[2] + 5.625 * sqrt21 * Y3_p2 * V[3];
    d[62] = -Y3_p1 * Y6_p4 * V[0] + ((573.0 / 572.0) * sqrt2 * Y7_p3 + (87.0 / 572.0) * sqrt22 * Y7_p5 + (5.0 / 13.0) * sqrt15 * Y5_p3 * R2 + (15.0 / 13.0) * sqrt3 * Y5_p5 * R2 + (5.0 / 22.0) * sqrt105 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt15 * Y5_p3 - 3.75 * sqrt3 * Y5_p5 - 1.25 * sqrt105 * Y3_p3 * R2) * V[2] + 1.875 * sqrt105 * Y3_p3 * V[3];
    d[63] = -Y3_p1 * Y6_p5 * V[0] + (1.5 * Y7_p4 + (-3.0 / 52.0) * sqrt26 * Y7_p6) * V[1];
    d[64] = -Y3_p1 * Y6_p6 * V[0] + ((9.0 / 13.0) * sqrt3 * Y7_p5 + (-3.0 / 26.0) * sqrt273 * Y7_p7 + (-15.0 / 26.0) * sqrt22 * Y5_p5 * R2) * V[1] + 1.875 * sqrt22 * Y5_p5 * V[2];
    d[57] = -Y3_p1 * Y6_n1 * V[0] + ((177.0 / 286.0) * sqrt3 * Y7_n2 + (5.0 / 26.0) * sqrt30 * Y5_n2 * R2 + (-5.0 / 44.0) * sqrt210 * Y3_n2 * R4) * V[1] + (-0.625 * sqrt30 * Y5_n2 + 0.625 * sqrt210 * Y3_n2 * R2) * V[2] - 0.9375 * sqrt210 * Y3_n2 * V[3];
    d[56] = -Y3_p1 * Y6_n2 * V[0] + ((48.0 / 143.0) * sqrt15 * Y7_n3 + (69.0 / 286.0) * sqrt5 * Y7_n1 + (35.0 / 26.0) * sqrt2 * Y5_n3 * R2 + (5.0 / 13.0) * sqrt21 * Y5_n1 * R2 + (-5.0 / 22.0) * sqrt14 * Y3_n3 * R4 + (5.0 / 22.0) * sqrt210 * Y3_n1 * R4) * V[1] + (-4.375 * sqrt2 * Y5_n3 - 1.25 * sqrt21 * Y5_n1 + 1.25 * sqrt14 * Y3_n3 * R2 - 1.25 * sqrt210 * Y3_n1 * R2) * V[2] + (-1.875 * sqrt14 * Y3_n3 + 1.875 * sqrt210 * Y3_n1) * V[3];
    d[55] = -Y3_p1 * Y6_n3 * V[0] + ((27.0 / 286.0) * sqrt165 * Y7_n4 + (111.0 / 572.0) * sqrt30 * Y7_n2 + (30.0 / 13.0) * Y5_n4 * R2 + (15.0 / 13.0) * sqrt3 * Y5_n2 * R2 + (15.0 / 22.0) * sqrt21 * Y3_n2 * R4) * V[1] + (-7.5 * Y5_n4 - 3.75 * sqrt3 * Y5_n2 - 3.75 * sqrt21 * Y3_n2 * R2) * V[2] + 5.625 * sqrt21 * Y3_n2 * V[3];
    d[54] = -Y3_p1 * Y6_n4 * V[0] + ((87.0 / 572.0) * sqrt22 * Y7_n5 + (573.0 / 572.0) * sqrt2 * Y7_n3 + (15.0 / 13.0) * sqrt3 * Y5_n5 * R2 + (5.0 / 13.0) * sqrt15 * Y5_n3 * R2 + (5.0 / 22.0) * sqrt105 * Y3_n3 * R4) * V[1] + (-3.75 * sqrt3 * Y5_n5 - 1.25 * sqrt15 * Y5_n3 - 1.25 * sqrt105 * Y3_n3 * R2) * V[2] + 1.875 * sqrt105 * Y3_n3 * V[3];
    d[53] = -Y3_p1 * Y6_n5 * V[0] + ((-3.0 / 52.0) * sqrt26 * Y7_n6 + 1.5 * Y7_n4) * V[1];
    d[52] = -Y3_p1 * Y6_n6 * V[0] + ((-3.0 / 26.0) * sqrt273 * Y7_n7 + (9.0 / 13.0) * sqrt3 * Y7_n5 + (-15.0 / 26.0) * sqrt22 * Y5_n5 * R2) * V[1] + 1.875 * sqrt22 * Y5_n5 * V[2];
    d[71] = -Y3_p2 * Y6_p0 * V[0] + ((-18.0 / 143.0) * sqrt70 * Y7_p2 + (-15.0 / 13.0) * sqrt7 * Y5_p2 * R2 + (15.0 / 11.0) * Y3_p2 * R4) * V[1] + (3.75 * sqrt7 * Y5_p2 - 7.5 * Y3_p2 * R2) * V[2] + 11.25 * Y3_p2 * V[3];
    d[72] = -Y3_p2 * Y6_p1 * V[0] + ((147.0 / 286.0) * sqrt5 * Y7_p1 + (-15.0 / 286.0) * sqrt15 * Y7_p3 + (5.0 / 13.0) * sqrt21 * Y5_p1 * R2 + (-20.0 / 13.0) * sqrt2 * Y5_p3 * R2 + (-5.0 / 44.0) * sqrt210 * Y3_p1 * R4 + (5.0 / 44.0) * sqrt14 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt21 * Y5_p1 + 5.0 * sqrt2 * Y5_p3 + 0.625 * sqrt210 * Y3_p1 * R2 - 0.625 * sqrt14 * Y3_p3 * R2) * V[2] + (-0.9375 * sqrt210 * Y3_p1 + 0.9375 * sqrt14 * Y3_p3) * V[3];
    d[73] = -Y3_p2 * Y6_p2 * V[0] + ((-75.0 / 143.0) * sqrt14 * Y7_p0 + (15.0 / 286.0) * sqrt66 * Y7_p4 + (-5.0 / 13.0) * sqrt14 * Y5_p0 * R2 + (-15.0 / 26.0) * sqrt10 * Y5_p4 * R2 + (10.0 / 11.0) * sqrt14 * Y3_p0 * R4) * V[1] + (1.25 * sqrt14 * Y5_p0 + 1.875 * sqrt10 * Y5_p4 - 5.0 * sqrt14 * Y3_p0 * R2) * V[2] + 7.5 * sqrt14 * Y3_p0 * V[3];
    d[74] = -Y3_p2 * Y6_p3 * V[0] + ((-45.0 / 44.0) * sqrt2 * Y7_p1 + (75.0 / 572.0) * sqrt66 * Y7_p5 + (-15.0 / 13.0) * Y5_p5 * R2 + (15.0 / 22.0) * sqrt21 * Y3_p1 * R4) * V[1] + (3.75 * Y5_p5 - 3.75 * sqrt21 * Y3_p1 * R2) * V[2] + 5.625 * sqrt21 * Y3_p1 * V[3];
    d[75] = -Y3_p2 * Y6_p4 * V[0] + ((-60.0 / 143.0) * sqrt10 * Y7_p2 + (6.0 / 143.0) * sqrt1430 * Y7_p6 + (15.0 / 13.0) * Y5_p2 * R2 + (15.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[1] + (-3.75 * Y5_p2 - 7.5 * sqrt7 * Y3_p2 * R2) * V[2] + 11.25 * sqrt7 * Y3_p2 * V[3];
    d[76] = -Y3_p2 * Y6_p5 * V[0] + ((-57.0 / 572.0) * sqrt110 * Y7_p3 + (3.0 / 52.0) * sqrt910 * Y7_p7 + (5.0 / 13.0) * sqrt33 * Y5_p3 * R2 + (5.0 / 22.0) * sqrt231 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt33 * Y5_p3 - 1.25 * sqrt231 * Y3_p3 * R2) * V[2] + 1.875 * sqrt231 * Y3_p3 * V[3];
    d[77] = -Y3_p2 * Y6_p6 * V[0] + ((-3.0 / 26.0) * sqrt30 * Y7_p4 + (15.0 / 26.0) * sqrt22 * Y5_p4 * R2) * V[1] - 1.875 * sqrt22 * Y5_p4 * V[2];
    d[70] = -Y3_p2 * Y6_n1 * V[0] + ((-15.0 / 286.0) * sqrt15 * Y7_n3 + (-147.0 / 286.0) * sqrt5 * Y7_n1 + (-20.0 / 13.0) * sqrt2 * Y5_n3 * R2 + (-5.0 / 13.0) * sqrt21 * Y5_n1 * R2 + (5.0 / 44.0) * sqrt14 * Y3_n3 * R4 + (5.0 / 44.0) * sqrt210 * Y3_n1 * R4) * V[1] + (5.0 * sqrt2 * Y5_n3 + 1.25 * sqrt21 * Y5_n1 - 0.625 * sqrt14 * Y3_n3 * R2 - 0.625 * sqrt210 * Y3_n1 * R2) * V[2] + (0.9375 * sqrt14 * Y3_n3 + 0.9375 * sqrt210 * Y3_n1) * V[3];
    d[69] = -Y3_p2 * Y6_n2 * V[0] + ((15.0 / 286.0) * sqrt66 * Y7_n4 + (-15.0 / 26.0) * sqrt10 * Y5_n4 * R2) * V[1] + 1.875 * sqrt10 * Y5_n4 * V[2];
    d[68] = -Y3_p2 * Y6_n3 * V[0] + ((75.0 / 572.0) * sqrt66 * Y7_n5 + (-45.0 / 44.0) * sqrt2 * Y7_n1 + (-15.0 / 13.0) * Y5_n5 * R2 + (15.0 / 22.0) * sqrt21 * Y3_n1 * R4) * V[1] + (3.75 * Y5_n5 - 3.75 * sqrt21 * Y3_n1 * R2) * V[2] + 5.625 * sqrt21 * Y3_n1 * V[3];
    d[67] = -Y3_p2 * Y6_n4 * V[0] + ((6.0 / 143.0) * sqrt1430 * Y7_n6 + (-60.0 / 143.0) * sqrt10 * Y7_n2 + (15.0 / 13.0) * Y5_n2 * R2 + (15.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[1] + (-3.75 * Y5_n2 - 7.5 * sqrt7 * Y3_n2 * R2) * V[2] + 11.25 * sqrt7 * Y3_n2 * V[3];
    d[66] = -Y3_p2 * Y6_n5 * V[0] + ((3.0 / 52.0) * sqrt910 * Y7_n7 + (-57.0 / 572.0) * sqrt110 * Y7_n3 + (5.0 / 13.0) * sqrt33 * Y5_n3 * R2 + (5.0 / 22.0) * sqrt231 * Y3_n3 * R4) * V[1] + (-1.25 * sqrt33 * Y5_n3 - 1.25 * sqrt231 * Y3_n3 * R2) * V[2] + 1.875 * sqrt231 * Y3_n3 * V[3];
    d[65] = -Y3_p2 * Y6_n6 * V[0] + ((-3.0 / 26.0) * sqrt30 * Y7_n4 + (15.0 / 26.0) * sqrt22 * Y5_n4 * R2) * V[1] - 1.875 * sqrt22 * Y5_n4 * V[2];
    d[84] = -Y3_p3 * Y6_p0 * V[0] + ((-45.0 / 286.0) * sqrt210 * Y7_p3 + (10.0 / 13.0) * sqrt7 * Y5_p3 * R2 + (-5.0 / 22.0) * Y3_p3 * R4) * V[1] + (-2.5 * sqrt7 * Y5_p3 + 1.25 * Y3_p3 * R2) * V[2] - 1.875 * Y3_p3 * V[3];
    d[85] = -Y3_p3 * Y6_p1 * V[0] + ((189.0 / 286.0) * sqrt5 * Y7_p2 + (-45.0 / 286.0) * sqrt110 * Y7_p4 + (-35.0 / 26.0) * sqrt2 * Y5_p2 * R2 + (5.0 / 13.0) * sqrt6 * Y5_p4 * R2 + (5.0 / 44.0) * sqrt14 * Y3_p2 * R4) * V[1] + (4.375 * sqrt2 * Y5_p2 - 1.25 * sqrt6 * Y5_p4 - 0.625 * sqrt14 * Y3_p2 * R2) * V[2] + 0.9375 * sqrt14 * Y3_p2 * V[3];
    d[86] = -Y3_p3 * Y6_p2 * V[0] + ((-105.0 / 143.0) * sqrt3 * Y7_p1 + (-135.0 / 286.0) * sqrt11 * Y7_p5 + (5.0 / 13.0) * sqrt35 * Y5_p1 * R2 + (5.0 / 26.0) * sqrt6 * Y5_p5 * R2 + (-5.0 / 22.0) * sqrt14 * Y3_p1 * R4) * V[1] + (-1.25 * sqrt35 * Y5_p1 - 0.625 * sqrt6 * Y5_p5 + 1.25 * sqrt14 * Y3_p1 * R2) * V[2] - 1.875 * sqrt14 * Y3_p1 * V[3];
    d[87] = -Y3_p3 * Y6_p3 * V[0] + ((45.0 / 143.0) * sqrt21 * Y7_p0 + (-45.0 / 572.0) * sqrt286 * Y7_p6 + (-10.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (5.0 / 11.0) * sqrt21 * Y3_p0 * R4) * V[1] + (2.5 * sqrt21 * Y5_p0 - 2.5 * sqrt21 * Y3_p0 * R2) * V[2] + 3.75 * sqrt21 * Y3_p0 * V[3];
    d[88] = -Y3_p3 * Y6_p4 * V[0] + ((135.0 / 572.0) * sqrt10 * Y7_p1 + (-3.0 / 572.0) * sqrt30030 * Y7_p7 + (-5.0 / 13.0) * sqrt42 * Y5_p1 * R2 + (5.0 / 22.0) * sqrt105 * Y3_p1 * R4) * V[1] + (1.25 * sqrt42 * Y5_p1 - 1.25 * sqrt105 * Y3_p1 * R2) * V[2] + 1.875 * sqrt105 * Y3_p1 * V[3];
    d[89] = -Y3_p3 * Y6_p5 * V[0] + ((15.0 / 572.0) * sqrt330 * Y7_p2 + (-5.0 / 13.0) * sqrt33 * Y5_p2 * R2 + (5.0 / 22.0) * sqrt231 * Y3_p2 * R4) * V[1] + (1.25 * sqrt33 * Y5_p2 - 1.25 * sqrt231 * Y3_p2 * R2) * V[2] + 1.875 * sqrt231 * Y3_p2 * V[3];
    d[90] = -Y3_p3 * Y6_p6 * V[0] + ((9.0 / 286.0) * sqrt55 * Y7_p3 + (-5.0 / 26.0) * sqrt66 * Y5_p3 * R2 + (5.0 / 22.0) * sqrt462 * Y3_p3 * R4) * V[1] + (0.625 * sqrt66 * Y5_p3 - 1.25 * sqrt462 * Y3_p3 * R2) * V[2] + 1.875 * sqrt462 * Y3_p3 * V[3];
    d[83] = -Y3_p3 * Y6_n1 * V[0] + ((-45.0 / 286.0) * sqrt110 * Y7_n4 + (-189.0 / 286.0) * sqrt5 * Y7_n2 + (5.0 / 13.0) * sqrt6 * Y5_n4 * R2 + (35.0 / 26.0) * sqrt2 * Y5_n2 * R2 + (-5.0 / 44.0) * sqrt14 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt6 * Y5_n4 - 4.375 * sqrt2 * Y5_n2 + 0.625 * sqrt14 * Y3_n2 * R2) * V[2] - 0.9375 * sqrt14 * Y3_n2 * V[3];
    d[82] = -Y3_p3 * Y6_n2 * V[0] + ((-135.0 / 286.0) * sqrt11 * Y7_n5 + (105.0 / 143.0) * sqrt3 * Y7_n1 + (5.0 / 26.0) * sqrt6 * Y5_n5 * R2 + (-5.0 / 13.0) * sqrt35 * Y5_n1 * R2 + (5.0 / 22.0) * sqrt14 * Y3_n1 * R4) * V[1] + (-0.625 * sqrt6 * Y5_n5 + 1.25 * sqrt35 * Y5_n1 - 1.25 * sqrt14 * Y3_n1 * R2) * V[2] + 1.875 * sqrt14 * Y3_n1 * V[3];
    d[81] = -Y3_p3 * Y6_n3 * V[0] + (-45.0 / 572.0) * sqrt286 * Y7_n6 * V[1];
    d[80] = -Y3_p3 * Y6_n4 * V[0] + ((-3.0 / 572.0) * sqrt30030 * Y7_n7 + (135.0 / 572.0) * sqrt10 * Y7_n1 + (-5.0 / 13.0) * sqrt42 * Y5_n1 * R2 + (5.0 / 22.0) * sqrt105 * Y3_n1 * R4) * V[1] + (1.25 * sqrt42 * Y5_n1 - 1.25 * sqrt105 * Y3_n1 * R2) * V[2] + 1.875 * sqrt105 * Y3_n1 * V[3];
    d[79] = -Y3_p3 * Y6_n5 * V[0] + ((15.0 / 572.0) * sqrt330 * Y7_n2 + (-5.0 / 13.0) * sqrt33 * Y5_n2 * R2 + (5.0 / 22.0) * sqrt231 * Y3_n2 * R4) * V[1] + (1.25 * sqrt33 * Y5_n2 - 1.25 * sqrt231 * Y3_n2 * R2) * V[2] + 1.875 * sqrt231 * Y3_n2 * V[3];
    d[78] = -Y3_p3 * Y6_n6 * V[0] + ((9.0 / 286.0) * sqrt55 * Y7_n3 + (-5.0 / 26.0) * sqrt66 * Y5_n3 * R2 + (5.0 / 22.0) * sqrt462 * Y3_n3 * R4) * V[1] + (0.625 * sqrt66 * Y5_n3 - 1.25 * sqrt462 * Y3_n3 * R2) * V[2] + 1.875 * sqrt462 * Y3_n3 * V[3];
    d[32] = -Y3_n1 * Y6_p0 * V[0] + ((3.0 / 22.0) * sqrt42 * Y7_n1 + (-75.0 / 22.0) * Y3_n1 * R4) * V[1] + 18.75 * Y3_n1 * R2 * V[2] - 28.125 * Y3_n1 * V[3];
    d[33] = -Y3_n1 * Y6_p1 * V[0] + ((177.0 / 286.0) * sqrt3 * Y7_n2 + (5.0 / 26.0) * sqrt30 * Y5_n2 * R2 + (-5.0 / 44.0) * sqrt210 * Y3_n2 * R4) * V[1] + (-0.625 * sqrt30 * Y5_n2 + 0.625 * sqrt210 * Y3_n2 * R2) * V[2] - 0.9375 * sqrt210 * Y3_n2 * V[3];
    d[34] = -Y3_n1 * Y6_p2 * V[0] + ((48.0 / 143.0) * sqrt15 * Y7_n3 + (-69.0 / 286.0) * sqrt5 * Y7_n1 + (35.0 / 26.0) * sqrt2 * Y5_n3 * R2 + (-5.0 / 13.0) * sqrt21 * Y5_n1 * R2 + (-5.0 / 22.0) * sqrt14 * Y3_n3 * R4 + (-5.0 / 22.0) * sqrt210 * Y3_n1 * R4) * V[1] + (-4.375 * sqrt2 * Y5_n3 + 1.25 * sqrt21 * Y5_n1 + 1.25 * sqrt14 * Y3_n3 * R2 + 1.25 * sqrt210 * Y3_n1 * R2) * V[2] + (-1.875 * sqrt14 * Y3_n3 - 1.875 * sqrt210 * Y3_n1) * V[3];
    d[35] = -Y3_n1 * Y6_p3 * V[0] + ((27.0 / 286.0) * sqrt165 * Y7_n4 + (-111.0 / 572.0) * sqrt30 * Y7_n2 + (30.0 / 13.0) * Y5_n4 * R2 + (-15.0 / 13.0) * sqrt3 * Y5_n2 * R2 + (-15.0 / 22.0) * sqrt21 * Y3_n2 * R4) * V[1] + (-7.5 * Y5_n4 + 3.75 * sqrt3 * Y5_n2 + 3.75 * sqrt21 * Y3_n2 * R2) * V[2] - 5.625 * sqrt21 * Y3_n2 * V[3];
    d[36] = -Y3_n1 * Y6_p4 * V[0] + ((87.0 / 572.0) * sqrt22 * Y7_n5 + (-573.0 / 572.0) * sqrt2 * Y7_n3 + (15.0 / 13.0) * sqrt3 * Y5_n5 * R2 + (-5.0 / 13.0) * sqrt15 * Y5_n3 * R2 + (-5.0 / 22.0) * sqrt105 * Y3_n3 * R4) * V[1] + (-3.75 * sqrt3 * Y5_n5 + 1.25 * sqrt15 * Y5_n3 + 1.25 * sqrt105 * Y3_n3 * R2) * V[2] - 1.875 * sqrt105 * Y3_n3 * V[3];
    d[37] = -Y3_n1 * Y6_p5 * V[0] + ((-3.0 / 52.0) * sqrt26 * Y7_n6 - 1.5 * Y7_n4) * V[1];
    d[38] = -Y3_n1 * Y6_p6 * V[0] + ((-3.0 / 26.0) * sqrt273 * Y7_n7 + (-9.0 / 13.0) * sqrt3 * Y7_n5 + (15.0 / 26.0) * sqrt22 * Y5_n5 * R2) * V[1] - 1.875 * sqrt22 * Y5_n5 * V[2];
    d[31] = -Y3_n1 * Y6_n1 * V[0] + ((-3.0 / 143.0) * sqrt14 * Y7_p0 + (-177.0 / 286.0) * sqrt3 * Y7_p2 + (5.0 / 13.0) * sqrt14 * Y5_p0 * R2 + (-5.0 / 26.0) * sqrt30 * Y5_p2 * R2 + (25.0 / 22.0) * sqrt14 * Y3_p0 * R4 + (5.0 / 44.0) * sqrt210 * Y3_p2 * R4) * V[1] + (-1.25 * sqrt14 * Y5_p0 + 0.625 * sqrt30 * Y5_p2 - 6.25 * sqrt14 * Y3_p0 * R2 - 0.625 * sqrt210 * Y3_p2 * R2) * V[2] + (9.375 * sqrt14 * Y3_p0 + 0.9375 * sqrt210 * Y3_p2) * V[3];
    d[30] = -Y3_n1 * Y6_n2 * V[0] + ((69.0 / 286.0) * sqrt5 * Y7_p1 + (-48.0 / 143.0) * sqrt15 * Y7_p3 + (5.0 / 13.0) * sqrt21 * Y5_p1 * R2 + (-35.0 / 26.0) * sqrt2 * Y5_p3 * R2 + (5.0 / 22.0) * sqrt210 * Y3_p1 * R4 + (5.0 / 22.0) * sqrt14 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt21 * Y5_p1 + 4.375 * sqrt2 * Y5_p3 - 1.25 * sqrt210 * Y3_p1 * R2 - 1.25 * sqrt14 * Y3_p3 * R2) * V[2] + (1.875 * sqrt210 * Y3_p1 + 1.875 * sqrt14 * Y3_p3) * V[3];
    d[29] = -Y3_n1 * Y6_n3 * V[0] + ((111.0 / 572.0) * sqrt30 * Y7_p2 + (-27.0 / 286.0) * sqrt165 * Y7_p4 + (15.0 / 13.0) * sqrt3 * Y5_p2 * R2 + (-30.0 / 13.0) * Y5_p4 * R2 + (15.0 / 22.0) * sqrt21 * Y3_p2 * R4) * V[1] + (-3.75 * sqrt3 * Y5_p2 + 7.5 * Y5_p4 - 3.75 * sqrt21 * Y3_p2 * R2) * V[2] + 5.625 * sqrt21 * Y3_p2 * V[3];
    d[28] = -Y3_n1 * Y6_n4 * V[0] + ((573.0 / 572.0) * sqrt2 * Y7_p3 + (-87.0 / 572.0) * sqrt22 * Y7_p5 + (5.0 / 13.0) * sqrt15 * Y5_p3 * R2 + (-15.0 / 13.0) * sqrt3 * Y5_p5 * R2 + (5.0 / 22.0) * sqrt105 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt15 * Y5_p3 + 3.75 * sqrt3 * Y5_p5 - 1.25 * sqrt105 * Y3_p3 * R2) * V[2] + 1.875 * sqrt105 * Y3_p3 * V[3];
    d[27] = -Y3_n1 * Y6_n5 * V[0] + (1.5 * Y7_p4 + (3.0 / 52.0) * sqrt26 * Y7_p6) * V[1];
    d[26] = -Y3_n1 * Y6_n6 * V[0] + ((9.0 / 13.0) * sqrt3 * Y7_p5 + (3.0 / 26.0) * sqrt273 * Y7_p7 + (-15.0 / 26.0) * sqrt22 * Y5_p5 * R2) * V[1] + 1.875 * sqrt22 * Y5_p5 * V[2];
    d[19] = -Y3_n2 * Y6_p0 * V[0] + ((-18.0 / 143.0) * sqrt70 * Y7_n2 + (-15.0 / 13.0) * sqrt7 * Y5_n2 * R2 + (15.0 / 11.0) * Y3_n2 * R4) * V[1] + (3.75 * sqrt7 * Y5_n2 - 7.5 * Y3_n2 * R2) * V[2] + 11.25 * Y3_n2 * V[3];
    d[20] = -Y3_n2 * Y6_p1 * V[0] + ((-15.0 / 286.0) * sqrt15 * Y7_n3 + (147.0 / 286.0) * sqrt5 * Y7_n1 + (-20.0 / 13.0) * sqrt2 * Y5_n3 * R2 + (5.0 / 13.0) * sqrt21 * Y5_n1 * R2 + (5.0 / 44.0) * sqrt14 * Y3_n3 * R4 + (-5.0 / 44.0) * sqrt210 * Y3_n1 * R4) * V[1] + (5.0 * sqrt2 * Y5_n3 - 1.25 * sqrt21 * Y5_n1 - 0.625 * sqrt14 * Y3_n3 * R2 + 0.625 * sqrt210 * Y3_n1 * R2) * V[2] + (0.9375 * sqrt14 * Y3_n3 - 0.9375 * sqrt210 * Y3_n1) * V[3];
    d[21] = -Y3_n2 * Y6_p2 * V[0] + ((15.0 / 286.0) * sqrt66 * Y7_n4 + (-15.0 / 26.0) * sqrt10 * Y5_n4 * R2) * V[1] + 1.875 * sqrt10 * Y5_n4 * V[2];
    d[22] = -Y3_n2 * Y6_p3 * V[0] + ((75.0 / 572.0) * sqrt66 * Y7_n5 + (45.0 / 44.0) * sqrt2 * Y7_n1 + (-15.0 / 13.0) * Y5_n5 * R2 + (-15.0 / 22.0) * sqrt21 * Y3_n1 * R4) * V[1] + (3.75 * Y5_n5 + 3.75 * sqrt21 * Y3_n1 * R2) * V[2] - 5.625 * sqrt21 * Y3_n1 * V[3];
    d[23] = -Y3_n2 * Y6_p4 * V[0] + ((6.0 / 143.0) * sqrt1430 * Y7_n6 + (60.0 / 143.0) * sqrt10 * Y7_n2 + (-15.0 / 13.0) * Y5_n2 * R2 + (-15.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[1] + (3.75 * Y5_n2 + 7.5 * sqrt7 * Y3_n2 * R2) * V[2] - 11.25 * sqrt7 * Y3_n2 * V[3];
    d[24] = -Y3_n2 * Y6_p5 * V[0] + ((3.0 / 52.0) * sqrt910 * Y7_n7 + (57.0 / 572.0) * sqrt110 * Y7_n3 + (-5.0 / 13.0) * sqrt33 * Y5_n3 * R2 + (-5.0 / 22.0) * sqrt231 * Y3_n3 * R4) * V[1] + (1.25 * sqrt33 * Y5_n3 + 1.25 * sqrt231 * Y3_n3 * R2) * V[2] - 1.875 * sqrt231 * Y3_n3 * V[3];
    d[25] = -Y3_n2 * Y6_p6 * V[0] + ((3.0 / 26.0) * sqrt30 * Y7_n4 + (-15.0 / 26.0) * sqrt22 * Y5_n4 * R2) * V[1] + 1.875 * sqrt22 * Y5_n4 * V[2];
    d[18] = -Y3_n2 * Y6_n1 * V[0] + ((147.0 / 286.0) * sqrt5 * Y7_p1 + (15.0 / 286.0) * sqrt15 * Y7_p3 + (5.0 / 13.0) * sqrt21 * Y5_p1 * R2 + (20.0 / 13.0) * sqrt2 * Y5_p3 * R2 + (-5.0 / 44.0) * sqrt210 * Y3_p1 * R4 + (-5.0 / 44.0) * sqrt14 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt21 * Y5_p1 - 5.0 * sqrt2 * Y5_p3 + 0.625 * sqrt210 * Y3_p1 * R2 + 0.625 * sqrt14 * Y3_p3 * R2) * V[2] + (-0.9375 * sqrt210 * Y3_p1 - 0.9375 * sqrt14 * Y3_p3) * V[3];
    d[17] = -Y3_n2 * Y6_n2 * V[0] + ((-75.0 / 143.0) * sqrt14 * Y7_p0 + (-15.0 / 286.0) * sqrt66 * Y7_p4 + (-5.0 / 13.0) * sqrt14 * Y5_p0 * R2 + (15.0 / 26.0) * sqrt10 * Y5_p4 * R2 + (10.0 / 11.0) * sqrt14 * Y3_p0 * R4) * V[1] + (1.25 * sqrt14 * Y5_p0 - 1.875 * sqrt10 * Y5_p4 - 5.0 * sqrt14 * Y3_p0 * R2) * V[2] + 7.5 * sqrt14 * Y3_p0 * V[3];
    d[16] = -Y3_n2 * Y6_n3 * V[0] + ((-45.0 / 44.0) * sqrt2 * Y7_p1 + (-75.0 / 572.0) * sqrt66 * Y7_p5 + (15.0 / 13.0) * Y5_p5 * R2 + (15.0 / 22.0) * sqrt21 * Y3_p1 * R4) * V[1] + (-3.75 * Y5_p5 - 3.75 * sqrt21 * Y3_p1 * R2) * V[2] + 5.625 * sqrt21 * Y3_p1 * V[3];
    d[15] = -Y3_n2 * Y6_n4 * V[0] + ((-60.0 / 143.0) * sqrt10 * Y7_p2 + (-6.0 / 143.0) * sqrt1430 * Y7_p6 + (15.0 / 13.0) * Y5_p2 * R2 + (15.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[1] + (-3.75 * Y5_p2 - 7.5 * sqrt7 * Y3_p2 * R2) * V[2] + 11.25 * sqrt7 * Y3_p2 * V[3];
    d[14] = -Y3_n2 * Y6_n5 * V[0] + ((-57.0 / 572.0) * sqrt110 * Y7_p3 + (-3.0 / 52.0) * sqrt910 * Y7_p7 + (5.0 / 13.0) * sqrt33 * Y5_p3 * R2 + (5.0 / 22.0) * sqrt231 * Y3_p3 * R4) * V[1] + (-1.25 * sqrt33 * Y5_p3 - 1.25 * sqrt231 * Y3_p3 * R2) * V[2] + 1.875 * sqrt231 * Y3_p3 * V[3];
    d[13] = -Y3_n2 * Y6_n6 * V[0] + ((-3.0 / 26.0) * sqrt30 * Y7_p4 + (15.0 / 26.0) * sqrt22 * Y5_p4 * R2) * V[1] - 1.875 * sqrt22 * Y5_p4 * V[2];
    d[6] = -Y3_n3 * Y6_p0 * V[0] + ((-45.0 / 286.0) * sqrt210 * Y7_n3 + (10.0 / 13.0) * sqrt7 * Y5_n3 * R2 + (-5.0 / 22.0) * Y3_n3 * R4) * V[1] + (-2.5 * sqrt7 * Y5_n3 + 1.25 * Y3_n3 * R2) * V[2] - 1.875 * Y3_n3 * V[3];
    d[7] = -Y3_n3 * Y6_p1 * V[0] + ((-45.0 / 286.0) * sqrt110 * Y7_n4 + (189.0 / 286.0) * sqrt5 * Y7_n2 + (5.0 / 13.0) * sqrt6 * Y5_n4 * R2 + (-35.0 / 26.0) * sqrt2 * Y5_n2 * R2 + (5.0 / 44.0) * sqrt14 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt6 * Y5_n4 + 4.375 * sqrt2 * Y5_n2 - 0.625 * sqrt14 * Y3_n2 * R2) * V[2] + 0.9375 * sqrt14 * Y3_n2 * V[3];
    d[8] = -Y3_n3 * Y6_p2 * V[0] + ((-135.0 / 286.0) * sqrt11 * Y7_n5 + (-105.0 / 143.0) * sqrt3 * Y7_n1 + (5.0 / 26.0) * sqrt6 * Y5_n5 * R2 + (5.0 / 13.0) * sqrt35 * Y5_n1 * R2 + (-5.0 / 22.0) * sqrt14 * Y3_n1 * R4) * V[1] + (-0.625 * sqrt6 * Y5_n5 - 1.25 * sqrt35 * Y5_n1 + 1.25 * sqrt14 * Y3_n1 * R2) * V[2] - 1.875 * sqrt14 * Y3_n1 * V[3];
    d[9] = -Y3_n3 * Y6_p3 * V[0] + (-45.0 / 572.0) * sqrt286 * Y7_n6 * V[1];
    d[10] = -Y3_n3 * Y6_p4 * V[0] + ((-3.0 / 572.0) * sqrt30030 * Y7_n7 + (-135.0 / 572.0) * sqrt10 * Y7_n1 + (5.0 / 13.0) * sqrt42 * Y5_n1 * R2 + (-5.0 / 22.0) * sqrt105 * Y3_n1 * R4) * V[1] + (-1.25 * sqrt42 * Y5_n1 + 1.25 * sqrt105 * Y3_n1 * R2) * V[2] - 1.875 * sqrt105 * Y3_n1 * V[3];
    d[11] = -Y3_n3 * Y6_p5 * V[0] + ((-15.0 / 572.0) * sqrt330 * Y7_n2 + (5.0 / 13.0) * sqrt33 * Y5_n2 * R2 + (-5.0 / 22.0) * sqrt231 * Y3_n2 * R4) * V[1] + (-1.25 * sqrt33 * Y5_n2 + 1.25 * sqrt231 * Y3_n2 * R2) * V[2] - 1.875 * sqrt231 * Y3_n2 * V[3];
    d[12] = -Y3_n3 * Y6_p6 * V[0] + ((-9.0 / 286.0) * sqrt55 * Y7_n3 + (5.0 / 26.0) * sqrt66 * Y5_n3 * R2 + (-5.0 / 22.0) * sqrt462 * Y3_n3 * R4) * V[1] + (-0.625 * sqrt66 * Y5_n3 + 1.25 * sqrt462 * Y3_n3 * R2) * V[2] - 1.875 * sqrt462 * Y3_n3 * V[3];
    d[5] = -Y3_n3 * Y6_n1 * V[0] + ((189.0 / 286.0) * sqrt5 * Y7_p2 + (45.0 / 286.0) * sqrt110 * Y7_p4 + (-35.0 / 26.0) * sqrt2 * Y5_p2 * R2 + (-5.0 / 13.0) * sqrt6 * Y5_p4 * R2 + (5.0 / 44.0) * sqrt14 * Y3_p2 * R4) * V[1] + (4.375 * sqrt2 * Y5_p2 + 1.25 * sqrt6 * Y5_p4 - 0.625 * sqrt14 * Y3_p2 * R2) * V[2] + 0.9375 * sqrt14 * Y3_p2 * V[3];
    d[4] = -Y3_n3 * Y6_n2 * V[0] + ((-105.0 / 143.0) * sqrt3 * Y7_p1 + (135.0 / 286.0) * sqrt11 * Y7_p5 + (5.0 / 13.0) * sqrt35 * Y5_p1 * R2 + (-5.0 / 26.0) * sqrt6 * Y5_p5 * R2 + (-5.0 / 22.0) * sqrt14 * Y3_p1 * R4) * V[1] + (-1.25 * sqrt35 * Y5_p1 + 0.625 * sqrt6 * Y5_p5 + 1.25 * sqrt14 * Y3_p1 * R2) * V[2] - 1.875 * sqrt14 * Y3_p1 * V[3];
    d[3] = -Y3_n3 * Y6_n3 * V[0] + ((45.0 / 143.0) * sqrt21 * Y7_p0 + (45.0 / 572.0) * sqrt286 * Y7_p6 + (-10.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (5.0 / 11.0) * sqrt21 * Y3_p0 * R4) * V[1] + (2.5 * sqrt21 * Y5_p0 - 2.5 * sqrt21 * Y3_p0 * R2) * V[2] + 3.75 * sqrt21 * Y3_p0 * V[3];
    d[2] = -Y3_n3 * Y6_n4 * V[0] + ((135.0 / 572.0) * sqrt10 * Y7_p1 + (3.0 / 572.0) * sqrt30030 * Y7_p7 + (-5.0 / 13.0) * sqrt42 * Y5_p1 * R2 + (5.0 / 22.0) * sqrt105 * Y3_p1 * R4) * V[1] + (1.25 * sqrt42 * Y5_p1 - 1.25 * sqrt105 * Y3_p1 * R2) * V[2] + 1.875 * sqrt105 * Y3_p1 * V[3];
    d[1] = -Y3_n3 * Y6_n5 * V[0] + ((15.0 / 572.0) * sqrt330 * Y7_p2 + (-5.0 / 13.0) * sqrt33 * Y5_p2 * R2 + (5.0 / 22.0) * sqrt231 * Y3_p2 * R4) * V[1] + (1.25 * sqrt33 * Y5_p2 - 1.25 * sqrt231 * Y3_p2 * R2) * V[2] + 1.875 * sqrt231 * Y3_p2 * V[3];
    d[0] = -Y3_n3 * Y6_n6 * V[0] + ((9.0 / 286.0) * sqrt55 * Y7_p3 + (-5.0 / 26.0) * sqrt66 * Y5_p3 * R2 + (5.0 / 22.0) * sqrt462 * Y3_p3 * R4) * V[1] + (0.625 * sqrt66 * Y5_p3 - 1.25 * sqrt462 * Y3_p3 * R2) * V[2] + 1.875 * sqrt462 * Y3_p3 * V[3];
}

}  // namespace ovlab