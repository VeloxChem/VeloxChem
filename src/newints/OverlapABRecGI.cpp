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

#include "OverlapABRecGI.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_g_i(
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
    const auto sqrt91 = std::sqrt(91.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt143 = std::sqrt(143.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt182 = std::sqrt(182.0);
    const auto sqrt195 = std::sqrt(195.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt231 = std::sqrt(231.0);
    const auto sqrt273 = std::sqrt(273.0);
    const auto sqrt286 = std::sqrt(286.0);
    const auto sqrt330 = std::sqrt(330.0);
    const auto sqrt385 = std::sqrt(385.0);
    const auto sqrt462 = std::sqrt(462.0);
    const auto sqrt546 = std::sqrt(546.0);
    const auto sqrt770 = std::sqrt(770.0);
    const auto sqrt858 = std::sqrt(858.0);
    const auto sqrt1155 = std::sqrt(1155.0);
    const auto sqrt1430 = std::sqrt(1430.0);
    const auto sqrt2310 = std::sqrt(2310.0);
    const auto sqrt2730 = std::sqrt(2730.0);
    const auto sqrt3003 = std::sqrt(3003.0);
    const auto sqrt4290 = std::sqrt(4290.0);
    const auto sqrt15015 = std::sqrt(15015.0);
    const auto sqrt30030 = std::sqrt(30030.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^6 · β^4 · p^{-10} · (s|s)
    //   V[1] ↔ α^5 · β^3 · p^{-9} · (s|s)
    //   V[2] ↔ α^4 · β^2 · p^{-8} · (s|s)
    //   V[3] ↔ α^3 · β · p^{-7} · (s|s)
    //   V[4] ↔ α^2 · p^{-6} · (s|s)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 5> V = {0.0, 0.0, 0.0, 0.0, 0.0};

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

            V[0] += cab_ss * alpha6 * beta4 * pinv10;
            V[1] += cab_ss * alpha5 * beta3 * pinv9;
            V[2] += cab_ss * alpha4 * beta2 * pinv8;
            V[3] += cab_ss * alpha3 * beta * pinv7;
            V[4] += cab_ss * alpha2 * pinv6;
        }
    }

    // ---- Phase 3: fused M·V → 9 × 13 spherical block ----
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
    const auto R4 = R2 * R2;
    const auto R6 = R4 * R2;
    auto *d = buffer;
    d[58] = Y4_p0 * Y6_p0 * V[0] + ((-252.0 / 143.0) * Y8_p0 + (-28.0 / 11.0) * Y6_p0 * R2 + (-450.0 / 143.0) * Y4_p0 * R4 + (-50.0 / 11.0) * Y2_p0 * R6) * V[1] + ((105.0 / 11.0) * Y6_p0 + (225.0 / 11.0) * Y4_p0 * R2 + 37.5 * Y2_p0 * R4) * V[2] + (-37.5 * Y4_p0 - 112.5 * Y2_p0 * R2) * V[3] + 98.4375 * Y2_p0 * V[4];
    d[59] = Y4_p0 * Y6_p1 * V[0] + ((-47.0 / 143.0) * sqrt21 * Y8_p1 + (-64.0 / 33.0) * Y6_p1 * R2 + (-45.0 / 286.0) * sqrt210 * Y4_p1 * R4 + (-50.0 / 33.0) * sqrt7 * Y2_p1 * R6) * V[1] + ((80.0 / 11.0) * Y6_p1 + (45.0 / 44.0) * sqrt210 * Y4_p1 * R2 + 12.5 * sqrt7 * Y2_p1 * R4) * V[2] + (-1.875 * sqrt210 * Y4_p1 - 37.5 * sqrt7 * Y2_p1 * R2) * V[3] + 32.8125 * sqrt7 * Y2_p1 * V[4];
    d[60] = Y4_p0 * Y6_p2 * V[0] + ((-5.0 / 11.0) * sqrt3 * Y8_p2 + (-1.0 / 3.0) * Y6_p2 * R2 + (-10.0 / 33.0) * sqrt70 * Y2_p2 * R6) * V[1] + (1.25 * Y6_p2 + 2.5 * sqrt70 * Y2_p2 * R4) * V[2] - 7.5 * sqrt70 * Y2_p2 * R2 * V[3] + 6.5625 * sqrt70 * Y2_p2 * V[4];
    d[61] = Y4_p0 * Y6_p3 * V[0] + ((15.0 / 286.0) * sqrt22 * Y8_p3 + (18.0 / 11.0) * Y6_p3 * R2 + (225.0 / 143.0) * sqrt3 * Y4_p3 * R4) * V[1] + ((-135.0 / 22.0) * Y6_p3 + (-225.0 / 22.0) * sqrt3 * Y4_p3 * R2) * V[2] + 18.75 * sqrt3 * Y4_p3 * V[3];
    d[62] = Y4_p0 * Y6_p4 * V[0] + ((58.0 / 143.0) * sqrt11 * Y8_p4 + (32.0 / 11.0) * Y6_p4 * R2 + (270.0 / 143.0) * sqrt5 * Y4_p4 * R4) * V[1] + ((-120.0 / 11.0) * Y6_p4 + (-135.0 / 11.0) * sqrt5 * Y4_p4 * R2) * V[2] + 22.5 * sqrt5 * Y4_p4 * V[3];
    d[63] = Y4_p0 * Y6_p5 * V[0] + ((11.0 / 26.0) * sqrt26 * Y8_p5 + 2.0 * Y6_p5 * R2) * V[1] - 7.5 * Y6_p5 * V[2];
    d[64] = Y4_p0 * Y6_p6 * V[0] + ((3.0 / 13.0) * sqrt91 * Y8_p6 - 3.0 * Y6_p6 * R2) * V[1] + 11.25 * Y6_p6 * V[2];
    d[57] = Y4_p0 * Y6_n1 * V[0] + ((-47.0 / 143.0) * sqrt21 * Y8_n1 + (-64.0 / 33.0) * Y6_n1 * R2 + (-45.0 / 286.0) * sqrt210 * Y4_n1 * R4 + (-50.0 / 33.0) * sqrt7 * Y2_n1 * R6) * V[1] + ((80.0 / 11.0) * Y6_n1 + (45.0 / 44.0) * sqrt210 * Y4_n1 * R2 + 12.5 * sqrt7 * Y2_n1 * R4) * V[2] + (-1.875 * sqrt210 * Y4_n1 - 37.5 * sqrt7 * Y2_n1 * R2) * V[3] + 32.8125 * sqrt7 * Y2_n1 * V[4];
    d[56] = Y4_p0 * Y6_n2 * V[0] + ((-5.0 / 11.0) * sqrt3 * Y8_n2 + (-1.0 / 3.0) * Y6_n2 * R2 + (-10.0 / 33.0) * sqrt70 * Y2_n2 * R6) * V[1] + (1.25 * Y6_n2 + 2.5 * sqrt70 * Y2_n2 * R4) * V[2] - 7.5 * sqrt70 * Y2_n2 * R2 * V[3] + 6.5625 * sqrt70 * Y2_n2 * V[4];
    d[55] = Y4_p0 * Y6_n3 * V[0] + ((15.0 / 286.0) * sqrt22 * Y8_n3 + (18.0 / 11.0) * Y6_n3 * R2 + (225.0 / 143.0) * sqrt3 * Y4_n3 * R4) * V[1] + ((-135.0 / 22.0) * Y6_n3 + (-225.0 / 22.0) * sqrt3 * Y4_n3 * R2) * V[2] + 18.75 * sqrt3 * Y4_n3 * V[3];
    d[54] = Y4_p0 * Y6_n4 * V[0] + ((58.0 / 143.0) * sqrt11 * Y8_n4 + (32.0 / 11.0) * Y6_n4 * R2 + (270.0 / 143.0) * sqrt5 * Y4_n4 * R4) * V[1] + ((-120.0 / 11.0) * Y6_n4 + (-135.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[2] + 22.5 * sqrt5 * Y4_n4 * V[3];
    d[53] = Y4_p0 * Y6_n5 * V[0] + ((11.0 / 26.0) * sqrt26 * Y8_n5 + 2.0 * Y6_n5 * R2) * V[1] - 7.5 * Y6_n5 * V[2];
    d[52] = Y4_p0 * Y6_n6 * V[0] + ((3.0 / 13.0) * sqrt91 * Y8_n6 - 3.0 * Y6_n6 * R2) * V[1] + 11.25 * Y6_n6 * V[2];
    d[71] = Y4_p1 * Y6_p0 * V[0] + ((-105.0 / 286.0) * sqrt10 * Y8_p1 + (-2.0 / 33.0) * sqrt210 * Y6_p1 * R2 + (45.0 / 286.0) * Y4_p1 * R4 + (20.0 / 33.0) * sqrt30 * Y2_p1 * R6) * V[1] + ((5.0 / 22.0) * sqrt210 * Y6_p1 + (-45.0 / 44.0) * Y4_p1 * R2 - 5.0 * sqrt30 * Y2_p1 * R4) * V[2] + (1.875 * Y4_p1 + 15.0 * sqrt30 * Y2_p1 * R2) * V[3] - 13.125 * sqrt30 * Y2_p1 * V[4];
    d[72] = Y4_p1 * Y6_p1 * V[0] + ((3.0 / 143.0) * sqrt210 * Y8_p0 + (-199.0 / 286.0) * sqrt3 * Y8_p2 + (-2.0 / 33.0) * sqrt210 * Y6_p0 * R2 + (-53.0 / 33.0) * Y6_p2 * R2 + (-45.0 / 286.0) * sqrt210 * Y4_p0 * R4 + (-135.0 / 572.0) * sqrt42 * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt210 * Y2_p0 * R6 + (5.0 / 33.0) * sqrt70 * Y2_p2 * R6) * V[1] + ((5.0 / 22.0) * sqrt210 * Y6_p0 + (265.0 / 44.0) * Y6_p2 + (45.0 / 44.0) * sqrt210 * Y4_p0 * R2 + (135.0 / 88.0) * sqrt42 * Y4_p2 * R2 + 2.5 * sqrt210 * Y2_p0 * R4 - 1.25 * sqrt70 * Y2_p2 * R4) * V[2] + (-1.875 * sqrt210 * Y4_p0 - 2.8125 * sqrt42 * Y4_p2 - 7.5 * sqrt210 * Y2_p0 * R2 + 3.75 * sqrt70 * Y2_p2 * R2) * V[3] + (6.5625 * sqrt210 * Y2_p0 - 3.28125 * sqrt70 * Y2_p2) * V[4];
    d[73] = Y4_p1 * Y6_p2 * V[0] + ((-29.0 / 286.0) * sqrt21 * Y8_p1 + (-24.0 / 143.0) * sqrt55 * Y8_p3 + (-53.0 / 33.0) * Y6_p1 * R2 + (-13.0 / 22.0) * sqrt10 * Y6_p3 * R2 + (-45.0 / 286.0) * sqrt210 * Y4_p1 * R4 + (-135.0 / 286.0) * sqrt30 * Y4_p3 * R4 + (-40.0 / 33.0) * sqrt7 * Y2_p1 * R6) * V[1] + ((265.0 / 44.0) * Y6_p1 + (195.0 / 88.0) * sqrt10 * Y6_p3 + (45.0 / 44.0) * sqrt210 * Y4_p1 * R2 + (135.0 / 44.0) * sqrt30 * Y4_p3 * R2 + 10.0 * sqrt7 * Y2_p1 * R4) * V[2] + (-1.875 * sqrt210 * Y4_p1 - 5.625 * sqrt30 * Y4_p3 - 30.0 * sqrt7 * Y2_p1 * R2) * V[3] + 26.25 * sqrt7 * Y2_p1 * V[4];
    d[74] = Y4_p1 * Y6_p3 * V[0] + ((-111.0 / 572.0) * sqrt30 * Y8_p2 + (-43.0 / 286.0) * sqrt33 * Y8_p4 + (-13.0 / 22.0) * sqrt10 * Y6_p2 * R2 + (-7.0 / 11.0) * sqrt3 * Y6_p4 * R2 + (-45.0 / 286.0) * sqrt105 * Y4_p2 * R4 + (-90.0 / 143.0) * sqrt15 * Y4_p4 * R4 + (-10.0 / 11.0) * sqrt7 * Y2_p2 * R6) * V[1] + ((195.0 / 88.0) * sqrt10 * Y6_p2 + (105.0 / 44.0) * sqrt3 * Y6_p4 + (45.0 / 44.0) * sqrt105 * Y4_p2 * R2 + (45.0 / 11.0) * sqrt15 * Y4_p4 * R2 + 7.5 * sqrt7 * Y2_p2 * R4) * V[2] + (-1.875 * sqrt105 * Y4_p2 - 7.5 * sqrt15 * Y4_p4 - 22.5 * sqrt7 * Y2_p2 * R2) * V[3] + 19.6875 * sqrt7 * Y2_p2 * V[4];
    d[75] = Y4_p1 * Y6_p4 * V[0] + ((-101.0 / 572.0) * sqrt66 * Y8_p3 + (-1.0 / 572.0) * sqrt1430 * Y8_p5 + (-7.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (1.0 / 11.0) * sqrt55 * Y6_p5 * R2 + (135.0 / 286.0) * Y4_p3 * R4) * V[1] + ((105.0 / 44.0) * sqrt3 * Y6_p3 + (-15.0 / 44.0) * sqrt55 * Y6_p5 + (-135.0 / 44.0) * Y4_p3 * R2) * V[2] + 5.625 * Y4_p3 * V[3];
    d[76] = Y4_p1 * Y6_p5 * V[0] + ((-17.0 / 26.0) * sqrt5 * Y8_p4 + (1.0 / 52.0) * sqrt2730 * Y8_p6 + (1.0 / 11.0) * sqrt55 * Y6_p4 * R2 + 0.5 * sqrt30 * Y6_p6 * R2 + (135.0 / 143.0) * sqrt11 * Y4_p4 * R4) * V[1] + ((-15.0 / 44.0) * sqrt55 * Y6_p4 - 1.875 * sqrt30 * Y6_p6 + (-135.0 / 22.0) * sqrt11 * Y4_p4 * R2) * V[2] + 11.25 * sqrt11 * Y4_p4 * V[3];
    d[77] = Y4_p1 * Y6_p6 * V[0] + ((-1.0 / 13.0) * sqrt195 * Y8_p5 + (3.0 / 26.0) * sqrt273 * Y8_p7 + 0.5 * sqrt30 * Y6_p5 * R2) * V[1] - 1.875 * sqrt30 * Y6_p5 * V[2];
    d[70] = Y4_p1 * Y6_n1 * V[0] + ((-199.0 / 286.0) * sqrt3 * Y8_n2 + (-53.0 / 33.0) * Y6_n2 * R2 + (-135.0 / 572.0) * sqrt42 * Y4_n2 * R4 + (5.0 / 33.0) * sqrt70 * Y2_n2 * R6) * V[1] + ((265.0 / 44.0) * Y6_n2 + (135.0 / 88.0) * sqrt42 * Y4_n2 * R2 - 1.25 * sqrt70 * Y2_n2 * R4) * V[2] + (-2.8125 * sqrt42 * Y4_n2 + 3.75 * sqrt70 * Y2_n2 * R2) * V[3] - 3.28125 * sqrt70 * Y2_n2 * V[4];
    d[69] = Y4_p1 * Y6_n2 * V[0] + ((-24.0 / 143.0) * sqrt55 * Y8_n3 + (-29.0 / 286.0) * sqrt21 * Y8_n1 + (-13.0 / 22.0) * sqrt10 * Y6_n3 * R2 + (-53.0 / 33.0) * Y6_n1 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_n3 * R4 + (-45.0 / 286.0) * sqrt210 * Y4_n1 * R4 + (-40.0 / 33.0) * sqrt7 * Y2_n1 * R6) * V[1] + ((195.0 / 88.0) * sqrt10 * Y6_n3 + (265.0 / 44.0) * Y6_n1 + (135.0 / 44.0) * sqrt30 * Y4_n3 * R2 + (45.0 / 44.0) * sqrt210 * Y4_n1 * R2 + 10.0 * sqrt7 * Y2_n1 * R4) * V[2] + (-5.625 * sqrt30 * Y4_n3 - 1.875 * sqrt210 * Y4_n1 - 30.0 * sqrt7 * Y2_n1 * R2) * V[3] + 26.25 * sqrt7 * Y2_n1 * V[4];
    d[68] = Y4_p1 * Y6_n3 * V[0] + ((-43.0 / 286.0) * sqrt33 * Y8_n4 + (-111.0 / 572.0) * sqrt30 * Y8_n2 + (-7.0 / 11.0) * sqrt3 * Y6_n4 * R2 + (-13.0 / 22.0) * sqrt10 * Y6_n2 * R2 + (-90.0 / 143.0) * sqrt15 * Y4_n4 * R4 + (-45.0 / 286.0) * sqrt105 * Y4_n2 * R4 + (-10.0 / 11.0) * sqrt7 * Y2_n2 * R6) * V[1] + ((105.0 / 44.0) * sqrt3 * Y6_n4 + (195.0 / 88.0) * sqrt10 * Y6_n2 + (45.0 / 11.0) * sqrt15 * Y4_n4 * R2 + (45.0 / 44.0) * sqrt105 * Y4_n2 * R2 + 7.5 * sqrt7 * Y2_n2 * R4) * V[2] + (-7.5 * sqrt15 * Y4_n4 - 1.875 * sqrt105 * Y4_n2 - 22.5 * sqrt7 * Y2_n2 * R2) * V[3] + 19.6875 * sqrt7 * Y2_n2 * V[4];
    d[67] = Y4_p1 * Y6_n4 * V[0] + ((-1.0 / 572.0) * sqrt1430 * Y8_n5 + (-101.0 / 572.0) * sqrt66 * Y8_n3 + (1.0 / 11.0) * sqrt55 * Y6_n5 * R2 + (-7.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (135.0 / 286.0) * Y4_n3 * R4) * V[1] + ((-15.0 / 44.0) * sqrt55 * Y6_n5 + (105.0 / 44.0) * sqrt3 * Y6_n3 + (-135.0 / 44.0) * Y4_n3 * R2) * V[2] + 5.625 * Y4_n3 * V[3];
    d[66] = Y4_p1 * Y6_n5 * V[0] + ((1.0 / 52.0) * sqrt2730 * Y8_n6 + (-17.0 / 26.0) * sqrt5 * Y8_n4 + 0.5 * sqrt30 * Y6_n6 * R2 + (1.0 / 11.0) * sqrt55 * Y6_n4 * R2 + (135.0 / 143.0) * sqrt11 * Y4_n4 * R4) * V[1] + (-1.875 * sqrt30 * Y6_n6 + (-15.0 / 44.0) * sqrt55 * Y6_n4 + (-135.0 / 22.0) * sqrt11 * Y4_n4 * R2) * V[2] + 11.25 * sqrt11 * Y4_n4 * V[3];
    d[65] = Y4_p1 * Y6_n6 * V[0] + ((3.0 / 26.0) * sqrt273 * Y8_n7 + (-1.0 / 13.0) * sqrt195 * Y8_n5 + 0.5 * sqrt30 * Y6_n5 * R2) * V[1] - 1.875 * sqrt30 * Y6_n5 * V[2];
    d[84] = Y4_p2 * Y6_p0 * V[0] + ((12.0 / 143.0) * sqrt14 * Y8_p2 + (1.0 / 3.0) * sqrt42 * Y6_p2 * R2 + (45.0 / 13.0) * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt15 * Y2_p2 * R6) * V[1] + (-1.25 * sqrt42 * Y6_p2 - 22.5 * Y4_p2 * R2 + 2.5 * sqrt15 * Y2_p2 * R4) * V[2] + (41.25 * Y4_p2 - 7.5 * sqrt15 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt15 * Y2_p2 * V[4];
    d[85] = Y4_p2 * Y6_p1 * V[0] + ((-23.0 / 286.0) * sqrt105 * Y8_p1 + (-3.0 / 22.0) * sqrt11 * Y8_p3 + (-28.0 / 33.0) * sqrt5 * Y6_p1 * R2 + (4.0 / 11.0) * sqrt2 * Y6_p3 * R2 + (-135.0 / 572.0) * sqrt42 * Y4_p1 * R4 + (45.0 / 44.0) * sqrt6 * Y4_p3 * R4 + (10.0 / 33.0) * sqrt35 * Y2_p1 * R6) * V[1] + ((35.0 / 11.0) * sqrt5 * Y6_p1 + (-15.0 / 11.0) * sqrt2 * Y6_p3 + (135.0 / 88.0) * sqrt42 * Y4_p1 * R2 + (-585.0 / 88.0) * sqrt6 * Y4_p3 * R2 - 2.5 * sqrt35 * Y2_p1 * R4) * V[2] + (-2.8125 * sqrt42 * Y4_p1 + 12.1875 * sqrt6 * Y4_p3 + 7.5 * sqrt35 * Y2_p1 * R2) * V[3] - 6.5625 * sqrt35 * Y2_p1 * V[4];
    d[86] = Y4_p2 * Y6_p2 * V[0] + ((3.0 / 11.0) * sqrt42 * Y8_p0 + (-37.0 / 286.0) * sqrt66 * Y8_p4 + (1.0 / 3.0) * sqrt42 * Y6_p0 * R2 + (-23.0 / 66.0) * sqrt6 * Y6_p4 * R2 + (45.0 / 143.0) * sqrt30 * Y4_p4 * R4 + (-20.0 / 33.0) * sqrt42 * Y2_p0 * R6) * V[1] + (-1.25 * sqrt42 * Y6_p0 + (115.0 / 88.0) * sqrt6 * Y6_p4 + (-45.0 / 22.0) * sqrt30 * Y4_p4 * R2 + 5.0 * sqrt42 * Y2_p0 * R4) * V[2] + (3.75 * sqrt30 * Y4_p4 - 15.0 * sqrt42 * Y2_p0 * R2) * V[3] + 13.125 * sqrt42 * Y2_p0 * V[4];
    d[87] = Y4_p2 * Y6_p3 * V[0] + ((127.0 / 572.0) * sqrt42 * Y8_p1 + (-27.0 / 572.0) * sqrt858 * Y8_p5 + (4.0 / 11.0) * sqrt2 * Y6_p1 * R2 + (-4.0 / 11.0) * sqrt33 * Y6_p5 * R2 + (-45.0 / 286.0) * sqrt105 * Y4_p1 * R4 + (-10.0 / 11.0) * sqrt14 * Y2_p1 * R6) * V[1] + ((-15.0 / 11.0) * sqrt2 * Y6_p1 + (15.0 / 11.0) * sqrt33 * Y6_p5 + (45.0 / 44.0) * sqrt105 * Y4_p1 * R2 + 7.5 * sqrt14 * Y2_p1 * R4) * V[2] + (-1.875 * sqrt105 * Y4_p1 - 22.5 * sqrt14 * Y2_p1 * R2) * V[3] + 19.6875 * sqrt14 * Y2_p1 * V[4];
    d[88] = Y4_p2 * Y6_p4 * V[0] + ((139.0 / 143.0) * sqrt2 * Y8_p2 + (-1.0 / 143.0) * sqrt30030 * Y8_p6 + (-23.0 / 66.0) * sqrt6 * Y6_p2 * R2 + (-3.0 / 22.0) * sqrt330 * Y6_p6 * R2 + (-135.0 / 143.0) * sqrt7 * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt105 * Y2_p2 * R6) * V[1] + ((115.0 / 88.0) * sqrt6 * Y6_p2 + (45.0 / 88.0) * sqrt330 * Y6_p6 + (135.0 / 22.0) * sqrt7 * Y4_p2 * R2 + 2.5 * sqrt105 * Y2_p2 * R4) * V[2] + (-11.25 * sqrt7 * Y4_p2 - 7.5 * sqrt105 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt105 * Y2_p2 * V[4];
    d[89] = Y4_p2 * Y6_p5 * V[0] + ((23.0 / 52.0) * sqrt6 * Y8_p3 + (-1.0 / 52.0) * sqrt182 * Y8_p7 + (-4.0 / 11.0) * sqrt33 * Y6_p3 * R2 + (-135.0 / 286.0) * sqrt11 * Y4_p3 * R4) * V[1] + ((15.0 / 11.0) * sqrt33 * Y6_p3 + (135.0 / 44.0) * sqrt11 * Y4_p3 * R2) * V[2] - 5.625 * sqrt11 * Y4_p3 * V[3];
    d[90] = Y4_p2 * Y6_p6 * V[0] + ((3.0 / 26.0) * sqrt30 * Y8_p4 + (1.0 / 13.0) * sqrt546 * Y8_p8 + (-3.0 / 22.0) * sqrt330 * Y6_p4 * R2 + (45.0 / 143.0) * sqrt66 * Y4_p4 * R4) * V[1] + ((45.0 / 88.0) * sqrt330 * Y6_p4 + (-45.0 / 22.0) * sqrt66 * Y4_p4 * R2) * V[2] + 3.75 * sqrt66 * Y4_p4 * V[3];
    d[83] = Y4_p2 * Y6_n1 * V[0] + ((-3.0 / 22.0) * sqrt11 * Y8_n3 + (23.0 / 286.0) * sqrt105 * Y8_n1 + (4.0 / 11.0) * sqrt2 * Y6_n3 * R2 + (28.0 / 33.0) * sqrt5 * Y6_n1 * R2 + (45.0 / 44.0) * sqrt6 * Y4_n3 * R4 + (135.0 / 572.0) * sqrt42 * Y4_n1 * R4 + (-10.0 / 33.0) * sqrt35 * Y2_n1 * R6) * V[1] + ((-15.0 / 11.0) * sqrt2 * Y6_n3 + (-35.0 / 11.0) * sqrt5 * Y6_n1 + (-585.0 / 88.0) * sqrt6 * Y4_n3 * R2 + (-135.0 / 88.0) * sqrt42 * Y4_n1 * R2 + 2.5 * sqrt35 * Y2_n1 * R4) * V[2] + (12.1875 * sqrt6 * Y4_n3 + 2.8125 * sqrt42 * Y4_n1 - 7.5 * sqrt35 * Y2_n1 * R2) * V[3] + 6.5625 * sqrt35 * Y2_n1 * V[4];
    d[82] = Y4_p2 * Y6_n2 * V[0] + ((-37.0 / 286.0) * sqrt66 * Y8_n4 + (-23.0 / 66.0) * sqrt6 * Y6_n4 * R2 + (45.0 / 143.0) * sqrt30 * Y4_n4 * R4) * V[1] + ((115.0 / 88.0) * sqrt6 * Y6_n4 + (-45.0 / 22.0) * sqrt30 * Y4_n4 * R2) * V[2] + 3.75 * sqrt30 * Y4_n4 * V[3];
    d[81] = Y4_p2 * Y6_n3 * V[0] + ((-27.0 / 572.0) * sqrt858 * Y8_n5 + (127.0 / 572.0) * sqrt42 * Y8_n1 + (-4.0 / 11.0) * sqrt33 * Y6_n5 * R2 + (4.0 / 11.0) * sqrt2 * Y6_n1 * R2 + (-45.0 / 286.0) * sqrt105 * Y4_n1 * R4 + (-10.0 / 11.0) * sqrt14 * Y2_n1 * R6) * V[1] + ((15.0 / 11.0) * sqrt33 * Y6_n5 + (-15.0 / 11.0) * sqrt2 * Y6_n1 + (45.0 / 44.0) * sqrt105 * Y4_n1 * R2 + 7.5 * sqrt14 * Y2_n1 * R4) * V[2] + (-1.875 * sqrt105 * Y4_n1 - 22.5 * sqrt14 * Y2_n1 * R2) * V[3] + 19.6875 * sqrt14 * Y2_n1 * V[4];
    d[80] = Y4_p2 * Y6_n4 * V[0] + ((-1.0 / 143.0) * sqrt30030 * Y8_n6 + (139.0 / 143.0) * sqrt2 * Y8_n2 + (-3.0 / 22.0) * sqrt330 * Y6_n6 * R2 + (-23.0 / 66.0) * sqrt6 * Y6_n2 * R2 + (-135.0 / 143.0) * sqrt7 * Y4_n2 * R4 + (-10.0 / 33.0) * sqrt105 * Y2_n2 * R6) * V[1] + ((45.0 / 88.0) * sqrt330 * Y6_n6 + (115.0 / 88.0) * sqrt6 * Y6_n2 + (135.0 / 22.0) * sqrt7 * Y4_n2 * R2 + 2.5 * sqrt105 * Y2_n2 * R4) * V[2] + (-11.25 * sqrt7 * Y4_n2 - 7.5 * sqrt105 * Y2_n2 * R2) * V[3] + 6.5625 * sqrt105 * Y2_n2 * V[4];
    d[79] = Y4_p2 * Y6_n5 * V[0] + ((-1.0 / 52.0) * sqrt182 * Y8_n7 + (23.0 / 52.0) * sqrt6 * Y8_n3 + (-4.0 / 11.0) * sqrt33 * Y6_n3 * R2 + (-135.0 / 286.0) * sqrt11 * Y4_n3 * R4) * V[1] + ((15.0 / 11.0) * sqrt33 * Y6_n3 + (135.0 / 44.0) * sqrt11 * Y4_n3 * R2) * V[2] - 5.625 * sqrt11 * Y4_n3 * V[3];
    d[78] = Y4_p2 * Y6_n6 * V[0] + ((1.0 / 13.0) * sqrt546 * Y8_n8 + (3.0 / 26.0) * sqrt30 * Y8_n4 + (-3.0 / 22.0) * sqrt330 * Y6_n4 * R2 + (45.0 / 143.0) * sqrt66 * Y4_n4 * R4) * V[1] + ((45.0 / 88.0) * sqrt330 * Y6_n4 + (-45.0 / 22.0) * sqrt66 * Y4_n4 * R2) * V[2] + 3.75 * sqrt66 * Y4_n4 * V[3];
    d[97] = Y4_p3 * Y6_p0 * V[0] + ((63.0 / 286.0) * sqrt66 * Y8_p3 + (14.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (-765.0 / 286.0) * Y4_p3 * R4) * V[1] + ((-105.0 / 22.0) * sqrt3 * Y6_p3 + (765.0 / 44.0) * Y4_p3 * R2) * V[2] - 31.875 * Y4_p3 * V[3];
    d[98] = Y4_p3 * Y6_p1 * V[0] + ((-7.0 / 22.0) * sqrt21 * Y8_p2 + (5.0 / 286.0) * sqrt2310 * Y8_p4 + (-7.0 / 33.0) * sqrt7 * Y6_p2 * R2 + (5.0 / 33.0) * sqrt210 * Y6_p4 * R2 + (45.0 / 44.0) * sqrt6 * Y4_p2 * R4 + (-45.0 / 286.0) * sqrt42 * Y4_p4 * R4 + (-5.0 / 33.0) * sqrt10 * Y2_p2 * R6) * V[1] + ((35.0 / 44.0) * sqrt7 * Y6_p2 + (-25.0 / 44.0) * sqrt210 * Y6_p4 + (-585.0 / 88.0) * sqrt6 * Y4_p2 * R2 + (45.0 / 44.0) * sqrt42 * Y4_p4 * R2 + 1.25 * sqrt10 * Y2_p2 * R4) * V[2] + (12.1875 * sqrt6 * Y4_p2 - 1.875 * sqrt42 * Y4_p4 - 3.75 * sqrt10 * Y2_p2 * R2) * V[3] + 3.28125 * sqrt10 * Y2_p2 * V[4];
    d[99] = Y4_p3 * Y6_p2 * V[0] + ((119.0 / 143.0) * sqrt3 * Y8_p1 + (1.0 / 286.0) * sqrt3003 * Y8_p5 + (-7.0 / 33.0) * sqrt7 * Y6_p1 * R2 + (7.0 / 66.0) * sqrt462 * Y6_p5 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_p1 * R4 + (40.0 / 33.0) * Y2_p1 * R6) * V[1] + ((35.0 / 44.0) * sqrt7 * Y6_p1 + (-35.0 / 88.0) * sqrt462 * Y6_p5 + (135.0 / 44.0) * sqrt30 * Y4_p1 * R2 - 10.0 * Y2_p1 * R4) * V[2] + (-5.625 * sqrt30 * Y4_p1 + 30.0 * Y2_p1 * R2) * V[3] - 26.25 * Y2_p1 * V[4];
    d[100] = Y4_p3 * Y6_p3 * V[0] + ((-147.0 / 143.0) * sqrt3 * Y8_p0 + (-21.0 / 572.0) * sqrt286 * Y8_p6 + (14.0 / 11.0) * sqrt3 * Y6_p0 * R2 + (3.0 / 22.0) * sqrt154 * Y6_p6 * R2 + (225.0 / 143.0) * sqrt3 * Y4_p0 * R4 + (-20.0 / 11.0) * sqrt3 * Y2_p0 * R6) * V[1] + ((-105.0 / 22.0) * sqrt3 * Y6_p0 + (-45.0 / 88.0) * sqrt154 * Y6_p6 + (-225.0 / 22.0) * sqrt3 * Y4_p0 * R2 + 15.0 * sqrt3 * Y2_p0 * R4) * V[2] + (18.75 * sqrt3 * Y4_p0 - 45.0 * sqrt3 * Y2_p0 * R2) * V[3] + 39.375 * sqrt3 * Y2_p0 * V[4];
    d[101] = Y4_p3 * Y6_p4 * V[0] + ((-175.0 / 572.0) * sqrt10 * Y8_p1 + (-49.0 / 572.0) * sqrt286 * Y8_p7 + (5.0 / 33.0) * sqrt210 * Y6_p1 * R2 + (135.0 / 286.0) * Y4_p1 * R4 + (-20.0 / 33.0) * sqrt30 * Y2_p1 * R6) * V[1] + ((-25.0 / 44.0) * sqrt210 * Y6_p1 + (-135.0 / 44.0) * Y4_p1 * R2 + 5.0 * sqrt30 * Y2_p1 * R4) * V[2] + (5.625 * Y4_p1 - 15.0 * sqrt30 * Y2_p1 * R2) * V[3] + 13.125 * sqrt30 * Y2_p1 * V[4];
    d[102] = Y4_p3 * Y6_p5 * V[0] + ((-29.0 / 572.0) * sqrt154 * Y8_p2 + (-7.0 / 13.0) * sqrt13 * Y8_p8 + (7.0 / 66.0) * sqrt462 * Y6_p2 * R2 + (-135.0 / 286.0) * sqrt11 * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt165 * Y2_p2 * R6) * V[1] + ((-35.0 / 88.0) * sqrt462 * Y6_p2 + (135.0 / 44.0) * sqrt11 * Y4_p2 * R2 + 2.5 * sqrt165 * Y2_p2 * R4) * V[2] + (-5.625 * sqrt11 * Y4_p2 - 7.5 * sqrt165 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt165 * Y2_p2 * V[4];
    d[103] = Y4_p3 * Y6_p6 * V[0] + ((-3.0 / 26.0) * sqrt7 * Y8_p3 + (3.0 / 22.0) * sqrt154 * Y6_p3 * R2 + (-45.0 / 286.0) * sqrt462 * Y4_p3 * R4) * V[1] + ((-45.0 / 88.0) * sqrt154 * Y6_p3 + (45.0 / 44.0) * sqrt462 * Y4_p3 * R2) * V[2] - 1.875 * sqrt462 * Y4_p3 * V[3];
    d[96] = Y4_p3 * Y6_n1 * V[0] + ((5.0 / 286.0) * sqrt2310 * Y8_n4 + (7.0 / 22.0) * sqrt21 * Y8_n2 + (5.0 / 33.0) * sqrt210 * Y6_n4 * R2 + (7.0 / 33.0) * sqrt7 * Y6_n2 * R2 + (-45.0 / 286.0) * sqrt42 * Y4_n4 * R4 + (-45.0 / 44.0) * sqrt6 * Y4_n2 * R4 + (5.0 / 33.0) * sqrt10 * Y2_n2 * R6) * V[1] + ((-25.0 / 44.0) * sqrt210 * Y6_n4 + (-35.0 / 44.0) * sqrt7 * Y6_n2 + (45.0 / 44.0) * sqrt42 * Y4_n4 * R2 + (585.0 / 88.0) * sqrt6 * Y4_n2 * R2 - 1.25 * sqrt10 * Y2_n2 * R4) * V[2] + (-1.875 * sqrt42 * Y4_n4 - 12.1875 * sqrt6 * Y4_n2 + 3.75 * sqrt10 * Y2_n2 * R2) * V[3] - 3.28125 * sqrt10 * Y2_n2 * V[4];
    d[95] = Y4_p3 * Y6_n2 * V[0] + ((1.0 / 286.0) * sqrt3003 * Y8_n5 + (-119.0 / 143.0) * sqrt3 * Y8_n1 + (7.0 / 66.0) * sqrt462 * Y6_n5 * R2 + (7.0 / 33.0) * sqrt7 * Y6_n1 * R2 + (135.0 / 286.0) * sqrt30 * Y4_n1 * R4 + (-40.0 / 33.0) * Y2_n1 * R6) * V[1] + ((-35.0 / 88.0) * sqrt462 * Y6_n5 + (-35.0 / 44.0) * sqrt7 * Y6_n1 + (-135.0 / 44.0) * sqrt30 * Y4_n1 * R2 + 10.0 * Y2_n1 * R4) * V[2] + (5.625 * sqrt30 * Y4_n1 - 30.0 * Y2_n1 * R2) * V[3] + 26.25 * Y2_n1 * V[4];
    d[94] = Y4_p3 * Y6_n3 * V[0] + ((-21.0 / 572.0) * sqrt286 * Y8_n6 + (3.0 / 22.0) * sqrt154 * Y6_n6 * R2) * V[1] + (-45.0 / 88.0) * sqrt154 * Y6_n6 * V[2];
    d[93] = Y4_p3 * Y6_n4 * V[0] + ((-49.0 / 572.0) * sqrt286 * Y8_n7 + (-175.0 / 572.0) * sqrt10 * Y8_n1 + (5.0 / 33.0) * sqrt210 * Y6_n1 * R2 + (135.0 / 286.0) * Y4_n1 * R4 + (-20.0 / 33.0) * sqrt30 * Y2_n1 * R6) * V[1] + ((-25.0 / 44.0) * sqrt210 * Y6_n1 + (-135.0 / 44.0) * Y4_n1 * R2 + 5.0 * sqrt30 * Y2_n1 * R4) * V[2] + (5.625 * Y4_n1 - 15.0 * sqrt30 * Y2_n1 * R2) * V[3] + 13.125 * sqrt30 * Y2_n1 * V[4];
    d[92] = Y4_p3 * Y6_n5 * V[0] + ((-7.0 / 13.0) * sqrt13 * Y8_n8 + (-29.0 / 572.0) * sqrt154 * Y8_n2 + (7.0 / 66.0) * sqrt462 * Y6_n2 * R2 + (-135.0 / 286.0) * sqrt11 * Y4_n2 * R4 + (-10.0 / 33.0) * sqrt165 * Y2_n2 * R6) * V[1] + ((-35.0 / 88.0) * sqrt462 * Y6_n2 + (135.0 / 44.0) * sqrt11 * Y4_n2 * R2 + 2.5 * sqrt165 * Y2_n2 * R4) * V[2] + (-5.625 * sqrt11 * Y4_n2 - 7.5 * sqrt165 * Y2_n2 * R2) * V[3] + 6.5625 * sqrt165 * Y2_n2 * V[4];
    d[91] = Y4_p3 * Y6_n6 * V[0] + ((-3.0 / 26.0) * sqrt7 * Y8_n3 + (3.0 / 22.0) * sqrt154 * Y6_n3 * R2 + (-45.0 / 286.0) * sqrt462 * Y4_n3 * R4) * V[1] + ((-45.0 / 88.0) * sqrt154 * Y6_n3 + (45.0 / 44.0) * sqrt462 * Y4_n3 * R2) * V[2] - 1.875 * sqrt462 * Y4_n3 * V[3];
    d[110] = Y4_p4 * Y6_p0 * V[0] + ((42.0 / 143.0) * sqrt55 * Y8_p4 + (-14.0 / 11.0) * sqrt5 * Y6_p4 * R2 + (90.0 / 143.0) * Y4_p4 * R4) * V[1] + ((105.0 / 22.0) * sqrt5 * Y6_p4 + (-45.0 / 11.0) * Y4_p4 * R2) * V[2] + 7.5 * Y4_p4 * V[3];
    d[111] = Y4_p4 * Y6_p1 * V[0] + ((-21.0 / 143.0) * sqrt77 * Y8_p3 + (2.0 / 143.0) * sqrt15015 * Y8_p5 + (7.0 / 11.0) * sqrt14 * Y6_p3 * R2 + (-1.0 / 33.0) * sqrt2310 * Y6_p5 * R2 + (-45.0 / 286.0) * sqrt42 * Y4_p3 * R4) * V[1] + ((-105.0 / 44.0) * sqrt14 * Y6_p3 + (5.0 / 44.0) * sqrt2310 * Y6_p5 + (45.0 / 44.0) * sqrt42 * Y4_p3 * R2) * V[2] - 1.875 * sqrt42 * Y4_p3 * V[3];
    d[112] = Y4_p4 * Y6_p2 * V[0] + ((14.0 / 143.0) * sqrt105 * Y8_p2 + (21.0 / 143.0) * sqrt143 * Y8_p6 + (-14.0 / 33.0) * sqrt35 * Y6_p2 * R2 + (-1.0 / 11.0) * sqrt77 * Y6_p6 * R2 + (45.0 / 143.0) * sqrt30 * Y4_p2 * R4 + (-5.0 / 33.0) * sqrt2 * Y2_p2 * R6) * V[1] + ((35.0 / 22.0) * sqrt35 * Y6_p2 + (15.0 / 44.0) * sqrt77 * Y6_p6 + (-45.0 / 22.0) * sqrt30 * Y4_p2 * R2 + 1.25 * sqrt2 * Y2_p2 * R4) * V[2] + (3.75 * sqrt30 * Y4_p2 - 3.75 * sqrt2 * Y2_p2 * R2) * V[3] + 3.28125 * sqrt2 * Y2_p2 * V[4];
    d[113] = Y4_p4 * Y6_p3 * V[0] + ((-42.0 / 143.0) * sqrt6 * Y8_p1 + (7.0 / 286.0) * sqrt4290 * Y8_p7 + (7.0 / 11.0) * sqrt14 * Y6_p1 * R2 + (-90.0 / 143.0) * sqrt15 * Y4_p1 * R4 + (5.0 / 11.0) * sqrt2 * Y2_p1 * R6) * V[1] + ((-105.0 / 44.0) * sqrt14 * Y6_p1 + (45.0 / 11.0) * sqrt15 * Y4_p1 * R2 - 3.75 * sqrt2 * Y2_p1 * R4) * V[2] + (-7.5 * sqrt15 * Y4_p1 + 11.25 * sqrt2 * Y2_p1 * R2) * V[3] - 9.84375 * sqrt2 * Y2_p1 * V[4];
    d[114] = Y4_p4 * Y6_p4 * V[0] + ((42.0 / 143.0) * sqrt5 * Y8_p0 + (14.0 / 143.0) * sqrt143 * Y8_p8 + (-14.0 / 11.0) * sqrt5 * Y6_p0 * R2 + (270.0 / 143.0) * sqrt5 * Y4_p0 * R4 + (-10.0 / 11.0) * sqrt5 * Y2_p0 * R6) * V[1] + ((105.0 / 22.0) * sqrt5 * Y6_p0 + (-135.0 / 11.0) * sqrt5 * Y4_p0 * R2 + 7.5 * sqrt5 * Y2_p0 * R4) * V[2] + (22.5 * sqrt5 * Y4_p0 - 22.5 * sqrt5 * Y2_p0 * R2) * V[3] + 19.6875 * sqrt5 * Y2_p0 * V[4];
    d[115] = Y4_p4 * Y6_p5 * V[0] + ((7.0 / 286.0) * sqrt110 * Y8_p1 + (-1.0 / 33.0) * sqrt2310 * Y6_p1 * R2 + (135.0 / 143.0) * sqrt11 * Y4_p1 * R4 + (-5.0 / 33.0) * sqrt330 * Y2_p1 * R6) * V[1] + ((5.0 / 44.0) * sqrt2310 * Y6_p1 + (-135.0 / 22.0) * sqrt11 * Y4_p1 * R2 + 1.25 * sqrt330 * Y2_p1 * R4) * V[2] + (11.25 * sqrt11 * Y4_p1 - 3.75 * sqrt330 * Y2_p1 * R2) * V[3] + 3.28125 * sqrt330 * Y2_p1 * V[4];
    d[116] = Y4_p4 * Y6_p6 * V[0] + ((1.0 / 143.0) * sqrt231 * Y8_p2 + (-1.0 / 11.0) * sqrt77 * Y6_p2 * R2 + (45.0 / 143.0) * sqrt66 * Y4_p2 * R4 + (-5.0 / 11.0) * sqrt110 * Y2_p2 * R6) * V[1] + ((15.0 / 44.0) * sqrt77 * Y6_p2 + (-45.0 / 22.0) * sqrt66 * Y4_p2 * R2 + 3.75 * sqrt110 * Y2_p2 * R4) * V[2] + (3.75 * sqrt66 * Y4_p2 - 11.25 * sqrt110 * Y2_p2 * R2) * V[3] + 9.84375 * sqrt110 * Y2_p2 * V[4];
    d[109] = Y4_p4 * Y6_n1 * V[0] + ((2.0 / 143.0) * sqrt15015 * Y8_n5 + (21.0 / 143.0) * sqrt77 * Y8_n3 + (-1.0 / 33.0) * sqrt2310 * Y6_n5 * R2 + (-7.0 / 11.0) * sqrt14 * Y6_n3 * R2 + (45.0 / 286.0) * sqrt42 * Y4_n3 * R4) * V[1] + ((5.0 / 44.0) * sqrt2310 * Y6_n5 + (105.0 / 44.0) * sqrt14 * Y6_n3 + (-45.0 / 44.0) * sqrt42 * Y4_n3 * R2) * V[2] + 1.875 * sqrt42 * Y4_n3 * V[3];
    d[108] = Y4_p4 * Y6_n2 * V[0] + ((21.0 / 143.0) * sqrt143 * Y8_n6 + (-14.0 / 143.0) * sqrt105 * Y8_n2 + (-1.0 / 11.0) * sqrt77 * Y6_n6 * R2 + (14.0 / 33.0) * sqrt35 * Y6_n2 * R2 + (-45.0 / 143.0) * sqrt30 * Y4_n2 * R4 + (5.0 / 33.0) * sqrt2 * Y2_n2 * R6) * V[1] + ((15.0 / 44.0) * sqrt77 * Y6_n6 + (-35.0 / 22.0) * sqrt35 * Y6_n2 + (45.0 / 22.0) * sqrt30 * Y4_n2 * R2 - 1.25 * sqrt2 * Y2_n2 * R4) * V[2] + (-3.75 * sqrt30 * Y4_n2 + 3.75 * sqrt2 * Y2_n2 * R2) * V[3] - 3.28125 * sqrt2 * Y2_n2 * V[4];
    d[107] = Y4_p4 * Y6_n3 * V[0] + ((7.0 / 286.0) * sqrt4290 * Y8_n7 + (42.0 / 143.0) * sqrt6 * Y8_n1 + (-7.0 / 11.0) * sqrt14 * Y6_n1 * R2 + (90.0 / 143.0) * sqrt15 * Y4_n1 * R4 + (-5.0 / 11.0) * sqrt2 * Y2_n1 * R6) * V[1] + ((105.0 / 44.0) * sqrt14 * Y6_n1 + (-45.0 / 11.0) * sqrt15 * Y4_n1 * R2 + 3.75 * sqrt2 * Y2_n1 * R4) * V[2] + (7.5 * sqrt15 * Y4_n1 - 11.25 * sqrt2 * Y2_n1 * R2) * V[3] + 9.84375 * sqrt2 * Y2_n1 * V[4];
    d[106] = Y4_p4 * Y6_n4 * V[0] + (14.0 / 143.0) * sqrt143 * Y8_n8 * V[1];
    d[105] = Y4_p4 * Y6_n5 * V[0] + ((7.0 / 286.0) * sqrt110 * Y8_n1 + (-1.0 / 33.0) * sqrt2310 * Y6_n1 * R2 + (135.0 / 143.0) * sqrt11 * Y4_n1 * R4 + (-5.0 / 33.0) * sqrt330 * Y2_n1 * R6) * V[1] + ((5.0 / 44.0) * sqrt2310 * Y6_n1 + (-135.0 / 22.0) * sqrt11 * Y4_n1 * R2 + 1.25 * sqrt330 * Y2_n1 * R4) * V[2] + (11.25 * sqrt11 * Y4_n1 - 3.75 * sqrt330 * Y2_n1 * R2) * V[3] + 3.28125 * sqrt330 * Y2_n1 * V[4];
    d[104] = Y4_p4 * Y6_n6 * V[0] + ((1.0 / 143.0) * sqrt231 * Y8_n2 + (-1.0 / 11.0) * sqrt77 * Y6_n2 * R2 + (45.0 / 143.0) * sqrt66 * Y4_n2 * R4 + (-5.0 / 11.0) * sqrt110 * Y2_n2 * R6) * V[1] + ((15.0 / 44.0) * sqrt77 * Y6_n2 + (-45.0 / 22.0) * sqrt66 * Y4_n2 * R2 + 3.75 * sqrt110 * Y2_n2 * R4) * V[2] + (3.75 * sqrt66 * Y4_n2 - 11.25 * sqrt110 * Y2_n2 * R2) * V[3] + 9.84375 * sqrt110 * Y2_n2 * V[4];
    d[45] = Y4_n1 * Y6_p0 * V[0] + ((-105.0 / 286.0) * sqrt10 * Y8_n1 + (-2.0 / 33.0) * sqrt210 * Y6_n1 * R2 + (45.0 / 286.0) * Y4_n1 * R4 + (20.0 / 33.0) * sqrt30 * Y2_n1 * R6) * V[1] + ((5.0 / 22.0) * sqrt210 * Y6_n1 + (-45.0 / 44.0) * Y4_n1 * R2 - 5.0 * sqrt30 * Y2_n1 * R4) * V[2] + (1.875 * Y4_n1 + 15.0 * sqrt30 * Y2_n1 * R2) * V[3] - 13.125 * sqrt30 * Y2_n1 * V[4];
    d[46] = Y4_n1 * Y6_p1 * V[0] + ((-199.0 / 286.0) * sqrt3 * Y8_n2 + (-53.0 / 33.0) * Y6_n2 * R2 + (-135.0 / 572.0) * sqrt42 * Y4_n2 * R4 + (5.0 / 33.0) * sqrt70 * Y2_n2 * R6) * V[1] + ((265.0 / 44.0) * Y6_n2 + (135.0 / 88.0) * sqrt42 * Y4_n2 * R2 - 1.25 * sqrt70 * Y2_n2 * R4) * V[2] + (-2.8125 * sqrt42 * Y4_n2 + 3.75 * sqrt70 * Y2_n2 * R2) * V[3] - 3.28125 * sqrt70 * Y2_n2 * V[4];
    d[47] = Y4_n1 * Y6_p2 * V[0] + ((-24.0 / 143.0) * sqrt55 * Y8_n3 + (29.0 / 286.0) * sqrt21 * Y8_n1 + (-13.0 / 22.0) * sqrt10 * Y6_n3 * R2 + (53.0 / 33.0) * Y6_n1 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_n3 * R4 + (45.0 / 286.0) * sqrt210 * Y4_n1 * R4 + (40.0 / 33.0) * sqrt7 * Y2_n1 * R6) * V[1] + ((195.0 / 88.0) * sqrt10 * Y6_n3 + (-265.0 / 44.0) * Y6_n1 + (135.0 / 44.0) * sqrt30 * Y4_n3 * R2 + (-45.0 / 44.0) * sqrt210 * Y4_n1 * R2 - 10.0 * sqrt7 * Y2_n1 * R4) * V[2] + (-5.625 * sqrt30 * Y4_n3 + 1.875 * sqrt210 * Y4_n1 + 30.0 * sqrt7 * Y2_n1 * R2) * V[3] - 26.25 * sqrt7 * Y2_n1 * V[4];
    d[48] = Y4_n1 * Y6_p3 * V[0] + ((-43.0 / 286.0) * sqrt33 * Y8_n4 + (111.0 / 572.0) * sqrt30 * Y8_n2 + (-7.0 / 11.0) * sqrt3 * Y6_n4 * R2 + (13.0 / 22.0) * sqrt10 * Y6_n2 * R2 + (-90.0 / 143.0) * sqrt15 * Y4_n4 * R4 + (45.0 / 286.0) * sqrt105 * Y4_n2 * R4 + (10.0 / 11.0) * sqrt7 * Y2_n2 * R6) * V[1] + ((105.0 / 44.0) * sqrt3 * Y6_n4 + (-195.0 / 88.0) * sqrt10 * Y6_n2 + (45.0 / 11.0) * sqrt15 * Y4_n4 * R2 + (-45.0 / 44.0) * sqrt105 * Y4_n2 * R2 - 7.5 * sqrt7 * Y2_n2 * R4) * V[2] + (-7.5 * sqrt15 * Y4_n4 + 1.875 * sqrt105 * Y4_n2 + 22.5 * sqrt7 * Y2_n2 * R2) * V[3] - 19.6875 * sqrt7 * Y2_n2 * V[4];
    d[49] = Y4_n1 * Y6_p4 * V[0] + ((-1.0 / 572.0) * sqrt1430 * Y8_n5 + (101.0 / 572.0) * sqrt66 * Y8_n3 + (1.0 / 11.0) * sqrt55 * Y6_n5 * R2 + (7.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (-135.0 / 286.0) * Y4_n3 * R4) * V[1] + ((-15.0 / 44.0) * sqrt55 * Y6_n5 + (-105.0 / 44.0) * sqrt3 * Y6_n3 + (135.0 / 44.0) * Y4_n3 * R2) * V[2] - 5.625 * Y4_n3 * V[3];
    d[50] = Y4_n1 * Y6_p5 * V[0] + ((1.0 / 52.0) * sqrt2730 * Y8_n6 + (17.0 / 26.0) * sqrt5 * Y8_n4 + 0.5 * sqrt30 * Y6_n6 * R2 + (-1.0 / 11.0) * sqrt55 * Y6_n4 * R2 + (-135.0 / 143.0) * sqrt11 * Y4_n4 * R4) * V[1] + (-1.875 * sqrt30 * Y6_n6 + (15.0 / 44.0) * sqrt55 * Y6_n4 + (135.0 / 22.0) * sqrt11 * Y4_n4 * R2) * V[2] - 11.25 * sqrt11 * Y4_n4 * V[3];
    d[51] = Y4_n1 * Y6_p6 * V[0] + ((3.0 / 26.0) * sqrt273 * Y8_n7 + (1.0 / 13.0) * sqrt195 * Y8_n5 - 0.5 * sqrt30 * Y6_n5 * R2) * V[1] + 1.875 * sqrt30 * Y6_n5 * V[2];
    d[44] = Y4_n1 * Y6_n1 * V[0] + ((3.0 / 143.0) * sqrt210 * Y8_p0 + (199.0 / 286.0) * sqrt3 * Y8_p2 + (-2.0 / 33.0) * sqrt210 * Y6_p0 * R2 + (53.0 / 33.0) * Y6_p2 * R2 + (-45.0 / 286.0) * sqrt210 * Y4_p0 * R4 + (135.0 / 572.0) * sqrt42 * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt210 * Y2_p0 * R6 + (-5.0 / 33.0) * sqrt70 * Y2_p2 * R6) * V[1] + ((5.0 / 22.0) * sqrt210 * Y6_p0 + (-265.0 / 44.0) * Y6_p2 + (45.0 / 44.0) * sqrt210 * Y4_p0 * R2 + (-135.0 / 88.0) * sqrt42 * Y4_p2 * R2 + 2.5 * sqrt210 * Y2_p0 * R4 + 1.25 * sqrt70 * Y2_p2 * R4) * V[2] + (-1.875 * sqrt210 * Y4_p0 + 2.8125 * sqrt42 * Y4_p2 - 7.5 * sqrt210 * Y2_p0 * R2 - 3.75 * sqrt70 * Y2_p2 * R2) * V[3] + (6.5625 * sqrt210 * Y2_p0 + 3.28125 * sqrt70 * Y2_p2) * V[4];
    d[43] = Y4_n1 * Y6_n2 * V[0] + ((-29.0 / 286.0) * sqrt21 * Y8_p1 + (24.0 / 143.0) * sqrt55 * Y8_p3 + (-53.0 / 33.0) * Y6_p1 * R2 + (13.0 / 22.0) * sqrt10 * Y6_p3 * R2 + (-45.0 / 286.0) * sqrt210 * Y4_p1 * R4 + (135.0 / 286.0) * sqrt30 * Y4_p3 * R4 + (-40.0 / 33.0) * sqrt7 * Y2_p1 * R6) * V[1] + ((265.0 / 44.0) * Y6_p1 + (-195.0 / 88.0) * sqrt10 * Y6_p3 + (45.0 / 44.0) * sqrt210 * Y4_p1 * R2 + (-135.0 / 44.0) * sqrt30 * Y4_p3 * R2 + 10.0 * sqrt7 * Y2_p1 * R4) * V[2] + (-1.875 * sqrt210 * Y4_p1 + 5.625 * sqrt30 * Y4_p3 - 30.0 * sqrt7 * Y2_p1 * R2) * V[3] + 26.25 * sqrt7 * Y2_p1 * V[4];
    d[42] = Y4_n1 * Y6_n3 * V[0] + ((-111.0 / 572.0) * sqrt30 * Y8_p2 + (43.0 / 286.0) * sqrt33 * Y8_p4 + (-13.0 / 22.0) * sqrt10 * Y6_p2 * R2 + (7.0 / 11.0) * sqrt3 * Y6_p4 * R2 + (-45.0 / 286.0) * sqrt105 * Y4_p2 * R4 + (90.0 / 143.0) * sqrt15 * Y4_p4 * R4 + (-10.0 / 11.0) * sqrt7 * Y2_p2 * R6) * V[1] + ((195.0 / 88.0) * sqrt10 * Y6_p2 + (-105.0 / 44.0) * sqrt3 * Y6_p4 + (45.0 / 44.0) * sqrt105 * Y4_p2 * R2 + (-45.0 / 11.0) * sqrt15 * Y4_p4 * R2 + 7.5 * sqrt7 * Y2_p2 * R4) * V[2] + (-1.875 * sqrt105 * Y4_p2 + 7.5 * sqrt15 * Y4_p4 - 22.5 * sqrt7 * Y2_p2 * R2) * V[3] + 19.6875 * sqrt7 * Y2_p2 * V[4];
    d[41] = Y4_n1 * Y6_n4 * V[0] + ((-101.0 / 572.0) * sqrt66 * Y8_p3 + (1.0 / 572.0) * sqrt1430 * Y8_p5 + (-7.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (-1.0 / 11.0) * sqrt55 * Y6_p5 * R2 + (135.0 / 286.0) * Y4_p3 * R4) * V[1] + ((105.0 / 44.0) * sqrt3 * Y6_p3 + (15.0 / 44.0) * sqrt55 * Y6_p5 + (-135.0 / 44.0) * Y4_p3 * R2) * V[2] + 5.625 * Y4_p3 * V[3];
    d[40] = Y4_n1 * Y6_n5 * V[0] + ((-17.0 / 26.0) * sqrt5 * Y8_p4 + (-1.0 / 52.0) * sqrt2730 * Y8_p6 + (1.0 / 11.0) * sqrt55 * Y6_p4 * R2 - 0.5 * sqrt30 * Y6_p6 * R2 + (135.0 / 143.0) * sqrt11 * Y4_p4 * R4) * V[1] + ((-15.0 / 44.0) * sqrt55 * Y6_p4 + 1.875 * sqrt30 * Y6_p6 + (-135.0 / 22.0) * sqrt11 * Y4_p4 * R2) * V[2] + 11.25 * sqrt11 * Y4_p4 * V[3];
    d[39] = Y4_n1 * Y6_n6 * V[0] + ((-1.0 / 13.0) * sqrt195 * Y8_p5 + (-3.0 / 26.0) * sqrt273 * Y8_p7 + 0.5 * sqrt30 * Y6_p5 * R2) * V[1] - 1.875 * sqrt30 * Y6_p5 * V[2];
    d[32] = Y4_n2 * Y6_p0 * V[0] + ((12.0 / 143.0) * sqrt14 * Y8_n2 + (1.0 / 3.0) * sqrt42 * Y6_n2 * R2 + (45.0 / 13.0) * Y4_n2 * R4 + (-10.0 / 33.0) * sqrt15 * Y2_n2 * R6) * V[1] + (-1.25 * sqrt42 * Y6_n2 - 22.5 * Y4_n2 * R2 + 2.5 * sqrt15 * Y2_n2 * R4) * V[2] + (41.25 * Y4_n2 - 7.5 * sqrt15 * Y2_n2 * R2) * V[3] + 6.5625 * sqrt15 * Y2_n2 * V[4];
    d[33] = Y4_n2 * Y6_p1 * V[0] + ((-3.0 / 22.0) * sqrt11 * Y8_n3 + (-23.0 / 286.0) * sqrt105 * Y8_n1 + (4.0 / 11.0) * sqrt2 * Y6_n3 * R2 + (-28.0 / 33.0) * sqrt5 * Y6_n1 * R2 + (45.0 / 44.0) * sqrt6 * Y4_n3 * R4 + (-135.0 / 572.0) * sqrt42 * Y4_n1 * R4 + (10.0 / 33.0) * sqrt35 * Y2_n1 * R6) * V[1] + ((-15.0 / 11.0) * sqrt2 * Y6_n3 + (35.0 / 11.0) * sqrt5 * Y6_n1 + (-585.0 / 88.0) * sqrt6 * Y4_n3 * R2 + (135.0 / 88.0) * sqrt42 * Y4_n1 * R2 - 2.5 * sqrt35 * Y2_n1 * R4) * V[2] + (12.1875 * sqrt6 * Y4_n3 - 2.8125 * sqrt42 * Y4_n1 + 7.5 * sqrt35 * Y2_n1 * R2) * V[3] - 6.5625 * sqrt35 * Y2_n1 * V[4];
    d[34] = Y4_n2 * Y6_p2 * V[0] + ((-37.0 / 286.0) * sqrt66 * Y8_n4 + (-23.0 / 66.0) * sqrt6 * Y6_n4 * R2 + (45.0 / 143.0) * sqrt30 * Y4_n4 * R4) * V[1] + ((115.0 / 88.0) * sqrt6 * Y6_n4 + (-45.0 / 22.0) * sqrt30 * Y4_n4 * R2) * V[2] + 3.75 * sqrt30 * Y4_n4 * V[3];
    d[35] = Y4_n2 * Y6_p3 * V[0] + ((-27.0 / 572.0) * sqrt858 * Y8_n5 + (-127.0 / 572.0) * sqrt42 * Y8_n1 + (-4.0 / 11.0) * sqrt33 * Y6_n5 * R2 + (-4.0 / 11.0) * sqrt2 * Y6_n1 * R2 + (45.0 / 286.0) * sqrt105 * Y4_n1 * R4 + (10.0 / 11.0) * sqrt14 * Y2_n1 * R6) * V[1] + ((15.0 / 11.0) * sqrt33 * Y6_n5 + (15.0 / 11.0) * sqrt2 * Y6_n1 + (-45.0 / 44.0) * sqrt105 * Y4_n1 * R2 - 7.5 * sqrt14 * Y2_n1 * R4) * V[2] + (1.875 * sqrt105 * Y4_n1 + 22.5 * sqrt14 * Y2_n1 * R2) * V[3] - 19.6875 * sqrt14 * Y2_n1 * V[4];
    d[36] = Y4_n2 * Y6_p4 * V[0] + ((-1.0 / 143.0) * sqrt30030 * Y8_n6 + (-139.0 / 143.0) * sqrt2 * Y8_n2 + (-3.0 / 22.0) * sqrt330 * Y6_n6 * R2 + (23.0 / 66.0) * sqrt6 * Y6_n2 * R2 + (135.0 / 143.0) * sqrt7 * Y4_n2 * R4 + (10.0 / 33.0) * sqrt105 * Y2_n2 * R6) * V[1] + ((45.0 / 88.0) * sqrt330 * Y6_n6 + (-115.0 / 88.0) * sqrt6 * Y6_n2 + (-135.0 / 22.0) * sqrt7 * Y4_n2 * R2 - 2.5 * sqrt105 * Y2_n2 * R4) * V[2] + (11.25 * sqrt7 * Y4_n2 + 7.5 * sqrt105 * Y2_n2 * R2) * V[3] - 6.5625 * sqrt105 * Y2_n2 * V[4];
    d[37] = Y4_n2 * Y6_p5 * V[0] + ((-1.0 / 52.0) * sqrt182 * Y8_n7 + (-23.0 / 52.0) * sqrt6 * Y8_n3 + (4.0 / 11.0) * sqrt33 * Y6_n3 * R2 + (135.0 / 286.0) * sqrt11 * Y4_n3 * R4) * V[1] + ((-15.0 / 11.0) * sqrt33 * Y6_n3 + (-135.0 / 44.0) * sqrt11 * Y4_n3 * R2) * V[2] + 5.625 * sqrt11 * Y4_n3 * V[3];
    d[38] = Y4_n2 * Y6_p6 * V[0] + ((1.0 / 13.0) * sqrt546 * Y8_n8 + (-3.0 / 26.0) * sqrt30 * Y8_n4 + (3.0 / 22.0) * sqrt330 * Y6_n4 * R2 + (-45.0 / 143.0) * sqrt66 * Y4_n4 * R4) * V[1] + ((-45.0 / 88.0) * sqrt330 * Y6_n4 + (45.0 / 22.0) * sqrt66 * Y4_n4 * R2) * V[2] - 3.75 * sqrt66 * Y4_n4 * V[3];
    d[31] = Y4_n2 * Y6_n1 * V[0] + ((-23.0 / 286.0) * sqrt105 * Y8_p1 + (3.0 / 22.0) * sqrt11 * Y8_p3 + (-28.0 / 33.0) * sqrt5 * Y6_p1 * R2 + (-4.0 / 11.0) * sqrt2 * Y6_p3 * R2 + (-135.0 / 572.0) * sqrt42 * Y4_p1 * R4 + (-45.0 / 44.0) * sqrt6 * Y4_p3 * R4 + (10.0 / 33.0) * sqrt35 * Y2_p1 * R6) * V[1] + ((35.0 / 11.0) * sqrt5 * Y6_p1 + (15.0 / 11.0) * sqrt2 * Y6_p3 + (135.0 / 88.0) * sqrt42 * Y4_p1 * R2 + (585.0 / 88.0) * sqrt6 * Y4_p3 * R2 - 2.5 * sqrt35 * Y2_p1 * R4) * V[2] + (-2.8125 * sqrt42 * Y4_p1 - 12.1875 * sqrt6 * Y4_p3 + 7.5 * sqrt35 * Y2_p1 * R2) * V[3] - 6.5625 * sqrt35 * Y2_p1 * V[4];
    d[30] = Y4_n2 * Y6_n2 * V[0] + ((3.0 / 11.0) * sqrt42 * Y8_p0 + (37.0 / 286.0) * sqrt66 * Y8_p4 + (1.0 / 3.0) * sqrt42 * Y6_p0 * R2 + (23.0 / 66.0) * sqrt6 * Y6_p4 * R2 + (-45.0 / 143.0) * sqrt30 * Y4_p4 * R4 + (-20.0 / 33.0) * sqrt42 * Y2_p0 * R6) * V[1] + (-1.25 * sqrt42 * Y6_p0 + (-115.0 / 88.0) * sqrt6 * Y6_p4 + (45.0 / 22.0) * sqrt30 * Y4_p4 * R2 + 5.0 * sqrt42 * Y2_p0 * R4) * V[2] + (-3.75 * sqrt30 * Y4_p4 - 15.0 * sqrt42 * Y2_p0 * R2) * V[3] + 13.125 * sqrt42 * Y2_p0 * V[4];
    d[29] = Y4_n2 * Y6_n3 * V[0] + ((127.0 / 572.0) * sqrt42 * Y8_p1 + (27.0 / 572.0) * sqrt858 * Y8_p5 + (4.0 / 11.0) * sqrt2 * Y6_p1 * R2 + (4.0 / 11.0) * sqrt33 * Y6_p5 * R2 + (-45.0 / 286.0) * sqrt105 * Y4_p1 * R4 + (-10.0 / 11.0) * sqrt14 * Y2_p1 * R6) * V[1] + ((-15.0 / 11.0) * sqrt2 * Y6_p1 + (-15.0 / 11.0) * sqrt33 * Y6_p5 + (45.0 / 44.0) * sqrt105 * Y4_p1 * R2 + 7.5 * sqrt14 * Y2_p1 * R4) * V[2] + (-1.875 * sqrt105 * Y4_p1 - 22.5 * sqrt14 * Y2_p1 * R2) * V[3] + 19.6875 * sqrt14 * Y2_p1 * V[4];
    d[28] = Y4_n2 * Y6_n4 * V[0] + ((139.0 / 143.0) * sqrt2 * Y8_p2 + (1.0 / 143.0) * sqrt30030 * Y8_p6 + (-23.0 / 66.0) * sqrt6 * Y6_p2 * R2 + (3.0 / 22.0) * sqrt330 * Y6_p6 * R2 + (-135.0 / 143.0) * sqrt7 * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt105 * Y2_p2 * R6) * V[1] + ((115.0 / 88.0) * sqrt6 * Y6_p2 + (-45.0 / 88.0) * sqrt330 * Y6_p6 + (135.0 / 22.0) * sqrt7 * Y4_p2 * R2 + 2.5 * sqrt105 * Y2_p2 * R4) * V[2] + (-11.25 * sqrt7 * Y4_p2 - 7.5 * sqrt105 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt105 * Y2_p2 * V[4];
    d[27] = Y4_n2 * Y6_n5 * V[0] + ((23.0 / 52.0) * sqrt6 * Y8_p3 + (1.0 / 52.0) * sqrt182 * Y8_p7 + (-4.0 / 11.0) * sqrt33 * Y6_p3 * R2 + (-135.0 / 286.0) * sqrt11 * Y4_p3 * R4) * V[1] + ((15.0 / 11.0) * sqrt33 * Y6_p3 + (135.0 / 44.0) * sqrt11 * Y4_p3 * R2) * V[2] - 5.625 * sqrt11 * Y4_p3 * V[3];
    d[26] = Y4_n2 * Y6_n6 * V[0] + ((3.0 / 26.0) * sqrt30 * Y8_p4 + (-1.0 / 13.0) * sqrt546 * Y8_p8 + (-3.0 / 22.0) * sqrt330 * Y6_p4 * R2 + (45.0 / 143.0) * sqrt66 * Y4_p4 * R4) * V[1] + ((45.0 / 88.0) * sqrt330 * Y6_p4 + (-45.0 / 22.0) * sqrt66 * Y4_p4 * R2) * V[2] + 3.75 * sqrt66 * Y4_p4 * V[3];
    d[19] = Y4_n3 * Y6_p0 * V[0] + ((63.0 / 286.0) * sqrt66 * Y8_n3 + (14.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (-765.0 / 286.0) * Y4_n3 * R4) * V[1] + ((-105.0 / 22.0) * sqrt3 * Y6_n3 + (765.0 / 44.0) * Y4_n3 * R2) * V[2] - 31.875 * Y4_n3 * V[3];
    d[20] = Y4_n3 * Y6_p1 * V[0] + ((5.0 / 286.0) * sqrt2310 * Y8_n4 + (-7.0 / 22.0) * sqrt21 * Y8_n2 + (5.0 / 33.0) * sqrt210 * Y6_n4 * R2 + (-7.0 / 33.0) * sqrt7 * Y6_n2 * R2 + (-45.0 / 286.0) * sqrt42 * Y4_n4 * R4 + (45.0 / 44.0) * sqrt6 * Y4_n2 * R4 + (-5.0 / 33.0) * sqrt10 * Y2_n2 * R6) * V[1] + ((-25.0 / 44.0) * sqrt210 * Y6_n4 + (35.0 / 44.0) * sqrt7 * Y6_n2 + (45.0 / 44.0) * sqrt42 * Y4_n4 * R2 + (-585.0 / 88.0) * sqrt6 * Y4_n2 * R2 + 1.25 * sqrt10 * Y2_n2 * R4) * V[2] + (-1.875 * sqrt42 * Y4_n4 + 12.1875 * sqrt6 * Y4_n2 - 3.75 * sqrt10 * Y2_n2 * R2) * V[3] + 3.28125 * sqrt10 * Y2_n2 * V[4];
    d[21] = Y4_n3 * Y6_p2 * V[0] + ((1.0 / 286.0) * sqrt3003 * Y8_n5 + (119.0 / 143.0) * sqrt3 * Y8_n1 + (7.0 / 66.0) * sqrt462 * Y6_n5 * R2 + (-7.0 / 33.0) * sqrt7 * Y6_n1 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_n1 * R4 + (40.0 / 33.0) * Y2_n1 * R6) * V[1] + ((-35.0 / 88.0) * sqrt462 * Y6_n5 + (35.0 / 44.0) * sqrt7 * Y6_n1 + (135.0 / 44.0) * sqrt30 * Y4_n1 * R2 - 10.0 * Y2_n1 * R4) * V[2] + (-5.625 * sqrt30 * Y4_n1 + 30.0 * Y2_n1 * R2) * V[3] - 26.25 * Y2_n1 * V[4];
    d[22] = Y4_n3 * Y6_p3 * V[0] + ((-21.0 / 572.0) * sqrt286 * Y8_n6 + (3.0 / 22.0) * sqrt154 * Y6_n6 * R2) * V[1] + (-45.0 / 88.0) * sqrt154 * Y6_n6 * V[2];
    d[23] = Y4_n3 * Y6_p4 * V[0] + ((-49.0 / 572.0) * sqrt286 * Y8_n7 + (175.0 / 572.0) * sqrt10 * Y8_n1 + (-5.0 / 33.0) * sqrt210 * Y6_n1 * R2 + (-135.0 / 286.0) * Y4_n1 * R4 + (20.0 / 33.0) * sqrt30 * Y2_n1 * R6) * V[1] + ((25.0 / 44.0) * sqrt210 * Y6_n1 + (135.0 / 44.0) * Y4_n1 * R2 - 5.0 * sqrt30 * Y2_n1 * R4) * V[2] + (-5.625 * Y4_n1 + 15.0 * sqrt30 * Y2_n1 * R2) * V[3] - 13.125 * sqrt30 * Y2_n1 * V[4];
    d[24] = Y4_n3 * Y6_p5 * V[0] + ((-7.0 / 13.0) * sqrt13 * Y8_n8 + (29.0 / 572.0) * sqrt154 * Y8_n2 + (-7.0 / 66.0) * sqrt462 * Y6_n2 * R2 + (135.0 / 286.0) * sqrt11 * Y4_n2 * R4 + (10.0 / 33.0) * sqrt165 * Y2_n2 * R6) * V[1] + ((35.0 / 88.0) * sqrt462 * Y6_n2 + (-135.0 / 44.0) * sqrt11 * Y4_n2 * R2 - 2.5 * sqrt165 * Y2_n2 * R4) * V[2] + (5.625 * sqrt11 * Y4_n2 + 7.5 * sqrt165 * Y2_n2 * R2) * V[3] - 6.5625 * sqrt165 * Y2_n2 * V[4];
    d[25] = Y4_n3 * Y6_p6 * V[0] + ((3.0 / 26.0) * sqrt7 * Y8_n3 + (-3.0 / 22.0) * sqrt154 * Y6_n3 * R2 + (45.0 / 286.0) * sqrt462 * Y4_n3 * R4) * V[1] + ((45.0 / 88.0) * sqrt154 * Y6_n3 + (-45.0 / 44.0) * sqrt462 * Y4_n3 * R2) * V[2] + 1.875 * sqrt462 * Y4_n3 * V[3];
    d[18] = Y4_n3 * Y6_n1 * V[0] + ((-7.0 / 22.0) * sqrt21 * Y8_p2 + (-5.0 / 286.0) * sqrt2310 * Y8_p4 + (-7.0 / 33.0) * sqrt7 * Y6_p2 * R2 + (-5.0 / 33.0) * sqrt210 * Y6_p4 * R2 + (45.0 / 44.0) * sqrt6 * Y4_p2 * R4 + (45.0 / 286.0) * sqrt42 * Y4_p4 * R4 + (-5.0 / 33.0) * sqrt10 * Y2_p2 * R6) * V[1] + ((35.0 / 44.0) * sqrt7 * Y6_p2 + (25.0 / 44.0) * sqrt210 * Y6_p4 + (-585.0 / 88.0) * sqrt6 * Y4_p2 * R2 + (-45.0 / 44.0) * sqrt42 * Y4_p4 * R2 + 1.25 * sqrt10 * Y2_p2 * R4) * V[2] + (12.1875 * sqrt6 * Y4_p2 + 1.875 * sqrt42 * Y4_p4 - 3.75 * sqrt10 * Y2_p2 * R2) * V[3] + 3.28125 * sqrt10 * Y2_p2 * V[4];
    d[17] = Y4_n3 * Y6_n2 * V[0] + ((119.0 / 143.0) * sqrt3 * Y8_p1 + (-1.0 / 286.0) * sqrt3003 * Y8_p5 + (-7.0 / 33.0) * sqrt7 * Y6_p1 * R2 + (-7.0 / 66.0) * sqrt462 * Y6_p5 * R2 + (-135.0 / 286.0) * sqrt30 * Y4_p1 * R4 + (40.0 / 33.0) * Y2_p1 * R6) * V[1] + ((35.0 / 44.0) * sqrt7 * Y6_p1 + (35.0 / 88.0) * sqrt462 * Y6_p5 + (135.0 / 44.0) * sqrt30 * Y4_p1 * R2 - 10.0 * Y2_p1 * R4) * V[2] + (-5.625 * sqrt30 * Y4_p1 + 30.0 * Y2_p1 * R2) * V[3] - 26.25 * Y2_p1 * V[4];
    d[16] = Y4_n3 * Y6_n3 * V[0] + ((-147.0 / 143.0) * sqrt3 * Y8_p0 + (21.0 / 572.0) * sqrt286 * Y8_p6 + (14.0 / 11.0) * sqrt3 * Y6_p0 * R2 + (-3.0 / 22.0) * sqrt154 * Y6_p6 * R2 + (225.0 / 143.0) * sqrt3 * Y4_p0 * R4 + (-20.0 / 11.0) * sqrt3 * Y2_p0 * R6) * V[1] + ((-105.0 / 22.0) * sqrt3 * Y6_p0 + (45.0 / 88.0) * sqrt154 * Y6_p6 + (-225.0 / 22.0) * sqrt3 * Y4_p0 * R2 + 15.0 * sqrt3 * Y2_p0 * R4) * V[2] + (18.75 * sqrt3 * Y4_p0 - 45.0 * sqrt3 * Y2_p0 * R2) * V[3] + 39.375 * sqrt3 * Y2_p0 * V[4];
    d[15] = Y4_n3 * Y6_n4 * V[0] + ((-175.0 / 572.0) * sqrt10 * Y8_p1 + (49.0 / 572.0) * sqrt286 * Y8_p7 + (5.0 / 33.0) * sqrt210 * Y6_p1 * R2 + (135.0 / 286.0) * Y4_p1 * R4 + (-20.0 / 33.0) * sqrt30 * Y2_p1 * R6) * V[1] + ((-25.0 / 44.0) * sqrt210 * Y6_p1 + (-135.0 / 44.0) * Y4_p1 * R2 + 5.0 * sqrt30 * Y2_p1 * R4) * V[2] + (5.625 * Y4_p1 - 15.0 * sqrt30 * Y2_p1 * R2) * V[3] + 13.125 * sqrt30 * Y2_p1 * V[4];
    d[14] = Y4_n3 * Y6_n5 * V[0] + ((-29.0 / 572.0) * sqrt154 * Y8_p2 + (7.0 / 13.0) * sqrt13 * Y8_p8 + (7.0 / 66.0) * sqrt462 * Y6_p2 * R2 + (-135.0 / 286.0) * sqrt11 * Y4_p2 * R4 + (-10.0 / 33.0) * sqrt165 * Y2_p2 * R6) * V[1] + ((-35.0 / 88.0) * sqrt462 * Y6_p2 + (135.0 / 44.0) * sqrt11 * Y4_p2 * R2 + 2.5 * sqrt165 * Y2_p2 * R4) * V[2] + (-5.625 * sqrt11 * Y4_p2 - 7.5 * sqrt165 * Y2_p2 * R2) * V[3] + 6.5625 * sqrt165 * Y2_p2 * V[4];
    d[13] = Y4_n3 * Y6_n6 * V[0] + ((-3.0 / 26.0) * sqrt7 * Y8_p3 + (3.0 / 22.0) * sqrt154 * Y6_p3 * R2 + (-45.0 / 286.0) * sqrt462 * Y4_p3 * R4) * V[1] + ((-45.0 / 88.0) * sqrt154 * Y6_p3 + (45.0 / 44.0) * sqrt462 * Y4_p3 * R2) * V[2] - 1.875 * sqrt462 * Y4_p3 * V[3];
    d[6] = Y4_n4 * Y6_p0 * V[0] + ((42.0 / 143.0) * sqrt55 * Y8_n4 + (-14.0 / 11.0) * sqrt5 * Y6_n4 * R2 + (90.0 / 143.0) * Y4_n4 * R4) * V[1] + ((105.0 / 22.0) * sqrt5 * Y6_n4 + (-45.0 / 11.0) * Y4_n4 * R2) * V[2] + 7.5 * Y4_n4 * V[3];
    d[7] = Y4_n4 * Y6_p1 * V[0] + ((2.0 / 143.0) * sqrt15015 * Y8_n5 + (-21.0 / 143.0) * sqrt77 * Y8_n3 + (-1.0 / 33.0) * sqrt2310 * Y6_n5 * R2 + (7.0 / 11.0) * sqrt14 * Y6_n3 * R2 + (-45.0 / 286.0) * sqrt42 * Y4_n3 * R4) * V[1] + ((5.0 / 44.0) * sqrt2310 * Y6_n5 + (-105.0 / 44.0) * sqrt14 * Y6_n3 + (45.0 / 44.0) * sqrt42 * Y4_n3 * R2) * V[2] - 1.875 * sqrt42 * Y4_n3 * V[3];
    d[8] = Y4_n4 * Y6_p2 * V[0] + ((21.0 / 143.0) * sqrt143 * Y8_n6 + (14.0 / 143.0) * sqrt105 * Y8_n2 + (-1.0 / 11.0) * sqrt77 * Y6_n6 * R2 + (-14.0 / 33.0) * sqrt35 * Y6_n2 * R2 + (45.0 / 143.0) * sqrt30 * Y4_n2 * R4 + (-5.0 / 33.0) * sqrt2 * Y2_n2 * R6) * V[1] + ((15.0 / 44.0) * sqrt77 * Y6_n6 + (35.0 / 22.0) * sqrt35 * Y6_n2 + (-45.0 / 22.0) * sqrt30 * Y4_n2 * R2 + 1.25 * sqrt2 * Y2_n2 * R4) * V[2] + (3.75 * sqrt30 * Y4_n2 - 3.75 * sqrt2 * Y2_n2 * R2) * V[3] + 3.28125 * sqrt2 * Y2_n2 * V[4];
    d[9] = Y4_n4 * Y6_p3 * V[0] + ((7.0 / 286.0) * sqrt4290 * Y8_n7 + (-42.0 / 143.0) * sqrt6 * Y8_n1 + (7.0 / 11.0) * sqrt14 * Y6_n1 * R2 + (-90.0 / 143.0) * sqrt15 * Y4_n1 * R4 + (5.0 / 11.0) * sqrt2 * Y2_n1 * R6) * V[1] + ((-105.0 / 44.0) * sqrt14 * Y6_n1 + (45.0 / 11.0) * sqrt15 * Y4_n1 * R2 - 3.75 * sqrt2 * Y2_n1 * R4) * V[2] + (-7.5 * sqrt15 * Y4_n1 + 11.25 * sqrt2 * Y2_n1 * R2) * V[3] - 9.84375 * sqrt2 * Y2_n1 * V[4];
    d[10] = Y4_n4 * Y6_p4 * V[0] + (14.0 / 143.0) * sqrt143 * Y8_n8 * V[1];
    d[11] = Y4_n4 * Y6_p5 * V[0] + ((-7.0 / 286.0) * sqrt110 * Y8_n1 + (1.0 / 33.0) * sqrt2310 * Y6_n1 * R2 + (-135.0 / 143.0) * sqrt11 * Y4_n1 * R4 + (5.0 / 33.0) * sqrt330 * Y2_n1 * R6) * V[1] + ((-5.0 / 44.0) * sqrt2310 * Y6_n1 + (135.0 / 22.0) * sqrt11 * Y4_n1 * R2 - 1.25 * sqrt330 * Y2_n1 * R4) * V[2] + (-11.25 * sqrt11 * Y4_n1 + 3.75 * sqrt330 * Y2_n1 * R2) * V[3] - 3.28125 * sqrt330 * Y2_n1 * V[4];
    d[12] = Y4_n4 * Y6_p6 * V[0] + ((-1.0 / 143.0) * sqrt231 * Y8_n2 + (1.0 / 11.0) * sqrt77 * Y6_n2 * R2 + (-45.0 / 143.0) * sqrt66 * Y4_n2 * R4 + (5.0 / 11.0) * sqrt110 * Y2_n2 * R6) * V[1] + ((-15.0 / 44.0) * sqrt77 * Y6_n2 + (45.0 / 22.0) * sqrt66 * Y4_n2 * R2 - 3.75 * sqrt110 * Y2_n2 * R4) * V[2] + (-3.75 * sqrt66 * Y4_n2 + 11.25 * sqrt110 * Y2_n2 * R2) * V[3] - 9.84375 * sqrt110 * Y2_n2 * V[4];
    d[5] = Y4_n4 * Y6_n1 * V[0] + ((-21.0 / 143.0) * sqrt77 * Y8_p3 + (-2.0 / 143.0) * sqrt15015 * Y8_p5 + (7.0 / 11.0) * sqrt14 * Y6_p3 * R2 + (1.0 / 33.0) * sqrt2310 * Y6_p5 * R2 + (-45.0 / 286.0) * sqrt42 * Y4_p3 * R4) * V[1] + ((-105.0 / 44.0) * sqrt14 * Y6_p3 + (-5.0 / 44.0) * sqrt2310 * Y6_p5 + (45.0 / 44.0) * sqrt42 * Y4_p3 * R2) * V[2] - 1.875 * sqrt42 * Y4_p3 * V[3];
    d[4] = Y4_n4 * Y6_n2 * V[0] + ((14.0 / 143.0) * sqrt105 * Y8_p2 + (-21.0 / 143.0) * sqrt143 * Y8_p6 + (-14.0 / 33.0) * sqrt35 * Y6_p2 * R2 + (1.0 / 11.0) * sqrt77 * Y6_p6 * R2 + (45.0 / 143.0) * sqrt30 * Y4_p2 * R4 + (-5.0 / 33.0) * sqrt2 * Y2_p2 * R6) * V[1] + ((35.0 / 22.0) * sqrt35 * Y6_p2 + (-15.0 / 44.0) * sqrt77 * Y6_p6 + (-45.0 / 22.0) * sqrt30 * Y4_p2 * R2 + 1.25 * sqrt2 * Y2_p2 * R4) * V[2] + (3.75 * sqrt30 * Y4_p2 - 3.75 * sqrt2 * Y2_p2 * R2) * V[3] + 3.28125 * sqrt2 * Y2_p2 * V[4];
    d[3] = Y4_n4 * Y6_n3 * V[0] + ((-42.0 / 143.0) * sqrt6 * Y8_p1 + (-7.0 / 286.0) * sqrt4290 * Y8_p7 + (7.0 / 11.0) * sqrt14 * Y6_p1 * R2 + (-90.0 / 143.0) * sqrt15 * Y4_p1 * R4 + (5.0 / 11.0) * sqrt2 * Y2_p1 * R6) * V[1] + ((-105.0 / 44.0) * sqrt14 * Y6_p1 + (45.0 / 11.0) * sqrt15 * Y4_p1 * R2 - 3.75 * sqrt2 * Y2_p1 * R4) * V[2] + (-7.5 * sqrt15 * Y4_p1 + 11.25 * sqrt2 * Y2_p1 * R2) * V[3] - 9.84375 * sqrt2 * Y2_p1 * V[4];
    d[2] = Y4_n4 * Y6_n4 * V[0] + ((42.0 / 143.0) * sqrt5 * Y8_p0 + (-14.0 / 143.0) * sqrt143 * Y8_p8 + (-14.0 / 11.0) * sqrt5 * Y6_p0 * R2 + (270.0 / 143.0) * sqrt5 * Y4_p0 * R4 + (-10.0 / 11.0) * sqrt5 * Y2_p0 * R6) * V[1] + ((105.0 / 22.0) * sqrt5 * Y6_p0 + (-135.0 / 11.0) * sqrt5 * Y4_p0 * R2 + 7.5 * sqrt5 * Y2_p0 * R4) * V[2] + (22.5 * sqrt5 * Y4_p0 - 22.5 * sqrt5 * Y2_p0 * R2) * V[3] + 19.6875 * sqrt5 * Y2_p0 * V[4];
    d[1] = Y4_n4 * Y6_n5 * V[0] + ((7.0 / 286.0) * sqrt110 * Y8_p1 + (-1.0 / 33.0) * sqrt2310 * Y6_p1 * R2 + (135.0 / 143.0) * sqrt11 * Y4_p1 * R4 + (-5.0 / 33.0) * sqrt330 * Y2_p1 * R6) * V[1] + ((5.0 / 44.0) * sqrt2310 * Y6_p1 + (-135.0 / 22.0) * sqrt11 * Y4_p1 * R2 + 1.25 * sqrt330 * Y2_p1 * R4) * V[2] + (11.25 * sqrt11 * Y4_p1 - 3.75 * sqrt330 * Y2_p1 * R2) * V[3] + 3.28125 * sqrt330 * Y2_p1 * V[4];
    d[0] = Y4_n4 * Y6_n6 * V[0] + ((1.0 / 143.0) * sqrt231 * Y8_p2 + (-1.0 / 11.0) * sqrt77 * Y6_p2 * R2 + (45.0 / 143.0) * sqrt66 * Y4_p2 * R4 + (-5.0 / 11.0) * sqrt110 * Y2_p2 * R6) * V[1] + ((15.0 / 44.0) * sqrt77 * Y6_p2 + (-45.0 / 22.0) * sqrt66 * Y4_p2 * R2 + 3.75 * sqrt110 * Y2_p2 * R4) * V[2] + (3.75 * sqrt66 * Y4_p2 - 11.25 * sqrt110 * Y2_p2 * R2) * V[3] + 9.84375 * sqrt110 * Y2_p2 * V[4];
}

}  // namespace ovlab