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

#include "OverlapABRecGG.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_g_g(
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
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt462 = std::sqrt(462.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^4 · β^4 · p^{-8} · (s|s)
    //   V[1] ↔ α^3 · β^3 · p^{-7} · (s|s)
    //   V[2] ↔ α^2 · β^2 · p^{-6} · (s|s)
    //   V[3] ↔ α · β · p^{-5} · (s|s)
    //   V[4] ↔ p^{-4} · (s|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto pinv7 = pinv6 * pinv;
            const auto pinv8 = pinv7 * pinv;

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha4 * beta4 * pinv8;
            V[1] += cab_ss * alpha3 * beta3 * pinv7;
            V[2] += cab_ss * alpha2 * beta2 * pinv6;
            V[3] += cab_ss * alpha * beta * pinv5;
            V[4] += cab_ss * pinv4;
        }
    }

    // ---- Phase 3: fused M·V → 9 × 9 spherical block ----
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
    const auto R4 = R2 * R2;
    const auto R6 = R4 * R2;
    newints::Block out{9, 9, std::vector<double>(81, 0.0)};
    auto *d = out.data.data();
    d[40] = Y4_p0 * Y4_p0 * V[0] + ((-50.0 / 33.0) * Y6_p0 + (-162.0 / 77.0) * Y4_p0 * R2 + (-50.0 / 21.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + ((81.0 / 14.0) * Y4_p0 + (75.0 / 7.0) * Y2_p0 * R2 + 10.5 * R4) * V[2] + (-12.5 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[41] = d[49] = Y4_p0 * Y4_p1 * V[0] + ((-5.0 / 66.0) * sqrt210 * Y6_p1 + (-81.0 / 77.0) * Y4_p1 * R2 + (-5.0 / 42.0) * sqrt30 * Y2_p1 * R4) * V[1] + ((81.0 / 28.0) * Y4_p1 + (15.0 / 28.0) * sqrt30 * Y2_p1 * R2) * V[2] - 0.625 * sqrt30 * Y2_p1 * V[3];
    d[42] = d[58] = Y4_p0 * Y4_p2 * V[0] + ((9.0 / 7.0) * Y4_p2 * R2 + (5.0 / 7.0) * sqrt15 * Y2_p2 * R4) * V[1] + ((-99.0 / 28.0) * Y4_p2 + (-45.0 / 14.0) * sqrt15 * Y2_p2 * R2) * V[2] + 3.75 * sqrt15 * Y2_p2 * V[3];
    d[43] = d[67] = Y4_p0 * Y4_p3 * V[0] + ((25.0 / 33.0) * sqrt3 * Y6_p3 + (27.0 / 11.0) * Y4_p3 * R2) * V[1] - 6.75 * Y4_p3 * V[2];
    d[44] = d[76] = Y4_p0 * Y4_p4 * V[0] + ((10.0 / 11.0) * sqrt5 * Y6_p4 + (-18.0 / 11.0) * Y4_p4 * R2) * V[1] + 4.5 * Y4_p4 * V[2];
    d[39] = d[31] = Y4_p0 * Y4_n1 * V[0] + ((-5.0 / 66.0) * sqrt210 * Y6_n1 + (-81.0 / 77.0) * Y4_n1 * R2 + (-5.0 / 42.0) * sqrt30 * Y2_n1 * R4) * V[1] + ((81.0 / 28.0) * Y4_n1 + (15.0 / 28.0) * sqrt30 * Y2_n1 * R2) * V[2] - 0.625 * sqrt30 * Y2_n1 * V[3];
    d[38] = d[22] = Y4_p0 * Y4_n2 * V[0] + ((9.0 / 7.0) * Y4_n2 * R2 + (5.0 / 7.0) * sqrt15 * Y2_n2 * R4) * V[1] + ((-99.0 / 28.0) * Y4_n2 + (-45.0 / 14.0) * sqrt15 * Y2_n2 * R2) * V[2] + 3.75 * sqrt15 * Y2_n2 * V[3];
    d[37] = d[13] = Y4_p0 * Y4_n3 * V[0] + ((25.0 / 33.0) * sqrt3 * Y6_n3 + (27.0 / 11.0) * Y4_n3 * R2) * V[1] - 6.75 * Y4_n3 * V[2];
    d[36] = d[4] = Y4_p0 * Y4_n4 * V[0] + ((10.0 / 11.0) * sqrt5 * Y6_n4 + (-18.0 / 11.0) * Y4_n4 * R2) * V[1] + 4.5 * Y4_n4 * V[2];
    d[50] = Y4_p1 * Y4_p1 * V[0] + ((5.0 / 66.0) * Y6_p0 + (-5.0 / 66.0) * sqrt210 * Y6_p2 + (-81.0 / 77.0) * Y4_p0 * R2 + (-54.0 / 77.0) * sqrt5 * Y4_p2 * R2 + (-85.0 / 42.0) * Y2_p0 * R4 + (-25.0 / 21.0) * sqrt3 * Y2_p2 * R4 - 2.0 * R6) * V[1] + ((81.0 / 28.0) * Y4_p0 + (27.0 / 14.0) * sqrt5 * Y4_p2 + (255.0 / 28.0) * Y2_p0 * R2 + (75.0 / 14.0) * sqrt3 * Y2_p2 * R2 + 10.5 * R4) * V[2] + (-10.625 * Y2_p0 - 6.25 * sqrt3 * Y2_p2 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[51] = d[59] = Y4_p1 * Y4_p2 * V[0] + ((-5.0 / 44.0) * sqrt42 * Y6_p1 + (-5.0 / 66.0) * sqrt105 * Y6_p3 + (-54.0 / 77.0) * sqrt5 * Y4_p1 * R2 + (-9.0 / 77.0) * sqrt35 * Y4_p3 * R2 + (-15.0 / 28.0) * sqrt6 * Y2_p1 * R4) * V[1] + ((27.0 / 14.0) * sqrt5 * Y4_p1 + (9.0 / 28.0) * sqrt35 * Y4_p3 + (135.0 / 56.0) * sqrt6 * Y2_p1 * R2) * V[2] - 2.8125 * sqrt6 * Y2_p1 * V[3];
    d[52] = d[68] = Y4_p1 * Y4_p3 * V[0] + ((-5.0 / 22.0) * sqrt30 * Y6_p2 + (5.0 / 22.0) * Y6_p4 + (-9.0 / 77.0) * sqrt35 * Y4_p2 * R2 + (9.0 / 11.0) * sqrt5 * Y4_p4 * R2 + (5.0 / 14.0) * sqrt21 * Y2_p2 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_p2 - 2.25 * sqrt5 * Y4_p4 + (-45.0 / 28.0) * sqrt21 * Y2_p2 * R2) * V[2] + 1.875 * sqrt21 * Y2_p2 * V[3];
    d[53] = d[77] = Y4_p1 * Y4_p4 * V[0] + ((-10.0 / 33.0) * sqrt15 * Y6_p3 + (5.0 / 11.0) * sqrt11 * Y6_p5 + (9.0 / 11.0) * sqrt5 * Y4_p3 * R2) * V[1] - 2.25 * sqrt5 * Y4_p3 * V[2];
    d[48] = d[32] = Y4_p1 * Y4_n1 * V[0] + ((-5.0 / 66.0) * sqrt210 * Y6_n2 + (-54.0 / 77.0) * sqrt5 * Y4_n2 * R2 + (-25.0 / 21.0) * sqrt3 * Y2_n2 * R4) * V[1] + ((27.0 / 14.0) * sqrt5 * Y4_n2 + (75.0 / 14.0) * sqrt3 * Y2_n2 * R2) * V[2] - 6.25 * sqrt3 * Y2_n2 * V[3];
    d[47] = d[23] = Y4_p1 * Y4_n2 * V[0] + ((-5.0 / 66.0) * sqrt105 * Y6_n3 + (-5.0 / 44.0) * sqrt42 * Y6_n1 + (-9.0 / 77.0) * sqrt35 * Y4_n3 * R2 + (-54.0 / 77.0) * sqrt5 * Y4_n1 * R2 + (-15.0 / 28.0) * sqrt6 * Y2_n1 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_n3 + (27.0 / 14.0) * sqrt5 * Y4_n1 + (135.0 / 56.0) * sqrt6 * Y2_n1 * R2) * V[2] - 2.8125 * sqrt6 * Y2_n1 * V[3];
    d[46] = d[14] = Y4_p1 * Y4_n3 * V[0] + ((5.0 / 22.0) * Y6_n4 + (-5.0 / 22.0) * sqrt30 * Y6_n2 + (9.0 / 11.0) * sqrt5 * Y4_n4 * R2 + (-9.0 / 77.0) * sqrt35 * Y4_n2 * R2 + (5.0 / 14.0) * sqrt21 * Y2_n2 * R4) * V[1] + (-2.25 * sqrt5 * Y4_n4 + (9.0 / 28.0) * sqrt35 * Y4_n2 + (-45.0 / 28.0) * sqrt21 * Y2_n2 * R2) * V[2] + 1.875 * sqrt21 * Y2_n2 * V[3];
    d[45] = d[5] = Y4_p1 * Y4_n4 * V[0] + ((5.0 / 11.0) * sqrt11 * Y6_n5 + (-10.0 / 33.0) * sqrt15 * Y6_n3 + (9.0 / 11.0) * sqrt5 * Y4_n3 * R2) * V[1] - 2.25 * sqrt5 * Y4_n3 * V[2];
    d[60] = Y4_p2 * Y4_p2 * V[0] + ((5.0 / 3.0) * Y6_p0 + (-5.0 / 11.0) * sqrt7 * Y6_p4 + (9.0 / 7.0) * Y4_p0 * R2 + (-27.0 / 77.0) * sqrt35 * Y4_p4 * R2 + (-20.0 / 21.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + ((-99.0 / 28.0) * Y4_p0 + (27.0 / 28.0) * sqrt35 * Y4_p4 + (30.0 / 7.0) * Y2_p0 * R2 + 10.5 * R4) * V[2] + (-5.0 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[61] = d[69] = Y4_p2 * Y4_p3 * V[0] + ((65.0 / 132.0) * sqrt6 * Y6_p1 + (-5.0 / 22.0) * sqrt11 * Y6_p5 + (-9.0 / 77.0) * sqrt35 * Y4_p1 * R2 + (-25.0 / 84.0) * sqrt42 * Y2_p1 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_p1 + (75.0 / 56.0) * sqrt42 * Y2_p1 * R2) * V[2] - 1.5625 * sqrt42 * Y2_p1 * V[3];
    d[62] = d[78] = Y4_p2 * Y4_p4 * V[0] + ((5.0 / 33.0) * sqrt30 * Y6_p2 + (5.0 / 33.0) * sqrt66 * Y6_p6 + (-27.0 / 77.0) * sqrt35 * Y4_p2 * R2 + (5.0 / 21.0) * sqrt21 * Y2_p2 * R4) * V[1] + ((27.0 / 28.0) * sqrt35 * Y4_p2 + (-15.0 / 14.0) * sqrt21 * Y2_p2 * R2) * V[2] + 1.25 * sqrt21 * Y2_p2 * V[3];
    d[57] = d[33] = Y4_p2 * Y4_n1 * V[0] + ((-5.0 / 66.0) * sqrt105 * Y6_n3 + (5.0 / 44.0) * sqrt42 * Y6_n1 + (-9.0 / 77.0) * sqrt35 * Y4_n3 * R2 + (54.0 / 77.0) * sqrt5 * Y4_n1 * R2 + (15.0 / 28.0) * sqrt6 * Y2_n1 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_n3 + (-27.0 / 14.0) * sqrt5 * Y4_n1 + (-135.0 / 56.0) * sqrt6 * Y2_n1 * R2) * V[2] + 2.8125 * sqrt6 * Y2_n1 * V[3];
    d[56] = d[24] = Y4_p2 * Y4_n2 * V[0] + ((-5.0 / 11.0) * sqrt7 * Y6_n4 + (-27.0 / 77.0) * sqrt35 * Y4_n4 * R2) * V[1] + (27.0 / 28.0) * sqrt35 * Y4_n4 * V[2];
    d[55] = d[15] = Y4_p2 * Y4_n3 * V[0] + ((-5.0 / 22.0) * sqrt11 * Y6_n5 + (65.0 / 132.0) * sqrt6 * Y6_n1 + (-9.0 / 77.0) * sqrt35 * Y4_n1 * R2 + (-25.0 / 84.0) * sqrt42 * Y2_n1 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_n1 + (75.0 / 56.0) * sqrt42 * Y2_n1 * R2) * V[2] - 1.5625 * sqrt42 * Y2_n1 * V[3];
    d[54] = d[6] = Y4_p2 * Y4_n4 * V[0] + ((5.0 / 33.0) * sqrt66 * Y6_n6 + (5.0 / 33.0) * sqrt30 * Y6_n2 + (-27.0 / 77.0) * sqrt35 * Y4_n2 * R2 + (5.0 / 21.0) * sqrt21 * Y2_n2 * R4) * V[1] + ((27.0 / 28.0) * sqrt35 * Y4_n2 + (-15.0 / 14.0) * sqrt21 * Y2_n2 * R2) * V[2] + 1.25 * sqrt21 * Y2_n2 * V[3];
    d[70] = Y4_p3 * Y4_p3 * V[0] + ((-85.0 / 66.0) * Y6_p0 + (-5.0 / 66.0) * sqrt462 * Y6_p6 + (27.0 / 11.0) * Y4_p0 * R2 + (5.0 / 6.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + (-6.75 * Y4_p0 - 3.75 * Y2_p0 * R2 + 10.5 * R4) * V[2] + (4.375 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[71] = d[79] = Y4_p3 * Y4_p4 * V[0] + ((-5.0 / 66.0) * sqrt42 * Y6_p1 + (9.0 / 11.0) * sqrt5 * Y4_p1 * R2 + (-5.0 / 6.0) * sqrt6 * Y2_p1 * R4) * V[1] + (-2.25 * sqrt5 * Y4_p1 + 3.75 * sqrt6 * Y2_p1 * R2) * V[2] - 4.375 * sqrt6 * Y2_p1 * V[3];
    d[66] = d[34] = Y4_p3 * Y4_n1 * V[0] + ((5.0 / 22.0) * Y6_n4 + (5.0 / 22.0) * sqrt30 * Y6_n2 + (9.0 / 11.0) * sqrt5 * Y4_n4 * R2 + (9.0 / 77.0) * sqrt35 * Y4_n2 * R2 + (-5.0 / 14.0) * sqrt21 * Y2_n2 * R4) * V[1] + (-2.25 * sqrt5 * Y4_n4 + (-9.0 / 28.0) * sqrt35 * Y4_n2 + (45.0 / 28.0) * sqrt21 * Y2_n2 * R2) * V[2] - 1.875 * sqrt21 * Y2_n2 * V[3];
    d[65] = d[25] = Y4_p3 * Y4_n2 * V[0] + ((-5.0 / 22.0) * sqrt11 * Y6_n5 + (-65.0 / 132.0) * sqrt6 * Y6_n1 + (9.0 / 77.0) * sqrt35 * Y4_n1 * R2 + (25.0 / 84.0) * sqrt42 * Y2_n1 * R4) * V[1] + ((-9.0 / 28.0) * sqrt35 * Y4_n1 + (-75.0 / 56.0) * sqrt42 * Y2_n1 * R2) * V[2] + 1.5625 * sqrt42 * Y2_n1 * V[3];
    d[64] = d[16] = Y4_p3 * Y4_n3 * V[0] + (-5.0 / 66.0) * sqrt462 * Y6_n6 * V[1];
    d[63] = d[7] = Y4_p3 * Y4_n4 * V[0] + ((-5.0 / 66.0) * sqrt42 * Y6_n1 + (9.0 / 11.0) * sqrt5 * Y4_n1 * R2 + (-5.0 / 6.0) * sqrt6 * Y2_n1 * R4) * V[1] + (-2.25 * sqrt5 * Y4_n1 + 3.75 * sqrt6 * Y2_n1 * R2) * V[2] - 4.375 * sqrt6 * Y2_n1 * V[3];
    d[80] = Y4_p4 * Y4_p4 * V[0] + ((10.0 / 33.0) * Y6_p0 + (-18.0 / 11.0) * Y4_p0 * R2 + (10.0 / 3.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + (4.5 * Y4_p0 - 15.0 * Y2_p0 * R2 + 10.5 * R4) * V[2] + (17.5 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[75] = d[35] = Y4_p4 * Y4_n1 * V[0] + ((5.0 / 11.0) * sqrt11 * Y6_n5 + (10.0 / 33.0) * sqrt15 * Y6_n3 + (-9.0 / 11.0) * sqrt5 * Y4_n3 * R2) * V[1] + 2.25 * sqrt5 * Y4_n3 * V[2];
    d[74] = d[26] = Y4_p4 * Y4_n2 * V[0] + ((5.0 / 33.0) * sqrt66 * Y6_n6 + (-5.0 / 33.0) * sqrt30 * Y6_n2 + (27.0 / 77.0) * sqrt35 * Y4_n2 * R2 + (-5.0 / 21.0) * sqrt21 * Y2_n2 * R4) * V[1] + ((-27.0 / 28.0) * sqrt35 * Y4_n2 + (15.0 / 14.0) * sqrt21 * Y2_n2 * R2) * V[2] - 1.25 * sqrt21 * Y2_n2 * V[3];
    d[73] = d[17] = Y4_p4 * Y4_n3 * V[0] + ((5.0 / 66.0) * sqrt42 * Y6_n1 + (-9.0 / 11.0) * sqrt5 * Y4_n1 * R2 + (5.0 / 6.0) * sqrt6 * Y2_n1 * R4) * V[1] + (2.25 * sqrt5 * Y4_n1 - 3.75 * sqrt6 * Y2_n1 * R2) * V[2] + 4.375 * sqrt6 * Y2_n1 * V[3];
    d[72] = d[8] = Y4_p4 * Y4_n4 * V[0];
    d[30] = Y4_n1 * Y4_n1 * V[0] + ((5.0 / 66.0) * Y6_p0 + (5.0 / 66.0) * sqrt210 * Y6_p2 + (-81.0 / 77.0) * Y4_p0 * R2 + (54.0 / 77.0) * sqrt5 * Y4_p2 * R2 + (-85.0 / 42.0) * Y2_p0 * R4 + (25.0 / 21.0) * sqrt3 * Y2_p2 * R4 - 2.0 * R6) * V[1] + ((81.0 / 28.0) * Y4_p0 + (-27.0 / 14.0) * sqrt5 * Y4_p2 + (255.0 / 28.0) * Y2_p0 * R2 + (-75.0 / 14.0) * sqrt3 * Y2_p2 * R2 + 10.5 * R4) * V[2] + (-10.625 * Y2_p0 + 6.25 * sqrt3 * Y2_p2 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[29] = d[21] = Y4_n1 * Y4_n2 * V[0] + ((-5.0 / 44.0) * sqrt42 * Y6_p1 + (5.0 / 66.0) * sqrt105 * Y6_p3 + (-54.0 / 77.0) * sqrt5 * Y4_p1 * R2 + (9.0 / 77.0) * sqrt35 * Y4_p3 * R2 + (-15.0 / 28.0) * sqrt6 * Y2_p1 * R4) * V[1] + ((27.0 / 14.0) * sqrt5 * Y4_p1 + (-9.0 / 28.0) * sqrt35 * Y4_p3 + (135.0 / 56.0) * sqrt6 * Y2_p1 * R2) * V[2] - 2.8125 * sqrt6 * Y2_p1 * V[3];
    d[28] = d[12] = Y4_n1 * Y4_n3 * V[0] + ((-5.0 / 22.0) * sqrt30 * Y6_p2 + (-5.0 / 22.0) * Y6_p4 + (-9.0 / 77.0) * sqrt35 * Y4_p2 * R2 + (-9.0 / 11.0) * sqrt5 * Y4_p4 * R2 + (5.0 / 14.0) * sqrt21 * Y2_p2 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_p2 + 2.25 * sqrt5 * Y4_p4 + (-45.0 / 28.0) * sqrt21 * Y2_p2 * R2) * V[2] + 1.875 * sqrt21 * Y2_p2 * V[3];
    d[27] = d[3] = Y4_n1 * Y4_n4 * V[0] + ((-10.0 / 33.0) * sqrt15 * Y6_p3 + (-5.0 / 11.0) * sqrt11 * Y6_p5 + (9.0 / 11.0) * sqrt5 * Y4_p3 * R2) * V[1] - 2.25 * sqrt5 * Y4_p3 * V[2];
    d[20] = Y4_n2 * Y4_n2 * V[0] + ((5.0 / 3.0) * Y6_p0 + (5.0 / 11.0) * sqrt7 * Y6_p4 + (9.0 / 7.0) * Y4_p0 * R2 + (27.0 / 77.0) * sqrt35 * Y4_p4 * R2 + (-20.0 / 21.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + ((-99.0 / 28.0) * Y4_p0 + (-27.0 / 28.0) * sqrt35 * Y4_p4 + (30.0 / 7.0) * Y2_p0 * R2 + 10.5 * R4) * V[2] + (-5.0 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[19] = d[11] = Y4_n2 * Y4_n3 * V[0] + ((65.0 / 132.0) * sqrt6 * Y6_p1 + (5.0 / 22.0) * sqrt11 * Y6_p5 + (-9.0 / 77.0) * sqrt35 * Y4_p1 * R2 + (-25.0 / 84.0) * sqrt42 * Y2_p1 * R4) * V[1] + ((9.0 / 28.0) * sqrt35 * Y4_p1 + (75.0 / 56.0) * sqrt42 * Y2_p1 * R2) * V[2] - 1.5625 * sqrt42 * Y2_p1 * V[3];
    d[18] = d[2] = Y4_n2 * Y4_n4 * V[0] + ((5.0 / 33.0) * sqrt30 * Y6_p2 + (-5.0 / 33.0) * sqrt66 * Y6_p6 + (-27.0 / 77.0) * sqrt35 * Y4_p2 * R2 + (5.0 / 21.0) * sqrt21 * Y2_p2 * R4) * V[1] + ((27.0 / 28.0) * sqrt35 * Y4_p2 + (-15.0 / 14.0) * sqrt21 * Y2_p2 * R2) * V[2] + 1.25 * sqrt21 * Y2_p2 * V[3];
    d[10] = Y4_n3 * Y4_n3 * V[0] + ((-85.0 / 66.0) * Y6_p0 + (5.0 / 66.0) * sqrt462 * Y6_p6 + (27.0 / 11.0) * Y4_p0 * R2 + (5.0 / 6.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + (-6.75 * Y4_p0 - 3.75 * Y2_p0 * R2 + 10.5 * R4) * V[2] + (4.375 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];
    d[9] = d[1] = Y4_n3 * Y4_n4 * V[0] + ((-5.0 / 66.0) * sqrt42 * Y6_p1 + (9.0 / 11.0) * sqrt5 * Y4_p1 * R2 + (-5.0 / 6.0) * sqrt6 * Y2_p1 * R4) * V[1] + (-2.25 * sqrt5 * Y4_p1 + 3.75 * sqrt6 * Y2_p1 * R2) * V[2] - 4.375 * sqrt6 * Y2_p1 * V[3];
    d[0] = Y4_n4 * Y4_n4 * V[0] + ((10.0 / 33.0) * Y6_p0 + (-18.0 / 11.0) * Y4_p0 * R2 + (10.0 / 3.0) * Y2_p0 * R4 - 2.0 * R6) * V[1] + (4.5 * Y4_p0 - 15.0 * Y2_p0 * R2 + 10.5 * R4) * V[2] + (17.5 * Y2_p0 - 17.5 * R2) * V[3] + 6.5625 * V[4];

    return out;
}

}  // namespace ovlab