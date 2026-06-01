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

#include "KineticEnergyABRecID.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_i_d(
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
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt22 = std::sqrt(22.0);
    const auto sqrt26 = std::sqrt(26.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt55 = std::sqrt(55.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt91 = std::sqrt(91.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt330 = std::sqrt(330.0);
    const auto sqrt429 = std::sqrt(429.0);
    const auto sqrt462 = std::sqrt(462.0);
    const auto sqrt910 = std::sqrt(910.0);
    const auto sqrt1430 = std::sqrt(1430.0);
    const auto sqrt2002 = std::sqrt(2002.0);
    const auto sqrt2730 = std::sqrt(2730.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^2 · β^6 · p^{-8} · (s|T|s)
    //   V[1] ↔ α · β^5 · p^{-7} · (s|T|s)
    //   V[2] ↔ β^4 · p^{-6} · (s|T|s)
    //   V[3] ↔ α^2 · β^6 · p^{-8} · ξ · (s|S|s)
    //   V[4] ↔ α · β^5 · p^{-7} · ξ · (s|S|s)
    //   V[5] ↔ β^4 · p^{-6} · ξ · (s|S|s)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 6> V = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < exps_a.size(); ++i)
    {
        const auto alpha = exps_a[i];
        const auto ca    = coefs_a[i];
        const auto alpha2 = alpha * alpha;

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
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha2 * beta6 * pinv8;
            V[1] += cab_tt * alpha * beta5 * pinv7;
            V[2] += cab_tt * beta4 * pinv6;
            V[3] += cab_ss * alpha2 * beta6 * pinv8 * xi;
            V[4] += cab_ss * alpha * beta5 * pinv7 * xi;
            V[5] += cab_ss * beta4 * pinv6 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 13 × 5 spherical block ----
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
    auto *d = buffer;
    d[32] = ((28.0 / 65.0) * Y8_p0 + (14.0 / 55.0) * Y6_p0 * R2 + (45.0 / 143.0) * Y4_p0 * R4) * V[0] + ((-21.0 / 11.0) * Y6_p0 + (-45.0 / 11.0) * Y4_p0 * R2) * V[1] + 11.25 * Y4_p0 * V[2] + ((448.0 / 65.0) * Y8_p0 + (224.0 / 55.0) * Y6_p0 * R2 + (720.0 / 143.0) * Y4_p0 * R4) * V[3] + ((-294.0 / 11.0) * Y6_p0 + (-630.0 / 11.0) * Y4_p0 * R2) * V[4] + 135.0 * Y4_p0 * V[5];
    d[33] = ((14.0 / 65.0) * sqrt3 * Y8_p1 + (1.0 / 55.0) * sqrt7 * Y6_p1 * R2 + (-6.0 / 143.0) * sqrt30 * Y4_p1 * R4) * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_p1 + (6.0 / 11.0) * sqrt30 * Y4_p1 * R2) * V[1] - 1.5 * sqrt30 * Y4_p1 * V[2] + ((224.0 / 65.0) * sqrt3 * Y8_p1 + (16.0 / 55.0) * sqrt7 * Y6_p1 * R2 + (-96.0 / 143.0) * sqrt30 * Y4_p1 * R4) * V[3] + ((-21.0 / 11.0) * sqrt7 * Y6_p1 + (84.0 / 11.0) * sqrt30 * Y4_p1 * R2) * V[4] - 18.0 * sqrt30 * Y4_p1 * V[5];
    d[34] = ((1.0 / 65.0) * sqrt210 * Y8_p2 + (-2.0 / 55.0) * sqrt70 * Y6_p2 * R2 + (3.0 / 143.0) * sqrt15 * Y4_p2 * R4) * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_p2 + (-3.0 / 11.0) * sqrt15 * Y4_p2 * R2) * V[1] + 0.75 * sqrt15 * Y4_p2 * V[2] + ((16.0 / 65.0) * sqrt210 * Y8_p2 + (-32.0 / 55.0) * sqrt70 * Y6_p2 * R2 + (48.0 / 143.0) * sqrt15 * Y4_p2 * R4) * V[3] + ((42.0 / 11.0) * sqrt70 * Y6_p2 + (-42.0 / 11.0) * sqrt15 * Y4_p2 * R2) * V[4] + 9.0 * sqrt15 * Y4_p2 * V[5];
    d[31] = ((14.0 / 65.0) * sqrt3 * Y8_n1 + (1.0 / 55.0) * sqrt7 * Y6_n1 * R2 + (-6.0 / 143.0) * sqrt30 * Y4_n1 * R4) * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_n1 + (6.0 / 11.0) * sqrt30 * Y4_n1 * R2) * V[1] - 1.5 * sqrt30 * Y4_n1 * V[2] + ((224.0 / 65.0) * sqrt3 * Y8_n1 + (16.0 / 55.0) * sqrt7 * Y6_n1 * R2 + (-96.0 / 143.0) * sqrt30 * Y4_n1 * R4) * V[3] + ((-21.0 / 11.0) * sqrt7 * Y6_n1 + (84.0 / 11.0) * sqrt30 * Y4_n1 * R2) * V[4] - 18.0 * sqrt30 * Y4_n1 * V[5];
    d[30] = ((1.0 / 65.0) * sqrt210 * Y8_n2 + (-2.0 / 55.0) * sqrt70 * Y6_n2 * R2 + (3.0 / 143.0) * sqrt15 * Y4_n2 * R4) * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_n2 + (-3.0 / 11.0) * sqrt15 * Y4_n2 * R2) * V[1] + 0.75 * sqrt15 * Y4_n2 * V[2] + ((16.0 / 65.0) * sqrt210 * Y8_n2 + (-32.0 / 55.0) * sqrt70 * Y6_n2 * R2 + (48.0 / 143.0) * sqrt15 * Y4_n2 * R4) * V[3] + ((42.0 / 11.0) * sqrt70 * Y6_n2 + (-42.0 / 11.0) * sqrt15 * Y4_n2 * R2) * V[4] + 9.0 * sqrt15 * Y4_n2 * V[5];
    d[37] = ((6.0 / 65.0) * sqrt21 * Y8_p1 + (13.0 / 55.0) * Y6_p1 * R2 + (3.0 / 143.0) * sqrt210 * Y4_p1 * R4) * V[0] + ((-39.0 / 22.0) * Y6_p1 + (-3.0 / 11.0) * sqrt210 * Y4_p1 * R2) * V[1] + 0.75 * sqrt210 * Y4_p1 * V[2] + ((96.0 / 65.0) * sqrt21 * Y8_p1 + (208.0 / 55.0) * Y6_p1 * R2 + (48.0 / 143.0) * sqrt210 * Y4_p1 * R4) * V[3] + ((-273.0 / 11.0) * Y6_p1 + (-42.0 / 11.0) * sqrt210 * Y4_p1 * R2) * V[4] + 9.0 * sqrt210 * Y4_p1 * V[5];
    d[38] = ((-8.0 / 65.0) * sqrt7 * Y8_p0 + (6.0 / 65.0) * sqrt10 * Y8_p2 + (1.0 / 55.0) * sqrt7 * Y6_p0 * R2 + (1.0 / 55.0) * sqrt30 * Y6_p2 * R2 + (15.0 / 143.0) * sqrt7 * Y4_p0 * R4 + (-3.0 / 143.0) * sqrt35 * Y4_p2 * R4) * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_p0 + (-3.0 / 22.0) * sqrt30 * Y6_p2 + (-15.0 / 11.0) * sqrt7 * Y4_p0 * R2 + (3.0 / 11.0) * sqrt35 * Y4_p2 * R2) * V[1] + (3.75 * sqrt7 * Y4_p0 - 0.75 * sqrt35 * Y4_p2) * V[2] + ((-128.0 / 65.0) * sqrt7 * Y8_p0 + (96.0 / 65.0) * sqrt10 * Y8_p2 + (16.0 / 55.0) * sqrt7 * Y6_p0 * R2 + (16.0 / 55.0) * sqrt30 * Y6_p2 * R2 + (240.0 / 143.0) * sqrt7 * Y4_p0 * R4 + (-48.0 / 143.0) * sqrt35 * Y4_p2 * R4) * V[3] + ((-21.0 / 11.0) * sqrt7 * Y6_p0 + (-21.0 / 11.0) * sqrt30 * Y6_p2 + (-210.0 / 11.0) * sqrt7 * Y4_p0 * R2 + (42.0 / 11.0) * sqrt35 * Y4_p2 * R2) * V[4] + (45.0 * sqrt7 * Y4_p0 - 9.0 * sqrt35 * Y4_p2) * V[5];
    d[39] = ((-3.0 / 65.0) * sqrt7 * Y8_p1 + (1.0 / 65.0) * sqrt165 * Y8_p3 + (7.0 / 55.0) * sqrt3 * Y6_p1 * R2 + (-2.0 / 55.0) * sqrt30 * Y6_p3 * R2 + (-3.0 / 286.0) * sqrt70 * Y4_p1 * R4 + (3.0 / 286.0) * sqrt10 * Y4_p3 * R4) * V[0] + ((-21.0 / 22.0) * sqrt3 * Y6_p1 + (3.0 / 11.0) * sqrt30 * Y6_p3 + (3.0 / 22.0) * sqrt70 * Y4_p1 * R2 + (-3.0 / 22.0) * sqrt10 * Y4_p3 * R2) * V[1] + (-0.375 * sqrt70 * Y4_p1 + 0.375 * sqrt10 * Y4_p3) * V[2] + ((-48.0 / 65.0) * sqrt7 * Y8_p1 + (16.0 / 65.0) * sqrt165 * Y8_p3 + (112.0 / 55.0) * sqrt3 * Y6_p1 * R2 + (-32.0 / 55.0) * sqrt30 * Y6_p3 * R2 + (-24.0 / 143.0) * sqrt70 * Y4_p1 * R4 + (24.0 / 143.0) * sqrt10 * Y4_p3 * R4) * V[3] + ((-147.0 / 11.0) * sqrt3 * Y6_p1 + (42.0 / 11.0) * sqrt30 * Y6_p3 + (21.0 / 11.0) * sqrt70 * Y4_p1 * R2 + (-21.0 / 11.0) * sqrt10 * Y4_p3 * R2) * V[4] + (-4.5 * sqrt70 * Y4_p1 + 4.5 * sqrt10 * Y4_p3) * V[5];
    d[36] = ((6.0 / 65.0) * sqrt10 * Y8_n2 + (1.0 / 55.0) * sqrt30 * Y6_n2 * R2 + (-3.0 / 143.0) * sqrt35 * Y4_n2 * R4) * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_n2 + (3.0 / 11.0) * sqrt35 * Y4_n2 * R2) * V[1] - 0.75 * sqrt35 * Y4_n2 * V[2] + ((96.0 / 65.0) * sqrt10 * Y8_n2 + (16.0 / 55.0) * sqrt30 * Y6_n2 * R2 + (-48.0 / 143.0) * sqrt35 * Y4_n2 * R4) * V[3] + ((-21.0 / 11.0) * sqrt30 * Y6_n2 + (42.0 / 11.0) * sqrt35 * Y4_n2 * R2) * V[4] - 9.0 * sqrt35 * Y4_n2 * V[5];
    d[35] = ((1.0 / 65.0) * sqrt165 * Y8_n3 + (-3.0 / 65.0) * sqrt7 * Y8_n1 + (-2.0 / 55.0) * sqrt30 * Y6_n3 * R2 + (7.0 / 55.0) * sqrt3 * Y6_n1 * R2 + (3.0 / 286.0) * sqrt10 * Y4_n3 * R4 + (-3.0 / 286.0) * sqrt70 * Y4_n1 * R4) * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_n3 + (-21.0 / 22.0) * sqrt3 * Y6_n1 + (-3.0 / 22.0) * sqrt10 * Y4_n3 * R2 + (3.0 / 22.0) * sqrt70 * Y4_n1 * R2) * V[1] + (0.375 * sqrt10 * Y4_n3 - 0.375 * sqrt70 * Y4_n1) * V[2] + ((16.0 / 65.0) * sqrt165 * Y8_n3 + (-48.0 / 65.0) * sqrt7 * Y8_n1 + (-32.0 / 55.0) * sqrt30 * Y6_n3 * R2 + (112.0 / 55.0) * sqrt3 * Y6_n1 * R2 + (24.0 / 143.0) * sqrt10 * Y4_n3 * R4 + (-24.0 / 143.0) * sqrt70 * Y4_n1 * R4) * V[3] + ((42.0 / 11.0) * sqrt30 * Y6_n3 + (-147.0 / 11.0) * sqrt3 * Y6_n1 + (-21.0 / 11.0) * sqrt10 * Y4_n3 * R2 + (21.0 / 11.0) * sqrt70 * Y4_n1 * R2) * V[4] + (4.5 * sqrt10 * Y4_n3 - 4.5 * sqrt70 * Y4_n1) * V[5];
    d[42] = ((3.0 / 13.0) * sqrt3 * Y8_p2 + (2.0 / 11.0) * Y6_p2 * R2 + (6.0 / 143.0) * sqrt42 * Y4_p2 * R4) * V[0] + ((-15.0 / 11.0) * Y6_p2 + (-6.0 / 11.0) * sqrt42 * Y4_p2 * R2) * V[1] + 1.5 * sqrt42 * Y4_p2 * V[2] + ((48.0 / 13.0) * sqrt3 * Y8_p2 + (32.0 / 11.0) * Y6_p2 * R2 + (96.0 / 143.0) * sqrt42 * Y4_p2 * R4) * V[3] + ((-210.0 / 11.0) * Y6_p2 + (-84.0 / 11.0) * sqrt42 * Y4_p2 * R2) * V[4] + 18.0 * sqrt42 * Y4_p2 * V[5];
    d[43] = ((-3.0 / 130.0) * sqrt70 * Y8_p1 + (1.0 / 26.0) * sqrt66 * Y8_p3 + (1.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (1.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (12.0 / 143.0) * sqrt7 * Y4_p1 * R4 + (-12.0 / 143.0) * Y4_p3 * R4) * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_p1 + (-15.0 / 22.0) * sqrt3 * Y6_p3 + (-12.0 / 11.0) * sqrt7 * Y4_p1 * R2 + (12.0 / 11.0) * Y4_p3 * R2) * V[1] + (3.0 * sqrt7 * Y4_p1 - 3.0 * Y4_p3) * V[2] + ((-24.0 / 65.0) * sqrt70 * Y8_p1 + (8.0 / 13.0) * sqrt66 * Y8_p3 + (16.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (16.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (192.0 / 143.0) * sqrt7 * Y4_p1 * R4 + (-192.0 / 143.0) * Y4_p3 * R4) * V[3] + ((-21.0 / 11.0) * sqrt30 * Y6_p1 + (-105.0 / 11.0) * sqrt3 * Y6_p3 + (-168.0 / 11.0) * sqrt7 * Y4_p1 * R2 + (168.0 / 11.0) * Y4_p3 * R2) * V[4] + (36.0 * sqrt7 * Y4_p1 - 36.0 * Y4_p3) * V[5];
    d[44] = ((1.0 / 65.0) * sqrt70 * Y8_p0 + (3.0 / 130.0) * sqrt110 * Y8_p4 + (-2.0 / 55.0) * sqrt70 * Y6_p0 * R2 + (-3.0 / 55.0) * sqrt10 * Y6_p4 * R2 + (3.0 / 143.0) * sqrt70 * Y4_p0 * R4 + (3.0 / 286.0) * sqrt2 * Y4_p4 * R4) * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_p0 + (9.0 / 22.0) * sqrt10 * Y6_p4 + (-3.0 / 11.0) * sqrt70 * Y4_p0 * R2 + (-3.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1] + (0.75 * sqrt70 * Y4_p0 + 0.375 * sqrt2 * Y4_p4) * V[2] + ((16.0 / 65.0) * sqrt70 * Y8_p0 + (24.0 / 65.0) * sqrt110 * Y8_p4 + (-32.0 / 55.0) * sqrt70 * Y6_p0 * R2 + (-48.0 / 55.0) * sqrt10 * Y6_p4 * R2 + (48.0 / 143.0) * sqrt70 * Y4_p0 * R4 + (24.0 / 143.0) * sqrt2 * Y4_p4 * R4) * V[3] + ((42.0 / 11.0) * sqrt70 * Y6_p0 + (63.0 / 11.0) * sqrt10 * Y6_p4 + (-42.0 / 11.0) * sqrt70 * Y4_p0 * R2 + (-21.0 / 11.0) * sqrt2 * Y4_p4 * R2) * V[4] + (9.0 * sqrt70 * Y4_p0 + 4.5 * sqrt2 * Y4_p4) * V[5];
    d[41] = ((1.0 / 26.0) * sqrt66 * Y8_n3 + (3.0 / 130.0) * sqrt70 * Y8_n1 + (1.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (-1.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (-12.0 / 143.0) * Y4_n3 * R4 + (-12.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_n3 + (3.0 / 22.0) * sqrt30 * Y6_n1 + (12.0 / 11.0) * Y4_n3 * R2 + (12.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1] + (-3.0 * Y4_n3 - 3.0 * sqrt7 * Y4_n1) * V[2] + ((8.0 / 13.0) * sqrt66 * Y8_n3 + (24.0 / 65.0) * sqrt70 * Y8_n1 + (16.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (-16.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (-192.0 / 143.0) * Y4_n3 * R4 + (-192.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[3] + ((-105.0 / 11.0) * sqrt3 * Y6_n3 + (21.0 / 11.0) * sqrt30 * Y6_n1 + (168.0 / 11.0) * Y4_n3 * R2 + (168.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[4] + (-36.0 * Y4_n3 - 36.0 * sqrt7 * Y4_n1) * V[5];
    d[40] = ((3.0 / 130.0) * sqrt110 * Y8_n4 + (-3.0 / 55.0) * sqrt10 * Y6_n4 * R2 + (3.0 / 286.0) * sqrt2 * Y4_n4 * R4) * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_n4 + (-3.0 / 22.0) * sqrt2 * Y4_n4 * R2) * V[1] + 0.375 * sqrt2 * Y4_n4 * V[2] + ((24.0 / 65.0) * sqrt110 * Y8_n4 + (-48.0 / 55.0) * sqrt10 * Y6_n4 * R2 + (24.0 / 143.0) * sqrt2 * Y4_n4 * R4) * V[3] + ((63.0 / 11.0) * sqrt10 * Y6_n4 + (-21.0 / 11.0) * sqrt2 * Y4_n4 * R2) * V[4] + 4.5 * sqrt2 * Y4_n4 * V[5];
    d[47] = ((1.0 / 13.0) * sqrt22 * Y8_p3 + (1.0 / 11.0) * Y6_p3 * R2 + (18.0 / 143.0) * sqrt3 * Y4_p3 * R4) * V[0] + ((-15.0 / 22.0) * Y6_p3 + (-18.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[1] + 4.5 * sqrt3 * Y4_p3 * V[2] + ((16.0 / 13.0) * sqrt22 * Y8_p3 + (16.0 / 11.0) * Y6_p3 * R2 + (288.0 / 143.0) * sqrt3 * Y4_p3 * R4) * V[3] + ((-105.0 / 11.0) * Y6_p3 + (-252.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[4] + 54.0 * sqrt3 * Y4_p3 * V[5];
    d[48] = ((-2.0 / 13.0) * Y8_p2 + (2.0 / 65.0) * sqrt110 * Y8_p4 + (1.0 / 11.0) * sqrt3 * Y6_p2 * R2 + (7.0 / 110.0) * sqrt10 * Y6_p4 * R2 + (9.0 / 143.0) * sqrt14 * Y4_p2 * R4 + (-9.0 / 286.0) * sqrt2 * Y4_p4 * R4) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_p2 + (-21.0 / 44.0) * sqrt10 * Y6_p4 + (-9.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (9.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1] + (2.25 * sqrt14 * Y4_p2 - 1.125 * sqrt2 * Y4_p4) * V[2] + ((-32.0 / 13.0) * Y8_p2 + (32.0 / 65.0) * sqrt110 * Y8_p4 + (16.0 / 11.0) * sqrt3 * Y6_p2 * R2 + (56.0 / 55.0) * sqrt10 * Y6_p4 * R2 + (144.0 / 143.0) * sqrt14 * Y4_p2 * R4 + (-72.0 / 143.0) * sqrt2 * Y4_p4 * R4) * V[3] + ((-105.0 / 11.0) * sqrt3 * Y6_p2 + (-147.0 / 22.0) * sqrt10 * Y6_p4 + (-126.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (63.0 / 11.0) * sqrt2 * Y4_p4 * R2) * V[4] + (27.0 * sqrt14 * Y4_p2 - 13.5 * sqrt2 * Y4_p4) * V[5];
    d[49] = ((1.0 / 130.0) * sqrt70 * Y8_p1 + (1.0 / 130.0) * sqrt1430 * Y8_p5 + (-2.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (-1.0 / 55.0) * sqrt55 * Y6_p5 * R2 + (9.0 / 143.0) * sqrt7 * Y4_p1 * R4) * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_p1 + (3.0 / 22.0) * sqrt55 * Y6_p5 + (-9.0 / 11.0) * sqrt7 * Y4_p1 * R2) * V[1] + 2.25 * sqrt7 * Y4_p1 * V[2] + ((8.0 / 65.0) * sqrt70 * Y8_p1 + (8.0 / 65.0) * sqrt1430 * Y8_p5 + (-32.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (-16.0 / 55.0) * sqrt55 * Y6_p5 * R2 + (144.0 / 143.0) * sqrt7 * Y4_p1 * R4) * V[3] + ((42.0 / 11.0) * sqrt30 * Y6_p1 + (21.0 / 11.0) * sqrt55 * Y6_p5 + (-126.0 / 11.0) * sqrt7 * Y4_p1 * R2) * V[4] + 27.0 * sqrt7 * Y4_p1 * V[5];
    d[46] = ((2.0 / 65.0) * sqrt110 * Y8_n4 + (2.0 / 13.0) * Y8_n2 + (7.0 / 110.0) * sqrt10 * Y6_n4 * R2 + (-1.0 / 11.0) * sqrt3 * Y6_n2 * R2 + (-9.0 / 286.0) * sqrt2 * Y4_n4 * R4 + (-9.0 / 143.0) * sqrt14 * Y4_n2 * R4) * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_n4 + (15.0 / 22.0) * sqrt3 * Y6_n2 + (9.0 / 22.0) * sqrt2 * Y4_n4 * R2 + (9.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[1] + (-1.125 * sqrt2 * Y4_n4 - 2.25 * sqrt14 * Y4_n2) * V[2] + ((32.0 / 65.0) * sqrt110 * Y8_n4 + (32.0 / 13.0) * Y8_n2 + (56.0 / 55.0) * sqrt10 * Y6_n4 * R2 + (-16.0 / 11.0) * sqrt3 * Y6_n2 * R2 + (-72.0 / 143.0) * sqrt2 * Y4_n4 * R4 + (-144.0 / 143.0) * sqrt14 * Y4_n2 * R4) * V[3] + ((-147.0 / 22.0) * sqrt10 * Y6_n4 + (105.0 / 11.0) * sqrt3 * Y6_n2 + (63.0 / 11.0) * sqrt2 * Y4_n4 * R2 + (126.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[4] + (-13.5 * sqrt2 * Y4_n4 - 27.0 * sqrt14 * Y4_n2) * V[5];
    d[45] = ((1.0 / 130.0) * sqrt1430 * Y8_n5 + (-1.0 / 130.0) * sqrt70 * Y8_n1 + (-1.0 / 55.0) * sqrt55 * Y6_n5 * R2 + (2.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (-9.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_n5 + (-3.0 / 11.0) * sqrt30 * Y6_n1 + (9.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1] - 2.25 * sqrt7 * Y4_n1 * V[2] + ((8.0 / 65.0) * sqrt1430 * Y8_n5 + (-8.0 / 65.0) * sqrt70 * Y8_n1 + (-16.0 / 55.0) * sqrt55 * Y6_n5 * R2 + (32.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (-144.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[3] + ((21.0 / 11.0) * sqrt55 * Y6_n5 + (-42.0 / 11.0) * sqrt30 * Y6_n1 + (126.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[4] - 27.0 * sqrt7 * Y4_n1 * V[5];
    d[52] = ((6.0 / 65.0) * sqrt11 * Y8_p4 + (-2.0 / 55.0) * Y6_p4 * R2 + (9.0 / 143.0) * sqrt5 * Y4_p4 * R4) * V[0] + ((3.0 / 11.0) * Y6_p4 + (-9.0 / 11.0) * sqrt5 * Y4_p4 * R2) * V[1] + 2.25 * sqrt5 * Y4_p4 * V[2] + ((96.0 / 65.0) * sqrt11 * Y8_p4 + (-32.0 / 55.0) * Y6_p4 * R2 + (144.0 / 143.0) * sqrt5 * Y4_p4 * R4) * V[3] + ((42.0 / 11.0) * Y6_p4 + (-126.0 / 11.0) * sqrt5 * Y4_p4 * R2) * V[4] + 27.0 * sqrt5 * Y4_p4 * V[5];
    d[53] = ((-1.0 / 65.0) * sqrt55 * Y8_p3 + (1.0 / 65.0) * sqrt429 * Y8_p5 + (7.0 / 110.0) * sqrt10 * Y6_p3 * R2 + (3.0 / 110.0) * sqrt66 * Y6_p5 * R2 + (6.0 / 143.0) * sqrt30 * Y4_p3 * R4) * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_p3 + (-9.0 / 44.0) * sqrt66 * Y6_p5 + (-6.0 / 11.0) * sqrt30 * Y4_p3 * R2) * V[1] + 1.5 * sqrt30 * Y4_p3 * V[2] + ((-16.0 / 65.0) * sqrt55 * Y8_p3 + (16.0 / 65.0) * sqrt429 * Y8_p5 + (56.0 / 55.0) * sqrt10 * Y6_p3 * R2 + (24.0 / 55.0) * sqrt66 * Y6_p5 * R2 + (96.0 / 143.0) * sqrt30 * Y4_p3 * R4) * V[3] + ((-147.0 / 22.0) * sqrt10 * Y6_p3 + (-63.0 / 22.0) * sqrt66 * Y6_p5 + (-84.0 / 11.0) * sqrt30 * Y4_p3 * R2) * V[4] + 18.0 * sqrt30 * Y4_p3 * V[5];
    d[54] = ((1.0 / 130.0) * sqrt30 * Y8_p2 + (1.0 / 130.0) * sqrt2002 * Y8_p6 + (-3.0 / 55.0) * sqrt10 * Y6_p2 * R2 + (-1.0 / 55.0) * sqrt22 * Y6_p6 * R2 + (3.0 / 143.0) * sqrt105 * Y4_p2 * R4) * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_p2 + (3.0 / 22.0) * sqrt22 * Y6_p6 + (-3.0 / 11.0) * sqrt105 * Y4_p2 * R2) * V[1] + 0.75 * sqrt105 * Y4_p2 * V[2] + ((8.0 / 65.0) * sqrt30 * Y8_p2 + (8.0 / 65.0) * sqrt2002 * Y8_p6 + (-48.0 / 55.0) * sqrt10 * Y6_p2 * R2 + (-16.0 / 55.0) * sqrt22 * Y6_p6 * R2 + (48.0 / 143.0) * sqrt105 * Y4_p2 * R4) * V[3] + ((63.0 / 11.0) * sqrt10 * Y6_p2 + (21.0 / 11.0) * sqrt22 * Y6_p6 + (-42.0 / 11.0) * sqrt105 * Y4_p2 * R2) * V[4] + 9.0 * sqrt105 * Y4_p2 * V[5];
    d[51] = ((1.0 / 65.0) * sqrt429 * Y8_n5 + (1.0 / 65.0) * sqrt55 * Y8_n3 + (3.0 / 110.0) * sqrt66 * Y6_n5 * R2 + (-7.0 / 110.0) * sqrt10 * Y6_n3 * R2 + (-6.0 / 143.0) * sqrt30 * Y4_n3 * R4) * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_n5 + (21.0 / 44.0) * sqrt10 * Y6_n3 + (6.0 / 11.0) * sqrt30 * Y4_n3 * R2) * V[1] - 1.5 * sqrt30 * Y4_n3 * V[2] + ((16.0 / 65.0) * sqrt429 * Y8_n5 + (16.0 / 65.0) * sqrt55 * Y8_n3 + (24.0 / 55.0) * sqrt66 * Y6_n5 * R2 + (-56.0 / 55.0) * sqrt10 * Y6_n3 * R2 + (-96.0 / 143.0) * sqrt30 * Y4_n3 * R4) * V[3] + ((-63.0 / 22.0) * sqrt66 * Y6_n5 + (147.0 / 22.0) * sqrt10 * Y6_n3 + (84.0 / 11.0) * sqrt30 * Y4_n3 * R2) * V[4] - 18.0 * sqrt30 * Y4_n3 * V[5];
    d[50] = ((1.0 / 130.0) * sqrt2002 * Y8_n6 + (-1.0 / 130.0) * sqrt30 * Y8_n2 + (-1.0 / 55.0) * sqrt22 * Y6_n6 * R2 + (3.0 / 55.0) * sqrt10 * Y6_n2 * R2 + (-3.0 / 143.0) * sqrt105 * Y4_n2 * R4) * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_n6 + (-9.0 / 22.0) * sqrt10 * Y6_n2 + (3.0 / 11.0) * sqrt105 * Y4_n2 * R2) * V[1] - 0.75 * sqrt105 * Y4_n2 * V[2] + ((8.0 / 65.0) * sqrt2002 * Y8_n6 + (-8.0 / 65.0) * sqrt30 * Y8_n2 + (-16.0 / 55.0) * sqrt22 * Y6_n6 * R2 + (48.0 / 55.0) * sqrt10 * Y6_n2 * R2 + (-48.0 / 143.0) * sqrt105 * Y4_n2 * R4) * V[3] + ((21.0 / 11.0) * sqrt22 * Y6_n6 + (-63.0 / 11.0) * sqrt10 * Y6_n2 + (42.0 / 11.0) * sqrt105 * Y4_n2 * R2) * V[4] - 9.0 * sqrt105 * Y4_n2 * V[5];
    d[57] = ((3.0 / 65.0) * sqrt26 * Y8_p5 - 0.2 * Y6_p5 * R2) * V[0] + 1.5 * Y6_p5 * V[1] + ((48.0 / 65.0) * sqrt26 * Y8_p5 - 3.2 * Y6_p5 * R2) * V[3] + 21.0 * Y6_p5 * V[4];
    d[58] = ((-2.0 / 65.0) * sqrt6 * Y8_p4 + (2.0 / 65.0) * sqrt91 * Y8_p6 + (3.0 / 110.0) * sqrt66 * Y6_p4 * R2 + 0.2 * Y6_p6 * R2 + (3.0 / 286.0) * sqrt330 * Y4_p4 * R4) * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_p4 - 1.5 * Y6_p6 + (-3.0 / 22.0) * sqrt330 * Y4_p4 * R2) * V[1] + 0.375 * sqrt330 * Y4_p4 * V[2] + ((-32.0 / 65.0) * sqrt6 * Y8_p4 + (32.0 / 65.0) * sqrt91 * Y8_p6 + (24.0 / 55.0) * sqrt66 * Y6_p4 * R2 + 3.2 * Y6_p6 * R2 + (24.0 / 143.0) * sqrt330 * Y4_p4 * R4) * V[3] + ((-63.0 / 22.0) * sqrt66 * Y6_p4 - 21.0 * Y6_p6 + (-21.0 / 11.0) * sqrt330 * Y4_p4 * R2) * V[4] + 4.5 * sqrt330 * Y4_p4 * V[5];
    d[59] = ((1.0 / 130.0) * sqrt10 * Y8_p3 + (1.0 / 130.0) * sqrt2730 * Y8_p7 + (-1.0 / 55.0) * sqrt55 * Y6_p3 * R2 + (3.0 / 143.0) * sqrt165 * Y4_p3 * R4) * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_p3 + (-3.0 / 11.0) * sqrt165 * Y4_p3 * R2) * V[1] + 0.75 * sqrt165 * Y4_p3 * V[2] + ((8.0 / 65.0) * sqrt10 * Y8_p3 + (8.0 / 65.0) * sqrt2730 * Y8_p7 + (-16.0 / 55.0) * sqrt55 * Y6_p3 * R2 + (48.0 / 143.0) * sqrt165 * Y4_p3 * R4) * V[3] + ((21.0 / 11.0) * sqrt55 * Y6_p3 + (-42.0 / 11.0) * sqrt165 * Y4_p3 * R2) * V[4] + 9.0 * sqrt165 * Y4_p3 * V[5];
    d[56] = ((2.0 / 65.0) * sqrt91 * Y8_n6 + (2.0 / 65.0) * sqrt6 * Y8_n4 + 0.2 * Y6_n6 * R2 + (-3.0 / 110.0) * sqrt66 * Y6_n4 * R2 + (-3.0 / 286.0) * sqrt330 * Y4_n4 * R4) * V[0] + (-1.5 * Y6_n6 + (9.0 / 44.0) * sqrt66 * Y6_n4 + (3.0 / 22.0) * sqrt330 * Y4_n4 * R2) * V[1] - 0.375 * sqrt330 * Y4_n4 * V[2] + ((32.0 / 65.0) * sqrt91 * Y8_n6 + (32.0 / 65.0) * sqrt6 * Y8_n4 + 3.2 * Y6_n6 * R2 + (-24.0 / 55.0) * sqrt66 * Y6_n4 * R2 + (-24.0 / 143.0) * sqrt330 * Y4_n4 * R4) * V[3] + (-21.0 * Y6_n6 + (63.0 / 22.0) * sqrt66 * Y6_n4 + (21.0 / 11.0) * sqrt330 * Y4_n4 * R2) * V[4] - 4.5 * sqrt330 * Y4_n4 * V[5];
    d[55] = ((1.0 / 130.0) * sqrt2730 * Y8_n7 + (-1.0 / 130.0) * sqrt10 * Y8_n3 + (1.0 / 55.0) * sqrt55 * Y6_n3 * R2 + (-3.0 / 143.0) * sqrt165 * Y4_n3 * R4) * V[0] + ((-3.0 / 22.0) * sqrt55 * Y6_n3 + (3.0 / 11.0) * sqrt165 * Y4_n3 * R2) * V[1] - 0.75 * sqrt165 * Y4_n3 * V[2] + ((8.0 / 65.0) * sqrt2730 * Y8_n7 + (-8.0 / 65.0) * sqrt10 * Y8_n3 + (16.0 / 55.0) * sqrt55 * Y6_n3 * R2 + (-48.0 / 143.0) * sqrt165 * Y4_n3 * R4) * V[3] + ((-21.0 / 11.0) * sqrt55 * Y6_n3 + (42.0 / 11.0) * sqrt165 * Y4_n3 * R2) * V[4] - 9.0 * sqrt165 * Y4_n3 * V[5];
    d[62] = ((1.0 / 65.0) * sqrt91 * Y8_p6 - 0.4 * Y6_p6 * R2) * V[0] + 3.0 * Y6_p6 * V[1] + ((16.0 / 65.0) * sqrt91 * Y8_p6 - 6.4 * Y6_p6 * R2) * V[3] + 42.0 * Y6_p6 * V[4];
    d[63] = ((-1.0 / 130.0) * sqrt26 * Y8_p5 + (1.0 / 130.0) * sqrt910 * Y8_p7 + 0.2 * Y6_p5 * R2) * V[0] - 1.5 * Y6_p5 * V[1] + ((-8.0 / 65.0) * sqrt26 * Y8_p5 + (8.0 / 65.0) * sqrt910 * Y8_p7 + 3.2 * Y6_p5 * R2) * V[3] - 21.0 * Y6_p5 * V[4];
    d[64] = ((1.0 / 130.0) * sqrt2 * Y8_p4 + (1.0 / 65.0) * sqrt910 * Y8_p8 + (-1.0 / 55.0) * sqrt22 * Y6_p4 * R2 + (9.0 / 286.0) * sqrt110 * Y4_p4 * R4) * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_p4 + (-9.0 / 22.0) * sqrt110 * Y4_p4 * R2) * V[1] + 1.125 * sqrt110 * Y4_p4 * V[2] + ((8.0 / 65.0) * sqrt2 * Y8_p4 + (16.0 / 65.0) * sqrt910 * Y8_p8 + (-16.0 / 55.0) * sqrt22 * Y6_p4 * R2 + (72.0 / 143.0) * sqrt110 * Y4_p4 * R4) * V[3] + ((21.0 / 11.0) * sqrt22 * Y6_p4 + (-63.0 / 11.0) * sqrt110 * Y4_p4 * R2) * V[4] + 13.5 * sqrt110 * Y4_p4 * V[5];
    d[61] = ((1.0 / 130.0) * sqrt910 * Y8_n7 + (1.0 / 130.0) * sqrt26 * Y8_n5 - 0.2 * Y6_n5 * R2) * V[0] + 1.5 * Y6_n5 * V[1] + ((8.0 / 65.0) * sqrt910 * Y8_n7 + (8.0 / 65.0) * sqrt26 * Y8_n5 - 3.2 * Y6_n5 * R2) * V[3] + 21.0 * Y6_n5 * V[4];
    d[60] = ((1.0 / 65.0) * sqrt910 * Y8_n8 + (-1.0 / 130.0) * sqrt2 * Y8_n4 + (1.0 / 55.0) * sqrt22 * Y6_n4 * R2 + (-9.0 / 286.0) * sqrt110 * Y4_n4 * R4) * V[0] + ((-3.0 / 22.0) * sqrt22 * Y6_n4 + (9.0 / 22.0) * sqrt110 * Y4_n4 * R2) * V[1] - 1.125 * sqrt110 * Y4_n4 * V[2] + ((16.0 / 65.0) * sqrt910 * Y8_n8 + (-8.0 / 65.0) * sqrt2 * Y8_n4 + (16.0 / 55.0) * sqrt22 * Y6_n4 * R2 + (-72.0 / 143.0) * sqrt110 * Y4_n4 * R4) * V[3] + ((-21.0 / 11.0) * sqrt22 * Y6_n4 + (63.0 / 11.0) * sqrt110 * Y4_n4 * R2) * V[4] - 13.5 * sqrt110 * Y4_n4 * V[5];
    d[27] = ((6.0 / 65.0) * sqrt21 * Y8_n1 + (13.0 / 55.0) * Y6_n1 * R2 + (3.0 / 143.0) * sqrt210 * Y4_n1 * R4) * V[0] + ((-39.0 / 22.0) * Y6_n1 + (-3.0 / 11.0) * sqrt210 * Y4_n1 * R2) * V[1] + 0.75 * sqrt210 * Y4_n1 * V[2] + ((96.0 / 65.0) * sqrt21 * Y8_n1 + (208.0 / 55.0) * Y6_n1 * R2 + (48.0 / 143.0) * sqrt210 * Y4_n1 * R4) * V[3] + ((-273.0 / 11.0) * Y6_n1 + (-42.0 / 11.0) * sqrt210 * Y4_n1 * R2) * V[4] + 9.0 * sqrt210 * Y4_n1 * V[5];
    d[28] = ((6.0 / 65.0) * sqrt10 * Y8_n2 + (1.0 / 55.0) * sqrt30 * Y6_n2 * R2 + (-3.0 / 143.0) * sqrt35 * Y4_n2 * R4) * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_n2 + (3.0 / 11.0) * sqrt35 * Y4_n2 * R2) * V[1] - 0.75 * sqrt35 * Y4_n2 * V[2] + ((96.0 / 65.0) * sqrt10 * Y8_n2 + (16.0 / 55.0) * sqrt30 * Y6_n2 * R2 + (-48.0 / 143.0) * sqrt35 * Y4_n2 * R4) * V[3] + ((-21.0 / 11.0) * sqrt30 * Y6_n2 + (42.0 / 11.0) * sqrt35 * Y4_n2 * R2) * V[4] - 9.0 * sqrt35 * Y4_n2 * V[5];
    d[29] = ((1.0 / 65.0) * sqrt165 * Y8_n3 + (3.0 / 65.0) * sqrt7 * Y8_n1 + (-2.0 / 55.0) * sqrt30 * Y6_n3 * R2 + (-7.0 / 55.0) * sqrt3 * Y6_n1 * R2 + (3.0 / 286.0) * sqrt10 * Y4_n3 * R4 + (3.0 / 286.0) * sqrt70 * Y4_n1 * R4) * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_n3 + (21.0 / 22.0) * sqrt3 * Y6_n1 + (-3.0 / 22.0) * sqrt10 * Y4_n3 * R2 + (-3.0 / 22.0) * sqrt70 * Y4_n1 * R2) * V[1] + (0.375 * sqrt10 * Y4_n3 + 0.375 * sqrt70 * Y4_n1) * V[2] + ((16.0 / 65.0) * sqrt165 * Y8_n3 + (48.0 / 65.0) * sqrt7 * Y8_n1 + (-32.0 / 55.0) * sqrt30 * Y6_n3 * R2 + (-112.0 / 55.0) * sqrt3 * Y6_n1 * R2 + (24.0 / 143.0) * sqrt10 * Y4_n3 * R4 + (24.0 / 143.0) * sqrt70 * Y4_n1 * R4) * V[3] + ((42.0 / 11.0) * sqrt30 * Y6_n3 + (147.0 / 11.0) * sqrt3 * Y6_n1 + (-21.0 / 11.0) * sqrt10 * Y4_n3 * R2 + (-21.0 / 11.0) * sqrt70 * Y4_n1 * R2) * V[4] + (4.5 * sqrt10 * Y4_n3 + 4.5 * sqrt70 * Y4_n1) * V[5];
    d[26] = ((-8.0 / 65.0) * sqrt7 * Y8_p0 + (-6.0 / 65.0) * sqrt10 * Y8_p2 + (1.0 / 55.0) * sqrt7 * Y6_p0 * R2 + (-1.0 / 55.0) * sqrt30 * Y6_p2 * R2 + (15.0 / 143.0) * sqrt7 * Y4_p0 * R4 + (3.0 / 143.0) * sqrt35 * Y4_p2 * R4) * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_p0 + (3.0 / 22.0) * sqrt30 * Y6_p2 + (-15.0 / 11.0) * sqrt7 * Y4_p0 * R2 + (-3.0 / 11.0) * sqrt35 * Y4_p2 * R2) * V[1] + (3.75 * sqrt7 * Y4_p0 + 0.75 * sqrt35 * Y4_p2) * V[2] + ((-128.0 / 65.0) * sqrt7 * Y8_p0 + (-96.0 / 65.0) * sqrt10 * Y8_p2 + (16.0 / 55.0) * sqrt7 * Y6_p0 * R2 + (-16.0 / 55.0) * sqrt30 * Y6_p2 * R2 + (240.0 / 143.0) * sqrt7 * Y4_p0 * R4 + (48.0 / 143.0) * sqrt35 * Y4_p2 * R4) * V[3] + ((-21.0 / 11.0) * sqrt7 * Y6_p0 + (21.0 / 11.0) * sqrt30 * Y6_p2 + (-210.0 / 11.0) * sqrt7 * Y4_p0 * R2 + (-42.0 / 11.0) * sqrt35 * Y4_p2 * R2) * V[4] + (45.0 * sqrt7 * Y4_p0 + 9.0 * sqrt35 * Y4_p2) * V[5];
    d[25] = ((-3.0 / 65.0) * sqrt7 * Y8_p1 + (-1.0 / 65.0) * sqrt165 * Y8_p3 + (7.0 / 55.0) * sqrt3 * Y6_p1 * R2 + (2.0 / 55.0) * sqrt30 * Y6_p3 * R2 + (-3.0 / 286.0) * sqrt70 * Y4_p1 * R4 + (-3.0 / 286.0) * sqrt10 * Y4_p3 * R4) * V[0] + ((-21.0 / 22.0) * sqrt3 * Y6_p1 + (-3.0 / 11.0) * sqrt30 * Y6_p3 + (3.0 / 22.0) * sqrt70 * Y4_p1 * R2 + (3.0 / 22.0) * sqrt10 * Y4_p3 * R2) * V[1] + (-0.375 * sqrt70 * Y4_p1 - 0.375 * sqrt10 * Y4_p3) * V[2] + ((-48.0 / 65.0) * sqrt7 * Y8_p1 + (-16.0 / 65.0) * sqrt165 * Y8_p3 + (112.0 / 55.0) * sqrt3 * Y6_p1 * R2 + (32.0 / 55.0) * sqrt30 * Y6_p3 * R2 + (-24.0 / 143.0) * sqrt70 * Y4_p1 * R4 + (-24.0 / 143.0) * sqrt10 * Y4_p3 * R4) * V[3] + ((-147.0 / 11.0) * sqrt3 * Y6_p1 + (-42.0 / 11.0) * sqrt30 * Y6_p3 + (21.0 / 11.0) * sqrt70 * Y4_p1 * R2 + (21.0 / 11.0) * sqrt10 * Y4_p3 * R2) * V[4] + (-4.5 * sqrt70 * Y4_p1 - 4.5 * sqrt10 * Y4_p3) * V[5];
    d[22] = ((3.0 / 13.0) * sqrt3 * Y8_n2 + (2.0 / 11.0) * Y6_n2 * R2 + (6.0 / 143.0) * sqrt42 * Y4_n2 * R4) * V[0] + ((-15.0 / 11.0) * Y6_n2 + (-6.0 / 11.0) * sqrt42 * Y4_n2 * R2) * V[1] + 1.5 * sqrt42 * Y4_n2 * V[2] + ((48.0 / 13.0) * sqrt3 * Y8_n2 + (32.0 / 11.0) * Y6_n2 * R2 + (96.0 / 143.0) * sqrt42 * Y4_n2 * R4) * V[3] + ((-210.0 / 11.0) * Y6_n2 + (-84.0 / 11.0) * sqrt42 * Y4_n2 * R2) * V[4] + 18.0 * sqrt42 * Y4_n2 * V[5];
    d[23] = ((1.0 / 26.0) * sqrt66 * Y8_n3 + (-3.0 / 130.0) * sqrt70 * Y8_n1 + (1.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (1.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (-12.0 / 143.0) * Y4_n3 * R4 + (12.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_n3 + (-3.0 / 22.0) * sqrt30 * Y6_n1 + (12.0 / 11.0) * Y4_n3 * R2 + (-12.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1] + (-3.0 * Y4_n3 + 3.0 * sqrt7 * Y4_n1) * V[2] + ((8.0 / 13.0) * sqrt66 * Y8_n3 + (-24.0 / 65.0) * sqrt70 * Y8_n1 + (16.0 / 11.0) * sqrt3 * Y6_n3 * R2 + (16.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (-192.0 / 143.0) * Y4_n3 * R4 + (192.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[3] + ((-105.0 / 11.0) * sqrt3 * Y6_n3 + (-21.0 / 11.0) * sqrt30 * Y6_n1 + (168.0 / 11.0) * Y4_n3 * R2 + (-168.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[4] + (-36.0 * Y4_n3 + 36.0 * sqrt7 * Y4_n1) * V[5];
    d[24] = ((3.0 / 130.0) * sqrt110 * Y8_n4 + (-3.0 / 55.0) * sqrt10 * Y6_n4 * R2 + (3.0 / 286.0) * sqrt2 * Y4_n4 * R4) * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_n4 + (-3.0 / 22.0) * sqrt2 * Y4_n4 * R2) * V[1] + 0.375 * sqrt2 * Y4_n4 * V[2] + ((24.0 / 65.0) * sqrt110 * Y8_n4 + (-48.0 / 55.0) * sqrt10 * Y6_n4 * R2 + (24.0 / 143.0) * sqrt2 * Y4_n4 * R4) * V[3] + ((63.0 / 11.0) * sqrt10 * Y6_n4 + (-21.0 / 11.0) * sqrt2 * Y4_n4 * R2) * V[4] + 4.5 * sqrt2 * Y4_n4 * V[5];
    d[21] = ((-3.0 / 130.0) * sqrt70 * Y8_p1 + (-1.0 / 26.0) * sqrt66 * Y8_p3 + (1.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (-1.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (12.0 / 143.0) * sqrt7 * Y4_p1 * R4 + (12.0 / 143.0) * Y4_p3 * R4) * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_p1 + (15.0 / 22.0) * sqrt3 * Y6_p3 + (-12.0 / 11.0) * sqrt7 * Y4_p1 * R2 + (-12.0 / 11.0) * Y4_p3 * R2) * V[1] + (3.0 * sqrt7 * Y4_p1 + 3.0 * Y4_p3) * V[2] + ((-24.0 / 65.0) * sqrt70 * Y8_p1 + (-8.0 / 13.0) * sqrt66 * Y8_p3 + (16.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (-16.0 / 11.0) * sqrt3 * Y6_p3 * R2 + (192.0 / 143.0) * sqrt7 * Y4_p1 * R4 + (192.0 / 143.0) * Y4_p3 * R4) * V[3] + ((-21.0 / 11.0) * sqrt30 * Y6_p1 + (105.0 / 11.0) * sqrt3 * Y6_p3 + (-168.0 / 11.0) * sqrt7 * Y4_p1 * R2 + (-168.0 / 11.0) * Y4_p3 * R2) * V[4] + (36.0 * sqrt7 * Y4_p1 + 36.0 * Y4_p3) * V[5];
    d[20] = ((1.0 / 65.0) * sqrt70 * Y8_p0 + (-3.0 / 130.0) * sqrt110 * Y8_p4 + (-2.0 / 55.0) * sqrt70 * Y6_p0 * R2 + (3.0 / 55.0) * sqrt10 * Y6_p4 * R2 + (3.0 / 143.0) * sqrt70 * Y4_p0 * R4 + (-3.0 / 286.0) * sqrt2 * Y4_p4 * R4) * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_p0 + (-9.0 / 22.0) * sqrt10 * Y6_p4 + (-3.0 / 11.0) * sqrt70 * Y4_p0 * R2 + (3.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1] + (0.75 * sqrt70 * Y4_p0 - 0.375 * sqrt2 * Y4_p4) * V[2] + ((16.0 / 65.0) * sqrt70 * Y8_p0 + (-24.0 / 65.0) * sqrt110 * Y8_p4 + (-32.0 / 55.0) * sqrt70 * Y6_p0 * R2 + (48.0 / 55.0) * sqrt10 * Y6_p4 * R2 + (48.0 / 143.0) * sqrt70 * Y4_p0 * R4 + (-24.0 / 143.0) * sqrt2 * Y4_p4 * R4) * V[3] + ((42.0 / 11.0) * sqrt70 * Y6_p0 + (-63.0 / 11.0) * sqrt10 * Y6_p4 + (-42.0 / 11.0) * sqrt70 * Y4_p0 * R2 + (21.0 / 11.0) * sqrt2 * Y4_p4 * R2) * V[4] + (9.0 * sqrt70 * Y4_p0 - 4.5 * sqrt2 * Y4_p4) * V[5];
    d[17] = ((1.0 / 13.0) * sqrt22 * Y8_n3 + (1.0 / 11.0) * Y6_n3 * R2 + (18.0 / 143.0) * sqrt3 * Y4_n3 * R4) * V[0] + ((-15.0 / 22.0) * Y6_n3 + (-18.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[1] + 4.5 * sqrt3 * Y4_n3 * V[2] + ((16.0 / 13.0) * sqrt22 * Y8_n3 + (16.0 / 11.0) * Y6_n3 * R2 + (288.0 / 143.0) * sqrt3 * Y4_n3 * R4) * V[3] + ((-105.0 / 11.0) * Y6_n3 + (-252.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[4] + 54.0 * sqrt3 * Y4_n3 * V[5];
    d[18] = ((2.0 / 65.0) * sqrt110 * Y8_n4 + (-2.0 / 13.0) * Y8_n2 + (7.0 / 110.0) * sqrt10 * Y6_n4 * R2 + (1.0 / 11.0) * sqrt3 * Y6_n2 * R2 + (-9.0 / 286.0) * sqrt2 * Y4_n4 * R4 + (9.0 / 143.0) * sqrt14 * Y4_n2 * R4) * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_n4 + (-15.0 / 22.0) * sqrt3 * Y6_n2 + (9.0 / 22.0) * sqrt2 * Y4_n4 * R2 + (-9.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[1] + (-1.125 * sqrt2 * Y4_n4 + 2.25 * sqrt14 * Y4_n2) * V[2] + ((32.0 / 65.0) * sqrt110 * Y8_n4 + (-32.0 / 13.0) * Y8_n2 + (56.0 / 55.0) * sqrt10 * Y6_n4 * R2 + (16.0 / 11.0) * sqrt3 * Y6_n2 * R2 + (-72.0 / 143.0) * sqrt2 * Y4_n4 * R4 + (144.0 / 143.0) * sqrt14 * Y4_n2 * R4) * V[3] + ((-147.0 / 22.0) * sqrt10 * Y6_n4 + (-105.0 / 11.0) * sqrt3 * Y6_n2 + (63.0 / 11.0) * sqrt2 * Y4_n4 * R2 + (-126.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[4] + (-13.5 * sqrt2 * Y4_n4 + 27.0 * sqrt14 * Y4_n2) * V[5];
    d[19] = ((1.0 / 130.0) * sqrt1430 * Y8_n5 + (1.0 / 130.0) * sqrt70 * Y8_n1 + (-1.0 / 55.0) * sqrt55 * Y6_n5 * R2 + (-2.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (9.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_n5 + (3.0 / 11.0) * sqrt30 * Y6_n1 + (-9.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1] + 2.25 * sqrt7 * Y4_n1 * V[2] + ((8.0 / 65.0) * sqrt1430 * Y8_n5 + (8.0 / 65.0) * sqrt70 * Y8_n1 + (-16.0 / 55.0) * sqrt55 * Y6_n5 * R2 + (-32.0 / 55.0) * sqrt30 * Y6_n1 * R2 + (144.0 / 143.0) * sqrt7 * Y4_n1 * R4) * V[3] + ((21.0 / 11.0) * sqrt55 * Y6_n5 + (42.0 / 11.0) * sqrt30 * Y6_n1 + (-126.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[4] + 27.0 * sqrt7 * Y4_n1 * V[5];
    d[16] = ((-2.0 / 13.0) * Y8_p2 + (-2.0 / 65.0) * sqrt110 * Y8_p4 + (1.0 / 11.0) * sqrt3 * Y6_p2 * R2 + (-7.0 / 110.0) * sqrt10 * Y6_p4 * R2 + (9.0 / 143.0) * sqrt14 * Y4_p2 * R4 + (9.0 / 286.0) * sqrt2 * Y4_p4 * R4) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_p2 + (21.0 / 44.0) * sqrt10 * Y6_p4 + (-9.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (-9.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1] + (2.25 * sqrt14 * Y4_p2 + 1.125 * sqrt2 * Y4_p4) * V[2] + ((-32.0 / 13.0) * Y8_p2 + (-32.0 / 65.0) * sqrt110 * Y8_p4 + (16.0 / 11.0) * sqrt3 * Y6_p2 * R2 + (-56.0 / 55.0) * sqrt10 * Y6_p4 * R2 + (144.0 / 143.0) * sqrt14 * Y4_p2 * R4 + (72.0 / 143.0) * sqrt2 * Y4_p4 * R4) * V[3] + ((-105.0 / 11.0) * sqrt3 * Y6_p2 + (147.0 / 22.0) * sqrt10 * Y6_p4 + (-126.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (-63.0 / 11.0) * sqrt2 * Y4_p4 * R2) * V[4] + (27.0 * sqrt14 * Y4_p2 + 13.5 * sqrt2 * Y4_p4) * V[5];
    d[15] = ((1.0 / 130.0) * sqrt70 * Y8_p1 + (-1.0 / 130.0) * sqrt1430 * Y8_p5 + (-2.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (1.0 / 55.0) * sqrt55 * Y6_p5 * R2 + (9.0 / 143.0) * sqrt7 * Y4_p1 * R4) * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_p1 + (-3.0 / 22.0) * sqrt55 * Y6_p5 + (-9.0 / 11.0) * sqrt7 * Y4_p1 * R2) * V[1] + 2.25 * sqrt7 * Y4_p1 * V[2] + ((8.0 / 65.0) * sqrt70 * Y8_p1 + (-8.0 / 65.0) * sqrt1430 * Y8_p5 + (-32.0 / 55.0) * sqrt30 * Y6_p1 * R2 + (16.0 / 55.0) * sqrt55 * Y6_p5 * R2 + (144.0 / 143.0) * sqrt7 * Y4_p1 * R4) * V[3] + ((42.0 / 11.0) * sqrt30 * Y6_p1 + (-21.0 / 11.0) * sqrt55 * Y6_p5 + (-126.0 / 11.0) * sqrt7 * Y4_p1 * R2) * V[4] + 27.0 * sqrt7 * Y4_p1 * V[5];
    d[12] = ((6.0 / 65.0) * sqrt11 * Y8_n4 + (-2.0 / 55.0) * Y6_n4 * R2 + (9.0 / 143.0) * sqrt5 * Y4_n4 * R4) * V[0] + ((3.0 / 11.0) * Y6_n4 + (-9.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[1] + 2.25 * sqrt5 * Y4_n4 * V[2] + ((96.0 / 65.0) * sqrt11 * Y8_n4 + (-32.0 / 55.0) * Y6_n4 * R2 + (144.0 / 143.0) * sqrt5 * Y4_n4 * R4) * V[3] + ((42.0 / 11.0) * Y6_n4 + (-126.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[4] + 27.0 * sqrt5 * Y4_n4 * V[5];
    d[13] = ((1.0 / 65.0) * sqrt429 * Y8_n5 + (-1.0 / 65.0) * sqrt55 * Y8_n3 + (3.0 / 110.0) * sqrt66 * Y6_n5 * R2 + (7.0 / 110.0) * sqrt10 * Y6_n3 * R2 + (6.0 / 143.0) * sqrt30 * Y4_n3 * R4) * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_n5 + (-21.0 / 44.0) * sqrt10 * Y6_n3 + (-6.0 / 11.0) * sqrt30 * Y4_n3 * R2) * V[1] + 1.5 * sqrt30 * Y4_n3 * V[2] + ((16.0 / 65.0) * sqrt429 * Y8_n5 + (-16.0 / 65.0) * sqrt55 * Y8_n3 + (24.0 / 55.0) * sqrt66 * Y6_n5 * R2 + (56.0 / 55.0) * sqrt10 * Y6_n3 * R2 + (96.0 / 143.0) * sqrt30 * Y4_n3 * R4) * V[3] + ((-63.0 / 22.0) * sqrt66 * Y6_n5 + (-147.0 / 22.0) * sqrt10 * Y6_n3 + (-84.0 / 11.0) * sqrt30 * Y4_n3 * R2) * V[4] + 18.0 * sqrt30 * Y4_n3 * V[5];
    d[14] = ((1.0 / 130.0) * sqrt2002 * Y8_n6 + (1.0 / 130.0) * sqrt30 * Y8_n2 + (-1.0 / 55.0) * sqrt22 * Y6_n6 * R2 + (-3.0 / 55.0) * sqrt10 * Y6_n2 * R2 + (3.0 / 143.0) * sqrt105 * Y4_n2 * R4) * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_n6 + (9.0 / 22.0) * sqrt10 * Y6_n2 + (-3.0 / 11.0) * sqrt105 * Y4_n2 * R2) * V[1] + 0.75 * sqrt105 * Y4_n2 * V[2] + ((8.0 / 65.0) * sqrt2002 * Y8_n6 + (8.0 / 65.0) * sqrt30 * Y8_n2 + (-16.0 / 55.0) * sqrt22 * Y6_n6 * R2 + (-48.0 / 55.0) * sqrt10 * Y6_n2 * R2 + (48.0 / 143.0) * sqrt105 * Y4_n2 * R4) * V[3] + ((21.0 / 11.0) * sqrt22 * Y6_n6 + (63.0 / 11.0) * sqrt10 * Y6_n2 + (-42.0 / 11.0) * sqrt105 * Y4_n2 * R2) * V[4] + 9.0 * sqrt105 * Y4_n2 * V[5];
    d[11] = ((-1.0 / 65.0) * sqrt55 * Y8_p3 + (-1.0 / 65.0) * sqrt429 * Y8_p5 + (7.0 / 110.0) * sqrt10 * Y6_p3 * R2 + (-3.0 / 110.0) * sqrt66 * Y6_p5 * R2 + (6.0 / 143.0) * sqrt30 * Y4_p3 * R4) * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_p3 + (9.0 / 44.0) * sqrt66 * Y6_p5 + (-6.0 / 11.0) * sqrt30 * Y4_p3 * R2) * V[1] + 1.5 * sqrt30 * Y4_p3 * V[2] + ((-16.0 / 65.0) * sqrt55 * Y8_p3 + (-16.0 / 65.0) * sqrt429 * Y8_p5 + (56.0 / 55.0) * sqrt10 * Y6_p3 * R2 + (-24.0 / 55.0) * sqrt66 * Y6_p5 * R2 + (96.0 / 143.0) * sqrt30 * Y4_p3 * R4) * V[3] + ((-147.0 / 22.0) * sqrt10 * Y6_p3 + (63.0 / 22.0) * sqrt66 * Y6_p5 + (-84.0 / 11.0) * sqrt30 * Y4_p3 * R2) * V[4] + 18.0 * sqrt30 * Y4_p3 * V[5];
    d[10] = ((1.0 / 130.0) * sqrt30 * Y8_p2 + (-1.0 / 130.0) * sqrt2002 * Y8_p6 + (-3.0 / 55.0) * sqrt10 * Y6_p2 * R2 + (1.0 / 55.0) * sqrt22 * Y6_p6 * R2 + (3.0 / 143.0) * sqrt105 * Y4_p2 * R4) * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_p2 + (-3.0 / 22.0) * sqrt22 * Y6_p6 + (-3.0 / 11.0) * sqrt105 * Y4_p2 * R2) * V[1] + 0.75 * sqrt105 * Y4_p2 * V[2] + ((8.0 / 65.0) * sqrt30 * Y8_p2 + (-8.0 / 65.0) * sqrt2002 * Y8_p6 + (-48.0 / 55.0) * sqrt10 * Y6_p2 * R2 + (16.0 / 55.0) * sqrt22 * Y6_p6 * R2 + (48.0 / 143.0) * sqrt105 * Y4_p2 * R4) * V[3] + ((63.0 / 11.0) * sqrt10 * Y6_p2 + (-21.0 / 11.0) * sqrt22 * Y6_p6 + (-42.0 / 11.0) * sqrt105 * Y4_p2 * R2) * V[4] + 9.0 * sqrt105 * Y4_p2 * V[5];
    d[7] = ((3.0 / 65.0) * sqrt26 * Y8_n5 - 0.2 * Y6_n5 * R2) * V[0] + 1.5 * Y6_n5 * V[1] + ((48.0 / 65.0) * sqrt26 * Y8_n5 - 3.2 * Y6_n5 * R2) * V[3] + 21.0 * Y6_n5 * V[4];
    d[8] = ((2.0 / 65.0) * sqrt91 * Y8_n6 + (-2.0 / 65.0) * sqrt6 * Y8_n4 + 0.2 * Y6_n6 * R2 + (3.0 / 110.0) * sqrt66 * Y6_n4 * R2 + (3.0 / 286.0) * sqrt330 * Y4_n4 * R4) * V[0] + (-1.5 * Y6_n6 + (-9.0 / 44.0) * sqrt66 * Y6_n4 + (-3.0 / 22.0) * sqrt330 * Y4_n4 * R2) * V[1] + 0.375 * sqrt330 * Y4_n4 * V[2] + ((32.0 / 65.0) * sqrt91 * Y8_n6 + (-32.0 / 65.0) * sqrt6 * Y8_n4 + 3.2 * Y6_n6 * R2 + (24.0 / 55.0) * sqrt66 * Y6_n4 * R2 + (24.0 / 143.0) * sqrt330 * Y4_n4 * R4) * V[3] + (-21.0 * Y6_n6 + (-63.0 / 22.0) * sqrt66 * Y6_n4 + (-21.0 / 11.0) * sqrt330 * Y4_n4 * R2) * V[4] + 4.5 * sqrt330 * Y4_n4 * V[5];
    d[9] = ((1.0 / 130.0) * sqrt2730 * Y8_n7 + (1.0 / 130.0) * sqrt10 * Y8_n3 + (-1.0 / 55.0) * sqrt55 * Y6_n3 * R2 + (3.0 / 143.0) * sqrt165 * Y4_n3 * R4) * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_n3 + (-3.0 / 11.0) * sqrt165 * Y4_n3 * R2) * V[1] + 0.75 * sqrt165 * Y4_n3 * V[2] + ((8.0 / 65.0) * sqrt2730 * Y8_n7 + (8.0 / 65.0) * sqrt10 * Y8_n3 + (-16.0 / 55.0) * sqrt55 * Y6_n3 * R2 + (48.0 / 143.0) * sqrt165 * Y4_n3 * R4) * V[3] + ((21.0 / 11.0) * sqrt55 * Y6_n3 + (-42.0 / 11.0) * sqrt165 * Y4_n3 * R2) * V[4] + 9.0 * sqrt165 * Y4_n3 * V[5];
    d[6] = ((-2.0 / 65.0) * sqrt6 * Y8_p4 + (-2.0 / 65.0) * sqrt91 * Y8_p6 + (3.0 / 110.0) * sqrt66 * Y6_p4 * R2 - 0.2 * Y6_p6 * R2 + (3.0 / 286.0) * sqrt330 * Y4_p4 * R4) * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_p4 + 1.5 * Y6_p6 + (-3.0 / 22.0) * sqrt330 * Y4_p4 * R2) * V[1] + 0.375 * sqrt330 * Y4_p4 * V[2] + ((-32.0 / 65.0) * sqrt6 * Y8_p4 + (-32.0 / 65.0) * sqrt91 * Y8_p6 + (24.0 / 55.0) * sqrt66 * Y6_p4 * R2 - 3.2 * Y6_p6 * R2 + (24.0 / 143.0) * sqrt330 * Y4_p4 * R4) * V[3] + ((-63.0 / 22.0) * sqrt66 * Y6_p4 + 21.0 * Y6_p6 + (-21.0 / 11.0) * sqrt330 * Y4_p4 * R2) * V[4] + 4.5 * sqrt330 * Y4_p4 * V[5];
    d[5] = ((1.0 / 130.0) * sqrt10 * Y8_p3 + (-1.0 / 130.0) * sqrt2730 * Y8_p7 + (-1.0 / 55.0) * sqrt55 * Y6_p3 * R2 + (3.0 / 143.0) * sqrt165 * Y4_p3 * R4) * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_p3 + (-3.0 / 11.0) * sqrt165 * Y4_p3 * R2) * V[1] + 0.75 * sqrt165 * Y4_p3 * V[2] + ((8.0 / 65.0) * sqrt10 * Y8_p3 + (-8.0 / 65.0) * sqrt2730 * Y8_p7 + (-16.0 / 55.0) * sqrt55 * Y6_p3 * R2 + (48.0 / 143.0) * sqrt165 * Y4_p3 * R4) * V[3] + ((21.0 / 11.0) * sqrt55 * Y6_p3 + (-42.0 / 11.0) * sqrt165 * Y4_p3 * R2) * V[4] + 9.0 * sqrt165 * Y4_p3 * V[5];
    d[2] = ((1.0 / 65.0) * sqrt91 * Y8_n6 - 0.4 * Y6_n6 * R2) * V[0] + 3.0 * Y6_n6 * V[1] + ((16.0 / 65.0) * sqrt91 * Y8_n6 - 6.4 * Y6_n6 * R2) * V[3] + 42.0 * Y6_n6 * V[4];
    d[3] = ((1.0 / 130.0) * sqrt910 * Y8_n7 + (-1.0 / 130.0) * sqrt26 * Y8_n5 + 0.2 * Y6_n5 * R2) * V[0] - 1.5 * Y6_n5 * V[1] + ((8.0 / 65.0) * sqrt910 * Y8_n7 + (-8.0 / 65.0) * sqrt26 * Y8_n5 + 3.2 * Y6_n5 * R2) * V[3] - 21.0 * Y6_n5 * V[4];
    d[4] = ((1.0 / 65.0) * sqrt910 * Y8_n8 + (1.0 / 130.0) * sqrt2 * Y8_n4 + (-1.0 / 55.0) * sqrt22 * Y6_n4 * R2 + (9.0 / 286.0) * sqrt110 * Y4_n4 * R4) * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_n4 + (-9.0 / 22.0) * sqrt110 * Y4_n4 * R2) * V[1] + 1.125 * sqrt110 * Y4_n4 * V[2] + ((16.0 / 65.0) * sqrt910 * Y8_n8 + (8.0 / 65.0) * sqrt2 * Y8_n4 + (-16.0 / 55.0) * sqrt22 * Y6_n4 * R2 + (72.0 / 143.0) * sqrt110 * Y4_n4 * R4) * V[3] + ((21.0 / 11.0) * sqrt22 * Y6_n4 + (-63.0 / 11.0) * sqrt110 * Y4_n4 * R2) * V[4] + 13.5 * sqrt110 * Y4_n4 * V[5];
    d[1] = ((-1.0 / 130.0) * sqrt26 * Y8_p5 + (-1.0 / 130.0) * sqrt910 * Y8_p7 + 0.2 * Y6_p5 * R2) * V[0] - 1.5 * Y6_p5 * V[1] + ((-8.0 / 65.0) * sqrt26 * Y8_p5 + (-8.0 / 65.0) * sqrt910 * Y8_p7 + 3.2 * Y6_p5 * R2) * V[3] - 21.0 * Y6_p5 * V[4];
    d[0] = ((1.0 / 130.0) * sqrt2 * Y8_p4 + (-1.0 / 65.0) * sqrt910 * Y8_p8 + (-1.0 / 55.0) * sqrt22 * Y6_p4 * R2 + (9.0 / 286.0) * sqrt110 * Y4_p4 * R4) * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_p4 + (-9.0 / 22.0) * sqrt110 * Y4_p4 * R2) * V[1] + 1.125 * sqrt110 * Y4_p4 * V[2] + ((8.0 / 65.0) * sqrt2 * Y8_p4 + (-16.0 / 65.0) * sqrt910 * Y8_p8 + (-16.0 / 55.0) * sqrt22 * Y6_p4 * R2 + (72.0 / 143.0) * sqrt110 * Y4_p4 * R4) * V[3] + ((21.0 / 11.0) * sqrt22 * Y6_p4 + (-63.0 / 11.0) * sqrt110 * Y4_p4 * R2) * V[4] + 13.5 * sqrt110 * Y4_p4 * V[5];
}

}  // namespace kinab