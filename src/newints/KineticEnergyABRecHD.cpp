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

#include "KineticEnergyABRecHD.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_h_d(
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
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt22 = std::sqrt(22.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt55 = std::sqrt(55.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt143 = std::sqrt(143.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt1430 = std::sqrt(1430.0);
    const auto sqrt2002 = std::sqrt(2002.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^2 · β^5 · p^{-7} · (s|T|s)
    //   V[1] ↔ α · β^4 · p^{-6} · (s|T|s)
    //   V[2] ↔ β^3 · p^{-5} · (s|T|s)
    //   V[3] ↔ α^2 · β^5 · p^{-7} · ξ · (s|S|s)
    //   V[4] ↔ α · β^4 · p^{-6} · ξ · (s|S|s)
    //   V[5] ↔ β^3 · p^{-5} · ξ · (s|S|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto pinv7 = pinv6 * pinv;
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha2 * beta5 * pinv7;
            V[1] += cab_tt * alpha * beta4 * pinv6;
            V[2] += cab_tt * beta3 * pinv5;
            V[3] += cab_ss * alpha2 * beta5 * pinv7 * xi;
            V[4] += cab_ss * alpha * beta4 * pinv6 * xi;
            V[5] += cab_ss * beta3 * pinv5 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 11 × 5 spherical block ----
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
    d[27] = ((-63.0 / 143.0) * Y7_p0 + (-10.0 / 39.0) * Y5_p0 * R2 + (-10.0 / 33.0) * Y3_p0 * R4) * V[0] + ((5.0 / 3.0) * Y5_p0 + (10.0 / 3.0) * Y3_p0 * R2) * V[1] - 7.5 * Y3_p0 * V[2] + ((-882.0 / 143.0) * Y7_p0 + (-140.0 / 39.0) * Y5_p0 * R2 + (-140.0 / 33.0) * Y3_p0 * R4) * V[3] + (20.0 * Y5_p0 + 40.0 * Y3_p0 * R2) * V[4] - 75.0 * Y3_p0 * V[5];
    d[28] = ((-12.0 / 143.0) * sqrt21 * Y7_p1 + (-1.0 / 39.0) * sqrt5 * Y5_p1 * R2 + (5.0 / 33.0) * sqrt2 * Y3_p1 * R4) * V[0] + ((1.0 / 6.0) * sqrt5 * Y5_p1 + (-5.0 / 3.0) * sqrt2 * Y3_p1 * R2) * V[1] + 3.75 * sqrt2 * Y3_p1 * V[2] + ((-168.0 / 143.0) * sqrt21 * Y7_p1 + (-14.0 / 39.0) * sqrt5 * Y5_p1 * R2 + (70.0 / 33.0) * sqrt2 * Y3_p1 * R4) * V[3] + (2.0 * sqrt5 * Y5_p1 - 20.0 * sqrt2 * Y3_p1 * R2) * V[4] + 37.5 * sqrt2 * Y3_p1 * V[5];
    d[29] = ((-9.0 / 143.0) * sqrt14 * Y7_p2 + (2.0 / 39.0) * sqrt35 * Y5_p2 * R2 + (-1.0 / 33.0) * sqrt5 * Y3_p2 * R4) * V[0] + ((-1.0 / 3.0) * sqrt35 * Y5_p2 + (1.0 / 3.0) * sqrt5 * Y3_p2 * R2) * V[1] - 0.75 * sqrt5 * Y3_p2 * V[2] + ((-126.0 / 143.0) * sqrt14 * Y7_p2 + (28.0 / 39.0) * sqrt35 * Y5_p2 * R2 + (-14.0 / 33.0) * sqrt5 * Y3_p2 * R4) * V[3] + (-4.0 * sqrt35 * Y5_p2 + 4.0 * sqrt5 * Y3_p2 * R2) * V[4] - 7.5 * sqrt5 * Y3_p2 * V[5];
    d[26] = ((-12.0 / 143.0) * sqrt21 * Y7_n1 + (-1.0 / 39.0) * sqrt5 * Y5_n1 * R2 + (5.0 / 33.0) * sqrt2 * Y3_n1 * R4) * V[0] + ((1.0 / 6.0) * sqrt5 * Y5_n1 + (-5.0 / 3.0) * sqrt2 * Y3_n1 * R2) * V[1] + 3.75 * sqrt2 * Y3_n1 * V[2] + ((-168.0 / 143.0) * sqrt21 * Y7_n1 + (-14.0 / 39.0) * sqrt5 * Y5_n1 * R2 + (70.0 / 33.0) * sqrt2 * Y3_n1 * R4) * V[3] + (2.0 * sqrt5 * Y5_n1 - 20.0 * sqrt2 * Y3_n1 * R2) * V[4] + 37.5 * sqrt2 * Y3_n1 * V[5];
    d[25] = ((-9.0 / 143.0) * sqrt14 * Y7_n2 + (2.0 / 39.0) * sqrt35 * Y5_n2 * R2 + (-1.0 / 33.0) * sqrt5 * Y3_n2 * R4) * V[0] + ((-1.0 / 3.0) * sqrt35 * Y5_n2 + (1.0 / 3.0) * sqrt5 * Y3_n2 * R2) * V[1] - 0.75 * sqrt5 * Y3_n2 * V[2] + ((-126.0 / 143.0) * sqrt14 * Y7_n2 + (28.0 / 39.0) * sqrt35 * Y5_n2 * R2 + (-14.0 / 33.0) * sqrt5 * Y3_n2 * R4) * V[3] + (-4.0 * sqrt35 * Y5_n2 + 4.0 * sqrt5 * Y3_n2 * R2) * V[4] - 7.5 * sqrt5 * Y3_n2 * V[5];
    d[32] = ((-6.0 / 143.0) * sqrt105 * Y7_p1 + (-3.0 / 13.0) * Y5_p1 * R2 + (-1.0 / 11.0) * sqrt10 * Y3_p1 * R4) * V[0] + (1.5 * Y5_p1 + sqrt10 * Y3_p1 * R2) * V[1] - 2.25 * sqrt10 * Y3_p1 * V[2] + ((-84.0 / 143.0) * sqrt105 * Y7_p1 + (-42.0 / 13.0) * Y5_p1 * R2 + (-14.0 / 11.0) * sqrt10 * Y3_p1 * R4) * V[3] + (18.0 * Y5_p1 + 12.0 * sqrt10 * Y3_p1 * R2) * V[4] - 22.5 * sqrt10 * Y3_p1 * V[5];
    d[33] = ((21.0 / 143.0) * sqrt5 * Y7_p0 + (-3.0 / 143.0) * sqrt210 * Y7_p2 + (-1.0 / 39.0) * sqrt5 * Y5_p0 * R2 + (-1.0 / 39.0) * sqrt21 * Y5_p2 * R2 + (-4.0 / 33.0) * sqrt5 * Y3_p0 * R4 + (2.0 / 33.0) * sqrt3 * Y3_p2 * R4) * V[0] + ((1.0 / 6.0) * sqrt5 * Y5_p0 + (1.0 / 6.0) * sqrt21 * Y5_p2 + (4.0 / 3.0) * sqrt5 * Y3_p0 * R2 + (-2.0 / 3.0) * sqrt3 * Y3_p2 * R2) * V[1] + (-3.0 * sqrt5 * Y3_p0 + 1.5 * sqrt3 * Y3_p2) * V[2] + ((294.0 / 143.0) * sqrt5 * Y7_p0 + (-42.0 / 143.0) * sqrt210 * Y7_p2 + (-14.0 / 39.0) * sqrt5 * Y5_p0 * R2 + (-14.0 / 39.0) * sqrt21 * Y5_p2 * R2 + (-56.0 / 33.0) * sqrt5 * Y3_p0 * R4 + (28.0 / 33.0) * sqrt3 * Y3_p2 * R4) * V[3] + (2.0 * sqrt5 * Y5_p0 + 2.0 * sqrt21 * Y5_p2 + 16.0 * sqrt5 * Y3_p0 * R2 - 8.0 * sqrt3 * Y3_p2 * R2) * V[4] + (-30.0 * sqrt5 * Y3_p0 + 15.0 * sqrt3 * Y3_p2) * V[5];
    d[34] = ((3.0 / 143.0) * sqrt35 * Y7_p1 + (-3.0 / 143.0) * sqrt105 * Y7_p3 + (-5.0 / 39.0) * sqrt3 * Y5_p1 * R2 + (2.0 / 39.0) * sqrt14 * Y5_p3 * R2 + (1.0 / 66.0) * sqrt30 * Y3_p1 * R4 + (-1.0 / 66.0) * sqrt2 * Y3_p3 * R4) * V[0] + ((5.0 / 6.0) * sqrt3 * Y5_p1 + (-1.0 / 3.0) * sqrt14 * Y5_p3 + (-1.0 / 6.0) * sqrt30 * Y3_p1 * R2 + (1.0 / 6.0) * sqrt2 * Y3_p3 * R2) * V[1] + (0.375 * sqrt30 * Y3_p1 - 0.375 * sqrt2 * Y3_p3) * V[2] + ((42.0 / 143.0) * sqrt35 * Y7_p1 + (-42.0 / 143.0) * sqrt105 * Y7_p3 + (-70.0 / 39.0) * sqrt3 * Y5_p1 * R2 + (28.0 / 39.0) * sqrt14 * Y5_p3 * R2 + (7.0 / 33.0) * sqrt30 * Y3_p1 * R4 + (-7.0 / 33.0) * sqrt2 * Y3_p3 * R4) * V[3] + (10.0 * sqrt3 * Y5_p1 - 4.0 * sqrt14 * Y5_p3 - 2.0 * sqrt30 * Y3_p1 * R2 + 2.0 * sqrt2 * Y3_p3 * R2) * V[4] + (3.75 * sqrt30 * Y3_p1 - 3.75 * sqrt2 * Y3_p3) * V[5];
    d[31] = ((-3.0 / 143.0) * sqrt210 * Y7_n2 + (-1.0 / 39.0) * sqrt21 * Y5_n2 * R2 + (2.0 / 33.0) * sqrt3 * Y3_n2 * R4) * V[0] + ((1.0 / 6.0) * sqrt21 * Y5_n2 + (-2.0 / 3.0) * sqrt3 * Y3_n2 * R2) * V[1] + 1.5 * sqrt3 * Y3_n2 * V[2] + ((-42.0 / 143.0) * sqrt210 * Y7_n2 + (-14.0 / 39.0) * sqrt21 * Y5_n2 * R2 + (28.0 / 33.0) * sqrt3 * Y3_n2 * R4) * V[3] + (2.0 * sqrt21 * Y5_n2 - 8.0 * sqrt3 * Y3_n2 * R2) * V[4] + 15.0 * sqrt3 * Y3_n2 * V[5];
    d[30] = ((-3.0 / 143.0) * sqrt105 * Y7_n3 + (3.0 / 143.0) * sqrt35 * Y7_n1 + (2.0 / 39.0) * sqrt14 * Y5_n3 * R2 + (-5.0 / 39.0) * sqrt3 * Y5_n1 * R2 + (-1.0 / 66.0) * sqrt2 * Y3_n3 * R4 + (1.0 / 66.0) * sqrt30 * Y3_n1 * R4) * V[0] + ((-1.0 / 3.0) * sqrt14 * Y5_n3 + (5.0 / 6.0) * sqrt3 * Y5_n1 + (1.0 / 6.0) * sqrt2 * Y3_n3 * R2 + (-1.0 / 6.0) * sqrt30 * Y3_n1 * R2) * V[1] + (-0.375 * sqrt2 * Y3_n3 + 0.375 * sqrt30 * Y3_n1) * V[2] + ((-42.0 / 143.0) * sqrt105 * Y7_n3 + (42.0 / 143.0) * sqrt35 * Y7_n1 + (28.0 / 39.0) * sqrt14 * Y5_n3 * R2 + (-70.0 / 39.0) * sqrt3 * Y5_n1 * R2 + (-7.0 / 33.0) * sqrt2 * Y3_n3 * R4 + (7.0 / 33.0) * sqrt30 * Y3_n1 * R4) * V[3] + (-4.0 * sqrt14 * Y5_n3 + 10.0 * sqrt3 * Y5_n1 + 2.0 * sqrt2 * Y3_n3 * R2 - 2.0 * sqrt30 * Y3_n1 * R2) * V[4] + (-3.75 * sqrt2 * Y3_n3 + 3.75 * sqrt30 * Y3_n1) * V[5];
    d[37] = ((-18.0 / 143.0) * sqrt10 * Y7_p2 + (-2.0 / 13.0) * Y5_p2 * R2 + (-1.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[0] + (Y5_p2 + sqrt7 * Y3_p2 * R2) * V[1] - 2.25 * sqrt7 * Y3_p2 * V[2] + ((-252.0 / 143.0) * sqrt10 * Y7_p2 + (-28.0 / 13.0) * Y5_p2 * R2 + (-14.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[3] + (12.0 * Y5_p2 + 12.0 * sqrt7 * Y3_p2 * R2) * V[4] - 22.5 * sqrt7 * Y3_p2 * V[5];
    d[38] = ((12.0 / 143.0) * sqrt5 * Y7_p1 + (-12.0 / 143.0) * sqrt15 * Y7_p3 + (-1.0 / 39.0) * sqrt21 * Y5_p1 * R2 + (-5.0 / 39.0) * sqrt2 * Y5_p3 * R2 + (-1.0 / 66.0) * sqrt210 * Y3_p1 * R4 + (1.0 / 66.0) * sqrt14 * Y3_p3 * R4) * V[0] + ((1.0 / 6.0) * sqrt21 * Y5_p1 + (5.0 / 6.0) * sqrt2 * Y5_p3 + (1.0 / 6.0) * sqrt210 * Y3_p1 * R2 + (-1.0 / 6.0) * sqrt14 * Y3_p3 * R2) * V[1] + (-0.375 * sqrt210 * Y3_p1 + 0.375 * sqrt14 * Y3_p3) * V[2] + ((168.0 / 143.0) * sqrt5 * Y7_p1 + (-168.0 / 143.0) * sqrt15 * Y7_p3 + (-14.0 / 39.0) * sqrt21 * Y5_p1 * R2 + (-70.0 / 39.0) * sqrt2 * Y5_p3 * R2 + (-7.0 / 33.0) * sqrt210 * Y3_p1 * R4 + (7.0 / 33.0) * sqrt14 * Y3_p3 * R4) * V[3] + (2.0 * sqrt21 * Y5_p1 + 10.0 * sqrt2 * Y5_p3 + 2.0 * sqrt210 * Y3_p1 * R2 - 2.0 * sqrt14 * Y3_p3 * R2) * V[4] + (-3.75 * sqrt210 * Y3_p1 + 3.75 * sqrt14 * Y3_p3) * V[5];
    d[39] = ((-3.0 / 143.0) * sqrt35 * Y7_p0 + (-3.0 / 143.0) * sqrt165 * Y7_p4 + (2.0 / 39.0) * sqrt35 * Y5_p0 * R2 + (2.0 / 13.0) * Y5_p4 * R2 + (-1.0 / 33.0) * sqrt35 * Y3_p0 * R4) * V[0] + ((-1.0 / 3.0) * sqrt35 * Y5_p0 - Y5_p4 + (1.0 / 3.0) * sqrt35 * Y3_p0 * R2) * V[1] - 0.75 * sqrt35 * Y3_p0 * V[2] + ((-42.0 / 143.0) * sqrt35 * Y7_p0 + (-42.0 / 143.0) * sqrt165 * Y7_p4 + (28.0 / 39.0) * sqrt35 * Y5_p0 * R2 + (28.0 / 13.0) * Y5_p4 * R2 + (-14.0 / 33.0) * sqrt35 * Y3_p0 * R4) * V[3] + (-4.0 * sqrt35 * Y5_p0 - 12.0 * Y5_p4 + 4.0 * sqrt35 * Y3_p0 * R2) * V[4] - 7.5 * sqrt35 * Y3_p0 * V[5];
    d[36] = ((-12.0 / 143.0) * sqrt15 * Y7_n3 + (-12.0 / 143.0) * sqrt5 * Y7_n1 + (-5.0 / 39.0) * sqrt2 * Y5_n3 * R2 + (1.0 / 39.0) * sqrt21 * Y5_n1 * R2 + (1.0 / 66.0) * sqrt14 * Y3_n3 * R4 + (1.0 / 66.0) * sqrt210 * Y3_n1 * R4) * V[0] + ((5.0 / 6.0) * sqrt2 * Y5_n3 + (-1.0 / 6.0) * sqrt21 * Y5_n1 + (-1.0 / 6.0) * sqrt14 * Y3_n3 * R2 + (-1.0 / 6.0) * sqrt210 * Y3_n1 * R2) * V[1] + (0.375 * sqrt14 * Y3_n3 + 0.375 * sqrt210 * Y3_n1) * V[2] + ((-168.0 / 143.0) * sqrt15 * Y7_n3 + (-168.0 / 143.0) * sqrt5 * Y7_n1 + (-70.0 / 39.0) * sqrt2 * Y5_n3 * R2 + (14.0 / 39.0) * sqrt21 * Y5_n1 * R2 + (7.0 / 33.0) * sqrt14 * Y3_n3 * R4 + (7.0 / 33.0) * sqrt210 * Y3_n1 * R4) * V[3] + (10.0 * sqrt2 * Y5_n3 - 2.0 * sqrt21 * Y5_n1 - 2.0 * sqrt14 * Y3_n3 * R2 - 2.0 * sqrt210 * Y3_n1 * R2) * V[4] + (3.75 * sqrt14 * Y3_n3 + 3.75 * sqrt210 * Y3_n1) * V[5];
    d[35] = ((-3.0 / 143.0) * sqrt165 * Y7_n4 + (2.0 / 13.0) * Y5_n4 * R2) * V[0] - Y5_n4 * V[1] + ((-42.0 / 143.0) * sqrt165 * Y7_n4 + (28.0 / 13.0) * Y5_n4 * R2) * V[3] - 12.0 * Y5_n4 * V[4];
    d[42] = ((-9.0 / 143.0) * sqrt30 * Y7_p3 + (-1.0 / 39.0) * Y5_p3 * R2 + (-2.0 / 33.0) * sqrt7 * Y3_p3 * R4) * V[0] + ((1.0 / 6.0) * Y5_p3 + (2.0 / 3.0) * sqrt7 * Y3_p3 * R2) * V[1] - 1.5 * sqrt7 * Y3_p3 * V[2] + ((-126.0 / 143.0) * sqrt30 * Y7_p3 + (-14.0 / 39.0) * Y5_p3 * R2 + (-28.0 / 33.0) * sqrt7 * Y3_p3 * R4) * V[3] + (2.0 * Y5_p3 + 8.0 * sqrt7 * Y3_p3 * R2) * V[4] - 15.0 * sqrt7 * Y3_p3 * V[5];
    d[43] = ((9.0 / 143.0) * sqrt5 * Y7_p2 + (-9.0 / 286.0) * sqrt110 * Y7_p4 + (-5.0 / 39.0) * sqrt2 * Y5_p2 * R2 + (-7.0 / 78.0) * sqrt6 * Y5_p4 * R2 + (-2.0 / 33.0) * sqrt14 * Y3_p2 * R4) * V[0] + ((5.0 / 6.0) * sqrt2 * Y5_p2 + (7.0 / 12.0) * sqrt6 * Y5_p4 + (2.0 / 3.0) * sqrt14 * Y3_p2 * R2) * V[1] - 1.5 * sqrt14 * Y3_p2 * V[2] + ((126.0 / 143.0) * sqrt5 * Y7_p2 + (-63.0 / 143.0) * sqrt110 * Y7_p4 + (-70.0 / 39.0) * sqrt2 * Y5_p2 * R2 + (-49.0 / 39.0) * sqrt6 * Y5_p4 * R2 + (-28.0 / 33.0) * sqrt14 * Y3_p2 * R4) * V[3] + (10.0 * sqrt2 * Y5_p2 + 7.0 * sqrt6 * Y5_p4 + 8.0 * sqrt14 * Y3_p2 * R2) * V[4] - 15.0 * sqrt14 * Y3_p2 * V[5];
    d[44] = ((-3.0 / 286.0) * sqrt30 * Y7_p1 + (-9.0 / 286.0) * sqrt110 * Y7_p5 + (2.0 / 39.0) * sqrt14 * Y5_p1 * R2 + (1.0 / 39.0) * sqrt15 * Y5_p5 * R2 + (-1.0 / 33.0) * sqrt35 * Y3_p1 * R4) * V[0] + ((-1.0 / 3.0) * sqrt14 * Y5_p1 + (-1.0 / 6.0) * sqrt15 * Y5_p5 + (1.0 / 3.0) * sqrt35 * Y3_p1 * R2) * V[1] - 0.75 * sqrt35 * Y3_p1 * V[2] + ((-21.0 / 143.0) * sqrt30 * Y7_p1 + (-63.0 / 143.0) * sqrt110 * Y7_p5 + (28.0 / 39.0) * sqrt14 * Y5_p1 * R2 + (14.0 / 39.0) * sqrt15 * Y5_p5 * R2 + (-14.0 / 33.0) * sqrt35 * Y3_p1 * R4) * V[3] + (-4.0 * sqrt14 * Y5_p1 - 2.0 * sqrt15 * Y5_p5 + 4.0 * sqrt35 * Y3_p1 * R2) * V[4] - 7.5 * sqrt35 * Y3_p1 * V[5];
    d[41] = ((-9.0 / 286.0) * sqrt110 * Y7_n4 + (-9.0 / 143.0) * sqrt5 * Y7_n2 + (-7.0 / 78.0) * sqrt6 * Y5_n4 * R2 + (5.0 / 39.0) * sqrt2 * Y5_n2 * R2 + (2.0 / 33.0) * sqrt14 * Y3_n2 * R4) * V[0] + ((7.0 / 12.0) * sqrt6 * Y5_n4 + (-5.0 / 6.0) * sqrt2 * Y5_n2 + (-2.0 / 3.0) * sqrt14 * Y3_n2 * R2) * V[1] + 1.5 * sqrt14 * Y3_n2 * V[2] + ((-63.0 / 143.0) * sqrt110 * Y7_n4 + (-126.0 / 143.0) * sqrt5 * Y7_n2 + (-49.0 / 39.0) * sqrt6 * Y5_n4 * R2 + (70.0 / 39.0) * sqrt2 * Y5_n2 * R2 + (28.0 / 33.0) * sqrt14 * Y3_n2 * R4) * V[3] + (7.0 * sqrt6 * Y5_n4 - 10.0 * sqrt2 * Y5_n2 - 8.0 * sqrt14 * Y3_n2 * R2) * V[4] + 15.0 * sqrt14 * Y3_n2 * V[5];
    d[40] = ((-9.0 / 286.0) * sqrt110 * Y7_n5 + (3.0 / 286.0) * sqrt30 * Y7_n1 + (1.0 / 39.0) * sqrt15 * Y5_n5 * R2 + (-2.0 / 39.0) * sqrt14 * Y5_n1 * R2 + (1.0 / 33.0) * sqrt35 * Y3_n1 * R4) * V[0] + ((-1.0 / 6.0) * sqrt15 * Y5_n5 + (1.0 / 3.0) * sqrt14 * Y5_n1 + (-1.0 / 3.0) * sqrt35 * Y3_n1 * R2) * V[1] + 0.75 * sqrt35 * Y3_n1 * V[2] + ((-63.0 / 143.0) * sqrt110 * Y7_n5 + (21.0 / 143.0) * sqrt30 * Y7_n1 + (14.0 / 39.0) * sqrt15 * Y5_n5 * R2 + (-28.0 / 39.0) * sqrt14 * Y5_n1 * R2 + (14.0 / 33.0) * sqrt35 * Y3_n1 * R4) * V[3] + (-2.0 * sqrt15 * Y5_n5 + 4.0 * sqrt14 * Y5_n1 - 4.0 * sqrt35 * Y3_n1 * R2) * V[4] + 7.5 * sqrt35 * Y3_n1 * V[5];
    d[47] = ((-3.0 / 143.0) * sqrt165 * Y7_p4 + (2.0 / 13.0) * Y5_p4 * R2) * V[0] - Y5_p4 * V[1] + ((-42.0 / 143.0) * sqrt165 * Y7_p4 + (28.0 / 13.0) * Y5_p4 * R2) * V[3] - 12.0 * Y5_p4 * V[4];
    d[48] = ((6.0 / 143.0) * sqrt5 * Y7_p3 + (-6.0 / 143.0) * sqrt55 * Y7_p5 + (-7.0 / 78.0) * sqrt6 * Y5_p3 * R2 + (-1.0 / 26.0) * sqrt30 * Y5_p5 * R2 + (-1.0 / 33.0) * sqrt42 * Y3_p3 * R4) * V[0] + ((7.0 / 12.0) * sqrt6 * Y5_p3 + 0.25 * sqrt30 * Y5_p5 + (1.0 / 3.0) * sqrt42 * Y3_p3 * R2) * V[1] - 0.75 * sqrt42 * Y3_p3 * V[2] + ((84.0 / 143.0) * sqrt5 * Y7_p3 + (-84.0 / 143.0) * sqrt55 * Y7_p5 + (-49.0 / 39.0) * sqrt6 * Y5_p3 * R2 + (-7.0 / 13.0) * sqrt30 * Y5_p5 * R2 + (-14.0 / 33.0) * sqrt42 * Y3_p3 * R4) * V[3] + (7.0 * sqrt6 * Y5_p3 + 3.0 * sqrt30 * Y5_p5 + 4.0 * sqrt42 * Y3_p3 * R2) * V[4] - 7.5 * sqrt42 * Y3_p3 * V[5];
    d[49] = ((-3.0 / 286.0) * sqrt10 * Y7_p2 + (-3.0 / 286.0) * sqrt1430 * Y7_p6 + (2.0 / 13.0) * Y5_p2 * R2 + (-1.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[0] + (-Y5_p2 + sqrt7 * Y3_p2 * R2) * V[1] - 2.25 * sqrt7 * Y3_p2 * V[2] + ((-21.0 / 143.0) * sqrt10 * Y7_p2 + (-21.0 / 143.0) * sqrt1430 * Y7_p6 + (28.0 / 13.0) * Y5_p2 * R2 + (-14.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[3] + (-12.0 * Y5_p2 + 12.0 * sqrt7 * Y3_p2 * R2) * V[4] - 22.5 * sqrt7 * Y3_p2 * V[5];
    d[46] = ((-6.0 / 143.0) * sqrt55 * Y7_n5 + (-6.0 / 143.0) * sqrt5 * Y7_n3 + (-1.0 / 26.0) * sqrt30 * Y5_n5 * R2 + (7.0 / 78.0) * sqrt6 * Y5_n3 * R2 + (1.0 / 33.0) * sqrt42 * Y3_n3 * R4) * V[0] + (0.25 * sqrt30 * Y5_n5 + (-7.0 / 12.0) * sqrt6 * Y5_n3 + (-1.0 / 3.0) * sqrt42 * Y3_n3 * R2) * V[1] + 0.75 * sqrt42 * Y3_n3 * V[2] + ((-84.0 / 143.0) * sqrt55 * Y7_n5 + (-84.0 / 143.0) * sqrt5 * Y7_n3 + (-7.0 / 13.0) * sqrt30 * Y5_n5 * R2 + (49.0 / 39.0) * sqrt6 * Y5_n3 * R2 + (14.0 / 33.0) * sqrt42 * Y3_n3 * R4) * V[3] + (3.0 * sqrt30 * Y5_n5 - 7.0 * sqrt6 * Y5_n3 - 4.0 * sqrt42 * Y3_n3 * R2) * V[4] + 7.5 * sqrt42 * Y3_n3 * V[5];
    d[45] = ((-3.0 / 286.0) * sqrt1430 * Y7_n6 + (3.0 / 286.0) * sqrt10 * Y7_n2 + (-2.0 / 13.0) * Y5_n2 * R2 + (1.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[0] + (Y5_n2 - sqrt7 * Y3_n2 * R2) * V[1] + 2.25 * sqrt7 * Y3_n2 * V[2] + ((-21.0 / 143.0) * sqrt1430 * Y7_n6 + (21.0 / 143.0) * sqrt10 * Y7_n2 + (-28.0 / 13.0) * Y5_n2 * R2 + (14.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[3] + (12.0 * Y5_n2 - 12.0 * sqrt7 * Y3_n2 * R2) * V[4] + 22.5 * sqrt7 * Y3_n2 * V[5];
    d[52] = ((-3.0 / 143.0) * sqrt66 * Y7_p5 + (5.0 / 13.0) * Y5_p5 * R2) * V[0] - 2.5 * Y5_p5 * V[1] + ((-42.0 / 143.0) * sqrt66 * Y7_p5 + (70.0 / 13.0) * Y5_p5 * R2) * V[3] - 30.0 * Y5_p5 * V[4];
    d[53] = ((3.0 / 286.0) * sqrt22 * Y7_p4 + (-3.0 / 143.0) * sqrt143 * Y7_p6 + (-1.0 / 26.0) * sqrt30 * Y5_p4 * R2) * V[0] + 0.25 * sqrt30 * Y5_p4 * V[1] + ((21.0 / 143.0) * sqrt22 * Y7_p4 + (-42.0 / 143.0) * sqrt143 * Y7_p6 + (-7.0 / 13.0) * sqrt30 * Y5_p4 * R2) * V[3] + 3.0 * sqrt30 * Y5_p4 * V[4];
    d[54] = ((-3.0 / 286.0) * sqrt2 * Y7_p3 + (-3.0 / 286.0) * sqrt2002 * Y7_p7 + (1.0 / 39.0) * sqrt15 * Y5_p3 * R2 + (-1.0 / 33.0) * sqrt105 * Y3_p3 * R4) * V[0] + ((-1.0 / 6.0) * sqrt15 * Y5_p3 + (1.0 / 3.0) * sqrt105 * Y3_p3 * R2) * V[1] - 0.75 * sqrt105 * Y3_p3 * V[2] + ((-21.0 / 143.0) * sqrt2 * Y7_p3 + (-21.0 / 143.0) * sqrt2002 * Y7_p7 + (14.0 / 39.0) * sqrt15 * Y5_p3 * R2 + (-14.0 / 33.0) * sqrt105 * Y3_p3 * R4) * V[3] + (-2.0 * sqrt15 * Y5_p3 + 4.0 * sqrt105 * Y3_p3 * R2) * V[4] - 7.5 * sqrt105 * Y3_p3 * V[5];
    d[51] = ((-3.0 / 143.0) * sqrt143 * Y7_n6 + (-3.0 / 286.0) * sqrt22 * Y7_n4 + (1.0 / 26.0) * sqrt30 * Y5_n4 * R2) * V[0] - 0.25 * sqrt30 * Y5_n4 * V[1] + ((-42.0 / 143.0) * sqrt143 * Y7_n6 + (-21.0 / 143.0) * sqrt22 * Y7_n4 + (7.0 / 13.0) * sqrt30 * Y5_n4 * R2) * V[3] - 3.0 * sqrt30 * Y5_n4 * V[4];
    d[50] = ((-3.0 / 286.0) * sqrt2002 * Y7_n7 + (3.0 / 286.0) * sqrt2 * Y7_n3 + (-1.0 / 39.0) * sqrt15 * Y5_n3 * R2 + (1.0 / 33.0) * sqrt105 * Y3_n3 * R4) * V[0] + ((1.0 / 6.0) * sqrt15 * Y5_n3 + (-1.0 / 3.0) * sqrt105 * Y3_n3 * R2) * V[1] + 0.75 * sqrt105 * Y3_n3 * V[2] + ((-21.0 / 143.0) * sqrt2002 * Y7_n7 + (21.0 / 143.0) * sqrt2 * Y7_n3 + (-14.0 / 39.0) * sqrt15 * Y5_n3 * R2 + (14.0 / 33.0) * sqrt105 * Y3_n3 * R4) * V[3] + (2.0 * sqrt15 * Y5_n3 - 4.0 * sqrt105 * Y3_n3 * R2) * V[4] + 7.5 * sqrt105 * Y3_n3 * V[5];
    d[22] = ((-6.0 / 143.0) * sqrt105 * Y7_n1 + (-3.0 / 13.0) * Y5_n1 * R2 + (-1.0 / 11.0) * sqrt10 * Y3_n1 * R4) * V[0] + (1.5 * Y5_n1 + sqrt10 * Y3_n1 * R2) * V[1] - 2.25 * sqrt10 * Y3_n1 * V[2] + ((-84.0 / 143.0) * sqrt105 * Y7_n1 + (-42.0 / 13.0) * Y5_n1 * R2 + (-14.0 / 11.0) * sqrt10 * Y3_n1 * R4) * V[3] + (18.0 * Y5_n1 + 12.0 * sqrt10 * Y3_n1 * R2) * V[4] - 22.5 * sqrt10 * Y3_n1 * V[5];
    d[23] = ((-3.0 / 143.0) * sqrt210 * Y7_n2 + (-1.0 / 39.0) * sqrt21 * Y5_n2 * R2 + (2.0 / 33.0) * sqrt3 * Y3_n2 * R4) * V[0] + ((1.0 / 6.0) * sqrt21 * Y5_n2 + (-2.0 / 3.0) * sqrt3 * Y3_n2 * R2) * V[1] + 1.5 * sqrt3 * Y3_n2 * V[2] + ((-42.0 / 143.0) * sqrt210 * Y7_n2 + (-14.0 / 39.0) * sqrt21 * Y5_n2 * R2 + (28.0 / 33.0) * sqrt3 * Y3_n2 * R4) * V[3] + (2.0 * sqrt21 * Y5_n2 - 8.0 * sqrt3 * Y3_n2 * R2) * V[4] + 15.0 * sqrt3 * Y3_n2 * V[5];
    d[24] = ((-3.0 / 143.0) * sqrt105 * Y7_n3 + (-3.0 / 143.0) * sqrt35 * Y7_n1 + (2.0 / 39.0) * sqrt14 * Y5_n3 * R2 + (5.0 / 39.0) * sqrt3 * Y5_n1 * R2 + (-1.0 / 66.0) * sqrt2 * Y3_n3 * R4 + (-1.0 / 66.0) * sqrt30 * Y3_n1 * R4) * V[0] + ((-1.0 / 3.0) * sqrt14 * Y5_n3 + (-5.0 / 6.0) * sqrt3 * Y5_n1 + (1.0 / 6.0) * sqrt2 * Y3_n3 * R2 + (1.0 / 6.0) * sqrt30 * Y3_n1 * R2) * V[1] + (-0.375 * sqrt2 * Y3_n3 - 0.375 * sqrt30 * Y3_n1) * V[2] + ((-42.0 / 143.0) * sqrt105 * Y7_n3 + (-42.0 / 143.0) * sqrt35 * Y7_n1 + (28.0 / 39.0) * sqrt14 * Y5_n3 * R2 + (70.0 / 39.0) * sqrt3 * Y5_n1 * R2 + (-7.0 / 33.0) * sqrt2 * Y3_n3 * R4 + (-7.0 / 33.0) * sqrt30 * Y3_n1 * R4) * V[3] + (-4.0 * sqrt14 * Y5_n3 - 10.0 * sqrt3 * Y5_n1 + 2.0 * sqrt2 * Y3_n3 * R2 + 2.0 * sqrt30 * Y3_n1 * R2) * V[4] + (-3.75 * sqrt2 * Y3_n3 - 3.75 * sqrt30 * Y3_n1) * V[5];
    d[21] = ((21.0 / 143.0) * sqrt5 * Y7_p0 + (3.0 / 143.0) * sqrt210 * Y7_p2 + (-1.0 / 39.0) * sqrt5 * Y5_p0 * R2 + (1.0 / 39.0) * sqrt21 * Y5_p2 * R2 + (-4.0 / 33.0) * sqrt5 * Y3_p0 * R4 + (-2.0 / 33.0) * sqrt3 * Y3_p2 * R4) * V[0] + ((1.0 / 6.0) * sqrt5 * Y5_p0 + (-1.0 / 6.0) * sqrt21 * Y5_p2 + (4.0 / 3.0) * sqrt5 * Y3_p0 * R2 + (2.0 / 3.0) * sqrt3 * Y3_p2 * R2) * V[1] + (-3.0 * sqrt5 * Y3_p0 - 1.5 * sqrt3 * Y3_p2) * V[2] + ((294.0 / 143.0) * sqrt5 * Y7_p0 + (42.0 / 143.0) * sqrt210 * Y7_p2 + (-14.0 / 39.0) * sqrt5 * Y5_p0 * R2 + (14.0 / 39.0) * sqrt21 * Y5_p2 * R2 + (-56.0 / 33.0) * sqrt5 * Y3_p0 * R4 + (-28.0 / 33.0) * sqrt3 * Y3_p2 * R4) * V[3] + (2.0 * sqrt5 * Y5_p0 - 2.0 * sqrt21 * Y5_p2 + 16.0 * sqrt5 * Y3_p0 * R2 + 8.0 * sqrt3 * Y3_p2 * R2) * V[4] + (-30.0 * sqrt5 * Y3_p0 - 15.0 * sqrt3 * Y3_p2) * V[5];
    d[20] = ((3.0 / 143.0) * sqrt35 * Y7_p1 + (3.0 / 143.0) * sqrt105 * Y7_p3 + (-5.0 / 39.0) * sqrt3 * Y5_p1 * R2 + (-2.0 / 39.0) * sqrt14 * Y5_p3 * R2 + (1.0 / 66.0) * sqrt30 * Y3_p1 * R4 + (1.0 / 66.0) * sqrt2 * Y3_p3 * R4) * V[0] + ((5.0 / 6.0) * sqrt3 * Y5_p1 + (1.0 / 3.0) * sqrt14 * Y5_p3 + (-1.0 / 6.0) * sqrt30 * Y3_p1 * R2 + (-1.0 / 6.0) * sqrt2 * Y3_p3 * R2) * V[1] + (0.375 * sqrt30 * Y3_p1 + 0.375 * sqrt2 * Y3_p3) * V[2] + ((42.0 / 143.0) * sqrt35 * Y7_p1 + (42.0 / 143.0) * sqrt105 * Y7_p3 + (-70.0 / 39.0) * sqrt3 * Y5_p1 * R2 + (-28.0 / 39.0) * sqrt14 * Y5_p3 * R2 + (7.0 / 33.0) * sqrt30 * Y3_p1 * R4 + (7.0 / 33.0) * sqrt2 * Y3_p3 * R4) * V[3] + (10.0 * sqrt3 * Y5_p1 + 4.0 * sqrt14 * Y5_p3 - 2.0 * sqrt30 * Y3_p1 * R2 - 2.0 * sqrt2 * Y3_p3 * R2) * V[4] + (3.75 * sqrt30 * Y3_p1 + 3.75 * sqrt2 * Y3_p3) * V[5];
    d[17] = ((-18.0 / 143.0) * sqrt10 * Y7_n2 + (-2.0 / 13.0) * Y5_n2 * R2 + (-1.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[0] + (Y5_n2 + sqrt7 * Y3_n2 * R2) * V[1] - 2.25 * sqrt7 * Y3_n2 * V[2] + ((-252.0 / 143.0) * sqrt10 * Y7_n2 + (-28.0 / 13.0) * Y5_n2 * R2 + (-14.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[3] + (12.0 * Y5_n2 + 12.0 * sqrt7 * Y3_n2 * R2) * V[4] - 22.5 * sqrt7 * Y3_n2 * V[5];
    d[18] = ((-12.0 / 143.0) * sqrt15 * Y7_n3 + (12.0 / 143.0) * sqrt5 * Y7_n1 + (-5.0 / 39.0) * sqrt2 * Y5_n3 * R2 + (-1.0 / 39.0) * sqrt21 * Y5_n1 * R2 + (1.0 / 66.0) * sqrt14 * Y3_n3 * R4 + (-1.0 / 66.0) * sqrt210 * Y3_n1 * R4) * V[0] + ((5.0 / 6.0) * sqrt2 * Y5_n3 + (1.0 / 6.0) * sqrt21 * Y5_n1 + (-1.0 / 6.0) * sqrt14 * Y3_n3 * R2 + (1.0 / 6.0) * sqrt210 * Y3_n1 * R2) * V[1] + (0.375 * sqrt14 * Y3_n3 - 0.375 * sqrt210 * Y3_n1) * V[2] + ((-168.0 / 143.0) * sqrt15 * Y7_n3 + (168.0 / 143.0) * sqrt5 * Y7_n1 + (-70.0 / 39.0) * sqrt2 * Y5_n3 * R2 + (-14.0 / 39.0) * sqrt21 * Y5_n1 * R2 + (7.0 / 33.0) * sqrt14 * Y3_n3 * R4 + (-7.0 / 33.0) * sqrt210 * Y3_n1 * R4) * V[3] + (10.0 * sqrt2 * Y5_n3 + 2.0 * sqrt21 * Y5_n1 - 2.0 * sqrt14 * Y3_n3 * R2 + 2.0 * sqrt210 * Y3_n1 * R2) * V[4] + (3.75 * sqrt14 * Y3_n3 - 3.75 * sqrt210 * Y3_n1) * V[5];
    d[19] = ((-3.0 / 143.0) * sqrt165 * Y7_n4 + (2.0 / 13.0) * Y5_n4 * R2) * V[0] - Y5_n4 * V[1] + ((-42.0 / 143.0) * sqrt165 * Y7_n4 + (28.0 / 13.0) * Y5_n4 * R2) * V[3] - 12.0 * Y5_n4 * V[4];
    d[16] = ((12.0 / 143.0) * sqrt5 * Y7_p1 + (12.0 / 143.0) * sqrt15 * Y7_p3 + (-1.0 / 39.0) * sqrt21 * Y5_p1 * R2 + (5.0 / 39.0) * sqrt2 * Y5_p3 * R2 + (-1.0 / 66.0) * sqrt210 * Y3_p1 * R4 + (-1.0 / 66.0) * sqrt14 * Y3_p3 * R4) * V[0] + ((1.0 / 6.0) * sqrt21 * Y5_p1 + (-5.0 / 6.0) * sqrt2 * Y5_p3 + (1.0 / 6.0) * sqrt210 * Y3_p1 * R2 + (1.0 / 6.0) * sqrt14 * Y3_p3 * R2) * V[1] + (-0.375 * sqrt210 * Y3_p1 - 0.375 * sqrt14 * Y3_p3) * V[2] + ((168.0 / 143.0) * sqrt5 * Y7_p1 + (168.0 / 143.0) * sqrt15 * Y7_p3 + (-14.0 / 39.0) * sqrt21 * Y5_p1 * R2 + (70.0 / 39.0) * sqrt2 * Y5_p3 * R2 + (-7.0 / 33.0) * sqrt210 * Y3_p1 * R4 + (-7.0 / 33.0) * sqrt14 * Y3_p3 * R4) * V[3] + (2.0 * sqrt21 * Y5_p1 - 10.0 * sqrt2 * Y5_p3 + 2.0 * sqrt210 * Y3_p1 * R2 + 2.0 * sqrt14 * Y3_p3 * R2) * V[4] + (-3.75 * sqrt210 * Y3_p1 - 3.75 * sqrt14 * Y3_p3) * V[5];
    d[15] = ((-3.0 / 143.0) * sqrt35 * Y7_p0 + (3.0 / 143.0) * sqrt165 * Y7_p4 + (2.0 / 39.0) * sqrt35 * Y5_p0 * R2 + (-2.0 / 13.0) * Y5_p4 * R2 + (-1.0 / 33.0) * sqrt35 * Y3_p0 * R4) * V[0] + ((-1.0 / 3.0) * sqrt35 * Y5_p0 + Y5_p4 + (1.0 / 3.0) * sqrt35 * Y3_p0 * R2) * V[1] - 0.75 * sqrt35 * Y3_p0 * V[2] + ((-42.0 / 143.0) * sqrt35 * Y7_p0 + (42.0 / 143.0) * sqrt165 * Y7_p4 + (28.0 / 39.0) * sqrt35 * Y5_p0 * R2 + (-28.0 / 13.0) * Y5_p4 * R2 + (-14.0 / 33.0) * sqrt35 * Y3_p0 * R4) * V[3] + (-4.0 * sqrt35 * Y5_p0 + 12.0 * Y5_p4 + 4.0 * sqrt35 * Y3_p0 * R2) * V[4] - 7.5 * sqrt35 * Y3_p0 * V[5];
    d[12] = ((-9.0 / 143.0) * sqrt30 * Y7_n3 + (-1.0 / 39.0) * Y5_n3 * R2 + (-2.0 / 33.0) * sqrt7 * Y3_n3 * R4) * V[0] + ((1.0 / 6.0) * Y5_n3 + (2.0 / 3.0) * sqrt7 * Y3_n3 * R2) * V[1] - 1.5 * sqrt7 * Y3_n3 * V[2] + ((-126.0 / 143.0) * sqrt30 * Y7_n3 + (-14.0 / 39.0) * Y5_n3 * R2 + (-28.0 / 33.0) * sqrt7 * Y3_n3 * R4) * V[3] + (2.0 * Y5_n3 + 8.0 * sqrt7 * Y3_n3 * R2) * V[4] - 15.0 * sqrt7 * Y3_n3 * V[5];
    d[13] = ((-9.0 / 286.0) * sqrt110 * Y7_n4 + (9.0 / 143.0) * sqrt5 * Y7_n2 + (-7.0 / 78.0) * sqrt6 * Y5_n4 * R2 + (-5.0 / 39.0) * sqrt2 * Y5_n2 * R2 + (-2.0 / 33.0) * sqrt14 * Y3_n2 * R4) * V[0] + ((7.0 / 12.0) * sqrt6 * Y5_n4 + (5.0 / 6.0) * sqrt2 * Y5_n2 + (2.0 / 3.0) * sqrt14 * Y3_n2 * R2) * V[1] - 1.5 * sqrt14 * Y3_n2 * V[2] + ((-63.0 / 143.0) * sqrt110 * Y7_n4 + (126.0 / 143.0) * sqrt5 * Y7_n2 + (-49.0 / 39.0) * sqrt6 * Y5_n4 * R2 + (-70.0 / 39.0) * sqrt2 * Y5_n2 * R2 + (-28.0 / 33.0) * sqrt14 * Y3_n2 * R4) * V[3] + (7.0 * sqrt6 * Y5_n4 + 10.0 * sqrt2 * Y5_n2 + 8.0 * sqrt14 * Y3_n2 * R2) * V[4] - 15.0 * sqrt14 * Y3_n2 * V[5];
    d[14] = ((-9.0 / 286.0) * sqrt110 * Y7_n5 + (-3.0 / 286.0) * sqrt30 * Y7_n1 + (1.0 / 39.0) * sqrt15 * Y5_n5 * R2 + (2.0 / 39.0) * sqrt14 * Y5_n1 * R2 + (-1.0 / 33.0) * sqrt35 * Y3_n1 * R4) * V[0] + ((-1.0 / 6.0) * sqrt15 * Y5_n5 + (-1.0 / 3.0) * sqrt14 * Y5_n1 + (1.0 / 3.0) * sqrt35 * Y3_n1 * R2) * V[1] - 0.75 * sqrt35 * Y3_n1 * V[2] + ((-63.0 / 143.0) * sqrt110 * Y7_n5 + (-21.0 / 143.0) * sqrt30 * Y7_n1 + (14.0 / 39.0) * sqrt15 * Y5_n5 * R2 + (28.0 / 39.0) * sqrt14 * Y5_n1 * R2 + (-14.0 / 33.0) * sqrt35 * Y3_n1 * R4) * V[3] + (-2.0 * sqrt15 * Y5_n5 - 4.0 * sqrt14 * Y5_n1 + 4.0 * sqrt35 * Y3_n1 * R2) * V[4] - 7.5 * sqrt35 * Y3_n1 * V[5];
    d[11] = ((9.0 / 143.0) * sqrt5 * Y7_p2 + (9.0 / 286.0) * sqrt110 * Y7_p4 + (-5.0 / 39.0) * sqrt2 * Y5_p2 * R2 + (7.0 / 78.0) * sqrt6 * Y5_p4 * R2 + (-2.0 / 33.0) * sqrt14 * Y3_p2 * R4) * V[0] + ((5.0 / 6.0) * sqrt2 * Y5_p2 + (-7.0 / 12.0) * sqrt6 * Y5_p4 + (2.0 / 3.0) * sqrt14 * Y3_p2 * R2) * V[1] - 1.5 * sqrt14 * Y3_p2 * V[2] + ((126.0 / 143.0) * sqrt5 * Y7_p2 + (63.0 / 143.0) * sqrt110 * Y7_p4 + (-70.0 / 39.0) * sqrt2 * Y5_p2 * R2 + (49.0 / 39.0) * sqrt6 * Y5_p4 * R2 + (-28.0 / 33.0) * sqrt14 * Y3_p2 * R4) * V[3] + (10.0 * sqrt2 * Y5_p2 - 7.0 * sqrt6 * Y5_p4 + 8.0 * sqrt14 * Y3_p2 * R2) * V[4] - 15.0 * sqrt14 * Y3_p2 * V[5];
    d[10] = ((-3.0 / 286.0) * sqrt30 * Y7_p1 + (9.0 / 286.0) * sqrt110 * Y7_p5 + (2.0 / 39.0) * sqrt14 * Y5_p1 * R2 + (-1.0 / 39.0) * sqrt15 * Y5_p5 * R2 + (-1.0 / 33.0) * sqrt35 * Y3_p1 * R4) * V[0] + ((-1.0 / 3.0) * sqrt14 * Y5_p1 + (1.0 / 6.0) * sqrt15 * Y5_p5 + (1.0 / 3.0) * sqrt35 * Y3_p1 * R2) * V[1] - 0.75 * sqrt35 * Y3_p1 * V[2] + ((-21.0 / 143.0) * sqrt30 * Y7_p1 + (63.0 / 143.0) * sqrt110 * Y7_p5 + (28.0 / 39.0) * sqrt14 * Y5_p1 * R2 + (-14.0 / 39.0) * sqrt15 * Y5_p5 * R2 + (-14.0 / 33.0) * sqrt35 * Y3_p1 * R4) * V[3] + (-4.0 * sqrt14 * Y5_p1 + 2.0 * sqrt15 * Y5_p5 + 4.0 * sqrt35 * Y3_p1 * R2) * V[4] - 7.5 * sqrt35 * Y3_p1 * V[5];
    d[7] = ((-3.0 / 143.0) * sqrt165 * Y7_n4 + (2.0 / 13.0) * Y5_n4 * R2) * V[0] - Y5_n4 * V[1] + ((-42.0 / 143.0) * sqrt165 * Y7_n4 + (28.0 / 13.0) * Y5_n4 * R2) * V[3] - 12.0 * Y5_n4 * V[4];
    d[8] = ((-6.0 / 143.0) * sqrt55 * Y7_n5 + (6.0 / 143.0) * sqrt5 * Y7_n3 + (-1.0 / 26.0) * sqrt30 * Y5_n5 * R2 + (-7.0 / 78.0) * sqrt6 * Y5_n3 * R2 + (-1.0 / 33.0) * sqrt42 * Y3_n3 * R4) * V[0] + (0.25 * sqrt30 * Y5_n5 + (7.0 / 12.0) * sqrt6 * Y5_n3 + (1.0 / 3.0) * sqrt42 * Y3_n3 * R2) * V[1] - 0.75 * sqrt42 * Y3_n3 * V[2] + ((-84.0 / 143.0) * sqrt55 * Y7_n5 + (84.0 / 143.0) * sqrt5 * Y7_n3 + (-7.0 / 13.0) * sqrt30 * Y5_n5 * R2 + (-49.0 / 39.0) * sqrt6 * Y5_n3 * R2 + (-14.0 / 33.0) * sqrt42 * Y3_n3 * R4) * V[3] + (3.0 * sqrt30 * Y5_n5 + 7.0 * sqrt6 * Y5_n3 + 4.0 * sqrt42 * Y3_n3 * R2) * V[4] - 7.5 * sqrt42 * Y3_n3 * V[5];
    d[9] = ((-3.0 / 286.0) * sqrt1430 * Y7_n6 + (-3.0 / 286.0) * sqrt10 * Y7_n2 + (2.0 / 13.0) * Y5_n2 * R2 + (-1.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[0] + (-Y5_n2 + sqrt7 * Y3_n2 * R2) * V[1] - 2.25 * sqrt7 * Y3_n2 * V[2] + ((-21.0 / 143.0) * sqrt1430 * Y7_n6 + (-21.0 / 143.0) * sqrt10 * Y7_n2 + (28.0 / 13.0) * Y5_n2 * R2 + (-14.0 / 11.0) * sqrt7 * Y3_n2 * R4) * V[3] + (-12.0 * Y5_n2 + 12.0 * sqrt7 * Y3_n2 * R2) * V[4] - 22.5 * sqrt7 * Y3_n2 * V[5];
    d[6] = ((6.0 / 143.0) * sqrt5 * Y7_p3 + (6.0 / 143.0) * sqrt55 * Y7_p5 + (-7.0 / 78.0) * sqrt6 * Y5_p3 * R2 + (1.0 / 26.0) * sqrt30 * Y5_p5 * R2 + (-1.0 / 33.0) * sqrt42 * Y3_p3 * R4) * V[0] + ((7.0 / 12.0) * sqrt6 * Y5_p3 - 0.25 * sqrt30 * Y5_p5 + (1.0 / 3.0) * sqrt42 * Y3_p3 * R2) * V[1] - 0.75 * sqrt42 * Y3_p3 * V[2] + ((84.0 / 143.0) * sqrt5 * Y7_p3 + (84.0 / 143.0) * sqrt55 * Y7_p5 + (-49.0 / 39.0) * sqrt6 * Y5_p3 * R2 + (7.0 / 13.0) * sqrt30 * Y5_p5 * R2 + (-14.0 / 33.0) * sqrt42 * Y3_p3 * R4) * V[3] + (7.0 * sqrt6 * Y5_p3 - 3.0 * sqrt30 * Y5_p5 + 4.0 * sqrt42 * Y3_p3 * R2) * V[4] - 7.5 * sqrt42 * Y3_p3 * V[5];
    d[5] = ((-3.0 / 286.0) * sqrt10 * Y7_p2 + (3.0 / 286.0) * sqrt1430 * Y7_p6 + (2.0 / 13.0) * Y5_p2 * R2 + (-1.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[0] + (-Y5_p2 + sqrt7 * Y3_p2 * R2) * V[1] - 2.25 * sqrt7 * Y3_p2 * V[2] + ((-21.0 / 143.0) * sqrt10 * Y7_p2 + (21.0 / 143.0) * sqrt1430 * Y7_p6 + (28.0 / 13.0) * Y5_p2 * R2 + (-14.0 / 11.0) * sqrt7 * Y3_p2 * R4) * V[3] + (-12.0 * Y5_p2 + 12.0 * sqrt7 * Y3_p2 * R2) * V[4] - 22.5 * sqrt7 * Y3_p2 * V[5];
    d[2] = ((-3.0 / 143.0) * sqrt66 * Y7_n5 + (5.0 / 13.0) * Y5_n5 * R2) * V[0] - 2.5 * Y5_n5 * V[1] + ((-42.0 / 143.0) * sqrt66 * Y7_n5 + (70.0 / 13.0) * Y5_n5 * R2) * V[3] - 30.0 * Y5_n5 * V[4];
    d[3] = ((-3.0 / 143.0) * sqrt143 * Y7_n6 + (3.0 / 286.0) * sqrt22 * Y7_n4 + (-1.0 / 26.0) * sqrt30 * Y5_n4 * R2) * V[0] + 0.25 * sqrt30 * Y5_n4 * V[1] + ((-42.0 / 143.0) * sqrt143 * Y7_n6 + (21.0 / 143.0) * sqrt22 * Y7_n4 + (-7.0 / 13.0) * sqrt30 * Y5_n4 * R2) * V[3] + 3.0 * sqrt30 * Y5_n4 * V[4];
    d[4] = ((-3.0 / 286.0) * sqrt2002 * Y7_n7 + (-3.0 / 286.0) * sqrt2 * Y7_n3 + (1.0 / 39.0) * sqrt15 * Y5_n3 * R2 + (-1.0 / 33.0) * sqrt105 * Y3_n3 * R4) * V[0] + ((-1.0 / 6.0) * sqrt15 * Y5_n3 + (1.0 / 3.0) * sqrt105 * Y3_n3 * R2) * V[1] - 0.75 * sqrt105 * Y3_n3 * V[2] + ((-21.0 / 143.0) * sqrt2002 * Y7_n7 + (-21.0 / 143.0) * sqrt2 * Y7_n3 + (14.0 / 39.0) * sqrt15 * Y5_n3 * R2 + (-14.0 / 33.0) * sqrt105 * Y3_n3 * R4) * V[3] + (-2.0 * sqrt15 * Y5_n3 + 4.0 * sqrt105 * Y3_n3 * R2) * V[4] - 7.5 * sqrt105 * Y3_n3 * V[5];
    d[1] = ((3.0 / 286.0) * sqrt22 * Y7_p4 + (3.0 / 143.0) * sqrt143 * Y7_p6 + (-1.0 / 26.0) * sqrt30 * Y5_p4 * R2) * V[0] + 0.25 * sqrt30 * Y5_p4 * V[1] + ((21.0 / 143.0) * sqrt22 * Y7_p4 + (42.0 / 143.0) * sqrt143 * Y7_p6 + (-7.0 / 13.0) * sqrt30 * Y5_p4 * R2) * V[3] + 3.0 * sqrt30 * Y5_p4 * V[4];
    d[0] = ((-3.0 / 286.0) * sqrt2 * Y7_p3 + (3.0 / 286.0) * sqrt2002 * Y7_p7 + (1.0 / 39.0) * sqrt15 * Y5_p3 * R2 + (-1.0 / 33.0) * sqrt105 * Y3_p3 * R4) * V[0] + ((-1.0 / 6.0) * sqrt15 * Y5_p3 + (1.0 / 3.0) * sqrt105 * Y3_p3 * R2) * V[1] - 0.75 * sqrt105 * Y3_p3 * V[2] + ((-21.0 / 143.0) * sqrt2 * Y7_p3 + (21.0 / 143.0) * sqrt2002 * Y7_p7 + (14.0 / 39.0) * sqrt15 * Y5_p3 * R2 + (-14.0 / 33.0) * sqrt105 * Y3_p3 * R4) * V[3] + (-2.0 * sqrt15 * Y5_p3 + 4.0 * sqrt105 * Y3_p3 * R2) * V[4] - 7.5 * sqrt105 * Y3_p3 * V[5];
}

}  // namespace kinab