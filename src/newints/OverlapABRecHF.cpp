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

#include "OverlapABRecHF.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_h_f(
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
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt66 = std::sqrt(66.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt210 = std::sqrt(210.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^3 · β^5 · p^{-8} · (s|s)
    //   V[1] ↔ α^2 · β^4 · p^{-7} · (s|s)
    //   V[2] ↔ α · β^3 · p^{-6} · (s|s)
    //   V[3] ↔ β^2 · p^{-5} · (s|s)
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

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha3 * beta5 * pinv8;
            V[1] += cab_ss * alpha2 * beta4 * pinv7;
            V[2] += cab_ss * alpha * beta3 * pinv6;
            V[3] += cab_ss * beta2 * pinv5;
        }
    }

    // ---- Phase 3: fused M·V → 11 × 7 spherical block ----
    const auto Y2_n2 = harm::Y_ll_2_m_n2(AB_x, AB_y, AB_z);
    const auto Y2_n1 = harm::Y_ll_2_m_n1(AB_x, AB_y, AB_z);
    const auto Y2_p0 = harm::Y_ll_2_m_p0(AB_x, AB_y, AB_z);
    const auto Y2_p1 = harm::Y_ll_2_m_p1(AB_x, AB_y, AB_z);
    const auto Y2_p2 = harm::Y_ll_2_m_p2(AB_x, AB_y, AB_z);
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
    auto *d = buffer;
    d[38] = -Y5_p0 * Y3_p0 * V[0] + ((35.0 / 22.0) * Y6_p0 + (180.0 / 77.0) * Y4_p0 * R2 + (25.0 / 7.0) * Y2_p0 * R4) * V[1] + ((-45.0 / 7.0) * Y4_p0 + (-225.0 / 14.0) * Y2_p0 * R2) * V[2] + 18.75 * Y2_p0 * V[3];
    d[39] = -Y5_p0 * Y3_p1 * V[0] + ((5.0 / 22.0) * sqrt14 * Y6_p1 + (3.0 / 77.0) * sqrt15 * Y4_p1 * R2 + (-25.0 / 14.0) * sqrt2 * Y2_p1 * R4) * V[1] + ((-3.0 / 28.0) * sqrt15 * Y4_p1 + (225.0 / 28.0) * sqrt2 * Y2_p1 * R2) * V[2] - 9.375 * sqrt2 * Y2_p1 * V[3];
    d[40] = -Y5_p0 * Y3_p2 * V[0] + ((-5.0 / 22.0) * sqrt14 * Y6_p2 + (-120.0 / 77.0) * sqrt3 * Y4_p2 * R2 + (5.0 / 14.0) * sqrt5 * Y2_p2 * R4) * V[1] + ((30.0 / 7.0) * sqrt3 * Y4_p2 + (-45.0 / 28.0) * sqrt5 * Y2_p2 * R2) * V[2] + 1.875 * sqrt5 * Y2_p2 * V[3];
    d[41] = -Y5_p0 * Y3_p3 * V[0] + ((-5.0 / 11.0) * sqrt21 * Y6_p3 + (45.0 / 77.0) * sqrt7 * Y4_p3 * R2) * V[1] + (-45.0 / 28.0) * sqrt7 * Y4_p3 * V[2];
    d[37] = -Y5_p0 * Y3_n1 * V[0] + ((5.0 / 22.0) * sqrt14 * Y6_n1 + (3.0 / 77.0) * sqrt15 * Y4_n1 * R2 + (-25.0 / 14.0) * sqrt2 * Y2_n1 * R4) * V[1] + ((-3.0 / 28.0) * sqrt15 * Y4_n1 + (225.0 / 28.0) * sqrt2 * Y2_n1 * R2) * V[2] - 9.375 * sqrt2 * Y2_n1 * V[3];
    d[36] = -Y5_p0 * Y3_n2 * V[0] + ((-5.0 / 22.0) * sqrt14 * Y6_n2 + (-120.0 / 77.0) * sqrt3 * Y4_n2 * R2 + (5.0 / 14.0) * sqrt5 * Y2_n2 * R4) * V[1] + ((30.0 / 7.0) * sqrt3 * Y4_n2 + (-45.0 / 28.0) * sqrt5 * Y2_n2 * R2) * V[2] + 1.875 * sqrt5 * Y2_n2 * V[3];
    d[35] = -Y5_p0 * Y3_n3 * V[0] + ((-5.0 / 11.0) * sqrt21 * Y6_n3 + (45.0 / 77.0) * sqrt7 * Y4_n3 * R2) * V[1] + (-45.0 / 28.0) * sqrt7 * Y4_n3 * V[2];
    d[45] = -Y5_p1 * Y3_p0 * V[0] + ((5.0 / 22.0) * sqrt35 * Y6_p1 + (57.0 / 77.0) * sqrt6 * Y4_p1 * R2 + (10.0 / 7.0) * sqrt5 * Y2_p1 * R4) * V[1] + ((-57.0 / 28.0) * sqrt6 * Y4_p1 + (-45.0 / 7.0) * sqrt5 * Y2_p1 * R2) * V[2] + 7.5 * sqrt5 * Y2_p1 * V[3];
    d[46] = -Y5_p1 * Y3_p1 * V[0] + ((5.0 / 22.0) * sqrt21 * Y6_p2 + (3.0 / 7.0) * sqrt10 * Y4_p0 * R2 + (129.0 / 154.0) * sqrt2 * Y4_p2 * R2 + (15.0 / 14.0) * sqrt10 * Y2_p0 * R4 + (-5.0 / 28.0) * sqrt30 * Y2_p2 * R4) * V[1] + ((-33.0 / 28.0) * sqrt10 * Y4_p0 + (-129.0 / 56.0) * sqrt2 * Y4_p2 + (-135.0 / 28.0) * sqrt10 * Y2_p0 * R2 + (45.0 / 56.0) * sqrt30 * Y2_p2 * R2) * V[2] + (5.625 * sqrt10 * Y2_p0 - 0.9375 * sqrt30 * Y2_p2) * V[3];
    d[47] = -Y5_p1 * Y3_p2 * V[0] + ((5.0 / 22.0) * sqrt21 * Y6_p1 + (75.0 / 154.0) * sqrt10 * Y4_p1 * R2 + (-3.0 / 14.0) * sqrt70 * Y4_p3 * R2 + (-5.0 / 7.0) * sqrt3 * Y2_p1 * R4) * V[1] + ((-75.0 / 56.0) * sqrt10 * Y4_p1 + (33.0 / 56.0) * sqrt70 * Y4_p3 + (45.0 / 14.0) * sqrt3 * Y2_p1 * R2) * V[2] - 3.75 * sqrt3 * Y2_p1 * V[3];
    d[48] = -Y5_p1 * Y3_p3 * V[0] + ((5.0 / 22.0) * sqrt35 * Y6_p2 + (-5.0 / 22.0) * sqrt42 * Y6_p4 + (-45.0 / 154.0) * sqrt30 * Y4_p2 * R2 + (3.0 / 77.0) * sqrt210 * Y4_p4 * R2 + (5.0 / 28.0) * sqrt2 * Y2_p2 * R4) * V[1] + ((45.0 / 56.0) * sqrt30 * Y4_p2 + (-3.0 / 28.0) * sqrt210 * Y4_p4 + (-45.0 / 56.0) * sqrt2 * Y2_p2 * R2) * V[2] + 0.9375 * sqrt2 * Y2_p2 * V[3];
    d[44] = -Y5_p1 * Y3_n1 * V[0] + ((5.0 / 22.0) * sqrt21 * Y6_n2 + (129.0 / 154.0) * sqrt2 * Y4_n2 * R2 + (-5.0 / 28.0) * sqrt30 * Y2_n2 * R4) * V[1] + ((-129.0 / 56.0) * sqrt2 * Y4_n2 + (45.0 / 56.0) * sqrt30 * Y2_n2 * R2) * V[2] - 0.9375 * sqrt30 * Y2_n2 * V[3];
    d[43] = -Y5_p1 * Y3_n2 * V[0] + ((5.0 / 22.0) * sqrt21 * Y6_n1 + (-3.0 / 14.0) * sqrt70 * Y4_n3 * R2 + (75.0 / 154.0) * sqrt10 * Y4_n1 * R2 + (-5.0 / 7.0) * sqrt3 * Y2_n1 * R4) * V[1] + ((33.0 / 56.0) * sqrt70 * Y4_n3 + (-75.0 / 56.0) * sqrt10 * Y4_n1 + (45.0 / 14.0) * sqrt3 * Y2_n1 * R2) * V[2] - 3.75 * sqrt3 * Y2_n1 * V[3];
    d[42] = -Y5_p1 * Y3_n3 * V[0] + ((-5.0 / 22.0) * sqrt42 * Y6_n4 + (5.0 / 22.0) * sqrt35 * Y6_n2 + (3.0 / 77.0) * sqrt210 * Y4_n4 * R2 + (-45.0 / 154.0) * sqrt30 * Y4_n2 * R2 + (5.0 / 28.0) * sqrt2 * Y2_n2 * R4) * V[1] + ((-3.0 / 28.0) * sqrt210 * Y4_n4 + (45.0 / 56.0) * sqrt30 * Y4_n2 + (-45.0 / 56.0) * sqrt2 * Y2_n2 * R2) * V[2] + 0.9375 * sqrt2 * Y2_n2 * V[3];
    d[52] = -Y5_p2 * Y3_p0 * V[0] + ((5.0 / 11.0) * sqrt2 * Y6_p2 + (6.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (5.0 / 14.0) * sqrt35 * Y2_p2 * R4) * V[1] + ((-3.0 / 14.0) * sqrt21 * Y4_p2 + (-45.0 / 28.0) * sqrt35 * Y2_p2 * R2) * V[2] + 1.875 * sqrt35 * Y2_p2 * V[3];
    d[53] = -Y5_p2 * Y3_p1 * V[0] + ((5.0 / 44.0) * sqrt30 * Y6_p1 + (15.0 / 22.0) * sqrt3 * Y6_p3 + (48.0 / 77.0) * sqrt7 * Y4_p1 * R2 + (21.0 / 11.0) * Y4_p3 * R2 + (5.0 / 28.0) * sqrt210 * Y2_p1 * R4) * V[1] + ((-12.0 / 7.0) * sqrt7 * Y4_p1 - 5.25 * Y4_p3 + (-45.0 / 56.0) * sqrt210 * Y2_p1 * R2) * V[2] + 0.9375 * sqrt210 * Y2_p1 * V[3];
    d[54] = -Y5_p2 * Y3_p2 * V[0] + ((-15.0 / 22.0) * sqrt7 * Y6_p0 + (15.0 / 22.0) * Y6_p4 + (-30.0 / 77.0) * sqrt7 * Y4_p0 * R2 + (-6.0 / 11.0) * sqrt5 * Y4_p4 * R2 + (15.0 / 14.0) * sqrt7 * Y2_p0 * R4) * V[1] + ((15.0 / 14.0) * sqrt7 * Y4_p0 + 1.5 * sqrt5 * Y4_p4 + (-135.0 / 28.0) * sqrt7 * Y2_p0 * R2) * V[2] + 5.625 * sqrt7 * Y2_p0 * V[3];
    d[55] = -Y5_p2 * Y3_p3 * V[0] + ((-35.0 / 44.0) * sqrt2 * Y6_p1 + (-5.0 / 22.0) * sqrt33 * Y6_p5 + (15.0 / 77.0) * sqrt105 * Y4_p1 * R2 + (-5.0 / 28.0) * sqrt14 * Y2_p1 * R4) * V[1] + ((-15.0 / 28.0) * sqrt105 * Y4_p1 + (45.0 / 56.0) * sqrt14 * Y2_p1 * R2) * V[2] - 0.9375 * sqrt14 * Y2_p1 * V[3];
    d[51] = -Y5_p2 * Y3_n1 * V[0] + ((15.0 / 22.0) * sqrt3 * Y6_n3 + (-5.0 / 44.0) * sqrt30 * Y6_n1 + (21.0 / 11.0) * Y4_n3 * R2 + (-48.0 / 77.0) * sqrt7 * Y4_n1 * R2 + (-5.0 / 28.0) * sqrt210 * Y2_n1 * R4) * V[1] + (-5.25 * Y4_n3 + (12.0 / 7.0) * sqrt7 * Y4_n1 + (45.0 / 56.0) * sqrt210 * Y2_n1 * R2) * V[2] - 0.9375 * sqrt210 * Y2_n1 * V[3];
    d[50] = -Y5_p2 * Y3_n2 * V[0] + ((15.0 / 22.0) * Y6_n4 + (-6.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[1] + 1.5 * sqrt5 * Y4_n4 * V[2];
    d[49] = -Y5_p2 * Y3_n3 * V[0] + ((-5.0 / 22.0) * sqrt33 * Y6_n5 + (-35.0 / 44.0) * sqrt2 * Y6_n1 + (15.0 / 77.0) * sqrt105 * Y4_n1 * R2 + (-5.0 / 28.0) * sqrt14 * Y2_n1 * R4) * V[1] + ((-15.0 / 28.0) * sqrt105 * Y4_n1 + (45.0 / 56.0) * sqrt14 * Y2_n1 * R2) * V[2] - 0.9375 * sqrt14 * Y2_n1 * V[3];
    d[59] = -Y5_p3 * Y3_p0 * V[0] + ((-5.0 / 22.0) * sqrt3 * Y6_p3 + (-18.0 / 11.0) * Y4_p3 * R2) * V[1] + 4.5 * Y4_p3 * V[2];
    d[60] = -Y5_p3 * Y3_p1 * V[0] + ((35.0 / 44.0) * sqrt2 * Y6_p2 + (5.0 / 22.0) * sqrt15 * Y6_p4 + (27.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (12.0 / 11.0) * sqrt3 * Y4_p4 * R2 + (5.0 / 14.0) * sqrt35 * Y2_p2 * R4) * V[1] + ((-27.0 / 28.0) * sqrt21 * Y4_p2 - 3.0 * sqrt3 * Y4_p4 + (-45.0 / 28.0) * sqrt35 * Y2_p2 * R2) * V[2] + 1.875 * sqrt35 * Y2_p2 * V[3];
    d[61] = -Y5_p3 * Y3_p2 * V[0] + ((-10.0 / 11.0) * sqrt2 * Y6_p1 + (5.0 / 22.0) * sqrt33 * Y6_p5 + (3.0 / 77.0) * sqrt105 * Y4_p1 * R2 + (5.0 / 7.0) * sqrt14 * Y2_p1 * R4) * V[1] + ((-3.0 / 28.0) * sqrt105 * Y4_p1 + (-45.0 / 14.0) * sqrt14 * Y2_p1 * R2) * V[2] + 3.75 * sqrt14 * Y2_p1 * V[3];
    d[62] = -Y5_p3 * Y3_p3 * V[0] + ((5.0 / 11.0) * sqrt7 * Y6_p0 + (-5.0 / 44.0) * sqrt66 * Y6_p6 + (-90.0 / 77.0) * sqrt7 * Y4_p0 * R2 + (5.0 / 7.0) * sqrt7 * Y2_p0 * R4) * V[1] + ((45.0 / 14.0) * sqrt7 * Y4_p0 + (-45.0 / 14.0) * sqrt7 * Y2_p0 * R2) * V[2] + 3.75 * sqrt7 * Y2_p0 * V[3];
    d[58] = -Y5_p3 * Y3_n1 * V[0] + ((5.0 / 22.0) * sqrt15 * Y6_n4 + (-35.0 / 44.0) * sqrt2 * Y6_n2 + (12.0 / 11.0) * sqrt3 * Y4_n4 * R2 + (-27.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (-5.0 / 14.0) * sqrt35 * Y2_n2 * R4) * V[1] + (-3.0 * sqrt3 * Y4_n4 + (27.0 / 28.0) * sqrt21 * Y4_n2 + (45.0 / 28.0) * sqrt35 * Y2_n2 * R2) * V[2] - 1.875 * sqrt35 * Y2_n2 * V[3];
    d[57] = -Y5_p3 * Y3_n2 * V[0] + ((5.0 / 22.0) * sqrt33 * Y6_n5 + (10.0 / 11.0) * sqrt2 * Y6_n1 + (-3.0 / 77.0) * sqrt105 * Y4_n1 * R2 + (-5.0 / 7.0) * sqrt14 * Y2_n1 * R4) * V[1] + ((3.0 / 28.0) * sqrt105 * Y4_n1 + (45.0 / 14.0) * sqrt14 * Y2_n1 * R2) * V[2] - 3.75 * sqrt14 * Y2_n1 * V[3];
    d[56] = -Y5_p3 * Y3_n3 * V[0] + (-5.0 / 44.0) * sqrt66 * Y6_n6 * V[1];
    d[66] = -Y5_p4 * Y3_p0 * V[0] + ((-15.0 / 22.0) * sqrt5 * Y6_p4 + (-36.0 / 11.0) * Y4_p4 * R2) * V[1] + 9.0 * Y4_p4 * V[2];
    d[67] = -Y5_p4 * Y3_p1 * V[0] + ((15.0 / 11.0) * Y6_p3 + (3.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[1] - 0.75 * sqrt3 * Y4_p3 * V[2];
    d[68] = -Y5_p4 * Y3_p2 * V[0] + ((-15.0 / 44.0) * sqrt10 * Y6_p2 + (15.0 / 44.0) * sqrt22 * Y6_p6 + (12.0 / 77.0) * sqrt105 * Y4_p2 * R2 + (15.0 / 14.0) * sqrt7 * Y2_p2 * R4) * V[1] + ((-3.0 / 7.0) * sqrt105 * Y4_p2 + (-135.0 / 28.0) * sqrt7 * Y2_p2 * R2) * V[2] + 5.625 * sqrt7 * Y2_p2 * V[3];
    d[69] = -Y5_p4 * Y3_p3 * V[0] + ((5.0 / 22.0) * sqrt6 * Y6_p1 + (-27.0 / 77.0) * sqrt35 * Y4_p1 * R2 + (5.0 / 14.0) * sqrt42 * Y2_p1 * R4) * V[1] + ((27.0 / 28.0) * sqrt35 * Y4_p1 + (-45.0 / 28.0) * sqrt42 * Y2_p1 * R2) * V[2] + 1.875 * sqrt42 * Y2_p1 * V[3];
    d[65] = -Y5_p4 * Y3_n1 * V[0] + ((-15.0 / 11.0) * Y6_n3 + (-3.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[1] + 0.75 * sqrt3 * Y4_n3 * V[2];
    d[64] = -Y5_p4 * Y3_n2 * V[0] + ((15.0 / 44.0) * sqrt22 * Y6_n6 + (15.0 / 44.0) * sqrt10 * Y6_n2 + (-12.0 / 77.0) * sqrt105 * Y4_n2 * R2 + (-15.0 / 14.0) * sqrt7 * Y2_n2 * R4) * V[1] + ((3.0 / 7.0) * sqrt105 * Y4_n2 + (135.0 / 28.0) * sqrt7 * Y2_n2 * R2) * V[2] - 5.625 * sqrt7 * Y2_n2 * V[3];
    d[63] = -Y5_p4 * Y3_n3 * V[0] + ((-5.0 / 22.0) * sqrt6 * Y6_n1 + (27.0 / 77.0) * sqrt35 * Y4_n1 * R2 + (-5.0 / 14.0) * sqrt42 * Y2_n1 * R4) * V[1] + ((-27.0 / 28.0) * sqrt35 * Y4_n1 + (45.0 / 28.0) * sqrt42 * Y2_n1 * R2) * V[2] - 1.875 * sqrt42 * Y2_n1 * V[3];
    d[73] = -Y5_p5 * Y3_p0 * V[0] + (-15.0 / 22.0) * sqrt11 * Y6_p5 * V[1];
    d[74] = -Y5_p5 * Y3_p1 * V[0] + ((15.0 / 22.0) * sqrt3 * Y6_p4 + (-15.0 / 44.0) * sqrt22 * Y6_p6 + (-6.0 / 11.0) * sqrt15 * Y4_p4 * R2) * V[1] + 1.5 * sqrt15 * Y4_p4 * V[2];
    d[75] = -Y5_p5 * Y3_p2 * V[0] + ((-15.0 / 22.0) * Y6_p3 + (15.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[1] - 3.75 * sqrt3 * Y4_p3 * V[2];
    d[76] = -Y5_p5 * Y3_p3 * V[0] + ((5.0 / 44.0) * sqrt6 * Y6_p2 + (-45.0 / 77.0) * sqrt7 * Y4_p2 * R2 + (5.0 / 14.0) * sqrt105 * Y2_p2 * R4) * V[1] + ((45.0 / 28.0) * sqrt7 * Y4_p2 + (-45.0 / 28.0) * sqrt105 * Y2_p2 * R2) * V[2] + 1.875 * sqrt105 * Y2_p2 * V[3];
    d[72] = -Y5_p5 * Y3_n1 * V[0] + ((-15.0 / 44.0) * sqrt22 * Y6_n6 + (-15.0 / 22.0) * sqrt3 * Y6_n4 + (6.0 / 11.0) * sqrt15 * Y4_n4 * R2) * V[1] - 1.5 * sqrt15 * Y4_n4 * V[2];
    d[71] = -Y5_p5 * Y3_n2 * V[0] + ((15.0 / 22.0) * Y6_n3 + (-15.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[1] + 3.75 * sqrt3 * Y4_n3 * V[2];
    d[70] = -Y5_p5 * Y3_n3 * V[0] + ((-5.0 / 44.0) * sqrt6 * Y6_n2 + (45.0 / 77.0) * sqrt7 * Y4_n2 * R2 + (-5.0 / 14.0) * sqrt105 * Y2_n2 * R4) * V[1] + ((-45.0 / 28.0) * sqrt7 * Y4_n2 + (45.0 / 28.0) * sqrt105 * Y2_n2 * R2) * V[2] - 1.875 * sqrt105 * Y2_n2 * V[3];
    d[31] = -Y5_n1 * Y3_p0 * V[0] + ((5.0 / 22.0) * sqrt35 * Y6_n1 + (57.0 / 77.0) * sqrt6 * Y4_n1 * R2 + (10.0 / 7.0) * sqrt5 * Y2_n1 * R4) * V[1] + ((-57.0 / 28.0) * sqrt6 * Y4_n1 + (-45.0 / 7.0) * sqrt5 * Y2_n1 * R2) * V[2] + 7.5 * sqrt5 * Y2_n1 * V[3];
    d[32] = -Y5_n1 * Y3_p1 * V[0] + ((5.0 / 22.0) * sqrt21 * Y6_n2 + (129.0 / 154.0) * sqrt2 * Y4_n2 * R2 + (-5.0 / 28.0) * sqrt30 * Y2_n2 * R4) * V[1] + ((-129.0 / 56.0) * sqrt2 * Y4_n2 + (45.0 / 56.0) * sqrt30 * Y2_n2 * R2) * V[2] - 0.9375 * sqrt30 * Y2_n2 * V[3];
    d[33] = -Y5_n1 * Y3_p2 * V[0] + ((-5.0 / 22.0) * sqrt21 * Y6_n1 + (-3.0 / 14.0) * sqrt70 * Y4_n3 * R2 + (-75.0 / 154.0) * sqrt10 * Y4_n1 * R2 + (5.0 / 7.0) * sqrt3 * Y2_n1 * R4) * V[1] + ((33.0 / 56.0) * sqrt70 * Y4_n3 + (75.0 / 56.0) * sqrt10 * Y4_n1 + (-45.0 / 14.0) * sqrt3 * Y2_n1 * R2) * V[2] + 3.75 * sqrt3 * Y2_n1 * V[3];
    d[34] = -Y5_n1 * Y3_p3 * V[0] + ((-5.0 / 22.0) * sqrt42 * Y6_n4 + (-5.0 / 22.0) * sqrt35 * Y6_n2 + (3.0 / 77.0) * sqrt210 * Y4_n4 * R2 + (45.0 / 154.0) * sqrt30 * Y4_n2 * R2 + (-5.0 / 28.0) * sqrt2 * Y2_n2 * R4) * V[1] + ((-3.0 / 28.0) * sqrt210 * Y4_n4 + (-45.0 / 56.0) * sqrt30 * Y4_n2 + (45.0 / 56.0) * sqrt2 * Y2_n2 * R2) * V[2] - 0.9375 * sqrt2 * Y2_n2 * V[3];
    d[30] = -Y5_n1 * Y3_n1 * V[0] + ((-5.0 / 22.0) * sqrt21 * Y6_p2 + (3.0 / 7.0) * sqrt10 * Y4_p0 * R2 + (-129.0 / 154.0) * sqrt2 * Y4_p2 * R2 + (15.0 / 14.0) * sqrt10 * Y2_p0 * R4 + (5.0 / 28.0) * sqrt30 * Y2_p2 * R4) * V[1] + ((-33.0 / 28.0) * sqrt10 * Y4_p0 + (129.0 / 56.0) * sqrt2 * Y4_p2 + (-135.0 / 28.0) * sqrt10 * Y2_p0 * R2 + (-45.0 / 56.0) * sqrt30 * Y2_p2 * R2) * V[2] + (5.625 * sqrt10 * Y2_p0 + 0.9375 * sqrt30 * Y2_p2) * V[3];
    d[29] = -Y5_n1 * Y3_n2 * V[0] + ((5.0 / 22.0) * sqrt21 * Y6_p1 + (75.0 / 154.0) * sqrt10 * Y4_p1 * R2 + (3.0 / 14.0) * sqrt70 * Y4_p3 * R2 + (-5.0 / 7.0) * sqrt3 * Y2_p1 * R4) * V[1] + ((-75.0 / 56.0) * sqrt10 * Y4_p1 + (-33.0 / 56.0) * sqrt70 * Y4_p3 + (45.0 / 14.0) * sqrt3 * Y2_p1 * R2) * V[2] - 3.75 * sqrt3 * Y2_p1 * V[3];
    d[28] = -Y5_n1 * Y3_n3 * V[0] + ((5.0 / 22.0) * sqrt35 * Y6_p2 + (5.0 / 22.0) * sqrt42 * Y6_p4 + (-45.0 / 154.0) * sqrt30 * Y4_p2 * R2 + (-3.0 / 77.0) * sqrt210 * Y4_p4 * R2 + (5.0 / 28.0) * sqrt2 * Y2_p2 * R4) * V[1] + ((45.0 / 56.0) * sqrt30 * Y4_p2 + (3.0 / 28.0) * sqrt210 * Y4_p4 + (-45.0 / 56.0) * sqrt2 * Y2_p2 * R2) * V[2] + 0.9375 * sqrt2 * Y2_p2 * V[3];
    d[24] = -Y5_n2 * Y3_p0 * V[0] + ((5.0 / 11.0) * sqrt2 * Y6_n2 + (6.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (5.0 / 14.0) * sqrt35 * Y2_n2 * R4) * V[1] + ((-3.0 / 14.0) * sqrt21 * Y4_n2 + (-45.0 / 28.0) * sqrt35 * Y2_n2 * R2) * V[2] + 1.875 * sqrt35 * Y2_n2 * V[3];
    d[25] = -Y5_n2 * Y3_p1 * V[0] + ((15.0 / 22.0) * sqrt3 * Y6_n3 + (5.0 / 44.0) * sqrt30 * Y6_n1 + (21.0 / 11.0) * Y4_n3 * R2 + (48.0 / 77.0) * sqrt7 * Y4_n1 * R2 + (5.0 / 28.0) * sqrt210 * Y2_n1 * R4) * V[1] + (-5.25 * Y4_n3 + (-12.0 / 7.0) * sqrt7 * Y4_n1 + (-45.0 / 56.0) * sqrt210 * Y2_n1 * R2) * V[2] + 0.9375 * sqrt210 * Y2_n1 * V[3];
    d[26] = -Y5_n2 * Y3_p2 * V[0] + ((15.0 / 22.0) * Y6_n4 + (-6.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[1] + 1.5 * sqrt5 * Y4_n4 * V[2];
    d[27] = -Y5_n2 * Y3_p3 * V[0] + ((-5.0 / 22.0) * sqrt33 * Y6_n5 + (35.0 / 44.0) * sqrt2 * Y6_n1 + (-15.0 / 77.0) * sqrt105 * Y4_n1 * R2 + (5.0 / 28.0) * sqrt14 * Y2_n1 * R4) * V[1] + ((15.0 / 28.0) * sqrt105 * Y4_n1 + (-45.0 / 56.0) * sqrt14 * Y2_n1 * R2) * V[2] + 0.9375 * sqrt14 * Y2_n1 * V[3];
    d[23] = -Y5_n2 * Y3_n1 * V[0] + ((5.0 / 44.0) * sqrt30 * Y6_p1 + (-15.0 / 22.0) * sqrt3 * Y6_p3 + (48.0 / 77.0) * sqrt7 * Y4_p1 * R2 + (-21.0 / 11.0) * Y4_p3 * R2 + (5.0 / 28.0) * sqrt210 * Y2_p1 * R4) * V[1] + ((-12.0 / 7.0) * sqrt7 * Y4_p1 + 5.25 * Y4_p3 + (-45.0 / 56.0) * sqrt210 * Y2_p1 * R2) * V[2] + 0.9375 * sqrt210 * Y2_p1 * V[3];
    d[22] = -Y5_n2 * Y3_n2 * V[0] + ((-15.0 / 22.0) * sqrt7 * Y6_p0 + (-15.0 / 22.0) * Y6_p4 + (-30.0 / 77.0) * sqrt7 * Y4_p0 * R2 + (6.0 / 11.0) * sqrt5 * Y4_p4 * R2 + (15.0 / 14.0) * sqrt7 * Y2_p0 * R4) * V[1] + ((15.0 / 14.0) * sqrt7 * Y4_p0 - 1.5 * sqrt5 * Y4_p4 + (-135.0 / 28.0) * sqrt7 * Y2_p0 * R2) * V[2] + 5.625 * sqrt7 * Y2_p0 * V[3];
    d[21] = -Y5_n2 * Y3_n3 * V[0] + ((-35.0 / 44.0) * sqrt2 * Y6_p1 + (5.0 / 22.0) * sqrt33 * Y6_p5 + (15.0 / 77.0) * sqrt105 * Y4_p1 * R2 + (-5.0 / 28.0) * sqrt14 * Y2_p1 * R4) * V[1] + ((-15.0 / 28.0) * sqrt105 * Y4_p1 + (45.0 / 56.0) * sqrt14 * Y2_p1 * R2) * V[2] - 0.9375 * sqrt14 * Y2_p1 * V[3];
    d[17] = -Y5_n3 * Y3_p0 * V[0] + ((-5.0 / 22.0) * sqrt3 * Y6_n3 + (-18.0 / 11.0) * Y4_n3 * R2) * V[1] + 4.5 * Y4_n3 * V[2];
    d[18] = -Y5_n3 * Y3_p1 * V[0] + ((5.0 / 22.0) * sqrt15 * Y6_n4 + (35.0 / 44.0) * sqrt2 * Y6_n2 + (12.0 / 11.0) * sqrt3 * Y4_n4 * R2 + (27.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (5.0 / 14.0) * sqrt35 * Y2_n2 * R4) * V[1] + (-3.0 * sqrt3 * Y4_n4 + (-27.0 / 28.0) * sqrt21 * Y4_n2 + (-45.0 / 28.0) * sqrt35 * Y2_n2 * R2) * V[2] + 1.875 * sqrt35 * Y2_n2 * V[3];
    d[19] = -Y5_n3 * Y3_p2 * V[0] + ((5.0 / 22.0) * sqrt33 * Y6_n5 + (-10.0 / 11.0) * sqrt2 * Y6_n1 + (3.0 / 77.0) * sqrt105 * Y4_n1 * R2 + (5.0 / 7.0) * sqrt14 * Y2_n1 * R4) * V[1] + ((-3.0 / 28.0) * sqrt105 * Y4_n1 + (-45.0 / 14.0) * sqrt14 * Y2_n1 * R2) * V[2] + 3.75 * sqrt14 * Y2_n1 * V[3];
    d[20] = -Y5_n3 * Y3_p3 * V[0] + (-5.0 / 44.0) * sqrt66 * Y6_n6 * V[1];
    d[16] = -Y5_n3 * Y3_n1 * V[0] + ((35.0 / 44.0) * sqrt2 * Y6_p2 + (-5.0 / 22.0) * sqrt15 * Y6_p4 + (27.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (-12.0 / 11.0) * sqrt3 * Y4_p4 * R2 + (5.0 / 14.0) * sqrt35 * Y2_p2 * R4) * V[1] + ((-27.0 / 28.0) * sqrt21 * Y4_p2 + 3.0 * sqrt3 * Y4_p4 + (-45.0 / 28.0) * sqrt35 * Y2_p2 * R2) * V[2] + 1.875 * sqrt35 * Y2_p2 * V[3];
    d[15] = -Y5_n3 * Y3_n2 * V[0] + ((-10.0 / 11.0) * sqrt2 * Y6_p1 + (-5.0 / 22.0) * sqrt33 * Y6_p5 + (3.0 / 77.0) * sqrt105 * Y4_p1 * R2 + (5.0 / 7.0) * sqrt14 * Y2_p1 * R4) * V[1] + ((-3.0 / 28.0) * sqrt105 * Y4_p1 + (-45.0 / 14.0) * sqrt14 * Y2_p1 * R2) * V[2] + 3.75 * sqrt14 * Y2_p1 * V[3];
    d[14] = -Y5_n3 * Y3_n3 * V[0] + ((5.0 / 11.0) * sqrt7 * Y6_p0 + (5.0 / 44.0) * sqrt66 * Y6_p6 + (-90.0 / 77.0) * sqrt7 * Y4_p0 * R2 + (5.0 / 7.0) * sqrt7 * Y2_p0 * R4) * V[1] + ((45.0 / 14.0) * sqrt7 * Y4_p0 + (-45.0 / 14.0) * sqrt7 * Y2_p0 * R2) * V[2] + 3.75 * sqrt7 * Y2_p0 * V[3];
    d[10] = -Y5_n4 * Y3_p0 * V[0] + ((-15.0 / 22.0) * sqrt5 * Y6_n4 + (-36.0 / 11.0) * Y4_n4 * R2) * V[1] + 9.0 * Y4_n4 * V[2];
    d[11] = -Y5_n4 * Y3_p1 * V[0] + ((15.0 / 11.0) * Y6_n3 + (3.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[1] - 0.75 * sqrt3 * Y4_n3 * V[2];
    d[12] = -Y5_n4 * Y3_p2 * V[0] + ((15.0 / 44.0) * sqrt22 * Y6_n6 + (-15.0 / 44.0) * sqrt10 * Y6_n2 + (12.0 / 77.0) * sqrt105 * Y4_n2 * R2 + (15.0 / 14.0) * sqrt7 * Y2_n2 * R4) * V[1] + ((-3.0 / 7.0) * sqrt105 * Y4_n2 + (-135.0 / 28.0) * sqrt7 * Y2_n2 * R2) * V[2] + 5.625 * sqrt7 * Y2_n2 * V[3];
    d[13] = -Y5_n4 * Y3_p3 * V[0] + ((5.0 / 22.0) * sqrt6 * Y6_n1 + (-27.0 / 77.0) * sqrt35 * Y4_n1 * R2 + (5.0 / 14.0) * sqrt42 * Y2_n1 * R4) * V[1] + ((27.0 / 28.0) * sqrt35 * Y4_n1 + (-45.0 / 28.0) * sqrt42 * Y2_n1 * R2) * V[2] + 1.875 * sqrt42 * Y2_n1 * V[3];
    d[9] = -Y5_n4 * Y3_n1 * V[0] + ((15.0 / 11.0) * Y6_p3 + (3.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[1] - 0.75 * sqrt3 * Y4_p3 * V[2];
    d[8] = -Y5_n4 * Y3_n2 * V[0] + ((-15.0 / 44.0) * sqrt10 * Y6_p2 + (-15.0 / 44.0) * sqrt22 * Y6_p6 + (12.0 / 77.0) * sqrt105 * Y4_p2 * R2 + (15.0 / 14.0) * sqrt7 * Y2_p2 * R4) * V[1] + ((-3.0 / 7.0) * sqrt105 * Y4_p2 + (-135.0 / 28.0) * sqrt7 * Y2_p2 * R2) * V[2] + 5.625 * sqrt7 * Y2_p2 * V[3];
    d[7] = -Y5_n4 * Y3_n3 * V[0] + ((5.0 / 22.0) * sqrt6 * Y6_p1 + (-27.0 / 77.0) * sqrt35 * Y4_p1 * R2 + (5.0 / 14.0) * sqrt42 * Y2_p1 * R4) * V[1] + ((27.0 / 28.0) * sqrt35 * Y4_p1 + (-45.0 / 28.0) * sqrt42 * Y2_p1 * R2) * V[2] + 1.875 * sqrt42 * Y2_p1 * V[3];
    d[3] = -Y5_n5 * Y3_p0 * V[0] + (-15.0 / 22.0) * sqrt11 * Y6_n5 * V[1];
    d[4] = -Y5_n5 * Y3_p1 * V[0] + ((-15.0 / 44.0) * sqrt22 * Y6_n6 + (15.0 / 22.0) * sqrt3 * Y6_n4 + (-6.0 / 11.0) * sqrt15 * Y4_n4 * R2) * V[1] + 1.5 * sqrt15 * Y4_n4 * V[2];
    d[5] = -Y5_n5 * Y3_p2 * V[0] + ((-15.0 / 22.0) * Y6_n3 + (15.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[1] - 3.75 * sqrt3 * Y4_n3 * V[2];
    d[6] = -Y5_n5 * Y3_p3 * V[0] + ((5.0 / 44.0) * sqrt6 * Y6_n2 + (-45.0 / 77.0) * sqrt7 * Y4_n2 * R2 + (5.0 / 14.0) * sqrt105 * Y2_n2 * R4) * V[1] + ((45.0 / 28.0) * sqrt7 * Y4_n2 + (-45.0 / 28.0) * sqrt105 * Y2_n2 * R2) * V[2] + 1.875 * sqrt105 * Y2_n2 * V[3];
    d[2] = -Y5_n5 * Y3_n1 * V[0] + ((15.0 / 22.0) * sqrt3 * Y6_p4 + (15.0 / 44.0) * sqrt22 * Y6_p6 + (-6.0 / 11.0) * sqrt15 * Y4_p4 * R2) * V[1] + 1.5 * sqrt15 * Y4_p4 * V[2];
    d[1] = -Y5_n5 * Y3_n2 * V[0] + ((-15.0 / 22.0) * Y6_p3 + (15.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[1] - 3.75 * sqrt3 * Y4_p3 * V[2];
    d[0] = -Y5_n5 * Y3_n3 * V[0] + ((5.0 / 44.0) * sqrt6 * Y6_p2 + (-45.0 / 77.0) * sqrt7 * Y4_p2 * R2 + (5.0 / 14.0) * sqrt105 * Y2_p2 * R4) * V[1] + ((45.0 / 28.0) * sqrt7 * Y4_p2 + (-45.0 / 28.0) * sqrt105 * Y2_p2 * R2) * V[2] + 1.875 * sqrt105 * Y2_p2 * V[3];
}

}  // namespace ovlab