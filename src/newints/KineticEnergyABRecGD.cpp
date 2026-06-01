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

#include "KineticEnergyABRecGD.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_g_d(
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
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt330 = std::sqrt(330.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^2 · β^4 · p^{-6} · (s|T|s)
    //   V[1] ↔ α · β^3 · p^{-5} · (s|T|s)
    //   V[2] ↔ β^2 · p^{-4} · (s|T|s)
    //   V[3] ↔ α^2 · β^4 · p^{-6} · ξ · (s|S|s)
    //   V[4] ↔ α · β^3 · p^{-5} · ξ · (s|S|s)
    //   V[5] ↔ β^2 · p^{-4} · ξ · (s|S|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha2 * beta4 * pinv6;
            V[1] += cab_tt * alpha * beta3 * pinv5;
            V[2] += cab_tt * beta2 * pinv4;
            V[3] += cab_ss * alpha2 * beta4 * pinv6 * xi;
            V[4] += cab_ss * alpha * beta3 * pinv5 * xi;
            V[5] += cab_ss * beta2 * pinv4 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 9 × 5 spherical block ----
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
    auto *d = buffer;
    d[22] = ((5.0 / 11.0) * Y6_p0 + (20.0 / 77.0) * Y4_p0 * R2 + (2.0 / 7.0) * Y2_p0 * R4) * V[0] + ((-10.0 / 7.0) * Y4_p0 + (-18.0 / 7.0) * Y2_p0 * R2) * V[1] + 4.5 * Y2_p0 * V[2] + ((60.0 / 11.0) * Y6_p0 + (240.0 / 77.0) * Y4_p0 * R2 + (24.0 / 7.0) * Y2_p0 * R4) * V[3] + ((-100.0 / 7.0) * Y4_p0 + (-180.0 / 7.0) * Y2_p0 * R2) * V[4] + 36.0 * Y2_p0 * V[5];
    d[23] = ((5.0 / 33.0) * sqrt7 * Y6_p1 + (1.0 / 77.0) * sqrt30 * Y4_p1 * R2 + (-4.0 / 21.0) * Y2_p1 * R4) * V[0] + ((-1.0 / 14.0) * sqrt30 * Y4_p1 + (12.0 / 7.0) * Y2_p1 * R2) * V[1] - 3.0 * Y2_p1 * V[2] + ((20.0 / 11.0) * sqrt7 * Y6_p1 + (12.0 / 77.0) * sqrt30 * Y4_p1 * R2 + (-16.0 / 7.0) * Y2_p1 * R4) * V[3] + ((-5.0 / 7.0) * sqrt30 * Y4_p1 + (120.0 / 7.0) * Y2_p1 * R2) * V[4] - 24.0 * Y2_p1 * V[5];
    d[24] = ((1.0 / 33.0) * sqrt70 * Y6_p2 + (-6.0 / 77.0) * sqrt15 * Y4_p2 * R2 + (1.0 / 21.0) * Y2_p2 * R4) * V[0] + ((3.0 / 7.0) * sqrt15 * Y4_p2 + (-3.0 / 7.0) * Y2_p2 * R2) * V[1] + 0.75 * Y2_p2 * V[2] + ((4.0 / 11.0) * sqrt70 * Y6_p2 + (-72.0 / 77.0) * sqrt15 * Y4_p2 * R2 + (4.0 / 7.0) * Y2_p2 * R4) * V[3] + ((30.0 / 7.0) * sqrt15 * Y4_p2 + (-30.0 / 7.0) * Y2_p2 * R2) * V[4] + 6.0 * Y2_p2 * V[5];
    d[21] = ((5.0 / 33.0) * sqrt7 * Y6_n1 + (1.0 / 77.0) * sqrt30 * Y4_n1 * R2 + (-4.0 / 21.0) * Y2_n1 * R4) * V[0] + ((-1.0 / 14.0) * sqrt30 * Y4_n1 + (12.0 / 7.0) * Y2_n1 * R2) * V[1] - 3.0 * Y2_n1 * V[2] + ((20.0 / 11.0) * sqrt7 * Y6_n1 + (12.0 / 77.0) * sqrt30 * Y4_n1 * R2 + (-16.0 / 7.0) * Y2_n1 * R4) * V[3] + ((-5.0 / 7.0) * sqrt30 * Y4_n1 + (120.0 / 7.0) * Y2_n1 * R2) * V[4] - 24.0 * Y2_n1 * V[5];
    d[20] = ((1.0 / 33.0) * sqrt70 * Y6_n2 + (-6.0 / 77.0) * sqrt15 * Y4_n2 * R2 + (1.0 / 21.0) * Y2_n2 * R4) * V[0] + ((3.0 / 7.0) * sqrt15 * Y4_n2 + (-3.0 / 7.0) * Y2_n2 * R2) * V[1] + 0.75 * Y2_n2 * V[2] + ((4.0 / 11.0) * sqrt70 * Y6_n2 + (-72.0 / 77.0) * sqrt15 * Y4_n2 * R2 + (4.0 / 7.0) * Y2_n2 * R4) * V[3] + ((30.0 / 7.0) * sqrt15 * Y4_n2 + (-30.0 / 7.0) * Y2_n2 * R2) * V[4] + 6.0 * Y2_n2 * V[5];
    d[27] = ((1.0 / 33.0) * sqrt210 * Y6_p1 + (17.0 / 77.0) * Y4_p1 * R2 + (1.0 / 21.0) * sqrt30 * Y2_p1 * R4) * V[0] + ((-17.0 / 14.0) * Y4_p1 + (-3.0 / 7.0) * sqrt30 * Y2_p1 * R2) * V[1] + 0.75 * sqrt30 * Y2_p1 * V[2] + ((4.0 / 11.0) * sqrt210 * Y6_p1 + (204.0 / 77.0) * Y4_p1 * R2 + (4.0 / 7.0) * sqrt30 * Y2_p1 * R4) * V[3] + ((-85.0 / 7.0) * Y4_p1 + (-30.0 / 7.0) * sqrt30 * Y2_p1 * R2) * V[4] + 6.0 * sqrt30 * Y2_p1 * V[5];
    d[28] = ((-2.0 / 33.0) * sqrt30 * Y6_p0 + (4.0 / 33.0) * sqrt7 * Y6_p2 + (1.0 / 77.0) * sqrt30 * Y4_p0 * R2 + (9.0 / 154.0) * sqrt6 * Y4_p2 * R2 + (1.0 / 21.0) * sqrt30 * Y2_p0 * R4 + (-1.0 / 42.0) * sqrt10 * Y2_p2 * R4) * V[0] + ((-1.0 / 14.0) * sqrt30 * Y4_p0 + (-9.0 / 28.0) * sqrt6 * Y4_p2 + (-3.0 / 7.0) * sqrt30 * Y2_p0 * R2 + (3.0 / 14.0) * sqrt10 * Y2_p2 * R2) * V[1] + (0.75 * sqrt30 * Y2_p0 - 0.375 * sqrt10 * Y2_p2) * V[2] + ((-8.0 / 11.0) * sqrt30 * Y6_p0 + (16.0 / 11.0) * sqrt7 * Y6_p2 + (12.0 / 77.0) * sqrt30 * Y4_p0 * R2 + (54.0 / 77.0) * sqrt6 * Y4_p2 * R2 + (4.0 / 7.0) * sqrt30 * Y2_p0 * R4 + (-2.0 / 7.0) * sqrt10 * Y2_p2 * R4) * V[3] + ((-5.0 / 7.0) * sqrt30 * Y4_p0 + (-45.0 / 14.0) * sqrt6 * Y4_p2 + (-30.0 / 7.0) * sqrt30 * Y2_p0 * R2 + (15.0 / 7.0) * sqrt10 * Y2_p2 * R2) * V[4] + (6.0 * sqrt30 * Y2_p0 - 3.0 * sqrt10 * Y2_p2) * V[5];
    d[29] = ((-1.0 / 66.0) * sqrt70 * Y6_p1 + (1.0 / 11.0) * sqrt7 * Y6_p3 + (10.0 / 77.0) * sqrt3 * Y4_p1 * R2 + (-3.0 / 77.0) * sqrt21 * Y4_p3 * R2 + (-1.0 / 42.0) * sqrt10 * Y2_p1 * R4) * V[0] + ((-5.0 / 7.0) * sqrt3 * Y4_p1 + (3.0 / 14.0) * sqrt21 * Y4_p3 + (3.0 / 14.0) * sqrt10 * Y2_p1 * R2) * V[1] - 0.375 * sqrt10 * Y2_p1 * V[2] + ((-2.0 / 11.0) * sqrt70 * Y6_p1 + (12.0 / 11.0) * sqrt7 * Y6_p3 + (120.0 / 77.0) * sqrt3 * Y4_p1 * R2 + (-36.0 / 77.0) * sqrt21 * Y4_p3 * R2 + (-2.0 / 7.0) * sqrt10 * Y2_p1 * R4) * V[3] + ((-50.0 / 7.0) * sqrt3 * Y4_p1 + (15.0 / 7.0) * sqrt21 * Y4_p3 + (15.0 / 7.0) * sqrt10 * Y2_p1 * R2) * V[4] - 3.0 * sqrt10 * Y2_p1 * V[5];
    d[26] = ((4.0 / 33.0) * sqrt7 * Y6_n2 + (9.0 / 154.0) * sqrt6 * Y4_n2 * R2 + (-1.0 / 42.0) * sqrt10 * Y2_n2 * R4) * V[0] + ((-9.0 / 28.0) * sqrt6 * Y4_n2 + (3.0 / 14.0) * sqrt10 * Y2_n2 * R2) * V[1] - 0.375 * sqrt10 * Y2_n2 * V[2] + ((16.0 / 11.0) * sqrt7 * Y6_n2 + (54.0 / 77.0) * sqrt6 * Y4_n2 * R2 + (-2.0 / 7.0) * sqrt10 * Y2_n2 * R4) * V[3] + ((-45.0 / 14.0) * sqrt6 * Y4_n2 + (15.0 / 7.0) * sqrt10 * Y2_n2 * R2) * V[4] - 3.0 * sqrt10 * Y2_n2 * V[5];
    d[25] = ((1.0 / 11.0) * sqrt7 * Y6_n3 + (-1.0 / 66.0) * sqrt70 * Y6_n1 + (-3.0 / 77.0) * sqrt21 * Y4_n3 * R2 + (10.0 / 77.0) * sqrt3 * Y4_n1 * R2 + (-1.0 / 42.0) * sqrt10 * Y2_n1 * R4) * V[0] + ((3.0 / 14.0) * sqrt21 * Y4_n3 + (-5.0 / 7.0) * sqrt3 * Y4_n1 + (3.0 / 14.0) * sqrt10 * Y2_n1 * R2) * V[1] - 0.375 * sqrt10 * Y2_n1 * V[2] + ((12.0 / 11.0) * sqrt7 * Y6_n3 + (-2.0 / 11.0) * sqrt70 * Y6_n1 + (-36.0 / 77.0) * sqrt21 * Y4_n3 * R2 + (120.0 / 77.0) * sqrt3 * Y4_n1 * R2 + (-2.0 / 7.0) * sqrt10 * Y2_n1 * R4) * V[3] + ((15.0 / 7.0) * sqrt21 * Y4_n3 + (-50.0 / 7.0) * sqrt3 * Y4_n1 + (15.0 / 7.0) * sqrt10 * Y2_n1 * R2) * V[4] - 3.0 * sqrt10 * Y2_n1 * V[5];
    d[32] = ((2.0 / 33.0) * sqrt42 * Y6_p2 + (8.0 / 77.0) * Y4_p2 * R2 + (1.0 / 21.0) * sqrt15 * Y2_p2 * R4) * V[0] + ((-4.0 / 7.0) * Y4_p2 + (-3.0 / 7.0) * sqrt15 * Y2_p2 * R2) * V[1] + 0.75 * sqrt15 * Y2_p2 * V[2] + ((8.0 / 11.0) * sqrt42 * Y6_p2 + (96.0 / 77.0) * Y4_p2 * R2 + (4.0 / 7.0) * sqrt15 * Y2_p2 * R4) * V[3] + ((-40.0 / 7.0) * Y4_p2 + (-30.0 / 7.0) * sqrt15 * Y2_p2 * R2) * V[4] + 6.0 * sqrt15 * Y2_p2 * V[5];
    d[33] = ((-1.0 / 33.0) * sqrt35 * Y6_p1 + (1.0 / 11.0) * sqrt14 * Y6_p3 + (9.0 / 154.0) * sqrt6 * Y4_p1 * R2 + (5.0 / 154.0) * sqrt42 * Y4_p3 * R2 + (2.0 / 21.0) * sqrt5 * Y2_p1 * R4) * V[0] + ((-9.0 / 28.0) * sqrt6 * Y4_p1 + (-5.0 / 28.0) * sqrt42 * Y4_p3 + (-6.0 / 7.0) * sqrt5 * Y2_p1 * R2) * V[1] + 1.5 * sqrt5 * Y2_p1 * V[2] + ((-4.0 / 11.0) * sqrt35 * Y6_p1 + (12.0 / 11.0) * sqrt14 * Y6_p3 + (54.0 / 77.0) * sqrt6 * Y4_p1 * R2 + (30.0 / 77.0) * sqrt42 * Y4_p3 * R2 + (8.0 / 7.0) * sqrt5 * Y2_p1 * R4) * V[3] + ((-45.0 / 14.0) * sqrt6 * Y4_p1 + (-25.0 / 14.0) * sqrt42 * Y4_p3 + (-60.0 / 7.0) * sqrt5 * Y2_p1 * R2) * V[4] + 12.0 * sqrt5 * Y2_p1 * V[5];
    d[34] = ((1.0 / 33.0) * sqrt15 * Y6_p0 + (1.0 / 33.0) * sqrt105 * Y6_p4 + (-6.0 / 77.0) * sqrt15 * Y4_p0 * R2 + (-2.0 / 77.0) * sqrt21 * Y4_p4 * R2 + (1.0 / 21.0) * sqrt15 * Y2_p0 * R4) * V[0] + ((3.0 / 7.0) * sqrt15 * Y4_p0 + (1.0 / 7.0) * sqrt21 * Y4_p4 + (-3.0 / 7.0) * sqrt15 * Y2_p0 * R2) * V[1] + 0.75 * sqrt15 * Y2_p0 * V[2] + ((4.0 / 11.0) * sqrt15 * Y6_p0 + (4.0 / 11.0) * sqrt105 * Y6_p4 + (-72.0 / 77.0) * sqrt15 * Y4_p0 * R2 + (-24.0 / 77.0) * sqrt21 * Y4_p4 * R2 + (4.0 / 7.0) * sqrt15 * Y2_p0 * R4) * V[3] + ((30.0 / 7.0) * sqrt15 * Y4_p0 + (10.0 / 7.0) * sqrt21 * Y4_p4 + (-30.0 / 7.0) * sqrt15 * Y2_p0 * R2) * V[4] + 6.0 * sqrt15 * Y2_p0 * V[5];
    d[31] = ((1.0 / 11.0) * sqrt14 * Y6_n3 + (1.0 / 33.0) * sqrt35 * Y6_n1 + (5.0 / 154.0) * sqrt42 * Y4_n3 * R2 + (-9.0 / 154.0) * sqrt6 * Y4_n1 * R2 + (-2.0 / 21.0) * sqrt5 * Y2_n1 * R4) * V[0] + ((-5.0 / 28.0) * sqrt42 * Y4_n3 + (9.0 / 28.0) * sqrt6 * Y4_n1 + (6.0 / 7.0) * sqrt5 * Y2_n1 * R2) * V[1] - 1.5 * sqrt5 * Y2_n1 * V[2] + ((12.0 / 11.0) * sqrt14 * Y6_n3 + (4.0 / 11.0) * sqrt35 * Y6_n1 + (30.0 / 77.0) * sqrt42 * Y4_n3 * R2 + (-54.0 / 77.0) * sqrt6 * Y4_n1 * R2 + (-8.0 / 7.0) * sqrt5 * Y2_n1 * R4) * V[3] + ((-25.0 / 14.0) * sqrt42 * Y4_n3 + (45.0 / 14.0) * sqrt6 * Y4_n1 + (60.0 / 7.0) * sqrt5 * Y2_n1 * R2) * V[4] - 12.0 * sqrt5 * Y2_n1 * V[5];
    d[30] = ((1.0 / 33.0) * sqrt105 * Y6_n4 + (-2.0 / 77.0) * sqrt21 * Y4_n4 * R2) * V[0] + (1.0 / 7.0) * sqrt21 * Y4_n4 * V[1] + ((4.0 / 11.0) * sqrt105 * Y6_n4 + (-24.0 / 77.0) * sqrt21 * Y4_n4 * R2) * V[3] + (10.0 / 7.0) * sqrt21 * Y4_n4 * V[4];
    d[37] = ((2.0 / 11.0) * sqrt3 * Y6_p3 + (-1.0 / 11.0) * Y4_p3 * R2) * V[0] + 0.5 * Y4_p3 * V[1] + ((24.0 / 11.0) * sqrt3 * Y6_p3 + (-12.0 / 11.0) * Y4_p3 * R2) * V[3] + 5.0 * Y4_p3 * V[4];
    d[38] = ((-4.0 / 33.0) * Y6_p2 + (2.0 / 33.0) * sqrt30 * Y6_p4 + (5.0 / 154.0) * sqrt42 * Y4_p2 * R2 + (1.0 / 11.0) * sqrt6 * Y4_p4 * R2 + (1.0 / 42.0) * sqrt70 * Y2_p2 * R4) * V[0] + ((-5.0 / 28.0) * sqrt42 * Y4_p2 - 0.5 * sqrt6 * Y4_p4 + (-3.0 / 14.0) * sqrt70 * Y2_p2 * R2) * V[1] + 0.375 * sqrt70 * Y2_p2 * V[2] + ((-16.0 / 11.0) * Y6_p2 + (8.0 / 11.0) * sqrt30 * Y6_p4 + (30.0 / 77.0) * sqrt42 * Y4_p2 * R2 + (12.0 / 11.0) * sqrt6 * Y4_p4 * R2 + (2.0 / 7.0) * sqrt70 * Y2_p2 * R4) * V[3] + ((-25.0 / 14.0) * sqrt42 * Y4_p2 - 5.0 * sqrt6 * Y4_p4 + (-15.0 / 7.0) * sqrt70 * Y2_p2 * R2) * V[4] + 3.0 * sqrt70 * Y2_p2 * V[5];
    d[39] = ((1.0 / 66.0) * sqrt10 * Y6_p1 + (1.0 / 33.0) * sqrt165 * Y6_p5 + (-3.0 / 77.0) * sqrt21 * Y4_p1 * R2 + (1.0 / 42.0) * sqrt70 * Y2_p1 * R4) * V[0] + ((3.0 / 14.0) * sqrt21 * Y4_p1 + (-3.0 / 14.0) * sqrt70 * Y2_p1 * R2) * V[1] + 0.375 * sqrt70 * Y2_p1 * V[2] + ((2.0 / 11.0) * sqrt10 * Y6_p1 + (4.0 / 11.0) * sqrt165 * Y6_p5 + (-36.0 / 77.0) * sqrt21 * Y4_p1 * R2 + (2.0 / 7.0) * sqrt70 * Y2_p1 * R4) * V[3] + ((15.0 / 7.0) * sqrt21 * Y4_p1 + (-15.0 / 7.0) * sqrt70 * Y2_p1 * R2) * V[4] + 3.0 * sqrt70 * Y2_p1 * V[5];
    d[36] = ((2.0 / 33.0) * sqrt30 * Y6_n4 + (4.0 / 33.0) * Y6_n2 + (1.0 / 11.0) * sqrt6 * Y4_n4 * R2 + (-5.0 / 154.0) * sqrt42 * Y4_n2 * R2 + (-1.0 / 42.0) * sqrt70 * Y2_n2 * R4) * V[0] + (-0.5 * sqrt6 * Y4_n4 + (5.0 / 28.0) * sqrt42 * Y4_n2 + (3.0 / 14.0) * sqrt70 * Y2_n2 * R2) * V[1] - 0.375 * sqrt70 * Y2_n2 * V[2] + ((8.0 / 11.0) * sqrt30 * Y6_n4 + (16.0 / 11.0) * Y6_n2 + (12.0 / 11.0) * sqrt6 * Y4_n4 * R2 + (-30.0 / 77.0) * sqrt42 * Y4_n2 * R2 + (-2.0 / 7.0) * sqrt70 * Y2_n2 * R4) * V[3] + (-5.0 * sqrt6 * Y4_n4 + (25.0 / 14.0) * sqrt42 * Y4_n2 + (15.0 / 7.0) * sqrt70 * Y2_n2 * R2) * V[4] - 3.0 * sqrt70 * Y2_n2 * V[5];
    d[35] = ((1.0 / 33.0) * sqrt165 * Y6_n5 + (-1.0 / 66.0) * sqrt10 * Y6_n1 + (3.0 / 77.0) * sqrt21 * Y4_n1 * R2 + (-1.0 / 42.0) * sqrt70 * Y2_n1 * R4) * V[0] + ((-3.0 / 14.0) * sqrt21 * Y4_n1 + (3.0 / 14.0) * sqrt70 * Y2_n1 * R2) * V[1] - 0.375 * sqrt70 * Y2_n1 * V[2] + ((4.0 / 11.0) * sqrt165 * Y6_n5 + (-2.0 / 11.0) * sqrt10 * Y6_n1 + (36.0 / 77.0) * sqrt21 * Y4_n1 * R2 + (-2.0 / 7.0) * sqrt70 * Y2_n1 * R4) * V[3] + ((-15.0 / 7.0) * sqrt21 * Y4_n1 + (15.0 / 7.0) * sqrt70 * Y2_n1 * R2) * V[4] - 3.0 * sqrt70 * Y2_n1 * V[5];
    d[42] = ((1.0 / 11.0) * sqrt5 * Y6_p4 + (-4.0 / 11.0) * Y4_p4 * R2) * V[0] + 2.0 * Y4_p4 * V[1] + ((12.0 / 11.0) * sqrt5 * Y6_p4 + (-48.0 / 11.0) * Y4_p4 * R2) * V[3] + 20.0 * Y4_p4 * V[4];
    d[43] = ((-1.0 / 22.0) * sqrt2 * Y6_p3 + (1.0 / 66.0) * sqrt330 * Y6_p5 + (1.0 / 11.0) * sqrt6 * Y4_p3 * R2) * V[0] - 0.5 * sqrt6 * Y4_p3 * V[1] + ((-6.0 / 11.0) * sqrt2 * Y6_p3 + (2.0 / 11.0) * sqrt330 * Y6_p5 + (12.0 / 11.0) * sqrt6 * Y4_p3 * R2) * V[3] - 5.0 * sqrt6 * Y4_p3 * V[4];
    d[44] = ((1.0 / 66.0) * sqrt2 * Y6_p2 + (1.0 / 22.0) * sqrt110 * Y6_p6 + (-2.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (1.0 / 21.0) * sqrt35 * Y2_p2 * R4) * V[0] + ((1.0 / 7.0) * sqrt21 * Y4_p2 + (-3.0 / 7.0) * sqrt35 * Y2_p2 * R2) * V[1] + 0.75 * sqrt35 * Y2_p2 * V[2] + ((2.0 / 11.0) * sqrt2 * Y6_p2 + (6.0 / 11.0) * sqrt110 * Y6_p6 + (-24.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (4.0 / 7.0) * sqrt35 * Y2_p2 * R4) * V[3] + ((10.0 / 7.0) * sqrt21 * Y4_p2 + (-30.0 / 7.0) * sqrt35 * Y2_p2 * R2) * V[4] + 6.0 * sqrt35 * Y2_p2 * V[5];
    d[41] = ((1.0 / 66.0) * sqrt330 * Y6_n5 + (1.0 / 22.0) * sqrt2 * Y6_n3 + (-1.0 / 11.0) * sqrt6 * Y4_n3 * R2) * V[0] + 0.5 * sqrt6 * Y4_n3 * V[1] + ((2.0 / 11.0) * sqrt330 * Y6_n5 + (6.0 / 11.0) * sqrt2 * Y6_n3 + (-12.0 / 11.0) * sqrt6 * Y4_n3 * R2) * V[3] + 5.0 * sqrt6 * Y4_n3 * V[4];
    d[40] = ((1.0 / 22.0) * sqrt110 * Y6_n6 + (-1.0 / 66.0) * sqrt2 * Y6_n2 + (2.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (-1.0 / 21.0) * sqrt35 * Y2_n2 * R4) * V[0] + ((-1.0 / 7.0) * sqrt21 * Y4_n2 + (3.0 / 7.0) * sqrt35 * Y2_n2 * R2) * V[1] - 0.75 * sqrt35 * Y2_n2 * V[2] + ((6.0 / 11.0) * sqrt110 * Y6_n6 + (-2.0 / 11.0) * sqrt2 * Y6_n2 + (24.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (-4.0 / 7.0) * sqrt35 * Y2_n2 * R4) * V[3] + ((-10.0 / 7.0) * sqrt21 * Y4_n2 + (30.0 / 7.0) * sqrt35 * Y2_n2 * R2) * V[4] - 6.0 * sqrt35 * Y2_n2 * V[5];
    d[17] = ((1.0 / 33.0) * sqrt210 * Y6_n1 + (17.0 / 77.0) * Y4_n1 * R2 + (1.0 / 21.0) * sqrt30 * Y2_n1 * R4) * V[0] + ((-17.0 / 14.0) * Y4_n1 + (-3.0 / 7.0) * sqrt30 * Y2_n1 * R2) * V[1] + 0.75 * sqrt30 * Y2_n1 * V[2] + ((4.0 / 11.0) * sqrt210 * Y6_n1 + (204.0 / 77.0) * Y4_n1 * R2 + (4.0 / 7.0) * sqrt30 * Y2_n1 * R4) * V[3] + ((-85.0 / 7.0) * Y4_n1 + (-30.0 / 7.0) * sqrt30 * Y2_n1 * R2) * V[4] + 6.0 * sqrt30 * Y2_n1 * V[5];
    d[18] = ((4.0 / 33.0) * sqrt7 * Y6_n2 + (9.0 / 154.0) * sqrt6 * Y4_n2 * R2 + (-1.0 / 42.0) * sqrt10 * Y2_n2 * R4) * V[0] + ((-9.0 / 28.0) * sqrt6 * Y4_n2 + (3.0 / 14.0) * sqrt10 * Y2_n2 * R2) * V[1] - 0.375 * sqrt10 * Y2_n2 * V[2] + ((16.0 / 11.0) * sqrt7 * Y6_n2 + (54.0 / 77.0) * sqrt6 * Y4_n2 * R2 + (-2.0 / 7.0) * sqrt10 * Y2_n2 * R4) * V[3] + ((-45.0 / 14.0) * sqrt6 * Y4_n2 + (15.0 / 7.0) * sqrt10 * Y2_n2 * R2) * V[4] - 3.0 * sqrt10 * Y2_n2 * V[5];
    d[19] = ((1.0 / 11.0) * sqrt7 * Y6_n3 + (1.0 / 66.0) * sqrt70 * Y6_n1 + (-3.0 / 77.0) * sqrt21 * Y4_n3 * R2 + (-10.0 / 77.0) * sqrt3 * Y4_n1 * R2 + (1.0 / 42.0) * sqrt10 * Y2_n1 * R4) * V[0] + ((3.0 / 14.0) * sqrt21 * Y4_n3 + (5.0 / 7.0) * sqrt3 * Y4_n1 + (-3.0 / 14.0) * sqrt10 * Y2_n1 * R2) * V[1] + 0.375 * sqrt10 * Y2_n1 * V[2] + ((12.0 / 11.0) * sqrt7 * Y6_n3 + (2.0 / 11.0) * sqrt70 * Y6_n1 + (-36.0 / 77.0) * sqrt21 * Y4_n3 * R2 + (-120.0 / 77.0) * sqrt3 * Y4_n1 * R2 + (2.0 / 7.0) * sqrt10 * Y2_n1 * R4) * V[3] + ((15.0 / 7.0) * sqrt21 * Y4_n3 + (50.0 / 7.0) * sqrt3 * Y4_n1 + (-15.0 / 7.0) * sqrt10 * Y2_n1 * R2) * V[4] + 3.0 * sqrt10 * Y2_n1 * V[5];
    d[16] = ((-2.0 / 33.0) * sqrt30 * Y6_p0 + (-4.0 / 33.0) * sqrt7 * Y6_p2 + (1.0 / 77.0) * sqrt30 * Y4_p0 * R2 + (-9.0 / 154.0) * sqrt6 * Y4_p2 * R2 + (1.0 / 21.0) * sqrt30 * Y2_p0 * R4 + (1.0 / 42.0) * sqrt10 * Y2_p2 * R4) * V[0] + ((-1.0 / 14.0) * sqrt30 * Y4_p0 + (9.0 / 28.0) * sqrt6 * Y4_p2 + (-3.0 / 7.0) * sqrt30 * Y2_p0 * R2 + (-3.0 / 14.0) * sqrt10 * Y2_p2 * R2) * V[1] + (0.75 * sqrt30 * Y2_p0 + 0.375 * sqrt10 * Y2_p2) * V[2] + ((-8.0 / 11.0) * sqrt30 * Y6_p0 + (-16.0 / 11.0) * sqrt7 * Y6_p2 + (12.0 / 77.0) * sqrt30 * Y4_p0 * R2 + (-54.0 / 77.0) * sqrt6 * Y4_p2 * R2 + (4.0 / 7.0) * sqrt30 * Y2_p0 * R4 + (2.0 / 7.0) * sqrt10 * Y2_p2 * R4) * V[3] + ((-5.0 / 7.0) * sqrt30 * Y4_p0 + (45.0 / 14.0) * sqrt6 * Y4_p2 + (-30.0 / 7.0) * sqrt30 * Y2_p0 * R2 + (-15.0 / 7.0) * sqrt10 * Y2_p2 * R2) * V[4] + (6.0 * sqrt30 * Y2_p0 + 3.0 * sqrt10 * Y2_p2) * V[5];
    d[15] = ((-1.0 / 66.0) * sqrt70 * Y6_p1 + (-1.0 / 11.0) * sqrt7 * Y6_p3 + (10.0 / 77.0) * sqrt3 * Y4_p1 * R2 + (3.0 / 77.0) * sqrt21 * Y4_p3 * R2 + (-1.0 / 42.0) * sqrt10 * Y2_p1 * R4) * V[0] + ((-5.0 / 7.0) * sqrt3 * Y4_p1 + (-3.0 / 14.0) * sqrt21 * Y4_p3 + (3.0 / 14.0) * sqrt10 * Y2_p1 * R2) * V[1] - 0.375 * sqrt10 * Y2_p1 * V[2] + ((-2.0 / 11.0) * sqrt70 * Y6_p1 + (-12.0 / 11.0) * sqrt7 * Y6_p3 + (120.0 / 77.0) * sqrt3 * Y4_p1 * R2 + (36.0 / 77.0) * sqrt21 * Y4_p3 * R2 + (-2.0 / 7.0) * sqrt10 * Y2_p1 * R4) * V[3] + ((-50.0 / 7.0) * sqrt3 * Y4_p1 + (-15.0 / 7.0) * sqrt21 * Y4_p3 + (15.0 / 7.0) * sqrt10 * Y2_p1 * R2) * V[4] - 3.0 * sqrt10 * Y2_p1 * V[5];
    d[12] = ((2.0 / 33.0) * sqrt42 * Y6_n2 + (8.0 / 77.0) * Y4_n2 * R2 + (1.0 / 21.0) * sqrt15 * Y2_n2 * R4) * V[0] + ((-4.0 / 7.0) * Y4_n2 + (-3.0 / 7.0) * sqrt15 * Y2_n2 * R2) * V[1] + 0.75 * sqrt15 * Y2_n2 * V[2] + ((8.0 / 11.0) * sqrt42 * Y6_n2 + (96.0 / 77.0) * Y4_n2 * R2 + (4.0 / 7.0) * sqrt15 * Y2_n2 * R4) * V[3] + ((-40.0 / 7.0) * Y4_n2 + (-30.0 / 7.0) * sqrt15 * Y2_n2 * R2) * V[4] + 6.0 * sqrt15 * Y2_n2 * V[5];
    d[13] = ((1.0 / 11.0) * sqrt14 * Y6_n3 + (-1.0 / 33.0) * sqrt35 * Y6_n1 + (5.0 / 154.0) * sqrt42 * Y4_n3 * R2 + (9.0 / 154.0) * sqrt6 * Y4_n1 * R2 + (2.0 / 21.0) * sqrt5 * Y2_n1 * R4) * V[0] + ((-5.0 / 28.0) * sqrt42 * Y4_n3 + (-9.0 / 28.0) * sqrt6 * Y4_n1 + (-6.0 / 7.0) * sqrt5 * Y2_n1 * R2) * V[1] + 1.5 * sqrt5 * Y2_n1 * V[2] + ((12.0 / 11.0) * sqrt14 * Y6_n3 + (-4.0 / 11.0) * sqrt35 * Y6_n1 + (30.0 / 77.0) * sqrt42 * Y4_n3 * R2 + (54.0 / 77.0) * sqrt6 * Y4_n1 * R2 + (8.0 / 7.0) * sqrt5 * Y2_n1 * R4) * V[3] + ((-25.0 / 14.0) * sqrt42 * Y4_n3 + (-45.0 / 14.0) * sqrt6 * Y4_n1 + (-60.0 / 7.0) * sqrt5 * Y2_n1 * R2) * V[4] + 12.0 * sqrt5 * Y2_n1 * V[5];
    d[14] = ((1.0 / 33.0) * sqrt105 * Y6_n4 + (-2.0 / 77.0) * sqrt21 * Y4_n4 * R2) * V[0] + (1.0 / 7.0) * sqrt21 * Y4_n4 * V[1] + ((4.0 / 11.0) * sqrt105 * Y6_n4 + (-24.0 / 77.0) * sqrt21 * Y4_n4 * R2) * V[3] + (10.0 / 7.0) * sqrt21 * Y4_n4 * V[4];
    d[11] = ((-1.0 / 33.0) * sqrt35 * Y6_p1 + (-1.0 / 11.0) * sqrt14 * Y6_p3 + (9.0 / 154.0) * sqrt6 * Y4_p1 * R2 + (-5.0 / 154.0) * sqrt42 * Y4_p3 * R2 + (2.0 / 21.0) * sqrt5 * Y2_p1 * R4) * V[0] + ((-9.0 / 28.0) * sqrt6 * Y4_p1 + (5.0 / 28.0) * sqrt42 * Y4_p3 + (-6.0 / 7.0) * sqrt5 * Y2_p1 * R2) * V[1] + 1.5 * sqrt5 * Y2_p1 * V[2] + ((-4.0 / 11.0) * sqrt35 * Y6_p1 + (-12.0 / 11.0) * sqrt14 * Y6_p3 + (54.0 / 77.0) * sqrt6 * Y4_p1 * R2 + (-30.0 / 77.0) * sqrt42 * Y4_p3 * R2 + (8.0 / 7.0) * sqrt5 * Y2_p1 * R4) * V[3] + ((-45.0 / 14.0) * sqrt6 * Y4_p1 + (25.0 / 14.0) * sqrt42 * Y4_p3 + (-60.0 / 7.0) * sqrt5 * Y2_p1 * R2) * V[4] + 12.0 * sqrt5 * Y2_p1 * V[5];
    d[10] = ((1.0 / 33.0) * sqrt15 * Y6_p0 + (-1.0 / 33.0) * sqrt105 * Y6_p4 + (-6.0 / 77.0) * sqrt15 * Y4_p0 * R2 + (2.0 / 77.0) * sqrt21 * Y4_p4 * R2 + (1.0 / 21.0) * sqrt15 * Y2_p0 * R4) * V[0] + ((3.0 / 7.0) * sqrt15 * Y4_p0 + (-1.0 / 7.0) * sqrt21 * Y4_p4 + (-3.0 / 7.0) * sqrt15 * Y2_p0 * R2) * V[1] + 0.75 * sqrt15 * Y2_p0 * V[2] + ((4.0 / 11.0) * sqrt15 * Y6_p0 + (-4.0 / 11.0) * sqrt105 * Y6_p4 + (-72.0 / 77.0) * sqrt15 * Y4_p0 * R2 + (24.0 / 77.0) * sqrt21 * Y4_p4 * R2 + (4.0 / 7.0) * sqrt15 * Y2_p0 * R4) * V[3] + ((30.0 / 7.0) * sqrt15 * Y4_p0 + (-10.0 / 7.0) * sqrt21 * Y4_p4 + (-30.0 / 7.0) * sqrt15 * Y2_p0 * R2) * V[4] + 6.0 * sqrt15 * Y2_p0 * V[5];
    d[7] = ((2.0 / 11.0) * sqrt3 * Y6_n3 + (-1.0 / 11.0) * Y4_n3 * R2) * V[0] + 0.5 * Y4_n3 * V[1] + ((24.0 / 11.0) * sqrt3 * Y6_n3 + (-12.0 / 11.0) * Y4_n3 * R2) * V[3] + 5.0 * Y4_n3 * V[4];
    d[8] = ((2.0 / 33.0) * sqrt30 * Y6_n4 + (-4.0 / 33.0) * Y6_n2 + (1.0 / 11.0) * sqrt6 * Y4_n4 * R2 + (5.0 / 154.0) * sqrt42 * Y4_n2 * R2 + (1.0 / 42.0) * sqrt70 * Y2_n2 * R4) * V[0] + (-0.5 * sqrt6 * Y4_n4 + (-5.0 / 28.0) * sqrt42 * Y4_n2 + (-3.0 / 14.0) * sqrt70 * Y2_n2 * R2) * V[1] + 0.375 * sqrt70 * Y2_n2 * V[2] + ((8.0 / 11.0) * sqrt30 * Y6_n4 + (-16.0 / 11.0) * Y6_n2 + (12.0 / 11.0) * sqrt6 * Y4_n4 * R2 + (30.0 / 77.0) * sqrt42 * Y4_n2 * R2 + (2.0 / 7.0) * sqrt70 * Y2_n2 * R4) * V[3] + (-5.0 * sqrt6 * Y4_n4 + (-25.0 / 14.0) * sqrt42 * Y4_n2 + (-15.0 / 7.0) * sqrt70 * Y2_n2 * R2) * V[4] + 3.0 * sqrt70 * Y2_n2 * V[5];
    d[9] = ((1.0 / 33.0) * sqrt165 * Y6_n5 + (1.0 / 66.0) * sqrt10 * Y6_n1 + (-3.0 / 77.0) * sqrt21 * Y4_n1 * R2 + (1.0 / 42.0) * sqrt70 * Y2_n1 * R4) * V[0] + ((3.0 / 14.0) * sqrt21 * Y4_n1 + (-3.0 / 14.0) * sqrt70 * Y2_n1 * R2) * V[1] + 0.375 * sqrt70 * Y2_n1 * V[2] + ((4.0 / 11.0) * sqrt165 * Y6_n5 + (2.0 / 11.0) * sqrt10 * Y6_n1 + (-36.0 / 77.0) * sqrt21 * Y4_n1 * R2 + (2.0 / 7.0) * sqrt70 * Y2_n1 * R4) * V[3] + ((15.0 / 7.0) * sqrt21 * Y4_n1 + (-15.0 / 7.0) * sqrt70 * Y2_n1 * R2) * V[4] + 3.0 * sqrt70 * Y2_n1 * V[5];
    d[6] = ((-4.0 / 33.0) * Y6_p2 + (-2.0 / 33.0) * sqrt30 * Y6_p4 + (5.0 / 154.0) * sqrt42 * Y4_p2 * R2 + (-1.0 / 11.0) * sqrt6 * Y4_p4 * R2 + (1.0 / 42.0) * sqrt70 * Y2_p2 * R4) * V[0] + ((-5.0 / 28.0) * sqrt42 * Y4_p2 + 0.5 * sqrt6 * Y4_p4 + (-3.0 / 14.0) * sqrt70 * Y2_p2 * R2) * V[1] + 0.375 * sqrt70 * Y2_p2 * V[2] + ((-16.0 / 11.0) * Y6_p2 + (-8.0 / 11.0) * sqrt30 * Y6_p4 + (30.0 / 77.0) * sqrt42 * Y4_p2 * R2 + (-12.0 / 11.0) * sqrt6 * Y4_p4 * R2 + (2.0 / 7.0) * sqrt70 * Y2_p2 * R4) * V[3] + ((-25.0 / 14.0) * sqrt42 * Y4_p2 + 5.0 * sqrt6 * Y4_p4 + (-15.0 / 7.0) * sqrt70 * Y2_p2 * R2) * V[4] + 3.0 * sqrt70 * Y2_p2 * V[5];
    d[5] = ((1.0 / 66.0) * sqrt10 * Y6_p1 + (-1.0 / 33.0) * sqrt165 * Y6_p5 + (-3.0 / 77.0) * sqrt21 * Y4_p1 * R2 + (1.0 / 42.0) * sqrt70 * Y2_p1 * R4) * V[0] + ((3.0 / 14.0) * sqrt21 * Y4_p1 + (-3.0 / 14.0) * sqrt70 * Y2_p1 * R2) * V[1] + 0.375 * sqrt70 * Y2_p1 * V[2] + ((2.0 / 11.0) * sqrt10 * Y6_p1 + (-4.0 / 11.0) * sqrt165 * Y6_p5 + (-36.0 / 77.0) * sqrt21 * Y4_p1 * R2 + (2.0 / 7.0) * sqrt70 * Y2_p1 * R4) * V[3] + ((15.0 / 7.0) * sqrt21 * Y4_p1 + (-15.0 / 7.0) * sqrt70 * Y2_p1 * R2) * V[4] + 3.0 * sqrt70 * Y2_p1 * V[5];
    d[2] = ((1.0 / 11.0) * sqrt5 * Y6_n4 + (-4.0 / 11.0) * Y4_n4 * R2) * V[0] + 2.0 * Y4_n4 * V[1] + ((12.0 / 11.0) * sqrt5 * Y6_n4 + (-48.0 / 11.0) * Y4_n4 * R2) * V[3] + 20.0 * Y4_n4 * V[4];
    d[3] = ((1.0 / 66.0) * sqrt330 * Y6_n5 + (-1.0 / 22.0) * sqrt2 * Y6_n3 + (1.0 / 11.0) * sqrt6 * Y4_n3 * R2) * V[0] - 0.5 * sqrt6 * Y4_n3 * V[1] + ((2.0 / 11.0) * sqrt330 * Y6_n5 + (-6.0 / 11.0) * sqrt2 * Y6_n3 + (12.0 / 11.0) * sqrt6 * Y4_n3 * R2) * V[3] - 5.0 * sqrt6 * Y4_n3 * V[4];
    d[4] = ((1.0 / 22.0) * sqrt110 * Y6_n6 + (1.0 / 66.0) * sqrt2 * Y6_n2 + (-2.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (1.0 / 21.0) * sqrt35 * Y2_n2 * R4) * V[0] + ((1.0 / 7.0) * sqrt21 * Y4_n2 + (-3.0 / 7.0) * sqrt35 * Y2_n2 * R2) * V[1] + 0.75 * sqrt35 * Y2_n2 * V[2] + ((6.0 / 11.0) * sqrt110 * Y6_n6 + (2.0 / 11.0) * sqrt2 * Y6_n2 + (-24.0 / 77.0) * sqrt21 * Y4_n2 * R2 + (4.0 / 7.0) * sqrt35 * Y2_n2 * R4) * V[3] + ((10.0 / 7.0) * sqrt21 * Y4_n2 + (-30.0 / 7.0) * sqrt35 * Y2_n2 * R2) * V[4] + 6.0 * sqrt35 * Y2_n2 * V[5];
    d[1] = ((-1.0 / 22.0) * sqrt2 * Y6_p3 + (-1.0 / 66.0) * sqrt330 * Y6_p5 + (1.0 / 11.0) * sqrt6 * Y4_p3 * R2) * V[0] - 0.5 * sqrt6 * Y4_p3 * V[1] + ((-6.0 / 11.0) * sqrt2 * Y6_p3 + (-2.0 / 11.0) * sqrt330 * Y6_p5 + (12.0 / 11.0) * sqrt6 * Y4_p3 * R2) * V[3] - 5.0 * sqrt6 * Y4_p3 * V[4];
    d[0] = ((1.0 / 66.0) * sqrt2 * Y6_p2 + (-1.0 / 22.0) * sqrt110 * Y6_p6 + (-2.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (1.0 / 21.0) * sqrt35 * Y2_p2 * R4) * V[0] + ((1.0 / 7.0) * sqrt21 * Y4_p2 + (-3.0 / 7.0) * sqrt35 * Y2_p2 * R2) * V[1] + 0.75 * sqrt35 * Y2_p2 * V[2] + ((2.0 / 11.0) * sqrt2 * Y6_p2 + (-6.0 / 11.0) * sqrt110 * Y6_p6 + (-24.0 / 77.0) * sqrt21 * Y4_p2 * R2 + (4.0 / 7.0) * sqrt35 * Y2_p2 * R4) * V[3] + ((10.0 / 7.0) * sqrt21 * Y4_p2 + (-30.0 / 7.0) * sqrt35 * Y2_p2 * R2) * V[4] + 6.0 * sqrt35 * Y2_p2 * V[5];
}

}  // namespace kinab