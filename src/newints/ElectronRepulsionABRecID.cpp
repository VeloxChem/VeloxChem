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

#include "ElectronRepulsionABRecID.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "BoysFunction.hpp"
#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace eri2cab {

auto eri2c_i_d(
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
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt165 = std::sqrt(165.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt330 = std::sqrt(330.0);
    const auto sqrt462 = std::sqrt(462.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ β^4 · p^{-6} · (s|s)^(6)
    //   V[1] ↔ α · β^5 · p^{-7} · (s|s)^(7)
    //   V[2] ↔ α^2 · β^6 · p^{-8} · (s|s)^(8)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi      = mathconst::pi_value();
    const auto two_pi52 = 2.0 * pi * pi * std::sqrt(pi);

    std::array<double, 3> V = {0.0, 0.0, 0.0};

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

            const auto p    = alpha + beta;
            const auto pinv = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto pinv7 = pinv6 * pinv;
            const auto pinv8 = pinv7 * pinv;
            const auto xi   = alpha * beta * pinv;

            const auto T = xi * R2;
            std::array<double, 9> F;
            newints::boys_function<9>(T, F);

            const auto pref     = two_pi52 / (alpha * beta * std::sqrt(p));
            const auto cab_pref = ca * cb * pref;

            V[0] += cab_pref * beta4 * pinv6 * F[6];
            V[1] += cab_pref * alpha * beta5 * pinv7 * F[7];
            V[2] += cab_pref * alpha2 * beta6 * pinv8 * F[8];
        }
    }

    // ---- Phase 3: fused M·V → 13 × 5 spherical block ----
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
    auto *d = buffer;
    d[32] = Y6_p0 * Y2_p0 * V[2] + 11.25 * Y4_p0 * V[0] + ((-21.0 / 11.0) * Y6_p0 + (-45.0 / 11.0) * Y4_p0 * R2) * V[1];
    d[33] = Y6_p0 * Y2_p1 * V[2] - 1.5 * sqrt30 * Y4_p1 * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_p1 + (6.0 / 11.0) * sqrt30 * Y4_p1 * R2) * V[1];
    d[34] = Y6_p0 * Y2_p2 * V[2] + 0.75 * sqrt15 * Y4_p2 * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_p2 + (-3.0 / 11.0) * sqrt15 * Y4_p2 * R2) * V[1];
    d[31] = Y6_p0 * Y2_n1 * V[2] - 1.5 * sqrt30 * Y4_n1 * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_n1 + (6.0 / 11.0) * sqrt30 * Y4_n1 * R2) * V[1];
    d[30] = Y6_p0 * Y2_n2 * V[2] + 0.75 * sqrt15 * Y4_n2 * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_n2 + (-3.0 / 11.0) * sqrt15 * Y4_n2 * R2) * V[1];
    d[37] = Y6_p1 * Y2_p0 * V[2] + 0.75 * sqrt210 * Y4_p1 * V[0] + ((-39.0 / 22.0) * Y6_p1 + (-3.0 / 11.0) * sqrt210 * Y4_p1 * R2) * V[1];
    d[38] = Y6_p1 * Y2_p1 * V[2] + (3.75 * sqrt7 * Y4_p0 - 0.75 * sqrt35 * Y4_p2) * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_p0 + (-3.0 / 22.0) * sqrt30 * Y6_p2 + (-15.0 / 11.0) * sqrt7 * Y4_p0 * R2 + (3.0 / 11.0) * sqrt35 * Y4_p2 * R2) * V[1];
    d[39] = Y6_p1 * Y2_p2 * V[2] + (-0.375 * sqrt70 * Y4_p1 + 0.375 * sqrt10 * Y4_p3) * V[0] + ((-21.0 / 22.0) * sqrt3 * Y6_p1 + (3.0 / 11.0) * sqrt30 * Y6_p3 + (3.0 / 22.0) * sqrt70 * Y4_p1 * R2 + (-3.0 / 22.0) * sqrt10 * Y4_p3 * R2) * V[1];
    d[36] = Y6_p1 * Y2_n1 * V[2] - 0.75 * sqrt35 * Y4_n2 * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_n2 + (3.0 / 11.0) * sqrt35 * Y4_n2 * R2) * V[1];
    d[35] = Y6_p1 * Y2_n2 * V[2] + (0.375 * sqrt10 * Y4_n3 - 0.375 * sqrt70 * Y4_n1) * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_n3 + (-21.0 / 22.0) * sqrt3 * Y6_n1 + (-3.0 / 22.0) * sqrt10 * Y4_n3 * R2 + (3.0 / 22.0) * sqrt70 * Y4_n1 * R2) * V[1];
    d[42] = Y6_p2 * Y2_p0 * V[2] + 1.5 * sqrt42 * Y4_p2 * V[0] + ((-15.0 / 11.0) * Y6_p2 + (-6.0 / 11.0) * sqrt42 * Y4_p2 * R2) * V[1];
    d[43] = Y6_p2 * Y2_p1 * V[2] + (3.0 * sqrt7 * Y4_p1 - 3.0 * Y4_p3) * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_p1 + (-15.0 / 22.0) * sqrt3 * Y6_p3 + (-12.0 / 11.0) * sqrt7 * Y4_p1 * R2 + (12.0 / 11.0) * Y4_p3 * R2) * V[1];
    d[44] = Y6_p2 * Y2_p2 * V[2] + (0.75 * sqrt70 * Y4_p0 + 0.375 * sqrt2 * Y4_p4) * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_p0 + (9.0 / 22.0) * sqrt10 * Y6_p4 + (-3.0 / 11.0) * sqrt70 * Y4_p0 * R2 + (-3.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1];
    d[41] = Y6_p2 * Y2_n1 * V[2] + (-3.0 * Y4_n3 - 3.0 * sqrt7 * Y4_n1) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_n3 + (3.0 / 22.0) * sqrt30 * Y6_n1 + (12.0 / 11.0) * Y4_n3 * R2 + (12.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1];
    d[40] = Y6_p2 * Y2_n2 * V[2] + 0.375 * sqrt2 * Y4_n4 * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_n4 + (-3.0 / 22.0) * sqrt2 * Y4_n4 * R2) * V[1];
    d[47] = Y6_p3 * Y2_p0 * V[2] + 4.5 * sqrt3 * Y4_p3 * V[0] + ((-15.0 / 22.0) * Y6_p3 + (-18.0 / 11.0) * sqrt3 * Y4_p3 * R2) * V[1];
    d[48] = Y6_p3 * Y2_p1 * V[2] + (2.25 * sqrt14 * Y4_p2 - 1.125 * sqrt2 * Y4_p4) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_p2 + (-21.0 / 44.0) * sqrt10 * Y6_p4 + (-9.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (9.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1];
    d[49] = Y6_p3 * Y2_p2 * V[2] + 2.25 * sqrt7 * Y4_p1 * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_p1 + (3.0 / 22.0) * sqrt55 * Y6_p5 + (-9.0 / 11.0) * sqrt7 * Y4_p1 * R2) * V[1];
    d[46] = Y6_p3 * Y2_n1 * V[2] + (-1.125 * sqrt2 * Y4_n4 - 2.25 * sqrt14 * Y4_n2) * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_n4 + (15.0 / 22.0) * sqrt3 * Y6_n2 + (9.0 / 22.0) * sqrt2 * Y4_n4 * R2 + (9.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[1];
    d[45] = Y6_p3 * Y2_n2 * V[2] - 2.25 * sqrt7 * Y4_n1 * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_n5 + (-3.0 / 11.0) * sqrt30 * Y6_n1 + (9.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1];
    d[52] = Y6_p4 * Y2_p0 * V[2] + 2.25 * sqrt5 * Y4_p4 * V[0] + ((3.0 / 11.0) * Y6_p4 + (-9.0 / 11.0) * sqrt5 * Y4_p4 * R2) * V[1];
    d[53] = Y6_p4 * Y2_p1 * V[2] + 1.5 * sqrt30 * Y4_p3 * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_p3 + (-9.0 / 44.0) * sqrt66 * Y6_p5 + (-6.0 / 11.0) * sqrt30 * Y4_p3 * R2) * V[1];
    d[54] = Y6_p4 * Y2_p2 * V[2] + 0.75 * sqrt105 * Y4_p2 * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_p2 + (3.0 / 22.0) * sqrt22 * Y6_p6 + (-3.0 / 11.0) * sqrt105 * Y4_p2 * R2) * V[1];
    d[51] = Y6_p4 * Y2_n1 * V[2] - 1.5 * sqrt30 * Y4_n3 * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_n5 + (21.0 / 44.0) * sqrt10 * Y6_n3 + (6.0 / 11.0) * sqrt30 * Y4_n3 * R2) * V[1];
    d[50] = Y6_p4 * Y2_n2 * V[2] - 0.75 * sqrt105 * Y4_n2 * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_n6 + (-9.0 / 22.0) * sqrt10 * Y6_n2 + (3.0 / 11.0) * sqrt105 * Y4_n2 * R2) * V[1];
    d[57] = Y6_p5 * Y2_p0 * V[2] + 1.5 * Y6_p5 * V[1];
    d[58] = Y6_p5 * Y2_p1 * V[2] + 0.375 * sqrt330 * Y4_p4 * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_p4 - 1.5 * Y6_p6 + (-3.0 / 22.0) * sqrt330 * Y4_p4 * R2) * V[1];
    d[59] = Y6_p5 * Y2_p2 * V[2] + 0.75 * sqrt165 * Y4_p3 * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_p3 + (-3.0 / 11.0) * sqrt165 * Y4_p3 * R2) * V[1];
    d[56] = Y6_p5 * Y2_n1 * V[2] - 0.375 * sqrt330 * Y4_n4 * V[0] + (-1.5 * Y6_n6 + (9.0 / 44.0) * sqrt66 * Y6_n4 + (3.0 / 22.0) * sqrt330 * Y4_n4 * R2) * V[1];
    d[55] = Y6_p5 * Y2_n2 * V[2] - 0.75 * sqrt165 * Y4_n3 * V[0] + ((-3.0 / 22.0) * sqrt55 * Y6_n3 + (3.0 / 11.0) * sqrt165 * Y4_n3 * R2) * V[1];
    d[62] = Y6_p6 * Y2_p0 * V[2] + 3.0 * Y6_p6 * V[1];
    d[63] = Y6_p6 * Y2_p1 * V[2] - 1.5 * Y6_p5 * V[1];
    d[64] = Y6_p6 * Y2_p2 * V[2] + 1.125 * sqrt110 * Y4_p4 * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_p4 + (-9.0 / 22.0) * sqrt110 * Y4_p4 * R2) * V[1];
    d[61] = Y6_p6 * Y2_n1 * V[2] + 1.5 * Y6_n5 * V[1];
    d[60] = Y6_p6 * Y2_n2 * V[2] - 1.125 * sqrt110 * Y4_n4 * V[0] + ((-3.0 / 22.0) * sqrt22 * Y6_n4 + (9.0 / 22.0) * sqrt110 * Y4_n4 * R2) * V[1];
    d[27] = Y6_n1 * Y2_p0 * V[2] + 0.75 * sqrt210 * Y4_n1 * V[0] + ((-39.0 / 22.0) * Y6_n1 + (-3.0 / 11.0) * sqrt210 * Y4_n1 * R2) * V[1];
    d[28] = Y6_n1 * Y2_p1 * V[2] - 0.75 * sqrt35 * Y4_n2 * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_n2 + (3.0 / 11.0) * sqrt35 * Y4_n2 * R2) * V[1];
    d[29] = Y6_n1 * Y2_p2 * V[2] + (0.375 * sqrt10 * Y4_n3 + 0.375 * sqrt70 * Y4_n1) * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_n3 + (21.0 / 22.0) * sqrt3 * Y6_n1 + (-3.0 / 22.0) * sqrt10 * Y4_n3 * R2 + (-3.0 / 22.0) * sqrt70 * Y4_n1 * R2) * V[1];
    d[26] = Y6_n1 * Y2_n1 * V[2] + (3.75 * sqrt7 * Y4_p0 + 0.75 * sqrt35 * Y4_p2) * V[0] + ((-3.0 / 22.0) * sqrt7 * Y6_p0 + (3.0 / 22.0) * sqrt30 * Y6_p2 + (-15.0 / 11.0) * sqrt7 * Y4_p0 * R2 + (-3.0 / 11.0) * sqrt35 * Y4_p2 * R2) * V[1];
    d[25] = Y6_n1 * Y2_n2 * V[2] + (-0.375 * sqrt70 * Y4_p1 - 0.375 * sqrt10 * Y4_p3) * V[0] + ((-21.0 / 22.0) * sqrt3 * Y6_p1 + (-3.0 / 11.0) * sqrt30 * Y6_p3 + (3.0 / 22.0) * sqrt70 * Y4_p1 * R2 + (3.0 / 22.0) * sqrt10 * Y4_p3 * R2) * V[1];
    d[22] = Y6_n2 * Y2_p0 * V[2] + 1.5 * sqrt42 * Y4_n2 * V[0] + ((-15.0 / 11.0) * Y6_n2 + (-6.0 / 11.0) * sqrt42 * Y4_n2 * R2) * V[1];
    d[23] = Y6_n2 * Y2_p1 * V[2] + (-3.0 * Y4_n3 + 3.0 * sqrt7 * Y4_n1) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_n3 + (-3.0 / 22.0) * sqrt30 * Y6_n1 + (12.0 / 11.0) * Y4_n3 * R2 + (-12.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1];
    d[24] = Y6_n2 * Y2_p2 * V[2] + 0.375 * sqrt2 * Y4_n4 * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_n4 + (-3.0 / 22.0) * sqrt2 * Y4_n4 * R2) * V[1];
    d[21] = Y6_n2 * Y2_n1 * V[2] + (3.0 * sqrt7 * Y4_p1 + 3.0 * Y4_p3) * V[0] + ((-3.0 / 22.0) * sqrt30 * Y6_p1 + (15.0 / 22.0) * sqrt3 * Y6_p3 + (-12.0 / 11.0) * sqrt7 * Y4_p1 * R2 + (-12.0 / 11.0) * Y4_p3 * R2) * V[1];
    d[20] = Y6_n2 * Y2_n2 * V[2] + (0.75 * sqrt70 * Y4_p0 - 0.375 * sqrt2 * Y4_p4) * V[0] + ((3.0 / 11.0) * sqrt70 * Y6_p0 + (-9.0 / 22.0) * sqrt10 * Y6_p4 + (-3.0 / 11.0) * sqrt70 * Y4_p0 * R2 + (3.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1];
    d[17] = Y6_n3 * Y2_p0 * V[2] + 4.5 * sqrt3 * Y4_n3 * V[0] + ((-15.0 / 22.0) * Y6_n3 + (-18.0 / 11.0) * sqrt3 * Y4_n3 * R2) * V[1];
    d[18] = Y6_n3 * Y2_p1 * V[2] + (-1.125 * sqrt2 * Y4_n4 + 2.25 * sqrt14 * Y4_n2) * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_n4 + (-15.0 / 22.0) * sqrt3 * Y6_n2 + (9.0 / 22.0) * sqrt2 * Y4_n4 * R2 + (-9.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[1];
    d[19] = Y6_n3 * Y2_p2 * V[2] + 2.25 * sqrt7 * Y4_n1 * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_n5 + (3.0 / 11.0) * sqrt30 * Y6_n1 + (-9.0 / 11.0) * sqrt7 * Y4_n1 * R2) * V[1];
    d[16] = Y6_n3 * Y2_n1 * V[2] + (2.25 * sqrt14 * Y4_p2 + 1.125 * sqrt2 * Y4_p4) * V[0] + ((-15.0 / 22.0) * sqrt3 * Y6_p2 + (21.0 / 44.0) * sqrt10 * Y6_p4 + (-9.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (-9.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[1];
    d[15] = Y6_n3 * Y2_n2 * V[2] + 2.25 * sqrt7 * Y4_p1 * V[0] + ((3.0 / 11.0) * sqrt30 * Y6_p1 + (-3.0 / 22.0) * sqrt55 * Y6_p5 + (-9.0 / 11.0) * sqrt7 * Y4_p1 * R2) * V[1];
    d[12] = Y6_n4 * Y2_p0 * V[2] + 2.25 * sqrt5 * Y4_n4 * V[0] + ((3.0 / 11.0) * Y6_n4 + (-9.0 / 11.0) * sqrt5 * Y4_n4 * R2) * V[1];
    d[13] = Y6_n4 * Y2_p1 * V[2] + 1.5 * sqrt30 * Y4_n3 * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_n5 + (-21.0 / 44.0) * sqrt10 * Y6_n3 + (-6.0 / 11.0) * sqrt30 * Y4_n3 * R2) * V[1];
    d[14] = Y6_n4 * Y2_p2 * V[2] + 0.75 * sqrt105 * Y4_n2 * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_n6 + (9.0 / 22.0) * sqrt10 * Y6_n2 + (-3.0 / 11.0) * sqrt105 * Y4_n2 * R2) * V[1];
    d[11] = Y6_n4 * Y2_n1 * V[2] + 1.5 * sqrt30 * Y4_p3 * V[0] + ((-21.0 / 44.0) * sqrt10 * Y6_p3 + (9.0 / 44.0) * sqrt66 * Y6_p5 + (-6.0 / 11.0) * sqrt30 * Y4_p3 * R2) * V[1];
    d[10] = Y6_n4 * Y2_n2 * V[2] + 0.75 * sqrt105 * Y4_p2 * V[0] + ((9.0 / 22.0) * sqrt10 * Y6_p2 + (-3.0 / 22.0) * sqrt22 * Y6_p6 + (-3.0 / 11.0) * sqrt105 * Y4_p2 * R2) * V[1];
    d[7] = Y6_n5 * Y2_p0 * V[2] + 1.5 * Y6_n5 * V[1];
    d[8] = Y6_n5 * Y2_p1 * V[2] + 0.375 * sqrt330 * Y4_n4 * V[0] + (-1.5 * Y6_n6 + (-9.0 / 44.0) * sqrt66 * Y6_n4 + (-3.0 / 22.0) * sqrt330 * Y4_n4 * R2) * V[1];
    d[9] = Y6_n5 * Y2_p2 * V[2] + 0.75 * sqrt165 * Y4_n3 * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_n3 + (-3.0 / 11.0) * sqrt165 * Y4_n3 * R2) * V[1];
    d[6] = Y6_n5 * Y2_n1 * V[2] + 0.375 * sqrt330 * Y4_p4 * V[0] + ((-9.0 / 44.0) * sqrt66 * Y6_p4 + 1.5 * Y6_p6 + (-3.0 / 22.0) * sqrt330 * Y4_p4 * R2) * V[1];
    d[5] = Y6_n5 * Y2_n2 * V[2] + 0.75 * sqrt165 * Y4_p3 * V[0] + ((3.0 / 22.0) * sqrt55 * Y6_p3 + (-3.0 / 11.0) * sqrt165 * Y4_p3 * R2) * V[1];
    d[2] = Y6_n6 * Y2_p0 * V[2] + 3.0 * Y6_n6 * V[1];
    d[3] = Y6_n6 * Y2_p1 * V[2] - 1.5 * Y6_n5 * V[1];
    d[4] = Y6_n6 * Y2_p2 * V[2] + 1.125 * sqrt110 * Y4_n4 * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_n4 + (-9.0 / 22.0) * sqrt110 * Y4_n4 * R2) * V[1];
    d[1] = Y6_n6 * Y2_n1 * V[2] - 1.5 * Y6_p5 * V[1];
    d[0] = Y6_n6 * Y2_n2 * V[2] + 1.125 * sqrt110 * Y4_p4 * V[0] + ((3.0 / 22.0) * sqrt22 * Y6_p4 + (-9.0 / 22.0) * sqrt110 * Y4_p4 * R2) * V[1];
}

}  // namespace eri2cab