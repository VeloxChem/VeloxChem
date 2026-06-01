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

#include "OverlapABRecGF.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_g_f(
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
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt210 = std::sqrt(210.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^3 · β^4 · p^{-7} · (s|s)
    //   V[1] ↔ α^2 · β^3 · p^{-6} · (s|s)
    //   V[2] ↔ α · β^2 · p^{-5} · (s|s)
    //   V[3] ↔ β · p^{-4} · (s|s)
    const auto exps_a  = bra.get_exponents();
    const auto coefs_a = bra.get_normalization_factors();
    const auto exps_b  = ket.get_exponents();
    const auto coefs_b = ket.get_normalization_factors();

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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto pinv7 = pinv6 * pinv;

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha3 * beta4 * pinv7;
            V[1] += cab_ss * alpha2 * beta3 * pinv6;
            V[2] += cab_ss * alpha * beta2 * pinv5;
            V[3] += cab_ss * beta * pinv4;
        }
    }

    // ---- Phase 3: fused M·V → 9 × 7 spherical block ----
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
    const auto R4 = R2 * R2;
    newints::Block out{9, 7, std::vector<double>(63, 0.0)};
    auto *d = out.data.data();
    d[31] = Y4_p0 * Y3_p0 * V[0] + ((-10.0 / 7.0) * Y5_p0 - 2.0 * Y3_p0 * R2 + (-18.0 / 7.0) * Y1_p0 * R4) * V[1] + (4.5 * Y3_p0 + 9.0 * Y1_p0 * R2) * V[2] - 7.5 * Y1_p0 * V[3];
    d[32] = Y4_p0 * Y3_p1 * V[0] + ((-11.0 / 42.0) * sqrt10 * Y5_p1 + (-1.0 / 3.0) * Y3_p1 * R2 + (9.0 / 14.0) * sqrt6 * Y1_p1 * R4) * V[1] + (0.75 * Y3_p1 - 2.25 * sqrt6 * Y1_p1 * R2) * V[2] + 1.875 * sqrt6 * Y1_p1 * V[3];
    d[33] = Y4_p0 * Y3_p2 * V[0] + ((5.0 / 21.0) * sqrt7 * Y5_p2 + (7.0 / 3.0) * Y3_p2 * R2) * V[1] - 5.25 * Y3_p2 * V[2];
    d[34] = Y4_p0 * Y3_p3 * V[0] + ((5.0 / 7.0) * sqrt7 * Y5_p3 - Y3_p3 * R2) * V[1] + 2.25 * Y3_p3 * V[2];
    d[30] = Y4_p0 * Y3_n1 * V[0] + ((-11.0 / 42.0) * sqrt10 * Y5_n1 + (-1.0 / 3.0) * Y3_n1 * R2 + (9.0 / 14.0) * sqrt6 * Y1_n1 * R4) * V[1] + (0.75 * Y3_n1 - 2.25 * sqrt6 * Y1_n1 * R2) * V[2] + 1.875 * sqrt6 * Y1_n1 * V[3];
    d[29] = Y4_p0 * Y3_n2 * V[0] + ((5.0 / 21.0) * sqrt7 * Y5_n2 + (7.0 / 3.0) * Y3_n2 * R2) * V[1] - 5.25 * Y3_n2 * V[2];
    d[28] = Y4_p0 * Y3_n3 * V[0] + ((5.0 / 7.0) * sqrt7 * Y5_n3 - Y3_n3 * R2) * V[1] + 2.25 * Y3_n3 * V[2];
    d[38] = Y4_p1 * Y3_p0 * V[0] + ((-19.0 / 42.0) * sqrt6 * Y5_p1 + (-1.0 / 3.0) * sqrt15 * Y3_p1 * R2 + (-9.0 / 14.0) * sqrt10 * Y1_p1 * R4) * V[1] + (0.75 * sqrt15 * Y3_p1 + 2.25 * sqrt10 * Y1_p1 * R2) * V[2] - 1.875 * sqrt10 * Y1_p1 * V[3];
    d[39] = Y4_p1 * Y3_p1 * V[0] + ((-1.0 / 42.0) * sqrt15 * Y5_p0 + (-8.0 / 21.0) * sqrt7 * Y5_p2 + (-1.0 / 3.0) * sqrt15 * Y3_p0 * R2 + (-4.0 / 3.0) * Y3_p2 * R2 + (-9.0 / 14.0) * sqrt15 * Y1_p0 * R4) * V[1] + (0.75 * sqrt15 * Y3_p0 + 3.0 * Y3_p2 + 2.25 * sqrt15 * Y1_p0 * R2) * V[2] - 1.875 * sqrt15 * Y1_p0 * V[3];
    d[40] = Y4_p1 * Y3_p2 * V[0] + ((-25.0 / 84.0) * sqrt10 * Y5_p1 + (-1.0 / 42.0) * sqrt105 * Y5_p3 + (-4.0 / 3.0) * Y3_p1 * R2 + (1.0 / 3.0) * sqrt15 * Y3_p3 * R2 + (9.0 / 28.0) * sqrt6 * Y1_p1 * R4) * V[1] + (3.0 * Y3_p1 - 0.75 * sqrt15 * Y3_p3 - 1.125 * sqrt6 * Y1_p1 * R2) * V[2] + 0.9375 * sqrt6 * Y1_p1 * V[3];
    d[41] = Y4_p1 * Y3_p3 * V[0] + ((-5.0 / 42.0) * sqrt105 * Y5_p2 + (3.0 / 14.0) * sqrt35 * Y5_p4 + (1.0 / 3.0) * sqrt15 * Y3_p2 * R2) * V[1] - 0.75 * sqrt15 * Y3_p2 * V[2];
    d[37] = Y4_p1 * Y3_n1 * V[0] + ((-8.0 / 21.0) * sqrt7 * Y5_n2 + (-4.0 / 3.0) * Y3_n2 * R2) * V[1] + 3.0 * Y3_n2 * V[2];
    d[36] = Y4_p1 * Y3_n2 * V[0] + ((-1.0 / 42.0) * sqrt105 * Y5_n3 + (-25.0 / 84.0) * sqrt10 * Y5_n1 + (1.0 / 3.0) * sqrt15 * Y3_n3 * R2 + (-4.0 / 3.0) * Y3_n1 * R2 + (9.0 / 28.0) * sqrt6 * Y1_n1 * R4) * V[1] + (-0.75 * sqrt15 * Y3_n3 + 3.0 * Y3_n1 - 1.125 * sqrt6 * Y1_n1 * R2) * V[2] + 0.9375 * sqrt6 * Y1_n1 * V[3];
    d[35] = Y4_p1 * Y3_n3 * V[0] + ((3.0 / 14.0) * sqrt35 * Y5_n4 + (-5.0 / 42.0) * sqrt105 * Y5_n2 + (1.0 / 3.0) * sqrt15 * Y3_n2 * R2) * V[1] - 0.75 * sqrt15 * Y3_n2 * V[2];
    d[45] = Y4_p2 * Y3_p0 * V[0] + ((-1.0 / 21.0) * sqrt21 * Y5_p2 + (1.0 / 3.0) * sqrt3 * Y3_p2 * R2) * V[1] - 0.75 * sqrt3 * Y3_p2 * V[2];
    d[46] = Y4_p2 * Y3_p1 * V[0] + ((-43.0 / 84.0) * sqrt2 * Y5_p1 + (-3.0 / 14.0) * sqrt21 * Y5_p3 + (-2.0 / 3.0) * sqrt5 * Y3_p1 * R2 - sqrt3 * Y3_p3 * R2 + (-9.0 / 28.0) * sqrt30 * Y1_p1 * R4) * V[1] + (1.5 * sqrt5 * Y3_p1 + 2.25 * sqrt3 * Y3_p3 + 1.125 * sqrt30 * Y1_p1 * R2) * V[2] - 0.9375 * sqrt30 * Y1_p1 * V[3];
    d[47] = Y4_p2 * Y3_p2 * V[0] + ((20.0 / 21.0) * sqrt3 * Y5_p0 + (-2.0 / 21.0) * sqrt105 * Y5_p4 + (1.0 / 3.0) * sqrt3 * Y3_p0 * R2 + (-9.0 / 7.0) * sqrt3 * Y1_p0 * R4) * V[1] + (-0.75 * sqrt3 * Y3_p0 + 4.5 * sqrt3 * Y1_p0 * R2) * V[2] - 3.75 * sqrt3 * Y1_p0 * V[3];
    d[48] = Y4_p2 * Y3_p3 * V[0] + ((5.0 / 28.0) * sqrt30 * Y5_p1 + (5.0 / 14.0) * sqrt7 * Y5_p5 - sqrt3 * Y3_p1 * R2 + (9.0 / 28.0) * sqrt2 * Y1_p1 * R4) * V[1] + (2.25 * sqrt3 * Y3_p1 - 1.125 * sqrt2 * Y1_p1 * R2) * V[2] + 0.9375 * sqrt2 * Y1_p1 * V[3];
    d[44] = Y4_p2 * Y3_n1 * V[0] + ((-3.0 / 14.0) * sqrt21 * Y5_n3 + (43.0 / 84.0) * sqrt2 * Y5_n1 - sqrt3 * Y3_n3 * R2 + (2.0 / 3.0) * sqrt5 * Y3_n1 * R2 + (9.0 / 28.0) * sqrt30 * Y1_n1 * R4) * V[1] + (2.25 * sqrt3 * Y3_n3 - 1.5 * sqrt5 * Y3_n1 - 1.125 * sqrt30 * Y1_n1 * R2) * V[2] + 0.9375 * sqrt30 * Y1_n1 * V[3];
    d[43] = Y4_p2 * Y3_n2 * V[0] + (-2.0 / 21.0) * sqrt105 * Y5_n4 * V[1];
    d[42] = Y4_p2 * Y3_n3 * V[0] + ((5.0 / 14.0) * sqrt7 * Y5_n5 + (5.0 / 28.0) * sqrt30 * Y5_n1 - sqrt3 * Y3_n1 * R2 + (9.0 / 28.0) * sqrt2 * Y1_n1 * R4) * V[1] + (2.25 * sqrt3 * Y3_n1 - 1.125 * sqrt2 * Y1_n1 * R2) * V[2] + 0.9375 * sqrt2 * Y1_n1 * V[3];
    d[52] = Y4_p3 * Y3_p0 * V[0] + (Y5_p3 + sqrt7 * Y3_p3 * R2) * V[1] - 2.25 * sqrt7 * Y3_p3 * V[2];
    d[53] = Y4_p3 * Y3_p1 * V[0] + ((-7.0 / 6.0) * Y5_p2 + (-1.0 / 6.0) * sqrt3 * Y5_p4 + (-1.0 / 3.0) * sqrt7 * Y3_p2 * R2) * V[1] + 0.75 * sqrt7 * Y3_p2 * V[2];
    d[54] = Y4_p3 * Y3_p2 * V[0] + ((11.0 / 84.0) * sqrt70 * Y5_p1 + (-5.0 / 6.0) * sqrt3 * Y5_p5 + (-1.0 / 3.0) * sqrt7 * Y3_p1 * R2 + (-9.0 / 28.0) * sqrt42 * Y1_p1 * R4) * V[1] + (0.75 * sqrt7 * Y3_p1 + 1.125 * sqrt42 * Y1_p1 * R2) * V[2] - 0.9375 * sqrt42 * Y1_p1 * V[3];
    d[55] = Y4_p3 * Y3_p3 * V[0] + ((-5.0 / 14.0) * sqrt7 * Y5_p0 + sqrt7 * Y3_p0 * R2 + (-9.0 / 14.0) * sqrt7 * Y1_p0 * R4) * V[1] + (-2.25 * sqrt7 * Y3_p0 + 2.25 * sqrt7 * Y1_p0 * R2) * V[2] - 1.875 * sqrt7 * Y1_p0 * V[3];
    d[51] = Y4_p3 * Y3_n1 * V[0] + ((-1.0 / 6.0) * sqrt3 * Y5_n4 + (7.0 / 6.0) * Y5_n2 + (1.0 / 3.0) * sqrt7 * Y3_n2 * R2) * V[1] - 0.75 * sqrt7 * Y3_n2 * V[2];
    d[50] = Y4_p3 * Y3_n2 * V[0] + ((-5.0 / 6.0) * sqrt3 * Y5_n5 + (-11.0 / 84.0) * sqrt70 * Y5_n1 + (1.0 / 3.0) * sqrt7 * Y3_n1 * R2 + (9.0 / 28.0) * sqrt42 * Y1_n1 * R4) * V[1] + (-0.75 * sqrt7 * Y3_n1 - 1.125 * sqrt42 * Y1_n1 * R2) * V[2] + 0.9375 * sqrt42 * Y1_n1 * V[3];
    d[49] = Y4_p3 * Y3_n3 * V[0];
    d[59] = Y4_p4 * Y3_p0 * V[0] + 2.0 * Y5_p4 * V[1];
    d[60] = Y4_p4 * Y3_p1 * V[0] + ((-2.0 / 3.0) * sqrt3 * Y5_p3 + (1.0 / 3.0) * sqrt15 * Y5_p5 + (1.0 / 3.0) * sqrt21 * Y3_p3 * R2) * V[1] - 0.75 * sqrt21 * Y3_p3 * V[2];
    d[61] = Y4_p4 * Y3_p2 * V[0] + ((1.0 / 3.0) * sqrt5 * Y5_p2 + (-1.0 / 3.0) * sqrt35 * Y3_p2 * R2) * V[1] + 0.75 * sqrt35 * Y3_p2 * V[2];
    d[62] = Y4_p4 * Y3_p3 * V[0] + ((-1.0 / 42.0) * sqrt210 * Y5_p1 + (1.0 / 3.0) * sqrt21 * Y3_p1 * R2 + (-9.0 / 14.0) * sqrt14 * Y1_p1 * R4) * V[1] + (-0.75 * sqrt21 * Y3_p1 + 2.25 * sqrt14 * Y1_p1 * R2) * V[2] - 1.875 * sqrt14 * Y1_p1 * V[3];
    d[58] = Y4_p4 * Y3_n1 * V[0] + ((1.0 / 3.0) * sqrt15 * Y5_n5 + (2.0 / 3.0) * sqrt3 * Y5_n3 + (-1.0 / 3.0) * sqrt21 * Y3_n3 * R2) * V[1] + 0.75 * sqrt21 * Y3_n3 * V[2];
    d[57] = Y4_p4 * Y3_n2 * V[0] + ((-1.0 / 3.0) * sqrt5 * Y5_n2 + (1.0 / 3.0) * sqrt35 * Y3_n2 * R2) * V[1] - 0.75 * sqrt35 * Y3_n2 * V[2];
    d[56] = Y4_p4 * Y3_n3 * V[0] + ((1.0 / 42.0) * sqrt210 * Y5_n1 + (-1.0 / 3.0) * sqrt21 * Y3_n1 * R2 + (9.0 / 14.0) * sqrt14 * Y1_n1 * R4) * V[1] + (0.75 * sqrt21 * Y3_n1 - 2.25 * sqrt14 * Y1_n1 * R2) * V[2] + 1.875 * sqrt14 * Y1_n1 * V[3];
    d[24] = Y4_n1 * Y3_p0 * V[0] + ((-19.0 / 42.0) * sqrt6 * Y5_n1 + (-1.0 / 3.0) * sqrt15 * Y3_n1 * R2 + (-9.0 / 14.0) * sqrt10 * Y1_n1 * R4) * V[1] + (0.75 * sqrt15 * Y3_n1 + 2.25 * sqrt10 * Y1_n1 * R2) * V[2] - 1.875 * sqrt10 * Y1_n1 * V[3];
    d[25] = Y4_n1 * Y3_p1 * V[0] + ((-8.0 / 21.0) * sqrt7 * Y5_n2 + (-4.0 / 3.0) * Y3_n2 * R2) * V[1] + 3.0 * Y3_n2 * V[2];
    d[26] = Y4_n1 * Y3_p2 * V[0] + ((-1.0 / 42.0) * sqrt105 * Y5_n3 + (25.0 / 84.0) * sqrt10 * Y5_n1 + (1.0 / 3.0) * sqrt15 * Y3_n3 * R2 + (4.0 / 3.0) * Y3_n1 * R2 + (-9.0 / 28.0) * sqrt6 * Y1_n1 * R4) * V[1] + (-0.75 * sqrt15 * Y3_n3 - 3.0 * Y3_n1 + 1.125 * sqrt6 * Y1_n1 * R2) * V[2] - 0.9375 * sqrt6 * Y1_n1 * V[3];
    d[27] = Y4_n1 * Y3_p3 * V[0] + ((3.0 / 14.0) * sqrt35 * Y5_n4 + (5.0 / 42.0) * sqrt105 * Y5_n2 + (-1.0 / 3.0) * sqrt15 * Y3_n2 * R2) * V[1] + 0.75 * sqrt15 * Y3_n2 * V[2];
    d[23] = Y4_n1 * Y3_n1 * V[0] + ((-1.0 / 42.0) * sqrt15 * Y5_p0 + (8.0 / 21.0) * sqrt7 * Y5_p2 + (-1.0 / 3.0) * sqrt15 * Y3_p0 * R2 + (4.0 / 3.0) * Y3_p2 * R2 + (-9.0 / 14.0) * sqrt15 * Y1_p0 * R4) * V[1] + (0.75 * sqrt15 * Y3_p0 - 3.0 * Y3_p2 + 2.25 * sqrt15 * Y1_p0 * R2) * V[2] - 1.875 * sqrt15 * Y1_p0 * V[3];
    d[22] = Y4_n1 * Y3_n2 * V[0] + ((-25.0 / 84.0) * sqrt10 * Y5_p1 + (1.0 / 42.0) * sqrt105 * Y5_p3 + (-4.0 / 3.0) * Y3_p1 * R2 + (-1.0 / 3.0) * sqrt15 * Y3_p3 * R2 + (9.0 / 28.0) * sqrt6 * Y1_p1 * R4) * V[1] + (3.0 * Y3_p1 + 0.75 * sqrt15 * Y3_p3 - 1.125 * sqrt6 * Y1_p1 * R2) * V[2] + 0.9375 * sqrt6 * Y1_p1 * V[3];
    d[21] = Y4_n1 * Y3_n3 * V[0] + ((-5.0 / 42.0) * sqrt105 * Y5_p2 + (-3.0 / 14.0) * sqrt35 * Y5_p4 + (1.0 / 3.0) * sqrt15 * Y3_p2 * R2) * V[1] - 0.75 * sqrt15 * Y3_p2 * V[2];
    d[17] = Y4_n2 * Y3_p0 * V[0] + ((-1.0 / 21.0) * sqrt21 * Y5_n2 + (1.0 / 3.0) * sqrt3 * Y3_n2 * R2) * V[1] - 0.75 * sqrt3 * Y3_n2 * V[2];
    d[18] = Y4_n2 * Y3_p1 * V[0] + ((-3.0 / 14.0) * sqrt21 * Y5_n3 + (-43.0 / 84.0) * sqrt2 * Y5_n1 - sqrt3 * Y3_n3 * R2 + (-2.0 / 3.0) * sqrt5 * Y3_n1 * R2 + (-9.0 / 28.0) * sqrt30 * Y1_n1 * R4) * V[1] + (2.25 * sqrt3 * Y3_n3 + 1.5 * sqrt5 * Y3_n1 + 1.125 * sqrt30 * Y1_n1 * R2) * V[2] - 0.9375 * sqrt30 * Y1_n1 * V[3];
    d[19] = Y4_n2 * Y3_p2 * V[0] + (-2.0 / 21.0) * sqrt105 * Y5_n4 * V[1];
    d[20] = Y4_n2 * Y3_p3 * V[0] + ((5.0 / 14.0) * sqrt7 * Y5_n5 + (-5.0 / 28.0) * sqrt30 * Y5_n1 + sqrt3 * Y3_n1 * R2 + (-9.0 / 28.0) * sqrt2 * Y1_n1 * R4) * V[1] + (-2.25 * sqrt3 * Y3_n1 + 1.125 * sqrt2 * Y1_n1 * R2) * V[2] - 0.9375 * sqrt2 * Y1_n1 * V[3];
    d[16] = Y4_n2 * Y3_n1 * V[0] + ((-43.0 / 84.0) * sqrt2 * Y5_p1 + (3.0 / 14.0) * sqrt21 * Y5_p3 + (-2.0 / 3.0) * sqrt5 * Y3_p1 * R2 + sqrt3 * Y3_p3 * R2 + (-9.0 / 28.0) * sqrt30 * Y1_p1 * R4) * V[1] + (1.5 * sqrt5 * Y3_p1 - 2.25 * sqrt3 * Y3_p3 + 1.125 * sqrt30 * Y1_p1 * R2) * V[2] - 0.9375 * sqrt30 * Y1_p1 * V[3];
    d[15] = Y4_n2 * Y3_n2 * V[0] + ((20.0 / 21.0) * sqrt3 * Y5_p0 + (2.0 / 21.0) * sqrt105 * Y5_p4 + (1.0 / 3.0) * sqrt3 * Y3_p0 * R2 + (-9.0 / 7.0) * sqrt3 * Y1_p0 * R4) * V[1] + (-0.75 * sqrt3 * Y3_p0 + 4.5 * sqrt3 * Y1_p0 * R2) * V[2] - 3.75 * sqrt3 * Y1_p0 * V[3];
    d[14] = Y4_n2 * Y3_n3 * V[0] + ((5.0 / 28.0) * sqrt30 * Y5_p1 + (-5.0 / 14.0) * sqrt7 * Y5_p5 - sqrt3 * Y3_p1 * R2 + (9.0 / 28.0) * sqrt2 * Y1_p1 * R4) * V[1] + (2.25 * sqrt3 * Y3_p1 - 1.125 * sqrt2 * Y1_p1 * R2) * V[2] + 0.9375 * sqrt2 * Y1_p1 * V[3];
    d[10] = Y4_n3 * Y3_p0 * V[0] + (Y5_n3 + sqrt7 * Y3_n3 * R2) * V[1] - 2.25 * sqrt7 * Y3_n3 * V[2];
    d[11] = Y4_n3 * Y3_p1 * V[0] + ((-1.0 / 6.0) * sqrt3 * Y5_n4 + (-7.0 / 6.0) * Y5_n2 + (-1.0 / 3.0) * sqrt7 * Y3_n2 * R2) * V[1] + 0.75 * sqrt7 * Y3_n2 * V[2];
    d[12] = Y4_n3 * Y3_p2 * V[0] + ((-5.0 / 6.0) * sqrt3 * Y5_n5 + (11.0 / 84.0) * sqrt70 * Y5_n1 + (-1.0 / 3.0) * sqrt7 * Y3_n1 * R2 + (-9.0 / 28.0) * sqrt42 * Y1_n1 * R4) * V[1] + (0.75 * sqrt7 * Y3_n1 + 1.125 * sqrt42 * Y1_n1 * R2) * V[2] - 0.9375 * sqrt42 * Y1_n1 * V[3];
    d[13] = Y4_n3 * Y3_p3 * V[0];
    d[9] = Y4_n3 * Y3_n1 * V[0] + ((-7.0 / 6.0) * Y5_p2 + (1.0 / 6.0) * sqrt3 * Y5_p4 + (-1.0 / 3.0) * sqrt7 * Y3_p2 * R2) * V[1] + 0.75 * sqrt7 * Y3_p2 * V[2];
    d[8] = Y4_n3 * Y3_n2 * V[0] + ((11.0 / 84.0) * sqrt70 * Y5_p1 + (5.0 / 6.0) * sqrt3 * Y5_p5 + (-1.0 / 3.0) * sqrt7 * Y3_p1 * R2 + (-9.0 / 28.0) * sqrt42 * Y1_p1 * R4) * V[1] + (0.75 * sqrt7 * Y3_p1 + 1.125 * sqrt42 * Y1_p1 * R2) * V[2] - 0.9375 * sqrt42 * Y1_p1 * V[3];
    d[7] = Y4_n3 * Y3_n3 * V[0] + ((-5.0 / 14.0) * sqrt7 * Y5_p0 + sqrt7 * Y3_p0 * R2 + (-9.0 / 14.0) * sqrt7 * Y1_p0 * R4) * V[1] + (-2.25 * sqrt7 * Y3_p0 + 2.25 * sqrt7 * Y1_p0 * R2) * V[2] - 1.875 * sqrt7 * Y1_p0 * V[3];
    d[3] = Y4_n4 * Y3_p0 * V[0] + 2.0 * Y5_n4 * V[1];
    d[4] = Y4_n4 * Y3_p1 * V[0] + ((1.0 / 3.0) * sqrt15 * Y5_n5 + (-2.0 / 3.0) * sqrt3 * Y5_n3 + (1.0 / 3.0) * sqrt21 * Y3_n3 * R2) * V[1] - 0.75 * sqrt21 * Y3_n3 * V[2];
    d[5] = Y4_n4 * Y3_p2 * V[0] + ((1.0 / 3.0) * sqrt5 * Y5_n2 + (-1.0 / 3.0) * sqrt35 * Y3_n2 * R2) * V[1] + 0.75 * sqrt35 * Y3_n2 * V[2];
    d[6] = Y4_n4 * Y3_p3 * V[0] + ((-1.0 / 42.0) * sqrt210 * Y5_n1 + (1.0 / 3.0) * sqrt21 * Y3_n1 * R2 + (-9.0 / 14.0) * sqrt14 * Y1_n1 * R4) * V[1] + (-0.75 * sqrt21 * Y3_n1 + 2.25 * sqrt14 * Y1_n1 * R2) * V[2] - 1.875 * sqrt14 * Y1_n1 * V[3];
    d[2] = Y4_n4 * Y3_n1 * V[0] + ((-2.0 / 3.0) * sqrt3 * Y5_p3 + (-1.0 / 3.0) * sqrt15 * Y5_p5 + (1.0 / 3.0) * sqrt21 * Y3_p3 * R2) * V[1] - 0.75 * sqrt21 * Y3_p3 * V[2];
    d[1] = Y4_n4 * Y3_n2 * V[0] + ((1.0 / 3.0) * sqrt5 * Y5_p2 + (-1.0 / 3.0) * sqrt35 * Y3_p2 * R2) * V[1] + 0.75 * sqrt35 * Y3_p2 * V[2];
    d[0] = Y4_n4 * Y3_n3 * V[0] + ((-1.0 / 42.0) * sqrt210 * Y5_p1 + (1.0 / 3.0) * sqrt21 * Y3_p1 * R2 + (-9.0 / 14.0) * sqrt14 * Y1_p1 * R4) * V[1] + (-0.75 * sqrt21 * Y3_p1 + 2.25 * sqrt14 * Y1_p1 * R2) * V[2] - 1.875 * sqrt14 * Y1_p1 * V[3];

    return out;
}

}  // namespace ovlab