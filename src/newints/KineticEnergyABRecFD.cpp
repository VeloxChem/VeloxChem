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

#include "KineticEnergyABRecFD.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_f_d(
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
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt210 = std::sqrt(210.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^2 · β^3 · p^{-5} · (s|T|s)
    //   V[1] ↔ α · β^2 · p^{-4} · (s|T|s)
    //   V[2] ↔ β · p^{-3} · (s|T|s)
    //   V[3] ↔ α^2 · β^3 · p^{-5} · ξ · (s|S|s)
    //   V[4] ↔ α · β^2 · p^{-4} · ξ · (s|S|s)
    //   V[5] ↔ β · p^{-3} · ξ · (s|S|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha2 * beta3 * pinv5;
            V[1] += cab_tt * alpha * beta2 * pinv4;
            V[2] += cab_tt * beta * pinv3;
            V[3] += cab_ss * alpha2 * beta3 * pinv5 * xi;
            V[4] += cab_ss * alpha * beta2 * pinv4 * xi;
            V[5] += cab_ss * beta * pinv3 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 7 × 5 spherical block ----
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
    auto *d = buffer;
    d[17] = ((-10.0 / 21.0) * Y5_p0 + (-4.0 / 15.0) * Y3_p0 * R2 + (-9.0 / 35.0) * Y1_p0 * R4) * V[0] + (1.2 * Y3_p0 + 1.8 * Y1_p0 * R2) * V[1] - 2.25 * Y1_p0 * V[2] + ((-100.0 / 21.0) * Y5_p0 + (-8.0 / 3.0) * Y3_p0 * R2 + (-18.0 / 7.0) * Y1_p0 * R4) * V[3] + (9.6 * Y3_p0 + 14.4 * Y1_p0 * R2) * V[4] - 13.5 * Y1_p0 * V[5];
    d[18] = ((-4.0 / 21.0) * sqrt5 * Y5_p1 + (-1.0 / 15.0) * sqrt2 * Y3_p1 * R2 + (3.0 / 35.0) * sqrt3 * Y1_p1 * R4) * V[0] + (0.3 * sqrt2 * Y3_p1 - 0.6 * sqrt3 * Y1_p1 * R2) * V[1] + 0.75 * sqrt3 * Y1_p1 * V[2] + ((-40.0 / 21.0) * sqrt5 * Y5_p1 + (-2.0 / 3.0) * sqrt2 * Y3_p1 * R2 + (6.0 / 7.0) * sqrt3 * Y1_p1 * R4) * V[3] + (2.4 * sqrt2 * Y3_p1 - 4.8 * sqrt3 * Y1_p1 * R2) * V[4] + 4.5 * sqrt3 * Y1_p1 * V[5];
    d[19] = ((-1.0 / 21.0) * sqrt35 * Y5_p2 + (2.0 / 15.0) * sqrt5 * Y3_p2 * R2) * V[0] - 0.6 * sqrt5 * Y3_p2 * V[1] + ((-10.0 / 21.0) * sqrt35 * Y5_p2 + (4.0 / 3.0) * sqrt5 * Y3_p2 * R2) * V[3] - 4.8 * sqrt5 * Y3_p2 * V[4];
    d[16] = ((-4.0 / 21.0) * sqrt5 * Y5_n1 + (-1.0 / 15.0) * sqrt2 * Y3_n1 * R2 + (3.0 / 35.0) * sqrt3 * Y1_n1 * R4) * V[0] + (0.3 * sqrt2 * Y3_n1 - 0.6 * sqrt3 * Y1_n1 * R2) * V[1] + 0.75 * sqrt3 * Y1_n1 * V[2] + ((-40.0 / 21.0) * sqrt5 * Y5_n1 + (-2.0 / 3.0) * sqrt2 * Y3_n1 * R2 + (6.0 / 7.0) * sqrt3 * Y1_n1 * R4) * V[3] + (2.4 * sqrt2 * Y3_n1 - 4.8 * sqrt3 * Y1_n1 * R2) * V[4] + 4.5 * sqrt3 * Y1_n1 * V[5];
    d[15] = ((-1.0 / 21.0) * sqrt35 * Y5_n2 + (2.0 / 15.0) * sqrt5 * Y3_n2 * R2) * V[0] - 0.6 * sqrt5 * Y3_n2 * V[1] + ((-10.0 / 21.0) * sqrt35 * Y5_n2 + (4.0 / 3.0) * sqrt5 * Y3_n2 * R2) * V[3] - 4.8 * sqrt5 * Y3_n2 * V[4];
    d[22] = ((-1.0 / 7.0) * sqrt10 * Y5_p1 - 0.2 * Y3_p1 * R2 + (-3.0 / 35.0) * sqrt6 * Y1_p1 * R4) * V[0] + (0.9 * Y3_p1 + 0.6 * sqrt6 * Y1_p1 * R2) * V[1] - 0.75 * sqrt6 * Y1_p1 * V[2] + ((-10.0 / 7.0) * sqrt10 * Y5_p1 - 2.0 * Y3_p1 * R2 + (-6.0 / 7.0) * sqrt6 * Y1_p1 * R4) * V[3] + (7.2 * Y3_p1 + 4.8 * sqrt6 * Y1_p1 * R2) * V[4] - 4.5 * sqrt6 * Y1_p1 * V[5];
    d[23] = ((5.0 / 21.0) * sqrt2 * Y5_p0 + (-1.0 / 42.0) * sqrt210 * Y5_p2 + (-1.0 / 15.0) * sqrt2 * Y3_p0 * R2 + (-1.0 / 30.0) * sqrt30 * Y3_p2 * R2 + (-6.0 / 35.0) * sqrt2 * Y1_p0 * R4) * V[0] + (0.3 * sqrt2 * Y3_p0 + 0.15 * sqrt30 * Y3_p2 + 1.2 * sqrt2 * Y1_p0 * R2) * V[1] - 1.5 * sqrt2 * Y1_p0 * V[2] + ((50.0 / 21.0) * sqrt2 * Y5_p0 + (-5.0 / 21.0) * sqrt210 * Y5_p2 + (-2.0 / 3.0) * sqrt2 * Y3_p0 * R2 + (-1.0 / 3.0) * sqrt30 * Y3_p2 * R2 + (-12.0 / 7.0) * sqrt2 * Y1_p0 * R4) * V[3] + (2.4 * sqrt2 * Y3_p0 + 1.2 * sqrt30 * Y3_p2 + 9.6 * sqrt2 * Y1_p0 * R2) * V[4] - 9.0 * sqrt2 * Y1_p0 * V[5];
    d[24] = ((1.0 / 42.0) * sqrt30 * Y5_p1 + (-1.0 / 21.0) * sqrt35 * Y5_p3 + (-2.0 / 15.0) * sqrt3 * Y3_p1 * R2 + (1.0 / 15.0) * sqrt5 * Y3_p3 * R2 + (3.0 / 70.0) * sqrt2 * Y1_p1 * R4) * V[0] + (0.6 * sqrt3 * Y3_p1 - 0.3 * sqrt5 * Y3_p3 - 0.3 * sqrt2 * Y1_p1 * R2) * V[1] + 0.375 * sqrt2 * Y1_p1 * V[2] + ((5.0 / 21.0) * sqrt30 * Y5_p1 + (-10.0 / 21.0) * sqrt35 * Y5_p3 + (-4.0 / 3.0) * sqrt3 * Y3_p1 * R2 + (2.0 / 3.0) * sqrt5 * Y3_p3 * R2 + (3.0 / 7.0) * sqrt2 * Y1_p1 * R4) * V[3] + (4.8 * sqrt3 * Y3_p1 - 2.4 * sqrt5 * Y3_p3 - 2.4 * sqrt2 * Y1_p1 * R2) * V[4] + 2.25 * sqrt2 * Y1_p1 * V[5];
    d[21] = ((-1.0 / 42.0) * sqrt210 * Y5_n2 + (-1.0 / 30.0) * sqrt30 * Y3_n2 * R2) * V[0] + 0.15 * sqrt30 * Y3_n2 * V[1] + ((-5.0 / 21.0) * sqrt210 * Y5_n2 + (-1.0 / 3.0) * sqrt30 * Y3_n2 * R2) * V[3] + 1.2 * sqrt30 * Y3_n2 * V[4];
    d[20] = ((-1.0 / 21.0) * sqrt35 * Y5_n3 + (1.0 / 42.0) * sqrt30 * Y5_n1 + (1.0 / 15.0) * sqrt5 * Y3_n3 * R2 + (-2.0 / 15.0) * sqrt3 * Y3_n1 * R2 + (3.0 / 70.0) * sqrt2 * Y1_n1 * R4) * V[0] + (-0.3 * sqrt5 * Y3_n3 + 0.6 * sqrt3 * Y3_n1 - 0.3 * sqrt2 * Y1_n1 * R2) * V[1] + 0.375 * sqrt2 * Y1_n1 * V[2] + ((-10.0 / 21.0) * sqrt35 * Y5_n3 + (5.0 / 21.0) * sqrt30 * Y5_n1 + (2.0 / 3.0) * sqrt5 * Y3_n3 * R2 + (-4.0 / 3.0) * sqrt3 * Y3_n1 * R2 + (3.0 / 7.0) * sqrt2 * Y1_n1 * R4) * V[3] + (-2.4 * sqrt5 * Y3_n3 + 4.8 * sqrt3 * Y3_n1 - 2.4 * sqrt2 * Y1_n1 * R2) * V[4] + 2.25 * sqrt2 * Y1_n1 * V[5];
    d[27] = (-1.0 / 7.0) * sqrt7 * Y5_p2 * V[0] + (-10.0 / 7.0) * sqrt7 * Y5_p2 * V[3];
    d[28] = ((2.0 / 21.0) * sqrt3 * Y5_p1 + (-2.0 / 21.0) * sqrt14 * Y5_p3 + (-1.0 / 30.0) * sqrt30 * Y3_p1 * R2 + (-1.0 / 6.0) * sqrt2 * Y3_p3 * R2 + (-3.0 / 35.0) * sqrt5 * Y1_p1 * R4) * V[0] + (0.15 * sqrt30 * Y3_p1 + 0.75 * sqrt2 * Y3_p3 + 0.6 * sqrt5 * Y1_p1 * R2) * V[1] - 0.75 * sqrt5 * Y1_p1 * V[2] + ((20.0 / 21.0) * sqrt3 * Y5_p1 + (-20.0 / 21.0) * sqrt14 * Y5_p3 + (-1.0 / 3.0) * sqrt30 * Y3_p1 * R2 + (-5.0 / 3.0) * sqrt2 * Y3_p3 * R2 + (-6.0 / 7.0) * sqrt5 * Y1_p1 * R4) * V[3] + (1.2 * sqrt30 * Y3_p1 + 6.0 * sqrt2 * Y3_p3 + 4.8 * sqrt5 * Y1_p1 * R2) * V[4] - 4.5 * sqrt5 * Y1_p1 * V[5];
    d[29] = ((-1.0 / 21.0) * sqrt5 * Y5_p0 + (-1.0 / 7.0) * sqrt7 * Y5_p4 + (2.0 / 15.0) * sqrt5 * Y3_p0 * R2 + (-3.0 / 35.0) * sqrt5 * Y1_p0 * R4) * V[0] + (-0.6 * sqrt5 * Y3_p0 + 0.6 * sqrt5 * Y1_p0 * R2) * V[1] - 0.75 * sqrt5 * Y1_p0 * V[2] + ((-10.0 / 21.0) * sqrt5 * Y5_p0 + (-10.0 / 7.0) * sqrt7 * Y5_p4 + (4.0 / 3.0) * sqrt5 * Y3_p0 * R2 + (-6.0 / 7.0) * sqrt5 * Y1_p0 * R4) * V[3] + (-4.8 * sqrt5 * Y3_p0 + 4.8 * sqrt5 * Y1_p0 * R2) * V[4] - 4.5 * sqrt5 * Y1_p0 * V[5];
    d[26] = ((-2.0 / 21.0) * sqrt14 * Y5_n3 + (-2.0 / 21.0) * sqrt3 * Y5_n1 + (-1.0 / 6.0) * sqrt2 * Y3_n3 * R2 + (1.0 / 30.0) * sqrt30 * Y3_n1 * R2 + (3.0 / 35.0) * sqrt5 * Y1_n1 * R4) * V[0] + (0.75 * sqrt2 * Y3_n3 - 0.15 * sqrt30 * Y3_n1 - 0.6 * sqrt5 * Y1_n1 * R2) * V[1] + 0.75 * sqrt5 * Y1_n1 * V[2] + ((-20.0 / 21.0) * sqrt14 * Y5_n3 + (-20.0 / 21.0) * sqrt3 * Y5_n1 + (-5.0 / 3.0) * sqrt2 * Y3_n3 * R2 + (1.0 / 3.0) * sqrt30 * Y3_n1 * R2 + (6.0 / 7.0) * sqrt5 * Y1_n1 * R4) * V[3] + (6.0 * sqrt2 * Y3_n3 - 1.2 * sqrt30 * Y3_n1 - 4.8 * sqrt5 * Y1_n1 * R2) * V[4] + 4.5 * sqrt5 * Y1_n1 * V[5];
    d[25] = (-1.0 / 7.0) * sqrt7 * Y5_n4 * V[0] + (-10.0 / 7.0) * sqrt7 * Y5_n4 * V[3];
    d[32] = ((-2.0 / 21.0) * sqrt7 * Y5_p3 + (1.0 / 3.0) * Y3_p3 * R2) * V[0] - 1.5 * Y3_p3 * V[1] + ((-20.0 / 21.0) * sqrt7 * Y5_p3 + (10.0 / 3.0) * Y3_p3 * R2) * V[3] - 12.0 * Y3_p3 * V[4];
    d[33] = ((1.0 / 42.0) * sqrt14 * Y5_p2 + (-1.0 / 21.0) * sqrt42 * Y5_p4 + (-1.0 / 6.0) * sqrt2 * Y3_p2 * R2) * V[0] + 0.75 * sqrt2 * Y3_p2 * V[1] + ((5.0 / 21.0) * sqrt14 * Y5_p2 + (-10.0 / 21.0) * sqrt42 * Y5_p4 + (-5.0 / 3.0) * sqrt2 * Y3_p2 * R2) * V[3] + 6.0 * sqrt2 * Y3_p2 * V[4];
    d[34] = ((-1.0 / 42.0) * sqrt2 * Y5_p1 + (-1.0 / 21.0) * sqrt105 * Y5_p5 + (1.0 / 15.0) * sqrt5 * Y3_p1 * R2 + (-3.0 / 70.0) * sqrt30 * Y1_p1 * R4) * V[0] + (-0.3 * sqrt5 * Y3_p1 + 0.3 * sqrt30 * Y1_p1 * R2) * V[1] - 0.375 * sqrt30 * Y1_p1 * V[2] + ((-5.0 / 21.0) * sqrt2 * Y5_p1 + (-10.0 / 21.0) * sqrt105 * Y5_p5 + (2.0 / 3.0) * sqrt5 * Y3_p1 * R2 + (-3.0 / 7.0) * sqrt30 * Y1_p1 * R4) * V[3] + (-2.4 * sqrt5 * Y3_p1 + 2.4 * sqrt30 * Y1_p1 * R2) * V[4] - 2.25 * sqrt30 * Y1_p1 * V[5];
    d[31] = ((-1.0 / 21.0) * sqrt42 * Y5_n4 + (-1.0 / 42.0) * sqrt14 * Y5_n2 + (1.0 / 6.0) * sqrt2 * Y3_n2 * R2) * V[0] - 0.75 * sqrt2 * Y3_n2 * V[1] + ((-10.0 / 21.0) * sqrt42 * Y5_n4 + (-5.0 / 21.0) * sqrt14 * Y5_n2 + (5.0 / 3.0) * sqrt2 * Y3_n2 * R2) * V[3] - 6.0 * sqrt2 * Y3_n2 * V[4];
    d[30] = ((-1.0 / 21.0) * sqrt105 * Y5_n5 + (1.0 / 42.0) * sqrt2 * Y5_n1 + (-1.0 / 15.0) * sqrt5 * Y3_n1 * R2 + (3.0 / 70.0) * sqrt30 * Y1_n1 * R4) * V[0] + (0.3 * sqrt5 * Y3_n1 - 0.3 * sqrt30 * Y1_n1 * R2) * V[1] + 0.375 * sqrt30 * Y1_n1 * V[2] + ((-10.0 / 21.0) * sqrt105 * Y5_n5 + (5.0 / 21.0) * sqrt2 * Y5_n1 + (-2.0 / 3.0) * sqrt5 * Y3_n1 * R2 + (3.0 / 7.0) * sqrt30 * Y1_n1 * R4) * V[3] + (2.4 * sqrt5 * Y3_n1 - 2.4 * sqrt30 * Y1_n1 * R2) * V[4] + 2.25 * sqrt30 * Y1_n1 * V[5];
    d[12] = ((-1.0 / 7.0) * sqrt10 * Y5_n1 - 0.2 * Y3_n1 * R2 + (-3.0 / 35.0) * sqrt6 * Y1_n1 * R4) * V[0] + (0.9 * Y3_n1 + 0.6 * sqrt6 * Y1_n1 * R2) * V[1] - 0.75 * sqrt6 * Y1_n1 * V[2] + ((-10.0 / 7.0) * sqrt10 * Y5_n1 - 2.0 * Y3_n1 * R2 + (-6.0 / 7.0) * sqrt6 * Y1_n1 * R4) * V[3] + (7.2 * Y3_n1 + 4.8 * sqrt6 * Y1_n1 * R2) * V[4] - 4.5 * sqrt6 * Y1_n1 * V[5];
    d[13] = ((-1.0 / 42.0) * sqrt210 * Y5_n2 + (-1.0 / 30.0) * sqrt30 * Y3_n2 * R2) * V[0] + 0.15 * sqrt30 * Y3_n2 * V[1] + ((-5.0 / 21.0) * sqrt210 * Y5_n2 + (-1.0 / 3.0) * sqrt30 * Y3_n2 * R2) * V[3] + 1.2 * sqrt30 * Y3_n2 * V[4];
    d[14] = ((-1.0 / 21.0) * sqrt35 * Y5_n3 + (-1.0 / 42.0) * sqrt30 * Y5_n1 + (1.0 / 15.0) * sqrt5 * Y3_n3 * R2 + (2.0 / 15.0) * sqrt3 * Y3_n1 * R2 + (-3.0 / 70.0) * sqrt2 * Y1_n1 * R4) * V[0] + (-0.3 * sqrt5 * Y3_n3 - 0.6 * sqrt3 * Y3_n1 + 0.3 * sqrt2 * Y1_n1 * R2) * V[1] - 0.375 * sqrt2 * Y1_n1 * V[2] + ((-10.0 / 21.0) * sqrt35 * Y5_n3 + (-5.0 / 21.0) * sqrt30 * Y5_n1 + (2.0 / 3.0) * sqrt5 * Y3_n3 * R2 + (4.0 / 3.0) * sqrt3 * Y3_n1 * R2 + (-3.0 / 7.0) * sqrt2 * Y1_n1 * R4) * V[3] + (-2.4 * sqrt5 * Y3_n3 - 4.8 * sqrt3 * Y3_n1 + 2.4 * sqrt2 * Y1_n1 * R2) * V[4] - 2.25 * sqrt2 * Y1_n1 * V[5];
    d[11] = ((5.0 / 21.0) * sqrt2 * Y5_p0 + (1.0 / 42.0) * sqrt210 * Y5_p2 + (-1.0 / 15.0) * sqrt2 * Y3_p0 * R2 + (1.0 / 30.0) * sqrt30 * Y3_p2 * R2 + (-6.0 / 35.0) * sqrt2 * Y1_p0 * R4) * V[0] + (0.3 * sqrt2 * Y3_p0 - 0.15 * sqrt30 * Y3_p2 + 1.2 * sqrt2 * Y1_p0 * R2) * V[1] - 1.5 * sqrt2 * Y1_p0 * V[2] + ((50.0 / 21.0) * sqrt2 * Y5_p0 + (5.0 / 21.0) * sqrt210 * Y5_p2 + (-2.0 / 3.0) * sqrt2 * Y3_p0 * R2 + (1.0 / 3.0) * sqrt30 * Y3_p2 * R2 + (-12.0 / 7.0) * sqrt2 * Y1_p0 * R4) * V[3] + (2.4 * sqrt2 * Y3_p0 - 1.2 * sqrt30 * Y3_p2 + 9.6 * sqrt2 * Y1_p0 * R2) * V[4] - 9.0 * sqrt2 * Y1_p0 * V[5];
    d[10] = ((1.0 / 42.0) * sqrt30 * Y5_p1 + (1.0 / 21.0) * sqrt35 * Y5_p3 + (-2.0 / 15.0) * sqrt3 * Y3_p1 * R2 + (-1.0 / 15.0) * sqrt5 * Y3_p3 * R2 + (3.0 / 70.0) * sqrt2 * Y1_p1 * R4) * V[0] + (0.6 * sqrt3 * Y3_p1 + 0.3 * sqrt5 * Y3_p3 - 0.3 * sqrt2 * Y1_p1 * R2) * V[1] + 0.375 * sqrt2 * Y1_p1 * V[2] + ((5.0 / 21.0) * sqrt30 * Y5_p1 + (10.0 / 21.0) * sqrt35 * Y5_p3 + (-4.0 / 3.0) * sqrt3 * Y3_p1 * R2 + (-2.0 / 3.0) * sqrt5 * Y3_p3 * R2 + (3.0 / 7.0) * sqrt2 * Y1_p1 * R4) * V[3] + (4.8 * sqrt3 * Y3_p1 + 2.4 * sqrt5 * Y3_p3 - 2.4 * sqrt2 * Y1_p1 * R2) * V[4] + 2.25 * sqrt2 * Y1_p1 * V[5];
    d[7] = (-1.0 / 7.0) * sqrt7 * Y5_n2 * V[0] + (-10.0 / 7.0) * sqrt7 * Y5_n2 * V[3];
    d[8] = ((-2.0 / 21.0) * sqrt14 * Y5_n3 + (2.0 / 21.0) * sqrt3 * Y5_n1 + (-1.0 / 6.0) * sqrt2 * Y3_n3 * R2 + (-1.0 / 30.0) * sqrt30 * Y3_n1 * R2 + (-3.0 / 35.0) * sqrt5 * Y1_n1 * R4) * V[0] + (0.75 * sqrt2 * Y3_n3 + 0.15 * sqrt30 * Y3_n1 + 0.6 * sqrt5 * Y1_n1 * R2) * V[1] - 0.75 * sqrt5 * Y1_n1 * V[2] + ((-20.0 / 21.0) * sqrt14 * Y5_n3 + (20.0 / 21.0) * sqrt3 * Y5_n1 + (-5.0 / 3.0) * sqrt2 * Y3_n3 * R2 + (-1.0 / 3.0) * sqrt30 * Y3_n1 * R2 + (-6.0 / 7.0) * sqrt5 * Y1_n1 * R4) * V[3] + (6.0 * sqrt2 * Y3_n3 + 1.2 * sqrt30 * Y3_n1 + 4.8 * sqrt5 * Y1_n1 * R2) * V[4] - 4.5 * sqrt5 * Y1_n1 * V[5];
    d[9] = (-1.0 / 7.0) * sqrt7 * Y5_n4 * V[0] + (-10.0 / 7.0) * sqrt7 * Y5_n4 * V[3];
    d[6] = ((2.0 / 21.0) * sqrt3 * Y5_p1 + (2.0 / 21.0) * sqrt14 * Y5_p3 + (-1.0 / 30.0) * sqrt30 * Y3_p1 * R2 + (1.0 / 6.0) * sqrt2 * Y3_p3 * R2 + (-3.0 / 35.0) * sqrt5 * Y1_p1 * R4) * V[0] + (0.15 * sqrt30 * Y3_p1 - 0.75 * sqrt2 * Y3_p3 + 0.6 * sqrt5 * Y1_p1 * R2) * V[1] - 0.75 * sqrt5 * Y1_p1 * V[2] + ((20.0 / 21.0) * sqrt3 * Y5_p1 + (20.0 / 21.0) * sqrt14 * Y5_p3 + (-1.0 / 3.0) * sqrt30 * Y3_p1 * R2 + (5.0 / 3.0) * sqrt2 * Y3_p3 * R2 + (-6.0 / 7.0) * sqrt5 * Y1_p1 * R4) * V[3] + (1.2 * sqrt30 * Y3_p1 - 6.0 * sqrt2 * Y3_p3 + 4.8 * sqrt5 * Y1_p1 * R2) * V[4] - 4.5 * sqrt5 * Y1_p1 * V[5];
    d[5] = ((-1.0 / 21.0) * sqrt5 * Y5_p0 + (1.0 / 7.0) * sqrt7 * Y5_p4 + (2.0 / 15.0) * sqrt5 * Y3_p0 * R2 + (-3.0 / 35.0) * sqrt5 * Y1_p0 * R4) * V[0] + (-0.6 * sqrt5 * Y3_p0 + 0.6 * sqrt5 * Y1_p0 * R2) * V[1] - 0.75 * sqrt5 * Y1_p0 * V[2] + ((-10.0 / 21.0) * sqrt5 * Y5_p0 + (10.0 / 7.0) * sqrt7 * Y5_p4 + (4.0 / 3.0) * sqrt5 * Y3_p0 * R2 + (-6.0 / 7.0) * sqrt5 * Y1_p0 * R4) * V[3] + (-4.8 * sqrt5 * Y3_p0 + 4.8 * sqrt5 * Y1_p0 * R2) * V[4] - 4.5 * sqrt5 * Y1_p0 * V[5];
    d[2] = ((-2.0 / 21.0) * sqrt7 * Y5_n3 + (1.0 / 3.0) * Y3_n3 * R2) * V[0] - 1.5 * Y3_n3 * V[1] + ((-20.0 / 21.0) * sqrt7 * Y5_n3 + (10.0 / 3.0) * Y3_n3 * R2) * V[3] - 12.0 * Y3_n3 * V[4];
    d[3] = ((-1.0 / 21.0) * sqrt42 * Y5_n4 + (1.0 / 42.0) * sqrt14 * Y5_n2 + (-1.0 / 6.0) * sqrt2 * Y3_n2 * R2) * V[0] + 0.75 * sqrt2 * Y3_n2 * V[1] + ((-10.0 / 21.0) * sqrt42 * Y5_n4 + (5.0 / 21.0) * sqrt14 * Y5_n2 + (-5.0 / 3.0) * sqrt2 * Y3_n2 * R2) * V[3] + 6.0 * sqrt2 * Y3_n2 * V[4];
    d[4] = ((-1.0 / 21.0) * sqrt105 * Y5_n5 + (-1.0 / 42.0) * sqrt2 * Y5_n1 + (1.0 / 15.0) * sqrt5 * Y3_n1 * R2 + (-3.0 / 70.0) * sqrt30 * Y1_n1 * R4) * V[0] + (-0.3 * sqrt5 * Y3_n1 + 0.3 * sqrt30 * Y1_n1 * R2) * V[1] - 0.375 * sqrt30 * Y1_n1 * V[2] + ((-10.0 / 21.0) * sqrt105 * Y5_n5 + (-5.0 / 21.0) * sqrt2 * Y5_n1 + (2.0 / 3.0) * sqrt5 * Y3_n1 * R2 + (-3.0 / 7.0) * sqrt30 * Y1_n1 * R4) * V[3] + (-2.4 * sqrt5 * Y3_n1 + 2.4 * sqrt30 * Y1_n1 * R2) * V[4] - 2.25 * sqrt30 * Y1_n1 * V[5];
    d[1] = ((1.0 / 42.0) * sqrt14 * Y5_p2 + (1.0 / 21.0) * sqrt42 * Y5_p4 + (-1.0 / 6.0) * sqrt2 * Y3_p2 * R2) * V[0] + 0.75 * sqrt2 * Y3_p2 * V[1] + ((5.0 / 21.0) * sqrt14 * Y5_p2 + (10.0 / 21.0) * sqrt42 * Y5_p4 + (-5.0 / 3.0) * sqrt2 * Y3_p2 * R2) * V[3] + 6.0 * sqrt2 * Y3_p2 * V[4];
    d[0] = ((-1.0 / 42.0) * sqrt2 * Y5_p1 + (1.0 / 21.0) * sqrt105 * Y5_p5 + (1.0 / 15.0) * sqrt5 * Y3_p1 * R2 + (-3.0 / 70.0) * sqrt30 * Y1_p1 * R4) * V[0] + (-0.3 * sqrt5 * Y3_p1 + 0.3 * sqrt30 * Y1_p1 * R2) * V[1] - 0.375 * sqrt30 * Y1_p1 * V[2] + ((-5.0 / 21.0) * sqrt2 * Y5_p1 + (10.0 / 21.0) * sqrt105 * Y5_p5 + (2.0 / 3.0) * sqrt5 * Y3_p1 * R2 + (-3.0 / 7.0) * sqrt30 * Y1_p1 * R4) * V[3] + (-2.4 * sqrt5 * Y3_p1 + 2.4 * sqrt30 * Y1_p1 * R2) * V[4] - 2.25 * sqrt30 * Y1_p1 * V[5];
}

}  // namespace kinab