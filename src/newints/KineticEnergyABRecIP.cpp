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

#include "KineticEnergyABRecIP.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_i_p(
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
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt39 = std::sqrt(39.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt182 = std::sqrt(182.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt462 = std::sqrt(462.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α · β^6 · p^{-7} · (s|T|s)
    //   V[1] ↔ β^5 · p^{-6} · (s|T|s)
    //   V[2] ↔ α · β^6 · p^{-7} · ξ · (s|S|s)
    //   V[3] ↔ β^5 · p^{-6} · ξ · (s|S|s)
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
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha * beta6 * pinv7;
            V[1] += cab_tt * beta5 * pinv6;
            V[2] += cab_ss * alpha * beta6 * pinv7 * xi;
            V[3] += cab_ss * beta5 * pinv6 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 13 × 3 spherical block ----
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
    auto *d = buffer;
    d[19] = ((7.0 / 13.0) * Y7_p0 + (6.0 / 13.0) * Y5_p0 * R2) * V[0] - 3.0 * Y5_p0 * V[1] + ((98.0 / 13.0) * Y7_p0 + (84.0 / 13.0) * Y5_p0 * R2) * V[2] - 36.0 * Y5_p0 * V[3];
    d[20] = ((2.0 / 13.0) * sqrt7 * Y7_p1 + (-1.0 / 13.0) * sqrt15 * Y5_p1 * R2) * V[0] + 0.5 * sqrt15 * Y5_p1 * V[1] + ((28.0 / 13.0) * sqrt7 * Y7_p1 + (-14.0 / 13.0) * sqrt15 * Y5_p1 * R2) * V[2] + 6.0 * sqrt15 * Y5_p1 * V[3];
    d[18] = ((2.0 / 13.0) * sqrt7 * Y7_n1 + (-1.0 / 13.0) * sqrt15 * Y5_n1 * R2) * V[0] + 0.5 * sqrt15 * Y5_n1 * V[1] + ((28.0 / 13.0) * sqrt7 * Y7_n1 + (-14.0 / 13.0) * sqrt15 * Y5_n1 * R2) * V[2] + 6.0 * sqrt15 * Y5_n1 * V[3];
    d[22] = ((4.0 / 13.0) * sqrt3 * Y7_p1 + (1.0 / 13.0) * sqrt35 * Y5_p1 * R2) * V[0] - 0.5 * sqrt35 * Y5_p1 * V[1] + ((56.0 / 13.0) * sqrt3 * Y7_p1 + (14.0 / 13.0) * sqrt35 * Y5_p1 * R2) * V[2] - 6.0 * sqrt35 * Y5_p1 * V[3];
    d[23] = ((-1.0 / 13.0) * sqrt21 * Y7_p0 + (3.0 / 13.0) * sqrt2 * Y7_p2 + (1.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (-1.0 / 13.0) * sqrt5 * Y5_p2 * R2) * V[0] + (-0.5 * sqrt21 * Y5_p0 + 0.5 * sqrt5 * Y5_p2) * V[1] + ((-14.0 / 13.0) * sqrt21 * Y7_p0 + (42.0 / 13.0) * sqrt2 * Y7_p2 + (14.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (-14.0 / 13.0) * sqrt5 * Y5_p2 * R2) * V[2] + (-6.0 * sqrt21 * Y5_p0 + 6.0 * sqrt5 * Y5_p2) * V[3];
    d[21] = ((3.0 / 13.0) * sqrt2 * Y7_n2 + (-1.0 / 13.0) * sqrt5 * Y5_n2 * R2) * V[0] + 0.5 * sqrt5 * Y5_n2 * V[1] + ((42.0 / 13.0) * sqrt2 * Y7_n2 + (-14.0 / 13.0) * sqrt5 * Y5_n2 * R2) * V[2] + 6.0 * sqrt5 * Y5_n2 * V[3];
    d[25] = ((3.0 / 13.0) * sqrt5 * Y7_p2 + (4.0 / 13.0) * sqrt2 * Y5_p2 * R2) * V[0] - 2.0 * sqrt2 * Y5_p2 * V[1] + ((42.0 / 13.0) * sqrt5 * Y7_p2 + (56.0 / 13.0) * sqrt2 * Y5_p2 * R2) * V[2] - 24.0 * sqrt2 * Y5_p2 * V[3];
    d[26] = ((-1.0 / 26.0) * sqrt30 * Y7_p1 + (3.0 / 26.0) * sqrt10 * Y7_p3 + (1.0 / 13.0) * sqrt14 * Y5_p1 * R2 + (-1.0 / 13.0) * sqrt3 * Y5_p3 * R2) * V[0] + (-0.5 * sqrt14 * Y5_p1 + 0.5 * sqrt3 * Y5_p3) * V[1] + ((-7.0 / 13.0) * sqrt30 * Y7_p1 + (21.0 / 13.0) * sqrt10 * Y7_p3 + (14.0 / 13.0) * sqrt14 * Y5_p1 * R2 + (-14.0 / 13.0) * sqrt3 * Y5_p3 * R2) * V[2] + (-6.0 * sqrt14 * Y5_p1 + 6.0 * sqrt3 * Y5_p3) * V[3];
    d[24] = ((3.0 / 26.0) * sqrt10 * Y7_n3 + (1.0 / 26.0) * sqrt30 * Y7_n1 + (-1.0 / 13.0) * sqrt3 * Y5_n3 * R2 + (-1.0 / 13.0) * sqrt14 * Y5_n1 * R2) * V[0] + (0.5 * sqrt3 * Y5_n3 + 0.5 * sqrt14 * Y5_n1) * V[1] + ((21.0 / 13.0) * sqrt10 * Y7_n3 + (7.0 / 13.0) * sqrt30 * Y7_n1 + (-14.0 / 13.0) * sqrt3 * Y5_n3 * R2 + (-14.0 / 13.0) * sqrt14 * Y5_n1 * R2) * V[2] + (6.0 * sqrt3 * Y5_n3 + 6.0 * sqrt14 * Y5_n1) * V[3];
    d[28] = ((2.0 / 13.0) * sqrt10 * Y7_p3 + (3.0 / 13.0) * sqrt3 * Y5_p3 * R2) * V[0] - 1.5 * sqrt3 * Y5_p3 * V[1] + ((28.0 / 13.0) * sqrt10 * Y7_p3 + (42.0 / 13.0) * sqrt3 * Y5_p3 * R2) * V[2] - 18.0 * sqrt3 * Y5_p3 * V[3];
    d[29] = ((-1.0 / 13.0) * sqrt5 * Y7_p2 + (1.0 / 26.0) * sqrt110 * Y7_p4 + (3.0 / 13.0) * sqrt2 * Y5_p2 * R2 + (-1.0 / 26.0) * sqrt6 * Y5_p4 * R2) * V[0] + (-1.5 * sqrt2 * Y5_p2 + 0.25 * sqrt6 * Y5_p4) * V[1] + ((-14.0 / 13.0) * sqrt5 * Y7_p2 + (7.0 / 13.0) * sqrt110 * Y7_p4 + (42.0 / 13.0) * sqrt2 * Y5_p2 * R2 + (-7.0 / 13.0) * sqrt6 * Y5_p4 * R2) * V[2] + (-18.0 * sqrt2 * Y5_p2 + 3.0 * sqrt6 * Y5_p4) * V[3];
    d[27] = ((1.0 / 26.0) * sqrt110 * Y7_n4 + (1.0 / 13.0) * sqrt5 * Y7_n2 + (-1.0 / 26.0) * sqrt6 * Y5_n4 * R2 + (-3.0 / 13.0) * sqrt2 * Y5_n2 * R2) * V[0] + (0.25 * sqrt6 * Y5_n4 + 1.5 * sqrt2 * Y5_n2) * V[1] + ((7.0 / 13.0) * sqrt110 * Y7_n4 + (14.0 / 13.0) * sqrt5 * Y7_n2 + (-7.0 / 13.0) * sqrt6 * Y5_n4 * R2 + (-42.0 / 13.0) * sqrt2 * Y5_n2 * R2) * V[2] + (3.0 * sqrt6 * Y5_n4 + 18.0 * sqrt2 * Y5_n2) * V[3];
    d[31] = ((1.0 / 13.0) * sqrt33 * Y7_p4 + (2.0 / 13.0) * sqrt5 * Y5_p4 * R2) * V[0] - sqrt5 * Y5_p4 * V[1] + ((14.0 / 13.0) * sqrt33 * Y7_p4 + (28.0 / 13.0) * sqrt5 * Y5_p4 * R2) * V[2] - 12.0 * sqrt5 * Y5_p4 * V[3];
    d[32] = ((-1.0 / 13.0) * sqrt3 * Y7_p3 + (1.0 / 13.0) * sqrt33 * Y7_p5 + (3.0 / 26.0) * sqrt10 * Y5_p3 * R2 + (-1.0 / 26.0) * sqrt2 * Y5_p5 * R2) * V[0] + (-0.75 * sqrt10 * Y5_p3 + 0.25 * sqrt2 * Y5_p5) * V[1] + ((-14.0 / 13.0) * sqrt3 * Y7_p3 + (14.0 / 13.0) * sqrt33 * Y7_p5 + (21.0 / 13.0) * sqrt10 * Y5_p3 * R2 + (-7.0 / 13.0) * sqrt2 * Y5_p5 * R2) * V[2] + (-9.0 * sqrt10 * Y5_p3 + 3.0 * sqrt2 * Y5_p5) * V[3];
    d[30] = ((1.0 / 13.0) * sqrt33 * Y7_n5 + (1.0 / 13.0) * sqrt3 * Y7_n3 + (-1.0 / 26.0) * sqrt2 * Y5_n5 * R2 + (-3.0 / 26.0) * sqrt10 * Y5_n3 * R2) * V[0] + (0.25 * sqrt2 * Y5_n5 + 0.75 * sqrt10 * Y5_n3) * V[1] + ((14.0 / 13.0) * sqrt33 * Y7_n5 + (14.0 / 13.0) * sqrt3 * Y7_n3 + (-7.0 / 13.0) * sqrt2 * Y5_n5 * R2 + (-21.0 / 13.0) * sqrt10 * Y5_n3 * R2) * V[2] + (3.0 * sqrt2 * Y5_n5 + 9.0 * sqrt10 * Y5_n3) * V[3];
    d[34] = ((2.0 / 13.0) * sqrt6 * Y7_p5 + (1.0 / 13.0) * sqrt11 * Y5_p5 * R2) * V[0] - 0.5 * sqrt11 * Y5_p5 * V[1] + ((28.0 / 13.0) * sqrt6 * Y7_p5 + (14.0 / 13.0) * sqrt11 * Y5_p5 * R2) * V[2] - 6.0 * sqrt11 * Y5_p5 * V[3];
    d[35] = ((-1.0 / 26.0) * sqrt6 * Y7_p4 + (1.0 / 13.0) * sqrt39 * Y7_p6 + (1.0 / 26.0) * sqrt110 * Y5_p4 * R2) * V[0] - 0.25 * sqrt110 * Y5_p4 * V[1] + ((-7.0 / 13.0) * sqrt6 * Y7_p4 + (14.0 / 13.0) * sqrt39 * Y7_p6 + (7.0 / 13.0) * sqrt110 * Y5_p4 * R2) * V[2] - 3.0 * sqrt110 * Y5_p4 * V[3];
    d[33] = ((1.0 / 13.0) * sqrt39 * Y7_n6 + (1.0 / 26.0) * sqrt6 * Y7_n4 + (-1.0 / 26.0) * sqrt110 * Y5_n4 * R2) * V[0] + 0.25 * sqrt110 * Y5_n4 * V[1] + ((14.0 / 13.0) * sqrt39 * Y7_n6 + (7.0 / 13.0) * sqrt6 * Y7_n4 + (-7.0 / 13.0) * sqrt110 * Y5_n4 * R2) * V[2] + 3.0 * sqrt110 * Y5_n4 * V[3];
    d[37] = (1.0 / 13.0) * sqrt13 * Y7_p6 * V[0] + (14.0 / 13.0) * sqrt13 * Y7_p6 * V[2];
    d[38] = ((-1.0 / 26.0) * sqrt2 * Y7_p5 + (1.0 / 26.0) * sqrt182 * Y7_p7 + (1.0 / 13.0) * sqrt33 * Y5_p5 * R2) * V[0] - 0.5 * sqrt33 * Y5_p5 * V[1] + ((-7.0 / 13.0) * sqrt2 * Y7_p5 + (7.0 / 13.0) * sqrt182 * Y7_p7 + (14.0 / 13.0) * sqrt33 * Y5_p5 * R2) * V[2] - 6.0 * sqrt33 * Y5_p5 * V[3];
    d[36] = ((1.0 / 26.0) * sqrt182 * Y7_n7 + (1.0 / 26.0) * sqrt2 * Y7_n5 + (-1.0 / 13.0) * sqrt33 * Y5_n5 * R2) * V[0] + 0.5 * sqrt33 * Y5_n5 * V[1] + ((7.0 / 13.0) * sqrt182 * Y7_n7 + (7.0 / 13.0) * sqrt2 * Y7_n5 + (-14.0 / 13.0) * sqrt33 * Y5_n5 * R2) * V[2] + 6.0 * sqrt33 * Y5_n5 * V[3];
    d[16] = ((4.0 / 13.0) * sqrt3 * Y7_n1 + (1.0 / 13.0) * sqrt35 * Y5_n1 * R2) * V[0] - 0.5 * sqrt35 * Y5_n1 * V[1] + ((56.0 / 13.0) * sqrt3 * Y7_n1 + (14.0 / 13.0) * sqrt35 * Y5_n1 * R2) * V[2] - 6.0 * sqrt35 * Y5_n1 * V[3];
    d[17] = ((3.0 / 13.0) * sqrt2 * Y7_n2 + (-1.0 / 13.0) * sqrt5 * Y5_n2 * R2) * V[0] + 0.5 * sqrt5 * Y5_n2 * V[1] + ((42.0 / 13.0) * sqrt2 * Y7_n2 + (-14.0 / 13.0) * sqrt5 * Y5_n2 * R2) * V[2] + 6.0 * sqrt5 * Y5_n2 * V[3];
    d[15] = ((-1.0 / 13.0) * sqrt21 * Y7_p0 + (-3.0 / 13.0) * sqrt2 * Y7_p2 + (1.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (1.0 / 13.0) * sqrt5 * Y5_p2 * R2) * V[0] + (-0.5 * sqrt21 * Y5_p0 - 0.5 * sqrt5 * Y5_p2) * V[1] + ((-14.0 / 13.0) * sqrt21 * Y7_p0 + (-42.0 / 13.0) * sqrt2 * Y7_p2 + (14.0 / 13.0) * sqrt21 * Y5_p0 * R2 + (14.0 / 13.0) * sqrt5 * Y5_p2 * R2) * V[2] + (-6.0 * sqrt21 * Y5_p0 - 6.0 * sqrt5 * Y5_p2) * V[3];
    d[13] = ((3.0 / 13.0) * sqrt5 * Y7_n2 + (4.0 / 13.0) * sqrt2 * Y5_n2 * R2) * V[0] - 2.0 * sqrt2 * Y5_n2 * V[1] + ((42.0 / 13.0) * sqrt5 * Y7_n2 + (56.0 / 13.0) * sqrt2 * Y5_n2 * R2) * V[2] - 24.0 * sqrt2 * Y5_n2 * V[3];
    d[14] = ((3.0 / 26.0) * sqrt10 * Y7_n3 + (-1.0 / 26.0) * sqrt30 * Y7_n1 + (-1.0 / 13.0) * sqrt3 * Y5_n3 * R2 + (1.0 / 13.0) * sqrt14 * Y5_n1 * R2) * V[0] + (0.5 * sqrt3 * Y5_n3 - 0.5 * sqrt14 * Y5_n1) * V[1] + ((21.0 / 13.0) * sqrt10 * Y7_n3 + (-7.0 / 13.0) * sqrt30 * Y7_n1 + (-14.0 / 13.0) * sqrt3 * Y5_n3 * R2 + (14.0 / 13.0) * sqrt14 * Y5_n1 * R2) * V[2] + (6.0 * sqrt3 * Y5_n3 - 6.0 * sqrt14 * Y5_n1) * V[3];
    d[12] = ((-1.0 / 26.0) * sqrt30 * Y7_p1 + (-3.0 / 26.0) * sqrt10 * Y7_p3 + (1.0 / 13.0) * sqrt14 * Y5_p1 * R2 + (1.0 / 13.0) * sqrt3 * Y5_p3 * R2) * V[0] + (-0.5 * sqrt14 * Y5_p1 - 0.5 * sqrt3 * Y5_p3) * V[1] + ((-7.0 / 13.0) * sqrt30 * Y7_p1 + (-21.0 / 13.0) * sqrt10 * Y7_p3 + (14.0 / 13.0) * sqrt14 * Y5_p1 * R2 + (14.0 / 13.0) * sqrt3 * Y5_p3 * R2) * V[2] + (-6.0 * sqrt14 * Y5_p1 - 6.0 * sqrt3 * Y5_p3) * V[3];
    d[10] = ((2.0 / 13.0) * sqrt10 * Y7_n3 + (3.0 / 13.0) * sqrt3 * Y5_n3 * R2) * V[0] - 1.5 * sqrt3 * Y5_n3 * V[1] + ((28.0 / 13.0) * sqrt10 * Y7_n3 + (42.0 / 13.0) * sqrt3 * Y5_n3 * R2) * V[2] - 18.0 * sqrt3 * Y5_n3 * V[3];
    d[11] = ((1.0 / 26.0) * sqrt110 * Y7_n4 + (-1.0 / 13.0) * sqrt5 * Y7_n2 + (-1.0 / 26.0) * sqrt6 * Y5_n4 * R2 + (3.0 / 13.0) * sqrt2 * Y5_n2 * R2) * V[0] + (0.25 * sqrt6 * Y5_n4 - 1.5 * sqrt2 * Y5_n2) * V[1] + ((7.0 / 13.0) * sqrt110 * Y7_n4 + (-14.0 / 13.0) * sqrt5 * Y7_n2 + (-7.0 / 13.0) * sqrt6 * Y5_n4 * R2 + (42.0 / 13.0) * sqrt2 * Y5_n2 * R2) * V[2] + (3.0 * sqrt6 * Y5_n4 - 18.0 * sqrt2 * Y5_n2) * V[3];
    d[9] = ((-1.0 / 13.0) * sqrt5 * Y7_p2 + (-1.0 / 26.0) * sqrt110 * Y7_p4 + (3.0 / 13.0) * sqrt2 * Y5_p2 * R2 + (1.0 / 26.0) * sqrt6 * Y5_p4 * R2) * V[0] + (-1.5 * sqrt2 * Y5_p2 - 0.25 * sqrt6 * Y5_p4) * V[1] + ((-14.0 / 13.0) * sqrt5 * Y7_p2 + (-7.0 / 13.0) * sqrt110 * Y7_p4 + (42.0 / 13.0) * sqrt2 * Y5_p2 * R2 + (7.0 / 13.0) * sqrt6 * Y5_p4 * R2) * V[2] + (-18.0 * sqrt2 * Y5_p2 - 3.0 * sqrt6 * Y5_p4) * V[3];
    d[7] = ((1.0 / 13.0) * sqrt33 * Y7_n4 + (2.0 / 13.0) * sqrt5 * Y5_n4 * R2) * V[0] - sqrt5 * Y5_n4 * V[1] + ((14.0 / 13.0) * sqrt33 * Y7_n4 + (28.0 / 13.0) * sqrt5 * Y5_n4 * R2) * V[2] - 12.0 * sqrt5 * Y5_n4 * V[3];
    d[8] = ((1.0 / 13.0) * sqrt33 * Y7_n5 + (-1.0 / 13.0) * sqrt3 * Y7_n3 + (-1.0 / 26.0) * sqrt2 * Y5_n5 * R2 + (3.0 / 26.0) * sqrt10 * Y5_n3 * R2) * V[0] + (0.25 * sqrt2 * Y5_n5 - 0.75 * sqrt10 * Y5_n3) * V[1] + ((14.0 / 13.0) * sqrt33 * Y7_n5 + (-14.0 / 13.0) * sqrt3 * Y7_n3 + (-7.0 / 13.0) * sqrt2 * Y5_n5 * R2 + (21.0 / 13.0) * sqrt10 * Y5_n3 * R2) * V[2] + (3.0 * sqrt2 * Y5_n5 - 9.0 * sqrt10 * Y5_n3) * V[3];
    d[6] = ((-1.0 / 13.0) * sqrt3 * Y7_p3 + (-1.0 / 13.0) * sqrt33 * Y7_p5 + (3.0 / 26.0) * sqrt10 * Y5_p3 * R2 + (1.0 / 26.0) * sqrt2 * Y5_p5 * R2) * V[0] + (-0.75 * sqrt10 * Y5_p3 - 0.25 * sqrt2 * Y5_p5) * V[1] + ((-14.0 / 13.0) * sqrt3 * Y7_p3 + (-14.0 / 13.0) * sqrt33 * Y7_p5 + (21.0 / 13.0) * sqrt10 * Y5_p3 * R2 + (7.0 / 13.0) * sqrt2 * Y5_p5 * R2) * V[2] + (-9.0 * sqrt10 * Y5_p3 - 3.0 * sqrt2 * Y5_p5) * V[3];
    d[4] = ((2.0 / 13.0) * sqrt6 * Y7_n5 + (1.0 / 13.0) * sqrt11 * Y5_n5 * R2) * V[0] - 0.5 * sqrt11 * Y5_n5 * V[1] + ((28.0 / 13.0) * sqrt6 * Y7_n5 + (14.0 / 13.0) * sqrt11 * Y5_n5 * R2) * V[2] - 6.0 * sqrt11 * Y5_n5 * V[3];
    d[5] = ((1.0 / 13.0) * sqrt39 * Y7_n6 + (-1.0 / 26.0) * sqrt6 * Y7_n4 + (1.0 / 26.0) * sqrt110 * Y5_n4 * R2) * V[0] - 0.25 * sqrt110 * Y5_n4 * V[1] + ((14.0 / 13.0) * sqrt39 * Y7_n6 + (-7.0 / 13.0) * sqrt6 * Y7_n4 + (7.0 / 13.0) * sqrt110 * Y5_n4 * R2) * V[2] - 3.0 * sqrt110 * Y5_n4 * V[3];
    d[3] = ((-1.0 / 26.0) * sqrt6 * Y7_p4 + (-1.0 / 13.0) * sqrt39 * Y7_p6 + (1.0 / 26.0) * sqrt110 * Y5_p4 * R2) * V[0] - 0.25 * sqrt110 * Y5_p4 * V[1] + ((-7.0 / 13.0) * sqrt6 * Y7_p4 + (-14.0 / 13.0) * sqrt39 * Y7_p6 + (7.0 / 13.0) * sqrt110 * Y5_p4 * R2) * V[2] - 3.0 * sqrt110 * Y5_p4 * V[3];
    d[1] = (1.0 / 13.0) * sqrt13 * Y7_n6 * V[0] + (14.0 / 13.0) * sqrt13 * Y7_n6 * V[2];
    d[2] = ((1.0 / 26.0) * sqrt182 * Y7_n7 + (-1.0 / 26.0) * sqrt2 * Y7_n5 + (1.0 / 13.0) * sqrt33 * Y5_n5 * R2) * V[0] - 0.5 * sqrt33 * Y5_n5 * V[1] + ((7.0 / 13.0) * sqrt182 * Y7_n7 + (-7.0 / 13.0) * sqrt2 * Y7_n5 + (14.0 / 13.0) * sqrt33 * Y5_n5 * R2) * V[2] - 6.0 * sqrt33 * Y5_n5 * V[3];
    d[0] = ((-1.0 / 26.0) * sqrt2 * Y7_p5 + (-1.0 / 26.0) * sqrt182 * Y7_p7 + (1.0 / 13.0) * sqrt33 * Y5_p5 * R2) * V[0] - 0.5 * sqrt33 * Y5_p5 * V[1] + ((-7.0 / 13.0) * sqrt2 * Y7_p5 + (-7.0 / 13.0) * sqrt182 * Y7_p7 + (14.0 / 13.0) * sqrt33 * Y5_p5 * R2) * V[2] - 6.0 * sqrt33 * Y5_p5 * V[3];
}

}  // namespace kinab