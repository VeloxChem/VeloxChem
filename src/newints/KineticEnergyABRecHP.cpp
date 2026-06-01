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

#include "KineticEnergyABRecHP.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_h_p(
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
    const auto sqrt10 = std::sqrt(10.0);
    const auto sqrt11 = std::sqrt(11.0);
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);
    const auto sqrt110 = std::sqrt(110.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α · β^5 · p^{-6} · (s|T|s)
    //   V[1] ↔ β^4 · p^{-5} · (s|T|s)
    //   V[2] ↔ α · β^5 · p^{-6} · ξ · (s|S|s)
    //   V[3] ↔ β^4 · p^{-5} · ξ · (s|S|s)
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

            V[0] += cab_tt * alpha * beta5 * pinv6;
            V[1] += cab_tt * beta4 * pinv5;
            V[2] += cab_ss * alpha * beta5 * pinv6 * xi;
            V[3] += cab_ss * beta4 * pinv5 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 11 × 3 spherical block ----
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
    d[16] = ((-6.0 / 11.0) * Y6_p0 + (-5.0 / 11.0) * Y4_p0 * R2) * V[0] + 2.5 * Y4_p0 * V[1] + ((-72.0 / 11.0) * Y6_p0 + (-60.0 / 11.0) * Y4_p0 * R2) * V[2] + 25.0 * Y4_p0 * V[3];
    d[17] = ((-1.0 / 11.0) * sqrt21 * Y6_p1 + (1.0 / 11.0) * sqrt10 * Y4_p1 * R2) * V[0] - 0.5 * sqrt10 * Y4_p1 * V[1] + ((-12.0 / 11.0) * sqrt21 * Y6_p1 + (12.0 / 11.0) * sqrt10 * Y4_p1 * R2) * V[2] - 5.0 * sqrt10 * Y4_p1 * V[3];
    d[15] = ((-1.0 / 11.0) * sqrt21 * Y6_n1 + (1.0 / 11.0) * sqrt10 * Y4_n1 * R2) * V[0] - 0.5 * sqrt10 * Y4_n1 * V[1] + ((-12.0 / 11.0) * sqrt21 * Y6_n1 + (12.0 / 11.0) * sqrt10 * Y4_n1 * R2) * V[2] - 5.0 * sqrt10 * Y4_n1 * V[3];
    d[19] = ((-1.0 / 11.0) * sqrt35 * Y6_p1 + (-2.0 / 11.0) * sqrt6 * Y4_p1 * R2) * V[0] + sqrt6 * Y4_p1 * V[1] + ((-12.0 / 11.0) * sqrt35 * Y6_p1 + (-24.0 / 11.0) * sqrt6 * Y4_p1 * R2) * V[2] + 10.0 * sqrt6 * Y4_p1 * V[3];
    d[20] = ((1.0 / 11.0) * sqrt15 * Y6_p0 + (-1.0 / 11.0) * sqrt14 * Y6_p2 + (-1.0 / 11.0) * sqrt15 * Y4_p0 * R2 + (1.0 / 11.0) * sqrt3 * Y4_p2 * R2) * V[0] + (0.5 * sqrt15 * Y4_p0 - 0.5 * sqrt3 * Y4_p2) * V[1] + ((12.0 / 11.0) * sqrt15 * Y6_p0 + (-12.0 / 11.0) * sqrt14 * Y6_p2 + (-12.0 / 11.0) * sqrt15 * Y4_p0 * R2 + (12.0 / 11.0) * sqrt3 * Y4_p2 * R2) * V[2] + (5.0 * sqrt15 * Y4_p0 - 5.0 * sqrt3 * Y4_p2) * V[3];
    d[18] = ((-1.0 / 11.0) * sqrt14 * Y6_n2 + (1.0 / 11.0) * sqrt3 * Y4_n2 * R2) * V[0] - 0.5 * sqrt3 * Y4_n2 * V[1] + ((-12.0 / 11.0) * sqrt14 * Y6_n2 + (12.0 / 11.0) * sqrt3 * Y4_n2 * R2) * V[2] - 5.0 * sqrt3 * Y4_n2 * V[3];
    d[22] = ((-4.0 / 11.0) * sqrt2 * Y6_p2 + (-1.0 / 11.0) * sqrt21 * Y4_p2 * R2) * V[0] + 0.5 * sqrt21 * Y4_p2 * V[1] + ((-48.0 / 11.0) * sqrt2 * Y6_p2 + (-12.0 / 11.0) * sqrt21 * Y4_p2 * R2) * V[2] + 5.0 * sqrt21 * Y4_p2 * V[3];
    d[23] = ((1.0 / 11.0) * sqrt5 * Y6_p1 + (-3.0 / 11.0) * sqrt2 * Y6_p3 + (-1.0 / 22.0) * sqrt42 * Y4_p1 * R2 + (1.0 / 22.0) * sqrt6 * Y4_p3 * R2) * V[0] + (0.25 * sqrt42 * Y4_p1 - 0.25 * sqrt6 * Y4_p3) * V[1] + ((12.0 / 11.0) * sqrt5 * Y6_p1 + (-36.0 / 11.0) * sqrt2 * Y6_p3 + (-6.0 / 11.0) * sqrt42 * Y4_p1 * R2 + (6.0 / 11.0) * sqrt6 * Y4_p3 * R2) * V[2] + (2.5 * sqrt42 * Y4_p1 - 2.5 * sqrt6 * Y4_p3) * V[3];
    d[21] = ((-3.0 / 11.0) * sqrt2 * Y6_n3 + (-1.0 / 11.0) * sqrt5 * Y6_n1 + (1.0 / 22.0) * sqrt6 * Y4_n3 * R2 + (1.0 / 22.0) * sqrt42 * Y4_n1 * R2) * V[0] + (-0.25 * sqrt6 * Y4_n3 - 0.25 * sqrt42 * Y4_n1) * V[1] + ((-36.0 / 11.0) * sqrt2 * Y6_n3 + (-12.0 / 11.0) * sqrt5 * Y6_n1 + (6.0 / 11.0) * sqrt6 * Y4_n3 * R2 + (6.0 / 11.0) * sqrt42 * Y4_n1 * R2) * V[2] + (-2.5 * sqrt6 * Y4_n3 - 2.5 * sqrt42 * Y4_n1) * V[3];
    d[25] = ((-3.0 / 11.0) * sqrt3 * Y6_p3 + (-4.0 / 11.0) * Y4_p3 * R2) * V[0] + 2.0 * Y4_p3 * V[1] + ((-36.0 / 11.0) * sqrt3 * Y6_p3 + (-48.0 / 11.0) * Y4_p3 * R2) * V[2] + 20.0 * Y4_p3 * V[3];
    d[26] = ((1.0 / 11.0) * sqrt3 * Y6_p2 + (-3.0 / 22.0) * sqrt10 * Y6_p4 + (-1.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (1.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[0] + (0.5 * sqrt14 * Y4_p2 - 0.25 * sqrt2 * Y4_p4) * V[1] + ((12.0 / 11.0) * sqrt3 * Y6_p2 + (-18.0 / 11.0) * sqrt10 * Y6_p4 + (-12.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (6.0 / 11.0) * sqrt2 * Y4_p4 * R2) * V[2] + (5.0 * sqrt14 * Y4_p2 - 2.5 * sqrt2 * Y4_p4) * V[3];
    d[24] = ((-3.0 / 22.0) * sqrt10 * Y6_n4 + (-1.0 / 11.0) * sqrt3 * Y6_n2 + (1.0 / 22.0) * sqrt2 * Y4_n4 * R2 + (1.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[0] + (-0.25 * sqrt2 * Y4_n4 - 0.5 * sqrt14 * Y4_n2) * V[1] + ((-18.0 / 11.0) * sqrt10 * Y6_n4 + (-12.0 / 11.0) * sqrt3 * Y6_n2 + (6.0 / 11.0) * sqrt2 * Y4_n4 * R2 + (12.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[2] + (-2.5 * sqrt2 * Y4_n4 - 5.0 * sqrt14 * Y4_n2) * V[3];
    d[28] = ((-2.0 / 11.0) * sqrt5 * Y6_p4 + (-3.0 / 11.0) * Y4_p4 * R2) * V[0] + 1.5 * Y4_p4 * V[1] + ((-24.0 / 11.0) * sqrt5 * Y6_p4 + (-36.0 / 11.0) * Y4_p4 * R2) * V[2] + 15.0 * Y4_p4 * V[3];
    d[29] = ((1.0 / 22.0) * sqrt6 * Y6_p3 + (-1.0 / 22.0) * sqrt110 * Y6_p5 + (-3.0 / 11.0) * sqrt2 * Y4_p3 * R2) * V[0] + 1.5 * sqrt2 * Y4_p3 * V[1] + ((6.0 / 11.0) * sqrt6 * Y6_p3 + (-6.0 / 11.0) * sqrt110 * Y6_p5 + (-36.0 / 11.0) * sqrt2 * Y4_p3 * R2) * V[2] + 15.0 * sqrt2 * Y4_p3 * V[3];
    d[27] = ((-1.0 / 22.0) * sqrt110 * Y6_n5 + (-1.0 / 22.0) * sqrt6 * Y6_n3 + (3.0 / 11.0) * sqrt2 * Y4_n3 * R2) * V[0] - 1.5 * sqrt2 * Y4_n3 * V[1] + ((-6.0 / 11.0) * sqrt110 * Y6_n5 + (-6.0 / 11.0) * sqrt6 * Y6_n3 + (36.0 / 11.0) * sqrt2 * Y4_n3 * R2) * V[2] - 15.0 * sqrt2 * Y4_n3 * V[3];
    d[31] = (-1.0 / 11.0) * sqrt11 * Y6_p5 * V[0] + (-12.0 / 11.0) * sqrt11 * Y6_p5 * V[2];
    d[32] = ((1.0 / 22.0) * sqrt2 * Y6_p4 + (-1.0 / 11.0) * sqrt33 * Y6_p6 + (-3.0 / 22.0) * sqrt10 * Y4_p4 * R2) * V[0] + 0.75 * sqrt10 * Y4_p4 * V[1] + ((6.0 / 11.0) * sqrt2 * Y6_p4 + (-12.0 / 11.0) * sqrt33 * Y6_p6 + (-18.0 / 11.0) * sqrt10 * Y4_p4 * R2) * V[2] + 7.5 * sqrt10 * Y4_p4 * V[3];
    d[30] = ((-1.0 / 11.0) * sqrt33 * Y6_n6 + (-1.0 / 22.0) * sqrt2 * Y6_n4 + (3.0 / 22.0) * sqrt10 * Y4_n4 * R2) * V[0] - 0.75 * sqrt10 * Y4_n4 * V[1] + ((-12.0 / 11.0) * sqrt33 * Y6_n6 + (-6.0 / 11.0) * sqrt2 * Y6_n4 + (18.0 / 11.0) * sqrt10 * Y4_n4 * R2) * V[2] - 7.5 * sqrt10 * Y4_n4 * V[3];
    d[13] = ((-1.0 / 11.0) * sqrt35 * Y6_n1 + (-2.0 / 11.0) * sqrt6 * Y4_n1 * R2) * V[0] + sqrt6 * Y4_n1 * V[1] + ((-12.0 / 11.0) * sqrt35 * Y6_n1 + (-24.0 / 11.0) * sqrt6 * Y4_n1 * R2) * V[2] + 10.0 * sqrt6 * Y4_n1 * V[3];
    d[14] = ((-1.0 / 11.0) * sqrt14 * Y6_n2 + (1.0 / 11.0) * sqrt3 * Y4_n2 * R2) * V[0] - 0.5 * sqrt3 * Y4_n2 * V[1] + ((-12.0 / 11.0) * sqrt14 * Y6_n2 + (12.0 / 11.0) * sqrt3 * Y4_n2 * R2) * V[2] - 5.0 * sqrt3 * Y4_n2 * V[3];
    d[12] = ((1.0 / 11.0) * sqrt15 * Y6_p0 + (1.0 / 11.0) * sqrt14 * Y6_p2 + (-1.0 / 11.0) * sqrt15 * Y4_p0 * R2 + (-1.0 / 11.0) * sqrt3 * Y4_p2 * R2) * V[0] + (0.5 * sqrt15 * Y4_p0 + 0.5 * sqrt3 * Y4_p2) * V[1] + ((12.0 / 11.0) * sqrt15 * Y6_p0 + (12.0 / 11.0) * sqrt14 * Y6_p2 + (-12.0 / 11.0) * sqrt15 * Y4_p0 * R2 + (-12.0 / 11.0) * sqrt3 * Y4_p2 * R2) * V[2] + (5.0 * sqrt15 * Y4_p0 + 5.0 * sqrt3 * Y4_p2) * V[3];
    d[10] = ((-4.0 / 11.0) * sqrt2 * Y6_n2 + (-1.0 / 11.0) * sqrt21 * Y4_n2 * R2) * V[0] + 0.5 * sqrt21 * Y4_n2 * V[1] + ((-48.0 / 11.0) * sqrt2 * Y6_n2 + (-12.0 / 11.0) * sqrt21 * Y4_n2 * R2) * V[2] + 5.0 * sqrt21 * Y4_n2 * V[3];
    d[11] = ((-3.0 / 11.0) * sqrt2 * Y6_n3 + (1.0 / 11.0) * sqrt5 * Y6_n1 + (1.0 / 22.0) * sqrt6 * Y4_n3 * R2 + (-1.0 / 22.0) * sqrt42 * Y4_n1 * R2) * V[0] + (-0.25 * sqrt6 * Y4_n3 + 0.25 * sqrt42 * Y4_n1) * V[1] + ((-36.0 / 11.0) * sqrt2 * Y6_n3 + (12.0 / 11.0) * sqrt5 * Y6_n1 + (6.0 / 11.0) * sqrt6 * Y4_n3 * R2 + (-6.0 / 11.0) * sqrt42 * Y4_n1 * R2) * V[2] + (-2.5 * sqrt6 * Y4_n3 + 2.5 * sqrt42 * Y4_n1) * V[3];
    d[9] = ((1.0 / 11.0) * sqrt5 * Y6_p1 + (3.0 / 11.0) * sqrt2 * Y6_p3 + (-1.0 / 22.0) * sqrt42 * Y4_p1 * R2 + (-1.0 / 22.0) * sqrt6 * Y4_p3 * R2) * V[0] + (0.25 * sqrt42 * Y4_p1 + 0.25 * sqrt6 * Y4_p3) * V[1] + ((12.0 / 11.0) * sqrt5 * Y6_p1 + (36.0 / 11.0) * sqrt2 * Y6_p3 + (-6.0 / 11.0) * sqrt42 * Y4_p1 * R2 + (-6.0 / 11.0) * sqrt6 * Y4_p3 * R2) * V[2] + (2.5 * sqrt42 * Y4_p1 + 2.5 * sqrt6 * Y4_p3) * V[3];
    d[7] = ((-3.0 / 11.0) * sqrt3 * Y6_n3 + (-4.0 / 11.0) * Y4_n3 * R2) * V[0] + 2.0 * Y4_n3 * V[1] + ((-36.0 / 11.0) * sqrt3 * Y6_n3 + (-48.0 / 11.0) * Y4_n3 * R2) * V[2] + 20.0 * Y4_n3 * V[3];
    d[8] = ((-3.0 / 22.0) * sqrt10 * Y6_n4 + (1.0 / 11.0) * sqrt3 * Y6_n2 + (1.0 / 22.0) * sqrt2 * Y4_n4 * R2 + (-1.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[0] + (-0.25 * sqrt2 * Y4_n4 + 0.5 * sqrt14 * Y4_n2) * V[1] + ((-18.0 / 11.0) * sqrt10 * Y6_n4 + (12.0 / 11.0) * sqrt3 * Y6_n2 + (6.0 / 11.0) * sqrt2 * Y4_n4 * R2 + (-12.0 / 11.0) * sqrt14 * Y4_n2 * R2) * V[2] + (-2.5 * sqrt2 * Y4_n4 + 5.0 * sqrt14 * Y4_n2) * V[3];
    d[6] = ((1.0 / 11.0) * sqrt3 * Y6_p2 + (3.0 / 22.0) * sqrt10 * Y6_p4 + (-1.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (-1.0 / 22.0) * sqrt2 * Y4_p4 * R2) * V[0] + (0.5 * sqrt14 * Y4_p2 + 0.25 * sqrt2 * Y4_p4) * V[1] + ((12.0 / 11.0) * sqrt3 * Y6_p2 + (18.0 / 11.0) * sqrt10 * Y6_p4 + (-12.0 / 11.0) * sqrt14 * Y4_p2 * R2 + (-6.0 / 11.0) * sqrt2 * Y4_p4 * R2) * V[2] + (5.0 * sqrt14 * Y4_p2 + 2.5 * sqrt2 * Y4_p4) * V[3];
    d[4] = ((-2.0 / 11.0) * sqrt5 * Y6_n4 + (-3.0 / 11.0) * Y4_n4 * R2) * V[0] + 1.5 * Y4_n4 * V[1] + ((-24.0 / 11.0) * sqrt5 * Y6_n4 + (-36.0 / 11.0) * Y4_n4 * R2) * V[2] + 15.0 * Y4_n4 * V[3];
    d[5] = ((-1.0 / 22.0) * sqrt110 * Y6_n5 + (1.0 / 22.0) * sqrt6 * Y6_n3 + (-3.0 / 11.0) * sqrt2 * Y4_n3 * R2) * V[0] + 1.5 * sqrt2 * Y4_n3 * V[1] + ((-6.0 / 11.0) * sqrt110 * Y6_n5 + (6.0 / 11.0) * sqrt6 * Y6_n3 + (-36.0 / 11.0) * sqrt2 * Y4_n3 * R2) * V[2] + 15.0 * sqrt2 * Y4_n3 * V[3];
    d[3] = ((1.0 / 22.0) * sqrt6 * Y6_p3 + (1.0 / 22.0) * sqrt110 * Y6_p5 + (-3.0 / 11.0) * sqrt2 * Y4_p3 * R2) * V[0] + 1.5 * sqrt2 * Y4_p3 * V[1] + ((6.0 / 11.0) * sqrt6 * Y6_p3 + (6.0 / 11.0) * sqrt110 * Y6_p5 + (-36.0 / 11.0) * sqrt2 * Y4_p3 * R2) * V[2] + 15.0 * sqrt2 * Y4_p3 * V[3];
    d[1] = (-1.0 / 11.0) * sqrt11 * Y6_n5 * V[0] + (-12.0 / 11.0) * sqrt11 * Y6_n5 * V[2];
    d[2] = ((-1.0 / 11.0) * sqrt33 * Y6_n6 + (1.0 / 22.0) * sqrt2 * Y6_n4 + (-3.0 / 22.0) * sqrt10 * Y4_n4 * R2) * V[0] + 0.75 * sqrt10 * Y4_n4 * V[1] + ((-12.0 / 11.0) * sqrt33 * Y6_n6 + (6.0 / 11.0) * sqrt2 * Y6_n4 + (-18.0 / 11.0) * sqrt10 * Y4_n4 * R2) * V[2] + 7.5 * sqrt10 * Y4_n4 * V[3];
    d[0] = ((1.0 / 22.0) * sqrt2 * Y6_p4 + (1.0 / 11.0) * sqrt33 * Y6_p6 + (-3.0 / 22.0) * sqrt10 * Y4_p4 * R2) * V[0] + 0.75 * sqrt10 * Y4_p4 * V[1] + ((6.0 / 11.0) * sqrt2 * Y6_p4 + (12.0 / 11.0) * sqrt33 * Y6_p6 + (-18.0 / 11.0) * sqrt10 * Y4_p4 * R2) * V[2] + 7.5 * sqrt10 * Y4_p4 * V[3];
}

}  // namespace kinab