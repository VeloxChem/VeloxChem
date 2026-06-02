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

#include "ElectronRepulsionABRecHP.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "BoysFunction.hpp"
#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace eri2cab {

auto eri2c_h_p(
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
    const auto sqrt6 = std::sqrt(6.0);
    const auto sqrt10 = std::sqrt(10.0);
    const auto sqrt14 = std::sqrt(14.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt42 = std::sqrt(42.0);
    const auto sqrt70 = std::sqrt(70.0);
    const auto sqrt105 = std::sqrt(105.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ β^4 · p^{-5} · (s|s)^(5)
    //   V[1] ↔ α · β^5 · p^{-6} · (s|s)^(6)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi      = mathconst::pi_value();
    const auto two_pi52 = 2.0 * pi * pi * std::sqrt(pi);

    std::array<double, 2> V = {0.0, 0.0};

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

            const auto p    = alpha + beta;
            const auto pinv = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;
            const auto xi   = alpha * beta * pinv;

            const auto T = xi * R2;
            std::array<double, 7> F;
            newints::boys_function<7>(T, F);

            const auto pref     = two_pi52 / (alpha * beta * std::sqrt(p));
            const auto cab_pref = ca * cb * pref;

            V[0] += cab_pref * beta4 * pinv5 * F[5];
            V[1] += cab_pref * alpha * beta5 * pinv6 * F[6];
        }
    }

    // ---- Phase 3: fused M·V → 11 × 3 spherical block ----
    const auto Y1_n1 = harm::Y_ll_1_m_n1(AB_x, AB_y, AB_z);
    const auto Y1_p0 = harm::Y_ll_1_m_p0(AB_x, AB_y, AB_z);
    const auto Y1_p1 = harm::Y_ll_1_m_p1(AB_x, AB_y, AB_z);
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
    auto *d = buffer;
    d[16] = -Y5_p0 * Y1_p0 * V[1] + 2.5 * Y4_p0 * V[0];
    d[17] = -Y5_p0 * Y1_p1 * V[1] - 0.5 * sqrt10 * Y4_p1 * V[0];
    d[15] = -Y5_p0 * Y1_n1 * V[1] - 0.5 * sqrt10 * Y4_n1 * V[0];
    d[19] = -Y5_p1 * Y1_p0 * V[1] + sqrt6 * Y4_p1 * V[0];
    d[20] = -Y5_p1 * Y1_p1 * V[1] + (0.5 * sqrt15 * Y4_p0 - 0.5 * sqrt3 * Y4_p2) * V[0];
    d[18] = -Y5_p1 * Y1_n1 * V[1] - 0.5 * sqrt3 * Y4_n2 * V[0];
    d[22] = -Y5_p2 * Y1_p0 * V[1] + 0.5 * sqrt21 * Y4_p2 * V[0];
    d[23] = -Y5_p2 * Y1_p1 * V[1] + (0.25 * sqrt42 * Y4_p1 - 0.25 * sqrt6 * Y4_p3) * V[0];
    d[21] = -Y5_p2 * Y1_n1 * V[1] + (-0.25 * sqrt6 * Y4_n3 - 0.25 * sqrt42 * Y4_n1) * V[0];
    d[25] = -Y5_p3 * Y1_p0 * V[1] + 2.0 * Y4_p3 * V[0];
    d[26] = -Y5_p3 * Y1_p1 * V[1] + (0.5 * sqrt14 * Y4_p2 - 0.25 * sqrt2 * Y4_p4) * V[0];
    d[24] = -Y5_p3 * Y1_n1 * V[1] + (-0.25 * sqrt2 * Y4_n4 - 0.5 * sqrt14 * Y4_n2) * V[0];
    d[28] = -Y5_p4 * Y1_p0 * V[1] + 1.5 * Y4_p4 * V[0];
    d[29] = -Y5_p4 * Y1_p1 * V[1] + 1.5 * sqrt2 * Y4_p3 * V[0];
    d[27] = -Y5_p4 * Y1_n1 * V[1] - 1.5 * sqrt2 * Y4_n3 * V[0];
    d[31] = -Y5_p5 * Y1_p0 * V[1];
    d[32] = -Y5_p5 * Y1_p1 * V[1] + 0.75 * sqrt10 * Y4_p4 * V[0];
    d[30] = -Y5_p5 * Y1_n1 * V[1] - 0.75 * sqrt10 * Y4_n4 * V[0];
    d[13] = -Y5_n1 * Y1_p0 * V[1] + sqrt6 * Y4_n1 * V[0];
    d[14] = -Y5_n1 * Y1_p1 * V[1] - 0.5 * sqrt3 * Y4_n2 * V[0];
    d[12] = -Y5_n1 * Y1_n1 * V[1] + (0.5 * sqrt15 * Y4_p0 + 0.5 * sqrt3 * Y4_p2) * V[0];
    d[10] = -Y5_n2 * Y1_p0 * V[1] + 0.5 * sqrt21 * Y4_n2 * V[0];
    d[11] = -Y5_n2 * Y1_p1 * V[1] + (-0.25 * sqrt6 * Y4_n3 + 0.25 * sqrt42 * Y4_n1) * V[0];
    d[9] = -Y5_n2 * Y1_n1 * V[1] + (0.25 * sqrt42 * Y4_p1 + 0.25 * sqrt6 * Y4_p3) * V[0];
    d[7] = -Y5_n3 * Y1_p0 * V[1] + 2.0 * Y4_n3 * V[0];
    d[8] = -Y5_n3 * Y1_p1 * V[1] + (-0.25 * sqrt2 * Y4_n4 + 0.5 * sqrt14 * Y4_n2) * V[0];
    d[6] = -Y5_n3 * Y1_n1 * V[1] + (0.5 * sqrt14 * Y4_p2 + 0.25 * sqrt2 * Y4_p4) * V[0];
    d[4] = -Y5_n4 * Y1_p0 * V[1] + 1.5 * Y4_n4 * V[0];
    d[5] = -Y5_n4 * Y1_p1 * V[1] + 1.5 * sqrt2 * Y4_n3 * V[0];
    d[3] = -Y5_n4 * Y1_n1 * V[1] + 1.5 * sqrt2 * Y4_p3 * V[0];
    d[1] = -Y5_n5 * Y1_p0 * V[1];
    d[2] = -Y5_n5 * Y1_p1 * V[1] + 0.75 * sqrt10 * Y4_n4 * V[0];
    d[0] = -Y5_n5 * Y1_n1 * V[1] + 0.75 * sqrt10 * Y4_p4 * V[0];
}

}  // namespace eri2cab