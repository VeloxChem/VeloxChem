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

#include "OverlapABRecPI.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_p_i(
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
    const auto sqrt33 = std::sqrt(33.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt110 = std::sqrt(110.0);
    const auto sqrt154 = std::sqrt(154.0);
    const auto sqrt210 = std::sqrt(210.0);
    const auto sqrt462 = std::sqrt(462.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^6 · β · p^{-7} · (s|s)
    //   V[1] ↔ α^5 · p^{-6} · (s|s)
    const auto &exps_a  = bra.exponents();
    const auto &coefs_a = bra.normalization_factors();
    const auto &exps_b  = ket.exponents();
    const auto &coefs_b = ket.normalization_factors();

    const auto pi = mathconst::pi_value();

    std::array<double, 2> V = {0.0, 0.0};

    for (std::size_t i = 0; i < exps_a.size(); ++i)
    {
        const auto alpha = exps_a[i];
        const auto ca    = coefs_a[i];
        const auto alpha2 = alpha * alpha;
        const auto alpha3 = alpha2 * alpha;
        const auto alpha4 = alpha3 * alpha;
        const auto alpha5 = alpha4 * alpha;
        const auto alpha6 = alpha5 * alpha;
        for (std::size_t j = 0; j < exps_b.size(); ++j)
        {
            const auto beta = exps_b[j];
            const auto cb   = coefs_b[j];

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

            V[0] += cab_ss * alpha6 * beta * pinv7;
            V[1] += cab_ss * alpha5 * pinv6;
        }
    }

    // ---- Phase 3: fused M·V → 3 × 13 spherical block ----
    const auto Y1_n1 = harm::Y_ll_1_m_n1(AB_x, AB_y, AB_z);
    const auto Y1_p0 = harm::Y_ll_1_m_p0(AB_x, AB_y, AB_z);
    const auto Y1_p1 = harm::Y_ll_1_m_p1(AB_x, AB_y, AB_z);
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
    auto *d = buffer;
    d[19] = -Y1_p0 * Y6_p0 * V[0] + 3.0 * Y5_p0 * V[1];
    d[20] = -Y1_p0 * Y6_p1 * V[0] + 0.5 * sqrt35 * Y5_p1 * V[1];
    d[21] = -Y1_p0 * Y6_p2 * V[0] + 2.0 * sqrt2 * Y5_p2 * V[1];
    d[22] = -Y1_p0 * Y6_p3 * V[0] + 1.5 * sqrt3 * Y5_p3 * V[1];
    d[23] = -Y1_p0 * Y6_p4 * V[0] + sqrt5 * Y5_p4 * V[1];
    d[24] = -Y1_p0 * Y6_p5 * V[0] + 0.5 * sqrt11 * Y5_p5 * V[1];
    d[25] = -Y1_p0 * Y6_p6 * V[0];
    d[18] = -Y1_p0 * Y6_n1 * V[0] + 0.5 * sqrt35 * Y5_n1 * V[1];
    d[17] = -Y1_p0 * Y6_n2 * V[0] + 2.0 * sqrt2 * Y5_n2 * V[1];
    d[16] = -Y1_p0 * Y6_n3 * V[0] + 1.5 * sqrt3 * Y5_n3 * V[1];
    d[15] = -Y1_p0 * Y6_n4 * V[0] + sqrt5 * Y5_n4 * V[1];
    d[14] = -Y1_p0 * Y6_n5 * V[0] + 0.5 * sqrt11 * Y5_n5 * V[1];
    d[13] = -Y1_p0 * Y6_n6 * V[0];
    d[32] = -Y1_p1 * Y6_p0 * V[0] - 0.5 * sqrt15 * Y5_p1 * V[1];
    d[33] = -Y1_p1 * Y6_p1 * V[0] + (0.5 * sqrt21 * Y5_p0 - 0.5 * sqrt5 * Y5_p2) * V[1];
    d[34] = -Y1_p1 * Y6_p2 * V[0] + (0.5 * sqrt14 * Y5_p1 - 0.5 * sqrt3 * Y5_p3) * V[1];
    d[35] = -Y1_p1 * Y6_p3 * V[0] + (1.5 * sqrt2 * Y5_p2 - 0.25 * sqrt6 * Y5_p4) * V[1];
    d[36] = -Y1_p1 * Y6_p4 * V[0] + (0.75 * sqrt10 * Y5_p3 - 0.25 * sqrt2 * Y5_p5) * V[1];
    d[37] = -Y1_p1 * Y6_p5 * V[0] + 0.25 * sqrt110 * Y5_p4 * V[1];
    d[38] = -Y1_p1 * Y6_p6 * V[0] + 0.5 * sqrt33 * Y5_p5 * V[1];
    d[31] = -Y1_p1 * Y6_n1 * V[0] - 0.5 * sqrt5 * Y5_n2 * V[1];
    d[30] = -Y1_p1 * Y6_n2 * V[0] + (-0.5 * sqrt3 * Y5_n3 + 0.5 * sqrt14 * Y5_n1) * V[1];
    d[29] = -Y1_p1 * Y6_n3 * V[0] + (-0.25 * sqrt6 * Y5_n4 + 1.5 * sqrt2 * Y5_n2) * V[1];
    d[28] = -Y1_p1 * Y6_n4 * V[0] + (-0.25 * sqrt2 * Y5_n5 + 0.75 * sqrt10 * Y5_n3) * V[1];
    d[27] = -Y1_p1 * Y6_n5 * V[0] + 0.25 * sqrt110 * Y5_n4 * V[1];
    d[26] = -Y1_p1 * Y6_n6 * V[0] + 0.5 * sqrt33 * Y5_n5 * V[1];
    d[6] = -Y1_n1 * Y6_p0 * V[0] - 0.5 * sqrt15 * Y5_n1 * V[1];
    d[7] = -Y1_n1 * Y6_p1 * V[0] - 0.5 * sqrt5 * Y5_n2 * V[1];
    d[8] = -Y1_n1 * Y6_p2 * V[0] + (-0.5 * sqrt3 * Y5_n3 - 0.5 * sqrt14 * Y5_n1) * V[1];
    d[9] = -Y1_n1 * Y6_p3 * V[0] + (-0.25 * sqrt6 * Y5_n4 - 1.5 * sqrt2 * Y5_n2) * V[1];
    d[10] = -Y1_n1 * Y6_p4 * V[0] + (-0.25 * sqrt2 * Y5_n5 - 0.75 * sqrt10 * Y5_n3) * V[1];
    d[11] = -Y1_n1 * Y6_p5 * V[0] - 0.25 * sqrt110 * Y5_n4 * V[1];
    d[12] = -Y1_n1 * Y6_p6 * V[0] - 0.5 * sqrt33 * Y5_n5 * V[1];
    d[5] = -Y1_n1 * Y6_n1 * V[0] + (0.5 * sqrt21 * Y5_p0 + 0.5 * sqrt5 * Y5_p2) * V[1];
    d[4] = -Y1_n1 * Y6_n2 * V[0] + (0.5 * sqrt14 * Y5_p1 + 0.5 * sqrt3 * Y5_p3) * V[1];
    d[3] = -Y1_n1 * Y6_n3 * V[0] + (1.5 * sqrt2 * Y5_p2 + 0.25 * sqrt6 * Y5_p4) * V[1];
    d[2] = -Y1_n1 * Y6_n4 * V[0] + (0.75 * sqrt10 * Y5_p3 + 0.25 * sqrt2 * Y5_p5) * V[1];
    d[1] = -Y1_n1 * Y6_n5 * V[0] + 0.25 * sqrt110 * Y5_p4 * V[1];
    d[0] = -Y1_n1 * Y6_n6 * V[0] + 0.5 * sqrt33 * Y5_p5 * V[1];
}

}  // namespace ovlab