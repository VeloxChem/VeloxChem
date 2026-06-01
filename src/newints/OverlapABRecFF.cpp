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

#include "OverlapABRecFF.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace ovlab {

auto overlap_f_f(
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
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt21 = std::sqrt(21.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^3 · β^3 · p^{-6} · (s|s)
    //   V[1] ↔ α^2 · β^2 · p^{-5} · (s|s)
    //   V[2] ↔ α · β · p^{-4} · (s|s)
    //   V[3] ↔ p^{-3} · (s|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto pinv6 = pinv5 * pinv;

            const auto fpi    = pi * pinv;
            const auto ss     = fpi * std::sqrt(fpi)
                              * std::exp(-alpha * beta * pinv * R2);
            const auto cab_ss = ca * cb * ss;

            V[0] += cab_ss * alpha3 * beta3 * pinv6;
            V[1] += cab_ss * alpha2 * beta2 * pinv5;
            V[2] += cab_ss * alpha * beta * pinv4;
            V[3] += cab_ss * pinv3;
        }
    }

    // ---- Phase 3: fused M·V → 7 × 7 spherical block ----
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
    const auto R4 = R2 * R2;
    newints::Block out{7, 7, std::vector<double>(49, 0.0)};
    auto *d = out.data.data();
    d[24] = -Y3_p0 * Y3_p0 * V[0] + ((9.0 / 7.0) * Y4_p0 + (12.0 / 7.0) * Y2_p0 * R2 + 1.5 * R4) * V[1] + (-3.0 * Y2_p0 - 3.75 * R2) * V[2] + 1.875 * V[3];
    d[25] = d[31] = -Y3_p0 * Y3_p1 * V[0] + ((3.0 / 14.0) * sqrt15 * Y4_p1 + (3.0 / 7.0) * sqrt2 * Y2_p1 * R2) * V[1] - 0.75 * sqrt2 * Y2_p1 * V[2];
    d[26] = d[38] = -Y3_p0 * Y3_p2 * V[0] + ((-3.0 / 14.0) * sqrt3 * Y4_p2 + (-6.0 / 7.0) * sqrt5 * Y2_p2 * R2) * V[1] + 1.5 * sqrt5 * Y2_p2 * V[2];
    d[27] = d[45] = -Y3_p0 * Y3_p3 * V[0] + (-9.0 / 14.0) * sqrt7 * Y4_p3 * V[1];
    d[23] = d[17] = -Y3_p0 * Y3_n1 * V[0] + ((3.0 / 14.0) * sqrt15 * Y4_n1 + (3.0 / 7.0) * sqrt2 * Y2_n1 * R2) * V[1] - 0.75 * sqrt2 * Y2_n1 * V[2];
    d[22] = d[10] = -Y3_p0 * Y3_n2 * V[0] + ((-3.0 / 14.0) * sqrt3 * Y4_n2 + (-6.0 / 7.0) * sqrt5 * Y2_n2 * R2) * V[1] + 1.5 * sqrt5 * Y2_n2 * V[2];
    d[21] = d[3] = -Y3_p0 * Y3_n3 * V[0] + (-9.0 / 14.0) * sqrt7 * Y4_n3 * V[1];
    d[32] = -Y3_p1 * Y3_p1 * V[0] + ((3.0 / 14.0) * Y4_p0 + (3.0 / 7.0) * sqrt5 * Y4_p2 + (9.0 / 7.0) * Y2_p0 * R2 + (6.0 / 7.0) * sqrt3 * Y2_p2 * R2 + 1.5 * R4) * V[1] + (-2.25 * Y2_p0 - 1.5 * sqrt3 * Y2_p2 - 3.75 * R2) * V[2] + 1.875 * V[3];
    d[33] = d[39] = -Y3_p1 * Y3_p2 * V[0] + ((6.0 / 7.0) * Y4_p1 + (3.0 / 14.0) * sqrt7 * Y4_p3 + (3.0 / 14.0) * sqrt30 * Y2_p1 * R2) * V[1] - 0.375 * sqrt30 * Y2_p1 * V[2];
    d[34] = d[46] = -Y3_p1 * Y3_p3 * V[0] + ((9.0 / 14.0) * sqrt3 * Y4_p2 + (-3.0 / 14.0) * sqrt21 * Y4_p4 + (-3.0 / 7.0) * sqrt5 * Y2_p2 * R2) * V[1] + 0.75 * sqrt5 * Y2_p2 * V[2];
    d[30] = d[18] = -Y3_p1 * Y3_n1 * V[0] + ((3.0 / 7.0) * sqrt5 * Y4_n2 + (6.0 / 7.0) * sqrt3 * Y2_n2 * R2) * V[1] - 1.5 * sqrt3 * Y2_n2 * V[2];
    d[29] = d[11] = -Y3_p1 * Y3_n2 * V[0] + ((3.0 / 14.0) * sqrt7 * Y4_n3 + (6.0 / 7.0) * Y4_n1 + (3.0 / 14.0) * sqrt30 * Y2_n1 * R2) * V[1] - 0.375 * sqrt30 * Y2_n1 * V[2];
    d[28] = d[4] = -Y3_p1 * Y3_n3 * V[0] + ((-3.0 / 14.0) * sqrt21 * Y4_n4 + (9.0 / 14.0) * sqrt3 * Y4_n2 + (-3.0 / 7.0) * sqrt5 * Y2_n2 * R2) * V[1] + 0.75 * sqrt5 * Y2_n2 * V[2];
    d[40] = -Y3_p2 * Y3_p2 * V[0] + (-1.5 * Y4_p0 + (3.0 / 14.0) * sqrt35 * Y4_p4 + 1.5 * R4) * V[1] - 3.75 * R2 * V[2] + 1.875 * V[3];
    d[41] = d[47] = -Y3_p2 * Y3_p3 * V[0] + ((-3.0 / 14.0) * sqrt15 * Y4_p1 + (15.0 / 14.0) * sqrt2 * Y2_p1 * R2) * V[1] - 1.875 * sqrt2 * Y2_p1 * V[2];
    d[37] = d[19] = -Y3_p2 * Y3_n1 * V[0] + ((3.0 / 14.0) * sqrt7 * Y4_n3 + (-6.0 / 7.0) * Y4_n1 + (-3.0 / 14.0) * sqrt30 * Y2_n1 * R2) * V[1] + 0.375 * sqrt30 * Y2_n1 * V[2];
    d[36] = d[12] = -Y3_p2 * Y3_n2 * V[0] + (3.0 / 14.0) * sqrt35 * Y4_n4 * V[1];
    d[35] = d[5] = -Y3_p2 * Y3_n3 * V[0] + ((-3.0 / 14.0) * sqrt15 * Y4_n1 + (15.0 / 14.0) * sqrt2 * Y2_n1 * R2) * V[1] - 1.875 * sqrt2 * Y2_n1 * V[2];
    d[48] = -Y3_p3 * Y3_p3 * V[0] + ((9.0 / 14.0) * Y4_p0 + (-15.0 / 7.0) * Y2_p0 * R2 + 1.5 * R4) * V[1] + (3.75 * Y2_p0 - 3.75 * R2) * V[2] + 1.875 * V[3];
    d[44] = d[20] = -Y3_p3 * Y3_n1 * V[0] + ((-3.0 / 14.0) * sqrt21 * Y4_n4 + (-9.0 / 14.0) * sqrt3 * Y4_n2 + (3.0 / 7.0) * sqrt5 * Y2_n2 * R2) * V[1] - 0.75 * sqrt5 * Y2_n2 * V[2];
    d[43] = d[13] = -Y3_p3 * Y3_n2 * V[0] + ((3.0 / 14.0) * sqrt15 * Y4_n1 + (-15.0 / 14.0) * sqrt2 * Y2_n1 * R2) * V[1] + 1.875 * sqrt2 * Y2_n1 * V[2];
    d[42] = d[6] = -Y3_p3 * Y3_n3 * V[0];
    d[16] = -Y3_n1 * Y3_n1 * V[0] + ((3.0 / 14.0) * Y4_p0 + (-3.0 / 7.0) * sqrt5 * Y4_p2 + (9.0 / 7.0) * Y2_p0 * R2 + (-6.0 / 7.0) * sqrt3 * Y2_p2 * R2 + 1.5 * R4) * V[1] + (-2.25 * Y2_p0 + 1.5 * sqrt3 * Y2_p2 - 3.75 * R2) * V[2] + 1.875 * V[3];
    d[15] = d[9] = -Y3_n1 * Y3_n2 * V[0] + ((6.0 / 7.0) * Y4_p1 + (-3.0 / 14.0) * sqrt7 * Y4_p3 + (3.0 / 14.0) * sqrt30 * Y2_p1 * R2) * V[1] - 0.375 * sqrt30 * Y2_p1 * V[2];
    d[14] = d[2] = -Y3_n1 * Y3_n3 * V[0] + ((9.0 / 14.0) * sqrt3 * Y4_p2 + (3.0 / 14.0) * sqrt21 * Y4_p4 + (-3.0 / 7.0) * sqrt5 * Y2_p2 * R2) * V[1] + 0.75 * sqrt5 * Y2_p2 * V[2];
    d[8] = -Y3_n2 * Y3_n2 * V[0] + (-1.5 * Y4_p0 + (-3.0 / 14.0) * sqrt35 * Y4_p4 + 1.5 * R4) * V[1] - 3.75 * R2 * V[2] + 1.875 * V[3];
    d[7] = d[1] = -Y3_n2 * Y3_n3 * V[0] + ((-3.0 / 14.0) * sqrt15 * Y4_p1 + (15.0 / 14.0) * sqrt2 * Y2_p1 * R2) * V[1] - 1.875 * sqrt2 * Y2_p1 * V[2];
    d[0] = -Y3_n3 * Y3_n3 * V[0] + ((9.0 / 14.0) * Y4_p0 + (-15.0 / 7.0) * Y2_p0 * R2 + 1.5 * R4) * V[1] + (3.75 * Y2_p0 - 3.75 * R2) * V[2] + 1.875 * V[3];

    return out;
}

}  // namespace ovlab