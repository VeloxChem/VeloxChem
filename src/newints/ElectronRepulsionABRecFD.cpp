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

#include "ElectronRepulsionABRecFD.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "BoysFunction.hpp"
#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace eri2cab {

auto eri2c_f_d(
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
    const auto sqrt30 = std::sqrt(30.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ β · p^{-3} · (s|s)^(3)
    //   V[1] ↔ α · β^2 · p^{-4} · (s|s)^(4)
    //   V[2] ↔ α^2 · β^3 · p^{-5} · (s|s)^(5)
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

            const auto p    = alpha + beta;
            const auto pinv = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto pinv5 = pinv4 * pinv;
            const auto xi   = alpha * beta * pinv;

            const auto T = xi * R2;
            std::array<double, 6> F;
            newints::boys_function<6>(T, F);

            const auto pref     = two_pi52 / (alpha * beta * std::sqrt(p));
            const auto cab_pref = ca * cb * pref;

            V[0] += cab_pref * beta * pinv3 * F[3];
            V[1] += cab_pref * alpha * beta2 * pinv4 * F[4];
            V[2] += cab_pref * alpha2 * beta3 * pinv5 * F[5];
        }
    }

    // ---- Phase 3: fused M·V → 7 × 5 spherical block ----
    const auto Y1_n1 = harm::Y_ll_1_m_n1(AB_x, AB_y, AB_z);
    const auto Y1_p0 = harm::Y_ll_1_m_p0(AB_x, AB_y, AB_z);
    const auto Y1_p1 = harm::Y_ll_1_m_p1(AB_x, AB_y, AB_z);
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
    auto *d = buffer;
    d[17] = -Y3_p0 * Y2_p0 * V[2] - 2.25 * Y1_p0 * V[0] + (1.2 * Y3_p0 + 1.8 * Y1_p0 * R2) * V[1];
    d[18] = -Y3_p0 * Y2_p1 * V[2] + 0.75 * sqrt3 * Y1_p1 * V[0] + (0.3 * sqrt2 * Y3_p1 - 0.6 * sqrt3 * Y1_p1 * R2) * V[1];
    d[19] = -Y3_p0 * Y2_p2 * V[2] - 0.6 * sqrt5 * Y3_p2 * V[1];
    d[16] = -Y3_p0 * Y2_n1 * V[2] + 0.75 * sqrt3 * Y1_n1 * V[0] + (0.3 * sqrt2 * Y3_n1 - 0.6 * sqrt3 * Y1_n1 * R2) * V[1];
    d[15] = -Y3_p0 * Y2_n2 * V[2] - 0.6 * sqrt5 * Y3_n2 * V[1];
    d[22] = -Y3_p1 * Y2_p0 * V[2] - 0.75 * sqrt6 * Y1_p1 * V[0] + (0.9 * Y3_p1 + 0.6 * sqrt6 * Y1_p1 * R2) * V[1];
    d[23] = -Y3_p1 * Y2_p1 * V[2] - 1.5 * sqrt2 * Y1_p0 * V[0] + (0.3 * sqrt2 * Y3_p0 + 0.15 * sqrt30 * Y3_p2 + 1.2 * sqrt2 * Y1_p0 * R2) * V[1];
    d[24] = -Y3_p1 * Y2_p2 * V[2] + 0.375 * sqrt2 * Y1_p1 * V[0] + (0.6 * sqrt3 * Y3_p1 - 0.3 * sqrt5 * Y3_p3 - 0.3 * sqrt2 * Y1_p1 * R2) * V[1];
    d[21] = -Y3_p1 * Y2_n1 * V[2] + 0.15 * sqrt30 * Y3_n2 * V[1];
    d[20] = -Y3_p1 * Y2_n2 * V[2] + 0.375 * sqrt2 * Y1_n1 * V[0] + (-0.3 * sqrt5 * Y3_n3 + 0.6 * sqrt3 * Y3_n1 - 0.3 * sqrt2 * Y1_n1 * R2) * V[1];
    d[27] = -Y3_p2 * Y2_p0 * V[2];
    d[28] = -Y3_p2 * Y2_p1 * V[2] - 0.75 * sqrt5 * Y1_p1 * V[0] + (0.15 * sqrt30 * Y3_p1 + 0.75 * sqrt2 * Y3_p3 + 0.6 * sqrt5 * Y1_p1 * R2) * V[1];
    d[29] = -Y3_p2 * Y2_p2 * V[2] - 0.75 * sqrt5 * Y1_p0 * V[0] + (-0.6 * sqrt5 * Y3_p0 + 0.6 * sqrt5 * Y1_p0 * R2) * V[1];
    d[26] = -Y3_p2 * Y2_n1 * V[2] + 0.75 * sqrt5 * Y1_n1 * V[0] + (0.75 * sqrt2 * Y3_n3 - 0.15 * sqrt30 * Y3_n1 - 0.6 * sqrt5 * Y1_n1 * R2) * V[1];
    d[25] = -Y3_p2 * Y2_n2 * V[2];
    d[32] = -Y3_p3 * Y2_p0 * V[2] - 1.5 * Y3_p3 * V[1];
    d[33] = -Y3_p3 * Y2_p1 * V[2] + 0.75 * sqrt2 * Y3_p2 * V[1];
    d[34] = -Y3_p3 * Y2_p2 * V[2] - 0.375 * sqrt30 * Y1_p1 * V[0] + (-0.3 * sqrt5 * Y3_p1 + 0.3 * sqrt30 * Y1_p1 * R2) * V[1];
    d[31] = -Y3_p3 * Y2_n1 * V[2] - 0.75 * sqrt2 * Y3_n2 * V[1];
    d[30] = -Y3_p3 * Y2_n2 * V[2] + 0.375 * sqrt30 * Y1_n1 * V[0] + (0.3 * sqrt5 * Y3_n1 - 0.3 * sqrt30 * Y1_n1 * R2) * V[1];
    d[12] = -Y3_n1 * Y2_p0 * V[2] - 0.75 * sqrt6 * Y1_n1 * V[0] + (0.9 * Y3_n1 + 0.6 * sqrt6 * Y1_n1 * R2) * V[1];
    d[13] = -Y3_n1 * Y2_p1 * V[2] + 0.15 * sqrt30 * Y3_n2 * V[1];
    d[14] = -Y3_n1 * Y2_p2 * V[2] - 0.375 * sqrt2 * Y1_n1 * V[0] + (-0.3 * sqrt5 * Y3_n3 - 0.6 * sqrt3 * Y3_n1 + 0.3 * sqrt2 * Y1_n1 * R2) * V[1];
    d[11] = -Y3_n1 * Y2_n1 * V[2] - 1.5 * sqrt2 * Y1_p0 * V[0] + (0.3 * sqrt2 * Y3_p0 - 0.15 * sqrt30 * Y3_p2 + 1.2 * sqrt2 * Y1_p0 * R2) * V[1];
    d[10] = -Y3_n1 * Y2_n2 * V[2] + 0.375 * sqrt2 * Y1_p1 * V[0] + (0.6 * sqrt3 * Y3_p1 + 0.3 * sqrt5 * Y3_p3 - 0.3 * sqrt2 * Y1_p1 * R2) * V[1];
    d[7] = -Y3_n2 * Y2_p0 * V[2];
    d[8] = -Y3_n2 * Y2_p1 * V[2] - 0.75 * sqrt5 * Y1_n1 * V[0] + (0.75 * sqrt2 * Y3_n3 + 0.15 * sqrt30 * Y3_n1 + 0.6 * sqrt5 * Y1_n1 * R2) * V[1];
    d[9] = -Y3_n2 * Y2_p2 * V[2];
    d[6] = -Y3_n2 * Y2_n1 * V[2] - 0.75 * sqrt5 * Y1_p1 * V[0] + (0.15 * sqrt30 * Y3_p1 - 0.75 * sqrt2 * Y3_p3 + 0.6 * sqrt5 * Y1_p1 * R2) * V[1];
    d[5] = -Y3_n2 * Y2_n2 * V[2] - 0.75 * sqrt5 * Y1_p0 * V[0] + (-0.6 * sqrt5 * Y3_p0 + 0.6 * sqrt5 * Y1_p0 * R2) * V[1];
    d[2] = -Y3_n3 * Y2_p0 * V[2] - 1.5 * Y3_n3 * V[1];
    d[3] = -Y3_n3 * Y2_p1 * V[2] + 0.75 * sqrt2 * Y3_n2 * V[1];
    d[4] = -Y3_n3 * Y2_p2 * V[2] - 0.375 * sqrt30 * Y1_n1 * V[0] + (-0.3 * sqrt5 * Y3_n1 + 0.3 * sqrt30 * Y1_n1 * R2) * V[1];
    d[1] = -Y3_n3 * Y2_n1 * V[2] + 0.75 * sqrt2 * Y3_p2 * V[1];
    d[0] = -Y3_n3 * Y2_n2 * V[2] - 0.375 * sqrt30 * Y1_p1 * V[0] + (-0.3 * sqrt5 * Y3_p1 + 0.3 * sqrt30 * Y1_p1 * R2) * V[1];
}

}  // namespace eri2cab