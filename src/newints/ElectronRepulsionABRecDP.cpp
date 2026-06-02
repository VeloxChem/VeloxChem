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

#include "ElectronRepulsionABRecDP.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "BoysFunction.hpp"
#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace eri2cab {

auto eri2c_d_p(
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

    const auto sqrt3 = std::sqrt(3.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ β · p^{-2} · (s|s)^(2)
    //   V[1] ↔ α · β^2 · p^{-3} · (s|s)^(3)
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

            const auto p    = alpha + beta;
            const auto pinv = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto xi   = alpha * beta * pinv;

            const auto T = xi * R2;
            std::array<double, 4> F;
            newints::boys_function<4>(T, F);

            const auto pref     = two_pi52 / (alpha * beta * std::sqrt(p));
            const auto cab_pref = ca * cb * pref;

            V[0] += cab_pref * beta * pinv2 * F[2];
            V[1] += cab_pref * alpha * beta2 * pinv3 * F[3];
        }
    }

    // ---- Phase 3: fused M·V → 5 × 3 spherical block ----
    const auto Y1_n1 = harm::Y_ll_1_m_n1(AB_x, AB_y, AB_z);
    const auto Y1_p0 = harm::Y_ll_1_m_p0(AB_x, AB_y, AB_z);
    const auto Y1_p1 = harm::Y_ll_1_m_p1(AB_x, AB_y, AB_z);
    const auto Y2_n2 = harm::Y_ll_2_m_n2(AB_x, AB_y, AB_z);
    const auto Y2_n1 = harm::Y_ll_2_m_n1(AB_x, AB_y, AB_z);
    const auto Y2_p0 = harm::Y_ll_2_m_p0(AB_x, AB_y, AB_z);
    const auto Y2_p1 = harm::Y_ll_2_m_p1(AB_x, AB_y, AB_z);
    const auto Y2_p2 = harm::Y_ll_2_m_p2(AB_x, AB_y, AB_z);
    auto *d = buffer;
    d[7] = Y2_p0 * Y1_p0 * V[1] - Y1_p0 * V[0];
    d[8] = Y2_p0 * Y1_p1 * V[1] + 0.5 * Y1_p1 * V[0];
    d[6] = Y2_p0 * Y1_n1 * V[1] + 0.5 * Y1_n1 * V[0];
    d[10] = Y2_p1 * Y1_p0 * V[1] - 0.5 * sqrt3 * Y1_p1 * V[0];
    d[11] = Y2_p1 * Y1_p1 * V[1] - 0.5 * sqrt3 * Y1_p0 * V[0];
    d[9] = Y2_p1 * Y1_n1 * V[1];
    d[13] = Y2_p2 * Y1_p0 * V[1];
    d[14] = Y2_p2 * Y1_p1 * V[1] - 0.5 * sqrt3 * Y1_p1 * V[0];
    d[12] = Y2_p2 * Y1_n1 * V[1] + 0.5 * sqrt3 * Y1_n1 * V[0];
    d[4] = Y2_n1 * Y1_p0 * V[1] - 0.5 * sqrt3 * Y1_n1 * V[0];
    d[5] = Y2_n1 * Y1_p1 * V[1];
    d[3] = Y2_n1 * Y1_n1 * V[1] - 0.5 * sqrt3 * Y1_p0 * V[0];
    d[1] = Y2_n2 * Y1_p0 * V[1];
    d[2] = Y2_n2 * Y1_p1 * V[1] - 0.5 * sqrt3 * Y1_n1 * V[0];
    d[0] = Y2_n2 * Y1_n1 * V[1] - 0.5 * sqrt3 * Y1_p1 * V[0];
}

}  // namespace eri2cab