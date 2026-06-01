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

#include "KineticEnergyABRecDP.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_d_p(
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
    const auto sqrt30 = std::sqrt(30.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α · β^2 · p^{-3} · (s|T|s)
    //   V[1] ↔ β · p^{-2} · (s|T|s)
    //   V[2] ↔ α · β^2 · p^{-3} · ξ · (s|S|s)
    //   V[3] ↔ β · p^{-2} · ξ · (s|S|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha * beta2 * pinv3;
            V[1] += cab_tt * beta * pinv2;
            V[2] += cab_ss * alpha * beta2 * pinv3 * xi;
            V[3] += cab_ss * beta * pinv2 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 5 × 3 spherical block ----
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
    auto *d = buffer;
    d[7] = (0.6 * Y3_p0 + 0.4 * Y1_p0 * R2) * V[0] - Y1_p0 * V[1] + (3.6 * Y3_p0 + 2.4 * Y1_p0 * R2) * V[2] - 4.0 * Y1_p0 * V[3];
    d[8] = (0.2 * sqrt6 * Y3_p1 - 0.2 * Y1_p1 * R2) * V[0] + 0.5 * Y1_p1 * V[1] + (1.2 * sqrt6 * Y3_p1 - 1.2 * Y1_p1 * R2) * V[2] + 2.0 * Y1_p1 * V[3];
    d[6] = (0.2 * sqrt6 * Y3_n1 - 0.2 * Y1_n1 * R2) * V[0] + 0.5 * Y1_n1 * V[1] + (1.2 * sqrt6 * Y3_n1 - 1.2 * Y1_n1 * R2) * V[2] + 2.0 * Y1_n1 * V[3];
    d[10] = (0.4 * sqrt2 * Y3_p1 + 0.2 * sqrt3 * Y1_p1 * R2) * V[0] - 0.5 * sqrt3 * Y1_p1 * V[1] + (2.4 * sqrt2 * Y3_p1 + 1.2 * sqrt3 * Y1_p1 * R2) * V[2] - 2.0 * sqrt3 * Y1_p1 * V[3];
    d[11] = (-0.2 * sqrt3 * Y3_p0 + 0.2 * sqrt5 * Y3_p2 + 0.2 * sqrt3 * Y1_p0 * R2) * V[0] - 0.5 * sqrt3 * Y1_p0 * V[1] + (-1.2 * sqrt3 * Y3_p0 + 1.2 * sqrt5 * Y3_p2 + 1.2 * sqrt3 * Y1_p0 * R2) * V[2] - 2.0 * sqrt3 * Y1_p0 * V[3];
    d[9] = 0.2 * sqrt5 * Y3_n2 * V[0] + 1.2 * sqrt5 * Y3_n2 * V[2];
    d[13] = 0.2 * sqrt5 * Y3_p2 * V[0] + 1.2 * sqrt5 * Y3_p2 * V[2];
    d[14] = (-0.1 * sqrt2 * Y3_p1 + 0.1 * sqrt30 * Y3_p3 + 0.2 * sqrt3 * Y1_p1 * R2) * V[0] - 0.5 * sqrt3 * Y1_p1 * V[1] + (-0.6 * sqrt2 * Y3_p1 + 0.6 * sqrt30 * Y3_p3 + 1.2 * sqrt3 * Y1_p1 * R2) * V[2] - 2.0 * sqrt3 * Y1_p1 * V[3];
    d[12] = (0.1 * sqrt30 * Y3_n3 + 0.1 * sqrt2 * Y3_n1 - 0.2 * sqrt3 * Y1_n1 * R2) * V[0] + 0.5 * sqrt3 * Y1_n1 * V[1] + (0.6 * sqrt30 * Y3_n3 + 0.6 * sqrt2 * Y3_n1 - 1.2 * sqrt3 * Y1_n1 * R2) * V[2] + 2.0 * sqrt3 * Y1_n1 * V[3];
    d[4] = (0.4 * sqrt2 * Y3_n1 + 0.2 * sqrt3 * Y1_n1 * R2) * V[0] - 0.5 * sqrt3 * Y1_n1 * V[1] + (2.4 * sqrt2 * Y3_n1 + 1.2 * sqrt3 * Y1_n1 * R2) * V[2] - 2.0 * sqrt3 * Y1_n1 * V[3];
    d[5] = 0.2 * sqrt5 * Y3_n2 * V[0] + 1.2 * sqrt5 * Y3_n2 * V[2];
    d[3] = (-0.2 * sqrt3 * Y3_p0 - 0.2 * sqrt5 * Y3_p2 + 0.2 * sqrt3 * Y1_p0 * R2) * V[0] - 0.5 * sqrt3 * Y1_p0 * V[1] + (-1.2 * sqrt3 * Y3_p0 - 1.2 * sqrt5 * Y3_p2 + 1.2 * sqrt3 * Y1_p0 * R2) * V[2] - 2.0 * sqrt3 * Y1_p0 * V[3];
    d[1] = 0.2 * sqrt5 * Y3_n2 * V[0] + 1.2 * sqrt5 * Y3_n2 * V[2];
    d[2] = (0.1 * sqrt30 * Y3_n3 - 0.1 * sqrt2 * Y3_n1 + 0.2 * sqrt3 * Y1_n1 * R2) * V[0] - 0.5 * sqrt3 * Y1_n1 * V[1] + (0.6 * sqrt30 * Y3_n3 - 0.6 * sqrt2 * Y3_n1 + 1.2 * sqrt3 * Y1_n1 * R2) * V[2] - 2.0 * sqrt3 * Y1_n1 * V[3];
    d[0] = (-0.1 * sqrt2 * Y3_p1 - 0.1 * sqrt30 * Y3_p3 + 0.2 * sqrt3 * Y1_p1 * R2) * V[0] - 0.5 * sqrt3 * Y1_p1 * V[1] + (-0.6 * sqrt2 * Y3_p1 - 0.6 * sqrt30 * Y3_p3 + 1.2 * sqrt3 * Y1_p1 * R2) * V[2] - 2.0 * sqrt3 * Y1_p1 * V[3];
}

}  // namespace kinab