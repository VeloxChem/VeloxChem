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

#include "KineticEnergyABRecDD.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "MathConst.hpp"
#include "RealSolidHarmonicAB.hpp"

namespace kinab {

auto kinetic_d_d(
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
    const auto sqrt5 = std::sqrt(5.0);
    const auto sqrt10 = std::sqrt(10.0);
    const auto sqrt15 = std::sqrt(15.0);
    const auto sqrt30 = std::sqrt(30.0);
    const auto sqrt35 = std::sqrt(35.0);
    const auto sqrt70 = std::sqrt(70.0);

    // ---- Phase 2: V evaluation + primitive contraction ----
    //   V[0] ↔ α^2 · β^2 · p^{-4} · (s|T|s)
    //   V[1] ↔ α · β · p^{-3} · (s|T|s)
    //   V[2] ↔ p^{-2} · (s|T|s)
    //   V[3] ↔ α^2 · β^2 · p^{-4} · ξ · (s|S|s)
    //   V[4] ↔ α · β · p^{-3} · ξ · (s|S|s)
    //   V[5] ↔ p^{-2} · ξ · (s|S|s)
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

            const auto p     = alpha + beta;
            const auto pinv  = 1.0 / p;
            const auto pinv2 = pinv * pinv;
            const auto pinv3 = pinv2 * pinv;
            const auto pinv4 = pinv3 * pinv;
            const auto fpi   = pi * pinv;
            const auto ss    = fpi * std::sqrt(fpi) * std::exp(-alpha * beta * pinv * R2);
            const auto xi    = alpha * beta * pinv;
            const auto tt    = xi * (3.0 - 2.0 * xi * R2) * ss;
            const auto cab_ss = ca * cb * ss;
            const auto cab_tt = ca * cb * tt;

            V[0] += cab_tt * alpha2 * beta2 * pinv4;
            V[1] += cab_tt * alpha * beta * pinv3;
            V[2] += cab_tt * pinv2;
            V[3] += cab_ss * alpha2 * beta2 * pinv4 * xi;
            V[4] += cab_ss * alpha * beta * pinv3 * xi;
            V[5] += cab_ss * pinv2 * xi;
        }
    }

    // ---- Phase 3: fused M·V → 5 × 5 spherical block ----
    const auto Y2_n2 = harm::Y_ll_2_m_n2(AB_x, AB_y, AB_z);
    const auto Y2_n1 = harm::Y_ll_2_m_n1(AB_x, AB_y, AB_z);
    const auto Y2_p0 = harm::Y_ll_2_m_p0(AB_x, AB_y, AB_z);
    const auto Y2_p1 = harm::Y_ll_2_m_p1(AB_x, AB_y, AB_z);
    const auto Y2_p2 = harm::Y_ll_2_m_p2(AB_x, AB_y, AB_z);
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
    auto *d = buffer;
    d[12] = ((18.0 / 35.0) * Y4_p0 + (2.0 / 7.0) * Y2_p0 * R2 + 0.2 * R4) * V[0] + (-Y2_p0 - R2) * V[1] + 0.75 * V[2] + ((144.0 / 35.0) * Y4_p0 + (16.0 / 7.0) * Y2_p0 * R2 + 1.6 * R4) * V[3] + (-6.0 * Y2_p0 - 6.0 * R2) * V[4] + 3.0 * V[5];
    d[13] = d[17] = ((3.0 / 35.0) * sqrt30 * Y4_p1 + (1.0 / 7.0) * Y2_p1 * R2) * V[0] - 0.5 * Y2_p1 * V[1] + ((24.0 / 35.0) * sqrt30 * Y4_p1 + (8.0 / 7.0) * Y2_p1 * R2) * V[3] - 3.0 * Y2_p1 * V[4];
    d[14] = d[22] = ((3.0 / 35.0) * sqrt15 * Y4_p2 + (-2.0 / 7.0) * Y2_p2 * R2) * V[0] + Y2_p2 * V[1] + ((24.0 / 35.0) * sqrt15 * Y4_p2 + (-16.0 / 7.0) * Y2_p2 * R2) * V[3] + 6.0 * Y2_p2 * V[4];
    d[11] = d[7] = ((3.0 / 35.0) * sqrt30 * Y4_n1 + (1.0 / 7.0) * Y2_n1 * R2) * V[0] - 0.5 * Y2_n1 * V[1] + ((24.0 / 35.0) * sqrt30 * Y4_n1 + (8.0 / 7.0) * Y2_n1 * R2) * V[3] - 3.0 * Y2_n1 * V[4];
    d[10] = d[2] = ((3.0 / 35.0) * sqrt15 * Y4_n2 + (-2.0 / 7.0) * Y2_n2 * R2) * V[0] + Y2_n2 * V[1] + ((24.0 / 35.0) * sqrt15 * Y4_n2 + (-16.0 / 7.0) * Y2_n2 * R2) * V[3] + 6.0 * Y2_n2 * V[4];
    d[18] = ((-12.0 / 35.0) * Y4_p0 + (6.0 / 35.0) * sqrt5 * Y4_p2 + (1.0 / 7.0) * Y2_p0 * R2 + (1.0 / 7.0) * sqrt3 * Y2_p2 * R2 + 0.2 * R4) * V[0] + (-0.5 * Y2_p0 - 0.5 * sqrt3 * Y2_p2 - R2) * V[1] + 0.75 * V[2] + ((-96.0 / 35.0) * Y4_p0 + (48.0 / 35.0) * sqrt5 * Y4_p2 + (8.0 / 7.0) * Y2_p0 * R2 + (8.0 / 7.0) * sqrt3 * Y2_p2 * R2 + 1.6 * R4) * V[3] + (-3.0 * Y2_p0 - 3.0 * sqrt3 * Y2_p2 - 6.0 * R2) * V[4] + 3.0 * V[5];
    d[19] = d[23] = ((-3.0 / 70.0) * sqrt10 * Y4_p1 + (3.0 / 70.0) * sqrt70 * Y4_p3 + (1.0 / 7.0) * sqrt3 * Y2_p1 * R2) * V[0] - 0.5 * sqrt3 * Y2_p1 * V[1] + ((-12.0 / 35.0) * sqrt10 * Y4_p1 + (12.0 / 35.0) * sqrt70 * Y4_p3 + (8.0 / 7.0) * sqrt3 * Y2_p1 * R2) * V[3] - 3.0 * sqrt3 * Y2_p1 * V[4];
    d[16] = d[8] = ((6.0 / 35.0) * sqrt5 * Y4_n2 + (1.0 / 7.0) * sqrt3 * Y2_n2 * R2) * V[0] - 0.5 * sqrt3 * Y2_n2 * V[1] + ((48.0 / 35.0) * sqrt5 * Y4_n2 + (8.0 / 7.0) * sqrt3 * Y2_n2 * R2) * V[3] - 3.0 * sqrt3 * Y2_n2 * V[4];
    d[15] = d[3] = ((3.0 / 70.0) * sqrt70 * Y4_n3 + (-3.0 / 70.0) * sqrt10 * Y4_n1 + (1.0 / 7.0) * sqrt3 * Y2_n1 * R2) * V[0] - 0.5 * sqrt3 * Y2_n1 * V[1] + ((12.0 / 35.0) * sqrt70 * Y4_n3 + (-12.0 / 35.0) * sqrt10 * Y4_n1 + (8.0 / 7.0) * sqrt3 * Y2_n1 * R2) * V[3] - 3.0 * sqrt3 * Y2_n1 * V[4];
    d[24] = ((3.0 / 35.0) * Y4_p0 + (3.0 / 35.0) * sqrt35 * Y4_p4 + (-2.0 / 7.0) * Y2_p0 * R2 + 0.2 * R4) * V[0] + (Y2_p0 - R2) * V[1] + 0.75 * V[2] + ((24.0 / 35.0) * Y4_p0 + (24.0 / 35.0) * sqrt35 * Y4_p4 + (-16.0 / 7.0) * Y2_p0 * R2 + 1.6 * R4) * V[3] + (6.0 * Y2_p0 - 6.0 * R2) * V[4] + 3.0 * V[5];
    d[21] = d[9] = ((3.0 / 70.0) * sqrt70 * Y4_n3 + (3.0 / 70.0) * sqrt10 * Y4_n1 + (-1.0 / 7.0) * sqrt3 * Y2_n1 * R2) * V[0] + 0.5 * sqrt3 * Y2_n1 * V[1] + ((12.0 / 35.0) * sqrt70 * Y4_n3 + (12.0 / 35.0) * sqrt10 * Y4_n1 + (-8.0 / 7.0) * sqrt3 * Y2_n1 * R2) * V[3] + 3.0 * sqrt3 * Y2_n1 * V[4];
    d[20] = d[4] = (3.0 / 35.0) * sqrt35 * Y4_n4 * V[0] + (24.0 / 35.0) * sqrt35 * Y4_n4 * V[3];
    d[6] = ((-12.0 / 35.0) * Y4_p0 + (-6.0 / 35.0) * sqrt5 * Y4_p2 + (1.0 / 7.0) * Y2_p0 * R2 + (-1.0 / 7.0) * sqrt3 * Y2_p2 * R2 + 0.2 * R4) * V[0] + (-0.5 * Y2_p0 + 0.5 * sqrt3 * Y2_p2 - R2) * V[1] + 0.75 * V[2] + ((-96.0 / 35.0) * Y4_p0 + (-48.0 / 35.0) * sqrt5 * Y4_p2 + (8.0 / 7.0) * Y2_p0 * R2 + (-8.0 / 7.0) * sqrt3 * Y2_p2 * R2 + 1.6 * R4) * V[3] + (-3.0 * Y2_p0 + 3.0 * sqrt3 * Y2_p2 - 6.0 * R2) * V[4] + 3.0 * V[5];
    d[5] = d[1] = ((-3.0 / 70.0) * sqrt10 * Y4_p1 + (-3.0 / 70.0) * sqrt70 * Y4_p3 + (1.0 / 7.0) * sqrt3 * Y2_p1 * R2) * V[0] - 0.5 * sqrt3 * Y2_p1 * V[1] + ((-12.0 / 35.0) * sqrt10 * Y4_p1 + (-12.0 / 35.0) * sqrt70 * Y4_p3 + (8.0 / 7.0) * sqrt3 * Y2_p1 * R2) * V[3] - 3.0 * sqrt3 * Y2_p1 * V[4];
    d[0] = ((3.0 / 35.0) * Y4_p0 + (-3.0 / 35.0) * sqrt35 * Y4_p4 + (-2.0 / 7.0) * Y2_p0 * R2 + 0.2 * R4) * V[0] + (Y2_p0 - R2) * V[1] + 0.75 * V[2] + ((24.0 / 35.0) * Y4_p0 + (-24.0 / 35.0) * sqrt35 * Y4_p4 + (-16.0 / 7.0) * Y2_p0 * R2 + 1.6 * R4) * V[3] + (6.0 * Y2_p0 - 6.0 * R2) * V[4] + 3.0 * V[5];
}

}  // namespace kinab