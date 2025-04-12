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

#include "PrimitiveOverlapFD_XYZ_T.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapFD_XYZ_T(TDoubleArray&       buffer_xx,
                             TDoubleArray&       buffer_xy,
                             TDoubleArray&       buffer_xz,
                             TDoubleArray&       buffer_yy,
                             TDoubleArray&       buffer_yz,
                             TDoubleArray&       buffer_zz,
                             const double        bra_exp,
                             const double        bra_norm,
                             const TPoint3D&     bra_coord,
                             const TDoubleArray& ket_exps,
                             const TDoubleArray& ket_norms,
                             const TDoubleArray& ket_coords_x,
                             const TDoubleArray& ket_coords_y,
                             const TDoubleArray& ket_coords_z,
                             const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

#pragma omp simd aligned(fints_xx, fints_xy, fints_xz, fints_yy, fints_yz, fints_zz, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x + fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_y * rpa_x * rpb_x * rpb_x);

        fints_xy[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x +
                              (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_y * rpa_x * rpb_y * rpb_x);

        fints_xz[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x +
                              (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_y * rpa_x * rpb_z * rpb_x);

        fints_yy[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x + fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_y * rpa_x * rpb_y * rpb_y);

        fints_yz[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y +
                              (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_y * rpa_x * rpb_z * rpb_y);

        fints_zz[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x + fe_0 * rpa_y * rpa_x * rpb_z + rpa_z * rpa_y * rpa_x * rpb_z * rpb_z);
    }
}

}  // namespace ovlrec
