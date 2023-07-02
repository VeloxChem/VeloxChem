#include "PrimitiveKineticEnergyFD_YZZ_T.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace kinrec {  // kinrec namespace

auto
compPrimitiveKineticEnergyFD_YZZ_T(TDoubleArray&       buffer_xx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xx[i] +=
            fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * fz_0 - fbe_0 * rpa_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * fz_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0 + 4.0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_xx[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * fz_0 +
                              10.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0);

        fints_xx[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x +
                              (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y * rpb_x * rpb_x);

        fints_xy[i] +=
            fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_y * rpb_x * fz_0 + 4.0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 +
                   4.0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xy[i] += fss * 10.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0;

        fints_xy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x +
                              (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_y * rpb_y * rpb_x);

        fints_xz[i] += fss * (-fbe_0 * rpa_y * rpb_z * rpb_x * fz_0 + 8.0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 +
                              4.0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 + 10.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0);

        fints_xz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x + rpa_z * rpa_z * rpa_y * rpb_z * rpb_x);

        fints_yy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * fz_0 - fbe_0 * fe_0 * rpb_y * fz_0 - fbe_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                              (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * fz_0 + 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0);

        fints_yy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 + 4.0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 +
                              (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * fz_0 + 3.0 * fe_0 * fe_0 * rpb_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * fz_0);

        fints_yy[i] += fss * 10.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0;

        fints_yy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y + fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y +
                              (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y;

        fints_yz[i] +=
            fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_y * fz_0 + 8.0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + 4.0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yz[i] +=
            fss * (3.0 * fe_0 * fe_0 * rpa_z * fz_0 + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * fz_0 + 10.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y +
                              (1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z);

        fints_yz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y;

        fints_zz[i] +=
            fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * fz_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 + 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0);

        fints_zz[i] += fss * (4.0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * fz_0 -
                              fke_0 * rpa_z * rpa_z * rpa_y * fz_0 + 10.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0);

        fints_zz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y +
                              (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y * rpb_z * rpb_z);
    }
}

}  // namespace kinrec
