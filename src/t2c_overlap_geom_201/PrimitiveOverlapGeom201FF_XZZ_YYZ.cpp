#include "PrimitiveOverlapGeom201FF_XZZ_YYZ.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapGeom201FF_XZZ_YYZ(TDoubleArray&       buffer_xx_x,
                                      TDoubleArray&       buffer_xx_y,
                                      TDoubleArray&       buffer_xx_z,
                                      TDoubleArray&       buffer_xy_x,
                                      TDoubleArray&       buffer_xy_y,
                                      TDoubleArray&       buffer_xy_z,
                                      TDoubleArray&       buffer_xz_x,
                                      TDoubleArray&       buffer_xz_y,
                                      TDoubleArray&       buffer_xz_z,
                                      TDoubleArray&       buffer_yy_x,
                                      TDoubleArray&       buffer_yy_y,
                                      TDoubleArray&       buffer_yy_z,
                                      TDoubleArray&       buffer_yz_x,
                                      TDoubleArray&       buffer_yz_y,
                                      TDoubleArray&       buffer_yz_z,
                                      TDoubleArray&       buffer_zz_x,
                                      TDoubleArray&       buffer_zz_y,
                                      TDoubleArray&       buffer_zz_z,
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

    auto fints_xx_x = buffer_xx_x.data();

    auto fints_xx_y = buffer_xx_y.data();

    auto fints_xx_z = buffer_xx_z.data();

    auto fints_xy_x = buffer_xy_x.data();

    auto fints_xy_y = buffer_xy_y.data();

    auto fints_xy_z = buffer_xy_z.data();

    auto fints_xz_x = buffer_xz_x.data();

    auto fints_xz_y = buffer_xz_y.data();

    auto fints_xz_z = buffer_xz_z.data();

    auto fints_yy_x = buffer_yy_x.data();

    auto fints_yy_y = buffer_yy_y.data();

    auto fints_yy_z = buffer_yy_z.data();

    auto fints_yz_x = buffer_yz_x.data();

    auto fints_yz_y = buffer_yz_y.data();

    auto fints_yz_z = buffer_yz_z.data();

    auto fints_zz_x = buffer_zz_x.data();

    auto fints_zz_y = buffer_zz_y.data();

    auto fints_zz_z = buffer_zz_z.data();

#pragma omp simd aligned(fints_xx_x,     \
                             fints_xx_y, \
                             fints_xx_z, \
                             fints_xy_x, \
                             fints_xy_y, \
                             fints_xy_z, \
                             fints_xz_x, \
                             fints_xz_y, \
                             fints_xz_z, \
                             fints_yy_x, \
                             fints_yy_y, \
                             fints_yy_z, \
                             fints_yz_x, \
                             fints_yz_y, \
                             fints_yz_z, \
                             fints_zz_x, \
                             fints_zz_y, \
                             fints_zz_z, \
                             ket_fe,     \
                             ket_fn,     \
                             ket_rx,     \
                             ket_ry,     \
                             ket_rz : 64)
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

        const auto tke_0 = ket_fe[i];

        const auto tbe_0 = bra_exp;

        fints_xx_x[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 - (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 +
                   fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
                   (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] +=
            fss *
            ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] +=
            fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_xx_x[i] +=
            fss * (-3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 +
                   (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] +=
            fss *
            (fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0);

        fints_xx_x[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
             6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
             6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] +=
            fss * (6.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_xx_y[i] +=
            fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 +
                   12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0);

        fints_xx_y[i] +=
            fss * (-3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0);

        fints_xx_y[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                                6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_y[i] += fss * (4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * tbe_0 + 3.0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 -
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0);

        fints_xx_z[i] +=
            fss * (-fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] +=
            fss * (-2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 3.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 -
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 6.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0);

        fints_xx_z[i] +=
            fss * (-3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
                   fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xx_z[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 + (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xx_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (-4.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 -
             6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 - 6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 -
             12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (-12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] +=
            fss *
            (fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tke_0);

        fints_xy_x[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 -
                                2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xy_x[i] += fss * (-4.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (fe_0 * fe_0 * rpa_z * tbe_0 + fe_0 * fe_0 * rpa_z * tbe_0 - fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_xy_y[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 + fe_0 * fe_0 * rpb_z * tbe_0 +
                                2.0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 + 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 +
                                2.0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0);

        fints_xy_y[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 - (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 -
                                fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0);

        fints_xy_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_xy_y[i] +=
            fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 +
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0);

        fints_xy_y[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 + 4.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 -
                                fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0);

        fints_xy_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0);

        fints_xy_y[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xy_y[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xy_y[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * tbe_0 + fe_0 * fe_0 * rpb_y * tbe_0 + fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 +
                                2.0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 - fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 -
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_xy_z[i] +=
            fss * (-fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 +
                   fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_xy_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 + fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 -
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (-fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 2.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 -
                                fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tke_0);

        fints_xy_z[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 -
                   fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xy_z[i] +=
            fss *
            ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 +
             fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xy_z[i] +=
            fss *
            (-2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * (-2.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_z[i] += fss * 8.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0;

        fints_xz_x[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] +=
            fss *
            ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * rpb_x * tke_0 -
                                fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (-fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0);

        fints_xz_x[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] +=
            fss *
            (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (2.0 * fe_0 * rpa_z * rpb_z * rpb_x * tke_0 + 2.0 * fe_0 * rpb_y * rpb_y * rpb_x * tke_0 -
                                fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xz_x[i] +=
            fss * (-fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xz_x[i] +=
            fss *
            (-2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tke_0 - 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (-4.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (-2.0 * fe_0 * rpb_y + 2.0 * fe_0 * fe_0 * rpb_y * tbe_0 + fe_0 * fe_0 * rpb_y * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpb_y * tbe_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0);

        fints_xz_y[i] += fss * (2.0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 +
                                4.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0);

        fints_xz_y[i] +=
            fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0);

        fints_xz_y[i] +=
            fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 - 4.0 * rpa_z * rpb_z * rpb_y + fe_0 * fe_0 * rpb_y * tke_0);

        fints_xz_y[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpb_y * tke_0 + 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 + 4.0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 -
                   fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 - fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (4.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (-fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 -
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 +
                   fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] +=
            fss *
            (fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 -
             8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0);

        fints_xz_y[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * rpa_z * rpb_z * rpb_y * tke_0 + 4.0 * fe_0 * rpa_z * rpb_z * rpb_y * tke_0 +
                                2.0 * fe_0 * rpb_y * rpb_y * rpb_y * tke_0 - fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] +=
            fss *
            (-2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 - 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tke_0 - 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (-4.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (-fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * tbe_0 + fe_0 * fe_0 * rpa_z * tbe_0 +
                                fe_0 * fe_0 * rpa_z * tbe_0 + fe_0 * rpa_z * rpa_z * rpa_z * tbe_0);

        fints_xz_z[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 -
                                fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 -
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0);

        fints_xz_z[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 +
                                fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] +=
            fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * rpa_z * rpb_y * rpb_y +
                                fe_0 * fe_0 * rpa_z * tke_0 + 2.0 * fe_0 * fe_0 * rpb_z * tke_0 + fe_0 * rpa_z * rpb_y * rpb_y * tbe_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 -
                                fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 -
                                fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 - fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 + 2.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 -
                                fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tke_0);

        fints_xz_z[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 + 4.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tke_0 -
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xz_z[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 + (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (-2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpb_y * rpb_y * tke_0 +
                   2.0 * fe_0 * rpa_z * rpb_z * rpb_z * tke_0 + 4.0 * fe_0 * rpb_z * rpb_y * rpb_y * tke_0);

        fints_xz_z[i] +=
            fss * (-fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_z[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_z[i] +=
            fss *
            (-2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 -
                                8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tke_0);

        fints_xz_z[i] += fss * (-2.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xz_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_z[i] += fss * 8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0;

        fints_yy_x[i] +=
            fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 +
                   (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (-fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0);

        fints_yy_x[i] +=
            fss *
            (-fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 +
                   2.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0);

        fints_yy_y[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_yy_y[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (2.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 + 4.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 -
                                fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yy_y[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
             3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_yy_y[i] +=
            fss *
            (-2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * (-2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_y[i] += fss * 8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0;

        fints_yy_z[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * tbe_0 +
                                fe_0 * rpa_z * rpa_z * rpa_x * tbe_0);

        fints_yy_z[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 -
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_yy_z[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_yy_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 + fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 -
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 2.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0);

        fints_yy_z[i] += fss * (-fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 -
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_yy_z[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                2.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 - fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0);

        fints_yy_z[i] +=
            fss *
            (-2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * (-2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_z[i] += fss * 8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0;

        fints_yz_x[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 +
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] +=
            fss *
            (fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] +=
            fss *
            (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tke_0);

        fints_yz_x[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                   fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_x[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_x * tbe_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0);

        fints_yz_y[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 + 4.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 +
                                4.0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 -
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_yz_y[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] +=
            fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + 8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yz_y[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (-4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_y[i] += fss * (12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 + 4.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 -
                                fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 -
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_yz_z[i] +=
            fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] +=
            fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 + 4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0);

        fints_yz_z[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_yz_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_yz_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] +=
            fss *
            (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 - 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 -
             8.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_yz_z[i] += fss * (4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (-(5.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 - (5.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] +=
            fss *
            ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
             (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 +
             fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * rpb_z * tke_0 -
                                (5.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_zz_x[i] +=
            fss * (-5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 -
                   5.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * tke_0);

        fints_zz_x[i] +=
            fss * (2.0 * fe_0 * rpb_z * rpb_y * rpb_y * tke_0 - 5.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
                   5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 10.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
                   10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (-10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                10.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_x[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tke_0 - 10.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                   20.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_zz_x[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] +=
            fss * (10.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 + 10.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_zz_y[i] +=
            fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_zz_y[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 4.0 * rpa_x * rpb_z * rpb_y +
                                10.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 + 20.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0);

        fints_zz_y[i] +=
            fss * (-5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                   5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 - 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0);

        fints_zz_y[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (-8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * tke_0);

        fints_zz_y[i] +=
            fss * (4.0 * fe_0 * rpa_x * rpb_z * rpb_y * tke_0 - 5.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 - 10.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                   20.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (-10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_y[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   4.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tke_0 - 10.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   20.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_zz_y[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (-fe_0 * rpa_x + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * tbe_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 -
                                5.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 - (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0);

        fints_zz_z[i] +=
            fss * (-fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 -
                   3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] +=
            fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * rpa_x * rpb_y * rpb_y +
                                fe_0 * fe_0 * rpa_x * tke_0);

        fints_zz_z[i] += fss * (5.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 - (5.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                                10.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 - 5.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tke_0 -
                                10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0);

        fints_zz_z[i] +=
            fss *
            (-10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 - 10.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
             3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_zz_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                4.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_zz_z[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpb_y * rpb_y * tke_0 +
                   2.0 * fe_0 * rpa_x * rpb_z * rpb_z * tke_0 - 5.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
                   5.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (-10.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
                                10.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 -
                                20.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                20.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tke_0);

        fints_zz_z[i] += fss * (-10.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                                20.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                                6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_z[i] += fss * (12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);
    }
}

}  // namespace ovlrec
