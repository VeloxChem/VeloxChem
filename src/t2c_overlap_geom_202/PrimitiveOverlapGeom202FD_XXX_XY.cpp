#include "PrimitiveOverlapGeom202FD_XXX_XY.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapGeom202FD_XXX_XY(TDoubleArray&       buffer_xx_xx,
                                     TDoubleArray&       buffer_xx_xy,
                                     TDoubleArray&       buffer_xx_xz,
                                     TDoubleArray&       buffer_xx_yy,
                                     TDoubleArray&       buffer_xx_yz,
                                     TDoubleArray&       buffer_xx_zz,
                                     TDoubleArray&       buffer_xy_xx,
                                     TDoubleArray&       buffer_xy_xy,
                                     TDoubleArray&       buffer_xy_xz,
                                     TDoubleArray&       buffer_xy_yy,
                                     TDoubleArray&       buffer_xy_yz,
                                     TDoubleArray&       buffer_xy_zz,
                                     TDoubleArray&       buffer_xz_xx,
                                     TDoubleArray&       buffer_xz_xy,
                                     TDoubleArray&       buffer_xz_xz,
                                     TDoubleArray&       buffer_xz_yy,
                                     TDoubleArray&       buffer_xz_yz,
                                     TDoubleArray&       buffer_xz_zz,
                                     TDoubleArray&       buffer_yy_xx,
                                     TDoubleArray&       buffer_yy_xy,
                                     TDoubleArray&       buffer_yy_xz,
                                     TDoubleArray&       buffer_yy_yy,
                                     TDoubleArray&       buffer_yy_yz,
                                     TDoubleArray&       buffer_yy_zz,
                                     TDoubleArray&       buffer_yz_xx,
                                     TDoubleArray&       buffer_yz_xy,
                                     TDoubleArray&       buffer_yz_xz,
                                     TDoubleArray&       buffer_yz_yy,
                                     TDoubleArray&       buffer_yz_yz,
                                     TDoubleArray&       buffer_yz_zz,
                                     TDoubleArray&       buffer_zz_xx,
                                     TDoubleArray&       buffer_zz_xy,
                                     TDoubleArray&       buffer_zz_xz,
                                     TDoubleArray&       buffer_zz_yy,
                                     TDoubleArray&       buffer_zz_yz,
                                     TDoubleArray&       buffer_zz_zz,
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

    auto fints_xx_xx = buffer_xx_xx.data();

    auto fints_xx_xy = buffer_xx_xy.data();

    auto fints_xx_xz = buffer_xx_xz.data();

    auto fints_xx_yy = buffer_xx_yy.data();

    auto fints_xx_yz = buffer_xx_yz.data();

    auto fints_xx_zz = buffer_xx_zz.data();

    auto fints_xy_xx = buffer_xy_xx.data();

    auto fints_xy_xy = buffer_xy_xy.data();

    auto fints_xy_xz = buffer_xy_xz.data();

    auto fints_xy_yy = buffer_xy_yy.data();

    auto fints_xy_yz = buffer_xy_yz.data();

    auto fints_xy_zz = buffer_xy_zz.data();

    auto fints_xz_xx = buffer_xz_xx.data();

    auto fints_xz_xy = buffer_xz_xy.data();

    auto fints_xz_xz = buffer_xz_xz.data();

    auto fints_xz_yy = buffer_xz_yy.data();

    auto fints_xz_yz = buffer_xz_yz.data();

    auto fints_xz_zz = buffer_xz_zz.data();

    auto fints_yy_xx = buffer_yy_xx.data();

    auto fints_yy_xy = buffer_yy_xy.data();

    auto fints_yy_xz = buffer_yy_xz.data();

    auto fints_yy_yy = buffer_yy_yy.data();

    auto fints_yy_yz = buffer_yy_yz.data();

    auto fints_yy_zz = buffer_yy_zz.data();

    auto fints_yz_xx = buffer_yz_xx.data();

    auto fints_yz_xy = buffer_yz_xy.data();

    auto fints_yz_xz = buffer_yz_xz.data();

    auto fints_yz_yy = buffer_yz_yy.data();

    auto fints_yz_yz = buffer_yz_yz.data();

    auto fints_yz_zz = buffer_yz_zz.data();

    auto fints_zz_xx = buffer_zz_xx.data();

    auto fints_zz_xy = buffer_zz_xy.data();

    auto fints_zz_xz = buffer_zz_xz.data();

    auto fints_zz_yy = buffer_zz_yy.data();

    auto fints_zz_yz = buffer_zz_yz.data();

    auto fints_zz_zz = buffer_zz_zz.data();

#pragma omp simd aligned(fints_xx_xx,     \
                             fints_xx_xy, \
                             fints_xx_xz, \
                             fints_xx_yy, \
                             fints_xx_yz, \
                             fints_xx_zz, \
                             fints_xy_xx, \
                             fints_xy_xy, \
                             fints_xy_xz, \
                             fints_xy_yy, \
                             fints_xy_yz, \
                             fints_xy_zz, \
                             fints_xz_xx, \
                             fints_xz_xy, \
                             fints_xz_xz, \
                             fints_xz_yy, \
                             fints_xz_yz, \
                             fints_xz_zz, \
                             fints_yy_xx, \
                             fints_yy_xy, \
                             fints_yy_xz, \
                             fints_yy_yy, \
                             fints_yy_yz, \
                             fints_yy_zz, \
                             fints_yz_xx, \
                             fints_yz_xy, \
                             fints_yz_xz, \
                             fints_yz_yy, \
                             fints_yz_yz, \
                             fints_yz_zz, \
                             fints_zz_xx, \
                             fints_zz_xy, \
                             fints_zz_xz, \
                             fints_zz_yy, \
                             fints_zz_yz, \
                             fints_zz_zz, \
                             ket_fe,      \
                             ket_fn,      \
                             ket_rx,      \
                             ket_ry,      \
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

        fints_xx_xx[i] +=
            fss * (-18.0 * fe_0 * rpb_y * tke_0 + 42.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 21.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                   42.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 42.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_xx_xx[i] +=
            fss * (42.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 - 42.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 -
                   24.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   9.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_xx[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_xx[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_xx[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_xx[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 36.0 * rpa_x * rpb_y * rpb_x * tke_0 + 18.0 * fe_0 * fe_0 * rpb_y * tke_0 * tke_0 +
                                 42.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 84.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xx_xx[i] +=
            fss * (-42.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 21.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                   84.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 42.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
                   42.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] +=
            fss *
            (-84.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             42.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 84.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             84.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xx[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (48.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (36.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (24.0 * fe_0 * rpa_x * rpb_y * rpb_x * tke_0 * tke_0 + 36.0 * fe_0 * rpb_y * rpb_x * rpb_x * tke_0 * tke_0 -
                                 14.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (-56.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 84.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 42.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (-84.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 84.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 84.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (36.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (-28.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * (32.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xx[i] += fss * 16.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xx_xy[i] += fss * (6.0 * rpa_x - 7.0 * fe_0 * rpa_x * tbe_0 - 14.0 * fe_0 * rpa_x * tbe_0 - 14.0 * rpa_x * rpa_x * rpa_x * tbe_0 +
                                 14.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss * (14.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 14.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                   3.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 + 4.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 + 8.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0);

        fints_xx_xy[i] += fss * (2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 + 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 + 4.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * rpa_x * tke_0 -
                                 6.0 * fe_0 * rpa_x * tke_0 - 12.0 * fe_0 * rpb_x * tke_0);

        fints_xx_xy[i] += fss * (7.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 14.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                                 7.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 14.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                                 28.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (14.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 +
                                 14.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 + 28.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 +
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss * (28.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 - 14.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                   14.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss * (-16.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] +=
            fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 - 12.0 * rpa_x * rpb_y * rpb_y * tke_0 -
                   12.0 * rpa_x * rpb_x * rpb_x * tke_0 + 6.0 * fe_0 * fe_0 * rpa_x * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (12.0 * fe_0 * fe_0 * rpb_x * tke_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                                 28.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 +
                                 28.0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss * (-7.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                   28.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 +
                   28.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (28.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 - 28.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-28.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 - 28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] +=
            fss *
            (-16.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 16.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_xy[i] += fss * (-8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * rpa_x * rpb_x * rpb_x * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (12.0 * fe_0 * rpa_x * rpb_y * rpb_y * tke_0 * tke_0 + 24.0 * fe_0 * rpb_y * rpb_y * rpb_x * tke_0 * tke_0 -
                                 14.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 14.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (-28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (-56.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (-28.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * (32.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xy[i] += fss * 16.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xx_xz[i] += fss * (-12.0 * rpa_x * rpb_z * rpb_y * tke_0 + 14.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 +
                                 28.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 28.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] +=
            fss *
            (-28.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 16.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * rpa_x * rpb_z * rpb_y * tke_0 * tke_0 + 24.0 * fe_0 * rpb_z * rpb_y * rpb_x * tke_0 * tke_0 -
                                 14.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (-28.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (-56.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (32.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (24.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_xz[i] += fss * (32.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] +=
            fss * (-18.0 * fe_0 * rpb_y * tke_0 + 42.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 21.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                   42.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 42.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_xx_yy[i] +=
            fss * (42.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 - 24.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] +=
            fss *
            (-24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 36.0 * rpa_x * rpb_y * rpb_x * tke_0 +
                                 6.0 * fe_0 * fe_0 * rpb_y * tke_0 * tke_0 + 12.0 * fe_0 * fe_0 * rpb_y * tke_0 * tke_0 +
                                 42.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xx_yy[i] +=
            fss * (84.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 14.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 -
                   28.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 7.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 -
                   14.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] +=
            fss *
            (84.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] +=
            fss *
            (-14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             24.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 - 48.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] +=
            fss * (12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tke_0 * tke_0 + 24.0 * fe_0 * rpa_x * rpb_y * rpb_x * tke_0 * tke_0 +
                   12.0 * fe_0 * rpb_y * rpb_y * rpb_y * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                   28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (-28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 14.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (-56.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (24.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yy[i] += fss * (32.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] +=
            fss * (-6.0 * fe_0 * rpb_z * tke_0 + 14.0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 + 7.0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 +
                   14.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0);

        fints_xx_yz[i] += fss * (14.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_yz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_yz[i] +=
            fss *
            (-8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_yz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xx_yz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 12.0 * rpa_x * rpb_z * rpb_x * tke_0 +
                                 6.0 * fe_0 * fe_0 * rpb_z * tke_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 +
                                 28.0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0);

        fints_xx_yz[i] +=
            fss * (-14.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 - 7.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 +
                   28.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 -
                   14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] +=
            fss *
            (-14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 - 16.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
             8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_yz[i] += fss * (-16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * rpa_x * rpb_z * rpb_x * tke_0 * tke_0);

        fints_xx_yz[i] +=
            fss *
            (12.0 * fe_0 * rpb_z * rpb_y * rpb_y * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 - 28.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             14.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (-28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (-56.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_yz[i] += fss * (16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] +=
            fss * (-6.0 * fe_0 * rpb_y * tke_0 + 14.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 7.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                   14.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_xx_zz[i] += fss * (14.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_zz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_zz[i] +=
            fss *
            (-8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_zz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xx_zz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * rpa_x * rpb_y * rpb_x * tke_0 +
                                 6.0 * fe_0 * fe_0 * rpb_y * tke_0 * tke_0 + 14.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 +
                                 28.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xx_zz[i] +=
            fss * (-14.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 7.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                   28.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
                   14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] +=
            fss *
            (-14.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 - 16.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
             8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xx_zz[i] += fss * (-16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tke_0 * tke_0);

        fints_xx_zz[i] +=
            fss *
            (12.0 * fe_0 * rpb_z * rpb_z * rpb_y * tke_0 * tke_0 - 14.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             28.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 28.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             14.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (-28.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tke_0 * tke_0 -
                                 28.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 56.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (-56.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xx_zz[i] += fss * (16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 32.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (9.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 9.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xx[i] +=
            fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xx[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 + 9.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 +
                                 18.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 + 18.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 18.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_xy_xx[i] +=
            fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                   18.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xx[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xx[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xx[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] +=
            fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                   18.0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 + 36.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xy_xx[i] +=
            fss *
            (-18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 36.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (-24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 36.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 36.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (36.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * (-24.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xx[i] += fss * 16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xy_xy[i] += fss * (-3.0 * fe_0 * rpa_y * tbe_0 - 6.0 * rpa_y * rpa_x * rpa_x * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_xy_xy[i] += fss * (4.0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] +=
            fss *
            (4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tke_0 +
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_xy_xy[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (-8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (6.0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 12.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xy_xy[i] +=
            fss *
            (12.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             24.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 24.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (-4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xy[i] += fss * (-8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (-12.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xy[i] += fss * (24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (3.0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xy_xz[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 +
                                 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xy_xz[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_xz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * (-24.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_xz[i] += fss * 16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xy_yy[i] += fss * (9.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 9.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_yy[i] +=
            fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_yy[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 9.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 +
                                 18.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 + 18.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_xy_yy[i] += fss * (18.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_yy[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_yy[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] +=
            fss *
            (36.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (-18.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (-36.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (24.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xy_yz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xy_yz[i] +=
            fss * (-4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tke_0 +
                   12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_yz[i] += fss * (-8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * (-24.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_yz[i] += fss * 16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xy_zz[i] += fss * (3.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_zz[i] +=
            fss * (-3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_zz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_xy_zz[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_zz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_zz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xy_zz[i] += fss * (-4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] +=
            fss *
            (12.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xy_zz[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xy_zz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (18.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_xx[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_xx[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 + 36.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (-36.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xx[i] += fss * (16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (-3.0 * fe_0 * rpa_z * tbe_0 - 6.0 * rpa_z * rpa_x * rpa_x * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_xz_xy[i] += fss * (4.0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] +=
            fss *
            (4.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tke_0);

        fints_xz_xy[i] += fss * (6.0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 +
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpb_x * rpb_x * tbe_0 * tke_0);

        fints_xz_xy[i] +=
            fss * (-3.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tke_0 * tke_0 + 12.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   12.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tke_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_xy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (3.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_xz[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xz_xz[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_xz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * (-24.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_xz[i] += fss * 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xz_yy[i] += fss * (18.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_yy[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_yy[i] +=
            fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 + 18.0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 +
                   36.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 18.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_yy[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * (-24.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yy[i] += fss * 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_xz_yz[i] += fss * (3.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_yz[i] +=
            fss * (-3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_yz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0);

        fints_xz_yz[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_yz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_yz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_xz_yz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] +=
            fss *
            (12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_yz[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_yz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_zz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xz_zz[i] +=
            fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 +
                   12.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_zz[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xz_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 24.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * (-24.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_xz_zz[i] += fss * 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_yy_xx[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 + 12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (-24.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (24.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (-fe_0 * rpa_x * tbe_0 - 2.0 * fe_0 * rpa_x * tbe_0 - 2.0 * rpa_x * rpa_x * rpa_x * tbe_0 +
                                 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 +
                   2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0);

        fints_yy_xy[i] += fss * (4.0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 + 4.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 + fe_0 * fe_0 * rpa_x * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 +
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 +
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (-8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (4.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (-8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_xy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] +=
            fss *
            (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (-8.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (-8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (-4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xy[i] += fss * (8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_xz[i] += fss * (-4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 +
                                 4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_yy_xz[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] +=
            fss *
            (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_yy[i] += fss * (6.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yy_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_yy[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_yy[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] +=
            fss * (-fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                   12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_yy[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (-8.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (2.0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 + fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0);

        fints_yy_yz[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_yz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yy_yz[i] += fss * (-4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (2.0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 +
                                 4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0);

        fints_yy_yz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (-2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_yz[i] += fss * 16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_yy_zz[i] +=
            fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                   fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (2.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (-4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] +=
            fss *
            (4.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (-8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yy_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] +=
            fss *
            (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yy_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xx[i] +=
            fss *
            (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_xx[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xx[i] += fss * (16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 + 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 +
                   4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xy[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_xy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (-fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] +=
            fss *
            (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_xz[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_xz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_xz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_yy[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_yy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yy[i] += fss * 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_yz_yz[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_yz[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_yz_yz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (-8.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_yz[i] += fss * (16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_zz[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yz_zz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yz_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (-8.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_yz_zz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (6.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_zz_xx[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_xx[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_xx[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_zz_xx[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                   12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xx[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xx[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xx[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] +=
            fss *
            (24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (-12.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (-8.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xx[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (-fe_0 * rpa_x * tbe_0 - 2.0 * fe_0 * rpa_x * tbe_0 - 2.0 * rpa_x * rpa_x * rpa_x * tbe_0 +
                                 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0);

        fints_zz_xy[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 +
                   2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0);

        fints_zz_xy[i] += fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_zz_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0);

        fints_zz_xy[i] +=
            fss * (fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 + 4.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0);

        fints_zz_xy[i] += fss * (2.0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 +
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xy[i] +=
            fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
             2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xy[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (2.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 -
                                 fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] +=
            fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tke_0 * tke_0 + 4.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   4.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0);

        fints_zz_xy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xy[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_xy[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] +=
            fss *
            (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (-8.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (-8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (-8.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xy[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_xz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 +
                                 4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_zz_xz[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] +=
            fss *
            (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_xz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (6.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_zz_yy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_yy[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_yy[i] += fss * (6.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 + 12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yy[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
             8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yy[i] += fss * 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;

        fints_zz_yz[i] +=
            fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (2.0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 + fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] +=
            fss *
            (4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] +=
            fss *
            (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 - 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_yz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (2.0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_zz_zz[i] += fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_zz[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zz_zz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (2.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tke_0 * tke_0 +
                                 4.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_zz_zz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tke_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zz_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (-2.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0 -
                                 8.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0);

        fints_zz_zz[i] += fss * 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 * tke_0;
    }
}

}  // namespace ovlrec
