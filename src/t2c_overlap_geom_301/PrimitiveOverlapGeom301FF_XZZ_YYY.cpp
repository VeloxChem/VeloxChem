#include "PrimitiveOverlapGeom301FF_XZZ_YYY.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapGeom301FF_XZZ_YYY(TDoubleArray&       buffer_xxx_x,
                                      TDoubleArray&       buffer_xxx_y,
                                      TDoubleArray&       buffer_xxx_z,
                                      TDoubleArray&       buffer_xxy_x,
                                      TDoubleArray&       buffer_xxy_y,
                                      TDoubleArray&       buffer_xxy_z,
                                      TDoubleArray&       buffer_xxz_x,
                                      TDoubleArray&       buffer_xxz_y,
                                      TDoubleArray&       buffer_xxz_z,
                                      TDoubleArray&       buffer_xyy_x,
                                      TDoubleArray&       buffer_xyy_y,
                                      TDoubleArray&       buffer_xyy_z,
                                      TDoubleArray&       buffer_xyz_x,
                                      TDoubleArray&       buffer_xyz_y,
                                      TDoubleArray&       buffer_xyz_z,
                                      TDoubleArray&       buffer_xzz_x,
                                      TDoubleArray&       buffer_xzz_y,
                                      TDoubleArray&       buffer_xzz_z,
                                      TDoubleArray&       buffer_yyy_x,
                                      TDoubleArray&       buffer_yyy_y,
                                      TDoubleArray&       buffer_yyy_z,
                                      TDoubleArray&       buffer_yyz_x,
                                      TDoubleArray&       buffer_yyz_y,
                                      TDoubleArray&       buffer_yyz_z,
                                      TDoubleArray&       buffer_yzz_x,
                                      TDoubleArray&       buffer_yzz_y,
                                      TDoubleArray&       buffer_yzz_z,
                                      TDoubleArray&       buffer_zzz_x,
                                      TDoubleArray&       buffer_zzz_y,
                                      TDoubleArray&       buffer_zzz_z,
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

    auto fints_xxx_x = buffer_xxx_x.data();

    auto fints_xxx_y = buffer_xxx_y.data();

    auto fints_xxx_z = buffer_xxx_z.data();

    auto fints_xxy_x = buffer_xxy_x.data();

    auto fints_xxy_y = buffer_xxy_y.data();

    auto fints_xxy_z = buffer_xxy_z.data();

    auto fints_xxz_x = buffer_xxz_x.data();

    auto fints_xxz_y = buffer_xxz_y.data();

    auto fints_xxz_z = buffer_xxz_z.data();

    auto fints_xyy_x = buffer_xyy_x.data();

    auto fints_xyy_y = buffer_xyy_y.data();

    auto fints_xyy_z = buffer_xyy_z.data();

    auto fints_xyz_x = buffer_xyz_x.data();

    auto fints_xyz_y = buffer_xyz_y.data();

    auto fints_xyz_z = buffer_xyz_z.data();

    auto fints_xzz_x = buffer_xzz_x.data();

    auto fints_xzz_y = buffer_xzz_y.data();

    auto fints_xzz_z = buffer_xzz_z.data();

    auto fints_yyy_x = buffer_yyy_x.data();

    auto fints_yyy_y = buffer_yyy_y.data();

    auto fints_yyy_z = buffer_yyy_z.data();

    auto fints_yyz_x = buffer_yyz_x.data();

    auto fints_yyz_y = buffer_yyz_y.data();

    auto fints_yyz_z = buffer_yyz_z.data();

    auto fints_yzz_x = buffer_yzz_x.data();

    auto fints_yzz_y = buffer_yzz_y.data();

    auto fints_yzz_z = buffer_yzz_z.data();

    auto fints_zzz_x = buffer_zzz_x.data();

    auto fints_zzz_y = buffer_zzz_y.data();

    auto fints_zzz_z = buffer_zzz_z.data();

#pragma omp simd aligned(fints_xxx_x,     \
                             fints_xxx_y, \
                             fints_xxx_z, \
                             fints_xxy_x, \
                             fints_xxy_y, \
                             fints_xxy_z, \
                             fints_xxz_x, \
                             fints_xxz_y, \
                             fints_xxz_z, \
                             fints_xyy_x, \
                             fints_xyy_y, \
                             fints_xyy_z, \
                             fints_xyz_x, \
                             fints_xyz_y, \
                             fints_xyz_z, \
                             fints_xzz_x, \
                             fints_xzz_y, \
                             fints_xzz_z, \
                             fints_yyy_x, \
                             fints_yyy_y, \
                             fints_yyy_z, \
                             fints_yyz_x, \
                             fints_yyz_y, \
                             fints_yyz_z, \
                             fints_yzz_x, \
                             fints_yzz_y, \
                             fints_yzz_z, \
                             fints_zzz_x, \
                             fints_zzz_y, \
                             fints_zzz_z, \
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

        fints_xxx_x[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xxx_x[i] +=
            fss *
            (6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (12.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] += fss * (-48.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_x[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (-(9.0 / 2.0) * fe_0 * fe_0 * tbe_0 - 9.0 * fe_0 * rpa_z * rpa_z * tbe_0 + 9.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 +
                                 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 + 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0);

        fints_xxx_y[i] +=
            fss *
            (-(9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
             9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxx_y[i] += fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxx_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 9.0 * fe_0 * rpb_y * rpb_y * tbe_0 + (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 -
                                 18.0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 + 9.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0);

        fints_xxx_y[i] +=
            fss * (18.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 - 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 +
                   36.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 -
                   18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (-18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 72.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xxx_y[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xxx_y[i] +=
            fss * (9.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (-12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] +=
            fss * (24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + 12.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   24.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_y[i] +=
            fss * (24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (3.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (-12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 +
                                 6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xxx_z[i] +=
            fss * (12.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                   6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] +=
            fss * (6.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + 12.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   12.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   24.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   24.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxx_z[i] += fss * 16.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0;

        fints_xxy_x[i] += fss * (-(9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_x[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] +=
            fss * (9.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 + 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxy_y[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xxy_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxy_y[i] +=
            fss *
            (-24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 +
             18.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (36.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xxy_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xxy_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xxy_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (-48.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_y[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxy_z[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (12.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] += fss * (-16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_x[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-18.0 * fe_0 * rpa_z * rpa_x * tbe_0 + 9.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 +
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0);

        fints_xxz_y[i] +=
            fss *
            (18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxz_y[i] += fss * (-12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxz_y[i] +=
            fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   36.0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 + 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tke_0 +
                   18.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 + 36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xxz_y[i] +=
            fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 +
                   24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xxz_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] +=
            fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                   24.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (-16.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_y[i] +=
            fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 12.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (12.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] += fss * (-16.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xxz_z[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] +=
            fss * (-fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xyy_x[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 + 4.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (-2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_x[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss * (3.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - (3.0 / 2.0) * fe_0 * fe_0 * tbe_0 -
                                 3.0 * fe_0 * rpa_z * rpa_z * tbe_0 + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 +
                                 (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] += fss * (3.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss * (3.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 -
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0);

        fints_xyy_y[i] +=
            fss * (12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] += fss * (-3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] +=
            fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                   24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * rpb_y * rpb_y * tbe_0 +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 - 6.0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0);

        fints_xyy_y[i] += fss * (3.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 + 3.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss *
            (6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss *
            (-8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
             6.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xyy_y[i] +=
            fss *
            (6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
             (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (12.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyy_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss * (3.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (4.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (-4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_y[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (-fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (-2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (-fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (2.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + 4.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] += fss * (-8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyy_z[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] +=
            fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xyz_x[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_x[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] +=
            fss * (-6.0 * fe_0 * rpa_z * rpa_y * tbe_0 - 12.0 * fe_0 * rpa_z * rpb_y * tbe_0 + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0);

        fints_xyz_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0);

        fints_xyz_y[i] +=
            fss *
            (12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyz_y[i] += fss * (-12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyz_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0);

        fints_xyz_y[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 + 8.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xyz_y[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] +=
            fss *
            (-16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 - 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (12.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xyz_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyz_y[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] +=
            fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (16.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (-32.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (8.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] += fss * (-16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_y[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] +=
            fss * (3.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 -
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0);

        fints_xyz_z[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             3.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xyz_z[i] +=
            fss *
            (8.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
             12.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xyz_z[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 5.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-10.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 5.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-10.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpb_y * rpb_x * tke_0 - 4.0 * fe_0 * rpb_y * rpb_x * tke_0 + 5.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (10.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 + 10.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 +
                                 20.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xzz_x[i] +=
            fss * (-3.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xzz_x[i] +=
            fss * (4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-20.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-40.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                   4.0 * rpb_y * rpb_y * rpb_y * rpb_x * tke_0 + 10.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 + 20.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-6.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] += fss * (-40.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_x[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (3.0 * fe_0 - (15.0 / 2.0) * fe_0 * fe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * tbe_0 - 15.0 * fe_0 * rpa_z * rpa_z * tbe_0 +
                                 (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0);

        fints_xzz_y[i] += fss * (-6.0 * fe_0 * rpa_x * rpa_x * tbe_0 + (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 +
                                 9.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0);

        fints_xzz_y[i] +=
            fss * (15.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 -
                   (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 +
                   30.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_xzz_y[i] +=
            fss *
            (-9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xzz_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 + 6.0 * rpb_y * rpb_y);

        fints_xzz_y[i] +=
            fss * (-3.0 * fe_0 * fe_0 * tke_0 - 15.0 * fe_0 * rpb_y * rpb_y * tbe_0 + (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * rpb_y * rpb_y * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-30.0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 - (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0);

        fints_xzz_y[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tke_0 + 15.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 -
                   (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xzz_y[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 -
                   9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 + 30.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_xzz_y[i] +=
            fss *
            (30.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 - 15.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
             15.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
             (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (12.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                 60.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xzz_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_xzz_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-24.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * rpb_y * rpb_y * tke_0 - 4.0 * fe_0 * rpb_y * rpb_y * tke_0 - 6.0 * fe_0 * rpb_y * rpb_y * tke_0);

        fints_xzz_y[i] += fss * (5.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 10.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 +
                                 15.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xzz_y[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 10.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                   20.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 + 30.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   9.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-5.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 15.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-10.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-40.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 60.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * rpb_y * rpb_y * rpb_y * rpb_y * tke_0);

        fints_xzz_y[i] +=
            fss * (10.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                   20.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                   8.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-10.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (-20.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_y[i] +=
            fss * (24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (5.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 +
                                 5.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             10.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             10.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * rpb_z * rpb_y * tke_0 - 4.0 * fe_0 * rpb_z * rpb_y * tke_0 + 5.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (10.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 + 10.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 +
                                 20.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss * (10.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + 10.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss * (8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 - 5.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-10.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-10.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-40.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                   4.0 * rpb_z * rpb_y * rpb_y * rpb_y * tke_0 + 10.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 + 20.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] += fss * (-40.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_xzz_z[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] += fss * (-24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_x[i] +=
            fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yyy_y[i] +=
            fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 + 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 +
                   18.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 +
                   36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0);

        fints_yyy_y[i] +=
            fss *
            (-3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
             12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yyy_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_yyy_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_yyy_y[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_yyy_y[i] += fss * (-9.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_yyy_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] += fss * (32.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_y[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] += fss * (-24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyy_z[i] +=
            fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 -
                                 fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] +=
            fss *
            (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
             4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 +
             4.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (8.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] += fss * (-16.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_x[i] +=
            fss * (16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * rpa_z * rpa_x * tbe_0);

        fints_yyz_y[i] += fss * (3.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0);

        fints_yyz_y[i] +=
            fss *
            (24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 -
             24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 - 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yyz_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yyz_y[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yyz_y[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tke_0);

        fints_yyz_y[i] +=
            fss * (6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                   3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 + 12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_yyz_y[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_yyz_y[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-24.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] +=
            fss * (32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-32.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 32.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (32.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] += fss * (-8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_y[i] +=
            fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yyz_z[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
             2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 16.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yyz_z[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 - (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tke_0 -
                                 15.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tke_0 +
                                 (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 15.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 5.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (-15.0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] +=
            fss * (4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (-20.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (-60.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] += fss * (24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_yzz_x[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] +=
            fss * (-6.0 * fe_0 * rpa_y * rpa_x * tbe_0 - 12.0 * fe_0 * rpa_x * rpb_y * tbe_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 +
                   30.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 + 30.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0);

        fints_yzz_y[i] +=
            fss * (60.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 - 9.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yzz_y[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yzz_y[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 16.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_yzz_y[i] +=
            fss *
            (30.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 - 15.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
             20.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 - 40.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
             60.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_yzz_y[i] += fss * (-30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 80.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_yzz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_yzz_y[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_yzz_y[i] +=
            fss * (12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 16.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   20.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                   30.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (-40.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 60.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 80.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] += fss * (8.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_y[i] +=
            fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (-15.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                                 15.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] +=
            fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   6.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tke_0 - 15.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (-20.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 10.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (-30.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 8.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 +
                   12.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tke_0 -
                   10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (-20.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 30.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (-20.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 60.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 20.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 40.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_yzz_z[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + 24.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (-36.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_zzz_x[i] +=
            fss * (48.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 + 24.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                   18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                   36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (-72.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (48.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_x[i] +=
            fss * (32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_x * tbe_0 + 27.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 +
                                 54.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 + 54.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_zzz_y[i] += fss * (-12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_zzz_y[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                 72.0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 + 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tke_0 +
                                 54.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0);

        fints_zzz_y[i] +=
            fss * (108.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 - 27.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 -
                   54.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0 +
                   108.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 -
                   54.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 +
                                 9.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (24.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 12.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 -
                                 48.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_zzz_y[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 -
                                 24.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0);

        fints_zzz_y[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 + 48.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   72.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 -
                   18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (-36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 54.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 108.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 108.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (18.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] += fss * (-72.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_y[i] +=
            fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 24.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (-18.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 + 48.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (-18.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (-36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (16.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 24.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 48.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tke_0 -
                                 36.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 -
                                 72.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] += fss * (-72.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tke_0 +
                                 12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 16.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 32.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                                 8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);

        fints_zzz_z[i] +=
            fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   24.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   32.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0 +
                   16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * tbe_0 * tbe_0 * tbe_0 * tke_0);
    }
}

}  // namespace ovlrec
