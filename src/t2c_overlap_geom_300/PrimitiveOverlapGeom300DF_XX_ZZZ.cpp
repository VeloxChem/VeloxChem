#include "PrimitiveOverlapGeom300DF_XX_ZZZ.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapGeom300DF_XX_ZZZ(TDoubleArray&       buffer_xxx,
                                     TDoubleArray&       buffer_xxy,
                                     TDoubleArray&       buffer_xxz,
                                     TDoubleArray&       buffer_xyy,
                                     TDoubleArray&       buffer_xyz,
                                     TDoubleArray&       buffer_xzz,
                                     TDoubleArray&       buffer_yyy,
                                     TDoubleArray&       buffer_yyz,
                                     TDoubleArray&       buffer_yzz,
                                     TDoubleArray&       buffer_zzz,
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

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
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

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto tbe_0 = bra_exp;

        fints_xxx[i] +=
            fss * (12.0 * fe_0 * rpa_x * rpb_z * tbe_0 + 24.0 * fe_0 * rpa_x * rpb_z * tbe_0 - 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   18.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - 18.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xxx[i] +=
            fss * (-36.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - 18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   36.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxx[i] +=
            fss *
            (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 16.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxx[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxx[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   16.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 24.0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0);

        fints_xxx[i] +=
            fss * (-18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 36.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   36.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxx[i] += fss * (16.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               16.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxx[i] += fss * 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0;

        fints_xxy[i] +=
            fss * (2.0 * fe_0 * rpa_y * rpb_z * tbe_0 + 4.0 * fe_0 * rpa_y * rpb_z * tbe_0 - 5.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 -
                   10.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 - 10.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xxy[i] +=
            fss * (-20.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxy[i] += fss * (8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 - 10.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                               20.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxy[i] += fss * (4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxz[i] += fss * (3.0 * fe_0 * fe_0 * tbe_0 - (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 -
                               15.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 + (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_z * rpb_z * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpb_z * tbe_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * rpb_z * rpb_z * tbe_0 - 5.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 -
                               10.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 - 15.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 -
                               10.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xxz[i] +=
            fss *
            (-20.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 - 30.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 +
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             9.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               18.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 - 10.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                               20.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyy[i] +=
            fss * (2.0 * fe_0 * rpa_x * rpb_z * tbe_0 + 4.0 * fe_0 * rpa_x * rpb_z * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xyy[i] += fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                               4.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 -
                               8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xyy[i] +=
            fss * (-2.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 +
                   fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyy[i] += fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyy[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 4.0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0);

        fints_xyy[i] +=
            fss * (-4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 8.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xyy[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyy[i] += fss * 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0;

        fints_xyz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                               4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xyz[i] +=
            fss * (-8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyz[i] += fss * (12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                               8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] += fss * (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 -
                               12.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] +=
            fss *
            (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 +
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_x * rpb_z * tbe_0 +
                               4.0 * fe_0 * rpa_x * rpb_z * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xzz[i] += fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                               2.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                               4.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xzz[i] +=
            fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   2.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_xzz[i] +=
            fss *
            (-4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 + fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 4.0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0);

        fints_xzz[i] +=
            fss * (-4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 8.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xzz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzz[i] += fss * 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0;

        fints_yyy[i] +=
            fss * (-3.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 -
                   6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 +
                   fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyy[i] +=
            fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyy[i] += fss * (4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                               fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 -
                               3.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_yyz[i] +=
            fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 - 6.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 +
                   fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   3.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                               2.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_yyz[i] += fss * (-4.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yzz[i] +=
            fss * (12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 - fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_yzz[i] +=
            fss * (fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzz[i] +=
            fss * (12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   2.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 4.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                   2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzz[i] += fss * (4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 -
                               9.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] +=
            fss *
            ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
             3.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 +
             3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] +=
            fss *
            (6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 +
             6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (-3.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 -
                               9.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 - 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 -
                               12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0);

        fints_zzz[i] +=
            fss *
            (-18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 + fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               6.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] +=
            fss * (12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 - 12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_zzz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                               4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzz[i] += fss * 8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0;
    }
}

}  // namespace ovlrec
