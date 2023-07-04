#include "PrimitiveOctupoleFG_YYY_XXZZ.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace mpol {  // mpol namespace

auto
compPrimitiveOctupoleFG_YYY_XXZZ(TDoubleArray&       buffer_xxx,
                                 TDoubleArray&       buffer_xxy,
                                 TDoubleArray&       buffer_xxz,
                                 TDoubleArray&       buffer_xyy,
                                 TDoubleArray&       buffer_xyz,
                                 TDoubleArray&       buffer_xzz,
                                 TDoubleArray&       buffer_yyy,
                                 TDoubleArray&       buffer_yyz,
                                 TDoubleArray&       buffer_yzz,
                                 TDoubleArray&       buffer_zzz,
                                 const TPoint3D&     point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = rpc_x * fss;

        const auto faa_y = rpc_y * fss;

        const auto faa_z = rpc_z * fss;

        const auto faa_xx = fss * (rpc_x * rpc_x + 0.5 * fe_0);

        const auto faa_xy = fss * rpc_x * rpc_y;

        const auto faa_xz = fss * rpc_x * rpc_z;

        const auto faa_yy = fss * (rpc_y * rpc_y + 0.5 * fe_0);

        const auto faa_yz = fss * rpc_x * rpc_z;

        const auto faa_zz = fss * (rpc_z * rpc_z + 0.5 * fe_0);

        const auto faa_xxx = fss * (rpc_x * rpc_x * rpc_x + 1.5 * fe_0 * rpc_x);

        const auto faa_xxy = fss * (rpc_x * rpc_x * rpc_y + 0.5 * fe_0 * rpc_y);

        const auto faa_xxz = fss * (rpc_x * rpc_x * rpc_z + 0.5 * fe_0 * rpc_z);

        const auto faa_xyy = fss * (rpc_x * rpc_y * rpc_y + 0.5 * fe_0 * rpc_x);

        const auto faa_xyz = fss * rpc_x * rpc_y * rpc_z;

        const auto faa_xzz = fss * (rpc_x * rpc_z * rpc_z + 0.5 * fe_0 * rpc_x);

        const auto faa_yyy = fss * (rpc_y * rpc_y * rpc_y + 1.5 * fe_0 * rpc_y);

        const auto faa_yyz = fss * (rpc_y * rpc_y * rpc_z + 0.5 * fe_0 * rpc_z);

        const auto faa_yzz = fss * (rpc_y * rpc_z * rpc_z + 0.5 * fe_0 * rpc_y);

        const auto faa_zzz = fss * (rpc_z * rpc_z * rpc_z + 1.5 * fe_0 * rpc_z);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            faa_x * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xxx[i] +=
            faa_xx * (3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                      (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxx[i] +=
            faa_xxx * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxx[i] += faa_xxx * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xxy[i] +=
            fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0);

        fints_xxy[i] +=
            faa_y * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                     (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xxy[i] +=
            faa_x * (3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                     (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xxy[i] += faa_xy * (2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                                  fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxy[i] +=
            faa_xx * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x +
                      (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxy[i] += faa_xx * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x +
                                  (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_xxy[i] +=
            faa_xxy * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxy[i] += faa_xxy * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] +=
            fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xxz[i] +=
            faa_z * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                     (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xxz[i] += faa_x * (2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x + 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xxz[i] += faa_xz * (2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                                  fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxz[i] += faa_xx * (fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                                  (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xxz[i] +=
            faa_xxz * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxz[i] += faa_xxz * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xyy[i] +=
            faa_y * (3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                     (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xyy[i] += faa_yy * (fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                                  (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xyy[i] +=
            faa_x * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] +=
            faa_xy * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                      (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xyy[i] += faa_xy * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x +
                                  (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_xyy[i] +=
            faa_xyy * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xyy[i] += faa_xyy * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xyz[i] +=
            faa_z * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xyz[i] += faa_y * (fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xyz[i] += faa_yz * (fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                                  (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xyz[i] +=
            faa_x * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_xyz[i] +=
            faa_xz * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x +
                      (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xyz[i] += faa_xz * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x +
                                  (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_xyz[i] += faa_xy * (fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                                  (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xyz[i] +=
            faa_xyz * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xyz[i] += faa_xyz * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xzz[i] +=
            fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xzz[i] += faa_z * (2.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x + 3.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xzz[i] += faa_zz * (fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                                  (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xzz[i] +=
            faa_x * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                     (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xzz[i] += faa_xz * (2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                                  fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xzz[i] +=
            faa_xzz * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xzz[i] += faa_xzz * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yyy[i] +=
            faa_y * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                     (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] +=
            faa_yy * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                      (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x +
                      (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yyy[i] += faa_yy * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x +
                                  (9.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yyy[i] +=
            faa_yyy * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_yyy[i] += faa_yyy * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_yyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_yyz[i] +=
            faa_z * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] +=
            faa_y * (3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                     (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_yyz[i] +=
            faa_yz * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                      (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yyz[i] += faa_yz * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x +
                                  (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yyz[i] += faa_yy * (fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                                  (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_yyz[i] +=
            faa_yyz * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_yyz[i] += faa_yyz * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_yzz[i] +=
            fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yzz[i] +=
            faa_z * (3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                     (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_yzz[i] +=
            faa_zz * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x +
                      (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yzz[i] += faa_zz * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x +
                                  (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yzz[i] +=
            faa_y * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                     (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yzz[i] += faa_yz * (2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                                  fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_yzz[i] +=
            faa_yzz * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_yzz[i] += faa_yzz * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_zzz[i] +=
            faa_z * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_zzz[i] +=
            faa_zz * (3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                      (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_zzz[i] +=
            faa_zzz * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                       (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_zzz[i] += faa_zzz * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                                   rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);
    }
}

}  // namespace mpol
