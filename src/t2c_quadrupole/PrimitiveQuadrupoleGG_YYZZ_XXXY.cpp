#include "PrimitiveQuadrupoleGG_YYZZ_XXXY.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace mpol {  // mpol namespace

auto
compPrimitiveQuadrupoleGG_YYZZ_XXXY(TDoubleArray&       buffer_xx,
                                    TDoubleArray&       buffer_xy,
                                    TDoubleArray&       buffer_xz,
                                    TDoubleArray&       buffer_yy,
                                    TDoubleArray&       buffer_yz,
                                    TDoubleArray&       buffer_zz,
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

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_xx[i] +=
            fss *
            ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x +
             (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xx[i] += fss * (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x;

        fints_xx[i] +=
            faa_x *
            (3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xx[i] += faa_x * ((3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                                (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                                (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x);

        fints_xx[i] += faa_x * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xx[i] +=
            faa_xx * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_xx[i] +=
            faa_xx * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xx[i] += faa_xx * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xy[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x +
                              (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x +
                              (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y +
                              (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y +
                              (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x);

        fints_xy[i] +=
            fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (9.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xy[i] += fss * ((9.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x + (9.0 / 32.0) * fe_0 * fe_0 * fe_0 * fe_0 * fe_0);

        fints_xy[i] +=
            faa_y *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xy[i] += faa_y * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                                (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                                (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x);

        fints_xy[i] += faa_y * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xy[i] += faa_x * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x +
                                (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
                                (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x +
                                (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x +
                                (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_xy[i] +=
            faa_x * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x +
                     (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_xy[i] += faa_x * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (9.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xy[i] +=
            faa_xy * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_xy[i] +=
            faa_xy * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xy[i] += faa_xy * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xz[i] +=
            fss *
            ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y +
             (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_xz[i] += fss * (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y;

        fints_xz[i] +=
            faa_z *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xz[i] += faa_z * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                                (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                                (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x);

        fints_xz[i] += faa_z * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xz[i] +=
            faa_x * (fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x +
                     (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                     (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_xz[i] += faa_x * (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x;

        fints_xz[i] +=
            faa_xz * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_xz[i] +=
            faa_xz * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xz[i] += faa_xz * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_yy[i] +=
            fss *
            (fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x +
             (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_yy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                              (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_yy[i] +=
            faa_y *
            (2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
             3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_yy[i] +=
            faa_y * (fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
                     (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yy[i] += faa_y * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_yy[i] +=
            faa_yy * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_yy[i] +=
            faa_yy * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_yy[i] += faa_yy * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_yz[i] +=
            fss *
            (fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x +
             (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_yz[i] += fss * (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x;

        fints_yz[i] += faa_z * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x +
                                (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
                                (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x +
                                (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x +
                                (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_yz[i] +=
            faa_z * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x +
                     (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                     (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yz[i] += faa_z * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (9.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_yz[i] +=
            faa_y * (fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x +
                     (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                     (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_yz[i] += faa_y * (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x;

        fints_yz[i] +=
            faa_yz * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_yz[i] +=
            faa_yz * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_yz[i] += faa_yz * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_zz[i] +=
            fss *
            ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x +
             (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x +
             (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_zz[i] += fss * (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x;

        fints_zz[i] +=
            faa_z * (2.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x +
                     3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x + fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                     3.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_zz[i] += faa_z * (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x;

        fints_zz[i] +=
            faa_zz * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x +
                      (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_zz[i] +=
            faa_zz * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x +
                      (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_zz[i] += faa_zz * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);
    }
}

}  // namespace mpol
