#include "PrimitiveNuclearPotentialGeom010GG_XXXY_ZZZZ.hpp"

#include <cmath>

#include "BoysFunc.hpp"
#include "MathConst.hpp"

namespace geom_npotrec {  // geom_npotrec namespace

auto
compPrimitiveNuclearPotentialGeom010GG_XXXY_ZZZZ(TDoubleArray&       buffer_x,
                                                 TDoubleArray&       buffer_y,
                                                 TDoubleArray&       buffer_z,
                                                 const TPoint3D&     dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<9> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<10> bf_values;

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto b7_vals = bf_values[7].data();

    auto b8_vals = bf_values[8].data();

    auto b9_vals = bf_values[9].data();

    auto targs = bf_args.data();

    // compute Boys function values

#pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<10>(bf_values, bf_args, ket_dim);

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
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

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = dip_x * 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = dip_y * 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = dip_z * 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] *
                      (9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z +
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_x[i] += dip_x * fss * b1_vals[i] * 3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_x[i] += dip_x * fss * b2_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x - 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z - 6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z -
                       (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_x[i] += dip_x * fss * b2_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x -
                       9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_x[i] += dip_x * fss * b2_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y);

        fints_x[i] += dip_x * fss * b2_vals[i] *
                      (-6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x - 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z -
                       3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] *
                      (36.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x + 18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x +
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z + 9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z +
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] *
                      (9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z +
                       18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y +
                       9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] *
                      (6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y +
                       9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x +
                       18.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] *
                      ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z +
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y);

        fints_x[i] +=
            dip_x * fss * b3_vals[i] *
            (9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + 9.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
             (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 24.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * fss * b3_vals[i] *
                      (18.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + 3.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x +
                       6.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x + 12.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-36.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x - 18.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x -
                       9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x -
                       36.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x -
                       18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y -
                       9.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x - 9.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x -
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-18.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 36.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
                       12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - 12.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
                       24.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] *
                      (-18.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - 3.0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b5_vals[i] *
                      (18.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x + 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x +
                       (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b5_vals[i] *
                      (36.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x + 18.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x +
                       9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + 18.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                       6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b5_vals[i] *
                      (9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x);

        fints_x[i] += dip_x * fss * b5_vals[i] *
                      (9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + 9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y +
                       24.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * fss * b5_vals[i] *
                      (3.0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + 18.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
                       36.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 12.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y +
                       12.0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b6_vals[i] *
                      (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z -
                       18.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - 18.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
                       6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b6_vals[i] *
                      (-9.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x -
                       6.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * fss * b6_vals[i] *
                      (-12.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 24.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
                       3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 18.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b7_vals[i] *
                      (9.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
                       3.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + 6.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       12.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b8_vals[i] * (-3.0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += faa_x * b1_vals[i] *
                      ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z +
                       3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_x[i] += faa_x * b1_vals[i] * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_x[i] +=
            faa_x * b2_vals[i] *
            (-6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z -
             9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z -
             3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z);

        fints_x[i] +=
            faa_x * b2_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y -
             9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z);

        fints_x[i] += faa_x * b2_vals[i] *
                      (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y);

        fints_x[i] += faa_x * b2_vals[i] *
                      (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - 3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
                       4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_x[i] += faa_x * b2_vals[i] * (-rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_x[i] += faa_x * b3_vals[i] *
                      (9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_x[i] +=
            faa_x * b3_vals[i] *
            (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z +
             6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x +
             6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_x[i] +=
            faa_x * b3_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y + 9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b3_vals[i] *
                      (18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_x[i] += faa_x * b3_vals[i] *
                      ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y +
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_x[i] += faa_x * b3_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x +
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_x[i] += faa_x * b3_vals[i] *
                      ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x +
                       3.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x +
                       12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x +
                       6.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_x[i] += faa_x * b3_vals[i] *
                      (3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
                       4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x -
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y - 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x -
                       9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y -
                       6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
                       (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x -
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - 18.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                       18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y - 9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x -
                       9.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b4_vals[i] *
                      (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x -
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x -
                       12.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_x[i] +=
            faa_x * b4_vals[i] *
            (-18.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
             4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b4_vals[i] * (-6.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_x[i] +=
            faa_x * b5_vals[i] *
            (18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z +
             9.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_x[i] += faa_x * b5_vals[i] *
                      (6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x +
                       18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += faa_x * b5_vals[i] *
                      (6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y +
                       9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x +
                       9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b5_vals[i] *
                      (3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z);

        fints_x[i] += faa_x * b5_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_x[i] += faa_x * b5_vals[i] *
                      (9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x +
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + 18.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b5_vals[i] *
                      ((9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x +
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y +
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_x[i] +=
            faa_x * b5_vals[i] *
            (18.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
             12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             4.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
             12.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] +=
            faa_x * b5_vals[i] *
            (18.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] +=
            faa_x * b6_vals[i] *
            (-9.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z -
             6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_x[i] +=
            faa_x * b6_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_x[i] += faa_x * b6_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x -
                       3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] += faa_x * b6_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                       9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_x[i] += faa_x * b6_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 12.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
                       3.0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_x[i] +=
            faa_x * b6_vals[i] *
            (-6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             18.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             4.0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] +=
            faa_x * b7_vals[i] *
            (3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x +
             9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] += faa_x * b7_vals[i] *
                      (6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                       (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] += faa_x * b7_vals[i] *
                      (3.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
                       4.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                       12.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                       3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       6.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] +=
            faa_x * b8_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             4.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] += faa_x * b9_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] *
                      ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x +
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_y[i] += dip_y * fss * b1_vals[i] * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_y[i] += dip_y * fss * b2_vals[i] *
                      (-6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z -
                       9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x - 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z -
                       3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z);

        fints_y[i] += dip_y * fss * b2_vals[i] *
                      (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x - 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z -
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_x -
                       (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x);

        fints_y[i] += dip_y * fss * b2_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_x - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
                       4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] *
                      (9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] *
                      (6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z +
                       6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x +
                       18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] *
                      ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z +
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_x +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x);

        fints_y[i] +=
            dip_y * fss * b3_vals[i] *
            (9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x + 9.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpc_x + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
             (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x + 3.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] *
                      (12.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + 6.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x -
                       18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z -
                       9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_x -
                       18.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x * rpc_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x -
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] *
                      (-12.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 18.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
                       4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b5_vals[i] *
                      (18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z +
                       9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
                       9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b5_vals[i] *
                      (6.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b5_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x * rpc_x + 9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x * rpc_x +
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x);

        fints_y[i] += dip_y * fss * b5_vals[i] *
                      (18.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 12.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x +
                       rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + 4.0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b6_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z -
                       6.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x -
                       3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b6_vals[i] *
                      (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x * rpc_x - 12.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
                       3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b6_vals[i] * (-6.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b7_vals[i] *
                      (3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x +
                       3.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + 4.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b8_vals[i] * (-rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += faa_y * b1_vals[i] *
                      ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z +
                       3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_y[i] += faa_y * b1_vals[i] * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_y[i] +=
            faa_y * b2_vals[i] *
            (-6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z -
             9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z -
             3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z);

        fints_y[i] +=
            faa_y * b2_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y -
             9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z);

        fints_y[i] += faa_y * b2_vals[i] *
                      (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y);

        fints_y[i] += faa_y * b2_vals[i] *
                      (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - 3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
                       4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_y[i] += faa_y * b2_vals[i] * (-rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_y[i] += faa_y * b3_vals[i] *
                      (9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_y[i] +=
            faa_y * b3_vals[i] *
            (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z +
             6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x +
             6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_y[i] +=
            faa_y * b3_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y + 9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b3_vals[i] *
                      (18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_y[i] += faa_y * b3_vals[i] *
                      ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y +
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_y[i] += faa_y * b3_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x +
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_y[i] += faa_y * b3_vals[i] *
                      ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x +
                       3.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x +
                       12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x +
                       6.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_y[i] += faa_y * b3_vals[i] *
                      (3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
                       4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x -
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y - 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x -
                       9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y -
                       6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
                       (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x -
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - 18.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                       18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y - 9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x -
                       9.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b4_vals[i] *
                      (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x -
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x -
                       12.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_y[i] +=
            faa_y * b4_vals[i] *
            (-18.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
             4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b4_vals[i] * (-6.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_y[i] +=
            faa_y * b5_vals[i] *
            (18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z +
             9.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_y[i] += faa_y * b5_vals[i] *
                      (6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x +
                       18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_y[i] += faa_y * b5_vals[i] *
                      (6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y +
                       9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x +
                       9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b5_vals[i] *
                      (3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z);

        fints_y[i] += faa_y * b5_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_y[i] += faa_y * b5_vals[i] *
                      (9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x +
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + 18.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b5_vals[i] *
                      ((9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x +
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y +
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_y[i] +=
            faa_y * b5_vals[i] *
            (18.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
             12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             4.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
             12.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_y[i] +=
            faa_y * b5_vals[i] *
            (18.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] +=
            faa_y * b6_vals[i] *
            (-9.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z -
             6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] +=
            faa_y * b6_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_y[i] += faa_y * b6_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x -
                       3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] += faa_y * b6_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                       9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_y[i] += faa_y * b6_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 12.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
                       3.0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_y[i] +=
            faa_y * b6_vals[i] *
            (-6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             18.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             4.0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] +=
            faa_y * b7_vals[i] *
            (3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x +
             9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] += faa_y * b7_vals[i] *
                      (6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                       (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] += faa_y * b7_vals[i] *
                      (3.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
                       4.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                       12.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                       3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       6.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] +=
            faa_y * b8_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             4.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] += faa_y * b9_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] *
                      (6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z +
                       9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z + 4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_z[i] += dip_z * fss * b2_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z - 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z -
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z -
                       6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z);

        fints_z[i] += dip_z * fss * b2_vals[i] *
                      (-6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x - 6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y -
                       6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_y - 18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z -
                       9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z);

        fints_z[i] += dip_z * fss * b2_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x - 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y -
                       12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x - 12.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z -
                       4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] *
                      (18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z + 18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x +
                       18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x +
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] *
                      (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z + 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x +
                       6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x + 18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y +
                       6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] *
                      (18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_y +
                       6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z);

        fints_z[i] +=
            dip_z * fss * b3_vals[i] *
            (18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z + 18.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x + 9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x +
             18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y + 9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] *
                      (9.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x + 12.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x +
                       36.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 12.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z +
                       12.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * 12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y;

        fints_z[i] += dip_z * fss * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z - 18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x -
                       18.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z -
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - 6.0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x -
                       18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x - 18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y -
                       18.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y - 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x -
                       18.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_y -
                       18.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_z[i] +=
            dip_z * fss * b4_vals[i] *
            (-6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z -
             9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x - 18.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x - 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y);

        fints_z[i] +=
            dip_z * fss * b4_vals[i] *
            (-18.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y - 18.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x - 9.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x -
             36.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 36.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] *
                      (-4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z - 4.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
                       12.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 36.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
                       12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b5_vals[i] *
                      (18.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z +
                       18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x + 6.0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b5_vals[i] *
                      (6.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + 18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y +
                       18.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x + 18.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b5_vals[i] *
                      (18.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x + 18.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       6.0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 18.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x);

        fints_z[i] +=
            dip_z * fss * b5_vals[i] *
            (9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y + 9.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x + 18.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x +
             36.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 12.0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b5_vals[i] *
                      (12.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 36.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                       36.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 4.0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y +
                       4.0 * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b6_vals[i] *
                      (-6.0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x -
                       18.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x - 6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y -
                       18.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b6_vals[i] *
                      (-6.0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 9.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x -
                       12.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b6_vals[i] *
                      (-12.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 36.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                       12.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 12.0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b7_vals[i] *
                      (6.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       4.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + 12.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                       12.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b8_vals[i] * (-4.0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += faa_z * b1_vals[i] *
                      ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z +
                       3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_z[i] += faa_z * b1_vals[i] * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_z[i] +=
            faa_z * b2_vals[i] *
            (-6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z -
             9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z -
             3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z);

        fints_z[i] +=
            faa_z * b2_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y -
             9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z);

        fints_z[i] += faa_z * b2_vals[i] *
                      (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y);

        fints_z[i] += faa_z * b2_vals[i] *
                      (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - 3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
                       4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_z[i] += faa_z * b2_vals[i] * (-rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_z[i] += faa_z * b3_vals[i] *
                      (9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x +
                       6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_z[i] +=
            faa_z * b3_vals[i] *
            (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z +
             6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x +
             6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_z[i] +=
            faa_z * b3_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y + 9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b3_vals[i] *
                      (18.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_z[i] += faa_z * b3_vals[i] *
                      ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y +
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_z[i] += faa_z * b3_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x +
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_z[i] += faa_z * b3_vals[i] *
                      ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x +
                       3.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x +
                       12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x +
                       6.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_z[i] += faa_z * b3_vals[i] *
                      (3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
                       4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x -
                       18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z -
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
                       6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y - 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x -
                       9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y -
                       6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
                       (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x -
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - 18.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                       18.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x -
                       (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y - 9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x -
                       9.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b4_vals[i] *
                      (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x -
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x -
                       12.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_z[i] +=
            faa_z * b4_vals[i] *
            (-18.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
             4.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b4_vals[i] * (-6.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_z[i] +=
            faa_z * b5_vals[i] *
            (18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z +
             9.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_z[i] += faa_z * b5_vals[i] *
                      (6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x +
                       9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x +
                       18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_z[i] += faa_z * b5_vals[i] *
                      (6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y +
                       9.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x +
                       9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b5_vals[i] *
                      (3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z);

        fints_z[i] += faa_z * b5_vals[i] *
                      ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x +
                       9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_z[i] += faa_z * b5_vals[i] *
                      (9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x +
                       (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + 18.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b5_vals[i] *
                      ((9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x +
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y +
                       (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_z[i] +=
            faa_z * b5_vals[i] *
            (18.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
             12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             4.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
             12.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_z[i] +=
            faa_z * b5_vals[i] *
            (18.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             4.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] +=
            faa_z * b6_vals[i] *
            (-9.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z -
             6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_z[i] +=
            faa_z * b6_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 18.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_z[i] += faa_z * b6_vals[i] *
                      (-9.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                       6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x -
                       3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += faa_z * b6_vals[i] *
                      (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                       (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                       9.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_z[i] += faa_z * b6_vals[i] *
                      (-9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x -
                       (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 12.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
                       3.0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_z[i] +=
            faa_z * b6_vals[i] *
            (-6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             18.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             4.0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] +=
            faa_z * b7_vals[i] *
            (3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x +
             9.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += faa_z * b7_vals[i] *
                      (6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                       (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += faa_z * b7_vals[i] *
                      (3.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
                       4.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                       12.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                       3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                       6.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] +=
            faa_z * b8_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             4.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += faa_z * b9_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;
    }
}

}  // namespace geom_npotrec
