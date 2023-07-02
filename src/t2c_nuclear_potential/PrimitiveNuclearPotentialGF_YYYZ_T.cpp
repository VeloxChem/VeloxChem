#include "PrimitiveNuclearPotentialGF_YYYZ_T.hpp"

#include <cmath>

#include "BoysFunc.hpp"
#include "MathConst.hpp"

namespace npotrec {  // npotrec namespace

auto
compPrimitiveNuclearPotentialGF_YYYZ_T(TDoubleArray&       buffer_xxx,
                                       TDoubleArray&       buffer_xxy,
                                       TDoubleArray&       buffer_xxz,
                                       TDoubleArray&       buffer_xyy,
                                       TDoubleArray&       buffer_xyz,
                                       TDoubleArray&       buffer_xzz,
                                       TDoubleArray&       buffer_yyy,
                                       TDoubleArray&       buffer_yyz,
                                       TDoubleArray&       buffer_yzz,
                                       TDoubleArray&       buffer_zzz,
                                       const double        charge,
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

    // set up Boys function variables

    const CBoysFunc<7> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<8> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto b7_vals = bf_values[7].data();

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

    bf_table.compute<8>(bf_values, bf_args, ket_dim);

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

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_x -
                         rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_xxx[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y +
                         9.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpc_x * rpc_x +
                         3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * 3.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x;

        fints_xxx[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x - 9.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x -
                         9.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] *
                        (-rpa_z * rpa_y * rpa_y * rpa_y * rpc_x * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - 9.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
                         3.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x + 9.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
                         3.0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] *
                        (3.0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + 9.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                         9.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpa_y * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x +
                         rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxx[i] += fss * b5_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x -
                         3.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 9.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                         3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x +
                         rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
                         3.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b7_vals[i] * (-rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_xxy[i] += fss * b0_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y -
                         3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_xxy[i] += fss * b1_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z -
                         (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z - 3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y -
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_x);

        fints_xxy[i] +=
            fss * b1_vals[i] * (-rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x + 9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x * rpc_x +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y +
                         3.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y +
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z +
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] *
                        ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_xxy[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x);

        fints_xxy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z +
             3.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y + 6.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y + rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpc_x * rpc_x +
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y +
                         2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y;

        fints_xxy[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x -
                         9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x * rpc_x -
                         3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-6.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y -
                         3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y -
                         3.0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y -
                         3.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x);

        fints_xxy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y -
             (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-6.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - 3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - 6.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x -
                         rpa_z * rpa_y * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] *
                        (-rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
                         6.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
                         rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y +
                         6.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x +
                         9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x * rpc_x +
                         3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         6.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x +
                         3.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z +
                         3.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 6.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
                         3.0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] *
                        (2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y +
                         6.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] *
                        (6.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x +
                         rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b5_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x -
                         3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b5_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x - 3.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
                         rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
                         3.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] *
                        (-6.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                         2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b6_vals[i] *
                        (3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
                         rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
                         rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * 2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xxy[i] += fss * b7_vals[i] * (-rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z +
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * b0_vals[i] *
                        ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z);

        fints_xxz[i] += fss * b1_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] *
                        (-fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y);

        fints_xxz[i] += fss * b1_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y -
                         (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y);

        fints_xxz[i] += fss * b1_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y -
                         (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y -
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_xxz[i] +=
            fss * b1_vals[i] * (-rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z - rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_z + 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z +
                         fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] *
                        ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z + 3.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * b2_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y +
                         (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y);

        fints_xxz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y + 6.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x +
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpc_x * rpc_x +
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] *
                        (3.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x +
                         rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x -
                         3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x -
                         3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y -
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x -
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x * rpc_x -
                         3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x -
             (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-6.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - 3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - 6.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x -
                         rpa_z * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] *
                        (-rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
                         6.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
                         rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y +
                         3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x +
                         3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] *
                        (3.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] *
                        (fe_0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_xxz[i] += fss * b4_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y +
                         (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
                         6.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
                         rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                         3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] *
                        (3.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 6.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                         rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxz[i] += fss * b5_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x -
                         (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y);

        fints_xxz[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
                         2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                         6.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
                         rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxz[i] += fss * b6_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
                         3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] +=
            fss * b6_vals[i] * (rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b7_vals[i] * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xyy[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyy[i] += fss * b1_vals[i] *
                        (-9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x -
                         3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x -
                         3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] *
                        (-2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_x -
                         rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] *
                        (9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y + 9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x + 9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y +
                         3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_x + 6.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x + 9.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x +
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x +
                         (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] *
                        ((15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] *
                        ((15.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y +
                         6.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x +
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] *
                        (rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y +
                         2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x - 9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y -
                         9.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x -
                         6.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-6.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x -
                         5.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y -
                         9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x - 9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x - 6.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y -
                         (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x -
                         (15.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x -
                         6.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - 6.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y - rpa_z * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x -
                         rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x -
                         rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] *
                        (9.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x +
                         5.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x +
                         9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] *
                        (9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y +
                         6.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + 5.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
                         (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] *
                        ((15.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x +
                         6.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] *
                        (2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x +
                         6.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
                         6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] *
                        (3.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x +
                         rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b5_vals[i] *
                        (-5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x - 9.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x -
                         6.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x - 5.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y -
                         5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] *
                        (-(15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x - 3.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
                         2.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y -
                         6.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
                         2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] *
                        (5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
                         3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
                         rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b7_vals[i] * (-rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x +
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * b0_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z);

        fints_xyz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x -
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_x -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y -
                         (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] *
                        (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x - 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y -
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_x - rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y -
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y +
                         (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xyz[i] += fss * b2_vals[i] *
                        ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x + 3.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x + 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] *
                        (rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x + rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_x +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y +
                         rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] +=
            fss * b2_vals[i] * (rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
                         3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_x -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z -
                         3.0 * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x -
                         (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x -
                         3.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-rpa_z * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
                         3.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x -
                         rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y +
                         3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] *
                        (3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] *
                        ((1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] *
                        (3.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x +
                         3.0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
                         3.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x +
                         rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] *
                        (3.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
                         3.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                         3.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] +=
            fss * b4_vals[i] * (rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b5_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
                         (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x -
                         3.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
                         rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] *
                        (-rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
                         3.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] *
                        (-rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b6_vals[i] *
                        (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
                         rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x +
                         rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] +=
            fss * b6_vals[i] * (rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b7_vals[i] * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x +
                         fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xzz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x;

        fints_xzz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y -
                         fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] *
                        (-fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] *
                        (-2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z - rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_x -
                         rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_x + 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z +
                         3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x +
                         3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y +
                         6.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y +
                         2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x;

        fints_xzz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y -
                         3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z -
                         3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x -
                         3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_x -
                         3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - fe_0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x - 6.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - 6.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x -
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] *
                        (-rpa_z * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_x - rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y -
                         3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] +=
            fss * b3_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x +
                         3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x +
                         3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x +
                         3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] *
                        (fe_0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                         fe_0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                         (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x + 6.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                         3.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x +
                         2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] *
                        (rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 6.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
                         3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                         3.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] +=
            fss * b4_vals[i] * (rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b5_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x -
                         3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] *
                        (-fe_0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
                         (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x - 3.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
                         2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         6.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
                         2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                         rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x +
                         2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_xzz[i] += fss * b7_vals[i] * (-rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_yyy[i] += fss * b0_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] *
                        (-(27.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y -
                         (27.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b1_vals[i] *
                        (-(27.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y - (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y - (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] *
                        (-(27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z -
                         (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z);

        fints_yyy[i] += fss * b1_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y -
                         rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] *
                        (27.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y + (27.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y +
                         (27.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y + 9.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] *
                        (9.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y +
                         (27.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z +
                         (27.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y +
                         (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yyy[i] += fss * b2_vals[i] *
                        ((45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y +
                         (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                         (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] *
                        ((27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z + (45.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z + (45.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] *
                        ((45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z +
                         3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y + 9.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] *
                        (3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] *
                        (-27.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y - 15.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y -
                         9.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y - 15.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y -
                         9.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] *
                        (-27.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - (27.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y -
                         (27.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] *
                        (-9.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y -
                         (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y - (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y -
                         (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] *
                        (-(27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z - (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z - (45.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] *
                        (-(45.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z -
                         (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z - 9.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y -
                         9.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] *
                        (-rpa_z * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 9.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] *
                        (15.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y + 15.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y +
                         (15.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y + 27.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y +
                         15.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] *
                        (9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y + 15.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         9.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y + (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y +
                         (45.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] *
                        ((45.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y + (45.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y +
                         (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z + 9.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] *
                        (3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + 9.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         9.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b5_vals[i] *
                        (-(15.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y - 15.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y -
                         15.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y - (15.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         (45.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y -
                         9.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b6_vals[i] *
                        ((15.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b7_vals[i] * (-rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y +
                         (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z - 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z -
                         3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] *
                        (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y -
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z -
                         rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] *
                        (9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y +
                         9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_z + 6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y +
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y +
                         3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z +
                         9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y +
                         3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z +
                         3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y +
                         (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y +
                         3.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] *
                        (rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y +
                         3.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y +
                         rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y -
                         9.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y -
                         6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-5.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y -
                         9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z -
                         9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-(15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z -
                         6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y -
                         (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-(5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y -
                         (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y -
                         6.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y -
                         6.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_z * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y -
                         rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y -
                         2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        (9.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y + 5.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y +
                         6.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y + 5.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        (9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y +
                         9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y +
                         5.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y +
                         fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y +
                         (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        ((15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z +
                         3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y +
                         5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        ((15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y +
                         6.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] *
                        (rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] +=
            fss * b4_vals[i] * (rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] *
                        (-5.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] *
                        (-fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y - 5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y -
                         (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y -
                         2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b6_vals[i] *
                        (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y +
                         rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
                         rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b6_vals[i] * 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yyz[i] += fss * b7_vals[i] * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += fss * b0_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                         (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z);

        fints_yzz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y -
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z -
                         fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z);

        fints_yzz[i] += fss * b1_vals[i] *
                        (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z -
                         3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_yzz[i] +=
            fss * b1_vals[i] *
            (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y -
             (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_yzz[i] += fss * b1_vals[i] *
                        (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z - 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y -
                         2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y -
                         rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y +
                         3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z +
                         3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y +
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z +
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y +
                         fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y +
                         3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        (9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z +
                         (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] *
                        (3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z +
                         3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] +=
            fss * b2_vals[i] *
            ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z +
             3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 6.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y +
                         2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y;

        fints_yzz[i] += fss * b3_vals[i] *
                        (-9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y -
                         3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y -
                         9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z -
                         3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - fe_0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z - (27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y -
                         3.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-6.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z -
             6.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-6.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y -
                         rpa_z * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y -
                         3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] *
                        (-6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y -
                         2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y +
                         6.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y +
                         (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] *
                        (3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y +
                         9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b4_vals[i] *
                        (3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y +
                         6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y +
                         3.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y +
                         (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y +
                         (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y + 9.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b4_vals[i] *
                        ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z + 6.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         3.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y +
                         2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] *
                        (rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y +
                         3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y +
                         3.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] +=
            fss * b4_vals[i] * (rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] *
                        (-3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y -
                         6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] *
                        (-fe_0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         3.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z -
                         3.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                         rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] *
                        (-6.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
                         3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b6_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
                         rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b6_vals[i] * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_yzz[i] += fss * b7_vals[i] * (-rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b0_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_zzz[i] += fss * b0_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                         rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z -
                         (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z -
                         (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y - (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] *
                        (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y -
                         (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y);

        fints_zzz[i] += fss * b1_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z -
                         rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z + (27.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z +
                         3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] *
                        ((27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] *
                        ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (27.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y +
                         (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                         (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] *
                        (3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y + 9.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y +
                         3.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y +
                         3.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y - (27.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y -
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y -
                         (27.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y -
                         3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y - (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z -
                         9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y -
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(27.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y -
                         (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y -
                         (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 9.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y -
                         9.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_z * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z -
                         rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] *
                        (-3.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y -
                         3.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] *
                        ((9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z +
                         (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y + (27.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] *
                        (9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y +
                         (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y +
                         (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] *
                        ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z +
                         (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (27.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y +
                         9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] *
                        ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y +
                         9.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y +
                         3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] *
                        (9.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 9.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y +
                         rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z + rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b5_vals[i] *
                        (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y -
                         9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z -
                         (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b5_vals[i] *
                        (-(9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y -
                         (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b5_vals[i] *
                        (-3.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                         9.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
                         3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b6_vals[i] *
                        (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
                         rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
                         3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b7_vals[i] * (-rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);
    }
}

}  // namespace npotrec
