#include "PrimitiveNuclearPotentialGG_XXYZ_T.hpp"

#include <cmath>

#include "BoysFunc.hpp"
#include "MathConst.hpp"

namespace npotrec {  // npotrec namespace

auto
compPrimitiveNuclearPotentialGG_XXYZ_T(TDoubleArray&       buffer_xxxx,
                                       TDoubleArray&       buffer_xxxy,
                                       TDoubleArray&       buffer_xxxz,
                                       TDoubleArray&       buffer_xxyy,
                                       TDoubleArray&       buffer_xxyz,
                                       TDoubleArray&       buffer_xxzz,
                                       TDoubleArray&       buffer_xyyy,
                                       TDoubleArray&       buffer_xyyz,
                                       TDoubleArray&       buffer_xyzz,
                                       TDoubleArray&       buffer_xzzz,
                                       TDoubleArray&       buffer_yyyy,
                                       TDoubleArray&       buffer_yyyz,
                                       TDoubleArray&       buffer_yyzz,
                                       TDoubleArray&       buffer_yzzz,
                                       TDoubleArray&       buffer_zzzz,
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

    auto fints_xxxx = buffer_xxxx.data();

    auto fints_xxxy = buffer_xxxy.data();

    auto fints_xxxz = buffer_xxxz.data();

    auto fints_xxyy = buffer_xxyy.data();

    auto fints_xxyz = buffer_xxyz.data();

    auto fints_xxzz = buffer_xxzz.data();

    auto fints_xyyy = buffer_xyyy.data();

    auto fints_xyyz = buffer_xyyz.data();

    auto fints_xyzz = buffer_xyzz.data();

    auto fints_xzzz = buffer_xzzz.data();

    auto fints_yyyy = buffer_yyyy.data();

    auto fints_yyyz = buffer_yyyz.data();

    auto fints_yyzz = buffer_yyzz.data();

    auto fints_yzzz = buffer_yzzz.data();

    auto fints_zzzz = buffer_zzzz.data();

    // set up Boys function variables

    const CBoysFunc<8> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<9> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto b7_vals = bf_values[7].data();

    auto b8_vals = bf_values[8].data();

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

    bf_table.compute<9>(bf_values, bf_args, ket_dim);

#pragma omp simd aligned(fints_xxxx,     \
                             fints_xxxy, \
                             fints_xxxz, \
                             fints_xxyy, \
                             fints_xxyz, \
                             fints_xxzz, \
                             fints_xyyy, \
                             fints_xyyz, \
                             fints_xyzz, \
                             fints_xzzz, \
                             fints_yyyy, \
                             fints_yyyz, \
                             fints_yyzz, \
                             fints_yzzz, \
                             fints_zzzz, \
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

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxxx[i] += fss * b0_vals[i] *
                         (4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x + 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x);

        fints_xxxx[i] += fss * b0_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xxxx[i] += fss * b1_vals[i] *
                         (-18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x - 4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x -
                          6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x -
                          6.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxxx[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x - 4.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y -
             3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * rpc_y -
             4.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxx[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z - 12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x -
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x);

        fints_xxxx[i] += fss * b1_vals[i] *
                         (-15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x - 9.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x -
                          6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y);

        fints_xxxx[i] += fss * b1_vals[i] *
                         (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxxx[i] +=
            fss * b1_vals[i] *
            (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z - 2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_x -
             4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_y -
             rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxx[i] += fss * b2_vals[i] *
                         (24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x + 18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x +
                          18.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x + 18.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x +
                          4.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxxx[i] +=
            fss * b2_vals[i] *
            (6.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * rpc_y +
             18.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + 4.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z +
             6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xxxx[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z + 6.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z +
                          4.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxx[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y + 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x +
                          15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x +
                          30.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x);

        fints_xxxx[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_xxxx[i] += fss * b2_vals[i] *
                         (15.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + 9.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y +
                          12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z);

        fints_xxxx[i] += fss * b2_vals[i] *
                         (15.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + 9.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z +
                          6.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxx[i] += fss * b2_vals[i] *
                         ((45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y +
                          8.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b2_vals[i] *
            (6.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x + rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x +
             2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxxx[i] +=
            fss * b2_vals[i] *
            (4.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x - 10.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x - 20.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x -
                          18.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-24.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - 18.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x -
                          6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x -
                          18.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x - 24.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x -
                          18.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-18.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x -
                          18.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - 4.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y -
                          6.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y - 6.0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y -
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - 15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-(45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y -
                          15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y -
                          30.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y - (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x -
                          6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - 15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-30.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z -
                          (45.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - 12.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y -
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - 15.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x -
                          9.0 * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxxx[i] += fss * b3_vals[i] *
                         (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y -
                          12.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b3_vals[i] *
            (-8.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             8.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b3_vals[i] *
            (-rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             4.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         (10.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x + 20.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x +
                          (15.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          24.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x + 10.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x + 20.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          18.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x + 24.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x +
                          10.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x + 20.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          18.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x + 24.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          18.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          18.0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + 6.0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x +
                          (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + 15.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x +
                          (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x +
                          15.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         ((45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y +
                          15.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                          30.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (45.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x +
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxxx[i] +=
            fss * b4_vals[i] *
            (8.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x +
             12.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] *
                         (4.0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          12.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          8.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b5_vals[i] *
            (-(15.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x - 10.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x -
             20.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             10.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b5_vals[i] *
            (-20.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             24.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x - 10.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b5_vals[i] *
                         (-20.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 18.0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x -
                          (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (45.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x -
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xxxx[i] += fss * b5_vals[i] *
                         (-15.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (45.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x -
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b5_vals[i] *
            (-8.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             8.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b5_vals[i] *
                         (-6.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
                          12.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          4.0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b6_vals[i] *
            ((15.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             (15.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + 10.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             20.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + (15.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b6_vals[i] *
                         ((45.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b6_vals[i] *
            (4.0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
             8.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             6.0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] +=
            fss * b7_vals[i] *
            (-(15.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             4.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b8_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxxy[i] +=
            fss * b0_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y);

        fints_xxxy[i] += fss * b0_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xxxy[i] +=
            fss * b0_vals[i] * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xxxy[i] +=
            fss * b1_vals[i] *
            (-9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x -
             3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x);

        fints_xxxy[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y -
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y);

        fints_xxxy[i] +=
            fss * b1_vals[i] *
            (-fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxxy[i] +=
            fss * b1_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z -
             3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y);

        fints_xxxy[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xxxy[i] += fss * b1_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xxxy[i] += fss * b1_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxxy[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x -
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xxxy[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z -
             2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxxy[i] +=
            fss * b1_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_y - rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxy[i] += fss * b2_vals[i] *
                         (9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x +
                          9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x +
             9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y +
                          9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            (9.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y +
                          9.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y);

        fints_xxxy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x +
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xxxy[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z +
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x + 6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            (3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxxy[i] +=
            fss * b2_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y +
             rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x -
             6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x -
             9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-5.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x -
                          9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-9.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x -
                          9.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - 9.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y -
                          9.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x -
             (9.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y -
                          6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y -
                          9.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z -
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - 9.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxy[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x -
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x - rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x -
             3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxxy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x +
                          9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_z * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x +
             5.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            ((9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x +
             6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x + 9.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x + 9.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x +
             5.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x +
             (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y +
             3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b4_vals[i] *
                         ((15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x +
                          5.0 * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_xxxy[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y +
                          6.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxxy[i] += fss * b4_vals[i] *
                         ((15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y +
                          9.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b4_vals[i] *
                         ((5.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 6.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            (3.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             6.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             3.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b4_vals[i] *
            (3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b5_vals[i] *
                         (-5.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_z * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x - 5.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x -
                          9.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x - 5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x -
             5.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b5_vals[i] *
                         (-9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x - 5.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-6.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x -
             9.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-5.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b5_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b5_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x -
                          5.0 * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-6.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-6.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b5_vals[i] *
            (-6.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b6_vals[i] *
            (5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
             5.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             fe_0 * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b6_vals[i] *
            (5.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + 5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b6_vals[i] *
            ((15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (5.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b6_vals[i] *
            (3.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b6_vals[i] *
            (6.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b7_vals[i] *
            (-5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] +=
            fss * b7_vals[i] *
            (-rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b8_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxxz[i] +=
            fss * b0_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z);

        fints_xxxz[i] += fss * b0_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xxxz[i] +=
            fss * b0_vals[i] * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_xxxz[i] +=
            fss * b1_vals[i] *
            (-9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x -
             3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_xxxz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z -
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxxz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z -
             fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z);

        fints_xxxz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z);

        fints_xxxz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y);

        fints_xxxz[i] += fss * b1_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x);

        fints_xxxz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxxz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x -
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xxxz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y -
             2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxz[i] +=
            fss * b1_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxz[i] += fss * b2_vals[i] *
                         (9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x +
                          9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x +
             9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z + 9.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y +
             (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b2_vals[i] *
                         (9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z +
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x +
                          fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            ((9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z +
                          9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z);

        fints_xxxz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xxxz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y +
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x + 6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            (3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxxz[i] +=
            fss * b2_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z +
             rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x -
             6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x -
             9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-5.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x -
                          9.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-9.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y -
             9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 9.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-9.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x -
                          9.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y -
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - 9.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxxz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x -
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x - rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x -
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x + 5.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x +
                          9.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x + 5.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x +
                          9.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + 9.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x +
             9.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (5.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x +
             (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y +
             (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x);

        fints_xxxz[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_xxxz[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x +
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + 5.0 * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y +
                          6.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxxz[i] += fss * b4_vals[i] *
                         ((15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y +
                          9.0 * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxxz[i] += fss * b4_vals[i] *
                         ((5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + 6.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (3.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             6.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             3.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b4_vals[i] *
            (3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b5_vals[i] *
                         (-5.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x -
                          5.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x -
                          5.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b5_vals[i] *
                         (-6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x -
                          5.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-5.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x -
             9.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-5.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b5_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x -
                          (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b5_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x -
                          5.0 * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-6.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-6.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b5_vals[i] *
            (-6.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             3.0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b6_vals[i] *
                         (5.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + 5.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b6_vals[i] *
            (5.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + 5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b6_vals[i] *
            ((15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b6_vals[i] *
            (3.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b6_vals[i] *
            (6.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b7_vals[i] *
            (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] +=
            fss * b7_vals[i] *
            (-rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxyy[i] +=
            fss * b0_vals[i] *
            (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x +
             fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxyy[i] += fss * b0_vals[i] *
                         (fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x +
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x);

        fints_xxyy[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y - 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x - fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y - fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x - fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y -
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y -
                          2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x - fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxyy[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * rpc_y -
             2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z - fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z -
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x - fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y -
                          fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_xxyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xxyy[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_x -
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y -
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_y);

        fints_xxyy[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y +
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y +
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x + fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x +
                          4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y + 4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x + fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y +
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y + fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y +
             fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * rpc_y +
             4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z + fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z +
                          2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxyy[i] +=
            fss * b2_vals[i] *
            (2.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y +
                          2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z);

        fints_xxyy[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z +
                          4.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y +
                          fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xxyy[i] += fss * b2_vals[i] *
                         ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxyy[i] +=
            fss * b2_vals[i] *
            (4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x +
             4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x +
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y + rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b2_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyy[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x -
             6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x -
             4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x -
                          6.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y - fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x -
             2.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x - fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y - 6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x - fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y -
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
                          4.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 4.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x - fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y - fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
             fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y -
                          9.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-9.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y -
                          fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z -
                          6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y -
                          6.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y -
                          3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y -
                          8.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b3_vals[i] *
                         (-2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-4.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y - rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxyy[i] +=
            fss * b3_vals[i] *
            (-rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
             2.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x +
             6.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x +
             6.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b4_vals[i] *
            (2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x + 2.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_z * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y +
             6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y + 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x +
                          4.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y + 4.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x + fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y +
                          6.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y +
                          9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x +
                          3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 6.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxyy[i] += fss * b4_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z +
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 4.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
             8.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (2.0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          8.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] *
                         (2.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          4.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b4_vals[i] *
            (4.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpb_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x -
             6.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-fe_0 * rpa_z * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x -
                          6.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-6.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x -
                          6.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - fe_0 * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y -
             (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x -
                          9.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-2.0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 4.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-8.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
                          2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b5_vals[i] *
                         (-2.0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b6_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + 6.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b6_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b6_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             4.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b6_vals[i] *
            (rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyy[i] +=
            fss * b7_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b7_vals[i] *
                         (-2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
                          2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyy[i] += fss * b8_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxyz[i] +=
            fss * b0_vals[i] *
            (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxyz[i] += fss * b0_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x +
                          fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x);

        fints_xxyz[i] += fss * b0_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxyz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x +
                          (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y - 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y);

        fints_xxyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z - 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y - fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y - fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x);

        fints_xxyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * rpc_y -
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z - fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x -
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z -
                          2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x - fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y - fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z -
                          fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_xxyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxyz[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y -
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z - rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_y -
             rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y +
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y + 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y +
                          fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y +
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y + fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x +
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z +
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x + fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x +
             fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y +
             3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x +
             2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y +
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z + fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x +
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y +
                          2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x +
                          2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x + fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z +
             fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y + fe_0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y +
                          fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y +
             rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y + rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x +
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] +=
            fss * b2_vals[i] *
            (rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z + rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y -
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y -
             2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z -
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x -
                          fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x - fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x -
                          fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x -
                          2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y - 2.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x -
                          3.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_x -
                          (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpb_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x -
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
             4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y -
             rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
             rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b3_vals[i] *
            (-rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y +
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x +
                          fe_0 * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y +
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             fe_0 * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y +
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y + 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x +
             fe_0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x +
             3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x +
             3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x + fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z +
                          fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (fe_0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + fe_0 * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + fe_0 * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] *
                         ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
             4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
             rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             2.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b4_vals[i] *
            (rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
             rpb_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-fe_0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x -
                          fe_0 * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-2.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - fe_0 * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x -
             3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
                          fe_0 * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y -
                          3.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             4.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             4.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            (fe_0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + fe_0 * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b6_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                          fe_0 * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            ((1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b7_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] +=
            fss * b7_vals[i] *
            (-2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxyz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxzz[i] +=
            fss * b0_vals[i] *
            (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x +
             fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xxzz[i] += fss * b0_vals[i] *
                         (fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x +
                          2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x);

        fints_xxzz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z - 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x - fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z - fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x - fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z -
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x);

        fints_xxzz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z -
                          2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x - fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z - fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y -
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y -
                          4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z -
                          fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y);

        fints_xxzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xxzz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y - 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_x -
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z -
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxzz[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (4.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x + 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z +
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z +
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x + fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z +
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x +
                          4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y +
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y +
                          fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * rpc_y +
             4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z + 4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x +
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x +
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z +
                          3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z + fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z +
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxzz[i] +=
            fss * b2_vals[i] *
            (2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y);

        fints_xxzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y + 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x +
                          6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x + 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z +
                          6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z +
                          4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y +
                          fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxzz[i] +=
            fss * b2_vals[i] *
            (4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x +
             4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x +
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z + rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxzz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z -
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x -
             6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x -
             4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y - fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z - 4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x -
                          6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x - fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z -
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z - 6.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x -
             2.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z - fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z -
                          4.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - 4.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x - fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y -
                          fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y - fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
             fe_0 * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x -
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z -
                          6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z -
                          9.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y -
                          6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y -
                          6.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y -
                          3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y -
                          8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b3_vals[i] *
                         (-2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-4.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
             rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxzz[i] +=
            fss * b3_vals[i] *
            (-rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y -
             2.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x +
             6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y +
             2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x +
             6.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x +
                          6.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z +
             6.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
                          fe_0 * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z + 4.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y +
             4.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x + fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y +
                          fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y +
                          6.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          fe_0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y + fe_0 * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z);

        fints_xxzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z +
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y +
                          9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x +
                          3.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + 6.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxzz[i] += fss * b4_vals[i] *
                         ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y +
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 4.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x + rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             8.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (2.0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                          8.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] *
                         (2.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x +
                          2.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
                          2.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x +
                          4.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          4.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b4_vals[i] *
            (4.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x - 6.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x -
                          fe_0 * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b5_vals[i] *
            (-2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-6.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - 2.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x -
                          6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
                          3.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x -
                          9.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             4.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-8.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
                          2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b5_vals[i] *
                         (-2.0 * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                          2.0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b6_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxzz[i] += fss * b6_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          3.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b6_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                          2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             4.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b6_vals[i] *
            (rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxzz[i] +=
            fss * b7_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b7_vals[i] *
                         (-2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
                          2.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxzz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xyyy[i] +=
            fss * b0_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y);

        fints_xyyy[i] += fss * b0_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xyyy[i] +=
            fss * b0_vals[i] * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_xyyy[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x);

        fints_xyyy[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xyyy[i] +=
            fss * b1_vals[i] *
            (-fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyyy[i] += fss * b1_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y);

        fints_xyyy[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xyyy[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xyyy[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_xyyy[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xyyy[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z -
             2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x);

        fints_xyyy[i] +=
            fss * b1_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_y - rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x + 9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y +
                          (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y +
                          fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y +
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x +
                          9.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y +
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xyyy[i] += fss * b2_vals[i] *
                         ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b2_vals[i] *
            (3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y + rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x +
             2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y;

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y -
             9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x - 9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y -
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x - 9.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y -
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(27.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x -
                          9.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(27.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x -
             6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y +
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + 9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x +
                          (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            ((9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x + 9.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y +
                          3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x +
                          9.0 * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y);

        fints_xyyy[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xyyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y +
                          (27.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b4_vals[i] *
                         ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y +
             3.0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             6.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
             3.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b4_vals[i] *
            (rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
             3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b4_vals[i] * rpb_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x -
             fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y -
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             9.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] += fss * b5_vals[i] *
                         (-6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y -
             (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
             (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x -
             (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x -
                          (27.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b5_vals[i] *
                         (-9.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
                          3.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 6.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
             3.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             3.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             6.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b6_vals[i] *
            (fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
             (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b6_vals[i] *
            (rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + 6.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
             3.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b6_vals[i] * 3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyyy[i] +=
            fss * b7_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
             rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] +=
            fss * b7_vals[i] *
            (-3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyy[i] += fss * b8_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyyz[i] += fss * b0_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x + fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyyz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x +
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y);

        fints_xyyz[i] += fss * b0_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x +
                          (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xyyz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y - fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x - fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x);

        fints_xyyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z - fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y -
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x - fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y - fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * rpc_y - fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z -
                          fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xyyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z - fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y -
                          fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z -
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y -
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y - fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x -
                          fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y -
                          fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xyyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xyyz[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_x -
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x -
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z - rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y + fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y + fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z + fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x + fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z +
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y + fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y +
                          fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y +
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z +
                          fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x + fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z +
                          fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y + fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y +
                          fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z +
                          fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z +
                          fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x +
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y +
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x +
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z +
                          2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                          fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x +
                          fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y + fe_0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_xyyz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x +
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x +
             rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x +
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyz[i] +=
            fss * b2_vals[i] *
            (rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z + rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y - fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y - fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x - fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y -
             3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y -
             2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x - fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y -
                          fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y -
             fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y -
             fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x - fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x -
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
             fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x - fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y -
                          fe_0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z -
                          fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - fe_0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y - fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z -
                          fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z -
                          2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z -
                          fe_0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y - fe_0 * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y -
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y -
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xyyz[i] += fss * b3_vals[i] *
                         (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x -
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y - 2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y -
             rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x -
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y -
             2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyyz[i] +=
            fss * b3_vals[i] *
            (-rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x - rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] * (-rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y + fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y +
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y +
             fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x + fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y +
                          fe_0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                          fe_0 * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + fe_0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x + fe_0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
             fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z +
                          fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b4_vals[i] *
                         (fe_0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x +
                          fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x +
                          3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x +
                          fe_0 * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y +
             2.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
             2.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (4.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
             rpb_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x -
             fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x -
             fe_0 * rpa_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x -
             fe_0 * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
                          fe_0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x -
             fe_0 * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-fe_0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
                          fe_0 * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - fe_0 * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             4.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             4.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
             2.0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b5_vals[i] *
            (-rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b6_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y +
                          fe_0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] += fss * b6_vals[i] *
                         (3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
                          fe_0 * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyyz[i] +=
            fss * b6_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x +
             rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b6_vals[i] *
            (rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyyz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b7_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyyz[i] +=
            fss * b7_vals[i] *
            (-2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyyz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyzz[i] += fss * b0_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x + fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xyzz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x);

        fints_xyzz[i] += fss * b0_vals[i] *
                         (fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x +
                          (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xyzz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z - fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y - fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x);

        fints_xyzz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y - fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y - fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y -
                          fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x - fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xyzz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * rpc_y -
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x - fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x -
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z -
                          fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y -
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y -
                          2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y - fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z -
                          fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z -
                          fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xyzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_x -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xyzz[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_x -
             2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x -
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xyzz[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z + 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y +
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y + fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y +
             2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y + fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y + fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z + fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y +
             fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x +
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z +
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z +
                          fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x +
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y + fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z + fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x + fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y +
                          fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y +
                          fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z + fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x +
                          2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z +
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x + fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x + fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y +
                          2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x + fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z +
                          fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_xyzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b2_vals[i] *
                         (4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z +
             rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x +
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyzz[i] +=
            fss * b2_vals[i] *
            (rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y + rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y - fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x - fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x - fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y -
             2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y -
                          fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y -
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x - fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y -
                          fe_0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y -
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y - fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x -
                          fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x - fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z -
                          fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y - fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x -
                          fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - fe_0 * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z -
                          fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y -
                          fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y -
                          2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x -
                          fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y - fe_0 * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] *
                         (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x - 4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x -
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - 4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x -
             rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x -
             2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b3_vals[i] *
            (-rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b3_vals[i] * (-rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y + fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x +
             2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y +
                          fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x + fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x +
                          fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x + fe_0 * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x +
                          fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                          fe_0 * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y + 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x +
                          2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y +
                          3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x +
                          fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + fe_0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y +
                          fe_0 * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                          fe_0 * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
                          fe_0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x +
             fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + fe_0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x +
                          fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x +
                          fe_0 * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y +
             2.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
                          4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            (rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y +
             2.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
             rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b4_vals[i] *
                         (4.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] +=
            fss * b4_vals[i] *
            (rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y +
             rpb_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x -
             fe_0 * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x - fe_0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x -
                          3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
                          fe_0 * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_x -
             2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             4.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             4.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b5_vals[i] *
            (-rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyzz[i] += fss * b6_vals[i] *
                         (3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            (fe_0 * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x +
             (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x +
             rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            (rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            (rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + 4.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyzz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x +
             rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b7_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyzz[i] +=
            fss * b7_vals[i] *
            (-2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyzz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xzzz[i] +=
            fss * b0_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z);

        fints_xzzz[i] += fss * b0_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xzzz[i] +=
            fss * b0_vals[i] * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_xzzz[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_xzzz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_xzzz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x -
             fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z);

        fints_xzzz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z);

        fints_xzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y);

        fints_xzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x);

        fints_xzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_xzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xzzz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y -
             2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x);

        fints_xzzz[i] +=
            fss * b1_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z +
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
             9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x + fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z +
             (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x + fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y +
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y +
                          9.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xzzz[i] += fss * b2_vals[i] *
                         ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y;

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z -
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y -
             3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x -
                          9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x - 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x - 9.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y -
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x -
                          6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x -
                          9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(27.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_x -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x -
             6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z +
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y +
                          9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z +
                          6.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z +
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            ((9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_xzzz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x + 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y +
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_xzzz[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z);

        fints_xzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + 9.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y +
                          (27.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] *
                         ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x * rpc_y + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y +
             3.0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x +
             6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z +
             3.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b4_vals[i] *
            (rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y +
             3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] += fss * b4_vals[i] * rpb_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x -
             fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z -
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] += fss * b5_vals[i] *
                         (-6.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y -
             (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                          (27.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b5_vals[i] *
                         (-9.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
                          3.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 6.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x -
             3.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             6.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             6.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x +
             6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b6_vals[i] *
            (fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzzz[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
             (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x +
             rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
             3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b6_vals[i] *
            (rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + 6.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
             3.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b6_vals[i] * 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xzzz[i] +=
            fss * b7_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] +=
            fss * b7_vals[i] *
            (-3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzzz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yyyy[i] +=
            fss * b0_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y +
             2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y);

        fints_yyyy[i] += fss * b0_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_yyyy[i] += fss * b1_vals[i] *
                         (-6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y - 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_yyyy[i] +=
            fss * b1_vals[i] *
            (-4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x - 9.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y -
             2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * rpc_y -
             3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z -
                          2.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y);

        fints_yyyy[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x -
                          6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_yyyy[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z -
                          fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b1_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y - (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y -
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yyyy[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_x -
             4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y - rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_y -
             rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y +
                          18.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 4.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x +
                          12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yyyy[i] +=
            fss * b2_vals[i] *
            (9.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y + 2.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y +
             2.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * rpc_y +
             6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_yyyy[i] +=
            fss * b2_vals[i] *
            (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z +
             2.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z +
             4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + 2.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z +
                          (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x);

        fints_yyyy[i] += fss * b2_vals[i] *
                         ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y +
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y + 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y +
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (9.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y + fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z +
                          3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x +
                          6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyyy[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_yyyy[i] +=
            fss * b2_vals[i] *
            ((9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y +
             8.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x +
             6.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y + rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyyy[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x +
             4.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y - 6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x -
                          24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - 18.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x -
                          12.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y - 2.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x -
                          6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y -
                          3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y -
                          18.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 4.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x -
                          12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyy[i] +=
            fss * b3_vals[i] *
            (-9.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y - 2.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y -
             2.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x -
                          15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y - 12.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y -
                          6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-(5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z - 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x -
                          12.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-(15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - 6.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y -
                          3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - 9.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y -
                          fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y - (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y -
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y - 12.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y -
                          4.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
                          8.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] +=
            fss * b3_vals[i] *
            (-6.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y - rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             8.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
             6.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b3_vals[i] *
                         (-2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
                          4.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x + 10.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x +
                          5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y +
                          12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
                          3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x +
                          24.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (18.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + 12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y +
                          5.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y + 2.0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x +
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x +
                          15.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_yyyy[i] += fss * b4_vals[i] *
                         ((15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + 12.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y +
                          6.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b4_vals[i] *
                         ((5.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyyy[i] +=
            fss * b4_vals[i] *
            (8.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y +
             6.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
             12.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] *
                         (4.0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          12.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
                          4.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
                          4.0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          8.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + rpb_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y -
             10.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x - 12.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] +=
            fss * b5_vals[i] *
            (-5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y -
             6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
             24.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x - 10.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             5.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b5_vals[i] *
                         (-12.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
                          9.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y -
                          (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b5_vals[i] *
                         (-(5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          6.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x -
                          5.0 * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y -
                          2.0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             8.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y -
             6.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             8.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] +=
            fss * b5_vals[i] *
            (-rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             12.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             4.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
             4.0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] +=
            fss * b6_vals[i] *
            (5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y +
             3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
             10.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] +=
            fss * b6_vals[i] *
            (12.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
             5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y +
             (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b6_vals[i] *
                         ((5.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
                          4.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] +=
            fss * b6_vals[i] *
            (4.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             8.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y +
             6.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] +=
            fss * b7_vals[i] *
            (-5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y -
             rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyyy[i] += fss * b7_vals[i] * (-4.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyy[i] += fss * b8_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yyyz[i] +=
            fss * b0_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y);

        fints_yyyz[i] += fss * b0_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y);

        fints_yyyz[i] += fss * b0_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y +
                          (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yyyz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y;

        fints_yyyz[i] +=
            fss * b1_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yyyz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_yyyz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y);

        fints_yyyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_yyyz[i] += fss * b1_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y);

        fints_yyyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_yyyz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yyyz[i] += fss * b1_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_yyyz[i] += fss * b1_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yyyz[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y - rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z -
             rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x +
             (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x +
             fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z + 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_yyyz[i] += fss * b2_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b2_vals[i] *
            (3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyyz[i] += fss * b2_vals[i] * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y;

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y -
                          9.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x -
                          9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y -
             3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y -
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x -
                          9.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_yyyz[i] += fss * b3_vals[i] *
                         (-(27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpb_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x -
             6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x -
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y +
             6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x + 9.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            ((9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yyyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                          9.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y +
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y +
             3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y);

        fints_yyyz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + 6.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yyyz[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x +
                          (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
             3.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             6.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y +
             3.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
             3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y +
             3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-fe_0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] += fss * b5_vals[i] *
                         (-9.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x -
                          6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y -
             (9.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y -
             3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b5_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b5_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             6.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyyz[i] +=
            fss * b5_vals[i] *
            (-rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
             3.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b6_vals[i] *
            (6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
             3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
             (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b6_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y +
                          rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             3.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b6_vals[i] *
            (3.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x +
             6.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y +
             3.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b6_vals[i] * 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yyyz[i] +=
            fss * b7_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] +=
            fss * b7_vals[i] *
            (-2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyyz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yyzz[i] +=
            fss * b0_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y +
             fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yyzz[i] += fss * b0_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y);

        fints_yyzz[i] += fss * b0_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y +
                          fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_yyzz[i] += fss * b0_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x - fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y - fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z -
                          fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yyzz[i] +=
            fss * b1_vals[i] *
            (-2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y -
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z - fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y -
                          fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y -
                          fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y -
                          fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z - fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y - fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x -
                          fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z -
                          fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x -
                          2.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z -
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_yyzz[i] +=
            fss * b1_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y +
                          2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y + fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z +
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y +
                          fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y +
                          4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x + fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x +
                          3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y +
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z +
                          3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z + fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x + fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z +
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z + 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y +
                          fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y + fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z +
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y + 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z +
                          2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (4.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x +
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y +
                          3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y);

        fints_yyzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
             fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z + fe_0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y + fe_0 * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y + fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x +
                          4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yyzz[i] +=
            fss * b2_vals[i] *
            (4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y +
             rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z + rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b2_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y +
                          2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b2_vals[i] * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y;

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x -
                          fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x - fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y -
                          fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y -
                          fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z - fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y -
                          2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
                          fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x -
                          2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x -
                          6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x - fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z -
                          fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z -
                          2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z - fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b3_vals[i] *
            (-fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x - 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x -
                          2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z -
                          2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y -
                          fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y - fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y - fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y -
                          fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y - fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z -
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z -
                          fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_yyzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z - fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y -
             fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - fe_0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z - 3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x -
                          6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - 2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z -
                          2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y -
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y -
                          (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x -
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y -
                          2.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x -
                          4.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b3_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y -
             rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - 4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] +=
            fss * b3_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y -
             rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b3_vals[i] *
                         (-2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y -
                          2.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x +
                          fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x +
             6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y +
                          2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y + fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x +
             2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z +
             2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y +
                          fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x + 4.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x +
                          6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y +
                          fe_0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y + fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y +
             2.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            (fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y +
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z +
                          3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + 9.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z +
                          fe_0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (fe_0 * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + fe_0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x +
                          fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_z + (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y +
             4.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            (4.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
             rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 8.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y +
                          2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y +
                          2.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (8.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                          2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y +
                          2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] *
                         (2.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
                          2.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
                          4.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
                          4.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x +
                          4.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b4_vals[i] *
            (rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y +
             rpb_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - fe_0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y -
                          fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x -
                          fe_0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x -
                          fe_0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x -
                          fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             6.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y -
                          2.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x -
                          2.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - fe_0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y -
             3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                          4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
             4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             4.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
             4.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b5_vals[i] *
            (-rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 8.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b5_vals[i] *
                         (-2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
                          2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                          2.0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b6_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
                          3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] += fss * b6_vals[i] *
                         (fe_0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
                          fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y +
             (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] +=
            fss * b6_vals[i] *
            (rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyzz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             4.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x +
             4.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyzz[i] +=
            fss * b6_vals[i] *
            (4.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x +
             rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] +=
            fss * b7_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y -
             rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b7_vals[i] *
                         (-2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x -
                          2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x -
                          2.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyzz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yzzz[i] +=
            fss * b0_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y);

        fints_yzzz[i] += fss * b0_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z);

        fints_yzzz[i] += fss * b0_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z +
                          (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yzzz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y;

        fints_yzzz[i] +=
            fss * b1_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzzz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y - fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_yzzz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_yzzz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y);

        fints_yzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_yzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z -
                          (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y);

        fints_yzzz[i] += fss * b1_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_yzzz[i] += fss * b1_vals[i] *
                         (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_yzzz[i] += fss * b1_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_yzzz[i] += fss * b1_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0);

        fints_yzzz[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z - rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y -
             rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_y - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x +
                          fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x + 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x +
             3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z +
             (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y +
             3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            (fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z +
                          (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x +
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * fe_0 +
                          6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z +
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x +
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b2_vals[i] *
            (3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y;

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x -
             3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z -
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x -
             9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x - 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x -
                          (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
             3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y -
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y -
                          3.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - 9.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x -
                          (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_yzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y -
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x -
             6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x -
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] +=
            fss * b3_vals[i] *
            (-3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x +
             3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x +
             6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y +
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z +
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] += fss * b4_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x + 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x +
                          3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y +
             3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x +
                          6.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x * rpc_x +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y * rpc_y +
                          (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yzzz[i] += fss * b4_vals[i] *
                         ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x +
                          (3.0 / 16.0) * fe_0 * fe_0 * fe_0 * fe_0 + 6.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
                          2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
             6.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y +
             3.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x +
             6.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z +
             3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x +
             3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y +
             3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x -
             fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x -
             (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 6.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x -
             (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z -
             3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yzzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x - 6.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x -
                          6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z -
             (9.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y -
             (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x -
             (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b5_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x -
                          (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y -
                          (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y -
                          3.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b5_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y * rpc_y -
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             3.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             6.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 3.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             6.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             6.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             3.0 * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b6_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzzz[i] +=
            fss * b6_vals[i] *
            (6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x +
             (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
             (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
             (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b6_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
                          (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x +
                          rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b6_vals[i] *
            (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x +
             3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
             2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
             3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b6_vals[i] *
            (rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + 6.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x +
             2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y +
             3.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b6_vals[i] * 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yzzz[i] +=
            fss * b7_vals[i] *
            (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y -
             (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
             rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] +=
            fss * b7_vals[i] *
            (-2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x -
             3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzzz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_zzzz[i] +=
            fss * b0_vals[i] *
            (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z +
             2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z);

        fints_zzzz[i] += fss * b0_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z +
                          (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z +
                          rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_zzzz[i] += fss * b1_vals[i] *
                         (-6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x - 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z - 2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z -
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_zzzz[i] +=
            fss * b1_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y -
             4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x - 9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z -
             2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_zzzz[i] += fss * b1_vals[i] *
                         (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z -
                          2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z);

        fints_zzzz[i] += fss * b1_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x -
                          6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z);

        fints_zzzz[i] += fss * b1_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z -
                          2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z - 3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y -
                          fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_zzzz[i] += fss * b1_vals[i] *
                         (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y -
                          (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_zzzz[i] +=
            fss * b1_vals[i] *
            (-2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x -
             4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y -
             rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x +
                          6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z +
                          3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x + 2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z +
                          6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y +
                          3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_zzzz[i] +=
            fss * b2_vals[i] *
            (2.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y +
             18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x +
             12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_zzzz[i] +=
            fss * b2_vals[i] *
            (9.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z + 2.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z +
             2.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z +
             4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + 2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y +
                          (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x +
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x);

        fints_zzzz[i] += fss * b2_vals[i] *
                         ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y +
                          3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y + 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x +
                          (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z +
                          6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x +
                          9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z + 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x +
                          6.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                          (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b2_vals[i] *
                         (2.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y +
                          (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_zzzz[i] +=
            fss * b2_vals[i] *
            ((9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y +
             8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x +
             6.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_zzzz[i] +=
            fss * b2_vals[i] *
            (2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y +
             2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x +
             4.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z + rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x - 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x -
                          3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z - 6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x -
                          2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z - 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x -
                          12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x -
                          6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y -
                          3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y -
                          24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x - 12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z -
                          5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
                          3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-2.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z - 2.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x -
                          18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 4.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x -
                          12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] +=
            fss * b3_vals[i] *
            (-9.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y - 2.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y -
             2.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x -
                          (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x -
                          6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z -
                          12.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x -
                          (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-(15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z -
                          12.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-(15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - 6.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y -
                          3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - 9.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y -
                          fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-(3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z -
                          (9.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y - 12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
                          4.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z -
                          4.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x -
                          8.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_zzzz[i] +=
            fss * b3_vals[i] *
            (-6.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x -
             8.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x -
             6.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_zzzz[i] += fss * b3_vals[i] *
                         (-2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x -
                          4.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x + 6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x +
                          (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x + 6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x +
                          3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y + 6.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x +
                          24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x + 10.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x +
                          5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (12.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z +
                          9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z +
                          24.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x + 12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y +
                          5.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          3.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (2.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y + 2.0 * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y +
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x +
                          6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x +
                          5.0 * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x +
                          15.0 * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_zzzz[i] += fss * b4_vals[i] *
                         ((15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + 12.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y +
                          6.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y +
                          (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzzz[i] += fss * b4_vals[i] *
                         ((5.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z * rpc_y +
                          (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z * rpc_y +
                          (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_zzzz[i] +=
            fss * b4_vals[i] *
            (8.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_z * rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z +
             6.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
             12.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x +
             4.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b4_vals[i] *
                         (4.0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          12.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x +
                          4.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z +
                          4.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
                          8.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzzz[i] +=
            fss * b4_vals[i] *
            (6.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z -
             6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x - 6.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x -
             2.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] +=
            fss * b5_vals[i] *
            (-3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             10.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x - 12.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x -
             2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzzz[i] +=
            fss * b5_vals[i] *
            (-5.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z -
             24.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x - 10.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             5.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b5_vals[i] *
                         (-12.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
                          9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y -
                          (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b5_vals[i] *
                         (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x -
                          (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                          6.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_zzzz[i] += fss * b5_vals[i] *
                         (-3.0 * fe_0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x -
                          5.0 * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_z * rpc_y -
                          2.0 * rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zzzz[i] +=
            fss * b5_vals[i] *
            (-4.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             8.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             6.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
             8.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zzzz[i] +=
            fss * b5_vals[i] *
            (-rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z - 6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             12.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x -
             4.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             4.0 * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzzz[i] +=
            fss * b6_vals[i] *
            (3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             5.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z +
             10.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzzz[i] +=
            fss * b6_vals[i] *
            (12.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzzz[i] += fss * b6_vals[i] *
                         ((5.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
                          2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x +
                          4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x +
                          2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zzzz[i] +=
            fss * b6_vals[i] *
            (4.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x +
             8.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y +
             6.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzzz[i] +=
            fss * b7_vals[i] *
            (-5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y -
             rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
             2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzzz[i] += fss * b7_vals[i] * (-4.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzzz[i] += fss * b8_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;
    }
}

}  // namespace npotrec
