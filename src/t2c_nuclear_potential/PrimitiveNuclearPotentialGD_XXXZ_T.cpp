#include "PrimitiveNuclearPotentialGD_XXXZ_T.hpp"

#include <cmath>

#include "BoysFunc.hpp"
#include "MathConst.hpp"

namespace npotrec {  // npotrec namespace

auto
compPrimitiveNuclearPotentialGD_XXXZ_T(TDoubleArray&       buffer_xx,
                                       TDoubleArray&       buffer_xy,
                                       TDoubleArray&       buffer_xz,
                                       TDoubleArray&       buffer_yy,
                                       TDoubleArray&       buffer_yz,
                                       TDoubleArray&       buffer_zz,
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

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

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

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

#pragma omp simd aligned(fints_xx, fints_xy, fints_xz, fints_yy, fints_yz, fints_zz, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
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

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xx[i] += fss * b0_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x +
                        (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                        (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xx[i] += fss * b0_vals[i] * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x;

        fints_xx[i] += fss * b1_vals[i] *
                       (-9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x -
                        3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x -
                        (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x);

        fints_xx[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z -
                        3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z -
                        (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xx[i] += fss * b1_vals[i] *
                       (-3.0 * fe_0 * fe_0 * rpa_z * rpb_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z -
                        (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - 3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xx[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xx[i] += fss * b2_vals[i] *
                       (9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + 9.0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x +
                        (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x + 6.0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_xx[i] += fss * b2_vals[i] *
                       (9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z +
                        3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x +
                        (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z);

        fints_xx[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xx[i] +=
            fss * b2_vals[i] *
            (3.0 * fe_0 * fe_0 * rpb_x * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x +
             6.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x + rpa_z * rpa_x * rpa_x * rpa_x * rpc_x * rpc_x);

        fints_xx[i] += fss * b2_vals[i] * (3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xx[i] +=
            fss * b3_vals[i] *
            (-9.0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - 6.0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x - 5.0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x -
             9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - 9.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * b3_vals[i] *
                       (-(9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - 6.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x -
                        (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xx[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x -
                        6.0 * rpa_z * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x -
                        rpa_z * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * b3_vals[i] *
                       (-3.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 6.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x -
                        rpa_x * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xx[i] +=
            fss * b4_vals[i] *
            (5.0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + 9.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + 6.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x +
             5.0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xx[i] += fss * b4_vals[i] *
                       (3.0 * rpa_z * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x +
                        6.0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x +
                        rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * b5_vals[i] *
                       (-5.0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - rpa_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x -
                        3.0 * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * b6_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * b0_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y +
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x);

        fints_xy[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_xy[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y -
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xy[i] += fss * b1_vals[i] *
                       (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x -
                        rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpc_x - rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpc_y -
                        rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z);

        fints_xy[i] += fss * b2_vals[i] *
                       ((9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y +
                        (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x);

        fints_xy[i] += fss * b2_vals[i] *
                       (3.0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xy[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y +
                        (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y +
                        (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xy[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y +
                        3.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x + 3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x +
                        3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * b2_vals[i] *
                       (rpa_z * rpa_x * rpa_x * rpa_x * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x +
                        rpa_x * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xy[i] += fss * b3_vals[i] *
                       (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x -
                        (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xy[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x -
                        3.0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xy[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y -
                        (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y -
                        3.0 * rpa_z * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xy[i] += fss * b3_vals[i] *
                       (-3.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x -
                        rpa_z * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x -
                        3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xy[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * b4_vals[i] *
                       (3.0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x +
                        3.0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x +
                        3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * b4_vals[i] *
                       ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x +
                        rpa_z * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x +
                        3.0 * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xy[i] += fss * b4_vals[i] *
                       (3.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x +
                        rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xy[i] += fss * b5_vals[i] *
                       (-3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x -
                        3.0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x -
                        rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xz[i] += fss * b0_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z +
                        (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xz[i] += fss * b0_vals[i] *
                       ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x);

        fints_xz[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_xz[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x -
                        (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x);

        fints_xz[i] +=
            fss * b1_vals[i] *
            (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x -
             (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x);

        fints_xz[i] += fss * b1_vals[i] *
                       (-rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpc_x - rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpc_z -
                        rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z);

        fints_xz[i] += fss * b2_vals[i] *
                       ((9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z +
                        (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x);

        fints_xz[i] += fss * b2_vals[i] *
                       (3.0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xz[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x);

        fints_xz[i] +=
            fss * b2_vals[i] *
            ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x +
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xz[i] +=
            fss * b2_vals[i] *
            ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z +
             (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x +
             3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * b2_vals[i] *
                       (rpa_z * rpa_x * rpa_x * rpa_x * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x +
                        rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xz[i] += fss * b3_vals[i] *
                       (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x -
                        (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_xz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x -
                        (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_xz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x -
                        3.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x -
                        (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xz[i] +=
            fss * b3_vals[i] *
            (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z -
             (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xz[i] += fss * b3_vals[i] *
                       (-3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_z * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x -
                        3.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * b3_vals[i] *
                       (-rpa_z * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x -
                        3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x -
                        rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * b4_vals[i] *
                       (3.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * b4_vals[i] *
                       ((1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x +
                        (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z +
                        (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_xz[i] += fss * b4_vals[i] *
                       (3.0 * rpa_z * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x +
                        rpa_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x +
                        3.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xz[i] += fss * b5_vals[i] *
                       (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x -
                        rpa_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x -
                        rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xz[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_yy[i] += fss * b0_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x +
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y);

        fints_yy[i] += fss * b1_vals[i] *
                       (-3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x);

        fints_yy[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z -
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_yy[i] += fss * b1_vals[i] *
                       (-3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpc_y -
                        rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yy[i] += fss * b2_vals[i] *
                       (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x +
                        3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x);

        fints_yy[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x +
                        (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z);

        fints_yy[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x +
             (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_yy[i] += fss * b2_vals[i] *
                       (3.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x + 6.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x +
                        rpa_z * rpa_x * rpa_x * rpa_x * rpc_y * rpc_y + 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x +
                        2.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_yy[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x -
                        3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x -
                        (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x);

        fints_yy[i] += fss * b3_vals[i] *
                       (-3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x -
                        3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yy[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x -
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x -
                        6.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yy[i] += fss * b3_vals[i] *
                       (-3.0 * rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x -
                        3.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 6.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x -
                        rpa_x * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * b4_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x +
                        3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yy[i] += fss * b4_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x +
                        (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x +
                        2.0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yy[i] += fss * b4_vals[i] *
                       (6.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x +
                        rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_yy[i] += fss * b5_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x -
                        rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x -
                        2.0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yz[i] += fss * b0_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y +
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y);

        fints_yz[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_yz[i] += fss * b1_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y -
                        (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y -
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_yz[i] += fss * b1_vals[i] *
                       (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x -
                        rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpc_y - rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpc_z -
                        rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_yz[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x);

        fints_yz[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_yz[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x +
                        (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x +
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yz[i] +=
            fss * b2_vals[i] *
            ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x +
             3.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x + 3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_yz[i] += fss * b2_vals[i] *
                       (3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpa_x * rpc_z * rpc_y +
                        3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y +
                        rpa_x * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_yz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_yz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_yz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x -
                        (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x -
                        (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_yz[i] += fss * b3_vals[i] *
                       (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x -
                        3.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 3.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x -
                        3.0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * b3_vals[i] *
                       (-rpa_z * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x -
                        3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x -
                        rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * b4_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y +
                        (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yz[i] += fss * b4_vals[i] *
                       ((1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x +
                        (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x +
                        3.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yz[i] += fss * b4_vals[i] *
                       (rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x +
                        3.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x +
                        3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * b4_vals[i] * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_yz[i] += fss * b5_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x -
                        rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x -
                        rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_yz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_zz[i] += fss * b0_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x +
                        fe_0 * rpa_x * rpa_x * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_zz[i] += fss * b0_vals[i] * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z;

        fints_zz[i] += fss * b1_vals[i] *
                       (-3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x -
                        (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x);

        fints_zz[i] +=
            fss * b1_vals[i] *
            (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z - 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x - fe_0 * rpa_x * rpa_x * rpa_x * rpb_z -
             (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_zz[i] += fss * b1_vals[i] *
                       (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - 3.0 * fe_0 * fe_0 * rpa_x * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z -
                        (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_zz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z - rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_zz[i] += fss * b2_vals[i] *
                       (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x +
                        3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x);

        fints_zz[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z +
                        3.0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z +
                        3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_zz[i] += fss * b2_vals[i] *
                       ((9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z +
                        (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                        (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_zz[i] += fss * b2_vals[i] *
                       ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + 3.0 * fe_0 * fe_0 * rpb_z * rpc_x +
                        (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_zz[i] += fss * b2_vals[i] *
                       (6.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpa_x * rpc_z * rpc_z +
                        3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_zz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x -
                        3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x -
                        (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x);

        fints_zz[i] += fss * b3_vals[i] *
                       (-3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - 3.0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x -
                        (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z -
                        (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_zz[i] +=
            fss * b3_vals[i] *
            (-3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x -
             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_zz[i] += fss * b3_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x -
                        6.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x -
                        rpa_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_zz[i] += fss * b3_vals[i] *
                       (-3.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x -
                        rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * b4_vals[i] *
                       ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x +
                        (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z +
                        3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * b4_vals[i] *
                       (fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x +
                        (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x +
                        3.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_zz[i] += fss * b4_vals[i] *
                       (2.0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x +
                        3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_zz[i] += fss * b5_vals[i] *
                       (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x -
                        rpa_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x -
                        2.0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_zz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;
    }
}

}  // namespace npotrec
