#include "PrimitiveNuclearPotentialGeom020FG_YYY_YYZZ.hpp"

#include <cmath>

#include "BoysFunc.hpp"
#include "MathConst.hpp"

namespace geom_npotrec { // geom_npotrec namespace

auto
compPrimitiveNuclearPotentialGeom020FG_YYY_YYZZ(      TDoubleArray& buffer_xx,
                                                      TDoubleArray& buffer_xy,
                                                      TDoubleArray& buffer_xz,
                                                      TDoubleArray& buffer_yy,
                                                      TDoubleArray& buffer_yz,
                                                      TDoubleArray& buffer_zz,
                              const T2Tensor& quadrupole,
                              const TPoint3D& point,
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

    // set up quadrupole components

    const auto qpol_xx = quadrupole[0];

    const auto qpol_xy = quadrupole[1];

    const auto qpol_xz = quadrupole[2];

    const auto qpol_yy = quadrupole[3];

    const auto qpol_yz = quadrupole[4];

    const auto qpol_zz = quadrupole[5];

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

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

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
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

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_xx = qpol_xx * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_x;

        const auto faa_xy = qpol_xy * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_y;

        const auto faa_xz = qpol_xz * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpcz;

        const auto faa_yy = qpol_yy * fss * 4.0 * fxi_0 * fxi_0 * rpc_y * rpc_y;

        const auto faa_yz = qpol_yz * fss * 4.0 * fxi_0 * fxi_0 * rpc_y * rpcz;

        const auto faa_zz = qpol_zz * fss * 4.0 * fxi_0 * fxi_0 * rpc_z * rpcz;

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto faa = -2.0 * fxi_0 * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += (faa_xx * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_xx[i] += (faa_xx * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xx[i] += (faa_xx * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xx[i] += (faa_xx * b3_vals[i] + faa * b2_vals[i]) * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_xx[i] += (faa_xx * b3_vals[i] + faa * b2_vals[i]) * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_xx[i] += (faa_xx * b3_vals[i] + faa * b2_vals[i]) * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y);

        fints_xx[i] += (faa_xx * b3_vals[i] + faa * b2_vals[i]) * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_xx[i] += (faa_xx * b3_vals[i] + faa * b2_vals[i]) * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xx[i] += (faa_xx * b3_vals[i] + faa * b2_vals[i]) * (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xx[i] += (faa_xx * b4_vals[i] + faa * b3_vals[i]) * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b6_vals[i] + faa * b5_vals[i]) * (18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b6_vals[i] + faa * b5_vals[i]) * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b6_vals[i] + faa * b5_vals[i]) * (6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_xx[i] += (faa_xx * b6_vals[i] + faa * b5_vals[i]) * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_xx[i] += (faa_xx * b6_vals[i] + faa * b5_vals[i]) * (5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b6_vals[i] + faa * b5_vals[i]) * (6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b7_vals[i] + faa * b6_vals[i]) * (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b7_vals[i] + faa * b6_vals[i]) * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b7_vals[i] + faa * b6_vals[i]) * (-6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b8_vals[i] + faa * b7_vals[i]) * (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += (faa_xx * b9_vals[i] + faa * b8_vals[i]) * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xy[i] += qpol_xy * faa_x * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xy[i] += qpol_xy * faa_x * b2_vals[i] * (3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_xy[i] += qpol_xy * faa_x * b3_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y - 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z);

        fints_xy[i] += qpol_xy * faa_x * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z - 6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y - fe_0 * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b3_vals[i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z - 12.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y - 9.0 * fe_0 * fe_0 * rpa_y * rpb_y - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_xy[i] += qpol_xy * faa_x * b3_vals[i] * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 12.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_xy[i] += qpol_xy * faa_x * b3_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + 36.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + 6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * (fe_0 * rpa_y * rpa_y * rpa_y * rpc_y + 24.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z + 12.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y + 15.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y + 18.0 * fe_0 * fe_0 * rpa_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * (15.0 * fe_0 * fe_0 * rpb_z * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + 12.0 * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * ((15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 + 12.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 18.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 24.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b4_vals[i] * (9.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b5_vals[i] * (-36.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - 18.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - 6.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - 24.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 30.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - 15.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b5_vals[i] * (-12.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - 4.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b5_vals[i] * (-(15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 15.0 * fe_0 * fe_0 * rpc_y * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xy[i] += qpol_xy * faa_x * b5_vals[i] * (-36.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 12.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 18.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b5_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - 6.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 8.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + 6.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 30.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + 12.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + 4.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b6_vals[i] * (15.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (5.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + 24.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b6_vals[i] * (18.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 9.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 16.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b7_vals[i] * (-15.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (5.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 10.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 8.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += qpol_xy * faa_x * b8_vals[i] * 5.0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_xy[i] += faa_xy * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_xy[i] += faa_xy * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xy[i] += faa_xy * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xy[i] += faa_xy * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_xy[i] += faa_xy * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_xy[i] += faa_xy * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y);

        fints_xy[i] += faa_xy * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_xy[i] += faa_xy * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xy[i] += faa_xy * b3_vals[i] * (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_xy[i] += faa_xy * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xy[i] += faa_xy * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_xy[i] += faa_xy * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_xy[i] += faa_xy * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_xy[i] += faa_xy * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b4_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_xy[i] += faa_xy * b4_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_xy[i] += faa_xy * b4_vals[i] * ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xy[i] += faa_xy * b4_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xy[i] += faa_xy * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b5_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b6_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b6_vals[i] * (6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_xy[i] += faa_xy * b6_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_xy[i] += faa_xy * b6_vals[i] * (5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b6_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b7_vals[i] * (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b7_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b7_vals[i] * (-6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b8_vals[i] * (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xy[i] += faa_xy * b9_vals[i] * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + 3.0 * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_xz[i] += qpol_xz * faa_x * b2_vals[i] * 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y;

        fints_xz[i] += qpol_xz * faa_x * b3_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b3_vals[i] * (-6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpc_z - 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - 9.0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xz[i] += qpol_xz * faa_x * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - 6.0 * fe_0 * fe_0 * rpb_z * rpb_y - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z - 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b3_vals[i] * (-4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_xz[i] += qpol_xz * faa_x * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + 18.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + 18.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b4_vals[i] * (6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpc_z + 12.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b4_vals[i] * (3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpc_z + 3.0 * fe_0 * fe_0 * rpb_z * rpb_y + 15.0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b4_vals[i] * (6.0 * fe_0 * fe_0 * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y + 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b4_vals[i] * (2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - 18.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b5_vals[i] * (-10.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b5_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpc_z * rpc_y - 12.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b5_vals[i] * (-12.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b6_vals[i] * (6.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 12.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 4.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b7_vals[i] * (-10.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 4.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += qpol_xz * faa_x * b8_vals[i] * 2.0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_xz[i] += faa_xz * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_xz[i] += faa_xz * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xz[i] += faa_xz * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xz[i] += faa_xz * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_xz[i] += faa_xz * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_xz[i] += faa_xz * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y);

        fints_xz[i] += faa_xz * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_xz[i] += faa_xz * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xz[i] += faa_xz * b3_vals[i] * (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_xz[i] += faa_xz * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xz[i] += faa_xz * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_xz[i] += faa_xz * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_xz[i] += faa_xz * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_xz[i] += faa_xz * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b4_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_xz[i] += faa_xz * b4_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_xz[i] += faa_xz * b4_vals[i] * ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xz[i] += faa_xz * b4_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_xz[i] += faa_xz * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b5_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b6_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b6_vals[i] * (6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_xz[i] += faa_xz * b6_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_xz[i] += faa_xz * b6_vals[i] * (5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b6_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b7_vals[i] * (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b7_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b7_vals[i] * (-6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b8_vals[i] * (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += faa_xz * b9_vals[i] * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b2_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_z + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y + 12.0 * fe_0 * rpb_z * rpb_z * rpb_y);

        fints_yy[i] += qpol_yy * fss * b2_vals[i] * (9.0 * fe_0 * fe_0 * rpa_y + 6.0 * fe_0 * fe_0 * rpb_y + 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 12.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z);

        fints_yy[i] += qpol_yy * fss * b3_vals[i] * (-36.0 * fe_0 * rpa_y * rpb_z * rpc_z - 18.0 * fe_0 * rpa_y * rpb_z * rpb_z - 18.0 * fe_0 * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_y);

        fints_yy[i] += qpol_yy * fss * b3_vals[i] * (-9.0 * fe_0 * rpa_y * rpa_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y - 24.0 * fe_0 * rpb_z * rpb_y * rpc_z - 12.0 * fe_0 * rpb_z * rpb_z * rpb_y - 30.0 * fe_0 * rpb_z * rpb_z * rpc_y);

        fints_yy[i] += qpol_yy * fss * b3_vals[i] * (-3.0 * fe_0 * rpb_y * rpb_y * rpc_y - 18.0 * fe_0 * fe_0 * rpa_y - 12.0 * fe_0 * fe_0 * rpb_y - 15.0 * fe_0 * fe_0 * rpc_y - 12.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yy[i] += qpol_yy * fss * b3_vals[i] * (-36.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 24.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 18.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - 6.0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b4_vals[i] * (36.0 * fe_0 * rpa_y * rpb_z * rpc_z + 18.0 * fe_0 * rpa_y * rpb_y * rpc_y + 18.0 * fe_0 * rpa_y * rpc_z * rpc_z + 18.0 * fe_0 * rpa_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b4_vals[i] * (24.0 * fe_0 * rpb_z * rpb_y * rpc_z + 60.0 * fe_0 * rpb_z * rpc_z * rpc_y + 30.0 * fe_0 * rpb_z * rpb_z * rpc_y + 12.0 * fe_0 * rpb_y * rpc_z * rpc_z + 12.0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b4_vals[i] * (3.0 * fe_0 * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y + 6.0 * fe_0 * fe_0 * rpb_y + 30.0 * fe_0 * fe_0 * rpc_y + 72.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yy[i] += qpol_yy * fss * b4_vals[i] * (36.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + 6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z + 36.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + 12.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + 2.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_yy[i] += qpol_yy * fss * b4_vals[i] * (12.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 24.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpc_z * rpc_z - 18.0 * fe_0 * rpa_y * rpc_y * rpc_y - 60.0 * fe_0 * rpb_z * rpc_z * rpc_y - 12.0 * fe_0 * rpb_y * rpc_z * rpc_z - 12.0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b5_vals[i] * (-30.0 * fe_0 * rpc_z * rpc_z * rpc_y - 10.0 * fe_0 * rpc_y * rpc_y * rpc_y - 15.0 * fe_0 * fe_0 * rpc_y - 72.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 36.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yy[i] += qpol_yy * fss * b5_vals[i] * (-18.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - 48.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 20.0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - 6.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yy[i] += qpol_yy * fss * b6_vals[i] * (30.0 * fe_0 * rpc_z * rpc_z * rpc_y + 10.0 * fe_0 * rpc_y * rpc_y * rpc_y + 36.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 40.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 24.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * fss * b7_vals[i] * (-20.0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b2_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + 3.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yy[i] += qpol_yy * faa_y * b2_vals[i] * (9.0 * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 4.0) * fe_0 * fe_0 * fe_0);

        fints_yy[i] += qpol_yy * faa_y * b2_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yy[i] += qpol_yy * faa_y * b3_vals[i] * (-36.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y - 36.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - 6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - 18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z);

        fints_yy[i] += qpol_yy * faa_y * b3_vals[i] * (-9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z - 12.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y - 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y - 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b3_vals[i] * (-6.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z - 24.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y - 18.0 * fe_0 * fe_0 * rpa_y * rpb_y - 18.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b3_vals[i] * (-9.0 * fe_0 * fe_0 * rpa_y * rpa_y - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpb_z - 12.0 * fe_0 * fe_0 * rpb_y * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_yy[i] += qpol_yy * faa_y * b3_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * fe_0 - 12.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 24.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 8.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_yy[i] += qpol_yy * faa_y * b3_vals[i] * (-4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (36.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + 72.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + 36.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + 18.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + 18.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + 18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + 12.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y + 48.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z + 24.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y + 30.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y + 36.0 * fe_0 * fe_0 * rpa_y * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (30.0 * fe_0 * fe_0 * rpb_z * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z + 24.0 * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (15.0 * fe_0 * fe_0 * rpc_y * rpc_y + (45.0 / 4.0) * fe_0 * fe_0 * fe_0 + 24.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 36.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 48.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b4_vals[i] * (18.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z + 8.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + 6.0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b5_vals[i] * (-72.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - 18.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - 36.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - 12.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b5_vals[i] * (-9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - 48.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 60.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - 30.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b5_vals[i] * (-24.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - 8.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - 18.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b5_vals[i] * (-15.0 * fe_0 * fe_0 * rpb_z * rpc_z - 12.0 * fe_0 * fe_0 * rpb_y * rpc_y - 15.0 * fe_0 * fe_0 * rpc_z * rpc_z - 30.0 * fe_0 * fe_0 * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * fe_0);

        fints_yy[i] += qpol_yy * faa_y * b5_vals[i] * (-72.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 24.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 36.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 24.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b5_vals[i] * (-4.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - 12.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 16.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b6_vals[i] * (36.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + 12.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 60.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + 24.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + 8.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b6_vals[i] * (30.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + 5.0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + 15.0 * fe_0 * fe_0 * rpc_y * rpc_y + 48.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b6_vals[i] * (36.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 18.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 32.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 10.0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 6.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b7_vals[i] * (-30.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - 5.0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 24.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 20.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 16.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += qpol_yy * faa_y * b8_vals[i] * 10.0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yy[i] += (faa_yy * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_yy[i] += (faa_yy * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yy[i] += (faa_yy * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yy[i] += (faa_yy * b3_vals[i] + faa * b2_vals[i]) * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yy[i] += (faa_yy * b3_vals[i] + faa * b2_vals[i]) * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_yy[i] += (faa_yy * b3_vals[i] + faa * b2_vals[i]) * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y);

        fints_yy[i] += (faa_yy * b3_vals[i] + faa * b2_vals[i]) * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_yy[i] += (faa_yy * b3_vals[i] + faa * b2_vals[i]) * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += (faa_yy * b3_vals[i] + faa * b2_vals[i]) * (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yy[i] += (faa_yy * b4_vals[i] + faa * b3_vals[i]) * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b6_vals[i] + faa * b5_vals[i]) * (18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b6_vals[i] + faa * b5_vals[i]) * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b6_vals[i] + faa * b5_vals[i]) * (6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_yy[i] += (faa_yy * b6_vals[i] + faa * b5_vals[i]) * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_yy[i] += (faa_yy * b6_vals[i] + faa * b5_vals[i]) * (5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b6_vals[i] + faa * b5_vals[i]) * (6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b7_vals[i] + faa * b6_vals[i]) * (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b7_vals[i] + faa * b6_vals[i]) * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b7_vals[i] + faa * b6_vals[i]) * (-6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b8_vals[i] + faa * b7_vals[i]) * (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += (faa_yy * b9_vals[i] + faa * b8_vals[i]) * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * fss * b2_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yz[i] += qpol_yz * fss * b2_vals[i] * 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y;

        fints_yz[i] += qpol_yz * fss * b3_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y - 36.0 * fe_0 * rpa_y * rpb_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z - 9.0 * fe_0 * rpa_y * rpa_y * rpc_z);

        fints_yz[i] += qpol_yz * fss * b3_vals[i] * (-24.0 * fe_0 * rpb_z * rpb_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpb_y - 3.0 * fe_0 * rpb_y * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z);

        fints_yz[i] += qpol_yz * fss * b3_vals[i] * (-12.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y - 24.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z - 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y - 4.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z);

        fints_yz[i] += qpol_yz * fss * b4_vals[i] * (36.0 * fe_0 * rpa_y * rpb_z * rpc_y + 18.0 * fe_0 * rpa_y * rpb_y * rpc_z + 36.0 * fe_0 * rpa_y * rpc_z * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpc_z + 24.0 * fe_0 * rpb_z * rpb_y * rpc_y);

        fints_yz[i] += qpol_yz * fss * b4_vals[i] * (30.0 * fe_0 * rpb_z * rpc_y * rpc_y + 24.0 * fe_0 * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpb_y * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z + 15.0 * fe_0 * fe_0 * rpc_z);

        fints_yz[i] += qpol_yz * fss * b4_vals[i] * (36.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 12.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y + 18.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 24.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_yz[i] += qpol_yz * fss * b4_vals[i] * 6.0 * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y;

        fints_yz[i] += qpol_yz * fss * b5_vals[i] * (-36.0 * fe_0 * rpa_y * rpc_z * rpc_y - 30.0 * fe_0 * rpb_z * rpc_y * rpc_y - 24.0 * fe_0 * rpb_y * rpc_z * rpc_y - 30.0 * fe_0 * rpc_z * rpc_y * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z);

        fints_yz[i] += qpol_yz * fss * b5_vals[i] * (-24.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y - 36.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - 18.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - 16.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 6.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * fss * b6_vals[i] * (30.0 * fe_0 * rpc_z * rpc_y * rpc_y + 24.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 10.0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 16.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * fss * b7_vals[i] * (-10.0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yz[i] += qpol_yz * faa_z * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_yz[i] += qpol_yz * faa_z * b2_vals[i] * (3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yz[i] += qpol_yz * faa_z * b3_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y - 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z);

        fints_yz[i] += qpol_yz * faa_z * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z - 6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y - fe_0 * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b3_vals[i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z - 12.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y - 9.0 * fe_0 * fe_0 * rpa_y * rpb_y - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_yz[i] += qpol_yz * faa_z * b3_vals[i] * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 12.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_yz[i] += qpol_yz * faa_z * b3_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + 36.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + 6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * (fe_0 * rpa_y * rpa_y * rpa_y * rpc_y + 24.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z + 12.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y + 15.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y + 18.0 * fe_0 * fe_0 * rpa_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * (15.0 * fe_0 * fe_0 * rpb_z * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + 12.0 * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * ((15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 + 12.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 18.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 24.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b4_vals[i] * (9.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b5_vals[i] * (-36.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - 18.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - 6.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - 24.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 30.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - 15.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b5_vals[i] * (-12.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - 4.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b5_vals[i] * (-(15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 15.0 * fe_0 * fe_0 * rpc_y * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_yz[i] += qpol_yz * faa_z * b5_vals[i] * (-36.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 12.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 18.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b5_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - 6.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 8.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + 6.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 30.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + 12.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + 4.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b6_vals[i] * (15.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (5.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + 24.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b6_vals[i] * (18.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 9.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 16.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b7_vals[i] * (-15.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (5.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 10.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 8.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_z * b8_vals[i] * 5.0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yz[i] += qpol_yz * faa_y * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + 3.0 * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yz[i] += qpol_yz * faa_y * b2_vals[i] * 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y;

        fints_yz[i] += qpol_yz * faa_y * b3_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b3_vals[i] * (-6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpc_z - 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - 9.0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_yz[i] += qpol_yz * faa_y * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - 6.0 * fe_0 * fe_0 * rpb_z * rpb_y - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z - 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b3_vals[i] * (-4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yz[i] += qpol_yz * faa_y * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + 18.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + 18.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b4_vals[i] * (6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpc_z + 12.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b4_vals[i] * (3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpc_z + 3.0 * fe_0 * fe_0 * rpb_z * rpb_y + 15.0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b4_vals[i] * (6.0 * fe_0 * fe_0 * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y + 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b4_vals[i] * (2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - 18.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b5_vals[i] * (-10.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b5_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpc_z * rpc_y - 12.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b5_vals[i] * (-12.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b6_vals[i] * (6.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 12.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 4.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b7_vals[i] * (-10.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 4.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += qpol_yz * faa_y * b8_vals[i] * 2.0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yz[i] += faa_yz * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_yz[i] += faa_yz * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yz[i] += faa_yz * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yz[i] += faa_yz * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yz[i] += faa_yz * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_yz[i] += faa_yz * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y);

        fints_yz[i] += faa_yz * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_yz[i] += faa_yz * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yz[i] += faa_yz * b3_vals[i] * (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_yz[i] += faa_yz * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yz[i] += faa_yz * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yz[i] += faa_yz * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_yz[i] += faa_yz * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_yz[i] += faa_yz * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b4_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_yz[i] += faa_yz * b4_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_yz[i] += faa_yz * b4_vals[i] * ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += faa_yz * b4_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yz[i] += faa_yz * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b5_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b6_vals[i] * (18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b6_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b6_vals[i] * (6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_yz[i] += faa_yz * b6_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += faa_yz * b6_vals[i] * (5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b6_vals[i] * (6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b7_vals[i] * (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b7_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b7_vals[i] * (-6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b8_vals[i] * (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += faa_yz * b9_vals[i] * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_y * rpb_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y + 3.0 * fe_0 * fe_0 * rpb_y);

        fints_zz[i] += qpol_zz * fss * b2_vals[i] * 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y;

        fints_zz[i] += qpol_zz * fss * b3_vals[i] * (-18.0 * fe_0 * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_y - 9.0 * fe_0 * rpa_y * rpa_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y);

        fints_zz[i] += qpol_zz * fss * b3_vals[i] * (-3.0 * fe_0 * rpb_y * rpb_y * rpc_y - 9.0 * fe_0 * fe_0 * rpa_y - 6.0 * fe_0 * fe_0 * rpb_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b3_vals[i] * (-4.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b4_vals[i] * (18.0 * fe_0 * rpa_y * rpb_y * rpc_y + 18.0 * fe_0 * rpa_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpc_y + 12.0 * fe_0 * rpb_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y + 3.0 * fe_0 * fe_0 * rpb_y + 15.0 * fe_0 * fe_0 * rpc_y + 6.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 12.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b4_vals[i] * 2.0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y;

        fints_zz[i] += qpol_zz * fss * b5_vals[i] * (-18.0 * fe_0 * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_y * rpc_y * rpc_y - 10.0 * fe_0 * rpc_y * rpc_y * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_y - 12.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b5_vals[i] * (-6.0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - 2.0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b6_vals[i] * (10.0 * fe_0 * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y + 4.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * fss * b7_vals[i] * (-2.0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b2_vals[i] * (6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + 12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_zz[i] += qpol_zz * faa_z * b2_vals[i] * 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y;

        fints_zz[i] += qpol_zz * faa_z * b3_vals[i] * (-36.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - 6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y - 18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b3_vals[i] * (-12.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z - 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z - 6.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - 18.0 * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_zz[i] += qpol_zz * faa_z * b3_vals[i] * (-9.0 * fe_0 * fe_0 * rpa_y * rpc_z - 12.0 * fe_0 * fe_0 * rpb_z * rpb_y - 15.0 * fe_0 * fe_0 * rpb_z * rpc_y - 6.0 * fe_0 * fe_0 * rpb_y * rpc_z - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b3_vals[i] * (-8.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - 4.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_zz[i] += qpol_zz * faa_z * b4_vals[i] * (36.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + 36.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + 36.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + 18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b4_vals[i] * (12.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + 18.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z + 24.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b4_vals[i] * (6.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z + 18.0 * fe_0 * fe_0 * rpa_y * rpc_z + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y + 30.0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b4_vals[i] * (12.0 * fe_0 * fe_0 * rpb_y * rpc_z + 15.0 * fe_0 * fe_0 * rpc_z * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 24.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 12.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b4_vals[i] * (4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 8.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b5_vals[i] * (-36.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - 36.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - 36.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - 18.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - 24.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b5_vals[i] * (-20.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - 24.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - 9.0 * fe_0 * fe_0 * rpa_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b5_vals[i] * (-6.0 * fe_0 * fe_0 * rpb_y * rpc_z - 30.0 * fe_0 * fe_0 * rpc_z * rpc_y - 24.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b5_vals[i] * (-24.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - 4.0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b6_vals[i] * (36.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + 20.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + 24.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 20.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + 15.0 * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b6_vals[i] * (12.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 24.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 12.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 8.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + 4.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b7_vals[i] * (-20.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 8.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += qpol_zz * faa_z * b8_vals[i] * 4.0 * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_zz[i] += (faa_zz * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_zz[i] += (faa_zz * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_zz[i] += (faa_zz * b2_vals[i] + faa * b1_vals[i]) * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_zz[i] += (faa_zz * b3_vals[i] + faa * b2_vals[i]) * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y - 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zz[i] += (faa_zz * b3_vals[i] + faa * b2_vals[i]) * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_zz[i] += (faa_zz * b3_vals[i] + faa * b2_vals[i]) * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y);

        fints_zz[i] += (faa_zz * b3_vals[i] + faa * b2_vals[i]) * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_zz[i] += (faa_zz * b3_vals[i] + faa * b2_vals[i]) * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_zz[i] += (faa_zz * b3_vals[i] + faa * b2_vals[i]) * (-(15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * (18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 6.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * ((27.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_zz[i] += (faa_zz * b4_vals[i] + faa * b3_vals[i]) * (6.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 4.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y - 18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y * rpc_y - 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpc_z - 15.0 * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z - 6.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-(45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 12.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b5_vals[i] + faa * b4_vals[i]) * (-3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b6_vals[i] + faa * b5_vals[i]) * (18.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b6_vals[i] + faa * b5_vals[i]) * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y + 12.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y + 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b6_vals[i] + faa * b5_vals[i]) * (6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_zz[i] += (faa_zz * b6_vals[i] + faa * b5_vals[i]) * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += (faa_zz * b6_vals[i] + faa * b5_vals[i]) * (5.0 * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpc_y + 12.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b6_vals[i] + faa * b5_vals[i]) * (6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b7_vals[i] + faa * b6_vals[i]) * (-9.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y - 10.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b7_vals[i] + faa * b6_vals[i]) * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z * rpc_y - (5.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b7_vals[i] + faa * b6_vals[i]) * (-6.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 4.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b8_vals[i] + faa * b7_vals[i]) * (5.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += (faa_zz * b9_vals[i] + faa * b8_vals[i]) * (-rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

    }
}

} // geom_npotrec namespace

