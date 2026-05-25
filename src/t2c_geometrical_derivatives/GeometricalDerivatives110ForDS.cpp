#include "GeometricalDerivatives110ForDS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_ds(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_ds,
                         const int idx_op_ss,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SS

    auto tr_0_0 = pbuffer.data(idx_op_ss);

    // Set up components of auxiliary buffer : PP

    auto tr_x_x = pbuffer.data(idx_op_pp);

    auto tr_x_y = pbuffer.data(idx_op_pp + 1);

    auto tr_x_z = pbuffer.data(idx_op_pp + 2);

    auto tr_y_x = pbuffer.data(idx_op_pp + 3);

    auto tr_y_y = pbuffer.data(idx_op_pp + 4);

    auto tr_y_z = pbuffer.data(idx_op_pp + 5);

    auto tr_z_x = pbuffer.data(idx_op_pp + 6);

    auto tr_z_y = pbuffer.data(idx_op_pp + 7);

    auto tr_z_z = pbuffer.data(idx_op_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto tr_xx_0 = pbuffer.data(idx_op_ds);

    auto tr_xy_0 = pbuffer.data(idx_op_ds + 1);

    auto tr_xz_0 = pbuffer.data(idx_op_ds + 2);

    auto tr_yy_0 = pbuffer.data(idx_op_ds + 3);

    auto tr_yz_0 = pbuffer.data(idx_op_ds + 4);

    auto tr_zz_0 = pbuffer.data(idx_op_ds + 5);

    // Set up components of auxiliary buffer : FP

    auto tr_xxx_x = pbuffer.data(idx_op_fp);

    auto tr_xxx_y = pbuffer.data(idx_op_fp + 1);

    auto tr_xxx_z = pbuffer.data(idx_op_fp + 2);

    auto tr_xxy_x = pbuffer.data(idx_op_fp + 3);

    auto tr_xxy_y = pbuffer.data(idx_op_fp + 4);

    auto tr_xxy_z = pbuffer.data(idx_op_fp + 5);

    auto tr_xxz_x = pbuffer.data(idx_op_fp + 6);

    auto tr_xxz_y = pbuffer.data(idx_op_fp + 7);

    auto tr_xxz_z = pbuffer.data(idx_op_fp + 8);

    auto tr_xyy_x = pbuffer.data(idx_op_fp + 9);

    auto tr_xyy_y = pbuffer.data(idx_op_fp + 10);

    auto tr_xyy_z = pbuffer.data(idx_op_fp + 11);

    auto tr_xyz_x = pbuffer.data(idx_op_fp + 12);

    auto tr_xyz_y = pbuffer.data(idx_op_fp + 13);

    auto tr_xyz_z = pbuffer.data(idx_op_fp + 14);

    auto tr_xzz_x = pbuffer.data(idx_op_fp + 15);

    auto tr_xzz_y = pbuffer.data(idx_op_fp + 16);

    auto tr_xzz_z = pbuffer.data(idx_op_fp + 17);

    auto tr_yyy_x = pbuffer.data(idx_op_fp + 18);

    auto tr_yyy_y = pbuffer.data(idx_op_fp + 19);

    auto tr_yyy_z = pbuffer.data(idx_op_fp + 20);

    auto tr_yyz_x = pbuffer.data(idx_op_fp + 21);

    auto tr_yyz_y = pbuffer.data(idx_op_fp + 22);

    auto tr_yyz_z = pbuffer.data(idx_op_fp + 23);

    auto tr_yzz_x = pbuffer.data(idx_op_fp + 24);

    auto tr_yzz_y = pbuffer.data(idx_op_fp + 25);

    auto tr_yzz_z = pbuffer.data(idx_op_fp + 26);

    auto tr_zzz_x = pbuffer.data(idx_op_fp + 27);

    auto tr_zzz_y = pbuffer.data(idx_op_fp + 28);

    auto tr_zzz_z = pbuffer.data(idx_op_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto tr_xxxx_0 = pbuffer.data(idx_op_gs);

    auto tr_xxxy_0 = pbuffer.data(idx_op_gs + 1);

    auto tr_xxxz_0 = pbuffer.data(idx_op_gs + 2);

    auto tr_xxyy_0 = pbuffer.data(idx_op_gs + 3);

    auto tr_xxyz_0 = pbuffer.data(idx_op_gs + 4);

    auto tr_xxzz_0 = pbuffer.data(idx_op_gs + 5);

    auto tr_xyyy_0 = pbuffer.data(idx_op_gs + 6);

    auto tr_xyyz_0 = pbuffer.data(idx_op_gs + 7);

    auto tr_xyzz_0 = pbuffer.data(idx_op_gs + 8);

    auto tr_xzzz_0 = pbuffer.data(idx_op_gs + 9);

    auto tr_yyyy_0 = pbuffer.data(idx_op_gs + 10);

    auto tr_yyyz_0 = pbuffer.data(idx_op_gs + 11);

    auto tr_yyzz_0 = pbuffer.data(idx_op_gs + 12);

    auto tr_yzzz_0 = pbuffer.data(idx_op_gs + 13);

    auto tr_zzzz_0 = pbuffer.data(idx_op_gs + 14);

    // Set up components of targeted buffer : DS

    auto tr_x_0_x_xx_0 = pbuffer.data(idx_op_geom_110_ds);

    auto tr_x_0_x_xy_0 = pbuffer.data(idx_op_geom_110_ds + 1);

    auto tr_x_0_x_xz_0 = pbuffer.data(idx_op_geom_110_ds + 2);

    auto tr_x_0_x_yy_0 = pbuffer.data(idx_op_geom_110_ds + 3);

    auto tr_x_0_x_yz_0 = pbuffer.data(idx_op_geom_110_ds + 4);

    auto tr_x_0_x_zz_0 = pbuffer.data(idx_op_geom_110_ds + 5);

    auto tr_x_0_y_xx_0 = pbuffer.data(idx_op_geom_110_ds + 6);

    auto tr_x_0_y_xy_0 = pbuffer.data(idx_op_geom_110_ds + 7);

    auto tr_x_0_y_xz_0 = pbuffer.data(idx_op_geom_110_ds + 8);

    auto tr_x_0_y_yy_0 = pbuffer.data(idx_op_geom_110_ds + 9);

    auto tr_x_0_y_yz_0 = pbuffer.data(idx_op_geom_110_ds + 10);

    auto tr_x_0_y_zz_0 = pbuffer.data(idx_op_geom_110_ds + 11);

    auto tr_x_0_z_xx_0 = pbuffer.data(idx_op_geom_110_ds + 12);

    auto tr_x_0_z_xy_0 = pbuffer.data(idx_op_geom_110_ds + 13);

    auto tr_x_0_z_xz_0 = pbuffer.data(idx_op_geom_110_ds + 14);

    auto tr_x_0_z_yy_0 = pbuffer.data(idx_op_geom_110_ds + 15);

    auto tr_x_0_z_yz_0 = pbuffer.data(idx_op_geom_110_ds + 16);

    auto tr_x_0_z_zz_0 = pbuffer.data(idx_op_geom_110_ds + 17);

    auto tr_y_0_x_xx_0 = pbuffer.data(idx_op_geom_110_ds + 18);

    auto tr_y_0_x_xy_0 = pbuffer.data(idx_op_geom_110_ds + 19);

    auto tr_y_0_x_xz_0 = pbuffer.data(idx_op_geom_110_ds + 20);

    auto tr_y_0_x_yy_0 = pbuffer.data(idx_op_geom_110_ds + 21);

    auto tr_y_0_x_yz_0 = pbuffer.data(idx_op_geom_110_ds + 22);

    auto tr_y_0_x_zz_0 = pbuffer.data(idx_op_geom_110_ds + 23);

    auto tr_y_0_y_xx_0 = pbuffer.data(idx_op_geom_110_ds + 24);

    auto tr_y_0_y_xy_0 = pbuffer.data(idx_op_geom_110_ds + 25);

    auto tr_y_0_y_xz_0 = pbuffer.data(idx_op_geom_110_ds + 26);

    auto tr_y_0_y_yy_0 = pbuffer.data(idx_op_geom_110_ds + 27);

    auto tr_y_0_y_yz_0 = pbuffer.data(idx_op_geom_110_ds + 28);

    auto tr_y_0_y_zz_0 = pbuffer.data(idx_op_geom_110_ds + 29);

    auto tr_y_0_z_xx_0 = pbuffer.data(idx_op_geom_110_ds + 30);

    auto tr_y_0_z_xy_0 = pbuffer.data(idx_op_geom_110_ds + 31);

    auto tr_y_0_z_xz_0 = pbuffer.data(idx_op_geom_110_ds + 32);

    auto tr_y_0_z_yy_0 = pbuffer.data(idx_op_geom_110_ds + 33);

    auto tr_y_0_z_yz_0 = pbuffer.data(idx_op_geom_110_ds + 34);

    auto tr_y_0_z_zz_0 = pbuffer.data(idx_op_geom_110_ds + 35);

    auto tr_z_0_x_xx_0 = pbuffer.data(idx_op_geom_110_ds + 36);

    auto tr_z_0_x_xy_0 = pbuffer.data(idx_op_geom_110_ds + 37);

    auto tr_z_0_x_xz_0 = pbuffer.data(idx_op_geom_110_ds + 38);

    auto tr_z_0_x_yy_0 = pbuffer.data(idx_op_geom_110_ds + 39);

    auto tr_z_0_x_yz_0 = pbuffer.data(idx_op_geom_110_ds + 40);

    auto tr_z_0_x_zz_0 = pbuffer.data(idx_op_geom_110_ds + 41);

    auto tr_z_0_y_xx_0 = pbuffer.data(idx_op_geom_110_ds + 42);

    auto tr_z_0_y_xy_0 = pbuffer.data(idx_op_geom_110_ds + 43);

    auto tr_z_0_y_xz_0 = pbuffer.data(idx_op_geom_110_ds + 44);

    auto tr_z_0_y_yy_0 = pbuffer.data(idx_op_geom_110_ds + 45);

    auto tr_z_0_y_yz_0 = pbuffer.data(idx_op_geom_110_ds + 46);

    auto tr_z_0_y_zz_0 = pbuffer.data(idx_op_geom_110_ds + 47);

    auto tr_z_0_z_xx_0 = pbuffer.data(idx_op_geom_110_ds + 48);

    auto tr_z_0_z_xy_0 = pbuffer.data(idx_op_geom_110_ds + 49);

    auto tr_z_0_z_xz_0 = pbuffer.data(idx_op_geom_110_ds + 50);

    auto tr_z_0_z_yy_0 = pbuffer.data(idx_op_geom_110_ds + 51);

    auto tr_z_0_z_yz_0 = pbuffer.data(idx_op_geom_110_ds + 52);

    auto tr_z_0_z_zz_0 = pbuffer.data(idx_op_geom_110_ds + 53);

    #pragma omp simd aligned(tr_0_0, tr_x_0_x_xx_0, tr_x_0_x_xy_0, tr_x_0_x_xz_0, tr_x_0_x_yy_0, tr_x_0_x_yz_0, tr_x_0_x_zz_0, tr_x_0_y_xx_0, tr_x_0_y_xy_0, tr_x_0_y_xz_0, tr_x_0_y_yy_0, tr_x_0_y_yz_0, tr_x_0_y_zz_0, tr_x_0_z_xx_0, tr_x_0_z_xy_0, tr_x_0_z_xz_0, tr_x_0_z_yy_0, tr_x_0_z_yz_0, tr_x_0_z_zz_0, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxy_0, tr_xxxz_0, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyz_0, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xy_0, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyz_0, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xz_0, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_y_0_x_xx_0, tr_y_0_x_xy_0, tr_y_0_x_xz_0, tr_y_0_x_yy_0, tr_y_0_x_yz_0, tr_y_0_x_zz_0, tr_y_0_y_xx_0, tr_y_0_y_xy_0, tr_y_0_y_xz_0, tr_y_0_y_yy_0, tr_y_0_y_yz_0, tr_y_0_y_zz_0, tr_y_0_z_xx_0, tr_y_0_z_xy_0, tr_y_0_z_xz_0, tr_y_0_z_yy_0, tr_y_0_z_yz_0, tr_y_0_z_zz_0, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyy_0, tr_yyyz_0, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yz_0, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_z_0_x_xx_0, tr_z_0_x_xy_0, tr_z_0_x_xz_0, tr_z_0_x_yy_0, tr_z_0_x_yz_0, tr_z_0_x_zz_0, tr_z_0_y_xx_0, tr_z_0_y_xy_0, tr_z_0_y_xz_0, tr_z_0_y_yy_0, tr_z_0_y_yz_0, tr_z_0_y_zz_0, tr_z_0_z_xx_0, tr_z_0_z_xy_0, tr_z_0_z_xz_0, tr_z_0_z_yy_0, tr_z_0_z_yz_0, tr_z_0_z_zz_0, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zzz_x, tr_zzz_y, tr_zzz_z, tr_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xx_0[i] = 2.0 * tr_0_0[i] - 4.0 * tr_x_x[i] * tke_0 - 10.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxx_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_0[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_0[i] = -2.0 * tr_y_x[i] * tke_0 - 6.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xxy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_0[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_0[i] = -2.0 * tr_z_x[i] * tke_0 - 6.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xxz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_0[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_0[i] = -2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_xyy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_0[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_0[i] = -2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_0[i] = -2.0 * tr_zz_0[i] * tbe_0 + 4.0 * tr_xzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_0[i] = -4.0 * tr_x_y[i] * tke_0 - 4.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xxx_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_0[i] = tr_0_0[i] - 2.0 * tr_y_y[i] * tke_0 - 2.0 * tr_yy_0[i] * tbe_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_0[i] = -2.0 * tr_z_y[i] * tke_0 - 2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_xxz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_0[i] = -4.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_0[i] = -2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_0[i] = 4.0 * tr_xzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_0[i] = -4.0 * tr_x_z[i] * tke_0 - 4.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xxx_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_0[i] = -2.0 * tr_y_z[i] * tke_0 - 2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_xxy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_0[i] = tr_0_0[i] - 2.0 * tr_z_z[i] * tke_0 - 2.0 * tr_zz_0[i] * tbe_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_0[i] = 4.0 * tr_xyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_0[i] = -2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_0[i] = -4.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_0[i] = -4.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xxy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_0[i] = tr_0_0[i] - 2.0 * tr_yy_0[i] * tbe_0 - 2.0 * tr_x_x[i] * tke_0 + 4.0 * tr_xyy_x[i] * tbe_0 * tke_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxyy_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_0[i] = -2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_0[i] = -4.0 * tr_y_x[i] * tke_0 + 4.0 * tr_yyy_x[i] * tbe_0 * tke_0 - 4.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyyy_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_0[i] = -2.0 * tr_z_x[i] * tke_0 + 4.0 * tr_yyz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_0[i] = 4.0 * tr_yzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_0[i] = -2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_0[i] = -2.0 * tr_x_y[i] * tke_0 - 6.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_0[i] = -2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_0[i] = 2.0 * tr_0_0[i] - 4.0 * tr_y_y[i] * tke_0 - 10.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_0[i] = -2.0 * tr_z_y[i] * tke_0 - 6.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_0[i] = -2.0 * tr_zz_0[i] * tbe_0 + 4.0 * tr_yzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_0[i] = 4.0 * tr_xxy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_0[i] = -2.0 * tr_x_z[i] * tke_0 - 2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_0[i] = -2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_0[i] = -4.0 * tr_y_z[i] * tke_0 - 4.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_0[i] = tr_0_0[i] - 2.0 * tr_z_z[i] * tke_0 - 2.0 * tr_zz_0[i] * tbe_0 - 2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_0[i] = -4.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_0[i] = -4.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xxz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_0[i] = -2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_0[i] = tr_0_0[i] - 2.0 * tr_zz_0[i] * tbe_0 - 2.0 * tr_x_x[i] * tke_0 + 4.0 * tr_xzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_0[i] = 4.0 * tr_yyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_0[i] = -2.0 * tr_y_x[i] * tke_0 + 4.0 * tr_yzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_0[i] = -4.0 * tr_z_x[i] * tke_0 + 4.0 * tr_zzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xzzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_0[i] = 4.0 * tr_xxz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_0[i] = -2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_0[i] = -2.0 * tr_x_y[i] * tke_0 + 4.0 * tr_xzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_0[i] = -4.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_0[i] = tr_0_0[i] - 2.0 * tr_zz_0[i] * tbe_0 - 2.0 * tr_y_y[i] * tke_0 + 4.0 * tr_yzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yyzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_0[i] = -4.0 * tr_z_y[i] * tke_0 + 4.0 * tr_zzz_y[i] * tbe_0 * tke_0 - 4.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yzzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_0[i] = -2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xxz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_0[i] = -2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_0[i] = -2.0 * tr_x_z[i] * tke_0 - 6.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_0[i] = -2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_0[i] = -2.0 * tr_y_z[i] * tke_0 - 6.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_0[i] = 2.0 * tr_0_0[i] - 4.0 * tr_z_z[i] * tke_0 - 10.0 * tr_zz_0[i] * tbe_0 + 4.0 * tr_zzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_0[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

