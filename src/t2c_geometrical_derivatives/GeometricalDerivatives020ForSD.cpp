#include "GeometricalDerivatives020ForSD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_sd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_sd,
                         const int idx_op_ss,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_pp,
                         const int idx_op_pf,
                         const int idx_op_dd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SS

    auto tr_0_0 = pbuffer.data(idx_op_ss);

    // Set up components of auxiliary buffer : SD

    auto tr_0_xx = pbuffer.data(idx_op_sd);

    auto tr_0_xy = pbuffer.data(idx_op_sd + 1);

    auto tr_0_xz = pbuffer.data(idx_op_sd + 2);

    auto tr_0_yy = pbuffer.data(idx_op_sd + 3);

    auto tr_0_yz = pbuffer.data(idx_op_sd + 4);

    auto tr_0_zz = pbuffer.data(idx_op_sd + 5);

    // Set up components of auxiliary buffer : SG

    auto tr_0_xxxx = pbuffer.data(idx_op_sg);

    auto tr_0_xxxy = pbuffer.data(idx_op_sg + 1);

    auto tr_0_xxxz = pbuffer.data(idx_op_sg + 2);

    auto tr_0_xxyy = pbuffer.data(idx_op_sg + 3);

    auto tr_0_xxyz = pbuffer.data(idx_op_sg + 4);

    auto tr_0_xxzz = pbuffer.data(idx_op_sg + 5);

    auto tr_0_xyyy = pbuffer.data(idx_op_sg + 6);

    auto tr_0_xyyz = pbuffer.data(idx_op_sg + 7);

    auto tr_0_xyzz = pbuffer.data(idx_op_sg + 8);

    auto tr_0_xzzz = pbuffer.data(idx_op_sg + 9);

    auto tr_0_yyyy = pbuffer.data(idx_op_sg + 10);

    auto tr_0_yyyz = pbuffer.data(idx_op_sg + 11);

    auto tr_0_yyzz = pbuffer.data(idx_op_sg + 12);

    auto tr_0_yzzz = pbuffer.data(idx_op_sg + 13);

    auto tr_0_zzzz = pbuffer.data(idx_op_sg + 14);

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

    // Set up components of auxiliary buffer : PF

    auto tr_x_xxx = pbuffer.data(idx_op_pf);

    auto tr_x_xxy = pbuffer.data(idx_op_pf + 1);

    auto tr_x_xxz = pbuffer.data(idx_op_pf + 2);

    auto tr_x_xyy = pbuffer.data(idx_op_pf + 3);

    auto tr_x_xyz = pbuffer.data(idx_op_pf + 4);

    auto tr_x_xzz = pbuffer.data(idx_op_pf + 5);

    auto tr_x_yyy = pbuffer.data(idx_op_pf + 6);

    auto tr_x_yyz = pbuffer.data(idx_op_pf + 7);

    auto tr_x_yzz = pbuffer.data(idx_op_pf + 8);

    auto tr_x_zzz = pbuffer.data(idx_op_pf + 9);

    auto tr_y_xxx = pbuffer.data(idx_op_pf + 10);

    auto tr_y_xxy = pbuffer.data(idx_op_pf + 11);

    auto tr_y_xxz = pbuffer.data(idx_op_pf + 12);

    auto tr_y_xyy = pbuffer.data(idx_op_pf + 13);

    auto tr_y_xyz = pbuffer.data(idx_op_pf + 14);

    auto tr_y_xzz = pbuffer.data(idx_op_pf + 15);

    auto tr_y_yyy = pbuffer.data(idx_op_pf + 16);

    auto tr_y_yyz = pbuffer.data(idx_op_pf + 17);

    auto tr_y_yzz = pbuffer.data(idx_op_pf + 18);

    auto tr_y_zzz = pbuffer.data(idx_op_pf + 19);

    auto tr_z_xxx = pbuffer.data(idx_op_pf + 20);

    auto tr_z_xxy = pbuffer.data(idx_op_pf + 21);

    auto tr_z_xxz = pbuffer.data(idx_op_pf + 22);

    auto tr_z_xyy = pbuffer.data(idx_op_pf + 23);

    auto tr_z_xyz = pbuffer.data(idx_op_pf + 24);

    auto tr_z_xzz = pbuffer.data(idx_op_pf + 25);

    auto tr_z_yyy = pbuffer.data(idx_op_pf + 26);

    auto tr_z_yyz = pbuffer.data(idx_op_pf + 27);

    auto tr_z_yzz = pbuffer.data(idx_op_pf + 28);

    auto tr_z_zzz = pbuffer.data(idx_op_pf + 29);

    // Set up components of auxiliary buffer : DD

    auto tr_xx_xx = pbuffer.data(idx_op_dd);

    auto tr_xx_xy = pbuffer.data(idx_op_dd + 1);

    auto tr_xx_xz = pbuffer.data(idx_op_dd + 2);

    auto tr_xx_yy = pbuffer.data(idx_op_dd + 3);

    auto tr_xx_yz = pbuffer.data(idx_op_dd + 4);

    auto tr_xx_zz = pbuffer.data(idx_op_dd + 5);

    auto tr_xy_xx = pbuffer.data(idx_op_dd + 6);

    auto tr_xy_xy = pbuffer.data(idx_op_dd + 7);

    auto tr_xy_xz = pbuffer.data(idx_op_dd + 8);

    auto tr_xy_yy = pbuffer.data(idx_op_dd + 9);

    auto tr_xy_yz = pbuffer.data(idx_op_dd + 10);

    auto tr_xy_zz = pbuffer.data(idx_op_dd + 11);

    auto tr_xz_xx = pbuffer.data(idx_op_dd + 12);

    auto tr_xz_xy = pbuffer.data(idx_op_dd + 13);

    auto tr_xz_xz = pbuffer.data(idx_op_dd + 14);

    auto tr_xz_yy = pbuffer.data(idx_op_dd + 15);

    auto tr_xz_yz = pbuffer.data(idx_op_dd + 16);

    auto tr_xz_zz = pbuffer.data(idx_op_dd + 17);

    auto tr_yy_xx = pbuffer.data(idx_op_dd + 18);

    auto tr_yy_xy = pbuffer.data(idx_op_dd + 19);

    auto tr_yy_xz = pbuffer.data(idx_op_dd + 20);

    auto tr_yy_yy = pbuffer.data(idx_op_dd + 21);

    auto tr_yy_yz = pbuffer.data(idx_op_dd + 22);

    auto tr_yy_zz = pbuffer.data(idx_op_dd + 23);

    auto tr_yz_xx = pbuffer.data(idx_op_dd + 24);

    auto tr_yz_xy = pbuffer.data(idx_op_dd + 25);

    auto tr_yz_xz = pbuffer.data(idx_op_dd + 26);

    auto tr_yz_yy = pbuffer.data(idx_op_dd + 27);

    auto tr_yz_yz = pbuffer.data(idx_op_dd + 28);

    auto tr_yz_zz = pbuffer.data(idx_op_dd + 29);

    auto tr_zz_xx = pbuffer.data(idx_op_dd + 30);

    auto tr_zz_xy = pbuffer.data(idx_op_dd + 31);

    auto tr_zz_xz = pbuffer.data(idx_op_dd + 32);

    auto tr_zz_yy = pbuffer.data(idx_op_dd + 33);

    auto tr_zz_yz = pbuffer.data(idx_op_dd + 34);

    auto tr_zz_zz = pbuffer.data(idx_op_dd + 35);

    // Set up components of targeted buffer : SD

    auto tr_0_0_xx_0_xx = pbuffer.data(idx_op_geom_020_sd);

    auto tr_0_0_xx_0_xy = pbuffer.data(idx_op_geom_020_sd + 1);

    auto tr_0_0_xx_0_xz = pbuffer.data(idx_op_geom_020_sd + 2);

    auto tr_0_0_xx_0_yy = pbuffer.data(idx_op_geom_020_sd + 3);

    auto tr_0_0_xx_0_yz = pbuffer.data(idx_op_geom_020_sd + 4);

    auto tr_0_0_xx_0_zz = pbuffer.data(idx_op_geom_020_sd + 5);

    auto tr_0_0_xy_0_xx = pbuffer.data(idx_op_geom_020_sd + 6);

    auto tr_0_0_xy_0_xy = pbuffer.data(idx_op_geom_020_sd + 7);

    auto tr_0_0_xy_0_xz = pbuffer.data(idx_op_geom_020_sd + 8);

    auto tr_0_0_xy_0_yy = pbuffer.data(idx_op_geom_020_sd + 9);

    auto tr_0_0_xy_0_yz = pbuffer.data(idx_op_geom_020_sd + 10);

    auto tr_0_0_xy_0_zz = pbuffer.data(idx_op_geom_020_sd + 11);

    auto tr_0_0_xz_0_xx = pbuffer.data(idx_op_geom_020_sd + 12);

    auto tr_0_0_xz_0_xy = pbuffer.data(idx_op_geom_020_sd + 13);

    auto tr_0_0_xz_0_xz = pbuffer.data(idx_op_geom_020_sd + 14);

    auto tr_0_0_xz_0_yy = pbuffer.data(idx_op_geom_020_sd + 15);

    auto tr_0_0_xz_0_yz = pbuffer.data(idx_op_geom_020_sd + 16);

    auto tr_0_0_xz_0_zz = pbuffer.data(idx_op_geom_020_sd + 17);

    auto tr_0_0_yy_0_xx = pbuffer.data(idx_op_geom_020_sd + 18);

    auto tr_0_0_yy_0_xy = pbuffer.data(idx_op_geom_020_sd + 19);

    auto tr_0_0_yy_0_xz = pbuffer.data(idx_op_geom_020_sd + 20);

    auto tr_0_0_yy_0_yy = pbuffer.data(idx_op_geom_020_sd + 21);

    auto tr_0_0_yy_0_yz = pbuffer.data(idx_op_geom_020_sd + 22);

    auto tr_0_0_yy_0_zz = pbuffer.data(idx_op_geom_020_sd + 23);

    auto tr_0_0_yz_0_xx = pbuffer.data(idx_op_geom_020_sd + 24);

    auto tr_0_0_yz_0_xy = pbuffer.data(idx_op_geom_020_sd + 25);

    auto tr_0_0_yz_0_xz = pbuffer.data(idx_op_geom_020_sd + 26);

    auto tr_0_0_yz_0_yy = pbuffer.data(idx_op_geom_020_sd + 27);

    auto tr_0_0_yz_0_yz = pbuffer.data(idx_op_geom_020_sd + 28);

    auto tr_0_0_yz_0_zz = pbuffer.data(idx_op_geom_020_sd + 29);

    auto tr_0_0_zz_0_xx = pbuffer.data(idx_op_geom_020_sd + 30);

    auto tr_0_0_zz_0_xy = pbuffer.data(idx_op_geom_020_sd + 31);

    auto tr_0_0_zz_0_xz = pbuffer.data(idx_op_geom_020_sd + 32);

    auto tr_0_0_zz_0_yy = pbuffer.data(idx_op_geom_020_sd + 33);

    auto tr_0_0_zz_0_yz = pbuffer.data(idx_op_geom_020_sd + 34);

    auto tr_0_0_zz_0_zz = pbuffer.data(idx_op_geom_020_sd + 35);

    #pragma omp simd aligned(tr_0_0, tr_0_0_xx_0_xx, tr_0_0_xx_0_xy, tr_0_0_xx_0_xz, tr_0_0_xx_0_yy, tr_0_0_xx_0_yz, tr_0_0_xx_0_zz, tr_0_0_xy_0_xx, tr_0_0_xy_0_xy, tr_0_0_xy_0_xz, tr_0_0_xy_0_yy, tr_0_0_xy_0_yz, tr_0_0_xy_0_zz, tr_0_0_xz_0_xx, tr_0_0_xz_0_xy, tr_0_0_xz_0_xz, tr_0_0_xz_0_yy, tr_0_0_xz_0_yz, tr_0_0_xz_0_zz, tr_0_0_yy_0_xx, tr_0_0_yy_0_xy, tr_0_0_yy_0_xz, tr_0_0_yy_0_yy, tr_0_0_yy_0_yz, tr_0_0_yy_0_zz, tr_0_0_yz_0_xx, tr_0_0_yz_0_xy, tr_0_0_yz_0_xz, tr_0_0_yz_0_yy, tr_0_0_yz_0_yz, tr_0_0_yz_0_zz, tr_0_0_zz_0_xx, tr_0_0_zz_0_xy, tr_0_0_zz_0_xz, tr_0_0_zz_0_yy, tr_0_0_zz_0_yz, tr_0_0_zz_0_zz, tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_x_x, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_y_x, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_z_x, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_0_xx[i] = 2.0 * tr_0_0[i] - 2.0 * tr_0_xx[i] * tbe_0 - 10.0 * tr_0_xx[i] * tke_0 + 4.0 * tr_0_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_x_x[i] * tbe_0 + 8.0 * tr_x_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xy[i] = -2.0 * tr_0_xy[i] * tbe_0 - 6.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_x_y[i] * tbe_0 + 8.0 * tr_x_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xz[i] = -2.0 * tr_0_xz[i] * tbe_0 - 6.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_x_z[i] * tbe_0 + 8.0 * tr_x_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_yy[i] = -2.0 * tr_0_yy[i] * tbe_0 - 2.0 * tr_0_yy[i] * tke_0 + 4.0 * tr_0_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_x_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_yz[i] = -2.0 * tr_0_yz[i] * tbe_0 - 2.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_x_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_zz[i] = -2.0 * tr_0_zz[i] * tbe_0 - 2.0 * tr_0_zz[i] * tke_0 + 4.0 * tr_0_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_x_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_zz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xx[i] = -4.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_y_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xy[i] = tr_0_0[i] - 2.0 * tr_0_yy[i] * tke_0 - 2.0 * tr_0_xx[i] * tke_0 + 4.0 * tr_0_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_y_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_x_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xz[i] = -2.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_y_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_x_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_yy[i] = -4.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_y_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_x_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_yz[i] = -2.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_y_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_x_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_zz[i] = 4.0 * tr_0_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_y_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_x_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_zz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xx[i] = -4.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_z_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xy[i] = -2.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_z_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_x_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xz[i] = tr_0_0[i] - 2.0 * tr_0_zz[i] * tke_0 - 2.0 * tr_0_xx[i] * tke_0 + 4.0 * tr_0_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_z_z[i] * tbe_0 + 4.0 * tr_z_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_x_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_yy[i] = 4.0 * tr_0_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_x_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_yz[i] = -2.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_x_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_zz[i] = -4.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_x_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_zz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xx[i] = -2.0 * tr_0_xx[i] * tbe_0 - 2.0 * tr_0_xx[i] * tke_0 + 4.0 * tr_0_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_y_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xy[i] = -2.0 * tr_0_xy[i] * tbe_0 - 6.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_y_x[i] * tbe_0 + 8.0 * tr_y_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xz[i] = -2.0 * tr_0_xz[i] * tbe_0 - 2.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_y_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_yy[i] = 2.0 * tr_0_0[i] - 2.0 * tr_0_yy[i] * tbe_0 - 10.0 * tr_0_yy[i] * tke_0 + 4.0 * tr_0_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_y_y[i] * tbe_0 + 8.0 * tr_y_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_yz[i] = -2.0 * tr_0_yz[i] * tbe_0 - 6.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_y_z[i] * tbe_0 + 8.0 * tr_y_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_zz[i] = -2.0 * tr_0_zz[i] * tbe_0 - 2.0 * tr_0_zz[i] * tke_0 + 4.0 * tr_0_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_y_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_zz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xx[i] = 4.0 * tr_0_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_z_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_y_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xy[i] = -2.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_z_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_y_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xz[i] = -2.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_y_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_yy[i] = -4.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_z_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_y_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_yz[i] = tr_0_0[i] - 2.0 * tr_0_zz[i] * tke_0 - 2.0 * tr_0_yy[i] * tke_0 + 4.0 * tr_0_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_z_z[i] * tbe_0 + 4.0 * tr_z_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_y_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_zz[i] = -4.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_y_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_zz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xx[i] = -2.0 * tr_0_xx[i] * tbe_0 - 2.0 * tr_0_xx[i] * tke_0 + 4.0 * tr_0_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_z_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xy[i] = -2.0 * tr_0_xy[i] * tbe_0 - 2.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_0_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_z_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xz[i] = -2.0 * tr_0_xz[i] * tbe_0 - 6.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_0_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_x[i] * tbe_0 + 8.0 * tr_z_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_yy[i] = -2.0 * tr_0_yy[i] * tbe_0 - 2.0 * tr_0_yy[i] * tke_0 + 4.0 * tr_0_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_z_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_yz[i] = -2.0 * tr_0_yz[i] * tbe_0 - 6.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_0_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_y[i] * tbe_0 + 8.0 * tr_z_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_zz[i] = 2.0 * tr_0_0[i] - 2.0 * tr_0_zz[i] * tbe_0 - 10.0 * tr_0_zz[i] * tke_0 + 4.0 * tr_0_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_z_z[i] * tbe_0 + 8.0 * tr_z_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_zz[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

