#include "GeometricalDerivatives110ForPP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_pp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_pp,
                         const int idx_op_ss,
                         const int idx_op_sd,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_fp,
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

    // Set up 0-3 components of targeted buffer : PP

    auto tr_x_0_x_x_x = pbuffer.data(idx_op_geom_110_pp);

    auto tr_x_0_x_x_y = pbuffer.data(idx_op_geom_110_pp + 1);

    auto tr_x_0_x_x_z = pbuffer.data(idx_op_geom_110_pp + 2);

    #pragma omp simd aligned(tr_0_0, tr_0_xx, tr_0_xy, tr_0_xz, tr_x_0_x_x_x, tr_x_0_x_x_y, tr_x_0_x_x_z, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_x_x[i] = tr_0_0[i] - 2.0 * tr_0_xx[i] * tke_0 - 6.0 * tr_x_x[i] * tbe_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xx_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_y[i] = -2.0 * tr_0_xy[i] * tke_0 - 6.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_xx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_z[i] = -2.0 * tr_0_xz[i] * tke_0 - 6.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_xx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_z[i] * tbe_0 * tbe_0;
    }

    // Set up 3-6 components of targeted buffer : PP

    auto tr_x_0_x_y_x = pbuffer.data(idx_op_geom_110_pp + 3);

    auto tr_x_0_x_y_y = pbuffer.data(idx_op_geom_110_pp + 4);

    auto tr_x_0_x_y_z = pbuffer.data(idx_op_geom_110_pp + 5);

    #pragma omp simd aligned(tr_x_0_x_y_x, tr_x_0_x_y_y, tr_x_0_x_y_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_y_x, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_y_x[i] = -2.0 * tr_y_x[i] * tbe_0 - 2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_y[i] = -2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_xy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_z[i] = -2.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_xy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 6-9 components of targeted buffer : PP

    auto tr_x_0_x_z_x = pbuffer.data(idx_op_geom_110_pp + 6);

    auto tr_x_0_x_z_y = pbuffer.data(idx_op_geom_110_pp + 7);

    auto tr_x_0_x_z_z = pbuffer.data(idx_op_geom_110_pp + 8);

    #pragma omp simd aligned(tr_x_0_x_z_x, tr_x_0_x_z_y, tr_x_0_x_z_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_z_x[i] = -2.0 * tr_z_x[i] * tbe_0 - 2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_y[i] = -2.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_xz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_z[i] = -2.0 * tr_z_z[i] * tbe_0 + 4.0 * tr_xz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 9-12 components of targeted buffer : PP

    auto tr_x_0_y_x_x = pbuffer.data(idx_op_geom_110_pp + 9);

    auto tr_x_0_y_x_y = pbuffer.data(idx_op_geom_110_pp + 10);

    auto tr_x_0_y_x_z = pbuffer.data(idx_op_geom_110_pp + 11);

    #pragma omp simd aligned(tr_0_0, tr_0_xy, tr_0_yy, tr_0_yz, tr_x_0_y_x_x, tr_x_0_y_x_y, tr_x_0_y_x_z, tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_y_x, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_x_x[i] = -2.0 * tr_0_xy[i] * tke_0 - 2.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_xx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_y[i] = tr_0_0[i] - 2.0 * tr_0_yy[i] * tke_0 - 2.0 * tr_y_y[i] * tbe_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xx_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_z[i] = -2.0 * tr_0_yz[i] * tke_0 - 2.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_xx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 12-15 components of targeted buffer : PP

    auto tr_x_0_y_y_x = pbuffer.data(idx_op_geom_110_pp + 12);

    auto tr_x_0_y_y_y = pbuffer.data(idx_op_geom_110_pp + 13);

    auto tr_x_0_y_y_z = pbuffer.data(idx_op_geom_110_pp + 14);

    #pragma omp simd aligned(tr_x_0_y_y_x, tr_x_0_y_y_y, tr_x_0_y_y_z, tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_y_x[i] = -2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_xy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_y[i] = -2.0 * tr_x_y[i] * tbe_0 - 2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_z[i] = -2.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_xy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 15-18 components of targeted buffer : PP

    auto tr_x_0_y_z_x = pbuffer.data(idx_op_geom_110_pp + 15);

    auto tr_x_0_y_z_y = pbuffer.data(idx_op_geom_110_pp + 16);

    auto tr_x_0_y_z_z = pbuffer.data(idx_op_geom_110_pp + 17);

    #pragma omp simd aligned(tr_x_0_y_z_x, tr_x_0_y_z_y, tr_x_0_y_z_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_z_x[i] = 4.0 * tr_xz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_y[i] = -2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_z[i] = 4.0 * tr_xz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 18-21 components of targeted buffer : PP

    auto tr_x_0_z_x_x = pbuffer.data(idx_op_geom_110_pp + 18);

    auto tr_x_0_z_x_y = pbuffer.data(idx_op_geom_110_pp + 19);

    auto tr_x_0_z_x_z = pbuffer.data(idx_op_geom_110_pp + 20);

    #pragma omp simd aligned(tr_0_0, tr_0_xz, tr_0_yz, tr_0_zz, tr_x_0_z_x_x, tr_x_0_z_x_y, tr_x_0_z_x_z, tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_x_x[i] = -2.0 * tr_0_xz[i] * tke_0 - 2.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_xx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_y[i] = -2.0 * tr_0_yz[i] * tke_0 - 2.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_xx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_z[i] = tr_0_0[i] - 2.0 * tr_0_zz[i] * tke_0 - 2.0 * tr_z_z[i] * tbe_0 - 2.0 * tr_xx_0[i] * tbe_0 + 4.0 * tr_xx_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 21-24 components of targeted buffer : PP

    auto tr_x_0_z_y_x = pbuffer.data(idx_op_geom_110_pp + 21);

    auto tr_x_0_z_y_y = pbuffer.data(idx_op_geom_110_pp + 22);

    auto tr_x_0_z_y_z = pbuffer.data(idx_op_geom_110_pp + 23);

    #pragma omp simd aligned(tr_x_0_z_y_x, tr_x_0_z_y_y, tr_x_0_z_y_z, tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_y_x[i] = 4.0 * tr_xy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_y[i] = 4.0 * tr_xy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_z[i] = -2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 24-27 components of targeted buffer : PP

    auto tr_x_0_z_z_x = pbuffer.data(idx_op_geom_110_pp + 24);

    auto tr_x_0_z_z_y = pbuffer.data(idx_op_geom_110_pp + 25);

    auto tr_x_0_z_z_z = pbuffer.data(idx_op_geom_110_pp + 26);

    #pragma omp simd aligned(tr_x_0_z_z_x, tr_x_0_z_z_y, tr_x_0_z_z_z, tr_x_x, tr_x_y, tr_x_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_z_x[i] = -2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_xz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_y[i] = -2.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_xz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_z[i] = -2.0 * tr_x_z[i] * tbe_0 - 2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 27-30 components of targeted buffer : PP

    auto tr_y_0_x_x_x = pbuffer.data(idx_op_geom_110_pp + 27);

    auto tr_y_0_x_x_y = pbuffer.data(idx_op_geom_110_pp + 28);

    auto tr_y_0_x_x_z = pbuffer.data(idx_op_geom_110_pp + 29);

    #pragma omp simd aligned(tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_y_0_x_x_x, tr_y_0_x_x_y, tr_y_0_x_x_z, tr_y_x, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_x_x[i] = -2.0 * tr_y_x[i] * tbe_0 - 2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_y[i] = -2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_xy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_z[i] = -2.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_xy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 30-33 components of targeted buffer : PP

    auto tr_y_0_x_y_x = pbuffer.data(idx_op_geom_110_pp + 30);

    auto tr_y_0_x_y_y = pbuffer.data(idx_op_geom_110_pp + 31);

    auto tr_y_0_x_y_z = pbuffer.data(idx_op_geom_110_pp + 32);

    #pragma omp simd aligned(tr_0_0, tr_0_xx, tr_0_xy, tr_0_xz, tr_x_x, tr_x_y, tr_x_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_y_0_x_y_x, tr_y_0_x_y_y, tr_y_0_x_y_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_y_x[i] = tr_0_0[i] - 2.0 * tr_0_xx[i] * tke_0 - 2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yy_xx[i] * tbe_0 * tke_0 - 2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_xyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_y[i] = -2.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_yy_xy[i] * tbe_0 * tke_0 - 2.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_xyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_z[i] = -2.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_yy_xz[i] * tbe_0 * tke_0 - 2.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_xyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 33-36 components of targeted buffer : PP

    auto tr_y_0_x_z_x = pbuffer.data(idx_op_geom_110_pp + 33);

    auto tr_y_0_x_z_y = pbuffer.data(idx_op_geom_110_pp + 34);

    auto tr_y_0_x_z_z = pbuffer.data(idx_op_geom_110_pp + 35);

    #pragma omp simd aligned(tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_y_0_x_z_x, tr_y_0_x_z_y, tr_y_0_x_z_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_z_x[i] = -2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_y[i] = 4.0 * tr_yz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_z[i] = 4.0 * tr_yz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 36-39 components of targeted buffer : PP

    auto tr_y_0_y_x_x = pbuffer.data(idx_op_geom_110_pp + 36);

    auto tr_y_0_y_x_y = pbuffer.data(idx_op_geom_110_pp + 37);

    auto tr_y_0_y_x_z = pbuffer.data(idx_op_geom_110_pp + 38);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_y_0_y_x_x, tr_y_0_y_x_y, tr_y_0_y_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_x_x[i] = -2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_xy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_y[i] = -2.0 * tr_x_y[i] * tbe_0 - 2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_z[i] = -2.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_xy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 39-42 components of targeted buffer : PP

    auto tr_y_0_y_y_x = pbuffer.data(idx_op_geom_110_pp + 39);

    auto tr_y_0_y_y_y = pbuffer.data(idx_op_geom_110_pp + 40);

    auto tr_y_0_y_y_z = pbuffer.data(idx_op_geom_110_pp + 41);

    #pragma omp simd aligned(tr_0_0, tr_0_xy, tr_0_yy, tr_0_yz, tr_y_0_y_y_x, tr_y_0_y_y_y, tr_y_0_y_y_z, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_y_x[i] = -2.0 * tr_0_xy[i] * tke_0 - 6.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_yy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_y[i] = tr_0_0[i] - 2.0 * tr_0_yy[i] * tke_0 - 6.0 * tr_y_y[i] * tbe_0 - 2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_z[i] = -2.0 * tr_0_yz[i] * tke_0 - 6.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_yy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 42-45 components of targeted buffer : PP

    auto tr_y_0_y_z_x = pbuffer.data(idx_op_geom_110_pp + 42);

    auto tr_y_0_y_z_y = pbuffer.data(idx_op_geom_110_pp + 43);

    auto tr_y_0_y_z_z = pbuffer.data(idx_op_geom_110_pp + 44);

    #pragma omp simd aligned(tr_y_0_y_z_x, tr_y_0_y_z_y, tr_y_0_y_z_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_z_x[i] = -2.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_yz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_y[i] = -2.0 * tr_z_y[i] * tbe_0 - 2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_z[i] = -2.0 * tr_z_z[i] * tbe_0 + 4.0 * tr_yz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 45-48 components of targeted buffer : PP

    auto tr_y_0_z_x_x = pbuffer.data(idx_op_geom_110_pp + 45);

    auto tr_y_0_z_x_y = pbuffer.data(idx_op_geom_110_pp + 46);

    auto tr_y_0_z_x_z = pbuffer.data(idx_op_geom_110_pp + 47);

    #pragma omp simd aligned(tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_y_0_z_x_x, tr_y_0_z_x_y, tr_y_0_z_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_x_x[i] = 4.0 * tr_xy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_y[i] = 4.0 * tr_xy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_z[i] = -2.0 * tr_xy_0[i] * tbe_0 + 4.0 * tr_xy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 48-51 components of targeted buffer : PP

    auto tr_y_0_z_y_x = pbuffer.data(idx_op_geom_110_pp + 48);

    auto tr_y_0_z_y_y = pbuffer.data(idx_op_geom_110_pp + 49);

    auto tr_y_0_z_y_z = pbuffer.data(idx_op_geom_110_pp + 50);

    #pragma omp simd aligned(tr_0_0, tr_0_xz, tr_0_yz, tr_0_zz, tr_y_0_z_y_x, tr_y_0_z_y_y, tr_y_0_z_y_z, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_y_x[i] = -2.0 * tr_0_xz[i] * tke_0 - 2.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_yy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_y[i] = -2.0 * tr_0_yz[i] * tke_0 - 2.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_yy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_z[i] = tr_0_0[i] - 2.0 * tr_0_zz[i] * tke_0 - 2.0 * tr_z_z[i] * tbe_0 - 2.0 * tr_yy_0[i] * tbe_0 + 4.0 * tr_yy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 51-54 components of targeted buffer : PP

    auto tr_y_0_z_z_x = pbuffer.data(idx_op_geom_110_pp + 51);

    auto tr_y_0_z_z_y = pbuffer.data(idx_op_geom_110_pp + 52);

    auto tr_y_0_z_z_z = pbuffer.data(idx_op_geom_110_pp + 53);

    #pragma omp simd aligned(tr_y_0_z_z_x, tr_y_0_z_z_y, tr_y_0_z_z_z, tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_z_x[i] = -2.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_yz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_y[i] = -2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_yz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_z[i] = -2.0 * tr_y_z[i] * tbe_0 - 2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 54-57 components of targeted buffer : PP

    auto tr_z_0_x_x_x = pbuffer.data(idx_op_geom_110_pp + 54);

    auto tr_z_0_x_x_y = pbuffer.data(idx_op_geom_110_pp + 55);

    auto tr_z_0_x_x_z = pbuffer.data(idx_op_geom_110_pp + 56);

    #pragma omp simd aligned(tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_z_0_x_x_x, tr_z_0_x_x_y, tr_z_0_x_x_z, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_x_x[i] = -2.0 * tr_z_x[i] * tbe_0 - 2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_y[i] = -2.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_xz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_z[i] = -2.0 * tr_z_z[i] * tbe_0 + 4.0 * tr_xz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 57-60 components of targeted buffer : PP

    auto tr_z_0_x_y_x = pbuffer.data(idx_op_geom_110_pp + 57);

    auto tr_z_0_x_y_y = pbuffer.data(idx_op_geom_110_pp + 58);

    auto tr_z_0_x_y_z = pbuffer.data(idx_op_geom_110_pp + 59);

    #pragma omp simd aligned(tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_z_0_x_y_x, tr_z_0_x_y_y, tr_z_0_x_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_y_x[i] = -2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_y[i] = 4.0 * tr_yz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_z[i] = 4.0 * tr_yz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 60-63 components of targeted buffer : PP

    auto tr_z_0_x_z_x = pbuffer.data(idx_op_geom_110_pp + 60);

    auto tr_z_0_x_z_y = pbuffer.data(idx_op_geom_110_pp + 61);

    auto tr_z_0_x_z_z = pbuffer.data(idx_op_geom_110_pp + 62);

    #pragma omp simd aligned(tr_0_0, tr_0_xx, tr_0_xy, tr_0_xz, tr_x_x, tr_x_y, tr_x_z, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_z_0_x_z_x, tr_z_0_x_z_y, tr_z_0_x_z_z, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_z_x[i] = tr_0_0[i] - 2.0 * tr_0_xx[i] * tke_0 - 2.0 * tr_zz_0[i] * tbe_0 + 4.0 * tr_zz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_xzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_y[i] = -2.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_zz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_xzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_z[i] = -2.0 * tr_0_xz[i] * tke_0 + 4.0 * tr_zz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_x_z[i] * tbe_0 + 4.0 * tr_xzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 63-66 components of targeted buffer : PP

    auto tr_z_0_y_x_x = pbuffer.data(idx_op_geom_110_pp + 63);

    auto tr_z_0_y_x_y = pbuffer.data(idx_op_geom_110_pp + 64);

    auto tr_z_0_y_x_z = pbuffer.data(idx_op_geom_110_pp + 65);

    #pragma omp simd aligned(tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, tr_z_0_y_x_x, tr_z_0_y_x_y, tr_z_0_y_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_x_x[i] = 4.0 * tr_xz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_y[i] = -2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_z[i] = 4.0 * tr_xz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 66-69 components of targeted buffer : PP

    auto tr_z_0_y_y_x = pbuffer.data(idx_op_geom_110_pp + 66);

    auto tr_z_0_y_y_y = pbuffer.data(idx_op_geom_110_pp + 67);

    auto tr_z_0_y_y_z = pbuffer.data(idx_op_geom_110_pp + 68);

    #pragma omp simd aligned(tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_0_y_y_x, tr_z_0_y_y_y, tr_z_0_y_y_z, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_y_x[i] = -2.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_yz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_y[i] = -2.0 * tr_z_y[i] * tbe_0 - 2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_z[i] = -2.0 * tr_z_z[i] * tbe_0 + 4.0 * tr_yz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 69-72 components of targeted buffer : PP

    auto tr_z_0_y_z_x = pbuffer.data(idx_op_geom_110_pp + 69);

    auto tr_z_0_y_z_y = pbuffer.data(idx_op_geom_110_pp + 70);

    auto tr_z_0_y_z_z = pbuffer.data(idx_op_geom_110_pp + 71);

    #pragma omp simd aligned(tr_0_0, tr_0_xy, tr_0_yy, tr_0_yz, tr_y_x, tr_y_y, tr_y_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_z_0_y_z_x, tr_z_0_y_z_y, tr_z_0_y_z_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_z_x[i] = -2.0 * tr_0_xy[i] * tke_0 + 4.0 * tr_zz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_yzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_y[i] = tr_0_0[i] - 2.0 * tr_0_yy[i] * tke_0 - 2.0 * tr_zz_0[i] * tbe_0 + 4.0 * tr_zz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_yzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_z[i] = -2.0 * tr_0_yz[i] * tke_0 + 4.0 * tr_zz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_y_z[i] * tbe_0 + 4.0 * tr_yzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 72-75 components of targeted buffer : PP

    auto tr_z_0_z_x_x = pbuffer.data(idx_op_geom_110_pp + 72);

    auto tr_z_0_z_x_y = pbuffer.data(idx_op_geom_110_pp + 73);

    auto tr_z_0_z_x_z = pbuffer.data(idx_op_geom_110_pp + 74);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_z_0_z_x_x, tr_z_0_z_x_y, tr_z_0_z_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_x_x[i] = -2.0 * tr_x_x[i] * tbe_0 + 4.0 * tr_xz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_y[i] = -2.0 * tr_x_y[i] * tbe_0 + 4.0 * tr_xz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_z[i] = -2.0 * tr_x_z[i] * tbe_0 - 2.0 * tr_xz_0[i] * tbe_0 + 4.0 * tr_xz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 75-78 components of targeted buffer : PP

    auto tr_z_0_z_y_x = pbuffer.data(idx_op_geom_110_pp + 75);

    auto tr_z_0_z_y_y = pbuffer.data(idx_op_geom_110_pp + 76);

    auto tr_z_0_z_y_z = pbuffer.data(idx_op_geom_110_pp + 77);

    #pragma omp simd aligned(tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_z_0_z_y_x, tr_z_0_z_y_y, tr_z_0_z_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_y_x[i] = -2.0 * tr_y_x[i] * tbe_0 + 4.0 * tr_yz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_y[i] = -2.0 * tr_y_y[i] * tbe_0 + 4.0 * tr_yz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_z[i] = -2.0 * tr_y_z[i] * tbe_0 - 2.0 * tr_yz_0[i] * tbe_0 + 4.0 * tr_yz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 78-81 components of targeted buffer : PP

    auto tr_z_0_z_z_x = pbuffer.data(idx_op_geom_110_pp + 78);

    auto tr_z_0_z_z_y = pbuffer.data(idx_op_geom_110_pp + 79);

    auto tr_z_0_z_z_z = pbuffer.data(idx_op_geom_110_pp + 80);

    #pragma omp simd aligned(tr_0_0, tr_0_xz, tr_0_yz, tr_0_zz, tr_z_0_z_z_x, tr_z_0_z_z_y, tr_z_0_z_z_z, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_z_x[i] = -2.0 * tr_0_xz[i] * tke_0 - 6.0 * tr_z_x[i] * tbe_0 + 4.0 * tr_zz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_y[i] = -2.0 * tr_0_yz[i] * tke_0 - 6.0 * tr_z_y[i] * tbe_0 + 4.0 * tr_zz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_z[i] = tr_0_0[i] - 2.0 * tr_0_zz[i] * tke_0 - 6.0 * tr_z_z[i] * tbe_0 - 2.0 * tr_zz_0[i] * tbe_0 + 4.0 * tr_zz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_z[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

