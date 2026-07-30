#include "GeometricalDerivatives010ForDP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_dp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_dp,
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

    // Set up 0-3 components of targeted buffer : DP

    auto tr_0_0_x_xx_x = pbuffer.data(idx_op_geom_010_dp);

    auto tr_0_0_x_xx_y = pbuffer.data(idx_op_geom_010_dp + 1);

    auto tr_0_0_x_xx_z = pbuffer.data(idx_op_geom_010_dp + 2);

    #pragma omp simd aligned(tr_0_0_x_xx_x, tr_0_0_x_xx_y, tr_0_0_x_xx_z, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xx_x[i] = 2.0 * tr_xxx_x[i] * tbe_0 + 2.0 * tr_xx_xx[i] * tke_0 - 2.0 * tr_x_x[i] - tr_xx_0[i];

        tr_0_0_x_xx_y[i] = 2.0 * tr_xxx_y[i] * tbe_0 + 2.0 * tr_xx_xy[i] * tke_0 - 2.0 * tr_x_y[i];

        tr_0_0_x_xx_z[i] = 2.0 * tr_xxx_z[i] * tbe_0 + 2.0 * tr_xx_xz[i] * tke_0 - 2.0 * tr_x_z[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto tr_0_0_x_xy_x = pbuffer.data(idx_op_geom_010_dp + 3);

    auto tr_0_0_x_xy_y = pbuffer.data(idx_op_geom_010_dp + 4);

    auto tr_0_0_x_xy_z = pbuffer.data(idx_op_geom_010_dp + 5);

    #pragma omp simd aligned(tr_0_0_x_xy_x, tr_0_0_x_xy_y, tr_0_0_x_xy_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_y_x, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xy_x[i] = 2.0 * tr_xxy_x[i] * tbe_0 + 2.0 * tr_xy_xx[i] * tke_0 - tr_y_x[i] - tr_xy_0[i];

        tr_0_0_x_xy_y[i] = 2.0 * tr_xxy_y[i] * tbe_0 + 2.0 * tr_xy_xy[i] * tke_0 - tr_y_y[i];

        tr_0_0_x_xy_z[i] = 2.0 * tr_xxy_z[i] * tbe_0 + 2.0 * tr_xy_xz[i] * tke_0 - tr_y_z[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto tr_0_0_x_xz_x = pbuffer.data(idx_op_geom_010_dp + 6);

    auto tr_0_0_x_xz_y = pbuffer.data(idx_op_geom_010_dp + 7);

    auto tr_0_0_x_xz_z = pbuffer.data(idx_op_geom_010_dp + 8);

    #pragma omp simd aligned(tr_0_0_x_xz_x, tr_0_0_x_xz_y, tr_0_0_x_xz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xz_x[i] = 2.0 * tr_xxz_x[i] * tbe_0 + 2.0 * tr_xz_xx[i] * tke_0 - tr_z_x[i] - tr_xz_0[i];

        tr_0_0_x_xz_y[i] = 2.0 * tr_xxz_y[i] * tbe_0 + 2.0 * tr_xz_xy[i] * tke_0 - tr_z_y[i];

        tr_0_0_x_xz_z[i] = 2.0 * tr_xxz_z[i] * tbe_0 + 2.0 * tr_xz_xz[i] * tke_0 - tr_z_z[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto tr_0_0_x_yy_x = pbuffer.data(idx_op_geom_010_dp + 9);

    auto tr_0_0_x_yy_y = pbuffer.data(idx_op_geom_010_dp + 10);

    auto tr_0_0_x_yy_z = pbuffer.data(idx_op_geom_010_dp + 11);

    #pragma omp simd aligned(tr_0_0_x_yy_x, tr_0_0_x_yy_y, tr_0_0_x_yy_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yy_x[i] = 2.0 * tr_xyy_x[i] * tbe_0 + 2.0 * tr_yy_xx[i] * tke_0 - tr_yy_0[i];

        tr_0_0_x_yy_y[i] = 2.0 * tr_xyy_y[i] * tbe_0 + 2.0 * tr_yy_xy[i] * tke_0;

        tr_0_0_x_yy_z[i] = 2.0 * tr_xyy_z[i] * tbe_0 + 2.0 * tr_yy_xz[i] * tke_0;
    }

    // Set up 12-15 components of targeted buffer : DP

    auto tr_0_0_x_yz_x = pbuffer.data(idx_op_geom_010_dp + 12);

    auto tr_0_0_x_yz_y = pbuffer.data(idx_op_geom_010_dp + 13);

    auto tr_0_0_x_yz_z = pbuffer.data(idx_op_geom_010_dp + 14);

    #pragma omp simd aligned(tr_0_0_x_yz_x, tr_0_0_x_yz_y, tr_0_0_x_yz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yz_x[i] = 2.0 * tr_xyz_x[i] * tbe_0 + 2.0 * tr_yz_xx[i] * tke_0 - tr_yz_0[i];

        tr_0_0_x_yz_y[i] = 2.0 * tr_xyz_y[i] * tbe_0 + 2.0 * tr_yz_xy[i] * tke_0;

        tr_0_0_x_yz_z[i] = 2.0 * tr_xyz_z[i] * tbe_0 + 2.0 * tr_yz_xz[i] * tke_0;
    }

    // Set up 15-18 components of targeted buffer : DP

    auto tr_0_0_x_zz_x = pbuffer.data(idx_op_geom_010_dp + 15);

    auto tr_0_0_x_zz_y = pbuffer.data(idx_op_geom_010_dp + 16);

    auto tr_0_0_x_zz_z = pbuffer.data(idx_op_geom_010_dp + 17);

    #pragma omp simd aligned(tr_0_0_x_zz_x, tr_0_0_x_zz_y, tr_0_0_x_zz_z, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zz_x[i] = 2.0 * tr_xzz_x[i] * tbe_0 + 2.0 * tr_zz_xx[i] * tke_0 - tr_zz_0[i];

        tr_0_0_x_zz_y[i] = 2.0 * tr_xzz_y[i] * tbe_0 + 2.0 * tr_zz_xy[i] * tke_0;

        tr_0_0_x_zz_z[i] = 2.0 * tr_xzz_z[i] * tbe_0 + 2.0 * tr_zz_xz[i] * tke_0;
    }

    // Set up 18-21 components of targeted buffer : DP

    auto tr_0_0_y_xx_x = pbuffer.data(idx_op_geom_010_dp + 18);

    auto tr_0_0_y_xx_y = pbuffer.data(idx_op_geom_010_dp + 19);

    auto tr_0_0_y_xx_z = pbuffer.data(idx_op_geom_010_dp + 20);

    #pragma omp simd aligned(tr_0_0_y_xx_x, tr_0_0_y_xx_y, tr_0_0_y_xx_z, tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxy_x, tr_xxy_y, tr_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xx_x[i] = 2.0 * tr_xxy_x[i] * tbe_0 + 2.0 * tr_xx_xy[i] * tke_0;

        tr_0_0_y_xx_y[i] = 2.0 * tr_xxy_y[i] * tbe_0 + 2.0 * tr_xx_yy[i] * tke_0 - tr_xx_0[i];

        tr_0_0_y_xx_z[i] = 2.0 * tr_xxy_z[i] * tbe_0 + 2.0 * tr_xx_yz[i] * tke_0;
    }

    // Set up 21-24 components of targeted buffer : DP

    auto tr_0_0_y_xy_x = pbuffer.data(idx_op_geom_010_dp + 21);

    auto tr_0_0_y_xy_y = pbuffer.data(idx_op_geom_010_dp + 22);

    auto tr_0_0_y_xy_z = pbuffer.data(idx_op_geom_010_dp + 23);

    #pragma omp simd aligned(tr_0_0_y_xy_x, tr_0_0_y_xy_y, tr_0_0_y_xy_z, tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xy_x[i] = 2.0 * tr_xyy_x[i] * tbe_0 + 2.0 * tr_xy_xy[i] * tke_0 - tr_x_x[i];

        tr_0_0_y_xy_y[i] = 2.0 * tr_xyy_y[i] * tbe_0 + 2.0 * tr_xy_yy[i] * tke_0 - tr_x_y[i] - tr_xy_0[i];

        tr_0_0_y_xy_z[i] = 2.0 * tr_xyy_z[i] * tbe_0 + 2.0 * tr_xy_yz[i] * tke_0 - tr_x_z[i];
    }

    // Set up 24-27 components of targeted buffer : DP

    auto tr_0_0_y_xz_x = pbuffer.data(idx_op_geom_010_dp + 24);

    auto tr_0_0_y_xz_y = pbuffer.data(idx_op_geom_010_dp + 25);

    auto tr_0_0_y_xz_z = pbuffer.data(idx_op_geom_010_dp + 26);

    #pragma omp simd aligned(tr_0_0_y_xz_x, tr_0_0_y_xz_y, tr_0_0_y_xz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xz_x[i] = 2.0 * tr_xyz_x[i] * tbe_0 + 2.0 * tr_xz_xy[i] * tke_0;

        tr_0_0_y_xz_y[i] = 2.0 * tr_xyz_y[i] * tbe_0 + 2.0 * tr_xz_yy[i] * tke_0 - tr_xz_0[i];

        tr_0_0_y_xz_z[i] = 2.0 * tr_xyz_z[i] * tbe_0 + 2.0 * tr_xz_yz[i] * tke_0;
    }

    // Set up 27-30 components of targeted buffer : DP

    auto tr_0_0_y_yy_x = pbuffer.data(idx_op_geom_010_dp + 27);

    auto tr_0_0_y_yy_y = pbuffer.data(idx_op_geom_010_dp + 28);

    auto tr_0_0_y_yy_z = pbuffer.data(idx_op_geom_010_dp + 29);

    #pragma omp simd aligned(tr_0_0_y_yy_x, tr_0_0_y_yy_y, tr_0_0_y_yy_z, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yy_x[i] = 2.0 * tr_yyy_x[i] * tbe_0 + 2.0 * tr_yy_xy[i] * tke_0 - 2.0 * tr_y_x[i];

        tr_0_0_y_yy_y[i] = 2.0 * tr_yyy_y[i] * tbe_0 + 2.0 * tr_yy_yy[i] * tke_0 - 2.0 * tr_y_y[i] - tr_yy_0[i];

        tr_0_0_y_yy_z[i] = 2.0 * tr_yyy_z[i] * tbe_0 + 2.0 * tr_yy_yz[i] * tke_0 - 2.0 * tr_y_z[i];
    }

    // Set up 30-33 components of targeted buffer : DP

    auto tr_0_0_y_yz_x = pbuffer.data(idx_op_geom_010_dp + 30);

    auto tr_0_0_y_yz_y = pbuffer.data(idx_op_geom_010_dp + 31);

    auto tr_0_0_y_yz_z = pbuffer.data(idx_op_geom_010_dp + 32);

    #pragma omp simd aligned(tr_0_0_y_yz_x, tr_0_0_y_yz_y, tr_0_0_y_yz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yz_x[i] = 2.0 * tr_yyz_x[i] * tbe_0 + 2.0 * tr_yz_xy[i] * tke_0 - tr_z_x[i];

        tr_0_0_y_yz_y[i] = 2.0 * tr_yyz_y[i] * tbe_0 + 2.0 * tr_yz_yy[i] * tke_0 - tr_z_y[i] - tr_yz_0[i];

        tr_0_0_y_yz_z[i] = 2.0 * tr_yyz_z[i] * tbe_0 + 2.0 * tr_yz_yz[i] * tke_0 - tr_z_z[i];
    }

    // Set up 33-36 components of targeted buffer : DP

    auto tr_0_0_y_zz_x = pbuffer.data(idx_op_geom_010_dp + 33);

    auto tr_0_0_y_zz_y = pbuffer.data(idx_op_geom_010_dp + 34);

    auto tr_0_0_y_zz_z = pbuffer.data(idx_op_geom_010_dp + 35);

    #pragma omp simd aligned(tr_0_0_y_zz_x, tr_0_0_y_zz_y, tr_0_0_y_zz_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zz_x[i] = 2.0 * tr_yzz_x[i] * tbe_0 + 2.0 * tr_zz_xy[i] * tke_0;

        tr_0_0_y_zz_y[i] = 2.0 * tr_yzz_y[i] * tbe_0 + 2.0 * tr_zz_yy[i] * tke_0 - tr_zz_0[i];

        tr_0_0_y_zz_z[i] = 2.0 * tr_yzz_z[i] * tbe_0 + 2.0 * tr_zz_yz[i] * tke_0;
    }

    // Set up 36-39 components of targeted buffer : DP

    auto tr_0_0_z_xx_x = pbuffer.data(idx_op_geom_010_dp + 36);

    auto tr_0_0_z_xx_y = pbuffer.data(idx_op_geom_010_dp + 37);

    auto tr_0_0_z_xx_z = pbuffer.data(idx_op_geom_010_dp + 38);

    #pragma omp simd aligned(tr_0_0_z_xx_x, tr_0_0_z_xx_y, tr_0_0_z_xx_z, tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_y, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xx_x[i] = 2.0 * tr_xxz_x[i] * tbe_0 + 2.0 * tr_xx_xz[i] * tke_0;

        tr_0_0_z_xx_y[i] = 2.0 * tr_xxz_y[i] * tbe_0 + 2.0 * tr_xx_yz[i] * tke_0;

        tr_0_0_z_xx_z[i] = 2.0 * tr_xxz_z[i] * tbe_0 + 2.0 * tr_xx_zz[i] * tke_0 - tr_xx_0[i];
    }

    // Set up 39-42 components of targeted buffer : DP

    auto tr_0_0_z_xy_x = pbuffer.data(idx_op_geom_010_dp + 39);

    auto tr_0_0_z_xy_y = pbuffer.data(idx_op_geom_010_dp + 40);

    auto tr_0_0_z_xy_z = pbuffer.data(idx_op_geom_010_dp + 41);

    #pragma omp simd aligned(tr_0_0_z_xy_x, tr_0_0_z_xy_y, tr_0_0_z_xy_z, tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xy_x[i] = 2.0 * tr_xyz_x[i] * tbe_0 + 2.0 * tr_xy_xz[i] * tke_0;

        tr_0_0_z_xy_y[i] = 2.0 * tr_xyz_y[i] * tbe_0 + 2.0 * tr_xy_yz[i] * tke_0;

        tr_0_0_z_xy_z[i] = 2.0 * tr_xyz_z[i] * tbe_0 + 2.0 * tr_xy_zz[i] * tke_0 - tr_xy_0[i];
    }

    // Set up 42-45 components of targeted buffer : DP

    auto tr_0_0_z_xz_x = pbuffer.data(idx_op_geom_010_dp + 42);

    auto tr_0_0_z_xz_y = pbuffer.data(idx_op_geom_010_dp + 43);

    auto tr_0_0_z_xz_z = pbuffer.data(idx_op_geom_010_dp + 44);

    #pragma omp simd aligned(tr_0_0_z_xz_x, tr_0_0_z_xz_y, tr_0_0_z_xz_z, tr_x_x, tr_x_y, tr_x_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xz_x[i] = 2.0 * tr_xzz_x[i] * tbe_0 + 2.0 * tr_xz_xz[i] * tke_0 - tr_x_x[i];

        tr_0_0_z_xz_y[i] = 2.0 * tr_xzz_y[i] * tbe_0 + 2.0 * tr_xz_yz[i] * tke_0 - tr_x_y[i];

        tr_0_0_z_xz_z[i] = 2.0 * tr_xzz_z[i] * tbe_0 + 2.0 * tr_xz_zz[i] * tke_0 - tr_x_z[i] - tr_xz_0[i];
    }

    // Set up 45-48 components of targeted buffer : DP

    auto tr_0_0_z_yy_x = pbuffer.data(idx_op_geom_010_dp + 45);

    auto tr_0_0_z_yy_y = pbuffer.data(idx_op_geom_010_dp + 46);

    auto tr_0_0_z_yy_z = pbuffer.data(idx_op_geom_010_dp + 47);

    #pragma omp simd aligned(tr_0_0_z_yy_x, tr_0_0_z_yy_y, tr_0_0_z_yy_z, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yy_x[i] = 2.0 * tr_yyz_x[i] * tbe_0 + 2.0 * tr_yy_xz[i] * tke_0;

        tr_0_0_z_yy_y[i] = 2.0 * tr_yyz_y[i] * tbe_0 + 2.0 * tr_yy_yz[i] * tke_0;

        tr_0_0_z_yy_z[i] = 2.0 * tr_yyz_z[i] * tbe_0 + 2.0 * tr_yy_zz[i] * tke_0 - tr_yy_0[i];
    }

    // Set up 48-51 components of targeted buffer : DP

    auto tr_0_0_z_yz_x = pbuffer.data(idx_op_geom_010_dp + 48);

    auto tr_0_0_z_yz_y = pbuffer.data(idx_op_geom_010_dp + 49);

    auto tr_0_0_z_yz_z = pbuffer.data(idx_op_geom_010_dp + 50);

    #pragma omp simd aligned(tr_0_0_z_yz_x, tr_0_0_z_yz_y, tr_0_0_z_yz_z, tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yz_x[i] = 2.0 * tr_yzz_x[i] * tbe_0 + 2.0 * tr_yz_xz[i] * tke_0 - tr_y_x[i];

        tr_0_0_z_yz_y[i] = 2.0 * tr_yzz_y[i] * tbe_0 + 2.0 * tr_yz_yz[i] * tke_0 - tr_y_y[i];

        tr_0_0_z_yz_z[i] = 2.0 * tr_yzz_z[i] * tbe_0 + 2.0 * tr_yz_zz[i] * tke_0 - tr_y_z[i] - tr_yz_0[i];
    }

    // Set up 51-54 components of targeted buffer : DP

    auto tr_0_0_z_zz_x = pbuffer.data(idx_op_geom_010_dp + 51);

    auto tr_0_0_z_zz_y = pbuffer.data(idx_op_geom_010_dp + 52);

    auto tr_0_0_z_zz_z = pbuffer.data(idx_op_geom_010_dp + 53);

    #pragma omp simd aligned(tr_0_0_z_zz_x, tr_0_0_z_zz_y, tr_0_0_z_zz_z, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zz_x[i] = 2.0 * tr_zzz_x[i] * tbe_0 + 2.0 * tr_zz_xz[i] * tke_0 - 2.0 * tr_z_x[i];

        tr_0_0_z_zz_y[i] = 2.0 * tr_zzz_y[i] * tbe_0 + 2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_z_y[i];

        tr_0_0_z_zz_z[i] = 2.0 * tr_zzz_z[i] * tbe_0 + 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_z_z[i] - tr_zz_0[i];
    }

}

} // t2cgeom namespace

