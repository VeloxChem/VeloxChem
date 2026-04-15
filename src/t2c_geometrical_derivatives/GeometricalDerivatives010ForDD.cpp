#include "GeometricalDerivatives010ForDD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_dd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_dd,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const int idx_op_df,
                         const int idx_op_fd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : PD

    auto tr_x_xx = pbuffer.data(idx_op_pd);

    auto tr_x_xy = pbuffer.data(idx_op_pd + 1);

    auto tr_x_xz = pbuffer.data(idx_op_pd + 2);

    auto tr_x_yy = pbuffer.data(idx_op_pd + 3);

    auto tr_x_yz = pbuffer.data(idx_op_pd + 4);

    auto tr_x_zz = pbuffer.data(idx_op_pd + 5);

    auto tr_y_xx = pbuffer.data(idx_op_pd + 6);

    auto tr_y_xy = pbuffer.data(idx_op_pd + 7);

    auto tr_y_xz = pbuffer.data(idx_op_pd + 8);

    auto tr_y_yy = pbuffer.data(idx_op_pd + 9);

    auto tr_y_yz = pbuffer.data(idx_op_pd + 10);

    auto tr_y_zz = pbuffer.data(idx_op_pd + 11);

    auto tr_z_xx = pbuffer.data(idx_op_pd + 12);

    auto tr_z_xy = pbuffer.data(idx_op_pd + 13);

    auto tr_z_xz = pbuffer.data(idx_op_pd + 14);

    auto tr_z_yy = pbuffer.data(idx_op_pd + 15);

    auto tr_z_yz = pbuffer.data(idx_op_pd + 16);

    auto tr_z_zz = pbuffer.data(idx_op_pd + 17);

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

    // Set up components of auxiliary buffer : DF

    auto tr_xx_xxx = pbuffer.data(idx_op_df);

    auto tr_xx_xxy = pbuffer.data(idx_op_df + 1);

    auto tr_xx_xxz = pbuffer.data(idx_op_df + 2);

    auto tr_xx_xyy = pbuffer.data(idx_op_df + 3);

    auto tr_xx_xyz = pbuffer.data(idx_op_df + 4);

    auto tr_xx_xzz = pbuffer.data(idx_op_df + 5);

    auto tr_xx_yyy = pbuffer.data(idx_op_df + 6);

    auto tr_xx_yyz = pbuffer.data(idx_op_df + 7);

    auto tr_xx_yzz = pbuffer.data(idx_op_df + 8);

    auto tr_xx_zzz = pbuffer.data(idx_op_df + 9);

    auto tr_xy_xxx = pbuffer.data(idx_op_df + 10);

    auto tr_xy_xxy = pbuffer.data(idx_op_df + 11);

    auto tr_xy_xxz = pbuffer.data(idx_op_df + 12);

    auto tr_xy_xyy = pbuffer.data(idx_op_df + 13);

    auto tr_xy_xyz = pbuffer.data(idx_op_df + 14);

    auto tr_xy_xzz = pbuffer.data(idx_op_df + 15);

    auto tr_xy_yyy = pbuffer.data(idx_op_df + 16);

    auto tr_xy_yyz = pbuffer.data(idx_op_df + 17);

    auto tr_xy_yzz = pbuffer.data(idx_op_df + 18);

    auto tr_xy_zzz = pbuffer.data(idx_op_df + 19);

    auto tr_xz_xxx = pbuffer.data(idx_op_df + 20);

    auto tr_xz_xxy = pbuffer.data(idx_op_df + 21);

    auto tr_xz_xxz = pbuffer.data(idx_op_df + 22);

    auto tr_xz_xyy = pbuffer.data(idx_op_df + 23);

    auto tr_xz_xyz = pbuffer.data(idx_op_df + 24);

    auto tr_xz_xzz = pbuffer.data(idx_op_df + 25);

    auto tr_xz_yyy = pbuffer.data(idx_op_df + 26);

    auto tr_xz_yyz = pbuffer.data(idx_op_df + 27);

    auto tr_xz_yzz = pbuffer.data(idx_op_df + 28);

    auto tr_xz_zzz = pbuffer.data(idx_op_df + 29);

    auto tr_yy_xxx = pbuffer.data(idx_op_df + 30);

    auto tr_yy_xxy = pbuffer.data(idx_op_df + 31);

    auto tr_yy_xxz = pbuffer.data(idx_op_df + 32);

    auto tr_yy_xyy = pbuffer.data(idx_op_df + 33);

    auto tr_yy_xyz = pbuffer.data(idx_op_df + 34);

    auto tr_yy_xzz = pbuffer.data(idx_op_df + 35);

    auto tr_yy_yyy = pbuffer.data(idx_op_df + 36);

    auto tr_yy_yyz = pbuffer.data(idx_op_df + 37);

    auto tr_yy_yzz = pbuffer.data(idx_op_df + 38);

    auto tr_yy_zzz = pbuffer.data(idx_op_df + 39);

    auto tr_yz_xxx = pbuffer.data(idx_op_df + 40);

    auto tr_yz_xxy = pbuffer.data(idx_op_df + 41);

    auto tr_yz_xxz = pbuffer.data(idx_op_df + 42);

    auto tr_yz_xyy = pbuffer.data(idx_op_df + 43);

    auto tr_yz_xyz = pbuffer.data(idx_op_df + 44);

    auto tr_yz_xzz = pbuffer.data(idx_op_df + 45);

    auto tr_yz_yyy = pbuffer.data(idx_op_df + 46);

    auto tr_yz_yyz = pbuffer.data(idx_op_df + 47);

    auto tr_yz_yzz = pbuffer.data(idx_op_df + 48);

    auto tr_yz_zzz = pbuffer.data(idx_op_df + 49);

    auto tr_zz_xxx = pbuffer.data(idx_op_df + 50);

    auto tr_zz_xxy = pbuffer.data(idx_op_df + 51);

    auto tr_zz_xxz = pbuffer.data(idx_op_df + 52);

    auto tr_zz_xyy = pbuffer.data(idx_op_df + 53);

    auto tr_zz_xyz = pbuffer.data(idx_op_df + 54);

    auto tr_zz_xzz = pbuffer.data(idx_op_df + 55);

    auto tr_zz_yyy = pbuffer.data(idx_op_df + 56);

    auto tr_zz_yyz = pbuffer.data(idx_op_df + 57);

    auto tr_zz_yzz = pbuffer.data(idx_op_df + 58);

    auto tr_zz_zzz = pbuffer.data(idx_op_df + 59);

    // Set up components of auxiliary buffer : FD

    auto tr_xxx_xx = pbuffer.data(idx_op_fd);

    auto tr_xxx_xy = pbuffer.data(idx_op_fd + 1);

    auto tr_xxx_xz = pbuffer.data(idx_op_fd + 2);

    auto tr_xxx_yy = pbuffer.data(idx_op_fd + 3);

    auto tr_xxx_yz = pbuffer.data(idx_op_fd + 4);

    auto tr_xxx_zz = pbuffer.data(idx_op_fd + 5);

    auto tr_xxy_xx = pbuffer.data(idx_op_fd + 6);

    auto tr_xxy_xy = pbuffer.data(idx_op_fd + 7);

    auto tr_xxy_xz = pbuffer.data(idx_op_fd + 8);

    auto tr_xxy_yy = pbuffer.data(idx_op_fd + 9);

    auto tr_xxy_yz = pbuffer.data(idx_op_fd + 10);

    auto tr_xxy_zz = pbuffer.data(idx_op_fd + 11);

    auto tr_xxz_xx = pbuffer.data(idx_op_fd + 12);

    auto tr_xxz_xy = pbuffer.data(idx_op_fd + 13);

    auto tr_xxz_xz = pbuffer.data(idx_op_fd + 14);

    auto tr_xxz_yy = pbuffer.data(idx_op_fd + 15);

    auto tr_xxz_yz = pbuffer.data(idx_op_fd + 16);

    auto tr_xxz_zz = pbuffer.data(idx_op_fd + 17);

    auto tr_xyy_xx = pbuffer.data(idx_op_fd + 18);

    auto tr_xyy_xy = pbuffer.data(idx_op_fd + 19);

    auto tr_xyy_xz = pbuffer.data(idx_op_fd + 20);

    auto tr_xyy_yy = pbuffer.data(idx_op_fd + 21);

    auto tr_xyy_yz = pbuffer.data(idx_op_fd + 22);

    auto tr_xyy_zz = pbuffer.data(idx_op_fd + 23);

    auto tr_xyz_xx = pbuffer.data(idx_op_fd + 24);

    auto tr_xyz_xy = pbuffer.data(idx_op_fd + 25);

    auto tr_xyz_xz = pbuffer.data(idx_op_fd + 26);

    auto tr_xyz_yy = pbuffer.data(idx_op_fd + 27);

    auto tr_xyz_yz = pbuffer.data(idx_op_fd + 28);

    auto tr_xyz_zz = pbuffer.data(idx_op_fd + 29);

    auto tr_xzz_xx = pbuffer.data(idx_op_fd + 30);

    auto tr_xzz_xy = pbuffer.data(idx_op_fd + 31);

    auto tr_xzz_xz = pbuffer.data(idx_op_fd + 32);

    auto tr_xzz_yy = pbuffer.data(idx_op_fd + 33);

    auto tr_xzz_yz = pbuffer.data(idx_op_fd + 34);

    auto tr_xzz_zz = pbuffer.data(idx_op_fd + 35);

    auto tr_yyy_xx = pbuffer.data(idx_op_fd + 36);

    auto tr_yyy_xy = pbuffer.data(idx_op_fd + 37);

    auto tr_yyy_xz = pbuffer.data(idx_op_fd + 38);

    auto tr_yyy_yy = pbuffer.data(idx_op_fd + 39);

    auto tr_yyy_yz = pbuffer.data(idx_op_fd + 40);

    auto tr_yyy_zz = pbuffer.data(idx_op_fd + 41);

    auto tr_yyz_xx = pbuffer.data(idx_op_fd + 42);

    auto tr_yyz_xy = pbuffer.data(idx_op_fd + 43);

    auto tr_yyz_xz = pbuffer.data(idx_op_fd + 44);

    auto tr_yyz_yy = pbuffer.data(idx_op_fd + 45);

    auto tr_yyz_yz = pbuffer.data(idx_op_fd + 46);

    auto tr_yyz_zz = pbuffer.data(idx_op_fd + 47);

    auto tr_yzz_xx = pbuffer.data(idx_op_fd + 48);

    auto tr_yzz_xy = pbuffer.data(idx_op_fd + 49);

    auto tr_yzz_xz = pbuffer.data(idx_op_fd + 50);

    auto tr_yzz_yy = pbuffer.data(idx_op_fd + 51);

    auto tr_yzz_yz = pbuffer.data(idx_op_fd + 52);

    auto tr_yzz_zz = pbuffer.data(idx_op_fd + 53);

    auto tr_zzz_xx = pbuffer.data(idx_op_fd + 54);

    auto tr_zzz_xy = pbuffer.data(idx_op_fd + 55);

    auto tr_zzz_xz = pbuffer.data(idx_op_fd + 56);

    auto tr_zzz_yy = pbuffer.data(idx_op_fd + 57);

    auto tr_zzz_yz = pbuffer.data(idx_op_fd + 58);

    auto tr_zzz_zz = pbuffer.data(idx_op_fd + 59);

    // Set up 0-6 components of targeted buffer : DD

    auto tr_0_0_x_xx_xx = pbuffer.data(idx_op_geom_010_dd);

    auto tr_0_0_x_xx_xy = pbuffer.data(idx_op_geom_010_dd + 1);

    auto tr_0_0_x_xx_xz = pbuffer.data(idx_op_geom_010_dd + 2);

    auto tr_0_0_x_xx_yy = pbuffer.data(idx_op_geom_010_dd + 3);

    auto tr_0_0_x_xx_yz = pbuffer.data(idx_op_geom_010_dd + 4);

    auto tr_0_0_x_xx_zz = pbuffer.data(idx_op_geom_010_dd + 5);

    #pragma omp simd aligned(tr_0_0_x_xx_xx, tr_0_0_x_xx_xy, tr_0_0_x_xx_xz, tr_0_0_x_xx_yy, tr_0_0_x_xx_yz, tr_0_0_x_xx_zz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_z, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xx_xx[i] = 2.0 * tr_xxx_xx[i] * tbe_0 + 2.0 * tr_xx_xxx[i] * tke_0 - 2.0 * tr_x_xx[i] - 2.0 * tr_xx_x[i];

        tr_0_0_x_xx_xy[i] = 2.0 * tr_xxx_xy[i] * tbe_0 + 2.0 * tr_xx_xxy[i] * tke_0 - 2.0 * tr_x_xy[i] - tr_xx_y[i];

        tr_0_0_x_xx_xz[i] = 2.0 * tr_xxx_xz[i] * tbe_0 + 2.0 * tr_xx_xxz[i] * tke_0 - 2.0 * tr_x_xz[i] - tr_xx_z[i];

        tr_0_0_x_xx_yy[i] = 2.0 * tr_xxx_yy[i] * tbe_0 + 2.0 * tr_xx_xyy[i] * tke_0 - 2.0 * tr_x_yy[i];

        tr_0_0_x_xx_yz[i] = 2.0 * tr_xxx_yz[i] * tbe_0 + 2.0 * tr_xx_xyz[i] * tke_0 - 2.0 * tr_x_yz[i];

        tr_0_0_x_xx_zz[i] = 2.0 * tr_xxx_zz[i] * tbe_0 + 2.0 * tr_xx_xzz[i] * tke_0 - 2.0 * tr_x_zz[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto tr_0_0_x_xy_xx = pbuffer.data(idx_op_geom_010_dd + 6);

    auto tr_0_0_x_xy_xy = pbuffer.data(idx_op_geom_010_dd + 7);

    auto tr_0_0_x_xy_xz = pbuffer.data(idx_op_geom_010_dd + 8);

    auto tr_0_0_x_xy_yy = pbuffer.data(idx_op_geom_010_dd + 9);

    auto tr_0_0_x_xy_yz = pbuffer.data(idx_op_geom_010_dd + 10);

    auto tr_0_0_x_xy_zz = pbuffer.data(idx_op_geom_010_dd + 11);

    #pragma omp simd aligned(tr_0_0_x_xy_xx, tr_0_0_x_xy_xy, tr_0_0_x_xy_xz, tr_0_0_x_xy_yy, tr_0_0_x_xy_yz, tr_0_0_x_xy_zz, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xy_x, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_z, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xy_xx[i] = 2.0 * tr_xxy_xx[i] * tbe_0 + 2.0 * tr_xy_xxx[i] * tke_0 - tr_y_xx[i] - 2.0 * tr_xy_x[i];

        tr_0_0_x_xy_xy[i] = 2.0 * tr_xxy_xy[i] * tbe_0 + 2.0 * tr_xy_xxy[i] * tke_0 - tr_y_xy[i] - tr_xy_y[i];

        tr_0_0_x_xy_xz[i] = 2.0 * tr_xxy_xz[i] * tbe_0 + 2.0 * tr_xy_xxz[i] * tke_0 - tr_y_xz[i] - tr_xy_z[i];

        tr_0_0_x_xy_yy[i] = 2.0 * tr_xxy_yy[i] * tbe_0 + 2.0 * tr_xy_xyy[i] * tke_0 - tr_y_yy[i];

        tr_0_0_x_xy_yz[i] = 2.0 * tr_xxy_yz[i] * tbe_0 + 2.0 * tr_xy_xyz[i] * tke_0 - tr_y_yz[i];

        tr_0_0_x_xy_zz[i] = 2.0 * tr_xxy_zz[i] * tbe_0 + 2.0 * tr_xy_xzz[i] * tke_0 - tr_y_zz[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto tr_0_0_x_xz_xx = pbuffer.data(idx_op_geom_010_dd + 12);

    auto tr_0_0_x_xz_xy = pbuffer.data(idx_op_geom_010_dd + 13);

    auto tr_0_0_x_xz_xz = pbuffer.data(idx_op_geom_010_dd + 14);

    auto tr_0_0_x_xz_yy = pbuffer.data(idx_op_geom_010_dd + 15);

    auto tr_0_0_x_xz_yz = pbuffer.data(idx_op_geom_010_dd + 16);

    auto tr_0_0_x_xz_zz = pbuffer.data(idx_op_geom_010_dd + 17);

    #pragma omp simd aligned(tr_0_0_x_xz_xx, tr_0_0_x_xz_xy, tr_0_0_x_xz_xz, tr_0_0_x_xz_yy, tr_0_0_x_xz_yz, tr_0_0_x_xz_zz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xz_x, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xz_xx[i] = 2.0 * tr_xxz_xx[i] * tbe_0 + 2.0 * tr_xz_xxx[i] * tke_0 - tr_z_xx[i] - 2.0 * tr_xz_x[i];

        tr_0_0_x_xz_xy[i] = 2.0 * tr_xxz_xy[i] * tbe_0 + 2.0 * tr_xz_xxy[i] * tke_0 - tr_z_xy[i] - tr_xz_y[i];

        tr_0_0_x_xz_xz[i] = 2.0 * tr_xxz_xz[i] * tbe_0 + 2.0 * tr_xz_xxz[i] * tke_0 - tr_z_xz[i] - tr_xz_z[i];

        tr_0_0_x_xz_yy[i] = 2.0 * tr_xxz_yy[i] * tbe_0 + 2.0 * tr_xz_xyy[i] * tke_0 - tr_z_yy[i];

        tr_0_0_x_xz_yz[i] = 2.0 * tr_xxz_yz[i] * tbe_0 + 2.0 * tr_xz_xyz[i] * tke_0 - tr_z_yz[i];

        tr_0_0_x_xz_zz[i] = 2.0 * tr_xxz_zz[i] * tbe_0 + 2.0 * tr_xz_xzz[i] * tke_0 - tr_z_zz[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto tr_0_0_x_yy_xx = pbuffer.data(idx_op_geom_010_dd + 18);

    auto tr_0_0_x_yy_xy = pbuffer.data(idx_op_geom_010_dd + 19);

    auto tr_0_0_x_yy_xz = pbuffer.data(idx_op_geom_010_dd + 20);

    auto tr_0_0_x_yy_yy = pbuffer.data(idx_op_geom_010_dd + 21);

    auto tr_0_0_x_yy_yz = pbuffer.data(idx_op_geom_010_dd + 22);

    auto tr_0_0_x_yy_zz = pbuffer.data(idx_op_geom_010_dd + 23);

    #pragma omp simd aligned(tr_0_0_x_yy_xx, tr_0_0_x_yy_xy, tr_0_0_x_yy_xz, tr_0_0_x_yy_yy, tr_0_0_x_yy_yz, tr_0_0_x_yy_zz, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yy_xx[i] = 2.0 * tr_xyy_xx[i] * tbe_0 + 2.0 * tr_yy_xxx[i] * tke_0 - 2.0 * tr_yy_x[i];

        tr_0_0_x_yy_xy[i] = 2.0 * tr_xyy_xy[i] * tbe_0 + 2.0 * tr_yy_xxy[i] * tke_0 - tr_yy_y[i];

        tr_0_0_x_yy_xz[i] = 2.0 * tr_xyy_xz[i] * tbe_0 + 2.0 * tr_yy_xxz[i] * tke_0 - tr_yy_z[i];

        tr_0_0_x_yy_yy[i] = 2.0 * tr_xyy_yy[i] * tbe_0 + 2.0 * tr_yy_xyy[i] * tke_0;

        tr_0_0_x_yy_yz[i] = 2.0 * tr_xyy_yz[i] * tbe_0 + 2.0 * tr_yy_xyz[i] * tke_0;

        tr_0_0_x_yy_zz[i] = 2.0 * tr_xyy_zz[i] * tbe_0 + 2.0 * tr_yy_xzz[i] * tke_0;
    }

    // Set up 24-30 components of targeted buffer : DD

    auto tr_0_0_x_yz_xx = pbuffer.data(idx_op_geom_010_dd + 24);

    auto tr_0_0_x_yz_xy = pbuffer.data(idx_op_geom_010_dd + 25);

    auto tr_0_0_x_yz_xz = pbuffer.data(idx_op_geom_010_dd + 26);

    auto tr_0_0_x_yz_yy = pbuffer.data(idx_op_geom_010_dd + 27);

    auto tr_0_0_x_yz_yz = pbuffer.data(idx_op_geom_010_dd + 28);

    auto tr_0_0_x_yz_zz = pbuffer.data(idx_op_geom_010_dd + 29);

    #pragma omp simd aligned(tr_0_0_x_yz_xx, tr_0_0_x_yz_xy, tr_0_0_x_yz_xz, tr_0_0_x_yz_yy, tr_0_0_x_yz_yz, tr_0_0_x_yz_zz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yz_xx[i] = 2.0 * tr_xyz_xx[i] * tbe_0 + 2.0 * tr_yz_xxx[i] * tke_0 - 2.0 * tr_yz_x[i];

        tr_0_0_x_yz_xy[i] = 2.0 * tr_xyz_xy[i] * tbe_0 + 2.0 * tr_yz_xxy[i] * tke_0 - tr_yz_y[i];

        tr_0_0_x_yz_xz[i] = 2.0 * tr_xyz_xz[i] * tbe_0 + 2.0 * tr_yz_xxz[i] * tke_0 - tr_yz_z[i];

        tr_0_0_x_yz_yy[i] = 2.0 * tr_xyz_yy[i] * tbe_0 + 2.0 * tr_yz_xyy[i] * tke_0;

        tr_0_0_x_yz_yz[i] = 2.0 * tr_xyz_yz[i] * tbe_0 + 2.0 * tr_yz_xyz[i] * tke_0;

        tr_0_0_x_yz_zz[i] = 2.0 * tr_xyz_zz[i] * tbe_0 + 2.0 * tr_yz_xzz[i] * tke_0;
    }

    // Set up 30-36 components of targeted buffer : DD

    auto tr_0_0_x_zz_xx = pbuffer.data(idx_op_geom_010_dd + 30);

    auto tr_0_0_x_zz_xy = pbuffer.data(idx_op_geom_010_dd + 31);

    auto tr_0_0_x_zz_xz = pbuffer.data(idx_op_geom_010_dd + 32);

    auto tr_0_0_x_zz_yy = pbuffer.data(idx_op_geom_010_dd + 33);

    auto tr_0_0_x_zz_yz = pbuffer.data(idx_op_geom_010_dd + 34);

    auto tr_0_0_x_zz_zz = pbuffer.data(idx_op_geom_010_dd + 35);

    #pragma omp simd aligned(tr_0_0_x_zz_xx, tr_0_0_x_zz_xy, tr_0_0_x_zz_xz, tr_0_0_x_zz_yy, tr_0_0_x_zz_yz, tr_0_0_x_zz_zz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zz_xx[i] = 2.0 * tr_xzz_xx[i] * tbe_0 + 2.0 * tr_zz_xxx[i] * tke_0 - 2.0 * tr_zz_x[i];

        tr_0_0_x_zz_xy[i] = 2.0 * tr_xzz_xy[i] * tbe_0 + 2.0 * tr_zz_xxy[i] * tke_0 - tr_zz_y[i];

        tr_0_0_x_zz_xz[i] = 2.0 * tr_xzz_xz[i] * tbe_0 + 2.0 * tr_zz_xxz[i] * tke_0 - tr_zz_z[i];

        tr_0_0_x_zz_yy[i] = 2.0 * tr_xzz_yy[i] * tbe_0 + 2.0 * tr_zz_xyy[i] * tke_0;

        tr_0_0_x_zz_yz[i] = 2.0 * tr_xzz_yz[i] * tbe_0 + 2.0 * tr_zz_xyz[i] * tke_0;

        tr_0_0_x_zz_zz[i] = 2.0 * tr_xzz_zz[i] * tbe_0 + 2.0 * tr_zz_xzz[i] * tke_0;
    }

    // Set up 36-42 components of targeted buffer : DD

    auto tr_0_0_y_xx_xx = pbuffer.data(idx_op_geom_010_dd + 36);

    auto tr_0_0_y_xx_xy = pbuffer.data(idx_op_geom_010_dd + 37);

    auto tr_0_0_y_xx_xz = pbuffer.data(idx_op_geom_010_dd + 38);

    auto tr_0_0_y_xx_yy = pbuffer.data(idx_op_geom_010_dd + 39);

    auto tr_0_0_y_xx_yz = pbuffer.data(idx_op_geom_010_dd + 40);

    auto tr_0_0_y_xx_zz = pbuffer.data(idx_op_geom_010_dd + 41);

    #pragma omp simd aligned(tr_0_0_y_xx_xx, tr_0_0_y_xx_xy, tr_0_0_y_xx_xz, tr_0_0_y_xx_yy, tr_0_0_y_xx_yz, tr_0_0_y_xx_zz, tr_xx_x, tr_xx_xxy, tr_xx_xyy, tr_xx_xyz, tr_xx_y, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xx_xx[i] = 2.0 * tr_xxy_xx[i] * tbe_0 + 2.0 * tr_xx_xxy[i] * tke_0;

        tr_0_0_y_xx_xy[i] = 2.0 * tr_xxy_xy[i] * tbe_0 + 2.0 * tr_xx_xyy[i] * tke_0 - tr_xx_x[i];

        tr_0_0_y_xx_xz[i] = 2.0 * tr_xxy_xz[i] * tbe_0 + 2.0 * tr_xx_xyz[i] * tke_0;

        tr_0_0_y_xx_yy[i] = 2.0 * tr_xxy_yy[i] * tbe_0 + 2.0 * tr_xx_yyy[i] * tke_0 - 2.0 * tr_xx_y[i];

        tr_0_0_y_xx_yz[i] = 2.0 * tr_xxy_yz[i] * tbe_0 + 2.0 * tr_xx_yyz[i] * tke_0 - tr_xx_z[i];

        tr_0_0_y_xx_zz[i] = 2.0 * tr_xxy_zz[i] * tbe_0 + 2.0 * tr_xx_yzz[i] * tke_0;
    }

    // Set up 42-48 components of targeted buffer : DD

    auto tr_0_0_y_xy_xx = pbuffer.data(idx_op_geom_010_dd + 42);

    auto tr_0_0_y_xy_xy = pbuffer.data(idx_op_geom_010_dd + 43);

    auto tr_0_0_y_xy_xz = pbuffer.data(idx_op_geom_010_dd + 44);

    auto tr_0_0_y_xy_yy = pbuffer.data(idx_op_geom_010_dd + 45);

    auto tr_0_0_y_xy_yz = pbuffer.data(idx_op_geom_010_dd + 46);

    auto tr_0_0_y_xy_zz = pbuffer.data(idx_op_geom_010_dd + 47);

    #pragma omp simd aligned(tr_0_0_y_xy_xx, tr_0_0_y_xy_xy, tr_0_0_y_xy_xz, tr_0_0_y_xy_yy, tr_0_0_y_xy_yz, tr_0_0_y_xy_zz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xy_x, tr_xy_xxy, tr_xy_xyy, tr_xy_xyz, tr_xy_y, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xy_xx[i] = 2.0 * tr_xyy_xx[i] * tbe_0 + 2.0 * tr_xy_xxy[i] * tke_0 - tr_x_xx[i];

        tr_0_0_y_xy_xy[i] = 2.0 * tr_xyy_xy[i] * tbe_0 + 2.0 * tr_xy_xyy[i] * tke_0 - tr_x_xy[i] - tr_xy_x[i];

        tr_0_0_y_xy_xz[i] = 2.0 * tr_xyy_xz[i] * tbe_0 + 2.0 * tr_xy_xyz[i] * tke_0 - tr_x_xz[i];

        tr_0_0_y_xy_yy[i] = 2.0 * tr_xyy_yy[i] * tbe_0 + 2.0 * tr_xy_yyy[i] * tke_0 - tr_x_yy[i] - 2.0 * tr_xy_y[i];

        tr_0_0_y_xy_yz[i] = 2.0 * tr_xyy_yz[i] * tbe_0 + 2.0 * tr_xy_yyz[i] * tke_0 - tr_x_yz[i] - tr_xy_z[i];

        tr_0_0_y_xy_zz[i] = 2.0 * tr_xyy_zz[i] * tbe_0 + 2.0 * tr_xy_yzz[i] * tke_0 - tr_x_zz[i];
    }

    // Set up 48-54 components of targeted buffer : DD

    auto tr_0_0_y_xz_xx = pbuffer.data(idx_op_geom_010_dd + 48);

    auto tr_0_0_y_xz_xy = pbuffer.data(idx_op_geom_010_dd + 49);

    auto tr_0_0_y_xz_xz = pbuffer.data(idx_op_geom_010_dd + 50);

    auto tr_0_0_y_xz_yy = pbuffer.data(idx_op_geom_010_dd + 51);

    auto tr_0_0_y_xz_yz = pbuffer.data(idx_op_geom_010_dd + 52);

    auto tr_0_0_y_xz_zz = pbuffer.data(idx_op_geom_010_dd + 53);

    #pragma omp simd aligned(tr_0_0_y_xz_xx, tr_0_0_y_xz_xy, tr_0_0_y_xz_xz, tr_0_0_y_xz_yy, tr_0_0_y_xz_yz, tr_0_0_y_xz_zz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xz_x, tr_xz_xxy, tr_xz_xyy, tr_xz_xyz, tr_xz_y, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xz_xx[i] = 2.0 * tr_xyz_xx[i] * tbe_0 + 2.0 * tr_xz_xxy[i] * tke_0;

        tr_0_0_y_xz_xy[i] = 2.0 * tr_xyz_xy[i] * tbe_0 + 2.0 * tr_xz_xyy[i] * tke_0 - tr_xz_x[i];

        tr_0_0_y_xz_xz[i] = 2.0 * tr_xyz_xz[i] * tbe_0 + 2.0 * tr_xz_xyz[i] * tke_0;

        tr_0_0_y_xz_yy[i] = 2.0 * tr_xyz_yy[i] * tbe_0 + 2.0 * tr_xz_yyy[i] * tke_0 - 2.0 * tr_xz_y[i];

        tr_0_0_y_xz_yz[i] = 2.0 * tr_xyz_yz[i] * tbe_0 + 2.0 * tr_xz_yyz[i] * tke_0 - tr_xz_z[i];

        tr_0_0_y_xz_zz[i] = 2.0 * tr_xyz_zz[i] * tbe_0 + 2.0 * tr_xz_yzz[i] * tke_0;
    }

    // Set up 54-60 components of targeted buffer : DD

    auto tr_0_0_y_yy_xx = pbuffer.data(idx_op_geom_010_dd + 54);

    auto tr_0_0_y_yy_xy = pbuffer.data(idx_op_geom_010_dd + 55);

    auto tr_0_0_y_yy_xz = pbuffer.data(idx_op_geom_010_dd + 56);

    auto tr_0_0_y_yy_yy = pbuffer.data(idx_op_geom_010_dd + 57);

    auto tr_0_0_y_yy_yz = pbuffer.data(idx_op_geom_010_dd + 58);

    auto tr_0_0_y_yy_zz = pbuffer.data(idx_op_geom_010_dd + 59);

    #pragma omp simd aligned(tr_0_0_y_yy_xx, tr_0_0_y_yy_xy, tr_0_0_y_yy_xz, tr_0_0_y_yy_yy, tr_0_0_y_yy_yz, tr_0_0_y_yy_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_xxy, tr_yy_xyy, tr_yy_xyz, tr_yy_y, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yy_xx[i] = 2.0 * tr_yyy_xx[i] * tbe_0 + 2.0 * tr_yy_xxy[i] * tke_0 - 2.0 * tr_y_xx[i];

        tr_0_0_y_yy_xy[i] = 2.0 * tr_yyy_xy[i] * tbe_0 + 2.0 * tr_yy_xyy[i] * tke_0 - 2.0 * tr_y_xy[i] - tr_yy_x[i];

        tr_0_0_y_yy_xz[i] = 2.0 * tr_yyy_xz[i] * tbe_0 + 2.0 * tr_yy_xyz[i] * tke_0 - 2.0 * tr_y_xz[i];

        tr_0_0_y_yy_yy[i] = 2.0 * tr_yyy_yy[i] * tbe_0 + 2.0 * tr_yy_yyy[i] * tke_0 - 2.0 * tr_y_yy[i] - 2.0 * tr_yy_y[i];

        tr_0_0_y_yy_yz[i] = 2.0 * tr_yyy_yz[i] * tbe_0 + 2.0 * tr_yy_yyz[i] * tke_0 - 2.0 * tr_y_yz[i] - tr_yy_z[i];

        tr_0_0_y_yy_zz[i] = 2.0 * tr_yyy_zz[i] * tbe_0 + 2.0 * tr_yy_yzz[i] * tke_0 - 2.0 * tr_y_zz[i];
    }

    // Set up 60-66 components of targeted buffer : DD

    auto tr_0_0_y_yz_xx = pbuffer.data(idx_op_geom_010_dd + 60);

    auto tr_0_0_y_yz_xy = pbuffer.data(idx_op_geom_010_dd + 61);

    auto tr_0_0_y_yz_xz = pbuffer.data(idx_op_geom_010_dd + 62);

    auto tr_0_0_y_yz_yy = pbuffer.data(idx_op_geom_010_dd + 63);

    auto tr_0_0_y_yz_yz = pbuffer.data(idx_op_geom_010_dd + 64);

    auto tr_0_0_y_yz_zz = pbuffer.data(idx_op_geom_010_dd + 65);

    #pragma omp simd aligned(tr_0_0_y_yz_xx, tr_0_0_y_yz_xy, tr_0_0_y_yz_xz, tr_0_0_y_yz_yy, tr_0_0_y_yz_yz, tr_0_0_y_yz_zz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yz_x, tr_yz_xxy, tr_yz_xyy, tr_yz_xyz, tr_yz_y, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yz_xx[i] = 2.0 * tr_yyz_xx[i] * tbe_0 + 2.0 * tr_yz_xxy[i] * tke_0 - tr_z_xx[i];

        tr_0_0_y_yz_xy[i] = 2.0 * tr_yyz_xy[i] * tbe_0 + 2.0 * tr_yz_xyy[i] * tke_0 - tr_z_xy[i] - tr_yz_x[i];

        tr_0_0_y_yz_xz[i] = 2.0 * tr_yyz_xz[i] * tbe_0 + 2.0 * tr_yz_xyz[i] * tke_0 - tr_z_xz[i];

        tr_0_0_y_yz_yy[i] = 2.0 * tr_yyz_yy[i] * tbe_0 + 2.0 * tr_yz_yyy[i] * tke_0 - tr_z_yy[i] - 2.0 * tr_yz_y[i];

        tr_0_0_y_yz_yz[i] = 2.0 * tr_yyz_yz[i] * tbe_0 + 2.0 * tr_yz_yyz[i] * tke_0 - tr_z_yz[i] - tr_yz_z[i];

        tr_0_0_y_yz_zz[i] = 2.0 * tr_yyz_zz[i] * tbe_0 + 2.0 * tr_yz_yzz[i] * tke_0 - tr_z_zz[i];
    }

    // Set up 66-72 components of targeted buffer : DD

    auto tr_0_0_y_zz_xx = pbuffer.data(idx_op_geom_010_dd + 66);

    auto tr_0_0_y_zz_xy = pbuffer.data(idx_op_geom_010_dd + 67);

    auto tr_0_0_y_zz_xz = pbuffer.data(idx_op_geom_010_dd + 68);

    auto tr_0_0_y_zz_yy = pbuffer.data(idx_op_geom_010_dd + 69);

    auto tr_0_0_y_zz_yz = pbuffer.data(idx_op_geom_010_dd + 70);

    auto tr_0_0_y_zz_zz = pbuffer.data(idx_op_geom_010_dd + 71);

    #pragma omp simd aligned(tr_0_0_y_zz_xx, tr_0_0_y_zz_xy, tr_0_0_y_zz_xz, tr_0_0_y_zz_yy, tr_0_0_y_zz_yz, tr_0_0_y_zz_zz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_zz_x, tr_zz_xxy, tr_zz_xyy, tr_zz_xyz, tr_zz_y, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zz_xx[i] = 2.0 * tr_yzz_xx[i] * tbe_0 + 2.0 * tr_zz_xxy[i] * tke_0;

        tr_0_0_y_zz_xy[i] = 2.0 * tr_yzz_xy[i] * tbe_0 + 2.0 * tr_zz_xyy[i] * tke_0 - tr_zz_x[i];

        tr_0_0_y_zz_xz[i] = 2.0 * tr_yzz_xz[i] * tbe_0 + 2.0 * tr_zz_xyz[i] * tke_0;

        tr_0_0_y_zz_yy[i] = 2.0 * tr_yzz_yy[i] * tbe_0 + 2.0 * tr_zz_yyy[i] * tke_0 - 2.0 * tr_zz_y[i];

        tr_0_0_y_zz_yz[i] = 2.0 * tr_yzz_yz[i] * tbe_0 + 2.0 * tr_zz_yyz[i] * tke_0 - tr_zz_z[i];

        tr_0_0_y_zz_zz[i] = 2.0 * tr_yzz_zz[i] * tbe_0 + 2.0 * tr_zz_yzz[i] * tke_0;
    }

    // Set up 72-78 components of targeted buffer : DD

    auto tr_0_0_z_xx_xx = pbuffer.data(idx_op_geom_010_dd + 72);

    auto tr_0_0_z_xx_xy = pbuffer.data(idx_op_geom_010_dd + 73);

    auto tr_0_0_z_xx_xz = pbuffer.data(idx_op_geom_010_dd + 74);

    auto tr_0_0_z_xx_yy = pbuffer.data(idx_op_geom_010_dd + 75);

    auto tr_0_0_z_xx_yz = pbuffer.data(idx_op_geom_010_dd + 76);

    auto tr_0_0_z_xx_zz = pbuffer.data(idx_op_geom_010_dd + 77);

    #pragma omp simd aligned(tr_0_0_z_xx_xx, tr_0_0_z_xx_xy, tr_0_0_z_xx_xz, tr_0_0_z_xx_yy, tr_0_0_z_xx_yz, tr_0_0_z_xx_zz, tr_xx_x, tr_xx_xxz, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xx_xx[i] = 2.0 * tr_xxz_xx[i] * tbe_0 + 2.0 * tr_xx_xxz[i] * tke_0;

        tr_0_0_z_xx_xy[i] = 2.0 * tr_xxz_xy[i] * tbe_0 + 2.0 * tr_xx_xyz[i] * tke_0;

        tr_0_0_z_xx_xz[i] = 2.0 * tr_xxz_xz[i] * tbe_0 + 2.0 * tr_xx_xzz[i] * tke_0 - tr_xx_x[i];

        tr_0_0_z_xx_yy[i] = 2.0 * tr_xxz_yy[i] * tbe_0 + 2.0 * tr_xx_yyz[i] * tke_0;

        tr_0_0_z_xx_yz[i] = 2.0 * tr_xxz_yz[i] * tbe_0 + 2.0 * tr_xx_yzz[i] * tke_0 - tr_xx_y[i];

        tr_0_0_z_xx_zz[i] = 2.0 * tr_xxz_zz[i] * tbe_0 + 2.0 * tr_xx_zzz[i] * tke_0 - 2.0 * tr_xx_z[i];
    }

    // Set up 78-84 components of targeted buffer : DD

    auto tr_0_0_z_xy_xx = pbuffer.data(idx_op_geom_010_dd + 78);

    auto tr_0_0_z_xy_xy = pbuffer.data(idx_op_geom_010_dd + 79);

    auto tr_0_0_z_xy_xz = pbuffer.data(idx_op_geom_010_dd + 80);

    auto tr_0_0_z_xy_yy = pbuffer.data(idx_op_geom_010_dd + 81);

    auto tr_0_0_z_xy_yz = pbuffer.data(idx_op_geom_010_dd + 82);

    auto tr_0_0_z_xy_zz = pbuffer.data(idx_op_geom_010_dd + 83);

    #pragma omp simd aligned(tr_0_0_z_xy_xx, tr_0_0_z_xy_xy, tr_0_0_z_xy_xz, tr_0_0_z_xy_yy, tr_0_0_z_xy_yz, tr_0_0_z_xy_zz, tr_xy_x, tr_xy_xxz, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xy_xx[i] = 2.0 * tr_xyz_xx[i] * tbe_0 + 2.0 * tr_xy_xxz[i] * tke_0;

        tr_0_0_z_xy_xy[i] = 2.0 * tr_xyz_xy[i] * tbe_0 + 2.0 * tr_xy_xyz[i] * tke_0;

        tr_0_0_z_xy_xz[i] = 2.0 * tr_xyz_xz[i] * tbe_0 + 2.0 * tr_xy_xzz[i] * tke_0 - tr_xy_x[i];

        tr_0_0_z_xy_yy[i] = 2.0 * tr_xyz_yy[i] * tbe_0 + 2.0 * tr_xy_yyz[i] * tke_0;

        tr_0_0_z_xy_yz[i] = 2.0 * tr_xyz_yz[i] * tbe_0 + 2.0 * tr_xy_yzz[i] * tke_0 - tr_xy_y[i];

        tr_0_0_z_xy_zz[i] = 2.0 * tr_xyz_zz[i] * tbe_0 + 2.0 * tr_xy_zzz[i] * tke_0 - 2.0 * tr_xy_z[i];
    }

    // Set up 84-90 components of targeted buffer : DD

    auto tr_0_0_z_xz_xx = pbuffer.data(idx_op_geom_010_dd + 84);

    auto tr_0_0_z_xz_xy = pbuffer.data(idx_op_geom_010_dd + 85);

    auto tr_0_0_z_xz_xz = pbuffer.data(idx_op_geom_010_dd + 86);

    auto tr_0_0_z_xz_yy = pbuffer.data(idx_op_geom_010_dd + 87);

    auto tr_0_0_z_xz_yz = pbuffer.data(idx_op_geom_010_dd + 88);

    auto tr_0_0_z_xz_zz = pbuffer.data(idx_op_geom_010_dd + 89);

    #pragma omp simd aligned(tr_0_0_z_xz_xx, tr_0_0_z_xz_xy, tr_0_0_z_xz_xz, tr_0_0_z_xz_yy, tr_0_0_z_xz_yz, tr_0_0_z_xz_zz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xz_x, tr_xz_xxz, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xz_xx[i] = 2.0 * tr_xzz_xx[i] * tbe_0 + 2.0 * tr_xz_xxz[i] * tke_0 - tr_x_xx[i];

        tr_0_0_z_xz_xy[i] = 2.0 * tr_xzz_xy[i] * tbe_0 + 2.0 * tr_xz_xyz[i] * tke_0 - tr_x_xy[i];

        tr_0_0_z_xz_xz[i] = 2.0 * tr_xzz_xz[i] * tbe_0 + 2.0 * tr_xz_xzz[i] * tke_0 - tr_x_xz[i] - tr_xz_x[i];

        tr_0_0_z_xz_yy[i] = 2.0 * tr_xzz_yy[i] * tbe_0 + 2.0 * tr_xz_yyz[i] * tke_0 - tr_x_yy[i];

        tr_0_0_z_xz_yz[i] = 2.0 * tr_xzz_yz[i] * tbe_0 + 2.0 * tr_xz_yzz[i] * tke_0 - tr_x_yz[i] - tr_xz_y[i];

        tr_0_0_z_xz_zz[i] = 2.0 * tr_xzz_zz[i] * tbe_0 + 2.0 * tr_xz_zzz[i] * tke_0 - tr_x_zz[i] - 2.0 * tr_xz_z[i];
    }

    // Set up 90-96 components of targeted buffer : DD

    auto tr_0_0_z_yy_xx = pbuffer.data(idx_op_geom_010_dd + 90);

    auto tr_0_0_z_yy_xy = pbuffer.data(idx_op_geom_010_dd + 91);

    auto tr_0_0_z_yy_xz = pbuffer.data(idx_op_geom_010_dd + 92);

    auto tr_0_0_z_yy_yy = pbuffer.data(idx_op_geom_010_dd + 93);

    auto tr_0_0_z_yy_yz = pbuffer.data(idx_op_geom_010_dd + 94);

    auto tr_0_0_z_yy_zz = pbuffer.data(idx_op_geom_010_dd + 95);

    #pragma omp simd aligned(tr_0_0_z_yy_xx, tr_0_0_z_yy_xy, tr_0_0_z_yy_xz, tr_0_0_z_yy_yy, tr_0_0_z_yy_yz, tr_0_0_z_yy_zz, tr_yy_x, tr_yy_xxz, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yy_xx[i] = 2.0 * tr_yyz_xx[i] * tbe_0 + 2.0 * tr_yy_xxz[i] * tke_0;

        tr_0_0_z_yy_xy[i] = 2.0 * tr_yyz_xy[i] * tbe_0 + 2.0 * tr_yy_xyz[i] * tke_0;

        tr_0_0_z_yy_xz[i] = 2.0 * tr_yyz_xz[i] * tbe_0 + 2.0 * tr_yy_xzz[i] * tke_0 - tr_yy_x[i];

        tr_0_0_z_yy_yy[i] = 2.0 * tr_yyz_yy[i] * tbe_0 + 2.0 * tr_yy_yyz[i] * tke_0;

        tr_0_0_z_yy_yz[i] = 2.0 * tr_yyz_yz[i] * tbe_0 + 2.0 * tr_yy_yzz[i] * tke_0 - tr_yy_y[i];

        tr_0_0_z_yy_zz[i] = 2.0 * tr_yyz_zz[i] * tbe_0 + 2.0 * tr_yy_zzz[i] * tke_0 - 2.0 * tr_yy_z[i];
    }

    // Set up 96-102 components of targeted buffer : DD

    auto tr_0_0_z_yz_xx = pbuffer.data(idx_op_geom_010_dd + 96);

    auto tr_0_0_z_yz_xy = pbuffer.data(idx_op_geom_010_dd + 97);

    auto tr_0_0_z_yz_xz = pbuffer.data(idx_op_geom_010_dd + 98);

    auto tr_0_0_z_yz_yy = pbuffer.data(idx_op_geom_010_dd + 99);

    auto tr_0_0_z_yz_yz = pbuffer.data(idx_op_geom_010_dd + 100);

    auto tr_0_0_z_yz_zz = pbuffer.data(idx_op_geom_010_dd + 101);

    #pragma omp simd aligned(tr_0_0_z_yz_xx, tr_0_0_z_yz_xy, tr_0_0_z_yz_xz, tr_0_0_z_yz_yy, tr_0_0_z_yz_yz, tr_0_0_z_yz_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_xxz, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yz_xx[i] = 2.0 * tr_yzz_xx[i] * tbe_0 + 2.0 * tr_yz_xxz[i] * tke_0 - tr_y_xx[i];

        tr_0_0_z_yz_xy[i] = 2.0 * tr_yzz_xy[i] * tbe_0 + 2.0 * tr_yz_xyz[i] * tke_0 - tr_y_xy[i];

        tr_0_0_z_yz_xz[i] = 2.0 * tr_yzz_xz[i] * tbe_0 + 2.0 * tr_yz_xzz[i] * tke_0 - tr_y_xz[i] - tr_yz_x[i];

        tr_0_0_z_yz_yy[i] = 2.0 * tr_yzz_yy[i] * tbe_0 + 2.0 * tr_yz_yyz[i] * tke_0 - tr_y_yy[i];

        tr_0_0_z_yz_yz[i] = 2.0 * tr_yzz_yz[i] * tbe_0 + 2.0 * tr_yz_yzz[i] * tke_0 - tr_y_yz[i] - tr_yz_y[i];

        tr_0_0_z_yz_zz[i] = 2.0 * tr_yzz_zz[i] * tbe_0 + 2.0 * tr_yz_zzz[i] * tke_0 - tr_y_zz[i] - 2.0 * tr_yz_z[i];
    }

    // Set up 102-108 components of targeted buffer : DD

    auto tr_0_0_z_zz_xx = pbuffer.data(idx_op_geom_010_dd + 102);

    auto tr_0_0_z_zz_xy = pbuffer.data(idx_op_geom_010_dd + 103);

    auto tr_0_0_z_zz_xz = pbuffer.data(idx_op_geom_010_dd + 104);

    auto tr_0_0_z_zz_yy = pbuffer.data(idx_op_geom_010_dd + 105);

    auto tr_0_0_z_zz_yz = pbuffer.data(idx_op_geom_010_dd + 106);

    auto tr_0_0_z_zz_zz = pbuffer.data(idx_op_geom_010_dd + 107);

    #pragma omp simd aligned(tr_0_0_z_zz_xx, tr_0_0_z_zz_xy, tr_0_0_z_zz_xz, tr_0_0_z_zz_yy, tr_0_0_z_zz_yz, tr_0_0_z_zz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_xxz, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zz_xx[i] = 2.0 * tr_zzz_xx[i] * tbe_0 + 2.0 * tr_zz_xxz[i] * tke_0 - 2.0 * tr_z_xx[i];

        tr_0_0_z_zz_xy[i] = 2.0 * tr_zzz_xy[i] * tbe_0 + 2.0 * tr_zz_xyz[i] * tke_0 - 2.0 * tr_z_xy[i];

        tr_0_0_z_zz_xz[i] = 2.0 * tr_zzz_xz[i] * tbe_0 + 2.0 * tr_zz_xzz[i] * tke_0 - 2.0 * tr_z_xz[i] - tr_zz_x[i];

        tr_0_0_z_zz_yy[i] = 2.0 * tr_zzz_yy[i] * tbe_0 + 2.0 * tr_zz_yyz[i] * tke_0 - 2.0 * tr_z_yy[i];

        tr_0_0_z_zz_yz[i] = 2.0 * tr_zzz_yz[i] * tbe_0 + 2.0 * tr_zz_yzz[i] * tke_0 - 2.0 * tr_z_yz[i] - tr_zz_y[i];

        tr_0_0_z_zz_zz[i] = 2.0 * tr_zzz_zz[i] * tbe_0 + 2.0 * tr_zz_zzz[i] * tke_0 - 2.0 * tr_z_zz[i] - 2.0 * tr_zz_z[i];
    }

}

} // t2cgeom namespace

