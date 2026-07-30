#include "GeometricalDerivatives010ForPD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_pd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pd,
                         const int idx_op_sd,
                         const int idx_op_pp,
                         const int idx_op_pf,
                         const int idx_op_dd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up 0-6 components of targeted buffer : PD

    auto tr_0_0_x_x_xx = pbuffer.data(idx_op_geom_010_pd);

    auto tr_0_0_x_x_xy = pbuffer.data(idx_op_geom_010_pd + 1);

    auto tr_0_0_x_x_xz = pbuffer.data(idx_op_geom_010_pd + 2);

    auto tr_0_0_x_x_yy = pbuffer.data(idx_op_geom_010_pd + 3);

    auto tr_0_0_x_x_yz = pbuffer.data(idx_op_geom_010_pd + 4);

    auto tr_0_0_x_x_zz = pbuffer.data(idx_op_geom_010_pd + 5);

    #pragma omp simd aligned(tr_0_0_x_x_xx, tr_0_0_x_x_xy, tr_0_0_x_x_xz, tr_0_0_x_x_yy, tr_0_0_x_x_yz, tr_0_0_x_x_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_x, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_z, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_x_xx[i] = 2.0 * tr_xx_xx[i] * tbe_0 + 2.0 * tr_x_xxx[i] * tke_0 - tr_0_xx[i] - 2.0 * tr_x_x[i];

        tr_0_0_x_x_xy[i] = 2.0 * tr_xx_xy[i] * tbe_0 + 2.0 * tr_x_xxy[i] * tke_0 - tr_0_xy[i] - tr_x_y[i];

        tr_0_0_x_x_xz[i] = 2.0 * tr_xx_xz[i] * tbe_0 + 2.0 * tr_x_xxz[i] * tke_0 - tr_0_xz[i] - tr_x_z[i];

        tr_0_0_x_x_yy[i] = 2.0 * tr_xx_yy[i] * tbe_0 + 2.0 * tr_x_xyy[i] * tke_0 - tr_0_yy[i];

        tr_0_0_x_x_yz[i] = 2.0 * tr_xx_yz[i] * tbe_0 + 2.0 * tr_x_xyz[i] * tke_0 - tr_0_yz[i];

        tr_0_0_x_x_zz[i] = 2.0 * tr_xx_zz[i] * tbe_0 + 2.0 * tr_x_xzz[i] * tke_0 - tr_0_zz[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto tr_0_0_x_y_xx = pbuffer.data(idx_op_geom_010_pd + 6);

    auto tr_0_0_x_y_xy = pbuffer.data(idx_op_geom_010_pd + 7);

    auto tr_0_0_x_y_xz = pbuffer.data(idx_op_geom_010_pd + 8);

    auto tr_0_0_x_y_yy = pbuffer.data(idx_op_geom_010_pd + 9);

    auto tr_0_0_x_y_yz = pbuffer.data(idx_op_geom_010_pd + 10);

    auto tr_0_0_x_y_zz = pbuffer.data(idx_op_geom_010_pd + 11);

    #pragma omp simd aligned(tr_0_0_x_y_xx, tr_0_0_x_y_xy, tr_0_0_x_y_xz, tr_0_0_x_y_yy, tr_0_0_x_y_yz, tr_0_0_x_y_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_x, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_y_xx[i] = 2.0 * tr_xy_xx[i] * tbe_0 + 2.0 * tr_y_xxx[i] * tke_0 - 2.0 * tr_y_x[i];

        tr_0_0_x_y_xy[i] = 2.0 * tr_xy_xy[i] * tbe_0 + 2.0 * tr_y_xxy[i] * tke_0 - tr_y_y[i];

        tr_0_0_x_y_xz[i] = 2.0 * tr_xy_xz[i] * tbe_0 + 2.0 * tr_y_xxz[i] * tke_0 - tr_y_z[i];

        tr_0_0_x_y_yy[i] = 2.0 * tr_xy_yy[i] * tbe_0 + 2.0 * tr_y_xyy[i] * tke_0;

        tr_0_0_x_y_yz[i] = 2.0 * tr_xy_yz[i] * tbe_0 + 2.0 * tr_y_xyz[i] * tke_0;

        tr_0_0_x_y_zz[i] = 2.0 * tr_xy_zz[i] * tbe_0 + 2.0 * tr_y_xzz[i] * tke_0;
    }

    // Set up 12-18 components of targeted buffer : PD

    auto tr_0_0_x_z_xx = pbuffer.data(idx_op_geom_010_pd + 12);

    auto tr_0_0_x_z_xy = pbuffer.data(idx_op_geom_010_pd + 13);

    auto tr_0_0_x_z_xz = pbuffer.data(idx_op_geom_010_pd + 14);

    auto tr_0_0_x_z_yy = pbuffer.data(idx_op_geom_010_pd + 15);

    auto tr_0_0_x_z_yz = pbuffer.data(idx_op_geom_010_pd + 16);

    auto tr_0_0_x_z_zz = pbuffer.data(idx_op_geom_010_pd + 17);

    #pragma omp simd aligned(tr_0_0_x_z_xx, tr_0_0_x_z_xy, tr_0_0_x_z_xz, tr_0_0_x_z_yy, tr_0_0_x_z_yz, tr_0_0_x_z_zz, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_x, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_z_xx[i] = 2.0 * tr_xz_xx[i] * tbe_0 + 2.0 * tr_z_xxx[i] * tke_0 - 2.0 * tr_z_x[i];

        tr_0_0_x_z_xy[i] = 2.0 * tr_xz_xy[i] * tbe_0 + 2.0 * tr_z_xxy[i] * tke_0 - tr_z_y[i];

        tr_0_0_x_z_xz[i] = 2.0 * tr_xz_xz[i] * tbe_0 + 2.0 * tr_z_xxz[i] * tke_0 - tr_z_z[i];

        tr_0_0_x_z_yy[i] = 2.0 * tr_xz_yy[i] * tbe_0 + 2.0 * tr_z_xyy[i] * tke_0;

        tr_0_0_x_z_yz[i] = 2.0 * tr_xz_yz[i] * tbe_0 + 2.0 * tr_z_xyz[i] * tke_0;

        tr_0_0_x_z_zz[i] = 2.0 * tr_xz_zz[i] * tbe_0 + 2.0 * tr_z_xzz[i] * tke_0;
    }

    // Set up 18-24 components of targeted buffer : PD

    auto tr_0_0_y_x_xx = pbuffer.data(idx_op_geom_010_pd + 18);

    auto tr_0_0_y_x_xy = pbuffer.data(idx_op_geom_010_pd + 19);

    auto tr_0_0_y_x_xz = pbuffer.data(idx_op_geom_010_pd + 20);

    auto tr_0_0_y_x_yy = pbuffer.data(idx_op_geom_010_pd + 21);

    auto tr_0_0_y_x_yz = pbuffer.data(idx_op_geom_010_pd + 22);

    auto tr_0_0_y_x_zz = pbuffer.data(idx_op_geom_010_pd + 23);

    #pragma omp simd aligned(tr_0_0_y_x_xx, tr_0_0_y_x_xy, tr_0_0_y_x_xz, tr_0_0_y_x_yy, tr_0_0_y_x_yz, tr_0_0_y_x_zz, tr_x_x, tr_x_xxy, tr_x_xyy, tr_x_xyz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_x_xx[i] = 2.0 * tr_xy_xx[i] * tbe_0 + 2.0 * tr_x_xxy[i] * tke_0;

        tr_0_0_y_x_xy[i] = 2.0 * tr_xy_xy[i] * tbe_0 + 2.0 * tr_x_xyy[i] * tke_0 - tr_x_x[i];

        tr_0_0_y_x_xz[i] = 2.0 * tr_xy_xz[i] * tbe_0 + 2.0 * tr_x_xyz[i] * tke_0;

        tr_0_0_y_x_yy[i] = 2.0 * tr_xy_yy[i] * tbe_0 + 2.0 * tr_x_yyy[i] * tke_0 - 2.0 * tr_x_y[i];

        tr_0_0_y_x_yz[i] = 2.0 * tr_xy_yz[i] * tbe_0 + 2.0 * tr_x_yyz[i] * tke_0 - tr_x_z[i];

        tr_0_0_y_x_zz[i] = 2.0 * tr_xy_zz[i] * tbe_0 + 2.0 * tr_x_yzz[i] * tke_0;
    }

    // Set up 24-30 components of targeted buffer : PD

    auto tr_0_0_y_y_xx = pbuffer.data(idx_op_geom_010_pd + 24);

    auto tr_0_0_y_y_xy = pbuffer.data(idx_op_geom_010_pd + 25);

    auto tr_0_0_y_y_xz = pbuffer.data(idx_op_geom_010_pd + 26);

    auto tr_0_0_y_y_yy = pbuffer.data(idx_op_geom_010_pd + 27);

    auto tr_0_0_y_y_yz = pbuffer.data(idx_op_geom_010_pd + 28);

    auto tr_0_0_y_y_zz = pbuffer.data(idx_op_geom_010_pd + 29);

    #pragma omp simd aligned(tr_0_0_y_y_xx, tr_0_0_y_y_xy, tr_0_0_y_y_xz, tr_0_0_y_y_yy, tr_0_0_y_y_yz, tr_0_0_y_y_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_y_x, tr_y_xxy, tr_y_xyy, tr_y_xyz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_y_xx[i] = 2.0 * tr_yy_xx[i] * tbe_0 + 2.0 * tr_y_xxy[i] * tke_0 - tr_0_xx[i];

        tr_0_0_y_y_xy[i] = 2.0 * tr_yy_xy[i] * tbe_0 + 2.0 * tr_y_xyy[i] * tke_0 - tr_0_xy[i] - tr_y_x[i];

        tr_0_0_y_y_xz[i] = 2.0 * tr_yy_xz[i] * tbe_0 + 2.0 * tr_y_xyz[i] * tke_0 - tr_0_xz[i];

        tr_0_0_y_y_yy[i] = 2.0 * tr_yy_yy[i] * tbe_0 + 2.0 * tr_y_yyy[i] * tke_0 - tr_0_yy[i] - 2.0 * tr_y_y[i];

        tr_0_0_y_y_yz[i] = 2.0 * tr_yy_yz[i] * tbe_0 + 2.0 * tr_y_yyz[i] * tke_0 - tr_0_yz[i] - tr_y_z[i];

        tr_0_0_y_y_zz[i] = 2.0 * tr_yy_zz[i] * tbe_0 + 2.0 * tr_y_yzz[i] * tke_0 - tr_0_zz[i];
    }

    // Set up 30-36 components of targeted buffer : PD

    auto tr_0_0_y_z_xx = pbuffer.data(idx_op_geom_010_pd + 30);

    auto tr_0_0_y_z_xy = pbuffer.data(idx_op_geom_010_pd + 31);

    auto tr_0_0_y_z_xz = pbuffer.data(idx_op_geom_010_pd + 32);

    auto tr_0_0_y_z_yy = pbuffer.data(idx_op_geom_010_pd + 33);

    auto tr_0_0_y_z_yz = pbuffer.data(idx_op_geom_010_pd + 34);

    auto tr_0_0_y_z_zz = pbuffer.data(idx_op_geom_010_pd + 35);

    #pragma omp simd aligned(tr_0_0_y_z_xx, tr_0_0_y_z_xy, tr_0_0_y_z_xz, tr_0_0_y_z_yy, tr_0_0_y_z_yz, tr_0_0_y_z_zz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_z_x, tr_z_xxy, tr_z_xyy, tr_z_xyz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_z_xx[i] = 2.0 * tr_yz_xx[i] * tbe_0 + 2.0 * tr_z_xxy[i] * tke_0;

        tr_0_0_y_z_xy[i] = 2.0 * tr_yz_xy[i] * tbe_0 + 2.0 * tr_z_xyy[i] * tke_0 - tr_z_x[i];

        tr_0_0_y_z_xz[i] = 2.0 * tr_yz_xz[i] * tbe_0 + 2.0 * tr_z_xyz[i] * tke_0;

        tr_0_0_y_z_yy[i] = 2.0 * tr_yz_yy[i] * tbe_0 + 2.0 * tr_z_yyy[i] * tke_0 - 2.0 * tr_z_y[i];

        tr_0_0_y_z_yz[i] = 2.0 * tr_yz_yz[i] * tbe_0 + 2.0 * tr_z_yyz[i] * tke_0 - tr_z_z[i];

        tr_0_0_y_z_zz[i] = 2.0 * tr_yz_zz[i] * tbe_0 + 2.0 * tr_z_yzz[i] * tke_0;
    }

    // Set up 36-42 components of targeted buffer : PD

    auto tr_0_0_z_x_xx = pbuffer.data(idx_op_geom_010_pd + 36);

    auto tr_0_0_z_x_xy = pbuffer.data(idx_op_geom_010_pd + 37);

    auto tr_0_0_z_x_xz = pbuffer.data(idx_op_geom_010_pd + 38);

    auto tr_0_0_z_x_yy = pbuffer.data(idx_op_geom_010_pd + 39);

    auto tr_0_0_z_x_yz = pbuffer.data(idx_op_geom_010_pd + 40);

    auto tr_0_0_z_x_zz = pbuffer.data(idx_op_geom_010_pd + 41);

    #pragma omp simd aligned(tr_0_0_z_x_xx, tr_0_0_z_x_xy, tr_0_0_z_x_xz, tr_0_0_z_x_yy, tr_0_0_z_x_yz, tr_0_0_z_x_zz, tr_x_x, tr_x_xxz, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_x_xx[i] = 2.0 * tr_xz_xx[i] * tbe_0 + 2.0 * tr_x_xxz[i] * tke_0;

        tr_0_0_z_x_xy[i] = 2.0 * tr_xz_xy[i] * tbe_0 + 2.0 * tr_x_xyz[i] * tke_0;

        tr_0_0_z_x_xz[i] = 2.0 * tr_xz_xz[i] * tbe_0 + 2.0 * tr_x_xzz[i] * tke_0 - tr_x_x[i];

        tr_0_0_z_x_yy[i] = 2.0 * tr_xz_yy[i] * tbe_0 + 2.0 * tr_x_yyz[i] * tke_0;

        tr_0_0_z_x_yz[i] = 2.0 * tr_xz_yz[i] * tbe_0 + 2.0 * tr_x_yzz[i] * tke_0 - tr_x_y[i];

        tr_0_0_z_x_zz[i] = 2.0 * tr_xz_zz[i] * tbe_0 + 2.0 * tr_x_zzz[i] * tke_0 - 2.0 * tr_x_z[i];
    }

    // Set up 42-48 components of targeted buffer : PD

    auto tr_0_0_z_y_xx = pbuffer.data(idx_op_geom_010_pd + 42);

    auto tr_0_0_z_y_xy = pbuffer.data(idx_op_geom_010_pd + 43);

    auto tr_0_0_z_y_xz = pbuffer.data(idx_op_geom_010_pd + 44);

    auto tr_0_0_z_y_yy = pbuffer.data(idx_op_geom_010_pd + 45);

    auto tr_0_0_z_y_yz = pbuffer.data(idx_op_geom_010_pd + 46);

    auto tr_0_0_z_y_zz = pbuffer.data(idx_op_geom_010_pd + 47);

    #pragma omp simd aligned(tr_0_0_z_y_xx, tr_0_0_z_y_xy, tr_0_0_z_y_xz, tr_0_0_z_y_yy, tr_0_0_z_y_yz, tr_0_0_z_y_zz, tr_y_x, tr_y_xxz, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_y_xx[i] = 2.0 * tr_yz_xx[i] * tbe_0 + 2.0 * tr_y_xxz[i] * tke_0;

        tr_0_0_z_y_xy[i] = 2.0 * tr_yz_xy[i] * tbe_0 + 2.0 * tr_y_xyz[i] * tke_0;

        tr_0_0_z_y_xz[i] = 2.0 * tr_yz_xz[i] * tbe_0 + 2.0 * tr_y_xzz[i] * tke_0 - tr_y_x[i];

        tr_0_0_z_y_yy[i] = 2.0 * tr_yz_yy[i] * tbe_0 + 2.0 * tr_y_yyz[i] * tke_0;

        tr_0_0_z_y_yz[i] = 2.0 * tr_yz_yz[i] * tbe_0 + 2.0 * tr_y_yzz[i] * tke_0 - tr_y_y[i];

        tr_0_0_z_y_zz[i] = 2.0 * tr_yz_zz[i] * tbe_0 + 2.0 * tr_y_zzz[i] * tke_0 - 2.0 * tr_y_z[i];
    }

    // Set up 48-54 components of targeted buffer : PD

    auto tr_0_0_z_z_xx = pbuffer.data(idx_op_geom_010_pd + 48);

    auto tr_0_0_z_z_xy = pbuffer.data(idx_op_geom_010_pd + 49);

    auto tr_0_0_z_z_xz = pbuffer.data(idx_op_geom_010_pd + 50);

    auto tr_0_0_z_z_yy = pbuffer.data(idx_op_geom_010_pd + 51);

    auto tr_0_0_z_z_yz = pbuffer.data(idx_op_geom_010_pd + 52);

    auto tr_0_0_z_z_zz = pbuffer.data(idx_op_geom_010_pd + 53);

    #pragma omp simd aligned(tr_0_0_z_z_xx, tr_0_0_z_z_xy, tr_0_0_z_z_xz, tr_0_0_z_z_yy, tr_0_0_z_z_yz, tr_0_0_z_z_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_z_x, tr_z_xxz, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_z_xx[i] = 2.0 * tr_zz_xx[i] * tbe_0 + 2.0 * tr_z_xxz[i] * tke_0 - tr_0_xx[i];

        tr_0_0_z_z_xy[i] = 2.0 * tr_zz_xy[i] * tbe_0 + 2.0 * tr_z_xyz[i] * tke_0 - tr_0_xy[i];

        tr_0_0_z_z_xz[i] = 2.0 * tr_zz_xz[i] * tbe_0 + 2.0 * tr_z_xzz[i] * tke_0 - tr_0_xz[i] - tr_z_x[i];

        tr_0_0_z_z_yy[i] = 2.0 * tr_zz_yy[i] * tbe_0 + 2.0 * tr_z_yyz[i] * tke_0 - tr_0_yy[i];

        tr_0_0_z_z_yz[i] = 2.0 * tr_zz_yz[i] * tbe_0 + 2.0 * tr_z_yzz[i] * tke_0 - tr_0_yz[i] - tr_z_y[i];

        tr_0_0_z_z_zz[i] = 2.0 * tr_zz_zz[i] * tbe_0 + 2.0 * tr_z_zzz[i] * tke_0 - tr_0_zz[i] - 2.0 * tr_z_z[i];
    }

}

} // t2cgeom namespace

