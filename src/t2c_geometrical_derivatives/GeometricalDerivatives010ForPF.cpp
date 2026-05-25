#include "GeometricalDerivatives010ForPF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_pf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pf,
                         const int idx_op_sf,
                         const int idx_op_pd,
                         const int idx_op_pg,
                         const int idx_op_df,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SF

    auto tr_0_xxx = pbuffer.data(idx_op_sf);

    auto tr_0_xxy = pbuffer.data(idx_op_sf + 1);

    auto tr_0_xxz = pbuffer.data(idx_op_sf + 2);

    auto tr_0_xyy = pbuffer.data(idx_op_sf + 3);

    auto tr_0_xyz = pbuffer.data(idx_op_sf + 4);

    auto tr_0_xzz = pbuffer.data(idx_op_sf + 5);

    auto tr_0_yyy = pbuffer.data(idx_op_sf + 6);

    auto tr_0_yyz = pbuffer.data(idx_op_sf + 7);

    auto tr_0_yzz = pbuffer.data(idx_op_sf + 8);

    auto tr_0_zzz = pbuffer.data(idx_op_sf + 9);

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

    // Set up components of auxiliary buffer : PG

    auto tr_x_xxxx = pbuffer.data(idx_op_pg);

    auto tr_x_xxxy = pbuffer.data(idx_op_pg + 1);

    auto tr_x_xxxz = pbuffer.data(idx_op_pg + 2);

    auto tr_x_xxyy = pbuffer.data(idx_op_pg + 3);

    auto tr_x_xxyz = pbuffer.data(idx_op_pg + 4);

    auto tr_x_xxzz = pbuffer.data(idx_op_pg + 5);

    auto tr_x_xyyy = pbuffer.data(idx_op_pg + 6);

    auto tr_x_xyyz = pbuffer.data(idx_op_pg + 7);

    auto tr_x_xyzz = pbuffer.data(idx_op_pg + 8);

    auto tr_x_xzzz = pbuffer.data(idx_op_pg + 9);

    auto tr_x_yyyy = pbuffer.data(idx_op_pg + 10);

    auto tr_x_yyyz = pbuffer.data(idx_op_pg + 11);

    auto tr_x_yyzz = pbuffer.data(idx_op_pg + 12);

    auto tr_x_yzzz = pbuffer.data(idx_op_pg + 13);

    auto tr_x_zzzz = pbuffer.data(idx_op_pg + 14);

    auto tr_y_xxxx = pbuffer.data(idx_op_pg + 15);

    auto tr_y_xxxy = pbuffer.data(idx_op_pg + 16);

    auto tr_y_xxxz = pbuffer.data(idx_op_pg + 17);

    auto tr_y_xxyy = pbuffer.data(idx_op_pg + 18);

    auto tr_y_xxyz = pbuffer.data(idx_op_pg + 19);

    auto tr_y_xxzz = pbuffer.data(idx_op_pg + 20);

    auto tr_y_xyyy = pbuffer.data(idx_op_pg + 21);

    auto tr_y_xyyz = pbuffer.data(idx_op_pg + 22);

    auto tr_y_xyzz = pbuffer.data(idx_op_pg + 23);

    auto tr_y_xzzz = pbuffer.data(idx_op_pg + 24);

    auto tr_y_yyyy = pbuffer.data(idx_op_pg + 25);

    auto tr_y_yyyz = pbuffer.data(idx_op_pg + 26);

    auto tr_y_yyzz = pbuffer.data(idx_op_pg + 27);

    auto tr_y_yzzz = pbuffer.data(idx_op_pg + 28);

    auto tr_y_zzzz = pbuffer.data(idx_op_pg + 29);

    auto tr_z_xxxx = pbuffer.data(idx_op_pg + 30);

    auto tr_z_xxxy = pbuffer.data(idx_op_pg + 31);

    auto tr_z_xxxz = pbuffer.data(idx_op_pg + 32);

    auto tr_z_xxyy = pbuffer.data(idx_op_pg + 33);

    auto tr_z_xxyz = pbuffer.data(idx_op_pg + 34);

    auto tr_z_xxzz = pbuffer.data(idx_op_pg + 35);

    auto tr_z_xyyy = pbuffer.data(idx_op_pg + 36);

    auto tr_z_xyyz = pbuffer.data(idx_op_pg + 37);

    auto tr_z_xyzz = pbuffer.data(idx_op_pg + 38);

    auto tr_z_xzzz = pbuffer.data(idx_op_pg + 39);

    auto tr_z_yyyy = pbuffer.data(idx_op_pg + 40);

    auto tr_z_yyyz = pbuffer.data(idx_op_pg + 41);

    auto tr_z_yyzz = pbuffer.data(idx_op_pg + 42);

    auto tr_z_yzzz = pbuffer.data(idx_op_pg + 43);

    auto tr_z_zzzz = pbuffer.data(idx_op_pg + 44);

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

    // Set up 0-10 components of targeted buffer : PF

    auto tr_0_0_x_x_xxx = pbuffer.data(idx_op_geom_010_pf);

    auto tr_0_0_x_x_xxy = pbuffer.data(idx_op_geom_010_pf + 1);

    auto tr_0_0_x_x_xxz = pbuffer.data(idx_op_geom_010_pf + 2);

    auto tr_0_0_x_x_xyy = pbuffer.data(idx_op_geom_010_pf + 3);

    auto tr_0_0_x_x_xyz = pbuffer.data(idx_op_geom_010_pf + 4);

    auto tr_0_0_x_x_xzz = pbuffer.data(idx_op_geom_010_pf + 5);

    auto tr_0_0_x_x_yyy = pbuffer.data(idx_op_geom_010_pf + 6);

    auto tr_0_0_x_x_yyz = pbuffer.data(idx_op_geom_010_pf + 7);

    auto tr_0_0_x_x_yzz = pbuffer.data(idx_op_geom_010_pf + 8);

    auto tr_0_0_x_x_zzz = pbuffer.data(idx_op_geom_010_pf + 9);

    #pragma omp simd aligned(tr_0_0_x_x_xxx, tr_0_0_x_x_xxy, tr_0_0_x_x_xxz, tr_0_0_x_x_xyy, tr_0_0_x_x_xyz, tr_0_0_x_x_xzz, tr_0_0_x_x_yyy, tr_0_0_x_x_yyz, tr_0_0_x_x_yzz, tr_0_0_x_x_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xx, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_x_xxx[i] = 2.0 * tr_xx_xxx[i] * tbe_0 + 2.0 * tr_x_xxxx[i] * tke_0 - tr_0_xxx[i] - 3.0 * tr_x_xx[i];

        tr_0_0_x_x_xxy[i] = 2.0 * tr_xx_xxy[i] * tbe_0 + 2.0 * tr_x_xxxy[i] * tke_0 - tr_0_xxy[i] - 2.0 * tr_x_xy[i];

        tr_0_0_x_x_xxz[i] = 2.0 * tr_xx_xxz[i] * tbe_0 + 2.0 * tr_x_xxxz[i] * tke_0 - tr_0_xxz[i] - 2.0 * tr_x_xz[i];

        tr_0_0_x_x_xyy[i] = 2.0 * tr_xx_xyy[i] * tbe_0 + 2.0 * tr_x_xxyy[i] * tke_0 - tr_0_xyy[i] - tr_x_yy[i];

        tr_0_0_x_x_xyz[i] = 2.0 * tr_xx_xyz[i] * tbe_0 + 2.0 * tr_x_xxyz[i] * tke_0 - tr_0_xyz[i] - tr_x_yz[i];

        tr_0_0_x_x_xzz[i] = 2.0 * tr_xx_xzz[i] * tbe_0 + 2.0 * tr_x_xxzz[i] * tke_0 - tr_0_xzz[i] - tr_x_zz[i];

        tr_0_0_x_x_yyy[i] = 2.0 * tr_xx_yyy[i] * tbe_0 + 2.0 * tr_x_xyyy[i] * tke_0 - tr_0_yyy[i];

        tr_0_0_x_x_yyz[i] = 2.0 * tr_xx_yyz[i] * tbe_0 + 2.0 * tr_x_xyyz[i] * tke_0 - tr_0_yyz[i];

        tr_0_0_x_x_yzz[i] = 2.0 * tr_xx_yzz[i] * tbe_0 + 2.0 * tr_x_xyzz[i] * tke_0 - tr_0_yzz[i];

        tr_0_0_x_x_zzz[i] = 2.0 * tr_xx_zzz[i] * tbe_0 + 2.0 * tr_x_xzzz[i] * tke_0 - tr_0_zzz[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto tr_0_0_x_y_xxx = pbuffer.data(idx_op_geom_010_pf + 10);

    auto tr_0_0_x_y_xxy = pbuffer.data(idx_op_geom_010_pf + 11);

    auto tr_0_0_x_y_xxz = pbuffer.data(idx_op_geom_010_pf + 12);

    auto tr_0_0_x_y_xyy = pbuffer.data(idx_op_geom_010_pf + 13);

    auto tr_0_0_x_y_xyz = pbuffer.data(idx_op_geom_010_pf + 14);

    auto tr_0_0_x_y_xzz = pbuffer.data(idx_op_geom_010_pf + 15);

    auto tr_0_0_x_y_yyy = pbuffer.data(idx_op_geom_010_pf + 16);

    auto tr_0_0_x_y_yyz = pbuffer.data(idx_op_geom_010_pf + 17);

    auto tr_0_0_x_y_yzz = pbuffer.data(idx_op_geom_010_pf + 18);

    auto tr_0_0_x_y_zzz = pbuffer.data(idx_op_geom_010_pf + 19);

    #pragma omp simd aligned(tr_0_0_x_y_xxx, tr_0_0_x_y_xxy, tr_0_0_x_y_xxz, tr_0_0_x_y_xyy, tr_0_0_x_y_xyz, tr_0_0_x_y_xzz, tr_0_0_x_y_yyy, tr_0_0_x_y_yyz, tr_0_0_x_y_yzz, tr_0_0_x_y_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_xx, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yz, tr_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_y_xxx[i] = 2.0 * tr_xy_xxx[i] * tbe_0 + 2.0 * tr_y_xxxx[i] * tke_0 - 3.0 * tr_y_xx[i];

        tr_0_0_x_y_xxy[i] = 2.0 * tr_xy_xxy[i] * tbe_0 + 2.0 * tr_y_xxxy[i] * tke_0 - 2.0 * tr_y_xy[i];

        tr_0_0_x_y_xxz[i] = 2.0 * tr_xy_xxz[i] * tbe_0 + 2.0 * tr_y_xxxz[i] * tke_0 - 2.0 * tr_y_xz[i];

        tr_0_0_x_y_xyy[i] = 2.0 * tr_xy_xyy[i] * tbe_0 + 2.0 * tr_y_xxyy[i] * tke_0 - tr_y_yy[i];

        tr_0_0_x_y_xyz[i] = 2.0 * tr_xy_xyz[i] * tbe_0 + 2.0 * tr_y_xxyz[i] * tke_0 - tr_y_yz[i];

        tr_0_0_x_y_xzz[i] = 2.0 * tr_xy_xzz[i] * tbe_0 + 2.0 * tr_y_xxzz[i] * tke_0 - tr_y_zz[i];

        tr_0_0_x_y_yyy[i] = 2.0 * tr_xy_yyy[i] * tbe_0 + 2.0 * tr_y_xyyy[i] * tke_0;

        tr_0_0_x_y_yyz[i] = 2.0 * tr_xy_yyz[i] * tbe_0 + 2.0 * tr_y_xyyz[i] * tke_0;

        tr_0_0_x_y_yzz[i] = 2.0 * tr_xy_yzz[i] * tbe_0 + 2.0 * tr_y_xyzz[i] * tke_0;

        tr_0_0_x_y_zzz[i] = 2.0 * tr_xy_zzz[i] * tbe_0 + 2.0 * tr_y_xzzz[i] * tke_0;
    }

    // Set up 20-30 components of targeted buffer : PF

    auto tr_0_0_x_z_xxx = pbuffer.data(idx_op_geom_010_pf + 20);

    auto tr_0_0_x_z_xxy = pbuffer.data(idx_op_geom_010_pf + 21);

    auto tr_0_0_x_z_xxz = pbuffer.data(idx_op_geom_010_pf + 22);

    auto tr_0_0_x_z_xyy = pbuffer.data(idx_op_geom_010_pf + 23);

    auto tr_0_0_x_z_xyz = pbuffer.data(idx_op_geom_010_pf + 24);

    auto tr_0_0_x_z_xzz = pbuffer.data(idx_op_geom_010_pf + 25);

    auto tr_0_0_x_z_yyy = pbuffer.data(idx_op_geom_010_pf + 26);

    auto tr_0_0_x_z_yyz = pbuffer.data(idx_op_geom_010_pf + 27);

    auto tr_0_0_x_z_yzz = pbuffer.data(idx_op_geom_010_pf + 28);

    auto tr_0_0_x_z_zzz = pbuffer.data(idx_op_geom_010_pf + 29);

    #pragma omp simd aligned(tr_0_0_x_z_xxx, tr_0_0_x_z_xxy, tr_0_0_x_z_xxz, tr_0_0_x_z_xyy, tr_0_0_x_z_xyz, tr_0_0_x_z_xzz, tr_0_0_x_z_yyy, tr_0_0_x_z_yyz, tr_0_0_x_z_yzz, tr_0_0_x_z_zzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_xx, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_z_xxx[i] = 2.0 * tr_xz_xxx[i] * tbe_0 + 2.0 * tr_z_xxxx[i] * tke_0 - 3.0 * tr_z_xx[i];

        tr_0_0_x_z_xxy[i] = 2.0 * tr_xz_xxy[i] * tbe_0 + 2.0 * tr_z_xxxy[i] * tke_0 - 2.0 * tr_z_xy[i];

        tr_0_0_x_z_xxz[i] = 2.0 * tr_xz_xxz[i] * tbe_0 + 2.0 * tr_z_xxxz[i] * tke_0 - 2.0 * tr_z_xz[i];

        tr_0_0_x_z_xyy[i] = 2.0 * tr_xz_xyy[i] * tbe_0 + 2.0 * tr_z_xxyy[i] * tke_0 - tr_z_yy[i];

        tr_0_0_x_z_xyz[i] = 2.0 * tr_xz_xyz[i] * tbe_0 + 2.0 * tr_z_xxyz[i] * tke_0 - tr_z_yz[i];

        tr_0_0_x_z_xzz[i] = 2.0 * tr_xz_xzz[i] * tbe_0 + 2.0 * tr_z_xxzz[i] * tke_0 - tr_z_zz[i];

        tr_0_0_x_z_yyy[i] = 2.0 * tr_xz_yyy[i] * tbe_0 + 2.0 * tr_z_xyyy[i] * tke_0;

        tr_0_0_x_z_yyz[i] = 2.0 * tr_xz_yyz[i] * tbe_0 + 2.0 * tr_z_xyyz[i] * tke_0;

        tr_0_0_x_z_yzz[i] = 2.0 * tr_xz_yzz[i] * tbe_0 + 2.0 * tr_z_xyzz[i] * tke_0;

        tr_0_0_x_z_zzz[i] = 2.0 * tr_xz_zzz[i] * tbe_0 + 2.0 * tr_z_xzzz[i] * tke_0;
    }

    // Set up 30-40 components of targeted buffer : PF

    auto tr_0_0_y_x_xxx = pbuffer.data(idx_op_geom_010_pf + 30);

    auto tr_0_0_y_x_xxy = pbuffer.data(idx_op_geom_010_pf + 31);

    auto tr_0_0_y_x_xxz = pbuffer.data(idx_op_geom_010_pf + 32);

    auto tr_0_0_y_x_xyy = pbuffer.data(idx_op_geom_010_pf + 33);

    auto tr_0_0_y_x_xyz = pbuffer.data(idx_op_geom_010_pf + 34);

    auto tr_0_0_y_x_xzz = pbuffer.data(idx_op_geom_010_pf + 35);

    auto tr_0_0_y_x_yyy = pbuffer.data(idx_op_geom_010_pf + 36);

    auto tr_0_0_y_x_yyz = pbuffer.data(idx_op_geom_010_pf + 37);

    auto tr_0_0_y_x_yzz = pbuffer.data(idx_op_geom_010_pf + 38);

    auto tr_0_0_y_x_zzz = pbuffer.data(idx_op_geom_010_pf + 39);

    #pragma omp simd aligned(tr_0_0_y_x_xxx, tr_0_0_y_x_xxy, tr_0_0_y_x_xxz, tr_0_0_y_x_xyy, tr_0_0_y_x_xyz, tr_0_0_y_x_xzz, tr_0_0_y_x_yyy, tr_0_0_y_x_yyz, tr_0_0_y_x_yzz, tr_0_0_y_x_zzz, tr_x_xx, tr_x_xxxy, tr_x_xxyy, tr_x_xxyz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_yy, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_x_xxx[i] = 2.0 * tr_xy_xxx[i] * tbe_0 + 2.0 * tr_x_xxxy[i] * tke_0;

        tr_0_0_y_x_xxy[i] = 2.0 * tr_xy_xxy[i] * tbe_0 + 2.0 * tr_x_xxyy[i] * tke_0 - tr_x_xx[i];

        tr_0_0_y_x_xxz[i] = 2.0 * tr_xy_xxz[i] * tbe_0 + 2.0 * tr_x_xxyz[i] * tke_0;

        tr_0_0_y_x_xyy[i] = 2.0 * tr_xy_xyy[i] * tbe_0 + 2.0 * tr_x_xyyy[i] * tke_0 - 2.0 * tr_x_xy[i];

        tr_0_0_y_x_xyz[i] = 2.0 * tr_xy_xyz[i] * tbe_0 + 2.0 * tr_x_xyyz[i] * tke_0 - tr_x_xz[i];

        tr_0_0_y_x_xzz[i] = 2.0 * tr_xy_xzz[i] * tbe_0 + 2.0 * tr_x_xyzz[i] * tke_0;

        tr_0_0_y_x_yyy[i] = 2.0 * tr_xy_yyy[i] * tbe_0 + 2.0 * tr_x_yyyy[i] * tke_0 - 3.0 * tr_x_yy[i];

        tr_0_0_y_x_yyz[i] = 2.0 * tr_xy_yyz[i] * tbe_0 + 2.0 * tr_x_yyyz[i] * tke_0 - 2.0 * tr_x_yz[i];

        tr_0_0_y_x_yzz[i] = 2.0 * tr_xy_yzz[i] * tbe_0 + 2.0 * tr_x_yyzz[i] * tke_0 - tr_x_zz[i];

        tr_0_0_y_x_zzz[i] = 2.0 * tr_xy_zzz[i] * tbe_0 + 2.0 * tr_x_yzzz[i] * tke_0;
    }

    // Set up 40-50 components of targeted buffer : PF

    auto tr_0_0_y_y_xxx = pbuffer.data(idx_op_geom_010_pf + 40);

    auto tr_0_0_y_y_xxy = pbuffer.data(idx_op_geom_010_pf + 41);

    auto tr_0_0_y_y_xxz = pbuffer.data(idx_op_geom_010_pf + 42);

    auto tr_0_0_y_y_xyy = pbuffer.data(idx_op_geom_010_pf + 43);

    auto tr_0_0_y_y_xyz = pbuffer.data(idx_op_geom_010_pf + 44);

    auto tr_0_0_y_y_xzz = pbuffer.data(idx_op_geom_010_pf + 45);

    auto tr_0_0_y_y_yyy = pbuffer.data(idx_op_geom_010_pf + 46);

    auto tr_0_0_y_y_yyz = pbuffer.data(idx_op_geom_010_pf + 47);

    auto tr_0_0_y_y_yzz = pbuffer.data(idx_op_geom_010_pf + 48);

    auto tr_0_0_y_y_zzz = pbuffer.data(idx_op_geom_010_pf + 49);

    #pragma omp simd aligned(tr_0_0_y_y_xxx, tr_0_0_y_y_xxy, tr_0_0_y_y_xxz, tr_0_0_y_y_xyy, tr_0_0_y_y_xyz, tr_0_0_y_y_xzz, tr_0_0_y_y_yyy, tr_0_0_y_y_yyz, tr_0_0_y_y_yzz, tr_0_0_y_y_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_y_xx, tr_y_xxxy, tr_y_xxyy, tr_y_xxyz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_yy, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_y_xxx[i] = 2.0 * tr_yy_xxx[i] * tbe_0 + 2.0 * tr_y_xxxy[i] * tke_0 - tr_0_xxx[i];

        tr_0_0_y_y_xxy[i] = 2.0 * tr_yy_xxy[i] * tbe_0 + 2.0 * tr_y_xxyy[i] * tke_0 - tr_0_xxy[i] - tr_y_xx[i];

        tr_0_0_y_y_xxz[i] = 2.0 * tr_yy_xxz[i] * tbe_0 + 2.0 * tr_y_xxyz[i] * tke_0 - tr_0_xxz[i];

        tr_0_0_y_y_xyy[i] = 2.0 * tr_yy_xyy[i] * tbe_0 + 2.0 * tr_y_xyyy[i] * tke_0 - tr_0_xyy[i] - 2.0 * tr_y_xy[i];

        tr_0_0_y_y_xyz[i] = 2.0 * tr_yy_xyz[i] * tbe_0 + 2.0 * tr_y_xyyz[i] * tke_0 - tr_0_xyz[i] - tr_y_xz[i];

        tr_0_0_y_y_xzz[i] = 2.0 * tr_yy_xzz[i] * tbe_0 + 2.0 * tr_y_xyzz[i] * tke_0 - tr_0_xzz[i];

        tr_0_0_y_y_yyy[i] = 2.0 * tr_yy_yyy[i] * tbe_0 + 2.0 * tr_y_yyyy[i] * tke_0 - tr_0_yyy[i] - 3.0 * tr_y_yy[i];

        tr_0_0_y_y_yyz[i] = 2.0 * tr_yy_yyz[i] * tbe_0 + 2.0 * tr_y_yyyz[i] * tke_0 - tr_0_yyz[i] - 2.0 * tr_y_yz[i];

        tr_0_0_y_y_yzz[i] = 2.0 * tr_yy_yzz[i] * tbe_0 + 2.0 * tr_y_yyzz[i] * tke_0 - tr_0_yzz[i] - tr_y_zz[i];

        tr_0_0_y_y_zzz[i] = 2.0 * tr_yy_zzz[i] * tbe_0 + 2.0 * tr_y_yzzz[i] * tke_0 - tr_0_zzz[i];
    }

    // Set up 50-60 components of targeted buffer : PF

    auto tr_0_0_y_z_xxx = pbuffer.data(idx_op_geom_010_pf + 50);

    auto tr_0_0_y_z_xxy = pbuffer.data(idx_op_geom_010_pf + 51);

    auto tr_0_0_y_z_xxz = pbuffer.data(idx_op_geom_010_pf + 52);

    auto tr_0_0_y_z_xyy = pbuffer.data(idx_op_geom_010_pf + 53);

    auto tr_0_0_y_z_xyz = pbuffer.data(idx_op_geom_010_pf + 54);

    auto tr_0_0_y_z_xzz = pbuffer.data(idx_op_geom_010_pf + 55);

    auto tr_0_0_y_z_yyy = pbuffer.data(idx_op_geom_010_pf + 56);

    auto tr_0_0_y_z_yyz = pbuffer.data(idx_op_geom_010_pf + 57);

    auto tr_0_0_y_z_yzz = pbuffer.data(idx_op_geom_010_pf + 58);

    auto tr_0_0_y_z_zzz = pbuffer.data(idx_op_geom_010_pf + 59);

    #pragma omp simd aligned(tr_0_0_y_z_xxx, tr_0_0_y_z_xxy, tr_0_0_y_z_xxz, tr_0_0_y_z_xyy, tr_0_0_y_z_xyz, tr_0_0_y_z_xzz, tr_0_0_y_z_yyy, tr_0_0_y_z_yyz, tr_0_0_y_z_yzz, tr_0_0_y_z_zzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_z_xx, tr_z_xxxy, tr_z_xxyy, tr_z_xxyz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_yy, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_z_xxx[i] = 2.0 * tr_yz_xxx[i] * tbe_0 + 2.0 * tr_z_xxxy[i] * tke_0;

        tr_0_0_y_z_xxy[i] = 2.0 * tr_yz_xxy[i] * tbe_0 + 2.0 * tr_z_xxyy[i] * tke_0 - tr_z_xx[i];

        tr_0_0_y_z_xxz[i] = 2.0 * tr_yz_xxz[i] * tbe_0 + 2.0 * tr_z_xxyz[i] * tke_0;

        tr_0_0_y_z_xyy[i] = 2.0 * tr_yz_xyy[i] * tbe_0 + 2.0 * tr_z_xyyy[i] * tke_0 - 2.0 * tr_z_xy[i];

        tr_0_0_y_z_xyz[i] = 2.0 * tr_yz_xyz[i] * tbe_0 + 2.0 * tr_z_xyyz[i] * tke_0 - tr_z_xz[i];

        tr_0_0_y_z_xzz[i] = 2.0 * tr_yz_xzz[i] * tbe_0 + 2.0 * tr_z_xyzz[i] * tke_0;

        tr_0_0_y_z_yyy[i] = 2.0 * tr_yz_yyy[i] * tbe_0 + 2.0 * tr_z_yyyy[i] * tke_0 - 3.0 * tr_z_yy[i];

        tr_0_0_y_z_yyz[i] = 2.0 * tr_yz_yyz[i] * tbe_0 + 2.0 * tr_z_yyyz[i] * tke_0 - 2.0 * tr_z_yz[i];

        tr_0_0_y_z_yzz[i] = 2.0 * tr_yz_yzz[i] * tbe_0 + 2.0 * tr_z_yyzz[i] * tke_0 - tr_z_zz[i];

        tr_0_0_y_z_zzz[i] = 2.0 * tr_yz_zzz[i] * tbe_0 + 2.0 * tr_z_yzzz[i] * tke_0;
    }

    // Set up 60-70 components of targeted buffer : PF

    auto tr_0_0_z_x_xxx = pbuffer.data(idx_op_geom_010_pf + 60);

    auto tr_0_0_z_x_xxy = pbuffer.data(idx_op_geom_010_pf + 61);

    auto tr_0_0_z_x_xxz = pbuffer.data(idx_op_geom_010_pf + 62);

    auto tr_0_0_z_x_xyy = pbuffer.data(idx_op_geom_010_pf + 63);

    auto tr_0_0_z_x_xyz = pbuffer.data(idx_op_geom_010_pf + 64);

    auto tr_0_0_z_x_xzz = pbuffer.data(idx_op_geom_010_pf + 65);

    auto tr_0_0_z_x_yyy = pbuffer.data(idx_op_geom_010_pf + 66);

    auto tr_0_0_z_x_yyz = pbuffer.data(idx_op_geom_010_pf + 67);

    auto tr_0_0_z_x_yzz = pbuffer.data(idx_op_geom_010_pf + 68);

    auto tr_0_0_z_x_zzz = pbuffer.data(idx_op_geom_010_pf + 69);

    #pragma omp simd aligned(tr_0_0_z_x_xxx, tr_0_0_z_x_xxy, tr_0_0_z_x_xxz, tr_0_0_z_x_xyy, tr_0_0_z_x_xyz, tr_0_0_z_x_xzz, tr_0_0_z_x_yyy, tr_0_0_z_x_yyz, tr_0_0_z_x_yzz, tr_0_0_z_x_zzz, tr_x_xx, tr_x_xxxz, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_x_xxx[i] = 2.0 * tr_xz_xxx[i] * tbe_0 + 2.0 * tr_x_xxxz[i] * tke_0;

        tr_0_0_z_x_xxy[i] = 2.0 * tr_xz_xxy[i] * tbe_0 + 2.0 * tr_x_xxyz[i] * tke_0;

        tr_0_0_z_x_xxz[i] = 2.0 * tr_xz_xxz[i] * tbe_0 + 2.0 * tr_x_xxzz[i] * tke_0 - tr_x_xx[i];

        tr_0_0_z_x_xyy[i] = 2.0 * tr_xz_xyy[i] * tbe_0 + 2.0 * tr_x_xyyz[i] * tke_0;

        tr_0_0_z_x_xyz[i] = 2.0 * tr_xz_xyz[i] * tbe_0 + 2.0 * tr_x_xyzz[i] * tke_0 - tr_x_xy[i];

        tr_0_0_z_x_xzz[i] = 2.0 * tr_xz_xzz[i] * tbe_0 + 2.0 * tr_x_xzzz[i] * tke_0 - 2.0 * tr_x_xz[i];

        tr_0_0_z_x_yyy[i] = 2.0 * tr_xz_yyy[i] * tbe_0 + 2.0 * tr_x_yyyz[i] * tke_0;

        tr_0_0_z_x_yyz[i] = 2.0 * tr_xz_yyz[i] * tbe_0 + 2.0 * tr_x_yyzz[i] * tke_0 - tr_x_yy[i];

        tr_0_0_z_x_yzz[i] = 2.0 * tr_xz_yzz[i] * tbe_0 + 2.0 * tr_x_yzzz[i] * tke_0 - 2.0 * tr_x_yz[i];

        tr_0_0_z_x_zzz[i] = 2.0 * tr_xz_zzz[i] * tbe_0 + 2.0 * tr_x_zzzz[i] * tke_0 - 3.0 * tr_x_zz[i];
    }

    // Set up 70-80 components of targeted buffer : PF

    auto tr_0_0_z_y_xxx = pbuffer.data(idx_op_geom_010_pf + 70);

    auto tr_0_0_z_y_xxy = pbuffer.data(idx_op_geom_010_pf + 71);

    auto tr_0_0_z_y_xxz = pbuffer.data(idx_op_geom_010_pf + 72);

    auto tr_0_0_z_y_xyy = pbuffer.data(idx_op_geom_010_pf + 73);

    auto tr_0_0_z_y_xyz = pbuffer.data(idx_op_geom_010_pf + 74);

    auto tr_0_0_z_y_xzz = pbuffer.data(idx_op_geom_010_pf + 75);

    auto tr_0_0_z_y_yyy = pbuffer.data(idx_op_geom_010_pf + 76);

    auto tr_0_0_z_y_yyz = pbuffer.data(idx_op_geom_010_pf + 77);

    auto tr_0_0_z_y_yzz = pbuffer.data(idx_op_geom_010_pf + 78);

    auto tr_0_0_z_y_zzz = pbuffer.data(idx_op_geom_010_pf + 79);

    #pragma omp simd aligned(tr_0_0_z_y_xxx, tr_0_0_z_y_xxy, tr_0_0_z_y_xxz, tr_0_0_z_y_xyy, tr_0_0_z_y_xyz, tr_0_0_z_y_xzz, tr_0_0_z_y_yyy, tr_0_0_z_y_yyz, tr_0_0_z_y_yzz, tr_0_0_z_y_zzz, tr_y_xx, tr_y_xxxz, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_y_xxx[i] = 2.0 * tr_yz_xxx[i] * tbe_0 + 2.0 * tr_y_xxxz[i] * tke_0;

        tr_0_0_z_y_xxy[i] = 2.0 * tr_yz_xxy[i] * tbe_0 + 2.0 * tr_y_xxyz[i] * tke_0;

        tr_0_0_z_y_xxz[i] = 2.0 * tr_yz_xxz[i] * tbe_0 + 2.0 * tr_y_xxzz[i] * tke_0 - tr_y_xx[i];

        tr_0_0_z_y_xyy[i] = 2.0 * tr_yz_xyy[i] * tbe_0 + 2.0 * tr_y_xyyz[i] * tke_0;

        tr_0_0_z_y_xyz[i] = 2.0 * tr_yz_xyz[i] * tbe_0 + 2.0 * tr_y_xyzz[i] * tke_0 - tr_y_xy[i];

        tr_0_0_z_y_xzz[i] = 2.0 * tr_yz_xzz[i] * tbe_0 + 2.0 * tr_y_xzzz[i] * tke_0 - 2.0 * tr_y_xz[i];

        tr_0_0_z_y_yyy[i] = 2.0 * tr_yz_yyy[i] * tbe_0 + 2.0 * tr_y_yyyz[i] * tke_0;

        tr_0_0_z_y_yyz[i] = 2.0 * tr_yz_yyz[i] * tbe_0 + 2.0 * tr_y_yyzz[i] * tke_0 - tr_y_yy[i];

        tr_0_0_z_y_yzz[i] = 2.0 * tr_yz_yzz[i] * tbe_0 + 2.0 * tr_y_yzzz[i] * tke_0 - 2.0 * tr_y_yz[i];

        tr_0_0_z_y_zzz[i] = 2.0 * tr_yz_zzz[i] * tbe_0 + 2.0 * tr_y_zzzz[i] * tke_0 - 3.0 * tr_y_zz[i];
    }

    // Set up 80-90 components of targeted buffer : PF

    auto tr_0_0_z_z_xxx = pbuffer.data(idx_op_geom_010_pf + 80);

    auto tr_0_0_z_z_xxy = pbuffer.data(idx_op_geom_010_pf + 81);

    auto tr_0_0_z_z_xxz = pbuffer.data(idx_op_geom_010_pf + 82);

    auto tr_0_0_z_z_xyy = pbuffer.data(idx_op_geom_010_pf + 83);

    auto tr_0_0_z_z_xyz = pbuffer.data(idx_op_geom_010_pf + 84);

    auto tr_0_0_z_z_xzz = pbuffer.data(idx_op_geom_010_pf + 85);

    auto tr_0_0_z_z_yyy = pbuffer.data(idx_op_geom_010_pf + 86);

    auto tr_0_0_z_z_yyz = pbuffer.data(idx_op_geom_010_pf + 87);

    auto tr_0_0_z_z_yzz = pbuffer.data(idx_op_geom_010_pf + 88);

    auto tr_0_0_z_z_zzz = pbuffer.data(idx_op_geom_010_pf + 89);

    #pragma omp simd aligned(tr_0_0_z_z_xxx, tr_0_0_z_z_xxy, tr_0_0_z_z_xxz, tr_0_0_z_z_xyy, tr_0_0_z_z_xyz, tr_0_0_z_z_xzz, tr_0_0_z_z_yyy, tr_0_0_z_z_yyz, tr_0_0_z_z_yzz, tr_0_0_z_z_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_z_xx, tr_z_xxxz, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_z_xxx[i] = 2.0 * tr_zz_xxx[i] * tbe_0 + 2.0 * tr_z_xxxz[i] * tke_0 - tr_0_xxx[i];

        tr_0_0_z_z_xxy[i] = 2.0 * tr_zz_xxy[i] * tbe_0 + 2.0 * tr_z_xxyz[i] * tke_0 - tr_0_xxy[i];

        tr_0_0_z_z_xxz[i] = 2.0 * tr_zz_xxz[i] * tbe_0 + 2.0 * tr_z_xxzz[i] * tke_0 - tr_0_xxz[i] - tr_z_xx[i];

        tr_0_0_z_z_xyy[i] = 2.0 * tr_zz_xyy[i] * tbe_0 + 2.0 * tr_z_xyyz[i] * tke_0 - tr_0_xyy[i];

        tr_0_0_z_z_xyz[i] = 2.0 * tr_zz_xyz[i] * tbe_0 + 2.0 * tr_z_xyzz[i] * tke_0 - tr_0_xyz[i] - tr_z_xy[i];

        tr_0_0_z_z_xzz[i] = 2.0 * tr_zz_xzz[i] * tbe_0 + 2.0 * tr_z_xzzz[i] * tke_0 - tr_0_xzz[i] - 2.0 * tr_z_xz[i];

        tr_0_0_z_z_yyy[i] = 2.0 * tr_zz_yyy[i] * tbe_0 + 2.0 * tr_z_yyyz[i] * tke_0 - tr_0_yyy[i];

        tr_0_0_z_z_yyz[i] = 2.0 * tr_zz_yyz[i] * tbe_0 + 2.0 * tr_z_yyzz[i] * tke_0 - tr_0_yyz[i] - tr_z_yy[i];

        tr_0_0_z_z_yzz[i] = 2.0 * tr_zz_yzz[i] * tbe_0 + 2.0 * tr_z_yzzz[i] * tke_0 - tr_0_yzz[i] - 2.0 * tr_z_yz[i];

        tr_0_0_z_z_zzz[i] = 2.0 * tr_zz_zzz[i] * tbe_0 + 2.0 * tr_z_zzzz[i] * tke_0 - tr_0_zzz[i] - 3.0 * tr_z_zz[i];
    }

}

} // t2cgeom namespace

