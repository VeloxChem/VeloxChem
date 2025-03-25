#include "ThreeCenterOverlapGradientPrimRecDF.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_df(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_df,
                              const size_t idx_pf,
                              const size_t idx_dd,
                              const size_t idx_df,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : PF

    auto ts_x_xxx = pbuffer.data(idx_pf);

    auto ts_x_xxy = pbuffer.data(idx_pf + 1);

    auto ts_x_xxz = pbuffer.data(idx_pf + 2);

    auto ts_x_xyy = pbuffer.data(idx_pf + 3);

    auto ts_x_xyz = pbuffer.data(idx_pf + 4);

    auto ts_x_xzz = pbuffer.data(idx_pf + 5);

    auto ts_x_yyy = pbuffer.data(idx_pf + 6);

    auto ts_x_yyz = pbuffer.data(idx_pf + 7);

    auto ts_x_yzz = pbuffer.data(idx_pf + 8);

    auto ts_x_zzz = pbuffer.data(idx_pf + 9);

    auto ts_y_xxx = pbuffer.data(idx_pf + 10);

    auto ts_y_xxy = pbuffer.data(idx_pf + 11);

    auto ts_y_xxz = pbuffer.data(idx_pf + 12);

    auto ts_y_xyy = pbuffer.data(idx_pf + 13);

    auto ts_y_xyz = pbuffer.data(idx_pf + 14);

    auto ts_y_xzz = pbuffer.data(idx_pf + 15);

    auto ts_y_yyy = pbuffer.data(idx_pf + 16);

    auto ts_y_yyz = pbuffer.data(idx_pf + 17);

    auto ts_y_yzz = pbuffer.data(idx_pf + 18);

    auto ts_y_zzz = pbuffer.data(idx_pf + 19);

    auto ts_z_xxx = pbuffer.data(idx_pf + 20);

    auto ts_z_xxy = pbuffer.data(idx_pf + 21);

    auto ts_z_xxz = pbuffer.data(idx_pf + 22);

    auto ts_z_xyy = pbuffer.data(idx_pf + 23);

    auto ts_z_xyz = pbuffer.data(idx_pf + 24);

    auto ts_z_xzz = pbuffer.data(idx_pf + 25);

    auto ts_z_yyy = pbuffer.data(idx_pf + 26);

    auto ts_z_yyz = pbuffer.data(idx_pf + 27);

    auto ts_z_yzz = pbuffer.data(idx_pf + 28);

    auto ts_z_zzz = pbuffer.data(idx_pf + 29);

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_dd);

    auto ts_xx_xy = pbuffer.data(idx_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_dd + 5);

    auto ts_xy_xx = pbuffer.data(idx_dd + 6);

    auto ts_xy_xy = pbuffer.data(idx_dd + 7);

    auto ts_xy_xz = pbuffer.data(idx_dd + 8);

    auto ts_xy_yy = pbuffer.data(idx_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_dd + 10);

    auto ts_xy_zz = pbuffer.data(idx_dd + 11);

    auto ts_xz_xx = pbuffer.data(idx_dd + 12);

    auto ts_xz_xy = pbuffer.data(idx_dd + 13);

    auto ts_xz_xz = pbuffer.data(idx_dd + 14);

    auto ts_xz_yy = pbuffer.data(idx_dd + 15);

    auto ts_xz_yz = pbuffer.data(idx_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_dd + 23);

    auto ts_yz_xx = pbuffer.data(idx_dd + 24);

    auto ts_yz_xy = pbuffer.data(idx_dd + 25);

    auto ts_yz_xz = pbuffer.data(idx_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_dd + 35);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_df);

    auto ts_xx_xxy = pbuffer.data(idx_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_df + 9);

    auto ts_xy_xxx = pbuffer.data(idx_df + 10);

    auto ts_xy_xxy = pbuffer.data(idx_df + 11);

    auto ts_xy_xxz = pbuffer.data(idx_df + 12);

    auto ts_xy_xyy = pbuffer.data(idx_df + 13);

    auto ts_xy_xyz = pbuffer.data(idx_df + 14);

    auto ts_xy_xzz = pbuffer.data(idx_df + 15);

    auto ts_xy_yyy = pbuffer.data(idx_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_df + 18);

    auto ts_xy_zzz = pbuffer.data(idx_df + 19);

    auto ts_xz_xxx = pbuffer.data(idx_df + 20);

    auto ts_xz_xxy = pbuffer.data(idx_df + 21);

    auto ts_xz_xxz = pbuffer.data(idx_df + 22);

    auto ts_xz_xyy = pbuffer.data(idx_df + 23);

    auto ts_xz_xyz = pbuffer.data(idx_df + 24);

    auto ts_xz_xzz = pbuffer.data(idx_df + 25);

    auto ts_xz_yyy = pbuffer.data(idx_df + 26);

    auto ts_xz_yyz = pbuffer.data(idx_df + 27);

    auto ts_xz_yzz = pbuffer.data(idx_df + 28);

    auto ts_xz_zzz = pbuffer.data(idx_df + 29);

    auto ts_yy_xxx = pbuffer.data(idx_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_df + 39);

    auto ts_yz_xxx = pbuffer.data(idx_df + 40);

    auto ts_yz_xxy = pbuffer.data(idx_df + 41);

    auto ts_yz_xxz = pbuffer.data(idx_df + 42);

    auto ts_yz_xyy = pbuffer.data(idx_df + 43);

    auto ts_yz_xyz = pbuffer.data(idx_df + 44);

    auto ts_yz_xzz = pbuffer.data(idx_df + 45);

    auto ts_yz_yyy = pbuffer.data(idx_df + 46);

    auto ts_yz_yyz = pbuffer.data(idx_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_df + 48);

    auto ts_yz_zzz = pbuffer.data(idx_df + 49);

    auto ts_zz_xxx = pbuffer.data(idx_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_df + 59);

    // Set up 0-10 components of targeted buffer : DF

    auto gs_x_xx_xxx = pbuffer.data(idx_g_df);

    auto gs_x_xx_xxy = pbuffer.data(idx_g_df + 1);

    auto gs_x_xx_xxz = pbuffer.data(idx_g_df + 2);

    auto gs_x_xx_xyy = pbuffer.data(idx_g_df + 3);

    auto gs_x_xx_xyz = pbuffer.data(idx_g_df + 4);

    auto gs_x_xx_xzz = pbuffer.data(idx_g_df + 5);

    auto gs_x_xx_yyy = pbuffer.data(idx_g_df + 6);

    auto gs_x_xx_yyz = pbuffer.data(idx_g_df + 7);

    auto gs_x_xx_yzz = pbuffer.data(idx_g_df + 8);

    auto gs_x_xx_zzz = pbuffer.data(idx_g_df + 9);

    #pragma omp simd aligned(gc_x, gs_x_xx_xxx, gs_x_xx_xxy, gs_x_xx_xxz, gs_x_xx_xyy, gs_x_xx_xyz, gs_x_xx_xzz, gs_x_xx_yyy, gs_x_xx_yyz, gs_x_xx_yzz, gs_x_xx_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xx_xxx[i] = 4.0 * ts_x_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxx[i] * gc_x[i] * tce_0;

        gs_x_xx_xxy[i] = 4.0 * ts_x_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxy[i] * gc_x[i] * tce_0;

        gs_x_xx_xxz[i] = 4.0 * ts_x_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxz[i] * gc_x[i] * tce_0;

        gs_x_xx_xyy[i] = 4.0 * ts_x_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyy[i] * gc_x[i] * tce_0;

        gs_x_xx_xyz[i] = 4.0 * ts_x_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyz[i] * gc_x[i] * tce_0;

        gs_x_xx_xzz[i] = 4.0 * ts_x_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xzz[i] * gc_x[i] * tce_0;

        gs_x_xx_yyy[i] = 4.0 * ts_x_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyy[i] * gc_x[i] * tce_0;

        gs_x_xx_yyz[i] = 4.0 * ts_x_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyz[i] * gc_x[i] * tce_0;

        gs_x_xx_yzz[i] = 4.0 * ts_x_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzz[i] * gc_x[i] * tce_0;

        gs_x_xx_zzz[i] = 4.0 * ts_x_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 10-20 components of targeted buffer : DF

    auto gs_x_xy_xxx = pbuffer.data(idx_g_df + 10);

    auto gs_x_xy_xxy = pbuffer.data(idx_g_df + 11);

    auto gs_x_xy_xxz = pbuffer.data(idx_g_df + 12);

    auto gs_x_xy_xyy = pbuffer.data(idx_g_df + 13);

    auto gs_x_xy_xyz = pbuffer.data(idx_g_df + 14);

    auto gs_x_xy_xzz = pbuffer.data(idx_g_df + 15);

    auto gs_x_xy_yyy = pbuffer.data(idx_g_df + 16);

    auto gs_x_xy_yyz = pbuffer.data(idx_g_df + 17);

    auto gs_x_xy_yzz = pbuffer.data(idx_g_df + 18);

    auto gs_x_xy_zzz = pbuffer.data(idx_g_df + 19);

    #pragma omp simd aligned(gc_x, gs_x_xy_xxx, gs_x_xy_xxy, gs_x_xy_xxz, gs_x_xy_xyy, gs_x_xy_xyz, gs_x_xy_xzz, gs_x_xy_yyy, gs_x_xy_yyz, gs_x_xy_yzz, gs_x_xy_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xy_xxx[i] = 2.0 * ts_y_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xy_xxy[i] = 2.0 * ts_y_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xy_xxz[i] = 2.0 * ts_y_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xy_xyy[i] = 2.0 * ts_y_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xy_xyz[i] = 2.0 * ts_y_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xy_xzz[i] = 2.0 * ts_y_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xy_yyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xy_yyz[i] = 2.0 * ts_y_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xy_yzz[i] = 2.0 * ts_y_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xy_zzz[i] = 2.0 * ts_y_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 20-30 components of targeted buffer : DF

    auto gs_x_xz_xxx = pbuffer.data(idx_g_df + 20);

    auto gs_x_xz_xxy = pbuffer.data(idx_g_df + 21);

    auto gs_x_xz_xxz = pbuffer.data(idx_g_df + 22);

    auto gs_x_xz_xyy = pbuffer.data(idx_g_df + 23);

    auto gs_x_xz_xyz = pbuffer.data(idx_g_df + 24);

    auto gs_x_xz_xzz = pbuffer.data(idx_g_df + 25);

    auto gs_x_xz_yyy = pbuffer.data(idx_g_df + 26);

    auto gs_x_xz_yyz = pbuffer.data(idx_g_df + 27);

    auto gs_x_xz_yzz = pbuffer.data(idx_g_df + 28);

    auto gs_x_xz_zzz = pbuffer.data(idx_g_df + 29);

    #pragma omp simd aligned(gc_x, gs_x_xz_xxx, gs_x_xz_xxy, gs_x_xz_xxz, gs_x_xz_xyy, gs_x_xz_xyz, gs_x_xz_xzz, gs_x_xz_yyy, gs_x_xz_yyz, gs_x_xz_yzz, gs_x_xz_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xz_xxx[i] = 2.0 * ts_z_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xz_xxy[i] = 2.0 * ts_z_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xz_xxz[i] = 2.0 * ts_z_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xz_xyy[i] = 2.0 * ts_z_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xz_xyz[i] = 2.0 * ts_z_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xz_xzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xz_yyy[i] = 2.0 * ts_z_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xz_yyz[i] = 2.0 * ts_z_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xz_yzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-40 components of targeted buffer : DF

    auto gs_x_yy_xxx = pbuffer.data(idx_g_df + 30);

    auto gs_x_yy_xxy = pbuffer.data(idx_g_df + 31);

    auto gs_x_yy_xxz = pbuffer.data(idx_g_df + 32);

    auto gs_x_yy_xyy = pbuffer.data(idx_g_df + 33);

    auto gs_x_yy_xyz = pbuffer.data(idx_g_df + 34);

    auto gs_x_yy_xzz = pbuffer.data(idx_g_df + 35);

    auto gs_x_yy_yyy = pbuffer.data(idx_g_df + 36);

    auto gs_x_yy_yyz = pbuffer.data(idx_g_df + 37);

    auto gs_x_yy_yzz = pbuffer.data(idx_g_df + 38);

    auto gs_x_yy_zzz = pbuffer.data(idx_g_df + 39);

    #pragma omp simd aligned(gc_x, gs_x_yy_xxx, gs_x_yy_xxy, gs_x_yy_xxz, gs_x_yy_xyy, gs_x_yy_xyz, gs_x_yy_xzz, gs_x_yy_yyy, gs_x_yy_yyz, gs_x_yy_yzz, gs_x_yy_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yy_xxx[i] = 6.0 * ts_yy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxx[i] * gc_x[i] * tce_0;

        gs_x_yy_xxy[i] = 4.0 * ts_yy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxy[i] * gc_x[i] * tce_0;

        gs_x_yy_xxz[i] = 4.0 * ts_yy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxz[i] * gc_x[i] * tce_0;

        gs_x_yy_xyy[i] = 2.0 * ts_yy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyy[i] * gc_x[i] * tce_0;

        gs_x_yy_xyz[i] = 2.0 * ts_yy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyz[i] * gc_x[i] * tce_0;

        gs_x_yy_xzz[i] = 2.0 * ts_yy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzz[i] * gc_x[i] * tce_0;

        gs_x_yy_yyy[i] = 2.0 * ts_yy_yyy[i] * gc_x[i] * tce_0;

        gs_x_yy_yyz[i] = 2.0 * ts_yy_yyz[i] * gc_x[i] * tce_0;

        gs_x_yy_yzz[i] = 2.0 * ts_yy_yzz[i] * gc_x[i] * tce_0;

        gs_x_yy_zzz[i] = 2.0 * ts_yy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 40-50 components of targeted buffer : DF

    auto gs_x_yz_xxx = pbuffer.data(idx_g_df + 40);

    auto gs_x_yz_xxy = pbuffer.data(idx_g_df + 41);

    auto gs_x_yz_xxz = pbuffer.data(idx_g_df + 42);

    auto gs_x_yz_xyy = pbuffer.data(idx_g_df + 43);

    auto gs_x_yz_xyz = pbuffer.data(idx_g_df + 44);

    auto gs_x_yz_xzz = pbuffer.data(idx_g_df + 45);

    auto gs_x_yz_yyy = pbuffer.data(idx_g_df + 46);

    auto gs_x_yz_yyz = pbuffer.data(idx_g_df + 47);

    auto gs_x_yz_yzz = pbuffer.data(idx_g_df + 48);

    auto gs_x_yz_zzz = pbuffer.data(idx_g_df + 49);

    #pragma omp simd aligned(gc_x, gs_x_yz_xxx, gs_x_yz_xxy, gs_x_yz_xxz, gs_x_yz_xyy, gs_x_yz_xyz, gs_x_yz_xzz, gs_x_yz_yyy, gs_x_yz_yyz, gs_x_yz_yzz, gs_x_yz_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yz_xxx[i] = 6.0 * ts_yz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yz_xxy[i] = 4.0 * ts_yz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yz_xxz[i] = 4.0 * ts_yz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yz_xyy[i] = 2.0 * ts_yz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yz_xyz[i] = 2.0 * ts_yz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yz_xzz[i] = 2.0 * ts_yz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yz_yyy[i] = 2.0 * ts_yz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yz_yyz[i] = 2.0 * ts_yz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yz_yzz[i] = 2.0 * ts_yz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yz_zzz[i] = 2.0 * ts_yz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 50-60 components of targeted buffer : DF

    auto gs_x_zz_xxx = pbuffer.data(idx_g_df + 50);

    auto gs_x_zz_xxy = pbuffer.data(idx_g_df + 51);

    auto gs_x_zz_xxz = pbuffer.data(idx_g_df + 52);

    auto gs_x_zz_xyy = pbuffer.data(idx_g_df + 53);

    auto gs_x_zz_xyz = pbuffer.data(idx_g_df + 54);

    auto gs_x_zz_xzz = pbuffer.data(idx_g_df + 55);

    auto gs_x_zz_yyy = pbuffer.data(idx_g_df + 56);

    auto gs_x_zz_yyz = pbuffer.data(idx_g_df + 57);

    auto gs_x_zz_yzz = pbuffer.data(idx_g_df + 58);

    auto gs_x_zz_zzz = pbuffer.data(idx_g_df + 59);

    #pragma omp simd aligned(gc_x, gs_x_zz_xxx, gs_x_zz_xxy, gs_x_zz_xxz, gs_x_zz_xyy, gs_x_zz_xyz, gs_x_zz_xzz, gs_x_zz_yyy, gs_x_zz_yyz, gs_x_zz_yzz, gs_x_zz_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zz_xxx[i] = 6.0 * ts_zz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxx[i] * gc_x[i] * tce_0;

        gs_x_zz_xxy[i] = 4.0 * ts_zz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxy[i] * gc_x[i] * tce_0;

        gs_x_zz_xxz[i] = 4.0 * ts_zz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxz[i] * gc_x[i] * tce_0;

        gs_x_zz_xyy[i] = 2.0 * ts_zz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyy[i] * gc_x[i] * tce_0;

        gs_x_zz_xyz[i] = 2.0 * ts_zz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyz[i] * gc_x[i] * tce_0;

        gs_x_zz_xzz[i] = 2.0 * ts_zz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xzz[i] * gc_x[i] * tce_0;

        gs_x_zz_yyy[i] = 2.0 * ts_zz_yyy[i] * gc_x[i] * tce_0;

        gs_x_zz_yyz[i] = 2.0 * ts_zz_yyz[i] * gc_x[i] * tce_0;

        gs_x_zz_yzz[i] = 2.0 * ts_zz_yzz[i] * gc_x[i] * tce_0;

        gs_x_zz_zzz[i] = 2.0 * ts_zz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-70 components of targeted buffer : DF

    auto gs_y_xx_xxx = pbuffer.data(idx_g_df + 60);

    auto gs_y_xx_xxy = pbuffer.data(idx_g_df + 61);

    auto gs_y_xx_xxz = pbuffer.data(idx_g_df + 62);

    auto gs_y_xx_xyy = pbuffer.data(idx_g_df + 63);

    auto gs_y_xx_xyz = pbuffer.data(idx_g_df + 64);

    auto gs_y_xx_xzz = pbuffer.data(idx_g_df + 65);

    auto gs_y_xx_yyy = pbuffer.data(idx_g_df + 66);

    auto gs_y_xx_yyz = pbuffer.data(idx_g_df + 67);

    auto gs_y_xx_yzz = pbuffer.data(idx_g_df + 68);

    auto gs_y_xx_zzz = pbuffer.data(idx_g_df + 69);

    #pragma omp simd aligned(gc_y, gs_y_xx_xxx, gs_y_xx_xxy, gs_y_xx_xxz, gs_y_xx_xyy, gs_y_xx_xyz, gs_y_xx_xzz, gs_y_xx_yyy, gs_y_xx_yyz, gs_y_xx_yzz, gs_y_xx_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xx_xxx[i] = 2.0 * ts_xx_xxx[i] * gc_y[i] * tce_0;

        gs_y_xx_xxy[i] = 2.0 * ts_xx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxy[i] * gc_y[i] * tce_0;

        gs_y_xx_xxz[i] = 2.0 * ts_xx_xxz[i] * gc_y[i] * tce_0;

        gs_y_xx_xyy[i] = 4.0 * ts_xx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyy[i] * gc_y[i] * tce_0;

        gs_y_xx_xyz[i] = 2.0 * ts_xx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyz[i] * gc_y[i] * tce_0;

        gs_y_xx_xzz[i] = 2.0 * ts_xx_xzz[i] * gc_y[i] * tce_0;

        gs_y_xx_yyy[i] = 6.0 * ts_xx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyy[i] * gc_y[i] * tce_0;

        gs_y_xx_yyz[i] = 4.0 * ts_xx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyz[i] * gc_y[i] * tce_0;

        gs_y_xx_yzz[i] = 2.0 * ts_xx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzz[i] * gc_y[i] * tce_0;

        gs_y_xx_zzz[i] = 2.0 * ts_xx_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 70-80 components of targeted buffer : DF

    auto gs_y_xy_xxx = pbuffer.data(idx_g_df + 70);

    auto gs_y_xy_xxy = pbuffer.data(idx_g_df + 71);

    auto gs_y_xy_xxz = pbuffer.data(idx_g_df + 72);

    auto gs_y_xy_xyy = pbuffer.data(idx_g_df + 73);

    auto gs_y_xy_xyz = pbuffer.data(idx_g_df + 74);

    auto gs_y_xy_xzz = pbuffer.data(idx_g_df + 75);

    auto gs_y_xy_yyy = pbuffer.data(idx_g_df + 76);

    auto gs_y_xy_yyz = pbuffer.data(idx_g_df + 77);

    auto gs_y_xy_yzz = pbuffer.data(idx_g_df + 78);

    auto gs_y_xy_zzz = pbuffer.data(idx_g_df + 79);

    #pragma omp simd aligned(gc_y, gs_y_xy_xxx, gs_y_xy_xxy, gs_y_xy_xxz, gs_y_xy_xyy, gs_y_xy_xyz, gs_y_xy_xzz, gs_y_xy_yyy, gs_y_xy_yyz, gs_y_xy_yzz, gs_y_xy_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xy_xxx[i] = 2.0 * ts_x_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xy_xxy[i] = 2.0 * ts_x_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xy_xxz[i] = 2.0 * ts_x_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xy_xyy[i] = 2.0 * ts_x_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xy_xyz[i] = 2.0 * ts_x_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xy_xzz[i] = 2.0 * ts_x_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xy_yyy[i] = 2.0 * ts_x_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xy_yyz[i] = 2.0 * ts_x_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xy_yzz[i] = 2.0 * ts_x_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xy_zzz[i] = 2.0 * ts_x_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 80-90 components of targeted buffer : DF

    auto gs_y_xz_xxx = pbuffer.data(idx_g_df + 80);

    auto gs_y_xz_xxy = pbuffer.data(idx_g_df + 81);

    auto gs_y_xz_xxz = pbuffer.data(idx_g_df + 82);

    auto gs_y_xz_xyy = pbuffer.data(idx_g_df + 83);

    auto gs_y_xz_xyz = pbuffer.data(idx_g_df + 84);

    auto gs_y_xz_xzz = pbuffer.data(idx_g_df + 85);

    auto gs_y_xz_yyy = pbuffer.data(idx_g_df + 86);

    auto gs_y_xz_yyz = pbuffer.data(idx_g_df + 87);

    auto gs_y_xz_yzz = pbuffer.data(idx_g_df + 88);

    auto gs_y_xz_zzz = pbuffer.data(idx_g_df + 89);

    #pragma omp simd aligned(gc_y, gs_y_xz_xxx, gs_y_xz_xxy, gs_y_xz_xxz, gs_y_xz_xyy, gs_y_xz_xyz, gs_y_xz_xzz, gs_y_xz_yyy, gs_y_xz_yyz, gs_y_xz_yzz, gs_y_xz_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xz_xxx[i] = 2.0 * ts_xz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xz_xxy[i] = 2.0 * ts_xz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xz_xxz[i] = 2.0 * ts_xz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xz_xyy[i] = 4.0 * ts_xz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xz_xyz[i] = 2.0 * ts_xz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xz_xzz[i] = 2.0 * ts_xz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xz_yyy[i] = 6.0 * ts_xz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xz_yyz[i] = 4.0 * ts_xz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xz_yzz[i] = 2.0 * ts_xz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xz_zzz[i] = 2.0 * ts_xz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 90-100 components of targeted buffer : DF

    auto gs_y_yy_xxx = pbuffer.data(idx_g_df + 90);

    auto gs_y_yy_xxy = pbuffer.data(idx_g_df + 91);

    auto gs_y_yy_xxz = pbuffer.data(idx_g_df + 92);

    auto gs_y_yy_xyy = pbuffer.data(idx_g_df + 93);

    auto gs_y_yy_xyz = pbuffer.data(idx_g_df + 94);

    auto gs_y_yy_xzz = pbuffer.data(idx_g_df + 95);

    auto gs_y_yy_yyy = pbuffer.data(idx_g_df + 96);

    auto gs_y_yy_yyz = pbuffer.data(idx_g_df + 97);

    auto gs_y_yy_yzz = pbuffer.data(idx_g_df + 98);

    auto gs_y_yy_zzz = pbuffer.data(idx_g_df + 99);

    #pragma omp simd aligned(gc_y, gs_y_yy_xxx, gs_y_yy_xxy, gs_y_yy_xxz, gs_y_yy_xyy, gs_y_yy_xyz, gs_y_yy_xzz, gs_y_yy_yyy, gs_y_yy_yyz, gs_y_yy_yzz, gs_y_yy_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yy_xxx[i] = 4.0 * ts_y_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxx[i] * gc_y[i] * tce_0;

        gs_y_yy_xxy[i] = 4.0 * ts_y_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxy[i] * gc_y[i] * tce_0;

        gs_y_yy_xxz[i] = 4.0 * ts_y_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxz[i] * gc_y[i] * tce_0;

        gs_y_yy_xyy[i] = 4.0 * ts_y_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyy[i] * gc_y[i] * tce_0;

        gs_y_yy_xyz[i] = 4.0 * ts_y_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyz[i] * gc_y[i] * tce_0;

        gs_y_yy_xzz[i] = 4.0 * ts_y_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzz[i] * gc_y[i] * tce_0;

        gs_y_yy_yyy[i] = 4.0 * ts_y_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyy[i] * gc_y[i] * tce_0;

        gs_y_yy_yyz[i] = 4.0 * ts_y_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyz[i] * gc_y[i] * tce_0;

        gs_y_yy_yzz[i] = 4.0 * ts_y_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yzz[i] * gc_y[i] * tce_0;

        gs_y_yy_zzz[i] = 4.0 * ts_y_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 100-110 components of targeted buffer : DF

    auto gs_y_yz_xxx = pbuffer.data(idx_g_df + 100);

    auto gs_y_yz_xxy = pbuffer.data(idx_g_df + 101);

    auto gs_y_yz_xxz = pbuffer.data(idx_g_df + 102);

    auto gs_y_yz_xyy = pbuffer.data(idx_g_df + 103);

    auto gs_y_yz_xyz = pbuffer.data(idx_g_df + 104);

    auto gs_y_yz_xzz = pbuffer.data(idx_g_df + 105);

    auto gs_y_yz_yyy = pbuffer.data(idx_g_df + 106);

    auto gs_y_yz_yyz = pbuffer.data(idx_g_df + 107);

    auto gs_y_yz_yzz = pbuffer.data(idx_g_df + 108);

    auto gs_y_yz_zzz = pbuffer.data(idx_g_df + 109);

    #pragma omp simd aligned(gc_y, gs_y_yz_xxx, gs_y_yz_xxy, gs_y_yz_xxz, gs_y_yz_xyy, gs_y_yz_xyz, gs_y_yz_xzz, gs_y_yz_yyy, gs_y_yz_yyz, gs_y_yz_yzz, gs_y_yz_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yz_xxx[i] = 2.0 * ts_z_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yz_xxy[i] = 2.0 * ts_z_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yz_xxz[i] = 2.0 * ts_z_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yz_xyy[i] = 2.0 * ts_z_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yz_xyz[i] = 2.0 * ts_z_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yz_xzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yz_yyy[i] = 2.0 * ts_z_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yz_yyz[i] = 2.0 * ts_z_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yz_yzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 110-120 components of targeted buffer : DF

    auto gs_y_zz_xxx = pbuffer.data(idx_g_df + 110);

    auto gs_y_zz_xxy = pbuffer.data(idx_g_df + 111);

    auto gs_y_zz_xxz = pbuffer.data(idx_g_df + 112);

    auto gs_y_zz_xyy = pbuffer.data(idx_g_df + 113);

    auto gs_y_zz_xyz = pbuffer.data(idx_g_df + 114);

    auto gs_y_zz_xzz = pbuffer.data(idx_g_df + 115);

    auto gs_y_zz_yyy = pbuffer.data(idx_g_df + 116);

    auto gs_y_zz_yyz = pbuffer.data(idx_g_df + 117);

    auto gs_y_zz_yzz = pbuffer.data(idx_g_df + 118);

    auto gs_y_zz_zzz = pbuffer.data(idx_g_df + 119);

    #pragma omp simd aligned(gc_y, gs_y_zz_xxx, gs_y_zz_xxy, gs_y_zz_xxz, gs_y_zz_xyy, gs_y_zz_xyz, gs_y_zz_xzz, gs_y_zz_yyy, gs_y_zz_yyz, gs_y_zz_yzz, gs_y_zz_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zz_xxx[i] = 2.0 * ts_zz_xxx[i] * gc_y[i] * tce_0;

        gs_y_zz_xxy[i] = 2.0 * ts_zz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxy[i] * gc_y[i] * tce_0;

        gs_y_zz_xxz[i] = 2.0 * ts_zz_xxz[i] * gc_y[i] * tce_0;

        gs_y_zz_xyy[i] = 4.0 * ts_zz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyy[i] * gc_y[i] * tce_0;

        gs_y_zz_xyz[i] = 2.0 * ts_zz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyz[i] * gc_y[i] * tce_0;

        gs_y_zz_xzz[i] = 2.0 * ts_zz_xzz[i] * gc_y[i] * tce_0;

        gs_y_zz_yyy[i] = 6.0 * ts_zz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyy[i] * gc_y[i] * tce_0;

        gs_y_zz_yyz[i] = 4.0 * ts_zz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyz[i] * gc_y[i] * tce_0;

        gs_y_zz_yzz[i] = 2.0 * ts_zz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yzz[i] * gc_y[i] * tce_0;

        gs_y_zz_zzz[i] = 2.0 * ts_zz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 120-130 components of targeted buffer : DF

    auto gs_z_xx_xxx = pbuffer.data(idx_g_df + 120);

    auto gs_z_xx_xxy = pbuffer.data(idx_g_df + 121);

    auto gs_z_xx_xxz = pbuffer.data(idx_g_df + 122);

    auto gs_z_xx_xyy = pbuffer.data(idx_g_df + 123);

    auto gs_z_xx_xyz = pbuffer.data(idx_g_df + 124);

    auto gs_z_xx_xzz = pbuffer.data(idx_g_df + 125);

    auto gs_z_xx_yyy = pbuffer.data(idx_g_df + 126);

    auto gs_z_xx_yyz = pbuffer.data(idx_g_df + 127);

    auto gs_z_xx_yzz = pbuffer.data(idx_g_df + 128);

    auto gs_z_xx_zzz = pbuffer.data(idx_g_df + 129);

    #pragma omp simd aligned(gc_z, gs_z_xx_xxx, gs_z_xx_xxy, gs_z_xx_xxz, gs_z_xx_xyy, gs_z_xx_xyz, gs_z_xx_xzz, gs_z_xx_yyy, gs_z_xx_yyz, gs_z_xx_yzz, gs_z_xx_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xx_xxx[i] = 2.0 * ts_xx_xxx[i] * gc_z[i] * tce_0;

        gs_z_xx_xxy[i] = 2.0 * ts_xx_xxy[i] * gc_z[i] * tce_0;

        gs_z_xx_xxz[i] = 2.0 * ts_xx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxz[i] * gc_z[i] * tce_0;

        gs_z_xx_xyy[i] = 2.0 * ts_xx_xyy[i] * gc_z[i] * tce_0;

        gs_z_xx_xyz[i] = 2.0 * ts_xx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyz[i] * gc_z[i] * tce_0;

        gs_z_xx_xzz[i] = 4.0 * ts_xx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xzz[i] * gc_z[i] * tce_0;

        gs_z_xx_yyy[i] = 2.0 * ts_xx_yyy[i] * gc_z[i] * tce_0;

        gs_z_xx_yyz[i] = 2.0 * ts_xx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyz[i] * gc_z[i] * tce_0;

        gs_z_xx_yzz[i] = 4.0 * ts_xx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzz[i] * gc_z[i] * tce_0;

        gs_z_xx_zzz[i] = 6.0 * ts_xx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 130-140 components of targeted buffer : DF

    auto gs_z_xy_xxx = pbuffer.data(idx_g_df + 130);

    auto gs_z_xy_xxy = pbuffer.data(idx_g_df + 131);

    auto gs_z_xy_xxz = pbuffer.data(idx_g_df + 132);

    auto gs_z_xy_xyy = pbuffer.data(idx_g_df + 133);

    auto gs_z_xy_xyz = pbuffer.data(idx_g_df + 134);

    auto gs_z_xy_xzz = pbuffer.data(idx_g_df + 135);

    auto gs_z_xy_yyy = pbuffer.data(idx_g_df + 136);

    auto gs_z_xy_yyz = pbuffer.data(idx_g_df + 137);

    auto gs_z_xy_yzz = pbuffer.data(idx_g_df + 138);

    auto gs_z_xy_zzz = pbuffer.data(idx_g_df + 139);

    #pragma omp simd aligned(gc_z, gs_z_xy_xxx, gs_z_xy_xxy, gs_z_xy_xxz, gs_z_xy_xyy, gs_z_xy_xyz, gs_z_xy_xzz, gs_z_xy_yyy, gs_z_xy_yyz, gs_z_xy_yzz, gs_z_xy_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xy_xxx[i] = 2.0 * ts_xy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xy_xxy[i] = 2.0 * ts_xy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xy_xxz[i] = 2.0 * ts_xy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xy_xyy[i] = 2.0 * ts_xy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xy_xyz[i] = 2.0 * ts_xy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xy_xzz[i] = 4.0 * ts_xy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xy_yyy[i] = 2.0 * ts_xy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xy_yyz[i] = 2.0 * ts_xy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xy_yzz[i] = 4.0 * ts_xy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xy_zzz[i] = 6.0 * ts_xy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 140-150 components of targeted buffer : DF

    auto gs_z_xz_xxx = pbuffer.data(idx_g_df + 140);

    auto gs_z_xz_xxy = pbuffer.data(idx_g_df + 141);

    auto gs_z_xz_xxz = pbuffer.data(idx_g_df + 142);

    auto gs_z_xz_xyy = pbuffer.data(idx_g_df + 143);

    auto gs_z_xz_xyz = pbuffer.data(idx_g_df + 144);

    auto gs_z_xz_xzz = pbuffer.data(idx_g_df + 145);

    auto gs_z_xz_yyy = pbuffer.data(idx_g_df + 146);

    auto gs_z_xz_yyz = pbuffer.data(idx_g_df + 147);

    auto gs_z_xz_yzz = pbuffer.data(idx_g_df + 148);

    auto gs_z_xz_zzz = pbuffer.data(idx_g_df + 149);

    #pragma omp simd aligned(gc_z, gs_z_xz_xxx, gs_z_xz_xxy, gs_z_xz_xxz, gs_z_xz_xyy, gs_z_xz_xyz, gs_z_xz_xzz, gs_z_xz_yyy, gs_z_xz_yyz, gs_z_xz_yzz, gs_z_xz_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xz_xxx[i] = 2.0 * ts_x_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xz_xxy[i] = 2.0 * ts_x_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xz_xxz[i] = 2.0 * ts_x_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xz_xyy[i] = 2.0 * ts_x_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xz_xyz[i] = 2.0 * ts_x_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xz_xzz[i] = 2.0 * ts_x_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xz_yyy[i] = 2.0 * ts_x_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xz_yyz[i] = 2.0 * ts_x_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xz_yzz[i] = 2.0 * ts_x_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xz_zzz[i] = 2.0 * ts_x_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 150-160 components of targeted buffer : DF

    auto gs_z_yy_xxx = pbuffer.data(idx_g_df + 150);

    auto gs_z_yy_xxy = pbuffer.data(idx_g_df + 151);

    auto gs_z_yy_xxz = pbuffer.data(idx_g_df + 152);

    auto gs_z_yy_xyy = pbuffer.data(idx_g_df + 153);

    auto gs_z_yy_xyz = pbuffer.data(idx_g_df + 154);

    auto gs_z_yy_xzz = pbuffer.data(idx_g_df + 155);

    auto gs_z_yy_yyy = pbuffer.data(idx_g_df + 156);

    auto gs_z_yy_yyz = pbuffer.data(idx_g_df + 157);

    auto gs_z_yy_yzz = pbuffer.data(idx_g_df + 158);

    auto gs_z_yy_zzz = pbuffer.data(idx_g_df + 159);

    #pragma omp simd aligned(gc_z, gs_z_yy_xxx, gs_z_yy_xxy, gs_z_yy_xxz, gs_z_yy_xyy, gs_z_yy_xyz, gs_z_yy_xzz, gs_z_yy_yyy, gs_z_yy_yyz, gs_z_yy_yzz, gs_z_yy_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yy_xxx[i] = 2.0 * ts_yy_xxx[i] * gc_z[i] * tce_0;

        gs_z_yy_xxy[i] = 2.0 * ts_yy_xxy[i] * gc_z[i] * tce_0;

        gs_z_yy_xxz[i] = 2.0 * ts_yy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxz[i] * gc_z[i] * tce_0;

        gs_z_yy_xyy[i] = 2.0 * ts_yy_xyy[i] * gc_z[i] * tce_0;

        gs_z_yy_xyz[i] = 2.0 * ts_yy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyz[i] * gc_z[i] * tce_0;

        gs_z_yy_xzz[i] = 4.0 * ts_yy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzz[i] * gc_z[i] * tce_0;

        gs_z_yy_yyy[i] = 2.0 * ts_yy_yyy[i] * gc_z[i] * tce_0;

        gs_z_yy_yyz[i] = 2.0 * ts_yy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyz[i] * gc_z[i] * tce_0;

        gs_z_yy_yzz[i] = 4.0 * ts_yy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yzz[i] * gc_z[i] * tce_0;

        gs_z_yy_zzz[i] = 6.0 * ts_yy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 160-170 components of targeted buffer : DF

    auto gs_z_yz_xxx = pbuffer.data(idx_g_df + 160);

    auto gs_z_yz_xxy = pbuffer.data(idx_g_df + 161);

    auto gs_z_yz_xxz = pbuffer.data(idx_g_df + 162);

    auto gs_z_yz_xyy = pbuffer.data(idx_g_df + 163);

    auto gs_z_yz_xyz = pbuffer.data(idx_g_df + 164);

    auto gs_z_yz_xzz = pbuffer.data(idx_g_df + 165);

    auto gs_z_yz_yyy = pbuffer.data(idx_g_df + 166);

    auto gs_z_yz_yyz = pbuffer.data(idx_g_df + 167);

    auto gs_z_yz_yzz = pbuffer.data(idx_g_df + 168);

    auto gs_z_yz_zzz = pbuffer.data(idx_g_df + 169);

    #pragma omp simd aligned(gc_z, gs_z_yz_xxx, gs_z_yz_xxy, gs_z_yz_xxz, gs_z_yz_xyy, gs_z_yz_xyz, gs_z_yz_xzz, gs_z_yz_yyy, gs_z_yz_yyz, gs_z_yz_yzz, gs_z_yz_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yz_xxx[i] = 2.0 * ts_y_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yz_xxy[i] = 2.0 * ts_y_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yz_xxz[i] = 2.0 * ts_y_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yz_xyy[i] = 2.0 * ts_y_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yz_xyz[i] = 2.0 * ts_y_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yz_xzz[i] = 2.0 * ts_y_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yz_yyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yz_yyz[i] = 2.0 * ts_y_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yz_yzz[i] = 2.0 * ts_y_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yz_zzz[i] = 2.0 * ts_y_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 170-180 components of targeted buffer : DF

    auto gs_z_zz_xxx = pbuffer.data(idx_g_df + 170);

    auto gs_z_zz_xxy = pbuffer.data(idx_g_df + 171);

    auto gs_z_zz_xxz = pbuffer.data(idx_g_df + 172);

    auto gs_z_zz_xyy = pbuffer.data(idx_g_df + 173);

    auto gs_z_zz_xyz = pbuffer.data(idx_g_df + 174);

    auto gs_z_zz_xzz = pbuffer.data(idx_g_df + 175);

    auto gs_z_zz_yyy = pbuffer.data(idx_g_df + 176);

    auto gs_z_zz_yyz = pbuffer.data(idx_g_df + 177);

    auto gs_z_zz_yzz = pbuffer.data(idx_g_df + 178);

    auto gs_z_zz_zzz = pbuffer.data(idx_g_df + 179);

    #pragma omp simd aligned(gc_z, gs_z_zz_xxx, gs_z_zz_xxy, gs_z_zz_xxz, gs_z_zz_xyy, gs_z_zz_xyz, gs_z_zz_xzz, gs_z_zz_yyy, gs_z_zz_yyz, gs_z_zz_yzz, gs_z_zz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zz_xxx[i] = 4.0 * ts_z_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxx[i] * gc_z[i] * tce_0;

        gs_z_zz_xxy[i] = 4.0 * ts_z_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxy[i] * gc_z[i] * tce_0;

        gs_z_zz_xxz[i] = 4.0 * ts_z_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxz[i] * gc_z[i] * tce_0;

        gs_z_zz_xyy[i] = 4.0 * ts_z_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyy[i] * gc_z[i] * tce_0;

        gs_z_zz_xyz[i] = 4.0 * ts_z_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyz[i] * gc_z[i] * tce_0;

        gs_z_zz_xzz[i] = 4.0 * ts_z_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xzz[i] * gc_z[i] * tce_0;

        gs_z_zz_yyy[i] = 4.0 * ts_z_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyy[i] * gc_z[i] * tce_0;

        gs_z_zz_yyz[i] = 4.0 * ts_z_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyz[i] * gc_z[i] * tce_0;

        gs_z_zz_yzz[i] = 4.0 * ts_z_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yzz[i] * gc_z[i] * tce_0;

        gs_z_zz_zzz[i] = 4.0 * ts_z_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_zz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_zzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

