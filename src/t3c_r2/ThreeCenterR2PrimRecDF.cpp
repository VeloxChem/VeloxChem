#include "ThreeCenterR2PrimRecDF.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_df(CSimdArray<double>& pbuffer, 
                const size_t idx_g_df,
                const size_t idx_sf,
                const size_t idx_pd,
                const size_t idx_pf,
                const size_t idx_dp,
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

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_sf);

    auto ts_0_xxy = pbuffer.data(idx_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_pd);

    auto ts_x_xy = pbuffer.data(idx_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_pd + 17);

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

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    auto ts_xy_x = pbuffer.data(idx_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_dp + 5);

    auto ts_xz_x = pbuffer.data(idx_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    auto ts_yz_x = pbuffer.data(idx_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

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

    auto gr_xx_xxx = pbuffer.data(idx_g_df);

    auto gr_xx_xxy = pbuffer.data(idx_g_df + 1);

    auto gr_xx_xxz = pbuffer.data(idx_g_df + 2);

    auto gr_xx_xyy = pbuffer.data(idx_g_df + 3);

    auto gr_xx_xyz = pbuffer.data(idx_g_df + 4);

    auto gr_xx_xzz = pbuffer.data(idx_g_df + 5);

    auto gr_xx_yyy = pbuffer.data(idx_g_df + 6);

    auto gr_xx_yyz = pbuffer.data(idx_g_df + 7);

    auto gr_xx_yzz = pbuffer.data(idx_g_df + 8);

    auto gr_xx_zzz = pbuffer.data(idx_g_df + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xyy, gr_xx_xyz, gr_xx_xzz, gr_xx_yyy, gr_xx_yyz, gr_xx_yzz, gr_xx_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, ts_xx_x, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_y, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_z, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xx_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 + 12.0 * ts_x_xx[i] * gfe_0 + 4.0 * ts_x_xxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_x[i] * gfe_0 + 6.0 * ts_xx_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xx_xxx[i] * gfe_0 + ts_xx_xxx[i] * rgc2_0;

        gr_xx_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 + 8.0 * ts_x_xy[i] * gfe_0 + 4.0 * ts_x_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_y[i] * gfe_0 + 4.0 * ts_xx_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_xxy[i] * gfe_0 + ts_xx_xxy[i] * rgc2_0;

        gr_xx_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 + 8.0 * ts_x_xz[i] * gfe_0 + 4.0 * ts_x_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 + 4.0 * ts_xx_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xxz[i] * gfe_0 + ts_xx_xxz[i] * rgc2_0;

        gr_xx_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 + 4.0 * ts_x_yy[i] * gfe_0 + 4.0 * ts_x_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 + 4.0 * ts_xx_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_xyy[i] * gfe_0 + ts_xx_xyy[i] * rgc2_0;

        gr_xx_xyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 + 4.0 * ts_x_yz[i] * gfe_0 + 4.0 * ts_x_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xyz[i] * gfe_0 + ts_xx_xyz[i] * rgc2_0;

        gr_xx_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 + 4.0 * ts_x_zz[i] * gfe_0 + 4.0 * ts_x_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 + 4.0 * ts_xx_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xzz[i] * gfe_0 + ts_xx_xzz[i] * rgc2_0;

        gr_xx_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 + 4.0 * ts_x_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_y[i] * gfe_0 + 6.0 * ts_xx_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_yyy[i] * gfe_0 + ts_xx_yyy[i] * rgc2_0;

        gr_xx_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 + 4.0 * ts_x_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 + 4.0 * ts_xx_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_yyz[i] * gfe_0 + ts_xx_yyz[i] * rgc2_0;

        gr_xx_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 + 4.0 * ts_x_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_y[i] * gfe_0 + 4.0 * ts_xx_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_yzz[i] * gfe_0 + ts_xx_yzz[i] * rgc2_0;

        gr_xx_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 + 4.0 * ts_x_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_z[i] * gfe_0 + 6.0 * ts_xx_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_zzz[i] * gfe_0 + ts_xx_zzz[i] * rgc2_0;
    }

    // Set up 10-20 components of targeted buffer : DF

    auto gr_xy_xxx = pbuffer.data(idx_g_df + 10);

    auto gr_xy_xxy = pbuffer.data(idx_g_df + 11);

    auto gr_xy_xxz = pbuffer.data(idx_g_df + 12);

    auto gr_xy_xyy = pbuffer.data(idx_g_df + 13);

    auto gr_xy_xyz = pbuffer.data(idx_g_df + 14);

    auto gr_xy_xzz = pbuffer.data(idx_g_df + 15);

    auto gr_xy_yyy = pbuffer.data(idx_g_df + 16);

    auto gr_xy_yyz = pbuffer.data(idx_g_df + 17);

    auto gr_xy_yzz = pbuffer.data(idx_g_df + 18);

    auto gr_xy_zzz = pbuffer.data(idx_g_df + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xyy, gr_xy_xyz, gr_xy_xzz, gr_xy_yyy, gr_xy_yyz, gr_xy_yzz, gr_xy_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, ts_xy_x, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_y, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_z, ts_xy_zz, ts_xy_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xy_xxx[i] = 6.0 * ts_y_xx[i] * gfe_0 + 2.0 * ts_y_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_x[i] * gfe_0 + 6.0 * ts_xy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xy_xxx[i] * gfe_0 + ts_xy_xxx[i] * rgc2_0;

        gr_xy_xxy[i] = 4.0 * ts_y_xy[i] * gfe_0 + 2.0 * ts_y_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 + 2.0 * ts_x_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe_0 + 4.0 * ts_xy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_xxy[i] * gfe_0 + ts_xy_xxy[i] * rgc2_0;

        gr_xy_xxz[i] = 4.0 * ts_y_xz[i] * gfe_0 + 2.0 * ts_y_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_z[i] * gfe_0 + 4.0 * ts_xy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xxz[i] * gfe_0 + ts_xy_xxz[i] * rgc2_0;

        gr_xy_xyy[i] = 2.0 * ts_y_yy[i] * gfe_0 + 2.0 * ts_y_xyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_xy[i] * gfe_0 + 2.0 * ts_x_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_x[i] * gfe_0 + 4.0 * ts_xy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_xyy[i] * gfe_0 + ts_xy_xyy[i] * rgc2_0;

        gr_xy_xyz[i] = 2.0 * ts_y_yz[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xz[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xyz[i] * gfe_0 + ts_xy_xyz[i] * rgc2_0;

        gr_xy_xzz[i] = 2.0 * ts_y_zz[i] * gfe_0 + 2.0 * ts_y_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_x[i] * gfe_0 + 4.0 * ts_xy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xzz[i] * gfe_0 + ts_xy_xzz[i] * rgc2_0;

        gr_xy_yyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_yy[i] * gfe_0 + 2.0 * ts_x_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_y[i] * gfe_0 + 6.0 * ts_xy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_yyy[i] * gfe_0 + ts_xy_yyy[i] * rgc2_0;

        gr_xy_yyz[i] = 2.0 * ts_y_yyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_yz[i] * gfe_0 + 2.0 * ts_x_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_z[i] * gfe_0 + 4.0 * ts_xy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_yyz[i] * gfe_0 + ts_xy_yyz[i] * rgc2_0;

        gr_xy_yzz[i] = 2.0 * ts_y_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe_0 + 2.0 * ts_x_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe_0 + 4.0 * ts_xy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_yzz[i] * gfe_0 + ts_xy_yzz[i] * rgc2_0;

        gr_xy_zzz[i] = 2.0 * ts_y_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_z[i] * gfe_0 + 6.0 * ts_xy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_zzz[i] * gfe_0 + ts_xy_zzz[i] * rgc2_0;
    }

    // Set up 20-30 components of targeted buffer : DF

    auto gr_xz_xxx = pbuffer.data(idx_g_df + 20);

    auto gr_xz_xxy = pbuffer.data(idx_g_df + 21);

    auto gr_xz_xxz = pbuffer.data(idx_g_df + 22);

    auto gr_xz_xyy = pbuffer.data(idx_g_df + 23);

    auto gr_xz_xyz = pbuffer.data(idx_g_df + 24);

    auto gr_xz_xzz = pbuffer.data(idx_g_df + 25);

    auto gr_xz_yyy = pbuffer.data(idx_g_df + 26);

    auto gr_xz_yyz = pbuffer.data(idx_g_df + 27);

    auto gr_xz_yzz = pbuffer.data(idx_g_df + 28);

    auto gr_xz_zzz = pbuffer.data(idx_g_df + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xyy, gr_xz_xyz, gr_xz_xzz, gr_xz_yyy, gr_xz_yyz, gr_xz_yzz, gr_xz_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, ts_xz_x, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_y, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_z, ts_xz_zz, ts_xz_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xz_xxx[i] = 6.0 * ts_z_xx[i] * gfe_0 + 2.0 * ts_z_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xz_x[i] * gfe_0 + 6.0 * ts_xz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xz_xxx[i] * gfe_0 + ts_xz_xxx[i] * rgc2_0;

        gr_xz_xxy[i] = 4.0 * ts_z_xy[i] * gfe_0 + 2.0 * ts_z_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_y[i] * gfe_0 + 4.0 * ts_xz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_xxy[i] * gfe_0 + ts_xz_xxy[i] * rgc2_0;

        gr_xz_xxz[i] = 4.0 * ts_z_xz[i] * gfe_0 + 2.0 * ts_z_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 + 2.0 * ts_x_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_z[i] * gfe_0 + 4.0 * ts_xz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xxz[i] * gfe_0 + ts_xz_xxz[i] * rgc2_0;

        gr_xz_xyy[i] = 2.0 * ts_z_yy[i] * gfe_0 + 2.0 * ts_z_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe_0 + 4.0 * ts_xz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_xyy[i] * gfe_0 + ts_xz_xyy[i] * rgc2_0;

        gr_xz_xyz[i] = 2.0 * ts_z_yz[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xy[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xyz[i] * gfe_0 + ts_xz_xyz[i] * rgc2_0;

        gr_xz_xzz[i] = 2.0 * ts_z_zz[i] * gfe_0 + 2.0 * ts_z_xzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_xz[i] * gfe_0 + 2.0 * ts_x_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe_0 + 4.0 * ts_xz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xzz[i] * gfe_0 + ts_xz_xzz[i] * rgc2_0;

        gr_xz_yyy[i] = 2.0 * ts_z_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xz_y[i] * gfe_0 + 6.0 * ts_xz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_yyy[i] * gfe_0 + ts_xz_yyy[i] * rgc2_0;

        gr_xz_yyz[i] = 2.0 * ts_z_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yy[i] * gfe_0 + 2.0 * ts_x_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_z[i] * gfe_0 + 4.0 * ts_xz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_yyz[i] * gfe_0 + ts_xz_yyz[i] * rgc2_0;

        gr_xz_yzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_yz[i] * gfe_0 + 2.0 * ts_x_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_y[i] * gfe_0 + 4.0 * ts_xz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_yzz[i] * gfe_0 + ts_xz_yzz[i] * rgc2_0;

        gr_xz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_zz[i] * gfe_0 + 2.0 * ts_x_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xz_z[i] * gfe_0 + 6.0 * ts_xz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_zzz[i] * gfe_0 + ts_xz_zzz[i] * rgc2_0;
    }

    // Set up 30-40 components of targeted buffer : DF

    auto gr_yy_xxx = pbuffer.data(idx_g_df + 30);

    auto gr_yy_xxy = pbuffer.data(idx_g_df + 31);

    auto gr_yy_xxz = pbuffer.data(idx_g_df + 32);

    auto gr_yy_xyy = pbuffer.data(idx_g_df + 33);

    auto gr_yy_xyz = pbuffer.data(idx_g_df + 34);

    auto gr_yy_xzz = pbuffer.data(idx_g_df + 35);

    auto gr_yy_yyy = pbuffer.data(idx_g_df + 36);

    auto gr_yy_yyz = pbuffer.data(idx_g_df + 37);

    auto gr_yy_yzz = pbuffer.data(idx_g_df + 38);

    auto gr_yy_zzz = pbuffer.data(idx_g_df + 39);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xyy, gr_yy_xyz, gr_yy_xzz, gr_yy_yyy, gr_yy_yyz, gr_yy_yzz, gr_yy_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, ts_yy_x, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_y, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_z, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yy_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 + 4.0 * ts_y_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_x[i] * gfe_0 + 6.0 * ts_yy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yy_xxx[i] * gfe_0 + ts_yy_xxx[i] * rgc2_0;

        gr_yy_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 + 4.0 * ts_y_xx[i] * gfe_0 + 4.0 * ts_y_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe_0 + 4.0 * ts_yy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_xxy[i] * gfe_0 + ts_yy_xxy[i] * rgc2_0;

        gr_yy_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 + 4.0 * ts_y_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_z[i] * gfe_0 + 4.0 * ts_yy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xxz[i] * gfe_0 + ts_yy_xxz[i] * rgc2_0;

        gr_yy_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 + 8.0 * ts_y_xy[i] * gfe_0 + 4.0 * ts_y_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_x[i] * gfe_0 + 4.0 * ts_yy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_xyy[i] * gfe_0 + ts_yy_xyy[i] * rgc2_0;

        gr_yy_xyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 + 4.0 * ts_y_xz[i] * gfe_0 + 4.0 * ts_y_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xyz[i] * gfe_0 + ts_yy_xyz[i] * rgc2_0;

        gr_yy_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 + 4.0 * ts_y_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_x[i] * gfe_0 + 4.0 * ts_yy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xzz[i] * gfe_0 + ts_yy_xzz[i] * rgc2_0;

        gr_yy_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 + 12.0 * ts_y_yy[i] * gfe_0 + 4.0 * ts_y_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_y[i] * gfe_0 + 6.0 * ts_yy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_yyy[i] * gfe_0 + ts_yy_yyy[i] * rgc2_0;

        gr_yy_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 + 8.0 * ts_y_yz[i] * gfe_0 + 4.0 * ts_y_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_z[i] * gfe_0 + 4.0 * ts_yy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_yyz[i] * gfe_0 + ts_yy_yyz[i] * rgc2_0;

        gr_yy_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 + 4.0 * ts_y_zz[i] * gfe_0 + 4.0 * ts_y_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe_0 + 4.0 * ts_yy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_yzz[i] * gfe_0 + ts_yy_yzz[i] * rgc2_0;

        gr_yy_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 + 4.0 * ts_y_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_z[i] * gfe_0 + 6.0 * ts_yy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_zzz[i] * gfe_0 + ts_yy_zzz[i] * rgc2_0;
    }

    // Set up 40-50 components of targeted buffer : DF

    auto gr_yz_xxx = pbuffer.data(idx_g_df + 40);

    auto gr_yz_xxy = pbuffer.data(idx_g_df + 41);

    auto gr_yz_xxz = pbuffer.data(idx_g_df + 42);

    auto gr_yz_xyy = pbuffer.data(idx_g_df + 43);

    auto gr_yz_xyz = pbuffer.data(idx_g_df + 44);

    auto gr_yz_xzz = pbuffer.data(idx_g_df + 45);

    auto gr_yz_yyy = pbuffer.data(idx_g_df + 46);

    auto gr_yz_yyz = pbuffer.data(idx_g_df + 47);

    auto gr_yz_yzz = pbuffer.data(idx_g_df + 48);

    auto gr_yz_zzz = pbuffer.data(idx_g_df + 49);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xyy, gr_yz_xyz, gr_yz_xzz, gr_yz_yyy, gr_yz_yyz, gr_yz_yzz, gr_yz_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, ts_yz_x, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_y, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_z, ts_yz_zz, ts_yz_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yz_xxx[i] = 2.0 * ts_z_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_yz_x[i] * gfe_0 + 6.0 * ts_yz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yz_xxx[i] * gfe_0 + ts_yz_xxx[i] * rgc2_0;

        gr_yz_xxy[i] = 2.0 * ts_z_xx[i] * gfe_0 + 2.0 * ts_z_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_y[i] * gfe_0 + 4.0 * ts_yz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_xxy[i] * gfe_0 + ts_yz_xxy[i] * rgc2_0;

        gr_yz_xxz[i] = 2.0 * ts_z_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xx[i] * gfe_0 + 2.0 * ts_y_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_z[i] * gfe_0 + 4.0 * ts_yz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xxz[i] * gfe_0 + ts_yz_xxz[i] * rgc2_0;

        gr_yz_xyy[i] = 4.0 * ts_z_xy[i] * gfe_0 + 2.0 * ts_z_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_x[i] * gfe_0 + 4.0 * ts_yz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_xyy[i] * gfe_0 + ts_yz_xyy[i] * rgc2_0;

        gr_yz_xyz[i] = 2.0 * ts_z_xz[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xy[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xyz[i] * gfe_0 + ts_yz_xyz[i] * rgc2_0;

        gr_yz_xzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_y_xz[i] * gfe_0 + 2.0 * ts_y_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_x[i] * gfe_0 + 4.0 * ts_yz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xzz[i] * gfe_0 + ts_yz_xzz[i] * rgc2_0;

        gr_yz_yyy[i] = 6.0 * ts_z_yy[i] * gfe_0 + 2.0 * ts_z_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yz_y[i] * gfe_0 + 6.0 * ts_yz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_yyy[i] * gfe_0 + ts_yz_yyy[i] * rgc2_0;

        gr_yz_yyz[i] = 4.0 * ts_z_yz[i] * gfe_0 + 2.0 * ts_z_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe_0 + 2.0 * ts_y_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_z[i] * gfe_0 + 4.0 * ts_yz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_yyz[i] * gfe_0 + ts_yz_yyz[i] * rgc2_0;

        gr_yz_yzz[i] = 2.0 * ts_z_zz[i] * gfe_0 + 2.0 * ts_z_yzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_y_yz[i] * gfe_0 + 2.0 * ts_y_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_y[i] * gfe_0 + 4.0 * ts_yz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_yzz[i] * gfe_0 + ts_yz_yzz[i] * rgc2_0;

        gr_yz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_zz[i] * gfe_0 + 2.0 * ts_y_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yz_z[i] * gfe_0 + 6.0 * ts_yz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_zzz[i] * gfe_0 + ts_yz_zzz[i] * rgc2_0;
    }

    // Set up 50-60 components of targeted buffer : DF

    auto gr_zz_xxx = pbuffer.data(idx_g_df + 50);

    auto gr_zz_xxy = pbuffer.data(idx_g_df + 51);

    auto gr_zz_xxz = pbuffer.data(idx_g_df + 52);

    auto gr_zz_xyy = pbuffer.data(idx_g_df + 53);

    auto gr_zz_xyz = pbuffer.data(idx_g_df + 54);

    auto gr_zz_xzz = pbuffer.data(idx_g_df + 55);

    auto gr_zz_yyy = pbuffer.data(idx_g_df + 56);

    auto gr_zz_yyz = pbuffer.data(idx_g_df + 57);

    auto gr_zz_yzz = pbuffer.data(idx_g_df + 58);

    auto gr_zz_zzz = pbuffer.data(idx_g_df + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xyy, gr_zz_xyz, gr_zz_xzz, gr_zz_yyy, gr_zz_yyz, gr_zz_yzz, gr_zz_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, ts_zz_x, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_y, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_z, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zz_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 + 4.0 * ts_z_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_zz_x[i] * gfe_0 + 6.0 * ts_zz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_zz_xxx[i] * gfe_0 + ts_zz_xxx[i] * rgc2_0;

        gr_zz_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 + 4.0 * ts_z_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_y[i] * gfe_0 + 4.0 * ts_zz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_xxy[i] * gfe_0 + ts_zz_xxy[i] * rgc2_0;

        gr_zz_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 + 4.0 * ts_z_xx[i] * gfe_0 + 4.0 * ts_z_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_z[i] * gfe_0 + 4.0 * ts_zz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xxz[i] * gfe_0 + ts_zz_xxz[i] * rgc2_0;

        gr_zz_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 + 4.0 * ts_z_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_x[i] * gfe_0 + 4.0 * ts_zz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_xyy[i] * gfe_0 + ts_zz_xyy[i] * rgc2_0;

        gr_zz_xyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 + 4.0 * ts_z_xy[i] * gfe_0 + 4.0 * ts_z_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xyz[i] * gfe_0 + ts_zz_xyz[i] * rgc2_0;

        gr_zz_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 + 8.0 * ts_z_xz[i] * gfe_0 + 4.0 * ts_z_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_x[i] * gfe_0 + 4.0 * ts_zz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xzz[i] * gfe_0 + ts_zz_xzz[i] * rgc2_0;

        gr_zz_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 + 4.0 * ts_z_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_zz_y[i] * gfe_0 + 6.0 * ts_zz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_yyy[i] * gfe_0 + ts_zz_yyy[i] * rgc2_0;

        gr_zz_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 + 4.0 * ts_z_yy[i] * gfe_0 + 4.0 * ts_z_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_z[i] * gfe_0 + 4.0 * ts_zz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_yyz[i] * gfe_0 + ts_zz_yyz[i] * rgc2_0;

        gr_zz_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 + 8.0 * ts_z_yz[i] * gfe_0 + 4.0 * ts_z_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_y[i] * gfe_0 + 4.0 * ts_zz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_yzz[i] * gfe_0 + ts_zz_yzz[i] * rgc2_0;

        gr_zz_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 + 12.0 * ts_z_zz[i] * gfe_0 + 4.0 * ts_z_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zz_z[i] * gfe_0 + 6.0 * ts_zz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_zzz[i] * gfe_0 + ts_zz_zzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

