#include "ThreeCenterRR2PrimRecDF.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_df(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_df,
                  const size_t idx_pf,
                  const size_t idx_g_pf,
                  const size_t idx_dd,
                  const size_t idx_g_dd,
                  const size_t idx_df,
                  const size_t idx_g_df,
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

    // Set up components of auxiliary buffer : PF

    auto gr_x_xxx = pbuffer.data(idx_g_pf);

    auto gr_x_xxy = pbuffer.data(idx_g_pf + 1);

    auto gr_x_xxz = pbuffer.data(idx_g_pf + 2);

    auto gr_x_xyy = pbuffer.data(idx_g_pf + 3);

    auto gr_x_xyz = pbuffer.data(idx_g_pf + 4);

    auto gr_x_xzz = pbuffer.data(idx_g_pf + 5);

    auto gr_x_yyy = pbuffer.data(idx_g_pf + 6);

    auto gr_x_yyz = pbuffer.data(idx_g_pf + 7);

    auto gr_x_yzz = pbuffer.data(idx_g_pf + 8);

    auto gr_x_zzz = pbuffer.data(idx_g_pf + 9);

    auto gr_y_xxx = pbuffer.data(idx_g_pf + 10);

    auto gr_y_xxy = pbuffer.data(idx_g_pf + 11);

    auto gr_y_xxz = pbuffer.data(idx_g_pf + 12);

    auto gr_y_xyy = pbuffer.data(idx_g_pf + 13);

    auto gr_y_xyz = pbuffer.data(idx_g_pf + 14);

    auto gr_y_xzz = pbuffer.data(idx_g_pf + 15);

    auto gr_y_yyy = pbuffer.data(idx_g_pf + 16);

    auto gr_y_yyz = pbuffer.data(idx_g_pf + 17);

    auto gr_y_yzz = pbuffer.data(idx_g_pf + 18);

    auto gr_y_zzz = pbuffer.data(idx_g_pf + 19);

    auto gr_z_xxx = pbuffer.data(idx_g_pf + 20);

    auto gr_z_xxy = pbuffer.data(idx_g_pf + 21);

    auto gr_z_xxz = pbuffer.data(idx_g_pf + 22);

    auto gr_z_xyy = pbuffer.data(idx_g_pf + 23);

    auto gr_z_xyz = pbuffer.data(idx_g_pf + 24);

    auto gr_z_xzz = pbuffer.data(idx_g_pf + 25);

    auto gr_z_yyy = pbuffer.data(idx_g_pf + 26);

    auto gr_z_yyz = pbuffer.data(idx_g_pf + 27);

    auto gr_z_yzz = pbuffer.data(idx_g_pf + 28);

    auto gr_z_zzz = pbuffer.data(idx_g_pf + 29);

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

    // Set up components of auxiliary buffer : DD

    auto gr_xx_xx = pbuffer.data(idx_g_dd);

    auto gr_xx_xy = pbuffer.data(idx_g_dd + 1);

    auto gr_xx_xz = pbuffer.data(idx_g_dd + 2);

    auto gr_xx_yy = pbuffer.data(idx_g_dd + 3);

    auto gr_xx_yz = pbuffer.data(idx_g_dd + 4);

    auto gr_xx_zz = pbuffer.data(idx_g_dd + 5);

    auto gr_xy_xx = pbuffer.data(idx_g_dd + 6);

    auto gr_xy_xy = pbuffer.data(idx_g_dd + 7);

    auto gr_xy_xz = pbuffer.data(idx_g_dd + 8);

    auto gr_xy_yy = pbuffer.data(idx_g_dd + 9);

    auto gr_xy_yz = pbuffer.data(idx_g_dd + 10);

    auto gr_xy_zz = pbuffer.data(idx_g_dd + 11);

    auto gr_xz_xx = pbuffer.data(idx_g_dd + 12);

    auto gr_xz_xy = pbuffer.data(idx_g_dd + 13);

    auto gr_xz_xz = pbuffer.data(idx_g_dd + 14);

    auto gr_xz_yy = pbuffer.data(idx_g_dd + 15);

    auto gr_xz_yz = pbuffer.data(idx_g_dd + 16);

    auto gr_xz_zz = pbuffer.data(idx_g_dd + 17);

    auto gr_yy_xx = pbuffer.data(idx_g_dd + 18);

    auto gr_yy_xy = pbuffer.data(idx_g_dd + 19);

    auto gr_yy_xz = pbuffer.data(idx_g_dd + 20);

    auto gr_yy_yy = pbuffer.data(idx_g_dd + 21);

    auto gr_yy_yz = pbuffer.data(idx_g_dd + 22);

    auto gr_yy_zz = pbuffer.data(idx_g_dd + 23);

    auto gr_yz_xx = pbuffer.data(idx_g_dd + 24);

    auto gr_yz_xy = pbuffer.data(idx_g_dd + 25);

    auto gr_yz_xz = pbuffer.data(idx_g_dd + 26);

    auto gr_yz_yy = pbuffer.data(idx_g_dd + 27);

    auto gr_yz_yz = pbuffer.data(idx_g_dd + 28);

    auto gr_yz_zz = pbuffer.data(idx_g_dd + 29);

    auto gr_zz_xx = pbuffer.data(idx_g_dd + 30);

    auto gr_zz_xy = pbuffer.data(idx_g_dd + 31);

    auto gr_zz_xz = pbuffer.data(idx_g_dd + 32);

    auto gr_zz_yy = pbuffer.data(idx_g_dd + 33);

    auto gr_zz_yz = pbuffer.data(idx_g_dd + 34);

    auto gr_zz_zz = pbuffer.data(idx_g_dd + 35);

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

    // Set up components of auxiliary buffer : DF

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

    // Set up 0-10 components of targeted buffer : DF

    auto grr_x_xx_xxx = pbuffer.data(idx_gr_df);

    auto grr_x_xx_xxy = pbuffer.data(idx_gr_df + 1);

    auto grr_x_xx_xxz = pbuffer.data(idx_gr_df + 2);

    auto grr_x_xx_xyy = pbuffer.data(idx_gr_df + 3);

    auto grr_x_xx_xyz = pbuffer.data(idx_gr_df + 4);

    auto grr_x_xx_xzz = pbuffer.data(idx_gr_df + 5);

    auto grr_x_xx_yyy = pbuffer.data(idx_gr_df + 6);

    auto grr_x_xx_yyz = pbuffer.data(idx_gr_df + 7);

    auto grr_x_xx_yzz = pbuffer.data(idx_gr_df + 8);

    auto grr_x_xx_zzz = pbuffer.data(idx_gr_df + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xyy, gr_x_xyz, gr_x_xzz, gr_x_yyy, gr_x_yyz, gr_x_yzz, gr_x_zzz, gr_xx_xx, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xy, gr_xx_xyy, gr_xx_xyz, gr_xx_xz, gr_xx_xzz, gr_xx_yy, gr_xx_yyy, gr_xx_yyz, gr_xx_yz, gr_xx_yzz, gr_xx_zz, gr_xx_zzz, grr_x_xx_xxx, grr_x_xx_xxy, grr_x_xx_xxz, grr_x_xx_xyy, grr_x_xx_xyz, grr_x_xx_xzz, grr_x_xx_yyy, grr_x_xx_yyz, grr_x_xx_yzz, grr_x_xx_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xx_xxx[i] = 2.0 * ts_x_xxx[i] * gfe2_0 + 2.0 * gr_x_xxx[i] * gfe_0 + 3.0 * ts_xx_xx[i] * gfe2_0 + 3.0 * gr_xx_xx[i] * gfe_0 + ts_xx_xxx[i] * gfe_0 * gc_x[i] + gr_xx_xxx[i] * gc_x[i];

        grr_x_xx_xxy[i] = 2.0 * ts_x_xxy[i] * gfe2_0 + 2.0 * gr_x_xxy[i] * gfe_0 + 2.0 * ts_xx_xy[i] * gfe2_0 + 2.0 * gr_xx_xy[i] * gfe_0 + ts_xx_xxy[i] * gfe_0 * gc_x[i] + gr_xx_xxy[i] * gc_x[i];

        grr_x_xx_xxz[i] = 2.0 * ts_x_xxz[i] * gfe2_0 + 2.0 * gr_x_xxz[i] * gfe_0 + 2.0 * ts_xx_xz[i] * gfe2_0 + 2.0 * gr_xx_xz[i] * gfe_0 + ts_xx_xxz[i] * gfe_0 * gc_x[i] + gr_xx_xxz[i] * gc_x[i];

        grr_x_xx_xyy[i] = 2.0 * ts_x_xyy[i] * gfe2_0 + 2.0 * gr_x_xyy[i] * gfe_0 + ts_xx_yy[i] * gfe2_0 + gr_xx_yy[i] * gfe_0 + ts_xx_xyy[i] * gfe_0 * gc_x[i] + gr_xx_xyy[i] * gc_x[i];

        grr_x_xx_xyz[i] = 2.0 * ts_x_xyz[i] * gfe2_0 + 2.0 * gr_x_xyz[i] * gfe_0 + ts_xx_yz[i] * gfe2_0 + gr_xx_yz[i] * gfe_0 + ts_xx_xyz[i] * gfe_0 * gc_x[i] + gr_xx_xyz[i] * gc_x[i];

        grr_x_xx_xzz[i] = 2.0 * ts_x_xzz[i] * gfe2_0 + 2.0 * gr_x_xzz[i] * gfe_0 + ts_xx_zz[i] * gfe2_0 + gr_xx_zz[i] * gfe_0 + ts_xx_xzz[i] * gfe_0 * gc_x[i] + gr_xx_xzz[i] * gc_x[i];

        grr_x_xx_yyy[i] = 2.0 * ts_x_yyy[i] * gfe2_0 + 2.0 * gr_x_yyy[i] * gfe_0 + ts_xx_yyy[i] * gfe_0 * gc_x[i] + gr_xx_yyy[i] * gc_x[i];

        grr_x_xx_yyz[i] = 2.0 * ts_x_yyz[i] * gfe2_0 + 2.0 * gr_x_yyz[i] * gfe_0 + ts_xx_yyz[i] * gfe_0 * gc_x[i] + gr_xx_yyz[i] * gc_x[i];

        grr_x_xx_yzz[i] = 2.0 * ts_x_yzz[i] * gfe2_0 + 2.0 * gr_x_yzz[i] * gfe_0 + ts_xx_yzz[i] * gfe_0 * gc_x[i] + gr_xx_yzz[i] * gc_x[i];

        grr_x_xx_zzz[i] = 2.0 * ts_x_zzz[i] * gfe2_0 + 2.0 * gr_x_zzz[i] * gfe_0 + ts_xx_zzz[i] * gfe_0 * gc_x[i] + gr_xx_zzz[i] * gc_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto grr_x_xy_xxx = pbuffer.data(idx_gr_df + 10);

    auto grr_x_xy_xxy = pbuffer.data(idx_gr_df + 11);

    auto grr_x_xy_xxz = pbuffer.data(idx_gr_df + 12);

    auto grr_x_xy_xyy = pbuffer.data(idx_gr_df + 13);

    auto grr_x_xy_xyz = pbuffer.data(idx_gr_df + 14);

    auto grr_x_xy_xzz = pbuffer.data(idx_gr_df + 15);

    auto grr_x_xy_yyy = pbuffer.data(idx_gr_df + 16);

    auto grr_x_xy_yyz = pbuffer.data(idx_gr_df + 17);

    auto grr_x_xy_yzz = pbuffer.data(idx_gr_df + 18);

    auto grr_x_xy_zzz = pbuffer.data(idx_gr_df + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xx, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xy, gr_xy_xyy, gr_xy_xyz, gr_xy_xz, gr_xy_xzz, gr_xy_yy, gr_xy_yyy, gr_xy_yyz, gr_xy_yz, gr_xy_yzz, gr_xy_zz, gr_xy_zzz, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xyy, gr_y_xyz, gr_y_xzz, gr_y_yyy, gr_y_yyz, gr_y_yzz, gr_y_zzz, grr_x_xy_xxx, grr_x_xy_xxy, grr_x_xy_xxz, grr_x_xy_xyy, grr_x_xy_xyz, grr_x_xy_xzz, grr_x_xy_yyy, grr_x_xy_yyz, grr_x_xy_yzz, grr_x_xy_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xy_xxx[i] = ts_y_xxx[i] * gfe2_0 + gr_y_xxx[i] * gfe_0 + 3.0 * ts_xy_xx[i] * gfe2_0 + 3.0 * gr_xy_xx[i] * gfe_0 + ts_xy_xxx[i] * gfe_0 * gc_x[i] + gr_xy_xxx[i] * gc_x[i];

        grr_x_xy_xxy[i] = ts_y_xxy[i] * gfe2_0 + gr_y_xxy[i] * gfe_0 + 2.0 * ts_xy_xy[i] * gfe2_0 + 2.0 * gr_xy_xy[i] * gfe_0 + ts_xy_xxy[i] * gfe_0 * gc_x[i] + gr_xy_xxy[i] * gc_x[i];

        grr_x_xy_xxz[i] = ts_y_xxz[i] * gfe2_0 + gr_y_xxz[i] * gfe_0 + 2.0 * ts_xy_xz[i] * gfe2_0 + 2.0 * gr_xy_xz[i] * gfe_0 + ts_xy_xxz[i] * gfe_0 * gc_x[i] + gr_xy_xxz[i] * gc_x[i];

        grr_x_xy_xyy[i] = ts_y_xyy[i] * gfe2_0 + gr_y_xyy[i] * gfe_0 + ts_xy_yy[i] * gfe2_0 + gr_xy_yy[i] * gfe_0 + ts_xy_xyy[i] * gfe_0 * gc_x[i] + gr_xy_xyy[i] * gc_x[i];

        grr_x_xy_xyz[i] = ts_y_xyz[i] * gfe2_0 + gr_y_xyz[i] * gfe_0 + ts_xy_yz[i] * gfe2_0 + gr_xy_yz[i] * gfe_0 + ts_xy_xyz[i] * gfe_0 * gc_x[i] + gr_xy_xyz[i] * gc_x[i];

        grr_x_xy_xzz[i] = ts_y_xzz[i] * gfe2_0 + gr_y_xzz[i] * gfe_0 + ts_xy_zz[i] * gfe2_0 + gr_xy_zz[i] * gfe_0 + ts_xy_xzz[i] * gfe_0 * gc_x[i] + gr_xy_xzz[i] * gc_x[i];

        grr_x_xy_yyy[i] = ts_y_yyy[i] * gfe2_0 + gr_y_yyy[i] * gfe_0 + ts_xy_yyy[i] * gfe_0 * gc_x[i] + gr_xy_yyy[i] * gc_x[i];

        grr_x_xy_yyz[i] = ts_y_yyz[i] * gfe2_0 + gr_y_yyz[i] * gfe_0 + ts_xy_yyz[i] * gfe_0 * gc_x[i] + gr_xy_yyz[i] * gc_x[i];

        grr_x_xy_yzz[i] = ts_y_yzz[i] * gfe2_0 + gr_y_yzz[i] * gfe_0 + ts_xy_yzz[i] * gfe_0 * gc_x[i] + gr_xy_yzz[i] * gc_x[i];

        grr_x_xy_zzz[i] = ts_y_zzz[i] * gfe2_0 + gr_y_zzz[i] * gfe_0 + ts_xy_zzz[i] * gfe_0 * gc_x[i] + gr_xy_zzz[i] * gc_x[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto grr_x_xz_xxx = pbuffer.data(idx_gr_df + 20);

    auto grr_x_xz_xxy = pbuffer.data(idx_gr_df + 21);

    auto grr_x_xz_xxz = pbuffer.data(idx_gr_df + 22);

    auto grr_x_xz_xyy = pbuffer.data(idx_gr_df + 23);

    auto grr_x_xz_xyz = pbuffer.data(idx_gr_df + 24);

    auto grr_x_xz_xzz = pbuffer.data(idx_gr_df + 25);

    auto grr_x_xz_yyy = pbuffer.data(idx_gr_df + 26);

    auto grr_x_xz_yyz = pbuffer.data(idx_gr_df + 27);

    auto grr_x_xz_yzz = pbuffer.data(idx_gr_df + 28);

    auto grr_x_xz_zzz = pbuffer.data(idx_gr_df + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xx, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xy, gr_xz_xyy, gr_xz_xyz, gr_xz_xz, gr_xz_xzz, gr_xz_yy, gr_xz_yyy, gr_xz_yyz, gr_xz_yz, gr_xz_yzz, gr_xz_zz, gr_xz_zzz, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xyy, gr_z_xyz, gr_z_xzz, gr_z_yyy, gr_z_yyz, gr_z_yzz, gr_z_zzz, grr_x_xz_xxx, grr_x_xz_xxy, grr_x_xz_xxz, grr_x_xz_xyy, grr_x_xz_xyz, grr_x_xz_xzz, grr_x_xz_yyy, grr_x_xz_yyz, grr_x_xz_yzz, grr_x_xz_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xz_xxx[i] = ts_z_xxx[i] * gfe2_0 + gr_z_xxx[i] * gfe_0 + 3.0 * ts_xz_xx[i] * gfe2_0 + 3.0 * gr_xz_xx[i] * gfe_0 + ts_xz_xxx[i] * gfe_0 * gc_x[i] + gr_xz_xxx[i] * gc_x[i];

        grr_x_xz_xxy[i] = ts_z_xxy[i] * gfe2_0 + gr_z_xxy[i] * gfe_0 + 2.0 * ts_xz_xy[i] * gfe2_0 + 2.0 * gr_xz_xy[i] * gfe_0 + ts_xz_xxy[i] * gfe_0 * gc_x[i] + gr_xz_xxy[i] * gc_x[i];

        grr_x_xz_xxz[i] = ts_z_xxz[i] * gfe2_0 + gr_z_xxz[i] * gfe_0 + 2.0 * ts_xz_xz[i] * gfe2_0 + 2.0 * gr_xz_xz[i] * gfe_0 + ts_xz_xxz[i] * gfe_0 * gc_x[i] + gr_xz_xxz[i] * gc_x[i];

        grr_x_xz_xyy[i] = ts_z_xyy[i] * gfe2_0 + gr_z_xyy[i] * gfe_0 + ts_xz_yy[i] * gfe2_0 + gr_xz_yy[i] * gfe_0 + ts_xz_xyy[i] * gfe_0 * gc_x[i] + gr_xz_xyy[i] * gc_x[i];

        grr_x_xz_xyz[i] = ts_z_xyz[i] * gfe2_0 + gr_z_xyz[i] * gfe_0 + ts_xz_yz[i] * gfe2_0 + gr_xz_yz[i] * gfe_0 + ts_xz_xyz[i] * gfe_0 * gc_x[i] + gr_xz_xyz[i] * gc_x[i];

        grr_x_xz_xzz[i] = ts_z_xzz[i] * gfe2_0 + gr_z_xzz[i] * gfe_0 + ts_xz_zz[i] * gfe2_0 + gr_xz_zz[i] * gfe_0 + ts_xz_xzz[i] * gfe_0 * gc_x[i] + gr_xz_xzz[i] * gc_x[i];

        grr_x_xz_yyy[i] = ts_z_yyy[i] * gfe2_0 + gr_z_yyy[i] * gfe_0 + ts_xz_yyy[i] * gfe_0 * gc_x[i] + gr_xz_yyy[i] * gc_x[i];

        grr_x_xz_yyz[i] = ts_z_yyz[i] * gfe2_0 + gr_z_yyz[i] * gfe_0 + ts_xz_yyz[i] * gfe_0 * gc_x[i] + gr_xz_yyz[i] * gc_x[i];

        grr_x_xz_yzz[i] = ts_z_yzz[i] * gfe2_0 + gr_z_yzz[i] * gfe_0 + ts_xz_yzz[i] * gfe_0 * gc_x[i] + gr_xz_yzz[i] * gc_x[i];

        grr_x_xz_zzz[i] = ts_z_zzz[i] * gfe2_0 + gr_z_zzz[i] * gfe_0 + ts_xz_zzz[i] * gfe_0 * gc_x[i] + gr_xz_zzz[i] * gc_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto grr_x_yy_xxx = pbuffer.data(idx_gr_df + 30);

    auto grr_x_yy_xxy = pbuffer.data(idx_gr_df + 31);

    auto grr_x_yy_xxz = pbuffer.data(idx_gr_df + 32);

    auto grr_x_yy_xyy = pbuffer.data(idx_gr_df + 33);

    auto grr_x_yy_xyz = pbuffer.data(idx_gr_df + 34);

    auto grr_x_yy_xzz = pbuffer.data(idx_gr_df + 35);

    auto grr_x_yy_yyy = pbuffer.data(idx_gr_df + 36);

    auto grr_x_yy_yyz = pbuffer.data(idx_gr_df + 37);

    auto grr_x_yy_yzz = pbuffer.data(idx_gr_df + 38);

    auto grr_x_yy_zzz = pbuffer.data(idx_gr_df + 39);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xx, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xy, gr_yy_xyy, gr_yy_xyz, gr_yy_xz, gr_yy_xzz, gr_yy_yy, gr_yy_yyy, gr_yy_yyz, gr_yy_yz, gr_yy_yzz, gr_yy_zz, gr_yy_zzz, grr_x_yy_xxx, grr_x_yy_xxy, grr_x_yy_xxz, grr_x_yy_xyy, grr_x_yy_xyz, grr_x_yy_xzz, grr_x_yy_yyy, grr_x_yy_yyz, grr_x_yy_yzz, grr_x_yy_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yy_xxx[i] = 3.0 * ts_yy_xx[i] * gfe2_0 + 3.0 * gr_yy_xx[i] * gfe_0 + ts_yy_xxx[i] * gfe_0 * gc_x[i] + gr_yy_xxx[i] * gc_x[i];

        grr_x_yy_xxy[i] = 2.0 * ts_yy_xy[i] * gfe2_0 + 2.0 * gr_yy_xy[i] * gfe_0 + ts_yy_xxy[i] * gfe_0 * gc_x[i] + gr_yy_xxy[i] * gc_x[i];

        grr_x_yy_xxz[i] = 2.0 * ts_yy_xz[i] * gfe2_0 + 2.0 * gr_yy_xz[i] * gfe_0 + ts_yy_xxz[i] * gfe_0 * gc_x[i] + gr_yy_xxz[i] * gc_x[i];

        grr_x_yy_xyy[i] = ts_yy_yy[i] * gfe2_0 + gr_yy_yy[i] * gfe_0 + ts_yy_xyy[i] * gfe_0 * gc_x[i] + gr_yy_xyy[i] * gc_x[i];

        grr_x_yy_xyz[i] = ts_yy_yz[i] * gfe2_0 + gr_yy_yz[i] * gfe_0 + ts_yy_xyz[i] * gfe_0 * gc_x[i] + gr_yy_xyz[i] * gc_x[i];

        grr_x_yy_xzz[i] = ts_yy_zz[i] * gfe2_0 + gr_yy_zz[i] * gfe_0 + ts_yy_xzz[i] * gfe_0 * gc_x[i] + gr_yy_xzz[i] * gc_x[i];

        grr_x_yy_yyy[i] = ts_yy_yyy[i] * gfe_0 * gc_x[i] + gr_yy_yyy[i] * gc_x[i];

        grr_x_yy_yyz[i] = ts_yy_yyz[i] * gfe_0 * gc_x[i] + gr_yy_yyz[i] * gc_x[i];

        grr_x_yy_yzz[i] = ts_yy_yzz[i] * gfe_0 * gc_x[i] + gr_yy_yzz[i] * gc_x[i];

        grr_x_yy_zzz[i] = ts_yy_zzz[i] * gfe_0 * gc_x[i] + gr_yy_zzz[i] * gc_x[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto grr_x_yz_xxx = pbuffer.data(idx_gr_df + 40);

    auto grr_x_yz_xxy = pbuffer.data(idx_gr_df + 41);

    auto grr_x_yz_xxz = pbuffer.data(idx_gr_df + 42);

    auto grr_x_yz_xyy = pbuffer.data(idx_gr_df + 43);

    auto grr_x_yz_xyz = pbuffer.data(idx_gr_df + 44);

    auto grr_x_yz_xzz = pbuffer.data(idx_gr_df + 45);

    auto grr_x_yz_yyy = pbuffer.data(idx_gr_df + 46);

    auto grr_x_yz_yyz = pbuffer.data(idx_gr_df + 47);

    auto grr_x_yz_yzz = pbuffer.data(idx_gr_df + 48);

    auto grr_x_yz_zzz = pbuffer.data(idx_gr_df + 49);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xx, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xy, gr_yz_xyy, gr_yz_xyz, gr_yz_xz, gr_yz_xzz, gr_yz_yy, gr_yz_yyy, gr_yz_yyz, gr_yz_yz, gr_yz_yzz, gr_yz_zz, gr_yz_zzz, grr_x_yz_xxx, grr_x_yz_xxy, grr_x_yz_xxz, grr_x_yz_xyy, grr_x_yz_xyz, grr_x_yz_xzz, grr_x_yz_yyy, grr_x_yz_yyz, grr_x_yz_yzz, grr_x_yz_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yz_xxx[i] = 3.0 * ts_yz_xx[i] * gfe2_0 + 3.0 * gr_yz_xx[i] * gfe_0 + ts_yz_xxx[i] * gfe_0 * gc_x[i] + gr_yz_xxx[i] * gc_x[i];

        grr_x_yz_xxy[i] = 2.0 * ts_yz_xy[i] * gfe2_0 + 2.0 * gr_yz_xy[i] * gfe_0 + ts_yz_xxy[i] * gfe_0 * gc_x[i] + gr_yz_xxy[i] * gc_x[i];

        grr_x_yz_xxz[i] = 2.0 * ts_yz_xz[i] * gfe2_0 + 2.0 * gr_yz_xz[i] * gfe_0 + ts_yz_xxz[i] * gfe_0 * gc_x[i] + gr_yz_xxz[i] * gc_x[i];

        grr_x_yz_xyy[i] = ts_yz_yy[i] * gfe2_0 + gr_yz_yy[i] * gfe_0 + ts_yz_xyy[i] * gfe_0 * gc_x[i] + gr_yz_xyy[i] * gc_x[i];

        grr_x_yz_xyz[i] = ts_yz_yz[i] * gfe2_0 + gr_yz_yz[i] * gfe_0 + ts_yz_xyz[i] * gfe_0 * gc_x[i] + gr_yz_xyz[i] * gc_x[i];

        grr_x_yz_xzz[i] = ts_yz_zz[i] * gfe2_0 + gr_yz_zz[i] * gfe_0 + ts_yz_xzz[i] * gfe_0 * gc_x[i] + gr_yz_xzz[i] * gc_x[i];

        grr_x_yz_yyy[i] = ts_yz_yyy[i] * gfe_0 * gc_x[i] + gr_yz_yyy[i] * gc_x[i];

        grr_x_yz_yyz[i] = ts_yz_yyz[i] * gfe_0 * gc_x[i] + gr_yz_yyz[i] * gc_x[i];

        grr_x_yz_yzz[i] = ts_yz_yzz[i] * gfe_0 * gc_x[i] + gr_yz_yzz[i] * gc_x[i];

        grr_x_yz_zzz[i] = ts_yz_zzz[i] * gfe_0 * gc_x[i] + gr_yz_zzz[i] * gc_x[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto grr_x_zz_xxx = pbuffer.data(idx_gr_df + 50);

    auto grr_x_zz_xxy = pbuffer.data(idx_gr_df + 51);

    auto grr_x_zz_xxz = pbuffer.data(idx_gr_df + 52);

    auto grr_x_zz_xyy = pbuffer.data(idx_gr_df + 53);

    auto grr_x_zz_xyz = pbuffer.data(idx_gr_df + 54);

    auto grr_x_zz_xzz = pbuffer.data(idx_gr_df + 55);

    auto grr_x_zz_yyy = pbuffer.data(idx_gr_df + 56);

    auto grr_x_zz_yyz = pbuffer.data(idx_gr_df + 57);

    auto grr_x_zz_yzz = pbuffer.data(idx_gr_df + 58);

    auto grr_x_zz_zzz = pbuffer.data(idx_gr_df + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xx, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xy, gr_zz_xyy, gr_zz_xyz, gr_zz_xz, gr_zz_xzz, gr_zz_yy, gr_zz_yyy, gr_zz_yyz, gr_zz_yz, gr_zz_yzz, gr_zz_zz, gr_zz_zzz, grr_x_zz_xxx, grr_x_zz_xxy, grr_x_zz_xxz, grr_x_zz_xyy, grr_x_zz_xyz, grr_x_zz_xzz, grr_x_zz_yyy, grr_x_zz_yyz, grr_x_zz_yzz, grr_x_zz_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zz_xxx[i] = 3.0 * ts_zz_xx[i] * gfe2_0 + 3.0 * gr_zz_xx[i] * gfe_0 + ts_zz_xxx[i] * gfe_0 * gc_x[i] + gr_zz_xxx[i] * gc_x[i];

        grr_x_zz_xxy[i] = 2.0 * ts_zz_xy[i] * gfe2_0 + 2.0 * gr_zz_xy[i] * gfe_0 + ts_zz_xxy[i] * gfe_0 * gc_x[i] + gr_zz_xxy[i] * gc_x[i];

        grr_x_zz_xxz[i] = 2.0 * ts_zz_xz[i] * gfe2_0 + 2.0 * gr_zz_xz[i] * gfe_0 + ts_zz_xxz[i] * gfe_0 * gc_x[i] + gr_zz_xxz[i] * gc_x[i];

        grr_x_zz_xyy[i] = ts_zz_yy[i] * gfe2_0 + gr_zz_yy[i] * gfe_0 + ts_zz_xyy[i] * gfe_0 * gc_x[i] + gr_zz_xyy[i] * gc_x[i];

        grr_x_zz_xyz[i] = ts_zz_yz[i] * gfe2_0 + gr_zz_yz[i] * gfe_0 + ts_zz_xyz[i] * gfe_0 * gc_x[i] + gr_zz_xyz[i] * gc_x[i];

        grr_x_zz_xzz[i] = ts_zz_zz[i] * gfe2_0 + gr_zz_zz[i] * gfe_0 + ts_zz_xzz[i] * gfe_0 * gc_x[i] + gr_zz_xzz[i] * gc_x[i];

        grr_x_zz_yyy[i] = ts_zz_yyy[i] * gfe_0 * gc_x[i] + gr_zz_yyy[i] * gc_x[i];

        grr_x_zz_yyz[i] = ts_zz_yyz[i] * gfe_0 * gc_x[i] + gr_zz_yyz[i] * gc_x[i];

        grr_x_zz_yzz[i] = ts_zz_yzz[i] * gfe_0 * gc_x[i] + gr_zz_yzz[i] * gc_x[i];

        grr_x_zz_zzz[i] = ts_zz_zzz[i] * gfe_0 * gc_x[i] + gr_zz_zzz[i] * gc_x[i];
    }

    // Set up 60-70 components of targeted buffer : DF

    auto grr_y_xx_xxx = pbuffer.data(idx_gr_df + 60);

    auto grr_y_xx_xxy = pbuffer.data(idx_gr_df + 61);

    auto grr_y_xx_xxz = pbuffer.data(idx_gr_df + 62);

    auto grr_y_xx_xyy = pbuffer.data(idx_gr_df + 63);

    auto grr_y_xx_xyz = pbuffer.data(idx_gr_df + 64);

    auto grr_y_xx_xzz = pbuffer.data(idx_gr_df + 65);

    auto grr_y_xx_yyy = pbuffer.data(idx_gr_df + 66);

    auto grr_y_xx_yyz = pbuffer.data(idx_gr_df + 67);

    auto grr_y_xx_yzz = pbuffer.data(idx_gr_df + 68);

    auto grr_y_xx_zzz = pbuffer.data(idx_gr_df + 69);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xx, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xy, gr_xx_xyy, gr_xx_xyz, gr_xx_xz, gr_xx_xzz, gr_xx_yy, gr_xx_yyy, gr_xx_yyz, gr_xx_yz, gr_xx_yzz, gr_xx_zz, gr_xx_zzz, grr_y_xx_xxx, grr_y_xx_xxy, grr_y_xx_xxz, grr_y_xx_xyy, grr_y_xx_xyz, grr_y_xx_xzz, grr_y_xx_yyy, grr_y_xx_yyz, grr_y_xx_yzz, grr_y_xx_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xx_xxx[i] = ts_xx_xxx[i] * gfe_0 * gc_y[i] + gr_xx_xxx[i] * gc_y[i];

        grr_y_xx_xxy[i] = ts_xx_xx[i] * gfe2_0 + gr_xx_xx[i] * gfe_0 + ts_xx_xxy[i] * gfe_0 * gc_y[i] + gr_xx_xxy[i] * gc_y[i];

        grr_y_xx_xxz[i] = ts_xx_xxz[i] * gfe_0 * gc_y[i] + gr_xx_xxz[i] * gc_y[i];

        grr_y_xx_xyy[i] = 2.0 * ts_xx_xy[i] * gfe2_0 + 2.0 * gr_xx_xy[i] * gfe_0 + ts_xx_xyy[i] * gfe_0 * gc_y[i] + gr_xx_xyy[i] * gc_y[i];

        grr_y_xx_xyz[i] = ts_xx_xz[i] * gfe2_0 + gr_xx_xz[i] * gfe_0 + ts_xx_xyz[i] * gfe_0 * gc_y[i] + gr_xx_xyz[i] * gc_y[i];

        grr_y_xx_xzz[i] = ts_xx_xzz[i] * gfe_0 * gc_y[i] + gr_xx_xzz[i] * gc_y[i];

        grr_y_xx_yyy[i] = 3.0 * ts_xx_yy[i] * gfe2_0 + 3.0 * gr_xx_yy[i] * gfe_0 + ts_xx_yyy[i] * gfe_0 * gc_y[i] + gr_xx_yyy[i] * gc_y[i];

        grr_y_xx_yyz[i] = 2.0 * ts_xx_yz[i] * gfe2_0 + 2.0 * gr_xx_yz[i] * gfe_0 + ts_xx_yyz[i] * gfe_0 * gc_y[i] + gr_xx_yyz[i] * gc_y[i];

        grr_y_xx_yzz[i] = ts_xx_zz[i] * gfe2_0 + gr_xx_zz[i] * gfe_0 + ts_xx_yzz[i] * gfe_0 * gc_y[i] + gr_xx_yzz[i] * gc_y[i];

        grr_y_xx_zzz[i] = ts_xx_zzz[i] * gfe_0 * gc_y[i] + gr_xx_zzz[i] * gc_y[i];
    }

    // Set up 70-80 components of targeted buffer : DF

    auto grr_y_xy_xxx = pbuffer.data(idx_gr_df + 70);

    auto grr_y_xy_xxy = pbuffer.data(idx_gr_df + 71);

    auto grr_y_xy_xxz = pbuffer.data(idx_gr_df + 72);

    auto grr_y_xy_xyy = pbuffer.data(idx_gr_df + 73);

    auto grr_y_xy_xyz = pbuffer.data(idx_gr_df + 74);

    auto grr_y_xy_xzz = pbuffer.data(idx_gr_df + 75);

    auto grr_y_xy_yyy = pbuffer.data(idx_gr_df + 76);

    auto grr_y_xy_yyz = pbuffer.data(idx_gr_df + 77);

    auto grr_y_xy_yzz = pbuffer.data(idx_gr_df + 78);

    auto grr_y_xy_zzz = pbuffer.data(idx_gr_df + 79);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xyy, gr_x_xyz, gr_x_xzz, gr_x_yyy, gr_x_yyz, gr_x_yzz, gr_x_zzz, gr_xy_xx, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xy, gr_xy_xyy, gr_xy_xyz, gr_xy_xz, gr_xy_xzz, gr_xy_yy, gr_xy_yyy, gr_xy_yyz, gr_xy_yz, gr_xy_yzz, gr_xy_zz, gr_xy_zzz, grr_y_xy_xxx, grr_y_xy_xxy, grr_y_xy_xxz, grr_y_xy_xyy, grr_y_xy_xyz, grr_y_xy_xzz, grr_y_xy_yyy, grr_y_xy_yyz, grr_y_xy_yzz, grr_y_xy_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xy_xxx[i] = ts_x_xxx[i] * gfe2_0 + gr_x_xxx[i] * gfe_0 + ts_xy_xxx[i] * gfe_0 * gc_y[i] + gr_xy_xxx[i] * gc_y[i];

        grr_y_xy_xxy[i] = ts_x_xxy[i] * gfe2_0 + gr_x_xxy[i] * gfe_0 + ts_xy_xx[i] * gfe2_0 + gr_xy_xx[i] * gfe_0 + ts_xy_xxy[i] * gfe_0 * gc_y[i] + gr_xy_xxy[i] * gc_y[i];

        grr_y_xy_xxz[i] = ts_x_xxz[i] * gfe2_0 + gr_x_xxz[i] * gfe_0 + ts_xy_xxz[i] * gfe_0 * gc_y[i] + gr_xy_xxz[i] * gc_y[i];

        grr_y_xy_xyy[i] = ts_x_xyy[i] * gfe2_0 + gr_x_xyy[i] * gfe_0 + 2.0 * ts_xy_xy[i] * gfe2_0 + 2.0 * gr_xy_xy[i] * gfe_0 + ts_xy_xyy[i] * gfe_0 * gc_y[i] + gr_xy_xyy[i] * gc_y[i];

        grr_y_xy_xyz[i] = ts_x_xyz[i] * gfe2_0 + gr_x_xyz[i] * gfe_0 + ts_xy_xz[i] * gfe2_0 + gr_xy_xz[i] * gfe_0 + ts_xy_xyz[i] * gfe_0 * gc_y[i] + gr_xy_xyz[i] * gc_y[i];

        grr_y_xy_xzz[i] = ts_x_xzz[i] * gfe2_0 + gr_x_xzz[i] * gfe_0 + ts_xy_xzz[i] * gfe_0 * gc_y[i] + gr_xy_xzz[i] * gc_y[i];

        grr_y_xy_yyy[i] = ts_x_yyy[i] * gfe2_0 + gr_x_yyy[i] * gfe_0 + 3.0 * ts_xy_yy[i] * gfe2_0 + 3.0 * gr_xy_yy[i] * gfe_0 + ts_xy_yyy[i] * gfe_0 * gc_y[i] + gr_xy_yyy[i] * gc_y[i];

        grr_y_xy_yyz[i] = ts_x_yyz[i] * gfe2_0 + gr_x_yyz[i] * gfe_0 + 2.0 * ts_xy_yz[i] * gfe2_0 + 2.0 * gr_xy_yz[i] * gfe_0 + ts_xy_yyz[i] * gfe_0 * gc_y[i] + gr_xy_yyz[i] * gc_y[i];

        grr_y_xy_yzz[i] = ts_x_yzz[i] * gfe2_0 + gr_x_yzz[i] * gfe_0 + ts_xy_zz[i] * gfe2_0 + gr_xy_zz[i] * gfe_0 + ts_xy_yzz[i] * gfe_0 * gc_y[i] + gr_xy_yzz[i] * gc_y[i];

        grr_y_xy_zzz[i] = ts_x_zzz[i] * gfe2_0 + gr_x_zzz[i] * gfe_0 + ts_xy_zzz[i] * gfe_0 * gc_y[i] + gr_xy_zzz[i] * gc_y[i];
    }

    // Set up 80-90 components of targeted buffer : DF

    auto grr_y_xz_xxx = pbuffer.data(idx_gr_df + 80);

    auto grr_y_xz_xxy = pbuffer.data(idx_gr_df + 81);

    auto grr_y_xz_xxz = pbuffer.data(idx_gr_df + 82);

    auto grr_y_xz_xyy = pbuffer.data(idx_gr_df + 83);

    auto grr_y_xz_xyz = pbuffer.data(idx_gr_df + 84);

    auto grr_y_xz_xzz = pbuffer.data(idx_gr_df + 85);

    auto grr_y_xz_yyy = pbuffer.data(idx_gr_df + 86);

    auto grr_y_xz_yyz = pbuffer.data(idx_gr_df + 87);

    auto grr_y_xz_yzz = pbuffer.data(idx_gr_df + 88);

    auto grr_y_xz_zzz = pbuffer.data(idx_gr_df + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xx, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xy, gr_xz_xyy, gr_xz_xyz, gr_xz_xz, gr_xz_xzz, gr_xz_yy, gr_xz_yyy, gr_xz_yyz, gr_xz_yz, gr_xz_yzz, gr_xz_zz, gr_xz_zzz, grr_y_xz_xxx, grr_y_xz_xxy, grr_y_xz_xxz, grr_y_xz_xyy, grr_y_xz_xyz, grr_y_xz_xzz, grr_y_xz_yyy, grr_y_xz_yyz, grr_y_xz_yzz, grr_y_xz_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xz_xxx[i] = ts_xz_xxx[i] * gfe_0 * gc_y[i] + gr_xz_xxx[i] * gc_y[i];

        grr_y_xz_xxy[i] = ts_xz_xx[i] * gfe2_0 + gr_xz_xx[i] * gfe_0 + ts_xz_xxy[i] * gfe_0 * gc_y[i] + gr_xz_xxy[i] * gc_y[i];

        grr_y_xz_xxz[i] = ts_xz_xxz[i] * gfe_0 * gc_y[i] + gr_xz_xxz[i] * gc_y[i];

        grr_y_xz_xyy[i] = 2.0 * ts_xz_xy[i] * gfe2_0 + 2.0 * gr_xz_xy[i] * gfe_0 + ts_xz_xyy[i] * gfe_0 * gc_y[i] + gr_xz_xyy[i] * gc_y[i];

        grr_y_xz_xyz[i] = ts_xz_xz[i] * gfe2_0 + gr_xz_xz[i] * gfe_0 + ts_xz_xyz[i] * gfe_0 * gc_y[i] + gr_xz_xyz[i] * gc_y[i];

        grr_y_xz_xzz[i] = ts_xz_xzz[i] * gfe_0 * gc_y[i] + gr_xz_xzz[i] * gc_y[i];

        grr_y_xz_yyy[i] = 3.0 * ts_xz_yy[i] * gfe2_0 + 3.0 * gr_xz_yy[i] * gfe_0 + ts_xz_yyy[i] * gfe_0 * gc_y[i] + gr_xz_yyy[i] * gc_y[i];

        grr_y_xz_yyz[i] = 2.0 * ts_xz_yz[i] * gfe2_0 + 2.0 * gr_xz_yz[i] * gfe_0 + ts_xz_yyz[i] * gfe_0 * gc_y[i] + gr_xz_yyz[i] * gc_y[i];

        grr_y_xz_yzz[i] = ts_xz_zz[i] * gfe2_0 + gr_xz_zz[i] * gfe_0 + ts_xz_yzz[i] * gfe_0 * gc_y[i] + gr_xz_yzz[i] * gc_y[i];

        grr_y_xz_zzz[i] = ts_xz_zzz[i] * gfe_0 * gc_y[i] + gr_xz_zzz[i] * gc_y[i];
    }

    // Set up 90-100 components of targeted buffer : DF

    auto grr_y_yy_xxx = pbuffer.data(idx_gr_df + 90);

    auto grr_y_yy_xxy = pbuffer.data(idx_gr_df + 91);

    auto grr_y_yy_xxz = pbuffer.data(idx_gr_df + 92);

    auto grr_y_yy_xyy = pbuffer.data(idx_gr_df + 93);

    auto grr_y_yy_xyz = pbuffer.data(idx_gr_df + 94);

    auto grr_y_yy_xzz = pbuffer.data(idx_gr_df + 95);

    auto grr_y_yy_yyy = pbuffer.data(idx_gr_df + 96);

    auto grr_y_yy_yyz = pbuffer.data(idx_gr_df + 97);

    auto grr_y_yy_yzz = pbuffer.data(idx_gr_df + 98);

    auto grr_y_yy_zzz = pbuffer.data(idx_gr_df + 99);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xyy, gr_y_xyz, gr_y_xzz, gr_y_yyy, gr_y_yyz, gr_y_yzz, gr_y_zzz, gr_yy_xx, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xy, gr_yy_xyy, gr_yy_xyz, gr_yy_xz, gr_yy_xzz, gr_yy_yy, gr_yy_yyy, gr_yy_yyz, gr_yy_yz, gr_yy_yzz, gr_yy_zz, gr_yy_zzz, grr_y_yy_xxx, grr_y_yy_xxy, grr_y_yy_xxz, grr_y_yy_xyy, grr_y_yy_xyz, grr_y_yy_xzz, grr_y_yy_yyy, grr_y_yy_yyz, grr_y_yy_yzz, grr_y_yy_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yy_xxx[i] = 2.0 * ts_y_xxx[i] * gfe2_0 + 2.0 * gr_y_xxx[i] * gfe_0 + ts_yy_xxx[i] * gfe_0 * gc_y[i] + gr_yy_xxx[i] * gc_y[i];

        grr_y_yy_xxy[i] = 2.0 * ts_y_xxy[i] * gfe2_0 + 2.0 * gr_y_xxy[i] * gfe_0 + ts_yy_xx[i] * gfe2_0 + gr_yy_xx[i] * gfe_0 + ts_yy_xxy[i] * gfe_0 * gc_y[i] + gr_yy_xxy[i] * gc_y[i];

        grr_y_yy_xxz[i] = 2.0 * ts_y_xxz[i] * gfe2_0 + 2.0 * gr_y_xxz[i] * gfe_0 + ts_yy_xxz[i] * gfe_0 * gc_y[i] + gr_yy_xxz[i] * gc_y[i];

        grr_y_yy_xyy[i] = 2.0 * ts_y_xyy[i] * gfe2_0 + 2.0 * gr_y_xyy[i] * gfe_0 + 2.0 * ts_yy_xy[i] * gfe2_0 + 2.0 * gr_yy_xy[i] * gfe_0 + ts_yy_xyy[i] * gfe_0 * gc_y[i] + gr_yy_xyy[i] * gc_y[i];

        grr_y_yy_xyz[i] = 2.0 * ts_y_xyz[i] * gfe2_0 + 2.0 * gr_y_xyz[i] * gfe_0 + ts_yy_xz[i] * gfe2_0 + gr_yy_xz[i] * gfe_0 + ts_yy_xyz[i] * gfe_0 * gc_y[i] + gr_yy_xyz[i] * gc_y[i];

        grr_y_yy_xzz[i] = 2.0 * ts_y_xzz[i] * gfe2_0 + 2.0 * gr_y_xzz[i] * gfe_0 + ts_yy_xzz[i] * gfe_0 * gc_y[i] + gr_yy_xzz[i] * gc_y[i];

        grr_y_yy_yyy[i] = 2.0 * ts_y_yyy[i] * gfe2_0 + 2.0 * gr_y_yyy[i] * gfe_0 + 3.0 * ts_yy_yy[i] * gfe2_0 + 3.0 * gr_yy_yy[i] * gfe_0 + ts_yy_yyy[i] * gfe_0 * gc_y[i] + gr_yy_yyy[i] * gc_y[i];

        grr_y_yy_yyz[i] = 2.0 * ts_y_yyz[i] * gfe2_0 + 2.0 * gr_y_yyz[i] * gfe_0 + 2.0 * ts_yy_yz[i] * gfe2_0 + 2.0 * gr_yy_yz[i] * gfe_0 + ts_yy_yyz[i] * gfe_0 * gc_y[i] + gr_yy_yyz[i] * gc_y[i];

        grr_y_yy_yzz[i] = 2.0 * ts_y_yzz[i] * gfe2_0 + 2.0 * gr_y_yzz[i] * gfe_0 + ts_yy_zz[i] * gfe2_0 + gr_yy_zz[i] * gfe_0 + ts_yy_yzz[i] * gfe_0 * gc_y[i] + gr_yy_yzz[i] * gc_y[i];

        grr_y_yy_zzz[i] = 2.0 * ts_y_zzz[i] * gfe2_0 + 2.0 * gr_y_zzz[i] * gfe_0 + ts_yy_zzz[i] * gfe_0 * gc_y[i] + gr_yy_zzz[i] * gc_y[i];
    }

    // Set up 100-110 components of targeted buffer : DF

    auto grr_y_yz_xxx = pbuffer.data(idx_gr_df + 100);

    auto grr_y_yz_xxy = pbuffer.data(idx_gr_df + 101);

    auto grr_y_yz_xxz = pbuffer.data(idx_gr_df + 102);

    auto grr_y_yz_xyy = pbuffer.data(idx_gr_df + 103);

    auto grr_y_yz_xyz = pbuffer.data(idx_gr_df + 104);

    auto grr_y_yz_xzz = pbuffer.data(idx_gr_df + 105);

    auto grr_y_yz_yyy = pbuffer.data(idx_gr_df + 106);

    auto grr_y_yz_yyz = pbuffer.data(idx_gr_df + 107);

    auto grr_y_yz_yzz = pbuffer.data(idx_gr_df + 108);

    auto grr_y_yz_zzz = pbuffer.data(idx_gr_df + 109);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xx, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xy, gr_yz_xyy, gr_yz_xyz, gr_yz_xz, gr_yz_xzz, gr_yz_yy, gr_yz_yyy, gr_yz_yyz, gr_yz_yz, gr_yz_yzz, gr_yz_zz, gr_yz_zzz, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xyy, gr_z_xyz, gr_z_xzz, gr_z_yyy, gr_z_yyz, gr_z_yzz, gr_z_zzz, grr_y_yz_xxx, grr_y_yz_xxy, grr_y_yz_xxz, grr_y_yz_xyy, grr_y_yz_xyz, grr_y_yz_xzz, grr_y_yz_yyy, grr_y_yz_yyz, grr_y_yz_yzz, grr_y_yz_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yz_xxx[i] = ts_z_xxx[i] * gfe2_0 + gr_z_xxx[i] * gfe_0 + ts_yz_xxx[i] * gfe_0 * gc_y[i] + gr_yz_xxx[i] * gc_y[i];

        grr_y_yz_xxy[i] = ts_z_xxy[i] * gfe2_0 + gr_z_xxy[i] * gfe_0 + ts_yz_xx[i] * gfe2_0 + gr_yz_xx[i] * gfe_0 + ts_yz_xxy[i] * gfe_0 * gc_y[i] + gr_yz_xxy[i] * gc_y[i];

        grr_y_yz_xxz[i] = ts_z_xxz[i] * gfe2_0 + gr_z_xxz[i] * gfe_0 + ts_yz_xxz[i] * gfe_0 * gc_y[i] + gr_yz_xxz[i] * gc_y[i];

        grr_y_yz_xyy[i] = ts_z_xyy[i] * gfe2_0 + gr_z_xyy[i] * gfe_0 + 2.0 * ts_yz_xy[i] * gfe2_0 + 2.0 * gr_yz_xy[i] * gfe_0 + ts_yz_xyy[i] * gfe_0 * gc_y[i] + gr_yz_xyy[i] * gc_y[i];

        grr_y_yz_xyz[i] = ts_z_xyz[i] * gfe2_0 + gr_z_xyz[i] * gfe_0 + ts_yz_xz[i] * gfe2_0 + gr_yz_xz[i] * gfe_0 + ts_yz_xyz[i] * gfe_0 * gc_y[i] + gr_yz_xyz[i] * gc_y[i];

        grr_y_yz_xzz[i] = ts_z_xzz[i] * gfe2_0 + gr_z_xzz[i] * gfe_0 + ts_yz_xzz[i] * gfe_0 * gc_y[i] + gr_yz_xzz[i] * gc_y[i];

        grr_y_yz_yyy[i] = ts_z_yyy[i] * gfe2_0 + gr_z_yyy[i] * gfe_0 + 3.0 * ts_yz_yy[i] * gfe2_0 + 3.0 * gr_yz_yy[i] * gfe_0 + ts_yz_yyy[i] * gfe_0 * gc_y[i] + gr_yz_yyy[i] * gc_y[i];

        grr_y_yz_yyz[i] = ts_z_yyz[i] * gfe2_0 + gr_z_yyz[i] * gfe_0 + 2.0 * ts_yz_yz[i] * gfe2_0 + 2.0 * gr_yz_yz[i] * gfe_0 + ts_yz_yyz[i] * gfe_0 * gc_y[i] + gr_yz_yyz[i] * gc_y[i];

        grr_y_yz_yzz[i] = ts_z_yzz[i] * gfe2_0 + gr_z_yzz[i] * gfe_0 + ts_yz_zz[i] * gfe2_0 + gr_yz_zz[i] * gfe_0 + ts_yz_yzz[i] * gfe_0 * gc_y[i] + gr_yz_yzz[i] * gc_y[i];

        grr_y_yz_zzz[i] = ts_z_zzz[i] * gfe2_0 + gr_z_zzz[i] * gfe_0 + ts_yz_zzz[i] * gfe_0 * gc_y[i] + gr_yz_zzz[i] * gc_y[i];
    }

    // Set up 110-120 components of targeted buffer : DF

    auto grr_y_zz_xxx = pbuffer.data(idx_gr_df + 110);

    auto grr_y_zz_xxy = pbuffer.data(idx_gr_df + 111);

    auto grr_y_zz_xxz = pbuffer.data(idx_gr_df + 112);

    auto grr_y_zz_xyy = pbuffer.data(idx_gr_df + 113);

    auto grr_y_zz_xyz = pbuffer.data(idx_gr_df + 114);

    auto grr_y_zz_xzz = pbuffer.data(idx_gr_df + 115);

    auto grr_y_zz_yyy = pbuffer.data(idx_gr_df + 116);

    auto grr_y_zz_yyz = pbuffer.data(idx_gr_df + 117);

    auto grr_y_zz_yzz = pbuffer.data(idx_gr_df + 118);

    auto grr_y_zz_zzz = pbuffer.data(idx_gr_df + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xx, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xy, gr_zz_xyy, gr_zz_xyz, gr_zz_xz, gr_zz_xzz, gr_zz_yy, gr_zz_yyy, gr_zz_yyz, gr_zz_yz, gr_zz_yzz, gr_zz_zz, gr_zz_zzz, grr_y_zz_xxx, grr_y_zz_xxy, grr_y_zz_xxz, grr_y_zz_xyy, grr_y_zz_xyz, grr_y_zz_xzz, grr_y_zz_yyy, grr_y_zz_yyz, grr_y_zz_yzz, grr_y_zz_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zz_xxx[i] = ts_zz_xxx[i] * gfe_0 * gc_y[i] + gr_zz_xxx[i] * gc_y[i];

        grr_y_zz_xxy[i] = ts_zz_xx[i] * gfe2_0 + gr_zz_xx[i] * gfe_0 + ts_zz_xxy[i] * gfe_0 * gc_y[i] + gr_zz_xxy[i] * gc_y[i];

        grr_y_zz_xxz[i] = ts_zz_xxz[i] * gfe_0 * gc_y[i] + gr_zz_xxz[i] * gc_y[i];

        grr_y_zz_xyy[i] = 2.0 * ts_zz_xy[i] * gfe2_0 + 2.0 * gr_zz_xy[i] * gfe_0 + ts_zz_xyy[i] * gfe_0 * gc_y[i] + gr_zz_xyy[i] * gc_y[i];

        grr_y_zz_xyz[i] = ts_zz_xz[i] * gfe2_0 + gr_zz_xz[i] * gfe_0 + ts_zz_xyz[i] * gfe_0 * gc_y[i] + gr_zz_xyz[i] * gc_y[i];

        grr_y_zz_xzz[i] = ts_zz_xzz[i] * gfe_0 * gc_y[i] + gr_zz_xzz[i] * gc_y[i];

        grr_y_zz_yyy[i] = 3.0 * ts_zz_yy[i] * gfe2_0 + 3.0 * gr_zz_yy[i] * gfe_0 + ts_zz_yyy[i] * gfe_0 * gc_y[i] + gr_zz_yyy[i] * gc_y[i];

        grr_y_zz_yyz[i] = 2.0 * ts_zz_yz[i] * gfe2_0 + 2.0 * gr_zz_yz[i] * gfe_0 + ts_zz_yyz[i] * gfe_0 * gc_y[i] + gr_zz_yyz[i] * gc_y[i];

        grr_y_zz_yzz[i] = ts_zz_zz[i] * gfe2_0 + gr_zz_zz[i] * gfe_0 + ts_zz_yzz[i] * gfe_0 * gc_y[i] + gr_zz_yzz[i] * gc_y[i];

        grr_y_zz_zzz[i] = ts_zz_zzz[i] * gfe_0 * gc_y[i] + gr_zz_zzz[i] * gc_y[i];
    }

    // Set up 120-130 components of targeted buffer : DF

    auto grr_z_xx_xxx = pbuffer.data(idx_gr_df + 120);

    auto grr_z_xx_xxy = pbuffer.data(idx_gr_df + 121);

    auto grr_z_xx_xxz = pbuffer.data(idx_gr_df + 122);

    auto grr_z_xx_xyy = pbuffer.data(idx_gr_df + 123);

    auto grr_z_xx_xyz = pbuffer.data(idx_gr_df + 124);

    auto grr_z_xx_xzz = pbuffer.data(idx_gr_df + 125);

    auto grr_z_xx_yyy = pbuffer.data(idx_gr_df + 126);

    auto grr_z_xx_yyz = pbuffer.data(idx_gr_df + 127);

    auto grr_z_xx_yzz = pbuffer.data(idx_gr_df + 128);

    auto grr_z_xx_zzz = pbuffer.data(idx_gr_df + 129);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xx, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xy, gr_xx_xyy, gr_xx_xyz, gr_xx_xz, gr_xx_xzz, gr_xx_yy, gr_xx_yyy, gr_xx_yyz, gr_xx_yz, gr_xx_yzz, gr_xx_zz, gr_xx_zzz, grr_z_xx_xxx, grr_z_xx_xxy, grr_z_xx_xxz, grr_z_xx_xyy, grr_z_xx_xyz, grr_z_xx_xzz, grr_z_xx_yyy, grr_z_xx_yyz, grr_z_xx_yzz, grr_z_xx_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xx_xxx[i] = ts_xx_xxx[i] * gfe_0 * gc_z[i] + gr_xx_xxx[i] * gc_z[i];

        grr_z_xx_xxy[i] = ts_xx_xxy[i] * gfe_0 * gc_z[i] + gr_xx_xxy[i] * gc_z[i];

        grr_z_xx_xxz[i] = ts_xx_xx[i] * gfe2_0 + gr_xx_xx[i] * gfe_0 + ts_xx_xxz[i] * gfe_0 * gc_z[i] + gr_xx_xxz[i] * gc_z[i];

        grr_z_xx_xyy[i] = ts_xx_xyy[i] * gfe_0 * gc_z[i] + gr_xx_xyy[i] * gc_z[i];

        grr_z_xx_xyz[i] = ts_xx_xy[i] * gfe2_0 + gr_xx_xy[i] * gfe_0 + ts_xx_xyz[i] * gfe_0 * gc_z[i] + gr_xx_xyz[i] * gc_z[i];

        grr_z_xx_xzz[i] = 2.0 * ts_xx_xz[i] * gfe2_0 + 2.0 * gr_xx_xz[i] * gfe_0 + ts_xx_xzz[i] * gfe_0 * gc_z[i] + gr_xx_xzz[i] * gc_z[i];

        grr_z_xx_yyy[i] = ts_xx_yyy[i] * gfe_0 * gc_z[i] + gr_xx_yyy[i] * gc_z[i];

        grr_z_xx_yyz[i] = ts_xx_yy[i] * gfe2_0 + gr_xx_yy[i] * gfe_0 + ts_xx_yyz[i] * gfe_0 * gc_z[i] + gr_xx_yyz[i] * gc_z[i];

        grr_z_xx_yzz[i] = 2.0 * ts_xx_yz[i] * gfe2_0 + 2.0 * gr_xx_yz[i] * gfe_0 + ts_xx_yzz[i] * gfe_0 * gc_z[i] + gr_xx_yzz[i] * gc_z[i];

        grr_z_xx_zzz[i] = 3.0 * ts_xx_zz[i] * gfe2_0 + 3.0 * gr_xx_zz[i] * gfe_0 + ts_xx_zzz[i] * gfe_0 * gc_z[i] + gr_xx_zzz[i] * gc_z[i];
    }

    // Set up 130-140 components of targeted buffer : DF

    auto grr_z_xy_xxx = pbuffer.data(idx_gr_df + 130);

    auto grr_z_xy_xxy = pbuffer.data(idx_gr_df + 131);

    auto grr_z_xy_xxz = pbuffer.data(idx_gr_df + 132);

    auto grr_z_xy_xyy = pbuffer.data(idx_gr_df + 133);

    auto grr_z_xy_xyz = pbuffer.data(idx_gr_df + 134);

    auto grr_z_xy_xzz = pbuffer.data(idx_gr_df + 135);

    auto grr_z_xy_yyy = pbuffer.data(idx_gr_df + 136);

    auto grr_z_xy_yyz = pbuffer.data(idx_gr_df + 137);

    auto grr_z_xy_yzz = pbuffer.data(idx_gr_df + 138);

    auto grr_z_xy_zzz = pbuffer.data(idx_gr_df + 139);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xx, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xy, gr_xy_xyy, gr_xy_xyz, gr_xy_xz, gr_xy_xzz, gr_xy_yy, gr_xy_yyy, gr_xy_yyz, gr_xy_yz, gr_xy_yzz, gr_xy_zz, gr_xy_zzz, grr_z_xy_xxx, grr_z_xy_xxy, grr_z_xy_xxz, grr_z_xy_xyy, grr_z_xy_xyz, grr_z_xy_xzz, grr_z_xy_yyy, grr_z_xy_yyz, grr_z_xy_yzz, grr_z_xy_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xy_xxx[i] = ts_xy_xxx[i] * gfe_0 * gc_z[i] + gr_xy_xxx[i] * gc_z[i];

        grr_z_xy_xxy[i] = ts_xy_xxy[i] * gfe_0 * gc_z[i] + gr_xy_xxy[i] * gc_z[i];

        grr_z_xy_xxz[i] = ts_xy_xx[i] * gfe2_0 + gr_xy_xx[i] * gfe_0 + ts_xy_xxz[i] * gfe_0 * gc_z[i] + gr_xy_xxz[i] * gc_z[i];

        grr_z_xy_xyy[i] = ts_xy_xyy[i] * gfe_0 * gc_z[i] + gr_xy_xyy[i] * gc_z[i];

        grr_z_xy_xyz[i] = ts_xy_xy[i] * gfe2_0 + gr_xy_xy[i] * gfe_0 + ts_xy_xyz[i] * gfe_0 * gc_z[i] + gr_xy_xyz[i] * gc_z[i];

        grr_z_xy_xzz[i] = 2.0 * ts_xy_xz[i] * gfe2_0 + 2.0 * gr_xy_xz[i] * gfe_0 + ts_xy_xzz[i] * gfe_0 * gc_z[i] + gr_xy_xzz[i] * gc_z[i];

        grr_z_xy_yyy[i] = ts_xy_yyy[i] * gfe_0 * gc_z[i] + gr_xy_yyy[i] * gc_z[i];

        grr_z_xy_yyz[i] = ts_xy_yy[i] * gfe2_0 + gr_xy_yy[i] * gfe_0 + ts_xy_yyz[i] * gfe_0 * gc_z[i] + gr_xy_yyz[i] * gc_z[i];

        grr_z_xy_yzz[i] = 2.0 * ts_xy_yz[i] * gfe2_0 + 2.0 * gr_xy_yz[i] * gfe_0 + ts_xy_yzz[i] * gfe_0 * gc_z[i] + gr_xy_yzz[i] * gc_z[i];

        grr_z_xy_zzz[i] = 3.0 * ts_xy_zz[i] * gfe2_0 + 3.0 * gr_xy_zz[i] * gfe_0 + ts_xy_zzz[i] * gfe_0 * gc_z[i] + gr_xy_zzz[i] * gc_z[i];
    }

    // Set up 140-150 components of targeted buffer : DF

    auto grr_z_xz_xxx = pbuffer.data(idx_gr_df + 140);

    auto grr_z_xz_xxy = pbuffer.data(idx_gr_df + 141);

    auto grr_z_xz_xxz = pbuffer.data(idx_gr_df + 142);

    auto grr_z_xz_xyy = pbuffer.data(idx_gr_df + 143);

    auto grr_z_xz_xyz = pbuffer.data(idx_gr_df + 144);

    auto grr_z_xz_xzz = pbuffer.data(idx_gr_df + 145);

    auto grr_z_xz_yyy = pbuffer.data(idx_gr_df + 146);

    auto grr_z_xz_yyz = pbuffer.data(idx_gr_df + 147);

    auto grr_z_xz_yzz = pbuffer.data(idx_gr_df + 148);

    auto grr_z_xz_zzz = pbuffer.data(idx_gr_df + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xyy, gr_x_xyz, gr_x_xzz, gr_x_yyy, gr_x_yyz, gr_x_yzz, gr_x_zzz, gr_xz_xx, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xy, gr_xz_xyy, gr_xz_xyz, gr_xz_xz, gr_xz_xzz, gr_xz_yy, gr_xz_yyy, gr_xz_yyz, gr_xz_yz, gr_xz_yzz, gr_xz_zz, gr_xz_zzz, grr_z_xz_xxx, grr_z_xz_xxy, grr_z_xz_xxz, grr_z_xz_xyy, grr_z_xz_xyz, grr_z_xz_xzz, grr_z_xz_yyy, grr_z_xz_yyz, grr_z_xz_yzz, grr_z_xz_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xz_xxx[i] = ts_x_xxx[i] * gfe2_0 + gr_x_xxx[i] * gfe_0 + ts_xz_xxx[i] * gfe_0 * gc_z[i] + gr_xz_xxx[i] * gc_z[i];

        grr_z_xz_xxy[i] = ts_x_xxy[i] * gfe2_0 + gr_x_xxy[i] * gfe_0 + ts_xz_xxy[i] * gfe_0 * gc_z[i] + gr_xz_xxy[i] * gc_z[i];

        grr_z_xz_xxz[i] = ts_x_xxz[i] * gfe2_0 + gr_x_xxz[i] * gfe_0 + ts_xz_xx[i] * gfe2_0 + gr_xz_xx[i] * gfe_0 + ts_xz_xxz[i] * gfe_0 * gc_z[i] + gr_xz_xxz[i] * gc_z[i];

        grr_z_xz_xyy[i] = ts_x_xyy[i] * gfe2_0 + gr_x_xyy[i] * gfe_0 + ts_xz_xyy[i] * gfe_0 * gc_z[i] + gr_xz_xyy[i] * gc_z[i];

        grr_z_xz_xyz[i] = ts_x_xyz[i] * gfe2_0 + gr_x_xyz[i] * gfe_0 + ts_xz_xy[i] * gfe2_0 + gr_xz_xy[i] * gfe_0 + ts_xz_xyz[i] * gfe_0 * gc_z[i] + gr_xz_xyz[i] * gc_z[i];

        grr_z_xz_xzz[i] = ts_x_xzz[i] * gfe2_0 + gr_x_xzz[i] * gfe_0 + 2.0 * ts_xz_xz[i] * gfe2_0 + 2.0 * gr_xz_xz[i] * gfe_0 + ts_xz_xzz[i] * gfe_0 * gc_z[i] + gr_xz_xzz[i] * gc_z[i];

        grr_z_xz_yyy[i] = ts_x_yyy[i] * gfe2_0 + gr_x_yyy[i] * gfe_0 + ts_xz_yyy[i] * gfe_0 * gc_z[i] + gr_xz_yyy[i] * gc_z[i];

        grr_z_xz_yyz[i] = ts_x_yyz[i] * gfe2_0 + gr_x_yyz[i] * gfe_0 + ts_xz_yy[i] * gfe2_0 + gr_xz_yy[i] * gfe_0 + ts_xz_yyz[i] * gfe_0 * gc_z[i] + gr_xz_yyz[i] * gc_z[i];

        grr_z_xz_yzz[i] = ts_x_yzz[i] * gfe2_0 + gr_x_yzz[i] * gfe_0 + 2.0 * ts_xz_yz[i] * gfe2_0 + 2.0 * gr_xz_yz[i] * gfe_0 + ts_xz_yzz[i] * gfe_0 * gc_z[i] + gr_xz_yzz[i] * gc_z[i];

        grr_z_xz_zzz[i] = ts_x_zzz[i] * gfe2_0 + gr_x_zzz[i] * gfe_0 + 3.0 * ts_xz_zz[i] * gfe2_0 + 3.0 * gr_xz_zz[i] * gfe_0 + ts_xz_zzz[i] * gfe_0 * gc_z[i] + gr_xz_zzz[i] * gc_z[i];
    }

    // Set up 150-160 components of targeted buffer : DF

    auto grr_z_yy_xxx = pbuffer.data(idx_gr_df + 150);

    auto grr_z_yy_xxy = pbuffer.data(idx_gr_df + 151);

    auto grr_z_yy_xxz = pbuffer.data(idx_gr_df + 152);

    auto grr_z_yy_xyy = pbuffer.data(idx_gr_df + 153);

    auto grr_z_yy_xyz = pbuffer.data(idx_gr_df + 154);

    auto grr_z_yy_xzz = pbuffer.data(idx_gr_df + 155);

    auto grr_z_yy_yyy = pbuffer.data(idx_gr_df + 156);

    auto grr_z_yy_yyz = pbuffer.data(idx_gr_df + 157);

    auto grr_z_yy_yzz = pbuffer.data(idx_gr_df + 158);

    auto grr_z_yy_zzz = pbuffer.data(idx_gr_df + 159);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xx, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xy, gr_yy_xyy, gr_yy_xyz, gr_yy_xz, gr_yy_xzz, gr_yy_yy, gr_yy_yyy, gr_yy_yyz, gr_yy_yz, gr_yy_yzz, gr_yy_zz, gr_yy_zzz, grr_z_yy_xxx, grr_z_yy_xxy, grr_z_yy_xxz, grr_z_yy_xyy, grr_z_yy_xyz, grr_z_yy_xzz, grr_z_yy_yyy, grr_z_yy_yyz, grr_z_yy_yzz, grr_z_yy_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yy_xxx[i] = ts_yy_xxx[i] * gfe_0 * gc_z[i] + gr_yy_xxx[i] * gc_z[i];

        grr_z_yy_xxy[i] = ts_yy_xxy[i] * gfe_0 * gc_z[i] + gr_yy_xxy[i] * gc_z[i];

        grr_z_yy_xxz[i] = ts_yy_xx[i] * gfe2_0 + gr_yy_xx[i] * gfe_0 + ts_yy_xxz[i] * gfe_0 * gc_z[i] + gr_yy_xxz[i] * gc_z[i];

        grr_z_yy_xyy[i] = ts_yy_xyy[i] * gfe_0 * gc_z[i] + gr_yy_xyy[i] * gc_z[i];

        grr_z_yy_xyz[i] = ts_yy_xy[i] * gfe2_0 + gr_yy_xy[i] * gfe_0 + ts_yy_xyz[i] * gfe_0 * gc_z[i] + gr_yy_xyz[i] * gc_z[i];

        grr_z_yy_xzz[i] = 2.0 * ts_yy_xz[i] * gfe2_0 + 2.0 * gr_yy_xz[i] * gfe_0 + ts_yy_xzz[i] * gfe_0 * gc_z[i] + gr_yy_xzz[i] * gc_z[i];

        grr_z_yy_yyy[i] = ts_yy_yyy[i] * gfe_0 * gc_z[i] + gr_yy_yyy[i] * gc_z[i];

        grr_z_yy_yyz[i] = ts_yy_yy[i] * gfe2_0 + gr_yy_yy[i] * gfe_0 + ts_yy_yyz[i] * gfe_0 * gc_z[i] + gr_yy_yyz[i] * gc_z[i];

        grr_z_yy_yzz[i] = 2.0 * ts_yy_yz[i] * gfe2_0 + 2.0 * gr_yy_yz[i] * gfe_0 + ts_yy_yzz[i] * gfe_0 * gc_z[i] + gr_yy_yzz[i] * gc_z[i];

        grr_z_yy_zzz[i] = 3.0 * ts_yy_zz[i] * gfe2_0 + 3.0 * gr_yy_zz[i] * gfe_0 + ts_yy_zzz[i] * gfe_0 * gc_z[i] + gr_yy_zzz[i] * gc_z[i];
    }

    // Set up 160-170 components of targeted buffer : DF

    auto grr_z_yz_xxx = pbuffer.data(idx_gr_df + 160);

    auto grr_z_yz_xxy = pbuffer.data(idx_gr_df + 161);

    auto grr_z_yz_xxz = pbuffer.data(idx_gr_df + 162);

    auto grr_z_yz_xyy = pbuffer.data(idx_gr_df + 163);

    auto grr_z_yz_xyz = pbuffer.data(idx_gr_df + 164);

    auto grr_z_yz_xzz = pbuffer.data(idx_gr_df + 165);

    auto grr_z_yz_yyy = pbuffer.data(idx_gr_df + 166);

    auto grr_z_yz_yyz = pbuffer.data(idx_gr_df + 167);

    auto grr_z_yz_yzz = pbuffer.data(idx_gr_df + 168);

    auto grr_z_yz_zzz = pbuffer.data(idx_gr_df + 169);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xyy, gr_y_xyz, gr_y_xzz, gr_y_yyy, gr_y_yyz, gr_y_yzz, gr_y_zzz, gr_yz_xx, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xy, gr_yz_xyy, gr_yz_xyz, gr_yz_xz, gr_yz_xzz, gr_yz_yy, gr_yz_yyy, gr_yz_yyz, gr_yz_yz, gr_yz_yzz, gr_yz_zz, gr_yz_zzz, grr_z_yz_xxx, grr_z_yz_xxy, grr_z_yz_xxz, grr_z_yz_xyy, grr_z_yz_xyz, grr_z_yz_xzz, grr_z_yz_yyy, grr_z_yz_yyz, grr_z_yz_yzz, grr_z_yz_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yz_xxx[i] = ts_y_xxx[i] * gfe2_0 + gr_y_xxx[i] * gfe_0 + ts_yz_xxx[i] * gfe_0 * gc_z[i] + gr_yz_xxx[i] * gc_z[i];

        grr_z_yz_xxy[i] = ts_y_xxy[i] * gfe2_0 + gr_y_xxy[i] * gfe_0 + ts_yz_xxy[i] * gfe_0 * gc_z[i] + gr_yz_xxy[i] * gc_z[i];

        grr_z_yz_xxz[i] = ts_y_xxz[i] * gfe2_0 + gr_y_xxz[i] * gfe_0 + ts_yz_xx[i] * gfe2_0 + gr_yz_xx[i] * gfe_0 + ts_yz_xxz[i] * gfe_0 * gc_z[i] + gr_yz_xxz[i] * gc_z[i];

        grr_z_yz_xyy[i] = ts_y_xyy[i] * gfe2_0 + gr_y_xyy[i] * gfe_0 + ts_yz_xyy[i] * gfe_0 * gc_z[i] + gr_yz_xyy[i] * gc_z[i];

        grr_z_yz_xyz[i] = ts_y_xyz[i] * gfe2_0 + gr_y_xyz[i] * gfe_0 + ts_yz_xy[i] * gfe2_0 + gr_yz_xy[i] * gfe_0 + ts_yz_xyz[i] * gfe_0 * gc_z[i] + gr_yz_xyz[i] * gc_z[i];

        grr_z_yz_xzz[i] = ts_y_xzz[i] * gfe2_0 + gr_y_xzz[i] * gfe_0 + 2.0 * ts_yz_xz[i] * gfe2_0 + 2.0 * gr_yz_xz[i] * gfe_0 + ts_yz_xzz[i] * gfe_0 * gc_z[i] + gr_yz_xzz[i] * gc_z[i];

        grr_z_yz_yyy[i] = ts_y_yyy[i] * gfe2_0 + gr_y_yyy[i] * gfe_0 + ts_yz_yyy[i] * gfe_0 * gc_z[i] + gr_yz_yyy[i] * gc_z[i];

        grr_z_yz_yyz[i] = ts_y_yyz[i] * gfe2_0 + gr_y_yyz[i] * gfe_0 + ts_yz_yy[i] * gfe2_0 + gr_yz_yy[i] * gfe_0 + ts_yz_yyz[i] * gfe_0 * gc_z[i] + gr_yz_yyz[i] * gc_z[i];

        grr_z_yz_yzz[i] = ts_y_yzz[i] * gfe2_0 + gr_y_yzz[i] * gfe_0 + 2.0 * ts_yz_yz[i] * gfe2_0 + 2.0 * gr_yz_yz[i] * gfe_0 + ts_yz_yzz[i] * gfe_0 * gc_z[i] + gr_yz_yzz[i] * gc_z[i];

        grr_z_yz_zzz[i] = ts_y_zzz[i] * gfe2_0 + gr_y_zzz[i] * gfe_0 + 3.0 * ts_yz_zz[i] * gfe2_0 + 3.0 * gr_yz_zz[i] * gfe_0 + ts_yz_zzz[i] * gfe_0 * gc_z[i] + gr_yz_zzz[i] * gc_z[i];
    }

    // Set up 170-180 components of targeted buffer : DF

    auto grr_z_zz_xxx = pbuffer.data(idx_gr_df + 170);

    auto grr_z_zz_xxy = pbuffer.data(idx_gr_df + 171);

    auto grr_z_zz_xxz = pbuffer.data(idx_gr_df + 172);

    auto grr_z_zz_xyy = pbuffer.data(idx_gr_df + 173);

    auto grr_z_zz_xyz = pbuffer.data(idx_gr_df + 174);

    auto grr_z_zz_xzz = pbuffer.data(idx_gr_df + 175);

    auto grr_z_zz_yyy = pbuffer.data(idx_gr_df + 176);

    auto grr_z_zz_yyz = pbuffer.data(idx_gr_df + 177);

    auto grr_z_zz_yzz = pbuffer.data(idx_gr_df + 178);

    auto grr_z_zz_zzz = pbuffer.data(idx_gr_df + 179);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xyy, gr_z_xyz, gr_z_xzz, gr_z_yyy, gr_z_yyz, gr_z_yzz, gr_z_zzz, gr_zz_xx, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xy, gr_zz_xyy, gr_zz_xyz, gr_zz_xz, gr_zz_xzz, gr_zz_yy, gr_zz_yyy, gr_zz_yyz, gr_zz_yz, gr_zz_yzz, gr_zz_zz, gr_zz_zzz, grr_z_zz_xxx, grr_z_zz_xxy, grr_z_zz_xxz, grr_z_zz_xyy, grr_z_zz_xyz, grr_z_zz_xzz, grr_z_zz_yyy, grr_z_zz_yyz, grr_z_zz_yzz, grr_z_zz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zz_xxx[i] = 2.0 * ts_z_xxx[i] * gfe2_0 + 2.0 * gr_z_xxx[i] * gfe_0 + ts_zz_xxx[i] * gfe_0 * gc_z[i] + gr_zz_xxx[i] * gc_z[i];

        grr_z_zz_xxy[i] = 2.0 * ts_z_xxy[i] * gfe2_0 + 2.0 * gr_z_xxy[i] * gfe_0 + ts_zz_xxy[i] * gfe_0 * gc_z[i] + gr_zz_xxy[i] * gc_z[i];

        grr_z_zz_xxz[i] = 2.0 * ts_z_xxz[i] * gfe2_0 + 2.0 * gr_z_xxz[i] * gfe_0 + ts_zz_xx[i] * gfe2_0 + gr_zz_xx[i] * gfe_0 + ts_zz_xxz[i] * gfe_0 * gc_z[i] + gr_zz_xxz[i] * gc_z[i];

        grr_z_zz_xyy[i] = 2.0 * ts_z_xyy[i] * gfe2_0 + 2.0 * gr_z_xyy[i] * gfe_0 + ts_zz_xyy[i] * gfe_0 * gc_z[i] + gr_zz_xyy[i] * gc_z[i];

        grr_z_zz_xyz[i] = 2.0 * ts_z_xyz[i] * gfe2_0 + 2.0 * gr_z_xyz[i] * gfe_0 + ts_zz_xy[i] * gfe2_0 + gr_zz_xy[i] * gfe_0 + ts_zz_xyz[i] * gfe_0 * gc_z[i] + gr_zz_xyz[i] * gc_z[i];

        grr_z_zz_xzz[i] = 2.0 * ts_z_xzz[i] * gfe2_0 + 2.0 * gr_z_xzz[i] * gfe_0 + 2.0 * ts_zz_xz[i] * gfe2_0 + 2.0 * gr_zz_xz[i] * gfe_0 + ts_zz_xzz[i] * gfe_0 * gc_z[i] + gr_zz_xzz[i] * gc_z[i];

        grr_z_zz_yyy[i] = 2.0 * ts_z_yyy[i] * gfe2_0 + 2.0 * gr_z_yyy[i] * gfe_0 + ts_zz_yyy[i] * gfe_0 * gc_z[i] + gr_zz_yyy[i] * gc_z[i];

        grr_z_zz_yyz[i] = 2.0 * ts_z_yyz[i] * gfe2_0 + 2.0 * gr_z_yyz[i] * gfe_0 + ts_zz_yy[i] * gfe2_0 + gr_zz_yy[i] * gfe_0 + ts_zz_yyz[i] * gfe_0 * gc_z[i] + gr_zz_yyz[i] * gc_z[i];

        grr_z_zz_yzz[i] = 2.0 * ts_z_yzz[i] * gfe2_0 + 2.0 * gr_z_yzz[i] * gfe_0 + 2.0 * ts_zz_yz[i] * gfe2_0 + 2.0 * gr_zz_yz[i] * gfe_0 + ts_zz_yzz[i] * gfe_0 * gc_z[i] + gr_zz_yzz[i] * gc_z[i];

        grr_z_zz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe2_0 + 2.0 * gr_z_zzz[i] * gfe_0 + 3.0 * ts_zz_zz[i] * gfe2_0 + 3.0 * gr_zz_zz[i] * gfe_0 + ts_zz_zzz[i] * gfe_0 * gc_z[i] + gr_zz_zzz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

