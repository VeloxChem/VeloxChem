#include "ThreeCenterOverlapGradientPrimRecFF.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_ff(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_ff,
                              const size_t idx_df,
                              const size_t idx_fd,
                              const size_t idx_ff,
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

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_fd);

    auto ts_xxx_xy = pbuffer.data(idx_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_fd + 10);

    auto ts_xxy_zz = pbuffer.data(idx_fd + 11);

    auto ts_xxz_xx = pbuffer.data(idx_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_fd + 14);

    auto ts_xxz_yy = pbuffer.data(idx_fd + 15);

    auto ts_xxz_yz = pbuffer.data(idx_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_fd + 17);

    auto ts_xyy_xx = pbuffer.data(idx_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_fd + 19);

    auto ts_xyy_xz = pbuffer.data(idx_fd + 20);

    auto ts_xyy_yy = pbuffer.data(idx_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_fd + 23);

    auto ts_xyz_xx = pbuffer.data(idx_fd + 24);

    auto ts_xyz_xy = pbuffer.data(idx_fd + 25);

    auto ts_xyz_xz = pbuffer.data(idx_fd + 26);

    auto ts_xyz_yy = pbuffer.data(idx_fd + 27);

    auto ts_xyz_yz = pbuffer.data(idx_fd + 28);

    auto ts_xyz_zz = pbuffer.data(idx_fd + 29);

    auto ts_xzz_xx = pbuffer.data(idx_fd + 30);

    auto ts_xzz_xy = pbuffer.data(idx_fd + 31);

    auto ts_xzz_xz = pbuffer.data(idx_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_fd + 41);

    auto ts_yyz_xx = pbuffer.data(idx_fd + 42);

    auto ts_yyz_xy = pbuffer.data(idx_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_fd + 47);

    auto ts_yzz_xx = pbuffer.data(idx_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ff + 9);

    auto ts_xxy_xxx = pbuffer.data(idx_ff + 10);

    auto ts_xxy_xxy = pbuffer.data(idx_ff + 11);

    auto ts_xxy_xxz = pbuffer.data(idx_ff + 12);

    auto ts_xxy_xyy = pbuffer.data(idx_ff + 13);

    auto ts_xxy_xyz = pbuffer.data(idx_ff + 14);

    auto ts_xxy_xzz = pbuffer.data(idx_ff + 15);

    auto ts_xxy_yyy = pbuffer.data(idx_ff + 16);

    auto ts_xxy_yyz = pbuffer.data(idx_ff + 17);

    auto ts_xxy_yzz = pbuffer.data(idx_ff + 18);

    auto ts_xxy_zzz = pbuffer.data(idx_ff + 19);

    auto ts_xxz_xxx = pbuffer.data(idx_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ff + 21);

    auto ts_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto ts_xxz_xyy = pbuffer.data(idx_ff + 23);

    auto ts_xxz_xyz = pbuffer.data(idx_ff + 24);

    auto ts_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto ts_xxz_yyy = pbuffer.data(idx_ff + 26);

    auto ts_xxz_yyz = pbuffer.data(idx_ff + 27);

    auto ts_xxz_yzz = pbuffer.data(idx_ff + 28);

    auto ts_xxz_zzz = pbuffer.data(idx_ff + 29);

    auto ts_xyy_xxx = pbuffer.data(idx_ff + 30);

    auto ts_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto ts_xyy_xxz = pbuffer.data(idx_ff + 32);

    auto ts_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto ts_xyy_xzz = pbuffer.data(idx_ff + 35);

    auto ts_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ff + 39);

    auto ts_xyz_xxx = pbuffer.data(idx_ff + 40);

    auto ts_xyz_xxy = pbuffer.data(idx_ff + 41);

    auto ts_xyz_xxz = pbuffer.data(idx_ff + 42);

    auto ts_xyz_xyy = pbuffer.data(idx_ff + 43);

    auto ts_xyz_xyz = pbuffer.data(idx_ff + 44);

    auto ts_xyz_xzz = pbuffer.data(idx_ff + 45);

    auto ts_xyz_yyy = pbuffer.data(idx_ff + 46);

    auto ts_xyz_yyz = pbuffer.data(idx_ff + 47);

    auto ts_xyz_yzz = pbuffer.data(idx_ff + 48);

    auto ts_xyz_zzz = pbuffer.data(idx_ff + 49);

    auto ts_xzz_xxx = pbuffer.data(idx_ff + 50);

    auto ts_xzz_xxy = pbuffer.data(idx_ff + 51);

    auto ts_xzz_xxz = pbuffer.data(idx_ff + 52);

    auto ts_xzz_xyy = pbuffer.data(idx_ff + 53);

    auto ts_xzz_xyz = pbuffer.data(idx_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ff + 55);

    auto ts_xzz_yyy = pbuffer.data(idx_ff + 56);

    auto ts_xzz_yyz = pbuffer.data(idx_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ff + 69);

    auto ts_yyz_xxx = pbuffer.data(idx_ff + 70);

    auto ts_yyz_xxy = pbuffer.data(idx_ff + 71);

    auto ts_yyz_xxz = pbuffer.data(idx_ff + 72);

    auto ts_yyz_xyy = pbuffer.data(idx_ff + 73);

    auto ts_yyz_xyz = pbuffer.data(idx_ff + 74);

    auto ts_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto ts_yyz_yyy = pbuffer.data(idx_ff + 76);

    auto ts_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto ts_yzz_xxx = pbuffer.data(idx_ff + 80);

    auto ts_yzz_xxy = pbuffer.data(idx_ff + 81);

    auto ts_yzz_xxz = pbuffer.data(idx_ff + 82);

    auto ts_yzz_xyy = pbuffer.data(idx_ff + 83);

    auto ts_yzz_xyz = pbuffer.data(idx_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ff + 99);

    // Set up 0-10 components of targeted buffer : FF

    auto gs_x_xxx_xxx = pbuffer.data(idx_g_ff);

    auto gs_x_xxx_xxy = pbuffer.data(idx_g_ff + 1);

    auto gs_x_xxx_xxz = pbuffer.data(idx_g_ff + 2);

    auto gs_x_xxx_xyy = pbuffer.data(idx_g_ff + 3);

    auto gs_x_xxx_xyz = pbuffer.data(idx_g_ff + 4);

    auto gs_x_xxx_xzz = pbuffer.data(idx_g_ff + 5);

    auto gs_x_xxx_yyy = pbuffer.data(idx_g_ff + 6);

    auto gs_x_xxx_yyz = pbuffer.data(idx_g_ff + 7);

    auto gs_x_xxx_yzz = pbuffer.data(idx_g_ff + 8);

    auto gs_x_xxx_zzz = pbuffer.data(idx_g_ff + 9);

    #pragma omp simd aligned(gc_x, gs_x_xxx_xxx, gs_x_xxx_xxy, gs_x_xxx_xxz, gs_x_xxx_xyy, gs_x_xxx_xyz, gs_x_xxx_xzz, gs_x_xxx_yyy, gs_x_xxx_yyz, gs_x_xxx_yzz, gs_x_xxx_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxx_xxx[i] = 6.0 * ts_xx_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxy[i] = 6.0 * ts_xx_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxz[i] = 6.0 * ts_xx_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyy[i] = 6.0 * ts_xx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyz[i] = 6.0 * ts_xx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xzz[i] = 6.0 * ts_xx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyy[i] = 6.0 * ts_xx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyz[i] = 6.0 * ts_xx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yzz[i] = 6.0 * ts_xx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_zzz[i] = 6.0 * ts_xx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 10-20 components of targeted buffer : FF

    auto gs_x_xxy_xxx = pbuffer.data(idx_g_ff + 10);

    auto gs_x_xxy_xxy = pbuffer.data(idx_g_ff + 11);

    auto gs_x_xxy_xxz = pbuffer.data(idx_g_ff + 12);

    auto gs_x_xxy_xyy = pbuffer.data(idx_g_ff + 13);

    auto gs_x_xxy_xyz = pbuffer.data(idx_g_ff + 14);

    auto gs_x_xxy_xzz = pbuffer.data(idx_g_ff + 15);

    auto gs_x_xxy_yyy = pbuffer.data(idx_g_ff + 16);

    auto gs_x_xxy_yyz = pbuffer.data(idx_g_ff + 17);

    auto gs_x_xxy_yzz = pbuffer.data(idx_g_ff + 18);

    auto gs_x_xxy_zzz = pbuffer.data(idx_g_ff + 19);

    #pragma omp simd aligned(gc_x, gs_x_xxy_xxx, gs_x_xxy_xxy, gs_x_xxy_xxz, gs_x_xxy_xyy, gs_x_xxy_xyz, gs_x_xxy_xzz, gs_x_xxy_yyy, gs_x_xxy_yyz, gs_x_xxy_yzz, gs_x_xxy_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxy_xxx[i] = 4.0 * ts_xy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxy[i] = 4.0 * ts_xy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxz[i] = 4.0 * ts_xy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyy[i] = 4.0 * ts_xy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyz[i] = 4.0 * ts_xy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xzz[i] = 4.0 * ts_xy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyy[i] = 4.0 * ts_xy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyz[i] = 4.0 * ts_xy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yzz[i] = 4.0 * ts_xy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_zzz[i] = 4.0 * ts_xy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 20-30 components of targeted buffer : FF

    auto gs_x_xxz_xxx = pbuffer.data(idx_g_ff + 20);

    auto gs_x_xxz_xxy = pbuffer.data(idx_g_ff + 21);

    auto gs_x_xxz_xxz = pbuffer.data(idx_g_ff + 22);

    auto gs_x_xxz_xyy = pbuffer.data(idx_g_ff + 23);

    auto gs_x_xxz_xyz = pbuffer.data(idx_g_ff + 24);

    auto gs_x_xxz_xzz = pbuffer.data(idx_g_ff + 25);

    auto gs_x_xxz_yyy = pbuffer.data(idx_g_ff + 26);

    auto gs_x_xxz_yyz = pbuffer.data(idx_g_ff + 27);

    auto gs_x_xxz_yzz = pbuffer.data(idx_g_ff + 28);

    auto gs_x_xxz_zzz = pbuffer.data(idx_g_ff + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxz_xxx, gs_x_xxz_xxy, gs_x_xxz_xxz, gs_x_xxz_xyy, gs_x_xxz_xyz, gs_x_xxz_xzz, gs_x_xxz_yyy, gs_x_xxz_yyz, gs_x_xxz_yzz, gs_x_xxz_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxz_xxx[i] = 4.0 * ts_xz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxy[i] = 4.0 * ts_xz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxz[i] = 4.0 * ts_xz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyy[i] = 4.0 * ts_xz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyz[i] = 4.0 * ts_xz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xzz[i] = 4.0 * ts_xz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyy[i] = 4.0 * ts_xz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyz[i] = 4.0 * ts_xz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yzz[i] = 4.0 * ts_xz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_zzz[i] = 4.0 * ts_xz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-40 components of targeted buffer : FF

    auto gs_x_xyy_xxx = pbuffer.data(idx_g_ff + 30);

    auto gs_x_xyy_xxy = pbuffer.data(idx_g_ff + 31);

    auto gs_x_xyy_xxz = pbuffer.data(idx_g_ff + 32);

    auto gs_x_xyy_xyy = pbuffer.data(idx_g_ff + 33);

    auto gs_x_xyy_xyz = pbuffer.data(idx_g_ff + 34);

    auto gs_x_xyy_xzz = pbuffer.data(idx_g_ff + 35);

    auto gs_x_xyy_yyy = pbuffer.data(idx_g_ff + 36);

    auto gs_x_xyy_yyz = pbuffer.data(idx_g_ff + 37);

    auto gs_x_xyy_yzz = pbuffer.data(idx_g_ff + 38);

    auto gs_x_xyy_zzz = pbuffer.data(idx_g_ff + 39);

    #pragma omp simd aligned(gc_x, gs_x_xyy_xxx, gs_x_xyy_xxy, gs_x_xyy_xxz, gs_x_xyy_xyy, gs_x_xyy_xyz, gs_x_xyy_xzz, gs_x_xyy_yyy, gs_x_xyy_yyz, gs_x_xyy_yzz, gs_x_xyy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyy_xxx[i] = 2.0 * ts_yy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxy[i] = 2.0 * ts_yy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxz[i] = 2.0 * ts_yy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyy[i] = 2.0 * ts_yy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyz[i] = 2.0 * ts_yy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xzz[i] = 2.0 * ts_yy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyy[i] = 2.0 * ts_yy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyz[i] = 2.0 * ts_yy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yzz[i] = 2.0 * ts_yy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_zzz[i] = 2.0 * ts_yy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 40-50 components of targeted buffer : FF

    auto gs_x_xyz_xxx = pbuffer.data(idx_g_ff + 40);

    auto gs_x_xyz_xxy = pbuffer.data(idx_g_ff + 41);

    auto gs_x_xyz_xxz = pbuffer.data(idx_g_ff + 42);

    auto gs_x_xyz_xyy = pbuffer.data(idx_g_ff + 43);

    auto gs_x_xyz_xyz = pbuffer.data(idx_g_ff + 44);

    auto gs_x_xyz_xzz = pbuffer.data(idx_g_ff + 45);

    auto gs_x_xyz_yyy = pbuffer.data(idx_g_ff + 46);

    auto gs_x_xyz_yyz = pbuffer.data(idx_g_ff + 47);

    auto gs_x_xyz_yzz = pbuffer.data(idx_g_ff + 48);

    auto gs_x_xyz_zzz = pbuffer.data(idx_g_ff + 49);

    #pragma omp simd aligned(gc_x, gs_x_xyz_xxx, gs_x_xyz_xxy, gs_x_xyz_xxz, gs_x_xyz_xyy, gs_x_xyz_xyz, gs_x_xyz_xzz, gs_x_xyz_yyy, gs_x_xyz_yyz, gs_x_xyz_yzz, gs_x_xyz_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyz_xxx[i] = 2.0 * ts_yz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxy[i] = 2.0 * ts_yz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxz[i] = 2.0 * ts_yz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyy[i] = 2.0 * ts_yz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyz[i] = 2.0 * ts_yz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xzz[i] = 2.0 * ts_yz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyy[i] = 2.0 * ts_yz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyz[i] = 2.0 * ts_yz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yzz[i] = 2.0 * ts_yz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_zzz[i] = 2.0 * ts_yz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 50-60 components of targeted buffer : FF

    auto gs_x_xzz_xxx = pbuffer.data(idx_g_ff + 50);

    auto gs_x_xzz_xxy = pbuffer.data(idx_g_ff + 51);

    auto gs_x_xzz_xxz = pbuffer.data(idx_g_ff + 52);

    auto gs_x_xzz_xyy = pbuffer.data(idx_g_ff + 53);

    auto gs_x_xzz_xyz = pbuffer.data(idx_g_ff + 54);

    auto gs_x_xzz_xzz = pbuffer.data(idx_g_ff + 55);

    auto gs_x_xzz_yyy = pbuffer.data(idx_g_ff + 56);

    auto gs_x_xzz_yyz = pbuffer.data(idx_g_ff + 57);

    auto gs_x_xzz_yzz = pbuffer.data(idx_g_ff + 58);

    auto gs_x_xzz_zzz = pbuffer.data(idx_g_ff + 59);

    #pragma omp simd aligned(gc_x, gs_x_xzz_xxx, gs_x_xzz_xxy, gs_x_xzz_xxz, gs_x_xzz_xyy, gs_x_xzz_xyz, gs_x_xzz_xzz, gs_x_xzz_yyy, gs_x_xzz_yyz, gs_x_xzz_yzz, gs_x_xzz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxy[i] = 2.0 * ts_zz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyy[i] = 2.0 * ts_zz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyz[i] = 2.0 * ts_zz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-70 components of targeted buffer : FF

    auto gs_x_yyy_xxx = pbuffer.data(idx_g_ff + 60);

    auto gs_x_yyy_xxy = pbuffer.data(idx_g_ff + 61);

    auto gs_x_yyy_xxz = pbuffer.data(idx_g_ff + 62);

    auto gs_x_yyy_xyy = pbuffer.data(idx_g_ff + 63);

    auto gs_x_yyy_xyz = pbuffer.data(idx_g_ff + 64);

    auto gs_x_yyy_xzz = pbuffer.data(idx_g_ff + 65);

    auto gs_x_yyy_yyy = pbuffer.data(idx_g_ff + 66);

    auto gs_x_yyy_yyz = pbuffer.data(idx_g_ff + 67);

    auto gs_x_yyy_yzz = pbuffer.data(idx_g_ff + 68);

    auto gs_x_yyy_zzz = pbuffer.data(idx_g_ff + 69);

    #pragma omp simd aligned(gc_x, gs_x_yyy_xxx, gs_x_yyy_xxy, gs_x_yyy_xxz, gs_x_yyy_xyy, gs_x_yyy_xyz, gs_x_yyy_xzz, gs_x_yyy_yyy, gs_x_yyy_yyz, gs_x_yyy_yzz, gs_x_yyy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyy_xxx[i] = 6.0 * ts_yyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxy[i] = 4.0 * ts_yyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxz[i] = 4.0 * ts_yyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyy[i] = 2.0 * ts_yyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyz[i] = 2.0 * ts_yyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xzz[i] = 2.0 * ts_yyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyy[i] = 2.0 * ts_yyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyz[i] = 2.0 * ts_yyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yzz[i] = 2.0 * ts_yyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_zzz[i] = 2.0 * ts_yyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 70-80 components of targeted buffer : FF

    auto gs_x_yyz_xxx = pbuffer.data(idx_g_ff + 70);

    auto gs_x_yyz_xxy = pbuffer.data(idx_g_ff + 71);

    auto gs_x_yyz_xxz = pbuffer.data(idx_g_ff + 72);

    auto gs_x_yyz_xyy = pbuffer.data(idx_g_ff + 73);

    auto gs_x_yyz_xyz = pbuffer.data(idx_g_ff + 74);

    auto gs_x_yyz_xzz = pbuffer.data(idx_g_ff + 75);

    auto gs_x_yyz_yyy = pbuffer.data(idx_g_ff + 76);

    auto gs_x_yyz_yyz = pbuffer.data(idx_g_ff + 77);

    auto gs_x_yyz_yzz = pbuffer.data(idx_g_ff + 78);

    auto gs_x_yyz_zzz = pbuffer.data(idx_g_ff + 79);

    #pragma omp simd aligned(gc_x, gs_x_yyz_xxx, gs_x_yyz_xxy, gs_x_yyz_xxz, gs_x_yyz_xyy, gs_x_yyz_xyz, gs_x_yyz_xzz, gs_x_yyz_yyy, gs_x_yyz_yyz, gs_x_yyz_yzz, gs_x_yyz_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyz_xxx[i] = 6.0 * ts_yyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxy[i] = 4.0 * ts_yyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxz[i] = 4.0 * ts_yyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyy[i] = 2.0 * ts_yyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyz[i] = 2.0 * ts_yyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xzz[i] = 2.0 * ts_yyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyy[i] = 2.0 * ts_yyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyz[i] = 2.0 * ts_yyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yzz[i] = 2.0 * ts_yyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_zzz[i] = 2.0 * ts_yyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 80-90 components of targeted buffer : FF

    auto gs_x_yzz_xxx = pbuffer.data(idx_g_ff + 80);

    auto gs_x_yzz_xxy = pbuffer.data(idx_g_ff + 81);

    auto gs_x_yzz_xxz = pbuffer.data(idx_g_ff + 82);

    auto gs_x_yzz_xyy = pbuffer.data(idx_g_ff + 83);

    auto gs_x_yzz_xyz = pbuffer.data(idx_g_ff + 84);

    auto gs_x_yzz_xzz = pbuffer.data(idx_g_ff + 85);

    auto gs_x_yzz_yyy = pbuffer.data(idx_g_ff + 86);

    auto gs_x_yzz_yyz = pbuffer.data(idx_g_ff + 87);

    auto gs_x_yzz_yzz = pbuffer.data(idx_g_ff + 88);

    auto gs_x_yzz_zzz = pbuffer.data(idx_g_ff + 89);

    #pragma omp simd aligned(gc_x, gs_x_yzz_xxx, gs_x_yzz_xxy, gs_x_yzz_xxz, gs_x_yzz_xyy, gs_x_yzz_xyz, gs_x_yzz_xzz, gs_x_yzz_yyy, gs_x_yzz_yyz, gs_x_yzz_yzz, gs_x_yzz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzz_xxx[i] = 6.0 * ts_yzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxy[i] = 4.0 * ts_yzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxz[i] = 4.0 * ts_yzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyy[i] = 2.0 * ts_yzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyz[i] = 2.0 * ts_yzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xzz[i] = 2.0 * ts_yzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyy[i] = 2.0 * ts_yzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyz[i] = 2.0 * ts_yzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yzz[i] = 2.0 * ts_yzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_zzz[i] = 2.0 * ts_yzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-100 components of targeted buffer : FF

    auto gs_x_zzz_xxx = pbuffer.data(idx_g_ff + 90);

    auto gs_x_zzz_xxy = pbuffer.data(idx_g_ff + 91);

    auto gs_x_zzz_xxz = pbuffer.data(idx_g_ff + 92);

    auto gs_x_zzz_xyy = pbuffer.data(idx_g_ff + 93);

    auto gs_x_zzz_xyz = pbuffer.data(idx_g_ff + 94);

    auto gs_x_zzz_xzz = pbuffer.data(idx_g_ff + 95);

    auto gs_x_zzz_yyy = pbuffer.data(idx_g_ff + 96);

    auto gs_x_zzz_yyz = pbuffer.data(idx_g_ff + 97);

    auto gs_x_zzz_yzz = pbuffer.data(idx_g_ff + 98);

    auto gs_x_zzz_zzz = pbuffer.data(idx_g_ff + 99);

    #pragma omp simd aligned(gc_x, gs_x_zzz_xxx, gs_x_zzz_xxy, gs_x_zzz_xxz, gs_x_zzz_xyy, gs_x_zzz_xyz, gs_x_zzz_xzz, gs_x_zzz_yyy, gs_x_zzz_yyz, gs_x_zzz_yzz, gs_x_zzz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzz_xxx[i] = 6.0 * ts_zzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxy[i] = 4.0 * ts_zzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxz[i] = 4.0 * ts_zzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyy[i] = 2.0 * ts_zzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyz[i] = 2.0 * ts_zzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xzz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyy[i] = 2.0 * ts_zzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyz[i] = 2.0 * ts_zzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yzz[i] = 2.0 * ts_zzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 100-110 components of targeted buffer : FF

    auto gs_y_xxx_xxx = pbuffer.data(idx_g_ff + 100);

    auto gs_y_xxx_xxy = pbuffer.data(idx_g_ff + 101);

    auto gs_y_xxx_xxz = pbuffer.data(idx_g_ff + 102);

    auto gs_y_xxx_xyy = pbuffer.data(idx_g_ff + 103);

    auto gs_y_xxx_xyz = pbuffer.data(idx_g_ff + 104);

    auto gs_y_xxx_xzz = pbuffer.data(idx_g_ff + 105);

    auto gs_y_xxx_yyy = pbuffer.data(idx_g_ff + 106);

    auto gs_y_xxx_yyz = pbuffer.data(idx_g_ff + 107);

    auto gs_y_xxx_yzz = pbuffer.data(idx_g_ff + 108);

    auto gs_y_xxx_zzz = pbuffer.data(idx_g_ff + 109);

    #pragma omp simd aligned(gc_y, gs_y_xxx_xxx, gs_y_xxx_xxy, gs_y_xxx_xxz, gs_y_xxx_xyy, gs_y_xxx_xyz, gs_y_xxx_xzz, gs_y_xxx_yyy, gs_y_xxx_yyz, gs_y_xxx_yzz, gs_y_xxx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxx_xxx[i] = 2.0 * ts_xxx_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxy[i] = 2.0 * ts_xxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxz[i] = 2.0 * ts_xxx_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyy[i] = 4.0 * ts_xxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyz[i] = 2.0 * ts_xxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xzz[i] = 2.0 * ts_xxx_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyy[i] = 6.0 * ts_xxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyz[i] = 4.0 * ts_xxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yzz[i] = 2.0 * ts_xxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_zzz[i] = 2.0 * ts_xxx_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 110-120 components of targeted buffer : FF

    auto gs_y_xxy_xxx = pbuffer.data(idx_g_ff + 110);

    auto gs_y_xxy_xxy = pbuffer.data(idx_g_ff + 111);

    auto gs_y_xxy_xxz = pbuffer.data(idx_g_ff + 112);

    auto gs_y_xxy_xyy = pbuffer.data(idx_g_ff + 113);

    auto gs_y_xxy_xyz = pbuffer.data(idx_g_ff + 114);

    auto gs_y_xxy_xzz = pbuffer.data(idx_g_ff + 115);

    auto gs_y_xxy_yyy = pbuffer.data(idx_g_ff + 116);

    auto gs_y_xxy_yyz = pbuffer.data(idx_g_ff + 117);

    auto gs_y_xxy_yzz = pbuffer.data(idx_g_ff + 118);

    auto gs_y_xxy_zzz = pbuffer.data(idx_g_ff + 119);

    #pragma omp simd aligned(gc_y, gs_y_xxy_xxx, gs_y_xxy_xxy, gs_y_xxy_xxz, gs_y_xxy_xyy, gs_y_xxy_xyz, gs_y_xxy_xzz, gs_y_xxy_yyy, gs_y_xxy_yyz, gs_y_xxy_yzz, gs_y_xxy_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxy_xxx[i] = 2.0 * ts_xx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxy[i] = 2.0 * ts_xx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxz[i] = 2.0 * ts_xx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyy[i] = 2.0 * ts_xx_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyz[i] = 2.0 * ts_xx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xzz[i] = 2.0 * ts_xx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyy[i] = 2.0 * ts_xx_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyz[i] = 2.0 * ts_xx_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yzz[i] = 2.0 * ts_xx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_zzz[i] = 2.0 * ts_xx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 120-130 components of targeted buffer : FF

    auto gs_y_xxz_xxx = pbuffer.data(idx_g_ff + 120);

    auto gs_y_xxz_xxy = pbuffer.data(idx_g_ff + 121);

    auto gs_y_xxz_xxz = pbuffer.data(idx_g_ff + 122);

    auto gs_y_xxz_xyy = pbuffer.data(idx_g_ff + 123);

    auto gs_y_xxz_xyz = pbuffer.data(idx_g_ff + 124);

    auto gs_y_xxz_xzz = pbuffer.data(idx_g_ff + 125);

    auto gs_y_xxz_yyy = pbuffer.data(idx_g_ff + 126);

    auto gs_y_xxz_yyz = pbuffer.data(idx_g_ff + 127);

    auto gs_y_xxz_yzz = pbuffer.data(idx_g_ff + 128);

    auto gs_y_xxz_zzz = pbuffer.data(idx_g_ff + 129);

    #pragma omp simd aligned(gc_y, gs_y_xxz_xxx, gs_y_xxz_xxy, gs_y_xxz_xxz, gs_y_xxz_xyy, gs_y_xxz_xyz, gs_y_xxz_xzz, gs_y_xxz_yyy, gs_y_xxz_yyz, gs_y_xxz_yzz, gs_y_xxz_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxz_xxx[i] = 2.0 * ts_xxz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxy[i] = 2.0 * ts_xxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxz[i] = 2.0 * ts_xxz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyy[i] = 4.0 * ts_xxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyz[i] = 2.0 * ts_xxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xzz[i] = 2.0 * ts_xxz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyy[i] = 6.0 * ts_xxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyz[i] = 4.0 * ts_xxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yzz[i] = 2.0 * ts_xxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_zzz[i] = 2.0 * ts_xxz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 130-140 components of targeted buffer : FF

    auto gs_y_xyy_xxx = pbuffer.data(idx_g_ff + 130);

    auto gs_y_xyy_xxy = pbuffer.data(idx_g_ff + 131);

    auto gs_y_xyy_xxz = pbuffer.data(idx_g_ff + 132);

    auto gs_y_xyy_xyy = pbuffer.data(idx_g_ff + 133);

    auto gs_y_xyy_xyz = pbuffer.data(idx_g_ff + 134);

    auto gs_y_xyy_xzz = pbuffer.data(idx_g_ff + 135);

    auto gs_y_xyy_yyy = pbuffer.data(idx_g_ff + 136);

    auto gs_y_xyy_yyz = pbuffer.data(idx_g_ff + 137);

    auto gs_y_xyy_yzz = pbuffer.data(idx_g_ff + 138);

    auto gs_y_xyy_zzz = pbuffer.data(idx_g_ff + 139);

    #pragma omp simd aligned(gc_y, gs_y_xyy_xxx, gs_y_xyy_xxy, gs_y_xyy_xxz, gs_y_xyy_xyy, gs_y_xyy_xyz, gs_y_xyy_xzz, gs_y_xyy_yyy, gs_y_xyy_yyz, gs_y_xyy_yzz, gs_y_xyy_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyy_xxx[i] = 4.0 * ts_xy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxy[i] = 4.0 * ts_xy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxz[i] = 4.0 * ts_xy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyy[i] = 4.0 * ts_xy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyz[i] = 4.0 * ts_xy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xzz[i] = 4.0 * ts_xy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyy[i] = 4.0 * ts_xy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyz[i] = 4.0 * ts_xy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yzz[i] = 4.0 * ts_xy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_zzz[i] = 4.0 * ts_xy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 140-150 components of targeted buffer : FF

    auto gs_y_xyz_xxx = pbuffer.data(idx_g_ff + 140);

    auto gs_y_xyz_xxy = pbuffer.data(idx_g_ff + 141);

    auto gs_y_xyz_xxz = pbuffer.data(idx_g_ff + 142);

    auto gs_y_xyz_xyy = pbuffer.data(idx_g_ff + 143);

    auto gs_y_xyz_xyz = pbuffer.data(idx_g_ff + 144);

    auto gs_y_xyz_xzz = pbuffer.data(idx_g_ff + 145);

    auto gs_y_xyz_yyy = pbuffer.data(idx_g_ff + 146);

    auto gs_y_xyz_yyz = pbuffer.data(idx_g_ff + 147);

    auto gs_y_xyz_yzz = pbuffer.data(idx_g_ff + 148);

    auto gs_y_xyz_zzz = pbuffer.data(idx_g_ff + 149);

    #pragma omp simd aligned(gc_y, gs_y_xyz_xxx, gs_y_xyz_xxy, gs_y_xyz_xxz, gs_y_xyz_xyy, gs_y_xyz_xyz, gs_y_xyz_xzz, gs_y_xyz_yyy, gs_y_xyz_yyz, gs_y_xyz_yzz, gs_y_xyz_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyz_xxx[i] = 2.0 * ts_xz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxy[i] = 2.0 * ts_xz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxz[i] = 2.0 * ts_xz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyy[i] = 2.0 * ts_xz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyz[i] = 2.0 * ts_xz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xzz[i] = 2.0 * ts_xz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyy[i] = 2.0 * ts_xz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyz[i] = 2.0 * ts_xz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yzz[i] = 2.0 * ts_xz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_zzz[i] = 2.0 * ts_xz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 150-160 components of targeted buffer : FF

    auto gs_y_xzz_xxx = pbuffer.data(idx_g_ff + 150);

    auto gs_y_xzz_xxy = pbuffer.data(idx_g_ff + 151);

    auto gs_y_xzz_xxz = pbuffer.data(idx_g_ff + 152);

    auto gs_y_xzz_xyy = pbuffer.data(idx_g_ff + 153);

    auto gs_y_xzz_xyz = pbuffer.data(idx_g_ff + 154);

    auto gs_y_xzz_xzz = pbuffer.data(idx_g_ff + 155);

    auto gs_y_xzz_yyy = pbuffer.data(idx_g_ff + 156);

    auto gs_y_xzz_yyz = pbuffer.data(idx_g_ff + 157);

    auto gs_y_xzz_yzz = pbuffer.data(idx_g_ff + 158);

    auto gs_y_xzz_zzz = pbuffer.data(idx_g_ff + 159);

    #pragma omp simd aligned(gc_y, gs_y_xzz_xxx, gs_y_xzz_xxy, gs_y_xzz_xxz, gs_y_xzz_xyy, gs_y_xzz_xyz, gs_y_xzz_xzz, gs_y_xzz_yyy, gs_y_xzz_yyz, gs_y_xzz_yzz, gs_y_xzz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzz_xxx[i] = 2.0 * ts_xzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxy[i] = 2.0 * ts_xzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxz[i] = 2.0 * ts_xzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyy[i] = 4.0 * ts_xzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyz[i] = 2.0 * ts_xzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xzz[i] = 2.0 * ts_xzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyy[i] = 6.0 * ts_xzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyz[i] = 4.0 * ts_xzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yzz[i] = 2.0 * ts_xzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_zzz[i] = 2.0 * ts_xzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 160-170 components of targeted buffer : FF

    auto gs_y_yyy_xxx = pbuffer.data(idx_g_ff + 160);

    auto gs_y_yyy_xxy = pbuffer.data(idx_g_ff + 161);

    auto gs_y_yyy_xxz = pbuffer.data(idx_g_ff + 162);

    auto gs_y_yyy_xyy = pbuffer.data(idx_g_ff + 163);

    auto gs_y_yyy_xyz = pbuffer.data(idx_g_ff + 164);

    auto gs_y_yyy_xzz = pbuffer.data(idx_g_ff + 165);

    auto gs_y_yyy_yyy = pbuffer.data(idx_g_ff + 166);

    auto gs_y_yyy_yyz = pbuffer.data(idx_g_ff + 167);

    auto gs_y_yyy_yzz = pbuffer.data(idx_g_ff + 168);

    auto gs_y_yyy_zzz = pbuffer.data(idx_g_ff + 169);

    #pragma omp simd aligned(gc_y, gs_y_yyy_xxx, gs_y_yyy_xxy, gs_y_yyy_xxz, gs_y_yyy_xyy, gs_y_yyy_xyz, gs_y_yyy_xzz, gs_y_yyy_yyy, gs_y_yyy_yyz, gs_y_yyy_yzz, gs_y_yyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyy_xxx[i] = 6.0 * ts_yy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxy[i] = 6.0 * ts_yy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxz[i] = 6.0 * ts_yy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyy[i] = 6.0 * ts_yy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyz[i] = 6.0 * ts_yy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xzz[i] = 6.0 * ts_yy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyy[i] = 6.0 * ts_yy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyz[i] = 6.0 * ts_yy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yzz[i] = 6.0 * ts_yy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_zzz[i] = 6.0 * ts_yy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 170-180 components of targeted buffer : FF

    auto gs_y_yyz_xxx = pbuffer.data(idx_g_ff + 170);

    auto gs_y_yyz_xxy = pbuffer.data(idx_g_ff + 171);

    auto gs_y_yyz_xxz = pbuffer.data(idx_g_ff + 172);

    auto gs_y_yyz_xyy = pbuffer.data(idx_g_ff + 173);

    auto gs_y_yyz_xyz = pbuffer.data(idx_g_ff + 174);

    auto gs_y_yyz_xzz = pbuffer.data(idx_g_ff + 175);

    auto gs_y_yyz_yyy = pbuffer.data(idx_g_ff + 176);

    auto gs_y_yyz_yyz = pbuffer.data(idx_g_ff + 177);

    auto gs_y_yyz_yzz = pbuffer.data(idx_g_ff + 178);

    auto gs_y_yyz_zzz = pbuffer.data(idx_g_ff + 179);

    #pragma omp simd aligned(gc_y, gs_y_yyz_xxx, gs_y_yyz_xxy, gs_y_yyz_xxz, gs_y_yyz_xyy, gs_y_yyz_xyz, gs_y_yyz_xzz, gs_y_yyz_yyy, gs_y_yyz_yyz, gs_y_yyz_yzz, gs_y_yyz_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyz_xxx[i] = 4.0 * ts_yz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxy[i] = 4.0 * ts_yz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxz[i] = 4.0 * ts_yz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyy[i] = 4.0 * ts_yz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyz[i] = 4.0 * ts_yz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xzz[i] = 4.0 * ts_yz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyy[i] = 4.0 * ts_yz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyz[i] = 4.0 * ts_yz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yzz[i] = 4.0 * ts_yz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_zzz[i] = 4.0 * ts_yz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 180-190 components of targeted buffer : FF

    auto gs_y_yzz_xxx = pbuffer.data(idx_g_ff + 180);

    auto gs_y_yzz_xxy = pbuffer.data(idx_g_ff + 181);

    auto gs_y_yzz_xxz = pbuffer.data(idx_g_ff + 182);

    auto gs_y_yzz_xyy = pbuffer.data(idx_g_ff + 183);

    auto gs_y_yzz_xyz = pbuffer.data(idx_g_ff + 184);

    auto gs_y_yzz_xzz = pbuffer.data(idx_g_ff + 185);

    auto gs_y_yzz_yyy = pbuffer.data(idx_g_ff + 186);

    auto gs_y_yzz_yyz = pbuffer.data(idx_g_ff + 187);

    auto gs_y_yzz_yzz = pbuffer.data(idx_g_ff + 188);

    auto gs_y_yzz_zzz = pbuffer.data(idx_g_ff + 189);

    #pragma omp simd aligned(gc_y, gs_y_yzz_xxx, gs_y_yzz_xxy, gs_y_yzz_xxz, gs_y_yzz_xyy, gs_y_yzz_xyz, gs_y_yzz_xzz, gs_y_yzz_yyy, gs_y_yzz_yyz, gs_y_yzz_yzz, gs_y_yzz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxy[i] = 2.0 * ts_zz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyy[i] = 2.0 * ts_zz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyz[i] = 2.0 * ts_zz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 190-200 components of targeted buffer : FF

    auto gs_y_zzz_xxx = pbuffer.data(idx_g_ff + 190);

    auto gs_y_zzz_xxy = pbuffer.data(idx_g_ff + 191);

    auto gs_y_zzz_xxz = pbuffer.data(idx_g_ff + 192);

    auto gs_y_zzz_xyy = pbuffer.data(idx_g_ff + 193);

    auto gs_y_zzz_xyz = pbuffer.data(idx_g_ff + 194);

    auto gs_y_zzz_xzz = pbuffer.data(idx_g_ff + 195);

    auto gs_y_zzz_yyy = pbuffer.data(idx_g_ff + 196);

    auto gs_y_zzz_yyz = pbuffer.data(idx_g_ff + 197);

    auto gs_y_zzz_yzz = pbuffer.data(idx_g_ff + 198);

    auto gs_y_zzz_zzz = pbuffer.data(idx_g_ff + 199);

    #pragma omp simd aligned(gc_y, gs_y_zzz_xxx, gs_y_zzz_xxy, gs_y_zzz_xxz, gs_y_zzz_xyy, gs_y_zzz_xyz, gs_y_zzz_xzz, gs_y_zzz_yyy, gs_y_zzz_yyz, gs_y_zzz_yzz, gs_y_zzz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzz_xxx[i] = 2.0 * ts_zzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxy[i] = 2.0 * ts_zzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxz[i] = 2.0 * ts_zzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyy[i] = 4.0 * ts_zzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyz[i] = 2.0 * ts_zzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xzz[i] = 2.0 * ts_zzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyy[i] = 6.0 * ts_zzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyz[i] = 4.0 * ts_zzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yzz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 200-210 components of targeted buffer : FF

    auto gs_z_xxx_xxx = pbuffer.data(idx_g_ff + 200);

    auto gs_z_xxx_xxy = pbuffer.data(idx_g_ff + 201);

    auto gs_z_xxx_xxz = pbuffer.data(idx_g_ff + 202);

    auto gs_z_xxx_xyy = pbuffer.data(idx_g_ff + 203);

    auto gs_z_xxx_xyz = pbuffer.data(idx_g_ff + 204);

    auto gs_z_xxx_xzz = pbuffer.data(idx_g_ff + 205);

    auto gs_z_xxx_yyy = pbuffer.data(idx_g_ff + 206);

    auto gs_z_xxx_yyz = pbuffer.data(idx_g_ff + 207);

    auto gs_z_xxx_yzz = pbuffer.data(idx_g_ff + 208);

    auto gs_z_xxx_zzz = pbuffer.data(idx_g_ff + 209);

    #pragma omp simd aligned(gc_z, gs_z_xxx_xxx, gs_z_xxx_xxy, gs_z_xxx_xxz, gs_z_xxx_xyy, gs_z_xxx_xyz, gs_z_xxx_xzz, gs_z_xxx_yyy, gs_z_xxx_yyz, gs_z_xxx_yzz, gs_z_xxx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxx_xxx[i] = 2.0 * ts_xxx_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxy[i] = 2.0 * ts_xxx_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxz[i] = 2.0 * ts_xxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyy[i] = 2.0 * ts_xxx_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyz[i] = 2.0 * ts_xxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xzz[i] = 4.0 * ts_xxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyy[i] = 2.0 * ts_xxx_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyz[i] = 2.0 * ts_xxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yzz[i] = 4.0 * ts_xxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_zzz[i] = 6.0 * ts_xxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 210-220 components of targeted buffer : FF

    auto gs_z_xxy_xxx = pbuffer.data(idx_g_ff + 210);

    auto gs_z_xxy_xxy = pbuffer.data(idx_g_ff + 211);

    auto gs_z_xxy_xxz = pbuffer.data(idx_g_ff + 212);

    auto gs_z_xxy_xyy = pbuffer.data(idx_g_ff + 213);

    auto gs_z_xxy_xyz = pbuffer.data(idx_g_ff + 214);

    auto gs_z_xxy_xzz = pbuffer.data(idx_g_ff + 215);

    auto gs_z_xxy_yyy = pbuffer.data(idx_g_ff + 216);

    auto gs_z_xxy_yyz = pbuffer.data(idx_g_ff + 217);

    auto gs_z_xxy_yzz = pbuffer.data(idx_g_ff + 218);

    auto gs_z_xxy_zzz = pbuffer.data(idx_g_ff + 219);

    #pragma omp simd aligned(gc_z, gs_z_xxy_xxx, gs_z_xxy_xxy, gs_z_xxy_xxz, gs_z_xxy_xyy, gs_z_xxy_xyz, gs_z_xxy_xzz, gs_z_xxy_yyy, gs_z_xxy_yyz, gs_z_xxy_yzz, gs_z_xxy_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxy_xxx[i] = 2.0 * ts_xxy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxy[i] = 2.0 * ts_xxy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxz[i] = 2.0 * ts_xxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyy[i] = 2.0 * ts_xxy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyz[i] = 2.0 * ts_xxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xzz[i] = 4.0 * ts_xxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyy[i] = 2.0 * ts_xxy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyz[i] = 2.0 * ts_xxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yzz[i] = 4.0 * ts_xxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_zzz[i] = 6.0 * ts_xxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 220-230 components of targeted buffer : FF

    auto gs_z_xxz_xxx = pbuffer.data(idx_g_ff + 220);

    auto gs_z_xxz_xxy = pbuffer.data(idx_g_ff + 221);

    auto gs_z_xxz_xxz = pbuffer.data(idx_g_ff + 222);

    auto gs_z_xxz_xyy = pbuffer.data(idx_g_ff + 223);

    auto gs_z_xxz_xyz = pbuffer.data(idx_g_ff + 224);

    auto gs_z_xxz_xzz = pbuffer.data(idx_g_ff + 225);

    auto gs_z_xxz_yyy = pbuffer.data(idx_g_ff + 226);

    auto gs_z_xxz_yyz = pbuffer.data(idx_g_ff + 227);

    auto gs_z_xxz_yzz = pbuffer.data(idx_g_ff + 228);

    auto gs_z_xxz_zzz = pbuffer.data(idx_g_ff + 229);

    #pragma omp simd aligned(gc_z, gs_z_xxz_xxx, gs_z_xxz_xxy, gs_z_xxz_xxz, gs_z_xxz_xyy, gs_z_xxz_xyz, gs_z_xxz_xzz, gs_z_xxz_yyy, gs_z_xxz_yyz, gs_z_xxz_yzz, gs_z_xxz_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxz_xxx[i] = 2.0 * ts_xx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxy[i] = 2.0 * ts_xx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxz[i] = 2.0 * ts_xx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyy[i] = 2.0 * ts_xx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyz[i] = 2.0 * ts_xx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xzz[i] = 2.0 * ts_xx_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyy[i] = 2.0 * ts_xx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyz[i] = 2.0 * ts_xx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yzz[i] = 2.0 * ts_xx_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_zzz[i] = 2.0 * ts_xx_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 230-240 components of targeted buffer : FF

    auto gs_z_xyy_xxx = pbuffer.data(idx_g_ff + 230);

    auto gs_z_xyy_xxy = pbuffer.data(idx_g_ff + 231);

    auto gs_z_xyy_xxz = pbuffer.data(idx_g_ff + 232);

    auto gs_z_xyy_xyy = pbuffer.data(idx_g_ff + 233);

    auto gs_z_xyy_xyz = pbuffer.data(idx_g_ff + 234);

    auto gs_z_xyy_xzz = pbuffer.data(idx_g_ff + 235);

    auto gs_z_xyy_yyy = pbuffer.data(idx_g_ff + 236);

    auto gs_z_xyy_yyz = pbuffer.data(idx_g_ff + 237);

    auto gs_z_xyy_yzz = pbuffer.data(idx_g_ff + 238);

    auto gs_z_xyy_zzz = pbuffer.data(idx_g_ff + 239);

    #pragma omp simd aligned(gc_z, gs_z_xyy_xxx, gs_z_xyy_xxy, gs_z_xyy_xxz, gs_z_xyy_xyy, gs_z_xyy_xyz, gs_z_xyy_xzz, gs_z_xyy_yyy, gs_z_xyy_yyz, gs_z_xyy_yzz, gs_z_xyy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyy_xxx[i] = 2.0 * ts_xyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxy[i] = 2.0 * ts_xyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxz[i] = 2.0 * ts_xyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyy[i] = 2.0 * ts_xyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyz[i] = 2.0 * ts_xyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xzz[i] = 4.0 * ts_xyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyy[i] = 2.0 * ts_xyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyz[i] = 2.0 * ts_xyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yzz[i] = 4.0 * ts_xyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_zzz[i] = 6.0 * ts_xyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 240-250 components of targeted buffer : FF

    auto gs_z_xyz_xxx = pbuffer.data(idx_g_ff + 240);

    auto gs_z_xyz_xxy = pbuffer.data(idx_g_ff + 241);

    auto gs_z_xyz_xxz = pbuffer.data(idx_g_ff + 242);

    auto gs_z_xyz_xyy = pbuffer.data(idx_g_ff + 243);

    auto gs_z_xyz_xyz = pbuffer.data(idx_g_ff + 244);

    auto gs_z_xyz_xzz = pbuffer.data(idx_g_ff + 245);

    auto gs_z_xyz_yyy = pbuffer.data(idx_g_ff + 246);

    auto gs_z_xyz_yyz = pbuffer.data(idx_g_ff + 247);

    auto gs_z_xyz_yzz = pbuffer.data(idx_g_ff + 248);

    auto gs_z_xyz_zzz = pbuffer.data(idx_g_ff + 249);

    #pragma omp simd aligned(gc_z, gs_z_xyz_xxx, gs_z_xyz_xxy, gs_z_xyz_xxz, gs_z_xyz_xyy, gs_z_xyz_xyz, gs_z_xyz_xzz, gs_z_xyz_yyy, gs_z_xyz_yyz, gs_z_xyz_yzz, gs_z_xyz_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyz_xxx[i] = 2.0 * ts_xy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxy[i] = 2.0 * ts_xy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxz[i] = 2.0 * ts_xy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyy[i] = 2.0 * ts_xy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyz[i] = 2.0 * ts_xy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xzz[i] = 2.0 * ts_xy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyy[i] = 2.0 * ts_xy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyz[i] = 2.0 * ts_xy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yzz[i] = 2.0 * ts_xy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_zzz[i] = 2.0 * ts_xy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 250-260 components of targeted buffer : FF

    auto gs_z_xzz_xxx = pbuffer.data(idx_g_ff + 250);

    auto gs_z_xzz_xxy = pbuffer.data(idx_g_ff + 251);

    auto gs_z_xzz_xxz = pbuffer.data(idx_g_ff + 252);

    auto gs_z_xzz_xyy = pbuffer.data(idx_g_ff + 253);

    auto gs_z_xzz_xyz = pbuffer.data(idx_g_ff + 254);

    auto gs_z_xzz_xzz = pbuffer.data(idx_g_ff + 255);

    auto gs_z_xzz_yyy = pbuffer.data(idx_g_ff + 256);

    auto gs_z_xzz_yyz = pbuffer.data(idx_g_ff + 257);

    auto gs_z_xzz_yzz = pbuffer.data(idx_g_ff + 258);

    auto gs_z_xzz_zzz = pbuffer.data(idx_g_ff + 259);

    #pragma omp simd aligned(gc_z, gs_z_xzz_xxx, gs_z_xzz_xxy, gs_z_xzz_xxz, gs_z_xzz_xyy, gs_z_xzz_xyz, gs_z_xzz_xzz, gs_z_xzz_yyy, gs_z_xzz_yyz, gs_z_xzz_yzz, gs_z_xzz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzz_xxx[i] = 4.0 * ts_xz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxy[i] = 4.0 * ts_xz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxz[i] = 4.0 * ts_xz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyy[i] = 4.0 * ts_xz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyz[i] = 4.0 * ts_xz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xzz[i] = 4.0 * ts_xz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyy[i] = 4.0 * ts_xz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyz[i] = 4.0 * ts_xz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yzz[i] = 4.0 * ts_xz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_zzz[i] = 4.0 * ts_xz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 260-270 components of targeted buffer : FF

    auto gs_z_yyy_xxx = pbuffer.data(idx_g_ff + 260);

    auto gs_z_yyy_xxy = pbuffer.data(idx_g_ff + 261);

    auto gs_z_yyy_xxz = pbuffer.data(idx_g_ff + 262);

    auto gs_z_yyy_xyy = pbuffer.data(idx_g_ff + 263);

    auto gs_z_yyy_xyz = pbuffer.data(idx_g_ff + 264);

    auto gs_z_yyy_xzz = pbuffer.data(idx_g_ff + 265);

    auto gs_z_yyy_yyy = pbuffer.data(idx_g_ff + 266);

    auto gs_z_yyy_yyz = pbuffer.data(idx_g_ff + 267);

    auto gs_z_yyy_yzz = pbuffer.data(idx_g_ff + 268);

    auto gs_z_yyy_zzz = pbuffer.data(idx_g_ff + 269);

    #pragma omp simd aligned(gc_z, gs_z_yyy_xxx, gs_z_yyy_xxy, gs_z_yyy_xxz, gs_z_yyy_xyy, gs_z_yyy_xyz, gs_z_yyy_xzz, gs_z_yyy_yyy, gs_z_yyy_yyz, gs_z_yyy_yzz, gs_z_yyy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyy_xxx[i] = 2.0 * ts_yyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxy[i] = 2.0 * ts_yyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxz[i] = 2.0 * ts_yyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyy[i] = 2.0 * ts_yyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyz[i] = 2.0 * ts_yyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xzz[i] = 4.0 * ts_yyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyy[i] = 2.0 * ts_yyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyz[i] = 2.0 * ts_yyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yzz[i] = 4.0 * ts_yyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_zzz[i] = 6.0 * ts_yyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 270-280 components of targeted buffer : FF

    auto gs_z_yyz_xxx = pbuffer.data(idx_g_ff + 270);

    auto gs_z_yyz_xxy = pbuffer.data(idx_g_ff + 271);

    auto gs_z_yyz_xxz = pbuffer.data(idx_g_ff + 272);

    auto gs_z_yyz_xyy = pbuffer.data(idx_g_ff + 273);

    auto gs_z_yyz_xyz = pbuffer.data(idx_g_ff + 274);

    auto gs_z_yyz_xzz = pbuffer.data(idx_g_ff + 275);

    auto gs_z_yyz_yyy = pbuffer.data(idx_g_ff + 276);

    auto gs_z_yyz_yyz = pbuffer.data(idx_g_ff + 277);

    auto gs_z_yyz_yzz = pbuffer.data(idx_g_ff + 278);

    auto gs_z_yyz_zzz = pbuffer.data(idx_g_ff + 279);

    #pragma omp simd aligned(gc_z, gs_z_yyz_xxx, gs_z_yyz_xxy, gs_z_yyz_xxz, gs_z_yyz_xyy, gs_z_yyz_xyz, gs_z_yyz_xzz, gs_z_yyz_yyy, gs_z_yyz_yyz, gs_z_yyz_yzz, gs_z_yyz_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyz_xxx[i] = 2.0 * ts_yy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxy[i] = 2.0 * ts_yy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxz[i] = 2.0 * ts_yy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyy[i] = 2.0 * ts_yy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyz[i] = 2.0 * ts_yy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xzz[i] = 2.0 * ts_yy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyy[i] = 2.0 * ts_yy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyz[i] = 2.0 * ts_yy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yzz[i] = 2.0 * ts_yy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_zzz[i] = 2.0 * ts_yy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 280-290 components of targeted buffer : FF

    auto gs_z_yzz_xxx = pbuffer.data(idx_g_ff + 280);

    auto gs_z_yzz_xxy = pbuffer.data(idx_g_ff + 281);

    auto gs_z_yzz_xxz = pbuffer.data(idx_g_ff + 282);

    auto gs_z_yzz_xyy = pbuffer.data(idx_g_ff + 283);

    auto gs_z_yzz_xyz = pbuffer.data(idx_g_ff + 284);

    auto gs_z_yzz_xzz = pbuffer.data(idx_g_ff + 285);

    auto gs_z_yzz_yyy = pbuffer.data(idx_g_ff + 286);

    auto gs_z_yzz_yyz = pbuffer.data(idx_g_ff + 287);

    auto gs_z_yzz_yzz = pbuffer.data(idx_g_ff + 288);

    auto gs_z_yzz_zzz = pbuffer.data(idx_g_ff + 289);

    #pragma omp simd aligned(gc_z, gs_z_yzz_xxx, gs_z_yzz_xxy, gs_z_yzz_xxz, gs_z_yzz_xyy, gs_z_yzz_xyz, gs_z_yzz_xzz, gs_z_yzz_yyy, gs_z_yzz_yyz, gs_z_yzz_yzz, gs_z_yzz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzz_xxx[i] = 4.0 * ts_yz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxy[i] = 4.0 * ts_yz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxz[i] = 4.0 * ts_yz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyy[i] = 4.0 * ts_yz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyz[i] = 4.0 * ts_yz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xzz[i] = 4.0 * ts_yz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyy[i] = 4.0 * ts_yz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyz[i] = 4.0 * ts_yz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yzz[i] = 4.0 * ts_yz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_zzz[i] = 4.0 * ts_yz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 290-300 components of targeted buffer : FF

    auto gs_z_zzz_xxx = pbuffer.data(idx_g_ff + 290);

    auto gs_z_zzz_xxy = pbuffer.data(idx_g_ff + 291);

    auto gs_z_zzz_xxz = pbuffer.data(idx_g_ff + 292);

    auto gs_z_zzz_xyy = pbuffer.data(idx_g_ff + 293);

    auto gs_z_zzz_xyz = pbuffer.data(idx_g_ff + 294);

    auto gs_z_zzz_xzz = pbuffer.data(idx_g_ff + 295);

    auto gs_z_zzz_yyy = pbuffer.data(idx_g_ff + 296);

    auto gs_z_zzz_yyz = pbuffer.data(idx_g_ff + 297);

    auto gs_z_zzz_yzz = pbuffer.data(idx_g_ff + 298);

    auto gs_z_zzz_zzz = pbuffer.data(idx_g_ff + 299);

    #pragma omp simd aligned(gc_z, gs_z_zzz_xxx, gs_z_zzz_xxy, gs_z_zzz_xxz, gs_z_zzz_xyy, gs_z_zzz_xyz, gs_z_zzz_xzz, gs_z_zzz_yyy, gs_z_zzz_yyz, gs_z_zzz_yzz, gs_z_zzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzz_xxx[i] = 6.0 * ts_zz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxy[i] = 6.0 * ts_zz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxz[i] = 6.0 * ts_zz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyy[i] = 6.0 * ts_zz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyz[i] = 6.0 * ts_zz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xzz[i] = 6.0 * ts_zz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyy[i] = 6.0 * ts_zz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyz[i] = 6.0 * ts_zz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yzz[i] = 6.0 * ts_zz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_zzz[i] = 6.0 * ts_zz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_zzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

