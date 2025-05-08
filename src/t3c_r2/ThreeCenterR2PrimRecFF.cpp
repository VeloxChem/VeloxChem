#include "ThreeCenterR2PrimRecFF.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_ff(CSimdArray<double>& pbuffer, 
                const size_t idx_g_ff,
                const size_t idx_pf,
                const size_t idx_dd,
                const size_t idx_df,
                const size_t idx_fp,
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

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_fp);

    auto ts_xxx_y = pbuffer.data(idx_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_fp + 4);

    auto ts_xxy_z = pbuffer.data(idx_fp + 5);

    auto ts_xxz_x = pbuffer.data(idx_fp + 6);

    auto ts_xxz_y = pbuffer.data(idx_fp + 7);

    auto ts_xxz_z = pbuffer.data(idx_fp + 8);

    auto ts_xyy_x = pbuffer.data(idx_fp + 9);

    auto ts_xyy_y = pbuffer.data(idx_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_fp + 11);

    auto ts_xyz_x = pbuffer.data(idx_fp + 12);

    auto ts_xyz_y = pbuffer.data(idx_fp + 13);

    auto ts_xyz_z = pbuffer.data(idx_fp + 14);

    auto ts_xzz_x = pbuffer.data(idx_fp + 15);

    auto ts_xzz_y = pbuffer.data(idx_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_fp + 20);

    auto ts_yyz_x = pbuffer.data(idx_fp + 21);

    auto ts_yyz_y = pbuffer.data(idx_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_fp + 23);

    auto ts_yzz_x = pbuffer.data(idx_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_fp + 29);

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

    auto gr_xxx_xxx = pbuffer.data(idx_g_ff);

    auto gr_xxx_xxy = pbuffer.data(idx_g_ff + 1);

    auto gr_xxx_xxz = pbuffer.data(idx_g_ff + 2);

    auto gr_xxx_xyy = pbuffer.data(idx_g_ff + 3);

    auto gr_xxx_xyz = pbuffer.data(idx_g_ff + 4);

    auto gr_xxx_xzz = pbuffer.data(idx_g_ff + 5);

    auto gr_xxx_yyy = pbuffer.data(idx_g_ff + 6);

    auto gr_xxx_yyz = pbuffer.data(idx_g_ff + 7);

    auto gr_xxx_yzz = pbuffer.data(idx_g_ff + 8);

    auto gr_xxx_zzz = pbuffer.data(idx_g_ff + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xxx, gr_xxx_xxy, gr_xxx_xxz, gr_xxx_xyy, gr_xxx_xyz, gr_xxx_xzz, gr_xxx_yyy, gr_xxx_yyz, gr_xxx_yzz, gr_xxx_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, ts_xxx_x, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_y, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_z, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxx_xxx[i] = 6.0 * ts_x_xxx[i] * gfe2_0 + 18.0 * ts_xx_xx[i] * gfe2_0 + 6.0 * ts_xx_xxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_x[i] * gfe2_0 + 6.0 * ts_xxx_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxx_xxx[i] * gfe_0 + ts_xxx_xxx[i] * rgc2_0;

        gr_xxx_xxy[i] = 6.0 * ts_x_xxy[i] * gfe2_0 + 12.0 * ts_xx_xy[i] * gfe2_0 + 6.0 * ts_xx_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_y[i] * gfe2_0 + 4.0 * ts_xxx_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_xxy[i] * gfe_0 + ts_xxx_xxy[i] * rgc2_0;

        gr_xxx_xxz[i] = 6.0 * ts_x_xxz[i] * gfe2_0 + 12.0 * ts_xx_xz[i] * gfe2_0 + 6.0 * ts_xx_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_z[i] * gfe2_0 + 4.0 * ts_xxx_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xxz[i] * gfe_0 + ts_xxx_xxz[i] * rgc2_0;

        gr_xxx_xyy[i] = 6.0 * ts_x_xyy[i] * gfe2_0 + 6.0 * ts_xx_yy[i] * gfe2_0 + 6.0 * ts_xx_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe2_0 + 4.0 * ts_xxx_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_xyy[i] * gfe_0 + ts_xxx_xyy[i] * rgc2_0;

        gr_xxx_xyz[i] = 6.0 * ts_x_xyz[i] * gfe2_0 + 6.0 * ts_xx_yz[i] * gfe2_0 + 6.0 * ts_xx_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xyz[i] * gfe_0 + ts_xxx_xyz[i] * rgc2_0;

        gr_xxx_xzz[i] = 6.0 * ts_x_xzz[i] * gfe2_0 + 6.0 * ts_xx_zz[i] * gfe2_0 + 6.0 * ts_xx_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe2_0 + 4.0 * ts_xxx_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xzz[i] * gfe_0 + ts_xxx_xzz[i] * rgc2_0;

        gr_xxx_yyy[i] = 6.0 * ts_x_yyy[i] * gfe2_0 + 6.0 * ts_xx_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_y[i] * gfe2_0 + 6.0 * ts_xxx_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_yyy[i] * gfe_0 + ts_xxx_yyy[i] * rgc2_0;

        gr_xxx_yyz[i] = 6.0 * ts_x_yyz[i] * gfe2_0 + 6.0 * ts_xx_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_z[i] * gfe2_0 + 4.0 * ts_xxx_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_yyz[i] * gfe_0 + ts_xxx_yyz[i] * rgc2_0;

        gr_xxx_yzz[i] = 6.0 * ts_x_yzz[i] * gfe2_0 + 6.0 * ts_xx_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_y[i] * gfe2_0 + 4.0 * ts_xxx_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_yzz[i] * gfe_0 + ts_xxx_yzz[i] * rgc2_0;

        gr_xxx_zzz[i] = 6.0 * ts_x_zzz[i] * gfe2_0 + 6.0 * ts_xx_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_z[i] * gfe2_0 + 6.0 * ts_xxx_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_zzz[i] * gfe_0 + ts_xxx_zzz[i] * rgc2_0;
    }

    // Set up 10-20 components of targeted buffer : FF

    auto gr_xxy_xxx = pbuffer.data(idx_g_ff + 10);

    auto gr_xxy_xxy = pbuffer.data(idx_g_ff + 11);

    auto gr_xxy_xxz = pbuffer.data(idx_g_ff + 12);

    auto gr_xxy_xyy = pbuffer.data(idx_g_ff + 13);

    auto gr_xxy_xyz = pbuffer.data(idx_g_ff + 14);

    auto gr_xxy_xzz = pbuffer.data(idx_g_ff + 15);

    auto gr_xxy_yyy = pbuffer.data(idx_g_ff + 16);

    auto gr_xxy_yyz = pbuffer.data(idx_g_ff + 17);

    auto gr_xxy_yzz = pbuffer.data(idx_g_ff + 18);

    auto gr_xxy_zzz = pbuffer.data(idx_g_ff + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xxx, gr_xxy_xxy, gr_xxy_xxz, gr_xxy_xyy, gr_xxy_xyz, gr_xxy_xzz, gr_xxy_yyy, gr_xxy_yyz, gr_xxy_yzz, gr_xxy_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, ts_xxy_x, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_y, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_z, ts_xxy_zz, ts_xxy_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxy_xxx[i] = 2.0 * ts_y_xxx[i] * gfe2_0 + 12.0 * ts_xy_xx[i] * gfe2_0 + 4.0 * ts_xy_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_x[i] * gfe2_0 + 6.0 * ts_xxy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxy_xxx[i] * gfe_0 + ts_xxy_xxx[i] * rgc2_0;

        gr_xxy_xxy[i] = 2.0 * ts_y_xxy[i] * gfe2_0 + 8.0 * ts_xy_xy[i] * gfe2_0 + 4.0 * ts_xy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe2_0 + 2.0 * ts_xx_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_y[i] * gfe2_0 + 4.0 * ts_xxy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_xxy[i] * gfe_0 + ts_xxy_xxy[i] * rgc2_0;

        gr_xxy_xxz[i] = 2.0 * ts_y_xxz[i] * gfe2_0 + 8.0 * ts_xy_xz[i] * gfe2_0 + 4.0 * ts_xy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_z[i] * gfe2_0 + 4.0 * ts_xxy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xxz[i] * gfe_0 + ts_xxy_xxz[i] * rgc2_0;

        gr_xxy_xyy[i] = 2.0 * ts_y_xyy[i] * gfe2_0 + 4.0 * ts_xy_yy[i] * gfe2_0 + 4.0 * ts_xy_xyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_xy[i] * gfe2_0 + 2.0 * ts_xx_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_x[i] * gfe2_0 + 4.0 * ts_xxy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_xyy[i] * gfe_0 + ts_xxy_xyy[i] * rgc2_0;

        gr_xxy_xyz[i] = 2.0 * ts_y_xyz[i] * gfe2_0 + 4.0 * ts_xy_yz[i] * gfe2_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xz[i] * gfe2_0 + 2.0 * ts_xx_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xyz[i] * gfe_0 + ts_xxy_xyz[i] * rgc2_0;

        gr_xxy_xzz[i] = 2.0 * ts_y_xzz[i] * gfe2_0 + 4.0 * ts_xy_zz[i] * gfe2_0 + 4.0 * ts_xy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_x[i] * gfe2_0 + 4.0 * ts_xxy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xzz[i] * gfe_0 + ts_xxy_xzz[i] * rgc2_0;

        gr_xxy_yyy[i] = 2.0 * ts_y_yyy[i] * gfe2_0 + 4.0 * ts_xy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_yy[i] * gfe2_0 + 2.0 * ts_xx_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_y[i] * gfe2_0 + 6.0 * ts_xxy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_yyy[i] * gfe_0 + ts_xxy_yyy[i] * rgc2_0;

        gr_xxy_yyz[i] = 2.0 * ts_y_yyz[i] * gfe2_0 + 4.0 * ts_xy_yyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_yz[i] * gfe2_0 + 2.0 * ts_xx_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_z[i] * gfe2_0 + 4.0 * ts_xxy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_yyz[i] * gfe_0 + ts_xxy_yyz[i] * rgc2_0;

        gr_xxy_yzz[i] = 2.0 * ts_y_yzz[i] * gfe2_0 + 4.0 * ts_xy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe2_0 + 2.0 * ts_xx_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_y[i] * gfe2_0 + 4.0 * ts_xxy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_yzz[i] * gfe_0 + ts_xxy_yzz[i] * rgc2_0;

        gr_xxy_zzz[i] = 2.0 * ts_y_zzz[i] * gfe2_0 + 4.0 * ts_xy_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_z[i] * gfe2_0 + 6.0 * ts_xxy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_zzz[i] * gfe_0 + ts_xxy_zzz[i] * rgc2_0;
    }

    // Set up 20-30 components of targeted buffer : FF

    auto gr_xxz_xxx = pbuffer.data(idx_g_ff + 20);

    auto gr_xxz_xxy = pbuffer.data(idx_g_ff + 21);

    auto gr_xxz_xxz = pbuffer.data(idx_g_ff + 22);

    auto gr_xxz_xyy = pbuffer.data(idx_g_ff + 23);

    auto gr_xxz_xyz = pbuffer.data(idx_g_ff + 24);

    auto gr_xxz_xzz = pbuffer.data(idx_g_ff + 25);

    auto gr_xxz_yyy = pbuffer.data(idx_g_ff + 26);

    auto gr_xxz_yyz = pbuffer.data(idx_g_ff + 27);

    auto gr_xxz_yzz = pbuffer.data(idx_g_ff + 28);

    auto gr_xxz_zzz = pbuffer.data(idx_g_ff + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xxx, gr_xxz_xxy, gr_xxz_xxz, gr_xxz_xyy, gr_xxz_xyz, gr_xxz_xzz, gr_xxz_yyy, gr_xxz_yyz, gr_xxz_yzz, gr_xxz_zzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xy, ts_xx_xyy, ts_xx_xyz, ts_xx_xz, ts_xx_xzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyz, ts_xx_yz, ts_xx_yzz, ts_xx_zz, ts_xx_zzz, ts_xxz_x, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_y, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_z, ts_xxz_zz, ts_xxz_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxz_xxx[i] = 2.0 * ts_z_xxx[i] * gfe2_0 + 12.0 * ts_xz_xx[i] * gfe2_0 + 4.0 * ts_xz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxz_x[i] * gfe2_0 + 6.0 * ts_xxz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxz_xxx[i] * gfe_0 + ts_xxz_xxx[i] * rgc2_0;

        gr_xxz_xxy[i] = 2.0 * ts_z_xxy[i] * gfe2_0 + 8.0 * ts_xz_xy[i] * gfe2_0 + 4.0 * ts_xz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_y[i] * gfe2_0 + 4.0 * ts_xxz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_xxy[i] * gfe_0 + ts_xxz_xxy[i] * rgc2_0;

        gr_xxz_xxz[i] = 2.0 * ts_z_xxz[i] * gfe2_0 + 8.0 * ts_xz_xz[i] * gfe2_0 + 4.0 * ts_xz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe2_0 + 2.0 * ts_xx_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_z[i] * gfe2_0 + 4.0 * ts_xxz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xxz[i] * gfe_0 + ts_xxz_xxz[i] * rgc2_0;

        gr_xxz_xyy[i] = 2.0 * ts_z_xyy[i] * gfe2_0 + 4.0 * ts_xz_yy[i] * gfe2_0 + 4.0 * ts_xz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_x[i] * gfe2_0 + 4.0 * ts_xxz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_xyy[i] * gfe_0 + ts_xxz_xyy[i] * rgc2_0;

        gr_xxz_xyz[i] = 2.0 * ts_z_xyz[i] * gfe2_0 + 4.0 * ts_xz_yz[i] * gfe2_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xy[i] * gfe2_0 + 2.0 * ts_xx_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xyz[i] * gfe_0 + ts_xxz_xyz[i] * rgc2_0;

        gr_xxz_xzz[i] = 2.0 * ts_z_xzz[i] * gfe2_0 + 4.0 * ts_xz_zz[i] * gfe2_0 + 4.0 * ts_xz_xzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_xz[i] * gfe2_0 + 2.0 * ts_xx_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_x[i] * gfe2_0 + 4.0 * ts_xxz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xzz[i] * gfe_0 + ts_xxz_xzz[i] * rgc2_0;

        gr_xxz_yyy[i] = 2.0 * ts_z_yyy[i] * gfe2_0 + 4.0 * ts_xz_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxz_y[i] * gfe2_0 + 6.0 * ts_xxz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_yyy[i] * gfe_0 + ts_xxz_yyy[i] * rgc2_0;

        gr_xxz_yyz[i] = 2.0 * ts_z_yyz[i] * gfe2_0 + 4.0 * ts_xz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yy[i] * gfe2_0 + 2.0 * ts_xx_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_z[i] * gfe2_0 + 4.0 * ts_xxz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_yyz[i] * gfe_0 + ts_xxz_yyz[i] * rgc2_0;

        gr_xxz_yzz[i] = 2.0 * ts_z_yzz[i] * gfe2_0 + 4.0 * ts_xz_yzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_yz[i] * gfe2_0 + 2.0 * ts_xx_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_y[i] * gfe2_0 + 4.0 * ts_xxz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_yzz[i] * gfe_0 + ts_xxz_yzz[i] * rgc2_0;

        gr_xxz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe2_0 + 4.0 * ts_xz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_zz[i] * gfe2_0 + 2.0 * ts_xx_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxz_z[i] * gfe2_0 + 6.0 * ts_xxz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_zzz[i] * gfe_0 + ts_xxz_zzz[i] * rgc2_0;
    }

    // Set up 30-40 components of targeted buffer : FF

    auto gr_xyy_xxx = pbuffer.data(idx_g_ff + 30);

    auto gr_xyy_xxy = pbuffer.data(idx_g_ff + 31);

    auto gr_xyy_xxz = pbuffer.data(idx_g_ff + 32);

    auto gr_xyy_xyy = pbuffer.data(idx_g_ff + 33);

    auto gr_xyy_xyz = pbuffer.data(idx_g_ff + 34);

    auto gr_xyy_xzz = pbuffer.data(idx_g_ff + 35);

    auto gr_xyy_yyy = pbuffer.data(idx_g_ff + 36);

    auto gr_xyy_yyz = pbuffer.data(idx_g_ff + 37);

    auto gr_xyy_yzz = pbuffer.data(idx_g_ff + 38);

    auto gr_xyy_zzz = pbuffer.data(idx_g_ff + 39);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xxx, gr_xyy_xxy, gr_xyy_xxz, gr_xyy_xyy, gr_xyy_xyz, gr_xyy_xzz, gr_xyy_yyy, gr_xyy_yyz, gr_xyy_yzz, gr_xyy_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, ts_xyy_x, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_y, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_z, ts_xyy_zz, ts_xyy_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyy_xxx[i] = 6.0 * ts_yy_xx[i] * gfe2_0 + 2.0 * ts_yy_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe2_0 + 4.0 * ts_xy_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_x[i] * gfe2_0 + 6.0 * ts_xyy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyy_xxx[i] * gfe_0 + ts_xyy_xxx[i] * rgc2_0;

        gr_xyy_xxy[i] = 4.0 * ts_yy_xy[i] * gfe2_0 + 2.0 * ts_yy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxy[i] * gfe2_0 + 4.0 * ts_xy_xx[i] * gfe2_0 + 4.0 * ts_xy_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_y[i] * gfe2_0 + 4.0 * ts_xyy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_xxy[i] * gfe_0 + ts_xyy_xxy[i] * rgc2_0;

        gr_xyy_xxz[i] = 4.0 * ts_yy_xz[i] * gfe2_0 + 2.0 * ts_yy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxz[i] * gfe2_0 + 4.0 * ts_xy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_z[i] * gfe2_0 + 4.0 * ts_xyy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xxz[i] * gfe_0 + ts_xyy_xxz[i] * rgc2_0;

        gr_xyy_xyy[i] = 2.0 * ts_yy_yy[i] * gfe2_0 + 2.0 * ts_yy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyy[i] * gfe2_0 + 8.0 * ts_xy_xy[i] * gfe2_0 + 4.0 * ts_xy_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_x[i] * gfe2_0 + 4.0 * ts_xyy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_xyy[i] * gfe_0 + ts_xyy_xyy[i] * rgc2_0;

        gr_xyy_xyz[i] = 2.0 * ts_yy_yz[i] * gfe2_0 + 2.0 * ts_yy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyz[i] * gfe2_0 + 4.0 * ts_xy_xz[i] * gfe2_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xyz[i] * gfe_0 + ts_xyy_xyz[i] * rgc2_0;

        gr_xyy_xzz[i] = 2.0 * ts_yy_zz[i] * gfe2_0 + 2.0 * ts_yy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzz[i] * gfe2_0 + 4.0 * ts_xy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_x[i] * gfe2_0 + 4.0 * ts_xyy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xzz[i] * gfe_0 + ts_xyy_xzz[i] * rgc2_0;

        gr_xyy_yyy[i] = 2.0 * ts_yy_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyy[i] * gfe2_0 + 12.0 * ts_xy_yy[i] * gfe2_0 + 4.0 * ts_xy_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_y[i] * gfe2_0 + 6.0 * ts_xyy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_yyy[i] * gfe_0 + ts_xyy_yyy[i] * rgc2_0;

        gr_xyy_yyz[i] = 2.0 * ts_yy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyz[i] * gfe2_0 + 8.0 * ts_xy_yz[i] * gfe2_0 + 4.0 * ts_xy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_z[i] * gfe2_0 + 4.0 * ts_xyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_yyz[i] * gfe_0 + ts_xyy_yyz[i] * rgc2_0;

        gr_xyy_yzz[i] = 2.0 * ts_yy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yzz[i] * gfe2_0 + 4.0 * ts_xy_zz[i] * gfe2_0 + 4.0 * ts_xy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_y[i] * gfe2_0 + 4.0 * ts_xyy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_yzz[i] * gfe_0 + ts_xyy_yzz[i] * rgc2_0;

        gr_xyy_zzz[i] = 2.0 * ts_yy_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzz[i] * gfe2_0 + 4.0 * ts_xy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_z[i] * gfe2_0 + 6.0 * ts_xyy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_zzz[i] * gfe_0 + ts_xyy_zzz[i] * rgc2_0;
    }

    // Set up 40-50 components of targeted buffer : FF

    auto gr_xyz_xxx = pbuffer.data(idx_g_ff + 40);

    auto gr_xyz_xxy = pbuffer.data(idx_g_ff + 41);

    auto gr_xyz_xxz = pbuffer.data(idx_g_ff + 42);

    auto gr_xyz_xyy = pbuffer.data(idx_g_ff + 43);

    auto gr_xyz_xyz = pbuffer.data(idx_g_ff + 44);

    auto gr_xyz_xzz = pbuffer.data(idx_g_ff + 45);

    auto gr_xyz_yyy = pbuffer.data(idx_g_ff + 46);

    auto gr_xyz_yyz = pbuffer.data(idx_g_ff + 47);

    auto gr_xyz_yzz = pbuffer.data(idx_g_ff + 48);

    auto gr_xyz_zzz = pbuffer.data(idx_g_ff + 49);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xxx, gr_xyz_xxy, gr_xyz_xxz, gr_xyz_xyy, gr_xyz_xyz, gr_xyz_xzz, gr_xyz_yyy, gr_xyz_yyz, gr_xyz_yzz, gr_xyz_zzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xy, ts_xy_xyy, ts_xy_xyz, ts_xy_xz, ts_xy_xzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyz, ts_xy_yz, ts_xy_yzz, ts_xy_zz, ts_xy_zzz, ts_xyz_x, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_y, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_z, ts_xyz_zz, ts_xyz_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyz_xxx[i] = 6.0 * ts_yz_xx[i] * gfe2_0 + 2.0 * ts_yz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyz_x[i] * gfe2_0 + 6.0 * ts_xyz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyz_xxx[i] * gfe_0 + ts_xyz_xxx[i] * rgc2_0;

        gr_xyz_xxy[i] = 4.0 * ts_yz_xy[i] * gfe2_0 + 2.0 * ts_yz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe2_0 + 2.0 * ts_xz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_y[i] * gfe2_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_xxy[i] * gfe_0 + ts_xyz_xxy[i] * rgc2_0;

        gr_xyz_xxz[i] = 4.0 * ts_yz_xz[i] * gfe2_0 + 2.0 * ts_yz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xx[i] * gfe2_0 + 2.0 * ts_xy_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_z[i] * gfe2_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xxz[i] * gfe_0 + ts_xyz_xxz[i] * rgc2_0;

        gr_xyz_xyy[i] = 2.0 * ts_yz_yy[i] * gfe2_0 + 2.0 * ts_yz_xyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xz_xy[i] * gfe2_0 + 2.0 * ts_xz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_x[i] * gfe2_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_xyy[i] * gfe_0 + ts_xyz_xyy[i] * rgc2_0;

        gr_xyz_xyz[i] = 2.0 * ts_yz_yz[i] * gfe2_0 + 2.0 * ts_yz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xz[i] * gfe2_0 + 2.0 * ts_xz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xy[i] * gfe2_0 + 2.0 * ts_xy_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xyz[i] * gfe_0 + ts_xyz_xyz[i] * rgc2_0;

        gr_xyz_xzz[i] = 2.0 * ts_yz_zz[i] * gfe2_0 + 2.0 * ts_yz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xy_xz[i] * gfe2_0 + 2.0 * ts_xy_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_x[i] * gfe2_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xzz[i] * gfe_0 + ts_xyz_xzz[i] * rgc2_0;

        gr_xyz_yyy[i] = 2.0 * ts_yz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yy[i] * gfe2_0 + 2.0 * ts_xz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyz_y[i] * gfe2_0 + 6.0 * ts_xyz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_yyy[i] * gfe_0 + ts_xyz_yyy[i] * rgc2_0;

        gr_xyz_yyz[i] = 2.0 * ts_yz_yyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xz_yz[i] * gfe2_0 + 2.0 * ts_xz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe2_0 + 2.0 * ts_xy_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_z[i] * gfe2_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_yyz[i] * gfe_0 + ts_xyz_yyz[i] * rgc2_0;

        gr_xyz_yzz[i] = 2.0 * ts_yz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zz[i] * gfe2_0 + 2.0 * ts_xz_yzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xy_yz[i] * gfe2_0 + 2.0 * ts_xy_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_y[i] * gfe2_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_yzz[i] * gfe_0 + ts_xyz_yzz[i] * rgc2_0;

        gr_xyz_zzz[i] = 2.0 * ts_yz_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_zz[i] * gfe2_0 + 2.0 * ts_xy_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyz_z[i] * gfe2_0 + 6.0 * ts_xyz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_zzz[i] * gfe_0 + ts_xyz_zzz[i] * rgc2_0;
    }

    // Set up 50-60 components of targeted buffer : FF

    auto gr_xzz_xxx = pbuffer.data(idx_g_ff + 50);

    auto gr_xzz_xxy = pbuffer.data(idx_g_ff + 51);

    auto gr_xzz_xxz = pbuffer.data(idx_g_ff + 52);

    auto gr_xzz_xyy = pbuffer.data(idx_g_ff + 53);

    auto gr_xzz_xyz = pbuffer.data(idx_g_ff + 54);

    auto gr_xzz_xzz = pbuffer.data(idx_g_ff + 55);

    auto gr_xzz_yyy = pbuffer.data(idx_g_ff + 56);

    auto gr_xzz_yyz = pbuffer.data(idx_g_ff + 57);

    auto gr_xzz_yzz = pbuffer.data(idx_g_ff + 58);

    auto gr_xzz_zzz = pbuffer.data(idx_g_ff + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xxx, gr_xzz_xxy, gr_xzz_xxz, gr_xzz_xyy, gr_xzz_xyz, gr_xzz_xzz, gr_xzz_yyy, gr_xzz_yyz, gr_xzz_yzz, gr_xzz_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xy, ts_xz_xyy, ts_xz_xyz, ts_xz_xz, ts_xz_xzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyz, ts_xz_yz, ts_xz_yzz, ts_xz_zz, ts_xz_zzz, ts_xzz_x, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_y, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_z, ts_xzz_zz, ts_xzz_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xzz_xxx[i] = 6.0 * ts_zz_xx[i] * gfe2_0 + 2.0 * ts_zz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe2_0 + 4.0 * ts_xz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzz_x[i] * gfe2_0 + 6.0 * ts_xzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzz_xxx[i] * gfe_0 + ts_xzz_xxx[i] * rgc2_0;

        gr_xzz_xxy[i] = 4.0 * ts_zz_xy[i] * gfe2_0 + 2.0 * ts_zz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxy[i] * gfe2_0 + 4.0 * ts_xz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_y[i] * gfe2_0 + 4.0 * ts_xzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_xxy[i] * gfe_0 + ts_xzz_xxy[i] * rgc2_0;

        gr_xzz_xxz[i] = 4.0 * ts_zz_xz[i] * gfe2_0 + 2.0 * ts_zz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxz[i] * gfe2_0 + 4.0 * ts_xz_xx[i] * gfe2_0 + 4.0 * ts_xz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_z[i] * gfe2_0 + 4.0 * ts_xzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xxz[i] * gfe_0 + ts_xzz_xxz[i] * rgc2_0;

        gr_xzz_xyy[i] = 2.0 * ts_zz_yy[i] * gfe2_0 + 2.0 * ts_zz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyy[i] * gfe2_0 + 4.0 * ts_xz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_x[i] * gfe2_0 + 4.0 * ts_xzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_xyy[i] * gfe_0 + ts_xzz_xyy[i] * rgc2_0;

        gr_xzz_xyz[i] = 2.0 * ts_zz_yz[i] * gfe2_0 + 2.0 * ts_zz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyz[i] * gfe2_0 + 4.0 * ts_xz_xy[i] * gfe2_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xyz[i] * gfe_0 + ts_xzz_xyz[i] * rgc2_0;

        gr_xzz_xzz[i] = 2.0 * ts_zz_zz[i] * gfe2_0 + 2.0 * ts_zz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzz[i] * gfe2_0 + 8.0 * ts_xz_xz[i] * gfe2_0 + 4.0 * ts_xz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_x[i] * gfe2_0 + 4.0 * ts_xzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xzz[i] * gfe_0 + ts_xzz_xzz[i] * rgc2_0;

        gr_xzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyy[i] * gfe2_0 + 4.0 * ts_xz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzz_y[i] * gfe2_0 + 6.0 * ts_xzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_yyy[i] * gfe_0 + ts_xzz_yyy[i] * rgc2_0;

        gr_xzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyz[i] * gfe2_0 + 4.0 * ts_xz_yy[i] * gfe2_0 + 4.0 * ts_xz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_z[i] * gfe2_0 + 4.0 * ts_xzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_yyz[i] * gfe_0 + ts_xzz_yyz[i] * rgc2_0;

        gr_xzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yzz[i] * gfe2_0 + 8.0 * ts_xz_yz[i] * gfe2_0 + 4.0 * ts_xz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_y[i] * gfe2_0 + 4.0 * ts_xzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_yzz[i] * gfe_0 + ts_xzz_yzz[i] * rgc2_0;

        gr_xzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzz[i] * gfe2_0 + 12.0 * ts_xz_zz[i] * gfe2_0 + 4.0 * ts_xz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzz_z[i] * gfe2_0 + 6.0 * ts_xzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_zzz[i] * gfe_0 + ts_xzz_zzz[i] * rgc2_0;
    }

    // Set up 60-70 components of targeted buffer : FF

    auto gr_yyy_xxx = pbuffer.data(idx_g_ff + 60);

    auto gr_yyy_xxy = pbuffer.data(idx_g_ff + 61);

    auto gr_yyy_xxz = pbuffer.data(idx_g_ff + 62);

    auto gr_yyy_xyy = pbuffer.data(idx_g_ff + 63);

    auto gr_yyy_xyz = pbuffer.data(idx_g_ff + 64);

    auto gr_yyy_xzz = pbuffer.data(idx_g_ff + 65);

    auto gr_yyy_yyy = pbuffer.data(idx_g_ff + 66);

    auto gr_yyy_yyz = pbuffer.data(idx_g_ff + 67);

    auto gr_yyy_yzz = pbuffer.data(idx_g_ff + 68);

    auto gr_yyy_zzz = pbuffer.data(idx_g_ff + 69);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xxx, gr_yyy_xxy, gr_yyy_xxz, gr_yyy_xyy, gr_yyy_xyz, gr_yyy_xzz, gr_yyy_yyy, gr_yyy_yyz, gr_yyy_yzz, gr_yyy_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, ts_yyy_x, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_y, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_z, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyy_xxx[i] = 6.0 * ts_y_xxx[i] * gfe2_0 + 6.0 * ts_yy_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_x[i] * gfe2_0 + 6.0 * ts_yyy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyy_xxx[i] * gfe_0 + ts_yyy_xxx[i] * rgc2_0;

        gr_yyy_xxy[i] = 6.0 * ts_y_xxy[i] * gfe2_0 + 6.0 * ts_yy_xx[i] * gfe2_0 + 6.0 * ts_yy_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_y[i] * gfe2_0 + 4.0 * ts_yyy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_xxy[i] * gfe_0 + ts_yyy_xxy[i] * rgc2_0;

        gr_yyy_xxz[i] = 6.0 * ts_y_xxz[i] * gfe2_0 + 6.0 * ts_yy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_z[i] * gfe2_0 + 4.0 * ts_yyy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xxz[i] * gfe_0 + ts_yyy_xxz[i] * rgc2_0;

        gr_yyy_xyy[i] = 6.0 * ts_y_xyy[i] * gfe2_0 + 12.0 * ts_yy_xy[i] * gfe2_0 + 6.0 * ts_yy_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_x[i] * gfe2_0 + 4.0 * ts_yyy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_xyy[i] * gfe_0 + ts_yyy_xyy[i] * rgc2_0;

        gr_yyy_xyz[i] = 6.0 * ts_y_xyz[i] * gfe2_0 + 6.0 * ts_yy_xz[i] * gfe2_0 + 6.0 * ts_yy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xyz[i] * gfe_0 + ts_yyy_xyz[i] * rgc2_0;

        gr_yyy_xzz[i] = 6.0 * ts_y_xzz[i] * gfe2_0 + 6.0 * ts_yy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_x[i] * gfe2_0 + 4.0 * ts_yyy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xzz[i] * gfe_0 + ts_yyy_xzz[i] * rgc2_0;

        gr_yyy_yyy[i] = 6.0 * ts_y_yyy[i] * gfe2_0 + 18.0 * ts_yy_yy[i] * gfe2_0 + 6.0 * ts_yy_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_y[i] * gfe2_0 + 6.0 * ts_yyy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_yyy[i] * gfe_0 + ts_yyy_yyy[i] * rgc2_0;

        gr_yyy_yyz[i] = 6.0 * ts_y_yyz[i] * gfe2_0 + 12.0 * ts_yy_yz[i] * gfe2_0 + 6.0 * ts_yy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_z[i] * gfe2_0 + 4.0 * ts_yyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_yyz[i] * gfe_0 + ts_yyy_yyz[i] * rgc2_0;

        gr_yyy_yzz[i] = 6.0 * ts_y_yzz[i] * gfe2_0 + 6.0 * ts_yy_zz[i] * gfe2_0 + 6.0 * ts_yy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_y[i] * gfe2_0 + 4.0 * ts_yyy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_yzz[i] * gfe_0 + ts_yyy_yzz[i] * rgc2_0;

        gr_yyy_zzz[i] = 6.0 * ts_y_zzz[i] * gfe2_0 + 6.0 * ts_yy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_z[i] * gfe2_0 + 6.0 * ts_yyy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_zzz[i] * gfe_0 + ts_yyy_zzz[i] * rgc2_0;
    }

    // Set up 70-80 components of targeted buffer : FF

    auto gr_yyz_xxx = pbuffer.data(idx_g_ff + 70);

    auto gr_yyz_xxy = pbuffer.data(idx_g_ff + 71);

    auto gr_yyz_xxz = pbuffer.data(idx_g_ff + 72);

    auto gr_yyz_xyy = pbuffer.data(idx_g_ff + 73);

    auto gr_yyz_xyz = pbuffer.data(idx_g_ff + 74);

    auto gr_yyz_xzz = pbuffer.data(idx_g_ff + 75);

    auto gr_yyz_yyy = pbuffer.data(idx_g_ff + 76);

    auto gr_yyz_yyz = pbuffer.data(idx_g_ff + 77);

    auto gr_yyz_yzz = pbuffer.data(idx_g_ff + 78);

    auto gr_yyz_zzz = pbuffer.data(idx_g_ff + 79);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xxx, gr_yyz_xxy, gr_yyz_xxz, gr_yyz_xyy, gr_yyz_xyz, gr_yyz_xzz, gr_yyz_yyy, gr_yyz_yyz, gr_yyz_yzz, gr_yyz_zzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xy, ts_yy_xyy, ts_yy_xyz, ts_yy_xz, ts_yy_xzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyz, ts_yy_yz, ts_yy_yzz, ts_yy_zz, ts_yy_zzz, ts_yyz_x, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_y, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_z, ts_yyz_zz, ts_yyz_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyz_xxx[i] = 2.0 * ts_z_xxx[i] * gfe2_0 + 4.0 * ts_yz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyz_x[i] * gfe2_0 + 6.0 * ts_yyz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyz_xxx[i] * gfe_0 + ts_yyz_xxx[i] * rgc2_0;

        gr_yyz_xxy[i] = 2.0 * ts_z_xxy[i] * gfe2_0 + 4.0 * ts_yz_xx[i] * gfe2_0 + 4.0 * ts_yz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_y[i] * gfe2_0 + 4.0 * ts_yyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_xxy[i] * gfe_0 + ts_yyz_xxy[i] * rgc2_0;

        gr_yyz_xxz[i] = 2.0 * ts_z_xxz[i] * gfe2_0 + 4.0 * ts_yz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xx[i] * gfe2_0 + 2.0 * ts_yy_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_z[i] * gfe2_0 + 4.0 * ts_yyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xxz[i] * gfe_0 + ts_yyz_xxz[i] * rgc2_0;

        gr_yyz_xyy[i] = 2.0 * ts_z_xyy[i] * gfe2_0 + 8.0 * ts_yz_xy[i] * gfe2_0 + 4.0 * ts_yz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_x[i] * gfe2_0 + 4.0 * ts_yyz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_xyy[i] * gfe_0 + ts_yyz_xyy[i] * rgc2_0;

        gr_yyz_xyz[i] = 2.0 * ts_z_xyz[i] * gfe2_0 + 4.0 * ts_yz_xz[i] * gfe2_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xy[i] * gfe2_0 + 2.0 * ts_yy_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xyz[i] * gfe_0 + ts_yyz_xyz[i] * rgc2_0;

        gr_yyz_xzz[i] = 2.0 * ts_z_xzz[i] * gfe2_0 + 4.0 * ts_yz_xzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yy_xz[i] * gfe2_0 + 2.0 * ts_yy_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_x[i] * gfe2_0 + 4.0 * ts_yyz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xzz[i] * gfe_0 + ts_yyz_xzz[i] * rgc2_0;

        gr_yyz_yyy[i] = 2.0 * ts_z_yyy[i] * gfe2_0 + 12.0 * ts_yz_yy[i] * gfe2_0 + 4.0 * ts_yz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyz_y[i] * gfe2_0 + 6.0 * ts_yyz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_yyy[i] * gfe_0 + ts_yyz_yyy[i] * rgc2_0;

        gr_yyz_yyz[i] = 2.0 * ts_z_yyz[i] * gfe2_0 + 8.0 * ts_yz_yz[i] * gfe2_0 + 4.0 * ts_yz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe2_0 + 2.0 * ts_yy_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_z[i] * gfe2_0 + 4.0 * ts_yyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_yyz[i] * gfe_0 + ts_yyz_yyz[i] * rgc2_0;

        gr_yyz_yzz[i] = 2.0 * ts_z_yzz[i] * gfe2_0 + 4.0 * ts_yz_zz[i] * gfe2_0 + 4.0 * ts_yz_yzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yy_yz[i] * gfe2_0 + 2.0 * ts_yy_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_y[i] * gfe2_0 + 4.0 * ts_yyz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_yzz[i] * gfe_0 + ts_yyz_yzz[i] * rgc2_0;

        gr_yyz_zzz[i] = 2.0 * ts_z_zzz[i] * gfe2_0 + 4.0 * ts_yz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_zz[i] * gfe2_0 + 2.0 * ts_yy_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyz_z[i] * gfe2_0 + 6.0 * ts_yyz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_zzz[i] * gfe_0 + ts_yyz_zzz[i] * rgc2_0;
    }

    // Set up 80-90 components of targeted buffer : FF

    auto gr_yzz_xxx = pbuffer.data(idx_g_ff + 80);

    auto gr_yzz_xxy = pbuffer.data(idx_g_ff + 81);

    auto gr_yzz_xxz = pbuffer.data(idx_g_ff + 82);

    auto gr_yzz_xyy = pbuffer.data(idx_g_ff + 83);

    auto gr_yzz_xyz = pbuffer.data(idx_g_ff + 84);

    auto gr_yzz_xzz = pbuffer.data(idx_g_ff + 85);

    auto gr_yzz_yyy = pbuffer.data(idx_g_ff + 86);

    auto gr_yzz_yyz = pbuffer.data(idx_g_ff + 87);

    auto gr_yzz_yzz = pbuffer.data(idx_g_ff + 88);

    auto gr_yzz_zzz = pbuffer.data(idx_g_ff + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xxx, gr_yzz_xxy, gr_yzz_xxz, gr_yzz_xyy, gr_yzz_xyz, gr_yzz_xzz, gr_yzz_yyy, gr_yzz_yyz, gr_yzz_yzz, gr_yzz_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xy, ts_yz_xyy, ts_yz_xyz, ts_yz_xz, ts_yz_xzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyz, ts_yz_yz, ts_yz_yzz, ts_yz_zz, ts_yz_zzz, ts_yzz_x, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_y, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_z, ts_yzz_zz, ts_yzz_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxx[i] * gfe2_0 + 4.0 * ts_yz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzz_x[i] * gfe2_0 + 6.0 * ts_yzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzz_xxx[i] * gfe_0 + ts_yzz_xxx[i] * rgc2_0;

        gr_yzz_xxy[i] = 2.0 * ts_zz_xx[i] * gfe2_0 + 2.0 * ts_zz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxy[i] * gfe2_0 + 4.0 * ts_yz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_y[i] * gfe2_0 + 4.0 * ts_yzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_xxy[i] * gfe_0 + ts_yzz_xxy[i] * rgc2_0;

        gr_yzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxz[i] * gfe2_0 + 4.0 * ts_yz_xx[i] * gfe2_0 + 4.0 * ts_yz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_z[i] * gfe2_0 + 4.0 * ts_yzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xxz[i] * gfe_0 + ts_yzz_xxz[i] * rgc2_0;

        gr_yzz_xyy[i] = 4.0 * ts_zz_xy[i] * gfe2_0 + 2.0 * ts_zz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyy[i] * gfe2_0 + 4.0 * ts_yz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_x[i] * gfe2_0 + 4.0 * ts_yzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_xyy[i] * gfe_0 + ts_yzz_xyy[i] * rgc2_0;

        gr_yzz_xyz[i] = 2.0 * ts_zz_xz[i] * gfe2_0 + 2.0 * ts_zz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyz[i] * gfe2_0 + 4.0 * ts_yz_xy[i] * gfe2_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xyz[i] * gfe_0 + ts_yzz_xyz[i] * rgc2_0;

        gr_yzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xzz[i] * gfe2_0 + 8.0 * ts_yz_xz[i] * gfe2_0 + 4.0 * ts_yz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_x[i] * gfe2_0 + 4.0 * ts_yzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xzz[i] * gfe_0 + ts_yzz_xzz[i] * rgc2_0;

        gr_yzz_yyy[i] = 6.0 * ts_zz_yy[i] * gfe2_0 + 2.0 * ts_zz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyy[i] * gfe2_0 + 4.0 * ts_yz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzz_y[i] * gfe2_0 + 6.0 * ts_yzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_yyy[i] * gfe_0 + ts_yzz_yyy[i] * rgc2_0;

        gr_yzz_yyz[i] = 4.0 * ts_zz_yz[i] * gfe2_0 + 2.0 * ts_zz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyz[i] * gfe2_0 + 4.0 * ts_yz_yy[i] * gfe2_0 + 4.0 * ts_yz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_z[i] * gfe2_0 + 4.0 * ts_yzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_yyz[i] * gfe_0 + ts_yzz_yyz[i] * rgc2_0;

        gr_yzz_yzz[i] = 2.0 * ts_zz_zz[i] * gfe2_0 + 2.0 * ts_zz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yzz[i] * gfe2_0 + 8.0 * ts_yz_yz[i] * gfe2_0 + 4.0 * ts_yz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_y[i] * gfe2_0 + 4.0 * ts_yzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_yzz[i] * gfe_0 + ts_yzz_yzz[i] * rgc2_0;

        gr_yzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zzz[i] * gfe2_0 + 12.0 * ts_yz_zz[i] * gfe2_0 + 4.0 * ts_yz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzz_z[i] * gfe2_0 + 6.0 * ts_yzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_zzz[i] * gfe_0 + ts_yzz_zzz[i] * rgc2_0;
    }

    // Set up 90-100 components of targeted buffer : FF

    auto gr_zzz_xxx = pbuffer.data(idx_g_ff + 90);

    auto gr_zzz_xxy = pbuffer.data(idx_g_ff + 91);

    auto gr_zzz_xxz = pbuffer.data(idx_g_ff + 92);

    auto gr_zzz_xyy = pbuffer.data(idx_g_ff + 93);

    auto gr_zzz_xyz = pbuffer.data(idx_g_ff + 94);

    auto gr_zzz_xzz = pbuffer.data(idx_g_ff + 95);

    auto gr_zzz_yyy = pbuffer.data(idx_g_ff + 96);

    auto gr_zzz_yyz = pbuffer.data(idx_g_ff + 97);

    auto gr_zzz_yzz = pbuffer.data(idx_g_ff + 98);

    auto gr_zzz_zzz = pbuffer.data(idx_g_ff + 99);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xxx, gr_zzz_xxy, gr_zzz_xxz, gr_zzz_xyy, gr_zzz_xyz, gr_zzz_xzz, gr_zzz_yyy, gr_zzz_yyz, gr_zzz_yzz, gr_zzz_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xy, ts_zz_xyy, ts_zz_xyz, ts_zz_xz, ts_zz_xzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyz, ts_zz_yz, ts_zz_yzz, ts_zz_zz, ts_zz_zzz, ts_zzz_x, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_y, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_z, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_zzz_xxx[i] = 6.0 * ts_z_xxx[i] * gfe2_0 + 6.0 * ts_zz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzz_x[i] * gfe2_0 + 6.0 * ts_zzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzz_xxx[i] * gfe_0 + ts_zzz_xxx[i] * rgc2_0;

        gr_zzz_xxy[i] = 6.0 * ts_z_xxy[i] * gfe2_0 + 6.0 * ts_zz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_y[i] * gfe2_0 + 4.0 * ts_zzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_xxy[i] * gfe_0 + ts_zzz_xxy[i] * rgc2_0;

        gr_zzz_xxz[i] = 6.0 * ts_z_xxz[i] * gfe2_0 + 6.0 * ts_zz_xx[i] * gfe2_0 + 6.0 * ts_zz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_z[i] * gfe2_0 + 4.0 * ts_zzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xxz[i] * gfe_0 + ts_zzz_xxz[i] * rgc2_0;

        gr_zzz_xyy[i] = 6.0 * ts_z_xyy[i] * gfe2_0 + 6.0 * ts_zz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_x[i] * gfe2_0 + 4.0 * ts_zzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_xyy[i] * gfe_0 + ts_zzz_xyy[i] * rgc2_0;

        gr_zzz_xyz[i] = 6.0 * ts_z_xyz[i] * gfe2_0 + 6.0 * ts_zz_xy[i] * gfe2_0 + 6.0 * ts_zz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xyz[i] * gfe_0 + ts_zzz_xyz[i] * rgc2_0;

        gr_zzz_xzz[i] = 6.0 * ts_z_xzz[i] * gfe2_0 + 12.0 * ts_zz_xz[i] * gfe2_0 + 6.0 * ts_zz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_x[i] * gfe2_0 + 4.0 * ts_zzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xzz[i] * gfe_0 + ts_zzz_xzz[i] * rgc2_0;

        gr_zzz_yyy[i] = 6.0 * ts_z_yyy[i] * gfe2_0 + 6.0 * ts_zz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzz_y[i] * gfe2_0 + 6.0 * ts_zzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_yyy[i] * gfe_0 + ts_zzz_yyy[i] * rgc2_0;

        gr_zzz_yyz[i] = 6.0 * ts_z_yyz[i] * gfe2_0 + 6.0 * ts_zz_yy[i] * gfe2_0 + 6.0 * ts_zz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_z[i] * gfe2_0 + 4.0 * ts_zzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_yyz[i] * gfe_0 + ts_zzz_yyz[i] * rgc2_0;

        gr_zzz_yzz[i] = 6.0 * ts_z_yzz[i] * gfe2_0 + 12.0 * ts_zz_yz[i] * gfe2_0 + 6.0 * ts_zz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_y[i] * gfe2_0 + 4.0 * ts_zzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_yzz[i] * gfe_0 + ts_zzz_yzz[i] * rgc2_0;

        gr_zzz_zzz[i] = 6.0 * ts_z_zzz[i] * gfe2_0 + 18.0 * ts_zz_zz[i] * gfe2_0 + 6.0 * ts_zz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzz_z[i] * gfe2_0 + 6.0 * ts_zzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_zzz[i] * gfe_0 + ts_zzz_zzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

