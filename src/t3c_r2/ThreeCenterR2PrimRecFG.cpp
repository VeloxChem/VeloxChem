#include "ThreeCenterR2PrimRecFG.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_fg(CSimdArray<double>& pbuffer, 
                const size_t idx_g_fg,
                const size_t idx_pg,
                const size_t idx_df,
                const size_t idx_dg,
                const size_t idx_fd,
                const size_t idx_ff,
                const size_t idx_fg,
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

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_pg);

    auto ts_x_xxxy = pbuffer.data(idx_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_pg + 44);

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

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto ts_xy_xxxx = pbuffer.data(idx_dg + 15);

    auto ts_xy_xxxy = pbuffer.data(idx_dg + 16);

    auto ts_xy_xxxz = pbuffer.data(idx_dg + 17);

    auto ts_xy_xxyy = pbuffer.data(idx_dg + 18);

    auto ts_xy_xxyz = pbuffer.data(idx_dg + 19);

    auto ts_xy_xxzz = pbuffer.data(idx_dg + 20);

    auto ts_xy_xyyy = pbuffer.data(idx_dg + 21);

    auto ts_xy_xyyz = pbuffer.data(idx_dg + 22);

    auto ts_xy_xyzz = pbuffer.data(idx_dg + 23);

    auto ts_xy_xzzz = pbuffer.data(idx_dg + 24);

    auto ts_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto ts_xy_zzzz = pbuffer.data(idx_dg + 29);

    auto ts_xz_xxxx = pbuffer.data(idx_dg + 30);

    auto ts_xz_xxxy = pbuffer.data(idx_dg + 31);

    auto ts_xz_xxxz = pbuffer.data(idx_dg + 32);

    auto ts_xz_xxyy = pbuffer.data(idx_dg + 33);

    auto ts_xz_xxyz = pbuffer.data(idx_dg + 34);

    auto ts_xz_xxzz = pbuffer.data(idx_dg + 35);

    auto ts_xz_xyyy = pbuffer.data(idx_dg + 36);

    auto ts_xz_xyyz = pbuffer.data(idx_dg + 37);

    auto ts_xz_xyzz = pbuffer.data(idx_dg + 38);

    auto ts_xz_xzzz = pbuffer.data(idx_dg + 39);

    auto ts_xz_yyyy = pbuffer.data(idx_dg + 40);

    auto ts_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto ts_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto ts_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto ts_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto ts_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto ts_yz_xxxx = pbuffer.data(idx_dg + 60);

    auto ts_yz_xxxy = pbuffer.data(idx_dg + 61);

    auto ts_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto ts_yz_xxyy = pbuffer.data(idx_dg + 63);

    auto ts_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto ts_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto ts_yz_xyyy = pbuffer.data(idx_dg + 66);

    auto ts_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto ts_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto ts_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto ts_yz_yyyy = pbuffer.data(idx_dg + 70);

    auto ts_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto ts_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto ts_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_dg + 89);

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

    // Set up components of auxiliary buffer : FG

    auto ts_xxx_xxxx = pbuffer.data(idx_fg);

    auto ts_xxx_xxxy = pbuffer.data(idx_fg + 1);

    auto ts_xxx_xxxz = pbuffer.data(idx_fg + 2);

    auto ts_xxx_xxyy = pbuffer.data(idx_fg + 3);

    auto ts_xxx_xxyz = pbuffer.data(idx_fg + 4);

    auto ts_xxx_xxzz = pbuffer.data(idx_fg + 5);

    auto ts_xxx_xyyy = pbuffer.data(idx_fg + 6);

    auto ts_xxx_xyyz = pbuffer.data(idx_fg + 7);

    auto ts_xxx_xyzz = pbuffer.data(idx_fg + 8);

    auto ts_xxx_xzzz = pbuffer.data(idx_fg + 9);

    auto ts_xxx_yyyy = pbuffer.data(idx_fg + 10);

    auto ts_xxx_yyyz = pbuffer.data(idx_fg + 11);

    auto ts_xxx_yyzz = pbuffer.data(idx_fg + 12);

    auto ts_xxx_yzzz = pbuffer.data(idx_fg + 13);

    auto ts_xxx_zzzz = pbuffer.data(idx_fg + 14);

    auto ts_xxy_xxxx = pbuffer.data(idx_fg + 15);

    auto ts_xxy_xxxy = pbuffer.data(idx_fg + 16);

    auto ts_xxy_xxxz = pbuffer.data(idx_fg + 17);

    auto ts_xxy_xxyy = pbuffer.data(idx_fg + 18);

    auto ts_xxy_xxyz = pbuffer.data(idx_fg + 19);

    auto ts_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto ts_xxy_xyyy = pbuffer.data(idx_fg + 21);

    auto ts_xxy_xyyz = pbuffer.data(idx_fg + 22);

    auto ts_xxy_xyzz = pbuffer.data(idx_fg + 23);

    auto ts_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto ts_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto ts_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto ts_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto ts_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto ts_xxy_zzzz = pbuffer.data(idx_fg + 29);

    auto ts_xxz_xxxx = pbuffer.data(idx_fg + 30);

    auto ts_xxz_xxxy = pbuffer.data(idx_fg + 31);

    auto ts_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto ts_xxz_xxyy = pbuffer.data(idx_fg + 33);

    auto ts_xxz_xxyz = pbuffer.data(idx_fg + 34);

    auto ts_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto ts_xxz_xyyy = pbuffer.data(idx_fg + 36);

    auto ts_xxz_xyyz = pbuffer.data(idx_fg + 37);

    auto ts_xxz_xyzz = pbuffer.data(idx_fg + 38);

    auto ts_xxz_xzzz = pbuffer.data(idx_fg + 39);

    auto ts_xxz_yyyy = pbuffer.data(idx_fg + 40);

    auto ts_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto ts_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto ts_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto ts_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto ts_xyy_xxxx = pbuffer.data(idx_fg + 45);

    auto ts_xyy_xxxy = pbuffer.data(idx_fg + 46);

    auto ts_xyy_xxxz = pbuffer.data(idx_fg + 47);

    auto ts_xyy_xxyy = pbuffer.data(idx_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_fg + 49);

    auto ts_xyy_xxzz = pbuffer.data(idx_fg + 50);

    auto ts_xyy_xyyy = pbuffer.data(idx_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_fg + 53);

    auto ts_xyy_xzzz = pbuffer.data(idx_fg + 54);

    auto ts_xyy_yyyy = pbuffer.data(idx_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_fg + 59);

    auto ts_xyz_xxxx = pbuffer.data(idx_fg + 60);

    auto ts_xyz_xxxy = pbuffer.data(idx_fg + 61);

    auto ts_xyz_xxxz = pbuffer.data(idx_fg + 62);

    auto ts_xyz_xxyy = pbuffer.data(idx_fg + 63);

    auto ts_xyz_xxyz = pbuffer.data(idx_fg + 64);

    auto ts_xyz_xxzz = pbuffer.data(idx_fg + 65);

    auto ts_xyz_xyyy = pbuffer.data(idx_fg + 66);

    auto ts_xyz_xyyz = pbuffer.data(idx_fg + 67);

    auto ts_xyz_xyzz = pbuffer.data(idx_fg + 68);

    auto ts_xyz_xzzz = pbuffer.data(idx_fg + 69);

    auto ts_xyz_yyyy = pbuffer.data(idx_fg + 70);

    auto ts_xyz_yyyz = pbuffer.data(idx_fg + 71);

    auto ts_xyz_yyzz = pbuffer.data(idx_fg + 72);

    auto ts_xyz_yzzz = pbuffer.data(idx_fg + 73);

    auto ts_xyz_zzzz = pbuffer.data(idx_fg + 74);

    auto ts_xzz_xxxx = pbuffer.data(idx_fg + 75);

    auto ts_xzz_xxxy = pbuffer.data(idx_fg + 76);

    auto ts_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto ts_xzz_xxyy = pbuffer.data(idx_fg + 78);

    auto ts_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto ts_xzz_xyyy = pbuffer.data(idx_fg + 81);

    auto ts_xzz_xyyz = pbuffer.data(idx_fg + 82);

    auto ts_xzz_xyzz = pbuffer.data(idx_fg + 83);

    auto ts_xzz_xzzz = pbuffer.data(idx_fg + 84);

    auto ts_xzz_yyyy = pbuffer.data(idx_fg + 85);

    auto ts_xzz_yyyz = pbuffer.data(idx_fg + 86);

    auto ts_xzz_yyzz = pbuffer.data(idx_fg + 87);

    auto ts_xzz_yzzz = pbuffer.data(idx_fg + 88);

    auto ts_xzz_zzzz = pbuffer.data(idx_fg + 89);

    auto ts_yyy_xxxx = pbuffer.data(idx_fg + 90);

    auto ts_yyy_xxxy = pbuffer.data(idx_fg + 91);

    auto ts_yyy_xxxz = pbuffer.data(idx_fg + 92);

    auto ts_yyy_xxyy = pbuffer.data(idx_fg + 93);

    auto ts_yyy_xxyz = pbuffer.data(idx_fg + 94);

    auto ts_yyy_xxzz = pbuffer.data(idx_fg + 95);

    auto ts_yyy_xyyy = pbuffer.data(idx_fg + 96);

    auto ts_yyy_xyyz = pbuffer.data(idx_fg + 97);

    auto ts_yyy_xyzz = pbuffer.data(idx_fg + 98);

    auto ts_yyy_xzzz = pbuffer.data(idx_fg + 99);

    auto ts_yyy_yyyy = pbuffer.data(idx_fg + 100);

    auto ts_yyy_yyyz = pbuffer.data(idx_fg + 101);

    auto ts_yyy_yyzz = pbuffer.data(idx_fg + 102);

    auto ts_yyy_yzzz = pbuffer.data(idx_fg + 103);

    auto ts_yyy_zzzz = pbuffer.data(idx_fg + 104);

    auto ts_yyz_xxxx = pbuffer.data(idx_fg + 105);

    auto ts_yyz_xxxy = pbuffer.data(idx_fg + 106);

    auto ts_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto ts_yyz_xxyy = pbuffer.data(idx_fg + 108);

    auto ts_yyz_xxyz = pbuffer.data(idx_fg + 109);

    auto ts_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto ts_yyz_xyyy = pbuffer.data(idx_fg + 111);

    auto ts_yyz_xyyz = pbuffer.data(idx_fg + 112);

    auto ts_yyz_xyzz = pbuffer.data(idx_fg + 113);

    auto ts_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto ts_yyz_yyyy = pbuffer.data(idx_fg + 115);

    auto ts_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_fg + 119);

    auto ts_yzz_xxxx = pbuffer.data(idx_fg + 120);

    auto ts_yzz_xxxy = pbuffer.data(idx_fg + 121);

    auto ts_yzz_xxxz = pbuffer.data(idx_fg + 122);

    auto ts_yzz_xxyy = pbuffer.data(idx_fg + 123);

    auto ts_yzz_xxyz = pbuffer.data(idx_fg + 124);

    auto ts_yzz_xxzz = pbuffer.data(idx_fg + 125);

    auto ts_yzz_xyyy = pbuffer.data(idx_fg + 126);

    auto ts_yzz_xyyz = pbuffer.data(idx_fg + 127);

    auto ts_yzz_xyzz = pbuffer.data(idx_fg + 128);

    auto ts_yzz_xzzz = pbuffer.data(idx_fg + 129);

    auto ts_yzz_yyyy = pbuffer.data(idx_fg + 130);

    auto ts_yzz_yyyz = pbuffer.data(idx_fg + 131);

    auto ts_yzz_yyzz = pbuffer.data(idx_fg + 132);

    auto ts_yzz_yzzz = pbuffer.data(idx_fg + 133);

    auto ts_yzz_zzzz = pbuffer.data(idx_fg + 134);

    auto ts_zzz_xxxx = pbuffer.data(idx_fg + 135);

    auto ts_zzz_xxxy = pbuffer.data(idx_fg + 136);

    auto ts_zzz_xxxz = pbuffer.data(idx_fg + 137);

    auto ts_zzz_xxyy = pbuffer.data(idx_fg + 138);

    auto ts_zzz_xxyz = pbuffer.data(idx_fg + 139);

    auto ts_zzz_xxzz = pbuffer.data(idx_fg + 140);

    auto ts_zzz_xyyy = pbuffer.data(idx_fg + 141);

    auto ts_zzz_xyyz = pbuffer.data(idx_fg + 142);

    auto ts_zzz_xyzz = pbuffer.data(idx_fg + 143);

    auto ts_zzz_xzzz = pbuffer.data(idx_fg + 144);

    auto ts_zzz_yyyy = pbuffer.data(idx_fg + 145);

    auto ts_zzz_yyyz = pbuffer.data(idx_fg + 146);

    auto ts_zzz_yyzz = pbuffer.data(idx_fg + 147);

    auto ts_zzz_yzzz = pbuffer.data(idx_fg + 148);

    auto ts_zzz_zzzz = pbuffer.data(idx_fg + 149);

    // Set up 0-15 components of targeted buffer : FG

    auto gr_xxx_xxxx = pbuffer.data(idx_g_fg);

    auto gr_xxx_xxxy = pbuffer.data(idx_g_fg + 1);

    auto gr_xxx_xxxz = pbuffer.data(idx_g_fg + 2);

    auto gr_xxx_xxyy = pbuffer.data(idx_g_fg + 3);

    auto gr_xxx_xxyz = pbuffer.data(idx_g_fg + 4);

    auto gr_xxx_xxzz = pbuffer.data(idx_g_fg + 5);

    auto gr_xxx_xyyy = pbuffer.data(idx_g_fg + 6);

    auto gr_xxx_xyyz = pbuffer.data(idx_g_fg + 7);

    auto gr_xxx_xyzz = pbuffer.data(idx_g_fg + 8);

    auto gr_xxx_xzzz = pbuffer.data(idx_g_fg + 9);

    auto gr_xxx_yyyy = pbuffer.data(idx_g_fg + 10);

    auto gr_xxx_yyyz = pbuffer.data(idx_g_fg + 11);

    auto gr_xxx_yyzz = pbuffer.data(idx_g_fg + 12);

    auto gr_xxx_yzzz = pbuffer.data(idx_g_fg + 13);

    auto gr_xxx_zzzz = pbuffer.data(idx_g_fg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xxxx, gr_xxx_xxxy, gr_xxx_xxxz, gr_xxx_xxyy, gr_xxx_xxyz, gr_xxx_xxzz, gr_xxx_xyyy, gr_xxx_xyyz, gr_xxx_xyzz, gr_xxx_xzzz, gr_xxx_yyyy, gr_xxx_yyyz, gr_xxx_yyzz, gr_xxx_yzzz, gr_xxx_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxx_xxxx[i] = 6.0 * ts_x_xxxx[i] * gfe2_0 + 24.0 * ts_xx_xxx[i] * gfe2_0 + 6.0 * ts_xx_xxxx[i] * gfe_0 * gc_x[i] + 12.0 * ts_xxx_xx[i] * gfe2_0 + 8.0 * ts_xxx_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxx_xxxx[i] * gfe_0 + ts_xxx_xxxx[i] * rgc2_0;

        gr_xxx_xxxy[i] = 6.0 * ts_x_xxxy[i] * gfe2_0 + 18.0 * ts_xx_xxy[i] * gfe2_0 + 6.0 * ts_xx_xxxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_xy[i] * gfe2_0 + 6.0 * ts_xxx_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_xxxy[i] * gfe_0 + ts_xxx_xxxy[i] * rgc2_0;

        gr_xxx_xxxz[i] = 6.0 * ts_x_xxxz[i] * gfe2_0 + 18.0 * ts_xx_xxz[i] * gfe2_0 + 6.0 * ts_xx_xxxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_xz[i] * gfe2_0 + 6.0 * ts_xxx_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xxxz[i] * gfe_0 + ts_xxx_xxxz[i] * rgc2_0;

        gr_xxx_xxyy[i] = 6.0 * ts_x_xxyy[i] * gfe2_0 + 12.0 * ts_xx_xyy[i] * gfe2_0 + 6.0 * ts_xx_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yy[i] * gfe2_0 + 4.0 * ts_xxx_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe2_0 + 4.0 * ts_xxx_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_xxyy[i] * gfe_0 + ts_xxx_xxyy[i] * rgc2_0;

        gr_xxx_xxyz[i] = 6.0 * ts_x_xxyz[i] * gfe2_0 + 12.0 * ts_xx_xyz[i] * gfe2_0 + 6.0 * ts_xx_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yz[i] * gfe2_0 + 4.0 * ts_xxx_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xxyz[i] * gfe_0 + ts_xxx_xxyz[i] * rgc2_0;

        gr_xxx_xxzz[i] = 6.0 * ts_x_xxzz[i] * gfe2_0 + 12.0 * ts_xx_xzz[i] * gfe2_0 + 6.0 * ts_xx_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zz[i] * gfe2_0 + 4.0 * ts_xxx_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe2_0 + 4.0 * ts_xxx_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xxzz[i] * gfe_0 + ts_xxx_xxzz[i] * rgc2_0;

        gr_xxx_xyyy[i] = 6.0 * ts_x_xyyy[i] * gfe2_0 + 6.0 * ts_xx_yyy[i] * gfe2_0 + 6.0 * ts_xx_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_xy[i] * gfe2_0 + 6.0 * ts_xxx_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_xyyy[i] * gfe_0 + ts_xxx_xyyy[i] * rgc2_0;

        gr_xxx_xyyz[i] = 6.0 * ts_x_xyyz[i] * gfe2_0 + 6.0 * ts_xx_yyz[i] * gfe2_0 + 6.0 * ts_xx_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xz[i] * gfe2_0 + 4.0 * ts_xxx_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xyyz[i] * gfe_0 + ts_xxx_xyyz[i] * rgc2_0;

        gr_xxx_xyzz[i] = 6.0 * ts_x_xyzz[i] * gfe2_0 + 6.0 * ts_xx_yzz[i] * gfe2_0 + 6.0 * ts_xx_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_xy[i] * gfe2_0 + 4.0 * ts_xxx_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xyzz[i] * gfe_0 + ts_xxx_xyzz[i] * rgc2_0;

        gr_xxx_xzzz[i] = 6.0 * ts_x_xzzz[i] * gfe2_0 + 6.0 * ts_xx_zzz[i] * gfe2_0 + 6.0 * ts_xx_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_xz[i] * gfe2_0 + 6.0 * ts_xxx_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xzzz[i] * gfe_0 + ts_xxx_xzzz[i] * rgc2_0;

        gr_xxx_yyyy[i] = 6.0 * ts_x_yyyy[i] * gfe2_0 + 6.0 * ts_xx_yyyy[i] * gfe_0 * gc_x[i] + 12.0 * ts_xxx_yy[i] * gfe2_0 + 8.0 * ts_xxx_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_yyyy[i] * gfe_0 + ts_xxx_yyyy[i] * rgc2_0;

        gr_xxx_yyyz[i] = 6.0 * ts_x_yyyz[i] * gfe2_0 + 6.0 * ts_xx_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_yz[i] * gfe2_0 + 6.0 * ts_xxx_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_yyyz[i] * gfe_0 + ts_xxx_yyyz[i] * rgc2_0;

        gr_xxx_yyzz[i] = 6.0 * ts_x_yyzz[i] * gfe2_0 + 6.0 * ts_xx_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zz[i] * gfe2_0 + 4.0 * ts_xxx_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_yy[i] * gfe2_0 + 4.0 * ts_xxx_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_yyzz[i] * gfe_0 + ts_xxx_yyzz[i] * rgc2_0;

        gr_xxx_yzzz[i] = 6.0 * ts_x_yzzz[i] * gfe2_0 + 6.0 * ts_xx_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxx_yz[i] * gfe2_0 + 6.0 * ts_xxx_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_yzzz[i] * gfe_0 + ts_xxx_yzzz[i] * rgc2_0;

        gr_xxx_zzzz[i] = 6.0 * ts_x_zzzz[i] * gfe2_0 + 6.0 * ts_xx_zzzz[i] * gfe_0 * gc_x[i] + 12.0 * ts_xxx_zz[i] * gfe2_0 + 8.0 * ts_xxx_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_zzzz[i] * gfe_0 + ts_xxx_zzzz[i] * rgc2_0;
    }

    // Set up 15-30 components of targeted buffer : FG

    auto gr_xxy_xxxx = pbuffer.data(idx_g_fg + 15);

    auto gr_xxy_xxxy = pbuffer.data(idx_g_fg + 16);

    auto gr_xxy_xxxz = pbuffer.data(idx_g_fg + 17);

    auto gr_xxy_xxyy = pbuffer.data(idx_g_fg + 18);

    auto gr_xxy_xxyz = pbuffer.data(idx_g_fg + 19);

    auto gr_xxy_xxzz = pbuffer.data(idx_g_fg + 20);

    auto gr_xxy_xyyy = pbuffer.data(idx_g_fg + 21);

    auto gr_xxy_xyyz = pbuffer.data(idx_g_fg + 22);

    auto gr_xxy_xyzz = pbuffer.data(idx_g_fg + 23);

    auto gr_xxy_xzzz = pbuffer.data(idx_g_fg + 24);

    auto gr_xxy_yyyy = pbuffer.data(idx_g_fg + 25);

    auto gr_xxy_yyyz = pbuffer.data(idx_g_fg + 26);

    auto gr_xxy_yyzz = pbuffer.data(idx_g_fg + 27);

    auto gr_xxy_yzzz = pbuffer.data(idx_g_fg + 28);

    auto gr_xxy_zzzz = pbuffer.data(idx_g_fg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xxxx, gr_xxy_xxxy, gr_xxy_xxxz, gr_xxy_xxyy, gr_xxy_xxyz, gr_xxy_xxzz, gr_xxy_xyyy, gr_xxy_xyyz, gr_xxy_xyzz, gr_xxy_xzzz, gr_xxy_yyyy, gr_xxy_yyyz, gr_xxy_yyzz, gr_xxy_yzzz, gr_xxy_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zz, ts_xxy_zzz, ts_xxy_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxy_xxxx[i] = 2.0 * ts_y_xxxx[i] * gfe2_0 + 16.0 * ts_xy_xxx[i] * gfe2_0 + 4.0 * ts_xy_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxy_xx[i] * gfe2_0 + 8.0 * ts_xxy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxy_xxxx[i] * gfe_0 + ts_xxy_xxxx[i] * rgc2_0;

        gr_xxy_xxxy[i] = 2.0 * ts_y_xxxy[i] * gfe2_0 + 12.0 * ts_xy_xxy[i] * gfe2_0 + 4.0 * ts_xy_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe2_0 + 2.0 * ts_xx_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_xy[i] * gfe2_0 + 6.0 * ts_xxy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_xxxy[i] * gfe_0 + ts_xxy_xxxy[i] * rgc2_0;

        gr_xxy_xxxz[i] = 2.0 * ts_y_xxxz[i] * gfe2_0 + 12.0 * ts_xy_xxz[i] * gfe2_0 + 4.0 * ts_xy_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_xz[i] * gfe2_0 + 6.0 * ts_xxy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xxxz[i] * gfe_0 + ts_xxy_xxxz[i] * rgc2_0;

        gr_xxy_xxyy[i] = 2.0 * ts_y_xxyy[i] * gfe2_0 + 8.0 * ts_xy_xyy[i] * gfe2_0 + 4.0 * ts_xy_xxyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_xxy[i] * gfe2_0 + 2.0 * ts_xx_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yy[i] * gfe2_0 + 4.0 * ts_xxy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xx[i] * gfe2_0 + 4.0 * ts_xxy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_xxyy[i] * gfe_0 + ts_xxy_xxyy[i] * rgc2_0;

        gr_xxy_xxyz[i] = 2.0 * ts_y_xxyz[i] * gfe2_0 + 8.0 * ts_xy_xyz[i] * gfe2_0 + 4.0 * ts_xy_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxz[i] * gfe2_0 + 2.0 * ts_xx_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yz[i] * gfe2_0 + 4.0 * ts_xxy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xxyz[i] * gfe_0 + ts_xxy_xxyz[i] * rgc2_0;

        gr_xxy_xxzz[i] = 2.0 * ts_y_xxzz[i] * gfe2_0 + 8.0 * ts_xy_xzz[i] * gfe2_0 + 4.0 * ts_xy_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_zz[i] * gfe2_0 + 4.0 * ts_xxy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xx[i] * gfe2_0 + 4.0 * ts_xxy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xxzz[i] * gfe_0 + ts_xxy_xxzz[i] * rgc2_0;

        gr_xxy_xyyy[i] = 2.0 * ts_y_xyyy[i] * gfe2_0 + 4.0 * ts_xy_yyy[i] * gfe2_0 + 4.0 * ts_xy_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_xyy[i] * gfe2_0 + 2.0 * ts_xx_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxy_xy[i] * gfe2_0 + 6.0 * ts_xxy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_xyyy[i] * gfe_0 + ts_xxy_xyyy[i] * rgc2_0;

        gr_xxy_xyyz[i] = 2.0 * ts_y_xyyz[i] * gfe2_0 + 4.0 * ts_xy_yyz[i] * gfe2_0 + 4.0 * ts_xy_xyyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_xyz[i] * gfe2_0 + 2.0 * ts_xx_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xz[i] * gfe2_0 + 4.0 * ts_xxy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xyyz[i] * gfe_0 + ts_xxy_xyyz[i] * rgc2_0;

        gr_xxy_xyzz[i] = 2.0 * ts_y_xyzz[i] * gfe2_0 + 4.0 * ts_xy_yzz[i] * gfe2_0 + 4.0 * ts_xy_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzz[i] * gfe2_0 + 2.0 * ts_xx_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xy[i] * gfe2_0 + 4.0 * ts_xxy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xyzz[i] * gfe_0 + ts_xxy_xyzz[i] * rgc2_0;

        gr_xxy_xzzz[i] = 2.0 * ts_y_xzzz[i] * gfe2_0 + 4.0 * ts_xy_zzz[i] * gfe2_0 + 4.0 * ts_xy_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxy_xz[i] * gfe2_0 + 6.0 * ts_xxy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xzzz[i] * gfe_0 + ts_xxy_xzzz[i] * rgc2_0;

        gr_xxy_yyyy[i] = 2.0 * ts_y_yyyy[i] * gfe2_0 + 4.0 * ts_xy_yyyy[i] * gfe_0 * gc_x[i] + 8.0 * ts_xx_yyy[i] * gfe2_0 + 2.0 * ts_xx_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxy_yy[i] * gfe2_0 + 8.0 * ts_xxy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_yyyy[i] * gfe_0 + ts_xxy_yyyy[i] * rgc2_0;

        gr_xxy_yyyz[i] = 2.0 * ts_y_yyyz[i] * gfe2_0 + 4.0 * ts_xy_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_yyz[i] * gfe2_0 + 2.0 * ts_xx_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_yz[i] * gfe2_0 + 6.0 * ts_xxy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_yyyz[i] * gfe_0 + ts_xxy_yyyz[i] * rgc2_0;

        gr_xxy_yyzz[i] = 2.0 * ts_y_yyzz[i] * gfe2_0 + 4.0 * ts_xy_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_yzz[i] * gfe2_0 + 2.0 * ts_xx_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_zz[i] * gfe2_0 + 4.0 * ts_xxy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yy[i] * gfe2_0 + 4.0 * ts_xxy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_yyzz[i] * gfe_0 + ts_xxy_yyzz[i] * rgc2_0;

        gr_xxy_yzzz[i] = 2.0 * ts_y_yzzz[i] * gfe2_0 + 4.0 * ts_xy_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzz[i] * gfe2_0 + 2.0 * ts_xx_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_yz[i] * gfe2_0 + 6.0 * ts_xxy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_yzzz[i] * gfe_0 + ts_xxy_yzzz[i] * rgc2_0;

        gr_xxy_zzzz[i] = 2.0 * ts_y_zzzz[i] * gfe2_0 + 4.0 * ts_xy_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxy_zz[i] * gfe2_0 + 8.0 * ts_xxy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_zzzz[i] * gfe_0 + ts_xxy_zzzz[i] * rgc2_0;
    }

    // Set up 30-45 components of targeted buffer : FG

    auto gr_xxz_xxxx = pbuffer.data(idx_g_fg + 30);

    auto gr_xxz_xxxy = pbuffer.data(idx_g_fg + 31);

    auto gr_xxz_xxxz = pbuffer.data(idx_g_fg + 32);

    auto gr_xxz_xxyy = pbuffer.data(idx_g_fg + 33);

    auto gr_xxz_xxyz = pbuffer.data(idx_g_fg + 34);

    auto gr_xxz_xxzz = pbuffer.data(idx_g_fg + 35);

    auto gr_xxz_xyyy = pbuffer.data(idx_g_fg + 36);

    auto gr_xxz_xyyz = pbuffer.data(idx_g_fg + 37);

    auto gr_xxz_xyzz = pbuffer.data(idx_g_fg + 38);

    auto gr_xxz_xzzz = pbuffer.data(idx_g_fg + 39);

    auto gr_xxz_yyyy = pbuffer.data(idx_g_fg + 40);

    auto gr_xxz_yyyz = pbuffer.data(idx_g_fg + 41);

    auto gr_xxz_yyzz = pbuffer.data(idx_g_fg + 42);

    auto gr_xxz_yzzz = pbuffer.data(idx_g_fg + 43);

    auto gr_xxz_zzzz = pbuffer.data(idx_g_fg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xxxx, gr_xxz_xxxy, gr_xxz_xxxz, gr_xxz_xxyy, gr_xxz_xxyz, gr_xxz_xxzz, gr_xxz_xyyy, gr_xxz_xyyz, gr_xxz_xyzz, gr_xxz_xzzz, gr_xxz_yyyy, gr_xxz_yyyz, gr_xxz_yyzz, gr_xxz_yzzz, gr_xxz_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zz, ts_xxz_zzz, ts_xxz_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe2_0 + 16.0 * ts_xz_xxx[i] * gfe2_0 + 4.0 * ts_xz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxz_xx[i] * gfe2_0 + 8.0 * ts_xxz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxz_xxxx[i] * gfe_0 + ts_xxz_xxxx[i] * rgc2_0;

        gr_xxz_xxxy[i] = 2.0 * ts_z_xxxy[i] * gfe2_0 + 12.0 * ts_xz_xxy[i] * gfe2_0 + 4.0 * ts_xz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxz_xy[i] * gfe2_0 + 6.0 * ts_xxz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_xxxy[i] * gfe_0 + ts_xxz_xxxy[i] * rgc2_0;

        gr_xxz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe2_0 + 12.0 * ts_xz_xxz[i] * gfe2_0 + 4.0 * ts_xz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe2_0 + 2.0 * ts_xx_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxz_xz[i] * gfe2_0 + 6.0 * ts_xxz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xxxz[i] * gfe_0 + ts_xxz_xxxz[i] * rgc2_0;

        gr_xxz_xxyy[i] = 2.0 * ts_z_xxyy[i] * gfe2_0 + 8.0 * ts_xz_xyy[i] * gfe2_0 + 4.0 * ts_xz_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yy[i] * gfe2_0 + 4.0 * ts_xxz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xx[i] * gfe2_0 + 4.0 * ts_xxz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_xxyy[i] * gfe_0 + ts_xxz_xxyy[i] * rgc2_0;

        gr_xxz_xxyz[i] = 2.0 * ts_z_xxyz[i] * gfe2_0 + 8.0 * ts_xz_xyz[i] * gfe2_0 + 4.0 * ts_xz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxy[i] * gfe2_0 + 2.0 * ts_xx_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yz[i] * gfe2_0 + 4.0 * ts_xxz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xxyz[i] * gfe_0 + ts_xxz_xxyz[i] * rgc2_0;

        gr_xxz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe2_0 + 8.0 * ts_xz_xzz[i] * gfe2_0 + 4.0 * ts_xz_xxzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_xxz[i] * gfe2_0 + 2.0 * ts_xx_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_zz[i] * gfe2_0 + 4.0 * ts_xxz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xx[i] * gfe2_0 + 4.0 * ts_xxz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xxzz[i] * gfe_0 + ts_xxz_xxzz[i] * rgc2_0;

        gr_xxz_xyyy[i] = 2.0 * ts_z_xyyy[i] * gfe2_0 + 4.0 * ts_xz_yyy[i] * gfe2_0 + 4.0 * ts_xz_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxz_xy[i] * gfe2_0 + 6.0 * ts_xxz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_xyyy[i] * gfe_0 + ts_xxz_xyyy[i] * rgc2_0;

        gr_xxz_xyyz[i] = 2.0 * ts_z_xyyz[i] * gfe2_0 + 4.0 * ts_xz_yyz[i] * gfe2_0 + 4.0 * ts_xz_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyy[i] * gfe2_0 + 2.0 * ts_xx_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xz[i] * gfe2_0 + 4.0 * ts_xxz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xyyz[i] * gfe_0 + ts_xxz_xyyz[i] * rgc2_0;

        gr_xxz_xyzz[i] = 2.0 * ts_z_xyzz[i] * gfe2_0 + 4.0 * ts_xz_yzz[i] * gfe2_0 + 4.0 * ts_xz_xyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_xyz[i] * gfe2_0 + 2.0 * ts_xx_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_xy[i] * gfe2_0 + 4.0 * ts_xxz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xyzz[i] * gfe_0 + ts_xxz_xyzz[i] * rgc2_0;

        gr_xxz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe2_0 + 4.0 * ts_xz_zzz[i] * gfe2_0 + 4.0 * ts_xz_xzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_xzz[i] * gfe2_0 + 2.0 * ts_xx_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxz_xz[i] * gfe2_0 + 6.0 * ts_xxz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xzzz[i] * gfe_0 + ts_xxz_xzzz[i] * rgc2_0;

        gr_xxz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe2_0 + 4.0 * ts_xz_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxz_yy[i] * gfe2_0 + 8.0 * ts_xxz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_yyyy[i] * gfe_0 + ts_xxz_yyyy[i] * rgc2_0;

        gr_xxz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe2_0 + 4.0 * ts_xz_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyy[i] * gfe2_0 + 2.0 * ts_xx_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxz_yz[i] * gfe2_0 + 6.0 * ts_xxz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_yyyz[i] * gfe_0 + ts_xxz_yyyz[i] * rgc2_0;

        gr_xxz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe2_0 + 4.0 * ts_xz_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_yyz[i] * gfe2_0 + 2.0 * ts_xx_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_zz[i] * gfe2_0 + 4.0 * ts_xxz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_yy[i] * gfe2_0 + 4.0 * ts_xxz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_yyzz[i] * gfe_0 + ts_xxz_yyzz[i] * rgc2_0;

        gr_xxz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe2_0 + 4.0 * ts_xz_yzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_yzz[i] * gfe2_0 + 2.0 * ts_xx_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxz_yz[i] * gfe2_0 + 6.0 * ts_xxz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_yzzz[i] * gfe_0 + ts_xxz_yzzz[i] * rgc2_0;

        gr_xxz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe2_0 + 4.0 * ts_xz_zzzz[i] * gfe_0 * gc_x[i] + 8.0 * ts_xx_zzz[i] * gfe2_0 + 2.0 * ts_xx_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxz_zz[i] * gfe2_0 + 8.0 * ts_xxz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_zzzz[i] * gfe_0 + ts_xxz_zzzz[i] * rgc2_0;
    }

    // Set up 45-60 components of targeted buffer : FG

    auto gr_xyy_xxxx = pbuffer.data(idx_g_fg + 45);

    auto gr_xyy_xxxy = pbuffer.data(idx_g_fg + 46);

    auto gr_xyy_xxxz = pbuffer.data(idx_g_fg + 47);

    auto gr_xyy_xxyy = pbuffer.data(idx_g_fg + 48);

    auto gr_xyy_xxyz = pbuffer.data(idx_g_fg + 49);

    auto gr_xyy_xxzz = pbuffer.data(idx_g_fg + 50);

    auto gr_xyy_xyyy = pbuffer.data(idx_g_fg + 51);

    auto gr_xyy_xyyz = pbuffer.data(idx_g_fg + 52);

    auto gr_xyy_xyzz = pbuffer.data(idx_g_fg + 53);

    auto gr_xyy_xzzz = pbuffer.data(idx_g_fg + 54);

    auto gr_xyy_yyyy = pbuffer.data(idx_g_fg + 55);

    auto gr_xyy_yyyz = pbuffer.data(idx_g_fg + 56);

    auto gr_xyy_yyzz = pbuffer.data(idx_g_fg + 57);

    auto gr_xyy_yzzz = pbuffer.data(idx_g_fg + 58);

    auto gr_xyy_zzzz = pbuffer.data(idx_g_fg + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xxxx, gr_xyy_xxxy, gr_xyy_xxxz, gr_xyy_xxyy, gr_xyy_xxyz, gr_xyy_xxzz, gr_xyy_xyyy, gr_xyy_xyyz, gr_xyy_xyzz, gr_xyy_xzzz, gr_xyy_yyyy, gr_xyy_yyyz, gr_xyy_yyzz, gr_xyy_yzzz, gr_xyy_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zz, ts_xyy_zzz, ts_xyy_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyy_xxxx[i] = 8.0 * ts_yy_xxx[i] * gfe2_0 + 2.0 * ts_yy_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxx[i] * gfe2_0 + 4.0 * ts_xy_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_xyy_xx[i] * gfe2_0 + 8.0 * ts_xyy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyy_xxxx[i] * gfe_0 + ts_xyy_xxxx[i] * rgc2_0;

        gr_xyy_xxxy[i] = 6.0 * ts_yy_xxy[i] * gfe2_0 + 2.0 * ts_yy_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxy[i] * gfe2_0 + 4.0 * ts_xy_xxx[i] * gfe2_0 + 4.0 * ts_xy_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_xy[i] * gfe2_0 + 6.0 * ts_xyy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_xxxy[i] * gfe_0 + ts_xyy_xxxy[i] * rgc2_0;

        gr_xyy_xxxz[i] = 6.0 * ts_yy_xxz[i] * gfe2_0 + 2.0 * ts_yy_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxz[i] * gfe2_0 + 4.0 * ts_xy_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_xz[i] * gfe2_0 + 6.0 * ts_xyy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xxxz[i] * gfe_0 + ts_xyy_xxxz[i] * rgc2_0;

        gr_xyy_xxyy[i] = 4.0 * ts_yy_xyy[i] * gfe2_0 + 2.0 * ts_yy_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxyy[i] * gfe2_0 + 8.0 * ts_xy_xxy[i] * gfe2_0 + 4.0 * ts_xy_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yy[i] * gfe2_0 + 4.0 * ts_xyy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xx[i] * gfe2_0 + 4.0 * ts_xyy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_xxyy[i] * gfe_0 + ts_xyy_xxyy[i] * rgc2_0;

        gr_xyy_xxyz[i] = 4.0 * ts_yy_xyz[i] * gfe2_0 + 2.0 * ts_yy_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxyz[i] * gfe2_0 + 4.0 * ts_xy_xxz[i] * gfe2_0 + 4.0 * ts_xy_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yz[i] * gfe2_0 + 4.0 * ts_xyy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xxyz[i] * gfe_0 + ts_xyy_xxyz[i] * rgc2_0;

        gr_xyy_xxzz[i] = 4.0 * ts_yy_xzz[i] * gfe2_0 + 2.0 * ts_yy_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxzz[i] * gfe2_0 + 4.0 * ts_xy_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_zz[i] * gfe2_0 + 4.0 * ts_xyy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xx[i] * gfe2_0 + 4.0 * ts_xyy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xxzz[i] * gfe_0 + ts_xyy_xxzz[i] * rgc2_0;

        gr_xyy_xyyy[i] = 2.0 * ts_yy_yyy[i] * gfe2_0 + 2.0 * ts_yy_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyyy[i] * gfe2_0 + 12.0 * ts_xy_xyy[i] * gfe2_0 + 4.0 * ts_xy_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyy_xy[i] * gfe2_0 + 6.0 * ts_xyy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_xyyy[i] * gfe_0 + ts_xyy_xyyy[i] * rgc2_0;

        gr_xyy_xyyz[i] = 2.0 * ts_yy_yyz[i] * gfe2_0 + 2.0 * ts_yy_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyyz[i] * gfe2_0 + 8.0 * ts_xy_xyz[i] * gfe2_0 + 4.0 * ts_xy_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xz[i] * gfe2_0 + 4.0 * ts_xyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xyyz[i] * gfe_0 + ts_xyy_xyyz[i] * rgc2_0;

        gr_xyy_xyzz[i] = 2.0 * ts_yy_yzz[i] * gfe2_0 + 2.0 * ts_yy_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyzz[i] * gfe2_0 + 4.0 * ts_xy_xzz[i] * gfe2_0 + 4.0 * ts_xy_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xy[i] * gfe2_0 + 4.0 * ts_xyy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xyzz[i] * gfe_0 + ts_xyy_xyzz[i] * rgc2_0;

        gr_xyy_xzzz[i] = 2.0 * ts_yy_zzz[i] * gfe2_0 + 2.0 * ts_yy_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzzz[i] * gfe2_0 + 4.0 * ts_xy_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyy_xz[i] * gfe2_0 + 6.0 * ts_xyy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xzzz[i] * gfe_0 + ts_xyy_xzzz[i] * rgc2_0;

        gr_xyy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyyy[i] * gfe2_0 + 16.0 * ts_xy_yyy[i] * gfe2_0 + 4.0 * ts_xy_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_xyy_yy[i] * gfe2_0 + 8.0 * ts_xyy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_yyyy[i] * gfe_0 + ts_xyy_yyyy[i] * rgc2_0;

        gr_xyy_yyyz[i] = 2.0 * ts_yy_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyyz[i] * gfe2_0 + 12.0 * ts_xy_yyz[i] * gfe2_0 + 4.0 * ts_xy_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_yz[i] * gfe2_0 + 6.0 * ts_xyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_yyyz[i] * gfe_0 + ts_xyy_yyyz[i] * rgc2_0;

        gr_xyy_yyzz[i] = 2.0 * ts_yy_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyzz[i] * gfe2_0 + 8.0 * ts_xy_yzz[i] * gfe2_0 + 4.0 * ts_xy_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_zz[i] * gfe2_0 + 4.0 * ts_xyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yy[i] * gfe2_0 + 4.0 * ts_xyy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_yyzz[i] * gfe_0 + ts_xyy_yyzz[i] * rgc2_0;

        gr_xyy_yzzz[i] = 2.0 * ts_yy_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yzzz[i] * gfe2_0 + 4.0 * ts_xy_zzz[i] * gfe2_0 + 4.0 * ts_xy_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_yz[i] * gfe2_0 + 6.0 * ts_xyy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_yzzz[i] * gfe_0 + ts_xyy_yzzz[i] * rgc2_0;

        gr_xyy_zzzz[i] = 2.0 * ts_yy_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzzz[i] * gfe2_0 + 4.0 * ts_xy_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_xyy_zz[i] * gfe2_0 + 8.0 * ts_xyy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_zzzz[i] * gfe_0 + ts_xyy_zzzz[i] * rgc2_0;
    }

    // Set up 60-75 components of targeted buffer : FG

    auto gr_xyz_xxxx = pbuffer.data(idx_g_fg + 60);

    auto gr_xyz_xxxy = pbuffer.data(idx_g_fg + 61);

    auto gr_xyz_xxxz = pbuffer.data(idx_g_fg + 62);

    auto gr_xyz_xxyy = pbuffer.data(idx_g_fg + 63);

    auto gr_xyz_xxyz = pbuffer.data(idx_g_fg + 64);

    auto gr_xyz_xxzz = pbuffer.data(idx_g_fg + 65);

    auto gr_xyz_xyyy = pbuffer.data(idx_g_fg + 66);

    auto gr_xyz_xyyz = pbuffer.data(idx_g_fg + 67);

    auto gr_xyz_xyzz = pbuffer.data(idx_g_fg + 68);

    auto gr_xyz_xzzz = pbuffer.data(idx_g_fg + 69);

    auto gr_xyz_yyyy = pbuffer.data(idx_g_fg + 70);

    auto gr_xyz_yyyz = pbuffer.data(idx_g_fg + 71);

    auto gr_xyz_yyzz = pbuffer.data(idx_g_fg + 72);

    auto gr_xyz_yzzz = pbuffer.data(idx_g_fg + 73);

    auto gr_xyz_zzzz = pbuffer.data(idx_g_fg + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xxxx, gr_xyz_xxxy, gr_xyz_xxxz, gr_xyz_xxyy, gr_xyz_xxyz, gr_xyz_xxzz, gr_xyz_xyyy, gr_xyz_xyyz, gr_xyz_xyzz, gr_xyz_xzzz, gr_xyz_yyyy, gr_xyz_yyyz, gr_xyz_yyzz, gr_xyz_yzzz, gr_xyz_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zz, ts_xyz_zzz, ts_xyz_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyz_xxxx[i] = 8.0 * ts_yz_xxx[i] * gfe2_0 + 2.0 * ts_yz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyz_xx[i] * gfe2_0 + 8.0 * ts_xyz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyz_xxxx[i] * gfe_0 + ts_xyz_xxxx[i] * rgc2_0;

        gr_xyz_xxxy[i] = 6.0 * ts_yz_xxy[i] * gfe2_0 + 2.0 * ts_yz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxx[i] * gfe2_0 + 2.0 * ts_xz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyz_xy[i] * gfe2_0 + 6.0 * ts_xyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_xxxy[i] * gfe_0 + ts_xyz_xxxy[i] * rgc2_0;

        gr_xyz_xxxz[i] = 6.0 * ts_yz_xxz[i] * gfe2_0 + 2.0 * ts_yz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxx[i] * gfe2_0 + 2.0 * ts_xy_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyz_xz[i] * gfe2_0 + 6.0 * ts_xyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xxxz[i] * gfe_0 + ts_xyz_xxxz[i] * rgc2_0;

        gr_xyz_xxyy[i] = 4.0 * ts_yz_xyy[i] * gfe2_0 + 2.0 * ts_yz_xxyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xz_xxy[i] * gfe2_0 + 2.0 * ts_xz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yy[i] * gfe2_0 + 4.0 * ts_xyz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xx[i] * gfe2_0 + 4.0 * ts_xyz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_xxyy[i] * gfe_0 + ts_xyz_xxyy[i] * rgc2_0;

        gr_xyz_xxyz[i] = 4.0 * ts_yz_xyz[i] * gfe2_0 + 2.0 * ts_yz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxz[i] * gfe2_0 + 2.0 * ts_xz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxy[i] * gfe2_0 + 2.0 * ts_xy_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yz[i] * gfe2_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xxyz[i] * gfe_0 + ts_xyz_xxyz[i] * rgc2_0;

        gr_xyz_xxzz[i] = 4.0 * ts_yz_xzz[i] * gfe2_0 + 2.0 * ts_yz_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xy_xxz[i] * gfe2_0 + 2.0 * ts_xy_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_zz[i] * gfe2_0 + 4.0 * ts_xyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xx[i] * gfe2_0 + 4.0 * ts_xyz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xxzz[i] * gfe_0 + ts_xyz_xxzz[i] * rgc2_0;

        gr_xyz_xyyy[i] = 2.0 * ts_yz_yyy[i] * gfe2_0 + 2.0 * ts_yz_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xyy[i] * gfe2_0 + 2.0 * ts_xz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyz_xy[i] * gfe2_0 + 6.0 * ts_xyz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_xyyy[i] * gfe_0 + ts_xyz_xyyy[i] * rgc2_0;

        gr_xyz_xyyz[i] = 2.0 * ts_yz_yyz[i] * gfe2_0 + 2.0 * ts_yz_xyyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xz_xyz[i] * gfe2_0 + 2.0 * ts_xz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyy[i] * gfe2_0 + 2.0 * ts_xy_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xz[i] * gfe2_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xyyz[i] * gfe_0 + ts_xyz_xyyz[i] * rgc2_0;

        gr_xyz_xyzz[i] = 2.0 * ts_yz_yzz[i] * gfe2_0 + 2.0 * ts_yz_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xzz[i] * gfe2_0 + 2.0 * ts_xz_xyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xy_xyz[i] * gfe2_0 + 2.0 * ts_xy_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_xy[i] * gfe2_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xyzz[i] * gfe_0 + ts_xyz_xyzz[i] * rgc2_0;

        gr_xyz_xzzz[i] = 2.0 * ts_yz_zzz[i] * gfe2_0 + 2.0 * ts_yz_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_xzz[i] * gfe2_0 + 2.0 * ts_xy_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyz_xz[i] * gfe2_0 + 6.0 * ts_xyz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xzzz[i] * gfe_0 + ts_xyz_xzzz[i] * rgc2_0;

        gr_xyz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe_0 * gc_x[i] + 8.0 * ts_xz_yyy[i] * gfe2_0 + 2.0 * ts_xz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyz_yy[i] * gfe2_0 + 8.0 * ts_xyz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_yyyy[i] * gfe_0 + ts_xyz_yyyy[i] * rgc2_0;

        gr_xyz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yyz[i] * gfe2_0 + 2.0 * ts_xz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyy[i] * gfe2_0 + 2.0 * ts_xy_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyz_yz[i] * gfe2_0 + 6.0 * ts_xyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_yyyz[i] * gfe_0 + ts_xyz_yyyz[i] * rgc2_0;

        gr_xyz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xz_yzz[i] * gfe2_0 + 2.0 * ts_xz_yyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xy_yyz[i] * gfe2_0 + 2.0 * ts_xy_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_zz[i] * gfe2_0 + 4.0 * ts_xyz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_yy[i] * gfe2_0 + 4.0 * ts_xyz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_yyzz[i] * gfe_0 + ts_xyz_yyzz[i] * rgc2_0;

        gr_xyz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zzz[i] * gfe2_0 + 2.0 * ts_xz_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_yzz[i] * gfe2_0 + 2.0 * ts_xy_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyz_yz[i] * gfe2_0 + 6.0 * ts_xyz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_yzzz[i] * gfe_0 + ts_xyz_yzzz[i] * rgc2_0;

        gr_xyz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zzzz[i] * gfe_0 * gc_y[i] + 8.0 * ts_xy_zzz[i] * gfe2_0 + 2.0 * ts_xy_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyz_zz[i] * gfe2_0 + 8.0 * ts_xyz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_zzzz[i] * gfe_0 + ts_xyz_zzzz[i] * rgc2_0;
    }

    // Set up 75-90 components of targeted buffer : FG

    auto gr_xzz_xxxx = pbuffer.data(idx_g_fg + 75);

    auto gr_xzz_xxxy = pbuffer.data(idx_g_fg + 76);

    auto gr_xzz_xxxz = pbuffer.data(idx_g_fg + 77);

    auto gr_xzz_xxyy = pbuffer.data(idx_g_fg + 78);

    auto gr_xzz_xxyz = pbuffer.data(idx_g_fg + 79);

    auto gr_xzz_xxzz = pbuffer.data(idx_g_fg + 80);

    auto gr_xzz_xyyy = pbuffer.data(idx_g_fg + 81);

    auto gr_xzz_xyyz = pbuffer.data(idx_g_fg + 82);

    auto gr_xzz_xyzz = pbuffer.data(idx_g_fg + 83);

    auto gr_xzz_xzzz = pbuffer.data(idx_g_fg + 84);

    auto gr_xzz_yyyy = pbuffer.data(idx_g_fg + 85);

    auto gr_xzz_yyyz = pbuffer.data(idx_g_fg + 86);

    auto gr_xzz_yyzz = pbuffer.data(idx_g_fg + 87);

    auto gr_xzz_yzzz = pbuffer.data(idx_g_fg + 88);

    auto gr_xzz_zzzz = pbuffer.data(idx_g_fg + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xxxx, gr_xzz_xxxy, gr_xzz_xxxz, gr_xzz_xxyy, gr_xzz_xxyz, gr_xzz_xxzz, gr_xzz_xyyy, gr_xzz_xyyz, gr_xzz_xyzz, gr_xzz_xzzz, gr_xzz_yyyy, gr_xzz_yyyz, gr_xzz_yyzz, gr_xzz_yzzz, gr_xzz_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zz, ts_xzz_zzz, ts_xzz_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xzz_xxxx[i] = 8.0 * ts_zz_xxx[i] * gfe2_0 + 2.0 * ts_zz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxx[i] * gfe2_0 + 4.0 * ts_xz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xzz_xx[i] * gfe2_0 + 8.0 * ts_xzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzz_xxxx[i] * gfe_0 + ts_xzz_xxxx[i] * rgc2_0;

        gr_xzz_xxxy[i] = 6.0 * ts_zz_xxy[i] * gfe2_0 + 2.0 * ts_zz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxy[i] * gfe2_0 + 4.0 * ts_xz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzz_xy[i] * gfe2_0 + 6.0 * ts_xzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_xxxy[i] * gfe_0 + ts_xzz_xxxy[i] * rgc2_0;

        gr_xzz_xxxz[i] = 6.0 * ts_zz_xxz[i] * gfe2_0 + 2.0 * ts_zz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxz[i] * gfe2_0 + 4.0 * ts_xz_xxx[i] * gfe2_0 + 4.0 * ts_xz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzz_xz[i] * gfe2_0 + 6.0 * ts_xzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xxxz[i] * gfe_0 + ts_xzz_xxxz[i] * rgc2_0;

        gr_xzz_xxyy[i] = 4.0 * ts_zz_xyy[i] * gfe2_0 + 2.0 * ts_zz_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxyy[i] * gfe2_0 + 4.0 * ts_xz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yy[i] * gfe2_0 + 4.0 * ts_xzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xx[i] * gfe2_0 + 4.0 * ts_xzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_xxyy[i] * gfe_0 + ts_xzz_xxyy[i] * rgc2_0;

        gr_xzz_xxyz[i] = 4.0 * ts_zz_xyz[i] * gfe2_0 + 2.0 * ts_zz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxyz[i] * gfe2_0 + 4.0 * ts_xz_xxy[i] * gfe2_0 + 4.0 * ts_xz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yz[i] * gfe2_0 + 4.0 * ts_xzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xxyz[i] * gfe_0 + ts_xzz_xxyz[i] * rgc2_0;

        gr_xzz_xxzz[i] = 4.0 * ts_zz_xzz[i] * gfe2_0 + 2.0 * ts_zz_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxzz[i] * gfe2_0 + 8.0 * ts_xz_xxz[i] * gfe2_0 + 4.0 * ts_xz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_zz[i] * gfe2_0 + 4.0 * ts_xzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xx[i] * gfe2_0 + 4.0 * ts_xzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xxzz[i] * gfe_0 + ts_xzz_xxzz[i] * rgc2_0;

        gr_xzz_xyyy[i] = 2.0 * ts_zz_yyy[i] * gfe2_0 + 2.0 * ts_zz_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyyy[i] * gfe2_0 + 4.0 * ts_xz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzz_xy[i] * gfe2_0 + 6.0 * ts_xzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_xyyy[i] * gfe_0 + ts_xzz_xyyy[i] * rgc2_0;

        gr_xzz_xyyz[i] = 2.0 * ts_zz_yyz[i] * gfe2_0 + 2.0 * ts_zz_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyyz[i] * gfe2_0 + 4.0 * ts_xz_xyy[i] * gfe2_0 + 4.0 * ts_xz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xz[i] * gfe2_0 + 4.0 * ts_xzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xyyz[i] * gfe_0 + ts_xzz_xyyz[i] * rgc2_0;

        gr_xzz_xyzz[i] = 2.0 * ts_zz_yzz[i] * gfe2_0 + 2.0 * ts_zz_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyzz[i] * gfe2_0 + 8.0 * ts_xz_xyz[i] * gfe2_0 + 4.0 * ts_xz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_xy[i] * gfe2_0 + 4.0 * ts_xzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xyzz[i] * gfe_0 + ts_xzz_xyzz[i] * rgc2_0;

        gr_xzz_xzzz[i] = 2.0 * ts_zz_zzz[i] * gfe2_0 + 2.0 * ts_zz_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzzz[i] * gfe2_0 + 12.0 * ts_xz_xzz[i] * gfe2_0 + 4.0 * ts_xz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzz_xz[i] * gfe2_0 + 6.0 * ts_xzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xzzz[i] * gfe_0 + ts_xzz_xzzz[i] * rgc2_0;

        gr_xzz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyyy[i] * gfe2_0 + 4.0 * ts_xz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xzz_yy[i] * gfe2_0 + 8.0 * ts_xzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_yyyy[i] * gfe_0 + ts_xzz_yyyy[i] * rgc2_0;

        gr_xzz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyyz[i] * gfe2_0 + 4.0 * ts_xz_yyy[i] * gfe2_0 + 4.0 * ts_xz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzz_yz[i] * gfe2_0 + 6.0 * ts_xzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_yyyz[i] * gfe_0 + ts_xzz_yyyz[i] * rgc2_0;

        gr_xzz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyzz[i] * gfe2_0 + 8.0 * ts_xz_yyz[i] * gfe2_0 + 4.0 * ts_xz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_zz[i] * gfe2_0 + 4.0 * ts_xzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_yy[i] * gfe2_0 + 4.0 * ts_xzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_yyzz[i] * gfe_0 + ts_xzz_yyzz[i] * rgc2_0;

        gr_xzz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yzzz[i] * gfe2_0 + 12.0 * ts_xz_yzz[i] * gfe2_0 + 4.0 * ts_xz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xzz_yz[i] * gfe2_0 + 6.0 * ts_xzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_yzzz[i] * gfe_0 + ts_xzz_yzzz[i] * rgc2_0;

        gr_xzz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzzz[i] * gfe2_0 + 16.0 * ts_xz_zzz[i] * gfe2_0 + 4.0 * ts_xz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xzz_zz[i] * gfe2_0 + 8.0 * ts_xzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_zzzz[i] * gfe_0 + ts_xzz_zzzz[i] * rgc2_0;
    }

    // Set up 90-105 components of targeted buffer : FG

    auto gr_yyy_xxxx = pbuffer.data(idx_g_fg + 90);

    auto gr_yyy_xxxy = pbuffer.data(idx_g_fg + 91);

    auto gr_yyy_xxxz = pbuffer.data(idx_g_fg + 92);

    auto gr_yyy_xxyy = pbuffer.data(idx_g_fg + 93);

    auto gr_yyy_xxyz = pbuffer.data(idx_g_fg + 94);

    auto gr_yyy_xxzz = pbuffer.data(idx_g_fg + 95);

    auto gr_yyy_xyyy = pbuffer.data(idx_g_fg + 96);

    auto gr_yyy_xyyz = pbuffer.data(idx_g_fg + 97);

    auto gr_yyy_xyzz = pbuffer.data(idx_g_fg + 98);

    auto gr_yyy_xzzz = pbuffer.data(idx_g_fg + 99);

    auto gr_yyy_yyyy = pbuffer.data(idx_g_fg + 100);

    auto gr_yyy_yyyz = pbuffer.data(idx_g_fg + 101);

    auto gr_yyy_yyzz = pbuffer.data(idx_g_fg + 102);

    auto gr_yyy_yzzz = pbuffer.data(idx_g_fg + 103);

    auto gr_yyy_zzzz = pbuffer.data(idx_g_fg + 104);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xxxx, gr_yyy_xxxy, gr_yyy_xxxz, gr_yyy_xxyy, gr_yyy_xxyz, gr_yyy_xxzz, gr_yyy_xyyy, gr_yyy_xyyz, gr_yyy_xyzz, gr_yyy_xzzz, gr_yyy_yyyy, gr_yyy_yyyz, gr_yyy_yyzz, gr_yyy_yzzz, gr_yyy_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyy_xxxx[i] = 6.0 * ts_y_xxxx[i] * gfe2_0 + 6.0 * ts_yy_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_yyy_xx[i] * gfe2_0 + 8.0 * ts_yyy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyy_xxxx[i] * gfe_0 + ts_yyy_xxxx[i] * rgc2_0;

        gr_yyy_xxxy[i] = 6.0 * ts_y_xxxy[i] * gfe2_0 + 6.0 * ts_yy_xxx[i] * gfe2_0 + 6.0 * ts_yy_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_xy[i] * gfe2_0 + 6.0 * ts_yyy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_xxxy[i] * gfe_0 + ts_yyy_xxxy[i] * rgc2_0;

        gr_yyy_xxxz[i] = 6.0 * ts_y_xxxz[i] * gfe2_0 + 6.0 * ts_yy_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_xz[i] * gfe2_0 + 6.0 * ts_yyy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xxxz[i] * gfe_0 + ts_yyy_xxxz[i] * rgc2_0;

        gr_yyy_xxyy[i] = 6.0 * ts_y_xxyy[i] * gfe2_0 + 12.0 * ts_yy_xxy[i] * gfe2_0 + 6.0 * ts_yy_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yy[i] * gfe2_0 + 4.0 * ts_yyy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xx[i] * gfe2_0 + 4.0 * ts_yyy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_xxyy[i] * gfe_0 + ts_yyy_xxyy[i] * rgc2_0;

        gr_yyy_xxyz[i] = 6.0 * ts_y_xxyz[i] * gfe2_0 + 6.0 * ts_yy_xxz[i] * gfe2_0 + 6.0 * ts_yy_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yz[i] * gfe2_0 + 4.0 * ts_yyy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xxyz[i] * gfe_0 + ts_yyy_xxyz[i] * rgc2_0;

        gr_yyy_xxzz[i] = 6.0 * ts_y_xxzz[i] * gfe2_0 + 6.0 * ts_yy_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_zz[i] * gfe2_0 + 4.0 * ts_yyy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xx[i] * gfe2_0 + 4.0 * ts_yyy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xxzz[i] * gfe_0 + ts_yyy_xxzz[i] * rgc2_0;

        gr_yyy_xyyy[i] = 6.0 * ts_y_xyyy[i] * gfe2_0 + 18.0 * ts_yy_xyy[i] * gfe2_0 + 6.0 * ts_yy_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyy_xy[i] * gfe2_0 + 6.0 * ts_yyy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_xyyy[i] * gfe_0 + ts_yyy_xyyy[i] * rgc2_0;

        gr_yyy_xyyz[i] = 6.0 * ts_y_xyyz[i] * gfe2_0 + 12.0 * ts_yy_xyz[i] * gfe2_0 + 6.0 * ts_yy_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xz[i] * gfe2_0 + 4.0 * ts_yyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xyyz[i] * gfe_0 + ts_yyy_xyyz[i] * rgc2_0;

        gr_yyy_xyzz[i] = 6.0 * ts_y_xyzz[i] * gfe2_0 + 6.0 * ts_yy_xzz[i] * gfe2_0 + 6.0 * ts_yy_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xy[i] * gfe2_0 + 4.0 * ts_yyy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xyzz[i] * gfe_0 + ts_yyy_xyzz[i] * rgc2_0;

        gr_yyy_xzzz[i] = 6.0 * ts_y_xzzz[i] * gfe2_0 + 6.0 * ts_yy_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyy_xz[i] * gfe2_0 + 6.0 * ts_yyy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xzzz[i] * gfe_0 + ts_yyy_xzzz[i] * rgc2_0;

        gr_yyy_yyyy[i] = 6.0 * ts_y_yyyy[i] * gfe2_0 + 24.0 * ts_yy_yyy[i] * gfe2_0 + 6.0 * ts_yy_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_yyy_yy[i] * gfe2_0 + 8.0 * ts_yyy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_yyyy[i] * gfe_0 + ts_yyy_yyyy[i] * rgc2_0;

        gr_yyy_yyyz[i] = 6.0 * ts_y_yyyz[i] * gfe2_0 + 18.0 * ts_yy_yyz[i] * gfe2_0 + 6.0 * ts_yy_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_yz[i] * gfe2_0 + 6.0 * ts_yyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_yyyz[i] * gfe_0 + ts_yyy_yyyz[i] * rgc2_0;

        gr_yyy_yyzz[i] = 6.0 * ts_y_yyzz[i] * gfe2_0 + 12.0 * ts_yy_yzz[i] * gfe2_0 + 6.0 * ts_yy_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_zz[i] * gfe2_0 + 4.0 * ts_yyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yy[i] * gfe2_0 + 4.0 * ts_yyy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_yyzz[i] * gfe_0 + ts_yyy_yyzz[i] * rgc2_0;

        gr_yyy_yzzz[i] = 6.0 * ts_y_yzzz[i] * gfe2_0 + 6.0 * ts_yy_zzz[i] * gfe2_0 + 6.0 * ts_yy_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_yz[i] * gfe2_0 + 6.0 * ts_yyy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_yzzz[i] * gfe_0 + ts_yyy_yzzz[i] * rgc2_0;

        gr_yyy_zzzz[i] = 6.0 * ts_y_zzzz[i] * gfe2_0 + 6.0 * ts_yy_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_yyy_zz[i] * gfe2_0 + 8.0 * ts_yyy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_zzzz[i] * gfe_0 + ts_yyy_zzzz[i] * rgc2_0;
    }

    // Set up 105-120 components of targeted buffer : FG

    auto gr_yyz_xxxx = pbuffer.data(idx_g_fg + 105);

    auto gr_yyz_xxxy = pbuffer.data(idx_g_fg + 106);

    auto gr_yyz_xxxz = pbuffer.data(idx_g_fg + 107);

    auto gr_yyz_xxyy = pbuffer.data(idx_g_fg + 108);

    auto gr_yyz_xxyz = pbuffer.data(idx_g_fg + 109);

    auto gr_yyz_xxzz = pbuffer.data(idx_g_fg + 110);

    auto gr_yyz_xyyy = pbuffer.data(idx_g_fg + 111);

    auto gr_yyz_xyyz = pbuffer.data(idx_g_fg + 112);

    auto gr_yyz_xyzz = pbuffer.data(idx_g_fg + 113);

    auto gr_yyz_xzzz = pbuffer.data(idx_g_fg + 114);

    auto gr_yyz_yyyy = pbuffer.data(idx_g_fg + 115);

    auto gr_yyz_yyyz = pbuffer.data(idx_g_fg + 116);

    auto gr_yyz_yyzz = pbuffer.data(idx_g_fg + 117);

    auto gr_yyz_yzzz = pbuffer.data(idx_g_fg + 118);

    auto gr_yyz_zzzz = pbuffer.data(idx_g_fg + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xxxx, gr_yyz_xxxy, gr_yyz_xxxz, gr_yyz_xxyy, gr_yyz_xxyz, gr_yyz_xxzz, gr_yyz_xyyy, gr_yyz_xyyz, gr_yyz_xyzz, gr_yyz_xzzz, gr_yyz_yyyy, gr_yyz_yyyz, gr_yyz_yyzz, gr_yyz_yzzz, gr_yyz_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zz, ts_yyz_zzz, ts_yyz_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe2_0 + 4.0 * ts_yz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyz_xx[i] * gfe2_0 + 8.0 * ts_yyz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyz_xxxx[i] * gfe_0 + ts_yyz_xxxx[i] * rgc2_0;

        gr_yyz_xxxy[i] = 2.0 * ts_z_xxxy[i] * gfe2_0 + 4.0 * ts_yz_xxx[i] * gfe2_0 + 4.0 * ts_yz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyz_xy[i] * gfe2_0 + 6.0 * ts_yyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_xxxy[i] * gfe_0 + ts_yyz_xxxy[i] * rgc2_0;

        gr_yyz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe2_0 + 4.0 * ts_yz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxx[i] * gfe2_0 + 2.0 * ts_yy_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyz_xz[i] * gfe2_0 + 6.0 * ts_yyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xxxz[i] * gfe_0 + ts_yyz_xxxz[i] * rgc2_0;

        gr_yyz_xxyy[i] = 2.0 * ts_z_xxyy[i] * gfe2_0 + 8.0 * ts_yz_xxy[i] * gfe2_0 + 4.0 * ts_yz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yy[i] * gfe2_0 + 4.0 * ts_yyz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xx[i] * gfe2_0 + 4.0 * ts_yyz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_xxyy[i] * gfe_0 + ts_yyz_xxyy[i] * rgc2_0;

        gr_yyz_xxyz[i] = 2.0 * ts_z_xxyz[i] * gfe2_0 + 4.0 * ts_yz_xxz[i] * gfe2_0 + 4.0 * ts_yz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxy[i] * gfe2_0 + 2.0 * ts_yy_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yz[i] * gfe2_0 + 4.0 * ts_yyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xxyz[i] * gfe_0 + ts_yyz_xxyz[i] * rgc2_0;

        gr_yyz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe2_0 + 4.0 * ts_yz_xxzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yy_xxz[i] * gfe2_0 + 2.0 * ts_yy_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_zz[i] * gfe2_0 + 4.0 * ts_yyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xx[i] * gfe2_0 + 4.0 * ts_yyz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xxzz[i] * gfe_0 + ts_yyz_xxzz[i] * rgc2_0;

        gr_yyz_xyyy[i] = 2.0 * ts_z_xyyy[i] * gfe2_0 + 12.0 * ts_yz_xyy[i] * gfe2_0 + 4.0 * ts_yz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyz_xy[i] * gfe2_0 + 6.0 * ts_yyz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_xyyy[i] * gfe_0 + ts_yyz_xyyy[i] * rgc2_0;

        gr_yyz_xyyz[i] = 2.0 * ts_z_xyyz[i] * gfe2_0 + 8.0 * ts_yz_xyz[i] * gfe2_0 + 4.0 * ts_yz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyy[i] * gfe2_0 + 2.0 * ts_yy_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xz[i] * gfe2_0 + 4.0 * ts_yyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xyyz[i] * gfe_0 + ts_yyz_xyyz[i] * rgc2_0;

        gr_yyz_xyzz[i] = 2.0 * ts_z_xyzz[i] * gfe2_0 + 4.0 * ts_yz_xzz[i] * gfe2_0 + 4.0 * ts_yz_xyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yy_xyz[i] * gfe2_0 + 2.0 * ts_yy_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_xy[i] * gfe2_0 + 4.0 * ts_yyz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xyzz[i] * gfe_0 + ts_yyz_xyzz[i] * rgc2_0;

        gr_yyz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe2_0 + 4.0 * ts_yz_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_xzz[i] * gfe2_0 + 2.0 * ts_yy_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyz_xz[i] * gfe2_0 + 6.0 * ts_yyz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xzzz[i] * gfe_0 + ts_yyz_xzzz[i] * rgc2_0;

        gr_yyz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe2_0 + 16.0 * ts_yz_yyy[i] * gfe2_0 + 4.0 * ts_yz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyz_yy[i] * gfe2_0 + 8.0 * ts_yyz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_yyyy[i] * gfe_0 + ts_yyz_yyyy[i] * rgc2_0;

        gr_yyz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe2_0 + 12.0 * ts_yz_yyz[i] * gfe2_0 + 4.0 * ts_yz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyy[i] * gfe2_0 + 2.0 * ts_yy_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyz_yz[i] * gfe2_0 + 6.0 * ts_yyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_yyyz[i] * gfe_0 + ts_yyz_yyyz[i] * rgc2_0;

        gr_yyz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe2_0 + 8.0 * ts_yz_yzz[i] * gfe2_0 + 4.0 * ts_yz_yyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yy_yyz[i] * gfe2_0 + 2.0 * ts_yy_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_zz[i] * gfe2_0 + 4.0 * ts_yyz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_yy[i] * gfe2_0 + 4.0 * ts_yyz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_yyzz[i] * gfe_0 + ts_yyz_yyzz[i] * rgc2_0;

        gr_yyz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe2_0 + 4.0 * ts_yz_zzz[i] * gfe2_0 + 4.0 * ts_yz_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_yzz[i] * gfe2_0 + 2.0 * ts_yy_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyz_yz[i] * gfe2_0 + 6.0 * ts_yyz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_yzzz[i] * gfe_0 + ts_yyz_yzzz[i] * rgc2_0;

        gr_yyz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe2_0 + 4.0 * ts_yz_zzzz[i] * gfe_0 * gc_y[i] + 8.0 * ts_yy_zzz[i] * gfe2_0 + 2.0 * ts_yy_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyz_zz[i] * gfe2_0 + 8.0 * ts_yyz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_zzzz[i] * gfe_0 + ts_yyz_zzzz[i] * rgc2_0;
    }

    // Set up 120-135 components of targeted buffer : FG

    auto gr_yzz_xxxx = pbuffer.data(idx_g_fg + 120);

    auto gr_yzz_xxxy = pbuffer.data(idx_g_fg + 121);

    auto gr_yzz_xxxz = pbuffer.data(idx_g_fg + 122);

    auto gr_yzz_xxyy = pbuffer.data(idx_g_fg + 123);

    auto gr_yzz_xxyz = pbuffer.data(idx_g_fg + 124);

    auto gr_yzz_xxzz = pbuffer.data(idx_g_fg + 125);

    auto gr_yzz_xyyy = pbuffer.data(idx_g_fg + 126);

    auto gr_yzz_xyyz = pbuffer.data(idx_g_fg + 127);

    auto gr_yzz_xyzz = pbuffer.data(idx_g_fg + 128);

    auto gr_yzz_xzzz = pbuffer.data(idx_g_fg + 129);

    auto gr_yzz_yyyy = pbuffer.data(idx_g_fg + 130);

    auto gr_yzz_yyyz = pbuffer.data(idx_g_fg + 131);

    auto gr_yzz_yyzz = pbuffer.data(idx_g_fg + 132);

    auto gr_yzz_yzzz = pbuffer.data(idx_g_fg + 133);

    auto gr_yzz_zzzz = pbuffer.data(idx_g_fg + 134);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xxxx, gr_yzz_xxxy, gr_yzz_xxxz, gr_yzz_xxyy, gr_yzz_xxyz, gr_yzz_xxzz, gr_yzz_xyyy, gr_yzz_xyyz, gr_yzz_xyzz, gr_yzz_xzzz, gr_yzz_yyyy, gr_yzz_yyyz, gr_yzz_yyzz, gr_yzz_yzzz, gr_yzz_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zz, ts_yzz_zzz, ts_yzz_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yzz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxxx[i] * gfe2_0 + 4.0 * ts_yz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_yzz_xx[i] * gfe2_0 + 8.0 * ts_yzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzz_xxxx[i] * gfe_0 + ts_yzz_xxxx[i] * rgc2_0;

        gr_yzz_xxxy[i] = 2.0 * ts_zz_xxx[i] * gfe2_0 + 2.0 * ts_zz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxxy[i] * gfe2_0 + 4.0 * ts_yz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzz_xy[i] * gfe2_0 + 6.0 * ts_yzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_xxxy[i] * gfe_0 + ts_yzz_xxxy[i] * rgc2_0;

        gr_yzz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxxz[i] * gfe2_0 + 4.0 * ts_yz_xxx[i] * gfe2_0 + 4.0 * ts_yz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzz_xz[i] * gfe2_0 + 6.0 * ts_yzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xxxz[i] * gfe_0 + ts_yzz_xxxz[i] * rgc2_0;

        gr_yzz_xxyy[i] = 4.0 * ts_zz_xxy[i] * gfe2_0 + 2.0 * ts_zz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxyy[i] * gfe2_0 + 4.0 * ts_yz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yy[i] * gfe2_0 + 4.0 * ts_yzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xx[i] * gfe2_0 + 4.0 * ts_yzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_xxyy[i] * gfe_0 + ts_yzz_xxyy[i] * rgc2_0;

        gr_yzz_xxyz[i] = 2.0 * ts_zz_xxz[i] * gfe2_0 + 2.0 * ts_zz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxyz[i] * gfe2_0 + 4.0 * ts_yz_xxy[i] * gfe2_0 + 4.0 * ts_yz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yz[i] * gfe2_0 + 4.0 * ts_yzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xxyz[i] * gfe_0 + ts_yzz_xxyz[i] * rgc2_0;

        gr_yzz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxzz[i] * gfe2_0 + 8.0 * ts_yz_xxz[i] * gfe2_0 + 4.0 * ts_yz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_zz[i] * gfe2_0 + 4.0 * ts_yzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xx[i] * gfe2_0 + 4.0 * ts_yzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xxzz[i] * gfe_0 + ts_yzz_xxzz[i] * rgc2_0;

        gr_yzz_xyyy[i] = 6.0 * ts_zz_xyy[i] * gfe2_0 + 2.0 * ts_zz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyyy[i] * gfe2_0 + 4.0 * ts_yz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yzz_xy[i] * gfe2_0 + 6.0 * ts_yzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_xyyy[i] * gfe_0 + ts_yzz_xyyy[i] * rgc2_0;

        gr_yzz_xyyz[i] = 4.0 * ts_zz_xyz[i] * gfe2_0 + 2.0 * ts_zz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyyz[i] * gfe2_0 + 4.0 * ts_yz_xyy[i] * gfe2_0 + 4.0 * ts_yz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xz[i] * gfe2_0 + 4.0 * ts_yzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xyyz[i] * gfe_0 + ts_yzz_xyyz[i] * rgc2_0;

        gr_yzz_xyzz[i] = 2.0 * ts_zz_xzz[i] * gfe2_0 + 2.0 * ts_zz_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyzz[i] * gfe2_0 + 8.0 * ts_yz_xyz[i] * gfe2_0 + 4.0 * ts_yz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_xy[i] * gfe2_0 + 4.0 * ts_yzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xyzz[i] * gfe_0 + ts_yzz_xyzz[i] * rgc2_0;

        gr_yzz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xzzz[i] * gfe2_0 + 12.0 * ts_yz_xzz[i] * gfe2_0 + 4.0 * ts_yz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yzz_xz[i] * gfe2_0 + 6.0 * ts_yzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xzzz[i] * gfe_0 + ts_yzz_xzzz[i] * rgc2_0;

        gr_yzz_yyyy[i] = 8.0 * ts_zz_yyy[i] * gfe2_0 + 2.0 * ts_zz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyyy[i] * gfe2_0 + 4.0 * ts_yz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_yzz_yy[i] * gfe2_0 + 8.0 * ts_yzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_yyyy[i] * gfe_0 + ts_yzz_yyyy[i] * rgc2_0;

        gr_yzz_yyyz[i] = 6.0 * ts_zz_yyz[i] * gfe2_0 + 2.0 * ts_zz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyyz[i] * gfe2_0 + 4.0 * ts_yz_yyy[i] * gfe2_0 + 4.0 * ts_yz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzz_yz[i] * gfe2_0 + 6.0 * ts_yzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_yyyz[i] * gfe_0 + ts_yzz_yyyz[i] * rgc2_0;

        gr_yzz_yyzz[i] = 4.0 * ts_zz_yzz[i] * gfe2_0 + 2.0 * ts_zz_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyzz[i] * gfe2_0 + 8.0 * ts_yz_yyz[i] * gfe2_0 + 4.0 * ts_yz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_zz[i] * gfe2_0 + 4.0 * ts_yzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_yy[i] * gfe2_0 + 4.0 * ts_yzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_yyzz[i] * gfe_0 + ts_yzz_yyzz[i] * rgc2_0;

        gr_yzz_yzzz[i] = 2.0 * ts_zz_zzz[i] * gfe2_0 + 2.0 * ts_zz_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yzzz[i] * gfe2_0 + 12.0 * ts_yz_yzz[i] * gfe2_0 + 4.0 * ts_yz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yzz_yz[i] * gfe2_0 + 6.0 * ts_yzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_yzzz[i] * gfe_0 + ts_yzz_yzzz[i] * rgc2_0;

        gr_yzz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zzzz[i] * gfe2_0 + 16.0 * ts_yz_zzz[i] * gfe2_0 + 4.0 * ts_yz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_yzz_zz[i] * gfe2_0 + 8.0 * ts_yzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_zzzz[i] * gfe_0 + ts_yzz_zzzz[i] * rgc2_0;
    }

    // Set up 135-150 components of targeted buffer : FG

    auto gr_zzz_xxxx = pbuffer.data(idx_g_fg + 135);

    auto gr_zzz_xxxy = pbuffer.data(idx_g_fg + 136);

    auto gr_zzz_xxxz = pbuffer.data(idx_g_fg + 137);

    auto gr_zzz_xxyy = pbuffer.data(idx_g_fg + 138);

    auto gr_zzz_xxyz = pbuffer.data(idx_g_fg + 139);

    auto gr_zzz_xxzz = pbuffer.data(idx_g_fg + 140);

    auto gr_zzz_xyyy = pbuffer.data(idx_g_fg + 141);

    auto gr_zzz_xyyz = pbuffer.data(idx_g_fg + 142);

    auto gr_zzz_xyzz = pbuffer.data(idx_g_fg + 143);

    auto gr_zzz_xzzz = pbuffer.data(idx_g_fg + 144);

    auto gr_zzz_yyyy = pbuffer.data(idx_g_fg + 145);

    auto gr_zzz_yyyz = pbuffer.data(idx_g_fg + 146);

    auto gr_zzz_yyzz = pbuffer.data(idx_g_fg + 147);

    auto gr_zzz_yzzz = pbuffer.data(idx_g_fg + 148);

    auto gr_zzz_zzzz = pbuffer.data(idx_g_fg + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xxxx, gr_zzz_xxxy, gr_zzz_xxxz, gr_zzz_xxyy, gr_zzz_xxyz, gr_zzz_xxzz, gr_zzz_xyyy, gr_zzz_xyyz, gr_zzz_xyzz, gr_zzz_xzzz, gr_zzz_yyyy, gr_zzz_yyyz, gr_zzz_yyzz, gr_zzz_yzzz, gr_zzz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_zzz_xxxx[i] = 6.0 * ts_z_xxxx[i] * gfe2_0 + 6.0 * ts_zz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_zzz_xx[i] * gfe2_0 + 8.0 * ts_zzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzz_xxxx[i] * gfe_0 + ts_zzz_xxxx[i] * rgc2_0;

        gr_zzz_xxxy[i] = 6.0 * ts_z_xxxy[i] * gfe2_0 + 6.0 * ts_zz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzz_xy[i] * gfe2_0 + 6.0 * ts_zzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_xxxy[i] * gfe_0 + ts_zzz_xxxy[i] * rgc2_0;

        gr_zzz_xxxz[i] = 6.0 * ts_z_xxxz[i] * gfe2_0 + 6.0 * ts_zz_xxx[i] * gfe2_0 + 6.0 * ts_zz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzz_xz[i] * gfe2_0 + 6.0 * ts_zzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xxxz[i] * gfe_0 + ts_zzz_xxxz[i] * rgc2_0;

        gr_zzz_xxyy[i] = 6.0 * ts_z_xxyy[i] * gfe2_0 + 6.0 * ts_zz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yy[i] * gfe2_0 + 4.0 * ts_zzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xx[i] * gfe2_0 + 4.0 * ts_zzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_xxyy[i] * gfe_0 + ts_zzz_xxyy[i] * rgc2_0;

        gr_zzz_xxyz[i] = 6.0 * ts_z_xxyz[i] * gfe2_0 + 6.0 * ts_zz_xxy[i] * gfe2_0 + 6.0 * ts_zz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yz[i] * gfe2_0 + 4.0 * ts_zzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xxyz[i] * gfe_0 + ts_zzz_xxyz[i] * rgc2_0;

        gr_zzz_xxzz[i] = 6.0 * ts_z_xxzz[i] * gfe2_0 + 12.0 * ts_zz_xxz[i] * gfe2_0 + 6.0 * ts_zz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_zz[i] * gfe2_0 + 4.0 * ts_zzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xx[i] * gfe2_0 + 4.0 * ts_zzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xxzz[i] * gfe_0 + ts_zzz_xxzz[i] * rgc2_0;

        gr_zzz_xyyy[i] = 6.0 * ts_z_xyyy[i] * gfe2_0 + 6.0 * ts_zz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_zzz_xy[i] * gfe2_0 + 6.0 * ts_zzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_xyyy[i] * gfe_0 + ts_zzz_xyyy[i] * rgc2_0;

        gr_zzz_xyyz[i] = 6.0 * ts_z_xyyz[i] * gfe2_0 + 6.0 * ts_zz_xyy[i] * gfe2_0 + 6.0 * ts_zz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xz[i] * gfe2_0 + 4.0 * ts_zzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xyyz[i] * gfe_0 + ts_zzz_xyyz[i] * rgc2_0;

        gr_zzz_xyzz[i] = 6.0 * ts_z_xyzz[i] * gfe2_0 + 12.0 * ts_zz_xyz[i] * gfe2_0 + 6.0 * ts_zz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_xy[i] * gfe2_0 + 4.0 * ts_zzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xyzz[i] * gfe_0 + ts_zzz_xyzz[i] * rgc2_0;

        gr_zzz_xzzz[i] = 6.0 * ts_z_xzzz[i] * gfe2_0 + 18.0 * ts_zz_xzz[i] * gfe2_0 + 6.0 * ts_zz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_zzz_xz[i] * gfe2_0 + 6.0 * ts_zzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xzzz[i] * gfe_0 + ts_zzz_xzzz[i] * rgc2_0;

        gr_zzz_yyyy[i] = 6.0 * ts_z_yyyy[i] * gfe2_0 + 6.0 * ts_zz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_zzz_yy[i] * gfe2_0 + 8.0 * ts_zzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_yyyy[i] * gfe_0 + ts_zzz_yyyy[i] * rgc2_0;

        gr_zzz_yyyz[i] = 6.0 * ts_z_yyyz[i] * gfe2_0 + 6.0 * ts_zz_yyy[i] * gfe2_0 + 6.0 * ts_zz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzz_yz[i] * gfe2_0 + 6.0 * ts_zzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_yyyz[i] * gfe_0 + ts_zzz_yyyz[i] * rgc2_0;

        gr_zzz_yyzz[i] = 6.0 * ts_z_yyzz[i] * gfe2_0 + 12.0 * ts_zz_yyz[i] * gfe2_0 + 6.0 * ts_zz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_zz[i] * gfe2_0 + 4.0 * ts_zzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_yy[i] * gfe2_0 + 4.0 * ts_zzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_yyzz[i] * gfe_0 + ts_zzz_yyzz[i] * rgc2_0;

        gr_zzz_yzzz[i] = 6.0 * ts_z_yzzz[i] * gfe2_0 + 18.0 * ts_zz_yzz[i] * gfe2_0 + 6.0 * ts_zz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_zzz_yz[i] * gfe2_0 + 6.0 * ts_zzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_yzzz[i] * gfe_0 + ts_zzz_yzzz[i] * rgc2_0;

        gr_zzz_zzzz[i] = 6.0 * ts_z_zzzz[i] * gfe2_0 + 24.0 * ts_zz_zzz[i] * gfe2_0 + 6.0 * ts_zz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_zzz_zz[i] * gfe2_0 + 8.0 * ts_zzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_zzzz[i] * gfe_0 + ts_zzz_zzzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

