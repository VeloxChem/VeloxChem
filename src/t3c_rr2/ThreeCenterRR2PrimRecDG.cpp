#include "ThreeCenterRR2PrimRecDG.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_dg(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_dg,
                  const size_t idx_pg,
                  const size_t idx_g_pg,
                  const size_t idx_df,
                  const size_t idx_g_df,
                  const size_t idx_dg,
                  const size_t idx_g_dg,
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

    // Set up components of auxiliary buffer : PG

    auto gr_x_xxxx = pbuffer.data(idx_g_pg);

    auto gr_x_xxxy = pbuffer.data(idx_g_pg + 1);

    auto gr_x_xxxz = pbuffer.data(idx_g_pg + 2);

    auto gr_x_xxyy = pbuffer.data(idx_g_pg + 3);

    auto gr_x_xxyz = pbuffer.data(idx_g_pg + 4);

    auto gr_x_xxzz = pbuffer.data(idx_g_pg + 5);

    auto gr_x_xyyy = pbuffer.data(idx_g_pg + 6);

    auto gr_x_xyyz = pbuffer.data(idx_g_pg + 7);

    auto gr_x_xyzz = pbuffer.data(idx_g_pg + 8);

    auto gr_x_xzzz = pbuffer.data(idx_g_pg + 9);

    auto gr_x_yyyy = pbuffer.data(idx_g_pg + 10);

    auto gr_x_yyyz = pbuffer.data(idx_g_pg + 11);

    auto gr_x_yyzz = pbuffer.data(idx_g_pg + 12);

    auto gr_x_yzzz = pbuffer.data(idx_g_pg + 13);

    auto gr_x_zzzz = pbuffer.data(idx_g_pg + 14);

    auto gr_y_xxxx = pbuffer.data(idx_g_pg + 15);

    auto gr_y_xxxy = pbuffer.data(idx_g_pg + 16);

    auto gr_y_xxxz = pbuffer.data(idx_g_pg + 17);

    auto gr_y_xxyy = pbuffer.data(idx_g_pg + 18);

    auto gr_y_xxyz = pbuffer.data(idx_g_pg + 19);

    auto gr_y_xxzz = pbuffer.data(idx_g_pg + 20);

    auto gr_y_xyyy = pbuffer.data(idx_g_pg + 21);

    auto gr_y_xyyz = pbuffer.data(idx_g_pg + 22);

    auto gr_y_xyzz = pbuffer.data(idx_g_pg + 23);

    auto gr_y_xzzz = pbuffer.data(idx_g_pg + 24);

    auto gr_y_yyyy = pbuffer.data(idx_g_pg + 25);

    auto gr_y_yyyz = pbuffer.data(idx_g_pg + 26);

    auto gr_y_yyzz = pbuffer.data(idx_g_pg + 27);

    auto gr_y_yzzz = pbuffer.data(idx_g_pg + 28);

    auto gr_y_zzzz = pbuffer.data(idx_g_pg + 29);

    auto gr_z_xxxx = pbuffer.data(idx_g_pg + 30);

    auto gr_z_xxxy = pbuffer.data(idx_g_pg + 31);

    auto gr_z_xxxz = pbuffer.data(idx_g_pg + 32);

    auto gr_z_xxyy = pbuffer.data(idx_g_pg + 33);

    auto gr_z_xxyz = pbuffer.data(idx_g_pg + 34);

    auto gr_z_xxzz = pbuffer.data(idx_g_pg + 35);

    auto gr_z_xyyy = pbuffer.data(idx_g_pg + 36);

    auto gr_z_xyyz = pbuffer.data(idx_g_pg + 37);

    auto gr_z_xyzz = pbuffer.data(idx_g_pg + 38);

    auto gr_z_xzzz = pbuffer.data(idx_g_pg + 39);

    auto gr_z_yyyy = pbuffer.data(idx_g_pg + 40);

    auto gr_z_yyyz = pbuffer.data(idx_g_pg + 41);

    auto gr_z_yyzz = pbuffer.data(idx_g_pg + 42);

    auto gr_z_yzzz = pbuffer.data(idx_g_pg + 43);

    auto gr_z_zzzz = pbuffer.data(idx_g_pg + 44);

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

    // Set up components of auxiliary buffer : DG

    auto gr_xx_xxxx = pbuffer.data(idx_g_dg);

    auto gr_xx_xxxy = pbuffer.data(idx_g_dg + 1);

    auto gr_xx_xxxz = pbuffer.data(idx_g_dg + 2);

    auto gr_xx_xxyy = pbuffer.data(idx_g_dg + 3);

    auto gr_xx_xxyz = pbuffer.data(idx_g_dg + 4);

    auto gr_xx_xxzz = pbuffer.data(idx_g_dg + 5);

    auto gr_xx_xyyy = pbuffer.data(idx_g_dg + 6);

    auto gr_xx_xyyz = pbuffer.data(idx_g_dg + 7);

    auto gr_xx_xyzz = pbuffer.data(idx_g_dg + 8);

    auto gr_xx_xzzz = pbuffer.data(idx_g_dg + 9);

    auto gr_xx_yyyy = pbuffer.data(idx_g_dg + 10);

    auto gr_xx_yyyz = pbuffer.data(idx_g_dg + 11);

    auto gr_xx_yyzz = pbuffer.data(idx_g_dg + 12);

    auto gr_xx_yzzz = pbuffer.data(idx_g_dg + 13);

    auto gr_xx_zzzz = pbuffer.data(idx_g_dg + 14);

    auto gr_xy_xxxx = pbuffer.data(idx_g_dg + 15);

    auto gr_xy_xxxy = pbuffer.data(idx_g_dg + 16);

    auto gr_xy_xxxz = pbuffer.data(idx_g_dg + 17);

    auto gr_xy_xxyy = pbuffer.data(idx_g_dg + 18);

    auto gr_xy_xxyz = pbuffer.data(idx_g_dg + 19);

    auto gr_xy_xxzz = pbuffer.data(idx_g_dg + 20);

    auto gr_xy_xyyy = pbuffer.data(idx_g_dg + 21);

    auto gr_xy_xyyz = pbuffer.data(idx_g_dg + 22);

    auto gr_xy_xyzz = pbuffer.data(idx_g_dg + 23);

    auto gr_xy_xzzz = pbuffer.data(idx_g_dg + 24);

    auto gr_xy_yyyy = pbuffer.data(idx_g_dg + 25);

    auto gr_xy_yyyz = pbuffer.data(idx_g_dg + 26);

    auto gr_xy_yyzz = pbuffer.data(idx_g_dg + 27);

    auto gr_xy_yzzz = pbuffer.data(idx_g_dg + 28);

    auto gr_xy_zzzz = pbuffer.data(idx_g_dg + 29);

    auto gr_xz_xxxx = pbuffer.data(idx_g_dg + 30);

    auto gr_xz_xxxy = pbuffer.data(idx_g_dg + 31);

    auto gr_xz_xxxz = pbuffer.data(idx_g_dg + 32);

    auto gr_xz_xxyy = pbuffer.data(idx_g_dg + 33);

    auto gr_xz_xxyz = pbuffer.data(idx_g_dg + 34);

    auto gr_xz_xxzz = pbuffer.data(idx_g_dg + 35);

    auto gr_xz_xyyy = pbuffer.data(idx_g_dg + 36);

    auto gr_xz_xyyz = pbuffer.data(idx_g_dg + 37);

    auto gr_xz_xyzz = pbuffer.data(idx_g_dg + 38);

    auto gr_xz_xzzz = pbuffer.data(idx_g_dg + 39);

    auto gr_xz_yyyy = pbuffer.data(idx_g_dg + 40);

    auto gr_xz_yyyz = pbuffer.data(idx_g_dg + 41);

    auto gr_xz_yyzz = pbuffer.data(idx_g_dg + 42);

    auto gr_xz_yzzz = pbuffer.data(idx_g_dg + 43);

    auto gr_xz_zzzz = pbuffer.data(idx_g_dg + 44);

    auto gr_yy_xxxx = pbuffer.data(idx_g_dg + 45);

    auto gr_yy_xxxy = pbuffer.data(idx_g_dg + 46);

    auto gr_yy_xxxz = pbuffer.data(idx_g_dg + 47);

    auto gr_yy_xxyy = pbuffer.data(idx_g_dg + 48);

    auto gr_yy_xxyz = pbuffer.data(idx_g_dg + 49);

    auto gr_yy_xxzz = pbuffer.data(idx_g_dg + 50);

    auto gr_yy_xyyy = pbuffer.data(idx_g_dg + 51);

    auto gr_yy_xyyz = pbuffer.data(idx_g_dg + 52);

    auto gr_yy_xyzz = pbuffer.data(idx_g_dg + 53);

    auto gr_yy_xzzz = pbuffer.data(idx_g_dg + 54);

    auto gr_yy_yyyy = pbuffer.data(idx_g_dg + 55);

    auto gr_yy_yyyz = pbuffer.data(idx_g_dg + 56);

    auto gr_yy_yyzz = pbuffer.data(idx_g_dg + 57);

    auto gr_yy_yzzz = pbuffer.data(idx_g_dg + 58);

    auto gr_yy_zzzz = pbuffer.data(idx_g_dg + 59);

    auto gr_yz_xxxx = pbuffer.data(idx_g_dg + 60);

    auto gr_yz_xxxy = pbuffer.data(idx_g_dg + 61);

    auto gr_yz_xxxz = pbuffer.data(idx_g_dg + 62);

    auto gr_yz_xxyy = pbuffer.data(idx_g_dg + 63);

    auto gr_yz_xxyz = pbuffer.data(idx_g_dg + 64);

    auto gr_yz_xxzz = pbuffer.data(idx_g_dg + 65);

    auto gr_yz_xyyy = pbuffer.data(idx_g_dg + 66);

    auto gr_yz_xyyz = pbuffer.data(idx_g_dg + 67);

    auto gr_yz_xyzz = pbuffer.data(idx_g_dg + 68);

    auto gr_yz_xzzz = pbuffer.data(idx_g_dg + 69);

    auto gr_yz_yyyy = pbuffer.data(idx_g_dg + 70);

    auto gr_yz_yyyz = pbuffer.data(idx_g_dg + 71);

    auto gr_yz_yyzz = pbuffer.data(idx_g_dg + 72);

    auto gr_yz_yzzz = pbuffer.data(idx_g_dg + 73);

    auto gr_yz_zzzz = pbuffer.data(idx_g_dg + 74);

    auto gr_zz_xxxx = pbuffer.data(idx_g_dg + 75);

    auto gr_zz_xxxy = pbuffer.data(idx_g_dg + 76);

    auto gr_zz_xxxz = pbuffer.data(idx_g_dg + 77);

    auto gr_zz_xxyy = pbuffer.data(idx_g_dg + 78);

    auto gr_zz_xxyz = pbuffer.data(idx_g_dg + 79);

    auto gr_zz_xxzz = pbuffer.data(idx_g_dg + 80);

    auto gr_zz_xyyy = pbuffer.data(idx_g_dg + 81);

    auto gr_zz_xyyz = pbuffer.data(idx_g_dg + 82);

    auto gr_zz_xyzz = pbuffer.data(idx_g_dg + 83);

    auto gr_zz_xzzz = pbuffer.data(idx_g_dg + 84);

    auto gr_zz_yyyy = pbuffer.data(idx_g_dg + 85);

    auto gr_zz_yyyz = pbuffer.data(idx_g_dg + 86);

    auto gr_zz_yyzz = pbuffer.data(idx_g_dg + 87);

    auto gr_zz_yzzz = pbuffer.data(idx_g_dg + 88);

    auto gr_zz_zzzz = pbuffer.data(idx_g_dg + 89);

    // Set up 0-15 components of targeted buffer : DG

    auto grr_x_xx_xxxx = pbuffer.data(idx_gr_dg);

    auto grr_x_xx_xxxy = pbuffer.data(idx_gr_dg + 1);

    auto grr_x_xx_xxxz = pbuffer.data(idx_gr_dg + 2);

    auto grr_x_xx_xxyy = pbuffer.data(idx_gr_dg + 3);

    auto grr_x_xx_xxyz = pbuffer.data(idx_gr_dg + 4);

    auto grr_x_xx_xxzz = pbuffer.data(idx_gr_dg + 5);

    auto grr_x_xx_xyyy = pbuffer.data(idx_gr_dg + 6);

    auto grr_x_xx_xyyz = pbuffer.data(idx_gr_dg + 7);

    auto grr_x_xx_xyzz = pbuffer.data(idx_gr_dg + 8);

    auto grr_x_xx_xzzz = pbuffer.data(idx_gr_dg + 9);

    auto grr_x_xx_yyyy = pbuffer.data(idx_gr_dg + 10);

    auto grr_x_xx_yyyz = pbuffer.data(idx_gr_dg + 11);

    auto grr_x_xx_yyzz = pbuffer.data(idx_gr_dg + 12);

    auto grr_x_xx_yzzz = pbuffer.data(idx_gr_dg + 13);

    auto grr_x_xx_zzzz = pbuffer.data(idx_gr_dg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxyy, gr_x_xxyz, gr_x_xxzz, gr_x_xyyy, gr_x_xyyz, gr_x_xyzz, gr_x_xzzz, gr_x_yyyy, gr_x_yyyz, gr_x_yyzz, gr_x_yzzz, gr_x_zzzz, gr_xx_xxx, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxy, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxz, gr_xx_xxzz, gr_xx_xyy, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyz, gr_xx_xyzz, gr_xx_xzz, gr_xx_xzzz, gr_xx_yyy, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyz, gr_xx_yyzz, gr_xx_yzz, gr_xx_yzzz, gr_xx_zzz, gr_xx_zzzz, grr_x_xx_xxxx, grr_x_xx_xxxy, grr_x_xx_xxxz, grr_x_xx_xxyy, grr_x_xx_xxyz, grr_x_xx_xxzz, grr_x_xx_xyyy, grr_x_xx_xyyz, grr_x_xx_xyzz, grr_x_xx_xzzz, grr_x_xx_yyyy, grr_x_xx_yyyz, grr_x_xx_yyzz, grr_x_xx_yzzz, grr_x_xx_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xx_xxxx[i] = 4.0 * ts_x_xxxx[i] * gfe2_0 + 2.0 * gr_x_xxxx[i] * gfe_0 + 8.0 * ts_xx_xxx[i] * gfe2_0 + 4.0 * gr_xx_xxx[i] * gfe_0 + 2.0 * ts_xx_xxxx[i] * gfe_0 * gc_x[i] + gr_xx_xxxx[i] * gc_x[i];

        grr_x_xx_xxxy[i] = 4.0 * ts_x_xxxy[i] * gfe2_0 + 2.0 * gr_x_xxxy[i] * gfe_0 + 6.0 * ts_xx_xxy[i] * gfe2_0 + 3.0 * gr_xx_xxy[i] * gfe_0 + 2.0 * ts_xx_xxxy[i] * gfe_0 * gc_x[i] + gr_xx_xxxy[i] * gc_x[i];

        grr_x_xx_xxxz[i] = 4.0 * ts_x_xxxz[i] * gfe2_0 + 2.0 * gr_x_xxxz[i] * gfe_0 + 6.0 * ts_xx_xxz[i] * gfe2_0 + 3.0 * gr_xx_xxz[i] * gfe_0 + 2.0 * ts_xx_xxxz[i] * gfe_0 * gc_x[i] + gr_xx_xxxz[i] * gc_x[i];

        grr_x_xx_xxyy[i] = 4.0 * ts_x_xxyy[i] * gfe2_0 + 2.0 * gr_x_xxyy[i] * gfe_0 + 4.0 * ts_xx_xyy[i] * gfe2_0 + 2.0 * gr_xx_xyy[i] * gfe_0 + 2.0 * ts_xx_xxyy[i] * gfe_0 * gc_x[i] + gr_xx_xxyy[i] * gc_x[i];

        grr_x_xx_xxyz[i] = 4.0 * ts_x_xxyz[i] * gfe2_0 + 2.0 * gr_x_xxyz[i] * gfe_0 + 4.0 * ts_xx_xyz[i] * gfe2_0 + 2.0 * gr_xx_xyz[i] * gfe_0 + 2.0 * ts_xx_xxyz[i] * gfe_0 * gc_x[i] + gr_xx_xxyz[i] * gc_x[i];

        grr_x_xx_xxzz[i] = 4.0 * ts_x_xxzz[i] * gfe2_0 + 2.0 * gr_x_xxzz[i] * gfe_0 + 4.0 * ts_xx_xzz[i] * gfe2_0 + 2.0 * gr_xx_xzz[i] * gfe_0 + 2.0 * ts_xx_xxzz[i] * gfe_0 * gc_x[i] + gr_xx_xxzz[i] * gc_x[i];

        grr_x_xx_xyyy[i] = 4.0 * ts_x_xyyy[i] * gfe2_0 + 2.0 * gr_x_xyyy[i] * gfe_0 + 2.0 * ts_xx_yyy[i] * gfe2_0 + gr_xx_yyy[i] * gfe_0 + 2.0 * ts_xx_xyyy[i] * gfe_0 * gc_x[i] + gr_xx_xyyy[i] * gc_x[i];

        grr_x_xx_xyyz[i] = 4.0 * ts_x_xyyz[i] * gfe2_0 + 2.0 * gr_x_xyyz[i] * gfe_0 + 2.0 * ts_xx_yyz[i] * gfe2_0 + gr_xx_yyz[i] * gfe_0 + 2.0 * ts_xx_xyyz[i] * gfe_0 * gc_x[i] + gr_xx_xyyz[i] * gc_x[i];

        grr_x_xx_xyzz[i] = 4.0 * ts_x_xyzz[i] * gfe2_0 + 2.0 * gr_x_xyzz[i] * gfe_0 + 2.0 * ts_xx_yzz[i] * gfe2_0 + gr_xx_yzz[i] * gfe_0 + 2.0 * ts_xx_xyzz[i] * gfe_0 * gc_x[i] + gr_xx_xyzz[i] * gc_x[i];

        grr_x_xx_xzzz[i] = 4.0 * ts_x_xzzz[i] * gfe2_0 + 2.0 * gr_x_xzzz[i] * gfe_0 + 2.0 * ts_xx_zzz[i] * gfe2_0 + gr_xx_zzz[i] * gfe_0 + 2.0 * ts_xx_xzzz[i] * gfe_0 * gc_x[i] + gr_xx_xzzz[i] * gc_x[i];

        grr_x_xx_yyyy[i] = 4.0 * ts_x_yyyy[i] * gfe2_0 + 2.0 * gr_x_yyyy[i] * gfe_0 + 2.0 * ts_xx_yyyy[i] * gfe_0 * gc_x[i] + gr_xx_yyyy[i] * gc_x[i];

        grr_x_xx_yyyz[i] = 4.0 * ts_x_yyyz[i] * gfe2_0 + 2.0 * gr_x_yyyz[i] * gfe_0 + 2.0 * ts_xx_yyyz[i] * gfe_0 * gc_x[i] + gr_xx_yyyz[i] * gc_x[i];

        grr_x_xx_yyzz[i] = 4.0 * ts_x_yyzz[i] * gfe2_0 + 2.0 * gr_x_yyzz[i] * gfe_0 + 2.0 * ts_xx_yyzz[i] * gfe_0 * gc_x[i] + gr_xx_yyzz[i] * gc_x[i];

        grr_x_xx_yzzz[i] = 4.0 * ts_x_yzzz[i] * gfe2_0 + 2.0 * gr_x_yzzz[i] * gfe_0 + 2.0 * ts_xx_yzzz[i] * gfe_0 * gc_x[i] + gr_xx_yzzz[i] * gc_x[i];

        grr_x_xx_zzzz[i] = 4.0 * ts_x_zzzz[i] * gfe2_0 + 2.0 * gr_x_zzzz[i] * gfe_0 + 2.0 * ts_xx_zzzz[i] * gfe_0 * gc_x[i] + gr_xx_zzzz[i] * gc_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto grr_x_xy_xxxx = pbuffer.data(idx_gr_dg + 15);

    auto grr_x_xy_xxxy = pbuffer.data(idx_gr_dg + 16);

    auto grr_x_xy_xxxz = pbuffer.data(idx_gr_dg + 17);

    auto grr_x_xy_xxyy = pbuffer.data(idx_gr_dg + 18);

    auto grr_x_xy_xxyz = pbuffer.data(idx_gr_dg + 19);

    auto grr_x_xy_xxzz = pbuffer.data(idx_gr_dg + 20);

    auto grr_x_xy_xyyy = pbuffer.data(idx_gr_dg + 21);

    auto grr_x_xy_xyyz = pbuffer.data(idx_gr_dg + 22);

    auto grr_x_xy_xyzz = pbuffer.data(idx_gr_dg + 23);

    auto grr_x_xy_xzzz = pbuffer.data(idx_gr_dg + 24);

    auto grr_x_xy_yyyy = pbuffer.data(idx_gr_dg + 25);

    auto grr_x_xy_yyyz = pbuffer.data(idx_gr_dg + 26);

    auto grr_x_xy_yyzz = pbuffer.data(idx_gr_dg + 27);

    auto grr_x_xy_yzzz = pbuffer.data(idx_gr_dg + 28);

    auto grr_x_xy_zzzz = pbuffer.data(idx_gr_dg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxx, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxy, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxz, gr_xy_xxzz, gr_xy_xyy, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyz, gr_xy_xyzz, gr_xy_xzz, gr_xy_xzzz, gr_xy_yyy, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyz, gr_xy_yyzz, gr_xy_yzz, gr_xy_yzzz, gr_xy_zzz, gr_xy_zzzz, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxyy, gr_y_xxyz, gr_y_xxzz, gr_y_xyyy, gr_y_xyyz, gr_y_xyzz, gr_y_xzzz, gr_y_yyyy, gr_y_yyyz, gr_y_yyzz, gr_y_yzzz, gr_y_zzzz, grr_x_xy_xxxx, grr_x_xy_xxxy, grr_x_xy_xxxz, grr_x_xy_xxyy, grr_x_xy_xxyz, grr_x_xy_xxzz, grr_x_xy_xyyy, grr_x_xy_xyyz, grr_x_xy_xyzz, grr_x_xy_xzzz, grr_x_xy_yyyy, grr_x_xy_yyyz, grr_x_xy_yyzz, grr_x_xy_yzzz, grr_x_xy_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xy_xxxx[i] = 2.0 * ts_y_xxxx[i] * gfe2_0 + gr_y_xxxx[i] * gfe_0 + 8.0 * ts_xy_xxx[i] * gfe2_0 + 4.0 * gr_xy_xxx[i] * gfe_0 + 2.0 * ts_xy_xxxx[i] * gfe_0 * gc_x[i] + gr_xy_xxxx[i] * gc_x[i];

        grr_x_xy_xxxy[i] = 2.0 * ts_y_xxxy[i] * gfe2_0 + gr_y_xxxy[i] * gfe_0 + 6.0 * ts_xy_xxy[i] * gfe2_0 + 3.0 * gr_xy_xxy[i] * gfe_0 + 2.0 * ts_xy_xxxy[i] * gfe_0 * gc_x[i] + gr_xy_xxxy[i] * gc_x[i];

        grr_x_xy_xxxz[i] = 2.0 * ts_y_xxxz[i] * gfe2_0 + gr_y_xxxz[i] * gfe_0 + 6.0 * ts_xy_xxz[i] * gfe2_0 + 3.0 * gr_xy_xxz[i] * gfe_0 + 2.0 * ts_xy_xxxz[i] * gfe_0 * gc_x[i] + gr_xy_xxxz[i] * gc_x[i];

        grr_x_xy_xxyy[i] = 2.0 * ts_y_xxyy[i] * gfe2_0 + gr_y_xxyy[i] * gfe_0 + 4.0 * ts_xy_xyy[i] * gfe2_0 + 2.0 * gr_xy_xyy[i] * gfe_0 + 2.0 * ts_xy_xxyy[i] * gfe_0 * gc_x[i] + gr_xy_xxyy[i] * gc_x[i];

        grr_x_xy_xxyz[i] = 2.0 * ts_y_xxyz[i] * gfe2_0 + gr_y_xxyz[i] * gfe_0 + 4.0 * ts_xy_xyz[i] * gfe2_0 + 2.0 * gr_xy_xyz[i] * gfe_0 + 2.0 * ts_xy_xxyz[i] * gfe_0 * gc_x[i] + gr_xy_xxyz[i] * gc_x[i];

        grr_x_xy_xxzz[i] = 2.0 * ts_y_xxzz[i] * gfe2_0 + gr_y_xxzz[i] * gfe_0 + 4.0 * ts_xy_xzz[i] * gfe2_0 + 2.0 * gr_xy_xzz[i] * gfe_0 + 2.0 * ts_xy_xxzz[i] * gfe_0 * gc_x[i] + gr_xy_xxzz[i] * gc_x[i];

        grr_x_xy_xyyy[i] = 2.0 * ts_y_xyyy[i] * gfe2_0 + gr_y_xyyy[i] * gfe_0 + 2.0 * ts_xy_yyy[i] * gfe2_0 + gr_xy_yyy[i] * gfe_0 + 2.0 * ts_xy_xyyy[i] * gfe_0 * gc_x[i] + gr_xy_xyyy[i] * gc_x[i];

        grr_x_xy_xyyz[i] = 2.0 * ts_y_xyyz[i] * gfe2_0 + gr_y_xyyz[i] * gfe_0 + 2.0 * ts_xy_yyz[i] * gfe2_0 + gr_xy_yyz[i] * gfe_0 + 2.0 * ts_xy_xyyz[i] * gfe_0 * gc_x[i] + gr_xy_xyyz[i] * gc_x[i];

        grr_x_xy_xyzz[i] = 2.0 * ts_y_xyzz[i] * gfe2_0 + gr_y_xyzz[i] * gfe_0 + 2.0 * ts_xy_yzz[i] * gfe2_0 + gr_xy_yzz[i] * gfe_0 + 2.0 * ts_xy_xyzz[i] * gfe_0 * gc_x[i] + gr_xy_xyzz[i] * gc_x[i];

        grr_x_xy_xzzz[i] = 2.0 * ts_y_xzzz[i] * gfe2_0 + gr_y_xzzz[i] * gfe_0 + 2.0 * ts_xy_zzz[i] * gfe2_0 + gr_xy_zzz[i] * gfe_0 + 2.0 * ts_xy_xzzz[i] * gfe_0 * gc_x[i] + gr_xy_xzzz[i] * gc_x[i];

        grr_x_xy_yyyy[i] = 2.0 * ts_y_yyyy[i] * gfe2_0 + gr_y_yyyy[i] * gfe_0 + 2.0 * ts_xy_yyyy[i] * gfe_0 * gc_x[i] + gr_xy_yyyy[i] * gc_x[i];

        grr_x_xy_yyyz[i] = 2.0 * ts_y_yyyz[i] * gfe2_0 + gr_y_yyyz[i] * gfe_0 + 2.0 * ts_xy_yyyz[i] * gfe_0 * gc_x[i] + gr_xy_yyyz[i] * gc_x[i];

        grr_x_xy_yyzz[i] = 2.0 * ts_y_yyzz[i] * gfe2_0 + gr_y_yyzz[i] * gfe_0 + 2.0 * ts_xy_yyzz[i] * gfe_0 * gc_x[i] + gr_xy_yyzz[i] * gc_x[i];

        grr_x_xy_yzzz[i] = 2.0 * ts_y_yzzz[i] * gfe2_0 + gr_y_yzzz[i] * gfe_0 + 2.0 * ts_xy_yzzz[i] * gfe_0 * gc_x[i] + gr_xy_yzzz[i] * gc_x[i];

        grr_x_xy_zzzz[i] = 2.0 * ts_y_zzzz[i] * gfe2_0 + gr_y_zzzz[i] * gfe_0 + 2.0 * ts_xy_zzzz[i] * gfe_0 * gc_x[i] + gr_xy_zzzz[i] * gc_x[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto grr_x_xz_xxxx = pbuffer.data(idx_gr_dg + 30);

    auto grr_x_xz_xxxy = pbuffer.data(idx_gr_dg + 31);

    auto grr_x_xz_xxxz = pbuffer.data(idx_gr_dg + 32);

    auto grr_x_xz_xxyy = pbuffer.data(idx_gr_dg + 33);

    auto grr_x_xz_xxyz = pbuffer.data(idx_gr_dg + 34);

    auto grr_x_xz_xxzz = pbuffer.data(idx_gr_dg + 35);

    auto grr_x_xz_xyyy = pbuffer.data(idx_gr_dg + 36);

    auto grr_x_xz_xyyz = pbuffer.data(idx_gr_dg + 37);

    auto grr_x_xz_xyzz = pbuffer.data(idx_gr_dg + 38);

    auto grr_x_xz_xzzz = pbuffer.data(idx_gr_dg + 39);

    auto grr_x_xz_yyyy = pbuffer.data(idx_gr_dg + 40);

    auto grr_x_xz_yyyz = pbuffer.data(idx_gr_dg + 41);

    auto grr_x_xz_yyzz = pbuffer.data(idx_gr_dg + 42);

    auto grr_x_xz_yzzz = pbuffer.data(idx_gr_dg + 43);

    auto grr_x_xz_zzzz = pbuffer.data(idx_gr_dg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xxx, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxy, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxz, gr_xz_xxzz, gr_xz_xyy, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyz, gr_xz_xyzz, gr_xz_xzz, gr_xz_xzzz, gr_xz_yyy, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyz, gr_xz_yyzz, gr_xz_yzz, gr_xz_yzzz, gr_xz_zzz, gr_xz_zzzz, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxyy, gr_z_xxyz, gr_z_xxzz, gr_z_xyyy, gr_z_xyyz, gr_z_xyzz, gr_z_xzzz, gr_z_yyyy, gr_z_yyyz, gr_z_yyzz, gr_z_yzzz, gr_z_zzzz, grr_x_xz_xxxx, grr_x_xz_xxxy, grr_x_xz_xxxz, grr_x_xz_xxyy, grr_x_xz_xxyz, grr_x_xz_xxzz, grr_x_xz_xyyy, grr_x_xz_xyyz, grr_x_xz_xyzz, grr_x_xz_xzzz, grr_x_xz_yyyy, grr_x_xz_yyyz, grr_x_xz_yyzz, grr_x_xz_yzzz, grr_x_xz_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe2_0 + gr_z_xxxx[i] * gfe_0 + 8.0 * ts_xz_xxx[i] * gfe2_0 + 4.0 * gr_xz_xxx[i] * gfe_0 + 2.0 * ts_xz_xxxx[i] * gfe_0 * gc_x[i] + gr_xz_xxxx[i] * gc_x[i];

        grr_x_xz_xxxy[i] = 2.0 * ts_z_xxxy[i] * gfe2_0 + gr_z_xxxy[i] * gfe_0 + 6.0 * ts_xz_xxy[i] * gfe2_0 + 3.0 * gr_xz_xxy[i] * gfe_0 + 2.0 * ts_xz_xxxy[i] * gfe_0 * gc_x[i] + gr_xz_xxxy[i] * gc_x[i];

        grr_x_xz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe2_0 + gr_z_xxxz[i] * gfe_0 + 6.0 * ts_xz_xxz[i] * gfe2_0 + 3.0 * gr_xz_xxz[i] * gfe_0 + 2.0 * ts_xz_xxxz[i] * gfe_0 * gc_x[i] + gr_xz_xxxz[i] * gc_x[i];

        grr_x_xz_xxyy[i] = 2.0 * ts_z_xxyy[i] * gfe2_0 + gr_z_xxyy[i] * gfe_0 + 4.0 * ts_xz_xyy[i] * gfe2_0 + 2.0 * gr_xz_xyy[i] * gfe_0 + 2.0 * ts_xz_xxyy[i] * gfe_0 * gc_x[i] + gr_xz_xxyy[i] * gc_x[i];

        grr_x_xz_xxyz[i] = 2.0 * ts_z_xxyz[i] * gfe2_0 + gr_z_xxyz[i] * gfe_0 + 4.0 * ts_xz_xyz[i] * gfe2_0 + 2.0 * gr_xz_xyz[i] * gfe_0 + 2.0 * ts_xz_xxyz[i] * gfe_0 * gc_x[i] + gr_xz_xxyz[i] * gc_x[i];

        grr_x_xz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe2_0 + gr_z_xxzz[i] * gfe_0 + 4.0 * ts_xz_xzz[i] * gfe2_0 + 2.0 * gr_xz_xzz[i] * gfe_0 + 2.0 * ts_xz_xxzz[i] * gfe_0 * gc_x[i] + gr_xz_xxzz[i] * gc_x[i];

        grr_x_xz_xyyy[i] = 2.0 * ts_z_xyyy[i] * gfe2_0 + gr_z_xyyy[i] * gfe_0 + 2.0 * ts_xz_yyy[i] * gfe2_0 + gr_xz_yyy[i] * gfe_0 + 2.0 * ts_xz_xyyy[i] * gfe_0 * gc_x[i] + gr_xz_xyyy[i] * gc_x[i];

        grr_x_xz_xyyz[i] = 2.0 * ts_z_xyyz[i] * gfe2_0 + gr_z_xyyz[i] * gfe_0 + 2.0 * ts_xz_yyz[i] * gfe2_0 + gr_xz_yyz[i] * gfe_0 + 2.0 * ts_xz_xyyz[i] * gfe_0 * gc_x[i] + gr_xz_xyyz[i] * gc_x[i];

        grr_x_xz_xyzz[i] = 2.0 * ts_z_xyzz[i] * gfe2_0 + gr_z_xyzz[i] * gfe_0 + 2.0 * ts_xz_yzz[i] * gfe2_0 + gr_xz_yzz[i] * gfe_0 + 2.0 * ts_xz_xyzz[i] * gfe_0 * gc_x[i] + gr_xz_xyzz[i] * gc_x[i];

        grr_x_xz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe2_0 + gr_z_xzzz[i] * gfe_0 + 2.0 * ts_xz_zzz[i] * gfe2_0 + gr_xz_zzz[i] * gfe_0 + 2.0 * ts_xz_xzzz[i] * gfe_0 * gc_x[i] + gr_xz_xzzz[i] * gc_x[i];

        grr_x_xz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe2_0 + gr_z_yyyy[i] * gfe_0 + 2.0 * ts_xz_yyyy[i] * gfe_0 * gc_x[i] + gr_xz_yyyy[i] * gc_x[i];

        grr_x_xz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe2_0 + gr_z_yyyz[i] * gfe_0 + 2.0 * ts_xz_yyyz[i] * gfe_0 * gc_x[i] + gr_xz_yyyz[i] * gc_x[i];

        grr_x_xz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe2_0 + gr_z_yyzz[i] * gfe_0 + 2.0 * ts_xz_yyzz[i] * gfe_0 * gc_x[i] + gr_xz_yyzz[i] * gc_x[i];

        grr_x_xz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe2_0 + gr_z_yzzz[i] * gfe_0 + 2.0 * ts_xz_yzzz[i] * gfe_0 * gc_x[i] + gr_xz_yzzz[i] * gc_x[i];

        grr_x_xz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe2_0 + gr_z_zzzz[i] * gfe_0 + 2.0 * ts_xz_zzzz[i] * gfe_0 * gc_x[i] + gr_xz_zzzz[i] * gc_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto grr_x_yy_xxxx = pbuffer.data(idx_gr_dg + 45);

    auto grr_x_yy_xxxy = pbuffer.data(idx_gr_dg + 46);

    auto grr_x_yy_xxxz = pbuffer.data(idx_gr_dg + 47);

    auto grr_x_yy_xxyy = pbuffer.data(idx_gr_dg + 48);

    auto grr_x_yy_xxyz = pbuffer.data(idx_gr_dg + 49);

    auto grr_x_yy_xxzz = pbuffer.data(idx_gr_dg + 50);

    auto grr_x_yy_xyyy = pbuffer.data(idx_gr_dg + 51);

    auto grr_x_yy_xyyz = pbuffer.data(idx_gr_dg + 52);

    auto grr_x_yy_xyzz = pbuffer.data(idx_gr_dg + 53);

    auto grr_x_yy_xzzz = pbuffer.data(idx_gr_dg + 54);

    auto grr_x_yy_yyyy = pbuffer.data(idx_gr_dg + 55);

    auto grr_x_yy_yyyz = pbuffer.data(idx_gr_dg + 56);

    auto grr_x_yy_yyzz = pbuffer.data(idx_gr_dg + 57);

    auto grr_x_yy_yzzz = pbuffer.data(idx_gr_dg + 58);

    auto grr_x_yy_zzzz = pbuffer.data(idx_gr_dg + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxx, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxy, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxz, gr_yy_xxzz, gr_yy_xyy, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyz, gr_yy_xyzz, gr_yy_xzz, gr_yy_xzzz, gr_yy_yyy, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyz, gr_yy_yyzz, gr_yy_yzz, gr_yy_yzzz, gr_yy_zzz, gr_yy_zzzz, grr_x_yy_xxxx, grr_x_yy_xxxy, grr_x_yy_xxxz, grr_x_yy_xxyy, grr_x_yy_xxyz, grr_x_yy_xxzz, grr_x_yy_xyyy, grr_x_yy_xyyz, grr_x_yy_xyzz, grr_x_yy_xzzz, grr_x_yy_yyyy, grr_x_yy_yyyz, grr_x_yy_yyzz, grr_x_yy_yzzz, grr_x_yy_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yy_xxxx[i] = 8.0 * ts_yy_xxx[i] * gfe2_0 + 4.0 * gr_yy_xxx[i] * gfe_0 + 2.0 * ts_yy_xxxx[i] * gfe_0 * gc_x[i] + gr_yy_xxxx[i] * gc_x[i];

        grr_x_yy_xxxy[i] = 6.0 * ts_yy_xxy[i] * gfe2_0 + 3.0 * gr_yy_xxy[i] * gfe_0 + 2.0 * ts_yy_xxxy[i] * gfe_0 * gc_x[i] + gr_yy_xxxy[i] * gc_x[i];

        grr_x_yy_xxxz[i] = 6.0 * ts_yy_xxz[i] * gfe2_0 + 3.0 * gr_yy_xxz[i] * gfe_0 + 2.0 * ts_yy_xxxz[i] * gfe_0 * gc_x[i] + gr_yy_xxxz[i] * gc_x[i];

        grr_x_yy_xxyy[i] = 4.0 * ts_yy_xyy[i] * gfe2_0 + 2.0 * gr_yy_xyy[i] * gfe_0 + 2.0 * ts_yy_xxyy[i] * gfe_0 * gc_x[i] + gr_yy_xxyy[i] * gc_x[i];

        grr_x_yy_xxyz[i] = 4.0 * ts_yy_xyz[i] * gfe2_0 + 2.0 * gr_yy_xyz[i] * gfe_0 + 2.0 * ts_yy_xxyz[i] * gfe_0 * gc_x[i] + gr_yy_xxyz[i] * gc_x[i];

        grr_x_yy_xxzz[i] = 4.0 * ts_yy_xzz[i] * gfe2_0 + 2.0 * gr_yy_xzz[i] * gfe_0 + 2.0 * ts_yy_xxzz[i] * gfe_0 * gc_x[i] + gr_yy_xxzz[i] * gc_x[i];

        grr_x_yy_xyyy[i] = 2.0 * ts_yy_yyy[i] * gfe2_0 + gr_yy_yyy[i] * gfe_0 + 2.0 * ts_yy_xyyy[i] * gfe_0 * gc_x[i] + gr_yy_xyyy[i] * gc_x[i];

        grr_x_yy_xyyz[i] = 2.0 * ts_yy_yyz[i] * gfe2_0 + gr_yy_yyz[i] * gfe_0 + 2.0 * ts_yy_xyyz[i] * gfe_0 * gc_x[i] + gr_yy_xyyz[i] * gc_x[i];

        grr_x_yy_xyzz[i] = 2.0 * ts_yy_yzz[i] * gfe2_0 + gr_yy_yzz[i] * gfe_0 + 2.0 * ts_yy_xyzz[i] * gfe_0 * gc_x[i] + gr_yy_xyzz[i] * gc_x[i];

        grr_x_yy_xzzz[i] = 2.0 * ts_yy_zzz[i] * gfe2_0 + gr_yy_zzz[i] * gfe_0 + 2.0 * ts_yy_xzzz[i] * gfe_0 * gc_x[i] + gr_yy_xzzz[i] * gc_x[i];

        grr_x_yy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * gc_x[i] + gr_yy_yyyy[i] * gc_x[i];

        grr_x_yy_yyyz[i] = 2.0 * ts_yy_yyyz[i] * gfe_0 * gc_x[i] + gr_yy_yyyz[i] * gc_x[i];

        grr_x_yy_yyzz[i] = 2.0 * ts_yy_yyzz[i] * gfe_0 * gc_x[i] + gr_yy_yyzz[i] * gc_x[i];

        grr_x_yy_yzzz[i] = 2.0 * ts_yy_yzzz[i] * gfe_0 * gc_x[i] + gr_yy_yzzz[i] * gc_x[i];

        grr_x_yy_zzzz[i] = 2.0 * ts_yy_zzzz[i] * gfe_0 * gc_x[i] + gr_yy_zzzz[i] * gc_x[i];
    }

    // Set up 60-75 components of targeted buffer : DG

    auto grr_x_yz_xxxx = pbuffer.data(idx_gr_dg + 60);

    auto grr_x_yz_xxxy = pbuffer.data(idx_gr_dg + 61);

    auto grr_x_yz_xxxz = pbuffer.data(idx_gr_dg + 62);

    auto grr_x_yz_xxyy = pbuffer.data(idx_gr_dg + 63);

    auto grr_x_yz_xxyz = pbuffer.data(idx_gr_dg + 64);

    auto grr_x_yz_xxzz = pbuffer.data(idx_gr_dg + 65);

    auto grr_x_yz_xyyy = pbuffer.data(idx_gr_dg + 66);

    auto grr_x_yz_xyyz = pbuffer.data(idx_gr_dg + 67);

    auto grr_x_yz_xyzz = pbuffer.data(idx_gr_dg + 68);

    auto grr_x_yz_xzzz = pbuffer.data(idx_gr_dg + 69);

    auto grr_x_yz_yyyy = pbuffer.data(idx_gr_dg + 70);

    auto grr_x_yz_yyyz = pbuffer.data(idx_gr_dg + 71);

    auto grr_x_yz_yyzz = pbuffer.data(idx_gr_dg + 72);

    auto grr_x_yz_yzzz = pbuffer.data(idx_gr_dg + 73);

    auto grr_x_yz_zzzz = pbuffer.data(idx_gr_dg + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xxx, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxy, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxz, gr_yz_xxzz, gr_yz_xyy, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyz, gr_yz_xyzz, gr_yz_xzz, gr_yz_xzzz, gr_yz_yyy, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyz, gr_yz_yyzz, gr_yz_yzz, gr_yz_yzzz, gr_yz_zzz, gr_yz_zzzz, grr_x_yz_xxxx, grr_x_yz_xxxy, grr_x_yz_xxxz, grr_x_yz_xxyy, grr_x_yz_xxyz, grr_x_yz_xxzz, grr_x_yz_xyyy, grr_x_yz_xyyz, grr_x_yz_xyzz, grr_x_yz_xzzz, grr_x_yz_yyyy, grr_x_yz_yyyz, grr_x_yz_yyzz, grr_x_yz_yzzz, grr_x_yz_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yz_xxxx[i] = 8.0 * ts_yz_xxx[i] * gfe2_0 + 4.0 * gr_yz_xxx[i] * gfe_0 + 2.0 * ts_yz_xxxx[i] * gfe_0 * gc_x[i] + gr_yz_xxxx[i] * gc_x[i];

        grr_x_yz_xxxy[i] = 6.0 * ts_yz_xxy[i] * gfe2_0 + 3.0 * gr_yz_xxy[i] * gfe_0 + 2.0 * ts_yz_xxxy[i] * gfe_0 * gc_x[i] + gr_yz_xxxy[i] * gc_x[i];

        grr_x_yz_xxxz[i] = 6.0 * ts_yz_xxz[i] * gfe2_0 + 3.0 * gr_yz_xxz[i] * gfe_0 + 2.0 * ts_yz_xxxz[i] * gfe_0 * gc_x[i] + gr_yz_xxxz[i] * gc_x[i];

        grr_x_yz_xxyy[i] = 4.0 * ts_yz_xyy[i] * gfe2_0 + 2.0 * gr_yz_xyy[i] * gfe_0 + 2.0 * ts_yz_xxyy[i] * gfe_0 * gc_x[i] + gr_yz_xxyy[i] * gc_x[i];

        grr_x_yz_xxyz[i] = 4.0 * ts_yz_xyz[i] * gfe2_0 + 2.0 * gr_yz_xyz[i] * gfe_0 + 2.0 * ts_yz_xxyz[i] * gfe_0 * gc_x[i] + gr_yz_xxyz[i] * gc_x[i];

        grr_x_yz_xxzz[i] = 4.0 * ts_yz_xzz[i] * gfe2_0 + 2.0 * gr_yz_xzz[i] * gfe_0 + 2.0 * ts_yz_xxzz[i] * gfe_0 * gc_x[i] + gr_yz_xxzz[i] * gc_x[i];

        grr_x_yz_xyyy[i] = 2.0 * ts_yz_yyy[i] * gfe2_0 + gr_yz_yyy[i] * gfe_0 + 2.0 * ts_yz_xyyy[i] * gfe_0 * gc_x[i] + gr_yz_xyyy[i] * gc_x[i];

        grr_x_yz_xyyz[i] = 2.0 * ts_yz_yyz[i] * gfe2_0 + gr_yz_yyz[i] * gfe_0 + 2.0 * ts_yz_xyyz[i] * gfe_0 * gc_x[i] + gr_yz_xyyz[i] * gc_x[i];

        grr_x_yz_xyzz[i] = 2.0 * ts_yz_yzz[i] * gfe2_0 + gr_yz_yzz[i] * gfe_0 + 2.0 * ts_yz_xyzz[i] * gfe_0 * gc_x[i] + gr_yz_xyzz[i] * gc_x[i];

        grr_x_yz_xzzz[i] = 2.0 * ts_yz_zzz[i] * gfe2_0 + gr_yz_zzz[i] * gfe_0 + 2.0 * ts_yz_xzzz[i] * gfe_0 * gc_x[i] + gr_yz_xzzz[i] * gc_x[i];

        grr_x_yz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe_0 * gc_x[i] + gr_yz_yyyy[i] * gc_x[i];

        grr_x_yz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe_0 * gc_x[i] + gr_yz_yyyz[i] * gc_x[i];

        grr_x_yz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe_0 * gc_x[i] + gr_yz_yyzz[i] * gc_x[i];

        grr_x_yz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe_0 * gc_x[i] + gr_yz_yzzz[i] * gc_x[i];

        grr_x_yz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe_0 * gc_x[i] + gr_yz_zzzz[i] * gc_x[i];
    }

    // Set up 75-90 components of targeted buffer : DG

    auto grr_x_zz_xxxx = pbuffer.data(idx_gr_dg + 75);

    auto grr_x_zz_xxxy = pbuffer.data(idx_gr_dg + 76);

    auto grr_x_zz_xxxz = pbuffer.data(idx_gr_dg + 77);

    auto grr_x_zz_xxyy = pbuffer.data(idx_gr_dg + 78);

    auto grr_x_zz_xxyz = pbuffer.data(idx_gr_dg + 79);

    auto grr_x_zz_xxzz = pbuffer.data(idx_gr_dg + 80);

    auto grr_x_zz_xyyy = pbuffer.data(idx_gr_dg + 81);

    auto grr_x_zz_xyyz = pbuffer.data(idx_gr_dg + 82);

    auto grr_x_zz_xyzz = pbuffer.data(idx_gr_dg + 83);

    auto grr_x_zz_xzzz = pbuffer.data(idx_gr_dg + 84);

    auto grr_x_zz_yyyy = pbuffer.data(idx_gr_dg + 85);

    auto grr_x_zz_yyyz = pbuffer.data(idx_gr_dg + 86);

    auto grr_x_zz_yyzz = pbuffer.data(idx_gr_dg + 87);

    auto grr_x_zz_yzzz = pbuffer.data(idx_gr_dg + 88);

    auto grr_x_zz_zzzz = pbuffer.data(idx_gr_dg + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xxx, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxy, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxz, gr_zz_xxzz, gr_zz_xyy, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyz, gr_zz_xyzz, gr_zz_xzz, gr_zz_xzzz, gr_zz_yyy, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyz, gr_zz_yyzz, gr_zz_yzz, gr_zz_yzzz, gr_zz_zzz, gr_zz_zzzz, grr_x_zz_xxxx, grr_x_zz_xxxy, grr_x_zz_xxxz, grr_x_zz_xxyy, grr_x_zz_xxyz, grr_x_zz_xxzz, grr_x_zz_xyyy, grr_x_zz_xyyz, grr_x_zz_xyzz, grr_x_zz_xzzz, grr_x_zz_yyyy, grr_x_zz_yyyz, grr_x_zz_yyzz, grr_x_zz_yzzz, grr_x_zz_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zz_xxxx[i] = 8.0 * ts_zz_xxx[i] * gfe2_0 + 4.0 * gr_zz_xxx[i] * gfe_0 + 2.0 * ts_zz_xxxx[i] * gfe_0 * gc_x[i] + gr_zz_xxxx[i] * gc_x[i];

        grr_x_zz_xxxy[i] = 6.0 * ts_zz_xxy[i] * gfe2_0 + 3.0 * gr_zz_xxy[i] * gfe_0 + 2.0 * ts_zz_xxxy[i] * gfe_0 * gc_x[i] + gr_zz_xxxy[i] * gc_x[i];

        grr_x_zz_xxxz[i] = 6.0 * ts_zz_xxz[i] * gfe2_0 + 3.0 * gr_zz_xxz[i] * gfe_0 + 2.0 * ts_zz_xxxz[i] * gfe_0 * gc_x[i] + gr_zz_xxxz[i] * gc_x[i];

        grr_x_zz_xxyy[i] = 4.0 * ts_zz_xyy[i] * gfe2_0 + 2.0 * gr_zz_xyy[i] * gfe_0 + 2.0 * ts_zz_xxyy[i] * gfe_0 * gc_x[i] + gr_zz_xxyy[i] * gc_x[i];

        grr_x_zz_xxyz[i] = 4.0 * ts_zz_xyz[i] * gfe2_0 + 2.0 * gr_zz_xyz[i] * gfe_0 + 2.0 * ts_zz_xxyz[i] * gfe_0 * gc_x[i] + gr_zz_xxyz[i] * gc_x[i];

        grr_x_zz_xxzz[i] = 4.0 * ts_zz_xzz[i] * gfe2_0 + 2.0 * gr_zz_xzz[i] * gfe_0 + 2.0 * ts_zz_xxzz[i] * gfe_0 * gc_x[i] + gr_zz_xxzz[i] * gc_x[i];

        grr_x_zz_xyyy[i] = 2.0 * ts_zz_yyy[i] * gfe2_0 + gr_zz_yyy[i] * gfe_0 + 2.0 * ts_zz_xyyy[i] * gfe_0 * gc_x[i] + gr_zz_xyyy[i] * gc_x[i];

        grr_x_zz_xyyz[i] = 2.0 * ts_zz_yyz[i] * gfe2_0 + gr_zz_yyz[i] * gfe_0 + 2.0 * ts_zz_xyyz[i] * gfe_0 * gc_x[i] + gr_zz_xyyz[i] * gc_x[i];

        grr_x_zz_xyzz[i] = 2.0 * ts_zz_yzz[i] * gfe2_0 + gr_zz_yzz[i] * gfe_0 + 2.0 * ts_zz_xyzz[i] * gfe_0 * gc_x[i] + gr_zz_xyzz[i] * gc_x[i];

        grr_x_zz_xzzz[i] = 2.0 * ts_zz_zzz[i] * gfe2_0 + gr_zz_zzz[i] * gfe_0 + 2.0 * ts_zz_xzzz[i] * gfe_0 * gc_x[i] + gr_zz_xzzz[i] * gc_x[i];

        grr_x_zz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe_0 * gc_x[i] + gr_zz_yyyy[i] * gc_x[i];

        grr_x_zz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe_0 * gc_x[i] + gr_zz_yyyz[i] * gc_x[i];

        grr_x_zz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe_0 * gc_x[i] + gr_zz_yyzz[i] * gc_x[i];

        grr_x_zz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe_0 * gc_x[i] + gr_zz_yzzz[i] * gc_x[i];

        grr_x_zz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * gc_x[i] + gr_zz_zzzz[i] * gc_x[i];
    }

    // Set up 90-105 components of targeted buffer : DG

    auto grr_y_xx_xxxx = pbuffer.data(idx_gr_dg + 90);

    auto grr_y_xx_xxxy = pbuffer.data(idx_gr_dg + 91);

    auto grr_y_xx_xxxz = pbuffer.data(idx_gr_dg + 92);

    auto grr_y_xx_xxyy = pbuffer.data(idx_gr_dg + 93);

    auto grr_y_xx_xxyz = pbuffer.data(idx_gr_dg + 94);

    auto grr_y_xx_xxzz = pbuffer.data(idx_gr_dg + 95);

    auto grr_y_xx_xyyy = pbuffer.data(idx_gr_dg + 96);

    auto grr_y_xx_xyyz = pbuffer.data(idx_gr_dg + 97);

    auto grr_y_xx_xyzz = pbuffer.data(idx_gr_dg + 98);

    auto grr_y_xx_xzzz = pbuffer.data(idx_gr_dg + 99);

    auto grr_y_xx_yyyy = pbuffer.data(idx_gr_dg + 100);

    auto grr_y_xx_yyyz = pbuffer.data(idx_gr_dg + 101);

    auto grr_y_xx_yyzz = pbuffer.data(idx_gr_dg + 102);

    auto grr_y_xx_yzzz = pbuffer.data(idx_gr_dg + 103);

    auto grr_y_xx_zzzz = pbuffer.data(idx_gr_dg + 104);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxx, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxy, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxz, gr_xx_xxzz, gr_xx_xyy, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyz, gr_xx_xyzz, gr_xx_xzz, gr_xx_xzzz, gr_xx_yyy, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyz, gr_xx_yyzz, gr_xx_yzz, gr_xx_yzzz, gr_xx_zzz, gr_xx_zzzz, grr_y_xx_xxxx, grr_y_xx_xxxy, grr_y_xx_xxxz, grr_y_xx_xxyy, grr_y_xx_xxyz, grr_y_xx_xxzz, grr_y_xx_xyyy, grr_y_xx_xyyz, grr_y_xx_xyzz, grr_y_xx_xzzz, grr_y_xx_yyyy, grr_y_xx_yyyz, grr_y_xx_yyzz, grr_y_xx_yzzz, grr_y_xx_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xx_xxxx[i] = 2.0 * ts_xx_xxxx[i] * gfe_0 * gc_y[i] + gr_xx_xxxx[i] * gc_y[i];

        grr_y_xx_xxxy[i] = 2.0 * ts_xx_xxx[i] * gfe2_0 + gr_xx_xxx[i] * gfe_0 + 2.0 * ts_xx_xxxy[i] * gfe_0 * gc_y[i] + gr_xx_xxxy[i] * gc_y[i];

        grr_y_xx_xxxz[i] = 2.0 * ts_xx_xxxz[i] * gfe_0 * gc_y[i] + gr_xx_xxxz[i] * gc_y[i];

        grr_y_xx_xxyy[i] = 4.0 * ts_xx_xxy[i] * gfe2_0 + 2.0 * gr_xx_xxy[i] * gfe_0 + 2.0 * ts_xx_xxyy[i] * gfe_0 * gc_y[i] + gr_xx_xxyy[i] * gc_y[i];

        grr_y_xx_xxyz[i] = 2.0 * ts_xx_xxz[i] * gfe2_0 + gr_xx_xxz[i] * gfe_0 + 2.0 * ts_xx_xxyz[i] * gfe_0 * gc_y[i] + gr_xx_xxyz[i] * gc_y[i];

        grr_y_xx_xxzz[i] = 2.0 * ts_xx_xxzz[i] * gfe_0 * gc_y[i] + gr_xx_xxzz[i] * gc_y[i];

        grr_y_xx_xyyy[i] = 6.0 * ts_xx_xyy[i] * gfe2_0 + 3.0 * gr_xx_xyy[i] * gfe_0 + 2.0 * ts_xx_xyyy[i] * gfe_0 * gc_y[i] + gr_xx_xyyy[i] * gc_y[i];

        grr_y_xx_xyyz[i] = 4.0 * ts_xx_xyz[i] * gfe2_0 + 2.0 * gr_xx_xyz[i] * gfe_0 + 2.0 * ts_xx_xyyz[i] * gfe_0 * gc_y[i] + gr_xx_xyyz[i] * gc_y[i];

        grr_y_xx_xyzz[i] = 2.0 * ts_xx_xzz[i] * gfe2_0 + gr_xx_xzz[i] * gfe_0 + 2.0 * ts_xx_xyzz[i] * gfe_0 * gc_y[i] + gr_xx_xyzz[i] * gc_y[i];

        grr_y_xx_xzzz[i] = 2.0 * ts_xx_xzzz[i] * gfe_0 * gc_y[i] + gr_xx_xzzz[i] * gc_y[i];

        grr_y_xx_yyyy[i] = 8.0 * ts_xx_yyy[i] * gfe2_0 + 4.0 * gr_xx_yyy[i] * gfe_0 + 2.0 * ts_xx_yyyy[i] * gfe_0 * gc_y[i] + gr_xx_yyyy[i] * gc_y[i];

        grr_y_xx_yyyz[i] = 6.0 * ts_xx_yyz[i] * gfe2_0 + 3.0 * gr_xx_yyz[i] * gfe_0 + 2.0 * ts_xx_yyyz[i] * gfe_0 * gc_y[i] + gr_xx_yyyz[i] * gc_y[i];

        grr_y_xx_yyzz[i] = 4.0 * ts_xx_yzz[i] * gfe2_0 + 2.0 * gr_xx_yzz[i] * gfe_0 + 2.0 * ts_xx_yyzz[i] * gfe_0 * gc_y[i] + gr_xx_yyzz[i] * gc_y[i];

        grr_y_xx_yzzz[i] = 2.0 * ts_xx_zzz[i] * gfe2_0 + gr_xx_zzz[i] * gfe_0 + 2.0 * ts_xx_yzzz[i] * gfe_0 * gc_y[i] + gr_xx_yzzz[i] * gc_y[i];

        grr_y_xx_zzzz[i] = 2.0 * ts_xx_zzzz[i] * gfe_0 * gc_y[i] + gr_xx_zzzz[i] * gc_y[i];
    }

    // Set up 105-120 components of targeted buffer : DG

    auto grr_y_xy_xxxx = pbuffer.data(idx_gr_dg + 105);

    auto grr_y_xy_xxxy = pbuffer.data(idx_gr_dg + 106);

    auto grr_y_xy_xxxz = pbuffer.data(idx_gr_dg + 107);

    auto grr_y_xy_xxyy = pbuffer.data(idx_gr_dg + 108);

    auto grr_y_xy_xxyz = pbuffer.data(idx_gr_dg + 109);

    auto grr_y_xy_xxzz = pbuffer.data(idx_gr_dg + 110);

    auto grr_y_xy_xyyy = pbuffer.data(idx_gr_dg + 111);

    auto grr_y_xy_xyyz = pbuffer.data(idx_gr_dg + 112);

    auto grr_y_xy_xyzz = pbuffer.data(idx_gr_dg + 113);

    auto grr_y_xy_xzzz = pbuffer.data(idx_gr_dg + 114);

    auto grr_y_xy_yyyy = pbuffer.data(idx_gr_dg + 115);

    auto grr_y_xy_yyyz = pbuffer.data(idx_gr_dg + 116);

    auto grr_y_xy_yyzz = pbuffer.data(idx_gr_dg + 117);

    auto grr_y_xy_yzzz = pbuffer.data(idx_gr_dg + 118);

    auto grr_y_xy_zzzz = pbuffer.data(idx_gr_dg + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxyy, gr_x_xxyz, gr_x_xxzz, gr_x_xyyy, gr_x_xyyz, gr_x_xyzz, gr_x_xzzz, gr_x_yyyy, gr_x_yyyz, gr_x_yyzz, gr_x_yzzz, gr_x_zzzz, gr_xy_xxx, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxy, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxz, gr_xy_xxzz, gr_xy_xyy, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyz, gr_xy_xyzz, gr_xy_xzz, gr_xy_xzzz, gr_xy_yyy, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyz, gr_xy_yyzz, gr_xy_yzz, gr_xy_yzzz, gr_xy_zzz, gr_xy_zzzz, grr_y_xy_xxxx, grr_y_xy_xxxy, grr_y_xy_xxxz, grr_y_xy_xxyy, grr_y_xy_xxyz, grr_y_xy_xxzz, grr_y_xy_xyyy, grr_y_xy_xyyz, grr_y_xy_xyzz, grr_y_xy_xzzz, grr_y_xy_yyyy, grr_y_xy_yyyz, grr_y_xy_yyzz, grr_y_xy_yzzz, grr_y_xy_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xy_xxxx[i] = 2.0 * ts_x_xxxx[i] * gfe2_0 + gr_x_xxxx[i] * gfe_0 + 2.0 * ts_xy_xxxx[i] * gfe_0 * gc_y[i] + gr_xy_xxxx[i] * gc_y[i];

        grr_y_xy_xxxy[i] = 2.0 * ts_x_xxxy[i] * gfe2_0 + gr_x_xxxy[i] * gfe_0 + 2.0 * ts_xy_xxx[i] * gfe2_0 + gr_xy_xxx[i] * gfe_0 + 2.0 * ts_xy_xxxy[i] * gfe_0 * gc_y[i] + gr_xy_xxxy[i] * gc_y[i];

        grr_y_xy_xxxz[i] = 2.0 * ts_x_xxxz[i] * gfe2_0 + gr_x_xxxz[i] * gfe_0 + 2.0 * ts_xy_xxxz[i] * gfe_0 * gc_y[i] + gr_xy_xxxz[i] * gc_y[i];

        grr_y_xy_xxyy[i] = 2.0 * ts_x_xxyy[i] * gfe2_0 + gr_x_xxyy[i] * gfe_0 + 4.0 * ts_xy_xxy[i] * gfe2_0 + 2.0 * gr_xy_xxy[i] * gfe_0 + 2.0 * ts_xy_xxyy[i] * gfe_0 * gc_y[i] + gr_xy_xxyy[i] * gc_y[i];

        grr_y_xy_xxyz[i] = 2.0 * ts_x_xxyz[i] * gfe2_0 + gr_x_xxyz[i] * gfe_0 + 2.0 * ts_xy_xxz[i] * gfe2_0 + gr_xy_xxz[i] * gfe_0 + 2.0 * ts_xy_xxyz[i] * gfe_0 * gc_y[i] + gr_xy_xxyz[i] * gc_y[i];

        grr_y_xy_xxzz[i] = 2.0 * ts_x_xxzz[i] * gfe2_0 + gr_x_xxzz[i] * gfe_0 + 2.0 * ts_xy_xxzz[i] * gfe_0 * gc_y[i] + gr_xy_xxzz[i] * gc_y[i];

        grr_y_xy_xyyy[i] = 2.0 * ts_x_xyyy[i] * gfe2_0 + gr_x_xyyy[i] * gfe_0 + 6.0 * ts_xy_xyy[i] * gfe2_0 + 3.0 * gr_xy_xyy[i] * gfe_0 + 2.0 * ts_xy_xyyy[i] * gfe_0 * gc_y[i] + gr_xy_xyyy[i] * gc_y[i];

        grr_y_xy_xyyz[i] = 2.0 * ts_x_xyyz[i] * gfe2_0 + gr_x_xyyz[i] * gfe_0 + 4.0 * ts_xy_xyz[i] * gfe2_0 + 2.0 * gr_xy_xyz[i] * gfe_0 + 2.0 * ts_xy_xyyz[i] * gfe_0 * gc_y[i] + gr_xy_xyyz[i] * gc_y[i];

        grr_y_xy_xyzz[i] = 2.0 * ts_x_xyzz[i] * gfe2_0 + gr_x_xyzz[i] * gfe_0 + 2.0 * ts_xy_xzz[i] * gfe2_0 + gr_xy_xzz[i] * gfe_0 + 2.0 * ts_xy_xyzz[i] * gfe_0 * gc_y[i] + gr_xy_xyzz[i] * gc_y[i];

        grr_y_xy_xzzz[i] = 2.0 * ts_x_xzzz[i] * gfe2_0 + gr_x_xzzz[i] * gfe_0 + 2.0 * ts_xy_xzzz[i] * gfe_0 * gc_y[i] + gr_xy_xzzz[i] * gc_y[i];

        grr_y_xy_yyyy[i] = 2.0 * ts_x_yyyy[i] * gfe2_0 + gr_x_yyyy[i] * gfe_0 + 8.0 * ts_xy_yyy[i] * gfe2_0 + 4.0 * gr_xy_yyy[i] * gfe_0 + 2.0 * ts_xy_yyyy[i] * gfe_0 * gc_y[i] + gr_xy_yyyy[i] * gc_y[i];

        grr_y_xy_yyyz[i] = 2.0 * ts_x_yyyz[i] * gfe2_0 + gr_x_yyyz[i] * gfe_0 + 6.0 * ts_xy_yyz[i] * gfe2_0 + 3.0 * gr_xy_yyz[i] * gfe_0 + 2.0 * ts_xy_yyyz[i] * gfe_0 * gc_y[i] + gr_xy_yyyz[i] * gc_y[i];

        grr_y_xy_yyzz[i] = 2.0 * ts_x_yyzz[i] * gfe2_0 + gr_x_yyzz[i] * gfe_0 + 4.0 * ts_xy_yzz[i] * gfe2_0 + 2.0 * gr_xy_yzz[i] * gfe_0 + 2.0 * ts_xy_yyzz[i] * gfe_0 * gc_y[i] + gr_xy_yyzz[i] * gc_y[i];

        grr_y_xy_yzzz[i] = 2.0 * ts_x_yzzz[i] * gfe2_0 + gr_x_yzzz[i] * gfe_0 + 2.0 * ts_xy_zzz[i] * gfe2_0 + gr_xy_zzz[i] * gfe_0 + 2.0 * ts_xy_yzzz[i] * gfe_0 * gc_y[i] + gr_xy_yzzz[i] * gc_y[i];

        grr_y_xy_zzzz[i] = 2.0 * ts_x_zzzz[i] * gfe2_0 + gr_x_zzzz[i] * gfe_0 + 2.0 * ts_xy_zzzz[i] * gfe_0 * gc_y[i] + gr_xy_zzzz[i] * gc_y[i];
    }

    // Set up 120-135 components of targeted buffer : DG

    auto grr_y_xz_xxxx = pbuffer.data(idx_gr_dg + 120);

    auto grr_y_xz_xxxy = pbuffer.data(idx_gr_dg + 121);

    auto grr_y_xz_xxxz = pbuffer.data(idx_gr_dg + 122);

    auto grr_y_xz_xxyy = pbuffer.data(idx_gr_dg + 123);

    auto grr_y_xz_xxyz = pbuffer.data(idx_gr_dg + 124);

    auto grr_y_xz_xxzz = pbuffer.data(idx_gr_dg + 125);

    auto grr_y_xz_xyyy = pbuffer.data(idx_gr_dg + 126);

    auto grr_y_xz_xyyz = pbuffer.data(idx_gr_dg + 127);

    auto grr_y_xz_xyzz = pbuffer.data(idx_gr_dg + 128);

    auto grr_y_xz_xzzz = pbuffer.data(idx_gr_dg + 129);

    auto grr_y_xz_yyyy = pbuffer.data(idx_gr_dg + 130);

    auto grr_y_xz_yyyz = pbuffer.data(idx_gr_dg + 131);

    auto grr_y_xz_yyzz = pbuffer.data(idx_gr_dg + 132);

    auto grr_y_xz_yzzz = pbuffer.data(idx_gr_dg + 133);

    auto grr_y_xz_zzzz = pbuffer.data(idx_gr_dg + 134);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xxx, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxy, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxz, gr_xz_xxzz, gr_xz_xyy, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyz, gr_xz_xyzz, gr_xz_xzz, gr_xz_xzzz, gr_xz_yyy, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyz, gr_xz_yyzz, gr_xz_yzz, gr_xz_yzzz, gr_xz_zzz, gr_xz_zzzz, grr_y_xz_xxxx, grr_y_xz_xxxy, grr_y_xz_xxxz, grr_y_xz_xxyy, grr_y_xz_xxyz, grr_y_xz_xxzz, grr_y_xz_xyyy, grr_y_xz_xyyz, grr_y_xz_xyzz, grr_y_xz_xzzz, grr_y_xz_yyyy, grr_y_xz_yyyz, grr_y_xz_yyzz, grr_y_xz_yzzz, grr_y_xz_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xz_xxxx[i] = 2.0 * ts_xz_xxxx[i] * gfe_0 * gc_y[i] + gr_xz_xxxx[i] * gc_y[i];

        grr_y_xz_xxxy[i] = 2.0 * ts_xz_xxx[i] * gfe2_0 + gr_xz_xxx[i] * gfe_0 + 2.0 * ts_xz_xxxy[i] * gfe_0 * gc_y[i] + gr_xz_xxxy[i] * gc_y[i];

        grr_y_xz_xxxz[i] = 2.0 * ts_xz_xxxz[i] * gfe_0 * gc_y[i] + gr_xz_xxxz[i] * gc_y[i];

        grr_y_xz_xxyy[i] = 4.0 * ts_xz_xxy[i] * gfe2_0 + 2.0 * gr_xz_xxy[i] * gfe_0 + 2.0 * ts_xz_xxyy[i] * gfe_0 * gc_y[i] + gr_xz_xxyy[i] * gc_y[i];

        grr_y_xz_xxyz[i] = 2.0 * ts_xz_xxz[i] * gfe2_0 + gr_xz_xxz[i] * gfe_0 + 2.0 * ts_xz_xxyz[i] * gfe_0 * gc_y[i] + gr_xz_xxyz[i] * gc_y[i];

        grr_y_xz_xxzz[i] = 2.0 * ts_xz_xxzz[i] * gfe_0 * gc_y[i] + gr_xz_xxzz[i] * gc_y[i];

        grr_y_xz_xyyy[i] = 6.0 * ts_xz_xyy[i] * gfe2_0 + 3.0 * gr_xz_xyy[i] * gfe_0 + 2.0 * ts_xz_xyyy[i] * gfe_0 * gc_y[i] + gr_xz_xyyy[i] * gc_y[i];

        grr_y_xz_xyyz[i] = 4.0 * ts_xz_xyz[i] * gfe2_0 + 2.0 * gr_xz_xyz[i] * gfe_0 + 2.0 * ts_xz_xyyz[i] * gfe_0 * gc_y[i] + gr_xz_xyyz[i] * gc_y[i];

        grr_y_xz_xyzz[i] = 2.0 * ts_xz_xzz[i] * gfe2_0 + gr_xz_xzz[i] * gfe_0 + 2.0 * ts_xz_xyzz[i] * gfe_0 * gc_y[i] + gr_xz_xyzz[i] * gc_y[i];

        grr_y_xz_xzzz[i] = 2.0 * ts_xz_xzzz[i] * gfe_0 * gc_y[i] + gr_xz_xzzz[i] * gc_y[i];

        grr_y_xz_yyyy[i] = 8.0 * ts_xz_yyy[i] * gfe2_0 + 4.0 * gr_xz_yyy[i] * gfe_0 + 2.0 * ts_xz_yyyy[i] * gfe_0 * gc_y[i] + gr_xz_yyyy[i] * gc_y[i];

        grr_y_xz_yyyz[i] = 6.0 * ts_xz_yyz[i] * gfe2_0 + 3.0 * gr_xz_yyz[i] * gfe_0 + 2.0 * ts_xz_yyyz[i] * gfe_0 * gc_y[i] + gr_xz_yyyz[i] * gc_y[i];

        grr_y_xz_yyzz[i] = 4.0 * ts_xz_yzz[i] * gfe2_0 + 2.0 * gr_xz_yzz[i] * gfe_0 + 2.0 * ts_xz_yyzz[i] * gfe_0 * gc_y[i] + gr_xz_yyzz[i] * gc_y[i];

        grr_y_xz_yzzz[i] = 2.0 * ts_xz_zzz[i] * gfe2_0 + gr_xz_zzz[i] * gfe_0 + 2.0 * ts_xz_yzzz[i] * gfe_0 * gc_y[i] + gr_xz_yzzz[i] * gc_y[i];

        grr_y_xz_zzzz[i] = 2.0 * ts_xz_zzzz[i] * gfe_0 * gc_y[i] + gr_xz_zzzz[i] * gc_y[i];
    }

    // Set up 135-150 components of targeted buffer : DG

    auto grr_y_yy_xxxx = pbuffer.data(idx_gr_dg + 135);

    auto grr_y_yy_xxxy = pbuffer.data(idx_gr_dg + 136);

    auto grr_y_yy_xxxz = pbuffer.data(idx_gr_dg + 137);

    auto grr_y_yy_xxyy = pbuffer.data(idx_gr_dg + 138);

    auto grr_y_yy_xxyz = pbuffer.data(idx_gr_dg + 139);

    auto grr_y_yy_xxzz = pbuffer.data(idx_gr_dg + 140);

    auto grr_y_yy_xyyy = pbuffer.data(idx_gr_dg + 141);

    auto grr_y_yy_xyyz = pbuffer.data(idx_gr_dg + 142);

    auto grr_y_yy_xyzz = pbuffer.data(idx_gr_dg + 143);

    auto grr_y_yy_xzzz = pbuffer.data(idx_gr_dg + 144);

    auto grr_y_yy_yyyy = pbuffer.data(idx_gr_dg + 145);

    auto grr_y_yy_yyyz = pbuffer.data(idx_gr_dg + 146);

    auto grr_y_yy_yyzz = pbuffer.data(idx_gr_dg + 147);

    auto grr_y_yy_yzzz = pbuffer.data(idx_gr_dg + 148);

    auto grr_y_yy_zzzz = pbuffer.data(idx_gr_dg + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxyy, gr_y_xxyz, gr_y_xxzz, gr_y_xyyy, gr_y_xyyz, gr_y_xyzz, gr_y_xzzz, gr_y_yyyy, gr_y_yyyz, gr_y_yyzz, gr_y_yzzz, gr_y_zzzz, gr_yy_xxx, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxy, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxz, gr_yy_xxzz, gr_yy_xyy, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyz, gr_yy_xyzz, gr_yy_xzz, gr_yy_xzzz, gr_yy_yyy, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyz, gr_yy_yyzz, gr_yy_yzz, gr_yy_yzzz, gr_yy_zzz, gr_yy_zzzz, grr_y_yy_xxxx, grr_y_yy_xxxy, grr_y_yy_xxxz, grr_y_yy_xxyy, grr_y_yy_xxyz, grr_y_yy_xxzz, grr_y_yy_xyyy, grr_y_yy_xyyz, grr_y_yy_xyzz, grr_y_yy_xzzz, grr_y_yy_yyyy, grr_y_yy_yyyz, grr_y_yy_yyzz, grr_y_yy_yzzz, grr_y_yy_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yy_xxxx[i] = 4.0 * ts_y_xxxx[i] * gfe2_0 + 2.0 * gr_y_xxxx[i] * gfe_0 + 2.0 * ts_yy_xxxx[i] * gfe_0 * gc_y[i] + gr_yy_xxxx[i] * gc_y[i];

        grr_y_yy_xxxy[i] = 4.0 * ts_y_xxxy[i] * gfe2_0 + 2.0 * gr_y_xxxy[i] * gfe_0 + 2.0 * ts_yy_xxx[i] * gfe2_0 + gr_yy_xxx[i] * gfe_0 + 2.0 * ts_yy_xxxy[i] * gfe_0 * gc_y[i] + gr_yy_xxxy[i] * gc_y[i];

        grr_y_yy_xxxz[i] = 4.0 * ts_y_xxxz[i] * gfe2_0 + 2.0 * gr_y_xxxz[i] * gfe_0 + 2.0 * ts_yy_xxxz[i] * gfe_0 * gc_y[i] + gr_yy_xxxz[i] * gc_y[i];

        grr_y_yy_xxyy[i] = 4.0 * ts_y_xxyy[i] * gfe2_0 + 2.0 * gr_y_xxyy[i] * gfe_0 + 4.0 * ts_yy_xxy[i] * gfe2_0 + 2.0 * gr_yy_xxy[i] * gfe_0 + 2.0 * ts_yy_xxyy[i] * gfe_0 * gc_y[i] + gr_yy_xxyy[i] * gc_y[i];

        grr_y_yy_xxyz[i] = 4.0 * ts_y_xxyz[i] * gfe2_0 + 2.0 * gr_y_xxyz[i] * gfe_0 + 2.0 * ts_yy_xxz[i] * gfe2_0 + gr_yy_xxz[i] * gfe_0 + 2.0 * ts_yy_xxyz[i] * gfe_0 * gc_y[i] + gr_yy_xxyz[i] * gc_y[i];

        grr_y_yy_xxzz[i] = 4.0 * ts_y_xxzz[i] * gfe2_0 + 2.0 * gr_y_xxzz[i] * gfe_0 + 2.0 * ts_yy_xxzz[i] * gfe_0 * gc_y[i] + gr_yy_xxzz[i] * gc_y[i];

        grr_y_yy_xyyy[i] = 4.0 * ts_y_xyyy[i] * gfe2_0 + 2.0 * gr_y_xyyy[i] * gfe_0 + 6.0 * ts_yy_xyy[i] * gfe2_0 + 3.0 * gr_yy_xyy[i] * gfe_0 + 2.0 * ts_yy_xyyy[i] * gfe_0 * gc_y[i] + gr_yy_xyyy[i] * gc_y[i];

        grr_y_yy_xyyz[i] = 4.0 * ts_y_xyyz[i] * gfe2_0 + 2.0 * gr_y_xyyz[i] * gfe_0 + 4.0 * ts_yy_xyz[i] * gfe2_0 + 2.0 * gr_yy_xyz[i] * gfe_0 + 2.0 * ts_yy_xyyz[i] * gfe_0 * gc_y[i] + gr_yy_xyyz[i] * gc_y[i];

        grr_y_yy_xyzz[i] = 4.0 * ts_y_xyzz[i] * gfe2_0 + 2.0 * gr_y_xyzz[i] * gfe_0 + 2.0 * ts_yy_xzz[i] * gfe2_0 + gr_yy_xzz[i] * gfe_0 + 2.0 * ts_yy_xyzz[i] * gfe_0 * gc_y[i] + gr_yy_xyzz[i] * gc_y[i];

        grr_y_yy_xzzz[i] = 4.0 * ts_y_xzzz[i] * gfe2_0 + 2.0 * gr_y_xzzz[i] * gfe_0 + 2.0 * ts_yy_xzzz[i] * gfe_0 * gc_y[i] + gr_yy_xzzz[i] * gc_y[i];

        grr_y_yy_yyyy[i] = 4.0 * ts_y_yyyy[i] * gfe2_0 + 2.0 * gr_y_yyyy[i] * gfe_0 + 8.0 * ts_yy_yyy[i] * gfe2_0 + 4.0 * gr_yy_yyy[i] * gfe_0 + 2.0 * ts_yy_yyyy[i] * gfe_0 * gc_y[i] + gr_yy_yyyy[i] * gc_y[i];

        grr_y_yy_yyyz[i] = 4.0 * ts_y_yyyz[i] * gfe2_0 + 2.0 * gr_y_yyyz[i] * gfe_0 + 6.0 * ts_yy_yyz[i] * gfe2_0 + 3.0 * gr_yy_yyz[i] * gfe_0 + 2.0 * ts_yy_yyyz[i] * gfe_0 * gc_y[i] + gr_yy_yyyz[i] * gc_y[i];

        grr_y_yy_yyzz[i] = 4.0 * ts_y_yyzz[i] * gfe2_0 + 2.0 * gr_y_yyzz[i] * gfe_0 + 4.0 * ts_yy_yzz[i] * gfe2_0 + 2.0 * gr_yy_yzz[i] * gfe_0 + 2.0 * ts_yy_yyzz[i] * gfe_0 * gc_y[i] + gr_yy_yyzz[i] * gc_y[i];

        grr_y_yy_yzzz[i] = 4.0 * ts_y_yzzz[i] * gfe2_0 + 2.0 * gr_y_yzzz[i] * gfe_0 + 2.0 * ts_yy_zzz[i] * gfe2_0 + gr_yy_zzz[i] * gfe_0 + 2.0 * ts_yy_yzzz[i] * gfe_0 * gc_y[i] + gr_yy_yzzz[i] * gc_y[i];

        grr_y_yy_zzzz[i] = 4.0 * ts_y_zzzz[i] * gfe2_0 + 2.0 * gr_y_zzzz[i] * gfe_0 + 2.0 * ts_yy_zzzz[i] * gfe_0 * gc_y[i] + gr_yy_zzzz[i] * gc_y[i];
    }

    // Set up 150-165 components of targeted buffer : DG

    auto grr_y_yz_xxxx = pbuffer.data(idx_gr_dg + 150);

    auto grr_y_yz_xxxy = pbuffer.data(idx_gr_dg + 151);

    auto grr_y_yz_xxxz = pbuffer.data(idx_gr_dg + 152);

    auto grr_y_yz_xxyy = pbuffer.data(idx_gr_dg + 153);

    auto grr_y_yz_xxyz = pbuffer.data(idx_gr_dg + 154);

    auto grr_y_yz_xxzz = pbuffer.data(idx_gr_dg + 155);

    auto grr_y_yz_xyyy = pbuffer.data(idx_gr_dg + 156);

    auto grr_y_yz_xyyz = pbuffer.data(idx_gr_dg + 157);

    auto grr_y_yz_xyzz = pbuffer.data(idx_gr_dg + 158);

    auto grr_y_yz_xzzz = pbuffer.data(idx_gr_dg + 159);

    auto grr_y_yz_yyyy = pbuffer.data(idx_gr_dg + 160);

    auto grr_y_yz_yyyz = pbuffer.data(idx_gr_dg + 161);

    auto grr_y_yz_yyzz = pbuffer.data(idx_gr_dg + 162);

    auto grr_y_yz_yzzz = pbuffer.data(idx_gr_dg + 163);

    auto grr_y_yz_zzzz = pbuffer.data(idx_gr_dg + 164);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xxx, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxy, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxz, gr_yz_xxzz, gr_yz_xyy, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyz, gr_yz_xyzz, gr_yz_xzz, gr_yz_xzzz, gr_yz_yyy, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyz, gr_yz_yyzz, gr_yz_yzz, gr_yz_yzzz, gr_yz_zzz, gr_yz_zzzz, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxyy, gr_z_xxyz, gr_z_xxzz, gr_z_xyyy, gr_z_xyyz, gr_z_xyzz, gr_z_xzzz, gr_z_yyyy, gr_z_yyyz, gr_z_yyzz, gr_z_yzzz, gr_z_zzzz, grr_y_yz_xxxx, grr_y_yz_xxxy, grr_y_yz_xxxz, grr_y_yz_xxyy, grr_y_yz_xxyz, grr_y_yz_xxzz, grr_y_yz_xyyy, grr_y_yz_xyyz, grr_y_yz_xyzz, grr_y_yz_xzzz, grr_y_yz_yyyy, grr_y_yz_yyyz, grr_y_yz_yyzz, grr_y_yz_yzzz, grr_y_yz_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe2_0 + gr_z_xxxx[i] * gfe_0 + 2.0 * ts_yz_xxxx[i] * gfe_0 * gc_y[i] + gr_yz_xxxx[i] * gc_y[i];

        grr_y_yz_xxxy[i] = 2.0 * ts_z_xxxy[i] * gfe2_0 + gr_z_xxxy[i] * gfe_0 + 2.0 * ts_yz_xxx[i] * gfe2_0 + gr_yz_xxx[i] * gfe_0 + 2.0 * ts_yz_xxxy[i] * gfe_0 * gc_y[i] + gr_yz_xxxy[i] * gc_y[i];

        grr_y_yz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe2_0 + gr_z_xxxz[i] * gfe_0 + 2.0 * ts_yz_xxxz[i] * gfe_0 * gc_y[i] + gr_yz_xxxz[i] * gc_y[i];

        grr_y_yz_xxyy[i] = 2.0 * ts_z_xxyy[i] * gfe2_0 + gr_z_xxyy[i] * gfe_0 + 4.0 * ts_yz_xxy[i] * gfe2_0 + 2.0 * gr_yz_xxy[i] * gfe_0 + 2.0 * ts_yz_xxyy[i] * gfe_0 * gc_y[i] + gr_yz_xxyy[i] * gc_y[i];

        grr_y_yz_xxyz[i] = 2.0 * ts_z_xxyz[i] * gfe2_0 + gr_z_xxyz[i] * gfe_0 + 2.0 * ts_yz_xxz[i] * gfe2_0 + gr_yz_xxz[i] * gfe_0 + 2.0 * ts_yz_xxyz[i] * gfe_0 * gc_y[i] + gr_yz_xxyz[i] * gc_y[i];

        grr_y_yz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe2_0 + gr_z_xxzz[i] * gfe_0 + 2.0 * ts_yz_xxzz[i] * gfe_0 * gc_y[i] + gr_yz_xxzz[i] * gc_y[i];

        grr_y_yz_xyyy[i] = 2.0 * ts_z_xyyy[i] * gfe2_0 + gr_z_xyyy[i] * gfe_0 + 6.0 * ts_yz_xyy[i] * gfe2_0 + 3.0 * gr_yz_xyy[i] * gfe_0 + 2.0 * ts_yz_xyyy[i] * gfe_0 * gc_y[i] + gr_yz_xyyy[i] * gc_y[i];

        grr_y_yz_xyyz[i] = 2.0 * ts_z_xyyz[i] * gfe2_0 + gr_z_xyyz[i] * gfe_0 + 4.0 * ts_yz_xyz[i] * gfe2_0 + 2.0 * gr_yz_xyz[i] * gfe_0 + 2.0 * ts_yz_xyyz[i] * gfe_0 * gc_y[i] + gr_yz_xyyz[i] * gc_y[i];

        grr_y_yz_xyzz[i] = 2.0 * ts_z_xyzz[i] * gfe2_0 + gr_z_xyzz[i] * gfe_0 + 2.0 * ts_yz_xzz[i] * gfe2_0 + gr_yz_xzz[i] * gfe_0 + 2.0 * ts_yz_xyzz[i] * gfe_0 * gc_y[i] + gr_yz_xyzz[i] * gc_y[i];

        grr_y_yz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe2_0 + gr_z_xzzz[i] * gfe_0 + 2.0 * ts_yz_xzzz[i] * gfe_0 * gc_y[i] + gr_yz_xzzz[i] * gc_y[i];

        grr_y_yz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe2_0 + gr_z_yyyy[i] * gfe_0 + 8.0 * ts_yz_yyy[i] * gfe2_0 + 4.0 * gr_yz_yyy[i] * gfe_0 + 2.0 * ts_yz_yyyy[i] * gfe_0 * gc_y[i] + gr_yz_yyyy[i] * gc_y[i];

        grr_y_yz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe2_0 + gr_z_yyyz[i] * gfe_0 + 6.0 * ts_yz_yyz[i] * gfe2_0 + 3.0 * gr_yz_yyz[i] * gfe_0 + 2.0 * ts_yz_yyyz[i] * gfe_0 * gc_y[i] + gr_yz_yyyz[i] * gc_y[i];

        grr_y_yz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe2_0 + gr_z_yyzz[i] * gfe_0 + 4.0 * ts_yz_yzz[i] * gfe2_0 + 2.0 * gr_yz_yzz[i] * gfe_0 + 2.0 * ts_yz_yyzz[i] * gfe_0 * gc_y[i] + gr_yz_yyzz[i] * gc_y[i];

        grr_y_yz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe2_0 + gr_z_yzzz[i] * gfe_0 + 2.0 * ts_yz_zzz[i] * gfe2_0 + gr_yz_zzz[i] * gfe_0 + 2.0 * ts_yz_yzzz[i] * gfe_0 * gc_y[i] + gr_yz_yzzz[i] * gc_y[i];

        grr_y_yz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe2_0 + gr_z_zzzz[i] * gfe_0 + 2.0 * ts_yz_zzzz[i] * gfe_0 * gc_y[i] + gr_yz_zzzz[i] * gc_y[i];
    }

    // Set up 165-180 components of targeted buffer : DG

    auto grr_y_zz_xxxx = pbuffer.data(idx_gr_dg + 165);

    auto grr_y_zz_xxxy = pbuffer.data(idx_gr_dg + 166);

    auto grr_y_zz_xxxz = pbuffer.data(idx_gr_dg + 167);

    auto grr_y_zz_xxyy = pbuffer.data(idx_gr_dg + 168);

    auto grr_y_zz_xxyz = pbuffer.data(idx_gr_dg + 169);

    auto grr_y_zz_xxzz = pbuffer.data(idx_gr_dg + 170);

    auto grr_y_zz_xyyy = pbuffer.data(idx_gr_dg + 171);

    auto grr_y_zz_xyyz = pbuffer.data(idx_gr_dg + 172);

    auto grr_y_zz_xyzz = pbuffer.data(idx_gr_dg + 173);

    auto grr_y_zz_xzzz = pbuffer.data(idx_gr_dg + 174);

    auto grr_y_zz_yyyy = pbuffer.data(idx_gr_dg + 175);

    auto grr_y_zz_yyyz = pbuffer.data(idx_gr_dg + 176);

    auto grr_y_zz_yyzz = pbuffer.data(idx_gr_dg + 177);

    auto grr_y_zz_yzzz = pbuffer.data(idx_gr_dg + 178);

    auto grr_y_zz_zzzz = pbuffer.data(idx_gr_dg + 179);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xxx, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxy, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxz, gr_zz_xxzz, gr_zz_xyy, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyz, gr_zz_xyzz, gr_zz_xzz, gr_zz_xzzz, gr_zz_yyy, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyz, gr_zz_yyzz, gr_zz_yzz, gr_zz_yzzz, gr_zz_zzz, gr_zz_zzzz, grr_y_zz_xxxx, grr_y_zz_xxxy, grr_y_zz_xxxz, grr_y_zz_xxyy, grr_y_zz_xxyz, grr_y_zz_xxzz, grr_y_zz_xyyy, grr_y_zz_xyyz, grr_y_zz_xyzz, grr_y_zz_xzzz, grr_y_zz_yyyy, grr_y_zz_yyyz, grr_y_zz_yyzz, grr_y_zz_yzzz, grr_y_zz_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gfe_0 * gc_y[i] + gr_zz_xxxx[i] * gc_y[i];

        grr_y_zz_xxxy[i] = 2.0 * ts_zz_xxx[i] * gfe2_0 + gr_zz_xxx[i] * gfe_0 + 2.0 * ts_zz_xxxy[i] * gfe_0 * gc_y[i] + gr_zz_xxxy[i] * gc_y[i];

        grr_y_zz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gfe_0 * gc_y[i] + gr_zz_xxxz[i] * gc_y[i];

        grr_y_zz_xxyy[i] = 4.0 * ts_zz_xxy[i] * gfe2_0 + 2.0 * gr_zz_xxy[i] * gfe_0 + 2.0 * ts_zz_xxyy[i] * gfe_0 * gc_y[i] + gr_zz_xxyy[i] * gc_y[i];

        grr_y_zz_xxyz[i] = 2.0 * ts_zz_xxz[i] * gfe2_0 + gr_zz_xxz[i] * gfe_0 + 2.0 * ts_zz_xxyz[i] * gfe_0 * gc_y[i] + gr_zz_xxyz[i] * gc_y[i];

        grr_y_zz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gfe_0 * gc_y[i] + gr_zz_xxzz[i] * gc_y[i];

        grr_y_zz_xyyy[i] = 6.0 * ts_zz_xyy[i] * gfe2_0 + 3.0 * gr_zz_xyy[i] * gfe_0 + 2.0 * ts_zz_xyyy[i] * gfe_0 * gc_y[i] + gr_zz_xyyy[i] * gc_y[i];

        grr_y_zz_xyyz[i] = 4.0 * ts_zz_xyz[i] * gfe2_0 + 2.0 * gr_zz_xyz[i] * gfe_0 + 2.0 * ts_zz_xyyz[i] * gfe_0 * gc_y[i] + gr_zz_xyyz[i] * gc_y[i];

        grr_y_zz_xyzz[i] = 2.0 * ts_zz_xzz[i] * gfe2_0 + gr_zz_xzz[i] * gfe_0 + 2.0 * ts_zz_xyzz[i] * gfe_0 * gc_y[i] + gr_zz_xyzz[i] * gc_y[i];

        grr_y_zz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe_0 * gc_y[i] + gr_zz_xzzz[i] * gc_y[i];

        grr_y_zz_yyyy[i] = 8.0 * ts_zz_yyy[i] * gfe2_0 + 4.0 * gr_zz_yyy[i] * gfe_0 + 2.0 * ts_zz_yyyy[i] * gfe_0 * gc_y[i] + gr_zz_yyyy[i] * gc_y[i];

        grr_y_zz_yyyz[i] = 6.0 * ts_zz_yyz[i] * gfe2_0 + 3.0 * gr_zz_yyz[i] * gfe_0 + 2.0 * ts_zz_yyyz[i] * gfe_0 * gc_y[i] + gr_zz_yyyz[i] * gc_y[i];

        grr_y_zz_yyzz[i] = 4.0 * ts_zz_yzz[i] * gfe2_0 + 2.0 * gr_zz_yzz[i] * gfe_0 + 2.0 * ts_zz_yyzz[i] * gfe_0 * gc_y[i] + gr_zz_yyzz[i] * gc_y[i];

        grr_y_zz_yzzz[i] = 2.0 * ts_zz_zzz[i] * gfe2_0 + gr_zz_zzz[i] * gfe_0 + 2.0 * ts_zz_yzzz[i] * gfe_0 * gc_y[i] + gr_zz_yzzz[i] * gc_y[i];

        grr_y_zz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * gc_y[i] + gr_zz_zzzz[i] * gc_y[i];
    }

    // Set up 180-195 components of targeted buffer : DG

    auto grr_z_xx_xxxx = pbuffer.data(idx_gr_dg + 180);

    auto grr_z_xx_xxxy = pbuffer.data(idx_gr_dg + 181);

    auto grr_z_xx_xxxz = pbuffer.data(idx_gr_dg + 182);

    auto grr_z_xx_xxyy = pbuffer.data(idx_gr_dg + 183);

    auto grr_z_xx_xxyz = pbuffer.data(idx_gr_dg + 184);

    auto grr_z_xx_xxzz = pbuffer.data(idx_gr_dg + 185);

    auto grr_z_xx_xyyy = pbuffer.data(idx_gr_dg + 186);

    auto grr_z_xx_xyyz = pbuffer.data(idx_gr_dg + 187);

    auto grr_z_xx_xyzz = pbuffer.data(idx_gr_dg + 188);

    auto grr_z_xx_xzzz = pbuffer.data(idx_gr_dg + 189);

    auto grr_z_xx_yyyy = pbuffer.data(idx_gr_dg + 190);

    auto grr_z_xx_yyyz = pbuffer.data(idx_gr_dg + 191);

    auto grr_z_xx_yyzz = pbuffer.data(idx_gr_dg + 192);

    auto grr_z_xx_yzzz = pbuffer.data(idx_gr_dg + 193);

    auto grr_z_xx_zzzz = pbuffer.data(idx_gr_dg + 194);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxx, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxy, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxz, gr_xx_xxzz, gr_xx_xyy, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyz, gr_xx_xyzz, gr_xx_xzz, gr_xx_xzzz, gr_xx_yyy, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyz, gr_xx_yyzz, gr_xx_yzz, gr_xx_yzzz, gr_xx_zzz, gr_xx_zzzz, grr_z_xx_xxxx, grr_z_xx_xxxy, grr_z_xx_xxxz, grr_z_xx_xxyy, grr_z_xx_xxyz, grr_z_xx_xxzz, grr_z_xx_xyyy, grr_z_xx_xyyz, grr_z_xx_xyzz, grr_z_xx_xzzz, grr_z_xx_yyyy, grr_z_xx_yyyz, grr_z_xx_yyzz, grr_z_xx_yzzz, grr_z_xx_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xx_xxxx[i] = 2.0 * ts_xx_xxxx[i] * gfe_0 * gc_z[i] + gr_xx_xxxx[i] * gc_z[i];

        grr_z_xx_xxxy[i] = 2.0 * ts_xx_xxxy[i] * gfe_0 * gc_z[i] + gr_xx_xxxy[i] * gc_z[i];

        grr_z_xx_xxxz[i] = 2.0 * ts_xx_xxx[i] * gfe2_0 + gr_xx_xxx[i] * gfe_0 + 2.0 * ts_xx_xxxz[i] * gfe_0 * gc_z[i] + gr_xx_xxxz[i] * gc_z[i];

        grr_z_xx_xxyy[i] = 2.0 * ts_xx_xxyy[i] * gfe_0 * gc_z[i] + gr_xx_xxyy[i] * gc_z[i];

        grr_z_xx_xxyz[i] = 2.0 * ts_xx_xxy[i] * gfe2_0 + gr_xx_xxy[i] * gfe_0 + 2.0 * ts_xx_xxyz[i] * gfe_0 * gc_z[i] + gr_xx_xxyz[i] * gc_z[i];

        grr_z_xx_xxzz[i] = 4.0 * ts_xx_xxz[i] * gfe2_0 + 2.0 * gr_xx_xxz[i] * gfe_0 + 2.0 * ts_xx_xxzz[i] * gfe_0 * gc_z[i] + gr_xx_xxzz[i] * gc_z[i];

        grr_z_xx_xyyy[i] = 2.0 * ts_xx_xyyy[i] * gfe_0 * gc_z[i] + gr_xx_xyyy[i] * gc_z[i];

        grr_z_xx_xyyz[i] = 2.0 * ts_xx_xyy[i] * gfe2_0 + gr_xx_xyy[i] * gfe_0 + 2.0 * ts_xx_xyyz[i] * gfe_0 * gc_z[i] + gr_xx_xyyz[i] * gc_z[i];

        grr_z_xx_xyzz[i] = 4.0 * ts_xx_xyz[i] * gfe2_0 + 2.0 * gr_xx_xyz[i] * gfe_0 + 2.0 * ts_xx_xyzz[i] * gfe_0 * gc_z[i] + gr_xx_xyzz[i] * gc_z[i];

        grr_z_xx_xzzz[i] = 6.0 * ts_xx_xzz[i] * gfe2_0 + 3.0 * gr_xx_xzz[i] * gfe_0 + 2.0 * ts_xx_xzzz[i] * gfe_0 * gc_z[i] + gr_xx_xzzz[i] * gc_z[i];

        grr_z_xx_yyyy[i] = 2.0 * ts_xx_yyyy[i] * gfe_0 * gc_z[i] + gr_xx_yyyy[i] * gc_z[i];

        grr_z_xx_yyyz[i] = 2.0 * ts_xx_yyy[i] * gfe2_0 + gr_xx_yyy[i] * gfe_0 + 2.0 * ts_xx_yyyz[i] * gfe_0 * gc_z[i] + gr_xx_yyyz[i] * gc_z[i];

        grr_z_xx_yyzz[i] = 4.0 * ts_xx_yyz[i] * gfe2_0 + 2.0 * gr_xx_yyz[i] * gfe_0 + 2.0 * ts_xx_yyzz[i] * gfe_0 * gc_z[i] + gr_xx_yyzz[i] * gc_z[i];

        grr_z_xx_yzzz[i] = 6.0 * ts_xx_yzz[i] * gfe2_0 + 3.0 * gr_xx_yzz[i] * gfe_0 + 2.0 * ts_xx_yzzz[i] * gfe_0 * gc_z[i] + gr_xx_yzzz[i] * gc_z[i];

        grr_z_xx_zzzz[i] = 8.0 * ts_xx_zzz[i] * gfe2_0 + 4.0 * gr_xx_zzz[i] * gfe_0 + 2.0 * ts_xx_zzzz[i] * gfe_0 * gc_z[i] + gr_xx_zzzz[i] * gc_z[i];
    }

    // Set up 195-210 components of targeted buffer : DG

    auto grr_z_xy_xxxx = pbuffer.data(idx_gr_dg + 195);

    auto grr_z_xy_xxxy = pbuffer.data(idx_gr_dg + 196);

    auto grr_z_xy_xxxz = pbuffer.data(idx_gr_dg + 197);

    auto grr_z_xy_xxyy = pbuffer.data(idx_gr_dg + 198);

    auto grr_z_xy_xxyz = pbuffer.data(idx_gr_dg + 199);

    auto grr_z_xy_xxzz = pbuffer.data(idx_gr_dg + 200);

    auto grr_z_xy_xyyy = pbuffer.data(idx_gr_dg + 201);

    auto grr_z_xy_xyyz = pbuffer.data(idx_gr_dg + 202);

    auto grr_z_xy_xyzz = pbuffer.data(idx_gr_dg + 203);

    auto grr_z_xy_xzzz = pbuffer.data(idx_gr_dg + 204);

    auto grr_z_xy_yyyy = pbuffer.data(idx_gr_dg + 205);

    auto grr_z_xy_yyyz = pbuffer.data(idx_gr_dg + 206);

    auto grr_z_xy_yyzz = pbuffer.data(idx_gr_dg + 207);

    auto grr_z_xy_yzzz = pbuffer.data(idx_gr_dg + 208);

    auto grr_z_xy_zzzz = pbuffer.data(idx_gr_dg + 209);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxx, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxy, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxz, gr_xy_xxzz, gr_xy_xyy, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyz, gr_xy_xyzz, gr_xy_xzz, gr_xy_xzzz, gr_xy_yyy, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyz, gr_xy_yyzz, gr_xy_yzz, gr_xy_yzzz, gr_xy_zzz, gr_xy_zzzz, grr_z_xy_xxxx, grr_z_xy_xxxy, grr_z_xy_xxxz, grr_z_xy_xxyy, grr_z_xy_xxyz, grr_z_xy_xxzz, grr_z_xy_xyyy, grr_z_xy_xyyz, grr_z_xy_xyzz, grr_z_xy_xzzz, grr_z_xy_yyyy, grr_z_xy_yyyz, grr_z_xy_yyzz, grr_z_xy_yzzz, grr_z_xy_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xy_xxxx[i] = 2.0 * ts_xy_xxxx[i] * gfe_0 * gc_z[i] + gr_xy_xxxx[i] * gc_z[i];

        grr_z_xy_xxxy[i] = 2.0 * ts_xy_xxxy[i] * gfe_0 * gc_z[i] + gr_xy_xxxy[i] * gc_z[i];

        grr_z_xy_xxxz[i] = 2.0 * ts_xy_xxx[i] * gfe2_0 + gr_xy_xxx[i] * gfe_0 + 2.0 * ts_xy_xxxz[i] * gfe_0 * gc_z[i] + gr_xy_xxxz[i] * gc_z[i];

        grr_z_xy_xxyy[i] = 2.0 * ts_xy_xxyy[i] * gfe_0 * gc_z[i] + gr_xy_xxyy[i] * gc_z[i];

        grr_z_xy_xxyz[i] = 2.0 * ts_xy_xxy[i] * gfe2_0 + gr_xy_xxy[i] * gfe_0 + 2.0 * ts_xy_xxyz[i] * gfe_0 * gc_z[i] + gr_xy_xxyz[i] * gc_z[i];

        grr_z_xy_xxzz[i] = 4.0 * ts_xy_xxz[i] * gfe2_0 + 2.0 * gr_xy_xxz[i] * gfe_0 + 2.0 * ts_xy_xxzz[i] * gfe_0 * gc_z[i] + gr_xy_xxzz[i] * gc_z[i];

        grr_z_xy_xyyy[i] = 2.0 * ts_xy_xyyy[i] * gfe_0 * gc_z[i] + gr_xy_xyyy[i] * gc_z[i];

        grr_z_xy_xyyz[i] = 2.0 * ts_xy_xyy[i] * gfe2_0 + gr_xy_xyy[i] * gfe_0 + 2.0 * ts_xy_xyyz[i] * gfe_0 * gc_z[i] + gr_xy_xyyz[i] * gc_z[i];

        grr_z_xy_xyzz[i] = 4.0 * ts_xy_xyz[i] * gfe2_0 + 2.0 * gr_xy_xyz[i] * gfe_0 + 2.0 * ts_xy_xyzz[i] * gfe_0 * gc_z[i] + gr_xy_xyzz[i] * gc_z[i];

        grr_z_xy_xzzz[i] = 6.0 * ts_xy_xzz[i] * gfe2_0 + 3.0 * gr_xy_xzz[i] * gfe_0 + 2.0 * ts_xy_xzzz[i] * gfe_0 * gc_z[i] + gr_xy_xzzz[i] * gc_z[i];

        grr_z_xy_yyyy[i] = 2.0 * ts_xy_yyyy[i] * gfe_0 * gc_z[i] + gr_xy_yyyy[i] * gc_z[i];

        grr_z_xy_yyyz[i] = 2.0 * ts_xy_yyy[i] * gfe2_0 + gr_xy_yyy[i] * gfe_0 + 2.0 * ts_xy_yyyz[i] * gfe_0 * gc_z[i] + gr_xy_yyyz[i] * gc_z[i];

        grr_z_xy_yyzz[i] = 4.0 * ts_xy_yyz[i] * gfe2_0 + 2.0 * gr_xy_yyz[i] * gfe_0 + 2.0 * ts_xy_yyzz[i] * gfe_0 * gc_z[i] + gr_xy_yyzz[i] * gc_z[i];

        grr_z_xy_yzzz[i] = 6.0 * ts_xy_yzz[i] * gfe2_0 + 3.0 * gr_xy_yzz[i] * gfe_0 + 2.0 * ts_xy_yzzz[i] * gfe_0 * gc_z[i] + gr_xy_yzzz[i] * gc_z[i];

        grr_z_xy_zzzz[i] = 8.0 * ts_xy_zzz[i] * gfe2_0 + 4.0 * gr_xy_zzz[i] * gfe_0 + 2.0 * ts_xy_zzzz[i] * gfe_0 * gc_z[i] + gr_xy_zzzz[i] * gc_z[i];
    }

    // Set up 210-225 components of targeted buffer : DG

    auto grr_z_xz_xxxx = pbuffer.data(idx_gr_dg + 210);

    auto grr_z_xz_xxxy = pbuffer.data(idx_gr_dg + 211);

    auto grr_z_xz_xxxz = pbuffer.data(idx_gr_dg + 212);

    auto grr_z_xz_xxyy = pbuffer.data(idx_gr_dg + 213);

    auto grr_z_xz_xxyz = pbuffer.data(idx_gr_dg + 214);

    auto grr_z_xz_xxzz = pbuffer.data(idx_gr_dg + 215);

    auto grr_z_xz_xyyy = pbuffer.data(idx_gr_dg + 216);

    auto grr_z_xz_xyyz = pbuffer.data(idx_gr_dg + 217);

    auto grr_z_xz_xyzz = pbuffer.data(idx_gr_dg + 218);

    auto grr_z_xz_xzzz = pbuffer.data(idx_gr_dg + 219);

    auto grr_z_xz_yyyy = pbuffer.data(idx_gr_dg + 220);

    auto grr_z_xz_yyyz = pbuffer.data(idx_gr_dg + 221);

    auto grr_z_xz_yyzz = pbuffer.data(idx_gr_dg + 222);

    auto grr_z_xz_yzzz = pbuffer.data(idx_gr_dg + 223);

    auto grr_z_xz_zzzz = pbuffer.data(idx_gr_dg + 224);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxyy, gr_x_xxyz, gr_x_xxzz, gr_x_xyyy, gr_x_xyyz, gr_x_xyzz, gr_x_xzzz, gr_x_yyyy, gr_x_yyyz, gr_x_yyzz, gr_x_yzzz, gr_x_zzzz, gr_xz_xxx, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxy, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxz, gr_xz_xxzz, gr_xz_xyy, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyz, gr_xz_xyzz, gr_xz_xzz, gr_xz_xzzz, gr_xz_yyy, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyz, gr_xz_yyzz, gr_xz_yzz, gr_xz_yzzz, gr_xz_zzz, gr_xz_zzzz, grr_z_xz_xxxx, grr_z_xz_xxxy, grr_z_xz_xxxz, grr_z_xz_xxyy, grr_z_xz_xxyz, grr_z_xz_xxzz, grr_z_xz_xyyy, grr_z_xz_xyyz, grr_z_xz_xyzz, grr_z_xz_xzzz, grr_z_xz_yyyy, grr_z_xz_yyyz, grr_z_xz_yyzz, grr_z_xz_yzzz, grr_z_xz_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xz_xxxx[i] = 2.0 * ts_x_xxxx[i] * gfe2_0 + gr_x_xxxx[i] * gfe_0 + 2.0 * ts_xz_xxxx[i] * gfe_0 * gc_z[i] + gr_xz_xxxx[i] * gc_z[i];

        grr_z_xz_xxxy[i] = 2.0 * ts_x_xxxy[i] * gfe2_0 + gr_x_xxxy[i] * gfe_0 + 2.0 * ts_xz_xxxy[i] * gfe_0 * gc_z[i] + gr_xz_xxxy[i] * gc_z[i];

        grr_z_xz_xxxz[i] = 2.0 * ts_x_xxxz[i] * gfe2_0 + gr_x_xxxz[i] * gfe_0 + 2.0 * ts_xz_xxx[i] * gfe2_0 + gr_xz_xxx[i] * gfe_0 + 2.0 * ts_xz_xxxz[i] * gfe_0 * gc_z[i] + gr_xz_xxxz[i] * gc_z[i];

        grr_z_xz_xxyy[i] = 2.0 * ts_x_xxyy[i] * gfe2_0 + gr_x_xxyy[i] * gfe_0 + 2.0 * ts_xz_xxyy[i] * gfe_0 * gc_z[i] + gr_xz_xxyy[i] * gc_z[i];

        grr_z_xz_xxyz[i] = 2.0 * ts_x_xxyz[i] * gfe2_0 + gr_x_xxyz[i] * gfe_0 + 2.0 * ts_xz_xxy[i] * gfe2_0 + gr_xz_xxy[i] * gfe_0 + 2.0 * ts_xz_xxyz[i] * gfe_0 * gc_z[i] + gr_xz_xxyz[i] * gc_z[i];

        grr_z_xz_xxzz[i] = 2.0 * ts_x_xxzz[i] * gfe2_0 + gr_x_xxzz[i] * gfe_0 + 4.0 * ts_xz_xxz[i] * gfe2_0 + 2.0 * gr_xz_xxz[i] * gfe_0 + 2.0 * ts_xz_xxzz[i] * gfe_0 * gc_z[i] + gr_xz_xxzz[i] * gc_z[i];

        grr_z_xz_xyyy[i] = 2.0 * ts_x_xyyy[i] * gfe2_0 + gr_x_xyyy[i] * gfe_0 + 2.0 * ts_xz_xyyy[i] * gfe_0 * gc_z[i] + gr_xz_xyyy[i] * gc_z[i];

        grr_z_xz_xyyz[i] = 2.0 * ts_x_xyyz[i] * gfe2_0 + gr_x_xyyz[i] * gfe_0 + 2.0 * ts_xz_xyy[i] * gfe2_0 + gr_xz_xyy[i] * gfe_0 + 2.0 * ts_xz_xyyz[i] * gfe_0 * gc_z[i] + gr_xz_xyyz[i] * gc_z[i];

        grr_z_xz_xyzz[i] = 2.0 * ts_x_xyzz[i] * gfe2_0 + gr_x_xyzz[i] * gfe_0 + 4.0 * ts_xz_xyz[i] * gfe2_0 + 2.0 * gr_xz_xyz[i] * gfe_0 + 2.0 * ts_xz_xyzz[i] * gfe_0 * gc_z[i] + gr_xz_xyzz[i] * gc_z[i];

        grr_z_xz_xzzz[i] = 2.0 * ts_x_xzzz[i] * gfe2_0 + gr_x_xzzz[i] * gfe_0 + 6.0 * ts_xz_xzz[i] * gfe2_0 + 3.0 * gr_xz_xzz[i] * gfe_0 + 2.0 * ts_xz_xzzz[i] * gfe_0 * gc_z[i] + gr_xz_xzzz[i] * gc_z[i];

        grr_z_xz_yyyy[i] = 2.0 * ts_x_yyyy[i] * gfe2_0 + gr_x_yyyy[i] * gfe_0 + 2.0 * ts_xz_yyyy[i] * gfe_0 * gc_z[i] + gr_xz_yyyy[i] * gc_z[i];

        grr_z_xz_yyyz[i] = 2.0 * ts_x_yyyz[i] * gfe2_0 + gr_x_yyyz[i] * gfe_0 + 2.0 * ts_xz_yyy[i] * gfe2_0 + gr_xz_yyy[i] * gfe_0 + 2.0 * ts_xz_yyyz[i] * gfe_0 * gc_z[i] + gr_xz_yyyz[i] * gc_z[i];

        grr_z_xz_yyzz[i] = 2.0 * ts_x_yyzz[i] * gfe2_0 + gr_x_yyzz[i] * gfe_0 + 4.0 * ts_xz_yyz[i] * gfe2_0 + 2.0 * gr_xz_yyz[i] * gfe_0 + 2.0 * ts_xz_yyzz[i] * gfe_0 * gc_z[i] + gr_xz_yyzz[i] * gc_z[i];

        grr_z_xz_yzzz[i] = 2.0 * ts_x_yzzz[i] * gfe2_0 + gr_x_yzzz[i] * gfe_0 + 6.0 * ts_xz_yzz[i] * gfe2_0 + 3.0 * gr_xz_yzz[i] * gfe_0 + 2.0 * ts_xz_yzzz[i] * gfe_0 * gc_z[i] + gr_xz_yzzz[i] * gc_z[i];

        grr_z_xz_zzzz[i] = 2.0 * ts_x_zzzz[i] * gfe2_0 + gr_x_zzzz[i] * gfe_0 + 8.0 * ts_xz_zzz[i] * gfe2_0 + 4.0 * gr_xz_zzz[i] * gfe_0 + 2.0 * ts_xz_zzzz[i] * gfe_0 * gc_z[i] + gr_xz_zzzz[i] * gc_z[i];
    }

    // Set up 225-240 components of targeted buffer : DG

    auto grr_z_yy_xxxx = pbuffer.data(idx_gr_dg + 225);

    auto grr_z_yy_xxxy = pbuffer.data(idx_gr_dg + 226);

    auto grr_z_yy_xxxz = pbuffer.data(idx_gr_dg + 227);

    auto grr_z_yy_xxyy = pbuffer.data(idx_gr_dg + 228);

    auto grr_z_yy_xxyz = pbuffer.data(idx_gr_dg + 229);

    auto grr_z_yy_xxzz = pbuffer.data(idx_gr_dg + 230);

    auto grr_z_yy_xyyy = pbuffer.data(idx_gr_dg + 231);

    auto grr_z_yy_xyyz = pbuffer.data(idx_gr_dg + 232);

    auto grr_z_yy_xyzz = pbuffer.data(idx_gr_dg + 233);

    auto grr_z_yy_xzzz = pbuffer.data(idx_gr_dg + 234);

    auto grr_z_yy_yyyy = pbuffer.data(idx_gr_dg + 235);

    auto grr_z_yy_yyyz = pbuffer.data(idx_gr_dg + 236);

    auto grr_z_yy_yyzz = pbuffer.data(idx_gr_dg + 237);

    auto grr_z_yy_yzzz = pbuffer.data(idx_gr_dg + 238);

    auto grr_z_yy_zzzz = pbuffer.data(idx_gr_dg + 239);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxx, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxy, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxz, gr_yy_xxzz, gr_yy_xyy, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyz, gr_yy_xyzz, gr_yy_xzz, gr_yy_xzzz, gr_yy_yyy, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyz, gr_yy_yyzz, gr_yy_yzz, gr_yy_yzzz, gr_yy_zzz, gr_yy_zzzz, grr_z_yy_xxxx, grr_z_yy_xxxy, grr_z_yy_xxxz, grr_z_yy_xxyy, grr_z_yy_xxyz, grr_z_yy_xxzz, grr_z_yy_xyyy, grr_z_yy_xyyz, grr_z_yy_xyzz, grr_z_yy_xzzz, grr_z_yy_yyyy, grr_z_yy_yyyz, grr_z_yy_yyzz, grr_z_yy_yzzz, grr_z_yy_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yy_xxxx[i] = 2.0 * ts_yy_xxxx[i] * gfe_0 * gc_z[i] + gr_yy_xxxx[i] * gc_z[i];

        grr_z_yy_xxxy[i] = 2.0 * ts_yy_xxxy[i] * gfe_0 * gc_z[i] + gr_yy_xxxy[i] * gc_z[i];

        grr_z_yy_xxxz[i] = 2.0 * ts_yy_xxx[i] * gfe2_0 + gr_yy_xxx[i] * gfe_0 + 2.0 * ts_yy_xxxz[i] * gfe_0 * gc_z[i] + gr_yy_xxxz[i] * gc_z[i];

        grr_z_yy_xxyy[i] = 2.0 * ts_yy_xxyy[i] * gfe_0 * gc_z[i] + gr_yy_xxyy[i] * gc_z[i];

        grr_z_yy_xxyz[i] = 2.0 * ts_yy_xxy[i] * gfe2_0 + gr_yy_xxy[i] * gfe_0 + 2.0 * ts_yy_xxyz[i] * gfe_0 * gc_z[i] + gr_yy_xxyz[i] * gc_z[i];

        grr_z_yy_xxzz[i] = 4.0 * ts_yy_xxz[i] * gfe2_0 + 2.0 * gr_yy_xxz[i] * gfe_0 + 2.0 * ts_yy_xxzz[i] * gfe_0 * gc_z[i] + gr_yy_xxzz[i] * gc_z[i];

        grr_z_yy_xyyy[i] = 2.0 * ts_yy_xyyy[i] * gfe_0 * gc_z[i] + gr_yy_xyyy[i] * gc_z[i];

        grr_z_yy_xyyz[i] = 2.0 * ts_yy_xyy[i] * gfe2_0 + gr_yy_xyy[i] * gfe_0 + 2.0 * ts_yy_xyyz[i] * gfe_0 * gc_z[i] + gr_yy_xyyz[i] * gc_z[i];

        grr_z_yy_xyzz[i] = 4.0 * ts_yy_xyz[i] * gfe2_0 + 2.0 * gr_yy_xyz[i] * gfe_0 + 2.0 * ts_yy_xyzz[i] * gfe_0 * gc_z[i] + gr_yy_xyzz[i] * gc_z[i];

        grr_z_yy_xzzz[i] = 6.0 * ts_yy_xzz[i] * gfe2_0 + 3.0 * gr_yy_xzz[i] * gfe_0 + 2.0 * ts_yy_xzzz[i] * gfe_0 * gc_z[i] + gr_yy_xzzz[i] * gc_z[i];

        grr_z_yy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * gc_z[i] + gr_yy_yyyy[i] * gc_z[i];

        grr_z_yy_yyyz[i] = 2.0 * ts_yy_yyy[i] * gfe2_0 + gr_yy_yyy[i] * gfe_0 + 2.0 * ts_yy_yyyz[i] * gfe_0 * gc_z[i] + gr_yy_yyyz[i] * gc_z[i];

        grr_z_yy_yyzz[i] = 4.0 * ts_yy_yyz[i] * gfe2_0 + 2.0 * gr_yy_yyz[i] * gfe_0 + 2.0 * ts_yy_yyzz[i] * gfe_0 * gc_z[i] + gr_yy_yyzz[i] * gc_z[i];

        grr_z_yy_yzzz[i] = 6.0 * ts_yy_yzz[i] * gfe2_0 + 3.0 * gr_yy_yzz[i] * gfe_0 + 2.0 * ts_yy_yzzz[i] * gfe_0 * gc_z[i] + gr_yy_yzzz[i] * gc_z[i];

        grr_z_yy_zzzz[i] = 8.0 * ts_yy_zzz[i] * gfe2_0 + 4.0 * gr_yy_zzz[i] * gfe_0 + 2.0 * ts_yy_zzzz[i] * gfe_0 * gc_z[i] + gr_yy_zzzz[i] * gc_z[i];
    }

    // Set up 240-255 components of targeted buffer : DG

    auto grr_z_yz_xxxx = pbuffer.data(idx_gr_dg + 240);

    auto grr_z_yz_xxxy = pbuffer.data(idx_gr_dg + 241);

    auto grr_z_yz_xxxz = pbuffer.data(idx_gr_dg + 242);

    auto grr_z_yz_xxyy = pbuffer.data(idx_gr_dg + 243);

    auto grr_z_yz_xxyz = pbuffer.data(idx_gr_dg + 244);

    auto grr_z_yz_xxzz = pbuffer.data(idx_gr_dg + 245);

    auto grr_z_yz_xyyy = pbuffer.data(idx_gr_dg + 246);

    auto grr_z_yz_xyyz = pbuffer.data(idx_gr_dg + 247);

    auto grr_z_yz_xyzz = pbuffer.data(idx_gr_dg + 248);

    auto grr_z_yz_xzzz = pbuffer.data(idx_gr_dg + 249);

    auto grr_z_yz_yyyy = pbuffer.data(idx_gr_dg + 250);

    auto grr_z_yz_yyyz = pbuffer.data(idx_gr_dg + 251);

    auto grr_z_yz_yyzz = pbuffer.data(idx_gr_dg + 252);

    auto grr_z_yz_yzzz = pbuffer.data(idx_gr_dg + 253);

    auto grr_z_yz_zzzz = pbuffer.data(idx_gr_dg + 254);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxyy, gr_y_xxyz, gr_y_xxzz, gr_y_xyyy, gr_y_xyyz, gr_y_xyzz, gr_y_xzzz, gr_y_yyyy, gr_y_yyyz, gr_y_yyzz, gr_y_yzzz, gr_y_zzzz, gr_yz_xxx, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxy, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxz, gr_yz_xxzz, gr_yz_xyy, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyz, gr_yz_xyzz, gr_yz_xzz, gr_yz_xzzz, gr_yz_yyy, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyz, gr_yz_yyzz, gr_yz_yzz, gr_yz_yzzz, gr_yz_zzz, gr_yz_zzzz, grr_z_yz_xxxx, grr_z_yz_xxxy, grr_z_yz_xxxz, grr_z_yz_xxyy, grr_z_yz_xxyz, grr_z_yz_xxzz, grr_z_yz_xyyy, grr_z_yz_xyyz, grr_z_yz_xyzz, grr_z_yz_xzzz, grr_z_yz_yyyy, grr_z_yz_yyyz, grr_z_yz_yyzz, grr_z_yz_yzzz, grr_z_yz_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yz_xxxx[i] = 2.0 * ts_y_xxxx[i] * gfe2_0 + gr_y_xxxx[i] * gfe_0 + 2.0 * ts_yz_xxxx[i] * gfe_0 * gc_z[i] + gr_yz_xxxx[i] * gc_z[i];

        grr_z_yz_xxxy[i] = 2.0 * ts_y_xxxy[i] * gfe2_0 + gr_y_xxxy[i] * gfe_0 + 2.0 * ts_yz_xxxy[i] * gfe_0 * gc_z[i] + gr_yz_xxxy[i] * gc_z[i];

        grr_z_yz_xxxz[i] = 2.0 * ts_y_xxxz[i] * gfe2_0 + gr_y_xxxz[i] * gfe_0 + 2.0 * ts_yz_xxx[i] * gfe2_0 + gr_yz_xxx[i] * gfe_0 + 2.0 * ts_yz_xxxz[i] * gfe_0 * gc_z[i] + gr_yz_xxxz[i] * gc_z[i];

        grr_z_yz_xxyy[i] = 2.0 * ts_y_xxyy[i] * gfe2_0 + gr_y_xxyy[i] * gfe_0 + 2.0 * ts_yz_xxyy[i] * gfe_0 * gc_z[i] + gr_yz_xxyy[i] * gc_z[i];

        grr_z_yz_xxyz[i] = 2.0 * ts_y_xxyz[i] * gfe2_0 + gr_y_xxyz[i] * gfe_0 + 2.0 * ts_yz_xxy[i] * gfe2_0 + gr_yz_xxy[i] * gfe_0 + 2.0 * ts_yz_xxyz[i] * gfe_0 * gc_z[i] + gr_yz_xxyz[i] * gc_z[i];

        grr_z_yz_xxzz[i] = 2.0 * ts_y_xxzz[i] * gfe2_0 + gr_y_xxzz[i] * gfe_0 + 4.0 * ts_yz_xxz[i] * gfe2_0 + 2.0 * gr_yz_xxz[i] * gfe_0 + 2.0 * ts_yz_xxzz[i] * gfe_0 * gc_z[i] + gr_yz_xxzz[i] * gc_z[i];

        grr_z_yz_xyyy[i] = 2.0 * ts_y_xyyy[i] * gfe2_0 + gr_y_xyyy[i] * gfe_0 + 2.0 * ts_yz_xyyy[i] * gfe_0 * gc_z[i] + gr_yz_xyyy[i] * gc_z[i];

        grr_z_yz_xyyz[i] = 2.0 * ts_y_xyyz[i] * gfe2_0 + gr_y_xyyz[i] * gfe_0 + 2.0 * ts_yz_xyy[i] * gfe2_0 + gr_yz_xyy[i] * gfe_0 + 2.0 * ts_yz_xyyz[i] * gfe_0 * gc_z[i] + gr_yz_xyyz[i] * gc_z[i];

        grr_z_yz_xyzz[i] = 2.0 * ts_y_xyzz[i] * gfe2_0 + gr_y_xyzz[i] * gfe_0 + 4.0 * ts_yz_xyz[i] * gfe2_0 + 2.0 * gr_yz_xyz[i] * gfe_0 + 2.0 * ts_yz_xyzz[i] * gfe_0 * gc_z[i] + gr_yz_xyzz[i] * gc_z[i];

        grr_z_yz_xzzz[i] = 2.0 * ts_y_xzzz[i] * gfe2_0 + gr_y_xzzz[i] * gfe_0 + 6.0 * ts_yz_xzz[i] * gfe2_0 + 3.0 * gr_yz_xzz[i] * gfe_0 + 2.0 * ts_yz_xzzz[i] * gfe_0 * gc_z[i] + gr_yz_xzzz[i] * gc_z[i];

        grr_z_yz_yyyy[i] = 2.0 * ts_y_yyyy[i] * gfe2_0 + gr_y_yyyy[i] * gfe_0 + 2.0 * ts_yz_yyyy[i] * gfe_0 * gc_z[i] + gr_yz_yyyy[i] * gc_z[i];

        grr_z_yz_yyyz[i] = 2.0 * ts_y_yyyz[i] * gfe2_0 + gr_y_yyyz[i] * gfe_0 + 2.0 * ts_yz_yyy[i] * gfe2_0 + gr_yz_yyy[i] * gfe_0 + 2.0 * ts_yz_yyyz[i] * gfe_0 * gc_z[i] + gr_yz_yyyz[i] * gc_z[i];

        grr_z_yz_yyzz[i] = 2.0 * ts_y_yyzz[i] * gfe2_0 + gr_y_yyzz[i] * gfe_0 + 4.0 * ts_yz_yyz[i] * gfe2_0 + 2.0 * gr_yz_yyz[i] * gfe_0 + 2.0 * ts_yz_yyzz[i] * gfe_0 * gc_z[i] + gr_yz_yyzz[i] * gc_z[i];

        grr_z_yz_yzzz[i] = 2.0 * ts_y_yzzz[i] * gfe2_0 + gr_y_yzzz[i] * gfe_0 + 6.0 * ts_yz_yzz[i] * gfe2_0 + 3.0 * gr_yz_yzz[i] * gfe_0 + 2.0 * ts_yz_yzzz[i] * gfe_0 * gc_z[i] + gr_yz_yzzz[i] * gc_z[i];

        grr_z_yz_zzzz[i] = 2.0 * ts_y_zzzz[i] * gfe2_0 + gr_y_zzzz[i] * gfe_0 + 8.0 * ts_yz_zzz[i] * gfe2_0 + 4.0 * gr_yz_zzz[i] * gfe_0 + 2.0 * ts_yz_zzzz[i] * gfe_0 * gc_z[i] + gr_yz_zzzz[i] * gc_z[i];
    }

    // Set up 255-270 components of targeted buffer : DG

    auto grr_z_zz_xxxx = pbuffer.data(idx_gr_dg + 255);

    auto grr_z_zz_xxxy = pbuffer.data(idx_gr_dg + 256);

    auto grr_z_zz_xxxz = pbuffer.data(idx_gr_dg + 257);

    auto grr_z_zz_xxyy = pbuffer.data(idx_gr_dg + 258);

    auto grr_z_zz_xxyz = pbuffer.data(idx_gr_dg + 259);

    auto grr_z_zz_xxzz = pbuffer.data(idx_gr_dg + 260);

    auto grr_z_zz_xyyy = pbuffer.data(idx_gr_dg + 261);

    auto grr_z_zz_xyyz = pbuffer.data(idx_gr_dg + 262);

    auto grr_z_zz_xyzz = pbuffer.data(idx_gr_dg + 263);

    auto grr_z_zz_xzzz = pbuffer.data(idx_gr_dg + 264);

    auto grr_z_zz_yyyy = pbuffer.data(idx_gr_dg + 265);

    auto grr_z_zz_yyyz = pbuffer.data(idx_gr_dg + 266);

    auto grr_z_zz_yyzz = pbuffer.data(idx_gr_dg + 267);

    auto grr_z_zz_yzzz = pbuffer.data(idx_gr_dg + 268);

    auto grr_z_zz_zzzz = pbuffer.data(idx_gr_dg + 269);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxyy, gr_z_xxyz, gr_z_xxzz, gr_z_xyyy, gr_z_xyyz, gr_z_xyzz, gr_z_xzzz, gr_z_yyyy, gr_z_yyyz, gr_z_yyzz, gr_z_yzzz, gr_z_zzzz, gr_zz_xxx, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxy, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxz, gr_zz_xxzz, gr_zz_xyy, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyz, gr_zz_xyzz, gr_zz_xzz, gr_zz_xzzz, gr_zz_yyy, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyz, gr_zz_yyzz, gr_zz_yzz, gr_zz_yzzz, gr_zz_zzz, gr_zz_zzzz, grr_z_zz_xxxx, grr_z_zz_xxxy, grr_z_zz_xxxz, grr_z_zz_xxyy, grr_z_zz_xxyz, grr_z_zz_xxzz, grr_z_zz_xyyy, grr_z_zz_xyyz, grr_z_zz_xyzz, grr_z_zz_xzzz, grr_z_zz_yyyy, grr_z_zz_yyyz, grr_z_zz_yyzz, grr_z_zz_yzzz, grr_z_zz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zz_xxxx[i] = 4.0 * ts_z_xxxx[i] * gfe2_0 + 2.0 * gr_z_xxxx[i] * gfe_0 + 2.0 * ts_zz_xxxx[i] * gfe_0 * gc_z[i] + gr_zz_xxxx[i] * gc_z[i];

        grr_z_zz_xxxy[i] = 4.0 * ts_z_xxxy[i] * gfe2_0 + 2.0 * gr_z_xxxy[i] * gfe_0 + 2.0 * ts_zz_xxxy[i] * gfe_0 * gc_z[i] + gr_zz_xxxy[i] * gc_z[i];

        grr_z_zz_xxxz[i] = 4.0 * ts_z_xxxz[i] * gfe2_0 + 2.0 * gr_z_xxxz[i] * gfe_0 + 2.0 * ts_zz_xxx[i] * gfe2_0 + gr_zz_xxx[i] * gfe_0 + 2.0 * ts_zz_xxxz[i] * gfe_0 * gc_z[i] + gr_zz_xxxz[i] * gc_z[i];

        grr_z_zz_xxyy[i] = 4.0 * ts_z_xxyy[i] * gfe2_0 + 2.0 * gr_z_xxyy[i] * gfe_0 + 2.0 * ts_zz_xxyy[i] * gfe_0 * gc_z[i] + gr_zz_xxyy[i] * gc_z[i];

        grr_z_zz_xxyz[i] = 4.0 * ts_z_xxyz[i] * gfe2_0 + 2.0 * gr_z_xxyz[i] * gfe_0 + 2.0 * ts_zz_xxy[i] * gfe2_0 + gr_zz_xxy[i] * gfe_0 + 2.0 * ts_zz_xxyz[i] * gfe_0 * gc_z[i] + gr_zz_xxyz[i] * gc_z[i];

        grr_z_zz_xxzz[i] = 4.0 * ts_z_xxzz[i] * gfe2_0 + 2.0 * gr_z_xxzz[i] * gfe_0 + 4.0 * ts_zz_xxz[i] * gfe2_0 + 2.0 * gr_zz_xxz[i] * gfe_0 + 2.0 * ts_zz_xxzz[i] * gfe_0 * gc_z[i] + gr_zz_xxzz[i] * gc_z[i];

        grr_z_zz_xyyy[i] = 4.0 * ts_z_xyyy[i] * gfe2_0 + 2.0 * gr_z_xyyy[i] * gfe_0 + 2.0 * ts_zz_xyyy[i] * gfe_0 * gc_z[i] + gr_zz_xyyy[i] * gc_z[i];

        grr_z_zz_xyyz[i] = 4.0 * ts_z_xyyz[i] * gfe2_0 + 2.0 * gr_z_xyyz[i] * gfe_0 + 2.0 * ts_zz_xyy[i] * gfe2_0 + gr_zz_xyy[i] * gfe_0 + 2.0 * ts_zz_xyyz[i] * gfe_0 * gc_z[i] + gr_zz_xyyz[i] * gc_z[i];

        grr_z_zz_xyzz[i] = 4.0 * ts_z_xyzz[i] * gfe2_0 + 2.0 * gr_z_xyzz[i] * gfe_0 + 4.0 * ts_zz_xyz[i] * gfe2_0 + 2.0 * gr_zz_xyz[i] * gfe_0 + 2.0 * ts_zz_xyzz[i] * gfe_0 * gc_z[i] + gr_zz_xyzz[i] * gc_z[i];

        grr_z_zz_xzzz[i] = 4.0 * ts_z_xzzz[i] * gfe2_0 + 2.0 * gr_z_xzzz[i] * gfe_0 + 6.0 * ts_zz_xzz[i] * gfe2_0 + 3.0 * gr_zz_xzz[i] * gfe_0 + 2.0 * ts_zz_xzzz[i] * gfe_0 * gc_z[i] + gr_zz_xzzz[i] * gc_z[i];

        grr_z_zz_yyyy[i] = 4.0 * ts_z_yyyy[i] * gfe2_0 + 2.0 * gr_z_yyyy[i] * gfe_0 + 2.0 * ts_zz_yyyy[i] * gfe_0 * gc_z[i] + gr_zz_yyyy[i] * gc_z[i];

        grr_z_zz_yyyz[i] = 4.0 * ts_z_yyyz[i] * gfe2_0 + 2.0 * gr_z_yyyz[i] * gfe_0 + 2.0 * ts_zz_yyy[i] * gfe2_0 + gr_zz_yyy[i] * gfe_0 + 2.0 * ts_zz_yyyz[i] * gfe_0 * gc_z[i] + gr_zz_yyyz[i] * gc_z[i];

        grr_z_zz_yyzz[i] = 4.0 * ts_z_yyzz[i] * gfe2_0 + 2.0 * gr_z_yyzz[i] * gfe_0 + 4.0 * ts_zz_yyz[i] * gfe2_0 + 2.0 * gr_zz_yyz[i] * gfe_0 + 2.0 * ts_zz_yyzz[i] * gfe_0 * gc_z[i] + gr_zz_yyzz[i] * gc_z[i];

        grr_z_zz_yzzz[i] = 4.0 * ts_z_yzzz[i] * gfe2_0 + 2.0 * gr_z_yzzz[i] * gfe_0 + 6.0 * ts_zz_yzz[i] * gfe2_0 + 3.0 * gr_zz_yzz[i] * gfe_0 + 2.0 * ts_zz_yzzz[i] * gfe_0 * gc_z[i] + gr_zz_yzzz[i] * gc_z[i];

        grr_z_zz_zzzz[i] = 4.0 * ts_z_zzzz[i] * gfe2_0 + 2.0 * gr_z_zzzz[i] * gfe_0 + 8.0 * ts_zz_zzz[i] * gfe2_0 + 4.0 * gr_zz_zzz[i] * gfe_0 + 2.0 * ts_zz_zzzz[i] * gfe_0 * gc_z[i] + gr_zz_zzzz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

