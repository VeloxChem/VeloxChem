#include "ThreeCenterRR2PrimRecFG.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_fg(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_fg,
                  const size_t idx_dg,
                  const size_t idx_g_dg,
                  const size_t idx_ff,
                  const size_t idx_g_ff,
                  const size_t idx_fg,
                  const size_t idx_g_fg,
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

    // Set up components of auxiliary buffer : FF

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

    // Set up components of auxiliary buffer : FG

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

    // Set up 0-15 components of targeted buffer : FG

    auto grr_x_xxx_xxxx = pbuffer.data(idx_gr_fg);

    auto grr_x_xxx_xxxy = pbuffer.data(idx_gr_fg + 1);

    auto grr_x_xxx_xxxz = pbuffer.data(idx_gr_fg + 2);

    auto grr_x_xxx_xxyy = pbuffer.data(idx_gr_fg + 3);

    auto grr_x_xxx_xxyz = pbuffer.data(idx_gr_fg + 4);

    auto grr_x_xxx_xxzz = pbuffer.data(idx_gr_fg + 5);

    auto grr_x_xxx_xyyy = pbuffer.data(idx_gr_fg + 6);

    auto grr_x_xxx_xyyz = pbuffer.data(idx_gr_fg + 7);

    auto grr_x_xxx_xyzz = pbuffer.data(idx_gr_fg + 8);

    auto grr_x_xxx_xzzz = pbuffer.data(idx_gr_fg + 9);

    auto grr_x_xxx_yyyy = pbuffer.data(idx_gr_fg + 10);

    auto grr_x_xxx_yyyz = pbuffer.data(idx_gr_fg + 11);

    auto grr_x_xxx_yyzz = pbuffer.data(idx_gr_fg + 12);

    auto grr_x_xxx_yzzz = pbuffer.data(idx_gr_fg + 13);

    auto grr_x_xxx_zzzz = pbuffer.data(idx_gr_fg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxzz, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyzz, gr_xx_xzzz, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyzz, gr_xx_yzzz, gr_xx_zzzz, gr_xxx_xxx, gr_xxx_xxxx, gr_xxx_xxxy, gr_xxx_xxxz, gr_xxx_xxy, gr_xxx_xxyy, gr_xxx_xxyz, gr_xxx_xxz, gr_xxx_xxzz, gr_xxx_xyy, gr_xxx_xyyy, gr_xxx_xyyz, gr_xxx_xyz, gr_xxx_xyzz, gr_xxx_xzz, gr_xxx_xzzz, gr_xxx_yyy, gr_xxx_yyyy, gr_xxx_yyyz, gr_xxx_yyz, gr_xxx_yyzz, gr_xxx_yzz, gr_xxx_yzzz, gr_xxx_zzz, gr_xxx_zzzz, grr_x_xxx_xxxx, grr_x_xxx_xxxy, grr_x_xxx_xxxz, grr_x_xxx_xxyy, grr_x_xxx_xxyz, grr_x_xxx_xxzz, grr_x_xxx_xyyy, grr_x_xxx_xyyz, grr_x_xxx_xyzz, grr_x_xxx_xzzz, grr_x_xxx_yyyy, grr_x_xxx_yyyz, grr_x_xxx_yyzz, grr_x_xxx_yzzz, grr_x_xxx_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxx_xxxx[i] = 3.0 * ts_xx_xxxx[i] * gfe2_0 + 3.0 * gr_xx_xxxx[i] * gfe_0 + 4.0 * ts_xxx_xxx[i] * gfe2_0 + 4.0 * gr_xxx_xxx[i] * gfe_0 + ts_xxx_xxxx[i] * gfe_0 * gc_x[i] + gr_xxx_xxxx[i] * gc_x[i];

        grr_x_xxx_xxxy[i] = 3.0 * ts_xx_xxxy[i] * gfe2_0 + 3.0 * gr_xx_xxxy[i] * gfe_0 + 3.0 * ts_xxx_xxy[i] * gfe2_0 + 3.0 * gr_xxx_xxy[i] * gfe_0 + ts_xxx_xxxy[i] * gfe_0 * gc_x[i] + gr_xxx_xxxy[i] * gc_x[i];

        grr_x_xxx_xxxz[i] = 3.0 * ts_xx_xxxz[i] * gfe2_0 + 3.0 * gr_xx_xxxz[i] * gfe_0 + 3.0 * ts_xxx_xxz[i] * gfe2_0 + 3.0 * gr_xxx_xxz[i] * gfe_0 + ts_xxx_xxxz[i] * gfe_0 * gc_x[i] + gr_xxx_xxxz[i] * gc_x[i];

        grr_x_xxx_xxyy[i] = 3.0 * ts_xx_xxyy[i] * gfe2_0 + 3.0 * gr_xx_xxyy[i] * gfe_0 + 2.0 * ts_xxx_xyy[i] * gfe2_0 + 2.0 * gr_xxx_xyy[i] * gfe_0 + ts_xxx_xxyy[i] * gfe_0 * gc_x[i] + gr_xxx_xxyy[i] * gc_x[i];

        grr_x_xxx_xxyz[i] = 3.0 * ts_xx_xxyz[i] * gfe2_0 + 3.0 * gr_xx_xxyz[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe2_0 + 2.0 * gr_xxx_xyz[i] * gfe_0 + ts_xxx_xxyz[i] * gfe_0 * gc_x[i] + gr_xxx_xxyz[i] * gc_x[i];

        grr_x_xxx_xxzz[i] = 3.0 * ts_xx_xxzz[i] * gfe2_0 + 3.0 * gr_xx_xxzz[i] * gfe_0 + 2.0 * ts_xxx_xzz[i] * gfe2_0 + 2.0 * gr_xxx_xzz[i] * gfe_0 + ts_xxx_xxzz[i] * gfe_0 * gc_x[i] + gr_xxx_xxzz[i] * gc_x[i];

        grr_x_xxx_xyyy[i] = 3.0 * ts_xx_xyyy[i] * gfe2_0 + 3.0 * gr_xx_xyyy[i] * gfe_0 + ts_xxx_yyy[i] * gfe2_0 + gr_xxx_yyy[i] * gfe_0 + ts_xxx_xyyy[i] * gfe_0 * gc_x[i] + gr_xxx_xyyy[i] * gc_x[i];

        grr_x_xxx_xyyz[i] = 3.0 * ts_xx_xyyz[i] * gfe2_0 + 3.0 * gr_xx_xyyz[i] * gfe_0 + ts_xxx_yyz[i] * gfe2_0 + gr_xxx_yyz[i] * gfe_0 + ts_xxx_xyyz[i] * gfe_0 * gc_x[i] + gr_xxx_xyyz[i] * gc_x[i];

        grr_x_xxx_xyzz[i] = 3.0 * ts_xx_xyzz[i] * gfe2_0 + 3.0 * gr_xx_xyzz[i] * gfe_0 + ts_xxx_yzz[i] * gfe2_0 + gr_xxx_yzz[i] * gfe_0 + ts_xxx_xyzz[i] * gfe_0 * gc_x[i] + gr_xxx_xyzz[i] * gc_x[i];

        grr_x_xxx_xzzz[i] = 3.0 * ts_xx_xzzz[i] * gfe2_0 + 3.0 * gr_xx_xzzz[i] * gfe_0 + ts_xxx_zzz[i] * gfe2_0 + gr_xxx_zzz[i] * gfe_0 + ts_xxx_xzzz[i] * gfe_0 * gc_x[i] + gr_xxx_xzzz[i] * gc_x[i];

        grr_x_xxx_yyyy[i] = 3.0 * ts_xx_yyyy[i] * gfe2_0 + 3.0 * gr_xx_yyyy[i] * gfe_0 + ts_xxx_yyyy[i] * gfe_0 * gc_x[i] + gr_xxx_yyyy[i] * gc_x[i];

        grr_x_xxx_yyyz[i] = 3.0 * ts_xx_yyyz[i] * gfe2_0 + 3.0 * gr_xx_yyyz[i] * gfe_0 + ts_xxx_yyyz[i] * gfe_0 * gc_x[i] + gr_xxx_yyyz[i] * gc_x[i];

        grr_x_xxx_yyzz[i] = 3.0 * ts_xx_yyzz[i] * gfe2_0 + 3.0 * gr_xx_yyzz[i] * gfe_0 + ts_xxx_yyzz[i] * gfe_0 * gc_x[i] + gr_xxx_yyzz[i] * gc_x[i];

        grr_x_xxx_yzzz[i] = 3.0 * ts_xx_yzzz[i] * gfe2_0 + 3.0 * gr_xx_yzzz[i] * gfe_0 + ts_xxx_yzzz[i] * gfe_0 * gc_x[i] + gr_xxx_yzzz[i] * gc_x[i];

        grr_x_xxx_zzzz[i] = 3.0 * ts_xx_zzzz[i] * gfe2_0 + 3.0 * gr_xx_zzzz[i] * gfe_0 + ts_xxx_zzzz[i] * gfe_0 * gc_x[i] + gr_xxx_zzzz[i] * gc_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto grr_x_xxy_xxxx = pbuffer.data(idx_gr_fg + 15);

    auto grr_x_xxy_xxxy = pbuffer.data(idx_gr_fg + 16);

    auto grr_x_xxy_xxxz = pbuffer.data(idx_gr_fg + 17);

    auto grr_x_xxy_xxyy = pbuffer.data(idx_gr_fg + 18);

    auto grr_x_xxy_xxyz = pbuffer.data(idx_gr_fg + 19);

    auto grr_x_xxy_xxzz = pbuffer.data(idx_gr_fg + 20);

    auto grr_x_xxy_xyyy = pbuffer.data(idx_gr_fg + 21);

    auto grr_x_xxy_xyyz = pbuffer.data(idx_gr_fg + 22);

    auto grr_x_xxy_xyzz = pbuffer.data(idx_gr_fg + 23);

    auto grr_x_xxy_xzzz = pbuffer.data(idx_gr_fg + 24);

    auto grr_x_xxy_yyyy = pbuffer.data(idx_gr_fg + 25);

    auto grr_x_xxy_yyyz = pbuffer.data(idx_gr_fg + 26);

    auto grr_x_xxy_yyzz = pbuffer.data(idx_gr_fg + 27);

    auto grr_x_xxy_yzzz = pbuffer.data(idx_gr_fg + 28);

    auto grr_x_xxy_zzzz = pbuffer.data(idx_gr_fg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xxx, gr_xxy_xxxx, gr_xxy_xxxy, gr_xxy_xxxz, gr_xxy_xxy, gr_xxy_xxyy, gr_xxy_xxyz, gr_xxy_xxz, gr_xxy_xxzz, gr_xxy_xyy, gr_xxy_xyyy, gr_xxy_xyyz, gr_xxy_xyz, gr_xxy_xyzz, gr_xxy_xzz, gr_xxy_xzzz, gr_xxy_yyy, gr_xxy_yyyy, gr_xxy_yyyz, gr_xxy_yyz, gr_xxy_yyzz, gr_xxy_yzz, gr_xxy_yzzz, gr_xxy_zzz, gr_xxy_zzzz, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxzz, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyzz, gr_xy_xzzz, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyzz, gr_xy_yzzz, gr_xy_zzzz, grr_x_xxy_xxxx, grr_x_xxy_xxxy, grr_x_xxy_xxxz, grr_x_xxy_xxyy, grr_x_xxy_xxyz, grr_x_xxy_xxzz, grr_x_xxy_xyyy, grr_x_xxy_xyyz, grr_x_xxy_xyzz, grr_x_xxy_xzzz, grr_x_xxy_yyyy, grr_x_xxy_yyyz, grr_x_xxy_yyzz, grr_x_xxy_yzzz, grr_x_xxy_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxy_xxxx[i] = 2.0 * ts_xy_xxxx[i] * gfe2_0 + 2.0 * gr_xy_xxxx[i] * gfe_0 + 4.0 * ts_xxy_xxx[i] * gfe2_0 + 4.0 * gr_xxy_xxx[i] * gfe_0 + ts_xxy_xxxx[i] * gfe_0 * gc_x[i] + gr_xxy_xxxx[i] * gc_x[i];

        grr_x_xxy_xxxy[i] = 2.0 * ts_xy_xxxy[i] * gfe2_0 + 2.0 * gr_xy_xxxy[i] * gfe_0 + 3.0 * ts_xxy_xxy[i] * gfe2_0 + 3.0 * gr_xxy_xxy[i] * gfe_0 + ts_xxy_xxxy[i] * gfe_0 * gc_x[i] + gr_xxy_xxxy[i] * gc_x[i];

        grr_x_xxy_xxxz[i] = 2.0 * ts_xy_xxxz[i] * gfe2_0 + 2.0 * gr_xy_xxxz[i] * gfe_0 + 3.0 * ts_xxy_xxz[i] * gfe2_0 + 3.0 * gr_xxy_xxz[i] * gfe_0 + ts_xxy_xxxz[i] * gfe_0 * gc_x[i] + gr_xxy_xxxz[i] * gc_x[i];

        grr_x_xxy_xxyy[i] = 2.0 * ts_xy_xxyy[i] * gfe2_0 + 2.0 * gr_xy_xxyy[i] * gfe_0 + 2.0 * ts_xxy_xyy[i] * gfe2_0 + 2.0 * gr_xxy_xyy[i] * gfe_0 + ts_xxy_xxyy[i] * gfe_0 * gc_x[i] + gr_xxy_xxyy[i] * gc_x[i];

        grr_x_xxy_xxyz[i] = 2.0 * ts_xy_xxyz[i] * gfe2_0 + 2.0 * gr_xy_xxyz[i] * gfe_0 + 2.0 * ts_xxy_xyz[i] * gfe2_0 + 2.0 * gr_xxy_xyz[i] * gfe_0 + ts_xxy_xxyz[i] * gfe_0 * gc_x[i] + gr_xxy_xxyz[i] * gc_x[i];

        grr_x_xxy_xxzz[i] = 2.0 * ts_xy_xxzz[i] * gfe2_0 + 2.0 * gr_xy_xxzz[i] * gfe_0 + 2.0 * ts_xxy_xzz[i] * gfe2_0 + 2.0 * gr_xxy_xzz[i] * gfe_0 + ts_xxy_xxzz[i] * gfe_0 * gc_x[i] + gr_xxy_xxzz[i] * gc_x[i];

        grr_x_xxy_xyyy[i] = 2.0 * ts_xy_xyyy[i] * gfe2_0 + 2.0 * gr_xy_xyyy[i] * gfe_0 + ts_xxy_yyy[i] * gfe2_0 + gr_xxy_yyy[i] * gfe_0 + ts_xxy_xyyy[i] * gfe_0 * gc_x[i] + gr_xxy_xyyy[i] * gc_x[i];

        grr_x_xxy_xyyz[i] = 2.0 * ts_xy_xyyz[i] * gfe2_0 + 2.0 * gr_xy_xyyz[i] * gfe_0 + ts_xxy_yyz[i] * gfe2_0 + gr_xxy_yyz[i] * gfe_0 + ts_xxy_xyyz[i] * gfe_0 * gc_x[i] + gr_xxy_xyyz[i] * gc_x[i];

        grr_x_xxy_xyzz[i] = 2.0 * ts_xy_xyzz[i] * gfe2_0 + 2.0 * gr_xy_xyzz[i] * gfe_0 + ts_xxy_yzz[i] * gfe2_0 + gr_xxy_yzz[i] * gfe_0 + ts_xxy_xyzz[i] * gfe_0 * gc_x[i] + gr_xxy_xyzz[i] * gc_x[i];

        grr_x_xxy_xzzz[i] = 2.0 * ts_xy_xzzz[i] * gfe2_0 + 2.0 * gr_xy_xzzz[i] * gfe_0 + ts_xxy_zzz[i] * gfe2_0 + gr_xxy_zzz[i] * gfe_0 + ts_xxy_xzzz[i] * gfe_0 * gc_x[i] + gr_xxy_xzzz[i] * gc_x[i];

        grr_x_xxy_yyyy[i] = 2.0 * ts_xy_yyyy[i] * gfe2_0 + 2.0 * gr_xy_yyyy[i] * gfe_0 + ts_xxy_yyyy[i] * gfe_0 * gc_x[i] + gr_xxy_yyyy[i] * gc_x[i];

        grr_x_xxy_yyyz[i] = 2.0 * ts_xy_yyyz[i] * gfe2_0 + 2.0 * gr_xy_yyyz[i] * gfe_0 + ts_xxy_yyyz[i] * gfe_0 * gc_x[i] + gr_xxy_yyyz[i] * gc_x[i];

        grr_x_xxy_yyzz[i] = 2.0 * ts_xy_yyzz[i] * gfe2_0 + 2.0 * gr_xy_yyzz[i] * gfe_0 + ts_xxy_yyzz[i] * gfe_0 * gc_x[i] + gr_xxy_yyzz[i] * gc_x[i];

        grr_x_xxy_yzzz[i] = 2.0 * ts_xy_yzzz[i] * gfe2_0 + 2.0 * gr_xy_yzzz[i] * gfe_0 + ts_xxy_yzzz[i] * gfe_0 * gc_x[i] + gr_xxy_yzzz[i] * gc_x[i];

        grr_x_xxy_zzzz[i] = 2.0 * ts_xy_zzzz[i] * gfe2_0 + 2.0 * gr_xy_zzzz[i] * gfe_0 + ts_xxy_zzzz[i] * gfe_0 * gc_x[i] + gr_xxy_zzzz[i] * gc_x[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto grr_x_xxz_xxxx = pbuffer.data(idx_gr_fg + 30);

    auto grr_x_xxz_xxxy = pbuffer.data(idx_gr_fg + 31);

    auto grr_x_xxz_xxxz = pbuffer.data(idx_gr_fg + 32);

    auto grr_x_xxz_xxyy = pbuffer.data(idx_gr_fg + 33);

    auto grr_x_xxz_xxyz = pbuffer.data(idx_gr_fg + 34);

    auto grr_x_xxz_xxzz = pbuffer.data(idx_gr_fg + 35);

    auto grr_x_xxz_xyyy = pbuffer.data(idx_gr_fg + 36);

    auto grr_x_xxz_xyyz = pbuffer.data(idx_gr_fg + 37);

    auto grr_x_xxz_xyzz = pbuffer.data(idx_gr_fg + 38);

    auto grr_x_xxz_xzzz = pbuffer.data(idx_gr_fg + 39);

    auto grr_x_xxz_yyyy = pbuffer.data(idx_gr_fg + 40);

    auto grr_x_xxz_yyyz = pbuffer.data(idx_gr_fg + 41);

    auto grr_x_xxz_yyzz = pbuffer.data(idx_gr_fg + 42);

    auto grr_x_xxz_yzzz = pbuffer.data(idx_gr_fg + 43);

    auto grr_x_xxz_zzzz = pbuffer.data(idx_gr_fg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xxx, gr_xxz_xxxx, gr_xxz_xxxy, gr_xxz_xxxz, gr_xxz_xxy, gr_xxz_xxyy, gr_xxz_xxyz, gr_xxz_xxz, gr_xxz_xxzz, gr_xxz_xyy, gr_xxz_xyyy, gr_xxz_xyyz, gr_xxz_xyz, gr_xxz_xyzz, gr_xxz_xzz, gr_xxz_xzzz, gr_xxz_yyy, gr_xxz_yyyy, gr_xxz_yyyz, gr_xxz_yyz, gr_xxz_yyzz, gr_xxz_yzz, gr_xxz_yzzz, gr_xxz_zzz, gr_xxz_zzzz, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxzz, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyzz, gr_xz_xzzz, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyzz, gr_xz_yzzz, gr_xz_zzzz, grr_x_xxz_xxxx, grr_x_xxz_xxxy, grr_x_xxz_xxxz, grr_x_xxz_xxyy, grr_x_xxz_xxyz, grr_x_xxz_xxzz, grr_x_xxz_xyyy, grr_x_xxz_xyyz, grr_x_xxz_xyzz, grr_x_xxz_xzzz, grr_x_xxz_yyyy, grr_x_xxz_yyyz, grr_x_xxz_yyzz, grr_x_xxz_yzzz, grr_x_xxz_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxz_xxxx[i] = 2.0 * ts_xz_xxxx[i] * gfe2_0 + 2.0 * gr_xz_xxxx[i] * gfe_0 + 4.0 * ts_xxz_xxx[i] * gfe2_0 + 4.0 * gr_xxz_xxx[i] * gfe_0 + ts_xxz_xxxx[i] * gfe_0 * gc_x[i] + gr_xxz_xxxx[i] * gc_x[i];

        grr_x_xxz_xxxy[i] = 2.0 * ts_xz_xxxy[i] * gfe2_0 + 2.0 * gr_xz_xxxy[i] * gfe_0 + 3.0 * ts_xxz_xxy[i] * gfe2_0 + 3.0 * gr_xxz_xxy[i] * gfe_0 + ts_xxz_xxxy[i] * gfe_0 * gc_x[i] + gr_xxz_xxxy[i] * gc_x[i];

        grr_x_xxz_xxxz[i] = 2.0 * ts_xz_xxxz[i] * gfe2_0 + 2.0 * gr_xz_xxxz[i] * gfe_0 + 3.0 * ts_xxz_xxz[i] * gfe2_0 + 3.0 * gr_xxz_xxz[i] * gfe_0 + ts_xxz_xxxz[i] * gfe_0 * gc_x[i] + gr_xxz_xxxz[i] * gc_x[i];

        grr_x_xxz_xxyy[i] = 2.0 * ts_xz_xxyy[i] * gfe2_0 + 2.0 * gr_xz_xxyy[i] * gfe_0 + 2.0 * ts_xxz_xyy[i] * gfe2_0 + 2.0 * gr_xxz_xyy[i] * gfe_0 + ts_xxz_xxyy[i] * gfe_0 * gc_x[i] + gr_xxz_xxyy[i] * gc_x[i];

        grr_x_xxz_xxyz[i] = 2.0 * ts_xz_xxyz[i] * gfe2_0 + 2.0 * gr_xz_xxyz[i] * gfe_0 + 2.0 * ts_xxz_xyz[i] * gfe2_0 + 2.0 * gr_xxz_xyz[i] * gfe_0 + ts_xxz_xxyz[i] * gfe_0 * gc_x[i] + gr_xxz_xxyz[i] * gc_x[i];

        grr_x_xxz_xxzz[i] = 2.0 * ts_xz_xxzz[i] * gfe2_0 + 2.0 * gr_xz_xxzz[i] * gfe_0 + 2.0 * ts_xxz_xzz[i] * gfe2_0 + 2.0 * gr_xxz_xzz[i] * gfe_0 + ts_xxz_xxzz[i] * gfe_0 * gc_x[i] + gr_xxz_xxzz[i] * gc_x[i];

        grr_x_xxz_xyyy[i] = 2.0 * ts_xz_xyyy[i] * gfe2_0 + 2.0 * gr_xz_xyyy[i] * gfe_0 + ts_xxz_yyy[i] * gfe2_0 + gr_xxz_yyy[i] * gfe_0 + ts_xxz_xyyy[i] * gfe_0 * gc_x[i] + gr_xxz_xyyy[i] * gc_x[i];

        grr_x_xxz_xyyz[i] = 2.0 * ts_xz_xyyz[i] * gfe2_0 + 2.0 * gr_xz_xyyz[i] * gfe_0 + ts_xxz_yyz[i] * gfe2_0 + gr_xxz_yyz[i] * gfe_0 + ts_xxz_xyyz[i] * gfe_0 * gc_x[i] + gr_xxz_xyyz[i] * gc_x[i];

        grr_x_xxz_xyzz[i] = 2.0 * ts_xz_xyzz[i] * gfe2_0 + 2.0 * gr_xz_xyzz[i] * gfe_0 + ts_xxz_yzz[i] * gfe2_0 + gr_xxz_yzz[i] * gfe_0 + ts_xxz_xyzz[i] * gfe_0 * gc_x[i] + gr_xxz_xyzz[i] * gc_x[i];

        grr_x_xxz_xzzz[i] = 2.0 * ts_xz_xzzz[i] * gfe2_0 + 2.0 * gr_xz_xzzz[i] * gfe_0 + ts_xxz_zzz[i] * gfe2_0 + gr_xxz_zzz[i] * gfe_0 + ts_xxz_xzzz[i] * gfe_0 * gc_x[i] + gr_xxz_xzzz[i] * gc_x[i];

        grr_x_xxz_yyyy[i] = 2.0 * ts_xz_yyyy[i] * gfe2_0 + 2.0 * gr_xz_yyyy[i] * gfe_0 + ts_xxz_yyyy[i] * gfe_0 * gc_x[i] + gr_xxz_yyyy[i] * gc_x[i];

        grr_x_xxz_yyyz[i] = 2.0 * ts_xz_yyyz[i] * gfe2_0 + 2.0 * gr_xz_yyyz[i] * gfe_0 + ts_xxz_yyyz[i] * gfe_0 * gc_x[i] + gr_xxz_yyyz[i] * gc_x[i];

        grr_x_xxz_yyzz[i] = 2.0 * ts_xz_yyzz[i] * gfe2_0 + 2.0 * gr_xz_yyzz[i] * gfe_0 + ts_xxz_yyzz[i] * gfe_0 * gc_x[i] + gr_xxz_yyzz[i] * gc_x[i];

        grr_x_xxz_yzzz[i] = 2.0 * ts_xz_yzzz[i] * gfe2_0 + 2.0 * gr_xz_yzzz[i] * gfe_0 + ts_xxz_yzzz[i] * gfe_0 * gc_x[i] + gr_xxz_yzzz[i] * gc_x[i];

        grr_x_xxz_zzzz[i] = 2.0 * ts_xz_zzzz[i] * gfe2_0 + 2.0 * gr_xz_zzzz[i] * gfe_0 + ts_xxz_zzzz[i] * gfe_0 * gc_x[i] + gr_xxz_zzzz[i] * gc_x[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto grr_x_xyy_xxxx = pbuffer.data(idx_gr_fg + 45);

    auto grr_x_xyy_xxxy = pbuffer.data(idx_gr_fg + 46);

    auto grr_x_xyy_xxxz = pbuffer.data(idx_gr_fg + 47);

    auto grr_x_xyy_xxyy = pbuffer.data(idx_gr_fg + 48);

    auto grr_x_xyy_xxyz = pbuffer.data(idx_gr_fg + 49);

    auto grr_x_xyy_xxzz = pbuffer.data(idx_gr_fg + 50);

    auto grr_x_xyy_xyyy = pbuffer.data(idx_gr_fg + 51);

    auto grr_x_xyy_xyyz = pbuffer.data(idx_gr_fg + 52);

    auto grr_x_xyy_xyzz = pbuffer.data(idx_gr_fg + 53);

    auto grr_x_xyy_xzzz = pbuffer.data(idx_gr_fg + 54);

    auto grr_x_xyy_yyyy = pbuffer.data(idx_gr_fg + 55);

    auto grr_x_xyy_yyyz = pbuffer.data(idx_gr_fg + 56);

    auto grr_x_xyy_yyzz = pbuffer.data(idx_gr_fg + 57);

    auto grr_x_xyy_yzzz = pbuffer.data(idx_gr_fg + 58);

    auto grr_x_xyy_zzzz = pbuffer.data(idx_gr_fg + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xxx, gr_xyy_xxxx, gr_xyy_xxxy, gr_xyy_xxxz, gr_xyy_xxy, gr_xyy_xxyy, gr_xyy_xxyz, gr_xyy_xxz, gr_xyy_xxzz, gr_xyy_xyy, gr_xyy_xyyy, gr_xyy_xyyz, gr_xyy_xyz, gr_xyy_xyzz, gr_xyy_xzz, gr_xyy_xzzz, gr_xyy_yyy, gr_xyy_yyyy, gr_xyy_yyyz, gr_xyy_yyz, gr_xyy_yyzz, gr_xyy_yzz, gr_xyy_yzzz, gr_xyy_zzz, gr_xyy_zzzz, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxzz, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyzz, gr_yy_xzzz, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyzz, gr_yy_yzzz, gr_yy_zzzz, grr_x_xyy_xxxx, grr_x_xyy_xxxy, grr_x_xyy_xxxz, grr_x_xyy_xxyy, grr_x_xyy_xxyz, grr_x_xyy_xxzz, grr_x_xyy_xyyy, grr_x_xyy_xyyz, grr_x_xyy_xyzz, grr_x_xyy_xzzz, grr_x_xyy_yyyy, grr_x_xyy_yyyz, grr_x_xyy_yyzz, grr_x_xyy_yzzz, grr_x_xyy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyy_xxxx[i] = ts_yy_xxxx[i] * gfe2_0 + gr_yy_xxxx[i] * gfe_0 + 4.0 * ts_xyy_xxx[i] * gfe2_0 + 4.0 * gr_xyy_xxx[i] * gfe_0 + ts_xyy_xxxx[i] * gfe_0 * gc_x[i] + gr_xyy_xxxx[i] * gc_x[i];

        grr_x_xyy_xxxy[i] = ts_yy_xxxy[i] * gfe2_0 + gr_yy_xxxy[i] * gfe_0 + 3.0 * ts_xyy_xxy[i] * gfe2_0 + 3.0 * gr_xyy_xxy[i] * gfe_0 + ts_xyy_xxxy[i] * gfe_0 * gc_x[i] + gr_xyy_xxxy[i] * gc_x[i];

        grr_x_xyy_xxxz[i] = ts_yy_xxxz[i] * gfe2_0 + gr_yy_xxxz[i] * gfe_0 + 3.0 * ts_xyy_xxz[i] * gfe2_0 + 3.0 * gr_xyy_xxz[i] * gfe_0 + ts_xyy_xxxz[i] * gfe_0 * gc_x[i] + gr_xyy_xxxz[i] * gc_x[i];

        grr_x_xyy_xxyy[i] = ts_yy_xxyy[i] * gfe2_0 + gr_yy_xxyy[i] * gfe_0 + 2.0 * ts_xyy_xyy[i] * gfe2_0 + 2.0 * gr_xyy_xyy[i] * gfe_0 + ts_xyy_xxyy[i] * gfe_0 * gc_x[i] + gr_xyy_xxyy[i] * gc_x[i];

        grr_x_xyy_xxyz[i] = ts_yy_xxyz[i] * gfe2_0 + gr_yy_xxyz[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe2_0 + 2.0 * gr_xyy_xyz[i] * gfe_0 + ts_xyy_xxyz[i] * gfe_0 * gc_x[i] + gr_xyy_xxyz[i] * gc_x[i];

        grr_x_xyy_xxzz[i] = ts_yy_xxzz[i] * gfe2_0 + gr_yy_xxzz[i] * gfe_0 + 2.0 * ts_xyy_xzz[i] * gfe2_0 + 2.0 * gr_xyy_xzz[i] * gfe_0 + ts_xyy_xxzz[i] * gfe_0 * gc_x[i] + gr_xyy_xxzz[i] * gc_x[i];

        grr_x_xyy_xyyy[i] = ts_yy_xyyy[i] * gfe2_0 + gr_yy_xyyy[i] * gfe_0 + ts_xyy_yyy[i] * gfe2_0 + gr_xyy_yyy[i] * gfe_0 + ts_xyy_xyyy[i] * gfe_0 * gc_x[i] + gr_xyy_xyyy[i] * gc_x[i];

        grr_x_xyy_xyyz[i] = ts_yy_xyyz[i] * gfe2_0 + gr_yy_xyyz[i] * gfe_0 + ts_xyy_yyz[i] * gfe2_0 + gr_xyy_yyz[i] * gfe_0 + ts_xyy_xyyz[i] * gfe_0 * gc_x[i] + gr_xyy_xyyz[i] * gc_x[i];

        grr_x_xyy_xyzz[i] = ts_yy_xyzz[i] * gfe2_0 + gr_yy_xyzz[i] * gfe_0 + ts_xyy_yzz[i] * gfe2_0 + gr_xyy_yzz[i] * gfe_0 + ts_xyy_xyzz[i] * gfe_0 * gc_x[i] + gr_xyy_xyzz[i] * gc_x[i];

        grr_x_xyy_xzzz[i] = ts_yy_xzzz[i] * gfe2_0 + gr_yy_xzzz[i] * gfe_0 + ts_xyy_zzz[i] * gfe2_0 + gr_xyy_zzz[i] * gfe_0 + ts_xyy_xzzz[i] * gfe_0 * gc_x[i] + gr_xyy_xzzz[i] * gc_x[i];

        grr_x_xyy_yyyy[i] = ts_yy_yyyy[i] * gfe2_0 + gr_yy_yyyy[i] * gfe_0 + ts_xyy_yyyy[i] * gfe_0 * gc_x[i] + gr_xyy_yyyy[i] * gc_x[i];

        grr_x_xyy_yyyz[i] = ts_yy_yyyz[i] * gfe2_0 + gr_yy_yyyz[i] * gfe_0 + ts_xyy_yyyz[i] * gfe_0 * gc_x[i] + gr_xyy_yyyz[i] * gc_x[i];

        grr_x_xyy_yyzz[i] = ts_yy_yyzz[i] * gfe2_0 + gr_yy_yyzz[i] * gfe_0 + ts_xyy_yyzz[i] * gfe_0 * gc_x[i] + gr_xyy_yyzz[i] * gc_x[i];

        grr_x_xyy_yzzz[i] = ts_yy_yzzz[i] * gfe2_0 + gr_yy_yzzz[i] * gfe_0 + ts_xyy_yzzz[i] * gfe_0 * gc_x[i] + gr_xyy_yzzz[i] * gc_x[i];

        grr_x_xyy_zzzz[i] = ts_yy_zzzz[i] * gfe2_0 + gr_yy_zzzz[i] * gfe_0 + ts_xyy_zzzz[i] * gfe_0 * gc_x[i] + gr_xyy_zzzz[i] * gc_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto grr_x_xyz_xxxx = pbuffer.data(idx_gr_fg + 60);

    auto grr_x_xyz_xxxy = pbuffer.data(idx_gr_fg + 61);

    auto grr_x_xyz_xxxz = pbuffer.data(idx_gr_fg + 62);

    auto grr_x_xyz_xxyy = pbuffer.data(idx_gr_fg + 63);

    auto grr_x_xyz_xxyz = pbuffer.data(idx_gr_fg + 64);

    auto grr_x_xyz_xxzz = pbuffer.data(idx_gr_fg + 65);

    auto grr_x_xyz_xyyy = pbuffer.data(idx_gr_fg + 66);

    auto grr_x_xyz_xyyz = pbuffer.data(idx_gr_fg + 67);

    auto grr_x_xyz_xyzz = pbuffer.data(idx_gr_fg + 68);

    auto grr_x_xyz_xzzz = pbuffer.data(idx_gr_fg + 69);

    auto grr_x_xyz_yyyy = pbuffer.data(idx_gr_fg + 70);

    auto grr_x_xyz_yyyz = pbuffer.data(idx_gr_fg + 71);

    auto grr_x_xyz_yyzz = pbuffer.data(idx_gr_fg + 72);

    auto grr_x_xyz_yzzz = pbuffer.data(idx_gr_fg + 73);

    auto grr_x_xyz_zzzz = pbuffer.data(idx_gr_fg + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xxx, gr_xyz_xxxx, gr_xyz_xxxy, gr_xyz_xxxz, gr_xyz_xxy, gr_xyz_xxyy, gr_xyz_xxyz, gr_xyz_xxz, gr_xyz_xxzz, gr_xyz_xyy, gr_xyz_xyyy, gr_xyz_xyyz, gr_xyz_xyz, gr_xyz_xyzz, gr_xyz_xzz, gr_xyz_xzzz, gr_xyz_yyy, gr_xyz_yyyy, gr_xyz_yyyz, gr_xyz_yyz, gr_xyz_yyzz, gr_xyz_yzz, gr_xyz_yzzz, gr_xyz_zzz, gr_xyz_zzzz, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxzz, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyzz, gr_yz_xzzz, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyzz, gr_yz_yzzz, gr_yz_zzzz, grr_x_xyz_xxxx, grr_x_xyz_xxxy, grr_x_xyz_xxxz, grr_x_xyz_xxyy, grr_x_xyz_xxyz, grr_x_xyz_xxzz, grr_x_xyz_xyyy, grr_x_xyz_xyyz, grr_x_xyz_xyzz, grr_x_xyz_xzzz, grr_x_xyz_yyyy, grr_x_xyz_yyyz, grr_x_xyz_yyzz, grr_x_xyz_yzzz, grr_x_xyz_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyz_xxxx[i] = ts_yz_xxxx[i] * gfe2_0 + gr_yz_xxxx[i] * gfe_0 + 4.0 * ts_xyz_xxx[i] * gfe2_0 + 4.0 * gr_xyz_xxx[i] * gfe_0 + ts_xyz_xxxx[i] * gfe_0 * gc_x[i] + gr_xyz_xxxx[i] * gc_x[i];

        grr_x_xyz_xxxy[i] = ts_yz_xxxy[i] * gfe2_0 + gr_yz_xxxy[i] * gfe_0 + 3.0 * ts_xyz_xxy[i] * gfe2_0 + 3.0 * gr_xyz_xxy[i] * gfe_0 + ts_xyz_xxxy[i] * gfe_0 * gc_x[i] + gr_xyz_xxxy[i] * gc_x[i];

        grr_x_xyz_xxxz[i] = ts_yz_xxxz[i] * gfe2_0 + gr_yz_xxxz[i] * gfe_0 + 3.0 * ts_xyz_xxz[i] * gfe2_0 + 3.0 * gr_xyz_xxz[i] * gfe_0 + ts_xyz_xxxz[i] * gfe_0 * gc_x[i] + gr_xyz_xxxz[i] * gc_x[i];

        grr_x_xyz_xxyy[i] = ts_yz_xxyy[i] * gfe2_0 + gr_yz_xxyy[i] * gfe_0 + 2.0 * ts_xyz_xyy[i] * gfe2_0 + 2.0 * gr_xyz_xyy[i] * gfe_0 + ts_xyz_xxyy[i] * gfe_0 * gc_x[i] + gr_xyz_xxyy[i] * gc_x[i];

        grr_x_xyz_xxyz[i] = ts_yz_xxyz[i] * gfe2_0 + gr_yz_xxyz[i] * gfe_0 + 2.0 * ts_xyz_xyz[i] * gfe2_0 + 2.0 * gr_xyz_xyz[i] * gfe_0 + ts_xyz_xxyz[i] * gfe_0 * gc_x[i] + gr_xyz_xxyz[i] * gc_x[i];

        grr_x_xyz_xxzz[i] = ts_yz_xxzz[i] * gfe2_0 + gr_yz_xxzz[i] * gfe_0 + 2.0 * ts_xyz_xzz[i] * gfe2_0 + 2.0 * gr_xyz_xzz[i] * gfe_0 + ts_xyz_xxzz[i] * gfe_0 * gc_x[i] + gr_xyz_xxzz[i] * gc_x[i];

        grr_x_xyz_xyyy[i] = ts_yz_xyyy[i] * gfe2_0 + gr_yz_xyyy[i] * gfe_0 + ts_xyz_yyy[i] * gfe2_0 + gr_xyz_yyy[i] * gfe_0 + ts_xyz_xyyy[i] * gfe_0 * gc_x[i] + gr_xyz_xyyy[i] * gc_x[i];

        grr_x_xyz_xyyz[i] = ts_yz_xyyz[i] * gfe2_0 + gr_yz_xyyz[i] * gfe_0 + ts_xyz_yyz[i] * gfe2_0 + gr_xyz_yyz[i] * gfe_0 + ts_xyz_xyyz[i] * gfe_0 * gc_x[i] + gr_xyz_xyyz[i] * gc_x[i];

        grr_x_xyz_xyzz[i] = ts_yz_xyzz[i] * gfe2_0 + gr_yz_xyzz[i] * gfe_0 + ts_xyz_yzz[i] * gfe2_0 + gr_xyz_yzz[i] * gfe_0 + ts_xyz_xyzz[i] * gfe_0 * gc_x[i] + gr_xyz_xyzz[i] * gc_x[i];

        grr_x_xyz_xzzz[i] = ts_yz_xzzz[i] * gfe2_0 + gr_yz_xzzz[i] * gfe_0 + ts_xyz_zzz[i] * gfe2_0 + gr_xyz_zzz[i] * gfe_0 + ts_xyz_xzzz[i] * gfe_0 * gc_x[i] + gr_xyz_xzzz[i] * gc_x[i];

        grr_x_xyz_yyyy[i] = ts_yz_yyyy[i] * gfe2_0 + gr_yz_yyyy[i] * gfe_0 + ts_xyz_yyyy[i] * gfe_0 * gc_x[i] + gr_xyz_yyyy[i] * gc_x[i];

        grr_x_xyz_yyyz[i] = ts_yz_yyyz[i] * gfe2_0 + gr_yz_yyyz[i] * gfe_0 + ts_xyz_yyyz[i] * gfe_0 * gc_x[i] + gr_xyz_yyyz[i] * gc_x[i];

        grr_x_xyz_yyzz[i] = ts_yz_yyzz[i] * gfe2_0 + gr_yz_yyzz[i] * gfe_0 + ts_xyz_yyzz[i] * gfe_0 * gc_x[i] + gr_xyz_yyzz[i] * gc_x[i];

        grr_x_xyz_yzzz[i] = ts_yz_yzzz[i] * gfe2_0 + gr_yz_yzzz[i] * gfe_0 + ts_xyz_yzzz[i] * gfe_0 * gc_x[i] + gr_xyz_yzzz[i] * gc_x[i];

        grr_x_xyz_zzzz[i] = ts_yz_zzzz[i] * gfe2_0 + gr_yz_zzzz[i] * gfe_0 + ts_xyz_zzzz[i] * gfe_0 * gc_x[i] + gr_xyz_zzzz[i] * gc_x[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto grr_x_xzz_xxxx = pbuffer.data(idx_gr_fg + 75);

    auto grr_x_xzz_xxxy = pbuffer.data(idx_gr_fg + 76);

    auto grr_x_xzz_xxxz = pbuffer.data(idx_gr_fg + 77);

    auto grr_x_xzz_xxyy = pbuffer.data(idx_gr_fg + 78);

    auto grr_x_xzz_xxyz = pbuffer.data(idx_gr_fg + 79);

    auto grr_x_xzz_xxzz = pbuffer.data(idx_gr_fg + 80);

    auto grr_x_xzz_xyyy = pbuffer.data(idx_gr_fg + 81);

    auto grr_x_xzz_xyyz = pbuffer.data(idx_gr_fg + 82);

    auto grr_x_xzz_xyzz = pbuffer.data(idx_gr_fg + 83);

    auto grr_x_xzz_xzzz = pbuffer.data(idx_gr_fg + 84);

    auto grr_x_xzz_yyyy = pbuffer.data(idx_gr_fg + 85);

    auto grr_x_xzz_yyyz = pbuffer.data(idx_gr_fg + 86);

    auto grr_x_xzz_yyzz = pbuffer.data(idx_gr_fg + 87);

    auto grr_x_xzz_yzzz = pbuffer.data(idx_gr_fg + 88);

    auto grr_x_xzz_zzzz = pbuffer.data(idx_gr_fg + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xxx, gr_xzz_xxxx, gr_xzz_xxxy, gr_xzz_xxxz, gr_xzz_xxy, gr_xzz_xxyy, gr_xzz_xxyz, gr_xzz_xxz, gr_xzz_xxzz, gr_xzz_xyy, gr_xzz_xyyy, gr_xzz_xyyz, gr_xzz_xyz, gr_xzz_xyzz, gr_xzz_xzz, gr_xzz_xzzz, gr_xzz_yyy, gr_xzz_yyyy, gr_xzz_yyyz, gr_xzz_yyz, gr_xzz_yyzz, gr_xzz_yzz, gr_xzz_yzzz, gr_xzz_zzz, gr_xzz_zzzz, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxzz, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyzz, gr_zz_xzzz, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyzz, gr_zz_yzzz, gr_zz_zzzz, grr_x_xzz_xxxx, grr_x_xzz_xxxy, grr_x_xzz_xxxz, grr_x_xzz_xxyy, grr_x_xzz_xxyz, grr_x_xzz_xxzz, grr_x_xzz_xyyy, grr_x_xzz_xyyz, grr_x_xzz_xyzz, grr_x_xzz_xzzz, grr_x_xzz_yyyy, grr_x_xzz_yyyz, grr_x_xzz_yyzz, grr_x_xzz_yzzz, grr_x_xzz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xzz_xxxx[i] = ts_zz_xxxx[i] * gfe2_0 + gr_zz_xxxx[i] * gfe_0 + 4.0 * ts_xzz_xxx[i] * gfe2_0 + 4.0 * gr_xzz_xxx[i] * gfe_0 + ts_xzz_xxxx[i] * gfe_0 * gc_x[i] + gr_xzz_xxxx[i] * gc_x[i];

        grr_x_xzz_xxxy[i] = ts_zz_xxxy[i] * gfe2_0 + gr_zz_xxxy[i] * gfe_0 + 3.0 * ts_xzz_xxy[i] * gfe2_0 + 3.0 * gr_xzz_xxy[i] * gfe_0 + ts_xzz_xxxy[i] * gfe_0 * gc_x[i] + gr_xzz_xxxy[i] * gc_x[i];

        grr_x_xzz_xxxz[i] = ts_zz_xxxz[i] * gfe2_0 + gr_zz_xxxz[i] * gfe_0 + 3.0 * ts_xzz_xxz[i] * gfe2_0 + 3.0 * gr_xzz_xxz[i] * gfe_0 + ts_xzz_xxxz[i] * gfe_0 * gc_x[i] + gr_xzz_xxxz[i] * gc_x[i];

        grr_x_xzz_xxyy[i] = ts_zz_xxyy[i] * gfe2_0 + gr_zz_xxyy[i] * gfe_0 + 2.0 * ts_xzz_xyy[i] * gfe2_0 + 2.0 * gr_xzz_xyy[i] * gfe_0 + ts_xzz_xxyy[i] * gfe_0 * gc_x[i] + gr_xzz_xxyy[i] * gc_x[i];

        grr_x_xzz_xxyz[i] = ts_zz_xxyz[i] * gfe2_0 + gr_zz_xxyz[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe2_0 + 2.0 * gr_xzz_xyz[i] * gfe_0 + ts_xzz_xxyz[i] * gfe_0 * gc_x[i] + gr_xzz_xxyz[i] * gc_x[i];

        grr_x_xzz_xxzz[i] = ts_zz_xxzz[i] * gfe2_0 + gr_zz_xxzz[i] * gfe_0 + 2.0 * ts_xzz_xzz[i] * gfe2_0 + 2.0 * gr_xzz_xzz[i] * gfe_0 + ts_xzz_xxzz[i] * gfe_0 * gc_x[i] + gr_xzz_xxzz[i] * gc_x[i];

        grr_x_xzz_xyyy[i] = ts_zz_xyyy[i] * gfe2_0 + gr_zz_xyyy[i] * gfe_0 + ts_xzz_yyy[i] * gfe2_0 + gr_xzz_yyy[i] * gfe_0 + ts_xzz_xyyy[i] * gfe_0 * gc_x[i] + gr_xzz_xyyy[i] * gc_x[i];

        grr_x_xzz_xyyz[i] = ts_zz_xyyz[i] * gfe2_0 + gr_zz_xyyz[i] * gfe_0 + ts_xzz_yyz[i] * gfe2_0 + gr_xzz_yyz[i] * gfe_0 + ts_xzz_xyyz[i] * gfe_0 * gc_x[i] + gr_xzz_xyyz[i] * gc_x[i];

        grr_x_xzz_xyzz[i] = ts_zz_xyzz[i] * gfe2_0 + gr_zz_xyzz[i] * gfe_0 + ts_xzz_yzz[i] * gfe2_0 + gr_xzz_yzz[i] * gfe_0 + ts_xzz_xyzz[i] * gfe_0 * gc_x[i] + gr_xzz_xyzz[i] * gc_x[i];

        grr_x_xzz_xzzz[i] = ts_zz_xzzz[i] * gfe2_0 + gr_zz_xzzz[i] * gfe_0 + ts_xzz_zzz[i] * gfe2_0 + gr_xzz_zzz[i] * gfe_0 + ts_xzz_xzzz[i] * gfe_0 * gc_x[i] + gr_xzz_xzzz[i] * gc_x[i];

        grr_x_xzz_yyyy[i] = ts_zz_yyyy[i] * gfe2_0 + gr_zz_yyyy[i] * gfe_0 + ts_xzz_yyyy[i] * gfe_0 * gc_x[i] + gr_xzz_yyyy[i] * gc_x[i];

        grr_x_xzz_yyyz[i] = ts_zz_yyyz[i] * gfe2_0 + gr_zz_yyyz[i] * gfe_0 + ts_xzz_yyyz[i] * gfe_0 * gc_x[i] + gr_xzz_yyyz[i] * gc_x[i];

        grr_x_xzz_yyzz[i] = ts_zz_yyzz[i] * gfe2_0 + gr_zz_yyzz[i] * gfe_0 + ts_xzz_yyzz[i] * gfe_0 * gc_x[i] + gr_xzz_yyzz[i] * gc_x[i];

        grr_x_xzz_yzzz[i] = ts_zz_yzzz[i] * gfe2_0 + gr_zz_yzzz[i] * gfe_0 + ts_xzz_yzzz[i] * gfe_0 * gc_x[i] + gr_xzz_yzzz[i] * gc_x[i];

        grr_x_xzz_zzzz[i] = ts_zz_zzzz[i] * gfe2_0 + gr_zz_zzzz[i] * gfe_0 + ts_xzz_zzzz[i] * gfe_0 * gc_x[i] + gr_xzz_zzzz[i] * gc_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto grr_x_yyy_xxxx = pbuffer.data(idx_gr_fg + 90);

    auto grr_x_yyy_xxxy = pbuffer.data(idx_gr_fg + 91);

    auto grr_x_yyy_xxxz = pbuffer.data(idx_gr_fg + 92);

    auto grr_x_yyy_xxyy = pbuffer.data(idx_gr_fg + 93);

    auto grr_x_yyy_xxyz = pbuffer.data(idx_gr_fg + 94);

    auto grr_x_yyy_xxzz = pbuffer.data(idx_gr_fg + 95);

    auto grr_x_yyy_xyyy = pbuffer.data(idx_gr_fg + 96);

    auto grr_x_yyy_xyyz = pbuffer.data(idx_gr_fg + 97);

    auto grr_x_yyy_xyzz = pbuffer.data(idx_gr_fg + 98);

    auto grr_x_yyy_xzzz = pbuffer.data(idx_gr_fg + 99);

    auto grr_x_yyy_yyyy = pbuffer.data(idx_gr_fg + 100);

    auto grr_x_yyy_yyyz = pbuffer.data(idx_gr_fg + 101);

    auto grr_x_yyy_yyzz = pbuffer.data(idx_gr_fg + 102);

    auto grr_x_yyy_yzzz = pbuffer.data(idx_gr_fg + 103);

    auto grr_x_yyy_zzzz = pbuffer.data(idx_gr_fg + 104);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xxx, gr_yyy_xxxx, gr_yyy_xxxy, gr_yyy_xxxz, gr_yyy_xxy, gr_yyy_xxyy, gr_yyy_xxyz, gr_yyy_xxz, gr_yyy_xxzz, gr_yyy_xyy, gr_yyy_xyyy, gr_yyy_xyyz, gr_yyy_xyz, gr_yyy_xyzz, gr_yyy_xzz, gr_yyy_xzzz, gr_yyy_yyy, gr_yyy_yyyy, gr_yyy_yyyz, gr_yyy_yyz, gr_yyy_yyzz, gr_yyy_yzz, gr_yyy_yzzz, gr_yyy_zzz, gr_yyy_zzzz, grr_x_yyy_xxxx, grr_x_yyy_xxxy, grr_x_yyy_xxxz, grr_x_yyy_xxyy, grr_x_yyy_xxyz, grr_x_yyy_xxzz, grr_x_yyy_xyyy, grr_x_yyy_xyyz, grr_x_yyy_xyzz, grr_x_yyy_xzzz, grr_x_yyy_yyyy, grr_x_yyy_yyyz, grr_x_yyy_yyzz, grr_x_yyy_yzzz, grr_x_yyy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyy_xxxx[i] = 4.0 * ts_yyy_xxx[i] * gfe2_0 + 4.0 * gr_yyy_xxx[i] * gfe_0 + ts_yyy_xxxx[i] * gfe_0 * gc_x[i] + gr_yyy_xxxx[i] * gc_x[i];

        grr_x_yyy_xxxy[i] = 3.0 * ts_yyy_xxy[i] * gfe2_0 + 3.0 * gr_yyy_xxy[i] * gfe_0 + ts_yyy_xxxy[i] * gfe_0 * gc_x[i] + gr_yyy_xxxy[i] * gc_x[i];

        grr_x_yyy_xxxz[i] = 3.0 * ts_yyy_xxz[i] * gfe2_0 + 3.0 * gr_yyy_xxz[i] * gfe_0 + ts_yyy_xxxz[i] * gfe_0 * gc_x[i] + gr_yyy_xxxz[i] * gc_x[i];

        grr_x_yyy_xxyy[i] = 2.0 * ts_yyy_xyy[i] * gfe2_0 + 2.0 * gr_yyy_xyy[i] * gfe_0 + ts_yyy_xxyy[i] * gfe_0 * gc_x[i] + gr_yyy_xxyy[i] * gc_x[i];

        grr_x_yyy_xxyz[i] = 2.0 * ts_yyy_xyz[i] * gfe2_0 + 2.0 * gr_yyy_xyz[i] * gfe_0 + ts_yyy_xxyz[i] * gfe_0 * gc_x[i] + gr_yyy_xxyz[i] * gc_x[i];

        grr_x_yyy_xxzz[i] = 2.0 * ts_yyy_xzz[i] * gfe2_0 + 2.0 * gr_yyy_xzz[i] * gfe_0 + ts_yyy_xxzz[i] * gfe_0 * gc_x[i] + gr_yyy_xxzz[i] * gc_x[i];

        grr_x_yyy_xyyy[i] = ts_yyy_yyy[i] * gfe2_0 + gr_yyy_yyy[i] * gfe_0 + ts_yyy_xyyy[i] * gfe_0 * gc_x[i] + gr_yyy_xyyy[i] * gc_x[i];

        grr_x_yyy_xyyz[i] = ts_yyy_yyz[i] * gfe2_0 + gr_yyy_yyz[i] * gfe_0 + ts_yyy_xyyz[i] * gfe_0 * gc_x[i] + gr_yyy_xyyz[i] * gc_x[i];

        grr_x_yyy_xyzz[i] = ts_yyy_yzz[i] * gfe2_0 + gr_yyy_yzz[i] * gfe_0 + ts_yyy_xyzz[i] * gfe_0 * gc_x[i] + gr_yyy_xyzz[i] * gc_x[i];

        grr_x_yyy_xzzz[i] = ts_yyy_zzz[i] * gfe2_0 + gr_yyy_zzz[i] * gfe_0 + ts_yyy_xzzz[i] * gfe_0 * gc_x[i] + gr_yyy_xzzz[i] * gc_x[i];

        grr_x_yyy_yyyy[i] = ts_yyy_yyyy[i] * gfe_0 * gc_x[i] + gr_yyy_yyyy[i] * gc_x[i];

        grr_x_yyy_yyyz[i] = ts_yyy_yyyz[i] * gfe_0 * gc_x[i] + gr_yyy_yyyz[i] * gc_x[i];

        grr_x_yyy_yyzz[i] = ts_yyy_yyzz[i] * gfe_0 * gc_x[i] + gr_yyy_yyzz[i] * gc_x[i];

        grr_x_yyy_yzzz[i] = ts_yyy_yzzz[i] * gfe_0 * gc_x[i] + gr_yyy_yzzz[i] * gc_x[i];

        grr_x_yyy_zzzz[i] = ts_yyy_zzzz[i] * gfe_0 * gc_x[i] + gr_yyy_zzzz[i] * gc_x[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto grr_x_yyz_xxxx = pbuffer.data(idx_gr_fg + 105);

    auto grr_x_yyz_xxxy = pbuffer.data(idx_gr_fg + 106);

    auto grr_x_yyz_xxxz = pbuffer.data(idx_gr_fg + 107);

    auto grr_x_yyz_xxyy = pbuffer.data(idx_gr_fg + 108);

    auto grr_x_yyz_xxyz = pbuffer.data(idx_gr_fg + 109);

    auto grr_x_yyz_xxzz = pbuffer.data(idx_gr_fg + 110);

    auto grr_x_yyz_xyyy = pbuffer.data(idx_gr_fg + 111);

    auto grr_x_yyz_xyyz = pbuffer.data(idx_gr_fg + 112);

    auto grr_x_yyz_xyzz = pbuffer.data(idx_gr_fg + 113);

    auto grr_x_yyz_xzzz = pbuffer.data(idx_gr_fg + 114);

    auto grr_x_yyz_yyyy = pbuffer.data(idx_gr_fg + 115);

    auto grr_x_yyz_yyyz = pbuffer.data(idx_gr_fg + 116);

    auto grr_x_yyz_yyzz = pbuffer.data(idx_gr_fg + 117);

    auto grr_x_yyz_yzzz = pbuffer.data(idx_gr_fg + 118);

    auto grr_x_yyz_zzzz = pbuffer.data(idx_gr_fg + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xxx, gr_yyz_xxxx, gr_yyz_xxxy, gr_yyz_xxxz, gr_yyz_xxy, gr_yyz_xxyy, gr_yyz_xxyz, gr_yyz_xxz, gr_yyz_xxzz, gr_yyz_xyy, gr_yyz_xyyy, gr_yyz_xyyz, gr_yyz_xyz, gr_yyz_xyzz, gr_yyz_xzz, gr_yyz_xzzz, gr_yyz_yyy, gr_yyz_yyyy, gr_yyz_yyyz, gr_yyz_yyz, gr_yyz_yyzz, gr_yyz_yzz, gr_yyz_yzzz, gr_yyz_zzz, gr_yyz_zzzz, grr_x_yyz_xxxx, grr_x_yyz_xxxy, grr_x_yyz_xxxz, grr_x_yyz_xxyy, grr_x_yyz_xxyz, grr_x_yyz_xxzz, grr_x_yyz_xyyy, grr_x_yyz_xyyz, grr_x_yyz_xyzz, grr_x_yyz_xzzz, grr_x_yyz_yyyy, grr_x_yyz_yyyz, grr_x_yyz_yyzz, grr_x_yyz_yzzz, grr_x_yyz_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyz_xxxx[i] = 4.0 * ts_yyz_xxx[i] * gfe2_0 + 4.0 * gr_yyz_xxx[i] * gfe_0 + ts_yyz_xxxx[i] * gfe_0 * gc_x[i] + gr_yyz_xxxx[i] * gc_x[i];

        grr_x_yyz_xxxy[i] = 3.0 * ts_yyz_xxy[i] * gfe2_0 + 3.0 * gr_yyz_xxy[i] * gfe_0 + ts_yyz_xxxy[i] * gfe_0 * gc_x[i] + gr_yyz_xxxy[i] * gc_x[i];

        grr_x_yyz_xxxz[i] = 3.0 * ts_yyz_xxz[i] * gfe2_0 + 3.0 * gr_yyz_xxz[i] * gfe_0 + ts_yyz_xxxz[i] * gfe_0 * gc_x[i] + gr_yyz_xxxz[i] * gc_x[i];

        grr_x_yyz_xxyy[i] = 2.0 * ts_yyz_xyy[i] * gfe2_0 + 2.0 * gr_yyz_xyy[i] * gfe_0 + ts_yyz_xxyy[i] * gfe_0 * gc_x[i] + gr_yyz_xxyy[i] * gc_x[i];

        grr_x_yyz_xxyz[i] = 2.0 * ts_yyz_xyz[i] * gfe2_0 + 2.0 * gr_yyz_xyz[i] * gfe_0 + ts_yyz_xxyz[i] * gfe_0 * gc_x[i] + gr_yyz_xxyz[i] * gc_x[i];

        grr_x_yyz_xxzz[i] = 2.0 * ts_yyz_xzz[i] * gfe2_0 + 2.0 * gr_yyz_xzz[i] * gfe_0 + ts_yyz_xxzz[i] * gfe_0 * gc_x[i] + gr_yyz_xxzz[i] * gc_x[i];

        grr_x_yyz_xyyy[i] = ts_yyz_yyy[i] * gfe2_0 + gr_yyz_yyy[i] * gfe_0 + ts_yyz_xyyy[i] * gfe_0 * gc_x[i] + gr_yyz_xyyy[i] * gc_x[i];

        grr_x_yyz_xyyz[i] = ts_yyz_yyz[i] * gfe2_0 + gr_yyz_yyz[i] * gfe_0 + ts_yyz_xyyz[i] * gfe_0 * gc_x[i] + gr_yyz_xyyz[i] * gc_x[i];

        grr_x_yyz_xyzz[i] = ts_yyz_yzz[i] * gfe2_0 + gr_yyz_yzz[i] * gfe_0 + ts_yyz_xyzz[i] * gfe_0 * gc_x[i] + gr_yyz_xyzz[i] * gc_x[i];

        grr_x_yyz_xzzz[i] = ts_yyz_zzz[i] * gfe2_0 + gr_yyz_zzz[i] * gfe_0 + ts_yyz_xzzz[i] * gfe_0 * gc_x[i] + gr_yyz_xzzz[i] * gc_x[i];

        grr_x_yyz_yyyy[i] = ts_yyz_yyyy[i] * gfe_0 * gc_x[i] + gr_yyz_yyyy[i] * gc_x[i];

        grr_x_yyz_yyyz[i] = ts_yyz_yyyz[i] * gfe_0 * gc_x[i] + gr_yyz_yyyz[i] * gc_x[i];

        grr_x_yyz_yyzz[i] = ts_yyz_yyzz[i] * gfe_0 * gc_x[i] + gr_yyz_yyzz[i] * gc_x[i];

        grr_x_yyz_yzzz[i] = ts_yyz_yzzz[i] * gfe_0 * gc_x[i] + gr_yyz_yzzz[i] * gc_x[i];

        grr_x_yyz_zzzz[i] = ts_yyz_zzzz[i] * gfe_0 * gc_x[i] + gr_yyz_zzzz[i] * gc_x[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto grr_x_yzz_xxxx = pbuffer.data(idx_gr_fg + 120);

    auto grr_x_yzz_xxxy = pbuffer.data(idx_gr_fg + 121);

    auto grr_x_yzz_xxxz = pbuffer.data(idx_gr_fg + 122);

    auto grr_x_yzz_xxyy = pbuffer.data(idx_gr_fg + 123);

    auto grr_x_yzz_xxyz = pbuffer.data(idx_gr_fg + 124);

    auto grr_x_yzz_xxzz = pbuffer.data(idx_gr_fg + 125);

    auto grr_x_yzz_xyyy = pbuffer.data(idx_gr_fg + 126);

    auto grr_x_yzz_xyyz = pbuffer.data(idx_gr_fg + 127);

    auto grr_x_yzz_xyzz = pbuffer.data(idx_gr_fg + 128);

    auto grr_x_yzz_xzzz = pbuffer.data(idx_gr_fg + 129);

    auto grr_x_yzz_yyyy = pbuffer.data(idx_gr_fg + 130);

    auto grr_x_yzz_yyyz = pbuffer.data(idx_gr_fg + 131);

    auto grr_x_yzz_yyzz = pbuffer.data(idx_gr_fg + 132);

    auto grr_x_yzz_yzzz = pbuffer.data(idx_gr_fg + 133);

    auto grr_x_yzz_zzzz = pbuffer.data(idx_gr_fg + 134);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xxx, gr_yzz_xxxx, gr_yzz_xxxy, gr_yzz_xxxz, gr_yzz_xxy, gr_yzz_xxyy, gr_yzz_xxyz, gr_yzz_xxz, gr_yzz_xxzz, gr_yzz_xyy, gr_yzz_xyyy, gr_yzz_xyyz, gr_yzz_xyz, gr_yzz_xyzz, gr_yzz_xzz, gr_yzz_xzzz, gr_yzz_yyy, gr_yzz_yyyy, gr_yzz_yyyz, gr_yzz_yyz, gr_yzz_yyzz, gr_yzz_yzz, gr_yzz_yzzz, gr_yzz_zzz, gr_yzz_zzzz, grr_x_yzz_xxxx, grr_x_yzz_xxxy, grr_x_yzz_xxxz, grr_x_yzz_xxyy, grr_x_yzz_xxyz, grr_x_yzz_xxzz, grr_x_yzz_xyyy, grr_x_yzz_xyyz, grr_x_yzz_xyzz, grr_x_yzz_xzzz, grr_x_yzz_yyyy, grr_x_yzz_yyyz, grr_x_yzz_yyzz, grr_x_yzz_yzzz, grr_x_yzz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yzz_xxxx[i] = 4.0 * ts_yzz_xxx[i] * gfe2_0 + 4.0 * gr_yzz_xxx[i] * gfe_0 + ts_yzz_xxxx[i] * gfe_0 * gc_x[i] + gr_yzz_xxxx[i] * gc_x[i];

        grr_x_yzz_xxxy[i] = 3.0 * ts_yzz_xxy[i] * gfe2_0 + 3.0 * gr_yzz_xxy[i] * gfe_0 + ts_yzz_xxxy[i] * gfe_0 * gc_x[i] + gr_yzz_xxxy[i] * gc_x[i];

        grr_x_yzz_xxxz[i] = 3.0 * ts_yzz_xxz[i] * gfe2_0 + 3.0 * gr_yzz_xxz[i] * gfe_0 + ts_yzz_xxxz[i] * gfe_0 * gc_x[i] + gr_yzz_xxxz[i] * gc_x[i];

        grr_x_yzz_xxyy[i] = 2.0 * ts_yzz_xyy[i] * gfe2_0 + 2.0 * gr_yzz_xyy[i] * gfe_0 + ts_yzz_xxyy[i] * gfe_0 * gc_x[i] + gr_yzz_xxyy[i] * gc_x[i];

        grr_x_yzz_xxyz[i] = 2.0 * ts_yzz_xyz[i] * gfe2_0 + 2.0 * gr_yzz_xyz[i] * gfe_0 + ts_yzz_xxyz[i] * gfe_0 * gc_x[i] + gr_yzz_xxyz[i] * gc_x[i];

        grr_x_yzz_xxzz[i] = 2.0 * ts_yzz_xzz[i] * gfe2_0 + 2.0 * gr_yzz_xzz[i] * gfe_0 + ts_yzz_xxzz[i] * gfe_0 * gc_x[i] + gr_yzz_xxzz[i] * gc_x[i];

        grr_x_yzz_xyyy[i] = ts_yzz_yyy[i] * gfe2_0 + gr_yzz_yyy[i] * gfe_0 + ts_yzz_xyyy[i] * gfe_0 * gc_x[i] + gr_yzz_xyyy[i] * gc_x[i];

        grr_x_yzz_xyyz[i] = ts_yzz_yyz[i] * gfe2_0 + gr_yzz_yyz[i] * gfe_0 + ts_yzz_xyyz[i] * gfe_0 * gc_x[i] + gr_yzz_xyyz[i] * gc_x[i];

        grr_x_yzz_xyzz[i] = ts_yzz_yzz[i] * gfe2_0 + gr_yzz_yzz[i] * gfe_0 + ts_yzz_xyzz[i] * gfe_0 * gc_x[i] + gr_yzz_xyzz[i] * gc_x[i];

        grr_x_yzz_xzzz[i] = ts_yzz_zzz[i] * gfe2_0 + gr_yzz_zzz[i] * gfe_0 + ts_yzz_xzzz[i] * gfe_0 * gc_x[i] + gr_yzz_xzzz[i] * gc_x[i];

        grr_x_yzz_yyyy[i] = ts_yzz_yyyy[i] * gfe_0 * gc_x[i] + gr_yzz_yyyy[i] * gc_x[i];

        grr_x_yzz_yyyz[i] = ts_yzz_yyyz[i] * gfe_0 * gc_x[i] + gr_yzz_yyyz[i] * gc_x[i];

        grr_x_yzz_yyzz[i] = ts_yzz_yyzz[i] * gfe_0 * gc_x[i] + gr_yzz_yyzz[i] * gc_x[i];

        grr_x_yzz_yzzz[i] = ts_yzz_yzzz[i] * gfe_0 * gc_x[i] + gr_yzz_yzzz[i] * gc_x[i];

        grr_x_yzz_zzzz[i] = ts_yzz_zzzz[i] * gfe_0 * gc_x[i] + gr_yzz_zzzz[i] * gc_x[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto grr_x_zzz_xxxx = pbuffer.data(idx_gr_fg + 135);

    auto grr_x_zzz_xxxy = pbuffer.data(idx_gr_fg + 136);

    auto grr_x_zzz_xxxz = pbuffer.data(idx_gr_fg + 137);

    auto grr_x_zzz_xxyy = pbuffer.data(idx_gr_fg + 138);

    auto grr_x_zzz_xxyz = pbuffer.data(idx_gr_fg + 139);

    auto grr_x_zzz_xxzz = pbuffer.data(idx_gr_fg + 140);

    auto grr_x_zzz_xyyy = pbuffer.data(idx_gr_fg + 141);

    auto grr_x_zzz_xyyz = pbuffer.data(idx_gr_fg + 142);

    auto grr_x_zzz_xyzz = pbuffer.data(idx_gr_fg + 143);

    auto grr_x_zzz_xzzz = pbuffer.data(idx_gr_fg + 144);

    auto grr_x_zzz_yyyy = pbuffer.data(idx_gr_fg + 145);

    auto grr_x_zzz_yyyz = pbuffer.data(idx_gr_fg + 146);

    auto grr_x_zzz_yyzz = pbuffer.data(idx_gr_fg + 147);

    auto grr_x_zzz_yzzz = pbuffer.data(idx_gr_fg + 148);

    auto grr_x_zzz_zzzz = pbuffer.data(idx_gr_fg + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xxx, gr_zzz_xxxx, gr_zzz_xxxy, gr_zzz_xxxz, gr_zzz_xxy, gr_zzz_xxyy, gr_zzz_xxyz, gr_zzz_xxz, gr_zzz_xxzz, gr_zzz_xyy, gr_zzz_xyyy, gr_zzz_xyyz, gr_zzz_xyz, gr_zzz_xyzz, gr_zzz_xzz, gr_zzz_xzzz, gr_zzz_yyy, gr_zzz_yyyy, gr_zzz_yyyz, gr_zzz_yyz, gr_zzz_yyzz, gr_zzz_yzz, gr_zzz_yzzz, gr_zzz_zzz, gr_zzz_zzzz, grr_x_zzz_xxxx, grr_x_zzz_xxxy, grr_x_zzz_xxxz, grr_x_zzz_xxyy, grr_x_zzz_xxyz, grr_x_zzz_xxzz, grr_x_zzz_xyyy, grr_x_zzz_xyyz, grr_x_zzz_xyzz, grr_x_zzz_xzzz, grr_x_zzz_yyyy, grr_x_zzz_yyyz, grr_x_zzz_yyzz, grr_x_zzz_yzzz, grr_x_zzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zzz_xxxx[i] = 4.0 * ts_zzz_xxx[i] * gfe2_0 + 4.0 * gr_zzz_xxx[i] * gfe_0 + ts_zzz_xxxx[i] * gfe_0 * gc_x[i] + gr_zzz_xxxx[i] * gc_x[i];

        grr_x_zzz_xxxy[i] = 3.0 * ts_zzz_xxy[i] * gfe2_0 + 3.0 * gr_zzz_xxy[i] * gfe_0 + ts_zzz_xxxy[i] * gfe_0 * gc_x[i] + gr_zzz_xxxy[i] * gc_x[i];

        grr_x_zzz_xxxz[i] = 3.0 * ts_zzz_xxz[i] * gfe2_0 + 3.0 * gr_zzz_xxz[i] * gfe_0 + ts_zzz_xxxz[i] * gfe_0 * gc_x[i] + gr_zzz_xxxz[i] * gc_x[i];

        grr_x_zzz_xxyy[i] = 2.0 * ts_zzz_xyy[i] * gfe2_0 + 2.0 * gr_zzz_xyy[i] * gfe_0 + ts_zzz_xxyy[i] * gfe_0 * gc_x[i] + gr_zzz_xxyy[i] * gc_x[i];

        grr_x_zzz_xxyz[i] = 2.0 * ts_zzz_xyz[i] * gfe2_0 + 2.0 * gr_zzz_xyz[i] * gfe_0 + ts_zzz_xxyz[i] * gfe_0 * gc_x[i] + gr_zzz_xxyz[i] * gc_x[i];

        grr_x_zzz_xxzz[i] = 2.0 * ts_zzz_xzz[i] * gfe2_0 + 2.0 * gr_zzz_xzz[i] * gfe_0 + ts_zzz_xxzz[i] * gfe_0 * gc_x[i] + gr_zzz_xxzz[i] * gc_x[i];

        grr_x_zzz_xyyy[i] = ts_zzz_yyy[i] * gfe2_0 + gr_zzz_yyy[i] * gfe_0 + ts_zzz_xyyy[i] * gfe_0 * gc_x[i] + gr_zzz_xyyy[i] * gc_x[i];

        grr_x_zzz_xyyz[i] = ts_zzz_yyz[i] * gfe2_0 + gr_zzz_yyz[i] * gfe_0 + ts_zzz_xyyz[i] * gfe_0 * gc_x[i] + gr_zzz_xyyz[i] * gc_x[i];

        grr_x_zzz_xyzz[i] = ts_zzz_yzz[i] * gfe2_0 + gr_zzz_yzz[i] * gfe_0 + ts_zzz_xyzz[i] * gfe_0 * gc_x[i] + gr_zzz_xyzz[i] * gc_x[i];

        grr_x_zzz_xzzz[i] = ts_zzz_zzz[i] * gfe2_0 + gr_zzz_zzz[i] * gfe_0 + ts_zzz_xzzz[i] * gfe_0 * gc_x[i] + gr_zzz_xzzz[i] * gc_x[i];

        grr_x_zzz_yyyy[i] = ts_zzz_yyyy[i] * gfe_0 * gc_x[i] + gr_zzz_yyyy[i] * gc_x[i];

        grr_x_zzz_yyyz[i] = ts_zzz_yyyz[i] * gfe_0 * gc_x[i] + gr_zzz_yyyz[i] * gc_x[i];

        grr_x_zzz_yyzz[i] = ts_zzz_yyzz[i] * gfe_0 * gc_x[i] + gr_zzz_yyzz[i] * gc_x[i];

        grr_x_zzz_yzzz[i] = ts_zzz_yzzz[i] * gfe_0 * gc_x[i] + gr_zzz_yzzz[i] * gc_x[i];

        grr_x_zzz_zzzz[i] = ts_zzz_zzzz[i] * gfe_0 * gc_x[i] + gr_zzz_zzzz[i] * gc_x[i];
    }

    // Set up 150-165 components of targeted buffer : FG

    auto grr_y_xxx_xxxx = pbuffer.data(idx_gr_fg + 150);

    auto grr_y_xxx_xxxy = pbuffer.data(idx_gr_fg + 151);

    auto grr_y_xxx_xxxz = pbuffer.data(idx_gr_fg + 152);

    auto grr_y_xxx_xxyy = pbuffer.data(idx_gr_fg + 153);

    auto grr_y_xxx_xxyz = pbuffer.data(idx_gr_fg + 154);

    auto grr_y_xxx_xxzz = pbuffer.data(idx_gr_fg + 155);

    auto grr_y_xxx_xyyy = pbuffer.data(idx_gr_fg + 156);

    auto grr_y_xxx_xyyz = pbuffer.data(idx_gr_fg + 157);

    auto grr_y_xxx_xyzz = pbuffer.data(idx_gr_fg + 158);

    auto grr_y_xxx_xzzz = pbuffer.data(idx_gr_fg + 159);

    auto grr_y_xxx_yyyy = pbuffer.data(idx_gr_fg + 160);

    auto grr_y_xxx_yyyz = pbuffer.data(idx_gr_fg + 161);

    auto grr_y_xxx_yyzz = pbuffer.data(idx_gr_fg + 162);

    auto grr_y_xxx_yzzz = pbuffer.data(idx_gr_fg + 163);

    auto grr_y_xxx_zzzz = pbuffer.data(idx_gr_fg + 164);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xxx, gr_xxx_xxxx, gr_xxx_xxxy, gr_xxx_xxxz, gr_xxx_xxy, gr_xxx_xxyy, gr_xxx_xxyz, gr_xxx_xxz, gr_xxx_xxzz, gr_xxx_xyy, gr_xxx_xyyy, gr_xxx_xyyz, gr_xxx_xyz, gr_xxx_xyzz, gr_xxx_xzz, gr_xxx_xzzz, gr_xxx_yyy, gr_xxx_yyyy, gr_xxx_yyyz, gr_xxx_yyz, gr_xxx_yyzz, gr_xxx_yzz, gr_xxx_yzzz, gr_xxx_zzz, gr_xxx_zzzz, grr_y_xxx_xxxx, grr_y_xxx_xxxy, grr_y_xxx_xxxz, grr_y_xxx_xxyy, grr_y_xxx_xxyz, grr_y_xxx_xxzz, grr_y_xxx_xyyy, grr_y_xxx_xyyz, grr_y_xxx_xyzz, grr_y_xxx_xzzz, grr_y_xxx_yyyy, grr_y_xxx_yyyz, grr_y_xxx_yyzz, grr_y_xxx_yzzz, grr_y_xxx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxx_xxxx[i] = ts_xxx_xxxx[i] * gfe_0 * gc_y[i] + gr_xxx_xxxx[i] * gc_y[i];

        grr_y_xxx_xxxy[i] = ts_xxx_xxx[i] * gfe2_0 + gr_xxx_xxx[i] * gfe_0 + ts_xxx_xxxy[i] * gfe_0 * gc_y[i] + gr_xxx_xxxy[i] * gc_y[i];

        grr_y_xxx_xxxz[i] = ts_xxx_xxxz[i] * gfe_0 * gc_y[i] + gr_xxx_xxxz[i] * gc_y[i];

        grr_y_xxx_xxyy[i] = 2.0 * ts_xxx_xxy[i] * gfe2_0 + 2.0 * gr_xxx_xxy[i] * gfe_0 + ts_xxx_xxyy[i] * gfe_0 * gc_y[i] + gr_xxx_xxyy[i] * gc_y[i];

        grr_y_xxx_xxyz[i] = ts_xxx_xxz[i] * gfe2_0 + gr_xxx_xxz[i] * gfe_0 + ts_xxx_xxyz[i] * gfe_0 * gc_y[i] + gr_xxx_xxyz[i] * gc_y[i];

        grr_y_xxx_xxzz[i] = ts_xxx_xxzz[i] * gfe_0 * gc_y[i] + gr_xxx_xxzz[i] * gc_y[i];

        grr_y_xxx_xyyy[i] = 3.0 * ts_xxx_xyy[i] * gfe2_0 + 3.0 * gr_xxx_xyy[i] * gfe_0 + ts_xxx_xyyy[i] * gfe_0 * gc_y[i] + gr_xxx_xyyy[i] * gc_y[i];

        grr_y_xxx_xyyz[i] = 2.0 * ts_xxx_xyz[i] * gfe2_0 + 2.0 * gr_xxx_xyz[i] * gfe_0 + ts_xxx_xyyz[i] * gfe_0 * gc_y[i] + gr_xxx_xyyz[i] * gc_y[i];

        grr_y_xxx_xyzz[i] = ts_xxx_xzz[i] * gfe2_0 + gr_xxx_xzz[i] * gfe_0 + ts_xxx_xyzz[i] * gfe_0 * gc_y[i] + gr_xxx_xyzz[i] * gc_y[i];

        grr_y_xxx_xzzz[i] = ts_xxx_xzzz[i] * gfe_0 * gc_y[i] + gr_xxx_xzzz[i] * gc_y[i];

        grr_y_xxx_yyyy[i] = 4.0 * ts_xxx_yyy[i] * gfe2_0 + 4.0 * gr_xxx_yyy[i] * gfe_0 + ts_xxx_yyyy[i] * gfe_0 * gc_y[i] + gr_xxx_yyyy[i] * gc_y[i];

        grr_y_xxx_yyyz[i] = 3.0 * ts_xxx_yyz[i] * gfe2_0 + 3.0 * gr_xxx_yyz[i] * gfe_0 + ts_xxx_yyyz[i] * gfe_0 * gc_y[i] + gr_xxx_yyyz[i] * gc_y[i];

        grr_y_xxx_yyzz[i] = 2.0 * ts_xxx_yzz[i] * gfe2_0 + 2.0 * gr_xxx_yzz[i] * gfe_0 + ts_xxx_yyzz[i] * gfe_0 * gc_y[i] + gr_xxx_yyzz[i] * gc_y[i];

        grr_y_xxx_yzzz[i] = ts_xxx_zzz[i] * gfe2_0 + gr_xxx_zzz[i] * gfe_0 + ts_xxx_yzzz[i] * gfe_0 * gc_y[i] + gr_xxx_yzzz[i] * gc_y[i];

        grr_y_xxx_zzzz[i] = ts_xxx_zzzz[i] * gfe_0 * gc_y[i] + gr_xxx_zzzz[i] * gc_y[i];
    }

    // Set up 165-180 components of targeted buffer : FG

    auto grr_y_xxy_xxxx = pbuffer.data(idx_gr_fg + 165);

    auto grr_y_xxy_xxxy = pbuffer.data(idx_gr_fg + 166);

    auto grr_y_xxy_xxxz = pbuffer.data(idx_gr_fg + 167);

    auto grr_y_xxy_xxyy = pbuffer.data(idx_gr_fg + 168);

    auto grr_y_xxy_xxyz = pbuffer.data(idx_gr_fg + 169);

    auto grr_y_xxy_xxzz = pbuffer.data(idx_gr_fg + 170);

    auto grr_y_xxy_xyyy = pbuffer.data(idx_gr_fg + 171);

    auto grr_y_xxy_xyyz = pbuffer.data(idx_gr_fg + 172);

    auto grr_y_xxy_xyzz = pbuffer.data(idx_gr_fg + 173);

    auto grr_y_xxy_xzzz = pbuffer.data(idx_gr_fg + 174);

    auto grr_y_xxy_yyyy = pbuffer.data(idx_gr_fg + 175);

    auto grr_y_xxy_yyyz = pbuffer.data(idx_gr_fg + 176);

    auto grr_y_xxy_yyzz = pbuffer.data(idx_gr_fg + 177);

    auto grr_y_xxy_yzzz = pbuffer.data(idx_gr_fg + 178);

    auto grr_y_xxy_zzzz = pbuffer.data(idx_gr_fg + 179);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxzz, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyzz, gr_xx_xzzz, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyzz, gr_xx_yzzz, gr_xx_zzzz, gr_xxy_xxx, gr_xxy_xxxx, gr_xxy_xxxy, gr_xxy_xxxz, gr_xxy_xxy, gr_xxy_xxyy, gr_xxy_xxyz, gr_xxy_xxz, gr_xxy_xxzz, gr_xxy_xyy, gr_xxy_xyyy, gr_xxy_xyyz, gr_xxy_xyz, gr_xxy_xyzz, gr_xxy_xzz, gr_xxy_xzzz, gr_xxy_yyy, gr_xxy_yyyy, gr_xxy_yyyz, gr_xxy_yyz, gr_xxy_yyzz, gr_xxy_yzz, gr_xxy_yzzz, gr_xxy_zzz, gr_xxy_zzzz, grr_y_xxy_xxxx, grr_y_xxy_xxxy, grr_y_xxy_xxxz, grr_y_xxy_xxyy, grr_y_xxy_xxyz, grr_y_xxy_xxzz, grr_y_xxy_xyyy, grr_y_xxy_xyyz, grr_y_xxy_xyzz, grr_y_xxy_xzzz, grr_y_xxy_yyyy, grr_y_xxy_yyyz, grr_y_xxy_yyzz, grr_y_xxy_yzzz, grr_y_xxy_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxy_xxxx[i] = ts_xx_xxxx[i] * gfe2_0 + gr_xx_xxxx[i] * gfe_0 + ts_xxy_xxxx[i] * gfe_0 * gc_y[i] + gr_xxy_xxxx[i] * gc_y[i];

        grr_y_xxy_xxxy[i] = ts_xx_xxxy[i] * gfe2_0 + gr_xx_xxxy[i] * gfe_0 + ts_xxy_xxx[i] * gfe2_0 + gr_xxy_xxx[i] * gfe_0 + ts_xxy_xxxy[i] * gfe_0 * gc_y[i] + gr_xxy_xxxy[i] * gc_y[i];

        grr_y_xxy_xxxz[i] = ts_xx_xxxz[i] * gfe2_0 + gr_xx_xxxz[i] * gfe_0 + ts_xxy_xxxz[i] * gfe_0 * gc_y[i] + gr_xxy_xxxz[i] * gc_y[i];

        grr_y_xxy_xxyy[i] = ts_xx_xxyy[i] * gfe2_0 + gr_xx_xxyy[i] * gfe_0 + 2.0 * ts_xxy_xxy[i] * gfe2_0 + 2.0 * gr_xxy_xxy[i] * gfe_0 + ts_xxy_xxyy[i] * gfe_0 * gc_y[i] + gr_xxy_xxyy[i] * gc_y[i];

        grr_y_xxy_xxyz[i] = ts_xx_xxyz[i] * gfe2_0 + gr_xx_xxyz[i] * gfe_0 + ts_xxy_xxz[i] * gfe2_0 + gr_xxy_xxz[i] * gfe_0 + ts_xxy_xxyz[i] * gfe_0 * gc_y[i] + gr_xxy_xxyz[i] * gc_y[i];

        grr_y_xxy_xxzz[i] = ts_xx_xxzz[i] * gfe2_0 + gr_xx_xxzz[i] * gfe_0 + ts_xxy_xxzz[i] * gfe_0 * gc_y[i] + gr_xxy_xxzz[i] * gc_y[i];

        grr_y_xxy_xyyy[i] = ts_xx_xyyy[i] * gfe2_0 + gr_xx_xyyy[i] * gfe_0 + 3.0 * ts_xxy_xyy[i] * gfe2_0 + 3.0 * gr_xxy_xyy[i] * gfe_0 + ts_xxy_xyyy[i] * gfe_0 * gc_y[i] + gr_xxy_xyyy[i] * gc_y[i];

        grr_y_xxy_xyyz[i] = ts_xx_xyyz[i] * gfe2_0 + gr_xx_xyyz[i] * gfe_0 + 2.0 * ts_xxy_xyz[i] * gfe2_0 + 2.0 * gr_xxy_xyz[i] * gfe_0 + ts_xxy_xyyz[i] * gfe_0 * gc_y[i] + gr_xxy_xyyz[i] * gc_y[i];

        grr_y_xxy_xyzz[i] = ts_xx_xyzz[i] * gfe2_0 + gr_xx_xyzz[i] * gfe_0 + ts_xxy_xzz[i] * gfe2_0 + gr_xxy_xzz[i] * gfe_0 + ts_xxy_xyzz[i] * gfe_0 * gc_y[i] + gr_xxy_xyzz[i] * gc_y[i];

        grr_y_xxy_xzzz[i] = ts_xx_xzzz[i] * gfe2_0 + gr_xx_xzzz[i] * gfe_0 + ts_xxy_xzzz[i] * gfe_0 * gc_y[i] + gr_xxy_xzzz[i] * gc_y[i];

        grr_y_xxy_yyyy[i] = ts_xx_yyyy[i] * gfe2_0 + gr_xx_yyyy[i] * gfe_0 + 4.0 * ts_xxy_yyy[i] * gfe2_0 + 4.0 * gr_xxy_yyy[i] * gfe_0 + ts_xxy_yyyy[i] * gfe_0 * gc_y[i] + gr_xxy_yyyy[i] * gc_y[i];

        grr_y_xxy_yyyz[i] = ts_xx_yyyz[i] * gfe2_0 + gr_xx_yyyz[i] * gfe_0 + 3.0 * ts_xxy_yyz[i] * gfe2_0 + 3.0 * gr_xxy_yyz[i] * gfe_0 + ts_xxy_yyyz[i] * gfe_0 * gc_y[i] + gr_xxy_yyyz[i] * gc_y[i];

        grr_y_xxy_yyzz[i] = ts_xx_yyzz[i] * gfe2_0 + gr_xx_yyzz[i] * gfe_0 + 2.0 * ts_xxy_yzz[i] * gfe2_0 + 2.0 * gr_xxy_yzz[i] * gfe_0 + ts_xxy_yyzz[i] * gfe_0 * gc_y[i] + gr_xxy_yyzz[i] * gc_y[i];

        grr_y_xxy_yzzz[i] = ts_xx_yzzz[i] * gfe2_0 + gr_xx_yzzz[i] * gfe_0 + ts_xxy_zzz[i] * gfe2_0 + gr_xxy_zzz[i] * gfe_0 + ts_xxy_yzzz[i] * gfe_0 * gc_y[i] + gr_xxy_yzzz[i] * gc_y[i];

        grr_y_xxy_zzzz[i] = ts_xx_zzzz[i] * gfe2_0 + gr_xx_zzzz[i] * gfe_0 + ts_xxy_zzzz[i] * gfe_0 * gc_y[i] + gr_xxy_zzzz[i] * gc_y[i];
    }

    // Set up 180-195 components of targeted buffer : FG

    auto grr_y_xxz_xxxx = pbuffer.data(idx_gr_fg + 180);

    auto grr_y_xxz_xxxy = pbuffer.data(idx_gr_fg + 181);

    auto grr_y_xxz_xxxz = pbuffer.data(idx_gr_fg + 182);

    auto grr_y_xxz_xxyy = pbuffer.data(idx_gr_fg + 183);

    auto grr_y_xxz_xxyz = pbuffer.data(idx_gr_fg + 184);

    auto grr_y_xxz_xxzz = pbuffer.data(idx_gr_fg + 185);

    auto grr_y_xxz_xyyy = pbuffer.data(idx_gr_fg + 186);

    auto grr_y_xxz_xyyz = pbuffer.data(idx_gr_fg + 187);

    auto grr_y_xxz_xyzz = pbuffer.data(idx_gr_fg + 188);

    auto grr_y_xxz_xzzz = pbuffer.data(idx_gr_fg + 189);

    auto grr_y_xxz_yyyy = pbuffer.data(idx_gr_fg + 190);

    auto grr_y_xxz_yyyz = pbuffer.data(idx_gr_fg + 191);

    auto grr_y_xxz_yyzz = pbuffer.data(idx_gr_fg + 192);

    auto grr_y_xxz_yzzz = pbuffer.data(idx_gr_fg + 193);

    auto grr_y_xxz_zzzz = pbuffer.data(idx_gr_fg + 194);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xxx, gr_xxz_xxxx, gr_xxz_xxxy, gr_xxz_xxxz, gr_xxz_xxy, gr_xxz_xxyy, gr_xxz_xxyz, gr_xxz_xxz, gr_xxz_xxzz, gr_xxz_xyy, gr_xxz_xyyy, gr_xxz_xyyz, gr_xxz_xyz, gr_xxz_xyzz, gr_xxz_xzz, gr_xxz_xzzz, gr_xxz_yyy, gr_xxz_yyyy, gr_xxz_yyyz, gr_xxz_yyz, gr_xxz_yyzz, gr_xxz_yzz, gr_xxz_yzzz, gr_xxz_zzz, gr_xxz_zzzz, grr_y_xxz_xxxx, grr_y_xxz_xxxy, grr_y_xxz_xxxz, grr_y_xxz_xxyy, grr_y_xxz_xxyz, grr_y_xxz_xxzz, grr_y_xxz_xyyy, grr_y_xxz_xyyz, grr_y_xxz_xyzz, grr_y_xxz_xzzz, grr_y_xxz_yyyy, grr_y_xxz_yyyz, grr_y_xxz_yyzz, grr_y_xxz_yzzz, grr_y_xxz_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxz_xxxx[i] = ts_xxz_xxxx[i] * gfe_0 * gc_y[i] + gr_xxz_xxxx[i] * gc_y[i];

        grr_y_xxz_xxxy[i] = ts_xxz_xxx[i] * gfe2_0 + gr_xxz_xxx[i] * gfe_0 + ts_xxz_xxxy[i] * gfe_0 * gc_y[i] + gr_xxz_xxxy[i] * gc_y[i];

        grr_y_xxz_xxxz[i] = ts_xxz_xxxz[i] * gfe_0 * gc_y[i] + gr_xxz_xxxz[i] * gc_y[i];

        grr_y_xxz_xxyy[i] = 2.0 * ts_xxz_xxy[i] * gfe2_0 + 2.0 * gr_xxz_xxy[i] * gfe_0 + ts_xxz_xxyy[i] * gfe_0 * gc_y[i] + gr_xxz_xxyy[i] * gc_y[i];

        grr_y_xxz_xxyz[i] = ts_xxz_xxz[i] * gfe2_0 + gr_xxz_xxz[i] * gfe_0 + ts_xxz_xxyz[i] * gfe_0 * gc_y[i] + gr_xxz_xxyz[i] * gc_y[i];

        grr_y_xxz_xxzz[i] = ts_xxz_xxzz[i] * gfe_0 * gc_y[i] + gr_xxz_xxzz[i] * gc_y[i];

        grr_y_xxz_xyyy[i] = 3.0 * ts_xxz_xyy[i] * gfe2_0 + 3.0 * gr_xxz_xyy[i] * gfe_0 + ts_xxz_xyyy[i] * gfe_0 * gc_y[i] + gr_xxz_xyyy[i] * gc_y[i];

        grr_y_xxz_xyyz[i] = 2.0 * ts_xxz_xyz[i] * gfe2_0 + 2.0 * gr_xxz_xyz[i] * gfe_0 + ts_xxz_xyyz[i] * gfe_0 * gc_y[i] + gr_xxz_xyyz[i] * gc_y[i];

        grr_y_xxz_xyzz[i] = ts_xxz_xzz[i] * gfe2_0 + gr_xxz_xzz[i] * gfe_0 + ts_xxz_xyzz[i] * gfe_0 * gc_y[i] + gr_xxz_xyzz[i] * gc_y[i];

        grr_y_xxz_xzzz[i] = ts_xxz_xzzz[i] * gfe_0 * gc_y[i] + gr_xxz_xzzz[i] * gc_y[i];

        grr_y_xxz_yyyy[i] = 4.0 * ts_xxz_yyy[i] * gfe2_0 + 4.0 * gr_xxz_yyy[i] * gfe_0 + ts_xxz_yyyy[i] * gfe_0 * gc_y[i] + gr_xxz_yyyy[i] * gc_y[i];

        grr_y_xxz_yyyz[i] = 3.0 * ts_xxz_yyz[i] * gfe2_0 + 3.0 * gr_xxz_yyz[i] * gfe_0 + ts_xxz_yyyz[i] * gfe_0 * gc_y[i] + gr_xxz_yyyz[i] * gc_y[i];

        grr_y_xxz_yyzz[i] = 2.0 * ts_xxz_yzz[i] * gfe2_0 + 2.0 * gr_xxz_yzz[i] * gfe_0 + ts_xxz_yyzz[i] * gfe_0 * gc_y[i] + gr_xxz_yyzz[i] * gc_y[i];

        grr_y_xxz_yzzz[i] = ts_xxz_zzz[i] * gfe2_0 + gr_xxz_zzz[i] * gfe_0 + ts_xxz_yzzz[i] * gfe_0 * gc_y[i] + gr_xxz_yzzz[i] * gc_y[i];

        grr_y_xxz_zzzz[i] = ts_xxz_zzzz[i] * gfe_0 * gc_y[i] + gr_xxz_zzzz[i] * gc_y[i];
    }

    // Set up 195-210 components of targeted buffer : FG

    auto grr_y_xyy_xxxx = pbuffer.data(idx_gr_fg + 195);

    auto grr_y_xyy_xxxy = pbuffer.data(idx_gr_fg + 196);

    auto grr_y_xyy_xxxz = pbuffer.data(idx_gr_fg + 197);

    auto grr_y_xyy_xxyy = pbuffer.data(idx_gr_fg + 198);

    auto grr_y_xyy_xxyz = pbuffer.data(idx_gr_fg + 199);

    auto grr_y_xyy_xxzz = pbuffer.data(idx_gr_fg + 200);

    auto grr_y_xyy_xyyy = pbuffer.data(idx_gr_fg + 201);

    auto grr_y_xyy_xyyz = pbuffer.data(idx_gr_fg + 202);

    auto grr_y_xyy_xyzz = pbuffer.data(idx_gr_fg + 203);

    auto grr_y_xyy_xzzz = pbuffer.data(idx_gr_fg + 204);

    auto grr_y_xyy_yyyy = pbuffer.data(idx_gr_fg + 205);

    auto grr_y_xyy_yyyz = pbuffer.data(idx_gr_fg + 206);

    auto grr_y_xyy_yyzz = pbuffer.data(idx_gr_fg + 207);

    auto grr_y_xyy_yzzz = pbuffer.data(idx_gr_fg + 208);

    auto grr_y_xyy_zzzz = pbuffer.data(idx_gr_fg + 209);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxzz, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyzz, gr_xy_xzzz, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyzz, gr_xy_yzzz, gr_xy_zzzz, gr_xyy_xxx, gr_xyy_xxxx, gr_xyy_xxxy, gr_xyy_xxxz, gr_xyy_xxy, gr_xyy_xxyy, gr_xyy_xxyz, gr_xyy_xxz, gr_xyy_xxzz, gr_xyy_xyy, gr_xyy_xyyy, gr_xyy_xyyz, gr_xyy_xyz, gr_xyy_xyzz, gr_xyy_xzz, gr_xyy_xzzz, gr_xyy_yyy, gr_xyy_yyyy, gr_xyy_yyyz, gr_xyy_yyz, gr_xyy_yyzz, gr_xyy_yzz, gr_xyy_yzzz, gr_xyy_zzz, gr_xyy_zzzz, grr_y_xyy_xxxx, grr_y_xyy_xxxy, grr_y_xyy_xxxz, grr_y_xyy_xxyy, grr_y_xyy_xxyz, grr_y_xyy_xxzz, grr_y_xyy_xyyy, grr_y_xyy_xyyz, grr_y_xyy_xyzz, grr_y_xyy_xzzz, grr_y_xyy_yyyy, grr_y_xyy_yyyz, grr_y_xyy_yyzz, grr_y_xyy_yzzz, grr_y_xyy_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyy_xxxx[i] = 2.0 * ts_xy_xxxx[i] * gfe2_0 + 2.0 * gr_xy_xxxx[i] * gfe_0 + ts_xyy_xxxx[i] * gfe_0 * gc_y[i] + gr_xyy_xxxx[i] * gc_y[i];

        grr_y_xyy_xxxy[i] = 2.0 * ts_xy_xxxy[i] * gfe2_0 + 2.0 * gr_xy_xxxy[i] * gfe_0 + ts_xyy_xxx[i] * gfe2_0 + gr_xyy_xxx[i] * gfe_0 + ts_xyy_xxxy[i] * gfe_0 * gc_y[i] + gr_xyy_xxxy[i] * gc_y[i];

        grr_y_xyy_xxxz[i] = 2.0 * ts_xy_xxxz[i] * gfe2_0 + 2.0 * gr_xy_xxxz[i] * gfe_0 + ts_xyy_xxxz[i] * gfe_0 * gc_y[i] + gr_xyy_xxxz[i] * gc_y[i];

        grr_y_xyy_xxyy[i] = 2.0 * ts_xy_xxyy[i] * gfe2_0 + 2.0 * gr_xy_xxyy[i] * gfe_0 + 2.0 * ts_xyy_xxy[i] * gfe2_0 + 2.0 * gr_xyy_xxy[i] * gfe_0 + ts_xyy_xxyy[i] * gfe_0 * gc_y[i] + gr_xyy_xxyy[i] * gc_y[i];

        grr_y_xyy_xxyz[i] = 2.0 * ts_xy_xxyz[i] * gfe2_0 + 2.0 * gr_xy_xxyz[i] * gfe_0 + ts_xyy_xxz[i] * gfe2_0 + gr_xyy_xxz[i] * gfe_0 + ts_xyy_xxyz[i] * gfe_0 * gc_y[i] + gr_xyy_xxyz[i] * gc_y[i];

        grr_y_xyy_xxzz[i] = 2.0 * ts_xy_xxzz[i] * gfe2_0 + 2.0 * gr_xy_xxzz[i] * gfe_0 + ts_xyy_xxzz[i] * gfe_0 * gc_y[i] + gr_xyy_xxzz[i] * gc_y[i];

        grr_y_xyy_xyyy[i] = 2.0 * ts_xy_xyyy[i] * gfe2_0 + 2.0 * gr_xy_xyyy[i] * gfe_0 + 3.0 * ts_xyy_xyy[i] * gfe2_0 + 3.0 * gr_xyy_xyy[i] * gfe_0 + ts_xyy_xyyy[i] * gfe_0 * gc_y[i] + gr_xyy_xyyy[i] * gc_y[i];

        grr_y_xyy_xyyz[i] = 2.0 * ts_xy_xyyz[i] * gfe2_0 + 2.0 * gr_xy_xyyz[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe2_0 + 2.0 * gr_xyy_xyz[i] * gfe_0 + ts_xyy_xyyz[i] * gfe_0 * gc_y[i] + gr_xyy_xyyz[i] * gc_y[i];

        grr_y_xyy_xyzz[i] = 2.0 * ts_xy_xyzz[i] * gfe2_0 + 2.0 * gr_xy_xyzz[i] * gfe_0 + ts_xyy_xzz[i] * gfe2_0 + gr_xyy_xzz[i] * gfe_0 + ts_xyy_xyzz[i] * gfe_0 * gc_y[i] + gr_xyy_xyzz[i] * gc_y[i];

        grr_y_xyy_xzzz[i] = 2.0 * ts_xy_xzzz[i] * gfe2_0 + 2.0 * gr_xy_xzzz[i] * gfe_0 + ts_xyy_xzzz[i] * gfe_0 * gc_y[i] + gr_xyy_xzzz[i] * gc_y[i];

        grr_y_xyy_yyyy[i] = 2.0 * ts_xy_yyyy[i] * gfe2_0 + 2.0 * gr_xy_yyyy[i] * gfe_0 + 4.0 * ts_xyy_yyy[i] * gfe2_0 + 4.0 * gr_xyy_yyy[i] * gfe_0 + ts_xyy_yyyy[i] * gfe_0 * gc_y[i] + gr_xyy_yyyy[i] * gc_y[i];

        grr_y_xyy_yyyz[i] = 2.0 * ts_xy_yyyz[i] * gfe2_0 + 2.0 * gr_xy_yyyz[i] * gfe_0 + 3.0 * ts_xyy_yyz[i] * gfe2_0 + 3.0 * gr_xyy_yyz[i] * gfe_0 + ts_xyy_yyyz[i] * gfe_0 * gc_y[i] + gr_xyy_yyyz[i] * gc_y[i];

        grr_y_xyy_yyzz[i] = 2.0 * ts_xy_yyzz[i] * gfe2_0 + 2.0 * gr_xy_yyzz[i] * gfe_0 + 2.0 * ts_xyy_yzz[i] * gfe2_0 + 2.0 * gr_xyy_yzz[i] * gfe_0 + ts_xyy_yyzz[i] * gfe_0 * gc_y[i] + gr_xyy_yyzz[i] * gc_y[i];

        grr_y_xyy_yzzz[i] = 2.0 * ts_xy_yzzz[i] * gfe2_0 + 2.0 * gr_xy_yzzz[i] * gfe_0 + ts_xyy_zzz[i] * gfe2_0 + gr_xyy_zzz[i] * gfe_0 + ts_xyy_yzzz[i] * gfe_0 * gc_y[i] + gr_xyy_yzzz[i] * gc_y[i];

        grr_y_xyy_zzzz[i] = 2.0 * ts_xy_zzzz[i] * gfe2_0 + 2.0 * gr_xy_zzzz[i] * gfe_0 + ts_xyy_zzzz[i] * gfe_0 * gc_y[i] + gr_xyy_zzzz[i] * gc_y[i];
    }

    // Set up 210-225 components of targeted buffer : FG

    auto grr_y_xyz_xxxx = pbuffer.data(idx_gr_fg + 210);

    auto grr_y_xyz_xxxy = pbuffer.data(idx_gr_fg + 211);

    auto grr_y_xyz_xxxz = pbuffer.data(idx_gr_fg + 212);

    auto grr_y_xyz_xxyy = pbuffer.data(idx_gr_fg + 213);

    auto grr_y_xyz_xxyz = pbuffer.data(idx_gr_fg + 214);

    auto grr_y_xyz_xxzz = pbuffer.data(idx_gr_fg + 215);

    auto grr_y_xyz_xyyy = pbuffer.data(idx_gr_fg + 216);

    auto grr_y_xyz_xyyz = pbuffer.data(idx_gr_fg + 217);

    auto grr_y_xyz_xyzz = pbuffer.data(idx_gr_fg + 218);

    auto grr_y_xyz_xzzz = pbuffer.data(idx_gr_fg + 219);

    auto grr_y_xyz_yyyy = pbuffer.data(idx_gr_fg + 220);

    auto grr_y_xyz_yyyz = pbuffer.data(idx_gr_fg + 221);

    auto grr_y_xyz_yyzz = pbuffer.data(idx_gr_fg + 222);

    auto grr_y_xyz_yzzz = pbuffer.data(idx_gr_fg + 223);

    auto grr_y_xyz_zzzz = pbuffer.data(idx_gr_fg + 224);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xxx, gr_xyz_xxxx, gr_xyz_xxxy, gr_xyz_xxxz, gr_xyz_xxy, gr_xyz_xxyy, gr_xyz_xxyz, gr_xyz_xxz, gr_xyz_xxzz, gr_xyz_xyy, gr_xyz_xyyy, gr_xyz_xyyz, gr_xyz_xyz, gr_xyz_xyzz, gr_xyz_xzz, gr_xyz_xzzz, gr_xyz_yyy, gr_xyz_yyyy, gr_xyz_yyyz, gr_xyz_yyz, gr_xyz_yyzz, gr_xyz_yzz, gr_xyz_yzzz, gr_xyz_zzz, gr_xyz_zzzz, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxzz, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyzz, gr_xz_xzzz, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyzz, gr_xz_yzzz, gr_xz_zzzz, grr_y_xyz_xxxx, grr_y_xyz_xxxy, grr_y_xyz_xxxz, grr_y_xyz_xxyy, grr_y_xyz_xxyz, grr_y_xyz_xxzz, grr_y_xyz_xyyy, grr_y_xyz_xyyz, grr_y_xyz_xyzz, grr_y_xyz_xzzz, grr_y_xyz_yyyy, grr_y_xyz_yyyz, grr_y_xyz_yyzz, grr_y_xyz_yzzz, grr_y_xyz_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyz_xxxx[i] = ts_xz_xxxx[i] * gfe2_0 + gr_xz_xxxx[i] * gfe_0 + ts_xyz_xxxx[i] * gfe_0 * gc_y[i] + gr_xyz_xxxx[i] * gc_y[i];

        grr_y_xyz_xxxy[i] = ts_xz_xxxy[i] * gfe2_0 + gr_xz_xxxy[i] * gfe_0 + ts_xyz_xxx[i] * gfe2_0 + gr_xyz_xxx[i] * gfe_0 + ts_xyz_xxxy[i] * gfe_0 * gc_y[i] + gr_xyz_xxxy[i] * gc_y[i];

        grr_y_xyz_xxxz[i] = ts_xz_xxxz[i] * gfe2_0 + gr_xz_xxxz[i] * gfe_0 + ts_xyz_xxxz[i] * gfe_0 * gc_y[i] + gr_xyz_xxxz[i] * gc_y[i];

        grr_y_xyz_xxyy[i] = ts_xz_xxyy[i] * gfe2_0 + gr_xz_xxyy[i] * gfe_0 + 2.0 * ts_xyz_xxy[i] * gfe2_0 + 2.0 * gr_xyz_xxy[i] * gfe_0 + ts_xyz_xxyy[i] * gfe_0 * gc_y[i] + gr_xyz_xxyy[i] * gc_y[i];

        grr_y_xyz_xxyz[i] = ts_xz_xxyz[i] * gfe2_0 + gr_xz_xxyz[i] * gfe_0 + ts_xyz_xxz[i] * gfe2_0 + gr_xyz_xxz[i] * gfe_0 + ts_xyz_xxyz[i] * gfe_0 * gc_y[i] + gr_xyz_xxyz[i] * gc_y[i];

        grr_y_xyz_xxzz[i] = ts_xz_xxzz[i] * gfe2_0 + gr_xz_xxzz[i] * gfe_0 + ts_xyz_xxzz[i] * gfe_0 * gc_y[i] + gr_xyz_xxzz[i] * gc_y[i];

        grr_y_xyz_xyyy[i] = ts_xz_xyyy[i] * gfe2_0 + gr_xz_xyyy[i] * gfe_0 + 3.0 * ts_xyz_xyy[i] * gfe2_0 + 3.0 * gr_xyz_xyy[i] * gfe_0 + ts_xyz_xyyy[i] * gfe_0 * gc_y[i] + gr_xyz_xyyy[i] * gc_y[i];

        grr_y_xyz_xyyz[i] = ts_xz_xyyz[i] * gfe2_0 + gr_xz_xyyz[i] * gfe_0 + 2.0 * ts_xyz_xyz[i] * gfe2_0 + 2.0 * gr_xyz_xyz[i] * gfe_0 + ts_xyz_xyyz[i] * gfe_0 * gc_y[i] + gr_xyz_xyyz[i] * gc_y[i];

        grr_y_xyz_xyzz[i] = ts_xz_xyzz[i] * gfe2_0 + gr_xz_xyzz[i] * gfe_0 + ts_xyz_xzz[i] * gfe2_0 + gr_xyz_xzz[i] * gfe_0 + ts_xyz_xyzz[i] * gfe_0 * gc_y[i] + gr_xyz_xyzz[i] * gc_y[i];

        grr_y_xyz_xzzz[i] = ts_xz_xzzz[i] * gfe2_0 + gr_xz_xzzz[i] * gfe_0 + ts_xyz_xzzz[i] * gfe_0 * gc_y[i] + gr_xyz_xzzz[i] * gc_y[i];

        grr_y_xyz_yyyy[i] = ts_xz_yyyy[i] * gfe2_0 + gr_xz_yyyy[i] * gfe_0 + 4.0 * ts_xyz_yyy[i] * gfe2_0 + 4.0 * gr_xyz_yyy[i] * gfe_0 + ts_xyz_yyyy[i] * gfe_0 * gc_y[i] + gr_xyz_yyyy[i] * gc_y[i];

        grr_y_xyz_yyyz[i] = ts_xz_yyyz[i] * gfe2_0 + gr_xz_yyyz[i] * gfe_0 + 3.0 * ts_xyz_yyz[i] * gfe2_0 + 3.0 * gr_xyz_yyz[i] * gfe_0 + ts_xyz_yyyz[i] * gfe_0 * gc_y[i] + gr_xyz_yyyz[i] * gc_y[i];

        grr_y_xyz_yyzz[i] = ts_xz_yyzz[i] * gfe2_0 + gr_xz_yyzz[i] * gfe_0 + 2.0 * ts_xyz_yzz[i] * gfe2_0 + 2.0 * gr_xyz_yzz[i] * gfe_0 + ts_xyz_yyzz[i] * gfe_0 * gc_y[i] + gr_xyz_yyzz[i] * gc_y[i];

        grr_y_xyz_yzzz[i] = ts_xz_yzzz[i] * gfe2_0 + gr_xz_yzzz[i] * gfe_0 + ts_xyz_zzz[i] * gfe2_0 + gr_xyz_zzz[i] * gfe_0 + ts_xyz_yzzz[i] * gfe_0 * gc_y[i] + gr_xyz_yzzz[i] * gc_y[i];

        grr_y_xyz_zzzz[i] = ts_xz_zzzz[i] * gfe2_0 + gr_xz_zzzz[i] * gfe_0 + ts_xyz_zzzz[i] * gfe_0 * gc_y[i] + gr_xyz_zzzz[i] * gc_y[i];
    }

    // Set up 225-240 components of targeted buffer : FG

    auto grr_y_xzz_xxxx = pbuffer.data(idx_gr_fg + 225);

    auto grr_y_xzz_xxxy = pbuffer.data(idx_gr_fg + 226);

    auto grr_y_xzz_xxxz = pbuffer.data(idx_gr_fg + 227);

    auto grr_y_xzz_xxyy = pbuffer.data(idx_gr_fg + 228);

    auto grr_y_xzz_xxyz = pbuffer.data(idx_gr_fg + 229);

    auto grr_y_xzz_xxzz = pbuffer.data(idx_gr_fg + 230);

    auto grr_y_xzz_xyyy = pbuffer.data(idx_gr_fg + 231);

    auto grr_y_xzz_xyyz = pbuffer.data(idx_gr_fg + 232);

    auto grr_y_xzz_xyzz = pbuffer.data(idx_gr_fg + 233);

    auto grr_y_xzz_xzzz = pbuffer.data(idx_gr_fg + 234);

    auto grr_y_xzz_yyyy = pbuffer.data(idx_gr_fg + 235);

    auto grr_y_xzz_yyyz = pbuffer.data(idx_gr_fg + 236);

    auto grr_y_xzz_yyzz = pbuffer.data(idx_gr_fg + 237);

    auto grr_y_xzz_yzzz = pbuffer.data(idx_gr_fg + 238);

    auto grr_y_xzz_zzzz = pbuffer.data(idx_gr_fg + 239);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xxx, gr_xzz_xxxx, gr_xzz_xxxy, gr_xzz_xxxz, gr_xzz_xxy, gr_xzz_xxyy, gr_xzz_xxyz, gr_xzz_xxz, gr_xzz_xxzz, gr_xzz_xyy, gr_xzz_xyyy, gr_xzz_xyyz, gr_xzz_xyz, gr_xzz_xyzz, gr_xzz_xzz, gr_xzz_xzzz, gr_xzz_yyy, gr_xzz_yyyy, gr_xzz_yyyz, gr_xzz_yyz, gr_xzz_yyzz, gr_xzz_yzz, gr_xzz_yzzz, gr_xzz_zzz, gr_xzz_zzzz, grr_y_xzz_xxxx, grr_y_xzz_xxxy, grr_y_xzz_xxxz, grr_y_xzz_xxyy, grr_y_xzz_xxyz, grr_y_xzz_xxzz, grr_y_xzz_xyyy, grr_y_xzz_xyyz, grr_y_xzz_xyzz, grr_y_xzz_xzzz, grr_y_xzz_yyyy, grr_y_xzz_yyyz, grr_y_xzz_yyzz, grr_y_xzz_yzzz, grr_y_xzz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xzz_xxxx[i] = ts_xzz_xxxx[i] * gfe_0 * gc_y[i] + gr_xzz_xxxx[i] * gc_y[i];

        grr_y_xzz_xxxy[i] = ts_xzz_xxx[i] * gfe2_0 + gr_xzz_xxx[i] * gfe_0 + ts_xzz_xxxy[i] * gfe_0 * gc_y[i] + gr_xzz_xxxy[i] * gc_y[i];

        grr_y_xzz_xxxz[i] = ts_xzz_xxxz[i] * gfe_0 * gc_y[i] + gr_xzz_xxxz[i] * gc_y[i];

        grr_y_xzz_xxyy[i] = 2.0 * ts_xzz_xxy[i] * gfe2_0 + 2.0 * gr_xzz_xxy[i] * gfe_0 + ts_xzz_xxyy[i] * gfe_0 * gc_y[i] + gr_xzz_xxyy[i] * gc_y[i];

        grr_y_xzz_xxyz[i] = ts_xzz_xxz[i] * gfe2_0 + gr_xzz_xxz[i] * gfe_0 + ts_xzz_xxyz[i] * gfe_0 * gc_y[i] + gr_xzz_xxyz[i] * gc_y[i];

        grr_y_xzz_xxzz[i] = ts_xzz_xxzz[i] * gfe_0 * gc_y[i] + gr_xzz_xxzz[i] * gc_y[i];

        grr_y_xzz_xyyy[i] = 3.0 * ts_xzz_xyy[i] * gfe2_0 + 3.0 * gr_xzz_xyy[i] * gfe_0 + ts_xzz_xyyy[i] * gfe_0 * gc_y[i] + gr_xzz_xyyy[i] * gc_y[i];

        grr_y_xzz_xyyz[i] = 2.0 * ts_xzz_xyz[i] * gfe2_0 + 2.0 * gr_xzz_xyz[i] * gfe_0 + ts_xzz_xyyz[i] * gfe_0 * gc_y[i] + gr_xzz_xyyz[i] * gc_y[i];

        grr_y_xzz_xyzz[i] = ts_xzz_xzz[i] * gfe2_0 + gr_xzz_xzz[i] * gfe_0 + ts_xzz_xyzz[i] * gfe_0 * gc_y[i] + gr_xzz_xyzz[i] * gc_y[i];

        grr_y_xzz_xzzz[i] = ts_xzz_xzzz[i] * gfe_0 * gc_y[i] + gr_xzz_xzzz[i] * gc_y[i];

        grr_y_xzz_yyyy[i] = 4.0 * ts_xzz_yyy[i] * gfe2_0 + 4.0 * gr_xzz_yyy[i] * gfe_0 + ts_xzz_yyyy[i] * gfe_0 * gc_y[i] + gr_xzz_yyyy[i] * gc_y[i];

        grr_y_xzz_yyyz[i] = 3.0 * ts_xzz_yyz[i] * gfe2_0 + 3.0 * gr_xzz_yyz[i] * gfe_0 + ts_xzz_yyyz[i] * gfe_0 * gc_y[i] + gr_xzz_yyyz[i] * gc_y[i];

        grr_y_xzz_yyzz[i] = 2.0 * ts_xzz_yzz[i] * gfe2_0 + 2.0 * gr_xzz_yzz[i] * gfe_0 + ts_xzz_yyzz[i] * gfe_0 * gc_y[i] + gr_xzz_yyzz[i] * gc_y[i];

        grr_y_xzz_yzzz[i] = ts_xzz_zzz[i] * gfe2_0 + gr_xzz_zzz[i] * gfe_0 + ts_xzz_yzzz[i] * gfe_0 * gc_y[i] + gr_xzz_yzzz[i] * gc_y[i];

        grr_y_xzz_zzzz[i] = ts_xzz_zzzz[i] * gfe_0 * gc_y[i] + gr_xzz_zzzz[i] * gc_y[i];
    }

    // Set up 240-255 components of targeted buffer : FG

    auto grr_y_yyy_xxxx = pbuffer.data(idx_gr_fg + 240);

    auto grr_y_yyy_xxxy = pbuffer.data(idx_gr_fg + 241);

    auto grr_y_yyy_xxxz = pbuffer.data(idx_gr_fg + 242);

    auto grr_y_yyy_xxyy = pbuffer.data(idx_gr_fg + 243);

    auto grr_y_yyy_xxyz = pbuffer.data(idx_gr_fg + 244);

    auto grr_y_yyy_xxzz = pbuffer.data(idx_gr_fg + 245);

    auto grr_y_yyy_xyyy = pbuffer.data(idx_gr_fg + 246);

    auto grr_y_yyy_xyyz = pbuffer.data(idx_gr_fg + 247);

    auto grr_y_yyy_xyzz = pbuffer.data(idx_gr_fg + 248);

    auto grr_y_yyy_xzzz = pbuffer.data(idx_gr_fg + 249);

    auto grr_y_yyy_yyyy = pbuffer.data(idx_gr_fg + 250);

    auto grr_y_yyy_yyyz = pbuffer.data(idx_gr_fg + 251);

    auto grr_y_yyy_yyzz = pbuffer.data(idx_gr_fg + 252);

    auto grr_y_yyy_yzzz = pbuffer.data(idx_gr_fg + 253);

    auto grr_y_yyy_zzzz = pbuffer.data(idx_gr_fg + 254);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxzz, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyzz, gr_yy_xzzz, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyzz, gr_yy_yzzz, gr_yy_zzzz, gr_yyy_xxx, gr_yyy_xxxx, gr_yyy_xxxy, gr_yyy_xxxz, gr_yyy_xxy, gr_yyy_xxyy, gr_yyy_xxyz, gr_yyy_xxz, gr_yyy_xxzz, gr_yyy_xyy, gr_yyy_xyyy, gr_yyy_xyyz, gr_yyy_xyz, gr_yyy_xyzz, gr_yyy_xzz, gr_yyy_xzzz, gr_yyy_yyy, gr_yyy_yyyy, gr_yyy_yyyz, gr_yyy_yyz, gr_yyy_yyzz, gr_yyy_yzz, gr_yyy_yzzz, gr_yyy_zzz, gr_yyy_zzzz, grr_y_yyy_xxxx, grr_y_yyy_xxxy, grr_y_yyy_xxxz, grr_y_yyy_xxyy, grr_y_yyy_xxyz, grr_y_yyy_xxzz, grr_y_yyy_xyyy, grr_y_yyy_xyyz, grr_y_yyy_xyzz, grr_y_yyy_xzzz, grr_y_yyy_yyyy, grr_y_yyy_yyyz, grr_y_yyy_yyzz, grr_y_yyy_yzzz, grr_y_yyy_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyy_xxxx[i] = 3.0 * ts_yy_xxxx[i] * gfe2_0 + 3.0 * gr_yy_xxxx[i] * gfe_0 + ts_yyy_xxxx[i] * gfe_0 * gc_y[i] + gr_yyy_xxxx[i] * gc_y[i];

        grr_y_yyy_xxxy[i] = 3.0 * ts_yy_xxxy[i] * gfe2_0 + 3.0 * gr_yy_xxxy[i] * gfe_0 + ts_yyy_xxx[i] * gfe2_0 + gr_yyy_xxx[i] * gfe_0 + ts_yyy_xxxy[i] * gfe_0 * gc_y[i] + gr_yyy_xxxy[i] * gc_y[i];

        grr_y_yyy_xxxz[i] = 3.0 * ts_yy_xxxz[i] * gfe2_0 + 3.0 * gr_yy_xxxz[i] * gfe_0 + ts_yyy_xxxz[i] * gfe_0 * gc_y[i] + gr_yyy_xxxz[i] * gc_y[i];

        grr_y_yyy_xxyy[i] = 3.0 * ts_yy_xxyy[i] * gfe2_0 + 3.0 * gr_yy_xxyy[i] * gfe_0 + 2.0 * ts_yyy_xxy[i] * gfe2_0 + 2.0 * gr_yyy_xxy[i] * gfe_0 + ts_yyy_xxyy[i] * gfe_0 * gc_y[i] + gr_yyy_xxyy[i] * gc_y[i];

        grr_y_yyy_xxyz[i] = 3.0 * ts_yy_xxyz[i] * gfe2_0 + 3.0 * gr_yy_xxyz[i] * gfe_0 + ts_yyy_xxz[i] * gfe2_0 + gr_yyy_xxz[i] * gfe_0 + ts_yyy_xxyz[i] * gfe_0 * gc_y[i] + gr_yyy_xxyz[i] * gc_y[i];

        grr_y_yyy_xxzz[i] = 3.0 * ts_yy_xxzz[i] * gfe2_0 + 3.0 * gr_yy_xxzz[i] * gfe_0 + ts_yyy_xxzz[i] * gfe_0 * gc_y[i] + gr_yyy_xxzz[i] * gc_y[i];

        grr_y_yyy_xyyy[i] = 3.0 * ts_yy_xyyy[i] * gfe2_0 + 3.0 * gr_yy_xyyy[i] * gfe_0 + 3.0 * ts_yyy_xyy[i] * gfe2_0 + 3.0 * gr_yyy_xyy[i] * gfe_0 + ts_yyy_xyyy[i] * gfe_0 * gc_y[i] + gr_yyy_xyyy[i] * gc_y[i];

        grr_y_yyy_xyyz[i] = 3.0 * ts_yy_xyyz[i] * gfe2_0 + 3.0 * gr_yy_xyyz[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe2_0 + 2.0 * gr_yyy_xyz[i] * gfe_0 + ts_yyy_xyyz[i] * gfe_0 * gc_y[i] + gr_yyy_xyyz[i] * gc_y[i];

        grr_y_yyy_xyzz[i] = 3.0 * ts_yy_xyzz[i] * gfe2_0 + 3.0 * gr_yy_xyzz[i] * gfe_0 + ts_yyy_xzz[i] * gfe2_0 + gr_yyy_xzz[i] * gfe_0 + ts_yyy_xyzz[i] * gfe_0 * gc_y[i] + gr_yyy_xyzz[i] * gc_y[i];

        grr_y_yyy_xzzz[i] = 3.0 * ts_yy_xzzz[i] * gfe2_0 + 3.0 * gr_yy_xzzz[i] * gfe_0 + ts_yyy_xzzz[i] * gfe_0 * gc_y[i] + gr_yyy_xzzz[i] * gc_y[i];

        grr_y_yyy_yyyy[i] = 3.0 * ts_yy_yyyy[i] * gfe2_0 + 3.0 * gr_yy_yyyy[i] * gfe_0 + 4.0 * ts_yyy_yyy[i] * gfe2_0 + 4.0 * gr_yyy_yyy[i] * gfe_0 + ts_yyy_yyyy[i] * gfe_0 * gc_y[i] + gr_yyy_yyyy[i] * gc_y[i];

        grr_y_yyy_yyyz[i] = 3.0 * ts_yy_yyyz[i] * gfe2_0 + 3.0 * gr_yy_yyyz[i] * gfe_0 + 3.0 * ts_yyy_yyz[i] * gfe2_0 + 3.0 * gr_yyy_yyz[i] * gfe_0 + ts_yyy_yyyz[i] * gfe_0 * gc_y[i] + gr_yyy_yyyz[i] * gc_y[i];

        grr_y_yyy_yyzz[i] = 3.0 * ts_yy_yyzz[i] * gfe2_0 + 3.0 * gr_yy_yyzz[i] * gfe_0 + 2.0 * ts_yyy_yzz[i] * gfe2_0 + 2.0 * gr_yyy_yzz[i] * gfe_0 + ts_yyy_yyzz[i] * gfe_0 * gc_y[i] + gr_yyy_yyzz[i] * gc_y[i];

        grr_y_yyy_yzzz[i] = 3.0 * ts_yy_yzzz[i] * gfe2_0 + 3.0 * gr_yy_yzzz[i] * gfe_0 + ts_yyy_zzz[i] * gfe2_0 + gr_yyy_zzz[i] * gfe_0 + ts_yyy_yzzz[i] * gfe_0 * gc_y[i] + gr_yyy_yzzz[i] * gc_y[i];

        grr_y_yyy_zzzz[i] = 3.0 * ts_yy_zzzz[i] * gfe2_0 + 3.0 * gr_yy_zzzz[i] * gfe_0 + ts_yyy_zzzz[i] * gfe_0 * gc_y[i] + gr_yyy_zzzz[i] * gc_y[i];
    }

    // Set up 255-270 components of targeted buffer : FG

    auto grr_y_yyz_xxxx = pbuffer.data(idx_gr_fg + 255);

    auto grr_y_yyz_xxxy = pbuffer.data(idx_gr_fg + 256);

    auto grr_y_yyz_xxxz = pbuffer.data(idx_gr_fg + 257);

    auto grr_y_yyz_xxyy = pbuffer.data(idx_gr_fg + 258);

    auto grr_y_yyz_xxyz = pbuffer.data(idx_gr_fg + 259);

    auto grr_y_yyz_xxzz = pbuffer.data(idx_gr_fg + 260);

    auto grr_y_yyz_xyyy = pbuffer.data(idx_gr_fg + 261);

    auto grr_y_yyz_xyyz = pbuffer.data(idx_gr_fg + 262);

    auto grr_y_yyz_xyzz = pbuffer.data(idx_gr_fg + 263);

    auto grr_y_yyz_xzzz = pbuffer.data(idx_gr_fg + 264);

    auto grr_y_yyz_yyyy = pbuffer.data(idx_gr_fg + 265);

    auto grr_y_yyz_yyyz = pbuffer.data(idx_gr_fg + 266);

    auto grr_y_yyz_yyzz = pbuffer.data(idx_gr_fg + 267);

    auto grr_y_yyz_yzzz = pbuffer.data(idx_gr_fg + 268);

    auto grr_y_yyz_zzzz = pbuffer.data(idx_gr_fg + 269);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xxx, gr_yyz_xxxx, gr_yyz_xxxy, gr_yyz_xxxz, gr_yyz_xxy, gr_yyz_xxyy, gr_yyz_xxyz, gr_yyz_xxz, gr_yyz_xxzz, gr_yyz_xyy, gr_yyz_xyyy, gr_yyz_xyyz, gr_yyz_xyz, gr_yyz_xyzz, gr_yyz_xzz, gr_yyz_xzzz, gr_yyz_yyy, gr_yyz_yyyy, gr_yyz_yyyz, gr_yyz_yyz, gr_yyz_yyzz, gr_yyz_yzz, gr_yyz_yzzz, gr_yyz_zzz, gr_yyz_zzzz, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxzz, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyzz, gr_yz_xzzz, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyzz, gr_yz_yzzz, gr_yz_zzzz, grr_y_yyz_xxxx, grr_y_yyz_xxxy, grr_y_yyz_xxxz, grr_y_yyz_xxyy, grr_y_yyz_xxyz, grr_y_yyz_xxzz, grr_y_yyz_xyyy, grr_y_yyz_xyyz, grr_y_yyz_xyzz, grr_y_yyz_xzzz, grr_y_yyz_yyyy, grr_y_yyz_yyyz, grr_y_yyz_yyzz, grr_y_yyz_yzzz, grr_y_yyz_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyz_xxxx[i] = 2.0 * ts_yz_xxxx[i] * gfe2_0 + 2.0 * gr_yz_xxxx[i] * gfe_0 + ts_yyz_xxxx[i] * gfe_0 * gc_y[i] + gr_yyz_xxxx[i] * gc_y[i];

        grr_y_yyz_xxxy[i] = 2.0 * ts_yz_xxxy[i] * gfe2_0 + 2.0 * gr_yz_xxxy[i] * gfe_0 + ts_yyz_xxx[i] * gfe2_0 + gr_yyz_xxx[i] * gfe_0 + ts_yyz_xxxy[i] * gfe_0 * gc_y[i] + gr_yyz_xxxy[i] * gc_y[i];

        grr_y_yyz_xxxz[i] = 2.0 * ts_yz_xxxz[i] * gfe2_0 + 2.0 * gr_yz_xxxz[i] * gfe_0 + ts_yyz_xxxz[i] * gfe_0 * gc_y[i] + gr_yyz_xxxz[i] * gc_y[i];

        grr_y_yyz_xxyy[i] = 2.0 * ts_yz_xxyy[i] * gfe2_0 + 2.0 * gr_yz_xxyy[i] * gfe_0 + 2.0 * ts_yyz_xxy[i] * gfe2_0 + 2.0 * gr_yyz_xxy[i] * gfe_0 + ts_yyz_xxyy[i] * gfe_0 * gc_y[i] + gr_yyz_xxyy[i] * gc_y[i];

        grr_y_yyz_xxyz[i] = 2.0 * ts_yz_xxyz[i] * gfe2_0 + 2.0 * gr_yz_xxyz[i] * gfe_0 + ts_yyz_xxz[i] * gfe2_0 + gr_yyz_xxz[i] * gfe_0 + ts_yyz_xxyz[i] * gfe_0 * gc_y[i] + gr_yyz_xxyz[i] * gc_y[i];

        grr_y_yyz_xxzz[i] = 2.0 * ts_yz_xxzz[i] * gfe2_0 + 2.0 * gr_yz_xxzz[i] * gfe_0 + ts_yyz_xxzz[i] * gfe_0 * gc_y[i] + gr_yyz_xxzz[i] * gc_y[i];

        grr_y_yyz_xyyy[i] = 2.0 * ts_yz_xyyy[i] * gfe2_0 + 2.0 * gr_yz_xyyy[i] * gfe_0 + 3.0 * ts_yyz_xyy[i] * gfe2_0 + 3.0 * gr_yyz_xyy[i] * gfe_0 + ts_yyz_xyyy[i] * gfe_0 * gc_y[i] + gr_yyz_xyyy[i] * gc_y[i];

        grr_y_yyz_xyyz[i] = 2.0 * ts_yz_xyyz[i] * gfe2_0 + 2.0 * gr_yz_xyyz[i] * gfe_0 + 2.0 * ts_yyz_xyz[i] * gfe2_0 + 2.0 * gr_yyz_xyz[i] * gfe_0 + ts_yyz_xyyz[i] * gfe_0 * gc_y[i] + gr_yyz_xyyz[i] * gc_y[i];

        grr_y_yyz_xyzz[i] = 2.0 * ts_yz_xyzz[i] * gfe2_0 + 2.0 * gr_yz_xyzz[i] * gfe_0 + ts_yyz_xzz[i] * gfe2_0 + gr_yyz_xzz[i] * gfe_0 + ts_yyz_xyzz[i] * gfe_0 * gc_y[i] + gr_yyz_xyzz[i] * gc_y[i];

        grr_y_yyz_xzzz[i] = 2.0 * ts_yz_xzzz[i] * gfe2_0 + 2.0 * gr_yz_xzzz[i] * gfe_0 + ts_yyz_xzzz[i] * gfe_0 * gc_y[i] + gr_yyz_xzzz[i] * gc_y[i];

        grr_y_yyz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe2_0 + 2.0 * gr_yz_yyyy[i] * gfe_0 + 4.0 * ts_yyz_yyy[i] * gfe2_0 + 4.0 * gr_yyz_yyy[i] * gfe_0 + ts_yyz_yyyy[i] * gfe_0 * gc_y[i] + gr_yyz_yyyy[i] * gc_y[i];

        grr_y_yyz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe2_0 + 2.0 * gr_yz_yyyz[i] * gfe_0 + 3.0 * ts_yyz_yyz[i] * gfe2_0 + 3.0 * gr_yyz_yyz[i] * gfe_0 + ts_yyz_yyyz[i] * gfe_0 * gc_y[i] + gr_yyz_yyyz[i] * gc_y[i];

        grr_y_yyz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe2_0 + 2.0 * gr_yz_yyzz[i] * gfe_0 + 2.0 * ts_yyz_yzz[i] * gfe2_0 + 2.0 * gr_yyz_yzz[i] * gfe_0 + ts_yyz_yyzz[i] * gfe_0 * gc_y[i] + gr_yyz_yyzz[i] * gc_y[i];

        grr_y_yyz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe2_0 + 2.0 * gr_yz_yzzz[i] * gfe_0 + ts_yyz_zzz[i] * gfe2_0 + gr_yyz_zzz[i] * gfe_0 + ts_yyz_yzzz[i] * gfe_0 * gc_y[i] + gr_yyz_yzzz[i] * gc_y[i];

        grr_y_yyz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe2_0 + 2.0 * gr_yz_zzzz[i] * gfe_0 + ts_yyz_zzzz[i] * gfe_0 * gc_y[i] + gr_yyz_zzzz[i] * gc_y[i];
    }

    // Set up 270-285 components of targeted buffer : FG

    auto grr_y_yzz_xxxx = pbuffer.data(idx_gr_fg + 270);

    auto grr_y_yzz_xxxy = pbuffer.data(idx_gr_fg + 271);

    auto grr_y_yzz_xxxz = pbuffer.data(idx_gr_fg + 272);

    auto grr_y_yzz_xxyy = pbuffer.data(idx_gr_fg + 273);

    auto grr_y_yzz_xxyz = pbuffer.data(idx_gr_fg + 274);

    auto grr_y_yzz_xxzz = pbuffer.data(idx_gr_fg + 275);

    auto grr_y_yzz_xyyy = pbuffer.data(idx_gr_fg + 276);

    auto grr_y_yzz_xyyz = pbuffer.data(idx_gr_fg + 277);

    auto grr_y_yzz_xyzz = pbuffer.data(idx_gr_fg + 278);

    auto grr_y_yzz_xzzz = pbuffer.data(idx_gr_fg + 279);

    auto grr_y_yzz_yyyy = pbuffer.data(idx_gr_fg + 280);

    auto grr_y_yzz_yyyz = pbuffer.data(idx_gr_fg + 281);

    auto grr_y_yzz_yyzz = pbuffer.data(idx_gr_fg + 282);

    auto grr_y_yzz_yzzz = pbuffer.data(idx_gr_fg + 283);

    auto grr_y_yzz_zzzz = pbuffer.data(idx_gr_fg + 284);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xxx, gr_yzz_xxxx, gr_yzz_xxxy, gr_yzz_xxxz, gr_yzz_xxy, gr_yzz_xxyy, gr_yzz_xxyz, gr_yzz_xxz, gr_yzz_xxzz, gr_yzz_xyy, gr_yzz_xyyy, gr_yzz_xyyz, gr_yzz_xyz, gr_yzz_xyzz, gr_yzz_xzz, gr_yzz_xzzz, gr_yzz_yyy, gr_yzz_yyyy, gr_yzz_yyyz, gr_yzz_yyz, gr_yzz_yyzz, gr_yzz_yzz, gr_yzz_yzzz, gr_yzz_zzz, gr_yzz_zzzz, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxzz, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyzz, gr_zz_xzzz, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyzz, gr_zz_yzzz, gr_zz_zzzz, grr_y_yzz_xxxx, grr_y_yzz_xxxy, grr_y_yzz_xxxz, grr_y_yzz_xxyy, grr_y_yzz_xxyz, grr_y_yzz_xxzz, grr_y_yzz_xyyy, grr_y_yzz_xyyz, grr_y_yzz_xyzz, grr_y_yzz_xzzz, grr_y_yzz_yyyy, grr_y_yzz_yyyz, grr_y_yzz_yyzz, grr_y_yzz_yzzz, grr_y_yzz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yzz_xxxx[i] = ts_zz_xxxx[i] * gfe2_0 + gr_zz_xxxx[i] * gfe_0 + ts_yzz_xxxx[i] * gfe_0 * gc_y[i] + gr_yzz_xxxx[i] * gc_y[i];

        grr_y_yzz_xxxy[i] = ts_zz_xxxy[i] * gfe2_0 + gr_zz_xxxy[i] * gfe_0 + ts_yzz_xxx[i] * gfe2_0 + gr_yzz_xxx[i] * gfe_0 + ts_yzz_xxxy[i] * gfe_0 * gc_y[i] + gr_yzz_xxxy[i] * gc_y[i];

        grr_y_yzz_xxxz[i] = ts_zz_xxxz[i] * gfe2_0 + gr_zz_xxxz[i] * gfe_0 + ts_yzz_xxxz[i] * gfe_0 * gc_y[i] + gr_yzz_xxxz[i] * gc_y[i];

        grr_y_yzz_xxyy[i] = ts_zz_xxyy[i] * gfe2_0 + gr_zz_xxyy[i] * gfe_0 + 2.0 * ts_yzz_xxy[i] * gfe2_0 + 2.0 * gr_yzz_xxy[i] * gfe_0 + ts_yzz_xxyy[i] * gfe_0 * gc_y[i] + gr_yzz_xxyy[i] * gc_y[i];

        grr_y_yzz_xxyz[i] = ts_zz_xxyz[i] * gfe2_0 + gr_zz_xxyz[i] * gfe_0 + ts_yzz_xxz[i] * gfe2_0 + gr_yzz_xxz[i] * gfe_0 + ts_yzz_xxyz[i] * gfe_0 * gc_y[i] + gr_yzz_xxyz[i] * gc_y[i];

        grr_y_yzz_xxzz[i] = ts_zz_xxzz[i] * gfe2_0 + gr_zz_xxzz[i] * gfe_0 + ts_yzz_xxzz[i] * gfe_0 * gc_y[i] + gr_yzz_xxzz[i] * gc_y[i];

        grr_y_yzz_xyyy[i] = ts_zz_xyyy[i] * gfe2_0 + gr_zz_xyyy[i] * gfe_0 + 3.0 * ts_yzz_xyy[i] * gfe2_0 + 3.0 * gr_yzz_xyy[i] * gfe_0 + ts_yzz_xyyy[i] * gfe_0 * gc_y[i] + gr_yzz_xyyy[i] * gc_y[i];

        grr_y_yzz_xyyz[i] = ts_zz_xyyz[i] * gfe2_0 + gr_zz_xyyz[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe2_0 + 2.0 * gr_yzz_xyz[i] * gfe_0 + ts_yzz_xyyz[i] * gfe_0 * gc_y[i] + gr_yzz_xyyz[i] * gc_y[i];

        grr_y_yzz_xyzz[i] = ts_zz_xyzz[i] * gfe2_0 + gr_zz_xyzz[i] * gfe_0 + ts_yzz_xzz[i] * gfe2_0 + gr_yzz_xzz[i] * gfe_0 + ts_yzz_xyzz[i] * gfe_0 * gc_y[i] + gr_yzz_xyzz[i] * gc_y[i];

        grr_y_yzz_xzzz[i] = ts_zz_xzzz[i] * gfe2_0 + gr_zz_xzzz[i] * gfe_0 + ts_yzz_xzzz[i] * gfe_0 * gc_y[i] + gr_yzz_xzzz[i] * gc_y[i];

        grr_y_yzz_yyyy[i] = ts_zz_yyyy[i] * gfe2_0 + gr_zz_yyyy[i] * gfe_0 + 4.0 * ts_yzz_yyy[i] * gfe2_0 + 4.0 * gr_yzz_yyy[i] * gfe_0 + ts_yzz_yyyy[i] * gfe_0 * gc_y[i] + gr_yzz_yyyy[i] * gc_y[i];

        grr_y_yzz_yyyz[i] = ts_zz_yyyz[i] * gfe2_0 + gr_zz_yyyz[i] * gfe_0 + 3.0 * ts_yzz_yyz[i] * gfe2_0 + 3.0 * gr_yzz_yyz[i] * gfe_0 + ts_yzz_yyyz[i] * gfe_0 * gc_y[i] + gr_yzz_yyyz[i] * gc_y[i];

        grr_y_yzz_yyzz[i] = ts_zz_yyzz[i] * gfe2_0 + gr_zz_yyzz[i] * gfe_0 + 2.0 * ts_yzz_yzz[i] * gfe2_0 + 2.0 * gr_yzz_yzz[i] * gfe_0 + ts_yzz_yyzz[i] * gfe_0 * gc_y[i] + gr_yzz_yyzz[i] * gc_y[i];

        grr_y_yzz_yzzz[i] = ts_zz_yzzz[i] * gfe2_0 + gr_zz_yzzz[i] * gfe_0 + ts_yzz_zzz[i] * gfe2_0 + gr_yzz_zzz[i] * gfe_0 + ts_yzz_yzzz[i] * gfe_0 * gc_y[i] + gr_yzz_yzzz[i] * gc_y[i];

        grr_y_yzz_zzzz[i] = ts_zz_zzzz[i] * gfe2_0 + gr_zz_zzzz[i] * gfe_0 + ts_yzz_zzzz[i] * gfe_0 * gc_y[i] + gr_yzz_zzzz[i] * gc_y[i];
    }

    // Set up 285-300 components of targeted buffer : FG

    auto grr_y_zzz_xxxx = pbuffer.data(idx_gr_fg + 285);

    auto grr_y_zzz_xxxy = pbuffer.data(idx_gr_fg + 286);

    auto grr_y_zzz_xxxz = pbuffer.data(idx_gr_fg + 287);

    auto grr_y_zzz_xxyy = pbuffer.data(idx_gr_fg + 288);

    auto grr_y_zzz_xxyz = pbuffer.data(idx_gr_fg + 289);

    auto grr_y_zzz_xxzz = pbuffer.data(idx_gr_fg + 290);

    auto grr_y_zzz_xyyy = pbuffer.data(idx_gr_fg + 291);

    auto grr_y_zzz_xyyz = pbuffer.data(idx_gr_fg + 292);

    auto grr_y_zzz_xyzz = pbuffer.data(idx_gr_fg + 293);

    auto grr_y_zzz_xzzz = pbuffer.data(idx_gr_fg + 294);

    auto grr_y_zzz_yyyy = pbuffer.data(idx_gr_fg + 295);

    auto grr_y_zzz_yyyz = pbuffer.data(idx_gr_fg + 296);

    auto grr_y_zzz_yyzz = pbuffer.data(idx_gr_fg + 297);

    auto grr_y_zzz_yzzz = pbuffer.data(idx_gr_fg + 298);

    auto grr_y_zzz_zzzz = pbuffer.data(idx_gr_fg + 299);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xxx, gr_zzz_xxxx, gr_zzz_xxxy, gr_zzz_xxxz, gr_zzz_xxy, gr_zzz_xxyy, gr_zzz_xxyz, gr_zzz_xxz, gr_zzz_xxzz, gr_zzz_xyy, gr_zzz_xyyy, gr_zzz_xyyz, gr_zzz_xyz, gr_zzz_xyzz, gr_zzz_xzz, gr_zzz_xzzz, gr_zzz_yyy, gr_zzz_yyyy, gr_zzz_yyyz, gr_zzz_yyz, gr_zzz_yyzz, gr_zzz_yzz, gr_zzz_yzzz, gr_zzz_zzz, gr_zzz_zzzz, grr_y_zzz_xxxx, grr_y_zzz_xxxy, grr_y_zzz_xxxz, grr_y_zzz_xxyy, grr_y_zzz_xxyz, grr_y_zzz_xxzz, grr_y_zzz_xyyy, grr_y_zzz_xyyz, grr_y_zzz_xyzz, grr_y_zzz_xzzz, grr_y_zzz_yyyy, grr_y_zzz_yyyz, grr_y_zzz_yyzz, grr_y_zzz_yzzz, grr_y_zzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zzz_xxxx[i] = ts_zzz_xxxx[i] * gfe_0 * gc_y[i] + gr_zzz_xxxx[i] * gc_y[i];

        grr_y_zzz_xxxy[i] = ts_zzz_xxx[i] * gfe2_0 + gr_zzz_xxx[i] * gfe_0 + ts_zzz_xxxy[i] * gfe_0 * gc_y[i] + gr_zzz_xxxy[i] * gc_y[i];

        grr_y_zzz_xxxz[i] = ts_zzz_xxxz[i] * gfe_0 * gc_y[i] + gr_zzz_xxxz[i] * gc_y[i];

        grr_y_zzz_xxyy[i] = 2.0 * ts_zzz_xxy[i] * gfe2_0 + 2.0 * gr_zzz_xxy[i] * gfe_0 + ts_zzz_xxyy[i] * gfe_0 * gc_y[i] + gr_zzz_xxyy[i] * gc_y[i];

        grr_y_zzz_xxyz[i] = ts_zzz_xxz[i] * gfe2_0 + gr_zzz_xxz[i] * gfe_0 + ts_zzz_xxyz[i] * gfe_0 * gc_y[i] + gr_zzz_xxyz[i] * gc_y[i];

        grr_y_zzz_xxzz[i] = ts_zzz_xxzz[i] * gfe_0 * gc_y[i] + gr_zzz_xxzz[i] * gc_y[i];

        grr_y_zzz_xyyy[i] = 3.0 * ts_zzz_xyy[i] * gfe2_0 + 3.0 * gr_zzz_xyy[i] * gfe_0 + ts_zzz_xyyy[i] * gfe_0 * gc_y[i] + gr_zzz_xyyy[i] * gc_y[i];

        grr_y_zzz_xyyz[i] = 2.0 * ts_zzz_xyz[i] * gfe2_0 + 2.0 * gr_zzz_xyz[i] * gfe_0 + ts_zzz_xyyz[i] * gfe_0 * gc_y[i] + gr_zzz_xyyz[i] * gc_y[i];

        grr_y_zzz_xyzz[i] = ts_zzz_xzz[i] * gfe2_0 + gr_zzz_xzz[i] * gfe_0 + ts_zzz_xyzz[i] * gfe_0 * gc_y[i] + gr_zzz_xyzz[i] * gc_y[i];

        grr_y_zzz_xzzz[i] = ts_zzz_xzzz[i] * gfe_0 * gc_y[i] + gr_zzz_xzzz[i] * gc_y[i];

        grr_y_zzz_yyyy[i] = 4.0 * ts_zzz_yyy[i] * gfe2_0 + 4.0 * gr_zzz_yyy[i] * gfe_0 + ts_zzz_yyyy[i] * gfe_0 * gc_y[i] + gr_zzz_yyyy[i] * gc_y[i];

        grr_y_zzz_yyyz[i] = 3.0 * ts_zzz_yyz[i] * gfe2_0 + 3.0 * gr_zzz_yyz[i] * gfe_0 + ts_zzz_yyyz[i] * gfe_0 * gc_y[i] + gr_zzz_yyyz[i] * gc_y[i];

        grr_y_zzz_yyzz[i] = 2.0 * ts_zzz_yzz[i] * gfe2_0 + 2.0 * gr_zzz_yzz[i] * gfe_0 + ts_zzz_yyzz[i] * gfe_0 * gc_y[i] + gr_zzz_yyzz[i] * gc_y[i];

        grr_y_zzz_yzzz[i] = ts_zzz_zzz[i] * gfe2_0 + gr_zzz_zzz[i] * gfe_0 + ts_zzz_yzzz[i] * gfe_0 * gc_y[i] + gr_zzz_yzzz[i] * gc_y[i];

        grr_y_zzz_zzzz[i] = ts_zzz_zzzz[i] * gfe_0 * gc_y[i] + gr_zzz_zzzz[i] * gc_y[i];
    }

    // Set up 300-315 components of targeted buffer : FG

    auto grr_z_xxx_xxxx = pbuffer.data(idx_gr_fg + 300);

    auto grr_z_xxx_xxxy = pbuffer.data(idx_gr_fg + 301);

    auto grr_z_xxx_xxxz = pbuffer.data(idx_gr_fg + 302);

    auto grr_z_xxx_xxyy = pbuffer.data(idx_gr_fg + 303);

    auto grr_z_xxx_xxyz = pbuffer.data(idx_gr_fg + 304);

    auto grr_z_xxx_xxzz = pbuffer.data(idx_gr_fg + 305);

    auto grr_z_xxx_xyyy = pbuffer.data(idx_gr_fg + 306);

    auto grr_z_xxx_xyyz = pbuffer.data(idx_gr_fg + 307);

    auto grr_z_xxx_xyzz = pbuffer.data(idx_gr_fg + 308);

    auto grr_z_xxx_xzzz = pbuffer.data(idx_gr_fg + 309);

    auto grr_z_xxx_yyyy = pbuffer.data(idx_gr_fg + 310);

    auto grr_z_xxx_yyyz = pbuffer.data(idx_gr_fg + 311);

    auto grr_z_xxx_yyzz = pbuffer.data(idx_gr_fg + 312);

    auto grr_z_xxx_yzzz = pbuffer.data(idx_gr_fg + 313);

    auto grr_z_xxx_zzzz = pbuffer.data(idx_gr_fg + 314);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xxx, gr_xxx_xxxx, gr_xxx_xxxy, gr_xxx_xxxz, gr_xxx_xxy, gr_xxx_xxyy, gr_xxx_xxyz, gr_xxx_xxz, gr_xxx_xxzz, gr_xxx_xyy, gr_xxx_xyyy, gr_xxx_xyyz, gr_xxx_xyz, gr_xxx_xyzz, gr_xxx_xzz, gr_xxx_xzzz, gr_xxx_yyy, gr_xxx_yyyy, gr_xxx_yyyz, gr_xxx_yyz, gr_xxx_yyzz, gr_xxx_yzz, gr_xxx_yzzz, gr_xxx_zzz, gr_xxx_zzzz, grr_z_xxx_xxxx, grr_z_xxx_xxxy, grr_z_xxx_xxxz, grr_z_xxx_xxyy, grr_z_xxx_xxyz, grr_z_xxx_xxzz, grr_z_xxx_xyyy, grr_z_xxx_xyyz, grr_z_xxx_xyzz, grr_z_xxx_xzzz, grr_z_xxx_yyyy, grr_z_xxx_yyyz, grr_z_xxx_yyzz, grr_z_xxx_yzzz, grr_z_xxx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxx_xxxx[i] = ts_xxx_xxxx[i] * gfe_0 * gc_z[i] + gr_xxx_xxxx[i] * gc_z[i];

        grr_z_xxx_xxxy[i] = ts_xxx_xxxy[i] * gfe_0 * gc_z[i] + gr_xxx_xxxy[i] * gc_z[i];

        grr_z_xxx_xxxz[i] = ts_xxx_xxx[i] * gfe2_0 + gr_xxx_xxx[i] * gfe_0 + ts_xxx_xxxz[i] * gfe_0 * gc_z[i] + gr_xxx_xxxz[i] * gc_z[i];

        grr_z_xxx_xxyy[i] = ts_xxx_xxyy[i] * gfe_0 * gc_z[i] + gr_xxx_xxyy[i] * gc_z[i];

        grr_z_xxx_xxyz[i] = ts_xxx_xxy[i] * gfe2_0 + gr_xxx_xxy[i] * gfe_0 + ts_xxx_xxyz[i] * gfe_0 * gc_z[i] + gr_xxx_xxyz[i] * gc_z[i];

        grr_z_xxx_xxzz[i] = 2.0 * ts_xxx_xxz[i] * gfe2_0 + 2.0 * gr_xxx_xxz[i] * gfe_0 + ts_xxx_xxzz[i] * gfe_0 * gc_z[i] + gr_xxx_xxzz[i] * gc_z[i];

        grr_z_xxx_xyyy[i] = ts_xxx_xyyy[i] * gfe_0 * gc_z[i] + gr_xxx_xyyy[i] * gc_z[i];

        grr_z_xxx_xyyz[i] = ts_xxx_xyy[i] * gfe2_0 + gr_xxx_xyy[i] * gfe_0 + ts_xxx_xyyz[i] * gfe_0 * gc_z[i] + gr_xxx_xyyz[i] * gc_z[i];

        grr_z_xxx_xyzz[i] = 2.0 * ts_xxx_xyz[i] * gfe2_0 + 2.0 * gr_xxx_xyz[i] * gfe_0 + ts_xxx_xyzz[i] * gfe_0 * gc_z[i] + gr_xxx_xyzz[i] * gc_z[i];

        grr_z_xxx_xzzz[i] = 3.0 * ts_xxx_xzz[i] * gfe2_0 + 3.0 * gr_xxx_xzz[i] * gfe_0 + ts_xxx_xzzz[i] * gfe_0 * gc_z[i] + gr_xxx_xzzz[i] * gc_z[i];

        grr_z_xxx_yyyy[i] = ts_xxx_yyyy[i] * gfe_0 * gc_z[i] + gr_xxx_yyyy[i] * gc_z[i];

        grr_z_xxx_yyyz[i] = ts_xxx_yyy[i] * gfe2_0 + gr_xxx_yyy[i] * gfe_0 + ts_xxx_yyyz[i] * gfe_0 * gc_z[i] + gr_xxx_yyyz[i] * gc_z[i];

        grr_z_xxx_yyzz[i] = 2.0 * ts_xxx_yyz[i] * gfe2_0 + 2.0 * gr_xxx_yyz[i] * gfe_0 + ts_xxx_yyzz[i] * gfe_0 * gc_z[i] + gr_xxx_yyzz[i] * gc_z[i];

        grr_z_xxx_yzzz[i] = 3.0 * ts_xxx_yzz[i] * gfe2_0 + 3.0 * gr_xxx_yzz[i] * gfe_0 + ts_xxx_yzzz[i] * gfe_0 * gc_z[i] + gr_xxx_yzzz[i] * gc_z[i];

        grr_z_xxx_zzzz[i] = 4.0 * ts_xxx_zzz[i] * gfe2_0 + 4.0 * gr_xxx_zzz[i] * gfe_0 + ts_xxx_zzzz[i] * gfe_0 * gc_z[i] + gr_xxx_zzzz[i] * gc_z[i];
    }

    // Set up 315-330 components of targeted buffer : FG

    auto grr_z_xxy_xxxx = pbuffer.data(idx_gr_fg + 315);

    auto grr_z_xxy_xxxy = pbuffer.data(idx_gr_fg + 316);

    auto grr_z_xxy_xxxz = pbuffer.data(idx_gr_fg + 317);

    auto grr_z_xxy_xxyy = pbuffer.data(idx_gr_fg + 318);

    auto grr_z_xxy_xxyz = pbuffer.data(idx_gr_fg + 319);

    auto grr_z_xxy_xxzz = pbuffer.data(idx_gr_fg + 320);

    auto grr_z_xxy_xyyy = pbuffer.data(idx_gr_fg + 321);

    auto grr_z_xxy_xyyz = pbuffer.data(idx_gr_fg + 322);

    auto grr_z_xxy_xyzz = pbuffer.data(idx_gr_fg + 323);

    auto grr_z_xxy_xzzz = pbuffer.data(idx_gr_fg + 324);

    auto grr_z_xxy_yyyy = pbuffer.data(idx_gr_fg + 325);

    auto grr_z_xxy_yyyz = pbuffer.data(idx_gr_fg + 326);

    auto grr_z_xxy_yyzz = pbuffer.data(idx_gr_fg + 327);

    auto grr_z_xxy_yzzz = pbuffer.data(idx_gr_fg + 328);

    auto grr_z_xxy_zzzz = pbuffer.data(idx_gr_fg + 329);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xxx, gr_xxy_xxxx, gr_xxy_xxxy, gr_xxy_xxxz, gr_xxy_xxy, gr_xxy_xxyy, gr_xxy_xxyz, gr_xxy_xxz, gr_xxy_xxzz, gr_xxy_xyy, gr_xxy_xyyy, gr_xxy_xyyz, gr_xxy_xyz, gr_xxy_xyzz, gr_xxy_xzz, gr_xxy_xzzz, gr_xxy_yyy, gr_xxy_yyyy, gr_xxy_yyyz, gr_xxy_yyz, gr_xxy_yyzz, gr_xxy_yzz, gr_xxy_yzzz, gr_xxy_zzz, gr_xxy_zzzz, grr_z_xxy_xxxx, grr_z_xxy_xxxy, grr_z_xxy_xxxz, grr_z_xxy_xxyy, grr_z_xxy_xxyz, grr_z_xxy_xxzz, grr_z_xxy_xyyy, grr_z_xxy_xyyz, grr_z_xxy_xyzz, grr_z_xxy_xzzz, grr_z_xxy_yyyy, grr_z_xxy_yyyz, grr_z_xxy_yyzz, grr_z_xxy_yzzz, grr_z_xxy_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxy_xxxx[i] = ts_xxy_xxxx[i] * gfe_0 * gc_z[i] + gr_xxy_xxxx[i] * gc_z[i];

        grr_z_xxy_xxxy[i] = ts_xxy_xxxy[i] * gfe_0 * gc_z[i] + gr_xxy_xxxy[i] * gc_z[i];

        grr_z_xxy_xxxz[i] = ts_xxy_xxx[i] * gfe2_0 + gr_xxy_xxx[i] * gfe_0 + ts_xxy_xxxz[i] * gfe_0 * gc_z[i] + gr_xxy_xxxz[i] * gc_z[i];

        grr_z_xxy_xxyy[i] = ts_xxy_xxyy[i] * gfe_0 * gc_z[i] + gr_xxy_xxyy[i] * gc_z[i];

        grr_z_xxy_xxyz[i] = ts_xxy_xxy[i] * gfe2_0 + gr_xxy_xxy[i] * gfe_0 + ts_xxy_xxyz[i] * gfe_0 * gc_z[i] + gr_xxy_xxyz[i] * gc_z[i];

        grr_z_xxy_xxzz[i] = 2.0 * ts_xxy_xxz[i] * gfe2_0 + 2.0 * gr_xxy_xxz[i] * gfe_0 + ts_xxy_xxzz[i] * gfe_0 * gc_z[i] + gr_xxy_xxzz[i] * gc_z[i];

        grr_z_xxy_xyyy[i] = ts_xxy_xyyy[i] * gfe_0 * gc_z[i] + gr_xxy_xyyy[i] * gc_z[i];

        grr_z_xxy_xyyz[i] = ts_xxy_xyy[i] * gfe2_0 + gr_xxy_xyy[i] * gfe_0 + ts_xxy_xyyz[i] * gfe_0 * gc_z[i] + gr_xxy_xyyz[i] * gc_z[i];

        grr_z_xxy_xyzz[i] = 2.0 * ts_xxy_xyz[i] * gfe2_0 + 2.0 * gr_xxy_xyz[i] * gfe_0 + ts_xxy_xyzz[i] * gfe_0 * gc_z[i] + gr_xxy_xyzz[i] * gc_z[i];

        grr_z_xxy_xzzz[i] = 3.0 * ts_xxy_xzz[i] * gfe2_0 + 3.0 * gr_xxy_xzz[i] * gfe_0 + ts_xxy_xzzz[i] * gfe_0 * gc_z[i] + gr_xxy_xzzz[i] * gc_z[i];

        grr_z_xxy_yyyy[i] = ts_xxy_yyyy[i] * gfe_0 * gc_z[i] + gr_xxy_yyyy[i] * gc_z[i];

        grr_z_xxy_yyyz[i] = ts_xxy_yyy[i] * gfe2_0 + gr_xxy_yyy[i] * gfe_0 + ts_xxy_yyyz[i] * gfe_0 * gc_z[i] + gr_xxy_yyyz[i] * gc_z[i];

        grr_z_xxy_yyzz[i] = 2.0 * ts_xxy_yyz[i] * gfe2_0 + 2.0 * gr_xxy_yyz[i] * gfe_0 + ts_xxy_yyzz[i] * gfe_0 * gc_z[i] + gr_xxy_yyzz[i] * gc_z[i];

        grr_z_xxy_yzzz[i] = 3.0 * ts_xxy_yzz[i] * gfe2_0 + 3.0 * gr_xxy_yzz[i] * gfe_0 + ts_xxy_yzzz[i] * gfe_0 * gc_z[i] + gr_xxy_yzzz[i] * gc_z[i];

        grr_z_xxy_zzzz[i] = 4.0 * ts_xxy_zzz[i] * gfe2_0 + 4.0 * gr_xxy_zzz[i] * gfe_0 + ts_xxy_zzzz[i] * gfe_0 * gc_z[i] + gr_xxy_zzzz[i] * gc_z[i];
    }

    // Set up 330-345 components of targeted buffer : FG

    auto grr_z_xxz_xxxx = pbuffer.data(idx_gr_fg + 330);

    auto grr_z_xxz_xxxy = pbuffer.data(idx_gr_fg + 331);

    auto grr_z_xxz_xxxz = pbuffer.data(idx_gr_fg + 332);

    auto grr_z_xxz_xxyy = pbuffer.data(idx_gr_fg + 333);

    auto grr_z_xxz_xxyz = pbuffer.data(idx_gr_fg + 334);

    auto grr_z_xxz_xxzz = pbuffer.data(idx_gr_fg + 335);

    auto grr_z_xxz_xyyy = pbuffer.data(idx_gr_fg + 336);

    auto grr_z_xxz_xyyz = pbuffer.data(idx_gr_fg + 337);

    auto grr_z_xxz_xyzz = pbuffer.data(idx_gr_fg + 338);

    auto grr_z_xxz_xzzz = pbuffer.data(idx_gr_fg + 339);

    auto grr_z_xxz_yyyy = pbuffer.data(idx_gr_fg + 340);

    auto grr_z_xxz_yyyz = pbuffer.data(idx_gr_fg + 341);

    auto grr_z_xxz_yyzz = pbuffer.data(idx_gr_fg + 342);

    auto grr_z_xxz_yzzz = pbuffer.data(idx_gr_fg + 343);

    auto grr_z_xxz_zzzz = pbuffer.data(idx_gr_fg + 344);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxzz, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyzz, gr_xx_xzzz, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyzz, gr_xx_yzzz, gr_xx_zzzz, gr_xxz_xxx, gr_xxz_xxxx, gr_xxz_xxxy, gr_xxz_xxxz, gr_xxz_xxy, gr_xxz_xxyy, gr_xxz_xxyz, gr_xxz_xxz, gr_xxz_xxzz, gr_xxz_xyy, gr_xxz_xyyy, gr_xxz_xyyz, gr_xxz_xyz, gr_xxz_xyzz, gr_xxz_xzz, gr_xxz_xzzz, gr_xxz_yyy, gr_xxz_yyyy, gr_xxz_yyyz, gr_xxz_yyz, gr_xxz_yyzz, gr_xxz_yzz, gr_xxz_yzzz, gr_xxz_zzz, gr_xxz_zzzz, grr_z_xxz_xxxx, grr_z_xxz_xxxy, grr_z_xxz_xxxz, grr_z_xxz_xxyy, grr_z_xxz_xxyz, grr_z_xxz_xxzz, grr_z_xxz_xyyy, grr_z_xxz_xyyz, grr_z_xxz_xyzz, grr_z_xxz_xzzz, grr_z_xxz_yyyy, grr_z_xxz_yyyz, grr_z_xxz_yyzz, grr_z_xxz_yzzz, grr_z_xxz_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxz_xxxx[i] = ts_xx_xxxx[i] * gfe2_0 + gr_xx_xxxx[i] * gfe_0 + ts_xxz_xxxx[i] * gfe_0 * gc_z[i] + gr_xxz_xxxx[i] * gc_z[i];

        grr_z_xxz_xxxy[i] = ts_xx_xxxy[i] * gfe2_0 + gr_xx_xxxy[i] * gfe_0 + ts_xxz_xxxy[i] * gfe_0 * gc_z[i] + gr_xxz_xxxy[i] * gc_z[i];

        grr_z_xxz_xxxz[i] = ts_xx_xxxz[i] * gfe2_0 + gr_xx_xxxz[i] * gfe_0 + ts_xxz_xxx[i] * gfe2_0 + gr_xxz_xxx[i] * gfe_0 + ts_xxz_xxxz[i] * gfe_0 * gc_z[i] + gr_xxz_xxxz[i] * gc_z[i];

        grr_z_xxz_xxyy[i] = ts_xx_xxyy[i] * gfe2_0 + gr_xx_xxyy[i] * gfe_0 + ts_xxz_xxyy[i] * gfe_0 * gc_z[i] + gr_xxz_xxyy[i] * gc_z[i];

        grr_z_xxz_xxyz[i] = ts_xx_xxyz[i] * gfe2_0 + gr_xx_xxyz[i] * gfe_0 + ts_xxz_xxy[i] * gfe2_0 + gr_xxz_xxy[i] * gfe_0 + ts_xxz_xxyz[i] * gfe_0 * gc_z[i] + gr_xxz_xxyz[i] * gc_z[i];

        grr_z_xxz_xxzz[i] = ts_xx_xxzz[i] * gfe2_0 + gr_xx_xxzz[i] * gfe_0 + 2.0 * ts_xxz_xxz[i] * gfe2_0 + 2.0 * gr_xxz_xxz[i] * gfe_0 + ts_xxz_xxzz[i] * gfe_0 * gc_z[i] + gr_xxz_xxzz[i] * gc_z[i];

        grr_z_xxz_xyyy[i] = ts_xx_xyyy[i] * gfe2_0 + gr_xx_xyyy[i] * gfe_0 + ts_xxz_xyyy[i] * gfe_0 * gc_z[i] + gr_xxz_xyyy[i] * gc_z[i];

        grr_z_xxz_xyyz[i] = ts_xx_xyyz[i] * gfe2_0 + gr_xx_xyyz[i] * gfe_0 + ts_xxz_xyy[i] * gfe2_0 + gr_xxz_xyy[i] * gfe_0 + ts_xxz_xyyz[i] * gfe_0 * gc_z[i] + gr_xxz_xyyz[i] * gc_z[i];

        grr_z_xxz_xyzz[i] = ts_xx_xyzz[i] * gfe2_0 + gr_xx_xyzz[i] * gfe_0 + 2.0 * ts_xxz_xyz[i] * gfe2_0 + 2.0 * gr_xxz_xyz[i] * gfe_0 + ts_xxz_xyzz[i] * gfe_0 * gc_z[i] + gr_xxz_xyzz[i] * gc_z[i];

        grr_z_xxz_xzzz[i] = ts_xx_xzzz[i] * gfe2_0 + gr_xx_xzzz[i] * gfe_0 + 3.0 * ts_xxz_xzz[i] * gfe2_0 + 3.0 * gr_xxz_xzz[i] * gfe_0 + ts_xxz_xzzz[i] * gfe_0 * gc_z[i] + gr_xxz_xzzz[i] * gc_z[i];

        grr_z_xxz_yyyy[i] = ts_xx_yyyy[i] * gfe2_0 + gr_xx_yyyy[i] * gfe_0 + ts_xxz_yyyy[i] * gfe_0 * gc_z[i] + gr_xxz_yyyy[i] * gc_z[i];

        grr_z_xxz_yyyz[i] = ts_xx_yyyz[i] * gfe2_0 + gr_xx_yyyz[i] * gfe_0 + ts_xxz_yyy[i] * gfe2_0 + gr_xxz_yyy[i] * gfe_0 + ts_xxz_yyyz[i] * gfe_0 * gc_z[i] + gr_xxz_yyyz[i] * gc_z[i];

        grr_z_xxz_yyzz[i] = ts_xx_yyzz[i] * gfe2_0 + gr_xx_yyzz[i] * gfe_0 + 2.0 * ts_xxz_yyz[i] * gfe2_0 + 2.0 * gr_xxz_yyz[i] * gfe_0 + ts_xxz_yyzz[i] * gfe_0 * gc_z[i] + gr_xxz_yyzz[i] * gc_z[i];

        grr_z_xxz_yzzz[i] = ts_xx_yzzz[i] * gfe2_0 + gr_xx_yzzz[i] * gfe_0 + 3.0 * ts_xxz_yzz[i] * gfe2_0 + 3.0 * gr_xxz_yzz[i] * gfe_0 + ts_xxz_yzzz[i] * gfe_0 * gc_z[i] + gr_xxz_yzzz[i] * gc_z[i];

        grr_z_xxz_zzzz[i] = ts_xx_zzzz[i] * gfe2_0 + gr_xx_zzzz[i] * gfe_0 + 4.0 * ts_xxz_zzz[i] * gfe2_0 + 4.0 * gr_xxz_zzz[i] * gfe_0 + ts_xxz_zzzz[i] * gfe_0 * gc_z[i] + gr_xxz_zzzz[i] * gc_z[i];
    }

    // Set up 345-360 components of targeted buffer : FG

    auto grr_z_xyy_xxxx = pbuffer.data(idx_gr_fg + 345);

    auto grr_z_xyy_xxxy = pbuffer.data(idx_gr_fg + 346);

    auto grr_z_xyy_xxxz = pbuffer.data(idx_gr_fg + 347);

    auto grr_z_xyy_xxyy = pbuffer.data(idx_gr_fg + 348);

    auto grr_z_xyy_xxyz = pbuffer.data(idx_gr_fg + 349);

    auto grr_z_xyy_xxzz = pbuffer.data(idx_gr_fg + 350);

    auto grr_z_xyy_xyyy = pbuffer.data(idx_gr_fg + 351);

    auto grr_z_xyy_xyyz = pbuffer.data(idx_gr_fg + 352);

    auto grr_z_xyy_xyzz = pbuffer.data(idx_gr_fg + 353);

    auto grr_z_xyy_xzzz = pbuffer.data(idx_gr_fg + 354);

    auto grr_z_xyy_yyyy = pbuffer.data(idx_gr_fg + 355);

    auto grr_z_xyy_yyyz = pbuffer.data(idx_gr_fg + 356);

    auto grr_z_xyy_yyzz = pbuffer.data(idx_gr_fg + 357);

    auto grr_z_xyy_yzzz = pbuffer.data(idx_gr_fg + 358);

    auto grr_z_xyy_zzzz = pbuffer.data(idx_gr_fg + 359);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xxx, gr_xyy_xxxx, gr_xyy_xxxy, gr_xyy_xxxz, gr_xyy_xxy, gr_xyy_xxyy, gr_xyy_xxyz, gr_xyy_xxz, gr_xyy_xxzz, gr_xyy_xyy, gr_xyy_xyyy, gr_xyy_xyyz, gr_xyy_xyz, gr_xyy_xyzz, gr_xyy_xzz, gr_xyy_xzzz, gr_xyy_yyy, gr_xyy_yyyy, gr_xyy_yyyz, gr_xyy_yyz, gr_xyy_yyzz, gr_xyy_yzz, gr_xyy_yzzz, gr_xyy_zzz, gr_xyy_zzzz, grr_z_xyy_xxxx, grr_z_xyy_xxxy, grr_z_xyy_xxxz, grr_z_xyy_xxyy, grr_z_xyy_xxyz, grr_z_xyy_xxzz, grr_z_xyy_xyyy, grr_z_xyy_xyyz, grr_z_xyy_xyzz, grr_z_xyy_xzzz, grr_z_xyy_yyyy, grr_z_xyy_yyyz, grr_z_xyy_yyzz, grr_z_xyy_yzzz, grr_z_xyy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyy_xxxx[i] = ts_xyy_xxxx[i] * gfe_0 * gc_z[i] + gr_xyy_xxxx[i] * gc_z[i];

        grr_z_xyy_xxxy[i] = ts_xyy_xxxy[i] * gfe_0 * gc_z[i] + gr_xyy_xxxy[i] * gc_z[i];

        grr_z_xyy_xxxz[i] = ts_xyy_xxx[i] * gfe2_0 + gr_xyy_xxx[i] * gfe_0 + ts_xyy_xxxz[i] * gfe_0 * gc_z[i] + gr_xyy_xxxz[i] * gc_z[i];

        grr_z_xyy_xxyy[i] = ts_xyy_xxyy[i] * gfe_0 * gc_z[i] + gr_xyy_xxyy[i] * gc_z[i];

        grr_z_xyy_xxyz[i] = ts_xyy_xxy[i] * gfe2_0 + gr_xyy_xxy[i] * gfe_0 + ts_xyy_xxyz[i] * gfe_0 * gc_z[i] + gr_xyy_xxyz[i] * gc_z[i];

        grr_z_xyy_xxzz[i] = 2.0 * ts_xyy_xxz[i] * gfe2_0 + 2.0 * gr_xyy_xxz[i] * gfe_0 + ts_xyy_xxzz[i] * gfe_0 * gc_z[i] + gr_xyy_xxzz[i] * gc_z[i];

        grr_z_xyy_xyyy[i] = ts_xyy_xyyy[i] * gfe_0 * gc_z[i] + gr_xyy_xyyy[i] * gc_z[i];

        grr_z_xyy_xyyz[i] = ts_xyy_xyy[i] * gfe2_0 + gr_xyy_xyy[i] * gfe_0 + ts_xyy_xyyz[i] * gfe_0 * gc_z[i] + gr_xyy_xyyz[i] * gc_z[i];

        grr_z_xyy_xyzz[i] = 2.0 * ts_xyy_xyz[i] * gfe2_0 + 2.0 * gr_xyy_xyz[i] * gfe_0 + ts_xyy_xyzz[i] * gfe_0 * gc_z[i] + gr_xyy_xyzz[i] * gc_z[i];

        grr_z_xyy_xzzz[i] = 3.0 * ts_xyy_xzz[i] * gfe2_0 + 3.0 * gr_xyy_xzz[i] * gfe_0 + ts_xyy_xzzz[i] * gfe_0 * gc_z[i] + gr_xyy_xzzz[i] * gc_z[i];

        grr_z_xyy_yyyy[i] = ts_xyy_yyyy[i] * gfe_0 * gc_z[i] + gr_xyy_yyyy[i] * gc_z[i];

        grr_z_xyy_yyyz[i] = ts_xyy_yyy[i] * gfe2_0 + gr_xyy_yyy[i] * gfe_0 + ts_xyy_yyyz[i] * gfe_0 * gc_z[i] + gr_xyy_yyyz[i] * gc_z[i];

        grr_z_xyy_yyzz[i] = 2.0 * ts_xyy_yyz[i] * gfe2_0 + 2.0 * gr_xyy_yyz[i] * gfe_0 + ts_xyy_yyzz[i] * gfe_0 * gc_z[i] + gr_xyy_yyzz[i] * gc_z[i];

        grr_z_xyy_yzzz[i] = 3.0 * ts_xyy_yzz[i] * gfe2_0 + 3.0 * gr_xyy_yzz[i] * gfe_0 + ts_xyy_yzzz[i] * gfe_0 * gc_z[i] + gr_xyy_yzzz[i] * gc_z[i];

        grr_z_xyy_zzzz[i] = 4.0 * ts_xyy_zzz[i] * gfe2_0 + 4.0 * gr_xyy_zzz[i] * gfe_0 + ts_xyy_zzzz[i] * gfe_0 * gc_z[i] + gr_xyy_zzzz[i] * gc_z[i];
    }

    // Set up 360-375 components of targeted buffer : FG

    auto grr_z_xyz_xxxx = pbuffer.data(idx_gr_fg + 360);

    auto grr_z_xyz_xxxy = pbuffer.data(idx_gr_fg + 361);

    auto grr_z_xyz_xxxz = pbuffer.data(idx_gr_fg + 362);

    auto grr_z_xyz_xxyy = pbuffer.data(idx_gr_fg + 363);

    auto grr_z_xyz_xxyz = pbuffer.data(idx_gr_fg + 364);

    auto grr_z_xyz_xxzz = pbuffer.data(idx_gr_fg + 365);

    auto grr_z_xyz_xyyy = pbuffer.data(idx_gr_fg + 366);

    auto grr_z_xyz_xyyz = pbuffer.data(idx_gr_fg + 367);

    auto grr_z_xyz_xyzz = pbuffer.data(idx_gr_fg + 368);

    auto grr_z_xyz_xzzz = pbuffer.data(idx_gr_fg + 369);

    auto grr_z_xyz_yyyy = pbuffer.data(idx_gr_fg + 370);

    auto grr_z_xyz_yyyz = pbuffer.data(idx_gr_fg + 371);

    auto grr_z_xyz_yyzz = pbuffer.data(idx_gr_fg + 372);

    auto grr_z_xyz_yzzz = pbuffer.data(idx_gr_fg + 373);

    auto grr_z_xyz_zzzz = pbuffer.data(idx_gr_fg + 374);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxzz, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyzz, gr_xy_xzzz, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyzz, gr_xy_yzzz, gr_xy_zzzz, gr_xyz_xxx, gr_xyz_xxxx, gr_xyz_xxxy, gr_xyz_xxxz, gr_xyz_xxy, gr_xyz_xxyy, gr_xyz_xxyz, gr_xyz_xxz, gr_xyz_xxzz, gr_xyz_xyy, gr_xyz_xyyy, gr_xyz_xyyz, gr_xyz_xyz, gr_xyz_xyzz, gr_xyz_xzz, gr_xyz_xzzz, gr_xyz_yyy, gr_xyz_yyyy, gr_xyz_yyyz, gr_xyz_yyz, gr_xyz_yyzz, gr_xyz_yzz, gr_xyz_yzzz, gr_xyz_zzz, gr_xyz_zzzz, grr_z_xyz_xxxx, grr_z_xyz_xxxy, grr_z_xyz_xxxz, grr_z_xyz_xxyy, grr_z_xyz_xxyz, grr_z_xyz_xxzz, grr_z_xyz_xyyy, grr_z_xyz_xyyz, grr_z_xyz_xyzz, grr_z_xyz_xzzz, grr_z_xyz_yyyy, grr_z_xyz_yyyz, grr_z_xyz_yyzz, grr_z_xyz_yzzz, grr_z_xyz_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyz_xxxx[i] = ts_xy_xxxx[i] * gfe2_0 + gr_xy_xxxx[i] * gfe_0 + ts_xyz_xxxx[i] * gfe_0 * gc_z[i] + gr_xyz_xxxx[i] * gc_z[i];

        grr_z_xyz_xxxy[i] = ts_xy_xxxy[i] * gfe2_0 + gr_xy_xxxy[i] * gfe_0 + ts_xyz_xxxy[i] * gfe_0 * gc_z[i] + gr_xyz_xxxy[i] * gc_z[i];

        grr_z_xyz_xxxz[i] = ts_xy_xxxz[i] * gfe2_0 + gr_xy_xxxz[i] * gfe_0 + ts_xyz_xxx[i] * gfe2_0 + gr_xyz_xxx[i] * gfe_0 + ts_xyz_xxxz[i] * gfe_0 * gc_z[i] + gr_xyz_xxxz[i] * gc_z[i];

        grr_z_xyz_xxyy[i] = ts_xy_xxyy[i] * gfe2_0 + gr_xy_xxyy[i] * gfe_0 + ts_xyz_xxyy[i] * gfe_0 * gc_z[i] + gr_xyz_xxyy[i] * gc_z[i];

        grr_z_xyz_xxyz[i] = ts_xy_xxyz[i] * gfe2_0 + gr_xy_xxyz[i] * gfe_0 + ts_xyz_xxy[i] * gfe2_0 + gr_xyz_xxy[i] * gfe_0 + ts_xyz_xxyz[i] * gfe_0 * gc_z[i] + gr_xyz_xxyz[i] * gc_z[i];

        grr_z_xyz_xxzz[i] = ts_xy_xxzz[i] * gfe2_0 + gr_xy_xxzz[i] * gfe_0 + 2.0 * ts_xyz_xxz[i] * gfe2_0 + 2.0 * gr_xyz_xxz[i] * gfe_0 + ts_xyz_xxzz[i] * gfe_0 * gc_z[i] + gr_xyz_xxzz[i] * gc_z[i];

        grr_z_xyz_xyyy[i] = ts_xy_xyyy[i] * gfe2_0 + gr_xy_xyyy[i] * gfe_0 + ts_xyz_xyyy[i] * gfe_0 * gc_z[i] + gr_xyz_xyyy[i] * gc_z[i];

        grr_z_xyz_xyyz[i] = ts_xy_xyyz[i] * gfe2_0 + gr_xy_xyyz[i] * gfe_0 + ts_xyz_xyy[i] * gfe2_0 + gr_xyz_xyy[i] * gfe_0 + ts_xyz_xyyz[i] * gfe_0 * gc_z[i] + gr_xyz_xyyz[i] * gc_z[i];

        grr_z_xyz_xyzz[i] = ts_xy_xyzz[i] * gfe2_0 + gr_xy_xyzz[i] * gfe_0 + 2.0 * ts_xyz_xyz[i] * gfe2_0 + 2.0 * gr_xyz_xyz[i] * gfe_0 + ts_xyz_xyzz[i] * gfe_0 * gc_z[i] + gr_xyz_xyzz[i] * gc_z[i];

        grr_z_xyz_xzzz[i] = ts_xy_xzzz[i] * gfe2_0 + gr_xy_xzzz[i] * gfe_0 + 3.0 * ts_xyz_xzz[i] * gfe2_0 + 3.0 * gr_xyz_xzz[i] * gfe_0 + ts_xyz_xzzz[i] * gfe_0 * gc_z[i] + gr_xyz_xzzz[i] * gc_z[i];

        grr_z_xyz_yyyy[i] = ts_xy_yyyy[i] * gfe2_0 + gr_xy_yyyy[i] * gfe_0 + ts_xyz_yyyy[i] * gfe_0 * gc_z[i] + gr_xyz_yyyy[i] * gc_z[i];

        grr_z_xyz_yyyz[i] = ts_xy_yyyz[i] * gfe2_0 + gr_xy_yyyz[i] * gfe_0 + ts_xyz_yyy[i] * gfe2_0 + gr_xyz_yyy[i] * gfe_0 + ts_xyz_yyyz[i] * gfe_0 * gc_z[i] + gr_xyz_yyyz[i] * gc_z[i];

        grr_z_xyz_yyzz[i] = ts_xy_yyzz[i] * gfe2_0 + gr_xy_yyzz[i] * gfe_0 + 2.0 * ts_xyz_yyz[i] * gfe2_0 + 2.0 * gr_xyz_yyz[i] * gfe_0 + ts_xyz_yyzz[i] * gfe_0 * gc_z[i] + gr_xyz_yyzz[i] * gc_z[i];

        grr_z_xyz_yzzz[i] = ts_xy_yzzz[i] * gfe2_0 + gr_xy_yzzz[i] * gfe_0 + 3.0 * ts_xyz_yzz[i] * gfe2_0 + 3.0 * gr_xyz_yzz[i] * gfe_0 + ts_xyz_yzzz[i] * gfe_0 * gc_z[i] + gr_xyz_yzzz[i] * gc_z[i];

        grr_z_xyz_zzzz[i] = ts_xy_zzzz[i] * gfe2_0 + gr_xy_zzzz[i] * gfe_0 + 4.0 * ts_xyz_zzz[i] * gfe2_0 + 4.0 * gr_xyz_zzz[i] * gfe_0 + ts_xyz_zzzz[i] * gfe_0 * gc_z[i] + gr_xyz_zzzz[i] * gc_z[i];
    }

    // Set up 375-390 components of targeted buffer : FG

    auto grr_z_xzz_xxxx = pbuffer.data(idx_gr_fg + 375);

    auto grr_z_xzz_xxxy = pbuffer.data(idx_gr_fg + 376);

    auto grr_z_xzz_xxxz = pbuffer.data(idx_gr_fg + 377);

    auto grr_z_xzz_xxyy = pbuffer.data(idx_gr_fg + 378);

    auto grr_z_xzz_xxyz = pbuffer.data(idx_gr_fg + 379);

    auto grr_z_xzz_xxzz = pbuffer.data(idx_gr_fg + 380);

    auto grr_z_xzz_xyyy = pbuffer.data(idx_gr_fg + 381);

    auto grr_z_xzz_xyyz = pbuffer.data(idx_gr_fg + 382);

    auto grr_z_xzz_xyzz = pbuffer.data(idx_gr_fg + 383);

    auto grr_z_xzz_xzzz = pbuffer.data(idx_gr_fg + 384);

    auto grr_z_xzz_yyyy = pbuffer.data(idx_gr_fg + 385);

    auto grr_z_xzz_yyyz = pbuffer.data(idx_gr_fg + 386);

    auto grr_z_xzz_yyzz = pbuffer.data(idx_gr_fg + 387);

    auto grr_z_xzz_yzzz = pbuffer.data(idx_gr_fg + 388);

    auto grr_z_xzz_zzzz = pbuffer.data(idx_gr_fg + 389);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxzz, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyzz, gr_xz_xzzz, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyzz, gr_xz_yzzz, gr_xz_zzzz, gr_xzz_xxx, gr_xzz_xxxx, gr_xzz_xxxy, gr_xzz_xxxz, gr_xzz_xxy, gr_xzz_xxyy, gr_xzz_xxyz, gr_xzz_xxz, gr_xzz_xxzz, gr_xzz_xyy, gr_xzz_xyyy, gr_xzz_xyyz, gr_xzz_xyz, gr_xzz_xyzz, gr_xzz_xzz, gr_xzz_xzzz, gr_xzz_yyy, gr_xzz_yyyy, gr_xzz_yyyz, gr_xzz_yyz, gr_xzz_yyzz, gr_xzz_yzz, gr_xzz_yzzz, gr_xzz_zzz, gr_xzz_zzzz, grr_z_xzz_xxxx, grr_z_xzz_xxxy, grr_z_xzz_xxxz, grr_z_xzz_xxyy, grr_z_xzz_xxyz, grr_z_xzz_xxzz, grr_z_xzz_xyyy, grr_z_xzz_xyyz, grr_z_xzz_xyzz, grr_z_xzz_xzzz, grr_z_xzz_yyyy, grr_z_xzz_yyyz, grr_z_xzz_yyzz, grr_z_xzz_yzzz, grr_z_xzz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xzz_xxxx[i] = 2.0 * ts_xz_xxxx[i] * gfe2_0 + 2.0 * gr_xz_xxxx[i] * gfe_0 + ts_xzz_xxxx[i] * gfe_0 * gc_z[i] + gr_xzz_xxxx[i] * gc_z[i];

        grr_z_xzz_xxxy[i] = 2.0 * ts_xz_xxxy[i] * gfe2_0 + 2.0 * gr_xz_xxxy[i] * gfe_0 + ts_xzz_xxxy[i] * gfe_0 * gc_z[i] + gr_xzz_xxxy[i] * gc_z[i];

        grr_z_xzz_xxxz[i] = 2.0 * ts_xz_xxxz[i] * gfe2_0 + 2.0 * gr_xz_xxxz[i] * gfe_0 + ts_xzz_xxx[i] * gfe2_0 + gr_xzz_xxx[i] * gfe_0 + ts_xzz_xxxz[i] * gfe_0 * gc_z[i] + gr_xzz_xxxz[i] * gc_z[i];

        grr_z_xzz_xxyy[i] = 2.0 * ts_xz_xxyy[i] * gfe2_0 + 2.0 * gr_xz_xxyy[i] * gfe_0 + ts_xzz_xxyy[i] * gfe_0 * gc_z[i] + gr_xzz_xxyy[i] * gc_z[i];

        grr_z_xzz_xxyz[i] = 2.0 * ts_xz_xxyz[i] * gfe2_0 + 2.0 * gr_xz_xxyz[i] * gfe_0 + ts_xzz_xxy[i] * gfe2_0 + gr_xzz_xxy[i] * gfe_0 + ts_xzz_xxyz[i] * gfe_0 * gc_z[i] + gr_xzz_xxyz[i] * gc_z[i];

        grr_z_xzz_xxzz[i] = 2.0 * ts_xz_xxzz[i] * gfe2_0 + 2.0 * gr_xz_xxzz[i] * gfe_0 + 2.0 * ts_xzz_xxz[i] * gfe2_0 + 2.0 * gr_xzz_xxz[i] * gfe_0 + ts_xzz_xxzz[i] * gfe_0 * gc_z[i] + gr_xzz_xxzz[i] * gc_z[i];

        grr_z_xzz_xyyy[i] = 2.0 * ts_xz_xyyy[i] * gfe2_0 + 2.0 * gr_xz_xyyy[i] * gfe_0 + ts_xzz_xyyy[i] * gfe_0 * gc_z[i] + gr_xzz_xyyy[i] * gc_z[i];

        grr_z_xzz_xyyz[i] = 2.0 * ts_xz_xyyz[i] * gfe2_0 + 2.0 * gr_xz_xyyz[i] * gfe_0 + ts_xzz_xyy[i] * gfe2_0 + gr_xzz_xyy[i] * gfe_0 + ts_xzz_xyyz[i] * gfe_0 * gc_z[i] + gr_xzz_xyyz[i] * gc_z[i];

        grr_z_xzz_xyzz[i] = 2.0 * ts_xz_xyzz[i] * gfe2_0 + 2.0 * gr_xz_xyzz[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe2_0 + 2.0 * gr_xzz_xyz[i] * gfe_0 + ts_xzz_xyzz[i] * gfe_0 * gc_z[i] + gr_xzz_xyzz[i] * gc_z[i];

        grr_z_xzz_xzzz[i] = 2.0 * ts_xz_xzzz[i] * gfe2_0 + 2.0 * gr_xz_xzzz[i] * gfe_0 + 3.0 * ts_xzz_xzz[i] * gfe2_0 + 3.0 * gr_xzz_xzz[i] * gfe_0 + ts_xzz_xzzz[i] * gfe_0 * gc_z[i] + gr_xzz_xzzz[i] * gc_z[i];

        grr_z_xzz_yyyy[i] = 2.0 * ts_xz_yyyy[i] * gfe2_0 + 2.0 * gr_xz_yyyy[i] * gfe_0 + ts_xzz_yyyy[i] * gfe_0 * gc_z[i] + gr_xzz_yyyy[i] * gc_z[i];

        grr_z_xzz_yyyz[i] = 2.0 * ts_xz_yyyz[i] * gfe2_0 + 2.0 * gr_xz_yyyz[i] * gfe_0 + ts_xzz_yyy[i] * gfe2_0 + gr_xzz_yyy[i] * gfe_0 + ts_xzz_yyyz[i] * gfe_0 * gc_z[i] + gr_xzz_yyyz[i] * gc_z[i];

        grr_z_xzz_yyzz[i] = 2.0 * ts_xz_yyzz[i] * gfe2_0 + 2.0 * gr_xz_yyzz[i] * gfe_0 + 2.0 * ts_xzz_yyz[i] * gfe2_0 + 2.0 * gr_xzz_yyz[i] * gfe_0 + ts_xzz_yyzz[i] * gfe_0 * gc_z[i] + gr_xzz_yyzz[i] * gc_z[i];

        grr_z_xzz_yzzz[i] = 2.0 * ts_xz_yzzz[i] * gfe2_0 + 2.0 * gr_xz_yzzz[i] * gfe_0 + 3.0 * ts_xzz_yzz[i] * gfe2_0 + 3.0 * gr_xzz_yzz[i] * gfe_0 + ts_xzz_yzzz[i] * gfe_0 * gc_z[i] + gr_xzz_yzzz[i] * gc_z[i];

        grr_z_xzz_zzzz[i] = 2.0 * ts_xz_zzzz[i] * gfe2_0 + 2.0 * gr_xz_zzzz[i] * gfe_0 + 4.0 * ts_xzz_zzz[i] * gfe2_0 + 4.0 * gr_xzz_zzz[i] * gfe_0 + ts_xzz_zzzz[i] * gfe_0 * gc_z[i] + gr_xzz_zzzz[i] * gc_z[i];
    }

    // Set up 390-405 components of targeted buffer : FG

    auto grr_z_yyy_xxxx = pbuffer.data(idx_gr_fg + 390);

    auto grr_z_yyy_xxxy = pbuffer.data(idx_gr_fg + 391);

    auto grr_z_yyy_xxxz = pbuffer.data(idx_gr_fg + 392);

    auto grr_z_yyy_xxyy = pbuffer.data(idx_gr_fg + 393);

    auto grr_z_yyy_xxyz = pbuffer.data(idx_gr_fg + 394);

    auto grr_z_yyy_xxzz = pbuffer.data(idx_gr_fg + 395);

    auto grr_z_yyy_xyyy = pbuffer.data(idx_gr_fg + 396);

    auto grr_z_yyy_xyyz = pbuffer.data(idx_gr_fg + 397);

    auto grr_z_yyy_xyzz = pbuffer.data(idx_gr_fg + 398);

    auto grr_z_yyy_xzzz = pbuffer.data(idx_gr_fg + 399);

    auto grr_z_yyy_yyyy = pbuffer.data(idx_gr_fg + 400);

    auto grr_z_yyy_yyyz = pbuffer.data(idx_gr_fg + 401);

    auto grr_z_yyy_yyzz = pbuffer.data(idx_gr_fg + 402);

    auto grr_z_yyy_yzzz = pbuffer.data(idx_gr_fg + 403);

    auto grr_z_yyy_zzzz = pbuffer.data(idx_gr_fg + 404);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xxx, gr_yyy_xxxx, gr_yyy_xxxy, gr_yyy_xxxz, gr_yyy_xxy, gr_yyy_xxyy, gr_yyy_xxyz, gr_yyy_xxz, gr_yyy_xxzz, gr_yyy_xyy, gr_yyy_xyyy, gr_yyy_xyyz, gr_yyy_xyz, gr_yyy_xyzz, gr_yyy_xzz, gr_yyy_xzzz, gr_yyy_yyy, gr_yyy_yyyy, gr_yyy_yyyz, gr_yyy_yyz, gr_yyy_yyzz, gr_yyy_yzz, gr_yyy_yzzz, gr_yyy_zzz, gr_yyy_zzzz, grr_z_yyy_xxxx, grr_z_yyy_xxxy, grr_z_yyy_xxxz, grr_z_yyy_xxyy, grr_z_yyy_xxyz, grr_z_yyy_xxzz, grr_z_yyy_xyyy, grr_z_yyy_xyyz, grr_z_yyy_xyzz, grr_z_yyy_xzzz, grr_z_yyy_yyyy, grr_z_yyy_yyyz, grr_z_yyy_yyzz, grr_z_yyy_yzzz, grr_z_yyy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyy_xxxx[i] = ts_yyy_xxxx[i] * gfe_0 * gc_z[i] + gr_yyy_xxxx[i] * gc_z[i];

        grr_z_yyy_xxxy[i] = ts_yyy_xxxy[i] * gfe_0 * gc_z[i] + gr_yyy_xxxy[i] * gc_z[i];

        grr_z_yyy_xxxz[i] = ts_yyy_xxx[i] * gfe2_0 + gr_yyy_xxx[i] * gfe_0 + ts_yyy_xxxz[i] * gfe_0 * gc_z[i] + gr_yyy_xxxz[i] * gc_z[i];

        grr_z_yyy_xxyy[i] = ts_yyy_xxyy[i] * gfe_0 * gc_z[i] + gr_yyy_xxyy[i] * gc_z[i];

        grr_z_yyy_xxyz[i] = ts_yyy_xxy[i] * gfe2_0 + gr_yyy_xxy[i] * gfe_0 + ts_yyy_xxyz[i] * gfe_0 * gc_z[i] + gr_yyy_xxyz[i] * gc_z[i];

        grr_z_yyy_xxzz[i] = 2.0 * ts_yyy_xxz[i] * gfe2_0 + 2.0 * gr_yyy_xxz[i] * gfe_0 + ts_yyy_xxzz[i] * gfe_0 * gc_z[i] + gr_yyy_xxzz[i] * gc_z[i];

        grr_z_yyy_xyyy[i] = ts_yyy_xyyy[i] * gfe_0 * gc_z[i] + gr_yyy_xyyy[i] * gc_z[i];

        grr_z_yyy_xyyz[i] = ts_yyy_xyy[i] * gfe2_0 + gr_yyy_xyy[i] * gfe_0 + ts_yyy_xyyz[i] * gfe_0 * gc_z[i] + gr_yyy_xyyz[i] * gc_z[i];

        grr_z_yyy_xyzz[i] = 2.0 * ts_yyy_xyz[i] * gfe2_0 + 2.0 * gr_yyy_xyz[i] * gfe_0 + ts_yyy_xyzz[i] * gfe_0 * gc_z[i] + gr_yyy_xyzz[i] * gc_z[i];

        grr_z_yyy_xzzz[i] = 3.0 * ts_yyy_xzz[i] * gfe2_0 + 3.0 * gr_yyy_xzz[i] * gfe_0 + ts_yyy_xzzz[i] * gfe_0 * gc_z[i] + gr_yyy_xzzz[i] * gc_z[i];

        grr_z_yyy_yyyy[i] = ts_yyy_yyyy[i] * gfe_0 * gc_z[i] + gr_yyy_yyyy[i] * gc_z[i];

        grr_z_yyy_yyyz[i] = ts_yyy_yyy[i] * gfe2_0 + gr_yyy_yyy[i] * gfe_0 + ts_yyy_yyyz[i] * gfe_0 * gc_z[i] + gr_yyy_yyyz[i] * gc_z[i];

        grr_z_yyy_yyzz[i] = 2.0 * ts_yyy_yyz[i] * gfe2_0 + 2.0 * gr_yyy_yyz[i] * gfe_0 + ts_yyy_yyzz[i] * gfe_0 * gc_z[i] + gr_yyy_yyzz[i] * gc_z[i];

        grr_z_yyy_yzzz[i] = 3.0 * ts_yyy_yzz[i] * gfe2_0 + 3.0 * gr_yyy_yzz[i] * gfe_0 + ts_yyy_yzzz[i] * gfe_0 * gc_z[i] + gr_yyy_yzzz[i] * gc_z[i];

        grr_z_yyy_zzzz[i] = 4.0 * ts_yyy_zzz[i] * gfe2_0 + 4.0 * gr_yyy_zzz[i] * gfe_0 + ts_yyy_zzzz[i] * gfe_0 * gc_z[i] + gr_yyy_zzzz[i] * gc_z[i];
    }

    // Set up 405-420 components of targeted buffer : FG

    auto grr_z_yyz_xxxx = pbuffer.data(idx_gr_fg + 405);

    auto grr_z_yyz_xxxy = pbuffer.data(idx_gr_fg + 406);

    auto grr_z_yyz_xxxz = pbuffer.data(idx_gr_fg + 407);

    auto grr_z_yyz_xxyy = pbuffer.data(idx_gr_fg + 408);

    auto grr_z_yyz_xxyz = pbuffer.data(idx_gr_fg + 409);

    auto grr_z_yyz_xxzz = pbuffer.data(idx_gr_fg + 410);

    auto grr_z_yyz_xyyy = pbuffer.data(idx_gr_fg + 411);

    auto grr_z_yyz_xyyz = pbuffer.data(idx_gr_fg + 412);

    auto grr_z_yyz_xyzz = pbuffer.data(idx_gr_fg + 413);

    auto grr_z_yyz_xzzz = pbuffer.data(idx_gr_fg + 414);

    auto grr_z_yyz_yyyy = pbuffer.data(idx_gr_fg + 415);

    auto grr_z_yyz_yyyz = pbuffer.data(idx_gr_fg + 416);

    auto grr_z_yyz_yyzz = pbuffer.data(idx_gr_fg + 417);

    auto grr_z_yyz_yzzz = pbuffer.data(idx_gr_fg + 418);

    auto grr_z_yyz_zzzz = pbuffer.data(idx_gr_fg + 419);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxzz, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyzz, gr_yy_xzzz, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyzz, gr_yy_yzzz, gr_yy_zzzz, gr_yyz_xxx, gr_yyz_xxxx, gr_yyz_xxxy, gr_yyz_xxxz, gr_yyz_xxy, gr_yyz_xxyy, gr_yyz_xxyz, gr_yyz_xxz, gr_yyz_xxzz, gr_yyz_xyy, gr_yyz_xyyy, gr_yyz_xyyz, gr_yyz_xyz, gr_yyz_xyzz, gr_yyz_xzz, gr_yyz_xzzz, gr_yyz_yyy, gr_yyz_yyyy, gr_yyz_yyyz, gr_yyz_yyz, gr_yyz_yyzz, gr_yyz_yzz, gr_yyz_yzzz, gr_yyz_zzz, gr_yyz_zzzz, grr_z_yyz_xxxx, grr_z_yyz_xxxy, grr_z_yyz_xxxz, grr_z_yyz_xxyy, grr_z_yyz_xxyz, grr_z_yyz_xxzz, grr_z_yyz_xyyy, grr_z_yyz_xyyz, grr_z_yyz_xyzz, grr_z_yyz_xzzz, grr_z_yyz_yyyy, grr_z_yyz_yyyz, grr_z_yyz_yyzz, grr_z_yyz_yzzz, grr_z_yyz_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyz_xxxx[i] = ts_yy_xxxx[i] * gfe2_0 + gr_yy_xxxx[i] * gfe_0 + ts_yyz_xxxx[i] * gfe_0 * gc_z[i] + gr_yyz_xxxx[i] * gc_z[i];

        grr_z_yyz_xxxy[i] = ts_yy_xxxy[i] * gfe2_0 + gr_yy_xxxy[i] * gfe_0 + ts_yyz_xxxy[i] * gfe_0 * gc_z[i] + gr_yyz_xxxy[i] * gc_z[i];

        grr_z_yyz_xxxz[i] = ts_yy_xxxz[i] * gfe2_0 + gr_yy_xxxz[i] * gfe_0 + ts_yyz_xxx[i] * gfe2_0 + gr_yyz_xxx[i] * gfe_0 + ts_yyz_xxxz[i] * gfe_0 * gc_z[i] + gr_yyz_xxxz[i] * gc_z[i];

        grr_z_yyz_xxyy[i] = ts_yy_xxyy[i] * gfe2_0 + gr_yy_xxyy[i] * gfe_0 + ts_yyz_xxyy[i] * gfe_0 * gc_z[i] + gr_yyz_xxyy[i] * gc_z[i];

        grr_z_yyz_xxyz[i] = ts_yy_xxyz[i] * gfe2_0 + gr_yy_xxyz[i] * gfe_0 + ts_yyz_xxy[i] * gfe2_0 + gr_yyz_xxy[i] * gfe_0 + ts_yyz_xxyz[i] * gfe_0 * gc_z[i] + gr_yyz_xxyz[i] * gc_z[i];

        grr_z_yyz_xxzz[i] = ts_yy_xxzz[i] * gfe2_0 + gr_yy_xxzz[i] * gfe_0 + 2.0 * ts_yyz_xxz[i] * gfe2_0 + 2.0 * gr_yyz_xxz[i] * gfe_0 + ts_yyz_xxzz[i] * gfe_0 * gc_z[i] + gr_yyz_xxzz[i] * gc_z[i];

        grr_z_yyz_xyyy[i] = ts_yy_xyyy[i] * gfe2_0 + gr_yy_xyyy[i] * gfe_0 + ts_yyz_xyyy[i] * gfe_0 * gc_z[i] + gr_yyz_xyyy[i] * gc_z[i];

        grr_z_yyz_xyyz[i] = ts_yy_xyyz[i] * gfe2_0 + gr_yy_xyyz[i] * gfe_0 + ts_yyz_xyy[i] * gfe2_0 + gr_yyz_xyy[i] * gfe_0 + ts_yyz_xyyz[i] * gfe_0 * gc_z[i] + gr_yyz_xyyz[i] * gc_z[i];

        grr_z_yyz_xyzz[i] = ts_yy_xyzz[i] * gfe2_0 + gr_yy_xyzz[i] * gfe_0 + 2.0 * ts_yyz_xyz[i] * gfe2_0 + 2.0 * gr_yyz_xyz[i] * gfe_0 + ts_yyz_xyzz[i] * gfe_0 * gc_z[i] + gr_yyz_xyzz[i] * gc_z[i];

        grr_z_yyz_xzzz[i] = ts_yy_xzzz[i] * gfe2_0 + gr_yy_xzzz[i] * gfe_0 + 3.0 * ts_yyz_xzz[i] * gfe2_0 + 3.0 * gr_yyz_xzz[i] * gfe_0 + ts_yyz_xzzz[i] * gfe_0 * gc_z[i] + gr_yyz_xzzz[i] * gc_z[i];

        grr_z_yyz_yyyy[i] = ts_yy_yyyy[i] * gfe2_0 + gr_yy_yyyy[i] * gfe_0 + ts_yyz_yyyy[i] * gfe_0 * gc_z[i] + gr_yyz_yyyy[i] * gc_z[i];

        grr_z_yyz_yyyz[i] = ts_yy_yyyz[i] * gfe2_0 + gr_yy_yyyz[i] * gfe_0 + ts_yyz_yyy[i] * gfe2_0 + gr_yyz_yyy[i] * gfe_0 + ts_yyz_yyyz[i] * gfe_0 * gc_z[i] + gr_yyz_yyyz[i] * gc_z[i];

        grr_z_yyz_yyzz[i] = ts_yy_yyzz[i] * gfe2_0 + gr_yy_yyzz[i] * gfe_0 + 2.0 * ts_yyz_yyz[i] * gfe2_0 + 2.0 * gr_yyz_yyz[i] * gfe_0 + ts_yyz_yyzz[i] * gfe_0 * gc_z[i] + gr_yyz_yyzz[i] * gc_z[i];

        grr_z_yyz_yzzz[i] = ts_yy_yzzz[i] * gfe2_0 + gr_yy_yzzz[i] * gfe_0 + 3.0 * ts_yyz_yzz[i] * gfe2_0 + 3.0 * gr_yyz_yzz[i] * gfe_0 + ts_yyz_yzzz[i] * gfe_0 * gc_z[i] + gr_yyz_yzzz[i] * gc_z[i];

        grr_z_yyz_zzzz[i] = ts_yy_zzzz[i] * gfe2_0 + gr_yy_zzzz[i] * gfe_0 + 4.0 * ts_yyz_zzz[i] * gfe2_0 + 4.0 * gr_yyz_zzz[i] * gfe_0 + ts_yyz_zzzz[i] * gfe_0 * gc_z[i] + gr_yyz_zzzz[i] * gc_z[i];
    }

    // Set up 420-435 components of targeted buffer : FG

    auto grr_z_yzz_xxxx = pbuffer.data(idx_gr_fg + 420);

    auto grr_z_yzz_xxxy = pbuffer.data(idx_gr_fg + 421);

    auto grr_z_yzz_xxxz = pbuffer.data(idx_gr_fg + 422);

    auto grr_z_yzz_xxyy = pbuffer.data(idx_gr_fg + 423);

    auto grr_z_yzz_xxyz = pbuffer.data(idx_gr_fg + 424);

    auto grr_z_yzz_xxzz = pbuffer.data(idx_gr_fg + 425);

    auto grr_z_yzz_xyyy = pbuffer.data(idx_gr_fg + 426);

    auto grr_z_yzz_xyyz = pbuffer.data(idx_gr_fg + 427);

    auto grr_z_yzz_xyzz = pbuffer.data(idx_gr_fg + 428);

    auto grr_z_yzz_xzzz = pbuffer.data(idx_gr_fg + 429);

    auto grr_z_yzz_yyyy = pbuffer.data(idx_gr_fg + 430);

    auto grr_z_yzz_yyyz = pbuffer.data(idx_gr_fg + 431);

    auto grr_z_yzz_yyzz = pbuffer.data(idx_gr_fg + 432);

    auto grr_z_yzz_yzzz = pbuffer.data(idx_gr_fg + 433);

    auto grr_z_yzz_zzzz = pbuffer.data(idx_gr_fg + 434);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxzz, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyzz, gr_yz_xzzz, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyzz, gr_yz_yzzz, gr_yz_zzzz, gr_yzz_xxx, gr_yzz_xxxx, gr_yzz_xxxy, gr_yzz_xxxz, gr_yzz_xxy, gr_yzz_xxyy, gr_yzz_xxyz, gr_yzz_xxz, gr_yzz_xxzz, gr_yzz_xyy, gr_yzz_xyyy, gr_yzz_xyyz, gr_yzz_xyz, gr_yzz_xyzz, gr_yzz_xzz, gr_yzz_xzzz, gr_yzz_yyy, gr_yzz_yyyy, gr_yzz_yyyz, gr_yzz_yyz, gr_yzz_yyzz, gr_yzz_yzz, gr_yzz_yzzz, gr_yzz_zzz, gr_yzz_zzzz, grr_z_yzz_xxxx, grr_z_yzz_xxxy, grr_z_yzz_xxxz, grr_z_yzz_xxyy, grr_z_yzz_xxyz, grr_z_yzz_xxzz, grr_z_yzz_xyyy, grr_z_yzz_xyyz, grr_z_yzz_xyzz, grr_z_yzz_xzzz, grr_z_yzz_yyyy, grr_z_yzz_yyyz, grr_z_yzz_yyzz, grr_z_yzz_yzzz, grr_z_yzz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yzz_xxxx[i] = 2.0 * ts_yz_xxxx[i] * gfe2_0 + 2.0 * gr_yz_xxxx[i] * gfe_0 + ts_yzz_xxxx[i] * gfe_0 * gc_z[i] + gr_yzz_xxxx[i] * gc_z[i];

        grr_z_yzz_xxxy[i] = 2.0 * ts_yz_xxxy[i] * gfe2_0 + 2.0 * gr_yz_xxxy[i] * gfe_0 + ts_yzz_xxxy[i] * gfe_0 * gc_z[i] + gr_yzz_xxxy[i] * gc_z[i];

        grr_z_yzz_xxxz[i] = 2.0 * ts_yz_xxxz[i] * gfe2_0 + 2.0 * gr_yz_xxxz[i] * gfe_0 + ts_yzz_xxx[i] * gfe2_0 + gr_yzz_xxx[i] * gfe_0 + ts_yzz_xxxz[i] * gfe_0 * gc_z[i] + gr_yzz_xxxz[i] * gc_z[i];

        grr_z_yzz_xxyy[i] = 2.0 * ts_yz_xxyy[i] * gfe2_0 + 2.0 * gr_yz_xxyy[i] * gfe_0 + ts_yzz_xxyy[i] * gfe_0 * gc_z[i] + gr_yzz_xxyy[i] * gc_z[i];

        grr_z_yzz_xxyz[i] = 2.0 * ts_yz_xxyz[i] * gfe2_0 + 2.0 * gr_yz_xxyz[i] * gfe_0 + ts_yzz_xxy[i] * gfe2_0 + gr_yzz_xxy[i] * gfe_0 + ts_yzz_xxyz[i] * gfe_0 * gc_z[i] + gr_yzz_xxyz[i] * gc_z[i];

        grr_z_yzz_xxzz[i] = 2.0 * ts_yz_xxzz[i] * gfe2_0 + 2.0 * gr_yz_xxzz[i] * gfe_0 + 2.0 * ts_yzz_xxz[i] * gfe2_0 + 2.0 * gr_yzz_xxz[i] * gfe_0 + ts_yzz_xxzz[i] * gfe_0 * gc_z[i] + gr_yzz_xxzz[i] * gc_z[i];

        grr_z_yzz_xyyy[i] = 2.0 * ts_yz_xyyy[i] * gfe2_0 + 2.0 * gr_yz_xyyy[i] * gfe_0 + ts_yzz_xyyy[i] * gfe_0 * gc_z[i] + gr_yzz_xyyy[i] * gc_z[i];

        grr_z_yzz_xyyz[i] = 2.0 * ts_yz_xyyz[i] * gfe2_0 + 2.0 * gr_yz_xyyz[i] * gfe_0 + ts_yzz_xyy[i] * gfe2_0 + gr_yzz_xyy[i] * gfe_0 + ts_yzz_xyyz[i] * gfe_0 * gc_z[i] + gr_yzz_xyyz[i] * gc_z[i];

        grr_z_yzz_xyzz[i] = 2.0 * ts_yz_xyzz[i] * gfe2_0 + 2.0 * gr_yz_xyzz[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe2_0 + 2.0 * gr_yzz_xyz[i] * gfe_0 + ts_yzz_xyzz[i] * gfe_0 * gc_z[i] + gr_yzz_xyzz[i] * gc_z[i];

        grr_z_yzz_xzzz[i] = 2.0 * ts_yz_xzzz[i] * gfe2_0 + 2.0 * gr_yz_xzzz[i] * gfe_0 + 3.0 * ts_yzz_xzz[i] * gfe2_0 + 3.0 * gr_yzz_xzz[i] * gfe_0 + ts_yzz_xzzz[i] * gfe_0 * gc_z[i] + gr_yzz_xzzz[i] * gc_z[i];

        grr_z_yzz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe2_0 + 2.0 * gr_yz_yyyy[i] * gfe_0 + ts_yzz_yyyy[i] * gfe_0 * gc_z[i] + gr_yzz_yyyy[i] * gc_z[i];

        grr_z_yzz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe2_0 + 2.0 * gr_yz_yyyz[i] * gfe_0 + ts_yzz_yyy[i] * gfe2_0 + gr_yzz_yyy[i] * gfe_0 + ts_yzz_yyyz[i] * gfe_0 * gc_z[i] + gr_yzz_yyyz[i] * gc_z[i];

        grr_z_yzz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe2_0 + 2.0 * gr_yz_yyzz[i] * gfe_0 + 2.0 * ts_yzz_yyz[i] * gfe2_0 + 2.0 * gr_yzz_yyz[i] * gfe_0 + ts_yzz_yyzz[i] * gfe_0 * gc_z[i] + gr_yzz_yyzz[i] * gc_z[i];

        grr_z_yzz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe2_0 + 2.0 * gr_yz_yzzz[i] * gfe_0 + 3.0 * ts_yzz_yzz[i] * gfe2_0 + 3.0 * gr_yzz_yzz[i] * gfe_0 + ts_yzz_yzzz[i] * gfe_0 * gc_z[i] + gr_yzz_yzzz[i] * gc_z[i];

        grr_z_yzz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe2_0 + 2.0 * gr_yz_zzzz[i] * gfe_0 + 4.0 * ts_yzz_zzz[i] * gfe2_0 + 4.0 * gr_yzz_zzz[i] * gfe_0 + ts_yzz_zzzz[i] * gfe_0 * gc_z[i] + gr_yzz_zzzz[i] * gc_z[i];
    }

    // Set up 435-450 components of targeted buffer : FG

    auto grr_z_zzz_xxxx = pbuffer.data(idx_gr_fg + 435);

    auto grr_z_zzz_xxxy = pbuffer.data(idx_gr_fg + 436);

    auto grr_z_zzz_xxxz = pbuffer.data(idx_gr_fg + 437);

    auto grr_z_zzz_xxyy = pbuffer.data(idx_gr_fg + 438);

    auto grr_z_zzz_xxyz = pbuffer.data(idx_gr_fg + 439);

    auto grr_z_zzz_xxzz = pbuffer.data(idx_gr_fg + 440);

    auto grr_z_zzz_xyyy = pbuffer.data(idx_gr_fg + 441);

    auto grr_z_zzz_xyyz = pbuffer.data(idx_gr_fg + 442);

    auto grr_z_zzz_xyzz = pbuffer.data(idx_gr_fg + 443);

    auto grr_z_zzz_xzzz = pbuffer.data(idx_gr_fg + 444);

    auto grr_z_zzz_yyyy = pbuffer.data(idx_gr_fg + 445);

    auto grr_z_zzz_yyyz = pbuffer.data(idx_gr_fg + 446);

    auto grr_z_zzz_yyzz = pbuffer.data(idx_gr_fg + 447);

    auto grr_z_zzz_yzzz = pbuffer.data(idx_gr_fg + 448);

    auto grr_z_zzz_zzzz = pbuffer.data(idx_gr_fg + 449);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxzz, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyzz, gr_zz_xzzz, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyzz, gr_zz_yzzz, gr_zz_zzzz, gr_zzz_xxx, gr_zzz_xxxx, gr_zzz_xxxy, gr_zzz_xxxz, gr_zzz_xxy, gr_zzz_xxyy, gr_zzz_xxyz, gr_zzz_xxz, gr_zzz_xxzz, gr_zzz_xyy, gr_zzz_xyyy, gr_zzz_xyyz, gr_zzz_xyz, gr_zzz_xyzz, gr_zzz_xzz, gr_zzz_xzzz, gr_zzz_yyy, gr_zzz_yyyy, gr_zzz_yyyz, gr_zzz_yyz, gr_zzz_yyzz, gr_zzz_yzz, gr_zzz_yzzz, gr_zzz_zzz, gr_zzz_zzzz, grr_z_zzz_xxxx, grr_z_zzz_xxxy, grr_z_zzz_xxxz, grr_z_zzz_xxyy, grr_z_zzz_xxyz, grr_z_zzz_xxzz, grr_z_zzz_xyyy, grr_z_zzz_xyyz, grr_z_zzz_xyzz, grr_z_zzz_xzzz, grr_z_zzz_yyyy, grr_z_zzz_yyyz, grr_z_zzz_yyzz, grr_z_zzz_yzzz, grr_z_zzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zzz_xxxx[i] = 3.0 * ts_zz_xxxx[i] * gfe2_0 + 3.0 * gr_zz_xxxx[i] * gfe_0 + ts_zzz_xxxx[i] * gfe_0 * gc_z[i] + gr_zzz_xxxx[i] * gc_z[i];

        grr_z_zzz_xxxy[i] = 3.0 * ts_zz_xxxy[i] * gfe2_0 + 3.0 * gr_zz_xxxy[i] * gfe_0 + ts_zzz_xxxy[i] * gfe_0 * gc_z[i] + gr_zzz_xxxy[i] * gc_z[i];

        grr_z_zzz_xxxz[i] = 3.0 * ts_zz_xxxz[i] * gfe2_0 + 3.0 * gr_zz_xxxz[i] * gfe_0 + ts_zzz_xxx[i] * gfe2_0 + gr_zzz_xxx[i] * gfe_0 + ts_zzz_xxxz[i] * gfe_0 * gc_z[i] + gr_zzz_xxxz[i] * gc_z[i];

        grr_z_zzz_xxyy[i] = 3.0 * ts_zz_xxyy[i] * gfe2_0 + 3.0 * gr_zz_xxyy[i] * gfe_0 + ts_zzz_xxyy[i] * gfe_0 * gc_z[i] + gr_zzz_xxyy[i] * gc_z[i];

        grr_z_zzz_xxyz[i] = 3.0 * ts_zz_xxyz[i] * gfe2_0 + 3.0 * gr_zz_xxyz[i] * gfe_0 + ts_zzz_xxy[i] * gfe2_0 + gr_zzz_xxy[i] * gfe_0 + ts_zzz_xxyz[i] * gfe_0 * gc_z[i] + gr_zzz_xxyz[i] * gc_z[i];

        grr_z_zzz_xxzz[i] = 3.0 * ts_zz_xxzz[i] * gfe2_0 + 3.0 * gr_zz_xxzz[i] * gfe_0 + 2.0 * ts_zzz_xxz[i] * gfe2_0 + 2.0 * gr_zzz_xxz[i] * gfe_0 + ts_zzz_xxzz[i] * gfe_0 * gc_z[i] + gr_zzz_xxzz[i] * gc_z[i];

        grr_z_zzz_xyyy[i] = 3.0 * ts_zz_xyyy[i] * gfe2_0 + 3.0 * gr_zz_xyyy[i] * gfe_0 + ts_zzz_xyyy[i] * gfe_0 * gc_z[i] + gr_zzz_xyyy[i] * gc_z[i];

        grr_z_zzz_xyyz[i] = 3.0 * ts_zz_xyyz[i] * gfe2_0 + 3.0 * gr_zz_xyyz[i] * gfe_0 + ts_zzz_xyy[i] * gfe2_0 + gr_zzz_xyy[i] * gfe_0 + ts_zzz_xyyz[i] * gfe_0 * gc_z[i] + gr_zzz_xyyz[i] * gc_z[i];

        grr_z_zzz_xyzz[i] = 3.0 * ts_zz_xyzz[i] * gfe2_0 + 3.0 * gr_zz_xyzz[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe2_0 + 2.0 * gr_zzz_xyz[i] * gfe_0 + ts_zzz_xyzz[i] * gfe_0 * gc_z[i] + gr_zzz_xyzz[i] * gc_z[i];

        grr_z_zzz_xzzz[i] = 3.0 * ts_zz_xzzz[i] * gfe2_0 + 3.0 * gr_zz_xzzz[i] * gfe_0 + 3.0 * ts_zzz_xzz[i] * gfe2_0 + 3.0 * gr_zzz_xzz[i] * gfe_0 + ts_zzz_xzzz[i] * gfe_0 * gc_z[i] + gr_zzz_xzzz[i] * gc_z[i];

        grr_z_zzz_yyyy[i] = 3.0 * ts_zz_yyyy[i] * gfe2_0 + 3.0 * gr_zz_yyyy[i] * gfe_0 + ts_zzz_yyyy[i] * gfe_0 * gc_z[i] + gr_zzz_yyyy[i] * gc_z[i];

        grr_z_zzz_yyyz[i] = 3.0 * ts_zz_yyyz[i] * gfe2_0 + 3.0 * gr_zz_yyyz[i] * gfe_0 + ts_zzz_yyy[i] * gfe2_0 + gr_zzz_yyy[i] * gfe_0 + ts_zzz_yyyz[i] * gfe_0 * gc_z[i] + gr_zzz_yyyz[i] * gc_z[i];

        grr_z_zzz_yyzz[i] = 3.0 * ts_zz_yyzz[i] * gfe2_0 + 3.0 * gr_zz_yyzz[i] * gfe_0 + 2.0 * ts_zzz_yyz[i] * gfe2_0 + 2.0 * gr_zzz_yyz[i] * gfe_0 + ts_zzz_yyzz[i] * gfe_0 * gc_z[i] + gr_zzz_yyzz[i] * gc_z[i];

        grr_z_zzz_yzzz[i] = 3.0 * ts_zz_yzzz[i] * gfe2_0 + 3.0 * gr_zz_yzzz[i] * gfe_0 + 3.0 * ts_zzz_yzz[i] * gfe2_0 + 3.0 * gr_zzz_yzz[i] * gfe_0 + ts_zzz_yzzz[i] * gfe_0 * gc_z[i] + gr_zzz_yzzz[i] * gc_z[i];

        grr_z_zzz_zzzz[i] = 3.0 * ts_zz_zzzz[i] * gfe2_0 + 3.0 * gr_zz_zzzz[i] * gfe_0 + 4.0 * ts_zzz_zzz[i] * gfe2_0 + 4.0 * gr_zzz_zzz[i] * gfe_0 + ts_zzz_zzzz[i] * gfe_0 * gc_z[i] + gr_zzz_zzzz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

