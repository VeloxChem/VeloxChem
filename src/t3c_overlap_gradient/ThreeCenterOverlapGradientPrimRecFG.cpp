#include "ThreeCenterOverlapGradientPrimRecFG.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_fg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_fg,
                              const size_t idx_dg,
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

    auto gs_x_xxx_xxxx = pbuffer.data(idx_g_fg);

    auto gs_x_xxx_xxxy = pbuffer.data(idx_g_fg + 1);

    auto gs_x_xxx_xxxz = pbuffer.data(idx_g_fg + 2);

    auto gs_x_xxx_xxyy = pbuffer.data(idx_g_fg + 3);

    auto gs_x_xxx_xxyz = pbuffer.data(idx_g_fg + 4);

    auto gs_x_xxx_xxzz = pbuffer.data(idx_g_fg + 5);

    auto gs_x_xxx_xyyy = pbuffer.data(idx_g_fg + 6);

    auto gs_x_xxx_xyyz = pbuffer.data(idx_g_fg + 7);

    auto gs_x_xxx_xyzz = pbuffer.data(idx_g_fg + 8);

    auto gs_x_xxx_xzzz = pbuffer.data(idx_g_fg + 9);

    auto gs_x_xxx_yyyy = pbuffer.data(idx_g_fg + 10);

    auto gs_x_xxx_yyyz = pbuffer.data(idx_g_fg + 11);

    auto gs_x_xxx_yyzz = pbuffer.data(idx_g_fg + 12);

    auto gs_x_xxx_yzzz = pbuffer.data(idx_g_fg + 13);

    auto gs_x_xxx_zzzz = pbuffer.data(idx_g_fg + 14);

    #pragma omp simd aligned(gc_x, gs_x_xxx_xxxx, gs_x_xxx_xxxy, gs_x_xxx_xxxz, gs_x_xxx_xxyy, gs_x_xxx_xxyz, gs_x_xxx_xxzz, gs_x_xxx_xyyy, gs_x_xxx_xyyz, gs_x_xxx_xyzz, gs_x_xxx_xzzz, gs_x_xxx_yyyy, gs_x_xxx_yyyz, gs_x_xxx_yyzz, gs_x_xxx_yzzz, gs_x_xxx_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxx_xxxx[i] = 6.0 * ts_xx_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxy[i] = 6.0 * ts_xx_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxz[i] = 6.0 * ts_xx_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxyy[i] = 6.0 * ts_xx_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxyz[i] = 6.0 * ts_xx_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxzz[i] = 6.0 * ts_xx_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyyy[i] = 6.0 * ts_xx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyyz[i] = 6.0 * ts_xx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyzz[i] = 6.0 * ts_xx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xzzz[i] = 6.0 * ts_xx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyyy[i] = 6.0 * ts_xx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyyz[i] = 6.0 * ts_xx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyzz[i] = 6.0 * ts_xx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yzzz[i] = 6.0 * ts_xx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_zzzz[i] = 6.0 * ts_xx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 15-30 components of targeted buffer : FG

    auto gs_x_xxy_xxxx = pbuffer.data(idx_g_fg + 15);

    auto gs_x_xxy_xxxy = pbuffer.data(idx_g_fg + 16);

    auto gs_x_xxy_xxxz = pbuffer.data(idx_g_fg + 17);

    auto gs_x_xxy_xxyy = pbuffer.data(idx_g_fg + 18);

    auto gs_x_xxy_xxyz = pbuffer.data(idx_g_fg + 19);

    auto gs_x_xxy_xxzz = pbuffer.data(idx_g_fg + 20);

    auto gs_x_xxy_xyyy = pbuffer.data(idx_g_fg + 21);

    auto gs_x_xxy_xyyz = pbuffer.data(idx_g_fg + 22);

    auto gs_x_xxy_xyzz = pbuffer.data(idx_g_fg + 23);

    auto gs_x_xxy_xzzz = pbuffer.data(idx_g_fg + 24);

    auto gs_x_xxy_yyyy = pbuffer.data(idx_g_fg + 25);

    auto gs_x_xxy_yyyz = pbuffer.data(idx_g_fg + 26);

    auto gs_x_xxy_yyzz = pbuffer.data(idx_g_fg + 27);

    auto gs_x_xxy_yzzz = pbuffer.data(idx_g_fg + 28);

    auto gs_x_xxy_zzzz = pbuffer.data(idx_g_fg + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxy_xxxx, gs_x_xxy_xxxy, gs_x_xxy_xxxz, gs_x_xxy_xxyy, gs_x_xxy_xxyz, gs_x_xxy_xxzz, gs_x_xxy_xyyy, gs_x_xxy_xyyz, gs_x_xxy_xyzz, gs_x_xxy_xzzz, gs_x_xxy_yyyy, gs_x_xxy_yyyz, gs_x_xxy_yyzz, gs_x_xxy_yzzz, gs_x_xxy_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxy_xxxx[i] = 4.0 * ts_xy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxy[i] = 4.0 * ts_xy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxz[i] = 4.0 * ts_xy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxyy[i] = 4.0 * ts_xy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxyz[i] = 4.0 * ts_xy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxzz[i] = 4.0 * ts_xy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyyy[i] = 4.0 * ts_xy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyyz[i] = 4.0 * ts_xy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyzz[i] = 4.0 * ts_xy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xzzz[i] = 4.0 * ts_xy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyyy[i] = 4.0 * ts_xy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyyz[i] = 4.0 * ts_xy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyzz[i] = 4.0 * ts_xy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yzzz[i] = 4.0 * ts_xy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_zzzz[i] = 4.0 * ts_xy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-45 components of targeted buffer : FG

    auto gs_x_xxz_xxxx = pbuffer.data(idx_g_fg + 30);

    auto gs_x_xxz_xxxy = pbuffer.data(idx_g_fg + 31);

    auto gs_x_xxz_xxxz = pbuffer.data(idx_g_fg + 32);

    auto gs_x_xxz_xxyy = pbuffer.data(idx_g_fg + 33);

    auto gs_x_xxz_xxyz = pbuffer.data(idx_g_fg + 34);

    auto gs_x_xxz_xxzz = pbuffer.data(idx_g_fg + 35);

    auto gs_x_xxz_xyyy = pbuffer.data(idx_g_fg + 36);

    auto gs_x_xxz_xyyz = pbuffer.data(idx_g_fg + 37);

    auto gs_x_xxz_xyzz = pbuffer.data(idx_g_fg + 38);

    auto gs_x_xxz_xzzz = pbuffer.data(idx_g_fg + 39);

    auto gs_x_xxz_yyyy = pbuffer.data(idx_g_fg + 40);

    auto gs_x_xxz_yyyz = pbuffer.data(idx_g_fg + 41);

    auto gs_x_xxz_yyzz = pbuffer.data(idx_g_fg + 42);

    auto gs_x_xxz_yzzz = pbuffer.data(idx_g_fg + 43);

    auto gs_x_xxz_zzzz = pbuffer.data(idx_g_fg + 44);

    #pragma omp simd aligned(gc_x, gs_x_xxz_xxxx, gs_x_xxz_xxxy, gs_x_xxz_xxxz, gs_x_xxz_xxyy, gs_x_xxz_xxyz, gs_x_xxz_xxzz, gs_x_xxz_xyyy, gs_x_xxz_xyyz, gs_x_xxz_xyzz, gs_x_xxz_xzzz, gs_x_xxz_yyyy, gs_x_xxz_yyyz, gs_x_xxz_yyzz, gs_x_xxz_yzzz, gs_x_xxz_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxz_xxxx[i] = 4.0 * ts_xz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxy[i] = 4.0 * ts_xz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxz[i] = 4.0 * ts_xz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxyy[i] = 4.0 * ts_xz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxyz[i] = 4.0 * ts_xz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxzz[i] = 4.0 * ts_xz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyyy[i] = 4.0 * ts_xz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyyz[i] = 4.0 * ts_xz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyzz[i] = 4.0 * ts_xz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xzzz[i] = 4.0 * ts_xz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyyy[i] = 4.0 * ts_xz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyyz[i] = 4.0 * ts_xz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyzz[i] = 4.0 * ts_xz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yzzz[i] = 4.0 * ts_xz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_zzzz[i] = 4.0 * ts_xz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 45-60 components of targeted buffer : FG

    auto gs_x_xyy_xxxx = pbuffer.data(idx_g_fg + 45);

    auto gs_x_xyy_xxxy = pbuffer.data(idx_g_fg + 46);

    auto gs_x_xyy_xxxz = pbuffer.data(idx_g_fg + 47);

    auto gs_x_xyy_xxyy = pbuffer.data(idx_g_fg + 48);

    auto gs_x_xyy_xxyz = pbuffer.data(idx_g_fg + 49);

    auto gs_x_xyy_xxzz = pbuffer.data(idx_g_fg + 50);

    auto gs_x_xyy_xyyy = pbuffer.data(idx_g_fg + 51);

    auto gs_x_xyy_xyyz = pbuffer.data(idx_g_fg + 52);

    auto gs_x_xyy_xyzz = pbuffer.data(idx_g_fg + 53);

    auto gs_x_xyy_xzzz = pbuffer.data(idx_g_fg + 54);

    auto gs_x_xyy_yyyy = pbuffer.data(idx_g_fg + 55);

    auto gs_x_xyy_yyyz = pbuffer.data(idx_g_fg + 56);

    auto gs_x_xyy_yyzz = pbuffer.data(idx_g_fg + 57);

    auto gs_x_xyy_yzzz = pbuffer.data(idx_g_fg + 58);

    auto gs_x_xyy_zzzz = pbuffer.data(idx_g_fg + 59);

    #pragma omp simd aligned(gc_x, gs_x_xyy_xxxx, gs_x_xyy_xxxy, gs_x_xyy_xxxz, gs_x_xyy_xxyy, gs_x_xyy_xxyz, gs_x_xyy_xxzz, gs_x_xyy_xyyy, gs_x_xyy_xyyz, gs_x_xyy_xyzz, gs_x_xyy_xzzz, gs_x_xyy_yyyy, gs_x_xyy_yyyz, gs_x_xyy_yyzz, gs_x_xyy_yzzz, gs_x_xyy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyy_xxxx[i] = 2.0 * ts_yy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxy[i] = 2.0 * ts_yy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxz[i] = 2.0 * ts_yy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxyy[i] = 2.0 * ts_yy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxyz[i] = 2.0 * ts_yy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxzz[i] = 2.0 * ts_yy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyyy[i] = 2.0 * ts_yy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyyz[i] = 2.0 * ts_yy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyzz[i] = 2.0 * ts_yy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xzzz[i] = 2.0 * ts_yy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyyz[i] = 2.0 * ts_yy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyzz[i] = 2.0 * ts_yy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yzzz[i] = 2.0 * ts_yy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_zzzz[i] = 2.0 * ts_yy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-75 components of targeted buffer : FG

    auto gs_x_xyz_xxxx = pbuffer.data(idx_g_fg + 60);

    auto gs_x_xyz_xxxy = pbuffer.data(idx_g_fg + 61);

    auto gs_x_xyz_xxxz = pbuffer.data(idx_g_fg + 62);

    auto gs_x_xyz_xxyy = pbuffer.data(idx_g_fg + 63);

    auto gs_x_xyz_xxyz = pbuffer.data(idx_g_fg + 64);

    auto gs_x_xyz_xxzz = pbuffer.data(idx_g_fg + 65);

    auto gs_x_xyz_xyyy = pbuffer.data(idx_g_fg + 66);

    auto gs_x_xyz_xyyz = pbuffer.data(idx_g_fg + 67);

    auto gs_x_xyz_xyzz = pbuffer.data(idx_g_fg + 68);

    auto gs_x_xyz_xzzz = pbuffer.data(idx_g_fg + 69);

    auto gs_x_xyz_yyyy = pbuffer.data(idx_g_fg + 70);

    auto gs_x_xyz_yyyz = pbuffer.data(idx_g_fg + 71);

    auto gs_x_xyz_yyzz = pbuffer.data(idx_g_fg + 72);

    auto gs_x_xyz_yzzz = pbuffer.data(idx_g_fg + 73);

    auto gs_x_xyz_zzzz = pbuffer.data(idx_g_fg + 74);

    #pragma omp simd aligned(gc_x, gs_x_xyz_xxxx, gs_x_xyz_xxxy, gs_x_xyz_xxxz, gs_x_xyz_xxyy, gs_x_xyz_xxyz, gs_x_xyz_xxzz, gs_x_xyz_xyyy, gs_x_xyz_xyyz, gs_x_xyz_xyzz, gs_x_xyz_xzzz, gs_x_xyz_yyyy, gs_x_xyz_yyyz, gs_x_xyz_yyzz, gs_x_xyz_yzzz, gs_x_xyz_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyz_xxxx[i] = 2.0 * ts_yz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxy[i] = 2.0 * ts_yz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxz[i] = 2.0 * ts_yz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxyy[i] = 2.0 * ts_yz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxyz[i] = 2.0 * ts_yz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxzz[i] = 2.0 * ts_yz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyyy[i] = 2.0 * ts_yz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyyz[i] = 2.0 * ts_yz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyzz[i] = 2.0 * ts_yz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xzzz[i] = 2.0 * ts_yz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 75-90 components of targeted buffer : FG

    auto gs_x_xzz_xxxx = pbuffer.data(idx_g_fg + 75);

    auto gs_x_xzz_xxxy = pbuffer.data(idx_g_fg + 76);

    auto gs_x_xzz_xxxz = pbuffer.data(idx_g_fg + 77);

    auto gs_x_xzz_xxyy = pbuffer.data(idx_g_fg + 78);

    auto gs_x_xzz_xxyz = pbuffer.data(idx_g_fg + 79);

    auto gs_x_xzz_xxzz = pbuffer.data(idx_g_fg + 80);

    auto gs_x_xzz_xyyy = pbuffer.data(idx_g_fg + 81);

    auto gs_x_xzz_xyyz = pbuffer.data(idx_g_fg + 82);

    auto gs_x_xzz_xyzz = pbuffer.data(idx_g_fg + 83);

    auto gs_x_xzz_xzzz = pbuffer.data(idx_g_fg + 84);

    auto gs_x_xzz_yyyy = pbuffer.data(idx_g_fg + 85);

    auto gs_x_xzz_yyyz = pbuffer.data(idx_g_fg + 86);

    auto gs_x_xzz_yyzz = pbuffer.data(idx_g_fg + 87);

    auto gs_x_xzz_yzzz = pbuffer.data(idx_g_fg + 88);

    auto gs_x_xzz_zzzz = pbuffer.data(idx_g_fg + 89);

    #pragma omp simd aligned(gc_x, gs_x_xzz_xxxx, gs_x_xzz_xxxy, gs_x_xzz_xxxz, gs_x_xzz_xxyy, gs_x_xzz_xxyz, gs_x_xzz_xxzz, gs_x_xzz_xyyy, gs_x_xzz_xyyz, gs_x_xzz_xyzz, gs_x_xzz_xzzz, gs_x_xzz_yyyy, gs_x_xzz_yyyz, gs_x_xzz_yyzz, gs_x_xzz_yzzz, gs_x_xzz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxy[i] = 2.0 * ts_zz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxyy[i] = 2.0 * ts_zz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxyz[i] = 2.0 * ts_zz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyyy[i] = 2.0 * ts_zz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyyz[i] = 2.0 * ts_zz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyzz[i] = 2.0 * ts_zz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-105 components of targeted buffer : FG

    auto gs_x_yyy_xxxx = pbuffer.data(idx_g_fg + 90);

    auto gs_x_yyy_xxxy = pbuffer.data(idx_g_fg + 91);

    auto gs_x_yyy_xxxz = pbuffer.data(idx_g_fg + 92);

    auto gs_x_yyy_xxyy = pbuffer.data(idx_g_fg + 93);

    auto gs_x_yyy_xxyz = pbuffer.data(idx_g_fg + 94);

    auto gs_x_yyy_xxzz = pbuffer.data(idx_g_fg + 95);

    auto gs_x_yyy_xyyy = pbuffer.data(idx_g_fg + 96);

    auto gs_x_yyy_xyyz = pbuffer.data(idx_g_fg + 97);

    auto gs_x_yyy_xyzz = pbuffer.data(idx_g_fg + 98);

    auto gs_x_yyy_xzzz = pbuffer.data(idx_g_fg + 99);

    auto gs_x_yyy_yyyy = pbuffer.data(idx_g_fg + 100);

    auto gs_x_yyy_yyyz = pbuffer.data(idx_g_fg + 101);

    auto gs_x_yyy_yyzz = pbuffer.data(idx_g_fg + 102);

    auto gs_x_yyy_yzzz = pbuffer.data(idx_g_fg + 103);

    auto gs_x_yyy_zzzz = pbuffer.data(idx_g_fg + 104);

    #pragma omp simd aligned(gc_x, gs_x_yyy_xxxx, gs_x_yyy_xxxy, gs_x_yyy_xxxz, gs_x_yyy_xxyy, gs_x_yyy_xxyz, gs_x_yyy_xxzz, gs_x_yyy_xyyy, gs_x_yyy_xyyz, gs_x_yyy_xyzz, gs_x_yyy_xzzz, gs_x_yyy_yyyy, gs_x_yyy_yyyz, gs_x_yyy_yyzz, gs_x_yyy_yzzz, gs_x_yyy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyy_xxxx[i] = 8.0 * ts_yyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxy[i] = 6.0 * ts_yyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxz[i] = 6.0 * ts_yyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxyy[i] = 4.0 * ts_yyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxyz[i] = 4.0 * ts_yyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxzz[i] = 4.0 * ts_yyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyyy[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyyz[i] = 2.0 * ts_yyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyzz[i] = 2.0 * ts_yyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xzzz[i] = 2.0 * ts_yyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyyy[i] = 2.0 * ts_yyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyyz[i] = 2.0 * ts_yyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyzz[i] = 2.0 * ts_yyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yzzz[i] = 2.0 * ts_yyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_zzzz[i] = 2.0 * ts_yyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-120 components of targeted buffer : FG

    auto gs_x_yyz_xxxx = pbuffer.data(idx_g_fg + 105);

    auto gs_x_yyz_xxxy = pbuffer.data(idx_g_fg + 106);

    auto gs_x_yyz_xxxz = pbuffer.data(idx_g_fg + 107);

    auto gs_x_yyz_xxyy = pbuffer.data(idx_g_fg + 108);

    auto gs_x_yyz_xxyz = pbuffer.data(idx_g_fg + 109);

    auto gs_x_yyz_xxzz = pbuffer.data(idx_g_fg + 110);

    auto gs_x_yyz_xyyy = pbuffer.data(idx_g_fg + 111);

    auto gs_x_yyz_xyyz = pbuffer.data(idx_g_fg + 112);

    auto gs_x_yyz_xyzz = pbuffer.data(idx_g_fg + 113);

    auto gs_x_yyz_xzzz = pbuffer.data(idx_g_fg + 114);

    auto gs_x_yyz_yyyy = pbuffer.data(idx_g_fg + 115);

    auto gs_x_yyz_yyyz = pbuffer.data(idx_g_fg + 116);

    auto gs_x_yyz_yyzz = pbuffer.data(idx_g_fg + 117);

    auto gs_x_yyz_yzzz = pbuffer.data(idx_g_fg + 118);

    auto gs_x_yyz_zzzz = pbuffer.data(idx_g_fg + 119);

    #pragma omp simd aligned(gc_x, gs_x_yyz_xxxx, gs_x_yyz_xxxy, gs_x_yyz_xxxz, gs_x_yyz_xxyy, gs_x_yyz_xxyz, gs_x_yyz_xxzz, gs_x_yyz_xyyy, gs_x_yyz_xyyz, gs_x_yyz_xyzz, gs_x_yyz_xzzz, gs_x_yyz_yyyy, gs_x_yyz_yyyz, gs_x_yyz_yyzz, gs_x_yyz_yzzz, gs_x_yyz_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyz_xxxx[i] = 8.0 * ts_yyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxy[i] = 6.0 * ts_yyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxz[i] = 6.0 * ts_yyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxyy[i] = 4.0 * ts_yyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxyz[i] = 4.0 * ts_yyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxzz[i] = 4.0 * ts_yyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyyy[i] = 2.0 * ts_yyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyyz[i] = 2.0 * ts_yyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyzz[i] = 2.0 * ts_yyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xzzz[i] = 2.0 * ts_yyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyyy[i] = 2.0 * ts_yyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyyz[i] = 2.0 * ts_yyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyzz[i] = 2.0 * ts_yyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yzzz[i] = 2.0 * ts_yyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_zzzz[i] = 2.0 * ts_yyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 120-135 components of targeted buffer : FG

    auto gs_x_yzz_xxxx = pbuffer.data(idx_g_fg + 120);

    auto gs_x_yzz_xxxy = pbuffer.data(idx_g_fg + 121);

    auto gs_x_yzz_xxxz = pbuffer.data(idx_g_fg + 122);

    auto gs_x_yzz_xxyy = pbuffer.data(idx_g_fg + 123);

    auto gs_x_yzz_xxyz = pbuffer.data(idx_g_fg + 124);

    auto gs_x_yzz_xxzz = pbuffer.data(idx_g_fg + 125);

    auto gs_x_yzz_xyyy = pbuffer.data(idx_g_fg + 126);

    auto gs_x_yzz_xyyz = pbuffer.data(idx_g_fg + 127);

    auto gs_x_yzz_xyzz = pbuffer.data(idx_g_fg + 128);

    auto gs_x_yzz_xzzz = pbuffer.data(idx_g_fg + 129);

    auto gs_x_yzz_yyyy = pbuffer.data(idx_g_fg + 130);

    auto gs_x_yzz_yyyz = pbuffer.data(idx_g_fg + 131);

    auto gs_x_yzz_yyzz = pbuffer.data(idx_g_fg + 132);

    auto gs_x_yzz_yzzz = pbuffer.data(idx_g_fg + 133);

    auto gs_x_yzz_zzzz = pbuffer.data(idx_g_fg + 134);

    #pragma omp simd aligned(gc_x, gs_x_yzz_xxxx, gs_x_yzz_xxxy, gs_x_yzz_xxxz, gs_x_yzz_xxyy, gs_x_yzz_xxyz, gs_x_yzz_xxzz, gs_x_yzz_xyyy, gs_x_yzz_xyyz, gs_x_yzz_xyzz, gs_x_yzz_xzzz, gs_x_yzz_yyyy, gs_x_yzz_yyyz, gs_x_yzz_yyzz, gs_x_yzz_yzzz, gs_x_yzz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzz_xxxx[i] = 8.0 * ts_yzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxy[i] = 6.0 * ts_yzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxz[i] = 6.0 * ts_yzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxyy[i] = 4.0 * ts_yzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxyz[i] = 4.0 * ts_yzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxzz[i] = 4.0 * ts_yzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyyy[i] = 2.0 * ts_yzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyyz[i] = 2.0 * ts_yzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyzz[i] = 2.0 * ts_yzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xzzz[i] = 2.0 * ts_yzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyyy[i] = 2.0 * ts_yzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyyz[i] = 2.0 * ts_yzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyzz[i] = 2.0 * ts_yzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yzzz[i] = 2.0 * ts_yzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_zzzz[i] = 2.0 * ts_yzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 135-150 components of targeted buffer : FG

    auto gs_x_zzz_xxxx = pbuffer.data(idx_g_fg + 135);

    auto gs_x_zzz_xxxy = pbuffer.data(idx_g_fg + 136);

    auto gs_x_zzz_xxxz = pbuffer.data(idx_g_fg + 137);

    auto gs_x_zzz_xxyy = pbuffer.data(idx_g_fg + 138);

    auto gs_x_zzz_xxyz = pbuffer.data(idx_g_fg + 139);

    auto gs_x_zzz_xxzz = pbuffer.data(idx_g_fg + 140);

    auto gs_x_zzz_xyyy = pbuffer.data(idx_g_fg + 141);

    auto gs_x_zzz_xyyz = pbuffer.data(idx_g_fg + 142);

    auto gs_x_zzz_xyzz = pbuffer.data(idx_g_fg + 143);

    auto gs_x_zzz_xzzz = pbuffer.data(idx_g_fg + 144);

    auto gs_x_zzz_yyyy = pbuffer.data(idx_g_fg + 145);

    auto gs_x_zzz_yyyz = pbuffer.data(idx_g_fg + 146);

    auto gs_x_zzz_yyzz = pbuffer.data(idx_g_fg + 147);

    auto gs_x_zzz_yzzz = pbuffer.data(idx_g_fg + 148);

    auto gs_x_zzz_zzzz = pbuffer.data(idx_g_fg + 149);

    #pragma omp simd aligned(gc_x, gs_x_zzz_xxxx, gs_x_zzz_xxxy, gs_x_zzz_xxxz, gs_x_zzz_xxyy, gs_x_zzz_xxyz, gs_x_zzz_xxzz, gs_x_zzz_xyyy, gs_x_zzz_xyyz, gs_x_zzz_xyzz, gs_x_zzz_xzzz, gs_x_zzz_yyyy, gs_x_zzz_yyyz, gs_x_zzz_yyzz, gs_x_zzz_yzzz, gs_x_zzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzz_xxxx[i] = 8.0 * ts_zzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxy[i] = 6.0 * ts_zzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxz[i] = 6.0 * ts_zzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxyy[i] = 4.0 * ts_zzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxyz[i] = 4.0 * ts_zzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxzz[i] = 4.0 * ts_zzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyyy[i] = 2.0 * ts_zzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyyz[i] = 2.0 * ts_zzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyzz[i] = 2.0 * ts_zzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xzzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyyy[i] = 2.0 * ts_zzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyyz[i] = 2.0 * ts_zzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyzz[i] = 2.0 * ts_zzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yzzz[i] = 2.0 * ts_zzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_zzzz[i] = 2.0 * ts_zzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 150-165 components of targeted buffer : FG

    auto gs_y_xxx_xxxx = pbuffer.data(idx_g_fg + 150);

    auto gs_y_xxx_xxxy = pbuffer.data(idx_g_fg + 151);

    auto gs_y_xxx_xxxz = pbuffer.data(idx_g_fg + 152);

    auto gs_y_xxx_xxyy = pbuffer.data(idx_g_fg + 153);

    auto gs_y_xxx_xxyz = pbuffer.data(idx_g_fg + 154);

    auto gs_y_xxx_xxzz = pbuffer.data(idx_g_fg + 155);

    auto gs_y_xxx_xyyy = pbuffer.data(idx_g_fg + 156);

    auto gs_y_xxx_xyyz = pbuffer.data(idx_g_fg + 157);

    auto gs_y_xxx_xyzz = pbuffer.data(idx_g_fg + 158);

    auto gs_y_xxx_xzzz = pbuffer.data(idx_g_fg + 159);

    auto gs_y_xxx_yyyy = pbuffer.data(idx_g_fg + 160);

    auto gs_y_xxx_yyyz = pbuffer.data(idx_g_fg + 161);

    auto gs_y_xxx_yyzz = pbuffer.data(idx_g_fg + 162);

    auto gs_y_xxx_yzzz = pbuffer.data(idx_g_fg + 163);

    auto gs_y_xxx_zzzz = pbuffer.data(idx_g_fg + 164);

    #pragma omp simd aligned(gc_y, gs_y_xxx_xxxx, gs_y_xxx_xxxy, gs_y_xxx_xxxz, gs_y_xxx_xxyy, gs_y_xxx_xxyz, gs_y_xxx_xxzz, gs_y_xxx_xyyy, gs_y_xxx_xyyz, gs_y_xxx_xyzz, gs_y_xxx_xzzz, gs_y_xxx_yyyy, gs_y_xxx_yyyz, gs_y_xxx_yyzz, gs_y_xxx_yzzz, gs_y_xxx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxx_xxxx[i] = 2.0 * ts_xxx_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxy[i] = 2.0 * ts_xxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxz[i] = 2.0 * ts_xxx_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxyy[i] = 4.0 * ts_xxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxyz[i] = 2.0 * ts_xxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxzz[i] = 2.0 * ts_xxx_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyyy[i] = 6.0 * ts_xxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyyz[i] = 4.0 * ts_xxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyzz[i] = 2.0 * ts_xxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xzzz[i] = 2.0 * ts_xxx_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyyy[i] = 8.0 * ts_xxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyyz[i] = 6.0 * ts_xxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyzz[i] = 4.0 * ts_xxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yzzz[i] = 2.0 * ts_xxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_zzzz[i] = 2.0 * ts_xxx_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 165-180 components of targeted buffer : FG

    auto gs_y_xxy_xxxx = pbuffer.data(idx_g_fg + 165);

    auto gs_y_xxy_xxxy = pbuffer.data(idx_g_fg + 166);

    auto gs_y_xxy_xxxz = pbuffer.data(idx_g_fg + 167);

    auto gs_y_xxy_xxyy = pbuffer.data(idx_g_fg + 168);

    auto gs_y_xxy_xxyz = pbuffer.data(idx_g_fg + 169);

    auto gs_y_xxy_xxzz = pbuffer.data(idx_g_fg + 170);

    auto gs_y_xxy_xyyy = pbuffer.data(idx_g_fg + 171);

    auto gs_y_xxy_xyyz = pbuffer.data(idx_g_fg + 172);

    auto gs_y_xxy_xyzz = pbuffer.data(idx_g_fg + 173);

    auto gs_y_xxy_xzzz = pbuffer.data(idx_g_fg + 174);

    auto gs_y_xxy_yyyy = pbuffer.data(idx_g_fg + 175);

    auto gs_y_xxy_yyyz = pbuffer.data(idx_g_fg + 176);

    auto gs_y_xxy_yyzz = pbuffer.data(idx_g_fg + 177);

    auto gs_y_xxy_yzzz = pbuffer.data(idx_g_fg + 178);

    auto gs_y_xxy_zzzz = pbuffer.data(idx_g_fg + 179);

    #pragma omp simd aligned(gc_y, gs_y_xxy_xxxx, gs_y_xxy_xxxy, gs_y_xxy_xxxz, gs_y_xxy_xxyy, gs_y_xxy_xxyz, gs_y_xxy_xxzz, gs_y_xxy_xyyy, gs_y_xxy_xyyz, gs_y_xxy_xyzz, gs_y_xxy_xzzz, gs_y_xxy_yyyy, gs_y_xxy_yyyz, gs_y_xxy_yyzz, gs_y_xxy_yzzz, gs_y_xxy_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxy_xxxx[i] = 2.0 * ts_xx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxy[i] = 2.0 * ts_xx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxz[i] = 2.0 * ts_xx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxyy[i] = 2.0 * ts_xx_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxyz[i] = 2.0 * ts_xx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxzz[i] = 2.0 * ts_xx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyyy[i] = 2.0 * ts_xx_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyyz[i] = 2.0 * ts_xx_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyzz[i] = 2.0 * ts_xx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xzzz[i] = 2.0 * ts_xx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyyy[i] = 2.0 * ts_xx_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyyz[i] = 2.0 * ts_xx_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyzz[i] = 2.0 * ts_xx_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yzzz[i] = 2.0 * ts_xx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_zzzz[i] = 2.0 * ts_xx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 180-195 components of targeted buffer : FG

    auto gs_y_xxz_xxxx = pbuffer.data(idx_g_fg + 180);

    auto gs_y_xxz_xxxy = pbuffer.data(idx_g_fg + 181);

    auto gs_y_xxz_xxxz = pbuffer.data(idx_g_fg + 182);

    auto gs_y_xxz_xxyy = pbuffer.data(idx_g_fg + 183);

    auto gs_y_xxz_xxyz = pbuffer.data(idx_g_fg + 184);

    auto gs_y_xxz_xxzz = pbuffer.data(idx_g_fg + 185);

    auto gs_y_xxz_xyyy = pbuffer.data(idx_g_fg + 186);

    auto gs_y_xxz_xyyz = pbuffer.data(idx_g_fg + 187);

    auto gs_y_xxz_xyzz = pbuffer.data(idx_g_fg + 188);

    auto gs_y_xxz_xzzz = pbuffer.data(idx_g_fg + 189);

    auto gs_y_xxz_yyyy = pbuffer.data(idx_g_fg + 190);

    auto gs_y_xxz_yyyz = pbuffer.data(idx_g_fg + 191);

    auto gs_y_xxz_yyzz = pbuffer.data(idx_g_fg + 192);

    auto gs_y_xxz_yzzz = pbuffer.data(idx_g_fg + 193);

    auto gs_y_xxz_zzzz = pbuffer.data(idx_g_fg + 194);

    #pragma omp simd aligned(gc_y, gs_y_xxz_xxxx, gs_y_xxz_xxxy, gs_y_xxz_xxxz, gs_y_xxz_xxyy, gs_y_xxz_xxyz, gs_y_xxz_xxzz, gs_y_xxz_xyyy, gs_y_xxz_xyyz, gs_y_xxz_xyzz, gs_y_xxz_xzzz, gs_y_xxz_yyyy, gs_y_xxz_yyyz, gs_y_xxz_yyzz, gs_y_xxz_yzzz, gs_y_xxz_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxz_xxxx[i] = 2.0 * ts_xxz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxy[i] = 2.0 * ts_xxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxz[i] = 2.0 * ts_xxz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxyy[i] = 4.0 * ts_xxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxyz[i] = 2.0 * ts_xxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxzz[i] = 2.0 * ts_xxz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyyy[i] = 6.0 * ts_xxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyyz[i] = 4.0 * ts_xxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyzz[i] = 2.0 * ts_xxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xzzz[i] = 2.0 * ts_xxz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyyy[i] = 8.0 * ts_xxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyyz[i] = 6.0 * ts_xxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyzz[i] = 4.0 * ts_xxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yzzz[i] = 2.0 * ts_xxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_zzzz[i] = 2.0 * ts_xxz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 195-210 components of targeted buffer : FG

    auto gs_y_xyy_xxxx = pbuffer.data(idx_g_fg + 195);

    auto gs_y_xyy_xxxy = pbuffer.data(idx_g_fg + 196);

    auto gs_y_xyy_xxxz = pbuffer.data(idx_g_fg + 197);

    auto gs_y_xyy_xxyy = pbuffer.data(idx_g_fg + 198);

    auto gs_y_xyy_xxyz = pbuffer.data(idx_g_fg + 199);

    auto gs_y_xyy_xxzz = pbuffer.data(idx_g_fg + 200);

    auto gs_y_xyy_xyyy = pbuffer.data(idx_g_fg + 201);

    auto gs_y_xyy_xyyz = pbuffer.data(idx_g_fg + 202);

    auto gs_y_xyy_xyzz = pbuffer.data(idx_g_fg + 203);

    auto gs_y_xyy_xzzz = pbuffer.data(idx_g_fg + 204);

    auto gs_y_xyy_yyyy = pbuffer.data(idx_g_fg + 205);

    auto gs_y_xyy_yyyz = pbuffer.data(idx_g_fg + 206);

    auto gs_y_xyy_yyzz = pbuffer.data(idx_g_fg + 207);

    auto gs_y_xyy_yzzz = pbuffer.data(idx_g_fg + 208);

    auto gs_y_xyy_zzzz = pbuffer.data(idx_g_fg + 209);

    #pragma omp simd aligned(gc_y, gs_y_xyy_xxxx, gs_y_xyy_xxxy, gs_y_xyy_xxxz, gs_y_xyy_xxyy, gs_y_xyy_xxyz, gs_y_xyy_xxzz, gs_y_xyy_xyyy, gs_y_xyy_xyyz, gs_y_xyy_xyzz, gs_y_xyy_xzzz, gs_y_xyy_yyyy, gs_y_xyy_yyyz, gs_y_xyy_yyzz, gs_y_xyy_yzzz, gs_y_xyy_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyy_xxxx[i] = 4.0 * ts_xy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxy[i] = 4.0 * ts_xy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxz[i] = 4.0 * ts_xy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxyy[i] = 4.0 * ts_xy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxyz[i] = 4.0 * ts_xy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxzz[i] = 4.0 * ts_xy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyyy[i] = 4.0 * ts_xy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyyz[i] = 4.0 * ts_xy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyzz[i] = 4.0 * ts_xy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xzzz[i] = 4.0 * ts_xy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyyy[i] = 4.0 * ts_xy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyyz[i] = 4.0 * ts_xy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyzz[i] = 4.0 * ts_xy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yzzz[i] = 4.0 * ts_xy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_zzzz[i] = 4.0 * ts_xy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 210-225 components of targeted buffer : FG

    auto gs_y_xyz_xxxx = pbuffer.data(idx_g_fg + 210);

    auto gs_y_xyz_xxxy = pbuffer.data(idx_g_fg + 211);

    auto gs_y_xyz_xxxz = pbuffer.data(idx_g_fg + 212);

    auto gs_y_xyz_xxyy = pbuffer.data(idx_g_fg + 213);

    auto gs_y_xyz_xxyz = pbuffer.data(idx_g_fg + 214);

    auto gs_y_xyz_xxzz = pbuffer.data(idx_g_fg + 215);

    auto gs_y_xyz_xyyy = pbuffer.data(idx_g_fg + 216);

    auto gs_y_xyz_xyyz = pbuffer.data(idx_g_fg + 217);

    auto gs_y_xyz_xyzz = pbuffer.data(idx_g_fg + 218);

    auto gs_y_xyz_xzzz = pbuffer.data(idx_g_fg + 219);

    auto gs_y_xyz_yyyy = pbuffer.data(idx_g_fg + 220);

    auto gs_y_xyz_yyyz = pbuffer.data(idx_g_fg + 221);

    auto gs_y_xyz_yyzz = pbuffer.data(idx_g_fg + 222);

    auto gs_y_xyz_yzzz = pbuffer.data(idx_g_fg + 223);

    auto gs_y_xyz_zzzz = pbuffer.data(idx_g_fg + 224);

    #pragma omp simd aligned(gc_y, gs_y_xyz_xxxx, gs_y_xyz_xxxy, gs_y_xyz_xxxz, gs_y_xyz_xxyy, gs_y_xyz_xxyz, gs_y_xyz_xxzz, gs_y_xyz_xyyy, gs_y_xyz_xyyz, gs_y_xyz_xyzz, gs_y_xyz_xzzz, gs_y_xyz_yyyy, gs_y_xyz_yyyz, gs_y_xyz_yyzz, gs_y_xyz_yzzz, gs_y_xyz_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyz_xxxx[i] = 2.0 * ts_xz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxy[i] = 2.0 * ts_xz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxz[i] = 2.0 * ts_xz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxyy[i] = 2.0 * ts_xz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxyz[i] = 2.0 * ts_xz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxzz[i] = 2.0 * ts_xz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyyy[i] = 2.0 * ts_xz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyyz[i] = 2.0 * ts_xz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyzz[i] = 2.0 * ts_xz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xzzz[i] = 2.0 * ts_xz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyyy[i] = 2.0 * ts_xz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyyz[i] = 2.0 * ts_xz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyzz[i] = 2.0 * ts_xz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yzzz[i] = 2.0 * ts_xz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_zzzz[i] = 2.0 * ts_xz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 225-240 components of targeted buffer : FG

    auto gs_y_xzz_xxxx = pbuffer.data(idx_g_fg + 225);

    auto gs_y_xzz_xxxy = pbuffer.data(idx_g_fg + 226);

    auto gs_y_xzz_xxxz = pbuffer.data(idx_g_fg + 227);

    auto gs_y_xzz_xxyy = pbuffer.data(idx_g_fg + 228);

    auto gs_y_xzz_xxyz = pbuffer.data(idx_g_fg + 229);

    auto gs_y_xzz_xxzz = pbuffer.data(idx_g_fg + 230);

    auto gs_y_xzz_xyyy = pbuffer.data(idx_g_fg + 231);

    auto gs_y_xzz_xyyz = pbuffer.data(idx_g_fg + 232);

    auto gs_y_xzz_xyzz = pbuffer.data(idx_g_fg + 233);

    auto gs_y_xzz_xzzz = pbuffer.data(idx_g_fg + 234);

    auto gs_y_xzz_yyyy = pbuffer.data(idx_g_fg + 235);

    auto gs_y_xzz_yyyz = pbuffer.data(idx_g_fg + 236);

    auto gs_y_xzz_yyzz = pbuffer.data(idx_g_fg + 237);

    auto gs_y_xzz_yzzz = pbuffer.data(idx_g_fg + 238);

    auto gs_y_xzz_zzzz = pbuffer.data(idx_g_fg + 239);

    #pragma omp simd aligned(gc_y, gs_y_xzz_xxxx, gs_y_xzz_xxxy, gs_y_xzz_xxxz, gs_y_xzz_xxyy, gs_y_xzz_xxyz, gs_y_xzz_xxzz, gs_y_xzz_xyyy, gs_y_xzz_xyyz, gs_y_xzz_xyzz, gs_y_xzz_xzzz, gs_y_xzz_yyyy, gs_y_xzz_yyyz, gs_y_xzz_yyzz, gs_y_xzz_yzzz, gs_y_xzz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzz_xxxx[i] = 2.0 * ts_xzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxy[i] = 2.0 * ts_xzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxz[i] = 2.0 * ts_xzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxyy[i] = 4.0 * ts_xzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxyz[i] = 2.0 * ts_xzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxzz[i] = 2.0 * ts_xzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyyy[i] = 6.0 * ts_xzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyyz[i] = 4.0 * ts_xzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyzz[i] = 2.0 * ts_xzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xzzz[i] = 2.0 * ts_xzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyyy[i] = 8.0 * ts_xzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyyz[i] = 6.0 * ts_xzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyzz[i] = 4.0 * ts_xzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yzzz[i] = 2.0 * ts_xzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_zzzz[i] = 2.0 * ts_xzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 240-255 components of targeted buffer : FG

    auto gs_y_yyy_xxxx = pbuffer.data(idx_g_fg + 240);

    auto gs_y_yyy_xxxy = pbuffer.data(idx_g_fg + 241);

    auto gs_y_yyy_xxxz = pbuffer.data(idx_g_fg + 242);

    auto gs_y_yyy_xxyy = pbuffer.data(idx_g_fg + 243);

    auto gs_y_yyy_xxyz = pbuffer.data(idx_g_fg + 244);

    auto gs_y_yyy_xxzz = pbuffer.data(idx_g_fg + 245);

    auto gs_y_yyy_xyyy = pbuffer.data(idx_g_fg + 246);

    auto gs_y_yyy_xyyz = pbuffer.data(idx_g_fg + 247);

    auto gs_y_yyy_xyzz = pbuffer.data(idx_g_fg + 248);

    auto gs_y_yyy_xzzz = pbuffer.data(idx_g_fg + 249);

    auto gs_y_yyy_yyyy = pbuffer.data(idx_g_fg + 250);

    auto gs_y_yyy_yyyz = pbuffer.data(idx_g_fg + 251);

    auto gs_y_yyy_yyzz = pbuffer.data(idx_g_fg + 252);

    auto gs_y_yyy_yzzz = pbuffer.data(idx_g_fg + 253);

    auto gs_y_yyy_zzzz = pbuffer.data(idx_g_fg + 254);

    #pragma omp simd aligned(gc_y, gs_y_yyy_xxxx, gs_y_yyy_xxxy, gs_y_yyy_xxxz, gs_y_yyy_xxyy, gs_y_yyy_xxyz, gs_y_yyy_xxzz, gs_y_yyy_xyyy, gs_y_yyy_xyyz, gs_y_yyy_xyzz, gs_y_yyy_xzzz, gs_y_yyy_yyyy, gs_y_yyy_yyyz, gs_y_yyy_yyzz, gs_y_yyy_yzzz, gs_y_yyy_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyy_xxxx[i] = 6.0 * ts_yy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxy[i] = 6.0 * ts_yy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxz[i] = 6.0 * ts_yy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxyy[i] = 6.0 * ts_yy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxyz[i] = 6.0 * ts_yy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxzz[i] = 6.0 * ts_yy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyyy[i] = 6.0 * ts_yy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyyz[i] = 6.0 * ts_yy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyzz[i] = 6.0 * ts_yy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xzzz[i] = 6.0 * ts_yy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyyy[i] = 6.0 * ts_yy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyyz[i] = 6.0 * ts_yy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyzz[i] = 6.0 * ts_yy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yzzz[i] = 6.0 * ts_yy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_zzzz[i] = 6.0 * ts_yy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 255-270 components of targeted buffer : FG

    auto gs_y_yyz_xxxx = pbuffer.data(idx_g_fg + 255);

    auto gs_y_yyz_xxxy = pbuffer.data(idx_g_fg + 256);

    auto gs_y_yyz_xxxz = pbuffer.data(idx_g_fg + 257);

    auto gs_y_yyz_xxyy = pbuffer.data(idx_g_fg + 258);

    auto gs_y_yyz_xxyz = pbuffer.data(idx_g_fg + 259);

    auto gs_y_yyz_xxzz = pbuffer.data(idx_g_fg + 260);

    auto gs_y_yyz_xyyy = pbuffer.data(idx_g_fg + 261);

    auto gs_y_yyz_xyyz = pbuffer.data(idx_g_fg + 262);

    auto gs_y_yyz_xyzz = pbuffer.data(idx_g_fg + 263);

    auto gs_y_yyz_xzzz = pbuffer.data(idx_g_fg + 264);

    auto gs_y_yyz_yyyy = pbuffer.data(idx_g_fg + 265);

    auto gs_y_yyz_yyyz = pbuffer.data(idx_g_fg + 266);

    auto gs_y_yyz_yyzz = pbuffer.data(idx_g_fg + 267);

    auto gs_y_yyz_yzzz = pbuffer.data(idx_g_fg + 268);

    auto gs_y_yyz_zzzz = pbuffer.data(idx_g_fg + 269);

    #pragma omp simd aligned(gc_y, gs_y_yyz_xxxx, gs_y_yyz_xxxy, gs_y_yyz_xxxz, gs_y_yyz_xxyy, gs_y_yyz_xxyz, gs_y_yyz_xxzz, gs_y_yyz_xyyy, gs_y_yyz_xyyz, gs_y_yyz_xyzz, gs_y_yyz_xzzz, gs_y_yyz_yyyy, gs_y_yyz_yyyz, gs_y_yyz_yyzz, gs_y_yyz_yzzz, gs_y_yyz_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyz_xxxx[i] = 4.0 * ts_yz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxy[i] = 4.0 * ts_yz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxz[i] = 4.0 * ts_yz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxyy[i] = 4.0 * ts_yz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxyz[i] = 4.0 * ts_yz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxzz[i] = 4.0 * ts_yz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyyy[i] = 4.0 * ts_yz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyyz[i] = 4.0 * ts_yz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyzz[i] = 4.0 * ts_yz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xzzz[i] = 4.0 * ts_yz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyyy[i] = 4.0 * ts_yz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyyz[i] = 4.0 * ts_yz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyzz[i] = 4.0 * ts_yz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yzzz[i] = 4.0 * ts_yz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_zzzz[i] = 4.0 * ts_yz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 270-285 components of targeted buffer : FG

    auto gs_y_yzz_xxxx = pbuffer.data(idx_g_fg + 270);

    auto gs_y_yzz_xxxy = pbuffer.data(idx_g_fg + 271);

    auto gs_y_yzz_xxxz = pbuffer.data(idx_g_fg + 272);

    auto gs_y_yzz_xxyy = pbuffer.data(idx_g_fg + 273);

    auto gs_y_yzz_xxyz = pbuffer.data(idx_g_fg + 274);

    auto gs_y_yzz_xxzz = pbuffer.data(idx_g_fg + 275);

    auto gs_y_yzz_xyyy = pbuffer.data(idx_g_fg + 276);

    auto gs_y_yzz_xyyz = pbuffer.data(idx_g_fg + 277);

    auto gs_y_yzz_xyzz = pbuffer.data(idx_g_fg + 278);

    auto gs_y_yzz_xzzz = pbuffer.data(idx_g_fg + 279);

    auto gs_y_yzz_yyyy = pbuffer.data(idx_g_fg + 280);

    auto gs_y_yzz_yyyz = pbuffer.data(idx_g_fg + 281);

    auto gs_y_yzz_yyzz = pbuffer.data(idx_g_fg + 282);

    auto gs_y_yzz_yzzz = pbuffer.data(idx_g_fg + 283);

    auto gs_y_yzz_zzzz = pbuffer.data(idx_g_fg + 284);

    #pragma omp simd aligned(gc_y, gs_y_yzz_xxxx, gs_y_yzz_xxxy, gs_y_yzz_xxxz, gs_y_yzz_xxyy, gs_y_yzz_xxyz, gs_y_yzz_xxzz, gs_y_yzz_xyyy, gs_y_yzz_xyyz, gs_y_yzz_xyzz, gs_y_yzz_xzzz, gs_y_yzz_yyyy, gs_y_yzz_yyyz, gs_y_yzz_yyzz, gs_y_yzz_yzzz, gs_y_yzz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxy[i] = 2.0 * ts_zz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxyy[i] = 2.0 * ts_zz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxyz[i] = 2.0 * ts_zz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyyy[i] = 2.0 * ts_zz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyyz[i] = 2.0 * ts_zz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyzz[i] = 2.0 * ts_zz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 285-300 components of targeted buffer : FG

    auto gs_y_zzz_xxxx = pbuffer.data(idx_g_fg + 285);

    auto gs_y_zzz_xxxy = pbuffer.data(idx_g_fg + 286);

    auto gs_y_zzz_xxxz = pbuffer.data(idx_g_fg + 287);

    auto gs_y_zzz_xxyy = pbuffer.data(idx_g_fg + 288);

    auto gs_y_zzz_xxyz = pbuffer.data(idx_g_fg + 289);

    auto gs_y_zzz_xxzz = pbuffer.data(idx_g_fg + 290);

    auto gs_y_zzz_xyyy = pbuffer.data(idx_g_fg + 291);

    auto gs_y_zzz_xyyz = pbuffer.data(idx_g_fg + 292);

    auto gs_y_zzz_xyzz = pbuffer.data(idx_g_fg + 293);

    auto gs_y_zzz_xzzz = pbuffer.data(idx_g_fg + 294);

    auto gs_y_zzz_yyyy = pbuffer.data(idx_g_fg + 295);

    auto gs_y_zzz_yyyz = pbuffer.data(idx_g_fg + 296);

    auto gs_y_zzz_yyzz = pbuffer.data(idx_g_fg + 297);

    auto gs_y_zzz_yzzz = pbuffer.data(idx_g_fg + 298);

    auto gs_y_zzz_zzzz = pbuffer.data(idx_g_fg + 299);

    #pragma omp simd aligned(gc_y, gs_y_zzz_xxxx, gs_y_zzz_xxxy, gs_y_zzz_xxxz, gs_y_zzz_xxyy, gs_y_zzz_xxyz, gs_y_zzz_xxzz, gs_y_zzz_xyyy, gs_y_zzz_xyyz, gs_y_zzz_xyzz, gs_y_zzz_xzzz, gs_y_zzz_yyyy, gs_y_zzz_yyyz, gs_y_zzz_yyzz, gs_y_zzz_yzzz, gs_y_zzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzz_xxxx[i] = 2.0 * ts_zzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxy[i] = 2.0 * ts_zzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxz[i] = 2.0 * ts_zzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxyy[i] = 4.0 * ts_zzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxyz[i] = 2.0 * ts_zzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxzz[i] = 2.0 * ts_zzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyyy[i] = 6.0 * ts_zzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyyz[i] = 4.0 * ts_zzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyzz[i] = 2.0 * ts_zzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xzzz[i] = 2.0 * ts_zzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyyy[i] = 8.0 * ts_zzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyyz[i] = 6.0 * ts_zzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyzz[i] = 4.0 * ts_zzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yzzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_zzzz[i] = 2.0 * ts_zzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 300-315 components of targeted buffer : FG

    auto gs_z_xxx_xxxx = pbuffer.data(idx_g_fg + 300);

    auto gs_z_xxx_xxxy = pbuffer.data(idx_g_fg + 301);

    auto gs_z_xxx_xxxz = pbuffer.data(idx_g_fg + 302);

    auto gs_z_xxx_xxyy = pbuffer.data(idx_g_fg + 303);

    auto gs_z_xxx_xxyz = pbuffer.data(idx_g_fg + 304);

    auto gs_z_xxx_xxzz = pbuffer.data(idx_g_fg + 305);

    auto gs_z_xxx_xyyy = pbuffer.data(idx_g_fg + 306);

    auto gs_z_xxx_xyyz = pbuffer.data(idx_g_fg + 307);

    auto gs_z_xxx_xyzz = pbuffer.data(idx_g_fg + 308);

    auto gs_z_xxx_xzzz = pbuffer.data(idx_g_fg + 309);

    auto gs_z_xxx_yyyy = pbuffer.data(idx_g_fg + 310);

    auto gs_z_xxx_yyyz = pbuffer.data(idx_g_fg + 311);

    auto gs_z_xxx_yyzz = pbuffer.data(idx_g_fg + 312);

    auto gs_z_xxx_yzzz = pbuffer.data(idx_g_fg + 313);

    auto gs_z_xxx_zzzz = pbuffer.data(idx_g_fg + 314);

    #pragma omp simd aligned(gc_z, gs_z_xxx_xxxx, gs_z_xxx_xxxy, gs_z_xxx_xxxz, gs_z_xxx_xxyy, gs_z_xxx_xxyz, gs_z_xxx_xxzz, gs_z_xxx_xyyy, gs_z_xxx_xyyz, gs_z_xxx_xyzz, gs_z_xxx_xzzz, gs_z_xxx_yyyy, gs_z_xxx_yyyz, gs_z_xxx_yyzz, gs_z_xxx_yzzz, gs_z_xxx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxx_xxxx[i] = 2.0 * ts_xxx_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxy[i] = 2.0 * ts_xxx_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxz[i] = 2.0 * ts_xxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxyy[i] = 2.0 * ts_xxx_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxyz[i] = 2.0 * ts_xxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxzz[i] = 4.0 * ts_xxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyyy[i] = 2.0 * ts_xxx_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyyz[i] = 2.0 * ts_xxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyzz[i] = 4.0 * ts_xxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xzzz[i] = 6.0 * ts_xxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyyy[i] = 2.0 * ts_xxx_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyyz[i] = 2.0 * ts_xxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyzz[i] = 4.0 * ts_xxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yzzz[i] = 6.0 * ts_xxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_zzzz[i] = 8.0 * ts_xxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 315-330 components of targeted buffer : FG

    auto gs_z_xxy_xxxx = pbuffer.data(idx_g_fg + 315);

    auto gs_z_xxy_xxxy = pbuffer.data(idx_g_fg + 316);

    auto gs_z_xxy_xxxz = pbuffer.data(idx_g_fg + 317);

    auto gs_z_xxy_xxyy = pbuffer.data(idx_g_fg + 318);

    auto gs_z_xxy_xxyz = pbuffer.data(idx_g_fg + 319);

    auto gs_z_xxy_xxzz = pbuffer.data(idx_g_fg + 320);

    auto gs_z_xxy_xyyy = pbuffer.data(idx_g_fg + 321);

    auto gs_z_xxy_xyyz = pbuffer.data(idx_g_fg + 322);

    auto gs_z_xxy_xyzz = pbuffer.data(idx_g_fg + 323);

    auto gs_z_xxy_xzzz = pbuffer.data(idx_g_fg + 324);

    auto gs_z_xxy_yyyy = pbuffer.data(idx_g_fg + 325);

    auto gs_z_xxy_yyyz = pbuffer.data(idx_g_fg + 326);

    auto gs_z_xxy_yyzz = pbuffer.data(idx_g_fg + 327);

    auto gs_z_xxy_yzzz = pbuffer.data(idx_g_fg + 328);

    auto gs_z_xxy_zzzz = pbuffer.data(idx_g_fg + 329);

    #pragma omp simd aligned(gc_z, gs_z_xxy_xxxx, gs_z_xxy_xxxy, gs_z_xxy_xxxz, gs_z_xxy_xxyy, gs_z_xxy_xxyz, gs_z_xxy_xxzz, gs_z_xxy_xyyy, gs_z_xxy_xyyz, gs_z_xxy_xyzz, gs_z_xxy_xzzz, gs_z_xxy_yyyy, gs_z_xxy_yyyz, gs_z_xxy_yyzz, gs_z_xxy_yzzz, gs_z_xxy_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxy_xxxx[i] = 2.0 * ts_xxy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxy[i] = 2.0 * ts_xxy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxz[i] = 2.0 * ts_xxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxyy[i] = 2.0 * ts_xxy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxyz[i] = 2.0 * ts_xxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxzz[i] = 4.0 * ts_xxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyyy[i] = 2.0 * ts_xxy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyyz[i] = 2.0 * ts_xxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyzz[i] = 4.0 * ts_xxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xzzz[i] = 6.0 * ts_xxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyyy[i] = 2.0 * ts_xxy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyyz[i] = 2.0 * ts_xxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyzz[i] = 4.0 * ts_xxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yzzz[i] = 6.0 * ts_xxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_zzzz[i] = 8.0 * ts_xxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 330-345 components of targeted buffer : FG

    auto gs_z_xxz_xxxx = pbuffer.data(idx_g_fg + 330);

    auto gs_z_xxz_xxxy = pbuffer.data(idx_g_fg + 331);

    auto gs_z_xxz_xxxz = pbuffer.data(idx_g_fg + 332);

    auto gs_z_xxz_xxyy = pbuffer.data(idx_g_fg + 333);

    auto gs_z_xxz_xxyz = pbuffer.data(idx_g_fg + 334);

    auto gs_z_xxz_xxzz = pbuffer.data(idx_g_fg + 335);

    auto gs_z_xxz_xyyy = pbuffer.data(idx_g_fg + 336);

    auto gs_z_xxz_xyyz = pbuffer.data(idx_g_fg + 337);

    auto gs_z_xxz_xyzz = pbuffer.data(idx_g_fg + 338);

    auto gs_z_xxz_xzzz = pbuffer.data(idx_g_fg + 339);

    auto gs_z_xxz_yyyy = pbuffer.data(idx_g_fg + 340);

    auto gs_z_xxz_yyyz = pbuffer.data(idx_g_fg + 341);

    auto gs_z_xxz_yyzz = pbuffer.data(idx_g_fg + 342);

    auto gs_z_xxz_yzzz = pbuffer.data(idx_g_fg + 343);

    auto gs_z_xxz_zzzz = pbuffer.data(idx_g_fg + 344);

    #pragma omp simd aligned(gc_z, gs_z_xxz_xxxx, gs_z_xxz_xxxy, gs_z_xxz_xxxz, gs_z_xxz_xxyy, gs_z_xxz_xxyz, gs_z_xxz_xxzz, gs_z_xxz_xyyy, gs_z_xxz_xyyz, gs_z_xxz_xyzz, gs_z_xxz_xzzz, gs_z_xxz_yyyy, gs_z_xxz_yyyz, gs_z_xxz_yyzz, gs_z_xxz_yzzz, gs_z_xxz_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxz_xxxx[i] = 2.0 * ts_xx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxy[i] = 2.0 * ts_xx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxz[i] = 2.0 * ts_xx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxyy[i] = 2.0 * ts_xx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxyz[i] = 2.0 * ts_xx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxzz[i] = 2.0 * ts_xx_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyyy[i] = 2.0 * ts_xx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyyz[i] = 2.0 * ts_xx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyzz[i] = 2.0 * ts_xx_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xzzz[i] = 2.0 * ts_xx_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyyy[i] = 2.0 * ts_xx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyyz[i] = 2.0 * ts_xx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyzz[i] = 2.0 * ts_xx_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yzzz[i] = 2.0 * ts_xx_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_zzzz[i] = 2.0 * ts_xx_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 345-360 components of targeted buffer : FG

    auto gs_z_xyy_xxxx = pbuffer.data(idx_g_fg + 345);

    auto gs_z_xyy_xxxy = pbuffer.data(idx_g_fg + 346);

    auto gs_z_xyy_xxxz = pbuffer.data(idx_g_fg + 347);

    auto gs_z_xyy_xxyy = pbuffer.data(idx_g_fg + 348);

    auto gs_z_xyy_xxyz = pbuffer.data(idx_g_fg + 349);

    auto gs_z_xyy_xxzz = pbuffer.data(idx_g_fg + 350);

    auto gs_z_xyy_xyyy = pbuffer.data(idx_g_fg + 351);

    auto gs_z_xyy_xyyz = pbuffer.data(idx_g_fg + 352);

    auto gs_z_xyy_xyzz = pbuffer.data(idx_g_fg + 353);

    auto gs_z_xyy_xzzz = pbuffer.data(idx_g_fg + 354);

    auto gs_z_xyy_yyyy = pbuffer.data(idx_g_fg + 355);

    auto gs_z_xyy_yyyz = pbuffer.data(idx_g_fg + 356);

    auto gs_z_xyy_yyzz = pbuffer.data(idx_g_fg + 357);

    auto gs_z_xyy_yzzz = pbuffer.data(idx_g_fg + 358);

    auto gs_z_xyy_zzzz = pbuffer.data(idx_g_fg + 359);

    #pragma omp simd aligned(gc_z, gs_z_xyy_xxxx, gs_z_xyy_xxxy, gs_z_xyy_xxxz, gs_z_xyy_xxyy, gs_z_xyy_xxyz, gs_z_xyy_xxzz, gs_z_xyy_xyyy, gs_z_xyy_xyyz, gs_z_xyy_xyzz, gs_z_xyy_xzzz, gs_z_xyy_yyyy, gs_z_xyy_yyyz, gs_z_xyy_yyzz, gs_z_xyy_yzzz, gs_z_xyy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyy_xxxx[i] = 2.0 * ts_xyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxy[i] = 2.0 * ts_xyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxz[i] = 2.0 * ts_xyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxyy[i] = 2.0 * ts_xyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxyz[i] = 2.0 * ts_xyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxzz[i] = 4.0 * ts_xyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyyy[i] = 2.0 * ts_xyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyyz[i] = 2.0 * ts_xyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyzz[i] = 4.0 * ts_xyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xzzz[i] = 6.0 * ts_xyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyyy[i] = 2.0 * ts_xyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyyz[i] = 2.0 * ts_xyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyzz[i] = 4.0 * ts_xyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yzzz[i] = 6.0 * ts_xyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_zzzz[i] = 8.0 * ts_xyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 360-375 components of targeted buffer : FG

    auto gs_z_xyz_xxxx = pbuffer.data(idx_g_fg + 360);

    auto gs_z_xyz_xxxy = pbuffer.data(idx_g_fg + 361);

    auto gs_z_xyz_xxxz = pbuffer.data(idx_g_fg + 362);

    auto gs_z_xyz_xxyy = pbuffer.data(idx_g_fg + 363);

    auto gs_z_xyz_xxyz = pbuffer.data(idx_g_fg + 364);

    auto gs_z_xyz_xxzz = pbuffer.data(idx_g_fg + 365);

    auto gs_z_xyz_xyyy = pbuffer.data(idx_g_fg + 366);

    auto gs_z_xyz_xyyz = pbuffer.data(idx_g_fg + 367);

    auto gs_z_xyz_xyzz = pbuffer.data(idx_g_fg + 368);

    auto gs_z_xyz_xzzz = pbuffer.data(idx_g_fg + 369);

    auto gs_z_xyz_yyyy = pbuffer.data(idx_g_fg + 370);

    auto gs_z_xyz_yyyz = pbuffer.data(idx_g_fg + 371);

    auto gs_z_xyz_yyzz = pbuffer.data(idx_g_fg + 372);

    auto gs_z_xyz_yzzz = pbuffer.data(idx_g_fg + 373);

    auto gs_z_xyz_zzzz = pbuffer.data(idx_g_fg + 374);

    #pragma omp simd aligned(gc_z, gs_z_xyz_xxxx, gs_z_xyz_xxxy, gs_z_xyz_xxxz, gs_z_xyz_xxyy, gs_z_xyz_xxyz, gs_z_xyz_xxzz, gs_z_xyz_xyyy, gs_z_xyz_xyyz, gs_z_xyz_xyzz, gs_z_xyz_xzzz, gs_z_xyz_yyyy, gs_z_xyz_yyyz, gs_z_xyz_yyzz, gs_z_xyz_yzzz, gs_z_xyz_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyz_xxxx[i] = 2.0 * ts_xy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxy[i] = 2.0 * ts_xy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxz[i] = 2.0 * ts_xy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxyy[i] = 2.0 * ts_xy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxyz[i] = 2.0 * ts_xy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxzz[i] = 2.0 * ts_xy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyyy[i] = 2.0 * ts_xy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyyz[i] = 2.0 * ts_xy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyzz[i] = 2.0 * ts_xy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xzzz[i] = 2.0 * ts_xy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyyy[i] = 2.0 * ts_xy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyyz[i] = 2.0 * ts_xy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyzz[i] = 2.0 * ts_xy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yzzz[i] = 2.0 * ts_xy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_zzzz[i] = 2.0 * ts_xy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 375-390 components of targeted buffer : FG

    auto gs_z_xzz_xxxx = pbuffer.data(idx_g_fg + 375);

    auto gs_z_xzz_xxxy = pbuffer.data(idx_g_fg + 376);

    auto gs_z_xzz_xxxz = pbuffer.data(idx_g_fg + 377);

    auto gs_z_xzz_xxyy = pbuffer.data(idx_g_fg + 378);

    auto gs_z_xzz_xxyz = pbuffer.data(idx_g_fg + 379);

    auto gs_z_xzz_xxzz = pbuffer.data(idx_g_fg + 380);

    auto gs_z_xzz_xyyy = pbuffer.data(idx_g_fg + 381);

    auto gs_z_xzz_xyyz = pbuffer.data(idx_g_fg + 382);

    auto gs_z_xzz_xyzz = pbuffer.data(idx_g_fg + 383);

    auto gs_z_xzz_xzzz = pbuffer.data(idx_g_fg + 384);

    auto gs_z_xzz_yyyy = pbuffer.data(idx_g_fg + 385);

    auto gs_z_xzz_yyyz = pbuffer.data(idx_g_fg + 386);

    auto gs_z_xzz_yyzz = pbuffer.data(idx_g_fg + 387);

    auto gs_z_xzz_yzzz = pbuffer.data(idx_g_fg + 388);

    auto gs_z_xzz_zzzz = pbuffer.data(idx_g_fg + 389);

    #pragma omp simd aligned(gc_z, gs_z_xzz_xxxx, gs_z_xzz_xxxy, gs_z_xzz_xxxz, gs_z_xzz_xxyy, gs_z_xzz_xxyz, gs_z_xzz_xxzz, gs_z_xzz_xyyy, gs_z_xzz_xyyz, gs_z_xzz_xyzz, gs_z_xzz_xzzz, gs_z_xzz_yyyy, gs_z_xzz_yyyz, gs_z_xzz_yyzz, gs_z_xzz_yzzz, gs_z_xzz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzz_xxxx[i] = 4.0 * ts_xz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxy[i] = 4.0 * ts_xz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxz[i] = 4.0 * ts_xz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxyy[i] = 4.0 * ts_xz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxyz[i] = 4.0 * ts_xz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxzz[i] = 4.0 * ts_xz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyyy[i] = 4.0 * ts_xz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyyz[i] = 4.0 * ts_xz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyzz[i] = 4.0 * ts_xz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xzzz[i] = 4.0 * ts_xz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyyy[i] = 4.0 * ts_xz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyyz[i] = 4.0 * ts_xz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyzz[i] = 4.0 * ts_xz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yzzz[i] = 4.0 * ts_xz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_zzzz[i] = 4.0 * ts_xz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 390-405 components of targeted buffer : FG

    auto gs_z_yyy_xxxx = pbuffer.data(idx_g_fg + 390);

    auto gs_z_yyy_xxxy = pbuffer.data(idx_g_fg + 391);

    auto gs_z_yyy_xxxz = pbuffer.data(idx_g_fg + 392);

    auto gs_z_yyy_xxyy = pbuffer.data(idx_g_fg + 393);

    auto gs_z_yyy_xxyz = pbuffer.data(idx_g_fg + 394);

    auto gs_z_yyy_xxzz = pbuffer.data(idx_g_fg + 395);

    auto gs_z_yyy_xyyy = pbuffer.data(idx_g_fg + 396);

    auto gs_z_yyy_xyyz = pbuffer.data(idx_g_fg + 397);

    auto gs_z_yyy_xyzz = pbuffer.data(idx_g_fg + 398);

    auto gs_z_yyy_xzzz = pbuffer.data(idx_g_fg + 399);

    auto gs_z_yyy_yyyy = pbuffer.data(idx_g_fg + 400);

    auto gs_z_yyy_yyyz = pbuffer.data(idx_g_fg + 401);

    auto gs_z_yyy_yyzz = pbuffer.data(idx_g_fg + 402);

    auto gs_z_yyy_yzzz = pbuffer.data(idx_g_fg + 403);

    auto gs_z_yyy_zzzz = pbuffer.data(idx_g_fg + 404);

    #pragma omp simd aligned(gc_z, gs_z_yyy_xxxx, gs_z_yyy_xxxy, gs_z_yyy_xxxz, gs_z_yyy_xxyy, gs_z_yyy_xxyz, gs_z_yyy_xxzz, gs_z_yyy_xyyy, gs_z_yyy_xyyz, gs_z_yyy_xyzz, gs_z_yyy_xzzz, gs_z_yyy_yyyy, gs_z_yyy_yyyz, gs_z_yyy_yyzz, gs_z_yyy_yzzz, gs_z_yyy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyy_xxxx[i] = 2.0 * ts_yyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxy[i] = 2.0 * ts_yyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxz[i] = 2.0 * ts_yyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxyy[i] = 2.0 * ts_yyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxyz[i] = 2.0 * ts_yyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxzz[i] = 4.0 * ts_yyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyyy[i] = 2.0 * ts_yyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyyz[i] = 2.0 * ts_yyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyzz[i] = 4.0 * ts_yyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xzzz[i] = 6.0 * ts_yyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyyy[i] = 2.0 * ts_yyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyyz[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyzz[i] = 4.0 * ts_yyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yzzz[i] = 6.0 * ts_yyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_zzzz[i] = 8.0 * ts_yyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 405-420 components of targeted buffer : FG

    auto gs_z_yyz_xxxx = pbuffer.data(idx_g_fg + 405);

    auto gs_z_yyz_xxxy = pbuffer.data(idx_g_fg + 406);

    auto gs_z_yyz_xxxz = pbuffer.data(idx_g_fg + 407);

    auto gs_z_yyz_xxyy = pbuffer.data(idx_g_fg + 408);

    auto gs_z_yyz_xxyz = pbuffer.data(idx_g_fg + 409);

    auto gs_z_yyz_xxzz = pbuffer.data(idx_g_fg + 410);

    auto gs_z_yyz_xyyy = pbuffer.data(idx_g_fg + 411);

    auto gs_z_yyz_xyyz = pbuffer.data(idx_g_fg + 412);

    auto gs_z_yyz_xyzz = pbuffer.data(idx_g_fg + 413);

    auto gs_z_yyz_xzzz = pbuffer.data(idx_g_fg + 414);

    auto gs_z_yyz_yyyy = pbuffer.data(idx_g_fg + 415);

    auto gs_z_yyz_yyyz = pbuffer.data(idx_g_fg + 416);

    auto gs_z_yyz_yyzz = pbuffer.data(idx_g_fg + 417);

    auto gs_z_yyz_yzzz = pbuffer.data(idx_g_fg + 418);

    auto gs_z_yyz_zzzz = pbuffer.data(idx_g_fg + 419);

    #pragma omp simd aligned(gc_z, gs_z_yyz_xxxx, gs_z_yyz_xxxy, gs_z_yyz_xxxz, gs_z_yyz_xxyy, gs_z_yyz_xxyz, gs_z_yyz_xxzz, gs_z_yyz_xyyy, gs_z_yyz_xyyz, gs_z_yyz_xyzz, gs_z_yyz_xzzz, gs_z_yyz_yyyy, gs_z_yyz_yyyz, gs_z_yyz_yyzz, gs_z_yyz_yzzz, gs_z_yyz_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyz_xxxx[i] = 2.0 * ts_yy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxy[i] = 2.0 * ts_yy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxz[i] = 2.0 * ts_yy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxyy[i] = 2.0 * ts_yy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxyz[i] = 2.0 * ts_yy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxzz[i] = 2.0 * ts_yy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyyy[i] = 2.0 * ts_yy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyyz[i] = 2.0 * ts_yy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyzz[i] = 2.0 * ts_yy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xzzz[i] = 2.0 * ts_yy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyyz[i] = 2.0 * ts_yy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyzz[i] = 2.0 * ts_yy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yzzz[i] = 2.0 * ts_yy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_zzzz[i] = 2.0 * ts_yy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 420-435 components of targeted buffer : FG

    auto gs_z_yzz_xxxx = pbuffer.data(idx_g_fg + 420);

    auto gs_z_yzz_xxxy = pbuffer.data(idx_g_fg + 421);

    auto gs_z_yzz_xxxz = pbuffer.data(idx_g_fg + 422);

    auto gs_z_yzz_xxyy = pbuffer.data(idx_g_fg + 423);

    auto gs_z_yzz_xxyz = pbuffer.data(idx_g_fg + 424);

    auto gs_z_yzz_xxzz = pbuffer.data(idx_g_fg + 425);

    auto gs_z_yzz_xyyy = pbuffer.data(idx_g_fg + 426);

    auto gs_z_yzz_xyyz = pbuffer.data(idx_g_fg + 427);

    auto gs_z_yzz_xyzz = pbuffer.data(idx_g_fg + 428);

    auto gs_z_yzz_xzzz = pbuffer.data(idx_g_fg + 429);

    auto gs_z_yzz_yyyy = pbuffer.data(idx_g_fg + 430);

    auto gs_z_yzz_yyyz = pbuffer.data(idx_g_fg + 431);

    auto gs_z_yzz_yyzz = pbuffer.data(idx_g_fg + 432);

    auto gs_z_yzz_yzzz = pbuffer.data(idx_g_fg + 433);

    auto gs_z_yzz_zzzz = pbuffer.data(idx_g_fg + 434);

    #pragma omp simd aligned(gc_z, gs_z_yzz_xxxx, gs_z_yzz_xxxy, gs_z_yzz_xxxz, gs_z_yzz_xxyy, gs_z_yzz_xxyz, gs_z_yzz_xxzz, gs_z_yzz_xyyy, gs_z_yzz_xyyz, gs_z_yzz_xyzz, gs_z_yzz_xzzz, gs_z_yzz_yyyy, gs_z_yzz_yyyz, gs_z_yzz_yyzz, gs_z_yzz_yzzz, gs_z_yzz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzz_xxxx[i] = 4.0 * ts_yz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxy[i] = 4.0 * ts_yz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxz[i] = 4.0 * ts_yz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxyy[i] = 4.0 * ts_yz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxyz[i] = 4.0 * ts_yz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxzz[i] = 4.0 * ts_yz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyyy[i] = 4.0 * ts_yz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyyz[i] = 4.0 * ts_yz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyzz[i] = 4.0 * ts_yz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xzzz[i] = 4.0 * ts_yz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyyy[i] = 4.0 * ts_yz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyyz[i] = 4.0 * ts_yz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyzz[i] = 4.0 * ts_yz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yzzz[i] = 4.0 * ts_yz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_zzzz[i] = 4.0 * ts_yz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 435-450 components of targeted buffer : FG

    auto gs_z_zzz_xxxx = pbuffer.data(idx_g_fg + 435);

    auto gs_z_zzz_xxxy = pbuffer.data(idx_g_fg + 436);

    auto gs_z_zzz_xxxz = pbuffer.data(idx_g_fg + 437);

    auto gs_z_zzz_xxyy = pbuffer.data(idx_g_fg + 438);

    auto gs_z_zzz_xxyz = pbuffer.data(idx_g_fg + 439);

    auto gs_z_zzz_xxzz = pbuffer.data(idx_g_fg + 440);

    auto gs_z_zzz_xyyy = pbuffer.data(idx_g_fg + 441);

    auto gs_z_zzz_xyyz = pbuffer.data(idx_g_fg + 442);

    auto gs_z_zzz_xyzz = pbuffer.data(idx_g_fg + 443);

    auto gs_z_zzz_xzzz = pbuffer.data(idx_g_fg + 444);

    auto gs_z_zzz_yyyy = pbuffer.data(idx_g_fg + 445);

    auto gs_z_zzz_yyyz = pbuffer.data(idx_g_fg + 446);

    auto gs_z_zzz_yyzz = pbuffer.data(idx_g_fg + 447);

    auto gs_z_zzz_yzzz = pbuffer.data(idx_g_fg + 448);

    auto gs_z_zzz_zzzz = pbuffer.data(idx_g_fg + 449);

    #pragma omp simd aligned(gc_z, gs_z_zzz_xxxx, gs_z_zzz_xxxy, gs_z_zzz_xxxz, gs_z_zzz_xxyy, gs_z_zzz_xxyz, gs_z_zzz_xxzz, gs_z_zzz_xyyy, gs_z_zzz_xyyz, gs_z_zzz_xyzz, gs_z_zzz_xzzz, gs_z_zzz_yyyy, gs_z_zzz_yyyz, gs_z_zzz_yyzz, gs_z_zzz_yzzz, gs_z_zzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzz_xxxx[i] = 6.0 * ts_zz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxy[i] = 6.0 * ts_zz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxz[i] = 6.0 * ts_zz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxyy[i] = 6.0 * ts_zz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxyz[i] = 6.0 * ts_zz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxzz[i] = 6.0 * ts_zz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyyy[i] = 6.0 * ts_zz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyyz[i] = 6.0 * ts_zz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyzz[i] = 6.0 * ts_zz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xzzz[i] = 6.0 * ts_zz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyyy[i] = 6.0 * ts_zz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyyz[i] = 6.0 * ts_zz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyzz[i] = 6.0 * ts_zz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yzzz[i] = 6.0 * ts_zz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_zzzz[i] = 6.0 * ts_zz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_zzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

