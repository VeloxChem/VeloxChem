#include "ThreeCenterR2PrimRecGG.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_gg(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gg,
                const size_t idx_dg,
                const size_t idx_ff,
                const size_t idx_fg,
                const size_t idx_gd,
                const size_t idx_gf,
                const size_t idx_gg,
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

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_gd + 11);

    auto ts_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto ts_xxyz_xx = pbuffer.data(idx_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_gd + 29);

    auto ts_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto ts_xyyz_xx = pbuffer.data(idx_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto ts_xyzz_xx = pbuffer.data(idx_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_gd + 53);

    auto ts_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto ts_yyyz_xx = pbuffer.data(idx_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_gf + 10);

    auto ts_xxxy_xxy = pbuffer.data(idx_gf + 11);

    auto ts_xxxy_xxz = pbuffer.data(idx_gf + 12);

    auto ts_xxxy_xyy = pbuffer.data(idx_gf + 13);

    auto ts_xxxy_xyz = pbuffer.data(idx_gf + 14);

    auto ts_xxxy_xzz = pbuffer.data(idx_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_gf + 18);

    auto ts_xxxy_zzz = pbuffer.data(idx_gf + 19);

    auto ts_xxxz_xxx = pbuffer.data(idx_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_gf + 23);

    auto ts_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_gf + 25);

    auto ts_xxxz_yyy = pbuffer.data(idx_gf + 26);

    auto ts_xxxz_yyz = pbuffer.data(idx_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_gf + 29);

    auto ts_xxyy_xxx = pbuffer.data(idx_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_gf + 39);

    auto ts_xxyz_xxx = pbuffer.data(idx_gf + 40);

    auto ts_xxyz_xxy = pbuffer.data(idx_gf + 41);

    auto ts_xxyz_xxz = pbuffer.data(idx_gf + 42);

    auto ts_xxyz_xyy = pbuffer.data(idx_gf + 43);

    auto ts_xxyz_xyz = pbuffer.data(idx_gf + 44);

    auto ts_xxyz_xzz = pbuffer.data(idx_gf + 45);

    auto ts_xxyz_yyy = pbuffer.data(idx_gf + 46);

    auto ts_xxyz_yyz = pbuffer.data(idx_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_gf + 48);

    auto ts_xxyz_zzz = pbuffer.data(idx_gf + 49);

    auto ts_xxzz_xxx = pbuffer.data(idx_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_gf + 59);

    auto ts_xyyy_xxx = pbuffer.data(idx_gf + 60);

    auto ts_xyyy_xxy = pbuffer.data(idx_gf + 61);

    auto ts_xyyy_xxz = pbuffer.data(idx_gf + 62);

    auto ts_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto ts_xyyy_xzz = pbuffer.data(idx_gf + 65);

    auto ts_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_gf + 69);

    auto ts_xyyz_xxx = pbuffer.data(idx_gf + 70);

    auto ts_xyyz_xxy = pbuffer.data(idx_gf + 71);

    auto ts_xyyz_xxz = pbuffer.data(idx_gf + 72);

    auto ts_xyyz_xyy = pbuffer.data(idx_gf + 73);

    auto ts_xyyz_xyz = pbuffer.data(idx_gf + 74);

    auto ts_xyyz_xzz = pbuffer.data(idx_gf + 75);

    auto ts_xyyz_yyy = pbuffer.data(idx_gf + 76);

    auto ts_xyyz_yyz = pbuffer.data(idx_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_gf + 79);

    auto ts_xyzz_xxx = pbuffer.data(idx_gf + 80);

    auto ts_xyzz_xxy = pbuffer.data(idx_gf + 81);

    auto ts_xyzz_xxz = pbuffer.data(idx_gf + 82);

    auto ts_xyzz_xyy = pbuffer.data(idx_gf + 83);

    auto ts_xyzz_xyz = pbuffer.data(idx_gf + 84);

    auto ts_xyzz_xzz = pbuffer.data(idx_gf + 85);

    auto ts_xyzz_yyy = pbuffer.data(idx_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_gf + 88);

    auto ts_xyzz_zzz = pbuffer.data(idx_gf + 89);

    auto ts_xzzz_xxx = pbuffer.data(idx_gf + 90);

    auto ts_xzzz_xxy = pbuffer.data(idx_gf + 91);

    auto ts_xzzz_xxz = pbuffer.data(idx_gf + 92);

    auto ts_xzzz_xyy = pbuffer.data(idx_gf + 93);

    auto ts_xzzz_xyz = pbuffer.data(idx_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_gf + 109);

    auto ts_yyyz_xxx = pbuffer.data(idx_gf + 110);

    auto ts_yyyz_xxy = pbuffer.data(idx_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_gf + 113);

    auto ts_yyyz_xyz = pbuffer.data(idx_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_gf + 119);

    auto ts_yyzz_xxx = pbuffer.data(idx_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_gf + 130);

    auto ts_yzzz_xxy = pbuffer.data(idx_gf + 131);

    auto ts_yzzz_xxz = pbuffer.data(idx_gf + 132);

    auto ts_yzzz_xyy = pbuffer.data(idx_gf + 133);

    auto ts_yzzz_xyz = pbuffer.data(idx_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_gf + 149);

    // Set up components of auxiliary buffer : GG

    auto ts_xxxx_xxxx = pbuffer.data(idx_gg);

    auto ts_xxxx_xxxy = pbuffer.data(idx_gg + 1);

    auto ts_xxxx_xxxz = pbuffer.data(idx_gg + 2);

    auto ts_xxxx_xxyy = pbuffer.data(idx_gg + 3);

    auto ts_xxxx_xxyz = pbuffer.data(idx_gg + 4);

    auto ts_xxxx_xxzz = pbuffer.data(idx_gg + 5);

    auto ts_xxxx_xyyy = pbuffer.data(idx_gg + 6);

    auto ts_xxxx_xyyz = pbuffer.data(idx_gg + 7);

    auto ts_xxxx_xyzz = pbuffer.data(idx_gg + 8);

    auto ts_xxxx_xzzz = pbuffer.data(idx_gg + 9);

    auto ts_xxxx_yyyy = pbuffer.data(idx_gg + 10);

    auto ts_xxxx_yyyz = pbuffer.data(idx_gg + 11);

    auto ts_xxxx_yyzz = pbuffer.data(idx_gg + 12);

    auto ts_xxxx_yzzz = pbuffer.data(idx_gg + 13);

    auto ts_xxxx_zzzz = pbuffer.data(idx_gg + 14);

    auto ts_xxxy_xxxx = pbuffer.data(idx_gg + 15);

    auto ts_xxxy_xxxy = pbuffer.data(idx_gg + 16);

    auto ts_xxxy_xxxz = pbuffer.data(idx_gg + 17);

    auto ts_xxxy_xxyy = pbuffer.data(idx_gg + 18);

    auto ts_xxxy_xxyz = pbuffer.data(idx_gg + 19);

    auto ts_xxxy_xxzz = pbuffer.data(idx_gg + 20);

    auto ts_xxxy_xyyy = pbuffer.data(idx_gg + 21);

    auto ts_xxxy_xyyz = pbuffer.data(idx_gg + 22);

    auto ts_xxxy_xyzz = pbuffer.data(idx_gg + 23);

    auto ts_xxxy_xzzz = pbuffer.data(idx_gg + 24);

    auto ts_xxxy_yyyy = pbuffer.data(idx_gg + 25);

    auto ts_xxxy_yyyz = pbuffer.data(idx_gg + 26);

    auto ts_xxxy_yyzz = pbuffer.data(idx_gg + 27);

    auto ts_xxxy_yzzz = pbuffer.data(idx_gg + 28);

    auto ts_xxxy_zzzz = pbuffer.data(idx_gg + 29);

    auto ts_xxxz_xxxx = pbuffer.data(idx_gg + 30);

    auto ts_xxxz_xxxy = pbuffer.data(idx_gg + 31);

    auto ts_xxxz_xxxz = pbuffer.data(idx_gg + 32);

    auto ts_xxxz_xxyy = pbuffer.data(idx_gg + 33);

    auto ts_xxxz_xxyz = pbuffer.data(idx_gg + 34);

    auto ts_xxxz_xxzz = pbuffer.data(idx_gg + 35);

    auto ts_xxxz_xyyy = pbuffer.data(idx_gg + 36);

    auto ts_xxxz_xyyz = pbuffer.data(idx_gg + 37);

    auto ts_xxxz_xyzz = pbuffer.data(idx_gg + 38);

    auto ts_xxxz_xzzz = pbuffer.data(idx_gg + 39);

    auto ts_xxxz_yyyy = pbuffer.data(idx_gg + 40);

    auto ts_xxxz_yyyz = pbuffer.data(idx_gg + 41);

    auto ts_xxxz_yyzz = pbuffer.data(idx_gg + 42);

    auto ts_xxxz_yzzz = pbuffer.data(idx_gg + 43);

    auto ts_xxxz_zzzz = pbuffer.data(idx_gg + 44);

    auto ts_xxyy_xxxx = pbuffer.data(idx_gg + 45);

    auto ts_xxyy_xxxy = pbuffer.data(idx_gg + 46);

    auto ts_xxyy_xxxz = pbuffer.data(idx_gg + 47);

    auto ts_xxyy_xxyy = pbuffer.data(idx_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_gg + 49);

    auto ts_xxyy_xxzz = pbuffer.data(idx_gg + 50);

    auto ts_xxyy_xyyy = pbuffer.data(idx_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_gg + 53);

    auto ts_xxyy_xzzz = pbuffer.data(idx_gg + 54);

    auto ts_xxyy_yyyy = pbuffer.data(idx_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_gg + 58);

    auto ts_xxyy_zzzz = pbuffer.data(idx_gg + 59);

    auto ts_xxyz_xxxx = pbuffer.data(idx_gg + 60);

    auto ts_xxyz_xxxy = pbuffer.data(idx_gg + 61);

    auto ts_xxyz_xxxz = pbuffer.data(idx_gg + 62);

    auto ts_xxyz_xxyy = pbuffer.data(idx_gg + 63);

    auto ts_xxyz_xxyz = pbuffer.data(idx_gg + 64);

    auto ts_xxyz_xxzz = pbuffer.data(idx_gg + 65);

    auto ts_xxyz_xyyy = pbuffer.data(idx_gg + 66);

    auto ts_xxyz_xyyz = pbuffer.data(idx_gg + 67);

    auto ts_xxyz_xyzz = pbuffer.data(idx_gg + 68);

    auto ts_xxyz_xzzz = pbuffer.data(idx_gg + 69);

    auto ts_xxyz_yyyy = pbuffer.data(idx_gg + 70);

    auto ts_xxyz_yyyz = pbuffer.data(idx_gg + 71);

    auto ts_xxyz_yyzz = pbuffer.data(idx_gg + 72);

    auto ts_xxyz_yzzz = pbuffer.data(idx_gg + 73);

    auto ts_xxyz_zzzz = pbuffer.data(idx_gg + 74);

    auto ts_xxzz_xxxx = pbuffer.data(idx_gg + 75);

    auto ts_xxzz_xxxy = pbuffer.data(idx_gg + 76);

    auto ts_xxzz_xxxz = pbuffer.data(idx_gg + 77);

    auto ts_xxzz_xxyy = pbuffer.data(idx_gg + 78);

    auto ts_xxzz_xxyz = pbuffer.data(idx_gg + 79);

    auto ts_xxzz_xxzz = pbuffer.data(idx_gg + 80);

    auto ts_xxzz_xyyy = pbuffer.data(idx_gg + 81);

    auto ts_xxzz_xyyz = pbuffer.data(idx_gg + 82);

    auto ts_xxzz_xyzz = pbuffer.data(idx_gg + 83);

    auto ts_xxzz_xzzz = pbuffer.data(idx_gg + 84);

    auto ts_xxzz_yyyy = pbuffer.data(idx_gg + 85);

    auto ts_xxzz_yyyz = pbuffer.data(idx_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_gg + 89);

    auto ts_xyyy_xxxx = pbuffer.data(idx_gg + 90);

    auto ts_xyyy_xxxy = pbuffer.data(idx_gg + 91);

    auto ts_xyyy_xxxz = pbuffer.data(idx_gg + 92);

    auto ts_xyyy_xxyy = pbuffer.data(idx_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_gg + 94);

    auto ts_xyyy_xxzz = pbuffer.data(idx_gg + 95);

    auto ts_xyyy_xyyy = pbuffer.data(idx_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_gg + 98);

    auto ts_xyyy_xzzz = pbuffer.data(idx_gg + 99);

    auto ts_xyyy_yyyy = pbuffer.data(idx_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_gg + 104);

    auto ts_xyyz_xxxx = pbuffer.data(idx_gg + 105);

    auto ts_xyyz_xxxy = pbuffer.data(idx_gg + 106);

    auto ts_xyyz_xxxz = pbuffer.data(idx_gg + 107);

    auto ts_xyyz_xxyy = pbuffer.data(idx_gg + 108);

    auto ts_xyyz_xxyz = pbuffer.data(idx_gg + 109);

    auto ts_xyyz_xxzz = pbuffer.data(idx_gg + 110);

    auto ts_xyyz_xyyy = pbuffer.data(idx_gg + 111);

    auto ts_xyyz_xyyz = pbuffer.data(idx_gg + 112);

    auto ts_xyyz_xyzz = pbuffer.data(idx_gg + 113);

    auto ts_xyyz_xzzz = pbuffer.data(idx_gg + 114);

    auto ts_xyyz_yyyy = pbuffer.data(idx_gg + 115);

    auto ts_xyyz_yyyz = pbuffer.data(idx_gg + 116);

    auto ts_xyyz_yyzz = pbuffer.data(idx_gg + 117);

    auto ts_xyyz_yzzz = pbuffer.data(idx_gg + 118);

    auto ts_xyyz_zzzz = pbuffer.data(idx_gg + 119);

    auto ts_xyzz_xxxx = pbuffer.data(idx_gg + 120);

    auto ts_xyzz_xxxy = pbuffer.data(idx_gg + 121);

    auto ts_xyzz_xxxz = pbuffer.data(idx_gg + 122);

    auto ts_xyzz_xxyy = pbuffer.data(idx_gg + 123);

    auto ts_xyzz_xxyz = pbuffer.data(idx_gg + 124);

    auto ts_xyzz_xxzz = pbuffer.data(idx_gg + 125);

    auto ts_xyzz_xyyy = pbuffer.data(idx_gg + 126);

    auto ts_xyzz_xyyz = pbuffer.data(idx_gg + 127);

    auto ts_xyzz_xyzz = pbuffer.data(idx_gg + 128);

    auto ts_xyzz_xzzz = pbuffer.data(idx_gg + 129);

    auto ts_xyzz_yyyy = pbuffer.data(idx_gg + 130);

    auto ts_xyzz_yyyz = pbuffer.data(idx_gg + 131);

    auto ts_xyzz_yyzz = pbuffer.data(idx_gg + 132);

    auto ts_xyzz_yzzz = pbuffer.data(idx_gg + 133);

    auto ts_xyzz_zzzz = pbuffer.data(idx_gg + 134);

    auto ts_xzzz_xxxx = pbuffer.data(idx_gg + 135);

    auto ts_xzzz_xxxy = pbuffer.data(idx_gg + 136);

    auto ts_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto ts_xzzz_xxyy = pbuffer.data(idx_gg + 138);

    auto ts_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_gg + 140);

    auto ts_xzzz_xyyy = pbuffer.data(idx_gg + 141);

    auto ts_xzzz_xyyz = pbuffer.data(idx_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_gg + 144);

    auto ts_xzzz_yyyy = pbuffer.data(idx_gg + 145);

    auto ts_xzzz_yyyz = pbuffer.data(idx_gg + 146);

    auto ts_xzzz_yyzz = pbuffer.data(idx_gg + 147);

    auto ts_xzzz_yzzz = pbuffer.data(idx_gg + 148);

    auto ts_xzzz_zzzz = pbuffer.data(idx_gg + 149);

    auto ts_yyyy_xxxx = pbuffer.data(idx_gg + 150);

    auto ts_yyyy_xxxy = pbuffer.data(idx_gg + 151);

    auto ts_yyyy_xxxz = pbuffer.data(idx_gg + 152);

    auto ts_yyyy_xxyy = pbuffer.data(idx_gg + 153);

    auto ts_yyyy_xxyz = pbuffer.data(idx_gg + 154);

    auto ts_yyyy_xxzz = pbuffer.data(idx_gg + 155);

    auto ts_yyyy_xyyy = pbuffer.data(idx_gg + 156);

    auto ts_yyyy_xyyz = pbuffer.data(idx_gg + 157);

    auto ts_yyyy_xyzz = pbuffer.data(idx_gg + 158);

    auto ts_yyyy_xzzz = pbuffer.data(idx_gg + 159);

    auto ts_yyyy_yyyy = pbuffer.data(idx_gg + 160);

    auto ts_yyyy_yyyz = pbuffer.data(idx_gg + 161);

    auto ts_yyyy_yyzz = pbuffer.data(idx_gg + 162);

    auto ts_yyyy_yzzz = pbuffer.data(idx_gg + 163);

    auto ts_yyyy_zzzz = pbuffer.data(idx_gg + 164);

    auto ts_yyyz_xxxx = pbuffer.data(idx_gg + 165);

    auto ts_yyyz_xxxy = pbuffer.data(idx_gg + 166);

    auto ts_yyyz_xxxz = pbuffer.data(idx_gg + 167);

    auto ts_yyyz_xxyy = pbuffer.data(idx_gg + 168);

    auto ts_yyyz_xxyz = pbuffer.data(idx_gg + 169);

    auto ts_yyyz_xxzz = pbuffer.data(idx_gg + 170);

    auto ts_yyyz_xyyy = pbuffer.data(idx_gg + 171);

    auto ts_yyyz_xyyz = pbuffer.data(idx_gg + 172);

    auto ts_yyyz_xyzz = pbuffer.data(idx_gg + 173);

    auto ts_yyyz_xzzz = pbuffer.data(idx_gg + 174);

    auto ts_yyyz_yyyy = pbuffer.data(idx_gg + 175);

    auto ts_yyyz_yyyz = pbuffer.data(idx_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_gg + 179);

    auto ts_yyzz_xxxx = pbuffer.data(idx_gg + 180);

    auto ts_yyzz_xxxy = pbuffer.data(idx_gg + 181);

    auto ts_yyzz_xxxz = pbuffer.data(idx_gg + 182);

    auto ts_yyzz_xxyy = pbuffer.data(idx_gg + 183);

    auto ts_yyzz_xxyz = pbuffer.data(idx_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_gg + 185);

    auto ts_yyzz_xyyy = pbuffer.data(idx_gg + 186);

    auto ts_yyzz_xyyz = pbuffer.data(idx_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_gg + 194);

    auto ts_yzzz_xxxx = pbuffer.data(idx_gg + 195);

    auto ts_yzzz_xxxy = pbuffer.data(idx_gg + 196);

    auto ts_yzzz_xxxz = pbuffer.data(idx_gg + 197);

    auto ts_yzzz_xxyy = pbuffer.data(idx_gg + 198);

    auto ts_yzzz_xxyz = pbuffer.data(idx_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_gg + 200);

    auto ts_yzzz_xyyy = pbuffer.data(idx_gg + 201);

    auto ts_yzzz_xyyz = pbuffer.data(idx_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_gg + 203);

    auto ts_yzzz_xzzz = pbuffer.data(idx_gg + 204);

    auto ts_yzzz_yyyy = pbuffer.data(idx_gg + 205);

    auto ts_yzzz_yyyz = pbuffer.data(idx_gg + 206);

    auto ts_yzzz_yyzz = pbuffer.data(idx_gg + 207);

    auto ts_yzzz_yzzz = pbuffer.data(idx_gg + 208);

    auto ts_yzzz_zzzz = pbuffer.data(idx_gg + 209);

    auto ts_zzzz_xxxx = pbuffer.data(idx_gg + 210);

    auto ts_zzzz_xxxy = pbuffer.data(idx_gg + 211);

    auto ts_zzzz_xxxz = pbuffer.data(idx_gg + 212);

    auto ts_zzzz_xxyy = pbuffer.data(idx_gg + 213);

    auto ts_zzzz_xxyz = pbuffer.data(idx_gg + 214);

    auto ts_zzzz_xxzz = pbuffer.data(idx_gg + 215);

    auto ts_zzzz_xyyy = pbuffer.data(idx_gg + 216);

    auto ts_zzzz_xyyz = pbuffer.data(idx_gg + 217);

    auto ts_zzzz_xyzz = pbuffer.data(idx_gg + 218);

    auto ts_zzzz_xzzz = pbuffer.data(idx_gg + 219);

    auto ts_zzzz_yyyy = pbuffer.data(idx_gg + 220);

    auto ts_zzzz_yyyz = pbuffer.data(idx_gg + 221);

    auto ts_zzzz_yyzz = pbuffer.data(idx_gg + 222);

    auto ts_zzzz_yzzz = pbuffer.data(idx_gg + 223);

    auto ts_zzzz_zzzz = pbuffer.data(idx_gg + 224);

    // Set up 0-15 components of targeted buffer : GG

    auto gr_xxxx_xxxx = pbuffer.data(idx_g_gg);

    auto gr_xxxx_xxxy = pbuffer.data(idx_g_gg + 1);

    auto gr_xxxx_xxxz = pbuffer.data(idx_g_gg + 2);

    auto gr_xxxx_xxyy = pbuffer.data(idx_g_gg + 3);

    auto gr_xxxx_xxyz = pbuffer.data(idx_g_gg + 4);

    auto gr_xxxx_xxzz = pbuffer.data(idx_g_gg + 5);

    auto gr_xxxx_xyyy = pbuffer.data(idx_g_gg + 6);

    auto gr_xxxx_xyyz = pbuffer.data(idx_g_gg + 7);

    auto gr_xxxx_xyzz = pbuffer.data(idx_g_gg + 8);

    auto gr_xxxx_xzzz = pbuffer.data(idx_g_gg + 9);

    auto gr_xxxx_yyyy = pbuffer.data(idx_g_gg + 10);

    auto gr_xxxx_yyyz = pbuffer.data(idx_g_gg + 11);

    auto gr_xxxx_yyzz = pbuffer.data(idx_g_gg + 12);

    auto gr_xxxx_yzzz = pbuffer.data(idx_g_gg + 13);

    auto gr_xxxx_zzzz = pbuffer.data(idx_g_gg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_xxxx, gr_xxxx_xxxy, gr_xxxx_xxxz, gr_xxxx_xxyy, gr_xxxx_xxyz, gr_xxxx_xxzz, gr_xxxx_xyyy, gr_xxxx_xyyz, gr_xxxx_xyzz, gr_xxxx_xzzz, gr_xxxx_yyyy, gr_xxxx_yyyz, gr_xxxx_yyzz, gr_xxxx_yzzz, gr_xxxx_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, ts_xxxx_xx, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xy, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_yy, ts_xxxx_yyy, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyz, ts_xxxx_yyzz, ts_xxxx_yz, ts_xxxx_yzz, ts_xxxx_yzzz, ts_xxxx_zz, ts_xxxx_zzz, ts_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxxx_xxxx[i] = 12.0 * ts_xx_xxxx[i] * gfe2_0 + 32.0 * ts_xxx_xxx[i] * gfe2_0 + 8.0 * ts_xxx_xxxx[i] * gfe_0 * gc_x[i] + 12.0 * ts_xxxx_xx[i] * gfe2_0 + 8.0 * ts_xxxx_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxx_xxxx[i] * gfe_0 + ts_xxxx_xxxx[i] * rgc2_0;

        gr_xxxx_xxxy[i] = 12.0 * ts_xx_xxxy[i] * gfe2_0 + 24.0 * ts_xxx_xxy[i] * gfe2_0 + 8.0 * ts_xxx_xxxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_xy[i] * gfe2_0 + 6.0 * ts_xxxx_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_xxxy[i] * gfe_0 + ts_xxxx_xxxy[i] * rgc2_0;

        gr_xxxx_xxxz[i] = 12.0 * ts_xx_xxxz[i] * gfe2_0 + 24.0 * ts_xxx_xxz[i] * gfe2_0 + 8.0 * ts_xxx_xxxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_xz[i] * gfe2_0 + 6.0 * ts_xxxx_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xxxz[i] * gfe_0 + ts_xxxx_xxxz[i] * rgc2_0;

        gr_xxxx_xxyy[i] = 12.0 * ts_xx_xxyy[i] * gfe2_0 + 16.0 * ts_xxx_xyy[i] * gfe2_0 + 8.0 * ts_xxx_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yy[i] * gfe2_0 + 4.0 * ts_xxxx_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xx[i] * gfe2_0 + 4.0 * ts_xxxx_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_xxyy[i] * gfe_0 + ts_xxxx_xxyy[i] * rgc2_0;

        gr_xxxx_xxyz[i] = 12.0 * ts_xx_xxyz[i] * gfe2_0 + 16.0 * ts_xxx_xyz[i] * gfe2_0 + 8.0 * ts_xxx_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yz[i] * gfe2_0 + 4.0 * ts_xxxx_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xxyz[i] * gfe_0 + ts_xxxx_xxyz[i] * rgc2_0;

        gr_xxxx_xxzz[i] = 12.0 * ts_xx_xxzz[i] * gfe2_0 + 16.0 * ts_xxx_xzz[i] * gfe2_0 + 8.0 * ts_xxx_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_zz[i] * gfe2_0 + 4.0 * ts_xxxx_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xx[i] * gfe2_0 + 4.0 * ts_xxxx_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xxzz[i] * gfe_0 + ts_xxxx_xxzz[i] * rgc2_0;

        gr_xxxx_xyyy[i] = 12.0 * ts_xx_xyyy[i] * gfe2_0 + 8.0 * ts_xxx_yyy[i] * gfe2_0 + 8.0 * ts_xxx_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_xy[i] * gfe2_0 + 6.0 * ts_xxxx_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_xyyy[i] * gfe_0 + ts_xxxx_xyyy[i] * rgc2_0;

        gr_xxxx_xyyz[i] = 12.0 * ts_xx_xyyz[i] * gfe2_0 + 8.0 * ts_xxx_yyz[i] * gfe2_0 + 8.0 * ts_xxx_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xz[i] * gfe2_0 + 4.0 * ts_xxxx_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xyyz[i] * gfe_0 + ts_xxxx_xyyz[i] * rgc2_0;

        gr_xxxx_xyzz[i] = 12.0 * ts_xx_xyzz[i] * gfe2_0 + 8.0 * ts_xxx_yzz[i] * gfe2_0 + 8.0 * ts_xxx_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_xy[i] * gfe2_0 + 4.0 * ts_xxxx_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xyzz[i] * gfe_0 + ts_xxxx_xyzz[i] * rgc2_0;

        gr_xxxx_xzzz[i] = 12.0 * ts_xx_xzzz[i] * gfe2_0 + 8.0 * ts_xxx_zzz[i] * gfe2_0 + 8.0 * ts_xxx_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_xz[i] * gfe2_0 + 6.0 * ts_xxxx_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xzzz[i] * gfe_0 + ts_xxxx_xzzz[i] * rgc2_0;

        gr_xxxx_yyyy[i] = 12.0 * ts_xx_yyyy[i] * gfe2_0 + 8.0 * ts_xxx_yyyy[i] * gfe_0 * gc_x[i] + 12.0 * ts_xxxx_yy[i] * gfe2_0 + 8.0 * ts_xxxx_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_yyyy[i] * gfe_0 + ts_xxxx_yyyy[i] * rgc2_0;

        gr_xxxx_yyyz[i] = 12.0 * ts_xx_yyyz[i] * gfe2_0 + 8.0 * ts_xxx_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_yz[i] * gfe2_0 + 6.0 * ts_xxxx_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_yyyz[i] * gfe_0 + ts_xxxx_yyyz[i] * rgc2_0;

        gr_xxxx_yyzz[i] = 12.0 * ts_xx_yyzz[i] * gfe2_0 + 8.0 * ts_xxx_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_zz[i] * gfe2_0 + 4.0 * ts_xxxx_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_yy[i] * gfe2_0 + 4.0 * ts_xxxx_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_yyzz[i] * gfe_0 + ts_xxxx_yyzz[i] * rgc2_0;

        gr_xxxx_yzzz[i] = 12.0 * ts_xx_yzzz[i] * gfe2_0 + 8.0 * ts_xxx_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxx_yz[i] * gfe2_0 + 6.0 * ts_xxxx_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_yzzz[i] * gfe_0 + ts_xxxx_yzzz[i] * rgc2_0;

        gr_xxxx_zzzz[i] = 12.0 * ts_xx_zzzz[i] * gfe2_0 + 8.0 * ts_xxx_zzzz[i] * gfe_0 * gc_x[i] + 12.0 * ts_xxxx_zz[i] * gfe2_0 + 8.0 * ts_xxxx_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_zzzz[i] * gfe_0 + ts_xxxx_zzzz[i] * rgc2_0;
    }

    // Set up 15-30 components of targeted buffer : GG

    auto gr_xxxy_xxxx = pbuffer.data(idx_g_gg + 15);

    auto gr_xxxy_xxxy = pbuffer.data(idx_g_gg + 16);

    auto gr_xxxy_xxxz = pbuffer.data(idx_g_gg + 17);

    auto gr_xxxy_xxyy = pbuffer.data(idx_g_gg + 18);

    auto gr_xxxy_xxyz = pbuffer.data(idx_g_gg + 19);

    auto gr_xxxy_xxzz = pbuffer.data(idx_g_gg + 20);

    auto gr_xxxy_xyyy = pbuffer.data(idx_g_gg + 21);

    auto gr_xxxy_xyyz = pbuffer.data(idx_g_gg + 22);

    auto gr_xxxy_xyzz = pbuffer.data(idx_g_gg + 23);

    auto gr_xxxy_xzzz = pbuffer.data(idx_g_gg + 24);

    auto gr_xxxy_yyyy = pbuffer.data(idx_g_gg + 25);

    auto gr_xxxy_yyyz = pbuffer.data(idx_g_gg + 26);

    auto gr_xxxy_yyzz = pbuffer.data(idx_g_gg + 27);

    auto gr_xxxy_yzzz = pbuffer.data(idx_g_gg + 28);

    auto gr_xxxy_zzzz = pbuffer.data(idx_g_gg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_xxxx, gr_xxxy_xxxy, gr_xxxy_xxxz, gr_xxxy_xxyy, gr_xxxy_xxyz, gr_xxxy_xxzz, gr_xxxy_xyyy, gr_xxxy_xyyz, gr_xxxy_xyzz, gr_xxxy_xzzz, gr_xxxy_yyyy, gr_xxxy_yyyz, gr_xxxy_yyzz, gr_xxxy_yzzz, gr_xxxy_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, ts_xxxy_xx, ts_xxxy_xxx, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxy, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxz, ts_xxxy_xxzz, ts_xxxy_xy, ts_xxxy_xyy, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyz, ts_xxxy_xyzz, ts_xxxy_xz, ts_xxxy_xzz, ts_xxxy_xzzz, ts_xxxy_yy, ts_xxxy_yyy, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyz, ts_xxxy_yyzz, ts_xxxy_yz, ts_xxxy_yzz, ts_xxxy_yzzz, ts_xxxy_zz, ts_xxxy_zzz, ts_xxxy_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxxy_xxxx[i] = 6.0 * ts_xy_xxxx[i] * gfe2_0 + 24.0 * ts_xxy_xxx[i] * gfe2_0 + 6.0 * ts_xxy_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxxy_xx[i] * gfe2_0 + 8.0 * ts_xxxy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxy_xxxx[i] * gfe_0 + ts_xxxy_xxxx[i] * rgc2_0;

        gr_xxxy_xxxy[i] = 6.0 * ts_xy_xxxy[i] * gfe2_0 + 18.0 * ts_xxy_xxy[i] * gfe2_0 + 6.0 * ts_xxy_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxx[i] * gfe2_0 + 2.0 * ts_xxx_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_xy[i] * gfe2_0 + 6.0 * ts_xxxy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_xxxy[i] * gfe_0 + ts_xxxy_xxxy[i] * rgc2_0;

        gr_xxxy_xxxz[i] = 6.0 * ts_xy_xxxz[i] * gfe2_0 + 18.0 * ts_xxy_xxz[i] * gfe2_0 + 6.0 * ts_xxy_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_xz[i] * gfe2_0 + 6.0 * ts_xxxy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xxxz[i] * gfe_0 + ts_xxxy_xxxz[i] * rgc2_0;

        gr_xxxy_xxyy[i] = 6.0 * ts_xy_xxyy[i] * gfe2_0 + 12.0 * ts_xxy_xyy[i] * gfe2_0 + 6.0 * ts_xxy_xxyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_xxy[i] * gfe2_0 + 2.0 * ts_xxx_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yy[i] * gfe2_0 + 4.0 * ts_xxxy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xx[i] * gfe2_0 + 4.0 * ts_xxxy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_xxyy[i] * gfe_0 + ts_xxxy_xxyy[i] * rgc2_0;

        gr_xxxy_xxyz[i] = 6.0 * ts_xy_xxyz[i] * gfe2_0 + 12.0 * ts_xxy_xyz[i] * gfe2_0 + 6.0 * ts_xxy_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxz[i] * gfe2_0 + 2.0 * ts_xxx_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yz[i] * gfe2_0 + 4.0 * ts_xxxy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xxyz[i] * gfe_0 + ts_xxxy_xxyz[i] * rgc2_0;

        gr_xxxy_xxzz[i] = 6.0 * ts_xy_xxzz[i] * gfe2_0 + 12.0 * ts_xxy_xzz[i] * gfe2_0 + 6.0 * ts_xxy_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_zz[i] * gfe2_0 + 4.0 * ts_xxxy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xx[i] * gfe2_0 + 4.0 * ts_xxxy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xxzz[i] * gfe_0 + ts_xxxy_xxzz[i] * rgc2_0;

        gr_xxxy_xyyy[i] = 6.0 * ts_xy_xyyy[i] * gfe2_0 + 6.0 * ts_xxy_yyy[i] * gfe2_0 + 6.0 * ts_xxy_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_xyy[i] * gfe2_0 + 2.0 * ts_xxx_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxy_xy[i] * gfe2_0 + 6.0 * ts_xxxy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_xyyy[i] * gfe_0 + ts_xxxy_xyyy[i] * rgc2_0;

        gr_xxxy_xyyz[i] = 6.0 * ts_xy_xyyz[i] * gfe2_0 + 6.0 * ts_xxy_yyz[i] * gfe2_0 + 6.0 * ts_xxy_xyyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_xyz[i] * gfe2_0 + 2.0 * ts_xxx_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xz[i] * gfe2_0 + 4.0 * ts_xxxy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xyyz[i] * gfe_0 + ts_xxxy_xyyz[i] * rgc2_0;

        gr_xxxy_xyzz[i] = 6.0 * ts_xy_xyzz[i] * gfe2_0 + 6.0 * ts_xxy_yzz[i] * gfe2_0 + 6.0 * ts_xxy_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xzz[i] * gfe2_0 + 2.0 * ts_xxx_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_xy[i] * gfe2_0 + 4.0 * ts_xxxy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xyzz[i] * gfe_0 + ts_xxxy_xyzz[i] * rgc2_0;

        gr_xxxy_xzzz[i] = 6.0 * ts_xy_xzzz[i] * gfe2_0 + 6.0 * ts_xxy_zzz[i] * gfe2_0 + 6.0 * ts_xxy_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxy_xz[i] * gfe2_0 + 6.0 * ts_xxxy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xzzz[i] * gfe_0 + ts_xxxy_xzzz[i] * rgc2_0;

        gr_xxxy_yyyy[i] = 6.0 * ts_xy_yyyy[i] * gfe2_0 + 6.0 * ts_xxy_yyyy[i] * gfe_0 * gc_x[i] + 8.0 * ts_xxx_yyy[i] * gfe2_0 + 2.0 * ts_xxx_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxxy_yy[i] * gfe2_0 + 8.0 * ts_xxxy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_yyyy[i] * gfe_0 + ts_xxxy_yyyy[i] * rgc2_0;

        gr_xxxy_yyyz[i] = 6.0 * ts_xy_yyyz[i] * gfe2_0 + 6.0 * ts_xxy_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_yyz[i] * gfe2_0 + 2.0 * ts_xxx_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_yz[i] * gfe2_0 + 6.0 * ts_xxxy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_yyyz[i] * gfe_0 + ts_xxxy_yyyz[i] * rgc2_0;

        gr_xxxy_yyzz[i] = 6.0 * ts_xy_yyzz[i] * gfe2_0 + 6.0 * ts_xxy_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_yzz[i] * gfe2_0 + 2.0 * ts_xxx_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_zz[i] * gfe2_0 + 4.0 * ts_xxxy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yy[i] * gfe2_0 + 4.0 * ts_xxxy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_yyzz[i] * gfe_0 + ts_xxxy_yyzz[i] * rgc2_0;

        gr_xxxy_yzzz[i] = 6.0 * ts_xy_yzzz[i] * gfe2_0 + 6.0 * ts_xxy_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zzz[i] * gfe2_0 + 2.0 * ts_xxx_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_yz[i] * gfe2_0 + 6.0 * ts_xxxy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_yzzz[i] * gfe_0 + ts_xxxy_yzzz[i] * rgc2_0;

        gr_xxxy_zzzz[i] = 6.0 * ts_xy_zzzz[i] * gfe2_0 + 6.0 * ts_xxy_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxxy_zz[i] * gfe2_0 + 8.0 * ts_xxxy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_zzzz[i] * gfe_0 + ts_xxxy_zzzz[i] * rgc2_0;
    }

    // Set up 30-45 components of targeted buffer : GG

    auto gr_xxxz_xxxx = pbuffer.data(idx_g_gg + 30);

    auto gr_xxxz_xxxy = pbuffer.data(idx_g_gg + 31);

    auto gr_xxxz_xxxz = pbuffer.data(idx_g_gg + 32);

    auto gr_xxxz_xxyy = pbuffer.data(idx_g_gg + 33);

    auto gr_xxxz_xxyz = pbuffer.data(idx_g_gg + 34);

    auto gr_xxxz_xxzz = pbuffer.data(idx_g_gg + 35);

    auto gr_xxxz_xyyy = pbuffer.data(idx_g_gg + 36);

    auto gr_xxxz_xyyz = pbuffer.data(idx_g_gg + 37);

    auto gr_xxxz_xyzz = pbuffer.data(idx_g_gg + 38);

    auto gr_xxxz_xzzz = pbuffer.data(idx_g_gg + 39);

    auto gr_xxxz_yyyy = pbuffer.data(idx_g_gg + 40);

    auto gr_xxxz_yyyz = pbuffer.data(idx_g_gg + 41);

    auto gr_xxxz_yyzz = pbuffer.data(idx_g_gg + 42);

    auto gr_xxxz_yzzz = pbuffer.data(idx_g_gg + 43);

    auto gr_xxxz_zzzz = pbuffer.data(idx_g_gg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_xxxx, gr_xxxz_xxxy, gr_xxxz_xxxz, gr_xxxz_xxyy, gr_xxxz_xxyz, gr_xxxz_xxzz, gr_xxxz_xyyy, gr_xxxz_xyyz, gr_xxxz_xyzz, gr_xxxz_xzzz, gr_xxxz_yyyy, gr_xxxz_yyyz, gr_xxxz_yyzz, gr_xxxz_yzzz, gr_xxxz_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, ts_xxxz_xx, ts_xxxz_xxx, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxy, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxz, ts_xxxz_xxzz, ts_xxxz_xy, ts_xxxz_xyy, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyz, ts_xxxz_xyzz, ts_xxxz_xz, ts_xxxz_xzz, ts_xxxz_xzzz, ts_xxxz_yy, ts_xxxz_yyy, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyz, ts_xxxz_yyzz, ts_xxxz_yz, ts_xxxz_yzz, ts_xxxz_yzzz, ts_xxxz_zz, ts_xxxz_zzz, ts_xxxz_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxxz_xxxx[i] = 6.0 * ts_xz_xxxx[i] * gfe2_0 + 24.0 * ts_xxz_xxx[i] * gfe2_0 + 6.0 * ts_xxz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxxz_xx[i] * gfe2_0 + 8.0 * ts_xxxz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxz_xxxx[i] * gfe_0 + ts_xxxz_xxxx[i] * rgc2_0;

        gr_xxxz_xxxy[i] = 6.0 * ts_xz_xxxy[i] * gfe2_0 + 18.0 * ts_xxz_xxy[i] * gfe2_0 + 6.0 * ts_xxz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxxz_xy[i] * gfe2_0 + 6.0 * ts_xxxz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_xxxy[i] * gfe_0 + ts_xxxz_xxxy[i] * rgc2_0;

        gr_xxxz_xxxz[i] = 6.0 * ts_xz_xxxz[i] * gfe2_0 + 18.0 * ts_xxz_xxz[i] * gfe2_0 + 6.0 * ts_xxz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxx[i] * gfe2_0 + 2.0 * ts_xxx_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxxz_xz[i] * gfe2_0 + 6.0 * ts_xxxz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xxxz[i] * gfe_0 + ts_xxxz_xxxz[i] * rgc2_0;

        gr_xxxz_xxyy[i] = 6.0 * ts_xz_xxyy[i] * gfe2_0 + 12.0 * ts_xxz_xyy[i] * gfe2_0 + 6.0 * ts_xxz_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yy[i] * gfe2_0 + 4.0 * ts_xxxz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xx[i] * gfe2_0 + 4.0 * ts_xxxz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_xxyy[i] * gfe_0 + ts_xxxz_xxyy[i] * rgc2_0;

        gr_xxxz_xxyz[i] = 6.0 * ts_xz_xxyz[i] * gfe2_0 + 12.0 * ts_xxz_xyz[i] * gfe2_0 + 6.0 * ts_xxz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxy[i] * gfe2_0 + 2.0 * ts_xxx_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yz[i] * gfe2_0 + 4.0 * ts_xxxz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xxyz[i] * gfe_0 + ts_xxxz_xxyz[i] * rgc2_0;

        gr_xxxz_xxzz[i] = 6.0 * ts_xz_xxzz[i] * gfe2_0 + 12.0 * ts_xxz_xzz[i] * gfe2_0 + 6.0 * ts_xxz_xxzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_xxz[i] * gfe2_0 + 2.0 * ts_xxx_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_zz[i] * gfe2_0 + 4.0 * ts_xxxz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xx[i] * gfe2_0 + 4.0 * ts_xxxz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xxzz[i] * gfe_0 + ts_xxxz_xxzz[i] * rgc2_0;

        gr_xxxz_xyyy[i] = 6.0 * ts_xz_xyyy[i] * gfe2_0 + 6.0 * ts_xxz_yyy[i] * gfe2_0 + 6.0 * ts_xxz_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxz_xy[i] * gfe2_0 + 6.0 * ts_xxxz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_xyyy[i] * gfe_0 + ts_xxxz_xyyy[i] * rgc2_0;

        gr_xxxz_xyyz[i] = 6.0 * ts_xz_xyyz[i] * gfe2_0 + 6.0 * ts_xxz_yyz[i] * gfe2_0 + 6.0 * ts_xxz_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xyy[i] * gfe2_0 + 2.0 * ts_xxx_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xz[i] * gfe2_0 + 4.0 * ts_xxxz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xyyz[i] * gfe_0 + ts_xxxz_xyyz[i] * rgc2_0;

        gr_xxxz_xyzz[i] = 6.0 * ts_xz_xyzz[i] * gfe2_0 + 6.0 * ts_xxz_yzz[i] * gfe2_0 + 6.0 * ts_xxz_xyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_xyz[i] * gfe2_0 + 2.0 * ts_xxx_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_xy[i] * gfe2_0 + 4.0 * ts_xxxz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xyzz[i] * gfe_0 + ts_xxxz_xyzz[i] * rgc2_0;

        gr_xxxz_xzzz[i] = 6.0 * ts_xz_xzzz[i] * gfe2_0 + 6.0 * ts_xxz_zzz[i] * gfe2_0 + 6.0 * ts_xxz_xzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_xzz[i] * gfe2_0 + 2.0 * ts_xxx_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxz_xz[i] * gfe2_0 + 6.0 * ts_xxxz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xzzz[i] * gfe_0 + ts_xxxz_xzzz[i] * rgc2_0;

        gr_xxxz_yyyy[i] = 6.0 * ts_xz_yyyy[i] * gfe2_0 + 6.0 * ts_xxz_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxxz_yy[i] * gfe2_0 + 8.0 * ts_xxxz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_yyyy[i] * gfe_0 + ts_xxxz_yyyy[i] * rgc2_0;

        gr_xxxz_yyyz[i] = 6.0 * ts_xz_yyyz[i] * gfe2_0 + 6.0 * ts_xxz_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yyy[i] * gfe2_0 + 2.0 * ts_xxx_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxxz_yz[i] * gfe2_0 + 6.0 * ts_xxxz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_yyyz[i] * gfe_0 + ts_xxxz_yyyz[i] * rgc2_0;

        gr_xxxz_yyzz[i] = 6.0 * ts_xz_yyzz[i] * gfe2_0 + 6.0 * ts_xxz_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_yyz[i] * gfe2_0 + 2.0 * ts_xxx_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_zz[i] * gfe2_0 + 4.0 * ts_xxxz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_yy[i] * gfe2_0 + 4.0 * ts_xxxz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_yyzz[i] * gfe_0 + ts_xxxz_yyzz[i] * rgc2_0;

        gr_xxxz_yzzz[i] = 6.0 * ts_xz_yzzz[i] * gfe2_0 + 6.0 * ts_xxz_yzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_yzz[i] * gfe2_0 + 2.0 * ts_xxx_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxz_yz[i] * gfe2_0 + 6.0 * ts_xxxz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_yzzz[i] * gfe_0 + ts_xxxz_yzzz[i] * rgc2_0;

        gr_xxxz_zzzz[i] = 6.0 * ts_xz_zzzz[i] * gfe2_0 + 6.0 * ts_xxz_zzzz[i] * gfe_0 * gc_x[i] + 8.0 * ts_xxx_zzz[i] * gfe2_0 + 2.0 * ts_xxx_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxxz_zz[i] * gfe2_0 + 8.0 * ts_xxxz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_zzzz[i] * gfe_0 + ts_xxxz_zzzz[i] * rgc2_0;
    }

    // Set up 45-60 components of targeted buffer : GG

    auto gr_xxyy_xxxx = pbuffer.data(idx_g_gg + 45);

    auto gr_xxyy_xxxy = pbuffer.data(idx_g_gg + 46);

    auto gr_xxyy_xxxz = pbuffer.data(idx_g_gg + 47);

    auto gr_xxyy_xxyy = pbuffer.data(idx_g_gg + 48);

    auto gr_xxyy_xxyz = pbuffer.data(idx_g_gg + 49);

    auto gr_xxyy_xxzz = pbuffer.data(idx_g_gg + 50);

    auto gr_xxyy_xyyy = pbuffer.data(idx_g_gg + 51);

    auto gr_xxyy_xyyz = pbuffer.data(idx_g_gg + 52);

    auto gr_xxyy_xyzz = pbuffer.data(idx_g_gg + 53);

    auto gr_xxyy_xzzz = pbuffer.data(idx_g_gg + 54);

    auto gr_xxyy_yyyy = pbuffer.data(idx_g_gg + 55);

    auto gr_xxyy_yyyz = pbuffer.data(idx_g_gg + 56);

    auto gr_xxyy_yyzz = pbuffer.data(idx_g_gg + 57);

    auto gr_xxyy_yzzz = pbuffer.data(idx_g_gg + 58);

    auto gr_xxyy_zzzz = pbuffer.data(idx_g_gg + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_xxxx, gr_xxyy_xxxy, gr_xxyy_xxxz, gr_xxyy_xxyy, gr_xxyy_xxyz, gr_xxyy_xxzz, gr_xxyy_xyyy, gr_xxyy_xyyz, gr_xxyy_xyzz, gr_xxyy_xzzz, gr_xxyy_yyyy, gr_xxyy_yyyz, gr_xxyy_yyzz, gr_xxyy_yzzz, gr_xxyy_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, ts_xxyy_xx, ts_xxyy_xxx, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxz, ts_xxyy_xxzz, ts_xxyy_xy, ts_xxyy_xyy, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyz, ts_xxyy_xyzz, ts_xxyy_xz, ts_xxyy_xzz, ts_xxyy_xzzz, ts_xxyy_yy, ts_xxyy_yyy, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyz, ts_xxyy_yyzz, ts_xxyy_yz, ts_xxyy_yzz, ts_xxyy_yzzz, ts_xxyy_zz, ts_xxyy_zzz, ts_xxyy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxyy_xxxx[i] = 2.0 * ts_yy_xxxx[i] * gfe2_0 + 16.0 * ts_xyy_xxx[i] * gfe2_0 + 4.0 * ts_xyy_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxx[i] * gfe2_0 + 4.0 * ts_xxy_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxyy_xx[i] * gfe2_0 + 8.0 * ts_xxyy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyy_xxxx[i] * gfe_0 + ts_xxyy_xxxx[i] * rgc2_0;

        gr_xxyy_xxxy[i] = 2.0 * ts_yy_xxxy[i] * gfe2_0 + 12.0 * ts_xyy_xxy[i] * gfe2_0 + 4.0 * ts_xyy_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxy[i] * gfe2_0 + 4.0 * ts_xxy_xxx[i] * gfe2_0 + 4.0 * ts_xxy_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_xy[i] * gfe2_0 + 6.0 * ts_xxyy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_xxxy[i] * gfe_0 + ts_xxyy_xxxy[i] * rgc2_0;

        gr_xxyy_xxxz[i] = 2.0 * ts_yy_xxxz[i] * gfe2_0 + 12.0 * ts_xyy_xxz[i] * gfe2_0 + 4.0 * ts_xyy_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxz[i] * gfe2_0 + 4.0 * ts_xxy_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_xz[i] * gfe2_0 + 6.0 * ts_xxyy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xxxz[i] * gfe_0 + ts_xxyy_xxxz[i] * rgc2_0;

        gr_xxyy_xxyy[i] = 2.0 * ts_yy_xxyy[i] * gfe2_0 + 8.0 * ts_xyy_xyy[i] * gfe2_0 + 4.0 * ts_xyy_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxyy[i] * gfe2_0 + 8.0 * ts_xxy_xxy[i] * gfe2_0 + 4.0 * ts_xxy_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yy[i] * gfe2_0 + 4.0 * ts_xxyy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xx[i] * gfe2_0 + 4.0 * ts_xxyy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_xxyy[i] * gfe_0 + ts_xxyy_xxyy[i] * rgc2_0;

        gr_xxyy_xxyz[i] = 2.0 * ts_yy_xxyz[i] * gfe2_0 + 8.0 * ts_xyy_xyz[i] * gfe2_0 + 4.0 * ts_xyy_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxyz[i] * gfe2_0 + 4.0 * ts_xxy_xxz[i] * gfe2_0 + 4.0 * ts_xxy_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yz[i] * gfe2_0 + 4.0 * ts_xxyy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xxyz[i] * gfe_0 + ts_xxyy_xxyz[i] * rgc2_0;

        gr_xxyy_xxzz[i] = 2.0 * ts_yy_xxzz[i] * gfe2_0 + 8.0 * ts_xyy_xzz[i] * gfe2_0 + 4.0 * ts_xyy_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxzz[i] * gfe2_0 + 4.0 * ts_xxy_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_zz[i] * gfe2_0 + 4.0 * ts_xxyy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xx[i] * gfe2_0 + 4.0 * ts_xxyy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xxzz[i] * gfe_0 + ts_xxyy_xxzz[i] * rgc2_0;

        gr_xxyy_xyyy[i] = 2.0 * ts_yy_xyyy[i] * gfe2_0 + 4.0 * ts_xyy_yyy[i] * gfe2_0 + 4.0 * ts_xyy_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyyy[i] * gfe2_0 + 12.0 * ts_xxy_xyy[i] * gfe2_0 + 4.0 * ts_xxy_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxyy_xy[i] * gfe2_0 + 6.0 * ts_xxyy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_xyyy[i] * gfe_0 + ts_xxyy_xyyy[i] * rgc2_0;

        gr_xxyy_xyyz[i] = 2.0 * ts_yy_xyyz[i] * gfe2_0 + 4.0 * ts_xyy_yyz[i] * gfe2_0 + 4.0 * ts_xyy_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyyz[i] * gfe2_0 + 8.0 * ts_xxy_xyz[i] * gfe2_0 + 4.0 * ts_xxy_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xz[i] * gfe2_0 + 4.0 * ts_xxyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xyyz[i] * gfe_0 + ts_xxyy_xyyz[i] * rgc2_0;

        gr_xxyy_xyzz[i] = 2.0 * ts_yy_xyzz[i] * gfe2_0 + 4.0 * ts_xyy_yzz[i] * gfe2_0 + 4.0 * ts_xyy_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyzz[i] * gfe2_0 + 4.0 * ts_xxy_xzz[i] * gfe2_0 + 4.0 * ts_xxy_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_xy[i] * gfe2_0 + 4.0 * ts_xxyy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xyzz[i] * gfe_0 + ts_xxyy_xyzz[i] * rgc2_0;

        gr_xxyy_xzzz[i] = 2.0 * ts_yy_xzzz[i] * gfe2_0 + 4.0 * ts_xyy_zzz[i] * gfe2_0 + 4.0 * ts_xyy_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzzz[i] * gfe2_0 + 4.0 * ts_xxy_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxyy_xz[i] * gfe2_0 + 6.0 * ts_xxyy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xzzz[i] * gfe_0 + ts_xxyy_xzzz[i] * rgc2_0;

        gr_xxyy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe2_0 + 4.0 * ts_xyy_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyyy[i] * gfe2_0 + 16.0 * ts_xxy_yyy[i] * gfe2_0 + 4.0 * ts_xxy_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxyy_yy[i] * gfe2_0 + 8.0 * ts_xxyy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_yyyy[i] * gfe_0 + ts_xxyy_yyyy[i] * rgc2_0;

        gr_xxyy_yyyz[i] = 2.0 * ts_yy_yyyz[i] * gfe2_0 + 4.0 * ts_xyy_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyyz[i] * gfe2_0 + 12.0 * ts_xxy_yyz[i] * gfe2_0 + 4.0 * ts_xxy_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_yz[i] * gfe2_0 + 6.0 * ts_xxyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_yyyz[i] * gfe_0 + ts_xxyy_yyyz[i] * rgc2_0;

        gr_xxyy_yyzz[i] = 2.0 * ts_yy_yyzz[i] * gfe2_0 + 4.0 * ts_xyy_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyzz[i] * gfe2_0 + 8.0 * ts_xxy_yzz[i] * gfe2_0 + 4.0 * ts_xxy_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_zz[i] * gfe2_0 + 4.0 * ts_xxyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yy[i] * gfe2_0 + 4.0 * ts_xxyy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_yyzz[i] * gfe_0 + ts_xxyy_yyzz[i] * rgc2_0;

        gr_xxyy_yzzz[i] = 2.0 * ts_yy_yzzz[i] * gfe2_0 + 4.0 * ts_xyy_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yzzz[i] * gfe2_0 + 4.0 * ts_xxy_zzz[i] * gfe2_0 + 4.0 * ts_xxy_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_yz[i] * gfe2_0 + 6.0 * ts_xxyy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_yzzz[i] * gfe_0 + ts_xxyy_yzzz[i] * rgc2_0;

        gr_xxyy_zzzz[i] = 2.0 * ts_yy_zzzz[i] * gfe2_0 + 4.0 * ts_xyy_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzzz[i] * gfe2_0 + 4.0 * ts_xxy_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_xxyy_zz[i] * gfe2_0 + 8.0 * ts_xxyy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_zzzz[i] * gfe_0 + ts_xxyy_zzzz[i] * rgc2_0;
    }

    // Set up 60-75 components of targeted buffer : GG

    auto gr_xxyz_xxxx = pbuffer.data(idx_g_gg + 60);

    auto gr_xxyz_xxxy = pbuffer.data(idx_g_gg + 61);

    auto gr_xxyz_xxxz = pbuffer.data(idx_g_gg + 62);

    auto gr_xxyz_xxyy = pbuffer.data(idx_g_gg + 63);

    auto gr_xxyz_xxyz = pbuffer.data(idx_g_gg + 64);

    auto gr_xxyz_xxzz = pbuffer.data(idx_g_gg + 65);

    auto gr_xxyz_xyyy = pbuffer.data(idx_g_gg + 66);

    auto gr_xxyz_xyyz = pbuffer.data(idx_g_gg + 67);

    auto gr_xxyz_xyzz = pbuffer.data(idx_g_gg + 68);

    auto gr_xxyz_xzzz = pbuffer.data(idx_g_gg + 69);

    auto gr_xxyz_yyyy = pbuffer.data(idx_g_gg + 70);

    auto gr_xxyz_yyyz = pbuffer.data(idx_g_gg + 71);

    auto gr_xxyz_yyzz = pbuffer.data(idx_g_gg + 72);

    auto gr_xxyz_yzzz = pbuffer.data(idx_g_gg + 73);

    auto gr_xxyz_zzzz = pbuffer.data(idx_g_gg + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_xxxx, gr_xxyz_xxxy, gr_xxyz_xxxz, gr_xxyz_xxyy, gr_xxyz_xxyz, gr_xxyz_xxzz, gr_xxyz_xyyy, gr_xxyz_xyyz, gr_xxyz_xyzz, gr_xxyz_xzzz, gr_xxyz_yyyy, gr_xxyz_yyyz, gr_xxyz_yyzz, gr_xxyz_yzzz, gr_xxyz_zzzz, ts_xxy_xxx, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxy, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxz, ts_xxy_xxzz, ts_xxy_xyy, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyz, ts_xxy_xyzz, ts_xxy_xzz, ts_xxy_xzzz, ts_xxy_yyy, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyz, ts_xxy_yyzz, ts_xxy_yzz, ts_xxy_yzzz, ts_xxy_zzz, ts_xxy_zzzz, ts_xxyz_xx, ts_xxyz_xxx, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxy, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxz, ts_xxyz_xxzz, ts_xxyz_xy, ts_xxyz_xyy, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyz, ts_xxyz_xyzz, ts_xxyz_xz, ts_xxyz_xzz, ts_xxyz_xzzz, ts_xxyz_yy, ts_xxyz_yyy, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyz, ts_xxyz_yyzz, ts_xxyz_yz, ts_xxyz_yzz, ts_xxyz_yzzz, ts_xxyz_zz, ts_xxyz_zzz, ts_xxyz_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxyz_xxxx[i] = 2.0 * ts_yz_xxxx[i] * gfe2_0 + 16.0 * ts_xyz_xxx[i] * gfe2_0 + 4.0 * ts_xyz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxyz_xx[i] * gfe2_0 + 8.0 * ts_xxyz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyz_xxxx[i] * gfe_0 + ts_xxyz_xxxx[i] * rgc2_0;

        gr_xxyz_xxxy[i] = 2.0 * ts_yz_xxxy[i] * gfe2_0 + 12.0 * ts_xyz_xxy[i] * gfe2_0 + 4.0 * ts_xyz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxx[i] * gfe2_0 + 2.0 * ts_xxz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxyz_xy[i] * gfe2_0 + 6.0 * ts_xxyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_xxxy[i] * gfe_0 + ts_xxyz_xxxy[i] * rgc2_0;

        gr_xxyz_xxxz[i] = 2.0 * ts_yz_xxxz[i] * gfe2_0 + 12.0 * ts_xyz_xxz[i] * gfe2_0 + 4.0 * ts_xyz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxx[i] * gfe2_0 + 2.0 * ts_xxy_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxyz_xz[i] * gfe2_0 + 6.0 * ts_xxyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xxxz[i] * gfe_0 + ts_xxyz_xxxz[i] * rgc2_0;

        gr_xxyz_xxyy[i] = 2.0 * ts_yz_xxyy[i] * gfe2_0 + 8.0 * ts_xyz_xyy[i] * gfe2_0 + 4.0 * ts_xyz_xxyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxz_xxy[i] * gfe2_0 + 2.0 * ts_xxz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yy[i] * gfe2_0 + 4.0 * ts_xxyz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xx[i] * gfe2_0 + 4.0 * ts_xxyz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_xxyy[i] * gfe_0 + ts_xxyz_xxyy[i] * rgc2_0;

        gr_xxyz_xxyz[i] = 2.0 * ts_yz_xxyz[i] * gfe2_0 + 8.0 * ts_xyz_xyz[i] * gfe2_0 + 4.0 * ts_xyz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxz[i] * gfe2_0 + 2.0 * ts_xxz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxy[i] * gfe2_0 + 2.0 * ts_xxy_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yz[i] * gfe2_0 + 4.0 * ts_xxyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xxyz[i] * gfe_0 + ts_xxyz_xxyz[i] * rgc2_0;

        gr_xxyz_xxzz[i] = 2.0 * ts_yz_xxzz[i] * gfe2_0 + 8.0 * ts_xyz_xzz[i] * gfe2_0 + 4.0 * ts_xyz_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xxy_xxz[i] * gfe2_0 + 2.0 * ts_xxy_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_zz[i] * gfe2_0 + 4.0 * ts_xxyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xx[i] * gfe2_0 + 4.0 * ts_xxyz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xxzz[i] * gfe_0 + ts_xxyz_xxzz[i] * rgc2_0;

        gr_xxyz_xyyy[i] = 2.0 * ts_yz_xyyy[i] * gfe2_0 + 4.0 * ts_xyz_yyy[i] * gfe2_0 + 4.0 * ts_xyz_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxz_xyy[i] * gfe2_0 + 2.0 * ts_xxz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxyz_xy[i] * gfe2_0 + 6.0 * ts_xxyz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_xyyy[i] * gfe_0 + ts_xxyz_xyyy[i] * rgc2_0;

        gr_xxyz_xyyz[i] = 2.0 * ts_yz_xyyz[i] * gfe2_0 + 4.0 * ts_xyz_yyz[i] * gfe2_0 + 4.0 * ts_xyz_xyyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxz_xyz[i] * gfe2_0 + 2.0 * ts_xxz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xyy[i] * gfe2_0 + 2.0 * ts_xxy_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xz[i] * gfe2_0 + 4.0 * ts_xxyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xyyz[i] * gfe_0 + ts_xxyz_xyyz[i] * rgc2_0;

        gr_xxyz_xyzz[i] = 2.0 * ts_yz_xyzz[i] * gfe2_0 + 4.0 * ts_xyz_yzz[i] * gfe2_0 + 4.0 * ts_xyz_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xzz[i] * gfe2_0 + 2.0 * ts_xxz_xyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xxy_xyz[i] * gfe2_0 + 2.0 * ts_xxy_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_xy[i] * gfe2_0 + 4.0 * ts_xxyz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xyzz[i] * gfe_0 + ts_xxyz_xyzz[i] * rgc2_0;

        gr_xxyz_xzzz[i] = 2.0 * ts_yz_xzzz[i] * gfe2_0 + 4.0 * ts_xyz_zzz[i] * gfe2_0 + 4.0 * ts_xyz_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_xzz[i] * gfe2_0 + 2.0 * ts_xxy_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxyz_xz[i] * gfe2_0 + 6.0 * ts_xxyz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xzzz[i] * gfe_0 + ts_xxyz_xzzz[i] * rgc2_0;

        gr_xxyz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe2_0 + 4.0 * ts_xyz_yyyy[i] * gfe_0 * gc_x[i] + 8.0 * ts_xxz_yyy[i] * gfe2_0 + 2.0 * ts_xxz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxyz_yy[i] * gfe2_0 + 8.0 * ts_xxyz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_yyyy[i] * gfe_0 + ts_xxyz_yyyy[i] * rgc2_0;

        gr_xxyz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe2_0 + 4.0 * ts_xyz_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxz_yyz[i] * gfe2_0 + 2.0 * ts_xxz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yyy[i] * gfe2_0 + 2.0 * ts_xxy_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxyz_yz[i] * gfe2_0 + 6.0 * ts_xxyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_yyyz[i] * gfe_0 + ts_xxyz_yyyz[i] * rgc2_0;

        gr_xxyz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe2_0 + 4.0 * ts_xyz_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxz_yzz[i] * gfe2_0 + 2.0 * ts_xxz_yyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xxy_yyz[i] * gfe2_0 + 2.0 * ts_xxy_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_zz[i] * gfe2_0 + 4.0 * ts_xxyz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_yy[i] * gfe2_0 + 4.0 * ts_xxyz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_yyzz[i] * gfe_0 + ts_xxyz_yyzz[i] * rgc2_0;

        gr_xxyz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe2_0 + 4.0 * ts_xyz_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_zzz[i] * gfe2_0 + 2.0 * ts_xxz_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_yzz[i] * gfe2_0 + 2.0 * ts_xxy_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyz_yz[i] * gfe2_0 + 6.0 * ts_xxyz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_yzzz[i] * gfe_0 + ts_xxyz_yzzz[i] * rgc2_0;

        gr_xxyz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe2_0 + 4.0 * ts_xyz_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_zzzz[i] * gfe_0 * gc_y[i] + 8.0 * ts_xxy_zzz[i] * gfe2_0 + 2.0 * ts_xxy_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxyz_zz[i] * gfe2_0 + 8.0 * ts_xxyz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_zzzz[i] * gfe_0 + ts_xxyz_zzzz[i] * rgc2_0;
    }

    // Set up 75-90 components of targeted buffer : GG

    auto gr_xxzz_xxxx = pbuffer.data(idx_g_gg + 75);

    auto gr_xxzz_xxxy = pbuffer.data(idx_g_gg + 76);

    auto gr_xxzz_xxxz = pbuffer.data(idx_g_gg + 77);

    auto gr_xxzz_xxyy = pbuffer.data(idx_g_gg + 78);

    auto gr_xxzz_xxyz = pbuffer.data(idx_g_gg + 79);

    auto gr_xxzz_xxzz = pbuffer.data(idx_g_gg + 80);

    auto gr_xxzz_xyyy = pbuffer.data(idx_g_gg + 81);

    auto gr_xxzz_xyyz = pbuffer.data(idx_g_gg + 82);

    auto gr_xxzz_xyzz = pbuffer.data(idx_g_gg + 83);

    auto gr_xxzz_xzzz = pbuffer.data(idx_g_gg + 84);

    auto gr_xxzz_yyyy = pbuffer.data(idx_g_gg + 85);

    auto gr_xxzz_yyyz = pbuffer.data(idx_g_gg + 86);

    auto gr_xxzz_yyzz = pbuffer.data(idx_g_gg + 87);

    auto gr_xxzz_yzzz = pbuffer.data(idx_g_gg + 88);

    auto gr_xxzz_zzzz = pbuffer.data(idx_g_gg + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_xxxx, gr_xxzz_xxxy, gr_xxzz_xxxz, gr_xxzz_xxyy, gr_xxzz_xxyz, gr_xxzz_xxzz, gr_xxzz_xyyy, gr_xxzz_xyyz, gr_xxzz_xyzz, gr_xxzz_xzzz, gr_xxzz_yyyy, gr_xxzz_yyyz, gr_xxzz_yyzz, gr_xxzz_yzzz, gr_xxzz_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxz_xxx, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxy, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyy, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_yyy, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyz, ts_xxz_yyzz, ts_xxz_yzz, ts_xxz_yzzz, ts_xxz_zzz, ts_xxz_zzzz, ts_xxzz_xx, ts_xxzz_xxx, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxy, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxz, ts_xxzz_xxzz, ts_xxzz_xy, ts_xxzz_xyy, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyz, ts_xxzz_xyzz, ts_xxzz_xz, ts_xxzz_xzz, ts_xxzz_xzzz, ts_xxzz_yy, ts_xxzz_yyy, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyz, ts_xxzz_yyzz, ts_xxzz_yz, ts_xxzz_yzz, ts_xxzz_yzzz, ts_xxzz_zz, ts_xxzz_zzz, ts_xxzz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxzz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gfe2_0 + 16.0 * ts_xzz_xxx[i] * gfe2_0 + 4.0 * ts_xzz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxx[i] * gfe2_0 + 4.0 * ts_xxz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxzz_xx[i] * gfe2_0 + 8.0 * ts_xxzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxzz_xxxx[i] * gfe_0 + ts_xxzz_xxxx[i] * rgc2_0;

        gr_xxzz_xxxy[i] = 2.0 * ts_zz_xxxy[i] * gfe2_0 + 12.0 * ts_xzz_xxy[i] * gfe2_0 + 4.0 * ts_xzz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxy[i] * gfe2_0 + 4.0 * ts_xxz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxzz_xy[i] * gfe2_0 + 6.0 * ts_xxzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_xxxy[i] * gfe_0 + ts_xxzz_xxxy[i] * rgc2_0;

        gr_xxzz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gfe2_0 + 12.0 * ts_xzz_xxz[i] * gfe2_0 + 4.0 * ts_xzz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxxz[i] * gfe2_0 + 4.0 * ts_xxz_xxx[i] * gfe2_0 + 4.0 * ts_xxz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxzz_xz[i] * gfe2_0 + 6.0 * ts_xxzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xxxz[i] * gfe_0 + ts_xxzz_xxxz[i] * rgc2_0;

        gr_xxzz_xxyy[i] = 2.0 * ts_zz_xxyy[i] * gfe2_0 + 8.0 * ts_xzz_xyy[i] * gfe2_0 + 4.0 * ts_xzz_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxyy[i] * gfe2_0 + 4.0 * ts_xxz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yy[i] * gfe2_0 + 4.0 * ts_xxzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xx[i] * gfe2_0 + 4.0 * ts_xxzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_xxyy[i] * gfe_0 + ts_xxzz_xxyy[i] * rgc2_0;

        gr_xxzz_xxyz[i] = 2.0 * ts_zz_xxyz[i] * gfe2_0 + 8.0 * ts_xzz_xyz[i] * gfe2_0 + 4.0 * ts_xzz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxyz[i] * gfe2_0 + 4.0 * ts_xxz_xxy[i] * gfe2_0 + 4.0 * ts_xxz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yz[i] * gfe2_0 + 4.0 * ts_xxzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xxyz[i] * gfe_0 + ts_xxzz_xxyz[i] * rgc2_0;

        gr_xxzz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gfe2_0 + 8.0 * ts_xzz_xzz[i] * gfe2_0 + 4.0 * ts_xzz_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxzz[i] * gfe2_0 + 8.0 * ts_xxz_xxz[i] * gfe2_0 + 4.0 * ts_xxz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_zz[i] * gfe2_0 + 4.0 * ts_xxzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xx[i] * gfe2_0 + 4.0 * ts_xxzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xxzz[i] * gfe_0 + ts_xxzz_xxzz[i] * rgc2_0;

        gr_xxzz_xyyy[i] = 2.0 * ts_zz_xyyy[i] * gfe2_0 + 4.0 * ts_xzz_yyy[i] * gfe2_0 + 4.0 * ts_xzz_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyyy[i] * gfe2_0 + 4.0 * ts_xxz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxzz_xy[i] * gfe2_0 + 6.0 * ts_xxzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_xyyy[i] * gfe_0 + ts_xxzz_xyyy[i] * rgc2_0;

        gr_xxzz_xyyz[i] = 2.0 * ts_zz_xyyz[i] * gfe2_0 + 4.0 * ts_xzz_yyz[i] * gfe2_0 + 4.0 * ts_xzz_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyyz[i] * gfe2_0 + 4.0 * ts_xxz_xyy[i] * gfe2_0 + 4.0 * ts_xxz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xz[i] * gfe2_0 + 4.0 * ts_xxzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xyyz[i] * gfe_0 + ts_xxzz_xyyz[i] * rgc2_0;

        gr_xxzz_xyzz[i] = 2.0 * ts_zz_xyzz[i] * gfe2_0 + 4.0 * ts_xzz_yzz[i] * gfe2_0 + 4.0 * ts_xzz_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyzz[i] * gfe2_0 + 8.0 * ts_xxz_xyz[i] * gfe2_0 + 4.0 * ts_xxz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_xy[i] * gfe2_0 + 4.0 * ts_xxzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xyzz[i] * gfe_0 + ts_xxzz_xyzz[i] * rgc2_0;

        gr_xxzz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe2_0 + 4.0 * ts_xzz_zzz[i] * gfe2_0 + 4.0 * ts_xzz_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzzz[i] * gfe2_0 + 12.0 * ts_xxz_xzz[i] * gfe2_0 + 4.0 * ts_xxz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxzz_xz[i] * gfe2_0 + 6.0 * ts_xxzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xzzz[i] * gfe_0 + ts_xxzz_xzzz[i] * rgc2_0;

        gr_xxzz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe2_0 + 4.0 * ts_xzz_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyyy[i] * gfe2_0 + 4.0 * ts_xxz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxzz_yy[i] * gfe2_0 + 8.0 * ts_xxzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_yyyy[i] * gfe_0 + ts_xxzz_yyyy[i] * rgc2_0;

        gr_xxzz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe2_0 + 4.0 * ts_xzz_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyyz[i] * gfe2_0 + 4.0 * ts_xxz_yyy[i] * gfe2_0 + 4.0 * ts_xxz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxzz_yz[i] * gfe2_0 + 6.0 * ts_xxzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_yyyz[i] * gfe_0 + ts_xxzz_yyyz[i] * rgc2_0;

        gr_xxzz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe2_0 + 4.0 * ts_xzz_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyzz[i] * gfe2_0 + 8.0 * ts_xxz_yyz[i] * gfe2_0 + 4.0 * ts_xxz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_zz[i] * gfe2_0 + 4.0 * ts_xxzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_yy[i] * gfe2_0 + 4.0 * ts_xxzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_yyzz[i] * gfe_0 + ts_xxzz_yyzz[i] * rgc2_0;

        gr_xxzz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe2_0 + 4.0 * ts_xzz_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yzzz[i] * gfe2_0 + 12.0 * ts_xxz_yzz[i] * gfe2_0 + 4.0 * ts_xxz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxzz_yz[i] * gfe2_0 + 6.0 * ts_xxzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_yzzz[i] * gfe_0 + ts_xxzz_yzzz[i] * rgc2_0;

        gr_xxzz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe2_0 + 4.0 * ts_xzz_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzzz[i] * gfe2_0 + 16.0 * ts_xxz_zzz[i] * gfe2_0 + 4.0 * ts_xxz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xxzz_zz[i] * gfe2_0 + 8.0 * ts_xxzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_zzzz[i] * gfe_0 + ts_xxzz_zzzz[i] * rgc2_0;
    }

    // Set up 90-105 components of targeted buffer : GG

    auto gr_xyyy_xxxx = pbuffer.data(idx_g_gg + 90);

    auto gr_xyyy_xxxy = pbuffer.data(idx_g_gg + 91);

    auto gr_xyyy_xxxz = pbuffer.data(idx_g_gg + 92);

    auto gr_xyyy_xxyy = pbuffer.data(idx_g_gg + 93);

    auto gr_xyyy_xxyz = pbuffer.data(idx_g_gg + 94);

    auto gr_xyyy_xxzz = pbuffer.data(idx_g_gg + 95);

    auto gr_xyyy_xyyy = pbuffer.data(idx_g_gg + 96);

    auto gr_xyyy_xyyz = pbuffer.data(idx_g_gg + 97);

    auto gr_xyyy_xyzz = pbuffer.data(idx_g_gg + 98);

    auto gr_xyyy_xzzz = pbuffer.data(idx_g_gg + 99);

    auto gr_xyyy_yyyy = pbuffer.data(idx_g_gg + 100);

    auto gr_xyyy_yyyz = pbuffer.data(idx_g_gg + 101);

    auto gr_xyyy_yyzz = pbuffer.data(idx_g_gg + 102);

    auto gr_xyyy_yzzz = pbuffer.data(idx_g_gg + 103);

    auto gr_xyyy_zzzz = pbuffer.data(idx_g_gg + 104);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_xxxx, gr_xyyy_xxxy, gr_xyyy_xxxz, gr_xyyy_xxyy, gr_xyyy_xxyz, gr_xyyy_xxzz, gr_xyyy_xyyy, gr_xyyy_xyyz, gr_xyyy_xyzz, gr_xyyy_xzzz, gr_xyyy_yyyy, gr_xyyy_yyyz, gr_xyyy_yyzz, gr_xyyy_yzzz, gr_xyyy_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, ts_xyyy_xx, ts_xyyy_xxx, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxy, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxz, ts_xyyy_xxzz, ts_xyyy_xy, ts_xyyy_xyy, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyz, ts_xyyy_xyzz, ts_xyyy_xz, ts_xyyy_xzz, ts_xyyy_xzzz, ts_xyyy_yy, ts_xyyy_yyy, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyz, ts_xyyy_yyzz, ts_xyyy_yz, ts_xyyy_yzz, ts_xyyy_yzzz, ts_xyyy_zz, ts_xyyy_zzz, ts_xyyy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyyy_xxxx[i] = 8.0 * ts_yyy_xxx[i] * gfe2_0 + 2.0 * ts_yyy_xxxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxxx[i] * gfe2_0 + 6.0 * ts_xyy_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_xyyy_xx[i] * gfe2_0 + 8.0 * ts_xyyy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyy_xxxx[i] * gfe_0 + ts_xyyy_xxxx[i] * rgc2_0;

        gr_xyyy_xxxy[i] = 6.0 * ts_yyy_xxy[i] * gfe2_0 + 2.0 * ts_yyy_xxxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxxy[i] * gfe2_0 + 6.0 * ts_xyy_xxx[i] * gfe2_0 + 6.0 * ts_xyy_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_xy[i] * gfe2_0 + 6.0 * ts_xyyy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_xxxy[i] * gfe_0 + ts_xyyy_xxxy[i] * rgc2_0;

        gr_xyyy_xxxz[i] = 6.0 * ts_yyy_xxz[i] * gfe2_0 + 2.0 * ts_yyy_xxxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxxz[i] * gfe2_0 + 6.0 * ts_xyy_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_xz[i] * gfe2_0 + 6.0 * ts_xyyy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xxxz[i] * gfe_0 + ts_xyyy_xxxz[i] * rgc2_0;

        gr_xyyy_xxyy[i] = 4.0 * ts_yyy_xyy[i] * gfe2_0 + 2.0 * ts_yyy_xxyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxyy[i] * gfe2_0 + 12.0 * ts_xyy_xxy[i] * gfe2_0 + 6.0 * ts_xyy_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yy[i] * gfe2_0 + 4.0 * ts_xyyy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xx[i] * gfe2_0 + 4.0 * ts_xyyy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_xxyy[i] * gfe_0 + ts_xyyy_xxyy[i] * rgc2_0;

        gr_xyyy_xxyz[i] = 4.0 * ts_yyy_xyz[i] * gfe2_0 + 2.0 * ts_yyy_xxyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxyz[i] * gfe2_0 + 6.0 * ts_xyy_xxz[i] * gfe2_0 + 6.0 * ts_xyy_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yz[i] * gfe2_0 + 4.0 * ts_xyyy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xxyz[i] * gfe_0 + ts_xyyy_xxyz[i] * rgc2_0;

        gr_xyyy_xxzz[i] = 4.0 * ts_yyy_xzz[i] * gfe2_0 + 2.0 * ts_yyy_xxzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxzz[i] * gfe2_0 + 6.0 * ts_xyy_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_zz[i] * gfe2_0 + 4.0 * ts_xyyy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xx[i] * gfe2_0 + 4.0 * ts_xyyy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xxzz[i] * gfe_0 + ts_xyyy_xxzz[i] * rgc2_0;

        gr_xyyy_xyyy[i] = 2.0 * ts_yyy_yyy[i] * gfe2_0 + 2.0 * ts_yyy_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xyyy[i] * gfe2_0 + 18.0 * ts_xyy_xyy[i] * gfe2_0 + 6.0 * ts_xyy_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyyy_xy[i] * gfe2_0 + 6.0 * ts_xyyy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_xyyy[i] * gfe_0 + ts_xyyy_xyyy[i] * rgc2_0;

        gr_xyyy_xyyz[i] = 2.0 * ts_yyy_yyz[i] * gfe2_0 + 2.0 * ts_yyy_xyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xyyz[i] * gfe2_0 + 12.0 * ts_xyy_xyz[i] * gfe2_0 + 6.0 * ts_xyy_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xz[i] * gfe2_0 + 4.0 * ts_xyyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xyyz[i] * gfe_0 + ts_xyyy_xyyz[i] * rgc2_0;

        gr_xyyy_xyzz[i] = 2.0 * ts_yyy_yzz[i] * gfe2_0 + 2.0 * ts_yyy_xyzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xyzz[i] * gfe2_0 + 6.0 * ts_xyy_xzz[i] * gfe2_0 + 6.0 * ts_xyy_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_xy[i] * gfe2_0 + 4.0 * ts_xyyy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xyzz[i] * gfe_0 + ts_xyyy_xyzz[i] * rgc2_0;

        gr_xyyy_xzzz[i] = 2.0 * ts_yyy_zzz[i] * gfe2_0 + 2.0 * ts_yyy_xzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xzzz[i] * gfe2_0 + 6.0 * ts_xyy_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyyy_xz[i] * gfe2_0 + 6.0 * ts_xyyy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xzzz[i] * gfe_0 + ts_xyyy_xzzz[i] * rgc2_0;

        gr_xyyy_yyyy[i] = 2.0 * ts_yyy_yyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yyyy[i] * gfe2_0 + 24.0 * ts_xyy_yyy[i] * gfe2_0 + 6.0 * ts_xyy_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_xyyy_yy[i] * gfe2_0 + 8.0 * ts_xyyy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_yyyy[i] * gfe_0 + ts_xyyy_yyyy[i] * rgc2_0;

        gr_xyyy_yyyz[i] = 2.0 * ts_yyy_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yyyz[i] * gfe2_0 + 18.0 * ts_xyy_yyz[i] * gfe2_0 + 6.0 * ts_xyy_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_yz[i] * gfe2_0 + 6.0 * ts_xyyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_yyyz[i] * gfe_0 + ts_xyyy_yyyz[i] * rgc2_0;

        gr_xyyy_yyzz[i] = 2.0 * ts_yyy_yyzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yyzz[i] * gfe2_0 + 12.0 * ts_xyy_yzz[i] * gfe2_0 + 6.0 * ts_xyy_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_zz[i] * gfe2_0 + 4.0 * ts_xyyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yy[i] * gfe2_0 + 4.0 * ts_xyyy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_yyzz[i] * gfe_0 + ts_xyyy_yyzz[i] * rgc2_0;

        gr_xyyy_yzzz[i] = 2.0 * ts_yyy_yzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yzzz[i] * gfe2_0 + 6.0 * ts_xyy_zzz[i] * gfe2_0 + 6.0 * ts_xyy_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_yz[i] * gfe2_0 + 6.0 * ts_xyyy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_yzzz[i] * gfe_0 + ts_xyyy_yzzz[i] * rgc2_0;

        gr_xyyy_zzzz[i] = 2.0 * ts_yyy_zzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_zzzz[i] * gfe2_0 + 6.0 * ts_xyy_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_xyyy_zz[i] * gfe2_0 + 8.0 * ts_xyyy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_zzzz[i] * gfe_0 + ts_xyyy_zzzz[i] * rgc2_0;
    }

    // Set up 105-120 components of targeted buffer : GG

    auto gr_xyyz_xxxx = pbuffer.data(idx_g_gg + 105);

    auto gr_xyyz_xxxy = pbuffer.data(idx_g_gg + 106);

    auto gr_xyyz_xxxz = pbuffer.data(idx_g_gg + 107);

    auto gr_xyyz_xxyy = pbuffer.data(idx_g_gg + 108);

    auto gr_xyyz_xxyz = pbuffer.data(idx_g_gg + 109);

    auto gr_xyyz_xxzz = pbuffer.data(idx_g_gg + 110);

    auto gr_xyyz_xyyy = pbuffer.data(idx_g_gg + 111);

    auto gr_xyyz_xyyz = pbuffer.data(idx_g_gg + 112);

    auto gr_xyyz_xyzz = pbuffer.data(idx_g_gg + 113);

    auto gr_xyyz_xzzz = pbuffer.data(idx_g_gg + 114);

    auto gr_xyyz_yyyy = pbuffer.data(idx_g_gg + 115);

    auto gr_xyyz_yyyz = pbuffer.data(idx_g_gg + 116);

    auto gr_xyyz_yyzz = pbuffer.data(idx_g_gg + 117);

    auto gr_xyyz_yzzz = pbuffer.data(idx_g_gg + 118);

    auto gr_xyyz_zzzz = pbuffer.data(idx_g_gg + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_xxxx, gr_xyyz_xxxy, gr_xyyz_xxxz, gr_xyyz_xxyy, gr_xyyz_xxyz, gr_xyyz_xxzz, gr_xyyz_xyyy, gr_xyyz_xyyz, gr_xyyz_xyzz, gr_xyyz_xzzz, gr_xyyz_yyyy, gr_xyyz_yyyz, gr_xyyz_yyzz, gr_xyyz_yzzz, gr_xyyz_zzzz, ts_xyy_xxx, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxz, ts_xyy_xxzz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_xzz, ts_xyy_xzzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzz, ts_xyy_zzzz, ts_xyyz_xx, ts_xyyz_xxx, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxy, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxz, ts_xyyz_xxzz, ts_xyyz_xy, ts_xyyz_xyy, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyz, ts_xyyz_xyzz, ts_xyyz_xz, ts_xyyz_xzz, ts_xyyz_xzzz, ts_xyyz_yy, ts_xyyz_yyy, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyz, ts_xyyz_yyzz, ts_xyyz_yz, ts_xyyz_yzz, ts_xyyz_yzzz, ts_xyyz_zz, ts_xyyz_zzz, ts_xyyz_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyyz_xxxx[i] = 8.0 * ts_yyz_xxx[i] * gfe2_0 + 2.0 * ts_yyz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxxx[i] * gfe2_0 + 4.0 * ts_xyz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyyz_xx[i] * gfe2_0 + 8.0 * ts_xyyz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyz_xxxx[i] * gfe_0 + ts_xyyz_xxxx[i] * rgc2_0;

        gr_xyyz_xxxy[i] = 6.0 * ts_yyz_xxy[i] * gfe2_0 + 2.0 * ts_yyz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxxy[i] * gfe2_0 + 4.0 * ts_xyz_xxx[i] * gfe2_0 + 4.0 * ts_xyz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyyz_xy[i] * gfe2_0 + 6.0 * ts_xyyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_xxxy[i] * gfe_0 + ts_xyyz_xxxy[i] * rgc2_0;

        gr_xyyz_xxxz[i] = 6.0 * ts_yyz_xxz[i] * gfe2_0 + 2.0 * ts_yyz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxxz[i] * gfe2_0 + 4.0 * ts_xyz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxx[i] * gfe2_0 + 2.0 * ts_xyy_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyyz_xz[i] * gfe2_0 + 6.0 * ts_xyyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xxxz[i] * gfe_0 + ts_xyyz_xxxz[i] * rgc2_0;

        gr_xyyz_xxyy[i] = 4.0 * ts_yyz_xyy[i] * gfe2_0 + 2.0 * ts_yyz_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxyy[i] * gfe2_0 + 8.0 * ts_xyz_xxy[i] * gfe2_0 + 4.0 * ts_xyz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yy[i] * gfe2_0 + 4.0 * ts_xyyz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xx[i] * gfe2_0 + 4.0 * ts_xyyz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_xxyy[i] * gfe_0 + ts_xyyz_xxyy[i] * rgc2_0;

        gr_xyyz_xxyz[i] = 4.0 * ts_yyz_xyz[i] * gfe2_0 + 2.0 * ts_yyz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxyz[i] * gfe2_0 + 4.0 * ts_xyz_xxz[i] * gfe2_0 + 4.0 * ts_xyz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxy[i] * gfe2_0 + 2.0 * ts_xyy_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yz[i] * gfe2_0 + 4.0 * ts_xyyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xxyz[i] * gfe_0 + ts_xyyz_xxyz[i] * rgc2_0;

        gr_xyyz_xxzz[i] = 4.0 * ts_yyz_xzz[i] * gfe2_0 + 2.0 * ts_yyz_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxzz[i] * gfe2_0 + 4.0 * ts_xyz_xxzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xyy_xxz[i] * gfe2_0 + 2.0 * ts_xyy_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_zz[i] * gfe2_0 + 4.0 * ts_xyyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xx[i] * gfe2_0 + 4.0 * ts_xyyz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xxzz[i] * gfe_0 + ts_xyyz_xxzz[i] * rgc2_0;

        gr_xyyz_xyyy[i] = 2.0 * ts_yyz_yyy[i] * gfe2_0 + 2.0 * ts_yyz_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xyyy[i] * gfe2_0 + 12.0 * ts_xyz_xyy[i] * gfe2_0 + 4.0 * ts_xyz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyyz_xy[i] * gfe2_0 + 6.0 * ts_xyyz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_xyyy[i] * gfe_0 + ts_xyyz_xyyy[i] * rgc2_0;

        gr_xyyz_xyyz[i] = 2.0 * ts_yyz_yyz[i] * gfe2_0 + 2.0 * ts_yyz_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xyyz[i] * gfe2_0 + 8.0 * ts_xyz_xyz[i] * gfe2_0 + 4.0 * ts_xyz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xyy[i] * gfe2_0 + 2.0 * ts_xyy_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xz[i] * gfe2_0 + 4.0 * ts_xyyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xyyz[i] * gfe_0 + ts_xyyz_xyyz[i] * rgc2_0;

        gr_xyyz_xyzz[i] = 2.0 * ts_yyz_yzz[i] * gfe2_0 + 2.0 * ts_yyz_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xyzz[i] * gfe2_0 + 4.0 * ts_xyz_xzz[i] * gfe2_0 + 4.0 * ts_xyz_xyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xyy_xyz[i] * gfe2_0 + 2.0 * ts_xyy_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_xy[i] * gfe2_0 + 4.0 * ts_xyyz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xyzz[i] * gfe_0 + ts_xyyz_xyzz[i] * rgc2_0;

        gr_xyyz_xzzz[i] = 2.0 * ts_yyz_zzz[i] * gfe2_0 + 2.0 * ts_yyz_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xzzz[i] * gfe2_0 + 4.0 * ts_xyz_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_xzz[i] * gfe2_0 + 2.0 * ts_xyy_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyyz_xz[i] * gfe2_0 + 6.0 * ts_xyyz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xzzz[i] * gfe_0 + ts_xyyz_xzzz[i] * rgc2_0;

        gr_xyyz_yyyy[i] = 2.0 * ts_yyz_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yyyy[i] * gfe2_0 + 16.0 * ts_xyz_yyy[i] * gfe2_0 + 4.0 * ts_xyz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyyz_yy[i] * gfe2_0 + 8.0 * ts_xyyz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_yyyy[i] * gfe_0 + ts_xyyz_yyyy[i] * rgc2_0;

        gr_xyyz_yyyz[i] = 2.0 * ts_yyz_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yyyz[i] * gfe2_0 + 12.0 * ts_xyz_yyz[i] * gfe2_0 + 4.0 * ts_xyz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yyy[i] * gfe2_0 + 2.0 * ts_xyy_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyyz_yz[i] * gfe2_0 + 6.0 * ts_xyyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_yyyz[i] * gfe_0 + ts_xyyz_yyyz[i] * rgc2_0;

        gr_xyyz_yyzz[i] = 2.0 * ts_yyz_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yyzz[i] * gfe2_0 + 8.0 * ts_xyz_yzz[i] * gfe2_0 + 4.0 * ts_xyz_yyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xyy_yyz[i] * gfe2_0 + 2.0 * ts_xyy_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_zz[i] * gfe2_0 + 4.0 * ts_xyyz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_yy[i] * gfe2_0 + 4.0 * ts_xyyz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_yyzz[i] * gfe_0 + ts_xyyz_yyzz[i] * rgc2_0;

        gr_xyyz_yzzz[i] = 2.0 * ts_yyz_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yzzz[i] * gfe2_0 + 4.0 * ts_xyz_zzz[i] * gfe2_0 + 4.0 * ts_xyz_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_yzz[i] * gfe2_0 + 2.0 * ts_xyy_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyz_yz[i] * gfe2_0 + 6.0 * ts_xyyz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_yzzz[i] * gfe_0 + ts_xyyz_yzzz[i] * rgc2_0;

        gr_xyyz_zzzz[i] = 2.0 * ts_yyz_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zzzz[i] * gfe2_0 + 4.0 * ts_xyz_zzzz[i] * gfe_0 * gc_y[i] + 8.0 * ts_xyy_zzz[i] * gfe2_0 + 2.0 * ts_xyy_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyyz_zz[i] * gfe2_0 + 8.0 * ts_xyyz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_zzzz[i] * gfe_0 + ts_xyyz_zzzz[i] * rgc2_0;
    }

    // Set up 120-135 components of targeted buffer : GG

    auto gr_xyzz_xxxx = pbuffer.data(idx_g_gg + 120);

    auto gr_xyzz_xxxy = pbuffer.data(idx_g_gg + 121);

    auto gr_xyzz_xxxz = pbuffer.data(idx_g_gg + 122);

    auto gr_xyzz_xxyy = pbuffer.data(idx_g_gg + 123);

    auto gr_xyzz_xxyz = pbuffer.data(idx_g_gg + 124);

    auto gr_xyzz_xxzz = pbuffer.data(idx_g_gg + 125);

    auto gr_xyzz_xyyy = pbuffer.data(idx_g_gg + 126);

    auto gr_xyzz_xyyz = pbuffer.data(idx_g_gg + 127);

    auto gr_xyzz_xyzz = pbuffer.data(idx_g_gg + 128);

    auto gr_xyzz_xzzz = pbuffer.data(idx_g_gg + 129);

    auto gr_xyzz_yyyy = pbuffer.data(idx_g_gg + 130);

    auto gr_xyzz_yyyz = pbuffer.data(idx_g_gg + 131);

    auto gr_xyzz_yyzz = pbuffer.data(idx_g_gg + 132);

    auto gr_xyzz_yzzz = pbuffer.data(idx_g_gg + 133);

    auto gr_xyzz_zzzz = pbuffer.data(idx_g_gg + 134);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_xxxx, gr_xyzz_xxxy, gr_xyzz_xxxz, gr_xyzz_xxyy, gr_xyzz_xxyz, gr_xyzz_xxzz, gr_xyzz_xyyy, gr_xyzz_xyyz, gr_xyzz_xyzz, gr_xyzz_xzzz, gr_xyzz_yyyy, gr_xyzz_yyyz, gr_xyzz_yyzz, gr_xyzz_yzzz, gr_xyzz_zzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_xyz_xxx, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxy, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxz, ts_xyz_xxzz, ts_xyz_xyy, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyz, ts_xyz_xyzz, ts_xyz_xzz, ts_xyz_xzzz, ts_xyz_yyy, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyz, ts_xyz_yyzz, ts_xyz_yzz, ts_xyz_yzzz, ts_xyz_zzz, ts_xyz_zzzz, ts_xyzz_xx, ts_xyzz_xxx, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxy, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxz, ts_xyzz_xxzz, ts_xyzz_xy, ts_xyzz_xyy, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyz, ts_xyzz_xyzz, ts_xyzz_xz, ts_xyzz_xzz, ts_xyzz_xzzz, ts_xyzz_yy, ts_xyzz_yyy, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyz, ts_xyzz_yyzz, ts_xyzz_yz, ts_xyzz_yzz, ts_xyzz_yzzz, ts_xyzz_zz, ts_xyzz_zzz, ts_xyzz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyzz_xxxx[i] = 8.0 * ts_yzz_xxx[i] * gfe2_0 + 2.0 * ts_yzz_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxxx[i] * gfe2_0 + 4.0 * ts_xyz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyzz_xx[i] * gfe2_0 + 8.0 * ts_xyzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyzz_xxxx[i] * gfe_0 + ts_xyzz_xxxx[i] * rgc2_0;

        gr_xyzz_xxxy[i] = 6.0 * ts_yzz_xxy[i] * gfe2_0 + 2.0 * ts_yzz_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxx[i] * gfe2_0 + 2.0 * ts_xzz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxxy[i] * gfe2_0 + 4.0 * ts_xyz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyzz_xy[i] * gfe2_0 + 6.0 * ts_xyzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_xxxy[i] * gfe_0 + ts_xyzz_xxxy[i] * rgc2_0;

        gr_xyzz_xxxz[i] = 6.0 * ts_yzz_xxz[i] * gfe2_0 + 2.0 * ts_yzz_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxxz[i] * gfe2_0 + 4.0 * ts_xyz_xxx[i] * gfe2_0 + 4.0 * ts_xyz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyzz_xz[i] * gfe2_0 + 6.0 * ts_xyzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xxxz[i] * gfe_0 + ts_xyzz_xxxz[i] * rgc2_0;

        gr_xyzz_xxyy[i] = 4.0 * ts_yzz_xyy[i] * gfe2_0 + 2.0 * ts_yzz_xxyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xzz_xxy[i] * gfe2_0 + 2.0 * ts_xzz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxyy[i] * gfe2_0 + 4.0 * ts_xyz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yy[i] * gfe2_0 + 4.0 * ts_xyzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xx[i] * gfe2_0 + 4.0 * ts_xyzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_xxyy[i] * gfe_0 + ts_xyzz_xxyy[i] * rgc2_0;

        gr_xyzz_xxyz[i] = 4.0 * ts_yzz_xyz[i] * gfe2_0 + 2.0 * ts_yzz_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxz[i] * gfe2_0 + 2.0 * ts_xzz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxyz[i] * gfe2_0 + 4.0 * ts_xyz_xxy[i] * gfe2_0 + 4.0 * ts_xyz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yz[i] * gfe2_0 + 4.0 * ts_xyzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xxyz[i] * gfe_0 + ts_xyzz_xxyz[i] * rgc2_0;

        gr_xyzz_xxzz[i] = 4.0 * ts_yzz_xzz[i] * gfe2_0 + 2.0 * ts_yzz_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxzz[i] * gfe2_0 + 8.0 * ts_xyz_xxz[i] * gfe2_0 + 4.0 * ts_xyz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_zz[i] * gfe2_0 + 4.0 * ts_xyzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xx[i] * gfe2_0 + 4.0 * ts_xyzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xxzz[i] * gfe_0 + ts_xyzz_xxzz[i] * rgc2_0;

        gr_xyzz_xyyy[i] = 2.0 * ts_yzz_yyy[i] * gfe2_0 + 2.0 * ts_yzz_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzz_xyy[i] * gfe2_0 + 2.0 * ts_xzz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyyy[i] * gfe2_0 + 4.0 * ts_xyz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyzz_xy[i] * gfe2_0 + 6.0 * ts_xyzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_xyyy[i] * gfe_0 + ts_xyzz_xyyy[i] * rgc2_0;

        gr_xyzz_xyyz[i] = 2.0 * ts_yzz_yyz[i] * gfe2_0 + 2.0 * ts_yzz_xyyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xzz_xyz[i] * gfe2_0 + 2.0 * ts_xzz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyyz[i] * gfe2_0 + 4.0 * ts_xyz_xyy[i] * gfe2_0 + 4.0 * ts_xyz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xz[i] * gfe2_0 + 4.0 * ts_xyzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xyyz[i] * gfe_0 + ts_xyzz_xyyz[i] * rgc2_0;

        gr_xyzz_xyzz[i] = 2.0 * ts_yzz_yzz[i] * gfe2_0 + 2.0 * ts_yzz_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xzz[i] * gfe2_0 + 2.0 * ts_xzz_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyzz[i] * gfe2_0 + 8.0 * ts_xyz_xyz[i] * gfe2_0 + 4.0 * ts_xyz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_xy[i] * gfe2_0 + 4.0 * ts_xyzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xyzz[i] * gfe_0 + ts_xyzz_xyzz[i] * rgc2_0;

        gr_xyzz_xzzz[i] = 2.0 * ts_yzz_zzz[i] * gfe2_0 + 2.0 * ts_yzz_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xzzz[i] * gfe2_0 + 12.0 * ts_xyz_xzz[i] * gfe2_0 + 4.0 * ts_xyz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xyzz_xz[i] * gfe2_0 + 6.0 * ts_xyzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xzzz[i] * gfe_0 + ts_xyzz_xzzz[i] * rgc2_0;

        gr_xyzz_yyyy[i] = 2.0 * ts_yzz_yyyy[i] * gfe_0 * gc_x[i] + 8.0 * ts_xzz_yyy[i] * gfe2_0 + 2.0 * ts_xzz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyyy[i] * gfe2_0 + 4.0 * ts_xyz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyzz_yy[i] * gfe2_0 + 8.0 * ts_xyzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_yyyy[i] * gfe_0 + ts_xyzz_yyyy[i] * rgc2_0;

        gr_xyzz_yyyz[i] = 2.0 * ts_yzz_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzz_yyz[i] * gfe2_0 + 2.0 * ts_xzz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyyz[i] * gfe2_0 + 4.0 * ts_xyz_yyy[i] * gfe2_0 + 4.0 * ts_xyz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyzz_yz[i] * gfe2_0 + 6.0 * ts_xyzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_yyyz[i] * gfe_0 + ts_xyzz_yyyz[i] * rgc2_0;

        gr_xyzz_yyzz[i] = 2.0 * ts_yzz_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xzz_yzz[i] * gfe2_0 + 2.0 * ts_xzz_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyzz[i] * gfe2_0 + 8.0 * ts_xyz_yyz[i] * gfe2_0 + 4.0 * ts_xyz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_zz[i] * gfe2_0 + 4.0 * ts_xyzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_yy[i] * gfe2_0 + 4.0 * ts_xyzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_yyzz[i] * gfe_0 + ts_xyzz_yyzz[i] * rgc2_0;

        gr_xyzz_yzzz[i] = 2.0 * ts_yzz_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_zzz[i] * gfe2_0 + 2.0 * ts_xzz_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yzzz[i] * gfe2_0 + 12.0 * ts_xyz_yzz[i] * gfe2_0 + 4.0 * ts_xyz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyzz_yz[i] * gfe2_0 + 6.0 * ts_xyzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_yzzz[i] * gfe_0 + ts_xyzz_yzzz[i] * rgc2_0;

        gr_xyzz_zzzz[i] = 2.0 * ts_yzz_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_zzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zzzz[i] * gfe2_0 + 16.0 * ts_xyz_zzz[i] * gfe2_0 + 4.0 * ts_xyz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xyzz_zz[i] * gfe2_0 + 8.0 * ts_xyzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_zzzz[i] * gfe_0 + ts_xyzz_zzzz[i] * rgc2_0;
    }

    // Set up 135-150 components of targeted buffer : GG

    auto gr_xzzz_xxxx = pbuffer.data(idx_g_gg + 135);

    auto gr_xzzz_xxxy = pbuffer.data(idx_g_gg + 136);

    auto gr_xzzz_xxxz = pbuffer.data(idx_g_gg + 137);

    auto gr_xzzz_xxyy = pbuffer.data(idx_g_gg + 138);

    auto gr_xzzz_xxyz = pbuffer.data(idx_g_gg + 139);

    auto gr_xzzz_xxzz = pbuffer.data(idx_g_gg + 140);

    auto gr_xzzz_xyyy = pbuffer.data(idx_g_gg + 141);

    auto gr_xzzz_xyyz = pbuffer.data(idx_g_gg + 142);

    auto gr_xzzz_xyzz = pbuffer.data(idx_g_gg + 143);

    auto gr_xzzz_xzzz = pbuffer.data(idx_g_gg + 144);

    auto gr_xzzz_yyyy = pbuffer.data(idx_g_gg + 145);

    auto gr_xzzz_yyyz = pbuffer.data(idx_g_gg + 146);

    auto gr_xzzz_yyzz = pbuffer.data(idx_g_gg + 147);

    auto gr_xzzz_yzzz = pbuffer.data(idx_g_gg + 148);

    auto gr_xzzz_zzzz = pbuffer.data(idx_g_gg + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_xxxx, gr_xzzz_xxxy, gr_xzzz_xxxz, gr_xzzz_xxyy, gr_xzzz_xxyz, gr_xzzz_xxzz, gr_xzzz_xyyy, gr_xzzz_xyyz, gr_xzzz_xyzz, gr_xzzz_xzzz, gr_xzzz_yyyy, gr_xzzz_yyyz, gr_xzzz_yyzz, gr_xzzz_yzzz, gr_xzzz_zzzz, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, ts_xzz_xxx, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxy, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyy, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyy, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, ts_xzzz_xx, ts_xzzz_xxx, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxy, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxz, ts_xzzz_xxzz, ts_xzzz_xy, ts_xzzz_xyy, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyz, ts_xzzz_xyzz, ts_xzzz_xz, ts_xzzz_xzz, ts_xzzz_xzzz, ts_xzzz_yy, ts_xzzz_yyy, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyz, ts_xzzz_yyzz, ts_xzzz_yz, ts_xzzz_yzz, ts_xzzz_yzzz, ts_xzzz_zz, ts_xzzz_zzz, ts_xzzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xzzz_xxxx[i] = 8.0 * ts_zzz_xxx[i] * gfe2_0 + 2.0 * ts_zzz_xxxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxxx[i] * gfe2_0 + 6.0 * ts_xzz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xzzz_xx[i] * gfe2_0 + 8.0 * ts_xzzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzzz_xxxx[i] * gfe_0 + ts_xzzz_xxxx[i] * rgc2_0;

        gr_xzzz_xxxy[i] = 6.0 * ts_zzz_xxy[i] * gfe2_0 + 2.0 * ts_zzz_xxxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxxy[i] * gfe2_0 + 6.0 * ts_xzz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzzz_xy[i] * gfe2_0 + 6.0 * ts_xzzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_xxxy[i] * gfe_0 + ts_xzzz_xxxy[i] * rgc2_0;

        gr_xzzz_xxxz[i] = 6.0 * ts_zzz_xxz[i] * gfe2_0 + 2.0 * ts_zzz_xxxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxxz[i] * gfe2_0 + 6.0 * ts_xzz_xxx[i] * gfe2_0 + 6.0 * ts_xzz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzzz_xz[i] * gfe2_0 + 6.0 * ts_xzzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xxxz[i] * gfe_0 + ts_xzzz_xxxz[i] * rgc2_0;

        gr_xzzz_xxyy[i] = 4.0 * ts_zzz_xyy[i] * gfe2_0 + 2.0 * ts_zzz_xxyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxyy[i] * gfe2_0 + 6.0 * ts_xzz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yy[i] * gfe2_0 + 4.0 * ts_xzzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xx[i] * gfe2_0 + 4.0 * ts_xzzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_xxyy[i] * gfe_0 + ts_xzzz_xxyy[i] * rgc2_0;

        gr_xzzz_xxyz[i] = 4.0 * ts_zzz_xyz[i] * gfe2_0 + 2.0 * ts_zzz_xxyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxyz[i] * gfe2_0 + 6.0 * ts_xzz_xxy[i] * gfe2_0 + 6.0 * ts_xzz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yz[i] * gfe2_0 + 4.0 * ts_xzzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xxyz[i] * gfe_0 + ts_xzzz_xxyz[i] * rgc2_0;

        gr_xzzz_xxzz[i] = 4.0 * ts_zzz_xzz[i] * gfe2_0 + 2.0 * ts_zzz_xxzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxzz[i] * gfe2_0 + 12.0 * ts_xzz_xxz[i] * gfe2_0 + 6.0 * ts_xzz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_zz[i] * gfe2_0 + 4.0 * ts_xzzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xx[i] * gfe2_0 + 4.0 * ts_xzzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xxzz[i] * gfe_0 + ts_xzzz_xxzz[i] * rgc2_0;

        gr_xzzz_xyyy[i] = 2.0 * ts_zzz_yyy[i] * gfe2_0 + 2.0 * ts_zzz_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xyyy[i] * gfe2_0 + 6.0 * ts_xzz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzzz_xy[i] * gfe2_0 + 6.0 * ts_xzzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_xyyy[i] * gfe_0 + ts_xzzz_xyyy[i] * rgc2_0;

        gr_xzzz_xyyz[i] = 2.0 * ts_zzz_yyz[i] * gfe2_0 + 2.0 * ts_zzz_xyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xyyz[i] * gfe2_0 + 6.0 * ts_xzz_xyy[i] * gfe2_0 + 6.0 * ts_xzz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xz[i] * gfe2_0 + 4.0 * ts_xzzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xyyz[i] * gfe_0 + ts_xzzz_xyyz[i] * rgc2_0;

        gr_xzzz_xyzz[i] = 2.0 * ts_zzz_yzz[i] * gfe2_0 + 2.0 * ts_zzz_xyzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xyzz[i] * gfe2_0 + 12.0 * ts_xzz_xyz[i] * gfe2_0 + 6.0 * ts_xzz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_xy[i] * gfe2_0 + 4.0 * ts_xzzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xyzz[i] * gfe_0 + ts_xzzz_xyzz[i] * rgc2_0;

        gr_xzzz_xzzz[i] = 2.0 * ts_zzz_zzz[i] * gfe2_0 + 2.0 * ts_zzz_xzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xzzz[i] * gfe2_0 + 18.0 * ts_xzz_xzz[i] * gfe2_0 + 6.0 * ts_xzz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzzz_xz[i] * gfe2_0 + 6.0 * ts_xzzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xzzz[i] * gfe_0 + ts_xzzz_xzzz[i] * rgc2_0;

        gr_xzzz_yyyy[i] = 2.0 * ts_zzz_yyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yyyy[i] * gfe2_0 + 6.0 * ts_xzz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xzzz_yy[i] * gfe2_0 + 8.0 * ts_xzzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_yyyy[i] * gfe_0 + ts_xzzz_yyyy[i] * rgc2_0;

        gr_xzzz_yyyz[i] = 2.0 * ts_zzz_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yyyz[i] * gfe2_0 + 6.0 * ts_xzz_yyy[i] * gfe2_0 + 6.0 * ts_xzz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzzz_yz[i] * gfe2_0 + 6.0 * ts_xzzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_yyyz[i] * gfe_0 + ts_xzzz_yyyz[i] * rgc2_0;

        gr_xzzz_yyzz[i] = 2.0 * ts_zzz_yyzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yyzz[i] * gfe2_0 + 12.0 * ts_xzz_yyz[i] * gfe2_0 + 6.0 * ts_xzz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_zz[i] * gfe2_0 + 4.0 * ts_xzzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_yy[i] * gfe2_0 + 4.0 * ts_xzzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_yyzz[i] * gfe_0 + ts_xzzz_yyzz[i] * rgc2_0;

        gr_xzzz_yzzz[i] = 2.0 * ts_zzz_yzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yzzz[i] * gfe2_0 + 18.0 * ts_xzz_yzz[i] * gfe2_0 + 6.0 * ts_xzz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xzzz_yz[i] * gfe2_0 + 6.0 * ts_xzzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_yzzz[i] * gfe_0 + ts_xzzz_yzzz[i] * rgc2_0;

        gr_xzzz_zzzz[i] = 2.0 * ts_zzz_zzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_zzzz[i] * gfe2_0 + 24.0 * ts_xzz_zzz[i] * gfe2_0 + 6.0 * ts_xzz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xzzz_zz[i] * gfe2_0 + 8.0 * ts_xzzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_zzzz[i] * gfe_0 + ts_xzzz_zzzz[i] * rgc2_0;
    }

    // Set up 150-165 components of targeted buffer : GG

    auto gr_yyyy_xxxx = pbuffer.data(idx_g_gg + 150);

    auto gr_yyyy_xxxy = pbuffer.data(idx_g_gg + 151);

    auto gr_yyyy_xxxz = pbuffer.data(idx_g_gg + 152);

    auto gr_yyyy_xxyy = pbuffer.data(idx_g_gg + 153);

    auto gr_yyyy_xxyz = pbuffer.data(idx_g_gg + 154);

    auto gr_yyyy_xxzz = pbuffer.data(idx_g_gg + 155);

    auto gr_yyyy_xyyy = pbuffer.data(idx_g_gg + 156);

    auto gr_yyyy_xyyz = pbuffer.data(idx_g_gg + 157);

    auto gr_yyyy_xyzz = pbuffer.data(idx_g_gg + 158);

    auto gr_yyyy_xzzz = pbuffer.data(idx_g_gg + 159);

    auto gr_yyyy_yyyy = pbuffer.data(idx_g_gg + 160);

    auto gr_yyyy_yyyz = pbuffer.data(idx_g_gg + 161);

    auto gr_yyyy_yyzz = pbuffer.data(idx_g_gg + 162);

    auto gr_yyyy_yzzz = pbuffer.data(idx_g_gg + 163);

    auto gr_yyyy_zzzz = pbuffer.data(idx_g_gg + 164);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_xxxx, gr_yyyy_xxxy, gr_yyyy_xxxz, gr_yyyy_xxyy, gr_yyyy_xxyz, gr_yyyy_xxzz, gr_yyyy_xyyy, gr_yyyy_xyyz, gr_yyyy_xyzz, gr_yyyy_xzzz, gr_yyyy_yyyy, gr_yyyy_yyyz, gr_yyyy_yyzz, gr_yyyy_yzzz, gr_yyyy_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, ts_yyyy_xx, ts_yyyy_xxx, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxz, ts_yyyy_xxzz, ts_yyyy_xy, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_xz, ts_yyyy_xzz, ts_yyyy_xzzz, ts_yyyy_yy, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyy_zz, ts_yyyy_zzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyyy_xxxx[i] = 12.0 * ts_yy_xxxx[i] * gfe2_0 + 8.0 * ts_yyy_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_yyyy_xx[i] * gfe2_0 + 8.0 * ts_yyyy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyy_xxxx[i] * gfe_0 + ts_yyyy_xxxx[i] * rgc2_0;

        gr_yyyy_xxxy[i] = 12.0 * ts_yy_xxxy[i] * gfe2_0 + 8.0 * ts_yyy_xxx[i] * gfe2_0 + 8.0 * ts_yyy_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_xy[i] * gfe2_0 + 6.0 * ts_yyyy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_xxxy[i] * gfe_0 + ts_yyyy_xxxy[i] * rgc2_0;

        gr_yyyy_xxxz[i] = 12.0 * ts_yy_xxxz[i] * gfe2_0 + 8.0 * ts_yyy_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_xz[i] * gfe2_0 + 6.0 * ts_yyyy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xxxz[i] * gfe_0 + ts_yyyy_xxxz[i] * rgc2_0;

        gr_yyyy_xxyy[i] = 12.0 * ts_yy_xxyy[i] * gfe2_0 + 16.0 * ts_yyy_xxy[i] * gfe2_0 + 8.0 * ts_yyy_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yy[i] * gfe2_0 + 4.0 * ts_yyyy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xx[i] * gfe2_0 + 4.0 * ts_yyyy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_xxyy[i] * gfe_0 + ts_yyyy_xxyy[i] * rgc2_0;

        gr_yyyy_xxyz[i] = 12.0 * ts_yy_xxyz[i] * gfe2_0 + 8.0 * ts_yyy_xxz[i] * gfe2_0 + 8.0 * ts_yyy_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yz[i] * gfe2_0 + 4.0 * ts_yyyy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xxyz[i] * gfe_0 + ts_yyyy_xxyz[i] * rgc2_0;

        gr_yyyy_xxzz[i] = 12.0 * ts_yy_xxzz[i] * gfe2_0 + 8.0 * ts_yyy_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_zz[i] * gfe2_0 + 4.0 * ts_yyyy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xx[i] * gfe2_0 + 4.0 * ts_yyyy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xxzz[i] * gfe_0 + ts_yyyy_xxzz[i] * rgc2_0;

        gr_yyyy_xyyy[i] = 12.0 * ts_yy_xyyy[i] * gfe2_0 + 24.0 * ts_yyy_xyy[i] * gfe2_0 + 8.0 * ts_yyy_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyyy_xy[i] * gfe2_0 + 6.0 * ts_yyyy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_xyyy[i] * gfe_0 + ts_yyyy_xyyy[i] * rgc2_0;

        gr_yyyy_xyyz[i] = 12.0 * ts_yy_xyyz[i] * gfe2_0 + 16.0 * ts_yyy_xyz[i] * gfe2_0 + 8.0 * ts_yyy_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xz[i] * gfe2_0 + 4.0 * ts_yyyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xyyz[i] * gfe_0 + ts_yyyy_xyyz[i] * rgc2_0;

        gr_yyyy_xyzz[i] = 12.0 * ts_yy_xyzz[i] * gfe2_0 + 8.0 * ts_yyy_xzz[i] * gfe2_0 + 8.0 * ts_yyy_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_xy[i] * gfe2_0 + 4.0 * ts_yyyy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xyzz[i] * gfe_0 + ts_yyyy_xyzz[i] * rgc2_0;

        gr_yyyy_xzzz[i] = 12.0 * ts_yy_xzzz[i] * gfe2_0 + 8.0 * ts_yyy_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyyy_xz[i] * gfe2_0 + 6.0 * ts_yyyy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xzzz[i] * gfe_0 + ts_yyyy_xzzz[i] * rgc2_0;

        gr_yyyy_yyyy[i] = 12.0 * ts_yy_yyyy[i] * gfe2_0 + 32.0 * ts_yyy_yyy[i] * gfe2_0 + 8.0 * ts_yyy_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_yyyy_yy[i] * gfe2_0 + 8.0 * ts_yyyy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_yyyy[i] * gfe_0 + ts_yyyy_yyyy[i] * rgc2_0;

        gr_yyyy_yyyz[i] = 12.0 * ts_yy_yyyz[i] * gfe2_0 + 24.0 * ts_yyy_yyz[i] * gfe2_0 + 8.0 * ts_yyy_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_yz[i] * gfe2_0 + 6.0 * ts_yyyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_yyyz[i] * gfe_0 + ts_yyyy_yyyz[i] * rgc2_0;

        gr_yyyy_yyzz[i] = 12.0 * ts_yy_yyzz[i] * gfe2_0 + 16.0 * ts_yyy_yzz[i] * gfe2_0 + 8.0 * ts_yyy_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_zz[i] * gfe2_0 + 4.0 * ts_yyyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yy[i] * gfe2_0 + 4.0 * ts_yyyy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_yyzz[i] * gfe_0 + ts_yyyy_yyzz[i] * rgc2_0;

        gr_yyyy_yzzz[i] = 12.0 * ts_yy_yzzz[i] * gfe2_0 + 8.0 * ts_yyy_zzz[i] * gfe2_0 + 8.0 * ts_yyy_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_yz[i] * gfe2_0 + 6.0 * ts_yyyy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_yzzz[i] * gfe_0 + ts_yyyy_yzzz[i] * rgc2_0;

        gr_yyyy_zzzz[i] = 12.0 * ts_yy_zzzz[i] * gfe2_0 + 8.0 * ts_yyy_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_yyyy_zz[i] * gfe2_0 + 8.0 * ts_yyyy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_zzzz[i] * gfe_0 + ts_yyyy_zzzz[i] * rgc2_0;
    }

    // Set up 165-180 components of targeted buffer : GG

    auto gr_yyyz_xxxx = pbuffer.data(idx_g_gg + 165);

    auto gr_yyyz_xxxy = pbuffer.data(idx_g_gg + 166);

    auto gr_yyyz_xxxz = pbuffer.data(idx_g_gg + 167);

    auto gr_yyyz_xxyy = pbuffer.data(idx_g_gg + 168);

    auto gr_yyyz_xxyz = pbuffer.data(idx_g_gg + 169);

    auto gr_yyyz_xxzz = pbuffer.data(idx_g_gg + 170);

    auto gr_yyyz_xyyy = pbuffer.data(idx_g_gg + 171);

    auto gr_yyyz_xyyz = pbuffer.data(idx_g_gg + 172);

    auto gr_yyyz_xyzz = pbuffer.data(idx_g_gg + 173);

    auto gr_yyyz_xzzz = pbuffer.data(idx_g_gg + 174);

    auto gr_yyyz_yyyy = pbuffer.data(idx_g_gg + 175);

    auto gr_yyyz_yyyz = pbuffer.data(idx_g_gg + 176);

    auto gr_yyyz_yyzz = pbuffer.data(idx_g_gg + 177);

    auto gr_yyyz_yzzz = pbuffer.data(idx_g_gg + 178);

    auto gr_yyyz_zzzz = pbuffer.data(idx_g_gg + 179);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_xxxx, gr_yyyz_xxxy, gr_yyyz_xxxz, gr_yyyz_xxyy, gr_yyyz_xxyz, gr_yyyz_xxzz, gr_yyyz_xyyy, gr_yyyz_xyyz, gr_yyyz_xyzz, gr_yyyz_xzzz, gr_yyyz_yyyy, gr_yyyz_yyyz, gr_yyyz_yyzz, gr_yyyz_yzzz, gr_yyyz_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, ts_yyyz_xx, ts_yyyz_xxx, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxy, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxz, ts_yyyz_xxzz, ts_yyyz_xy, ts_yyyz_xyy, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyz, ts_yyyz_xyzz, ts_yyyz_xz, ts_yyyz_xzz, ts_yyyz_xzzz, ts_yyyz_yy, ts_yyyz_yyy, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyz, ts_yyyz_yyzz, ts_yyyz_yz, ts_yyyz_yzz, ts_yyyz_yzzz, ts_yyyz_zz, ts_yyyz_zzz, ts_yyyz_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyyz_xxxx[i] = 6.0 * ts_yz_xxxx[i] * gfe2_0 + 6.0 * ts_yyz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyyz_xx[i] * gfe2_0 + 8.0 * ts_yyyz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyz_xxxx[i] * gfe_0 + ts_yyyz_xxxx[i] * rgc2_0;

        gr_yyyz_xxxy[i] = 6.0 * ts_yz_xxxy[i] * gfe2_0 + 6.0 * ts_yyz_xxx[i] * gfe2_0 + 6.0 * ts_yyz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyyz_xy[i] * gfe2_0 + 6.0 * ts_yyyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_xxxy[i] * gfe_0 + ts_yyyz_xxxy[i] * rgc2_0;

        gr_yyyz_xxxz[i] = 6.0 * ts_yz_xxxz[i] * gfe2_0 + 6.0 * ts_yyz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxx[i] * gfe2_0 + 2.0 * ts_yyy_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyyz_xz[i] * gfe2_0 + 6.0 * ts_yyyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xxxz[i] * gfe_0 + ts_yyyz_xxxz[i] * rgc2_0;

        gr_yyyz_xxyy[i] = 6.0 * ts_yz_xxyy[i] * gfe2_0 + 12.0 * ts_yyz_xxy[i] * gfe2_0 + 6.0 * ts_yyz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yy[i] * gfe2_0 + 4.0 * ts_yyyz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xx[i] * gfe2_0 + 4.0 * ts_yyyz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_xxyy[i] * gfe_0 + ts_yyyz_xxyy[i] * rgc2_0;

        gr_yyyz_xxyz[i] = 6.0 * ts_yz_xxyz[i] * gfe2_0 + 6.0 * ts_yyz_xxz[i] * gfe2_0 + 6.0 * ts_yyz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxy[i] * gfe2_0 + 2.0 * ts_yyy_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yz[i] * gfe2_0 + 4.0 * ts_yyyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xxyz[i] * gfe_0 + ts_yyyz_xxyz[i] * rgc2_0;

        gr_yyyz_xxzz[i] = 6.0 * ts_yz_xxzz[i] * gfe2_0 + 6.0 * ts_yyz_xxzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yyy_xxz[i] * gfe2_0 + 2.0 * ts_yyy_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_zz[i] * gfe2_0 + 4.0 * ts_yyyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xx[i] * gfe2_0 + 4.0 * ts_yyyz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xxzz[i] * gfe_0 + ts_yyyz_xxzz[i] * rgc2_0;

        gr_yyyz_xyyy[i] = 6.0 * ts_yz_xyyy[i] * gfe2_0 + 18.0 * ts_yyz_xyy[i] * gfe2_0 + 6.0 * ts_yyz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyyz_xy[i] * gfe2_0 + 6.0 * ts_yyyz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_xyyy[i] * gfe_0 + ts_yyyz_xyyy[i] * rgc2_0;

        gr_yyyz_xyyz[i] = 6.0 * ts_yz_xyyz[i] * gfe2_0 + 12.0 * ts_yyz_xyz[i] * gfe2_0 + 6.0 * ts_yyz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xyy[i] * gfe2_0 + 2.0 * ts_yyy_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xz[i] * gfe2_0 + 4.0 * ts_yyyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xyyz[i] * gfe_0 + ts_yyyz_xyyz[i] * rgc2_0;

        gr_yyyz_xyzz[i] = 6.0 * ts_yz_xyzz[i] * gfe2_0 + 6.0 * ts_yyz_xzz[i] * gfe2_0 + 6.0 * ts_yyz_xyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yyy_xyz[i] * gfe2_0 + 2.0 * ts_yyy_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_xy[i] * gfe2_0 + 4.0 * ts_yyyz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xyzz[i] * gfe_0 + ts_yyyz_xyzz[i] * rgc2_0;

        gr_yyyz_xzzz[i] = 6.0 * ts_yz_xzzz[i] * gfe2_0 + 6.0 * ts_yyz_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_xzz[i] * gfe2_0 + 2.0 * ts_yyy_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyyz_xz[i] * gfe2_0 + 6.0 * ts_yyyz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xzzz[i] * gfe_0 + ts_yyyz_xzzz[i] * rgc2_0;

        gr_yyyz_yyyy[i] = 6.0 * ts_yz_yyyy[i] * gfe2_0 + 24.0 * ts_yyz_yyy[i] * gfe2_0 + 6.0 * ts_yyz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyyz_yy[i] * gfe2_0 + 8.0 * ts_yyyz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_yyyy[i] * gfe_0 + ts_yyyz_yyyy[i] * rgc2_0;

        gr_yyyz_yyyz[i] = 6.0 * ts_yz_yyyz[i] * gfe2_0 + 18.0 * ts_yyz_yyz[i] * gfe2_0 + 6.0 * ts_yyz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yyy[i] * gfe2_0 + 2.0 * ts_yyy_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyyz_yz[i] * gfe2_0 + 6.0 * ts_yyyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_yyyz[i] * gfe_0 + ts_yyyz_yyyz[i] * rgc2_0;

        gr_yyyz_yyzz[i] = 6.0 * ts_yz_yyzz[i] * gfe2_0 + 12.0 * ts_yyz_yzz[i] * gfe2_0 + 6.0 * ts_yyz_yyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yyy_yyz[i] * gfe2_0 + 2.0 * ts_yyy_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_zz[i] * gfe2_0 + 4.0 * ts_yyyz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_yy[i] * gfe2_0 + 4.0 * ts_yyyz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_yyzz[i] * gfe_0 + ts_yyyz_yyzz[i] * rgc2_0;

        gr_yyyz_yzzz[i] = 6.0 * ts_yz_yzzz[i] * gfe2_0 + 6.0 * ts_yyz_zzz[i] * gfe2_0 + 6.0 * ts_yyz_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_yzz[i] * gfe2_0 + 2.0 * ts_yyy_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyz_yz[i] * gfe2_0 + 6.0 * ts_yyyz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_yzzz[i] * gfe_0 + ts_yyyz_yzzz[i] * rgc2_0;

        gr_yyyz_zzzz[i] = 6.0 * ts_yz_zzzz[i] * gfe2_0 + 6.0 * ts_yyz_zzzz[i] * gfe_0 * gc_y[i] + 8.0 * ts_yyy_zzz[i] * gfe2_0 + 2.0 * ts_yyy_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyyz_zz[i] * gfe2_0 + 8.0 * ts_yyyz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_zzzz[i] * gfe_0 + ts_yyyz_zzzz[i] * rgc2_0;
    }

    // Set up 180-195 components of targeted buffer : GG

    auto gr_yyzz_xxxx = pbuffer.data(idx_g_gg + 180);

    auto gr_yyzz_xxxy = pbuffer.data(idx_g_gg + 181);

    auto gr_yyzz_xxxz = pbuffer.data(idx_g_gg + 182);

    auto gr_yyzz_xxyy = pbuffer.data(idx_g_gg + 183);

    auto gr_yyzz_xxyz = pbuffer.data(idx_g_gg + 184);

    auto gr_yyzz_xxzz = pbuffer.data(idx_g_gg + 185);

    auto gr_yyzz_xyyy = pbuffer.data(idx_g_gg + 186);

    auto gr_yyzz_xyyz = pbuffer.data(idx_g_gg + 187);

    auto gr_yyzz_xyzz = pbuffer.data(idx_g_gg + 188);

    auto gr_yyzz_xzzz = pbuffer.data(idx_g_gg + 189);

    auto gr_yyzz_yyyy = pbuffer.data(idx_g_gg + 190);

    auto gr_yyzz_yyyz = pbuffer.data(idx_g_gg + 191);

    auto gr_yyzz_yyzz = pbuffer.data(idx_g_gg + 192);

    auto gr_yyzz_yzzz = pbuffer.data(idx_g_gg + 193);

    auto gr_yyzz_zzzz = pbuffer.data(idx_g_gg + 194);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_xxxx, gr_yyzz_xxxy, gr_yyzz_xxxz, gr_yyzz_xxyy, gr_yyzz_xxyz, gr_yyzz_xxzz, gr_yyzz_xyyy, gr_yyzz_xyyz, gr_yyzz_xyzz, gr_yyzz_xzzz, gr_yyzz_yyyy, gr_yyzz_yyyz, gr_yyzz_yyzz, gr_yyzz_yzzz, gr_yyzz_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyz_xxx, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxy, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyy, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyy, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, ts_yyzz_xx, ts_yyzz_xxx, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxy, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxz, ts_yyzz_xxzz, ts_yyzz_xy, ts_yyzz_xyy, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyz, ts_yyzz_xyzz, ts_yyzz_xz, ts_yyzz_xzz, ts_yyzz_xzzz, ts_yyzz_yy, ts_yyzz_yyy, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyz, ts_yyzz_yyzz, ts_yyzz_yz, ts_yyzz_yzz, ts_yyzz_yzzz, ts_yyzz_zz, ts_yyzz_zzz, ts_yyzz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyzz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gfe2_0 + 4.0 * ts_yzz_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxxx[i] * gfe2_0 + 4.0 * ts_yyz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyzz_xx[i] * gfe2_0 + 8.0 * ts_yyzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyzz_xxxx[i] * gfe_0 + ts_yyzz_xxxx[i] * rgc2_0;

        gr_yyzz_xxxy[i] = 2.0 * ts_zz_xxxy[i] * gfe2_0 + 4.0 * ts_yzz_xxx[i] * gfe2_0 + 4.0 * ts_yzz_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxxy[i] * gfe2_0 + 4.0 * ts_yyz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyzz_xy[i] * gfe2_0 + 6.0 * ts_yyzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_xxxy[i] * gfe_0 + ts_yyzz_xxxy[i] * rgc2_0;

        gr_yyzz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gfe2_0 + 4.0 * ts_yzz_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxxz[i] * gfe2_0 + 4.0 * ts_yyz_xxx[i] * gfe2_0 + 4.0 * ts_yyz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyzz_xz[i] * gfe2_0 + 6.0 * ts_yyzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xxxz[i] * gfe_0 + ts_yyzz_xxxz[i] * rgc2_0;

        gr_yyzz_xxyy[i] = 2.0 * ts_zz_xxyy[i] * gfe2_0 + 8.0 * ts_yzz_xxy[i] * gfe2_0 + 4.0 * ts_yzz_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxyy[i] * gfe2_0 + 4.0 * ts_yyz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yy[i] * gfe2_0 + 4.0 * ts_yyzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xx[i] * gfe2_0 + 4.0 * ts_yyzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_xxyy[i] * gfe_0 + ts_yyzz_xxyy[i] * rgc2_0;

        gr_yyzz_xxyz[i] = 2.0 * ts_zz_xxyz[i] * gfe2_0 + 4.0 * ts_yzz_xxz[i] * gfe2_0 + 4.0 * ts_yzz_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxyz[i] * gfe2_0 + 4.0 * ts_yyz_xxy[i] * gfe2_0 + 4.0 * ts_yyz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yz[i] * gfe2_0 + 4.0 * ts_yyzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xxyz[i] * gfe_0 + ts_yyzz_xxyz[i] * rgc2_0;

        gr_yyzz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gfe2_0 + 4.0 * ts_yzz_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxzz[i] * gfe2_0 + 8.0 * ts_yyz_xxz[i] * gfe2_0 + 4.0 * ts_yyz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_zz[i] * gfe2_0 + 4.0 * ts_yyzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xx[i] * gfe2_0 + 4.0 * ts_yyzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xxzz[i] * gfe_0 + ts_yyzz_xxzz[i] * rgc2_0;

        gr_yyzz_xyyy[i] = 2.0 * ts_zz_xyyy[i] * gfe2_0 + 12.0 * ts_yzz_xyy[i] * gfe2_0 + 4.0 * ts_yzz_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyyy[i] * gfe2_0 + 4.0 * ts_yyz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyzz_xy[i] * gfe2_0 + 6.0 * ts_yyzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_xyyy[i] * gfe_0 + ts_yyzz_xyyy[i] * rgc2_0;

        gr_yyzz_xyyz[i] = 2.0 * ts_zz_xyyz[i] * gfe2_0 + 8.0 * ts_yzz_xyz[i] * gfe2_0 + 4.0 * ts_yzz_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyyz[i] * gfe2_0 + 4.0 * ts_yyz_xyy[i] * gfe2_0 + 4.0 * ts_yyz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xz[i] * gfe2_0 + 4.0 * ts_yyzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xyyz[i] * gfe_0 + ts_yyzz_xyyz[i] * rgc2_0;

        gr_yyzz_xyzz[i] = 2.0 * ts_zz_xyzz[i] * gfe2_0 + 4.0 * ts_yzz_xzz[i] * gfe2_0 + 4.0 * ts_yzz_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyzz[i] * gfe2_0 + 8.0 * ts_yyz_xyz[i] * gfe2_0 + 4.0 * ts_yyz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_xy[i] * gfe2_0 + 4.0 * ts_yyzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xyzz[i] * gfe_0 + ts_yyzz_xyzz[i] * rgc2_0;

        gr_yyzz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe2_0 + 4.0 * ts_yzz_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xzzz[i] * gfe2_0 + 12.0 * ts_yyz_xzz[i] * gfe2_0 + 4.0 * ts_yyz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yyzz_xz[i] * gfe2_0 + 6.0 * ts_yyzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xzzz[i] * gfe_0 + ts_yyzz_xzzz[i] * rgc2_0;

        gr_yyzz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe2_0 + 16.0 * ts_yzz_yyy[i] * gfe2_0 + 4.0 * ts_yzz_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyyy[i] * gfe2_0 + 4.0 * ts_yyz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyzz_yy[i] * gfe2_0 + 8.0 * ts_yyzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_yyyy[i] * gfe_0 + ts_yyzz_yyyy[i] * rgc2_0;

        gr_yyzz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe2_0 + 12.0 * ts_yzz_yyz[i] * gfe2_0 + 4.0 * ts_yzz_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyyz[i] * gfe2_0 + 4.0 * ts_yyz_yyy[i] * gfe2_0 + 4.0 * ts_yyz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyzz_yz[i] * gfe2_0 + 6.0 * ts_yyzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_yyyz[i] * gfe_0 + ts_yyzz_yyyz[i] * rgc2_0;

        gr_yyzz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe2_0 + 8.0 * ts_yzz_yzz[i] * gfe2_0 + 4.0 * ts_yzz_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyzz[i] * gfe2_0 + 8.0 * ts_yyz_yyz[i] * gfe2_0 + 4.0 * ts_yyz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_zz[i] * gfe2_0 + 4.0 * ts_yyzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_yy[i] * gfe2_0 + 4.0 * ts_yyzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_yyzz[i] * gfe_0 + ts_yyzz_yyzz[i] * rgc2_0;

        gr_yyzz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe2_0 + 4.0 * ts_yzz_zzz[i] * gfe2_0 + 4.0 * ts_yzz_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yzzz[i] * gfe2_0 + 12.0 * ts_yyz_yzz[i] * gfe2_0 + 4.0 * ts_yyz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyzz_yz[i] * gfe2_0 + 6.0 * ts_yyzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_yzzz[i] * gfe_0 + ts_yyzz_yzzz[i] * rgc2_0;

        gr_yyzz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe2_0 + 4.0 * ts_yzz_zzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zzzz[i] * gfe2_0 + 16.0 * ts_yyz_zzz[i] * gfe2_0 + 4.0 * ts_yyz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_yyzz_zz[i] * gfe2_0 + 8.0 * ts_yyzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_zzzz[i] * gfe_0 + ts_yyzz_zzzz[i] * rgc2_0;
    }

    // Set up 195-210 components of targeted buffer : GG

    auto gr_yzzz_xxxx = pbuffer.data(idx_g_gg + 195);

    auto gr_yzzz_xxxy = pbuffer.data(idx_g_gg + 196);

    auto gr_yzzz_xxxz = pbuffer.data(idx_g_gg + 197);

    auto gr_yzzz_xxyy = pbuffer.data(idx_g_gg + 198);

    auto gr_yzzz_xxyz = pbuffer.data(idx_g_gg + 199);

    auto gr_yzzz_xxzz = pbuffer.data(idx_g_gg + 200);

    auto gr_yzzz_xyyy = pbuffer.data(idx_g_gg + 201);

    auto gr_yzzz_xyyz = pbuffer.data(idx_g_gg + 202);

    auto gr_yzzz_xyzz = pbuffer.data(idx_g_gg + 203);

    auto gr_yzzz_xzzz = pbuffer.data(idx_g_gg + 204);

    auto gr_yzzz_yyyy = pbuffer.data(idx_g_gg + 205);

    auto gr_yzzz_yyyz = pbuffer.data(idx_g_gg + 206);

    auto gr_yzzz_yyzz = pbuffer.data(idx_g_gg + 207);

    auto gr_yzzz_yzzz = pbuffer.data(idx_g_gg + 208);

    auto gr_yzzz_zzzz = pbuffer.data(idx_g_gg + 209);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_xxxx, gr_yzzz_xxxy, gr_yzzz_xxxz, gr_yzzz_xxyy, gr_yzzz_xxyz, gr_yzzz_xxzz, gr_yzzz_xyyy, gr_yzzz_xyyz, gr_yzzz_xyzz, gr_yzzz_xzzz, gr_yzzz_yyyy, gr_yzzz_yyyz, gr_yzzz_yyzz, gr_yzzz_yzzz, gr_yzzz_zzzz, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, ts_yzz_xxx, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, ts_yzzz_xx, ts_yzzz_xxx, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxy, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxz, ts_yzzz_xxzz, ts_yzzz_xy, ts_yzzz_xyy, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyz, ts_yzzz_xyzz, ts_yzzz_xz, ts_yzzz_xzz, ts_yzzz_xzzz, ts_yzzz_yy, ts_yzzz_yyy, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyz, ts_yzzz_yyzz, ts_yzzz_yz, ts_yzzz_yzz, ts_yzzz_yzzz, ts_yzzz_zz, ts_yzzz_zzz, ts_yzzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yzzz_xxxx[i] = 2.0 * ts_zzz_xxxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxxx[i] * gfe2_0 + 6.0 * ts_yzz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_yzzz_xx[i] * gfe2_0 + 8.0 * ts_yzzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzzz_xxxx[i] * gfe_0 + ts_yzzz_xxxx[i] * rgc2_0;

        gr_yzzz_xxxy[i] = 2.0 * ts_zzz_xxx[i] * gfe2_0 + 2.0 * ts_zzz_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxxy[i] * gfe2_0 + 6.0 * ts_yzz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzzz_xy[i] * gfe2_0 + 6.0 * ts_yzzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_xxxy[i] * gfe_0 + ts_yzzz_xxxy[i] * rgc2_0;

        gr_yzzz_xxxz[i] = 2.0 * ts_zzz_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxxz[i] * gfe2_0 + 6.0 * ts_yzz_xxx[i] * gfe2_0 + 6.0 * ts_yzz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzzz_xz[i] * gfe2_0 + 6.0 * ts_yzzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xxxz[i] * gfe_0 + ts_yzzz_xxxz[i] * rgc2_0;

        gr_yzzz_xxyy[i] = 4.0 * ts_zzz_xxy[i] * gfe2_0 + 2.0 * ts_zzz_xxyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxyy[i] * gfe2_0 + 6.0 * ts_yzz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yy[i] * gfe2_0 + 4.0 * ts_yzzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xx[i] * gfe2_0 + 4.0 * ts_yzzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_xxyy[i] * gfe_0 + ts_yzzz_xxyy[i] * rgc2_0;

        gr_yzzz_xxyz[i] = 2.0 * ts_zzz_xxz[i] * gfe2_0 + 2.0 * ts_zzz_xxyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxyz[i] * gfe2_0 + 6.0 * ts_yzz_xxy[i] * gfe2_0 + 6.0 * ts_yzz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yz[i] * gfe2_0 + 4.0 * ts_yzzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xxyz[i] * gfe_0 + ts_yzzz_xxyz[i] * rgc2_0;

        gr_yzzz_xxzz[i] = 2.0 * ts_zzz_xxzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxzz[i] * gfe2_0 + 12.0 * ts_yzz_xxz[i] * gfe2_0 + 6.0 * ts_yzz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_zz[i] * gfe2_0 + 4.0 * ts_yzzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xx[i] * gfe2_0 + 4.0 * ts_yzzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xxzz[i] * gfe_0 + ts_yzzz_xxzz[i] * rgc2_0;

        gr_yzzz_xyyy[i] = 6.0 * ts_zzz_xyy[i] * gfe2_0 + 2.0 * ts_zzz_xyyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xyyy[i] * gfe2_0 + 6.0 * ts_yzz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yzzz_xy[i] * gfe2_0 + 6.0 * ts_yzzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_xyyy[i] * gfe_0 + ts_yzzz_xyyy[i] * rgc2_0;

        gr_yzzz_xyyz[i] = 4.0 * ts_zzz_xyz[i] * gfe2_0 + 2.0 * ts_zzz_xyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xyyz[i] * gfe2_0 + 6.0 * ts_yzz_xyy[i] * gfe2_0 + 6.0 * ts_yzz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xz[i] * gfe2_0 + 4.0 * ts_yzzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xyyz[i] * gfe_0 + ts_yzzz_xyyz[i] * rgc2_0;

        gr_yzzz_xyzz[i] = 2.0 * ts_zzz_xzz[i] * gfe2_0 + 2.0 * ts_zzz_xyzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xyzz[i] * gfe2_0 + 12.0 * ts_yzz_xyz[i] * gfe2_0 + 6.0 * ts_yzz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_xy[i] * gfe2_0 + 4.0 * ts_yzzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xyzz[i] * gfe_0 + ts_yzzz_xyzz[i] * rgc2_0;

        gr_yzzz_xzzz[i] = 2.0 * ts_zzz_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xzzz[i] * gfe2_0 + 18.0 * ts_yzz_xzz[i] * gfe2_0 + 6.0 * ts_yzz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yzzz_xz[i] * gfe2_0 + 6.0 * ts_yzzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xzzz[i] * gfe_0 + ts_yzzz_xzzz[i] * rgc2_0;

        gr_yzzz_yyyy[i] = 8.0 * ts_zzz_yyy[i] * gfe2_0 + 2.0 * ts_zzz_yyyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yyyy[i] * gfe2_0 + 6.0 * ts_yzz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_yzzz_yy[i] * gfe2_0 + 8.0 * ts_yzzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_yyyy[i] * gfe_0 + ts_yzzz_yyyy[i] * rgc2_0;

        gr_yzzz_yyyz[i] = 6.0 * ts_zzz_yyz[i] * gfe2_0 + 2.0 * ts_zzz_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yyyz[i] * gfe2_0 + 6.0 * ts_yzz_yyy[i] * gfe2_0 + 6.0 * ts_yzz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzzz_yz[i] * gfe2_0 + 6.0 * ts_yzzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_yyyz[i] * gfe_0 + ts_yzzz_yyyz[i] * rgc2_0;

        gr_yzzz_yyzz[i] = 4.0 * ts_zzz_yzz[i] * gfe2_0 + 2.0 * ts_zzz_yyzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yyzz[i] * gfe2_0 + 12.0 * ts_yzz_yyz[i] * gfe2_0 + 6.0 * ts_yzz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_zz[i] * gfe2_0 + 4.0 * ts_yzzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_yy[i] * gfe2_0 + 4.0 * ts_yzzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_yyzz[i] * gfe_0 + ts_yzzz_yyzz[i] * rgc2_0;

        gr_yzzz_yzzz[i] = 2.0 * ts_zzz_zzz[i] * gfe2_0 + 2.0 * ts_zzz_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yzzz[i] * gfe2_0 + 18.0 * ts_yzz_yzz[i] * gfe2_0 + 6.0 * ts_yzz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yzzz_yz[i] * gfe2_0 + 6.0 * ts_yzzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_yzzz[i] * gfe_0 + ts_yzzz_yzzz[i] * rgc2_0;

        gr_yzzz_zzzz[i] = 2.0 * ts_zzz_zzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_zzzz[i] * gfe2_0 + 24.0 * ts_yzz_zzz[i] * gfe2_0 + 6.0 * ts_yzz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_yzzz_zz[i] * gfe2_0 + 8.0 * ts_yzzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_zzzz[i] * gfe_0 + ts_yzzz_zzzz[i] * rgc2_0;
    }

    // Set up 210-225 components of targeted buffer : GG

    auto gr_zzzz_xxxx = pbuffer.data(idx_g_gg + 210);

    auto gr_zzzz_xxxy = pbuffer.data(idx_g_gg + 211);

    auto gr_zzzz_xxxz = pbuffer.data(idx_g_gg + 212);

    auto gr_zzzz_xxyy = pbuffer.data(idx_g_gg + 213);

    auto gr_zzzz_xxyz = pbuffer.data(idx_g_gg + 214);

    auto gr_zzzz_xxzz = pbuffer.data(idx_g_gg + 215);

    auto gr_zzzz_xyyy = pbuffer.data(idx_g_gg + 216);

    auto gr_zzzz_xyyz = pbuffer.data(idx_g_gg + 217);

    auto gr_zzzz_xyzz = pbuffer.data(idx_g_gg + 218);

    auto gr_zzzz_xzzz = pbuffer.data(idx_g_gg + 219);

    auto gr_zzzz_yyyy = pbuffer.data(idx_g_gg + 220);

    auto gr_zzzz_yyyz = pbuffer.data(idx_g_gg + 221);

    auto gr_zzzz_yyzz = pbuffer.data(idx_g_gg + 222);

    auto gr_zzzz_yzzz = pbuffer.data(idx_g_gg + 223);

    auto gr_zzzz_zzzz = pbuffer.data(idx_g_gg + 224);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_xxxx, gr_zzzz_xxxy, gr_zzzz_xxxz, gr_zzzz_xxyy, gr_zzzz_xxyz, gr_zzzz_xxzz, gr_zzzz_xyyy, gr_zzzz_xyyz, gr_zzzz_xyzz, gr_zzzz_xzzz, gr_zzzz_yyyy, gr_zzzz_yyyz, gr_zzzz_yyzz, gr_zzzz_yzzz, gr_zzzz_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, ts_zzzz_xx, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xy, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yy, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zz, ts_zzzz_zzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_zzzz_xxxx[i] = 12.0 * ts_zz_xxxx[i] * gfe2_0 + 8.0 * ts_zzz_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_zzzz_xx[i] * gfe2_0 + 8.0 * ts_zzzz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzzz_xxxx[i] * gfe_0 + ts_zzzz_xxxx[i] * rgc2_0;

        gr_zzzz_xxxy[i] = 12.0 * ts_zz_xxxy[i] * gfe2_0 + 8.0 * ts_zzz_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzzz_xy[i] * gfe2_0 + 6.0 * ts_zzzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_xxxy[i] * gfe_0 + ts_zzzz_xxxy[i] * rgc2_0;

        gr_zzzz_xxxz[i] = 12.0 * ts_zz_xxxz[i] * gfe2_0 + 8.0 * ts_zzz_xxx[i] * gfe2_0 + 8.0 * ts_zzz_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzzz_xz[i] * gfe2_0 + 6.0 * ts_zzzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xxxz[i] * gfe_0 + ts_zzzz_xxxz[i] * rgc2_0;

        gr_zzzz_xxyy[i] = 12.0 * ts_zz_xxyy[i] * gfe2_0 + 8.0 * ts_zzz_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yy[i] * gfe2_0 + 4.0 * ts_zzzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xx[i] * gfe2_0 + 4.0 * ts_zzzz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_xxyy[i] * gfe_0 + ts_zzzz_xxyy[i] * rgc2_0;

        gr_zzzz_xxyz[i] = 12.0 * ts_zz_xxyz[i] * gfe2_0 + 8.0 * ts_zzz_xxy[i] * gfe2_0 + 8.0 * ts_zzz_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yz[i] * gfe2_0 + 4.0 * ts_zzzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xxyz[i] * gfe_0 + ts_zzzz_xxyz[i] * rgc2_0;

        gr_zzzz_xxzz[i] = 12.0 * ts_zz_xxzz[i] * gfe2_0 + 16.0 * ts_zzz_xxz[i] * gfe2_0 + 8.0 * ts_zzz_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_zz[i] * gfe2_0 + 4.0 * ts_zzzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xx[i] * gfe2_0 + 4.0 * ts_zzzz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xxzz[i] * gfe_0 + ts_zzzz_xxzz[i] * rgc2_0;

        gr_zzzz_xyyy[i] = 12.0 * ts_zz_xyyy[i] * gfe2_0 + 8.0 * ts_zzz_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_zzzz_xy[i] * gfe2_0 + 6.0 * ts_zzzz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_xyyy[i] * gfe_0 + ts_zzzz_xyyy[i] * rgc2_0;

        gr_zzzz_xyyz[i] = 12.0 * ts_zz_xyyz[i] * gfe2_0 + 8.0 * ts_zzz_xyy[i] * gfe2_0 + 8.0 * ts_zzz_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xz[i] * gfe2_0 + 4.0 * ts_zzzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xyyz[i] * gfe_0 + ts_zzzz_xyyz[i] * rgc2_0;

        gr_zzzz_xyzz[i] = 12.0 * ts_zz_xyzz[i] * gfe2_0 + 16.0 * ts_zzz_xyz[i] * gfe2_0 + 8.0 * ts_zzz_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_xy[i] * gfe2_0 + 4.0 * ts_zzzz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xyzz[i] * gfe_0 + ts_zzzz_xyzz[i] * rgc2_0;

        gr_zzzz_xzzz[i] = 12.0 * ts_zz_xzzz[i] * gfe2_0 + 24.0 * ts_zzz_xzz[i] * gfe2_0 + 8.0 * ts_zzz_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_zzzz_xz[i] * gfe2_0 + 6.0 * ts_zzzz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xzzz[i] * gfe_0 + ts_zzzz_xzzz[i] * rgc2_0;

        gr_zzzz_yyyy[i] = 12.0 * ts_zz_yyyy[i] * gfe2_0 + 8.0 * ts_zzz_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_zzzz_yy[i] * gfe2_0 + 8.0 * ts_zzzz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_yyyy[i] * gfe_0 + ts_zzzz_yyyy[i] * rgc2_0;

        gr_zzzz_yyyz[i] = 12.0 * ts_zz_yyyz[i] * gfe2_0 + 8.0 * ts_zzz_yyy[i] * gfe2_0 + 8.0 * ts_zzz_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzzz_yz[i] * gfe2_0 + 6.0 * ts_zzzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_yyyz[i] * gfe_0 + ts_zzzz_yyyz[i] * rgc2_0;

        gr_zzzz_yyzz[i] = 12.0 * ts_zz_yyzz[i] * gfe2_0 + 16.0 * ts_zzz_yyz[i] * gfe2_0 + 8.0 * ts_zzz_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_zz[i] * gfe2_0 + 4.0 * ts_zzzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_yy[i] * gfe2_0 + 4.0 * ts_zzzz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_yyzz[i] * gfe_0 + ts_zzzz_yyzz[i] * rgc2_0;

        gr_zzzz_yzzz[i] = 12.0 * ts_zz_yzzz[i] * gfe2_0 + 24.0 * ts_zzz_yzz[i] * gfe2_0 + 8.0 * ts_zzz_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_zzzz_yz[i] * gfe2_0 + 6.0 * ts_zzzz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_yzzz[i] * gfe_0 + ts_zzzz_yzzz[i] * rgc2_0;

        gr_zzzz_zzzz[i] = 12.0 * ts_zz_zzzz[i] * gfe2_0 + 32.0 * ts_zzz_zzz[i] * gfe2_0 + 8.0 * ts_zzz_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_zzzz_zz[i] * gfe2_0 + 8.0 * ts_zzzz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_zzzz[i] * gfe_0 + ts_zzzz_zzzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

