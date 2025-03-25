#include "ThreeCenterOverlapGradientPrimRecGF.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_gf(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_gf,
                              const size_t idx_ff,
                              const size_t idx_gd,
                              const size_t idx_gf,
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

    // Set up 0-10 components of targeted buffer : GF

    auto gs_x_xxxx_xxx = pbuffer.data(idx_g_gf);

    auto gs_x_xxxx_xxy = pbuffer.data(idx_g_gf + 1);

    auto gs_x_xxxx_xxz = pbuffer.data(idx_g_gf + 2);

    auto gs_x_xxxx_xyy = pbuffer.data(idx_g_gf + 3);

    auto gs_x_xxxx_xyz = pbuffer.data(idx_g_gf + 4);

    auto gs_x_xxxx_xzz = pbuffer.data(idx_g_gf + 5);

    auto gs_x_xxxx_yyy = pbuffer.data(idx_g_gf + 6);

    auto gs_x_xxxx_yyz = pbuffer.data(idx_g_gf + 7);

    auto gs_x_xxxx_yzz = pbuffer.data(idx_g_gf + 8);

    auto gs_x_xxxx_zzz = pbuffer.data(idx_g_gf + 9);

    #pragma omp simd aligned(gc_x, gs_x_xxxx_xxx, gs_x_xxxx_xxy, gs_x_xxxx_xxz, gs_x_xxxx_xyy, gs_x_xxxx_xyz, gs_x_xxxx_xzz, gs_x_xxxx_yyy, gs_x_xxxx_yyz, gs_x_xxxx_yzz, gs_x_xxxx_zzz, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xzz, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yzz, ts_xxx_zzz, ts_xxxx_xx, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xy, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xz, ts_xxxx_xzz, ts_xxxx_yy, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yz, ts_xxxx_yzz, ts_xxxx_zz, ts_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxx_xxx[i] = 8.0 * ts_xxx_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxy[i] = 8.0 * ts_xxx_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxz[i] = 8.0 * ts_xxx_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyy[i] = 8.0 * ts_xxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyz[i] = 8.0 * ts_xxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xzz[i] = 8.0 * ts_xxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyy[i] = 8.0 * ts_xxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyz[i] = 8.0 * ts_xxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yzz[i] = 8.0 * ts_xxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_zzz[i] = 8.0 * ts_xxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 10-20 components of targeted buffer : GF

    auto gs_x_xxxy_xxx = pbuffer.data(idx_g_gf + 10);

    auto gs_x_xxxy_xxy = pbuffer.data(idx_g_gf + 11);

    auto gs_x_xxxy_xxz = pbuffer.data(idx_g_gf + 12);

    auto gs_x_xxxy_xyy = pbuffer.data(idx_g_gf + 13);

    auto gs_x_xxxy_xyz = pbuffer.data(idx_g_gf + 14);

    auto gs_x_xxxy_xzz = pbuffer.data(idx_g_gf + 15);

    auto gs_x_xxxy_yyy = pbuffer.data(idx_g_gf + 16);

    auto gs_x_xxxy_yyz = pbuffer.data(idx_g_gf + 17);

    auto gs_x_xxxy_yzz = pbuffer.data(idx_g_gf + 18);

    auto gs_x_xxxy_zzz = pbuffer.data(idx_g_gf + 19);

    #pragma omp simd aligned(gc_x, gs_x_xxxy_xxx, gs_x_xxxy_xxy, gs_x_xxxy_xxz, gs_x_xxxy_xyy, gs_x_xxxy_xyz, gs_x_xxxy_xzz, gs_x_xxxy_yyy, gs_x_xxxy_yyz, gs_x_xxxy_yzz, gs_x_xxxy_zzz, ts_xxxy_xx, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xy, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xz, ts_xxxy_xzz, ts_xxxy_yy, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yz, ts_xxxy_yzz, ts_xxxy_zz, ts_xxxy_zzz, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xzz, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yzz, ts_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxy_xxx[i] = 6.0 * ts_xxy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxy[i] = 6.0 * ts_xxy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxz[i] = 6.0 * ts_xxy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyy[i] = 6.0 * ts_xxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyz[i] = 6.0 * ts_xxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xzz[i] = 6.0 * ts_xxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyy[i] = 6.0 * ts_xxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyz[i] = 6.0 * ts_xxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yzz[i] = 6.0 * ts_xxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_zzz[i] = 6.0 * ts_xxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 20-30 components of targeted buffer : GF

    auto gs_x_xxxz_xxx = pbuffer.data(idx_g_gf + 20);

    auto gs_x_xxxz_xxy = pbuffer.data(idx_g_gf + 21);

    auto gs_x_xxxz_xxz = pbuffer.data(idx_g_gf + 22);

    auto gs_x_xxxz_xyy = pbuffer.data(idx_g_gf + 23);

    auto gs_x_xxxz_xyz = pbuffer.data(idx_g_gf + 24);

    auto gs_x_xxxz_xzz = pbuffer.data(idx_g_gf + 25);

    auto gs_x_xxxz_yyy = pbuffer.data(idx_g_gf + 26);

    auto gs_x_xxxz_yyz = pbuffer.data(idx_g_gf + 27);

    auto gs_x_xxxz_yzz = pbuffer.data(idx_g_gf + 28);

    auto gs_x_xxxz_zzz = pbuffer.data(idx_g_gf + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxxz_xxx, gs_x_xxxz_xxy, gs_x_xxxz_xxz, gs_x_xxxz_xyy, gs_x_xxxz_xyz, gs_x_xxxz_xzz, gs_x_xxxz_yyy, gs_x_xxxz_yyz, gs_x_xxxz_yzz, gs_x_xxxz_zzz, ts_xxxz_xx, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xy, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xz, ts_xxxz_xzz, ts_xxxz_yy, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yz, ts_xxxz_yzz, ts_xxxz_zz, ts_xxxz_zzz, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xzz, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yzz, ts_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxz_xxx[i] = 6.0 * ts_xxz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxy[i] = 6.0 * ts_xxz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxz[i] = 6.0 * ts_xxz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyy[i] = 6.0 * ts_xxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyz[i] = 6.0 * ts_xxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xzz[i] = 6.0 * ts_xxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyy[i] = 6.0 * ts_xxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyz[i] = 6.0 * ts_xxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yzz[i] = 6.0 * ts_xxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_zzz[i] = 6.0 * ts_xxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-40 components of targeted buffer : GF

    auto gs_x_xxyy_xxx = pbuffer.data(idx_g_gf + 30);

    auto gs_x_xxyy_xxy = pbuffer.data(idx_g_gf + 31);

    auto gs_x_xxyy_xxz = pbuffer.data(idx_g_gf + 32);

    auto gs_x_xxyy_xyy = pbuffer.data(idx_g_gf + 33);

    auto gs_x_xxyy_xyz = pbuffer.data(idx_g_gf + 34);

    auto gs_x_xxyy_xzz = pbuffer.data(idx_g_gf + 35);

    auto gs_x_xxyy_yyy = pbuffer.data(idx_g_gf + 36);

    auto gs_x_xxyy_yyz = pbuffer.data(idx_g_gf + 37);

    auto gs_x_xxyy_yzz = pbuffer.data(idx_g_gf + 38);

    auto gs_x_xxyy_zzz = pbuffer.data(idx_g_gf + 39);

    #pragma omp simd aligned(gc_x, gs_x_xxyy_xxx, gs_x_xxyy_xxy, gs_x_xxyy_xxz, gs_x_xxyy_xyy, gs_x_xxyy_xyz, gs_x_xxyy_xzz, gs_x_xxyy_yyy, gs_x_xxyy_yyz, gs_x_xxyy_yzz, gs_x_xxyy_zzz, ts_xxyy_xx, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xy, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xz, ts_xxyy_xzz, ts_xxyy_yy, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yz, ts_xxyy_yzz, ts_xxyy_zz, ts_xxyy_zzz, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xzz, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yzz, ts_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyy_xxx[i] = 4.0 * ts_xyy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxy[i] = 4.0 * ts_xyy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxz[i] = 4.0 * ts_xyy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyy[i] = 4.0 * ts_xyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyz[i] = 4.0 * ts_xyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xzz[i] = 4.0 * ts_xyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyy[i] = 4.0 * ts_xyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyz[i] = 4.0 * ts_xyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yzz[i] = 4.0 * ts_xyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_zzz[i] = 4.0 * ts_xyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 40-50 components of targeted buffer : GF

    auto gs_x_xxyz_xxx = pbuffer.data(idx_g_gf + 40);

    auto gs_x_xxyz_xxy = pbuffer.data(idx_g_gf + 41);

    auto gs_x_xxyz_xxz = pbuffer.data(idx_g_gf + 42);

    auto gs_x_xxyz_xyy = pbuffer.data(idx_g_gf + 43);

    auto gs_x_xxyz_xyz = pbuffer.data(idx_g_gf + 44);

    auto gs_x_xxyz_xzz = pbuffer.data(idx_g_gf + 45);

    auto gs_x_xxyz_yyy = pbuffer.data(idx_g_gf + 46);

    auto gs_x_xxyz_yyz = pbuffer.data(idx_g_gf + 47);

    auto gs_x_xxyz_yzz = pbuffer.data(idx_g_gf + 48);

    auto gs_x_xxyz_zzz = pbuffer.data(idx_g_gf + 49);

    #pragma omp simd aligned(gc_x, gs_x_xxyz_xxx, gs_x_xxyz_xxy, gs_x_xxyz_xxz, gs_x_xxyz_xyy, gs_x_xxyz_xyz, gs_x_xxyz_xzz, gs_x_xxyz_yyy, gs_x_xxyz_yyz, gs_x_xxyz_yzz, gs_x_xxyz_zzz, ts_xxyz_xx, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xy, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xz, ts_xxyz_xzz, ts_xxyz_yy, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yz, ts_xxyz_yzz, ts_xxyz_zz, ts_xxyz_zzz, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xzz, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yzz, ts_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyz_xxx[i] = 4.0 * ts_xyz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxy[i] = 4.0 * ts_xyz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxz[i] = 4.0 * ts_xyz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyy[i] = 4.0 * ts_xyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyz[i] = 4.0 * ts_xyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xzz[i] = 4.0 * ts_xyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyy[i] = 4.0 * ts_xyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyz[i] = 4.0 * ts_xyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yzz[i] = 4.0 * ts_xyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_zzz[i] = 4.0 * ts_xyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 50-60 components of targeted buffer : GF

    auto gs_x_xxzz_xxx = pbuffer.data(idx_g_gf + 50);

    auto gs_x_xxzz_xxy = pbuffer.data(idx_g_gf + 51);

    auto gs_x_xxzz_xxz = pbuffer.data(idx_g_gf + 52);

    auto gs_x_xxzz_xyy = pbuffer.data(idx_g_gf + 53);

    auto gs_x_xxzz_xyz = pbuffer.data(idx_g_gf + 54);

    auto gs_x_xxzz_xzz = pbuffer.data(idx_g_gf + 55);

    auto gs_x_xxzz_yyy = pbuffer.data(idx_g_gf + 56);

    auto gs_x_xxzz_yyz = pbuffer.data(idx_g_gf + 57);

    auto gs_x_xxzz_yzz = pbuffer.data(idx_g_gf + 58);

    auto gs_x_xxzz_zzz = pbuffer.data(idx_g_gf + 59);

    #pragma omp simd aligned(gc_x, gs_x_xxzz_xxx, gs_x_xxzz_xxy, gs_x_xxzz_xxz, gs_x_xxzz_xyy, gs_x_xxzz_xyz, gs_x_xxzz_xzz, gs_x_xxzz_yyy, gs_x_xxzz_yyz, gs_x_xxzz_yzz, gs_x_xxzz_zzz, ts_xxzz_xx, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xy, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xz, ts_xxzz_xzz, ts_xxzz_yy, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yz, ts_xxzz_yzz, ts_xxzz_zz, ts_xxzz_zzz, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xzz, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yzz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzz_xxx[i] = 4.0 * ts_xzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxy[i] = 4.0 * ts_xzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxz[i] = 4.0 * ts_xzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyy[i] = 4.0 * ts_xzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyz[i] = 4.0 * ts_xzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xzz[i] = 4.0 * ts_xzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyy[i] = 4.0 * ts_xzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyz[i] = 4.0 * ts_xzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yzz[i] = 4.0 * ts_xzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_zzz[i] = 4.0 * ts_xzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-70 components of targeted buffer : GF

    auto gs_x_xyyy_xxx = pbuffer.data(idx_g_gf + 60);

    auto gs_x_xyyy_xxy = pbuffer.data(idx_g_gf + 61);

    auto gs_x_xyyy_xxz = pbuffer.data(idx_g_gf + 62);

    auto gs_x_xyyy_xyy = pbuffer.data(idx_g_gf + 63);

    auto gs_x_xyyy_xyz = pbuffer.data(idx_g_gf + 64);

    auto gs_x_xyyy_xzz = pbuffer.data(idx_g_gf + 65);

    auto gs_x_xyyy_yyy = pbuffer.data(idx_g_gf + 66);

    auto gs_x_xyyy_yyz = pbuffer.data(idx_g_gf + 67);

    auto gs_x_xyyy_yzz = pbuffer.data(idx_g_gf + 68);

    auto gs_x_xyyy_zzz = pbuffer.data(idx_g_gf + 69);

    #pragma omp simd aligned(gc_x, gs_x_xyyy_xxx, gs_x_xyyy_xxy, gs_x_xyyy_xxz, gs_x_xyyy_xyy, gs_x_xyyy_xyz, gs_x_xyyy_xzz, gs_x_xyyy_yyy, gs_x_xyyy_yyz, gs_x_xyyy_yzz, gs_x_xyyy_zzz, ts_xyyy_xx, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xy, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xz, ts_xyyy_xzz, ts_xyyy_yy, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yz, ts_xyyy_yzz, ts_xyyy_zz, ts_xyyy_zzz, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xzz, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yzz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyy_xxx[i] = 2.0 * ts_yyy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxy[i] = 2.0 * ts_yyy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxz[i] = 2.0 * ts_yyy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyy[i] = 2.0 * ts_yyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyz[i] = 2.0 * ts_yyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xzz[i] = 2.0 * ts_yyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyy[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyz[i] = 2.0 * ts_yyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yzz[i] = 2.0 * ts_yyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_zzz[i] = 2.0 * ts_yyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 70-80 components of targeted buffer : GF

    auto gs_x_xyyz_xxx = pbuffer.data(idx_g_gf + 70);

    auto gs_x_xyyz_xxy = pbuffer.data(idx_g_gf + 71);

    auto gs_x_xyyz_xxz = pbuffer.data(idx_g_gf + 72);

    auto gs_x_xyyz_xyy = pbuffer.data(idx_g_gf + 73);

    auto gs_x_xyyz_xyz = pbuffer.data(idx_g_gf + 74);

    auto gs_x_xyyz_xzz = pbuffer.data(idx_g_gf + 75);

    auto gs_x_xyyz_yyy = pbuffer.data(idx_g_gf + 76);

    auto gs_x_xyyz_yyz = pbuffer.data(idx_g_gf + 77);

    auto gs_x_xyyz_yzz = pbuffer.data(idx_g_gf + 78);

    auto gs_x_xyyz_zzz = pbuffer.data(idx_g_gf + 79);

    #pragma omp simd aligned(gc_x, gs_x_xyyz_xxx, gs_x_xyyz_xxy, gs_x_xyyz_xxz, gs_x_xyyz_xyy, gs_x_xyyz_xyz, gs_x_xyyz_xzz, gs_x_xyyz_yyy, gs_x_xyyz_yyz, gs_x_xyyz_yzz, gs_x_xyyz_zzz, ts_xyyz_xx, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xy, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xz, ts_xyyz_xzz, ts_xyyz_yy, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yz, ts_xyyz_yzz, ts_xyyz_zz, ts_xyyz_zzz, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xzz, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yzz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyz_xxx[i] = 2.0 * ts_yyz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxy[i] = 2.0 * ts_yyz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxz[i] = 2.0 * ts_yyz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyy[i] = 2.0 * ts_yyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyz[i] = 2.0 * ts_yyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xzz[i] = 2.0 * ts_yyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyy[i] = 2.0 * ts_yyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyz[i] = 2.0 * ts_yyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yzz[i] = 2.0 * ts_yyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_zzz[i] = 2.0 * ts_yyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 80-90 components of targeted buffer : GF

    auto gs_x_xyzz_xxx = pbuffer.data(idx_g_gf + 80);

    auto gs_x_xyzz_xxy = pbuffer.data(idx_g_gf + 81);

    auto gs_x_xyzz_xxz = pbuffer.data(idx_g_gf + 82);

    auto gs_x_xyzz_xyy = pbuffer.data(idx_g_gf + 83);

    auto gs_x_xyzz_xyz = pbuffer.data(idx_g_gf + 84);

    auto gs_x_xyzz_xzz = pbuffer.data(idx_g_gf + 85);

    auto gs_x_xyzz_yyy = pbuffer.data(idx_g_gf + 86);

    auto gs_x_xyzz_yyz = pbuffer.data(idx_g_gf + 87);

    auto gs_x_xyzz_yzz = pbuffer.data(idx_g_gf + 88);

    auto gs_x_xyzz_zzz = pbuffer.data(idx_g_gf + 89);

    #pragma omp simd aligned(gc_x, gs_x_xyzz_xxx, gs_x_xyzz_xxy, gs_x_xyzz_xxz, gs_x_xyzz_xyy, gs_x_xyzz_xyz, gs_x_xyzz_xzz, gs_x_xyzz_yyy, gs_x_xyzz_yyz, gs_x_xyzz_yzz, gs_x_xyzz_zzz, ts_xyzz_xx, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xy, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xz, ts_xyzz_xzz, ts_xyzz_yy, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yz, ts_xyzz_yzz, ts_xyzz_zz, ts_xyzz_zzz, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xzz, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yzz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzz_xxx[i] = 2.0 * ts_yzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxy[i] = 2.0 * ts_yzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxz[i] = 2.0 * ts_yzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyy[i] = 2.0 * ts_yzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyz[i] = 2.0 * ts_yzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xzz[i] = 2.0 * ts_yzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyy[i] = 2.0 * ts_yzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyz[i] = 2.0 * ts_yzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yzz[i] = 2.0 * ts_yzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_zzz[i] = 2.0 * ts_yzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-100 components of targeted buffer : GF

    auto gs_x_xzzz_xxx = pbuffer.data(idx_g_gf + 90);

    auto gs_x_xzzz_xxy = pbuffer.data(idx_g_gf + 91);

    auto gs_x_xzzz_xxz = pbuffer.data(idx_g_gf + 92);

    auto gs_x_xzzz_xyy = pbuffer.data(idx_g_gf + 93);

    auto gs_x_xzzz_xyz = pbuffer.data(idx_g_gf + 94);

    auto gs_x_xzzz_xzz = pbuffer.data(idx_g_gf + 95);

    auto gs_x_xzzz_yyy = pbuffer.data(idx_g_gf + 96);

    auto gs_x_xzzz_yyz = pbuffer.data(idx_g_gf + 97);

    auto gs_x_xzzz_yzz = pbuffer.data(idx_g_gf + 98);

    auto gs_x_xzzz_zzz = pbuffer.data(idx_g_gf + 99);

    #pragma omp simd aligned(gc_x, gs_x_xzzz_xxx, gs_x_xzzz_xxy, gs_x_xzzz_xxz, gs_x_xzzz_xyy, gs_x_xzzz_xyz, gs_x_xzzz_xzz, gs_x_xzzz_yyy, gs_x_xzzz_yyz, gs_x_xzzz_yzz, gs_x_xzzz_zzz, ts_xzzz_xx, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xy, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xz, ts_xzzz_xzz, ts_xzzz_yy, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yz, ts_xzzz_yzz, ts_xzzz_zz, ts_xzzz_zzz, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xzz, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yzz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzz_xxx[i] = 2.0 * ts_zzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxy[i] = 2.0 * ts_zzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxz[i] = 2.0 * ts_zzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyy[i] = 2.0 * ts_zzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyz[i] = 2.0 * ts_zzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xzz[i] = 2.0 * ts_zzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyy[i] = 2.0 * ts_zzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyz[i] = 2.0 * ts_zzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yzz[i] = 2.0 * ts_zzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 100-110 components of targeted buffer : GF

    auto gs_x_yyyy_xxx = pbuffer.data(idx_g_gf + 100);

    auto gs_x_yyyy_xxy = pbuffer.data(idx_g_gf + 101);

    auto gs_x_yyyy_xxz = pbuffer.data(idx_g_gf + 102);

    auto gs_x_yyyy_xyy = pbuffer.data(idx_g_gf + 103);

    auto gs_x_yyyy_xyz = pbuffer.data(idx_g_gf + 104);

    auto gs_x_yyyy_xzz = pbuffer.data(idx_g_gf + 105);

    auto gs_x_yyyy_yyy = pbuffer.data(idx_g_gf + 106);

    auto gs_x_yyyy_yyz = pbuffer.data(idx_g_gf + 107);

    auto gs_x_yyyy_yzz = pbuffer.data(idx_g_gf + 108);

    auto gs_x_yyyy_zzz = pbuffer.data(idx_g_gf + 109);

    #pragma omp simd aligned(gc_x, gs_x_yyyy_xxx, gs_x_yyyy_xxy, gs_x_yyyy_xxz, gs_x_yyyy_xyy, gs_x_yyyy_xyz, gs_x_yyyy_xzz, gs_x_yyyy_yyy, gs_x_yyyy_yyz, gs_x_yyyy_yzz, gs_x_yyyy_zzz, ts_yyyy_xx, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xy, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xz, ts_yyyy_xzz, ts_yyyy_yy, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yz, ts_yyyy_yzz, ts_yyyy_zz, ts_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyy_xxx[i] = 6.0 * ts_yyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxy[i] = 4.0 * ts_yyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxz[i] = 4.0 * ts_yyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyy[i] = 2.0 * ts_yyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyz[i] = 2.0 * ts_yyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xzz[i] = 2.0 * ts_yyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyy[i] = 2.0 * ts_yyyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyz[i] = 2.0 * ts_yyyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yzz[i] = 2.0 * ts_yyyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_zzz[i] = 2.0 * ts_yyyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 110-120 components of targeted buffer : GF

    auto gs_x_yyyz_xxx = pbuffer.data(idx_g_gf + 110);

    auto gs_x_yyyz_xxy = pbuffer.data(idx_g_gf + 111);

    auto gs_x_yyyz_xxz = pbuffer.data(idx_g_gf + 112);

    auto gs_x_yyyz_xyy = pbuffer.data(idx_g_gf + 113);

    auto gs_x_yyyz_xyz = pbuffer.data(idx_g_gf + 114);

    auto gs_x_yyyz_xzz = pbuffer.data(idx_g_gf + 115);

    auto gs_x_yyyz_yyy = pbuffer.data(idx_g_gf + 116);

    auto gs_x_yyyz_yyz = pbuffer.data(idx_g_gf + 117);

    auto gs_x_yyyz_yzz = pbuffer.data(idx_g_gf + 118);

    auto gs_x_yyyz_zzz = pbuffer.data(idx_g_gf + 119);

    #pragma omp simd aligned(gc_x, gs_x_yyyz_xxx, gs_x_yyyz_xxy, gs_x_yyyz_xxz, gs_x_yyyz_xyy, gs_x_yyyz_xyz, gs_x_yyyz_xzz, gs_x_yyyz_yyy, gs_x_yyyz_yyz, gs_x_yyyz_yzz, gs_x_yyyz_zzz, ts_yyyz_xx, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xy, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xz, ts_yyyz_xzz, ts_yyyz_yy, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yz, ts_yyyz_yzz, ts_yyyz_zz, ts_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyz_xxx[i] = 6.0 * ts_yyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxy[i] = 4.0 * ts_yyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxz[i] = 4.0 * ts_yyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyy[i] = 2.0 * ts_yyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyz[i] = 2.0 * ts_yyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xzz[i] = 2.0 * ts_yyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyy[i] = 2.0 * ts_yyyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyz[i] = 2.0 * ts_yyyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yzz[i] = 2.0 * ts_yyyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_zzz[i] = 2.0 * ts_yyyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 120-130 components of targeted buffer : GF

    auto gs_x_yyzz_xxx = pbuffer.data(idx_g_gf + 120);

    auto gs_x_yyzz_xxy = pbuffer.data(idx_g_gf + 121);

    auto gs_x_yyzz_xxz = pbuffer.data(idx_g_gf + 122);

    auto gs_x_yyzz_xyy = pbuffer.data(idx_g_gf + 123);

    auto gs_x_yyzz_xyz = pbuffer.data(idx_g_gf + 124);

    auto gs_x_yyzz_xzz = pbuffer.data(idx_g_gf + 125);

    auto gs_x_yyzz_yyy = pbuffer.data(idx_g_gf + 126);

    auto gs_x_yyzz_yyz = pbuffer.data(idx_g_gf + 127);

    auto gs_x_yyzz_yzz = pbuffer.data(idx_g_gf + 128);

    auto gs_x_yyzz_zzz = pbuffer.data(idx_g_gf + 129);

    #pragma omp simd aligned(gc_x, gs_x_yyzz_xxx, gs_x_yyzz_xxy, gs_x_yyzz_xxz, gs_x_yyzz_xyy, gs_x_yyzz_xyz, gs_x_yyzz_xzz, gs_x_yyzz_yyy, gs_x_yyzz_yyz, gs_x_yyzz_yzz, gs_x_yyzz_zzz, ts_yyzz_xx, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xy, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xz, ts_yyzz_xzz, ts_yyzz_yy, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yz, ts_yyzz_yzz, ts_yyzz_zz, ts_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzz_xxx[i] = 6.0 * ts_yyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxy[i] = 4.0 * ts_yyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxz[i] = 4.0 * ts_yyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyy[i] = 2.0 * ts_yyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyz[i] = 2.0 * ts_yyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xzz[i] = 2.0 * ts_yyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyy[i] = 2.0 * ts_yyzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyz[i] = 2.0 * ts_yyzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yzz[i] = 2.0 * ts_yyzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_zzz[i] = 2.0 * ts_yyzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 130-140 components of targeted buffer : GF

    auto gs_x_yzzz_xxx = pbuffer.data(idx_g_gf + 130);

    auto gs_x_yzzz_xxy = pbuffer.data(idx_g_gf + 131);

    auto gs_x_yzzz_xxz = pbuffer.data(idx_g_gf + 132);

    auto gs_x_yzzz_xyy = pbuffer.data(idx_g_gf + 133);

    auto gs_x_yzzz_xyz = pbuffer.data(idx_g_gf + 134);

    auto gs_x_yzzz_xzz = pbuffer.data(idx_g_gf + 135);

    auto gs_x_yzzz_yyy = pbuffer.data(idx_g_gf + 136);

    auto gs_x_yzzz_yyz = pbuffer.data(idx_g_gf + 137);

    auto gs_x_yzzz_yzz = pbuffer.data(idx_g_gf + 138);

    auto gs_x_yzzz_zzz = pbuffer.data(idx_g_gf + 139);

    #pragma omp simd aligned(gc_x, gs_x_yzzz_xxx, gs_x_yzzz_xxy, gs_x_yzzz_xxz, gs_x_yzzz_xyy, gs_x_yzzz_xyz, gs_x_yzzz_xzz, gs_x_yzzz_yyy, gs_x_yzzz_yyz, gs_x_yzzz_yzz, gs_x_yzzz_zzz, ts_yzzz_xx, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xy, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xz, ts_yzzz_xzz, ts_yzzz_yy, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yz, ts_yzzz_yzz, ts_yzzz_zz, ts_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzz_xxx[i] = 6.0 * ts_yzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxy[i] = 4.0 * ts_yzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxz[i] = 4.0 * ts_yzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyy[i] = 2.0 * ts_yzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyz[i] = 2.0 * ts_yzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xzz[i] = 2.0 * ts_yzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyy[i] = 2.0 * ts_yzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyz[i] = 2.0 * ts_yzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yzz[i] = 2.0 * ts_yzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_zzz[i] = 2.0 * ts_yzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 140-150 components of targeted buffer : GF

    auto gs_x_zzzz_xxx = pbuffer.data(idx_g_gf + 140);

    auto gs_x_zzzz_xxy = pbuffer.data(idx_g_gf + 141);

    auto gs_x_zzzz_xxz = pbuffer.data(idx_g_gf + 142);

    auto gs_x_zzzz_xyy = pbuffer.data(idx_g_gf + 143);

    auto gs_x_zzzz_xyz = pbuffer.data(idx_g_gf + 144);

    auto gs_x_zzzz_xzz = pbuffer.data(idx_g_gf + 145);

    auto gs_x_zzzz_yyy = pbuffer.data(idx_g_gf + 146);

    auto gs_x_zzzz_yyz = pbuffer.data(idx_g_gf + 147);

    auto gs_x_zzzz_yzz = pbuffer.data(idx_g_gf + 148);

    auto gs_x_zzzz_zzz = pbuffer.data(idx_g_gf + 149);

    #pragma omp simd aligned(gc_x, gs_x_zzzz_xxx, gs_x_zzzz_xxy, gs_x_zzzz_xxz, gs_x_zzzz_xyy, gs_x_zzzz_xyz, gs_x_zzzz_xzz, gs_x_zzzz_yyy, gs_x_zzzz_yyz, gs_x_zzzz_yzz, gs_x_zzzz_zzz, ts_zzzz_xx, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xy, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xz, ts_zzzz_xzz, ts_zzzz_yy, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yz, ts_zzzz_yzz, ts_zzzz_zz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzz_xxx[i] = 6.0 * ts_zzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxy[i] = 4.0 * ts_zzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxz[i] = 4.0 * ts_zzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyy[i] = 2.0 * ts_zzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyz[i] = 2.0 * ts_zzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xzz[i] = 2.0 * ts_zzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyy[i] = 2.0 * ts_zzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyz[i] = 2.0 * ts_zzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yzz[i] = 2.0 * ts_zzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_zzz[i] = 2.0 * ts_zzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 150-160 components of targeted buffer : GF

    auto gs_y_xxxx_xxx = pbuffer.data(idx_g_gf + 150);

    auto gs_y_xxxx_xxy = pbuffer.data(idx_g_gf + 151);

    auto gs_y_xxxx_xxz = pbuffer.data(idx_g_gf + 152);

    auto gs_y_xxxx_xyy = pbuffer.data(idx_g_gf + 153);

    auto gs_y_xxxx_xyz = pbuffer.data(idx_g_gf + 154);

    auto gs_y_xxxx_xzz = pbuffer.data(idx_g_gf + 155);

    auto gs_y_xxxx_yyy = pbuffer.data(idx_g_gf + 156);

    auto gs_y_xxxx_yyz = pbuffer.data(idx_g_gf + 157);

    auto gs_y_xxxx_yzz = pbuffer.data(idx_g_gf + 158);

    auto gs_y_xxxx_zzz = pbuffer.data(idx_g_gf + 159);

    #pragma omp simd aligned(gc_y, gs_y_xxxx_xxx, gs_y_xxxx_xxy, gs_y_xxxx_xxz, gs_y_xxxx_xyy, gs_y_xxxx_xyz, gs_y_xxxx_xzz, gs_y_xxxx_yyy, gs_y_xxxx_yyz, gs_y_xxxx_yzz, gs_y_xxxx_zzz, ts_xxxx_xx, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xy, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xz, ts_xxxx_xzz, ts_xxxx_yy, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yz, ts_xxxx_yzz, ts_xxxx_zz, ts_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxx_xxx[i] = 2.0 * ts_xxxx_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxy[i] = 2.0 * ts_xxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxz[i] = 2.0 * ts_xxxx_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyy[i] = 4.0 * ts_xxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyz[i] = 2.0 * ts_xxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xzz[i] = 2.0 * ts_xxxx_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyy[i] = 6.0 * ts_xxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyz[i] = 4.0 * ts_xxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yzz[i] = 2.0 * ts_xxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_zzz[i] = 2.0 * ts_xxxx_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 160-170 components of targeted buffer : GF

    auto gs_y_xxxy_xxx = pbuffer.data(idx_g_gf + 160);

    auto gs_y_xxxy_xxy = pbuffer.data(idx_g_gf + 161);

    auto gs_y_xxxy_xxz = pbuffer.data(idx_g_gf + 162);

    auto gs_y_xxxy_xyy = pbuffer.data(idx_g_gf + 163);

    auto gs_y_xxxy_xyz = pbuffer.data(idx_g_gf + 164);

    auto gs_y_xxxy_xzz = pbuffer.data(idx_g_gf + 165);

    auto gs_y_xxxy_yyy = pbuffer.data(idx_g_gf + 166);

    auto gs_y_xxxy_yyz = pbuffer.data(idx_g_gf + 167);

    auto gs_y_xxxy_yzz = pbuffer.data(idx_g_gf + 168);

    auto gs_y_xxxy_zzz = pbuffer.data(idx_g_gf + 169);

    #pragma omp simd aligned(gc_y, gs_y_xxxy_xxx, gs_y_xxxy_xxy, gs_y_xxxy_xxz, gs_y_xxxy_xyy, gs_y_xxxy_xyz, gs_y_xxxy_xzz, gs_y_xxxy_yyy, gs_y_xxxy_yyz, gs_y_xxxy_yzz, gs_y_xxxy_zzz, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xzz, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yzz, ts_xxx_zzz, ts_xxxy_xx, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xy, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xz, ts_xxxy_xzz, ts_xxxy_yy, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yz, ts_xxxy_yzz, ts_xxxy_zz, ts_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxy_xxx[i] = 2.0 * ts_xxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxy[i] = 2.0 * ts_xxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxz[i] = 2.0 * ts_xxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyy[i] = 2.0 * ts_xxx_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyz[i] = 2.0 * ts_xxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xzz[i] = 2.0 * ts_xxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyy[i] = 2.0 * ts_xxx_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyz[i] = 2.0 * ts_xxx_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yzz[i] = 2.0 * ts_xxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_zzz[i] = 2.0 * ts_xxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 170-180 components of targeted buffer : GF

    auto gs_y_xxxz_xxx = pbuffer.data(idx_g_gf + 170);

    auto gs_y_xxxz_xxy = pbuffer.data(idx_g_gf + 171);

    auto gs_y_xxxz_xxz = pbuffer.data(idx_g_gf + 172);

    auto gs_y_xxxz_xyy = pbuffer.data(idx_g_gf + 173);

    auto gs_y_xxxz_xyz = pbuffer.data(idx_g_gf + 174);

    auto gs_y_xxxz_xzz = pbuffer.data(idx_g_gf + 175);

    auto gs_y_xxxz_yyy = pbuffer.data(idx_g_gf + 176);

    auto gs_y_xxxz_yyz = pbuffer.data(idx_g_gf + 177);

    auto gs_y_xxxz_yzz = pbuffer.data(idx_g_gf + 178);

    auto gs_y_xxxz_zzz = pbuffer.data(idx_g_gf + 179);

    #pragma omp simd aligned(gc_y, gs_y_xxxz_xxx, gs_y_xxxz_xxy, gs_y_xxxz_xxz, gs_y_xxxz_xyy, gs_y_xxxz_xyz, gs_y_xxxz_xzz, gs_y_xxxz_yyy, gs_y_xxxz_yyz, gs_y_xxxz_yzz, gs_y_xxxz_zzz, ts_xxxz_xx, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xy, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xz, ts_xxxz_xzz, ts_xxxz_yy, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yz, ts_xxxz_yzz, ts_xxxz_zz, ts_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxz_xxx[i] = 2.0 * ts_xxxz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxy[i] = 2.0 * ts_xxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxz[i] = 2.0 * ts_xxxz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyy[i] = 4.0 * ts_xxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyz[i] = 2.0 * ts_xxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xzz[i] = 2.0 * ts_xxxz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyy[i] = 6.0 * ts_xxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyz[i] = 4.0 * ts_xxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yzz[i] = 2.0 * ts_xxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_zzz[i] = 2.0 * ts_xxxz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 180-190 components of targeted buffer : GF

    auto gs_y_xxyy_xxx = pbuffer.data(idx_g_gf + 180);

    auto gs_y_xxyy_xxy = pbuffer.data(idx_g_gf + 181);

    auto gs_y_xxyy_xxz = pbuffer.data(idx_g_gf + 182);

    auto gs_y_xxyy_xyy = pbuffer.data(idx_g_gf + 183);

    auto gs_y_xxyy_xyz = pbuffer.data(idx_g_gf + 184);

    auto gs_y_xxyy_xzz = pbuffer.data(idx_g_gf + 185);

    auto gs_y_xxyy_yyy = pbuffer.data(idx_g_gf + 186);

    auto gs_y_xxyy_yyz = pbuffer.data(idx_g_gf + 187);

    auto gs_y_xxyy_yzz = pbuffer.data(idx_g_gf + 188);

    auto gs_y_xxyy_zzz = pbuffer.data(idx_g_gf + 189);

    #pragma omp simd aligned(gc_y, gs_y_xxyy_xxx, gs_y_xxyy_xxy, gs_y_xxyy_xxz, gs_y_xxyy_xyy, gs_y_xxyy_xyz, gs_y_xxyy_xzz, gs_y_xxyy_yyy, gs_y_xxyy_yyz, gs_y_xxyy_yzz, gs_y_xxyy_zzz, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xzz, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yzz, ts_xxy_zzz, ts_xxyy_xx, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xy, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xz, ts_xxyy_xzz, ts_xxyy_yy, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yz, ts_xxyy_yzz, ts_xxyy_zz, ts_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyy_xxx[i] = 4.0 * ts_xxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxy[i] = 4.0 * ts_xxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxz[i] = 4.0 * ts_xxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyy[i] = 4.0 * ts_xxy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyz[i] = 4.0 * ts_xxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xzz[i] = 4.0 * ts_xxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyy[i] = 4.0 * ts_xxy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyz[i] = 4.0 * ts_xxy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yzz[i] = 4.0 * ts_xxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_zzz[i] = 4.0 * ts_xxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 190-200 components of targeted buffer : GF

    auto gs_y_xxyz_xxx = pbuffer.data(idx_g_gf + 190);

    auto gs_y_xxyz_xxy = pbuffer.data(idx_g_gf + 191);

    auto gs_y_xxyz_xxz = pbuffer.data(idx_g_gf + 192);

    auto gs_y_xxyz_xyy = pbuffer.data(idx_g_gf + 193);

    auto gs_y_xxyz_xyz = pbuffer.data(idx_g_gf + 194);

    auto gs_y_xxyz_xzz = pbuffer.data(idx_g_gf + 195);

    auto gs_y_xxyz_yyy = pbuffer.data(idx_g_gf + 196);

    auto gs_y_xxyz_yyz = pbuffer.data(idx_g_gf + 197);

    auto gs_y_xxyz_yzz = pbuffer.data(idx_g_gf + 198);

    auto gs_y_xxyz_zzz = pbuffer.data(idx_g_gf + 199);

    #pragma omp simd aligned(gc_y, gs_y_xxyz_xxx, gs_y_xxyz_xxy, gs_y_xxyz_xxz, gs_y_xxyz_xyy, gs_y_xxyz_xyz, gs_y_xxyz_xzz, gs_y_xxyz_yyy, gs_y_xxyz_yyz, gs_y_xxyz_yzz, gs_y_xxyz_zzz, ts_xxyz_xx, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xy, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xz, ts_xxyz_xzz, ts_xxyz_yy, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yz, ts_xxyz_yzz, ts_xxyz_zz, ts_xxyz_zzz, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xzz, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yzz, ts_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyz_xxx[i] = 2.0 * ts_xxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxy[i] = 2.0 * ts_xxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxz[i] = 2.0 * ts_xxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyy[i] = 2.0 * ts_xxz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyz[i] = 2.0 * ts_xxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xzz[i] = 2.0 * ts_xxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyy[i] = 2.0 * ts_xxz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyz[i] = 2.0 * ts_xxz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yzz[i] = 2.0 * ts_xxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_zzz[i] = 2.0 * ts_xxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 200-210 components of targeted buffer : GF

    auto gs_y_xxzz_xxx = pbuffer.data(idx_g_gf + 200);

    auto gs_y_xxzz_xxy = pbuffer.data(idx_g_gf + 201);

    auto gs_y_xxzz_xxz = pbuffer.data(idx_g_gf + 202);

    auto gs_y_xxzz_xyy = pbuffer.data(idx_g_gf + 203);

    auto gs_y_xxzz_xyz = pbuffer.data(idx_g_gf + 204);

    auto gs_y_xxzz_xzz = pbuffer.data(idx_g_gf + 205);

    auto gs_y_xxzz_yyy = pbuffer.data(idx_g_gf + 206);

    auto gs_y_xxzz_yyz = pbuffer.data(idx_g_gf + 207);

    auto gs_y_xxzz_yzz = pbuffer.data(idx_g_gf + 208);

    auto gs_y_xxzz_zzz = pbuffer.data(idx_g_gf + 209);

    #pragma omp simd aligned(gc_y, gs_y_xxzz_xxx, gs_y_xxzz_xxy, gs_y_xxzz_xxz, gs_y_xxzz_xyy, gs_y_xxzz_xyz, gs_y_xxzz_xzz, gs_y_xxzz_yyy, gs_y_xxzz_yyz, gs_y_xxzz_yzz, gs_y_xxzz_zzz, ts_xxzz_xx, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xy, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xz, ts_xxzz_xzz, ts_xxzz_yy, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yz, ts_xxzz_yzz, ts_xxzz_zz, ts_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzz_xxx[i] = 2.0 * ts_xxzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxy[i] = 2.0 * ts_xxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxz[i] = 2.0 * ts_xxzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyy[i] = 4.0 * ts_xxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyz[i] = 2.0 * ts_xxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xzz[i] = 2.0 * ts_xxzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyy[i] = 6.0 * ts_xxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyz[i] = 4.0 * ts_xxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yzz[i] = 2.0 * ts_xxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_zzz[i] = 2.0 * ts_xxzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 210-220 components of targeted buffer : GF

    auto gs_y_xyyy_xxx = pbuffer.data(idx_g_gf + 210);

    auto gs_y_xyyy_xxy = pbuffer.data(idx_g_gf + 211);

    auto gs_y_xyyy_xxz = pbuffer.data(idx_g_gf + 212);

    auto gs_y_xyyy_xyy = pbuffer.data(idx_g_gf + 213);

    auto gs_y_xyyy_xyz = pbuffer.data(idx_g_gf + 214);

    auto gs_y_xyyy_xzz = pbuffer.data(idx_g_gf + 215);

    auto gs_y_xyyy_yyy = pbuffer.data(idx_g_gf + 216);

    auto gs_y_xyyy_yyz = pbuffer.data(idx_g_gf + 217);

    auto gs_y_xyyy_yzz = pbuffer.data(idx_g_gf + 218);

    auto gs_y_xyyy_zzz = pbuffer.data(idx_g_gf + 219);

    #pragma omp simd aligned(gc_y, gs_y_xyyy_xxx, gs_y_xyyy_xxy, gs_y_xyyy_xxz, gs_y_xyyy_xyy, gs_y_xyyy_xyz, gs_y_xyyy_xzz, gs_y_xyyy_yyy, gs_y_xyyy_yyz, gs_y_xyyy_yzz, gs_y_xyyy_zzz, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xzz, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yzz, ts_xyy_zzz, ts_xyyy_xx, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xy, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xz, ts_xyyy_xzz, ts_xyyy_yy, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yz, ts_xyyy_yzz, ts_xyyy_zz, ts_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyy_xxx[i] = 6.0 * ts_xyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxy[i] = 6.0 * ts_xyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxz[i] = 6.0 * ts_xyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyy[i] = 6.0 * ts_xyy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyz[i] = 6.0 * ts_xyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xzz[i] = 6.0 * ts_xyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyy[i] = 6.0 * ts_xyy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyz[i] = 6.0 * ts_xyy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yzz[i] = 6.0 * ts_xyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_zzz[i] = 6.0 * ts_xyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 220-230 components of targeted buffer : GF

    auto gs_y_xyyz_xxx = pbuffer.data(idx_g_gf + 220);

    auto gs_y_xyyz_xxy = pbuffer.data(idx_g_gf + 221);

    auto gs_y_xyyz_xxz = pbuffer.data(idx_g_gf + 222);

    auto gs_y_xyyz_xyy = pbuffer.data(idx_g_gf + 223);

    auto gs_y_xyyz_xyz = pbuffer.data(idx_g_gf + 224);

    auto gs_y_xyyz_xzz = pbuffer.data(idx_g_gf + 225);

    auto gs_y_xyyz_yyy = pbuffer.data(idx_g_gf + 226);

    auto gs_y_xyyz_yyz = pbuffer.data(idx_g_gf + 227);

    auto gs_y_xyyz_yzz = pbuffer.data(idx_g_gf + 228);

    auto gs_y_xyyz_zzz = pbuffer.data(idx_g_gf + 229);

    #pragma omp simd aligned(gc_y, gs_y_xyyz_xxx, gs_y_xyyz_xxy, gs_y_xyyz_xxz, gs_y_xyyz_xyy, gs_y_xyyz_xyz, gs_y_xyyz_xzz, gs_y_xyyz_yyy, gs_y_xyyz_yyz, gs_y_xyyz_yzz, gs_y_xyyz_zzz, ts_xyyz_xx, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xy, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xz, ts_xyyz_xzz, ts_xyyz_yy, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yz, ts_xyyz_yzz, ts_xyyz_zz, ts_xyyz_zzz, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xzz, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yzz, ts_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyz_xxx[i] = 4.0 * ts_xyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxy[i] = 4.0 * ts_xyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxz[i] = 4.0 * ts_xyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyy[i] = 4.0 * ts_xyz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyz[i] = 4.0 * ts_xyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xzz[i] = 4.0 * ts_xyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyy[i] = 4.0 * ts_xyz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyz[i] = 4.0 * ts_xyz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yzz[i] = 4.0 * ts_xyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_zzz[i] = 4.0 * ts_xyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 230-240 components of targeted buffer : GF

    auto gs_y_xyzz_xxx = pbuffer.data(idx_g_gf + 230);

    auto gs_y_xyzz_xxy = pbuffer.data(idx_g_gf + 231);

    auto gs_y_xyzz_xxz = pbuffer.data(idx_g_gf + 232);

    auto gs_y_xyzz_xyy = pbuffer.data(idx_g_gf + 233);

    auto gs_y_xyzz_xyz = pbuffer.data(idx_g_gf + 234);

    auto gs_y_xyzz_xzz = pbuffer.data(idx_g_gf + 235);

    auto gs_y_xyzz_yyy = pbuffer.data(idx_g_gf + 236);

    auto gs_y_xyzz_yyz = pbuffer.data(idx_g_gf + 237);

    auto gs_y_xyzz_yzz = pbuffer.data(idx_g_gf + 238);

    auto gs_y_xyzz_zzz = pbuffer.data(idx_g_gf + 239);

    #pragma omp simd aligned(gc_y, gs_y_xyzz_xxx, gs_y_xyzz_xxy, gs_y_xyzz_xxz, gs_y_xyzz_xyy, gs_y_xyzz_xyz, gs_y_xyzz_xzz, gs_y_xyzz_yyy, gs_y_xyzz_yyz, gs_y_xyzz_yzz, gs_y_xyzz_zzz, ts_xyzz_xx, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xy, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xz, ts_xyzz_xzz, ts_xyzz_yy, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yz, ts_xyzz_yzz, ts_xyzz_zz, ts_xyzz_zzz, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xzz, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yzz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzz_xxx[i] = 2.0 * ts_xzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxy[i] = 2.0 * ts_xzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxz[i] = 2.0 * ts_xzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyy[i] = 2.0 * ts_xzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyz[i] = 2.0 * ts_xzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xzz[i] = 2.0 * ts_xzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyy[i] = 2.0 * ts_xzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyz[i] = 2.0 * ts_xzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yzz[i] = 2.0 * ts_xzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_zzz[i] = 2.0 * ts_xzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 240-250 components of targeted buffer : GF

    auto gs_y_xzzz_xxx = pbuffer.data(idx_g_gf + 240);

    auto gs_y_xzzz_xxy = pbuffer.data(idx_g_gf + 241);

    auto gs_y_xzzz_xxz = pbuffer.data(idx_g_gf + 242);

    auto gs_y_xzzz_xyy = pbuffer.data(idx_g_gf + 243);

    auto gs_y_xzzz_xyz = pbuffer.data(idx_g_gf + 244);

    auto gs_y_xzzz_xzz = pbuffer.data(idx_g_gf + 245);

    auto gs_y_xzzz_yyy = pbuffer.data(idx_g_gf + 246);

    auto gs_y_xzzz_yyz = pbuffer.data(idx_g_gf + 247);

    auto gs_y_xzzz_yzz = pbuffer.data(idx_g_gf + 248);

    auto gs_y_xzzz_zzz = pbuffer.data(idx_g_gf + 249);

    #pragma omp simd aligned(gc_y, gs_y_xzzz_xxx, gs_y_xzzz_xxy, gs_y_xzzz_xxz, gs_y_xzzz_xyy, gs_y_xzzz_xyz, gs_y_xzzz_xzz, gs_y_xzzz_yyy, gs_y_xzzz_yyz, gs_y_xzzz_yzz, gs_y_xzzz_zzz, ts_xzzz_xx, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xy, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xz, ts_xzzz_xzz, ts_xzzz_yy, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yz, ts_xzzz_yzz, ts_xzzz_zz, ts_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzz_xxx[i] = 2.0 * ts_xzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxy[i] = 2.0 * ts_xzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxz[i] = 2.0 * ts_xzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyy[i] = 4.0 * ts_xzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyz[i] = 2.0 * ts_xzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xzz[i] = 2.0 * ts_xzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyy[i] = 6.0 * ts_xzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyz[i] = 4.0 * ts_xzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yzz[i] = 2.0 * ts_xzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_zzz[i] = 2.0 * ts_xzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 250-260 components of targeted buffer : GF

    auto gs_y_yyyy_xxx = pbuffer.data(idx_g_gf + 250);

    auto gs_y_yyyy_xxy = pbuffer.data(idx_g_gf + 251);

    auto gs_y_yyyy_xxz = pbuffer.data(idx_g_gf + 252);

    auto gs_y_yyyy_xyy = pbuffer.data(idx_g_gf + 253);

    auto gs_y_yyyy_xyz = pbuffer.data(idx_g_gf + 254);

    auto gs_y_yyyy_xzz = pbuffer.data(idx_g_gf + 255);

    auto gs_y_yyyy_yyy = pbuffer.data(idx_g_gf + 256);

    auto gs_y_yyyy_yyz = pbuffer.data(idx_g_gf + 257);

    auto gs_y_yyyy_yzz = pbuffer.data(idx_g_gf + 258);

    auto gs_y_yyyy_zzz = pbuffer.data(idx_g_gf + 259);

    #pragma omp simd aligned(gc_y, gs_y_yyyy_xxx, gs_y_yyyy_xxy, gs_y_yyyy_xxz, gs_y_yyyy_xyy, gs_y_yyyy_xyz, gs_y_yyyy_xzz, gs_y_yyyy_yyy, gs_y_yyyy_yyz, gs_y_yyyy_yzz, gs_y_yyyy_zzz, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xzz, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yzz, ts_yyy_zzz, ts_yyyy_xx, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xy, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xz, ts_yyyy_xzz, ts_yyyy_yy, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yz, ts_yyyy_yzz, ts_yyyy_zz, ts_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyy_xxx[i] = 8.0 * ts_yyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxy[i] = 8.0 * ts_yyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxz[i] = 8.0 * ts_yyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyy[i] = 8.0 * ts_yyy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyz[i] = 8.0 * ts_yyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xzz[i] = 8.0 * ts_yyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyy[i] = 8.0 * ts_yyy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyz[i] = 8.0 * ts_yyy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yzz[i] = 8.0 * ts_yyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_zzz[i] = 8.0 * ts_yyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 260-270 components of targeted buffer : GF

    auto gs_y_yyyz_xxx = pbuffer.data(idx_g_gf + 260);

    auto gs_y_yyyz_xxy = pbuffer.data(idx_g_gf + 261);

    auto gs_y_yyyz_xxz = pbuffer.data(idx_g_gf + 262);

    auto gs_y_yyyz_xyy = pbuffer.data(idx_g_gf + 263);

    auto gs_y_yyyz_xyz = pbuffer.data(idx_g_gf + 264);

    auto gs_y_yyyz_xzz = pbuffer.data(idx_g_gf + 265);

    auto gs_y_yyyz_yyy = pbuffer.data(idx_g_gf + 266);

    auto gs_y_yyyz_yyz = pbuffer.data(idx_g_gf + 267);

    auto gs_y_yyyz_yzz = pbuffer.data(idx_g_gf + 268);

    auto gs_y_yyyz_zzz = pbuffer.data(idx_g_gf + 269);

    #pragma omp simd aligned(gc_y, gs_y_yyyz_xxx, gs_y_yyyz_xxy, gs_y_yyyz_xxz, gs_y_yyyz_xyy, gs_y_yyyz_xyz, gs_y_yyyz_xzz, gs_y_yyyz_yyy, gs_y_yyyz_yyz, gs_y_yyyz_yzz, gs_y_yyyz_zzz, ts_yyyz_xx, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xy, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xz, ts_yyyz_xzz, ts_yyyz_yy, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yz, ts_yyyz_yzz, ts_yyyz_zz, ts_yyyz_zzz, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xzz, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yzz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyz_xxx[i] = 6.0 * ts_yyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxy[i] = 6.0 * ts_yyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxz[i] = 6.0 * ts_yyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyy[i] = 6.0 * ts_yyz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyz[i] = 6.0 * ts_yyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xzz[i] = 6.0 * ts_yyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyy[i] = 6.0 * ts_yyz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyz[i] = 6.0 * ts_yyz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yzz[i] = 6.0 * ts_yyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_zzz[i] = 6.0 * ts_yyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 270-280 components of targeted buffer : GF

    auto gs_y_yyzz_xxx = pbuffer.data(idx_g_gf + 270);

    auto gs_y_yyzz_xxy = pbuffer.data(idx_g_gf + 271);

    auto gs_y_yyzz_xxz = pbuffer.data(idx_g_gf + 272);

    auto gs_y_yyzz_xyy = pbuffer.data(idx_g_gf + 273);

    auto gs_y_yyzz_xyz = pbuffer.data(idx_g_gf + 274);

    auto gs_y_yyzz_xzz = pbuffer.data(idx_g_gf + 275);

    auto gs_y_yyzz_yyy = pbuffer.data(idx_g_gf + 276);

    auto gs_y_yyzz_yyz = pbuffer.data(idx_g_gf + 277);

    auto gs_y_yyzz_yzz = pbuffer.data(idx_g_gf + 278);

    auto gs_y_yyzz_zzz = pbuffer.data(idx_g_gf + 279);

    #pragma omp simd aligned(gc_y, gs_y_yyzz_xxx, gs_y_yyzz_xxy, gs_y_yyzz_xxz, gs_y_yyzz_xyy, gs_y_yyzz_xyz, gs_y_yyzz_xzz, gs_y_yyzz_yyy, gs_y_yyzz_yyz, gs_y_yyzz_yzz, gs_y_yyzz_zzz, ts_yyzz_xx, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xy, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xz, ts_yyzz_xzz, ts_yyzz_yy, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yz, ts_yyzz_yzz, ts_yyzz_zz, ts_yyzz_zzz, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xzz, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yzz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzz_xxx[i] = 4.0 * ts_yzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxy[i] = 4.0 * ts_yzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxz[i] = 4.0 * ts_yzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyy[i] = 4.0 * ts_yzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyz[i] = 4.0 * ts_yzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xzz[i] = 4.0 * ts_yzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyy[i] = 4.0 * ts_yzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyz[i] = 4.0 * ts_yzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yzz[i] = 4.0 * ts_yzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_zzz[i] = 4.0 * ts_yzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 280-290 components of targeted buffer : GF

    auto gs_y_yzzz_xxx = pbuffer.data(idx_g_gf + 280);

    auto gs_y_yzzz_xxy = pbuffer.data(idx_g_gf + 281);

    auto gs_y_yzzz_xxz = pbuffer.data(idx_g_gf + 282);

    auto gs_y_yzzz_xyy = pbuffer.data(idx_g_gf + 283);

    auto gs_y_yzzz_xyz = pbuffer.data(idx_g_gf + 284);

    auto gs_y_yzzz_xzz = pbuffer.data(idx_g_gf + 285);

    auto gs_y_yzzz_yyy = pbuffer.data(idx_g_gf + 286);

    auto gs_y_yzzz_yyz = pbuffer.data(idx_g_gf + 287);

    auto gs_y_yzzz_yzz = pbuffer.data(idx_g_gf + 288);

    auto gs_y_yzzz_zzz = pbuffer.data(idx_g_gf + 289);

    #pragma omp simd aligned(gc_y, gs_y_yzzz_xxx, gs_y_yzzz_xxy, gs_y_yzzz_xxz, gs_y_yzzz_xyy, gs_y_yzzz_xyz, gs_y_yzzz_xzz, gs_y_yzzz_yyy, gs_y_yzzz_yyz, gs_y_yzzz_yzz, gs_y_yzzz_zzz, ts_yzzz_xx, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xy, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xz, ts_yzzz_xzz, ts_yzzz_yy, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yz, ts_yzzz_yzz, ts_yzzz_zz, ts_yzzz_zzz, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xzz, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yzz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzz_xxx[i] = 2.0 * ts_zzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxy[i] = 2.0 * ts_zzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxz[i] = 2.0 * ts_zzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyy[i] = 2.0 * ts_zzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyz[i] = 2.0 * ts_zzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xzz[i] = 2.0 * ts_zzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyy[i] = 2.0 * ts_zzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyz[i] = 2.0 * ts_zzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yzz[i] = 2.0 * ts_zzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 290-300 components of targeted buffer : GF

    auto gs_y_zzzz_xxx = pbuffer.data(idx_g_gf + 290);

    auto gs_y_zzzz_xxy = pbuffer.data(idx_g_gf + 291);

    auto gs_y_zzzz_xxz = pbuffer.data(idx_g_gf + 292);

    auto gs_y_zzzz_xyy = pbuffer.data(idx_g_gf + 293);

    auto gs_y_zzzz_xyz = pbuffer.data(idx_g_gf + 294);

    auto gs_y_zzzz_xzz = pbuffer.data(idx_g_gf + 295);

    auto gs_y_zzzz_yyy = pbuffer.data(idx_g_gf + 296);

    auto gs_y_zzzz_yyz = pbuffer.data(idx_g_gf + 297);

    auto gs_y_zzzz_yzz = pbuffer.data(idx_g_gf + 298);

    auto gs_y_zzzz_zzz = pbuffer.data(idx_g_gf + 299);

    #pragma omp simd aligned(gc_y, gs_y_zzzz_xxx, gs_y_zzzz_xxy, gs_y_zzzz_xxz, gs_y_zzzz_xyy, gs_y_zzzz_xyz, gs_y_zzzz_xzz, gs_y_zzzz_yyy, gs_y_zzzz_yyz, gs_y_zzzz_yzz, gs_y_zzzz_zzz, ts_zzzz_xx, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xy, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xz, ts_zzzz_xzz, ts_zzzz_yy, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yz, ts_zzzz_yzz, ts_zzzz_zz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzz_xxx[i] = 2.0 * ts_zzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxy[i] = 2.0 * ts_zzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxz[i] = 2.0 * ts_zzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyy[i] = 4.0 * ts_zzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyz[i] = 2.0 * ts_zzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xzz[i] = 2.0 * ts_zzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyy[i] = 6.0 * ts_zzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyz[i] = 4.0 * ts_zzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yzz[i] = 2.0 * ts_zzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_zzz[i] = 2.0 * ts_zzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 300-310 components of targeted buffer : GF

    auto gs_z_xxxx_xxx = pbuffer.data(idx_g_gf + 300);

    auto gs_z_xxxx_xxy = pbuffer.data(idx_g_gf + 301);

    auto gs_z_xxxx_xxz = pbuffer.data(idx_g_gf + 302);

    auto gs_z_xxxx_xyy = pbuffer.data(idx_g_gf + 303);

    auto gs_z_xxxx_xyz = pbuffer.data(idx_g_gf + 304);

    auto gs_z_xxxx_xzz = pbuffer.data(idx_g_gf + 305);

    auto gs_z_xxxx_yyy = pbuffer.data(idx_g_gf + 306);

    auto gs_z_xxxx_yyz = pbuffer.data(idx_g_gf + 307);

    auto gs_z_xxxx_yzz = pbuffer.data(idx_g_gf + 308);

    auto gs_z_xxxx_zzz = pbuffer.data(idx_g_gf + 309);

    #pragma omp simd aligned(gc_z, gs_z_xxxx_xxx, gs_z_xxxx_xxy, gs_z_xxxx_xxz, gs_z_xxxx_xyy, gs_z_xxxx_xyz, gs_z_xxxx_xzz, gs_z_xxxx_yyy, gs_z_xxxx_yyz, gs_z_xxxx_yzz, gs_z_xxxx_zzz, ts_xxxx_xx, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xy, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xz, ts_xxxx_xzz, ts_xxxx_yy, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yz, ts_xxxx_yzz, ts_xxxx_zz, ts_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxx_xxx[i] = 2.0 * ts_xxxx_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxy[i] = 2.0 * ts_xxxx_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxz[i] = 2.0 * ts_xxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyy[i] = 2.0 * ts_xxxx_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyz[i] = 2.0 * ts_xxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xzz[i] = 4.0 * ts_xxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyy[i] = 2.0 * ts_xxxx_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyz[i] = 2.0 * ts_xxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yzz[i] = 4.0 * ts_xxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_zzz[i] = 6.0 * ts_xxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 310-320 components of targeted buffer : GF

    auto gs_z_xxxy_xxx = pbuffer.data(idx_g_gf + 310);

    auto gs_z_xxxy_xxy = pbuffer.data(idx_g_gf + 311);

    auto gs_z_xxxy_xxz = pbuffer.data(idx_g_gf + 312);

    auto gs_z_xxxy_xyy = pbuffer.data(idx_g_gf + 313);

    auto gs_z_xxxy_xyz = pbuffer.data(idx_g_gf + 314);

    auto gs_z_xxxy_xzz = pbuffer.data(idx_g_gf + 315);

    auto gs_z_xxxy_yyy = pbuffer.data(idx_g_gf + 316);

    auto gs_z_xxxy_yyz = pbuffer.data(idx_g_gf + 317);

    auto gs_z_xxxy_yzz = pbuffer.data(idx_g_gf + 318);

    auto gs_z_xxxy_zzz = pbuffer.data(idx_g_gf + 319);

    #pragma omp simd aligned(gc_z, gs_z_xxxy_xxx, gs_z_xxxy_xxy, gs_z_xxxy_xxz, gs_z_xxxy_xyy, gs_z_xxxy_xyz, gs_z_xxxy_xzz, gs_z_xxxy_yyy, gs_z_xxxy_yyz, gs_z_xxxy_yzz, gs_z_xxxy_zzz, ts_xxxy_xx, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xy, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xz, ts_xxxy_xzz, ts_xxxy_yy, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yz, ts_xxxy_yzz, ts_xxxy_zz, ts_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxy_xxx[i] = 2.0 * ts_xxxy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxy[i] = 2.0 * ts_xxxy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxz[i] = 2.0 * ts_xxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyy[i] = 2.0 * ts_xxxy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyz[i] = 2.0 * ts_xxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xzz[i] = 4.0 * ts_xxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyy[i] = 2.0 * ts_xxxy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyz[i] = 2.0 * ts_xxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yzz[i] = 4.0 * ts_xxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_zzz[i] = 6.0 * ts_xxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 320-330 components of targeted buffer : GF

    auto gs_z_xxxz_xxx = pbuffer.data(idx_g_gf + 320);

    auto gs_z_xxxz_xxy = pbuffer.data(idx_g_gf + 321);

    auto gs_z_xxxz_xxz = pbuffer.data(idx_g_gf + 322);

    auto gs_z_xxxz_xyy = pbuffer.data(idx_g_gf + 323);

    auto gs_z_xxxz_xyz = pbuffer.data(idx_g_gf + 324);

    auto gs_z_xxxz_xzz = pbuffer.data(idx_g_gf + 325);

    auto gs_z_xxxz_yyy = pbuffer.data(idx_g_gf + 326);

    auto gs_z_xxxz_yyz = pbuffer.data(idx_g_gf + 327);

    auto gs_z_xxxz_yzz = pbuffer.data(idx_g_gf + 328);

    auto gs_z_xxxz_zzz = pbuffer.data(idx_g_gf + 329);

    #pragma omp simd aligned(gc_z, gs_z_xxxz_xxx, gs_z_xxxz_xxy, gs_z_xxxz_xxz, gs_z_xxxz_xyy, gs_z_xxxz_xyz, gs_z_xxxz_xzz, gs_z_xxxz_yyy, gs_z_xxxz_yyz, gs_z_xxxz_yzz, gs_z_xxxz_zzz, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xzz, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yzz, ts_xxx_zzz, ts_xxxz_xx, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xy, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xz, ts_xxxz_xzz, ts_xxxz_yy, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yz, ts_xxxz_yzz, ts_xxxz_zz, ts_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxz_xxx[i] = 2.0 * ts_xxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxy[i] = 2.0 * ts_xxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxz[i] = 2.0 * ts_xxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyy[i] = 2.0 * ts_xxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyz[i] = 2.0 * ts_xxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xzz[i] = 2.0 * ts_xxx_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyy[i] = 2.0 * ts_xxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyz[i] = 2.0 * ts_xxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yzz[i] = 2.0 * ts_xxx_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_zzz[i] = 2.0 * ts_xxx_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 330-340 components of targeted buffer : GF

    auto gs_z_xxyy_xxx = pbuffer.data(idx_g_gf + 330);

    auto gs_z_xxyy_xxy = pbuffer.data(idx_g_gf + 331);

    auto gs_z_xxyy_xxz = pbuffer.data(idx_g_gf + 332);

    auto gs_z_xxyy_xyy = pbuffer.data(idx_g_gf + 333);

    auto gs_z_xxyy_xyz = pbuffer.data(idx_g_gf + 334);

    auto gs_z_xxyy_xzz = pbuffer.data(idx_g_gf + 335);

    auto gs_z_xxyy_yyy = pbuffer.data(idx_g_gf + 336);

    auto gs_z_xxyy_yyz = pbuffer.data(idx_g_gf + 337);

    auto gs_z_xxyy_yzz = pbuffer.data(idx_g_gf + 338);

    auto gs_z_xxyy_zzz = pbuffer.data(idx_g_gf + 339);

    #pragma omp simd aligned(gc_z, gs_z_xxyy_xxx, gs_z_xxyy_xxy, gs_z_xxyy_xxz, gs_z_xxyy_xyy, gs_z_xxyy_xyz, gs_z_xxyy_xzz, gs_z_xxyy_yyy, gs_z_xxyy_yyz, gs_z_xxyy_yzz, gs_z_xxyy_zzz, ts_xxyy_xx, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xy, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xz, ts_xxyy_xzz, ts_xxyy_yy, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yz, ts_xxyy_yzz, ts_xxyy_zz, ts_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyy_xxx[i] = 2.0 * ts_xxyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxy[i] = 2.0 * ts_xxyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxz[i] = 2.0 * ts_xxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyy[i] = 2.0 * ts_xxyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyz[i] = 2.0 * ts_xxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xzz[i] = 4.0 * ts_xxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyy[i] = 2.0 * ts_xxyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyz[i] = 2.0 * ts_xxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yzz[i] = 4.0 * ts_xxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_zzz[i] = 6.0 * ts_xxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 340-350 components of targeted buffer : GF

    auto gs_z_xxyz_xxx = pbuffer.data(idx_g_gf + 340);

    auto gs_z_xxyz_xxy = pbuffer.data(idx_g_gf + 341);

    auto gs_z_xxyz_xxz = pbuffer.data(idx_g_gf + 342);

    auto gs_z_xxyz_xyy = pbuffer.data(idx_g_gf + 343);

    auto gs_z_xxyz_xyz = pbuffer.data(idx_g_gf + 344);

    auto gs_z_xxyz_xzz = pbuffer.data(idx_g_gf + 345);

    auto gs_z_xxyz_yyy = pbuffer.data(idx_g_gf + 346);

    auto gs_z_xxyz_yyz = pbuffer.data(idx_g_gf + 347);

    auto gs_z_xxyz_yzz = pbuffer.data(idx_g_gf + 348);

    auto gs_z_xxyz_zzz = pbuffer.data(idx_g_gf + 349);

    #pragma omp simd aligned(gc_z, gs_z_xxyz_xxx, gs_z_xxyz_xxy, gs_z_xxyz_xxz, gs_z_xxyz_xyy, gs_z_xxyz_xyz, gs_z_xxyz_xzz, gs_z_xxyz_yyy, gs_z_xxyz_yyz, gs_z_xxyz_yzz, gs_z_xxyz_zzz, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xzz, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yzz, ts_xxy_zzz, ts_xxyz_xx, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xy, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xz, ts_xxyz_xzz, ts_xxyz_yy, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yz, ts_xxyz_yzz, ts_xxyz_zz, ts_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyz_xxx[i] = 2.0 * ts_xxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxy[i] = 2.0 * ts_xxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxz[i] = 2.0 * ts_xxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyy[i] = 2.0 * ts_xxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyz[i] = 2.0 * ts_xxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xzz[i] = 2.0 * ts_xxy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyy[i] = 2.0 * ts_xxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyz[i] = 2.0 * ts_xxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yzz[i] = 2.0 * ts_xxy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_zzz[i] = 2.0 * ts_xxy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 350-360 components of targeted buffer : GF

    auto gs_z_xxzz_xxx = pbuffer.data(idx_g_gf + 350);

    auto gs_z_xxzz_xxy = pbuffer.data(idx_g_gf + 351);

    auto gs_z_xxzz_xxz = pbuffer.data(idx_g_gf + 352);

    auto gs_z_xxzz_xyy = pbuffer.data(idx_g_gf + 353);

    auto gs_z_xxzz_xyz = pbuffer.data(idx_g_gf + 354);

    auto gs_z_xxzz_xzz = pbuffer.data(idx_g_gf + 355);

    auto gs_z_xxzz_yyy = pbuffer.data(idx_g_gf + 356);

    auto gs_z_xxzz_yyz = pbuffer.data(idx_g_gf + 357);

    auto gs_z_xxzz_yzz = pbuffer.data(idx_g_gf + 358);

    auto gs_z_xxzz_zzz = pbuffer.data(idx_g_gf + 359);

    #pragma omp simd aligned(gc_z, gs_z_xxzz_xxx, gs_z_xxzz_xxy, gs_z_xxzz_xxz, gs_z_xxzz_xyy, gs_z_xxzz_xyz, gs_z_xxzz_xzz, gs_z_xxzz_yyy, gs_z_xxzz_yyz, gs_z_xxzz_yzz, gs_z_xxzz_zzz, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xzz, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yzz, ts_xxz_zzz, ts_xxzz_xx, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xy, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xz, ts_xxzz_xzz, ts_xxzz_yy, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yz, ts_xxzz_yzz, ts_xxzz_zz, ts_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzz_xxx[i] = 4.0 * ts_xxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxy[i] = 4.0 * ts_xxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxz[i] = 4.0 * ts_xxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyy[i] = 4.0 * ts_xxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyz[i] = 4.0 * ts_xxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xzz[i] = 4.0 * ts_xxz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyy[i] = 4.0 * ts_xxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyz[i] = 4.0 * ts_xxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yzz[i] = 4.0 * ts_xxz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_zzz[i] = 4.0 * ts_xxz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 360-370 components of targeted buffer : GF

    auto gs_z_xyyy_xxx = pbuffer.data(idx_g_gf + 360);

    auto gs_z_xyyy_xxy = pbuffer.data(idx_g_gf + 361);

    auto gs_z_xyyy_xxz = pbuffer.data(idx_g_gf + 362);

    auto gs_z_xyyy_xyy = pbuffer.data(idx_g_gf + 363);

    auto gs_z_xyyy_xyz = pbuffer.data(idx_g_gf + 364);

    auto gs_z_xyyy_xzz = pbuffer.data(idx_g_gf + 365);

    auto gs_z_xyyy_yyy = pbuffer.data(idx_g_gf + 366);

    auto gs_z_xyyy_yyz = pbuffer.data(idx_g_gf + 367);

    auto gs_z_xyyy_yzz = pbuffer.data(idx_g_gf + 368);

    auto gs_z_xyyy_zzz = pbuffer.data(idx_g_gf + 369);

    #pragma omp simd aligned(gc_z, gs_z_xyyy_xxx, gs_z_xyyy_xxy, gs_z_xyyy_xxz, gs_z_xyyy_xyy, gs_z_xyyy_xyz, gs_z_xyyy_xzz, gs_z_xyyy_yyy, gs_z_xyyy_yyz, gs_z_xyyy_yzz, gs_z_xyyy_zzz, ts_xyyy_xx, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xy, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xz, ts_xyyy_xzz, ts_xyyy_yy, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yz, ts_xyyy_yzz, ts_xyyy_zz, ts_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyy_xxx[i] = 2.0 * ts_xyyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxy[i] = 2.0 * ts_xyyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxz[i] = 2.0 * ts_xyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyy[i] = 2.0 * ts_xyyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyz[i] = 2.0 * ts_xyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xzz[i] = 4.0 * ts_xyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyy[i] = 2.0 * ts_xyyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyz[i] = 2.0 * ts_xyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yzz[i] = 4.0 * ts_xyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_zzz[i] = 6.0 * ts_xyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 370-380 components of targeted buffer : GF

    auto gs_z_xyyz_xxx = pbuffer.data(idx_g_gf + 370);

    auto gs_z_xyyz_xxy = pbuffer.data(idx_g_gf + 371);

    auto gs_z_xyyz_xxz = pbuffer.data(idx_g_gf + 372);

    auto gs_z_xyyz_xyy = pbuffer.data(idx_g_gf + 373);

    auto gs_z_xyyz_xyz = pbuffer.data(idx_g_gf + 374);

    auto gs_z_xyyz_xzz = pbuffer.data(idx_g_gf + 375);

    auto gs_z_xyyz_yyy = pbuffer.data(idx_g_gf + 376);

    auto gs_z_xyyz_yyz = pbuffer.data(idx_g_gf + 377);

    auto gs_z_xyyz_yzz = pbuffer.data(idx_g_gf + 378);

    auto gs_z_xyyz_zzz = pbuffer.data(idx_g_gf + 379);

    #pragma omp simd aligned(gc_z, gs_z_xyyz_xxx, gs_z_xyyz_xxy, gs_z_xyyz_xxz, gs_z_xyyz_xyy, gs_z_xyyz_xyz, gs_z_xyyz_xzz, gs_z_xyyz_yyy, gs_z_xyyz_yyz, gs_z_xyyz_yzz, gs_z_xyyz_zzz, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xzz, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yzz, ts_xyy_zzz, ts_xyyz_xx, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xy, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xz, ts_xyyz_xzz, ts_xyyz_yy, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yz, ts_xyyz_yzz, ts_xyyz_zz, ts_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyz_xxx[i] = 2.0 * ts_xyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxy[i] = 2.0 * ts_xyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxz[i] = 2.0 * ts_xyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyy[i] = 2.0 * ts_xyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyz[i] = 2.0 * ts_xyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xzz[i] = 2.0 * ts_xyy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyy[i] = 2.0 * ts_xyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyz[i] = 2.0 * ts_xyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yzz[i] = 2.0 * ts_xyy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_zzz[i] = 2.0 * ts_xyy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 380-390 components of targeted buffer : GF

    auto gs_z_xyzz_xxx = pbuffer.data(idx_g_gf + 380);

    auto gs_z_xyzz_xxy = pbuffer.data(idx_g_gf + 381);

    auto gs_z_xyzz_xxz = pbuffer.data(idx_g_gf + 382);

    auto gs_z_xyzz_xyy = pbuffer.data(idx_g_gf + 383);

    auto gs_z_xyzz_xyz = pbuffer.data(idx_g_gf + 384);

    auto gs_z_xyzz_xzz = pbuffer.data(idx_g_gf + 385);

    auto gs_z_xyzz_yyy = pbuffer.data(idx_g_gf + 386);

    auto gs_z_xyzz_yyz = pbuffer.data(idx_g_gf + 387);

    auto gs_z_xyzz_yzz = pbuffer.data(idx_g_gf + 388);

    auto gs_z_xyzz_zzz = pbuffer.data(idx_g_gf + 389);

    #pragma omp simd aligned(gc_z, gs_z_xyzz_xxx, gs_z_xyzz_xxy, gs_z_xyzz_xxz, gs_z_xyzz_xyy, gs_z_xyzz_xyz, gs_z_xyzz_xzz, gs_z_xyzz_yyy, gs_z_xyzz_yyz, gs_z_xyzz_yzz, gs_z_xyzz_zzz, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xzz, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yzz, ts_xyz_zzz, ts_xyzz_xx, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xy, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xz, ts_xyzz_xzz, ts_xyzz_yy, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yz, ts_xyzz_yzz, ts_xyzz_zz, ts_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzz_xxx[i] = 4.0 * ts_xyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxy[i] = 4.0 * ts_xyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxz[i] = 4.0 * ts_xyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyy[i] = 4.0 * ts_xyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyz[i] = 4.0 * ts_xyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xzz[i] = 4.0 * ts_xyz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyy[i] = 4.0 * ts_xyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyz[i] = 4.0 * ts_xyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yzz[i] = 4.0 * ts_xyz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_zzz[i] = 4.0 * ts_xyz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 390-400 components of targeted buffer : GF

    auto gs_z_xzzz_xxx = pbuffer.data(idx_g_gf + 390);

    auto gs_z_xzzz_xxy = pbuffer.data(idx_g_gf + 391);

    auto gs_z_xzzz_xxz = pbuffer.data(idx_g_gf + 392);

    auto gs_z_xzzz_xyy = pbuffer.data(idx_g_gf + 393);

    auto gs_z_xzzz_xyz = pbuffer.data(idx_g_gf + 394);

    auto gs_z_xzzz_xzz = pbuffer.data(idx_g_gf + 395);

    auto gs_z_xzzz_yyy = pbuffer.data(idx_g_gf + 396);

    auto gs_z_xzzz_yyz = pbuffer.data(idx_g_gf + 397);

    auto gs_z_xzzz_yzz = pbuffer.data(idx_g_gf + 398);

    auto gs_z_xzzz_zzz = pbuffer.data(idx_g_gf + 399);

    #pragma omp simd aligned(gc_z, gs_z_xzzz_xxx, gs_z_xzzz_xxy, gs_z_xzzz_xxz, gs_z_xzzz_xyy, gs_z_xzzz_xyz, gs_z_xzzz_xzz, gs_z_xzzz_yyy, gs_z_xzzz_yyz, gs_z_xzzz_yzz, gs_z_xzzz_zzz, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xzz, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yzz, ts_xzz_zzz, ts_xzzz_xx, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xy, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xz, ts_xzzz_xzz, ts_xzzz_yy, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yz, ts_xzzz_yzz, ts_xzzz_zz, ts_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzz_xxx[i] = 6.0 * ts_xzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxy[i] = 6.0 * ts_xzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxz[i] = 6.0 * ts_xzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyy[i] = 6.0 * ts_xzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyz[i] = 6.0 * ts_xzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xzz[i] = 6.0 * ts_xzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyy[i] = 6.0 * ts_xzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyz[i] = 6.0 * ts_xzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yzz[i] = 6.0 * ts_xzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_zzz[i] = 6.0 * ts_xzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 400-410 components of targeted buffer : GF

    auto gs_z_yyyy_xxx = pbuffer.data(idx_g_gf + 400);

    auto gs_z_yyyy_xxy = pbuffer.data(idx_g_gf + 401);

    auto gs_z_yyyy_xxz = pbuffer.data(idx_g_gf + 402);

    auto gs_z_yyyy_xyy = pbuffer.data(idx_g_gf + 403);

    auto gs_z_yyyy_xyz = pbuffer.data(idx_g_gf + 404);

    auto gs_z_yyyy_xzz = pbuffer.data(idx_g_gf + 405);

    auto gs_z_yyyy_yyy = pbuffer.data(idx_g_gf + 406);

    auto gs_z_yyyy_yyz = pbuffer.data(idx_g_gf + 407);

    auto gs_z_yyyy_yzz = pbuffer.data(idx_g_gf + 408);

    auto gs_z_yyyy_zzz = pbuffer.data(idx_g_gf + 409);

    #pragma omp simd aligned(gc_z, gs_z_yyyy_xxx, gs_z_yyyy_xxy, gs_z_yyyy_xxz, gs_z_yyyy_xyy, gs_z_yyyy_xyz, gs_z_yyyy_xzz, gs_z_yyyy_yyy, gs_z_yyyy_yyz, gs_z_yyyy_yzz, gs_z_yyyy_zzz, ts_yyyy_xx, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xy, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xz, ts_yyyy_xzz, ts_yyyy_yy, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yz, ts_yyyy_yzz, ts_yyyy_zz, ts_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyy_xxx[i] = 2.0 * ts_yyyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxy[i] = 2.0 * ts_yyyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxz[i] = 2.0 * ts_yyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyy[i] = 2.0 * ts_yyyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyz[i] = 2.0 * ts_yyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xzz[i] = 4.0 * ts_yyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyy[i] = 2.0 * ts_yyyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyz[i] = 2.0 * ts_yyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yzz[i] = 4.0 * ts_yyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_zzz[i] = 6.0 * ts_yyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 410-420 components of targeted buffer : GF

    auto gs_z_yyyz_xxx = pbuffer.data(idx_g_gf + 410);

    auto gs_z_yyyz_xxy = pbuffer.data(idx_g_gf + 411);

    auto gs_z_yyyz_xxz = pbuffer.data(idx_g_gf + 412);

    auto gs_z_yyyz_xyy = pbuffer.data(idx_g_gf + 413);

    auto gs_z_yyyz_xyz = pbuffer.data(idx_g_gf + 414);

    auto gs_z_yyyz_xzz = pbuffer.data(idx_g_gf + 415);

    auto gs_z_yyyz_yyy = pbuffer.data(idx_g_gf + 416);

    auto gs_z_yyyz_yyz = pbuffer.data(idx_g_gf + 417);

    auto gs_z_yyyz_yzz = pbuffer.data(idx_g_gf + 418);

    auto gs_z_yyyz_zzz = pbuffer.data(idx_g_gf + 419);

    #pragma omp simd aligned(gc_z, gs_z_yyyz_xxx, gs_z_yyyz_xxy, gs_z_yyyz_xxz, gs_z_yyyz_xyy, gs_z_yyyz_xyz, gs_z_yyyz_xzz, gs_z_yyyz_yyy, gs_z_yyyz_yyz, gs_z_yyyz_yzz, gs_z_yyyz_zzz, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xzz, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yzz, ts_yyy_zzz, ts_yyyz_xx, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xy, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xz, ts_yyyz_xzz, ts_yyyz_yy, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yz, ts_yyyz_yzz, ts_yyyz_zz, ts_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyz_xxx[i] = 2.0 * ts_yyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxy[i] = 2.0 * ts_yyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxz[i] = 2.0 * ts_yyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyy[i] = 2.0 * ts_yyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyz[i] = 2.0 * ts_yyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xzz[i] = 2.0 * ts_yyy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyy[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyz[i] = 2.0 * ts_yyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yzz[i] = 2.0 * ts_yyy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_zzz[i] = 2.0 * ts_yyy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 420-430 components of targeted buffer : GF

    auto gs_z_yyzz_xxx = pbuffer.data(idx_g_gf + 420);

    auto gs_z_yyzz_xxy = pbuffer.data(idx_g_gf + 421);

    auto gs_z_yyzz_xxz = pbuffer.data(idx_g_gf + 422);

    auto gs_z_yyzz_xyy = pbuffer.data(idx_g_gf + 423);

    auto gs_z_yyzz_xyz = pbuffer.data(idx_g_gf + 424);

    auto gs_z_yyzz_xzz = pbuffer.data(idx_g_gf + 425);

    auto gs_z_yyzz_yyy = pbuffer.data(idx_g_gf + 426);

    auto gs_z_yyzz_yyz = pbuffer.data(idx_g_gf + 427);

    auto gs_z_yyzz_yzz = pbuffer.data(idx_g_gf + 428);

    auto gs_z_yyzz_zzz = pbuffer.data(idx_g_gf + 429);

    #pragma omp simd aligned(gc_z, gs_z_yyzz_xxx, gs_z_yyzz_xxy, gs_z_yyzz_xxz, gs_z_yyzz_xyy, gs_z_yyzz_xyz, gs_z_yyzz_xzz, gs_z_yyzz_yyy, gs_z_yyzz_yyz, gs_z_yyzz_yzz, gs_z_yyzz_zzz, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xzz, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yzz, ts_yyz_zzz, ts_yyzz_xx, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xy, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xz, ts_yyzz_xzz, ts_yyzz_yy, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yz, ts_yyzz_yzz, ts_yyzz_zz, ts_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzz_xxx[i] = 4.0 * ts_yyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxy[i] = 4.0 * ts_yyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxz[i] = 4.0 * ts_yyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyy[i] = 4.0 * ts_yyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyz[i] = 4.0 * ts_yyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xzz[i] = 4.0 * ts_yyz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyy[i] = 4.0 * ts_yyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyz[i] = 4.0 * ts_yyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yzz[i] = 4.0 * ts_yyz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_zzz[i] = 4.0 * ts_yyz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 430-440 components of targeted buffer : GF

    auto gs_z_yzzz_xxx = pbuffer.data(idx_g_gf + 430);

    auto gs_z_yzzz_xxy = pbuffer.data(idx_g_gf + 431);

    auto gs_z_yzzz_xxz = pbuffer.data(idx_g_gf + 432);

    auto gs_z_yzzz_xyy = pbuffer.data(idx_g_gf + 433);

    auto gs_z_yzzz_xyz = pbuffer.data(idx_g_gf + 434);

    auto gs_z_yzzz_xzz = pbuffer.data(idx_g_gf + 435);

    auto gs_z_yzzz_yyy = pbuffer.data(idx_g_gf + 436);

    auto gs_z_yzzz_yyz = pbuffer.data(idx_g_gf + 437);

    auto gs_z_yzzz_yzz = pbuffer.data(idx_g_gf + 438);

    auto gs_z_yzzz_zzz = pbuffer.data(idx_g_gf + 439);

    #pragma omp simd aligned(gc_z, gs_z_yzzz_xxx, gs_z_yzzz_xxy, gs_z_yzzz_xxz, gs_z_yzzz_xyy, gs_z_yzzz_xyz, gs_z_yzzz_xzz, gs_z_yzzz_yyy, gs_z_yzzz_yyz, gs_z_yzzz_yzz, gs_z_yzzz_zzz, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xzz, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yzz, ts_yzz_zzz, ts_yzzz_xx, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xy, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xz, ts_yzzz_xzz, ts_yzzz_yy, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yz, ts_yzzz_yzz, ts_yzzz_zz, ts_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzz_xxx[i] = 6.0 * ts_yzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxy[i] = 6.0 * ts_yzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxz[i] = 6.0 * ts_yzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyy[i] = 6.0 * ts_yzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyz[i] = 6.0 * ts_yzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xzz[i] = 6.0 * ts_yzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyy[i] = 6.0 * ts_yzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyz[i] = 6.0 * ts_yzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yzz[i] = 6.0 * ts_yzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_zzz[i] = 6.0 * ts_yzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 440-450 components of targeted buffer : GF

    auto gs_z_zzzz_xxx = pbuffer.data(idx_g_gf + 440);

    auto gs_z_zzzz_xxy = pbuffer.data(idx_g_gf + 441);

    auto gs_z_zzzz_xxz = pbuffer.data(idx_g_gf + 442);

    auto gs_z_zzzz_xyy = pbuffer.data(idx_g_gf + 443);

    auto gs_z_zzzz_xyz = pbuffer.data(idx_g_gf + 444);

    auto gs_z_zzzz_xzz = pbuffer.data(idx_g_gf + 445);

    auto gs_z_zzzz_yyy = pbuffer.data(idx_g_gf + 446);

    auto gs_z_zzzz_yyz = pbuffer.data(idx_g_gf + 447);

    auto gs_z_zzzz_yzz = pbuffer.data(idx_g_gf + 448);

    auto gs_z_zzzz_zzz = pbuffer.data(idx_g_gf + 449);

    #pragma omp simd aligned(gc_z, gs_z_zzzz_xxx, gs_z_zzzz_xxy, gs_z_zzzz_xxz, gs_z_zzzz_xyy, gs_z_zzzz_xyz, gs_z_zzzz_xzz, gs_z_zzzz_yyy, gs_z_zzzz_yyz, gs_z_zzzz_yzz, gs_z_zzzz_zzz, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xzz, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yzz, ts_zzz_zzz, ts_zzzz_xx, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xy, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xz, ts_zzzz_xzz, ts_zzzz_yy, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yz, ts_zzzz_yzz, ts_zzzz_zz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzz_xxx[i] = 8.0 * ts_zzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxy[i] = 8.0 * ts_zzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxz[i] = 8.0 * ts_zzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyy[i] = 8.0 * ts_zzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyz[i] = 8.0 * ts_zzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xzz[i] = 8.0 * ts_zzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyy[i] = 8.0 * ts_zzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyz[i] = 8.0 * ts_zzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yzz[i] = 8.0 * ts_zzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_zzz[i] = 8.0 * ts_zzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_zzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

