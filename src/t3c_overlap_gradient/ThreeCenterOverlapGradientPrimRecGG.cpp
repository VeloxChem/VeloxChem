#include "ThreeCenterOverlapGradientPrimRecGG.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_gg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_gg,
                              const size_t idx_fg,
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

    auto gs_x_xxxx_xxxx = pbuffer.data(idx_g_gg);

    auto gs_x_xxxx_xxxy = pbuffer.data(idx_g_gg + 1);

    auto gs_x_xxxx_xxxz = pbuffer.data(idx_g_gg + 2);

    auto gs_x_xxxx_xxyy = pbuffer.data(idx_g_gg + 3);

    auto gs_x_xxxx_xxyz = pbuffer.data(idx_g_gg + 4);

    auto gs_x_xxxx_xxzz = pbuffer.data(idx_g_gg + 5);

    auto gs_x_xxxx_xyyy = pbuffer.data(idx_g_gg + 6);

    auto gs_x_xxxx_xyyz = pbuffer.data(idx_g_gg + 7);

    auto gs_x_xxxx_xyzz = pbuffer.data(idx_g_gg + 8);

    auto gs_x_xxxx_xzzz = pbuffer.data(idx_g_gg + 9);

    auto gs_x_xxxx_yyyy = pbuffer.data(idx_g_gg + 10);

    auto gs_x_xxxx_yyyz = pbuffer.data(idx_g_gg + 11);

    auto gs_x_xxxx_yyzz = pbuffer.data(idx_g_gg + 12);

    auto gs_x_xxxx_yzzz = pbuffer.data(idx_g_gg + 13);

    auto gs_x_xxxx_zzzz = pbuffer.data(idx_g_gg + 14);

    #pragma omp simd aligned(gc_x, gs_x_xxxx_xxxx, gs_x_xxxx_xxxy, gs_x_xxxx_xxxz, gs_x_xxxx_xxyy, gs_x_xxxx_xxyz, gs_x_xxxx_xxzz, gs_x_xxxx_xyyy, gs_x_xxxx_xyyz, gs_x_xxxx_xyzz, gs_x_xxxx_xzzz, gs_x_xxxx_yyyy, gs_x_xxxx_yyyz, gs_x_xxxx_yyzz, gs_x_xxxx_yzzz, gs_x_xxxx_zzzz, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxzz, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyzz, ts_xxx_xzzz, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyzz, ts_xxx_yzzz, ts_xxx_zzzz, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_yyy, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyz, ts_xxxx_yyzz, ts_xxxx_yzz, ts_xxxx_yzzz, ts_xxxx_zzz, ts_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxx_xxxx[i] = 8.0 * ts_xxx_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxy[i] = 8.0 * ts_xxx_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxxz[i] = 8.0 * ts_xxx_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxyy[i] = 8.0 * ts_xxx_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxyz[i] = 8.0 * ts_xxx_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xxzz[i] = 8.0 * ts_xxx_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyyy[i] = 8.0 * ts_xxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyyz[i] = 8.0 * ts_xxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xyzz[i] = 8.0 * ts_xxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xzzz[i] = 8.0 * ts_xxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyyy[i] = 8.0 * ts_xxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyyz[i] = 8.0 * ts_xxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yyzz[i] = 8.0 * ts_xxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yzzz[i] = 8.0 * ts_xxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_zzzz[i] = 8.0 * ts_xxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 15-30 components of targeted buffer : GG

    auto gs_x_xxxy_xxxx = pbuffer.data(idx_g_gg + 15);

    auto gs_x_xxxy_xxxy = pbuffer.data(idx_g_gg + 16);

    auto gs_x_xxxy_xxxz = pbuffer.data(idx_g_gg + 17);

    auto gs_x_xxxy_xxyy = pbuffer.data(idx_g_gg + 18);

    auto gs_x_xxxy_xxyz = pbuffer.data(idx_g_gg + 19);

    auto gs_x_xxxy_xxzz = pbuffer.data(idx_g_gg + 20);

    auto gs_x_xxxy_xyyy = pbuffer.data(idx_g_gg + 21);

    auto gs_x_xxxy_xyyz = pbuffer.data(idx_g_gg + 22);

    auto gs_x_xxxy_xyzz = pbuffer.data(idx_g_gg + 23);

    auto gs_x_xxxy_xzzz = pbuffer.data(idx_g_gg + 24);

    auto gs_x_xxxy_yyyy = pbuffer.data(idx_g_gg + 25);

    auto gs_x_xxxy_yyyz = pbuffer.data(idx_g_gg + 26);

    auto gs_x_xxxy_yyzz = pbuffer.data(idx_g_gg + 27);

    auto gs_x_xxxy_yzzz = pbuffer.data(idx_g_gg + 28);

    auto gs_x_xxxy_zzzz = pbuffer.data(idx_g_gg + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxxy_xxxx, gs_x_xxxy_xxxy, gs_x_xxxy_xxxz, gs_x_xxxy_xxyy, gs_x_xxxy_xxyz, gs_x_xxxy_xxzz, gs_x_xxxy_xyyy, gs_x_xxxy_xyyz, gs_x_xxxy_xyzz, gs_x_xxxy_xzzz, gs_x_xxxy_yyyy, gs_x_xxxy_yyyz, gs_x_xxxy_yyzz, gs_x_xxxy_yzzz, gs_x_xxxy_zzzz, ts_xxxy_xxx, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxy, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxz, ts_xxxy_xxzz, ts_xxxy_xyy, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyz, ts_xxxy_xyzz, ts_xxxy_xzz, ts_xxxy_xzzz, ts_xxxy_yyy, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyz, ts_xxxy_yyzz, ts_xxxy_yzz, ts_xxxy_yzzz, ts_xxxy_zzz, ts_xxxy_zzzz, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxzz, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyzz, ts_xxy_xzzz, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyzz, ts_xxy_yzzz, ts_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxy_xxxx[i] = 6.0 * ts_xxy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxy[i] = 6.0 * ts_xxy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxxz[i] = 6.0 * ts_xxy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxyy[i] = 6.0 * ts_xxy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxyz[i] = 6.0 * ts_xxy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xxzz[i] = 6.0 * ts_xxy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyyy[i] = 6.0 * ts_xxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyyz[i] = 6.0 * ts_xxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xyzz[i] = 6.0 * ts_xxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xzzz[i] = 6.0 * ts_xxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyyy[i] = 6.0 * ts_xxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyyz[i] = 6.0 * ts_xxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yyzz[i] = 6.0 * ts_xxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yzzz[i] = 6.0 * ts_xxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_zzzz[i] = 6.0 * ts_xxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-45 components of targeted buffer : GG

    auto gs_x_xxxz_xxxx = pbuffer.data(idx_g_gg + 30);

    auto gs_x_xxxz_xxxy = pbuffer.data(idx_g_gg + 31);

    auto gs_x_xxxz_xxxz = pbuffer.data(idx_g_gg + 32);

    auto gs_x_xxxz_xxyy = pbuffer.data(idx_g_gg + 33);

    auto gs_x_xxxz_xxyz = pbuffer.data(idx_g_gg + 34);

    auto gs_x_xxxz_xxzz = pbuffer.data(idx_g_gg + 35);

    auto gs_x_xxxz_xyyy = pbuffer.data(idx_g_gg + 36);

    auto gs_x_xxxz_xyyz = pbuffer.data(idx_g_gg + 37);

    auto gs_x_xxxz_xyzz = pbuffer.data(idx_g_gg + 38);

    auto gs_x_xxxz_xzzz = pbuffer.data(idx_g_gg + 39);

    auto gs_x_xxxz_yyyy = pbuffer.data(idx_g_gg + 40);

    auto gs_x_xxxz_yyyz = pbuffer.data(idx_g_gg + 41);

    auto gs_x_xxxz_yyzz = pbuffer.data(idx_g_gg + 42);

    auto gs_x_xxxz_yzzz = pbuffer.data(idx_g_gg + 43);

    auto gs_x_xxxz_zzzz = pbuffer.data(idx_g_gg + 44);

    #pragma omp simd aligned(gc_x, gs_x_xxxz_xxxx, gs_x_xxxz_xxxy, gs_x_xxxz_xxxz, gs_x_xxxz_xxyy, gs_x_xxxz_xxyz, gs_x_xxxz_xxzz, gs_x_xxxz_xyyy, gs_x_xxxz_xyyz, gs_x_xxxz_xyzz, gs_x_xxxz_xzzz, gs_x_xxxz_yyyy, gs_x_xxxz_yyyz, gs_x_xxxz_yyzz, gs_x_xxxz_yzzz, gs_x_xxxz_zzzz, ts_xxxz_xxx, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxy, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxz, ts_xxxz_xxzz, ts_xxxz_xyy, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyz, ts_xxxz_xyzz, ts_xxxz_xzz, ts_xxxz_xzzz, ts_xxxz_yyy, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyz, ts_xxxz_yyzz, ts_xxxz_yzz, ts_xxxz_yzzz, ts_xxxz_zzz, ts_xxxz_zzzz, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxzz, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyzz, ts_xxz_xzzz, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyzz, ts_xxz_yzzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxz_xxxx[i] = 6.0 * ts_xxz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxy[i] = 6.0 * ts_xxz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxxz[i] = 6.0 * ts_xxz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxyy[i] = 6.0 * ts_xxz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxyz[i] = 6.0 * ts_xxz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xxzz[i] = 6.0 * ts_xxz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyyy[i] = 6.0 * ts_xxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyyz[i] = 6.0 * ts_xxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xyzz[i] = 6.0 * ts_xxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xzzz[i] = 6.0 * ts_xxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyyy[i] = 6.0 * ts_xxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyyz[i] = 6.0 * ts_xxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yyzz[i] = 6.0 * ts_xxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yzzz[i] = 6.0 * ts_xxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_zzzz[i] = 6.0 * ts_xxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 45-60 components of targeted buffer : GG

    auto gs_x_xxyy_xxxx = pbuffer.data(idx_g_gg + 45);

    auto gs_x_xxyy_xxxy = pbuffer.data(idx_g_gg + 46);

    auto gs_x_xxyy_xxxz = pbuffer.data(idx_g_gg + 47);

    auto gs_x_xxyy_xxyy = pbuffer.data(idx_g_gg + 48);

    auto gs_x_xxyy_xxyz = pbuffer.data(idx_g_gg + 49);

    auto gs_x_xxyy_xxzz = pbuffer.data(idx_g_gg + 50);

    auto gs_x_xxyy_xyyy = pbuffer.data(idx_g_gg + 51);

    auto gs_x_xxyy_xyyz = pbuffer.data(idx_g_gg + 52);

    auto gs_x_xxyy_xyzz = pbuffer.data(idx_g_gg + 53);

    auto gs_x_xxyy_xzzz = pbuffer.data(idx_g_gg + 54);

    auto gs_x_xxyy_yyyy = pbuffer.data(idx_g_gg + 55);

    auto gs_x_xxyy_yyyz = pbuffer.data(idx_g_gg + 56);

    auto gs_x_xxyy_yyzz = pbuffer.data(idx_g_gg + 57);

    auto gs_x_xxyy_yzzz = pbuffer.data(idx_g_gg + 58);

    auto gs_x_xxyy_zzzz = pbuffer.data(idx_g_gg + 59);

    #pragma omp simd aligned(gc_x, gs_x_xxyy_xxxx, gs_x_xxyy_xxxy, gs_x_xxyy_xxxz, gs_x_xxyy_xxyy, gs_x_xxyy_xxyz, gs_x_xxyy_xxzz, gs_x_xxyy_xyyy, gs_x_xxyy_xyyz, gs_x_xxyy_xyzz, gs_x_xxyy_xzzz, gs_x_xxyy_yyyy, gs_x_xxyy_yyyz, gs_x_xxyy_yyzz, gs_x_xxyy_yzzz, gs_x_xxyy_zzzz, ts_xxyy_xxx, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxz, ts_xxyy_xxzz, ts_xxyy_xyy, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyz, ts_xxyy_xyzz, ts_xxyy_xzz, ts_xxyy_xzzz, ts_xxyy_yyy, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyz, ts_xxyy_yyzz, ts_xxyy_yzz, ts_xxyy_yzzz, ts_xxyy_zzz, ts_xxyy_zzzz, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxzz, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyzz, ts_xyy_xzzz, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyzz, ts_xyy_yzzz, ts_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyy_xxxx[i] = 4.0 * ts_xyy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxy[i] = 4.0 * ts_xyy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxxz[i] = 4.0 * ts_xyy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxyy[i] = 4.0 * ts_xyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxyz[i] = 4.0 * ts_xyy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xxzz[i] = 4.0 * ts_xyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyyy[i] = 4.0 * ts_xyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyyz[i] = 4.0 * ts_xyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xyzz[i] = 4.0 * ts_xyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xzzz[i] = 4.0 * ts_xyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyyy[i] = 4.0 * ts_xyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyyz[i] = 4.0 * ts_xyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yyzz[i] = 4.0 * ts_xyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yzzz[i] = 4.0 * ts_xyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_zzzz[i] = 4.0 * ts_xyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-75 components of targeted buffer : GG

    auto gs_x_xxyz_xxxx = pbuffer.data(idx_g_gg + 60);

    auto gs_x_xxyz_xxxy = pbuffer.data(idx_g_gg + 61);

    auto gs_x_xxyz_xxxz = pbuffer.data(idx_g_gg + 62);

    auto gs_x_xxyz_xxyy = pbuffer.data(idx_g_gg + 63);

    auto gs_x_xxyz_xxyz = pbuffer.data(idx_g_gg + 64);

    auto gs_x_xxyz_xxzz = pbuffer.data(idx_g_gg + 65);

    auto gs_x_xxyz_xyyy = pbuffer.data(idx_g_gg + 66);

    auto gs_x_xxyz_xyyz = pbuffer.data(idx_g_gg + 67);

    auto gs_x_xxyz_xyzz = pbuffer.data(idx_g_gg + 68);

    auto gs_x_xxyz_xzzz = pbuffer.data(idx_g_gg + 69);

    auto gs_x_xxyz_yyyy = pbuffer.data(idx_g_gg + 70);

    auto gs_x_xxyz_yyyz = pbuffer.data(idx_g_gg + 71);

    auto gs_x_xxyz_yyzz = pbuffer.data(idx_g_gg + 72);

    auto gs_x_xxyz_yzzz = pbuffer.data(idx_g_gg + 73);

    auto gs_x_xxyz_zzzz = pbuffer.data(idx_g_gg + 74);

    #pragma omp simd aligned(gc_x, gs_x_xxyz_xxxx, gs_x_xxyz_xxxy, gs_x_xxyz_xxxz, gs_x_xxyz_xxyy, gs_x_xxyz_xxyz, gs_x_xxyz_xxzz, gs_x_xxyz_xyyy, gs_x_xxyz_xyyz, gs_x_xxyz_xyzz, gs_x_xxyz_xzzz, gs_x_xxyz_yyyy, gs_x_xxyz_yyyz, gs_x_xxyz_yyzz, gs_x_xxyz_yzzz, gs_x_xxyz_zzzz, ts_xxyz_xxx, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxy, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxz, ts_xxyz_xxzz, ts_xxyz_xyy, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyz, ts_xxyz_xyzz, ts_xxyz_xzz, ts_xxyz_xzzz, ts_xxyz_yyy, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyz, ts_xxyz_yyzz, ts_xxyz_yzz, ts_xxyz_yzzz, ts_xxyz_zzz, ts_xxyz_zzzz, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxzz, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyzz, ts_xyz_xzzz, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyzz, ts_xyz_yzzz, ts_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyz_xxxx[i] = 4.0 * ts_xyz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxy[i] = 4.0 * ts_xyz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxxz[i] = 4.0 * ts_xyz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxyy[i] = 4.0 * ts_xyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxyz[i] = 4.0 * ts_xyz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xxzz[i] = 4.0 * ts_xyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyyy[i] = 4.0 * ts_xyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyyz[i] = 4.0 * ts_xyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xyzz[i] = 4.0 * ts_xyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xzzz[i] = 4.0 * ts_xyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyyy[i] = 4.0 * ts_xyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyyz[i] = 4.0 * ts_xyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yyzz[i] = 4.0 * ts_xyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yzzz[i] = 4.0 * ts_xyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_zzzz[i] = 4.0 * ts_xyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 75-90 components of targeted buffer : GG

    auto gs_x_xxzz_xxxx = pbuffer.data(idx_g_gg + 75);

    auto gs_x_xxzz_xxxy = pbuffer.data(idx_g_gg + 76);

    auto gs_x_xxzz_xxxz = pbuffer.data(idx_g_gg + 77);

    auto gs_x_xxzz_xxyy = pbuffer.data(idx_g_gg + 78);

    auto gs_x_xxzz_xxyz = pbuffer.data(idx_g_gg + 79);

    auto gs_x_xxzz_xxzz = pbuffer.data(idx_g_gg + 80);

    auto gs_x_xxzz_xyyy = pbuffer.data(idx_g_gg + 81);

    auto gs_x_xxzz_xyyz = pbuffer.data(idx_g_gg + 82);

    auto gs_x_xxzz_xyzz = pbuffer.data(idx_g_gg + 83);

    auto gs_x_xxzz_xzzz = pbuffer.data(idx_g_gg + 84);

    auto gs_x_xxzz_yyyy = pbuffer.data(idx_g_gg + 85);

    auto gs_x_xxzz_yyyz = pbuffer.data(idx_g_gg + 86);

    auto gs_x_xxzz_yyzz = pbuffer.data(idx_g_gg + 87);

    auto gs_x_xxzz_yzzz = pbuffer.data(idx_g_gg + 88);

    auto gs_x_xxzz_zzzz = pbuffer.data(idx_g_gg + 89);

    #pragma omp simd aligned(gc_x, gs_x_xxzz_xxxx, gs_x_xxzz_xxxy, gs_x_xxzz_xxxz, gs_x_xxzz_xxyy, gs_x_xxzz_xxyz, gs_x_xxzz_xxzz, gs_x_xxzz_xyyy, gs_x_xxzz_xyyz, gs_x_xxzz_xyzz, gs_x_xxzz_xzzz, gs_x_xxzz_yyyy, gs_x_xxzz_yyyz, gs_x_xxzz_yyzz, gs_x_xxzz_yzzz, gs_x_xxzz_zzzz, ts_xxzz_xxx, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxy, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxz, ts_xxzz_xxzz, ts_xxzz_xyy, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyz, ts_xxzz_xyzz, ts_xxzz_xzz, ts_xxzz_xzzz, ts_xxzz_yyy, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyz, ts_xxzz_yyzz, ts_xxzz_yzz, ts_xxzz_yzzz, ts_xxzz_zzz, ts_xxzz_zzzz, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxzz, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyzz, ts_xzz_xzzz, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyzz, ts_xzz_yzzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzz_xxxx[i] = 4.0 * ts_xzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxy[i] = 4.0 * ts_xzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxxz[i] = 4.0 * ts_xzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxyy[i] = 4.0 * ts_xzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxyz[i] = 4.0 * ts_xzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xxzz[i] = 4.0 * ts_xzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyyy[i] = 4.0 * ts_xzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyyz[i] = 4.0 * ts_xzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xyzz[i] = 4.0 * ts_xzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xzzz[i] = 4.0 * ts_xzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyyy[i] = 4.0 * ts_xzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyyz[i] = 4.0 * ts_xzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yyzz[i] = 4.0 * ts_xzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yzzz[i] = 4.0 * ts_xzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_zzzz[i] = 4.0 * ts_xzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-105 components of targeted buffer : GG

    auto gs_x_xyyy_xxxx = pbuffer.data(idx_g_gg + 90);

    auto gs_x_xyyy_xxxy = pbuffer.data(idx_g_gg + 91);

    auto gs_x_xyyy_xxxz = pbuffer.data(idx_g_gg + 92);

    auto gs_x_xyyy_xxyy = pbuffer.data(idx_g_gg + 93);

    auto gs_x_xyyy_xxyz = pbuffer.data(idx_g_gg + 94);

    auto gs_x_xyyy_xxzz = pbuffer.data(idx_g_gg + 95);

    auto gs_x_xyyy_xyyy = pbuffer.data(idx_g_gg + 96);

    auto gs_x_xyyy_xyyz = pbuffer.data(idx_g_gg + 97);

    auto gs_x_xyyy_xyzz = pbuffer.data(idx_g_gg + 98);

    auto gs_x_xyyy_xzzz = pbuffer.data(idx_g_gg + 99);

    auto gs_x_xyyy_yyyy = pbuffer.data(idx_g_gg + 100);

    auto gs_x_xyyy_yyyz = pbuffer.data(idx_g_gg + 101);

    auto gs_x_xyyy_yyzz = pbuffer.data(idx_g_gg + 102);

    auto gs_x_xyyy_yzzz = pbuffer.data(idx_g_gg + 103);

    auto gs_x_xyyy_zzzz = pbuffer.data(idx_g_gg + 104);

    #pragma omp simd aligned(gc_x, gs_x_xyyy_xxxx, gs_x_xyyy_xxxy, gs_x_xyyy_xxxz, gs_x_xyyy_xxyy, gs_x_xyyy_xxyz, gs_x_xyyy_xxzz, gs_x_xyyy_xyyy, gs_x_xyyy_xyyz, gs_x_xyyy_xyzz, gs_x_xyyy_xzzz, gs_x_xyyy_yyyy, gs_x_xyyy_yyyz, gs_x_xyyy_yyzz, gs_x_xyyy_yzzz, gs_x_xyyy_zzzz, ts_xyyy_xxx, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxy, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxz, ts_xyyy_xxzz, ts_xyyy_xyy, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyz, ts_xyyy_xyzz, ts_xyyy_xzz, ts_xyyy_xzzz, ts_xyyy_yyy, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyz, ts_xyyy_yyzz, ts_xyyy_yzz, ts_xyyy_yzzz, ts_xyyy_zzz, ts_xyyy_zzzz, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxzz, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyzz, ts_yyy_xzzz, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyzz, ts_yyy_yzzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyy_xxxx[i] = 2.0 * ts_yyy_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxy[i] = 2.0 * ts_yyy_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxxz[i] = 2.0 * ts_yyy_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxyy[i] = 2.0 * ts_yyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxyz[i] = 2.0 * ts_yyy_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xxzz[i] = 2.0 * ts_yyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyyy[i] = 2.0 * ts_yyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyyz[i] = 2.0 * ts_yyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xyzz[i] = 2.0 * ts_yyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xzzz[i] = 2.0 * ts_yyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyyy[i] = 2.0 * ts_yyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyyz[i] = 2.0 * ts_yyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yyzz[i] = 2.0 * ts_yyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yzzz[i] = 2.0 * ts_yyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_zzzz[i] = 2.0 * ts_yyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-120 components of targeted buffer : GG

    auto gs_x_xyyz_xxxx = pbuffer.data(idx_g_gg + 105);

    auto gs_x_xyyz_xxxy = pbuffer.data(idx_g_gg + 106);

    auto gs_x_xyyz_xxxz = pbuffer.data(idx_g_gg + 107);

    auto gs_x_xyyz_xxyy = pbuffer.data(idx_g_gg + 108);

    auto gs_x_xyyz_xxyz = pbuffer.data(idx_g_gg + 109);

    auto gs_x_xyyz_xxzz = pbuffer.data(idx_g_gg + 110);

    auto gs_x_xyyz_xyyy = pbuffer.data(idx_g_gg + 111);

    auto gs_x_xyyz_xyyz = pbuffer.data(idx_g_gg + 112);

    auto gs_x_xyyz_xyzz = pbuffer.data(idx_g_gg + 113);

    auto gs_x_xyyz_xzzz = pbuffer.data(idx_g_gg + 114);

    auto gs_x_xyyz_yyyy = pbuffer.data(idx_g_gg + 115);

    auto gs_x_xyyz_yyyz = pbuffer.data(idx_g_gg + 116);

    auto gs_x_xyyz_yyzz = pbuffer.data(idx_g_gg + 117);

    auto gs_x_xyyz_yzzz = pbuffer.data(idx_g_gg + 118);

    auto gs_x_xyyz_zzzz = pbuffer.data(idx_g_gg + 119);

    #pragma omp simd aligned(gc_x, gs_x_xyyz_xxxx, gs_x_xyyz_xxxy, gs_x_xyyz_xxxz, gs_x_xyyz_xxyy, gs_x_xyyz_xxyz, gs_x_xyyz_xxzz, gs_x_xyyz_xyyy, gs_x_xyyz_xyyz, gs_x_xyyz_xyzz, gs_x_xyyz_xzzz, gs_x_xyyz_yyyy, gs_x_xyyz_yyyz, gs_x_xyyz_yyzz, gs_x_xyyz_yzzz, gs_x_xyyz_zzzz, ts_xyyz_xxx, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxy, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxz, ts_xyyz_xxzz, ts_xyyz_xyy, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyz, ts_xyyz_xyzz, ts_xyyz_xzz, ts_xyyz_xzzz, ts_xyyz_yyy, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyz, ts_xyyz_yyzz, ts_xyyz_yzz, ts_xyyz_yzzz, ts_xyyz_zzz, ts_xyyz_zzzz, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxzz, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyzz, ts_yyz_xzzz, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyzz, ts_yyz_yzzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyz_xxxx[i] = 2.0 * ts_yyz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxy[i] = 2.0 * ts_yyz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxxz[i] = 2.0 * ts_yyz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxyy[i] = 2.0 * ts_yyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxyz[i] = 2.0 * ts_yyz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xxzz[i] = 2.0 * ts_yyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyyy[i] = 2.0 * ts_yyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyyz[i] = 2.0 * ts_yyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xyzz[i] = 2.0 * ts_yyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xzzz[i] = 2.0 * ts_yyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyyy[i] = 2.0 * ts_yyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyyz[i] = 2.0 * ts_yyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yyzz[i] = 2.0 * ts_yyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yzzz[i] = 2.0 * ts_yyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_zzzz[i] = 2.0 * ts_yyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 120-135 components of targeted buffer : GG

    auto gs_x_xyzz_xxxx = pbuffer.data(idx_g_gg + 120);

    auto gs_x_xyzz_xxxy = pbuffer.data(idx_g_gg + 121);

    auto gs_x_xyzz_xxxz = pbuffer.data(idx_g_gg + 122);

    auto gs_x_xyzz_xxyy = pbuffer.data(idx_g_gg + 123);

    auto gs_x_xyzz_xxyz = pbuffer.data(idx_g_gg + 124);

    auto gs_x_xyzz_xxzz = pbuffer.data(idx_g_gg + 125);

    auto gs_x_xyzz_xyyy = pbuffer.data(idx_g_gg + 126);

    auto gs_x_xyzz_xyyz = pbuffer.data(idx_g_gg + 127);

    auto gs_x_xyzz_xyzz = pbuffer.data(idx_g_gg + 128);

    auto gs_x_xyzz_xzzz = pbuffer.data(idx_g_gg + 129);

    auto gs_x_xyzz_yyyy = pbuffer.data(idx_g_gg + 130);

    auto gs_x_xyzz_yyyz = pbuffer.data(idx_g_gg + 131);

    auto gs_x_xyzz_yyzz = pbuffer.data(idx_g_gg + 132);

    auto gs_x_xyzz_yzzz = pbuffer.data(idx_g_gg + 133);

    auto gs_x_xyzz_zzzz = pbuffer.data(idx_g_gg + 134);

    #pragma omp simd aligned(gc_x, gs_x_xyzz_xxxx, gs_x_xyzz_xxxy, gs_x_xyzz_xxxz, gs_x_xyzz_xxyy, gs_x_xyzz_xxyz, gs_x_xyzz_xxzz, gs_x_xyzz_xyyy, gs_x_xyzz_xyyz, gs_x_xyzz_xyzz, gs_x_xyzz_xzzz, gs_x_xyzz_yyyy, gs_x_xyzz_yyyz, gs_x_xyzz_yyzz, gs_x_xyzz_yzzz, gs_x_xyzz_zzzz, ts_xyzz_xxx, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxy, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxz, ts_xyzz_xxzz, ts_xyzz_xyy, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyz, ts_xyzz_xyzz, ts_xyzz_xzz, ts_xyzz_xzzz, ts_xyzz_yyy, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyz, ts_xyzz_yyzz, ts_xyzz_yzz, ts_xyzz_yzzz, ts_xyzz_zzz, ts_xyzz_zzzz, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxzz, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyzz, ts_yzz_xzzz, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyzz, ts_yzz_yzzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzz_xxxx[i] = 2.0 * ts_yzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxy[i] = 2.0 * ts_yzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxxz[i] = 2.0 * ts_yzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxyy[i] = 2.0 * ts_yzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxyz[i] = 2.0 * ts_yzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xxzz[i] = 2.0 * ts_yzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyyy[i] = 2.0 * ts_yzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyyz[i] = 2.0 * ts_yzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xyzz[i] = 2.0 * ts_yzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xzzz[i] = 2.0 * ts_yzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyyy[i] = 2.0 * ts_yzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyyz[i] = 2.0 * ts_yzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yyzz[i] = 2.0 * ts_yzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yzzz[i] = 2.0 * ts_yzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_zzzz[i] = 2.0 * ts_yzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 135-150 components of targeted buffer : GG

    auto gs_x_xzzz_xxxx = pbuffer.data(idx_g_gg + 135);

    auto gs_x_xzzz_xxxy = pbuffer.data(idx_g_gg + 136);

    auto gs_x_xzzz_xxxz = pbuffer.data(idx_g_gg + 137);

    auto gs_x_xzzz_xxyy = pbuffer.data(idx_g_gg + 138);

    auto gs_x_xzzz_xxyz = pbuffer.data(idx_g_gg + 139);

    auto gs_x_xzzz_xxzz = pbuffer.data(idx_g_gg + 140);

    auto gs_x_xzzz_xyyy = pbuffer.data(idx_g_gg + 141);

    auto gs_x_xzzz_xyyz = pbuffer.data(idx_g_gg + 142);

    auto gs_x_xzzz_xyzz = pbuffer.data(idx_g_gg + 143);

    auto gs_x_xzzz_xzzz = pbuffer.data(idx_g_gg + 144);

    auto gs_x_xzzz_yyyy = pbuffer.data(idx_g_gg + 145);

    auto gs_x_xzzz_yyyz = pbuffer.data(idx_g_gg + 146);

    auto gs_x_xzzz_yyzz = pbuffer.data(idx_g_gg + 147);

    auto gs_x_xzzz_yzzz = pbuffer.data(idx_g_gg + 148);

    auto gs_x_xzzz_zzzz = pbuffer.data(idx_g_gg + 149);

    #pragma omp simd aligned(gc_x, gs_x_xzzz_xxxx, gs_x_xzzz_xxxy, gs_x_xzzz_xxxz, gs_x_xzzz_xxyy, gs_x_xzzz_xxyz, gs_x_xzzz_xxzz, gs_x_xzzz_xyyy, gs_x_xzzz_xyyz, gs_x_xzzz_xyzz, gs_x_xzzz_xzzz, gs_x_xzzz_yyyy, gs_x_xzzz_yyyz, gs_x_xzzz_yyzz, gs_x_xzzz_yzzz, gs_x_xzzz_zzzz, ts_xzzz_xxx, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxy, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxz, ts_xzzz_xxzz, ts_xzzz_xyy, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyz, ts_xzzz_xyzz, ts_xzzz_xzz, ts_xzzz_xzzz, ts_xzzz_yyy, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyz, ts_xzzz_yyzz, ts_xzzz_yzz, ts_xzzz_yzzz, ts_xzzz_zzz, ts_xzzz_zzzz, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxzz, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyzz, ts_zzz_xzzz, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyzz, ts_zzz_yzzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzz_xxxx[i] = 2.0 * ts_zzz_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxy[i] = 2.0 * ts_zzz_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxxz[i] = 2.0 * ts_zzz_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxyy[i] = 2.0 * ts_zzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxyz[i] = 2.0 * ts_zzz_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xxzz[i] = 2.0 * ts_zzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyyy[i] = 2.0 * ts_zzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyyz[i] = 2.0 * ts_zzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xyzz[i] = 2.0 * ts_zzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xzzz[i] = 2.0 * ts_zzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyyy[i] = 2.0 * ts_zzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyyz[i] = 2.0 * ts_zzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yyzz[i] = 2.0 * ts_zzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yzzz[i] = 2.0 * ts_zzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_zzzz[i] = 2.0 * ts_zzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 150-165 components of targeted buffer : GG

    auto gs_x_yyyy_xxxx = pbuffer.data(idx_g_gg + 150);

    auto gs_x_yyyy_xxxy = pbuffer.data(idx_g_gg + 151);

    auto gs_x_yyyy_xxxz = pbuffer.data(idx_g_gg + 152);

    auto gs_x_yyyy_xxyy = pbuffer.data(idx_g_gg + 153);

    auto gs_x_yyyy_xxyz = pbuffer.data(idx_g_gg + 154);

    auto gs_x_yyyy_xxzz = pbuffer.data(idx_g_gg + 155);

    auto gs_x_yyyy_xyyy = pbuffer.data(idx_g_gg + 156);

    auto gs_x_yyyy_xyyz = pbuffer.data(idx_g_gg + 157);

    auto gs_x_yyyy_xyzz = pbuffer.data(idx_g_gg + 158);

    auto gs_x_yyyy_xzzz = pbuffer.data(idx_g_gg + 159);

    auto gs_x_yyyy_yyyy = pbuffer.data(idx_g_gg + 160);

    auto gs_x_yyyy_yyyz = pbuffer.data(idx_g_gg + 161);

    auto gs_x_yyyy_yyzz = pbuffer.data(idx_g_gg + 162);

    auto gs_x_yyyy_yzzz = pbuffer.data(idx_g_gg + 163);

    auto gs_x_yyyy_zzzz = pbuffer.data(idx_g_gg + 164);

    #pragma omp simd aligned(gc_x, gs_x_yyyy_xxxx, gs_x_yyyy_xxxy, gs_x_yyyy_xxxz, gs_x_yyyy_xxyy, gs_x_yyyy_xxyz, gs_x_yyyy_xxzz, gs_x_yyyy_xyyy, gs_x_yyyy_xyyz, gs_x_yyyy_xyzz, gs_x_yyyy_xzzz, gs_x_yyyy_yyyy, gs_x_yyyy_yyyz, gs_x_yyyy_yyzz, gs_x_yyyy_yzzz, gs_x_yyyy_zzzz, ts_yyyy_xxx, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxz, ts_yyyy_xxzz, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_xzz, ts_yyyy_xzzz, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyy_zzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyy_xxxx[i] = 8.0 * ts_yyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxy[i] = 6.0 * ts_yyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxxz[i] = 6.0 * ts_yyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxyy[i] = 4.0 * ts_yyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxyz[i] = 4.0 * ts_yyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xxzz[i] = 4.0 * ts_yyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyyy[i] = 2.0 * ts_yyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyyz[i] = 2.0 * ts_yyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xyzz[i] = 2.0 * ts_yyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xzzz[i] = 2.0 * ts_yyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyyy[i] = 2.0 * ts_yyyy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyyz[i] = 2.0 * ts_yyyy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yyzz[i] = 2.0 * ts_yyyy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yzzz[i] = 2.0 * ts_yyyy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_zzzz[i] = 2.0 * ts_yyyy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 165-180 components of targeted buffer : GG

    auto gs_x_yyyz_xxxx = pbuffer.data(idx_g_gg + 165);

    auto gs_x_yyyz_xxxy = pbuffer.data(idx_g_gg + 166);

    auto gs_x_yyyz_xxxz = pbuffer.data(idx_g_gg + 167);

    auto gs_x_yyyz_xxyy = pbuffer.data(idx_g_gg + 168);

    auto gs_x_yyyz_xxyz = pbuffer.data(idx_g_gg + 169);

    auto gs_x_yyyz_xxzz = pbuffer.data(idx_g_gg + 170);

    auto gs_x_yyyz_xyyy = pbuffer.data(idx_g_gg + 171);

    auto gs_x_yyyz_xyyz = pbuffer.data(idx_g_gg + 172);

    auto gs_x_yyyz_xyzz = pbuffer.data(idx_g_gg + 173);

    auto gs_x_yyyz_xzzz = pbuffer.data(idx_g_gg + 174);

    auto gs_x_yyyz_yyyy = pbuffer.data(idx_g_gg + 175);

    auto gs_x_yyyz_yyyz = pbuffer.data(idx_g_gg + 176);

    auto gs_x_yyyz_yyzz = pbuffer.data(idx_g_gg + 177);

    auto gs_x_yyyz_yzzz = pbuffer.data(idx_g_gg + 178);

    auto gs_x_yyyz_zzzz = pbuffer.data(idx_g_gg + 179);

    #pragma omp simd aligned(gc_x, gs_x_yyyz_xxxx, gs_x_yyyz_xxxy, gs_x_yyyz_xxxz, gs_x_yyyz_xxyy, gs_x_yyyz_xxyz, gs_x_yyyz_xxzz, gs_x_yyyz_xyyy, gs_x_yyyz_xyyz, gs_x_yyyz_xyzz, gs_x_yyyz_xzzz, gs_x_yyyz_yyyy, gs_x_yyyz_yyyz, gs_x_yyyz_yyzz, gs_x_yyyz_yzzz, gs_x_yyyz_zzzz, ts_yyyz_xxx, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxy, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxz, ts_yyyz_xxzz, ts_yyyz_xyy, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyz, ts_yyyz_xyzz, ts_yyyz_xzz, ts_yyyz_xzzz, ts_yyyz_yyy, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyz, ts_yyyz_yyzz, ts_yyyz_yzz, ts_yyyz_yzzz, ts_yyyz_zzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyz_xxxx[i] = 8.0 * ts_yyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxy[i] = 6.0 * ts_yyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxxz[i] = 6.0 * ts_yyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxyy[i] = 4.0 * ts_yyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxyz[i] = 4.0 * ts_yyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xxzz[i] = 4.0 * ts_yyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyyy[i] = 2.0 * ts_yyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyyz[i] = 2.0 * ts_yyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xyzz[i] = 2.0 * ts_yyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xzzz[i] = 2.0 * ts_yyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyyy[i] = 2.0 * ts_yyyz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyyz[i] = 2.0 * ts_yyyz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yyzz[i] = 2.0 * ts_yyyz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yzzz[i] = 2.0 * ts_yyyz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_zzzz[i] = 2.0 * ts_yyyz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 180-195 components of targeted buffer : GG

    auto gs_x_yyzz_xxxx = pbuffer.data(idx_g_gg + 180);

    auto gs_x_yyzz_xxxy = pbuffer.data(idx_g_gg + 181);

    auto gs_x_yyzz_xxxz = pbuffer.data(idx_g_gg + 182);

    auto gs_x_yyzz_xxyy = pbuffer.data(idx_g_gg + 183);

    auto gs_x_yyzz_xxyz = pbuffer.data(idx_g_gg + 184);

    auto gs_x_yyzz_xxzz = pbuffer.data(idx_g_gg + 185);

    auto gs_x_yyzz_xyyy = pbuffer.data(idx_g_gg + 186);

    auto gs_x_yyzz_xyyz = pbuffer.data(idx_g_gg + 187);

    auto gs_x_yyzz_xyzz = pbuffer.data(idx_g_gg + 188);

    auto gs_x_yyzz_xzzz = pbuffer.data(idx_g_gg + 189);

    auto gs_x_yyzz_yyyy = pbuffer.data(idx_g_gg + 190);

    auto gs_x_yyzz_yyyz = pbuffer.data(idx_g_gg + 191);

    auto gs_x_yyzz_yyzz = pbuffer.data(idx_g_gg + 192);

    auto gs_x_yyzz_yzzz = pbuffer.data(idx_g_gg + 193);

    auto gs_x_yyzz_zzzz = pbuffer.data(idx_g_gg + 194);

    #pragma omp simd aligned(gc_x, gs_x_yyzz_xxxx, gs_x_yyzz_xxxy, gs_x_yyzz_xxxz, gs_x_yyzz_xxyy, gs_x_yyzz_xxyz, gs_x_yyzz_xxzz, gs_x_yyzz_xyyy, gs_x_yyzz_xyyz, gs_x_yyzz_xyzz, gs_x_yyzz_xzzz, gs_x_yyzz_yyyy, gs_x_yyzz_yyyz, gs_x_yyzz_yyzz, gs_x_yyzz_yzzz, gs_x_yyzz_zzzz, ts_yyzz_xxx, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxy, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxz, ts_yyzz_xxzz, ts_yyzz_xyy, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyz, ts_yyzz_xyzz, ts_yyzz_xzz, ts_yyzz_xzzz, ts_yyzz_yyy, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyz, ts_yyzz_yyzz, ts_yyzz_yzz, ts_yyzz_yzzz, ts_yyzz_zzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzz_xxxx[i] = 8.0 * ts_yyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxy[i] = 6.0 * ts_yyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxxz[i] = 6.0 * ts_yyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxyy[i] = 4.0 * ts_yyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxyz[i] = 4.0 * ts_yyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xxzz[i] = 4.0 * ts_yyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyyy[i] = 2.0 * ts_yyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyyz[i] = 2.0 * ts_yyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xyzz[i] = 2.0 * ts_yyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xzzz[i] = 2.0 * ts_yyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyyy[i] = 2.0 * ts_yyzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyyz[i] = 2.0 * ts_yyzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yyzz[i] = 2.0 * ts_yyzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yzzz[i] = 2.0 * ts_yyzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_zzzz[i] = 2.0 * ts_yyzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 195-210 components of targeted buffer : GG

    auto gs_x_yzzz_xxxx = pbuffer.data(idx_g_gg + 195);

    auto gs_x_yzzz_xxxy = pbuffer.data(idx_g_gg + 196);

    auto gs_x_yzzz_xxxz = pbuffer.data(idx_g_gg + 197);

    auto gs_x_yzzz_xxyy = pbuffer.data(idx_g_gg + 198);

    auto gs_x_yzzz_xxyz = pbuffer.data(idx_g_gg + 199);

    auto gs_x_yzzz_xxzz = pbuffer.data(idx_g_gg + 200);

    auto gs_x_yzzz_xyyy = pbuffer.data(idx_g_gg + 201);

    auto gs_x_yzzz_xyyz = pbuffer.data(idx_g_gg + 202);

    auto gs_x_yzzz_xyzz = pbuffer.data(idx_g_gg + 203);

    auto gs_x_yzzz_xzzz = pbuffer.data(idx_g_gg + 204);

    auto gs_x_yzzz_yyyy = pbuffer.data(idx_g_gg + 205);

    auto gs_x_yzzz_yyyz = pbuffer.data(idx_g_gg + 206);

    auto gs_x_yzzz_yyzz = pbuffer.data(idx_g_gg + 207);

    auto gs_x_yzzz_yzzz = pbuffer.data(idx_g_gg + 208);

    auto gs_x_yzzz_zzzz = pbuffer.data(idx_g_gg + 209);

    #pragma omp simd aligned(gc_x, gs_x_yzzz_xxxx, gs_x_yzzz_xxxy, gs_x_yzzz_xxxz, gs_x_yzzz_xxyy, gs_x_yzzz_xxyz, gs_x_yzzz_xxzz, gs_x_yzzz_xyyy, gs_x_yzzz_xyyz, gs_x_yzzz_xyzz, gs_x_yzzz_xzzz, gs_x_yzzz_yyyy, gs_x_yzzz_yyyz, gs_x_yzzz_yyzz, gs_x_yzzz_yzzz, gs_x_yzzz_zzzz, ts_yzzz_xxx, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxy, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxz, ts_yzzz_xxzz, ts_yzzz_xyy, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyz, ts_yzzz_xyzz, ts_yzzz_xzz, ts_yzzz_xzzz, ts_yzzz_yyy, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyz, ts_yzzz_yyzz, ts_yzzz_yzz, ts_yzzz_yzzz, ts_yzzz_zzz, ts_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzz_xxxx[i] = 8.0 * ts_yzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxy[i] = 6.0 * ts_yzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxxz[i] = 6.0 * ts_yzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxyy[i] = 4.0 * ts_yzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxyz[i] = 4.0 * ts_yzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xxzz[i] = 4.0 * ts_yzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyyy[i] = 2.0 * ts_yzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyyz[i] = 2.0 * ts_yzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xyzz[i] = 2.0 * ts_yzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xzzz[i] = 2.0 * ts_yzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyyy[i] = 2.0 * ts_yzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyyz[i] = 2.0 * ts_yzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yyzz[i] = 2.0 * ts_yzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yzzz[i] = 2.0 * ts_yzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_zzzz[i] = 2.0 * ts_yzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 210-225 components of targeted buffer : GG

    auto gs_x_zzzz_xxxx = pbuffer.data(idx_g_gg + 210);

    auto gs_x_zzzz_xxxy = pbuffer.data(idx_g_gg + 211);

    auto gs_x_zzzz_xxxz = pbuffer.data(idx_g_gg + 212);

    auto gs_x_zzzz_xxyy = pbuffer.data(idx_g_gg + 213);

    auto gs_x_zzzz_xxyz = pbuffer.data(idx_g_gg + 214);

    auto gs_x_zzzz_xxzz = pbuffer.data(idx_g_gg + 215);

    auto gs_x_zzzz_xyyy = pbuffer.data(idx_g_gg + 216);

    auto gs_x_zzzz_xyyz = pbuffer.data(idx_g_gg + 217);

    auto gs_x_zzzz_xyzz = pbuffer.data(idx_g_gg + 218);

    auto gs_x_zzzz_xzzz = pbuffer.data(idx_g_gg + 219);

    auto gs_x_zzzz_yyyy = pbuffer.data(idx_g_gg + 220);

    auto gs_x_zzzz_yyyz = pbuffer.data(idx_g_gg + 221);

    auto gs_x_zzzz_yyzz = pbuffer.data(idx_g_gg + 222);

    auto gs_x_zzzz_yzzz = pbuffer.data(idx_g_gg + 223);

    auto gs_x_zzzz_zzzz = pbuffer.data(idx_g_gg + 224);

    #pragma omp simd aligned(gc_x, gs_x_zzzz_xxxx, gs_x_zzzz_xxxy, gs_x_zzzz_xxxz, gs_x_zzzz_xxyy, gs_x_zzzz_xxyz, gs_x_zzzz_xxzz, gs_x_zzzz_xyyy, gs_x_zzzz_xyyz, gs_x_zzzz_xyzz, gs_x_zzzz_xzzz, gs_x_zzzz_yyyy, gs_x_zzzz_yyyz, gs_x_zzzz_yyzz, gs_x_zzzz_yzzz, gs_x_zzzz_zzzz, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzz_xxxx[i] = 8.0 * ts_zzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxy[i] = 6.0 * ts_zzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxxz[i] = 6.0 * ts_zzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxyy[i] = 4.0 * ts_zzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxyz[i] = 4.0 * ts_zzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xxzz[i] = 4.0 * ts_zzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyyy[i] = 2.0 * ts_zzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyyz[i] = 2.0 * ts_zzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xyzz[i] = 2.0 * ts_zzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xzzz[i] = 2.0 * ts_zzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyyy[i] = 2.0 * ts_zzzz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyyz[i] = 2.0 * ts_zzzz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yyzz[i] = 2.0 * ts_zzzz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yzzz[i] = 2.0 * ts_zzzz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_zzzz[i] = 2.0 * ts_zzzz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 225-240 components of targeted buffer : GG

    auto gs_y_xxxx_xxxx = pbuffer.data(idx_g_gg + 225);

    auto gs_y_xxxx_xxxy = pbuffer.data(idx_g_gg + 226);

    auto gs_y_xxxx_xxxz = pbuffer.data(idx_g_gg + 227);

    auto gs_y_xxxx_xxyy = pbuffer.data(idx_g_gg + 228);

    auto gs_y_xxxx_xxyz = pbuffer.data(idx_g_gg + 229);

    auto gs_y_xxxx_xxzz = pbuffer.data(idx_g_gg + 230);

    auto gs_y_xxxx_xyyy = pbuffer.data(idx_g_gg + 231);

    auto gs_y_xxxx_xyyz = pbuffer.data(idx_g_gg + 232);

    auto gs_y_xxxx_xyzz = pbuffer.data(idx_g_gg + 233);

    auto gs_y_xxxx_xzzz = pbuffer.data(idx_g_gg + 234);

    auto gs_y_xxxx_yyyy = pbuffer.data(idx_g_gg + 235);

    auto gs_y_xxxx_yyyz = pbuffer.data(idx_g_gg + 236);

    auto gs_y_xxxx_yyzz = pbuffer.data(idx_g_gg + 237);

    auto gs_y_xxxx_yzzz = pbuffer.data(idx_g_gg + 238);

    auto gs_y_xxxx_zzzz = pbuffer.data(idx_g_gg + 239);

    #pragma omp simd aligned(gc_y, gs_y_xxxx_xxxx, gs_y_xxxx_xxxy, gs_y_xxxx_xxxz, gs_y_xxxx_xxyy, gs_y_xxxx_xxyz, gs_y_xxxx_xxzz, gs_y_xxxx_xyyy, gs_y_xxxx_xyyz, gs_y_xxxx_xyzz, gs_y_xxxx_xzzz, gs_y_xxxx_yyyy, gs_y_xxxx_yyyz, gs_y_xxxx_yyzz, gs_y_xxxx_yzzz, gs_y_xxxx_zzzz, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_yyy, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyz, ts_xxxx_yyzz, ts_xxxx_yzz, ts_xxxx_yzzz, ts_xxxx_zzz, ts_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxx_xxxx[i] = 2.0 * ts_xxxx_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxy[i] = 2.0 * ts_xxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxxz[i] = 2.0 * ts_xxxx_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxyy[i] = 4.0 * ts_xxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxyz[i] = 2.0 * ts_xxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xxzz[i] = 2.0 * ts_xxxx_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyyy[i] = 6.0 * ts_xxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyyz[i] = 4.0 * ts_xxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xyzz[i] = 2.0 * ts_xxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xzzz[i] = 2.0 * ts_xxxx_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyyy[i] = 8.0 * ts_xxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyyz[i] = 6.0 * ts_xxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yyzz[i] = 4.0 * ts_xxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yzzz[i] = 2.0 * ts_xxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_zzzz[i] = 2.0 * ts_xxxx_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 240-255 components of targeted buffer : GG

    auto gs_y_xxxy_xxxx = pbuffer.data(idx_g_gg + 240);

    auto gs_y_xxxy_xxxy = pbuffer.data(idx_g_gg + 241);

    auto gs_y_xxxy_xxxz = pbuffer.data(idx_g_gg + 242);

    auto gs_y_xxxy_xxyy = pbuffer.data(idx_g_gg + 243);

    auto gs_y_xxxy_xxyz = pbuffer.data(idx_g_gg + 244);

    auto gs_y_xxxy_xxzz = pbuffer.data(idx_g_gg + 245);

    auto gs_y_xxxy_xyyy = pbuffer.data(idx_g_gg + 246);

    auto gs_y_xxxy_xyyz = pbuffer.data(idx_g_gg + 247);

    auto gs_y_xxxy_xyzz = pbuffer.data(idx_g_gg + 248);

    auto gs_y_xxxy_xzzz = pbuffer.data(idx_g_gg + 249);

    auto gs_y_xxxy_yyyy = pbuffer.data(idx_g_gg + 250);

    auto gs_y_xxxy_yyyz = pbuffer.data(idx_g_gg + 251);

    auto gs_y_xxxy_yyzz = pbuffer.data(idx_g_gg + 252);

    auto gs_y_xxxy_yzzz = pbuffer.data(idx_g_gg + 253);

    auto gs_y_xxxy_zzzz = pbuffer.data(idx_g_gg + 254);

    #pragma omp simd aligned(gc_y, gs_y_xxxy_xxxx, gs_y_xxxy_xxxy, gs_y_xxxy_xxxz, gs_y_xxxy_xxyy, gs_y_xxxy_xxyz, gs_y_xxxy_xxzz, gs_y_xxxy_xyyy, gs_y_xxxy_xyyz, gs_y_xxxy_xyzz, gs_y_xxxy_xzzz, gs_y_xxxy_yyyy, gs_y_xxxy_yyyz, gs_y_xxxy_yyzz, gs_y_xxxy_yzzz, gs_y_xxxy_zzzz, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxzz, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyzz, ts_xxx_xzzz, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyzz, ts_xxx_yzzz, ts_xxx_zzzz, ts_xxxy_xxx, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxy, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxz, ts_xxxy_xxzz, ts_xxxy_xyy, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyz, ts_xxxy_xyzz, ts_xxxy_xzz, ts_xxxy_xzzz, ts_xxxy_yyy, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyz, ts_xxxy_yyzz, ts_xxxy_yzz, ts_xxxy_yzzz, ts_xxxy_zzz, ts_xxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxy_xxxx[i] = 2.0 * ts_xxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxy[i] = 2.0 * ts_xxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxxz[i] = 2.0 * ts_xxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxyy[i] = 2.0 * ts_xxx_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxyz[i] = 2.0 * ts_xxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xxzz[i] = 2.0 * ts_xxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyyy[i] = 2.0 * ts_xxx_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyyz[i] = 2.0 * ts_xxx_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xyzz[i] = 2.0 * ts_xxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xzzz[i] = 2.0 * ts_xxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyyy[i] = 2.0 * ts_xxx_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyyz[i] = 2.0 * ts_xxx_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yyzz[i] = 2.0 * ts_xxx_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yzzz[i] = 2.0 * ts_xxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_zzzz[i] = 2.0 * ts_xxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 255-270 components of targeted buffer : GG

    auto gs_y_xxxz_xxxx = pbuffer.data(idx_g_gg + 255);

    auto gs_y_xxxz_xxxy = pbuffer.data(idx_g_gg + 256);

    auto gs_y_xxxz_xxxz = pbuffer.data(idx_g_gg + 257);

    auto gs_y_xxxz_xxyy = pbuffer.data(idx_g_gg + 258);

    auto gs_y_xxxz_xxyz = pbuffer.data(idx_g_gg + 259);

    auto gs_y_xxxz_xxzz = pbuffer.data(idx_g_gg + 260);

    auto gs_y_xxxz_xyyy = pbuffer.data(idx_g_gg + 261);

    auto gs_y_xxxz_xyyz = pbuffer.data(idx_g_gg + 262);

    auto gs_y_xxxz_xyzz = pbuffer.data(idx_g_gg + 263);

    auto gs_y_xxxz_xzzz = pbuffer.data(idx_g_gg + 264);

    auto gs_y_xxxz_yyyy = pbuffer.data(idx_g_gg + 265);

    auto gs_y_xxxz_yyyz = pbuffer.data(idx_g_gg + 266);

    auto gs_y_xxxz_yyzz = pbuffer.data(idx_g_gg + 267);

    auto gs_y_xxxz_yzzz = pbuffer.data(idx_g_gg + 268);

    auto gs_y_xxxz_zzzz = pbuffer.data(idx_g_gg + 269);

    #pragma omp simd aligned(gc_y, gs_y_xxxz_xxxx, gs_y_xxxz_xxxy, gs_y_xxxz_xxxz, gs_y_xxxz_xxyy, gs_y_xxxz_xxyz, gs_y_xxxz_xxzz, gs_y_xxxz_xyyy, gs_y_xxxz_xyyz, gs_y_xxxz_xyzz, gs_y_xxxz_xzzz, gs_y_xxxz_yyyy, gs_y_xxxz_yyyz, gs_y_xxxz_yyzz, gs_y_xxxz_yzzz, gs_y_xxxz_zzzz, ts_xxxz_xxx, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxy, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxz, ts_xxxz_xxzz, ts_xxxz_xyy, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyz, ts_xxxz_xyzz, ts_xxxz_xzz, ts_xxxz_xzzz, ts_xxxz_yyy, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyz, ts_xxxz_yyzz, ts_xxxz_yzz, ts_xxxz_yzzz, ts_xxxz_zzz, ts_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxz_xxxx[i] = 2.0 * ts_xxxz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxy[i] = 2.0 * ts_xxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxxz[i] = 2.0 * ts_xxxz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxyy[i] = 4.0 * ts_xxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxyz[i] = 2.0 * ts_xxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xxzz[i] = 2.0 * ts_xxxz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyyy[i] = 6.0 * ts_xxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyyz[i] = 4.0 * ts_xxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xyzz[i] = 2.0 * ts_xxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xzzz[i] = 2.0 * ts_xxxz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyyy[i] = 8.0 * ts_xxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyyz[i] = 6.0 * ts_xxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yyzz[i] = 4.0 * ts_xxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yzzz[i] = 2.0 * ts_xxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_zzzz[i] = 2.0 * ts_xxxz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 270-285 components of targeted buffer : GG

    auto gs_y_xxyy_xxxx = pbuffer.data(idx_g_gg + 270);

    auto gs_y_xxyy_xxxy = pbuffer.data(idx_g_gg + 271);

    auto gs_y_xxyy_xxxz = pbuffer.data(idx_g_gg + 272);

    auto gs_y_xxyy_xxyy = pbuffer.data(idx_g_gg + 273);

    auto gs_y_xxyy_xxyz = pbuffer.data(idx_g_gg + 274);

    auto gs_y_xxyy_xxzz = pbuffer.data(idx_g_gg + 275);

    auto gs_y_xxyy_xyyy = pbuffer.data(idx_g_gg + 276);

    auto gs_y_xxyy_xyyz = pbuffer.data(idx_g_gg + 277);

    auto gs_y_xxyy_xyzz = pbuffer.data(idx_g_gg + 278);

    auto gs_y_xxyy_xzzz = pbuffer.data(idx_g_gg + 279);

    auto gs_y_xxyy_yyyy = pbuffer.data(idx_g_gg + 280);

    auto gs_y_xxyy_yyyz = pbuffer.data(idx_g_gg + 281);

    auto gs_y_xxyy_yyzz = pbuffer.data(idx_g_gg + 282);

    auto gs_y_xxyy_yzzz = pbuffer.data(idx_g_gg + 283);

    auto gs_y_xxyy_zzzz = pbuffer.data(idx_g_gg + 284);

    #pragma omp simd aligned(gc_y, gs_y_xxyy_xxxx, gs_y_xxyy_xxxy, gs_y_xxyy_xxxz, gs_y_xxyy_xxyy, gs_y_xxyy_xxyz, gs_y_xxyy_xxzz, gs_y_xxyy_xyyy, gs_y_xxyy_xyyz, gs_y_xxyy_xyzz, gs_y_xxyy_xzzz, gs_y_xxyy_yyyy, gs_y_xxyy_yyyz, gs_y_xxyy_yyzz, gs_y_xxyy_yzzz, gs_y_xxyy_zzzz, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxzz, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyzz, ts_xxy_xzzz, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyzz, ts_xxy_yzzz, ts_xxy_zzzz, ts_xxyy_xxx, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxz, ts_xxyy_xxzz, ts_xxyy_xyy, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyz, ts_xxyy_xyzz, ts_xxyy_xzz, ts_xxyy_xzzz, ts_xxyy_yyy, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyz, ts_xxyy_yyzz, ts_xxyy_yzz, ts_xxyy_yzzz, ts_xxyy_zzz, ts_xxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyy_xxxx[i] = 4.0 * ts_xxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxy[i] = 4.0 * ts_xxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxxz[i] = 4.0 * ts_xxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxyy[i] = 4.0 * ts_xxy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxyz[i] = 4.0 * ts_xxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xxzz[i] = 4.0 * ts_xxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyyy[i] = 4.0 * ts_xxy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyyz[i] = 4.0 * ts_xxy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xyzz[i] = 4.0 * ts_xxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xzzz[i] = 4.0 * ts_xxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyyy[i] = 4.0 * ts_xxy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyyz[i] = 4.0 * ts_xxy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yyzz[i] = 4.0 * ts_xxy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yzzz[i] = 4.0 * ts_xxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_zzzz[i] = 4.0 * ts_xxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 285-300 components of targeted buffer : GG

    auto gs_y_xxyz_xxxx = pbuffer.data(idx_g_gg + 285);

    auto gs_y_xxyz_xxxy = pbuffer.data(idx_g_gg + 286);

    auto gs_y_xxyz_xxxz = pbuffer.data(idx_g_gg + 287);

    auto gs_y_xxyz_xxyy = pbuffer.data(idx_g_gg + 288);

    auto gs_y_xxyz_xxyz = pbuffer.data(idx_g_gg + 289);

    auto gs_y_xxyz_xxzz = pbuffer.data(idx_g_gg + 290);

    auto gs_y_xxyz_xyyy = pbuffer.data(idx_g_gg + 291);

    auto gs_y_xxyz_xyyz = pbuffer.data(idx_g_gg + 292);

    auto gs_y_xxyz_xyzz = pbuffer.data(idx_g_gg + 293);

    auto gs_y_xxyz_xzzz = pbuffer.data(idx_g_gg + 294);

    auto gs_y_xxyz_yyyy = pbuffer.data(idx_g_gg + 295);

    auto gs_y_xxyz_yyyz = pbuffer.data(idx_g_gg + 296);

    auto gs_y_xxyz_yyzz = pbuffer.data(idx_g_gg + 297);

    auto gs_y_xxyz_yzzz = pbuffer.data(idx_g_gg + 298);

    auto gs_y_xxyz_zzzz = pbuffer.data(idx_g_gg + 299);

    #pragma omp simd aligned(gc_y, gs_y_xxyz_xxxx, gs_y_xxyz_xxxy, gs_y_xxyz_xxxz, gs_y_xxyz_xxyy, gs_y_xxyz_xxyz, gs_y_xxyz_xxzz, gs_y_xxyz_xyyy, gs_y_xxyz_xyyz, gs_y_xxyz_xyzz, gs_y_xxyz_xzzz, gs_y_xxyz_yyyy, gs_y_xxyz_yyyz, gs_y_xxyz_yyzz, gs_y_xxyz_yzzz, gs_y_xxyz_zzzz, ts_xxyz_xxx, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxy, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxz, ts_xxyz_xxzz, ts_xxyz_xyy, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyz, ts_xxyz_xyzz, ts_xxyz_xzz, ts_xxyz_xzzz, ts_xxyz_yyy, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyz, ts_xxyz_yyzz, ts_xxyz_yzz, ts_xxyz_yzzz, ts_xxyz_zzz, ts_xxyz_zzzz, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxzz, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyzz, ts_xxz_xzzz, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyzz, ts_xxz_yzzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyz_xxxx[i] = 2.0 * ts_xxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxy[i] = 2.0 * ts_xxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxxz[i] = 2.0 * ts_xxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxyy[i] = 2.0 * ts_xxz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxyz[i] = 2.0 * ts_xxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xxzz[i] = 2.0 * ts_xxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyyy[i] = 2.0 * ts_xxz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyyz[i] = 2.0 * ts_xxz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xyzz[i] = 2.0 * ts_xxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xzzz[i] = 2.0 * ts_xxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyyy[i] = 2.0 * ts_xxz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyyz[i] = 2.0 * ts_xxz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yyzz[i] = 2.0 * ts_xxz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yzzz[i] = 2.0 * ts_xxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_zzzz[i] = 2.0 * ts_xxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 300-315 components of targeted buffer : GG

    auto gs_y_xxzz_xxxx = pbuffer.data(idx_g_gg + 300);

    auto gs_y_xxzz_xxxy = pbuffer.data(idx_g_gg + 301);

    auto gs_y_xxzz_xxxz = pbuffer.data(idx_g_gg + 302);

    auto gs_y_xxzz_xxyy = pbuffer.data(idx_g_gg + 303);

    auto gs_y_xxzz_xxyz = pbuffer.data(idx_g_gg + 304);

    auto gs_y_xxzz_xxzz = pbuffer.data(idx_g_gg + 305);

    auto gs_y_xxzz_xyyy = pbuffer.data(idx_g_gg + 306);

    auto gs_y_xxzz_xyyz = pbuffer.data(idx_g_gg + 307);

    auto gs_y_xxzz_xyzz = pbuffer.data(idx_g_gg + 308);

    auto gs_y_xxzz_xzzz = pbuffer.data(idx_g_gg + 309);

    auto gs_y_xxzz_yyyy = pbuffer.data(idx_g_gg + 310);

    auto gs_y_xxzz_yyyz = pbuffer.data(idx_g_gg + 311);

    auto gs_y_xxzz_yyzz = pbuffer.data(idx_g_gg + 312);

    auto gs_y_xxzz_yzzz = pbuffer.data(idx_g_gg + 313);

    auto gs_y_xxzz_zzzz = pbuffer.data(idx_g_gg + 314);

    #pragma omp simd aligned(gc_y, gs_y_xxzz_xxxx, gs_y_xxzz_xxxy, gs_y_xxzz_xxxz, gs_y_xxzz_xxyy, gs_y_xxzz_xxyz, gs_y_xxzz_xxzz, gs_y_xxzz_xyyy, gs_y_xxzz_xyyz, gs_y_xxzz_xyzz, gs_y_xxzz_xzzz, gs_y_xxzz_yyyy, gs_y_xxzz_yyyz, gs_y_xxzz_yyzz, gs_y_xxzz_yzzz, gs_y_xxzz_zzzz, ts_xxzz_xxx, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxy, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxz, ts_xxzz_xxzz, ts_xxzz_xyy, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyz, ts_xxzz_xyzz, ts_xxzz_xzz, ts_xxzz_xzzz, ts_xxzz_yyy, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyz, ts_xxzz_yyzz, ts_xxzz_yzz, ts_xxzz_yzzz, ts_xxzz_zzz, ts_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzz_xxxx[i] = 2.0 * ts_xxzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxy[i] = 2.0 * ts_xxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxxz[i] = 2.0 * ts_xxzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxyy[i] = 4.0 * ts_xxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxyz[i] = 2.0 * ts_xxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xxzz[i] = 2.0 * ts_xxzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyyy[i] = 6.0 * ts_xxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyyz[i] = 4.0 * ts_xxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xyzz[i] = 2.0 * ts_xxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xzzz[i] = 2.0 * ts_xxzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyyy[i] = 8.0 * ts_xxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyyz[i] = 6.0 * ts_xxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yyzz[i] = 4.0 * ts_xxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yzzz[i] = 2.0 * ts_xxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_zzzz[i] = 2.0 * ts_xxzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 315-330 components of targeted buffer : GG

    auto gs_y_xyyy_xxxx = pbuffer.data(idx_g_gg + 315);

    auto gs_y_xyyy_xxxy = pbuffer.data(idx_g_gg + 316);

    auto gs_y_xyyy_xxxz = pbuffer.data(idx_g_gg + 317);

    auto gs_y_xyyy_xxyy = pbuffer.data(idx_g_gg + 318);

    auto gs_y_xyyy_xxyz = pbuffer.data(idx_g_gg + 319);

    auto gs_y_xyyy_xxzz = pbuffer.data(idx_g_gg + 320);

    auto gs_y_xyyy_xyyy = pbuffer.data(idx_g_gg + 321);

    auto gs_y_xyyy_xyyz = pbuffer.data(idx_g_gg + 322);

    auto gs_y_xyyy_xyzz = pbuffer.data(idx_g_gg + 323);

    auto gs_y_xyyy_xzzz = pbuffer.data(idx_g_gg + 324);

    auto gs_y_xyyy_yyyy = pbuffer.data(idx_g_gg + 325);

    auto gs_y_xyyy_yyyz = pbuffer.data(idx_g_gg + 326);

    auto gs_y_xyyy_yyzz = pbuffer.data(idx_g_gg + 327);

    auto gs_y_xyyy_yzzz = pbuffer.data(idx_g_gg + 328);

    auto gs_y_xyyy_zzzz = pbuffer.data(idx_g_gg + 329);

    #pragma omp simd aligned(gc_y, gs_y_xyyy_xxxx, gs_y_xyyy_xxxy, gs_y_xyyy_xxxz, gs_y_xyyy_xxyy, gs_y_xyyy_xxyz, gs_y_xyyy_xxzz, gs_y_xyyy_xyyy, gs_y_xyyy_xyyz, gs_y_xyyy_xyzz, gs_y_xyyy_xzzz, gs_y_xyyy_yyyy, gs_y_xyyy_yyyz, gs_y_xyyy_yyzz, gs_y_xyyy_yzzz, gs_y_xyyy_zzzz, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxzz, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyzz, ts_xyy_xzzz, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyzz, ts_xyy_yzzz, ts_xyy_zzzz, ts_xyyy_xxx, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxy, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxz, ts_xyyy_xxzz, ts_xyyy_xyy, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyz, ts_xyyy_xyzz, ts_xyyy_xzz, ts_xyyy_xzzz, ts_xyyy_yyy, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyz, ts_xyyy_yyzz, ts_xyyy_yzz, ts_xyyy_yzzz, ts_xyyy_zzz, ts_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyy_xxxx[i] = 6.0 * ts_xyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxy[i] = 6.0 * ts_xyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxxz[i] = 6.0 * ts_xyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxyy[i] = 6.0 * ts_xyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxyz[i] = 6.0 * ts_xyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xxzz[i] = 6.0 * ts_xyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyyy[i] = 6.0 * ts_xyy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyyz[i] = 6.0 * ts_xyy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xyzz[i] = 6.0 * ts_xyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xzzz[i] = 6.0 * ts_xyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyyy[i] = 6.0 * ts_xyy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyyz[i] = 6.0 * ts_xyy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yyzz[i] = 6.0 * ts_xyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yzzz[i] = 6.0 * ts_xyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_zzzz[i] = 6.0 * ts_xyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 330-345 components of targeted buffer : GG

    auto gs_y_xyyz_xxxx = pbuffer.data(idx_g_gg + 330);

    auto gs_y_xyyz_xxxy = pbuffer.data(idx_g_gg + 331);

    auto gs_y_xyyz_xxxz = pbuffer.data(idx_g_gg + 332);

    auto gs_y_xyyz_xxyy = pbuffer.data(idx_g_gg + 333);

    auto gs_y_xyyz_xxyz = pbuffer.data(idx_g_gg + 334);

    auto gs_y_xyyz_xxzz = pbuffer.data(idx_g_gg + 335);

    auto gs_y_xyyz_xyyy = pbuffer.data(idx_g_gg + 336);

    auto gs_y_xyyz_xyyz = pbuffer.data(idx_g_gg + 337);

    auto gs_y_xyyz_xyzz = pbuffer.data(idx_g_gg + 338);

    auto gs_y_xyyz_xzzz = pbuffer.data(idx_g_gg + 339);

    auto gs_y_xyyz_yyyy = pbuffer.data(idx_g_gg + 340);

    auto gs_y_xyyz_yyyz = pbuffer.data(idx_g_gg + 341);

    auto gs_y_xyyz_yyzz = pbuffer.data(idx_g_gg + 342);

    auto gs_y_xyyz_yzzz = pbuffer.data(idx_g_gg + 343);

    auto gs_y_xyyz_zzzz = pbuffer.data(idx_g_gg + 344);

    #pragma omp simd aligned(gc_y, gs_y_xyyz_xxxx, gs_y_xyyz_xxxy, gs_y_xyyz_xxxz, gs_y_xyyz_xxyy, gs_y_xyyz_xxyz, gs_y_xyyz_xxzz, gs_y_xyyz_xyyy, gs_y_xyyz_xyyz, gs_y_xyyz_xyzz, gs_y_xyyz_xzzz, gs_y_xyyz_yyyy, gs_y_xyyz_yyyz, gs_y_xyyz_yyzz, gs_y_xyyz_yzzz, gs_y_xyyz_zzzz, ts_xyyz_xxx, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxy, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxz, ts_xyyz_xxzz, ts_xyyz_xyy, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyz, ts_xyyz_xyzz, ts_xyyz_xzz, ts_xyyz_xzzz, ts_xyyz_yyy, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyz, ts_xyyz_yyzz, ts_xyyz_yzz, ts_xyyz_yzzz, ts_xyyz_zzz, ts_xyyz_zzzz, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxzz, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyzz, ts_xyz_xzzz, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyzz, ts_xyz_yzzz, ts_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyz_xxxx[i] = 4.0 * ts_xyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxy[i] = 4.0 * ts_xyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxxz[i] = 4.0 * ts_xyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxyy[i] = 4.0 * ts_xyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxyz[i] = 4.0 * ts_xyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xxzz[i] = 4.0 * ts_xyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyyy[i] = 4.0 * ts_xyz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyyz[i] = 4.0 * ts_xyz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xyzz[i] = 4.0 * ts_xyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xzzz[i] = 4.0 * ts_xyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyyy[i] = 4.0 * ts_xyz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyyz[i] = 4.0 * ts_xyz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yyzz[i] = 4.0 * ts_xyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yzzz[i] = 4.0 * ts_xyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_zzzz[i] = 4.0 * ts_xyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 345-360 components of targeted buffer : GG

    auto gs_y_xyzz_xxxx = pbuffer.data(idx_g_gg + 345);

    auto gs_y_xyzz_xxxy = pbuffer.data(idx_g_gg + 346);

    auto gs_y_xyzz_xxxz = pbuffer.data(idx_g_gg + 347);

    auto gs_y_xyzz_xxyy = pbuffer.data(idx_g_gg + 348);

    auto gs_y_xyzz_xxyz = pbuffer.data(idx_g_gg + 349);

    auto gs_y_xyzz_xxzz = pbuffer.data(idx_g_gg + 350);

    auto gs_y_xyzz_xyyy = pbuffer.data(idx_g_gg + 351);

    auto gs_y_xyzz_xyyz = pbuffer.data(idx_g_gg + 352);

    auto gs_y_xyzz_xyzz = pbuffer.data(idx_g_gg + 353);

    auto gs_y_xyzz_xzzz = pbuffer.data(idx_g_gg + 354);

    auto gs_y_xyzz_yyyy = pbuffer.data(idx_g_gg + 355);

    auto gs_y_xyzz_yyyz = pbuffer.data(idx_g_gg + 356);

    auto gs_y_xyzz_yyzz = pbuffer.data(idx_g_gg + 357);

    auto gs_y_xyzz_yzzz = pbuffer.data(idx_g_gg + 358);

    auto gs_y_xyzz_zzzz = pbuffer.data(idx_g_gg + 359);

    #pragma omp simd aligned(gc_y, gs_y_xyzz_xxxx, gs_y_xyzz_xxxy, gs_y_xyzz_xxxz, gs_y_xyzz_xxyy, gs_y_xyzz_xxyz, gs_y_xyzz_xxzz, gs_y_xyzz_xyyy, gs_y_xyzz_xyyz, gs_y_xyzz_xyzz, gs_y_xyzz_xzzz, gs_y_xyzz_yyyy, gs_y_xyzz_yyyz, gs_y_xyzz_yyzz, gs_y_xyzz_yzzz, gs_y_xyzz_zzzz, ts_xyzz_xxx, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxy, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxz, ts_xyzz_xxzz, ts_xyzz_xyy, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyz, ts_xyzz_xyzz, ts_xyzz_xzz, ts_xyzz_xzzz, ts_xyzz_yyy, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyz, ts_xyzz_yyzz, ts_xyzz_yzz, ts_xyzz_yzzz, ts_xyzz_zzz, ts_xyzz_zzzz, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxzz, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyzz, ts_xzz_xzzz, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyzz, ts_xzz_yzzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzz_xxxx[i] = 2.0 * ts_xzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxy[i] = 2.0 * ts_xzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxxz[i] = 2.0 * ts_xzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxyy[i] = 2.0 * ts_xzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxyz[i] = 2.0 * ts_xzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xxzz[i] = 2.0 * ts_xzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyyy[i] = 2.0 * ts_xzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyyz[i] = 2.0 * ts_xzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xyzz[i] = 2.0 * ts_xzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xzzz[i] = 2.0 * ts_xzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyyy[i] = 2.0 * ts_xzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyyz[i] = 2.0 * ts_xzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yyzz[i] = 2.0 * ts_xzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yzzz[i] = 2.0 * ts_xzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_zzzz[i] = 2.0 * ts_xzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 360-375 components of targeted buffer : GG

    auto gs_y_xzzz_xxxx = pbuffer.data(idx_g_gg + 360);

    auto gs_y_xzzz_xxxy = pbuffer.data(idx_g_gg + 361);

    auto gs_y_xzzz_xxxz = pbuffer.data(idx_g_gg + 362);

    auto gs_y_xzzz_xxyy = pbuffer.data(idx_g_gg + 363);

    auto gs_y_xzzz_xxyz = pbuffer.data(idx_g_gg + 364);

    auto gs_y_xzzz_xxzz = pbuffer.data(idx_g_gg + 365);

    auto gs_y_xzzz_xyyy = pbuffer.data(idx_g_gg + 366);

    auto gs_y_xzzz_xyyz = pbuffer.data(idx_g_gg + 367);

    auto gs_y_xzzz_xyzz = pbuffer.data(idx_g_gg + 368);

    auto gs_y_xzzz_xzzz = pbuffer.data(idx_g_gg + 369);

    auto gs_y_xzzz_yyyy = pbuffer.data(idx_g_gg + 370);

    auto gs_y_xzzz_yyyz = pbuffer.data(idx_g_gg + 371);

    auto gs_y_xzzz_yyzz = pbuffer.data(idx_g_gg + 372);

    auto gs_y_xzzz_yzzz = pbuffer.data(idx_g_gg + 373);

    auto gs_y_xzzz_zzzz = pbuffer.data(idx_g_gg + 374);

    #pragma omp simd aligned(gc_y, gs_y_xzzz_xxxx, gs_y_xzzz_xxxy, gs_y_xzzz_xxxz, gs_y_xzzz_xxyy, gs_y_xzzz_xxyz, gs_y_xzzz_xxzz, gs_y_xzzz_xyyy, gs_y_xzzz_xyyz, gs_y_xzzz_xyzz, gs_y_xzzz_xzzz, gs_y_xzzz_yyyy, gs_y_xzzz_yyyz, gs_y_xzzz_yyzz, gs_y_xzzz_yzzz, gs_y_xzzz_zzzz, ts_xzzz_xxx, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxy, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxz, ts_xzzz_xxzz, ts_xzzz_xyy, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyz, ts_xzzz_xyzz, ts_xzzz_xzz, ts_xzzz_xzzz, ts_xzzz_yyy, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyz, ts_xzzz_yyzz, ts_xzzz_yzz, ts_xzzz_yzzz, ts_xzzz_zzz, ts_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzz_xxxx[i] = 2.0 * ts_xzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxy[i] = 2.0 * ts_xzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxxz[i] = 2.0 * ts_xzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxyy[i] = 4.0 * ts_xzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxyz[i] = 2.0 * ts_xzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xxzz[i] = 2.0 * ts_xzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyyy[i] = 6.0 * ts_xzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyyz[i] = 4.0 * ts_xzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xyzz[i] = 2.0 * ts_xzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xzzz[i] = 2.0 * ts_xzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyyy[i] = 8.0 * ts_xzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyyz[i] = 6.0 * ts_xzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yyzz[i] = 4.0 * ts_xzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yzzz[i] = 2.0 * ts_xzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_zzzz[i] = 2.0 * ts_xzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 375-390 components of targeted buffer : GG

    auto gs_y_yyyy_xxxx = pbuffer.data(idx_g_gg + 375);

    auto gs_y_yyyy_xxxy = pbuffer.data(idx_g_gg + 376);

    auto gs_y_yyyy_xxxz = pbuffer.data(idx_g_gg + 377);

    auto gs_y_yyyy_xxyy = pbuffer.data(idx_g_gg + 378);

    auto gs_y_yyyy_xxyz = pbuffer.data(idx_g_gg + 379);

    auto gs_y_yyyy_xxzz = pbuffer.data(idx_g_gg + 380);

    auto gs_y_yyyy_xyyy = pbuffer.data(idx_g_gg + 381);

    auto gs_y_yyyy_xyyz = pbuffer.data(idx_g_gg + 382);

    auto gs_y_yyyy_xyzz = pbuffer.data(idx_g_gg + 383);

    auto gs_y_yyyy_xzzz = pbuffer.data(idx_g_gg + 384);

    auto gs_y_yyyy_yyyy = pbuffer.data(idx_g_gg + 385);

    auto gs_y_yyyy_yyyz = pbuffer.data(idx_g_gg + 386);

    auto gs_y_yyyy_yyzz = pbuffer.data(idx_g_gg + 387);

    auto gs_y_yyyy_yzzz = pbuffer.data(idx_g_gg + 388);

    auto gs_y_yyyy_zzzz = pbuffer.data(idx_g_gg + 389);

    #pragma omp simd aligned(gc_y, gs_y_yyyy_xxxx, gs_y_yyyy_xxxy, gs_y_yyyy_xxxz, gs_y_yyyy_xxyy, gs_y_yyyy_xxyz, gs_y_yyyy_xxzz, gs_y_yyyy_xyyy, gs_y_yyyy_xyyz, gs_y_yyyy_xyzz, gs_y_yyyy_xzzz, gs_y_yyyy_yyyy, gs_y_yyyy_yyyz, gs_y_yyyy_yyzz, gs_y_yyyy_yzzz, gs_y_yyyy_zzzz, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxzz, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyzz, ts_yyy_xzzz, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyzz, ts_yyy_yzzz, ts_yyy_zzzz, ts_yyyy_xxx, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxz, ts_yyyy_xxzz, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_xzz, ts_yyyy_xzzz, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyy_zzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyy_xxxx[i] = 8.0 * ts_yyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxy[i] = 8.0 * ts_yyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxxz[i] = 8.0 * ts_yyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxyy[i] = 8.0 * ts_yyy_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxyz[i] = 8.0 * ts_yyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xxzz[i] = 8.0 * ts_yyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyyy[i] = 8.0 * ts_yyy_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyyz[i] = 8.0 * ts_yyy_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xyzz[i] = 8.0 * ts_yyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xzzz[i] = 8.0 * ts_yyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyyy[i] = 8.0 * ts_yyy_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyyz[i] = 8.0 * ts_yyy_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yyzz[i] = 8.0 * ts_yyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yzzz[i] = 8.0 * ts_yyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_zzzz[i] = 8.0 * ts_yyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 390-405 components of targeted buffer : GG

    auto gs_y_yyyz_xxxx = pbuffer.data(idx_g_gg + 390);

    auto gs_y_yyyz_xxxy = pbuffer.data(idx_g_gg + 391);

    auto gs_y_yyyz_xxxz = pbuffer.data(idx_g_gg + 392);

    auto gs_y_yyyz_xxyy = pbuffer.data(idx_g_gg + 393);

    auto gs_y_yyyz_xxyz = pbuffer.data(idx_g_gg + 394);

    auto gs_y_yyyz_xxzz = pbuffer.data(idx_g_gg + 395);

    auto gs_y_yyyz_xyyy = pbuffer.data(idx_g_gg + 396);

    auto gs_y_yyyz_xyyz = pbuffer.data(idx_g_gg + 397);

    auto gs_y_yyyz_xyzz = pbuffer.data(idx_g_gg + 398);

    auto gs_y_yyyz_xzzz = pbuffer.data(idx_g_gg + 399);

    auto gs_y_yyyz_yyyy = pbuffer.data(idx_g_gg + 400);

    auto gs_y_yyyz_yyyz = pbuffer.data(idx_g_gg + 401);

    auto gs_y_yyyz_yyzz = pbuffer.data(idx_g_gg + 402);

    auto gs_y_yyyz_yzzz = pbuffer.data(idx_g_gg + 403);

    auto gs_y_yyyz_zzzz = pbuffer.data(idx_g_gg + 404);

    #pragma omp simd aligned(gc_y, gs_y_yyyz_xxxx, gs_y_yyyz_xxxy, gs_y_yyyz_xxxz, gs_y_yyyz_xxyy, gs_y_yyyz_xxyz, gs_y_yyyz_xxzz, gs_y_yyyz_xyyy, gs_y_yyyz_xyyz, gs_y_yyyz_xyzz, gs_y_yyyz_xzzz, gs_y_yyyz_yyyy, gs_y_yyyz_yyyz, gs_y_yyyz_yyzz, gs_y_yyyz_yzzz, gs_y_yyyz_zzzz, ts_yyyz_xxx, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxy, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxz, ts_yyyz_xxzz, ts_yyyz_xyy, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyz, ts_yyyz_xyzz, ts_yyyz_xzz, ts_yyyz_xzzz, ts_yyyz_yyy, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyz, ts_yyyz_yyzz, ts_yyyz_yzz, ts_yyyz_yzzz, ts_yyyz_zzz, ts_yyyz_zzzz, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxzz, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyzz, ts_yyz_xzzz, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyzz, ts_yyz_yzzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyz_xxxx[i] = 6.0 * ts_yyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxy[i] = 6.0 * ts_yyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxxz[i] = 6.0 * ts_yyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxyy[i] = 6.0 * ts_yyz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxyz[i] = 6.0 * ts_yyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xxzz[i] = 6.0 * ts_yyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyyy[i] = 6.0 * ts_yyz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyyz[i] = 6.0 * ts_yyz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xyzz[i] = 6.0 * ts_yyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xzzz[i] = 6.0 * ts_yyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyyy[i] = 6.0 * ts_yyz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyyz[i] = 6.0 * ts_yyz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yyzz[i] = 6.0 * ts_yyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yzzz[i] = 6.0 * ts_yyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_zzzz[i] = 6.0 * ts_yyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 405-420 components of targeted buffer : GG

    auto gs_y_yyzz_xxxx = pbuffer.data(idx_g_gg + 405);

    auto gs_y_yyzz_xxxy = pbuffer.data(idx_g_gg + 406);

    auto gs_y_yyzz_xxxz = pbuffer.data(idx_g_gg + 407);

    auto gs_y_yyzz_xxyy = pbuffer.data(idx_g_gg + 408);

    auto gs_y_yyzz_xxyz = pbuffer.data(idx_g_gg + 409);

    auto gs_y_yyzz_xxzz = pbuffer.data(idx_g_gg + 410);

    auto gs_y_yyzz_xyyy = pbuffer.data(idx_g_gg + 411);

    auto gs_y_yyzz_xyyz = pbuffer.data(idx_g_gg + 412);

    auto gs_y_yyzz_xyzz = pbuffer.data(idx_g_gg + 413);

    auto gs_y_yyzz_xzzz = pbuffer.data(idx_g_gg + 414);

    auto gs_y_yyzz_yyyy = pbuffer.data(idx_g_gg + 415);

    auto gs_y_yyzz_yyyz = pbuffer.data(idx_g_gg + 416);

    auto gs_y_yyzz_yyzz = pbuffer.data(idx_g_gg + 417);

    auto gs_y_yyzz_yzzz = pbuffer.data(idx_g_gg + 418);

    auto gs_y_yyzz_zzzz = pbuffer.data(idx_g_gg + 419);

    #pragma omp simd aligned(gc_y, gs_y_yyzz_xxxx, gs_y_yyzz_xxxy, gs_y_yyzz_xxxz, gs_y_yyzz_xxyy, gs_y_yyzz_xxyz, gs_y_yyzz_xxzz, gs_y_yyzz_xyyy, gs_y_yyzz_xyyz, gs_y_yyzz_xyzz, gs_y_yyzz_xzzz, gs_y_yyzz_yyyy, gs_y_yyzz_yyyz, gs_y_yyzz_yyzz, gs_y_yyzz_yzzz, gs_y_yyzz_zzzz, ts_yyzz_xxx, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxy, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxz, ts_yyzz_xxzz, ts_yyzz_xyy, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyz, ts_yyzz_xyzz, ts_yyzz_xzz, ts_yyzz_xzzz, ts_yyzz_yyy, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyz, ts_yyzz_yyzz, ts_yyzz_yzz, ts_yyzz_yzzz, ts_yyzz_zzz, ts_yyzz_zzzz, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxzz, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyzz, ts_yzz_xzzz, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyzz, ts_yzz_yzzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzz_xxxx[i] = 4.0 * ts_yzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxy[i] = 4.0 * ts_yzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxxz[i] = 4.0 * ts_yzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxyy[i] = 4.0 * ts_yzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxyz[i] = 4.0 * ts_yzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xxzz[i] = 4.0 * ts_yzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyyy[i] = 4.0 * ts_yzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyyz[i] = 4.0 * ts_yzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xyzz[i] = 4.0 * ts_yzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xzzz[i] = 4.0 * ts_yzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyyy[i] = 4.0 * ts_yzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyyz[i] = 4.0 * ts_yzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yyzz[i] = 4.0 * ts_yzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yzzz[i] = 4.0 * ts_yzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_zzzz[i] = 4.0 * ts_yzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 420-435 components of targeted buffer : GG

    auto gs_y_yzzz_xxxx = pbuffer.data(idx_g_gg + 420);

    auto gs_y_yzzz_xxxy = pbuffer.data(idx_g_gg + 421);

    auto gs_y_yzzz_xxxz = pbuffer.data(idx_g_gg + 422);

    auto gs_y_yzzz_xxyy = pbuffer.data(idx_g_gg + 423);

    auto gs_y_yzzz_xxyz = pbuffer.data(idx_g_gg + 424);

    auto gs_y_yzzz_xxzz = pbuffer.data(idx_g_gg + 425);

    auto gs_y_yzzz_xyyy = pbuffer.data(idx_g_gg + 426);

    auto gs_y_yzzz_xyyz = pbuffer.data(idx_g_gg + 427);

    auto gs_y_yzzz_xyzz = pbuffer.data(idx_g_gg + 428);

    auto gs_y_yzzz_xzzz = pbuffer.data(idx_g_gg + 429);

    auto gs_y_yzzz_yyyy = pbuffer.data(idx_g_gg + 430);

    auto gs_y_yzzz_yyyz = pbuffer.data(idx_g_gg + 431);

    auto gs_y_yzzz_yyzz = pbuffer.data(idx_g_gg + 432);

    auto gs_y_yzzz_yzzz = pbuffer.data(idx_g_gg + 433);

    auto gs_y_yzzz_zzzz = pbuffer.data(idx_g_gg + 434);

    #pragma omp simd aligned(gc_y, gs_y_yzzz_xxxx, gs_y_yzzz_xxxy, gs_y_yzzz_xxxz, gs_y_yzzz_xxyy, gs_y_yzzz_xxyz, gs_y_yzzz_xxzz, gs_y_yzzz_xyyy, gs_y_yzzz_xyyz, gs_y_yzzz_xyzz, gs_y_yzzz_xzzz, gs_y_yzzz_yyyy, gs_y_yzzz_yyyz, gs_y_yzzz_yyzz, gs_y_yzzz_yzzz, gs_y_yzzz_zzzz, ts_yzzz_xxx, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxy, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxz, ts_yzzz_xxzz, ts_yzzz_xyy, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyz, ts_yzzz_xyzz, ts_yzzz_xzz, ts_yzzz_xzzz, ts_yzzz_yyy, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyz, ts_yzzz_yyzz, ts_yzzz_yzz, ts_yzzz_yzzz, ts_yzzz_zzz, ts_yzzz_zzzz, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxzz, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyzz, ts_zzz_xzzz, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyzz, ts_zzz_yzzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzz_xxxx[i] = 2.0 * ts_zzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxy[i] = 2.0 * ts_zzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxxz[i] = 2.0 * ts_zzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxyy[i] = 2.0 * ts_zzz_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxyz[i] = 2.0 * ts_zzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xxzz[i] = 2.0 * ts_zzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyyy[i] = 2.0 * ts_zzz_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyyz[i] = 2.0 * ts_zzz_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xyzz[i] = 2.0 * ts_zzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xzzz[i] = 2.0 * ts_zzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyyy[i] = 2.0 * ts_zzz_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyyz[i] = 2.0 * ts_zzz_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yyzz[i] = 2.0 * ts_zzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yzzz[i] = 2.0 * ts_zzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_zzzz[i] = 2.0 * ts_zzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 435-450 components of targeted buffer : GG

    auto gs_y_zzzz_xxxx = pbuffer.data(idx_g_gg + 435);

    auto gs_y_zzzz_xxxy = pbuffer.data(idx_g_gg + 436);

    auto gs_y_zzzz_xxxz = pbuffer.data(idx_g_gg + 437);

    auto gs_y_zzzz_xxyy = pbuffer.data(idx_g_gg + 438);

    auto gs_y_zzzz_xxyz = pbuffer.data(idx_g_gg + 439);

    auto gs_y_zzzz_xxzz = pbuffer.data(idx_g_gg + 440);

    auto gs_y_zzzz_xyyy = pbuffer.data(idx_g_gg + 441);

    auto gs_y_zzzz_xyyz = pbuffer.data(idx_g_gg + 442);

    auto gs_y_zzzz_xyzz = pbuffer.data(idx_g_gg + 443);

    auto gs_y_zzzz_xzzz = pbuffer.data(idx_g_gg + 444);

    auto gs_y_zzzz_yyyy = pbuffer.data(idx_g_gg + 445);

    auto gs_y_zzzz_yyyz = pbuffer.data(idx_g_gg + 446);

    auto gs_y_zzzz_yyzz = pbuffer.data(idx_g_gg + 447);

    auto gs_y_zzzz_yzzz = pbuffer.data(idx_g_gg + 448);

    auto gs_y_zzzz_zzzz = pbuffer.data(idx_g_gg + 449);

    #pragma omp simd aligned(gc_y, gs_y_zzzz_xxxx, gs_y_zzzz_xxxy, gs_y_zzzz_xxxz, gs_y_zzzz_xxyy, gs_y_zzzz_xxyz, gs_y_zzzz_xxzz, gs_y_zzzz_xyyy, gs_y_zzzz_xyyz, gs_y_zzzz_xyzz, gs_y_zzzz_xzzz, gs_y_zzzz_yyyy, gs_y_zzzz_yyyz, gs_y_zzzz_yyzz, gs_y_zzzz_yzzz, gs_y_zzzz_zzzz, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzz_xxxx[i] = 2.0 * ts_zzzz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxy[i] = 2.0 * ts_zzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxxz[i] = 2.0 * ts_zzzz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxyy[i] = 4.0 * ts_zzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxyz[i] = 2.0 * ts_zzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xxzz[i] = 2.0 * ts_zzzz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyyy[i] = 6.0 * ts_zzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyyz[i] = 4.0 * ts_zzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xyzz[i] = 2.0 * ts_zzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xzzz[i] = 2.0 * ts_zzzz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyyy[i] = 8.0 * ts_zzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyyz[i] = 6.0 * ts_zzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yyzz[i] = 4.0 * ts_zzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yzzz[i] = 2.0 * ts_zzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_zzzz[i] = 2.0 * ts_zzzz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 450-465 components of targeted buffer : GG

    auto gs_z_xxxx_xxxx = pbuffer.data(idx_g_gg + 450);

    auto gs_z_xxxx_xxxy = pbuffer.data(idx_g_gg + 451);

    auto gs_z_xxxx_xxxz = pbuffer.data(idx_g_gg + 452);

    auto gs_z_xxxx_xxyy = pbuffer.data(idx_g_gg + 453);

    auto gs_z_xxxx_xxyz = pbuffer.data(idx_g_gg + 454);

    auto gs_z_xxxx_xxzz = pbuffer.data(idx_g_gg + 455);

    auto gs_z_xxxx_xyyy = pbuffer.data(idx_g_gg + 456);

    auto gs_z_xxxx_xyyz = pbuffer.data(idx_g_gg + 457);

    auto gs_z_xxxx_xyzz = pbuffer.data(idx_g_gg + 458);

    auto gs_z_xxxx_xzzz = pbuffer.data(idx_g_gg + 459);

    auto gs_z_xxxx_yyyy = pbuffer.data(idx_g_gg + 460);

    auto gs_z_xxxx_yyyz = pbuffer.data(idx_g_gg + 461);

    auto gs_z_xxxx_yyzz = pbuffer.data(idx_g_gg + 462);

    auto gs_z_xxxx_yzzz = pbuffer.data(idx_g_gg + 463);

    auto gs_z_xxxx_zzzz = pbuffer.data(idx_g_gg + 464);

    #pragma omp simd aligned(gc_z, gs_z_xxxx_xxxx, gs_z_xxxx_xxxy, gs_z_xxxx_xxxz, gs_z_xxxx_xxyy, gs_z_xxxx_xxyz, gs_z_xxxx_xxzz, gs_z_xxxx_xyyy, gs_z_xxxx_xyyz, gs_z_xxxx_xyzz, gs_z_xxxx_xzzz, gs_z_xxxx_yyyy, gs_z_xxxx_yyyz, gs_z_xxxx_yyzz, gs_z_xxxx_yzzz, gs_z_xxxx_zzzz, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_yyy, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyz, ts_xxxx_yyzz, ts_xxxx_yzz, ts_xxxx_yzzz, ts_xxxx_zzz, ts_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxx_xxxx[i] = 2.0 * ts_xxxx_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxy[i] = 2.0 * ts_xxxx_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxxz[i] = 2.0 * ts_xxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxyy[i] = 2.0 * ts_xxxx_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxyz[i] = 2.0 * ts_xxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xxzz[i] = 4.0 * ts_xxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyyy[i] = 2.0 * ts_xxxx_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyyz[i] = 2.0 * ts_xxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xyzz[i] = 4.0 * ts_xxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xzzz[i] = 6.0 * ts_xxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyyy[i] = 2.0 * ts_xxxx_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyyz[i] = 2.0 * ts_xxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yyzz[i] = 4.0 * ts_xxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yzzz[i] = 6.0 * ts_xxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_zzzz[i] = 8.0 * ts_xxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 465-480 components of targeted buffer : GG

    auto gs_z_xxxy_xxxx = pbuffer.data(idx_g_gg + 465);

    auto gs_z_xxxy_xxxy = pbuffer.data(idx_g_gg + 466);

    auto gs_z_xxxy_xxxz = pbuffer.data(idx_g_gg + 467);

    auto gs_z_xxxy_xxyy = pbuffer.data(idx_g_gg + 468);

    auto gs_z_xxxy_xxyz = pbuffer.data(idx_g_gg + 469);

    auto gs_z_xxxy_xxzz = pbuffer.data(idx_g_gg + 470);

    auto gs_z_xxxy_xyyy = pbuffer.data(idx_g_gg + 471);

    auto gs_z_xxxy_xyyz = pbuffer.data(idx_g_gg + 472);

    auto gs_z_xxxy_xyzz = pbuffer.data(idx_g_gg + 473);

    auto gs_z_xxxy_xzzz = pbuffer.data(idx_g_gg + 474);

    auto gs_z_xxxy_yyyy = pbuffer.data(idx_g_gg + 475);

    auto gs_z_xxxy_yyyz = pbuffer.data(idx_g_gg + 476);

    auto gs_z_xxxy_yyzz = pbuffer.data(idx_g_gg + 477);

    auto gs_z_xxxy_yzzz = pbuffer.data(idx_g_gg + 478);

    auto gs_z_xxxy_zzzz = pbuffer.data(idx_g_gg + 479);

    #pragma omp simd aligned(gc_z, gs_z_xxxy_xxxx, gs_z_xxxy_xxxy, gs_z_xxxy_xxxz, gs_z_xxxy_xxyy, gs_z_xxxy_xxyz, gs_z_xxxy_xxzz, gs_z_xxxy_xyyy, gs_z_xxxy_xyyz, gs_z_xxxy_xyzz, gs_z_xxxy_xzzz, gs_z_xxxy_yyyy, gs_z_xxxy_yyyz, gs_z_xxxy_yyzz, gs_z_xxxy_yzzz, gs_z_xxxy_zzzz, ts_xxxy_xxx, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxy, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxz, ts_xxxy_xxzz, ts_xxxy_xyy, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyz, ts_xxxy_xyzz, ts_xxxy_xzz, ts_xxxy_xzzz, ts_xxxy_yyy, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyz, ts_xxxy_yyzz, ts_xxxy_yzz, ts_xxxy_yzzz, ts_xxxy_zzz, ts_xxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxy_xxxx[i] = 2.0 * ts_xxxy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxy[i] = 2.0 * ts_xxxy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxxz[i] = 2.0 * ts_xxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxyy[i] = 2.0 * ts_xxxy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxyz[i] = 2.0 * ts_xxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xxzz[i] = 4.0 * ts_xxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyyy[i] = 2.0 * ts_xxxy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyyz[i] = 2.0 * ts_xxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xyzz[i] = 4.0 * ts_xxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xzzz[i] = 6.0 * ts_xxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyyy[i] = 2.0 * ts_xxxy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyyz[i] = 2.0 * ts_xxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yyzz[i] = 4.0 * ts_xxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yzzz[i] = 6.0 * ts_xxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_zzzz[i] = 8.0 * ts_xxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 480-495 components of targeted buffer : GG

    auto gs_z_xxxz_xxxx = pbuffer.data(idx_g_gg + 480);

    auto gs_z_xxxz_xxxy = pbuffer.data(idx_g_gg + 481);

    auto gs_z_xxxz_xxxz = pbuffer.data(idx_g_gg + 482);

    auto gs_z_xxxz_xxyy = pbuffer.data(idx_g_gg + 483);

    auto gs_z_xxxz_xxyz = pbuffer.data(idx_g_gg + 484);

    auto gs_z_xxxz_xxzz = pbuffer.data(idx_g_gg + 485);

    auto gs_z_xxxz_xyyy = pbuffer.data(idx_g_gg + 486);

    auto gs_z_xxxz_xyyz = pbuffer.data(idx_g_gg + 487);

    auto gs_z_xxxz_xyzz = pbuffer.data(idx_g_gg + 488);

    auto gs_z_xxxz_xzzz = pbuffer.data(idx_g_gg + 489);

    auto gs_z_xxxz_yyyy = pbuffer.data(idx_g_gg + 490);

    auto gs_z_xxxz_yyyz = pbuffer.data(idx_g_gg + 491);

    auto gs_z_xxxz_yyzz = pbuffer.data(idx_g_gg + 492);

    auto gs_z_xxxz_yzzz = pbuffer.data(idx_g_gg + 493);

    auto gs_z_xxxz_zzzz = pbuffer.data(idx_g_gg + 494);

    #pragma omp simd aligned(gc_z, gs_z_xxxz_xxxx, gs_z_xxxz_xxxy, gs_z_xxxz_xxxz, gs_z_xxxz_xxyy, gs_z_xxxz_xxyz, gs_z_xxxz_xxzz, gs_z_xxxz_xyyy, gs_z_xxxz_xyyz, gs_z_xxxz_xyzz, gs_z_xxxz_xzzz, gs_z_xxxz_yyyy, gs_z_xxxz_yyyz, gs_z_xxxz_yyzz, gs_z_xxxz_yzzz, gs_z_xxxz_zzzz, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxzz, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyzz, ts_xxx_xzzz, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyzz, ts_xxx_yzzz, ts_xxx_zzzz, ts_xxxz_xxx, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxy, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxz, ts_xxxz_xxzz, ts_xxxz_xyy, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyz, ts_xxxz_xyzz, ts_xxxz_xzz, ts_xxxz_xzzz, ts_xxxz_yyy, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyz, ts_xxxz_yyzz, ts_xxxz_yzz, ts_xxxz_yzzz, ts_xxxz_zzz, ts_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxz_xxxx[i] = 2.0 * ts_xxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxy[i] = 2.0 * ts_xxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxxz[i] = 2.0 * ts_xxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxyy[i] = 2.0 * ts_xxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxyz[i] = 2.0 * ts_xxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xxzz[i] = 2.0 * ts_xxx_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyyy[i] = 2.0 * ts_xxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyyz[i] = 2.0 * ts_xxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xyzz[i] = 2.0 * ts_xxx_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xzzz[i] = 2.0 * ts_xxx_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyyy[i] = 2.0 * ts_xxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyyz[i] = 2.0 * ts_xxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yyzz[i] = 2.0 * ts_xxx_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yzzz[i] = 2.0 * ts_xxx_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_zzzz[i] = 2.0 * ts_xxx_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 495-510 components of targeted buffer : GG

    auto gs_z_xxyy_xxxx = pbuffer.data(idx_g_gg + 495);

    auto gs_z_xxyy_xxxy = pbuffer.data(idx_g_gg + 496);

    auto gs_z_xxyy_xxxz = pbuffer.data(idx_g_gg + 497);

    auto gs_z_xxyy_xxyy = pbuffer.data(idx_g_gg + 498);

    auto gs_z_xxyy_xxyz = pbuffer.data(idx_g_gg + 499);

    auto gs_z_xxyy_xxzz = pbuffer.data(idx_g_gg + 500);

    auto gs_z_xxyy_xyyy = pbuffer.data(idx_g_gg + 501);

    auto gs_z_xxyy_xyyz = pbuffer.data(idx_g_gg + 502);

    auto gs_z_xxyy_xyzz = pbuffer.data(idx_g_gg + 503);

    auto gs_z_xxyy_xzzz = pbuffer.data(idx_g_gg + 504);

    auto gs_z_xxyy_yyyy = pbuffer.data(idx_g_gg + 505);

    auto gs_z_xxyy_yyyz = pbuffer.data(idx_g_gg + 506);

    auto gs_z_xxyy_yyzz = pbuffer.data(idx_g_gg + 507);

    auto gs_z_xxyy_yzzz = pbuffer.data(idx_g_gg + 508);

    auto gs_z_xxyy_zzzz = pbuffer.data(idx_g_gg + 509);

    #pragma omp simd aligned(gc_z, gs_z_xxyy_xxxx, gs_z_xxyy_xxxy, gs_z_xxyy_xxxz, gs_z_xxyy_xxyy, gs_z_xxyy_xxyz, gs_z_xxyy_xxzz, gs_z_xxyy_xyyy, gs_z_xxyy_xyyz, gs_z_xxyy_xyzz, gs_z_xxyy_xzzz, gs_z_xxyy_yyyy, gs_z_xxyy_yyyz, gs_z_xxyy_yyzz, gs_z_xxyy_yzzz, gs_z_xxyy_zzzz, ts_xxyy_xxx, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxz, ts_xxyy_xxzz, ts_xxyy_xyy, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyz, ts_xxyy_xyzz, ts_xxyy_xzz, ts_xxyy_xzzz, ts_xxyy_yyy, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyz, ts_xxyy_yyzz, ts_xxyy_yzz, ts_xxyy_yzzz, ts_xxyy_zzz, ts_xxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyy_xxxx[i] = 2.0 * ts_xxyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxy[i] = 2.0 * ts_xxyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxxz[i] = 2.0 * ts_xxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxyy[i] = 2.0 * ts_xxyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxyz[i] = 2.0 * ts_xxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xxzz[i] = 4.0 * ts_xxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyyy[i] = 2.0 * ts_xxyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyyz[i] = 2.0 * ts_xxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xyzz[i] = 4.0 * ts_xxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xzzz[i] = 6.0 * ts_xxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyyy[i] = 2.0 * ts_xxyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyyz[i] = 2.0 * ts_xxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yyzz[i] = 4.0 * ts_xxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yzzz[i] = 6.0 * ts_xxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_zzzz[i] = 8.0 * ts_xxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 510-525 components of targeted buffer : GG

    auto gs_z_xxyz_xxxx = pbuffer.data(idx_g_gg + 510);

    auto gs_z_xxyz_xxxy = pbuffer.data(idx_g_gg + 511);

    auto gs_z_xxyz_xxxz = pbuffer.data(idx_g_gg + 512);

    auto gs_z_xxyz_xxyy = pbuffer.data(idx_g_gg + 513);

    auto gs_z_xxyz_xxyz = pbuffer.data(idx_g_gg + 514);

    auto gs_z_xxyz_xxzz = pbuffer.data(idx_g_gg + 515);

    auto gs_z_xxyz_xyyy = pbuffer.data(idx_g_gg + 516);

    auto gs_z_xxyz_xyyz = pbuffer.data(idx_g_gg + 517);

    auto gs_z_xxyz_xyzz = pbuffer.data(idx_g_gg + 518);

    auto gs_z_xxyz_xzzz = pbuffer.data(idx_g_gg + 519);

    auto gs_z_xxyz_yyyy = pbuffer.data(idx_g_gg + 520);

    auto gs_z_xxyz_yyyz = pbuffer.data(idx_g_gg + 521);

    auto gs_z_xxyz_yyzz = pbuffer.data(idx_g_gg + 522);

    auto gs_z_xxyz_yzzz = pbuffer.data(idx_g_gg + 523);

    auto gs_z_xxyz_zzzz = pbuffer.data(idx_g_gg + 524);

    #pragma omp simd aligned(gc_z, gs_z_xxyz_xxxx, gs_z_xxyz_xxxy, gs_z_xxyz_xxxz, gs_z_xxyz_xxyy, gs_z_xxyz_xxyz, gs_z_xxyz_xxzz, gs_z_xxyz_xyyy, gs_z_xxyz_xyyz, gs_z_xxyz_xyzz, gs_z_xxyz_xzzz, gs_z_xxyz_yyyy, gs_z_xxyz_yyyz, gs_z_xxyz_yyzz, gs_z_xxyz_yzzz, gs_z_xxyz_zzzz, ts_xxy_xxxx, ts_xxy_xxxy, ts_xxy_xxxz, ts_xxy_xxyy, ts_xxy_xxyz, ts_xxy_xxzz, ts_xxy_xyyy, ts_xxy_xyyz, ts_xxy_xyzz, ts_xxy_xzzz, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyzz, ts_xxy_yzzz, ts_xxy_zzzz, ts_xxyz_xxx, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxy, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxz, ts_xxyz_xxzz, ts_xxyz_xyy, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyz, ts_xxyz_xyzz, ts_xxyz_xzz, ts_xxyz_xzzz, ts_xxyz_yyy, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyz, ts_xxyz_yyzz, ts_xxyz_yzz, ts_xxyz_yzzz, ts_xxyz_zzz, ts_xxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyz_xxxx[i] = 2.0 * ts_xxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxy[i] = 2.0 * ts_xxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxxz[i] = 2.0 * ts_xxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxyy[i] = 2.0 * ts_xxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxyz[i] = 2.0 * ts_xxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xxzz[i] = 2.0 * ts_xxy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyyy[i] = 2.0 * ts_xxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyyz[i] = 2.0 * ts_xxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xyzz[i] = 2.0 * ts_xxy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xzzz[i] = 2.0 * ts_xxy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyyy[i] = 2.0 * ts_xxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyyz[i] = 2.0 * ts_xxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yyzz[i] = 2.0 * ts_xxy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yzzz[i] = 2.0 * ts_xxy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_zzzz[i] = 2.0 * ts_xxy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 525-540 components of targeted buffer : GG

    auto gs_z_xxzz_xxxx = pbuffer.data(idx_g_gg + 525);

    auto gs_z_xxzz_xxxy = pbuffer.data(idx_g_gg + 526);

    auto gs_z_xxzz_xxxz = pbuffer.data(idx_g_gg + 527);

    auto gs_z_xxzz_xxyy = pbuffer.data(idx_g_gg + 528);

    auto gs_z_xxzz_xxyz = pbuffer.data(idx_g_gg + 529);

    auto gs_z_xxzz_xxzz = pbuffer.data(idx_g_gg + 530);

    auto gs_z_xxzz_xyyy = pbuffer.data(idx_g_gg + 531);

    auto gs_z_xxzz_xyyz = pbuffer.data(idx_g_gg + 532);

    auto gs_z_xxzz_xyzz = pbuffer.data(idx_g_gg + 533);

    auto gs_z_xxzz_xzzz = pbuffer.data(idx_g_gg + 534);

    auto gs_z_xxzz_yyyy = pbuffer.data(idx_g_gg + 535);

    auto gs_z_xxzz_yyyz = pbuffer.data(idx_g_gg + 536);

    auto gs_z_xxzz_yyzz = pbuffer.data(idx_g_gg + 537);

    auto gs_z_xxzz_yzzz = pbuffer.data(idx_g_gg + 538);

    auto gs_z_xxzz_zzzz = pbuffer.data(idx_g_gg + 539);

    #pragma omp simd aligned(gc_z, gs_z_xxzz_xxxx, gs_z_xxzz_xxxy, gs_z_xxzz_xxxz, gs_z_xxzz_xxyy, gs_z_xxzz_xxyz, gs_z_xxzz_xxzz, gs_z_xxzz_xyyy, gs_z_xxzz_xyyz, gs_z_xxzz_xyzz, gs_z_xxzz_xzzz, gs_z_xxzz_yyyy, gs_z_xxzz_yyyz, gs_z_xxzz_yyzz, gs_z_xxzz_yzzz, gs_z_xxzz_zzzz, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxxz, ts_xxz_xxyy, ts_xxz_xxyz, ts_xxz_xxzz, ts_xxz_xyyy, ts_xxz_xyyz, ts_xxz_xyzz, ts_xxz_xzzz, ts_xxz_yyyy, ts_xxz_yyyz, ts_xxz_yyzz, ts_xxz_yzzz, ts_xxz_zzzz, ts_xxzz_xxx, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxy, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxz, ts_xxzz_xxzz, ts_xxzz_xyy, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyz, ts_xxzz_xyzz, ts_xxzz_xzz, ts_xxzz_xzzz, ts_xxzz_yyy, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyz, ts_xxzz_yyzz, ts_xxzz_yzz, ts_xxzz_yzzz, ts_xxzz_zzz, ts_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzz_xxxx[i] = 4.0 * ts_xxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxy[i] = 4.0 * ts_xxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxxz[i] = 4.0 * ts_xxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxyy[i] = 4.0 * ts_xxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxyz[i] = 4.0 * ts_xxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xxzz[i] = 4.0 * ts_xxz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyyy[i] = 4.0 * ts_xxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyyz[i] = 4.0 * ts_xxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xyzz[i] = 4.0 * ts_xxz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xzzz[i] = 4.0 * ts_xxz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyyy[i] = 4.0 * ts_xxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyyz[i] = 4.0 * ts_xxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yyzz[i] = 4.0 * ts_xxz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yzzz[i] = 4.0 * ts_xxz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_zzzz[i] = 4.0 * ts_xxz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 540-555 components of targeted buffer : GG

    auto gs_z_xyyy_xxxx = pbuffer.data(idx_g_gg + 540);

    auto gs_z_xyyy_xxxy = pbuffer.data(idx_g_gg + 541);

    auto gs_z_xyyy_xxxz = pbuffer.data(idx_g_gg + 542);

    auto gs_z_xyyy_xxyy = pbuffer.data(idx_g_gg + 543);

    auto gs_z_xyyy_xxyz = pbuffer.data(idx_g_gg + 544);

    auto gs_z_xyyy_xxzz = pbuffer.data(idx_g_gg + 545);

    auto gs_z_xyyy_xyyy = pbuffer.data(idx_g_gg + 546);

    auto gs_z_xyyy_xyyz = pbuffer.data(idx_g_gg + 547);

    auto gs_z_xyyy_xyzz = pbuffer.data(idx_g_gg + 548);

    auto gs_z_xyyy_xzzz = pbuffer.data(idx_g_gg + 549);

    auto gs_z_xyyy_yyyy = pbuffer.data(idx_g_gg + 550);

    auto gs_z_xyyy_yyyz = pbuffer.data(idx_g_gg + 551);

    auto gs_z_xyyy_yyzz = pbuffer.data(idx_g_gg + 552);

    auto gs_z_xyyy_yzzz = pbuffer.data(idx_g_gg + 553);

    auto gs_z_xyyy_zzzz = pbuffer.data(idx_g_gg + 554);

    #pragma omp simd aligned(gc_z, gs_z_xyyy_xxxx, gs_z_xyyy_xxxy, gs_z_xyyy_xxxz, gs_z_xyyy_xxyy, gs_z_xyyy_xxyz, gs_z_xyyy_xxzz, gs_z_xyyy_xyyy, gs_z_xyyy_xyyz, gs_z_xyyy_xyzz, gs_z_xyyy_xzzz, gs_z_xyyy_yyyy, gs_z_xyyy_yyyz, gs_z_xyyy_yyzz, gs_z_xyyy_yzzz, gs_z_xyyy_zzzz, ts_xyyy_xxx, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxy, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxz, ts_xyyy_xxzz, ts_xyyy_xyy, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyz, ts_xyyy_xyzz, ts_xyyy_xzz, ts_xyyy_xzzz, ts_xyyy_yyy, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyz, ts_xyyy_yyzz, ts_xyyy_yzz, ts_xyyy_yzzz, ts_xyyy_zzz, ts_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyy_xxxx[i] = 2.0 * ts_xyyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxy[i] = 2.0 * ts_xyyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxxz[i] = 2.0 * ts_xyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxyy[i] = 2.0 * ts_xyyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxyz[i] = 2.0 * ts_xyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xxzz[i] = 4.0 * ts_xyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyyy[i] = 2.0 * ts_xyyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyyz[i] = 2.0 * ts_xyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xyzz[i] = 4.0 * ts_xyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xzzz[i] = 6.0 * ts_xyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyyy[i] = 2.0 * ts_xyyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyyz[i] = 2.0 * ts_xyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yyzz[i] = 4.0 * ts_xyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yzzz[i] = 6.0 * ts_xyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_zzzz[i] = 8.0 * ts_xyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 555-570 components of targeted buffer : GG

    auto gs_z_xyyz_xxxx = pbuffer.data(idx_g_gg + 555);

    auto gs_z_xyyz_xxxy = pbuffer.data(idx_g_gg + 556);

    auto gs_z_xyyz_xxxz = pbuffer.data(idx_g_gg + 557);

    auto gs_z_xyyz_xxyy = pbuffer.data(idx_g_gg + 558);

    auto gs_z_xyyz_xxyz = pbuffer.data(idx_g_gg + 559);

    auto gs_z_xyyz_xxzz = pbuffer.data(idx_g_gg + 560);

    auto gs_z_xyyz_xyyy = pbuffer.data(idx_g_gg + 561);

    auto gs_z_xyyz_xyyz = pbuffer.data(idx_g_gg + 562);

    auto gs_z_xyyz_xyzz = pbuffer.data(idx_g_gg + 563);

    auto gs_z_xyyz_xzzz = pbuffer.data(idx_g_gg + 564);

    auto gs_z_xyyz_yyyy = pbuffer.data(idx_g_gg + 565);

    auto gs_z_xyyz_yyyz = pbuffer.data(idx_g_gg + 566);

    auto gs_z_xyyz_yyzz = pbuffer.data(idx_g_gg + 567);

    auto gs_z_xyyz_yzzz = pbuffer.data(idx_g_gg + 568);

    auto gs_z_xyyz_zzzz = pbuffer.data(idx_g_gg + 569);

    #pragma omp simd aligned(gc_z, gs_z_xyyz_xxxx, gs_z_xyyz_xxxy, gs_z_xyyz_xxxz, gs_z_xyyz_xxyy, gs_z_xyyz_xxyz, gs_z_xyyz_xxzz, gs_z_xyyz_xyyy, gs_z_xyyz_xyyz, gs_z_xyyz_xyzz, gs_z_xyyz_xzzz, gs_z_xyyz_yyyy, gs_z_xyyz_yyyz, gs_z_xyyz_yyzz, gs_z_xyyz_yzzz, gs_z_xyyz_zzzz, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxxz, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xxzz, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyzz, ts_xyy_xzzz, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyzz, ts_xyy_yzzz, ts_xyy_zzzz, ts_xyyz_xxx, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxy, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxz, ts_xyyz_xxzz, ts_xyyz_xyy, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyz, ts_xyyz_xyzz, ts_xyyz_xzz, ts_xyyz_xzzz, ts_xyyz_yyy, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyz, ts_xyyz_yyzz, ts_xyyz_yzz, ts_xyyz_yzzz, ts_xyyz_zzz, ts_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyz_xxxx[i] = 2.0 * ts_xyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxy[i] = 2.0 * ts_xyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxxz[i] = 2.0 * ts_xyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxyy[i] = 2.0 * ts_xyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxyz[i] = 2.0 * ts_xyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xxzz[i] = 2.0 * ts_xyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyyy[i] = 2.0 * ts_xyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyyz[i] = 2.0 * ts_xyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xyzz[i] = 2.0 * ts_xyy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xzzz[i] = 2.0 * ts_xyy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyyy[i] = 2.0 * ts_xyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyyz[i] = 2.0 * ts_xyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yyzz[i] = 2.0 * ts_xyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yzzz[i] = 2.0 * ts_xyy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_zzzz[i] = 2.0 * ts_xyy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 570-585 components of targeted buffer : GG

    auto gs_z_xyzz_xxxx = pbuffer.data(idx_g_gg + 570);

    auto gs_z_xyzz_xxxy = pbuffer.data(idx_g_gg + 571);

    auto gs_z_xyzz_xxxz = pbuffer.data(idx_g_gg + 572);

    auto gs_z_xyzz_xxyy = pbuffer.data(idx_g_gg + 573);

    auto gs_z_xyzz_xxyz = pbuffer.data(idx_g_gg + 574);

    auto gs_z_xyzz_xxzz = pbuffer.data(idx_g_gg + 575);

    auto gs_z_xyzz_xyyy = pbuffer.data(idx_g_gg + 576);

    auto gs_z_xyzz_xyyz = pbuffer.data(idx_g_gg + 577);

    auto gs_z_xyzz_xyzz = pbuffer.data(idx_g_gg + 578);

    auto gs_z_xyzz_xzzz = pbuffer.data(idx_g_gg + 579);

    auto gs_z_xyzz_yyyy = pbuffer.data(idx_g_gg + 580);

    auto gs_z_xyzz_yyyz = pbuffer.data(idx_g_gg + 581);

    auto gs_z_xyzz_yyzz = pbuffer.data(idx_g_gg + 582);

    auto gs_z_xyzz_yzzz = pbuffer.data(idx_g_gg + 583);

    auto gs_z_xyzz_zzzz = pbuffer.data(idx_g_gg + 584);

    #pragma omp simd aligned(gc_z, gs_z_xyzz_xxxx, gs_z_xyzz_xxxy, gs_z_xyzz_xxxz, gs_z_xyzz_xxyy, gs_z_xyzz_xxyz, gs_z_xyzz_xxzz, gs_z_xyzz_xyyy, gs_z_xyzz_xyyz, gs_z_xyzz_xyzz, gs_z_xyzz_xzzz, gs_z_xyzz_yyyy, gs_z_xyzz_yyyz, gs_z_xyzz_yyzz, gs_z_xyzz_yzzz, gs_z_xyzz_zzzz, ts_xyz_xxxx, ts_xyz_xxxy, ts_xyz_xxxz, ts_xyz_xxyy, ts_xyz_xxyz, ts_xyz_xxzz, ts_xyz_xyyy, ts_xyz_xyyz, ts_xyz_xyzz, ts_xyz_xzzz, ts_xyz_yyyy, ts_xyz_yyyz, ts_xyz_yyzz, ts_xyz_yzzz, ts_xyz_zzzz, ts_xyzz_xxx, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxy, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxz, ts_xyzz_xxzz, ts_xyzz_xyy, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyz, ts_xyzz_xyzz, ts_xyzz_xzz, ts_xyzz_xzzz, ts_xyzz_yyy, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyz, ts_xyzz_yyzz, ts_xyzz_yzz, ts_xyzz_yzzz, ts_xyzz_zzz, ts_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzz_xxxx[i] = 4.0 * ts_xyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxy[i] = 4.0 * ts_xyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxxz[i] = 4.0 * ts_xyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxyy[i] = 4.0 * ts_xyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxyz[i] = 4.0 * ts_xyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xxzz[i] = 4.0 * ts_xyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyyy[i] = 4.0 * ts_xyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyyz[i] = 4.0 * ts_xyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xyzz[i] = 4.0 * ts_xyz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xzzz[i] = 4.0 * ts_xyz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyyy[i] = 4.0 * ts_xyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyyz[i] = 4.0 * ts_xyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yyzz[i] = 4.0 * ts_xyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yzzz[i] = 4.0 * ts_xyz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_zzzz[i] = 4.0 * ts_xyz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 585-600 components of targeted buffer : GG

    auto gs_z_xzzz_xxxx = pbuffer.data(idx_g_gg + 585);

    auto gs_z_xzzz_xxxy = pbuffer.data(idx_g_gg + 586);

    auto gs_z_xzzz_xxxz = pbuffer.data(idx_g_gg + 587);

    auto gs_z_xzzz_xxyy = pbuffer.data(idx_g_gg + 588);

    auto gs_z_xzzz_xxyz = pbuffer.data(idx_g_gg + 589);

    auto gs_z_xzzz_xxzz = pbuffer.data(idx_g_gg + 590);

    auto gs_z_xzzz_xyyy = pbuffer.data(idx_g_gg + 591);

    auto gs_z_xzzz_xyyz = pbuffer.data(idx_g_gg + 592);

    auto gs_z_xzzz_xyzz = pbuffer.data(idx_g_gg + 593);

    auto gs_z_xzzz_xzzz = pbuffer.data(idx_g_gg + 594);

    auto gs_z_xzzz_yyyy = pbuffer.data(idx_g_gg + 595);

    auto gs_z_xzzz_yyyz = pbuffer.data(idx_g_gg + 596);

    auto gs_z_xzzz_yyzz = pbuffer.data(idx_g_gg + 597);

    auto gs_z_xzzz_yzzz = pbuffer.data(idx_g_gg + 598);

    auto gs_z_xzzz_zzzz = pbuffer.data(idx_g_gg + 599);

    #pragma omp simd aligned(gc_z, gs_z_xzzz_xxxx, gs_z_xzzz_xxxy, gs_z_xzzz_xxxz, gs_z_xzzz_xxyy, gs_z_xzzz_xxyz, gs_z_xzzz_xxzz, gs_z_xzzz_xyyy, gs_z_xzzz_xyyz, gs_z_xzzz_xyzz, gs_z_xzzz_xzzz, gs_z_xzzz_yyyy, gs_z_xzzz_yyyz, gs_z_xzzz_yyzz, gs_z_xzzz_yzzz, gs_z_xzzz_zzzz, ts_xzz_xxxx, ts_xzz_xxxy, ts_xzz_xxxz, ts_xzz_xxyy, ts_xzz_xxyz, ts_xzz_xxzz, ts_xzz_xyyy, ts_xzz_xyyz, ts_xzz_xyzz, ts_xzz_xzzz, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyzz, ts_xzz_yzzz, ts_xzz_zzzz, ts_xzzz_xxx, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxy, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxz, ts_xzzz_xxzz, ts_xzzz_xyy, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyz, ts_xzzz_xyzz, ts_xzzz_xzz, ts_xzzz_xzzz, ts_xzzz_yyy, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyz, ts_xzzz_yyzz, ts_xzzz_yzz, ts_xzzz_yzzz, ts_xzzz_zzz, ts_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzz_xxxx[i] = 6.0 * ts_xzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxy[i] = 6.0 * ts_xzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxxz[i] = 6.0 * ts_xzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxyy[i] = 6.0 * ts_xzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxyz[i] = 6.0 * ts_xzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xxzz[i] = 6.0 * ts_xzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyyy[i] = 6.0 * ts_xzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyyz[i] = 6.0 * ts_xzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xyzz[i] = 6.0 * ts_xzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xzzz[i] = 6.0 * ts_xzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyyy[i] = 6.0 * ts_xzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyyz[i] = 6.0 * ts_xzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yyzz[i] = 6.0 * ts_xzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yzzz[i] = 6.0 * ts_xzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_zzzz[i] = 6.0 * ts_xzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 600-615 components of targeted buffer : GG

    auto gs_z_yyyy_xxxx = pbuffer.data(idx_g_gg + 600);

    auto gs_z_yyyy_xxxy = pbuffer.data(idx_g_gg + 601);

    auto gs_z_yyyy_xxxz = pbuffer.data(idx_g_gg + 602);

    auto gs_z_yyyy_xxyy = pbuffer.data(idx_g_gg + 603);

    auto gs_z_yyyy_xxyz = pbuffer.data(idx_g_gg + 604);

    auto gs_z_yyyy_xxzz = pbuffer.data(idx_g_gg + 605);

    auto gs_z_yyyy_xyyy = pbuffer.data(idx_g_gg + 606);

    auto gs_z_yyyy_xyyz = pbuffer.data(idx_g_gg + 607);

    auto gs_z_yyyy_xyzz = pbuffer.data(idx_g_gg + 608);

    auto gs_z_yyyy_xzzz = pbuffer.data(idx_g_gg + 609);

    auto gs_z_yyyy_yyyy = pbuffer.data(idx_g_gg + 610);

    auto gs_z_yyyy_yyyz = pbuffer.data(idx_g_gg + 611);

    auto gs_z_yyyy_yyzz = pbuffer.data(idx_g_gg + 612);

    auto gs_z_yyyy_yzzz = pbuffer.data(idx_g_gg + 613);

    auto gs_z_yyyy_zzzz = pbuffer.data(idx_g_gg + 614);

    #pragma omp simd aligned(gc_z, gs_z_yyyy_xxxx, gs_z_yyyy_xxxy, gs_z_yyyy_xxxz, gs_z_yyyy_xxyy, gs_z_yyyy_xxyz, gs_z_yyyy_xxzz, gs_z_yyyy_xyyy, gs_z_yyyy_xyyz, gs_z_yyyy_xyzz, gs_z_yyyy_xzzz, gs_z_yyyy_yyyy, gs_z_yyyy_yyyz, gs_z_yyyy_yyzz, gs_z_yyyy_yzzz, gs_z_yyyy_zzzz, ts_yyyy_xxx, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxz, ts_yyyy_xxzz, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_xzz, ts_yyyy_xzzz, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyy_zzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyy_xxxx[i] = 2.0 * ts_yyyy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxy[i] = 2.0 * ts_yyyy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxxz[i] = 2.0 * ts_yyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxyy[i] = 2.0 * ts_yyyy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxyz[i] = 2.0 * ts_yyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xxzz[i] = 4.0 * ts_yyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyyy[i] = 2.0 * ts_yyyy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyyz[i] = 2.0 * ts_yyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xyzz[i] = 4.0 * ts_yyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xzzz[i] = 6.0 * ts_yyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyyy[i] = 2.0 * ts_yyyy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyyz[i] = 2.0 * ts_yyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yyzz[i] = 4.0 * ts_yyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yzzz[i] = 6.0 * ts_yyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_zzzz[i] = 8.0 * ts_yyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 615-630 components of targeted buffer : GG

    auto gs_z_yyyz_xxxx = pbuffer.data(idx_g_gg + 615);

    auto gs_z_yyyz_xxxy = pbuffer.data(idx_g_gg + 616);

    auto gs_z_yyyz_xxxz = pbuffer.data(idx_g_gg + 617);

    auto gs_z_yyyz_xxyy = pbuffer.data(idx_g_gg + 618);

    auto gs_z_yyyz_xxyz = pbuffer.data(idx_g_gg + 619);

    auto gs_z_yyyz_xxzz = pbuffer.data(idx_g_gg + 620);

    auto gs_z_yyyz_xyyy = pbuffer.data(idx_g_gg + 621);

    auto gs_z_yyyz_xyyz = pbuffer.data(idx_g_gg + 622);

    auto gs_z_yyyz_xyzz = pbuffer.data(idx_g_gg + 623);

    auto gs_z_yyyz_xzzz = pbuffer.data(idx_g_gg + 624);

    auto gs_z_yyyz_yyyy = pbuffer.data(idx_g_gg + 625);

    auto gs_z_yyyz_yyyz = pbuffer.data(idx_g_gg + 626);

    auto gs_z_yyyz_yyzz = pbuffer.data(idx_g_gg + 627);

    auto gs_z_yyyz_yzzz = pbuffer.data(idx_g_gg + 628);

    auto gs_z_yyyz_zzzz = pbuffer.data(idx_g_gg + 629);

    #pragma omp simd aligned(gc_z, gs_z_yyyz_xxxx, gs_z_yyyz_xxxy, gs_z_yyyz_xxxz, gs_z_yyyz_xxyy, gs_z_yyyz_xxyz, gs_z_yyyz_xxzz, gs_z_yyyz_xyyy, gs_z_yyyz_xyyz, gs_z_yyyz_xyzz, gs_z_yyyz_xzzz, gs_z_yyyz_yyyy, gs_z_yyyz_yyyz, gs_z_yyyz_yyzz, gs_z_yyyz_yzzz, gs_z_yyyz_zzzz, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxzz, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyzz, ts_yyy_xzzz, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyzz, ts_yyy_yzzz, ts_yyy_zzzz, ts_yyyz_xxx, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxy, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxz, ts_yyyz_xxzz, ts_yyyz_xyy, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyz, ts_yyyz_xyzz, ts_yyyz_xzz, ts_yyyz_xzzz, ts_yyyz_yyy, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyz, ts_yyyz_yyzz, ts_yyyz_yzz, ts_yyyz_yzzz, ts_yyyz_zzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyz_xxxx[i] = 2.0 * ts_yyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxy[i] = 2.0 * ts_yyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxxz[i] = 2.0 * ts_yyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxyy[i] = 2.0 * ts_yyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxyz[i] = 2.0 * ts_yyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xxzz[i] = 2.0 * ts_yyy_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyyy[i] = 2.0 * ts_yyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyyz[i] = 2.0 * ts_yyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xyzz[i] = 2.0 * ts_yyy_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xzzz[i] = 2.0 * ts_yyy_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyyy[i] = 2.0 * ts_yyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyyz[i] = 2.0 * ts_yyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yyzz[i] = 2.0 * ts_yyy_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yzzz[i] = 2.0 * ts_yyy_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_zzzz[i] = 2.0 * ts_yyy_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 630-645 components of targeted buffer : GG

    auto gs_z_yyzz_xxxx = pbuffer.data(idx_g_gg + 630);

    auto gs_z_yyzz_xxxy = pbuffer.data(idx_g_gg + 631);

    auto gs_z_yyzz_xxxz = pbuffer.data(idx_g_gg + 632);

    auto gs_z_yyzz_xxyy = pbuffer.data(idx_g_gg + 633);

    auto gs_z_yyzz_xxyz = pbuffer.data(idx_g_gg + 634);

    auto gs_z_yyzz_xxzz = pbuffer.data(idx_g_gg + 635);

    auto gs_z_yyzz_xyyy = pbuffer.data(idx_g_gg + 636);

    auto gs_z_yyzz_xyyz = pbuffer.data(idx_g_gg + 637);

    auto gs_z_yyzz_xyzz = pbuffer.data(idx_g_gg + 638);

    auto gs_z_yyzz_xzzz = pbuffer.data(idx_g_gg + 639);

    auto gs_z_yyzz_yyyy = pbuffer.data(idx_g_gg + 640);

    auto gs_z_yyzz_yyyz = pbuffer.data(idx_g_gg + 641);

    auto gs_z_yyzz_yyzz = pbuffer.data(idx_g_gg + 642);

    auto gs_z_yyzz_yzzz = pbuffer.data(idx_g_gg + 643);

    auto gs_z_yyzz_zzzz = pbuffer.data(idx_g_gg + 644);

    #pragma omp simd aligned(gc_z, gs_z_yyzz_xxxx, gs_z_yyzz_xxxy, gs_z_yyzz_xxxz, gs_z_yyzz_xxyy, gs_z_yyzz_xxyz, gs_z_yyzz_xxzz, gs_z_yyzz_xyyy, gs_z_yyzz_xyyz, gs_z_yyzz_xyzz, gs_z_yyzz_xzzz, gs_z_yyzz_yyyy, gs_z_yyzz_yyyz, gs_z_yyzz_yyzz, gs_z_yyzz_yzzz, gs_z_yyzz_zzzz, ts_yyz_xxxx, ts_yyz_xxxy, ts_yyz_xxxz, ts_yyz_xxyy, ts_yyz_xxyz, ts_yyz_xxzz, ts_yyz_xyyy, ts_yyz_xyyz, ts_yyz_xyzz, ts_yyz_xzzz, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyzz, ts_yyz_yzzz, ts_yyz_zzzz, ts_yyzz_xxx, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxy, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxz, ts_yyzz_xxzz, ts_yyzz_xyy, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyz, ts_yyzz_xyzz, ts_yyzz_xzz, ts_yyzz_xzzz, ts_yyzz_yyy, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyz, ts_yyzz_yyzz, ts_yyzz_yzz, ts_yyzz_yzzz, ts_yyzz_zzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzz_xxxx[i] = 4.0 * ts_yyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxy[i] = 4.0 * ts_yyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxxz[i] = 4.0 * ts_yyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxyy[i] = 4.0 * ts_yyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxyz[i] = 4.0 * ts_yyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xxzz[i] = 4.0 * ts_yyz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyyy[i] = 4.0 * ts_yyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyyz[i] = 4.0 * ts_yyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xyzz[i] = 4.0 * ts_yyz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xzzz[i] = 4.0 * ts_yyz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyyy[i] = 4.0 * ts_yyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyyz[i] = 4.0 * ts_yyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yyzz[i] = 4.0 * ts_yyz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yzzz[i] = 4.0 * ts_yyz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_zzzz[i] = 4.0 * ts_yyz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 645-660 components of targeted buffer : GG

    auto gs_z_yzzz_xxxx = pbuffer.data(idx_g_gg + 645);

    auto gs_z_yzzz_xxxy = pbuffer.data(idx_g_gg + 646);

    auto gs_z_yzzz_xxxz = pbuffer.data(idx_g_gg + 647);

    auto gs_z_yzzz_xxyy = pbuffer.data(idx_g_gg + 648);

    auto gs_z_yzzz_xxyz = pbuffer.data(idx_g_gg + 649);

    auto gs_z_yzzz_xxzz = pbuffer.data(idx_g_gg + 650);

    auto gs_z_yzzz_xyyy = pbuffer.data(idx_g_gg + 651);

    auto gs_z_yzzz_xyyz = pbuffer.data(idx_g_gg + 652);

    auto gs_z_yzzz_xyzz = pbuffer.data(idx_g_gg + 653);

    auto gs_z_yzzz_xzzz = pbuffer.data(idx_g_gg + 654);

    auto gs_z_yzzz_yyyy = pbuffer.data(idx_g_gg + 655);

    auto gs_z_yzzz_yyyz = pbuffer.data(idx_g_gg + 656);

    auto gs_z_yzzz_yyzz = pbuffer.data(idx_g_gg + 657);

    auto gs_z_yzzz_yzzz = pbuffer.data(idx_g_gg + 658);

    auto gs_z_yzzz_zzzz = pbuffer.data(idx_g_gg + 659);

    #pragma omp simd aligned(gc_z, gs_z_yzzz_xxxx, gs_z_yzzz_xxxy, gs_z_yzzz_xxxz, gs_z_yzzz_xxyy, gs_z_yzzz_xxyz, gs_z_yzzz_xxzz, gs_z_yzzz_xyyy, gs_z_yzzz_xyyz, gs_z_yzzz_xyzz, gs_z_yzzz_xzzz, gs_z_yzzz_yyyy, gs_z_yzzz_yyyz, gs_z_yzzz_yyzz, gs_z_yzzz_yzzz, gs_z_yzzz_zzzz, ts_yzz_xxxx, ts_yzz_xxxy, ts_yzz_xxxz, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xxzz, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyzz, ts_yzz_xzzz, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyzz, ts_yzz_yzzz, ts_yzz_zzzz, ts_yzzz_xxx, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxy, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxz, ts_yzzz_xxzz, ts_yzzz_xyy, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyz, ts_yzzz_xyzz, ts_yzzz_xzz, ts_yzzz_xzzz, ts_yzzz_yyy, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyz, ts_yzzz_yyzz, ts_yzzz_yzz, ts_yzzz_yzzz, ts_yzzz_zzz, ts_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzz_xxxx[i] = 6.0 * ts_yzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxy[i] = 6.0 * ts_yzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxxz[i] = 6.0 * ts_yzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxyy[i] = 6.0 * ts_yzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxyz[i] = 6.0 * ts_yzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xxzz[i] = 6.0 * ts_yzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyyy[i] = 6.0 * ts_yzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyyz[i] = 6.0 * ts_yzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xyzz[i] = 6.0 * ts_yzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xzzz[i] = 6.0 * ts_yzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyyy[i] = 6.0 * ts_yzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyyz[i] = 6.0 * ts_yzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yyzz[i] = 6.0 * ts_yzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yzzz[i] = 6.0 * ts_yzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_zzzz[i] = 6.0 * ts_yzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 660-675 components of targeted buffer : GG

    auto gs_z_zzzz_xxxx = pbuffer.data(idx_g_gg + 660);

    auto gs_z_zzzz_xxxy = pbuffer.data(idx_g_gg + 661);

    auto gs_z_zzzz_xxxz = pbuffer.data(idx_g_gg + 662);

    auto gs_z_zzzz_xxyy = pbuffer.data(idx_g_gg + 663);

    auto gs_z_zzzz_xxyz = pbuffer.data(idx_g_gg + 664);

    auto gs_z_zzzz_xxzz = pbuffer.data(idx_g_gg + 665);

    auto gs_z_zzzz_xyyy = pbuffer.data(idx_g_gg + 666);

    auto gs_z_zzzz_xyyz = pbuffer.data(idx_g_gg + 667);

    auto gs_z_zzzz_xyzz = pbuffer.data(idx_g_gg + 668);

    auto gs_z_zzzz_xzzz = pbuffer.data(idx_g_gg + 669);

    auto gs_z_zzzz_yyyy = pbuffer.data(idx_g_gg + 670);

    auto gs_z_zzzz_yyyz = pbuffer.data(idx_g_gg + 671);

    auto gs_z_zzzz_yyzz = pbuffer.data(idx_g_gg + 672);

    auto gs_z_zzzz_yzzz = pbuffer.data(idx_g_gg + 673);

    auto gs_z_zzzz_zzzz = pbuffer.data(idx_g_gg + 674);

    #pragma omp simd aligned(gc_z, gs_z_zzzz_xxxx, gs_z_zzzz_xxxy, gs_z_zzzz_xxxz, gs_z_zzzz_xxyy, gs_z_zzzz_xxyz, gs_z_zzzz_xxzz, gs_z_zzzz_xyyy, gs_z_zzzz_xyyz, gs_z_zzzz_xyzz, gs_z_zzzz_xzzz, gs_z_zzzz_yyyy, gs_z_zzzz_yyyz, gs_z_zzzz_yyzz, gs_z_zzzz_yzzz, gs_z_zzzz_zzzz, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxzz, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyzz, ts_zzz_xzzz, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyzz, ts_zzz_yzzz, ts_zzz_zzzz, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzz_xxxx[i] = 8.0 * ts_zzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxy[i] = 8.0 * ts_zzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxxz[i] = 8.0 * ts_zzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxyy[i] = 8.0 * ts_zzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxyz[i] = 8.0 * ts_zzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xxzz[i] = 8.0 * ts_zzz_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyyy[i] = 8.0 * ts_zzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyyz[i] = 8.0 * ts_zzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xyzz[i] = 8.0 * ts_zzz_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xzzz[i] = 8.0 * ts_zzz_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyyy[i] = 8.0 * ts_zzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyyz[i] = 8.0 * ts_zzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yyzz[i] = 8.0 * ts_zzz_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yzzz[i] = 8.0 * ts_zzz_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_zzzz[i] = 8.0 * ts_zzz_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_zzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

