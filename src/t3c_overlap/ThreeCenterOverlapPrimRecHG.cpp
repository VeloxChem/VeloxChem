#include "ThreeCenterOverlapPrimRecHG.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_hg(CSimdArray<double>& pbuffer, 
                     const size_t idx_hg,
                     const size_t idx_fg,
                     const size_t idx_gf,
                     const size_t idx_gg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

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

    auto ts_xxy_xxxz = pbuffer.data(idx_fg + 17);

    auto ts_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto ts_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto ts_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto ts_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto ts_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto ts_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto ts_xxz_xxxx = pbuffer.data(idx_fg + 30);

    auto ts_xxz_xxxy = pbuffer.data(idx_fg + 31);

    auto ts_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto ts_xxz_xxyy = pbuffer.data(idx_fg + 33);

    auto ts_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto ts_xxz_xyyy = pbuffer.data(idx_fg + 36);

    auto ts_xxz_xzzz = pbuffer.data(idx_fg + 39);

    auto ts_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto ts_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto ts_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto ts_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto ts_xyy_xxxy = pbuffer.data(idx_fg + 46);

    auto ts_xyy_xxyy = pbuffer.data(idx_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_fg + 49);

    auto ts_xyy_xyyy = pbuffer.data(idx_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_fg + 53);

    auto ts_xyy_yyyy = pbuffer.data(idx_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_fg + 59);

    auto ts_xyz_yyyz = pbuffer.data(idx_fg + 71);

    auto ts_xyz_yyzz = pbuffer.data(idx_fg + 72);

    auto ts_xyz_yzzz = pbuffer.data(idx_fg + 73);

    auto ts_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto ts_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_fg + 80);

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

    auto ts_yyz_xxxy = pbuffer.data(idx_fg + 106);

    auto ts_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto ts_yyz_xxyy = pbuffer.data(idx_fg + 108);

    auto ts_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto ts_yyz_xyyy = pbuffer.data(idx_fg + 111);

    auto ts_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto ts_yyz_yyyy = pbuffer.data(idx_fg + 115);

    auto ts_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_fg + 119);

    auto ts_yzz_xxxx = pbuffer.data(idx_fg + 120);

    auto ts_yzz_xxxz = pbuffer.data(idx_fg + 122);

    auto ts_yzz_xxyz = pbuffer.data(idx_fg + 124);

    auto ts_yzz_xxzz = pbuffer.data(idx_fg + 125);

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

    auto ts_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto ts_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_gf + 25);

    auto ts_xxyy_xxy = pbuffer.data(idx_gf + 31);

    auto ts_xxyy_xyy = pbuffer.data(idx_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_gf + 34);

    auto ts_xxyy_yyy = pbuffer.data(idx_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_gf + 38);

    auto ts_xxzz_xxx = pbuffer.data(idx_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_gf + 55);

    auto ts_xxzz_yyz = pbuffer.data(idx_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_gf + 59);

    auto ts_xyyy_xxy = pbuffer.data(idx_gf + 61);

    auto ts_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto ts_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto ts_xzzz_xxz = pbuffer.data(idx_gf + 92);

    auto ts_xzzz_xyz = pbuffer.data(idx_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_gf + 95);

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

    auto ts_yyyz_xxz = pbuffer.data(idx_gf + 112);

    auto ts_yyyz_xyz = pbuffer.data(idx_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_gf + 115);

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

    auto ts_xxxy_xxzz = pbuffer.data(idx_gg + 20);

    auto ts_xxxy_xyyy = pbuffer.data(idx_gg + 21);

    auto ts_xxxy_xzzz = pbuffer.data(idx_gg + 24);

    auto ts_xxxy_yyyy = pbuffer.data(idx_gg + 25);

    auto ts_xxxy_yyyz = pbuffer.data(idx_gg + 26);

    auto ts_xxxy_yyzz = pbuffer.data(idx_gg + 27);

    auto ts_xxxy_yzzz = pbuffer.data(idx_gg + 28);

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

    auto ts_xxyz_xxxz = pbuffer.data(idx_gg + 62);

    auto ts_xxyz_xxzz = pbuffer.data(idx_gg + 65);

    auto ts_xxyz_xzzz = pbuffer.data(idx_gg + 69);

    auto ts_xxyz_yyyz = pbuffer.data(idx_gg + 71);

    auto ts_xxyz_yyzz = pbuffer.data(idx_gg + 72);

    auto ts_xxyz_yzzz = pbuffer.data(idx_gg + 73);

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

    auto ts_xyyy_xxyy = pbuffer.data(idx_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_gg + 94);

    auto ts_xyyy_xyyy = pbuffer.data(idx_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_gg + 98);

    auto ts_xyyy_yyyy = pbuffer.data(idx_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_gg + 104);

    auto ts_xyyz_yyyz = pbuffer.data(idx_gg + 116);

    auto ts_xyyz_yyzz = pbuffer.data(idx_gg + 117);

    auto ts_xyyz_yzzz = pbuffer.data(idx_gg + 118);

    auto ts_xyyz_zzzz = pbuffer.data(idx_gg + 119);

    auto ts_xyzz_yyyy = pbuffer.data(idx_gg + 130);

    auto ts_xyzz_yyyz = pbuffer.data(idx_gg + 131);

    auto ts_xyzz_yyzz = pbuffer.data(idx_gg + 132);

    auto ts_xyzz_yzzz = pbuffer.data(idx_gg + 133);

    auto ts_xzzz_xxxx = pbuffer.data(idx_gg + 135);

    auto ts_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto ts_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_gg + 140);

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

    // Set up 0-15 components of targeted buffer : HG

    auto ts_xxxxx_xxxx = pbuffer.data(idx_hg);

    auto ts_xxxxx_xxxy = pbuffer.data(idx_hg + 1);

    auto ts_xxxxx_xxxz = pbuffer.data(idx_hg + 2);

    auto ts_xxxxx_xxyy = pbuffer.data(idx_hg + 3);

    auto ts_xxxxx_xxyz = pbuffer.data(idx_hg + 4);

    auto ts_xxxxx_xxzz = pbuffer.data(idx_hg + 5);

    auto ts_xxxxx_xyyy = pbuffer.data(idx_hg + 6);

    auto ts_xxxxx_xyyz = pbuffer.data(idx_hg + 7);

    auto ts_xxxxx_xyzz = pbuffer.data(idx_hg + 8);

    auto ts_xxxxx_xzzz = pbuffer.data(idx_hg + 9);

    auto ts_xxxxx_yyyy = pbuffer.data(idx_hg + 10);

    auto ts_xxxxx_yyyz = pbuffer.data(idx_hg + 11);

    auto ts_xxxxx_yyzz = pbuffer.data(idx_hg + 12);

    auto ts_xxxxx_yzzz = pbuffer.data(idx_hg + 13);

    auto ts_xxxxx_zzzz = pbuffer.data(idx_hg + 14);

    #pragma omp simd aligned(ga_x, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxzz, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyzz, ts_xxx_xzzz, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyzz, ts_xxx_yzzz, ts_xxx_zzzz, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_yyy, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyz, ts_xxxx_yyzz, ts_xxxx_yzz, ts_xxxx_yzzz, ts_xxxx_zzz, ts_xxxx_zzzz, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxzz, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyzz, ts_xxxxx_xzzz, ts_xxxxx_yyyy, ts_xxxxx_yyyz, ts_xxxxx_yyzz, ts_xxxxx_yzzz, ts_xxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxx_xxxx[i] = 4.0 * ts_xxx_xxxx[i] * gfe_0 + 4.0 * ts_xxxx_xxx[i] * gfe_0 + ts_xxxx_xxxx[i] * ga_x[i];

        ts_xxxxx_xxxy[i] = 4.0 * ts_xxx_xxxy[i] * gfe_0 + 3.0 * ts_xxxx_xxy[i] * gfe_0 + ts_xxxx_xxxy[i] * ga_x[i];

        ts_xxxxx_xxxz[i] = 4.0 * ts_xxx_xxxz[i] * gfe_0 + 3.0 * ts_xxxx_xxz[i] * gfe_0 + ts_xxxx_xxxz[i] * ga_x[i];

        ts_xxxxx_xxyy[i] = 4.0 * ts_xxx_xxyy[i] * gfe_0 + 2.0 * ts_xxxx_xyy[i] * gfe_0 + ts_xxxx_xxyy[i] * ga_x[i];

        ts_xxxxx_xxyz[i] = 4.0 * ts_xxx_xxyz[i] * gfe_0 + 2.0 * ts_xxxx_xyz[i] * gfe_0 + ts_xxxx_xxyz[i] * ga_x[i];

        ts_xxxxx_xxzz[i] = 4.0 * ts_xxx_xxzz[i] * gfe_0 + 2.0 * ts_xxxx_xzz[i] * gfe_0 + ts_xxxx_xxzz[i] * ga_x[i];

        ts_xxxxx_xyyy[i] = 4.0 * ts_xxx_xyyy[i] * gfe_0 + ts_xxxx_yyy[i] * gfe_0 + ts_xxxx_xyyy[i] * ga_x[i];

        ts_xxxxx_xyyz[i] = 4.0 * ts_xxx_xyyz[i] * gfe_0 + ts_xxxx_yyz[i] * gfe_0 + ts_xxxx_xyyz[i] * ga_x[i];

        ts_xxxxx_xyzz[i] = 4.0 * ts_xxx_xyzz[i] * gfe_0 + ts_xxxx_yzz[i] * gfe_0 + ts_xxxx_xyzz[i] * ga_x[i];

        ts_xxxxx_xzzz[i] = 4.0 * ts_xxx_xzzz[i] * gfe_0 + ts_xxxx_zzz[i] * gfe_0 + ts_xxxx_xzzz[i] * ga_x[i];

        ts_xxxxx_yyyy[i] = 4.0 * ts_xxx_yyyy[i] * gfe_0 + ts_xxxx_yyyy[i] * ga_x[i];

        ts_xxxxx_yyyz[i] = 4.0 * ts_xxx_yyyz[i] * gfe_0 + ts_xxxx_yyyz[i] * ga_x[i];

        ts_xxxxx_yyzz[i] = 4.0 * ts_xxx_yyzz[i] * gfe_0 + ts_xxxx_yyzz[i] * ga_x[i];

        ts_xxxxx_yzzz[i] = 4.0 * ts_xxx_yzzz[i] * gfe_0 + ts_xxxx_yzzz[i] * ga_x[i];

        ts_xxxxx_zzzz[i] = 4.0 * ts_xxx_zzzz[i] * gfe_0 + ts_xxxx_zzzz[i] * ga_x[i];
    }

    // Set up 15-30 components of targeted buffer : HG

    auto ts_xxxxy_xxxx = pbuffer.data(idx_hg + 15);

    auto ts_xxxxy_xxxy = pbuffer.data(idx_hg + 16);

    auto ts_xxxxy_xxxz = pbuffer.data(idx_hg + 17);

    auto ts_xxxxy_xxyy = pbuffer.data(idx_hg + 18);

    auto ts_xxxxy_xxyz = pbuffer.data(idx_hg + 19);

    auto ts_xxxxy_xxzz = pbuffer.data(idx_hg + 20);

    auto ts_xxxxy_xyyy = pbuffer.data(idx_hg + 21);

    auto ts_xxxxy_xyyz = pbuffer.data(idx_hg + 22);

    auto ts_xxxxy_xyzz = pbuffer.data(idx_hg + 23);

    auto ts_xxxxy_xzzz = pbuffer.data(idx_hg + 24);

    auto ts_xxxxy_yyyy = pbuffer.data(idx_hg + 25);

    auto ts_xxxxy_yyyz = pbuffer.data(idx_hg + 26);

    auto ts_xxxxy_yyzz = pbuffer.data(idx_hg + 27);

    auto ts_xxxxy_yzzz = pbuffer.data(idx_hg + 28);

    auto ts_xxxxy_zzzz = pbuffer.data(idx_hg + 29);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_zzzz, ts_xxxxy_xxxx, ts_xxxxy_xxxy, ts_xxxxy_xxxz, ts_xxxxy_xxyy, ts_xxxxy_xxyz, ts_xxxxy_xxzz, ts_xxxxy_xyyy, ts_xxxxy_xyyz, ts_xxxxy_xyzz, ts_xxxxy_xzzz, ts_xxxxy_yyyy, ts_xxxxy_yyyz, ts_xxxxy_yyzz, ts_xxxxy_yzzz, ts_xxxxy_zzzz, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyzz, ts_xxxy_yzzz, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyzz, ts_xxy_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxy_xxxx[i] = ts_xxxx_xxxx[i] * ga_y[i];

        ts_xxxxy_xxxy[i] = ts_xxxx_xxx[i] * gfe_0 + ts_xxxx_xxxy[i] * ga_y[i];

        ts_xxxxy_xxxz[i] = ts_xxxx_xxxz[i] * ga_y[i];

        ts_xxxxy_xxyy[i] = 2.0 * ts_xxxx_xxy[i] * gfe_0 + ts_xxxx_xxyy[i] * ga_y[i];

        ts_xxxxy_xxyz[i] = ts_xxxx_xxz[i] * gfe_0 + ts_xxxx_xxyz[i] * ga_y[i];

        ts_xxxxy_xxzz[i] = ts_xxxx_xxzz[i] * ga_y[i];

        ts_xxxxy_xyyy[i] = 3.0 * ts_xxxx_xyy[i] * gfe_0 + ts_xxxx_xyyy[i] * ga_y[i];

        ts_xxxxy_xyyz[i] = 2.0 * ts_xxxx_xyz[i] * gfe_0 + ts_xxxx_xyyz[i] * ga_y[i];

        ts_xxxxy_xyzz[i] = ts_xxxx_xzz[i] * gfe_0 + ts_xxxx_xyzz[i] * ga_y[i];

        ts_xxxxy_xzzz[i] = ts_xxxx_xzzz[i] * ga_y[i];

        ts_xxxxy_yyyy[i] = 3.0 * ts_xxy_yyyy[i] * gfe_0 + ts_xxxy_yyyy[i] * ga_x[i];

        ts_xxxxy_yyyz[i] = 3.0 * ts_xxy_yyyz[i] * gfe_0 + ts_xxxy_yyyz[i] * ga_x[i];

        ts_xxxxy_yyzz[i] = 3.0 * ts_xxy_yyzz[i] * gfe_0 + ts_xxxy_yyzz[i] * ga_x[i];

        ts_xxxxy_yzzz[i] = 3.0 * ts_xxy_yzzz[i] * gfe_0 + ts_xxxy_yzzz[i] * ga_x[i];

        ts_xxxxy_zzzz[i] = ts_xxxx_zzzz[i] * ga_y[i];
    }

    // Set up 30-45 components of targeted buffer : HG

    auto ts_xxxxz_xxxx = pbuffer.data(idx_hg + 30);

    auto ts_xxxxz_xxxy = pbuffer.data(idx_hg + 31);

    auto ts_xxxxz_xxxz = pbuffer.data(idx_hg + 32);

    auto ts_xxxxz_xxyy = pbuffer.data(idx_hg + 33);

    auto ts_xxxxz_xxyz = pbuffer.data(idx_hg + 34);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_hg + 35);

    auto ts_xxxxz_xyyy = pbuffer.data(idx_hg + 36);

    auto ts_xxxxz_xyyz = pbuffer.data(idx_hg + 37);

    auto ts_xxxxz_xyzz = pbuffer.data(idx_hg + 38);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_hg + 39);

    auto ts_xxxxz_yyyy = pbuffer.data(idx_hg + 40);

    auto ts_xxxxz_yyyz = pbuffer.data(idx_hg + 41);

    auto ts_xxxxz_yyzz = pbuffer.data(idx_hg + 42);

    auto ts_xxxxz_yzzz = pbuffer.data(idx_hg + 43);

    auto ts_xxxxz_zzzz = pbuffer.data(idx_hg + 44);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxxx_xxx, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxy, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxz, ts_xxxx_xxzz, ts_xxxx_xyy, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyz, ts_xxxx_xyzz, ts_xxxx_xzz, ts_xxxx_xzzz, ts_xxxx_yyyy, ts_xxxxz_xxxx, ts_xxxxz_xxxy, ts_xxxxz_xxxz, ts_xxxxz_xxyy, ts_xxxxz_xxyz, ts_xxxxz_xxzz, ts_xxxxz_xyyy, ts_xxxxz_xyyz, ts_xxxxz_xyzz, ts_xxxxz_xzzz, ts_xxxxz_yyyy, ts_xxxxz_yyyz, ts_xxxxz_yyzz, ts_xxxxz_yzzz, ts_xxxxz_zzzz, ts_xxxz_yyyz, ts_xxxz_yyzz, ts_xxxz_yzzz, ts_xxxz_zzzz, ts_xxz_yyyz, ts_xxz_yyzz, ts_xxz_yzzz, ts_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxz_xxxx[i] = ts_xxxx_xxxx[i] * ga_z[i];

        ts_xxxxz_xxxy[i] = ts_xxxx_xxxy[i] * ga_z[i];

        ts_xxxxz_xxxz[i] = ts_xxxx_xxx[i] * gfe_0 + ts_xxxx_xxxz[i] * ga_z[i];

        ts_xxxxz_xxyy[i] = ts_xxxx_xxyy[i] * ga_z[i];

        ts_xxxxz_xxyz[i] = ts_xxxx_xxy[i] * gfe_0 + ts_xxxx_xxyz[i] * ga_z[i];

        ts_xxxxz_xxzz[i] = 2.0 * ts_xxxx_xxz[i] * gfe_0 + ts_xxxx_xxzz[i] * ga_z[i];

        ts_xxxxz_xyyy[i] = ts_xxxx_xyyy[i] * ga_z[i];

        ts_xxxxz_xyyz[i] = ts_xxxx_xyy[i] * gfe_0 + ts_xxxx_xyyz[i] * ga_z[i];

        ts_xxxxz_xyzz[i] = 2.0 * ts_xxxx_xyz[i] * gfe_0 + ts_xxxx_xyzz[i] * ga_z[i];

        ts_xxxxz_xzzz[i] = 3.0 * ts_xxxx_xzz[i] * gfe_0 + ts_xxxx_xzzz[i] * ga_z[i];

        ts_xxxxz_yyyy[i] = ts_xxxx_yyyy[i] * ga_z[i];

        ts_xxxxz_yyyz[i] = 3.0 * ts_xxz_yyyz[i] * gfe_0 + ts_xxxz_yyyz[i] * ga_x[i];

        ts_xxxxz_yyzz[i] = 3.0 * ts_xxz_yyzz[i] * gfe_0 + ts_xxxz_yyzz[i] * ga_x[i];

        ts_xxxxz_yzzz[i] = 3.0 * ts_xxz_yzzz[i] * gfe_0 + ts_xxxz_yzzz[i] * ga_x[i];

        ts_xxxxz_zzzz[i] = 3.0 * ts_xxz_zzzz[i] * gfe_0 + ts_xxxz_zzzz[i] * ga_x[i];
    }

    // Set up 45-60 components of targeted buffer : HG

    auto ts_xxxyy_xxxx = pbuffer.data(idx_hg + 45);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_hg + 46);

    auto ts_xxxyy_xxxz = pbuffer.data(idx_hg + 47);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_hg + 48);

    auto ts_xxxyy_xxyz = pbuffer.data(idx_hg + 49);

    auto ts_xxxyy_xxzz = pbuffer.data(idx_hg + 50);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_hg + 51);

    auto ts_xxxyy_xyyz = pbuffer.data(idx_hg + 52);

    auto ts_xxxyy_xyzz = pbuffer.data(idx_hg + 53);

    auto ts_xxxyy_xzzz = pbuffer.data(idx_hg + 54);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_hg + 58);

    auto ts_xxxyy_zzzz = pbuffer.data(idx_hg + 59);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxx_xxxx, ts_xxx_xxxz, ts_xxx_xxzz, ts_xxx_xzzz, ts_xxxy_xxxx, ts_xxxy_xxxz, ts_xxxy_xxzz, ts_xxxy_xzzz, ts_xxxyy_xxxx, ts_xxxyy_xxxy, ts_xxxyy_xxxz, ts_xxxyy_xxyy, ts_xxxyy_xxyz, ts_xxxyy_xxzz, ts_xxxyy_xyyy, ts_xxxyy_xyyz, ts_xxxyy_xyzz, ts_xxxyy_xzzz, ts_xxxyy_yyyy, ts_xxxyy_yyyz, ts_xxxyy_yyzz, ts_xxxyy_yzzz, ts_xxxyy_zzzz, ts_xxyy_xxxy, ts_xxyy_xxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xyy, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyz, ts_xxyy_xyzz, ts_xxyy_yyy, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyz, ts_xxyy_yyzz, ts_xxyy_yzz, ts_xxyy_yzzz, ts_xxyy_zzzz, ts_xyy_xxxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyzz, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyzz, ts_xyy_yzzz, ts_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxyy_xxxx[i] = ts_xxx_xxxx[i] * gfe_0 + ts_xxxy_xxxx[i] * ga_y[i];

        ts_xxxyy_xxxy[i] = 2.0 * ts_xyy_xxxy[i] * gfe_0 + 3.0 * ts_xxyy_xxy[i] * gfe_0 + ts_xxyy_xxxy[i] * ga_x[i];

        ts_xxxyy_xxxz[i] = ts_xxx_xxxz[i] * gfe_0 + ts_xxxy_xxxz[i] * ga_y[i];

        ts_xxxyy_xxyy[i] = 2.0 * ts_xyy_xxyy[i] * gfe_0 + 2.0 * ts_xxyy_xyy[i] * gfe_0 + ts_xxyy_xxyy[i] * ga_x[i];

        ts_xxxyy_xxyz[i] = 2.0 * ts_xyy_xxyz[i] * gfe_0 + 2.0 * ts_xxyy_xyz[i] * gfe_0 + ts_xxyy_xxyz[i] * ga_x[i];

        ts_xxxyy_xxzz[i] = ts_xxx_xxzz[i] * gfe_0 + ts_xxxy_xxzz[i] * ga_y[i];

        ts_xxxyy_xyyy[i] = 2.0 * ts_xyy_xyyy[i] * gfe_0 + ts_xxyy_yyy[i] * gfe_0 + ts_xxyy_xyyy[i] * ga_x[i];

        ts_xxxyy_xyyz[i] = 2.0 * ts_xyy_xyyz[i] * gfe_0 + ts_xxyy_yyz[i] * gfe_0 + ts_xxyy_xyyz[i] * ga_x[i];

        ts_xxxyy_xyzz[i] = 2.0 * ts_xyy_xyzz[i] * gfe_0 + ts_xxyy_yzz[i] * gfe_0 + ts_xxyy_xyzz[i] * ga_x[i];

        ts_xxxyy_xzzz[i] = ts_xxx_xzzz[i] * gfe_0 + ts_xxxy_xzzz[i] * ga_y[i];

        ts_xxxyy_yyyy[i] = 2.0 * ts_xyy_yyyy[i] * gfe_0 + ts_xxyy_yyyy[i] * ga_x[i];

        ts_xxxyy_yyyz[i] = 2.0 * ts_xyy_yyyz[i] * gfe_0 + ts_xxyy_yyyz[i] * ga_x[i];

        ts_xxxyy_yyzz[i] = 2.0 * ts_xyy_yyzz[i] * gfe_0 + ts_xxyy_yyzz[i] * ga_x[i];

        ts_xxxyy_yzzz[i] = 2.0 * ts_xyy_yzzz[i] * gfe_0 + ts_xxyy_yzzz[i] * ga_x[i];

        ts_xxxyy_zzzz[i] = 2.0 * ts_xyy_zzzz[i] * gfe_0 + ts_xxyy_zzzz[i] * ga_x[i];
    }

    // Set up 60-75 components of targeted buffer : HG

    auto ts_xxxyz_xxxx = pbuffer.data(idx_hg + 60);

    auto ts_xxxyz_xxxy = pbuffer.data(idx_hg + 61);

    auto ts_xxxyz_xxxz = pbuffer.data(idx_hg + 62);

    auto ts_xxxyz_xxyy = pbuffer.data(idx_hg + 63);

    auto ts_xxxyz_xxyz = pbuffer.data(idx_hg + 64);

    auto ts_xxxyz_xxzz = pbuffer.data(idx_hg + 65);

    auto ts_xxxyz_xyyy = pbuffer.data(idx_hg + 66);

    auto ts_xxxyz_xyyz = pbuffer.data(idx_hg + 67);

    auto ts_xxxyz_xyzz = pbuffer.data(idx_hg + 68);

    auto ts_xxxyz_xzzz = pbuffer.data(idx_hg + 69);

    auto ts_xxxyz_yyyy = pbuffer.data(idx_hg + 70);

    auto ts_xxxyz_yyyz = pbuffer.data(idx_hg + 71);

    auto ts_xxxyz_yyzz = pbuffer.data(idx_hg + 72);

    auto ts_xxxyz_yzzz = pbuffer.data(idx_hg + 73);

    auto ts_xxxyz_zzzz = pbuffer.data(idx_hg + 74);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_xxxy_xxxy, ts_xxxy_xxyy, ts_xxxy_xyyy, ts_xxxy_yyyy, ts_xxxyz_xxxx, ts_xxxyz_xxxy, ts_xxxyz_xxxz, ts_xxxyz_xxyy, ts_xxxyz_xxyz, ts_xxxyz_xxzz, ts_xxxyz_xyyy, ts_xxxyz_xyyz, ts_xxxyz_xyzz, ts_xxxyz_xzzz, ts_xxxyz_yyyy, ts_xxxyz_yyyz, ts_xxxyz_yyzz, ts_xxxyz_yzzz, ts_xxxyz_zzzz, ts_xxxz_xxxx, ts_xxxz_xxxz, ts_xxxz_xxyz, ts_xxxz_xxz, ts_xxxz_xxzz, ts_xxxz_xyyz, ts_xxxz_xyz, ts_xxxz_xyzz, ts_xxxz_xzz, ts_xxxz_xzzz, ts_xxxz_zzzz, ts_xxyz_yyyz, ts_xxyz_yyzz, ts_xxyz_yzzz, ts_xyz_yyyz, ts_xyz_yyzz, ts_xyz_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxyz_xxxx[i] = ts_xxxz_xxxx[i] * ga_y[i];

        ts_xxxyz_xxxy[i] = ts_xxxy_xxxy[i] * ga_z[i];

        ts_xxxyz_xxxz[i] = ts_xxxz_xxxz[i] * ga_y[i];

        ts_xxxyz_xxyy[i] = ts_xxxy_xxyy[i] * ga_z[i];

        ts_xxxyz_xxyz[i] = ts_xxxz_xxz[i] * gfe_0 + ts_xxxz_xxyz[i] * ga_y[i];

        ts_xxxyz_xxzz[i] = ts_xxxz_xxzz[i] * ga_y[i];

        ts_xxxyz_xyyy[i] = ts_xxxy_xyyy[i] * ga_z[i];

        ts_xxxyz_xyyz[i] = 2.0 * ts_xxxz_xyz[i] * gfe_0 + ts_xxxz_xyyz[i] * ga_y[i];

        ts_xxxyz_xyzz[i] = ts_xxxz_xzz[i] * gfe_0 + ts_xxxz_xyzz[i] * ga_y[i];

        ts_xxxyz_xzzz[i] = ts_xxxz_xzzz[i] * ga_y[i];

        ts_xxxyz_yyyy[i] = ts_xxxy_yyyy[i] * ga_z[i];

        ts_xxxyz_yyyz[i] = 2.0 * ts_xyz_yyyz[i] * gfe_0 + ts_xxyz_yyyz[i] * ga_x[i];

        ts_xxxyz_yyzz[i] = 2.0 * ts_xyz_yyzz[i] * gfe_0 + ts_xxyz_yyzz[i] * ga_x[i];

        ts_xxxyz_yzzz[i] = 2.0 * ts_xyz_yzzz[i] * gfe_0 + ts_xxyz_yzzz[i] * ga_x[i];

        ts_xxxyz_zzzz[i] = ts_xxxz_zzzz[i] * ga_y[i];
    }

    // Set up 75-90 components of targeted buffer : HG

    auto ts_xxxzz_xxxx = pbuffer.data(idx_hg + 75);

    auto ts_xxxzz_xxxy = pbuffer.data(idx_hg + 76);

    auto ts_xxxzz_xxxz = pbuffer.data(idx_hg + 77);

    auto ts_xxxzz_xxyy = pbuffer.data(idx_hg + 78);

    auto ts_xxxzz_xxyz = pbuffer.data(idx_hg + 79);

    auto ts_xxxzz_xxzz = pbuffer.data(idx_hg + 80);

    auto ts_xxxzz_xyyy = pbuffer.data(idx_hg + 81);

    auto ts_xxxzz_xyyz = pbuffer.data(idx_hg + 82);

    auto ts_xxxzz_xyzz = pbuffer.data(idx_hg + 83);

    auto ts_xxxzz_xzzz = pbuffer.data(idx_hg + 84);

    auto ts_xxxzz_yyyy = pbuffer.data(idx_hg + 85);

    auto ts_xxxzz_yyyz = pbuffer.data(idx_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_hg + 89);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxyy, ts_xxx_xyyy, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxyy, ts_xxxz_xyyy, ts_xxxzz_xxxx, ts_xxxzz_xxxy, ts_xxxzz_xxxz, ts_xxxzz_xxyy, ts_xxxzz_xxyz, ts_xxxzz_xxzz, ts_xxxzz_xyyy, ts_xxxzz_xyyz, ts_xxxzz_xyzz, ts_xxxzz_xzzz, ts_xxxzz_yyyy, ts_xxxzz_yyyz, ts_xxxzz_yyzz, ts_xxxzz_yzzz, ts_xxxzz_zzzz, ts_xxzz_xxxz, ts_xxzz_xxyz, ts_xxzz_xxz, ts_xxzz_xxzz, ts_xxzz_xyyz, ts_xxzz_xyz, ts_xxzz_xyzz, ts_xxzz_xzz, ts_xxzz_xzzz, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyz, ts_xxzz_yyzz, ts_xxzz_yzz, ts_xxzz_yzzz, ts_xxzz_zzz, ts_xxzz_zzzz, ts_xzz_xxxz, ts_xzz_xxyz, ts_xzz_xxzz, ts_xzz_xyyz, ts_xzz_xyzz, ts_xzz_xzzz, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyzz, ts_xzz_yzzz, ts_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxzz_xxxx[i] = ts_xxx_xxxx[i] * gfe_0 + ts_xxxz_xxxx[i] * ga_z[i];

        ts_xxxzz_xxxy[i] = ts_xxx_xxxy[i] * gfe_0 + ts_xxxz_xxxy[i] * ga_z[i];

        ts_xxxzz_xxxz[i] = 2.0 * ts_xzz_xxxz[i] * gfe_0 + 3.0 * ts_xxzz_xxz[i] * gfe_0 + ts_xxzz_xxxz[i] * ga_x[i];

        ts_xxxzz_xxyy[i] = ts_xxx_xxyy[i] * gfe_0 + ts_xxxz_xxyy[i] * ga_z[i];

        ts_xxxzz_xxyz[i] = 2.0 * ts_xzz_xxyz[i] * gfe_0 + 2.0 * ts_xxzz_xyz[i] * gfe_0 + ts_xxzz_xxyz[i] * ga_x[i];

        ts_xxxzz_xxzz[i] = 2.0 * ts_xzz_xxzz[i] * gfe_0 + 2.0 * ts_xxzz_xzz[i] * gfe_0 + ts_xxzz_xxzz[i] * ga_x[i];

        ts_xxxzz_xyyy[i] = ts_xxx_xyyy[i] * gfe_0 + ts_xxxz_xyyy[i] * ga_z[i];

        ts_xxxzz_xyyz[i] = 2.0 * ts_xzz_xyyz[i] * gfe_0 + ts_xxzz_yyz[i] * gfe_0 + ts_xxzz_xyyz[i] * ga_x[i];

        ts_xxxzz_xyzz[i] = 2.0 * ts_xzz_xyzz[i] * gfe_0 + ts_xxzz_yzz[i] * gfe_0 + ts_xxzz_xyzz[i] * ga_x[i];

        ts_xxxzz_xzzz[i] = 2.0 * ts_xzz_xzzz[i] * gfe_0 + ts_xxzz_zzz[i] * gfe_0 + ts_xxzz_xzzz[i] * ga_x[i];

        ts_xxxzz_yyyy[i] = 2.0 * ts_xzz_yyyy[i] * gfe_0 + ts_xxzz_yyyy[i] * ga_x[i];

        ts_xxxzz_yyyz[i] = 2.0 * ts_xzz_yyyz[i] * gfe_0 + ts_xxzz_yyyz[i] * ga_x[i];

        ts_xxxzz_yyzz[i] = 2.0 * ts_xzz_yyzz[i] * gfe_0 + ts_xxzz_yyzz[i] * ga_x[i];

        ts_xxxzz_yzzz[i] = 2.0 * ts_xzz_yzzz[i] * gfe_0 + ts_xxzz_yzzz[i] * ga_x[i];

        ts_xxxzz_zzzz[i] = 2.0 * ts_xzz_zzzz[i] * gfe_0 + ts_xxzz_zzzz[i] * ga_x[i];
    }

    // Set up 90-105 components of targeted buffer : HG

    auto ts_xxyyy_xxxx = pbuffer.data(idx_hg + 90);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_hg + 91);

    auto ts_xxyyy_xxxz = pbuffer.data(idx_hg + 92);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_hg + 93);

    auto ts_xxyyy_xxyz = pbuffer.data(idx_hg + 94);

    auto ts_xxyyy_xxzz = pbuffer.data(idx_hg + 95);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_hg + 96);

    auto ts_xxyyy_xyyz = pbuffer.data(idx_hg + 97);

    auto ts_xxyyy_xyzz = pbuffer.data(idx_hg + 98);

    auto ts_xxyyy_xzzz = pbuffer.data(idx_hg + 99);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_hg + 103);

    auto ts_xxyyy_zzzz = pbuffer.data(idx_hg + 104);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxy_xxxx, ts_xxy_xxxz, ts_xxy_xxzz, ts_xxy_xzzz, ts_xxyy_xxxx, ts_xxyy_xxxz, ts_xxyy_xxzz, ts_xxyy_xzzz, ts_xxyyy_xxxx, ts_xxyyy_xxxy, ts_xxyyy_xxxz, ts_xxyyy_xxyy, ts_xxyyy_xxyz, ts_xxyyy_xxzz, ts_xxyyy_xyyy, ts_xxyyy_xyyz, ts_xxyyy_xyzz, ts_xxyyy_xzzz, ts_xxyyy_yyyy, ts_xxyyy_yyyz, ts_xxyyy_yyzz, ts_xxyyy_yzzz, ts_xxyyy_zzzz, ts_xyyy_xxxy, ts_xyyy_xxy, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xyy, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyz, ts_xyyy_xyzz, ts_xyyy_yyy, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyz, ts_xyyy_yyzz, ts_xyyy_yzz, ts_xyyy_yzzz, ts_xyyy_zzzz, ts_yyy_xxxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyzz, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyzz, ts_yyy_yzzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyyy_xxxx[i] = 2.0 * ts_xxy_xxxx[i] * gfe_0 + ts_xxyy_xxxx[i] * ga_y[i];

        ts_xxyyy_xxxy[i] = ts_yyy_xxxy[i] * gfe_0 + 3.0 * ts_xyyy_xxy[i] * gfe_0 + ts_xyyy_xxxy[i] * ga_x[i];

        ts_xxyyy_xxxz[i] = 2.0 * ts_xxy_xxxz[i] * gfe_0 + ts_xxyy_xxxz[i] * ga_y[i];

        ts_xxyyy_xxyy[i] = ts_yyy_xxyy[i] * gfe_0 + 2.0 * ts_xyyy_xyy[i] * gfe_0 + ts_xyyy_xxyy[i] * ga_x[i];

        ts_xxyyy_xxyz[i] = ts_yyy_xxyz[i] * gfe_0 + 2.0 * ts_xyyy_xyz[i] * gfe_0 + ts_xyyy_xxyz[i] * ga_x[i];

        ts_xxyyy_xxzz[i] = 2.0 * ts_xxy_xxzz[i] * gfe_0 + ts_xxyy_xxzz[i] * ga_y[i];

        ts_xxyyy_xyyy[i] = ts_yyy_xyyy[i] * gfe_0 + ts_xyyy_yyy[i] * gfe_0 + ts_xyyy_xyyy[i] * ga_x[i];

        ts_xxyyy_xyyz[i] = ts_yyy_xyyz[i] * gfe_0 + ts_xyyy_yyz[i] * gfe_0 + ts_xyyy_xyyz[i] * ga_x[i];

        ts_xxyyy_xyzz[i] = ts_yyy_xyzz[i] * gfe_0 + ts_xyyy_yzz[i] * gfe_0 + ts_xyyy_xyzz[i] * ga_x[i];

        ts_xxyyy_xzzz[i] = 2.0 * ts_xxy_xzzz[i] * gfe_0 + ts_xxyy_xzzz[i] * ga_y[i];

        ts_xxyyy_yyyy[i] = ts_yyy_yyyy[i] * gfe_0 + ts_xyyy_yyyy[i] * ga_x[i];

        ts_xxyyy_yyyz[i] = ts_yyy_yyyz[i] * gfe_0 + ts_xyyy_yyyz[i] * ga_x[i];

        ts_xxyyy_yyzz[i] = ts_yyy_yyzz[i] * gfe_0 + ts_xyyy_yyzz[i] * ga_x[i];

        ts_xxyyy_yzzz[i] = ts_yyy_yzzz[i] * gfe_0 + ts_xyyy_yzzz[i] * ga_x[i];

        ts_xxyyy_zzzz[i] = ts_yyy_zzzz[i] * gfe_0 + ts_xyyy_zzzz[i] * ga_x[i];
    }

    // Set up 105-120 components of targeted buffer : HG

    auto ts_xxyyz_xxxx = pbuffer.data(idx_hg + 105);

    auto ts_xxyyz_xxxy = pbuffer.data(idx_hg + 106);

    auto ts_xxyyz_xxxz = pbuffer.data(idx_hg + 107);

    auto ts_xxyyz_xxyy = pbuffer.data(idx_hg + 108);

    auto ts_xxyyz_xxyz = pbuffer.data(idx_hg + 109);

    auto ts_xxyyz_xxzz = pbuffer.data(idx_hg + 110);

    auto ts_xxyyz_xyyy = pbuffer.data(idx_hg + 111);

    auto ts_xxyyz_xyyz = pbuffer.data(idx_hg + 112);

    auto ts_xxyyz_xyzz = pbuffer.data(idx_hg + 113);

    auto ts_xxyyz_xzzz = pbuffer.data(idx_hg + 114);

    auto ts_xxyyz_yyyy = pbuffer.data(idx_hg + 115);

    auto ts_xxyyz_yyyz = pbuffer.data(idx_hg + 116);

    auto ts_xxyyz_yyzz = pbuffer.data(idx_hg + 117);

    auto ts_xxyyz_yzzz = pbuffer.data(idx_hg + 118);

    auto ts_xxyyz_zzzz = pbuffer.data(idx_hg + 119);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xyy, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyz, ts_xxyy_xyzz, ts_xxyy_yyyy, ts_xxyyz_xxxx, ts_xxyyz_xxxy, ts_xxyyz_xxxz, ts_xxyyz_xxyy, ts_xxyyz_xxyz, ts_xxyyz_xxzz, ts_xxyyz_xyyy, ts_xxyyz_xyyz, ts_xxyyz_xyzz, ts_xxyyz_xzzz, ts_xxyyz_yyyy, ts_xxyyz_yyyz, ts_xxyyz_yyzz, ts_xxyyz_yzzz, ts_xxyyz_zzzz, ts_xxyz_xxxz, ts_xxyz_xxzz, ts_xxyz_xzzz, ts_xxz_xxxz, ts_xxz_xxzz, ts_xxz_xzzz, ts_xyyz_yyyz, ts_xyyz_yyzz, ts_xyyz_yzzz, ts_xyyz_zzzz, ts_yyz_yyyz, ts_yyz_yyzz, ts_yyz_yzzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyyz_xxxx[i] = ts_xxyy_xxxx[i] * ga_z[i];

        ts_xxyyz_xxxy[i] = ts_xxyy_xxxy[i] * ga_z[i];

        ts_xxyyz_xxxz[i] = ts_xxz_xxxz[i] * gfe_0 + ts_xxyz_xxxz[i] * ga_y[i];

        ts_xxyyz_xxyy[i] = ts_xxyy_xxyy[i] * ga_z[i];

        ts_xxyyz_xxyz[i] = ts_xxyy_xxy[i] * gfe_0 + ts_xxyy_xxyz[i] * ga_z[i];

        ts_xxyyz_xxzz[i] = ts_xxz_xxzz[i] * gfe_0 + ts_xxyz_xxzz[i] * ga_y[i];

        ts_xxyyz_xyyy[i] = ts_xxyy_xyyy[i] * ga_z[i];

        ts_xxyyz_xyyz[i] = ts_xxyy_xyy[i] * gfe_0 + ts_xxyy_xyyz[i] * ga_z[i];

        ts_xxyyz_xyzz[i] = 2.0 * ts_xxyy_xyz[i] * gfe_0 + ts_xxyy_xyzz[i] * ga_z[i];

        ts_xxyyz_xzzz[i] = ts_xxz_xzzz[i] * gfe_0 + ts_xxyz_xzzz[i] * ga_y[i];

        ts_xxyyz_yyyy[i] = ts_xxyy_yyyy[i] * ga_z[i];

        ts_xxyyz_yyyz[i] = ts_yyz_yyyz[i] * gfe_0 + ts_xyyz_yyyz[i] * ga_x[i];

        ts_xxyyz_yyzz[i] = ts_yyz_yyzz[i] * gfe_0 + ts_xyyz_yyzz[i] * ga_x[i];

        ts_xxyyz_yzzz[i] = ts_yyz_yzzz[i] * gfe_0 + ts_xyyz_yzzz[i] * ga_x[i];

        ts_xxyyz_zzzz[i] = ts_yyz_zzzz[i] * gfe_0 + ts_xyyz_zzzz[i] * ga_x[i];
    }

    // Set up 120-135 components of targeted buffer : HG

    auto ts_xxyzz_xxxx = pbuffer.data(idx_hg + 120);

    auto ts_xxyzz_xxxy = pbuffer.data(idx_hg + 121);

    auto ts_xxyzz_xxxz = pbuffer.data(idx_hg + 122);

    auto ts_xxyzz_xxyy = pbuffer.data(idx_hg + 123);

    auto ts_xxyzz_xxyz = pbuffer.data(idx_hg + 124);

    auto ts_xxyzz_xxzz = pbuffer.data(idx_hg + 125);

    auto ts_xxyzz_xyyy = pbuffer.data(idx_hg + 126);

    auto ts_xxyzz_xyyz = pbuffer.data(idx_hg + 127);

    auto ts_xxyzz_xyzz = pbuffer.data(idx_hg + 128);

    auto ts_xxyzz_xzzz = pbuffer.data(idx_hg + 129);

    auto ts_xxyzz_yyyy = pbuffer.data(idx_hg + 130);

    auto ts_xxyzz_yyyz = pbuffer.data(idx_hg + 131);

    auto ts_xxyzz_yyzz = pbuffer.data(idx_hg + 132);

    auto ts_xxyzz_yzzz = pbuffer.data(idx_hg + 133);

    auto ts_xxyzz_zzzz = pbuffer.data(idx_hg + 134);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxyzz_xxxx, ts_xxyzz_xxxy, ts_xxyzz_xxxz, ts_xxyzz_xxyy, ts_xxyzz_xxyz, ts_xxyzz_xxzz, ts_xxyzz_xyyy, ts_xxyzz_xyyz, ts_xxyzz_xyzz, ts_xxyzz_xzzz, ts_xxyzz_yyyy, ts_xxyzz_yyyz, ts_xxyzz_yyzz, ts_xxyzz_yzzz, ts_xxyzz_zzzz, ts_xxzz_xxx, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxy, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxz, ts_xxzz_xxzz, ts_xxzz_xyy, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyz, ts_xxzz_xyzz, ts_xxzz_xzz, ts_xxzz_xzzz, ts_xxzz_zzzz, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyzz, ts_xyzz_yzzz, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyzz, ts_yzz_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyzz_xxxx[i] = ts_xxzz_xxxx[i] * ga_y[i];

        ts_xxyzz_xxxy[i] = ts_xxzz_xxx[i] * gfe_0 + ts_xxzz_xxxy[i] * ga_y[i];

        ts_xxyzz_xxxz[i] = ts_xxzz_xxxz[i] * ga_y[i];

        ts_xxyzz_xxyy[i] = 2.0 * ts_xxzz_xxy[i] * gfe_0 + ts_xxzz_xxyy[i] * ga_y[i];

        ts_xxyzz_xxyz[i] = ts_xxzz_xxz[i] * gfe_0 + ts_xxzz_xxyz[i] * ga_y[i];

        ts_xxyzz_xxzz[i] = ts_xxzz_xxzz[i] * ga_y[i];

        ts_xxyzz_xyyy[i] = 3.0 * ts_xxzz_xyy[i] * gfe_0 + ts_xxzz_xyyy[i] * ga_y[i];

        ts_xxyzz_xyyz[i] = 2.0 * ts_xxzz_xyz[i] * gfe_0 + ts_xxzz_xyyz[i] * ga_y[i];

        ts_xxyzz_xyzz[i] = ts_xxzz_xzz[i] * gfe_0 + ts_xxzz_xyzz[i] * ga_y[i];

        ts_xxyzz_xzzz[i] = ts_xxzz_xzzz[i] * ga_y[i];

        ts_xxyzz_yyyy[i] = ts_yzz_yyyy[i] * gfe_0 + ts_xyzz_yyyy[i] * ga_x[i];

        ts_xxyzz_yyyz[i] = ts_yzz_yyyz[i] * gfe_0 + ts_xyzz_yyyz[i] * ga_x[i];

        ts_xxyzz_yyzz[i] = ts_yzz_yyzz[i] * gfe_0 + ts_xyzz_yyzz[i] * ga_x[i];

        ts_xxyzz_yzzz[i] = ts_yzz_yzzz[i] * gfe_0 + ts_xyzz_yzzz[i] * ga_x[i];

        ts_xxyzz_zzzz[i] = ts_xxzz_zzzz[i] * ga_y[i];
    }

    // Set up 135-150 components of targeted buffer : HG

    auto ts_xxzzz_xxxx = pbuffer.data(idx_hg + 135);

    auto ts_xxzzz_xxxy = pbuffer.data(idx_hg + 136);

    auto ts_xxzzz_xxxz = pbuffer.data(idx_hg + 137);

    auto ts_xxzzz_xxyy = pbuffer.data(idx_hg + 138);

    auto ts_xxzzz_xxyz = pbuffer.data(idx_hg + 139);

    auto ts_xxzzz_xxzz = pbuffer.data(idx_hg + 140);

    auto ts_xxzzz_xyyy = pbuffer.data(idx_hg + 141);

    auto ts_xxzzz_xyyz = pbuffer.data(idx_hg + 142);

    auto ts_xxzzz_xyzz = pbuffer.data(idx_hg + 143);

    auto ts_xxzzz_xzzz = pbuffer.data(idx_hg + 144);

    auto ts_xxzzz_yyyy = pbuffer.data(idx_hg + 145);

    auto ts_xxzzz_yyyz = pbuffer.data(idx_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_hg + 149);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxyy, ts_xxz_xyyy, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxyy, ts_xxzz_xyyy, ts_xxzzz_xxxx, ts_xxzzz_xxxy, ts_xxzzz_xxxz, ts_xxzzz_xxyy, ts_xxzzz_xxyz, ts_xxzzz_xxzz, ts_xxzzz_xyyy, ts_xxzzz_xyyz, ts_xxzzz_xyzz, ts_xxzzz_xzzz, ts_xxzzz_yyyy, ts_xxzzz_yyyz, ts_xxzzz_yyzz, ts_xxzzz_yzzz, ts_xxzzz_zzzz, ts_xzzz_xxxz, ts_xzzz_xxyz, ts_xzzz_xxz, ts_xzzz_xxzz, ts_xzzz_xyyz, ts_xzzz_xyz, ts_xzzz_xyzz, ts_xzzz_xzz, ts_xzzz_xzzz, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyz, ts_xzzz_yyzz, ts_xzzz_yzz, ts_xzzz_yzzz, ts_xzzz_zzz, ts_xzzz_zzzz, ts_zzz_xxxz, ts_zzz_xxyz, ts_zzz_xxzz, ts_zzz_xyyz, ts_zzz_xyzz, ts_zzz_xzzz, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyzz, ts_zzz_yzzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxzzz_xxxx[i] = 2.0 * ts_xxz_xxxx[i] * gfe_0 + ts_xxzz_xxxx[i] * ga_z[i];

        ts_xxzzz_xxxy[i] = 2.0 * ts_xxz_xxxy[i] * gfe_0 + ts_xxzz_xxxy[i] * ga_z[i];

        ts_xxzzz_xxxz[i] = ts_zzz_xxxz[i] * gfe_0 + 3.0 * ts_xzzz_xxz[i] * gfe_0 + ts_xzzz_xxxz[i] * ga_x[i];

        ts_xxzzz_xxyy[i] = 2.0 * ts_xxz_xxyy[i] * gfe_0 + ts_xxzz_xxyy[i] * ga_z[i];

        ts_xxzzz_xxyz[i] = ts_zzz_xxyz[i] * gfe_0 + 2.0 * ts_xzzz_xyz[i] * gfe_0 + ts_xzzz_xxyz[i] * ga_x[i];

        ts_xxzzz_xxzz[i] = ts_zzz_xxzz[i] * gfe_0 + 2.0 * ts_xzzz_xzz[i] * gfe_0 + ts_xzzz_xxzz[i] * ga_x[i];

        ts_xxzzz_xyyy[i] = 2.0 * ts_xxz_xyyy[i] * gfe_0 + ts_xxzz_xyyy[i] * ga_z[i];

        ts_xxzzz_xyyz[i] = ts_zzz_xyyz[i] * gfe_0 + ts_xzzz_yyz[i] * gfe_0 + ts_xzzz_xyyz[i] * ga_x[i];

        ts_xxzzz_xyzz[i] = ts_zzz_xyzz[i] * gfe_0 + ts_xzzz_yzz[i] * gfe_0 + ts_xzzz_xyzz[i] * ga_x[i];

        ts_xxzzz_xzzz[i] = ts_zzz_xzzz[i] * gfe_0 + ts_xzzz_zzz[i] * gfe_0 + ts_xzzz_xzzz[i] * ga_x[i];

        ts_xxzzz_yyyy[i] = ts_zzz_yyyy[i] * gfe_0 + ts_xzzz_yyyy[i] * ga_x[i];

        ts_xxzzz_yyyz[i] = ts_zzz_yyyz[i] * gfe_0 + ts_xzzz_yyyz[i] * ga_x[i];

        ts_xxzzz_yyzz[i] = ts_zzz_yyzz[i] * gfe_0 + ts_xzzz_yyzz[i] * ga_x[i];

        ts_xxzzz_yzzz[i] = ts_zzz_yzzz[i] * gfe_0 + ts_xzzz_yzzz[i] * ga_x[i];

        ts_xxzzz_zzzz[i] = ts_zzz_zzzz[i] * gfe_0 + ts_xzzz_zzzz[i] * ga_x[i];
    }

    // Set up 150-165 components of targeted buffer : HG

    auto ts_xyyyy_xxxx = pbuffer.data(idx_hg + 150);

    auto ts_xyyyy_xxxy = pbuffer.data(idx_hg + 151);

    auto ts_xyyyy_xxxz = pbuffer.data(idx_hg + 152);

    auto ts_xyyyy_xxyy = pbuffer.data(idx_hg + 153);

    auto ts_xyyyy_xxyz = pbuffer.data(idx_hg + 154);

    auto ts_xyyyy_xxzz = pbuffer.data(idx_hg + 155);

    auto ts_xyyyy_xyyy = pbuffer.data(idx_hg + 156);

    auto ts_xyyyy_xyyz = pbuffer.data(idx_hg + 157);

    auto ts_xyyyy_xyzz = pbuffer.data(idx_hg + 158);

    auto ts_xyyyy_xzzz = pbuffer.data(idx_hg + 159);

    auto ts_xyyyy_yyyy = pbuffer.data(idx_hg + 160);

    auto ts_xyyyy_yyyz = pbuffer.data(idx_hg + 161);

    auto ts_xyyyy_yyzz = pbuffer.data(idx_hg + 162);

    auto ts_xyyyy_yzzz = pbuffer.data(idx_hg + 163);

    auto ts_xyyyy_zzzz = pbuffer.data(idx_hg + 164);

    #pragma omp simd aligned(ga_x, ts_xyyyy_xxxx, ts_xyyyy_xxxy, ts_xyyyy_xxxz, ts_xyyyy_xxyy, ts_xyyyy_xxyz, ts_xyyyy_xxzz, ts_xyyyy_xyyy, ts_xyyyy_xyyz, ts_xyyyy_xyzz, ts_xyyyy_xzzz, ts_xyyyy_yyyy, ts_xyyyy_yyyz, ts_xyyyy_yyzz, ts_xyyyy_yzzz, ts_xyyyy_zzzz, ts_yyyy_xxx, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxz, ts_yyyy_xxzz, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_xzz, ts_yyyy_xzzz, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyy_zzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyyy_xxxx[i] = 4.0 * ts_yyyy_xxx[i] * gfe_0 + ts_yyyy_xxxx[i] * ga_x[i];

        ts_xyyyy_xxxy[i] = 3.0 * ts_yyyy_xxy[i] * gfe_0 + ts_yyyy_xxxy[i] * ga_x[i];

        ts_xyyyy_xxxz[i] = 3.0 * ts_yyyy_xxz[i] * gfe_0 + ts_yyyy_xxxz[i] * ga_x[i];

        ts_xyyyy_xxyy[i] = 2.0 * ts_yyyy_xyy[i] * gfe_0 + ts_yyyy_xxyy[i] * ga_x[i];

        ts_xyyyy_xxyz[i] = 2.0 * ts_yyyy_xyz[i] * gfe_0 + ts_yyyy_xxyz[i] * ga_x[i];

        ts_xyyyy_xxzz[i] = 2.0 * ts_yyyy_xzz[i] * gfe_0 + ts_yyyy_xxzz[i] * ga_x[i];

        ts_xyyyy_xyyy[i] = ts_yyyy_yyy[i] * gfe_0 + ts_yyyy_xyyy[i] * ga_x[i];

        ts_xyyyy_xyyz[i] = ts_yyyy_yyz[i] * gfe_0 + ts_yyyy_xyyz[i] * ga_x[i];

        ts_xyyyy_xyzz[i] = ts_yyyy_yzz[i] * gfe_0 + ts_yyyy_xyzz[i] * ga_x[i];

        ts_xyyyy_xzzz[i] = ts_yyyy_zzz[i] * gfe_0 + ts_yyyy_xzzz[i] * ga_x[i];

        ts_xyyyy_yyyy[i] = ts_yyyy_yyyy[i] * ga_x[i];

        ts_xyyyy_yyyz[i] = ts_yyyy_yyyz[i] * ga_x[i];

        ts_xyyyy_yyzz[i] = ts_yyyy_yyzz[i] * ga_x[i];

        ts_xyyyy_yzzz[i] = ts_yyyy_yzzz[i] * ga_x[i];

        ts_xyyyy_zzzz[i] = ts_yyyy_zzzz[i] * ga_x[i];
    }

    // Set up 165-180 components of targeted buffer : HG

    auto ts_xyyyz_xxxx = pbuffer.data(idx_hg + 165);

    auto ts_xyyyz_xxxy = pbuffer.data(idx_hg + 166);

    auto ts_xyyyz_xxxz = pbuffer.data(idx_hg + 167);

    auto ts_xyyyz_xxyy = pbuffer.data(idx_hg + 168);

    auto ts_xyyyz_xxyz = pbuffer.data(idx_hg + 169);

    auto ts_xyyyz_xxzz = pbuffer.data(idx_hg + 170);

    auto ts_xyyyz_xyyy = pbuffer.data(idx_hg + 171);

    auto ts_xyyyz_xyyz = pbuffer.data(idx_hg + 172);

    auto ts_xyyyz_xyzz = pbuffer.data(idx_hg + 173);

    auto ts_xyyyz_xzzz = pbuffer.data(idx_hg + 174);

    auto ts_xyyyz_yyyy = pbuffer.data(idx_hg + 175);

    auto ts_xyyyz_yyyz = pbuffer.data(idx_hg + 176);

    auto ts_xyyyz_yyzz = pbuffer.data(idx_hg + 177);

    auto ts_xyyyz_yzzz = pbuffer.data(idx_hg + 178);

    auto ts_xyyyz_zzzz = pbuffer.data(idx_hg + 179);

    #pragma omp simd aligned(ga_x, ga_z, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxyy, ts_xyyy_xyyy, ts_xyyyz_xxxx, ts_xyyyz_xxxy, ts_xyyyz_xxxz, ts_xyyyz_xxyy, ts_xyyyz_xxyz, ts_xyyyz_xxzz, ts_xyyyz_xyyy, ts_xyyyz_xyyz, ts_xyyyz_xyzz, ts_xyyyz_xzzz, ts_xyyyz_yyyy, ts_xyyyz_yyyz, ts_xyyyz_yyzz, ts_xyyyz_yzzz, ts_xyyyz_zzzz, ts_yyyz_xxxz, ts_yyyz_xxyz, ts_yyyz_xxz, ts_yyyz_xxzz, ts_yyyz_xyyz, ts_yyyz_xyz, ts_yyyz_xyzz, ts_yyyz_xzz, ts_yyyz_xzzz, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyz, ts_yyyz_yyzz, ts_yyyz_yzz, ts_yyyz_yzzz, ts_yyyz_zzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyyz_xxxx[i] = ts_xyyy_xxxx[i] * ga_z[i];

        ts_xyyyz_xxxy[i] = ts_xyyy_xxxy[i] * ga_z[i];

        ts_xyyyz_xxxz[i] = 3.0 * ts_yyyz_xxz[i] * gfe_0 + ts_yyyz_xxxz[i] * ga_x[i];

        ts_xyyyz_xxyy[i] = ts_xyyy_xxyy[i] * ga_z[i];

        ts_xyyyz_xxyz[i] = 2.0 * ts_yyyz_xyz[i] * gfe_0 + ts_yyyz_xxyz[i] * ga_x[i];

        ts_xyyyz_xxzz[i] = 2.0 * ts_yyyz_xzz[i] * gfe_0 + ts_yyyz_xxzz[i] * ga_x[i];

        ts_xyyyz_xyyy[i] = ts_xyyy_xyyy[i] * ga_z[i];

        ts_xyyyz_xyyz[i] = ts_yyyz_yyz[i] * gfe_0 + ts_yyyz_xyyz[i] * ga_x[i];

        ts_xyyyz_xyzz[i] = ts_yyyz_yzz[i] * gfe_0 + ts_yyyz_xyzz[i] * ga_x[i];

        ts_xyyyz_xzzz[i] = ts_yyyz_zzz[i] * gfe_0 + ts_yyyz_xzzz[i] * ga_x[i];

        ts_xyyyz_yyyy[i] = ts_yyyz_yyyy[i] * ga_x[i];

        ts_xyyyz_yyyz[i] = ts_yyyz_yyyz[i] * ga_x[i];

        ts_xyyyz_yyzz[i] = ts_yyyz_yyzz[i] * ga_x[i];

        ts_xyyyz_yzzz[i] = ts_yyyz_yzzz[i] * ga_x[i];

        ts_xyyyz_zzzz[i] = ts_yyyz_zzzz[i] * ga_x[i];
    }

    // Set up 180-195 components of targeted buffer : HG

    auto ts_xyyzz_xxxx = pbuffer.data(idx_hg + 180);

    auto ts_xyyzz_xxxy = pbuffer.data(idx_hg + 181);

    auto ts_xyyzz_xxxz = pbuffer.data(idx_hg + 182);

    auto ts_xyyzz_xxyy = pbuffer.data(idx_hg + 183);

    auto ts_xyyzz_xxyz = pbuffer.data(idx_hg + 184);

    auto ts_xyyzz_xxzz = pbuffer.data(idx_hg + 185);

    auto ts_xyyzz_xyyy = pbuffer.data(idx_hg + 186);

    auto ts_xyyzz_xyyz = pbuffer.data(idx_hg + 187);

    auto ts_xyyzz_xyzz = pbuffer.data(idx_hg + 188);

    auto ts_xyyzz_xzzz = pbuffer.data(idx_hg + 189);

    auto ts_xyyzz_yyyy = pbuffer.data(idx_hg + 190);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_hg + 193);

    auto ts_xyyzz_zzzz = pbuffer.data(idx_hg + 194);

    #pragma omp simd aligned(ga_x, ts_xyyzz_xxxx, ts_xyyzz_xxxy, ts_xyyzz_xxxz, ts_xyyzz_xxyy, ts_xyyzz_xxyz, ts_xyyzz_xxzz, ts_xyyzz_xyyy, ts_xyyzz_xyyz, ts_xyyzz_xyzz, ts_xyyzz_xzzz, ts_xyyzz_yyyy, ts_xyyzz_yyyz, ts_xyyzz_yyzz, ts_xyyzz_yzzz, ts_xyyzz_zzzz, ts_yyzz_xxx, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxy, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxz, ts_yyzz_xxzz, ts_yyzz_xyy, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyz, ts_yyzz_xyzz, ts_yyzz_xzz, ts_yyzz_xzzz, ts_yyzz_yyy, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyz, ts_yyzz_yyzz, ts_yyzz_yzz, ts_yyzz_yzzz, ts_yyzz_zzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyzz_xxxx[i] = 4.0 * ts_yyzz_xxx[i] * gfe_0 + ts_yyzz_xxxx[i] * ga_x[i];

        ts_xyyzz_xxxy[i] = 3.0 * ts_yyzz_xxy[i] * gfe_0 + ts_yyzz_xxxy[i] * ga_x[i];

        ts_xyyzz_xxxz[i] = 3.0 * ts_yyzz_xxz[i] * gfe_0 + ts_yyzz_xxxz[i] * ga_x[i];

        ts_xyyzz_xxyy[i] = 2.0 * ts_yyzz_xyy[i] * gfe_0 + ts_yyzz_xxyy[i] * ga_x[i];

        ts_xyyzz_xxyz[i] = 2.0 * ts_yyzz_xyz[i] * gfe_0 + ts_yyzz_xxyz[i] * ga_x[i];

        ts_xyyzz_xxzz[i] = 2.0 * ts_yyzz_xzz[i] * gfe_0 + ts_yyzz_xxzz[i] * ga_x[i];

        ts_xyyzz_xyyy[i] = ts_yyzz_yyy[i] * gfe_0 + ts_yyzz_xyyy[i] * ga_x[i];

        ts_xyyzz_xyyz[i] = ts_yyzz_yyz[i] * gfe_0 + ts_yyzz_xyyz[i] * ga_x[i];

        ts_xyyzz_xyzz[i] = ts_yyzz_yzz[i] * gfe_0 + ts_yyzz_xyzz[i] * ga_x[i];

        ts_xyyzz_xzzz[i] = ts_yyzz_zzz[i] * gfe_0 + ts_yyzz_xzzz[i] * ga_x[i];

        ts_xyyzz_yyyy[i] = ts_yyzz_yyyy[i] * ga_x[i];

        ts_xyyzz_yyyz[i] = ts_yyzz_yyyz[i] * ga_x[i];

        ts_xyyzz_yyzz[i] = ts_yyzz_yyzz[i] * ga_x[i];

        ts_xyyzz_yzzz[i] = ts_yyzz_yzzz[i] * ga_x[i];

        ts_xyyzz_zzzz[i] = ts_yyzz_zzzz[i] * ga_x[i];
    }

    // Set up 195-210 components of targeted buffer : HG

    auto ts_xyzzz_xxxx = pbuffer.data(idx_hg + 195);

    auto ts_xyzzz_xxxy = pbuffer.data(idx_hg + 196);

    auto ts_xyzzz_xxxz = pbuffer.data(idx_hg + 197);

    auto ts_xyzzz_xxyy = pbuffer.data(idx_hg + 198);

    auto ts_xyzzz_xxyz = pbuffer.data(idx_hg + 199);

    auto ts_xyzzz_xxzz = pbuffer.data(idx_hg + 200);

    auto ts_xyzzz_xyyy = pbuffer.data(idx_hg + 201);

    auto ts_xyzzz_xyyz = pbuffer.data(idx_hg + 202);

    auto ts_xyzzz_xyzz = pbuffer.data(idx_hg + 203);

    auto ts_xyzzz_xzzz = pbuffer.data(idx_hg + 204);

    auto ts_xyzzz_yyyy = pbuffer.data(idx_hg + 205);

    auto ts_xyzzz_yyyz = pbuffer.data(idx_hg + 206);

    auto ts_xyzzz_yyzz = pbuffer.data(idx_hg + 207);

    auto ts_xyzzz_yzzz = pbuffer.data(idx_hg + 208);

    auto ts_xyzzz_zzzz = pbuffer.data(idx_hg + 209);

    #pragma omp simd aligned(ga_x, ga_y, ts_xyzzz_xxxx, ts_xyzzz_xxxy, ts_xyzzz_xxxz, ts_xyzzz_xxyy, ts_xyzzz_xxyz, ts_xyzzz_xxzz, ts_xyzzz_xyyy, ts_xyzzz_xyyz, ts_xyzzz_xyzz, ts_xyzzz_xzzz, ts_xyzzz_yyyy, ts_xyzzz_yyyz, ts_xyzzz_yyzz, ts_xyzzz_yzzz, ts_xyzzz_zzzz, ts_xzzz_xxxx, ts_xzzz_xxxz, ts_xzzz_xxzz, ts_xzzz_xzzz, ts_yzzz_xxxy, ts_yzzz_xxy, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xyy, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyz, ts_yzzz_xyzz, ts_yzzz_yyy, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyz, ts_yzzz_yyzz, ts_yzzz_yzz, ts_yzzz_yzzz, ts_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyzzz_xxxx[i] = ts_xzzz_xxxx[i] * ga_y[i];

        ts_xyzzz_xxxy[i] = 3.0 * ts_yzzz_xxy[i] * gfe_0 + ts_yzzz_xxxy[i] * ga_x[i];

        ts_xyzzz_xxxz[i] = ts_xzzz_xxxz[i] * ga_y[i];

        ts_xyzzz_xxyy[i] = 2.0 * ts_yzzz_xyy[i] * gfe_0 + ts_yzzz_xxyy[i] * ga_x[i];

        ts_xyzzz_xxyz[i] = 2.0 * ts_yzzz_xyz[i] * gfe_0 + ts_yzzz_xxyz[i] * ga_x[i];

        ts_xyzzz_xxzz[i] = ts_xzzz_xxzz[i] * ga_y[i];

        ts_xyzzz_xyyy[i] = ts_yzzz_yyy[i] * gfe_0 + ts_yzzz_xyyy[i] * ga_x[i];

        ts_xyzzz_xyyz[i] = ts_yzzz_yyz[i] * gfe_0 + ts_yzzz_xyyz[i] * ga_x[i];

        ts_xyzzz_xyzz[i] = ts_yzzz_yzz[i] * gfe_0 + ts_yzzz_xyzz[i] * ga_x[i];

        ts_xyzzz_xzzz[i] = ts_xzzz_xzzz[i] * ga_y[i];

        ts_xyzzz_yyyy[i] = ts_yzzz_yyyy[i] * ga_x[i];

        ts_xyzzz_yyyz[i] = ts_yzzz_yyyz[i] * ga_x[i];

        ts_xyzzz_yyzz[i] = ts_yzzz_yyzz[i] * ga_x[i];

        ts_xyzzz_yzzz[i] = ts_yzzz_yzzz[i] * ga_x[i];

        ts_xyzzz_zzzz[i] = ts_yzzz_zzzz[i] * ga_x[i];
    }

    // Set up 210-225 components of targeted buffer : HG

    auto ts_xzzzz_xxxx = pbuffer.data(idx_hg + 210);

    auto ts_xzzzz_xxxy = pbuffer.data(idx_hg + 211);

    auto ts_xzzzz_xxxz = pbuffer.data(idx_hg + 212);

    auto ts_xzzzz_xxyy = pbuffer.data(idx_hg + 213);

    auto ts_xzzzz_xxyz = pbuffer.data(idx_hg + 214);

    auto ts_xzzzz_xxzz = pbuffer.data(idx_hg + 215);

    auto ts_xzzzz_xyyy = pbuffer.data(idx_hg + 216);

    auto ts_xzzzz_xyyz = pbuffer.data(idx_hg + 217);

    auto ts_xzzzz_xyzz = pbuffer.data(idx_hg + 218);

    auto ts_xzzzz_xzzz = pbuffer.data(idx_hg + 219);

    auto ts_xzzzz_yyyy = pbuffer.data(idx_hg + 220);

    auto ts_xzzzz_yyyz = pbuffer.data(idx_hg + 221);

    auto ts_xzzzz_yyzz = pbuffer.data(idx_hg + 222);

    auto ts_xzzzz_yzzz = pbuffer.data(idx_hg + 223);

    auto ts_xzzzz_zzzz = pbuffer.data(idx_hg + 224);

    #pragma omp simd aligned(ga_x, ts_xzzzz_xxxx, ts_xzzzz_xxxy, ts_xzzzz_xxxz, ts_xzzzz_xxyy, ts_xzzzz_xxyz, ts_xzzzz_xxzz, ts_xzzzz_xyyy, ts_xzzzz_xyyz, ts_xzzzz_xyzz, ts_xzzzz_xzzz, ts_xzzzz_yyyy, ts_xzzzz_yyyz, ts_xzzzz_yyzz, ts_xzzzz_yzzz, ts_xzzzz_zzzz, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xzzzz_xxxx[i] = 4.0 * ts_zzzz_xxx[i] * gfe_0 + ts_zzzz_xxxx[i] * ga_x[i];

        ts_xzzzz_xxxy[i] = 3.0 * ts_zzzz_xxy[i] * gfe_0 + ts_zzzz_xxxy[i] * ga_x[i];

        ts_xzzzz_xxxz[i] = 3.0 * ts_zzzz_xxz[i] * gfe_0 + ts_zzzz_xxxz[i] * ga_x[i];

        ts_xzzzz_xxyy[i] = 2.0 * ts_zzzz_xyy[i] * gfe_0 + ts_zzzz_xxyy[i] * ga_x[i];

        ts_xzzzz_xxyz[i] = 2.0 * ts_zzzz_xyz[i] * gfe_0 + ts_zzzz_xxyz[i] * ga_x[i];

        ts_xzzzz_xxzz[i] = 2.0 * ts_zzzz_xzz[i] * gfe_0 + ts_zzzz_xxzz[i] * ga_x[i];

        ts_xzzzz_xyyy[i] = ts_zzzz_yyy[i] * gfe_0 + ts_zzzz_xyyy[i] * ga_x[i];

        ts_xzzzz_xyyz[i] = ts_zzzz_yyz[i] * gfe_0 + ts_zzzz_xyyz[i] * ga_x[i];

        ts_xzzzz_xyzz[i] = ts_zzzz_yzz[i] * gfe_0 + ts_zzzz_xyzz[i] * ga_x[i];

        ts_xzzzz_xzzz[i] = ts_zzzz_zzz[i] * gfe_0 + ts_zzzz_xzzz[i] * ga_x[i];

        ts_xzzzz_yyyy[i] = ts_zzzz_yyyy[i] * ga_x[i];

        ts_xzzzz_yyyz[i] = ts_zzzz_yyyz[i] * ga_x[i];

        ts_xzzzz_yyzz[i] = ts_zzzz_yyzz[i] * ga_x[i];

        ts_xzzzz_yzzz[i] = ts_zzzz_yzzz[i] * ga_x[i];

        ts_xzzzz_zzzz[i] = ts_zzzz_zzzz[i] * ga_x[i];
    }

    // Set up 225-240 components of targeted buffer : HG

    auto ts_yyyyy_xxxx = pbuffer.data(idx_hg + 225);

    auto ts_yyyyy_xxxy = pbuffer.data(idx_hg + 226);

    auto ts_yyyyy_xxxz = pbuffer.data(idx_hg + 227);

    auto ts_yyyyy_xxyy = pbuffer.data(idx_hg + 228);

    auto ts_yyyyy_xxyz = pbuffer.data(idx_hg + 229);

    auto ts_yyyyy_xxzz = pbuffer.data(idx_hg + 230);

    auto ts_yyyyy_xyyy = pbuffer.data(idx_hg + 231);

    auto ts_yyyyy_xyyz = pbuffer.data(idx_hg + 232);

    auto ts_yyyyy_xyzz = pbuffer.data(idx_hg + 233);

    auto ts_yyyyy_xzzz = pbuffer.data(idx_hg + 234);

    auto ts_yyyyy_yyyy = pbuffer.data(idx_hg + 235);

    auto ts_yyyyy_yyyz = pbuffer.data(idx_hg + 236);

    auto ts_yyyyy_yyzz = pbuffer.data(idx_hg + 237);

    auto ts_yyyyy_yzzz = pbuffer.data(idx_hg + 238);

    auto ts_yyyyy_zzzz = pbuffer.data(idx_hg + 239);

    #pragma omp simd aligned(ga_y, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxzz, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyzz, ts_yyy_xzzz, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyzz, ts_yyy_yzzz, ts_yyy_zzzz, ts_yyyy_xxx, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxz, ts_yyyy_xxzz, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_xzz, ts_yyyy_xzzz, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyy_zzz, ts_yyyy_zzzz, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxxz, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xxzz, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyzz, ts_yyyyy_xzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyzz, ts_yyyyy_yzzz, ts_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyyy_xxxx[i] = 4.0 * ts_yyy_xxxx[i] * gfe_0 + ts_yyyy_xxxx[i] * ga_y[i];

        ts_yyyyy_xxxy[i] = 4.0 * ts_yyy_xxxy[i] * gfe_0 + ts_yyyy_xxx[i] * gfe_0 + ts_yyyy_xxxy[i] * ga_y[i];

        ts_yyyyy_xxxz[i] = 4.0 * ts_yyy_xxxz[i] * gfe_0 + ts_yyyy_xxxz[i] * ga_y[i];

        ts_yyyyy_xxyy[i] = 4.0 * ts_yyy_xxyy[i] * gfe_0 + 2.0 * ts_yyyy_xxy[i] * gfe_0 + ts_yyyy_xxyy[i] * ga_y[i];

        ts_yyyyy_xxyz[i] = 4.0 * ts_yyy_xxyz[i] * gfe_0 + ts_yyyy_xxz[i] * gfe_0 + ts_yyyy_xxyz[i] * ga_y[i];

        ts_yyyyy_xxzz[i] = 4.0 * ts_yyy_xxzz[i] * gfe_0 + ts_yyyy_xxzz[i] * ga_y[i];

        ts_yyyyy_xyyy[i] = 4.0 * ts_yyy_xyyy[i] * gfe_0 + 3.0 * ts_yyyy_xyy[i] * gfe_0 + ts_yyyy_xyyy[i] * ga_y[i];

        ts_yyyyy_xyyz[i] = 4.0 * ts_yyy_xyyz[i] * gfe_0 + 2.0 * ts_yyyy_xyz[i] * gfe_0 + ts_yyyy_xyyz[i] * ga_y[i];

        ts_yyyyy_xyzz[i] = 4.0 * ts_yyy_xyzz[i] * gfe_0 + ts_yyyy_xzz[i] * gfe_0 + ts_yyyy_xyzz[i] * ga_y[i];

        ts_yyyyy_xzzz[i] = 4.0 * ts_yyy_xzzz[i] * gfe_0 + ts_yyyy_xzzz[i] * ga_y[i];

        ts_yyyyy_yyyy[i] = 4.0 * ts_yyy_yyyy[i] * gfe_0 + 4.0 * ts_yyyy_yyy[i] * gfe_0 + ts_yyyy_yyyy[i] * ga_y[i];

        ts_yyyyy_yyyz[i] = 4.0 * ts_yyy_yyyz[i] * gfe_0 + 3.0 * ts_yyyy_yyz[i] * gfe_0 + ts_yyyy_yyyz[i] * ga_y[i];

        ts_yyyyy_yyzz[i] = 4.0 * ts_yyy_yyzz[i] * gfe_0 + 2.0 * ts_yyyy_yzz[i] * gfe_0 + ts_yyyy_yyzz[i] * ga_y[i];

        ts_yyyyy_yzzz[i] = 4.0 * ts_yyy_yzzz[i] * gfe_0 + ts_yyyy_zzz[i] * gfe_0 + ts_yyyy_yzzz[i] * ga_y[i];

        ts_yyyyy_zzzz[i] = 4.0 * ts_yyy_zzzz[i] * gfe_0 + ts_yyyy_zzzz[i] * ga_y[i];
    }

    // Set up 240-255 components of targeted buffer : HG

    auto ts_yyyyz_xxxx = pbuffer.data(idx_hg + 240);

    auto ts_yyyyz_xxxy = pbuffer.data(idx_hg + 241);

    auto ts_yyyyz_xxxz = pbuffer.data(idx_hg + 242);

    auto ts_yyyyz_xxyy = pbuffer.data(idx_hg + 243);

    auto ts_yyyyz_xxyz = pbuffer.data(idx_hg + 244);

    auto ts_yyyyz_xxzz = pbuffer.data(idx_hg + 245);

    auto ts_yyyyz_xyyy = pbuffer.data(idx_hg + 246);

    auto ts_yyyyz_xyyz = pbuffer.data(idx_hg + 247);

    auto ts_yyyyz_xyzz = pbuffer.data(idx_hg + 248);

    auto ts_yyyyz_xzzz = pbuffer.data(idx_hg + 249);

    auto ts_yyyyz_yyyy = pbuffer.data(idx_hg + 250);

    auto ts_yyyyz_yyyz = pbuffer.data(idx_hg + 251);

    auto ts_yyyyz_yyzz = pbuffer.data(idx_hg + 252);

    auto ts_yyyyz_yzzz = pbuffer.data(idx_hg + 253);

    auto ts_yyyyz_zzzz = pbuffer.data(idx_hg + 254);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xyy, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyz, ts_yyyy_xyzz, ts_yyyy_yyy, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyz, ts_yyyy_yyzz, ts_yyyy_yzz, ts_yyyy_yzzz, ts_yyyyz_xxxx, ts_yyyyz_xxxy, ts_yyyyz_xxxz, ts_yyyyz_xxyy, ts_yyyyz_xxyz, ts_yyyyz_xxzz, ts_yyyyz_xyyy, ts_yyyyz_xyyz, ts_yyyyz_xyzz, ts_yyyyz_xzzz, ts_yyyyz_yyyy, ts_yyyyz_yyyz, ts_yyyyz_yyzz, ts_yyyyz_yzzz, ts_yyyyz_zzzz, ts_yyyz_xxxz, ts_yyyz_xxzz, ts_yyyz_xzzz, ts_yyyz_zzzz, ts_yyz_xxxz, ts_yyz_xxzz, ts_yyz_xzzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyyz_xxxx[i] = ts_yyyy_xxxx[i] * ga_z[i];

        ts_yyyyz_xxxy[i] = ts_yyyy_xxxy[i] * ga_z[i];

        ts_yyyyz_xxxz[i] = 3.0 * ts_yyz_xxxz[i] * gfe_0 + ts_yyyz_xxxz[i] * ga_y[i];

        ts_yyyyz_xxyy[i] = ts_yyyy_xxyy[i] * ga_z[i];

        ts_yyyyz_xxyz[i] = ts_yyyy_xxy[i] * gfe_0 + ts_yyyy_xxyz[i] * ga_z[i];

        ts_yyyyz_xxzz[i] = 3.0 * ts_yyz_xxzz[i] * gfe_0 + ts_yyyz_xxzz[i] * ga_y[i];

        ts_yyyyz_xyyy[i] = ts_yyyy_xyyy[i] * ga_z[i];

        ts_yyyyz_xyyz[i] = ts_yyyy_xyy[i] * gfe_0 + ts_yyyy_xyyz[i] * ga_z[i];

        ts_yyyyz_xyzz[i] = 2.0 * ts_yyyy_xyz[i] * gfe_0 + ts_yyyy_xyzz[i] * ga_z[i];

        ts_yyyyz_xzzz[i] = 3.0 * ts_yyz_xzzz[i] * gfe_0 + ts_yyyz_xzzz[i] * ga_y[i];

        ts_yyyyz_yyyy[i] = ts_yyyy_yyyy[i] * ga_z[i];

        ts_yyyyz_yyyz[i] = ts_yyyy_yyy[i] * gfe_0 + ts_yyyy_yyyz[i] * ga_z[i];

        ts_yyyyz_yyzz[i] = 2.0 * ts_yyyy_yyz[i] * gfe_0 + ts_yyyy_yyzz[i] * ga_z[i];

        ts_yyyyz_yzzz[i] = 3.0 * ts_yyyy_yzz[i] * gfe_0 + ts_yyyy_yzzz[i] * ga_z[i];

        ts_yyyyz_zzzz[i] = 3.0 * ts_yyz_zzzz[i] * gfe_0 + ts_yyyz_zzzz[i] * ga_y[i];
    }

    // Set up 255-270 components of targeted buffer : HG

    auto ts_yyyzz_xxxx = pbuffer.data(idx_hg + 255);

    auto ts_yyyzz_xxxy = pbuffer.data(idx_hg + 256);

    auto ts_yyyzz_xxxz = pbuffer.data(idx_hg + 257);

    auto ts_yyyzz_xxyy = pbuffer.data(idx_hg + 258);

    auto ts_yyyzz_xxyz = pbuffer.data(idx_hg + 259);

    auto ts_yyyzz_xxzz = pbuffer.data(idx_hg + 260);

    auto ts_yyyzz_xyyy = pbuffer.data(idx_hg + 261);

    auto ts_yyyzz_xyyz = pbuffer.data(idx_hg + 262);

    auto ts_yyyzz_xyzz = pbuffer.data(idx_hg + 263);

    auto ts_yyyzz_xzzz = pbuffer.data(idx_hg + 264);

    auto ts_yyyzz_yyyy = pbuffer.data(idx_hg + 265);

    auto ts_yyyzz_yyyz = pbuffer.data(idx_hg + 266);

    auto ts_yyyzz_yyzz = pbuffer.data(idx_hg + 267);

    auto ts_yyyzz_yzzz = pbuffer.data(idx_hg + 268);

    auto ts_yyyzz_zzzz = pbuffer.data(idx_hg + 269);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyy_xxxy, ts_yyy_xxyy, ts_yyy_xyyy, ts_yyy_yyyy, ts_yyyz_xxxy, ts_yyyz_xxyy, ts_yyyz_xyyy, ts_yyyz_yyyy, ts_yyyzz_xxxx, ts_yyyzz_xxxy, ts_yyyzz_xxxz, ts_yyyzz_xxyy, ts_yyyzz_xxyz, ts_yyyzz_xxzz, ts_yyyzz_xyyy, ts_yyyzz_xyyz, ts_yyyzz_xyzz, ts_yyyzz_xzzz, ts_yyyzz_yyyy, ts_yyyzz_yyyz, ts_yyyzz_yyzz, ts_yyyzz_yzzz, ts_yyyzz_zzzz, ts_yyzz_xxxx, ts_yyzz_xxxz, ts_yyzz_xxyz, ts_yyzz_xxz, ts_yyzz_xxzz, ts_yyzz_xyyz, ts_yyzz_xyz, ts_yyzz_xyzz, ts_yyzz_xzz, ts_yyzz_xzzz, ts_yyzz_yyyz, ts_yyzz_yyz, ts_yyzz_yyzz, ts_yyzz_yzz, ts_yyzz_yzzz, ts_yyzz_zzz, ts_yyzz_zzzz, ts_yzz_xxxx, ts_yzz_xxxz, ts_yzz_xxyz, ts_yzz_xxzz, ts_yzz_xyyz, ts_yzz_xyzz, ts_yzz_xzzz, ts_yzz_yyyz, ts_yzz_yyzz, ts_yzz_yzzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyzz_xxxx[i] = 2.0 * ts_yzz_xxxx[i] * gfe_0 + ts_yyzz_xxxx[i] * ga_y[i];

        ts_yyyzz_xxxy[i] = ts_yyy_xxxy[i] * gfe_0 + ts_yyyz_xxxy[i] * ga_z[i];

        ts_yyyzz_xxxz[i] = 2.0 * ts_yzz_xxxz[i] * gfe_0 + ts_yyzz_xxxz[i] * ga_y[i];

        ts_yyyzz_xxyy[i] = ts_yyy_xxyy[i] * gfe_0 + ts_yyyz_xxyy[i] * ga_z[i];

        ts_yyyzz_xxyz[i] = 2.0 * ts_yzz_xxyz[i] * gfe_0 + ts_yyzz_xxz[i] * gfe_0 + ts_yyzz_xxyz[i] * ga_y[i];

        ts_yyyzz_xxzz[i] = 2.0 * ts_yzz_xxzz[i] * gfe_0 + ts_yyzz_xxzz[i] * ga_y[i];

        ts_yyyzz_xyyy[i] = ts_yyy_xyyy[i] * gfe_0 + ts_yyyz_xyyy[i] * ga_z[i];

        ts_yyyzz_xyyz[i] = 2.0 * ts_yzz_xyyz[i] * gfe_0 + 2.0 * ts_yyzz_xyz[i] * gfe_0 + ts_yyzz_xyyz[i] * ga_y[i];

        ts_yyyzz_xyzz[i] = 2.0 * ts_yzz_xyzz[i] * gfe_0 + ts_yyzz_xzz[i] * gfe_0 + ts_yyzz_xyzz[i] * ga_y[i];

        ts_yyyzz_xzzz[i] = 2.0 * ts_yzz_xzzz[i] * gfe_0 + ts_yyzz_xzzz[i] * ga_y[i];

        ts_yyyzz_yyyy[i] = ts_yyy_yyyy[i] * gfe_0 + ts_yyyz_yyyy[i] * ga_z[i];

        ts_yyyzz_yyyz[i] = 2.0 * ts_yzz_yyyz[i] * gfe_0 + 3.0 * ts_yyzz_yyz[i] * gfe_0 + ts_yyzz_yyyz[i] * ga_y[i];

        ts_yyyzz_yyzz[i] = 2.0 * ts_yzz_yyzz[i] * gfe_0 + 2.0 * ts_yyzz_yzz[i] * gfe_0 + ts_yyzz_yyzz[i] * ga_y[i];

        ts_yyyzz_yzzz[i] = 2.0 * ts_yzz_yzzz[i] * gfe_0 + ts_yyzz_zzz[i] * gfe_0 + ts_yyzz_yzzz[i] * ga_y[i];

        ts_yyyzz_zzzz[i] = 2.0 * ts_yzz_zzzz[i] * gfe_0 + ts_yyzz_zzzz[i] * ga_y[i];
    }

    // Set up 270-285 components of targeted buffer : HG

    auto ts_yyzzz_xxxx = pbuffer.data(idx_hg + 270);

    auto ts_yyzzz_xxxy = pbuffer.data(idx_hg + 271);

    auto ts_yyzzz_xxxz = pbuffer.data(idx_hg + 272);

    auto ts_yyzzz_xxyy = pbuffer.data(idx_hg + 273);

    auto ts_yyzzz_xxyz = pbuffer.data(idx_hg + 274);

    auto ts_yyzzz_xxzz = pbuffer.data(idx_hg + 275);

    auto ts_yyzzz_xyyy = pbuffer.data(idx_hg + 276);

    auto ts_yyzzz_xyyz = pbuffer.data(idx_hg + 277);

    auto ts_yyzzz_xyzz = pbuffer.data(idx_hg + 278);

    auto ts_yyzzz_xzzz = pbuffer.data(idx_hg + 279);

    auto ts_yyzzz_yyyy = pbuffer.data(idx_hg + 280);

    auto ts_yyzzz_yyyz = pbuffer.data(idx_hg + 281);

    auto ts_yyzzz_yyzz = pbuffer.data(idx_hg + 282);

    auto ts_yyzzz_yzzz = pbuffer.data(idx_hg + 283);

    auto ts_yyzzz_zzzz = pbuffer.data(idx_hg + 284);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyz_xxxy, ts_yyz_xxyy, ts_yyz_xyyy, ts_yyz_yyyy, ts_yyzz_xxxy, ts_yyzz_xxyy, ts_yyzz_xyyy, ts_yyzz_yyyy, ts_yyzzz_xxxx, ts_yyzzz_xxxy, ts_yyzzz_xxxz, ts_yyzzz_xxyy, ts_yyzzz_xxyz, ts_yyzzz_xxzz, ts_yyzzz_xyyy, ts_yyzzz_xyyz, ts_yyzzz_xyzz, ts_yyzzz_xzzz, ts_yyzzz_yyyy, ts_yyzzz_yyyz, ts_yyzzz_yyzz, ts_yyzzz_yzzz, ts_yyzzz_zzzz, ts_yzzz_xxxx, ts_yzzz_xxxz, ts_yzzz_xxyz, ts_yzzz_xxz, ts_yzzz_xxzz, ts_yzzz_xyyz, ts_yzzz_xyz, ts_yzzz_xyzz, ts_yzzz_xzz, ts_yzzz_xzzz, ts_yzzz_yyyz, ts_yzzz_yyz, ts_yzzz_yyzz, ts_yzzz_yzz, ts_yzzz_yzzz, ts_yzzz_zzz, ts_yzzz_zzzz, ts_zzz_xxxx, ts_zzz_xxxz, ts_zzz_xxyz, ts_zzz_xxzz, ts_zzz_xyyz, ts_zzz_xyzz, ts_zzz_xzzz, ts_zzz_yyyz, ts_zzz_yyzz, ts_zzz_yzzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyzzz_xxxx[i] = ts_zzz_xxxx[i] * gfe_0 + ts_yzzz_xxxx[i] * ga_y[i];

        ts_yyzzz_xxxy[i] = 2.0 * ts_yyz_xxxy[i] * gfe_0 + ts_yyzz_xxxy[i] * ga_z[i];

        ts_yyzzz_xxxz[i] = ts_zzz_xxxz[i] * gfe_0 + ts_yzzz_xxxz[i] * ga_y[i];

        ts_yyzzz_xxyy[i] = 2.0 * ts_yyz_xxyy[i] * gfe_0 + ts_yyzz_xxyy[i] * ga_z[i];

        ts_yyzzz_xxyz[i] = ts_zzz_xxyz[i] * gfe_0 + ts_yzzz_xxz[i] * gfe_0 + ts_yzzz_xxyz[i] * ga_y[i];

        ts_yyzzz_xxzz[i] = ts_zzz_xxzz[i] * gfe_0 + ts_yzzz_xxzz[i] * ga_y[i];

        ts_yyzzz_xyyy[i] = 2.0 * ts_yyz_xyyy[i] * gfe_0 + ts_yyzz_xyyy[i] * ga_z[i];

        ts_yyzzz_xyyz[i] = ts_zzz_xyyz[i] * gfe_0 + 2.0 * ts_yzzz_xyz[i] * gfe_0 + ts_yzzz_xyyz[i] * ga_y[i];

        ts_yyzzz_xyzz[i] = ts_zzz_xyzz[i] * gfe_0 + ts_yzzz_xzz[i] * gfe_0 + ts_yzzz_xyzz[i] * ga_y[i];

        ts_yyzzz_xzzz[i] = ts_zzz_xzzz[i] * gfe_0 + ts_yzzz_xzzz[i] * ga_y[i];

        ts_yyzzz_yyyy[i] = 2.0 * ts_yyz_yyyy[i] * gfe_0 + ts_yyzz_yyyy[i] * ga_z[i];

        ts_yyzzz_yyyz[i] = ts_zzz_yyyz[i] * gfe_0 + 3.0 * ts_yzzz_yyz[i] * gfe_0 + ts_yzzz_yyyz[i] * ga_y[i];

        ts_yyzzz_yyzz[i] = ts_zzz_yyzz[i] * gfe_0 + 2.0 * ts_yzzz_yzz[i] * gfe_0 + ts_yzzz_yyzz[i] * ga_y[i];

        ts_yyzzz_yzzz[i] = ts_zzz_yzzz[i] * gfe_0 + ts_yzzz_zzz[i] * gfe_0 + ts_yzzz_yzzz[i] * ga_y[i];

        ts_yyzzz_zzzz[i] = ts_zzz_zzzz[i] * gfe_0 + ts_yzzz_zzzz[i] * ga_y[i];
    }

    // Set up 285-300 components of targeted buffer : HG

    auto ts_yzzzz_xxxx = pbuffer.data(idx_hg + 285);

    auto ts_yzzzz_xxxy = pbuffer.data(idx_hg + 286);

    auto ts_yzzzz_xxxz = pbuffer.data(idx_hg + 287);

    auto ts_yzzzz_xxyy = pbuffer.data(idx_hg + 288);

    auto ts_yzzzz_xxyz = pbuffer.data(idx_hg + 289);

    auto ts_yzzzz_xxzz = pbuffer.data(idx_hg + 290);

    auto ts_yzzzz_xyyy = pbuffer.data(idx_hg + 291);

    auto ts_yzzzz_xyyz = pbuffer.data(idx_hg + 292);

    auto ts_yzzzz_xyzz = pbuffer.data(idx_hg + 293);

    auto ts_yzzzz_xzzz = pbuffer.data(idx_hg + 294);

    auto ts_yzzzz_yyyy = pbuffer.data(idx_hg + 295);

    auto ts_yzzzz_yyyz = pbuffer.data(idx_hg + 296);

    auto ts_yzzzz_yyzz = pbuffer.data(idx_hg + 297);

    auto ts_yzzzz_yzzz = pbuffer.data(idx_hg + 298);

    auto ts_yzzzz_zzzz = pbuffer.data(idx_hg + 299);

    #pragma omp simd aligned(ga_y, ts_yzzzz_xxxx, ts_yzzzz_xxxy, ts_yzzzz_xxxz, ts_yzzzz_xxyy, ts_yzzzz_xxyz, ts_yzzzz_xxzz, ts_yzzzz_xyyy, ts_yzzzz_xyyz, ts_yzzzz_xyzz, ts_yzzzz_xzzz, ts_yzzzz_yyyy, ts_yzzzz_yyyz, ts_yzzzz_yyzz, ts_yzzzz_yzzz, ts_yzzzz_zzzz, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yzzzz_xxxx[i] = ts_zzzz_xxxx[i] * ga_y[i];

        ts_yzzzz_xxxy[i] = ts_zzzz_xxx[i] * gfe_0 + ts_zzzz_xxxy[i] * ga_y[i];

        ts_yzzzz_xxxz[i] = ts_zzzz_xxxz[i] * ga_y[i];

        ts_yzzzz_xxyy[i] = 2.0 * ts_zzzz_xxy[i] * gfe_0 + ts_zzzz_xxyy[i] * ga_y[i];

        ts_yzzzz_xxyz[i] = ts_zzzz_xxz[i] * gfe_0 + ts_zzzz_xxyz[i] * ga_y[i];

        ts_yzzzz_xxzz[i] = ts_zzzz_xxzz[i] * ga_y[i];

        ts_yzzzz_xyyy[i] = 3.0 * ts_zzzz_xyy[i] * gfe_0 + ts_zzzz_xyyy[i] * ga_y[i];

        ts_yzzzz_xyyz[i] = 2.0 * ts_zzzz_xyz[i] * gfe_0 + ts_zzzz_xyyz[i] * ga_y[i];

        ts_yzzzz_xyzz[i] = ts_zzzz_xzz[i] * gfe_0 + ts_zzzz_xyzz[i] * ga_y[i];

        ts_yzzzz_xzzz[i] = ts_zzzz_xzzz[i] * ga_y[i];

        ts_yzzzz_yyyy[i] = 4.0 * ts_zzzz_yyy[i] * gfe_0 + ts_zzzz_yyyy[i] * ga_y[i];

        ts_yzzzz_yyyz[i] = 3.0 * ts_zzzz_yyz[i] * gfe_0 + ts_zzzz_yyyz[i] * ga_y[i];

        ts_yzzzz_yyzz[i] = 2.0 * ts_zzzz_yzz[i] * gfe_0 + ts_zzzz_yyzz[i] * ga_y[i];

        ts_yzzzz_yzzz[i] = ts_zzzz_zzz[i] * gfe_0 + ts_zzzz_yzzz[i] * ga_y[i];

        ts_yzzzz_zzzz[i] = ts_zzzz_zzzz[i] * ga_y[i];
    }

    // Set up 300-315 components of targeted buffer : HG

    auto ts_zzzzz_xxxx = pbuffer.data(idx_hg + 300);

    auto ts_zzzzz_xxxy = pbuffer.data(idx_hg + 301);

    auto ts_zzzzz_xxxz = pbuffer.data(idx_hg + 302);

    auto ts_zzzzz_xxyy = pbuffer.data(idx_hg + 303);

    auto ts_zzzzz_xxyz = pbuffer.data(idx_hg + 304);

    auto ts_zzzzz_xxzz = pbuffer.data(idx_hg + 305);

    auto ts_zzzzz_xyyy = pbuffer.data(idx_hg + 306);

    auto ts_zzzzz_xyyz = pbuffer.data(idx_hg + 307);

    auto ts_zzzzz_xyzz = pbuffer.data(idx_hg + 308);

    auto ts_zzzzz_xzzz = pbuffer.data(idx_hg + 309);

    auto ts_zzzzz_yyyy = pbuffer.data(idx_hg + 310);

    auto ts_zzzzz_yyyz = pbuffer.data(idx_hg + 311);

    auto ts_zzzzz_yyzz = pbuffer.data(idx_hg + 312);

    auto ts_zzzzz_yzzz = pbuffer.data(idx_hg + 313);

    auto ts_zzzzz_zzzz = pbuffer.data(idx_hg + 314);

    #pragma omp simd aligned(ga_z, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxzz, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyzz, ts_zzz_xzzz, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyzz, ts_zzz_yzzz, ts_zzz_zzzz, ts_zzzz_xxx, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxy, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxz, ts_zzzz_xxzz, ts_zzzz_xyy, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyz, ts_zzzz_xyzz, ts_zzzz_xzz, ts_zzzz_xzzz, ts_zzzz_yyy, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyz, ts_zzzz_yyzz, ts_zzzz_yzz, ts_zzzz_yzzz, ts_zzzz_zzz, ts_zzzz_zzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxzz, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyzz, ts_zzzzz_xzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyzz, ts_zzzzz_yzzz, ts_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zzzzz_xxxx[i] = 4.0 * ts_zzz_xxxx[i] * gfe_0 + ts_zzzz_xxxx[i] * ga_z[i];

        ts_zzzzz_xxxy[i] = 4.0 * ts_zzz_xxxy[i] * gfe_0 + ts_zzzz_xxxy[i] * ga_z[i];

        ts_zzzzz_xxxz[i] = 4.0 * ts_zzz_xxxz[i] * gfe_0 + ts_zzzz_xxx[i] * gfe_0 + ts_zzzz_xxxz[i] * ga_z[i];

        ts_zzzzz_xxyy[i] = 4.0 * ts_zzz_xxyy[i] * gfe_0 + ts_zzzz_xxyy[i] * ga_z[i];

        ts_zzzzz_xxyz[i] = 4.0 * ts_zzz_xxyz[i] * gfe_0 + ts_zzzz_xxy[i] * gfe_0 + ts_zzzz_xxyz[i] * ga_z[i];

        ts_zzzzz_xxzz[i] = 4.0 * ts_zzz_xxzz[i] * gfe_0 + 2.0 * ts_zzzz_xxz[i] * gfe_0 + ts_zzzz_xxzz[i] * ga_z[i];

        ts_zzzzz_xyyy[i] = 4.0 * ts_zzz_xyyy[i] * gfe_0 + ts_zzzz_xyyy[i] * ga_z[i];

        ts_zzzzz_xyyz[i] = 4.0 * ts_zzz_xyyz[i] * gfe_0 + ts_zzzz_xyy[i] * gfe_0 + ts_zzzz_xyyz[i] * ga_z[i];

        ts_zzzzz_xyzz[i] = 4.0 * ts_zzz_xyzz[i] * gfe_0 + 2.0 * ts_zzzz_xyz[i] * gfe_0 + ts_zzzz_xyzz[i] * ga_z[i];

        ts_zzzzz_xzzz[i] = 4.0 * ts_zzz_xzzz[i] * gfe_0 + 3.0 * ts_zzzz_xzz[i] * gfe_0 + ts_zzzz_xzzz[i] * ga_z[i];

        ts_zzzzz_yyyy[i] = 4.0 * ts_zzz_yyyy[i] * gfe_0 + ts_zzzz_yyyy[i] * ga_z[i];

        ts_zzzzz_yyyz[i] = 4.0 * ts_zzz_yyyz[i] * gfe_0 + ts_zzzz_yyy[i] * gfe_0 + ts_zzzz_yyyz[i] * ga_z[i];

        ts_zzzzz_yyzz[i] = 4.0 * ts_zzz_yyzz[i] * gfe_0 + 2.0 * ts_zzzz_yyz[i] * gfe_0 + ts_zzzz_yyzz[i] * ga_z[i];

        ts_zzzzz_yzzz[i] = 4.0 * ts_zzz_yzzz[i] * gfe_0 + 3.0 * ts_zzzz_yzz[i] * gfe_0 + ts_zzzz_yzzz[i] * ga_z[i];

        ts_zzzzz_zzzz[i] = 4.0 * ts_zzz_zzzz[i] * gfe_0 + 4.0 * ts_zzzz_zzz[i] * gfe_0 + ts_zzzz_zzzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

