#include "LocalCorePotentialPrimRecHG.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_hg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hg,
                                  const size_t idx_fg,
                                  const size_t idx_gf,
                                  const size_t idx_gg,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : FG

    auto tg_xxx_xxxx = pbuffer.data(idx_fg);

    auto tg_xxx_xxxy = pbuffer.data(idx_fg + 1);

    auto tg_xxx_xxxz = pbuffer.data(idx_fg + 2);

    auto tg_xxx_xxyy = pbuffer.data(idx_fg + 3);

    auto tg_xxx_xxyz = pbuffer.data(idx_fg + 4);

    auto tg_xxx_xxzz = pbuffer.data(idx_fg + 5);

    auto tg_xxx_xyyy = pbuffer.data(idx_fg + 6);

    auto tg_xxx_xyyz = pbuffer.data(idx_fg + 7);

    auto tg_xxx_xyzz = pbuffer.data(idx_fg + 8);

    auto tg_xxx_xzzz = pbuffer.data(idx_fg + 9);

    auto tg_xxx_yyyy = pbuffer.data(idx_fg + 10);

    auto tg_xxx_yyyz = pbuffer.data(idx_fg + 11);

    auto tg_xxx_yyzz = pbuffer.data(idx_fg + 12);

    auto tg_xxx_yzzz = pbuffer.data(idx_fg + 13);

    auto tg_xxx_zzzz = pbuffer.data(idx_fg + 14);

    auto tg_xxy_xxxx = pbuffer.data(idx_fg + 15);

    auto tg_xxy_xxxz = pbuffer.data(idx_fg + 17);

    auto tg_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto tg_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto tg_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto tg_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto tg_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto tg_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto tg_xxz_xxxx = pbuffer.data(idx_fg + 30);

    auto tg_xxz_xxxy = pbuffer.data(idx_fg + 31);

    auto tg_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto tg_xxz_xxyy = pbuffer.data(idx_fg + 33);

    auto tg_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto tg_xxz_xyyy = pbuffer.data(idx_fg + 36);

    auto tg_xxz_xzzz = pbuffer.data(idx_fg + 39);

    auto tg_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto tg_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto tg_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto tg_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto tg_xyy_xxxy = pbuffer.data(idx_fg + 46);

    auto tg_xyy_xxyy = pbuffer.data(idx_fg + 48);

    auto tg_xyy_xxyz = pbuffer.data(idx_fg + 49);

    auto tg_xyy_xyyy = pbuffer.data(idx_fg + 51);

    auto tg_xyy_xyyz = pbuffer.data(idx_fg + 52);

    auto tg_xyy_xyzz = pbuffer.data(idx_fg + 53);

    auto tg_xyy_yyyy = pbuffer.data(idx_fg + 55);

    auto tg_xyy_yyyz = pbuffer.data(idx_fg + 56);

    auto tg_xyy_yyzz = pbuffer.data(idx_fg + 57);

    auto tg_xyy_yzzz = pbuffer.data(idx_fg + 58);

    auto tg_xyy_zzzz = pbuffer.data(idx_fg + 59);

    auto tg_xyz_yyyz = pbuffer.data(idx_fg + 71);

    auto tg_xyz_yyzz = pbuffer.data(idx_fg + 72);

    auto tg_xyz_yzzz = pbuffer.data(idx_fg + 73);

    auto tg_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto tg_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto tg_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto tg_xzz_xyyz = pbuffer.data(idx_fg + 82);

    auto tg_xzz_xyzz = pbuffer.data(idx_fg + 83);

    auto tg_xzz_xzzz = pbuffer.data(idx_fg + 84);

    auto tg_xzz_yyyy = pbuffer.data(idx_fg + 85);

    auto tg_xzz_yyyz = pbuffer.data(idx_fg + 86);

    auto tg_xzz_yyzz = pbuffer.data(idx_fg + 87);

    auto tg_xzz_yzzz = pbuffer.data(idx_fg + 88);

    auto tg_xzz_zzzz = pbuffer.data(idx_fg + 89);

    auto tg_yyy_xxxx = pbuffer.data(idx_fg + 90);

    auto tg_yyy_xxxy = pbuffer.data(idx_fg + 91);

    auto tg_yyy_xxxz = pbuffer.data(idx_fg + 92);

    auto tg_yyy_xxyy = pbuffer.data(idx_fg + 93);

    auto tg_yyy_xxyz = pbuffer.data(idx_fg + 94);

    auto tg_yyy_xxzz = pbuffer.data(idx_fg + 95);

    auto tg_yyy_xyyy = pbuffer.data(idx_fg + 96);

    auto tg_yyy_xyyz = pbuffer.data(idx_fg + 97);

    auto tg_yyy_xyzz = pbuffer.data(idx_fg + 98);

    auto tg_yyy_xzzz = pbuffer.data(idx_fg + 99);

    auto tg_yyy_yyyy = pbuffer.data(idx_fg + 100);

    auto tg_yyy_yyyz = pbuffer.data(idx_fg + 101);

    auto tg_yyy_yyzz = pbuffer.data(idx_fg + 102);

    auto tg_yyy_yzzz = pbuffer.data(idx_fg + 103);

    auto tg_yyy_zzzz = pbuffer.data(idx_fg + 104);

    auto tg_yyz_xxxy = pbuffer.data(idx_fg + 106);

    auto tg_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto tg_yyz_xxyy = pbuffer.data(idx_fg + 108);

    auto tg_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto tg_yyz_xyyy = pbuffer.data(idx_fg + 111);

    auto tg_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto tg_yyz_yyyy = pbuffer.data(idx_fg + 115);

    auto tg_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto tg_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto tg_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto tg_yyz_zzzz = pbuffer.data(idx_fg + 119);

    auto tg_yzz_xxxx = pbuffer.data(idx_fg + 120);

    auto tg_yzz_xxxz = pbuffer.data(idx_fg + 122);

    auto tg_yzz_xxyz = pbuffer.data(idx_fg + 124);

    auto tg_yzz_xxzz = pbuffer.data(idx_fg + 125);

    auto tg_yzz_xyyz = pbuffer.data(idx_fg + 127);

    auto tg_yzz_xyzz = pbuffer.data(idx_fg + 128);

    auto tg_yzz_xzzz = pbuffer.data(idx_fg + 129);

    auto tg_yzz_yyyy = pbuffer.data(idx_fg + 130);

    auto tg_yzz_yyyz = pbuffer.data(idx_fg + 131);

    auto tg_yzz_yyzz = pbuffer.data(idx_fg + 132);

    auto tg_yzz_yzzz = pbuffer.data(idx_fg + 133);

    auto tg_yzz_zzzz = pbuffer.data(idx_fg + 134);

    auto tg_zzz_xxxx = pbuffer.data(idx_fg + 135);

    auto tg_zzz_xxxy = pbuffer.data(idx_fg + 136);

    auto tg_zzz_xxxz = pbuffer.data(idx_fg + 137);

    auto tg_zzz_xxyy = pbuffer.data(idx_fg + 138);

    auto tg_zzz_xxyz = pbuffer.data(idx_fg + 139);

    auto tg_zzz_xxzz = pbuffer.data(idx_fg + 140);

    auto tg_zzz_xyyy = pbuffer.data(idx_fg + 141);

    auto tg_zzz_xyyz = pbuffer.data(idx_fg + 142);

    auto tg_zzz_xyzz = pbuffer.data(idx_fg + 143);

    auto tg_zzz_xzzz = pbuffer.data(idx_fg + 144);

    auto tg_zzz_yyyy = pbuffer.data(idx_fg + 145);

    auto tg_zzz_yyyz = pbuffer.data(idx_fg + 146);

    auto tg_zzz_yyzz = pbuffer.data(idx_fg + 147);

    auto tg_zzz_yzzz = pbuffer.data(idx_fg + 148);

    auto tg_zzz_zzzz = pbuffer.data(idx_fg + 149);

    // Set up components of auxiliary buffer : GF

    auto tg_xxxx_xxx = pbuffer.data(idx_gf);

    auto tg_xxxx_xxy = pbuffer.data(idx_gf + 1);

    auto tg_xxxx_xxz = pbuffer.data(idx_gf + 2);

    auto tg_xxxx_xyy = pbuffer.data(idx_gf + 3);

    auto tg_xxxx_xyz = pbuffer.data(idx_gf + 4);

    auto tg_xxxx_xzz = pbuffer.data(idx_gf + 5);

    auto tg_xxxx_yyy = pbuffer.data(idx_gf + 6);

    auto tg_xxxx_yyz = pbuffer.data(idx_gf + 7);

    auto tg_xxxx_yzz = pbuffer.data(idx_gf + 8);

    auto tg_xxxx_zzz = pbuffer.data(idx_gf + 9);

    auto tg_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto tg_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto tg_xxxz_xzz = pbuffer.data(idx_gf + 25);

    auto tg_xxyy_xxy = pbuffer.data(idx_gf + 31);

    auto tg_xxyy_xyy = pbuffer.data(idx_gf + 33);

    auto tg_xxyy_xyz = pbuffer.data(idx_gf + 34);

    auto tg_xxyy_yyy = pbuffer.data(idx_gf + 36);

    auto tg_xxyy_yyz = pbuffer.data(idx_gf + 37);

    auto tg_xxyy_yzz = pbuffer.data(idx_gf + 38);

    auto tg_xxzz_xxx = pbuffer.data(idx_gf + 50);

    auto tg_xxzz_xxy = pbuffer.data(idx_gf + 51);

    auto tg_xxzz_xxz = pbuffer.data(idx_gf + 52);

    auto tg_xxzz_xyy = pbuffer.data(idx_gf + 53);

    auto tg_xxzz_xyz = pbuffer.data(idx_gf + 54);

    auto tg_xxzz_xzz = pbuffer.data(idx_gf + 55);

    auto tg_xxzz_yyz = pbuffer.data(idx_gf + 57);

    auto tg_xxzz_yzz = pbuffer.data(idx_gf + 58);

    auto tg_xxzz_zzz = pbuffer.data(idx_gf + 59);

    auto tg_xyyy_xxy = pbuffer.data(idx_gf + 61);

    auto tg_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto tg_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto tg_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto tg_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto tg_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto tg_xzzz_xxz = pbuffer.data(idx_gf + 92);

    auto tg_xzzz_xyz = pbuffer.data(idx_gf + 94);

    auto tg_xzzz_xzz = pbuffer.data(idx_gf + 95);

    auto tg_xzzz_yyz = pbuffer.data(idx_gf + 97);

    auto tg_xzzz_yzz = pbuffer.data(idx_gf + 98);

    auto tg_xzzz_zzz = pbuffer.data(idx_gf + 99);

    auto tg_yyyy_xxx = pbuffer.data(idx_gf + 100);

    auto tg_yyyy_xxy = pbuffer.data(idx_gf + 101);

    auto tg_yyyy_xxz = pbuffer.data(idx_gf + 102);

    auto tg_yyyy_xyy = pbuffer.data(idx_gf + 103);

    auto tg_yyyy_xyz = pbuffer.data(idx_gf + 104);

    auto tg_yyyy_xzz = pbuffer.data(idx_gf + 105);

    auto tg_yyyy_yyy = pbuffer.data(idx_gf + 106);

    auto tg_yyyy_yyz = pbuffer.data(idx_gf + 107);

    auto tg_yyyy_yzz = pbuffer.data(idx_gf + 108);

    auto tg_yyyy_zzz = pbuffer.data(idx_gf + 109);

    auto tg_yyyz_xxz = pbuffer.data(idx_gf + 112);

    auto tg_yyyz_xyz = pbuffer.data(idx_gf + 114);

    auto tg_yyyz_xzz = pbuffer.data(idx_gf + 115);

    auto tg_yyyz_yyz = pbuffer.data(idx_gf + 117);

    auto tg_yyyz_yzz = pbuffer.data(idx_gf + 118);

    auto tg_yyyz_zzz = pbuffer.data(idx_gf + 119);

    auto tg_yyzz_xxx = pbuffer.data(idx_gf + 120);

    auto tg_yyzz_xxy = pbuffer.data(idx_gf + 121);

    auto tg_yyzz_xxz = pbuffer.data(idx_gf + 122);

    auto tg_yyzz_xyy = pbuffer.data(idx_gf + 123);

    auto tg_yyzz_xyz = pbuffer.data(idx_gf + 124);

    auto tg_yyzz_xzz = pbuffer.data(idx_gf + 125);

    auto tg_yyzz_yyy = pbuffer.data(idx_gf + 126);

    auto tg_yyzz_yyz = pbuffer.data(idx_gf + 127);

    auto tg_yyzz_yzz = pbuffer.data(idx_gf + 128);

    auto tg_yyzz_zzz = pbuffer.data(idx_gf + 129);

    auto tg_yzzz_xxy = pbuffer.data(idx_gf + 131);

    auto tg_yzzz_xxz = pbuffer.data(idx_gf + 132);

    auto tg_yzzz_xyy = pbuffer.data(idx_gf + 133);

    auto tg_yzzz_xyz = pbuffer.data(idx_gf + 134);

    auto tg_yzzz_xzz = pbuffer.data(idx_gf + 135);

    auto tg_yzzz_yyy = pbuffer.data(idx_gf + 136);

    auto tg_yzzz_yyz = pbuffer.data(idx_gf + 137);

    auto tg_yzzz_yzz = pbuffer.data(idx_gf + 138);

    auto tg_yzzz_zzz = pbuffer.data(idx_gf + 139);

    auto tg_zzzz_xxx = pbuffer.data(idx_gf + 140);

    auto tg_zzzz_xxy = pbuffer.data(idx_gf + 141);

    auto tg_zzzz_xxz = pbuffer.data(idx_gf + 142);

    auto tg_zzzz_xyy = pbuffer.data(idx_gf + 143);

    auto tg_zzzz_xyz = pbuffer.data(idx_gf + 144);

    auto tg_zzzz_xzz = pbuffer.data(idx_gf + 145);

    auto tg_zzzz_yyy = pbuffer.data(idx_gf + 146);

    auto tg_zzzz_yyz = pbuffer.data(idx_gf + 147);

    auto tg_zzzz_yzz = pbuffer.data(idx_gf + 148);

    auto tg_zzzz_zzz = pbuffer.data(idx_gf + 149);

    // Set up components of auxiliary buffer : GG

    auto tg_xxxx_xxxx = pbuffer.data(idx_gg);

    auto tg_xxxx_xxxy = pbuffer.data(idx_gg + 1);

    auto tg_xxxx_xxxz = pbuffer.data(idx_gg + 2);

    auto tg_xxxx_xxyy = pbuffer.data(idx_gg + 3);

    auto tg_xxxx_xxyz = pbuffer.data(idx_gg + 4);

    auto tg_xxxx_xxzz = pbuffer.data(idx_gg + 5);

    auto tg_xxxx_xyyy = pbuffer.data(idx_gg + 6);

    auto tg_xxxx_xyyz = pbuffer.data(idx_gg + 7);

    auto tg_xxxx_xyzz = pbuffer.data(idx_gg + 8);

    auto tg_xxxx_xzzz = pbuffer.data(idx_gg + 9);

    auto tg_xxxx_yyyy = pbuffer.data(idx_gg + 10);

    auto tg_xxxx_yyyz = pbuffer.data(idx_gg + 11);

    auto tg_xxxx_yyzz = pbuffer.data(idx_gg + 12);

    auto tg_xxxx_yzzz = pbuffer.data(idx_gg + 13);

    auto tg_xxxx_zzzz = pbuffer.data(idx_gg + 14);

    auto tg_xxxy_xxxx = pbuffer.data(idx_gg + 15);

    auto tg_xxxy_xxxy = pbuffer.data(idx_gg + 16);

    auto tg_xxxy_xxxz = pbuffer.data(idx_gg + 17);

    auto tg_xxxy_xxyy = pbuffer.data(idx_gg + 18);

    auto tg_xxxy_xxzz = pbuffer.data(idx_gg + 20);

    auto tg_xxxy_xyyy = pbuffer.data(idx_gg + 21);

    auto tg_xxxy_xzzz = pbuffer.data(idx_gg + 24);

    auto tg_xxxy_yyyy = pbuffer.data(idx_gg + 25);

    auto tg_xxxy_yyyz = pbuffer.data(idx_gg + 26);

    auto tg_xxxy_yyzz = pbuffer.data(idx_gg + 27);

    auto tg_xxxy_yzzz = pbuffer.data(idx_gg + 28);

    auto tg_xxxz_xxxx = pbuffer.data(idx_gg + 30);

    auto tg_xxxz_xxxy = pbuffer.data(idx_gg + 31);

    auto tg_xxxz_xxxz = pbuffer.data(idx_gg + 32);

    auto tg_xxxz_xxyy = pbuffer.data(idx_gg + 33);

    auto tg_xxxz_xxyz = pbuffer.data(idx_gg + 34);

    auto tg_xxxz_xxzz = pbuffer.data(idx_gg + 35);

    auto tg_xxxz_xyyy = pbuffer.data(idx_gg + 36);

    auto tg_xxxz_xyyz = pbuffer.data(idx_gg + 37);

    auto tg_xxxz_xyzz = pbuffer.data(idx_gg + 38);

    auto tg_xxxz_xzzz = pbuffer.data(idx_gg + 39);

    auto tg_xxxz_yyyz = pbuffer.data(idx_gg + 41);

    auto tg_xxxz_yyzz = pbuffer.data(idx_gg + 42);

    auto tg_xxxz_yzzz = pbuffer.data(idx_gg + 43);

    auto tg_xxxz_zzzz = pbuffer.data(idx_gg + 44);

    auto tg_xxyy_xxxx = pbuffer.data(idx_gg + 45);

    auto tg_xxyy_xxxy = pbuffer.data(idx_gg + 46);

    auto tg_xxyy_xxxz = pbuffer.data(idx_gg + 47);

    auto tg_xxyy_xxyy = pbuffer.data(idx_gg + 48);

    auto tg_xxyy_xxyz = pbuffer.data(idx_gg + 49);

    auto tg_xxyy_xxzz = pbuffer.data(idx_gg + 50);

    auto tg_xxyy_xyyy = pbuffer.data(idx_gg + 51);

    auto tg_xxyy_xyyz = pbuffer.data(idx_gg + 52);

    auto tg_xxyy_xyzz = pbuffer.data(idx_gg + 53);

    auto tg_xxyy_xzzz = pbuffer.data(idx_gg + 54);

    auto tg_xxyy_yyyy = pbuffer.data(idx_gg + 55);

    auto tg_xxyy_yyyz = pbuffer.data(idx_gg + 56);

    auto tg_xxyy_yyzz = pbuffer.data(idx_gg + 57);

    auto tg_xxyy_yzzz = pbuffer.data(idx_gg + 58);

    auto tg_xxyy_zzzz = pbuffer.data(idx_gg + 59);

    auto tg_xxyz_xxxz = pbuffer.data(idx_gg + 62);

    auto tg_xxyz_xxzz = pbuffer.data(idx_gg + 65);

    auto tg_xxyz_xzzz = pbuffer.data(idx_gg + 69);

    auto tg_xxyz_yyyz = pbuffer.data(idx_gg + 71);

    auto tg_xxyz_yyzz = pbuffer.data(idx_gg + 72);

    auto tg_xxyz_yzzz = pbuffer.data(idx_gg + 73);

    auto tg_xxzz_xxxx = pbuffer.data(idx_gg + 75);

    auto tg_xxzz_xxxy = pbuffer.data(idx_gg + 76);

    auto tg_xxzz_xxxz = pbuffer.data(idx_gg + 77);

    auto tg_xxzz_xxyy = pbuffer.data(idx_gg + 78);

    auto tg_xxzz_xxyz = pbuffer.data(idx_gg + 79);

    auto tg_xxzz_xxzz = pbuffer.data(idx_gg + 80);

    auto tg_xxzz_xyyy = pbuffer.data(idx_gg + 81);

    auto tg_xxzz_xyyz = pbuffer.data(idx_gg + 82);

    auto tg_xxzz_xyzz = pbuffer.data(idx_gg + 83);

    auto tg_xxzz_xzzz = pbuffer.data(idx_gg + 84);

    auto tg_xxzz_yyyy = pbuffer.data(idx_gg + 85);

    auto tg_xxzz_yyyz = pbuffer.data(idx_gg + 86);

    auto tg_xxzz_yyzz = pbuffer.data(idx_gg + 87);

    auto tg_xxzz_yzzz = pbuffer.data(idx_gg + 88);

    auto tg_xxzz_zzzz = pbuffer.data(idx_gg + 89);

    auto tg_xyyy_xxxx = pbuffer.data(idx_gg + 90);

    auto tg_xyyy_xxxy = pbuffer.data(idx_gg + 91);

    auto tg_xyyy_xxyy = pbuffer.data(idx_gg + 93);

    auto tg_xyyy_xxyz = pbuffer.data(idx_gg + 94);

    auto tg_xyyy_xyyy = pbuffer.data(idx_gg + 96);

    auto tg_xyyy_xyyz = pbuffer.data(idx_gg + 97);

    auto tg_xyyy_xyzz = pbuffer.data(idx_gg + 98);

    auto tg_xyyy_yyyy = pbuffer.data(idx_gg + 100);

    auto tg_xyyy_yyyz = pbuffer.data(idx_gg + 101);

    auto tg_xyyy_yyzz = pbuffer.data(idx_gg + 102);

    auto tg_xyyy_yzzz = pbuffer.data(idx_gg + 103);

    auto tg_xyyy_zzzz = pbuffer.data(idx_gg + 104);

    auto tg_xyyz_yyyz = pbuffer.data(idx_gg + 116);

    auto tg_xyyz_yyzz = pbuffer.data(idx_gg + 117);

    auto tg_xyyz_yzzz = pbuffer.data(idx_gg + 118);

    auto tg_xyyz_zzzz = pbuffer.data(idx_gg + 119);

    auto tg_xyzz_yyyy = pbuffer.data(idx_gg + 130);

    auto tg_xyzz_yyyz = pbuffer.data(idx_gg + 131);

    auto tg_xyzz_yyzz = pbuffer.data(idx_gg + 132);

    auto tg_xyzz_yzzz = pbuffer.data(idx_gg + 133);

    auto tg_xzzz_xxxx = pbuffer.data(idx_gg + 135);

    auto tg_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto tg_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto tg_xzzz_xxzz = pbuffer.data(idx_gg + 140);

    auto tg_xzzz_xyyz = pbuffer.data(idx_gg + 142);

    auto tg_xzzz_xyzz = pbuffer.data(idx_gg + 143);

    auto tg_xzzz_xzzz = pbuffer.data(idx_gg + 144);

    auto tg_xzzz_yyyy = pbuffer.data(idx_gg + 145);

    auto tg_xzzz_yyyz = pbuffer.data(idx_gg + 146);

    auto tg_xzzz_yyzz = pbuffer.data(idx_gg + 147);

    auto tg_xzzz_yzzz = pbuffer.data(idx_gg + 148);

    auto tg_xzzz_zzzz = pbuffer.data(idx_gg + 149);

    auto tg_yyyy_xxxx = pbuffer.data(idx_gg + 150);

    auto tg_yyyy_xxxy = pbuffer.data(idx_gg + 151);

    auto tg_yyyy_xxxz = pbuffer.data(idx_gg + 152);

    auto tg_yyyy_xxyy = pbuffer.data(idx_gg + 153);

    auto tg_yyyy_xxyz = pbuffer.data(idx_gg + 154);

    auto tg_yyyy_xxzz = pbuffer.data(idx_gg + 155);

    auto tg_yyyy_xyyy = pbuffer.data(idx_gg + 156);

    auto tg_yyyy_xyyz = pbuffer.data(idx_gg + 157);

    auto tg_yyyy_xyzz = pbuffer.data(idx_gg + 158);

    auto tg_yyyy_xzzz = pbuffer.data(idx_gg + 159);

    auto tg_yyyy_yyyy = pbuffer.data(idx_gg + 160);

    auto tg_yyyy_yyyz = pbuffer.data(idx_gg + 161);

    auto tg_yyyy_yyzz = pbuffer.data(idx_gg + 162);

    auto tg_yyyy_yzzz = pbuffer.data(idx_gg + 163);

    auto tg_yyyy_zzzz = pbuffer.data(idx_gg + 164);

    auto tg_yyyz_xxxy = pbuffer.data(idx_gg + 166);

    auto tg_yyyz_xxxz = pbuffer.data(idx_gg + 167);

    auto tg_yyyz_xxyy = pbuffer.data(idx_gg + 168);

    auto tg_yyyz_xxyz = pbuffer.data(idx_gg + 169);

    auto tg_yyyz_xxzz = pbuffer.data(idx_gg + 170);

    auto tg_yyyz_xyyy = pbuffer.data(idx_gg + 171);

    auto tg_yyyz_xyyz = pbuffer.data(idx_gg + 172);

    auto tg_yyyz_xyzz = pbuffer.data(idx_gg + 173);

    auto tg_yyyz_xzzz = pbuffer.data(idx_gg + 174);

    auto tg_yyyz_yyyy = pbuffer.data(idx_gg + 175);

    auto tg_yyyz_yyyz = pbuffer.data(idx_gg + 176);

    auto tg_yyyz_yyzz = pbuffer.data(idx_gg + 177);

    auto tg_yyyz_yzzz = pbuffer.data(idx_gg + 178);

    auto tg_yyyz_zzzz = pbuffer.data(idx_gg + 179);

    auto tg_yyzz_xxxx = pbuffer.data(idx_gg + 180);

    auto tg_yyzz_xxxy = pbuffer.data(idx_gg + 181);

    auto tg_yyzz_xxxz = pbuffer.data(idx_gg + 182);

    auto tg_yyzz_xxyy = pbuffer.data(idx_gg + 183);

    auto tg_yyzz_xxyz = pbuffer.data(idx_gg + 184);

    auto tg_yyzz_xxzz = pbuffer.data(idx_gg + 185);

    auto tg_yyzz_xyyy = pbuffer.data(idx_gg + 186);

    auto tg_yyzz_xyyz = pbuffer.data(idx_gg + 187);

    auto tg_yyzz_xyzz = pbuffer.data(idx_gg + 188);

    auto tg_yyzz_xzzz = pbuffer.data(idx_gg + 189);

    auto tg_yyzz_yyyy = pbuffer.data(idx_gg + 190);

    auto tg_yyzz_yyyz = pbuffer.data(idx_gg + 191);

    auto tg_yyzz_yyzz = pbuffer.data(idx_gg + 192);

    auto tg_yyzz_yzzz = pbuffer.data(idx_gg + 193);

    auto tg_yyzz_zzzz = pbuffer.data(idx_gg + 194);

    auto tg_yzzz_xxxx = pbuffer.data(idx_gg + 195);

    auto tg_yzzz_xxxy = pbuffer.data(idx_gg + 196);

    auto tg_yzzz_xxxz = pbuffer.data(idx_gg + 197);

    auto tg_yzzz_xxyy = pbuffer.data(idx_gg + 198);

    auto tg_yzzz_xxyz = pbuffer.data(idx_gg + 199);

    auto tg_yzzz_xxzz = pbuffer.data(idx_gg + 200);

    auto tg_yzzz_xyyy = pbuffer.data(idx_gg + 201);

    auto tg_yzzz_xyyz = pbuffer.data(idx_gg + 202);

    auto tg_yzzz_xyzz = pbuffer.data(idx_gg + 203);

    auto tg_yzzz_xzzz = pbuffer.data(idx_gg + 204);

    auto tg_yzzz_yyyy = pbuffer.data(idx_gg + 205);

    auto tg_yzzz_yyyz = pbuffer.data(idx_gg + 206);

    auto tg_yzzz_yyzz = pbuffer.data(idx_gg + 207);

    auto tg_yzzz_yzzz = pbuffer.data(idx_gg + 208);

    auto tg_yzzz_zzzz = pbuffer.data(idx_gg + 209);

    auto tg_zzzz_xxxx = pbuffer.data(idx_gg + 210);

    auto tg_zzzz_xxxy = pbuffer.data(idx_gg + 211);

    auto tg_zzzz_xxxz = pbuffer.data(idx_gg + 212);

    auto tg_zzzz_xxyy = pbuffer.data(idx_gg + 213);

    auto tg_zzzz_xxyz = pbuffer.data(idx_gg + 214);

    auto tg_zzzz_xxzz = pbuffer.data(idx_gg + 215);

    auto tg_zzzz_xyyy = pbuffer.data(idx_gg + 216);

    auto tg_zzzz_xyyz = pbuffer.data(idx_gg + 217);

    auto tg_zzzz_xyzz = pbuffer.data(idx_gg + 218);

    auto tg_zzzz_xzzz = pbuffer.data(idx_gg + 219);

    auto tg_zzzz_yyyy = pbuffer.data(idx_gg + 220);

    auto tg_zzzz_yyyz = pbuffer.data(idx_gg + 221);

    auto tg_zzzz_yyzz = pbuffer.data(idx_gg + 222);

    auto tg_zzzz_yzzz = pbuffer.data(idx_gg + 223);

    auto tg_zzzz_zzzz = pbuffer.data(idx_gg + 224);

    // Set up components of targeted buffer : HG

    auto tg_xxxxx_xxxx = pbuffer.data(idx_hg);

    auto tg_xxxxx_xxxy = pbuffer.data(idx_hg + 1);

    auto tg_xxxxx_xxxz = pbuffer.data(idx_hg + 2);

    auto tg_xxxxx_xxyy = pbuffer.data(idx_hg + 3);

    auto tg_xxxxx_xxyz = pbuffer.data(idx_hg + 4);

    auto tg_xxxxx_xxzz = pbuffer.data(idx_hg + 5);

    auto tg_xxxxx_xyyy = pbuffer.data(idx_hg + 6);

    auto tg_xxxxx_xyyz = pbuffer.data(idx_hg + 7);

    auto tg_xxxxx_xyzz = pbuffer.data(idx_hg + 8);

    auto tg_xxxxx_xzzz = pbuffer.data(idx_hg + 9);

    auto tg_xxxxx_yyyy = pbuffer.data(idx_hg + 10);

    auto tg_xxxxx_yyyz = pbuffer.data(idx_hg + 11);

    auto tg_xxxxx_yyzz = pbuffer.data(idx_hg + 12);

    auto tg_xxxxx_yzzz = pbuffer.data(idx_hg + 13);

    auto tg_xxxxx_zzzz = pbuffer.data(idx_hg + 14);

    auto tg_xxxxy_xxxx = pbuffer.data(idx_hg + 15);

    auto tg_xxxxy_xxxy = pbuffer.data(idx_hg + 16);

    auto tg_xxxxy_xxxz = pbuffer.data(idx_hg + 17);

    auto tg_xxxxy_xxyy = pbuffer.data(idx_hg + 18);

    auto tg_xxxxy_xxyz = pbuffer.data(idx_hg + 19);

    auto tg_xxxxy_xxzz = pbuffer.data(idx_hg + 20);

    auto tg_xxxxy_xyyy = pbuffer.data(idx_hg + 21);

    auto tg_xxxxy_xyyz = pbuffer.data(idx_hg + 22);

    auto tg_xxxxy_xyzz = pbuffer.data(idx_hg + 23);

    auto tg_xxxxy_xzzz = pbuffer.data(idx_hg + 24);

    auto tg_xxxxy_yyyy = pbuffer.data(idx_hg + 25);

    auto tg_xxxxy_yyyz = pbuffer.data(idx_hg + 26);

    auto tg_xxxxy_yyzz = pbuffer.data(idx_hg + 27);

    auto tg_xxxxy_yzzz = pbuffer.data(idx_hg + 28);

    auto tg_xxxxy_zzzz = pbuffer.data(idx_hg + 29);

    auto tg_xxxxz_xxxx = pbuffer.data(idx_hg + 30);

    auto tg_xxxxz_xxxy = pbuffer.data(idx_hg + 31);

    auto tg_xxxxz_xxxz = pbuffer.data(idx_hg + 32);

    auto tg_xxxxz_xxyy = pbuffer.data(idx_hg + 33);

    auto tg_xxxxz_xxyz = pbuffer.data(idx_hg + 34);

    auto tg_xxxxz_xxzz = pbuffer.data(idx_hg + 35);

    auto tg_xxxxz_xyyy = pbuffer.data(idx_hg + 36);

    auto tg_xxxxz_xyyz = pbuffer.data(idx_hg + 37);

    auto tg_xxxxz_xyzz = pbuffer.data(idx_hg + 38);

    auto tg_xxxxz_xzzz = pbuffer.data(idx_hg + 39);

    auto tg_xxxxz_yyyy = pbuffer.data(idx_hg + 40);

    auto tg_xxxxz_yyyz = pbuffer.data(idx_hg + 41);

    auto tg_xxxxz_yyzz = pbuffer.data(idx_hg + 42);

    auto tg_xxxxz_yzzz = pbuffer.data(idx_hg + 43);

    auto tg_xxxxz_zzzz = pbuffer.data(idx_hg + 44);

    auto tg_xxxyy_xxxx = pbuffer.data(idx_hg + 45);

    auto tg_xxxyy_xxxy = pbuffer.data(idx_hg + 46);

    auto tg_xxxyy_xxxz = pbuffer.data(idx_hg + 47);

    auto tg_xxxyy_xxyy = pbuffer.data(idx_hg + 48);

    auto tg_xxxyy_xxyz = pbuffer.data(idx_hg + 49);

    auto tg_xxxyy_xxzz = pbuffer.data(idx_hg + 50);

    auto tg_xxxyy_xyyy = pbuffer.data(idx_hg + 51);

    auto tg_xxxyy_xyyz = pbuffer.data(idx_hg + 52);

    auto tg_xxxyy_xyzz = pbuffer.data(idx_hg + 53);

    auto tg_xxxyy_xzzz = pbuffer.data(idx_hg + 54);

    auto tg_xxxyy_yyyy = pbuffer.data(idx_hg + 55);

    auto tg_xxxyy_yyyz = pbuffer.data(idx_hg + 56);

    auto tg_xxxyy_yyzz = pbuffer.data(idx_hg + 57);

    auto tg_xxxyy_yzzz = pbuffer.data(idx_hg + 58);

    auto tg_xxxyy_zzzz = pbuffer.data(idx_hg + 59);

    auto tg_xxxyz_xxxx = pbuffer.data(idx_hg + 60);

    auto tg_xxxyz_xxxy = pbuffer.data(idx_hg + 61);

    auto tg_xxxyz_xxxz = pbuffer.data(idx_hg + 62);

    auto tg_xxxyz_xxyy = pbuffer.data(idx_hg + 63);

    auto tg_xxxyz_xxyz = pbuffer.data(idx_hg + 64);

    auto tg_xxxyz_xxzz = pbuffer.data(idx_hg + 65);

    auto tg_xxxyz_xyyy = pbuffer.data(idx_hg + 66);

    auto tg_xxxyz_xyyz = pbuffer.data(idx_hg + 67);

    auto tg_xxxyz_xyzz = pbuffer.data(idx_hg + 68);

    auto tg_xxxyz_xzzz = pbuffer.data(idx_hg + 69);

    auto tg_xxxyz_yyyy = pbuffer.data(idx_hg + 70);

    auto tg_xxxyz_yyyz = pbuffer.data(idx_hg + 71);

    auto tg_xxxyz_yyzz = pbuffer.data(idx_hg + 72);

    auto tg_xxxyz_yzzz = pbuffer.data(idx_hg + 73);

    auto tg_xxxyz_zzzz = pbuffer.data(idx_hg + 74);

    auto tg_xxxzz_xxxx = pbuffer.data(idx_hg + 75);

    auto tg_xxxzz_xxxy = pbuffer.data(idx_hg + 76);

    auto tg_xxxzz_xxxz = pbuffer.data(idx_hg + 77);

    auto tg_xxxzz_xxyy = pbuffer.data(idx_hg + 78);

    auto tg_xxxzz_xxyz = pbuffer.data(idx_hg + 79);

    auto tg_xxxzz_xxzz = pbuffer.data(idx_hg + 80);

    auto tg_xxxzz_xyyy = pbuffer.data(idx_hg + 81);

    auto tg_xxxzz_xyyz = pbuffer.data(idx_hg + 82);

    auto tg_xxxzz_xyzz = pbuffer.data(idx_hg + 83);

    auto tg_xxxzz_xzzz = pbuffer.data(idx_hg + 84);

    auto tg_xxxzz_yyyy = pbuffer.data(idx_hg + 85);

    auto tg_xxxzz_yyyz = pbuffer.data(idx_hg + 86);

    auto tg_xxxzz_yyzz = pbuffer.data(idx_hg + 87);

    auto tg_xxxzz_yzzz = pbuffer.data(idx_hg + 88);

    auto tg_xxxzz_zzzz = pbuffer.data(idx_hg + 89);

    auto tg_xxyyy_xxxx = pbuffer.data(idx_hg + 90);

    auto tg_xxyyy_xxxy = pbuffer.data(idx_hg + 91);

    auto tg_xxyyy_xxxz = pbuffer.data(idx_hg + 92);

    auto tg_xxyyy_xxyy = pbuffer.data(idx_hg + 93);

    auto tg_xxyyy_xxyz = pbuffer.data(idx_hg + 94);

    auto tg_xxyyy_xxzz = pbuffer.data(idx_hg + 95);

    auto tg_xxyyy_xyyy = pbuffer.data(idx_hg + 96);

    auto tg_xxyyy_xyyz = pbuffer.data(idx_hg + 97);

    auto tg_xxyyy_xyzz = pbuffer.data(idx_hg + 98);

    auto tg_xxyyy_xzzz = pbuffer.data(idx_hg + 99);

    auto tg_xxyyy_yyyy = pbuffer.data(idx_hg + 100);

    auto tg_xxyyy_yyyz = pbuffer.data(idx_hg + 101);

    auto tg_xxyyy_yyzz = pbuffer.data(idx_hg + 102);

    auto tg_xxyyy_yzzz = pbuffer.data(idx_hg + 103);

    auto tg_xxyyy_zzzz = pbuffer.data(idx_hg + 104);

    auto tg_xxyyz_xxxx = pbuffer.data(idx_hg + 105);

    auto tg_xxyyz_xxxy = pbuffer.data(idx_hg + 106);

    auto tg_xxyyz_xxxz = pbuffer.data(idx_hg + 107);

    auto tg_xxyyz_xxyy = pbuffer.data(idx_hg + 108);

    auto tg_xxyyz_xxyz = pbuffer.data(idx_hg + 109);

    auto tg_xxyyz_xxzz = pbuffer.data(idx_hg + 110);

    auto tg_xxyyz_xyyy = pbuffer.data(idx_hg + 111);

    auto tg_xxyyz_xyyz = pbuffer.data(idx_hg + 112);

    auto tg_xxyyz_xyzz = pbuffer.data(idx_hg + 113);

    auto tg_xxyyz_xzzz = pbuffer.data(idx_hg + 114);

    auto tg_xxyyz_yyyy = pbuffer.data(idx_hg + 115);

    auto tg_xxyyz_yyyz = pbuffer.data(idx_hg + 116);

    auto tg_xxyyz_yyzz = pbuffer.data(idx_hg + 117);

    auto tg_xxyyz_yzzz = pbuffer.data(idx_hg + 118);

    auto tg_xxyyz_zzzz = pbuffer.data(idx_hg + 119);

    auto tg_xxyzz_xxxx = pbuffer.data(idx_hg + 120);

    auto tg_xxyzz_xxxy = pbuffer.data(idx_hg + 121);

    auto tg_xxyzz_xxxz = pbuffer.data(idx_hg + 122);

    auto tg_xxyzz_xxyy = pbuffer.data(idx_hg + 123);

    auto tg_xxyzz_xxyz = pbuffer.data(idx_hg + 124);

    auto tg_xxyzz_xxzz = pbuffer.data(idx_hg + 125);

    auto tg_xxyzz_xyyy = pbuffer.data(idx_hg + 126);

    auto tg_xxyzz_xyyz = pbuffer.data(idx_hg + 127);

    auto tg_xxyzz_xyzz = pbuffer.data(idx_hg + 128);

    auto tg_xxyzz_xzzz = pbuffer.data(idx_hg + 129);

    auto tg_xxyzz_yyyy = pbuffer.data(idx_hg + 130);

    auto tg_xxyzz_yyyz = pbuffer.data(idx_hg + 131);

    auto tg_xxyzz_yyzz = pbuffer.data(idx_hg + 132);

    auto tg_xxyzz_yzzz = pbuffer.data(idx_hg + 133);

    auto tg_xxyzz_zzzz = pbuffer.data(idx_hg + 134);

    auto tg_xxzzz_xxxx = pbuffer.data(idx_hg + 135);

    auto tg_xxzzz_xxxy = pbuffer.data(idx_hg + 136);

    auto tg_xxzzz_xxxz = pbuffer.data(idx_hg + 137);

    auto tg_xxzzz_xxyy = pbuffer.data(idx_hg + 138);

    auto tg_xxzzz_xxyz = pbuffer.data(idx_hg + 139);

    auto tg_xxzzz_xxzz = pbuffer.data(idx_hg + 140);

    auto tg_xxzzz_xyyy = pbuffer.data(idx_hg + 141);

    auto tg_xxzzz_xyyz = pbuffer.data(idx_hg + 142);

    auto tg_xxzzz_xyzz = pbuffer.data(idx_hg + 143);

    auto tg_xxzzz_xzzz = pbuffer.data(idx_hg + 144);

    auto tg_xxzzz_yyyy = pbuffer.data(idx_hg + 145);

    auto tg_xxzzz_yyyz = pbuffer.data(idx_hg + 146);

    auto tg_xxzzz_yyzz = pbuffer.data(idx_hg + 147);

    auto tg_xxzzz_yzzz = pbuffer.data(idx_hg + 148);

    auto tg_xxzzz_zzzz = pbuffer.data(idx_hg + 149);

    auto tg_xyyyy_xxxx = pbuffer.data(idx_hg + 150);

    auto tg_xyyyy_xxxy = pbuffer.data(idx_hg + 151);

    auto tg_xyyyy_xxxz = pbuffer.data(idx_hg + 152);

    auto tg_xyyyy_xxyy = pbuffer.data(idx_hg + 153);

    auto tg_xyyyy_xxyz = pbuffer.data(idx_hg + 154);

    auto tg_xyyyy_xxzz = pbuffer.data(idx_hg + 155);

    auto tg_xyyyy_xyyy = pbuffer.data(idx_hg + 156);

    auto tg_xyyyy_xyyz = pbuffer.data(idx_hg + 157);

    auto tg_xyyyy_xyzz = pbuffer.data(idx_hg + 158);

    auto tg_xyyyy_xzzz = pbuffer.data(idx_hg + 159);

    auto tg_xyyyy_yyyy = pbuffer.data(idx_hg + 160);

    auto tg_xyyyy_yyyz = pbuffer.data(idx_hg + 161);

    auto tg_xyyyy_yyzz = pbuffer.data(idx_hg + 162);

    auto tg_xyyyy_yzzz = pbuffer.data(idx_hg + 163);

    auto tg_xyyyy_zzzz = pbuffer.data(idx_hg + 164);

    auto tg_xyyyz_xxxx = pbuffer.data(idx_hg + 165);

    auto tg_xyyyz_xxxy = pbuffer.data(idx_hg + 166);

    auto tg_xyyyz_xxxz = pbuffer.data(idx_hg + 167);

    auto tg_xyyyz_xxyy = pbuffer.data(idx_hg + 168);

    auto tg_xyyyz_xxyz = pbuffer.data(idx_hg + 169);

    auto tg_xyyyz_xxzz = pbuffer.data(idx_hg + 170);

    auto tg_xyyyz_xyyy = pbuffer.data(idx_hg + 171);

    auto tg_xyyyz_xyyz = pbuffer.data(idx_hg + 172);

    auto tg_xyyyz_xyzz = pbuffer.data(idx_hg + 173);

    auto tg_xyyyz_xzzz = pbuffer.data(idx_hg + 174);

    auto tg_xyyyz_yyyy = pbuffer.data(idx_hg + 175);

    auto tg_xyyyz_yyyz = pbuffer.data(idx_hg + 176);

    auto tg_xyyyz_yyzz = pbuffer.data(idx_hg + 177);

    auto tg_xyyyz_yzzz = pbuffer.data(idx_hg + 178);

    auto tg_xyyyz_zzzz = pbuffer.data(idx_hg + 179);

    auto tg_xyyzz_xxxx = pbuffer.data(idx_hg + 180);

    auto tg_xyyzz_xxxy = pbuffer.data(idx_hg + 181);

    auto tg_xyyzz_xxxz = pbuffer.data(idx_hg + 182);

    auto tg_xyyzz_xxyy = pbuffer.data(idx_hg + 183);

    auto tg_xyyzz_xxyz = pbuffer.data(idx_hg + 184);

    auto tg_xyyzz_xxzz = pbuffer.data(idx_hg + 185);

    auto tg_xyyzz_xyyy = pbuffer.data(idx_hg + 186);

    auto tg_xyyzz_xyyz = pbuffer.data(idx_hg + 187);

    auto tg_xyyzz_xyzz = pbuffer.data(idx_hg + 188);

    auto tg_xyyzz_xzzz = pbuffer.data(idx_hg + 189);

    auto tg_xyyzz_yyyy = pbuffer.data(idx_hg + 190);

    auto tg_xyyzz_yyyz = pbuffer.data(idx_hg + 191);

    auto tg_xyyzz_yyzz = pbuffer.data(idx_hg + 192);

    auto tg_xyyzz_yzzz = pbuffer.data(idx_hg + 193);

    auto tg_xyyzz_zzzz = pbuffer.data(idx_hg + 194);

    auto tg_xyzzz_xxxx = pbuffer.data(idx_hg + 195);

    auto tg_xyzzz_xxxy = pbuffer.data(idx_hg + 196);

    auto tg_xyzzz_xxxz = pbuffer.data(idx_hg + 197);

    auto tg_xyzzz_xxyy = pbuffer.data(idx_hg + 198);

    auto tg_xyzzz_xxyz = pbuffer.data(idx_hg + 199);

    auto tg_xyzzz_xxzz = pbuffer.data(idx_hg + 200);

    auto tg_xyzzz_xyyy = pbuffer.data(idx_hg + 201);

    auto tg_xyzzz_xyyz = pbuffer.data(idx_hg + 202);

    auto tg_xyzzz_xyzz = pbuffer.data(idx_hg + 203);

    auto tg_xyzzz_xzzz = pbuffer.data(idx_hg + 204);

    auto tg_xyzzz_yyyy = pbuffer.data(idx_hg + 205);

    auto tg_xyzzz_yyyz = pbuffer.data(idx_hg + 206);

    auto tg_xyzzz_yyzz = pbuffer.data(idx_hg + 207);

    auto tg_xyzzz_yzzz = pbuffer.data(idx_hg + 208);

    auto tg_xyzzz_zzzz = pbuffer.data(idx_hg + 209);

    auto tg_xzzzz_xxxx = pbuffer.data(idx_hg + 210);

    auto tg_xzzzz_xxxy = pbuffer.data(idx_hg + 211);

    auto tg_xzzzz_xxxz = pbuffer.data(idx_hg + 212);

    auto tg_xzzzz_xxyy = pbuffer.data(idx_hg + 213);

    auto tg_xzzzz_xxyz = pbuffer.data(idx_hg + 214);

    auto tg_xzzzz_xxzz = pbuffer.data(idx_hg + 215);

    auto tg_xzzzz_xyyy = pbuffer.data(idx_hg + 216);

    auto tg_xzzzz_xyyz = pbuffer.data(idx_hg + 217);

    auto tg_xzzzz_xyzz = pbuffer.data(idx_hg + 218);

    auto tg_xzzzz_xzzz = pbuffer.data(idx_hg + 219);

    auto tg_xzzzz_yyyy = pbuffer.data(idx_hg + 220);

    auto tg_xzzzz_yyyz = pbuffer.data(idx_hg + 221);

    auto tg_xzzzz_yyzz = pbuffer.data(idx_hg + 222);

    auto tg_xzzzz_yzzz = pbuffer.data(idx_hg + 223);

    auto tg_xzzzz_zzzz = pbuffer.data(idx_hg + 224);

    auto tg_yyyyy_xxxx = pbuffer.data(idx_hg + 225);

    auto tg_yyyyy_xxxy = pbuffer.data(idx_hg + 226);

    auto tg_yyyyy_xxxz = pbuffer.data(idx_hg + 227);

    auto tg_yyyyy_xxyy = pbuffer.data(idx_hg + 228);

    auto tg_yyyyy_xxyz = pbuffer.data(idx_hg + 229);

    auto tg_yyyyy_xxzz = pbuffer.data(idx_hg + 230);

    auto tg_yyyyy_xyyy = pbuffer.data(idx_hg + 231);

    auto tg_yyyyy_xyyz = pbuffer.data(idx_hg + 232);

    auto tg_yyyyy_xyzz = pbuffer.data(idx_hg + 233);

    auto tg_yyyyy_xzzz = pbuffer.data(idx_hg + 234);

    auto tg_yyyyy_yyyy = pbuffer.data(idx_hg + 235);

    auto tg_yyyyy_yyyz = pbuffer.data(idx_hg + 236);

    auto tg_yyyyy_yyzz = pbuffer.data(idx_hg + 237);

    auto tg_yyyyy_yzzz = pbuffer.data(idx_hg + 238);

    auto tg_yyyyy_zzzz = pbuffer.data(idx_hg + 239);

    auto tg_yyyyz_xxxx = pbuffer.data(idx_hg + 240);

    auto tg_yyyyz_xxxy = pbuffer.data(idx_hg + 241);

    auto tg_yyyyz_xxxz = pbuffer.data(idx_hg + 242);

    auto tg_yyyyz_xxyy = pbuffer.data(idx_hg + 243);

    auto tg_yyyyz_xxyz = pbuffer.data(idx_hg + 244);

    auto tg_yyyyz_xxzz = pbuffer.data(idx_hg + 245);

    auto tg_yyyyz_xyyy = pbuffer.data(idx_hg + 246);

    auto tg_yyyyz_xyyz = pbuffer.data(idx_hg + 247);

    auto tg_yyyyz_xyzz = pbuffer.data(idx_hg + 248);

    auto tg_yyyyz_xzzz = pbuffer.data(idx_hg + 249);

    auto tg_yyyyz_yyyy = pbuffer.data(idx_hg + 250);

    auto tg_yyyyz_yyyz = pbuffer.data(idx_hg + 251);

    auto tg_yyyyz_yyzz = pbuffer.data(idx_hg + 252);

    auto tg_yyyyz_yzzz = pbuffer.data(idx_hg + 253);

    auto tg_yyyyz_zzzz = pbuffer.data(idx_hg + 254);

    auto tg_yyyzz_xxxx = pbuffer.data(idx_hg + 255);

    auto tg_yyyzz_xxxy = pbuffer.data(idx_hg + 256);

    auto tg_yyyzz_xxxz = pbuffer.data(idx_hg + 257);

    auto tg_yyyzz_xxyy = pbuffer.data(idx_hg + 258);

    auto tg_yyyzz_xxyz = pbuffer.data(idx_hg + 259);

    auto tg_yyyzz_xxzz = pbuffer.data(idx_hg + 260);

    auto tg_yyyzz_xyyy = pbuffer.data(idx_hg + 261);

    auto tg_yyyzz_xyyz = pbuffer.data(idx_hg + 262);

    auto tg_yyyzz_xyzz = pbuffer.data(idx_hg + 263);

    auto tg_yyyzz_xzzz = pbuffer.data(idx_hg + 264);

    auto tg_yyyzz_yyyy = pbuffer.data(idx_hg + 265);

    auto tg_yyyzz_yyyz = pbuffer.data(idx_hg + 266);

    auto tg_yyyzz_yyzz = pbuffer.data(idx_hg + 267);

    auto tg_yyyzz_yzzz = pbuffer.data(idx_hg + 268);

    auto tg_yyyzz_zzzz = pbuffer.data(idx_hg + 269);

    auto tg_yyzzz_xxxx = pbuffer.data(idx_hg + 270);

    auto tg_yyzzz_xxxy = pbuffer.data(idx_hg + 271);

    auto tg_yyzzz_xxxz = pbuffer.data(idx_hg + 272);

    auto tg_yyzzz_xxyy = pbuffer.data(idx_hg + 273);

    auto tg_yyzzz_xxyz = pbuffer.data(idx_hg + 274);

    auto tg_yyzzz_xxzz = pbuffer.data(idx_hg + 275);

    auto tg_yyzzz_xyyy = pbuffer.data(idx_hg + 276);

    auto tg_yyzzz_xyyz = pbuffer.data(idx_hg + 277);

    auto tg_yyzzz_xyzz = pbuffer.data(idx_hg + 278);

    auto tg_yyzzz_xzzz = pbuffer.data(idx_hg + 279);

    auto tg_yyzzz_yyyy = pbuffer.data(idx_hg + 280);

    auto tg_yyzzz_yyyz = pbuffer.data(idx_hg + 281);

    auto tg_yyzzz_yyzz = pbuffer.data(idx_hg + 282);

    auto tg_yyzzz_yzzz = pbuffer.data(idx_hg + 283);

    auto tg_yyzzz_zzzz = pbuffer.data(idx_hg + 284);

    auto tg_yzzzz_xxxx = pbuffer.data(idx_hg + 285);

    auto tg_yzzzz_xxxy = pbuffer.data(idx_hg + 286);

    auto tg_yzzzz_xxxz = pbuffer.data(idx_hg + 287);

    auto tg_yzzzz_xxyy = pbuffer.data(idx_hg + 288);

    auto tg_yzzzz_xxyz = pbuffer.data(idx_hg + 289);

    auto tg_yzzzz_xxzz = pbuffer.data(idx_hg + 290);

    auto tg_yzzzz_xyyy = pbuffer.data(idx_hg + 291);

    auto tg_yzzzz_xyyz = pbuffer.data(idx_hg + 292);

    auto tg_yzzzz_xyzz = pbuffer.data(idx_hg + 293);

    auto tg_yzzzz_xzzz = pbuffer.data(idx_hg + 294);

    auto tg_yzzzz_yyyy = pbuffer.data(idx_hg + 295);

    auto tg_yzzzz_yyyz = pbuffer.data(idx_hg + 296);

    auto tg_yzzzz_yyzz = pbuffer.data(idx_hg + 297);

    auto tg_yzzzz_yzzz = pbuffer.data(idx_hg + 298);

    auto tg_yzzzz_zzzz = pbuffer.data(idx_hg + 299);

    auto tg_zzzzz_xxxx = pbuffer.data(idx_hg + 300);

    auto tg_zzzzz_xxxy = pbuffer.data(idx_hg + 301);

    auto tg_zzzzz_xxxz = pbuffer.data(idx_hg + 302);

    auto tg_zzzzz_xxyy = pbuffer.data(idx_hg + 303);

    auto tg_zzzzz_xxyz = pbuffer.data(idx_hg + 304);

    auto tg_zzzzz_xxzz = pbuffer.data(idx_hg + 305);

    auto tg_zzzzz_xyyy = pbuffer.data(idx_hg + 306);

    auto tg_zzzzz_xyyz = pbuffer.data(idx_hg + 307);

    auto tg_zzzzz_xyzz = pbuffer.data(idx_hg + 308);

    auto tg_zzzzz_xzzz = pbuffer.data(idx_hg + 309);

    auto tg_zzzzz_yyyy = pbuffer.data(idx_hg + 310);

    auto tg_zzzzz_yyyz = pbuffer.data(idx_hg + 311);

    auto tg_zzzzz_yyzz = pbuffer.data(idx_hg + 312);

    auto tg_zzzzz_yzzz = pbuffer.data(idx_hg + 313);

    auto tg_zzzzz_zzzz = pbuffer.data(idx_hg + 314);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxx_xxxx, tg_xxx_xxxy, tg_xxx_xxxz, tg_xxx_xxyy, tg_xxx_xxyz, tg_xxx_xxzz, tg_xxx_xyyy, tg_xxx_xyyz, tg_xxx_xyzz, tg_xxx_xzzz, tg_xxx_yyyy, tg_xxx_yyyz, tg_xxx_yyzz, tg_xxx_yzzz, tg_xxx_zzzz, tg_xxxx_xxx, tg_xxxx_xxxx, tg_xxxx_xxxy, tg_xxxx_xxxz, tg_xxxx_xxy, tg_xxxx_xxyy, tg_xxxx_xxyz, tg_xxxx_xxz, tg_xxxx_xxzz, tg_xxxx_xyy, tg_xxxx_xyyy, tg_xxxx_xyyz, tg_xxxx_xyz, tg_xxxx_xyzz, tg_xxxx_xzz, tg_xxxx_xzzz, tg_xxxx_yyy, tg_xxxx_yyyy, tg_xxxx_yyyz, tg_xxxx_yyz, tg_xxxx_yyzz, tg_xxxx_yzz, tg_xxxx_yzzz, tg_xxxx_zzz, tg_xxxx_zzzz, tg_xxxxx_xxxx, tg_xxxxx_xxxy, tg_xxxxx_xxxz, tg_xxxxx_xxyy, tg_xxxxx_xxyz, tg_xxxxx_xxzz, tg_xxxxx_xyyy, tg_xxxxx_xyyz, tg_xxxxx_xyzz, tg_xxxxx_xzzz, tg_xxxxx_yyyy, tg_xxxxx_yyyz, tg_xxxxx_yyzz, tg_xxxxx_yzzz, tg_xxxxx_zzzz, tg_xxxxy_xxxx, tg_xxxxy_xxxy, tg_xxxxy_xxxz, tg_xxxxy_xxyy, tg_xxxxy_xxyz, tg_xxxxy_xxzz, tg_xxxxy_xyyy, tg_xxxxy_xyyz, tg_xxxxy_xyzz, tg_xxxxy_xzzz, tg_xxxxy_yyyy, tg_xxxxy_yyyz, tg_xxxxy_yyzz, tg_xxxxy_yzzz, tg_xxxxy_zzzz, tg_xxxxz_xxxx, tg_xxxxz_xxxy, tg_xxxxz_xxxz, tg_xxxxz_xxyy, tg_xxxxz_xxyz, tg_xxxxz_xxzz, tg_xxxxz_xyyy, tg_xxxxz_xyyz, tg_xxxxz_xyzz, tg_xxxxz_xzzz, tg_xxxxz_yyyy, tg_xxxxz_yyyz, tg_xxxxz_yyzz, tg_xxxxz_yzzz, tg_xxxxz_zzzz, tg_xxxy_xxxx, tg_xxxy_xxxy, tg_xxxy_xxxz, tg_xxxy_xxyy, tg_xxxy_xxzz, tg_xxxy_xyyy, tg_xxxy_xzzz, tg_xxxy_yyyy, tg_xxxy_yyyz, tg_xxxy_yyzz, tg_xxxy_yzzz, tg_xxxyy_xxxx, tg_xxxyy_xxxy, tg_xxxyy_xxxz, tg_xxxyy_xxyy, tg_xxxyy_xxyz, tg_xxxyy_xxzz, tg_xxxyy_xyyy, tg_xxxyy_xyyz, tg_xxxyy_xyzz, tg_xxxyy_xzzz, tg_xxxyy_yyyy, tg_xxxyy_yyyz, tg_xxxyy_yyzz, tg_xxxyy_yzzz, tg_xxxyy_zzzz, tg_xxxyz_xxxx, tg_xxxyz_xxxy, tg_xxxyz_xxxz, tg_xxxyz_xxyy, tg_xxxyz_xxyz, tg_xxxyz_xxzz, tg_xxxyz_xyyy, tg_xxxyz_xyyz, tg_xxxyz_xyzz, tg_xxxyz_xzzz, tg_xxxyz_yyyy, tg_xxxyz_yyyz, tg_xxxyz_yyzz, tg_xxxyz_yzzz, tg_xxxyz_zzzz, tg_xxxz_xxxx, tg_xxxz_xxxy, tg_xxxz_xxxz, tg_xxxz_xxyy, tg_xxxz_xxyz, tg_xxxz_xxz, tg_xxxz_xxzz, tg_xxxz_xyyy, tg_xxxz_xyyz, tg_xxxz_xyz, tg_xxxz_xyzz, tg_xxxz_xzz, tg_xxxz_xzzz, tg_xxxz_yyyz, tg_xxxz_yyzz, tg_xxxz_yzzz, tg_xxxz_zzzz, tg_xxxzz_xxxx, tg_xxxzz_xxxy, tg_xxxzz_xxxz, tg_xxxzz_xxyy, tg_xxxzz_xxyz, tg_xxxzz_xxzz, tg_xxxzz_xyyy, tg_xxxzz_xyyz, tg_xxxzz_xyzz, tg_xxxzz_xzzz, tg_xxxzz_yyyy, tg_xxxzz_yyyz, tg_xxxzz_yyzz, tg_xxxzz_yzzz, tg_xxxzz_zzzz, tg_xxy_xxxx, tg_xxy_xxxz, tg_xxy_xxzz, tg_xxy_xzzz, tg_xxy_yyyy, tg_xxy_yyyz, tg_xxy_yyzz, tg_xxy_yzzz, tg_xxyy_xxxx, tg_xxyy_xxxy, tg_xxyy_xxxz, tg_xxyy_xxy, tg_xxyy_xxyy, tg_xxyy_xxyz, tg_xxyy_xxzz, tg_xxyy_xyy, tg_xxyy_xyyy, tg_xxyy_xyyz, tg_xxyy_xyz, tg_xxyy_xyzz, tg_xxyy_xzzz, tg_xxyy_yyy, tg_xxyy_yyyy, tg_xxyy_yyyz, tg_xxyy_yyz, tg_xxyy_yyzz, tg_xxyy_yzz, tg_xxyy_yzzz, tg_xxyy_zzzz, tg_xxyyy_xxxx, tg_xxyyy_xxxy, tg_xxyyy_xxxz, tg_xxyyy_xxyy, tg_xxyyy_xxyz, tg_xxyyy_xxzz, tg_xxyyy_xyyy, tg_xxyyy_xyyz, tg_xxyyy_xyzz, tg_xxyyy_xzzz, tg_xxyyy_yyyy, tg_xxyyy_yyyz, tg_xxyyy_yyzz, tg_xxyyy_yzzz, tg_xxyyy_zzzz, tg_xxyyz_xxxx, tg_xxyyz_xxxy, tg_xxyyz_xxxz, tg_xxyyz_xxyy, tg_xxyyz_xxyz, tg_xxyyz_xxzz, tg_xxyyz_xyyy, tg_xxyyz_xyyz, tg_xxyyz_xyzz, tg_xxyyz_xzzz, tg_xxyyz_yyyy, tg_xxyyz_yyyz, tg_xxyyz_yyzz, tg_xxyyz_yzzz, tg_xxyyz_zzzz, tg_xxyz_xxxz, tg_xxyz_xxzz, tg_xxyz_xzzz, tg_xxyz_yyyz, tg_xxyz_yyzz, tg_xxyz_yzzz, tg_xxyzz_xxxx, tg_xxyzz_xxxy, tg_xxyzz_xxxz, tg_xxyzz_xxyy, tg_xxyzz_xxyz, tg_xxyzz_xxzz, tg_xxyzz_xyyy, tg_xxyzz_xyyz, tg_xxyzz_xyzz, tg_xxyzz_xzzz, tg_xxyzz_yyyy, tg_xxyzz_yyyz, tg_xxyzz_yyzz, tg_xxyzz_yzzz, tg_xxyzz_zzzz, tg_xxz_xxxx, tg_xxz_xxxy, tg_xxz_xxxz, tg_xxz_xxyy, tg_xxz_xxzz, tg_xxz_xyyy, tg_xxz_xzzz, tg_xxz_yyyz, tg_xxz_yyzz, tg_xxz_yzzz, tg_xxz_zzzz, tg_xxzz_xxx, tg_xxzz_xxxx, tg_xxzz_xxxy, tg_xxzz_xxxz, tg_xxzz_xxy, tg_xxzz_xxyy, tg_xxzz_xxyz, tg_xxzz_xxz, tg_xxzz_xxzz, tg_xxzz_xyy, tg_xxzz_xyyy, tg_xxzz_xyyz, tg_xxzz_xyz, tg_xxzz_xyzz, tg_xxzz_xzz, tg_xxzz_xzzz, tg_xxzz_yyyy, tg_xxzz_yyyz, tg_xxzz_yyz, tg_xxzz_yyzz, tg_xxzz_yzz, tg_xxzz_yzzz, tg_xxzz_zzz, tg_xxzz_zzzz, tg_xxzzz_xxxx, tg_xxzzz_xxxy, tg_xxzzz_xxxz, tg_xxzzz_xxyy, tg_xxzzz_xxyz, tg_xxzzz_xxzz, tg_xxzzz_xyyy, tg_xxzzz_xyyz, tg_xxzzz_xyzz, tg_xxzzz_xzzz, tg_xxzzz_yyyy, tg_xxzzz_yyyz, tg_xxzzz_yyzz, tg_xxzzz_yzzz, tg_xxzzz_zzzz, tg_xyy_xxxy, tg_xyy_xxyy, tg_xyy_xxyz, tg_xyy_xyyy, tg_xyy_xyyz, tg_xyy_xyzz, tg_xyy_yyyy, tg_xyy_yyyz, tg_xyy_yyzz, tg_xyy_yzzz, tg_xyy_zzzz, tg_xyyy_xxxx, tg_xyyy_xxxy, tg_xyyy_xxy, tg_xyyy_xxyy, tg_xyyy_xxyz, tg_xyyy_xyy, tg_xyyy_xyyy, tg_xyyy_xyyz, tg_xyyy_xyz, tg_xyyy_xyzz, tg_xyyy_yyy, tg_xyyy_yyyy, tg_xyyy_yyyz, tg_xyyy_yyz, tg_xyyy_yyzz, tg_xyyy_yzz, tg_xyyy_yzzz, tg_xyyy_zzzz, tg_xyyyy_xxxx, tg_xyyyy_xxxy, tg_xyyyy_xxxz, tg_xyyyy_xxyy, tg_xyyyy_xxyz, tg_xyyyy_xxzz, tg_xyyyy_xyyy, tg_xyyyy_xyyz, tg_xyyyy_xyzz, tg_xyyyy_xzzz, tg_xyyyy_yyyy, tg_xyyyy_yyyz, tg_xyyyy_yyzz, tg_xyyyy_yzzz, tg_xyyyy_zzzz, tg_xyyyz_xxxx, tg_xyyyz_xxxy, tg_xyyyz_xxxz, tg_xyyyz_xxyy, tg_xyyyz_xxyz, tg_xyyyz_xxzz, tg_xyyyz_xyyy, tg_xyyyz_xyyz, tg_xyyyz_xyzz, tg_xyyyz_xzzz, tg_xyyyz_yyyy, tg_xyyyz_yyyz, tg_xyyyz_yyzz, tg_xyyyz_yzzz, tg_xyyyz_zzzz, tg_xyyz_yyyz, tg_xyyz_yyzz, tg_xyyz_yzzz, tg_xyyz_zzzz, tg_xyyzz_xxxx, tg_xyyzz_xxxy, tg_xyyzz_xxxz, tg_xyyzz_xxyy, tg_xyyzz_xxyz, tg_xyyzz_xxzz, tg_xyyzz_xyyy, tg_xyyzz_xyyz, tg_xyyzz_xyzz, tg_xyyzz_xzzz, tg_xyyzz_yyyy, tg_xyyzz_yyyz, tg_xyyzz_yyzz, tg_xyyzz_yzzz, tg_xyyzz_zzzz, tg_xyz_yyyz, tg_xyz_yyzz, tg_xyz_yzzz, tg_xyzz_yyyy, tg_xyzz_yyyz, tg_xyzz_yyzz, tg_xyzz_yzzz, tg_xyzzz_xxxx, tg_xyzzz_xxxy, tg_xyzzz_xxxz, tg_xyzzz_xxyy, tg_xyzzz_xxyz, tg_xyzzz_xxzz, tg_xyzzz_xyyy, tg_xyzzz_xyyz, tg_xyzzz_xyzz, tg_xyzzz_xzzz, tg_xyzzz_yyyy, tg_xyzzz_yyyz, tg_xyzzz_yyzz, tg_xyzzz_yzzz, tg_xyzzz_zzzz, tg_xzz_xxxz, tg_xzz_xxyz, tg_xzz_xxzz, tg_xzz_xyyz, tg_xzz_xyzz, tg_xzz_xzzz, tg_xzz_yyyy, tg_xzz_yyyz, tg_xzz_yyzz, tg_xzz_yzzz, tg_xzz_zzzz, tg_xzzz_xxxx, tg_xzzz_xxxz, tg_xzzz_xxyz, tg_xzzz_xxz, tg_xzzz_xxzz, tg_xzzz_xyyz, tg_xzzz_xyz, tg_xzzz_xyzz, tg_xzzz_xzz, tg_xzzz_xzzz, tg_xzzz_yyyy, tg_xzzz_yyyz, tg_xzzz_yyz, tg_xzzz_yyzz, tg_xzzz_yzz, tg_xzzz_yzzz, tg_xzzz_zzz, tg_xzzz_zzzz, tg_xzzzz_xxxx, tg_xzzzz_xxxy, tg_xzzzz_xxxz, tg_xzzzz_xxyy, tg_xzzzz_xxyz, tg_xzzzz_xxzz, tg_xzzzz_xyyy, tg_xzzzz_xyyz, tg_xzzzz_xyzz, tg_xzzzz_xzzz, tg_xzzzz_yyyy, tg_xzzzz_yyyz, tg_xzzzz_yyzz, tg_xzzzz_yzzz, tg_xzzzz_zzzz, tg_yyy_xxxx, tg_yyy_xxxy, tg_yyy_xxxz, tg_yyy_xxyy, tg_yyy_xxyz, tg_yyy_xxzz, tg_yyy_xyyy, tg_yyy_xyyz, tg_yyy_xyzz, tg_yyy_xzzz, tg_yyy_yyyy, tg_yyy_yyyz, tg_yyy_yyzz, tg_yyy_yzzz, tg_yyy_zzzz, tg_yyyy_xxx, tg_yyyy_xxxx, tg_yyyy_xxxy, tg_yyyy_xxxz, tg_yyyy_xxy, tg_yyyy_xxyy, tg_yyyy_xxyz, tg_yyyy_xxz, tg_yyyy_xxzz, tg_yyyy_xyy, tg_yyyy_xyyy, tg_yyyy_xyyz, tg_yyyy_xyz, tg_yyyy_xyzz, tg_yyyy_xzz, tg_yyyy_xzzz, tg_yyyy_yyy, tg_yyyy_yyyy, tg_yyyy_yyyz, tg_yyyy_yyz, tg_yyyy_yyzz, tg_yyyy_yzz, tg_yyyy_yzzz, tg_yyyy_zzz, tg_yyyy_zzzz, tg_yyyyy_xxxx, tg_yyyyy_xxxy, tg_yyyyy_xxxz, tg_yyyyy_xxyy, tg_yyyyy_xxyz, tg_yyyyy_xxzz, tg_yyyyy_xyyy, tg_yyyyy_xyyz, tg_yyyyy_xyzz, tg_yyyyy_xzzz, tg_yyyyy_yyyy, tg_yyyyy_yyyz, tg_yyyyy_yyzz, tg_yyyyy_yzzz, tg_yyyyy_zzzz, tg_yyyyz_xxxx, tg_yyyyz_xxxy, tg_yyyyz_xxxz, tg_yyyyz_xxyy, tg_yyyyz_xxyz, tg_yyyyz_xxzz, tg_yyyyz_xyyy, tg_yyyyz_xyyz, tg_yyyyz_xyzz, tg_yyyyz_xzzz, tg_yyyyz_yyyy, tg_yyyyz_yyyz, tg_yyyyz_yyzz, tg_yyyyz_yzzz, tg_yyyyz_zzzz, tg_yyyz_xxxy, tg_yyyz_xxxz, tg_yyyz_xxyy, tg_yyyz_xxyz, tg_yyyz_xxz, tg_yyyz_xxzz, tg_yyyz_xyyy, tg_yyyz_xyyz, tg_yyyz_xyz, tg_yyyz_xyzz, tg_yyyz_xzz, tg_yyyz_xzzz, tg_yyyz_yyyy, tg_yyyz_yyyz, tg_yyyz_yyz, tg_yyyz_yyzz, tg_yyyz_yzz, tg_yyyz_yzzz, tg_yyyz_zzz, tg_yyyz_zzzz, tg_yyyzz_xxxx, tg_yyyzz_xxxy, tg_yyyzz_xxxz, tg_yyyzz_xxyy, tg_yyyzz_xxyz, tg_yyyzz_xxzz, tg_yyyzz_xyyy, tg_yyyzz_xyyz, tg_yyyzz_xyzz, tg_yyyzz_xzzz, tg_yyyzz_yyyy, tg_yyyzz_yyyz, tg_yyyzz_yyzz, tg_yyyzz_yzzz, tg_yyyzz_zzzz, tg_yyz_xxxy, tg_yyz_xxxz, tg_yyz_xxyy, tg_yyz_xxzz, tg_yyz_xyyy, tg_yyz_xzzz, tg_yyz_yyyy, tg_yyz_yyyz, tg_yyz_yyzz, tg_yyz_yzzz, tg_yyz_zzzz, tg_yyzz_xxx, tg_yyzz_xxxx, tg_yyzz_xxxy, tg_yyzz_xxxz, tg_yyzz_xxy, tg_yyzz_xxyy, tg_yyzz_xxyz, tg_yyzz_xxz, tg_yyzz_xxzz, tg_yyzz_xyy, tg_yyzz_xyyy, tg_yyzz_xyyz, tg_yyzz_xyz, tg_yyzz_xyzz, tg_yyzz_xzz, tg_yyzz_xzzz, tg_yyzz_yyy, tg_yyzz_yyyy, tg_yyzz_yyyz, tg_yyzz_yyz, tg_yyzz_yyzz, tg_yyzz_yzz, tg_yyzz_yzzz, tg_yyzz_zzz, tg_yyzz_zzzz, tg_yyzzz_xxxx, tg_yyzzz_xxxy, tg_yyzzz_xxxz, tg_yyzzz_xxyy, tg_yyzzz_xxyz, tg_yyzzz_xxzz, tg_yyzzz_xyyy, tg_yyzzz_xyyz, tg_yyzzz_xyzz, tg_yyzzz_xzzz, tg_yyzzz_yyyy, tg_yyzzz_yyyz, tg_yyzzz_yyzz, tg_yyzzz_yzzz, tg_yyzzz_zzzz, tg_yzz_xxxx, tg_yzz_xxxz, tg_yzz_xxyz, tg_yzz_xxzz, tg_yzz_xyyz, tg_yzz_xyzz, tg_yzz_xzzz, tg_yzz_yyyy, tg_yzz_yyyz, tg_yzz_yyzz, tg_yzz_yzzz, tg_yzz_zzzz, tg_yzzz_xxxx, tg_yzzz_xxxy, tg_yzzz_xxxz, tg_yzzz_xxy, tg_yzzz_xxyy, tg_yzzz_xxyz, tg_yzzz_xxz, tg_yzzz_xxzz, tg_yzzz_xyy, tg_yzzz_xyyy, tg_yzzz_xyyz, tg_yzzz_xyz, tg_yzzz_xyzz, tg_yzzz_xzz, tg_yzzz_xzzz, tg_yzzz_yyy, tg_yzzz_yyyy, tg_yzzz_yyyz, tg_yzzz_yyz, tg_yzzz_yyzz, tg_yzzz_yzz, tg_yzzz_yzzz, tg_yzzz_zzz, tg_yzzz_zzzz, tg_yzzzz_xxxx, tg_yzzzz_xxxy, tg_yzzzz_xxxz, tg_yzzzz_xxyy, tg_yzzzz_xxyz, tg_yzzzz_xxzz, tg_yzzzz_xyyy, tg_yzzzz_xyyz, tg_yzzzz_xyzz, tg_yzzzz_xzzz, tg_yzzzz_yyyy, tg_yzzzz_yyyz, tg_yzzzz_yyzz, tg_yzzzz_yzzz, tg_yzzzz_zzzz, tg_zzz_xxxx, tg_zzz_xxxy, tg_zzz_xxxz, tg_zzz_xxyy, tg_zzz_xxyz, tg_zzz_xxzz, tg_zzz_xyyy, tg_zzz_xyyz, tg_zzz_xyzz, tg_zzz_xzzz, tg_zzz_yyyy, tg_zzz_yyyz, tg_zzz_yyzz, tg_zzz_yzzz, tg_zzz_zzzz, tg_zzzz_xxx, tg_zzzz_xxxx, tg_zzzz_xxxy, tg_zzzz_xxxz, tg_zzzz_xxy, tg_zzzz_xxyy, tg_zzzz_xxyz, tg_zzzz_xxz, tg_zzzz_xxzz, tg_zzzz_xyy, tg_zzzz_xyyy, tg_zzzz_xyyz, tg_zzzz_xyz, tg_zzzz_xyzz, tg_zzzz_xzz, tg_zzzz_xzzz, tg_zzzz_yyy, tg_zzzz_yyyy, tg_zzzz_yyyz, tg_zzzz_yyz, tg_zzzz_yyzz, tg_zzzz_yzz, tg_zzzz_yzzz, tg_zzzz_zzz, tg_zzzz_zzzz, tg_zzzzz_xxxx, tg_zzzzz_xxxy, tg_zzzzz_xxxz, tg_zzzzz_xxyy, tg_zzzzz_xxyz, tg_zzzzz_xxzz, tg_zzzzz_xyyy, tg_zzzzz_xyyz, tg_zzzzz_xyzz, tg_zzzzz_xzzz, tg_zzzzz_yyyy, tg_zzzzz_yyyz, tg_zzzzz_yyzz, tg_zzzzz_yzzz, tg_zzzzz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxx_xxxx[i] = 4.0 * tg_xxx_xxxx[i] * fxi[i] + 4.0 * tg_xxxx_xxx[i] * fxi[i] + tg_xxxx_xxxx[i] * ra_x[i];

        tg_xxxxx_xxxy[i] = 4.0 * tg_xxx_xxxy[i] * fxi[i] + 3.0 * tg_xxxx_xxy[i] * fxi[i] + tg_xxxx_xxxy[i] * ra_x[i];

        tg_xxxxx_xxxz[i] = 4.0 * tg_xxx_xxxz[i] * fxi[i] + 3.0 * tg_xxxx_xxz[i] * fxi[i] + tg_xxxx_xxxz[i] * ra_x[i];

        tg_xxxxx_xxyy[i] = 4.0 * tg_xxx_xxyy[i] * fxi[i] + 2.0 * tg_xxxx_xyy[i] * fxi[i] + tg_xxxx_xxyy[i] * ra_x[i];

        tg_xxxxx_xxyz[i] = 4.0 * tg_xxx_xxyz[i] * fxi[i] + 2.0 * tg_xxxx_xyz[i] * fxi[i] + tg_xxxx_xxyz[i] * ra_x[i];

        tg_xxxxx_xxzz[i] = 4.0 * tg_xxx_xxzz[i] * fxi[i] + 2.0 * tg_xxxx_xzz[i] * fxi[i] + tg_xxxx_xxzz[i] * ra_x[i];

        tg_xxxxx_xyyy[i] = 4.0 * tg_xxx_xyyy[i] * fxi[i] + tg_xxxx_yyy[i] * fxi[i] + tg_xxxx_xyyy[i] * ra_x[i];

        tg_xxxxx_xyyz[i] = 4.0 * tg_xxx_xyyz[i] * fxi[i] + tg_xxxx_yyz[i] * fxi[i] + tg_xxxx_xyyz[i] * ra_x[i];

        tg_xxxxx_xyzz[i] = 4.0 * tg_xxx_xyzz[i] * fxi[i] + tg_xxxx_yzz[i] * fxi[i] + tg_xxxx_xyzz[i] * ra_x[i];

        tg_xxxxx_xzzz[i] = 4.0 * tg_xxx_xzzz[i] * fxi[i] + tg_xxxx_zzz[i] * fxi[i] + tg_xxxx_xzzz[i] * ra_x[i];

        tg_xxxxx_yyyy[i] = 4.0 * tg_xxx_yyyy[i] * fxi[i] + tg_xxxx_yyyy[i] * ra_x[i];

        tg_xxxxx_yyyz[i] = 4.0 * tg_xxx_yyyz[i] * fxi[i] + tg_xxxx_yyyz[i] * ra_x[i];

        tg_xxxxx_yyzz[i] = 4.0 * tg_xxx_yyzz[i] * fxi[i] + tg_xxxx_yyzz[i] * ra_x[i];

        tg_xxxxx_yzzz[i] = 4.0 * tg_xxx_yzzz[i] * fxi[i] + tg_xxxx_yzzz[i] * ra_x[i];

        tg_xxxxx_zzzz[i] = 4.0 * tg_xxx_zzzz[i] * fxi[i] + tg_xxxx_zzzz[i] * ra_x[i];

        tg_xxxxy_xxxx[i] = tg_xxxx_xxxx[i] * ra_y[i];

        tg_xxxxy_xxxy[i] = tg_xxxx_xxx[i] * fxi[i] + tg_xxxx_xxxy[i] * ra_y[i];

        tg_xxxxy_xxxz[i] = tg_xxxx_xxxz[i] * ra_y[i];

        tg_xxxxy_xxyy[i] = 2.0 * tg_xxxx_xxy[i] * fxi[i] + tg_xxxx_xxyy[i] * ra_y[i];

        tg_xxxxy_xxyz[i] = tg_xxxx_xxz[i] * fxi[i] + tg_xxxx_xxyz[i] * ra_y[i];

        tg_xxxxy_xxzz[i] = tg_xxxx_xxzz[i] * ra_y[i];

        tg_xxxxy_xyyy[i] = 3.0 * tg_xxxx_xyy[i] * fxi[i] + tg_xxxx_xyyy[i] * ra_y[i];

        tg_xxxxy_xyyz[i] = 2.0 * tg_xxxx_xyz[i] * fxi[i] + tg_xxxx_xyyz[i] * ra_y[i];

        tg_xxxxy_xyzz[i] = tg_xxxx_xzz[i] * fxi[i] + tg_xxxx_xyzz[i] * ra_y[i];

        tg_xxxxy_xzzz[i] = tg_xxxx_xzzz[i] * ra_y[i];

        tg_xxxxy_yyyy[i] = 3.0 * tg_xxy_yyyy[i] * fxi[i] + tg_xxxy_yyyy[i] * ra_x[i];

        tg_xxxxy_yyyz[i] = 3.0 * tg_xxy_yyyz[i] * fxi[i] + tg_xxxy_yyyz[i] * ra_x[i];

        tg_xxxxy_yyzz[i] = 3.0 * tg_xxy_yyzz[i] * fxi[i] + tg_xxxy_yyzz[i] * ra_x[i];

        tg_xxxxy_yzzz[i] = 3.0 * tg_xxy_yzzz[i] * fxi[i] + tg_xxxy_yzzz[i] * ra_x[i];

        tg_xxxxy_zzzz[i] = tg_xxxx_zzzz[i] * ra_y[i];

        tg_xxxxz_xxxx[i] = tg_xxxx_xxxx[i] * ra_z[i];

        tg_xxxxz_xxxy[i] = tg_xxxx_xxxy[i] * ra_z[i];

        tg_xxxxz_xxxz[i] = tg_xxxx_xxx[i] * fxi[i] + tg_xxxx_xxxz[i] * ra_z[i];

        tg_xxxxz_xxyy[i] = tg_xxxx_xxyy[i] * ra_z[i];

        tg_xxxxz_xxyz[i] = tg_xxxx_xxy[i] * fxi[i] + tg_xxxx_xxyz[i] * ra_z[i];

        tg_xxxxz_xxzz[i] = 2.0 * tg_xxxx_xxz[i] * fxi[i] + tg_xxxx_xxzz[i] * ra_z[i];

        tg_xxxxz_xyyy[i] = tg_xxxx_xyyy[i] * ra_z[i];

        tg_xxxxz_xyyz[i] = tg_xxxx_xyy[i] * fxi[i] + tg_xxxx_xyyz[i] * ra_z[i];

        tg_xxxxz_xyzz[i] = 2.0 * tg_xxxx_xyz[i] * fxi[i] + tg_xxxx_xyzz[i] * ra_z[i];

        tg_xxxxz_xzzz[i] = 3.0 * tg_xxxx_xzz[i] * fxi[i] + tg_xxxx_xzzz[i] * ra_z[i];

        tg_xxxxz_yyyy[i] = tg_xxxx_yyyy[i] * ra_z[i];

        tg_xxxxz_yyyz[i] = 3.0 * tg_xxz_yyyz[i] * fxi[i] + tg_xxxz_yyyz[i] * ra_x[i];

        tg_xxxxz_yyzz[i] = 3.0 * tg_xxz_yyzz[i] * fxi[i] + tg_xxxz_yyzz[i] * ra_x[i];

        tg_xxxxz_yzzz[i] = 3.0 * tg_xxz_yzzz[i] * fxi[i] + tg_xxxz_yzzz[i] * ra_x[i];

        tg_xxxxz_zzzz[i] = 3.0 * tg_xxz_zzzz[i] * fxi[i] + tg_xxxz_zzzz[i] * ra_x[i];

        tg_xxxyy_xxxx[i] = tg_xxx_xxxx[i] * fxi[i] + tg_xxxy_xxxx[i] * ra_y[i];

        tg_xxxyy_xxxy[i] = 2.0 * tg_xyy_xxxy[i] * fxi[i] + 3.0 * tg_xxyy_xxy[i] * fxi[i] + tg_xxyy_xxxy[i] * ra_x[i];

        tg_xxxyy_xxxz[i] = tg_xxx_xxxz[i] * fxi[i] + tg_xxxy_xxxz[i] * ra_y[i];

        tg_xxxyy_xxyy[i] = 2.0 * tg_xyy_xxyy[i] * fxi[i] + 2.0 * tg_xxyy_xyy[i] * fxi[i] + tg_xxyy_xxyy[i] * ra_x[i];

        tg_xxxyy_xxyz[i] = 2.0 * tg_xyy_xxyz[i] * fxi[i] + 2.0 * tg_xxyy_xyz[i] * fxi[i] + tg_xxyy_xxyz[i] * ra_x[i];

        tg_xxxyy_xxzz[i] = tg_xxx_xxzz[i] * fxi[i] + tg_xxxy_xxzz[i] * ra_y[i];

        tg_xxxyy_xyyy[i] = 2.0 * tg_xyy_xyyy[i] * fxi[i] + tg_xxyy_yyy[i] * fxi[i] + tg_xxyy_xyyy[i] * ra_x[i];

        tg_xxxyy_xyyz[i] = 2.0 * tg_xyy_xyyz[i] * fxi[i] + tg_xxyy_yyz[i] * fxi[i] + tg_xxyy_xyyz[i] * ra_x[i];

        tg_xxxyy_xyzz[i] = 2.0 * tg_xyy_xyzz[i] * fxi[i] + tg_xxyy_yzz[i] * fxi[i] + tg_xxyy_xyzz[i] * ra_x[i];

        tg_xxxyy_xzzz[i] = tg_xxx_xzzz[i] * fxi[i] + tg_xxxy_xzzz[i] * ra_y[i];

        tg_xxxyy_yyyy[i] = 2.0 * tg_xyy_yyyy[i] * fxi[i] + tg_xxyy_yyyy[i] * ra_x[i];

        tg_xxxyy_yyyz[i] = 2.0 * tg_xyy_yyyz[i] * fxi[i] + tg_xxyy_yyyz[i] * ra_x[i];

        tg_xxxyy_yyzz[i] = 2.0 * tg_xyy_yyzz[i] * fxi[i] + tg_xxyy_yyzz[i] * ra_x[i];

        tg_xxxyy_yzzz[i] = 2.0 * tg_xyy_yzzz[i] * fxi[i] + tg_xxyy_yzzz[i] * ra_x[i];

        tg_xxxyy_zzzz[i] = 2.0 * tg_xyy_zzzz[i] * fxi[i] + tg_xxyy_zzzz[i] * ra_x[i];

        tg_xxxyz_xxxx[i] = tg_xxxz_xxxx[i] * ra_y[i];

        tg_xxxyz_xxxy[i] = tg_xxxy_xxxy[i] * ra_z[i];

        tg_xxxyz_xxxz[i] = tg_xxxz_xxxz[i] * ra_y[i];

        tg_xxxyz_xxyy[i] = tg_xxxy_xxyy[i] * ra_z[i];

        tg_xxxyz_xxyz[i] = tg_xxxz_xxz[i] * fxi[i] + tg_xxxz_xxyz[i] * ra_y[i];

        tg_xxxyz_xxzz[i] = tg_xxxz_xxzz[i] * ra_y[i];

        tg_xxxyz_xyyy[i] = tg_xxxy_xyyy[i] * ra_z[i];

        tg_xxxyz_xyyz[i] = 2.0 * tg_xxxz_xyz[i] * fxi[i] + tg_xxxz_xyyz[i] * ra_y[i];

        tg_xxxyz_xyzz[i] = tg_xxxz_xzz[i] * fxi[i] + tg_xxxz_xyzz[i] * ra_y[i];

        tg_xxxyz_xzzz[i] = tg_xxxz_xzzz[i] * ra_y[i];

        tg_xxxyz_yyyy[i] = tg_xxxy_yyyy[i] * ra_z[i];

        tg_xxxyz_yyyz[i] = 2.0 * tg_xyz_yyyz[i] * fxi[i] + tg_xxyz_yyyz[i] * ra_x[i];

        tg_xxxyz_yyzz[i] = 2.0 * tg_xyz_yyzz[i] * fxi[i] + tg_xxyz_yyzz[i] * ra_x[i];

        tg_xxxyz_yzzz[i] = 2.0 * tg_xyz_yzzz[i] * fxi[i] + tg_xxyz_yzzz[i] * ra_x[i];

        tg_xxxyz_zzzz[i] = tg_xxxz_zzzz[i] * ra_y[i];

        tg_xxxzz_xxxx[i] = tg_xxx_xxxx[i] * fxi[i] + tg_xxxz_xxxx[i] * ra_z[i];

        tg_xxxzz_xxxy[i] = tg_xxx_xxxy[i] * fxi[i] + tg_xxxz_xxxy[i] * ra_z[i];

        tg_xxxzz_xxxz[i] = 2.0 * tg_xzz_xxxz[i] * fxi[i] + 3.0 * tg_xxzz_xxz[i] * fxi[i] + tg_xxzz_xxxz[i] * ra_x[i];

        tg_xxxzz_xxyy[i] = tg_xxx_xxyy[i] * fxi[i] + tg_xxxz_xxyy[i] * ra_z[i];

        tg_xxxzz_xxyz[i] = 2.0 * tg_xzz_xxyz[i] * fxi[i] + 2.0 * tg_xxzz_xyz[i] * fxi[i] + tg_xxzz_xxyz[i] * ra_x[i];

        tg_xxxzz_xxzz[i] = 2.0 * tg_xzz_xxzz[i] * fxi[i] + 2.0 * tg_xxzz_xzz[i] * fxi[i] + tg_xxzz_xxzz[i] * ra_x[i];

        tg_xxxzz_xyyy[i] = tg_xxx_xyyy[i] * fxi[i] + tg_xxxz_xyyy[i] * ra_z[i];

        tg_xxxzz_xyyz[i] = 2.0 * tg_xzz_xyyz[i] * fxi[i] + tg_xxzz_yyz[i] * fxi[i] + tg_xxzz_xyyz[i] * ra_x[i];

        tg_xxxzz_xyzz[i] = 2.0 * tg_xzz_xyzz[i] * fxi[i] + tg_xxzz_yzz[i] * fxi[i] + tg_xxzz_xyzz[i] * ra_x[i];

        tg_xxxzz_xzzz[i] = 2.0 * tg_xzz_xzzz[i] * fxi[i] + tg_xxzz_zzz[i] * fxi[i] + tg_xxzz_xzzz[i] * ra_x[i];

        tg_xxxzz_yyyy[i] = 2.0 * tg_xzz_yyyy[i] * fxi[i] + tg_xxzz_yyyy[i] * ra_x[i];

        tg_xxxzz_yyyz[i] = 2.0 * tg_xzz_yyyz[i] * fxi[i] + tg_xxzz_yyyz[i] * ra_x[i];

        tg_xxxzz_yyzz[i] = 2.0 * tg_xzz_yyzz[i] * fxi[i] + tg_xxzz_yyzz[i] * ra_x[i];

        tg_xxxzz_yzzz[i] = 2.0 * tg_xzz_yzzz[i] * fxi[i] + tg_xxzz_yzzz[i] * ra_x[i];

        tg_xxxzz_zzzz[i] = 2.0 * tg_xzz_zzzz[i] * fxi[i] + tg_xxzz_zzzz[i] * ra_x[i];

        tg_xxyyy_xxxx[i] = 2.0 * tg_xxy_xxxx[i] * fxi[i] + tg_xxyy_xxxx[i] * ra_y[i];

        tg_xxyyy_xxxy[i] = tg_yyy_xxxy[i] * fxi[i] + 3.0 * tg_xyyy_xxy[i] * fxi[i] + tg_xyyy_xxxy[i] * ra_x[i];

        tg_xxyyy_xxxz[i] = 2.0 * tg_xxy_xxxz[i] * fxi[i] + tg_xxyy_xxxz[i] * ra_y[i];

        tg_xxyyy_xxyy[i] = tg_yyy_xxyy[i] * fxi[i] + 2.0 * tg_xyyy_xyy[i] * fxi[i] + tg_xyyy_xxyy[i] * ra_x[i];

        tg_xxyyy_xxyz[i] = tg_yyy_xxyz[i] * fxi[i] + 2.0 * tg_xyyy_xyz[i] * fxi[i] + tg_xyyy_xxyz[i] * ra_x[i];

        tg_xxyyy_xxzz[i] = 2.0 * tg_xxy_xxzz[i] * fxi[i] + tg_xxyy_xxzz[i] * ra_y[i];

        tg_xxyyy_xyyy[i] = tg_yyy_xyyy[i] * fxi[i] + tg_xyyy_yyy[i] * fxi[i] + tg_xyyy_xyyy[i] * ra_x[i];

        tg_xxyyy_xyyz[i] = tg_yyy_xyyz[i] * fxi[i] + tg_xyyy_yyz[i] * fxi[i] + tg_xyyy_xyyz[i] * ra_x[i];

        tg_xxyyy_xyzz[i] = tg_yyy_xyzz[i] * fxi[i] + tg_xyyy_yzz[i] * fxi[i] + tg_xyyy_xyzz[i] * ra_x[i];

        tg_xxyyy_xzzz[i] = 2.0 * tg_xxy_xzzz[i] * fxi[i] + tg_xxyy_xzzz[i] * ra_y[i];

        tg_xxyyy_yyyy[i] = tg_yyy_yyyy[i] * fxi[i] + tg_xyyy_yyyy[i] * ra_x[i];

        tg_xxyyy_yyyz[i] = tg_yyy_yyyz[i] * fxi[i] + tg_xyyy_yyyz[i] * ra_x[i];

        tg_xxyyy_yyzz[i] = tg_yyy_yyzz[i] * fxi[i] + tg_xyyy_yyzz[i] * ra_x[i];

        tg_xxyyy_yzzz[i] = tg_yyy_yzzz[i] * fxi[i] + tg_xyyy_yzzz[i] * ra_x[i];

        tg_xxyyy_zzzz[i] = tg_yyy_zzzz[i] * fxi[i] + tg_xyyy_zzzz[i] * ra_x[i];

        tg_xxyyz_xxxx[i] = tg_xxyy_xxxx[i] * ra_z[i];

        tg_xxyyz_xxxy[i] = tg_xxyy_xxxy[i] * ra_z[i];

        tg_xxyyz_xxxz[i] = tg_xxz_xxxz[i] * fxi[i] + tg_xxyz_xxxz[i] * ra_y[i];

        tg_xxyyz_xxyy[i] = tg_xxyy_xxyy[i] * ra_z[i];

        tg_xxyyz_xxyz[i] = tg_xxyy_xxy[i] * fxi[i] + tg_xxyy_xxyz[i] * ra_z[i];

        tg_xxyyz_xxzz[i] = tg_xxz_xxzz[i] * fxi[i] + tg_xxyz_xxzz[i] * ra_y[i];

        tg_xxyyz_xyyy[i] = tg_xxyy_xyyy[i] * ra_z[i];

        tg_xxyyz_xyyz[i] = tg_xxyy_xyy[i] * fxi[i] + tg_xxyy_xyyz[i] * ra_z[i];

        tg_xxyyz_xyzz[i] = 2.0 * tg_xxyy_xyz[i] * fxi[i] + tg_xxyy_xyzz[i] * ra_z[i];

        tg_xxyyz_xzzz[i] = tg_xxz_xzzz[i] * fxi[i] + tg_xxyz_xzzz[i] * ra_y[i];

        tg_xxyyz_yyyy[i] = tg_xxyy_yyyy[i] * ra_z[i];

        tg_xxyyz_yyyz[i] = tg_yyz_yyyz[i] * fxi[i] + tg_xyyz_yyyz[i] * ra_x[i];

        tg_xxyyz_yyzz[i] = tg_yyz_yyzz[i] * fxi[i] + tg_xyyz_yyzz[i] * ra_x[i];

        tg_xxyyz_yzzz[i] = tg_yyz_yzzz[i] * fxi[i] + tg_xyyz_yzzz[i] * ra_x[i];

        tg_xxyyz_zzzz[i] = tg_yyz_zzzz[i] * fxi[i] + tg_xyyz_zzzz[i] * ra_x[i];

        tg_xxyzz_xxxx[i] = tg_xxzz_xxxx[i] * ra_y[i];

        tg_xxyzz_xxxy[i] = tg_xxzz_xxx[i] * fxi[i] + tg_xxzz_xxxy[i] * ra_y[i];

        tg_xxyzz_xxxz[i] = tg_xxzz_xxxz[i] * ra_y[i];

        tg_xxyzz_xxyy[i] = 2.0 * tg_xxzz_xxy[i] * fxi[i] + tg_xxzz_xxyy[i] * ra_y[i];

        tg_xxyzz_xxyz[i] = tg_xxzz_xxz[i] * fxi[i] + tg_xxzz_xxyz[i] * ra_y[i];

        tg_xxyzz_xxzz[i] = tg_xxzz_xxzz[i] * ra_y[i];

        tg_xxyzz_xyyy[i] = 3.0 * tg_xxzz_xyy[i] * fxi[i] + tg_xxzz_xyyy[i] * ra_y[i];

        tg_xxyzz_xyyz[i] = 2.0 * tg_xxzz_xyz[i] * fxi[i] + tg_xxzz_xyyz[i] * ra_y[i];

        tg_xxyzz_xyzz[i] = tg_xxzz_xzz[i] * fxi[i] + tg_xxzz_xyzz[i] * ra_y[i];

        tg_xxyzz_xzzz[i] = tg_xxzz_xzzz[i] * ra_y[i];

        tg_xxyzz_yyyy[i] = tg_yzz_yyyy[i] * fxi[i] + tg_xyzz_yyyy[i] * ra_x[i];

        tg_xxyzz_yyyz[i] = tg_yzz_yyyz[i] * fxi[i] + tg_xyzz_yyyz[i] * ra_x[i];

        tg_xxyzz_yyzz[i] = tg_yzz_yyzz[i] * fxi[i] + tg_xyzz_yyzz[i] * ra_x[i];

        tg_xxyzz_yzzz[i] = tg_yzz_yzzz[i] * fxi[i] + tg_xyzz_yzzz[i] * ra_x[i];

        tg_xxyzz_zzzz[i] = tg_xxzz_zzzz[i] * ra_y[i];

        tg_xxzzz_xxxx[i] = 2.0 * tg_xxz_xxxx[i] * fxi[i] + tg_xxzz_xxxx[i] * ra_z[i];

        tg_xxzzz_xxxy[i] = 2.0 * tg_xxz_xxxy[i] * fxi[i] + tg_xxzz_xxxy[i] * ra_z[i];

        tg_xxzzz_xxxz[i] = tg_zzz_xxxz[i] * fxi[i] + 3.0 * tg_xzzz_xxz[i] * fxi[i] + tg_xzzz_xxxz[i] * ra_x[i];

        tg_xxzzz_xxyy[i] = 2.0 * tg_xxz_xxyy[i] * fxi[i] + tg_xxzz_xxyy[i] * ra_z[i];

        tg_xxzzz_xxyz[i] = tg_zzz_xxyz[i] * fxi[i] + 2.0 * tg_xzzz_xyz[i] * fxi[i] + tg_xzzz_xxyz[i] * ra_x[i];

        tg_xxzzz_xxzz[i] = tg_zzz_xxzz[i] * fxi[i] + 2.0 * tg_xzzz_xzz[i] * fxi[i] + tg_xzzz_xxzz[i] * ra_x[i];

        tg_xxzzz_xyyy[i] = 2.0 * tg_xxz_xyyy[i] * fxi[i] + tg_xxzz_xyyy[i] * ra_z[i];

        tg_xxzzz_xyyz[i] = tg_zzz_xyyz[i] * fxi[i] + tg_xzzz_yyz[i] * fxi[i] + tg_xzzz_xyyz[i] * ra_x[i];

        tg_xxzzz_xyzz[i] = tg_zzz_xyzz[i] * fxi[i] + tg_xzzz_yzz[i] * fxi[i] + tg_xzzz_xyzz[i] * ra_x[i];

        tg_xxzzz_xzzz[i] = tg_zzz_xzzz[i] * fxi[i] + tg_xzzz_zzz[i] * fxi[i] + tg_xzzz_xzzz[i] * ra_x[i];

        tg_xxzzz_yyyy[i] = tg_zzz_yyyy[i] * fxi[i] + tg_xzzz_yyyy[i] * ra_x[i];

        tg_xxzzz_yyyz[i] = tg_zzz_yyyz[i] * fxi[i] + tg_xzzz_yyyz[i] * ra_x[i];

        tg_xxzzz_yyzz[i] = tg_zzz_yyzz[i] * fxi[i] + tg_xzzz_yyzz[i] * ra_x[i];

        tg_xxzzz_yzzz[i] = tg_zzz_yzzz[i] * fxi[i] + tg_xzzz_yzzz[i] * ra_x[i];

        tg_xxzzz_zzzz[i] = tg_zzz_zzzz[i] * fxi[i] + tg_xzzz_zzzz[i] * ra_x[i];

        tg_xyyyy_xxxx[i] = 4.0 * tg_yyyy_xxx[i] * fxi[i] + tg_yyyy_xxxx[i] * ra_x[i];

        tg_xyyyy_xxxy[i] = 3.0 * tg_yyyy_xxy[i] * fxi[i] + tg_yyyy_xxxy[i] * ra_x[i];

        tg_xyyyy_xxxz[i] = 3.0 * tg_yyyy_xxz[i] * fxi[i] + tg_yyyy_xxxz[i] * ra_x[i];

        tg_xyyyy_xxyy[i] = 2.0 * tg_yyyy_xyy[i] * fxi[i] + tg_yyyy_xxyy[i] * ra_x[i];

        tg_xyyyy_xxyz[i] = 2.0 * tg_yyyy_xyz[i] * fxi[i] + tg_yyyy_xxyz[i] * ra_x[i];

        tg_xyyyy_xxzz[i] = 2.0 * tg_yyyy_xzz[i] * fxi[i] + tg_yyyy_xxzz[i] * ra_x[i];

        tg_xyyyy_xyyy[i] = tg_yyyy_yyy[i] * fxi[i] + tg_yyyy_xyyy[i] * ra_x[i];

        tg_xyyyy_xyyz[i] = tg_yyyy_yyz[i] * fxi[i] + tg_yyyy_xyyz[i] * ra_x[i];

        tg_xyyyy_xyzz[i] = tg_yyyy_yzz[i] * fxi[i] + tg_yyyy_xyzz[i] * ra_x[i];

        tg_xyyyy_xzzz[i] = tg_yyyy_zzz[i] * fxi[i] + tg_yyyy_xzzz[i] * ra_x[i];

        tg_xyyyy_yyyy[i] = tg_yyyy_yyyy[i] * ra_x[i];

        tg_xyyyy_yyyz[i] = tg_yyyy_yyyz[i] * ra_x[i];

        tg_xyyyy_yyzz[i] = tg_yyyy_yyzz[i] * ra_x[i];

        tg_xyyyy_yzzz[i] = tg_yyyy_yzzz[i] * ra_x[i];

        tg_xyyyy_zzzz[i] = tg_yyyy_zzzz[i] * ra_x[i];

        tg_xyyyz_xxxx[i] = tg_xyyy_xxxx[i] * ra_z[i];

        tg_xyyyz_xxxy[i] = tg_xyyy_xxxy[i] * ra_z[i];

        tg_xyyyz_xxxz[i] = 3.0 * tg_yyyz_xxz[i] * fxi[i] + tg_yyyz_xxxz[i] * ra_x[i];

        tg_xyyyz_xxyy[i] = tg_xyyy_xxyy[i] * ra_z[i];

        tg_xyyyz_xxyz[i] = 2.0 * tg_yyyz_xyz[i] * fxi[i] + tg_yyyz_xxyz[i] * ra_x[i];

        tg_xyyyz_xxzz[i] = 2.0 * tg_yyyz_xzz[i] * fxi[i] + tg_yyyz_xxzz[i] * ra_x[i];

        tg_xyyyz_xyyy[i] = tg_xyyy_xyyy[i] * ra_z[i];

        tg_xyyyz_xyyz[i] = tg_yyyz_yyz[i] * fxi[i] + tg_yyyz_xyyz[i] * ra_x[i];

        tg_xyyyz_xyzz[i] = tg_yyyz_yzz[i] * fxi[i] + tg_yyyz_xyzz[i] * ra_x[i];

        tg_xyyyz_xzzz[i] = tg_yyyz_zzz[i] * fxi[i] + tg_yyyz_xzzz[i] * ra_x[i];

        tg_xyyyz_yyyy[i] = tg_yyyz_yyyy[i] * ra_x[i];

        tg_xyyyz_yyyz[i] = tg_yyyz_yyyz[i] * ra_x[i];

        tg_xyyyz_yyzz[i] = tg_yyyz_yyzz[i] * ra_x[i];

        tg_xyyyz_yzzz[i] = tg_yyyz_yzzz[i] * ra_x[i];

        tg_xyyyz_zzzz[i] = tg_yyyz_zzzz[i] * ra_x[i];

        tg_xyyzz_xxxx[i] = 4.0 * tg_yyzz_xxx[i] * fxi[i] + tg_yyzz_xxxx[i] * ra_x[i];

        tg_xyyzz_xxxy[i] = 3.0 * tg_yyzz_xxy[i] * fxi[i] + tg_yyzz_xxxy[i] * ra_x[i];

        tg_xyyzz_xxxz[i] = 3.0 * tg_yyzz_xxz[i] * fxi[i] + tg_yyzz_xxxz[i] * ra_x[i];

        tg_xyyzz_xxyy[i] = 2.0 * tg_yyzz_xyy[i] * fxi[i] + tg_yyzz_xxyy[i] * ra_x[i];

        tg_xyyzz_xxyz[i] = 2.0 * tg_yyzz_xyz[i] * fxi[i] + tg_yyzz_xxyz[i] * ra_x[i];

        tg_xyyzz_xxzz[i] = 2.0 * tg_yyzz_xzz[i] * fxi[i] + tg_yyzz_xxzz[i] * ra_x[i];

        tg_xyyzz_xyyy[i] = tg_yyzz_yyy[i] * fxi[i] + tg_yyzz_xyyy[i] * ra_x[i];

        tg_xyyzz_xyyz[i] = tg_yyzz_yyz[i] * fxi[i] + tg_yyzz_xyyz[i] * ra_x[i];

        tg_xyyzz_xyzz[i] = tg_yyzz_yzz[i] * fxi[i] + tg_yyzz_xyzz[i] * ra_x[i];

        tg_xyyzz_xzzz[i] = tg_yyzz_zzz[i] * fxi[i] + tg_yyzz_xzzz[i] * ra_x[i];

        tg_xyyzz_yyyy[i] = tg_yyzz_yyyy[i] * ra_x[i];

        tg_xyyzz_yyyz[i] = tg_yyzz_yyyz[i] * ra_x[i];

        tg_xyyzz_yyzz[i] = tg_yyzz_yyzz[i] * ra_x[i];

        tg_xyyzz_yzzz[i] = tg_yyzz_yzzz[i] * ra_x[i];

        tg_xyyzz_zzzz[i] = tg_yyzz_zzzz[i] * ra_x[i];

        tg_xyzzz_xxxx[i] = tg_xzzz_xxxx[i] * ra_y[i];

        tg_xyzzz_xxxy[i] = 3.0 * tg_yzzz_xxy[i] * fxi[i] + tg_yzzz_xxxy[i] * ra_x[i];

        tg_xyzzz_xxxz[i] = tg_xzzz_xxxz[i] * ra_y[i];

        tg_xyzzz_xxyy[i] = 2.0 * tg_yzzz_xyy[i] * fxi[i] + tg_yzzz_xxyy[i] * ra_x[i];

        tg_xyzzz_xxyz[i] = 2.0 * tg_yzzz_xyz[i] * fxi[i] + tg_yzzz_xxyz[i] * ra_x[i];

        tg_xyzzz_xxzz[i] = tg_xzzz_xxzz[i] * ra_y[i];

        tg_xyzzz_xyyy[i] = tg_yzzz_yyy[i] * fxi[i] + tg_yzzz_xyyy[i] * ra_x[i];

        tg_xyzzz_xyyz[i] = tg_yzzz_yyz[i] * fxi[i] + tg_yzzz_xyyz[i] * ra_x[i];

        tg_xyzzz_xyzz[i] = tg_yzzz_yzz[i] * fxi[i] + tg_yzzz_xyzz[i] * ra_x[i];

        tg_xyzzz_xzzz[i] = tg_xzzz_xzzz[i] * ra_y[i];

        tg_xyzzz_yyyy[i] = tg_yzzz_yyyy[i] * ra_x[i];

        tg_xyzzz_yyyz[i] = tg_yzzz_yyyz[i] * ra_x[i];

        tg_xyzzz_yyzz[i] = tg_yzzz_yyzz[i] * ra_x[i];

        tg_xyzzz_yzzz[i] = tg_yzzz_yzzz[i] * ra_x[i];

        tg_xyzzz_zzzz[i] = tg_yzzz_zzzz[i] * ra_x[i];

        tg_xzzzz_xxxx[i] = 4.0 * tg_zzzz_xxx[i] * fxi[i] + tg_zzzz_xxxx[i] * ra_x[i];

        tg_xzzzz_xxxy[i] = 3.0 * tg_zzzz_xxy[i] * fxi[i] + tg_zzzz_xxxy[i] * ra_x[i];

        tg_xzzzz_xxxz[i] = 3.0 * tg_zzzz_xxz[i] * fxi[i] + tg_zzzz_xxxz[i] * ra_x[i];

        tg_xzzzz_xxyy[i] = 2.0 * tg_zzzz_xyy[i] * fxi[i] + tg_zzzz_xxyy[i] * ra_x[i];

        tg_xzzzz_xxyz[i] = 2.0 * tg_zzzz_xyz[i] * fxi[i] + tg_zzzz_xxyz[i] * ra_x[i];

        tg_xzzzz_xxzz[i] = 2.0 * tg_zzzz_xzz[i] * fxi[i] + tg_zzzz_xxzz[i] * ra_x[i];

        tg_xzzzz_xyyy[i] = tg_zzzz_yyy[i] * fxi[i] + tg_zzzz_xyyy[i] * ra_x[i];

        tg_xzzzz_xyyz[i] = tg_zzzz_yyz[i] * fxi[i] + tg_zzzz_xyyz[i] * ra_x[i];

        tg_xzzzz_xyzz[i] = tg_zzzz_yzz[i] * fxi[i] + tg_zzzz_xyzz[i] * ra_x[i];

        tg_xzzzz_xzzz[i] = tg_zzzz_zzz[i] * fxi[i] + tg_zzzz_xzzz[i] * ra_x[i];

        tg_xzzzz_yyyy[i] = tg_zzzz_yyyy[i] * ra_x[i];

        tg_xzzzz_yyyz[i] = tg_zzzz_yyyz[i] * ra_x[i];

        tg_xzzzz_yyzz[i] = tg_zzzz_yyzz[i] * ra_x[i];

        tg_xzzzz_yzzz[i] = tg_zzzz_yzzz[i] * ra_x[i];

        tg_xzzzz_zzzz[i] = tg_zzzz_zzzz[i] * ra_x[i];

        tg_yyyyy_xxxx[i] = 4.0 * tg_yyy_xxxx[i] * fxi[i] + tg_yyyy_xxxx[i] * ra_y[i];

        tg_yyyyy_xxxy[i] = 4.0 * tg_yyy_xxxy[i] * fxi[i] + tg_yyyy_xxx[i] * fxi[i] + tg_yyyy_xxxy[i] * ra_y[i];

        tg_yyyyy_xxxz[i] = 4.0 * tg_yyy_xxxz[i] * fxi[i] + tg_yyyy_xxxz[i] * ra_y[i];

        tg_yyyyy_xxyy[i] = 4.0 * tg_yyy_xxyy[i] * fxi[i] + 2.0 * tg_yyyy_xxy[i] * fxi[i] + tg_yyyy_xxyy[i] * ra_y[i];

        tg_yyyyy_xxyz[i] = 4.0 * tg_yyy_xxyz[i] * fxi[i] + tg_yyyy_xxz[i] * fxi[i] + tg_yyyy_xxyz[i] * ra_y[i];

        tg_yyyyy_xxzz[i] = 4.0 * tg_yyy_xxzz[i] * fxi[i] + tg_yyyy_xxzz[i] * ra_y[i];

        tg_yyyyy_xyyy[i] = 4.0 * tg_yyy_xyyy[i] * fxi[i] + 3.0 * tg_yyyy_xyy[i] * fxi[i] + tg_yyyy_xyyy[i] * ra_y[i];

        tg_yyyyy_xyyz[i] = 4.0 * tg_yyy_xyyz[i] * fxi[i] + 2.0 * tg_yyyy_xyz[i] * fxi[i] + tg_yyyy_xyyz[i] * ra_y[i];

        tg_yyyyy_xyzz[i] = 4.0 * tg_yyy_xyzz[i] * fxi[i] + tg_yyyy_xzz[i] * fxi[i] + tg_yyyy_xyzz[i] * ra_y[i];

        tg_yyyyy_xzzz[i] = 4.0 * tg_yyy_xzzz[i] * fxi[i] + tg_yyyy_xzzz[i] * ra_y[i];

        tg_yyyyy_yyyy[i] = 4.0 * tg_yyy_yyyy[i] * fxi[i] + 4.0 * tg_yyyy_yyy[i] * fxi[i] + tg_yyyy_yyyy[i] * ra_y[i];

        tg_yyyyy_yyyz[i] = 4.0 * tg_yyy_yyyz[i] * fxi[i] + 3.0 * tg_yyyy_yyz[i] * fxi[i] + tg_yyyy_yyyz[i] * ra_y[i];

        tg_yyyyy_yyzz[i] = 4.0 * tg_yyy_yyzz[i] * fxi[i] + 2.0 * tg_yyyy_yzz[i] * fxi[i] + tg_yyyy_yyzz[i] * ra_y[i];

        tg_yyyyy_yzzz[i] = 4.0 * tg_yyy_yzzz[i] * fxi[i] + tg_yyyy_zzz[i] * fxi[i] + tg_yyyy_yzzz[i] * ra_y[i];

        tg_yyyyy_zzzz[i] = 4.0 * tg_yyy_zzzz[i] * fxi[i] + tg_yyyy_zzzz[i] * ra_y[i];

        tg_yyyyz_xxxx[i] = tg_yyyy_xxxx[i] * ra_z[i];

        tg_yyyyz_xxxy[i] = tg_yyyy_xxxy[i] * ra_z[i];

        tg_yyyyz_xxxz[i] = 3.0 * tg_yyz_xxxz[i] * fxi[i] + tg_yyyz_xxxz[i] * ra_y[i];

        tg_yyyyz_xxyy[i] = tg_yyyy_xxyy[i] * ra_z[i];

        tg_yyyyz_xxyz[i] = tg_yyyy_xxy[i] * fxi[i] + tg_yyyy_xxyz[i] * ra_z[i];

        tg_yyyyz_xxzz[i] = 3.0 * tg_yyz_xxzz[i] * fxi[i] + tg_yyyz_xxzz[i] * ra_y[i];

        tg_yyyyz_xyyy[i] = tg_yyyy_xyyy[i] * ra_z[i];

        tg_yyyyz_xyyz[i] = tg_yyyy_xyy[i] * fxi[i] + tg_yyyy_xyyz[i] * ra_z[i];

        tg_yyyyz_xyzz[i] = 2.0 * tg_yyyy_xyz[i] * fxi[i] + tg_yyyy_xyzz[i] * ra_z[i];

        tg_yyyyz_xzzz[i] = 3.0 * tg_yyz_xzzz[i] * fxi[i] + tg_yyyz_xzzz[i] * ra_y[i];

        tg_yyyyz_yyyy[i] = tg_yyyy_yyyy[i] * ra_z[i];

        tg_yyyyz_yyyz[i] = tg_yyyy_yyy[i] * fxi[i] + tg_yyyy_yyyz[i] * ra_z[i];

        tg_yyyyz_yyzz[i] = 2.0 * tg_yyyy_yyz[i] * fxi[i] + tg_yyyy_yyzz[i] * ra_z[i];

        tg_yyyyz_yzzz[i] = 3.0 * tg_yyyy_yzz[i] * fxi[i] + tg_yyyy_yzzz[i] * ra_z[i];

        tg_yyyyz_zzzz[i] = 3.0 * tg_yyz_zzzz[i] * fxi[i] + tg_yyyz_zzzz[i] * ra_y[i];

        tg_yyyzz_xxxx[i] = 2.0 * tg_yzz_xxxx[i] * fxi[i] + tg_yyzz_xxxx[i] * ra_y[i];

        tg_yyyzz_xxxy[i] = tg_yyy_xxxy[i] * fxi[i] + tg_yyyz_xxxy[i] * ra_z[i];

        tg_yyyzz_xxxz[i] = 2.0 * tg_yzz_xxxz[i] * fxi[i] + tg_yyzz_xxxz[i] * ra_y[i];

        tg_yyyzz_xxyy[i] = tg_yyy_xxyy[i] * fxi[i] + tg_yyyz_xxyy[i] * ra_z[i];

        tg_yyyzz_xxyz[i] = 2.0 * tg_yzz_xxyz[i] * fxi[i] + tg_yyzz_xxz[i] * fxi[i] + tg_yyzz_xxyz[i] * ra_y[i];

        tg_yyyzz_xxzz[i] = 2.0 * tg_yzz_xxzz[i] * fxi[i] + tg_yyzz_xxzz[i] * ra_y[i];

        tg_yyyzz_xyyy[i] = tg_yyy_xyyy[i] * fxi[i] + tg_yyyz_xyyy[i] * ra_z[i];

        tg_yyyzz_xyyz[i] = 2.0 * tg_yzz_xyyz[i] * fxi[i] + 2.0 * tg_yyzz_xyz[i] * fxi[i] + tg_yyzz_xyyz[i] * ra_y[i];

        tg_yyyzz_xyzz[i] = 2.0 * tg_yzz_xyzz[i] * fxi[i] + tg_yyzz_xzz[i] * fxi[i] + tg_yyzz_xyzz[i] * ra_y[i];

        tg_yyyzz_xzzz[i] = 2.0 * tg_yzz_xzzz[i] * fxi[i] + tg_yyzz_xzzz[i] * ra_y[i];

        tg_yyyzz_yyyy[i] = tg_yyy_yyyy[i] * fxi[i] + tg_yyyz_yyyy[i] * ra_z[i];

        tg_yyyzz_yyyz[i] = 2.0 * tg_yzz_yyyz[i] * fxi[i] + 3.0 * tg_yyzz_yyz[i] * fxi[i] + tg_yyzz_yyyz[i] * ra_y[i];

        tg_yyyzz_yyzz[i] = 2.0 * tg_yzz_yyzz[i] * fxi[i] + 2.0 * tg_yyzz_yzz[i] * fxi[i] + tg_yyzz_yyzz[i] * ra_y[i];

        tg_yyyzz_yzzz[i] = 2.0 * tg_yzz_yzzz[i] * fxi[i] + tg_yyzz_zzz[i] * fxi[i] + tg_yyzz_yzzz[i] * ra_y[i];

        tg_yyyzz_zzzz[i] = 2.0 * tg_yzz_zzzz[i] * fxi[i] + tg_yyzz_zzzz[i] * ra_y[i];

        tg_yyzzz_xxxx[i] = tg_zzz_xxxx[i] * fxi[i] + tg_yzzz_xxxx[i] * ra_y[i];

        tg_yyzzz_xxxy[i] = 2.0 * tg_yyz_xxxy[i] * fxi[i] + tg_yyzz_xxxy[i] * ra_z[i];

        tg_yyzzz_xxxz[i] = tg_zzz_xxxz[i] * fxi[i] + tg_yzzz_xxxz[i] * ra_y[i];

        tg_yyzzz_xxyy[i] = 2.0 * tg_yyz_xxyy[i] * fxi[i] + tg_yyzz_xxyy[i] * ra_z[i];

        tg_yyzzz_xxyz[i] = tg_zzz_xxyz[i] * fxi[i] + tg_yzzz_xxz[i] * fxi[i] + tg_yzzz_xxyz[i] * ra_y[i];

        tg_yyzzz_xxzz[i] = tg_zzz_xxzz[i] * fxi[i] + tg_yzzz_xxzz[i] * ra_y[i];

        tg_yyzzz_xyyy[i] = 2.0 * tg_yyz_xyyy[i] * fxi[i] + tg_yyzz_xyyy[i] * ra_z[i];

        tg_yyzzz_xyyz[i] = tg_zzz_xyyz[i] * fxi[i] + 2.0 * tg_yzzz_xyz[i] * fxi[i] + tg_yzzz_xyyz[i] * ra_y[i];

        tg_yyzzz_xyzz[i] = tg_zzz_xyzz[i] * fxi[i] + tg_yzzz_xzz[i] * fxi[i] + tg_yzzz_xyzz[i] * ra_y[i];

        tg_yyzzz_xzzz[i] = tg_zzz_xzzz[i] * fxi[i] + tg_yzzz_xzzz[i] * ra_y[i];

        tg_yyzzz_yyyy[i] = 2.0 * tg_yyz_yyyy[i] * fxi[i] + tg_yyzz_yyyy[i] * ra_z[i];

        tg_yyzzz_yyyz[i] = tg_zzz_yyyz[i] * fxi[i] + 3.0 * tg_yzzz_yyz[i] * fxi[i] + tg_yzzz_yyyz[i] * ra_y[i];

        tg_yyzzz_yyzz[i] = tg_zzz_yyzz[i] * fxi[i] + 2.0 * tg_yzzz_yzz[i] * fxi[i] + tg_yzzz_yyzz[i] * ra_y[i];

        tg_yyzzz_yzzz[i] = tg_zzz_yzzz[i] * fxi[i] + tg_yzzz_zzz[i] * fxi[i] + tg_yzzz_yzzz[i] * ra_y[i];

        tg_yyzzz_zzzz[i] = tg_zzz_zzzz[i] * fxi[i] + tg_yzzz_zzzz[i] * ra_y[i];

        tg_yzzzz_xxxx[i] = tg_zzzz_xxxx[i] * ra_y[i];

        tg_yzzzz_xxxy[i] = tg_zzzz_xxx[i] * fxi[i] + tg_zzzz_xxxy[i] * ra_y[i];

        tg_yzzzz_xxxz[i] = tg_zzzz_xxxz[i] * ra_y[i];

        tg_yzzzz_xxyy[i] = 2.0 * tg_zzzz_xxy[i] * fxi[i] + tg_zzzz_xxyy[i] * ra_y[i];

        tg_yzzzz_xxyz[i] = tg_zzzz_xxz[i] * fxi[i] + tg_zzzz_xxyz[i] * ra_y[i];

        tg_yzzzz_xxzz[i] = tg_zzzz_xxzz[i] * ra_y[i];

        tg_yzzzz_xyyy[i] = 3.0 * tg_zzzz_xyy[i] * fxi[i] + tg_zzzz_xyyy[i] * ra_y[i];

        tg_yzzzz_xyyz[i] = 2.0 * tg_zzzz_xyz[i] * fxi[i] + tg_zzzz_xyyz[i] * ra_y[i];

        tg_yzzzz_xyzz[i] = tg_zzzz_xzz[i] * fxi[i] + tg_zzzz_xyzz[i] * ra_y[i];

        tg_yzzzz_xzzz[i] = tg_zzzz_xzzz[i] * ra_y[i];

        tg_yzzzz_yyyy[i] = 4.0 * tg_zzzz_yyy[i] * fxi[i] + tg_zzzz_yyyy[i] * ra_y[i];

        tg_yzzzz_yyyz[i] = 3.0 * tg_zzzz_yyz[i] * fxi[i] + tg_zzzz_yyyz[i] * ra_y[i];

        tg_yzzzz_yyzz[i] = 2.0 * tg_zzzz_yzz[i] * fxi[i] + tg_zzzz_yyzz[i] * ra_y[i];

        tg_yzzzz_yzzz[i] = tg_zzzz_zzz[i] * fxi[i] + tg_zzzz_yzzz[i] * ra_y[i];

        tg_yzzzz_zzzz[i] = tg_zzzz_zzzz[i] * ra_y[i];

        tg_zzzzz_xxxx[i] = 4.0 * tg_zzz_xxxx[i] * fxi[i] + tg_zzzz_xxxx[i] * ra_z[i];

        tg_zzzzz_xxxy[i] = 4.0 * tg_zzz_xxxy[i] * fxi[i] + tg_zzzz_xxxy[i] * ra_z[i];

        tg_zzzzz_xxxz[i] = 4.0 * tg_zzz_xxxz[i] * fxi[i] + tg_zzzz_xxx[i] * fxi[i] + tg_zzzz_xxxz[i] * ra_z[i];

        tg_zzzzz_xxyy[i] = 4.0 * tg_zzz_xxyy[i] * fxi[i] + tg_zzzz_xxyy[i] * ra_z[i];

        tg_zzzzz_xxyz[i] = 4.0 * tg_zzz_xxyz[i] * fxi[i] + tg_zzzz_xxy[i] * fxi[i] + tg_zzzz_xxyz[i] * ra_z[i];

        tg_zzzzz_xxzz[i] = 4.0 * tg_zzz_xxzz[i] * fxi[i] + 2.0 * tg_zzzz_xxz[i] * fxi[i] + tg_zzzz_xxzz[i] * ra_z[i];

        tg_zzzzz_xyyy[i] = 4.0 * tg_zzz_xyyy[i] * fxi[i] + tg_zzzz_xyyy[i] * ra_z[i];

        tg_zzzzz_xyyz[i] = 4.0 * tg_zzz_xyyz[i] * fxi[i] + tg_zzzz_xyy[i] * fxi[i] + tg_zzzz_xyyz[i] * ra_z[i];

        tg_zzzzz_xyzz[i] = 4.0 * tg_zzz_xyzz[i] * fxi[i] + 2.0 * tg_zzzz_xyz[i] * fxi[i] + tg_zzzz_xyzz[i] * ra_z[i];

        tg_zzzzz_xzzz[i] = 4.0 * tg_zzz_xzzz[i] * fxi[i] + 3.0 * tg_zzzz_xzz[i] * fxi[i] + tg_zzzz_xzzz[i] * ra_z[i];

        tg_zzzzz_yyyy[i] = 4.0 * tg_zzz_yyyy[i] * fxi[i] + tg_zzzz_yyyy[i] * ra_z[i];

        tg_zzzzz_yyyz[i] = 4.0 * tg_zzz_yyyz[i] * fxi[i] + tg_zzzz_yyy[i] * fxi[i] + tg_zzzz_yyyz[i] * ra_z[i];

        tg_zzzzz_yyzz[i] = 4.0 * tg_zzz_yyzz[i] * fxi[i] + 2.0 * tg_zzzz_yyz[i] * fxi[i] + tg_zzzz_yyzz[i] * ra_z[i];

        tg_zzzzz_yzzz[i] = 4.0 * tg_zzz_yzzz[i] * fxi[i] + 3.0 * tg_zzzz_yzz[i] * fxi[i] + tg_zzzz_yzzz[i] * ra_z[i];

        tg_zzzzz_zzzz[i] = 4.0 * tg_zzz_zzzz[i] * fxi[i] + 4.0 * tg_zzzz_zzz[i] * fxi[i] + tg_zzzz_zzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

