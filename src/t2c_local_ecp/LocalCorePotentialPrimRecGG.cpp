#include "LocalCorePotentialPrimRecGG.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_gg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gg,
                                  const size_t idx_dg,
                                  const size_t idx_ff,
                                  const size_t idx_fg,
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

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx = pbuffer.data(idx_dg);

    auto tg_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto tg_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto tg_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto tg_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto tg_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto tg_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto tg_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto tg_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto tg_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto tg_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto tg_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto tg_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto tg_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto tg_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto tg_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto tg_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto tg_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto tg_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto tg_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto tg_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto tg_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto tg_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto tg_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto tg_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto tg_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto tg_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto tg_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto tg_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto tg_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto tg_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto tg_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto tg_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto tg_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto tg_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto tg_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto tg_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto tg_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto tg_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto tg_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto tg_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto tg_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto tg_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto tg_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto tg_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto tg_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto tg_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto tg_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto tg_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto tg_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto tg_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto tg_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto tg_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto tg_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto tg_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto tg_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto tg_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto tg_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto tg_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto tg_zz_zzzz = pbuffer.data(idx_dg + 89);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx = pbuffer.data(idx_ff);

    auto tg_xxx_xxy = pbuffer.data(idx_ff + 1);

    auto tg_xxx_xxz = pbuffer.data(idx_ff + 2);

    auto tg_xxx_xyy = pbuffer.data(idx_ff + 3);

    auto tg_xxx_xyz = pbuffer.data(idx_ff + 4);

    auto tg_xxx_xzz = pbuffer.data(idx_ff + 5);

    auto tg_xxx_yyy = pbuffer.data(idx_ff + 6);

    auto tg_xxx_yyz = pbuffer.data(idx_ff + 7);

    auto tg_xxx_yzz = pbuffer.data(idx_ff + 8);

    auto tg_xxx_zzz = pbuffer.data(idx_ff + 9);

    auto tg_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto tg_xxz_xyz = pbuffer.data(idx_ff + 24);

    auto tg_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto tg_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto tg_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto tg_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto tg_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto tg_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto tg_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto tg_xzz_xxz = pbuffer.data(idx_ff + 52);

    auto tg_xzz_xyz = pbuffer.data(idx_ff + 54);

    auto tg_xzz_xzz = pbuffer.data(idx_ff + 55);

    auto tg_xzz_yyz = pbuffer.data(idx_ff + 57);

    auto tg_xzz_yzz = pbuffer.data(idx_ff + 58);

    auto tg_xzz_zzz = pbuffer.data(idx_ff + 59);

    auto tg_yyy_xxx = pbuffer.data(idx_ff + 60);

    auto tg_yyy_xxy = pbuffer.data(idx_ff + 61);

    auto tg_yyy_xxz = pbuffer.data(idx_ff + 62);

    auto tg_yyy_xyy = pbuffer.data(idx_ff + 63);

    auto tg_yyy_xyz = pbuffer.data(idx_ff + 64);

    auto tg_yyy_xzz = pbuffer.data(idx_ff + 65);

    auto tg_yyy_yyy = pbuffer.data(idx_ff + 66);

    auto tg_yyy_yyz = pbuffer.data(idx_ff + 67);

    auto tg_yyy_yzz = pbuffer.data(idx_ff + 68);

    auto tg_yyy_zzz = pbuffer.data(idx_ff + 69);

    auto tg_yyz_xxz = pbuffer.data(idx_ff + 72);

    auto tg_yyz_xyz = pbuffer.data(idx_ff + 74);

    auto tg_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto tg_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto tg_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto tg_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto tg_yzz_xxy = pbuffer.data(idx_ff + 81);

    auto tg_yzz_xxz = pbuffer.data(idx_ff + 82);

    auto tg_yzz_xyy = pbuffer.data(idx_ff + 83);

    auto tg_yzz_xyz = pbuffer.data(idx_ff + 84);

    auto tg_yzz_xzz = pbuffer.data(idx_ff + 85);

    auto tg_yzz_yyy = pbuffer.data(idx_ff + 86);

    auto tg_yzz_yyz = pbuffer.data(idx_ff + 87);

    auto tg_yzz_yzz = pbuffer.data(idx_ff + 88);

    auto tg_yzz_zzz = pbuffer.data(idx_ff + 89);

    auto tg_zzz_xxx = pbuffer.data(idx_ff + 90);

    auto tg_zzz_xxy = pbuffer.data(idx_ff + 91);

    auto tg_zzz_xxz = pbuffer.data(idx_ff + 92);

    auto tg_zzz_xyy = pbuffer.data(idx_ff + 93);

    auto tg_zzz_xyz = pbuffer.data(idx_ff + 94);

    auto tg_zzz_xzz = pbuffer.data(idx_ff + 95);

    auto tg_zzz_yyy = pbuffer.data(idx_ff + 96);

    auto tg_zzz_yyz = pbuffer.data(idx_ff + 97);

    auto tg_zzz_yzz = pbuffer.data(idx_ff + 98);

    auto tg_zzz_zzz = pbuffer.data(idx_ff + 99);

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

    auto tg_xxy_xxxy = pbuffer.data(idx_fg + 16);

    auto tg_xxy_xxxz = pbuffer.data(idx_fg + 17);

    auto tg_xxy_xxyy = pbuffer.data(idx_fg + 18);

    auto tg_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto tg_xxy_xyyy = pbuffer.data(idx_fg + 21);

    auto tg_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto tg_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto tg_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto tg_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto tg_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto tg_xxz_xxxx = pbuffer.data(idx_fg + 30);

    auto tg_xxz_xxxy = pbuffer.data(idx_fg + 31);

    auto tg_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto tg_xxz_xxyy = pbuffer.data(idx_fg + 33);

    auto tg_xxz_xxyz = pbuffer.data(idx_fg + 34);

    auto tg_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto tg_xxz_xyyy = pbuffer.data(idx_fg + 36);

    auto tg_xxz_xyyz = pbuffer.data(idx_fg + 37);

    auto tg_xxz_xyzz = pbuffer.data(idx_fg + 38);

    auto tg_xxz_xzzz = pbuffer.data(idx_fg + 39);

    auto tg_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto tg_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto tg_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto tg_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto tg_xyy_xxxx = pbuffer.data(idx_fg + 45);

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

    auto tg_xzz_xxxx = pbuffer.data(idx_fg + 75);

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

    auto tg_yyz_xxyz = pbuffer.data(idx_fg + 109);

    auto tg_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto tg_yyz_xyyy = pbuffer.data(idx_fg + 111);

    auto tg_yyz_xyyz = pbuffer.data(idx_fg + 112);

    auto tg_yyz_xyzz = pbuffer.data(idx_fg + 113);

    auto tg_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto tg_yyz_yyyy = pbuffer.data(idx_fg + 115);

    auto tg_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto tg_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto tg_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto tg_yyz_zzzz = pbuffer.data(idx_fg + 119);

    auto tg_yzz_xxxx = pbuffer.data(idx_fg + 120);

    auto tg_yzz_xxxy = pbuffer.data(idx_fg + 121);

    auto tg_yzz_xxxz = pbuffer.data(idx_fg + 122);

    auto tg_yzz_xxyy = pbuffer.data(idx_fg + 123);

    auto tg_yzz_xxyz = pbuffer.data(idx_fg + 124);

    auto tg_yzz_xxzz = pbuffer.data(idx_fg + 125);

    auto tg_yzz_xyyy = pbuffer.data(idx_fg + 126);

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

    // Set up components of targeted buffer : GG

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

    auto tg_xxxy_xxyz = pbuffer.data(idx_gg + 19);

    auto tg_xxxy_xxzz = pbuffer.data(idx_gg + 20);

    auto tg_xxxy_xyyy = pbuffer.data(idx_gg + 21);

    auto tg_xxxy_xyyz = pbuffer.data(idx_gg + 22);

    auto tg_xxxy_xyzz = pbuffer.data(idx_gg + 23);

    auto tg_xxxy_xzzz = pbuffer.data(idx_gg + 24);

    auto tg_xxxy_yyyy = pbuffer.data(idx_gg + 25);

    auto tg_xxxy_yyyz = pbuffer.data(idx_gg + 26);

    auto tg_xxxy_yyzz = pbuffer.data(idx_gg + 27);

    auto tg_xxxy_yzzz = pbuffer.data(idx_gg + 28);

    auto tg_xxxy_zzzz = pbuffer.data(idx_gg + 29);

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

    auto tg_xxxz_yyyy = pbuffer.data(idx_gg + 40);

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

    auto tg_xxyz_xxxx = pbuffer.data(idx_gg + 60);

    auto tg_xxyz_xxxy = pbuffer.data(idx_gg + 61);

    auto tg_xxyz_xxxz = pbuffer.data(idx_gg + 62);

    auto tg_xxyz_xxyy = pbuffer.data(idx_gg + 63);

    auto tg_xxyz_xxyz = pbuffer.data(idx_gg + 64);

    auto tg_xxyz_xxzz = pbuffer.data(idx_gg + 65);

    auto tg_xxyz_xyyy = pbuffer.data(idx_gg + 66);

    auto tg_xxyz_xyyz = pbuffer.data(idx_gg + 67);

    auto tg_xxyz_xyzz = pbuffer.data(idx_gg + 68);

    auto tg_xxyz_xzzz = pbuffer.data(idx_gg + 69);

    auto tg_xxyz_yyyy = pbuffer.data(idx_gg + 70);

    auto tg_xxyz_yyyz = pbuffer.data(idx_gg + 71);

    auto tg_xxyz_yyzz = pbuffer.data(idx_gg + 72);

    auto tg_xxyz_yzzz = pbuffer.data(idx_gg + 73);

    auto tg_xxyz_zzzz = pbuffer.data(idx_gg + 74);

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

    auto tg_xyyy_xxxz = pbuffer.data(idx_gg + 92);

    auto tg_xyyy_xxyy = pbuffer.data(idx_gg + 93);

    auto tg_xyyy_xxyz = pbuffer.data(idx_gg + 94);

    auto tg_xyyy_xxzz = pbuffer.data(idx_gg + 95);

    auto tg_xyyy_xyyy = pbuffer.data(idx_gg + 96);

    auto tg_xyyy_xyyz = pbuffer.data(idx_gg + 97);

    auto tg_xyyy_xyzz = pbuffer.data(idx_gg + 98);

    auto tg_xyyy_xzzz = pbuffer.data(idx_gg + 99);

    auto tg_xyyy_yyyy = pbuffer.data(idx_gg + 100);

    auto tg_xyyy_yyyz = pbuffer.data(idx_gg + 101);

    auto tg_xyyy_yyzz = pbuffer.data(idx_gg + 102);

    auto tg_xyyy_yzzz = pbuffer.data(idx_gg + 103);

    auto tg_xyyy_zzzz = pbuffer.data(idx_gg + 104);

    auto tg_xyyz_xxxx = pbuffer.data(idx_gg + 105);

    auto tg_xyyz_xxxy = pbuffer.data(idx_gg + 106);

    auto tg_xyyz_xxxz = pbuffer.data(idx_gg + 107);

    auto tg_xyyz_xxyy = pbuffer.data(idx_gg + 108);

    auto tg_xyyz_xxyz = pbuffer.data(idx_gg + 109);

    auto tg_xyyz_xxzz = pbuffer.data(idx_gg + 110);

    auto tg_xyyz_xyyy = pbuffer.data(idx_gg + 111);

    auto tg_xyyz_xyyz = pbuffer.data(idx_gg + 112);

    auto tg_xyyz_xyzz = pbuffer.data(idx_gg + 113);

    auto tg_xyyz_xzzz = pbuffer.data(idx_gg + 114);

    auto tg_xyyz_yyyy = pbuffer.data(idx_gg + 115);

    auto tg_xyyz_yyyz = pbuffer.data(idx_gg + 116);

    auto tg_xyyz_yyzz = pbuffer.data(idx_gg + 117);

    auto tg_xyyz_yzzz = pbuffer.data(idx_gg + 118);

    auto tg_xyyz_zzzz = pbuffer.data(idx_gg + 119);

    auto tg_xyzz_xxxx = pbuffer.data(idx_gg + 120);

    auto tg_xyzz_xxxy = pbuffer.data(idx_gg + 121);

    auto tg_xyzz_xxxz = pbuffer.data(idx_gg + 122);

    auto tg_xyzz_xxyy = pbuffer.data(idx_gg + 123);

    auto tg_xyzz_xxyz = pbuffer.data(idx_gg + 124);

    auto tg_xyzz_xxzz = pbuffer.data(idx_gg + 125);

    auto tg_xyzz_xyyy = pbuffer.data(idx_gg + 126);

    auto tg_xyzz_xyyz = pbuffer.data(idx_gg + 127);

    auto tg_xyzz_xyzz = pbuffer.data(idx_gg + 128);

    auto tg_xyzz_xzzz = pbuffer.data(idx_gg + 129);

    auto tg_xyzz_yyyy = pbuffer.data(idx_gg + 130);

    auto tg_xyzz_yyyz = pbuffer.data(idx_gg + 131);

    auto tg_xyzz_yyzz = pbuffer.data(idx_gg + 132);

    auto tg_xyzz_yzzz = pbuffer.data(idx_gg + 133);

    auto tg_xyzz_zzzz = pbuffer.data(idx_gg + 134);

    auto tg_xzzz_xxxx = pbuffer.data(idx_gg + 135);

    auto tg_xzzz_xxxy = pbuffer.data(idx_gg + 136);

    auto tg_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto tg_xzzz_xxyy = pbuffer.data(idx_gg + 138);

    auto tg_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto tg_xzzz_xxzz = pbuffer.data(idx_gg + 140);

    auto tg_xzzz_xyyy = pbuffer.data(idx_gg + 141);

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

    auto tg_yyyz_xxxx = pbuffer.data(idx_gg + 165);

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

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xx_xxxx, tg_xx_xxxy, tg_xx_xxxz, tg_xx_xxyy, tg_xx_xxyz, tg_xx_xxzz, tg_xx_xyyy, tg_xx_xyyz, tg_xx_xyzz, tg_xx_xzzz, tg_xx_yyyy, tg_xx_yyyz, tg_xx_yyzz, tg_xx_yzzz, tg_xx_zzzz, tg_xxx_xxx, tg_xxx_xxxx, tg_xxx_xxxy, tg_xxx_xxxz, tg_xxx_xxy, tg_xxx_xxyy, tg_xxx_xxyz, tg_xxx_xxz, tg_xxx_xxzz, tg_xxx_xyy, tg_xxx_xyyy, tg_xxx_xyyz, tg_xxx_xyz, tg_xxx_xyzz, tg_xxx_xzz, tg_xxx_xzzz, tg_xxx_yyy, tg_xxx_yyyy, tg_xxx_yyyz, tg_xxx_yyz, tg_xxx_yyzz, tg_xxx_yzz, tg_xxx_yzzz, tg_xxx_zzz, tg_xxx_zzzz, tg_xxxx_xxxx, tg_xxxx_xxxy, tg_xxxx_xxxz, tg_xxxx_xxyy, tg_xxxx_xxyz, tg_xxxx_xxzz, tg_xxxx_xyyy, tg_xxxx_xyyz, tg_xxxx_xyzz, tg_xxxx_xzzz, tg_xxxx_yyyy, tg_xxxx_yyyz, tg_xxxx_yyzz, tg_xxxx_yzzz, tg_xxxx_zzzz, tg_xxxy_xxxx, tg_xxxy_xxxy, tg_xxxy_xxxz, tg_xxxy_xxyy, tg_xxxy_xxyz, tg_xxxy_xxzz, tg_xxxy_xyyy, tg_xxxy_xyyz, tg_xxxy_xyzz, tg_xxxy_xzzz, tg_xxxy_yyyy, tg_xxxy_yyyz, tg_xxxy_yyzz, tg_xxxy_yzzz, tg_xxxy_zzzz, tg_xxxz_xxxx, tg_xxxz_xxxy, tg_xxxz_xxxz, tg_xxxz_xxyy, tg_xxxz_xxyz, tg_xxxz_xxzz, tg_xxxz_xyyy, tg_xxxz_xyyz, tg_xxxz_xyzz, tg_xxxz_xzzz, tg_xxxz_yyyy, tg_xxxz_yyyz, tg_xxxz_yyzz, tg_xxxz_yzzz, tg_xxxz_zzzz, tg_xxy_xxxx, tg_xxy_xxxy, tg_xxy_xxxz, tg_xxy_xxyy, tg_xxy_xxzz, tg_xxy_xyyy, tg_xxy_xzzz, tg_xxy_yyyy, tg_xxy_yyyz, tg_xxy_yyzz, tg_xxy_yzzz, tg_xxyy_xxxx, tg_xxyy_xxxy, tg_xxyy_xxxz, tg_xxyy_xxyy, tg_xxyy_xxyz, tg_xxyy_xxzz, tg_xxyy_xyyy, tg_xxyy_xyyz, tg_xxyy_xyzz, tg_xxyy_xzzz, tg_xxyy_yyyy, tg_xxyy_yyyz, tg_xxyy_yyzz, tg_xxyy_yzzz, tg_xxyy_zzzz, tg_xxyz_xxxx, tg_xxyz_xxxy, tg_xxyz_xxxz, tg_xxyz_xxyy, tg_xxyz_xxyz, tg_xxyz_xxzz, tg_xxyz_xyyy, tg_xxyz_xyyz, tg_xxyz_xyzz, tg_xxyz_xzzz, tg_xxyz_yyyy, tg_xxyz_yyyz, tg_xxyz_yyzz, tg_xxyz_yzzz, tg_xxyz_zzzz, tg_xxz_xxxx, tg_xxz_xxxy, tg_xxz_xxxz, tg_xxz_xxyy, tg_xxz_xxyz, tg_xxz_xxz, tg_xxz_xxzz, tg_xxz_xyyy, tg_xxz_xyyz, tg_xxz_xyz, tg_xxz_xyzz, tg_xxz_xzz, tg_xxz_xzzz, tg_xxz_yyyz, tg_xxz_yyzz, tg_xxz_yzzz, tg_xxz_zzzz, tg_xxzz_xxxx, tg_xxzz_xxxy, tg_xxzz_xxxz, tg_xxzz_xxyy, tg_xxzz_xxyz, tg_xxzz_xxzz, tg_xxzz_xyyy, tg_xxzz_xyyz, tg_xxzz_xyzz, tg_xxzz_xzzz, tg_xxzz_yyyy, tg_xxzz_yyyz, tg_xxzz_yyzz, tg_xxzz_yzzz, tg_xxzz_zzzz, tg_xy_yyyy, tg_xy_yyyz, tg_xy_yyzz, tg_xy_yzzz, tg_xyy_xxxx, tg_xyy_xxxy, tg_xyy_xxy, tg_xyy_xxyy, tg_xyy_xxyz, tg_xyy_xyy, tg_xyy_xyyy, tg_xyy_xyyz, tg_xyy_xyz, tg_xyy_xyzz, tg_xyy_yyy, tg_xyy_yyyy, tg_xyy_yyyz, tg_xyy_yyz, tg_xyy_yyzz, tg_xyy_yzz, tg_xyy_yzzz, tg_xyy_zzzz, tg_xyyy_xxxx, tg_xyyy_xxxy, tg_xyyy_xxxz, tg_xyyy_xxyy, tg_xyyy_xxyz, tg_xyyy_xxzz, tg_xyyy_xyyy, tg_xyyy_xyyz, tg_xyyy_xyzz, tg_xyyy_xzzz, tg_xyyy_yyyy, tg_xyyy_yyyz, tg_xyyy_yyzz, tg_xyyy_yzzz, tg_xyyy_zzzz, tg_xyyz_xxxx, tg_xyyz_xxxy, tg_xyyz_xxxz, tg_xyyz_xxyy, tg_xyyz_xxyz, tg_xyyz_xxzz, tg_xyyz_xyyy, tg_xyyz_xyyz, tg_xyyz_xyzz, tg_xyyz_xzzz, tg_xyyz_yyyy, tg_xyyz_yyyz, tg_xyyz_yyzz, tg_xyyz_yzzz, tg_xyyz_zzzz, tg_xyz_yyyz, tg_xyz_yyzz, tg_xyz_yzzz, tg_xyzz_xxxx, tg_xyzz_xxxy, tg_xyzz_xxxz, tg_xyzz_xxyy, tg_xyzz_xxyz, tg_xyzz_xxzz, tg_xyzz_xyyy, tg_xyzz_xyyz, tg_xyzz_xyzz, tg_xyzz_xzzz, tg_xyzz_yyyy, tg_xyzz_yyyz, tg_xyzz_yyzz, tg_xyzz_yzzz, tg_xyzz_zzzz, tg_xz_yyyz, tg_xz_yyzz, tg_xz_yzzz, tg_xz_zzzz, tg_xzz_xxxx, tg_xzz_xxxz, tg_xzz_xxyz, tg_xzz_xxz, tg_xzz_xxzz, tg_xzz_xyyz, tg_xzz_xyz, tg_xzz_xyzz, tg_xzz_xzz, tg_xzz_xzzz, tg_xzz_yyyy, tg_xzz_yyyz, tg_xzz_yyz, tg_xzz_yyzz, tg_xzz_yzz, tg_xzz_yzzz, tg_xzz_zzz, tg_xzz_zzzz, tg_xzzz_xxxx, tg_xzzz_xxxy, tg_xzzz_xxxz, tg_xzzz_xxyy, tg_xzzz_xxyz, tg_xzzz_xxzz, tg_xzzz_xyyy, tg_xzzz_xyyz, tg_xzzz_xyzz, tg_xzzz_xzzz, tg_xzzz_yyyy, tg_xzzz_yyyz, tg_xzzz_yyzz, tg_xzzz_yzzz, tg_xzzz_zzzz, tg_yy_xxxx, tg_yy_xxxy, tg_yy_xxxz, tg_yy_xxyy, tg_yy_xxyz, tg_yy_xxzz, tg_yy_xyyy, tg_yy_xyyz, tg_yy_xyzz, tg_yy_xzzz, tg_yy_yyyy, tg_yy_yyyz, tg_yy_yyzz, tg_yy_yzzz, tg_yy_zzzz, tg_yyy_xxx, tg_yyy_xxxx, tg_yyy_xxxy, tg_yyy_xxxz, tg_yyy_xxy, tg_yyy_xxyy, tg_yyy_xxyz, tg_yyy_xxz, tg_yyy_xxzz, tg_yyy_xyy, tg_yyy_xyyy, tg_yyy_xyyz, tg_yyy_xyz, tg_yyy_xyzz, tg_yyy_xzz, tg_yyy_xzzz, tg_yyy_yyy, tg_yyy_yyyy, tg_yyy_yyyz, tg_yyy_yyz, tg_yyy_yyzz, tg_yyy_yzz, tg_yyy_yzzz, tg_yyy_zzz, tg_yyy_zzzz, tg_yyyy_xxxx, tg_yyyy_xxxy, tg_yyyy_xxxz, tg_yyyy_xxyy, tg_yyyy_xxyz, tg_yyyy_xxzz, tg_yyyy_xyyy, tg_yyyy_xyyz, tg_yyyy_xyzz, tg_yyyy_xzzz, tg_yyyy_yyyy, tg_yyyy_yyyz, tg_yyyy_yyzz, tg_yyyy_yzzz, tg_yyyy_zzzz, tg_yyyz_xxxx, tg_yyyz_xxxy, tg_yyyz_xxxz, tg_yyyz_xxyy, tg_yyyz_xxyz, tg_yyyz_xxzz, tg_yyyz_xyyy, tg_yyyz_xyyz, tg_yyyz_xyzz, tg_yyyz_xzzz, tg_yyyz_yyyy, tg_yyyz_yyyz, tg_yyyz_yyzz, tg_yyyz_yzzz, tg_yyyz_zzzz, tg_yyz_xxxy, tg_yyz_xxxz, tg_yyz_xxyy, tg_yyz_xxyz, tg_yyz_xxz, tg_yyz_xxzz, tg_yyz_xyyy, tg_yyz_xyyz, tg_yyz_xyz, tg_yyz_xyzz, tg_yyz_xzz, tg_yyz_xzzz, tg_yyz_yyyy, tg_yyz_yyyz, tg_yyz_yyz, tg_yyz_yyzz, tg_yyz_yzz, tg_yyz_yzzz, tg_yyz_zzz, tg_yyz_zzzz, tg_yyzz_xxxx, tg_yyzz_xxxy, tg_yyzz_xxxz, tg_yyzz_xxyy, tg_yyzz_xxyz, tg_yyzz_xxzz, tg_yyzz_xyyy, tg_yyzz_xyyz, tg_yyzz_xyzz, tg_yyzz_xzzz, tg_yyzz_yyyy, tg_yyzz_yyyz, tg_yyzz_yyzz, tg_yyzz_yzzz, tg_yyzz_zzzz, tg_yz_xxxz, tg_yz_xxzz, tg_yz_xzzz, tg_yz_yyyz, tg_yz_yyzz, tg_yz_yzzz, tg_yz_zzzz, tg_yzz_xxxx, tg_yzz_xxxy, tg_yzz_xxxz, tg_yzz_xxy, tg_yzz_xxyy, tg_yzz_xxyz, tg_yzz_xxz, tg_yzz_xxzz, tg_yzz_xyy, tg_yzz_xyyy, tg_yzz_xyyz, tg_yzz_xyz, tg_yzz_xyzz, tg_yzz_xzz, tg_yzz_xzzz, tg_yzz_yyy, tg_yzz_yyyy, tg_yzz_yyyz, tg_yzz_yyz, tg_yzz_yyzz, tg_yzz_yzz, tg_yzz_yzzz, tg_yzz_zzz, tg_yzz_zzzz, tg_yzzz_xxxx, tg_yzzz_xxxy, tg_yzzz_xxxz, tg_yzzz_xxyy, tg_yzzz_xxyz, tg_yzzz_xxzz, tg_yzzz_xyyy, tg_yzzz_xyyz, tg_yzzz_xyzz, tg_yzzz_xzzz, tg_yzzz_yyyy, tg_yzzz_yyyz, tg_yzzz_yyzz, tg_yzzz_yzzz, tg_yzzz_zzzz, tg_zz_xxxx, tg_zz_xxxy, tg_zz_xxxz, tg_zz_xxyy, tg_zz_xxyz, tg_zz_xxzz, tg_zz_xyyy, tg_zz_xyyz, tg_zz_xyzz, tg_zz_xzzz, tg_zz_yyyy, tg_zz_yyyz, tg_zz_yyzz, tg_zz_yzzz, tg_zz_zzzz, tg_zzz_xxx, tg_zzz_xxxx, tg_zzz_xxxy, tg_zzz_xxxz, tg_zzz_xxy, tg_zzz_xxyy, tg_zzz_xxyz, tg_zzz_xxz, tg_zzz_xxzz, tg_zzz_xyy, tg_zzz_xyyy, tg_zzz_xyyz, tg_zzz_xyz, tg_zzz_xyzz, tg_zzz_xzz, tg_zzz_xzzz, tg_zzz_yyy, tg_zzz_yyyy, tg_zzz_yyyz, tg_zzz_yyz, tg_zzz_yyzz, tg_zzz_yzz, tg_zzz_yzzz, tg_zzz_zzz, tg_zzz_zzzz, tg_zzzz_xxxx, tg_zzzz_xxxy, tg_zzzz_xxxz, tg_zzzz_xxyy, tg_zzzz_xxyz, tg_zzzz_xxzz, tg_zzzz_xyyy, tg_zzzz_xyyz, tg_zzzz_xyzz, tg_zzzz_xzzz, tg_zzzz_yyyy, tg_zzzz_yyyz, tg_zzzz_yyzz, tg_zzzz_yzzz, tg_zzzz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxx_xxxx[i] = 3.0 * tg_xx_xxxx[i] * fxi[i] + 4.0 * tg_xxx_xxx[i] * fxi[i] + tg_xxx_xxxx[i] * ra_x[i];

        tg_xxxx_xxxy[i] = 3.0 * tg_xx_xxxy[i] * fxi[i] + 3.0 * tg_xxx_xxy[i] * fxi[i] + tg_xxx_xxxy[i] * ra_x[i];

        tg_xxxx_xxxz[i] = 3.0 * tg_xx_xxxz[i] * fxi[i] + 3.0 * tg_xxx_xxz[i] * fxi[i] + tg_xxx_xxxz[i] * ra_x[i];

        tg_xxxx_xxyy[i] = 3.0 * tg_xx_xxyy[i] * fxi[i] + 2.0 * tg_xxx_xyy[i] * fxi[i] + tg_xxx_xxyy[i] * ra_x[i];

        tg_xxxx_xxyz[i] = 3.0 * tg_xx_xxyz[i] * fxi[i] + 2.0 * tg_xxx_xyz[i] * fxi[i] + tg_xxx_xxyz[i] * ra_x[i];

        tg_xxxx_xxzz[i] = 3.0 * tg_xx_xxzz[i] * fxi[i] + 2.0 * tg_xxx_xzz[i] * fxi[i] + tg_xxx_xxzz[i] * ra_x[i];

        tg_xxxx_xyyy[i] = 3.0 * tg_xx_xyyy[i] * fxi[i] + tg_xxx_yyy[i] * fxi[i] + tg_xxx_xyyy[i] * ra_x[i];

        tg_xxxx_xyyz[i] = 3.0 * tg_xx_xyyz[i] * fxi[i] + tg_xxx_yyz[i] * fxi[i] + tg_xxx_xyyz[i] * ra_x[i];

        tg_xxxx_xyzz[i] = 3.0 * tg_xx_xyzz[i] * fxi[i] + tg_xxx_yzz[i] * fxi[i] + tg_xxx_xyzz[i] * ra_x[i];

        tg_xxxx_xzzz[i] = 3.0 * tg_xx_xzzz[i] * fxi[i] + tg_xxx_zzz[i] * fxi[i] + tg_xxx_xzzz[i] * ra_x[i];

        tg_xxxx_yyyy[i] = 3.0 * tg_xx_yyyy[i] * fxi[i] + tg_xxx_yyyy[i] * ra_x[i];

        tg_xxxx_yyyz[i] = 3.0 * tg_xx_yyyz[i] * fxi[i] + tg_xxx_yyyz[i] * ra_x[i];

        tg_xxxx_yyzz[i] = 3.0 * tg_xx_yyzz[i] * fxi[i] + tg_xxx_yyzz[i] * ra_x[i];

        tg_xxxx_yzzz[i] = 3.0 * tg_xx_yzzz[i] * fxi[i] + tg_xxx_yzzz[i] * ra_x[i];

        tg_xxxx_zzzz[i] = 3.0 * tg_xx_zzzz[i] * fxi[i] + tg_xxx_zzzz[i] * ra_x[i];

        tg_xxxy_xxxx[i] = tg_xxx_xxxx[i] * ra_y[i];

        tg_xxxy_xxxy[i] = tg_xxx_xxx[i] * fxi[i] + tg_xxx_xxxy[i] * ra_y[i];

        tg_xxxy_xxxz[i] = tg_xxx_xxxz[i] * ra_y[i];

        tg_xxxy_xxyy[i] = 2.0 * tg_xxx_xxy[i] * fxi[i] + tg_xxx_xxyy[i] * ra_y[i];

        tg_xxxy_xxyz[i] = tg_xxx_xxz[i] * fxi[i] + tg_xxx_xxyz[i] * ra_y[i];

        tg_xxxy_xxzz[i] = tg_xxx_xxzz[i] * ra_y[i];

        tg_xxxy_xyyy[i] = 3.0 * tg_xxx_xyy[i] * fxi[i] + tg_xxx_xyyy[i] * ra_y[i];

        tg_xxxy_xyyz[i] = 2.0 * tg_xxx_xyz[i] * fxi[i] + tg_xxx_xyyz[i] * ra_y[i];

        tg_xxxy_xyzz[i] = tg_xxx_xzz[i] * fxi[i] + tg_xxx_xyzz[i] * ra_y[i];

        tg_xxxy_xzzz[i] = tg_xxx_xzzz[i] * ra_y[i];

        tg_xxxy_yyyy[i] = 2.0 * tg_xy_yyyy[i] * fxi[i] + tg_xxy_yyyy[i] * ra_x[i];

        tg_xxxy_yyyz[i] = 2.0 * tg_xy_yyyz[i] * fxi[i] + tg_xxy_yyyz[i] * ra_x[i];

        tg_xxxy_yyzz[i] = 2.0 * tg_xy_yyzz[i] * fxi[i] + tg_xxy_yyzz[i] * ra_x[i];

        tg_xxxy_yzzz[i] = 2.0 * tg_xy_yzzz[i] * fxi[i] + tg_xxy_yzzz[i] * ra_x[i];

        tg_xxxy_zzzz[i] = tg_xxx_zzzz[i] * ra_y[i];

        tg_xxxz_xxxx[i] = tg_xxx_xxxx[i] * ra_z[i];

        tg_xxxz_xxxy[i] = tg_xxx_xxxy[i] * ra_z[i];

        tg_xxxz_xxxz[i] = tg_xxx_xxx[i] * fxi[i] + tg_xxx_xxxz[i] * ra_z[i];

        tg_xxxz_xxyy[i] = tg_xxx_xxyy[i] * ra_z[i];

        tg_xxxz_xxyz[i] = tg_xxx_xxy[i] * fxi[i] + tg_xxx_xxyz[i] * ra_z[i];

        tg_xxxz_xxzz[i] = 2.0 * tg_xxx_xxz[i] * fxi[i] + tg_xxx_xxzz[i] * ra_z[i];

        tg_xxxz_xyyy[i] = tg_xxx_xyyy[i] * ra_z[i];

        tg_xxxz_xyyz[i] = tg_xxx_xyy[i] * fxi[i] + tg_xxx_xyyz[i] * ra_z[i];

        tg_xxxz_xyzz[i] = 2.0 * tg_xxx_xyz[i] * fxi[i] + tg_xxx_xyzz[i] * ra_z[i];

        tg_xxxz_xzzz[i] = 3.0 * tg_xxx_xzz[i] * fxi[i] + tg_xxx_xzzz[i] * ra_z[i];

        tg_xxxz_yyyy[i] = tg_xxx_yyyy[i] * ra_z[i];

        tg_xxxz_yyyz[i] = 2.0 * tg_xz_yyyz[i] * fxi[i] + tg_xxz_yyyz[i] * ra_x[i];

        tg_xxxz_yyzz[i] = 2.0 * tg_xz_yyzz[i] * fxi[i] + tg_xxz_yyzz[i] * ra_x[i];

        tg_xxxz_yzzz[i] = 2.0 * tg_xz_yzzz[i] * fxi[i] + tg_xxz_yzzz[i] * ra_x[i];

        tg_xxxz_zzzz[i] = 2.0 * tg_xz_zzzz[i] * fxi[i] + tg_xxz_zzzz[i] * ra_x[i];

        tg_xxyy_xxxx[i] = tg_xx_xxxx[i] * fxi[i] + tg_xxy_xxxx[i] * ra_y[i];

        tg_xxyy_xxxy[i] = tg_yy_xxxy[i] * fxi[i] + 3.0 * tg_xyy_xxy[i] * fxi[i] + tg_xyy_xxxy[i] * ra_x[i];

        tg_xxyy_xxxz[i] = tg_xx_xxxz[i] * fxi[i] + tg_xxy_xxxz[i] * ra_y[i];

        tg_xxyy_xxyy[i] = tg_yy_xxyy[i] * fxi[i] + 2.0 * tg_xyy_xyy[i] * fxi[i] + tg_xyy_xxyy[i] * ra_x[i];

        tg_xxyy_xxyz[i] = tg_yy_xxyz[i] * fxi[i] + 2.0 * tg_xyy_xyz[i] * fxi[i] + tg_xyy_xxyz[i] * ra_x[i];

        tg_xxyy_xxzz[i] = tg_xx_xxzz[i] * fxi[i] + tg_xxy_xxzz[i] * ra_y[i];

        tg_xxyy_xyyy[i] = tg_yy_xyyy[i] * fxi[i] + tg_xyy_yyy[i] * fxi[i] + tg_xyy_xyyy[i] * ra_x[i];

        tg_xxyy_xyyz[i] = tg_yy_xyyz[i] * fxi[i] + tg_xyy_yyz[i] * fxi[i] + tg_xyy_xyyz[i] * ra_x[i];

        tg_xxyy_xyzz[i] = tg_yy_xyzz[i] * fxi[i] + tg_xyy_yzz[i] * fxi[i] + tg_xyy_xyzz[i] * ra_x[i];

        tg_xxyy_xzzz[i] = tg_xx_xzzz[i] * fxi[i] + tg_xxy_xzzz[i] * ra_y[i];

        tg_xxyy_yyyy[i] = tg_yy_yyyy[i] * fxi[i] + tg_xyy_yyyy[i] * ra_x[i];

        tg_xxyy_yyyz[i] = tg_yy_yyyz[i] * fxi[i] + tg_xyy_yyyz[i] * ra_x[i];

        tg_xxyy_yyzz[i] = tg_yy_yyzz[i] * fxi[i] + tg_xyy_yyzz[i] * ra_x[i];

        tg_xxyy_yzzz[i] = tg_yy_yzzz[i] * fxi[i] + tg_xyy_yzzz[i] * ra_x[i];

        tg_xxyy_zzzz[i] = tg_yy_zzzz[i] * fxi[i] + tg_xyy_zzzz[i] * ra_x[i];

        tg_xxyz_xxxx[i] = tg_xxz_xxxx[i] * ra_y[i];

        tg_xxyz_xxxy[i] = tg_xxy_xxxy[i] * ra_z[i];

        tg_xxyz_xxxz[i] = tg_xxz_xxxz[i] * ra_y[i];

        tg_xxyz_xxyy[i] = tg_xxy_xxyy[i] * ra_z[i];

        tg_xxyz_xxyz[i] = tg_xxz_xxz[i] * fxi[i] + tg_xxz_xxyz[i] * ra_y[i];

        tg_xxyz_xxzz[i] = tg_xxz_xxzz[i] * ra_y[i];

        tg_xxyz_xyyy[i] = tg_xxy_xyyy[i] * ra_z[i];

        tg_xxyz_xyyz[i] = 2.0 * tg_xxz_xyz[i] * fxi[i] + tg_xxz_xyyz[i] * ra_y[i];

        tg_xxyz_xyzz[i] = tg_xxz_xzz[i] * fxi[i] + tg_xxz_xyzz[i] * ra_y[i];

        tg_xxyz_xzzz[i] = tg_xxz_xzzz[i] * ra_y[i];

        tg_xxyz_yyyy[i] = tg_xxy_yyyy[i] * ra_z[i];

        tg_xxyz_yyyz[i] = tg_yz_yyyz[i] * fxi[i] + tg_xyz_yyyz[i] * ra_x[i];

        tg_xxyz_yyzz[i] = tg_yz_yyzz[i] * fxi[i] + tg_xyz_yyzz[i] * ra_x[i];

        tg_xxyz_yzzz[i] = tg_yz_yzzz[i] * fxi[i] + tg_xyz_yzzz[i] * ra_x[i];

        tg_xxyz_zzzz[i] = tg_xxz_zzzz[i] * ra_y[i];

        tg_xxzz_xxxx[i] = tg_xx_xxxx[i] * fxi[i] + tg_xxz_xxxx[i] * ra_z[i];

        tg_xxzz_xxxy[i] = tg_xx_xxxy[i] * fxi[i] + tg_xxz_xxxy[i] * ra_z[i];

        tg_xxzz_xxxz[i] = tg_zz_xxxz[i] * fxi[i] + 3.0 * tg_xzz_xxz[i] * fxi[i] + tg_xzz_xxxz[i] * ra_x[i];

        tg_xxzz_xxyy[i] = tg_xx_xxyy[i] * fxi[i] + tg_xxz_xxyy[i] * ra_z[i];

        tg_xxzz_xxyz[i] = tg_zz_xxyz[i] * fxi[i] + 2.0 * tg_xzz_xyz[i] * fxi[i] + tg_xzz_xxyz[i] * ra_x[i];

        tg_xxzz_xxzz[i] = tg_zz_xxzz[i] * fxi[i] + 2.0 * tg_xzz_xzz[i] * fxi[i] + tg_xzz_xxzz[i] * ra_x[i];

        tg_xxzz_xyyy[i] = tg_xx_xyyy[i] * fxi[i] + tg_xxz_xyyy[i] * ra_z[i];

        tg_xxzz_xyyz[i] = tg_zz_xyyz[i] * fxi[i] + tg_xzz_yyz[i] * fxi[i] + tg_xzz_xyyz[i] * ra_x[i];

        tg_xxzz_xyzz[i] = tg_zz_xyzz[i] * fxi[i] + tg_xzz_yzz[i] * fxi[i] + tg_xzz_xyzz[i] * ra_x[i];

        tg_xxzz_xzzz[i] = tg_zz_xzzz[i] * fxi[i] + tg_xzz_zzz[i] * fxi[i] + tg_xzz_xzzz[i] * ra_x[i];

        tg_xxzz_yyyy[i] = tg_zz_yyyy[i] * fxi[i] + tg_xzz_yyyy[i] * ra_x[i];

        tg_xxzz_yyyz[i] = tg_zz_yyyz[i] * fxi[i] + tg_xzz_yyyz[i] * ra_x[i];

        tg_xxzz_yyzz[i] = tg_zz_yyzz[i] * fxi[i] + tg_xzz_yyzz[i] * ra_x[i];

        tg_xxzz_yzzz[i] = tg_zz_yzzz[i] * fxi[i] + tg_xzz_yzzz[i] * ra_x[i];

        tg_xxzz_zzzz[i] = tg_zz_zzzz[i] * fxi[i] + tg_xzz_zzzz[i] * ra_x[i];

        tg_xyyy_xxxx[i] = 4.0 * tg_yyy_xxx[i] * fxi[i] + tg_yyy_xxxx[i] * ra_x[i];

        tg_xyyy_xxxy[i] = 3.0 * tg_yyy_xxy[i] * fxi[i] + tg_yyy_xxxy[i] * ra_x[i];

        tg_xyyy_xxxz[i] = 3.0 * tg_yyy_xxz[i] * fxi[i] + tg_yyy_xxxz[i] * ra_x[i];

        tg_xyyy_xxyy[i] = 2.0 * tg_yyy_xyy[i] * fxi[i] + tg_yyy_xxyy[i] * ra_x[i];

        tg_xyyy_xxyz[i] = 2.0 * tg_yyy_xyz[i] * fxi[i] + tg_yyy_xxyz[i] * ra_x[i];

        tg_xyyy_xxzz[i] = 2.0 * tg_yyy_xzz[i] * fxi[i] + tg_yyy_xxzz[i] * ra_x[i];

        tg_xyyy_xyyy[i] = tg_yyy_yyy[i] * fxi[i] + tg_yyy_xyyy[i] * ra_x[i];

        tg_xyyy_xyyz[i] = tg_yyy_yyz[i] * fxi[i] + tg_yyy_xyyz[i] * ra_x[i];

        tg_xyyy_xyzz[i] = tg_yyy_yzz[i] * fxi[i] + tg_yyy_xyzz[i] * ra_x[i];

        tg_xyyy_xzzz[i] = tg_yyy_zzz[i] * fxi[i] + tg_yyy_xzzz[i] * ra_x[i];

        tg_xyyy_yyyy[i] = tg_yyy_yyyy[i] * ra_x[i];

        tg_xyyy_yyyz[i] = tg_yyy_yyyz[i] * ra_x[i];

        tg_xyyy_yyzz[i] = tg_yyy_yyzz[i] * ra_x[i];

        tg_xyyy_yzzz[i] = tg_yyy_yzzz[i] * ra_x[i];

        tg_xyyy_zzzz[i] = tg_yyy_zzzz[i] * ra_x[i];

        tg_xyyz_xxxx[i] = tg_xyy_xxxx[i] * ra_z[i];

        tg_xyyz_xxxy[i] = tg_xyy_xxxy[i] * ra_z[i];

        tg_xyyz_xxxz[i] = 3.0 * tg_yyz_xxz[i] * fxi[i] + tg_yyz_xxxz[i] * ra_x[i];

        tg_xyyz_xxyy[i] = tg_xyy_xxyy[i] * ra_z[i];

        tg_xyyz_xxyz[i] = 2.0 * tg_yyz_xyz[i] * fxi[i] + tg_yyz_xxyz[i] * ra_x[i];

        tg_xyyz_xxzz[i] = 2.0 * tg_yyz_xzz[i] * fxi[i] + tg_yyz_xxzz[i] * ra_x[i];

        tg_xyyz_xyyy[i] = tg_xyy_xyyy[i] * ra_z[i];

        tg_xyyz_xyyz[i] = tg_yyz_yyz[i] * fxi[i] + tg_yyz_xyyz[i] * ra_x[i];

        tg_xyyz_xyzz[i] = tg_yyz_yzz[i] * fxi[i] + tg_yyz_xyzz[i] * ra_x[i];

        tg_xyyz_xzzz[i] = tg_yyz_zzz[i] * fxi[i] + tg_yyz_xzzz[i] * ra_x[i];

        tg_xyyz_yyyy[i] = tg_yyz_yyyy[i] * ra_x[i];

        tg_xyyz_yyyz[i] = tg_yyz_yyyz[i] * ra_x[i];

        tg_xyyz_yyzz[i] = tg_yyz_yyzz[i] * ra_x[i];

        tg_xyyz_yzzz[i] = tg_yyz_yzzz[i] * ra_x[i];

        tg_xyyz_zzzz[i] = tg_yyz_zzzz[i] * ra_x[i];

        tg_xyzz_xxxx[i] = tg_xzz_xxxx[i] * ra_y[i];

        tg_xyzz_xxxy[i] = 3.0 * tg_yzz_xxy[i] * fxi[i] + tg_yzz_xxxy[i] * ra_x[i];

        tg_xyzz_xxxz[i] = tg_xzz_xxxz[i] * ra_y[i];

        tg_xyzz_xxyy[i] = 2.0 * tg_yzz_xyy[i] * fxi[i] + tg_yzz_xxyy[i] * ra_x[i];

        tg_xyzz_xxyz[i] = 2.0 * tg_yzz_xyz[i] * fxi[i] + tg_yzz_xxyz[i] * ra_x[i];

        tg_xyzz_xxzz[i] = tg_xzz_xxzz[i] * ra_y[i];

        tg_xyzz_xyyy[i] = tg_yzz_yyy[i] * fxi[i] + tg_yzz_xyyy[i] * ra_x[i];

        tg_xyzz_xyyz[i] = tg_yzz_yyz[i] * fxi[i] + tg_yzz_xyyz[i] * ra_x[i];

        tg_xyzz_xyzz[i] = tg_yzz_yzz[i] * fxi[i] + tg_yzz_xyzz[i] * ra_x[i];

        tg_xyzz_xzzz[i] = tg_xzz_xzzz[i] * ra_y[i];

        tg_xyzz_yyyy[i] = tg_yzz_yyyy[i] * ra_x[i];

        tg_xyzz_yyyz[i] = tg_yzz_yyyz[i] * ra_x[i];

        tg_xyzz_yyzz[i] = tg_yzz_yyzz[i] * ra_x[i];

        tg_xyzz_yzzz[i] = tg_yzz_yzzz[i] * ra_x[i];

        tg_xyzz_zzzz[i] = tg_yzz_zzzz[i] * ra_x[i];

        tg_xzzz_xxxx[i] = 4.0 * tg_zzz_xxx[i] * fxi[i] + tg_zzz_xxxx[i] * ra_x[i];

        tg_xzzz_xxxy[i] = 3.0 * tg_zzz_xxy[i] * fxi[i] + tg_zzz_xxxy[i] * ra_x[i];

        tg_xzzz_xxxz[i] = 3.0 * tg_zzz_xxz[i] * fxi[i] + tg_zzz_xxxz[i] * ra_x[i];

        tg_xzzz_xxyy[i] = 2.0 * tg_zzz_xyy[i] * fxi[i] + tg_zzz_xxyy[i] * ra_x[i];

        tg_xzzz_xxyz[i] = 2.0 * tg_zzz_xyz[i] * fxi[i] + tg_zzz_xxyz[i] * ra_x[i];

        tg_xzzz_xxzz[i] = 2.0 * tg_zzz_xzz[i] * fxi[i] + tg_zzz_xxzz[i] * ra_x[i];

        tg_xzzz_xyyy[i] = tg_zzz_yyy[i] * fxi[i] + tg_zzz_xyyy[i] * ra_x[i];

        tg_xzzz_xyyz[i] = tg_zzz_yyz[i] * fxi[i] + tg_zzz_xyyz[i] * ra_x[i];

        tg_xzzz_xyzz[i] = tg_zzz_yzz[i] * fxi[i] + tg_zzz_xyzz[i] * ra_x[i];

        tg_xzzz_xzzz[i] = tg_zzz_zzz[i] * fxi[i] + tg_zzz_xzzz[i] * ra_x[i];

        tg_xzzz_yyyy[i] = tg_zzz_yyyy[i] * ra_x[i];

        tg_xzzz_yyyz[i] = tg_zzz_yyyz[i] * ra_x[i];

        tg_xzzz_yyzz[i] = tg_zzz_yyzz[i] * ra_x[i];

        tg_xzzz_yzzz[i] = tg_zzz_yzzz[i] * ra_x[i];

        tg_xzzz_zzzz[i] = tg_zzz_zzzz[i] * ra_x[i];

        tg_yyyy_xxxx[i] = 3.0 * tg_yy_xxxx[i] * fxi[i] + tg_yyy_xxxx[i] * ra_y[i];

        tg_yyyy_xxxy[i] = 3.0 * tg_yy_xxxy[i] * fxi[i] + tg_yyy_xxx[i] * fxi[i] + tg_yyy_xxxy[i] * ra_y[i];

        tg_yyyy_xxxz[i] = 3.0 * tg_yy_xxxz[i] * fxi[i] + tg_yyy_xxxz[i] * ra_y[i];

        tg_yyyy_xxyy[i] = 3.0 * tg_yy_xxyy[i] * fxi[i] + 2.0 * tg_yyy_xxy[i] * fxi[i] + tg_yyy_xxyy[i] * ra_y[i];

        tg_yyyy_xxyz[i] = 3.0 * tg_yy_xxyz[i] * fxi[i] + tg_yyy_xxz[i] * fxi[i] + tg_yyy_xxyz[i] * ra_y[i];

        tg_yyyy_xxzz[i] = 3.0 * tg_yy_xxzz[i] * fxi[i] + tg_yyy_xxzz[i] * ra_y[i];

        tg_yyyy_xyyy[i] = 3.0 * tg_yy_xyyy[i] * fxi[i] + 3.0 * tg_yyy_xyy[i] * fxi[i] + tg_yyy_xyyy[i] * ra_y[i];

        tg_yyyy_xyyz[i] = 3.0 * tg_yy_xyyz[i] * fxi[i] + 2.0 * tg_yyy_xyz[i] * fxi[i] + tg_yyy_xyyz[i] * ra_y[i];

        tg_yyyy_xyzz[i] = 3.0 * tg_yy_xyzz[i] * fxi[i] + tg_yyy_xzz[i] * fxi[i] + tg_yyy_xyzz[i] * ra_y[i];

        tg_yyyy_xzzz[i] = 3.0 * tg_yy_xzzz[i] * fxi[i] + tg_yyy_xzzz[i] * ra_y[i];

        tg_yyyy_yyyy[i] = 3.0 * tg_yy_yyyy[i] * fxi[i] + 4.0 * tg_yyy_yyy[i] * fxi[i] + tg_yyy_yyyy[i] * ra_y[i];

        tg_yyyy_yyyz[i] = 3.0 * tg_yy_yyyz[i] * fxi[i] + 3.0 * tg_yyy_yyz[i] * fxi[i] + tg_yyy_yyyz[i] * ra_y[i];

        tg_yyyy_yyzz[i] = 3.0 * tg_yy_yyzz[i] * fxi[i] + 2.0 * tg_yyy_yzz[i] * fxi[i] + tg_yyy_yyzz[i] * ra_y[i];

        tg_yyyy_yzzz[i] = 3.0 * tg_yy_yzzz[i] * fxi[i] + tg_yyy_zzz[i] * fxi[i] + tg_yyy_yzzz[i] * ra_y[i];

        tg_yyyy_zzzz[i] = 3.0 * tg_yy_zzzz[i] * fxi[i] + tg_yyy_zzzz[i] * ra_y[i];

        tg_yyyz_xxxx[i] = tg_yyy_xxxx[i] * ra_z[i];

        tg_yyyz_xxxy[i] = tg_yyy_xxxy[i] * ra_z[i];

        tg_yyyz_xxxz[i] = 2.0 * tg_yz_xxxz[i] * fxi[i] + tg_yyz_xxxz[i] * ra_y[i];

        tg_yyyz_xxyy[i] = tg_yyy_xxyy[i] * ra_z[i];

        tg_yyyz_xxyz[i] = tg_yyy_xxy[i] * fxi[i] + tg_yyy_xxyz[i] * ra_z[i];

        tg_yyyz_xxzz[i] = 2.0 * tg_yz_xxzz[i] * fxi[i] + tg_yyz_xxzz[i] * ra_y[i];

        tg_yyyz_xyyy[i] = tg_yyy_xyyy[i] * ra_z[i];

        tg_yyyz_xyyz[i] = tg_yyy_xyy[i] * fxi[i] + tg_yyy_xyyz[i] * ra_z[i];

        tg_yyyz_xyzz[i] = 2.0 * tg_yyy_xyz[i] * fxi[i] + tg_yyy_xyzz[i] * ra_z[i];

        tg_yyyz_xzzz[i] = 2.0 * tg_yz_xzzz[i] * fxi[i] + tg_yyz_xzzz[i] * ra_y[i];

        tg_yyyz_yyyy[i] = tg_yyy_yyyy[i] * ra_z[i];

        tg_yyyz_yyyz[i] = tg_yyy_yyy[i] * fxi[i] + tg_yyy_yyyz[i] * ra_z[i];

        tg_yyyz_yyzz[i] = 2.0 * tg_yyy_yyz[i] * fxi[i] + tg_yyy_yyzz[i] * ra_z[i];

        tg_yyyz_yzzz[i] = 3.0 * tg_yyy_yzz[i] * fxi[i] + tg_yyy_yzzz[i] * ra_z[i];

        tg_yyyz_zzzz[i] = 2.0 * tg_yz_zzzz[i] * fxi[i] + tg_yyz_zzzz[i] * ra_y[i];

        tg_yyzz_xxxx[i] = tg_zz_xxxx[i] * fxi[i] + tg_yzz_xxxx[i] * ra_y[i];

        tg_yyzz_xxxy[i] = tg_yy_xxxy[i] * fxi[i] + tg_yyz_xxxy[i] * ra_z[i];

        tg_yyzz_xxxz[i] = tg_zz_xxxz[i] * fxi[i] + tg_yzz_xxxz[i] * ra_y[i];

        tg_yyzz_xxyy[i] = tg_yy_xxyy[i] * fxi[i] + tg_yyz_xxyy[i] * ra_z[i];

        tg_yyzz_xxyz[i] = tg_zz_xxyz[i] * fxi[i] + tg_yzz_xxz[i] * fxi[i] + tg_yzz_xxyz[i] * ra_y[i];

        tg_yyzz_xxzz[i] = tg_zz_xxzz[i] * fxi[i] + tg_yzz_xxzz[i] * ra_y[i];

        tg_yyzz_xyyy[i] = tg_yy_xyyy[i] * fxi[i] + tg_yyz_xyyy[i] * ra_z[i];

        tg_yyzz_xyyz[i] = tg_zz_xyyz[i] * fxi[i] + 2.0 * tg_yzz_xyz[i] * fxi[i] + tg_yzz_xyyz[i] * ra_y[i];

        tg_yyzz_xyzz[i] = tg_zz_xyzz[i] * fxi[i] + tg_yzz_xzz[i] * fxi[i] + tg_yzz_xyzz[i] * ra_y[i];

        tg_yyzz_xzzz[i] = tg_zz_xzzz[i] * fxi[i] + tg_yzz_xzzz[i] * ra_y[i];

        tg_yyzz_yyyy[i] = tg_yy_yyyy[i] * fxi[i] + tg_yyz_yyyy[i] * ra_z[i];

        tg_yyzz_yyyz[i] = tg_zz_yyyz[i] * fxi[i] + 3.0 * tg_yzz_yyz[i] * fxi[i] + tg_yzz_yyyz[i] * ra_y[i];

        tg_yyzz_yyzz[i] = tg_zz_yyzz[i] * fxi[i] + 2.0 * tg_yzz_yzz[i] * fxi[i] + tg_yzz_yyzz[i] * ra_y[i];

        tg_yyzz_yzzz[i] = tg_zz_yzzz[i] * fxi[i] + tg_yzz_zzz[i] * fxi[i] + tg_yzz_yzzz[i] * ra_y[i];

        tg_yyzz_zzzz[i] = tg_zz_zzzz[i] * fxi[i] + tg_yzz_zzzz[i] * ra_y[i];

        tg_yzzz_xxxx[i] = tg_zzz_xxxx[i] * ra_y[i];

        tg_yzzz_xxxy[i] = tg_zzz_xxx[i] * fxi[i] + tg_zzz_xxxy[i] * ra_y[i];

        tg_yzzz_xxxz[i] = tg_zzz_xxxz[i] * ra_y[i];

        tg_yzzz_xxyy[i] = 2.0 * tg_zzz_xxy[i] * fxi[i] + tg_zzz_xxyy[i] * ra_y[i];

        tg_yzzz_xxyz[i] = tg_zzz_xxz[i] * fxi[i] + tg_zzz_xxyz[i] * ra_y[i];

        tg_yzzz_xxzz[i] = tg_zzz_xxzz[i] * ra_y[i];

        tg_yzzz_xyyy[i] = 3.0 * tg_zzz_xyy[i] * fxi[i] + tg_zzz_xyyy[i] * ra_y[i];

        tg_yzzz_xyyz[i] = 2.0 * tg_zzz_xyz[i] * fxi[i] + tg_zzz_xyyz[i] * ra_y[i];

        tg_yzzz_xyzz[i] = tg_zzz_xzz[i] * fxi[i] + tg_zzz_xyzz[i] * ra_y[i];

        tg_yzzz_xzzz[i] = tg_zzz_xzzz[i] * ra_y[i];

        tg_yzzz_yyyy[i] = 4.0 * tg_zzz_yyy[i] * fxi[i] + tg_zzz_yyyy[i] * ra_y[i];

        tg_yzzz_yyyz[i] = 3.0 * tg_zzz_yyz[i] * fxi[i] + tg_zzz_yyyz[i] * ra_y[i];

        tg_yzzz_yyzz[i] = 2.0 * tg_zzz_yzz[i] * fxi[i] + tg_zzz_yyzz[i] * ra_y[i];

        tg_yzzz_yzzz[i] = tg_zzz_zzz[i] * fxi[i] + tg_zzz_yzzz[i] * ra_y[i];

        tg_yzzz_zzzz[i] = tg_zzz_zzzz[i] * ra_y[i];

        tg_zzzz_xxxx[i] = 3.0 * tg_zz_xxxx[i] * fxi[i] + tg_zzz_xxxx[i] * ra_z[i];

        tg_zzzz_xxxy[i] = 3.0 * tg_zz_xxxy[i] * fxi[i] + tg_zzz_xxxy[i] * ra_z[i];

        tg_zzzz_xxxz[i] = 3.0 * tg_zz_xxxz[i] * fxi[i] + tg_zzz_xxx[i] * fxi[i] + tg_zzz_xxxz[i] * ra_z[i];

        tg_zzzz_xxyy[i] = 3.0 * tg_zz_xxyy[i] * fxi[i] + tg_zzz_xxyy[i] * ra_z[i];

        tg_zzzz_xxyz[i] = 3.0 * tg_zz_xxyz[i] * fxi[i] + tg_zzz_xxy[i] * fxi[i] + tg_zzz_xxyz[i] * ra_z[i];

        tg_zzzz_xxzz[i] = 3.0 * tg_zz_xxzz[i] * fxi[i] + 2.0 * tg_zzz_xxz[i] * fxi[i] + tg_zzz_xxzz[i] * ra_z[i];

        tg_zzzz_xyyy[i] = 3.0 * tg_zz_xyyy[i] * fxi[i] + tg_zzz_xyyy[i] * ra_z[i];

        tg_zzzz_xyyz[i] = 3.0 * tg_zz_xyyz[i] * fxi[i] + tg_zzz_xyy[i] * fxi[i] + tg_zzz_xyyz[i] * ra_z[i];

        tg_zzzz_xyzz[i] = 3.0 * tg_zz_xyzz[i] * fxi[i] + 2.0 * tg_zzz_xyz[i] * fxi[i] + tg_zzz_xyzz[i] * ra_z[i];

        tg_zzzz_xzzz[i] = 3.0 * tg_zz_xzzz[i] * fxi[i] + 3.0 * tg_zzz_xzz[i] * fxi[i] + tg_zzz_xzzz[i] * ra_z[i];

        tg_zzzz_yyyy[i] = 3.0 * tg_zz_yyyy[i] * fxi[i] + tg_zzz_yyyy[i] * ra_z[i];

        tg_zzzz_yyyz[i] = 3.0 * tg_zz_yyyz[i] * fxi[i] + tg_zzz_yyy[i] * fxi[i] + tg_zzz_yyyz[i] * ra_z[i];

        tg_zzzz_yyzz[i] = 3.0 * tg_zz_yyzz[i] * fxi[i] + 2.0 * tg_zzz_yyz[i] * fxi[i] + tg_zzz_yyzz[i] * ra_z[i];

        tg_zzzz_yzzz[i] = 3.0 * tg_zz_yzzz[i] * fxi[i] + 3.0 * tg_zzz_yzz[i] * fxi[i] + tg_zzz_yzzz[i] * ra_z[i];

        tg_zzzz_zzzz[i] = 3.0 * tg_zz_zzzz[i] * fxi[i] + 4.0 * tg_zzz_zzz[i] * fxi[i] + tg_zzz_zzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

