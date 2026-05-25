#include "LocalCorePotentialPrimRecGH.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_gh(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gh,
                                  const size_t idx_dh,
                                  const size_t idx_fg,
                                  const size_t idx_fh,
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

    // Set up components of auxiliary buffer : DH

    auto tg_xx_xxxxx = pbuffer.data(idx_dh);

    auto tg_xx_xxxxy = pbuffer.data(idx_dh + 1);

    auto tg_xx_xxxxz = pbuffer.data(idx_dh + 2);

    auto tg_xx_xxxyy = pbuffer.data(idx_dh + 3);

    auto tg_xx_xxxyz = pbuffer.data(idx_dh + 4);

    auto tg_xx_xxxzz = pbuffer.data(idx_dh + 5);

    auto tg_xx_xxyyy = pbuffer.data(idx_dh + 6);

    auto tg_xx_xxyyz = pbuffer.data(idx_dh + 7);

    auto tg_xx_xxyzz = pbuffer.data(idx_dh + 8);

    auto tg_xx_xxzzz = pbuffer.data(idx_dh + 9);

    auto tg_xx_xyyyy = pbuffer.data(idx_dh + 10);

    auto tg_xx_xyyyz = pbuffer.data(idx_dh + 11);

    auto tg_xx_xyyzz = pbuffer.data(idx_dh + 12);

    auto tg_xx_xyzzz = pbuffer.data(idx_dh + 13);

    auto tg_xx_xzzzz = pbuffer.data(idx_dh + 14);

    auto tg_xx_yyyyy = pbuffer.data(idx_dh + 15);

    auto tg_xx_yyyyz = pbuffer.data(idx_dh + 16);

    auto tg_xx_yyyzz = pbuffer.data(idx_dh + 17);

    auto tg_xx_yyzzz = pbuffer.data(idx_dh + 18);

    auto tg_xx_yzzzz = pbuffer.data(idx_dh + 19);

    auto tg_xx_zzzzz = pbuffer.data(idx_dh + 20);

    auto tg_xy_yyyyy = pbuffer.data(idx_dh + 36);

    auto tg_xy_yyyyz = pbuffer.data(idx_dh + 37);

    auto tg_xy_yyyzz = pbuffer.data(idx_dh + 38);

    auto tg_xy_yyzzz = pbuffer.data(idx_dh + 39);

    auto tg_xy_yzzzz = pbuffer.data(idx_dh + 40);

    auto tg_xz_yyyyz = pbuffer.data(idx_dh + 58);

    auto tg_xz_yyyzz = pbuffer.data(idx_dh + 59);

    auto tg_xz_yyzzz = pbuffer.data(idx_dh + 60);

    auto tg_xz_yzzzz = pbuffer.data(idx_dh + 61);

    auto tg_xz_zzzzz = pbuffer.data(idx_dh + 62);

    auto tg_yy_xxxxx = pbuffer.data(idx_dh + 63);

    auto tg_yy_xxxxy = pbuffer.data(idx_dh + 64);

    auto tg_yy_xxxxz = pbuffer.data(idx_dh + 65);

    auto tg_yy_xxxyy = pbuffer.data(idx_dh + 66);

    auto tg_yy_xxxyz = pbuffer.data(idx_dh + 67);

    auto tg_yy_xxxzz = pbuffer.data(idx_dh + 68);

    auto tg_yy_xxyyy = pbuffer.data(idx_dh + 69);

    auto tg_yy_xxyyz = pbuffer.data(idx_dh + 70);

    auto tg_yy_xxyzz = pbuffer.data(idx_dh + 71);

    auto tg_yy_xxzzz = pbuffer.data(idx_dh + 72);

    auto tg_yy_xyyyy = pbuffer.data(idx_dh + 73);

    auto tg_yy_xyyyz = pbuffer.data(idx_dh + 74);

    auto tg_yy_xyyzz = pbuffer.data(idx_dh + 75);

    auto tg_yy_xyzzz = pbuffer.data(idx_dh + 76);

    auto tg_yy_xzzzz = pbuffer.data(idx_dh + 77);

    auto tg_yy_yyyyy = pbuffer.data(idx_dh + 78);

    auto tg_yy_yyyyz = pbuffer.data(idx_dh + 79);

    auto tg_yy_yyyzz = pbuffer.data(idx_dh + 80);

    auto tg_yy_yyzzz = pbuffer.data(idx_dh + 81);

    auto tg_yy_yzzzz = pbuffer.data(idx_dh + 82);

    auto tg_yy_zzzzz = pbuffer.data(idx_dh + 83);

    auto tg_yz_xxxxz = pbuffer.data(idx_dh + 86);

    auto tg_yz_xxxzz = pbuffer.data(idx_dh + 89);

    auto tg_yz_xxzzz = pbuffer.data(idx_dh + 93);

    auto tg_yz_xzzzz = pbuffer.data(idx_dh + 98);

    auto tg_yz_yyyyz = pbuffer.data(idx_dh + 100);

    auto tg_yz_yyyzz = pbuffer.data(idx_dh + 101);

    auto tg_yz_yyzzz = pbuffer.data(idx_dh + 102);

    auto tg_yz_yzzzz = pbuffer.data(idx_dh + 103);

    auto tg_yz_zzzzz = pbuffer.data(idx_dh + 104);

    auto tg_zz_xxxxx = pbuffer.data(idx_dh + 105);

    auto tg_zz_xxxxy = pbuffer.data(idx_dh + 106);

    auto tg_zz_xxxxz = pbuffer.data(idx_dh + 107);

    auto tg_zz_xxxyy = pbuffer.data(idx_dh + 108);

    auto tg_zz_xxxyz = pbuffer.data(idx_dh + 109);

    auto tg_zz_xxxzz = pbuffer.data(idx_dh + 110);

    auto tg_zz_xxyyy = pbuffer.data(idx_dh + 111);

    auto tg_zz_xxyyz = pbuffer.data(idx_dh + 112);

    auto tg_zz_xxyzz = pbuffer.data(idx_dh + 113);

    auto tg_zz_xxzzz = pbuffer.data(idx_dh + 114);

    auto tg_zz_xyyyy = pbuffer.data(idx_dh + 115);

    auto tg_zz_xyyyz = pbuffer.data(idx_dh + 116);

    auto tg_zz_xyyzz = pbuffer.data(idx_dh + 117);

    auto tg_zz_xyzzz = pbuffer.data(idx_dh + 118);

    auto tg_zz_xzzzz = pbuffer.data(idx_dh + 119);

    auto tg_zz_yyyyy = pbuffer.data(idx_dh + 120);

    auto tg_zz_yyyyz = pbuffer.data(idx_dh + 121);

    auto tg_zz_yyyzz = pbuffer.data(idx_dh + 122);

    auto tg_zz_yyzzz = pbuffer.data(idx_dh + 123);

    auto tg_zz_yzzzz = pbuffer.data(idx_dh + 124);

    auto tg_zz_zzzzz = pbuffer.data(idx_dh + 125);

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

    auto tg_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto tg_xxz_xxyz = pbuffer.data(idx_fg + 34);

    auto tg_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto tg_xxz_xyyz = pbuffer.data(idx_fg + 37);

    auto tg_xxz_xyzz = pbuffer.data(idx_fg + 38);

    auto tg_xxz_xzzz = pbuffer.data(idx_fg + 39);

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

    auto tg_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto tg_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto tg_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto tg_xzz_xyyz = pbuffer.data(idx_fg + 82);

    auto tg_xzz_xyzz = pbuffer.data(idx_fg + 83);

    auto tg_xzz_xzzz = pbuffer.data(idx_fg + 84);

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

    auto tg_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto tg_yyz_xxyz = pbuffer.data(idx_fg + 109);

    auto tg_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto tg_yyz_xyyz = pbuffer.data(idx_fg + 112);

    auto tg_yyz_xyzz = pbuffer.data(idx_fg + 113);

    auto tg_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto tg_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto tg_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto tg_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto tg_yyz_zzzz = pbuffer.data(idx_fg + 119);

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

    // Set up components of auxiliary buffer : FH

    auto tg_xxx_xxxxx = pbuffer.data(idx_fh);

    auto tg_xxx_xxxxy = pbuffer.data(idx_fh + 1);

    auto tg_xxx_xxxxz = pbuffer.data(idx_fh + 2);

    auto tg_xxx_xxxyy = pbuffer.data(idx_fh + 3);

    auto tg_xxx_xxxyz = pbuffer.data(idx_fh + 4);

    auto tg_xxx_xxxzz = pbuffer.data(idx_fh + 5);

    auto tg_xxx_xxyyy = pbuffer.data(idx_fh + 6);

    auto tg_xxx_xxyyz = pbuffer.data(idx_fh + 7);

    auto tg_xxx_xxyzz = pbuffer.data(idx_fh + 8);

    auto tg_xxx_xxzzz = pbuffer.data(idx_fh + 9);

    auto tg_xxx_xyyyy = pbuffer.data(idx_fh + 10);

    auto tg_xxx_xyyyz = pbuffer.data(idx_fh + 11);

    auto tg_xxx_xyyzz = pbuffer.data(idx_fh + 12);

    auto tg_xxx_xyzzz = pbuffer.data(idx_fh + 13);

    auto tg_xxx_xzzzz = pbuffer.data(idx_fh + 14);

    auto tg_xxx_yyyyy = pbuffer.data(idx_fh + 15);

    auto tg_xxx_yyyyz = pbuffer.data(idx_fh + 16);

    auto tg_xxx_yyyzz = pbuffer.data(idx_fh + 17);

    auto tg_xxx_yyzzz = pbuffer.data(idx_fh + 18);

    auto tg_xxx_yzzzz = pbuffer.data(idx_fh + 19);

    auto tg_xxx_zzzzz = pbuffer.data(idx_fh + 20);

    auto tg_xxy_xxxxx = pbuffer.data(idx_fh + 21);

    auto tg_xxy_xxxxy = pbuffer.data(idx_fh + 22);

    auto tg_xxy_xxxxz = pbuffer.data(idx_fh + 23);

    auto tg_xxy_xxxyy = pbuffer.data(idx_fh + 24);

    auto tg_xxy_xxxzz = pbuffer.data(idx_fh + 26);

    auto tg_xxy_xxyyy = pbuffer.data(idx_fh + 27);

    auto tg_xxy_xxzzz = pbuffer.data(idx_fh + 30);

    auto tg_xxy_xyyyy = pbuffer.data(idx_fh + 31);

    auto tg_xxy_xzzzz = pbuffer.data(idx_fh + 35);

    auto tg_xxy_yyyyy = pbuffer.data(idx_fh + 36);

    auto tg_xxy_yyyyz = pbuffer.data(idx_fh + 37);

    auto tg_xxy_yyyzz = pbuffer.data(idx_fh + 38);

    auto tg_xxy_yyzzz = pbuffer.data(idx_fh + 39);

    auto tg_xxy_yzzzz = pbuffer.data(idx_fh + 40);

    auto tg_xxz_xxxxx = pbuffer.data(idx_fh + 42);

    auto tg_xxz_xxxxy = pbuffer.data(idx_fh + 43);

    auto tg_xxz_xxxxz = pbuffer.data(idx_fh + 44);

    auto tg_xxz_xxxyy = pbuffer.data(idx_fh + 45);

    auto tg_xxz_xxxyz = pbuffer.data(idx_fh + 46);

    auto tg_xxz_xxxzz = pbuffer.data(idx_fh + 47);

    auto tg_xxz_xxyyy = pbuffer.data(idx_fh + 48);

    auto tg_xxz_xxyyz = pbuffer.data(idx_fh + 49);

    auto tg_xxz_xxyzz = pbuffer.data(idx_fh + 50);

    auto tg_xxz_xxzzz = pbuffer.data(idx_fh + 51);

    auto tg_xxz_xyyyy = pbuffer.data(idx_fh + 52);

    auto tg_xxz_xyyyz = pbuffer.data(idx_fh + 53);

    auto tg_xxz_xyyzz = pbuffer.data(idx_fh + 54);

    auto tg_xxz_xyzzz = pbuffer.data(idx_fh + 55);

    auto tg_xxz_xzzzz = pbuffer.data(idx_fh + 56);

    auto tg_xxz_yyyyz = pbuffer.data(idx_fh + 58);

    auto tg_xxz_yyyzz = pbuffer.data(idx_fh + 59);

    auto tg_xxz_yyzzz = pbuffer.data(idx_fh + 60);

    auto tg_xxz_yzzzz = pbuffer.data(idx_fh + 61);

    auto tg_xxz_zzzzz = pbuffer.data(idx_fh + 62);

    auto tg_xyy_xxxxx = pbuffer.data(idx_fh + 63);

    auto tg_xyy_xxxxy = pbuffer.data(idx_fh + 64);

    auto tg_xyy_xxxyy = pbuffer.data(idx_fh + 66);

    auto tg_xyy_xxxyz = pbuffer.data(idx_fh + 67);

    auto tg_xyy_xxyyy = pbuffer.data(idx_fh + 69);

    auto tg_xyy_xxyyz = pbuffer.data(idx_fh + 70);

    auto tg_xyy_xxyzz = pbuffer.data(idx_fh + 71);

    auto tg_xyy_xyyyy = pbuffer.data(idx_fh + 73);

    auto tg_xyy_xyyyz = pbuffer.data(idx_fh + 74);

    auto tg_xyy_xyyzz = pbuffer.data(idx_fh + 75);

    auto tg_xyy_xyzzz = pbuffer.data(idx_fh + 76);

    auto tg_xyy_yyyyy = pbuffer.data(idx_fh + 78);

    auto tg_xyy_yyyyz = pbuffer.data(idx_fh + 79);

    auto tg_xyy_yyyzz = pbuffer.data(idx_fh + 80);

    auto tg_xyy_yyzzz = pbuffer.data(idx_fh + 81);

    auto tg_xyy_yzzzz = pbuffer.data(idx_fh + 82);

    auto tg_xyy_zzzzz = pbuffer.data(idx_fh + 83);

    auto tg_xyz_yyyyz = pbuffer.data(idx_fh + 100);

    auto tg_xyz_yyyzz = pbuffer.data(idx_fh + 101);

    auto tg_xyz_yyzzz = pbuffer.data(idx_fh + 102);

    auto tg_xyz_yzzzz = pbuffer.data(idx_fh + 103);

    auto tg_xzz_xxxxx = pbuffer.data(idx_fh + 105);

    auto tg_xzz_xxxxz = pbuffer.data(idx_fh + 107);

    auto tg_xzz_xxxyz = pbuffer.data(idx_fh + 109);

    auto tg_xzz_xxxzz = pbuffer.data(idx_fh + 110);

    auto tg_xzz_xxyyz = pbuffer.data(idx_fh + 112);

    auto tg_xzz_xxyzz = pbuffer.data(idx_fh + 113);

    auto tg_xzz_xxzzz = pbuffer.data(idx_fh + 114);

    auto tg_xzz_xyyyz = pbuffer.data(idx_fh + 116);

    auto tg_xzz_xyyzz = pbuffer.data(idx_fh + 117);

    auto tg_xzz_xyzzz = pbuffer.data(idx_fh + 118);

    auto tg_xzz_xzzzz = pbuffer.data(idx_fh + 119);

    auto tg_xzz_yyyyy = pbuffer.data(idx_fh + 120);

    auto tg_xzz_yyyyz = pbuffer.data(idx_fh + 121);

    auto tg_xzz_yyyzz = pbuffer.data(idx_fh + 122);

    auto tg_xzz_yyzzz = pbuffer.data(idx_fh + 123);

    auto tg_xzz_yzzzz = pbuffer.data(idx_fh + 124);

    auto tg_xzz_zzzzz = pbuffer.data(idx_fh + 125);

    auto tg_yyy_xxxxx = pbuffer.data(idx_fh + 126);

    auto tg_yyy_xxxxy = pbuffer.data(idx_fh + 127);

    auto tg_yyy_xxxxz = pbuffer.data(idx_fh + 128);

    auto tg_yyy_xxxyy = pbuffer.data(idx_fh + 129);

    auto tg_yyy_xxxyz = pbuffer.data(idx_fh + 130);

    auto tg_yyy_xxxzz = pbuffer.data(idx_fh + 131);

    auto tg_yyy_xxyyy = pbuffer.data(idx_fh + 132);

    auto tg_yyy_xxyyz = pbuffer.data(idx_fh + 133);

    auto tg_yyy_xxyzz = pbuffer.data(idx_fh + 134);

    auto tg_yyy_xxzzz = pbuffer.data(idx_fh + 135);

    auto tg_yyy_xyyyy = pbuffer.data(idx_fh + 136);

    auto tg_yyy_xyyyz = pbuffer.data(idx_fh + 137);

    auto tg_yyy_xyyzz = pbuffer.data(idx_fh + 138);

    auto tg_yyy_xyzzz = pbuffer.data(idx_fh + 139);

    auto tg_yyy_xzzzz = pbuffer.data(idx_fh + 140);

    auto tg_yyy_yyyyy = pbuffer.data(idx_fh + 141);

    auto tg_yyy_yyyyz = pbuffer.data(idx_fh + 142);

    auto tg_yyy_yyyzz = pbuffer.data(idx_fh + 143);

    auto tg_yyy_yyzzz = pbuffer.data(idx_fh + 144);

    auto tg_yyy_yzzzz = pbuffer.data(idx_fh + 145);

    auto tg_yyy_zzzzz = pbuffer.data(idx_fh + 146);

    auto tg_yyz_xxxxy = pbuffer.data(idx_fh + 148);

    auto tg_yyz_xxxxz = pbuffer.data(idx_fh + 149);

    auto tg_yyz_xxxyy = pbuffer.data(idx_fh + 150);

    auto tg_yyz_xxxyz = pbuffer.data(idx_fh + 151);

    auto tg_yyz_xxxzz = pbuffer.data(idx_fh + 152);

    auto tg_yyz_xxyyy = pbuffer.data(idx_fh + 153);

    auto tg_yyz_xxyyz = pbuffer.data(idx_fh + 154);

    auto tg_yyz_xxyzz = pbuffer.data(idx_fh + 155);

    auto tg_yyz_xxzzz = pbuffer.data(idx_fh + 156);

    auto tg_yyz_xyyyy = pbuffer.data(idx_fh + 157);

    auto tg_yyz_xyyyz = pbuffer.data(idx_fh + 158);

    auto tg_yyz_xyyzz = pbuffer.data(idx_fh + 159);

    auto tg_yyz_xyzzz = pbuffer.data(idx_fh + 160);

    auto tg_yyz_xzzzz = pbuffer.data(idx_fh + 161);

    auto tg_yyz_yyyyy = pbuffer.data(idx_fh + 162);

    auto tg_yyz_yyyyz = pbuffer.data(idx_fh + 163);

    auto tg_yyz_yyyzz = pbuffer.data(idx_fh + 164);

    auto tg_yyz_yyzzz = pbuffer.data(idx_fh + 165);

    auto tg_yyz_yzzzz = pbuffer.data(idx_fh + 166);

    auto tg_yyz_zzzzz = pbuffer.data(idx_fh + 167);

    auto tg_yzz_xxxxx = pbuffer.data(idx_fh + 168);

    auto tg_yzz_xxxxy = pbuffer.data(idx_fh + 169);

    auto tg_yzz_xxxxz = pbuffer.data(idx_fh + 170);

    auto tg_yzz_xxxyy = pbuffer.data(idx_fh + 171);

    auto tg_yzz_xxxyz = pbuffer.data(idx_fh + 172);

    auto tg_yzz_xxxzz = pbuffer.data(idx_fh + 173);

    auto tg_yzz_xxyyy = pbuffer.data(idx_fh + 174);

    auto tg_yzz_xxyyz = pbuffer.data(idx_fh + 175);

    auto tg_yzz_xxyzz = pbuffer.data(idx_fh + 176);

    auto tg_yzz_xxzzz = pbuffer.data(idx_fh + 177);

    auto tg_yzz_xyyyy = pbuffer.data(idx_fh + 178);

    auto tg_yzz_xyyyz = pbuffer.data(idx_fh + 179);

    auto tg_yzz_xyyzz = pbuffer.data(idx_fh + 180);

    auto tg_yzz_xyzzz = pbuffer.data(idx_fh + 181);

    auto tg_yzz_xzzzz = pbuffer.data(idx_fh + 182);

    auto tg_yzz_yyyyy = pbuffer.data(idx_fh + 183);

    auto tg_yzz_yyyyz = pbuffer.data(idx_fh + 184);

    auto tg_yzz_yyyzz = pbuffer.data(idx_fh + 185);

    auto tg_yzz_yyzzz = pbuffer.data(idx_fh + 186);

    auto tg_yzz_yzzzz = pbuffer.data(idx_fh + 187);

    auto tg_yzz_zzzzz = pbuffer.data(idx_fh + 188);

    auto tg_zzz_xxxxx = pbuffer.data(idx_fh + 189);

    auto tg_zzz_xxxxy = pbuffer.data(idx_fh + 190);

    auto tg_zzz_xxxxz = pbuffer.data(idx_fh + 191);

    auto tg_zzz_xxxyy = pbuffer.data(idx_fh + 192);

    auto tg_zzz_xxxyz = pbuffer.data(idx_fh + 193);

    auto tg_zzz_xxxzz = pbuffer.data(idx_fh + 194);

    auto tg_zzz_xxyyy = pbuffer.data(idx_fh + 195);

    auto tg_zzz_xxyyz = pbuffer.data(idx_fh + 196);

    auto tg_zzz_xxyzz = pbuffer.data(idx_fh + 197);

    auto tg_zzz_xxzzz = pbuffer.data(idx_fh + 198);

    auto tg_zzz_xyyyy = pbuffer.data(idx_fh + 199);

    auto tg_zzz_xyyyz = pbuffer.data(idx_fh + 200);

    auto tg_zzz_xyyzz = pbuffer.data(idx_fh + 201);

    auto tg_zzz_xyzzz = pbuffer.data(idx_fh + 202);

    auto tg_zzz_xzzzz = pbuffer.data(idx_fh + 203);

    auto tg_zzz_yyyyy = pbuffer.data(idx_fh + 204);

    auto tg_zzz_yyyyz = pbuffer.data(idx_fh + 205);

    auto tg_zzz_yyyzz = pbuffer.data(idx_fh + 206);

    auto tg_zzz_yyzzz = pbuffer.data(idx_fh + 207);

    auto tg_zzz_yzzzz = pbuffer.data(idx_fh + 208);

    auto tg_zzz_zzzzz = pbuffer.data(idx_fh + 209);

    // Set up components of targeted buffer : GH

    auto tg_xxxx_xxxxx = pbuffer.data(idx_gh);

    auto tg_xxxx_xxxxy = pbuffer.data(idx_gh + 1);

    auto tg_xxxx_xxxxz = pbuffer.data(idx_gh + 2);

    auto tg_xxxx_xxxyy = pbuffer.data(idx_gh + 3);

    auto tg_xxxx_xxxyz = pbuffer.data(idx_gh + 4);

    auto tg_xxxx_xxxzz = pbuffer.data(idx_gh + 5);

    auto tg_xxxx_xxyyy = pbuffer.data(idx_gh + 6);

    auto tg_xxxx_xxyyz = pbuffer.data(idx_gh + 7);

    auto tg_xxxx_xxyzz = pbuffer.data(idx_gh + 8);

    auto tg_xxxx_xxzzz = pbuffer.data(idx_gh + 9);

    auto tg_xxxx_xyyyy = pbuffer.data(idx_gh + 10);

    auto tg_xxxx_xyyyz = pbuffer.data(idx_gh + 11);

    auto tg_xxxx_xyyzz = pbuffer.data(idx_gh + 12);

    auto tg_xxxx_xyzzz = pbuffer.data(idx_gh + 13);

    auto tg_xxxx_xzzzz = pbuffer.data(idx_gh + 14);

    auto tg_xxxx_yyyyy = pbuffer.data(idx_gh + 15);

    auto tg_xxxx_yyyyz = pbuffer.data(idx_gh + 16);

    auto tg_xxxx_yyyzz = pbuffer.data(idx_gh + 17);

    auto tg_xxxx_yyzzz = pbuffer.data(idx_gh + 18);

    auto tg_xxxx_yzzzz = pbuffer.data(idx_gh + 19);

    auto tg_xxxx_zzzzz = pbuffer.data(idx_gh + 20);

    auto tg_xxxy_xxxxx = pbuffer.data(idx_gh + 21);

    auto tg_xxxy_xxxxy = pbuffer.data(idx_gh + 22);

    auto tg_xxxy_xxxxz = pbuffer.data(idx_gh + 23);

    auto tg_xxxy_xxxyy = pbuffer.data(idx_gh + 24);

    auto tg_xxxy_xxxyz = pbuffer.data(idx_gh + 25);

    auto tg_xxxy_xxxzz = pbuffer.data(idx_gh + 26);

    auto tg_xxxy_xxyyy = pbuffer.data(idx_gh + 27);

    auto tg_xxxy_xxyyz = pbuffer.data(idx_gh + 28);

    auto tg_xxxy_xxyzz = pbuffer.data(idx_gh + 29);

    auto tg_xxxy_xxzzz = pbuffer.data(idx_gh + 30);

    auto tg_xxxy_xyyyy = pbuffer.data(idx_gh + 31);

    auto tg_xxxy_xyyyz = pbuffer.data(idx_gh + 32);

    auto tg_xxxy_xyyzz = pbuffer.data(idx_gh + 33);

    auto tg_xxxy_xyzzz = pbuffer.data(idx_gh + 34);

    auto tg_xxxy_xzzzz = pbuffer.data(idx_gh + 35);

    auto tg_xxxy_yyyyy = pbuffer.data(idx_gh + 36);

    auto tg_xxxy_yyyyz = pbuffer.data(idx_gh + 37);

    auto tg_xxxy_yyyzz = pbuffer.data(idx_gh + 38);

    auto tg_xxxy_yyzzz = pbuffer.data(idx_gh + 39);

    auto tg_xxxy_yzzzz = pbuffer.data(idx_gh + 40);

    auto tg_xxxy_zzzzz = pbuffer.data(idx_gh + 41);

    auto tg_xxxz_xxxxx = pbuffer.data(idx_gh + 42);

    auto tg_xxxz_xxxxy = pbuffer.data(idx_gh + 43);

    auto tg_xxxz_xxxxz = pbuffer.data(idx_gh + 44);

    auto tg_xxxz_xxxyy = pbuffer.data(idx_gh + 45);

    auto tg_xxxz_xxxyz = pbuffer.data(idx_gh + 46);

    auto tg_xxxz_xxxzz = pbuffer.data(idx_gh + 47);

    auto tg_xxxz_xxyyy = pbuffer.data(idx_gh + 48);

    auto tg_xxxz_xxyyz = pbuffer.data(idx_gh + 49);

    auto tg_xxxz_xxyzz = pbuffer.data(idx_gh + 50);

    auto tg_xxxz_xxzzz = pbuffer.data(idx_gh + 51);

    auto tg_xxxz_xyyyy = pbuffer.data(idx_gh + 52);

    auto tg_xxxz_xyyyz = pbuffer.data(idx_gh + 53);

    auto tg_xxxz_xyyzz = pbuffer.data(idx_gh + 54);

    auto tg_xxxz_xyzzz = pbuffer.data(idx_gh + 55);

    auto tg_xxxz_xzzzz = pbuffer.data(idx_gh + 56);

    auto tg_xxxz_yyyyy = pbuffer.data(idx_gh + 57);

    auto tg_xxxz_yyyyz = pbuffer.data(idx_gh + 58);

    auto tg_xxxz_yyyzz = pbuffer.data(idx_gh + 59);

    auto tg_xxxz_yyzzz = pbuffer.data(idx_gh + 60);

    auto tg_xxxz_yzzzz = pbuffer.data(idx_gh + 61);

    auto tg_xxxz_zzzzz = pbuffer.data(idx_gh + 62);

    auto tg_xxyy_xxxxx = pbuffer.data(idx_gh + 63);

    auto tg_xxyy_xxxxy = pbuffer.data(idx_gh + 64);

    auto tg_xxyy_xxxxz = pbuffer.data(idx_gh + 65);

    auto tg_xxyy_xxxyy = pbuffer.data(idx_gh + 66);

    auto tg_xxyy_xxxyz = pbuffer.data(idx_gh + 67);

    auto tg_xxyy_xxxzz = pbuffer.data(idx_gh + 68);

    auto tg_xxyy_xxyyy = pbuffer.data(idx_gh + 69);

    auto tg_xxyy_xxyyz = pbuffer.data(idx_gh + 70);

    auto tg_xxyy_xxyzz = pbuffer.data(idx_gh + 71);

    auto tg_xxyy_xxzzz = pbuffer.data(idx_gh + 72);

    auto tg_xxyy_xyyyy = pbuffer.data(idx_gh + 73);

    auto tg_xxyy_xyyyz = pbuffer.data(idx_gh + 74);

    auto tg_xxyy_xyyzz = pbuffer.data(idx_gh + 75);

    auto tg_xxyy_xyzzz = pbuffer.data(idx_gh + 76);

    auto tg_xxyy_xzzzz = pbuffer.data(idx_gh + 77);

    auto tg_xxyy_yyyyy = pbuffer.data(idx_gh + 78);

    auto tg_xxyy_yyyyz = pbuffer.data(idx_gh + 79);

    auto tg_xxyy_yyyzz = pbuffer.data(idx_gh + 80);

    auto tg_xxyy_yyzzz = pbuffer.data(idx_gh + 81);

    auto tg_xxyy_yzzzz = pbuffer.data(idx_gh + 82);

    auto tg_xxyy_zzzzz = pbuffer.data(idx_gh + 83);

    auto tg_xxyz_xxxxx = pbuffer.data(idx_gh + 84);

    auto tg_xxyz_xxxxy = pbuffer.data(idx_gh + 85);

    auto tg_xxyz_xxxxz = pbuffer.data(idx_gh + 86);

    auto tg_xxyz_xxxyy = pbuffer.data(idx_gh + 87);

    auto tg_xxyz_xxxyz = pbuffer.data(idx_gh + 88);

    auto tg_xxyz_xxxzz = pbuffer.data(idx_gh + 89);

    auto tg_xxyz_xxyyy = pbuffer.data(idx_gh + 90);

    auto tg_xxyz_xxyyz = pbuffer.data(idx_gh + 91);

    auto tg_xxyz_xxyzz = pbuffer.data(idx_gh + 92);

    auto tg_xxyz_xxzzz = pbuffer.data(idx_gh + 93);

    auto tg_xxyz_xyyyy = pbuffer.data(idx_gh + 94);

    auto tg_xxyz_xyyyz = pbuffer.data(idx_gh + 95);

    auto tg_xxyz_xyyzz = pbuffer.data(idx_gh + 96);

    auto tg_xxyz_xyzzz = pbuffer.data(idx_gh + 97);

    auto tg_xxyz_xzzzz = pbuffer.data(idx_gh + 98);

    auto tg_xxyz_yyyyy = pbuffer.data(idx_gh + 99);

    auto tg_xxyz_yyyyz = pbuffer.data(idx_gh + 100);

    auto tg_xxyz_yyyzz = pbuffer.data(idx_gh + 101);

    auto tg_xxyz_yyzzz = pbuffer.data(idx_gh + 102);

    auto tg_xxyz_yzzzz = pbuffer.data(idx_gh + 103);

    auto tg_xxyz_zzzzz = pbuffer.data(idx_gh + 104);

    auto tg_xxzz_xxxxx = pbuffer.data(idx_gh + 105);

    auto tg_xxzz_xxxxy = pbuffer.data(idx_gh + 106);

    auto tg_xxzz_xxxxz = pbuffer.data(idx_gh + 107);

    auto tg_xxzz_xxxyy = pbuffer.data(idx_gh + 108);

    auto tg_xxzz_xxxyz = pbuffer.data(idx_gh + 109);

    auto tg_xxzz_xxxzz = pbuffer.data(idx_gh + 110);

    auto tg_xxzz_xxyyy = pbuffer.data(idx_gh + 111);

    auto tg_xxzz_xxyyz = pbuffer.data(idx_gh + 112);

    auto tg_xxzz_xxyzz = pbuffer.data(idx_gh + 113);

    auto tg_xxzz_xxzzz = pbuffer.data(idx_gh + 114);

    auto tg_xxzz_xyyyy = pbuffer.data(idx_gh + 115);

    auto tg_xxzz_xyyyz = pbuffer.data(idx_gh + 116);

    auto tg_xxzz_xyyzz = pbuffer.data(idx_gh + 117);

    auto tg_xxzz_xyzzz = pbuffer.data(idx_gh + 118);

    auto tg_xxzz_xzzzz = pbuffer.data(idx_gh + 119);

    auto tg_xxzz_yyyyy = pbuffer.data(idx_gh + 120);

    auto tg_xxzz_yyyyz = pbuffer.data(idx_gh + 121);

    auto tg_xxzz_yyyzz = pbuffer.data(idx_gh + 122);

    auto tg_xxzz_yyzzz = pbuffer.data(idx_gh + 123);

    auto tg_xxzz_yzzzz = pbuffer.data(idx_gh + 124);

    auto tg_xxzz_zzzzz = pbuffer.data(idx_gh + 125);

    auto tg_xyyy_xxxxx = pbuffer.data(idx_gh + 126);

    auto tg_xyyy_xxxxy = pbuffer.data(idx_gh + 127);

    auto tg_xyyy_xxxxz = pbuffer.data(idx_gh + 128);

    auto tg_xyyy_xxxyy = pbuffer.data(idx_gh + 129);

    auto tg_xyyy_xxxyz = pbuffer.data(idx_gh + 130);

    auto tg_xyyy_xxxzz = pbuffer.data(idx_gh + 131);

    auto tg_xyyy_xxyyy = pbuffer.data(idx_gh + 132);

    auto tg_xyyy_xxyyz = pbuffer.data(idx_gh + 133);

    auto tg_xyyy_xxyzz = pbuffer.data(idx_gh + 134);

    auto tg_xyyy_xxzzz = pbuffer.data(idx_gh + 135);

    auto tg_xyyy_xyyyy = pbuffer.data(idx_gh + 136);

    auto tg_xyyy_xyyyz = pbuffer.data(idx_gh + 137);

    auto tg_xyyy_xyyzz = pbuffer.data(idx_gh + 138);

    auto tg_xyyy_xyzzz = pbuffer.data(idx_gh + 139);

    auto tg_xyyy_xzzzz = pbuffer.data(idx_gh + 140);

    auto tg_xyyy_yyyyy = pbuffer.data(idx_gh + 141);

    auto tg_xyyy_yyyyz = pbuffer.data(idx_gh + 142);

    auto tg_xyyy_yyyzz = pbuffer.data(idx_gh + 143);

    auto tg_xyyy_yyzzz = pbuffer.data(idx_gh + 144);

    auto tg_xyyy_yzzzz = pbuffer.data(idx_gh + 145);

    auto tg_xyyy_zzzzz = pbuffer.data(idx_gh + 146);

    auto tg_xyyz_xxxxx = pbuffer.data(idx_gh + 147);

    auto tg_xyyz_xxxxy = pbuffer.data(idx_gh + 148);

    auto tg_xyyz_xxxxz = pbuffer.data(idx_gh + 149);

    auto tg_xyyz_xxxyy = pbuffer.data(idx_gh + 150);

    auto tg_xyyz_xxxyz = pbuffer.data(idx_gh + 151);

    auto tg_xyyz_xxxzz = pbuffer.data(idx_gh + 152);

    auto tg_xyyz_xxyyy = pbuffer.data(idx_gh + 153);

    auto tg_xyyz_xxyyz = pbuffer.data(idx_gh + 154);

    auto tg_xyyz_xxyzz = pbuffer.data(idx_gh + 155);

    auto tg_xyyz_xxzzz = pbuffer.data(idx_gh + 156);

    auto tg_xyyz_xyyyy = pbuffer.data(idx_gh + 157);

    auto tg_xyyz_xyyyz = pbuffer.data(idx_gh + 158);

    auto tg_xyyz_xyyzz = pbuffer.data(idx_gh + 159);

    auto tg_xyyz_xyzzz = pbuffer.data(idx_gh + 160);

    auto tg_xyyz_xzzzz = pbuffer.data(idx_gh + 161);

    auto tg_xyyz_yyyyy = pbuffer.data(idx_gh + 162);

    auto tg_xyyz_yyyyz = pbuffer.data(idx_gh + 163);

    auto tg_xyyz_yyyzz = pbuffer.data(idx_gh + 164);

    auto tg_xyyz_yyzzz = pbuffer.data(idx_gh + 165);

    auto tg_xyyz_yzzzz = pbuffer.data(idx_gh + 166);

    auto tg_xyyz_zzzzz = pbuffer.data(idx_gh + 167);

    auto tg_xyzz_xxxxx = pbuffer.data(idx_gh + 168);

    auto tg_xyzz_xxxxy = pbuffer.data(idx_gh + 169);

    auto tg_xyzz_xxxxz = pbuffer.data(idx_gh + 170);

    auto tg_xyzz_xxxyy = pbuffer.data(idx_gh + 171);

    auto tg_xyzz_xxxyz = pbuffer.data(idx_gh + 172);

    auto tg_xyzz_xxxzz = pbuffer.data(idx_gh + 173);

    auto tg_xyzz_xxyyy = pbuffer.data(idx_gh + 174);

    auto tg_xyzz_xxyyz = pbuffer.data(idx_gh + 175);

    auto tg_xyzz_xxyzz = pbuffer.data(idx_gh + 176);

    auto tg_xyzz_xxzzz = pbuffer.data(idx_gh + 177);

    auto tg_xyzz_xyyyy = pbuffer.data(idx_gh + 178);

    auto tg_xyzz_xyyyz = pbuffer.data(idx_gh + 179);

    auto tg_xyzz_xyyzz = pbuffer.data(idx_gh + 180);

    auto tg_xyzz_xyzzz = pbuffer.data(idx_gh + 181);

    auto tg_xyzz_xzzzz = pbuffer.data(idx_gh + 182);

    auto tg_xyzz_yyyyy = pbuffer.data(idx_gh + 183);

    auto tg_xyzz_yyyyz = pbuffer.data(idx_gh + 184);

    auto tg_xyzz_yyyzz = pbuffer.data(idx_gh + 185);

    auto tg_xyzz_yyzzz = pbuffer.data(idx_gh + 186);

    auto tg_xyzz_yzzzz = pbuffer.data(idx_gh + 187);

    auto tg_xyzz_zzzzz = pbuffer.data(idx_gh + 188);

    auto tg_xzzz_xxxxx = pbuffer.data(idx_gh + 189);

    auto tg_xzzz_xxxxy = pbuffer.data(idx_gh + 190);

    auto tg_xzzz_xxxxz = pbuffer.data(idx_gh + 191);

    auto tg_xzzz_xxxyy = pbuffer.data(idx_gh + 192);

    auto tg_xzzz_xxxyz = pbuffer.data(idx_gh + 193);

    auto tg_xzzz_xxxzz = pbuffer.data(idx_gh + 194);

    auto tg_xzzz_xxyyy = pbuffer.data(idx_gh + 195);

    auto tg_xzzz_xxyyz = pbuffer.data(idx_gh + 196);

    auto tg_xzzz_xxyzz = pbuffer.data(idx_gh + 197);

    auto tg_xzzz_xxzzz = pbuffer.data(idx_gh + 198);

    auto tg_xzzz_xyyyy = pbuffer.data(idx_gh + 199);

    auto tg_xzzz_xyyyz = pbuffer.data(idx_gh + 200);

    auto tg_xzzz_xyyzz = pbuffer.data(idx_gh + 201);

    auto tg_xzzz_xyzzz = pbuffer.data(idx_gh + 202);

    auto tg_xzzz_xzzzz = pbuffer.data(idx_gh + 203);

    auto tg_xzzz_yyyyy = pbuffer.data(idx_gh + 204);

    auto tg_xzzz_yyyyz = pbuffer.data(idx_gh + 205);

    auto tg_xzzz_yyyzz = pbuffer.data(idx_gh + 206);

    auto tg_xzzz_yyzzz = pbuffer.data(idx_gh + 207);

    auto tg_xzzz_yzzzz = pbuffer.data(idx_gh + 208);

    auto tg_xzzz_zzzzz = pbuffer.data(idx_gh + 209);

    auto tg_yyyy_xxxxx = pbuffer.data(idx_gh + 210);

    auto tg_yyyy_xxxxy = pbuffer.data(idx_gh + 211);

    auto tg_yyyy_xxxxz = pbuffer.data(idx_gh + 212);

    auto tg_yyyy_xxxyy = pbuffer.data(idx_gh + 213);

    auto tg_yyyy_xxxyz = pbuffer.data(idx_gh + 214);

    auto tg_yyyy_xxxzz = pbuffer.data(idx_gh + 215);

    auto tg_yyyy_xxyyy = pbuffer.data(idx_gh + 216);

    auto tg_yyyy_xxyyz = pbuffer.data(idx_gh + 217);

    auto tg_yyyy_xxyzz = pbuffer.data(idx_gh + 218);

    auto tg_yyyy_xxzzz = pbuffer.data(idx_gh + 219);

    auto tg_yyyy_xyyyy = pbuffer.data(idx_gh + 220);

    auto tg_yyyy_xyyyz = pbuffer.data(idx_gh + 221);

    auto tg_yyyy_xyyzz = pbuffer.data(idx_gh + 222);

    auto tg_yyyy_xyzzz = pbuffer.data(idx_gh + 223);

    auto tg_yyyy_xzzzz = pbuffer.data(idx_gh + 224);

    auto tg_yyyy_yyyyy = pbuffer.data(idx_gh + 225);

    auto tg_yyyy_yyyyz = pbuffer.data(idx_gh + 226);

    auto tg_yyyy_yyyzz = pbuffer.data(idx_gh + 227);

    auto tg_yyyy_yyzzz = pbuffer.data(idx_gh + 228);

    auto tg_yyyy_yzzzz = pbuffer.data(idx_gh + 229);

    auto tg_yyyy_zzzzz = pbuffer.data(idx_gh + 230);

    auto tg_yyyz_xxxxx = pbuffer.data(idx_gh + 231);

    auto tg_yyyz_xxxxy = pbuffer.data(idx_gh + 232);

    auto tg_yyyz_xxxxz = pbuffer.data(idx_gh + 233);

    auto tg_yyyz_xxxyy = pbuffer.data(idx_gh + 234);

    auto tg_yyyz_xxxyz = pbuffer.data(idx_gh + 235);

    auto tg_yyyz_xxxzz = pbuffer.data(idx_gh + 236);

    auto tg_yyyz_xxyyy = pbuffer.data(idx_gh + 237);

    auto tg_yyyz_xxyyz = pbuffer.data(idx_gh + 238);

    auto tg_yyyz_xxyzz = pbuffer.data(idx_gh + 239);

    auto tg_yyyz_xxzzz = pbuffer.data(idx_gh + 240);

    auto tg_yyyz_xyyyy = pbuffer.data(idx_gh + 241);

    auto tg_yyyz_xyyyz = pbuffer.data(idx_gh + 242);

    auto tg_yyyz_xyyzz = pbuffer.data(idx_gh + 243);

    auto tg_yyyz_xyzzz = pbuffer.data(idx_gh + 244);

    auto tg_yyyz_xzzzz = pbuffer.data(idx_gh + 245);

    auto tg_yyyz_yyyyy = pbuffer.data(idx_gh + 246);

    auto tg_yyyz_yyyyz = pbuffer.data(idx_gh + 247);

    auto tg_yyyz_yyyzz = pbuffer.data(idx_gh + 248);

    auto tg_yyyz_yyzzz = pbuffer.data(idx_gh + 249);

    auto tg_yyyz_yzzzz = pbuffer.data(idx_gh + 250);

    auto tg_yyyz_zzzzz = pbuffer.data(idx_gh + 251);

    auto tg_yyzz_xxxxx = pbuffer.data(idx_gh + 252);

    auto tg_yyzz_xxxxy = pbuffer.data(idx_gh + 253);

    auto tg_yyzz_xxxxz = pbuffer.data(idx_gh + 254);

    auto tg_yyzz_xxxyy = pbuffer.data(idx_gh + 255);

    auto tg_yyzz_xxxyz = pbuffer.data(idx_gh + 256);

    auto tg_yyzz_xxxzz = pbuffer.data(idx_gh + 257);

    auto tg_yyzz_xxyyy = pbuffer.data(idx_gh + 258);

    auto tg_yyzz_xxyyz = pbuffer.data(idx_gh + 259);

    auto tg_yyzz_xxyzz = pbuffer.data(idx_gh + 260);

    auto tg_yyzz_xxzzz = pbuffer.data(idx_gh + 261);

    auto tg_yyzz_xyyyy = pbuffer.data(idx_gh + 262);

    auto tg_yyzz_xyyyz = pbuffer.data(idx_gh + 263);

    auto tg_yyzz_xyyzz = pbuffer.data(idx_gh + 264);

    auto tg_yyzz_xyzzz = pbuffer.data(idx_gh + 265);

    auto tg_yyzz_xzzzz = pbuffer.data(idx_gh + 266);

    auto tg_yyzz_yyyyy = pbuffer.data(idx_gh + 267);

    auto tg_yyzz_yyyyz = pbuffer.data(idx_gh + 268);

    auto tg_yyzz_yyyzz = pbuffer.data(idx_gh + 269);

    auto tg_yyzz_yyzzz = pbuffer.data(idx_gh + 270);

    auto tg_yyzz_yzzzz = pbuffer.data(idx_gh + 271);

    auto tg_yyzz_zzzzz = pbuffer.data(idx_gh + 272);

    auto tg_yzzz_xxxxx = pbuffer.data(idx_gh + 273);

    auto tg_yzzz_xxxxy = pbuffer.data(idx_gh + 274);

    auto tg_yzzz_xxxxz = pbuffer.data(idx_gh + 275);

    auto tg_yzzz_xxxyy = pbuffer.data(idx_gh + 276);

    auto tg_yzzz_xxxyz = pbuffer.data(idx_gh + 277);

    auto tg_yzzz_xxxzz = pbuffer.data(idx_gh + 278);

    auto tg_yzzz_xxyyy = pbuffer.data(idx_gh + 279);

    auto tg_yzzz_xxyyz = pbuffer.data(idx_gh + 280);

    auto tg_yzzz_xxyzz = pbuffer.data(idx_gh + 281);

    auto tg_yzzz_xxzzz = pbuffer.data(idx_gh + 282);

    auto tg_yzzz_xyyyy = pbuffer.data(idx_gh + 283);

    auto tg_yzzz_xyyyz = pbuffer.data(idx_gh + 284);

    auto tg_yzzz_xyyzz = pbuffer.data(idx_gh + 285);

    auto tg_yzzz_xyzzz = pbuffer.data(idx_gh + 286);

    auto tg_yzzz_xzzzz = pbuffer.data(idx_gh + 287);

    auto tg_yzzz_yyyyy = pbuffer.data(idx_gh + 288);

    auto tg_yzzz_yyyyz = pbuffer.data(idx_gh + 289);

    auto tg_yzzz_yyyzz = pbuffer.data(idx_gh + 290);

    auto tg_yzzz_yyzzz = pbuffer.data(idx_gh + 291);

    auto tg_yzzz_yzzzz = pbuffer.data(idx_gh + 292);

    auto tg_yzzz_zzzzz = pbuffer.data(idx_gh + 293);

    auto tg_zzzz_xxxxx = pbuffer.data(idx_gh + 294);

    auto tg_zzzz_xxxxy = pbuffer.data(idx_gh + 295);

    auto tg_zzzz_xxxxz = pbuffer.data(idx_gh + 296);

    auto tg_zzzz_xxxyy = pbuffer.data(idx_gh + 297);

    auto tg_zzzz_xxxyz = pbuffer.data(idx_gh + 298);

    auto tg_zzzz_xxxzz = pbuffer.data(idx_gh + 299);

    auto tg_zzzz_xxyyy = pbuffer.data(idx_gh + 300);

    auto tg_zzzz_xxyyz = pbuffer.data(idx_gh + 301);

    auto tg_zzzz_xxyzz = pbuffer.data(idx_gh + 302);

    auto tg_zzzz_xxzzz = pbuffer.data(idx_gh + 303);

    auto tg_zzzz_xyyyy = pbuffer.data(idx_gh + 304);

    auto tg_zzzz_xyyyz = pbuffer.data(idx_gh + 305);

    auto tg_zzzz_xyyzz = pbuffer.data(idx_gh + 306);

    auto tg_zzzz_xyzzz = pbuffer.data(idx_gh + 307);

    auto tg_zzzz_xzzzz = pbuffer.data(idx_gh + 308);

    auto tg_zzzz_yyyyy = pbuffer.data(idx_gh + 309);

    auto tg_zzzz_yyyyz = pbuffer.data(idx_gh + 310);

    auto tg_zzzz_yyyzz = pbuffer.data(idx_gh + 311);

    auto tg_zzzz_yyzzz = pbuffer.data(idx_gh + 312);

    auto tg_zzzz_yzzzz = pbuffer.data(idx_gh + 313);

    auto tg_zzzz_zzzzz = pbuffer.data(idx_gh + 314);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xx_xxxxx, tg_xx_xxxxy, tg_xx_xxxxz, tg_xx_xxxyy, tg_xx_xxxyz, tg_xx_xxxzz, tg_xx_xxyyy, tg_xx_xxyyz, tg_xx_xxyzz, tg_xx_xxzzz, tg_xx_xyyyy, tg_xx_xyyyz, tg_xx_xyyzz, tg_xx_xyzzz, tg_xx_xzzzz, tg_xx_yyyyy, tg_xx_yyyyz, tg_xx_yyyzz, tg_xx_yyzzz, tg_xx_yzzzz, tg_xx_zzzzz, tg_xxx_xxxx, tg_xxx_xxxxx, tg_xxx_xxxxy, tg_xxx_xxxxz, tg_xxx_xxxy, tg_xxx_xxxyy, tg_xxx_xxxyz, tg_xxx_xxxz, tg_xxx_xxxzz, tg_xxx_xxyy, tg_xxx_xxyyy, tg_xxx_xxyyz, tg_xxx_xxyz, tg_xxx_xxyzz, tg_xxx_xxzz, tg_xxx_xxzzz, tg_xxx_xyyy, tg_xxx_xyyyy, tg_xxx_xyyyz, tg_xxx_xyyz, tg_xxx_xyyzz, tg_xxx_xyzz, tg_xxx_xyzzz, tg_xxx_xzzz, tg_xxx_xzzzz, tg_xxx_yyyy, tg_xxx_yyyyy, tg_xxx_yyyyz, tg_xxx_yyyz, tg_xxx_yyyzz, tg_xxx_yyzz, tg_xxx_yyzzz, tg_xxx_yzzz, tg_xxx_yzzzz, tg_xxx_zzzz, tg_xxx_zzzzz, tg_xxxx_xxxxx, tg_xxxx_xxxxy, tg_xxxx_xxxxz, tg_xxxx_xxxyy, tg_xxxx_xxxyz, tg_xxxx_xxxzz, tg_xxxx_xxyyy, tg_xxxx_xxyyz, tg_xxxx_xxyzz, tg_xxxx_xxzzz, tg_xxxx_xyyyy, tg_xxxx_xyyyz, tg_xxxx_xyyzz, tg_xxxx_xyzzz, tg_xxxx_xzzzz, tg_xxxx_yyyyy, tg_xxxx_yyyyz, tg_xxxx_yyyzz, tg_xxxx_yyzzz, tg_xxxx_yzzzz, tg_xxxx_zzzzz, tg_xxxy_xxxxx, tg_xxxy_xxxxy, tg_xxxy_xxxxz, tg_xxxy_xxxyy, tg_xxxy_xxxyz, tg_xxxy_xxxzz, tg_xxxy_xxyyy, tg_xxxy_xxyyz, tg_xxxy_xxyzz, tg_xxxy_xxzzz, tg_xxxy_xyyyy, tg_xxxy_xyyyz, tg_xxxy_xyyzz, tg_xxxy_xyzzz, tg_xxxy_xzzzz, tg_xxxy_yyyyy, tg_xxxy_yyyyz, tg_xxxy_yyyzz, tg_xxxy_yyzzz, tg_xxxy_yzzzz, tg_xxxy_zzzzz, tg_xxxz_xxxxx, tg_xxxz_xxxxy, tg_xxxz_xxxxz, tg_xxxz_xxxyy, tg_xxxz_xxxyz, tg_xxxz_xxxzz, tg_xxxz_xxyyy, tg_xxxz_xxyyz, tg_xxxz_xxyzz, tg_xxxz_xxzzz, tg_xxxz_xyyyy, tg_xxxz_xyyyz, tg_xxxz_xyyzz, tg_xxxz_xyzzz, tg_xxxz_xzzzz, tg_xxxz_yyyyy, tg_xxxz_yyyyz, tg_xxxz_yyyzz, tg_xxxz_yyzzz, tg_xxxz_yzzzz, tg_xxxz_zzzzz, tg_xxy_xxxxx, tg_xxy_xxxxy, tg_xxy_xxxxz, tg_xxy_xxxyy, tg_xxy_xxxzz, tg_xxy_xxyyy, tg_xxy_xxzzz, tg_xxy_xyyyy, tg_xxy_xzzzz, tg_xxy_yyyyy, tg_xxy_yyyyz, tg_xxy_yyyzz, tg_xxy_yyzzz, tg_xxy_yzzzz, tg_xxyy_xxxxx, tg_xxyy_xxxxy, tg_xxyy_xxxxz, tg_xxyy_xxxyy, tg_xxyy_xxxyz, tg_xxyy_xxxzz, tg_xxyy_xxyyy, tg_xxyy_xxyyz, tg_xxyy_xxyzz, tg_xxyy_xxzzz, tg_xxyy_xyyyy, tg_xxyy_xyyyz, tg_xxyy_xyyzz, tg_xxyy_xyzzz, tg_xxyy_xzzzz, tg_xxyy_yyyyy, tg_xxyy_yyyyz, tg_xxyy_yyyzz, tg_xxyy_yyzzz, tg_xxyy_yzzzz, tg_xxyy_zzzzz, tg_xxyz_xxxxx, tg_xxyz_xxxxy, tg_xxyz_xxxxz, tg_xxyz_xxxyy, tg_xxyz_xxxyz, tg_xxyz_xxxzz, tg_xxyz_xxyyy, tg_xxyz_xxyyz, tg_xxyz_xxyzz, tg_xxyz_xxzzz, tg_xxyz_xyyyy, tg_xxyz_xyyyz, tg_xxyz_xyyzz, tg_xxyz_xyzzz, tg_xxyz_xzzzz, tg_xxyz_yyyyy, tg_xxyz_yyyyz, tg_xxyz_yyyzz, tg_xxyz_yyzzz, tg_xxyz_yzzzz, tg_xxyz_zzzzz, tg_xxz_xxxxx, tg_xxz_xxxxy, tg_xxz_xxxxz, tg_xxz_xxxyy, tg_xxz_xxxyz, tg_xxz_xxxz, tg_xxz_xxxzz, tg_xxz_xxyyy, tg_xxz_xxyyz, tg_xxz_xxyz, tg_xxz_xxyzz, tg_xxz_xxzz, tg_xxz_xxzzz, tg_xxz_xyyyy, tg_xxz_xyyyz, tg_xxz_xyyz, tg_xxz_xyyzz, tg_xxz_xyzz, tg_xxz_xyzzz, tg_xxz_xzzz, tg_xxz_xzzzz, tg_xxz_yyyyz, tg_xxz_yyyzz, tg_xxz_yyzzz, tg_xxz_yzzzz, tg_xxz_zzzzz, tg_xxzz_xxxxx, tg_xxzz_xxxxy, tg_xxzz_xxxxz, tg_xxzz_xxxyy, tg_xxzz_xxxyz, tg_xxzz_xxxzz, tg_xxzz_xxyyy, tg_xxzz_xxyyz, tg_xxzz_xxyzz, tg_xxzz_xxzzz, tg_xxzz_xyyyy, tg_xxzz_xyyyz, tg_xxzz_xyyzz, tg_xxzz_xyzzz, tg_xxzz_xzzzz, tg_xxzz_yyyyy, tg_xxzz_yyyyz, tg_xxzz_yyyzz, tg_xxzz_yyzzz, tg_xxzz_yzzzz, tg_xxzz_zzzzz, tg_xy_yyyyy, tg_xy_yyyyz, tg_xy_yyyzz, tg_xy_yyzzz, tg_xy_yzzzz, tg_xyy_xxxxx, tg_xyy_xxxxy, tg_xyy_xxxy, tg_xyy_xxxyy, tg_xyy_xxxyz, tg_xyy_xxyy, tg_xyy_xxyyy, tg_xyy_xxyyz, tg_xyy_xxyz, tg_xyy_xxyzz, tg_xyy_xyyy, tg_xyy_xyyyy, tg_xyy_xyyyz, tg_xyy_xyyz, tg_xyy_xyyzz, tg_xyy_xyzz, tg_xyy_xyzzz, tg_xyy_yyyy, tg_xyy_yyyyy, tg_xyy_yyyyz, tg_xyy_yyyz, tg_xyy_yyyzz, tg_xyy_yyzz, tg_xyy_yyzzz, tg_xyy_yzzz, tg_xyy_yzzzz, tg_xyy_zzzzz, tg_xyyy_xxxxx, tg_xyyy_xxxxy, tg_xyyy_xxxxz, tg_xyyy_xxxyy, tg_xyyy_xxxyz, tg_xyyy_xxxzz, tg_xyyy_xxyyy, tg_xyyy_xxyyz, tg_xyyy_xxyzz, tg_xyyy_xxzzz, tg_xyyy_xyyyy, tg_xyyy_xyyyz, tg_xyyy_xyyzz, tg_xyyy_xyzzz, tg_xyyy_xzzzz, tg_xyyy_yyyyy, tg_xyyy_yyyyz, tg_xyyy_yyyzz, tg_xyyy_yyzzz, tg_xyyy_yzzzz, tg_xyyy_zzzzz, tg_xyyz_xxxxx, tg_xyyz_xxxxy, tg_xyyz_xxxxz, tg_xyyz_xxxyy, tg_xyyz_xxxyz, tg_xyyz_xxxzz, tg_xyyz_xxyyy, tg_xyyz_xxyyz, tg_xyyz_xxyzz, tg_xyyz_xxzzz, tg_xyyz_xyyyy, tg_xyyz_xyyyz, tg_xyyz_xyyzz, tg_xyyz_xyzzz, tg_xyyz_xzzzz, tg_xyyz_yyyyy, tg_xyyz_yyyyz, tg_xyyz_yyyzz, tg_xyyz_yyzzz, tg_xyyz_yzzzz, tg_xyyz_zzzzz, tg_xyz_yyyyz, tg_xyz_yyyzz, tg_xyz_yyzzz, tg_xyz_yzzzz, tg_xyzz_xxxxx, tg_xyzz_xxxxy, tg_xyzz_xxxxz, tg_xyzz_xxxyy, tg_xyzz_xxxyz, tg_xyzz_xxxzz, tg_xyzz_xxyyy, tg_xyzz_xxyyz, tg_xyzz_xxyzz, tg_xyzz_xxzzz, tg_xyzz_xyyyy, tg_xyzz_xyyyz, tg_xyzz_xyyzz, tg_xyzz_xyzzz, tg_xyzz_xzzzz, tg_xyzz_yyyyy, tg_xyzz_yyyyz, tg_xyzz_yyyzz, tg_xyzz_yyzzz, tg_xyzz_yzzzz, tg_xyzz_zzzzz, tg_xz_yyyyz, tg_xz_yyyzz, tg_xz_yyzzz, tg_xz_yzzzz, tg_xz_zzzzz, tg_xzz_xxxxx, tg_xzz_xxxxz, tg_xzz_xxxyz, tg_xzz_xxxz, tg_xzz_xxxzz, tg_xzz_xxyyz, tg_xzz_xxyz, tg_xzz_xxyzz, tg_xzz_xxzz, tg_xzz_xxzzz, tg_xzz_xyyyz, tg_xzz_xyyz, tg_xzz_xyyzz, tg_xzz_xyzz, tg_xzz_xyzzz, tg_xzz_xzzz, tg_xzz_xzzzz, tg_xzz_yyyyy, tg_xzz_yyyyz, tg_xzz_yyyz, tg_xzz_yyyzz, tg_xzz_yyzz, tg_xzz_yyzzz, tg_xzz_yzzz, tg_xzz_yzzzz, tg_xzz_zzzz, tg_xzz_zzzzz, tg_xzzz_xxxxx, tg_xzzz_xxxxy, tg_xzzz_xxxxz, tg_xzzz_xxxyy, tg_xzzz_xxxyz, tg_xzzz_xxxzz, tg_xzzz_xxyyy, tg_xzzz_xxyyz, tg_xzzz_xxyzz, tg_xzzz_xxzzz, tg_xzzz_xyyyy, tg_xzzz_xyyyz, tg_xzzz_xyyzz, tg_xzzz_xyzzz, tg_xzzz_xzzzz, tg_xzzz_yyyyy, tg_xzzz_yyyyz, tg_xzzz_yyyzz, tg_xzzz_yyzzz, tg_xzzz_yzzzz, tg_xzzz_zzzzz, tg_yy_xxxxx, tg_yy_xxxxy, tg_yy_xxxxz, tg_yy_xxxyy, tg_yy_xxxyz, tg_yy_xxxzz, tg_yy_xxyyy, tg_yy_xxyyz, tg_yy_xxyzz, tg_yy_xxzzz, tg_yy_xyyyy, tg_yy_xyyyz, tg_yy_xyyzz, tg_yy_xyzzz, tg_yy_xzzzz, tg_yy_yyyyy, tg_yy_yyyyz, tg_yy_yyyzz, tg_yy_yyzzz, tg_yy_yzzzz, tg_yy_zzzzz, tg_yyy_xxxx, tg_yyy_xxxxx, tg_yyy_xxxxy, tg_yyy_xxxxz, tg_yyy_xxxy, tg_yyy_xxxyy, tg_yyy_xxxyz, tg_yyy_xxxz, tg_yyy_xxxzz, tg_yyy_xxyy, tg_yyy_xxyyy, tg_yyy_xxyyz, tg_yyy_xxyz, tg_yyy_xxyzz, tg_yyy_xxzz, tg_yyy_xxzzz, tg_yyy_xyyy, tg_yyy_xyyyy, tg_yyy_xyyyz, tg_yyy_xyyz, tg_yyy_xyyzz, tg_yyy_xyzz, tg_yyy_xyzzz, tg_yyy_xzzz, tg_yyy_xzzzz, tg_yyy_yyyy, tg_yyy_yyyyy, tg_yyy_yyyyz, tg_yyy_yyyz, tg_yyy_yyyzz, tg_yyy_yyzz, tg_yyy_yyzzz, tg_yyy_yzzz, tg_yyy_yzzzz, tg_yyy_zzzz, tg_yyy_zzzzz, tg_yyyy_xxxxx, tg_yyyy_xxxxy, tg_yyyy_xxxxz, tg_yyyy_xxxyy, tg_yyyy_xxxyz, tg_yyyy_xxxzz, tg_yyyy_xxyyy, tg_yyyy_xxyyz, tg_yyyy_xxyzz, tg_yyyy_xxzzz, tg_yyyy_xyyyy, tg_yyyy_xyyyz, tg_yyyy_xyyzz, tg_yyyy_xyzzz, tg_yyyy_xzzzz, tg_yyyy_yyyyy, tg_yyyy_yyyyz, tg_yyyy_yyyzz, tg_yyyy_yyzzz, tg_yyyy_yzzzz, tg_yyyy_zzzzz, tg_yyyz_xxxxx, tg_yyyz_xxxxy, tg_yyyz_xxxxz, tg_yyyz_xxxyy, tg_yyyz_xxxyz, tg_yyyz_xxxzz, tg_yyyz_xxyyy, tg_yyyz_xxyyz, tg_yyyz_xxyzz, tg_yyyz_xxzzz, tg_yyyz_xyyyy, tg_yyyz_xyyyz, tg_yyyz_xyyzz, tg_yyyz_xyzzz, tg_yyyz_xzzzz, tg_yyyz_yyyyy, tg_yyyz_yyyyz, tg_yyyz_yyyzz, tg_yyyz_yyzzz, tg_yyyz_yzzzz, tg_yyyz_zzzzz, tg_yyz_xxxxy, tg_yyz_xxxxz, tg_yyz_xxxyy, tg_yyz_xxxyz, tg_yyz_xxxz, tg_yyz_xxxzz, tg_yyz_xxyyy, tg_yyz_xxyyz, tg_yyz_xxyz, tg_yyz_xxyzz, tg_yyz_xxzz, tg_yyz_xxzzz, tg_yyz_xyyyy, tg_yyz_xyyyz, tg_yyz_xyyz, tg_yyz_xyyzz, tg_yyz_xyzz, tg_yyz_xyzzz, tg_yyz_xzzz, tg_yyz_xzzzz, tg_yyz_yyyyy, tg_yyz_yyyyz, tg_yyz_yyyz, tg_yyz_yyyzz, tg_yyz_yyzz, tg_yyz_yyzzz, tg_yyz_yzzz, tg_yyz_yzzzz, tg_yyz_zzzz, tg_yyz_zzzzz, tg_yyzz_xxxxx, tg_yyzz_xxxxy, tg_yyzz_xxxxz, tg_yyzz_xxxyy, tg_yyzz_xxxyz, tg_yyzz_xxxzz, tg_yyzz_xxyyy, tg_yyzz_xxyyz, tg_yyzz_xxyzz, tg_yyzz_xxzzz, tg_yyzz_xyyyy, tg_yyzz_xyyyz, tg_yyzz_xyyzz, tg_yyzz_xyzzz, tg_yyzz_xzzzz, tg_yyzz_yyyyy, tg_yyzz_yyyyz, tg_yyzz_yyyzz, tg_yyzz_yyzzz, tg_yyzz_yzzzz, tg_yyzz_zzzzz, tg_yz_xxxxz, tg_yz_xxxzz, tg_yz_xxzzz, tg_yz_xzzzz, tg_yz_yyyyz, tg_yz_yyyzz, tg_yz_yyzzz, tg_yz_yzzzz, tg_yz_zzzzz, tg_yzz_xxxxx, tg_yzz_xxxxy, tg_yzz_xxxxz, tg_yzz_xxxy, tg_yzz_xxxyy, tg_yzz_xxxyz, tg_yzz_xxxz, tg_yzz_xxxzz, tg_yzz_xxyy, tg_yzz_xxyyy, tg_yzz_xxyyz, tg_yzz_xxyz, tg_yzz_xxyzz, tg_yzz_xxzz, tg_yzz_xxzzz, tg_yzz_xyyy, tg_yzz_xyyyy, tg_yzz_xyyyz, tg_yzz_xyyz, tg_yzz_xyyzz, tg_yzz_xyzz, tg_yzz_xyzzz, tg_yzz_xzzz, tg_yzz_xzzzz, tg_yzz_yyyy, tg_yzz_yyyyy, tg_yzz_yyyyz, tg_yzz_yyyz, tg_yzz_yyyzz, tg_yzz_yyzz, tg_yzz_yyzzz, tg_yzz_yzzz, tg_yzz_yzzzz, tg_yzz_zzzz, tg_yzz_zzzzz, tg_yzzz_xxxxx, tg_yzzz_xxxxy, tg_yzzz_xxxxz, tg_yzzz_xxxyy, tg_yzzz_xxxyz, tg_yzzz_xxxzz, tg_yzzz_xxyyy, tg_yzzz_xxyyz, tg_yzzz_xxyzz, tg_yzzz_xxzzz, tg_yzzz_xyyyy, tg_yzzz_xyyyz, tg_yzzz_xyyzz, tg_yzzz_xyzzz, tg_yzzz_xzzzz, tg_yzzz_yyyyy, tg_yzzz_yyyyz, tg_yzzz_yyyzz, tg_yzzz_yyzzz, tg_yzzz_yzzzz, tg_yzzz_zzzzz, tg_zz_xxxxx, tg_zz_xxxxy, tg_zz_xxxxz, tg_zz_xxxyy, tg_zz_xxxyz, tg_zz_xxxzz, tg_zz_xxyyy, tg_zz_xxyyz, tg_zz_xxyzz, tg_zz_xxzzz, tg_zz_xyyyy, tg_zz_xyyyz, tg_zz_xyyzz, tg_zz_xyzzz, tg_zz_xzzzz, tg_zz_yyyyy, tg_zz_yyyyz, tg_zz_yyyzz, tg_zz_yyzzz, tg_zz_yzzzz, tg_zz_zzzzz, tg_zzz_xxxx, tg_zzz_xxxxx, tg_zzz_xxxxy, tg_zzz_xxxxz, tg_zzz_xxxy, tg_zzz_xxxyy, tg_zzz_xxxyz, tg_zzz_xxxz, tg_zzz_xxxzz, tg_zzz_xxyy, tg_zzz_xxyyy, tg_zzz_xxyyz, tg_zzz_xxyz, tg_zzz_xxyzz, tg_zzz_xxzz, tg_zzz_xxzzz, tg_zzz_xyyy, tg_zzz_xyyyy, tg_zzz_xyyyz, tg_zzz_xyyz, tg_zzz_xyyzz, tg_zzz_xyzz, tg_zzz_xyzzz, tg_zzz_xzzz, tg_zzz_xzzzz, tg_zzz_yyyy, tg_zzz_yyyyy, tg_zzz_yyyyz, tg_zzz_yyyz, tg_zzz_yyyzz, tg_zzz_yyzz, tg_zzz_yyzzz, tg_zzz_yzzz, tg_zzz_yzzzz, tg_zzz_zzzz, tg_zzz_zzzzz, tg_zzzz_xxxxx, tg_zzzz_xxxxy, tg_zzzz_xxxxz, tg_zzzz_xxxyy, tg_zzzz_xxxyz, tg_zzzz_xxxzz, tg_zzzz_xxyyy, tg_zzzz_xxyyz, tg_zzzz_xxyzz, tg_zzzz_xxzzz, tg_zzzz_xyyyy, tg_zzzz_xyyyz, tg_zzzz_xyyzz, tg_zzzz_xyzzz, tg_zzzz_xzzzz, tg_zzzz_yyyyy, tg_zzzz_yyyyz, tg_zzzz_yyyzz, tg_zzzz_yyzzz, tg_zzzz_yzzzz, tg_zzzz_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxx_xxxxx[i] = 3.0 * tg_xx_xxxxx[i] * fxi[i] + 5.0 * tg_xxx_xxxx[i] * fxi[i] + tg_xxx_xxxxx[i] * ra_x[i];

        tg_xxxx_xxxxy[i] = 3.0 * tg_xx_xxxxy[i] * fxi[i] + 4.0 * tg_xxx_xxxy[i] * fxi[i] + tg_xxx_xxxxy[i] * ra_x[i];

        tg_xxxx_xxxxz[i] = 3.0 * tg_xx_xxxxz[i] * fxi[i] + 4.0 * tg_xxx_xxxz[i] * fxi[i] + tg_xxx_xxxxz[i] * ra_x[i];

        tg_xxxx_xxxyy[i] = 3.0 * tg_xx_xxxyy[i] * fxi[i] + 3.0 * tg_xxx_xxyy[i] * fxi[i] + tg_xxx_xxxyy[i] * ra_x[i];

        tg_xxxx_xxxyz[i] = 3.0 * tg_xx_xxxyz[i] * fxi[i] + 3.0 * tg_xxx_xxyz[i] * fxi[i] + tg_xxx_xxxyz[i] * ra_x[i];

        tg_xxxx_xxxzz[i] = 3.0 * tg_xx_xxxzz[i] * fxi[i] + 3.0 * tg_xxx_xxzz[i] * fxi[i] + tg_xxx_xxxzz[i] * ra_x[i];

        tg_xxxx_xxyyy[i] = 3.0 * tg_xx_xxyyy[i] * fxi[i] + 2.0 * tg_xxx_xyyy[i] * fxi[i] + tg_xxx_xxyyy[i] * ra_x[i];

        tg_xxxx_xxyyz[i] = 3.0 * tg_xx_xxyyz[i] * fxi[i] + 2.0 * tg_xxx_xyyz[i] * fxi[i] + tg_xxx_xxyyz[i] * ra_x[i];

        tg_xxxx_xxyzz[i] = 3.0 * tg_xx_xxyzz[i] * fxi[i] + 2.0 * tg_xxx_xyzz[i] * fxi[i] + tg_xxx_xxyzz[i] * ra_x[i];

        tg_xxxx_xxzzz[i] = 3.0 * tg_xx_xxzzz[i] * fxi[i] + 2.0 * tg_xxx_xzzz[i] * fxi[i] + tg_xxx_xxzzz[i] * ra_x[i];

        tg_xxxx_xyyyy[i] = 3.0 * tg_xx_xyyyy[i] * fxi[i] + tg_xxx_yyyy[i] * fxi[i] + tg_xxx_xyyyy[i] * ra_x[i];

        tg_xxxx_xyyyz[i] = 3.0 * tg_xx_xyyyz[i] * fxi[i] + tg_xxx_yyyz[i] * fxi[i] + tg_xxx_xyyyz[i] * ra_x[i];

        tg_xxxx_xyyzz[i] = 3.0 * tg_xx_xyyzz[i] * fxi[i] + tg_xxx_yyzz[i] * fxi[i] + tg_xxx_xyyzz[i] * ra_x[i];

        tg_xxxx_xyzzz[i] = 3.0 * tg_xx_xyzzz[i] * fxi[i] + tg_xxx_yzzz[i] * fxi[i] + tg_xxx_xyzzz[i] * ra_x[i];

        tg_xxxx_xzzzz[i] = 3.0 * tg_xx_xzzzz[i] * fxi[i] + tg_xxx_zzzz[i] * fxi[i] + tg_xxx_xzzzz[i] * ra_x[i];

        tg_xxxx_yyyyy[i] = 3.0 * tg_xx_yyyyy[i] * fxi[i] + tg_xxx_yyyyy[i] * ra_x[i];

        tg_xxxx_yyyyz[i] = 3.0 * tg_xx_yyyyz[i] * fxi[i] + tg_xxx_yyyyz[i] * ra_x[i];

        tg_xxxx_yyyzz[i] = 3.0 * tg_xx_yyyzz[i] * fxi[i] + tg_xxx_yyyzz[i] * ra_x[i];

        tg_xxxx_yyzzz[i] = 3.0 * tg_xx_yyzzz[i] * fxi[i] + tg_xxx_yyzzz[i] * ra_x[i];

        tg_xxxx_yzzzz[i] = 3.0 * tg_xx_yzzzz[i] * fxi[i] + tg_xxx_yzzzz[i] * ra_x[i];

        tg_xxxx_zzzzz[i] = 3.0 * tg_xx_zzzzz[i] * fxi[i] + tg_xxx_zzzzz[i] * ra_x[i];

        tg_xxxy_xxxxx[i] = tg_xxx_xxxxx[i] * ra_y[i];

        tg_xxxy_xxxxy[i] = tg_xxx_xxxx[i] * fxi[i] + tg_xxx_xxxxy[i] * ra_y[i];

        tg_xxxy_xxxxz[i] = tg_xxx_xxxxz[i] * ra_y[i];

        tg_xxxy_xxxyy[i] = 2.0 * tg_xxx_xxxy[i] * fxi[i] + tg_xxx_xxxyy[i] * ra_y[i];

        tg_xxxy_xxxyz[i] = tg_xxx_xxxz[i] * fxi[i] + tg_xxx_xxxyz[i] * ra_y[i];

        tg_xxxy_xxxzz[i] = tg_xxx_xxxzz[i] * ra_y[i];

        tg_xxxy_xxyyy[i] = 3.0 * tg_xxx_xxyy[i] * fxi[i] + tg_xxx_xxyyy[i] * ra_y[i];

        tg_xxxy_xxyyz[i] = 2.0 * tg_xxx_xxyz[i] * fxi[i] + tg_xxx_xxyyz[i] * ra_y[i];

        tg_xxxy_xxyzz[i] = tg_xxx_xxzz[i] * fxi[i] + tg_xxx_xxyzz[i] * ra_y[i];

        tg_xxxy_xxzzz[i] = tg_xxx_xxzzz[i] * ra_y[i];

        tg_xxxy_xyyyy[i] = 4.0 * tg_xxx_xyyy[i] * fxi[i] + tg_xxx_xyyyy[i] * ra_y[i];

        tg_xxxy_xyyyz[i] = 3.0 * tg_xxx_xyyz[i] * fxi[i] + tg_xxx_xyyyz[i] * ra_y[i];

        tg_xxxy_xyyzz[i] = 2.0 * tg_xxx_xyzz[i] * fxi[i] + tg_xxx_xyyzz[i] * ra_y[i];

        tg_xxxy_xyzzz[i] = tg_xxx_xzzz[i] * fxi[i] + tg_xxx_xyzzz[i] * ra_y[i];

        tg_xxxy_xzzzz[i] = tg_xxx_xzzzz[i] * ra_y[i];

        tg_xxxy_yyyyy[i] = 2.0 * tg_xy_yyyyy[i] * fxi[i] + tg_xxy_yyyyy[i] * ra_x[i];

        tg_xxxy_yyyyz[i] = 2.0 * tg_xy_yyyyz[i] * fxi[i] + tg_xxy_yyyyz[i] * ra_x[i];

        tg_xxxy_yyyzz[i] = 2.0 * tg_xy_yyyzz[i] * fxi[i] + tg_xxy_yyyzz[i] * ra_x[i];

        tg_xxxy_yyzzz[i] = 2.0 * tg_xy_yyzzz[i] * fxi[i] + tg_xxy_yyzzz[i] * ra_x[i];

        tg_xxxy_yzzzz[i] = 2.0 * tg_xy_yzzzz[i] * fxi[i] + tg_xxy_yzzzz[i] * ra_x[i];

        tg_xxxy_zzzzz[i] = tg_xxx_zzzzz[i] * ra_y[i];

        tg_xxxz_xxxxx[i] = tg_xxx_xxxxx[i] * ra_z[i];

        tg_xxxz_xxxxy[i] = tg_xxx_xxxxy[i] * ra_z[i];

        tg_xxxz_xxxxz[i] = tg_xxx_xxxx[i] * fxi[i] + tg_xxx_xxxxz[i] * ra_z[i];

        tg_xxxz_xxxyy[i] = tg_xxx_xxxyy[i] * ra_z[i];

        tg_xxxz_xxxyz[i] = tg_xxx_xxxy[i] * fxi[i] + tg_xxx_xxxyz[i] * ra_z[i];

        tg_xxxz_xxxzz[i] = 2.0 * tg_xxx_xxxz[i] * fxi[i] + tg_xxx_xxxzz[i] * ra_z[i];

        tg_xxxz_xxyyy[i] = tg_xxx_xxyyy[i] * ra_z[i];

        tg_xxxz_xxyyz[i] = tg_xxx_xxyy[i] * fxi[i] + tg_xxx_xxyyz[i] * ra_z[i];

        tg_xxxz_xxyzz[i] = 2.0 * tg_xxx_xxyz[i] * fxi[i] + tg_xxx_xxyzz[i] * ra_z[i];

        tg_xxxz_xxzzz[i] = 3.0 * tg_xxx_xxzz[i] * fxi[i] + tg_xxx_xxzzz[i] * ra_z[i];

        tg_xxxz_xyyyy[i] = tg_xxx_xyyyy[i] * ra_z[i];

        tg_xxxz_xyyyz[i] = tg_xxx_xyyy[i] * fxi[i] + tg_xxx_xyyyz[i] * ra_z[i];

        tg_xxxz_xyyzz[i] = 2.0 * tg_xxx_xyyz[i] * fxi[i] + tg_xxx_xyyzz[i] * ra_z[i];

        tg_xxxz_xyzzz[i] = 3.0 * tg_xxx_xyzz[i] * fxi[i] + tg_xxx_xyzzz[i] * ra_z[i];

        tg_xxxz_xzzzz[i] = 4.0 * tg_xxx_xzzz[i] * fxi[i] + tg_xxx_xzzzz[i] * ra_z[i];

        tg_xxxz_yyyyy[i] = tg_xxx_yyyyy[i] * ra_z[i];

        tg_xxxz_yyyyz[i] = 2.0 * tg_xz_yyyyz[i] * fxi[i] + tg_xxz_yyyyz[i] * ra_x[i];

        tg_xxxz_yyyzz[i] = 2.0 * tg_xz_yyyzz[i] * fxi[i] + tg_xxz_yyyzz[i] * ra_x[i];

        tg_xxxz_yyzzz[i] = 2.0 * tg_xz_yyzzz[i] * fxi[i] + tg_xxz_yyzzz[i] * ra_x[i];

        tg_xxxz_yzzzz[i] = 2.0 * tg_xz_yzzzz[i] * fxi[i] + tg_xxz_yzzzz[i] * ra_x[i];

        tg_xxxz_zzzzz[i] = 2.0 * tg_xz_zzzzz[i] * fxi[i] + tg_xxz_zzzzz[i] * ra_x[i];

        tg_xxyy_xxxxx[i] = tg_xx_xxxxx[i] * fxi[i] + tg_xxy_xxxxx[i] * ra_y[i];

        tg_xxyy_xxxxy[i] = tg_yy_xxxxy[i] * fxi[i] + 4.0 * tg_xyy_xxxy[i] * fxi[i] + tg_xyy_xxxxy[i] * ra_x[i];

        tg_xxyy_xxxxz[i] = tg_xx_xxxxz[i] * fxi[i] + tg_xxy_xxxxz[i] * ra_y[i];

        tg_xxyy_xxxyy[i] = tg_yy_xxxyy[i] * fxi[i] + 3.0 * tg_xyy_xxyy[i] * fxi[i] + tg_xyy_xxxyy[i] * ra_x[i];

        tg_xxyy_xxxyz[i] = tg_yy_xxxyz[i] * fxi[i] + 3.0 * tg_xyy_xxyz[i] * fxi[i] + tg_xyy_xxxyz[i] * ra_x[i];

        tg_xxyy_xxxzz[i] = tg_xx_xxxzz[i] * fxi[i] + tg_xxy_xxxzz[i] * ra_y[i];

        tg_xxyy_xxyyy[i] = tg_yy_xxyyy[i] * fxi[i] + 2.0 * tg_xyy_xyyy[i] * fxi[i] + tg_xyy_xxyyy[i] * ra_x[i];

        tg_xxyy_xxyyz[i] = tg_yy_xxyyz[i] * fxi[i] + 2.0 * tg_xyy_xyyz[i] * fxi[i] + tg_xyy_xxyyz[i] * ra_x[i];

        tg_xxyy_xxyzz[i] = tg_yy_xxyzz[i] * fxi[i] + 2.0 * tg_xyy_xyzz[i] * fxi[i] + tg_xyy_xxyzz[i] * ra_x[i];

        tg_xxyy_xxzzz[i] = tg_xx_xxzzz[i] * fxi[i] + tg_xxy_xxzzz[i] * ra_y[i];

        tg_xxyy_xyyyy[i] = tg_yy_xyyyy[i] * fxi[i] + tg_xyy_yyyy[i] * fxi[i] + tg_xyy_xyyyy[i] * ra_x[i];

        tg_xxyy_xyyyz[i] = tg_yy_xyyyz[i] * fxi[i] + tg_xyy_yyyz[i] * fxi[i] + tg_xyy_xyyyz[i] * ra_x[i];

        tg_xxyy_xyyzz[i] = tg_yy_xyyzz[i] * fxi[i] + tg_xyy_yyzz[i] * fxi[i] + tg_xyy_xyyzz[i] * ra_x[i];

        tg_xxyy_xyzzz[i] = tg_yy_xyzzz[i] * fxi[i] + tg_xyy_yzzz[i] * fxi[i] + tg_xyy_xyzzz[i] * ra_x[i];

        tg_xxyy_xzzzz[i] = tg_xx_xzzzz[i] * fxi[i] + tg_xxy_xzzzz[i] * ra_y[i];

        tg_xxyy_yyyyy[i] = tg_yy_yyyyy[i] * fxi[i] + tg_xyy_yyyyy[i] * ra_x[i];

        tg_xxyy_yyyyz[i] = tg_yy_yyyyz[i] * fxi[i] + tg_xyy_yyyyz[i] * ra_x[i];

        tg_xxyy_yyyzz[i] = tg_yy_yyyzz[i] * fxi[i] + tg_xyy_yyyzz[i] * ra_x[i];

        tg_xxyy_yyzzz[i] = tg_yy_yyzzz[i] * fxi[i] + tg_xyy_yyzzz[i] * ra_x[i];

        tg_xxyy_yzzzz[i] = tg_yy_yzzzz[i] * fxi[i] + tg_xyy_yzzzz[i] * ra_x[i];

        tg_xxyy_zzzzz[i] = tg_yy_zzzzz[i] * fxi[i] + tg_xyy_zzzzz[i] * ra_x[i];

        tg_xxyz_xxxxx[i] = tg_xxz_xxxxx[i] * ra_y[i];

        tg_xxyz_xxxxy[i] = tg_xxy_xxxxy[i] * ra_z[i];

        tg_xxyz_xxxxz[i] = tg_xxz_xxxxz[i] * ra_y[i];

        tg_xxyz_xxxyy[i] = tg_xxy_xxxyy[i] * ra_z[i];

        tg_xxyz_xxxyz[i] = tg_xxz_xxxz[i] * fxi[i] + tg_xxz_xxxyz[i] * ra_y[i];

        tg_xxyz_xxxzz[i] = tg_xxz_xxxzz[i] * ra_y[i];

        tg_xxyz_xxyyy[i] = tg_xxy_xxyyy[i] * ra_z[i];

        tg_xxyz_xxyyz[i] = 2.0 * tg_xxz_xxyz[i] * fxi[i] + tg_xxz_xxyyz[i] * ra_y[i];

        tg_xxyz_xxyzz[i] = tg_xxz_xxzz[i] * fxi[i] + tg_xxz_xxyzz[i] * ra_y[i];

        tg_xxyz_xxzzz[i] = tg_xxz_xxzzz[i] * ra_y[i];

        tg_xxyz_xyyyy[i] = tg_xxy_xyyyy[i] * ra_z[i];

        tg_xxyz_xyyyz[i] = 3.0 * tg_xxz_xyyz[i] * fxi[i] + tg_xxz_xyyyz[i] * ra_y[i];

        tg_xxyz_xyyzz[i] = 2.0 * tg_xxz_xyzz[i] * fxi[i] + tg_xxz_xyyzz[i] * ra_y[i];

        tg_xxyz_xyzzz[i] = tg_xxz_xzzz[i] * fxi[i] + tg_xxz_xyzzz[i] * ra_y[i];

        tg_xxyz_xzzzz[i] = tg_xxz_xzzzz[i] * ra_y[i];

        tg_xxyz_yyyyy[i] = tg_xxy_yyyyy[i] * ra_z[i];

        tg_xxyz_yyyyz[i] = tg_yz_yyyyz[i] * fxi[i] + tg_xyz_yyyyz[i] * ra_x[i];

        tg_xxyz_yyyzz[i] = tg_yz_yyyzz[i] * fxi[i] + tg_xyz_yyyzz[i] * ra_x[i];

        tg_xxyz_yyzzz[i] = tg_yz_yyzzz[i] * fxi[i] + tg_xyz_yyzzz[i] * ra_x[i];

        tg_xxyz_yzzzz[i] = tg_yz_yzzzz[i] * fxi[i] + tg_xyz_yzzzz[i] * ra_x[i];

        tg_xxyz_zzzzz[i] = tg_xxz_zzzzz[i] * ra_y[i];

        tg_xxzz_xxxxx[i] = tg_xx_xxxxx[i] * fxi[i] + tg_xxz_xxxxx[i] * ra_z[i];

        tg_xxzz_xxxxy[i] = tg_xx_xxxxy[i] * fxi[i] + tg_xxz_xxxxy[i] * ra_z[i];

        tg_xxzz_xxxxz[i] = tg_zz_xxxxz[i] * fxi[i] + 4.0 * tg_xzz_xxxz[i] * fxi[i] + tg_xzz_xxxxz[i] * ra_x[i];

        tg_xxzz_xxxyy[i] = tg_xx_xxxyy[i] * fxi[i] + tg_xxz_xxxyy[i] * ra_z[i];

        tg_xxzz_xxxyz[i] = tg_zz_xxxyz[i] * fxi[i] + 3.0 * tg_xzz_xxyz[i] * fxi[i] + tg_xzz_xxxyz[i] * ra_x[i];

        tg_xxzz_xxxzz[i] = tg_zz_xxxzz[i] * fxi[i] + 3.0 * tg_xzz_xxzz[i] * fxi[i] + tg_xzz_xxxzz[i] * ra_x[i];

        tg_xxzz_xxyyy[i] = tg_xx_xxyyy[i] * fxi[i] + tg_xxz_xxyyy[i] * ra_z[i];

        tg_xxzz_xxyyz[i] = tg_zz_xxyyz[i] * fxi[i] + 2.0 * tg_xzz_xyyz[i] * fxi[i] + tg_xzz_xxyyz[i] * ra_x[i];

        tg_xxzz_xxyzz[i] = tg_zz_xxyzz[i] * fxi[i] + 2.0 * tg_xzz_xyzz[i] * fxi[i] + tg_xzz_xxyzz[i] * ra_x[i];

        tg_xxzz_xxzzz[i] = tg_zz_xxzzz[i] * fxi[i] + 2.0 * tg_xzz_xzzz[i] * fxi[i] + tg_xzz_xxzzz[i] * ra_x[i];

        tg_xxzz_xyyyy[i] = tg_xx_xyyyy[i] * fxi[i] + tg_xxz_xyyyy[i] * ra_z[i];

        tg_xxzz_xyyyz[i] = tg_zz_xyyyz[i] * fxi[i] + tg_xzz_yyyz[i] * fxi[i] + tg_xzz_xyyyz[i] * ra_x[i];

        tg_xxzz_xyyzz[i] = tg_zz_xyyzz[i] * fxi[i] + tg_xzz_yyzz[i] * fxi[i] + tg_xzz_xyyzz[i] * ra_x[i];

        tg_xxzz_xyzzz[i] = tg_zz_xyzzz[i] * fxi[i] + tg_xzz_yzzz[i] * fxi[i] + tg_xzz_xyzzz[i] * ra_x[i];

        tg_xxzz_xzzzz[i] = tg_zz_xzzzz[i] * fxi[i] + tg_xzz_zzzz[i] * fxi[i] + tg_xzz_xzzzz[i] * ra_x[i];

        tg_xxzz_yyyyy[i] = tg_zz_yyyyy[i] * fxi[i] + tg_xzz_yyyyy[i] * ra_x[i];

        tg_xxzz_yyyyz[i] = tg_zz_yyyyz[i] * fxi[i] + tg_xzz_yyyyz[i] * ra_x[i];

        tg_xxzz_yyyzz[i] = tg_zz_yyyzz[i] * fxi[i] + tg_xzz_yyyzz[i] * ra_x[i];

        tg_xxzz_yyzzz[i] = tg_zz_yyzzz[i] * fxi[i] + tg_xzz_yyzzz[i] * ra_x[i];

        tg_xxzz_yzzzz[i] = tg_zz_yzzzz[i] * fxi[i] + tg_xzz_yzzzz[i] * ra_x[i];

        tg_xxzz_zzzzz[i] = tg_zz_zzzzz[i] * fxi[i] + tg_xzz_zzzzz[i] * ra_x[i];

        tg_xyyy_xxxxx[i] = 5.0 * tg_yyy_xxxx[i] * fxi[i] + tg_yyy_xxxxx[i] * ra_x[i];

        tg_xyyy_xxxxy[i] = 4.0 * tg_yyy_xxxy[i] * fxi[i] + tg_yyy_xxxxy[i] * ra_x[i];

        tg_xyyy_xxxxz[i] = 4.0 * tg_yyy_xxxz[i] * fxi[i] + tg_yyy_xxxxz[i] * ra_x[i];

        tg_xyyy_xxxyy[i] = 3.0 * tg_yyy_xxyy[i] * fxi[i] + tg_yyy_xxxyy[i] * ra_x[i];

        tg_xyyy_xxxyz[i] = 3.0 * tg_yyy_xxyz[i] * fxi[i] + tg_yyy_xxxyz[i] * ra_x[i];

        tg_xyyy_xxxzz[i] = 3.0 * tg_yyy_xxzz[i] * fxi[i] + tg_yyy_xxxzz[i] * ra_x[i];

        tg_xyyy_xxyyy[i] = 2.0 * tg_yyy_xyyy[i] * fxi[i] + tg_yyy_xxyyy[i] * ra_x[i];

        tg_xyyy_xxyyz[i] = 2.0 * tg_yyy_xyyz[i] * fxi[i] + tg_yyy_xxyyz[i] * ra_x[i];

        tg_xyyy_xxyzz[i] = 2.0 * tg_yyy_xyzz[i] * fxi[i] + tg_yyy_xxyzz[i] * ra_x[i];

        tg_xyyy_xxzzz[i] = 2.0 * tg_yyy_xzzz[i] * fxi[i] + tg_yyy_xxzzz[i] * ra_x[i];

        tg_xyyy_xyyyy[i] = tg_yyy_yyyy[i] * fxi[i] + tg_yyy_xyyyy[i] * ra_x[i];

        tg_xyyy_xyyyz[i] = tg_yyy_yyyz[i] * fxi[i] + tg_yyy_xyyyz[i] * ra_x[i];

        tg_xyyy_xyyzz[i] = tg_yyy_yyzz[i] * fxi[i] + tg_yyy_xyyzz[i] * ra_x[i];

        tg_xyyy_xyzzz[i] = tg_yyy_yzzz[i] * fxi[i] + tg_yyy_xyzzz[i] * ra_x[i];

        tg_xyyy_xzzzz[i] = tg_yyy_zzzz[i] * fxi[i] + tg_yyy_xzzzz[i] * ra_x[i];

        tg_xyyy_yyyyy[i] = tg_yyy_yyyyy[i] * ra_x[i];

        tg_xyyy_yyyyz[i] = tg_yyy_yyyyz[i] * ra_x[i];

        tg_xyyy_yyyzz[i] = tg_yyy_yyyzz[i] * ra_x[i];

        tg_xyyy_yyzzz[i] = tg_yyy_yyzzz[i] * ra_x[i];

        tg_xyyy_yzzzz[i] = tg_yyy_yzzzz[i] * ra_x[i];

        tg_xyyy_zzzzz[i] = tg_yyy_zzzzz[i] * ra_x[i];

        tg_xyyz_xxxxx[i] = tg_xyy_xxxxx[i] * ra_z[i];

        tg_xyyz_xxxxy[i] = tg_xyy_xxxxy[i] * ra_z[i];

        tg_xyyz_xxxxz[i] = 4.0 * tg_yyz_xxxz[i] * fxi[i] + tg_yyz_xxxxz[i] * ra_x[i];

        tg_xyyz_xxxyy[i] = tg_xyy_xxxyy[i] * ra_z[i];

        tg_xyyz_xxxyz[i] = 3.0 * tg_yyz_xxyz[i] * fxi[i] + tg_yyz_xxxyz[i] * ra_x[i];

        tg_xyyz_xxxzz[i] = 3.0 * tg_yyz_xxzz[i] * fxi[i] + tg_yyz_xxxzz[i] * ra_x[i];

        tg_xyyz_xxyyy[i] = tg_xyy_xxyyy[i] * ra_z[i];

        tg_xyyz_xxyyz[i] = 2.0 * tg_yyz_xyyz[i] * fxi[i] + tg_yyz_xxyyz[i] * ra_x[i];

        tg_xyyz_xxyzz[i] = 2.0 * tg_yyz_xyzz[i] * fxi[i] + tg_yyz_xxyzz[i] * ra_x[i];

        tg_xyyz_xxzzz[i] = 2.0 * tg_yyz_xzzz[i] * fxi[i] + tg_yyz_xxzzz[i] * ra_x[i];

        tg_xyyz_xyyyy[i] = tg_xyy_xyyyy[i] * ra_z[i];

        tg_xyyz_xyyyz[i] = tg_yyz_yyyz[i] * fxi[i] + tg_yyz_xyyyz[i] * ra_x[i];

        tg_xyyz_xyyzz[i] = tg_yyz_yyzz[i] * fxi[i] + tg_yyz_xyyzz[i] * ra_x[i];

        tg_xyyz_xyzzz[i] = tg_yyz_yzzz[i] * fxi[i] + tg_yyz_xyzzz[i] * ra_x[i];

        tg_xyyz_xzzzz[i] = tg_yyz_zzzz[i] * fxi[i] + tg_yyz_xzzzz[i] * ra_x[i];

        tg_xyyz_yyyyy[i] = tg_yyz_yyyyy[i] * ra_x[i];

        tg_xyyz_yyyyz[i] = tg_yyz_yyyyz[i] * ra_x[i];

        tg_xyyz_yyyzz[i] = tg_yyz_yyyzz[i] * ra_x[i];

        tg_xyyz_yyzzz[i] = tg_yyz_yyzzz[i] * ra_x[i];

        tg_xyyz_yzzzz[i] = tg_yyz_yzzzz[i] * ra_x[i];

        tg_xyyz_zzzzz[i] = tg_yyz_zzzzz[i] * ra_x[i];

        tg_xyzz_xxxxx[i] = tg_xzz_xxxxx[i] * ra_y[i];

        tg_xyzz_xxxxy[i] = 4.0 * tg_yzz_xxxy[i] * fxi[i] + tg_yzz_xxxxy[i] * ra_x[i];

        tg_xyzz_xxxxz[i] = tg_xzz_xxxxz[i] * ra_y[i];

        tg_xyzz_xxxyy[i] = 3.0 * tg_yzz_xxyy[i] * fxi[i] + tg_yzz_xxxyy[i] * ra_x[i];

        tg_xyzz_xxxyz[i] = 3.0 * tg_yzz_xxyz[i] * fxi[i] + tg_yzz_xxxyz[i] * ra_x[i];

        tg_xyzz_xxxzz[i] = tg_xzz_xxxzz[i] * ra_y[i];

        tg_xyzz_xxyyy[i] = 2.0 * tg_yzz_xyyy[i] * fxi[i] + tg_yzz_xxyyy[i] * ra_x[i];

        tg_xyzz_xxyyz[i] = 2.0 * tg_yzz_xyyz[i] * fxi[i] + tg_yzz_xxyyz[i] * ra_x[i];

        tg_xyzz_xxyzz[i] = 2.0 * tg_yzz_xyzz[i] * fxi[i] + tg_yzz_xxyzz[i] * ra_x[i];

        tg_xyzz_xxzzz[i] = tg_xzz_xxzzz[i] * ra_y[i];

        tg_xyzz_xyyyy[i] = tg_yzz_yyyy[i] * fxi[i] + tg_yzz_xyyyy[i] * ra_x[i];

        tg_xyzz_xyyyz[i] = tg_yzz_yyyz[i] * fxi[i] + tg_yzz_xyyyz[i] * ra_x[i];

        tg_xyzz_xyyzz[i] = tg_yzz_yyzz[i] * fxi[i] + tg_yzz_xyyzz[i] * ra_x[i];

        tg_xyzz_xyzzz[i] = tg_yzz_yzzz[i] * fxi[i] + tg_yzz_xyzzz[i] * ra_x[i];

        tg_xyzz_xzzzz[i] = tg_xzz_xzzzz[i] * ra_y[i];

        tg_xyzz_yyyyy[i] = tg_yzz_yyyyy[i] * ra_x[i];

        tg_xyzz_yyyyz[i] = tg_yzz_yyyyz[i] * ra_x[i];

        tg_xyzz_yyyzz[i] = tg_yzz_yyyzz[i] * ra_x[i];

        tg_xyzz_yyzzz[i] = tg_yzz_yyzzz[i] * ra_x[i];

        tg_xyzz_yzzzz[i] = tg_yzz_yzzzz[i] * ra_x[i];

        tg_xyzz_zzzzz[i] = tg_yzz_zzzzz[i] * ra_x[i];

        tg_xzzz_xxxxx[i] = 5.0 * tg_zzz_xxxx[i] * fxi[i] + tg_zzz_xxxxx[i] * ra_x[i];

        tg_xzzz_xxxxy[i] = 4.0 * tg_zzz_xxxy[i] * fxi[i] + tg_zzz_xxxxy[i] * ra_x[i];

        tg_xzzz_xxxxz[i] = 4.0 * tg_zzz_xxxz[i] * fxi[i] + tg_zzz_xxxxz[i] * ra_x[i];

        tg_xzzz_xxxyy[i] = 3.0 * tg_zzz_xxyy[i] * fxi[i] + tg_zzz_xxxyy[i] * ra_x[i];

        tg_xzzz_xxxyz[i] = 3.0 * tg_zzz_xxyz[i] * fxi[i] + tg_zzz_xxxyz[i] * ra_x[i];

        tg_xzzz_xxxzz[i] = 3.0 * tg_zzz_xxzz[i] * fxi[i] + tg_zzz_xxxzz[i] * ra_x[i];

        tg_xzzz_xxyyy[i] = 2.0 * tg_zzz_xyyy[i] * fxi[i] + tg_zzz_xxyyy[i] * ra_x[i];

        tg_xzzz_xxyyz[i] = 2.0 * tg_zzz_xyyz[i] * fxi[i] + tg_zzz_xxyyz[i] * ra_x[i];

        tg_xzzz_xxyzz[i] = 2.0 * tg_zzz_xyzz[i] * fxi[i] + tg_zzz_xxyzz[i] * ra_x[i];

        tg_xzzz_xxzzz[i] = 2.0 * tg_zzz_xzzz[i] * fxi[i] + tg_zzz_xxzzz[i] * ra_x[i];

        tg_xzzz_xyyyy[i] = tg_zzz_yyyy[i] * fxi[i] + tg_zzz_xyyyy[i] * ra_x[i];

        tg_xzzz_xyyyz[i] = tg_zzz_yyyz[i] * fxi[i] + tg_zzz_xyyyz[i] * ra_x[i];

        tg_xzzz_xyyzz[i] = tg_zzz_yyzz[i] * fxi[i] + tg_zzz_xyyzz[i] * ra_x[i];

        tg_xzzz_xyzzz[i] = tg_zzz_yzzz[i] * fxi[i] + tg_zzz_xyzzz[i] * ra_x[i];

        tg_xzzz_xzzzz[i] = tg_zzz_zzzz[i] * fxi[i] + tg_zzz_xzzzz[i] * ra_x[i];

        tg_xzzz_yyyyy[i] = tg_zzz_yyyyy[i] * ra_x[i];

        tg_xzzz_yyyyz[i] = tg_zzz_yyyyz[i] * ra_x[i];

        tg_xzzz_yyyzz[i] = tg_zzz_yyyzz[i] * ra_x[i];

        tg_xzzz_yyzzz[i] = tg_zzz_yyzzz[i] * ra_x[i];

        tg_xzzz_yzzzz[i] = tg_zzz_yzzzz[i] * ra_x[i];

        tg_xzzz_zzzzz[i] = tg_zzz_zzzzz[i] * ra_x[i];

        tg_yyyy_xxxxx[i] = 3.0 * tg_yy_xxxxx[i] * fxi[i] + tg_yyy_xxxxx[i] * ra_y[i];

        tg_yyyy_xxxxy[i] = 3.0 * tg_yy_xxxxy[i] * fxi[i] + tg_yyy_xxxx[i] * fxi[i] + tg_yyy_xxxxy[i] * ra_y[i];

        tg_yyyy_xxxxz[i] = 3.0 * tg_yy_xxxxz[i] * fxi[i] + tg_yyy_xxxxz[i] * ra_y[i];

        tg_yyyy_xxxyy[i] = 3.0 * tg_yy_xxxyy[i] * fxi[i] + 2.0 * tg_yyy_xxxy[i] * fxi[i] + tg_yyy_xxxyy[i] * ra_y[i];

        tg_yyyy_xxxyz[i] = 3.0 * tg_yy_xxxyz[i] * fxi[i] + tg_yyy_xxxz[i] * fxi[i] + tg_yyy_xxxyz[i] * ra_y[i];

        tg_yyyy_xxxzz[i] = 3.0 * tg_yy_xxxzz[i] * fxi[i] + tg_yyy_xxxzz[i] * ra_y[i];

        tg_yyyy_xxyyy[i] = 3.0 * tg_yy_xxyyy[i] * fxi[i] + 3.0 * tg_yyy_xxyy[i] * fxi[i] + tg_yyy_xxyyy[i] * ra_y[i];

        tg_yyyy_xxyyz[i] = 3.0 * tg_yy_xxyyz[i] * fxi[i] + 2.0 * tg_yyy_xxyz[i] * fxi[i] + tg_yyy_xxyyz[i] * ra_y[i];

        tg_yyyy_xxyzz[i] = 3.0 * tg_yy_xxyzz[i] * fxi[i] + tg_yyy_xxzz[i] * fxi[i] + tg_yyy_xxyzz[i] * ra_y[i];

        tg_yyyy_xxzzz[i] = 3.0 * tg_yy_xxzzz[i] * fxi[i] + tg_yyy_xxzzz[i] * ra_y[i];

        tg_yyyy_xyyyy[i] = 3.0 * tg_yy_xyyyy[i] * fxi[i] + 4.0 * tg_yyy_xyyy[i] * fxi[i] + tg_yyy_xyyyy[i] * ra_y[i];

        tg_yyyy_xyyyz[i] = 3.0 * tg_yy_xyyyz[i] * fxi[i] + 3.0 * tg_yyy_xyyz[i] * fxi[i] + tg_yyy_xyyyz[i] * ra_y[i];

        tg_yyyy_xyyzz[i] = 3.0 * tg_yy_xyyzz[i] * fxi[i] + 2.0 * tg_yyy_xyzz[i] * fxi[i] + tg_yyy_xyyzz[i] * ra_y[i];

        tg_yyyy_xyzzz[i] = 3.0 * tg_yy_xyzzz[i] * fxi[i] + tg_yyy_xzzz[i] * fxi[i] + tg_yyy_xyzzz[i] * ra_y[i];

        tg_yyyy_xzzzz[i] = 3.0 * tg_yy_xzzzz[i] * fxi[i] + tg_yyy_xzzzz[i] * ra_y[i];

        tg_yyyy_yyyyy[i] = 3.0 * tg_yy_yyyyy[i] * fxi[i] + 5.0 * tg_yyy_yyyy[i] * fxi[i] + tg_yyy_yyyyy[i] * ra_y[i];

        tg_yyyy_yyyyz[i] = 3.0 * tg_yy_yyyyz[i] * fxi[i] + 4.0 * tg_yyy_yyyz[i] * fxi[i] + tg_yyy_yyyyz[i] * ra_y[i];

        tg_yyyy_yyyzz[i] = 3.0 * tg_yy_yyyzz[i] * fxi[i] + 3.0 * tg_yyy_yyzz[i] * fxi[i] + tg_yyy_yyyzz[i] * ra_y[i];

        tg_yyyy_yyzzz[i] = 3.0 * tg_yy_yyzzz[i] * fxi[i] + 2.0 * tg_yyy_yzzz[i] * fxi[i] + tg_yyy_yyzzz[i] * ra_y[i];

        tg_yyyy_yzzzz[i] = 3.0 * tg_yy_yzzzz[i] * fxi[i] + tg_yyy_zzzz[i] * fxi[i] + tg_yyy_yzzzz[i] * ra_y[i];

        tg_yyyy_zzzzz[i] = 3.0 * tg_yy_zzzzz[i] * fxi[i] + tg_yyy_zzzzz[i] * ra_y[i];

        tg_yyyz_xxxxx[i] = tg_yyy_xxxxx[i] * ra_z[i];

        tg_yyyz_xxxxy[i] = tg_yyy_xxxxy[i] * ra_z[i];

        tg_yyyz_xxxxz[i] = 2.0 * tg_yz_xxxxz[i] * fxi[i] + tg_yyz_xxxxz[i] * ra_y[i];

        tg_yyyz_xxxyy[i] = tg_yyy_xxxyy[i] * ra_z[i];

        tg_yyyz_xxxyz[i] = tg_yyy_xxxy[i] * fxi[i] + tg_yyy_xxxyz[i] * ra_z[i];

        tg_yyyz_xxxzz[i] = 2.0 * tg_yz_xxxzz[i] * fxi[i] + tg_yyz_xxxzz[i] * ra_y[i];

        tg_yyyz_xxyyy[i] = tg_yyy_xxyyy[i] * ra_z[i];

        tg_yyyz_xxyyz[i] = tg_yyy_xxyy[i] * fxi[i] + tg_yyy_xxyyz[i] * ra_z[i];

        tg_yyyz_xxyzz[i] = 2.0 * tg_yyy_xxyz[i] * fxi[i] + tg_yyy_xxyzz[i] * ra_z[i];

        tg_yyyz_xxzzz[i] = 2.0 * tg_yz_xxzzz[i] * fxi[i] + tg_yyz_xxzzz[i] * ra_y[i];

        tg_yyyz_xyyyy[i] = tg_yyy_xyyyy[i] * ra_z[i];

        tg_yyyz_xyyyz[i] = tg_yyy_xyyy[i] * fxi[i] + tg_yyy_xyyyz[i] * ra_z[i];

        tg_yyyz_xyyzz[i] = 2.0 * tg_yyy_xyyz[i] * fxi[i] + tg_yyy_xyyzz[i] * ra_z[i];

        tg_yyyz_xyzzz[i] = 3.0 * tg_yyy_xyzz[i] * fxi[i] + tg_yyy_xyzzz[i] * ra_z[i];

        tg_yyyz_xzzzz[i] = 2.0 * tg_yz_xzzzz[i] * fxi[i] + tg_yyz_xzzzz[i] * ra_y[i];

        tg_yyyz_yyyyy[i] = tg_yyy_yyyyy[i] * ra_z[i];

        tg_yyyz_yyyyz[i] = tg_yyy_yyyy[i] * fxi[i] + tg_yyy_yyyyz[i] * ra_z[i];

        tg_yyyz_yyyzz[i] = 2.0 * tg_yyy_yyyz[i] * fxi[i] + tg_yyy_yyyzz[i] * ra_z[i];

        tg_yyyz_yyzzz[i] = 3.0 * tg_yyy_yyzz[i] * fxi[i] + tg_yyy_yyzzz[i] * ra_z[i];

        tg_yyyz_yzzzz[i] = 4.0 * tg_yyy_yzzz[i] * fxi[i] + tg_yyy_yzzzz[i] * ra_z[i];

        tg_yyyz_zzzzz[i] = 2.0 * tg_yz_zzzzz[i] * fxi[i] + tg_yyz_zzzzz[i] * ra_y[i];

        tg_yyzz_xxxxx[i] = tg_zz_xxxxx[i] * fxi[i] + tg_yzz_xxxxx[i] * ra_y[i];

        tg_yyzz_xxxxy[i] = tg_yy_xxxxy[i] * fxi[i] + tg_yyz_xxxxy[i] * ra_z[i];

        tg_yyzz_xxxxz[i] = tg_zz_xxxxz[i] * fxi[i] + tg_yzz_xxxxz[i] * ra_y[i];

        tg_yyzz_xxxyy[i] = tg_yy_xxxyy[i] * fxi[i] + tg_yyz_xxxyy[i] * ra_z[i];

        tg_yyzz_xxxyz[i] = tg_zz_xxxyz[i] * fxi[i] + tg_yzz_xxxz[i] * fxi[i] + tg_yzz_xxxyz[i] * ra_y[i];

        tg_yyzz_xxxzz[i] = tg_zz_xxxzz[i] * fxi[i] + tg_yzz_xxxzz[i] * ra_y[i];

        tg_yyzz_xxyyy[i] = tg_yy_xxyyy[i] * fxi[i] + tg_yyz_xxyyy[i] * ra_z[i];

        tg_yyzz_xxyyz[i] = tg_zz_xxyyz[i] * fxi[i] + 2.0 * tg_yzz_xxyz[i] * fxi[i] + tg_yzz_xxyyz[i] * ra_y[i];

        tg_yyzz_xxyzz[i] = tg_zz_xxyzz[i] * fxi[i] + tg_yzz_xxzz[i] * fxi[i] + tg_yzz_xxyzz[i] * ra_y[i];

        tg_yyzz_xxzzz[i] = tg_zz_xxzzz[i] * fxi[i] + tg_yzz_xxzzz[i] * ra_y[i];

        tg_yyzz_xyyyy[i] = tg_yy_xyyyy[i] * fxi[i] + tg_yyz_xyyyy[i] * ra_z[i];

        tg_yyzz_xyyyz[i] = tg_zz_xyyyz[i] * fxi[i] + 3.0 * tg_yzz_xyyz[i] * fxi[i] + tg_yzz_xyyyz[i] * ra_y[i];

        tg_yyzz_xyyzz[i] = tg_zz_xyyzz[i] * fxi[i] + 2.0 * tg_yzz_xyzz[i] * fxi[i] + tg_yzz_xyyzz[i] * ra_y[i];

        tg_yyzz_xyzzz[i] = tg_zz_xyzzz[i] * fxi[i] + tg_yzz_xzzz[i] * fxi[i] + tg_yzz_xyzzz[i] * ra_y[i];

        tg_yyzz_xzzzz[i] = tg_zz_xzzzz[i] * fxi[i] + tg_yzz_xzzzz[i] * ra_y[i];

        tg_yyzz_yyyyy[i] = tg_yy_yyyyy[i] * fxi[i] + tg_yyz_yyyyy[i] * ra_z[i];

        tg_yyzz_yyyyz[i] = tg_zz_yyyyz[i] * fxi[i] + 4.0 * tg_yzz_yyyz[i] * fxi[i] + tg_yzz_yyyyz[i] * ra_y[i];

        tg_yyzz_yyyzz[i] = tg_zz_yyyzz[i] * fxi[i] + 3.0 * tg_yzz_yyzz[i] * fxi[i] + tg_yzz_yyyzz[i] * ra_y[i];

        tg_yyzz_yyzzz[i] = tg_zz_yyzzz[i] * fxi[i] + 2.0 * tg_yzz_yzzz[i] * fxi[i] + tg_yzz_yyzzz[i] * ra_y[i];

        tg_yyzz_yzzzz[i] = tg_zz_yzzzz[i] * fxi[i] + tg_yzz_zzzz[i] * fxi[i] + tg_yzz_yzzzz[i] * ra_y[i];

        tg_yyzz_zzzzz[i] = tg_zz_zzzzz[i] * fxi[i] + tg_yzz_zzzzz[i] * ra_y[i];

        tg_yzzz_xxxxx[i] = tg_zzz_xxxxx[i] * ra_y[i];

        tg_yzzz_xxxxy[i] = tg_zzz_xxxx[i] * fxi[i] + tg_zzz_xxxxy[i] * ra_y[i];

        tg_yzzz_xxxxz[i] = tg_zzz_xxxxz[i] * ra_y[i];

        tg_yzzz_xxxyy[i] = 2.0 * tg_zzz_xxxy[i] * fxi[i] + tg_zzz_xxxyy[i] * ra_y[i];

        tg_yzzz_xxxyz[i] = tg_zzz_xxxz[i] * fxi[i] + tg_zzz_xxxyz[i] * ra_y[i];

        tg_yzzz_xxxzz[i] = tg_zzz_xxxzz[i] * ra_y[i];

        tg_yzzz_xxyyy[i] = 3.0 * tg_zzz_xxyy[i] * fxi[i] + tg_zzz_xxyyy[i] * ra_y[i];

        tg_yzzz_xxyyz[i] = 2.0 * tg_zzz_xxyz[i] * fxi[i] + tg_zzz_xxyyz[i] * ra_y[i];

        tg_yzzz_xxyzz[i] = tg_zzz_xxzz[i] * fxi[i] + tg_zzz_xxyzz[i] * ra_y[i];

        tg_yzzz_xxzzz[i] = tg_zzz_xxzzz[i] * ra_y[i];

        tg_yzzz_xyyyy[i] = 4.0 * tg_zzz_xyyy[i] * fxi[i] + tg_zzz_xyyyy[i] * ra_y[i];

        tg_yzzz_xyyyz[i] = 3.0 * tg_zzz_xyyz[i] * fxi[i] + tg_zzz_xyyyz[i] * ra_y[i];

        tg_yzzz_xyyzz[i] = 2.0 * tg_zzz_xyzz[i] * fxi[i] + tg_zzz_xyyzz[i] * ra_y[i];

        tg_yzzz_xyzzz[i] = tg_zzz_xzzz[i] * fxi[i] + tg_zzz_xyzzz[i] * ra_y[i];

        tg_yzzz_xzzzz[i] = tg_zzz_xzzzz[i] * ra_y[i];

        tg_yzzz_yyyyy[i] = 5.0 * tg_zzz_yyyy[i] * fxi[i] + tg_zzz_yyyyy[i] * ra_y[i];

        tg_yzzz_yyyyz[i] = 4.0 * tg_zzz_yyyz[i] * fxi[i] + tg_zzz_yyyyz[i] * ra_y[i];

        tg_yzzz_yyyzz[i] = 3.0 * tg_zzz_yyzz[i] * fxi[i] + tg_zzz_yyyzz[i] * ra_y[i];

        tg_yzzz_yyzzz[i] = 2.0 * tg_zzz_yzzz[i] * fxi[i] + tg_zzz_yyzzz[i] * ra_y[i];

        tg_yzzz_yzzzz[i] = tg_zzz_zzzz[i] * fxi[i] + tg_zzz_yzzzz[i] * ra_y[i];

        tg_yzzz_zzzzz[i] = tg_zzz_zzzzz[i] * ra_y[i];

        tg_zzzz_xxxxx[i] = 3.0 * tg_zz_xxxxx[i] * fxi[i] + tg_zzz_xxxxx[i] * ra_z[i];

        tg_zzzz_xxxxy[i] = 3.0 * tg_zz_xxxxy[i] * fxi[i] + tg_zzz_xxxxy[i] * ra_z[i];

        tg_zzzz_xxxxz[i] = 3.0 * tg_zz_xxxxz[i] * fxi[i] + tg_zzz_xxxx[i] * fxi[i] + tg_zzz_xxxxz[i] * ra_z[i];

        tg_zzzz_xxxyy[i] = 3.0 * tg_zz_xxxyy[i] * fxi[i] + tg_zzz_xxxyy[i] * ra_z[i];

        tg_zzzz_xxxyz[i] = 3.0 * tg_zz_xxxyz[i] * fxi[i] + tg_zzz_xxxy[i] * fxi[i] + tg_zzz_xxxyz[i] * ra_z[i];

        tg_zzzz_xxxzz[i] = 3.0 * tg_zz_xxxzz[i] * fxi[i] + 2.0 * tg_zzz_xxxz[i] * fxi[i] + tg_zzz_xxxzz[i] * ra_z[i];

        tg_zzzz_xxyyy[i] = 3.0 * tg_zz_xxyyy[i] * fxi[i] + tg_zzz_xxyyy[i] * ra_z[i];

        tg_zzzz_xxyyz[i] = 3.0 * tg_zz_xxyyz[i] * fxi[i] + tg_zzz_xxyy[i] * fxi[i] + tg_zzz_xxyyz[i] * ra_z[i];

        tg_zzzz_xxyzz[i] = 3.0 * tg_zz_xxyzz[i] * fxi[i] + 2.0 * tg_zzz_xxyz[i] * fxi[i] + tg_zzz_xxyzz[i] * ra_z[i];

        tg_zzzz_xxzzz[i] = 3.0 * tg_zz_xxzzz[i] * fxi[i] + 3.0 * tg_zzz_xxzz[i] * fxi[i] + tg_zzz_xxzzz[i] * ra_z[i];

        tg_zzzz_xyyyy[i] = 3.0 * tg_zz_xyyyy[i] * fxi[i] + tg_zzz_xyyyy[i] * ra_z[i];

        tg_zzzz_xyyyz[i] = 3.0 * tg_zz_xyyyz[i] * fxi[i] + tg_zzz_xyyy[i] * fxi[i] + tg_zzz_xyyyz[i] * ra_z[i];

        tg_zzzz_xyyzz[i] = 3.0 * tg_zz_xyyzz[i] * fxi[i] + 2.0 * tg_zzz_xyyz[i] * fxi[i] + tg_zzz_xyyzz[i] * ra_z[i];

        tg_zzzz_xyzzz[i] = 3.0 * tg_zz_xyzzz[i] * fxi[i] + 3.0 * tg_zzz_xyzz[i] * fxi[i] + tg_zzz_xyzzz[i] * ra_z[i];

        tg_zzzz_xzzzz[i] = 3.0 * tg_zz_xzzzz[i] * fxi[i] + 4.0 * tg_zzz_xzzz[i] * fxi[i] + tg_zzz_xzzzz[i] * ra_z[i];

        tg_zzzz_yyyyy[i] = 3.0 * tg_zz_yyyyy[i] * fxi[i] + tg_zzz_yyyyy[i] * ra_z[i];

        tg_zzzz_yyyyz[i] = 3.0 * tg_zz_yyyyz[i] * fxi[i] + tg_zzz_yyyy[i] * fxi[i] + tg_zzz_yyyyz[i] * ra_z[i];

        tg_zzzz_yyyzz[i] = 3.0 * tg_zz_yyyzz[i] * fxi[i] + 2.0 * tg_zzz_yyyz[i] * fxi[i] + tg_zzz_yyyzz[i] * ra_z[i];

        tg_zzzz_yyzzz[i] = 3.0 * tg_zz_yyzzz[i] * fxi[i] + 3.0 * tg_zzz_yyzz[i] * fxi[i] + tg_zzz_yyzzz[i] * ra_z[i];

        tg_zzzz_yzzzz[i] = 3.0 * tg_zz_yzzzz[i] * fxi[i] + 4.0 * tg_zzz_yzzz[i] * fxi[i] + tg_zzz_yzzzz[i] * ra_z[i];

        tg_zzzz_zzzzz[i] = 3.0 * tg_zz_zzzzz[i] * fxi[i] + 5.0 * tg_zzz_zzzz[i] * fxi[i] + tg_zzz_zzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

