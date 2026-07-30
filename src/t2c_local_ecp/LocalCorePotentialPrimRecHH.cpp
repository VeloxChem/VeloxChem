#include "LocalCorePotentialPrimRecHH.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_hh(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hh,
                                  const size_t idx_fh,
                                  const size_t idx_gg,
                                  const size_t idx_gh,
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

    auto tg_xxy_xxxxz = pbuffer.data(idx_fh + 23);

    auto tg_xxy_xxxzz = pbuffer.data(idx_fh + 26);

    auto tg_xxy_xxzzz = pbuffer.data(idx_fh + 30);

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

    auto tg_xxz_xxxzz = pbuffer.data(idx_fh + 47);

    auto tg_xxz_xxyyy = pbuffer.data(idx_fh + 48);

    auto tg_xxz_xxzzz = pbuffer.data(idx_fh + 51);

    auto tg_xxz_xyyyy = pbuffer.data(idx_fh + 52);

    auto tg_xxz_xzzzz = pbuffer.data(idx_fh + 56);

    auto tg_xxz_yyyyz = pbuffer.data(idx_fh + 58);

    auto tg_xxz_yyyzz = pbuffer.data(idx_fh + 59);

    auto tg_xxz_yyzzz = pbuffer.data(idx_fh + 60);

    auto tg_xxz_yzzzz = pbuffer.data(idx_fh + 61);

    auto tg_xxz_zzzzz = pbuffer.data(idx_fh + 62);

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

    auto tg_yyz_xxxzz = pbuffer.data(idx_fh + 152);

    auto tg_yyz_xxyyy = pbuffer.data(idx_fh + 153);

    auto tg_yyz_xxzzz = pbuffer.data(idx_fh + 156);

    auto tg_yyz_xyyyy = pbuffer.data(idx_fh + 157);

    auto tg_yyz_xzzzz = pbuffer.data(idx_fh + 161);

    auto tg_yyz_yyyyy = pbuffer.data(idx_fh + 162);

    auto tg_yyz_yyyyz = pbuffer.data(idx_fh + 163);

    auto tg_yyz_yyyzz = pbuffer.data(idx_fh + 164);

    auto tg_yyz_yyzzz = pbuffer.data(idx_fh + 165);

    auto tg_yyz_yzzzz = pbuffer.data(idx_fh + 166);

    auto tg_yyz_zzzzz = pbuffer.data(idx_fh + 167);

    auto tg_yzz_xxxxx = pbuffer.data(idx_fh + 168);

    auto tg_yzz_xxxxz = pbuffer.data(idx_fh + 170);

    auto tg_yzz_xxxyz = pbuffer.data(idx_fh + 172);

    auto tg_yzz_xxxzz = pbuffer.data(idx_fh + 173);

    auto tg_yzz_xxyyz = pbuffer.data(idx_fh + 175);

    auto tg_yzz_xxyzz = pbuffer.data(idx_fh + 176);

    auto tg_yzz_xxzzz = pbuffer.data(idx_fh + 177);

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

    auto tg_xxxz_xxxz = pbuffer.data(idx_gg + 32);

    auto tg_xxxz_xxyz = pbuffer.data(idx_gg + 34);

    auto tg_xxxz_xxzz = pbuffer.data(idx_gg + 35);

    auto tg_xxxz_xyyz = pbuffer.data(idx_gg + 37);

    auto tg_xxxz_xyzz = pbuffer.data(idx_gg + 38);

    auto tg_xxxz_xzzz = pbuffer.data(idx_gg + 39);

    auto tg_xxyy_xxxy = pbuffer.data(idx_gg + 46);

    auto tg_xxyy_xxyy = pbuffer.data(idx_gg + 48);

    auto tg_xxyy_xxyz = pbuffer.data(idx_gg + 49);

    auto tg_xxyy_xyyy = pbuffer.data(idx_gg + 51);

    auto tg_xxyy_xyyz = pbuffer.data(idx_gg + 52);

    auto tg_xxyy_xyzz = pbuffer.data(idx_gg + 53);

    auto tg_xxyy_yyyy = pbuffer.data(idx_gg + 55);

    auto tg_xxyy_yyyz = pbuffer.data(idx_gg + 56);

    auto tg_xxyy_yyzz = pbuffer.data(idx_gg + 57);

    auto tg_xxyy_yzzz = pbuffer.data(idx_gg + 58);

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

    auto tg_xxzz_yyyz = pbuffer.data(idx_gg + 86);

    auto tg_xxzz_yyzz = pbuffer.data(idx_gg + 87);

    auto tg_xxzz_yzzz = pbuffer.data(idx_gg + 88);

    auto tg_xxzz_zzzz = pbuffer.data(idx_gg + 89);

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

    auto tg_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto tg_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto tg_xzzz_xxzz = pbuffer.data(idx_gg + 140);

    auto tg_xzzz_xyyz = pbuffer.data(idx_gg + 142);

    auto tg_xzzz_xyzz = pbuffer.data(idx_gg + 143);

    auto tg_xzzz_xzzz = pbuffer.data(idx_gg + 144);

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

    auto tg_yyyz_xxxz = pbuffer.data(idx_gg + 167);

    auto tg_yyyz_xxyz = pbuffer.data(idx_gg + 169);

    auto tg_yyyz_xxzz = pbuffer.data(idx_gg + 170);

    auto tg_yyyz_xyyz = pbuffer.data(idx_gg + 172);

    auto tg_yyyz_xyzz = pbuffer.data(idx_gg + 173);

    auto tg_yyyz_xzzz = pbuffer.data(idx_gg + 174);

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

    // Set up components of auxiliary buffer : GH

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

    auto tg_xxxy_xxxzz = pbuffer.data(idx_gh + 26);

    auto tg_xxxy_xxyyy = pbuffer.data(idx_gh + 27);

    auto tg_xxxy_xxzzz = pbuffer.data(idx_gh + 30);

    auto tg_xxxy_xyyyy = pbuffer.data(idx_gh + 31);

    auto tg_xxxy_xzzzz = pbuffer.data(idx_gh + 35);

    auto tg_xxxy_yyyyy = pbuffer.data(idx_gh + 36);

    auto tg_xxxy_yyyyz = pbuffer.data(idx_gh + 37);

    auto tg_xxxy_yyyzz = pbuffer.data(idx_gh + 38);

    auto tg_xxxy_yyzzz = pbuffer.data(idx_gh + 39);

    auto tg_xxxy_yzzzz = pbuffer.data(idx_gh + 40);

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

    auto tg_xxyz_xxxxz = pbuffer.data(idx_gh + 86);

    auto tg_xxyz_xxxzz = pbuffer.data(idx_gh + 89);

    auto tg_xxyz_xxzzz = pbuffer.data(idx_gh + 93);

    auto tg_xxyz_xzzzz = pbuffer.data(idx_gh + 98);

    auto tg_xxyz_yyyyz = pbuffer.data(idx_gh + 100);

    auto tg_xxyz_yyyzz = pbuffer.data(idx_gh + 101);

    auto tg_xxyz_yyzzz = pbuffer.data(idx_gh + 102);

    auto tg_xxyz_yzzzz = pbuffer.data(idx_gh + 103);

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

    auto tg_xyyy_xxxyy = pbuffer.data(idx_gh + 129);

    auto tg_xyyy_xxxyz = pbuffer.data(idx_gh + 130);

    auto tg_xyyy_xxyyy = pbuffer.data(idx_gh + 132);

    auto tg_xyyy_xxyyz = pbuffer.data(idx_gh + 133);

    auto tg_xyyy_xxyzz = pbuffer.data(idx_gh + 134);

    auto tg_xyyy_xyyyy = pbuffer.data(idx_gh + 136);

    auto tg_xyyy_xyyyz = pbuffer.data(idx_gh + 137);

    auto tg_xyyy_xyyzz = pbuffer.data(idx_gh + 138);

    auto tg_xyyy_xyzzz = pbuffer.data(idx_gh + 139);

    auto tg_xyyy_yyyyy = pbuffer.data(idx_gh + 141);

    auto tg_xyyy_yyyyz = pbuffer.data(idx_gh + 142);

    auto tg_xyyy_yyyzz = pbuffer.data(idx_gh + 143);

    auto tg_xyyy_yyzzz = pbuffer.data(idx_gh + 144);

    auto tg_xyyy_yzzzz = pbuffer.data(idx_gh + 145);

    auto tg_xyyy_zzzzz = pbuffer.data(idx_gh + 146);

    auto tg_xyyz_yyyyz = pbuffer.data(idx_gh + 163);

    auto tg_xyyz_yyyzz = pbuffer.data(idx_gh + 164);

    auto tg_xyyz_yyzzz = pbuffer.data(idx_gh + 165);

    auto tg_xyyz_yzzzz = pbuffer.data(idx_gh + 166);

    auto tg_xyyz_zzzzz = pbuffer.data(idx_gh + 167);

    auto tg_xyzz_yyyyy = pbuffer.data(idx_gh + 183);

    auto tg_xyzz_yyyyz = pbuffer.data(idx_gh + 184);

    auto tg_xyzz_yyyzz = pbuffer.data(idx_gh + 185);

    auto tg_xyzz_yyzzz = pbuffer.data(idx_gh + 186);

    auto tg_xyzz_yzzzz = pbuffer.data(idx_gh + 187);

    auto tg_xzzz_xxxxx = pbuffer.data(idx_gh + 189);

    auto tg_xzzz_xxxxz = pbuffer.data(idx_gh + 191);

    auto tg_xzzz_xxxyz = pbuffer.data(idx_gh + 193);

    auto tg_xzzz_xxxzz = pbuffer.data(idx_gh + 194);

    auto tg_xzzz_xxyyz = pbuffer.data(idx_gh + 196);

    auto tg_xzzz_xxyzz = pbuffer.data(idx_gh + 197);

    auto tg_xzzz_xxzzz = pbuffer.data(idx_gh + 198);

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

    // Set up components of targeted buffer : HH

    auto tg_xxxxx_xxxxx = pbuffer.data(idx_hh);

    auto tg_xxxxx_xxxxy = pbuffer.data(idx_hh + 1);

    auto tg_xxxxx_xxxxz = pbuffer.data(idx_hh + 2);

    auto tg_xxxxx_xxxyy = pbuffer.data(idx_hh + 3);

    auto tg_xxxxx_xxxyz = pbuffer.data(idx_hh + 4);

    auto tg_xxxxx_xxxzz = pbuffer.data(idx_hh + 5);

    auto tg_xxxxx_xxyyy = pbuffer.data(idx_hh + 6);

    auto tg_xxxxx_xxyyz = pbuffer.data(idx_hh + 7);

    auto tg_xxxxx_xxyzz = pbuffer.data(idx_hh + 8);

    auto tg_xxxxx_xxzzz = pbuffer.data(idx_hh + 9);

    auto tg_xxxxx_xyyyy = pbuffer.data(idx_hh + 10);

    auto tg_xxxxx_xyyyz = pbuffer.data(idx_hh + 11);

    auto tg_xxxxx_xyyzz = pbuffer.data(idx_hh + 12);

    auto tg_xxxxx_xyzzz = pbuffer.data(idx_hh + 13);

    auto tg_xxxxx_xzzzz = pbuffer.data(idx_hh + 14);

    auto tg_xxxxx_yyyyy = pbuffer.data(idx_hh + 15);

    auto tg_xxxxx_yyyyz = pbuffer.data(idx_hh + 16);

    auto tg_xxxxx_yyyzz = pbuffer.data(idx_hh + 17);

    auto tg_xxxxx_yyzzz = pbuffer.data(idx_hh + 18);

    auto tg_xxxxx_yzzzz = pbuffer.data(idx_hh + 19);

    auto tg_xxxxx_zzzzz = pbuffer.data(idx_hh + 20);

    auto tg_xxxxy_xxxxx = pbuffer.data(idx_hh + 21);

    auto tg_xxxxy_xxxxy = pbuffer.data(idx_hh + 22);

    auto tg_xxxxy_xxxxz = pbuffer.data(idx_hh + 23);

    auto tg_xxxxy_xxxyy = pbuffer.data(idx_hh + 24);

    auto tg_xxxxy_xxxyz = pbuffer.data(idx_hh + 25);

    auto tg_xxxxy_xxxzz = pbuffer.data(idx_hh + 26);

    auto tg_xxxxy_xxyyy = pbuffer.data(idx_hh + 27);

    auto tg_xxxxy_xxyyz = pbuffer.data(idx_hh + 28);

    auto tg_xxxxy_xxyzz = pbuffer.data(idx_hh + 29);

    auto tg_xxxxy_xxzzz = pbuffer.data(idx_hh + 30);

    auto tg_xxxxy_xyyyy = pbuffer.data(idx_hh + 31);

    auto tg_xxxxy_xyyyz = pbuffer.data(idx_hh + 32);

    auto tg_xxxxy_xyyzz = pbuffer.data(idx_hh + 33);

    auto tg_xxxxy_xyzzz = pbuffer.data(idx_hh + 34);

    auto tg_xxxxy_xzzzz = pbuffer.data(idx_hh + 35);

    auto tg_xxxxy_yyyyy = pbuffer.data(idx_hh + 36);

    auto tg_xxxxy_yyyyz = pbuffer.data(idx_hh + 37);

    auto tg_xxxxy_yyyzz = pbuffer.data(idx_hh + 38);

    auto tg_xxxxy_yyzzz = pbuffer.data(idx_hh + 39);

    auto tg_xxxxy_yzzzz = pbuffer.data(idx_hh + 40);

    auto tg_xxxxy_zzzzz = pbuffer.data(idx_hh + 41);

    auto tg_xxxxz_xxxxx = pbuffer.data(idx_hh + 42);

    auto tg_xxxxz_xxxxy = pbuffer.data(idx_hh + 43);

    auto tg_xxxxz_xxxxz = pbuffer.data(idx_hh + 44);

    auto tg_xxxxz_xxxyy = pbuffer.data(idx_hh + 45);

    auto tg_xxxxz_xxxyz = pbuffer.data(idx_hh + 46);

    auto tg_xxxxz_xxxzz = pbuffer.data(idx_hh + 47);

    auto tg_xxxxz_xxyyy = pbuffer.data(idx_hh + 48);

    auto tg_xxxxz_xxyyz = pbuffer.data(idx_hh + 49);

    auto tg_xxxxz_xxyzz = pbuffer.data(idx_hh + 50);

    auto tg_xxxxz_xxzzz = pbuffer.data(idx_hh + 51);

    auto tg_xxxxz_xyyyy = pbuffer.data(idx_hh + 52);

    auto tg_xxxxz_xyyyz = pbuffer.data(idx_hh + 53);

    auto tg_xxxxz_xyyzz = pbuffer.data(idx_hh + 54);

    auto tg_xxxxz_xyzzz = pbuffer.data(idx_hh + 55);

    auto tg_xxxxz_xzzzz = pbuffer.data(idx_hh + 56);

    auto tg_xxxxz_yyyyy = pbuffer.data(idx_hh + 57);

    auto tg_xxxxz_yyyyz = pbuffer.data(idx_hh + 58);

    auto tg_xxxxz_yyyzz = pbuffer.data(idx_hh + 59);

    auto tg_xxxxz_yyzzz = pbuffer.data(idx_hh + 60);

    auto tg_xxxxz_yzzzz = pbuffer.data(idx_hh + 61);

    auto tg_xxxxz_zzzzz = pbuffer.data(idx_hh + 62);

    auto tg_xxxyy_xxxxx = pbuffer.data(idx_hh + 63);

    auto tg_xxxyy_xxxxy = pbuffer.data(idx_hh + 64);

    auto tg_xxxyy_xxxxz = pbuffer.data(idx_hh + 65);

    auto tg_xxxyy_xxxyy = pbuffer.data(idx_hh + 66);

    auto tg_xxxyy_xxxyz = pbuffer.data(idx_hh + 67);

    auto tg_xxxyy_xxxzz = pbuffer.data(idx_hh + 68);

    auto tg_xxxyy_xxyyy = pbuffer.data(idx_hh + 69);

    auto tg_xxxyy_xxyyz = pbuffer.data(idx_hh + 70);

    auto tg_xxxyy_xxyzz = pbuffer.data(idx_hh + 71);

    auto tg_xxxyy_xxzzz = pbuffer.data(idx_hh + 72);

    auto tg_xxxyy_xyyyy = pbuffer.data(idx_hh + 73);

    auto tg_xxxyy_xyyyz = pbuffer.data(idx_hh + 74);

    auto tg_xxxyy_xyyzz = pbuffer.data(idx_hh + 75);

    auto tg_xxxyy_xyzzz = pbuffer.data(idx_hh + 76);

    auto tg_xxxyy_xzzzz = pbuffer.data(idx_hh + 77);

    auto tg_xxxyy_yyyyy = pbuffer.data(idx_hh + 78);

    auto tg_xxxyy_yyyyz = pbuffer.data(idx_hh + 79);

    auto tg_xxxyy_yyyzz = pbuffer.data(idx_hh + 80);

    auto tg_xxxyy_yyzzz = pbuffer.data(idx_hh + 81);

    auto tg_xxxyy_yzzzz = pbuffer.data(idx_hh + 82);

    auto tg_xxxyy_zzzzz = pbuffer.data(idx_hh + 83);

    auto tg_xxxyz_xxxxx = pbuffer.data(idx_hh + 84);

    auto tg_xxxyz_xxxxy = pbuffer.data(idx_hh + 85);

    auto tg_xxxyz_xxxxz = pbuffer.data(idx_hh + 86);

    auto tg_xxxyz_xxxyy = pbuffer.data(idx_hh + 87);

    auto tg_xxxyz_xxxyz = pbuffer.data(idx_hh + 88);

    auto tg_xxxyz_xxxzz = pbuffer.data(idx_hh + 89);

    auto tg_xxxyz_xxyyy = pbuffer.data(idx_hh + 90);

    auto tg_xxxyz_xxyyz = pbuffer.data(idx_hh + 91);

    auto tg_xxxyz_xxyzz = pbuffer.data(idx_hh + 92);

    auto tg_xxxyz_xxzzz = pbuffer.data(idx_hh + 93);

    auto tg_xxxyz_xyyyy = pbuffer.data(idx_hh + 94);

    auto tg_xxxyz_xyyyz = pbuffer.data(idx_hh + 95);

    auto tg_xxxyz_xyyzz = pbuffer.data(idx_hh + 96);

    auto tg_xxxyz_xyzzz = pbuffer.data(idx_hh + 97);

    auto tg_xxxyz_xzzzz = pbuffer.data(idx_hh + 98);

    auto tg_xxxyz_yyyyy = pbuffer.data(idx_hh + 99);

    auto tg_xxxyz_yyyyz = pbuffer.data(idx_hh + 100);

    auto tg_xxxyz_yyyzz = pbuffer.data(idx_hh + 101);

    auto tg_xxxyz_yyzzz = pbuffer.data(idx_hh + 102);

    auto tg_xxxyz_yzzzz = pbuffer.data(idx_hh + 103);

    auto tg_xxxyz_zzzzz = pbuffer.data(idx_hh + 104);

    auto tg_xxxzz_xxxxx = pbuffer.data(idx_hh + 105);

    auto tg_xxxzz_xxxxy = pbuffer.data(idx_hh + 106);

    auto tg_xxxzz_xxxxz = pbuffer.data(idx_hh + 107);

    auto tg_xxxzz_xxxyy = pbuffer.data(idx_hh + 108);

    auto tg_xxxzz_xxxyz = pbuffer.data(idx_hh + 109);

    auto tg_xxxzz_xxxzz = pbuffer.data(idx_hh + 110);

    auto tg_xxxzz_xxyyy = pbuffer.data(idx_hh + 111);

    auto tg_xxxzz_xxyyz = pbuffer.data(idx_hh + 112);

    auto tg_xxxzz_xxyzz = pbuffer.data(idx_hh + 113);

    auto tg_xxxzz_xxzzz = pbuffer.data(idx_hh + 114);

    auto tg_xxxzz_xyyyy = pbuffer.data(idx_hh + 115);

    auto tg_xxxzz_xyyyz = pbuffer.data(idx_hh + 116);

    auto tg_xxxzz_xyyzz = pbuffer.data(idx_hh + 117);

    auto tg_xxxzz_xyzzz = pbuffer.data(idx_hh + 118);

    auto tg_xxxzz_xzzzz = pbuffer.data(idx_hh + 119);

    auto tg_xxxzz_yyyyy = pbuffer.data(idx_hh + 120);

    auto tg_xxxzz_yyyyz = pbuffer.data(idx_hh + 121);

    auto tg_xxxzz_yyyzz = pbuffer.data(idx_hh + 122);

    auto tg_xxxzz_yyzzz = pbuffer.data(idx_hh + 123);

    auto tg_xxxzz_yzzzz = pbuffer.data(idx_hh + 124);

    auto tg_xxxzz_zzzzz = pbuffer.data(idx_hh + 125);

    auto tg_xxyyy_xxxxx = pbuffer.data(idx_hh + 126);

    auto tg_xxyyy_xxxxy = pbuffer.data(idx_hh + 127);

    auto tg_xxyyy_xxxxz = pbuffer.data(idx_hh + 128);

    auto tg_xxyyy_xxxyy = pbuffer.data(idx_hh + 129);

    auto tg_xxyyy_xxxyz = pbuffer.data(idx_hh + 130);

    auto tg_xxyyy_xxxzz = pbuffer.data(idx_hh + 131);

    auto tg_xxyyy_xxyyy = pbuffer.data(idx_hh + 132);

    auto tg_xxyyy_xxyyz = pbuffer.data(idx_hh + 133);

    auto tg_xxyyy_xxyzz = pbuffer.data(idx_hh + 134);

    auto tg_xxyyy_xxzzz = pbuffer.data(idx_hh + 135);

    auto tg_xxyyy_xyyyy = pbuffer.data(idx_hh + 136);

    auto tg_xxyyy_xyyyz = pbuffer.data(idx_hh + 137);

    auto tg_xxyyy_xyyzz = pbuffer.data(idx_hh + 138);

    auto tg_xxyyy_xyzzz = pbuffer.data(idx_hh + 139);

    auto tg_xxyyy_xzzzz = pbuffer.data(idx_hh + 140);

    auto tg_xxyyy_yyyyy = pbuffer.data(idx_hh + 141);

    auto tg_xxyyy_yyyyz = pbuffer.data(idx_hh + 142);

    auto tg_xxyyy_yyyzz = pbuffer.data(idx_hh + 143);

    auto tg_xxyyy_yyzzz = pbuffer.data(idx_hh + 144);

    auto tg_xxyyy_yzzzz = pbuffer.data(idx_hh + 145);

    auto tg_xxyyy_zzzzz = pbuffer.data(idx_hh + 146);

    auto tg_xxyyz_xxxxx = pbuffer.data(idx_hh + 147);

    auto tg_xxyyz_xxxxy = pbuffer.data(idx_hh + 148);

    auto tg_xxyyz_xxxxz = pbuffer.data(idx_hh + 149);

    auto tg_xxyyz_xxxyy = pbuffer.data(idx_hh + 150);

    auto tg_xxyyz_xxxyz = pbuffer.data(idx_hh + 151);

    auto tg_xxyyz_xxxzz = pbuffer.data(idx_hh + 152);

    auto tg_xxyyz_xxyyy = pbuffer.data(idx_hh + 153);

    auto tg_xxyyz_xxyyz = pbuffer.data(idx_hh + 154);

    auto tg_xxyyz_xxyzz = pbuffer.data(idx_hh + 155);

    auto tg_xxyyz_xxzzz = pbuffer.data(idx_hh + 156);

    auto tg_xxyyz_xyyyy = pbuffer.data(idx_hh + 157);

    auto tg_xxyyz_xyyyz = pbuffer.data(idx_hh + 158);

    auto tg_xxyyz_xyyzz = pbuffer.data(idx_hh + 159);

    auto tg_xxyyz_xyzzz = pbuffer.data(idx_hh + 160);

    auto tg_xxyyz_xzzzz = pbuffer.data(idx_hh + 161);

    auto tg_xxyyz_yyyyy = pbuffer.data(idx_hh + 162);

    auto tg_xxyyz_yyyyz = pbuffer.data(idx_hh + 163);

    auto tg_xxyyz_yyyzz = pbuffer.data(idx_hh + 164);

    auto tg_xxyyz_yyzzz = pbuffer.data(idx_hh + 165);

    auto tg_xxyyz_yzzzz = pbuffer.data(idx_hh + 166);

    auto tg_xxyyz_zzzzz = pbuffer.data(idx_hh + 167);

    auto tg_xxyzz_xxxxx = pbuffer.data(idx_hh + 168);

    auto tg_xxyzz_xxxxy = pbuffer.data(idx_hh + 169);

    auto tg_xxyzz_xxxxz = pbuffer.data(idx_hh + 170);

    auto tg_xxyzz_xxxyy = pbuffer.data(idx_hh + 171);

    auto tg_xxyzz_xxxyz = pbuffer.data(idx_hh + 172);

    auto tg_xxyzz_xxxzz = pbuffer.data(idx_hh + 173);

    auto tg_xxyzz_xxyyy = pbuffer.data(idx_hh + 174);

    auto tg_xxyzz_xxyyz = pbuffer.data(idx_hh + 175);

    auto tg_xxyzz_xxyzz = pbuffer.data(idx_hh + 176);

    auto tg_xxyzz_xxzzz = pbuffer.data(idx_hh + 177);

    auto tg_xxyzz_xyyyy = pbuffer.data(idx_hh + 178);

    auto tg_xxyzz_xyyyz = pbuffer.data(idx_hh + 179);

    auto tg_xxyzz_xyyzz = pbuffer.data(idx_hh + 180);

    auto tg_xxyzz_xyzzz = pbuffer.data(idx_hh + 181);

    auto tg_xxyzz_xzzzz = pbuffer.data(idx_hh + 182);

    auto tg_xxyzz_yyyyy = pbuffer.data(idx_hh + 183);

    auto tg_xxyzz_yyyyz = pbuffer.data(idx_hh + 184);

    auto tg_xxyzz_yyyzz = pbuffer.data(idx_hh + 185);

    auto tg_xxyzz_yyzzz = pbuffer.data(idx_hh + 186);

    auto tg_xxyzz_yzzzz = pbuffer.data(idx_hh + 187);

    auto tg_xxyzz_zzzzz = pbuffer.data(idx_hh + 188);

    auto tg_xxzzz_xxxxx = pbuffer.data(idx_hh + 189);

    auto tg_xxzzz_xxxxy = pbuffer.data(idx_hh + 190);

    auto tg_xxzzz_xxxxz = pbuffer.data(idx_hh + 191);

    auto tg_xxzzz_xxxyy = pbuffer.data(idx_hh + 192);

    auto tg_xxzzz_xxxyz = pbuffer.data(idx_hh + 193);

    auto tg_xxzzz_xxxzz = pbuffer.data(idx_hh + 194);

    auto tg_xxzzz_xxyyy = pbuffer.data(idx_hh + 195);

    auto tg_xxzzz_xxyyz = pbuffer.data(idx_hh + 196);

    auto tg_xxzzz_xxyzz = pbuffer.data(idx_hh + 197);

    auto tg_xxzzz_xxzzz = pbuffer.data(idx_hh + 198);

    auto tg_xxzzz_xyyyy = pbuffer.data(idx_hh + 199);

    auto tg_xxzzz_xyyyz = pbuffer.data(idx_hh + 200);

    auto tg_xxzzz_xyyzz = pbuffer.data(idx_hh + 201);

    auto tg_xxzzz_xyzzz = pbuffer.data(idx_hh + 202);

    auto tg_xxzzz_xzzzz = pbuffer.data(idx_hh + 203);

    auto tg_xxzzz_yyyyy = pbuffer.data(idx_hh + 204);

    auto tg_xxzzz_yyyyz = pbuffer.data(idx_hh + 205);

    auto tg_xxzzz_yyyzz = pbuffer.data(idx_hh + 206);

    auto tg_xxzzz_yyzzz = pbuffer.data(idx_hh + 207);

    auto tg_xxzzz_yzzzz = pbuffer.data(idx_hh + 208);

    auto tg_xxzzz_zzzzz = pbuffer.data(idx_hh + 209);

    auto tg_xyyyy_xxxxx = pbuffer.data(idx_hh + 210);

    auto tg_xyyyy_xxxxy = pbuffer.data(idx_hh + 211);

    auto tg_xyyyy_xxxxz = pbuffer.data(idx_hh + 212);

    auto tg_xyyyy_xxxyy = pbuffer.data(idx_hh + 213);

    auto tg_xyyyy_xxxyz = pbuffer.data(idx_hh + 214);

    auto tg_xyyyy_xxxzz = pbuffer.data(idx_hh + 215);

    auto tg_xyyyy_xxyyy = pbuffer.data(idx_hh + 216);

    auto tg_xyyyy_xxyyz = pbuffer.data(idx_hh + 217);

    auto tg_xyyyy_xxyzz = pbuffer.data(idx_hh + 218);

    auto tg_xyyyy_xxzzz = pbuffer.data(idx_hh + 219);

    auto tg_xyyyy_xyyyy = pbuffer.data(idx_hh + 220);

    auto tg_xyyyy_xyyyz = pbuffer.data(idx_hh + 221);

    auto tg_xyyyy_xyyzz = pbuffer.data(idx_hh + 222);

    auto tg_xyyyy_xyzzz = pbuffer.data(idx_hh + 223);

    auto tg_xyyyy_xzzzz = pbuffer.data(idx_hh + 224);

    auto tg_xyyyy_yyyyy = pbuffer.data(idx_hh + 225);

    auto tg_xyyyy_yyyyz = pbuffer.data(idx_hh + 226);

    auto tg_xyyyy_yyyzz = pbuffer.data(idx_hh + 227);

    auto tg_xyyyy_yyzzz = pbuffer.data(idx_hh + 228);

    auto tg_xyyyy_yzzzz = pbuffer.data(idx_hh + 229);

    auto tg_xyyyy_zzzzz = pbuffer.data(idx_hh + 230);

    auto tg_xyyyz_xxxxx = pbuffer.data(idx_hh + 231);

    auto tg_xyyyz_xxxxy = pbuffer.data(idx_hh + 232);

    auto tg_xyyyz_xxxxz = pbuffer.data(idx_hh + 233);

    auto tg_xyyyz_xxxyy = pbuffer.data(idx_hh + 234);

    auto tg_xyyyz_xxxyz = pbuffer.data(idx_hh + 235);

    auto tg_xyyyz_xxxzz = pbuffer.data(idx_hh + 236);

    auto tg_xyyyz_xxyyy = pbuffer.data(idx_hh + 237);

    auto tg_xyyyz_xxyyz = pbuffer.data(idx_hh + 238);

    auto tg_xyyyz_xxyzz = pbuffer.data(idx_hh + 239);

    auto tg_xyyyz_xxzzz = pbuffer.data(idx_hh + 240);

    auto tg_xyyyz_xyyyy = pbuffer.data(idx_hh + 241);

    auto tg_xyyyz_xyyyz = pbuffer.data(idx_hh + 242);

    auto tg_xyyyz_xyyzz = pbuffer.data(idx_hh + 243);

    auto tg_xyyyz_xyzzz = pbuffer.data(idx_hh + 244);

    auto tg_xyyyz_xzzzz = pbuffer.data(idx_hh + 245);

    auto tg_xyyyz_yyyyy = pbuffer.data(idx_hh + 246);

    auto tg_xyyyz_yyyyz = pbuffer.data(idx_hh + 247);

    auto tg_xyyyz_yyyzz = pbuffer.data(idx_hh + 248);

    auto tg_xyyyz_yyzzz = pbuffer.data(idx_hh + 249);

    auto tg_xyyyz_yzzzz = pbuffer.data(idx_hh + 250);

    auto tg_xyyyz_zzzzz = pbuffer.data(idx_hh + 251);

    auto tg_xyyzz_xxxxx = pbuffer.data(idx_hh + 252);

    auto tg_xyyzz_xxxxy = pbuffer.data(idx_hh + 253);

    auto tg_xyyzz_xxxxz = pbuffer.data(idx_hh + 254);

    auto tg_xyyzz_xxxyy = pbuffer.data(idx_hh + 255);

    auto tg_xyyzz_xxxyz = pbuffer.data(idx_hh + 256);

    auto tg_xyyzz_xxxzz = pbuffer.data(idx_hh + 257);

    auto tg_xyyzz_xxyyy = pbuffer.data(idx_hh + 258);

    auto tg_xyyzz_xxyyz = pbuffer.data(idx_hh + 259);

    auto tg_xyyzz_xxyzz = pbuffer.data(idx_hh + 260);

    auto tg_xyyzz_xxzzz = pbuffer.data(idx_hh + 261);

    auto tg_xyyzz_xyyyy = pbuffer.data(idx_hh + 262);

    auto tg_xyyzz_xyyyz = pbuffer.data(idx_hh + 263);

    auto tg_xyyzz_xyyzz = pbuffer.data(idx_hh + 264);

    auto tg_xyyzz_xyzzz = pbuffer.data(idx_hh + 265);

    auto tg_xyyzz_xzzzz = pbuffer.data(idx_hh + 266);

    auto tg_xyyzz_yyyyy = pbuffer.data(idx_hh + 267);

    auto tg_xyyzz_yyyyz = pbuffer.data(idx_hh + 268);

    auto tg_xyyzz_yyyzz = pbuffer.data(idx_hh + 269);

    auto tg_xyyzz_yyzzz = pbuffer.data(idx_hh + 270);

    auto tg_xyyzz_yzzzz = pbuffer.data(idx_hh + 271);

    auto tg_xyyzz_zzzzz = pbuffer.data(idx_hh + 272);

    auto tg_xyzzz_xxxxx = pbuffer.data(idx_hh + 273);

    auto tg_xyzzz_xxxxy = pbuffer.data(idx_hh + 274);

    auto tg_xyzzz_xxxxz = pbuffer.data(idx_hh + 275);

    auto tg_xyzzz_xxxyy = pbuffer.data(idx_hh + 276);

    auto tg_xyzzz_xxxyz = pbuffer.data(idx_hh + 277);

    auto tg_xyzzz_xxxzz = pbuffer.data(idx_hh + 278);

    auto tg_xyzzz_xxyyy = pbuffer.data(idx_hh + 279);

    auto tg_xyzzz_xxyyz = pbuffer.data(idx_hh + 280);

    auto tg_xyzzz_xxyzz = pbuffer.data(idx_hh + 281);

    auto tg_xyzzz_xxzzz = pbuffer.data(idx_hh + 282);

    auto tg_xyzzz_xyyyy = pbuffer.data(idx_hh + 283);

    auto tg_xyzzz_xyyyz = pbuffer.data(idx_hh + 284);

    auto tg_xyzzz_xyyzz = pbuffer.data(idx_hh + 285);

    auto tg_xyzzz_xyzzz = pbuffer.data(idx_hh + 286);

    auto tg_xyzzz_xzzzz = pbuffer.data(idx_hh + 287);

    auto tg_xyzzz_yyyyy = pbuffer.data(idx_hh + 288);

    auto tg_xyzzz_yyyyz = pbuffer.data(idx_hh + 289);

    auto tg_xyzzz_yyyzz = pbuffer.data(idx_hh + 290);

    auto tg_xyzzz_yyzzz = pbuffer.data(idx_hh + 291);

    auto tg_xyzzz_yzzzz = pbuffer.data(idx_hh + 292);

    auto tg_xyzzz_zzzzz = pbuffer.data(idx_hh + 293);

    auto tg_xzzzz_xxxxx = pbuffer.data(idx_hh + 294);

    auto tg_xzzzz_xxxxy = pbuffer.data(idx_hh + 295);

    auto tg_xzzzz_xxxxz = pbuffer.data(idx_hh + 296);

    auto tg_xzzzz_xxxyy = pbuffer.data(idx_hh + 297);

    auto tg_xzzzz_xxxyz = pbuffer.data(idx_hh + 298);

    auto tg_xzzzz_xxxzz = pbuffer.data(idx_hh + 299);

    auto tg_xzzzz_xxyyy = pbuffer.data(idx_hh + 300);

    auto tg_xzzzz_xxyyz = pbuffer.data(idx_hh + 301);

    auto tg_xzzzz_xxyzz = pbuffer.data(idx_hh + 302);

    auto tg_xzzzz_xxzzz = pbuffer.data(idx_hh + 303);

    auto tg_xzzzz_xyyyy = pbuffer.data(idx_hh + 304);

    auto tg_xzzzz_xyyyz = pbuffer.data(idx_hh + 305);

    auto tg_xzzzz_xyyzz = pbuffer.data(idx_hh + 306);

    auto tg_xzzzz_xyzzz = pbuffer.data(idx_hh + 307);

    auto tg_xzzzz_xzzzz = pbuffer.data(idx_hh + 308);

    auto tg_xzzzz_yyyyy = pbuffer.data(idx_hh + 309);

    auto tg_xzzzz_yyyyz = pbuffer.data(idx_hh + 310);

    auto tg_xzzzz_yyyzz = pbuffer.data(idx_hh + 311);

    auto tg_xzzzz_yyzzz = pbuffer.data(idx_hh + 312);

    auto tg_xzzzz_yzzzz = pbuffer.data(idx_hh + 313);

    auto tg_xzzzz_zzzzz = pbuffer.data(idx_hh + 314);

    auto tg_yyyyy_xxxxx = pbuffer.data(idx_hh + 315);

    auto tg_yyyyy_xxxxy = pbuffer.data(idx_hh + 316);

    auto tg_yyyyy_xxxxz = pbuffer.data(idx_hh + 317);

    auto tg_yyyyy_xxxyy = pbuffer.data(idx_hh + 318);

    auto tg_yyyyy_xxxyz = pbuffer.data(idx_hh + 319);

    auto tg_yyyyy_xxxzz = pbuffer.data(idx_hh + 320);

    auto tg_yyyyy_xxyyy = pbuffer.data(idx_hh + 321);

    auto tg_yyyyy_xxyyz = pbuffer.data(idx_hh + 322);

    auto tg_yyyyy_xxyzz = pbuffer.data(idx_hh + 323);

    auto tg_yyyyy_xxzzz = pbuffer.data(idx_hh + 324);

    auto tg_yyyyy_xyyyy = pbuffer.data(idx_hh + 325);

    auto tg_yyyyy_xyyyz = pbuffer.data(idx_hh + 326);

    auto tg_yyyyy_xyyzz = pbuffer.data(idx_hh + 327);

    auto tg_yyyyy_xyzzz = pbuffer.data(idx_hh + 328);

    auto tg_yyyyy_xzzzz = pbuffer.data(idx_hh + 329);

    auto tg_yyyyy_yyyyy = pbuffer.data(idx_hh + 330);

    auto tg_yyyyy_yyyyz = pbuffer.data(idx_hh + 331);

    auto tg_yyyyy_yyyzz = pbuffer.data(idx_hh + 332);

    auto tg_yyyyy_yyzzz = pbuffer.data(idx_hh + 333);

    auto tg_yyyyy_yzzzz = pbuffer.data(idx_hh + 334);

    auto tg_yyyyy_zzzzz = pbuffer.data(idx_hh + 335);

    auto tg_yyyyz_xxxxx = pbuffer.data(idx_hh + 336);

    auto tg_yyyyz_xxxxy = pbuffer.data(idx_hh + 337);

    auto tg_yyyyz_xxxxz = pbuffer.data(idx_hh + 338);

    auto tg_yyyyz_xxxyy = pbuffer.data(idx_hh + 339);

    auto tg_yyyyz_xxxyz = pbuffer.data(idx_hh + 340);

    auto tg_yyyyz_xxxzz = pbuffer.data(idx_hh + 341);

    auto tg_yyyyz_xxyyy = pbuffer.data(idx_hh + 342);

    auto tg_yyyyz_xxyyz = pbuffer.data(idx_hh + 343);

    auto tg_yyyyz_xxyzz = pbuffer.data(idx_hh + 344);

    auto tg_yyyyz_xxzzz = pbuffer.data(idx_hh + 345);

    auto tg_yyyyz_xyyyy = pbuffer.data(idx_hh + 346);

    auto tg_yyyyz_xyyyz = pbuffer.data(idx_hh + 347);

    auto tg_yyyyz_xyyzz = pbuffer.data(idx_hh + 348);

    auto tg_yyyyz_xyzzz = pbuffer.data(idx_hh + 349);

    auto tg_yyyyz_xzzzz = pbuffer.data(idx_hh + 350);

    auto tg_yyyyz_yyyyy = pbuffer.data(idx_hh + 351);

    auto tg_yyyyz_yyyyz = pbuffer.data(idx_hh + 352);

    auto tg_yyyyz_yyyzz = pbuffer.data(idx_hh + 353);

    auto tg_yyyyz_yyzzz = pbuffer.data(idx_hh + 354);

    auto tg_yyyyz_yzzzz = pbuffer.data(idx_hh + 355);

    auto tg_yyyyz_zzzzz = pbuffer.data(idx_hh + 356);

    auto tg_yyyzz_xxxxx = pbuffer.data(idx_hh + 357);

    auto tg_yyyzz_xxxxy = pbuffer.data(idx_hh + 358);

    auto tg_yyyzz_xxxxz = pbuffer.data(idx_hh + 359);

    auto tg_yyyzz_xxxyy = pbuffer.data(idx_hh + 360);

    auto tg_yyyzz_xxxyz = pbuffer.data(idx_hh + 361);

    auto tg_yyyzz_xxxzz = pbuffer.data(idx_hh + 362);

    auto tg_yyyzz_xxyyy = pbuffer.data(idx_hh + 363);

    auto tg_yyyzz_xxyyz = pbuffer.data(idx_hh + 364);

    auto tg_yyyzz_xxyzz = pbuffer.data(idx_hh + 365);

    auto tg_yyyzz_xxzzz = pbuffer.data(idx_hh + 366);

    auto tg_yyyzz_xyyyy = pbuffer.data(idx_hh + 367);

    auto tg_yyyzz_xyyyz = pbuffer.data(idx_hh + 368);

    auto tg_yyyzz_xyyzz = pbuffer.data(idx_hh + 369);

    auto tg_yyyzz_xyzzz = pbuffer.data(idx_hh + 370);

    auto tg_yyyzz_xzzzz = pbuffer.data(idx_hh + 371);

    auto tg_yyyzz_yyyyy = pbuffer.data(idx_hh + 372);

    auto tg_yyyzz_yyyyz = pbuffer.data(idx_hh + 373);

    auto tg_yyyzz_yyyzz = pbuffer.data(idx_hh + 374);

    auto tg_yyyzz_yyzzz = pbuffer.data(idx_hh + 375);

    auto tg_yyyzz_yzzzz = pbuffer.data(idx_hh + 376);

    auto tg_yyyzz_zzzzz = pbuffer.data(idx_hh + 377);

    auto tg_yyzzz_xxxxx = pbuffer.data(idx_hh + 378);

    auto tg_yyzzz_xxxxy = pbuffer.data(idx_hh + 379);

    auto tg_yyzzz_xxxxz = pbuffer.data(idx_hh + 380);

    auto tg_yyzzz_xxxyy = pbuffer.data(idx_hh + 381);

    auto tg_yyzzz_xxxyz = pbuffer.data(idx_hh + 382);

    auto tg_yyzzz_xxxzz = pbuffer.data(idx_hh + 383);

    auto tg_yyzzz_xxyyy = pbuffer.data(idx_hh + 384);

    auto tg_yyzzz_xxyyz = pbuffer.data(idx_hh + 385);

    auto tg_yyzzz_xxyzz = pbuffer.data(idx_hh + 386);

    auto tg_yyzzz_xxzzz = pbuffer.data(idx_hh + 387);

    auto tg_yyzzz_xyyyy = pbuffer.data(idx_hh + 388);

    auto tg_yyzzz_xyyyz = pbuffer.data(idx_hh + 389);

    auto tg_yyzzz_xyyzz = pbuffer.data(idx_hh + 390);

    auto tg_yyzzz_xyzzz = pbuffer.data(idx_hh + 391);

    auto tg_yyzzz_xzzzz = pbuffer.data(idx_hh + 392);

    auto tg_yyzzz_yyyyy = pbuffer.data(idx_hh + 393);

    auto tg_yyzzz_yyyyz = pbuffer.data(idx_hh + 394);

    auto tg_yyzzz_yyyzz = pbuffer.data(idx_hh + 395);

    auto tg_yyzzz_yyzzz = pbuffer.data(idx_hh + 396);

    auto tg_yyzzz_yzzzz = pbuffer.data(idx_hh + 397);

    auto tg_yyzzz_zzzzz = pbuffer.data(idx_hh + 398);

    auto tg_yzzzz_xxxxx = pbuffer.data(idx_hh + 399);

    auto tg_yzzzz_xxxxy = pbuffer.data(idx_hh + 400);

    auto tg_yzzzz_xxxxz = pbuffer.data(idx_hh + 401);

    auto tg_yzzzz_xxxyy = pbuffer.data(idx_hh + 402);

    auto tg_yzzzz_xxxyz = pbuffer.data(idx_hh + 403);

    auto tg_yzzzz_xxxzz = pbuffer.data(idx_hh + 404);

    auto tg_yzzzz_xxyyy = pbuffer.data(idx_hh + 405);

    auto tg_yzzzz_xxyyz = pbuffer.data(idx_hh + 406);

    auto tg_yzzzz_xxyzz = pbuffer.data(idx_hh + 407);

    auto tg_yzzzz_xxzzz = pbuffer.data(idx_hh + 408);

    auto tg_yzzzz_xyyyy = pbuffer.data(idx_hh + 409);

    auto tg_yzzzz_xyyyz = pbuffer.data(idx_hh + 410);

    auto tg_yzzzz_xyyzz = pbuffer.data(idx_hh + 411);

    auto tg_yzzzz_xyzzz = pbuffer.data(idx_hh + 412);

    auto tg_yzzzz_xzzzz = pbuffer.data(idx_hh + 413);

    auto tg_yzzzz_yyyyy = pbuffer.data(idx_hh + 414);

    auto tg_yzzzz_yyyyz = pbuffer.data(idx_hh + 415);

    auto tg_yzzzz_yyyzz = pbuffer.data(idx_hh + 416);

    auto tg_yzzzz_yyzzz = pbuffer.data(idx_hh + 417);

    auto tg_yzzzz_yzzzz = pbuffer.data(idx_hh + 418);

    auto tg_yzzzz_zzzzz = pbuffer.data(idx_hh + 419);

    auto tg_zzzzz_xxxxx = pbuffer.data(idx_hh + 420);

    auto tg_zzzzz_xxxxy = pbuffer.data(idx_hh + 421);

    auto tg_zzzzz_xxxxz = pbuffer.data(idx_hh + 422);

    auto tg_zzzzz_xxxyy = pbuffer.data(idx_hh + 423);

    auto tg_zzzzz_xxxyz = pbuffer.data(idx_hh + 424);

    auto tg_zzzzz_xxxzz = pbuffer.data(idx_hh + 425);

    auto tg_zzzzz_xxyyy = pbuffer.data(idx_hh + 426);

    auto tg_zzzzz_xxyyz = pbuffer.data(idx_hh + 427);

    auto tg_zzzzz_xxyzz = pbuffer.data(idx_hh + 428);

    auto tg_zzzzz_xxzzz = pbuffer.data(idx_hh + 429);

    auto tg_zzzzz_xyyyy = pbuffer.data(idx_hh + 430);

    auto tg_zzzzz_xyyyz = pbuffer.data(idx_hh + 431);

    auto tg_zzzzz_xyyzz = pbuffer.data(idx_hh + 432);

    auto tg_zzzzz_xyzzz = pbuffer.data(idx_hh + 433);

    auto tg_zzzzz_xzzzz = pbuffer.data(idx_hh + 434);

    auto tg_zzzzz_yyyyy = pbuffer.data(idx_hh + 435);

    auto tg_zzzzz_yyyyz = pbuffer.data(idx_hh + 436);

    auto tg_zzzzz_yyyzz = pbuffer.data(idx_hh + 437);

    auto tg_zzzzz_yyzzz = pbuffer.data(idx_hh + 438);

    auto tg_zzzzz_yzzzz = pbuffer.data(idx_hh + 439);

    auto tg_zzzzz_zzzzz = pbuffer.data(idx_hh + 440);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxx_xxxxx, tg_xxx_xxxxy, tg_xxx_xxxxz, tg_xxx_xxxyy, tg_xxx_xxxyz, tg_xxx_xxxzz, tg_xxx_xxyyy, tg_xxx_xxyyz, tg_xxx_xxyzz, tg_xxx_xxzzz, tg_xxx_xyyyy, tg_xxx_xyyyz, tg_xxx_xyyzz, tg_xxx_xyzzz, tg_xxx_xzzzz, tg_xxx_yyyyy, tg_xxx_yyyyz, tg_xxx_yyyzz, tg_xxx_yyzzz, tg_xxx_yzzzz, tg_xxx_zzzzz, tg_xxxx_xxxx, tg_xxxx_xxxxx, tg_xxxx_xxxxy, tg_xxxx_xxxxz, tg_xxxx_xxxy, tg_xxxx_xxxyy, tg_xxxx_xxxyz, tg_xxxx_xxxz, tg_xxxx_xxxzz, tg_xxxx_xxyy, tg_xxxx_xxyyy, tg_xxxx_xxyyz, tg_xxxx_xxyz, tg_xxxx_xxyzz, tg_xxxx_xxzz, tg_xxxx_xxzzz, tg_xxxx_xyyy, tg_xxxx_xyyyy, tg_xxxx_xyyyz, tg_xxxx_xyyz, tg_xxxx_xyyzz, tg_xxxx_xyzz, tg_xxxx_xyzzz, tg_xxxx_xzzz, tg_xxxx_xzzzz, tg_xxxx_yyyy, tg_xxxx_yyyyy, tg_xxxx_yyyyz, tg_xxxx_yyyz, tg_xxxx_yyyzz, tg_xxxx_yyzz, tg_xxxx_yyzzz, tg_xxxx_yzzz, tg_xxxx_yzzzz, tg_xxxx_zzzz, tg_xxxx_zzzzz, tg_xxxxx_xxxxx, tg_xxxxx_xxxxy, tg_xxxxx_xxxxz, tg_xxxxx_xxxyy, tg_xxxxx_xxxyz, tg_xxxxx_xxxzz, tg_xxxxx_xxyyy, tg_xxxxx_xxyyz, tg_xxxxx_xxyzz, tg_xxxxx_xxzzz, tg_xxxxx_xyyyy, tg_xxxxx_xyyyz, tg_xxxxx_xyyzz, tg_xxxxx_xyzzz, tg_xxxxx_xzzzz, tg_xxxxx_yyyyy, tg_xxxxx_yyyyz, tg_xxxxx_yyyzz, tg_xxxxx_yyzzz, tg_xxxxx_yzzzz, tg_xxxxx_zzzzz, tg_xxxxy_xxxxx, tg_xxxxy_xxxxy, tg_xxxxy_xxxxz, tg_xxxxy_xxxyy, tg_xxxxy_xxxyz, tg_xxxxy_xxxzz, tg_xxxxy_xxyyy, tg_xxxxy_xxyyz, tg_xxxxy_xxyzz, tg_xxxxy_xxzzz, tg_xxxxy_xyyyy, tg_xxxxy_xyyyz, tg_xxxxy_xyyzz, tg_xxxxy_xyzzz, tg_xxxxy_xzzzz, tg_xxxxy_yyyyy, tg_xxxxy_yyyyz, tg_xxxxy_yyyzz, tg_xxxxy_yyzzz, tg_xxxxy_yzzzz, tg_xxxxy_zzzzz, tg_xxxxz_xxxxx, tg_xxxxz_xxxxy, tg_xxxxz_xxxxz, tg_xxxxz_xxxyy, tg_xxxxz_xxxyz, tg_xxxxz_xxxzz, tg_xxxxz_xxyyy, tg_xxxxz_xxyyz, tg_xxxxz_xxyzz, tg_xxxxz_xxzzz, tg_xxxxz_xyyyy, tg_xxxxz_xyyyz, tg_xxxxz_xyyzz, tg_xxxxz_xyzzz, tg_xxxxz_xzzzz, tg_xxxxz_yyyyy, tg_xxxxz_yyyyz, tg_xxxxz_yyyzz, tg_xxxxz_yyzzz, tg_xxxxz_yzzzz, tg_xxxxz_zzzzz, tg_xxxy_xxxxx, tg_xxxy_xxxxy, tg_xxxy_xxxxz, tg_xxxy_xxxyy, tg_xxxy_xxxzz, tg_xxxy_xxyyy, tg_xxxy_xxzzz, tg_xxxy_xyyyy, tg_xxxy_xzzzz, tg_xxxy_yyyyy, tg_xxxy_yyyyz, tg_xxxy_yyyzz, tg_xxxy_yyzzz, tg_xxxy_yzzzz, tg_xxxyy_xxxxx, tg_xxxyy_xxxxy, tg_xxxyy_xxxxz, tg_xxxyy_xxxyy, tg_xxxyy_xxxyz, tg_xxxyy_xxxzz, tg_xxxyy_xxyyy, tg_xxxyy_xxyyz, tg_xxxyy_xxyzz, tg_xxxyy_xxzzz, tg_xxxyy_xyyyy, tg_xxxyy_xyyyz, tg_xxxyy_xyyzz, tg_xxxyy_xyzzz, tg_xxxyy_xzzzz, tg_xxxyy_yyyyy, tg_xxxyy_yyyyz, tg_xxxyy_yyyzz, tg_xxxyy_yyzzz, tg_xxxyy_yzzzz, tg_xxxyy_zzzzz, tg_xxxyz_xxxxx, tg_xxxyz_xxxxy, tg_xxxyz_xxxxz, tg_xxxyz_xxxyy, tg_xxxyz_xxxyz, tg_xxxyz_xxxzz, tg_xxxyz_xxyyy, tg_xxxyz_xxyyz, tg_xxxyz_xxyzz, tg_xxxyz_xxzzz, tg_xxxyz_xyyyy, tg_xxxyz_xyyyz, tg_xxxyz_xyyzz, tg_xxxyz_xyzzz, tg_xxxyz_xzzzz, tg_xxxyz_yyyyy, tg_xxxyz_yyyyz, tg_xxxyz_yyyzz, tg_xxxyz_yyzzz, tg_xxxyz_yzzzz, tg_xxxyz_zzzzz, tg_xxxz_xxxxx, tg_xxxz_xxxxy, tg_xxxz_xxxxz, tg_xxxz_xxxyy, tg_xxxz_xxxyz, tg_xxxz_xxxz, tg_xxxz_xxxzz, tg_xxxz_xxyyy, tg_xxxz_xxyyz, tg_xxxz_xxyz, tg_xxxz_xxyzz, tg_xxxz_xxzz, tg_xxxz_xxzzz, tg_xxxz_xyyyy, tg_xxxz_xyyyz, tg_xxxz_xyyz, tg_xxxz_xyyzz, tg_xxxz_xyzz, tg_xxxz_xyzzz, tg_xxxz_xzzz, tg_xxxz_xzzzz, tg_xxxz_yyyyz, tg_xxxz_yyyzz, tg_xxxz_yyzzz, tg_xxxz_yzzzz, tg_xxxz_zzzzz, tg_xxxzz_xxxxx, tg_xxxzz_xxxxy, tg_xxxzz_xxxxz, tg_xxxzz_xxxyy, tg_xxxzz_xxxyz, tg_xxxzz_xxxzz, tg_xxxzz_xxyyy, tg_xxxzz_xxyyz, tg_xxxzz_xxyzz, tg_xxxzz_xxzzz, tg_xxxzz_xyyyy, tg_xxxzz_xyyyz, tg_xxxzz_xyyzz, tg_xxxzz_xyzzz, tg_xxxzz_xzzzz, tg_xxxzz_yyyyy, tg_xxxzz_yyyyz, tg_xxxzz_yyyzz, tg_xxxzz_yyzzz, tg_xxxzz_yzzzz, tg_xxxzz_zzzzz, tg_xxy_xxxxx, tg_xxy_xxxxz, tg_xxy_xxxzz, tg_xxy_xxzzz, tg_xxy_xzzzz, tg_xxy_yyyyy, tg_xxy_yyyyz, tg_xxy_yyyzz, tg_xxy_yyzzz, tg_xxy_yzzzz, tg_xxyy_xxxxx, tg_xxyy_xxxxy, tg_xxyy_xxxxz, tg_xxyy_xxxy, tg_xxyy_xxxyy, tg_xxyy_xxxyz, tg_xxyy_xxxzz, tg_xxyy_xxyy, tg_xxyy_xxyyy, tg_xxyy_xxyyz, tg_xxyy_xxyz, tg_xxyy_xxyzz, tg_xxyy_xxzzz, tg_xxyy_xyyy, tg_xxyy_xyyyy, tg_xxyy_xyyyz, tg_xxyy_xyyz, tg_xxyy_xyyzz, tg_xxyy_xyzz, tg_xxyy_xyzzz, tg_xxyy_xzzzz, tg_xxyy_yyyy, tg_xxyy_yyyyy, tg_xxyy_yyyyz, tg_xxyy_yyyz, tg_xxyy_yyyzz, tg_xxyy_yyzz, tg_xxyy_yyzzz, tg_xxyy_yzzz, tg_xxyy_yzzzz, tg_xxyy_zzzzz, tg_xxyyy_xxxxx, tg_xxyyy_xxxxy, tg_xxyyy_xxxxz, tg_xxyyy_xxxyy, tg_xxyyy_xxxyz, tg_xxyyy_xxxzz, tg_xxyyy_xxyyy, tg_xxyyy_xxyyz, tg_xxyyy_xxyzz, tg_xxyyy_xxzzz, tg_xxyyy_xyyyy, tg_xxyyy_xyyyz, tg_xxyyy_xyyzz, tg_xxyyy_xyzzz, tg_xxyyy_xzzzz, tg_xxyyy_yyyyy, tg_xxyyy_yyyyz, tg_xxyyy_yyyzz, tg_xxyyy_yyzzz, tg_xxyyy_yzzzz, tg_xxyyy_zzzzz, tg_xxyyz_xxxxx, tg_xxyyz_xxxxy, tg_xxyyz_xxxxz, tg_xxyyz_xxxyy, tg_xxyyz_xxxyz, tg_xxyyz_xxxzz, tg_xxyyz_xxyyy, tg_xxyyz_xxyyz, tg_xxyyz_xxyzz, tg_xxyyz_xxzzz, tg_xxyyz_xyyyy, tg_xxyyz_xyyyz, tg_xxyyz_xyyzz, tg_xxyyz_xyzzz, tg_xxyyz_xzzzz, tg_xxyyz_yyyyy, tg_xxyyz_yyyyz, tg_xxyyz_yyyzz, tg_xxyyz_yyzzz, tg_xxyyz_yzzzz, tg_xxyyz_zzzzz, tg_xxyz_xxxxz, tg_xxyz_xxxzz, tg_xxyz_xxzzz, tg_xxyz_xzzzz, tg_xxyz_yyyyz, tg_xxyz_yyyzz, tg_xxyz_yyzzz, tg_xxyz_yzzzz, tg_xxyzz_xxxxx, tg_xxyzz_xxxxy, tg_xxyzz_xxxxz, tg_xxyzz_xxxyy, tg_xxyzz_xxxyz, tg_xxyzz_xxxzz, tg_xxyzz_xxyyy, tg_xxyzz_xxyyz, tg_xxyzz_xxyzz, tg_xxyzz_xxzzz, tg_xxyzz_xyyyy, tg_xxyzz_xyyyz, tg_xxyzz_xyyzz, tg_xxyzz_xyzzz, tg_xxyzz_xzzzz, tg_xxyzz_yyyyy, tg_xxyzz_yyyyz, tg_xxyzz_yyyzz, tg_xxyzz_yyzzz, tg_xxyzz_yzzzz, tg_xxyzz_zzzzz, tg_xxz_xxxxx, tg_xxz_xxxxy, tg_xxz_xxxxz, tg_xxz_xxxyy, tg_xxz_xxxzz, tg_xxz_xxyyy, tg_xxz_xxzzz, tg_xxz_xyyyy, tg_xxz_xzzzz, tg_xxz_yyyyz, tg_xxz_yyyzz, tg_xxz_yyzzz, tg_xxz_yzzzz, tg_xxz_zzzzz, tg_xxzz_xxxx, tg_xxzz_xxxxx, tg_xxzz_xxxxy, tg_xxzz_xxxxz, tg_xxzz_xxxy, tg_xxzz_xxxyy, tg_xxzz_xxxyz, tg_xxzz_xxxz, tg_xxzz_xxxzz, tg_xxzz_xxyy, tg_xxzz_xxyyy, tg_xxzz_xxyyz, tg_xxzz_xxyz, tg_xxzz_xxyzz, tg_xxzz_xxzz, tg_xxzz_xxzzz, tg_xxzz_xyyy, tg_xxzz_xyyyy, tg_xxzz_xyyyz, tg_xxzz_xyyz, tg_xxzz_xyyzz, tg_xxzz_xyzz, tg_xxzz_xyzzz, tg_xxzz_xzzz, tg_xxzz_xzzzz, tg_xxzz_yyyyy, tg_xxzz_yyyyz, tg_xxzz_yyyz, tg_xxzz_yyyzz, tg_xxzz_yyzz, tg_xxzz_yyzzz, tg_xxzz_yzzz, tg_xxzz_yzzzz, tg_xxzz_zzzz, tg_xxzz_zzzzz, tg_xxzzz_xxxxx, tg_xxzzz_xxxxy, tg_xxzzz_xxxxz, tg_xxzzz_xxxyy, tg_xxzzz_xxxyz, tg_xxzzz_xxxzz, tg_xxzzz_xxyyy, tg_xxzzz_xxyyz, tg_xxzzz_xxyzz, tg_xxzzz_xxzzz, tg_xxzzz_xyyyy, tg_xxzzz_xyyyz, tg_xxzzz_xyyzz, tg_xxzzz_xyzzz, tg_xxzzz_xzzzz, tg_xxzzz_yyyyy, tg_xxzzz_yyyyz, tg_xxzzz_yyyzz, tg_xxzzz_yyzzz, tg_xxzzz_yzzzz, tg_xxzzz_zzzzz, tg_xyy_xxxxy, tg_xyy_xxxyy, tg_xyy_xxxyz, tg_xyy_xxyyy, tg_xyy_xxyyz, tg_xyy_xxyzz, tg_xyy_xyyyy, tg_xyy_xyyyz, tg_xyy_xyyzz, tg_xyy_xyzzz, tg_xyy_yyyyy, tg_xyy_yyyyz, tg_xyy_yyyzz, tg_xyy_yyzzz, tg_xyy_yzzzz, tg_xyy_zzzzz, tg_xyyy_xxxxx, tg_xyyy_xxxxy, tg_xyyy_xxxy, tg_xyyy_xxxyy, tg_xyyy_xxxyz, tg_xyyy_xxyy, tg_xyyy_xxyyy, tg_xyyy_xxyyz, tg_xyyy_xxyz, tg_xyyy_xxyzz, tg_xyyy_xyyy, tg_xyyy_xyyyy, tg_xyyy_xyyyz, tg_xyyy_xyyz, tg_xyyy_xyyzz, tg_xyyy_xyzz, tg_xyyy_xyzzz, tg_xyyy_yyyy, tg_xyyy_yyyyy, tg_xyyy_yyyyz, tg_xyyy_yyyz, tg_xyyy_yyyzz, tg_xyyy_yyzz, tg_xyyy_yyzzz, tg_xyyy_yzzz, tg_xyyy_yzzzz, tg_xyyy_zzzzz, tg_xyyyy_xxxxx, tg_xyyyy_xxxxy, tg_xyyyy_xxxxz, tg_xyyyy_xxxyy, tg_xyyyy_xxxyz, tg_xyyyy_xxxzz, tg_xyyyy_xxyyy, tg_xyyyy_xxyyz, tg_xyyyy_xxyzz, tg_xyyyy_xxzzz, tg_xyyyy_xyyyy, tg_xyyyy_xyyyz, tg_xyyyy_xyyzz, tg_xyyyy_xyzzz, tg_xyyyy_xzzzz, tg_xyyyy_yyyyy, tg_xyyyy_yyyyz, tg_xyyyy_yyyzz, tg_xyyyy_yyzzz, tg_xyyyy_yzzzz, tg_xyyyy_zzzzz, tg_xyyyz_xxxxx, tg_xyyyz_xxxxy, tg_xyyyz_xxxxz, tg_xyyyz_xxxyy, tg_xyyyz_xxxyz, tg_xyyyz_xxxzz, tg_xyyyz_xxyyy, tg_xyyyz_xxyyz, tg_xyyyz_xxyzz, tg_xyyyz_xxzzz, tg_xyyyz_xyyyy, tg_xyyyz_xyyyz, tg_xyyyz_xyyzz, tg_xyyyz_xyzzz, tg_xyyyz_xzzzz, tg_xyyyz_yyyyy, tg_xyyyz_yyyyz, tg_xyyyz_yyyzz, tg_xyyyz_yyzzz, tg_xyyyz_yzzzz, tg_xyyyz_zzzzz, tg_xyyz_yyyyz, tg_xyyz_yyyzz, tg_xyyz_yyzzz, tg_xyyz_yzzzz, tg_xyyz_zzzzz, tg_xyyzz_xxxxx, tg_xyyzz_xxxxy, tg_xyyzz_xxxxz, tg_xyyzz_xxxyy, tg_xyyzz_xxxyz, tg_xyyzz_xxxzz, tg_xyyzz_xxyyy, tg_xyyzz_xxyyz, tg_xyyzz_xxyzz, tg_xyyzz_xxzzz, tg_xyyzz_xyyyy, tg_xyyzz_xyyyz, tg_xyyzz_xyyzz, tg_xyyzz_xyzzz, tg_xyyzz_xzzzz, tg_xyyzz_yyyyy, tg_xyyzz_yyyyz, tg_xyyzz_yyyzz, tg_xyyzz_yyzzz, tg_xyyzz_yzzzz, tg_xyyzz_zzzzz, tg_xyz_yyyyz, tg_xyz_yyyzz, tg_xyz_yyzzz, tg_xyz_yzzzz, tg_xyzz_yyyyy, tg_xyzz_yyyyz, tg_xyzz_yyyzz, tg_xyzz_yyzzz, tg_xyzz_yzzzz, tg_xyzzz_xxxxx, tg_xyzzz_xxxxy, tg_xyzzz_xxxxz, tg_xyzzz_xxxyy, tg_xyzzz_xxxyz, tg_xyzzz_xxxzz, tg_xyzzz_xxyyy, tg_xyzzz_xxyyz, tg_xyzzz_xxyzz, tg_xyzzz_xxzzz, tg_xyzzz_xyyyy, tg_xyzzz_xyyyz, tg_xyzzz_xyyzz, tg_xyzzz_xyzzz, tg_xyzzz_xzzzz, tg_xyzzz_yyyyy, tg_xyzzz_yyyyz, tg_xyzzz_yyyzz, tg_xyzzz_yyzzz, tg_xyzzz_yzzzz, tg_xyzzz_zzzzz, tg_xzz_xxxxz, tg_xzz_xxxyz, tg_xzz_xxxzz, tg_xzz_xxyyz, tg_xzz_xxyzz, tg_xzz_xxzzz, tg_xzz_xyyyz, tg_xzz_xyyzz, tg_xzz_xyzzz, tg_xzz_xzzzz, tg_xzz_yyyyy, tg_xzz_yyyyz, tg_xzz_yyyzz, tg_xzz_yyzzz, tg_xzz_yzzzz, tg_xzz_zzzzz, tg_xzzz_xxxxx, tg_xzzz_xxxxz, tg_xzzz_xxxyz, tg_xzzz_xxxz, tg_xzzz_xxxzz, tg_xzzz_xxyyz, tg_xzzz_xxyz, tg_xzzz_xxyzz, tg_xzzz_xxzz, tg_xzzz_xxzzz, tg_xzzz_xyyyz, tg_xzzz_xyyz, tg_xzzz_xyyzz, tg_xzzz_xyzz, tg_xzzz_xyzzz, tg_xzzz_xzzz, tg_xzzz_xzzzz, tg_xzzz_yyyyy, tg_xzzz_yyyyz, tg_xzzz_yyyz, tg_xzzz_yyyzz, tg_xzzz_yyzz, tg_xzzz_yyzzz, tg_xzzz_yzzz, tg_xzzz_yzzzz, tg_xzzz_zzzz, tg_xzzz_zzzzz, tg_xzzzz_xxxxx, tg_xzzzz_xxxxy, tg_xzzzz_xxxxz, tg_xzzzz_xxxyy, tg_xzzzz_xxxyz, tg_xzzzz_xxxzz, tg_xzzzz_xxyyy, tg_xzzzz_xxyyz, tg_xzzzz_xxyzz, tg_xzzzz_xxzzz, tg_xzzzz_xyyyy, tg_xzzzz_xyyyz, tg_xzzzz_xyyzz, tg_xzzzz_xyzzz, tg_xzzzz_xzzzz, tg_xzzzz_yyyyy, tg_xzzzz_yyyyz, tg_xzzzz_yyyzz, tg_xzzzz_yyzzz, tg_xzzzz_yzzzz, tg_xzzzz_zzzzz, tg_yyy_xxxxx, tg_yyy_xxxxy, tg_yyy_xxxxz, tg_yyy_xxxyy, tg_yyy_xxxyz, tg_yyy_xxxzz, tg_yyy_xxyyy, tg_yyy_xxyyz, tg_yyy_xxyzz, tg_yyy_xxzzz, tg_yyy_xyyyy, tg_yyy_xyyyz, tg_yyy_xyyzz, tg_yyy_xyzzz, tg_yyy_xzzzz, tg_yyy_yyyyy, tg_yyy_yyyyz, tg_yyy_yyyzz, tg_yyy_yyzzz, tg_yyy_yzzzz, tg_yyy_zzzzz, tg_yyyy_xxxx, tg_yyyy_xxxxx, tg_yyyy_xxxxy, tg_yyyy_xxxxz, tg_yyyy_xxxy, tg_yyyy_xxxyy, tg_yyyy_xxxyz, tg_yyyy_xxxz, tg_yyyy_xxxzz, tg_yyyy_xxyy, tg_yyyy_xxyyy, tg_yyyy_xxyyz, tg_yyyy_xxyz, tg_yyyy_xxyzz, tg_yyyy_xxzz, tg_yyyy_xxzzz, tg_yyyy_xyyy, tg_yyyy_xyyyy, tg_yyyy_xyyyz, tg_yyyy_xyyz, tg_yyyy_xyyzz, tg_yyyy_xyzz, tg_yyyy_xyzzz, tg_yyyy_xzzz, tg_yyyy_xzzzz, tg_yyyy_yyyy, tg_yyyy_yyyyy, tg_yyyy_yyyyz, tg_yyyy_yyyz, tg_yyyy_yyyzz, tg_yyyy_yyzz, tg_yyyy_yyzzz, tg_yyyy_yzzz, tg_yyyy_yzzzz, tg_yyyy_zzzz, tg_yyyy_zzzzz, tg_yyyyy_xxxxx, tg_yyyyy_xxxxy, tg_yyyyy_xxxxz, tg_yyyyy_xxxyy, tg_yyyyy_xxxyz, tg_yyyyy_xxxzz, tg_yyyyy_xxyyy, tg_yyyyy_xxyyz, tg_yyyyy_xxyzz, tg_yyyyy_xxzzz, tg_yyyyy_xyyyy, tg_yyyyy_xyyyz, tg_yyyyy_xyyzz, tg_yyyyy_xyzzz, tg_yyyyy_xzzzz, tg_yyyyy_yyyyy, tg_yyyyy_yyyyz, tg_yyyyy_yyyzz, tg_yyyyy_yyzzz, tg_yyyyy_yzzzz, tg_yyyyy_zzzzz, tg_yyyyz_xxxxx, tg_yyyyz_xxxxy, tg_yyyyz_xxxxz, tg_yyyyz_xxxyy, tg_yyyyz_xxxyz, tg_yyyyz_xxxzz, tg_yyyyz_xxyyy, tg_yyyyz_xxyyz, tg_yyyyz_xxyzz, tg_yyyyz_xxzzz, tg_yyyyz_xyyyy, tg_yyyyz_xyyyz, tg_yyyyz_xyyzz, tg_yyyyz_xyzzz, tg_yyyyz_xzzzz, tg_yyyyz_yyyyy, tg_yyyyz_yyyyz, tg_yyyyz_yyyzz, tg_yyyyz_yyzzz, tg_yyyyz_yzzzz, tg_yyyyz_zzzzz, tg_yyyz_xxxxy, tg_yyyz_xxxxz, tg_yyyz_xxxyy, tg_yyyz_xxxyz, tg_yyyz_xxxz, tg_yyyz_xxxzz, tg_yyyz_xxyyy, tg_yyyz_xxyyz, tg_yyyz_xxyz, tg_yyyz_xxyzz, tg_yyyz_xxzz, tg_yyyz_xxzzz, tg_yyyz_xyyyy, tg_yyyz_xyyyz, tg_yyyz_xyyz, tg_yyyz_xyyzz, tg_yyyz_xyzz, tg_yyyz_xyzzz, tg_yyyz_xzzz, tg_yyyz_xzzzz, tg_yyyz_yyyyy, tg_yyyz_yyyyz, tg_yyyz_yyyz, tg_yyyz_yyyzz, tg_yyyz_yyzz, tg_yyyz_yyzzz, tg_yyyz_yzzz, tg_yyyz_yzzzz, tg_yyyz_zzzz, tg_yyyz_zzzzz, tg_yyyzz_xxxxx, tg_yyyzz_xxxxy, tg_yyyzz_xxxxz, tg_yyyzz_xxxyy, tg_yyyzz_xxxyz, tg_yyyzz_xxxzz, tg_yyyzz_xxyyy, tg_yyyzz_xxyyz, tg_yyyzz_xxyzz, tg_yyyzz_xxzzz, tg_yyyzz_xyyyy, tg_yyyzz_xyyyz, tg_yyyzz_xyyzz, tg_yyyzz_xyzzz, tg_yyyzz_xzzzz, tg_yyyzz_yyyyy, tg_yyyzz_yyyyz, tg_yyyzz_yyyzz, tg_yyyzz_yyzzz, tg_yyyzz_yzzzz, tg_yyyzz_zzzzz, tg_yyz_xxxxy, tg_yyz_xxxxz, tg_yyz_xxxyy, tg_yyz_xxxzz, tg_yyz_xxyyy, tg_yyz_xxzzz, tg_yyz_xyyyy, tg_yyz_xzzzz, tg_yyz_yyyyy, tg_yyz_yyyyz, tg_yyz_yyyzz, tg_yyz_yyzzz, tg_yyz_yzzzz, tg_yyz_zzzzz, tg_yyzz_xxxx, tg_yyzz_xxxxx, tg_yyzz_xxxxy, tg_yyzz_xxxxz, tg_yyzz_xxxy, tg_yyzz_xxxyy, tg_yyzz_xxxyz, tg_yyzz_xxxz, tg_yyzz_xxxzz, tg_yyzz_xxyy, tg_yyzz_xxyyy, tg_yyzz_xxyyz, tg_yyzz_xxyz, tg_yyzz_xxyzz, tg_yyzz_xxzz, tg_yyzz_xxzzz, tg_yyzz_xyyy, tg_yyzz_xyyyy, tg_yyzz_xyyyz, tg_yyzz_xyyz, tg_yyzz_xyyzz, tg_yyzz_xyzz, tg_yyzz_xyzzz, tg_yyzz_xzzz, tg_yyzz_xzzzz, tg_yyzz_yyyy, tg_yyzz_yyyyy, tg_yyzz_yyyyz, tg_yyzz_yyyz, tg_yyzz_yyyzz, tg_yyzz_yyzz, tg_yyzz_yyzzz, tg_yyzz_yzzz, tg_yyzz_yzzzz, tg_yyzz_zzzz, tg_yyzz_zzzzz, tg_yyzzz_xxxxx, tg_yyzzz_xxxxy, tg_yyzzz_xxxxz, tg_yyzzz_xxxyy, tg_yyzzz_xxxyz, tg_yyzzz_xxxzz, tg_yyzzz_xxyyy, tg_yyzzz_xxyyz, tg_yyzzz_xxyzz, tg_yyzzz_xxzzz, tg_yyzzz_xyyyy, tg_yyzzz_xyyyz, tg_yyzzz_xyyzz, tg_yyzzz_xyzzz, tg_yyzzz_xzzzz, tg_yyzzz_yyyyy, tg_yyzzz_yyyyz, tg_yyzzz_yyyzz, tg_yyzzz_yyzzz, tg_yyzzz_yzzzz, tg_yyzzz_zzzzz, tg_yzz_xxxxx, tg_yzz_xxxxz, tg_yzz_xxxyz, tg_yzz_xxxzz, tg_yzz_xxyyz, tg_yzz_xxyzz, tg_yzz_xxzzz, tg_yzz_xyyyz, tg_yzz_xyyzz, tg_yzz_xyzzz, tg_yzz_xzzzz, tg_yzz_yyyyy, tg_yzz_yyyyz, tg_yzz_yyyzz, tg_yzz_yyzzz, tg_yzz_yzzzz, tg_yzz_zzzzz, tg_yzzz_xxxxx, tg_yzzz_xxxxy, tg_yzzz_xxxxz, tg_yzzz_xxxy, tg_yzzz_xxxyy, tg_yzzz_xxxyz, tg_yzzz_xxxz, tg_yzzz_xxxzz, tg_yzzz_xxyy, tg_yzzz_xxyyy, tg_yzzz_xxyyz, tg_yzzz_xxyz, tg_yzzz_xxyzz, tg_yzzz_xxzz, tg_yzzz_xxzzz, tg_yzzz_xyyy, tg_yzzz_xyyyy, tg_yzzz_xyyyz, tg_yzzz_xyyz, tg_yzzz_xyyzz, tg_yzzz_xyzz, tg_yzzz_xyzzz, tg_yzzz_xzzz, tg_yzzz_xzzzz, tg_yzzz_yyyy, tg_yzzz_yyyyy, tg_yzzz_yyyyz, tg_yzzz_yyyz, tg_yzzz_yyyzz, tg_yzzz_yyzz, tg_yzzz_yyzzz, tg_yzzz_yzzz, tg_yzzz_yzzzz, tg_yzzz_zzzz, tg_yzzz_zzzzz, tg_yzzzz_xxxxx, tg_yzzzz_xxxxy, tg_yzzzz_xxxxz, tg_yzzzz_xxxyy, tg_yzzzz_xxxyz, tg_yzzzz_xxxzz, tg_yzzzz_xxyyy, tg_yzzzz_xxyyz, tg_yzzzz_xxyzz, tg_yzzzz_xxzzz, tg_yzzzz_xyyyy, tg_yzzzz_xyyyz, tg_yzzzz_xyyzz, tg_yzzzz_xyzzz, tg_yzzzz_xzzzz, tg_yzzzz_yyyyy, tg_yzzzz_yyyyz, tg_yzzzz_yyyzz, tg_yzzzz_yyzzz, tg_yzzzz_yzzzz, tg_yzzzz_zzzzz, tg_zzz_xxxxx, tg_zzz_xxxxy, tg_zzz_xxxxz, tg_zzz_xxxyy, tg_zzz_xxxyz, tg_zzz_xxxzz, tg_zzz_xxyyy, tg_zzz_xxyyz, tg_zzz_xxyzz, tg_zzz_xxzzz, tg_zzz_xyyyy, tg_zzz_xyyyz, tg_zzz_xyyzz, tg_zzz_xyzzz, tg_zzz_xzzzz, tg_zzz_yyyyy, tg_zzz_yyyyz, tg_zzz_yyyzz, tg_zzz_yyzzz, tg_zzz_yzzzz, tg_zzz_zzzzz, tg_zzzz_xxxx, tg_zzzz_xxxxx, tg_zzzz_xxxxy, tg_zzzz_xxxxz, tg_zzzz_xxxy, tg_zzzz_xxxyy, tg_zzzz_xxxyz, tg_zzzz_xxxz, tg_zzzz_xxxzz, tg_zzzz_xxyy, tg_zzzz_xxyyy, tg_zzzz_xxyyz, tg_zzzz_xxyz, tg_zzzz_xxyzz, tg_zzzz_xxzz, tg_zzzz_xxzzz, tg_zzzz_xyyy, tg_zzzz_xyyyy, tg_zzzz_xyyyz, tg_zzzz_xyyz, tg_zzzz_xyyzz, tg_zzzz_xyzz, tg_zzzz_xyzzz, tg_zzzz_xzzz, tg_zzzz_xzzzz, tg_zzzz_yyyy, tg_zzzz_yyyyy, tg_zzzz_yyyyz, tg_zzzz_yyyz, tg_zzzz_yyyzz, tg_zzzz_yyzz, tg_zzzz_yyzzz, tg_zzzz_yzzz, tg_zzzz_yzzzz, tg_zzzz_zzzz, tg_zzzz_zzzzz, tg_zzzzz_xxxxx, tg_zzzzz_xxxxy, tg_zzzzz_xxxxz, tg_zzzzz_xxxyy, tg_zzzzz_xxxyz, tg_zzzzz_xxxzz, tg_zzzzz_xxyyy, tg_zzzzz_xxyyz, tg_zzzzz_xxyzz, tg_zzzzz_xxzzz, tg_zzzzz_xyyyy, tg_zzzzz_xyyyz, tg_zzzzz_xyyzz, tg_zzzzz_xyzzz, tg_zzzzz_xzzzz, tg_zzzzz_yyyyy, tg_zzzzz_yyyyz, tg_zzzzz_yyyzz, tg_zzzzz_yyzzz, tg_zzzzz_yzzzz, tg_zzzzz_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxx_xxxxx[i] = 4.0 * tg_xxx_xxxxx[i] * fxi[i] + 5.0 * tg_xxxx_xxxx[i] * fxi[i] + tg_xxxx_xxxxx[i] * ra_x[i];

        tg_xxxxx_xxxxy[i] = 4.0 * tg_xxx_xxxxy[i] * fxi[i] + 4.0 * tg_xxxx_xxxy[i] * fxi[i] + tg_xxxx_xxxxy[i] * ra_x[i];

        tg_xxxxx_xxxxz[i] = 4.0 * tg_xxx_xxxxz[i] * fxi[i] + 4.0 * tg_xxxx_xxxz[i] * fxi[i] + tg_xxxx_xxxxz[i] * ra_x[i];

        tg_xxxxx_xxxyy[i] = 4.0 * tg_xxx_xxxyy[i] * fxi[i] + 3.0 * tg_xxxx_xxyy[i] * fxi[i] + tg_xxxx_xxxyy[i] * ra_x[i];

        tg_xxxxx_xxxyz[i] = 4.0 * tg_xxx_xxxyz[i] * fxi[i] + 3.0 * tg_xxxx_xxyz[i] * fxi[i] + tg_xxxx_xxxyz[i] * ra_x[i];

        tg_xxxxx_xxxzz[i] = 4.0 * tg_xxx_xxxzz[i] * fxi[i] + 3.0 * tg_xxxx_xxzz[i] * fxi[i] + tg_xxxx_xxxzz[i] * ra_x[i];

        tg_xxxxx_xxyyy[i] = 4.0 * tg_xxx_xxyyy[i] * fxi[i] + 2.0 * tg_xxxx_xyyy[i] * fxi[i] + tg_xxxx_xxyyy[i] * ra_x[i];

        tg_xxxxx_xxyyz[i] = 4.0 * tg_xxx_xxyyz[i] * fxi[i] + 2.0 * tg_xxxx_xyyz[i] * fxi[i] + tg_xxxx_xxyyz[i] * ra_x[i];

        tg_xxxxx_xxyzz[i] = 4.0 * tg_xxx_xxyzz[i] * fxi[i] + 2.0 * tg_xxxx_xyzz[i] * fxi[i] + tg_xxxx_xxyzz[i] * ra_x[i];

        tg_xxxxx_xxzzz[i] = 4.0 * tg_xxx_xxzzz[i] * fxi[i] + 2.0 * tg_xxxx_xzzz[i] * fxi[i] + tg_xxxx_xxzzz[i] * ra_x[i];

        tg_xxxxx_xyyyy[i] = 4.0 * tg_xxx_xyyyy[i] * fxi[i] + tg_xxxx_yyyy[i] * fxi[i] + tg_xxxx_xyyyy[i] * ra_x[i];

        tg_xxxxx_xyyyz[i] = 4.0 * tg_xxx_xyyyz[i] * fxi[i] + tg_xxxx_yyyz[i] * fxi[i] + tg_xxxx_xyyyz[i] * ra_x[i];

        tg_xxxxx_xyyzz[i] = 4.0 * tg_xxx_xyyzz[i] * fxi[i] + tg_xxxx_yyzz[i] * fxi[i] + tg_xxxx_xyyzz[i] * ra_x[i];

        tg_xxxxx_xyzzz[i] = 4.0 * tg_xxx_xyzzz[i] * fxi[i] + tg_xxxx_yzzz[i] * fxi[i] + tg_xxxx_xyzzz[i] * ra_x[i];

        tg_xxxxx_xzzzz[i] = 4.0 * tg_xxx_xzzzz[i] * fxi[i] + tg_xxxx_zzzz[i] * fxi[i] + tg_xxxx_xzzzz[i] * ra_x[i];

        tg_xxxxx_yyyyy[i] = 4.0 * tg_xxx_yyyyy[i] * fxi[i] + tg_xxxx_yyyyy[i] * ra_x[i];

        tg_xxxxx_yyyyz[i] = 4.0 * tg_xxx_yyyyz[i] * fxi[i] + tg_xxxx_yyyyz[i] * ra_x[i];

        tg_xxxxx_yyyzz[i] = 4.0 * tg_xxx_yyyzz[i] * fxi[i] + tg_xxxx_yyyzz[i] * ra_x[i];

        tg_xxxxx_yyzzz[i] = 4.0 * tg_xxx_yyzzz[i] * fxi[i] + tg_xxxx_yyzzz[i] * ra_x[i];

        tg_xxxxx_yzzzz[i] = 4.0 * tg_xxx_yzzzz[i] * fxi[i] + tg_xxxx_yzzzz[i] * ra_x[i];

        tg_xxxxx_zzzzz[i] = 4.0 * tg_xxx_zzzzz[i] * fxi[i] + tg_xxxx_zzzzz[i] * ra_x[i];

        tg_xxxxy_xxxxx[i] = tg_xxxx_xxxxx[i] * ra_y[i];

        tg_xxxxy_xxxxy[i] = tg_xxxx_xxxx[i] * fxi[i] + tg_xxxx_xxxxy[i] * ra_y[i];

        tg_xxxxy_xxxxz[i] = tg_xxxx_xxxxz[i] * ra_y[i];

        tg_xxxxy_xxxyy[i] = 2.0 * tg_xxxx_xxxy[i] * fxi[i] + tg_xxxx_xxxyy[i] * ra_y[i];

        tg_xxxxy_xxxyz[i] = tg_xxxx_xxxz[i] * fxi[i] + tg_xxxx_xxxyz[i] * ra_y[i];

        tg_xxxxy_xxxzz[i] = tg_xxxx_xxxzz[i] * ra_y[i];

        tg_xxxxy_xxyyy[i] = 3.0 * tg_xxxx_xxyy[i] * fxi[i] + tg_xxxx_xxyyy[i] * ra_y[i];

        tg_xxxxy_xxyyz[i] = 2.0 * tg_xxxx_xxyz[i] * fxi[i] + tg_xxxx_xxyyz[i] * ra_y[i];

        tg_xxxxy_xxyzz[i] = tg_xxxx_xxzz[i] * fxi[i] + tg_xxxx_xxyzz[i] * ra_y[i];

        tg_xxxxy_xxzzz[i] = tg_xxxx_xxzzz[i] * ra_y[i];

        tg_xxxxy_xyyyy[i] = 4.0 * tg_xxxx_xyyy[i] * fxi[i] + tg_xxxx_xyyyy[i] * ra_y[i];

        tg_xxxxy_xyyyz[i] = 3.0 * tg_xxxx_xyyz[i] * fxi[i] + tg_xxxx_xyyyz[i] * ra_y[i];

        tg_xxxxy_xyyzz[i] = 2.0 * tg_xxxx_xyzz[i] * fxi[i] + tg_xxxx_xyyzz[i] * ra_y[i];

        tg_xxxxy_xyzzz[i] = tg_xxxx_xzzz[i] * fxi[i] + tg_xxxx_xyzzz[i] * ra_y[i];

        tg_xxxxy_xzzzz[i] = tg_xxxx_xzzzz[i] * ra_y[i];

        tg_xxxxy_yyyyy[i] = 3.0 * tg_xxy_yyyyy[i] * fxi[i] + tg_xxxy_yyyyy[i] * ra_x[i];

        tg_xxxxy_yyyyz[i] = 3.0 * tg_xxy_yyyyz[i] * fxi[i] + tg_xxxy_yyyyz[i] * ra_x[i];

        tg_xxxxy_yyyzz[i] = 3.0 * tg_xxy_yyyzz[i] * fxi[i] + tg_xxxy_yyyzz[i] * ra_x[i];

        tg_xxxxy_yyzzz[i] = 3.0 * tg_xxy_yyzzz[i] * fxi[i] + tg_xxxy_yyzzz[i] * ra_x[i];

        tg_xxxxy_yzzzz[i] = 3.0 * tg_xxy_yzzzz[i] * fxi[i] + tg_xxxy_yzzzz[i] * ra_x[i];

        tg_xxxxy_zzzzz[i] = tg_xxxx_zzzzz[i] * ra_y[i];

        tg_xxxxz_xxxxx[i] = tg_xxxx_xxxxx[i] * ra_z[i];

        tg_xxxxz_xxxxy[i] = tg_xxxx_xxxxy[i] * ra_z[i];

        tg_xxxxz_xxxxz[i] = tg_xxxx_xxxx[i] * fxi[i] + tg_xxxx_xxxxz[i] * ra_z[i];

        tg_xxxxz_xxxyy[i] = tg_xxxx_xxxyy[i] * ra_z[i];

        tg_xxxxz_xxxyz[i] = tg_xxxx_xxxy[i] * fxi[i] + tg_xxxx_xxxyz[i] * ra_z[i];

        tg_xxxxz_xxxzz[i] = 2.0 * tg_xxxx_xxxz[i] * fxi[i] + tg_xxxx_xxxzz[i] * ra_z[i];

        tg_xxxxz_xxyyy[i] = tg_xxxx_xxyyy[i] * ra_z[i];

        tg_xxxxz_xxyyz[i] = tg_xxxx_xxyy[i] * fxi[i] + tg_xxxx_xxyyz[i] * ra_z[i];

        tg_xxxxz_xxyzz[i] = 2.0 * tg_xxxx_xxyz[i] * fxi[i] + tg_xxxx_xxyzz[i] * ra_z[i];

        tg_xxxxz_xxzzz[i] = 3.0 * tg_xxxx_xxzz[i] * fxi[i] + tg_xxxx_xxzzz[i] * ra_z[i];

        tg_xxxxz_xyyyy[i] = tg_xxxx_xyyyy[i] * ra_z[i];

        tg_xxxxz_xyyyz[i] = tg_xxxx_xyyy[i] * fxi[i] + tg_xxxx_xyyyz[i] * ra_z[i];

        tg_xxxxz_xyyzz[i] = 2.0 * tg_xxxx_xyyz[i] * fxi[i] + tg_xxxx_xyyzz[i] * ra_z[i];

        tg_xxxxz_xyzzz[i] = 3.0 * tg_xxxx_xyzz[i] * fxi[i] + tg_xxxx_xyzzz[i] * ra_z[i];

        tg_xxxxz_xzzzz[i] = 4.0 * tg_xxxx_xzzz[i] * fxi[i] + tg_xxxx_xzzzz[i] * ra_z[i];

        tg_xxxxz_yyyyy[i] = tg_xxxx_yyyyy[i] * ra_z[i];

        tg_xxxxz_yyyyz[i] = 3.0 * tg_xxz_yyyyz[i] * fxi[i] + tg_xxxz_yyyyz[i] * ra_x[i];

        tg_xxxxz_yyyzz[i] = 3.0 * tg_xxz_yyyzz[i] * fxi[i] + tg_xxxz_yyyzz[i] * ra_x[i];

        tg_xxxxz_yyzzz[i] = 3.0 * tg_xxz_yyzzz[i] * fxi[i] + tg_xxxz_yyzzz[i] * ra_x[i];

        tg_xxxxz_yzzzz[i] = 3.0 * tg_xxz_yzzzz[i] * fxi[i] + tg_xxxz_yzzzz[i] * ra_x[i];

        tg_xxxxz_zzzzz[i] = 3.0 * tg_xxz_zzzzz[i] * fxi[i] + tg_xxxz_zzzzz[i] * ra_x[i];

        tg_xxxyy_xxxxx[i] = tg_xxx_xxxxx[i] * fxi[i] + tg_xxxy_xxxxx[i] * ra_y[i];

        tg_xxxyy_xxxxy[i] = 2.0 * tg_xyy_xxxxy[i] * fxi[i] + 4.0 * tg_xxyy_xxxy[i] * fxi[i] + tg_xxyy_xxxxy[i] * ra_x[i];

        tg_xxxyy_xxxxz[i] = tg_xxx_xxxxz[i] * fxi[i] + tg_xxxy_xxxxz[i] * ra_y[i];

        tg_xxxyy_xxxyy[i] = 2.0 * tg_xyy_xxxyy[i] * fxi[i] + 3.0 * tg_xxyy_xxyy[i] * fxi[i] + tg_xxyy_xxxyy[i] * ra_x[i];

        tg_xxxyy_xxxyz[i] = 2.0 * tg_xyy_xxxyz[i] * fxi[i] + 3.0 * tg_xxyy_xxyz[i] * fxi[i] + tg_xxyy_xxxyz[i] * ra_x[i];

        tg_xxxyy_xxxzz[i] = tg_xxx_xxxzz[i] * fxi[i] + tg_xxxy_xxxzz[i] * ra_y[i];

        tg_xxxyy_xxyyy[i] = 2.0 * tg_xyy_xxyyy[i] * fxi[i] + 2.0 * tg_xxyy_xyyy[i] * fxi[i] + tg_xxyy_xxyyy[i] * ra_x[i];

        tg_xxxyy_xxyyz[i] = 2.0 * tg_xyy_xxyyz[i] * fxi[i] + 2.0 * tg_xxyy_xyyz[i] * fxi[i] + tg_xxyy_xxyyz[i] * ra_x[i];

        tg_xxxyy_xxyzz[i] = 2.0 * tg_xyy_xxyzz[i] * fxi[i] + 2.0 * tg_xxyy_xyzz[i] * fxi[i] + tg_xxyy_xxyzz[i] * ra_x[i];

        tg_xxxyy_xxzzz[i] = tg_xxx_xxzzz[i] * fxi[i] + tg_xxxy_xxzzz[i] * ra_y[i];

        tg_xxxyy_xyyyy[i] = 2.0 * tg_xyy_xyyyy[i] * fxi[i] + tg_xxyy_yyyy[i] * fxi[i] + tg_xxyy_xyyyy[i] * ra_x[i];

        tg_xxxyy_xyyyz[i] = 2.0 * tg_xyy_xyyyz[i] * fxi[i] + tg_xxyy_yyyz[i] * fxi[i] + tg_xxyy_xyyyz[i] * ra_x[i];

        tg_xxxyy_xyyzz[i] = 2.0 * tg_xyy_xyyzz[i] * fxi[i] + tg_xxyy_yyzz[i] * fxi[i] + tg_xxyy_xyyzz[i] * ra_x[i];

        tg_xxxyy_xyzzz[i] = 2.0 * tg_xyy_xyzzz[i] * fxi[i] + tg_xxyy_yzzz[i] * fxi[i] + tg_xxyy_xyzzz[i] * ra_x[i];

        tg_xxxyy_xzzzz[i] = tg_xxx_xzzzz[i] * fxi[i] + tg_xxxy_xzzzz[i] * ra_y[i];

        tg_xxxyy_yyyyy[i] = 2.0 * tg_xyy_yyyyy[i] * fxi[i] + tg_xxyy_yyyyy[i] * ra_x[i];

        tg_xxxyy_yyyyz[i] = 2.0 * tg_xyy_yyyyz[i] * fxi[i] + tg_xxyy_yyyyz[i] * ra_x[i];

        tg_xxxyy_yyyzz[i] = 2.0 * tg_xyy_yyyzz[i] * fxi[i] + tg_xxyy_yyyzz[i] * ra_x[i];

        tg_xxxyy_yyzzz[i] = 2.0 * tg_xyy_yyzzz[i] * fxi[i] + tg_xxyy_yyzzz[i] * ra_x[i];

        tg_xxxyy_yzzzz[i] = 2.0 * tg_xyy_yzzzz[i] * fxi[i] + tg_xxyy_yzzzz[i] * ra_x[i];

        tg_xxxyy_zzzzz[i] = 2.0 * tg_xyy_zzzzz[i] * fxi[i] + tg_xxyy_zzzzz[i] * ra_x[i];

        tg_xxxyz_xxxxx[i] = tg_xxxz_xxxxx[i] * ra_y[i];

        tg_xxxyz_xxxxy[i] = tg_xxxy_xxxxy[i] * ra_z[i];

        tg_xxxyz_xxxxz[i] = tg_xxxz_xxxxz[i] * ra_y[i];

        tg_xxxyz_xxxyy[i] = tg_xxxy_xxxyy[i] * ra_z[i];

        tg_xxxyz_xxxyz[i] = tg_xxxz_xxxz[i] * fxi[i] + tg_xxxz_xxxyz[i] * ra_y[i];

        tg_xxxyz_xxxzz[i] = tg_xxxz_xxxzz[i] * ra_y[i];

        tg_xxxyz_xxyyy[i] = tg_xxxy_xxyyy[i] * ra_z[i];

        tg_xxxyz_xxyyz[i] = 2.0 * tg_xxxz_xxyz[i] * fxi[i] + tg_xxxz_xxyyz[i] * ra_y[i];

        tg_xxxyz_xxyzz[i] = tg_xxxz_xxzz[i] * fxi[i] + tg_xxxz_xxyzz[i] * ra_y[i];

        tg_xxxyz_xxzzz[i] = tg_xxxz_xxzzz[i] * ra_y[i];

        tg_xxxyz_xyyyy[i] = tg_xxxy_xyyyy[i] * ra_z[i];

        tg_xxxyz_xyyyz[i] = 3.0 * tg_xxxz_xyyz[i] * fxi[i] + tg_xxxz_xyyyz[i] * ra_y[i];

        tg_xxxyz_xyyzz[i] = 2.0 * tg_xxxz_xyzz[i] * fxi[i] + tg_xxxz_xyyzz[i] * ra_y[i];

        tg_xxxyz_xyzzz[i] = tg_xxxz_xzzz[i] * fxi[i] + tg_xxxz_xyzzz[i] * ra_y[i];

        tg_xxxyz_xzzzz[i] = tg_xxxz_xzzzz[i] * ra_y[i];

        tg_xxxyz_yyyyy[i] = tg_xxxy_yyyyy[i] * ra_z[i];

        tg_xxxyz_yyyyz[i] = 2.0 * tg_xyz_yyyyz[i] * fxi[i] + tg_xxyz_yyyyz[i] * ra_x[i];

        tg_xxxyz_yyyzz[i] = 2.0 * tg_xyz_yyyzz[i] * fxi[i] + tg_xxyz_yyyzz[i] * ra_x[i];

        tg_xxxyz_yyzzz[i] = 2.0 * tg_xyz_yyzzz[i] * fxi[i] + tg_xxyz_yyzzz[i] * ra_x[i];

        tg_xxxyz_yzzzz[i] = 2.0 * tg_xyz_yzzzz[i] * fxi[i] + tg_xxyz_yzzzz[i] * ra_x[i];

        tg_xxxyz_zzzzz[i] = tg_xxxz_zzzzz[i] * ra_y[i];

        tg_xxxzz_xxxxx[i] = tg_xxx_xxxxx[i] * fxi[i] + tg_xxxz_xxxxx[i] * ra_z[i];

        tg_xxxzz_xxxxy[i] = tg_xxx_xxxxy[i] * fxi[i] + tg_xxxz_xxxxy[i] * ra_z[i];

        tg_xxxzz_xxxxz[i] = 2.0 * tg_xzz_xxxxz[i] * fxi[i] + 4.0 * tg_xxzz_xxxz[i] * fxi[i] + tg_xxzz_xxxxz[i] * ra_x[i];

        tg_xxxzz_xxxyy[i] = tg_xxx_xxxyy[i] * fxi[i] + tg_xxxz_xxxyy[i] * ra_z[i];

        tg_xxxzz_xxxyz[i] = 2.0 * tg_xzz_xxxyz[i] * fxi[i] + 3.0 * tg_xxzz_xxyz[i] * fxi[i] + tg_xxzz_xxxyz[i] * ra_x[i];

        tg_xxxzz_xxxzz[i] = 2.0 * tg_xzz_xxxzz[i] * fxi[i] + 3.0 * tg_xxzz_xxzz[i] * fxi[i] + tg_xxzz_xxxzz[i] * ra_x[i];

        tg_xxxzz_xxyyy[i] = tg_xxx_xxyyy[i] * fxi[i] + tg_xxxz_xxyyy[i] * ra_z[i];

        tg_xxxzz_xxyyz[i] = 2.0 * tg_xzz_xxyyz[i] * fxi[i] + 2.0 * tg_xxzz_xyyz[i] * fxi[i] + tg_xxzz_xxyyz[i] * ra_x[i];

        tg_xxxzz_xxyzz[i] = 2.0 * tg_xzz_xxyzz[i] * fxi[i] + 2.0 * tg_xxzz_xyzz[i] * fxi[i] + tg_xxzz_xxyzz[i] * ra_x[i];

        tg_xxxzz_xxzzz[i] = 2.0 * tg_xzz_xxzzz[i] * fxi[i] + 2.0 * tg_xxzz_xzzz[i] * fxi[i] + tg_xxzz_xxzzz[i] * ra_x[i];

        tg_xxxzz_xyyyy[i] = tg_xxx_xyyyy[i] * fxi[i] + tg_xxxz_xyyyy[i] * ra_z[i];

        tg_xxxzz_xyyyz[i] = 2.0 * tg_xzz_xyyyz[i] * fxi[i] + tg_xxzz_yyyz[i] * fxi[i] + tg_xxzz_xyyyz[i] * ra_x[i];

        tg_xxxzz_xyyzz[i] = 2.0 * tg_xzz_xyyzz[i] * fxi[i] + tg_xxzz_yyzz[i] * fxi[i] + tg_xxzz_xyyzz[i] * ra_x[i];

        tg_xxxzz_xyzzz[i] = 2.0 * tg_xzz_xyzzz[i] * fxi[i] + tg_xxzz_yzzz[i] * fxi[i] + tg_xxzz_xyzzz[i] * ra_x[i];

        tg_xxxzz_xzzzz[i] = 2.0 * tg_xzz_xzzzz[i] * fxi[i] + tg_xxzz_zzzz[i] * fxi[i] + tg_xxzz_xzzzz[i] * ra_x[i];

        tg_xxxzz_yyyyy[i] = 2.0 * tg_xzz_yyyyy[i] * fxi[i] + tg_xxzz_yyyyy[i] * ra_x[i];

        tg_xxxzz_yyyyz[i] = 2.0 * tg_xzz_yyyyz[i] * fxi[i] + tg_xxzz_yyyyz[i] * ra_x[i];

        tg_xxxzz_yyyzz[i] = 2.0 * tg_xzz_yyyzz[i] * fxi[i] + tg_xxzz_yyyzz[i] * ra_x[i];

        tg_xxxzz_yyzzz[i] = 2.0 * tg_xzz_yyzzz[i] * fxi[i] + tg_xxzz_yyzzz[i] * ra_x[i];

        tg_xxxzz_yzzzz[i] = 2.0 * tg_xzz_yzzzz[i] * fxi[i] + tg_xxzz_yzzzz[i] * ra_x[i];

        tg_xxxzz_zzzzz[i] = 2.0 * tg_xzz_zzzzz[i] * fxi[i] + tg_xxzz_zzzzz[i] * ra_x[i];

        tg_xxyyy_xxxxx[i] = 2.0 * tg_xxy_xxxxx[i] * fxi[i] + tg_xxyy_xxxxx[i] * ra_y[i];

        tg_xxyyy_xxxxy[i] = tg_yyy_xxxxy[i] * fxi[i] + 4.0 * tg_xyyy_xxxy[i] * fxi[i] + tg_xyyy_xxxxy[i] * ra_x[i];

        tg_xxyyy_xxxxz[i] = 2.0 * tg_xxy_xxxxz[i] * fxi[i] + tg_xxyy_xxxxz[i] * ra_y[i];

        tg_xxyyy_xxxyy[i] = tg_yyy_xxxyy[i] * fxi[i] + 3.0 * tg_xyyy_xxyy[i] * fxi[i] + tg_xyyy_xxxyy[i] * ra_x[i];

        tg_xxyyy_xxxyz[i] = tg_yyy_xxxyz[i] * fxi[i] + 3.0 * tg_xyyy_xxyz[i] * fxi[i] + tg_xyyy_xxxyz[i] * ra_x[i];

        tg_xxyyy_xxxzz[i] = 2.0 * tg_xxy_xxxzz[i] * fxi[i] + tg_xxyy_xxxzz[i] * ra_y[i];

        tg_xxyyy_xxyyy[i] = tg_yyy_xxyyy[i] * fxi[i] + 2.0 * tg_xyyy_xyyy[i] * fxi[i] + tg_xyyy_xxyyy[i] * ra_x[i];

        tg_xxyyy_xxyyz[i] = tg_yyy_xxyyz[i] * fxi[i] + 2.0 * tg_xyyy_xyyz[i] * fxi[i] + tg_xyyy_xxyyz[i] * ra_x[i];

        tg_xxyyy_xxyzz[i] = tg_yyy_xxyzz[i] * fxi[i] + 2.0 * tg_xyyy_xyzz[i] * fxi[i] + tg_xyyy_xxyzz[i] * ra_x[i];

        tg_xxyyy_xxzzz[i] = 2.0 * tg_xxy_xxzzz[i] * fxi[i] + tg_xxyy_xxzzz[i] * ra_y[i];

        tg_xxyyy_xyyyy[i] = tg_yyy_xyyyy[i] * fxi[i] + tg_xyyy_yyyy[i] * fxi[i] + tg_xyyy_xyyyy[i] * ra_x[i];

        tg_xxyyy_xyyyz[i] = tg_yyy_xyyyz[i] * fxi[i] + tg_xyyy_yyyz[i] * fxi[i] + tg_xyyy_xyyyz[i] * ra_x[i];

        tg_xxyyy_xyyzz[i] = tg_yyy_xyyzz[i] * fxi[i] + tg_xyyy_yyzz[i] * fxi[i] + tg_xyyy_xyyzz[i] * ra_x[i];

        tg_xxyyy_xyzzz[i] = tg_yyy_xyzzz[i] * fxi[i] + tg_xyyy_yzzz[i] * fxi[i] + tg_xyyy_xyzzz[i] * ra_x[i];

        tg_xxyyy_xzzzz[i] = 2.0 * tg_xxy_xzzzz[i] * fxi[i] + tg_xxyy_xzzzz[i] * ra_y[i];

        tg_xxyyy_yyyyy[i] = tg_yyy_yyyyy[i] * fxi[i] + tg_xyyy_yyyyy[i] * ra_x[i];

        tg_xxyyy_yyyyz[i] = tg_yyy_yyyyz[i] * fxi[i] + tg_xyyy_yyyyz[i] * ra_x[i];

        tg_xxyyy_yyyzz[i] = tg_yyy_yyyzz[i] * fxi[i] + tg_xyyy_yyyzz[i] * ra_x[i];

        tg_xxyyy_yyzzz[i] = tg_yyy_yyzzz[i] * fxi[i] + tg_xyyy_yyzzz[i] * ra_x[i];

        tg_xxyyy_yzzzz[i] = tg_yyy_yzzzz[i] * fxi[i] + tg_xyyy_yzzzz[i] * ra_x[i];

        tg_xxyyy_zzzzz[i] = tg_yyy_zzzzz[i] * fxi[i] + tg_xyyy_zzzzz[i] * ra_x[i];

        tg_xxyyz_xxxxx[i] = tg_xxyy_xxxxx[i] * ra_z[i];

        tg_xxyyz_xxxxy[i] = tg_xxyy_xxxxy[i] * ra_z[i];

        tg_xxyyz_xxxxz[i] = tg_xxz_xxxxz[i] * fxi[i] + tg_xxyz_xxxxz[i] * ra_y[i];

        tg_xxyyz_xxxyy[i] = tg_xxyy_xxxyy[i] * ra_z[i];

        tg_xxyyz_xxxyz[i] = tg_xxyy_xxxy[i] * fxi[i] + tg_xxyy_xxxyz[i] * ra_z[i];

        tg_xxyyz_xxxzz[i] = tg_xxz_xxxzz[i] * fxi[i] + tg_xxyz_xxxzz[i] * ra_y[i];

        tg_xxyyz_xxyyy[i] = tg_xxyy_xxyyy[i] * ra_z[i];

        tg_xxyyz_xxyyz[i] = tg_xxyy_xxyy[i] * fxi[i] + tg_xxyy_xxyyz[i] * ra_z[i];

        tg_xxyyz_xxyzz[i] = 2.0 * tg_xxyy_xxyz[i] * fxi[i] + tg_xxyy_xxyzz[i] * ra_z[i];

        tg_xxyyz_xxzzz[i] = tg_xxz_xxzzz[i] * fxi[i] + tg_xxyz_xxzzz[i] * ra_y[i];

        tg_xxyyz_xyyyy[i] = tg_xxyy_xyyyy[i] * ra_z[i];

        tg_xxyyz_xyyyz[i] = tg_xxyy_xyyy[i] * fxi[i] + tg_xxyy_xyyyz[i] * ra_z[i];

        tg_xxyyz_xyyzz[i] = 2.0 * tg_xxyy_xyyz[i] * fxi[i] + tg_xxyy_xyyzz[i] * ra_z[i];

        tg_xxyyz_xyzzz[i] = 3.0 * tg_xxyy_xyzz[i] * fxi[i] + tg_xxyy_xyzzz[i] * ra_z[i];

        tg_xxyyz_xzzzz[i] = tg_xxz_xzzzz[i] * fxi[i] + tg_xxyz_xzzzz[i] * ra_y[i];

        tg_xxyyz_yyyyy[i] = tg_xxyy_yyyyy[i] * ra_z[i];

        tg_xxyyz_yyyyz[i] = tg_yyz_yyyyz[i] * fxi[i] + tg_xyyz_yyyyz[i] * ra_x[i];

        tg_xxyyz_yyyzz[i] = tg_yyz_yyyzz[i] * fxi[i] + tg_xyyz_yyyzz[i] * ra_x[i];

        tg_xxyyz_yyzzz[i] = tg_yyz_yyzzz[i] * fxi[i] + tg_xyyz_yyzzz[i] * ra_x[i];

        tg_xxyyz_yzzzz[i] = tg_yyz_yzzzz[i] * fxi[i] + tg_xyyz_yzzzz[i] * ra_x[i];

        tg_xxyyz_zzzzz[i] = tg_yyz_zzzzz[i] * fxi[i] + tg_xyyz_zzzzz[i] * ra_x[i];

        tg_xxyzz_xxxxx[i] = tg_xxzz_xxxxx[i] * ra_y[i];

        tg_xxyzz_xxxxy[i] = tg_xxzz_xxxx[i] * fxi[i] + tg_xxzz_xxxxy[i] * ra_y[i];

        tg_xxyzz_xxxxz[i] = tg_xxzz_xxxxz[i] * ra_y[i];

        tg_xxyzz_xxxyy[i] = 2.0 * tg_xxzz_xxxy[i] * fxi[i] + tg_xxzz_xxxyy[i] * ra_y[i];

        tg_xxyzz_xxxyz[i] = tg_xxzz_xxxz[i] * fxi[i] + tg_xxzz_xxxyz[i] * ra_y[i];

        tg_xxyzz_xxxzz[i] = tg_xxzz_xxxzz[i] * ra_y[i];

        tg_xxyzz_xxyyy[i] = 3.0 * tg_xxzz_xxyy[i] * fxi[i] + tg_xxzz_xxyyy[i] * ra_y[i];

        tg_xxyzz_xxyyz[i] = 2.0 * tg_xxzz_xxyz[i] * fxi[i] + tg_xxzz_xxyyz[i] * ra_y[i];

        tg_xxyzz_xxyzz[i] = tg_xxzz_xxzz[i] * fxi[i] + tg_xxzz_xxyzz[i] * ra_y[i];

        tg_xxyzz_xxzzz[i] = tg_xxzz_xxzzz[i] * ra_y[i];

        tg_xxyzz_xyyyy[i] = 4.0 * tg_xxzz_xyyy[i] * fxi[i] + tg_xxzz_xyyyy[i] * ra_y[i];

        tg_xxyzz_xyyyz[i] = 3.0 * tg_xxzz_xyyz[i] * fxi[i] + tg_xxzz_xyyyz[i] * ra_y[i];

        tg_xxyzz_xyyzz[i] = 2.0 * tg_xxzz_xyzz[i] * fxi[i] + tg_xxzz_xyyzz[i] * ra_y[i];

        tg_xxyzz_xyzzz[i] = tg_xxzz_xzzz[i] * fxi[i] + tg_xxzz_xyzzz[i] * ra_y[i];

        tg_xxyzz_xzzzz[i] = tg_xxzz_xzzzz[i] * ra_y[i];

        tg_xxyzz_yyyyy[i] = tg_yzz_yyyyy[i] * fxi[i] + tg_xyzz_yyyyy[i] * ra_x[i];

        tg_xxyzz_yyyyz[i] = tg_yzz_yyyyz[i] * fxi[i] + tg_xyzz_yyyyz[i] * ra_x[i];

        tg_xxyzz_yyyzz[i] = tg_yzz_yyyzz[i] * fxi[i] + tg_xyzz_yyyzz[i] * ra_x[i];

        tg_xxyzz_yyzzz[i] = tg_yzz_yyzzz[i] * fxi[i] + tg_xyzz_yyzzz[i] * ra_x[i];

        tg_xxyzz_yzzzz[i] = tg_yzz_yzzzz[i] * fxi[i] + tg_xyzz_yzzzz[i] * ra_x[i];

        tg_xxyzz_zzzzz[i] = tg_xxzz_zzzzz[i] * ra_y[i];

        tg_xxzzz_xxxxx[i] = 2.0 * tg_xxz_xxxxx[i] * fxi[i] + tg_xxzz_xxxxx[i] * ra_z[i];

        tg_xxzzz_xxxxy[i] = 2.0 * tg_xxz_xxxxy[i] * fxi[i] + tg_xxzz_xxxxy[i] * ra_z[i];

        tg_xxzzz_xxxxz[i] = tg_zzz_xxxxz[i] * fxi[i] + 4.0 * tg_xzzz_xxxz[i] * fxi[i] + tg_xzzz_xxxxz[i] * ra_x[i];

        tg_xxzzz_xxxyy[i] = 2.0 * tg_xxz_xxxyy[i] * fxi[i] + tg_xxzz_xxxyy[i] * ra_z[i];

        tg_xxzzz_xxxyz[i] = tg_zzz_xxxyz[i] * fxi[i] + 3.0 * tg_xzzz_xxyz[i] * fxi[i] + tg_xzzz_xxxyz[i] * ra_x[i];

        tg_xxzzz_xxxzz[i] = tg_zzz_xxxzz[i] * fxi[i] + 3.0 * tg_xzzz_xxzz[i] * fxi[i] + tg_xzzz_xxxzz[i] * ra_x[i];

        tg_xxzzz_xxyyy[i] = 2.0 * tg_xxz_xxyyy[i] * fxi[i] + tg_xxzz_xxyyy[i] * ra_z[i];

        tg_xxzzz_xxyyz[i] = tg_zzz_xxyyz[i] * fxi[i] + 2.0 * tg_xzzz_xyyz[i] * fxi[i] + tg_xzzz_xxyyz[i] * ra_x[i];

        tg_xxzzz_xxyzz[i] = tg_zzz_xxyzz[i] * fxi[i] + 2.0 * tg_xzzz_xyzz[i] * fxi[i] + tg_xzzz_xxyzz[i] * ra_x[i];

        tg_xxzzz_xxzzz[i] = tg_zzz_xxzzz[i] * fxi[i] + 2.0 * tg_xzzz_xzzz[i] * fxi[i] + tg_xzzz_xxzzz[i] * ra_x[i];

        tg_xxzzz_xyyyy[i] = 2.0 * tg_xxz_xyyyy[i] * fxi[i] + tg_xxzz_xyyyy[i] * ra_z[i];

        tg_xxzzz_xyyyz[i] = tg_zzz_xyyyz[i] * fxi[i] + tg_xzzz_yyyz[i] * fxi[i] + tg_xzzz_xyyyz[i] * ra_x[i];

        tg_xxzzz_xyyzz[i] = tg_zzz_xyyzz[i] * fxi[i] + tg_xzzz_yyzz[i] * fxi[i] + tg_xzzz_xyyzz[i] * ra_x[i];

        tg_xxzzz_xyzzz[i] = tg_zzz_xyzzz[i] * fxi[i] + tg_xzzz_yzzz[i] * fxi[i] + tg_xzzz_xyzzz[i] * ra_x[i];

        tg_xxzzz_xzzzz[i] = tg_zzz_xzzzz[i] * fxi[i] + tg_xzzz_zzzz[i] * fxi[i] + tg_xzzz_xzzzz[i] * ra_x[i];

        tg_xxzzz_yyyyy[i] = tg_zzz_yyyyy[i] * fxi[i] + tg_xzzz_yyyyy[i] * ra_x[i];

        tg_xxzzz_yyyyz[i] = tg_zzz_yyyyz[i] * fxi[i] + tg_xzzz_yyyyz[i] * ra_x[i];

        tg_xxzzz_yyyzz[i] = tg_zzz_yyyzz[i] * fxi[i] + tg_xzzz_yyyzz[i] * ra_x[i];

        tg_xxzzz_yyzzz[i] = tg_zzz_yyzzz[i] * fxi[i] + tg_xzzz_yyzzz[i] * ra_x[i];

        tg_xxzzz_yzzzz[i] = tg_zzz_yzzzz[i] * fxi[i] + tg_xzzz_yzzzz[i] * ra_x[i];

        tg_xxzzz_zzzzz[i] = tg_zzz_zzzzz[i] * fxi[i] + tg_xzzz_zzzzz[i] * ra_x[i];

        tg_xyyyy_xxxxx[i] = 5.0 * tg_yyyy_xxxx[i] * fxi[i] + tg_yyyy_xxxxx[i] * ra_x[i];

        tg_xyyyy_xxxxy[i] = 4.0 * tg_yyyy_xxxy[i] * fxi[i] + tg_yyyy_xxxxy[i] * ra_x[i];

        tg_xyyyy_xxxxz[i] = 4.0 * tg_yyyy_xxxz[i] * fxi[i] + tg_yyyy_xxxxz[i] * ra_x[i];

        tg_xyyyy_xxxyy[i] = 3.0 * tg_yyyy_xxyy[i] * fxi[i] + tg_yyyy_xxxyy[i] * ra_x[i];

        tg_xyyyy_xxxyz[i] = 3.0 * tg_yyyy_xxyz[i] * fxi[i] + tg_yyyy_xxxyz[i] * ra_x[i];

        tg_xyyyy_xxxzz[i] = 3.0 * tg_yyyy_xxzz[i] * fxi[i] + tg_yyyy_xxxzz[i] * ra_x[i];

        tg_xyyyy_xxyyy[i] = 2.0 * tg_yyyy_xyyy[i] * fxi[i] + tg_yyyy_xxyyy[i] * ra_x[i];

        tg_xyyyy_xxyyz[i] = 2.0 * tg_yyyy_xyyz[i] * fxi[i] + tg_yyyy_xxyyz[i] * ra_x[i];

        tg_xyyyy_xxyzz[i] = 2.0 * tg_yyyy_xyzz[i] * fxi[i] + tg_yyyy_xxyzz[i] * ra_x[i];

        tg_xyyyy_xxzzz[i] = 2.0 * tg_yyyy_xzzz[i] * fxi[i] + tg_yyyy_xxzzz[i] * ra_x[i];

        tg_xyyyy_xyyyy[i] = tg_yyyy_yyyy[i] * fxi[i] + tg_yyyy_xyyyy[i] * ra_x[i];

        tg_xyyyy_xyyyz[i] = tg_yyyy_yyyz[i] * fxi[i] + tg_yyyy_xyyyz[i] * ra_x[i];

        tg_xyyyy_xyyzz[i] = tg_yyyy_yyzz[i] * fxi[i] + tg_yyyy_xyyzz[i] * ra_x[i];

        tg_xyyyy_xyzzz[i] = tg_yyyy_yzzz[i] * fxi[i] + tg_yyyy_xyzzz[i] * ra_x[i];

        tg_xyyyy_xzzzz[i] = tg_yyyy_zzzz[i] * fxi[i] + tg_yyyy_xzzzz[i] * ra_x[i];

        tg_xyyyy_yyyyy[i] = tg_yyyy_yyyyy[i] * ra_x[i];

        tg_xyyyy_yyyyz[i] = tg_yyyy_yyyyz[i] * ra_x[i];

        tg_xyyyy_yyyzz[i] = tg_yyyy_yyyzz[i] * ra_x[i];

        tg_xyyyy_yyzzz[i] = tg_yyyy_yyzzz[i] * ra_x[i];

        tg_xyyyy_yzzzz[i] = tg_yyyy_yzzzz[i] * ra_x[i];

        tg_xyyyy_zzzzz[i] = tg_yyyy_zzzzz[i] * ra_x[i];

        tg_xyyyz_xxxxx[i] = tg_xyyy_xxxxx[i] * ra_z[i];

        tg_xyyyz_xxxxy[i] = tg_xyyy_xxxxy[i] * ra_z[i];

        tg_xyyyz_xxxxz[i] = 4.0 * tg_yyyz_xxxz[i] * fxi[i] + tg_yyyz_xxxxz[i] * ra_x[i];

        tg_xyyyz_xxxyy[i] = tg_xyyy_xxxyy[i] * ra_z[i];

        tg_xyyyz_xxxyz[i] = 3.0 * tg_yyyz_xxyz[i] * fxi[i] + tg_yyyz_xxxyz[i] * ra_x[i];

        tg_xyyyz_xxxzz[i] = 3.0 * tg_yyyz_xxzz[i] * fxi[i] + tg_yyyz_xxxzz[i] * ra_x[i];

        tg_xyyyz_xxyyy[i] = tg_xyyy_xxyyy[i] * ra_z[i];

        tg_xyyyz_xxyyz[i] = 2.0 * tg_yyyz_xyyz[i] * fxi[i] + tg_yyyz_xxyyz[i] * ra_x[i];

        tg_xyyyz_xxyzz[i] = 2.0 * tg_yyyz_xyzz[i] * fxi[i] + tg_yyyz_xxyzz[i] * ra_x[i];

        tg_xyyyz_xxzzz[i] = 2.0 * tg_yyyz_xzzz[i] * fxi[i] + tg_yyyz_xxzzz[i] * ra_x[i];

        tg_xyyyz_xyyyy[i] = tg_xyyy_xyyyy[i] * ra_z[i];

        tg_xyyyz_xyyyz[i] = tg_yyyz_yyyz[i] * fxi[i] + tg_yyyz_xyyyz[i] * ra_x[i];

        tg_xyyyz_xyyzz[i] = tg_yyyz_yyzz[i] * fxi[i] + tg_yyyz_xyyzz[i] * ra_x[i];

        tg_xyyyz_xyzzz[i] = tg_yyyz_yzzz[i] * fxi[i] + tg_yyyz_xyzzz[i] * ra_x[i];

        tg_xyyyz_xzzzz[i] = tg_yyyz_zzzz[i] * fxi[i] + tg_yyyz_xzzzz[i] * ra_x[i];

        tg_xyyyz_yyyyy[i] = tg_yyyz_yyyyy[i] * ra_x[i];

        tg_xyyyz_yyyyz[i] = tg_yyyz_yyyyz[i] * ra_x[i];

        tg_xyyyz_yyyzz[i] = tg_yyyz_yyyzz[i] * ra_x[i];

        tg_xyyyz_yyzzz[i] = tg_yyyz_yyzzz[i] * ra_x[i];

        tg_xyyyz_yzzzz[i] = tg_yyyz_yzzzz[i] * ra_x[i];

        tg_xyyyz_zzzzz[i] = tg_yyyz_zzzzz[i] * ra_x[i];

        tg_xyyzz_xxxxx[i] = 5.0 * tg_yyzz_xxxx[i] * fxi[i] + tg_yyzz_xxxxx[i] * ra_x[i];

        tg_xyyzz_xxxxy[i] = 4.0 * tg_yyzz_xxxy[i] * fxi[i] + tg_yyzz_xxxxy[i] * ra_x[i];

        tg_xyyzz_xxxxz[i] = 4.0 * tg_yyzz_xxxz[i] * fxi[i] + tg_yyzz_xxxxz[i] * ra_x[i];

        tg_xyyzz_xxxyy[i] = 3.0 * tg_yyzz_xxyy[i] * fxi[i] + tg_yyzz_xxxyy[i] * ra_x[i];

        tg_xyyzz_xxxyz[i] = 3.0 * tg_yyzz_xxyz[i] * fxi[i] + tg_yyzz_xxxyz[i] * ra_x[i];

        tg_xyyzz_xxxzz[i] = 3.0 * tg_yyzz_xxzz[i] * fxi[i] + tg_yyzz_xxxzz[i] * ra_x[i];

        tg_xyyzz_xxyyy[i] = 2.0 * tg_yyzz_xyyy[i] * fxi[i] + tg_yyzz_xxyyy[i] * ra_x[i];

        tg_xyyzz_xxyyz[i] = 2.0 * tg_yyzz_xyyz[i] * fxi[i] + tg_yyzz_xxyyz[i] * ra_x[i];

        tg_xyyzz_xxyzz[i] = 2.0 * tg_yyzz_xyzz[i] * fxi[i] + tg_yyzz_xxyzz[i] * ra_x[i];

        tg_xyyzz_xxzzz[i] = 2.0 * tg_yyzz_xzzz[i] * fxi[i] + tg_yyzz_xxzzz[i] * ra_x[i];

        tg_xyyzz_xyyyy[i] = tg_yyzz_yyyy[i] * fxi[i] + tg_yyzz_xyyyy[i] * ra_x[i];

        tg_xyyzz_xyyyz[i] = tg_yyzz_yyyz[i] * fxi[i] + tg_yyzz_xyyyz[i] * ra_x[i];

        tg_xyyzz_xyyzz[i] = tg_yyzz_yyzz[i] * fxi[i] + tg_yyzz_xyyzz[i] * ra_x[i];

        tg_xyyzz_xyzzz[i] = tg_yyzz_yzzz[i] * fxi[i] + tg_yyzz_xyzzz[i] * ra_x[i];

        tg_xyyzz_xzzzz[i] = tg_yyzz_zzzz[i] * fxi[i] + tg_yyzz_xzzzz[i] * ra_x[i];

        tg_xyyzz_yyyyy[i] = tg_yyzz_yyyyy[i] * ra_x[i];

        tg_xyyzz_yyyyz[i] = tg_yyzz_yyyyz[i] * ra_x[i];

        tg_xyyzz_yyyzz[i] = tg_yyzz_yyyzz[i] * ra_x[i];

        tg_xyyzz_yyzzz[i] = tg_yyzz_yyzzz[i] * ra_x[i];

        tg_xyyzz_yzzzz[i] = tg_yyzz_yzzzz[i] * ra_x[i];

        tg_xyyzz_zzzzz[i] = tg_yyzz_zzzzz[i] * ra_x[i];

        tg_xyzzz_xxxxx[i] = tg_xzzz_xxxxx[i] * ra_y[i];

        tg_xyzzz_xxxxy[i] = 4.0 * tg_yzzz_xxxy[i] * fxi[i] + tg_yzzz_xxxxy[i] * ra_x[i];

        tg_xyzzz_xxxxz[i] = tg_xzzz_xxxxz[i] * ra_y[i];

        tg_xyzzz_xxxyy[i] = 3.0 * tg_yzzz_xxyy[i] * fxi[i] + tg_yzzz_xxxyy[i] * ra_x[i];

        tg_xyzzz_xxxyz[i] = 3.0 * tg_yzzz_xxyz[i] * fxi[i] + tg_yzzz_xxxyz[i] * ra_x[i];

        tg_xyzzz_xxxzz[i] = tg_xzzz_xxxzz[i] * ra_y[i];

        tg_xyzzz_xxyyy[i] = 2.0 * tg_yzzz_xyyy[i] * fxi[i] + tg_yzzz_xxyyy[i] * ra_x[i];

        tg_xyzzz_xxyyz[i] = 2.0 * tg_yzzz_xyyz[i] * fxi[i] + tg_yzzz_xxyyz[i] * ra_x[i];

        tg_xyzzz_xxyzz[i] = 2.0 * tg_yzzz_xyzz[i] * fxi[i] + tg_yzzz_xxyzz[i] * ra_x[i];

        tg_xyzzz_xxzzz[i] = tg_xzzz_xxzzz[i] * ra_y[i];

        tg_xyzzz_xyyyy[i] = tg_yzzz_yyyy[i] * fxi[i] + tg_yzzz_xyyyy[i] * ra_x[i];

        tg_xyzzz_xyyyz[i] = tg_yzzz_yyyz[i] * fxi[i] + tg_yzzz_xyyyz[i] * ra_x[i];

        tg_xyzzz_xyyzz[i] = tg_yzzz_yyzz[i] * fxi[i] + tg_yzzz_xyyzz[i] * ra_x[i];

        tg_xyzzz_xyzzz[i] = tg_yzzz_yzzz[i] * fxi[i] + tg_yzzz_xyzzz[i] * ra_x[i];

        tg_xyzzz_xzzzz[i] = tg_xzzz_xzzzz[i] * ra_y[i];

        tg_xyzzz_yyyyy[i] = tg_yzzz_yyyyy[i] * ra_x[i];

        tg_xyzzz_yyyyz[i] = tg_yzzz_yyyyz[i] * ra_x[i];

        tg_xyzzz_yyyzz[i] = tg_yzzz_yyyzz[i] * ra_x[i];

        tg_xyzzz_yyzzz[i] = tg_yzzz_yyzzz[i] * ra_x[i];

        tg_xyzzz_yzzzz[i] = tg_yzzz_yzzzz[i] * ra_x[i];

        tg_xyzzz_zzzzz[i] = tg_yzzz_zzzzz[i] * ra_x[i];

        tg_xzzzz_xxxxx[i] = 5.0 * tg_zzzz_xxxx[i] * fxi[i] + tg_zzzz_xxxxx[i] * ra_x[i];

        tg_xzzzz_xxxxy[i] = 4.0 * tg_zzzz_xxxy[i] * fxi[i] + tg_zzzz_xxxxy[i] * ra_x[i];

        tg_xzzzz_xxxxz[i] = 4.0 * tg_zzzz_xxxz[i] * fxi[i] + tg_zzzz_xxxxz[i] * ra_x[i];

        tg_xzzzz_xxxyy[i] = 3.0 * tg_zzzz_xxyy[i] * fxi[i] + tg_zzzz_xxxyy[i] * ra_x[i];

        tg_xzzzz_xxxyz[i] = 3.0 * tg_zzzz_xxyz[i] * fxi[i] + tg_zzzz_xxxyz[i] * ra_x[i];

        tg_xzzzz_xxxzz[i] = 3.0 * tg_zzzz_xxzz[i] * fxi[i] + tg_zzzz_xxxzz[i] * ra_x[i];

        tg_xzzzz_xxyyy[i] = 2.0 * tg_zzzz_xyyy[i] * fxi[i] + tg_zzzz_xxyyy[i] * ra_x[i];

        tg_xzzzz_xxyyz[i] = 2.0 * tg_zzzz_xyyz[i] * fxi[i] + tg_zzzz_xxyyz[i] * ra_x[i];

        tg_xzzzz_xxyzz[i] = 2.0 * tg_zzzz_xyzz[i] * fxi[i] + tg_zzzz_xxyzz[i] * ra_x[i];

        tg_xzzzz_xxzzz[i] = 2.0 * tg_zzzz_xzzz[i] * fxi[i] + tg_zzzz_xxzzz[i] * ra_x[i];

        tg_xzzzz_xyyyy[i] = tg_zzzz_yyyy[i] * fxi[i] + tg_zzzz_xyyyy[i] * ra_x[i];

        tg_xzzzz_xyyyz[i] = tg_zzzz_yyyz[i] * fxi[i] + tg_zzzz_xyyyz[i] * ra_x[i];

        tg_xzzzz_xyyzz[i] = tg_zzzz_yyzz[i] * fxi[i] + tg_zzzz_xyyzz[i] * ra_x[i];

        tg_xzzzz_xyzzz[i] = tg_zzzz_yzzz[i] * fxi[i] + tg_zzzz_xyzzz[i] * ra_x[i];

        tg_xzzzz_xzzzz[i] = tg_zzzz_zzzz[i] * fxi[i] + tg_zzzz_xzzzz[i] * ra_x[i];

        tg_xzzzz_yyyyy[i] = tg_zzzz_yyyyy[i] * ra_x[i];

        tg_xzzzz_yyyyz[i] = tg_zzzz_yyyyz[i] * ra_x[i];

        tg_xzzzz_yyyzz[i] = tg_zzzz_yyyzz[i] * ra_x[i];

        tg_xzzzz_yyzzz[i] = tg_zzzz_yyzzz[i] * ra_x[i];

        tg_xzzzz_yzzzz[i] = tg_zzzz_yzzzz[i] * ra_x[i];

        tg_xzzzz_zzzzz[i] = tg_zzzz_zzzzz[i] * ra_x[i];

        tg_yyyyy_xxxxx[i] = 4.0 * tg_yyy_xxxxx[i] * fxi[i] + tg_yyyy_xxxxx[i] * ra_y[i];

        tg_yyyyy_xxxxy[i] = 4.0 * tg_yyy_xxxxy[i] * fxi[i] + tg_yyyy_xxxx[i] * fxi[i] + tg_yyyy_xxxxy[i] * ra_y[i];

        tg_yyyyy_xxxxz[i] = 4.0 * tg_yyy_xxxxz[i] * fxi[i] + tg_yyyy_xxxxz[i] * ra_y[i];

        tg_yyyyy_xxxyy[i] = 4.0 * tg_yyy_xxxyy[i] * fxi[i] + 2.0 * tg_yyyy_xxxy[i] * fxi[i] + tg_yyyy_xxxyy[i] * ra_y[i];

        tg_yyyyy_xxxyz[i] = 4.0 * tg_yyy_xxxyz[i] * fxi[i] + tg_yyyy_xxxz[i] * fxi[i] + tg_yyyy_xxxyz[i] * ra_y[i];

        tg_yyyyy_xxxzz[i] = 4.0 * tg_yyy_xxxzz[i] * fxi[i] + tg_yyyy_xxxzz[i] * ra_y[i];

        tg_yyyyy_xxyyy[i] = 4.0 * tg_yyy_xxyyy[i] * fxi[i] + 3.0 * tg_yyyy_xxyy[i] * fxi[i] + tg_yyyy_xxyyy[i] * ra_y[i];

        tg_yyyyy_xxyyz[i] = 4.0 * tg_yyy_xxyyz[i] * fxi[i] + 2.0 * tg_yyyy_xxyz[i] * fxi[i] + tg_yyyy_xxyyz[i] * ra_y[i];

        tg_yyyyy_xxyzz[i] = 4.0 * tg_yyy_xxyzz[i] * fxi[i] + tg_yyyy_xxzz[i] * fxi[i] + tg_yyyy_xxyzz[i] * ra_y[i];

        tg_yyyyy_xxzzz[i] = 4.0 * tg_yyy_xxzzz[i] * fxi[i] + tg_yyyy_xxzzz[i] * ra_y[i];

        tg_yyyyy_xyyyy[i] = 4.0 * tg_yyy_xyyyy[i] * fxi[i] + 4.0 * tg_yyyy_xyyy[i] * fxi[i] + tg_yyyy_xyyyy[i] * ra_y[i];

        tg_yyyyy_xyyyz[i] = 4.0 * tg_yyy_xyyyz[i] * fxi[i] + 3.0 * tg_yyyy_xyyz[i] * fxi[i] + tg_yyyy_xyyyz[i] * ra_y[i];

        tg_yyyyy_xyyzz[i] = 4.0 * tg_yyy_xyyzz[i] * fxi[i] + 2.0 * tg_yyyy_xyzz[i] * fxi[i] + tg_yyyy_xyyzz[i] * ra_y[i];

        tg_yyyyy_xyzzz[i] = 4.0 * tg_yyy_xyzzz[i] * fxi[i] + tg_yyyy_xzzz[i] * fxi[i] + tg_yyyy_xyzzz[i] * ra_y[i];

        tg_yyyyy_xzzzz[i] = 4.0 * tg_yyy_xzzzz[i] * fxi[i] + tg_yyyy_xzzzz[i] * ra_y[i];

        tg_yyyyy_yyyyy[i] = 4.0 * tg_yyy_yyyyy[i] * fxi[i] + 5.0 * tg_yyyy_yyyy[i] * fxi[i] + tg_yyyy_yyyyy[i] * ra_y[i];

        tg_yyyyy_yyyyz[i] = 4.0 * tg_yyy_yyyyz[i] * fxi[i] + 4.0 * tg_yyyy_yyyz[i] * fxi[i] + tg_yyyy_yyyyz[i] * ra_y[i];

        tg_yyyyy_yyyzz[i] = 4.0 * tg_yyy_yyyzz[i] * fxi[i] + 3.0 * tg_yyyy_yyzz[i] * fxi[i] + tg_yyyy_yyyzz[i] * ra_y[i];

        tg_yyyyy_yyzzz[i] = 4.0 * tg_yyy_yyzzz[i] * fxi[i] + 2.0 * tg_yyyy_yzzz[i] * fxi[i] + tg_yyyy_yyzzz[i] * ra_y[i];

        tg_yyyyy_yzzzz[i] = 4.0 * tg_yyy_yzzzz[i] * fxi[i] + tg_yyyy_zzzz[i] * fxi[i] + tg_yyyy_yzzzz[i] * ra_y[i];

        tg_yyyyy_zzzzz[i] = 4.0 * tg_yyy_zzzzz[i] * fxi[i] + tg_yyyy_zzzzz[i] * ra_y[i];

        tg_yyyyz_xxxxx[i] = tg_yyyy_xxxxx[i] * ra_z[i];

        tg_yyyyz_xxxxy[i] = tg_yyyy_xxxxy[i] * ra_z[i];

        tg_yyyyz_xxxxz[i] = 3.0 * tg_yyz_xxxxz[i] * fxi[i] + tg_yyyz_xxxxz[i] * ra_y[i];

        tg_yyyyz_xxxyy[i] = tg_yyyy_xxxyy[i] * ra_z[i];

        tg_yyyyz_xxxyz[i] = tg_yyyy_xxxy[i] * fxi[i] + tg_yyyy_xxxyz[i] * ra_z[i];

        tg_yyyyz_xxxzz[i] = 3.0 * tg_yyz_xxxzz[i] * fxi[i] + tg_yyyz_xxxzz[i] * ra_y[i];

        tg_yyyyz_xxyyy[i] = tg_yyyy_xxyyy[i] * ra_z[i];

        tg_yyyyz_xxyyz[i] = tg_yyyy_xxyy[i] * fxi[i] + tg_yyyy_xxyyz[i] * ra_z[i];

        tg_yyyyz_xxyzz[i] = 2.0 * tg_yyyy_xxyz[i] * fxi[i] + tg_yyyy_xxyzz[i] * ra_z[i];

        tg_yyyyz_xxzzz[i] = 3.0 * tg_yyz_xxzzz[i] * fxi[i] + tg_yyyz_xxzzz[i] * ra_y[i];

        tg_yyyyz_xyyyy[i] = tg_yyyy_xyyyy[i] * ra_z[i];

        tg_yyyyz_xyyyz[i] = tg_yyyy_xyyy[i] * fxi[i] + tg_yyyy_xyyyz[i] * ra_z[i];

        tg_yyyyz_xyyzz[i] = 2.0 * tg_yyyy_xyyz[i] * fxi[i] + tg_yyyy_xyyzz[i] * ra_z[i];

        tg_yyyyz_xyzzz[i] = 3.0 * tg_yyyy_xyzz[i] * fxi[i] + tg_yyyy_xyzzz[i] * ra_z[i];

        tg_yyyyz_xzzzz[i] = 3.0 * tg_yyz_xzzzz[i] * fxi[i] + tg_yyyz_xzzzz[i] * ra_y[i];

        tg_yyyyz_yyyyy[i] = tg_yyyy_yyyyy[i] * ra_z[i];

        tg_yyyyz_yyyyz[i] = tg_yyyy_yyyy[i] * fxi[i] + tg_yyyy_yyyyz[i] * ra_z[i];

        tg_yyyyz_yyyzz[i] = 2.0 * tg_yyyy_yyyz[i] * fxi[i] + tg_yyyy_yyyzz[i] * ra_z[i];

        tg_yyyyz_yyzzz[i] = 3.0 * tg_yyyy_yyzz[i] * fxi[i] + tg_yyyy_yyzzz[i] * ra_z[i];

        tg_yyyyz_yzzzz[i] = 4.0 * tg_yyyy_yzzz[i] * fxi[i] + tg_yyyy_yzzzz[i] * ra_z[i];

        tg_yyyyz_zzzzz[i] = 3.0 * tg_yyz_zzzzz[i] * fxi[i] + tg_yyyz_zzzzz[i] * ra_y[i];

        tg_yyyzz_xxxxx[i] = 2.0 * tg_yzz_xxxxx[i] * fxi[i] + tg_yyzz_xxxxx[i] * ra_y[i];

        tg_yyyzz_xxxxy[i] = tg_yyy_xxxxy[i] * fxi[i] + tg_yyyz_xxxxy[i] * ra_z[i];

        tg_yyyzz_xxxxz[i] = 2.0 * tg_yzz_xxxxz[i] * fxi[i] + tg_yyzz_xxxxz[i] * ra_y[i];

        tg_yyyzz_xxxyy[i] = tg_yyy_xxxyy[i] * fxi[i] + tg_yyyz_xxxyy[i] * ra_z[i];

        tg_yyyzz_xxxyz[i] = 2.0 * tg_yzz_xxxyz[i] * fxi[i] + tg_yyzz_xxxz[i] * fxi[i] + tg_yyzz_xxxyz[i] * ra_y[i];

        tg_yyyzz_xxxzz[i] = 2.0 * tg_yzz_xxxzz[i] * fxi[i] + tg_yyzz_xxxzz[i] * ra_y[i];

        tg_yyyzz_xxyyy[i] = tg_yyy_xxyyy[i] * fxi[i] + tg_yyyz_xxyyy[i] * ra_z[i];

        tg_yyyzz_xxyyz[i] = 2.0 * tg_yzz_xxyyz[i] * fxi[i] + 2.0 * tg_yyzz_xxyz[i] * fxi[i] + tg_yyzz_xxyyz[i] * ra_y[i];

        tg_yyyzz_xxyzz[i] = 2.0 * tg_yzz_xxyzz[i] * fxi[i] + tg_yyzz_xxzz[i] * fxi[i] + tg_yyzz_xxyzz[i] * ra_y[i];

        tg_yyyzz_xxzzz[i] = 2.0 * tg_yzz_xxzzz[i] * fxi[i] + tg_yyzz_xxzzz[i] * ra_y[i];

        tg_yyyzz_xyyyy[i] = tg_yyy_xyyyy[i] * fxi[i] + tg_yyyz_xyyyy[i] * ra_z[i];

        tg_yyyzz_xyyyz[i] = 2.0 * tg_yzz_xyyyz[i] * fxi[i] + 3.0 * tg_yyzz_xyyz[i] * fxi[i] + tg_yyzz_xyyyz[i] * ra_y[i];

        tg_yyyzz_xyyzz[i] = 2.0 * tg_yzz_xyyzz[i] * fxi[i] + 2.0 * tg_yyzz_xyzz[i] * fxi[i] + tg_yyzz_xyyzz[i] * ra_y[i];

        tg_yyyzz_xyzzz[i] = 2.0 * tg_yzz_xyzzz[i] * fxi[i] + tg_yyzz_xzzz[i] * fxi[i] + tg_yyzz_xyzzz[i] * ra_y[i];

        tg_yyyzz_xzzzz[i] = 2.0 * tg_yzz_xzzzz[i] * fxi[i] + tg_yyzz_xzzzz[i] * ra_y[i];

        tg_yyyzz_yyyyy[i] = tg_yyy_yyyyy[i] * fxi[i] + tg_yyyz_yyyyy[i] * ra_z[i];

        tg_yyyzz_yyyyz[i] = 2.0 * tg_yzz_yyyyz[i] * fxi[i] + 4.0 * tg_yyzz_yyyz[i] * fxi[i] + tg_yyzz_yyyyz[i] * ra_y[i];

        tg_yyyzz_yyyzz[i] = 2.0 * tg_yzz_yyyzz[i] * fxi[i] + 3.0 * tg_yyzz_yyzz[i] * fxi[i] + tg_yyzz_yyyzz[i] * ra_y[i];

        tg_yyyzz_yyzzz[i] = 2.0 * tg_yzz_yyzzz[i] * fxi[i] + 2.0 * tg_yyzz_yzzz[i] * fxi[i] + tg_yyzz_yyzzz[i] * ra_y[i];

        tg_yyyzz_yzzzz[i] = 2.0 * tg_yzz_yzzzz[i] * fxi[i] + tg_yyzz_zzzz[i] * fxi[i] + tg_yyzz_yzzzz[i] * ra_y[i];

        tg_yyyzz_zzzzz[i] = 2.0 * tg_yzz_zzzzz[i] * fxi[i] + tg_yyzz_zzzzz[i] * ra_y[i];

        tg_yyzzz_xxxxx[i] = tg_zzz_xxxxx[i] * fxi[i] + tg_yzzz_xxxxx[i] * ra_y[i];

        tg_yyzzz_xxxxy[i] = 2.0 * tg_yyz_xxxxy[i] * fxi[i] + tg_yyzz_xxxxy[i] * ra_z[i];

        tg_yyzzz_xxxxz[i] = tg_zzz_xxxxz[i] * fxi[i] + tg_yzzz_xxxxz[i] * ra_y[i];

        tg_yyzzz_xxxyy[i] = 2.0 * tg_yyz_xxxyy[i] * fxi[i] + tg_yyzz_xxxyy[i] * ra_z[i];

        tg_yyzzz_xxxyz[i] = tg_zzz_xxxyz[i] * fxi[i] + tg_yzzz_xxxz[i] * fxi[i] + tg_yzzz_xxxyz[i] * ra_y[i];

        tg_yyzzz_xxxzz[i] = tg_zzz_xxxzz[i] * fxi[i] + tg_yzzz_xxxzz[i] * ra_y[i];

        tg_yyzzz_xxyyy[i] = 2.0 * tg_yyz_xxyyy[i] * fxi[i] + tg_yyzz_xxyyy[i] * ra_z[i];

        tg_yyzzz_xxyyz[i] = tg_zzz_xxyyz[i] * fxi[i] + 2.0 * tg_yzzz_xxyz[i] * fxi[i] + tg_yzzz_xxyyz[i] * ra_y[i];

        tg_yyzzz_xxyzz[i] = tg_zzz_xxyzz[i] * fxi[i] + tg_yzzz_xxzz[i] * fxi[i] + tg_yzzz_xxyzz[i] * ra_y[i];

        tg_yyzzz_xxzzz[i] = tg_zzz_xxzzz[i] * fxi[i] + tg_yzzz_xxzzz[i] * ra_y[i];

        tg_yyzzz_xyyyy[i] = 2.0 * tg_yyz_xyyyy[i] * fxi[i] + tg_yyzz_xyyyy[i] * ra_z[i];

        tg_yyzzz_xyyyz[i] = tg_zzz_xyyyz[i] * fxi[i] + 3.0 * tg_yzzz_xyyz[i] * fxi[i] + tg_yzzz_xyyyz[i] * ra_y[i];

        tg_yyzzz_xyyzz[i] = tg_zzz_xyyzz[i] * fxi[i] + 2.0 * tg_yzzz_xyzz[i] * fxi[i] + tg_yzzz_xyyzz[i] * ra_y[i];

        tg_yyzzz_xyzzz[i] = tg_zzz_xyzzz[i] * fxi[i] + tg_yzzz_xzzz[i] * fxi[i] + tg_yzzz_xyzzz[i] * ra_y[i];

        tg_yyzzz_xzzzz[i] = tg_zzz_xzzzz[i] * fxi[i] + tg_yzzz_xzzzz[i] * ra_y[i];

        tg_yyzzz_yyyyy[i] = 2.0 * tg_yyz_yyyyy[i] * fxi[i] + tg_yyzz_yyyyy[i] * ra_z[i];

        tg_yyzzz_yyyyz[i] = tg_zzz_yyyyz[i] * fxi[i] + 4.0 * tg_yzzz_yyyz[i] * fxi[i] + tg_yzzz_yyyyz[i] * ra_y[i];

        tg_yyzzz_yyyzz[i] = tg_zzz_yyyzz[i] * fxi[i] + 3.0 * tg_yzzz_yyzz[i] * fxi[i] + tg_yzzz_yyyzz[i] * ra_y[i];

        tg_yyzzz_yyzzz[i] = tg_zzz_yyzzz[i] * fxi[i] + 2.0 * tg_yzzz_yzzz[i] * fxi[i] + tg_yzzz_yyzzz[i] * ra_y[i];

        tg_yyzzz_yzzzz[i] = tg_zzz_yzzzz[i] * fxi[i] + tg_yzzz_zzzz[i] * fxi[i] + tg_yzzz_yzzzz[i] * ra_y[i];

        tg_yyzzz_zzzzz[i] = tg_zzz_zzzzz[i] * fxi[i] + tg_yzzz_zzzzz[i] * ra_y[i];

        tg_yzzzz_xxxxx[i] = tg_zzzz_xxxxx[i] * ra_y[i];

        tg_yzzzz_xxxxy[i] = tg_zzzz_xxxx[i] * fxi[i] + tg_zzzz_xxxxy[i] * ra_y[i];

        tg_yzzzz_xxxxz[i] = tg_zzzz_xxxxz[i] * ra_y[i];

        tg_yzzzz_xxxyy[i] = 2.0 * tg_zzzz_xxxy[i] * fxi[i] + tg_zzzz_xxxyy[i] * ra_y[i];

        tg_yzzzz_xxxyz[i] = tg_zzzz_xxxz[i] * fxi[i] + tg_zzzz_xxxyz[i] * ra_y[i];

        tg_yzzzz_xxxzz[i] = tg_zzzz_xxxzz[i] * ra_y[i];

        tg_yzzzz_xxyyy[i] = 3.0 * tg_zzzz_xxyy[i] * fxi[i] + tg_zzzz_xxyyy[i] * ra_y[i];

        tg_yzzzz_xxyyz[i] = 2.0 * tg_zzzz_xxyz[i] * fxi[i] + tg_zzzz_xxyyz[i] * ra_y[i];

        tg_yzzzz_xxyzz[i] = tg_zzzz_xxzz[i] * fxi[i] + tg_zzzz_xxyzz[i] * ra_y[i];

        tg_yzzzz_xxzzz[i] = tg_zzzz_xxzzz[i] * ra_y[i];

        tg_yzzzz_xyyyy[i] = 4.0 * tg_zzzz_xyyy[i] * fxi[i] + tg_zzzz_xyyyy[i] * ra_y[i];

        tg_yzzzz_xyyyz[i] = 3.0 * tg_zzzz_xyyz[i] * fxi[i] + tg_zzzz_xyyyz[i] * ra_y[i];

        tg_yzzzz_xyyzz[i] = 2.0 * tg_zzzz_xyzz[i] * fxi[i] + tg_zzzz_xyyzz[i] * ra_y[i];

        tg_yzzzz_xyzzz[i] = tg_zzzz_xzzz[i] * fxi[i] + tg_zzzz_xyzzz[i] * ra_y[i];

        tg_yzzzz_xzzzz[i] = tg_zzzz_xzzzz[i] * ra_y[i];

        tg_yzzzz_yyyyy[i] = 5.0 * tg_zzzz_yyyy[i] * fxi[i] + tg_zzzz_yyyyy[i] * ra_y[i];

        tg_yzzzz_yyyyz[i] = 4.0 * tg_zzzz_yyyz[i] * fxi[i] + tg_zzzz_yyyyz[i] * ra_y[i];

        tg_yzzzz_yyyzz[i] = 3.0 * tg_zzzz_yyzz[i] * fxi[i] + tg_zzzz_yyyzz[i] * ra_y[i];

        tg_yzzzz_yyzzz[i] = 2.0 * tg_zzzz_yzzz[i] * fxi[i] + tg_zzzz_yyzzz[i] * ra_y[i];

        tg_yzzzz_yzzzz[i] = tg_zzzz_zzzz[i] * fxi[i] + tg_zzzz_yzzzz[i] * ra_y[i];

        tg_yzzzz_zzzzz[i] = tg_zzzz_zzzzz[i] * ra_y[i];

        tg_zzzzz_xxxxx[i] = 4.0 * tg_zzz_xxxxx[i] * fxi[i] + tg_zzzz_xxxxx[i] * ra_z[i];

        tg_zzzzz_xxxxy[i] = 4.0 * tg_zzz_xxxxy[i] * fxi[i] + tg_zzzz_xxxxy[i] * ra_z[i];

        tg_zzzzz_xxxxz[i] = 4.0 * tg_zzz_xxxxz[i] * fxi[i] + tg_zzzz_xxxx[i] * fxi[i] + tg_zzzz_xxxxz[i] * ra_z[i];

        tg_zzzzz_xxxyy[i] = 4.0 * tg_zzz_xxxyy[i] * fxi[i] + tg_zzzz_xxxyy[i] * ra_z[i];

        tg_zzzzz_xxxyz[i] = 4.0 * tg_zzz_xxxyz[i] * fxi[i] + tg_zzzz_xxxy[i] * fxi[i] + tg_zzzz_xxxyz[i] * ra_z[i];

        tg_zzzzz_xxxzz[i] = 4.0 * tg_zzz_xxxzz[i] * fxi[i] + 2.0 * tg_zzzz_xxxz[i] * fxi[i] + tg_zzzz_xxxzz[i] * ra_z[i];

        tg_zzzzz_xxyyy[i] = 4.0 * tg_zzz_xxyyy[i] * fxi[i] + tg_zzzz_xxyyy[i] * ra_z[i];

        tg_zzzzz_xxyyz[i] = 4.0 * tg_zzz_xxyyz[i] * fxi[i] + tg_zzzz_xxyy[i] * fxi[i] + tg_zzzz_xxyyz[i] * ra_z[i];

        tg_zzzzz_xxyzz[i] = 4.0 * tg_zzz_xxyzz[i] * fxi[i] + 2.0 * tg_zzzz_xxyz[i] * fxi[i] + tg_zzzz_xxyzz[i] * ra_z[i];

        tg_zzzzz_xxzzz[i] = 4.0 * tg_zzz_xxzzz[i] * fxi[i] + 3.0 * tg_zzzz_xxzz[i] * fxi[i] + tg_zzzz_xxzzz[i] * ra_z[i];

        tg_zzzzz_xyyyy[i] = 4.0 * tg_zzz_xyyyy[i] * fxi[i] + tg_zzzz_xyyyy[i] * ra_z[i];

        tg_zzzzz_xyyyz[i] = 4.0 * tg_zzz_xyyyz[i] * fxi[i] + tg_zzzz_xyyy[i] * fxi[i] + tg_zzzz_xyyyz[i] * ra_z[i];

        tg_zzzzz_xyyzz[i] = 4.0 * tg_zzz_xyyzz[i] * fxi[i] + 2.0 * tg_zzzz_xyyz[i] * fxi[i] + tg_zzzz_xyyzz[i] * ra_z[i];

        tg_zzzzz_xyzzz[i] = 4.0 * tg_zzz_xyzzz[i] * fxi[i] + 3.0 * tg_zzzz_xyzz[i] * fxi[i] + tg_zzzz_xyzzz[i] * ra_z[i];

        tg_zzzzz_xzzzz[i] = 4.0 * tg_zzz_xzzzz[i] * fxi[i] + 4.0 * tg_zzzz_xzzz[i] * fxi[i] + tg_zzzz_xzzzz[i] * ra_z[i];

        tg_zzzzz_yyyyy[i] = 4.0 * tg_zzz_yyyyy[i] * fxi[i] + tg_zzzz_yyyyy[i] * ra_z[i];

        tg_zzzzz_yyyyz[i] = 4.0 * tg_zzz_yyyyz[i] * fxi[i] + tg_zzzz_yyyy[i] * fxi[i] + tg_zzzz_yyyyz[i] * ra_z[i];

        tg_zzzzz_yyyzz[i] = 4.0 * tg_zzz_yyyzz[i] * fxi[i] + 2.0 * tg_zzzz_yyyz[i] * fxi[i] + tg_zzzz_yyyzz[i] * ra_z[i];

        tg_zzzzz_yyzzz[i] = 4.0 * tg_zzz_yyzzz[i] * fxi[i] + 3.0 * tg_zzzz_yyzz[i] * fxi[i] + tg_zzzz_yyzzz[i] * ra_z[i];

        tg_zzzzz_yzzzz[i] = 4.0 * tg_zzz_yzzzz[i] * fxi[i] + 4.0 * tg_zzzz_yzzz[i] * fxi[i] + tg_zzzz_yzzzz[i] * ra_z[i];

        tg_zzzzz_zzzzz[i] = 4.0 * tg_zzz_zzzzz[i] * fxi[i] + 5.0 * tg_zzzz_zzzz[i] * fxi[i] + tg_zzzz_zzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

