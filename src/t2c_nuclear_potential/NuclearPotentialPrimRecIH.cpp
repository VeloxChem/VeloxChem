#include "NuclearPotentialPrimRecIH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_ih(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_ih,
                               const size_t              idx_npot_0_gh,
                               const size_t              idx_npot_1_gh,
                               const size_t              idx_npot_0_hg,
                               const size_t              idx_npot_1_hg,
                               const size_t              idx_npot_0_hh,
                               const size_t              idx_npot_1_hh,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpa,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : GH

    auto ta_xxxx_xxxxx_0 = pbuffer.data(idx_npot_0_gh);

    auto ta_xxxx_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 1);

    auto ta_xxxx_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 2);

    auto ta_xxxx_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 3);

    auto ta_xxxx_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 4);

    auto ta_xxxx_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 5);

    auto ta_xxxx_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 6);

    auto ta_xxxx_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 7);

    auto ta_xxxx_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 8);

    auto ta_xxxx_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 9);

    auto ta_xxxx_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 10);

    auto ta_xxxx_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 11);

    auto ta_xxxx_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 12);

    auto ta_xxxx_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 13);

    auto ta_xxxx_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 14);

    auto ta_xxxx_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 15);

    auto ta_xxxx_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 16);

    auto ta_xxxx_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 17);

    auto ta_xxxx_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 18);

    auto ta_xxxx_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 19);

    auto ta_xxxx_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 20);

    auto ta_xxxy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 21);

    auto ta_xxxy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 23);

    auto ta_xxxy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 26);

    auto ta_xxxy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 30);

    auto ta_xxxy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 35);

    auto ta_xxxy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 36);

    auto ta_xxxy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 37);

    auto ta_xxxy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 38);

    auto ta_xxxy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 39);

    auto ta_xxxy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 40);

    auto ta_xxxz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 42);

    auto ta_xxxz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 43);

    auto ta_xxxz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 44);

    auto ta_xxxz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 45);

    auto ta_xxxz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 47);

    auto ta_xxxz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 48);

    auto ta_xxxz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 51);

    auto ta_xxxz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 52);

    auto ta_xxxz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 56);

    auto ta_xxxz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 58);

    auto ta_xxxz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 59);

    auto ta_xxxz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 60);

    auto ta_xxxz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 61);

    auto ta_xxxz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 62);

    auto ta_xxyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 63);

    auto ta_xxyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 64);

    auto ta_xxyy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 65);

    auto ta_xxyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 66);

    auto ta_xxyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 67);

    auto ta_xxyy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 68);

    auto ta_xxyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 69);

    auto ta_xxyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 70);

    auto ta_xxyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 71);

    auto ta_xxyy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 72);

    auto ta_xxyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 73);

    auto ta_xxyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 74);

    auto ta_xxyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 75);

    auto ta_xxyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 76);

    auto ta_xxyy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 77);

    auto ta_xxyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 78);

    auto ta_xxyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 79);

    auto ta_xxyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 80);

    auto ta_xxyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 81);

    auto ta_xxyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 82);

    auto ta_xxyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 83);

    auto ta_xxyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 86);

    auto ta_xxyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 89);

    auto ta_xxyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 93);

    auto ta_xxyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 98);

    auto ta_xxyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 100);

    auto ta_xxyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 101);

    auto ta_xxyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 102);

    auto ta_xxyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 103);

    auto ta_xxzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 105);

    auto ta_xxzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 106);

    auto ta_xxzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 107);

    auto ta_xxzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 108);

    auto ta_xxzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 109);

    auto ta_xxzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 110);

    auto ta_xxzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 111);

    auto ta_xxzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 112);

    auto ta_xxzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 113);

    auto ta_xxzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 114);

    auto ta_xxzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 115);

    auto ta_xxzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 116);

    auto ta_xxzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 117);

    auto ta_xxzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 118);

    auto ta_xxzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 119);

    auto ta_xxzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 120);

    auto ta_xxzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 121);

    auto ta_xxzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 122);

    auto ta_xxzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 123);

    auto ta_xxzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 124);

    auto ta_xxzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 125);

    auto ta_xyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 127);

    auto ta_xyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 129);

    auto ta_xyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 130);

    auto ta_xyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 132);

    auto ta_xyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 133);

    auto ta_xyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 134);

    auto ta_xyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 136);

    auto ta_xyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 137);

    auto ta_xyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 138);

    auto ta_xyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 139);

    auto ta_xyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 141);

    auto ta_xyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 142);

    auto ta_xyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 143);

    auto ta_xyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 144);

    auto ta_xyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 145);

    auto ta_xyyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 146);

    auto ta_xyyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 163);

    auto ta_xyyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 164);

    auto ta_xyyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 165);

    auto ta_xyyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 166);

    auto ta_xyyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 167);

    auto ta_xyzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 183);

    auto ta_xyzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 184);

    auto ta_xyzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 185);

    auto ta_xyzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 186);

    auto ta_xyzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 187);

    auto ta_xzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 191);

    auto ta_xzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 193);

    auto ta_xzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 194);

    auto ta_xzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 196);

    auto ta_xzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 197);

    auto ta_xzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 198);

    auto ta_xzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 200);

    auto ta_xzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 201);

    auto ta_xzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 202);

    auto ta_xzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 203);

    auto ta_xzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 204);

    auto ta_xzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 205);

    auto ta_xzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 206);

    auto ta_xzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 207);

    auto ta_xzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 208);

    auto ta_xzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 209);

    auto ta_yyyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 210);

    auto ta_yyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 211);

    auto ta_yyyy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 212);

    auto ta_yyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 213);

    auto ta_yyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 214);

    auto ta_yyyy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 215);

    auto ta_yyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 216);

    auto ta_yyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 217);

    auto ta_yyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 218);

    auto ta_yyyy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 219);

    auto ta_yyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 220);

    auto ta_yyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 221);

    auto ta_yyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 222);

    auto ta_yyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 223);

    auto ta_yyyy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 224);

    auto ta_yyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 225);

    auto ta_yyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 226);

    auto ta_yyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 227);

    auto ta_yyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 228);

    auto ta_yyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 229);

    auto ta_yyyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 230);

    auto ta_yyyz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 232);

    auto ta_yyyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 233);

    auto ta_yyyz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 234);

    auto ta_yyyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 236);

    auto ta_yyyz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 237);

    auto ta_yyyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 240);

    auto ta_yyyz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 241);

    auto ta_yyyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 245);

    auto ta_yyyz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 246);

    auto ta_yyyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 247);

    auto ta_yyyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 248);

    auto ta_yyyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 249);

    auto ta_yyyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 250);

    auto ta_yyyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 251);

    auto ta_yyzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 252);

    auto ta_yyzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 253);

    auto ta_yyzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 254);

    auto ta_yyzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 255);

    auto ta_yyzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 256);

    auto ta_yyzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 257);

    auto ta_yyzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 258);

    auto ta_yyzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 259);

    auto ta_yyzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 260);

    auto ta_yyzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 261);

    auto ta_yyzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 262);

    auto ta_yyzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 263);

    auto ta_yyzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 264);

    auto ta_yyzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 265);

    auto ta_yyzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 266);

    auto ta_yyzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 267);

    auto ta_yyzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 268);

    auto ta_yyzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 269);

    auto ta_yyzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 270);

    auto ta_yyzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 271);

    auto ta_yyzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 272);

    auto ta_yzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 273);

    auto ta_yzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 275);

    auto ta_yzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 277);

    auto ta_yzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 278);

    auto ta_yzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 280);

    auto ta_yzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 281);

    auto ta_yzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 282);

    auto ta_yzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 284);

    auto ta_yzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 285);

    auto ta_yzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 286);

    auto ta_yzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 287);

    auto ta_yzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 288);

    auto ta_yzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 289);

    auto ta_yzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 290);

    auto ta_yzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 291);

    auto ta_yzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 292);

    auto ta_yzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 293);

    auto ta_zzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 294);

    auto ta_zzzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 295);

    auto ta_zzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 296);

    auto ta_zzzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 297);

    auto ta_zzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 298);

    auto ta_zzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 299);

    auto ta_zzzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 300);

    auto ta_zzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 301);

    auto ta_zzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 302);

    auto ta_zzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 303);

    auto ta_zzzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 304);

    auto ta_zzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 305);

    auto ta_zzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 306);

    auto ta_zzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 307);

    auto ta_zzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 308);

    auto ta_zzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 309);

    auto ta_zzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 310);

    auto ta_zzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 311);

    auto ta_zzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 312);

    auto ta_zzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 313);

    auto ta_zzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 314);

    // Set up components of auxiliary buffer : GH

    auto ta_xxxx_xxxxx_1 = pbuffer.data(idx_npot_1_gh);

    auto ta_xxxx_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 1);

    auto ta_xxxx_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 2);

    auto ta_xxxx_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 3);

    auto ta_xxxx_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 4);

    auto ta_xxxx_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 5);

    auto ta_xxxx_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 6);

    auto ta_xxxx_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 7);

    auto ta_xxxx_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 8);

    auto ta_xxxx_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 9);

    auto ta_xxxx_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 10);

    auto ta_xxxx_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 11);

    auto ta_xxxx_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 12);

    auto ta_xxxx_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 13);

    auto ta_xxxx_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 14);

    auto ta_xxxx_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 15);

    auto ta_xxxx_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 16);

    auto ta_xxxx_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 17);

    auto ta_xxxx_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 18);

    auto ta_xxxx_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 19);

    auto ta_xxxx_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 20);

    auto ta_xxxy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 21);

    auto ta_xxxy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 23);

    auto ta_xxxy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 26);

    auto ta_xxxy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 30);

    auto ta_xxxy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 35);

    auto ta_xxxy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 36);

    auto ta_xxxy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 37);

    auto ta_xxxy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 38);

    auto ta_xxxy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 39);

    auto ta_xxxy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 40);

    auto ta_xxxz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 42);

    auto ta_xxxz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 43);

    auto ta_xxxz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 44);

    auto ta_xxxz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 45);

    auto ta_xxxz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 47);

    auto ta_xxxz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 48);

    auto ta_xxxz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 51);

    auto ta_xxxz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 52);

    auto ta_xxxz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 56);

    auto ta_xxxz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 58);

    auto ta_xxxz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 59);

    auto ta_xxxz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 60);

    auto ta_xxxz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 61);

    auto ta_xxxz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 62);

    auto ta_xxyy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 63);

    auto ta_xxyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 64);

    auto ta_xxyy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 65);

    auto ta_xxyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 66);

    auto ta_xxyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 67);

    auto ta_xxyy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 68);

    auto ta_xxyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 69);

    auto ta_xxyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 70);

    auto ta_xxyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 71);

    auto ta_xxyy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 72);

    auto ta_xxyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 73);

    auto ta_xxyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 74);

    auto ta_xxyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 75);

    auto ta_xxyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 76);

    auto ta_xxyy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 77);

    auto ta_xxyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 78);

    auto ta_xxyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 79);

    auto ta_xxyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 80);

    auto ta_xxyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 81);

    auto ta_xxyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 82);

    auto ta_xxyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 83);

    auto ta_xxyz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 86);

    auto ta_xxyz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 89);

    auto ta_xxyz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 93);

    auto ta_xxyz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 98);

    auto ta_xxyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 100);

    auto ta_xxyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 101);

    auto ta_xxyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 102);

    auto ta_xxyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 103);

    auto ta_xxzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 105);

    auto ta_xxzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 106);

    auto ta_xxzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 107);

    auto ta_xxzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 108);

    auto ta_xxzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 109);

    auto ta_xxzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 110);

    auto ta_xxzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 111);

    auto ta_xxzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 112);

    auto ta_xxzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 113);

    auto ta_xxzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 114);

    auto ta_xxzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 115);

    auto ta_xxzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 116);

    auto ta_xxzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 117);

    auto ta_xxzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 118);

    auto ta_xxzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 119);

    auto ta_xxzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 120);

    auto ta_xxzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 121);

    auto ta_xxzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 122);

    auto ta_xxzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 123);

    auto ta_xxzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 124);

    auto ta_xxzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 125);

    auto ta_xyyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 127);

    auto ta_xyyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 129);

    auto ta_xyyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 130);

    auto ta_xyyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 132);

    auto ta_xyyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 133);

    auto ta_xyyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 134);

    auto ta_xyyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 136);

    auto ta_xyyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 137);

    auto ta_xyyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 138);

    auto ta_xyyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 139);

    auto ta_xyyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 141);

    auto ta_xyyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 142);

    auto ta_xyyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 143);

    auto ta_xyyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 144);

    auto ta_xyyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 145);

    auto ta_xyyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 146);

    auto ta_xyyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 163);

    auto ta_xyyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 164);

    auto ta_xyyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 165);

    auto ta_xyyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 166);

    auto ta_xyyz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 167);

    auto ta_xyzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 183);

    auto ta_xyzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 184);

    auto ta_xyzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 185);

    auto ta_xyzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 186);

    auto ta_xyzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 187);

    auto ta_xzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 191);

    auto ta_xzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 193);

    auto ta_xzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 194);

    auto ta_xzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 196);

    auto ta_xzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 197);

    auto ta_xzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 198);

    auto ta_xzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 200);

    auto ta_xzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 201);

    auto ta_xzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 202);

    auto ta_xzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 203);

    auto ta_xzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 204);

    auto ta_xzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 205);

    auto ta_xzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 206);

    auto ta_xzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 207);

    auto ta_xzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 208);

    auto ta_xzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 209);

    auto ta_yyyy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 210);

    auto ta_yyyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 211);

    auto ta_yyyy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 212);

    auto ta_yyyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 213);

    auto ta_yyyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 214);

    auto ta_yyyy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 215);

    auto ta_yyyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 216);

    auto ta_yyyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 217);

    auto ta_yyyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 218);

    auto ta_yyyy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 219);

    auto ta_yyyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 220);

    auto ta_yyyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 221);

    auto ta_yyyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 222);

    auto ta_yyyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 223);

    auto ta_yyyy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 224);

    auto ta_yyyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 225);

    auto ta_yyyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 226);

    auto ta_yyyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 227);

    auto ta_yyyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 228);

    auto ta_yyyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 229);

    auto ta_yyyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 230);

    auto ta_yyyz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 232);

    auto ta_yyyz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 233);

    auto ta_yyyz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 234);

    auto ta_yyyz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 236);

    auto ta_yyyz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 237);

    auto ta_yyyz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 240);

    auto ta_yyyz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 241);

    auto ta_yyyz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 245);

    auto ta_yyyz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 246);

    auto ta_yyyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 247);

    auto ta_yyyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 248);

    auto ta_yyyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 249);

    auto ta_yyyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 250);

    auto ta_yyyz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 251);

    auto ta_yyzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 252);

    auto ta_yyzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 253);

    auto ta_yyzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 254);

    auto ta_yyzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 255);

    auto ta_yyzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 256);

    auto ta_yyzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 257);

    auto ta_yyzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 258);

    auto ta_yyzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 259);

    auto ta_yyzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 260);

    auto ta_yyzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 261);

    auto ta_yyzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 262);

    auto ta_yyzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 263);

    auto ta_yyzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 264);

    auto ta_yyzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 265);

    auto ta_yyzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 266);

    auto ta_yyzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 267);

    auto ta_yyzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 268);

    auto ta_yyzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 269);

    auto ta_yyzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 270);

    auto ta_yyzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 271);

    auto ta_yyzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 272);

    auto ta_yzzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 273);

    auto ta_yzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 275);

    auto ta_yzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 277);

    auto ta_yzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 278);

    auto ta_yzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 280);

    auto ta_yzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 281);

    auto ta_yzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 282);

    auto ta_yzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 284);

    auto ta_yzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 285);

    auto ta_yzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 286);

    auto ta_yzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 287);

    auto ta_yzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 288);

    auto ta_yzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 289);

    auto ta_yzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 290);

    auto ta_yzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 291);

    auto ta_yzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 292);

    auto ta_yzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 293);

    auto ta_zzzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 294);

    auto ta_zzzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 295);

    auto ta_zzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 296);

    auto ta_zzzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 297);

    auto ta_zzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 298);

    auto ta_zzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 299);

    auto ta_zzzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 300);

    auto ta_zzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 301);

    auto ta_zzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 302);

    auto ta_zzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 303);

    auto ta_zzzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 304);

    auto ta_zzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 305);

    auto ta_zzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 306);

    auto ta_zzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 307);

    auto ta_zzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 308);

    auto ta_zzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 309);

    auto ta_zzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 310);

    auto ta_zzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 311);

    auto ta_zzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 312);

    auto ta_zzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 313);

    auto ta_zzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 314);

    // Set up components of auxiliary buffer : HG

    auto ta_xxxxx_xxxx_0 = pbuffer.data(idx_npot_0_hg);

    auto ta_xxxxx_xxxy_0 = pbuffer.data(idx_npot_0_hg + 1);

    auto ta_xxxxx_xxxz_0 = pbuffer.data(idx_npot_0_hg + 2);

    auto ta_xxxxx_xxyy_0 = pbuffer.data(idx_npot_0_hg + 3);

    auto ta_xxxxx_xxyz_0 = pbuffer.data(idx_npot_0_hg + 4);

    auto ta_xxxxx_xxzz_0 = pbuffer.data(idx_npot_0_hg + 5);

    auto ta_xxxxx_xyyy_0 = pbuffer.data(idx_npot_0_hg + 6);

    auto ta_xxxxx_xyyz_0 = pbuffer.data(idx_npot_0_hg + 7);

    auto ta_xxxxx_xyzz_0 = pbuffer.data(idx_npot_0_hg + 8);

    auto ta_xxxxx_xzzz_0 = pbuffer.data(idx_npot_0_hg + 9);

    auto ta_xxxxx_yyyy_0 = pbuffer.data(idx_npot_0_hg + 10);

    auto ta_xxxxx_yyyz_0 = pbuffer.data(idx_npot_0_hg + 11);

    auto ta_xxxxx_yyzz_0 = pbuffer.data(idx_npot_0_hg + 12);

    auto ta_xxxxx_yzzz_0 = pbuffer.data(idx_npot_0_hg + 13);

    auto ta_xxxxx_zzzz_0 = pbuffer.data(idx_npot_0_hg + 14);

    auto ta_xxxxz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 32);

    auto ta_xxxxz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 34);

    auto ta_xxxxz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 35);

    auto ta_xxxxz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 37);

    auto ta_xxxxz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 38);

    auto ta_xxxxz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 39);

    auto ta_xxxyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 46);

    auto ta_xxxyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 48);

    auto ta_xxxyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 49);

    auto ta_xxxyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 51);

    auto ta_xxxyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 52);

    auto ta_xxxyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 53);

    auto ta_xxxyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 55);

    auto ta_xxxyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 56);

    auto ta_xxxyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 57);

    auto ta_xxxyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 58);

    auto ta_xxxzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 75);

    auto ta_xxxzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 76);

    auto ta_xxxzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 77);

    auto ta_xxxzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 78);

    auto ta_xxxzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 79);

    auto ta_xxxzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 80);

    auto ta_xxxzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 81);

    auto ta_xxxzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 82);

    auto ta_xxxzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 83);

    auto ta_xxxzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 84);

    auto ta_xxxzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 86);

    auto ta_xxxzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 87);

    auto ta_xxxzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 88);

    auto ta_xxxzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 89);

    auto ta_xxyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 91);

    auto ta_xxyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 93);

    auto ta_xxyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 94);

    auto ta_xxyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 96);

    auto ta_xxyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 97);

    auto ta_xxyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 98);

    auto ta_xxyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 100);

    auto ta_xxyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 101);

    auto ta_xxyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 102);

    auto ta_xxyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 103);

    auto ta_xxzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 135);

    auto ta_xxzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 136);

    auto ta_xxzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 137);

    auto ta_xxzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 138);

    auto ta_xxzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 139);

    auto ta_xxzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 140);

    auto ta_xxzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 141);

    auto ta_xxzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 142);

    auto ta_xxzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 143);

    auto ta_xxzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 144);

    auto ta_xxzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 146);

    auto ta_xxzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 147);

    auto ta_xxzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 148);

    auto ta_xxzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 149);

    auto ta_xyyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 151);

    auto ta_xyyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 153);

    auto ta_xyyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 154);

    auto ta_xyyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 156);

    auto ta_xyyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 157);

    auto ta_xyyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 158);

    auto ta_xyyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 160);

    auto ta_xyyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 161);

    auto ta_xyyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 162);

    auto ta_xyyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 163);

    auto ta_xyyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 184);

    auto ta_xyyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 187);

    auto ta_xyyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 188);

    auto ta_xyyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 191);

    auto ta_xyyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 192);

    auto ta_xyyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 193);

    auto ta_xzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 212);

    auto ta_xzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 214);

    auto ta_xzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 215);

    auto ta_xzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 217);

    auto ta_xzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 218);

    auto ta_xzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 219);

    auto ta_xzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 221);

    auto ta_xzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 222);

    auto ta_xzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 223);

    auto ta_xzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 224);

    auto ta_yyyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 225);

    auto ta_yyyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 226);

    auto ta_yyyyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 227);

    auto ta_yyyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 228);

    auto ta_yyyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 229);

    auto ta_yyyyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 230);

    auto ta_yyyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 231);

    auto ta_yyyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 232);

    auto ta_yyyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 233);

    auto ta_yyyyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 234);

    auto ta_yyyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 235);

    auto ta_yyyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 236);

    auto ta_yyyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 237);

    auto ta_yyyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 238);

    auto ta_yyyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 239);

    auto ta_yyyyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 242);

    auto ta_yyyyz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 244);

    auto ta_yyyyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 245);

    auto ta_yyyyz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 247);

    auto ta_yyyyz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 248);

    auto ta_yyyyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 249);

    auto ta_yyyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 251);

    auto ta_yyyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 252);

    auto ta_yyyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 253);

    auto ta_yyyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 254);

    auto ta_yyyzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 255);

    auto ta_yyyzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 256);

    auto ta_yyyzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 257);

    auto ta_yyyzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 258);

    auto ta_yyyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 259);

    auto ta_yyyzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 260);

    auto ta_yyyzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 261);

    auto ta_yyyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 262);

    auto ta_yyyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 263);

    auto ta_yyyzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 264);

    auto ta_yyyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 265);

    auto ta_yyyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 266);

    auto ta_yyyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 267);

    auto ta_yyyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 268);

    auto ta_yyyzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 269);

    auto ta_yyzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 270);

    auto ta_yyzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 271);

    auto ta_yyzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 272);

    auto ta_yyzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 273);

    auto ta_yyzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 274);

    auto ta_yyzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 275);

    auto ta_yyzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 276);

    auto ta_yyzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 277);

    auto ta_yyzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 278);

    auto ta_yyzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 279);

    auto ta_yyzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 280);

    auto ta_yyzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 281);

    auto ta_yyzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 282);

    auto ta_yyzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 283);

    auto ta_yyzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 284);

    auto ta_yzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 286);

    auto ta_yzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 287);

    auto ta_yzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 288);

    auto ta_yzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 289);

    auto ta_yzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 290);

    auto ta_yzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 291);

    auto ta_yzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 292);

    auto ta_yzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 293);

    auto ta_yzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 294);

    auto ta_yzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 295);

    auto ta_yzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 296);

    auto ta_yzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 297);

    auto ta_yzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 298);

    auto ta_yzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 299);

    auto ta_zzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 300);

    auto ta_zzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 301);

    auto ta_zzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 302);

    auto ta_zzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 303);

    auto ta_zzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 304);

    auto ta_zzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 305);

    auto ta_zzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 306);

    auto ta_zzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 307);

    auto ta_zzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 308);

    auto ta_zzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 309);

    auto ta_zzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 310);

    auto ta_zzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 311);

    auto ta_zzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 312);

    auto ta_zzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 313);

    auto ta_zzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 314);

    // Set up components of auxiliary buffer : HG

    auto ta_xxxxx_xxxx_1 = pbuffer.data(idx_npot_1_hg);

    auto ta_xxxxx_xxxy_1 = pbuffer.data(idx_npot_1_hg + 1);

    auto ta_xxxxx_xxxz_1 = pbuffer.data(idx_npot_1_hg + 2);

    auto ta_xxxxx_xxyy_1 = pbuffer.data(idx_npot_1_hg + 3);

    auto ta_xxxxx_xxyz_1 = pbuffer.data(idx_npot_1_hg + 4);

    auto ta_xxxxx_xxzz_1 = pbuffer.data(idx_npot_1_hg + 5);

    auto ta_xxxxx_xyyy_1 = pbuffer.data(idx_npot_1_hg + 6);

    auto ta_xxxxx_xyyz_1 = pbuffer.data(idx_npot_1_hg + 7);

    auto ta_xxxxx_xyzz_1 = pbuffer.data(idx_npot_1_hg + 8);

    auto ta_xxxxx_xzzz_1 = pbuffer.data(idx_npot_1_hg + 9);

    auto ta_xxxxx_yyyy_1 = pbuffer.data(idx_npot_1_hg + 10);

    auto ta_xxxxx_yyyz_1 = pbuffer.data(idx_npot_1_hg + 11);

    auto ta_xxxxx_yyzz_1 = pbuffer.data(idx_npot_1_hg + 12);

    auto ta_xxxxx_yzzz_1 = pbuffer.data(idx_npot_1_hg + 13);

    auto ta_xxxxx_zzzz_1 = pbuffer.data(idx_npot_1_hg + 14);

    auto ta_xxxxz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 32);

    auto ta_xxxxz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 34);

    auto ta_xxxxz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 35);

    auto ta_xxxxz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 37);

    auto ta_xxxxz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 38);

    auto ta_xxxxz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 39);

    auto ta_xxxyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 46);

    auto ta_xxxyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 48);

    auto ta_xxxyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 49);

    auto ta_xxxyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 51);

    auto ta_xxxyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 52);

    auto ta_xxxyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 53);

    auto ta_xxxyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 55);

    auto ta_xxxyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 56);

    auto ta_xxxyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 57);

    auto ta_xxxyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 58);

    auto ta_xxxzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 75);

    auto ta_xxxzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 76);

    auto ta_xxxzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 77);

    auto ta_xxxzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 78);

    auto ta_xxxzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 79);

    auto ta_xxxzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 80);

    auto ta_xxxzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 81);

    auto ta_xxxzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 82);

    auto ta_xxxzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 83);

    auto ta_xxxzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 84);

    auto ta_xxxzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 86);

    auto ta_xxxzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 87);

    auto ta_xxxzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 88);

    auto ta_xxxzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 89);

    auto ta_xxyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 91);

    auto ta_xxyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 93);

    auto ta_xxyyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 94);

    auto ta_xxyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 96);

    auto ta_xxyyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 97);

    auto ta_xxyyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 98);

    auto ta_xxyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 100);

    auto ta_xxyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 101);

    auto ta_xxyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 102);

    auto ta_xxyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 103);

    auto ta_xxzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 135);

    auto ta_xxzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 136);

    auto ta_xxzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 137);

    auto ta_xxzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 138);

    auto ta_xxzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 139);

    auto ta_xxzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 140);

    auto ta_xxzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 141);

    auto ta_xxzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 142);

    auto ta_xxzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 143);

    auto ta_xxzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 144);

    auto ta_xxzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 146);

    auto ta_xxzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 147);

    auto ta_xxzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 148);

    auto ta_xxzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 149);

    auto ta_xyyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 151);

    auto ta_xyyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 153);

    auto ta_xyyyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 154);

    auto ta_xyyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 156);

    auto ta_xyyyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 157);

    auto ta_xyyyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 158);

    auto ta_xyyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 160);

    auto ta_xyyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 161);

    auto ta_xyyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 162);

    auto ta_xyyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 163);

    auto ta_xyyzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 184);

    auto ta_xyyzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 187);

    auto ta_xyyzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 188);

    auto ta_xyyzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 191);

    auto ta_xyyzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 192);

    auto ta_xyyzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 193);

    auto ta_xzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 212);

    auto ta_xzzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 214);

    auto ta_xzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 215);

    auto ta_xzzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 217);

    auto ta_xzzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 218);

    auto ta_xzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 219);

    auto ta_xzzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 221);

    auto ta_xzzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 222);

    auto ta_xzzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 223);

    auto ta_xzzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 224);

    auto ta_yyyyy_xxxx_1 = pbuffer.data(idx_npot_1_hg + 225);

    auto ta_yyyyy_xxxy_1 = pbuffer.data(idx_npot_1_hg + 226);

    auto ta_yyyyy_xxxz_1 = pbuffer.data(idx_npot_1_hg + 227);

    auto ta_yyyyy_xxyy_1 = pbuffer.data(idx_npot_1_hg + 228);

    auto ta_yyyyy_xxyz_1 = pbuffer.data(idx_npot_1_hg + 229);

    auto ta_yyyyy_xxzz_1 = pbuffer.data(idx_npot_1_hg + 230);

    auto ta_yyyyy_xyyy_1 = pbuffer.data(idx_npot_1_hg + 231);

    auto ta_yyyyy_xyyz_1 = pbuffer.data(idx_npot_1_hg + 232);

    auto ta_yyyyy_xyzz_1 = pbuffer.data(idx_npot_1_hg + 233);

    auto ta_yyyyy_xzzz_1 = pbuffer.data(idx_npot_1_hg + 234);

    auto ta_yyyyy_yyyy_1 = pbuffer.data(idx_npot_1_hg + 235);

    auto ta_yyyyy_yyyz_1 = pbuffer.data(idx_npot_1_hg + 236);

    auto ta_yyyyy_yyzz_1 = pbuffer.data(idx_npot_1_hg + 237);

    auto ta_yyyyy_yzzz_1 = pbuffer.data(idx_npot_1_hg + 238);

    auto ta_yyyyy_zzzz_1 = pbuffer.data(idx_npot_1_hg + 239);

    auto ta_yyyyz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 242);

    auto ta_yyyyz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 244);

    auto ta_yyyyz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 245);

    auto ta_yyyyz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 247);

    auto ta_yyyyz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 248);

    auto ta_yyyyz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 249);

    auto ta_yyyyz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 251);

    auto ta_yyyyz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 252);

    auto ta_yyyyz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 253);

    auto ta_yyyyz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 254);

    auto ta_yyyzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 255);

    auto ta_yyyzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 256);

    auto ta_yyyzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 257);

    auto ta_yyyzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 258);

    auto ta_yyyzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 259);

    auto ta_yyyzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 260);

    auto ta_yyyzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 261);

    auto ta_yyyzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 262);

    auto ta_yyyzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 263);

    auto ta_yyyzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 264);

    auto ta_yyyzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 265);

    auto ta_yyyzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 266);

    auto ta_yyyzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 267);

    auto ta_yyyzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 268);

    auto ta_yyyzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 269);

    auto ta_yyzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 270);

    auto ta_yyzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 271);

    auto ta_yyzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 272);

    auto ta_yyzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 273);

    auto ta_yyzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 274);

    auto ta_yyzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 275);

    auto ta_yyzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 276);

    auto ta_yyzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 277);

    auto ta_yyzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 278);

    auto ta_yyzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 279);

    auto ta_yyzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 280);

    auto ta_yyzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 281);

    auto ta_yyzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 282);

    auto ta_yyzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 283);

    auto ta_yyzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 284);

    auto ta_yzzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 286);

    auto ta_yzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 287);

    auto ta_yzzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 288);

    auto ta_yzzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 289);

    auto ta_yzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 290);

    auto ta_yzzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 291);

    auto ta_yzzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 292);

    auto ta_yzzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 293);

    auto ta_yzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 294);

    auto ta_yzzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 295);

    auto ta_yzzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 296);

    auto ta_yzzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 297);

    auto ta_yzzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 298);

    auto ta_yzzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 299);

    auto ta_zzzzz_xxxx_1 = pbuffer.data(idx_npot_1_hg + 300);

    auto ta_zzzzz_xxxy_1 = pbuffer.data(idx_npot_1_hg + 301);

    auto ta_zzzzz_xxxz_1 = pbuffer.data(idx_npot_1_hg + 302);

    auto ta_zzzzz_xxyy_1 = pbuffer.data(idx_npot_1_hg + 303);

    auto ta_zzzzz_xxyz_1 = pbuffer.data(idx_npot_1_hg + 304);

    auto ta_zzzzz_xxzz_1 = pbuffer.data(idx_npot_1_hg + 305);

    auto ta_zzzzz_xyyy_1 = pbuffer.data(idx_npot_1_hg + 306);

    auto ta_zzzzz_xyyz_1 = pbuffer.data(idx_npot_1_hg + 307);

    auto ta_zzzzz_xyzz_1 = pbuffer.data(idx_npot_1_hg + 308);

    auto ta_zzzzz_xzzz_1 = pbuffer.data(idx_npot_1_hg + 309);

    auto ta_zzzzz_yyyy_1 = pbuffer.data(idx_npot_1_hg + 310);

    auto ta_zzzzz_yyyz_1 = pbuffer.data(idx_npot_1_hg + 311);

    auto ta_zzzzz_yyzz_1 = pbuffer.data(idx_npot_1_hg + 312);

    auto ta_zzzzz_yzzz_1 = pbuffer.data(idx_npot_1_hg + 313);

    auto ta_zzzzz_zzzz_1 = pbuffer.data(idx_npot_1_hg + 314);

    // Set up components of auxiliary buffer : HH

    auto ta_xxxxx_xxxxx_0 = pbuffer.data(idx_npot_0_hh);

    auto ta_xxxxx_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 1);

    auto ta_xxxxx_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 2);

    auto ta_xxxxx_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 3);

    auto ta_xxxxx_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 4);

    auto ta_xxxxx_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 5);

    auto ta_xxxxx_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 6);

    auto ta_xxxxx_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 7);

    auto ta_xxxxx_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 8);

    auto ta_xxxxx_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 9);

    auto ta_xxxxx_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 10);

    auto ta_xxxxx_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 11);

    auto ta_xxxxx_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 12);

    auto ta_xxxxx_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 13);

    auto ta_xxxxx_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 14);

    auto ta_xxxxx_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 15);

    auto ta_xxxxx_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 16);

    auto ta_xxxxx_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 17);

    auto ta_xxxxx_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 18);

    auto ta_xxxxx_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 19);

    auto ta_xxxxx_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 20);

    auto ta_xxxxy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 21);

    auto ta_xxxxy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 22);

    auto ta_xxxxy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 23);

    auto ta_xxxxy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 24);

    auto ta_xxxxy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 26);

    auto ta_xxxxy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 27);

    auto ta_xxxxy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 30);

    auto ta_xxxxy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 31);

    auto ta_xxxxy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 35);

    auto ta_xxxxy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 36);

    auto ta_xxxxy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 37);

    auto ta_xxxxy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 38);

    auto ta_xxxxy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 39);

    auto ta_xxxxy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 40);

    auto ta_xxxxz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 42);

    auto ta_xxxxz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 43);

    auto ta_xxxxz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 44);

    auto ta_xxxxz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 45);

    auto ta_xxxxz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 46);

    auto ta_xxxxz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 47);

    auto ta_xxxxz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 48);

    auto ta_xxxxz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 49);

    auto ta_xxxxz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 50);

    auto ta_xxxxz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 51);

    auto ta_xxxxz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 52);

    auto ta_xxxxz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 53);

    auto ta_xxxxz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 54);

    auto ta_xxxxz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 55);

    auto ta_xxxxz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 56);

    auto ta_xxxxz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 58);

    auto ta_xxxxz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 59);

    auto ta_xxxxz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 60);

    auto ta_xxxxz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 61);

    auto ta_xxxxz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 62);

    auto ta_xxxyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 63);

    auto ta_xxxyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 64);

    auto ta_xxxyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 65);

    auto ta_xxxyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 66);

    auto ta_xxxyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 67);

    auto ta_xxxyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 68);

    auto ta_xxxyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 69);

    auto ta_xxxyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 70);

    auto ta_xxxyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 71);

    auto ta_xxxyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 72);

    auto ta_xxxyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 73);

    auto ta_xxxyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 74);

    auto ta_xxxyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 75);

    auto ta_xxxyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 76);

    auto ta_xxxyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 77);

    auto ta_xxxyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 78);

    auto ta_xxxyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 79);

    auto ta_xxxyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 80);

    auto ta_xxxyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 81);

    auto ta_xxxyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 82);

    auto ta_xxxyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 83);

    auto ta_xxxyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 86);

    auto ta_xxxyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 89);

    auto ta_xxxyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 93);

    auto ta_xxxyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 98);

    auto ta_xxxyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 100);

    auto ta_xxxyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 101);

    auto ta_xxxyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 102);

    auto ta_xxxyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 103);

    auto ta_xxxzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 105);

    auto ta_xxxzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 106);

    auto ta_xxxzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 107);

    auto ta_xxxzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 108);

    auto ta_xxxzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 109);

    auto ta_xxxzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 110);

    auto ta_xxxzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 111);

    auto ta_xxxzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 112);

    auto ta_xxxzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 113);

    auto ta_xxxzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 114);

    auto ta_xxxzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 115);

    auto ta_xxxzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 116);

    auto ta_xxxzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 117);

    auto ta_xxxzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 118);

    auto ta_xxxzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 119);

    auto ta_xxxzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 120);

    auto ta_xxxzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 121);

    auto ta_xxxzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 122);

    auto ta_xxxzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 123);

    auto ta_xxxzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 124);

    auto ta_xxxzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 125);

    auto ta_xxyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 126);

    auto ta_xxyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 127);

    auto ta_xxyyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 128);

    auto ta_xxyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 129);

    auto ta_xxyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 130);

    auto ta_xxyyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 131);

    auto ta_xxyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 132);

    auto ta_xxyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 133);

    auto ta_xxyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 134);

    auto ta_xxyyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 135);

    auto ta_xxyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 136);

    auto ta_xxyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 137);

    auto ta_xxyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 138);

    auto ta_xxyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 139);

    auto ta_xxyyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 140);

    auto ta_xxyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 141);

    auto ta_xxyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 142);

    auto ta_xxyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 143);

    auto ta_xxyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 144);

    auto ta_xxyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 145);

    auto ta_xxyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 146);

    auto ta_xxyyz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 148);

    auto ta_xxyyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 149);

    auto ta_xxyyz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 150);

    auto ta_xxyyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 152);

    auto ta_xxyyz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 153);

    auto ta_xxyyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 156);

    auto ta_xxyyz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 157);

    auto ta_xxyyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 161);

    auto ta_xxyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 163);

    auto ta_xxyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 164);

    auto ta_xxyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 165);

    auto ta_xxyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 166);

    auto ta_xxyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 167);

    auto ta_xxyzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 168);

    auto ta_xxyzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 170);

    auto ta_xxyzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 173);

    auto ta_xxyzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 177);

    auto ta_xxyzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 182);

    auto ta_xxyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 183);

    auto ta_xxyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 184);

    auto ta_xxyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 185);

    auto ta_xxyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 186);

    auto ta_xxyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 187);

    auto ta_xxzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 189);

    auto ta_xxzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 190);

    auto ta_xxzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 191);

    auto ta_xxzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 192);

    auto ta_xxzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 193);

    auto ta_xxzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 194);

    auto ta_xxzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 195);

    auto ta_xxzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 196);

    auto ta_xxzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 197);

    auto ta_xxzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 198);

    auto ta_xxzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 199);

    auto ta_xxzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 200);

    auto ta_xxzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 201);

    auto ta_xxzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 202);

    auto ta_xxzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 203);

    auto ta_xxzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 204);

    auto ta_xxzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 205);

    auto ta_xxzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 206);

    auto ta_xxzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 207);

    auto ta_xxzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 208);

    auto ta_xxzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 209);

    auto ta_xyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 210);

    auto ta_xyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 211);

    auto ta_xyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 213);

    auto ta_xyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 214);

    auto ta_xyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 216);

    auto ta_xyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 217);

    auto ta_xyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 218);

    auto ta_xyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 220);

    auto ta_xyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 221);

    auto ta_xyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 222);

    auto ta_xyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 223);

    auto ta_xyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 225);

    auto ta_xyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 226);

    auto ta_xyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 227);

    auto ta_xyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 228);

    auto ta_xyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 229);

    auto ta_xyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 230);

    auto ta_xyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 247);

    auto ta_xyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 248);

    auto ta_xyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 249);

    auto ta_xyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 250);

    auto ta_xyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 251);

    auto ta_xyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 256);

    auto ta_xyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 259);

    auto ta_xyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 260);

    auto ta_xyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 263);

    auto ta_xyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 264);

    auto ta_xyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 265);

    auto ta_xyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 267);

    auto ta_xyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 268);

    auto ta_xyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 269);

    auto ta_xyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 270);

    auto ta_xyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 271);

    auto ta_xyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 272);

    auto ta_xyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 288);

    auto ta_xyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 289);

    auto ta_xyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 290);

    auto ta_xyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 291);

    auto ta_xyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 292);

    auto ta_xzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 294);

    auto ta_xzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 296);

    auto ta_xzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 298);

    auto ta_xzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 299);

    auto ta_xzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 301);

    auto ta_xzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 302);

    auto ta_xzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 303);

    auto ta_xzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 305);

    auto ta_xzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 306);

    auto ta_xzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 307);

    auto ta_xzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 308);

    auto ta_xzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 309);

    auto ta_xzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 310);

    auto ta_xzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 311);

    auto ta_xzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 312);

    auto ta_xzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 313);

    auto ta_xzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 314);

    auto ta_yyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 315);

    auto ta_yyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 316);

    auto ta_yyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 317);

    auto ta_yyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 318);

    auto ta_yyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 319);

    auto ta_yyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 320);

    auto ta_yyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 321);

    auto ta_yyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 322);

    auto ta_yyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 323);

    auto ta_yyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 324);

    auto ta_yyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 325);

    auto ta_yyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 326);

    auto ta_yyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 327);

    auto ta_yyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 328);

    auto ta_yyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 329);

    auto ta_yyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 330);

    auto ta_yyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 331);

    auto ta_yyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 332);

    auto ta_yyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 333);

    auto ta_yyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 334);

    auto ta_yyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 335);

    auto ta_yyyyz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 337);

    auto ta_yyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 338);

    auto ta_yyyyz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 339);

    auto ta_yyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 340);

    auto ta_yyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 341);

    auto ta_yyyyz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 342);

    auto ta_yyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 343);

    auto ta_yyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 344);

    auto ta_yyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 345);

    auto ta_yyyyz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 346);

    auto ta_yyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 347);

    auto ta_yyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 348);

    auto ta_yyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 349);

    auto ta_yyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 350);

    auto ta_yyyyz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 351);

    auto ta_yyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 352);

    auto ta_yyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 353);

    auto ta_yyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 354);

    auto ta_yyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 355);

    auto ta_yyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 356);

    auto ta_yyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 357);

    auto ta_yyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 358);

    auto ta_yyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 359);

    auto ta_yyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 360);

    auto ta_yyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 361);

    auto ta_yyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 362);

    auto ta_yyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 363);

    auto ta_yyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 364);

    auto ta_yyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 365);

    auto ta_yyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 366);

    auto ta_yyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 367);

    auto ta_yyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 368);

    auto ta_yyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 369);

    auto ta_yyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 370);

    auto ta_yyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 371);

    auto ta_yyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 372);

    auto ta_yyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 373);

    auto ta_yyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 374);

    auto ta_yyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 375);

    auto ta_yyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 376);

    auto ta_yyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 377);

    auto ta_yyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 378);

    auto ta_yyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 379);

    auto ta_yyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 380);

    auto ta_yyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 381);

    auto ta_yyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 382);

    auto ta_yyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 383);

    auto ta_yyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 384);

    auto ta_yyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 385);

    auto ta_yyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 386);

    auto ta_yyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 387);

    auto ta_yyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 388);

    auto ta_yyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 389);

    auto ta_yyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 390);

    auto ta_yyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 391);

    auto ta_yyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 392);

    auto ta_yyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 393);

    auto ta_yyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 394);

    auto ta_yyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 395);

    auto ta_yyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 396);

    auto ta_yyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 397);

    auto ta_yyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 398);

    auto ta_yzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 399);

    auto ta_yzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 400);

    auto ta_yzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 401);

    auto ta_yzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 402);

    auto ta_yzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 403);

    auto ta_yzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 404);

    auto ta_yzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 405);

    auto ta_yzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 406);

    auto ta_yzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 407);

    auto ta_yzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 408);

    auto ta_yzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 409);

    auto ta_yzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 410);

    auto ta_yzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 411);

    auto ta_yzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 412);

    auto ta_yzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 413);

    auto ta_yzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 414);

    auto ta_yzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 415);

    auto ta_yzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 416);

    auto ta_yzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 417);

    auto ta_yzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 418);

    auto ta_yzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 419);

    auto ta_zzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_hh + 420);

    auto ta_zzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_hh + 421);

    auto ta_zzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_hh + 422);

    auto ta_zzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_hh + 423);

    auto ta_zzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_hh + 424);

    auto ta_zzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_hh + 425);

    auto ta_zzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_hh + 426);

    auto ta_zzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_hh + 427);

    auto ta_zzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_hh + 428);

    auto ta_zzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_hh + 429);

    auto ta_zzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_hh + 430);

    auto ta_zzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_hh + 431);

    auto ta_zzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_hh + 432);

    auto ta_zzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_hh + 433);

    auto ta_zzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_hh + 434);

    auto ta_zzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_hh + 435);

    auto ta_zzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_hh + 436);

    auto ta_zzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_hh + 437);

    auto ta_zzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_hh + 438);

    auto ta_zzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_hh + 439);

    auto ta_zzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_hh + 440);

    // Set up components of auxiliary buffer : HH

    auto ta_xxxxx_xxxxx_1 = pbuffer.data(idx_npot_1_hh);

    auto ta_xxxxx_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 1);

    auto ta_xxxxx_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 2);

    auto ta_xxxxx_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 3);

    auto ta_xxxxx_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 4);

    auto ta_xxxxx_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 5);

    auto ta_xxxxx_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 6);

    auto ta_xxxxx_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 7);

    auto ta_xxxxx_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 8);

    auto ta_xxxxx_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 9);

    auto ta_xxxxx_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 10);

    auto ta_xxxxx_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 11);

    auto ta_xxxxx_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 12);

    auto ta_xxxxx_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 13);

    auto ta_xxxxx_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 14);

    auto ta_xxxxx_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 15);

    auto ta_xxxxx_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 16);

    auto ta_xxxxx_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 17);

    auto ta_xxxxx_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 18);

    auto ta_xxxxx_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 19);

    auto ta_xxxxx_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 20);

    auto ta_xxxxy_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 21);

    auto ta_xxxxy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 22);

    auto ta_xxxxy_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 23);

    auto ta_xxxxy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 24);

    auto ta_xxxxy_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 26);

    auto ta_xxxxy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 27);

    auto ta_xxxxy_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 30);

    auto ta_xxxxy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 31);

    auto ta_xxxxy_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 35);

    auto ta_xxxxy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 36);

    auto ta_xxxxy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 37);

    auto ta_xxxxy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 38);

    auto ta_xxxxy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 39);

    auto ta_xxxxy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 40);

    auto ta_xxxxz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 42);

    auto ta_xxxxz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 43);

    auto ta_xxxxz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 44);

    auto ta_xxxxz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 45);

    auto ta_xxxxz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 46);

    auto ta_xxxxz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 47);

    auto ta_xxxxz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 48);

    auto ta_xxxxz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 49);

    auto ta_xxxxz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 50);

    auto ta_xxxxz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 51);

    auto ta_xxxxz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 52);

    auto ta_xxxxz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 53);

    auto ta_xxxxz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 54);

    auto ta_xxxxz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 55);

    auto ta_xxxxz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 56);

    auto ta_xxxxz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 58);

    auto ta_xxxxz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 59);

    auto ta_xxxxz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 60);

    auto ta_xxxxz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 61);

    auto ta_xxxxz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 62);

    auto ta_xxxyy_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 63);

    auto ta_xxxyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 64);

    auto ta_xxxyy_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 65);

    auto ta_xxxyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 66);

    auto ta_xxxyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 67);

    auto ta_xxxyy_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 68);

    auto ta_xxxyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 69);

    auto ta_xxxyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 70);

    auto ta_xxxyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 71);

    auto ta_xxxyy_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 72);

    auto ta_xxxyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 73);

    auto ta_xxxyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 74);

    auto ta_xxxyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 75);

    auto ta_xxxyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 76);

    auto ta_xxxyy_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 77);

    auto ta_xxxyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 78);

    auto ta_xxxyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 79);

    auto ta_xxxyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 80);

    auto ta_xxxyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 81);

    auto ta_xxxyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 82);

    auto ta_xxxyy_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 83);

    auto ta_xxxyz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 86);

    auto ta_xxxyz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 89);

    auto ta_xxxyz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 93);

    auto ta_xxxyz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 98);

    auto ta_xxxyz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 100);

    auto ta_xxxyz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 101);

    auto ta_xxxyz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 102);

    auto ta_xxxyz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 103);

    auto ta_xxxzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 105);

    auto ta_xxxzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 106);

    auto ta_xxxzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 107);

    auto ta_xxxzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 108);

    auto ta_xxxzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 109);

    auto ta_xxxzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 110);

    auto ta_xxxzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 111);

    auto ta_xxxzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 112);

    auto ta_xxxzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 113);

    auto ta_xxxzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 114);

    auto ta_xxxzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 115);

    auto ta_xxxzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 116);

    auto ta_xxxzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 117);

    auto ta_xxxzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 118);

    auto ta_xxxzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 119);

    auto ta_xxxzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 120);

    auto ta_xxxzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 121);

    auto ta_xxxzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 122);

    auto ta_xxxzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 123);

    auto ta_xxxzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 124);

    auto ta_xxxzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 125);

    auto ta_xxyyy_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 126);

    auto ta_xxyyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 127);

    auto ta_xxyyy_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 128);

    auto ta_xxyyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 129);

    auto ta_xxyyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 130);

    auto ta_xxyyy_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 131);

    auto ta_xxyyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 132);

    auto ta_xxyyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 133);

    auto ta_xxyyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 134);

    auto ta_xxyyy_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 135);

    auto ta_xxyyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 136);

    auto ta_xxyyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 137);

    auto ta_xxyyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 138);

    auto ta_xxyyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 139);

    auto ta_xxyyy_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 140);

    auto ta_xxyyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 141);

    auto ta_xxyyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 142);

    auto ta_xxyyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 143);

    auto ta_xxyyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 144);

    auto ta_xxyyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 145);

    auto ta_xxyyy_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 146);

    auto ta_xxyyz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 148);

    auto ta_xxyyz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 149);

    auto ta_xxyyz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 150);

    auto ta_xxyyz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 152);

    auto ta_xxyyz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 153);

    auto ta_xxyyz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 156);

    auto ta_xxyyz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 157);

    auto ta_xxyyz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 161);

    auto ta_xxyyz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 163);

    auto ta_xxyyz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 164);

    auto ta_xxyyz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 165);

    auto ta_xxyyz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 166);

    auto ta_xxyyz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 167);

    auto ta_xxyzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 168);

    auto ta_xxyzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 170);

    auto ta_xxyzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 173);

    auto ta_xxyzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 177);

    auto ta_xxyzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 182);

    auto ta_xxyzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 183);

    auto ta_xxyzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 184);

    auto ta_xxyzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 185);

    auto ta_xxyzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 186);

    auto ta_xxyzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 187);

    auto ta_xxzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 189);

    auto ta_xxzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 190);

    auto ta_xxzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 191);

    auto ta_xxzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 192);

    auto ta_xxzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 193);

    auto ta_xxzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 194);

    auto ta_xxzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 195);

    auto ta_xxzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 196);

    auto ta_xxzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 197);

    auto ta_xxzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 198);

    auto ta_xxzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 199);

    auto ta_xxzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 200);

    auto ta_xxzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 201);

    auto ta_xxzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 202);

    auto ta_xxzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 203);

    auto ta_xxzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 204);

    auto ta_xxzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 205);

    auto ta_xxzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 206);

    auto ta_xxzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 207);

    auto ta_xxzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 208);

    auto ta_xxzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 209);

    auto ta_xyyyy_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 210);

    auto ta_xyyyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 211);

    auto ta_xyyyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 213);

    auto ta_xyyyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 214);

    auto ta_xyyyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 216);

    auto ta_xyyyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 217);

    auto ta_xyyyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 218);

    auto ta_xyyyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 220);

    auto ta_xyyyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 221);

    auto ta_xyyyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 222);

    auto ta_xyyyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 223);

    auto ta_xyyyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 225);

    auto ta_xyyyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 226);

    auto ta_xyyyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 227);

    auto ta_xyyyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 228);

    auto ta_xyyyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 229);

    auto ta_xyyyy_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 230);

    auto ta_xyyyz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 247);

    auto ta_xyyyz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 248);

    auto ta_xyyyz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 249);

    auto ta_xyyyz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 250);

    auto ta_xyyyz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 251);

    auto ta_xyyzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 256);

    auto ta_xyyzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 259);

    auto ta_xyyzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 260);

    auto ta_xyyzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 263);

    auto ta_xyyzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 264);

    auto ta_xyyzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 265);

    auto ta_xyyzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 267);

    auto ta_xyyzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 268);

    auto ta_xyyzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 269);

    auto ta_xyyzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 270);

    auto ta_xyyzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 271);

    auto ta_xyyzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 272);

    auto ta_xyzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 288);

    auto ta_xyzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 289);

    auto ta_xyzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 290);

    auto ta_xyzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 291);

    auto ta_xyzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 292);

    auto ta_xzzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 294);

    auto ta_xzzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 296);

    auto ta_xzzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 298);

    auto ta_xzzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 299);

    auto ta_xzzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 301);

    auto ta_xzzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 302);

    auto ta_xzzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 303);

    auto ta_xzzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 305);

    auto ta_xzzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 306);

    auto ta_xzzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 307);

    auto ta_xzzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 308);

    auto ta_xzzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 309);

    auto ta_xzzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 310);

    auto ta_xzzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 311);

    auto ta_xzzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 312);

    auto ta_xzzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 313);

    auto ta_xzzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 314);

    auto ta_yyyyy_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 315);

    auto ta_yyyyy_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 316);

    auto ta_yyyyy_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 317);

    auto ta_yyyyy_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 318);

    auto ta_yyyyy_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 319);

    auto ta_yyyyy_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 320);

    auto ta_yyyyy_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 321);

    auto ta_yyyyy_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 322);

    auto ta_yyyyy_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 323);

    auto ta_yyyyy_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 324);

    auto ta_yyyyy_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 325);

    auto ta_yyyyy_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 326);

    auto ta_yyyyy_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 327);

    auto ta_yyyyy_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 328);

    auto ta_yyyyy_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 329);

    auto ta_yyyyy_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 330);

    auto ta_yyyyy_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 331);

    auto ta_yyyyy_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 332);

    auto ta_yyyyy_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 333);

    auto ta_yyyyy_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 334);

    auto ta_yyyyy_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 335);

    auto ta_yyyyz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 337);

    auto ta_yyyyz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 338);

    auto ta_yyyyz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 339);

    auto ta_yyyyz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 340);

    auto ta_yyyyz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 341);

    auto ta_yyyyz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 342);

    auto ta_yyyyz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 343);

    auto ta_yyyyz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 344);

    auto ta_yyyyz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 345);

    auto ta_yyyyz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 346);

    auto ta_yyyyz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 347);

    auto ta_yyyyz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 348);

    auto ta_yyyyz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 349);

    auto ta_yyyyz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 350);

    auto ta_yyyyz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 351);

    auto ta_yyyyz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 352);

    auto ta_yyyyz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 353);

    auto ta_yyyyz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 354);

    auto ta_yyyyz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 355);

    auto ta_yyyyz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 356);

    auto ta_yyyzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 357);

    auto ta_yyyzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 358);

    auto ta_yyyzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 359);

    auto ta_yyyzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 360);

    auto ta_yyyzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 361);

    auto ta_yyyzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 362);

    auto ta_yyyzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 363);

    auto ta_yyyzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 364);

    auto ta_yyyzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 365);

    auto ta_yyyzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 366);

    auto ta_yyyzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 367);

    auto ta_yyyzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 368);

    auto ta_yyyzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 369);

    auto ta_yyyzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 370);

    auto ta_yyyzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 371);

    auto ta_yyyzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 372);

    auto ta_yyyzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 373);

    auto ta_yyyzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 374);

    auto ta_yyyzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 375);

    auto ta_yyyzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 376);

    auto ta_yyyzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 377);

    auto ta_yyzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 378);

    auto ta_yyzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 379);

    auto ta_yyzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 380);

    auto ta_yyzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 381);

    auto ta_yyzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 382);

    auto ta_yyzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 383);

    auto ta_yyzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 384);

    auto ta_yyzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 385);

    auto ta_yyzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 386);

    auto ta_yyzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 387);

    auto ta_yyzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 388);

    auto ta_yyzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 389);

    auto ta_yyzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 390);

    auto ta_yyzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 391);

    auto ta_yyzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 392);

    auto ta_yyzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 393);

    auto ta_yyzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 394);

    auto ta_yyzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 395);

    auto ta_yyzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 396);

    auto ta_yyzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 397);

    auto ta_yyzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 398);

    auto ta_yzzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 399);

    auto ta_yzzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 400);

    auto ta_yzzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 401);

    auto ta_yzzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 402);

    auto ta_yzzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 403);

    auto ta_yzzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 404);

    auto ta_yzzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 405);

    auto ta_yzzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 406);

    auto ta_yzzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 407);

    auto ta_yzzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 408);

    auto ta_yzzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 409);

    auto ta_yzzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 410);

    auto ta_yzzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 411);

    auto ta_yzzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 412);

    auto ta_yzzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 413);

    auto ta_yzzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 414);

    auto ta_yzzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 415);

    auto ta_yzzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 416);

    auto ta_yzzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 417);

    auto ta_yzzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 418);

    auto ta_yzzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 419);

    auto ta_zzzzz_xxxxx_1 = pbuffer.data(idx_npot_1_hh + 420);

    auto ta_zzzzz_xxxxy_1 = pbuffer.data(idx_npot_1_hh + 421);

    auto ta_zzzzz_xxxxz_1 = pbuffer.data(idx_npot_1_hh + 422);

    auto ta_zzzzz_xxxyy_1 = pbuffer.data(idx_npot_1_hh + 423);

    auto ta_zzzzz_xxxyz_1 = pbuffer.data(idx_npot_1_hh + 424);

    auto ta_zzzzz_xxxzz_1 = pbuffer.data(idx_npot_1_hh + 425);

    auto ta_zzzzz_xxyyy_1 = pbuffer.data(idx_npot_1_hh + 426);

    auto ta_zzzzz_xxyyz_1 = pbuffer.data(idx_npot_1_hh + 427);

    auto ta_zzzzz_xxyzz_1 = pbuffer.data(idx_npot_1_hh + 428);

    auto ta_zzzzz_xxzzz_1 = pbuffer.data(idx_npot_1_hh + 429);

    auto ta_zzzzz_xyyyy_1 = pbuffer.data(idx_npot_1_hh + 430);

    auto ta_zzzzz_xyyyz_1 = pbuffer.data(idx_npot_1_hh + 431);

    auto ta_zzzzz_xyyzz_1 = pbuffer.data(idx_npot_1_hh + 432);

    auto ta_zzzzz_xyzzz_1 = pbuffer.data(idx_npot_1_hh + 433);

    auto ta_zzzzz_xzzzz_1 = pbuffer.data(idx_npot_1_hh + 434);

    auto ta_zzzzz_yyyyy_1 = pbuffer.data(idx_npot_1_hh + 435);

    auto ta_zzzzz_yyyyz_1 = pbuffer.data(idx_npot_1_hh + 436);

    auto ta_zzzzz_yyyzz_1 = pbuffer.data(idx_npot_1_hh + 437);

    auto ta_zzzzz_yyzzz_1 = pbuffer.data(idx_npot_1_hh + 438);

    auto ta_zzzzz_yzzzz_1 = pbuffer.data(idx_npot_1_hh + 439);

    auto ta_zzzzz_zzzzz_1 = pbuffer.data(idx_npot_1_hh + 440);

    // Set up 0-21 components of targeted buffer : IH

    auto ta_xxxxxx_xxxxx_0 = pbuffer.data(idx_npot_0_ih);

    auto ta_xxxxxx_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 1);

    auto ta_xxxxxx_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 2);

    auto ta_xxxxxx_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 3);

    auto ta_xxxxxx_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 4);

    auto ta_xxxxxx_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 5);

    auto ta_xxxxxx_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 6);

    auto ta_xxxxxx_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 7);

    auto ta_xxxxxx_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 8);

    auto ta_xxxxxx_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 9);

    auto ta_xxxxxx_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 10);

    auto ta_xxxxxx_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 11);

    auto ta_xxxxxx_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 12);

    auto ta_xxxxxx_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 13);

    auto ta_xxxxxx_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 14);

    auto ta_xxxxxx_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 15);

    auto ta_xxxxxx_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 16);

    auto ta_xxxxxx_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 17);

    auto ta_xxxxxx_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 18);

    auto ta_xxxxxx_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 19);

    auto ta_xxxxxx_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 20);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xxxx_xxxxx_0,   \
                             ta_xxxx_xxxxx_1,   \
                             ta_xxxx_xxxxy_0,   \
                             ta_xxxx_xxxxy_1,   \
                             ta_xxxx_xxxxz_0,   \
                             ta_xxxx_xxxxz_1,   \
                             ta_xxxx_xxxyy_0,   \
                             ta_xxxx_xxxyy_1,   \
                             ta_xxxx_xxxyz_0,   \
                             ta_xxxx_xxxyz_1,   \
                             ta_xxxx_xxxzz_0,   \
                             ta_xxxx_xxxzz_1,   \
                             ta_xxxx_xxyyy_0,   \
                             ta_xxxx_xxyyy_1,   \
                             ta_xxxx_xxyyz_0,   \
                             ta_xxxx_xxyyz_1,   \
                             ta_xxxx_xxyzz_0,   \
                             ta_xxxx_xxyzz_1,   \
                             ta_xxxx_xxzzz_0,   \
                             ta_xxxx_xxzzz_1,   \
                             ta_xxxx_xyyyy_0,   \
                             ta_xxxx_xyyyy_1,   \
                             ta_xxxx_xyyyz_0,   \
                             ta_xxxx_xyyyz_1,   \
                             ta_xxxx_xyyzz_0,   \
                             ta_xxxx_xyyzz_1,   \
                             ta_xxxx_xyzzz_0,   \
                             ta_xxxx_xyzzz_1,   \
                             ta_xxxx_xzzzz_0,   \
                             ta_xxxx_xzzzz_1,   \
                             ta_xxxx_yyyyy_0,   \
                             ta_xxxx_yyyyy_1,   \
                             ta_xxxx_yyyyz_0,   \
                             ta_xxxx_yyyyz_1,   \
                             ta_xxxx_yyyzz_0,   \
                             ta_xxxx_yyyzz_1,   \
                             ta_xxxx_yyzzz_0,   \
                             ta_xxxx_yyzzz_1,   \
                             ta_xxxx_yzzzz_0,   \
                             ta_xxxx_yzzzz_1,   \
                             ta_xxxx_zzzzz_0,   \
                             ta_xxxx_zzzzz_1,   \
                             ta_xxxxx_xxxx_0,   \
                             ta_xxxxx_xxxx_1,   \
                             ta_xxxxx_xxxxx_0,  \
                             ta_xxxxx_xxxxx_1,  \
                             ta_xxxxx_xxxxy_0,  \
                             ta_xxxxx_xxxxy_1,  \
                             ta_xxxxx_xxxxz_0,  \
                             ta_xxxxx_xxxxz_1,  \
                             ta_xxxxx_xxxy_0,   \
                             ta_xxxxx_xxxy_1,   \
                             ta_xxxxx_xxxyy_0,  \
                             ta_xxxxx_xxxyy_1,  \
                             ta_xxxxx_xxxyz_0,  \
                             ta_xxxxx_xxxyz_1,  \
                             ta_xxxxx_xxxz_0,   \
                             ta_xxxxx_xxxz_1,   \
                             ta_xxxxx_xxxzz_0,  \
                             ta_xxxxx_xxxzz_1,  \
                             ta_xxxxx_xxyy_0,   \
                             ta_xxxxx_xxyy_1,   \
                             ta_xxxxx_xxyyy_0,  \
                             ta_xxxxx_xxyyy_1,  \
                             ta_xxxxx_xxyyz_0,  \
                             ta_xxxxx_xxyyz_1,  \
                             ta_xxxxx_xxyz_0,   \
                             ta_xxxxx_xxyz_1,   \
                             ta_xxxxx_xxyzz_0,  \
                             ta_xxxxx_xxyzz_1,  \
                             ta_xxxxx_xxzz_0,   \
                             ta_xxxxx_xxzz_1,   \
                             ta_xxxxx_xxzzz_0,  \
                             ta_xxxxx_xxzzz_1,  \
                             ta_xxxxx_xyyy_0,   \
                             ta_xxxxx_xyyy_1,   \
                             ta_xxxxx_xyyyy_0,  \
                             ta_xxxxx_xyyyy_1,  \
                             ta_xxxxx_xyyyz_0,  \
                             ta_xxxxx_xyyyz_1,  \
                             ta_xxxxx_xyyz_0,   \
                             ta_xxxxx_xyyz_1,   \
                             ta_xxxxx_xyyzz_0,  \
                             ta_xxxxx_xyyzz_1,  \
                             ta_xxxxx_xyzz_0,   \
                             ta_xxxxx_xyzz_1,   \
                             ta_xxxxx_xyzzz_0,  \
                             ta_xxxxx_xyzzz_1,  \
                             ta_xxxxx_xzzz_0,   \
                             ta_xxxxx_xzzz_1,   \
                             ta_xxxxx_xzzzz_0,  \
                             ta_xxxxx_xzzzz_1,  \
                             ta_xxxxx_yyyy_0,   \
                             ta_xxxxx_yyyy_1,   \
                             ta_xxxxx_yyyyy_0,  \
                             ta_xxxxx_yyyyy_1,  \
                             ta_xxxxx_yyyyz_0,  \
                             ta_xxxxx_yyyyz_1,  \
                             ta_xxxxx_yyyz_0,   \
                             ta_xxxxx_yyyz_1,   \
                             ta_xxxxx_yyyzz_0,  \
                             ta_xxxxx_yyyzz_1,  \
                             ta_xxxxx_yyzz_0,   \
                             ta_xxxxx_yyzz_1,   \
                             ta_xxxxx_yyzzz_0,  \
                             ta_xxxxx_yyzzz_1,  \
                             ta_xxxxx_yzzz_0,   \
                             ta_xxxxx_yzzz_1,   \
                             ta_xxxxx_yzzzz_0,  \
                             ta_xxxxx_yzzzz_1,  \
                             ta_xxxxx_zzzz_0,   \
                             ta_xxxxx_zzzz_1,   \
                             ta_xxxxx_zzzzz_0,  \
                             ta_xxxxx_zzzzz_1,  \
                             ta_xxxxxx_xxxxx_0, \
                             ta_xxxxxx_xxxxy_0, \
                             ta_xxxxxx_xxxxz_0, \
                             ta_xxxxxx_xxxyy_0, \
                             ta_xxxxxx_xxxyz_0, \
                             ta_xxxxxx_xxxzz_0, \
                             ta_xxxxxx_xxyyy_0, \
                             ta_xxxxxx_xxyyz_0, \
                             ta_xxxxxx_xxyzz_0, \
                             ta_xxxxxx_xxzzz_0, \
                             ta_xxxxxx_xyyyy_0, \
                             ta_xxxxxx_xyyyz_0, \
                             ta_xxxxxx_xyyzz_0, \
                             ta_xxxxxx_xyzzz_0, \
                             ta_xxxxxx_xzzzz_0, \
                             ta_xxxxxx_yyyyy_0, \
                             ta_xxxxxx_yyyyz_0, \
                             ta_xxxxxx_yyyzz_0, \
                             ta_xxxxxx_yyzzz_0, \
                             ta_xxxxxx_yzzzz_0, \
                             ta_xxxxxx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_xxxxx_0[i] = 5.0 * ta_xxxx_xxxxx_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxx_1[i] * fe_0 + 5.0 * ta_xxxxx_xxxx_0[i] * fe_0 -
                               5.0 * ta_xxxxx_xxxx_1[i] * fe_0 + ta_xxxxx_xxxxx_0[i] * pa_x[i] - ta_xxxxx_xxxxx_1[i] * pc_x[i];

        ta_xxxxxx_xxxxy_0[i] = 5.0 * ta_xxxx_xxxxy_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxy_1[i] * fe_0 + 4.0 * ta_xxxxx_xxxy_0[i] * fe_0 -
                               4.0 * ta_xxxxx_xxxy_1[i] * fe_0 + ta_xxxxx_xxxxy_0[i] * pa_x[i] - ta_xxxxx_xxxxy_1[i] * pc_x[i];

        ta_xxxxxx_xxxxz_0[i] = 5.0 * ta_xxxx_xxxxz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxxz_1[i] * fe_0 + 4.0 * ta_xxxxx_xxxz_0[i] * fe_0 -
                               4.0 * ta_xxxxx_xxxz_1[i] * fe_0 + ta_xxxxx_xxxxz_0[i] * pa_x[i] - ta_xxxxx_xxxxz_1[i] * pc_x[i];

        ta_xxxxxx_xxxyy_0[i] = 5.0 * ta_xxxx_xxxyy_0[i] * fe_0 - 5.0 * ta_xxxx_xxxyy_1[i] * fe_0 + 3.0 * ta_xxxxx_xxyy_0[i] * fe_0 -
                               3.0 * ta_xxxxx_xxyy_1[i] * fe_0 + ta_xxxxx_xxxyy_0[i] * pa_x[i] - ta_xxxxx_xxxyy_1[i] * pc_x[i];

        ta_xxxxxx_xxxyz_0[i] = 5.0 * ta_xxxx_xxxyz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxyz_1[i] * fe_0 + 3.0 * ta_xxxxx_xxyz_0[i] * fe_0 -
                               3.0 * ta_xxxxx_xxyz_1[i] * fe_0 + ta_xxxxx_xxxyz_0[i] * pa_x[i] - ta_xxxxx_xxxyz_1[i] * pc_x[i];

        ta_xxxxxx_xxxzz_0[i] = 5.0 * ta_xxxx_xxxzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxxzz_1[i] * fe_0 + 3.0 * ta_xxxxx_xxzz_0[i] * fe_0 -
                               3.0 * ta_xxxxx_xxzz_1[i] * fe_0 + ta_xxxxx_xxxzz_0[i] * pa_x[i] - ta_xxxxx_xxxzz_1[i] * pc_x[i];

        ta_xxxxxx_xxyyy_0[i] = 5.0 * ta_xxxx_xxyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xxyyy_1[i] * fe_0 + 2.0 * ta_xxxxx_xyyy_0[i] * fe_0 -
                               2.0 * ta_xxxxx_xyyy_1[i] * fe_0 + ta_xxxxx_xxyyy_0[i] * pa_x[i] - ta_xxxxx_xxyyy_1[i] * pc_x[i];

        ta_xxxxxx_xxyyz_0[i] = 5.0 * ta_xxxx_xxyyz_0[i] * fe_0 - 5.0 * ta_xxxx_xxyyz_1[i] * fe_0 + 2.0 * ta_xxxxx_xyyz_0[i] * fe_0 -
                               2.0 * ta_xxxxx_xyyz_1[i] * fe_0 + ta_xxxxx_xxyyz_0[i] * pa_x[i] - ta_xxxxx_xxyyz_1[i] * pc_x[i];

        ta_xxxxxx_xxyzz_0[i] = 5.0 * ta_xxxx_xxyzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxyzz_1[i] * fe_0 + 2.0 * ta_xxxxx_xyzz_0[i] * fe_0 -
                               2.0 * ta_xxxxx_xyzz_1[i] * fe_0 + ta_xxxxx_xxyzz_0[i] * pa_x[i] - ta_xxxxx_xxyzz_1[i] * pc_x[i];

        ta_xxxxxx_xxzzz_0[i] = 5.0 * ta_xxxx_xxzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xxzzz_1[i] * fe_0 + 2.0 * ta_xxxxx_xzzz_0[i] * fe_0 -
                               2.0 * ta_xxxxx_xzzz_1[i] * fe_0 + ta_xxxxx_xxzzz_0[i] * pa_x[i] - ta_xxxxx_xxzzz_1[i] * pc_x[i];

        ta_xxxxxx_xyyyy_0[i] = 5.0 * ta_xxxx_xyyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xyyyy_1[i] * fe_0 + ta_xxxxx_yyyy_0[i] * fe_0 -
                               ta_xxxxx_yyyy_1[i] * fe_0 + ta_xxxxx_xyyyy_0[i] * pa_x[i] - ta_xxxxx_xyyyy_1[i] * pc_x[i];

        ta_xxxxxx_xyyyz_0[i] = 5.0 * ta_xxxx_xyyyz_0[i] * fe_0 - 5.0 * ta_xxxx_xyyyz_1[i] * fe_0 + ta_xxxxx_yyyz_0[i] * fe_0 -
                               ta_xxxxx_yyyz_1[i] * fe_0 + ta_xxxxx_xyyyz_0[i] * pa_x[i] - ta_xxxxx_xyyyz_1[i] * pc_x[i];

        ta_xxxxxx_xyyzz_0[i] = 5.0 * ta_xxxx_xyyzz_0[i] * fe_0 - 5.0 * ta_xxxx_xyyzz_1[i] * fe_0 + ta_xxxxx_yyzz_0[i] * fe_0 -
                               ta_xxxxx_yyzz_1[i] * fe_0 + ta_xxxxx_xyyzz_0[i] * pa_x[i] - ta_xxxxx_xyyzz_1[i] * pc_x[i];

        ta_xxxxxx_xyzzz_0[i] = 5.0 * ta_xxxx_xyzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xyzzz_1[i] * fe_0 + ta_xxxxx_yzzz_0[i] * fe_0 -
                               ta_xxxxx_yzzz_1[i] * fe_0 + ta_xxxxx_xyzzz_0[i] * pa_x[i] - ta_xxxxx_xyzzz_1[i] * pc_x[i];

        ta_xxxxxx_xzzzz_0[i] = 5.0 * ta_xxxx_xzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xzzzz_1[i] * fe_0 + ta_xxxxx_zzzz_0[i] * fe_0 -
                               ta_xxxxx_zzzz_1[i] * fe_0 + ta_xxxxx_xzzzz_0[i] * pa_x[i] - ta_xxxxx_xzzzz_1[i] * pc_x[i];

        ta_xxxxxx_yyyyy_0[i] =
            5.0 * ta_xxxx_yyyyy_0[i] * fe_0 - 5.0 * ta_xxxx_yyyyy_1[i] * fe_0 + ta_xxxxx_yyyyy_0[i] * pa_x[i] - ta_xxxxx_yyyyy_1[i] * pc_x[i];

        ta_xxxxxx_yyyyz_0[i] =
            5.0 * ta_xxxx_yyyyz_0[i] * fe_0 - 5.0 * ta_xxxx_yyyyz_1[i] * fe_0 + ta_xxxxx_yyyyz_0[i] * pa_x[i] - ta_xxxxx_yyyyz_1[i] * pc_x[i];

        ta_xxxxxx_yyyzz_0[i] =
            5.0 * ta_xxxx_yyyzz_0[i] * fe_0 - 5.0 * ta_xxxx_yyyzz_1[i] * fe_0 + ta_xxxxx_yyyzz_0[i] * pa_x[i] - ta_xxxxx_yyyzz_1[i] * pc_x[i];

        ta_xxxxxx_yyzzz_0[i] =
            5.0 * ta_xxxx_yyzzz_0[i] * fe_0 - 5.0 * ta_xxxx_yyzzz_1[i] * fe_0 + ta_xxxxx_yyzzz_0[i] * pa_x[i] - ta_xxxxx_yyzzz_1[i] * pc_x[i];

        ta_xxxxxx_yzzzz_0[i] =
            5.0 * ta_xxxx_yzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_yzzzz_1[i] * fe_0 + ta_xxxxx_yzzzz_0[i] * pa_x[i] - ta_xxxxx_yzzzz_1[i] * pc_x[i];

        ta_xxxxxx_zzzzz_0[i] =
            5.0 * ta_xxxx_zzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_zzzzz_1[i] * fe_0 + ta_xxxxx_zzzzz_0[i] * pa_x[i] - ta_xxxxx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : IH

    auto ta_xxxxxy_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 21);

    auto ta_xxxxxy_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 22);

    auto ta_xxxxxy_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 23);

    auto ta_xxxxxy_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 24);

    auto ta_xxxxxy_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 25);

    auto ta_xxxxxy_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 26);

    auto ta_xxxxxy_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 27);

    auto ta_xxxxxy_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 28);

    auto ta_xxxxxy_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 29);

    auto ta_xxxxxy_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 30);

    auto ta_xxxxxy_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 31);

    auto ta_xxxxxy_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 32);

    auto ta_xxxxxy_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 33);

    auto ta_xxxxxy_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 34);

    auto ta_xxxxxy_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 35);

    auto ta_xxxxxy_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 36);

    auto ta_xxxxxy_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 37);

    auto ta_xxxxxy_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 38);

    auto ta_xxxxxy_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 39);

    auto ta_xxxxxy_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 40);

    auto ta_xxxxxy_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 41);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxxxx_xxxx_0,   \
                             ta_xxxxx_xxxx_1,   \
                             ta_xxxxx_xxxxx_0,  \
                             ta_xxxxx_xxxxx_1,  \
                             ta_xxxxx_xxxxy_0,  \
                             ta_xxxxx_xxxxy_1,  \
                             ta_xxxxx_xxxxz_0,  \
                             ta_xxxxx_xxxxz_1,  \
                             ta_xxxxx_xxxy_0,   \
                             ta_xxxxx_xxxy_1,   \
                             ta_xxxxx_xxxyy_0,  \
                             ta_xxxxx_xxxyy_1,  \
                             ta_xxxxx_xxxyz_0,  \
                             ta_xxxxx_xxxyz_1,  \
                             ta_xxxxx_xxxz_0,   \
                             ta_xxxxx_xxxz_1,   \
                             ta_xxxxx_xxxzz_0,  \
                             ta_xxxxx_xxxzz_1,  \
                             ta_xxxxx_xxyy_0,   \
                             ta_xxxxx_xxyy_1,   \
                             ta_xxxxx_xxyyy_0,  \
                             ta_xxxxx_xxyyy_1,  \
                             ta_xxxxx_xxyyz_0,  \
                             ta_xxxxx_xxyyz_1,  \
                             ta_xxxxx_xxyz_0,   \
                             ta_xxxxx_xxyz_1,   \
                             ta_xxxxx_xxyzz_0,  \
                             ta_xxxxx_xxyzz_1,  \
                             ta_xxxxx_xxzz_0,   \
                             ta_xxxxx_xxzz_1,   \
                             ta_xxxxx_xxzzz_0,  \
                             ta_xxxxx_xxzzz_1,  \
                             ta_xxxxx_xyyy_0,   \
                             ta_xxxxx_xyyy_1,   \
                             ta_xxxxx_xyyyy_0,  \
                             ta_xxxxx_xyyyy_1,  \
                             ta_xxxxx_xyyyz_0,  \
                             ta_xxxxx_xyyyz_1,  \
                             ta_xxxxx_xyyz_0,   \
                             ta_xxxxx_xyyz_1,   \
                             ta_xxxxx_xyyzz_0,  \
                             ta_xxxxx_xyyzz_1,  \
                             ta_xxxxx_xyzz_0,   \
                             ta_xxxxx_xyzz_1,   \
                             ta_xxxxx_xyzzz_0,  \
                             ta_xxxxx_xyzzz_1,  \
                             ta_xxxxx_xzzz_0,   \
                             ta_xxxxx_xzzz_1,   \
                             ta_xxxxx_xzzzz_0,  \
                             ta_xxxxx_xzzzz_1,  \
                             ta_xxxxx_zzzzz_0,  \
                             ta_xxxxx_zzzzz_1,  \
                             ta_xxxxxy_xxxxx_0, \
                             ta_xxxxxy_xxxxy_0, \
                             ta_xxxxxy_xxxxz_0, \
                             ta_xxxxxy_xxxyy_0, \
                             ta_xxxxxy_xxxyz_0, \
                             ta_xxxxxy_xxxzz_0, \
                             ta_xxxxxy_xxyyy_0, \
                             ta_xxxxxy_xxyyz_0, \
                             ta_xxxxxy_xxyzz_0, \
                             ta_xxxxxy_xxzzz_0, \
                             ta_xxxxxy_xyyyy_0, \
                             ta_xxxxxy_xyyyz_0, \
                             ta_xxxxxy_xyyzz_0, \
                             ta_xxxxxy_xyzzz_0, \
                             ta_xxxxxy_xzzzz_0, \
                             ta_xxxxxy_yyyyy_0, \
                             ta_xxxxxy_yyyyz_0, \
                             ta_xxxxxy_yyyzz_0, \
                             ta_xxxxxy_yyzzz_0, \
                             ta_xxxxxy_yzzzz_0, \
                             ta_xxxxxy_zzzzz_0, \
                             ta_xxxxy_yyyyy_0,  \
                             ta_xxxxy_yyyyy_1,  \
                             ta_xxxxy_yyyyz_0,  \
                             ta_xxxxy_yyyyz_1,  \
                             ta_xxxxy_yyyzz_0,  \
                             ta_xxxxy_yyyzz_1,  \
                             ta_xxxxy_yyzzz_0,  \
                             ta_xxxxy_yyzzz_1,  \
                             ta_xxxxy_yzzzz_0,  \
                             ta_xxxxy_yzzzz_1,  \
                             ta_xxxy_yyyyy_0,   \
                             ta_xxxy_yyyyy_1,   \
                             ta_xxxy_yyyyz_0,   \
                             ta_xxxy_yyyyz_1,   \
                             ta_xxxy_yyyzz_0,   \
                             ta_xxxy_yyyzz_1,   \
                             ta_xxxy_yyzzz_0,   \
                             ta_xxxy_yyzzz_1,   \
                             ta_xxxy_yzzzz_0,   \
                             ta_xxxy_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxy_xxxxx_0[i] = ta_xxxxx_xxxxx_0[i] * pa_y[i] - ta_xxxxx_xxxxx_1[i] * pc_y[i];

        ta_xxxxxy_xxxxy_0[i] = ta_xxxxx_xxxx_0[i] * fe_0 - ta_xxxxx_xxxx_1[i] * fe_0 + ta_xxxxx_xxxxy_0[i] * pa_y[i] - ta_xxxxx_xxxxy_1[i] * pc_y[i];

        ta_xxxxxy_xxxxz_0[i] = ta_xxxxx_xxxxz_0[i] * pa_y[i] - ta_xxxxx_xxxxz_1[i] * pc_y[i];

        ta_xxxxxy_xxxyy_0[i] =
            2.0 * ta_xxxxx_xxxy_0[i] * fe_0 - 2.0 * ta_xxxxx_xxxy_1[i] * fe_0 + ta_xxxxx_xxxyy_0[i] * pa_y[i] - ta_xxxxx_xxxyy_1[i] * pc_y[i];

        ta_xxxxxy_xxxyz_0[i] = ta_xxxxx_xxxz_0[i] * fe_0 - ta_xxxxx_xxxz_1[i] * fe_0 + ta_xxxxx_xxxyz_0[i] * pa_y[i] - ta_xxxxx_xxxyz_1[i] * pc_y[i];

        ta_xxxxxy_xxxzz_0[i] = ta_xxxxx_xxxzz_0[i] * pa_y[i] - ta_xxxxx_xxxzz_1[i] * pc_y[i];

        ta_xxxxxy_xxyyy_0[i] =
            3.0 * ta_xxxxx_xxyy_0[i] * fe_0 - 3.0 * ta_xxxxx_xxyy_1[i] * fe_0 + ta_xxxxx_xxyyy_0[i] * pa_y[i] - ta_xxxxx_xxyyy_1[i] * pc_y[i];

        ta_xxxxxy_xxyyz_0[i] =
            2.0 * ta_xxxxx_xxyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxyz_1[i] * fe_0 + ta_xxxxx_xxyyz_0[i] * pa_y[i] - ta_xxxxx_xxyyz_1[i] * pc_y[i];

        ta_xxxxxy_xxyzz_0[i] = ta_xxxxx_xxzz_0[i] * fe_0 - ta_xxxxx_xxzz_1[i] * fe_0 + ta_xxxxx_xxyzz_0[i] * pa_y[i] - ta_xxxxx_xxyzz_1[i] * pc_y[i];

        ta_xxxxxy_xxzzz_0[i] = ta_xxxxx_xxzzz_0[i] * pa_y[i] - ta_xxxxx_xxzzz_1[i] * pc_y[i];

        ta_xxxxxy_xyyyy_0[i] =
            4.0 * ta_xxxxx_xyyy_0[i] * fe_0 - 4.0 * ta_xxxxx_xyyy_1[i] * fe_0 + ta_xxxxx_xyyyy_0[i] * pa_y[i] - ta_xxxxx_xyyyy_1[i] * pc_y[i];

        ta_xxxxxy_xyyyz_0[i] =
            3.0 * ta_xxxxx_xyyz_0[i] * fe_0 - 3.0 * ta_xxxxx_xyyz_1[i] * fe_0 + ta_xxxxx_xyyyz_0[i] * pa_y[i] - ta_xxxxx_xyyyz_1[i] * pc_y[i];

        ta_xxxxxy_xyyzz_0[i] =
            2.0 * ta_xxxxx_xyzz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyzz_1[i] * fe_0 + ta_xxxxx_xyyzz_0[i] * pa_y[i] - ta_xxxxx_xyyzz_1[i] * pc_y[i];

        ta_xxxxxy_xyzzz_0[i] = ta_xxxxx_xzzz_0[i] * fe_0 - ta_xxxxx_xzzz_1[i] * fe_0 + ta_xxxxx_xyzzz_0[i] * pa_y[i] - ta_xxxxx_xyzzz_1[i] * pc_y[i];

        ta_xxxxxy_xzzzz_0[i] = ta_xxxxx_xzzzz_0[i] * pa_y[i] - ta_xxxxx_xzzzz_1[i] * pc_y[i];

        ta_xxxxxy_yyyyy_0[i] =
            4.0 * ta_xxxy_yyyyy_0[i] * fe_0 - 4.0 * ta_xxxy_yyyyy_1[i] * fe_0 + ta_xxxxy_yyyyy_0[i] * pa_x[i] - ta_xxxxy_yyyyy_1[i] * pc_x[i];

        ta_xxxxxy_yyyyz_0[i] =
            4.0 * ta_xxxy_yyyyz_0[i] * fe_0 - 4.0 * ta_xxxy_yyyyz_1[i] * fe_0 + ta_xxxxy_yyyyz_0[i] * pa_x[i] - ta_xxxxy_yyyyz_1[i] * pc_x[i];

        ta_xxxxxy_yyyzz_0[i] =
            4.0 * ta_xxxy_yyyzz_0[i] * fe_0 - 4.0 * ta_xxxy_yyyzz_1[i] * fe_0 + ta_xxxxy_yyyzz_0[i] * pa_x[i] - ta_xxxxy_yyyzz_1[i] * pc_x[i];

        ta_xxxxxy_yyzzz_0[i] =
            4.0 * ta_xxxy_yyzzz_0[i] * fe_0 - 4.0 * ta_xxxy_yyzzz_1[i] * fe_0 + ta_xxxxy_yyzzz_0[i] * pa_x[i] - ta_xxxxy_yyzzz_1[i] * pc_x[i];

        ta_xxxxxy_yzzzz_0[i] =
            4.0 * ta_xxxy_yzzzz_0[i] * fe_0 - 4.0 * ta_xxxy_yzzzz_1[i] * fe_0 + ta_xxxxy_yzzzz_0[i] * pa_x[i] - ta_xxxxy_yzzzz_1[i] * pc_x[i];

        ta_xxxxxy_zzzzz_0[i] = ta_xxxxx_zzzzz_0[i] * pa_y[i] - ta_xxxxx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : IH

    auto ta_xxxxxz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 42);

    auto ta_xxxxxz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 43);

    auto ta_xxxxxz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 44);

    auto ta_xxxxxz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 45);

    auto ta_xxxxxz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 46);

    auto ta_xxxxxz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 47);

    auto ta_xxxxxz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 48);

    auto ta_xxxxxz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 49);

    auto ta_xxxxxz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 50);

    auto ta_xxxxxz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 51);

    auto ta_xxxxxz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 52);

    auto ta_xxxxxz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 53);

    auto ta_xxxxxz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 54);

    auto ta_xxxxxz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 55);

    auto ta_xxxxxz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 56);

    auto ta_xxxxxz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 57);

    auto ta_xxxxxz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 58);

    auto ta_xxxxxz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 59);

    auto ta_xxxxxz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 60);

    auto ta_xxxxxz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 61);

    auto ta_xxxxxz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 62);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxxxx_xxxx_0,   \
                             ta_xxxxx_xxxx_1,   \
                             ta_xxxxx_xxxxx_0,  \
                             ta_xxxxx_xxxxx_1,  \
                             ta_xxxxx_xxxxy_0,  \
                             ta_xxxxx_xxxxy_1,  \
                             ta_xxxxx_xxxxz_0,  \
                             ta_xxxxx_xxxxz_1,  \
                             ta_xxxxx_xxxy_0,   \
                             ta_xxxxx_xxxy_1,   \
                             ta_xxxxx_xxxyy_0,  \
                             ta_xxxxx_xxxyy_1,  \
                             ta_xxxxx_xxxyz_0,  \
                             ta_xxxxx_xxxyz_1,  \
                             ta_xxxxx_xxxz_0,   \
                             ta_xxxxx_xxxz_1,   \
                             ta_xxxxx_xxxzz_0,  \
                             ta_xxxxx_xxxzz_1,  \
                             ta_xxxxx_xxyy_0,   \
                             ta_xxxxx_xxyy_1,   \
                             ta_xxxxx_xxyyy_0,  \
                             ta_xxxxx_xxyyy_1,  \
                             ta_xxxxx_xxyyz_0,  \
                             ta_xxxxx_xxyyz_1,  \
                             ta_xxxxx_xxyz_0,   \
                             ta_xxxxx_xxyz_1,   \
                             ta_xxxxx_xxyzz_0,  \
                             ta_xxxxx_xxyzz_1,  \
                             ta_xxxxx_xxzz_0,   \
                             ta_xxxxx_xxzz_1,   \
                             ta_xxxxx_xxzzz_0,  \
                             ta_xxxxx_xxzzz_1,  \
                             ta_xxxxx_xyyy_0,   \
                             ta_xxxxx_xyyy_1,   \
                             ta_xxxxx_xyyyy_0,  \
                             ta_xxxxx_xyyyy_1,  \
                             ta_xxxxx_xyyyz_0,  \
                             ta_xxxxx_xyyyz_1,  \
                             ta_xxxxx_xyyz_0,   \
                             ta_xxxxx_xyyz_1,   \
                             ta_xxxxx_xyyzz_0,  \
                             ta_xxxxx_xyyzz_1,  \
                             ta_xxxxx_xyzz_0,   \
                             ta_xxxxx_xyzz_1,   \
                             ta_xxxxx_xyzzz_0,  \
                             ta_xxxxx_xyzzz_1,  \
                             ta_xxxxx_xzzz_0,   \
                             ta_xxxxx_xzzz_1,   \
                             ta_xxxxx_xzzzz_0,  \
                             ta_xxxxx_xzzzz_1,  \
                             ta_xxxxx_yyyyy_0,  \
                             ta_xxxxx_yyyyy_1,  \
                             ta_xxxxxz_xxxxx_0, \
                             ta_xxxxxz_xxxxy_0, \
                             ta_xxxxxz_xxxxz_0, \
                             ta_xxxxxz_xxxyy_0, \
                             ta_xxxxxz_xxxyz_0, \
                             ta_xxxxxz_xxxzz_0, \
                             ta_xxxxxz_xxyyy_0, \
                             ta_xxxxxz_xxyyz_0, \
                             ta_xxxxxz_xxyzz_0, \
                             ta_xxxxxz_xxzzz_0, \
                             ta_xxxxxz_xyyyy_0, \
                             ta_xxxxxz_xyyyz_0, \
                             ta_xxxxxz_xyyzz_0, \
                             ta_xxxxxz_xyzzz_0, \
                             ta_xxxxxz_xzzzz_0, \
                             ta_xxxxxz_yyyyy_0, \
                             ta_xxxxxz_yyyyz_0, \
                             ta_xxxxxz_yyyzz_0, \
                             ta_xxxxxz_yyzzz_0, \
                             ta_xxxxxz_yzzzz_0, \
                             ta_xxxxxz_zzzzz_0, \
                             ta_xxxxz_yyyyz_0,  \
                             ta_xxxxz_yyyyz_1,  \
                             ta_xxxxz_yyyzz_0,  \
                             ta_xxxxz_yyyzz_1,  \
                             ta_xxxxz_yyzzz_0,  \
                             ta_xxxxz_yyzzz_1,  \
                             ta_xxxxz_yzzzz_0,  \
                             ta_xxxxz_yzzzz_1,  \
                             ta_xxxxz_zzzzz_0,  \
                             ta_xxxxz_zzzzz_1,  \
                             ta_xxxz_yyyyz_0,   \
                             ta_xxxz_yyyyz_1,   \
                             ta_xxxz_yyyzz_0,   \
                             ta_xxxz_yyyzz_1,   \
                             ta_xxxz_yyzzz_0,   \
                             ta_xxxz_yyzzz_1,   \
                             ta_xxxz_yzzzz_0,   \
                             ta_xxxz_yzzzz_1,   \
                             ta_xxxz_zzzzz_0,   \
                             ta_xxxz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxz_xxxxx_0[i] = ta_xxxxx_xxxxx_0[i] * pa_z[i] - ta_xxxxx_xxxxx_1[i] * pc_z[i];

        ta_xxxxxz_xxxxy_0[i] = ta_xxxxx_xxxxy_0[i] * pa_z[i] - ta_xxxxx_xxxxy_1[i] * pc_z[i];

        ta_xxxxxz_xxxxz_0[i] = ta_xxxxx_xxxx_0[i] * fe_0 - ta_xxxxx_xxxx_1[i] * fe_0 + ta_xxxxx_xxxxz_0[i] * pa_z[i] - ta_xxxxx_xxxxz_1[i] * pc_z[i];

        ta_xxxxxz_xxxyy_0[i] = ta_xxxxx_xxxyy_0[i] * pa_z[i] - ta_xxxxx_xxxyy_1[i] * pc_z[i];

        ta_xxxxxz_xxxyz_0[i] = ta_xxxxx_xxxy_0[i] * fe_0 - ta_xxxxx_xxxy_1[i] * fe_0 + ta_xxxxx_xxxyz_0[i] * pa_z[i] - ta_xxxxx_xxxyz_1[i] * pc_z[i];

        ta_xxxxxz_xxxzz_0[i] =
            2.0 * ta_xxxxx_xxxz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxxz_1[i] * fe_0 + ta_xxxxx_xxxzz_0[i] * pa_z[i] - ta_xxxxx_xxxzz_1[i] * pc_z[i];

        ta_xxxxxz_xxyyy_0[i] = ta_xxxxx_xxyyy_0[i] * pa_z[i] - ta_xxxxx_xxyyy_1[i] * pc_z[i];

        ta_xxxxxz_xxyyz_0[i] = ta_xxxxx_xxyy_0[i] * fe_0 - ta_xxxxx_xxyy_1[i] * fe_0 + ta_xxxxx_xxyyz_0[i] * pa_z[i] - ta_xxxxx_xxyyz_1[i] * pc_z[i];

        ta_xxxxxz_xxyzz_0[i] =
            2.0 * ta_xxxxx_xxyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xxyz_1[i] * fe_0 + ta_xxxxx_xxyzz_0[i] * pa_z[i] - ta_xxxxx_xxyzz_1[i] * pc_z[i];

        ta_xxxxxz_xxzzz_0[i] =
            3.0 * ta_xxxxx_xxzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xxzz_1[i] * fe_0 + ta_xxxxx_xxzzz_0[i] * pa_z[i] - ta_xxxxx_xxzzz_1[i] * pc_z[i];

        ta_xxxxxz_xyyyy_0[i] = ta_xxxxx_xyyyy_0[i] * pa_z[i] - ta_xxxxx_xyyyy_1[i] * pc_z[i];

        ta_xxxxxz_xyyyz_0[i] = ta_xxxxx_xyyy_0[i] * fe_0 - ta_xxxxx_xyyy_1[i] * fe_0 + ta_xxxxx_xyyyz_0[i] * pa_z[i] - ta_xxxxx_xyyyz_1[i] * pc_z[i];

        ta_xxxxxz_xyyzz_0[i] =
            2.0 * ta_xxxxx_xyyz_0[i] * fe_0 - 2.0 * ta_xxxxx_xyyz_1[i] * fe_0 + ta_xxxxx_xyyzz_0[i] * pa_z[i] - ta_xxxxx_xyyzz_1[i] * pc_z[i];

        ta_xxxxxz_xyzzz_0[i] =
            3.0 * ta_xxxxx_xyzz_0[i] * fe_0 - 3.0 * ta_xxxxx_xyzz_1[i] * fe_0 + ta_xxxxx_xyzzz_0[i] * pa_z[i] - ta_xxxxx_xyzzz_1[i] * pc_z[i];

        ta_xxxxxz_xzzzz_0[i] =
            4.0 * ta_xxxxx_xzzz_0[i] * fe_0 - 4.0 * ta_xxxxx_xzzz_1[i] * fe_0 + ta_xxxxx_xzzzz_0[i] * pa_z[i] - ta_xxxxx_xzzzz_1[i] * pc_z[i];

        ta_xxxxxz_yyyyy_0[i] = ta_xxxxx_yyyyy_0[i] * pa_z[i] - ta_xxxxx_yyyyy_1[i] * pc_z[i];

        ta_xxxxxz_yyyyz_0[i] =
            4.0 * ta_xxxz_yyyyz_0[i] * fe_0 - 4.0 * ta_xxxz_yyyyz_1[i] * fe_0 + ta_xxxxz_yyyyz_0[i] * pa_x[i] - ta_xxxxz_yyyyz_1[i] * pc_x[i];

        ta_xxxxxz_yyyzz_0[i] =
            4.0 * ta_xxxz_yyyzz_0[i] * fe_0 - 4.0 * ta_xxxz_yyyzz_1[i] * fe_0 + ta_xxxxz_yyyzz_0[i] * pa_x[i] - ta_xxxxz_yyyzz_1[i] * pc_x[i];

        ta_xxxxxz_yyzzz_0[i] =
            4.0 * ta_xxxz_yyzzz_0[i] * fe_0 - 4.0 * ta_xxxz_yyzzz_1[i] * fe_0 + ta_xxxxz_yyzzz_0[i] * pa_x[i] - ta_xxxxz_yyzzz_1[i] * pc_x[i];

        ta_xxxxxz_yzzzz_0[i] =
            4.0 * ta_xxxz_yzzzz_0[i] * fe_0 - 4.0 * ta_xxxz_yzzzz_1[i] * fe_0 + ta_xxxxz_yzzzz_0[i] * pa_x[i] - ta_xxxxz_yzzzz_1[i] * pc_x[i];

        ta_xxxxxz_zzzzz_0[i] =
            4.0 * ta_xxxz_zzzzz_0[i] * fe_0 - 4.0 * ta_xxxz_zzzzz_1[i] * fe_0 + ta_xxxxz_zzzzz_0[i] * pa_x[i] - ta_xxxxz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 63-84 components of targeted buffer : IH

    auto ta_xxxxyy_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 63);

    auto ta_xxxxyy_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 64);

    auto ta_xxxxyy_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 65);

    auto ta_xxxxyy_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 66);

    auto ta_xxxxyy_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 67);

    auto ta_xxxxyy_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 68);

    auto ta_xxxxyy_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 69);

    auto ta_xxxxyy_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 70);

    auto ta_xxxxyy_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 71);

    auto ta_xxxxyy_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 72);

    auto ta_xxxxyy_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 73);

    auto ta_xxxxyy_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 74);

    auto ta_xxxxyy_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 75);

    auto ta_xxxxyy_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 76);

    auto ta_xxxxyy_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 77);

    auto ta_xxxxyy_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 78);

    auto ta_xxxxyy_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 79);

    auto ta_xxxxyy_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 80);

    auto ta_xxxxyy_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 81);

    auto ta_xxxxyy_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 82);

    auto ta_xxxxyy_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 83);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxxx_xxxxx_0,   \
                             ta_xxxx_xxxxx_1,   \
                             ta_xxxx_xxxxz_0,   \
                             ta_xxxx_xxxxz_1,   \
                             ta_xxxx_xxxzz_0,   \
                             ta_xxxx_xxxzz_1,   \
                             ta_xxxx_xxzzz_0,   \
                             ta_xxxx_xxzzz_1,   \
                             ta_xxxx_xzzzz_0,   \
                             ta_xxxx_xzzzz_1,   \
                             ta_xxxxy_xxxxx_0,  \
                             ta_xxxxy_xxxxx_1,  \
                             ta_xxxxy_xxxxz_0,  \
                             ta_xxxxy_xxxxz_1,  \
                             ta_xxxxy_xxxzz_0,  \
                             ta_xxxxy_xxxzz_1,  \
                             ta_xxxxy_xxzzz_0,  \
                             ta_xxxxy_xxzzz_1,  \
                             ta_xxxxy_xzzzz_0,  \
                             ta_xxxxy_xzzzz_1,  \
                             ta_xxxxyy_xxxxx_0, \
                             ta_xxxxyy_xxxxy_0, \
                             ta_xxxxyy_xxxxz_0, \
                             ta_xxxxyy_xxxyy_0, \
                             ta_xxxxyy_xxxyz_0, \
                             ta_xxxxyy_xxxzz_0, \
                             ta_xxxxyy_xxyyy_0, \
                             ta_xxxxyy_xxyyz_0, \
                             ta_xxxxyy_xxyzz_0, \
                             ta_xxxxyy_xxzzz_0, \
                             ta_xxxxyy_xyyyy_0, \
                             ta_xxxxyy_xyyyz_0, \
                             ta_xxxxyy_xyyzz_0, \
                             ta_xxxxyy_xyzzz_0, \
                             ta_xxxxyy_xzzzz_0, \
                             ta_xxxxyy_yyyyy_0, \
                             ta_xxxxyy_yyyyz_0, \
                             ta_xxxxyy_yyyzz_0, \
                             ta_xxxxyy_yyzzz_0, \
                             ta_xxxxyy_yzzzz_0, \
                             ta_xxxxyy_zzzzz_0, \
                             ta_xxxyy_xxxxy_0,  \
                             ta_xxxyy_xxxxy_1,  \
                             ta_xxxyy_xxxy_0,   \
                             ta_xxxyy_xxxy_1,   \
                             ta_xxxyy_xxxyy_0,  \
                             ta_xxxyy_xxxyy_1,  \
                             ta_xxxyy_xxxyz_0,  \
                             ta_xxxyy_xxxyz_1,  \
                             ta_xxxyy_xxyy_0,   \
                             ta_xxxyy_xxyy_1,   \
                             ta_xxxyy_xxyyy_0,  \
                             ta_xxxyy_xxyyy_1,  \
                             ta_xxxyy_xxyyz_0,  \
                             ta_xxxyy_xxyyz_1,  \
                             ta_xxxyy_xxyz_0,   \
                             ta_xxxyy_xxyz_1,   \
                             ta_xxxyy_xxyzz_0,  \
                             ta_xxxyy_xxyzz_1,  \
                             ta_xxxyy_xyyy_0,   \
                             ta_xxxyy_xyyy_1,   \
                             ta_xxxyy_xyyyy_0,  \
                             ta_xxxyy_xyyyy_1,  \
                             ta_xxxyy_xyyyz_0,  \
                             ta_xxxyy_xyyyz_1,  \
                             ta_xxxyy_xyyz_0,   \
                             ta_xxxyy_xyyz_1,   \
                             ta_xxxyy_xyyzz_0,  \
                             ta_xxxyy_xyyzz_1,  \
                             ta_xxxyy_xyzz_0,   \
                             ta_xxxyy_xyzz_1,   \
                             ta_xxxyy_xyzzz_0,  \
                             ta_xxxyy_xyzzz_1,  \
                             ta_xxxyy_yyyy_0,   \
                             ta_xxxyy_yyyy_1,   \
                             ta_xxxyy_yyyyy_0,  \
                             ta_xxxyy_yyyyy_1,  \
                             ta_xxxyy_yyyyz_0,  \
                             ta_xxxyy_yyyyz_1,  \
                             ta_xxxyy_yyyz_0,   \
                             ta_xxxyy_yyyz_1,   \
                             ta_xxxyy_yyyzz_0,  \
                             ta_xxxyy_yyyzz_1,  \
                             ta_xxxyy_yyzz_0,   \
                             ta_xxxyy_yyzz_1,   \
                             ta_xxxyy_yyzzz_0,  \
                             ta_xxxyy_yyzzz_1,  \
                             ta_xxxyy_yzzz_0,   \
                             ta_xxxyy_yzzz_1,   \
                             ta_xxxyy_yzzzz_0,  \
                             ta_xxxyy_yzzzz_1,  \
                             ta_xxxyy_zzzzz_0,  \
                             ta_xxxyy_zzzzz_1,  \
                             ta_xxyy_xxxxy_0,   \
                             ta_xxyy_xxxxy_1,   \
                             ta_xxyy_xxxyy_0,   \
                             ta_xxyy_xxxyy_1,   \
                             ta_xxyy_xxxyz_0,   \
                             ta_xxyy_xxxyz_1,   \
                             ta_xxyy_xxyyy_0,   \
                             ta_xxyy_xxyyy_1,   \
                             ta_xxyy_xxyyz_0,   \
                             ta_xxyy_xxyyz_1,   \
                             ta_xxyy_xxyzz_0,   \
                             ta_xxyy_xxyzz_1,   \
                             ta_xxyy_xyyyy_0,   \
                             ta_xxyy_xyyyy_1,   \
                             ta_xxyy_xyyyz_0,   \
                             ta_xxyy_xyyyz_1,   \
                             ta_xxyy_xyyzz_0,   \
                             ta_xxyy_xyyzz_1,   \
                             ta_xxyy_xyzzz_0,   \
                             ta_xxyy_xyzzz_1,   \
                             ta_xxyy_yyyyy_0,   \
                             ta_xxyy_yyyyy_1,   \
                             ta_xxyy_yyyyz_0,   \
                             ta_xxyy_yyyyz_1,   \
                             ta_xxyy_yyyzz_0,   \
                             ta_xxyy_yyyzz_1,   \
                             ta_xxyy_yyzzz_0,   \
                             ta_xxyy_yyzzz_1,   \
                             ta_xxyy_yzzzz_0,   \
                             ta_xxyy_yzzzz_1,   \
                             ta_xxyy_zzzzz_0,   \
                             ta_xxyy_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyy_xxxxx_0[i] = ta_xxxx_xxxxx_0[i] * fe_0 - ta_xxxx_xxxxx_1[i] * fe_0 + ta_xxxxy_xxxxx_0[i] * pa_y[i] - ta_xxxxy_xxxxx_1[i] * pc_y[i];

        ta_xxxxyy_xxxxy_0[i] = 3.0 * ta_xxyy_xxxxy_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxy_1[i] * fe_0 + 4.0 * ta_xxxyy_xxxy_0[i] * fe_0 -
                               4.0 * ta_xxxyy_xxxy_1[i] * fe_0 + ta_xxxyy_xxxxy_0[i] * pa_x[i] - ta_xxxyy_xxxxy_1[i] * pc_x[i];

        ta_xxxxyy_xxxxz_0[i] = ta_xxxx_xxxxz_0[i] * fe_0 - ta_xxxx_xxxxz_1[i] * fe_0 + ta_xxxxy_xxxxz_0[i] * pa_y[i] - ta_xxxxy_xxxxz_1[i] * pc_y[i];

        ta_xxxxyy_xxxyy_0[i] = 3.0 * ta_xxyy_xxxyy_0[i] * fe_0 - 3.0 * ta_xxyy_xxxyy_1[i] * fe_0 + 3.0 * ta_xxxyy_xxyy_0[i] * fe_0 -
                               3.0 * ta_xxxyy_xxyy_1[i] * fe_0 + ta_xxxyy_xxxyy_0[i] * pa_x[i] - ta_xxxyy_xxxyy_1[i] * pc_x[i];

        ta_xxxxyy_xxxyz_0[i] = 3.0 * ta_xxyy_xxxyz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxyz_1[i] * fe_0 + 3.0 * ta_xxxyy_xxyz_0[i] * fe_0 -
                               3.0 * ta_xxxyy_xxyz_1[i] * fe_0 + ta_xxxyy_xxxyz_0[i] * pa_x[i] - ta_xxxyy_xxxyz_1[i] * pc_x[i];

        ta_xxxxyy_xxxzz_0[i] = ta_xxxx_xxxzz_0[i] * fe_0 - ta_xxxx_xxxzz_1[i] * fe_0 + ta_xxxxy_xxxzz_0[i] * pa_y[i] - ta_xxxxy_xxxzz_1[i] * pc_y[i];

        ta_xxxxyy_xxyyy_0[i] = 3.0 * ta_xxyy_xxyyy_0[i] * fe_0 - 3.0 * ta_xxyy_xxyyy_1[i] * fe_0 + 2.0 * ta_xxxyy_xyyy_0[i] * fe_0 -
                               2.0 * ta_xxxyy_xyyy_1[i] * fe_0 + ta_xxxyy_xxyyy_0[i] * pa_x[i] - ta_xxxyy_xxyyy_1[i] * pc_x[i];

        ta_xxxxyy_xxyyz_0[i] = 3.0 * ta_xxyy_xxyyz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyyz_1[i] * fe_0 + 2.0 * ta_xxxyy_xyyz_0[i] * fe_0 -
                               2.0 * ta_xxxyy_xyyz_1[i] * fe_0 + ta_xxxyy_xxyyz_0[i] * pa_x[i] - ta_xxxyy_xxyyz_1[i] * pc_x[i];

        ta_xxxxyy_xxyzz_0[i] = 3.0 * ta_xxyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyzz_1[i] * fe_0 + 2.0 * ta_xxxyy_xyzz_0[i] * fe_0 -
                               2.0 * ta_xxxyy_xyzz_1[i] * fe_0 + ta_xxxyy_xxyzz_0[i] * pa_x[i] - ta_xxxyy_xxyzz_1[i] * pc_x[i];

        ta_xxxxyy_xxzzz_0[i] = ta_xxxx_xxzzz_0[i] * fe_0 - ta_xxxx_xxzzz_1[i] * fe_0 + ta_xxxxy_xxzzz_0[i] * pa_y[i] - ta_xxxxy_xxzzz_1[i] * pc_y[i];

        ta_xxxxyy_xyyyy_0[i] = 3.0 * ta_xxyy_xyyyy_0[i] * fe_0 - 3.0 * ta_xxyy_xyyyy_1[i] * fe_0 + ta_xxxyy_yyyy_0[i] * fe_0 -
                               ta_xxxyy_yyyy_1[i] * fe_0 + ta_xxxyy_xyyyy_0[i] * pa_x[i] - ta_xxxyy_xyyyy_1[i] * pc_x[i];

        ta_xxxxyy_xyyyz_0[i] = 3.0 * ta_xxyy_xyyyz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyyz_1[i] * fe_0 + ta_xxxyy_yyyz_0[i] * fe_0 -
                               ta_xxxyy_yyyz_1[i] * fe_0 + ta_xxxyy_xyyyz_0[i] * pa_x[i] - ta_xxxyy_xyyyz_1[i] * pc_x[i];

        ta_xxxxyy_xyyzz_0[i] = 3.0 * ta_xxyy_xyyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyzz_1[i] * fe_0 + ta_xxxyy_yyzz_0[i] * fe_0 -
                               ta_xxxyy_yyzz_1[i] * fe_0 + ta_xxxyy_xyyzz_0[i] * pa_x[i] - ta_xxxyy_xyyzz_1[i] * pc_x[i];

        ta_xxxxyy_xyzzz_0[i] = 3.0 * ta_xxyy_xyzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyzzz_1[i] * fe_0 + ta_xxxyy_yzzz_0[i] * fe_0 -
                               ta_xxxyy_yzzz_1[i] * fe_0 + ta_xxxyy_xyzzz_0[i] * pa_x[i] - ta_xxxyy_xyzzz_1[i] * pc_x[i];

        ta_xxxxyy_xzzzz_0[i] = ta_xxxx_xzzzz_0[i] * fe_0 - ta_xxxx_xzzzz_1[i] * fe_0 + ta_xxxxy_xzzzz_0[i] * pa_y[i] - ta_xxxxy_xzzzz_1[i] * pc_y[i];

        ta_xxxxyy_yyyyy_0[i] =
            3.0 * ta_xxyy_yyyyy_0[i] * fe_0 - 3.0 * ta_xxyy_yyyyy_1[i] * fe_0 + ta_xxxyy_yyyyy_0[i] * pa_x[i] - ta_xxxyy_yyyyy_1[i] * pc_x[i];

        ta_xxxxyy_yyyyz_0[i] =
            3.0 * ta_xxyy_yyyyz_0[i] * fe_0 - 3.0 * ta_xxyy_yyyyz_1[i] * fe_0 + ta_xxxyy_yyyyz_0[i] * pa_x[i] - ta_xxxyy_yyyyz_1[i] * pc_x[i];

        ta_xxxxyy_yyyzz_0[i] =
            3.0 * ta_xxyy_yyyzz_0[i] * fe_0 - 3.0 * ta_xxyy_yyyzz_1[i] * fe_0 + ta_xxxyy_yyyzz_0[i] * pa_x[i] - ta_xxxyy_yyyzz_1[i] * pc_x[i];

        ta_xxxxyy_yyzzz_0[i] =
            3.0 * ta_xxyy_yyzzz_0[i] * fe_0 - 3.0 * ta_xxyy_yyzzz_1[i] * fe_0 + ta_xxxyy_yyzzz_0[i] * pa_x[i] - ta_xxxyy_yyzzz_1[i] * pc_x[i];

        ta_xxxxyy_yzzzz_0[i] =
            3.0 * ta_xxyy_yzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_yzzzz_1[i] * fe_0 + ta_xxxyy_yzzzz_0[i] * pa_x[i] - ta_xxxyy_yzzzz_1[i] * pc_x[i];

        ta_xxxxyy_zzzzz_0[i] =
            3.0 * ta_xxyy_zzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_zzzzz_1[i] * fe_0 + ta_xxxyy_zzzzz_0[i] * pa_x[i] - ta_xxxyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 84-105 components of targeted buffer : IH

    auto ta_xxxxyz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 84);

    auto ta_xxxxyz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 85);

    auto ta_xxxxyz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 86);

    auto ta_xxxxyz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 87);

    auto ta_xxxxyz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 88);

    auto ta_xxxxyz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 89);

    auto ta_xxxxyz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 90);

    auto ta_xxxxyz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 91);

    auto ta_xxxxyz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 92);

    auto ta_xxxxyz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 93);

    auto ta_xxxxyz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 94);

    auto ta_xxxxyz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 95);

    auto ta_xxxxyz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 96);

    auto ta_xxxxyz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 97);

    auto ta_xxxxyz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 98);

    auto ta_xxxxyz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 99);

    auto ta_xxxxyz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 100);

    auto ta_xxxxyz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 101);

    auto ta_xxxxyz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 102);

    auto ta_xxxxyz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 103);

    auto ta_xxxxyz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 104);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta_xxxxy_xxxxy_0,  \
                             ta_xxxxy_xxxxy_1,  \
                             ta_xxxxy_xxxyy_0,  \
                             ta_xxxxy_xxxyy_1,  \
                             ta_xxxxy_xxyyy_0,  \
                             ta_xxxxy_xxyyy_1,  \
                             ta_xxxxy_xyyyy_0,  \
                             ta_xxxxy_xyyyy_1,  \
                             ta_xxxxy_yyyyy_0,  \
                             ta_xxxxy_yyyyy_1,  \
                             ta_xxxxyz_xxxxx_0, \
                             ta_xxxxyz_xxxxy_0, \
                             ta_xxxxyz_xxxxz_0, \
                             ta_xxxxyz_xxxyy_0, \
                             ta_xxxxyz_xxxyz_0, \
                             ta_xxxxyz_xxxzz_0, \
                             ta_xxxxyz_xxyyy_0, \
                             ta_xxxxyz_xxyyz_0, \
                             ta_xxxxyz_xxyzz_0, \
                             ta_xxxxyz_xxzzz_0, \
                             ta_xxxxyz_xyyyy_0, \
                             ta_xxxxyz_xyyyz_0, \
                             ta_xxxxyz_xyyzz_0, \
                             ta_xxxxyz_xyzzz_0, \
                             ta_xxxxyz_xzzzz_0, \
                             ta_xxxxyz_yyyyy_0, \
                             ta_xxxxyz_yyyyz_0, \
                             ta_xxxxyz_yyyzz_0, \
                             ta_xxxxyz_yyzzz_0, \
                             ta_xxxxyz_yzzzz_0, \
                             ta_xxxxyz_zzzzz_0, \
                             ta_xxxxz_xxxxx_0,  \
                             ta_xxxxz_xxxxx_1,  \
                             ta_xxxxz_xxxxz_0,  \
                             ta_xxxxz_xxxxz_1,  \
                             ta_xxxxz_xxxyz_0,  \
                             ta_xxxxz_xxxyz_1,  \
                             ta_xxxxz_xxxz_0,   \
                             ta_xxxxz_xxxz_1,   \
                             ta_xxxxz_xxxzz_0,  \
                             ta_xxxxz_xxxzz_1,  \
                             ta_xxxxz_xxyyz_0,  \
                             ta_xxxxz_xxyyz_1,  \
                             ta_xxxxz_xxyz_0,   \
                             ta_xxxxz_xxyz_1,   \
                             ta_xxxxz_xxyzz_0,  \
                             ta_xxxxz_xxyzz_1,  \
                             ta_xxxxz_xxzz_0,   \
                             ta_xxxxz_xxzz_1,   \
                             ta_xxxxz_xxzzz_0,  \
                             ta_xxxxz_xxzzz_1,  \
                             ta_xxxxz_xyyyz_0,  \
                             ta_xxxxz_xyyyz_1,  \
                             ta_xxxxz_xyyz_0,   \
                             ta_xxxxz_xyyz_1,   \
                             ta_xxxxz_xyyzz_0,  \
                             ta_xxxxz_xyyzz_1,  \
                             ta_xxxxz_xyzz_0,   \
                             ta_xxxxz_xyzz_1,   \
                             ta_xxxxz_xyzzz_0,  \
                             ta_xxxxz_xyzzz_1,  \
                             ta_xxxxz_xzzz_0,   \
                             ta_xxxxz_xzzz_1,   \
                             ta_xxxxz_xzzzz_0,  \
                             ta_xxxxz_xzzzz_1,  \
                             ta_xxxxz_zzzzz_0,  \
                             ta_xxxxz_zzzzz_1,  \
                             ta_xxxyz_yyyyz_0,  \
                             ta_xxxyz_yyyyz_1,  \
                             ta_xxxyz_yyyzz_0,  \
                             ta_xxxyz_yyyzz_1,  \
                             ta_xxxyz_yyzzz_0,  \
                             ta_xxxyz_yyzzz_1,  \
                             ta_xxxyz_yzzzz_0,  \
                             ta_xxxyz_yzzzz_1,  \
                             ta_xxyz_yyyyz_0,   \
                             ta_xxyz_yyyyz_1,   \
                             ta_xxyz_yyyzz_0,   \
                             ta_xxyz_yyyzz_1,   \
                             ta_xxyz_yyzzz_0,   \
                             ta_xxyz_yyzzz_1,   \
                             ta_xxyz_yzzzz_0,   \
                             ta_xxyz_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyz_xxxxx_0[i] = ta_xxxxz_xxxxx_0[i] * pa_y[i] - ta_xxxxz_xxxxx_1[i] * pc_y[i];

        ta_xxxxyz_xxxxy_0[i] = ta_xxxxy_xxxxy_0[i] * pa_z[i] - ta_xxxxy_xxxxy_1[i] * pc_z[i];

        ta_xxxxyz_xxxxz_0[i] = ta_xxxxz_xxxxz_0[i] * pa_y[i] - ta_xxxxz_xxxxz_1[i] * pc_y[i];

        ta_xxxxyz_xxxyy_0[i] = ta_xxxxy_xxxyy_0[i] * pa_z[i] - ta_xxxxy_xxxyy_1[i] * pc_z[i];

        ta_xxxxyz_xxxyz_0[i] = ta_xxxxz_xxxz_0[i] * fe_0 - ta_xxxxz_xxxz_1[i] * fe_0 + ta_xxxxz_xxxyz_0[i] * pa_y[i] - ta_xxxxz_xxxyz_1[i] * pc_y[i];

        ta_xxxxyz_xxxzz_0[i] = ta_xxxxz_xxxzz_0[i] * pa_y[i] - ta_xxxxz_xxxzz_1[i] * pc_y[i];

        ta_xxxxyz_xxyyy_0[i] = ta_xxxxy_xxyyy_0[i] * pa_z[i] - ta_xxxxy_xxyyy_1[i] * pc_z[i];

        ta_xxxxyz_xxyyz_0[i] =
            2.0 * ta_xxxxz_xxyz_0[i] * fe_0 - 2.0 * ta_xxxxz_xxyz_1[i] * fe_0 + ta_xxxxz_xxyyz_0[i] * pa_y[i] - ta_xxxxz_xxyyz_1[i] * pc_y[i];

        ta_xxxxyz_xxyzz_0[i] = ta_xxxxz_xxzz_0[i] * fe_0 - ta_xxxxz_xxzz_1[i] * fe_0 + ta_xxxxz_xxyzz_0[i] * pa_y[i] - ta_xxxxz_xxyzz_1[i] * pc_y[i];

        ta_xxxxyz_xxzzz_0[i] = ta_xxxxz_xxzzz_0[i] * pa_y[i] - ta_xxxxz_xxzzz_1[i] * pc_y[i];

        ta_xxxxyz_xyyyy_0[i] = ta_xxxxy_xyyyy_0[i] * pa_z[i] - ta_xxxxy_xyyyy_1[i] * pc_z[i];

        ta_xxxxyz_xyyyz_0[i] =
            3.0 * ta_xxxxz_xyyz_0[i] * fe_0 - 3.0 * ta_xxxxz_xyyz_1[i] * fe_0 + ta_xxxxz_xyyyz_0[i] * pa_y[i] - ta_xxxxz_xyyyz_1[i] * pc_y[i];

        ta_xxxxyz_xyyzz_0[i] =
            2.0 * ta_xxxxz_xyzz_0[i] * fe_0 - 2.0 * ta_xxxxz_xyzz_1[i] * fe_0 + ta_xxxxz_xyyzz_0[i] * pa_y[i] - ta_xxxxz_xyyzz_1[i] * pc_y[i];

        ta_xxxxyz_xyzzz_0[i] = ta_xxxxz_xzzz_0[i] * fe_0 - ta_xxxxz_xzzz_1[i] * fe_0 + ta_xxxxz_xyzzz_0[i] * pa_y[i] - ta_xxxxz_xyzzz_1[i] * pc_y[i];

        ta_xxxxyz_xzzzz_0[i] = ta_xxxxz_xzzzz_0[i] * pa_y[i] - ta_xxxxz_xzzzz_1[i] * pc_y[i];

        ta_xxxxyz_yyyyy_0[i] = ta_xxxxy_yyyyy_0[i] * pa_z[i] - ta_xxxxy_yyyyy_1[i] * pc_z[i];

        ta_xxxxyz_yyyyz_0[i] =
            3.0 * ta_xxyz_yyyyz_0[i] * fe_0 - 3.0 * ta_xxyz_yyyyz_1[i] * fe_0 + ta_xxxyz_yyyyz_0[i] * pa_x[i] - ta_xxxyz_yyyyz_1[i] * pc_x[i];

        ta_xxxxyz_yyyzz_0[i] =
            3.0 * ta_xxyz_yyyzz_0[i] * fe_0 - 3.0 * ta_xxyz_yyyzz_1[i] * fe_0 + ta_xxxyz_yyyzz_0[i] * pa_x[i] - ta_xxxyz_yyyzz_1[i] * pc_x[i];

        ta_xxxxyz_yyzzz_0[i] =
            3.0 * ta_xxyz_yyzzz_0[i] * fe_0 - 3.0 * ta_xxyz_yyzzz_1[i] * fe_0 + ta_xxxyz_yyzzz_0[i] * pa_x[i] - ta_xxxyz_yyzzz_1[i] * pc_x[i];

        ta_xxxxyz_yzzzz_0[i] =
            3.0 * ta_xxyz_yzzzz_0[i] * fe_0 - 3.0 * ta_xxyz_yzzzz_1[i] * fe_0 + ta_xxxyz_yzzzz_0[i] * pa_x[i] - ta_xxxyz_yzzzz_1[i] * pc_x[i];

        ta_xxxxyz_zzzzz_0[i] = ta_xxxxz_zzzzz_0[i] * pa_y[i] - ta_xxxxz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : IH

    auto ta_xxxxzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 105);

    auto ta_xxxxzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 106);

    auto ta_xxxxzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 107);

    auto ta_xxxxzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 108);

    auto ta_xxxxzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 109);

    auto ta_xxxxzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 110);

    auto ta_xxxxzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 111);

    auto ta_xxxxzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 112);

    auto ta_xxxxzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 113);

    auto ta_xxxxzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 114);

    auto ta_xxxxzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 115);

    auto ta_xxxxzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 116);

    auto ta_xxxxzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 117);

    auto ta_xxxxzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 118);

    auto ta_xxxxzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 119);

    auto ta_xxxxzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 120);

    auto ta_xxxxzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 121);

    auto ta_xxxxzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 122);

    auto ta_xxxxzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 123);

    auto ta_xxxxzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 124);

    auto ta_xxxxzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 125);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxxx_xxxxx_0,   \
                             ta_xxxx_xxxxx_1,   \
                             ta_xxxx_xxxxy_0,   \
                             ta_xxxx_xxxxy_1,   \
                             ta_xxxx_xxxyy_0,   \
                             ta_xxxx_xxxyy_1,   \
                             ta_xxxx_xxyyy_0,   \
                             ta_xxxx_xxyyy_1,   \
                             ta_xxxx_xyyyy_0,   \
                             ta_xxxx_xyyyy_1,   \
                             ta_xxxxz_xxxxx_0,  \
                             ta_xxxxz_xxxxx_1,  \
                             ta_xxxxz_xxxxy_0,  \
                             ta_xxxxz_xxxxy_1,  \
                             ta_xxxxz_xxxyy_0,  \
                             ta_xxxxz_xxxyy_1,  \
                             ta_xxxxz_xxyyy_0,  \
                             ta_xxxxz_xxyyy_1,  \
                             ta_xxxxz_xyyyy_0,  \
                             ta_xxxxz_xyyyy_1,  \
                             ta_xxxxzz_xxxxx_0, \
                             ta_xxxxzz_xxxxy_0, \
                             ta_xxxxzz_xxxxz_0, \
                             ta_xxxxzz_xxxyy_0, \
                             ta_xxxxzz_xxxyz_0, \
                             ta_xxxxzz_xxxzz_0, \
                             ta_xxxxzz_xxyyy_0, \
                             ta_xxxxzz_xxyyz_0, \
                             ta_xxxxzz_xxyzz_0, \
                             ta_xxxxzz_xxzzz_0, \
                             ta_xxxxzz_xyyyy_0, \
                             ta_xxxxzz_xyyyz_0, \
                             ta_xxxxzz_xyyzz_0, \
                             ta_xxxxzz_xyzzz_0, \
                             ta_xxxxzz_xzzzz_0, \
                             ta_xxxxzz_yyyyy_0, \
                             ta_xxxxzz_yyyyz_0, \
                             ta_xxxxzz_yyyzz_0, \
                             ta_xxxxzz_yyzzz_0, \
                             ta_xxxxzz_yzzzz_0, \
                             ta_xxxxzz_zzzzz_0, \
                             ta_xxxzz_xxxxz_0,  \
                             ta_xxxzz_xxxxz_1,  \
                             ta_xxxzz_xxxyz_0,  \
                             ta_xxxzz_xxxyz_1,  \
                             ta_xxxzz_xxxz_0,   \
                             ta_xxxzz_xxxz_1,   \
                             ta_xxxzz_xxxzz_0,  \
                             ta_xxxzz_xxxzz_1,  \
                             ta_xxxzz_xxyyz_0,  \
                             ta_xxxzz_xxyyz_1,  \
                             ta_xxxzz_xxyz_0,   \
                             ta_xxxzz_xxyz_1,   \
                             ta_xxxzz_xxyzz_0,  \
                             ta_xxxzz_xxyzz_1,  \
                             ta_xxxzz_xxzz_0,   \
                             ta_xxxzz_xxzz_1,   \
                             ta_xxxzz_xxzzz_0,  \
                             ta_xxxzz_xxzzz_1,  \
                             ta_xxxzz_xyyyz_0,  \
                             ta_xxxzz_xyyyz_1,  \
                             ta_xxxzz_xyyz_0,   \
                             ta_xxxzz_xyyz_1,   \
                             ta_xxxzz_xyyzz_0,  \
                             ta_xxxzz_xyyzz_1,  \
                             ta_xxxzz_xyzz_0,   \
                             ta_xxxzz_xyzz_1,   \
                             ta_xxxzz_xyzzz_0,  \
                             ta_xxxzz_xyzzz_1,  \
                             ta_xxxzz_xzzz_0,   \
                             ta_xxxzz_xzzz_1,   \
                             ta_xxxzz_xzzzz_0,  \
                             ta_xxxzz_xzzzz_1,  \
                             ta_xxxzz_yyyyy_0,  \
                             ta_xxxzz_yyyyy_1,  \
                             ta_xxxzz_yyyyz_0,  \
                             ta_xxxzz_yyyyz_1,  \
                             ta_xxxzz_yyyz_0,   \
                             ta_xxxzz_yyyz_1,   \
                             ta_xxxzz_yyyzz_0,  \
                             ta_xxxzz_yyyzz_1,  \
                             ta_xxxzz_yyzz_0,   \
                             ta_xxxzz_yyzz_1,   \
                             ta_xxxzz_yyzzz_0,  \
                             ta_xxxzz_yyzzz_1,  \
                             ta_xxxzz_yzzz_0,   \
                             ta_xxxzz_yzzz_1,   \
                             ta_xxxzz_yzzzz_0,  \
                             ta_xxxzz_yzzzz_1,  \
                             ta_xxxzz_zzzz_0,   \
                             ta_xxxzz_zzzz_1,   \
                             ta_xxxzz_zzzzz_0,  \
                             ta_xxxzz_zzzzz_1,  \
                             ta_xxzz_xxxxz_0,   \
                             ta_xxzz_xxxxz_1,   \
                             ta_xxzz_xxxyz_0,   \
                             ta_xxzz_xxxyz_1,   \
                             ta_xxzz_xxxzz_0,   \
                             ta_xxzz_xxxzz_1,   \
                             ta_xxzz_xxyyz_0,   \
                             ta_xxzz_xxyyz_1,   \
                             ta_xxzz_xxyzz_0,   \
                             ta_xxzz_xxyzz_1,   \
                             ta_xxzz_xxzzz_0,   \
                             ta_xxzz_xxzzz_1,   \
                             ta_xxzz_xyyyz_0,   \
                             ta_xxzz_xyyyz_1,   \
                             ta_xxzz_xyyzz_0,   \
                             ta_xxzz_xyyzz_1,   \
                             ta_xxzz_xyzzz_0,   \
                             ta_xxzz_xyzzz_1,   \
                             ta_xxzz_xzzzz_0,   \
                             ta_xxzz_xzzzz_1,   \
                             ta_xxzz_yyyyy_0,   \
                             ta_xxzz_yyyyy_1,   \
                             ta_xxzz_yyyyz_0,   \
                             ta_xxzz_yyyyz_1,   \
                             ta_xxzz_yyyzz_0,   \
                             ta_xxzz_yyyzz_1,   \
                             ta_xxzz_yyzzz_0,   \
                             ta_xxzz_yyzzz_1,   \
                             ta_xxzz_yzzzz_0,   \
                             ta_xxzz_yzzzz_1,   \
                             ta_xxzz_zzzzz_0,   \
                             ta_xxzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxzz_xxxxx_0[i] = ta_xxxx_xxxxx_0[i] * fe_0 - ta_xxxx_xxxxx_1[i] * fe_0 + ta_xxxxz_xxxxx_0[i] * pa_z[i] - ta_xxxxz_xxxxx_1[i] * pc_z[i];

        ta_xxxxzz_xxxxy_0[i] = ta_xxxx_xxxxy_0[i] * fe_0 - ta_xxxx_xxxxy_1[i] * fe_0 + ta_xxxxz_xxxxy_0[i] * pa_z[i] - ta_xxxxz_xxxxy_1[i] * pc_z[i];

        ta_xxxxzz_xxxxz_0[i] = 3.0 * ta_xxzz_xxxxz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxz_1[i] * fe_0 + 4.0 * ta_xxxzz_xxxz_0[i] * fe_0 -
                               4.0 * ta_xxxzz_xxxz_1[i] * fe_0 + ta_xxxzz_xxxxz_0[i] * pa_x[i] - ta_xxxzz_xxxxz_1[i] * pc_x[i];

        ta_xxxxzz_xxxyy_0[i] = ta_xxxx_xxxyy_0[i] * fe_0 - ta_xxxx_xxxyy_1[i] * fe_0 + ta_xxxxz_xxxyy_0[i] * pa_z[i] - ta_xxxxz_xxxyy_1[i] * pc_z[i];

        ta_xxxxzz_xxxyz_0[i] = 3.0 * ta_xxzz_xxxyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxyz_1[i] * fe_0 + 3.0 * ta_xxxzz_xxyz_0[i] * fe_0 -
                               3.0 * ta_xxxzz_xxyz_1[i] * fe_0 + ta_xxxzz_xxxyz_0[i] * pa_x[i] - ta_xxxzz_xxxyz_1[i] * pc_x[i];

        ta_xxxxzz_xxxzz_0[i] = 3.0 * ta_xxzz_xxxzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxxzz_1[i] * fe_0 + 3.0 * ta_xxxzz_xxzz_0[i] * fe_0 -
                               3.0 * ta_xxxzz_xxzz_1[i] * fe_0 + ta_xxxzz_xxxzz_0[i] * pa_x[i] - ta_xxxzz_xxxzz_1[i] * pc_x[i];

        ta_xxxxzz_xxyyy_0[i] = ta_xxxx_xxyyy_0[i] * fe_0 - ta_xxxx_xxyyy_1[i] * fe_0 + ta_xxxxz_xxyyy_0[i] * pa_z[i] - ta_xxxxz_xxyyy_1[i] * pc_z[i];

        ta_xxxxzz_xxyyz_0[i] = 3.0 * ta_xxzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyyz_1[i] * fe_0 + 2.0 * ta_xxxzz_xyyz_0[i] * fe_0 -
                               2.0 * ta_xxxzz_xyyz_1[i] * fe_0 + ta_xxxzz_xxyyz_0[i] * pa_x[i] - ta_xxxzz_xxyyz_1[i] * pc_x[i];

        ta_xxxxzz_xxyzz_0[i] = 3.0 * ta_xxzz_xxyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyzz_1[i] * fe_0 + 2.0 * ta_xxxzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_xxxzz_xyzz_1[i] * fe_0 + ta_xxxzz_xxyzz_0[i] * pa_x[i] - ta_xxxzz_xxyzz_1[i] * pc_x[i];

        ta_xxxxzz_xxzzz_0[i] = 3.0 * ta_xxzz_xxzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xxzzz_1[i] * fe_0 + 2.0 * ta_xxxzz_xzzz_0[i] * fe_0 -
                               2.0 * ta_xxxzz_xzzz_1[i] * fe_0 + ta_xxxzz_xxzzz_0[i] * pa_x[i] - ta_xxxzz_xxzzz_1[i] * pc_x[i];

        ta_xxxxzz_xyyyy_0[i] = ta_xxxx_xyyyy_0[i] * fe_0 - ta_xxxx_xyyyy_1[i] * fe_0 + ta_xxxxz_xyyyy_0[i] * pa_z[i] - ta_xxxxz_xyyyy_1[i] * pc_z[i];

        ta_xxxxzz_xyyyz_0[i] = 3.0 * ta_xxzz_xyyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyyz_1[i] * fe_0 + ta_xxxzz_yyyz_0[i] * fe_0 -
                               ta_xxxzz_yyyz_1[i] * fe_0 + ta_xxxzz_xyyyz_0[i] * pa_x[i] - ta_xxxzz_xyyyz_1[i] * pc_x[i];

        ta_xxxxzz_xyyzz_0[i] = 3.0 * ta_xxzz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyzz_1[i] * fe_0 + ta_xxxzz_yyzz_0[i] * fe_0 -
                               ta_xxxzz_yyzz_1[i] * fe_0 + ta_xxxzz_xyyzz_0[i] * pa_x[i] - ta_xxxzz_xyyzz_1[i] * pc_x[i];

        ta_xxxxzz_xyzzz_0[i] = 3.0 * ta_xxzz_xyzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyzzz_1[i] * fe_0 + ta_xxxzz_yzzz_0[i] * fe_0 -
                               ta_xxxzz_yzzz_1[i] * fe_0 + ta_xxxzz_xyzzz_0[i] * pa_x[i] - ta_xxxzz_xyzzz_1[i] * pc_x[i];

        ta_xxxxzz_xzzzz_0[i] = 3.0 * ta_xxzz_xzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_xzzzz_1[i] * fe_0 + ta_xxxzz_zzzz_0[i] * fe_0 -
                               ta_xxxzz_zzzz_1[i] * fe_0 + ta_xxxzz_xzzzz_0[i] * pa_x[i] - ta_xxxzz_xzzzz_1[i] * pc_x[i];

        ta_xxxxzz_yyyyy_0[i] =
            3.0 * ta_xxzz_yyyyy_0[i] * fe_0 - 3.0 * ta_xxzz_yyyyy_1[i] * fe_0 + ta_xxxzz_yyyyy_0[i] * pa_x[i] - ta_xxxzz_yyyyy_1[i] * pc_x[i];

        ta_xxxxzz_yyyyz_0[i] =
            3.0 * ta_xxzz_yyyyz_0[i] * fe_0 - 3.0 * ta_xxzz_yyyyz_1[i] * fe_0 + ta_xxxzz_yyyyz_0[i] * pa_x[i] - ta_xxxzz_yyyyz_1[i] * pc_x[i];

        ta_xxxxzz_yyyzz_0[i] =
            3.0 * ta_xxzz_yyyzz_0[i] * fe_0 - 3.0 * ta_xxzz_yyyzz_1[i] * fe_0 + ta_xxxzz_yyyzz_0[i] * pa_x[i] - ta_xxxzz_yyyzz_1[i] * pc_x[i];

        ta_xxxxzz_yyzzz_0[i] =
            3.0 * ta_xxzz_yyzzz_0[i] * fe_0 - 3.0 * ta_xxzz_yyzzz_1[i] * fe_0 + ta_xxxzz_yyzzz_0[i] * pa_x[i] - ta_xxxzz_yyzzz_1[i] * pc_x[i];

        ta_xxxxzz_yzzzz_0[i] =
            3.0 * ta_xxzz_yzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_yzzzz_1[i] * fe_0 + ta_xxxzz_yzzzz_0[i] * pa_x[i] - ta_xxxzz_yzzzz_1[i] * pc_x[i];

        ta_xxxxzz_zzzzz_0[i] =
            3.0 * ta_xxzz_zzzzz_0[i] * fe_0 - 3.0 * ta_xxzz_zzzzz_1[i] * fe_0 + ta_xxxzz_zzzzz_0[i] * pa_x[i] - ta_xxxzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 126-147 components of targeted buffer : IH

    auto ta_xxxyyy_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 126);

    auto ta_xxxyyy_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 127);

    auto ta_xxxyyy_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 128);

    auto ta_xxxyyy_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 129);

    auto ta_xxxyyy_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 130);

    auto ta_xxxyyy_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 131);

    auto ta_xxxyyy_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 132);

    auto ta_xxxyyy_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 133);

    auto ta_xxxyyy_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 134);

    auto ta_xxxyyy_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 135);

    auto ta_xxxyyy_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 136);

    auto ta_xxxyyy_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 137);

    auto ta_xxxyyy_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 138);

    auto ta_xxxyyy_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 139);

    auto ta_xxxyyy_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 140);

    auto ta_xxxyyy_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 141);

    auto ta_xxxyyy_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 142);

    auto ta_xxxyyy_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 143);

    auto ta_xxxyyy_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 144);

    auto ta_xxxyyy_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 145);

    auto ta_xxxyyy_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 146);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxxy_xxxxx_0,   \
                             ta_xxxy_xxxxx_1,   \
                             ta_xxxy_xxxxz_0,   \
                             ta_xxxy_xxxxz_1,   \
                             ta_xxxy_xxxzz_0,   \
                             ta_xxxy_xxxzz_1,   \
                             ta_xxxy_xxzzz_0,   \
                             ta_xxxy_xxzzz_1,   \
                             ta_xxxy_xzzzz_0,   \
                             ta_xxxy_xzzzz_1,   \
                             ta_xxxyy_xxxxx_0,  \
                             ta_xxxyy_xxxxx_1,  \
                             ta_xxxyy_xxxxz_0,  \
                             ta_xxxyy_xxxxz_1,  \
                             ta_xxxyy_xxxzz_0,  \
                             ta_xxxyy_xxxzz_1,  \
                             ta_xxxyy_xxzzz_0,  \
                             ta_xxxyy_xxzzz_1,  \
                             ta_xxxyy_xzzzz_0,  \
                             ta_xxxyy_xzzzz_1,  \
                             ta_xxxyyy_xxxxx_0, \
                             ta_xxxyyy_xxxxy_0, \
                             ta_xxxyyy_xxxxz_0, \
                             ta_xxxyyy_xxxyy_0, \
                             ta_xxxyyy_xxxyz_0, \
                             ta_xxxyyy_xxxzz_0, \
                             ta_xxxyyy_xxyyy_0, \
                             ta_xxxyyy_xxyyz_0, \
                             ta_xxxyyy_xxyzz_0, \
                             ta_xxxyyy_xxzzz_0, \
                             ta_xxxyyy_xyyyy_0, \
                             ta_xxxyyy_xyyyz_0, \
                             ta_xxxyyy_xyyzz_0, \
                             ta_xxxyyy_xyzzz_0, \
                             ta_xxxyyy_xzzzz_0, \
                             ta_xxxyyy_yyyyy_0, \
                             ta_xxxyyy_yyyyz_0, \
                             ta_xxxyyy_yyyzz_0, \
                             ta_xxxyyy_yyzzz_0, \
                             ta_xxxyyy_yzzzz_0, \
                             ta_xxxyyy_zzzzz_0, \
                             ta_xxyyy_xxxxy_0,  \
                             ta_xxyyy_xxxxy_1,  \
                             ta_xxyyy_xxxy_0,   \
                             ta_xxyyy_xxxy_1,   \
                             ta_xxyyy_xxxyy_0,  \
                             ta_xxyyy_xxxyy_1,  \
                             ta_xxyyy_xxxyz_0,  \
                             ta_xxyyy_xxxyz_1,  \
                             ta_xxyyy_xxyy_0,   \
                             ta_xxyyy_xxyy_1,   \
                             ta_xxyyy_xxyyy_0,  \
                             ta_xxyyy_xxyyy_1,  \
                             ta_xxyyy_xxyyz_0,  \
                             ta_xxyyy_xxyyz_1,  \
                             ta_xxyyy_xxyz_0,   \
                             ta_xxyyy_xxyz_1,   \
                             ta_xxyyy_xxyzz_0,  \
                             ta_xxyyy_xxyzz_1,  \
                             ta_xxyyy_xyyy_0,   \
                             ta_xxyyy_xyyy_1,   \
                             ta_xxyyy_xyyyy_0,  \
                             ta_xxyyy_xyyyy_1,  \
                             ta_xxyyy_xyyyz_0,  \
                             ta_xxyyy_xyyyz_1,  \
                             ta_xxyyy_xyyz_0,   \
                             ta_xxyyy_xyyz_1,   \
                             ta_xxyyy_xyyzz_0,  \
                             ta_xxyyy_xyyzz_1,  \
                             ta_xxyyy_xyzz_0,   \
                             ta_xxyyy_xyzz_1,   \
                             ta_xxyyy_xyzzz_0,  \
                             ta_xxyyy_xyzzz_1,  \
                             ta_xxyyy_yyyy_0,   \
                             ta_xxyyy_yyyy_1,   \
                             ta_xxyyy_yyyyy_0,  \
                             ta_xxyyy_yyyyy_1,  \
                             ta_xxyyy_yyyyz_0,  \
                             ta_xxyyy_yyyyz_1,  \
                             ta_xxyyy_yyyz_0,   \
                             ta_xxyyy_yyyz_1,   \
                             ta_xxyyy_yyyzz_0,  \
                             ta_xxyyy_yyyzz_1,  \
                             ta_xxyyy_yyzz_0,   \
                             ta_xxyyy_yyzz_1,   \
                             ta_xxyyy_yyzzz_0,  \
                             ta_xxyyy_yyzzz_1,  \
                             ta_xxyyy_yzzz_0,   \
                             ta_xxyyy_yzzz_1,   \
                             ta_xxyyy_yzzzz_0,  \
                             ta_xxyyy_yzzzz_1,  \
                             ta_xxyyy_zzzzz_0,  \
                             ta_xxyyy_zzzzz_1,  \
                             ta_xyyy_xxxxy_0,   \
                             ta_xyyy_xxxxy_1,   \
                             ta_xyyy_xxxyy_0,   \
                             ta_xyyy_xxxyy_1,   \
                             ta_xyyy_xxxyz_0,   \
                             ta_xyyy_xxxyz_1,   \
                             ta_xyyy_xxyyy_0,   \
                             ta_xyyy_xxyyy_1,   \
                             ta_xyyy_xxyyz_0,   \
                             ta_xyyy_xxyyz_1,   \
                             ta_xyyy_xxyzz_0,   \
                             ta_xyyy_xxyzz_1,   \
                             ta_xyyy_xyyyy_0,   \
                             ta_xyyy_xyyyy_1,   \
                             ta_xyyy_xyyyz_0,   \
                             ta_xyyy_xyyyz_1,   \
                             ta_xyyy_xyyzz_0,   \
                             ta_xyyy_xyyzz_1,   \
                             ta_xyyy_xyzzz_0,   \
                             ta_xyyy_xyzzz_1,   \
                             ta_xyyy_yyyyy_0,   \
                             ta_xyyy_yyyyy_1,   \
                             ta_xyyy_yyyyz_0,   \
                             ta_xyyy_yyyyz_1,   \
                             ta_xyyy_yyyzz_0,   \
                             ta_xyyy_yyyzz_1,   \
                             ta_xyyy_yyzzz_0,   \
                             ta_xyyy_yyzzz_1,   \
                             ta_xyyy_yzzzz_0,   \
                             ta_xyyy_yzzzz_1,   \
                             ta_xyyy_zzzzz_0,   \
                             ta_xyyy_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyy_xxxxx_0[i] =
            2.0 * ta_xxxy_xxxxx_0[i] * fe_0 - 2.0 * ta_xxxy_xxxxx_1[i] * fe_0 + ta_xxxyy_xxxxx_0[i] * pa_y[i] - ta_xxxyy_xxxxx_1[i] * pc_y[i];

        ta_xxxyyy_xxxxy_0[i] = 2.0 * ta_xyyy_xxxxy_0[i] * fe_0 - 2.0 * ta_xyyy_xxxxy_1[i] * fe_0 + 4.0 * ta_xxyyy_xxxy_0[i] * fe_0 -
                               4.0 * ta_xxyyy_xxxy_1[i] * fe_0 + ta_xxyyy_xxxxy_0[i] * pa_x[i] - ta_xxyyy_xxxxy_1[i] * pc_x[i];

        ta_xxxyyy_xxxxz_0[i] =
            2.0 * ta_xxxy_xxxxz_0[i] * fe_0 - 2.0 * ta_xxxy_xxxxz_1[i] * fe_0 + ta_xxxyy_xxxxz_0[i] * pa_y[i] - ta_xxxyy_xxxxz_1[i] * pc_y[i];

        ta_xxxyyy_xxxyy_0[i] = 2.0 * ta_xyyy_xxxyy_0[i] * fe_0 - 2.0 * ta_xyyy_xxxyy_1[i] * fe_0 + 3.0 * ta_xxyyy_xxyy_0[i] * fe_0 -
                               3.0 * ta_xxyyy_xxyy_1[i] * fe_0 + ta_xxyyy_xxxyy_0[i] * pa_x[i] - ta_xxyyy_xxxyy_1[i] * pc_x[i];

        ta_xxxyyy_xxxyz_0[i] = 2.0 * ta_xyyy_xxxyz_0[i] * fe_0 - 2.0 * ta_xyyy_xxxyz_1[i] * fe_0 + 3.0 * ta_xxyyy_xxyz_0[i] * fe_0 -
                               3.0 * ta_xxyyy_xxyz_1[i] * fe_0 + ta_xxyyy_xxxyz_0[i] * pa_x[i] - ta_xxyyy_xxxyz_1[i] * pc_x[i];

        ta_xxxyyy_xxxzz_0[i] =
            2.0 * ta_xxxy_xxxzz_0[i] * fe_0 - 2.0 * ta_xxxy_xxxzz_1[i] * fe_0 + ta_xxxyy_xxxzz_0[i] * pa_y[i] - ta_xxxyy_xxxzz_1[i] * pc_y[i];

        ta_xxxyyy_xxyyy_0[i] = 2.0 * ta_xyyy_xxyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xxyyy_1[i] * fe_0 + 2.0 * ta_xxyyy_xyyy_0[i] * fe_0 -
                               2.0 * ta_xxyyy_xyyy_1[i] * fe_0 + ta_xxyyy_xxyyy_0[i] * pa_x[i] - ta_xxyyy_xxyyy_1[i] * pc_x[i];

        ta_xxxyyy_xxyyz_0[i] = 2.0 * ta_xyyy_xxyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xxyyz_1[i] * fe_0 + 2.0 * ta_xxyyy_xyyz_0[i] * fe_0 -
                               2.0 * ta_xxyyy_xyyz_1[i] * fe_0 + ta_xxyyy_xxyyz_0[i] * pa_x[i] - ta_xxyyy_xxyyz_1[i] * pc_x[i];

        ta_xxxyyy_xxyzz_0[i] = 2.0 * ta_xyyy_xxyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xxyzz_1[i] * fe_0 + 2.0 * ta_xxyyy_xyzz_0[i] * fe_0 -
                               2.0 * ta_xxyyy_xyzz_1[i] * fe_0 + ta_xxyyy_xxyzz_0[i] * pa_x[i] - ta_xxyyy_xxyzz_1[i] * pc_x[i];

        ta_xxxyyy_xxzzz_0[i] =
            2.0 * ta_xxxy_xxzzz_0[i] * fe_0 - 2.0 * ta_xxxy_xxzzz_1[i] * fe_0 + ta_xxxyy_xxzzz_0[i] * pa_y[i] - ta_xxxyy_xxzzz_1[i] * pc_y[i];

        ta_xxxyyy_xyyyy_0[i] = 2.0 * ta_xyyy_xyyyy_0[i] * fe_0 - 2.0 * ta_xyyy_xyyyy_1[i] * fe_0 + ta_xxyyy_yyyy_0[i] * fe_0 -
                               ta_xxyyy_yyyy_1[i] * fe_0 + ta_xxyyy_xyyyy_0[i] * pa_x[i] - ta_xxyyy_xyyyy_1[i] * pc_x[i];

        ta_xxxyyy_xyyyz_0[i] = 2.0 * ta_xyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyyz_1[i] * fe_0 + ta_xxyyy_yyyz_0[i] * fe_0 -
                               ta_xxyyy_yyyz_1[i] * fe_0 + ta_xxyyy_xyyyz_0[i] * pa_x[i] - ta_xxyyy_xyyyz_1[i] * pc_x[i];

        ta_xxxyyy_xyyzz_0[i] = 2.0 * ta_xyyy_xyyzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyyzz_1[i] * fe_0 + ta_xxyyy_yyzz_0[i] * fe_0 -
                               ta_xxyyy_yyzz_1[i] * fe_0 + ta_xxyyy_xyyzz_0[i] * pa_x[i] - ta_xxyyy_xyyzz_1[i] * pc_x[i];

        ta_xxxyyy_xyzzz_0[i] = 2.0 * ta_xyyy_xyzzz_0[i] * fe_0 - 2.0 * ta_xyyy_xyzzz_1[i] * fe_0 + ta_xxyyy_yzzz_0[i] * fe_0 -
                               ta_xxyyy_yzzz_1[i] * fe_0 + ta_xxyyy_xyzzz_0[i] * pa_x[i] - ta_xxyyy_xyzzz_1[i] * pc_x[i];

        ta_xxxyyy_xzzzz_0[i] =
            2.0 * ta_xxxy_xzzzz_0[i] * fe_0 - 2.0 * ta_xxxy_xzzzz_1[i] * fe_0 + ta_xxxyy_xzzzz_0[i] * pa_y[i] - ta_xxxyy_xzzzz_1[i] * pc_y[i];

        ta_xxxyyy_yyyyy_0[i] =
            2.0 * ta_xyyy_yyyyy_0[i] * fe_0 - 2.0 * ta_xyyy_yyyyy_1[i] * fe_0 + ta_xxyyy_yyyyy_0[i] * pa_x[i] - ta_xxyyy_yyyyy_1[i] * pc_x[i];

        ta_xxxyyy_yyyyz_0[i] =
            2.0 * ta_xyyy_yyyyz_0[i] * fe_0 - 2.0 * ta_xyyy_yyyyz_1[i] * fe_0 + ta_xxyyy_yyyyz_0[i] * pa_x[i] - ta_xxyyy_yyyyz_1[i] * pc_x[i];

        ta_xxxyyy_yyyzz_0[i] =
            2.0 * ta_xyyy_yyyzz_0[i] * fe_0 - 2.0 * ta_xyyy_yyyzz_1[i] * fe_0 + ta_xxyyy_yyyzz_0[i] * pa_x[i] - ta_xxyyy_yyyzz_1[i] * pc_x[i];

        ta_xxxyyy_yyzzz_0[i] =
            2.0 * ta_xyyy_yyzzz_0[i] * fe_0 - 2.0 * ta_xyyy_yyzzz_1[i] * fe_0 + ta_xxyyy_yyzzz_0[i] * pa_x[i] - ta_xxyyy_yyzzz_1[i] * pc_x[i];

        ta_xxxyyy_yzzzz_0[i] =
            2.0 * ta_xyyy_yzzzz_0[i] * fe_0 - 2.0 * ta_xyyy_yzzzz_1[i] * fe_0 + ta_xxyyy_yzzzz_0[i] * pa_x[i] - ta_xxyyy_yzzzz_1[i] * pc_x[i];

        ta_xxxyyy_zzzzz_0[i] =
            2.0 * ta_xyyy_zzzzz_0[i] * fe_0 - 2.0 * ta_xyyy_zzzzz_1[i] * fe_0 + ta_xxyyy_zzzzz_0[i] * pa_x[i] - ta_xxyyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 147-168 components of targeted buffer : IH

    auto ta_xxxyyz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 147);

    auto ta_xxxyyz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 148);

    auto ta_xxxyyz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 149);

    auto ta_xxxyyz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 150);

    auto ta_xxxyyz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 151);

    auto ta_xxxyyz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 152);

    auto ta_xxxyyz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 153);

    auto ta_xxxyyz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 154);

    auto ta_xxxyyz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 155);

    auto ta_xxxyyz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 156);

    auto ta_xxxyyz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 157);

    auto ta_xxxyyz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 158);

    auto ta_xxxyyz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 159);

    auto ta_xxxyyz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 160);

    auto ta_xxxyyz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 161);

    auto ta_xxxyyz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 162);

    auto ta_xxxyyz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 163);

    auto ta_xxxyyz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 164);

    auto ta_xxxyyz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 165);

    auto ta_xxxyyz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 166);

    auto ta_xxxyyz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 167);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta_xxxyy_xxxxx_0,  \
                             ta_xxxyy_xxxxx_1,  \
                             ta_xxxyy_xxxxy_0,  \
                             ta_xxxyy_xxxxy_1,  \
                             ta_xxxyy_xxxy_0,   \
                             ta_xxxyy_xxxy_1,   \
                             ta_xxxyy_xxxyy_0,  \
                             ta_xxxyy_xxxyy_1,  \
                             ta_xxxyy_xxxyz_0,  \
                             ta_xxxyy_xxxyz_1,  \
                             ta_xxxyy_xxyy_0,   \
                             ta_xxxyy_xxyy_1,   \
                             ta_xxxyy_xxyyy_0,  \
                             ta_xxxyy_xxyyy_1,  \
                             ta_xxxyy_xxyyz_0,  \
                             ta_xxxyy_xxyyz_1,  \
                             ta_xxxyy_xxyz_0,   \
                             ta_xxxyy_xxyz_1,   \
                             ta_xxxyy_xxyzz_0,  \
                             ta_xxxyy_xxyzz_1,  \
                             ta_xxxyy_xyyy_0,   \
                             ta_xxxyy_xyyy_1,   \
                             ta_xxxyy_xyyyy_0,  \
                             ta_xxxyy_xyyyy_1,  \
                             ta_xxxyy_xyyyz_0,  \
                             ta_xxxyy_xyyyz_1,  \
                             ta_xxxyy_xyyz_0,   \
                             ta_xxxyy_xyyz_1,   \
                             ta_xxxyy_xyyzz_0,  \
                             ta_xxxyy_xyyzz_1,  \
                             ta_xxxyy_xyzz_0,   \
                             ta_xxxyy_xyzz_1,   \
                             ta_xxxyy_xyzzz_0,  \
                             ta_xxxyy_xyzzz_1,  \
                             ta_xxxyy_yyyyy_0,  \
                             ta_xxxyy_yyyyy_1,  \
                             ta_xxxyyz_xxxxx_0, \
                             ta_xxxyyz_xxxxy_0, \
                             ta_xxxyyz_xxxxz_0, \
                             ta_xxxyyz_xxxyy_0, \
                             ta_xxxyyz_xxxyz_0, \
                             ta_xxxyyz_xxxzz_0, \
                             ta_xxxyyz_xxyyy_0, \
                             ta_xxxyyz_xxyyz_0, \
                             ta_xxxyyz_xxyzz_0, \
                             ta_xxxyyz_xxzzz_0, \
                             ta_xxxyyz_xyyyy_0, \
                             ta_xxxyyz_xyyyz_0, \
                             ta_xxxyyz_xyyzz_0, \
                             ta_xxxyyz_xyzzz_0, \
                             ta_xxxyyz_xzzzz_0, \
                             ta_xxxyyz_yyyyy_0, \
                             ta_xxxyyz_yyyyz_0, \
                             ta_xxxyyz_yyyzz_0, \
                             ta_xxxyyz_yyzzz_0, \
                             ta_xxxyyz_yzzzz_0, \
                             ta_xxxyyz_zzzzz_0, \
                             ta_xxxyz_xxxxz_0,  \
                             ta_xxxyz_xxxxz_1,  \
                             ta_xxxyz_xxxzz_0,  \
                             ta_xxxyz_xxxzz_1,  \
                             ta_xxxyz_xxzzz_0,  \
                             ta_xxxyz_xxzzz_1,  \
                             ta_xxxyz_xzzzz_0,  \
                             ta_xxxyz_xzzzz_1,  \
                             ta_xxxz_xxxxz_0,   \
                             ta_xxxz_xxxxz_1,   \
                             ta_xxxz_xxxzz_0,   \
                             ta_xxxz_xxxzz_1,   \
                             ta_xxxz_xxzzz_0,   \
                             ta_xxxz_xxzzz_1,   \
                             ta_xxxz_xzzzz_0,   \
                             ta_xxxz_xzzzz_1,   \
                             ta_xxyyz_yyyyz_0,  \
                             ta_xxyyz_yyyyz_1,  \
                             ta_xxyyz_yyyzz_0,  \
                             ta_xxyyz_yyyzz_1,  \
                             ta_xxyyz_yyzzz_0,  \
                             ta_xxyyz_yyzzz_1,  \
                             ta_xxyyz_yzzzz_0,  \
                             ta_xxyyz_yzzzz_1,  \
                             ta_xxyyz_zzzzz_0,  \
                             ta_xxyyz_zzzzz_1,  \
                             ta_xyyz_yyyyz_0,   \
                             ta_xyyz_yyyyz_1,   \
                             ta_xyyz_yyyzz_0,   \
                             ta_xyyz_yyyzz_1,   \
                             ta_xyyz_yyzzz_0,   \
                             ta_xyyz_yyzzz_1,   \
                             ta_xyyz_yzzzz_0,   \
                             ta_xyyz_yzzzz_1,   \
                             ta_xyyz_zzzzz_0,   \
                             ta_xyyz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyz_xxxxx_0[i] = ta_xxxyy_xxxxx_0[i] * pa_z[i] - ta_xxxyy_xxxxx_1[i] * pc_z[i];

        ta_xxxyyz_xxxxy_0[i] = ta_xxxyy_xxxxy_0[i] * pa_z[i] - ta_xxxyy_xxxxy_1[i] * pc_z[i];

        ta_xxxyyz_xxxxz_0[i] = ta_xxxz_xxxxz_0[i] * fe_0 - ta_xxxz_xxxxz_1[i] * fe_0 + ta_xxxyz_xxxxz_0[i] * pa_y[i] - ta_xxxyz_xxxxz_1[i] * pc_y[i];

        ta_xxxyyz_xxxyy_0[i] = ta_xxxyy_xxxyy_0[i] * pa_z[i] - ta_xxxyy_xxxyy_1[i] * pc_z[i];

        ta_xxxyyz_xxxyz_0[i] = ta_xxxyy_xxxy_0[i] * fe_0 - ta_xxxyy_xxxy_1[i] * fe_0 + ta_xxxyy_xxxyz_0[i] * pa_z[i] - ta_xxxyy_xxxyz_1[i] * pc_z[i];

        ta_xxxyyz_xxxzz_0[i] = ta_xxxz_xxxzz_0[i] * fe_0 - ta_xxxz_xxxzz_1[i] * fe_0 + ta_xxxyz_xxxzz_0[i] * pa_y[i] - ta_xxxyz_xxxzz_1[i] * pc_y[i];

        ta_xxxyyz_xxyyy_0[i] = ta_xxxyy_xxyyy_0[i] * pa_z[i] - ta_xxxyy_xxyyy_1[i] * pc_z[i];

        ta_xxxyyz_xxyyz_0[i] = ta_xxxyy_xxyy_0[i] * fe_0 - ta_xxxyy_xxyy_1[i] * fe_0 + ta_xxxyy_xxyyz_0[i] * pa_z[i] - ta_xxxyy_xxyyz_1[i] * pc_z[i];

        ta_xxxyyz_xxyzz_0[i] =
            2.0 * ta_xxxyy_xxyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xxyz_1[i] * fe_0 + ta_xxxyy_xxyzz_0[i] * pa_z[i] - ta_xxxyy_xxyzz_1[i] * pc_z[i];

        ta_xxxyyz_xxzzz_0[i] = ta_xxxz_xxzzz_0[i] * fe_0 - ta_xxxz_xxzzz_1[i] * fe_0 + ta_xxxyz_xxzzz_0[i] * pa_y[i] - ta_xxxyz_xxzzz_1[i] * pc_y[i];

        ta_xxxyyz_xyyyy_0[i] = ta_xxxyy_xyyyy_0[i] * pa_z[i] - ta_xxxyy_xyyyy_1[i] * pc_z[i];

        ta_xxxyyz_xyyyz_0[i] = ta_xxxyy_xyyy_0[i] * fe_0 - ta_xxxyy_xyyy_1[i] * fe_0 + ta_xxxyy_xyyyz_0[i] * pa_z[i] - ta_xxxyy_xyyyz_1[i] * pc_z[i];

        ta_xxxyyz_xyyzz_0[i] =
            2.0 * ta_xxxyy_xyyz_0[i] * fe_0 - 2.0 * ta_xxxyy_xyyz_1[i] * fe_0 + ta_xxxyy_xyyzz_0[i] * pa_z[i] - ta_xxxyy_xyyzz_1[i] * pc_z[i];

        ta_xxxyyz_xyzzz_0[i] =
            3.0 * ta_xxxyy_xyzz_0[i] * fe_0 - 3.0 * ta_xxxyy_xyzz_1[i] * fe_0 + ta_xxxyy_xyzzz_0[i] * pa_z[i] - ta_xxxyy_xyzzz_1[i] * pc_z[i];

        ta_xxxyyz_xzzzz_0[i] = ta_xxxz_xzzzz_0[i] * fe_0 - ta_xxxz_xzzzz_1[i] * fe_0 + ta_xxxyz_xzzzz_0[i] * pa_y[i] - ta_xxxyz_xzzzz_1[i] * pc_y[i];

        ta_xxxyyz_yyyyy_0[i] = ta_xxxyy_yyyyy_0[i] * pa_z[i] - ta_xxxyy_yyyyy_1[i] * pc_z[i];

        ta_xxxyyz_yyyyz_0[i] =
            2.0 * ta_xyyz_yyyyz_0[i] * fe_0 - 2.0 * ta_xyyz_yyyyz_1[i] * fe_0 + ta_xxyyz_yyyyz_0[i] * pa_x[i] - ta_xxyyz_yyyyz_1[i] * pc_x[i];

        ta_xxxyyz_yyyzz_0[i] =
            2.0 * ta_xyyz_yyyzz_0[i] * fe_0 - 2.0 * ta_xyyz_yyyzz_1[i] * fe_0 + ta_xxyyz_yyyzz_0[i] * pa_x[i] - ta_xxyyz_yyyzz_1[i] * pc_x[i];

        ta_xxxyyz_yyzzz_0[i] =
            2.0 * ta_xyyz_yyzzz_0[i] * fe_0 - 2.0 * ta_xyyz_yyzzz_1[i] * fe_0 + ta_xxyyz_yyzzz_0[i] * pa_x[i] - ta_xxyyz_yyzzz_1[i] * pc_x[i];

        ta_xxxyyz_yzzzz_0[i] =
            2.0 * ta_xyyz_yzzzz_0[i] * fe_0 - 2.0 * ta_xyyz_yzzzz_1[i] * fe_0 + ta_xxyyz_yzzzz_0[i] * pa_x[i] - ta_xxyyz_yzzzz_1[i] * pc_x[i];

        ta_xxxyyz_zzzzz_0[i] =
            2.0 * ta_xyyz_zzzzz_0[i] * fe_0 - 2.0 * ta_xyyz_zzzzz_1[i] * fe_0 + ta_xxyyz_zzzzz_0[i] * pa_x[i] - ta_xxyyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 168-189 components of targeted buffer : IH

    auto ta_xxxyzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 168);

    auto ta_xxxyzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 169);

    auto ta_xxxyzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 170);

    auto ta_xxxyzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 171);

    auto ta_xxxyzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 172);

    auto ta_xxxyzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 173);

    auto ta_xxxyzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 174);

    auto ta_xxxyzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 175);

    auto ta_xxxyzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 176);

    auto ta_xxxyzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 177);

    auto ta_xxxyzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 178);

    auto ta_xxxyzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 179);

    auto ta_xxxyzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 180);

    auto ta_xxxyzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 181);

    auto ta_xxxyzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 182);

    auto ta_xxxyzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 183);

    auto ta_xxxyzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 184);

    auto ta_xxxyzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 185);

    auto ta_xxxyzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 186);

    auto ta_xxxyzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 187);

    auto ta_xxxyzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 188);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxxyzz_xxxxx_0, \
                             ta_xxxyzz_xxxxy_0, \
                             ta_xxxyzz_xxxxz_0, \
                             ta_xxxyzz_xxxyy_0, \
                             ta_xxxyzz_xxxyz_0, \
                             ta_xxxyzz_xxxzz_0, \
                             ta_xxxyzz_xxyyy_0, \
                             ta_xxxyzz_xxyyz_0, \
                             ta_xxxyzz_xxyzz_0, \
                             ta_xxxyzz_xxzzz_0, \
                             ta_xxxyzz_xyyyy_0, \
                             ta_xxxyzz_xyyyz_0, \
                             ta_xxxyzz_xyyzz_0, \
                             ta_xxxyzz_xyzzz_0, \
                             ta_xxxyzz_xzzzz_0, \
                             ta_xxxyzz_yyyyy_0, \
                             ta_xxxyzz_yyyyz_0, \
                             ta_xxxyzz_yyyzz_0, \
                             ta_xxxyzz_yyzzz_0, \
                             ta_xxxyzz_yzzzz_0, \
                             ta_xxxyzz_zzzzz_0, \
                             ta_xxxzz_xxxx_0,   \
                             ta_xxxzz_xxxx_1,   \
                             ta_xxxzz_xxxxx_0,  \
                             ta_xxxzz_xxxxx_1,  \
                             ta_xxxzz_xxxxy_0,  \
                             ta_xxxzz_xxxxy_1,  \
                             ta_xxxzz_xxxxz_0,  \
                             ta_xxxzz_xxxxz_1,  \
                             ta_xxxzz_xxxy_0,   \
                             ta_xxxzz_xxxy_1,   \
                             ta_xxxzz_xxxyy_0,  \
                             ta_xxxzz_xxxyy_1,  \
                             ta_xxxzz_xxxyz_0,  \
                             ta_xxxzz_xxxyz_1,  \
                             ta_xxxzz_xxxz_0,   \
                             ta_xxxzz_xxxz_1,   \
                             ta_xxxzz_xxxzz_0,  \
                             ta_xxxzz_xxxzz_1,  \
                             ta_xxxzz_xxyy_0,   \
                             ta_xxxzz_xxyy_1,   \
                             ta_xxxzz_xxyyy_0,  \
                             ta_xxxzz_xxyyy_1,  \
                             ta_xxxzz_xxyyz_0,  \
                             ta_xxxzz_xxyyz_1,  \
                             ta_xxxzz_xxyz_0,   \
                             ta_xxxzz_xxyz_1,   \
                             ta_xxxzz_xxyzz_0,  \
                             ta_xxxzz_xxyzz_1,  \
                             ta_xxxzz_xxzz_0,   \
                             ta_xxxzz_xxzz_1,   \
                             ta_xxxzz_xxzzz_0,  \
                             ta_xxxzz_xxzzz_1,  \
                             ta_xxxzz_xyyy_0,   \
                             ta_xxxzz_xyyy_1,   \
                             ta_xxxzz_xyyyy_0,  \
                             ta_xxxzz_xyyyy_1,  \
                             ta_xxxzz_xyyyz_0,  \
                             ta_xxxzz_xyyyz_1,  \
                             ta_xxxzz_xyyz_0,   \
                             ta_xxxzz_xyyz_1,   \
                             ta_xxxzz_xyyzz_0,  \
                             ta_xxxzz_xyyzz_1,  \
                             ta_xxxzz_xyzz_0,   \
                             ta_xxxzz_xyzz_1,   \
                             ta_xxxzz_xyzzz_0,  \
                             ta_xxxzz_xyzzz_1,  \
                             ta_xxxzz_xzzz_0,   \
                             ta_xxxzz_xzzz_1,   \
                             ta_xxxzz_xzzzz_0,  \
                             ta_xxxzz_xzzzz_1,  \
                             ta_xxxzz_zzzzz_0,  \
                             ta_xxxzz_zzzzz_1,  \
                             ta_xxyzz_yyyyy_0,  \
                             ta_xxyzz_yyyyy_1,  \
                             ta_xxyzz_yyyyz_0,  \
                             ta_xxyzz_yyyyz_1,  \
                             ta_xxyzz_yyyzz_0,  \
                             ta_xxyzz_yyyzz_1,  \
                             ta_xxyzz_yyzzz_0,  \
                             ta_xxyzz_yyzzz_1,  \
                             ta_xxyzz_yzzzz_0,  \
                             ta_xxyzz_yzzzz_1,  \
                             ta_xyzz_yyyyy_0,   \
                             ta_xyzz_yyyyy_1,   \
                             ta_xyzz_yyyyz_0,   \
                             ta_xyzz_yyyyz_1,   \
                             ta_xyzz_yyyzz_0,   \
                             ta_xyzz_yyyzz_1,   \
                             ta_xyzz_yyzzz_0,   \
                             ta_xyzz_yyzzz_1,   \
                             ta_xyzz_yzzzz_0,   \
                             ta_xyzz_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyzz_xxxxx_0[i] = ta_xxxzz_xxxxx_0[i] * pa_y[i] - ta_xxxzz_xxxxx_1[i] * pc_y[i];

        ta_xxxyzz_xxxxy_0[i] = ta_xxxzz_xxxx_0[i] * fe_0 - ta_xxxzz_xxxx_1[i] * fe_0 + ta_xxxzz_xxxxy_0[i] * pa_y[i] - ta_xxxzz_xxxxy_1[i] * pc_y[i];

        ta_xxxyzz_xxxxz_0[i] = ta_xxxzz_xxxxz_0[i] * pa_y[i] - ta_xxxzz_xxxxz_1[i] * pc_y[i];

        ta_xxxyzz_xxxyy_0[i] =
            2.0 * ta_xxxzz_xxxy_0[i] * fe_0 - 2.0 * ta_xxxzz_xxxy_1[i] * fe_0 + ta_xxxzz_xxxyy_0[i] * pa_y[i] - ta_xxxzz_xxxyy_1[i] * pc_y[i];

        ta_xxxyzz_xxxyz_0[i] = ta_xxxzz_xxxz_0[i] * fe_0 - ta_xxxzz_xxxz_1[i] * fe_0 + ta_xxxzz_xxxyz_0[i] * pa_y[i] - ta_xxxzz_xxxyz_1[i] * pc_y[i];

        ta_xxxyzz_xxxzz_0[i] = ta_xxxzz_xxxzz_0[i] * pa_y[i] - ta_xxxzz_xxxzz_1[i] * pc_y[i];

        ta_xxxyzz_xxyyy_0[i] =
            3.0 * ta_xxxzz_xxyy_0[i] * fe_0 - 3.0 * ta_xxxzz_xxyy_1[i] * fe_0 + ta_xxxzz_xxyyy_0[i] * pa_y[i] - ta_xxxzz_xxyyy_1[i] * pc_y[i];

        ta_xxxyzz_xxyyz_0[i] =
            2.0 * ta_xxxzz_xxyz_0[i] * fe_0 - 2.0 * ta_xxxzz_xxyz_1[i] * fe_0 + ta_xxxzz_xxyyz_0[i] * pa_y[i] - ta_xxxzz_xxyyz_1[i] * pc_y[i];

        ta_xxxyzz_xxyzz_0[i] = ta_xxxzz_xxzz_0[i] * fe_0 - ta_xxxzz_xxzz_1[i] * fe_0 + ta_xxxzz_xxyzz_0[i] * pa_y[i] - ta_xxxzz_xxyzz_1[i] * pc_y[i];

        ta_xxxyzz_xxzzz_0[i] = ta_xxxzz_xxzzz_0[i] * pa_y[i] - ta_xxxzz_xxzzz_1[i] * pc_y[i];

        ta_xxxyzz_xyyyy_0[i] =
            4.0 * ta_xxxzz_xyyy_0[i] * fe_0 - 4.0 * ta_xxxzz_xyyy_1[i] * fe_0 + ta_xxxzz_xyyyy_0[i] * pa_y[i] - ta_xxxzz_xyyyy_1[i] * pc_y[i];

        ta_xxxyzz_xyyyz_0[i] =
            3.0 * ta_xxxzz_xyyz_0[i] * fe_0 - 3.0 * ta_xxxzz_xyyz_1[i] * fe_0 + ta_xxxzz_xyyyz_0[i] * pa_y[i] - ta_xxxzz_xyyyz_1[i] * pc_y[i];

        ta_xxxyzz_xyyzz_0[i] =
            2.0 * ta_xxxzz_xyzz_0[i] * fe_0 - 2.0 * ta_xxxzz_xyzz_1[i] * fe_0 + ta_xxxzz_xyyzz_0[i] * pa_y[i] - ta_xxxzz_xyyzz_1[i] * pc_y[i];

        ta_xxxyzz_xyzzz_0[i] = ta_xxxzz_xzzz_0[i] * fe_0 - ta_xxxzz_xzzz_1[i] * fe_0 + ta_xxxzz_xyzzz_0[i] * pa_y[i] - ta_xxxzz_xyzzz_1[i] * pc_y[i];

        ta_xxxyzz_xzzzz_0[i] = ta_xxxzz_xzzzz_0[i] * pa_y[i] - ta_xxxzz_xzzzz_1[i] * pc_y[i];

        ta_xxxyzz_yyyyy_0[i] =
            2.0 * ta_xyzz_yyyyy_0[i] * fe_0 - 2.0 * ta_xyzz_yyyyy_1[i] * fe_0 + ta_xxyzz_yyyyy_0[i] * pa_x[i] - ta_xxyzz_yyyyy_1[i] * pc_x[i];

        ta_xxxyzz_yyyyz_0[i] =
            2.0 * ta_xyzz_yyyyz_0[i] * fe_0 - 2.0 * ta_xyzz_yyyyz_1[i] * fe_0 + ta_xxyzz_yyyyz_0[i] * pa_x[i] - ta_xxyzz_yyyyz_1[i] * pc_x[i];

        ta_xxxyzz_yyyzz_0[i] =
            2.0 * ta_xyzz_yyyzz_0[i] * fe_0 - 2.0 * ta_xyzz_yyyzz_1[i] * fe_0 + ta_xxyzz_yyyzz_0[i] * pa_x[i] - ta_xxyzz_yyyzz_1[i] * pc_x[i];

        ta_xxxyzz_yyzzz_0[i] =
            2.0 * ta_xyzz_yyzzz_0[i] * fe_0 - 2.0 * ta_xyzz_yyzzz_1[i] * fe_0 + ta_xxyzz_yyzzz_0[i] * pa_x[i] - ta_xxyzz_yyzzz_1[i] * pc_x[i];

        ta_xxxyzz_yzzzz_0[i] =
            2.0 * ta_xyzz_yzzzz_0[i] * fe_0 - 2.0 * ta_xyzz_yzzzz_1[i] * fe_0 + ta_xxyzz_yzzzz_0[i] * pa_x[i] - ta_xxyzz_yzzzz_1[i] * pc_x[i];

        ta_xxxyzz_zzzzz_0[i] = ta_xxxzz_zzzzz_0[i] * pa_y[i] - ta_xxxzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 189-210 components of targeted buffer : IH

    auto ta_xxxzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 189);

    auto ta_xxxzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 190);

    auto ta_xxxzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 191);

    auto ta_xxxzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 192);

    auto ta_xxxzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 193);

    auto ta_xxxzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 194);

    auto ta_xxxzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 195);

    auto ta_xxxzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 196);

    auto ta_xxxzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 197);

    auto ta_xxxzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 198);

    auto ta_xxxzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 199);

    auto ta_xxxzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 200);

    auto ta_xxxzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 201);

    auto ta_xxxzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 202);

    auto ta_xxxzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 203);

    auto ta_xxxzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 204);

    auto ta_xxxzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 205);

    auto ta_xxxzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 206);

    auto ta_xxxzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 207);

    auto ta_xxxzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 208);

    auto ta_xxxzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 209);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxxz_xxxxx_0,   \
                             ta_xxxz_xxxxx_1,   \
                             ta_xxxz_xxxxy_0,   \
                             ta_xxxz_xxxxy_1,   \
                             ta_xxxz_xxxyy_0,   \
                             ta_xxxz_xxxyy_1,   \
                             ta_xxxz_xxyyy_0,   \
                             ta_xxxz_xxyyy_1,   \
                             ta_xxxz_xyyyy_0,   \
                             ta_xxxz_xyyyy_1,   \
                             ta_xxxzz_xxxxx_0,  \
                             ta_xxxzz_xxxxx_1,  \
                             ta_xxxzz_xxxxy_0,  \
                             ta_xxxzz_xxxxy_1,  \
                             ta_xxxzz_xxxyy_0,  \
                             ta_xxxzz_xxxyy_1,  \
                             ta_xxxzz_xxyyy_0,  \
                             ta_xxxzz_xxyyy_1,  \
                             ta_xxxzz_xyyyy_0,  \
                             ta_xxxzz_xyyyy_1,  \
                             ta_xxxzzz_xxxxx_0, \
                             ta_xxxzzz_xxxxy_0, \
                             ta_xxxzzz_xxxxz_0, \
                             ta_xxxzzz_xxxyy_0, \
                             ta_xxxzzz_xxxyz_0, \
                             ta_xxxzzz_xxxzz_0, \
                             ta_xxxzzz_xxyyy_0, \
                             ta_xxxzzz_xxyyz_0, \
                             ta_xxxzzz_xxyzz_0, \
                             ta_xxxzzz_xxzzz_0, \
                             ta_xxxzzz_xyyyy_0, \
                             ta_xxxzzz_xyyyz_0, \
                             ta_xxxzzz_xyyzz_0, \
                             ta_xxxzzz_xyzzz_0, \
                             ta_xxxzzz_xzzzz_0, \
                             ta_xxxzzz_yyyyy_0, \
                             ta_xxxzzz_yyyyz_0, \
                             ta_xxxzzz_yyyzz_0, \
                             ta_xxxzzz_yyzzz_0, \
                             ta_xxxzzz_yzzzz_0, \
                             ta_xxxzzz_zzzzz_0, \
                             ta_xxzzz_xxxxz_0,  \
                             ta_xxzzz_xxxxz_1,  \
                             ta_xxzzz_xxxyz_0,  \
                             ta_xxzzz_xxxyz_1,  \
                             ta_xxzzz_xxxz_0,   \
                             ta_xxzzz_xxxz_1,   \
                             ta_xxzzz_xxxzz_0,  \
                             ta_xxzzz_xxxzz_1,  \
                             ta_xxzzz_xxyyz_0,  \
                             ta_xxzzz_xxyyz_1,  \
                             ta_xxzzz_xxyz_0,   \
                             ta_xxzzz_xxyz_1,   \
                             ta_xxzzz_xxyzz_0,  \
                             ta_xxzzz_xxyzz_1,  \
                             ta_xxzzz_xxzz_0,   \
                             ta_xxzzz_xxzz_1,   \
                             ta_xxzzz_xxzzz_0,  \
                             ta_xxzzz_xxzzz_1,  \
                             ta_xxzzz_xyyyz_0,  \
                             ta_xxzzz_xyyyz_1,  \
                             ta_xxzzz_xyyz_0,   \
                             ta_xxzzz_xyyz_1,   \
                             ta_xxzzz_xyyzz_0,  \
                             ta_xxzzz_xyyzz_1,  \
                             ta_xxzzz_xyzz_0,   \
                             ta_xxzzz_xyzz_1,   \
                             ta_xxzzz_xyzzz_0,  \
                             ta_xxzzz_xyzzz_1,  \
                             ta_xxzzz_xzzz_0,   \
                             ta_xxzzz_xzzz_1,   \
                             ta_xxzzz_xzzzz_0,  \
                             ta_xxzzz_xzzzz_1,  \
                             ta_xxzzz_yyyyy_0,  \
                             ta_xxzzz_yyyyy_1,  \
                             ta_xxzzz_yyyyz_0,  \
                             ta_xxzzz_yyyyz_1,  \
                             ta_xxzzz_yyyz_0,   \
                             ta_xxzzz_yyyz_1,   \
                             ta_xxzzz_yyyzz_0,  \
                             ta_xxzzz_yyyzz_1,  \
                             ta_xxzzz_yyzz_0,   \
                             ta_xxzzz_yyzz_1,   \
                             ta_xxzzz_yyzzz_0,  \
                             ta_xxzzz_yyzzz_1,  \
                             ta_xxzzz_yzzz_0,   \
                             ta_xxzzz_yzzz_1,   \
                             ta_xxzzz_yzzzz_0,  \
                             ta_xxzzz_yzzzz_1,  \
                             ta_xxzzz_zzzz_0,   \
                             ta_xxzzz_zzzz_1,   \
                             ta_xxzzz_zzzzz_0,  \
                             ta_xxzzz_zzzzz_1,  \
                             ta_xzzz_xxxxz_0,   \
                             ta_xzzz_xxxxz_1,   \
                             ta_xzzz_xxxyz_0,   \
                             ta_xzzz_xxxyz_1,   \
                             ta_xzzz_xxxzz_0,   \
                             ta_xzzz_xxxzz_1,   \
                             ta_xzzz_xxyyz_0,   \
                             ta_xzzz_xxyyz_1,   \
                             ta_xzzz_xxyzz_0,   \
                             ta_xzzz_xxyzz_1,   \
                             ta_xzzz_xxzzz_0,   \
                             ta_xzzz_xxzzz_1,   \
                             ta_xzzz_xyyyz_0,   \
                             ta_xzzz_xyyyz_1,   \
                             ta_xzzz_xyyzz_0,   \
                             ta_xzzz_xyyzz_1,   \
                             ta_xzzz_xyzzz_0,   \
                             ta_xzzz_xyzzz_1,   \
                             ta_xzzz_xzzzz_0,   \
                             ta_xzzz_xzzzz_1,   \
                             ta_xzzz_yyyyy_0,   \
                             ta_xzzz_yyyyy_1,   \
                             ta_xzzz_yyyyz_0,   \
                             ta_xzzz_yyyyz_1,   \
                             ta_xzzz_yyyzz_0,   \
                             ta_xzzz_yyyzz_1,   \
                             ta_xzzz_yyzzz_0,   \
                             ta_xzzz_yyzzz_1,   \
                             ta_xzzz_yzzzz_0,   \
                             ta_xzzz_yzzzz_1,   \
                             ta_xzzz_zzzzz_0,   \
                             ta_xzzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzzz_xxxxx_0[i] =
            2.0 * ta_xxxz_xxxxx_0[i] * fe_0 - 2.0 * ta_xxxz_xxxxx_1[i] * fe_0 + ta_xxxzz_xxxxx_0[i] * pa_z[i] - ta_xxxzz_xxxxx_1[i] * pc_z[i];

        ta_xxxzzz_xxxxy_0[i] =
            2.0 * ta_xxxz_xxxxy_0[i] * fe_0 - 2.0 * ta_xxxz_xxxxy_1[i] * fe_0 + ta_xxxzz_xxxxy_0[i] * pa_z[i] - ta_xxxzz_xxxxy_1[i] * pc_z[i];

        ta_xxxzzz_xxxxz_0[i] = 2.0 * ta_xzzz_xxxxz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxxz_1[i] * fe_0 + 4.0 * ta_xxzzz_xxxz_0[i] * fe_0 -
                               4.0 * ta_xxzzz_xxxz_1[i] * fe_0 + ta_xxzzz_xxxxz_0[i] * pa_x[i] - ta_xxzzz_xxxxz_1[i] * pc_x[i];

        ta_xxxzzz_xxxyy_0[i] =
            2.0 * ta_xxxz_xxxyy_0[i] * fe_0 - 2.0 * ta_xxxz_xxxyy_1[i] * fe_0 + ta_xxxzz_xxxyy_0[i] * pa_z[i] - ta_xxxzz_xxxyy_1[i] * pc_z[i];

        ta_xxxzzz_xxxyz_0[i] = 2.0 * ta_xzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxyz_1[i] * fe_0 + 3.0 * ta_xxzzz_xxyz_0[i] * fe_0 -
                               3.0 * ta_xxzzz_xxyz_1[i] * fe_0 + ta_xxzzz_xxxyz_0[i] * pa_x[i] - ta_xxzzz_xxxyz_1[i] * pc_x[i];

        ta_xxxzzz_xxxzz_0[i] = 2.0 * ta_xzzz_xxxzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxxzz_1[i] * fe_0 + 3.0 * ta_xxzzz_xxzz_0[i] * fe_0 -
                               3.0 * ta_xxzzz_xxzz_1[i] * fe_0 + ta_xxzzz_xxxzz_0[i] * pa_x[i] - ta_xxzzz_xxxzz_1[i] * pc_x[i];

        ta_xxxzzz_xxyyy_0[i] =
            2.0 * ta_xxxz_xxyyy_0[i] * fe_0 - 2.0 * ta_xxxz_xxyyy_1[i] * fe_0 + ta_xxxzz_xxyyy_0[i] * pa_z[i] - ta_xxxzz_xxyyy_1[i] * pc_z[i];

        ta_xxxzzz_xxyyz_0[i] = 2.0 * ta_xzzz_xxyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xxyyz_1[i] * fe_0 + 2.0 * ta_xxzzz_xyyz_0[i] * fe_0 -
                               2.0 * ta_xxzzz_xyyz_1[i] * fe_0 + ta_xxzzz_xxyyz_0[i] * pa_x[i] - ta_xxzzz_xxyyz_1[i] * pc_x[i];

        ta_xxxzzz_xxyzz_0[i] = 2.0 * ta_xzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxyzz_1[i] * fe_0 + 2.0 * ta_xxzzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_xxzzz_xyzz_1[i] * fe_0 + ta_xxzzz_xxyzz_0[i] * pa_x[i] - ta_xxzzz_xxyzz_1[i] * pc_x[i];

        ta_xxxzzz_xxzzz_0[i] = 2.0 * ta_xzzz_xxzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xxzzz_1[i] * fe_0 + 2.0 * ta_xxzzz_xzzz_0[i] * fe_0 -
                               2.0 * ta_xxzzz_xzzz_1[i] * fe_0 + ta_xxzzz_xxzzz_0[i] * pa_x[i] - ta_xxzzz_xxzzz_1[i] * pc_x[i];

        ta_xxxzzz_xyyyy_0[i] =
            2.0 * ta_xxxz_xyyyy_0[i] * fe_0 - 2.0 * ta_xxxz_xyyyy_1[i] * fe_0 + ta_xxxzz_xyyyy_0[i] * pa_z[i] - ta_xxxzz_xyyyy_1[i] * pc_z[i];

        ta_xxxzzz_xyyyz_0[i] = 2.0 * ta_xzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyyz_1[i] * fe_0 + ta_xxzzz_yyyz_0[i] * fe_0 -
                               ta_xxzzz_yyyz_1[i] * fe_0 + ta_xxzzz_xyyyz_0[i] * pa_x[i] - ta_xxzzz_xyyyz_1[i] * pc_x[i];

        ta_xxxzzz_xyyzz_0[i] = 2.0 * ta_xzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyyzz_1[i] * fe_0 + ta_xxzzz_yyzz_0[i] * fe_0 -
                               ta_xxzzz_yyzz_1[i] * fe_0 + ta_xxzzz_xyyzz_0[i] * pa_x[i] - ta_xxzzz_xyyzz_1[i] * pc_x[i];

        ta_xxxzzz_xyzzz_0[i] = 2.0 * ta_xzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xyzzz_1[i] * fe_0 + ta_xxzzz_yzzz_0[i] * fe_0 -
                               ta_xxzzz_yzzz_1[i] * fe_0 + ta_xxzzz_xyzzz_0[i] * pa_x[i] - ta_xxzzz_xyzzz_1[i] * pc_x[i];

        ta_xxxzzz_xzzzz_0[i] = 2.0 * ta_xzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_xzzzz_1[i] * fe_0 + ta_xxzzz_zzzz_0[i] * fe_0 -
                               ta_xxzzz_zzzz_1[i] * fe_0 + ta_xxzzz_xzzzz_0[i] * pa_x[i] - ta_xxzzz_xzzzz_1[i] * pc_x[i];

        ta_xxxzzz_yyyyy_0[i] =
            2.0 * ta_xzzz_yyyyy_0[i] * fe_0 - 2.0 * ta_xzzz_yyyyy_1[i] * fe_0 + ta_xxzzz_yyyyy_0[i] * pa_x[i] - ta_xxzzz_yyyyy_1[i] * pc_x[i];

        ta_xxxzzz_yyyyz_0[i] =
            2.0 * ta_xzzz_yyyyz_0[i] * fe_0 - 2.0 * ta_xzzz_yyyyz_1[i] * fe_0 + ta_xxzzz_yyyyz_0[i] * pa_x[i] - ta_xxzzz_yyyyz_1[i] * pc_x[i];

        ta_xxxzzz_yyyzz_0[i] =
            2.0 * ta_xzzz_yyyzz_0[i] * fe_0 - 2.0 * ta_xzzz_yyyzz_1[i] * fe_0 + ta_xxzzz_yyyzz_0[i] * pa_x[i] - ta_xxzzz_yyyzz_1[i] * pc_x[i];

        ta_xxxzzz_yyzzz_0[i] =
            2.0 * ta_xzzz_yyzzz_0[i] * fe_0 - 2.0 * ta_xzzz_yyzzz_1[i] * fe_0 + ta_xxzzz_yyzzz_0[i] * pa_x[i] - ta_xxzzz_yyzzz_1[i] * pc_x[i];

        ta_xxxzzz_yzzzz_0[i] =
            2.0 * ta_xzzz_yzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_yzzzz_1[i] * fe_0 + ta_xxzzz_yzzzz_0[i] * pa_x[i] - ta_xxzzz_yzzzz_1[i] * pc_x[i];

        ta_xxxzzz_zzzzz_0[i] =
            2.0 * ta_xzzz_zzzzz_0[i] * fe_0 - 2.0 * ta_xzzz_zzzzz_1[i] * fe_0 + ta_xxzzz_zzzzz_0[i] * pa_x[i] - ta_xxzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 210-231 components of targeted buffer : IH

    auto ta_xxyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 210);

    auto ta_xxyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 211);

    auto ta_xxyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 212);

    auto ta_xxyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 213);

    auto ta_xxyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 214);

    auto ta_xxyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 215);

    auto ta_xxyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 216);

    auto ta_xxyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 217);

    auto ta_xxyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 218);

    auto ta_xxyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 219);

    auto ta_xxyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 220);

    auto ta_xxyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 221);

    auto ta_xxyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 222);

    auto ta_xxyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 223);

    auto ta_xxyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 224);

    auto ta_xxyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 225);

    auto ta_xxyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 226);

    auto ta_xxyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 227);

    auto ta_xxyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 228);

    auto ta_xxyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 229);

    auto ta_xxyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 230);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxyy_xxxxx_0,   \
                             ta_xxyy_xxxxx_1,   \
                             ta_xxyy_xxxxz_0,   \
                             ta_xxyy_xxxxz_1,   \
                             ta_xxyy_xxxzz_0,   \
                             ta_xxyy_xxxzz_1,   \
                             ta_xxyy_xxzzz_0,   \
                             ta_xxyy_xxzzz_1,   \
                             ta_xxyy_xzzzz_0,   \
                             ta_xxyy_xzzzz_1,   \
                             ta_xxyyy_xxxxx_0,  \
                             ta_xxyyy_xxxxx_1,  \
                             ta_xxyyy_xxxxz_0,  \
                             ta_xxyyy_xxxxz_1,  \
                             ta_xxyyy_xxxzz_0,  \
                             ta_xxyyy_xxxzz_1,  \
                             ta_xxyyy_xxzzz_0,  \
                             ta_xxyyy_xxzzz_1,  \
                             ta_xxyyy_xzzzz_0,  \
                             ta_xxyyy_xzzzz_1,  \
                             ta_xxyyyy_xxxxx_0, \
                             ta_xxyyyy_xxxxy_0, \
                             ta_xxyyyy_xxxxz_0, \
                             ta_xxyyyy_xxxyy_0, \
                             ta_xxyyyy_xxxyz_0, \
                             ta_xxyyyy_xxxzz_0, \
                             ta_xxyyyy_xxyyy_0, \
                             ta_xxyyyy_xxyyz_0, \
                             ta_xxyyyy_xxyzz_0, \
                             ta_xxyyyy_xxzzz_0, \
                             ta_xxyyyy_xyyyy_0, \
                             ta_xxyyyy_xyyyz_0, \
                             ta_xxyyyy_xyyzz_0, \
                             ta_xxyyyy_xyzzz_0, \
                             ta_xxyyyy_xzzzz_0, \
                             ta_xxyyyy_yyyyy_0, \
                             ta_xxyyyy_yyyyz_0, \
                             ta_xxyyyy_yyyzz_0, \
                             ta_xxyyyy_yyzzz_0, \
                             ta_xxyyyy_yzzzz_0, \
                             ta_xxyyyy_zzzzz_0, \
                             ta_xyyyy_xxxxy_0,  \
                             ta_xyyyy_xxxxy_1,  \
                             ta_xyyyy_xxxy_0,   \
                             ta_xyyyy_xxxy_1,   \
                             ta_xyyyy_xxxyy_0,  \
                             ta_xyyyy_xxxyy_1,  \
                             ta_xyyyy_xxxyz_0,  \
                             ta_xyyyy_xxxyz_1,  \
                             ta_xyyyy_xxyy_0,   \
                             ta_xyyyy_xxyy_1,   \
                             ta_xyyyy_xxyyy_0,  \
                             ta_xyyyy_xxyyy_1,  \
                             ta_xyyyy_xxyyz_0,  \
                             ta_xyyyy_xxyyz_1,  \
                             ta_xyyyy_xxyz_0,   \
                             ta_xyyyy_xxyz_1,   \
                             ta_xyyyy_xxyzz_0,  \
                             ta_xyyyy_xxyzz_1,  \
                             ta_xyyyy_xyyy_0,   \
                             ta_xyyyy_xyyy_1,   \
                             ta_xyyyy_xyyyy_0,  \
                             ta_xyyyy_xyyyy_1,  \
                             ta_xyyyy_xyyyz_0,  \
                             ta_xyyyy_xyyyz_1,  \
                             ta_xyyyy_xyyz_0,   \
                             ta_xyyyy_xyyz_1,   \
                             ta_xyyyy_xyyzz_0,  \
                             ta_xyyyy_xyyzz_1,  \
                             ta_xyyyy_xyzz_0,   \
                             ta_xyyyy_xyzz_1,   \
                             ta_xyyyy_xyzzz_0,  \
                             ta_xyyyy_xyzzz_1,  \
                             ta_xyyyy_yyyy_0,   \
                             ta_xyyyy_yyyy_1,   \
                             ta_xyyyy_yyyyy_0,  \
                             ta_xyyyy_yyyyy_1,  \
                             ta_xyyyy_yyyyz_0,  \
                             ta_xyyyy_yyyyz_1,  \
                             ta_xyyyy_yyyz_0,   \
                             ta_xyyyy_yyyz_1,   \
                             ta_xyyyy_yyyzz_0,  \
                             ta_xyyyy_yyyzz_1,  \
                             ta_xyyyy_yyzz_0,   \
                             ta_xyyyy_yyzz_1,   \
                             ta_xyyyy_yyzzz_0,  \
                             ta_xyyyy_yyzzz_1,  \
                             ta_xyyyy_yzzz_0,   \
                             ta_xyyyy_yzzz_1,   \
                             ta_xyyyy_yzzzz_0,  \
                             ta_xyyyy_yzzzz_1,  \
                             ta_xyyyy_zzzzz_0,  \
                             ta_xyyyy_zzzzz_1,  \
                             ta_yyyy_xxxxy_0,   \
                             ta_yyyy_xxxxy_1,   \
                             ta_yyyy_xxxyy_0,   \
                             ta_yyyy_xxxyy_1,   \
                             ta_yyyy_xxxyz_0,   \
                             ta_yyyy_xxxyz_1,   \
                             ta_yyyy_xxyyy_0,   \
                             ta_yyyy_xxyyy_1,   \
                             ta_yyyy_xxyyz_0,   \
                             ta_yyyy_xxyyz_1,   \
                             ta_yyyy_xxyzz_0,   \
                             ta_yyyy_xxyzz_1,   \
                             ta_yyyy_xyyyy_0,   \
                             ta_yyyy_xyyyy_1,   \
                             ta_yyyy_xyyyz_0,   \
                             ta_yyyy_xyyyz_1,   \
                             ta_yyyy_xyyzz_0,   \
                             ta_yyyy_xyyzz_1,   \
                             ta_yyyy_xyzzz_0,   \
                             ta_yyyy_xyzzz_1,   \
                             ta_yyyy_yyyyy_0,   \
                             ta_yyyy_yyyyy_1,   \
                             ta_yyyy_yyyyz_0,   \
                             ta_yyyy_yyyyz_1,   \
                             ta_yyyy_yyyzz_0,   \
                             ta_yyyy_yyyzz_1,   \
                             ta_yyyy_yyzzz_0,   \
                             ta_yyyy_yyzzz_1,   \
                             ta_yyyy_yzzzz_0,   \
                             ta_yyyy_yzzzz_1,   \
                             ta_yyyy_zzzzz_0,   \
                             ta_yyyy_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyy_xxxxx_0[i] =
            3.0 * ta_xxyy_xxxxx_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxx_1[i] * fe_0 + ta_xxyyy_xxxxx_0[i] * pa_y[i] - ta_xxyyy_xxxxx_1[i] * pc_y[i];

        ta_xxyyyy_xxxxy_0[i] = ta_yyyy_xxxxy_0[i] * fe_0 - ta_yyyy_xxxxy_1[i] * fe_0 + 4.0 * ta_xyyyy_xxxy_0[i] * fe_0 -
                               4.0 * ta_xyyyy_xxxy_1[i] * fe_0 + ta_xyyyy_xxxxy_0[i] * pa_x[i] - ta_xyyyy_xxxxy_1[i] * pc_x[i];

        ta_xxyyyy_xxxxz_0[i] =
            3.0 * ta_xxyy_xxxxz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxxz_1[i] * fe_0 + ta_xxyyy_xxxxz_0[i] * pa_y[i] - ta_xxyyy_xxxxz_1[i] * pc_y[i];

        ta_xxyyyy_xxxyy_0[i] = ta_yyyy_xxxyy_0[i] * fe_0 - ta_yyyy_xxxyy_1[i] * fe_0 + 3.0 * ta_xyyyy_xxyy_0[i] * fe_0 -
                               3.0 * ta_xyyyy_xxyy_1[i] * fe_0 + ta_xyyyy_xxxyy_0[i] * pa_x[i] - ta_xyyyy_xxxyy_1[i] * pc_x[i];

        ta_xxyyyy_xxxyz_0[i] = ta_yyyy_xxxyz_0[i] * fe_0 - ta_yyyy_xxxyz_1[i] * fe_0 + 3.0 * ta_xyyyy_xxyz_0[i] * fe_0 -
                               3.0 * ta_xyyyy_xxyz_1[i] * fe_0 + ta_xyyyy_xxxyz_0[i] * pa_x[i] - ta_xyyyy_xxxyz_1[i] * pc_x[i];

        ta_xxyyyy_xxxzz_0[i] =
            3.0 * ta_xxyy_xxxzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxxzz_1[i] * fe_0 + ta_xxyyy_xxxzz_0[i] * pa_y[i] - ta_xxyyy_xxxzz_1[i] * pc_y[i];

        ta_xxyyyy_xxyyy_0[i] = ta_yyyy_xxyyy_0[i] * fe_0 - ta_yyyy_xxyyy_1[i] * fe_0 + 2.0 * ta_xyyyy_xyyy_0[i] * fe_0 -
                               2.0 * ta_xyyyy_xyyy_1[i] * fe_0 + ta_xyyyy_xxyyy_0[i] * pa_x[i] - ta_xyyyy_xxyyy_1[i] * pc_x[i];

        ta_xxyyyy_xxyyz_0[i] = ta_yyyy_xxyyz_0[i] * fe_0 - ta_yyyy_xxyyz_1[i] * fe_0 + 2.0 * ta_xyyyy_xyyz_0[i] * fe_0 -
                               2.0 * ta_xyyyy_xyyz_1[i] * fe_0 + ta_xyyyy_xxyyz_0[i] * pa_x[i] - ta_xyyyy_xxyyz_1[i] * pc_x[i];

        ta_xxyyyy_xxyzz_0[i] = ta_yyyy_xxyzz_0[i] * fe_0 - ta_yyyy_xxyzz_1[i] * fe_0 + 2.0 * ta_xyyyy_xyzz_0[i] * fe_0 -
                               2.0 * ta_xyyyy_xyzz_1[i] * fe_0 + ta_xyyyy_xxyzz_0[i] * pa_x[i] - ta_xyyyy_xxyzz_1[i] * pc_x[i];

        ta_xxyyyy_xxzzz_0[i] =
            3.0 * ta_xxyy_xxzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxzzz_1[i] * fe_0 + ta_xxyyy_xxzzz_0[i] * pa_y[i] - ta_xxyyy_xxzzz_1[i] * pc_y[i];

        ta_xxyyyy_xyyyy_0[i] = ta_yyyy_xyyyy_0[i] * fe_0 - ta_yyyy_xyyyy_1[i] * fe_0 + ta_xyyyy_yyyy_0[i] * fe_0 - ta_xyyyy_yyyy_1[i] * fe_0 +
                               ta_xyyyy_xyyyy_0[i] * pa_x[i] - ta_xyyyy_xyyyy_1[i] * pc_x[i];

        ta_xxyyyy_xyyyz_0[i] = ta_yyyy_xyyyz_0[i] * fe_0 - ta_yyyy_xyyyz_1[i] * fe_0 + ta_xyyyy_yyyz_0[i] * fe_0 - ta_xyyyy_yyyz_1[i] * fe_0 +
                               ta_xyyyy_xyyyz_0[i] * pa_x[i] - ta_xyyyy_xyyyz_1[i] * pc_x[i];

        ta_xxyyyy_xyyzz_0[i] = ta_yyyy_xyyzz_0[i] * fe_0 - ta_yyyy_xyyzz_1[i] * fe_0 + ta_xyyyy_yyzz_0[i] * fe_0 - ta_xyyyy_yyzz_1[i] * fe_0 +
                               ta_xyyyy_xyyzz_0[i] * pa_x[i] - ta_xyyyy_xyyzz_1[i] * pc_x[i];

        ta_xxyyyy_xyzzz_0[i] = ta_yyyy_xyzzz_0[i] * fe_0 - ta_yyyy_xyzzz_1[i] * fe_0 + ta_xyyyy_yzzz_0[i] * fe_0 - ta_xyyyy_yzzz_1[i] * fe_0 +
                               ta_xyyyy_xyzzz_0[i] * pa_x[i] - ta_xyyyy_xyzzz_1[i] * pc_x[i];

        ta_xxyyyy_xzzzz_0[i] =
            3.0 * ta_xxyy_xzzzz_0[i] * fe_0 - 3.0 * ta_xxyy_xzzzz_1[i] * fe_0 + ta_xxyyy_xzzzz_0[i] * pa_y[i] - ta_xxyyy_xzzzz_1[i] * pc_y[i];

        ta_xxyyyy_yyyyy_0[i] = ta_yyyy_yyyyy_0[i] * fe_0 - ta_yyyy_yyyyy_1[i] * fe_0 + ta_xyyyy_yyyyy_0[i] * pa_x[i] - ta_xyyyy_yyyyy_1[i] * pc_x[i];

        ta_xxyyyy_yyyyz_0[i] = ta_yyyy_yyyyz_0[i] * fe_0 - ta_yyyy_yyyyz_1[i] * fe_0 + ta_xyyyy_yyyyz_0[i] * pa_x[i] - ta_xyyyy_yyyyz_1[i] * pc_x[i];

        ta_xxyyyy_yyyzz_0[i] = ta_yyyy_yyyzz_0[i] * fe_0 - ta_yyyy_yyyzz_1[i] * fe_0 + ta_xyyyy_yyyzz_0[i] * pa_x[i] - ta_xyyyy_yyyzz_1[i] * pc_x[i];

        ta_xxyyyy_yyzzz_0[i] = ta_yyyy_yyzzz_0[i] * fe_0 - ta_yyyy_yyzzz_1[i] * fe_0 + ta_xyyyy_yyzzz_0[i] * pa_x[i] - ta_xyyyy_yyzzz_1[i] * pc_x[i];

        ta_xxyyyy_yzzzz_0[i] = ta_yyyy_yzzzz_0[i] * fe_0 - ta_yyyy_yzzzz_1[i] * fe_0 + ta_xyyyy_yzzzz_0[i] * pa_x[i] - ta_xyyyy_yzzzz_1[i] * pc_x[i];

        ta_xxyyyy_zzzzz_0[i] = ta_yyyy_zzzzz_0[i] * fe_0 - ta_yyyy_zzzzz_1[i] * fe_0 + ta_xyyyy_zzzzz_0[i] * pa_x[i] - ta_xyyyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 231-252 components of targeted buffer : IH

    auto ta_xxyyyz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 231);

    auto ta_xxyyyz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 232);

    auto ta_xxyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 233);

    auto ta_xxyyyz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 234);

    auto ta_xxyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 235);

    auto ta_xxyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 236);

    auto ta_xxyyyz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 237);

    auto ta_xxyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 238);

    auto ta_xxyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 239);

    auto ta_xxyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 240);

    auto ta_xxyyyz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 241);

    auto ta_xxyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 242);

    auto ta_xxyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 243);

    auto ta_xxyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 244);

    auto ta_xxyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 245);

    auto ta_xxyyyz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 246);

    auto ta_xxyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 247);

    auto ta_xxyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 248);

    auto ta_xxyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 249);

    auto ta_xxyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 250);

    auto ta_xxyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 251);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta_xxyyy_xxxxx_0,  \
                             ta_xxyyy_xxxxx_1,  \
                             ta_xxyyy_xxxxy_0,  \
                             ta_xxyyy_xxxxy_1,  \
                             ta_xxyyy_xxxy_0,   \
                             ta_xxyyy_xxxy_1,   \
                             ta_xxyyy_xxxyy_0,  \
                             ta_xxyyy_xxxyy_1,  \
                             ta_xxyyy_xxxyz_0,  \
                             ta_xxyyy_xxxyz_1,  \
                             ta_xxyyy_xxyy_0,   \
                             ta_xxyyy_xxyy_1,   \
                             ta_xxyyy_xxyyy_0,  \
                             ta_xxyyy_xxyyy_1,  \
                             ta_xxyyy_xxyyz_0,  \
                             ta_xxyyy_xxyyz_1,  \
                             ta_xxyyy_xxyz_0,   \
                             ta_xxyyy_xxyz_1,   \
                             ta_xxyyy_xxyzz_0,  \
                             ta_xxyyy_xxyzz_1,  \
                             ta_xxyyy_xyyy_0,   \
                             ta_xxyyy_xyyy_1,   \
                             ta_xxyyy_xyyyy_0,  \
                             ta_xxyyy_xyyyy_1,  \
                             ta_xxyyy_xyyyz_0,  \
                             ta_xxyyy_xyyyz_1,  \
                             ta_xxyyy_xyyz_0,   \
                             ta_xxyyy_xyyz_1,   \
                             ta_xxyyy_xyyzz_0,  \
                             ta_xxyyy_xyyzz_1,  \
                             ta_xxyyy_xyzz_0,   \
                             ta_xxyyy_xyzz_1,   \
                             ta_xxyyy_xyzzz_0,  \
                             ta_xxyyy_xyzzz_1,  \
                             ta_xxyyy_yyyyy_0,  \
                             ta_xxyyy_yyyyy_1,  \
                             ta_xxyyyz_xxxxx_0, \
                             ta_xxyyyz_xxxxy_0, \
                             ta_xxyyyz_xxxxz_0, \
                             ta_xxyyyz_xxxyy_0, \
                             ta_xxyyyz_xxxyz_0, \
                             ta_xxyyyz_xxxzz_0, \
                             ta_xxyyyz_xxyyy_0, \
                             ta_xxyyyz_xxyyz_0, \
                             ta_xxyyyz_xxyzz_0, \
                             ta_xxyyyz_xxzzz_0, \
                             ta_xxyyyz_xyyyy_0, \
                             ta_xxyyyz_xyyyz_0, \
                             ta_xxyyyz_xyyzz_0, \
                             ta_xxyyyz_xyzzz_0, \
                             ta_xxyyyz_xzzzz_0, \
                             ta_xxyyyz_yyyyy_0, \
                             ta_xxyyyz_yyyyz_0, \
                             ta_xxyyyz_yyyzz_0, \
                             ta_xxyyyz_yyzzz_0, \
                             ta_xxyyyz_yzzzz_0, \
                             ta_xxyyyz_zzzzz_0, \
                             ta_xxyyz_xxxxz_0,  \
                             ta_xxyyz_xxxxz_1,  \
                             ta_xxyyz_xxxzz_0,  \
                             ta_xxyyz_xxxzz_1,  \
                             ta_xxyyz_xxzzz_0,  \
                             ta_xxyyz_xxzzz_1,  \
                             ta_xxyyz_xzzzz_0,  \
                             ta_xxyyz_xzzzz_1,  \
                             ta_xxyz_xxxxz_0,   \
                             ta_xxyz_xxxxz_1,   \
                             ta_xxyz_xxxzz_0,   \
                             ta_xxyz_xxxzz_1,   \
                             ta_xxyz_xxzzz_0,   \
                             ta_xxyz_xxzzz_1,   \
                             ta_xxyz_xzzzz_0,   \
                             ta_xxyz_xzzzz_1,   \
                             ta_xyyyz_yyyyz_0,  \
                             ta_xyyyz_yyyyz_1,  \
                             ta_xyyyz_yyyzz_0,  \
                             ta_xyyyz_yyyzz_1,  \
                             ta_xyyyz_yyzzz_0,  \
                             ta_xyyyz_yyzzz_1,  \
                             ta_xyyyz_yzzzz_0,  \
                             ta_xyyyz_yzzzz_1,  \
                             ta_xyyyz_zzzzz_0,  \
                             ta_xyyyz_zzzzz_1,  \
                             ta_yyyz_yyyyz_0,   \
                             ta_yyyz_yyyyz_1,   \
                             ta_yyyz_yyyzz_0,   \
                             ta_yyyz_yyyzz_1,   \
                             ta_yyyz_yyzzz_0,   \
                             ta_yyyz_yyzzz_1,   \
                             ta_yyyz_yzzzz_0,   \
                             ta_yyyz_yzzzz_1,   \
                             ta_yyyz_zzzzz_0,   \
                             ta_yyyz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyz_xxxxx_0[i] = ta_xxyyy_xxxxx_0[i] * pa_z[i] - ta_xxyyy_xxxxx_1[i] * pc_z[i];

        ta_xxyyyz_xxxxy_0[i] = ta_xxyyy_xxxxy_0[i] * pa_z[i] - ta_xxyyy_xxxxy_1[i] * pc_z[i];

        ta_xxyyyz_xxxxz_0[i] =
            2.0 * ta_xxyz_xxxxz_0[i] * fe_0 - 2.0 * ta_xxyz_xxxxz_1[i] * fe_0 + ta_xxyyz_xxxxz_0[i] * pa_y[i] - ta_xxyyz_xxxxz_1[i] * pc_y[i];

        ta_xxyyyz_xxxyy_0[i] = ta_xxyyy_xxxyy_0[i] * pa_z[i] - ta_xxyyy_xxxyy_1[i] * pc_z[i];

        ta_xxyyyz_xxxyz_0[i] = ta_xxyyy_xxxy_0[i] * fe_0 - ta_xxyyy_xxxy_1[i] * fe_0 + ta_xxyyy_xxxyz_0[i] * pa_z[i] - ta_xxyyy_xxxyz_1[i] * pc_z[i];

        ta_xxyyyz_xxxzz_0[i] =
            2.0 * ta_xxyz_xxxzz_0[i] * fe_0 - 2.0 * ta_xxyz_xxxzz_1[i] * fe_0 + ta_xxyyz_xxxzz_0[i] * pa_y[i] - ta_xxyyz_xxxzz_1[i] * pc_y[i];

        ta_xxyyyz_xxyyy_0[i] = ta_xxyyy_xxyyy_0[i] * pa_z[i] - ta_xxyyy_xxyyy_1[i] * pc_z[i];

        ta_xxyyyz_xxyyz_0[i] = ta_xxyyy_xxyy_0[i] * fe_0 - ta_xxyyy_xxyy_1[i] * fe_0 + ta_xxyyy_xxyyz_0[i] * pa_z[i] - ta_xxyyy_xxyyz_1[i] * pc_z[i];

        ta_xxyyyz_xxyzz_0[i] =
            2.0 * ta_xxyyy_xxyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xxyz_1[i] * fe_0 + ta_xxyyy_xxyzz_0[i] * pa_z[i] - ta_xxyyy_xxyzz_1[i] * pc_z[i];

        ta_xxyyyz_xxzzz_0[i] =
            2.0 * ta_xxyz_xxzzz_0[i] * fe_0 - 2.0 * ta_xxyz_xxzzz_1[i] * fe_0 + ta_xxyyz_xxzzz_0[i] * pa_y[i] - ta_xxyyz_xxzzz_1[i] * pc_y[i];

        ta_xxyyyz_xyyyy_0[i] = ta_xxyyy_xyyyy_0[i] * pa_z[i] - ta_xxyyy_xyyyy_1[i] * pc_z[i];

        ta_xxyyyz_xyyyz_0[i] = ta_xxyyy_xyyy_0[i] * fe_0 - ta_xxyyy_xyyy_1[i] * fe_0 + ta_xxyyy_xyyyz_0[i] * pa_z[i] - ta_xxyyy_xyyyz_1[i] * pc_z[i];

        ta_xxyyyz_xyyzz_0[i] =
            2.0 * ta_xxyyy_xyyz_0[i] * fe_0 - 2.0 * ta_xxyyy_xyyz_1[i] * fe_0 + ta_xxyyy_xyyzz_0[i] * pa_z[i] - ta_xxyyy_xyyzz_1[i] * pc_z[i];

        ta_xxyyyz_xyzzz_0[i] =
            3.0 * ta_xxyyy_xyzz_0[i] * fe_0 - 3.0 * ta_xxyyy_xyzz_1[i] * fe_0 + ta_xxyyy_xyzzz_0[i] * pa_z[i] - ta_xxyyy_xyzzz_1[i] * pc_z[i];

        ta_xxyyyz_xzzzz_0[i] =
            2.0 * ta_xxyz_xzzzz_0[i] * fe_0 - 2.0 * ta_xxyz_xzzzz_1[i] * fe_0 + ta_xxyyz_xzzzz_0[i] * pa_y[i] - ta_xxyyz_xzzzz_1[i] * pc_y[i];

        ta_xxyyyz_yyyyy_0[i] = ta_xxyyy_yyyyy_0[i] * pa_z[i] - ta_xxyyy_yyyyy_1[i] * pc_z[i];

        ta_xxyyyz_yyyyz_0[i] = ta_yyyz_yyyyz_0[i] * fe_0 - ta_yyyz_yyyyz_1[i] * fe_0 + ta_xyyyz_yyyyz_0[i] * pa_x[i] - ta_xyyyz_yyyyz_1[i] * pc_x[i];

        ta_xxyyyz_yyyzz_0[i] = ta_yyyz_yyyzz_0[i] * fe_0 - ta_yyyz_yyyzz_1[i] * fe_0 + ta_xyyyz_yyyzz_0[i] * pa_x[i] - ta_xyyyz_yyyzz_1[i] * pc_x[i];

        ta_xxyyyz_yyzzz_0[i] = ta_yyyz_yyzzz_0[i] * fe_0 - ta_yyyz_yyzzz_1[i] * fe_0 + ta_xyyyz_yyzzz_0[i] * pa_x[i] - ta_xyyyz_yyzzz_1[i] * pc_x[i];

        ta_xxyyyz_yzzzz_0[i] = ta_yyyz_yzzzz_0[i] * fe_0 - ta_yyyz_yzzzz_1[i] * fe_0 + ta_xyyyz_yzzzz_0[i] * pa_x[i] - ta_xyyyz_yzzzz_1[i] * pc_x[i];

        ta_xxyyyz_zzzzz_0[i] = ta_yyyz_zzzzz_0[i] * fe_0 - ta_yyyz_zzzzz_1[i] * fe_0 + ta_xyyyz_zzzzz_0[i] * pa_x[i] - ta_xyyyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 252-273 components of targeted buffer : IH

    auto ta_xxyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 252);

    auto ta_xxyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 253);

    auto ta_xxyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 254);

    auto ta_xxyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 255);

    auto ta_xxyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 256);

    auto ta_xxyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 257);

    auto ta_xxyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 258);

    auto ta_xxyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 259);

    auto ta_xxyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 260);

    auto ta_xxyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 261);

    auto ta_xxyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 262);

    auto ta_xxyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 263);

    auto ta_xxyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 264);

    auto ta_xxyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 265);

    auto ta_xxyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 266);

    auto ta_xxyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 267);

    auto ta_xxyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 268);

    auto ta_xxyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 269);

    auto ta_xxyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 270);

    auto ta_xxyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 271);

    auto ta_xxyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 272);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta_xxyy_xxxxy_0,   \
                             ta_xxyy_xxxxy_1,   \
                             ta_xxyy_xxxyy_0,   \
                             ta_xxyy_xxxyy_1,   \
                             ta_xxyy_xxyyy_0,   \
                             ta_xxyy_xxyyy_1,   \
                             ta_xxyy_xyyyy_0,   \
                             ta_xxyy_xyyyy_1,   \
                             ta_xxyyz_xxxxy_0,  \
                             ta_xxyyz_xxxxy_1,  \
                             ta_xxyyz_xxxyy_0,  \
                             ta_xxyyz_xxxyy_1,  \
                             ta_xxyyz_xxyyy_0,  \
                             ta_xxyyz_xxyyy_1,  \
                             ta_xxyyz_xyyyy_0,  \
                             ta_xxyyz_xyyyy_1,  \
                             ta_xxyyzz_xxxxx_0, \
                             ta_xxyyzz_xxxxy_0, \
                             ta_xxyyzz_xxxxz_0, \
                             ta_xxyyzz_xxxyy_0, \
                             ta_xxyyzz_xxxyz_0, \
                             ta_xxyyzz_xxxzz_0, \
                             ta_xxyyzz_xxyyy_0, \
                             ta_xxyyzz_xxyyz_0, \
                             ta_xxyyzz_xxyzz_0, \
                             ta_xxyyzz_xxzzz_0, \
                             ta_xxyyzz_xyyyy_0, \
                             ta_xxyyzz_xyyyz_0, \
                             ta_xxyyzz_xyyzz_0, \
                             ta_xxyyzz_xyzzz_0, \
                             ta_xxyyzz_xzzzz_0, \
                             ta_xxyyzz_yyyyy_0, \
                             ta_xxyyzz_yyyyz_0, \
                             ta_xxyyzz_yyyzz_0, \
                             ta_xxyyzz_yyzzz_0, \
                             ta_xxyyzz_yzzzz_0, \
                             ta_xxyyzz_zzzzz_0, \
                             ta_xxyzz_xxxxx_0,  \
                             ta_xxyzz_xxxxx_1,  \
                             ta_xxyzz_xxxxz_0,  \
                             ta_xxyzz_xxxxz_1,  \
                             ta_xxyzz_xxxzz_0,  \
                             ta_xxyzz_xxxzz_1,  \
                             ta_xxyzz_xxzzz_0,  \
                             ta_xxyzz_xxzzz_1,  \
                             ta_xxyzz_xzzzz_0,  \
                             ta_xxyzz_xzzzz_1,  \
                             ta_xxzz_xxxxx_0,   \
                             ta_xxzz_xxxxx_1,   \
                             ta_xxzz_xxxxz_0,   \
                             ta_xxzz_xxxxz_1,   \
                             ta_xxzz_xxxzz_0,   \
                             ta_xxzz_xxxzz_1,   \
                             ta_xxzz_xxzzz_0,   \
                             ta_xxzz_xxzzz_1,   \
                             ta_xxzz_xzzzz_0,   \
                             ta_xxzz_xzzzz_1,   \
                             ta_xyyzz_xxxyz_0,  \
                             ta_xyyzz_xxxyz_1,  \
                             ta_xyyzz_xxyyz_0,  \
                             ta_xyyzz_xxyyz_1,  \
                             ta_xyyzz_xxyz_0,   \
                             ta_xyyzz_xxyz_1,   \
                             ta_xyyzz_xxyzz_0,  \
                             ta_xyyzz_xxyzz_1,  \
                             ta_xyyzz_xyyyz_0,  \
                             ta_xyyzz_xyyyz_1,  \
                             ta_xyyzz_xyyz_0,   \
                             ta_xyyzz_xyyz_1,   \
                             ta_xyyzz_xyyzz_0,  \
                             ta_xyyzz_xyyzz_1,  \
                             ta_xyyzz_xyzz_0,   \
                             ta_xyyzz_xyzz_1,   \
                             ta_xyyzz_xyzzz_0,  \
                             ta_xyyzz_xyzzz_1,  \
                             ta_xyyzz_yyyyy_0,  \
                             ta_xyyzz_yyyyy_1,  \
                             ta_xyyzz_yyyyz_0,  \
                             ta_xyyzz_yyyyz_1,  \
                             ta_xyyzz_yyyz_0,   \
                             ta_xyyzz_yyyz_1,   \
                             ta_xyyzz_yyyzz_0,  \
                             ta_xyyzz_yyyzz_1,  \
                             ta_xyyzz_yyzz_0,   \
                             ta_xyyzz_yyzz_1,   \
                             ta_xyyzz_yyzzz_0,  \
                             ta_xyyzz_yyzzz_1,  \
                             ta_xyyzz_yzzz_0,   \
                             ta_xyyzz_yzzz_1,   \
                             ta_xyyzz_yzzzz_0,  \
                             ta_xyyzz_yzzzz_1,  \
                             ta_xyyzz_zzzzz_0,  \
                             ta_xyyzz_zzzzz_1,  \
                             ta_yyzz_xxxyz_0,   \
                             ta_yyzz_xxxyz_1,   \
                             ta_yyzz_xxyyz_0,   \
                             ta_yyzz_xxyyz_1,   \
                             ta_yyzz_xxyzz_0,   \
                             ta_yyzz_xxyzz_1,   \
                             ta_yyzz_xyyyz_0,   \
                             ta_yyzz_xyyyz_1,   \
                             ta_yyzz_xyyzz_0,   \
                             ta_yyzz_xyyzz_1,   \
                             ta_yyzz_xyzzz_0,   \
                             ta_yyzz_xyzzz_1,   \
                             ta_yyzz_yyyyy_0,   \
                             ta_yyzz_yyyyy_1,   \
                             ta_yyzz_yyyyz_0,   \
                             ta_yyzz_yyyyz_1,   \
                             ta_yyzz_yyyzz_0,   \
                             ta_yyzz_yyyzz_1,   \
                             ta_yyzz_yyzzz_0,   \
                             ta_yyzz_yyzzz_1,   \
                             ta_yyzz_yzzzz_0,   \
                             ta_yyzz_yzzzz_1,   \
                             ta_yyzz_zzzzz_0,   \
                             ta_yyzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyzz_xxxxx_0[i] = ta_xxzz_xxxxx_0[i] * fe_0 - ta_xxzz_xxxxx_1[i] * fe_0 + ta_xxyzz_xxxxx_0[i] * pa_y[i] - ta_xxyzz_xxxxx_1[i] * pc_y[i];

        ta_xxyyzz_xxxxy_0[i] = ta_xxyy_xxxxy_0[i] * fe_0 - ta_xxyy_xxxxy_1[i] * fe_0 + ta_xxyyz_xxxxy_0[i] * pa_z[i] - ta_xxyyz_xxxxy_1[i] * pc_z[i];

        ta_xxyyzz_xxxxz_0[i] = ta_xxzz_xxxxz_0[i] * fe_0 - ta_xxzz_xxxxz_1[i] * fe_0 + ta_xxyzz_xxxxz_0[i] * pa_y[i] - ta_xxyzz_xxxxz_1[i] * pc_y[i];

        ta_xxyyzz_xxxyy_0[i] = ta_xxyy_xxxyy_0[i] * fe_0 - ta_xxyy_xxxyy_1[i] * fe_0 + ta_xxyyz_xxxyy_0[i] * pa_z[i] - ta_xxyyz_xxxyy_1[i] * pc_z[i];

        ta_xxyyzz_xxxyz_0[i] = ta_yyzz_xxxyz_0[i] * fe_0 - ta_yyzz_xxxyz_1[i] * fe_0 + 3.0 * ta_xyyzz_xxyz_0[i] * fe_0 -
                               3.0 * ta_xyyzz_xxyz_1[i] * fe_0 + ta_xyyzz_xxxyz_0[i] * pa_x[i] - ta_xyyzz_xxxyz_1[i] * pc_x[i];

        ta_xxyyzz_xxxzz_0[i] = ta_xxzz_xxxzz_0[i] * fe_0 - ta_xxzz_xxxzz_1[i] * fe_0 + ta_xxyzz_xxxzz_0[i] * pa_y[i] - ta_xxyzz_xxxzz_1[i] * pc_y[i];

        ta_xxyyzz_xxyyy_0[i] = ta_xxyy_xxyyy_0[i] * fe_0 - ta_xxyy_xxyyy_1[i] * fe_0 + ta_xxyyz_xxyyy_0[i] * pa_z[i] - ta_xxyyz_xxyyy_1[i] * pc_z[i];

        ta_xxyyzz_xxyyz_0[i] = ta_yyzz_xxyyz_0[i] * fe_0 - ta_yyzz_xxyyz_1[i] * fe_0 + 2.0 * ta_xyyzz_xyyz_0[i] * fe_0 -
                               2.0 * ta_xyyzz_xyyz_1[i] * fe_0 + ta_xyyzz_xxyyz_0[i] * pa_x[i] - ta_xyyzz_xxyyz_1[i] * pc_x[i];

        ta_xxyyzz_xxyzz_0[i] = ta_yyzz_xxyzz_0[i] * fe_0 - ta_yyzz_xxyzz_1[i] * fe_0 + 2.0 * ta_xyyzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_xyyzz_xyzz_1[i] * fe_0 + ta_xyyzz_xxyzz_0[i] * pa_x[i] - ta_xyyzz_xxyzz_1[i] * pc_x[i];

        ta_xxyyzz_xxzzz_0[i] = ta_xxzz_xxzzz_0[i] * fe_0 - ta_xxzz_xxzzz_1[i] * fe_0 + ta_xxyzz_xxzzz_0[i] * pa_y[i] - ta_xxyzz_xxzzz_1[i] * pc_y[i];

        ta_xxyyzz_xyyyy_0[i] = ta_xxyy_xyyyy_0[i] * fe_0 - ta_xxyy_xyyyy_1[i] * fe_0 + ta_xxyyz_xyyyy_0[i] * pa_z[i] - ta_xxyyz_xyyyy_1[i] * pc_z[i];

        ta_xxyyzz_xyyyz_0[i] = ta_yyzz_xyyyz_0[i] * fe_0 - ta_yyzz_xyyyz_1[i] * fe_0 + ta_xyyzz_yyyz_0[i] * fe_0 - ta_xyyzz_yyyz_1[i] * fe_0 +
                               ta_xyyzz_xyyyz_0[i] * pa_x[i] - ta_xyyzz_xyyyz_1[i] * pc_x[i];

        ta_xxyyzz_xyyzz_0[i] = ta_yyzz_xyyzz_0[i] * fe_0 - ta_yyzz_xyyzz_1[i] * fe_0 + ta_xyyzz_yyzz_0[i] * fe_0 - ta_xyyzz_yyzz_1[i] * fe_0 +
                               ta_xyyzz_xyyzz_0[i] * pa_x[i] - ta_xyyzz_xyyzz_1[i] * pc_x[i];

        ta_xxyyzz_xyzzz_0[i] = ta_yyzz_xyzzz_0[i] * fe_0 - ta_yyzz_xyzzz_1[i] * fe_0 + ta_xyyzz_yzzz_0[i] * fe_0 - ta_xyyzz_yzzz_1[i] * fe_0 +
                               ta_xyyzz_xyzzz_0[i] * pa_x[i] - ta_xyyzz_xyzzz_1[i] * pc_x[i];

        ta_xxyyzz_xzzzz_0[i] = ta_xxzz_xzzzz_0[i] * fe_0 - ta_xxzz_xzzzz_1[i] * fe_0 + ta_xxyzz_xzzzz_0[i] * pa_y[i] - ta_xxyzz_xzzzz_1[i] * pc_y[i];

        ta_xxyyzz_yyyyy_0[i] = ta_yyzz_yyyyy_0[i] * fe_0 - ta_yyzz_yyyyy_1[i] * fe_0 + ta_xyyzz_yyyyy_0[i] * pa_x[i] - ta_xyyzz_yyyyy_1[i] * pc_x[i];

        ta_xxyyzz_yyyyz_0[i] = ta_yyzz_yyyyz_0[i] * fe_0 - ta_yyzz_yyyyz_1[i] * fe_0 + ta_xyyzz_yyyyz_0[i] * pa_x[i] - ta_xyyzz_yyyyz_1[i] * pc_x[i];

        ta_xxyyzz_yyyzz_0[i] = ta_yyzz_yyyzz_0[i] * fe_0 - ta_yyzz_yyyzz_1[i] * fe_0 + ta_xyyzz_yyyzz_0[i] * pa_x[i] - ta_xyyzz_yyyzz_1[i] * pc_x[i];

        ta_xxyyzz_yyzzz_0[i] = ta_yyzz_yyzzz_0[i] * fe_0 - ta_yyzz_yyzzz_1[i] * fe_0 + ta_xyyzz_yyzzz_0[i] * pa_x[i] - ta_xyyzz_yyzzz_1[i] * pc_x[i];

        ta_xxyyzz_yzzzz_0[i] = ta_yyzz_yzzzz_0[i] * fe_0 - ta_yyzz_yzzzz_1[i] * fe_0 + ta_xyyzz_yzzzz_0[i] * pa_x[i] - ta_xyyzz_yzzzz_1[i] * pc_x[i];

        ta_xxyyzz_zzzzz_0[i] = ta_yyzz_zzzzz_0[i] * fe_0 - ta_yyzz_zzzzz_1[i] * fe_0 + ta_xyyzz_zzzzz_0[i] * pa_x[i] - ta_xyyzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 273-294 components of targeted buffer : IH

    auto ta_xxyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 273);

    auto ta_xxyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 274);

    auto ta_xxyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 275);

    auto ta_xxyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 276);

    auto ta_xxyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 277);

    auto ta_xxyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 278);

    auto ta_xxyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 279);

    auto ta_xxyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 280);

    auto ta_xxyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 281);

    auto ta_xxyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 282);

    auto ta_xxyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 283);

    auto ta_xxyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 284);

    auto ta_xxyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 285);

    auto ta_xxyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 286);

    auto ta_xxyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 287);

    auto ta_xxyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 288);

    auto ta_xxyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 289);

    auto ta_xxyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 290);

    auto ta_xxyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 291);

    auto ta_xxyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 292);

    auto ta_xxyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 293);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxyzzz_xxxxx_0, \
                             ta_xxyzzz_xxxxy_0, \
                             ta_xxyzzz_xxxxz_0, \
                             ta_xxyzzz_xxxyy_0, \
                             ta_xxyzzz_xxxyz_0, \
                             ta_xxyzzz_xxxzz_0, \
                             ta_xxyzzz_xxyyy_0, \
                             ta_xxyzzz_xxyyz_0, \
                             ta_xxyzzz_xxyzz_0, \
                             ta_xxyzzz_xxzzz_0, \
                             ta_xxyzzz_xyyyy_0, \
                             ta_xxyzzz_xyyyz_0, \
                             ta_xxyzzz_xyyzz_0, \
                             ta_xxyzzz_xyzzz_0, \
                             ta_xxyzzz_xzzzz_0, \
                             ta_xxyzzz_yyyyy_0, \
                             ta_xxyzzz_yyyyz_0, \
                             ta_xxyzzz_yyyzz_0, \
                             ta_xxyzzz_yyzzz_0, \
                             ta_xxyzzz_yzzzz_0, \
                             ta_xxyzzz_zzzzz_0, \
                             ta_xxzzz_xxxx_0,   \
                             ta_xxzzz_xxxx_1,   \
                             ta_xxzzz_xxxxx_0,  \
                             ta_xxzzz_xxxxx_1,  \
                             ta_xxzzz_xxxxy_0,  \
                             ta_xxzzz_xxxxy_1,  \
                             ta_xxzzz_xxxxz_0,  \
                             ta_xxzzz_xxxxz_1,  \
                             ta_xxzzz_xxxy_0,   \
                             ta_xxzzz_xxxy_1,   \
                             ta_xxzzz_xxxyy_0,  \
                             ta_xxzzz_xxxyy_1,  \
                             ta_xxzzz_xxxyz_0,  \
                             ta_xxzzz_xxxyz_1,  \
                             ta_xxzzz_xxxz_0,   \
                             ta_xxzzz_xxxz_1,   \
                             ta_xxzzz_xxxzz_0,  \
                             ta_xxzzz_xxxzz_1,  \
                             ta_xxzzz_xxyy_0,   \
                             ta_xxzzz_xxyy_1,   \
                             ta_xxzzz_xxyyy_0,  \
                             ta_xxzzz_xxyyy_1,  \
                             ta_xxzzz_xxyyz_0,  \
                             ta_xxzzz_xxyyz_1,  \
                             ta_xxzzz_xxyz_0,   \
                             ta_xxzzz_xxyz_1,   \
                             ta_xxzzz_xxyzz_0,  \
                             ta_xxzzz_xxyzz_1,  \
                             ta_xxzzz_xxzz_0,   \
                             ta_xxzzz_xxzz_1,   \
                             ta_xxzzz_xxzzz_0,  \
                             ta_xxzzz_xxzzz_1,  \
                             ta_xxzzz_xyyy_0,   \
                             ta_xxzzz_xyyy_1,   \
                             ta_xxzzz_xyyyy_0,  \
                             ta_xxzzz_xyyyy_1,  \
                             ta_xxzzz_xyyyz_0,  \
                             ta_xxzzz_xyyyz_1,  \
                             ta_xxzzz_xyyz_0,   \
                             ta_xxzzz_xyyz_1,   \
                             ta_xxzzz_xyyzz_0,  \
                             ta_xxzzz_xyyzz_1,  \
                             ta_xxzzz_xyzz_0,   \
                             ta_xxzzz_xyzz_1,   \
                             ta_xxzzz_xyzzz_0,  \
                             ta_xxzzz_xyzzz_1,  \
                             ta_xxzzz_xzzz_0,   \
                             ta_xxzzz_xzzz_1,   \
                             ta_xxzzz_xzzzz_0,  \
                             ta_xxzzz_xzzzz_1,  \
                             ta_xxzzz_zzzzz_0,  \
                             ta_xxzzz_zzzzz_1,  \
                             ta_xyzzz_yyyyy_0,  \
                             ta_xyzzz_yyyyy_1,  \
                             ta_xyzzz_yyyyz_0,  \
                             ta_xyzzz_yyyyz_1,  \
                             ta_xyzzz_yyyzz_0,  \
                             ta_xyzzz_yyyzz_1,  \
                             ta_xyzzz_yyzzz_0,  \
                             ta_xyzzz_yyzzz_1,  \
                             ta_xyzzz_yzzzz_0,  \
                             ta_xyzzz_yzzzz_1,  \
                             ta_yzzz_yyyyy_0,   \
                             ta_yzzz_yyyyy_1,   \
                             ta_yzzz_yyyyz_0,   \
                             ta_yzzz_yyyyz_1,   \
                             ta_yzzz_yyyzz_0,   \
                             ta_yzzz_yyyzz_1,   \
                             ta_yzzz_yyzzz_0,   \
                             ta_yzzz_yyzzz_1,   \
                             ta_yzzz_yzzzz_0,   \
                             ta_yzzz_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzzz_xxxxx_0[i] = ta_xxzzz_xxxxx_0[i] * pa_y[i] - ta_xxzzz_xxxxx_1[i] * pc_y[i];

        ta_xxyzzz_xxxxy_0[i] = ta_xxzzz_xxxx_0[i] * fe_0 - ta_xxzzz_xxxx_1[i] * fe_0 + ta_xxzzz_xxxxy_0[i] * pa_y[i] - ta_xxzzz_xxxxy_1[i] * pc_y[i];

        ta_xxyzzz_xxxxz_0[i] = ta_xxzzz_xxxxz_0[i] * pa_y[i] - ta_xxzzz_xxxxz_1[i] * pc_y[i];

        ta_xxyzzz_xxxyy_0[i] =
            2.0 * ta_xxzzz_xxxy_0[i] * fe_0 - 2.0 * ta_xxzzz_xxxy_1[i] * fe_0 + ta_xxzzz_xxxyy_0[i] * pa_y[i] - ta_xxzzz_xxxyy_1[i] * pc_y[i];

        ta_xxyzzz_xxxyz_0[i] = ta_xxzzz_xxxz_0[i] * fe_0 - ta_xxzzz_xxxz_1[i] * fe_0 + ta_xxzzz_xxxyz_0[i] * pa_y[i] - ta_xxzzz_xxxyz_1[i] * pc_y[i];

        ta_xxyzzz_xxxzz_0[i] = ta_xxzzz_xxxzz_0[i] * pa_y[i] - ta_xxzzz_xxxzz_1[i] * pc_y[i];

        ta_xxyzzz_xxyyy_0[i] =
            3.0 * ta_xxzzz_xxyy_0[i] * fe_0 - 3.0 * ta_xxzzz_xxyy_1[i] * fe_0 + ta_xxzzz_xxyyy_0[i] * pa_y[i] - ta_xxzzz_xxyyy_1[i] * pc_y[i];

        ta_xxyzzz_xxyyz_0[i] =
            2.0 * ta_xxzzz_xxyz_0[i] * fe_0 - 2.0 * ta_xxzzz_xxyz_1[i] * fe_0 + ta_xxzzz_xxyyz_0[i] * pa_y[i] - ta_xxzzz_xxyyz_1[i] * pc_y[i];

        ta_xxyzzz_xxyzz_0[i] = ta_xxzzz_xxzz_0[i] * fe_0 - ta_xxzzz_xxzz_1[i] * fe_0 + ta_xxzzz_xxyzz_0[i] * pa_y[i] - ta_xxzzz_xxyzz_1[i] * pc_y[i];

        ta_xxyzzz_xxzzz_0[i] = ta_xxzzz_xxzzz_0[i] * pa_y[i] - ta_xxzzz_xxzzz_1[i] * pc_y[i];

        ta_xxyzzz_xyyyy_0[i] =
            4.0 * ta_xxzzz_xyyy_0[i] * fe_0 - 4.0 * ta_xxzzz_xyyy_1[i] * fe_0 + ta_xxzzz_xyyyy_0[i] * pa_y[i] - ta_xxzzz_xyyyy_1[i] * pc_y[i];

        ta_xxyzzz_xyyyz_0[i] =
            3.0 * ta_xxzzz_xyyz_0[i] * fe_0 - 3.0 * ta_xxzzz_xyyz_1[i] * fe_0 + ta_xxzzz_xyyyz_0[i] * pa_y[i] - ta_xxzzz_xyyyz_1[i] * pc_y[i];

        ta_xxyzzz_xyyzz_0[i] =
            2.0 * ta_xxzzz_xyzz_0[i] * fe_0 - 2.0 * ta_xxzzz_xyzz_1[i] * fe_0 + ta_xxzzz_xyyzz_0[i] * pa_y[i] - ta_xxzzz_xyyzz_1[i] * pc_y[i];

        ta_xxyzzz_xyzzz_0[i] = ta_xxzzz_xzzz_0[i] * fe_0 - ta_xxzzz_xzzz_1[i] * fe_0 + ta_xxzzz_xyzzz_0[i] * pa_y[i] - ta_xxzzz_xyzzz_1[i] * pc_y[i];

        ta_xxyzzz_xzzzz_0[i] = ta_xxzzz_xzzzz_0[i] * pa_y[i] - ta_xxzzz_xzzzz_1[i] * pc_y[i];

        ta_xxyzzz_yyyyy_0[i] = ta_yzzz_yyyyy_0[i] * fe_0 - ta_yzzz_yyyyy_1[i] * fe_0 + ta_xyzzz_yyyyy_0[i] * pa_x[i] - ta_xyzzz_yyyyy_1[i] * pc_x[i];

        ta_xxyzzz_yyyyz_0[i] = ta_yzzz_yyyyz_0[i] * fe_0 - ta_yzzz_yyyyz_1[i] * fe_0 + ta_xyzzz_yyyyz_0[i] * pa_x[i] - ta_xyzzz_yyyyz_1[i] * pc_x[i];

        ta_xxyzzz_yyyzz_0[i] = ta_yzzz_yyyzz_0[i] * fe_0 - ta_yzzz_yyyzz_1[i] * fe_0 + ta_xyzzz_yyyzz_0[i] * pa_x[i] - ta_xyzzz_yyyzz_1[i] * pc_x[i];

        ta_xxyzzz_yyzzz_0[i] = ta_yzzz_yyzzz_0[i] * fe_0 - ta_yzzz_yyzzz_1[i] * fe_0 + ta_xyzzz_yyzzz_0[i] * pa_x[i] - ta_xyzzz_yyzzz_1[i] * pc_x[i];

        ta_xxyzzz_yzzzz_0[i] = ta_yzzz_yzzzz_0[i] * fe_0 - ta_yzzz_yzzzz_1[i] * fe_0 + ta_xyzzz_yzzzz_0[i] * pa_x[i] - ta_xyzzz_yzzzz_1[i] * pc_x[i];

        ta_xxyzzz_zzzzz_0[i] = ta_xxzzz_zzzzz_0[i] * pa_y[i] - ta_xxzzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 294-315 components of targeted buffer : IH

    auto ta_xxzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 294);

    auto ta_xxzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 295);

    auto ta_xxzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 296);

    auto ta_xxzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 297);

    auto ta_xxzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 298);

    auto ta_xxzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 299);

    auto ta_xxzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 300);

    auto ta_xxzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 301);

    auto ta_xxzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 302);

    auto ta_xxzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 303);

    auto ta_xxzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 304);

    auto ta_xxzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 305);

    auto ta_xxzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 306);

    auto ta_xxzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 307);

    auto ta_xxzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 308);

    auto ta_xxzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 309);

    auto ta_xxzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 310);

    auto ta_xxzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 311);

    auto ta_xxzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 312);

    auto ta_xxzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 313);

    auto ta_xxzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 314);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxzz_xxxxx_0,   \
                             ta_xxzz_xxxxx_1,   \
                             ta_xxzz_xxxxy_0,   \
                             ta_xxzz_xxxxy_1,   \
                             ta_xxzz_xxxyy_0,   \
                             ta_xxzz_xxxyy_1,   \
                             ta_xxzz_xxyyy_0,   \
                             ta_xxzz_xxyyy_1,   \
                             ta_xxzz_xyyyy_0,   \
                             ta_xxzz_xyyyy_1,   \
                             ta_xxzzz_xxxxx_0,  \
                             ta_xxzzz_xxxxx_1,  \
                             ta_xxzzz_xxxxy_0,  \
                             ta_xxzzz_xxxxy_1,  \
                             ta_xxzzz_xxxyy_0,  \
                             ta_xxzzz_xxxyy_1,  \
                             ta_xxzzz_xxyyy_0,  \
                             ta_xxzzz_xxyyy_1,  \
                             ta_xxzzz_xyyyy_0,  \
                             ta_xxzzz_xyyyy_1,  \
                             ta_xxzzzz_xxxxx_0, \
                             ta_xxzzzz_xxxxy_0, \
                             ta_xxzzzz_xxxxz_0, \
                             ta_xxzzzz_xxxyy_0, \
                             ta_xxzzzz_xxxyz_0, \
                             ta_xxzzzz_xxxzz_0, \
                             ta_xxzzzz_xxyyy_0, \
                             ta_xxzzzz_xxyyz_0, \
                             ta_xxzzzz_xxyzz_0, \
                             ta_xxzzzz_xxzzz_0, \
                             ta_xxzzzz_xyyyy_0, \
                             ta_xxzzzz_xyyyz_0, \
                             ta_xxzzzz_xyyzz_0, \
                             ta_xxzzzz_xyzzz_0, \
                             ta_xxzzzz_xzzzz_0, \
                             ta_xxzzzz_yyyyy_0, \
                             ta_xxzzzz_yyyyz_0, \
                             ta_xxzzzz_yyyzz_0, \
                             ta_xxzzzz_yyzzz_0, \
                             ta_xxzzzz_yzzzz_0, \
                             ta_xxzzzz_zzzzz_0, \
                             ta_xzzzz_xxxxz_0,  \
                             ta_xzzzz_xxxxz_1,  \
                             ta_xzzzz_xxxyz_0,  \
                             ta_xzzzz_xxxyz_1,  \
                             ta_xzzzz_xxxz_0,   \
                             ta_xzzzz_xxxz_1,   \
                             ta_xzzzz_xxxzz_0,  \
                             ta_xzzzz_xxxzz_1,  \
                             ta_xzzzz_xxyyz_0,  \
                             ta_xzzzz_xxyyz_1,  \
                             ta_xzzzz_xxyz_0,   \
                             ta_xzzzz_xxyz_1,   \
                             ta_xzzzz_xxyzz_0,  \
                             ta_xzzzz_xxyzz_1,  \
                             ta_xzzzz_xxzz_0,   \
                             ta_xzzzz_xxzz_1,   \
                             ta_xzzzz_xxzzz_0,  \
                             ta_xzzzz_xxzzz_1,  \
                             ta_xzzzz_xyyyz_0,  \
                             ta_xzzzz_xyyyz_1,  \
                             ta_xzzzz_xyyz_0,   \
                             ta_xzzzz_xyyz_1,   \
                             ta_xzzzz_xyyzz_0,  \
                             ta_xzzzz_xyyzz_1,  \
                             ta_xzzzz_xyzz_0,   \
                             ta_xzzzz_xyzz_1,   \
                             ta_xzzzz_xyzzz_0,  \
                             ta_xzzzz_xyzzz_1,  \
                             ta_xzzzz_xzzz_0,   \
                             ta_xzzzz_xzzz_1,   \
                             ta_xzzzz_xzzzz_0,  \
                             ta_xzzzz_xzzzz_1,  \
                             ta_xzzzz_yyyyy_0,  \
                             ta_xzzzz_yyyyy_1,  \
                             ta_xzzzz_yyyyz_0,  \
                             ta_xzzzz_yyyyz_1,  \
                             ta_xzzzz_yyyz_0,   \
                             ta_xzzzz_yyyz_1,   \
                             ta_xzzzz_yyyzz_0,  \
                             ta_xzzzz_yyyzz_1,  \
                             ta_xzzzz_yyzz_0,   \
                             ta_xzzzz_yyzz_1,   \
                             ta_xzzzz_yyzzz_0,  \
                             ta_xzzzz_yyzzz_1,  \
                             ta_xzzzz_yzzz_0,   \
                             ta_xzzzz_yzzz_1,   \
                             ta_xzzzz_yzzzz_0,  \
                             ta_xzzzz_yzzzz_1,  \
                             ta_xzzzz_zzzz_0,   \
                             ta_xzzzz_zzzz_1,   \
                             ta_xzzzz_zzzzz_0,  \
                             ta_xzzzz_zzzzz_1,  \
                             ta_zzzz_xxxxz_0,   \
                             ta_zzzz_xxxxz_1,   \
                             ta_zzzz_xxxyz_0,   \
                             ta_zzzz_xxxyz_1,   \
                             ta_zzzz_xxxzz_0,   \
                             ta_zzzz_xxxzz_1,   \
                             ta_zzzz_xxyyz_0,   \
                             ta_zzzz_xxyyz_1,   \
                             ta_zzzz_xxyzz_0,   \
                             ta_zzzz_xxyzz_1,   \
                             ta_zzzz_xxzzz_0,   \
                             ta_zzzz_xxzzz_1,   \
                             ta_zzzz_xyyyz_0,   \
                             ta_zzzz_xyyyz_1,   \
                             ta_zzzz_xyyzz_0,   \
                             ta_zzzz_xyyzz_1,   \
                             ta_zzzz_xyzzz_0,   \
                             ta_zzzz_xyzzz_1,   \
                             ta_zzzz_xzzzz_0,   \
                             ta_zzzz_xzzzz_1,   \
                             ta_zzzz_yyyyy_0,   \
                             ta_zzzz_yyyyy_1,   \
                             ta_zzzz_yyyyz_0,   \
                             ta_zzzz_yyyyz_1,   \
                             ta_zzzz_yyyzz_0,   \
                             ta_zzzz_yyyzz_1,   \
                             ta_zzzz_yyzzz_0,   \
                             ta_zzzz_yyzzz_1,   \
                             ta_zzzz_yzzzz_0,   \
                             ta_zzzz_yzzzz_1,   \
                             ta_zzzz_zzzzz_0,   \
                             ta_zzzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzzz_xxxxx_0[i] =
            3.0 * ta_xxzz_xxxxx_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxx_1[i] * fe_0 + ta_xxzzz_xxxxx_0[i] * pa_z[i] - ta_xxzzz_xxxxx_1[i] * pc_z[i];

        ta_xxzzzz_xxxxy_0[i] =
            3.0 * ta_xxzz_xxxxy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxxy_1[i] * fe_0 + ta_xxzzz_xxxxy_0[i] * pa_z[i] - ta_xxzzz_xxxxy_1[i] * pc_z[i];

        ta_xxzzzz_xxxxz_0[i] = ta_zzzz_xxxxz_0[i] * fe_0 - ta_zzzz_xxxxz_1[i] * fe_0 + 4.0 * ta_xzzzz_xxxz_0[i] * fe_0 -
                               4.0 * ta_xzzzz_xxxz_1[i] * fe_0 + ta_xzzzz_xxxxz_0[i] * pa_x[i] - ta_xzzzz_xxxxz_1[i] * pc_x[i];

        ta_xxzzzz_xxxyy_0[i] =
            3.0 * ta_xxzz_xxxyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxyy_1[i] * fe_0 + ta_xxzzz_xxxyy_0[i] * pa_z[i] - ta_xxzzz_xxxyy_1[i] * pc_z[i];

        ta_xxzzzz_xxxyz_0[i] = ta_zzzz_xxxyz_0[i] * fe_0 - ta_zzzz_xxxyz_1[i] * fe_0 + 3.0 * ta_xzzzz_xxyz_0[i] * fe_0 -
                               3.0 * ta_xzzzz_xxyz_1[i] * fe_0 + ta_xzzzz_xxxyz_0[i] * pa_x[i] - ta_xzzzz_xxxyz_1[i] * pc_x[i];

        ta_xxzzzz_xxxzz_0[i] = ta_zzzz_xxxzz_0[i] * fe_0 - ta_zzzz_xxxzz_1[i] * fe_0 + 3.0 * ta_xzzzz_xxzz_0[i] * fe_0 -
                               3.0 * ta_xzzzz_xxzz_1[i] * fe_0 + ta_xzzzz_xxxzz_0[i] * pa_x[i] - ta_xzzzz_xxxzz_1[i] * pc_x[i];

        ta_xxzzzz_xxyyy_0[i] =
            3.0 * ta_xxzz_xxyyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxyyy_1[i] * fe_0 + ta_xxzzz_xxyyy_0[i] * pa_z[i] - ta_xxzzz_xxyyy_1[i] * pc_z[i];

        ta_xxzzzz_xxyyz_0[i] = ta_zzzz_xxyyz_0[i] * fe_0 - ta_zzzz_xxyyz_1[i] * fe_0 + 2.0 * ta_xzzzz_xyyz_0[i] * fe_0 -
                               2.0 * ta_xzzzz_xyyz_1[i] * fe_0 + ta_xzzzz_xxyyz_0[i] * pa_x[i] - ta_xzzzz_xxyyz_1[i] * pc_x[i];

        ta_xxzzzz_xxyzz_0[i] = ta_zzzz_xxyzz_0[i] * fe_0 - ta_zzzz_xxyzz_1[i] * fe_0 + 2.0 * ta_xzzzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_xzzzz_xyzz_1[i] * fe_0 + ta_xzzzz_xxyzz_0[i] * pa_x[i] - ta_xzzzz_xxyzz_1[i] * pc_x[i];

        ta_xxzzzz_xxzzz_0[i] = ta_zzzz_xxzzz_0[i] * fe_0 - ta_zzzz_xxzzz_1[i] * fe_0 + 2.0 * ta_xzzzz_xzzz_0[i] * fe_0 -
                               2.0 * ta_xzzzz_xzzz_1[i] * fe_0 + ta_xzzzz_xxzzz_0[i] * pa_x[i] - ta_xzzzz_xxzzz_1[i] * pc_x[i];

        ta_xxzzzz_xyyyy_0[i] =
            3.0 * ta_xxzz_xyyyy_0[i] * fe_0 - 3.0 * ta_xxzz_xyyyy_1[i] * fe_0 + ta_xxzzz_xyyyy_0[i] * pa_z[i] - ta_xxzzz_xyyyy_1[i] * pc_z[i];

        ta_xxzzzz_xyyyz_0[i] = ta_zzzz_xyyyz_0[i] * fe_0 - ta_zzzz_xyyyz_1[i] * fe_0 + ta_xzzzz_yyyz_0[i] * fe_0 - ta_xzzzz_yyyz_1[i] * fe_0 +
                               ta_xzzzz_xyyyz_0[i] * pa_x[i] - ta_xzzzz_xyyyz_1[i] * pc_x[i];

        ta_xxzzzz_xyyzz_0[i] = ta_zzzz_xyyzz_0[i] * fe_0 - ta_zzzz_xyyzz_1[i] * fe_0 + ta_xzzzz_yyzz_0[i] * fe_0 - ta_xzzzz_yyzz_1[i] * fe_0 +
                               ta_xzzzz_xyyzz_0[i] * pa_x[i] - ta_xzzzz_xyyzz_1[i] * pc_x[i];

        ta_xxzzzz_xyzzz_0[i] = ta_zzzz_xyzzz_0[i] * fe_0 - ta_zzzz_xyzzz_1[i] * fe_0 + ta_xzzzz_yzzz_0[i] * fe_0 - ta_xzzzz_yzzz_1[i] * fe_0 +
                               ta_xzzzz_xyzzz_0[i] * pa_x[i] - ta_xzzzz_xyzzz_1[i] * pc_x[i];

        ta_xxzzzz_xzzzz_0[i] = ta_zzzz_xzzzz_0[i] * fe_0 - ta_zzzz_xzzzz_1[i] * fe_0 + ta_xzzzz_zzzz_0[i] * fe_0 - ta_xzzzz_zzzz_1[i] * fe_0 +
                               ta_xzzzz_xzzzz_0[i] * pa_x[i] - ta_xzzzz_xzzzz_1[i] * pc_x[i];

        ta_xxzzzz_yyyyy_0[i] = ta_zzzz_yyyyy_0[i] * fe_0 - ta_zzzz_yyyyy_1[i] * fe_0 + ta_xzzzz_yyyyy_0[i] * pa_x[i] - ta_xzzzz_yyyyy_1[i] * pc_x[i];

        ta_xxzzzz_yyyyz_0[i] = ta_zzzz_yyyyz_0[i] * fe_0 - ta_zzzz_yyyyz_1[i] * fe_0 + ta_xzzzz_yyyyz_0[i] * pa_x[i] - ta_xzzzz_yyyyz_1[i] * pc_x[i];

        ta_xxzzzz_yyyzz_0[i] = ta_zzzz_yyyzz_0[i] * fe_0 - ta_zzzz_yyyzz_1[i] * fe_0 + ta_xzzzz_yyyzz_0[i] * pa_x[i] - ta_xzzzz_yyyzz_1[i] * pc_x[i];

        ta_xxzzzz_yyzzz_0[i] = ta_zzzz_yyzzz_0[i] * fe_0 - ta_zzzz_yyzzz_1[i] * fe_0 + ta_xzzzz_yyzzz_0[i] * pa_x[i] - ta_xzzzz_yyzzz_1[i] * pc_x[i];

        ta_xxzzzz_yzzzz_0[i] = ta_zzzz_yzzzz_0[i] * fe_0 - ta_zzzz_yzzzz_1[i] * fe_0 + ta_xzzzz_yzzzz_0[i] * pa_x[i] - ta_xzzzz_yzzzz_1[i] * pc_x[i];

        ta_xxzzzz_zzzzz_0[i] = ta_zzzz_zzzzz_0[i] * fe_0 - ta_zzzz_zzzzz_1[i] * fe_0 + ta_xzzzz_zzzzz_0[i] * pa_x[i] - ta_xzzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 315-336 components of targeted buffer : IH

    auto ta_xyyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 315);

    auto ta_xyyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 316);

    auto ta_xyyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 317);

    auto ta_xyyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 318);

    auto ta_xyyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 319);

    auto ta_xyyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 320);

    auto ta_xyyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 321);

    auto ta_xyyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 322);

    auto ta_xyyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 323);

    auto ta_xyyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 324);

    auto ta_xyyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 325);

    auto ta_xyyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 326);

    auto ta_xyyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 327);

    auto ta_xyyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 328);

    auto ta_xyyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 329);

    auto ta_xyyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 330);

    auto ta_xyyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 331);

    auto ta_xyyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 332);

    auto ta_xyyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 333);

    auto ta_xyyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 334);

    auto ta_xyyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 335);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xyyyyy_xxxxx_0, \
                             ta_xyyyyy_xxxxy_0, \
                             ta_xyyyyy_xxxxz_0, \
                             ta_xyyyyy_xxxyy_0, \
                             ta_xyyyyy_xxxyz_0, \
                             ta_xyyyyy_xxxzz_0, \
                             ta_xyyyyy_xxyyy_0, \
                             ta_xyyyyy_xxyyz_0, \
                             ta_xyyyyy_xxyzz_0, \
                             ta_xyyyyy_xxzzz_0, \
                             ta_xyyyyy_xyyyy_0, \
                             ta_xyyyyy_xyyyz_0, \
                             ta_xyyyyy_xyyzz_0, \
                             ta_xyyyyy_xyzzz_0, \
                             ta_xyyyyy_xzzzz_0, \
                             ta_xyyyyy_yyyyy_0, \
                             ta_xyyyyy_yyyyz_0, \
                             ta_xyyyyy_yyyzz_0, \
                             ta_xyyyyy_yyzzz_0, \
                             ta_xyyyyy_yzzzz_0, \
                             ta_xyyyyy_zzzzz_0, \
                             ta_yyyyy_xxxx_0,   \
                             ta_yyyyy_xxxx_1,   \
                             ta_yyyyy_xxxxx_0,  \
                             ta_yyyyy_xxxxx_1,  \
                             ta_yyyyy_xxxxy_0,  \
                             ta_yyyyy_xxxxy_1,  \
                             ta_yyyyy_xxxxz_0,  \
                             ta_yyyyy_xxxxz_1,  \
                             ta_yyyyy_xxxy_0,   \
                             ta_yyyyy_xxxy_1,   \
                             ta_yyyyy_xxxyy_0,  \
                             ta_yyyyy_xxxyy_1,  \
                             ta_yyyyy_xxxyz_0,  \
                             ta_yyyyy_xxxyz_1,  \
                             ta_yyyyy_xxxz_0,   \
                             ta_yyyyy_xxxz_1,   \
                             ta_yyyyy_xxxzz_0,  \
                             ta_yyyyy_xxxzz_1,  \
                             ta_yyyyy_xxyy_0,   \
                             ta_yyyyy_xxyy_1,   \
                             ta_yyyyy_xxyyy_0,  \
                             ta_yyyyy_xxyyy_1,  \
                             ta_yyyyy_xxyyz_0,  \
                             ta_yyyyy_xxyyz_1,  \
                             ta_yyyyy_xxyz_0,   \
                             ta_yyyyy_xxyz_1,   \
                             ta_yyyyy_xxyzz_0,  \
                             ta_yyyyy_xxyzz_1,  \
                             ta_yyyyy_xxzz_0,   \
                             ta_yyyyy_xxzz_1,   \
                             ta_yyyyy_xxzzz_0,  \
                             ta_yyyyy_xxzzz_1,  \
                             ta_yyyyy_xyyy_0,   \
                             ta_yyyyy_xyyy_1,   \
                             ta_yyyyy_xyyyy_0,  \
                             ta_yyyyy_xyyyy_1,  \
                             ta_yyyyy_xyyyz_0,  \
                             ta_yyyyy_xyyyz_1,  \
                             ta_yyyyy_xyyz_0,   \
                             ta_yyyyy_xyyz_1,   \
                             ta_yyyyy_xyyzz_0,  \
                             ta_yyyyy_xyyzz_1,  \
                             ta_yyyyy_xyzz_0,   \
                             ta_yyyyy_xyzz_1,   \
                             ta_yyyyy_xyzzz_0,  \
                             ta_yyyyy_xyzzz_1,  \
                             ta_yyyyy_xzzz_0,   \
                             ta_yyyyy_xzzz_1,   \
                             ta_yyyyy_xzzzz_0,  \
                             ta_yyyyy_xzzzz_1,  \
                             ta_yyyyy_yyyy_0,   \
                             ta_yyyyy_yyyy_1,   \
                             ta_yyyyy_yyyyy_0,  \
                             ta_yyyyy_yyyyy_1,  \
                             ta_yyyyy_yyyyz_0,  \
                             ta_yyyyy_yyyyz_1,  \
                             ta_yyyyy_yyyz_0,   \
                             ta_yyyyy_yyyz_1,   \
                             ta_yyyyy_yyyzz_0,  \
                             ta_yyyyy_yyyzz_1,  \
                             ta_yyyyy_yyzz_0,   \
                             ta_yyyyy_yyzz_1,   \
                             ta_yyyyy_yyzzz_0,  \
                             ta_yyyyy_yyzzz_1,  \
                             ta_yyyyy_yzzz_0,   \
                             ta_yyyyy_yzzz_1,   \
                             ta_yyyyy_yzzzz_0,  \
                             ta_yyyyy_yzzzz_1,  \
                             ta_yyyyy_zzzz_0,   \
                             ta_yyyyy_zzzz_1,   \
                             ta_yyyyy_zzzzz_0,  \
                             ta_yyyyy_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyy_xxxxx_0[i] =
            5.0 * ta_yyyyy_xxxx_0[i] * fe_0 - 5.0 * ta_yyyyy_xxxx_1[i] * fe_0 + ta_yyyyy_xxxxx_0[i] * pa_x[i] - ta_yyyyy_xxxxx_1[i] * pc_x[i];

        ta_xyyyyy_xxxxy_0[i] =
            4.0 * ta_yyyyy_xxxy_0[i] * fe_0 - 4.0 * ta_yyyyy_xxxy_1[i] * fe_0 + ta_yyyyy_xxxxy_0[i] * pa_x[i] - ta_yyyyy_xxxxy_1[i] * pc_x[i];

        ta_xyyyyy_xxxxz_0[i] =
            4.0 * ta_yyyyy_xxxz_0[i] * fe_0 - 4.0 * ta_yyyyy_xxxz_1[i] * fe_0 + ta_yyyyy_xxxxz_0[i] * pa_x[i] - ta_yyyyy_xxxxz_1[i] * pc_x[i];

        ta_xyyyyy_xxxyy_0[i] =
            3.0 * ta_yyyyy_xxyy_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyy_1[i] * fe_0 + ta_yyyyy_xxxyy_0[i] * pa_x[i] - ta_yyyyy_xxxyy_1[i] * pc_x[i];

        ta_xyyyyy_xxxyz_0[i] =
            3.0 * ta_yyyyy_xxyz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxyz_1[i] * fe_0 + ta_yyyyy_xxxyz_0[i] * pa_x[i] - ta_yyyyy_xxxyz_1[i] * pc_x[i];

        ta_xyyyyy_xxxzz_0[i] =
            3.0 * ta_yyyyy_xxzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xxzz_1[i] * fe_0 + ta_yyyyy_xxxzz_0[i] * pa_x[i] - ta_yyyyy_xxxzz_1[i] * pc_x[i];

        ta_xyyyyy_xxyyy_0[i] =
            2.0 * ta_yyyyy_xyyy_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyy_1[i] * fe_0 + ta_yyyyy_xxyyy_0[i] * pa_x[i] - ta_yyyyy_xxyyy_1[i] * pc_x[i];

        ta_xyyyyy_xxyyz_0[i] =
            2.0 * ta_yyyyy_xyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyz_1[i] * fe_0 + ta_yyyyy_xxyyz_0[i] * pa_x[i] - ta_yyyyy_xxyyz_1[i] * pc_x[i];

        ta_xyyyyy_xxyzz_0[i] =
            2.0 * ta_yyyyy_xyzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyzz_1[i] * fe_0 + ta_yyyyy_xxyzz_0[i] * pa_x[i] - ta_yyyyy_xxyzz_1[i] * pc_x[i];

        ta_xyyyyy_xxzzz_0[i] =
            2.0 * ta_yyyyy_xzzz_0[i] * fe_0 - 2.0 * ta_yyyyy_xzzz_1[i] * fe_0 + ta_yyyyy_xxzzz_0[i] * pa_x[i] - ta_yyyyy_xxzzz_1[i] * pc_x[i];

        ta_xyyyyy_xyyyy_0[i] = ta_yyyyy_yyyy_0[i] * fe_0 - ta_yyyyy_yyyy_1[i] * fe_0 + ta_yyyyy_xyyyy_0[i] * pa_x[i] - ta_yyyyy_xyyyy_1[i] * pc_x[i];

        ta_xyyyyy_xyyyz_0[i] = ta_yyyyy_yyyz_0[i] * fe_0 - ta_yyyyy_yyyz_1[i] * fe_0 + ta_yyyyy_xyyyz_0[i] * pa_x[i] - ta_yyyyy_xyyyz_1[i] * pc_x[i];

        ta_xyyyyy_xyyzz_0[i] = ta_yyyyy_yyzz_0[i] * fe_0 - ta_yyyyy_yyzz_1[i] * fe_0 + ta_yyyyy_xyyzz_0[i] * pa_x[i] - ta_yyyyy_xyyzz_1[i] * pc_x[i];

        ta_xyyyyy_xyzzz_0[i] = ta_yyyyy_yzzz_0[i] * fe_0 - ta_yyyyy_yzzz_1[i] * fe_0 + ta_yyyyy_xyzzz_0[i] * pa_x[i] - ta_yyyyy_xyzzz_1[i] * pc_x[i];

        ta_xyyyyy_xzzzz_0[i] = ta_yyyyy_zzzz_0[i] * fe_0 - ta_yyyyy_zzzz_1[i] * fe_0 + ta_yyyyy_xzzzz_0[i] * pa_x[i] - ta_yyyyy_xzzzz_1[i] * pc_x[i];

        ta_xyyyyy_yyyyy_0[i] = ta_yyyyy_yyyyy_0[i] * pa_x[i] - ta_yyyyy_yyyyy_1[i] * pc_x[i];

        ta_xyyyyy_yyyyz_0[i] = ta_yyyyy_yyyyz_0[i] * pa_x[i] - ta_yyyyy_yyyyz_1[i] * pc_x[i];

        ta_xyyyyy_yyyzz_0[i] = ta_yyyyy_yyyzz_0[i] * pa_x[i] - ta_yyyyy_yyyzz_1[i] * pc_x[i];

        ta_xyyyyy_yyzzz_0[i] = ta_yyyyy_yyzzz_0[i] * pa_x[i] - ta_yyyyy_yyzzz_1[i] * pc_x[i];

        ta_xyyyyy_yzzzz_0[i] = ta_yyyyy_yzzzz_0[i] * pa_x[i] - ta_yyyyy_yzzzz_1[i] * pc_x[i];

        ta_xyyyyy_zzzzz_0[i] = ta_yyyyy_zzzzz_0[i] * pa_x[i] - ta_yyyyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 336-357 components of targeted buffer : IH

    auto ta_xyyyyz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 336);

    auto ta_xyyyyz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 337);

    auto ta_xyyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 338);

    auto ta_xyyyyz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 339);

    auto ta_xyyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 340);

    auto ta_xyyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 341);

    auto ta_xyyyyz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 342);

    auto ta_xyyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 343);

    auto ta_xyyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 344);

    auto ta_xyyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 345);

    auto ta_xyyyyz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 346);

    auto ta_xyyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 347);

    auto ta_xyyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 348);

    auto ta_xyyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 349);

    auto ta_xyyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 350);

    auto ta_xyyyyz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 351);

    auto ta_xyyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 352);

    auto ta_xyyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 353);

    auto ta_xyyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 354);

    auto ta_xyyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 355);

    auto ta_xyyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 356);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xyyyy_xxxxx_0,  \
                             ta_xyyyy_xxxxx_1,  \
                             ta_xyyyy_xxxxy_0,  \
                             ta_xyyyy_xxxxy_1,  \
                             ta_xyyyy_xxxyy_0,  \
                             ta_xyyyy_xxxyy_1,  \
                             ta_xyyyy_xxyyy_0,  \
                             ta_xyyyy_xxyyy_1,  \
                             ta_xyyyy_xyyyy_0,  \
                             ta_xyyyy_xyyyy_1,  \
                             ta_xyyyyz_xxxxx_0, \
                             ta_xyyyyz_xxxxy_0, \
                             ta_xyyyyz_xxxxz_0, \
                             ta_xyyyyz_xxxyy_0, \
                             ta_xyyyyz_xxxyz_0, \
                             ta_xyyyyz_xxxzz_0, \
                             ta_xyyyyz_xxyyy_0, \
                             ta_xyyyyz_xxyyz_0, \
                             ta_xyyyyz_xxyzz_0, \
                             ta_xyyyyz_xxzzz_0, \
                             ta_xyyyyz_xyyyy_0, \
                             ta_xyyyyz_xyyyz_0, \
                             ta_xyyyyz_xyyzz_0, \
                             ta_xyyyyz_xyzzz_0, \
                             ta_xyyyyz_xzzzz_0, \
                             ta_xyyyyz_yyyyy_0, \
                             ta_xyyyyz_yyyyz_0, \
                             ta_xyyyyz_yyyzz_0, \
                             ta_xyyyyz_yyzzz_0, \
                             ta_xyyyyz_yzzzz_0, \
                             ta_xyyyyz_zzzzz_0, \
                             ta_yyyyz_xxxxz_0,  \
                             ta_yyyyz_xxxxz_1,  \
                             ta_yyyyz_xxxyz_0,  \
                             ta_yyyyz_xxxyz_1,  \
                             ta_yyyyz_xxxz_0,   \
                             ta_yyyyz_xxxz_1,   \
                             ta_yyyyz_xxxzz_0,  \
                             ta_yyyyz_xxxzz_1,  \
                             ta_yyyyz_xxyyz_0,  \
                             ta_yyyyz_xxyyz_1,  \
                             ta_yyyyz_xxyz_0,   \
                             ta_yyyyz_xxyz_1,   \
                             ta_yyyyz_xxyzz_0,  \
                             ta_yyyyz_xxyzz_1,  \
                             ta_yyyyz_xxzz_0,   \
                             ta_yyyyz_xxzz_1,   \
                             ta_yyyyz_xxzzz_0,  \
                             ta_yyyyz_xxzzz_1,  \
                             ta_yyyyz_xyyyz_0,  \
                             ta_yyyyz_xyyyz_1,  \
                             ta_yyyyz_xyyz_0,   \
                             ta_yyyyz_xyyz_1,   \
                             ta_yyyyz_xyyzz_0,  \
                             ta_yyyyz_xyyzz_1,  \
                             ta_yyyyz_xyzz_0,   \
                             ta_yyyyz_xyzz_1,   \
                             ta_yyyyz_xyzzz_0,  \
                             ta_yyyyz_xyzzz_1,  \
                             ta_yyyyz_xzzz_0,   \
                             ta_yyyyz_xzzz_1,   \
                             ta_yyyyz_xzzzz_0,  \
                             ta_yyyyz_xzzzz_1,  \
                             ta_yyyyz_yyyyy_0,  \
                             ta_yyyyz_yyyyy_1,  \
                             ta_yyyyz_yyyyz_0,  \
                             ta_yyyyz_yyyyz_1,  \
                             ta_yyyyz_yyyz_0,   \
                             ta_yyyyz_yyyz_1,   \
                             ta_yyyyz_yyyzz_0,  \
                             ta_yyyyz_yyyzz_1,  \
                             ta_yyyyz_yyzz_0,   \
                             ta_yyyyz_yyzz_1,   \
                             ta_yyyyz_yyzzz_0,  \
                             ta_yyyyz_yyzzz_1,  \
                             ta_yyyyz_yzzz_0,   \
                             ta_yyyyz_yzzz_1,   \
                             ta_yyyyz_yzzzz_0,  \
                             ta_yyyyz_yzzzz_1,  \
                             ta_yyyyz_zzzz_0,   \
                             ta_yyyyz_zzzz_1,   \
                             ta_yyyyz_zzzzz_0,  \
                             ta_yyyyz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyz_xxxxx_0[i] = ta_xyyyy_xxxxx_0[i] * pa_z[i] - ta_xyyyy_xxxxx_1[i] * pc_z[i];

        ta_xyyyyz_xxxxy_0[i] = ta_xyyyy_xxxxy_0[i] * pa_z[i] - ta_xyyyy_xxxxy_1[i] * pc_z[i];

        ta_xyyyyz_xxxxz_0[i] =
            4.0 * ta_yyyyz_xxxz_0[i] * fe_0 - 4.0 * ta_yyyyz_xxxz_1[i] * fe_0 + ta_yyyyz_xxxxz_0[i] * pa_x[i] - ta_yyyyz_xxxxz_1[i] * pc_x[i];

        ta_xyyyyz_xxxyy_0[i] = ta_xyyyy_xxxyy_0[i] * pa_z[i] - ta_xyyyy_xxxyy_1[i] * pc_z[i];

        ta_xyyyyz_xxxyz_0[i] =
            3.0 * ta_yyyyz_xxyz_0[i] * fe_0 - 3.0 * ta_yyyyz_xxyz_1[i] * fe_0 + ta_yyyyz_xxxyz_0[i] * pa_x[i] - ta_yyyyz_xxxyz_1[i] * pc_x[i];

        ta_xyyyyz_xxxzz_0[i] =
            3.0 * ta_yyyyz_xxzz_0[i] * fe_0 - 3.0 * ta_yyyyz_xxzz_1[i] * fe_0 + ta_yyyyz_xxxzz_0[i] * pa_x[i] - ta_yyyyz_xxxzz_1[i] * pc_x[i];

        ta_xyyyyz_xxyyy_0[i] = ta_xyyyy_xxyyy_0[i] * pa_z[i] - ta_xyyyy_xxyyy_1[i] * pc_z[i];

        ta_xyyyyz_xxyyz_0[i] =
            2.0 * ta_yyyyz_xyyz_0[i] * fe_0 - 2.0 * ta_yyyyz_xyyz_1[i] * fe_0 + ta_yyyyz_xxyyz_0[i] * pa_x[i] - ta_yyyyz_xxyyz_1[i] * pc_x[i];

        ta_xyyyyz_xxyzz_0[i] =
            2.0 * ta_yyyyz_xyzz_0[i] * fe_0 - 2.0 * ta_yyyyz_xyzz_1[i] * fe_0 + ta_yyyyz_xxyzz_0[i] * pa_x[i] - ta_yyyyz_xxyzz_1[i] * pc_x[i];

        ta_xyyyyz_xxzzz_0[i] =
            2.0 * ta_yyyyz_xzzz_0[i] * fe_0 - 2.0 * ta_yyyyz_xzzz_1[i] * fe_0 + ta_yyyyz_xxzzz_0[i] * pa_x[i] - ta_yyyyz_xxzzz_1[i] * pc_x[i];

        ta_xyyyyz_xyyyy_0[i] = ta_xyyyy_xyyyy_0[i] * pa_z[i] - ta_xyyyy_xyyyy_1[i] * pc_z[i];

        ta_xyyyyz_xyyyz_0[i] = ta_yyyyz_yyyz_0[i] * fe_0 - ta_yyyyz_yyyz_1[i] * fe_0 + ta_yyyyz_xyyyz_0[i] * pa_x[i] - ta_yyyyz_xyyyz_1[i] * pc_x[i];

        ta_xyyyyz_xyyzz_0[i] = ta_yyyyz_yyzz_0[i] * fe_0 - ta_yyyyz_yyzz_1[i] * fe_0 + ta_yyyyz_xyyzz_0[i] * pa_x[i] - ta_yyyyz_xyyzz_1[i] * pc_x[i];

        ta_xyyyyz_xyzzz_0[i] = ta_yyyyz_yzzz_0[i] * fe_0 - ta_yyyyz_yzzz_1[i] * fe_0 + ta_yyyyz_xyzzz_0[i] * pa_x[i] - ta_yyyyz_xyzzz_1[i] * pc_x[i];

        ta_xyyyyz_xzzzz_0[i] = ta_yyyyz_zzzz_0[i] * fe_0 - ta_yyyyz_zzzz_1[i] * fe_0 + ta_yyyyz_xzzzz_0[i] * pa_x[i] - ta_yyyyz_xzzzz_1[i] * pc_x[i];

        ta_xyyyyz_yyyyy_0[i] = ta_yyyyz_yyyyy_0[i] * pa_x[i] - ta_yyyyz_yyyyy_1[i] * pc_x[i];

        ta_xyyyyz_yyyyz_0[i] = ta_yyyyz_yyyyz_0[i] * pa_x[i] - ta_yyyyz_yyyyz_1[i] * pc_x[i];

        ta_xyyyyz_yyyzz_0[i] = ta_yyyyz_yyyzz_0[i] * pa_x[i] - ta_yyyyz_yyyzz_1[i] * pc_x[i];

        ta_xyyyyz_yyzzz_0[i] = ta_yyyyz_yyzzz_0[i] * pa_x[i] - ta_yyyyz_yyzzz_1[i] * pc_x[i];

        ta_xyyyyz_yzzzz_0[i] = ta_yyyyz_yzzzz_0[i] * pa_x[i] - ta_yyyyz_yzzzz_1[i] * pc_x[i];

        ta_xyyyyz_zzzzz_0[i] = ta_yyyyz_zzzzz_0[i] * pa_x[i] - ta_yyyyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 357-378 components of targeted buffer : IH

    auto ta_xyyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 357);

    auto ta_xyyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 358);

    auto ta_xyyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 359);

    auto ta_xyyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 360);

    auto ta_xyyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 361);

    auto ta_xyyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 362);

    auto ta_xyyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 363);

    auto ta_xyyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 364);

    auto ta_xyyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 365);

    auto ta_xyyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 366);

    auto ta_xyyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 367);

    auto ta_xyyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 368);

    auto ta_xyyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 369);

    auto ta_xyyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 370);

    auto ta_xyyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 371);

    auto ta_xyyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 372);

    auto ta_xyyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 373);

    auto ta_xyyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 374);

    auto ta_xyyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 375);

    auto ta_xyyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 376);

    auto ta_xyyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 377);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xyyyzz_xxxxx_0, \
                             ta_xyyyzz_xxxxy_0, \
                             ta_xyyyzz_xxxxz_0, \
                             ta_xyyyzz_xxxyy_0, \
                             ta_xyyyzz_xxxyz_0, \
                             ta_xyyyzz_xxxzz_0, \
                             ta_xyyyzz_xxyyy_0, \
                             ta_xyyyzz_xxyyz_0, \
                             ta_xyyyzz_xxyzz_0, \
                             ta_xyyyzz_xxzzz_0, \
                             ta_xyyyzz_xyyyy_0, \
                             ta_xyyyzz_xyyyz_0, \
                             ta_xyyyzz_xyyzz_0, \
                             ta_xyyyzz_xyzzz_0, \
                             ta_xyyyzz_xzzzz_0, \
                             ta_xyyyzz_yyyyy_0, \
                             ta_xyyyzz_yyyyz_0, \
                             ta_xyyyzz_yyyzz_0, \
                             ta_xyyyzz_yyzzz_0, \
                             ta_xyyyzz_yzzzz_0, \
                             ta_xyyyzz_zzzzz_0, \
                             ta_yyyzz_xxxx_0,   \
                             ta_yyyzz_xxxx_1,   \
                             ta_yyyzz_xxxxx_0,  \
                             ta_yyyzz_xxxxx_1,  \
                             ta_yyyzz_xxxxy_0,  \
                             ta_yyyzz_xxxxy_1,  \
                             ta_yyyzz_xxxxz_0,  \
                             ta_yyyzz_xxxxz_1,  \
                             ta_yyyzz_xxxy_0,   \
                             ta_yyyzz_xxxy_1,   \
                             ta_yyyzz_xxxyy_0,  \
                             ta_yyyzz_xxxyy_1,  \
                             ta_yyyzz_xxxyz_0,  \
                             ta_yyyzz_xxxyz_1,  \
                             ta_yyyzz_xxxz_0,   \
                             ta_yyyzz_xxxz_1,   \
                             ta_yyyzz_xxxzz_0,  \
                             ta_yyyzz_xxxzz_1,  \
                             ta_yyyzz_xxyy_0,   \
                             ta_yyyzz_xxyy_1,   \
                             ta_yyyzz_xxyyy_0,  \
                             ta_yyyzz_xxyyy_1,  \
                             ta_yyyzz_xxyyz_0,  \
                             ta_yyyzz_xxyyz_1,  \
                             ta_yyyzz_xxyz_0,   \
                             ta_yyyzz_xxyz_1,   \
                             ta_yyyzz_xxyzz_0,  \
                             ta_yyyzz_xxyzz_1,  \
                             ta_yyyzz_xxzz_0,   \
                             ta_yyyzz_xxzz_1,   \
                             ta_yyyzz_xxzzz_0,  \
                             ta_yyyzz_xxzzz_1,  \
                             ta_yyyzz_xyyy_0,   \
                             ta_yyyzz_xyyy_1,   \
                             ta_yyyzz_xyyyy_0,  \
                             ta_yyyzz_xyyyy_1,  \
                             ta_yyyzz_xyyyz_0,  \
                             ta_yyyzz_xyyyz_1,  \
                             ta_yyyzz_xyyz_0,   \
                             ta_yyyzz_xyyz_1,   \
                             ta_yyyzz_xyyzz_0,  \
                             ta_yyyzz_xyyzz_1,  \
                             ta_yyyzz_xyzz_0,   \
                             ta_yyyzz_xyzz_1,   \
                             ta_yyyzz_xyzzz_0,  \
                             ta_yyyzz_xyzzz_1,  \
                             ta_yyyzz_xzzz_0,   \
                             ta_yyyzz_xzzz_1,   \
                             ta_yyyzz_xzzzz_0,  \
                             ta_yyyzz_xzzzz_1,  \
                             ta_yyyzz_yyyy_0,   \
                             ta_yyyzz_yyyy_1,   \
                             ta_yyyzz_yyyyy_0,  \
                             ta_yyyzz_yyyyy_1,  \
                             ta_yyyzz_yyyyz_0,  \
                             ta_yyyzz_yyyyz_1,  \
                             ta_yyyzz_yyyz_0,   \
                             ta_yyyzz_yyyz_1,   \
                             ta_yyyzz_yyyzz_0,  \
                             ta_yyyzz_yyyzz_1,  \
                             ta_yyyzz_yyzz_0,   \
                             ta_yyyzz_yyzz_1,   \
                             ta_yyyzz_yyzzz_0,  \
                             ta_yyyzz_yyzzz_1,  \
                             ta_yyyzz_yzzz_0,   \
                             ta_yyyzz_yzzz_1,   \
                             ta_yyyzz_yzzzz_0,  \
                             ta_yyyzz_yzzzz_1,  \
                             ta_yyyzz_zzzz_0,   \
                             ta_yyyzz_zzzz_1,   \
                             ta_yyyzz_zzzzz_0,  \
                             ta_yyyzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyzz_xxxxx_0[i] =
            5.0 * ta_yyyzz_xxxx_0[i] * fe_0 - 5.0 * ta_yyyzz_xxxx_1[i] * fe_0 + ta_yyyzz_xxxxx_0[i] * pa_x[i] - ta_yyyzz_xxxxx_1[i] * pc_x[i];

        ta_xyyyzz_xxxxy_0[i] =
            4.0 * ta_yyyzz_xxxy_0[i] * fe_0 - 4.0 * ta_yyyzz_xxxy_1[i] * fe_0 + ta_yyyzz_xxxxy_0[i] * pa_x[i] - ta_yyyzz_xxxxy_1[i] * pc_x[i];

        ta_xyyyzz_xxxxz_0[i] =
            4.0 * ta_yyyzz_xxxz_0[i] * fe_0 - 4.0 * ta_yyyzz_xxxz_1[i] * fe_0 + ta_yyyzz_xxxxz_0[i] * pa_x[i] - ta_yyyzz_xxxxz_1[i] * pc_x[i];

        ta_xyyyzz_xxxyy_0[i] =
            3.0 * ta_yyyzz_xxyy_0[i] * fe_0 - 3.0 * ta_yyyzz_xxyy_1[i] * fe_0 + ta_yyyzz_xxxyy_0[i] * pa_x[i] - ta_yyyzz_xxxyy_1[i] * pc_x[i];

        ta_xyyyzz_xxxyz_0[i] =
            3.0 * ta_yyyzz_xxyz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxyz_1[i] * fe_0 + ta_yyyzz_xxxyz_0[i] * pa_x[i] - ta_yyyzz_xxxyz_1[i] * pc_x[i];

        ta_xyyyzz_xxxzz_0[i] =
            3.0 * ta_yyyzz_xxzz_0[i] * fe_0 - 3.0 * ta_yyyzz_xxzz_1[i] * fe_0 + ta_yyyzz_xxxzz_0[i] * pa_x[i] - ta_yyyzz_xxxzz_1[i] * pc_x[i];

        ta_xyyyzz_xxyyy_0[i] =
            2.0 * ta_yyyzz_xyyy_0[i] * fe_0 - 2.0 * ta_yyyzz_xyyy_1[i] * fe_0 + ta_yyyzz_xxyyy_0[i] * pa_x[i] - ta_yyyzz_xxyyy_1[i] * pc_x[i];

        ta_xyyyzz_xxyyz_0[i] =
            2.0 * ta_yyyzz_xyyz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyyz_1[i] * fe_0 + ta_yyyzz_xxyyz_0[i] * pa_x[i] - ta_yyyzz_xxyyz_1[i] * pc_x[i];

        ta_xyyyzz_xxyzz_0[i] =
            2.0 * ta_yyyzz_xyzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xyzz_1[i] * fe_0 + ta_yyyzz_xxyzz_0[i] * pa_x[i] - ta_yyyzz_xxyzz_1[i] * pc_x[i];

        ta_xyyyzz_xxzzz_0[i] =
            2.0 * ta_yyyzz_xzzz_0[i] * fe_0 - 2.0 * ta_yyyzz_xzzz_1[i] * fe_0 + ta_yyyzz_xxzzz_0[i] * pa_x[i] - ta_yyyzz_xxzzz_1[i] * pc_x[i];

        ta_xyyyzz_xyyyy_0[i] = ta_yyyzz_yyyy_0[i] * fe_0 - ta_yyyzz_yyyy_1[i] * fe_0 + ta_yyyzz_xyyyy_0[i] * pa_x[i] - ta_yyyzz_xyyyy_1[i] * pc_x[i];

        ta_xyyyzz_xyyyz_0[i] = ta_yyyzz_yyyz_0[i] * fe_0 - ta_yyyzz_yyyz_1[i] * fe_0 + ta_yyyzz_xyyyz_0[i] * pa_x[i] - ta_yyyzz_xyyyz_1[i] * pc_x[i];

        ta_xyyyzz_xyyzz_0[i] = ta_yyyzz_yyzz_0[i] * fe_0 - ta_yyyzz_yyzz_1[i] * fe_0 + ta_yyyzz_xyyzz_0[i] * pa_x[i] - ta_yyyzz_xyyzz_1[i] * pc_x[i];

        ta_xyyyzz_xyzzz_0[i] = ta_yyyzz_yzzz_0[i] * fe_0 - ta_yyyzz_yzzz_1[i] * fe_0 + ta_yyyzz_xyzzz_0[i] * pa_x[i] - ta_yyyzz_xyzzz_1[i] * pc_x[i];

        ta_xyyyzz_xzzzz_0[i] = ta_yyyzz_zzzz_0[i] * fe_0 - ta_yyyzz_zzzz_1[i] * fe_0 + ta_yyyzz_xzzzz_0[i] * pa_x[i] - ta_yyyzz_xzzzz_1[i] * pc_x[i];

        ta_xyyyzz_yyyyy_0[i] = ta_yyyzz_yyyyy_0[i] * pa_x[i] - ta_yyyzz_yyyyy_1[i] * pc_x[i];

        ta_xyyyzz_yyyyz_0[i] = ta_yyyzz_yyyyz_0[i] * pa_x[i] - ta_yyyzz_yyyyz_1[i] * pc_x[i];

        ta_xyyyzz_yyyzz_0[i] = ta_yyyzz_yyyzz_0[i] * pa_x[i] - ta_yyyzz_yyyzz_1[i] * pc_x[i];

        ta_xyyyzz_yyzzz_0[i] = ta_yyyzz_yyzzz_0[i] * pa_x[i] - ta_yyyzz_yyzzz_1[i] * pc_x[i];

        ta_xyyyzz_yzzzz_0[i] = ta_yyyzz_yzzzz_0[i] * pa_x[i] - ta_yyyzz_yzzzz_1[i] * pc_x[i];

        ta_xyyyzz_zzzzz_0[i] = ta_yyyzz_zzzzz_0[i] * pa_x[i] - ta_yyyzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 378-399 components of targeted buffer : IH

    auto ta_xyyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 378);

    auto ta_xyyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 379);

    auto ta_xyyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 380);

    auto ta_xyyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 381);

    auto ta_xyyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 382);

    auto ta_xyyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 383);

    auto ta_xyyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 384);

    auto ta_xyyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 385);

    auto ta_xyyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 386);

    auto ta_xyyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 387);

    auto ta_xyyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 388);

    auto ta_xyyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 389);

    auto ta_xyyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 390);

    auto ta_xyyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 391);

    auto ta_xyyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 392);

    auto ta_xyyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 393);

    auto ta_xyyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 394);

    auto ta_xyyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 395);

    auto ta_xyyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 396);

    auto ta_xyyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 397);

    auto ta_xyyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 398);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xyyzzz_xxxxx_0, \
                             ta_xyyzzz_xxxxy_0, \
                             ta_xyyzzz_xxxxz_0, \
                             ta_xyyzzz_xxxyy_0, \
                             ta_xyyzzz_xxxyz_0, \
                             ta_xyyzzz_xxxzz_0, \
                             ta_xyyzzz_xxyyy_0, \
                             ta_xyyzzz_xxyyz_0, \
                             ta_xyyzzz_xxyzz_0, \
                             ta_xyyzzz_xxzzz_0, \
                             ta_xyyzzz_xyyyy_0, \
                             ta_xyyzzz_xyyyz_0, \
                             ta_xyyzzz_xyyzz_0, \
                             ta_xyyzzz_xyzzz_0, \
                             ta_xyyzzz_xzzzz_0, \
                             ta_xyyzzz_yyyyy_0, \
                             ta_xyyzzz_yyyyz_0, \
                             ta_xyyzzz_yyyzz_0, \
                             ta_xyyzzz_yyzzz_0, \
                             ta_xyyzzz_yzzzz_0, \
                             ta_xyyzzz_zzzzz_0, \
                             ta_yyzzz_xxxx_0,   \
                             ta_yyzzz_xxxx_1,   \
                             ta_yyzzz_xxxxx_0,  \
                             ta_yyzzz_xxxxx_1,  \
                             ta_yyzzz_xxxxy_0,  \
                             ta_yyzzz_xxxxy_1,  \
                             ta_yyzzz_xxxxz_0,  \
                             ta_yyzzz_xxxxz_1,  \
                             ta_yyzzz_xxxy_0,   \
                             ta_yyzzz_xxxy_1,   \
                             ta_yyzzz_xxxyy_0,  \
                             ta_yyzzz_xxxyy_1,  \
                             ta_yyzzz_xxxyz_0,  \
                             ta_yyzzz_xxxyz_1,  \
                             ta_yyzzz_xxxz_0,   \
                             ta_yyzzz_xxxz_1,   \
                             ta_yyzzz_xxxzz_0,  \
                             ta_yyzzz_xxxzz_1,  \
                             ta_yyzzz_xxyy_0,   \
                             ta_yyzzz_xxyy_1,   \
                             ta_yyzzz_xxyyy_0,  \
                             ta_yyzzz_xxyyy_1,  \
                             ta_yyzzz_xxyyz_0,  \
                             ta_yyzzz_xxyyz_1,  \
                             ta_yyzzz_xxyz_0,   \
                             ta_yyzzz_xxyz_1,   \
                             ta_yyzzz_xxyzz_0,  \
                             ta_yyzzz_xxyzz_1,  \
                             ta_yyzzz_xxzz_0,   \
                             ta_yyzzz_xxzz_1,   \
                             ta_yyzzz_xxzzz_0,  \
                             ta_yyzzz_xxzzz_1,  \
                             ta_yyzzz_xyyy_0,   \
                             ta_yyzzz_xyyy_1,   \
                             ta_yyzzz_xyyyy_0,  \
                             ta_yyzzz_xyyyy_1,  \
                             ta_yyzzz_xyyyz_0,  \
                             ta_yyzzz_xyyyz_1,  \
                             ta_yyzzz_xyyz_0,   \
                             ta_yyzzz_xyyz_1,   \
                             ta_yyzzz_xyyzz_0,  \
                             ta_yyzzz_xyyzz_1,  \
                             ta_yyzzz_xyzz_0,   \
                             ta_yyzzz_xyzz_1,   \
                             ta_yyzzz_xyzzz_0,  \
                             ta_yyzzz_xyzzz_1,  \
                             ta_yyzzz_xzzz_0,   \
                             ta_yyzzz_xzzz_1,   \
                             ta_yyzzz_xzzzz_0,  \
                             ta_yyzzz_xzzzz_1,  \
                             ta_yyzzz_yyyy_0,   \
                             ta_yyzzz_yyyy_1,   \
                             ta_yyzzz_yyyyy_0,  \
                             ta_yyzzz_yyyyy_1,  \
                             ta_yyzzz_yyyyz_0,  \
                             ta_yyzzz_yyyyz_1,  \
                             ta_yyzzz_yyyz_0,   \
                             ta_yyzzz_yyyz_1,   \
                             ta_yyzzz_yyyzz_0,  \
                             ta_yyzzz_yyyzz_1,  \
                             ta_yyzzz_yyzz_0,   \
                             ta_yyzzz_yyzz_1,   \
                             ta_yyzzz_yyzzz_0,  \
                             ta_yyzzz_yyzzz_1,  \
                             ta_yyzzz_yzzz_0,   \
                             ta_yyzzz_yzzz_1,   \
                             ta_yyzzz_yzzzz_0,  \
                             ta_yyzzz_yzzzz_1,  \
                             ta_yyzzz_zzzz_0,   \
                             ta_yyzzz_zzzz_1,   \
                             ta_yyzzz_zzzzz_0,  \
                             ta_yyzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzzz_xxxxx_0[i] =
            5.0 * ta_yyzzz_xxxx_0[i] * fe_0 - 5.0 * ta_yyzzz_xxxx_1[i] * fe_0 + ta_yyzzz_xxxxx_0[i] * pa_x[i] - ta_yyzzz_xxxxx_1[i] * pc_x[i];

        ta_xyyzzz_xxxxy_0[i] =
            4.0 * ta_yyzzz_xxxy_0[i] * fe_0 - 4.0 * ta_yyzzz_xxxy_1[i] * fe_0 + ta_yyzzz_xxxxy_0[i] * pa_x[i] - ta_yyzzz_xxxxy_1[i] * pc_x[i];

        ta_xyyzzz_xxxxz_0[i] =
            4.0 * ta_yyzzz_xxxz_0[i] * fe_0 - 4.0 * ta_yyzzz_xxxz_1[i] * fe_0 + ta_yyzzz_xxxxz_0[i] * pa_x[i] - ta_yyzzz_xxxxz_1[i] * pc_x[i];

        ta_xyyzzz_xxxyy_0[i] =
            3.0 * ta_yyzzz_xxyy_0[i] * fe_0 - 3.0 * ta_yyzzz_xxyy_1[i] * fe_0 + ta_yyzzz_xxxyy_0[i] * pa_x[i] - ta_yyzzz_xxxyy_1[i] * pc_x[i];

        ta_xyyzzz_xxxyz_0[i] =
            3.0 * ta_yyzzz_xxyz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxyz_1[i] * fe_0 + ta_yyzzz_xxxyz_0[i] * pa_x[i] - ta_yyzzz_xxxyz_1[i] * pc_x[i];

        ta_xyyzzz_xxxzz_0[i] =
            3.0 * ta_yyzzz_xxzz_0[i] * fe_0 - 3.0 * ta_yyzzz_xxzz_1[i] * fe_0 + ta_yyzzz_xxxzz_0[i] * pa_x[i] - ta_yyzzz_xxxzz_1[i] * pc_x[i];

        ta_xyyzzz_xxyyy_0[i] =
            2.0 * ta_yyzzz_xyyy_0[i] * fe_0 - 2.0 * ta_yyzzz_xyyy_1[i] * fe_0 + ta_yyzzz_xxyyy_0[i] * pa_x[i] - ta_yyzzz_xxyyy_1[i] * pc_x[i];

        ta_xyyzzz_xxyyz_0[i] =
            2.0 * ta_yyzzz_xyyz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyyz_1[i] * fe_0 + ta_yyzzz_xxyyz_0[i] * pa_x[i] - ta_yyzzz_xxyyz_1[i] * pc_x[i];

        ta_xyyzzz_xxyzz_0[i] =
            2.0 * ta_yyzzz_xyzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xyzz_1[i] * fe_0 + ta_yyzzz_xxyzz_0[i] * pa_x[i] - ta_yyzzz_xxyzz_1[i] * pc_x[i];

        ta_xyyzzz_xxzzz_0[i] =
            2.0 * ta_yyzzz_xzzz_0[i] * fe_0 - 2.0 * ta_yyzzz_xzzz_1[i] * fe_0 + ta_yyzzz_xxzzz_0[i] * pa_x[i] - ta_yyzzz_xxzzz_1[i] * pc_x[i];

        ta_xyyzzz_xyyyy_0[i] = ta_yyzzz_yyyy_0[i] * fe_0 - ta_yyzzz_yyyy_1[i] * fe_0 + ta_yyzzz_xyyyy_0[i] * pa_x[i] - ta_yyzzz_xyyyy_1[i] * pc_x[i];

        ta_xyyzzz_xyyyz_0[i] = ta_yyzzz_yyyz_0[i] * fe_0 - ta_yyzzz_yyyz_1[i] * fe_0 + ta_yyzzz_xyyyz_0[i] * pa_x[i] - ta_yyzzz_xyyyz_1[i] * pc_x[i];

        ta_xyyzzz_xyyzz_0[i] = ta_yyzzz_yyzz_0[i] * fe_0 - ta_yyzzz_yyzz_1[i] * fe_0 + ta_yyzzz_xyyzz_0[i] * pa_x[i] - ta_yyzzz_xyyzz_1[i] * pc_x[i];

        ta_xyyzzz_xyzzz_0[i] = ta_yyzzz_yzzz_0[i] * fe_0 - ta_yyzzz_yzzz_1[i] * fe_0 + ta_yyzzz_xyzzz_0[i] * pa_x[i] - ta_yyzzz_xyzzz_1[i] * pc_x[i];

        ta_xyyzzz_xzzzz_0[i] = ta_yyzzz_zzzz_0[i] * fe_0 - ta_yyzzz_zzzz_1[i] * fe_0 + ta_yyzzz_xzzzz_0[i] * pa_x[i] - ta_yyzzz_xzzzz_1[i] * pc_x[i];

        ta_xyyzzz_yyyyy_0[i] = ta_yyzzz_yyyyy_0[i] * pa_x[i] - ta_yyzzz_yyyyy_1[i] * pc_x[i];

        ta_xyyzzz_yyyyz_0[i] = ta_yyzzz_yyyyz_0[i] * pa_x[i] - ta_yyzzz_yyyyz_1[i] * pc_x[i];

        ta_xyyzzz_yyyzz_0[i] = ta_yyzzz_yyyzz_0[i] * pa_x[i] - ta_yyzzz_yyyzz_1[i] * pc_x[i];

        ta_xyyzzz_yyzzz_0[i] = ta_yyzzz_yyzzz_0[i] * pa_x[i] - ta_yyzzz_yyzzz_1[i] * pc_x[i];

        ta_xyyzzz_yzzzz_0[i] = ta_yyzzz_yzzzz_0[i] * pa_x[i] - ta_yyzzz_yzzzz_1[i] * pc_x[i];

        ta_xyyzzz_zzzzz_0[i] = ta_yyzzz_zzzzz_0[i] * pa_x[i] - ta_yyzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 399-420 components of targeted buffer : IH

    auto ta_xyzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 399);

    auto ta_xyzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 400);

    auto ta_xyzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 401);

    auto ta_xyzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 402);

    auto ta_xyzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 403);

    auto ta_xyzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 404);

    auto ta_xyzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 405);

    auto ta_xyzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 406);

    auto ta_xyzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 407);

    auto ta_xyzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 408);

    auto ta_xyzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 409);

    auto ta_xyzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 410);

    auto ta_xyzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 411);

    auto ta_xyzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 412);

    auto ta_xyzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 413);

    auto ta_xyzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 414);

    auto ta_xyzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 415);

    auto ta_xyzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 416);

    auto ta_xyzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 417);

    auto ta_xyzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 418);

    auto ta_xyzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 419);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xyzzzz_xxxxx_0, \
                             ta_xyzzzz_xxxxy_0, \
                             ta_xyzzzz_xxxxz_0, \
                             ta_xyzzzz_xxxyy_0, \
                             ta_xyzzzz_xxxyz_0, \
                             ta_xyzzzz_xxxzz_0, \
                             ta_xyzzzz_xxyyy_0, \
                             ta_xyzzzz_xxyyz_0, \
                             ta_xyzzzz_xxyzz_0, \
                             ta_xyzzzz_xxzzz_0, \
                             ta_xyzzzz_xyyyy_0, \
                             ta_xyzzzz_xyyyz_0, \
                             ta_xyzzzz_xyyzz_0, \
                             ta_xyzzzz_xyzzz_0, \
                             ta_xyzzzz_xzzzz_0, \
                             ta_xyzzzz_yyyyy_0, \
                             ta_xyzzzz_yyyyz_0, \
                             ta_xyzzzz_yyyzz_0, \
                             ta_xyzzzz_yyzzz_0, \
                             ta_xyzzzz_yzzzz_0, \
                             ta_xyzzzz_zzzzz_0, \
                             ta_xzzzz_xxxxx_0,  \
                             ta_xzzzz_xxxxx_1,  \
                             ta_xzzzz_xxxxz_0,  \
                             ta_xzzzz_xxxxz_1,  \
                             ta_xzzzz_xxxzz_0,  \
                             ta_xzzzz_xxxzz_1,  \
                             ta_xzzzz_xxzzz_0,  \
                             ta_xzzzz_xxzzz_1,  \
                             ta_xzzzz_xzzzz_0,  \
                             ta_xzzzz_xzzzz_1,  \
                             ta_yzzzz_xxxxy_0,  \
                             ta_yzzzz_xxxxy_1,  \
                             ta_yzzzz_xxxy_0,   \
                             ta_yzzzz_xxxy_1,   \
                             ta_yzzzz_xxxyy_0,  \
                             ta_yzzzz_xxxyy_1,  \
                             ta_yzzzz_xxxyz_0,  \
                             ta_yzzzz_xxxyz_1,  \
                             ta_yzzzz_xxyy_0,   \
                             ta_yzzzz_xxyy_1,   \
                             ta_yzzzz_xxyyy_0,  \
                             ta_yzzzz_xxyyy_1,  \
                             ta_yzzzz_xxyyz_0,  \
                             ta_yzzzz_xxyyz_1,  \
                             ta_yzzzz_xxyz_0,   \
                             ta_yzzzz_xxyz_1,   \
                             ta_yzzzz_xxyzz_0,  \
                             ta_yzzzz_xxyzz_1,  \
                             ta_yzzzz_xyyy_0,   \
                             ta_yzzzz_xyyy_1,   \
                             ta_yzzzz_xyyyy_0,  \
                             ta_yzzzz_xyyyy_1,  \
                             ta_yzzzz_xyyyz_0,  \
                             ta_yzzzz_xyyyz_1,  \
                             ta_yzzzz_xyyz_0,   \
                             ta_yzzzz_xyyz_1,   \
                             ta_yzzzz_xyyzz_0,  \
                             ta_yzzzz_xyyzz_1,  \
                             ta_yzzzz_xyzz_0,   \
                             ta_yzzzz_xyzz_1,   \
                             ta_yzzzz_xyzzz_0,  \
                             ta_yzzzz_xyzzz_1,  \
                             ta_yzzzz_yyyy_0,   \
                             ta_yzzzz_yyyy_1,   \
                             ta_yzzzz_yyyyy_0,  \
                             ta_yzzzz_yyyyy_1,  \
                             ta_yzzzz_yyyyz_0,  \
                             ta_yzzzz_yyyyz_1,  \
                             ta_yzzzz_yyyz_0,   \
                             ta_yzzzz_yyyz_1,   \
                             ta_yzzzz_yyyzz_0,  \
                             ta_yzzzz_yyyzz_1,  \
                             ta_yzzzz_yyzz_0,   \
                             ta_yzzzz_yyzz_1,   \
                             ta_yzzzz_yyzzz_0,  \
                             ta_yzzzz_yyzzz_1,  \
                             ta_yzzzz_yzzz_0,   \
                             ta_yzzzz_yzzz_1,   \
                             ta_yzzzz_yzzzz_0,  \
                             ta_yzzzz_yzzzz_1,  \
                             ta_yzzzz_zzzzz_0,  \
                             ta_yzzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzzz_xxxxx_0[i] = ta_xzzzz_xxxxx_0[i] * pa_y[i] - ta_xzzzz_xxxxx_1[i] * pc_y[i];

        ta_xyzzzz_xxxxy_0[i] =
            4.0 * ta_yzzzz_xxxy_0[i] * fe_0 - 4.0 * ta_yzzzz_xxxy_1[i] * fe_0 + ta_yzzzz_xxxxy_0[i] * pa_x[i] - ta_yzzzz_xxxxy_1[i] * pc_x[i];

        ta_xyzzzz_xxxxz_0[i] = ta_xzzzz_xxxxz_0[i] * pa_y[i] - ta_xzzzz_xxxxz_1[i] * pc_y[i];

        ta_xyzzzz_xxxyy_0[i] =
            3.0 * ta_yzzzz_xxyy_0[i] * fe_0 - 3.0 * ta_yzzzz_xxyy_1[i] * fe_0 + ta_yzzzz_xxxyy_0[i] * pa_x[i] - ta_yzzzz_xxxyy_1[i] * pc_x[i];

        ta_xyzzzz_xxxyz_0[i] =
            3.0 * ta_yzzzz_xxyz_0[i] * fe_0 - 3.0 * ta_yzzzz_xxyz_1[i] * fe_0 + ta_yzzzz_xxxyz_0[i] * pa_x[i] - ta_yzzzz_xxxyz_1[i] * pc_x[i];

        ta_xyzzzz_xxxzz_0[i] = ta_xzzzz_xxxzz_0[i] * pa_y[i] - ta_xzzzz_xxxzz_1[i] * pc_y[i];

        ta_xyzzzz_xxyyy_0[i] =
            2.0 * ta_yzzzz_xyyy_0[i] * fe_0 - 2.0 * ta_yzzzz_xyyy_1[i] * fe_0 + ta_yzzzz_xxyyy_0[i] * pa_x[i] - ta_yzzzz_xxyyy_1[i] * pc_x[i];

        ta_xyzzzz_xxyyz_0[i] =
            2.0 * ta_yzzzz_xyyz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyyz_1[i] * fe_0 + ta_yzzzz_xxyyz_0[i] * pa_x[i] - ta_yzzzz_xxyyz_1[i] * pc_x[i];

        ta_xyzzzz_xxyzz_0[i] =
            2.0 * ta_yzzzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzzzz_xyzz_1[i] * fe_0 + ta_yzzzz_xxyzz_0[i] * pa_x[i] - ta_yzzzz_xxyzz_1[i] * pc_x[i];

        ta_xyzzzz_xxzzz_0[i] = ta_xzzzz_xxzzz_0[i] * pa_y[i] - ta_xzzzz_xxzzz_1[i] * pc_y[i];

        ta_xyzzzz_xyyyy_0[i] = ta_yzzzz_yyyy_0[i] * fe_0 - ta_yzzzz_yyyy_1[i] * fe_0 + ta_yzzzz_xyyyy_0[i] * pa_x[i] - ta_yzzzz_xyyyy_1[i] * pc_x[i];

        ta_xyzzzz_xyyyz_0[i] = ta_yzzzz_yyyz_0[i] * fe_0 - ta_yzzzz_yyyz_1[i] * fe_0 + ta_yzzzz_xyyyz_0[i] * pa_x[i] - ta_yzzzz_xyyyz_1[i] * pc_x[i];

        ta_xyzzzz_xyyzz_0[i] = ta_yzzzz_yyzz_0[i] * fe_0 - ta_yzzzz_yyzz_1[i] * fe_0 + ta_yzzzz_xyyzz_0[i] * pa_x[i] - ta_yzzzz_xyyzz_1[i] * pc_x[i];

        ta_xyzzzz_xyzzz_0[i] = ta_yzzzz_yzzz_0[i] * fe_0 - ta_yzzzz_yzzz_1[i] * fe_0 + ta_yzzzz_xyzzz_0[i] * pa_x[i] - ta_yzzzz_xyzzz_1[i] * pc_x[i];

        ta_xyzzzz_xzzzz_0[i] = ta_xzzzz_xzzzz_0[i] * pa_y[i] - ta_xzzzz_xzzzz_1[i] * pc_y[i];

        ta_xyzzzz_yyyyy_0[i] = ta_yzzzz_yyyyy_0[i] * pa_x[i] - ta_yzzzz_yyyyy_1[i] * pc_x[i];

        ta_xyzzzz_yyyyz_0[i] = ta_yzzzz_yyyyz_0[i] * pa_x[i] - ta_yzzzz_yyyyz_1[i] * pc_x[i];

        ta_xyzzzz_yyyzz_0[i] = ta_yzzzz_yyyzz_0[i] * pa_x[i] - ta_yzzzz_yyyzz_1[i] * pc_x[i];

        ta_xyzzzz_yyzzz_0[i] = ta_yzzzz_yyzzz_0[i] * pa_x[i] - ta_yzzzz_yyzzz_1[i] * pc_x[i];

        ta_xyzzzz_yzzzz_0[i] = ta_yzzzz_yzzzz_0[i] * pa_x[i] - ta_yzzzz_yzzzz_1[i] * pc_x[i];

        ta_xyzzzz_zzzzz_0[i] = ta_yzzzz_zzzzz_0[i] * pa_x[i] - ta_yzzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 420-441 components of targeted buffer : IH

    auto ta_xzzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 420);

    auto ta_xzzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 421);

    auto ta_xzzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 422);

    auto ta_xzzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 423);

    auto ta_xzzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 424);

    auto ta_xzzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 425);

    auto ta_xzzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 426);

    auto ta_xzzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 427);

    auto ta_xzzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 428);

    auto ta_xzzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 429);

    auto ta_xzzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 430);

    auto ta_xzzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 431);

    auto ta_xzzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 432);

    auto ta_xzzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 433);

    auto ta_xzzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 434);

    auto ta_xzzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 435);

    auto ta_xzzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 436);

    auto ta_xzzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 437);

    auto ta_xzzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 438);

    auto ta_xzzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 439);

    auto ta_xzzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 440);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xzzzzz_xxxxx_0, \
                             ta_xzzzzz_xxxxy_0, \
                             ta_xzzzzz_xxxxz_0, \
                             ta_xzzzzz_xxxyy_0, \
                             ta_xzzzzz_xxxyz_0, \
                             ta_xzzzzz_xxxzz_0, \
                             ta_xzzzzz_xxyyy_0, \
                             ta_xzzzzz_xxyyz_0, \
                             ta_xzzzzz_xxyzz_0, \
                             ta_xzzzzz_xxzzz_0, \
                             ta_xzzzzz_xyyyy_0, \
                             ta_xzzzzz_xyyyz_0, \
                             ta_xzzzzz_xyyzz_0, \
                             ta_xzzzzz_xyzzz_0, \
                             ta_xzzzzz_xzzzz_0, \
                             ta_xzzzzz_yyyyy_0, \
                             ta_xzzzzz_yyyyz_0, \
                             ta_xzzzzz_yyyzz_0, \
                             ta_xzzzzz_yyzzz_0, \
                             ta_xzzzzz_yzzzz_0, \
                             ta_xzzzzz_zzzzz_0, \
                             ta_zzzzz_xxxx_0,   \
                             ta_zzzzz_xxxx_1,   \
                             ta_zzzzz_xxxxx_0,  \
                             ta_zzzzz_xxxxx_1,  \
                             ta_zzzzz_xxxxy_0,  \
                             ta_zzzzz_xxxxy_1,  \
                             ta_zzzzz_xxxxz_0,  \
                             ta_zzzzz_xxxxz_1,  \
                             ta_zzzzz_xxxy_0,   \
                             ta_zzzzz_xxxy_1,   \
                             ta_zzzzz_xxxyy_0,  \
                             ta_zzzzz_xxxyy_1,  \
                             ta_zzzzz_xxxyz_0,  \
                             ta_zzzzz_xxxyz_1,  \
                             ta_zzzzz_xxxz_0,   \
                             ta_zzzzz_xxxz_1,   \
                             ta_zzzzz_xxxzz_0,  \
                             ta_zzzzz_xxxzz_1,  \
                             ta_zzzzz_xxyy_0,   \
                             ta_zzzzz_xxyy_1,   \
                             ta_zzzzz_xxyyy_0,  \
                             ta_zzzzz_xxyyy_1,  \
                             ta_zzzzz_xxyyz_0,  \
                             ta_zzzzz_xxyyz_1,  \
                             ta_zzzzz_xxyz_0,   \
                             ta_zzzzz_xxyz_1,   \
                             ta_zzzzz_xxyzz_0,  \
                             ta_zzzzz_xxyzz_1,  \
                             ta_zzzzz_xxzz_0,   \
                             ta_zzzzz_xxzz_1,   \
                             ta_zzzzz_xxzzz_0,  \
                             ta_zzzzz_xxzzz_1,  \
                             ta_zzzzz_xyyy_0,   \
                             ta_zzzzz_xyyy_1,   \
                             ta_zzzzz_xyyyy_0,  \
                             ta_zzzzz_xyyyy_1,  \
                             ta_zzzzz_xyyyz_0,  \
                             ta_zzzzz_xyyyz_1,  \
                             ta_zzzzz_xyyz_0,   \
                             ta_zzzzz_xyyz_1,   \
                             ta_zzzzz_xyyzz_0,  \
                             ta_zzzzz_xyyzz_1,  \
                             ta_zzzzz_xyzz_0,   \
                             ta_zzzzz_xyzz_1,   \
                             ta_zzzzz_xyzzz_0,  \
                             ta_zzzzz_xyzzz_1,  \
                             ta_zzzzz_xzzz_0,   \
                             ta_zzzzz_xzzz_1,   \
                             ta_zzzzz_xzzzz_0,  \
                             ta_zzzzz_xzzzz_1,  \
                             ta_zzzzz_yyyy_0,   \
                             ta_zzzzz_yyyy_1,   \
                             ta_zzzzz_yyyyy_0,  \
                             ta_zzzzz_yyyyy_1,  \
                             ta_zzzzz_yyyyz_0,  \
                             ta_zzzzz_yyyyz_1,  \
                             ta_zzzzz_yyyz_0,   \
                             ta_zzzzz_yyyz_1,   \
                             ta_zzzzz_yyyzz_0,  \
                             ta_zzzzz_yyyzz_1,  \
                             ta_zzzzz_yyzz_0,   \
                             ta_zzzzz_yyzz_1,   \
                             ta_zzzzz_yyzzz_0,  \
                             ta_zzzzz_yyzzz_1,  \
                             ta_zzzzz_yzzz_0,   \
                             ta_zzzzz_yzzz_1,   \
                             ta_zzzzz_yzzzz_0,  \
                             ta_zzzzz_yzzzz_1,  \
                             ta_zzzzz_zzzz_0,   \
                             ta_zzzzz_zzzz_1,   \
                             ta_zzzzz_zzzzz_0,  \
                             ta_zzzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzzz_xxxxx_0[i] =
            5.0 * ta_zzzzz_xxxx_0[i] * fe_0 - 5.0 * ta_zzzzz_xxxx_1[i] * fe_0 + ta_zzzzz_xxxxx_0[i] * pa_x[i] - ta_zzzzz_xxxxx_1[i] * pc_x[i];

        ta_xzzzzz_xxxxy_0[i] =
            4.0 * ta_zzzzz_xxxy_0[i] * fe_0 - 4.0 * ta_zzzzz_xxxy_1[i] * fe_0 + ta_zzzzz_xxxxy_0[i] * pa_x[i] - ta_zzzzz_xxxxy_1[i] * pc_x[i];

        ta_xzzzzz_xxxxz_0[i] =
            4.0 * ta_zzzzz_xxxz_0[i] * fe_0 - 4.0 * ta_zzzzz_xxxz_1[i] * fe_0 + ta_zzzzz_xxxxz_0[i] * pa_x[i] - ta_zzzzz_xxxxz_1[i] * pc_x[i];

        ta_xzzzzz_xxxyy_0[i] =
            3.0 * ta_zzzzz_xxyy_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyy_1[i] * fe_0 + ta_zzzzz_xxxyy_0[i] * pa_x[i] - ta_zzzzz_xxxyy_1[i] * pc_x[i];

        ta_xzzzzz_xxxyz_0[i] =
            3.0 * ta_zzzzz_xxyz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyz_1[i] * fe_0 + ta_zzzzz_xxxyz_0[i] * pa_x[i] - ta_zzzzz_xxxyz_1[i] * pc_x[i];

        ta_xzzzzz_xxxzz_0[i] =
            3.0 * ta_zzzzz_xxzz_0[i] * fe_0 - 3.0 * ta_zzzzz_xxzz_1[i] * fe_0 + ta_zzzzz_xxxzz_0[i] * pa_x[i] - ta_zzzzz_xxxzz_1[i] * pc_x[i];

        ta_xzzzzz_xxyyy_0[i] =
            2.0 * ta_zzzzz_xyyy_0[i] * fe_0 - 2.0 * ta_zzzzz_xyyy_1[i] * fe_0 + ta_zzzzz_xxyyy_0[i] * pa_x[i] - ta_zzzzz_xxyyy_1[i] * pc_x[i];

        ta_xzzzzz_xxyyz_0[i] =
            2.0 * ta_zzzzz_xyyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyyz_1[i] * fe_0 + ta_zzzzz_xxyyz_0[i] * pa_x[i] - ta_zzzzz_xxyyz_1[i] * pc_x[i];

        ta_xzzzzz_xxyzz_0[i] =
            2.0 * ta_zzzzz_xyzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyzz_1[i] * fe_0 + ta_zzzzz_xxyzz_0[i] * pa_x[i] - ta_zzzzz_xxyzz_1[i] * pc_x[i];

        ta_xzzzzz_xxzzz_0[i] =
            2.0 * ta_zzzzz_xzzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xzzz_1[i] * fe_0 + ta_zzzzz_xxzzz_0[i] * pa_x[i] - ta_zzzzz_xxzzz_1[i] * pc_x[i];

        ta_xzzzzz_xyyyy_0[i] = ta_zzzzz_yyyy_0[i] * fe_0 - ta_zzzzz_yyyy_1[i] * fe_0 + ta_zzzzz_xyyyy_0[i] * pa_x[i] - ta_zzzzz_xyyyy_1[i] * pc_x[i];

        ta_xzzzzz_xyyyz_0[i] = ta_zzzzz_yyyz_0[i] * fe_0 - ta_zzzzz_yyyz_1[i] * fe_0 + ta_zzzzz_xyyyz_0[i] * pa_x[i] - ta_zzzzz_xyyyz_1[i] * pc_x[i];

        ta_xzzzzz_xyyzz_0[i] = ta_zzzzz_yyzz_0[i] * fe_0 - ta_zzzzz_yyzz_1[i] * fe_0 + ta_zzzzz_xyyzz_0[i] * pa_x[i] - ta_zzzzz_xyyzz_1[i] * pc_x[i];

        ta_xzzzzz_xyzzz_0[i] = ta_zzzzz_yzzz_0[i] * fe_0 - ta_zzzzz_yzzz_1[i] * fe_0 + ta_zzzzz_xyzzz_0[i] * pa_x[i] - ta_zzzzz_xyzzz_1[i] * pc_x[i];

        ta_xzzzzz_xzzzz_0[i] = ta_zzzzz_zzzz_0[i] * fe_0 - ta_zzzzz_zzzz_1[i] * fe_0 + ta_zzzzz_xzzzz_0[i] * pa_x[i] - ta_zzzzz_xzzzz_1[i] * pc_x[i];

        ta_xzzzzz_yyyyy_0[i] = ta_zzzzz_yyyyy_0[i] * pa_x[i] - ta_zzzzz_yyyyy_1[i] * pc_x[i];

        ta_xzzzzz_yyyyz_0[i] = ta_zzzzz_yyyyz_0[i] * pa_x[i] - ta_zzzzz_yyyyz_1[i] * pc_x[i];

        ta_xzzzzz_yyyzz_0[i] = ta_zzzzz_yyyzz_0[i] * pa_x[i] - ta_zzzzz_yyyzz_1[i] * pc_x[i];

        ta_xzzzzz_yyzzz_0[i] = ta_zzzzz_yyzzz_0[i] * pa_x[i] - ta_zzzzz_yyzzz_1[i] * pc_x[i];

        ta_xzzzzz_yzzzz_0[i] = ta_zzzzz_yzzzz_0[i] * pa_x[i] - ta_zzzzz_yzzzz_1[i] * pc_x[i];

        ta_xzzzzz_zzzzz_0[i] = ta_zzzzz_zzzzz_0[i] * pa_x[i] - ta_zzzzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 441-462 components of targeted buffer : IH

    auto ta_yyyyyy_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 441);

    auto ta_yyyyyy_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 442);

    auto ta_yyyyyy_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 443);

    auto ta_yyyyyy_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 444);

    auto ta_yyyyyy_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 445);

    auto ta_yyyyyy_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 446);

    auto ta_yyyyyy_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 447);

    auto ta_yyyyyy_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 448);

    auto ta_yyyyyy_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 449);

    auto ta_yyyyyy_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 450);

    auto ta_yyyyyy_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 451);

    auto ta_yyyyyy_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 452);

    auto ta_yyyyyy_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 453);

    auto ta_yyyyyy_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 454);

    auto ta_yyyyyy_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 455);

    auto ta_yyyyyy_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 456);

    auto ta_yyyyyy_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 457);

    auto ta_yyyyyy_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 458);

    auto ta_yyyyyy_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 459);

    auto ta_yyyyyy_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 460);

    auto ta_yyyyyy_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 461);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta_yyyy_xxxxx_0,   \
                             ta_yyyy_xxxxx_1,   \
                             ta_yyyy_xxxxy_0,   \
                             ta_yyyy_xxxxy_1,   \
                             ta_yyyy_xxxxz_0,   \
                             ta_yyyy_xxxxz_1,   \
                             ta_yyyy_xxxyy_0,   \
                             ta_yyyy_xxxyy_1,   \
                             ta_yyyy_xxxyz_0,   \
                             ta_yyyy_xxxyz_1,   \
                             ta_yyyy_xxxzz_0,   \
                             ta_yyyy_xxxzz_1,   \
                             ta_yyyy_xxyyy_0,   \
                             ta_yyyy_xxyyy_1,   \
                             ta_yyyy_xxyyz_0,   \
                             ta_yyyy_xxyyz_1,   \
                             ta_yyyy_xxyzz_0,   \
                             ta_yyyy_xxyzz_1,   \
                             ta_yyyy_xxzzz_0,   \
                             ta_yyyy_xxzzz_1,   \
                             ta_yyyy_xyyyy_0,   \
                             ta_yyyy_xyyyy_1,   \
                             ta_yyyy_xyyyz_0,   \
                             ta_yyyy_xyyyz_1,   \
                             ta_yyyy_xyyzz_0,   \
                             ta_yyyy_xyyzz_1,   \
                             ta_yyyy_xyzzz_0,   \
                             ta_yyyy_xyzzz_1,   \
                             ta_yyyy_xzzzz_0,   \
                             ta_yyyy_xzzzz_1,   \
                             ta_yyyy_yyyyy_0,   \
                             ta_yyyy_yyyyy_1,   \
                             ta_yyyy_yyyyz_0,   \
                             ta_yyyy_yyyyz_1,   \
                             ta_yyyy_yyyzz_0,   \
                             ta_yyyy_yyyzz_1,   \
                             ta_yyyy_yyzzz_0,   \
                             ta_yyyy_yyzzz_1,   \
                             ta_yyyy_yzzzz_0,   \
                             ta_yyyy_yzzzz_1,   \
                             ta_yyyy_zzzzz_0,   \
                             ta_yyyy_zzzzz_1,   \
                             ta_yyyyy_xxxx_0,   \
                             ta_yyyyy_xxxx_1,   \
                             ta_yyyyy_xxxxx_0,  \
                             ta_yyyyy_xxxxx_1,  \
                             ta_yyyyy_xxxxy_0,  \
                             ta_yyyyy_xxxxy_1,  \
                             ta_yyyyy_xxxxz_0,  \
                             ta_yyyyy_xxxxz_1,  \
                             ta_yyyyy_xxxy_0,   \
                             ta_yyyyy_xxxy_1,   \
                             ta_yyyyy_xxxyy_0,  \
                             ta_yyyyy_xxxyy_1,  \
                             ta_yyyyy_xxxyz_0,  \
                             ta_yyyyy_xxxyz_1,  \
                             ta_yyyyy_xxxz_0,   \
                             ta_yyyyy_xxxz_1,   \
                             ta_yyyyy_xxxzz_0,  \
                             ta_yyyyy_xxxzz_1,  \
                             ta_yyyyy_xxyy_0,   \
                             ta_yyyyy_xxyy_1,   \
                             ta_yyyyy_xxyyy_0,  \
                             ta_yyyyy_xxyyy_1,  \
                             ta_yyyyy_xxyyz_0,  \
                             ta_yyyyy_xxyyz_1,  \
                             ta_yyyyy_xxyz_0,   \
                             ta_yyyyy_xxyz_1,   \
                             ta_yyyyy_xxyzz_0,  \
                             ta_yyyyy_xxyzz_1,  \
                             ta_yyyyy_xxzz_0,   \
                             ta_yyyyy_xxzz_1,   \
                             ta_yyyyy_xxzzz_0,  \
                             ta_yyyyy_xxzzz_1,  \
                             ta_yyyyy_xyyy_0,   \
                             ta_yyyyy_xyyy_1,   \
                             ta_yyyyy_xyyyy_0,  \
                             ta_yyyyy_xyyyy_1,  \
                             ta_yyyyy_xyyyz_0,  \
                             ta_yyyyy_xyyyz_1,  \
                             ta_yyyyy_xyyz_0,   \
                             ta_yyyyy_xyyz_1,   \
                             ta_yyyyy_xyyzz_0,  \
                             ta_yyyyy_xyyzz_1,  \
                             ta_yyyyy_xyzz_0,   \
                             ta_yyyyy_xyzz_1,   \
                             ta_yyyyy_xyzzz_0,  \
                             ta_yyyyy_xyzzz_1,  \
                             ta_yyyyy_xzzz_0,   \
                             ta_yyyyy_xzzz_1,   \
                             ta_yyyyy_xzzzz_0,  \
                             ta_yyyyy_xzzzz_1,  \
                             ta_yyyyy_yyyy_0,   \
                             ta_yyyyy_yyyy_1,   \
                             ta_yyyyy_yyyyy_0,  \
                             ta_yyyyy_yyyyy_1,  \
                             ta_yyyyy_yyyyz_0,  \
                             ta_yyyyy_yyyyz_1,  \
                             ta_yyyyy_yyyz_0,   \
                             ta_yyyyy_yyyz_1,   \
                             ta_yyyyy_yyyzz_0,  \
                             ta_yyyyy_yyyzz_1,  \
                             ta_yyyyy_yyzz_0,   \
                             ta_yyyyy_yyzz_1,   \
                             ta_yyyyy_yyzzz_0,  \
                             ta_yyyyy_yyzzz_1,  \
                             ta_yyyyy_yzzz_0,   \
                             ta_yyyyy_yzzz_1,   \
                             ta_yyyyy_yzzzz_0,  \
                             ta_yyyyy_yzzzz_1,  \
                             ta_yyyyy_zzzz_0,   \
                             ta_yyyyy_zzzz_1,   \
                             ta_yyyyy_zzzzz_0,  \
                             ta_yyyyy_zzzzz_1,  \
                             ta_yyyyyy_xxxxx_0, \
                             ta_yyyyyy_xxxxy_0, \
                             ta_yyyyyy_xxxxz_0, \
                             ta_yyyyyy_xxxyy_0, \
                             ta_yyyyyy_xxxyz_0, \
                             ta_yyyyyy_xxxzz_0, \
                             ta_yyyyyy_xxyyy_0, \
                             ta_yyyyyy_xxyyz_0, \
                             ta_yyyyyy_xxyzz_0, \
                             ta_yyyyyy_xxzzz_0, \
                             ta_yyyyyy_xyyyy_0, \
                             ta_yyyyyy_xyyyz_0, \
                             ta_yyyyyy_xyyzz_0, \
                             ta_yyyyyy_xyzzz_0, \
                             ta_yyyyyy_xzzzz_0, \
                             ta_yyyyyy_yyyyy_0, \
                             ta_yyyyyy_yyyyz_0, \
                             ta_yyyyyy_yyyzz_0, \
                             ta_yyyyyy_yyzzz_0, \
                             ta_yyyyyy_yzzzz_0, \
                             ta_yyyyyy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyy_xxxxx_0[i] =
            5.0 * ta_yyyy_xxxxx_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxx_1[i] * fe_0 + ta_yyyyy_xxxxx_0[i] * pa_y[i] - ta_yyyyy_xxxxx_1[i] * pc_y[i];

        ta_yyyyyy_xxxxy_0[i] = 5.0 * ta_yyyy_xxxxy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxy_1[i] * fe_0 + ta_yyyyy_xxxx_0[i] * fe_0 -
                               ta_yyyyy_xxxx_1[i] * fe_0 + ta_yyyyy_xxxxy_0[i] * pa_y[i] - ta_yyyyy_xxxxy_1[i] * pc_y[i];

        ta_yyyyyy_xxxxz_0[i] =
            5.0 * ta_yyyy_xxxxz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxz_1[i] * fe_0 + ta_yyyyy_xxxxz_0[i] * pa_y[i] - ta_yyyyy_xxxxz_1[i] * pc_y[i];

        ta_yyyyyy_xxxyy_0[i] = 5.0 * ta_yyyy_xxxyy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxyy_1[i] * fe_0 + 2.0 * ta_yyyyy_xxxy_0[i] * fe_0 -
                               2.0 * ta_yyyyy_xxxy_1[i] * fe_0 + ta_yyyyy_xxxyy_0[i] * pa_y[i] - ta_yyyyy_xxxyy_1[i] * pc_y[i];

        ta_yyyyyy_xxxyz_0[i] = 5.0 * ta_yyyy_xxxyz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxyz_1[i] * fe_0 + ta_yyyyy_xxxz_0[i] * fe_0 -
                               ta_yyyyy_xxxz_1[i] * fe_0 + ta_yyyyy_xxxyz_0[i] * pa_y[i] - ta_yyyyy_xxxyz_1[i] * pc_y[i];

        ta_yyyyyy_xxxzz_0[i] =
            5.0 * ta_yyyy_xxxzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxzz_1[i] * fe_0 + ta_yyyyy_xxxzz_0[i] * pa_y[i] - ta_yyyyy_xxxzz_1[i] * pc_y[i];

        ta_yyyyyy_xxyyy_0[i] = 5.0 * ta_yyyy_xxyyy_0[i] * fe_0 - 5.0 * ta_yyyy_xxyyy_1[i] * fe_0 + 3.0 * ta_yyyyy_xxyy_0[i] * fe_0 -
                               3.0 * ta_yyyyy_xxyy_1[i] * fe_0 + ta_yyyyy_xxyyy_0[i] * pa_y[i] - ta_yyyyy_xxyyy_1[i] * pc_y[i];

        ta_yyyyyy_xxyyz_0[i] = 5.0 * ta_yyyy_xxyyz_0[i] * fe_0 - 5.0 * ta_yyyy_xxyyz_1[i] * fe_0 + 2.0 * ta_yyyyy_xxyz_0[i] * fe_0 -
                               2.0 * ta_yyyyy_xxyz_1[i] * fe_0 + ta_yyyyy_xxyyz_0[i] * pa_y[i] - ta_yyyyy_xxyyz_1[i] * pc_y[i];

        ta_yyyyyy_xxyzz_0[i] = 5.0 * ta_yyyy_xxyzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxyzz_1[i] * fe_0 + ta_yyyyy_xxzz_0[i] * fe_0 -
                               ta_yyyyy_xxzz_1[i] * fe_0 + ta_yyyyy_xxyzz_0[i] * pa_y[i] - ta_yyyyy_xxyzz_1[i] * pc_y[i];

        ta_yyyyyy_xxzzz_0[i] =
            5.0 * ta_yyyy_xxzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xxzzz_1[i] * fe_0 + ta_yyyyy_xxzzz_0[i] * pa_y[i] - ta_yyyyy_xxzzz_1[i] * pc_y[i];

        ta_yyyyyy_xyyyy_0[i] = 5.0 * ta_yyyy_xyyyy_0[i] * fe_0 - 5.0 * ta_yyyy_xyyyy_1[i] * fe_0 + 4.0 * ta_yyyyy_xyyy_0[i] * fe_0 -
                               4.0 * ta_yyyyy_xyyy_1[i] * fe_0 + ta_yyyyy_xyyyy_0[i] * pa_y[i] - ta_yyyyy_xyyyy_1[i] * pc_y[i];

        ta_yyyyyy_xyyyz_0[i] = 5.0 * ta_yyyy_xyyyz_0[i] * fe_0 - 5.0 * ta_yyyy_xyyyz_1[i] * fe_0 + 3.0 * ta_yyyyy_xyyz_0[i] * fe_0 -
                               3.0 * ta_yyyyy_xyyz_1[i] * fe_0 + ta_yyyyy_xyyyz_0[i] * pa_y[i] - ta_yyyyy_xyyyz_1[i] * pc_y[i];

        ta_yyyyyy_xyyzz_0[i] = 5.0 * ta_yyyy_xyyzz_0[i] * fe_0 - 5.0 * ta_yyyy_xyyzz_1[i] * fe_0 + 2.0 * ta_yyyyy_xyzz_0[i] * fe_0 -
                               2.0 * ta_yyyyy_xyzz_1[i] * fe_0 + ta_yyyyy_xyyzz_0[i] * pa_y[i] - ta_yyyyy_xyyzz_1[i] * pc_y[i];

        ta_yyyyyy_xyzzz_0[i] = 5.0 * ta_yyyy_xyzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xyzzz_1[i] * fe_0 + ta_yyyyy_xzzz_0[i] * fe_0 -
                               ta_yyyyy_xzzz_1[i] * fe_0 + ta_yyyyy_xyzzz_0[i] * pa_y[i] - ta_yyyyy_xyzzz_1[i] * pc_y[i];

        ta_yyyyyy_xzzzz_0[i] =
            5.0 * ta_yyyy_xzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_xzzzz_1[i] * fe_0 + ta_yyyyy_xzzzz_0[i] * pa_y[i] - ta_yyyyy_xzzzz_1[i] * pc_y[i];

        ta_yyyyyy_yyyyy_0[i] = 5.0 * ta_yyyy_yyyyy_0[i] * fe_0 - 5.0 * ta_yyyy_yyyyy_1[i] * fe_0 + 5.0 * ta_yyyyy_yyyy_0[i] * fe_0 -
                               5.0 * ta_yyyyy_yyyy_1[i] * fe_0 + ta_yyyyy_yyyyy_0[i] * pa_y[i] - ta_yyyyy_yyyyy_1[i] * pc_y[i];

        ta_yyyyyy_yyyyz_0[i] = 5.0 * ta_yyyy_yyyyz_0[i] * fe_0 - 5.0 * ta_yyyy_yyyyz_1[i] * fe_0 + 4.0 * ta_yyyyy_yyyz_0[i] * fe_0 -
                               4.0 * ta_yyyyy_yyyz_1[i] * fe_0 + ta_yyyyy_yyyyz_0[i] * pa_y[i] - ta_yyyyy_yyyyz_1[i] * pc_y[i];

        ta_yyyyyy_yyyzz_0[i] = 5.0 * ta_yyyy_yyyzz_0[i] * fe_0 - 5.0 * ta_yyyy_yyyzz_1[i] * fe_0 + 3.0 * ta_yyyyy_yyzz_0[i] * fe_0 -
                               3.0 * ta_yyyyy_yyzz_1[i] * fe_0 + ta_yyyyy_yyyzz_0[i] * pa_y[i] - ta_yyyyy_yyyzz_1[i] * pc_y[i];

        ta_yyyyyy_yyzzz_0[i] = 5.0 * ta_yyyy_yyzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yyzzz_1[i] * fe_0 + 2.0 * ta_yyyyy_yzzz_0[i] * fe_0 -
                               2.0 * ta_yyyyy_yzzz_1[i] * fe_0 + ta_yyyyy_yyzzz_0[i] * pa_y[i] - ta_yyyyy_yyzzz_1[i] * pc_y[i];

        ta_yyyyyy_yzzzz_0[i] = 5.0 * ta_yyyy_yzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yzzzz_1[i] * fe_0 + ta_yyyyy_zzzz_0[i] * fe_0 -
                               ta_yyyyy_zzzz_1[i] * fe_0 + ta_yyyyy_yzzzz_0[i] * pa_y[i] - ta_yyyyy_yzzzz_1[i] * pc_y[i];

        ta_yyyyyy_zzzzz_0[i] =
            5.0 * ta_yyyy_zzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_zzzzz_1[i] * fe_0 + ta_yyyyy_zzzzz_0[i] * pa_y[i] - ta_yyyyy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 462-483 components of targeted buffer : IH

    auto ta_yyyyyz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 462);

    auto ta_yyyyyz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 463);

    auto ta_yyyyyz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 464);

    auto ta_yyyyyz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 465);

    auto ta_yyyyyz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 466);

    auto ta_yyyyyz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 467);

    auto ta_yyyyyz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 468);

    auto ta_yyyyyz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 469);

    auto ta_yyyyyz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 470);

    auto ta_yyyyyz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 471);

    auto ta_yyyyyz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 472);

    auto ta_yyyyyz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 473);

    auto ta_yyyyyz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 474);

    auto ta_yyyyyz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 475);

    auto ta_yyyyyz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 476);

    auto ta_yyyyyz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 477);

    auto ta_yyyyyz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 478);

    auto ta_yyyyyz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 479);

    auto ta_yyyyyz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 480);

    auto ta_yyyyyz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 481);

    auto ta_yyyyyz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 482);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyyyy_xxxxx_0,  \
                             ta_yyyyy_xxxxx_1,  \
                             ta_yyyyy_xxxxy_0,  \
                             ta_yyyyy_xxxxy_1,  \
                             ta_yyyyy_xxxy_0,   \
                             ta_yyyyy_xxxy_1,   \
                             ta_yyyyy_xxxyy_0,  \
                             ta_yyyyy_xxxyy_1,  \
                             ta_yyyyy_xxxyz_0,  \
                             ta_yyyyy_xxxyz_1,  \
                             ta_yyyyy_xxyy_0,   \
                             ta_yyyyy_xxyy_1,   \
                             ta_yyyyy_xxyyy_0,  \
                             ta_yyyyy_xxyyy_1,  \
                             ta_yyyyy_xxyyz_0,  \
                             ta_yyyyy_xxyyz_1,  \
                             ta_yyyyy_xxyz_0,   \
                             ta_yyyyy_xxyz_1,   \
                             ta_yyyyy_xxyzz_0,  \
                             ta_yyyyy_xxyzz_1,  \
                             ta_yyyyy_xyyy_0,   \
                             ta_yyyyy_xyyy_1,   \
                             ta_yyyyy_xyyyy_0,  \
                             ta_yyyyy_xyyyy_1,  \
                             ta_yyyyy_xyyyz_0,  \
                             ta_yyyyy_xyyyz_1,  \
                             ta_yyyyy_xyyz_0,   \
                             ta_yyyyy_xyyz_1,   \
                             ta_yyyyy_xyyzz_0,  \
                             ta_yyyyy_xyyzz_1,  \
                             ta_yyyyy_xyzz_0,   \
                             ta_yyyyy_xyzz_1,   \
                             ta_yyyyy_xyzzz_0,  \
                             ta_yyyyy_xyzzz_1,  \
                             ta_yyyyy_yyyy_0,   \
                             ta_yyyyy_yyyy_1,   \
                             ta_yyyyy_yyyyy_0,  \
                             ta_yyyyy_yyyyy_1,  \
                             ta_yyyyy_yyyyz_0,  \
                             ta_yyyyy_yyyyz_1,  \
                             ta_yyyyy_yyyz_0,   \
                             ta_yyyyy_yyyz_1,   \
                             ta_yyyyy_yyyzz_0,  \
                             ta_yyyyy_yyyzz_1,  \
                             ta_yyyyy_yyzz_0,   \
                             ta_yyyyy_yyzz_1,   \
                             ta_yyyyy_yyzzz_0,  \
                             ta_yyyyy_yyzzz_1,  \
                             ta_yyyyy_yzzz_0,   \
                             ta_yyyyy_yzzz_1,   \
                             ta_yyyyy_yzzzz_0,  \
                             ta_yyyyy_yzzzz_1,  \
                             ta_yyyyyz_xxxxx_0, \
                             ta_yyyyyz_xxxxy_0, \
                             ta_yyyyyz_xxxxz_0, \
                             ta_yyyyyz_xxxyy_0, \
                             ta_yyyyyz_xxxyz_0, \
                             ta_yyyyyz_xxxzz_0, \
                             ta_yyyyyz_xxyyy_0, \
                             ta_yyyyyz_xxyyz_0, \
                             ta_yyyyyz_xxyzz_0, \
                             ta_yyyyyz_xxzzz_0, \
                             ta_yyyyyz_xyyyy_0, \
                             ta_yyyyyz_xyyyz_0, \
                             ta_yyyyyz_xyyzz_0, \
                             ta_yyyyyz_xyzzz_0, \
                             ta_yyyyyz_xzzzz_0, \
                             ta_yyyyyz_yyyyy_0, \
                             ta_yyyyyz_yyyyz_0, \
                             ta_yyyyyz_yyyzz_0, \
                             ta_yyyyyz_yyzzz_0, \
                             ta_yyyyyz_yzzzz_0, \
                             ta_yyyyyz_zzzzz_0, \
                             ta_yyyyz_xxxxz_0,  \
                             ta_yyyyz_xxxxz_1,  \
                             ta_yyyyz_xxxzz_0,  \
                             ta_yyyyz_xxxzz_1,  \
                             ta_yyyyz_xxzzz_0,  \
                             ta_yyyyz_xxzzz_1,  \
                             ta_yyyyz_xzzzz_0,  \
                             ta_yyyyz_xzzzz_1,  \
                             ta_yyyyz_zzzzz_0,  \
                             ta_yyyyz_zzzzz_1,  \
                             ta_yyyz_xxxxz_0,   \
                             ta_yyyz_xxxxz_1,   \
                             ta_yyyz_xxxzz_0,   \
                             ta_yyyz_xxxzz_1,   \
                             ta_yyyz_xxzzz_0,   \
                             ta_yyyz_xxzzz_1,   \
                             ta_yyyz_xzzzz_0,   \
                             ta_yyyz_xzzzz_1,   \
                             ta_yyyz_zzzzz_0,   \
                             ta_yyyz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyz_xxxxx_0[i] = ta_yyyyy_xxxxx_0[i] * pa_z[i] - ta_yyyyy_xxxxx_1[i] * pc_z[i];

        ta_yyyyyz_xxxxy_0[i] = ta_yyyyy_xxxxy_0[i] * pa_z[i] - ta_yyyyy_xxxxy_1[i] * pc_z[i];

        ta_yyyyyz_xxxxz_0[i] =
            4.0 * ta_yyyz_xxxxz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxxz_1[i] * fe_0 + ta_yyyyz_xxxxz_0[i] * pa_y[i] - ta_yyyyz_xxxxz_1[i] * pc_y[i];

        ta_yyyyyz_xxxyy_0[i] = ta_yyyyy_xxxyy_0[i] * pa_z[i] - ta_yyyyy_xxxyy_1[i] * pc_z[i];

        ta_yyyyyz_xxxyz_0[i] = ta_yyyyy_xxxy_0[i] * fe_0 - ta_yyyyy_xxxy_1[i] * fe_0 + ta_yyyyy_xxxyz_0[i] * pa_z[i] - ta_yyyyy_xxxyz_1[i] * pc_z[i];

        ta_yyyyyz_xxxzz_0[i] =
            4.0 * ta_yyyz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxzz_1[i] * fe_0 + ta_yyyyz_xxxzz_0[i] * pa_y[i] - ta_yyyyz_xxxzz_1[i] * pc_y[i];

        ta_yyyyyz_xxyyy_0[i] = ta_yyyyy_xxyyy_0[i] * pa_z[i] - ta_yyyyy_xxyyy_1[i] * pc_z[i];

        ta_yyyyyz_xxyyz_0[i] = ta_yyyyy_xxyy_0[i] * fe_0 - ta_yyyyy_xxyy_1[i] * fe_0 + ta_yyyyy_xxyyz_0[i] * pa_z[i] - ta_yyyyy_xxyyz_1[i] * pc_z[i];

        ta_yyyyyz_xxyzz_0[i] =
            2.0 * ta_yyyyy_xxyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xxyz_1[i] * fe_0 + ta_yyyyy_xxyzz_0[i] * pa_z[i] - ta_yyyyy_xxyzz_1[i] * pc_z[i];

        ta_yyyyyz_xxzzz_0[i] =
            4.0 * ta_yyyz_xxzzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxzzz_1[i] * fe_0 + ta_yyyyz_xxzzz_0[i] * pa_y[i] - ta_yyyyz_xxzzz_1[i] * pc_y[i];

        ta_yyyyyz_xyyyy_0[i] = ta_yyyyy_xyyyy_0[i] * pa_z[i] - ta_yyyyy_xyyyy_1[i] * pc_z[i];

        ta_yyyyyz_xyyyz_0[i] = ta_yyyyy_xyyy_0[i] * fe_0 - ta_yyyyy_xyyy_1[i] * fe_0 + ta_yyyyy_xyyyz_0[i] * pa_z[i] - ta_yyyyy_xyyyz_1[i] * pc_z[i];

        ta_yyyyyz_xyyzz_0[i] =
            2.0 * ta_yyyyy_xyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_xyyz_1[i] * fe_0 + ta_yyyyy_xyyzz_0[i] * pa_z[i] - ta_yyyyy_xyyzz_1[i] * pc_z[i];

        ta_yyyyyz_xyzzz_0[i] =
            3.0 * ta_yyyyy_xyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_xyzz_1[i] * fe_0 + ta_yyyyy_xyzzz_0[i] * pa_z[i] - ta_yyyyy_xyzzz_1[i] * pc_z[i];

        ta_yyyyyz_xzzzz_0[i] =
            4.0 * ta_yyyz_xzzzz_0[i] * fe_0 - 4.0 * ta_yyyz_xzzzz_1[i] * fe_0 + ta_yyyyz_xzzzz_0[i] * pa_y[i] - ta_yyyyz_xzzzz_1[i] * pc_y[i];

        ta_yyyyyz_yyyyy_0[i] = ta_yyyyy_yyyyy_0[i] * pa_z[i] - ta_yyyyy_yyyyy_1[i] * pc_z[i];

        ta_yyyyyz_yyyyz_0[i] = ta_yyyyy_yyyy_0[i] * fe_0 - ta_yyyyy_yyyy_1[i] * fe_0 + ta_yyyyy_yyyyz_0[i] * pa_z[i] - ta_yyyyy_yyyyz_1[i] * pc_z[i];

        ta_yyyyyz_yyyzz_0[i] =
            2.0 * ta_yyyyy_yyyz_0[i] * fe_0 - 2.0 * ta_yyyyy_yyyz_1[i] * fe_0 + ta_yyyyy_yyyzz_0[i] * pa_z[i] - ta_yyyyy_yyyzz_1[i] * pc_z[i];

        ta_yyyyyz_yyzzz_0[i] =
            3.0 * ta_yyyyy_yyzz_0[i] * fe_0 - 3.0 * ta_yyyyy_yyzz_1[i] * fe_0 + ta_yyyyy_yyzzz_0[i] * pa_z[i] - ta_yyyyy_yyzzz_1[i] * pc_z[i];

        ta_yyyyyz_yzzzz_0[i] =
            4.0 * ta_yyyyy_yzzz_0[i] * fe_0 - 4.0 * ta_yyyyy_yzzz_1[i] * fe_0 + ta_yyyyy_yzzzz_0[i] * pa_z[i] - ta_yyyyy_yzzzz_1[i] * pc_z[i];

        ta_yyyyyz_zzzzz_0[i] =
            4.0 * ta_yyyz_zzzzz_0[i] * fe_0 - 4.0 * ta_yyyz_zzzzz_1[i] * fe_0 + ta_yyyyz_zzzzz_0[i] * pa_y[i] - ta_yyyyz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 483-504 components of targeted buffer : IH

    auto ta_yyyyzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 483);

    auto ta_yyyyzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 484);

    auto ta_yyyyzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 485);

    auto ta_yyyyzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 486);

    auto ta_yyyyzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 487);

    auto ta_yyyyzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 488);

    auto ta_yyyyzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 489);

    auto ta_yyyyzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 490);

    auto ta_yyyyzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 491);

    auto ta_yyyyzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 492);

    auto ta_yyyyzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 493);

    auto ta_yyyyzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 494);

    auto ta_yyyyzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 495);

    auto ta_yyyyzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 496);

    auto ta_yyyyzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 497);

    auto ta_yyyyzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 498);

    auto ta_yyyyzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 499);

    auto ta_yyyyzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 500);

    auto ta_yyyyzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 501);

    auto ta_yyyyzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 502);

    auto ta_yyyyzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 503);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyyy_xxxxy_0,   \
                             ta_yyyy_xxxxy_1,   \
                             ta_yyyy_xxxyy_0,   \
                             ta_yyyy_xxxyy_1,   \
                             ta_yyyy_xxyyy_0,   \
                             ta_yyyy_xxyyy_1,   \
                             ta_yyyy_xyyyy_0,   \
                             ta_yyyy_xyyyy_1,   \
                             ta_yyyy_yyyyy_0,   \
                             ta_yyyy_yyyyy_1,   \
                             ta_yyyyz_xxxxy_0,  \
                             ta_yyyyz_xxxxy_1,  \
                             ta_yyyyz_xxxyy_0,  \
                             ta_yyyyz_xxxyy_1,  \
                             ta_yyyyz_xxyyy_0,  \
                             ta_yyyyz_xxyyy_1,  \
                             ta_yyyyz_xyyyy_0,  \
                             ta_yyyyz_xyyyy_1,  \
                             ta_yyyyz_yyyyy_0,  \
                             ta_yyyyz_yyyyy_1,  \
                             ta_yyyyzz_xxxxx_0, \
                             ta_yyyyzz_xxxxy_0, \
                             ta_yyyyzz_xxxxz_0, \
                             ta_yyyyzz_xxxyy_0, \
                             ta_yyyyzz_xxxyz_0, \
                             ta_yyyyzz_xxxzz_0, \
                             ta_yyyyzz_xxyyy_0, \
                             ta_yyyyzz_xxyyz_0, \
                             ta_yyyyzz_xxyzz_0, \
                             ta_yyyyzz_xxzzz_0, \
                             ta_yyyyzz_xyyyy_0, \
                             ta_yyyyzz_xyyyz_0, \
                             ta_yyyyzz_xyyzz_0, \
                             ta_yyyyzz_xyzzz_0, \
                             ta_yyyyzz_xzzzz_0, \
                             ta_yyyyzz_yyyyy_0, \
                             ta_yyyyzz_yyyyz_0, \
                             ta_yyyyzz_yyyzz_0, \
                             ta_yyyyzz_yyzzz_0, \
                             ta_yyyyzz_yzzzz_0, \
                             ta_yyyyzz_zzzzz_0, \
                             ta_yyyzz_xxxxx_0,  \
                             ta_yyyzz_xxxxx_1,  \
                             ta_yyyzz_xxxxz_0,  \
                             ta_yyyzz_xxxxz_1,  \
                             ta_yyyzz_xxxyz_0,  \
                             ta_yyyzz_xxxyz_1,  \
                             ta_yyyzz_xxxz_0,   \
                             ta_yyyzz_xxxz_1,   \
                             ta_yyyzz_xxxzz_0,  \
                             ta_yyyzz_xxxzz_1,  \
                             ta_yyyzz_xxyyz_0,  \
                             ta_yyyzz_xxyyz_1,  \
                             ta_yyyzz_xxyz_0,   \
                             ta_yyyzz_xxyz_1,   \
                             ta_yyyzz_xxyzz_0,  \
                             ta_yyyzz_xxyzz_1,  \
                             ta_yyyzz_xxzz_0,   \
                             ta_yyyzz_xxzz_1,   \
                             ta_yyyzz_xxzzz_0,  \
                             ta_yyyzz_xxzzz_1,  \
                             ta_yyyzz_xyyyz_0,  \
                             ta_yyyzz_xyyyz_1,  \
                             ta_yyyzz_xyyz_0,   \
                             ta_yyyzz_xyyz_1,   \
                             ta_yyyzz_xyyzz_0,  \
                             ta_yyyzz_xyyzz_1,  \
                             ta_yyyzz_xyzz_0,   \
                             ta_yyyzz_xyzz_1,   \
                             ta_yyyzz_xyzzz_0,  \
                             ta_yyyzz_xyzzz_1,  \
                             ta_yyyzz_xzzz_0,   \
                             ta_yyyzz_xzzz_1,   \
                             ta_yyyzz_xzzzz_0,  \
                             ta_yyyzz_xzzzz_1,  \
                             ta_yyyzz_yyyyz_0,  \
                             ta_yyyzz_yyyyz_1,  \
                             ta_yyyzz_yyyz_0,   \
                             ta_yyyzz_yyyz_1,   \
                             ta_yyyzz_yyyzz_0,  \
                             ta_yyyzz_yyyzz_1,  \
                             ta_yyyzz_yyzz_0,   \
                             ta_yyyzz_yyzz_1,   \
                             ta_yyyzz_yyzzz_0,  \
                             ta_yyyzz_yyzzz_1,  \
                             ta_yyyzz_yzzz_0,   \
                             ta_yyyzz_yzzz_1,   \
                             ta_yyyzz_yzzzz_0,  \
                             ta_yyyzz_yzzzz_1,  \
                             ta_yyyzz_zzzz_0,   \
                             ta_yyyzz_zzzz_1,   \
                             ta_yyyzz_zzzzz_0,  \
                             ta_yyyzz_zzzzz_1,  \
                             ta_yyzz_xxxxx_0,   \
                             ta_yyzz_xxxxx_1,   \
                             ta_yyzz_xxxxz_0,   \
                             ta_yyzz_xxxxz_1,   \
                             ta_yyzz_xxxyz_0,   \
                             ta_yyzz_xxxyz_1,   \
                             ta_yyzz_xxxzz_0,   \
                             ta_yyzz_xxxzz_1,   \
                             ta_yyzz_xxyyz_0,   \
                             ta_yyzz_xxyyz_1,   \
                             ta_yyzz_xxyzz_0,   \
                             ta_yyzz_xxyzz_1,   \
                             ta_yyzz_xxzzz_0,   \
                             ta_yyzz_xxzzz_1,   \
                             ta_yyzz_xyyyz_0,   \
                             ta_yyzz_xyyyz_1,   \
                             ta_yyzz_xyyzz_0,   \
                             ta_yyzz_xyyzz_1,   \
                             ta_yyzz_xyzzz_0,   \
                             ta_yyzz_xyzzz_1,   \
                             ta_yyzz_xzzzz_0,   \
                             ta_yyzz_xzzzz_1,   \
                             ta_yyzz_yyyyz_0,   \
                             ta_yyzz_yyyyz_1,   \
                             ta_yyzz_yyyzz_0,   \
                             ta_yyzz_yyyzz_1,   \
                             ta_yyzz_yyzzz_0,   \
                             ta_yyzz_yyzzz_1,   \
                             ta_yyzz_yzzzz_0,   \
                             ta_yyzz_yzzzz_1,   \
                             ta_yyzz_zzzzz_0,   \
                             ta_yyzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyzz_xxxxx_0[i] =
            3.0 * ta_yyzz_xxxxx_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxx_1[i] * fe_0 + ta_yyyzz_xxxxx_0[i] * pa_y[i] - ta_yyyzz_xxxxx_1[i] * pc_y[i];

        ta_yyyyzz_xxxxy_0[i] = ta_yyyy_xxxxy_0[i] * fe_0 - ta_yyyy_xxxxy_1[i] * fe_0 + ta_yyyyz_xxxxy_0[i] * pa_z[i] - ta_yyyyz_xxxxy_1[i] * pc_z[i];

        ta_yyyyzz_xxxxz_0[i] =
            3.0 * ta_yyzz_xxxxz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxz_1[i] * fe_0 + ta_yyyzz_xxxxz_0[i] * pa_y[i] - ta_yyyzz_xxxxz_1[i] * pc_y[i];

        ta_yyyyzz_xxxyy_0[i] = ta_yyyy_xxxyy_0[i] * fe_0 - ta_yyyy_xxxyy_1[i] * fe_0 + ta_yyyyz_xxxyy_0[i] * pa_z[i] - ta_yyyyz_xxxyy_1[i] * pc_z[i];

        ta_yyyyzz_xxxyz_0[i] = 3.0 * ta_yyzz_xxxyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxyz_1[i] * fe_0 + ta_yyyzz_xxxz_0[i] * fe_0 -
                               ta_yyyzz_xxxz_1[i] * fe_0 + ta_yyyzz_xxxyz_0[i] * pa_y[i] - ta_yyyzz_xxxyz_1[i] * pc_y[i];

        ta_yyyyzz_xxxzz_0[i] =
            3.0 * ta_yyzz_xxxzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxxzz_1[i] * fe_0 + ta_yyyzz_xxxzz_0[i] * pa_y[i] - ta_yyyzz_xxxzz_1[i] * pc_y[i];

        ta_yyyyzz_xxyyy_0[i] = ta_yyyy_xxyyy_0[i] * fe_0 - ta_yyyy_xxyyy_1[i] * fe_0 + ta_yyyyz_xxyyy_0[i] * pa_z[i] - ta_yyyyz_xxyyy_1[i] * pc_z[i];

        ta_yyyyzz_xxyyz_0[i] = 3.0 * ta_yyzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyz_1[i] * fe_0 + 2.0 * ta_yyyzz_xxyz_0[i] * fe_0 -
                               2.0 * ta_yyyzz_xxyz_1[i] * fe_0 + ta_yyyzz_xxyyz_0[i] * pa_y[i] - ta_yyyzz_xxyyz_1[i] * pc_y[i];

        ta_yyyyzz_xxyzz_0[i] = 3.0 * ta_yyzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyzz_1[i] * fe_0 + ta_yyyzz_xxzz_0[i] * fe_0 -
                               ta_yyyzz_xxzz_1[i] * fe_0 + ta_yyyzz_xxyzz_0[i] * pa_y[i] - ta_yyyzz_xxyzz_1[i] * pc_y[i];

        ta_yyyyzz_xxzzz_0[i] =
            3.0 * ta_yyzz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxzzz_1[i] * fe_0 + ta_yyyzz_xxzzz_0[i] * pa_y[i] - ta_yyyzz_xxzzz_1[i] * pc_y[i];

        ta_yyyyzz_xyyyy_0[i] = ta_yyyy_xyyyy_0[i] * fe_0 - ta_yyyy_xyyyy_1[i] * fe_0 + ta_yyyyz_xyyyy_0[i] * pa_z[i] - ta_yyyyz_xyyyy_1[i] * pc_z[i];

        ta_yyyyzz_xyyyz_0[i] = 3.0 * ta_yyzz_xyyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xyyyz_1[i] * fe_0 + 3.0 * ta_yyyzz_xyyz_0[i] * fe_0 -
                               3.0 * ta_yyyzz_xyyz_1[i] * fe_0 + ta_yyyzz_xyyyz_0[i] * pa_y[i] - ta_yyyzz_xyyyz_1[i] * pc_y[i];

        ta_yyyyzz_xyyzz_0[i] = 3.0 * ta_yyzz_xyyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xyyzz_1[i] * fe_0 + 2.0 * ta_yyyzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_yyyzz_xyzz_1[i] * fe_0 + ta_yyyzz_xyyzz_0[i] * pa_y[i] - ta_yyyzz_xyyzz_1[i] * pc_y[i];

        ta_yyyyzz_xyzzz_0[i] = 3.0 * ta_yyzz_xyzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xyzzz_1[i] * fe_0 + ta_yyyzz_xzzz_0[i] * fe_0 -
                               ta_yyyzz_xzzz_1[i] * fe_0 + ta_yyyzz_xyzzz_0[i] * pa_y[i] - ta_yyyzz_xyzzz_1[i] * pc_y[i];

        ta_yyyyzz_xzzzz_0[i] =
            3.0 * ta_yyzz_xzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xzzzz_1[i] * fe_0 + ta_yyyzz_xzzzz_0[i] * pa_y[i] - ta_yyyzz_xzzzz_1[i] * pc_y[i];

        ta_yyyyzz_yyyyy_0[i] = ta_yyyy_yyyyy_0[i] * fe_0 - ta_yyyy_yyyyy_1[i] * fe_0 + ta_yyyyz_yyyyy_0[i] * pa_z[i] - ta_yyyyz_yyyyy_1[i] * pc_z[i];

        ta_yyyyzz_yyyyz_0[i] = 3.0 * ta_yyzz_yyyyz_0[i] * fe_0 - 3.0 * ta_yyzz_yyyyz_1[i] * fe_0 + 4.0 * ta_yyyzz_yyyz_0[i] * fe_0 -
                               4.0 * ta_yyyzz_yyyz_1[i] * fe_0 + ta_yyyzz_yyyyz_0[i] * pa_y[i] - ta_yyyzz_yyyyz_1[i] * pc_y[i];

        ta_yyyyzz_yyyzz_0[i] = 3.0 * ta_yyzz_yyyzz_0[i] * fe_0 - 3.0 * ta_yyzz_yyyzz_1[i] * fe_0 + 3.0 * ta_yyyzz_yyzz_0[i] * fe_0 -
                               3.0 * ta_yyyzz_yyzz_1[i] * fe_0 + ta_yyyzz_yyyzz_0[i] * pa_y[i] - ta_yyyzz_yyyzz_1[i] * pc_y[i];

        ta_yyyyzz_yyzzz_0[i] = 3.0 * ta_yyzz_yyzzz_0[i] * fe_0 - 3.0 * ta_yyzz_yyzzz_1[i] * fe_0 + 2.0 * ta_yyyzz_yzzz_0[i] * fe_0 -
                               2.0 * ta_yyyzz_yzzz_1[i] * fe_0 + ta_yyyzz_yyzzz_0[i] * pa_y[i] - ta_yyyzz_yyzzz_1[i] * pc_y[i];

        ta_yyyyzz_yzzzz_0[i] = 3.0 * ta_yyzz_yzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_yzzzz_1[i] * fe_0 + ta_yyyzz_zzzz_0[i] * fe_0 -
                               ta_yyyzz_zzzz_1[i] * fe_0 + ta_yyyzz_yzzzz_0[i] * pa_y[i] - ta_yyyzz_yzzzz_1[i] * pc_y[i];

        ta_yyyyzz_zzzzz_0[i] =
            3.0 * ta_yyzz_zzzzz_0[i] * fe_0 - 3.0 * ta_yyzz_zzzzz_1[i] * fe_0 + ta_yyyzz_zzzzz_0[i] * pa_y[i] - ta_yyyzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 504-525 components of targeted buffer : IH

    auto ta_yyyzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 504);

    auto ta_yyyzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 505);

    auto ta_yyyzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 506);

    auto ta_yyyzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 507);

    auto ta_yyyzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 508);

    auto ta_yyyzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 509);

    auto ta_yyyzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 510);

    auto ta_yyyzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 511);

    auto ta_yyyzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 512);

    auto ta_yyyzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 513);

    auto ta_yyyzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 514);

    auto ta_yyyzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 515);

    auto ta_yyyzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 516);

    auto ta_yyyzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 517);

    auto ta_yyyzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 518);

    auto ta_yyyzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 519);

    auto ta_yyyzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 520);

    auto ta_yyyzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 521);

    auto ta_yyyzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 522);

    auto ta_yyyzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 523);

    auto ta_yyyzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 524);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyyz_xxxxy_0,   \
                             ta_yyyz_xxxxy_1,   \
                             ta_yyyz_xxxyy_0,   \
                             ta_yyyz_xxxyy_1,   \
                             ta_yyyz_xxyyy_0,   \
                             ta_yyyz_xxyyy_1,   \
                             ta_yyyz_xyyyy_0,   \
                             ta_yyyz_xyyyy_1,   \
                             ta_yyyz_yyyyy_0,   \
                             ta_yyyz_yyyyy_1,   \
                             ta_yyyzz_xxxxy_0,  \
                             ta_yyyzz_xxxxy_1,  \
                             ta_yyyzz_xxxyy_0,  \
                             ta_yyyzz_xxxyy_1,  \
                             ta_yyyzz_xxyyy_0,  \
                             ta_yyyzz_xxyyy_1,  \
                             ta_yyyzz_xyyyy_0,  \
                             ta_yyyzz_xyyyy_1,  \
                             ta_yyyzz_yyyyy_0,  \
                             ta_yyyzz_yyyyy_1,  \
                             ta_yyyzzz_xxxxx_0, \
                             ta_yyyzzz_xxxxy_0, \
                             ta_yyyzzz_xxxxz_0, \
                             ta_yyyzzz_xxxyy_0, \
                             ta_yyyzzz_xxxyz_0, \
                             ta_yyyzzz_xxxzz_0, \
                             ta_yyyzzz_xxyyy_0, \
                             ta_yyyzzz_xxyyz_0, \
                             ta_yyyzzz_xxyzz_0, \
                             ta_yyyzzz_xxzzz_0, \
                             ta_yyyzzz_xyyyy_0, \
                             ta_yyyzzz_xyyyz_0, \
                             ta_yyyzzz_xyyzz_0, \
                             ta_yyyzzz_xyzzz_0, \
                             ta_yyyzzz_xzzzz_0, \
                             ta_yyyzzz_yyyyy_0, \
                             ta_yyyzzz_yyyyz_0, \
                             ta_yyyzzz_yyyzz_0, \
                             ta_yyyzzz_yyzzz_0, \
                             ta_yyyzzz_yzzzz_0, \
                             ta_yyyzzz_zzzzz_0, \
                             ta_yyzzz_xxxxx_0,  \
                             ta_yyzzz_xxxxx_1,  \
                             ta_yyzzz_xxxxz_0,  \
                             ta_yyzzz_xxxxz_1,  \
                             ta_yyzzz_xxxyz_0,  \
                             ta_yyzzz_xxxyz_1,  \
                             ta_yyzzz_xxxz_0,   \
                             ta_yyzzz_xxxz_1,   \
                             ta_yyzzz_xxxzz_0,  \
                             ta_yyzzz_xxxzz_1,  \
                             ta_yyzzz_xxyyz_0,  \
                             ta_yyzzz_xxyyz_1,  \
                             ta_yyzzz_xxyz_0,   \
                             ta_yyzzz_xxyz_1,   \
                             ta_yyzzz_xxyzz_0,  \
                             ta_yyzzz_xxyzz_1,  \
                             ta_yyzzz_xxzz_0,   \
                             ta_yyzzz_xxzz_1,   \
                             ta_yyzzz_xxzzz_0,  \
                             ta_yyzzz_xxzzz_1,  \
                             ta_yyzzz_xyyyz_0,  \
                             ta_yyzzz_xyyyz_1,  \
                             ta_yyzzz_xyyz_0,   \
                             ta_yyzzz_xyyz_1,   \
                             ta_yyzzz_xyyzz_0,  \
                             ta_yyzzz_xyyzz_1,  \
                             ta_yyzzz_xyzz_0,   \
                             ta_yyzzz_xyzz_1,   \
                             ta_yyzzz_xyzzz_0,  \
                             ta_yyzzz_xyzzz_1,  \
                             ta_yyzzz_xzzz_0,   \
                             ta_yyzzz_xzzz_1,   \
                             ta_yyzzz_xzzzz_0,  \
                             ta_yyzzz_xzzzz_1,  \
                             ta_yyzzz_yyyyz_0,  \
                             ta_yyzzz_yyyyz_1,  \
                             ta_yyzzz_yyyz_0,   \
                             ta_yyzzz_yyyz_1,   \
                             ta_yyzzz_yyyzz_0,  \
                             ta_yyzzz_yyyzz_1,  \
                             ta_yyzzz_yyzz_0,   \
                             ta_yyzzz_yyzz_1,   \
                             ta_yyzzz_yyzzz_0,  \
                             ta_yyzzz_yyzzz_1,  \
                             ta_yyzzz_yzzz_0,   \
                             ta_yyzzz_yzzz_1,   \
                             ta_yyzzz_yzzzz_0,  \
                             ta_yyzzz_yzzzz_1,  \
                             ta_yyzzz_zzzz_0,   \
                             ta_yyzzz_zzzz_1,   \
                             ta_yyzzz_zzzzz_0,  \
                             ta_yyzzz_zzzzz_1,  \
                             ta_yzzz_xxxxx_0,   \
                             ta_yzzz_xxxxx_1,   \
                             ta_yzzz_xxxxz_0,   \
                             ta_yzzz_xxxxz_1,   \
                             ta_yzzz_xxxyz_0,   \
                             ta_yzzz_xxxyz_1,   \
                             ta_yzzz_xxxzz_0,   \
                             ta_yzzz_xxxzz_1,   \
                             ta_yzzz_xxyyz_0,   \
                             ta_yzzz_xxyyz_1,   \
                             ta_yzzz_xxyzz_0,   \
                             ta_yzzz_xxyzz_1,   \
                             ta_yzzz_xxzzz_0,   \
                             ta_yzzz_xxzzz_1,   \
                             ta_yzzz_xyyyz_0,   \
                             ta_yzzz_xyyyz_1,   \
                             ta_yzzz_xyyzz_0,   \
                             ta_yzzz_xyyzz_1,   \
                             ta_yzzz_xyzzz_0,   \
                             ta_yzzz_xyzzz_1,   \
                             ta_yzzz_xzzzz_0,   \
                             ta_yzzz_xzzzz_1,   \
                             ta_yzzz_yyyyz_0,   \
                             ta_yzzz_yyyyz_1,   \
                             ta_yzzz_yyyzz_0,   \
                             ta_yzzz_yyyzz_1,   \
                             ta_yzzz_yyzzz_0,   \
                             ta_yzzz_yyzzz_1,   \
                             ta_yzzz_yzzzz_0,   \
                             ta_yzzz_yzzzz_1,   \
                             ta_yzzz_zzzzz_0,   \
                             ta_yzzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzzz_xxxxx_0[i] =
            2.0 * ta_yzzz_xxxxx_0[i] * fe_0 - 2.0 * ta_yzzz_xxxxx_1[i] * fe_0 + ta_yyzzz_xxxxx_0[i] * pa_y[i] - ta_yyzzz_xxxxx_1[i] * pc_y[i];

        ta_yyyzzz_xxxxy_0[i] =
            2.0 * ta_yyyz_xxxxy_0[i] * fe_0 - 2.0 * ta_yyyz_xxxxy_1[i] * fe_0 + ta_yyyzz_xxxxy_0[i] * pa_z[i] - ta_yyyzz_xxxxy_1[i] * pc_z[i];

        ta_yyyzzz_xxxxz_0[i] =
            2.0 * ta_yzzz_xxxxz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxxz_1[i] * fe_0 + ta_yyzzz_xxxxz_0[i] * pa_y[i] - ta_yyzzz_xxxxz_1[i] * pc_y[i];

        ta_yyyzzz_xxxyy_0[i] =
            2.0 * ta_yyyz_xxxyy_0[i] * fe_0 - 2.0 * ta_yyyz_xxxyy_1[i] * fe_0 + ta_yyyzz_xxxyy_0[i] * pa_z[i] - ta_yyyzz_xxxyy_1[i] * pc_z[i];

        ta_yyyzzz_xxxyz_0[i] = 2.0 * ta_yzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxyz_1[i] * fe_0 + ta_yyzzz_xxxz_0[i] * fe_0 -
                               ta_yyzzz_xxxz_1[i] * fe_0 + ta_yyzzz_xxxyz_0[i] * pa_y[i] - ta_yyzzz_xxxyz_1[i] * pc_y[i];

        ta_yyyzzz_xxxzz_0[i] =
            2.0 * ta_yzzz_xxxzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxxzz_1[i] * fe_0 + ta_yyzzz_xxxzz_0[i] * pa_y[i] - ta_yyzzz_xxxzz_1[i] * pc_y[i];

        ta_yyyzzz_xxyyy_0[i] =
            2.0 * ta_yyyz_xxyyy_0[i] * fe_0 - 2.0 * ta_yyyz_xxyyy_1[i] * fe_0 + ta_yyyzz_xxyyy_0[i] * pa_z[i] - ta_yyyzz_xxyyy_1[i] * pc_z[i];

        ta_yyyzzz_xxyyz_0[i] = 2.0 * ta_yzzz_xxyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyyz_1[i] * fe_0 + 2.0 * ta_yyzzz_xxyz_0[i] * fe_0 -
                               2.0 * ta_yyzzz_xxyz_1[i] * fe_0 + ta_yyzzz_xxyyz_0[i] * pa_y[i] - ta_yyzzz_xxyyz_1[i] * pc_y[i];

        ta_yyyzzz_xxyzz_0[i] = 2.0 * ta_yzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxyzz_1[i] * fe_0 + ta_yyzzz_xxzz_0[i] * fe_0 -
                               ta_yyzzz_xxzz_1[i] * fe_0 + ta_yyzzz_xxyzz_0[i] * pa_y[i] - ta_yyzzz_xxyzz_1[i] * pc_y[i];

        ta_yyyzzz_xxzzz_0[i] =
            2.0 * ta_yzzz_xxzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xxzzz_1[i] * fe_0 + ta_yyzzz_xxzzz_0[i] * pa_y[i] - ta_yyzzz_xxzzz_1[i] * pc_y[i];

        ta_yyyzzz_xyyyy_0[i] =
            2.0 * ta_yyyz_xyyyy_0[i] * fe_0 - 2.0 * ta_yyyz_xyyyy_1[i] * fe_0 + ta_yyyzz_xyyyy_0[i] * pa_z[i] - ta_yyyzz_xyyyy_1[i] * pc_z[i];

        ta_yyyzzz_xyyyz_0[i] = 2.0 * ta_yzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyyz_1[i] * fe_0 + 3.0 * ta_yyzzz_xyyz_0[i] * fe_0 -
                               3.0 * ta_yyzzz_xyyz_1[i] * fe_0 + ta_yyzzz_xyyyz_0[i] * pa_y[i] - ta_yyzzz_xyyyz_1[i] * pc_y[i];

        ta_yyyzzz_xyyzz_0[i] = 2.0 * ta_yzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyzz_1[i] * fe_0 + 2.0 * ta_yyzzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_yyzzz_xyzz_1[i] * fe_0 + ta_yyzzz_xyyzz_0[i] * pa_y[i] - ta_yyzzz_xyyzz_1[i] * pc_y[i];

        ta_yyyzzz_xyzzz_0[i] = 2.0 * ta_yzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyzzz_1[i] * fe_0 + ta_yyzzz_xzzz_0[i] * fe_0 -
                               ta_yyzzz_xzzz_1[i] * fe_0 + ta_yyzzz_xyzzz_0[i] * pa_y[i] - ta_yyzzz_xyzzz_1[i] * pc_y[i];

        ta_yyyzzz_xzzzz_0[i] =
            2.0 * ta_yzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xzzzz_1[i] * fe_0 + ta_yyzzz_xzzzz_0[i] * pa_y[i] - ta_yyzzz_xzzzz_1[i] * pc_y[i];

        ta_yyyzzz_yyyyy_0[i] =
            2.0 * ta_yyyz_yyyyy_0[i] * fe_0 - 2.0 * ta_yyyz_yyyyy_1[i] * fe_0 + ta_yyyzz_yyyyy_0[i] * pa_z[i] - ta_yyyzz_yyyyy_1[i] * pc_z[i];

        ta_yyyzzz_yyyyz_0[i] = 2.0 * ta_yzzz_yyyyz_0[i] * fe_0 - 2.0 * ta_yzzz_yyyyz_1[i] * fe_0 + 4.0 * ta_yyzzz_yyyz_0[i] * fe_0 -
                               4.0 * ta_yyzzz_yyyz_1[i] * fe_0 + ta_yyzzz_yyyyz_0[i] * pa_y[i] - ta_yyzzz_yyyyz_1[i] * pc_y[i];

        ta_yyyzzz_yyyzz_0[i] = 2.0 * ta_yzzz_yyyzz_0[i] * fe_0 - 2.0 * ta_yzzz_yyyzz_1[i] * fe_0 + 3.0 * ta_yyzzz_yyzz_0[i] * fe_0 -
                               3.0 * ta_yyzzz_yyzz_1[i] * fe_0 + ta_yyzzz_yyyzz_0[i] * pa_y[i] - ta_yyzzz_yyyzz_1[i] * pc_y[i];

        ta_yyyzzz_yyzzz_0[i] = 2.0 * ta_yzzz_yyzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yyzzz_1[i] * fe_0 + 2.0 * ta_yyzzz_yzzz_0[i] * fe_0 -
                               2.0 * ta_yyzzz_yzzz_1[i] * fe_0 + ta_yyzzz_yyzzz_0[i] * pa_y[i] - ta_yyzzz_yyzzz_1[i] * pc_y[i];

        ta_yyyzzz_yzzzz_0[i] = 2.0 * ta_yzzz_yzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_yzzzz_1[i] * fe_0 + ta_yyzzz_zzzz_0[i] * fe_0 -
                               ta_yyzzz_zzzz_1[i] * fe_0 + ta_yyzzz_yzzzz_0[i] * pa_y[i] - ta_yyzzz_yzzzz_1[i] * pc_y[i];

        ta_yyyzzz_zzzzz_0[i] =
            2.0 * ta_yzzz_zzzzz_0[i] * fe_0 - 2.0 * ta_yzzz_zzzzz_1[i] * fe_0 + ta_yyzzz_zzzzz_0[i] * pa_y[i] - ta_yyzzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 525-546 components of targeted buffer : IH

    auto ta_yyzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 525);

    auto ta_yyzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 526);

    auto ta_yyzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 527);

    auto ta_yyzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 528);

    auto ta_yyzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 529);

    auto ta_yyzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 530);

    auto ta_yyzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 531);

    auto ta_yyzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 532);

    auto ta_yyzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 533);

    auto ta_yyzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 534);

    auto ta_yyzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 535);

    auto ta_yyzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 536);

    auto ta_yyzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 537);

    auto ta_yyzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 538);

    auto ta_yyzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 539);

    auto ta_yyzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 540);

    auto ta_yyzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 541);

    auto ta_yyzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 542);

    auto ta_yyzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 543);

    auto ta_yyzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 544);

    auto ta_yyzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 545);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyzz_xxxxy_0,   \
                             ta_yyzz_xxxxy_1,   \
                             ta_yyzz_xxxyy_0,   \
                             ta_yyzz_xxxyy_1,   \
                             ta_yyzz_xxyyy_0,   \
                             ta_yyzz_xxyyy_1,   \
                             ta_yyzz_xyyyy_0,   \
                             ta_yyzz_xyyyy_1,   \
                             ta_yyzz_yyyyy_0,   \
                             ta_yyzz_yyyyy_1,   \
                             ta_yyzzz_xxxxy_0,  \
                             ta_yyzzz_xxxxy_1,  \
                             ta_yyzzz_xxxyy_0,  \
                             ta_yyzzz_xxxyy_1,  \
                             ta_yyzzz_xxyyy_0,  \
                             ta_yyzzz_xxyyy_1,  \
                             ta_yyzzz_xyyyy_0,  \
                             ta_yyzzz_xyyyy_1,  \
                             ta_yyzzz_yyyyy_0,  \
                             ta_yyzzz_yyyyy_1,  \
                             ta_yyzzzz_xxxxx_0, \
                             ta_yyzzzz_xxxxy_0, \
                             ta_yyzzzz_xxxxz_0, \
                             ta_yyzzzz_xxxyy_0, \
                             ta_yyzzzz_xxxyz_0, \
                             ta_yyzzzz_xxxzz_0, \
                             ta_yyzzzz_xxyyy_0, \
                             ta_yyzzzz_xxyyz_0, \
                             ta_yyzzzz_xxyzz_0, \
                             ta_yyzzzz_xxzzz_0, \
                             ta_yyzzzz_xyyyy_0, \
                             ta_yyzzzz_xyyyz_0, \
                             ta_yyzzzz_xyyzz_0, \
                             ta_yyzzzz_xyzzz_0, \
                             ta_yyzzzz_xzzzz_0, \
                             ta_yyzzzz_yyyyy_0, \
                             ta_yyzzzz_yyyyz_0, \
                             ta_yyzzzz_yyyzz_0, \
                             ta_yyzzzz_yyzzz_0, \
                             ta_yyzzzz_yzzzz_0, \
                             ta_yyzzzz_zzzzz_0, \
                             ta_yzzzz_xxxxx_0,  \
                             ta_yzzzz_xxxxx_1,  \
                             ta_yzzzz_xxxxz_0,  \
                             ta_yzzzz_xxxxz_1,  \
                             ta_yzzzz_xxxyz_0,  \
                             ta_yzzzz_xxxyz_1,  \
                             ta_yzzzz_xxxz_0,   \
                             ta_yzzzz_xxxz_1,   \
                             ta_yzzzz_xxxzz_0,  \
                             ta_yzzzz_xxxzz_1,  \
                             ta_yzzzz_xxyyz_0,  \
                             ta_yzzzz_xxyyz_1,  \
                             ta_yzzzz_xxyz_0,   \
                             ta_yzzzz_xxyz_1,   \
                             ta_yzzzz_xxyzz_0,  \
                             ta_yzzzz_xxyzz_1,  \
                             ta_yzzzz_xxzz_0,   \
                             ta_yzzzz_xxzz_1,   \
                             ta_yzzzz_xxzzz_0,  \
                             ta_yzzzz_xxzzz_1,  \
                             ta_yzzzz_xyyyz_0,  \
                             ta_yzzzz_xyyyz_1,  \
                             ta_yzzzz_xyyz_0,   \
                             ta_yzzzz_xyyz_1,   \
                             ta_yzzzz_xyyzz_0,  \
                             ta_yzzzz_xyyzz_1,  \
                             ta_yzzzz_xyzz_0,   \
                             ta_yzzzz_xyzz_1,   \
                             ta_yzzzz_xyzzz_0,  \
                             ta_yzzzz_xyzzz_1,  \
                             ta_yzzzz_xzzz_0,   \
                             ta_yzzzz_xzzz_1,   \
                             ta_yzzzz_xzzzz_0,  \
                             ta_yzzzz_xzzzz_1,  \
                             ta_yzzzz_yyyyz_0,  \
                             ta_yzzzz_yyyyz_1,  \
                             ta_yzzzz_yyyz_0,   \
                             ta_yzzzz_yyyz_1,   \
                             ta_yzzzz_yyyzz_0,  \
                             ta_yzzzz_yyyzz_1,  \
                             ta_yzzzz_yyzz_0,   \
                             ta_yzzzz_yyzz_1,   \
                             ta_yzzzz_yyzzz_0,  \
                             ta_yzzzz_yyzzz_1,  \
                             ta_yzzzz_yzzz_0,   \
                             ta_yzzzz_yzzz_1,   \
                             ta_yzzzz_yzzzz_0,  \
                             ta_yzzzz_yzzzz_1,  \
                             ta_yzzzz_zzzz_0,   \
                             ta_yzzzz_zzzz_1,   \
                             ta_yzzzz_zzzzz_0,  \
                             ta_yzzzz_zzzzz_1,  \
                             ta_zzzz_xxxxx_0,   \
                             ta_zzzz_xxxxx_1,   \
                             ta_zzzz_xxxxz_0,   \
                             ta_zzzz_xxxxz_1,   \
                             ta_zzzz_xxxyz_0,   \
                             ta_zzzz_xxxyz_1,   \
                             ta_zzzz_xxxzz_0,   \
                             ta_zzzz_xxxzz_1,   \
                             ta_zzzz_xxyyz_0,   \
                             ta_zzzz_xxyyz_1,   \
                             ta_zzzz_xxyzz_0,   \
                             ta_zzzz_xxyzz_1,   \
                             ta_zzzz_xxzzz_0,   \
                             ta_zzzz_xxzzz_1,   \
                             ta_zzzz_xyyyz_0,   \
                             ta_zzzz_xyyyz_1,   \
                             ta_zzzz_xyyzz_0,   \
                             ta_zzzz_xyyzz_1,   \
                             ta_zzzz_xyzzz_0,   \
                             ta_zzzz_xyzzz_1,   \
                             ta_zzzz_xzzzz_0,   \
                             ta_zzzz_xzzzz_1,   \
                             ta_zzzz_yyyyz_0,   \
                             ta_zzzz_yyyyz_1,   \
                             ta_zzzz_yyyzz_0,   \
                             ta_zzzz_yyyzz_1,   \
                             ta_zzzz_yyzzz_0,   \
                             ta_zzzz_yyzzz_1,   \
                             ta_zzzz_yzzzz_0,   \
                             ta_zzzz_yzzzz_1,   \
                             ta_zzzz_zzzzz_0,   \
                             ta_zzzz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzzz_xxxxx_0[i] = ta_zzzz_xxxxx_0[i] * fe_0 - ta_zzzz_xxxxx_1[i] * fe_0 + ta_yzzzz_xxxxx_0[i] * pa_y[i] - ta_yzzzz_xxxxx_1[i] * pc_y[i];

        ta_yyzzzz_xxxxy_0[i] =
            3.0 * ta_yyzz_xxxxy_0[i] * fe_0 - 3.0 * ta_yyzz_xxxxy_1[i] * fe_0 + ta_yyzzz_xxxxy_0[i] * pa_z[i] - ta_yyzzz_xxxxy_1[i] * pc_z[i];

        ta_yyzzzz_xxxxz_0[i] = ta_zzzz_xxxxz_0[i] * fe_0 - ta_zzzz_xxxxz_1[i] * fe_0 + ta_yzzzz_xxxxz_0[i] * pa_y[i] - ta_yzzzz_xxxxz_1[i] * pc_y[i];

        ta_yyzzzz_xxxyy_0[i] =
            3.0 * ta_yyzz_xxxyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxxyy_1[i] * fe_0 + ta_yyzzz_xxxyy_0[i] * pa_z[i] - ta_yyzzz_xxxyy_1[i] * pc_z[i];

        ta_yyzzzz_xxxyz_0[i] = ta_zzzz_xxxyz_0[i] * fe_0 - ta_zzzz_xxxyz_1[i] * fe_0 + ta_yzzzz_xxxz_0[i] * fe_0 - ta_yzzzz_xxxz_1[i] * fe_0 +
                               ta_yzzzz_xxxyz_0[i] * pa_y[i] - ta_yzzzz_xxxyz_1[i] * pc_y[i];

        ta_yyzzzz_xxxzz_0[i] = ta_zzzz_xxxzz_0[i] * fe_0 - ta_zzzz_xxxzz_1[i] * fe_0 + ta_yzzzz_xxxzz_0[i] * pa_y[i] - ta_yzzzz_xxxzz_1[i] * pc_y[i];

        ta_yyzzzz_xxyyy_0[i] =
            3.0 * ta_yyzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyy_1[i] * fe_0 + ta_yyzzz_xxyyy_0[i] * pa_z[i] - ta_yyzzz_xxyyy_1[i] * pc_z[i];

        ta_yyzzzz_xxyyz_0[i] = ta_zzzz_xxyyz_0[i] * fe_0 - ta_zzzz_xxyyz_1[i] * fe_0 + 2.0 * ta_yzzzz_xxyz_0[i] * fe_0 -
                               2.0 * ta_yzzzz_xxyz_1[i] * fe_0 + ta_yzzzz_xxyyz_0[i] * pa_y[i] - ta_yzzzz_xxyyz_1[i] * pc_y[i];

        ta_yyzzzz_xxyzz_0[i] = ta_zzzz_xxyzz_0[i] * fe_0 - ta_zzzz_xxyzz_1[i] * fe_0 + ta_yzzzz_xxzz_0[i] * fe_0 - ta_yzzzz_xxzz_1[i] * fe_0 +
                               ta_yzzzz_xxyzz_0[i] * pa_y[i] - ta_yzzzz_xxyzz_1[i] * pc_y[i];

        ta_yyzzzz_xxzzz_0[i] = ta_zzzz_xxzzz_0[i] * fe_0 - ta_zzzz_xxzzz_1[i] * fe_0 + ta_yzzzz_xxzzz_0[i] * pa_y[i] - ta_yzzzz_xxzzz_1[i] * pc_y[i];

        ta_yyzzzz_xyyyy_0[i] =
            3.0 * ta_yyzz_xyyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xyyyy_1[i] * fe_0 + ta_yyzzz_xyyyy_0[i] * pa_z[i] - ta_yyzzz_xyyyy_1[i] * pc_z[i];

        ta_yyzzzz_xyyyz_0[i] = ta_zzzz_xyyyz_0[i] * fe_0 - ta_zzzz_xyyyz_1[i] * fe_0 + 3.0 * ta_yzzzz_xyyz_0[i] * fe_0 -
                               3.0 * ta_yzzzz_xyyz_1[i] * fe_0 + ta_yzzzz_xyyyz_0[i] * pa_y[i] - ta_yzzzz_xyyyz_1[i] * pc_y[i];

        ta_yyzzzz_xyyzz_0[i] = ta_zzzz_xyyzz_0[i] * fe_0 - ta_zzzz_xyyzz_1[i] * fe_0 + 2.0 * ta_yzzzz_xyzz_0[i] * fe_0 -
                               2.0 * ta_yzzzz_xyzz_1[i] * fe_0 + ta_yzzzz_xyyzz_0[i] * pa_y[i] - ta_yzzzz_xyyzz_1[i] * pc_y[i];

        ta_yyzzzz_xyzzz_0[i] = ta_zzzz_xyzzz_0[i] * fe_0 - ta_zzzz_xyzzz_1[i] * fe_0 + ta_yzzzz_xzzz_0[i] * fe_0 - ta_yzzzz_xzzz_1[i] * fe_0 +
                               ta_yzzzz_xyzzz_0[i] * pa_y[i] - ta_yzzzz_xyzzz_1[i] * pc_y[i];

        ta_yyzzzz_xzzzz_0[i] = ta_zzzz_xzzzz_0[i] * fe_0 - ta_zzzz_xzzzz_1[i] * fe_0 + ta_yzzzz_xzzzz_0[i] * pa_y[i] - ta_yzzzz_xzzzz_1[i] * pc_y[i];

        ta_yyzzzz_yyyyy_0[i] =
            3.0 * ta_yyzz_yyyyy_0[i] * fe_0 - 3.0 * ta_yyzz_yyyyy_1[i] * fe_0 + ta_yyzzz_yyyyy_0[i] * pa_z[i] - ta_yyzzz_yyyyy_1[i] * pc_z[i];

        ta_yyzzzz_yyyyz_0[i] = ta_zzzz_yyyyz_0[i] * fe_0 - ta_zzzz_yyyyz_1[i] * fe_0 + 4.0 * ta_yzzzz_yyyz_0[i] * fe_0 -
                               4.0 * ta_yzzzz_yyyz_1[i] * fe_0 + ta_yzzzz_yyyyz_0[i] * pa_y[i] - ta_yzzzz_yyyyz_1[i] * pc_y[i];

        ta_yyzzzz_yyyzz_0[i] = ta_zzzz_yyyzz_0[i] * fe_0 - ta_zzzz_yyyzz_1[i] * fe_0 + 3.0 * ta_yzzzz_yyzz_0[i] * fe_0 -
                               3.0 * ta_yzzzz_yyzz_1[i] * fe_0 + ta_yzzzz_yyyzz_0[i] * pa_y[i] - ta_yzzzz_yyyzz_1[i] * pc_y[i];

        ta_yyzzzz_yyzzz_0[i] = ta_zzzz_yyzzz_0[i] * fe_0 - ta_zzzz_yyzzz_1[i] * fe_0 + 2.0 * ta_yzzzz_yzzz_0[i] * fe_0 -
                               2.0 * ta_yzzzz_yzzz_1[i] * fe_0 + ta_yzzzz_yyzzz_0[i] * pa_y[i] - ta_yzzzz_yyzzz_1[i] * pc_y[i];

        ta_yyzzzz_yzzzz_0[i] = ta_zzzz_yzzzz_0[i] * fe_0 - ta_zzzz_yzzzz_1[i] * fe_0 + ta_yzzzz_zzzz_0[i] * fe_0 - ta_yzzzz_zzzz_1[i] * fe_0 +
                               ta_yzzzz_yzzzz_0[i] * pa_y[i] - ta_yzzzz_yzzzz_1[i] * pc_y[i];

        ta_yyzzzz_zzzzz_0[i] = ta_zzzz_zzzzz_0[i] * fe_0 - ta_zzzz_zzzzz_1[i] * fe_0 + ta_yzzzz_zzzzz_0[i] * pa_y[i] - ta_yzzzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 546-567 components of targeted buffer : IH

    auto ta_yzzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 546);

    auto ta_yzzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 547);

    auto ta_yzzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 548);

    auto ta_yzzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 549);

    auto ta_yzzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 550);

    auto ta_yzzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 551);

    auto ta_yzzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 552);

    auto ta_yzzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 553);

    auto ta_yzzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 554);

    auto ta_yzzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 555);

    auto ta_yzzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 556);

    auto ta_yzzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 557);

    auto ta_yzzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 558);

    auto ta_yzzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 559);

    auto ta_yzzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 560);

    auto ta_yzzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 561);

    auto ta_yzzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 562);

    auto ta_yzzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 563);

    auto ta_yzzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 564);

    auto ta_yzzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 565);

    auto ta_yzzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 566);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta_yzzzzz_xxxxx_0, \
                             ta_yzzzzz_xxxxy_0, \
                             ta_yzzzzz_xxxxz_0, \
                             ta_yzzzzz_xxxyy_0, \
                             ta_yzzzzz_xxxyz_0, \
                             ta_yzzzzz_xxxzz_0, \
                             ta_yzzzzz_xxyyy_0, \
                             ta_yzzzzz_xxyyz_0, \
                             ta_yzzzzz_xxyzz_0, \
                             ta_yzzzzz_xxzzz_0, \
                             ta_yzzzzz_xyyyy_0, \
                             ta_yzzzzz_xyyyz_0, \
                             ta_yzzzzz_xyyzz_0, \
                             ta_yzzzzz_xyzzz_0, \
                             ta_yzzzzz_xzzzz_0, \
                             ta_yzzzzz_yyyyy_0, \
                             ta_yzzzzz_yyyyz_0, \
                             ta_yzzzzz_yyyzz_0, \
                             ta_yzzzzz_yyzzz_0, \
                             ta_yzzzzz_yzzzz_0, \
                             ta_yzzzzz_zzzzz_0, \
                             ta_zzzzz_xxxx_0,   \
                             ta_zzzzz_xxxx_1,   \
                             ta_zzzzz_xxxxx_0,  \
                             ta_zzzzz_xxxxx_1,  \
                             ta_zzzzz_xxxxy_0,  \
                             ta_zzzzz_xxxxy_1,  \
                             ta_zzzzz_xxxxz_0,  \
                             ta_zzzzz_xxxxz_1,  \
                             ta_zzzzz_xxxy_0,   \
                             ta_zzzzz_xxxy_1,   \
                             ta_zzzzz_xxxyy_0,  \
                             ta_zzzzz_xxxyy_1,  \
                             ta_zzzzz_xxxyz_0,  \
                             ta_zzzzz_xxxyz_1,  \
                             ta_zzzzz_xxxz_0,   \
                             ta_zzzzz_xxxz_1,   \
                             ta_zzzzz_xxxzz_0,  \
                             ta_zzzzz_xxxzz_1,  \
                             ta_zzzzz_xxyy_0,   \
                             ta_zzzzz_xxyy_1,   \
                             ta_zzzzz_xxyyy_0,  \
                             ta_zzzzz_xxyyy_1,  \
                             ta_zzzzz_xxyyz_0,  \
                             ta_zzzzz_xxyyz_1,  \
                             ta_zzzzz_xxyz_0,   \
                             ta_zzzzz_xxyz_1,   \
                             ta_zzzzz_xxyzz_0,  \
                             ta_zzzzz_xxyzz_1,  \
                             ta_zzzzz_xxzz_0,   \
                             ta_zzzzz_xxzz_1,   \
                             ta_zzzzz_xxzzz_0,  \
                             ta_zzzzz_xxzzz_1,  \
                             ta_zzzzz_xyyy_0,   \
                             ta_zzzzz_xyyy_1,   \
                             ta_zzzzz_xyyyy_0,  \
                             ta_zzzzz_xyyyy_1,  \
                             ta_zzzzz_xyyyz_0,  \
                             ta_zzzzz_xyyyz_1,  \
                             ta_zzzzz_xyyz_0,   \
                             ta_zzzzz_xyyz_1,   \
                             ta_zzzzz_xyyzz_0,  \
                             ta_zzzzz_xyyzz_1,  \
                             ta_zzzzz_xyzz_0,   \
                             ta_zzzzz_xyzz_1,   \
                             ta_zzzzz_xyzzz_0,  \
                             ta_zzzzz_xyzzz_1,  \
                             ta_zzzzz_xzzz_0,   \
                             ta_zzzzz_xzzz_1,   \
                             ta_zzzzz_xzzzz_0,  \
                             ta_zzzzz_xzzzz_1,  \
                             ta_zzzzz_yyyy_0,   \
                             ta_zzzzz_yyyy_1,   \
                             ta_zzzzz_yyyyy_0,  \
                             ta_zzzzz_yyyyy_1,  \
                             ta_zzzzz_yyyyz_0,  \
                             ta_zzzzz_yyyyz_1,  \
                             ta_zzzzz_yyyz_0,   \
                             ta_zzzzz_yyyz_1,   \
                             ta_zzzzz_yyyzz_0,  \
                             ta_zzzzz_yyyzz_1,  \
                             ta_zzzzz_yyzz_0,   \
                             ta_zzzzz_yyzz_1,   \
                             ta_zzzzz_yyzzz_0,  \
                             ta_zzzzz_yyzzz_1,  \
                             ta_zzzzz_yzzz_0,   \
                             ta_zzzzz_yzzz_1,   \
                             ta_zzzzz_yzzzz_0,  \
                             ta_zzzzz_yzzzz_1,  \
                             ta_zzzzz_zzzz_0,   \
                             ta_zzzzz_zzzz_1,   \
                             ta_zzzzz_zzzzz_0,  \
                             ta_zzzzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzzz_xxxxx_0[i] = ta_zzzzz_xxxxx_0[i] * pa_y[i] - ta_zzzzz_xxxxx_1[i] * pc_y[i];

        ta_yzzzzz_xxxxy_0[i] = ta_zzzzz_xxxx_0[i] * fe_0 - ta_zzzzz_xxxx_1[i] * fe_0 + ta_zzzzz_xxxxy_0[i] * pa_y[i] - ta_zzzzz_xxxxy_1[i] * pc_y[i];

        ta_yzzzzz_xxxxz_0[i] = ta_zzzzz_xxxxz_0[i] * pa_y[i] - ta_zzzzz_xxxxz_1[i] * pc_y[i];

        ta_yzzzzz_xxxyy_0[i] =
            2.0 * ta_zzzzz_xxxy_0[i] * fe_0 - 2.0 * ta_zzzzz_xxxy_1[i] * fe_0 + ta_zzzzz_xxxyy_0[i] * pa_y[i] - ta_zzzzz_xxxyy_1[i] * pc_y[i];

        ta_yzzzzz_xxxyz_0[i] = ta_zzzzz_xxxz_0[i] * fe_0 - ta_zzzzz_xxxz_1[i] * fe_0 + ta_zzzzz_xxxyz_0[i] * pa_y[i] - ta_zzzzz_xxxyz_1[i] * pc_y[i];

        ta_yzzzzz_xxxzz_0[i] = ta_zzzzz_xxxzz_0[i] * pa_y[i] - ta_zzzzz_xxxzz_1[i] * pc_y[i];

        ta_yzzzzz_xxyyy_0[i] =
            3.0 * ta_zzzzz_xxyy_0[i] * fe_0 - 3.0 * ta_zzzzz_xxyy_1[i] * fe_0 + ta_zzzzz_xxyyy_0[i] * pa_y[i] - ta_zzzzz_xxyyy_1[i] * pc_y[i];

        ta_yzzzzz_xxyyz_0[i] =
            2.0 * ta_zzzzz_xxyz_0[i] * fe_0 - 2.0 * ta_zzzzz_xxyz_1[i] * fe_0 + ta_zzzzz_xxyyz_0[i] * pa_y[i] - ta_zzzzz_xxyyz_1[i] * pc_y[i];

        ta_yzzzzz_xxyzz_0[i] = ta_zzzzz_xxzz_0[i] * fe_0 - ta_zzzzz_xxzz_1[i] * fe_0 + ta_zzzzz_xxyzz_0[i] * pa_y[i] - ta_zzzzz_xxyzz_1[i] * pc_y[i];

        ta_yzzzzz_xxzzz_0[i] = ta_zzzzz_xxzzz_0[i] * pa_y[i] - ta_zzzzz_xxzzz_1[i] * pc_y[i];

        ta_yzzzzz_xyyyy_0[i] =
            4.0 * ta_zzzzz_xyyy_0[i] * fe_0 - 4.0 * ta_zzzzz_xyyy_1[i] * fe_0 + ta_zzzzz_xyyyy_0[i] * pa_y[i] - ta_zzzzz_xyyyy_1[i] * pc_y[i];

        ta_yzzzzz_xyyyz_0[i] =
            3.0 * ta_zzzzz_xyyz_0[i] * fe_0 - 3.0 * ta_zzzzz_xyyz_1[i] * fe_0 + ta_zzzzz_xyyyz_0[i] * pa_y[i] - ta_zzzzz_xyyyz_1[i] * pc_y[i];

        ta_yzzzzz_xyyzz_0[i] =
            2.0 * ta_zzzzz_xyzz_0[i] * fe_0 - 2.0 * ta_zzzzz_xyzz_1[i] * fe_0 + ta_zzzzz_xyyzz_0[i] * pa_y[i] - ta_zzzzz_xyyzz_1[i] * pc_y[i];

        ta_yzzzzz_xyzzz_0[i] = ta_zzzzz_xzzz_0[i] * fe_0 - ta_zzzzz_xzzz_1[i] * fe_0 + ta_zzzzz_xyzzz_0[i] * pa_y[i] - ta_zzzzz_xyzzz_1[i] * pc_y[i];

        ta_yzzzzz_xzzzz_0[i] = ta_zzzzz_xzzzz_0[i] * pa_y[i] - ta_zzzzz_xzzzz_1[i] * pc_y[i];

        ta_yzzzzz_yyyyy_0[i] =
            5.0 * ta_zzzzz_yyyy_0[i] * fe_0 - 5.0 * ta_zzzzz_yyyy_1[i] * fe_0 + ta_zzzzz_yyyyy_0[i] * pa_y[i] - ta_zzzzz_yyyyy_1[i] * pc_y[i];

        ta_yzzzzz_yyyyz_0[i] =
            4.0 * ta_zzzzz_yyyz_0[i] * fe_0 - 4.0 * ta_zzzzz_yyyz_1[i] * fe_0 + ta_zzzzz_yyyyz_0[i] * pa_y[i] - ta_zzzzz_yyyyz_1[i] * pc_y[i];

        ta_yzzzzz_yyyzz_0[i] =
            3.0 * ta_zzzzz_yyzz_0[i] * fe_0 - 3.0 * ta_zzzzz_yyzz_1[i] * fe_0 + ta_zzzzz_yyyzz_0[i] * pa_y[i] - ta_zzzzz_yyyzz_1[i] * pc_y[i];

        ta_yzzzzz_yyzzz_0[i] =
            2.0 * ta_zzzzz_yzzz_0[i] * fe_0 - 2.0 * ta_zzzzz_yzzz_1[i] * fe_0 + ta_zzzzz_yyzzz_0[i] * pa_y[i] - ta_zzzzz_yyzzz_1[i] * pc_y[i];

        ta_yzzzzz_yzzzz_0[i] = ta_zzzzz_zzzz_0[i] * fe_0 - ta_zzzzz_zzzz_1[i] * fe_0 + ta_zzzzz_yzzzz_0[i] * pa_y[i] - ta_zzzzz_yzzzz_1[i] * pc_y[i];

        ta_yzzzzz_zzzzz_0[i] = ta_zzzzz_zzzzz_0[i] * pa_y[i] - ta_zzzzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 567-588 components of targeted buffer : IH

    auto ta_zzzzzz_xxxxx_0 = pbuffer.data(idx_npot_0_ih + 567);

    auto ta_zzzzzz_xxxxy_0 = pbuffer.data(idx_npot_0_ih + 568);

    auto ta_zzzzzz_xxxxz_0 = pbuffer.data(idx_npot_0_ih + 569);

    auto ta_zzzzzz_xxxyy_0 = pbuffer.data(idx_npot_0_ih + 570);

    auto ta_zzzzzz_xxxyz_0 = pbuffer.data(idx_npot_0_ih + 571);

    auto ta_zzzzzz_xxxzz_0 = pbuffer.data(idx_npot_0_ih + 572);

    auto ta_zzzzzz_xxyyy_0 = pbuffer.data(idx_npot_0_ih + 573);

    auto ta_zzzzzz_xxyyz_0 = pbuffer.data(idx_npot_0_ih + 574);

    auto ta_zzzzzz_xxyzz_0 = pbuffer.data(idx_npot_0_ih + 575);

    auto ta_zzzzzz_xxzzz_0 = pbuffer.data(idx_npot_0_ih + 576);

    auto ta_zzzzzz_xyyyy_0 = pbuffer.data(idx_npot_0_ih + 577);

    auto ta_zzzzzz_xyyyz_0 = pbuffer.data(idx_npot_0_ih + 578);

    auto ta_zzzzzz_xyyzz_0 = pbuffer.data(idx_npot_0_ih + 579);

    auto ta_zzzzzz_xyzzz_0 = pbuffer.data(idx_npot_0_ih + 580);

    auto ta_zzzzzz_xzzzz_0 = pbuffer.data(idx_npot_0_ih + 581);

    auto ta_zzzzzz_yyyyy_0 = pbuffer.data(idx_npot_0_ih + 582);

    auto ta_zzzzzz_yyyyz_0 = pbuffer.data(idx_npot_0_ih + 583);

    auto ta_zzzzzz_yyyzz_0 = pbuffer.data(idx_npot_0_ih + 584);

    auto ta_zzzzzz_yyzzz_0 = pbuffer.data(idx_npot_0_ih + 585);

    auto ta_zzzzzz_yzzzz_0 = pbuffer.data(idx_npot_0_ih + 586);

    auto ta_zzzzzz_zzzzz_0 = pbuffer.data(idx_npot_0_ih + 587);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta_zzzz_xxxxx_0,   \
                             ta_zzzz_xxxxx_1,   \
                             ta_zzzz_xxxxy_0,   \
                             ta_zzzz_xxxxy_1,   \
                             ta_zzzz_xxxxz_0,   \
                             ta_zzzz_xxxxz_1,   \
                             ta_zzzz_xxxyy_0,   \
                             ta_zzzz_xxxyy_1,   \
                             ta_zzzz_xxxyz_0,   \
                             ta_zzzz_xxxyz_1,   \
                             ta_zzzz_xxxzz_0,   \
                             ta_zzzz_xxxzz_1,   \
                             ta_zzzz_xxyyy_0,   \
                             ta_zzzz_xxyyy_1,   \
                             ta_zzzz_xxyyz_0,   \
                             ta_zzzz_xxyyz_1,   \
                             ta_zzzz_xxyzz_0,   \
                             ta_zzzz_xxyzz_1,   \
                             ta_zzzz_xxzzz_0,   \
                             ta_zzzz_xxzzz_1,   \
                             ta_zzzz_xyyyy_0,   \
                             ta_zzzz_xyyyy_1,   \
                             ta_zzzz_xyyyz_0,   \
                             ta_zzzz_xyyyz_1,   \
                             ta_zzzz_xyyzz_0,   \
                             ta_zzzz_xyyzz_1,   \
                             ta_zzzz_xyzzz_0,   \
                             ta_zzzz_xyzzz_1,   \
                             ta_zzzz_xzzzz_0,   \
                             ta_zzzz_xzzzz_1,   \
                             ta_zzzz_yyyyy_0,   \
                             ta_zzzz_yyyyy_1,   \
                             ta_zzzz_yyyyz_0,   \
                             ta_zzzz_yyyyz_1,   \
                             ta_zzzz_yyyzz_0,   \
                             ta_zzzz_yyyzz_1,   \
                             ta_zzzz_yyzzz_0,   \
                             ta_zzzz_yyzzz_1,   \
                             ta_zzzz_yzzzz_0,   \
                             ta_zzzz_yzzzz_1,   \
                             ta_zzzz_zzzzz_0,   \
                             ta_zzzz_zzzzz_1,   \
                             ta_zzzzz_xxxx_0,   \
                             ta_zzzzz_xxxx_1,   \
                             ta_zzzzz_xxxxx_0,  \
                             ta_zzzzz_xxxxx_1,  \
                             ta_zzzzz_xxxxy_0,  \
                             ta_zzzzz_xxxxy_1,  \
                             ta_zzzzz_xxxxz_0,  \
                             ta_zzzzz_xxxxz_1,  \
                             ta_zzzzz_xxxy_0,   \
                             ta_zzzzz_xxxy_1,   \
                             ta_zzzzz_xxxyy_0,  \
                             ta_zzzzz_xxxyy_1,  \
                             ta_zzzzz_xxxyz_0,  \
                             ta_zzzzz_xxxyz_1,  \
                             ta_zzzzz_xxxz_0,   \
                             ta_zzzzz_xxxz_1,   \
                             ta_zzzzz_xxxzz_0,  \
                             ta_zzzzz_xxxzz_1,  \
                             ta_zzzzz_xxyy_0,   \
                             ta_zzzzz_xxyy_1,   \
                             ta_zzzzz_xxyyy_0,  \
                             ta_zzzzz_xxyyy_1,  \
                             ta_zzzzz_xxyyz_0,  \
                             ta_zzzzz_xxyyz_1,  \
                             ta_zzzzz_xxyz_0,   \
                             ta_zzzzz_xxyz_1,   \
                             ta_zzzzz_xxyzz_0,  \
                             ta_zzzzz_xxyzz_1,  \
                             ta_zzzzz_xxzz_0,   \
                             ta_zzzzz_xxzz_1,   \
                             ta_zzzzz_xxzzz_0,  \
                             ta_zzzzz_xxzzz_1,  \
                             ta_zzzzz_xyyy_0,   \
                             ta_zzzzz_xyyy_1,   \
                             ta_zzzzz_xyyyy_0,  \
                             ta_zzzzz_xyyyy_1,  \
                             ta_zzzzz_xyyyz_0,  \
                             ta_zzzzz_xyyyz_1,  \
                             ta_zzzzz_xyyz_0,   \
                             ta_zzzzz_xyyz_1,   \
                             ta_zzzzz_xyyzz_0,  \
                             ta_zzzzz_xyyzz_1,  \
                             ta_zzzzz_xyzz_0,   \
                             ta_zzzzz_xyzz_1,   \
                             ta_zzzzz_xyzzz_0,  \
                             ta_zzzzz_xyzzz_1,  \
                             ta_zzzzz_xzzz_0,   \
                             ta_zzzzz_xzzz_1,   \
                             ta_zzzzz_xzzzz_0,  \
                             ta_zzzzz_xzzzz_1,  \
                             ta_zzzzz_yyyy_0,   \
                             ta_zzzzz_yyyy_1,   \
                             ta_zzzzz_yyyyy_0,  \
                             ta_zzzzz_yyyyy_1,  \
                             ta_zzzzz_yyyyz_0,  \
                             ta_zzzzz_yyyyz_1,  \
                             ta_zzzzz_yyyz_0,   \
                             ta_zzzzz_yyyz_1,   \
                             ta_zzzzz_yyyzz_0,  \
                             ta_zzzzz_yyyzz_1,  \
                             ta_zzzzz_yyzz_0,   \
                             ta_zzzzz_yyzz_1,   \
                             ta_zzzzz_yyzzz_0,  \
                             ta_zzzzz_yyzzz_1,  \
                             ta_zzzzz_yzzz_0,   \
                             ta_zzzzz_yzzz_1,   \
                             ta_zzzzz_yzzzz_0,  \
                             ta_zzzzz_yzzzz_1,  \
                             ta_zzzzz_zzzz_0,   \
                             ta_zzzzz_zzzz_1,   \
                             ta_zzzzz_zzzzz_0,  \
                             ta_zzzzz_zzzzz_1,  \
                             ta_zzzzzz_xxxxx_0, \
                             ta_zzzzzz_xxxxy_0, \
                             ta_zzzzzz_xxxxz_0, \
                             ta_zzzzzz_xxxyy_0, \
                             ta_zzzzzz_xxxyz_0, \
                             ta_zzzzzz_xxxzz_0, \
                             ta_zzzzzz_xxyyy_0, \
                             ta_zzzzzz_xxyyz_0, \
                             ta_zzzzzz_xxyzz_0, \
                             ta_zzzzzz_xxzzz_0, \
                             ta_zzzzzz_xyyyy_0, \
                             ta_zzzzzz_xyyyz_0, \
                             ta_zzzzzz_xyyzz_0, \
                             ta_zzzzzz_xyzzz_0, \
                             ta_zzzzzz_xzzzz_0, \
                             ta_zzzzzz_yyyyy_0, \
                             ta_zzzzzz_yyyyz_0, \
                             ta_zzzzzz_yyyzz_0, \
                             ta_zzzzzz_yyzzz_0, \
                             ta_zzzzzz_yzzzz_0, \
                             ta_zzzzzz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzzz_xxxxx_0[i] =
            5.0 * ta_zzzz_xxxxx_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxx_1[i] * fe_0 + ta_zzzzz_xxxxx_0[i] * pa_z[i] - ta_zzzzz_xxxxx_1[i] * pc_z[i];

        ta_zzzzzz_xxxxy_0[i] =
            5.0 * ta_zzzz_xxxxy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxy_1[i] * fe_0 + ta_zzzzz_xxxxy_0[i] * pa_z[i] - ta_zzzzz_xxxxy_1[i] * pc_z[i];

        ta_zzzzzz_xxxxz_0[i] = 5.0 * ta_zzzz_xxxxz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxz_1[i] * fe_0 + ta_zzzzz_xxxx_0[i] * fe_0 -
                               ta_zzzzz_xxxx_1[i] * fe_0 + ta_zzzzz_xxxxz_0[i] * pa_z[i] - ta_zzzzz_xxxxz_1[i] * pc_z[i];

        ta_zzzzzz_xxxyy_0[i] =
            5.0 * ta_zzzz_xxxyy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxyy_1[i] * fe_0 + ta_zzzzz_xxxyy_0[i] * pa_z[i] - ta_zzzzz_xxxyy_1[i] * pc_z[i];

        ta_zzzzzz_xxxyz_0[i] = 5.0 * ta_zzzz_xxxyz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxyz_1[i] * fe_0 + ta_zzzzz_xxxy_0[i] * fe_0 -
                               ta_zzzzz_xxxy_1[i] * fe_0 + ta_zzzzz_xxxyz_0[i] * pa_z[i] - ta_zzzzz_xxxyz_1[i] * pc_z[i];

        ta_zzzzzz_xxxzz_0[i] = 5.0 * ta_zzzz_xxxzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xxxz_0[i] * fe_0 -
                               2.0 * ta_zzzzz_xxxz_1[i] * fe_0 + ta_zzzzz_xxxzz_0[i] * pa_z[i] - ta_zzzzz_xxxzz_1[i] * pc_z[i];

        ta_zzzzzz_xxyyy_0[i] =
            5.0 * ta_zzzz_xxyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xxyyy_1[i] * fe_0 + ta_zzzzz_xxyyy_0[i] * pa_z[i] - ta_zzzzz_xxyyy_1[i] * pc_z[i];

        ta_zzzzzz_xxyyz_0[i] = 5.0 * ta_zzzz_xxyyz_0[i] * fe_0 - 5.0 * ta_zzzz_xxyyz_1[i] * fe_0 + ta_zzzzz_xxyy_0[i] * fe_0 -
                               ta_zzzzz_xxyy_1[i] * fe_0 + ta_zzzzz_xxyyz_0[i] * pa_z[i] - ta_zzzzz_xxyyz_1[i] * pc_z[i];

        ta_zzzzzz_xxyzz_0[i] = 5.0 * ta_zzzz_xxyzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xxyz_0[i] * fe_0 -
                               2.0 * ta_zzzzz_xxyz_1[i] * fe_0 + ta_zzzzz_xxyzz_0[i] * pa_z[i] - ta_zzzzz_xxyzz_1[i] * pc_z[i];

        ta_zzzzzz_xxzzz_0[i] = 5.0 * ta_zzzz_xxzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xxzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_xxzz_0[i] * fe_0 -
                               3.0 * ta_zzzzz_xxzz_1[i] * fe_0 + ta_zzzzz_xxzzz_0[i] * pa_z[i] - ta_zzzzz_xxzzz_1[i] * pc_z[i];

        ta_zzzzzz_xyyyy_0[i] =
            5.0 * ta_zzzz_xyyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xyyyy_1[i] * fe_0 + ta_zzzzz_xyyyy_0[i] * pa_z[i] - ta_zzzzz_xyyyy_1[i] * pc_z[i];

        ta_zzzzzz_xyyyz_0[i] = 5.0 * ta_zzzz_xyyyz_0[i] * fe_0 - 5.0 * ta_zzzz_xyyyz_1[i] * fe_0 + ta_zzzzz_xyyy_0[i] * fe_0 -
                               ta_zzzzz_xyyy_1[i] * fe_0 + ta_zzzzz_xyyyz_0[i] * pa_z[i] - ta_zzzzz_xyyyz_1[i] * pc_z[i];

        ta_zzzzzz_xyyzz_0[i] = 5.0 * ta_zzzz_xyyzz_0[i] * fe_0 - 5.0 * ta_zzzz_xyyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xyyz_0[i] * fe_0 -
                               2.0 * ta_zzzzz_xyyz_1[i] * fe_0 + ta_zzzzz_xyyzz_0[i] * pa_z[i] - ta_zzzzz_xyyzz_1[i] * pc_z[i];

        ta_zzzzzz_xyzzz_0[i] = 5.0 * ta_zzzz_xyzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xyzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_xyzz_0[i] * fe_0 -
                               3.0 * ta_zzzzz_xyzz_1[i] * fe_0 + ta_zzzzz_xyzzz_0[i] * pa_z[i] - ta_zzzzz_xyzzz_1[i] * pc_z[i];

        ta_zzzzzz_xzzzz_0[i] = 5.0 * ta_zzzz_xzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_xzzzz_1[i] * fe_0 + 4.0 * ta_zzzzz_xzzz_0[i] * fe_0 -
                               4.0 * ta_zzzzz_xzzz_1[i] * fe_0 + ta_zzzzz_xzzzz_0[i] * pa_z[i] - ta_zzzzz_xzzzz_1[i] * pc_z[i];

        ta_zzzzzz_yyyyy_0[i] =
            5.0 * ta_zzzz_yyyyy_0[i] * fe_0 - 5.0 * ta_zzzz_yyyyy_1[i] * fe_0 + ta_zzzzz_yyyyy_0[i] * pa_z[i] - ta_zzzzz_yyyyy_1[i] * pc_z[i];

        ta_zzzzzz_yyyyz_0[i] = 5.0 * ta_zzzz_yyyyz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyyz_1[i] * fe_0 + ta_zzzzz_yyyy_0[i] * fe_0 -
                               ta_zzzzz_yyyy_1[i] * fe_0 + ta_zzzzz_yyyyz_0[i] * pa_z[i] - ta_zzzzz_yyyyz_1[i] * pc_z[i];

        ta_zzzzzz_yyyzz_0[i] = 5.0 * ta_zzzz_yyyzz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyzz_1[i] * fe_0 + 2.0 * ta_zzzzz_yyyz_0[i] * fe_0 -
                               2.0 * ta_zzzzz_yyyz_1[i] * fe_0 + ta_zzzzz_yyyzz_0[i] * pa_z[i] - ta_zzzzz_yyyzz_1[i] * pc_z[i];

        ta_zzzzzz_yyzzz_0[i] = 5.0 * ta_zzzz_yyzzz_0[i] * fe_0 - 5.0 * ta_zzzz_yyzzz_1[i] * fe_0 + 3.0 * ta_zzzzz_yyzz_0[i] * fe_0 -
                               3.0 * ta_zzzzz_yyzz_1[i] * fe_0 + ta_zzzzz_yyzzz_0[i] * pa_z[i] - ta_zzzzz_yyzzz_1[i] * pc_z[i];

        ta_zzzzzz_yzzzz_0[i] = 5.0 * ta_zzzz_yzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_yzzzz_1[i] * fe_0 + 4.0 * ta_zzzzz_yzzz_0[i] * fe_0 -
                               4.0 * ta_zzzzz_yzzz_1[i] * fe_0 + ta_zzzzz_yzzzz_0[i] * pa_z[i] - ta_zzzzz_yzzzz_1[i] * pc_z[i];

        ta_zzzzzz_zzzzz_0[i] = 5.0 * ta_zzzz_zzzzz_0[i] * fe_0 - 5.0 * ta_zzzz_zzzzz_1[i] * fe_0 + 5.0 * ta_zzzzz_zzzz_0[i] * fe_0 -
                               5.0 * ta_zzzzz_zzzz_1[i] * fe_0 + ta_zzzzz_zzzzz_0[i] * pa_z[i] - ta_zzzzz_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
