#include "ThreeCenterElectronRepulsionPrimRecISH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ish(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ish,
                                 size_t idx_eri_0_gsh,
                                 size_t idx_eri_1_gsh,
                                 size_t idx_eri_1_hsg,
                                 size_t idx_eri_1_hsh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : GSH

    auto g_xxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh);

    auto g_xxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 1);

    auto g_xxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 2);

    auto g_xxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 3);

    auto g_xxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 4);

    auto g_xxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 5);

    auto g_xxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 6);

    auto g_xxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 7);

    auto g_xxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 8);

    auto g_xxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 9);

    auto g_xxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 10);

    auto g_xxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 11);

    auto g_xxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 12);

    auto g_xxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 13);

    auto g_xxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 14);

    auto g_xxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 15);

    auto g_xxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 16);

    auto g_xxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 17);

    auto g_xxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 18);

    auto g_xxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 19);

    auto g_xxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 20);

    auto g_xxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 21);

    auto g_xxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 23);

    auto g_xxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 26);

    auto g_xxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 30);

    auto g_xxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 35);

    auto g_xxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 42);

    auto g_xxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 43);

    auto g_xxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 45);

    auto g_xxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 48);

    auto g_xxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 52);

    auto g_xxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 63);

    auto g_xxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 64);

    auto g_xxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 65);

    auto g_xxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 66);

    auto g_xxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 67);

    auto g_xxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 68);

    auto g_xxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 69);

    auto g_xxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 70);

    auto g_xxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 71);

    auto g_xxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 72);

    auto g_xxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 73);

    auto g_xxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 74);

    auto g_xxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 75);

    auto g_xxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 76);

    auto g_xxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 77);

    auto g_xxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 78);

    auto g_xxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 79);

    auto g_xxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 80);

    auto g_xxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 81);

    auto g_xxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 82);

    auto g_xxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 83);

    auto g_xxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 105);

    auto g_xxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 106);

    auto g_xxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 107);

    auto g_xxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 108);

    auto g_xxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 109);

    auto g_xxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 110);

    auto g_xxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 111);

    auto g_xxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 112);

    auto g_xxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 113);

    auto g_xxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 114);

    auto g_xxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 115);

    auto g_xxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 116);

    auto g_xxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 117);

    auto g_xxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 118);

    auto g_xxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 119);

    auto g_xxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 120);

    auto g_xxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 121);

    auto g_xxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 122);

    auto g_xxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 123);

    auto g_xxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 124);

    auto g_xxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 125);

    auto g_xyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 127);

    auto g_xyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 129);

    auto g_xyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 130);

    auto g_xyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 132);

    auto g_xyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 133);

    auto g_xyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 134);

    auto g_xyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 136);

    auto g_xyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 137);

    auto g_xyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 138);

    auto g_xyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 139);

    auto g_xyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 141);

    auto g_xyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 142);

    auto g_xyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 143);

    auto g_xyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 144);

    auto g_xyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 145);

    auto g_xyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 146);

    auto g_xzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 191);

    auto g_xzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 193);

    auto g_xzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 194);

    auto g_xzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 196);

    auto g_xzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 197);

    auto g_xzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 198);

    auto g_xzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 200);

    auto g_xzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 201);

    auto g_xzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 202);

    auto g_xzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 203);

    auto g_xzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 204);

    auto g_xzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 205);

    auto g_xzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 206);

    auto g_xzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 207);

    auto g_xzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 208);

    auto g_xzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 209);

    auto g_yyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 210);

    auto g_yyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 211);

    auto g_yyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 212);

    auto g_yyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 213);

    auto g_yyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 214);

    auto g_yyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 215);

    auto g_yyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 216);

    auto g_yyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 217);

    auto g_yyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 218);

    auto g_yyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 219);

    auto g_yyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 220);

    auto g_yyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 221);

    auto g_yyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 222);

    auto g_yyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 223);

    auto g_yyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 224);

    auto g_yyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 225);

    auto g_yyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 226);

    auto g_yyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 227);

    auto g_yyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 228);

    auto g_yyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 229);

    auto g_yyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 230);

    auto g_yyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 232);

    auto g_yyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 234);

    auto g_yyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 237);

    auto g_yyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 241);

    auto g_yyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 246);

    auto g_yyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 252);

    auto g_yyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 253);

    auto g_yyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 254);

    auto g_yyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 255);

    auto g_yyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 256);

    auto g_yyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 257);

    auto g_yyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 258);

    auto g_yyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 259);

    auto g_yyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 260);

    auto g_yyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 261);

    auto g_yyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 262);

    auto g_yyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 263);

    auto g_yyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 264);

    auto g_yyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 265);

    auto g_yyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 266);

    auto g_yyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 267);

    auto g_yyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 268);

    auto g_yyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 269);

    auto g_yyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 270);

    auto g_yyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 271);

    auto g_yyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 272);

    auto g_yzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 273);

    auto g_yzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 275);

    auto g_yzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 277);

    auto g_yzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 278);

    auto g_yzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 280);

    auto g_yzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 281);

    auto g_yzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 282);

    auto g_yzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 284);

    auto g_yzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 285);

    auto g_yzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 286);

    auto g_yzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 287);

    auto g_yzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 289);

    auto g_yzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 290);

    auto g_yzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 291);

    auto g_yzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 292);

    auto g_yzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 293);

    auto g_zzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 294);

    auto g_zzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 295);

    auto g_zzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 296);

    auto g_zzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 297);

    auto g_zzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 298);

    auto g_zzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 299);

    auto g_zzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 300);

    auto g_zzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 301);

    auto g_zzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 302);

    auto g_zzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 303);

    auto g_zzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 304);

    auto g_zzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 305);

    auto g_zzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 306);

    auto g_zzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 307);

    auto g_zzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 308);

    auto g_zzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 309);

    auto g_zzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 310);

    auto g_zzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 311);

    auto g_zzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 312);

    auto g_zzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 313);

    auto g_zzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 314);

    /// Set up components of auxilary buffer : GSH

    auto g_xxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh);

    auto g_xxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 1);

    auto g_xxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 2);

    auto g_xxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 3);

    auto g_xxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 4);

    auto g_xxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 5);

    auto g_xxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 6);

    auto g_xxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 7);

    auto g_xxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 8);

    auto g_xxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 9);

    auto g_xxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 10);

    auto g_xxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 11);

    auto g_xxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 12);

    auto g_xxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 13);

    auto g_xxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 14);

    auto g_xxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 15);

    auto g_xxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 16);

    auto g_xxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 17);

    auto g_xxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 18);

    auto g_xxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 19);

    auto g_xxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 20);

    auto g_xxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 21);

    auto g_xxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 23);

    auto g_xxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 26);

    auto g_xxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 30);

    auto g_xxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 35);

    auto g_xxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 42);

    auto g_xxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 43);

    auto g_xxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 45);

    auto g_xxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 48);

    auto g_xxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 52);

    auto g_xxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 63);

    auto g_xxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 64);

    auto g_xxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 65);

    auto g_xxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 66);

    auto g_xxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 67);

    auto g_xxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 68);

    auto g_xxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 69);

    auto g_xxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 70);

    auto g_xxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 71);

    auto g_xxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 72);

    auto g_xxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 73);

    auto g_xxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 74);

    auto g_xxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 75);

    auto g_xxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 76);

    auto g_xxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 77);

    auto g_xxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 78);

    auto g_xxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 79);

    auto g_xxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 80);

    auto g_xxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 81);

    auto g_xxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 82);

    auto g_xxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 83);

    auto g_xxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 105);

    auto g_xxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 106);

    auto g_xxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 107);

    auto g_xxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 108);

    auto g_xxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 109);

    auto g_xxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 110);

    auto g_xxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 111);

    auto g_xxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 112);

    auto g_xxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 113);

    auto g_xxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 114);

    auto g_xxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 115);

    auto g_xxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 116);

    auto g_xxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 117);

    auto g_xxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 118);

    auto g_xxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 119);

    auto g_xxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 120);

    auto g_xxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 121);

    auto g_xxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 122);

    auto g_xxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 123);

    auto g_xxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 124);

    auto g_xxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 125);

    auto g_xyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 127);

    auto g_xyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 129);

    auto g_xyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 130);

    auto g_xyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 132);

    auto g_xyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 133);

    auto g_xyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 134);

    auto g_xyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 136);

    auto g_xyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 137);

    auto g_xyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 138);

    auto g_xyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 139);

    auto g_xyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 141);

    auto g_xyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 142);

    auto g_xyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 143);

    auto g_xyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 144);

    auto g_xyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 145);

    auto g_xyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 146);

    auto g_xzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 191);

    auto g_xzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 193);

    auto g_xzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 194);

    auto g_xzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 196);

    auto g_xzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 197);

    auto g_xzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 198);

    auto g_xzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 200);

    auto g_xzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 201);

    auto g_xzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 202);

    auto g_xzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 203);

    auto g_xzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 204);

    auto g_xzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 205);

    auto g_xzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 206);

    auto g_xzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 207);

    auto g_xzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 208);

    auto g_xzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 209);

    auto g_yyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 210);

    auto g_yyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 211);

    auto g_yyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 212);

    auto g_yyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 213);

    auto g_yyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 214);

    auto g_yyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 215);

    auto g_yyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 216);

    auto g_yyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 217);

    auto g_yyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 218);

    auto g_yyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 219);

    auto g_yyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 220);

    auto g_yyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 221);

    auto g_yyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 222);

    auto g_yyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 223);

    auto g_yyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 224);

    auto g_yyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 225);

    auto g_yyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 226);

    auto g_yyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 227);

    auto g_yyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 228);

    auto g_yyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 229);

    auto g_yyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 230);

    auto g_yyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 232);

    auto g_yyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 234);

    auto g_yyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 237);

    auto g_yyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 241);

    auto g_yyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 246);

    auto g_yyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 252);

    auto g_yyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 253);

    auto g_yyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 254);

    auto g_yyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 255);

    auto g_yyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 256);

    auto g_yyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 257);

    auto g_yyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 258);

    auto g_yyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 259);

    auto g_yyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 260);

    auto g_yyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 261);

    auto g_yyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 262);

    auto g_yyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 263);

    auto g_yyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 264);

    auto g_yyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 265);

    auto g_yyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 266);

    auto g_yyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 267);

    auto g_yyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 268);

    auto g_yyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 269);

    auto g_yyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 270);

    auto g_yyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 271);

    auto g_yyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 272);

    auto g_yzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 273);

    auto g_yzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 275);

    auto g_yzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 277);

    auto g_yzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 278);

    auto g_yzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 280);

    auto g_yzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 281);

    auto g_yzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 282);

    auto g_yzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 284);

    auto g_yzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 285);

    auto g_yzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 286);

    auto g_yzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 287);

    auto g_yzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 289);

    auto g_yzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 290);

    auto g_yzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 291);

    auto g_yzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 292);

    auto g_yzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 293);

    auto g_zzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 294);

    auto g_zzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 295);

    auto g_zzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 296);

    auto g_zzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 297);

    auto g_zzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 298);

    auto g_zzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 299);

    auto g_zzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 300);

    auto g_zzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 301);

    auto g_zzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 302);

    auto g_zzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 303);

    auto g_zzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 304);

    auto g_zzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 305);

    auto g_zzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 306);

    auto g_zzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 307);

    auto g_zzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 308);

    auto g_zzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 309);

    auto g_zzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 310);

    auto g_zzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 311);

    auto g_zzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 312);

    auto g_zzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 313);

    auto g_zzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 314);

    /// Set up components of auxilary buffer : HSG

    auto g_xxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg);

    auto g_xxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 1);

    auto g_xxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 2);

    auto g_xxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 3);

    auto g_xxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 4);

    auto g_xxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 5);

    auto g_xxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 6);

    auto g_xxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 7);

    auto g_xxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 8);

    auto g_xxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 9);

    auto g_xxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 10);

    auto g_xxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 11);

    auto g_xxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 12);

    auto g_xxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 13);

    auto g_xxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 14);

    auto g_xxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 32);

    auto g_xxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 34);

    auto g_xxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 35);

    auto g_xxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 37);

    auto g_xxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 38);

    auto g_xxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 39);

    auto g_xxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 41);

    auto g_xxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 42);

    auto g_xxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 43);

    auto g_xxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 44);

    auto g_xxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 45);

    auto g_xxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 46);

    auto g_xxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 47);

    auto g_xxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 48);

    auto g_xxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 49);

    auto g_xxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 50);

    auto g_xxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 51);

    auto g_xxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 52);

    auto g_xxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 53);

    auto g_xxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 54);

    auto g_xxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 55);

    auto g_xxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 56);

    auto g_xxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 57);

    auto g_xxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 58);

    auto g_xxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 59);

    auto g_xxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 75);

    auto g_xxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 76);

    auto g_xxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 77);

    auto g_xxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 78);

    auto g_xxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 79);

    auto g_xxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 80);

    auto g_xxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 81);

    auto g_xxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 82);

    auto g_xxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 83);

    auto g_xxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 84);

    auto g_xxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 85);

    auto g_xxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 86);

    auto g_xxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 87);

    auto g_xxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 88);

    auto g_xxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 89);

    auto g_xxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 90);

    auto g_xxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 91);

    auto g_xxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 92);

    auto g_xxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 93);

    auto g_xxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 94);

    auto g_xxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 95);

    auto g_xxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 96);

    auto g_xxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 97);

    auto g_xxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 98);

    auto g_xxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 99);

    auto g_xxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 100);

    auto g_xxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 101);

    auto g_xxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 102);

    auto g_xxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 103);

    auto g_xxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 104);

    auto g_xxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 135);

    auto g_xxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 136);

    auto g_xxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 137);

    auto g_xxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 138);

    auto g_xxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 139);

    auto g_xxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 140);

    auto g_xxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 141);

    auto g_xxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 142);

    auto g_xxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 143);

    auto g_xxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 144);

    auto g_xxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 145);

    auto g_xxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 146);

    auto g_xxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 147);

    auto g_xxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 148);

    auto g_xxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 149);

    auto g_xyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 151);

    auto g_xyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 153);

    auto g_xyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 154);

    auto g_xyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 156);

    auto g_xyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 157);

    auto g_xyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 158);

    auto g_xyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 160);

    auto g_xyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 161);

    auto g_xyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 162);

    auto g_xyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 163);

    auto g_xyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 184);

    auto g_xyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 187);

    auto g_xyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 188);

    auto g_xyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 191);

    auto g_xyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 192);

    auto g_xyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 193);

    auto g_xzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 212);

    auto g_xzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 214);

    auto g_xzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 215);

    auto g_xzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 217);

    auto g_xzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 218);

    auto g_xzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 219);

    auto g_xzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 221);

    auto g_xzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 222);

    auto g_xzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 223);

    auto g_xzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 224);

    auto g_yyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 225);

    auto g_yyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 226);

    auto g_yyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 227);

    auto g_yyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 228);

    auto g_yyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 229);

    auto g_yyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 230);

    auto g_yyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 231);

    auto g_yyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 232);

    auto g_yyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 233);

    auto g_yyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 234);

    auto g_yyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 235);

    auto g_yyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 236);

    auto g_yyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 237);

    auto g_yyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 238);

    auto g_yyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 239);

    auto g_yyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 242);

    auto g_yyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 244);

    auto g_yyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 245);

    auto g_yyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 247);

    auto g_yyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 248);

    auto g_yyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 249);

    auto g_yyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 251);

    auto g_yyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 252);

    auto g_yyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 253);

    auto g_yyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 254);

    auto g_yyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 255);

    auto g_yyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 256);

    auto g_yyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 257);

    auto g_yyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 258);

    auto g_yyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 259);

    auto g_yyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 260);

    auto g_yyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 261);

    auto g_yyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 262);

    auto g_yyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 263);

    auto g_yyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 264);

    auto g_yyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 265);

    auto g_yyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 266);

    auto g_yyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 267);

    auto g_yyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 268);

    auto g_yyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 269);

    auto g_yyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 270);

    auto g_yyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 271);

    auto g_yyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 272);

    auto g_yyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 273);

    auto g_yyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 274);

    auto g_yyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 275);

    auto g_yyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 276);

    auto g_yyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 277);

    auto g_yyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 278);

    auto g_yyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 279);

    auto g_yyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 280);

    auto g_yyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 281);

    auto g_yyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 282);

    auto g_yyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 283);

    auto g_yyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 284);

    auto g_yzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 286);

    auto g_yzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 287);

    auto g_yzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 288);

    auto g_yzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 289);

    auto g_yzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 290);

    auto g_yzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 291);

    auto g_yzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 292);

    auto g_yzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 293);

    auto g_yzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 294);

    auto g_yzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 295);

    auto g_yzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 296);

    auto g_yzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 297);

    auto g_yzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 298);

    auto g_yzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 299);

    auto g_zzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 300);

    auto g_zzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 301);

    auto g_zzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 302);

    auto g_zzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 303);

    auto g_zzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 304);

    auto g_zzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 305);

    auto g_zzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 306);

    auto g_zzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 307);

    auto g_zzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 308);

    auto g_zzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 309);

    auto g_zzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 310);

    auto g_zzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 311);

    auto g_zzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 312);

    auto g_zzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 313);

    auto g_zzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 314);

    /// Set up components of auxilary buffer : HSH

    auto g_xxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh);

    auto g_xxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 1);

    auto g_xxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 2);

    auto g_xxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 3);

    auto g_xxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 4);

    auto g_xxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 5);

    auto g_xxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 6);

    auto g_xxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 7);

    auto g_xxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 8);

    auto g_xxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 9);

    auto g_xxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 10);

    auto g_xxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 11);

    auto g_xxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 12);

    auto g_xxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 13);

    auto g_xxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 14);

    auto g_xxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 15);

    auto g_xxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 16);

    auto g_xxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 17);

    auto g_xxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 18);

    auto g_xxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 19);

    auto g_xxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 20);

    auto g_xxxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 21);

    auto g_xxxxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 22);

    auto g_xxxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 23);

    auto g_xxxxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 24);

    auto g_xxxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 26);

    auto g_xxxxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 27);

    auto g_xxxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 30);

    auto g_xxxxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 31);

    auto g_xxxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 35);

    auto g_xxxxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 36);

    auto g_xxxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 42);

    auto g_xxxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 43);

    auto g_xxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 44);

    auto g_xxxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 45);

    auto g_xxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 46);

    auto g_xxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 47);

    auto g_xxxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 48);

    auto g_xxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 49);

    auto g_xxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 50);

    auto g_xxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 51);

    auto g_xxxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 52);

    auto g_xxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 53);

    auto g_xxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 54);

    auto g_xxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 55);

    auto g_xxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 56);

    auto g_xxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 58);

    auto g_xxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 59);

    auto g_xxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 60);

    auto g_xxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 61);

    auto g_xxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 62);

    auto g_xxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 63);

    auto g_xxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 64);

    auto g_xxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 65);

    auto g_xxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 66);

    auto g_xxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 67);

    auto g_xxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 68);

    auto g_xxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 69);

    auto g_xxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 70);

    auto g_xxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 71);

    auto g_xxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 72);

    auto g_xxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 73);

    auto g_xxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 74);

    auto g_xxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 75);

    auto g_xxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 76);

    auto g_xxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 77);

    auto g_xxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 78);

    auto g_xxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 79);

    auto g_xxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 80);

    auto g_xxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 81);

    auto g_xxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 82);

    auto g_xxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 83);

    auto g_xxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 105);

    auto g_xxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 106);

    auto g_xxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 107);

    auto g_xxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 108);

    auto g_xxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 109);

    auto g_xxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 110);

    auto g_xxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 111);

    auto g_xxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 112);

    auto g_xxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 113);

    auto g_xxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 114);

    auto g_xxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 115);

    auto g_xxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 116);

    auto g_xxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 117);

    auto g_xxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 118);

    auto g_xxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 119);

    auto g_xxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 120);

    auto g_xxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 121);

    auto g_xxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 122);

    auto g_xxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 123);

    auto g_xxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 124);

    auto g_xxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 125);

    auto g_xxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 126);

    auto g_xxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 127);

    auto g_xxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 128);

    auto g_xxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 129);

    auto g_xxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 130);

    auto g_xxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 131);

    auto g_xxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 132);

    auto g_xxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 133);

    auto g_xxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 134);

    auto g_xxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 135);

    auto g_xxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 136);

    auto g_xxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 137);

    auto g_xxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 138);

    auto g_xxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 139);

    auto g_xxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 140);

    auto g_xxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 141);

    auto g_xxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 142);

    auto g_xxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 143);

    auto g_xxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 144);

    auto g_xxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 145);

    auto g_xxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 146);

    auto g_xxyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 148);

    auto g_xxyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 150);

    auto g_xxyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 153);

    auto g_xxyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 157);

    auto g_xxyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 168);

    auto g_xxyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 170);

    auto g_xxyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 173);

    auto g_xxyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 177);

    auto g_xxyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 182);

    auto g_xxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 189);

    auto g_xxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 190);

    auto g_xxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 191);

    auto g_xxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 192);

    auto g_xxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 193);

    auto g_xxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 194);

    auto g_xxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 195);

    auto g_xxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 196);

    auto g_xxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 197);

    auto g_xxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 198);

    auto g_xxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 199);

    auto g_xxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 200);

    auto g_xxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 201);

    auto g_xxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 202);

    auto g_xxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 203);

    auto g_xxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 204);

    auto g_xxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 205);

    auto g_xxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 206);

    auto g_xxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 207);

    auto g_xxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 208);

    auto g_xxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 209);

    auto g_xyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 210);

    auto g_xyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 211);

    auto g_xyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 213);

    auto g_xyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 214);

    auto g_xyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 216);

    auto g_xyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 217);

    auto g_xyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 218);

    auto g_xyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 220);

    auto g_xyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 221);

    auto g_xyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 222);

    auto g_xyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 223);

    auto g_xyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 225);

    auto g_xyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 226);

    auto g_xyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 227);

    auto g_xyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 228);

    auto g_xyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 229);

    auto g_xyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 230);

    auto g_xyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 256);

    auto g_xyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 259);

    auto g_xyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 260);

    auto g_xyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 263);

    auto g_xyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 264);

    auto g_xyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 265);

    auto g_xyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 267);

    auto g_xyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 268);

    auto g_xyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 269);

    auto g_xyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 270);

    auto g_xyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 271);

    auto g_xyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 272);

    auto g_xzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 294);

    auto g_xzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 296);

    auto g_xzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 298);

    auto g_xzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 299);

    auto g_xzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 301);

    auto g_xzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 302);

    auto g_xzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 303);

    auto g_xzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 305);

    auto g_xzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 306);

    auto g_xzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 307);

    auto g_xzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 308);

    auto g_xzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 309);

    auto g_xzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 310);

    auto g_xzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 311);

    auto g_xzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 312);

    auto g_xzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 313);

    auto g_xzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 314);

    auto g_yyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 315);

    auto g_yyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 316);

    auto g_yyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 317);

    auto g_yyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 318);

    auto g_yyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 319);

    auto g_yyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 320);

    auto g_yyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 321);

    auto g_yyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 322);

    auto g_yyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 323);

    auto g_yyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 324);

    auto g_yyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 325);

    auto g_yyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 326);

    auto g_yyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 327);

    auto g_yyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 328);

    auto g_yyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 329);

    auto g_yyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 330);

    auto g_yyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 331);

    auto g_yyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 332);

    auto g_yyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 333);

    auto g_yyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 334);

    auto g_yyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 335);

    auto g_yyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 337);

    auto g_yyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 338);

    auto g_yyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 339);

    auto g_yyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 340);

    auto g_yyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 341);

    auto g_yyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 342);

    auto g_yyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 343);

    auto g_yyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 344);

    auto g_yyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 345);

    auto g_yyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 346);

    auto g_yyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 347);

    auto g_yyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 348);

    auto g_yyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 349);

    auto g_yyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 350);

    auto g_yyyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 351);

    auto g_yyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 352);

    auto g_yyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 353);

    auto g_yyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 354);

    auto g_yyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 355);

    auto g_yyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 356);

    auto g_yyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 357);

    auto g_yyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 358);

    auto g_yyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 359);

    auto g_yyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 360);

    auto g_yyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 361);

    auto g_yyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 362);

    auto g_yyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 363);

    auto g_yyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 364);

    auto g_yyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 365);

    auto g_yyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 366);

    auto g_yyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 367);

    auto g_yyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 368);

    auto g_yyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 369);

    auto g_yyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 370);

    auto g_yyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 371);

    auto g_yyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 372);

    auto g_yyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 373);

    auto g_yyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 374);

    auto g_yyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 375);

    auto g_yyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 376);

    auto g_yyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 377);

    auto g_yyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 378);

    auto g_yyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 379);

    auto g_yyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 380);

    auto g_yyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 381);

    auto g_yyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 382);

    auto g_yyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 383);

    auto g_yyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 384);

    auto g_yyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 385);

    auto g_yyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 386);

    auto g_yyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 387);

    auto g_yyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 388);

    auto g_yyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 389);

    auto g_yyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 390);

    auto g_yyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 391);

    auto g_yyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 392);

    auto g_yyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 393);

    auto g_yyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 394);

    auto g_yyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 395);

    auto g_yyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 396);

    auto g_yyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 397);

    auto g_yyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 398);

    auto g_yzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 399);

    auto g_yzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 400);

    auto g_yzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 401);

    auto g_yzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 402);

    auto g_yzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 403);

    auto g_yzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 404);

    auto g_yzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 405);

    auto g_yzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 406);

    auto g_yzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 407);

    auto g_yzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 408);

    auto g_yzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 409);

    auto g_yzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 410);

    auto g_yzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 411);

    auto g_yzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 412);

    auto g_yzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 413);

    auto g_yzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 414);

    auto g_yzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 415);

    auto g_yzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 416);

    auto g_yzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 417);

    auto g_yzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 418);

    auto g_yzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 419);

    auto g_zzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 420);

    auto g_zzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 421);

    auto g_zzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 422);

    auto g_zzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 423);

    auto g_zzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 424);

    auto g_zzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 425);

    auto g_zzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 426);

    auto g_zzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 427);

    auto g_zzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 428);

    auto g_zzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 429);

    auto g_zzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 430);

    auto g_zzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 431);

    auto g_zzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 432);

    auto g_zzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 433);

    auto g_zzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 434);

    auto g_zzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 435);

    auto g_zzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 436);

    auto g_zzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 437);

    auto g_zzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 438);

    auto g_zzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 439);

    auto g_zzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 440);

    /// Set up 0-21 components of targeted buffer : ISH

    auto g_xxxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish);

    auto g_xxxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 1);

    auto g_xxxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 2);

    auto g_xxxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 3);

    auto g_xxxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 4);

    auto g_xxxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 5);

    auto g_xxxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 6);

    auto g_xxxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 7);

    auto g_xxxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 8);

    auto g_xxxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 9);

    auto g_xxxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 10);

    auto g_xxxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 11);

    auto g_xxxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 12);

    auto g_xxxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 13);

    auto g_xxxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 14);

    auto g_xxxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 15);

    auto g_xxxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 16);

    auto g_xxxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 17);

    auto g_xxxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 18);

    auto g_xxxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 19);

    auto g_xxxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 20);

    #pragma omp simd aligned(g_xxxx_0_xxxxx_0, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxy_0, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxz_0, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxyy_0, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyz_0, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxzz_0, g_xxxx_0_xxxzz_1, g_xxxx_0_xxyyy_0, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyz_0, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyzz_0, g_xxxx_0_xxyzz_1, g_xxxx_0_xxzzz_0, g_xxxx_0_xxzzz_1, g_xxxx_0_xyyyy_0, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyz_0, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyzz_0, g_xxxx_0_xyyzz_1, g_xxxx_0_xyzzz_0, g_xxxx_0_xyzzz_1, g_xxxx_0_xzzzz_0, g_xxxx_0_xzzzz_1, g_xxxx_0_yyyyy_0, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyz_0, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyzz_0, g_xxxx_0_yyyzz_1, g_xxxx_0_yyzzz_0, g_xxxx_0_yyzzz_1, g_xxxx_0_yzzzz_0, g_xxxx_0_yzzzz_1, g_xxxx_0_zzzzz_0, g_xxxx_0_zzzzz_1, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxzz_1, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyzz_1, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xzzz_1, g_xxxxx_0_xzzzz_1, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyzz_1, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yzzz_1, g_xxxxx_0_yzzzz_1, g_xxxxx_0_zzzz_1, g_xxxxx_0_zzzzz_1, g_xxxxxx_0_xxxxx_0, g_xxxxxx_0_xxxxy_0, g_xxxxxx_0_xxxxz_0, g_xxxxxx_0_xxxyy_0, g_xxxxxx_0_xxxyz_0, g_xxxxxx_0_xxxzz_0, g_xxxxxx_0_xxyyy_0, g_xxxxxx_0_xxyyz_0, g_xxxxxx_0_xxyzz_0, g_xxxxxx_0_xxzzz_0, g_xxxxxx_0_xyyyy_0, g_xxxxxx_0_xyyyz_0, g_xxxxxx_0_xyyzz_0, g_xxxxxx_0_xyzzz_0, g_xxxxxx_0_xzzzz_0, g_xxxxxx_0_yyyyy_0, g_xxxxxx_0_yyyyz_0, g_xxxxxx_0_yyyzz_0, g_xxxxxx_0_yyzzz_0, g_xxxxxx_0_yzzzz_0, g_xxxxxx_0_zzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xxxxx_0[i] = 5.0 * g_xxxx_0_xxxxx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxx_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxy_0[i] = 5.0 * g_xxxx_0_xxxxy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxz_0[i] = 5.0 * g_xxxx_0_xxxxz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyy_0[i] = 5.0 * g_xxxx_0_xxxyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyz_0[i] = 5.0 * g_xxxx_0_xxxyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxzz_0[i] = 5.0 * g_xxxx_0_xxxzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyy_0[i] = 5.0 * g_xxxx_0_xxyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyz_0[i] = 5.0 * g_xxxx_0_xxyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyzz_0[i] = 5.0 * g_xxxx_0_xxyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxzzz_0[i] = 5.0 * g_xxxx_0_xxzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyy_0[i] = 5.0 * g_xxxx_0_xyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyz_0[i] = 5.0 * g_xxxx_0_xyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyzz_0[i] = 5.0 * g_xxxx_0_xyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyzzz_0[i] = 5.0 * g_xxxx_0_xyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xzzzz_0[i] = 5.0 * g_xxxx_0_xzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyy_0[i] = 5.0 * g_xxxx_0_yyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyz_0[i] = 5.0 * g_xxxx_0_yyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyzz_0[i] = 5.0 * g_xxxx_0_yyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyzzz_0[i] = 5.0 * g_xxxx_0_yyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yzzzz_0[i] = 5.0 * g_xxxx_0_yzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_zzzzz_0[i] = 5.0 * g_xxxx_0_zzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : ISH

    auto g_xxxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 21);

    auto g_xxxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 22);

    auto g_xxxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 23);

    auto g_xxxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 24);

    auto g_xxxxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 25);

    auto g_xxxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 26);

    auto g_xxxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 27);

    auto g_xxxxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 28);

    auto g_xxxxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 29);

    auto g_xxxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 30);

    auto g_xxxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 31);

    auto g_xxxxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 32);

    auto g_xxxxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 33);

    auto g_xxxxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 34);

    auto g_xxxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 35);

    auto g_xxxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 36);

    auto g_xxxxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 37);

    auto g_xxxxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 38);

    auto g_xxxxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 39);

    auto g_xxxxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 40);

    auto g_xxxxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 41);

    #pragma omp simd aligned(g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxzz_1, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyzz_1, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xzzz_1, g_xxxxx_0_xzzzz_1, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyzz_1, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yzzz_1, g_xxxxx_0_yzzzz_1, g_xxxxx_0_zzzz_1, g_xxxxx_0_zzzzz_1, g_xxxxxy_0_xxxxx_0, g_xxxxxy_0_xxxxy_0, g_xxxxxy_0_xxxxz_0, g_xxxxxy_0_xxxyy_0, g_xxxxxy_0_xxxyz_0, g_xxxxxy_0_xxxzz_0, g_xxxxxy_0_xxyyy_0, g_xxxxxy_0_xxyyz_0, g_xxxxxy_0_xxyzz_0, g_xxxxxy_0_xxzzz_0, g_xxxxxy_0_xyyyy_0, g_xxxxxy_0_xyyyz_0, g_xxxxxy_0_xyyzz_0, g_xxxxxy_0_xyzzz_0, g_xxxxxy_0_xzzzz_0, g_xxxxxy_0_yyyyy_0, g_xxxxxy_0_yyyyz_0, g_xxxxxy_0_yyyzz_0, g_xxxxxy_0_yyzzz_0, g_xxxxxy_0_yzzzz_0, g_xxxxxy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xxxxx_0[i] = g_xxxxx_0_xxxxx_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxy_0[i] = g_xxxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxz_0[i] = g_xxxxx_0_xxxxz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyy_0[i] = 2.0 * g_xxxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyz_0[i] = g_xxxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxzz_0[i] = g_xxxxx_0_xxxzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyy_0[i] = 3.0 * g_xxxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyz_0[i] = 2.0 * g_xxxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyzz_0[i] = g_xxxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxzzz_0[i] = g_xxxxx_0_xxzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyy_0[i] = 4.0 * g_xxxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyz_0[i] = 3.0 * g_xxxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyzz_0[i] = 2.0 * g_xxxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyzzz_0[i] = g_xxxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xzzzz_0[i] = g_xxxxx_0_xzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyy_0[i] = 5.0 * g_xxxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyz_0[i] = 4.0 * g_xxxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyzz_0[i] = 3.0 * g_xxxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyzzz_0[i] = 2.0 * g_xxxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yzzzz_0[i] = g_xxxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_zzzzz_0[i] = g_xxxxx_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 42-63 components of targeted buffer : ISH

    auto g_xxxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 42);

    auto g_xxxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 43);

    auto g_xxxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 44);

    auto g_xxxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 45);

    auto g_xxxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 46);

    auto g_xxxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 47);

    auto g_xxxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 48);

    auto g_xxxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 49);

    auto g_xxxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 50);

    auto g_xxxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 51);

    auto g_xxxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 52);

    auto g_xxxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 53);

    auto g_xxxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 54);

    auto g_xxxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 55);

    auto g_xxxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 56);

    auto g_xxxxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 57);

    auto g_xxxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 58);

    auto g_xxxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 59);

    auto g_xxxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 60);

    auto g_xxxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 61);

    auto g_xxxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 62);

    #pragma omp simd aligned(g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxzz_1, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyzz_1, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xzzz_1, g_xxxxx_0_xzzzz_1, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyzz_1, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yzzz_1, g_xxxxx_0_yzzzz_1, g_xxxxx_0_zzzz_1, g_xxxxx_0_zzzzz_1, g_xxxxxz_0_xxxxx_0, g_xxxxxz_0_xxxxy_0, g_xxxxxz_0_xxxxz_0, g_xxxxxz_0_xxxyy_0, g_xxxxxz_0_xxxyz_0, g_xxxxxz_0_xxxzz_0, g_xxxxxz_0_xxyyy_0, g_xxxxxz_0_xxyyz_0, g_xxxxxz_0_xxyzz_0, g_xxxxxz_0_xxzzz_0, g_xxxxxz_0_xyyyy_0, g_xxxxxz_0_xyyyz_0, g_xxxxxz_0_xyyzz_0, g_xxxxxz_0_xyzzz_0, g_xxxxxz_0_xzzzz_0, g_xxxxxz_0_yyyyy_0, g_xxxxxz_0_yyyyz_0, g_xxxxxz_0_yyyzz_0, g_xxxxxz_0_yyzzz_0, g_xxxxxz_0_yzzzz_0, g_xxxxxz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xxxxx_0[i] = g_xxxxx_0_xxxxx_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxy_0[i] = g_xxxxx_0_xxxxy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxz_0[i] = g_xxxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyy_0[i] = g_xxxxx_0_xxxyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyz_0[i] = g_xxxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxzz_0[i] = 2.0 * g_xxxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyy_0[i] = g_xxxxx_0_xxyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyz_0[i] = g_xxxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyzz_0[i] = 2.0 * g_xxxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxzzz_0[i] = 3.0 * g_xxxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyy_0[i] = g_xxxxx_0_xyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyz_0[i] = g_xxxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyzz_0[i] = 2.0 * g_xxxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyzzz_0[i] = 3.0 * g_xxxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xzzzz_0[i] = 4.0 * g_xxxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyy_0[i] = g_xxxxx_0_yyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyz_0[i] = g_xxxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyzz_0[i] = 2.0 * g_xxxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyzzz_0[i] = 3.0 * g_xxxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yzzzz_0[i] = 4.0 * g_xxxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_zzzzz_0[i] = 5.0 * g_xxxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxxx_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 63-84 components of targeted buffer : ISH

    auto g_xxxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 63);

    auto g_xxxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 64);

    auto g_xxxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 65);

    auto g_xxxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 66);

    auto g_xxxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 67);

    auto g_xxxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 68);

    auto g_xxxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 69);

    auto g_xxxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 70);

    auto g_xxxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 71);

    auto g_xxxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 72);

    auto g_xxxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 73);

    auto g_xxxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 74);

    auto g_xxxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 75);

    auto g_xxxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 76);

    auto g_xxxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 77);

    auto g_xxxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 78);

    auto g_xxxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 79);

    auto g_xxxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 80);

    auto g_xxxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 81);

    auto g_xxxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 82);

    auto g_xxxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 83);

    #pragma omp simd aligned(g_xxxx_0_xxxxx_0, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxz_0, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxzz_0, g_xxxx_0_xxxzz_1, g_xxxx_0_xxzzz_0, g_xxxx_0_xxzzz_1, g_xxxx_0_xzzzz_0, g_xxxx_0_xzzzz_1, g_xxxxy_0_xxxxx_1, g_xxxxy_0_xxxxz_1, g_xxxxy_0_xxxzz_1, g_xxxxy_0_xxzzz_1, g_xxxxy_0_xzzzz_1, g_xxxxyy_0_xxxxx_0, g_xxxxyy_0_xxxxy_0, g_xxxxyy_0_xxxxz_0, g_xxxxyy_0_xxxyy_0, g_xxxxyy_0_xxxyz_0, g_xxxxyy_0_xxxzz_0, g_xxxxyy_0_xxyyy_0, g_xxxxyy_0_xxyyz_0, g_xxxxyy_0_xxyzz_0, g_xxxxyy_0_xxzzz_0, g_xxxxyy_0_xyyyy_0, g_xxxxyy_0_xyyyz_0, g_xxxxyy_0_xyyzz_0, g_xxxxyy_0_xyzzz_0, g_xxxxyy_0_xzzzz_0, g_xxxxyy_0_yyyyy_0, g_xxxxyy_0_yyyyz_0, g_xxxxyy_0_yyyzz_0, g_xxxxyy_0_yyzzz_0, g_xxxxyy_0_yzzzz_0, g_xxxxyy_0_zzzzz_0, g_xxxyy_0_xxxxy_1, g_xxxyy_0_xxxy_1, g_xxxyy_0_xxxyy_1, g_xxxyy_0_xxxyz_1, g_xxxyy_0_xxyy_1, g_xxxyy_0_xxyyy_1, g_xxxyy_0_xxyyz_1, g_xxxyy_0_xxyz_1, g_xxxyy_0_xxyzz_1, g_xxxyy_0_xyyy_1, g_xxxyy_0_xyyyy_1, g_xxxyy_0_xyyyz_1, g_xxxyy_0_xyyz_1, g_xxxyy_0_xyyzz_1, g_xxxyy_0_xyzz_1, g_xxxyy_0_xyzzz_1, g_xxxyy_0_yyyy_1, g_xxxyy_0_yyyyy_1, g_xxxyy_0_yyyyz_1, g_xxxyy_0_yyyz_1, g_xxxyy_0_yyyzz_1, g_xxxyy_0_yyzz_1, g_xxxyy_0_yyzzz_1, g_xxxyy_0_yzzz_1, g_xxxyy_0_yzzzz_1, g_xxxyy_0_zzzzz_1, g_xxyy_0_xxxxy_0, g_xxyy_0_xxxxy_1, g_xxyy_0_xxxyy_0, g_xxyy_0_xxxyy_1, g_xxyy_0_xxxyz_0, g_xxyy_0_xxxyz_1, g_xxyy_0_xxyyy_0, g_xxyy_0_xxyyy_1, g_xxyy_0_xxyyz_0, g_xxyy_0_xxyyz_1, g_xxyy_0_xxyzz_0, g_xxyy_0_xxyzz_1, g_xxyy_0_xyyyy_0, g_xxyy_0_xyyyy_1, g_xxyy_0_xyyyz_0, g_xxyy_0_xyyyz_1, g_xxyy_0_xyyzz_0, g_xxyy_0_xyyzz_1, g_xxyy_0_xyzzz_0, g_xxyy_0_xyzzz_1, g_xxyy_0_yyyyy_0, g_xxyy_0_yyyyy_1, g_xxyy_0_yyyyz_0, g_xxyy_0_yyyyz_1, g_xxyy_0_yyyzz_0, g_xxyy_0_yyyzz_1, g_xxyy_0_yyzzz_0, g_xxyy_0_yyzzz_1, g_xxyy_0_yzzzz_0, g_xxyy_0_yzzzz_1, g_xxyy_0_zzzzz_0, g_xxyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xxxxx_0[i] = g_xxxx_0_xxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxx_1[i] * fz_be_0 + g_xxxxy_0_xxxxx_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxy_0[i] = 3.0 * g_xxyy_0_xxxxy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxz_0[i] = g_xxxx_0_xxxxz_0[i] * fbe_0 - g_xxxx_0_xxxxz_1[i] * fz_be_0 + g_xxxxy_0_xxxxz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxyy_0[i] = 3.0 * g_xxyy_0_xxxyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyz_0[i] = 3.0 * g_xxyy_0_xxxyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxzz_0[i] = g_xxxx_0_xxxzz_0[i] * fbe_0 - g_xxxx_0_xxxzz_1[i] * fz_be_0 + g_xxxxy_0_xxxzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxyyy_0[i] = 3.0 * g_xxyy_0_xxyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyz_0[i] = 3.0 * g_xxyy_0_xxyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyzz_0[i] = 3.0 * g_xxyy_0_xxyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxzzz_0[i] = g_xxxx_0_xxzzz_0[i] * fbe_0 - g_xxxx_0_xxzzz_1[i] * fz_be_0 + g_xxxxy_0_xxzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xyyyy_0[i] = 3.0 * g_xxyy_0_xyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyz_0[i] = 3.0 * g_xxyy_0_xyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyzz_0[i] = 3.0 * g_xxyy_0_xyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyzzz_0[i] = 3.0 * g_xxyy_0_xyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xzzzz_0[i] = g_xxxx_0_xzzzz_0[i] * fbe_0 - g_xxxx_0_xzzzz_1[i] * fz_be_0 + g_xxxxy_0_xzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_yyyyy_0[i] = 3.0 * g_xxyy_0_yyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyz_0[i] = 3.0 * g_xxyy_0_yyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyzz_0[i] = 3.0 * g_xxyy_0_yyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyzzz_0[i] = 3.0 * g_xxyy_0_yyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yzzzz_0[i] = 3.0 * g_xxyy_0_yzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_zzzzz_0[i] = 3.0 * g_xxyy_0_zzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zzzzz_1[i] * fz_be_0 + g_xxxyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 84-105 components of targeted buffer : ISH

    auto g_xxxxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 84);

    auto g_xxxxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 85);

    auto g_xxxxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 86);

    auto g_xxxxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 87);

    auto g_xxxxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 88);

    auto g_xxxxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 89);

    auto g_xxxxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 90);

    auto g_xxxxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 91);

    auto g_xxxxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 92);

    auto g_xxxxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 93);

    auto g_xxxxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 94);

    auto g_xxxxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 95);

    auto g_xxxxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 96);

    auto g_xxxxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 97);

    auto g_xxxxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 98);

    auto g_xxxxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 99);

    auto g_xxxxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 100);

    auto g_xxxxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 101);

    auto g_xxxxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 102);

    auto g_xxxxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 103);

    auto g_xxxxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 104);

    #pragma omp simd aligned(g_xxxxy_0_xxxxy_1, g_xxxxy_0_xxxyy_1, g_xxxxy_0_xxyyy_1, g_xxxxy_0_xyyyy_1, g_xxxxy_0_yyyyy_1, g_xxxxyz_0_xxxxx_0, g_xxxxyz_0_xxxxy_0, g_xxxxyz_0_xxxxz_0, g_xxxxyz_0_xxxyy_0, g_xxxxyz_0_xxxyz_0, g_xxxxyz_0_xxxzz_0, g_xxxxyz_0_xxyyy_0, g_xxxxyz_0_xxyyz_0, g_xxxxyz_0_xxyzz_0, g_xxxxyz_0_xxzzz_0, g_xxxxyz_0_xyyyy_0, g_xxxxyz_0_xyyyz_0, g_xxxxyz_0_xyyzz_0, g_xxxxyz_0_xyzzz_0, g_xxxxyz_0_xzzzz_0, g_xxxxyz_0_yyyyy_0, g_xxxxyz_0_yyyyz_0, g_xxxxyz_0_yyyzz_0, g_xxxxyz_0_yyzzz_0, g_xxxxyz_0_yzzzz_0, g_xxxxyz_0_zzzzz_0, g_xxxxz_0_xxxxx_1, g_xxxxz_0_xxxxz_1, g_xxxxz_0_xxxyz_1, g_xxxxz_0_xxxz_1, g_xxxxz_0_xxxzz_1, g_xxxxz_0_xxyyz_1, g_xxxxz_0_xxyz_1, g_xxxxz_0_xxyzz_1, g_xxxxz_0_xxzz_1, g_xxxxz_0_xxzzz_1, g_xxxxz_0_xyyyz_1, g_xxxxz_0_xyyz_1, g_xxxxz_0_xyyzz_1, g_xxxxz_0_xyzz_1, g_xxxxz_0_xyzzz_1, g_xxxxz_0_xzzz_1, g_xxxxz_0_xzzzz_1, g_xxxxz_0_yyyyz_1, g_xxxxz_0_yyyz_1, g_xxxxz_0_yyyzz_1, g_xxxxz_0_yyzz_1, g_xxxxz_0_yyzzz_1, g_xxxxz_0_yzzz_1, g_xxxxz_0_yzzzz_1, g_xxxxz_0_zzzz_1, g_xxxxz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xxxxx_0[i] = g_xxxxz_0_xxxxx_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxy_0[i] = g_xxxxy_0_xxxxy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxz_0[i] = g_xxxxz_0_xxxxz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyy_0[i] = g_xxxxy_0_xxxyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxyz_0[i] = g_xxxxz_0_xxxz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxzz_0[i] = g_xxxxz_0_xxxzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyy_0[i] = g_xxxxy_0_xxyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxyyz_0[i] = 2.0 * g_xxxxz_0_xxyz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyzz_0[i] = g_xxxxz_0_xxzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxzzz_0[i] = g_xxxxz_0_xxzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyy_0[i] = g_xxxxy_0_xyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xyyyz_0[i] = 3.0 * g_xxxxz_0_xyyz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyzz_0[i] = 2.0 * g_xxxxz_0_xyzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyzzz_0[i] = g_xxxxz_0_xzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xzzzz_0[i] = g_xxxxz_0_xzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyy_0[i] = g_xxxxy_0_yyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_yyyyz_0[i] = 4.0 * g_xxxxz_0_yyyz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyzz_0[i] = 3.0 * g_xxxxz_0_yyzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyzzz_0[i] = 2.0 * g_xxxxz_0_yzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yzzzz_0[i] = g_xxxxz_0_zzzz_1[i] * fi_acd_0 + g_xxxxz_0_yzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_zzzzz_0[i] = g_xxxxz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 105-126 components of targeted buffer : ISH

    auto g_xxxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 105);

    auto g_xxxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 106);

    auto g_xxxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 107);

    auto g_xxxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 108);

    auto g_xxxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 109);

    auto g_xxxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 110);

    auto g_xxxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 111);

    auto g_xxxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 112);

    auto g_xxxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 113);

    auto g_xxxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 114);

    auto g_xxxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 115);

    auto g_xxxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 116);

    auto g_xxxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 117);

    auto g_xxxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 118);

    auto g_xxxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 119);

    auto g_xxxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 120);

    auto g_xxxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 121);

    auto g_xxxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 122);

    auto g_xxxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 123);

    auto g_xxxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 124);

    auto g_xxxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 125);

    #pragma omp simd aligned(g_xxxx_0_xxxxx_0, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxy_0, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxyy_0, g_xxxx_0_xxxyy_1, g_xxxx_0_xxyyy_0, g_xxxx_0_xxyyy_1, g_xxxx_0_xyyyy_0, g_xxxx_0_xyyyy_1, g_xxxxz_0_xxxxx_1, g_xxxxz_0_xxxxy_1, g_xxxxz_0_xxxyy_1, g_xxxxz_0_xxyyy_1, g_xxxxz_0_xyyyy_1, g_xxxxzz_0_xxxxx_0, g_xxxxzz_0_xxxxy_0, g_xxxxzz_0_xxxxz_0, g_xxxxzz_0_xxxyy_0, g_xxxxzz_0_xxxyz_0, g_xxxxzz_0_xxxzz_0, g_xxxxzz_0_xxyyy_0, g_xxxxzz_0_xxyyz_0, g_xxxxzz_0_xxyzz_0, g_xxxxzz_0_xxzzz_0, g_xxxxzz_0_xyyyy_0, g_xxxxzz_0_xyyyz_0, g_xxxxzz_0_xyyzz_0, g_xxxxzz_0_xyzzz_0, g_xxxxzz_0_xzzzz_0, g_xxxxzz_0_yyyyy_0, g_xxxxzz_0_yyyyz_0, g_xxxxzz_0_yyyzz_0, g_xxxxzz_0_yyzzz_0, g_xxxxzz_0_yzzzz_0, g_xxxxzz_0_zzzzz_0, g_xxxzz_0_xxxxz_1, g_xxxzz_0_xxxyz_1, g_xxxzz_0_xxxz_1, g_xxxzz_0_xxxzz_1, g_xxxzz_0_xxyyz_1, g_xxxzz_0_xxyz_1, g_xxxzz_0_xxyzz_1, g_xxxzz_0_xxzz_1, g_xxxzz_0_xxzzz_1, g_xxxzz_0_xyyyz_1, g_xxxzz_0_xyyz_1, g_xxxzz_0_xyyzz_1, g_xxxzz_0_xyzz_1, g_xxxzz_0_xyzzz_1, g_xxxzz_0_xzzz_1, g_xxxzz_0_xzzzz_1, g_xxxzz_0_yyyyy_1, g_xxxzz_0_yyyyz_1, g_xxxzz_0_yyyz_1, g_xxxzz_0_yyyzz_1, g_xxxzz_0_yyzz_1, g_xxxzz_0_yyzzz_1, g_xxxzz_0_yzzz_1, g_xxxzz_0_yzzzz_1, g_xxxzz_0_zzzz_1, g_xxxzz_0_zzzzz_1, g_xxzz_0_xxxxz_0, g_xxzz_0_xxxxz_1, g_xxzz_0_xxxyz_0, g_xxzz_0_xxxyz_1, g_xxzz_0_xxxzz_0, g_xxzz_0_xxxzz_1, g_xxzz_0_xxyyz_0, g_xxzz_0_xxyyz_1, g_xxzz_0_xxyzz_0, g_xxzz_0_xxyzz_1, g_xxzz_0_xxzzz_0, g_xxzz_0_xxzzz_1, g_xxzz_0_xyyyz_0, g_xxzz_0_xyyyz_1, g_xxzz_0_xyyzz_0, g_xxzz_0_xyyzz_1, g_xxzz_0_xyzzz_0, g_xxzz_0_xyzzz_1, g_xxzz_0_xzzzz_0, g_xxzz_0_xzzzz_1, g_xxzz_0_yyyyy_0, g_xxzz_0_yyyyy_1, g_xxzz_0_yyyyz_0, g_xxzz_0_yyyyz_1, g_xxzz_0_yyyzz_0, g_xxzz_0_yyyzz_1, g_xxzz_0_yyzzz_0, g_xxzz_0_yyzzz_1, g_xxzz_0_yzzzz_0, g_xxzz_0_yzzzz_1, g_xxzz_0_zzzzz_0, g_xxzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xxxxx_0[i] = g_xxxx_0_xxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxx_1[i] * fz_be_0 + g_xxxxz_0_xxxxx_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxy_0[i] = g_xxxx_0_xxxxy_0[i] * fbe_0 - g_xxxx_0_xxxxy_1[i] * fz_be_0 + g_xxxxz_0_xxxxy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxz_0[i] = 3.0 * g_xxzz_0_xxxxz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyy_0[i] = g_xxxx_0_xxxyy_0[i] * fbe_0 - g_xxxx_0_xxxyy_1[i] * fz_be_0 + g_xxxxz_0_xxxyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxyz_0[i] = 3.0 * g_xxzz_0_xxxyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxzz_0[i] = 3.0 * g_xxzz_0_xxxzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyy_0[i] = g_xxxx_0_xxyyy_0[i] * fbe_0 - g_xxxx_0_xxyyy_1[i] * fz_be_0 + g_xxxxz_0_xxyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxyyz_0[i] = 3.0 * g_xxzz_0_xxyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyzz_0[i] = 3.0 * g_xxzz_0_xxyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxzzz_0[i] = 3.0 * g_xxzz_0_xxzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyy_0[i] = g_xxxx_0_xyyyy_0[i] * fbe_0 - g_xxxx_0_xyyyy_1[i] * fz_be_0 + g_xxxxz_0_xyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xyyyz_0[i] = 3.0 * g_xxzz_0_xyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyzz_0[i] = 3.0 * g_xxzz_0_xyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyzzz_0[i] = 3.0 * g_xxzz_0_xyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xzzzz_0[i] = 3.0 * g_xxzz_0_xzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzz_1[i] * fi_acd_0 + g_xxxzz_0_xzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyy_0[i] = 3.0 * g_xxzz_0_yyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyy_1[i] * fz_be_0 + g_xxxzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyz_0[i] = 3.0 * g_xxzz_0_yyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyzz_0[i] = 3.0 * g_xxzz_0_yyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyzzz_0[i] = 3.0 * g_xxzz_0_yyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yzzzz_0[i] = 3.0 * g_xxzz_0_yzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_zzzzz_0[i] = 3.0 * g_xxzz_0_zzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 126-147 components of targeted buffer : ISH

    auto g_xxxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 126);

    auto g_xxxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 127);

    auto g_xxxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 128);

    auto g_xxxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 129);

    auto g_xxxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 130);

    auto g_xxxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 131);

    auto g_xxxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 132);

    auto g_xxxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 133);

    auto g_xxxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 134);

    auto g_xxxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 135);

    auto g_xxxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 136);

    auto g_xxxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 137);

    auto g_xxxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 138);

    auto g_xxxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 139);

    auto g_xxxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 140);

    auto g_xxxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 141);

    auto g_xxxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 142);

    auto g_xxxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 143);

    auto g_xxxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 144);

    auto g_xxxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 145);

    auto g_xxxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 146);

    #pragma omp simd aligned(g_xxxy_0_xxxxx_0, g_xxxy_0_xxxxx_1, g_xxxy_0_xxxxz_0, g_xxxy_0_xxxxz_1, g_xxxy_0_xxxzz_0, g_xxxy_0_xxxzz_1, g_xxxy_0_xxzzz_0, g_xxxy_0_xxzzz_1, g_xxxy_0_xzzzz_0, g_xxxy_0_xzzzz_1, g_xxxyy_0_xxxxx_1, g_xxxyy_0_xxxxz_1, g_xxxyy_0_xxxzz_1, g_xxxyy_0_xxzzz_1, g_xxxyy_0_xzzzz_1, g_xxxyyy_0_xxxxx_0, g_xxxyyy_0_xxxxy_0, g_xxxyyy_0_xxxxz_0, g_xxxyyy_0_xxxyy_0, g_xxxyyy_0_xxxyz_0, g_xxxyyy_0_xxxzz_0, g_xxxyyy_0_xxyyy_0, g_xxxyyy_0_xxyyz_0, g_xxxyyy_0_xxyzz_0, g_xxxyyy_0_xxzzz_0, g_xxxyyy_0_xyyyy_0, g_xxxyyy_0_xyyyz_0, g_xxxyyy_0_xyyzz_0, g_xxxyyy_0_xyzzz_0, g_xxxyyy_0_xzzzz_0, g_xxxyyy_0_yyyyy_0, g_xxxyyy_0_yyyyz_0, g_xxxyyy_0_yyyzz_0, g_xxxyyy_0_yyzzz_0, g_xxxyyy_0_yzzzz_0, g_xxxyyy_0_zzzzz_0, g_xxyyy_0_xxxxy_1, g_xxyyy_0_xxxy_1, g_xxyyy_0_xxxyy_1, g_xxyyy_0_xxxyz_1, g_xxyyy_0_xxyy_1, g_xxyyy_0_xxyyy_1, g_xxyyy_0_xxyyz_1, g_xxyyy_0_xxyz_1, g_xxyyy_0_xxyzz_1, g_xxyyy_0_xyyy_1, g_xxyyy_0_xyyyy_1, g_xxyyy_0_xyyyz_1, g_xxyyy_0_xyyz_1, g_xxyyy_0_xyyzz_1, g_xxyyy_0_xyzz_1, g_xxyyy_0_xyzzz_1, g_xxyyy_0_yyyy_1, g_xxyyy_0_yyyyy_1, g_xxyyy_0_yyyyz_1, g_xxyyy_0_yyyz_1, g_xxyyy_0_yyyzz_1, g_xxyyy_0_yyzz_1, g_xxyyy_0_yyzzz_1, g_xxyyy_0_yzzz_1, g_xxyyy_0_yzzzz_1, g_xxyyy_0_zzzzz_1, g_xyyy_0_xxxxy_0, g_xyyy_0_xxxxy_1, g_xyyy_0_xxxyy_0, g_xyyy_0_xxxyy_1, g_xyyy_0_xxxyz_0, g_xyyy_0_xxxyz_1, g_xyyy_0_xxyyy_0, g_xyyy_0_xxyyy_1, g_xyyy_0_xxyyz_0, g_xyyy_0_xxyyz_1, g_xyyy_0_xxyzz_0, g_xyyy_0_xxyzz_1, g_xyyy_0_xyyyy_0, g_xyyy_0_xyyyy_1, g_xyyy_0_xyyyz_0, g_xyyy_0_xyyyz_1, g_xyyy_0_xyyzz_0, g_xyyy_0_xyyzz_1, g_xyyy_0_xyzzz_0, g_xyyy_0_xyzzz_1, g_xyyy_0_yyyyy_0, g_xyyy_0_yyyyy_1, g_xyyy_0_yyyyz_0, g_xyyy_0_yyyyz_1, g_xyyy_0_yyyzz_0, g_xyyy_0_yyyzz_1, g_xyyy_0_yyzzz_0, g_xyyy_0_yyzzz_1, g_xyyy_0_yzzzz_0, g_xyyy_0_yzzzz_1, g_xyyy_0_zzzzz_0, g_xyyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xxxxx_0[i] = 2.0 * g_xxxy_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxx_1[i] * fz_be_0 + g_xxxyy_0_xxxxx_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxy_0[i] = 2.0 * g_xyyy_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxz_0[i] = 2.0 * g_xxxy_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxz_1[i] * fz_be_0 + g_xxxyy_0_xxxxz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxyy_0[i] = 2.0 * g_xyyy_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyz_0[i] = 2.0 * g_xyyy_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxzz_0[i] = 2.0 * g_xxxy_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxzz_1[i] * fz_be_0 + g_xxxyy_0_xxxzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxyyy_0[i] = 2.0 * g_xyyy_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyz_0[i] = 2.0 * g_xyyy_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyzz_0[i] = 2.0 * g_xyyy_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxzzz_0[i] = 2.0 * g_xxxy_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxzzz_1[i] * fz_be_0 + g_xxxyy_0_xxzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xyyyy_0[i] = 2.0 * g_xyyy_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyz_0[i] = 2.0 * g_xyyy_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyzz_0[i] = 2.0 * g_xyyy_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyzzz_0[i] = 2.0 * g_xyyy_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xzzzz_0[i] = 2.0 * g_xxxy_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xzzzz_1[i] * fz_be_0 + g_xxxyy_0_xzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_yyyyy_0[i] = 2.0 * g_xyyy_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyz_0[i] = 2.0 * g_xyyy_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyzz_0[i] = 2.0 * g_xyyy_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyzzz_0[i] = 2.0 * g_xyyy_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yzzzz_0[i] = 2.0 * g_xyyy_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_zzzzz_0[i] = 2.0 * g_xyyy_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zzzzz_1[i] * fz_be_0 + g_xxyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 147-168 components of targeted buffer : ISH

    auto g_xxxyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 147);

    auto g_xxxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 148);

    auto g_xxxyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 149);

    auto g_xxxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 150);

    auto g_xxxyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 151);

    auto g_xxxyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 152);

    auto g_xxxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 153);

    auto g_xxxyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 154);

    auto g_xxxyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 155);

    auto g_xxxyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 156);

    auto g_xxxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 157);

    auto g_xxxyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 158);

    auto g_xxxyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 159);

    auto g_xxxyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 160);

    auto g_xxxyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 161);

    auto g_xxxyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 162);

    auto g_xxxyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 163);

    auto g_xxxyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 164);

    auto g_xxxyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 165);

    auto g_xxxyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 166);

    auto g_xxxyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 167);

    #pragma omp simd aligned(g_xxxyy_0_xxxx_1, g_xxxyy_0_xxxxx_1, g_xxxyy_0_xxxxy_1, g_xxxyy_0_xxxxz_1, g_xxxyy_0_xxxy_1, g_xxxyy_0_xxxyy_1, g_xxxyy_0_xxxyz_1, g_xxxyy_0_xxxz_1, g_xxxyy_0_xxxzz_1, g_xxxyy_0_xxyy_1, g_xxxyy_0_xxyyy_1, g_xxxyy_0_xxyyz_1, g_xxxyy_0_xxyz_1, g_xxxyy_0_xxyzz_1, g_xxxyy_0_xxzz_1, g_xxxyy_0_xxzzz_1, g_xxxyy_0_xyyy_1, g_xxxyy_0_xyyyy_1, g_xxxyy_0_xyyyz_1, g_xxxyy_0_xyyz_1, g_xxxyy_0_xyyzz_1, g_xxxyy_0_xyzz_1, g_xxxyy_0_xyzzz_1, g_xxxyy_0_xzzz_1, g_xxxyy_0_xzzzz_1, g_xxxyy_0_yyyy_1, g_xxxyy_0_yyyyy_1, g_xxxyy_0_yyyyz_1, g_xxxyy_0_yyyz_1, g_xxxyy_0_yyyzz_1, g_xxxyy_0_yyzz_1, g_xxxyy_0_yyzzz_1, g_xxxyy_0_yzzz_1, g_xxxyy_0_yzzzz_1, g_xxxyy_0_zzzz_1, g_xxxyy_0_zzzzz_1, g_xxxyyz_0_xxxxx_0, g_xxxyyz_0_xxxxy_0, g_xxxyyz_0_xxxxz_0, g_xxxyyz_0_xxxyy_0, g_xxxyyz_0_xxxyz_0, g_xxxyyz_0_xxxzz_0, g_xxxyyz_0_xxyyy_0, g_xxxyyz_0_xxyyz_0, g_xxxyyz_0_xxyzz_0, g_xxxyyz_0_xxzzz_0, g_xxxyyz_0_xyyyy_0, g_xxxyyz_0_xyyyz_0, g_xxxyyz_0_xyyzz_0, g_xxxyyz_0_xyzzz_0, g_xxxyyz_0_xzzzz_0, g_xxxyyz_0_yyyyy_0, g_xxxyyz_0_yyyyz_0, g_xxxyyz_0_yyyzz_0, g_xxxyyz_0_yyzzz_0, g_xxxyyz_0_yzzzz_0, g_xxxyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xxxxx_0[i] = g_xxxyy_0_xxxxx_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxy_0[i] = g_xxxyy_0_xxxxy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxz_0[i] = g_xxxyy_0_xxxx_1[i] * fi_acd_0 + g_xxxyy_0_xxxxz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyy_0[i] = g_xxxyy_0_xxxyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyz_0[i] = g_xxxyy_0_xxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxzz_0[i] = 2.0 * g_xxxyy_0_xxxz_1[i] * fi_acd_0 + g_xxxyy_0_xxxzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyy_0[i] = g_xxxyy_0_xxyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyz_0[i] = g_xxxyy_0_xxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyzz_0[i] = 2.0 * g_xxxyy_0_xxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxzzz_0[i] = 3.0 * g_xxxyy_0_xxzz_1[i] * fi_acd_0 + g_xxxyy_0_xxzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyy_0[i] = g_xxxyy_0_xyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyz_0[i] = g_xxxyy_0_xyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyzz_0[i] = 2.0 * g_xxxyy_0_xyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyzzz_0[i] = 3.0 * g_xxxyy_0_xyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xzzzz_0[i] = 4.0 * g_xxxyy_0_xzzz_1[i] * fi_acd_0 + g_xxxyy_0_xzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyy_0[i] = g_xxxyy_0_yyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyz_0[i] = g_xxxyy_0_yyyy_1[i] * fi_acd_0 + g_xxxyy_0_yyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyzz_0[i] = 2.0 * g_xxxyy_0_yyyz_1[i] * fi_acd_0 + g_xxxyy_0_yyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyzzz_0[i] = 3.0 * g_xxxyy_0_yyzz_1[i] * fi_acd_0 + g_xxxyy_0_yyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yzzzz_0[i] = 4.0 * g_xxxyy_0_yzzz_1[i] * fi_acd_0 + g_xxxyy_0_yzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_zzzzz_0[i] = 5.0 * g_xxxyy_0_zzzz_1[i] * fi_acd_0 + g_xxxyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 168-189 components of targeted buffer : ISH

    auto g_xxxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 168);

    auto g_xxxyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 169);

    auto g_xxxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 170);

    auto g_xxxyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 171);

    auto g_xxxyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 172);

    auto g_xxxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 173);

    auto g_xxxyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 174);

    auto g_xxxyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 175);

    auto g_xxxyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 176);

    auto g_xxxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 177);

    auto g_xxxyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 178);

    auto g_xxxyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 179);

    auto g_xxxyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 180);

    auto g_xxxyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 181);

    auto g_xxxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 182);

    auto g_xxxyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 183);

    auto g_xxxyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 184);

    auto g_xxxyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 185);

    auto g_xxxyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 186);

    auto g_xxxyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 187);

    auto g_xxxyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 188);

    #pragma omp simd aligned(g_xxxyzz_0_xxxxx_0, g_xxxyzz_0_xxxxy_0, g_xxxyzz_0_xxxxz_0, g_xxxyzz_0_xxxyy_0, g_xxxyzz_0_xxxyz_0, g_xxxyzz_0_xxxzz_0, g_xxxyzz_0_xxyyy_0, g_xxxyzz_0_xxyyz_0, g_xxxyzz_0_xxyzz_0, g_xxxyzz_0_xxzzz_0, g_xxxyzz_0_xyyyy_0, g_xxxyzz_0_xyyyz_0, g_xxxyzz_0_xyyzz_0, g_xxxyzz_0_xyzzz_0, g_xxxyzz_0_xzzzz_0, g_xxxyzz_0_yyyyy_0, g_xxxyzz_0_yyyyz_0, g_xxxyzz_0_yyyzz_0, g_xxxyzz_0_yyzzz_0, g_xxxyzz_0_yzzzz_0, g_xxxyzz_0_zzzzz_0, g_xxxzz_0_xxxx_1, g_xxxzz_0_xxxxx_1, g_xxxzz_0_xxxxy_1, g_xxxzz_0_xxxxz_1, g_xxxzz_0_xxxy_1, g_xxxzz_0_xxxyy_1, g_xxxzz_0_xxxyz_1, g_xxxzz_0_xxxz_1, g_xxxzz_0_xxxzz_1, g_xxxzz_0_xxyy_1, g_xxxzz_0_xxyyy_1, g_xxxzz_0_xxyyz_1, g_xxxzz_0_xxyz_1, g_xxxzz_0_xxyzz_1, g_xxxzz_0_xxzz_1, g_xxxzz_0_xxzzz_1, g_xxxzz_0_xyyy_1, g_xxxzz_0_xyyyy_1, g_xxxzz_0_xyyyz_1, g_xxxzz_0_xyyz_1, g_xxxzz_0_xyyzz_1, g_xxxzz_0_xyzz_1, g_xxxzz_0_xyzzz_1, g_xxxzz_0_xzzz_1, g_xxxzz_0_xzzzz_1, g_xxxzz_0_yyyy_1, g_xxxzz_0_yyyyy_1, g_xxxzz_0_yyyyz_1, g_xxxzz_0_yyyz_1, g_xxxzz_0_yyyzz_1, g_xxxzz_0_yyzz_1, g_xxxzz_0_yyzzz_1, g_xxxzz_0_yzzz_1, g_xxxzz_0_yzzzz_1, g_xxxzz_0_zzzz_1, g_xxxzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xxxxx_0[i] = g_xxxzz_0_xxxxx_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxy_0[i] = g_xxxzz_0_xxxx_1[i] * fi_acd_0 + g_xxxzz_0_xxxxy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxz_0[i] = g_xxxzz_0_xxxxz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyy_0[i] = 2.0 * g_xxxzz_0_xxxy_1[i] * fi_acd_0 + g_xxxzz_0_xxxyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyz_0[i] = g_xxxzz_0_xxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxzz_0[i] = g_xxxzz_0_xxxzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyy_0[i] = 3.0 * g_xxxzz_0_xxyy_1[i] * fi_acd_0 + g_xxxzz_0_xxyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyz_0[i] = 2.0 * g_xxxzz_0_xxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyzz_0[i] = g_xxxzz_0_xxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxzzz_0[i] = g_xxxzz_0_xxzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyy_0[i] = 4.0 * g_xxxzz_0_xyyy_1[i] * fi_acd_0 + g_xxxzz_0_xyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyz_0[i] = 3.0 * g_xxxzz_0_xyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyzz_0[i] = 2.0 * g_xxxzz_0_xyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyzzz_0[i] = g_xxxzz_0_xzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xzzzz_0[i] = g_xxxzz_0_xzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyy_0[i] = 5.0 * g_xxxzz_0_yyyy_1[i] * fi_acd_0 + g_xxxzz_0_yyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyz_0[i] = 4.0 * g_xxxzz_0_yyyz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyzz_0[i] = 3.0 * g_xxxzz_0_yyzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyzzz_0[i] = 2.0 * g_xxxzz_0_yzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yzzzz_0[i] = g_xxxzz_0_zzzz_1[i] * fi_acd_0 + g_xxxzz_0_yzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_zzzzz_0[i] = g_xxxzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 189-210 components of targeted buffer : ISH

    auto g_xxxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 189);

    auto g_xxxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 190);

    auto g_xxxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 191);

    auto g_xxxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 192);

    auto g_xxxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 193);

    auto g_xxxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 194);

    auto g_xxxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 195);

    auto g_xxxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 196);

    auto g_xxxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 197);

    auto g_xxxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 198);

    auto g_xxxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 199);

    auto g_xxxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 200);

    auto g_xxxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 201);

    auto g_xxxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 202);

    auto g_xxxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 203);

    auto g_xxxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 204);

    auto g_xxxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 205);

    auto g_xxxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 206);

    auto g_xxxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 207);

    auto g_xxxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 208);

    auto g_xxxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 209);

    #pragma omp simd aligned(g_xxxz_0_xxxxx_0, g_xxxz_0_xxxxx_1, g_xxxz_0_xxxxy_0, g_xxxz_0_xxxxy_1, g_xxxz_0_xxxyy_0, g_xxxz_0_xxxyy_1, g_xxxz_0_xxyyy_0, g_xxxz_0_xxyyy_1, g_xxxz_0_xyyyy_0, g_xxxz_0_xyyyy_1, g_xxxzz_0_xxxxx_1, g_xxxzz_0_xxxxy_1, g_xxxzz_0_xxxyy_1, g_xxxzz_0_xxyyy_1, g_xxxzz_0_xyyyy_1, g_xxxzzz_0_xxxxx_0, g_xxxzzz_0_xxxxy_0, g_xxxzzz_0_xxxxz_0, g_xxxzzz_0_xxxyy_0, g_xxxzzz_0_xxxyz_0, g_xxxzzz_0_xxxzz_0, g_xxxzzz_0_xxyyy_0, g_xxxzzz_0_xxyyz_0, g_xxxzzz_0_xxyzz_0, g_xxxzzz_0_xxzzz_0, g_xxxzzz_0_xyyyy_0, g_xxxzzz_0_xyyyz_0, g_xxxzzz_0_xyyzz_0, g_xxxzzz_0_xyzzz_0, g_xxxzzz_0_xzzzz_0, g_xxxzzz_0_yyyyy_0, g_xxxzzz_0_yyyyz_0, g_xxxzzz_0_yyyzz_0, g_xxxzzz_0_yyzzz_0, g_xxxzzz_0_yzzzz_0, g_xxxzzz_0_zzzzz_0, g_xxzzz_0_xxxxz_1, g_xxzzz_0_xxxyz_1, g_xxzzz_0_xxxz_1, g_xxzzz_0_xxxzz_1, g_xxzzz_0_xxyyz_1, g_xxzzz_0_xxyz_1, g_xxzzz_0_xxyzz_1, g_xxzzz_0_xxzz_1, g_xxzzz_0_xxzzz_1, g_xxzzz_0_xyyyz_1, g_xxzzz_0_xyyz_1, g_xxzzz_0_xyyzz_1, g_xxzzz_0_xyzz_1, g_xxzzz_0_xyzzz_1, g_xxzzz_0_xzzz_1, g_xxzzz_0_xzzzz_1, g_xxzzz_0_yyyyy_1, g_xxzzz_0_yyyyz_1, g_xxzzz_0_yyyz_1, g_xxzzz_0_yyyzz_1, g_xxzzz_0_yyzz_1, g_xxzzz_0_yyzzz_1, g_xxzzz_0_yzzz_1, g_xxzzz_0_yzzzz_1, g_xxzzz_0_zzzz_1, g_xxzzz_0_zzzzz_1, g_xzzz_0_xxxxz_0, g_xzzz_0_xxxxz_1, g_xzzz_0_xxxyz_0, g_xzzz_0_xxxyz_1, g_xzzz_0_xxxzz_0, g_xzzz_0_xxxzz_1, g_xzzz_0_xxyyz_0, g_xzzz_0_xxyyz_1, g_xzzz_0_xxyzz_0, g_xzzz_0_xxyzz_1, g_xzzz_0_xxzzz_0, g_xzzz_0_xxzzz_1, g_xzzz_0_xyyyz_0, g_xzzz_0_xyyyz_1, g_xzzz_0_xyyzz_0, g_xzzz_0_xyyzz_1, g_xzzz_0_xyzzz_0, g_xzzz_0_xyzzz_1, g_xzzz_0_xzzzz_0, g_xzzz_0_xzzzz_1, g_xzzz_0_yyyyy_0, g_xzzz_0_yyyyy_1, g_xzzz_0_yyyyz_0, g_xzzz_0_yyyyz_1, g_xzzz_0_yyyzz_0, g_xzzz_0_yyyzz_1, g_xzzz_0_yyzzz_0, g_xzzz_0_yyzzz_1, g_xzzz_0_yzzzz_0, g_xzzz_0_yzzzz_1, g_xzzz_0_zzzzz_0, g_xzzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xxxxx_0[i] = 2.0 * g_xxxz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxx_1[i] * fz_be_0 + g_xxxzz_0_xxxxx_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxy_0[i] = 2.0 * g_xxxz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxy_1[i] * fz_be_0 + g_xxxzz_0_xxxxy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxz_0[i] = 2.0 * g_xzzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyy_0[i] = 2.0 * g_xxxz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxyy_1[i] * fz_be_0 + g_xxxzz_0_xxxyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxyz_0[i] = 2.0 * g_xzzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxzz_0[i] = 2.0 * g_xzzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyy_0[i] = 2.0 * g_xxxz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxyyy_1[i] * fz_be_0 + g_xxxzz_0_xxyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxyyz_0[i] = 2.0 * g_xzzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyzz_0[i] = 2.0 * g_xzzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxzzz_0[i] = 2.0 * g_xzzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyy_0[i] = 2.0 * g_xxxz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xyyyy_1[i] * fz_be_0 + g_xxxzz_0_xyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xyyyz_0[i] = 2.0 * g_xzzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyzz_0[i] = 2.0 * g_xzzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyzzz_0[i] = 2.0 * g_xzzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xzzzz_0[i] = 2.0 * g_xzzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzz_1[i] * fi_acd_0 + g_xxzzz_0_xzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyy_0[i] = 2.0 * g_xzzz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyy_1[i] * fz_be_0 + g_xxzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyz_0[i] = 2.0 * g_xzzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyzz_0[i] = 2.0 * g_xzzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyzzz_0[i] = 2.0 * g_xzzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yzzzz_0[i] = 2.0 * g_xzzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_zzzzz_0[i] = 2.0 * g_xzzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 210-231 components of targeted buffer : ISH

    auto g_xxyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 210);

    auto g_xxyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 211);

    auto g_xxyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 212);

    auto g_xxyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 213);

    auto g_xxyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 214);

    auto g_xxyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 215);

    auto g_xxyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 216);

    auto g_xxyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 217);

    auto g_xxyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 218);

    auto g_xxyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 219);

    auto g_xxyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 220);

    auto g_xxyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 221);

    auto g_xxyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 222);

    auto g_xxyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 223);

    auto g_xxyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 224);

    auto g_xxyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 225);

    auto g_xxyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 226);

    auto g_xxyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 227);

    auto g_xxyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 228);

    auto g_xxyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 229);

    auto g_xxyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 230);

    #pragma omp simd aligned(g_xxyy_0_xxxxx_0, g_xxyy_0_xxxxx_1, g_xxyy_0_xxxxz_0, g_xxyy_0_xxxxz_1, g_xxyy_0_xxxzz_0, g_xxyy_0_xxxzz_1, g_xxyy_0_xxzzz_0, g_xxyy_0_xxzzz_1, g_xxyy_0_xzzzz_0, g_xxyy_0_xzzzz_1, g_xxyyy_0_xxxxx_1, g_xxyyy_0_xxxxz_1, g_xxyyy_0_xxxzz_1, g_xxyyy_0_xxzzz_1, g_xxyyy_0_xzzzz_1, g_xxyyyy_0_xxxxx_0, g_xxyyyy_0_xxxxy_0, g_xxyyyy_0_xxxxz_0, g_xxyyyy_0_xxxyy_0, g_xxyyyy_0_xxxyz_0, g_xxyyyy_0_xxxzz_0, g_xxyyyy_0_xxyyy_0, g_xxyyyy_0_xxyyz_0, g_xxyyyy_0_xxyzz_0, g_xxyyyy_0_xxzzz_0, g_xxyyyy_0_xyyyy_0, g_xxyyyy_0_xyyyz_0, g_xxyyyy_0_xyyzz_0, g_xxyyyy_0_xyzzz_0, g_xxyyyy_0_xzzzz_0, g_xxyyyy_0_yyyyy_0, g_xxyyyy_0_yyyyz_0, g_xxyyyy_0_yyyzz_0, g_xxyyyy_0_yyzzz_0, g_xxyyyy_0_yzzzz_0, g_xxyyyy_0_zzzzz_0, g_xyyyy_0_xxxxy_1, g_xyyyy_0_xxxy_1, g_xyyyy_0_xxxyy_1, g_xyyyy_0_xxxyz_1, g_xyyyy_0_xxyy_1, g_xyyyy_0_xxyyy_1, g_xyyyy_0_xxyyz_1, g_xyyyy_0_xxyz_1, g_xyyyy_0_xxyzz_1, g_xyyyy_0_xyyy_1, g_xyyyy_0_xyyyy_1, g_xyyyy_0_xyyyz_1, g_xyyyy_0_xyyz_1, g_xyyyy_0_xyyzz_1, g_xyyyy_0_xyzz_1, g_xyyyy_0_xyzzz_1, g_xyyyy_0_yyyy_1, g_xyyyy_0_yyyyy_1, g_xyyyy_0_yyyyz_1, g_xyyyy_0_yyyz_1, g_xyyyy_0_yyyzz_1, g_xyyyy_0_yyzz_1, g_xyyyy_0_yyzzz_1, g_xyyyy_0_yzzz_1, g_xyyyy_0_yzzzz_1, g_xyyyy_0_zzzzz_1, g_yyyy_0_xxxxy_0, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxyy_0, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyz_0, g_yyyy_0_xxxyz_1, g_yyyy_0_xxyyy_0, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyz_0, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyzz_0, g_yyyy_0_xxyzz_1, g_yyyy_0_xyyyy_0, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyz_0, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyzz_0, g_yyyy_0_xyyzz_1, g_yyyy_0_xyzzz_0, g_yyyy_0_xyzzz_1, g_yyyy_0_yyyyy_0, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyz_0, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyzz_0, g_yyyy_0_yyyzz_1, g_yyyy_0_yyzzz_0, g_yyyy_0_yyzzz_1, g_yyyy_0_yzzzz_0, g_yyyy_0_yzzzz_1, g_yyyy_0_zzzzz_0, g_yyyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xxxxx_0[i] = 3.0 * g_xxyy_0_xxxxx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxx_1[i] * fz_be_0 + g_xxyyy_0_xxxxx_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxy_0[i] = g_yyyy_0_xxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxz_0[i] = 3.0 * g_xxyy_0_xxxxz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxz_1[i] * fz_be_0 + g_xxyyy_0_xxxxz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxyy_0[i] = g_yyyy_0_xxxyy_0[i] * fbe_0 - g_yyyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyz_0[i] = g_yyyy_0_xxxyz_0[i] * fbe_0 - g_yyyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxzz_0[i] = 3.0 * g_xxyy_0_xxxzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxzz_1[i] * fz_be_0 + g_xxyyy_0_xxxzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxyyy_0[i] = g_yyyy_0_xxyyy_0[i] * fbe_0 - g_yyyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyz_0[i] = g_yyyy_0_xxyyz_0[i] * fbe_0 - g_yyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyzz_0[i] = g_yyyy_0_xxyzz_0[i] * fbe_0 - g_yyyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxzzz_0[i] = 3.0 * g_xxyy_0_xxzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxzzz_1[i] * fz_be_0 + g_xxyyy_0_xxzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xyyyy_0[i] = g_yyyy_0_xyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyy_1[i] * fi_acd_0 + g_xyyyy_0_xyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyz_0[i] = g_yyyy_0_xyyyz_0[i] * fbe_0 - g_yyyy_0_xyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyzz_0[i] = g_yyyy_0_xyyzz_0[i] * fbe_0 - g_yyyy_0_xyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyzzz_0[i] = g_yyyy_0_xyzzz_0[i] * fbe_0 - g_yyyy_0_xyzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xzzzz_0[i] = 3.0 * g_xxyy_0_xzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xzzzz_1[i] * fz_be_0 + g_xxyyy_0_xzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_yyyyy_0[i] = g_yyyy_0_yyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyz_0[i] = g_yyyy_0_yyyyz_0[i] * fbe_0 - g_yyyy_0_yyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyzz_0[i] = g_yyyy_0_yyyzz_0[i] * fbe_0 - g_yyyy_0_yyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyzzz_0[i] = g_yyyy_0_yyzzz_0[i] * fbe_0 - g_yyyy_0_yyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yzzzz_0[i] = g_yyyy_0_yzzzz_0[i] * fbe_0 - g_yyyy_0_yzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_zzzzz_0[i] = g_yyyy_0_zzzzz_0[i] * fbe_0 - g_yyyy_0_zzzzz_1[i] * fz_be_0 + g_xyyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 231-252 components of targeted buffer : ISH

    auto g_xxyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 231);

    auto g_xxyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 232);

    auto g_xxyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 233);

    auto g_xxyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 234);

    auto g_xxyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 235);

    auto g_xxyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 236);

    auto g_xxyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 237);

    auto g_xxyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 238);

    auto g_xxyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 239);

    auto g_xxyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 240);

    auto g_xxyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 241);

    auto g_xxyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 242);

    auto g_xxyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 243);

    auto g_xxyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 244);

    auto g_xxyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 245);

    auto g_xxyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 246);

    auto g_xxyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 247);

    auto g_xxyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 248);

    auto g_xxyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 249);

    auto g_xxyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 250);

    auto g_xxyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 251);

    #pragma omp simd aligned(g_xxyyy_0_xxxx_1, g_xxyyy_0_xxxxx_1, g_xxyyy_0_xxxxy_1, g_xxyyy_0_xxxxz_1, g_xxyyy_0_xxxy_1, g_xxyyy_0_xxxyy_1, g_xxyyy_0_xxxyz_1, g_xxyyy_0_xxxz_1, g_xxyyy_0_xxxzz_1, g_xxyyy_0_xxyy_1, g_xxyyy_0_xxyyy_1, g_xxyyy_0_xxyyz_1, g_xxyyy_0_xxyz_1, g_xxyyy_0_xxyzz_1, g_xxyyy_0_xxzz_1, g_xxyyy_0_xxzzz_1, g_xxyyy_0_xyyy_1, g_xxyyy_0_xyyyy_1, g_xxyyy_0_xyyyz_1, g_xxyyy_0_xyyz_1, g_xxyyy_0_xyyzz_1, g_xxyyy_0_xyzz_1, g_xxyyy_0_xyzzz_1, g_xxyyy_0_xzzz_1, g_xxyyy_0_xzzzz_1, g_xxyyy_0_yyyy_1, g_xxyyy_0_yyyyy_1, g_xxyyy_0_yyyyz_1, g_xxyyy_0_yyyz_1, g_xxyyy_0_yyyzz_1, g_xxyyy_0_yyzz_1, g_xxyyy_0_yyzzz_1, g_xxyyy_0_yzzz_1, g_xxyyy_0_yzzzz_1, g_xxyyy_0_zzzz_1, g_xxyyy_0_zzzzz_1, g_xxyyyz_0_xxxxx_0, g_xxyyyz_0_xxxxy_0, g_xxyyyz_0_xxxxz_0, g_xxyyyz_0_xxxyy_0, g_xxyyyz_0_xxxyz_0, g_xxyyyz_0_xxxzz_0, g_xxyyyz_0_xxyyy_0, g_xxyyyz_0_xxyyz_0, g_xxyyyz_0_xxyzz_0, g_xxyyyz_0_xxzzz_0, g_xxyyyz_0_xyyyy_0, g_xxyyyz_0_xyyyz_0, g_xxyyyz_0_xyyzz_0, g_xxyyyz_0_xyzzz_0, g_xxyyyz_0_xzzzz_0, g_xxyyyz_0_yyyyy_0, g_xxyyyz_0_yyyyz_0, g_xxyyyz_0_yyyzz_0, g_xxyyyz_0_yyzzz_0, g_xxyyyz_0_yzzzz_0, g_xxyyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xxxxx_0[i] = g_xxyyy_0_xxxxx_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxy_0[i] = g_xxyyy_0_xxxxy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxz_0[i] = g_xxyyy_0_xxxx_1[i] * fi_acd_0 + g_xxyyy_0_xxxxz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyy_0[i] = g_xxyyy_0_xxxyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyz_0[i] = g_xxyyy_0_xxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxzz_0[i] = 2.0 * g_xxyyy_0_xxxz_1[i] * fi_acd_0 + g_xxyyy_0_xxxzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyy_0[i] = g_xxyyy_0_xxyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyz_0[i] = g_xxyyy_0_xxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyzz_0[i] = 2.0 * g_xxyyy_0_xxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxzzz_0[i] = 3.0 * g_xxyyy_0_xxzz_1[i] * fi_acd_0 + g_xxyyy_0_xxzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyy_0[i] = g_xxyyy_0_xyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyz_0[i] = g_xxyyy_0_xyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyzz_0[i] = 2.0 * g_xxyyy_0_xyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyzzz_0[i] = 3.0 * g_xxyyy_0_xyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xzzzz_0[i] = 4.0 * g_xxyyy_0_xzzz_1[i] * fi_acd_0 + g_xxyyy_0_xzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyy_0[i] = g_xxyyy_0_yyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyz_0[i] = g_xxyyy_0_yyyy_1[i] * fi_acd_0 + g_xxyyy_0_yyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyzz_0[i] = 2.0 * g_xxyyy_0_yyyz_1[i] * fi_acd_0 + g_xxyyy_0_yyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyzzz_0[i] = 3.0 * g_xxyyy_0_yyzz_1[i] * fi_acd_0 + g_xxyyy_0_yyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yzzzz_0[i] = 4.0 * g_xxyyy_0_yzzz_1[i] * fi_acd_0 + g_xxyyy_0_yzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_zzzzz_0[i] = 5.0 * g_xxyyy_0_zzzz_1[i] * fi_acd_0 + g_xxyyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 252-273 components of targeted buffer : ISH

    auto g_xxyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 252);

    auto g_xxyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 253);

    auto g_xxyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 254);

    auto g_xxyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 255);

    auto g_xxyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 256);

    auto g_xxyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 257);

    auto g_xxyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 258);

    auto g_xxyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 259);

    auto g_xxyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 260);

    auto g_xxyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 261);

    auto g_xxyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 262);

    auto g_xxyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 263);

    auto g_xxyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 264);

    auto g_xxyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 265);

    auto g_xxyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 266);

    auto g_xxyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 267);

    auto g_xxyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 268);

    auto g_xxyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 269);

    auto g_xxyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 270);

    auto g_xxyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 271);

    auto g_xxyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 272);

    #pragma omp simd aligned(g_xxyy_0_xxxxy_0, g_xxyy_0_xxxxy_1, g_xxyy_0_xxxyy_0, g_xxyy_0_xxxyy_1, g_xxyy_0_xxyyy_0, g_xxyy_0_xxyyy_1, g_xxyy_0_xyyyy_0, g_xxyy_0_xyyyy_1, g_xxyyz_0_xxxxy_1, g_xxyyz_0_xxxyy_1, g_xxyyz_0_xxyyy_1, g_xxyyz_0_xyyyy_1, g_xxyyzz_0_xxxxx_0, g_xxyyzz_0_xxxxy_0, g_xxyyzz_0_xxxxz_0, g_xxyyzz_0_xxxyy_0, g_xxyyzz_0_xxxyz_0, g_xxyyzz_0_xxxzz_0, g_xxyyzz_0_xxyyy_0, g_xxyyzz_0_xxyyz_0, g_xxyyzz_0_xxyzz_0, g_xxyyzz_0_xxzzz_0, g_xxyyzz_0_xyyyy_0, g_xxyyzz_0_xyyyz_0, g_xxyyzz_0_xyyzz_0, g_xxyyzz_0_xyzzz_0, g_xxyyzz_0_xzzzz_0, g_xxyyzz_0_yyyyy_0, g_xxyyzz_0_yyyyz_0, g_xxyyzz_0_yyyzz_0, g_xxyyzz_0_yyzzz_0, g_xxyyzz_0_yzzzz_0, g_xxyyzz_0_zzzzz_0, g_xxyzz_0_xxxxx_1, g_xxyzz_0_xxxxz_1, g_xxyzz_0_xxxzz_1, g_xxyzz_0_xxzzz_1, g_xxyzz_0_xzzzz_1, g_xxzz_0_xxxxx_0, g_xxzz_0_xxxxx_1, g_xxzz_0_xxxxz_0, g_xxzz_0_xxxxz_1, g_xxzz_0_xxxzz_0, g_xxzz_0_xxxzz_1, g_xxzz_0_xxzzz_0, g_xxzz_0_xxzzz_1, g_xxzz_0_xzzzz_0, g_xxzz_0_xzzzz_1, g_xyyzz_0_xxxyz_1, g_xyyzz_0_xxyyz_1, g_xyyzz_0_xxyz_1, g_xyyzz_0_xxyzz_1, g_xyyzz_0_xyyyz_1, g_xyyzz_0_xyyz_1, g_xyyzz_0_xyyzz_1, g_xyyzz_0_xyzz_1, g_xyyzz_0_xyzzz_1, g_xyyzz_0_yyyyy_1, g_xyyzz_0_yyyyz_1, g_xyyzz_0_yyyz_1, g_xyyzz_0_yyyzz_1, g_xyyzz_0_yyzz_1, g_xyyzz_0_yyzzz_1, g_xyyzz_0_yzzz_1, g_xyyzz_0_yzzzz_1, g_xyyzz_0_zzzzz_1, g_yyzz_0_xxxyz_0, g_yyzz_0_xxxyz_1, g_yyzz_0_xxyyz_0, g_yyzz_0_xxyyz_1, g_yyzz_0_xxyzz_0, g_yyzz_0_xxyzz_1, g_yyzz_0_xyyyz_0, g_yyzz_0_xyyyz_1, g_yyzz_0_xyyzz_0, g_yyzz_0_xyyzz_1, g_yyzz_0_xyzzz_0, g_yyzz_0_xyzzz_1, g_yyzz_0_yyyyy_0, g_yyzz_0_yyyyy_1, g_yyzz_0_yyyyz_0, g_yyzz_0_yyyyz_1, g_yyzz_0_yyyzz_0, g_yyzz_0_yyyzz_1, g_yyzz_0_yyzzz_0, g_yyzz_0_yyzzz_1, g_yyzz_0_yzzzz_0, g_yyzz_0_yzzzz_1, g_yyzz_0_zzzzz_0, g_yyzz_0_zzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzz_0_xxxxx_0[i] = g_xxzz_0_xxxxx_0[i] * fbe_0 - g_xxzz_0_xxxxx_1[i] * fz_be_0 + g_xxyzz_0_xxxxx_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxy_0[i] = g_xxyy_0_xxxxy_0[i] * fbe_0 - g_xxyy_0_xxxxy_1[i] * fz_be_0 + g_xxyyz_0_xxxxy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxz_0[i] = g_xxzz_0_xxxxz_0[i] * fbe_0 - g_xxzz_0_xxxxz_1[i] * fz_be_0 + g_xxyzz_0_xxxxz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxyy_0[i] = g_xxyy_0_xxxyy_0[i] * fbe_0 - g_xxyy_0_xxxyy_1[i] * fz_be_0 + g_xxyyz_0_xxxyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxyz_0[i] = g_yyzz_0_xxxyz_0[i] * fbe_0 - g_yyzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxzz_0[i] = g_xxzz_0_xxxzz_0[i] * fbe_0 - g_xxzz_0_xxxzz_1[i] * fz_be_0 + g_xxyzz_0_xxxzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxyyy_0[i] = g_xxyy_0_xxyyy_0[i] * fbe_0 - g_xxyy_0_xxyyy_1[i] * fz_be_0 + g_xxyyz_0_xxyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxyyz_0[i] = g_yyzz_0_xxyyz_0[i] * fbe_0 - g_yyzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyzz_0[i] = g_yyzz_0_xxyzz_0[i] * fbe_0 - g_yyzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxzzz_0[i] = g_xxzz_0_xxzzz_0[i] * fbe_0 - g_xxzz_0_xxzzz_1[i] * fz_be_0 + g_xxyzz_0_xxzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xyyyy_0[i] = g_xxyy_0_xyyyy_0[i] * fbe_0 - g_xxyy_0_xyyyy_1[i] * fz_be_0 + g_xxyyz_0_xyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xyyyz_0[i] = g_yyzz_0_xyyyz_0[i] * fbe_0 - g_yyzz_0_xyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyzz_0[i] = g_yyzz_0_xyyzz_0[i] * fbe_0 - g_yyzz_0_xyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyzzz_0[i] = g_yyzz_0_xyzzz_0[i] * fbe_0 - g_yyzz_0_xyzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xzzzz_0[i] = g_xxzz_0_xzzzz_0[i] * fbe_0 - g_xxzz_0_xzzzz_1[i] * fz_be_0 + g_xxyzz_0_xzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_yyyyy_0[i] = g_yyzz_0_yyyyy_0[i] * fbe_0 - g_yyzz_0_yyyyy_1[i] * fz_be_0 + g_xyyzz_0_yyyyy_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyz_0[i] = g_yyzz_0_yyyyz_0[i] * fbe_0 - g_yyzz_0_yyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyzz_0[i] = g_yyzz_0_yyyzz_0[i] * fbe_0 - g_yyzz_0_yyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyzzz_0[i] = g_yyzz_0_yyzzz_0[i] * fbe_0 - g_yyzz_0_yyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yzzzz_0[i] = g_yyzz_0_yzzzz_0[i] * fbe_0 - g_yyzz_0_yzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_zzzzz_0[i] = g_yyzz_0_zzzzz_0[i] * fbe_0 - g_yyzz_0_zzzzz_1[i] * fz_be_0 + g_xyyzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 273-294 components of targeted buffer : ISH

    auto g_xxyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 273);

    auto g_xxyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 274);

    auto g_xxyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 275);

    auto g_xxyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 276);

    auto g_xxyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 277);

    auto g_xxyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 278);

    auto g_xxyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 279);

    auto g_xxyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 280);

    auto g_xxyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 281);

    auto g_xxyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 282);

    auto g_xxyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 283);

    auto g_xxyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 284);

    auto g_xxyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 285);

    auto g_xxyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 286);

    auto g_xxyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 287);

    auto g_xxyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 288);

    auto g_xxyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 289);

    auto g_xxyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 290);

    auto g_xxyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 291);

    auto g_xxyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 292);

    auto g_xxyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 293);

    #pragma omp simd aligned(g_xxyzzz_0_xxxxx_0, g_xxyzzz_0_xxxxy_0, g_xxyzzz_0_xxxxz_0, g_xxyzzz_0_xxxyy_0, g_xxyzzz_0_xxxyz_0, g_xxyzzz_0_xxxzz_0, g_xxyzzz_0_xxyyy_0, g_xxyzzz_0_xxyyz_0, g_xxyzzz_0_xxyzz_0, g_xxyzzz_0_xxzzz_0, g_xxyzzz_0_xyyyy_0, g_xxyzzz_0_xyyyz_0, g_xxyzzz_0_xyyzz_0, g_xxyzzz_0_xyzzz_0, g_xxyzzz_0_xzzzz_0, g_xxyzzz_0_yyyyy_0, g_xxyzzz_0_yyyyz_0, g_xxyzzz_0_yyyzz_0, g_xxyzzz_0_yyzzz_0, g_xxyzzz_0_yzzzz_0, g_xxyzzz_0_zzzzz_0, g_xxzzz_0_xxxx_1, g_xxzzz_0_xxxxx_1, g_xxzzz_0_xxxxy_1, g_xxzzz_0_xxxxz_1, g_xxzzz_0_xxxy_1, g_xxzzz_0_xxxyy_1, g_xxzzz_0_xxxyz_1, g_xxzzz_0_xxxz_1, g_xxzzz_0_xxxzz_1, g_xxzzz_0_xxyy_1, g_xxzzz_0_xxyyy_1, g_xxzzz_0_xxyyz_1, g_xxzzz_0_xxyz_1, g_xxzzz_0_xxyzz_1, g_xxzzz_0_xxzz_1, g_xxzzz_0_xxzzz_1, g_xxzzz_0_xyyy_1, g_xxzzz_0_xyyyy_1, g_xxzzz_0_xyyyz_1, g_xxzzz_0_xyyz_1, g_xxzzz_0_xyyzz_1, g_xxzzz_0_xyzz_1, g_xxzzz_0_xyzzz_1, g_xxzzz_0_xzzz_1, g_xxzzz_0_xzzzz_1, g_xxzzz_0_yyyy_1, g_xxzzz_0_yyyyy_1, g_xxzzz_0_yyyyz_1, g_xxzzz_0_yyyz_1, g_xxzzz_0_yyyzz_1, g_xxzzz_0_yyzz_1, g_xxzzz_0_yyzzz_1, g_xxzzz_0_yzzz_1, g_xxzzz_0_yzzzz_1, g_xxzzz_0_zzzz_1, g_xxzzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xxxxx_0[i] = g_xxzzz_0_xxxxx_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxy_0[i] = g_xxzzz_0_xxxx_1[i] * fi_acd_0 + g_xxzzz_0_xxxxy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxz_0[i] = g_xxzzz_0_xxxxz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyy_0[i] = 2.0 * g_xxzzz_0_xxxy_1[i] * fi_acd_0 + g_xxzzz_0_xxxyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyz_0[i] = g_xxzzz_0_xxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxzz_0[i] = g_xxzzz_0_xxxzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyy_0[i] = 3.0 * g_xxzzz_0_xxyy_1[i] * fi_acd_0 + g_xxzzz_0_xxyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyz_0[i] = 2.0 * g_xxzzz_0_xxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyzz_0[i] = g_xxzzz_0_xxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxzzz_0[i] = g_xxzzz_0_xxzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyy_0[i] = 4.0 * g_xxzzz_0_xyyy_1[i] * fi_acd_0 + g_xxzzz_0_xyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyz_0[i] = 3.0 * g_xxzzz_0_xyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyzz_0[i] = 2.0 * g_xxzzz_0_xyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyzzz_0[i] = g_xxzzz_0_xzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xzzzz_0[i] = g_xxzzz_0_xzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyy_0[i] = 5.0 * g_xxzzz_0_yyyy_1[i] * fi_acd_0 + g_xxzzz_0_yyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyz_0[i] = 4.0 * g_xxzzz_0_yyyz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyzz_0[i] = 3.0 * g_xxzzz_0_yyzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyzzz_0[i] = 2.0 * g_xxzzz_0_yzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yzzzz_0[i] = g_xxzzz_0_zzzz_1[i] * fi_acd_0 + g_xxzzz_0_yzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_zzzzz_0[i] = g_xxzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 294-315 components of targeted buffer : ISH

    auto g_xxzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 294);

    auto g_xxzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 295);

    auto g_xxzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 296);

    auto g_xxzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 297);

    auto g_xxzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 298);

    auto g_xxzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 299);

    auto g_xxzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 300);

    auto g_xxzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 301);

    auto g_xxzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 302);

    auto g_xxzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 303);

    auto g_xxzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 304);

    auto g_xxzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 305);

    auto g_xxzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 306);

    auto g_xxzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 307);

    auto g_xxzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 308);

    auto g_xxzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 309);

    auto g_xxzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 310);

    auto g_xxzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 311);

    auto g_xxzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 312);

    auto g_xxzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 313);

    auto g_xxzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 314);

    #pragma omp simd aligned(g_xxzz_0_xxxxx_0, g_xxzz_0_xxxxx_1, g_xxzz_0_xxxxy_0, g_xxzz_0_xxxxy_1, g_xxzz_0_xxxyy_0, g_xxzz_0_xxxyy_1, g_xxzz_0_xxyyy_0, g_xxzz_0_xxyyy_1, g_xxzz_0_xyyyy_0, g_xxzz_0_xyyyy_1, g_xxzzz_0_xxxxx_1, g_xxzzz_0_xxxxy_1, g_xxzzz_0_xxxyy_1, g_xxzzz_0_xxyyy_1, g_xxzzz_0_xyyyy_1, g_xxzzzz_0_xxxxx_0, g_xxzzzz_0_xxxxy_0, g_xxzzzz_0_xxxxz_0, g_xxzzzz_0_xxxyy_0, g_xxzzzz_0_xxxyz_0, g_xxzzzz_0_xxxzz_0, g_xxzzzz_0_xxyyy_0, g_xxzzzz_0_xxyyz_0, g_xxzzzz_0_xxyzz_0, g_xxzzzz_0_xxzzz_0, g_xxzzzz_0_xyyyy_0, g_xxzzzz_0_xyyyz_0, g_xxzzzz_0_xyyzz_0, g_xxzzzz_0_xyzzz_0, g_xxzzzz_0_xzzzz_0, g_xxzzzz_0_yyyyy_0, g_xxzzzz_0_yyyyz_0, g_xxzzzz_0_yyyzz_0, g_xxzzzz_0_yyzzz_0, g_xxzzzz_0_yzzzz_0, g_xxzzzz_0_zzzzz_0, g_xzzzz_0_xxxxz_1, g_xzzzz_0_xxxyz_1, g_xzzzz_0_xxxz_1, g_xzzzz_0_xxxzz_1, g_xzzzz_0_xxyyz_1, g_xzzzz_0_xxyz_1, g_xzzzz_0_xxyzz_1, g_xzzzz_0_xxzz_1, g_xzzzz_0_xxzzz_1, g_xzzzz_0_xyyyz_1, g_xzzzz_0_xyyz_1, g_xzzzz_0_xyyzz_1, g_xzzzz_0_xyzz_1, g_xzzzz_0_xyzzz_1, g_xzzzz_0_xzzz_1, g_xzzzz_0_xzzzz_1, g_xzzzz_0_yyyyy_1, g_xzzzz_0_yyyyz_1, g_xzzzz_0_yyyz_1, g_xzzzz_0_yyyzz_1, g_xzzzz_0_yyzz_1, g_xzzzz_0_yyzzz_1, g_xzzzz_0_yzzz_1, g_xzzzz_0_yzzzz_1, g_xzzzz_0_zzzz_1, g_xzzzz_0_zzzzz_1, g_zzzz_0_xxxxz_0, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxyz_0, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxzz_0, g_zzzz_0_xxxzz_1, g_zzzz_0_xxyyz_0, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyzz_0, g_zzzz_0_xxyzz_1, g_zzzz_0_xxzzz_0, g_zzzz_0_xxzzz_1, g_zzzz_0_xyyyz_0, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyzz_0, g_zzzz_0_xyyzz_1, g_zzzz_0_xyzzz_0, g_zzzz_0_xyzzz_1, g_zzzz_0_xzzzz_0, g_zzzz_0_xzzzz_1, g_zzzz_0_yyyyy_0, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyz_0, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyzz_0, g_zzzz_0_yyyzz_1, g_zzzz_0_yyzzz_0, g_zzzz_0_yyzzz_1, g_zzzz_0_yzzzz_0, g_zzzz_0_yzzzz_1, g_zzzz_0_zzzzz_0, g_zzzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xxxxx_0[i] = 3.0 * g_xxzz_0_xxxxx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxx_1[i] * fz_be_0 + g_xxzzz_0_xxxxx_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxy_0[i] = 3.0 * g_xxzz_0_xxxxy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxy_1[i] * fz_be_0 + g_xxzzz_0_xxxxy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxz_0[i] = g_zzzz_0_xxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyy_0[i] = 3.0 * g_xxzz_0_xxxyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyy_1[i] * fz_be_0 + g_xxzzz_0_xxxyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxyz_0[i] = g_zzzz_0_xxxyz_0[i] * fbe_0 - g_zzzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxzz_0[i] = g_zzzz_0_xxxzz_0[i] * fbe_0 - g_zzzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyy_0[i] = 3.0 * g_xxzz_0_xxyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyy_1[i] * fz_be_0 + g_xxzzz_0_xxyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxyyz_0[i] = g_zzzz_0_xxyyz_0[i] * fbe_0 - g_zzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyzz_0[i] = g_zzzz_0_xxyzz_0[i] * fbe_0 - g_zzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxzzz_0[i] = g_zzzz_0_xxzzz_0[i] * fbe_0 - g_zzzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyy_0[i] = 3.0 * g_xxzz_0_xyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyy_1[i] * fz_be_0 + g_xxzzz_0_xyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xyyyz_0[i] = g_zzzz_0_xyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyzz_0[i] = g_zzzz_0_xyyzz_0[i] * fbe_0 - g_zzzz_0_xyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyzzz_0[i] = g_zzzz_0_xyzzz_0[i] * fbe_0 - g_zzzz_0_xyzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xzzzz_0[i] = g_zzzz_0_xzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzz_1[i] * fi_acd_0 + g_xzzzz_0_xzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyy_0[i] = g_zzzz_0_yyyyy_0[i] * fbe_0 - g_zzzz_0_yyyyy_1[i] * fz_be_0 + g_xzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyz_0[i] = g_zzzz_0_yyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyzz_0[i] = g_zzzz_0_yyyzz_0[i] * fbe_0 - g_zzzz_0_yyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyzzz_0[i] = g_zzzz_0_yyzzz_0[i] * fbe_0 - g_zzzz_0_yyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yzzzz_0[i] = g_zzzz_0_yzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_zzzzz_0[i] = g_zzzz_0_zzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-336 components of targeted buffer : ISH

    auto g_xyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 315);

    auto g_xyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 316);

    auto g_xyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 317);

    auto g_xyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 318);

    auto g_xyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 319);

    auto g_xyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 320);

    auto g_xyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 321);

    auto g_xyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 322);

    auto g_xyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 323);

    auto g_xyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 324);

    auto g_xyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 325);

    auto g_xyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 326);

    auto g_xyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 327);

    auto g_xyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 328);

    auto g_xyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 329);

    auto g_xyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 330);

    auto g_xyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 331);

    auto g_xyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 332);

    auto g_xyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 333);

    auto g_xyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 334);

    auto g_xyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 335);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxx_0, g_xyyyyy_0_xxxxy_0, g_xyyyyy_0_xxxxz_0, g_xyyyyy_0_xxxyy_0, g_xyyyyy_0_xxxyz_0, g_xyyyyy_0_xxxzz_0, g_xyyyyy_0_xxyyy_0, g_xyyyyy_0_xxyyz_0, g_xyyyyy_0_xxyzz_0, g_xyyyyy_0_xxzzz_0, g_xyyyyy_0_xyyyy_0, g_xyyyyy_0_xyyyz_0, g_xyyyyy_0_xyyzz_0, g_xyyyyy_0_xyzzz_0, g_xyyyyy_0_xzzzz_0, g_xyyyyy_0_yyyyy_0, g_xyyyyy_0_yyyyz_0, g_xyyyyy_0_yyyzz_0, g_xyyyyy_0_yyzzz_0, g_xyyyyy_0_yzzzz_0, g_xyyyyy_0_zzzzz_0, g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxzz_1, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyzz_1, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xzzz_1, g_yyyyy_0_xzzzz_1, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyzz_1, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yzzz_1, g_yyyyy_0_yzzzz_1, g_yyyyy_0_zzzz_1, g_yyyyy_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xxxxx_0[i] = 5.0 * g_yyyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxx_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxy_0[i] = 4.0 * g_yyyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxz_0[i] = 4.0 * g_yyyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyy_0[i] = 3.0 * g_yyyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyz_0[i] = 3.0 * g_yyyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxzz_0[i] = 3.0 * g_yyyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyy_0[i] = 2.0 * g_yyyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyz_0[i] = 2.0 * g_yyyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyzz_0[i] = 2.0 * g_yyyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxzzz_0[i] = 2.0 * g_yyyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyy_0[i] = g_yyyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyz_0[i] = g_yyyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyzz_0[i] = g_yyyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyzzz_0[i] = g_yyyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xzzzz_0[i] = g_yyyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyy_0[i] = g_yyyyy_0_yyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyz_0[i] = g_yyyyy_0_yyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyzz_0[i] = g_yyyyy_0_yyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyzzz_0[i] = g_yyyyy_0_yyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yzzzz_0[i] = g_yyyyy_0_yzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_zzzzz_0[i] = g_yyyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 336-357 components of targeted buffer : ISH

    auto g_xyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 336);

    auto g_xyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 337);

    auto g_xyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 338);

    auto g_xyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 339);

    auto g_xyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 340);

    auto g_xyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 341);

    auto g_xyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 342);

    auto g_xyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 343);

    auto g_xyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 344);

    auto g_xyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 345);

    auto g_xyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 346);

    auto g_xyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 347);

    auto g_xyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 348);

    auto g_xyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 349);

    auto g_xyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 350);

    auto g_xyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 351);

    auto g_xyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 352);

    auto g_xyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 353);

    auto g_xyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 354);

    auto g_xyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 355);

    auto g_xyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 356);

    #pragma omp simd aligned(g_xyyyy_0_xxxxx_1, g_xyyyy_0_xxxxy_1, g_xyyyy_0_xxxyy_1, g_xyyyy_0_xxyyy_1, g_xyyyy_0_xyyyy_1, g_xyyyyz_0_xxxxx_0, g_xyyyyz_0_xxxxy_0, g_xyyyyz_0_xxxxz_0, g_xyyyyz_0_xxxyy_0, g_xyyyyz_0_xxxyz_0, g_xyyyyz_0_xxxzz_0, g_xyyyyz_0_xxyyy_0, g_xyyyyz_0_xxyyz_0, g_xyyyyz_0_xxyzz_0, g_xyyyyz_0_xxzzz_0, g_xyyyyz_0_xyyyy_0, g_xyyyyz_0_xyyyz_0, g_xyyyyz_0_xyyzz_0, g_xyyyyz_0_xyzzz_0, g_xyyyyz_0_xzzzz_0, g_xyyyyz_0_yyyyy_0, g_xyyyyz_0_yyyyz_0, g_xyyyyz_0_yyyzz_0, g_xyyyyz_0_yyzzz_0, g_xyyyyz_0_yzzzz_0, g_xyyyyz_0_zzzzz_0, g_yyyyz_0_xxxxz_1, g_yyyyz_0_xxxyz_1, g_yyyyz_0_xxxz_1, g_yyyyz_0_xxxzz_1, g_yyyyz_0_xxyyz_1, g_yyyyz_0_xxyz_1, g_yyyyz_0_xxyzz_1, g_yyyyz_0_xxzz_1, g_yyyyz_0_xxzzz_1, g_yyyyz_0_xyyyz_1, g_yyyyz_0_xyyz_1, g_yyyyz_0_xyyzz_1, g_yyyyz_0_xyzz_1, g_yyyyz_0_xyzzz_1, g_yyyyz_0_xzzz_1, g_yyyyz_0_xzzzz_1, g_yyyyz_0_yyyyy_1, g_yyyyz_0_yyyyz_1, g_yyyyz_0_yyyz_1, g_yyyyz_0_yyyzz_1, g_yyyyz_0_yyzz_1, g_yyyyz_0_yyzzz_1, g_yyyyz_0_yzzz_1, g_yyyyz_0_yzzzz_1, g_yyyyz_0_zzzz_1, g_yyyyz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xxxxx_0[i] = g_xyyyy_0_xxxxx_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxy_0[i] = g_xyyyy_0_xxxxy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxz_0[i] = 4.0 * g_yyyyz_0_xxxz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyy_0[i] = g_xyyyy_0_xxxyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxyz_0[i] = 3.0 * g_yyyyz_0_xxyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxzz_0[i] = 3.0 * g_yyyyz_0_xxzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyy_0[i] = g_xyyyy_0_xxyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxyyz_0[i] = 2.0 * g_yyyyz_0_xyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyzz_0[i] = 2.0 * g_yyyyz_0_xyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxzzz_0[i] = 2.0 * g_yyyyz_0_xzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyy_0[i] = g_xyyyy_0_xyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xyyyz_0[i] = g_yyyyz_0_yyyz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyzz_0[i] = g_yyyyz_0_yyzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyzzz_0[i] = g_yyyyz_0_yzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xzzzz_0[i] = g_yyyyz_0_zzzz_1[i] * fi_acd_0 + g_yyyyz_0_xzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyy_0[i] = g_yyyyz_0_yyyyy_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyz_0[i] = g_yyyyz_0_yyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyzz_0[i] = g_yyyyz_0_yyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyzzz_0[i] = g_yyyyz_0_yyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yzzzz_0[i] = g_yyyyz_0_yzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_zzzzz_0[i] = g_yyyyz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 357-378 components of targeted buffer : ISH

    auto g_xyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 357);

    auto g_xyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 358);

    auto g_xyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 359);

    auto g_xyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 360);

    auto g_xyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 361);

    auto g_xyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 362);

    auto g_xyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 363);

    auto g_xyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 364);

    auto g_xyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 365);

    auto g_xyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 366);

    auto g_xyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 367);

    auto g_xyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 368);

    auto g_xyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 369);

    auto g_xyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 370);

    auto g_xyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 371);

    auto g_xyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 372);

    auto g_xyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 373);

    auto g_xyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 374);

    auto g_xyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 375);

    auto g_xyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 376);

    auto g_xyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 377);

    #pragma omp simd aligned(g_xyyyzz_0_xxxxx_0, g_xyyyzz_0_xxxxy_0, g_xyyyzz_0_xxxxz_0, g_xyyyzz_0_xxxyy_0, g_xyyyzz_0_xxxyz_0, g_xyyyzz_0_xxxzz_0, g_xyyyzz_0_xxyyy_0, g_xyyyzz_0_xxyyz_0, g_xyyyzz_0_xxyzz_0, g_xyyyzz_0_xxzzz_0, g_xyyyzz_0_xyyyy_0, g_xyyyzz_0_xyyyz_0, g_xyyyzz_0_xyyzz_0, g_xyyyzz_0_xyzzz_0, g_xyyyzz_0_xzzzz_0, g_xyyyzz_0_yyyyy_0, g_xyyyzz_0_yyyyz_0, g_xyyyzz_0_yyyzz_0, g_xyyyzz_0_yyzzz_0, g_xyyyzz_0_yzzzz_0, g_xyyyzz_0_zzzzz_0, g_yyyzz_0_xxxx_1, g_yyyzz_0_xxxxx_1, g_yyyzz_0_xxxxy_1, g_yyyzz_0_xxxxz_1, g_yyyzz_0_xxxy_1, g_yyyzz_0_xxxyy_1, g_yyyzz_0_xxxyz_1, g_yyyzz_0_xxxz_1, g_yyyzz_0_xxxzz_1, g_yyyzz_0_xxyy_1, g_yyyzz_0_xxyyy_1, g_yyyzz_0_xxyyz_1, g_yyyzz_0_xxyz_1, g_yyyzz_0_xxyzz_1, g_yyyzz_0_xxzz_1, g_yyyzz_0_xxzzz_1, g_yyyzz_0_xyyy_1, g_yyyzz_0_xyyyy_1, g_yyyzz_0_xyyyz_1, g_yyyzz_0_xyyz_1, g_yyyzz_0_xyyzz_1, g_yyyzz_0_xyzz_1, g_yyyzz_0_xyzzz_1, g_yyyzz_0_xzzz_1, g_yyyzz_0_xzzzz_1, g_yyyzz_0_yyyy_1, g_yyyzz_0_yyyyy_1, g_yyyzz_0_yyyyz_1, g_yyyzz_0_yyyz_1, g_yyyzz_0_yyyzz_1, g_yyyzz_0_yyzz_1, g_yyyzz_0_yyzzz_1, g_yyyzz_0_yzzz_1, g_yyyzz_0_yzzzz_1, g_yyyzz_0_zzzz_1, g_yyyzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xxxxx_0[i] = 5.0 * g_yyyzz_0_xxxx_1[i] * fi_acd_0 + g_yyyzz_0_xxxxx_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxy_0[i] = 4.0 * g_yyyzz_0_xxxy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxz_0[i] = 4.0 * g_yyyzz_0_xxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyy_0[i] = 3.0 * g_yyyzz_0_xxyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyz_0[i] = 3.0 * g_yyyzz_0_xxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxzz_0[i] = 3.0 * g_yyyzz_0_xxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyy_0[i] = 2.0 * g_yyyzz_0_xyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyz_0[i] = 2.0 * g_yyyzz_0_xyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyzz_0[i] = 2.0 * g_yyyzz_0_xyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxzzz_0[i] = 2.0 * g_yyyzz_0_xzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyy_0[i] = g_yyyzz_0_yyyy_1[i] * fi_acd_0 + g_yyyzz_0_xyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyz_0[i] = g_yyyzz_0_yyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyzz_0[i] = g_yyyzz_0_yyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyzzz_0[i] = g_yyyzz_0_yzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xzzzz_0[i] = g_yyyzz_0_zzzz_1[i] * fi_acd_0 + g_yyyzz_0_xzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyy_0[i] = g_yyyzz_0_yyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyz_0[i] = g_yyyzz_0_yyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyzz_0[i] = g_yyyzz_0_yyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyzzz_0[i] = g_yyyzz_0_yyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yzzzz_0[i] = g_yyyzz_0_yzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_zzzzz_0[i] = g_yyyzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 378-399 components of targeted buffer : ISH

    auto g_xyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 378);

    auto g_xyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 379);

    auto g_xyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 380);

    auto g_xyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 381);

    auto g_xyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 382);

    auto g_xyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 383);

    auto g_xyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 384);

    auto g_xyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 385);

    auto g_xyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 386);

    auto g_xyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 387);

    auto g_xyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 388);

    auto g_xyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 389);

    auto g_xyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 390);

    auto g_xyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 391);

    auto g_xyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 392);

    auto g_xyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 393);

    auto g_xyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 394);

    auto g_xyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 395);

    auto g_xyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 396);

    auto g_xyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 397);

    auto g_xyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 398);

    #pragma omp simd aligned(g_xyyzzz_0_xxxxx_0, g_xyyzzz_0_xxxxy_0, g_xyyzzz_0_xxxxz_0, g_xyyzzz_0_xxxyy_0, g_xyyzzz_0_xxxyz_0, g_xyyzzz_0_xxxzz_0, g_xyyzzz_0_xxyyy_0, g_xyyzzz_0_xxyyz_0, g_xyyzzz_0_xxyzz_0, g_xyyzzz_0_xxzzz_0, g_xyyzzz_0_xyyyy_0, g_xyyzzz_0_xyyyz_0, g_xyyzzz_0_xyyzz_0, g_xyyzzz_0_xyzzz_0, g_xyyzzz_0_xzzzz_0, g_xyyzzz_0_yyyyy_0, g_xyyzzz_0_yyyyz_0, g_xyyzzz_0_yyyzz_0, g_xyyzzz_0_yyzzz_0, g_xyyzzz_0_yzzzz_0, g_xyyzzz_0_zzzzz_0, g_yyzzz_0_xxxx_1, g_yyzzz_0_xxxxx_1, g_yyzzz_0_xxxxy_1, g_yyzzz_0_xxxxz_1, g_yyzzz_0_xxxy_1, g_yyzzz_0_xxxyy_1, g_yyzzz_0_xxxyz_1, g_yyzzz_0_xxxz_1, g_yyzzz_0_xxxzz_1, g_yyzzz_0_xxyy_1, g_yyzzz_0_xxyyy_1, g_yyzzz_0_xxyyz_1, g_yyzzz_0_xxyz_1, g_yyzzz_0_xxyzz_1, g_yyzzz_0_xxzz_1, g_yyzzz_0_xxzzz_1, g_yyzzz_0_xyyy_1, g_yyzzz_0_xyyyy_1, g_yyzzz_0_xyyyz_1, g_yyzzz_0_xyyz_1, g_yyzzz_0_xyyzz_1, g_yyzzz_0_xyzz_1, g_yyzzz_0_xyzzz_1, g_yyzzz_0_xzzz_1, g_yyzzz_0_xzzzz_1, g_yyzzz_0_yyyy_1, g_yyzzz_0_yyyyy_1, g_yyzzz_0_yyyyz_1, g_yyzzz_0_yyyz_1, g_yyzzz_0_yyyzz_1, g_yyzzz_0_yyzz_1, g_yyzzz_0_yyzzz_1, g_yyzzz_0_yzzz_1, g_yyzzz_0_yzzzz_1, g_yyzzz_0_zzzz_1, g_yyzzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xxxxx_0[i] = 5.0 * g_yyzzz_0_xxxx_1[i] * fi_acd_0 + g_yyzzz_0_xxxxx_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxy_0[i] = 4.0 * g_yyzzz_0_xxxy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxz_0[i] = 4.0 * g_yyzzz_0_xxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyy_0[i] = 3.0 * g_yyzzz_0_xxyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyz_0[i] = 3.0 * g_yyzzz_0_xxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxzz_0[i] = 3.0 * g_yyzzz_0_xxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyy_0[i] = 2.0 * g_yyzzz_0_xyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyz_0[i] = 2.0 * g_yyzzz_0_xyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyzz_0[i] = 2.0 * g_yyzzz_0_xyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxzzz_0[i] = 2.0 * g_yyzzz_0_xzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyy_0[i] = g_yyzzz_0_yyyy_1[i] * fi_acd_0 + g_yyzzz_0_xyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyz_0[i] = g_yyzzz_0_yyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyzz_0[i] = g_yyzzz_0_yyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyzzz_0[i] = g_yyzzz_0_yzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xzzzz_0[i] = g_yyzzz_0_zzzz_1[i] * fi_acd_0 + g_yyzzz_0_xzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyy_0[i] = g_yyzzz_0_yyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyz_0[i] = g_yyzzz_0_yyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyzz_0[i] = g_yyzzz_0_yyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyzzz_0[i] = g_yyzzz_0_yyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yzzzz_0[i] = g_yyzzz_0_yzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_zzzzz_0[i] = g_yyzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 399-420 components of targeted buffer : ISH

    auto g_xyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 399);

    auto g_xyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 400);

    auto g_xyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 401);

    auto g_xyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 402);

    auto g_xyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 403);

    auto g_xyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 404);

    auto g_xyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 405);

    auto g_xyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 406);

    auto g_xyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 407);

    auto g_xyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 408);

    auto g_xyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 409);

    auto g_xyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 410);

    auto g_xyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 411);

    auto g_xyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 412);

    auto g_xyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 413);

    auto g_xyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 414);

    auto g_xyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 415);

    auto g_xyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 416);

    auto g_xyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 417);

    auto g_xyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 418);

    auto g_xyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 419);

    #pragma omp simd aligned(g_xyzzzz_0_xxxxx_0, g_xyzzzz_0_xxxxy_0, g_xyzzzz_0_xxxxz_0, g_xyzzzz_0_xxxyy_0, g_xyzzzz_0_xxxyz_0, g_xyzzzz_0_xxxzz_0, g_xyzzzz_0_xxyyy_0, g_xyzzzz_0_xxyyz_0, g_xyzzzz_0_xxyzz_0, g_xyzzzz_0_xxzzz_0, g_xyzzzz_0_xyyyy_0, g_xyzzzz_0_xyyyz_0, g_xyzzzz_0_xyyzz_0, g_xyzzzz_0_xyzzz_0, g_xyzzzz_0_xzzzz_0, g_xyzzzz_0_yyyyy_0, g_xyzzzz_0_yyyyz_0, g_xyzzzz_0_yyyzz_0, g_xyzzzz_0_yyzzz_0, g_xyzzzz_0_yzzzz_0, g_xyzzzz_0_zzzzz_0, g_xzzzz_0_xxxxx_1, g_xzzzz_0_xxxxz_1, g_xzzzz_0_xxxzz_1, g_xzzzz_0_xxzzz_1, g_xzzzz_0_xzzzz_1, g_yzzzz_0_xxxxy_1, g_yzzzz_0_xxxy_1, g_yzzzz_0_xxxyy_1, g_yzzzz_0_xxxyz_1, g_yzzzz_0_xxyy_1, g_yzzzz_0_xxyyy_1, g_yzzzz_0_xxyyz_1, g_yzzzz_0_xxyz_1, g_yzzzz_0_xxyzz_1, g_yzzzz_0_xyyy_1, g_yzzzz_0_xyyyy_1, g_yzzzz_0_xyyyz_1, g_yzzzz_0_xyyz_1, g_yzzzz_0_xyyzz_1, g_yzzzz_0_xyzz_1, g_yzzzz_0_xyzzz_1, g_yzzzz_0_yyyy_1, g_yzzzz_0_yyyyy_1, g_yzzzz_0_yyyyz_1, g_yzzzz_0_yyyz_1, g_yzzzz_0_yyyzz_1, g_yzzzz_0_yyzz_1, g_yzzzz_0_yyzzz_1, g_yzzzz_0_yzzz_1, g_yzzzz_0_yzzzz_1, g_yzzzz_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xxxxx_0[i] = g_xzzzz_0_xxxxx_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxy_0[i] = 4.0 * g_yzzzz_0_xxxy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxz_0[i] = g_xzzzz_0_xxxxz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxyy_0[i] = 3.0 * g_yzzzz_0_xxyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyz_0[i] = 3.0 * g_yzzzz_0_xxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxzz_0[i] = g_xzzzz_0_xxxzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxyyy_0[i] = 2.0 * g_yzzzz_0_xyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyz_0[i] = 2.0 * g_yzzzz_0_xyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyzz_0[i] = 2.0 * g_yzzzz_0_xyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxzzz_0[i] = g_xzzzz_0_xxzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xyyyy_0[i] = g_yzzzz_0_yyyy_1[i] * fi_acd_0 + g_yzzzz_0_xyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyz_0[i] = g_yzzzz_0_yyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyzz_0[i] = g_yzzzz_0_yyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyzzz_0[i] = g_yzzzz_0_yzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xzzzz_0[i] = g_xzzzz_0_xzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_yyyyy_0[i] = g_yzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyz_0[i] = g_yzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyzz_0[i] = g_yzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyzzz_0[i] = g_yzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yzzzz_0[i] = g_yzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_zzzzz_0[i] = g_yzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 420-441 components of targeted buffer : ISH

    auto g_xzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 420);

    auto g_xzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 421);

    auto g_xzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 422);

    auto g_xzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 423);

    auto g_xzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 424);

    auto g_xzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 425);

    auto g_xzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 426);

    auto g_xzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 427);

    auto g_xzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 428);

    auto g_xzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 429);

    auto g_xzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 430);

    auto g_xzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 431);

    auto g_xzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 432);

    auto g_xzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 433);

    auto g_xzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 434);

    auto g_xzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 435);

    auto g_xzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 436);

    auto g_xzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 437);

    auto g_xzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 438);

    auto g_xzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 439);

    auto g_xzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 440);

    #pragma omp simd aligned(g_xzzzzz_0_xxxxx_0, g_xzzzzz_0_xxxxy_0, g_xzzzzz_0_xxxxz_0, g_xzzzzz_0_xxxyy_0, g_xzzzzz_0_xxxyz_0, g_xzzzzz_0_xxxzz_0, g_xzzzzz_0_xxyyy_0, g_xzzzzz_0_xxyyz_0, g_xzzzzz_0_xxyzz_0, g_xzzzzz_0_xxzzz_0, g_xzzzzz_0_xyyyy_0, g_xzzzzz_0_xyyyz_0, g_xzzzzz_0_xyyzz_0, g_xzzzzz_0_xyzzz_0, g_xzzzzz_0_xzzzz_0, g_xzzzzz_0_yyyyy_0, g_xzzzzz_0_yyyyz_0, g_xzzzzz_0_yyyzz_0, g_xzzzzz_0_yyzzz_0, g_xzzzzz_0_yzzzz_0, g_xzzzzz_0_zzzzz_0, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxzz_1, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyzz_1, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xzzz_1, g_zzzzz_0_xzzzz_1, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyzz_1, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yzzz_1, g_zzzzz_0_yzzzz_1, g_zzzzz_0_zzzz_1, g_zzzzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xxxxx_0[i] = 5.0 * g_zzzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxx_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxy_0[i] = 4.0 * g_zzzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxz_0[i] = 4.0 * g_zzzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyy_0[i] = 3.0 * g_zzzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyz_0[i] = 3.0 * g_zzzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxzz_0[i] = 3.0 * g_zzzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyy_0[i] = 2.0 * g_zzzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyz_0[i] = 2.0 * g_zzzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyzz_0[i] = 2.0 * g_zzzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxzzz_0[i] = 2.0 * g_zzzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyy_0[i] = g_zzzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyz_0[i] = g_zzzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyzz_0[i] = g_zzzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyzzz_0[i] = g_zzzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xzzzz_0[i] = g_zzzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyy_0[i] = g_zzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyz_0[i] = g_zzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyzz_0[i] = g_zzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyzzz_0[i] = g_zzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yzzzz_0[i] = g_zzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_zzzzz_0[i] = g_zzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 441-462 components of targeted buffer : ISH

    auto g_yyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 441);

    auto g_yyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 442);

    auto g_yyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 443);

    auto g_yyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 444);

    auto g_yyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 445);

    auto g_yyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 446);

    auto g_yyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 447);

    auto g_yyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 448);

    auto g_yyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 449);

    auto g_yyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 450);

    auto g_yyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 451);

    auto g_yyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 452);

    auto g_yyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 453);

    auto g_yyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 454);

    auto g_yyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 455);

    auto g_yyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 456);

    auto g_yyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 457);

    auto g_yyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 458);

    auto g_yyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 459);

    auto g_yyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 460);

    auto g_yyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 461);

    #pragma omp simd aligned(g_yyyy_0_xxxxx_0, g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxy_0, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxz_0, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxyy_0, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyz_0, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxzz_0, g_yyyy_0_xxxzz_1, g_yyyy_0_xxyyy_0, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyz_0, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyzz_0, g_yyyy_0_xxyzz_1, g_yyyy_0_xxzzz_0, g_yyyy_0_xxzzz_1, g_yyyy_0_xyyyy_0, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyz_0, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyzz_0, g_yyyy_0_xyyzz_1, g_yyyy_0_xyzzz_0, g_yyyy_0_xyzzz_1, g_yyyy_0_xzzzz_0, g_yyyy_0_xzzzz_1, g_yyyy_0_yyyyy_0, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyz_0, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyzz_0, g_yyyy_0_yyyzz_1, g_yyyy_0_yyzzz_0, g_yyyy_0_yyzzz_1, g_yyyy_0_yzzzz_0, g_yyyy_0_yzzzz_1, g_yyyy_0_zzzzz_0, g_yyyy_0_zzzzz_1, g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxzz_1, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyzz_1, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xzzz_1, g_yyyyy_0_xzzzz_1, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyzz_1, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yzzz_1, g_yyyyy_0_yzzzz_1, g_yyyyy_0_zzzz_1, g_yyyyy_0_zzzzz_1, g_yyyyyy_0_xxxxx_0, g_yyyyyy_0_xxxxy_0, g_yyyyyy_0_xxxxz_0, g_yyyyyy_0_xxxyy_0, g_yyyyyy_0_xxxyz_0, g_yyyyyy_0_xxxzz_0, g_yyyyyy_0_xxyyy_0, g_yyyyyy_0_xxyyz_0, g_yyyyyy_0_xxyzz_0, g_yyyyyy_0_xxzzz_0, g_yyyyyy_0_xyyyy_0, g_yyyyyy_0_xyyyz_0, g_yyyyyy_0_xyyzz_0, g_yyyyyy_0_xyzzz_0, g_yyyyyy_0_xzzzz_0, g_yyyyyy_0_yyyyy_0, g_yyyyyy_0_yyyyz_0, g_yyyyyy_0_yyyzz_0, g_yyyyyy_0_yyzzz_0, g_yyyyyy_0_yzzzz_0, g_yyyyyy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xxxxx_0[i] = 5.0 * g_yyyy_0_xxxxx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxx_1[i] * fz_be_0 + g_yyyyy_0_xxxxx_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxy_0[i] = 5.0 * g_yyyy_0_xxxxy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxy_1[i] * fz_be_0 + g_yyyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxz_0[i] = 5.0 * g_yyyy_0_xxxxz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxz_1[i] * fz_be_0 + g_yyyyy_0_xxxxz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyy_0[i] = 5.0 * g_yyyy_0_xxxyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyz_0[i] = 5.0 * g_yyyy_0_xxxyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyz_1[i] * fz_be_0 + g_yyyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxzz_0[i] = 5.0 * g_yyyy_0_xxxzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyy_0[i] = 5.0 * g_yyyy_0_xxyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyz_0[i] = 5.0 * g_yyyy_0_xxyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyzz_0[i] = 5.0 * g_yyyy_0_xxyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyzz_1[i] * fz_be_0 + g_yyyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxzzz_0[i] = 5.0 * g_yyyy_0_xxzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyy_0[i] = 5.0 * g_yyyy_0_xyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyz_0[i] = 5.0 * g_yyyy_0_xyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyzz_0[i] = 5.0 * g_yyyy_0_xyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyzzz_0[i] = 5.0 * g_yyyy_0_xyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xzzzz_0[i] = 5.0 * g_yyyy_0_xzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyy_0[i] = 5.0 * g_yyyy_0_yyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyz_0[i] = 5.0 * g_yyyy_0_yyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyzz_0[i] = 5.0 * g_yyyy_0_yyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyzzz_0[i] = 5.0 * g_yyyy_0_yyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yzzzz_0[i] = 5.0 * g_yyyy_0_yzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_zzzzz_0[i] = 5.0 * g_yyyy_0_zzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 462-483 components of targeted buffer : ISH

    auto g_yyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 462);

    auto g_yyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 463);

    auto g_yyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 464);

    auto g_yyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 465);

    auto g_yyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 466);

    auto g_yyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 467);

    auto g_yyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 468);

    auto g_yyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 469);

    auto g_yyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 470);

    auto g_yyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 471);

    auto g_yyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 472);

    auto g_yyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 473);

    auto g_yyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 474);

    auto g_yyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 475);

    auto g_yyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 476);

    auto g_yyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 477);

    auto g_yyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 478);

    auto g_yyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 479);

    auto g_yyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 480);

    auto g_yyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 481);

    auto g_yyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 482);

    #pragma omp simd aligned(g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxzz_1, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyzz_1, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xzzz_1, g_yyyyy_0_xzzzz_1, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyzz_1, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yzzz_1, g_yyyyy_0_yzzzz_1, g_yyyyy_0_zzzz_1, g_yyyyy_0_zzzzz_1, g_yyyyyz_0_xxxxx_0, g_yyyyyz_0_xxxxy_0, g_yyyyyz_0_xxxxz_0, g_yyyyyz_0_xxxyy_0, g_yyyyyz_0_xxxyz_0, g_yyyyyz_0_xxxzz_0, g_yyyyyz_0_xxyyy_0, g_yyyyyz_0_xxyyz_0, g_yyyyyz_0_xxyzz_0, g_yyyyyz_0_xxzzz_0, g_yyyyyz_0_xyyyy_0, g_yyyyyz_0_xyyyz_0, g_yyyyyz_0_xyyzz_0, g_yyyyyz_0_xyzzz_0, g_yyyyyz_0_xzzzz_0, g_yyyyyz_0_yyyyy_0, g_yyyyyz_0_yyyyz_0, g_yyyyyz_0_yyyzz_0, g_yyyyyz_0_yyzzz_0, g_yyyyyz_0_yzzzz_0, g_yyyyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xxxxx_0[i] = g_yyyyy_0_xxxxx_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxy_0[i] = g_yyyyy_0_xxxxy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxz_0[i] = g_yyyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyy_0[i] = g_yyyyy_0_xxxyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyz_0[i] = g_yyyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxzz_0[i] = 2.0 * g_yyyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyy_0[i] = g_yyyyy_0_xxyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyz_0[i] = g_yyyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyzz_0[i] = 2.0 * g_yyyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxzzz_0[i] = 3.0 * g_yyyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyy_0[i] = g_yyyyy_0_xyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyz_0[i] = g_yyyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyzz_0[i] = 2.0 * g_yyyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyzzz_0[i] = 3.0 * g_yyyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xzzzz_0[i] = 4.0 * g_yyyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyy_0[i] = g_yyyyy_0_yyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyz_0[i] = g_yyyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyzz_0[i] = 2.0 * g_yyyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyzzz_0[i] = 3.0 * g_yyyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yzzzz_0[i] = 4.0 * g_yyyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_zzzzz_0[i] = 5.0 * g_yyyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 483-504 components of targeted buffer : ISH

    auto g_yyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 483);

    auto g_yyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 484);

    auto g_yyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 485);

    auto g_yyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 486);

    auto g_yyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 487);

    auto g_yyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 488);

    auto g_yyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 489);

    auto g_yyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 490);

    auto g_yyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 491);

    auto g_yyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 492);

    auto g_yyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 493);

    auto g_yyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 494);

    auto g_yyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 495);

    auto g_yyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 496);

    auto g_yyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 497);

    auto g_yyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 498);

    auto g_yyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 499);

    auto g_yyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 500);

    auto g_yyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 501);

    auto g_yyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 502);

    auto g_yyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 503);

    #pragma omp simd aligned(g_yyyy_0_xxxxy_0, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxyy_0, g_yyyy_0_xxxyy_1, g_yyyy_0_xxyyy_0, g_yyyy_0_xxyyy_1, g_yyyy_0_xyyyy_0, g_yyyy_0_xyyyy_1, g_yyyy_0_yyyyy_0, g_yyyy_0_yyyyy_1, g_yyyyz_0_xxxxy_1, g_yyyyz_0_xxxyy_1, g_yyyyz_0_xxyyy_1, g_yyyyz_0_xyyyy_1, g_yyyyz_0_yyyyy_1, g_yyyyzz_0_xxxxx_0, g_yyyyzz_0_xxxxy_0, g_yyyyzz_0_xxxxz_0, g_yyyyzz_0_xxxyy_0, g_yyyyzz_0_xxxyz_0, g_yyyyzz_0_xxxzz_0, g_yyyyzz_0_xxyyy_0, g_yyyyzz_0_xxyyz_0, g_yyyyzz_0_xxyzz_0, g_yyyyzz_0_xxzzz_0, g_yyyyzz_0_xyyyy_0, g_yyyyzz_0_xyyyz_0, g_yyyyzz_0_xyyzz_0, g_yyyyzz_0_xyzzz_0, g_yyyyzz_0_xzzzz_0, g_yyyyzz_0_yyyyy_0, g_yyyyzz_0_yyyyz_0, g_yyyyzz_0_yyyzz_0, g_yyyyzz_0_yyzzz_0, g_yyyyzz_0_yzzzz_0, g_yyyyzz_0_zzzzz_0, g_yyyzz_0_xxxxx_1, g_yyyzz_0_xxxxz_1, g_yyyzz_0_xxxyz_1, g_yyyzz_0_xxxz_1, g_yyyzz_0_xxxzz_1, g_yyyzz_0_xxyyz_1, g_yyyzz_0_xxyz_1, g_yyyzz_0_xxyzz_1, g_yyyzz_0_xxzz_1, g_yyyzz_0_xxzzz_1, g_yyyzz_0_xyyyz_1, g_yyyzz_0_xyyz_1, g_yyyzz_0_xyyzz_1, g_yyyzz_0_xyzz_1, g_yyyzz_0_xyzzz_1, g_yyyzz_0_xzzz_1, g_yyyzz_0_xzzzz_1, g_yyyzz_0_yyyyz_1, g_yyyzz_0_yyyz_1, g_yyyzz_0_yyyzz_1, g_yyyzz_0_yyzz_1, g_yyyzz_0_yyzzz_1, g_yyyzz_0_yzzz_1, g_yyyzz_0_yzzzz_1, g_yyyzz_0_zzzz_1, g_yyyzz_0_zzzzz_1, g_yyzz_0_xxxxx_0, g_yyzz_0_xxxxx_1, g_yyzz_0_xxxxz_0, g_yyzz_0_xxxxz_1, g_yyzz_0_xxxyz_0, g_yyzz_0_xxxyz_1, g_yyzz_0_xxxzz_0, g_yyzz_0_xxxzz_1, g_yyzz_0_xxyyz_0, g_yyzz_0_xxyyz_1, g_yyzz_0_xxyzz_0, g_yyzz_0_xxyzz_1, g_yyzz_0_xxzzz_0, g_yyzz_0_xxzzz_1, g_yyzz_0_xyyyz_0, g_yyzz_0_xyyyz_1, g_yyzz_0_xyyzz_0, g_yyzz_0_xyyzz_1, g_yyzz_0_xyzzz_0, g_yyzz_0_xyzzz_1, g_yyzz_0_xzzzz_0, g_yyzz_0_xzzzz_1, g_yyzz_0_yyyyz_0, g_yyzz_0_yyyyz_1, g_yyzz_0_yyyzz_0, g_yyzz_0_yyyzz_1, g_yyzz_0_yyzzz_0, g_yyzz_0_yyzzz_1, g_yyzz_0_yzzzz_0, g_yyzz_0_yzzzz_1, g_yyzz_0_zzzzz_0, g_yyzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xxxxx_0[i] = 3.0 * g_yyzz_0_xxxxx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxx_1[i] * fz_be_0 + g_yyyzz_0_xxxxx_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxy_0[i] = g_yyyy_0_xxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxy_1[i] * fz_be_0 + g_yyyyz_0_xxxxy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxz_0[i] = 3.0 * g_yyzz_0_xxxxz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxz_1[i] * fz_be_0 + g_yyyzz_0_xxxxz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyy_0[i] = g_yyyy_0_xxxyy_0[i] * fbe_0 - g_yyyy_0_xxxyy_1[i] * fz_be_0 + g_yyyyz_0_xxxyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxyz_0[i] = 3.0 * g_yyzz_0_xxxyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyz_1[i] * fz_be_0 + g_yyyzz_0_xxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxzz_0[i] = 3.0 * g_yyzz_0_xxxzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyy_0[i] = g_yyyy_0_xxyyy_0[i] * fbe_0 - g_yyyy_0_xxyyy_1[i] * fz_be_0 + g_yyyyz_0_xxyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxyyz_0[i] = 3.0 * g_yyzz_0_xxyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyzz_0[i] = 3.0 * g_yyzz_0_xxyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyzz_1[i] * fz_be_0 + g_yyyzz_0_xxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxzzz_0[i] = 3.0 * g_yyzz_0_xxzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyy_0[i] = g_yyyy_0_xyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyy_1[i] * fz_be_0 + g_yyyyz_0_xyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xyyyz_0[i] = 3.0 * g_yyzz_0_xyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyzz_0[i] = 3.0 * g_yyzz_0_xyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyzzz_0[i] = 3.0 * g_yyzz_0_xyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xzzzz_0[i] = 3.0 * g_yyzz_0_xzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyy_0[i] = g_yyyy_0_yyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyy_1[i] * fz_be_0 + g_yyyyz_0_yyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_yyyyz_0[i] = 3.0 * g_yyzz_0_yyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_yyyz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyzz_0[i] = 3.0 * g_yyzz_0_yyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_yyzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyzzz_0[i] = 3.0 * g_yyzz_0_yyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_yzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yzzzz_0[i] = 3.0 * g_yyzz_0_yzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzz_1[i] * fi_acd_0 + g_yyyzz_0_yzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_zzzzz_0[i] = 3.0 * g_yyzz_0_zzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 504-525 components of targeted buffer : ISH

    auto g_yyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 504);

    auto g_yyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 505);

    auto g_yyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 506);

    auto g_yyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 507);

    auto g_yyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 508);

    auto g_yyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 509);

    auto g_yyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 510);

    auto g_yyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 511);

    auto g_yyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 512);

    auto g_yyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 513);

    auto g_yyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 514);

    auto g_yyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 515);

    auto g_yyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 516);

    auto g_yyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 517);

    auto g_yyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 518);

    auto g_yyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 519);

    auto g_yyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 520);

    auto g_yyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 521);

    auto g_yyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 522);

    auto g_yyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 523);

    auto g_yyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 524);

    #pragma omp simd aligned(g_yyyz_0_xxxxy_0, g_yyyz_0_xxxxy_1, g_yyyz_0_xxxyy_0, g_yyyz_0_xxxyy_1, g_yyyz_0_xxyyy_0, g_yyyz_0_xxyyy_1, g_yyyz_0_xyyyy_0, g_yyyz_0_xyyyy_1, g_yyyz_0_yyyyy_0, g_yyyz_0_yyyyy_1, g_yyyzz_0_xxxxy_1, g_yyyzz_0_xxxyy_1, g_yyyzz_0_xxyyy_1, g_yyyzz_0_xyyyy_1, g_yyyzz_0_yyyyy_1, g_yyyzzz_0_xxxxx_0, g_yyyzzz_0_xxxxy_0, g_yyyzzz_0_xxxxz_0, g_yyyzzz_0_xxxyy_0, g_yyyzzz_0_xxxyz_0, g_yyyzzz_0_xxxzz_0, g_yyyzzz_0_xxyyy_0, g_yyyzzz_0_xxyyz_0, g_yyyzzz_0_xxyzz_0, g_yyyzzz_0_xxzzz_0, g_yyyzzz_0_xyyyy_0, g_yyyzzz_0_xyyyz_0, g_yyyzzz_0_xyyzz_0, g_yyyzzz_0_xyzzz_0, g_yyyzzz_0_xzzzz_0, g_yyyzzz_0_yyyyy_0, g_yyyzzz_0_yyyyz_0, g_yyyzzz_0_yyyzz_0, g_yyyzzz_0_yyzzz_0, g_yyyzzz_0_yzzzz_0, g_yyyzzz_0_zzzzz_0, g_yyzzz_0_xxxxx_1, g_yyzzz_0_xxxxz_1, g_yyzzz_0_xxxyz_1, g_yyzzz_0_xxxz_1, g_yyzzz_0_xxxzz_1, g_yyzzz_0_xxyyz_1, g_yyzzz_0_xxyz_1, g_yyzzz_0_xxyzz_1, g_yyzzz_0_xxzz_1, g_yyzzz_0_xxzzz_1, g_yyzzz_0_xyyyz_1, g_yyzzz_0_xyyz_1, g_yyzzz_0_xyyzz_1, g_yyzzz_0_xyzz_1, g_yyzzz_0_xyzzz_1, g_yyzzz_0_xzzz_1, g_yyzzz_0_xzzzz_1, g_yyzzz_0_yyyyz_1, g_yyzzz_0_yyyz_1, g_yyzzz_0_yyyzz_1, g_yyzzz_0_yyzz_1, g_yyzzz_0_yyzzz_1, g_yyzzz_0_yzzz_1, g_yyzzz_0_yzzzz_1, g_yyzzz_0_zzzz_1, g_yyzzz_0_zzzzz_1, g_yzzz_0_xxxxx_0, g_yzzz_0_xxxxx_1, g_yzzz_0_xxxxz_0, g_yzzz_0_xxxxz_1, g_yzzz_0_xxxyz_0, g_yzzz_0_xxxyz_1, g_yzzz_0_xxxzz_0, g_yzzz_0_xxxzz_1, g_yzzz_0_xxyyz_0, g_yzzz_0_xxyyz_1, g_yzzz_0_xxyzz_0, g_yzzz_0_xxyzz_1, g_yzzz_0_xxzzz_0, g_yzzz_0_xxzzz_1, g_yzzz_0_xyyyz_0, g_yzzz_0_xyyyz_1, g_yzzz_0_xyyzz_0, g_yzzz_0_xyyzz_1, g_yzzz_0_xyzzz_0, g_yzzz_0_xyzzz_1, g_yzzz_0_xzzzz_0, g_yzzz_0_xzzzz_1, g_yzzz_0_yyyyz_0, g_yzzz_0_yyyyz_1, g_yzzz_0_yyyzz_0, g_yzzz_0_yyyzz_1, g_yzzz_0_yyzzz_0, g_yzzz_0_yyzzz_1, g_yzzz_0_yzzzz_0, g_yzzz_0_yzzzz_1, g_yzzz_0_zzzzz_0, g_yzzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xxxxx_0[i] = 2.0 * g_yzzz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxx_1[i] * fz_be_0 + g_yyzzz_0_xxxxx_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxy_0[i] = 2.0 * g_yyyz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxy_1[i] * fz_be_0 + g_yyyzz_0_xxxxy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxz_0[i] = 2.0 * g_yzzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxz_1[i] * fz_be_0 + g_yyzzz_0_xxxxz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyy_0[i] = 2.0 * g_yyyz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxyy_1[i] * fz_be_0 + g_yyyzz_0_xxxyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxyz_0[i] = 2.0 * g_yzzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyz_1[i] * fz_be_0 + g_yyzzz_0_xxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxzz_0[i] = 2.0 * g_yzzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyy_0[i] = 2.0 * g_yyyz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxyyy_1[i] * fz_be_0 + g_yyyzz_0_xxyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxyyz_0[i] = 2.0 * g_yzzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyzz_0[i] = 2.0 * g_yzzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyzz_1[i] * fz_be_0 + g_yyzzz_0_xxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxzzz_0[i] = 2.0 * g_yzzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyy_0[i] = 2.0 * g_yyyz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xyyyy_1[i] * fz_be_0 + g_yyyzz_0_xyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xyyyz_0[i] = 2.0 * g_yzzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyzz_0[i] = 2.0 * g_yzzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyzzz_0[i] = 2.0 * g_yzzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xzzzz_0[i] = 2.0 * g_yzzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyy_0[i] = 2.0 * g_yyyz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yyyyy_1[i] * fz_be_0 + g_yyyzz_0_yyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_yyyyz_0[i] = 2.0 * g_yzzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_yyyz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyzz_0[i] = 2.0 * g_yzzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_yyzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyzzz_0[i] = 2.0 * g_yzzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_yzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yzzzz_0[i] = 2.0 * g_yzzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzz_1[i] * fi_acd_0 + g_yyzzz_0_yzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_zzzzz_0[i] = 2.0 * g_yzzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 525-546 components of targeted buffer : ISH

    auto g_yyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 525);

    auto g_yyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 526);

    auto g_yyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 527);

    auto g_yyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 528);

    auto g_yyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 529);

    auto g_yyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 530);

    auto g_yyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 531);

    auto g_yyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 532);

    auto g_yyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 533);

    auto g_yyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 534);

    auto g_yyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 535);

    auto g_yyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 536);

    auto g_yyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 537);

    auto g_yyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 538);

    auto g_yyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 539);

    auto g_yyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 540);

    auto g_yyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 541);

    auto g_yyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 542);

    auto g_yyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 543);

    auto g_yyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 544);

    auto g_yyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 545);

    #pragma omp simd aligned(g_yyzz_0_xxxxy_0, g_yyzz_0_xxxxy_1, g_yyzz_0_xxxyy_0, g_yyzz_0_xxxyy_1, g_yyzz_0_xxyyy_0, g_yyzz_0_xxyyy_1, g_yyzz_0_xyyyy_0, g_yyzz_0_xyyyy_1, g_yyzz_0_yyyyy_0, g_yyzz_0_yyyyy_1, g_yyzzz_0_xxxxy_1, g_yyzzz_0_xxxyy_1, g_yyzzz_0_xxyyy_1, g_yyzzz_0_xyyyy_1, g_yyzzz_0_yyyyy_1, g_yyzzzz_0_xxxxx_0, g_yyzzzz_0_xxxxy_0, g_yyzzzz_0_xxxxz_0, g_yyzzzz_0_xxxyy_0, g_yyzzzz_0_xxxyz_0, g_yyzzzz_0_xxxzz_0, g_yyzzzz_0_xxyyy_0, g_yyzzzz_0_xxyyz_0, g_yyzzzz_0_xxyzz_0, g_yyzzzz_0_xxzzz_0, g_yyzzzz_0_xyyyy_0, g_yyzzzz_0_xyyyz_0, g_yyzzzz_0_xyyzz_0, g_yyzzzz_0_xyzzz_0, g_yyzzzz_0_xzzzz_0, g_yyzzzz_0_yyyyy_0, g_yyzzzz_0_yyyyz_0, g_yyzzzz_0_yyyzz_0, g_yyzzzz_0_yyzzz_0, g_yyzzzz_0_yzzzz_0, g_yyzzzz_0_zzzzz_0, g_yzzzz_0_xxxxx_1, g_yzzzz_0_xxxxz_1, g_yzzzz_0_xxxyz_1, g_yzzzz_0_xxxz_1, g_yzzzz_0_xxxzz_1, g_yzzzz_0_xxyyz_1, g_yzzzz_0_xxyz_1, g_yzzzz_0_xxyzz_1, g_yzzzz_0_xxzz_1, g_yzzzz_0_xxzzz_1, g_yzzzz_0_xyyyz_1, g_yzzzz_0_xyyz_1, g_yzzzz_0_xyyzz_1, g_yzzzz_0_xyzz_1, g_yzzzz_0_xyzzz_1, g_yzzzz_0_xzzz_1, g_yzzzz_0_xzzzz_1, g_yzzzz_0_yyyyz_1, g_yzzzz_0_yyyz_1, g_yzzzz_0_yyyzz_1, g_yzzzz_0_yyzz_1, g_yzzzz_0_yyzzz_1, g_yzzzz_0_yzzz_1, g_yzzzz_0_yzzzz_1, g_yzzzz_0_zzzz_1, g_yzzzz_0_zzzzz_1, g_zzzz_0_xxxxx_0, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxz_0, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxyz_0, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxzz_0, g_zzzz_0_xxxzz_1, g_zzzz_0_xxyyz_0, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyzz_0, g_zzzz_0_xxyzz_1, g_zzzz_0_xxzzz_0, g_zzzz_0_xxzzz_1, g_zzzz_0_xyyyz_0, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyzz_0, g_zzzz_0_xyyzz_1, g_zzzz_0_xyzzz_0, g_zzzz_0_xyzzz_1, g_zzzz_0_xzzzz_0, g_zzzz_0_xzzzz_1, g_zzzz_0_yyyyz_0, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyzz_0, g_zzzz_0_yyyzz_1, g_zzzz_0_yyzzz_0, g_zzzz_0_yyzzz_1, g_zzzz_0_yzzzz_0, g_zzzz_0_yzzzz_1, g_zzzz_0_zzzzz_0, g_zzzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xxxxx_0[i] = g_zzzz_0_xxxxx_0[i] * fbe_0 - g_zzzz_0_xxxxx_1[i] * fz_be_0 + g_yzzzz_0_xxxxx_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxy_0[i] = 3.0 * g_yyzz_0_xxxxy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxy_1[i] * fz_be_0 + g_yyzzz_0_xxxxy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxz_0[i] = g_zzzz_0_xxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxz_1[i] * fz_be_0 + g_yzzzz_0_xxxxz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyy_0[i] = 3.0 * g_yyzz_0_xxxyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyy_1[i] * fz_be_0 + g_yyzzz_0_xxxyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxyz_0[i] = g_zzzz_0_xxxyz_0[i] * fbe_0 - g_zzzz_0_xxxyz_1[i] * fz_be_0 + g_yzzzz_0_xxxz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxzz_0[i] = g_zzzz_0_xxxzz_0[i] * fbe_0 - g_zzzz_0_xxxzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyy_0[i] = 3.0 * g_yyzz_0_xxyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyy_1[i] * fz_be_0 + g_yyzzz_0_xxyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxyyz_0[i] = g_zzzz_0_xxyyz_0[i] * fbe_0 - g_zzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyzz_0[i] = g_zzzz_0_xxyzz_0[i] * fbe_0 - g_zzzz_0_xxyzz_1[i] * fz_be_0 + g_yzzzz_0_xxzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxzzz_0[i] = g_zzzz_0_xxzzz_0[i] * fbe_0 - g_zzzz_0_xxzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyy_0[i] = 3.0 * g_yyzz_0_xyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyy_1[i] * fz_be_0 + g_yyzzz_0_xyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xyyyz_0[i] = g_zzzz_0_xyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyzz_0[i] = g_zzzz_0_xyyzz_0[i] * fbe_0 - g_zzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyzzz_0[i] = g_zzzz_0_xyzzz_0[i] * fbe_0 - g_zzzz_0_xyzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xzzzz_0[i] = g_zzzz_0_xzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyy_0[i] = 3.0 * g_yyzz_0_yyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyy_1[i] * fz_be_0 + g_yyzzz_0_yyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_yyyyz_0[i] = g_zzzz_0_yyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_yyyz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyzz_0[i] = g_zzzz_0_yyyzz_0[i] * fbe_0 - g_zzzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_yyzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyzzz_0[i] = g_zzzz_0_yyzzz_0[i] * fbe_0 - g_zzzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_yzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yzzzz_0[i] = g_zzzz_0_yzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzz_1[i] * fi_acd_0 + g_yzzzz_0_yzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_zzzzz_0[i] = g_zzzz_0_zzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 546-567 components of targeted buffer : ISH

    auto g_yzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 546);

    auto g_yzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 547);

    auto g_yzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 548);

    auto g_yzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 549);

    auto g_yzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 550);

    auto g_yzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 551);

    auto g_yzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 552);

    auto g_yzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 553);

    auto g_yzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 554);

    auto g_yzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 555);

    auto g_yzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 556);

    auto g_yzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 557);

    auto g_yzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 558);

    auto g_yzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 559);

    auto g_yzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 560);

    auto g_yzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 561);

    auto g_yzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 562);

    auto g_yzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 563);

    auto g_yzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 564);

    auto g_yzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 565);

    auto g_yzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 566);

    #pragma omp simd aligned(g_yzzzzz_0_xxxxx_0, g_yzzzzz_0_xxxxy_0, g_yzzzzz_0_xxxxz_0, g_yzzzzz_0_xxxyy_0, g_yzzzzz_0_xxxyz_0, g_yzzzzz_0_xxxzz_0, g_yzzzzz_0_xxyyy_0, g_yzzzzz_0_xxyyz_0, g_yzzzzz_0_xxyzz_0, g_yzzzzz_0_xxzzz_0, g_yzzzzz_0_xyyyy_0, g_yzzzzz_0_xyyyz_0, g_yzzzzz_0_xyyzz_0, g_yzzzzz_0_xyzzz_0, g_yzzzzz_0_xzzzz_0, g_yzzzzz_0_yyyyy_0, g_yzzzzz_0_yyyyz_0, g_yzzzzz_0_yyyzz_0, g_yzzzzz_0_yyzzz_0, g_yzzzzz_0_yzzzz_0, g_yzzzzz_0_zzzzz_0, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxzz_1, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyzz_1, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xzzz_1, g_zzzzz_0_xzzzz_1, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyzz_1, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yzzz_1, g_zzzzz_0_yzzzz_1, g_zzzzz_0_zzzz_1, g_zzzzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xxxxx_0[i] = g_zzzzz_0_xxxxx_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxy_0[i] = g_zzzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxz_0[i] = g_zzzzz_0_xxxxz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyy_0[i] = 2.0 * g_zzzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyz_0[i] = g_zzzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxzz_0[i] = g_zzzzz_0_xxxzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyy_0[i] = 3.0 * g_zzzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyz_0[i] = 2.0 * g_zzzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyzz_0[i] = g_zzzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxzzz_0[i] = g_zzzzz_0_xxzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyy_0[i] = 4.0 * g_zzzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyz_0[i] = 3.0 * g_zzzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyzz_0[i] = 2.0 * g_zzzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyzzz_0[i] = g_zzzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xzzzz_0[i] = g_zzzzz_0_xzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyy_0[i] = 5.0 * g_zzzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyz_0[i] = 4.0 * g_zzzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyzz_0[i] = 3.0 * g_zzzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyzzz_0[i] = 2.0 * g_zzzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yzzzz_0[i] = g_zzzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_zzzzz_0[i] = g_zzzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 567-588 components of targeted buffer : ISH

    auto g_zzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ish + 567);

    auto g_zzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ish + 568);

    auto g_zzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ish + 569);

    auto g_zzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ish + 570);

    auto g_zzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ish + 571);

    auto g_zzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ish + 572);

    auto g_zzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ish + 573);

    auto g_zzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ish + 574);

    auto g_zzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ish + 575);

    auto g_zzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ish + 576);

    auto g_zzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ish + 577);

    auto g_zzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ish + 578);

    auto g_zzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ish + 579);

    auto g_zzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ish + 580);

    auto g_zzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ish + 581);

    auto g_zzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ish + 582);

    auto g_zzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ish + 583);

    auto g_zzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ish + 584);

    auto g_zzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ish + 585);

    auto g_zzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ish + 586);

    auto g_zzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ish + 587);

    #pragma omp simd aligned(g_zzzz_0_xxxxx_0, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxy_0, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxz_0, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxyy_0, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyz_0, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxzz_0, g_zzzz_0_xxxzz_1, g_zzzz_0_xxyyy_0, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyz_0, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyzz_0, g_zzzz_0_xxyzz_1, g_zzzz_0_xxzzz_0, g_zzzz_0_xxzzz_1, g_zzzz_0_xyyyy_0, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyz_0, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyzz_0, g_zzzz_0_xyyzz_1, g_zzzz_0_xyzzz_0, g_zzzz_0_xyzzz_1, g_zzzz_0_xzzzz_0, g_zzzz_0_xzzzz_1, g_zzzz_0_yyyyy_0, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyz_0, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyzz_0, g_zzzz_0_yyyzz_1, g_zzzz_0_yyzzz_0, g_zzzz_0_yyzzz_1, g_zzzz_0_yzzzz_0, g_zzzz_0_yzzzz_1, g_zzzz_0_zzzzz_0, g_zzzz_0_zzzzz_1, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxzz_1, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyzz_1, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xzzz_1, g_zzzzz_0_xzzzz_1, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyzz_1, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yzzz_1, g_zzzzz_0_yzzzz_1, g_zzzzz_0_zzzz_1, g_zzzzz_0_zzzzz_1, g_zzzzzz_0_xxxxx_0, g_zzzzzz_0_xxxxy_0, g_zzzzzz_0_xxxxz_0, g_zzzzzz_0_xxxyy_0, g_zzzzzz_0_xxxyz_0, g_zzzzzz_0_xxxzz_0, g_zzzzzz_0_xxyyy_0, g_zzzzzz_0_xxyyz_0, g_zzzzzz_0_xxyzz_0, g_zzzzzz_0_xxzzz_0, g_zzzzzz_0_xyyyy_0, g_zzzzzz_0_xyyyz_0, g_zzzzzz_0_xyyzz_0, g_zzzzzz_0_xyzzz_0, g_zzzzzz_0_xzzzz_0, g_zzzzzz_0_yyyyy_0, g_zzzzzz_0_yyyyz_0, g_zzzzzz_0_yyyzz_0, g_zzzzzz_0_yyzzz_0, g_zzzzzz_0_yzzzz_0, g_zzzzzz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xxxxx_0[i] = 5.0 * g_zzzz_0_xxxxx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxx_1[i] * fz_be_0 + g_zzzzz_0_xxxxx_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxy_0[i] = 5.0 * g_zzzz_0_xxxxy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxy_1[i] * fz_be_0 + g_zzzzz_0_xxxxy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxz_0[i] = 5.0 * g_zzzz_0_xxxxz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxz_1[i] * fz_be_0 + g_zzzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyy_0[i] = 5.0 * g_zzzz_0_xxxyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyy_1[i] * fz_be_0 + g_zzzzz_0_xxxyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyz_0[i] = 5.0 * g_zzzz_0_xxxyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyz_1[i] * fz_be_0 + g_zzzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxzz_0[i] = 5.0 * g_zzzz_0_xxxzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyy_0[i] = 5.0 * g_zzzz_0_xxyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyy_1[i] * fz_be_0 + g_zzzzz_0_xxyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyz_0[i] = 5.0 * g_zzzz_0_xxyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyz_1[i] * fz_be_0 + g_zzzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyzz_0[i] = 5.0 * g_zzzz_0_xxyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxzzz_0[i] = 5.0 * g_zzzz_0_xxzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyy_0[i] = 5.0 * g_zzzz_0_xyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyy_1[i] * fz_be_0 + g_zzzzz_0_xyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyz_0[i] = 5.0 * g_zzzz_0_xyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyz_1[i] * fz_be_0 + g_zzzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyzz_0[i] = 5.0 * g_zzzz_0_xyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyzzz_0[i] = 5.0 * g_zzzz_0_xyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xzzzz_0[i] = 5.0 * g_zzzz_0_xzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyy_0[i] = 5.0 * g_zzzz_0_yyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyy_1[i] * fz_be_0 + g_zzzzz_0_yyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyz_0[i] = 5.0 * g_zzzz_0_yyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyz_1[i] * fz_be_0 + g_zzzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyzz_0[i] = 5.0 * g_zzzz_0_yyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyzzz_0[i] = 5.0 * g_zzzz_0_yyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yzzzz_0[i] = 5.0 * g_zzzz_0_yzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_zzzzz_0[i] = 5.0 * g_zzzz_0_zzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzzz_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

