#include "ElectronRepulsionPrimRecSISH.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sish(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sish,
                                  size_t                idx_eri_0_sgsh,
                                  size_t                idx_eri_1_sgsh,
                                  size_t                idx_eri_1_shsg,
                                  size_t                idx_eri_0_shsh,
                                  size_t                idx_eri_1_shsh,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SGSH

    auto g_0_xxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh);

    auto g_0_xxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 1);

    auto g_0_xxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 2);

    auto g_0_xxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 3);

    auto g_0_xxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 4);

    auto g_0_xxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 5);

    auto g_0_xxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 6);

    auto g_0_xxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 7);

    auto g_0_xxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 8);

    auto g_0_xxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 9);

    auto g_0_xxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 10);

    auto g_0_xxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 11);

    auto g_0_xxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 12);

    auto g_0_xxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 13);

    auto g_0_xxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 14);

    auto g_0_xxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 15);

    auto g_0_xxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 16);

    auto g_0_xxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 17);

    auto g_0_xxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 18);

    auto g_0_xxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 19);

    auto g_0_xxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 20);

    auto g_0_xxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 21);

    auto g_0_xxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 23);

    auto g_0_xxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 26);

    auto g_0_xxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 30);

    auto g_0_xxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 35);

    auto g_0_xxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 42);

    auto g_0_xxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 43);

    auto g_0_xxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 45);

    auto g_0_xxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 48);

    auto g_0_xxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 52);

    auto g_0_xxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 63);

    auto g_0_xxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 64);

    auto g_0_xxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 65);

    auto g_0_xxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 66);

    auto g_0_xxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 67);

    auto g_0_xxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 68);

    auto g_0_xxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 69);

    auto g_0_xxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 70);

    auto g_0_xxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 71);

    auto g_0_xxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 72);

    auto g_0_xxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 73);

    auto g_0_xxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 74);

    auto g_0_xxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 75);

    auto g_0_xxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 76);

    auto g_0_xxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 77);

    auto g_0_xxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 78);

    auto g_0_xxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 79);

    auto g_0_xxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 80);

    auto g_0_xxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 81);

    auto g_0_xxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 82);

    auto g_0_xxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 83);

    auto g_0_xxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 105);

    auto g_0_xxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 106);

    auto g_0_xxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 107);

    auto g_0_xxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 108);

    auto g_0_xxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 109);

    auto g_0_xxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 110);

    auto g_0_xxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 111);

    auto g_0_xxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 112);

    auto g_0_xxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 113);

    auto g_0_xxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 114);

    auto g_0_xxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 115);

    auto g_0_xxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 116);

    auto g_0_xxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 117);

    auto g_0_xxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 118);

    auto g_0_xxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 119);

    auto g_0_xxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 120);

    auto g_0_xxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 121);

    auto g_0_xxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 122);

    auto g_0_xxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 123);

    auto g_0_xxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 124);

    auto g_0_xxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 125);

    auto g_0_xyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 127);

    auto g_0_xyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 129);

    auto g_0_xyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 130);

    auto g_0_xyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 132);

    auto g_0_xyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 133);

    auto g_0_xyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 134);

    auto g_0_xyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 136);

    auto g_0_xyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 137);

    auto g_0_xyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 138);

    auto g_0_xyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 139);

    auto g_0_xyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 141);

    auto g_0_xyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 142);

    auto g_0_xyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 143);

    auto g_0_xyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 144);

    auto g_0_xyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 145);

    auto g_0_xyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 146);

    auto g_0_xzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 191);

    auto g_0_xzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 193);

    auto g_0_xzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 194);

    auto g_0_xzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 196);

    auto g_0_xzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 197);

    auto g_0_xzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 198);

    auto g_0_xzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 200);

    auto g_0_xzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 201);

    auto g_0_xzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 202);

    auto g_0_xzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 203);

    auto g_0_xzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 204);

    auto g_0_xzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 205);

    auto g_0_xzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 206);

    auto g_0_xzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 207);

    auto g_0_xzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 208);

    auto g_0_xzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 209);

    auto g_0_yyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 210);

    auto g_0_yyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 211);

    auto g_0_yyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 212);

    auto g_0_yyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 213);

    auto g_0_yyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 214);

    auto g_0_yyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 215);

    auto g_0_yyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 216);

    auto g_0_yyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 217);

    auto g_0_yyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 218);

    auto g_0_yyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 219);

    auto g_0_yyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 220);

    auto g_0_yyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 221);

    auto g_0_yyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 222);

    auto g_0_yyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 223);

    auto g_0_yyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 224);

    auto g_0_yyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 225);

    auto g_0_yyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 226);

    auto g_0_yyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 227);

    auto g_0_yyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 228);

    auto g_0_yyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 229);

    auto g_0_yyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 230);

    auto g_0_yyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 232);

    auto g_0_yyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 234);

    auto g_0_yyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 237);

    auto g_0_yyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 241);

    auto g_0_yyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 246);

    auto g_0_yyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 252);

    auto g_0_yyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 253);

    auto g_0_yyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 254);

    auto g_0_yyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 255);

    auto g_0_yyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 256);

    auto g_0_yyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 257);

    auto g_0_yyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 258);

    auto g_0_yyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 259);

    auto g_0_yyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 260);

    auto g_0_yyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 261);

    auto g_0_yyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 262);

    auto g_0_yyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 263);

    auto g_0_yyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 264);

    auto g_0_yyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 265);

    auto g_0_yyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 266);

    auto g_0_yyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 267);

    auto g_0_yyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 268);

    auto g_0_yyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 269);

    auto g_0_yyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 270);

    auto g_0_yyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 271);

    auto g_0_yyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 272);

    auto g_0_yzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 273);

    auto g_0_yzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 275);

    auto g_0_yzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 277);

    auto g_0_yzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 278);

    auto g_0_yzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 280);

    auto g_0_yzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 281);

    auto g_0_yzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 282);

    auto g_0_yzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 284);

    auto g_0_yzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 285);

    auto g_0_yzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 286);

    auto g_0_yzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 287);

    auto g_0_yzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 289);

    auto g_0_yzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 290);

    auto g_0_yzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 291);

    auto g_0_yzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 292);

    auto g_0_yzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 293);

    auto g_0_zzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 294);

    auto g_0_zzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 295);

    auto g_0_zzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 296);

    auto g_0_zzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 297);

    auto g_0_zzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 298);

    auto g_0_zzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 299);

    auto g_0_zzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 300);

    auto g_0_zzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 301);

    auto g_0_zzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 302);

    auto g_0_zzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 303);

    auto g_0_zzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 304);

    auto g_0_zzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 305);

    auto g_0_zzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 306);

    auto g_0_zzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 307);

    auto g_0_zzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 308);

    auto g_0_zzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 309);

    auto g_0_zzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 310);

    auto g_0_zzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 311);

    auto g_0_zzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 312);

    auto g_0_zzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 313);

    auto g_0_zzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 314);

    /// Set up components of auxilary buffer : SGSH

    auto g_0_xxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh);

    auto g_0_xxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 1);

    auto g_0_xxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 2);

    auto g_0_xxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 3);

    auto g_0_xxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 4);

    auto g_0_xxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 5);

    auto g_0_xxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 6);

    auto g_0_xxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 7);

    auto g_0_xxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 8);

    auto g_0_xxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 9);

    auto g_0_xxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 10);

    auto g_0_xxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 11);

    auto g_0_xxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 12);

    auto g_0_xxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 13);

    auto g_0_xxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 14);

    auto g_0_xxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 15);

    auto g_0_xxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 16);

    auto g_0_xxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 17);

    auto g_0_xxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 18);

    auto g_0_xxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 19);

    auto g_0_xxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 20);

    auto g_0_xxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 21);

    auto g_0_xxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 23);

    auto g_0_xxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 26);

    auto g_0_xxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 30);

    auto g_0_xxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 35);

    auto g_0_xxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 42);

    auto g_0_xxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 43);

    auto g_0_xxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 45);

    auto g_0_xxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 48);

    auto g_0_xxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 52);

    auto g_0_xxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 63);

    auto g_0_xxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 64);

    auto g_0_xxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 65);

    auto g_0_xxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 66);

    auto g_0_xxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 67);

    auto g_0_xxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 68);

    auto g_0_xxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 69);

    auto g_0_xxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 70);

    auto g_0_xxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 71);

    auto g_0_xxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 72);

    auto g_0_xxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 73);

    auto g_0_xxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 74);

    auto g_0_xxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 75);

    auto g_0_xxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 76);

    auto g_0_xxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 77);

    auto g_0_xxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 78);

    auto g_0_xxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 79);

    auto g_0_xxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 80);

    auto g_0_xxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 81);

    auto g_0_xxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 82);

    auto g_0_xxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 83);

    auto g_0_xxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 105);

    auto g_0_xxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 106);

    auto g_0_xxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 107);

    auto g_0_xxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 108);

    auto g_0_xxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 109);

    auto g_0_xxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 110);

    auto g_0_xxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 111);

    auto g_0_xxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 112);

    auto g_0_xxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 113);

    auto g_0_xxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 114);

    auto g_0_xxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 115);

    auto g_0_xxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 116);

    auto g_0_xxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 117);

    auto g_0_xxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 118);

    auto g_0_xxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 119);

    auto g_0_xxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 120);

    auto g_0_xxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 121);

    auto g_0_xxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 122);

    auto g_0_xxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 123);

    auto g_0_xxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 124);

    auto g_0_xxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 125);

    auto g_0_xyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 127);

    auto g_0_xyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 129);

    auto g_0_xyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 130);

    auto g_0_xyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 132);

    auto g_0_xyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 133);

    auto g_0_xyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 134);

    auto g_0_xyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 136);

    auto g_0_xyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 137);

    auto g_0_xyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 138);

    auto g_0_xyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 139);

    auto g_0_xyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 141);

    auto g_0_xyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 142);

    auto g_0_xyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 143);

    auto g_0_xyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 144);

    auto g_0_xyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 145);

    auto g_0_xyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 146);

    auto g_0_xzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 191);

    auto g_0_xzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 193);

    auto g_0_xzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 194);

    auto g_0_xzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 196);

    auto g_0_xzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 197);

    auto g_0_xzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 198);

    auto g_0_xzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 200);

    auto g_0_xzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 201);

    auto g_0_xzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 202);

    auto g_0_xzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 203);

    auto g_0_xzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 204);

    auto g_0_xzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 205);

    auto g_0_xzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 206);

    auto g_0_xzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 207);

    auto g_0_xzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 208);

    auto g_0_xzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 209);

    auto g_0_yyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 210);

    auto g_0_yyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 211);

    auto g_0_yyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 212);

    auto g_0_yyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 213);

    auto g_0_yyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 214);

    auto g_0_yyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 215);

    auto g_0_yyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 216);

    auto g_0_yyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 217);

    auto g_0_yyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 218);

    auto g_0_yyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 219);

    auto g_0_yyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 220);

    auto g_0_yyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 221);

    auto g_0_yyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 222);

    auto g_0_yyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 223);

    auto g_0_yyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 224);

    auto g_0_yyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 225);

    auto g_0_yyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 226);

    auto g_0_yyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 227);

    auto g_0_yyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 228);

    auto g_0_yyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 229);

    auto g_0_yyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 230);

    auto g_0_yyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 232);

    auto g_0_yyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 234);

    auto g_0_yyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 237);

    auto g_0_yyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 241);

    auto g_0_yyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 246);

    auto g_0_yyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 252);

    auto g_0_yyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 253);

    auto g_0_yyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 254);

    auto g_0_yyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 255);

    auto g_0_yyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 256);

    auto g_0_yyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 257);

    auto g_0_yyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 258);

    auto g_0_yyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 259);

    auto g_0_yyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 260);

    auto g_0_yyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 261);

    auto g_0_yyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 262);

    auto g_0_yyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 263);

    auto g_0_yyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 264);

    auto g_0_yyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 265);

    auto g_0_yyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 266);

    auto g_0_yyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 267);

    auto g_0_yyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 268);

    auto g_0_yyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 269);

    auto g_0_yyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 270);

    auto g_0_yyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 271);

    auto g_0_yyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 272);

    auto g_0_yzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 273);

    auto g_0_yzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 275);

    auto g_0_yzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 277);

    auto g_0_yzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 278);

    auto g_0_yzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 280);

    auto g_0_yzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 281);

    auto g_0_yzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 282);

    auto g_0_yzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 284);

    auto g_0_yzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 285);

    auto g_0_yzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 286);

    auto g_0_yzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 287);

    auto g_0_yzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 289);

    auto g_0_yzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 290);

    auto g_0_yzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 291);

    auto g_0_yzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 292);

    auto g_0_yzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 293);

    auto g_0_zzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 294);

    auto g_0_zzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 295);

    auto g_0_zzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 296);

    auto g_0_zzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 297);

    auto g_0_zzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 298);

    auto g_0_zzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 299);

    auto g_0_zzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 300);

    auto g_0_zzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 301);

    auto g_0_zzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 302);

    auto g_0_zzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 303);

    auto g_0_zzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 304);

    auto g_0_zzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 305);

    auto g_0_zzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 306);

    auto g_0_zzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 307);

    auto g_0_zzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 308);

    auto g_0_zzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 309);

    auto g_0_zzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 310);

    auto g_0_zzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 311);

    auto g_0_zzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 312);

    auto g_0_zzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 313);

    auto g_0_zzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 314);

    /// Set up components of auxilary buffer : SHSG

    auto g_0_xxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg);

    auto g_0_xxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 1);

    auto g_0_xxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 2);

    auto g_0_xxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 3);

    auto g_0_xxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 4);

    auto g_0_xxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 5);

    auto g_0_xxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 6);

    auto g_0_xxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 7);

    auto g_0_xxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 8);

    auto g_0_xxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 9);

    auto g_0_xxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 10);

    auto g_0_xxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 11);

    auto g_0_xxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 12);

    auto g_0_xxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 13);

    auto g_0_xxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 14);

    auto g_0_xxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 32);

    auto g_0_xxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 34);

    auto g_0_xxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 35);

    auto g_0_xxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 37);

    auto g_0_xxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 38);

    auto g_0_xxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 39);

    auto g_0_xxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 41);

    auto g_0_xxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 42);

    auto g_0_xxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 43);

    auto g_0_xxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 44);

    auto g_0_xxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 45);

    auto g_0_xxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 46);

    auto g_0_xxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 47);

    auto g_0_xxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 48);

    auto g_0_xxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 49);

    auto g_0_xxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 50);

    auto g_0_xxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 51);

    auto g_0_xxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 52);

    auto g_0_xxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 53);

    auto g_0_xxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 54);

    auto g_0_xxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 55);

    auto g_0_xxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 56);

    auto g_0_xxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 57);

    auto g_0_xxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 58);

    auto g_0_xxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 59);

    auto g_0_xxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 75);

    auto g_0_xxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 76);

    auto g_0_xxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 77);

    auto g_0_xxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 78);

    auto g_0_xxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 79);

    auto g_0_xxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 80);

    auto g_0_xxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 81);

    auto g_0_xxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 82);

    auto g_0_xxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 83);

    auto g_0_xxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 84);

    auto g_0_xxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 85);

    auto g_0_xxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 86);

    auto g_0_xxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 87);

    auto g_0_xxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 88);

    auto g_0_xxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 89);

    auto g_0_xxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 90);

    auto g_0_xxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 91);

    auto g_0_xxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 92);

    auto g_0_xxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 93);

    auto g_0_xxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 94);

    auto g_0_xxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 95);

    auto g_0_xxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 96);

    auto g_0_xxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 97);

    auto g_0_xxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 98);

    auto g_0_xxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 99);

    auto g_0_xxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 100);

    auto g_0_xxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 101);

    auto g_0_xxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 102);

    auto g_0_xxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 103);

    auto g_0_xxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 104);

    auto g_0_xxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 135);

    auto g_0_xxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 136);

    auto g_0_xxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 137);

    auto g_0_xxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 138);

    auto g_0_xxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 139);

    auto g_0_xxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 140);

    auto g_0_xxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 141);

    auto g_0_xxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 142);

    auto g_0_xxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 143);

    auto g_0_xxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 144);

    auto g_0_xxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 145);

    auto g_0_xxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 146);

    auto g_0_xxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 147);

    auto g_0_xxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 148);

    auto g_0_xxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 149);

    auto g_0_xyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 151);

    auto g_0_xyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 153);

    auto g_0_xyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 154);

    auto g_0_xyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 156);

    auto g_0_xyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 157);

    auto g_0_xyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 158);

    auto g_0_xyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 160);

    auto g_0_xyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 161);

    auto g_0_xyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 162);

    auto g_0_xyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 163);

    auto g_0_xyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 184);

    auto g_0_xyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 187);

    auto g_0_xyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 188);

    auto g_0_xyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 191);

    auto g_0_xyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 192);

    auto g_0_xyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 193);

    auto g_0_xzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 212);

    auto g_0_xzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 214);

    auto g_0_xzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 215);

    auto g_0_xzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 217);

    auto g_0_xzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 218);

    auto g_0_xzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 219);

    auto g_0_xzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 221);

    auto g_0_xzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 222);

    auto g_0_xzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 223);

    auto g_0_xzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 224);

    auto g_0_yyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 225);

    auto g_0_yyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 226);

    auto g_0_yyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 227);

    auto g_0_yyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 228);

    auto g_0_yyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 229);

    auto g_0_yyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 230);

    auto g_0_yyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 231);

    auto g_0_yyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 232);

    auto g_0_yyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 233);

    auto g_0_yyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 234);

    auto g_0_yyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 235);

    auto g_0_yyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 236);

    auto g_0_yyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 237);

    auto g_0_yyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 238);

    auto g_0_yyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 239);

    auto g_0_yyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 242);

    auto g_0_yyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 244);

    auto g_0_yyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 245);

    auto g_0_yyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 247);

    auto g_0_yyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 248);

    auto g_0_yyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 249);

    auto g_0_yyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 251);

    auto g_0_yyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 252);

    auto g_0_yyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 253);

    auto g_0_yyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 254);

    auto g_0_yyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 255);

    auto g_0_yyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 256);

    auto g_0_yyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 257);

    auto g_0_yyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 258);

    auto g_0_yyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 259);

    auto g_0_yyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 260);

    auto g_0_yyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 261);

    auto g_0_yyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 262);

    auto g_0_yyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 263);

    auto g_0_yyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 264);

    auto g_0_yyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 265);

    auto g_0_yyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 266);

    auto g_0_yyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 267);

    auto g_0_yyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 268);

    auto g_0_yyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 269);

    auto g_0_yyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 270);

    auto g_0_yyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 271);

    auto g_0_yyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 272);

    auto g_0_yyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 273);

    auto g_0_yyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 274);

    auto g_0_yyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 275);

    auto g_0_yyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 276);

    auto g_0_yyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 277);

    auto g_0_yyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 278);

    auto g_0_yyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 279);

    auto g_0_yyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 280);

    auto g_0_yyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 281);

    auto g_0_yyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 282);

    auto g_0_yyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 283);

    auto g_0_yyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 284);

    auto g_0_yzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 286);

    auto g_0_yzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 287);

    auto g_0_yzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 288);

    auto g_0_yzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 289);

    auto g_0_yzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 290);

    auto g_0_yzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 291);

    auto g_0_yzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 292);

    auto g_0_yzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 293);

    auto g_0_yzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 294);

    auto g_0_yzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 295);

    auto g_0_yzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 296);

    auto g_0_yzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 297);

    auto g_0_yzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 298);

    auto g_0_yzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 299);

    auto g_0_zzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 300);

    auto g_0_zzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 301);

    auto g_0_zzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 302);

    auto g_0_zzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 303);

    auto g_0_zzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 304);

    auto g_0_zzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 305);

    auto g_0_zzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 306);

    auto g_0_zzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 307);

    auto g_0_zzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 308);

    auto g_0_zzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 309);

    auto g_0_zzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 310);

    auto g_0_zzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 311);

    auto g_0_zzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 312);

    auto g_0_zzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 313);

    auto g_0_zzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 314);

    /// Set up components of auxilary buffer : SHSH

    auto g_0_xxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh);

    auto g_0_xxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 1);

    auto g_0_xxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 2);

    auto g_0_xxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 3);

    auto g_0_xxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 4);

    auto g_0_xxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 5);

    auto g_0_xxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 6);

    auto g_0_xxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 7);

    auto g_0_xxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 8);

    auto g_0_xxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 9);

    auto g_0_xxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 10);

    auto g_0_xxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 11);

    auto g_0_xxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 12);

    auto g_0_xxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 13);

    auto g_0_xxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 14);

    auto g_0_xxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 15);

    auto g_0_xxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 16);

    auto g_0_xxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 17);

    auto g_0_xxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 18);

    auto g_0_xxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 19);

    auto g_0_xxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 20);

    auto g_0_xxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 21);

    auto g_0_xxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 22);

    auto g_0_xxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 23);

    auto g_0_xxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 24);

    auto g_0_xxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 26);

    auto g_0_xxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 27);

    auto g_0_xxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 30);

    auto g_0_xxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 31);

    auto g_0_xxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 35);

    auto g_0_xxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 36);

    auto g_0_xxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 42);

    auto g_0_xxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 43);

    auto g_0_xxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 44);

    auto g_0_xxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 45);

    auto g_0_xxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 46);

    auto g_0_xxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 47);

    auto g_0_xxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 48);

    auto g_0_xxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 49);

    auto g_0_xxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 50);

    auto g_0_xxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 51);

    auto g_0_xxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 52);

    auto g_0_xxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 53);

    auto g_0_xxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 54);

    auto g_0_xxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 55);

    auto g_0_xxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 56);

    auto g_0_xxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 58);

    auto g_0_xxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 59);

    auto g_0_xxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 60);

    auto g_0_xxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 61);

    auto g_0_xxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 62);

    auto g_0_xxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 63);

    auto g_0_xxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 64);

    auto g_0_xxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 65);

    auto g_0_xxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 66);

    auto g_0_xxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 67);

    auto g_0_xxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 68);

    auto g_0_xxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 69);

    auto g_0_xxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 70);

    auto g_0_xxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 71);

    auto g_0_xxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 72);

    auto g_0_xxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 73);

    auto g_0_xxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 74);

    auto g_0_xxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 75);

    auto g_0_xxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 76);

    auto g_0_xxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 77);

    auto g_0_xxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 78);

    auto g_0_xxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 79);

    auto g_0_xxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 80);

    auto g_0_xxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 81);

    auto g_0_xxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 82);

    auto g_0_xxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 83);

    auto g_0_xxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 105);

    auto g_0_xxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 106);

    auto g_0_xxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 107);

    auto g_0_xxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 108);

    auto g_0_xxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 109);

    auto g_0_xxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 110);

    auto g_0_xxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 111);

    auto g_0_xxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 112);

    auto g_0_xxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 113);

    auto g_0_xxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 114);

    auto g_0_xxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 115);

    auto g_0_xxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 116);

    auto g_0_xxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 117);

    auto g_0_xxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 118);

    auto g_0_xxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 119);

    auto g_0_xxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 120);

    auto g_0_xxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 121);

    auto g_0_xxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 122);

    auto g_0_xxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 123);

    auto g_0_xxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 124);

    auto g_0_xxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 125);

    auto g_0_xxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 126);

    auto g_0_xxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 127);

    auto g_0_xxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 128);

    auto g_0_xxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 129);

    auto g_0_xxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 130);

    auto g_0_xxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 131);

    auto g_0_xxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 132);

    auto g_0_xxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 133);

    auto g_0_xxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 134);

    auto g_0_xxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 135);

    auto g_0_xxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 136);

    auto g_0_xxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 137);

    auto g_0_xxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 138);

    auto g_0_xxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 139);

    auto g_0_xxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 140);

    auto g_0_xxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 141);

    auto g_0_xxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 142);

    auto g_0_xxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 143);

    auto g_0_xxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 144);

    auto g_0_xxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 145);

    auto g_0_xxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 146);

    auto g_0_xxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 148);

    auto g_0_xxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 150);

    auto g_0_xxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 153);

    auto g_0_xxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 157);

    auto g_0_xxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 168);

    auto g_0_xxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 170);

    auto g_0_xxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 173);

    auto g_0_xxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 177);

    auto g_0_xxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 182);

    auto g_0_xxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 189);

    auto g_0_xxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 190);

    auto g_0_xxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 191);

    auto g_0_xxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 192);

    auto g_0_xxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 193);

    auto g_0_xxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 194);

    auto g_0_xxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 195);

    auto g_0_xxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 196);

    auto g_0_xxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 197);

    auto g_0_xxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 198);

    auto g_0_xxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 199);

    auto g_0_xxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 200);

    auto g_0_xxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 201);

    auto g_0_xxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 202);

    auto g_0_xxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 203);

    auto g_0_xxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 204);

    auto g_0_xxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 205);

    auto g_0_xxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 206);

    auto g_0_xxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 207);

    auto g_0_xxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 208);

    auto g_0_xxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 209);

    auto g_0_xyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 210);

    auto g_0_xyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 211);

    auto g_0_xyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 213);

    auto g_0_xyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 214);

    auto g_0_xyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 216);

    auto g_0_xyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 217);

    auto g_0_xyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 218);

    auto g_0_xyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 220);

    auto g_0_xyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 221);

    auto g_0_xyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 222);

    auto g_0_xyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 223);

    auto g_0_xyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 225);

    auto g_0_xyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 226);

    auto g_0_xyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 227);

    auto g_0_xyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 228);

    auto g_0_xyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 229);

    auto g_0_xyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 230);

    auto g_0_xyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 256);

    auto g_0_xyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 259);

    auto g_0_xyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 260);

    auto g_0_xyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 263);

    auto g_0_xyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 264);

    auto g_0_xyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 265);

    auto g_0_xyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 267);

    auto g_0_xyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 268);

    auto g_0_xyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 269);

    auto g_0_xyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 270);

    auto g_0_xyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 271);

    auto g_0_xyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 272);

    auto g_0_xzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 294);

    auto g_0_xzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 296);

    auto g_0_xzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 298);

    auto g_0_xzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 299);

    auto g_0_xzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 301);

    auto g_0_xzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 302);

    auto g_0_xzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 303);

    auto g_0_xzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 305);

    auto g_0_xzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 306);

    auto g_0_xzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 307);

    auto g_0_xzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 308);

    auto g_0_xzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 309);

    auto g_0_xzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 310);

    auto g_0_xzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 311);

    auto g_0_xzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 312);

    auto g_0_xzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 313);

    auto g_0_xzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 314);

    auto g_0_yyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 315);

    auto g_0_yyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 316);

    auto g_0_yyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 317);

    auto g_0_yyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 318);

    auto g_0_yyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 319);

    auto g_0_yyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 320);

    auto g_0_yyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 321);

    auto g_0_yyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 322);

    auto g_0_yyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 323);

    auto g_0_yyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 324);

    auto g_0_yyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 325);

    auto g_0_yyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 326);

    auto g_0_yyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 327);

    auto g_0_yyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 328);

    auto g_0_yyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 329);

    auto g_0_yyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 330);

    auto g_0_yyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 331);

    auto g_0_yyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 332);

    auto g_0_yyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 333);

    auto g_0_yyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 334);

    auto g_0_yyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 335);

    auto g_0_yyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 337);

    auto g_0_yyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 338);

    auto g_0_yyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 339);

    auto g_0_yyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 340);

    auto g_0_yyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 341);

    auto g_0_yyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 342);

    auto g_0_yyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 343);

    auto g_0_yyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 344);

    auto g_0_yyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 345);

    auto g_0_yyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 346);

    auto g_0_yyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 347);

    auto g_0_yyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 348);

    auto g_0_yyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 349);

    auto g_0_yyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 350);

    auto g_0_yyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 351);

    auto g_0_yyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 352);

    auto g_0_yyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 353);

    auto g_0_yyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 354);

    auto g_0_yyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 355);

    auto g_0_yyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 356);

    auto g_0_yyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 357);

    auto g_0_yyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 358);

    auto g_0_yyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 359);

    auto g_0_yyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 360);

    auto g_0_yyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 361);

    auto g_0_yyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 362);

    auto g_0_yyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 363);

    auto g_0_yyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 364);

    auto g_0_yyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 365);

    auto g_0_yyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 366);

    auto g_0_yyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 367);

    auto g_0_yyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 368);

    auto g_0_yyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 369);

    auto g_0_yyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 370);

    auto g_0_yyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 371);

    auto g_0_yyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 372);

    auto g_0_yyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 373);

    auto g_0_yyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 374);

    auto g_0_yyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 375);

    auto g_0_yyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 376);

    auto g_0_yyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 377);

    auto g_0_yyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 378);

    auto g_0_yyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 379);

    auto g_0_yyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 380);

    auto g_0_yyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 381);

    auto g_0_yyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 382);

    auto g_0_yyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 383);

    auto g_0_yyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 384);

    auto g_0_yyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 385);

    auto g_0_yyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 386);

    auto g_0_yyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 387);

    auto g_0_yyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 388);

    auto g_0_yyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 389);

    auto g_0_yyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 390);

    auto g_0_yyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 391);

    auto g_0_yyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 392);

    auto g_0_yyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 393);

    auto g_0_yyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 394);

    auto g_0_yyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 395);

    auto g_0_yyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 396);

    auto g_0_yyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 397);

    auto g_0_yyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 398);

    auto g_0_yzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 399);

    auto g_0_yzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 400);

    auto g_0_yzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 401);

    auto g_0_yzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 402);

    auto g_0_yzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 403);

    auto g_0_yzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 404);

    auto g_0_yzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 405);

    auto g_0_yzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 406);

    auto g_0_yzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 407);

    auto g_0_yzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 408);

    auto g_0_yzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 409);

    auto g_0_yzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 410);

    auto g_0_yzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 411);

    auto g_0_yzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 412);

    auto g_0_yzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 413);

    auto g_0_yzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 414);

    auto g_0_yzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 415);

    auto g_0_yzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 416);

    auto g_0_yzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 417);

    auto g_0_yzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 418);

    auto g_0_yzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 419);

    auto g_0_zzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 420);

    auto g_0_zzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 421);

    auto g_0_zzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 422);

    auto g_0_zzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 423);

    auto g_0_zzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 424);

    auto g_0_zzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 425);

    auto g_0_zzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 426);

    auto g_0_zzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 427);

    auto g_0_zzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 428);

    auto g_0_zzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 429);

    auto g_0_zzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 430);

    auto g_0_zzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 431);

    auto g_0_zzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 432);

    auto g_0_zzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 433);

    auto g_0_zzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 434);

    auto g_0_zzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 435);

    auto g_0_zzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 436);

    auto g_0_zzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 437);

    auto g_0_zzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 438);

    auto g_0_zzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 439);

    auto g_0_zzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 440);

    /// Set up components of auxilary buffer : SHSH

    auto g_0_xxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh);

    auto g_0_xxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 1);

    auto g_0_xxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 2);

    auto g_0_xxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 3);

    auto g_0_xxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 4);

    auto g_0_xxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 5);

    auto g_0_xxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 6);

    auto g_0_xxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 7);

    auto g_0_xxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 8);

    auto g_0_xxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 9);

    auto g_0_xxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 10);

    auto g_0_xxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 11);

    auto g_0_xxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 12);

    auto g_0_xxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 13);

    auto g_0_xxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 14);

    auto g_0_xxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 15);

    auto g_0_xxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 16);

    auto g_0_xxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 17);

    auto g_0_xxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 18);

    auto g_0_xxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 19);

    auto g_0_xxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 20);

    auto g_0_xxxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 21);

    auto g_0_xxxxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 22);

    auto g_0_xxxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 23);

    auto g_0_xxxxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 24);

    auto g_0_xxxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 26);

    auto g_0_xxxxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 27);

    auto g_0_xxxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 30);

    auto g_0_xxxxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 31);

    auto g_0_xxxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 35);

    auto g_0_xxxxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 36);

    auto g_0_xxxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 42);

    auto g_0_xxxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 43);

    auto g_0_xxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 44);

    auto g_0_xxxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 45);

    auto g_0_xxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 46);

    auto g_0_xxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 47);

    auto g_0_xxxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 48);

    auto g_0_xxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 49);

    auto g_0_xxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 50);

    auto g_0_xxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 51);

    auto g_0_xxxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 52);

    auto g_0_xxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 53);

    auto g_0_xxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 54);

    auto g_0_xxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 55);

    auto g_0_xxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 56);

    auto g_0_xxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 58);

    auto g_0_xxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 59);

    auto g_0_xxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 60);

    auto g_0_xxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 61);

    auto g_0_xxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 62);

    auto g_0_xxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 63);

    auto g_0_xxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 64);

    auto g_0_xxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 65);

    auto g_0_xxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 66);

    auto g_0_xxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 67);

    auto g_0_xxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 68);

    auto g_0_xxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 69);

    auto g_0_xxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 70);

    auto g_0_xxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 71);

    auto g_0_xxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 72);

    auto g_0_xxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 73);

    auto g_0_xxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 74);

    auto g_0_xxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 75);

    auto g_0_xxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 76);

    auto g_0_xxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 77);

    auto g_0_xxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 78);

    auto g_0_xxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 79);

    auto g_0_xxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 80);

    auto g_0_xxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 81);

    auto g_0_xxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 82);

    auto g_0_xxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 83);

    auto g_0_xxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 105);

    auto g_0_xxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 106);

    auto g_0_xxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 107);

    auto g_0_xxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 108);

    auto g_0_xxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 109);

    auto g_0_xxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 110);

    auto g_0_xxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 111);

    auto g_0_xxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 112);

    auto g_0_xxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 113);

    auto g_0_xxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 114);

    auto g_0_xxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 115);

    auto g_0_xxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 116);

    auto g_0_xxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 117);

    auto g_0_xxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 118);

    auto g_0_xxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 119);

    auto g_0_xxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 120);

    auto g_0_xxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 121);

    auto g_0_xxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 122);

    auto g_0_xxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 123);

    auto g_0_xxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 124);

    auto g_0_xxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 125);

    auto g_0_xxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 126);

    auto g_0_xxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 127);

    auto g_0_xxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 128);

    auto g_0_xxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 129);

    auto g_0_xxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 130);

    auto g_0_xxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 131);

    auto g_0_xxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 132);

    auto g_0_xxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 133);

    auto g_0_xxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 134);

    auto g_0_xxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 135);

    auto g_0_xxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 136);

    auto g_0_xxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 137);

    auto g_0_xxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 138);

    auto g_0_xxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 139);

    auto g_0_xxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 140);

    auto g_0_xxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 141);

    auto g_0_xxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 142);

    auto g_0_xxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 143);

    auto g_0_xxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 144);

    auto g_0_xxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 145);

    auto g_0_xxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 146);

    auto g_0_xxyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 148);

    auto g_0_xxyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 150);

    auto g_0_xxyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 153);

    auto g_0_xxyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 157);

    auto g_0_xxyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 168);

    auto g_0_xxyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 170);

    auto g_0_xxyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 173);

    auto g_0_xxyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 177);

    auto g_0_xxyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 182);

    auto g_0_xxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 189);

    auto g_0_xxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 190);

    auto g_0_xxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 191);

    auto g_0_xxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 192);

    auto g_0_xxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 193);

    auto g_0_xxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 194);

    auto g_0_xxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 195);

    auto g_0_xxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 196);

    auto g_0_xxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 197);

    auto g_0_xxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 198);

    auto g_0_xxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 199);

    auto g_0_xxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 200);

    auto g_0_xxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 201);

    auto g_0_xxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 202);

    auto g_0_xxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 203);

    auto g_0_xxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 204);

    auto g_0_xxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 205);

    auto g_0_xxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 206);

    auto g_0_xxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 207);

    auto g_0_xxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 208);

    auto g_0_xxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 209);

    auto g_0_xyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 210);

    auto g_0_xyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 211);

    auto g_0_xyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 213);

    auto g_0_xyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 214);

    auto g_0_xyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 216);

    auto g_0_xyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 217);

    auto g_0_xyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 218);

    auto g_0_xyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 220);

    auto g_0_xyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 221);

    auto g_0_xyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 222);

    auto g_0_xyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 223);

    auto g_0_xyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 225);

    auto g_0_xyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 226);

    auto g_0_xyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 227);

    auto g_0_xyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 228);

    auto g_0_xyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 229);

    auto g_0_xyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 230);

    auto g_0_xyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 256);

    auto g_0_xyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 259);

    auto g_0_xyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 260);

    auto g_0_xyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 263);

    auto g_0_xyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 264);

    auto g_0_xyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 265);

    auto g_0_xyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 267);

    auto g_0_xyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 268);

    auto g_0_xyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 269);

    auto g_0_xyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 270);

    auto g_0_xyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 271);

    auto g_0_xyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 272);

    auto g_0_xzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 294);

    auto g_0_xzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 296);

    auto g_0_xzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 298);

    auto g_0_xzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 299);

    auto g_0_xzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 301);

    auto g_0_xzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 302);

    auto g_0_xzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 303);

    auto g_0_xzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 305);

    auto g_0_xzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 306);

    auto g_0_xzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 307);

    auto g_0_xzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 308);

    auto g_0_xzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 309);

    auto g_0_xzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 310);

    auto g_0_xzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 311);

    auto g_0_xzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 312);

    auto g_0_xzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 313);

    auto g_0_xzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 314);

    auto g_0_yyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 315);

    auto g_0_yyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 316);

    auto g_0_yyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 317);

    auto g_0_yyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 318);

    auto g_0_yyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 319);

    auto g_0_yyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 320);

    auto g_0_yyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 321);

    auto g_0_yyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 322);

    auto g_0_yyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 323);

    auto g_0_yyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 324);

    auto g_0_yyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 325);

    auto g_0_yyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 326);

    auto g_0_yyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 327);

    auto g_0_yyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 328);

    auto g_0_yyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 329);

    auto g_0_yyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 330);

    auto g_0_yyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 331);

    auto g_0_yyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 332);

    auto g_0_yyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 333);

    auto g_0_yyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 334);

    auto g_0_yyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 335);

    auto g_0_yyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 337);

    auto g_0_yyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 338);

    auto g_0_yyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 339);

    auto g_0_yyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 340);

    auto g_0_yyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 341);

    auto g_0_yyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 342);

    auto g_0_yyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 343);

    auto g_0_yyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 344);

    auto g_0_yyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 345);

    auto g_0_yyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 346);

    auto g_0_yyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 347);

    auto g_0_yyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 348);

    auto g_0_yyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 349);

    auto g_0_yyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 350);

    auto g_0_yyyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 351);

    auto g_0_yyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 352);

    auto g_0_yyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 353);

    auto g_0_yyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 354);

    auto g_0_yyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 355);

    auto g_0_yyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 356);

    auto g_0_yyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 357);

    auto g_0_yyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 358);

    auto g_0_yyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 359);

    auto g_0_yyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 360);

    auto g_0_yyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 361);

    auto g_0_yyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 362);

    auto g_0_yyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 363);

    auto g_0_yyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 364);

    auto g_0_yyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 365);

    auto g_0_yyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 366);

    auto g_0_yyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 367);

    auto g_0_yyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 368);

    auto g_0_yyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 369);

    auto g_0_yyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 370);

    auto g_0_yyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 371);

    auto g_0_yyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 372);

    auto g_0_yyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 373);

    auto g_0_yyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 374);

    auto g_0_yyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 375);

    auto g_0_yyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 376);

    auto g_0_yyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 377);

    auto g_0_yyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 378);

    auto g_0_yyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 379);

    auto g_0_yyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 380);

    auto g_0_yyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 381);

    auto g_0_yyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 382);

    auto g_0_yyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 383);

    auto g_0_yyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 384);

    auto g_0_yyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 385);

    auto g_0_yyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 386);

    auto g_0_yyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 387);

    auto g_0_yyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 388);

    auto g_0_yyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 389);

    auto g_0_yyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 390);

    auto g_0_yyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 391);

    auto g_0_yyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 392);

    auto g_0_yyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 393);

    auto g_0_yyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 394);

    auto g_0_yyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 395);

    auto g_0_yyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 396);

    auto g_0_yyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 397);

    auto g_0_yyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 398);

    auto g_0_yzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 399);

    auto g_0_yzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 400);

    auto g_0_yzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 401);

    auto g_0_yzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 402);

    auto g_0_yzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 403);

    auto g_0_yzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 404);

    auto g_0_yzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 405);

    auto g_0_yzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 406);

    auto g_0_yzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 407);

    auto g_0_yzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 408);

    auto g_0_yzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 409);

    auto g_0_yzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 410);

    auto g_0_yzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 411);

    auto g_0_yzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 412);

    auto g_0_yzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 413);

    auto g_0_yzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 414);

    auto g_0_yzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 415);

    auto g_0_yzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 416);

    auto g_0_yzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 417);

    auto g_0_yzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 418);

    auto g_0_yzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 419);

    auto g_0_zzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 420);

    auto g_0_zzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 421);

    auto g_0_zzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 422);

    auto g_0_zzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 423);

    auto g_0_zzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 424);

    auto g_0_zzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 425);

    auto g_0_zzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 426);

    auto g_0_zzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 427);

    auto g_0_zzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 428);

    auto g_0_zzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 429);

    auto g_0_zzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 430);

    auto g_0_zzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 431);

    auto g_0_zzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 432);

    auto g_0_zzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 433);

    auto g_0_zzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 434);

    auto g_0_zzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 435);

    auto g_0_zzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 436);

    auto g_0_zzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 437);

    auto g_0_zzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 438);

    auto g_0_zzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 439);

    auto g_0_zzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 440);

    /// Set up 0-21 components of targeted buffer : SISH

    auto g_0_xxxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish);

    auto g_0_xxxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 1);

    auto g_0_xxxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 2);

    auto g_0_xxxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 3);

    auto g_0_xxxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 4);

    auto g_0_xxxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 5);

    auto g_0_xxxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 6);

    auto g_0_xxxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 7);

    auto g_0_xxxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 8);

    auto g_0_xxxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 9);

    auto g_0_xxxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 10);

    auto g_0_xxxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 11);

    auto g_0_xxxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 12);

    auto g_0_xxxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 13);

    auto g_0_xxxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 14);

    auto g_0_xxxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 15);

    auto g_0_xxxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 16);

    auto g_0_xxxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 17);

    auto g_0_xxxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 18);

    auto g_0_xxxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 19);

    auto g_0_xxxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 20);

#pragma omp simd aligned(g_0_xxxx_0_xxxxx_0,       \
                             g_0_xxxx_0_xxxxx_1,   \
                             g_0_xxxx_0_xxxxy_0,   \
                             g_0_xxxx_0_xxxxy_1,   \
                             g_0_xxxx_0_xxxxz_0,   \
                             g_0_xxxx_0_xxxxz_1,   \
                             g_0_xxxx_0_xxxyy_0,   \
                             g_0_xxxx_0_xxxyy_1,   \
                             g_0_xxxx_0_xxxyz_0,   \
                             g_0_xxxx_0_xxxyz_1,   \
                             g_0_xxxx_0_xxxzz_0,   \
                             g_0_xxxx_0_xxxzz_1,   \
                             g_0_xxxx_0_xxyyy_0,   \
                             g_0_xxxx_0_xxyyy_1,   \
                             g_0_xxxx_0_xxyyz_0,   \
                             g_0_xxxx_0_xxyyz_1,   \
                             g_0_xxxx_0_xxyzz_0,   \
                             g_0_xxxx_0_xxyzz_1,   \
                             g_0_xxxx_0_xxzzz_0,   \
                             g_0_xxxx_0_xxzzz_1,   \
                             g_0_xxxx_0_xyyyy_0,   \
                             g_0_xxxx_0_xyyyy_1,   \
                             g_0_xxxx_0_xyyyz_0,   \
                             g_0_xxxx_0_xyyyz_1,   \
                             g_0_xxxx_0_xyyzz_0,   \
                             g_0_xxxx_0_xyyzz_1,   \
                             g_0_xxxx_0_xyzzz_0,   \
                             g_0_xxxx_0_xyzzz_1,   \
                             g_0_xxxx_0_xzzzz_0,   \
                             g_0_xxxx_0_xzzzz_1,   \
                             g_0_xxxx_0_yyyyy_0,   \
                             g_0_xxxx_0_yyyyy_1,   \
                             g_0_xxxx_0_yyyyz_0,   \
                             g_0_xxxx_0_yyyyz_1,   \
                             g_0_xxxx_0_yyyzz_0,   \
                             g_0_xxxx_0_yyyzz_1,   \
                             g_0_xxxx_0_yyzzz_0,   \
                             g_0_xxxx_0_yyzzz_1,   \
                             g_0_xxxx_0_yzzzz_0,   \
                             g_0_xxxx_0_yzzzz_1,   \
                             g_0_xxxx_0_zzzzz_0,   \
                             g_0_xxxx_0_zzzzz_1,   \
                             g_0_xxxxx_0_xxxx_1,   \
                             g_0_xxxxx_0_xxxxx_0,  \
                             g_0_xxxxx_0_xxxxx_1,  \
                             g_0_xxxxx_0_xxxxy_0,  \
                             g_0_xxxxx_0_xxxxy_1,  \
                             g_0_xxxxx_0_xxxxz_0,  \
                             g_0_xxxxx_0_xxxxz_1,  \
                             g_0_xxxxx_0_xxxy_1,   \
                             g_0_xxxxx_0_xxxyy_0,  \
                             g_0_xxxxx_0_xxxyy_1,  \
                             g_0_xxxxx_0_xxxyz_0,  \
                             g_0_xxxxx_0_xxxyz_1,  \
                             g_0_xxxxx_0_xxxz_1,   \
                             g_0_xxxxx_0_xxxzz_0,  \
                             g_0_xxxxx_0_xxxzz_1,  \
                             g_0_xxxxx_0_xxyy_1,   \
                             g_0_xxxxx_0_xxyyy_0,  \
                             g_0_xxxxx_0_xxyyy_1,  \
                             g_0_xxxxx_0_xxyyz_0,  \
                             g_0_xxxxx_0_xxyyz_1,  \
                             g_0_xxxxx_0_xxyz_1,   \
                             g_0_xxxxx_0_xxyzz_0,  \
                             g_0_xxxxx_0_xxyzz_1,  \
                             g_0_xxxxx_0_xxzz_1,   \
                             g_0_xxxxx_0_xxzzz_0,  \
                             g_0_xxxxx_0_xxzzz_1,  \
                             g_0_xxxxx_0_xyyy_1,   \
                             g_0_xxxxx_0_xyyyy_0,  \
                             g_0_xxxxx_0_xyyyy_1,  \
                             g_0_xxxxx_0_xyyyz_0,  \
                             g_0_xxxxx_0_xyyyz_1,  \
                             g_0_xxxxx_0_xyyz_1,   \
                             g_0_xxxxx_0_xyyzz_0,  \
                             g_0_xxxxx_0_xyyzz_1,  \
                             g_0_xxxxx_0_xyzz_1,   \
                             g_0_xxxxx_0_xyzzz_0,  \
                             g_0_xxxxx_0_xyzzz_1,  \
                             g_0_xxxxx_0_xzzz_1,   \
                             g_0_xxxxx_0_xzzzz_0,  \
                             g_0_xxxxx_0_xzzzz_1,  \
                             g_0_xxxxx_0_yyyy_1,   \
                             g_0_xxxxx_0_yyyyy_0,  \
                             g_0_xxxxx_0_yyyyy_1,  \
                             g_0_xxxxx_0_yyyyz_0,  \
                             g_0_xxxxx_0_yyyyz_1,  \
                             g_0_xxxxx_0_yyyz_1,   \
                             g_0_xxxxx_0_yyyzz_0,  \
                             g_0_xxxxx_0_yyyzz_1,  \
                             g_0_xxxxx_0_yyzz_1,   \
                             g_0_xxxxx_0_yyzzz_0,  \
                             g_0_xxxxx_0_yyzzz_1,  \
                             g_0_xxxxx_0_yzzz_1,   \
                             g_0_xxxxx_0_yzzzz_0,  \
                             g_0_xxxxx_0_yzzzz_1,  \
                             g_0_xxxxx_0_zzzz_1,   \
                             g_0_xxxxx_0_zzzzz_0,  \
                             g_0_xxxxx_0_zzzzz_1,  \
                             g_0_xxxxxx_0_xxxxx_0, \
                             g_0_xxxxxx_0_xxxxy_0, \
                             g_0_xxxxxx_0_xxxxz_0, \
                             g_0_xxxxxx_0_xxxyy_0, \
                             g_0_xxxxxx_0_xxxyz_0, \
                             g_0_xxxxxx_0_xxxzz_0, \
                             g_0_xxxxxx_0_xxyyy_0, \
                             g_0_xxxxxx_0_xxyyz_0, \
                             g_0_xxxxxx_0_xxyzz_0, \
                             g_0_xxxxxx_0_xxzzz_0, \
                             g_0_xxxxxx_0_xyyyy_0, \
                             g_0_xxxxxx_0_xyyyz_0, \
                             g_0_xxxxxx_0_xyyzz_0, \
                             g_0_xxxxxx_0_xyzzz_0, \
                             g_0_xxxxxx_0_xzzzz_0, \
                             g_0_xxxxxx_0_yyyyy_0, \
                             g_0_xxxxxx_0_yyyyz_0, \
                             g_0_xxxxxx_0_yyyzz_0, \
                             g_0_xxxxxx_0_yyzzz_0, \
                             g_0_xxxxxx_0_yzzzz_0, \
                             g_0_xxxxxx_0_zzzzz_0, \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xxxxx_0[i] = 5.0 * g_0_xxxx_0_xxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxx_1[i] * fti_ab_0 +
                                  5.0 * g_0_xxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxx_0[i] * pb_x + g_0_xxxxx_0_xxxxx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxy_0[i] = 5.0 * g_0_xxxx_0_xxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxy_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxy_0[i] * pb_x + g_0_xxxxx_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxz_0[i] = 5.0 * g_0_xxxx_0_xxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxz_0[i] * pb_x + g_0_xxxxx_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyy_0[i] = 5.0 * g_0_xxxx_0_xxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyy_0[i] * pb_x + g_0_xxxxx_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyz_0[i] = 5.0 * g_0_xxxx_0_xxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyz_0[i] * pb_x + g_0_xxxxx_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxzz_0[i] = 5.0 * g_0_xxxx_0_xxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzz_0[i] * pb_x + g_0_xxxxx_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyy_0[i] = 5.0 * g_0_xxxx_0_xxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyy_0[i] * pb_x + g_0_xxxxx_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyz_0[i] = 5.0 * g_0_xxxx_0_xxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyz_0[i] * pb_x + g_0_xxxxx_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyzz_0[i] = 5.0 * g_0_xxxx_0_xxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzz_0[i] * pb_x + g_0_xxxxx_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxzzz_0[i] = 5.0 * g_0_xxxx_0_xxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzz_0[i] * pb_x + g_0_xxxxx_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyy_0[i] = 5.0 * g_0_xxxx_0_xyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyy_1[i] * fi_abcd_0 +
                                  g_0_xxxxx_0_xyyyy_0[i] * pb_x + g_0_xxxxx_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyz_0[i] = 5.0 * g_0_xxxx_0_xyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xxxxx_0_xyyyz_0[i] * pb_x + g_0_xxxxx_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyzz_0[i] = 5.0 * g_0_xxxx_0_xyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxx_0_xyyzz_0[i] * pb_x + g_0_xxxxx_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyzzz_0[i] = 5.0 * g_0_xxxx_0_xyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxx_0_xyzzz_0[i] * pb_x + g_0_xxxxx_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xzzzz_0[i] = 5.0 * g_0_xxxx_0_xzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxx_0_xzzzz_0[i] * pb_x + g_0_xxxxx_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyy_0[i] = 5.0 * g_0_xxxx_0_yyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyy_0[i] * pb_x +
                                  g_0_xxxxx_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyz_0[i] = 5.0 * g_0_xxxx_0_yyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyz_0[i] * pb_x +
                                  g_0_xxxxx_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyzz_0[i] = 5.0 * g_0_xxxx_0_yyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyzz_0[i] * pb_x +
                                  g_0_xxxxx_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyzzz_0[i] = 5.0 * g_0_xxxx_0_yyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyzzz_0[i] * pb_x +
                                  g_0_xxxxx_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yzzzz_0[i] = 5.0 * g_0_xxxx_0_yzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzzzz_0[i] * pb_x +
                                  g_0_xxxxx_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zzzzz_0[i] = 5.0 * g_0_xxxx_0_zzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzzzz_0[i] * pb_x +
                                  g_0_xxxxx_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : SISH

    auto g_0_xxxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 21);

    auto g_0_xxxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 22);

    auto g_0_xxxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 23);

    auto g_0_xxxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 24);

    auto g_0_xxxxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 25);

    auto g_0_xxxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 26);

    auto g_0_xxxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 27);

    auto g_0_xxxxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 28);

    auto g_0_xxxxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 29);

    auto g_0_xxxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 30);

    auto g_0_xxxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 31);

    auto g_0_xxxxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 32);

    auto g_0_xxxxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 33);

    auto g_0_xxxxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 34);

    auto g_0_xxxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 35);

    auto g_0_xxxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 36);

    auto g_0_xxxxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 37);

    auto g_0_xxxxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 38);

    auto g_0_xxxxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 39);

    auto g_0_xxxxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 40);

    auto g_0_xxxxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 41);

#pragma omp simd aligned(g_0_xxxxx_0_xxxx_1,       \
                             g_0_xxxxx_0_xxxxx_0,  \
                             g_0_xxxxx_0_xxxxx_1,  \
                             g_0_xxxxx_0_xxxxy_0,  \
                             g_0_xxxxx_0_xxxxy_1,  \
                             g_0_xxxxx_0_xxxxz_0,  \
                             g_0_xxxxx_0_xxxxz_1,  \
                             g_0_xxxxx_0_xxxy_1,   \
                             g_0_xxxxx_0_xxxyy_0,  \
                             g_0_xxxxx_0_xxxyy_1,  \
                             g_0_xxxxx_0_xxxyz_0,  \
                             g_0_xxxxx_0_xxxyz_1,  \
                             g_0_xxxxx_0_xxxz_1,   \
                             g_0_xxxxx_0_xxxzz_0,  \
                             g_0_xxxxx_0_xxxzz_1,  \
                             g_0_xxxxx_0_xxyy_1,   \
                             g_0_xxxxx_0_xxyyy_0,  \
                             g_0_xxxxx_0_xxyyy_1,  \
                             g_0_xxxxx_0_xxyyz_0,  \
                             g_0_xxxxx_0_xxyyz_1,  \
                             g_0_xxxxx_0_xxyz_1,   \
                             g_0_xxxxx_0_xxyzz_0,  \
                             g_0_xxxxx_0_xxyzz_1,  \
                             g_0_xxxxx_0_xxzz_1,   \
                             g_0_xxxxx_0_xxzzz_0,  \
                             g_0_xxxxx_0_xxzzz_1,  \
                             g_0_xxxxx_0_xyyy_1,   \
                             g_0_xxxxx_0_xyyyy_0,  \
                             g_0_xxxxx_0_xyyyy_1,  \
                             g_0_xxxxx_0_xyyyz_0,  \
                             g_0_xxxxx_0_xyyyz_1,  \
                             g_0_xxxxx_0_xyyz_1,   \
                             g_0_xxxxx_0_xyyzz_0,  \
                             g_0_xxxxx_0_xyyzz_1,  \
                             g_0_xxxxx_0_xyzz_1,   \
                             g_0_xxxxx_0_xyzzz_0,  \
                             g_0_xxxxx_0_xyzzz_1,  \
                             g_0_xxxxx_0_xzzz_1,   \
                             g_0_xxxxx_0_xzzzz_0,  \
                             g_0_xxxxx_0_xzzzz_1,  \
                             g_0_xxxxx_0_yyyy_1,   \
                             g_0_xxxxx_0_yyyyy_0,  \
                             g_0_xxxxx_0_yyyyy_1,  \
                             g_0_xxxxx_0_yyyyz_0,  \
                             g_0_xxxxx_0_yyyyz_1,  \
                             g_0_xxxxx_0_yyyz_1,   \
                             g_0_xxxxx_0_yyyzz_0,  \
                             g_0_xxxxx_0_yyyzz_1,  \
                             g_0_xxxxx_0_yyzz_1,   \
                             g_0_xxxxx_0_yyzzz_0,  \
                             g_0_xxxxx_0_yyzzz_1,  \
                             g_0_xxxxx_0_yzzz_1,   \
                             g_0_xxxxx_0_yzzzz_0,  \
                             g_0_xxxxx_0_yzzzz_1,  \
                             g_0_xxxxx_0_zzzz_1,   \
                             g_0_xxxxx_0_zzzzz_0,  \
                             g_0_xxxxx_0_zzzzz_1,  \
                             g_0_xxxxxy_0_xxxxx_0, \
                             g_0_xxxxxy_0_xxxxy_0, \
                             g_0_xxxxxy_0_xxxxz_0, \
                             g_0_xxxxxy_0_xxxyy_0, \
                             g_0_xxxxxy_0_xxxyz_0, \
                             g_0_xxxxxy_0_xxxzz_0, \
                             g_0_xxxxxy_0_xxyyy_0, \
                             g_0_xxxxxy_0_xxyyz_0, \
                             g_0_xxxxxy_0_xxyzz_0, \
                             g_0_xxxxxy_0_xxzzz_0, \
                             g_0_xxxxxy_0_xyyyy_0, \
                             g_0_xxxxxy_0_xyyyz_0, \
                             g_0_xxxxxy_0_xyyzz_0, \
                             g_0_xxxxxy_0_xyzzz_0, \
                             g_0_xxxxxy_0_xzzzz_0, \
                             g_0_xxxxxy_0_yyyyy_0, \
                             g_0_xxxxxy_0_yyyyz_0, \
                             g_0_xxxxxy_0_yyyzz_0, \
                             g_0_xxxxxy_0_yyzzz_0, \
                             g_0_xxxxxy_0_yzzzz_0, \
                             g_0_xxxxxy_0_zzzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xxxxx_0[i] = g_0_xxxxx_0_xxxxx_0[i] * pb_y + g_0_xxxxx_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxy_0[i] = g_0_xxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxy_0[i] * pb_y + g_0_xxxxx_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxz_0[i] = g_0_xxxxx_0_xxxxz_0[i] * pb_y + g_0_xxxxx_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyy_0[i] = 2.0 * g_0_xxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyy_0[i] * pb_y + g_0_xxxxx_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyz_0[i] = g_0_xxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyz_0[i] * pb_y + g_0_xxxxx_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxzz_0[i] = g_0_xxxxx_0_xxxzz_0[i] * pb_y + g_0_xxxxx_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyy_0[i] = 3.0 * g_0_xxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyy_0[i] * pb_y + g_0_xxxxx_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyz_0[i] = 2.0 * g_0_xxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyz_0[i] * pb_y + g_0_xxxxx_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyzz_0[i] = g_0_xxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzz_0[i] * pb_y + g_0_xxxxx_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxzzz_0[i] = g_0_xxxxx_0_xxzzz_0[i] * pb_y + g_0_xxxxx_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyy_0[i] = 4.0 * g_0_xxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyy_0[i] * pb_y + g_0_xxxxx_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyz_0[i] = 3.0 * g_0_xxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyz_0[i] * pb_y + g_0_xxxxx_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyzz_0[i] = 2.0 * g_0_xxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzz_0[i] * pb_y + g_0_xxxxx_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyzzz_0[i] = g_0_xxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzz_0[i] * pb_y + g_0_xxxxx_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xzzzz_0[i] = g_0_xxxxx_0_xzzzz_0[i] * pb_y + g_0_xxxxx_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyy_0[i] = 5.0 * g_0_xxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyy_0[i] * pb_y + g_0_xxxxx_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyz_0[i] = 4.0 * g_0_xxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyz_0[i] * pb_y + g_0_xxxxx_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyzz_0[i] = 3.0 * g_0_xxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzz_0[i] * pb_y + g_0_xxxxx_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyzzz_0[i] = 2.0 * g_0_xxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzz_0[i] * pb_y + g_0_xxxxx_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yzzzz_0[i] = g_0_xxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzz_0[i] * pb_y + g_0_xxxxx_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zzzzz_0[i] = g_0_xxxxx_0_zzzzz_0[i] * pb_y + g_0_xxxxx_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : SISH

    auto g_0_xxxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 42);

    auto g_0_xxxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 43);

    auto g_0_xxxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 44);

    auto g_0_xxxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 45);

    auto g_0_xxxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 46);

    auto g_0_xxxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 47);

    auto g_0_xxxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 48);

    auto g_0_xxxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 49);

    auto g_0_xxxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 50);

    auto g_0_xxxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 51);

    auto g_0_xxxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 52);

    auto g_0_xxxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 53);

    auto g_0_xxxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 54);

    auto g_0_xxxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 55);

    auto g_0_xxxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 56);

    auto g_0_xxxxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 57);

    auto g_0_xxxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 58);

    auto g_0_xxxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 59);

    auto g_0_xxxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 60);

    auto g_0_xxxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 61);

    auto g_0_xxxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 62);

#pragma omp simd aligned(g_0_xxxxx_0_xxxx_1,       \
                             g_0_xxxxx_0_xxxxx_0,  \
                             g_0_xxxxx_0_xxxxx_1,  \
                             g_0_xxxxx_0_xxxxy_0,  \
                             g_0_xxxxx_0_xxxxy_1,  \
                             g_0_xxxxx_0_xxxxz_0,  \
                             g_0_xxxxx_0_xxxxz_1,  \
                             g_0_xxxxx_0_xxxy_1,   \
                             g_0_xxxxx_0_xxxyy_0,  \
                             g_0_xxxxx_0_xxxyy_1,  \
                             g_0_xxxxx_0_xxxyz_0,  \
                             g_0_xxxxx_0_xxxyz_1,  \
                             g_0_xxxxx_0_xxxz_1,   \
                             g_0_xxxxx_0_xxxzz_0,  \
                             g_0_xxxxx_0_xxxzz_1,  \
                             g_0_xxxxx_0_xxyy_1,   \
                             g_0_xxxxx_0_xxyyy_0,  \
                             g_0_xxxxx_0_xxyyy_1,  \
                             g_0_xxxxx_0_xxyyz_0,  \
                             g_0_xxxxx_0_xxyyz_1,  \
                             g_0_xxxxx_0_xxyz_1,   \
                             g_0_xxxxx_0_xxyzz_0,  \
                             g_0_xxxxx_0_xxyzz_1,  \
                             g_0_xxxxx_0_xxzz_1,   \
                             g_0_xxxxx_0_xxzzz_0,  \
                             g_0_xxxxx_0_xxzzz_1,  \
                             g_0_xxxxx_0_xyyy_1,   \
                             g_0_xxxxx_0_xyyyy_0,  \
                             g_0_xxxxx_0_xyyyy_1,  \
                             g_0_xxxxx_0_xyyyz_0,  \
                             g_0_xxxxx_0_xyyyz_1,  \
                             g_0_xxxxx_0_xyyz_1,   \
                             g_0_xxxxx_0_xyyzz_0,  \
                             g_0_xxxxx_0_xyyzz_1,  \
                             g_0_xxxxx_0_xyzz_1,   \
                             g_0_xxxxx_0_xyzzz_0,  \
                             g_0_xxxxx_0_xyzzz_1,  \
                             g_0_xxxxx_0_xzzz_1,   \
                             g_0_xxxxx_0_xzzzz_0,  \
                             g_0_xxxxx_0_xzzzz_1,  \
                             g_0_xxxxx_0_yyyy_1,   \
                             g_0_xxxxx_0_yyyyy_0,  \
                             g_0_xxxxx_0_yyyyy_1,  \
                             g_0_xxxxx_0_yyyyz_0,  \
                             g_0_xxxxx_0_yyyyz_1,  \
                             g_0_xxxxx_0_yyyz_1,   \
                             g_0_xxxxx_0_yyyzz_0,  \
                             g_0_xxxxx_0_yyyzz_1,  \
                             g_0_xxxxx_0_yyzz_1,   \
                             g_0_xxxxx_0_yyzzz_0,  \
                             g_0_xxxxx_0_yyzzz_1,  \
                             g_0_xxxxx_0_yzzz_1,   \
                             g_0_xxxxx_0_yzzzz_0,  \
                             g_0_xxxxx_0_yzzzz_1,  \
                             g_0_xxxxx_0_zzzz_1,   \
                             g_0_xxxxx_0_zzzzz_0,  \
                             g_0_xxxxx_0_zzzzz_1,  \
                             g_0_xxxxxz_0_xxxxx_0, \
                             g_0_xxxxxz_0_xxxxy_0, \
                             g_0_xxxxxz_0_xxxxz_0, \
                             g_0_xxxxxz_0_xxxyy_0, \
                             g_0_xxxxxz_0_xxxyz_0, \
                             g_0_xxxxxz_0_xxxzz_0, \
                             g_0_xxxxxz_0_xxyyy_0, \
                             g_0_xxxxxz_0_xxyyz_0, \
                             g_0_xxxxxz_0_xxyzz_0, \
                             g_0_xxxxxz_0_xxzzz_0, \
                             g_0_xxxxxz_0_xyyyy_0, \
                             g_0_xxxxxz_0_xyyyz_0, \
                             g_0_xxxxxz_0_xyyzz_0, \
                             g_0_xxxxxz_0_xyzzz_0, \
                             g_0_xxxxxz_0_xzzzz_0, \
                             g_0_xxxxxz_0_yyyyy_0, \
                             g_0_xxxxxz_0_yyyyz_0, \
                             g_0_xxxxxz_0_yyyzz_0, \
                             g_0_xxxxxz_0_yyzzz_0, \
                             g_0_xxxxxz_0_yzzzz_0, \
                             g_0_xxxxxz_0_zzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xxxxx_0[i] = g_0_xxxxx_0_xxxxx_0[i] * pb_z + g_0_xxxxx_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxy_0[i] = g_0_xxxxx_0_xxxxy_0[i] * pb_z + g_0_xxxxx_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxz_0[i] = g_0_xxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxz_0[i] * pb_z + g_0_xxxxx_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyy_0[i] = g_0_xxxxx_0_xxxyy_0[i] * pb_z + g_0_xxxxx_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyz_0[i] = g_0_xxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyz_0[i] * pb_z + g_0_xxxxx_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxzz_0[i] = 2.0 * g_0_xxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzz_0[i] * pb_z + g_0_xxxxx_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyy_0[i] = g_0_xxxxx_0_xxyyy_0[i] * pb_z + g_0_xxxxx_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyz_0[i] = g_0_xxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyz_0[i] * pb_z + g_0_xxxxx_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyzz_0[i] = 2.0 * g_0_xxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzz_0[i] * pb_z + g_0_xxxxx_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxzzz_0[i] = 3.0 * g_0_xxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzz_0[i] * pb_z + g_0_xxxxx_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyy_0[i] = g_0_xxxxx_0_xyyyy_0[i] * pb_z + g_0_xxxxx_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyz_0[i] = g_0_xxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyz_0[i] * pb_z + g_0_xxxxx_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyzz_0[i] = 2.0 * g_0_xxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzz_0[i] * pb_z + g_0_xxxxx_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyzzz_0[i] = 3.0 * g_0_xxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzz_0[i] * pb_z + g_0_xxxxx_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xzzzz_0[i] = 4.0 * g_0_xxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzz_0[i] * pb_z + g_0_xxxxx_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyy_0[i] = g_0_xxxxx_0_yyyyy_0[i] * pb_z + g_0_xxxxx_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyz_0[i] = g_0_xxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyz_0[i] * pb_z + g_0_xxxxx_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyzz_0[i] = 2.0 * g_0_xxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzz_0[i] * pb_z + g_0_xxxxx_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyzzz_0[i] = 3.0 * g_0_xxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzz_0[i] * pb_z + g_0_xxxxx_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yzzzz_0[i] = 4.0 * g_0_xxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzz_0[i] * pb_z + g_0_xxxxx_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zzzzz_0[i] = 5.0 * g_0_xxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_zzzzz_0[i] * pb_z + g_0_xxxxx_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 63-84 components of targeted buffer : SISH

    auto g_0_xxxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 63);

    auto g_0_xxxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 64);

    auto g_0_xxxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 65);

    auto g_0_xxxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 66);

    auto g_0_xxxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 67);

    auto g_0_xxxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 68);

    auto g_0_xxxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 69);

    auto g_0_xxxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 70);

    auto g_0_xxxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 71);

    auto g_0_xxxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 72);

    auto g_0_xxxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 73);

    auto g_0_xxxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 74);

    auto g_0_xxxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 75);

    auto g_0_xxxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 76);

    auto g_0_xxxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 77);

    auto g_0_xxxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 78);

    auto g_0_xxxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 79);

    auto g_0_xxxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 80);

    auto g_0_xxxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 81);

    auto g_0_xxxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 82);

    auto g_0_xxxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 83);

#pragma omp simd aligned(g_0_xxxx_0_xxxxx_0,       \
                             g_0_xxxx_0_xxxxx_1,   \
                             g_0_xxxx_0_xxxxz_0,   \
                             g_0_xxxx_0_xxxxz_1,   \
                             g_0_xxxx_0_xxxzz_0,   \
                             g_0_xxxx_0_xxxzz_1,   \
                             g_0_xxxx_0_xxzzz_0,   \
                             g_0_xxxx_0_xxzzz_1,   \
                             g_0_xxxx_0_xzzzz_0,   \
                             g_0_xxxx_0_xzzzz_1,   \
                             g_0_xxxxy_0_xxxxx_0,  \
                             g_0_xxxxy_0_xxxxx_1,  \
                             g_0_xxxxy_0_xxxxz_0,  \
                             g_0_xxxxy_0_xxxxz_1,  \
                             g_0_xxxxy_0_xxxzz_0,  \
                             g_0_xxxxy_0_xxxzz_1,  \
                             g_0_xxxxy_0_xxzzz_0,  \
                             g_0_xxxxy_0_xxzzz_1,  \
                             g_0_xxxxy_0_xzzzz_0,  \
                             g_0_xxxxy_0_xzzzz_1,  \
                             g_0_xxxxyy_0_xxxxx_0, \
                             g_0_xxxxyy_0_xxxxy_0, \
                             g_0_xxxxyy_0_xxxxz_0, \
                             g_0_xxxxyy_0_xxxyy_0, \
                             g_0_xxxxyy_0_xxxyz_0, \
                             g_0_xxxxyy_0_xxxzz_0, \
                             g_0_xxxxyy_0_xxyyy_0, \
                             g_0_xxxxyy_0_xxyyz_0, \
                             g_0_xxxxyy_0_xxyzz_0, \
                             g_0_xxxxyy_0_xxzzz_0, \
                             g_0_xxxxyy_0_xyyyy_0, \
                             g_0_xxxxyy_0_xyyyz_0, \
                             g_0_xxxxyy_0_xyyzz_0, \
                             g_0_xxxxyy_0_xyzzz_0, \
                             g_0_xxxxyy_0_xzzzz_0, \
                             g_0_xxxxyy_0_yyyyy_0, \
                             g_0_xxxxyy_0_yyyyz_0, \
                             g_0_xxxxyy_0_yyyzz_0, \
                             g_0_xxxxyy_0_yyzzz_0, \
                             g_0_xxxxyy_0_yzzzz_0, \
                             g_0_xxxxyy_0_zzzzz_0, \
                             g_0_xxxyy_0_xxxxy_0,  \
                             g_0_xxxyy_0_xxxxy_1,  \
                             g_0_xxxyy_0_xxxy_1,   \
                             g_0_xxxyy_0_xxxyy_0,  \
                             g_0_xxxyy_0_xxxyy_1,  \
                             g_0_xxxyy_0_xxxyz_0,  \
                             g_0_xxxyy_0_xxxyz_1,  \
                             g_0_xxxyy_0_xxyy_1,   \
                             g_0_xxxyy_0_xxyyy_0,  \
                             g_0_xxxyy_0_xxyyy_1,  \
                             g_0_xxxyy_0_xxyyz_0,  \
                             g_0_xxxyy_0_xxyyz_1,  \
                             g_0_xxxyy_0_xxyz_1,   \
                             g_0_xxxyy_0_xxyzz_0,  \
                             g_0_xxxyy_0_xxyzz_1,  \
                             g_0_xxxyy_0_xyyy_1,   \
                             g_0_xxxyy_0_xyyyy_0,  \
                             g_0_xxxyy_0_xyyyy_1,  \
                             g_0_xxxyy_0_xyyyz_0,  \
                             g_0_xxxyy_0_xyyyz_1,  \
                             g_0_xxxyy_0_xyyz_1,   \
                             g_0_xxxyy_0_xyyzz_0,  \
                             g_0_xxxyy_0_xyyzz_1,  \
                             g_0_xxxyy_0_xyzz_1,   \
                             g_0_xxxyy_0_xyzzz_0,  \
                             g_0_xxxyy_0_xyzzz_1,  \
                             g_0_xxxyy_0_yyyy_1,   \
                             g_0_xxxyy_0_yyyyy_0,  \
                             g_0_xxxyy_0_yyyyy_1,  \
                             g_0_xxxyy_0_yyyyz_0,  \
                             g_0_xxxyy_0_yyyyz_1,  \
                             g_0_xxxyy_0_yyyz_1,   \
                             g_0_xxxyy_0_yyyzz_0,  \
                             g_0_xxxyy_0_yyyzz_1,  \
                             g_0_xxxyy_0_yyzz_1,   \
                             g_0_xxxyy_0_yyzzz_0,  \
                             g_0_xxxyy_0_yyzzz_1,  \
                             g_0_xxxyy_0_yzzz_1,   \
                             g_0_xxxyy_0_yzzzz_0,  \
                             g_0_xxxyy_0_yzzzz_1,  \
                             g_0_xxxyy_0_zzzzz_0,  \
                             g_0_xxxyy_0_zzzzz_1,  \
                             g_0_xxyy_0_xxxxy_0,   \
                             g_0_xxyy_0_xxxxy_1,   \
                             g_0_xxyy_0_xxxyy_0,   \
                             g_0_xxyy_0_xxxyy_1,   \
                             g_0_xxyy_0_xxxyz_0,   \
                             g_0_xxyy_0_xxxyz_1,   \
                             g_0_xxyy_0_xxyyy_0,   \
                             g_0_xxyy_0_xxyyy_1,   \
                             g_0_xxyy_0_xxyyz_0,   \
                             g_0_xxyy_0_xxyyz_1,   \
                             g_0_xxyy_0_xxyzz_0,   \
                             g_0_xxyy_0_xxyzz_1,   \
                             g_0_xxyy_0_xyyyy_0,   \
                             g_0_xxyy_0_xyyyy_1,   \
                             g_0_xxyy_0_xyyyz_0,   \
                             g_0_xxyy_0_xyyyz_1,   \
                             g_0_xxyy_0_xyyzz_0,   \
                             g_0_xxyy_0_xyyzz_1,   \
                             g_0_xxyy_0_xyzzz_0,   \
                             g_0_xxyy_0_xyzzz_1,   \
                             g_0_xxyy_0_yyyyy_0,   \
                             g_0_xxyy_0_yyyyy_1,   \
                             g_0_xxyy_0_yyyyz_0,   \
                             g_0_xxyy_0_yyyyz_1,   \
                             g_0_xxyy_0_yyyzz_0,   \
                             g_0_xxyy_0_yyyzz_1,   \
                             g_0_xxyy_0_yyzzz_0,   \
                             g_0_xxyy_0_yyzzz_1,   \
                             g_0_xxyy_0_yzzzz_0,   \
                             g_0_xxyy_0_yzzzz_1,   \
                             g_0_xxyy_0_zzzzz_0,   \
                             g_0_xxyy_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xxxxx_0[i] =
            g_0_xxxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxx_0[i] * pb_y + g_0_xxxxy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxy_0[i] = 3.0 * g_0_xxyy_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxy_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxy_0[i] * pb_x + g_0_xxxyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxz_0[i] =
            g_0_xxxx_0_xxxxz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxz_0[i] * pb_y + g_0_xxxxy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxyy_0[i] = 3.0 * g_0_xxyy_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyy_0[i] * pb_x + g_0_xxxyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyz_0[i] = 3.0 * g_0_xxyy_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyz_0[i] * pb_x + g_0_xxxyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxzz_0[i] =
            g_0_xxxx_0_xxxzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxzz_0[i] * pb_y + g_0_xxxxy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxyyy_0[i] = 3.0 * g_0_xxyy_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyy_0[i] * pb_x + g_0_xxxyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyz_0[i] = 3.0 * g_0_xxyy_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyz_0[i] * pb_x + g_0_xxxyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyzz_0[i] = 3.0 * g_0_xxyy_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzz_0[i] * pb_x + g_0_xxxyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxzzz_0[i] =
            g_0_xxxx_0_xxzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxzzz_0[i] * pb_y + g_0_xxxxy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xyyyy_0[i] = 3.0 * g_0_xxyy_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyy_1[i] * fi_abcd_0 +
                                  g_0_xxxyy_0_xyyyy_0[i] * pb_x + g_0_xxxyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyz_0[i] = 3.0 * g_0_xxyy_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xxxyy_0_xyyyz_0[i] * pb_x + g_0_xxxyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyzz_0[i] = 3.0 * g_0_xxyy_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xxxyy_0_xyyzz_0[i] * pb_x + g_0_xxxyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyzzz_0[i] = 3.0 * g_0_xxyy_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxyy_0_xyzzz_0[i] * pb_x + g_0_xxxyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xzzzz_0[i] =
            g_0_xxxx_0_xzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xzzzz_0[i] * pb_y + g_0_xxxxy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yyyyy_0[i] = 3.0 * g_0_xxyy_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyy_0[i] * pb_x +
                                  g_0_xxxyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyz_0[i] = 3.0 * g_0_xxyy_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyz_0[i] * pb_x +
                                  g_0_xxxyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyzz_0[i] = 3.0 * g_0_xxyy_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyzz_0[i] * pb_x +
                                  g_0_xxxyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyzzz_0[i] = 3.0 * g_0_xxyy_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyzzz_0[i] * pb_x +
                                  g_0_xxxyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yzzzz_0[i] = 3.0 * g_0_xxyy_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzzzz_0[i] * pb_x +
                                  g_0_xxxyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zzzzz_0[i] = 3.0 * g_0_xxyy_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_zzzzz_0[i] * pb_x +
                                  g_0_xxxyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-105 components of targeted buffer : SISH

    auto g_0_xxxxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 84);

    auto g_0_xxxxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 85);

    auto g_0_xxxxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 86);

    auto g_0_xxxxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 87);

    auto g_0_xxxxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 88);

    auto g_0_xxxxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 89);

    auto g_0_xxxxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 90);

    auto g_0_xxxxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 91);

    auto g_0_xxxxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 92);

    auto g_0_xxxxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 93);

    auto g_0_xxxxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 94);

    auto g_0_xxxxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 95);

    auto g_0_xxxxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 96);

    auto g_0_xxxxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 97);

    auto g_0_xxxxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 98);

    auto g_0_xxxxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 99);

    auto g_0_xxxxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 100);

    auto g_0_xxxxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 101);

    auto g_0_xxxxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 102);

    auto g_0_xxxxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 103);

    auto g_0_xxxxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 104);

#pragma omp simd aligned(g_0_xxxxy_0_xxxxy_0,      \
                             g_0_xxxxy_0_xxxxy_1,  \
                             g_0_xxxxy_0_xxxyy_0,  \
                             g_0_xxxxy_0_xxxyy_1,  \
                             g_0_xxxxy_0_xxyyy_0,  \
                             g_0_xxxxy_0_xxyyy_1,  \
                             g_0_xxxxy_0_xyyyy_0,  \
                             g_0_xxxxy_0_xyyyy_1,  \
                             g_0_xxxxy_0_yyyyy_0,  \
                             g_0_xxxxy_0_yyyyy_1,  \
                             g_0_xxxxyz_0_xxxxx_0, \
                             g_0_xxxxyz_0_xxxxy_0, \
                             g_0_xxxxyz_0_xxxxz_0, \
                             g_0_xxxxyz_0_xxxyy_0, \
                             g_0_xxxxyz_0_xxxyz_0, \
                             g_0_xxxxyz_0_xxxzz_0, \
                             g_0_xxxxyz_0_xxyyy_0, \
                             g_0_xxxxyz_0_xxyyz_0, \
                             g_0_xxxxyz_0_xxyzz_0, \
                             g_0_xxxxyz_0_xxzzz_0, \
                             g_0_xxxxyz_0_xyyyy_0, \
                             g_0_xxxxyz_0_xyyyz_0, \
                             g_0_xxxxyz_0_xyyzz_0, \
                             g_0_xxxxyz_0_xyzzz_0, \
                             g_0_xxxxyz_0_xzzzz_0, \
                             g_0_xxxxyz_0_yyyyy_0, \
                             g_0_xxxxyz_0_yyyyz_0, \
                             g_0_xxxxyz_0_yyyzz_0, \
                             g_0_xxxxyz_0_yyzzz_0, \
                             g_0_xxxxyz_0_yzzzz_0, \
                             g_0_xxxxyz_0_zzzzz_0, \
                             g_0_xxxxz_0_xxxxx_0,  \
                             g_0_xxxxz_0_xxxxx_1,  \
                             g_0_xxxxz_0_xxxxz_0,  \
                             g_0_xxxxz_0_xxxxz_1,  \
                             g_0_xxxxz_0_xxxyz_0,  \
                             g_0_xxxxz_0_xxxyz_1,  \
                             g_0_xxxxz_0_xxxz_1,   \
                             g_0_xxxxz_0_xxxzz_0,  \
                             g_0_xxxxz_0_xxxzz_1,  \
                             g_0_xxxxz_0_xxyyz_0,  \
                             g_0_xxxxz_0_xxyyz_1,  \
                             g_0_xxxxz_0_xxyz_1,   \
                             g_0_xxxxz_0_xxyzz_0,  \
                             g_0_xxxxz_0_xxyzz_1,  \
                             g_0_xxxxz_0_xxzz_1,   \
                             g_0_xxxxz_0_xxzzz_0,  \
                             g_0_xxxxz_0_xxzzz_1,  \
                             g_0_xxxxz_0_xyyyz_0,  \
                             g_0_xxxxz_0_xyyyz_1,  \
                             g_0_xxxxz_0_xyyz_1,   \
                             g_0_xxxxz_0_xyyzz_0,  \
                             g_0_xxxxz_0_xyyzz_1,  \
                             g_0_xxxxz_0_xyzz_1,   \
                             g_0_xxxxz_0_xyzzz_0,  \
                             g_0_xxxxz_0_xyzzz_1,  \
                             g_0_xxxxz_0_xzzz_1,   \
                             g_0_xxxxz_0_xzzzz_0,  \
                             g_0_xxxxz_0_xzzzz_1,  \
                             g_0_xxxxz_0_yyyyz_0,  \
                             g_0_xxxxz_0_yyyyz_1,  \
                             g_0_xxxxz_0_yyyz_1,   \
                             g_0_xxxxz_0_yyyzz_0,  \
                             g_0_xxxxz_0_yyyzz_1,  \
                             g_0_xxxxz_0_yyzz_1,   \
                             g_0_xxxxz_0_yyzzz_0,  \
                             g_0_xxxxz_0_yyzzz_1,  \
                             g_0_xxxxz_0_yzzz_1,   \
                             g_0_xxxxz_0_yzzzz_0,  \
                             g_0_xxxxz_0_yzzzz_1,  \
                             g_0_xxxxz_0_zzzz_1,   \
                             g_0_xxxxz_0_zzzzz_0,  \
                             g_0_xxxxz_0_zzzzz_1,  \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xxxxx_0[i] = g_0_xxxxz_0_xxxxx_0[i] * pb_y + g_0_xxxxz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxy_0[i] = g_0_xxxxy_0_xxxxy_0[i] * pb_z + g_0_xxxxy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxz_0[i] = g_0_xxxxz_0_xxxxz_0[i] * pb_y + g_0_xxxxz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyy_0[i] = g_0_xxxxy_0_xxxyy_0[i] * pb_z + g_0_xxxxy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxyz_0[i] = g_0_xxxxz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyz_0[i] * pb_y + g_0_xxxxz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxzz_0[i] = g_0_xxxxz_0_xxxzz_0[i] * pb_y + g_0_xxxxz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyy_0[i] = g_0_xxxxy_0_xxyyy_0[i] * pb_z + g_0_xxxxy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxyyz_0[i] = 2.0 * g_0_xxxxz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyz_0[i] * pb_y + g_0_xxxxz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyzz_0[i] = g_0_xxxxz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyzz_0[i] * pb_y + g_0_xxxxz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxzzz_0[i] = g_0_xxxxz_0_xxzzz_0[i] * pb_y + g_0_xxxxz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyy_0[i] = g_0_xxxxy_0_xyyyy_0[i] * pb_z + g_0_xxxxy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xyyyz_0[i] = 3.0 * g_0_xxxxz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyz_0[i] * pb_y + g_0_xxxxz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyzz_0[i] = 2.0 * g_0_xxxxz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyzz_0[i] * pb_y + g_0_xxxxz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyzzz_0[i] = g_0_xxxxz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyzzz_0[i] * pb_y + g_0_xxxxz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xzzzz_0[i] = g_0_xxxxz_0_xzzzz_0[i] * pb_y + g_0_xxxxz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyy_0[i] = g_0_xxxxy_0_yyyyy_0[i] * pb_z + g_0_xxxxy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yyyyz_0[i] = 4.0 * g_0_xxxxz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyz_0[i] * pb_y + g_0_xxxxz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyzz_0[i] = 3.0 * g_0_xxxxz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyzz_0[i] * pb_y + g_0_xxxxz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyzzz_0[i] = 2.0 * g_0_xxxxz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyzzz_0[i] * pb_y + g_0_xxxxz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yzzzz_0[i] = g_0_xxxxz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yzzzz_0[i] * pb_y + g_0_xxxxz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zzzzz_0[i] = g_0_xxxxz_0_zzzzz_0[i] * pb_y + g_0_xxxxz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 105-126 components of targeted buffer : SISH

    auto g_0_xxxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 105);

    auto g_0_xxxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 106);

    auto g_0_xxxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 107);

    auto g_0_xxxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 108);

    auto g_0_xxxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 109);

    auto g_0_xxxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 110);

    auto g_0_xxxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 111);

    auto g_0_xxxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 112);

    auto g_0_xxxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 113);

    auto g_0_xxxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 114);

    auto g_0_xxxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 115);

    auto g_0_xxxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 116);

    auto g_0_xxxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 117);

    auto g_0_xxxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 118);

    auto g_0_xxxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 119);

    auto g_0_xxxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 120);

    auto g_0_xxxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 121);

    auto g_0_xxxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 122);

    auto g_0_xxxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 123);

    auto g_0_xxxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 124);

    auto g_0_xxxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 125);

#pragma omp simd aligned(g_0_xxxx_0_xxxxx_0,       \
                             g_0_xxxx_0_xxxxx_1,   \
                             g_0_xxxx_0_xxxxy_0,   \
                             g_0_xxxx_0_xxxxy_1,   \
                             g_0_xxxx_0_xxxyy_0,   \
                             g_0_xxxx_0_xxxyy_1,   \
                             g_0_xxxx_0_xxyyy_0,   \
                             g_0_xxxx_0_xxyyy_1,   \
                             g_0_xxxx_0_xyyyy_0,   \
                             g_0_xxxx_0_xyyyy_1,   \
                             g_0_xxxxz_0_xxxxx_0,  \
                             g_0_xxxxz_0_xxxxx_1,  \
                             g_0_xxxxz_0_xxxxy_0,  \
                             g_0_xxxxz_0_xxxxy_1,  \
                             g_0_xxxxz_0_xxxyy_0,  \
                             g_0_xxxxz_0_xxxyy_1,  \
                             g_0_xxxxz_0_xxyyy_0,  \
                             g_0_xxxxz_0_xxyyy_1,  \
                             g_0_xxxxz_0_xyyyy_0,  \
                             g_0_xxxxz_0_xyyyy_1,  \
                             g_0_xxxxzz_0_xxxxx_0, \
                             g_0_xxxxzz_0_xxxxy_0, \
                             g_0_xxxxzz_0_xxxxz_0, \
                             g_0_xxxxzz_0_xxxyy_0, \
                             g_0_xxxxzz_0_xxxyz_0, \
                             g_0_xxxxzz_0_xxxzz_0, \
                             g_0_xxxxzz_0_xxyyy_0, \
                             g_0_xxxxzz_0_xxyyz_0, \
                             g_0_xxxxzz_0_xxyzz_0, \
                             g_0_xxxxzz_0_xxzzz_0, \
                             g_0_xxxxzz_0_xyyyy_0, \
                             g_0_xxxxzz_0_xyyyz_0, \
                             g_0_xxxxzz_0_xyyzz_0, \
                             g_0_xxxxzz_0_xyzzz_0, \
                             g_0_xxxxzz_0_xzzzz_0, \
                             g_0_xxxxzz_0_yyyyy_0, \
                             g_0_xxxxzz_0_yyyyz_0, \
                             g_0_xxxxzz_0_yyyzz_0, \
                             g_0_xxxxzz_0_yyzzz_0, \
                             g_0_xxxxzz_0_yzzzz_0, \
                             g_0_xxxxzz_0_zzzzz_0, \
                             g_0_xxxzz_0_xxxxz_0,  \
                             g_0_xxxzz_0_xxxxz_1,  \
                             g_0_xxxzz_0_xxxyz_0,  \
                             g_0_xxxzz_0_xxxyz_1,  \
                             g_0_xxxzz_0_xxxz_1,   \
                             g_0_xxxzz_0_xxxzz_0,  \
                             g_0_xxxzz_0_xxxzz_1,  \
                             g_0_xxxzz_0_xxyyz_0,  \
                             g_0_xxxzz_0_xxyyz_1,  \
                             g_0_xxxzz_0_xxyz_1,   \
                             g_0_xxxzz_0_xxyzz_0,  \
                             g_0_xxxzz_0_xxyzz_1,  \
                             g_0_xxxzz_0_xxzz_1,   \
                             g_0_xxxzz_0_xxzzz_0,  \
                             g_0_xxxzz_0_xxzzz_1,  \
                             g_0_xxxzz_0_xyyyz_0,  \
                             g_0_xxxzz_0_xyyyz_1,  \
                             g_0_xxxzz_0_xyyz_1,   \
                             g_0_xxxzz_0_xyyzz_0,  \
                             g_0_xxxzz_0_xyyzz_1,  \
                             g_0_xxxzz_0_xyzz_1,   \
                             g_0_xxxzz_0_xyzzz_0,  \
                             g_0_xxxzz_0_xyzzz_1,  \
                             g_0_xxxzz_0_xzzz_1,   \
                             g_0_xxxzz_0_xzzzz_0,  \
                             g_0_xxxzz_0_xzzzz_1,  \
                             g_0_xxxzz_0_yyyyy_0,  \
                             g_0_xxxzz_0_yyyyy_1,  \
                             g_0_xxxzz_0_yyyyz_0,  \
                             g_0_xxxzz_0_yyyyz_1,  \
                             g_0_xxxzz_0_yyyz_1,   \
                             g_0_xxxzz_0_yyyzz_0,  \
                             g_0_xxxzz_0_yyyzz_1,  \
                             g_0_xxxzz_0_yyzz_1,   \
                             g_0_xxxzz_0_yyzzz_0,  \
                             g_0_xxxzz_0_yyzzz_1,  \
                             g_0_xxxzz_0_yzzz_1,   \
                             g_0_xxxzz_0_yzzzz_0,  \
                             g_0_xxxzz_0_yzzzz_1,  \
                             g_0_xxxzz_0_zzzz_1,   \
                             g_0_xxxzz_0_zzzzz_0,  \
                             g_0_xxxzz_0_zzzzz_1,  \
                             g_0_xxzz_0_xxxxz_0,   \
                             g_0_xxzz_0_xxxxz_1,   \
                             g_0_xxzz_0_xxxyz_0,   \
                             g_0_xxzz_0_xxxyz_1,   \
                             g_0_xxzz_0_xxxzz_0,   \
                             g_0_xxzz_0_xxxzz_1,   \
                             g_0_xxzz_0_xxyyz_0,   \
                             g_0_xxzz_0_xxyyz_1,   \
                             g_0_xxzz_0_xxyzz_0,   \
                             g_0_xxzz_0_xxyzz_1,   \
                             g_0_xxzz_0_xxzzz_0,   \
                             g_0_xxzz_0_xxzzz_1,   \
                             g_0_xxzz_0_xyyyz_0,   \
                             g_0_xxzz_0_xyyyz_1,   \
                             g_0_xxzz_0_xyyzz_0,   \
                             g_0_xxzz_0_xyyzz_1,   \
                             g_0_xxzz_0_xyzzz_0,   \
                             g_0_xxzz_0_xyzzz_1,   \
                             g_0_xxzz_0_xzzzz_0,   \
                             g_0_xxzz_0_xzzzz_1,   \
                             g_0_xxzz_0_yyyyy_0,   \
                             g_0_xxzz_0_yyyyy_1,   \
                             g_0_xxzz_0_yyyyz_0,   \
                             g_0_xxzz_0_yyyyz_1,   \
                             g_0_xxzz_0_yyyzz_0,   \
                             g_0_xxzz_0_yyyzz_1,   \
                             g_0_xxzz_0_yyzzz_0,   \
                             g_0_xxzz_0_yyzzz_1,   \
                             g_0_xxzz_0_yzzzz_0,   \
                             g_0_xxzz_0_yzzzz_1,   \
                             g_0_xxzz_0_zzzzz_0,   \
                             g_0_xxzz_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xxxxx_0[i] =
            g_0_xxxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxx_0[i] * pb_z + g_0_xxxxz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxy_0[i] =
            g_0_xxxx_0_xxxxy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxy_0[i] * pb_z + g_0_xxxxz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxz_0[i] = 3.0 * g_0_xxzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxz_0[i] * pb_x + g_0_xxxzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyy_0[i] =
            g_0_xxxx_0_xxxyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxyy_0[i] * pb_z + g_0_xxxxz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxyz_0[i] = 3.0 * g_0_xxzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyz_0[i] * pb_x + g_0_xxxzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxzz_0[i] = 3.0 * g_0_xxzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxzz_0[i] * pb_x + g_0_xxxzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyy_0[i] =
            g_0_xxxx_0_xxyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxyyy_0[i] * pb_z + g_0_xxxxz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxyyz_0[i] = 3.0 * g_0_xxzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyz_0[i] * pb_x + g_0_xxxzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyzz_0[i] = 3.0 * g_0_xxzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzz_0[i] * pb_x + g_0_xxxzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxzzz_0[i] = 3.0 * g_0_xxzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxzzz_0[i] * pb_x + g_0_xxxzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyy_0[i] =
            g_0_xxxx_0_xyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xyyyy_0[i] * pb_z + g_0_xxxxz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xyyyz_0[i] = 3.0 * g_0_xxzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xxxzz_0_xyyyz_0[i] * pb_x + g_0_xxxzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyzz_0[i] = 3.0 * g_0_xxzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xxxzz_0_xyyzz_0[i] * pb_x + g_0_xxxzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyzzz_0[i] = 3.0 * g_0_xxzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxzz_0_xyzzz_0[i] * pb_x + g_0_xxxzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xzzzz_0[i] = 3.0 * g_0_xxzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxzz_0_xzzzz_0[i] * pb_x + g_0_xxxzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyy_0[i] = 3.0 * g_0_xxzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyy_0[i] * pb_x +
                                  g_0_xxxzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyz_0[i] = 3.0 * g_0_xxzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyz_0[i] * pb_x +
                                  g_0_xxxzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyzz_0[i] = 3.0 * g_0_xxzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyzz_0[i] * pb_x +
                                  g_0_xxxzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyzzz_0[i] * pb_x +
                                  g_0_xxxzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yzzzz_0[i] = 3.0 * g_0_xxzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzzzz_0[i] * pb_x +
                                  g_0_xxxzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zzzzz_0[i] = 3.0 * g_0_xxzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzzzz_0[i] * pb_x +
                                  g_0_xxxzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 126-147 components of targeted buffer : SISH

    auto g_0_xxxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 126);

    auto g_0_xxxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 127);

    auto g_0_xxxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 128);

    auto g_0_xxxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 129);

    auto g_0_xxxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 130);

    auto g_0_xxxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 131);

    auto g_0_xxxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 132);

    auto g_0_xxxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 133);

    auto g_0_xxxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 134);

    auto g_0_xxxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 135);

    auto g_0_xxxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 136);

    auto g_0_xxxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 137);

    auto g_0_xxxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 138);

    auto g_0_xxxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 139);

    auto g_0_xxxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 140);

    auto g_0_xxxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 141);

    auto g_0_xxxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 142);

    auto g_0_xxxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 143);

    auto g_0_xxxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 144);

    auto g_0_xxxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 145);

    auto g_0_xxxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 146);

#pragma omp simd aligned(g_0_xxxy_0_xxxxx_0,       \
                             g_0_xxxy_0_xxxxx_1,   \
                             g_0_xxxy_0_xxxxz_0,   \
                             g_0_xxxy_0_xxxxz_1,   \
                             g_0_xxxy_0_xxxzz_0,   \
                             g_0_xxxy_0_xxxzz_1,   \
                             g_0_xxxy_0_xxzzz_0,   \
                             g_0_xxxy_0_xxzzz_1,   \
                             g_0_xxxy_0_xzzzz_0,   \
                             g_0_xxxy_0_xzzzz_1,   \
                             g_0_xxxyy_0_xxxxx_0,  \
                             g_0_xxxyy_0_xxxxx_1,  \
                             g_0_xxxyy_0_xxxxz_0,  \
                             g_0_xxxyy_0_xxxxz_1,  \
                             g_0_xxxyy_0_xxxzz_0,  \
                             g_0_xxxyy_0_xxxzz_1,  \
                             g_0_xxxyy_0_xxzzz_0,  \
                             g_0_xxxyy_0_xxzzz_1,  \
                             g_0_xxxyy_0_xzzzz_0,  \
                             g_0_xxxyy_0_xzzzz_1,  \
                             g_0_xxxyyy_0_xxxxx_0, \
                             g_0_xxxyyy_0_xxxxy_0, \
                             g_0_xxxyyy_0_xxxxz_0, \
                             g_0_xxxyyy_0_xxxyy_0, \
                             g_0_xxxyyy_0_xxxyz_0, \
                             g_0_xxxyyy_0_xxxzz_0, \
                             g_0_xxxyyy_0_xxyyy_0, \
                             g_0_xxxyyy_0_xxyyz_0, \
                             g_0_xxxyyy_0_xxyzz_0, \
                             g_0_xxxyyy_0_xxzzz_0, \
                             g_0_xxxyyy_0_xyyyy_0, \
                             g_0_xxxyyy_0_xyyyz_0, \
                             g_0_xxxyyy_0_xyyzz_0, \
                             g_0_xxxyyy_0_xyzzz_0, \
                             g_0_xxxyyy_0_xzzzz_0, \
                             g_0_xxxyyy_0_yyyyy_0, \
                             g_0_xxxyyy_0_yyyyz_0, \
                             g_0_xxxyyy_0_yyyzz_0, \
                             g_0_xxxyyy_0_yyzzz_0, \
                             g_0_xxxyyy_0_yzzzz_0, \
                             g_0_xxxyyy_0_zzzzz_0, \
                             g_0_xxyyy_0_xxxxy_0,  \
                             g_0_xxyyy_0_xxxxy_1,  \
                             g_0_xxyyy_0_xxxy_1,   \
                             g_0_xxyyy_0_xxxyy_0,  \
                             g_0_xxyyy_0_xxxyy_1,  \
                             g_0_xxyyy_0_xxxyz_0,  \
                             g_0_xxyyy_0_xxxyz_1,  \
                             g_0_xxyyy_0_xxyy_1,   \
                             g_0_xxyyy_0_xxyyy_0,  \
                             g_0_xxyyy_0_xxyyy_1,  \
                             g_0_xxyyy_0_xxyyz_0,  \
                             g_0_xxyyy_0_xxyyz_1,  \
                             g_0_xxyyy_0_xxyz_1,   \
                             g_0_xxyyy_0_xxyzz_0,  \
                             g_0_xxyyy_0_xxyzz_1,  \
                             g_0_xxyyy_0_xyyy_1,   \
                             g_0_xxyyy_0_xyyyy_0,  \
                             g_0_xxyyy_0_xyyyy_1,  \
                             g_0_xxyyy_0_xyyyz_0,  \
                             g_0_xxyyy_0_xyyyz_1,  \
                             g_0_xxyyy_0_xyyz_1,   \
                             g_0_xxyyy_0_xyyzz_0,  \
                             g_0_xxyyy_0_xyyzz_1,  \
                             g_0_xxyyy_0_xyzz_1,   \
                             g_0_xxyyy_0_xyzzz_0,  \
                             g_0_xxyyy_0_xyzzz_1,  \
                             g_0_xxyyy_0_yyyy_1,   \
                             g_0_xxyyy_0_yyyyy_0,  \
                             g_0_xxyyy_0_yyyyy_1,  \
                             g_0_xxyyy_0_yyyyz_0,  \
                             g_0_xxyyy_0_yyyyz_1,  \
                             g_0_xxyyy_0_yyyz_1,   \
                             g_0_xxyyy_0_yyyzz_0,  \
                             g_0_xxyyy_0_yyyzz_1,  \
                             g_0_xxyyy_0_yyzz_1,   \
                             g_0_xxyyy_0_yyzzz_0,  \
                             g_0_xxyyy_0_yyzzz_1,  \
                             g_0_xxyyy_0_yzzz_1,   \
                             g_0_xxyyy_0_yzzzz_0,  \
                             g_0_xxyyy_0_yzzzz_1,  \
                             g_0_xxyyy_0_zzzzz_0,  \
                             g_0_xxyyy_0_zzzzz_1,  \
                             g_0_xyyy_0_xxxxy_0,   \
                             g_0_xyyy_0_xxxxy_1,   \
                             g_0_xyyy_0_xxxyy_0,   \
                             g_0_xyyy_0_xxxyy_1,   \
                             g_0_xyyy_0_xxxyz_0,   \
                             g_0_xyyy_0_xxxyz_1,   \
                             g_0_xyyy_0_xxyyy_0,   \
                             g_0_xyyy_0_xxyyy_1,   \
                             g_0_xyyy_0_xxyyz_0,   \
                             g_0_xyyy_0_xxyyz_1,   \
                             g_0_xyyy_0_xxyzz_0,   \
                             g_0_xyyy_0_xxyzz_1,   \
                             g_0_xyyy_0_xyyyy_0,   \
                             g_0_xyyy_0_xyyyy_1,   \
                             g_0_xyyy_0_xyyyz_0,   \
                             g_0_xyyy_0_xyyyz_1,   \
                             g_0_xyyy_0_xyyzz_0,   \
                             g_0_xyyy_0_xyyzz_1,   \
                             g_0_xyyy_0_xyzzz_0,   \
                             g_0_xyyy_0_xyzzz_1,   \
                             g_0_xyyy_0_yyyyy_0,   \
                             g_0_xyyy_0_yyyyy_1,   \
                             g_0_xyyy_0_yyyyz_0,   \
                             g_0_xyyy_0_yyyyz_1,   \
                             g_0_xyyy_0_yyyzz_0,   \
                             g_0_xyyy_0_yyyzz_1,   \
                             g_0_xyyy_0_yyzzz_0,   \
                             g_0_xyyy_0_yyzzz_1,   \
                             g_0_xyyy_0_yzzzz_0,   \
                             g_0_xyyy_0_yzzzz_1,   \
                             g_0_xyyy_0_zzzzz_0,   \
                             g_0_xyyy_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xxxxx_0[i] = 2.0 * g_0_xxxy_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxx_0[i] * pb_y +
                                  g_0_xxxyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxy_0[i] = 2.0 * g_0_xyyy_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxy_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxy_0[i] * pb_x + g_0_xxyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxz_0[i] = 2.0 * g_0_xxxy_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxz_0[i] * pb_y +
                                  g_0_xxxyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxyy_0[i] = 2.0 * g_0_xyyy_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyy_0[i] * pb_x + g_0_xxyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyz_0[i] = 2.0 * g_0_xyyy_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyz_0[i] * pb_x + g_0_xxyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxzz_0[i] = 2.0 * g_0_xxxy_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxzz_0[i] * pb_y +
                                  g_0_xxxyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxyyy_0[i] = 2.0 * g_0_xyyy_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyy_0[i] * pb_x + g_0_xxyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyz_0[i] = 2.0 * g_0_xyyy_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyz_0[i] * pb_x + g_0_xxyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyzz_0[i] = 2.0 * g_0_xyyy_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzz_0[i] * pb_x + g_0_xxyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxzzz_0[i] = 2.0 * g_0_xxxy_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxzzz_0[i] * pb_y +
                                  g_0_xxxyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xyyyy_0[i] = 2.0 * g_0_xyyy_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyy_1[i] * fi_abcd_0 +
                                  g_0_xxyyy_0_xyyyy_0[i] * pb_x + g_0_xxyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyz_0[i] = 2.0 * g_0_xyyy_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xxyyy_0_xyyyz_0[i] * pb_x + g_0_xxyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyzz_0[i] = 2.0 * g_0_xyyy_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xxyyy_0_xyyzz_0[i] * pb_x + g_0_xxyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyzzz_0[i] = 2.0 * g_0_xyyy_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xxyyy_0_xyzzz_0[i] * pb_x + g_0_xxyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xzzzz_0[i] = 2.0 * g_0_xxxy_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xzzzz_0[i] * pb_y +
                                  g_0_xxxyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yyyyy_0[i] = 2.0 * g_0_xyyy_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyy_0[i] * pb_x +
                                  g_0_xxyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyz_0[i] = 2.0 * g_0_xyyy_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyz_0[i] * pb_x +
                                  g_0_xxyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyzz_0[i] = 2.0 * g_0_xyyy_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyzz_0[i] * pb_x +
                                  g_0_xxyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyzzz_0[i] = 2.0 * g_0_xyyy_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyzzz_0[i] * pb_x +
                                  g_0_xxyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yzzzz_0[i] = 2.0 * g_0_xyyy_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzzzz_0[i] * pb_x +
                                  g_0_xxyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zzzzz_0[i] = 2.0 * g_0_xyyy_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_zzzzz_0[i] * pb_x +
                                  g_0_xxyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 147-168 components of targeted buffer : SISH

    auto g_0_xxxyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 147);

    auto g_0_xxxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 148);

    auto g_0_xxxyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 149);

    auto g_0_xxxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 150);

    auto g_0_xxxyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 151);

    auto g_0_xxxyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 152);

    auto g_0_xxxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 153);

    auto g_0_xxxyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 154);

    auto g_0_xxxyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 155);

    auto g_0_xxxyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 156);

    auto g_0_xxxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 157);

    auto g_0_xxxyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 158);

    auto g_0_xxxyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 159);

    auto g_0_xxxyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 160);

    auto g_0_xxxyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 161);

    auto g_0_xxxyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 162);

    auto g_0_xxxyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 163);

    auto g_0_xxxyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 164);

    auto g_0_xxxyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 165);

    auto g_0_xxxyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 166);

    auto g_0_xxxyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 167);

#pragma omp simd aligned(g_0_xxxyy_0_xxxx_1,       \
                             g_0_xxxyy_0_xxxxx_0,  \
                             g_0_xxxyy_0_xxxxx_1,  \
                             g_0_xxxyy_0_xxxxy_0,  \
                             g_0_xxxyy_0_xxxxy_1,  \
                             g_0_xxxyy_0_xxxxz_0,  \
                             g_0_xxxyy_0_xxxxz_1,  \
                             g_0_xxxyy_0_xxxy_1,   \
                             g_0_xxxyy_0_xxxyy_0,  \
                             g_0_xxxyy_0_xxxyy_1,  \
                             g_0_xxxyy_0_xxxyz_0,  \
                             g_0_xxxyy_0_xxxyz_1,  \
                             g_0_xxxyy_0_xxxz_1,   \
                             g_0_xxxyy_0_xxxzz_0,  \
                             g_0_xxxyy_0_xxxzz_1,  \
                             g_0_xxxyy_0_xxyy_1,   \
                             g_0_xxxyy_0_xxyyy_0,  \
                             g_0_xxxyy_0_xxyyy_1,  \
                             g_0_xxxyy_0_xxyyz_0,  \
                             g_0_xxxyy_0_xxyyz_1,  \
                             g_0_xxxyy_0_xxyz_1,   \
                             g_0_xxxyy_0_xxyzz_0,  \
                             g_0_xxxyy_0_xxyzz_1,  \
                             g_0_xxxyy_0_xxzz_1,   \
                             g_0_xxxyy_0_xxzzz_0,  \
                             g_0_xxxyy_0_xxzzz_1,  \
                             g_0_xxxyy_0_xyyy_1,   \
                             g_0_xxxyy_0_xyyyy_0,  \
                             g_0_xxxyy_0_xyyyy_1,  \
                             g_0_xxxyy_0_xyyyz_0,  \
                             g_0_xxxyy_0_xyyyz_1,  \
                             g_0_xxxyy_0_xyyz_1,   \
                             g_0_xxxyy_0_xyyzz_0,  \
                             g_0_xxxyy_0_xyyzz_1,  \
                             g_0_xxxyy_0_xyzz_1,   \
                             g_0_xxxyy_0_xyzzz_0,  \
                             g_0_xxxyy_0_xyzzz_1,  \
                             g_0_xxxyy_0_xzzz_1,   \
                             g_0_xxxyy_0_xzzzz_0,  \
                             g_0_xxxyy_0_xzzzz_1,  \
                             g_0_xxxyy_0_yyyy_1,   \
                             g_0_xxxyy_0_yyyyy_0,  \
                             g_0_xxxyy_0_yyyyy_1,  \
                             g_0_xxxyy_0_yyyyz_0,  \
                             g_0_xxxyy_0_yyyyz_1,  \
                             g_0_xxxyy_0_yyyz_1,   \
                             g_0_xxxyy_0_yyyzz_0,  \
                             g_0_xxxyy_0_yyyzz_1,  \
                             g_0_xxxyy_0_yyzz_1,   \
                             g_0_xxxyy_0_yyzzz_0,  \
                             g_0_xxxyy_0_yyzzz_1,  \
                             g_0_xxxyy_0_yzzz_1,   \
                             g_0_xxxyy_0_yzzzz_0,  \
                             g_0_xxxyy_0_yzzzz_1,  \
                             g_0_xxxyy_0_zzzz_1,   \
                             g_0_xxxyy_0_zzzzz_0,  \
                             g_0_xxxyy_0_zzzzz_1,  \
                             g_0_xxxyyz_0_xxxxx_0, \
                             g_0_xxxyyz_0_xxxxy_0, \
                             g_0_xxxyyz_0_xxxxz_0, \
                             g_0_xxxyyz_0_xxxyy_0, \
                             g_0_xxxyyz_0_xxxyz_0, \
                             g_0_xxxyyz_0_xxxzz_0, \
                             g_0_xxxyyz_0_xxyyy_0, \
                             g_0_xxxyyz_0_xxyyz_0, \
                             g_0_xxxyyz_0_xxyzz_0, \
                             g_0_xxxyyz_0_xxzzz_0, \
                             g_0_xxxyyz_0_xyyyy_0, \
                             g_0_xxxyyz_0_xyyyz_0, \
                             g_0_xxxyyz_0_xyyzz_0, \
                             g_0_xxxyyz_0_xyzzz_0, \
                             g_0_xxxyyz_0_xzzzz_0, \
                             g_0_xxxyyz_0_yyyyy_0, \
                             g_0_xxxyyz_0_yyyyz_0, \
                             g_0_xxxyyz_0_yyyzz_0, \
                             g_0_xxxyyz_0_yyzzz_0, \
                             g_0_xxxyyz_0_yzzzz_0, \
                             g_0_xxxyyz_0_zzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xxxxx_0[i] = g_0_xxxyy_0_xxxxx_0[i] * pb_z + g_0_xxxyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxy_0[i] = g_0_xxxyy_0_xxxxy_0[i] * pb_z + g_0_xxxyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxz_0[i] = g_0_xxxyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxz_0[i] * pb_z + g_0_xxxyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyy_0[i] = g_0_xxxyy_0_xxxyy_0[i] * pb_z + g_0_xxxyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyz_0[i] = g_0_xxxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyz_0[i] * pb_z + g_0_xxxyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxzz_0[i] = 2.0 * g_0_xxxyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxzz_0[i] * pb_z + g_0_xxxyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyy_0[i] = g_0_xxxyy_0_xxyyy_0[i] * pb_z + g_0_xxxyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyz_0[i] = g_0_xxxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyz_0[i] * pb_z + g_0_xxxyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyzz_0[i] = 2.0 * g_0_xxxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzz_0[i] * pb_z + g_0_xxxyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxzzz_0[i] = 3.0 * g_0_xxxyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxzzz_0[i] * pb_z + g_0_xxxyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyy_0[i] = g_0_xxxyy_0_xyyyy_0[i] * pb_z + g_0_xxxyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyz_0[i] = g_0_xxxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyz_0[i] * pb_z + g_0_xxxyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyzz_0[i] = 2.0 * g_0_xxxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzz_0[i] * pb_z + g_0_xxxyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyzzz_0[i] = 3.0 * g_0_xxxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzz_0[i] * pb_z + g_0_xxxyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xzzzz_0[i] = 4.0 * g_0_xxxyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xzzzz_0[i] * pb_z + g_0_xxxyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyy_0[i] = g_0_xxxyy_0_yyyyy_0[i] * pb_z + g_0_xxxyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyz_0[i] = g_0_xxxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyz_0[i] * pb_z + g_0_xxxyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyzz_0[i] = 2.0 * g_0_xxxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyzz_0[i] * pb_z + g_0_xxxyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyzzz_0[i] = 3.0 * g_0_xxxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyzzz_0[i] * pb_z + g_0_xxxyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yzzzz_0[i] = 4.0 * g_0_xxxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yzzzz_0[i] * pb_z + g_0_xxxyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zzzzz_0[i] = 5.0 * g_0_xxxyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_zzzzz_0[i] * pb_z + g_0_xxxyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 168-189 components of targeted buffer : SISH

    auto g_0_xxxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 168);

    auto g_0_xxxyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 169);

    auto g_0_xxxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 170);

    auto g_0_xxxyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 171);

    auto g_0_xxxyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 172);

    auto g_0_xxxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 173);

    auto g_0_xxxyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 174);

    auto g_0_xxxyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 175);

    auto g_0_xxxyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 176);

    auto g_0_xxxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 177);

    auto g_0_xxxyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 178);

    auto g_0_xxxyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 179);

    auto g_0_xxxyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 180);

    auto g_0_xxxyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 181);

    auto g_0_xxxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 182);

    auto g_0_xxxyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 183);

    auto g_0_xxxyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 184);

    auto g_0_xxxyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 185);

    auto g_0_xxxyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 186);

    auto g_0_xxxyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 187);

    auto g_0_xxxyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 188);

#pragma omp simd aligned(g_0_xxxyzz_0_xxxxx_0,     \
                             g_0_xxxyzz_0_xxxxy_0, \
                             g_0_xxxyzz_0_xxxxz_0, \
                             g_0_xxxyzz_0_xxxyy_0, \
                             g_0_xxxyzz_0_xxxyz_0, \
                             g_0_xxxyzz_0_xxxzz_0, \
                             g_0_xxxyzz_0_xxyyy_0, \
                             g_0_xxxyzz_0_xxyyz_0, \
                             g_0_xxxyzz_0_xxyzz_0, \
                             g_0_xxxyzz_0_xxzzz_0, \
                             g_0_xxxyzz_0_xyyyy_0, \
                             g_0_xxxyzz_0_xyyyz_0, \
                             g_0_xxxyzz_0_xyyzz_0, \
                             g_0_xxxyzz_0_xyzzz_0, \
                             g_0_xxxyzz_0_xzzzz_0, \
                             g_0_xxxyzz_0_yyyyy_0, \
                             g_0_xxxyzz_0_yyyyz_0, \
                             g_0_xxxyzz_0_yyyzz_0, \
                             g_0_xxxyzz_0_yyzzz_0, \
                             g_0_xxxyzz_0_yzzzz_0, \
                             g_0_xxxyzz_0_zzzzz_0, \
                             g_0_xxxzz_0_xxxx_1,   \
                             g_0_xxxzz_0_xxxxx_0,  \
                             g_0_xxxzz_0_xxxxx_1,  \
                             g_0_xxxzz_0_xxxxy_0,  \
                             g_0_xxxzz_0_xxxxy_1,  \
                             g_0_xxxzz_0_xxxxz_0,  \
                             g_0_xxxzz_0_xxxxz_1,  \
                             g_0_xxxzz_0_xxxy_1,   \
                             g_0_xxxzz_0_xxxyy_0,  \
                             g_0_xxxzz_0_xxxyy_1,  \
                             g_0_xxxzz_0_xxxyz_0,  \
                             g_0_xxxzz_0_xxxyz_1,  \
                             g_0_xxxzz_0_xxxz_1,   \
                             g_0_xxxzz_0_xxxzz_0,  \
                             g_0_xxxzz_0_xxxzz_1,  \
                             g_0_xxxzz_0_xxyy_1,   \
                             g_0_xxxzz_0_xxyyy_0,  \
                             g_0_xxxzz_0_xxyyy_1,  \
                             g_0_xxxzz_0_xxyyz_0,  \
                             g_0_xxxzz_0_xxyyz_1,  \
                             g_0_xxxzz_0_xxyz_1,   \
                             g_0_xxxzz_0_xxyzz_0,  \
                             g_0_xxxzz_0_xxyzz_1,  \
                             g_0_xxxzz_0_xxzz_1,   \
                             g_0_xxxzz_0_xxzzz_0,  \
                             g_0_xxxzz_0_xxzzz_1,  \
                             g_0_xxxzz_0_xyyy_1,   \
                             g_0_xxxzz_0_xyyyy_0,  \
                             g_0_xxxzz_0_xyyyy_1,  \
                             g_0_xxxzz_0_xyyyz_0,  \
                             g_0_xxxzz_0_xyyyz_1,  \
                             g_0_xxxzz_0_xyyz_1,   \
                             g_0_xxxzz_0_xyyzz_0,  \
                             g_0_xxxzz_0_xyyzz_1,  \
                             g_0_xxxzz_0_xyzz_1,   \
                             g_0_xxxzz_0_xyzzz_0,  \
                             g_0_xxxzz_0_xyzzz_1,  \
                             g_0_xxxzz_0_xzzz_1,   \
                             g_0_xxxzz_0_xzzzz_0,  \
                             g_0_xxxzz_0_xzzzz_1,  \
                             g_0_xxxzz_0_yyyy_1,   \
                             g_0_xxxzz_0_yyyyy_0,  \
                             g_0_xxxzz_0_yyyyy_1,  \
                             g_0_xxxzz_0_yyyyz_0,  \
                             g_0_xxxzz_0_yyyyz_1,  \
                             g_0_xxxzz_0_yyyz_1,   \
                             g_0_xxxzz_0_yyyzz_0,  \
                             g_0_xxxzz_0_yyyzz_1,  \
                             g_0_xxxzz_0_yyzz_1,   \
                             g_0_xxxzz_0_yyzzz_0,  \
                             g_0_xxxzz_0_yyzzz_1,  \
                             g_0_xxxzz_0_yzzz_1,   \
                             g_0_xxxzz_0_yzzzz_0,  \
                             g_0_xxxzz_0_yzzzz_1,  \
                             g_0_xxxzz_0_zzzz_1,   \
                             g_0_xxxzz_0_zzzzz_0,  \
                             g_0_xxxzz_0_zzzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xxxxx_0[i] = g_0_xxxzz_0_xxxxx_0[i] * pb_y + g_0_xxxzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxy_0[i] = g_0_xxxzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxy_0[i] * pb_y + g_0_xxxzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxz_0[i] = g_0_xxxzz_0_xxxxz_0[i] * pb_y + g_0_xxxzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyy_0[i] = 2.0 * g_0_xxxzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyy_0[i] * pb_y + g_0_xxxzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyz_0[i] = g_0_xxxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyz_0[i] * pb_y + g_0_xxxzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxzz_0[i] = g_0_xxxzz_0_xxxzz_0[i] * pb_y + g_0_xxxzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyy_0[i] = 3.0 * g_0_xxxzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyy_0[i] * pb_y + g_0_xxxzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyz_0[i] = 2.0 * g_0_xxxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyz_0[i] * pb_y + g_0_xxxzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyzz_0[i] = g_0_xxxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzz_0[i] * pb_y + g_0_xxxzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxzzz_0[i] = g_0_xxxzz_0_xxzzz_0[i] * pb_y + g_0_xxxzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyy_0[i] = 4.0 * g_0_xxxzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyy_0[i] * pb_y + g_0_xxxzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyz_0[i] = 3.0 * g_0_xxxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyz_0[i] * pb_y + g_0_xxxzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyzz_0[i] = 2.0 * g_0_xxxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzz_0[i] * pb_y + g_0_xxxzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyzzz_0[i] = g_0_xxxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzz_0[i] * pb_y + g_0_xxxzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xzzzz_0[i] = g_0_xxxzz_0_xzzzz_0[i] * pb_y + g_0_xxxzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyy_0[i] = 5.0 * g_0_xxxzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyy_0[i] * pb_y + g_0_xxxzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyz_0[i] = 4.0 * g_0_xxxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyz_0[i] * pb_y + g_0_xxxzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyzz_0[i] = 3.0 * g_0_xxxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyzz_0[i] * pb_y + g_0_xxxzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyzzz_0[i] = 2.0 * g_0_xxxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyzzz_0[i] * pb_y + g_0_xxxzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yzzzz_0[i] = g_0_xxxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yzzzz_0[i] * pb_y + g_0_xxxzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zzzzz_0[i] = g_0_xxxzz_0_zzzzz_0[i] * pb_y + g_0_xxxzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 189-210 components of targeted buffer : SISH

    auto g_0_xxxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 189);

    auto g_0_xxxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 190);

    auto g_0_xxxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 191);

    auto g_0_xxxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 192);

    auto g_0_xxxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 193);

    auto g_0_xxxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 194);

    auto g_0_xxxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 195);

    auto g_0_xxxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 196);

    auto g_0_xxxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 197);

    auto g_0_xxxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 198);

    auto g_0_xxxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 199);

    auto g_0_xxxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 200);

    auto g_0_xxxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 201);

    auto g_0_xxxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 202);

    auto g_0_xxxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 203);

    auto g_0_xxxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 204);

    auto g_0_xxxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 205);

    auto g_0_xxxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 206);

    auto g_0_xxxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 207);

    auto g_0_xxxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 208);

    auto g_0_xxxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 209);

#pragma omp simd aligned(g_0_xxxz_0_xxxxx_0,       \
                             g_0_xxxz_0_xxxxx_1,   \
                             g_0_xxxz_0_xxxxy_0,   \
                             g_0_xxxz_0_xxxxy_1,   \
                             g_0_xxxz_0_xxxyy_0,   \
                             g_0_xxxz_0_xxxyy_1,   \
                             g_0_xxxz_0_xxyyy_0,   \
                             g_0_xxxz_0_xxyyy_1,   \
                             g_0_xxxz_0_xyyyy_0,   \
                             g_0_xxxz_0_xyyyy_1,   \
                             g_0_xxxzz_0_xxxxx_0,  \
                             g_0_xxxzz_0_xxxxx_1,  \
                             g_0_xxxzz_0_xxxxy_0,  \
                             g_0_xxxzz_0_xxxxy_1,  \
                             g_0_xxxzz_0_xxxyy_0,  \
                             g_0_xxxzz_0_xxxyy_1,  \
                             g_0_xxxzz_0_xxyyy_0,  \
                             g_0_xxxzz_0_xxyyy_1,  \
                             g_0_xxxzz_0_xyyyy_0,  \
                             g_0_xxxzz_0_xyyyy_1,  \
                             g_0_xxxzzz_0_xxxxx_0, \
                             g_0_xxxzzz_0_xxxxy_0, \
                             g_0_xxxzzz_0_xxxxz_0, \
                             g_0_xxxzzz_0_xxxyy_0, \
                             g_0_xxxzzz_0_xxxyz_0, \
                             g_0_xxxzzz_0_xxxzz_0, \
                             g_0_xxxzzz_0_xxyyy_0, \
                             g_0_xxxzzz_0_xxyyz_0, \
                             g_0_xxxzzz_0_xxyzz_0, \
                             g_0_xxxzzz_0_xxzzz_0, \
                             g_0_xxxzzz_0_xyyyy_0, \
                             g_0_xxxzzz_0_xyyyz_0, \
                             g_0_xxxzzz_0_xyyzz_0, \
                             g_0_xxxzzz_0_xyzzz_0, \
                             g_0_xxxzzz_0_xzzzz_0, \
                             g_0_xxxzzz_0_yyyyy_0, \
                             g_0_xxxzzz_0_yyyyz_0, \
                             g_0_xxxzzz_0_yyyzz_0, \
                             g_0_xxxzzz_0_yyzzz_0, \
                             g_0_xxxzzz_0_yzzzz_0, \
                             g_0_xxxzzz_0_zzzzz_0, \
                             g_0_xxzzz_0_xxxxz_0,  \
                             g_0_xxzzz_0_xxxxz_1,  \
                             g_0_xxzzz_0_xxxyz_0,  \
                             g_0_xxzzz_0_xxxyz_1,  \
                             g_0_xxzzz_0_xxxz_1,   \
                             g_0_xxzzz_0_xxxzz_0,  \
                             g_0_xxzzz_0_xxxzz_1,  \
                             g_0_xxzzz_0_xxyyz_0,  \
                             g_0_xxzzz_0_xxyyz_1,  \
                             g_0_xxzzz_0_xxyz_1,   \
                             g_0_xxzzz_0_xxyzz_0,  \
                             g_0_xxzzz_0_xxyzz_1,  \
                             g_0_xxzzz_0_xxzz_1,   \
                             g_0_xxzzz_0_xxzzz_0,  \
                             g_0_xxzzz_0_xxzzz_1,  \
                             g_0_xxzzz_0_xyyyz_0,  \
                             g_0_xxzzz_0_xyyyz_1,  \
                             g_0_xxzzz_0_xyyz_1,   \
                             g_0_xxzzz_0_xyyzz_0,  \
                             g_0_xxzzz_0_xyyzz_1,  \
                             g_0_xxzzz_0_xyzz_1,   \
                             g_0_xxzzz_0_xyzzz_0,  \
                             g_0_xxzzz_0_xyzzz_1,  \
                             g_0_xxzzz_0_xzzz_1,   \
                             g_0_xxzzz_0_xzzzz_0,  \
                             g_0_xxzzz_0_xzzzz_1,  \
                             g_0_xxzzz_0_yyyyy_0,  \
                             g_0_xxzzz_0_yyyyy_1,  \
                             g_0_xxzzz_0_yyyyz_0,  \
                             g_0_xxzzz_0_yyyyz_1,  \
                             g_0_xxzzz_0_yyyz_1,   \
                             g_0_xxzzz_0_yyyzz_0,  \
                             g_0_xxzzz_0_yyyzz_1,  \
                             g_0_xxzzz_0_yyzz_1,   \
                             g_0_xxzzz_0_yyzzz_0,  \
                             g_0_xxzzz_0_yyzzz_1,  \
                             g_0_xxzzz_0_yzzz_1,   \
                             g_0_xxzzz_0_yzzzz_0,  \
                             g_0_xxzzz_0_yzzzz_1,  \
                             g_0_xxzzz_0_zzzz_1,   \
                             g_0_xxzzz_0_zzzzz_0,  \
                             g_0_xxzzz_0_zzzzz_1,  \
                             g_0_xzzz_0_xxxxz_0,   \
                             g_0_xzzz_0_xxxxz_1,   \
                             g_0_xzzz_0_xxxyz_0,   \
                             g_0_xzzz_0_xxxyz_1,   \
                             g_0_xzzz_0_xxxzz_0,   \
                             g_0_xzzz_0_xxxzz_1,   \
                             g_0_xzzz_0_xxyyz_0,   \
                             g_0_xzzz_0_xxyyz_1,   \
                             g_0_xzzz_0_xxyzz_0,   \
                             g_0_xzzz_0_xxyzz_1,   \
                             g_0_xzzz_0_xxzzz_0,   \
                             g_0_xzzz_0_xxzzz_1,   \
                             g_0_xzzz_0_xyyyz_0,   \
                             g_0_xzzz_0_xyyyz_1,   \
                             g_0_xzzz_0_xyyzz_0,   \
                             g_0_xzzz_0_xyyzz_1,   \
                             g_0_xzzz_0_xyzzz_0,   \
                             g_0_xzzz_0_xyzzz_1,   \
                             g_0_xzzz_0_xzzzz_0,   \
                             g_0_xzzz_0_xzzzz_1,   \
                             g_0_xzzz_0_yyyyy_0,   \
                             g_0_xzzz_0_yyyyy_1,   \
                             g_0_xzzz_0_yyyyz_0,   \
                             g_0_xzzz_0_yyyyz_1,   \
                             g_0_xzzz_0_yyyzz_0,   \
                             g_0_xzzz_0_yyyzz_1,   \
                             g_0_xzzz_0_yyzzz_0,   \
                             g_0_xzzz_0_yyzzz_1,   \
                             g_0_xzzz_0_yzzzz_0,   \
                             g_0_xzzz_0_yzzzz_1,   \
                             g_0_xzzz_0_zzzzz_0,   \
                             g_0_xzzz_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xxxxx_0[i] = 2.0 * g_0_xxxz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxx_0[i] * pb_z +
                                  g_0_xxxzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxy_0[i] = 2.0 * g_0_xxxz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxy_0[i] * pb_z +
                                  g_0_xxxzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxz_0[i] = 2.0 * g_0_xzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxz_0[i] * pb_x + g_0_xxzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxyy_0[i] * pb_z +
                                  g_0_xxxzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxyz_0[i] = 2.0 * g_0_xzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyz_0[i] * pb_x + g_0_xxzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxzz_0[i] = 2.0 * g_0_xzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxzz_0[i] * pb_x + g_0_xxzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyy_0[i] = 2.0 * g_0_xxxz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxyyy_0[i] * pb_z +
                                  g_0_xxxzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxyyz_0[i] = 2.0 * g_0_xzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyz_0[i] * pb_x + g_0_xxzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyzz_0[i] = 2.0 * g_0_xzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzz_0[i] * pb_x + g_0_xxzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxzzz_0[i] = 2.0 * g_0_xzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxzzz_0[i] * pb_x + g_0_xxzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyy_0[i] = 2.0 * g_0_xxxz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xyyyy_0[i] * pb_z +
                                  g_0_xxxzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xyyyz_0[i] = 2.0 * g_0_xzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xxzzz_0_xyyyz_0[i] * pb_x + g_0_xxzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyzz_0[i] = 2.0 * g_0_xzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xxzzz_0_xyyzz_0[i] * pb_x + g_0_xxzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyzzz_0[i] = 2.0 * g_0_xzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xxzzz_0_xyzzz_0[i] * pb_x + g_0_xxzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xzzzz_0[i] = 2.0 * g_0_xzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_xxzzz_0_xzzzz_0[i] * pb_x + g_0_xxzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyy_0[i] = 2.0 * g_0_xzzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyy_0[i] * pb_x +
                                  g_0_xxzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyz_0[i] = 2.0 * g_0_xzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyz_0[i] * pb_x +
                                  g_0_xxzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyzz_0[i] = 2.0 * g_0_xzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyzz_0[i] * pb_x +
                                  g_0_xxzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyzzz_0[i] = 2.0 * g_0_xzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyzzz_0[i] * pb_x +
                                  g_0_xxzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yzzzz_0[i] = 2.0 * g_0_xzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzzzz_0[i] * pb_x +
                                  g_0_xxzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zzzzz_0[i] = 2.0 * g_0_xzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzzzz_0[i] * pb_x +
                                  g_0_xxzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 210-231 components of targeted buffer : SISH

    auto g_0_xxyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 210);

    auto g_0_xxyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 211);

    auto g_0_xxyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 212);

    auto g_0_xxyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 213);

    auto g_0_xxyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 214);

    auto g_0_xxyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 215);

    auto g_0_xxyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 216);

    auto g_0_xxyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 217);

    auto g_0_xxyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 218);

    auto g_0_xxyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 219);

    auto g_0_xxyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 220);

    auto g_0_xxyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 221);

    auto g_0_xxyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 222);

    auto g_0_xxyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 223);

    auto g_0_xxyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 224);

    auto g_0_xxyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 225);

    auto g_0_xxyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 226);

    auto g_0_xxyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 227);

    auto g_0_xxyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 228);

    auto g_0_xxyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 229);

    auto g_0_xxyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 230);

#pragma omp simd aligned(g_0_xxyy_0_xxxxx_0,       \
                             g_0_xxyy_0_xxxxx_1,   \
                             g_0_xxyy_0_xxxxz_0,   \
                             g_0_xxyy_0_xxxxz_1,   \
                             g_0_xxyy_0_xxxzz_0,   \
                             g_0_xxyy_0_xxxzz_1,   \
                             g_0_xxyy_0_xxzzz_0,   \
                             g_0_xxyy_0_xxzzz_1,   \
                             g_0_xxyy_0_xzzzz_0,   \
                             g_0_xxyy_0_xzzzz_1,   \
                             g_0_xxyyy_0_xxxxx_0,  \
                             g_0_xxyyy_0_xxxxx_1,  \
                             g_0_xxyyy_0_xxxxz_0,  \
                             g_0_xxyyy_0_xxxxz_1,  \
                             g_0_xxyyy_0_xxxzz_0,  \
                             g_0_xxyyy_0_xxxzz_1,  \
                             g_0_xxyyy_0_xxzzz_0,  \
                             g_0_xxyyy_0_xxzzz_1,  \
                             g_0_xxyyy_0_xzzzz_0,  \
                             g_0_xxyyy_0_xzzzz_1,  \
                             g_0_xxyyyy_0_xxxxx_0, \
                             g_0_xxyyyy_0_xxxxy_0, \
                             g_0_xxyyyy_0_xxxxz_0, \
                             g_0_xxyyyy_0_xxxyy_0, \
                             g_0_xxyyyy_0_xxxyz_0, \
                             g_0_xxyyyy_0_xxxzz_0, \
                             g_0_xxyyyy_0_xxyyy_0, \
                             g_0_xxyyyy_0_xxyyz_0, \
                             g_0_xxyyyy_0_xxyzz_0, \
                             g_0_xxyyyy_0_xxzzz_0, \
                             g_0_xxyyyy_0_xyyyy_0, \
                             g_0_xxyyyy_0_xyyyz_0, \
                             g_0_xxyyyy_0_xyyzz_0, \
                             g_0_xxyyyy_0_xyzzz_0, \
                             g_0_xxyyyy_0_xzzzz_0, \
                             g_0_xxyyyy_0_yyyyy_0, \
                             g_0_xxyyyy_0_yyyyz_0, \
                             g_0_xxyyyy_0_yyyzz_0, \
                             g_0_xxyyyy_0_yyzzz_0, \
                             g_0_xxyyyy_0_yzzzz_0, \
                             g_0_xxyyyy_0_zzzzz_0, \
                             g_0_xyyyy_0_xxxxy_0,  \
                             g_0_xyyyy_0_xxxxy_1,  \
                             g_0_xyyyy_0_xxxy_1,   \
                             g_0_xyyyy_0_xxxyy_0,  \
                             g_0_xyyyy_0_xxxyy_1,  \
                             g_0_xyyyy_0_xxxyz_0,  \
                             g_0_xyyyy_0_xxxyz_1,  \
                             g_0_xyyyy_0_xxyy_1,   \
                             g_0_xyyyy_0_xxyyy_0,  \
                             g_0_xyyyy_0_xxyyy_1,  \
                             g_0_xyyyy_0_xxyyz_0,  \
                             g_0_xyyyy_0_xxyyz_1,  \
                             g_0_xyyyy_0_xxyz_1,   \
                             g_0_xyyyy_0_xxyzz_0,  \
                             g_0_xyyyy_0_xxyzz_1,  \
                             g_0_xyyyy_0_xyyy_1,   \
                             g_0_xyyyy_0_xyyyy_0,  \
                             g_0_xyyyy_0_xyyyy_1,  \
                             g_0_xyyyy_0_xyyyz_0,  \
                             g_0_xyyyy_0_xyyyz_1,  \
                             g_0_xyyyy_0_xyyz_1,   \
                             g_0_xyyyy_0_xyyzz_0,  \
                             g_0_xyyyy_0_xyyzz_1,  \
                             g_0_xyyyy_0_xyzz_1,   \
                             g_0_xyyyy_0_xyzzz_0,  \
                             g_0_xyyyy_0_xyzzz_1,  \
                             g_0_xyyyy_0_yyyy_1,   \
                             g_0_xyyyy_0_yyyyy_0,  \
                             g_0_xyyyy_0_yyyyy_1,  \
                             g_0_xyyyy_0_yyyyz_0,  \
                             g_0_xyyyy_0_yyyyz_1,  \
                             g_0_xyyyy_0_yyyz_1,   \
                             g_0_xyyyy_0_yyyzz_0,  \
                             g_0_xyyyy_0_yyyzz_1,  \
                             g_0_xyyyy_0_yyzz_1,   \
                             g_0_xyyyy_0_yyzzz_0,  \
                             g_0_xyyyy_0_yyzzz_1,  \
                             g_0_xyyyy_0_yzzz_1,   \
                             g_0_xyyyy_0_yzzzz_0,  \
                             g_0_xyyyy_0_yzzzz_1,  \
                             g_0_xyyyy_0_zzzzz_0,  \
                             g_0_xyyyy_0_zzzzz_1,  \
                             g_0_yyyy_0_xxxxy_0,   \
                             g_0_yyyy_0_xxxxy_1,   \
                             g_0_yyyy_0_xxxyy_0,   \
                             g_0_yyyy_0_xxxyy_1,   \
                             g_0_yyyy_0_xxxyz_0,   \
                             g_0_yyyy_0_xxxyz_1,   \
                             g_0_yyyy_0_xxyyy_0,   \
                             g_0_yyyy_0_xxyyy_1,   \
                             g_0_yyyy_0_xxyyz_0,   \
                             g_0_yyyy_0_xxyyz_1,   \
                             g_0_yyyy_0_xxyzz_0,   \
                             g_0_yyyy_0_xxyzz_1,   \
                             g_0_yyyy_0_xyyyy_0,   \
                             g_0_yyyy_0_xyyyy_1,   \
                             g_0_yyyy_0_xyyyz_0,   \
                             g_0_yyyy_0_xyyyz_1,   \
                             g_0_yyyy_0_xyyzz_0,   \
                             g_0_yyyy_0_xyyzz_1,   \
                             g_0_yyyy_0_xyzzz_0,   \
                             g_0_yyyy_0_xyzzz_1,   \
                             g_0_yyyy_0_yyyyy_0,   \
                             g_0_yyyy_0_yyyyy_1,   \
                             g_0_yyyy_0_yyyyz_0,   \
                             g_0_yyyy_0_yyyyz_1,   \
                             g_0_yyyy_0_yyyzz_0,   \
                             g_0_yyyy_0_yyyzz_1,   \
                             g_0_yyyy_0_yyzzz_0,   \
                             g_0_yyyy_0_yyzzz_1,   \
                             g_0_yyyy_0_yzzzz_0,   \
                             g_0_yyyy_0_yzzzz_1,   \
                             g_0_yyyy_0_zzzzz_0,   \
                             g_0_yyyy_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xxxxx_0[i] = 3.0 * g_0_xxyy_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxx_0[i] * pb_y +
                                  g_0_xxyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxy_0[i] = g_0_yyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyy_0_xxxy_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xxxxy_0[i] * pb_x + g_0_xyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxz_0[i] = 3.0 * g_0_xxyy_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxz_0[i] * pb_y +
                                  g_0_xxyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxyy_0[i] = g_0_yyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyy_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xxxyy_0[i] * pb_x + g_0_xyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyz_0[i] = g_0_yyyy_0_xxxyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xxxyz_0[i] * pb_x + g_0_xyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxzz_0[i] = 3.0 * g_0_xxyy_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxzz_0[i] * pb_y +
                                  g_0_xxyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxyyy_0[i] = g_0_yyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyy_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xxyyy_0[i] * pb_x + g_0_xyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyz_0[i] = g_0_yyyy_0_xxyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xxyyz_0[i] * pb_x + g_0_xyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyzz_0[i] = g_0_yyyy_0_xxyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xxyzz_0[i] * pb_x + g_0_xyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxzzz_0[i] = 3.0 * g_0_xxyy_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxzzz_0[i] * pb_y +
                                  g_0_xxyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xyyyy_0[i] = g_0_yyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyy_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xyyyy_0[i] * pb_x + g_0_xyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyz_0[i] = g_0_yyyy_0_xyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xyyyz_0[i] * pb_x + g_0_xyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyzz_0[i] = g_0_yyyy_0_xyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xyyzz_0[i] * pb_x + g_0_xyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyzzz_0[i] = g_0_yyyy_0_xyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xyyyy_0_xyzzz_0[i] * pb_x + g_0_xyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xzzzz_0[i] = 3.0 * g_0_xxyy_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xzzzz_0[i] * pb_y +
                                  g_0_xxyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yyyyy_0[i] =
            g_0_yyyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyy_0[i] * pb_x + g_0_xyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyz_0[i] =
            g_0_yyyy_0_yyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyz_0[i] * pb_x + g_0_xyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyzz_0[i] =
            g_0_yyyy_0_yyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzz_0[i] * pb_x + g_0_xyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyzzz_0[i] =
            g_0_yyyy_0_yyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzz_0[i] * pb_x + g_0_xyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yzzzz_0[i] =
            g_0_yyyy_0_yzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzz_0[i] * pb_x + g_0_xyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zzzzz_0[i] =
            g_0_yyyy_0_zzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_zzzzz_0[i] * pb_x + g_0_xyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 231-252 components of targeted buffer : SISH

    auto g_0_xxyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 231);

    auto g_0_xxyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 232);

    auto g_0_xxyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 233);

    auto g_0_xxyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 234);

    auto g_0_xxyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 235);

    auto g_0_xxyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 236);

    auto g_0_xxyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 237);

    auto g_0_xxyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 238);

    auto g_0_xxyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 239);

    auto g_0_xxyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 240);

    auto g_0_xxyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 241);

    auto g_0_xxyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 242);

    auto g_0_xxyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 243);

    auto g_0_xxyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 244);

    auto g_0_xxyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 245);

    auto g_0_xxyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 246);

    auto g_0_xxyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 247);

    auto g_0_xxyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 248);

    auto g_0_xxyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 249);

    auto g_0_xxyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 250);

    auto g_0_xxyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 251);

#pragma omp simd aligned(g_0_xxyyy_0_xxxx_1,       \
                             g_0_xxyyy_0_xxxxx_0,  \
                             g_0_xxyyy_0_xxxxx_1,  \
                             g_0_xxyyy_0_xxxxy_0,  \
                             g_0_xxyyy_0_xxxxy_1,  \
                             g_0_xxyyy_0_xxxxz_0,  \
                             g_0_xxyyy_0_xxxxz_1,  \
                             g_0_xxyyy_0_xxxy_1,   \
                             g_0_xxyyy_0_xxxyy_0,  \
                             g_0_xxyyy_0_xxxyy_1,  \
                             g_0_xxyyy_0_xxxyz_0,  \
                             g_0_xxyyy_0_xxxyz_1,  \
                             g_0_xxyyy_0_xxxz_1,   \
                             g_0_xxyyy_0_xxxzz_0,  \
                             g_0_xxyyy_0_xxxzz_1,  \
                             g_0_xxyyy_0_xxyy_1,   \
                             g_0_xxyyy_0_xxyyy_0,  \
                             g_0_xxyyy_0_xxyyy_1,  \
                             g_0_xxyyy_0_xxyyz_0,  \
                             g_0_xxyyy_0_xxyyz_1,  \
                             g_0_xxyyy_0_xxyz_1,   \
                             g_0_xxyyy_0_xxyzz_0,  \
                             g_0_xxyyy_0_xxyzz_1,  \
                             g_0_xxyyy_0_xxzz_1,   \
                             g_0_xxyyy_0_xxzzz_0,  \
                             g_0_xxyyy_0_xxzzz_1,  \
                             g_0_xxyyy_0_xyyy_1,   \
                             g_0_xxyyy_0_xyyyy_0,  \
                             g_0_xxyyy_0_xyyyy_1,  \
                             g_0_xxyyy_0_xyyyz_0,  \
                             g_0_xxyyy_0_xyyyz_1,  \
                             g_0_xxyyy_0_xyyz_1,   \
                             g_0_xxyyy_0_xyyzz_0,  \
                             g_0_xxyyy_0_xyyzz_1,  \
                             g_0_xxyyy_0_xyzz_1,   \
                             g_0_xxyyy_0_xyzzz_0,  \
                             g_0_xxyyy_0_xyzzz_1,  \
                             g_0_xxyyy_0_xzzz_1,   \
                             g_0_xxyyy_0_xzzzz_0,  \
                             g_0_xxyyy_0_xzzzz_1,  \
                             g_0_xxyyy_0_yyyy_1,   \
                             g_0_xxyyy_0_yyyyy_0,  \
                             g_0_xxyyy_0_yyyyy_1,  \
                             g_0_xxyyy_0_yyyyz_0,  \
                             g_0_xxyyy_0_yyyyz_1,  \
                             g_0_xxyyy_0_yyyz_1,   \
                             g_0_xxyyy_0_yyyzz_0,  \
                             g_0_xxyyy_0_yyyzz_1,  \
                             g_0_xxyyy_0_yyzz_1,   \
                             g_0_xxyyy_0_yyzzz_0,  \
                             g_0_xxyyy_0_yyzzz_1,  \
                             g_0_xxyyy_0_yzzz_1,   \
                             g_0_xxyyy_0_yzzzz_0,  \
                             g_0_xxyyy_0_yzzzz_1,  \
                             g_0_xxyyy_0_zzzz_1,   \
                             g_0_xxyyy_0_zzzzz_0,  \
                             g_0_xxyyy_0_zzzzz_1,  \
                             g_0_xxyyyz_0_xxxxx_0, \
                             g_0_xxyyyz_0_xxxxy_0, \
                             g_0_xxyyyz_0_xxxxz_0, \
                             g_0_xxyyyz_0_xxxyy_0, \
                             g_0_xxyyyz_0_xxxyz_0, \
                             g_0_xxyyyz_0_xxxzz_0, \
                             g_0_xxyyyz_0_xxyyy_0, \
                             g_0_xxyyyz_0_xxyyz_0, \
                             g_0_xxyyyz_0_xxyzz_0, \
                             g_0_xxyyyz_0_xxzzz_0, \
                             g_0_xxyyyz_0_xyyyy_0, \
                             g_0_xxyyyz_0_xyyyz_0, \
                             g_0_xxyyyz_0_xyyzz_0, \
                             g_0_xxyyyz_0_xyzzz_0, \
                             g_0_xxyyyz_0_xzzzz_0, \
                             g_0_xxyyyz_0_yyyyy_0, \
                             g_0_xxyyyz_0_yyyyz_0, \
                             g_0_xxyyyz_0_yyyzz_0, \
                             g_0_xxyyyz_0_yyzzz_0, \
                             g_0_xxyyyz_0_yzzzz_0, \
                             g_0_xxyyyz_0_zzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xxxxx_0[i] = g_0_xxyyy_0_xxxxx_0[i] * pb_z + g_0_xxyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxy_0[i] = g_0_xxyyy_0_xxxxy_0[i] * pb_z + g_0_xxyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxz_0[i] = g_0_xxyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxz_0[i] * pb_z + g_0_xxyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyy_0[i] = g_0_xxyyy_0_xxxyy_0[i] * pb_z + g_0_xxyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyz_0[i] = g_0_xxyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyz_0[i] * pb_z + g_0_xxyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxzz_0[i] = 2.0 * g_0_xxyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxzz_0[i] * pb_z + g_0_xxyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyy_0[i] = g_0_xxyyy_0_xxyyy_0[i] * pb_z + g_0_xxyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyz_0[i] = g_0_xxyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyz_0[i] * pb_z + g_0_xxyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyzz_0[i] = 2.0 * g_0_xxyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzz_0[i] * pb_z + g_0_xxyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxzzz_0[i] = 3.0 * g_0_xxyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxzzz_0[i] * pb_z + g_0_xxyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyy_0[i] = g_0_xxyyy_0_xyyyy_0[i] * pb_z + g_0_xxyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyz_0[i] = g_0_xxyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyz_0[i] * pb_z + g_0_xxyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyzz_0[i] = 2.0 * g_0_xxyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzz_0[i] * pb_z + g_0_xxyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyzzz_0[i] = 3.0 * g_0_xxyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzz_0[i] * pb_z + g_0_xxyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xzzzz_0[i] = 4.0 * g_0_xxyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xzzzz_0[i] * pb_z + g_0_xxyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyy_0[i] = g_0_xxyyy_0_yyyyy_0[i] * pb_z + g_0_xxyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyz_0[i] = g_0_xxyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyz_0[i] * pb_z + g_0_xxyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyzz_0[i] = 2.0 * g_0_xxyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyzz_0[i] * pb_z + g_0_xxyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyzzz_0[i] = 3.0 * g_0_xxyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyzzz_0[i] * pb_z + g_0_xxyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yzzzz_0[i] = 4.0 * g_0_xxyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yzzzz_0[i] * pb_z + g_0_xxyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zzzzz_0[i] = 5.0 * g_0_xxyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_zzzzz_0[i] * pb_z + g_0_xxyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 252-273 components of targeted buffer : SISH

    auto g_0_xxyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 252);

    auto g_0_xxyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 253);

    auto g_0_xxyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 254);

    auto g_0_xxyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 255);

    auto g_0_xxyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 256);

    auto g_0_xxyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 257);

    auto g_0_xxyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 258);

    auto g_0_xxyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 259);

    auto g_0_xxyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 260);

    auto g_0_xxyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 261);

    auto g_0_xxyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 262);

    auto g_0_xxyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 263);

    auto g_0_xxyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 264);

    auto g_0_xxyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 265);

    auto g_0_xxyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 266);

    auto g_0_xxyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 267);

    auto g_0_xxyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 268);

    auto g_0_xxyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 269);

    auto g_0_xxyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 270);

    auto g_0_xxyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 271);

    auto g_0_xxyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 272);

#pragma omp simd aligned(g_0_xxyy_0_xxxxy_0,       \
                             g_0_xxyy_0_xxxxy_1,   \
                             g_0_xxyy_0_xxxyy_0,   \
                             g_0_xxyy_0_xxxyy_1,   \
                             g_0_xxyy_0_xxyyy_0,   \
                             g_0_xxyy_0_xxyyy_1,   \
                             g_0_xxyy_0_xyyyy_0,   \
                             g_0_xxyy_0_xyyyy_1,   \
                             g_0_xxyyz_0_xxxxy_0,  \
                             g_0_xxyyz_0_xxxxy_1,  \
                             g_0_xxyyz_0_xxxyy_0,  \
                             g_0_xxyyz_0_xxxyy_1,  \
                             g_0_xxyyz_0_xxyyy_0,  \
                             g_0_xxyyz_0_xxyyy_1,  \
                             g_0_xxyyz_0_xyyyy_0,  \
                             g_0_xxyyz_0_xyyyy_1,  \
                             g_0_xxyyzz_0_xxxxx_0, \
                             g_0_xxyyzz_0_xxxxy_0, \
                             g_0_xxyyzz_0_xxxxz_0, \
                             g_0_xxyyzz_0_xxxyy_0, \
                             g_0_xxyyzz_0_xxxyz_0, \
                             g_0_xxyyzz_0_xxxzz_0, \
                             g_0_xxyyzz_0_xxyyy_0, \
                             g_0_xxyyzz_0_xxyyz_0, \
                             g_0_xxyyzz_0_xxyzz_0, \
                             g_0_xxyyzz_0_xxzzz_0, \
                             g_0_xxyyzz_0_xyyyy_0, \
                             g_0_xxyyzz_0_xyyyz_0, \
                             g_0_xxyyzz_0_xyyzz_0, \
                             g_0_xxyyzz_0_xyzzz_0, \
                             g_0_xxyyzz_0_xzzzz_0, \
                             g_0_xxyyzz_0_yyyyy_0, \
                             g_0_xxyyzz_0_yyyyz_0, \
                             g_0_xxyyzz_0_yyyzz_0, \
                             g_0_xxyyzz_0_yyzzz_0, \
                             g_0_xxyyzz_0_yzzzz_0, \
                             g_0_xxyyzz_0_zzzzz_0, \
                             g_0_xxyzz_0_xxxxx_0,  \
                             g_0_xxyzz_0_xxxxx_1,  \
                             g_0_xxyzz_0_xxxxz_0,  \
                             g_0_xxyzz_0_xxxxz_1,  \
                             g_0_xxyzz_0_xxxzz_0,  \
                             g_0_xxyzz_0_xxxzz_1,  \
                             g_0_xxyzz_0_xxzzz_0,  \
                             g_0_xxyzz_0_xxzzz_1,  \
                             g_0_xxyzz_0_xzzzz_0,  \
                             g_0_xxyzz_0_xzzzz_1,  \
                             g_0_xxzz_0_xxxxx_0,   \
                             g_0_xxzz_0_xxxxx_1,   \
                             g_0_xxzz_0_xxxxz_0,   \
                             g_0_xxzz_0_xxxxz_1,   \
                             g_0_xxzz_0_xxxzz_0,   \
                             g_0_xxzz_0_xxxzz_1,   \
                             g_0_xxzz_0_xxzzz_0,   \
                             g_0_xxzz_0_xxzzz_1,   \
                             g_0_xxzz_0_xzzzz_0,   \
                             g_0_xxzz_0_xzzzz_1,   \
                             g_0_xyyzz_0_xxxyz_0,  \
                             g_0_xyyzz_0_xxxyz_1,  \
                             g_0_xyyzz_0_xxyyz_0,  \
                             g_0_xyyzz_0_xxyyz_1,  \
                             g_0_xyyzz_0_xxyz_1,   \
                             g_0_xyyzz_0_xxyzz_0,  \
                             g_0_xyyzz_0_xxyzz_1,  \
                             g_0_xyyzz_0_xyyyz_0,  \
                             g_0_xyyzz_0_xyyyz_1,  \
                             g_0_xyyzz_0_xyyz_1,   \
                             g_0_xyyzz_0_xyyzz_0,  \
                             g_0_xyyzz_0_xyyzz_1,  \
                             g_0_xyyzz_0_xyzz_1,   \
                             g_0_xyyzz_0_xyzzz_0,  \
                             g_0_xyyzz_0_xyzzz_1,  \
                             g_0_xyyzz_0_yyyyy_0,  \
                             g_0_xyyzz_0_yyyyy_1,  \
                             g_0_xyyzz_0_yyyyz_0,  \
                             g_0_xyyzz_0_yyyyz_1,  \
                             g_0_xyyzz_0_yyyz_1,   \
                             g_0_xyyzz_0_yyyzz_0,  \
                             g_0_xyyzz_0_yyyzz_1,  \
                             g_0_xyyzz_0_yyzz_1,   \
                             g_0_xyyzz_0_yyzzz_0,  \
                             g_0_xyyzz_0_yyzzz_1,  \
                             g_0_xyyzz_0_yzzz_1,   \
                             g_0_xyyzz_0_yzzzz_0,  \
                             g_0_xyyzz_0_yzzzz_1,  \
                             g_0_xyyzz_0_zzzzz_0,  \
                             g_0_xyyzz_0_zzzzz_1,  \
                             g_0_yyzz_0_xxxyz_0,   \
                             g_0_yyzz_0_xxxyz_1,   \
                             g_0_yyzz_0_xxyyz_0,   \
                             g_0_yyzz_0_xxyyz_1,   \
                             g_0_yyzz_0_xxyzz_0,   \
                             g_0_yyzz_0_xxyzz_1,   \
                             g_0_yyzz_0_xyyyz_0,   \
                             g_0_yyzz_0_xyyyz_1,   \
                             g_0_yyzz_0_xyyzz_0,   \
                             g_0_yyzz_0_xyyzz_1,   \
                             g_0_yyzz_0_xyzzz_0,   \
                             g_0_yyzz_0_xyzzz_1,   \
                             g_0_yyzz_0_yyyyy_0,   \
                             g_0_yyzz_0_yyyyy_1,   \
                             g_0_yyzz_0_yyyyz_0,   \
                             g_0_yyzz_0_yyyyz_1,   \
                             g_0_yyzz_0_yyyzz_0,   \
                             g_0_yyzz_0_yyyzz_1,   \
                             g_0_yyzz_0_yyzzz_0,   \
                             g_0_yyzz_0_yyzzz_1,   \
                             g_0_yyzz_0_yzzzz_0,   \
                             g_0_yyzz_0_yzzzz_1,   \
                             g_0_yyzz_0_zzzzz_0,   \
                             g_0_yyzz_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xxxxx_0[i] =
            g_0_xxzz_0_xxxxx_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxx_0[i] * pb_y + g_0_xxyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxy_0[i] =
            g_0_xxyy_0_xxxxy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxy_0[i] * pb_z + g_0_xxyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxz_0[i] =
            g_0_xxzz_0_xxxxz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxz_0[i] * pb_y + g_0_xxyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxyy_0[i] =
            g_0_xxyy_0_xxxyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxyy_0[i] * pb_z + g_0_xxyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxyz_0[i] = g_0_yyzz_0_xxxyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzz_0_xxyz_1[i] * fi_abcd_0 +
                                  g_0_xyyzz_0_xxxyz_0[i] * pb_x + g_0_xyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxzz_0[i] =
            g_0_xxzz_0_xxxzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxzz_0[i] * pb_y + g_0_xxyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxyyy_0[i] =
            g_0_xxyy_0_xxyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxyyy_0[i] * pb_z + g_0_xxyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxyyz_0[i] = g_0_yyzz_0_xxyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyzz_0_xxyyz_0[i] * pb_x + g_0_xyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyzz_0[i] = g_0_yyzz_0_xxyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyzz_0_xxyzz_0[i] * pb_x + g_0_xyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxzzz_0[i] =
            g_0_xxzz_0_xxzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxzzz_0[i] * pb_y + g_0_xxyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xyyyy_0[i] =
            g_0_xxyy_0_xyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xyyyy_0[i] * pb_z + g_0_xxyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xyyyz_0[i] = g_0_yyzz_0_xyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyzz_0_xyyyz_0[i] * pb_x + g_0_xyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyzz_0[i] = g_0_yyzz_0_xyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyzz_0_xyyzz_0[i] * pb_x + g_0_xyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyzzz_0[i] = g_0_yyzz_0_xyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xyyzz_0_xyzzz_0[i] * pb_x + g_0_xyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xzzzz_0[i] =
            g_0_xxzz_0_xzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xzzzz_0[i] * pb_y + g_0_xxyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yyyyy_0[i] =
            g_0_yyzz_0_yyyyy_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyy_0[i] * pb_x + g_0_xyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyz_0[i] =
            g_0_yyzz_0_yyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyz_0[i] * pb_x + g_0_xyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyzz_0[i] =
            g_0_yyzz_0_yyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzz_0[i] * pb_x + g_0_xyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyzzz_0[i] =
            g_0_yyzz_0_yyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzz_0[i] * pb_x + g_0_xyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yzzzz_0[i] =
            g_0_yyzz_0_yzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzz_0[i] * pb_x + g_0_xyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zzzzz_0[i] =
            g_0_yyzz_0_zzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_zzzzz_0[i] * pb_x + g_0_xyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 273-294 components of targeted buffer : SISH

    auto g_0_xxyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 273);

    auto g_0_xxyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 274);

    auto g_0_xxyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 275);

    auto g_0_xxyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 276);

    auto g_0_xxyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 277);

    auto g_0_xxyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 278);

    auto g_0_xxyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 279);

    auto g_0_xxyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 280);

    auto g_0_xxyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 281);

    auto g_0_xxyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 282);

    auto g_0_xxyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 283);

    auto g_0_xxyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 284);

    auto g_0_xxyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 285);

    auto g_0_xxyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 286);

    auto g_0_xxyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 287);

    auto g_0_xxyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 288);

    auto g_0_xxyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 289);

    auto g_0_xxyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 290);

    auto g_0_xxyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 291);

    auto g_0_xxyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 292);

    auto g_0_xxyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 293);

#pragma omp simd aligned(g_0_xxyzzz_0_xxxxx_0,     \
                             g_0_xxyzzz_0_xxxxy_0, \
                             g_0_xxyzzz_0_xxxxz_0, \
                             g_0_xxyzzz_0_xxxyy_0, \
                             g_0_xxyzzz_0_xxxyz_0, \
                             g_0_xxyzzz_0_xxxzz_0, \
                             g_0_xxyzzz_0_xxyyy_0, \
                             g_0_xxyzzz_0_xxyyz_0, \
                             g_0_xxyzzz_0_xxyzz_0, \
                             g_0_xxyzzz_0_xxzzz_0, \
                             g_0_xxyzzz_0_xyyyy_0, \
                             g_0_xxyzzz_0_xyyyz_0, \
                             g_0_xxyzzz_0_xyyzz_0, \
                             g_0_xxyzzz_0_xyzzz_0, \
                             g_0_xxyzzz_0_xzzzz_0, \
                             g_0_xxyzzz_0_yyyyy_0, \
                             g_0_xxyzzz_0_yyyyz_0, \
                             g_0_xxyzzz_0_yyyzz_0, \
                             g_0_xxyzzz_0_yyzzz_0, \
                             g_0_xxyzzz_0_yzzzz_0, \
                             g_0_xxyzzz_0_zzzzz_0, \
                             g_0_xxzzz_0_xxxx_1,   \
                             g_0_xxzzz_0_xxxxx_0,  \
                             g_0_xxzzz_0_xxxxx_1,  \
                             g_0_xxzzz_0_xxxxy_0,  \
                             g_0_xxzzz_0_xxxxy_1,  \
                             g_0_xxzzz_0_xxxxz_0,  \
                             g_0_xxzzz_0_xxxxz_1,  \
                             g_0_xxzzz_0_xxxy_1,   \
                             g_0_xxzzz_0_xxxyy_0,  \
                             g_0_xxzzz_0_xxxyy_1,  \
                             g_0_xxzzz_0_xxxyz_0,  \
                             g_0_xxzzz_0_xxxyz_1,  \
                             g_0_xxzzz_0_xxxz_1,   \
                             g_0_xxzzz_0_xxxzz_0,  \
                             g_0_xxzzz_0_xxxzz_1,  \
                             g_0_xxzzz_0_xxyy_1,   \
                             g_0_xxzzz_0_xxyyy_0,  \
                             g_0_xxzzz_0_xxyyy_1,  \
                             g_0_xxzzz_0_xxyyz_0,  \
                             g_0_xxzzz_0_xxyyz_1,  \
                             g_0_xxzzz_0_xxyz_1,   \
                             g_0_xxzzz_0_xxyzz_0,  \
                             g_0_xxzzz_0_xxyzz_1,  \
                             g_0_xxzzz_0_xxzz_1,   \
                             g_0_xxzzz_0_xxzzz_0,  \
                             g_0_xxzzz_0_xxzzz_1,  \
                             g_0_xxzzz_0_xyyy_1,   \
                             g_0_xxzzz_0_xyyyy_0,  \
                             g_0_xxzzz_0_xyyyy_1,  \
                             g_0_xxzzz_0_xyyyz_0,  \
                             g_0_xxzzz_0_xyyyz_1,  \
                             g_0_xxzzz_0_xyyz_1,   \
                             g_0_xxzzz_0_xyyzz_0,  \
                             g_0_xxzzz_0_xyyzz_1,  \
                             g_0_xxzzz_0_xyzz_1,   \
                             g_0_xxzzz_0_xyzzz_0,  \
                             g_0_xxzzz_0_xyzzz_1,  \
                             g_0_xxzzz_0_xzzz_1,   \
                             g_0_xxzzz_0_xzzzz_0,  \
                             g_0_xxzzz_0_xzzzz_1,  \
                             g_0_xxzzz_0_yyyy_1,   \
                             g_0_xxzzz_0_yyyyy_0,  \
                             g_0_xxzzz_0_yyyyy_1,  \
                             g_0_xxzzz_0_yyyyz_0,  \
                             g_0_xxzzz_0_yyyyz_1,  \
                             g_0_xxzzz_0_yyyz_1,   \
                             g_0_xxzzz_0_yyyzz_0,  \
                             g_0_xxzzz_0_yyyzz_1,  \
                             g_0_xxzzz_0_yyzz_1,   \
                             g_0_xxzzz_0_yyzzz_0,  \
                             g_0_xxzzz_0_yyzzz_1,  \
                             g_0_xxzzz_0_yzzz_1,   \
                             g_0_xxzzz_0_yzzzz_0,  \
                             g_0_xxzzz_0_yzzzz_1,  \
                             g_0_xxzzz_0_zzzz_1,   \
                             g_0_xxzzz_0_zzzzz_0,  \
                             g_0_xxzzz_0_zzzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xxxxx_0[i] = g_0_xxzzz_0_xxxxx_0[i] * pb_y + g_0_xxzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxy_0[i] = g_0_xxzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxy_0[i] * pb_y + g_0_xxzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxz_0[i] = g_0_xxzzz_0_xxxxz_0[i] * pb_y + g_0_xxzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyy_0[i] = 2.0 * g_0_xxzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyy_0[i] * pb_y + g_0_xxzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyz_0[i] = g_0_xxzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyz_0[i] * pb_y + g_0_xxzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxzz_0[i] = g_0_xxzzz_0_xxxzz_0[i] * pb_y + g_0_xxzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyy_0[i] = 3.0 * g_0_xxzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyy_0[i] * pb_y + g_0_xxzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyz_0[i] = 2.0 * g_0_xxzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyz_0[i] * pb_y + g_0_xxzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyzz_0[i] = g_0_xxzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzz_0[i] * pb_y + g_0_xxzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxzzz_0[i] = g_0_xxzzz_0_xxzzz_0[i] * pb_y + g_0_xxzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyy_0[i] = 4.0 * g_0_xxzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyy_0[i] * pb_y + g_0_xxzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyz_0[i] = 3.0 * g_0_xxzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyz_0[i] * pb_y + g_0_xxzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyzz_0[i] = 2.0 * g_0_xxzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzz_0[i] * pb_y + g_0_xxzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyzzz_0[i] = g_0_xxzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzz_0[i] * pb_y + g_0_xxzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xzzzz_0[i] = g_0_xxzzz_0_xzzzz_0[i] * pb_y + g_0_xxzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyy_0[i] = 5.0 * g_0_xxzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyy_0[i] * pb_y + g_0_xxzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyz_0[i] = 4.0 * g_0_xxzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyz_0[i] * pb_y + g_0_xxzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyzz_0[i] = 3.0 * g_0_xxzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyzz_0[i] * pb_y + g_0_xxzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyzzz_0[i] = 2.0 * g_0_xxzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyzzz_0[i] * pb_y + g_0_xxzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yzzzz_0[i] = g_0_xxzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yzzzz_0[i] * pb_y + g_0_xxzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zzzzz_0[i] = g_0_xxzzz_0_zzzzz_0[i] * pb_y + g_0_xxzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 294-315 components of targeted buffer : SISH

    auto g_0_xxzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 294);

    auto g_0_xxzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 295);

    auto g_0_xxzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 296);

    auto g_0_xxzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 297);

    auto g_0_xxzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 298);

    auto g_0_xxzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 299);

    auto g_0_xxzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 300);

    auto g_0_xxzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 301);

    auto g_0_xxzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 302);

    auto g_0_xxzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 303);

    auto g_0_xxzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 304);

    auto g_0_xxzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 305);

    auto g_0_xxzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 306);

    auto g_0_xxzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 307);

    auto g_0_xxzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 308);

    auto g_0_xxzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 309);

    auto g_0_xxzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 310);

    auto g_0_xxzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 311);

    auto g_0_xxzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 312);

    auto g_0_xxzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 313);

    auto g_0_xxzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 314);

#pragma omp simd aligned(g_0_xxzz_0_xxxxx_0,       \
                             g_0_xxzz_0_xxxxx_1,   \
                             g_0_xxzz_0_xxxxy_0,   \
                             g_0_xxzz_0_xxxxy_1,   \
                             g_0_xxzz_0_xxxyy_0,   \
                             g_0_xxzz_0_xxxyy_1,   \
                             g_0_xxzz_0_xxyyy_0,   \
                             g_0_xxzz_0_xxyyy_1,   \
                             g_0_xxzz_0_xyyyy_0,   \
                             g_0_xxzz_0_xyyyy_1,   \
                             g_0_xxzzz_0_xxxxx_0,  \
                             g_0_xxzzz_0_xxxxx_1,  \
                             g_0_xxzzz_0_xxxxy_0,  \
                             g_0_xxzzz_0_xxxxy_1,  \
                             g_0_xxzzz_0_xxxyy_0,  \
                             g_0_xxzzz_0_xxxyy_1,  \
                             g_0_xxzzz_0_xxyyy_0,  \
                             g_0_xxzzz_0_xxyyy_1,  \
                             g_0_xxzzz_0_xyyyy_0,  \
                             g_0_xxzzz_0_xyyyy_1,  \
                             g_0_xxzzzz_0_xxxxx_0, \
                             g_0_xxzzzz_0_xxxxy_0, \
                             g_0_xxzzzz_0_xxxxz_0, \
                             g_0_xxzzzz_0_xxxyy_0, \
                             g_0_xxzzzz_0_xxxyz_0, \
                             g_0_xxzzzz_0_xxxzz_0, \
                             g_0_xxzzzz_0_xxyyy_0, \
                             g_0_xxzzzz_0_xxyyz_0, \
                             g_0_xxzzzz_0_xxyzz_0, \
                             g_0_xxzzzz_0_xxzzz_0, \
                             g_0_xxzzzz_0_xyyyy_0, \
                             g_0_xxzzzz_0_xyyyz_0, \
                             g_0_xxzzzz_0_xyyzz_0, \
                             g_0_xxzzzz_0_xyzzz_0, \
                             g_0_xxzzzz_0_xzzzz_0, \
                             g_0_xxzzzz_0_yyyyy_0, \
                             g_0_xxzzzz_0_yyyyz_0, \
                             g_0_xxzzzz_0_yyyzz_0, \
                             g_0_xxzzzz_0_yyzzz_0, \
                             g_0_xxzzzz_0_yzzzz_0, \
                             g_0_xxzzzz_0_zzzzz_0, \
                             g_0_xzzzz_0_xxxxz_0,  \
                             g_0_xzzzz_0_xxxxz_1,  \
                             g_0_xzzzz_0_xxxyz_0,  \
                             g_0_xzzzz_0_xxxyz_1,  \
                             g_0_xzzzz_0_xxxz_1,   \
                             g_0_xzzzz_0_xxxzz_0,  \
                             g_0_xzzzz_0_xxxzz_1,  \
                             g_0_xzzzz_0_xxyyz_0,  \
                             g_0_xzzzz_0_xxyyz_1,  \
                             g_0_xzzzz_0_xxyz_1,   \
                             g_0_xzzzz_0_xxyzz_0,  \
                             g_0_xzzzz_0_xxyzz_1,  \
                             g_0_xzzzz_0_xxzz_1,   \
                             g_0_xzzzz_0_xxzzz_0,  \
                             g_0_xzzzz_0_xxzzz_1,  \
                             g_0_xzzzz_0_xyyyz_0,  \
                             g_0_xzzzz_0_xyyyz_1,  \
                             g_0_xzzzz_0_xyyz_1,   \
                             g_0_xzzzz_0_xyyzz_0,  \
                             g_0_xzzzz_0_xyyzz_1,  \
                             g_0_xzzzz_0_xyzz_1,   \
                             g_0_xzzzz_0_xyzzz_0,  \
                             g_0_xzzzz_0_xyzzz_1,  \
                             g_0_xzzzz_0_xzzz_1,   \
                             g_0_xzzzz_0_xzzzz_0,  \
                             g_0_xzzzz_0_xzzzz_1,  \
                             g_0_xzzzz_0_yyyyy_0,  \
                             g_0_xzzzz_0_yyyyy_1,  \
                             g_0_xzzzz_0_yyyyz_0,  \
                             g_0_xzzzz_0_yyyyz_1,  \
                             g_0_xzzzz_0_yyyz_1,   \
                             g_0_xzzzz_0_yyyzz_0,  \
                             g_0_xzzzz_0_yyyzz_1,  \
                             g_0_xzzzz_0_yyzz_1,   \
                             g_0_xzzzz_0_yyzzz_0,  \
                             g_0_xzzzz_0_yyzzz_1,  \
                             g_0_xzzzz_0_yzzz_1,   \
                             g_0_xzzzz_0_yzzzz_0,  \
                             g_0_xzzzz_0_yzzzz_1,  \
                             g_0_xzzzz_0_zzzz_1,   \
                             g_0_xzzzz_0_zzzzz_0,  \
                             g_0_xzzzz_0_zzzzz_1,  \
                             g_0_zzzz_0_xxxxz_0,   \
                             g_0_zzzz_0_xxxxz_1,   \
                             g_0_zzzz_0_xxxyz_0,   \
                             g_0_zzzz_0_xxxyz_1,   \
                             g_0_zzzz_0_xxxzz_0,   \
                             g_0_zzzz_0_xxxzz_1,   \
                             g_0_zzzz_0_xxyyz_0,   \
                             g_0_zzzz_0_xxyyz_1,   \
                             g_0_zzzz_0_xxyzz_0,   \
                             g_0_zzzz_0_xxyzz_1,   \
                             g_0_zzzz_0_xxzzz_0,   \
                             g_0_zzzz_0_xxzzz_1,   \
                             g_0_zzzz_0_xyyyz_0,   \
                             g_0_zzzz_0_xyyyz_1,   \
                             g_0_zzzz_0_xyyzz_0,   \
                             g_0_zzzz_0_xyyzz_1,   \
                             g_0_zzzz_0_xyzzz_0,   \
                             g_0_zzzz_0_xyzzz_1,   \
                             g_0_zzzz_0_xzzzz_0,   \
                             g_0_zzzz_0_xzzzz_1,   \
                             g_0_zzzz_0_yyyyy_0,   \
                             g_0_zzzz_0_yyyyy_1,   \
                             g_0_zzzz_0_yyyyz_0,   \
                             g_0_zzzz_0_yyyyz_1,   \
                             g_0_zzzz_0_yyyzz_0,   \
                             g_0_zzzz_0_yyyzz_1,   \
                             g_0_zzzz_0_yyzzz_0,   \
                             g_0_zzzz_0_yyzzz_1,   \
                             g_0_zzzz_0_yzzzz_0,   \
                             g_0_zzzz_0_yzzzz_1,   \
                             g_0_zzzz_0_zzzzz_0,   \
                             g_0_zzzz_0_zzzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xxxxx_0[i] = 3.0 * g_0_xxzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxx_0[i] * pb_z +
                                  g_0_xxzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxy_0[i] = 3.0 * g_0_xxzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxy_0[i] * pb_z +
                                  g_0_xxzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxz_0[i] = g_0_zzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzz_0_xxxz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xxxxz_0[i] * pb_x + g_0_xzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyy_0[i] = 3.0 * g_0_xxzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxyy_0[i] * pb_z +
                                  g_0_xxzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxyz_0[i] = g_0_zzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxyz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xxxyz_0[i] * pb_x + g_0_xzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxzz_0[i] = g_0_zzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xxxzz_0[i] * pb_x + g_0_xzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxyyy_0[i] * pb_z +
                                  g_0_xxzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxyyz_0[i] = g_0_zzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyyz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xxyyz_0[i] * pb_x + g_0_xzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyzz_0[i] = g_0_zzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xxyzz_0[i] * pb_x + g_0_xzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxzzz_0[i] = g_0_zzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xxzzz_0[i] * pb_x + g_0_xzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyy_0[i] = 3.0 * g_0_xxzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xyyyy_0[i] * pb_z +
                                  g_0_xxzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xyyyz_0[i] = g_0_zzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xyyyz_0[i] * pb_x + g_0_xzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyzz_0[i] = g_0_zzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xyyzz_0[i] * pb_x + g_0_xzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyzzz_0[i] = g_0_zzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xyzzz_0[i] * pb_x + g_0_xzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xzzzz_0[i] = g_0_zzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzz_0_xzzzz_0[i] * pb_x + g_0_xzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyy_0[i] =
            g_0_zzzz_0_yyyyy_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyy_0[i] * pb_x + g_0_xzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyz_0[i] =
            g_0_zzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyz_0[i] * pb_x + g_0_xzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyzz_0[i] =
            g_0_zzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzz_0[i] * pb_x + g_0_xzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyzzz_0[i] =
            g_0_zzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzz_0[i] * pb_x + g_0_xzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yzzzz_0[i] =
            g_0_zzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzz_0[i] * pb_x + g_0_xzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zzzzz_0[i] =
            g_0_zzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzz_0[i] * pb_x + g_0_xzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-336 components of targeted buffer : SISH

    auto g_0_xyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 315);

    auto g_0_xyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 316);

    auto g_0_xyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 317);

    auto g_0_xyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 318);

    auto g_0_xyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 319);

    auto g_0_xyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 320);

    auto g_0_xyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 321);

    auto g_0_xyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 322);

    auto g_0_xyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 323);

    auto g_0_xyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 324);

    auto g_0_xyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 325);

    auto g_0_xyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 326);

    auto g_0_xyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 327);

    auto g_0_xyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 328);

    auto g_0_xyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 329);

    auto g_0_xyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 330);

    auto g_0_xyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 331);

    auto g_0_xyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 332);

    auto g_0_xyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 333);

    auto g_0_xyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 334);

    auto g_0_xyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 335);

#pragma omp simd aligned(g_0_xyyyyy_0_xxxxx_0,     \
                             g_0_xyyyyy_0_xxxxy_0, \
                             g_0_xyyyyy_0_xxxxz_0, \
                             g_0_xyyyyy_0_xxxyy_0, \
                             g_0_xyyyyy_0_xxxyz_0, \
                             g_0_xyyyyy_0_xxxzz_0, \
                             g_0_xyyyyy_0_xxyyy_0, \
                             g_0_xyyyyy_0_xxyyz_0, \
                             g_0_xyyyyy_0_xxyzz_0, \
                             g_0_xyyyyy_0_xxzzz_0, \
                             g_0_xyyyyy_0_xyyyy_0, \
                             g_0_xyyyyy_0_xyyyz_0, \
                             g_0_xyyyyy_0_xyyzz_0, \
                             g_0_xyyyyy_0_xyzzz_0, \
                             g_0_xyyyyy_0_xzzzz_0, \
                             g_0_xyyyyy_0_yyyyy_0, \
                             g_0_xyyyyy_0_yyyyz_0, \
                             g_0_xyyyyy_0_yyyzz_0, \
                             g_0_xyyyyy_0_yyzzz_0, \
                             g_0_xyyyyy_0_yzzzz_0, \
                             g_0_xyyyyy_0_zzzzz_0, \
                             g_0_yyyyy_0_xxxx_1,   \
                             g_0_yyyyy_0_xxxxx_0,  \
                             g_0_yyyyy_0_xxxxx_1,  \
                             g_0_yyyyy_0_xxxxy_0,  \
                             g_0_yyyyy_0_xxxxy_1,  \
                             g_0_yyyyy_0_xxxxz_0,  \
                             g_0_yyyyy_0_xxxxz_1,  \
                             g_0_yyyyy_0_xxxy_1,   \
                             g_0_yyyyy_0_xxxyy_0,  \
                             g_0_yyyyy_0_xxxyy_1,  \
                             g_0_yyyyy_0_xxxyz_0,  \
                             g_0_yyyyy_0_xxxyz_1,  \
                             g_0_yyyyy_0_xxxz_1,   \
                             g_0_yyyyy_0_xxxzz_0,  \
                             g_0_yyyyy_0_xxxzz_1,  \
                             g_0_yyyyy_0_xxyy_1,   \
                             g_0_yyyyy_0_xxyyy_0,  \
                             g_0_yyyyy_0_xxyyy_1,  \
                             g_0_yyyyy_0_xxyyz_0,  \
                             g_0_yyyyy_0_xxyyz_1,  \
                             g_0_yyyyy_0_xxyz_1,   \
                             g_0_yyyyy_0_xxyzz_0,  \
                             g_0_yyyyy_0_xxyzz_1,  \
                             g_0_yyyyy_0_xxzz_1,   \
                             g_0_yyyyy_0_xxzzz_0,  \
                             g_0_yyyyy_0_xxzzz_1,  \
                             g_0_yyyyy_0_xyyy_1,   \
                             g_0_yyyyy_0_xyyyy_0,  \
                             g_0_yyyyy_0_xyyyy_1,  \
                             g_0_yyyyy_0_xyyyz_0,  \
                             g_0_yyyyy_0_xyyyz_1,  \
                             g_0_yyyyy_0_xyyz_1,   \
                             g_0_yyyyy_0_xyyzz_0,  \
                             g_0_yyyyy_0_xyyzz_1,  \
                             g_0_yyyyy_0_xyzz_1,   \
                             g_0_yyyyy_0_xyzzz_0,  \
                             g_0_yyyyy_0_xyzzz_1,  \
                             g_0_yyyyy_0_xzzz_1,   \
                             g_0_yyyyy_0_xzzzz_0,  \
                             g_0_yyyyy_0_xzzzz_1,  \
                             g_0_yyyyy_0_yyyy_1,   \
                             g_0_yyyyy_0_yyyyy_0,  \
                             g_0_yyyyy_0_yyyyy_1,  \
                             g_0_yyyyy_0_yyyyz_0,  \
                             g_0_yyyyy_0_yyyyz_1,  \
                             g_0_yyyyy_0_yyyz_1,   \
                             g_0_yyyyy_0_yyyzz_0,  \
                             g_0_yyyyy_0_yyyzz_1,  \
                             g_0_yyyyy_0_yyzz_1,   \
                             g_0_yyyyy_0_yyzzz_0,  \
                             g_0_yyyyy_0_yyzzz_1,  \
                             g_0_yyyyy_0_yzzz_1,   \
                             g_0_yyyyy_0_yzzzz_0,  \
                             g_0_yyyyy_0_yzzzz_1,  \
                             g_0_yyyyy_0_zzzz_1,   \
                             g_0_yyyyy_0_zzzzz_0,  \
                             g_0_yyyyy_0_zzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xxxxx_0[i] = 5.0 * g_0_yyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxx_0[i] * pb_x + g_0_yyyyy_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxy_0[i] = 4.0 * g_0_yyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxy_0[i] * pb_x + g_0_yyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxz_0[i] = 4.0 * g_0_yyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxz_0[i] * pb_x + g_0_yyyyy_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyy_0[i] = 3.0 * g_0_yyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyy_0[i] * pb_x + g_0_yyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyz_0[i] = 3.0 * g_0_yyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyz_0[i] * pb_x + g_0_yyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxzz_0[i] = 3.0 * g_0_yyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzz_0[i] * pb_x + g_0_yyyyy_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyy_0[i] = 2.0 * g_0_yyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyy_0[i] * pb_x + g_0_yyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyz_0[i] = 2.0 * g_0_yyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyz_0[i] * pb_x + g_0_yyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyzz_0[i] = 2.0 * g_0_yyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzz_0[i] * pb_x + g_0_yyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxzzz_0[i] = 2.0 * g_0_yyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzz_0[i] * pb_x + g_0_yyyyy_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyy_0[i] = g_0_yyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyy_0[i] * pb_x + g_0_yyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyz_0[i] = g_0_yyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyz_0[i] * pb_x + g_0_yyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyzz_0[i] = g_0_yyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzz_0[i] * pb_x + g_0_yyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyzzz_0[i] = g_0_yyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzz_0[i] * pb_x + g_0_yyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xzzzz_0[i] = g_0_yyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzz_0[i] * pb_x + g_0_yyyyy_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyy_0[i] = g_0_yyyyy_0_yyyyy_0[i] * pb_x + g_0_yyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyz_0[i] = g_0_yyyyy_0_yyyyz_0[i] * pb_x + g_0_yyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyzz_0[i] = g_0_yyyyy_0_yyyzz_0[i] * pb_x + g_0_yyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyzzz_0[i] = g_0_yyyyy_0_yyzzz_0[i] * pb_x + g_0_yyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yzzzz_0[i] = g_0_yyyyy_0_yzzzz_0[i] * pb_x + g_0_yyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zzzzz_0[i] = g_0_yyyyy_0_zzzzz_0[i] * pb_x + g_0_yyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 336-357 components of targeted buffer : SISH

    auto g_0_xyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 336);

    auto g_0_xyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 337);

    auto g_0_xyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 338);

    auto g_0_xyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 339);

    auto g_0_xyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 340);

    auto g_0_xyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 341);

    auto g_0_xyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 342);

    auto g_0_xyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 343);

    auto g_0_xyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 344);

    auto g_0_xyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 345);

    auto g_0_xyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 346);

    auto g_0_xyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 347);

    auto g_0_xyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 348);

    auto g_0_xyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 349);

    auto g_0_xyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 350);

    auto g_0_xyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 351);

    auto g_0_xyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 352);

    auto g_0_xyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 353);

    auto g_0_xyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 354);

    auto g_0_xyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 355);

    auto g_0_xyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 356);

#pragma omp simd aligned(g_0_xyyyy_0_xxxxx_0,      \
                             g_0_xyyyy_0_xxxxx_1,  \
                             g_0_xyyyy_0_xxxxy_0,  \
                             g_0_xyyyy_0_xxxxy_1,  \
                             g_0_xyyyy_0_xxxyy_0,  \
                             g_0_xyyyy_0_xxxyy_1,  \
                             g_0_xyyyy_0_xxyyy_0,  \
                             g_0_xyyyy_0_xxyyy_1,  \
                             g_0_xyyyy_0_xyyyy_0,  \
                             g_0_xyyyy_0_xyyyy_1,  \
                             g_0_xyyyyz_0_xxxxx_0, \
                             g_0_xyyyyz_0_xxxxy_0, \
                             g_0_xyyyyz_0_xxxxz_0, \
                             g_0_xyyyyz_0_xxxyy_0, \
                             g_0_xyyyyz_0_xxxyz_0, \
                             g_0_xyyyyz_0_xxxzz_0, \
                             g_0_xyyyyz_0_xxyyy_0, \
                             g_0_xyyyyz_0_xxyyz_0, \
                             g_0_xyyyyz_0_xxyzz_0, \
                             g_0_xyyyyz_0_xxzzz_0, \
                             g_0_xyyyyz_0_xyyyy_0, \
                             g_0_xyyyyz_0_xyyyz_0, \
                             g_0_xyyyyz_0_xyyzz_0, \
                             g_0_xyyyyz_0_xyzzz_0, \
                             g_0_xyyyyz_0_xzzzz_0, \
                             g_0_xyyyyz_0_yyyyy_0, \
                             g_0_xyyyyz_0_yyyyz_0, \
                             g_0_xyyyyz_0_yyyzz_0, \
                             g_0_xyyyyz_0_yyzzz_0, \
                             g_0_xyyyyz_0_yzzzz_0, \
                             g_0_xyyyyz_0_zzzzz_0, \
                             g_0_yyyyz_0_xxxxz_0,  \
                             g_0_yyyyz_0_xxxxz_1,  \
                             g_0_yyyyz_0_xxxyz_0,  \
                             g_0_yyyyz_0_xxxyz_1,  \
                             g_0_yyyyz_0_xxxz_1,   \
                             g_0_yyyyz_0_xxxzz_0,  \
                             g_0_yyyyz_0_xxxzz_1,  \
                             g_0_yyyyz_0_xxyyz_0,  \
                             g_0_yyyyz_0_xxyyz_1,  \
                             g_0_yyyyz_0_xxyz_1,   \
                             g_0_yyyyz_0_xxyzz_0,  \
                             g_0_yyyyz_0_xxyzz_1,  \
                             g_0_yyyyz_0_xxzz_1,   \
                             g_0_yyyyz_0_xxzzz_0,  \
                             g_0_yyyyz_0_xxzzz_1,  \
                             g_0_yyyyz_0_xyyyz_0,  \
                             g_0_yyyyz_0_xyyyz_1,  \
                             g_0_yyyyz_0_xyyz_1,   \
                             g_0_yyyyz_0_xyyzz_0,  \
                             g_0_yyyyz_0_xyyzz_1,  \
                             g_0_yyyyz_0_xyzz_1,   \
                             g_0_yyyyz_0_xyzzz_0,  \
                             g_0_yyyyz_0_xyzzz_1,  \
                             g_0_yyyyz_0_xzzz_1,   \
                             g_0_yyyyz_0_xzzzz_0,  \
                             g_0_yyyyz_0_xzzzz_1,  \
                             g_0_yyyyz_0_yyyyy_0,  \
                             g_0_yyyyz_0_yyyyy_1,  \
                             g_0_yyyyz_0_yyyyz_0,  \
                             g_0_yyyyz_0_yyyyz_1,  \
                             g_0_yyyyz_0_yyyz_1,   \
                             g_0_yyyyz_0_yyyzz_0,  \
                             g_0_yyyyz_0_yyyzz_1,  \
                             g_0_yyyyz_0_yyzz_1,   \
                             g_0_yyyyz_0_yyzzz_0,  \
                             g_0_yyyyz_0_yyzzz_1,  \
                             g_0_yyyyz_0_yzzz_1,   \
                             g_0_yyyyz_0_yzzzz_0,  \
                             g_0_yyyyz_0_yzzzz_1,  \
                             g_0_yyyyz_0_zzzz_1,   \
                             g_0_yyyyz_0_zzzzz_0,  \
                             g_0_yyyyz_0_zzzzz_1,  \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xxxxx_0[i] = g_0_xyyyy_0_xxxxx_0[i] * pb_z + g_0_xyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxy_0[i] = g_0_xyyyy_0_xxxxy_0[i] * pb_z + g_0_xyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxz_0[i] = 4.0 * g_0_yyyyz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxz_0[i] * pb_x + g_0_yyyyz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyy_0[i] = g_0_xyyyy_0_xxxyy_0[i] * pb_z + g_0_xyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxyz_0[i] = 3.0 * g_0_yyyyz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyz_0[i] * pb_x + g_0_yyyyz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxzz_0[i] = 3.0 * g_0_yyyyz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxzz_0[i] * pb_x + g_0_yyyyz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyy_0[i] = g_0_xyyyy_0_xxyyy_0[i] * pb_z + g_0_xyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxyyz_0[i] = 2.0 * g_0_yyyyz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyz_0[i] * pb_x + g_0_yyyyz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyyz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyzz_0[i] * pb_x + g_0_yyyyz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxzzz_0[i] = 2.0 * g_0_yyyyz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxzzz_0[i] * pb_x + g_0_yyyyz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyy_0[i] = g_0_xyyyy_0_xyyyy_0[i] * pb_z + g_0_xyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xyyyz_0[i] = g_0_yyyyz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyz_0[i] * pb_x + g_0_yyyyz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyzz_0[i] = g_0_yyyyz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyzz_0[i] * pb_x + g_0_yyyyz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyzzz_0[i] = g_0_yyyyz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyzzz_0[i] * pb_x + g_0_yyyyz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xzzzz_0[i] = g_0_yyyyz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xzzzz_0[i] * pb_x + g_0_yyyyz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyy_0[i] = g_0_yyyyz_0_yyyyy_0[i] * pb_x + g_0_yyyyz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyz_0[i] = g_0_yyyyz_0_yyyyz_0[i] * pb_x + g_0_yyyyz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyzz_0[i] = g_0_yyyyz_0_yyyzz_0[i] * pb_x + g_0_yyyyz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyzzz_0[i] = g_0_yyyyz_0_yyzzz_0[i] * pb_x + g_0_yyyyz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yzzzz_0[i] = g_0_yyyyz_0_yzzzz_0[i] * pb_x + g_0_yyyyz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zzzzz_0[i] = g_0_yyyyz_0_zzzzz_0[i] * pb_x + g_0_yyyyz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 357-378 components of targeted buffer : SISH

    auto g_0_xyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 357);

    auto g_0_xyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 358);

    auto g_0_xyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 359);

    auto g_0_xyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 360);

    auto g_0_xyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 361);

    auto g_0_xyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 362);

    auto g_0_xyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 363);

    auto g_0_xyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 364);

    auto g_0_xyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 365);

    auto g_0_xyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 366);

    auto g_0_xyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 367);

    auto g_0_xyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 368);

    auto g_0_xyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 369);

    auto g_0_xyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 370);

    auto g_0_xyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 371);

    auto g_0_xyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 372);

    auto g_0_xyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 373);

    auto g_0_xyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 374);

    auto g_0_xyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 375);

    auto g_0_xyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 376);

    auto g_0_xyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 377);

#pragma omp simd aligned(g_0_xyyyzz_0_xxxxx_0,     \
                             g_0_xyyyzz_0_xxxxy_0, \
                             g_0_xyyyzz_0_xxxxz_0, \
                             g_0_xyyyzz_0_xxxyy_0, \
                             g_0_xyyyzz_0_xxxyz_0, \
                             g_0_xyyyzz_0_xxxzz_0, \
                             g_0_xyyyzz_0_xxyyy_0, \
                             g_0_xyyyzz_0_xxyyz_0, \
                             g_0_xyyyzz_0_xxyzz_0, \
                             g_0_xyyyzz_0_xxzzz_0, \
                             g_0_xyyyzz_0_xyyyy_0, \
                             g_0_xyyyzz_0_xyyyz_0, \
                             g_0_xyyyzz_0_xyyzz_0, \
                             g_0_xyyyzz_0_xyzzz_0, \
                             g_0_xyyyzz_0_xzzzz_0, \
                             g_0_xyyyzz_0_yyyyy_0, \
                             g_0_xyyyzz_0_yyyyz_0, \
                             g_0_xyyyzz_0_yyyzz_0, \
                             g_0_xyyyzz_0_yyzzz_0, \
                             g_0_xyyyzz_0_yzzzz_0, \
                             g_0_xyyyzz_0_zzzzz_0, \
                             g_0_yyyzz_0_xxxx_1,   \
                             g_0_yyyzz_0_xxxxx_0,  \
                             g_0_yyyzz_0_xxxxx_1,  \
                             g_0_yyyzz_0_xxxxy_0,  \
                             g_0_yyyzz_0_xxxxy_1,  \
                             g_0_yyyzz_0_xxxxz_0,  \
                             g_0_yyyzz_0_xxxxz_1,  \
                             g_0_yyyzz_0_xxxy_1,   \
                             g_0_yyyzz_0_xxxyy_0,  \
                             g_0_yyyzz_0_xxxyy_1,  \
                             g_0_yyyzz_0_xxxyz_0,  \
                             g_0_yyyzz_0_xxxyz_1,  \
                             g_0_yyyzz_0_xxxz_1,   \
                             g_0_yyyzz_0_xxxzz_0,  \
                             g_0_yyyzz_0_xxxzz_1,  \
                             g_0_yyyzz_0_xxyy_1,   \
                             g_0_yyyzz_0_xxyyy_0,  \
                             g_0_yyyzz_0_xxyyy_1,  \
                             g_0_yyyzz_0_xxyyz_0,  \
                             g_0_yyyzz_0_xxyyz_1,  \
                             g_0_yyyzz_0_xxyz_1,   \
                             g_0_yyyzz_0_xxyzz_0,  \
                             g_0_yyyzz_0_xxyzz_1,  \
                             g_0_yyyzz_0_xxzz_1,   \
                             g_0_yyyzz_0_xxzzz_0,  \
                             g_0_yyyzz_0_xxzzz_1,  \
                             g_0_yyyzz_0_xyyy_1,   \
                             g_0_yyyzz_0_xyyyy_0,  \
                             g_0_yyyzz_0_xyyyy_1,  \
                             g_0_yyyzz_0_xyyyz_0,  \
                             g_0_yyyzz_0_xyyyz_1,  \
                             g_0_yyyzz_0_xyyz_1,   \
                             g_0_yyyzz_0_xyyzz_0,  \
                             g_0_yyyzz_0_xyyzz_1,  \
                             g_0_yyyzz_0_xyzz_1,   \
                             g_0_yyyzz_0_xyzzz_0,  \
                             g_0_yyyzz_0_xyzzz_1,  \
                             g_0_yyyzz_0_xzzz_1,   \
                             g_0_yyyzz_0_xzzzz_0,  \
                             g_0_yyyzz_0_xzzzz_1,  \
                             g_0_yyyzz_0_yyyy_1,   \
                             g_0_yyyzz_0_yyyyy_0,  \
                             g_0_yyyzz_0_yyyyy_1,  \
                             g_0_yyyzz_0_yyyyz_0,  \
                             g_0_yyyzz_0_yyyyz_1,  \
                             g_0_yyyzz_0_yyyz_1,   \
                             g_0_yyyzz_0_yyyzz_0,  \
                             g_0_yyyzz_0_yyyzz_1,  \
                             g_0_yyyzz_0_yyzz_1,   \
                             g_0_yyyzz_0_yyzzz_0,  \
                             g_0_yyyzz_0_yyzzz_1,  \
                             g_0_yyyzz_0_yzzz_1,   \
                             g_0_yyyzz_0_yzzzz_0,  \
                             g_0_yyyzz_0_yzzzz_1,  \
                             g_0_yyyzz_0_zzzz_1,   \
                             g_0_yyyzz_0_zzzzz_0,  \
                             g_0_yyyzz_0_zzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xxxxx_0[i] = 5.0 * g_0_yyyzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxx_0[i] * pb_x + g_0_yyyzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxy_0[i] = 4.0 * g_0_yyyzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxy_0[i] * pb_x + g_0_yyyzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxz_0[i] = 4.0 * g_0_yyyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxz_0[i] * pb_x + g_0_yyyzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyy_0[i] = 3.0 * g_0_yyyzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyy_0[i] * pb_x + g_0_yyyzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyz_0[i] = 3.0 * g_0_yyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyz_0[i] * pb_x + g_0_yyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxzz_0[i] = 3.0 * g_0_yyyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxzz_0[i] * pb_x + g_0_yyyzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyy_0[i] = 2.0 * g_0_yyyzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyy_0[i] * pb_x + g_0_yyyzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyz_0[i] = 2.0 * g_0_yyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyz_0[i] * pb_x + g_0_yyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyzz_0[i] = 2.0 * g_0_yyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzz_0[i] * pb_x + g_0_yyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxzzz_0[i] = 2.0 * g_0_yyyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxzzz_0[i] * pb_x + g_0_yyyzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyy_0[i] = g_0_yyyzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyy_0[i] * pb_x + g_0_yyyzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyz_0[i] = g_0_yyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyz_0[i] * pb_x + g_0_yyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyzz_0[i] = g_0_yyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzz_0[i] * pb_x + g_0_yyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyzzz_0[i] = g_0_yyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzz_0[i] * pb_x + g_0_yyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xzzzz_0[i] = g_0_yyyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xzzzz_0[i] * pb_x + g_0_yyyzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyy_0[i] = g_0_yyyzz_0_yyyyy_0[i] * pb_x + g_0_yyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyz_0[i] = g_0_yyyzz_0_yyyyz_0[i] * pb_x + g_0_yyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyzz_0[i] = g_0_yyyzz_0_yyyzz_0[i] * pb_x + g_0_yyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyzzz_0[i] = g_0_yyyzz_0_yyzzz_0[i] * pb_x + g_0_yyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yzzzz_0[i] = g_0_yyyzz_0_yzzzz_0[i] * pb_x + g_0_yyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zzzzz_0[i] = g_0_yyyzz_0_zzzzz_0[i] * pb_x + g_0_yyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 378-399 components of targeted buffer : SISH

    auto g_0_xyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 378);

    auto g_0_xyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 379);

    auto g_0_xyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 380);

    auto g_0_xyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 381);

    auto g_0_xyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 382);

    auto g_0_xyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 383);

    auto g_0_xyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 384);

    auto g_0_xyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 385);

    auto g_0_xyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 386);

    auto g_0_xyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 387);

    auto g_0_xyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 388);

    auto g_0_xyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 389);

    auto g_0_xyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 390);

    auto g_0_xyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 391);

    auto g_0_xyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 392);

    auto g_0_xyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 393);

    auto g_0_xyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 394);

    auto g_0_xyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 395);

    auto g_0_xyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 396);

    auto g_0_xyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 397);

    auto g_0_xyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 398);

#pragma omp simd aligned(g_0_xyyzzz_0_xxxxx_0,     \
                             g_0_xyyzzz_0_xxxxy_0, \
                             g_0_xyyzzz_0_xxxxz_0, \
                             g_0_xyyzzz_0_xxxyy_0, \
                             g_0_xyyzzz_0_xxxyz_0, \
                             g_0_xyyzzz_0_xxxzz_0, \
                             g_0_xyyzzz_0_xxyyy_0, \
                             g_0_xyyzzz_0_xxyyz_0, \
                             g_0_xyyzzz_0_xxyzz_0, \
                             g_0_xyyzzz_0_xxzzz_0, \
                             g_0_xyyzzz_0_xyyyy_0, \
                             g_0_xyyzzz_0_xyyyz_0, \
                             g_0_xyyzzz_0_xyyzz_0, \
                             g_0_xyyzzz_0_xyzzz_0, \
                             g_0_xyyzzz_0_xzzzz_0, \
                             g_0_xyyzzz_0_yyyyy_0, \
                             g_0_xyyzzz_0_yyyyz_0, \
                             g_0_xyyzzz_0_yyyzz_0, \
                             g_0_xyyzzz_0_yyzzz_0, \
                             g_0_xyyzzz_0_yzzzz_0, \
                             g_0_xyyzzz_0_zzzzz_0, \
                             g_0_yyzzz_0_xxxx_1,   \
                             g_0_yyzzz_0_xxxxx_0,  \
                             g_0_yyzzz_0_xxxxx_1,  \
                             g_0_yyzzz_0_xxxxy_0,  \
                             g_0_yyzzz_0_xxxxy_1,  \
                             g_0_yyzzz_0_xxxxz_0,  \
                             g_0_yyzzz_0_xxxxz_1,  \
                             g_0_yyzzz_0_xxxy_1,   \
                             g_0_yyzzz_0_xxxyy_0,  \
                             g_0_yyzzz_0_xxxyy_1,  \
                             g_0_yyzzz_0_xxxyz_0,  \
                             g_0_yyzzz_0_xxxyz_1,  \
                             g_0_yyzzz_0_xxxz_1,   \
                             g_0_yyzzz_0_xxxzz_0,  \
                             g_0_yyzzz_0_xxxzz_1,  \
                             g_0_yyzzz_0_xxyy_1,   \
                             g_0_yyzzz_0_xxyyy_0,  \
                             g_0_yyzzz_0_xxyyy_1,  \
                             g_0_yyzzz_0_xxyyz_0,  \
                             g_0_yyzzz_0_xxyyz_1,  \
                             g_0_yyzzz_0_xxyz_1,   \
                             g_0_yyzzz_0_xxyzz_0,  \
                             g_0_yyzzz_0_xxyzz_1,  \
                             g_0_yyzzz_0_xxzz_1,   \
                             g_0_yyzzz_0_xxzzz_0,  \
                             g_0_yyzzz_0_xxzzz_1,  \
                             g_0_yyzzz_0_xyyy_1,   \
                             g_0_yyzzz_0_xyyyy_0,  \
                             g_0_yyzzz_0_xyyyy_1,  \
                             g_0_yyzzz_0_xyyyz_0,  \
                             g_0_yyzzz_0_xyyyz_1,  \
                             g_0_yyzzz_0_xyyz_1,   \
                             g_0_yyzzz_0_xyyzz_0,  \
                             g_0_yyzzz_0_xyyzz_1,  \
                             g_0_yyzzz_0_xyzz_1,   \
                             g_0_yyzzz_0_xyzzz_0,  \
                             g_0_yyzzz_0_xyzzz_1,  \
                             g_0_yyzzz_0_xzzz_1,   \
                             g_0_yyzzz_0_xzzzz_0,  \
                             g_0_yyzzz_0_xzzzz_1,  \
                             g_0_yyzzz_0_yyyy_1,   \
                             g_0_yyzzz_0_yyyyy_0,  \
                             g_0_yyzzz_0_yyyyy_1,  \
                             g_0_yyzzz_0_yyyyz_0,  \
                             g_0_yyzzz_0_yyyyz_1,  \
                             g_0_yyzzz_0_yyyz_1,   \
                             g_0_yyzzz_0_yyyzz_0,  \
                             g_0_yyzzz_0_yyyzz_1,  \
                             g_0_yyzzz_0_yyzz_1,   \
                             g_0_yyzzz_0_yyzzz_0,  \
                             g_0_yyzzz_0_yyzzz_1,  \
                             g_0_yyzzz_0_yzzz_1,   \
                             g_0_yyzzz_0_yzzzz_0,  \
                             g_0_yyzzz_0_yzzzz_1,  \
                             g_0_yyzzz_0_zzzz_1,   \
                             g_0_yyzzz_0_zzzzz_0,  \
                             g_0_yyzzz_0_zzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xxxxx_0[i] = 5.0 * g_0_yyzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxx_0[i] * pb_x + g_0_yyzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxy_0[i] = 4.0 * g_0_yyzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxy_0[i] * pb_x + g_0_yyzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxz_0[i] = 4.0 * g_0_yyzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxz_0[i] * pb_x + g_0_yyzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyy_0[i] = 3.0 * g_0_yyzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyy_0[i] * pb_x + g_0_yyzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyz_0[i] = 3.0 * g_0_yyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyz_0[i] * pb_x + g_0_yyzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxzz_0[i] = 3.0 * g_0_yyzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxzz_0[i] * pb_x + g_0_yyzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyy_0[i] * pb_x + g_0_yyzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyz_0[i] = 2.0 * g_0_yyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyz_0[i] * pb_x + g_0_yyzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyzz_0[i] = 2.0 * g_0_yyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzz_0[i] * pb_x + g_0_yyzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxzzz_0[i] = 2.0 * g_0_yyzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxzzz_0[i] * pb_x + g_0_yyzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyy_0[i] = g_0_yyzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyy_0[i] * pb_x + g_0_yyzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyz_0[i] = g_0_yyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyz_0[i] * pb_x + g_0_yyzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyzz_0[i] = g_0_yyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzz_0[i] * pb_x + g_0_yyzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyzzz_0[i] = g_0_yyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzz_0[i] * pb_x + g_0_yyzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xzzzz_0[i] = g_0_yyzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xzzzz_0[i] * pb_x + g_0_yyzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyy_0[i] = g_0_yyzzz_0_yyyyy_0[i] * pb_x + g_0_yyzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyz_0[i] = g_0_yyzzz_0_yyyyz_0[i] * pb_x + g_0_yyzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyzz_0[i] = g_0_yyzzz_0_yyyzz_0[i] * pb_x + g_0_yyzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyzzz_0[i] = g_0_yyzzz_0_yyzzz_0[i] * pb_x + g_0_yyzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yzzzz_0[i] = g_0_yyzzz_0_yzzzz_0[i] * pb_x + g_0_yyzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zzzzz_0[i] = g_0_yyzzz_0_zzzzz_0[i] * pb_x + g_0_yyzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 399-420 components of targeted buffer : SISH

    auto g_0_xyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 399);

    auto g_0_xyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 400);

    auto g_0_xyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 401);

    auto g_0_xyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 402);

    auto g_0_xyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 403);

    auto g_0_xyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 404);

    auto g_0_xyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 405);

    auto g_0_xyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 406);

    auto g_0_xyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 407);

    auto g_0_xyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 408);

    auto g_0_xyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 409);

    auto g_0_xyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 410);

    auto g_0_xyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 411);

    auto g_0_xyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 412);

    auto g_0_xyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 413);

    auto g_0_xyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 414);

    auto g_0_xyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 415);

    auto g_0_xyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 416);

    auto g_0_xyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 417);

    auto g_0_xyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 418);

    auto g_0_xyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 419);

#pragma omp simd aligned(g_0_xyzzzz_0_xxxxx_0,     \
                             g_0_xyzzzz_0_xxxxy_0, \
                             g_0_xyzzzz_0_xxxxz_0, \
                             g_0_xyzzzz_0_xxxyy_0, \
                             g_0_xyzzzz_0_xxxyz_0, \
                             g_0_xyzzzz_0_xxxzz_0, \
                             g_0_xyzzzz_0_xxyyy_0, \
                             g_0_xyzzzz_0_xxyyz_0, \
                             g_0_xyzzzz_0_xxyzz_0, \
                             g_0_xyzzzz_0_xxzzz_0, \
                             g_0_xyzzzz_0_xyyyy_0, \
                             g_0_xyzzzz_0_xyyyz_0, \
                             g_0_xyzzzz_0_xyyzz_0, \
                             g_0_xyzzzz_0_xyzzz_0, \
                             g_0_xyzzzz_0_xzzzz_0, \
                             g_0_xyzzzz_0_yyyyy_0, \
                             g_0_xyzzzz_0_yyyyz_0, \
                             g_0_xyzzzz_0_yyyzz_0, \
                             g_0_xyzzzz_0_yyzzz_0, \
                             g_0_xyzzzz_0_yzzzz_0, \
                             g_0_xyzzzz_0_zzzzz_0, \
                             g_0_xzzzz_0_xxxxx_0,  \
                             g_0_xzzzz_0_xxxxx_1,  \
                             g_0_xzzzz_0_xxxxz_0,  \
                             g_0_xzzzz_0_xxxxz_1,  \
                             g_0_xzzzz_0_xxxzz_0,  \
                             g_0_xzzzz_0_xxxzz_1,  \
                             g_0_xzzzz_0_xxzzz_0,  \
                             g_0_xzzzz_0_xxzzz_1,  \
                             g_0_xzzzz_0_xzzzz_0,  \
                             g_0_xzzzz_0_xzzzz_1,  \
                             g_0_yzzzz_0_xxxxy_0,  \
                             g_0_yzzzz_0_xxxxy_1,  \
                             g_0_yzzzz_0_xxxy_1,   \
                             g_0_yzzzz_0_xxxyy_0,  \
                             g_0_yzzzz_0_xxxyy_1,  \
                             g_0_yzzzz_0_xxxyz_0,  \
                             g_0_yzzzz_0_xxxyz_1,  \
                             g_0_yzzzz_0_xxyy_1,   \
                             g_0_yzzzz_0_xxyyy_0,  \
                             g_0_yzzzz_0_xxyyy_1,  \
                             g_0_yzzzz_0_xxyyz_0,  \
                             g_0_yzzzz_0_xxyyz_1,  \
                             g_0_yzzzz_0_xxyz_1,   \
                             g_0_yzzzz_0_xxyzz_0,  \
                             g_0_yzzzz_0_xxyzz_1,  \
                             g_0_yzzzz_0_xyyy_1,   \
                             g_0_yzzzz_0_xyyyy_0,  \
                             g_0_yzzzz_0_xyyyy_1,  \
                             g_0_yzzzz_0_xyyyz_0,  \
                             g_0_yzzzz_0_xyyyz_1,  \
                             g_0_yzzzz_0_xyyz_1,   \
                             g_0_yzzzz_0_xyyzz_0,  \
                             g_0_yzzzz_0_xyyzz_1,  \
                             g_0_yzzzz_0_xyzz_1,   \
                             g_0_yzzzz_0_xyzzz_0,  \
                             g_0_yzzzz_0_xyzzz_1,  \
                             g_0_yzzzz_0_yyyy_1,   \
                             g_0_yzzzz_0_yyyyy_0,  \
                             g_0_yzzzz_0_yyyyy_1,  \
                             g_0_yzzzz_0_yyyyz_0,  \
                             g_0_yzzzz_0_yyyyz_1,  \
                             g_0_yzzzz_0_yyyz_1,   \
                             g_0_yzzzz_0_yyyzz_0,  \
                             g_0_yzzzz_0_yyyzz_1,  \
                             g_0_yzzzz_0_yyzz_1,   \
                             g_0_yzzzz_0_yyzzz_0,  \
                             g_0_yzzzz_0_yyzzz_1,  \
                             g_0_yzzzz_0_yzzz_1,   \
                             g_0_yzzzz_0_yzzzz_0,  \
                             g_0_yzzzz_0_yzzzz_1,  \
                             g_0_yzzzz_0_zzzzz_0,  \
                             g_0_yzzzz_0_zzzzz_1,  \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xxxxx_0[i] = g_0_xzzzz_0_xxxxx_0[i] * pb_y + g_0_xzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxy_0[i] = 4.0 * g_0_yzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxy_0[i] * pb_x + g_0_yzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxz_0[i] = g_0_xzzzz_0_xxxxz_0[i] * pb_y + g_0_xzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxyy_0[i] = 3.0 * g_0_yzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyy_0[i] * pb_x + g_0_yzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyz_0[i] = 3.0 * g_0_yzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyz_0[i] * pb_x + g_0_yzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxzz_0[i] = g_0_xzzzz_0_xxxzz_0[i] * pb_y + g_0_xzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxyyy_0[i] = 2.0 * g_0_yzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyy_0[i] * pb_x + g_0_yzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyz_0[i] * pb_x + g_0_yzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyzz_0[i] * pb_x + g_0_yzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxzzz_0[i] = g_0_xzzzz_0_xxzzz_0[i] * pb_y + g_0_xzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xyyyy_0[i] = g_0_yzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyy_0[i] * pb_x + g_0_yzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyz_0[i] = g_0_yzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyz_0[i] * pb_x + g_0_yzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyzz_0[i] = g_0_yzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyzz_0[i] * pb_x + g_0_yzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyzzz_0[i] = g_0_yzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyzzz_0[i] * pb_x + g_0_yzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xzzzz_0[i] = g_0_xzzzz_0_xzzzz_0[i] * pb_y + g_0_xzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yyyyy_0[i] = g_0_yzzzz_0_yyyyy_0[i] * pb_x + g_0_yzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyz_0[i] = g_0_yzzzz_0_yyyyz_0[i] * pb_x + g_0_yzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyzz_0[i] = g_0_yzzzz_0_yyyzz_0[i] * pb_x + g_0_yzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyzzz_0[i] = g_0_yzzzz_0_yyzzz_0[i] * pb_x + g_0_yzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yzzzz_0[i] = g_0_yzzzz_0_yzzzz_0[i] * pb_x + g_0_yzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zzzzz_0[i] = g_0_yzzzz_0_zzzzz_0[i] * pb_x + g_0_yzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 420-441 components of targeted buffer : SISH

    auto g_0_xzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 420);

    auto g_0_xzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 421);

    auto g_0_xzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 422);

    auto g_0_xzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 423);

    auto g_0_xzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 424);

    auto g_0_xzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 425);

    auto g_0_xzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 426);

    auto g_0_xzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 427);

    auto g_0_xzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 428);

    auto g_0_xzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 429);

    auto g_0_xzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 430);

    auto g_0_xzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 431);

    auto g_0_xzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 432);

    auto g_0_xzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 433);

    auto g_0_xzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 434);

    auto g_0_xzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 435);

    auto g_0_xzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 436);

    auto g_0_xzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 437);

    auto g_0_xzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 438);

    auto g_0_xzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 439);

    auto g_0_xzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 440);

#pragma omp simd aligned(g_0_xzzzzz_0_xxxxx_0,     \
                             g_0_xzzzzz_0_xxxxy_0, \
                             g_0_xzzzzz_0_xxxxz_0, \
                             g_0_xzzzzz_0_xxxyy_0, \
                             g_0_xzzzzz_0_xxxyz_0, \
                             g_0_xzzzzz_0_xxxzz_0, \
                             g_0_xzzzzz_0_xxyyy_0, \
                             g_0_xzzzzz_0_xxyyz_0, \
                             g_0_xzzzzz_0_xxyzz_0, \
                             g_0_xzzzzz_0_xxzzz_0, \
                             g_0_xzzzzz_0_xyyyy_0, \
                             g_0_xzzzzz_0_xyyyz_0, \
                             g_0_xzzzzz_0_xyyzz_0, \
                             g_0_xzzzzz_0_xyzzz_0, \
                             g_0_xzzzzz_0_xzzzz_0, \
                             g_0_xzzzzz_0_yyyyy_0, \
                             g_0_xzzzzz_0_yyyyz_0, \
                             g_0_xzzzzz_0_yyyzz_0, \
                             g_0_xzzzzz_0_yyzzz_0, \
                             g_0_xzzzzz_0_yzzzz_0, \
                             g_0_xzzzzz_0_zzzzz_0, \
                             g_0_zzzzz_0_xxxx_1,   \
                             g_0_zzzzz_0_xxxxx_0,  \
                             g_0_zzzzz_0_xxxxx_1,  \
                             g_0_zzzzz_0_xxxxy_0,  \
                             g_0_zzzzz_0_xxxxy_1,  \
                             g_0_zzzzz_0_xxxxz_0,  \
                             g_0_zzzzz_0_xxxxz_1,  \
                             g_0_zzzzz_0_xxxy_1,   \
                             g_0_zzzzz_0_xxxyy_0,  \
                             g_0_zzzzz_0_xxxyy_1,  \
                             g_0_zzzzz_0_xxxyz_0,  \
                             g_0_zzzzz_0_xxxyz_1,  \
                             g_0_zzzzz_0_xxxz_1,   \
                             g_0_zzzzz_0_xxxzz_0,  \
                             g_0_zzzzz_0_xxxzz_1,  \
                             g_0_zzzzz_0_xxyy_1,   \
                             g_0_zzzzz_0_xxyyy_0,  \
                             g_0_zzzzz_0_xxyyy_1,  \
                             g_0_zzzzz_0_xxyyz_0,  \
                             g_0_zzzzz_0_xxyyz_1,  \
                             g_0_zzzzz_0_xxyz_1,   \
                             g_0_zzzzz_0_xxyzz_0,  \
                             g_0_zzzzz_0_xxyzz_1,  \
                             g_0_zzzzz_0_xxzz_1,   \
                             g_0_zzzzz_0_xxzzz_0,  \
                             g_0_zzzzz_0_xxzzz_1,  \
                             g_0_zzzzz_0_xyyy_1,   \
                             g_0_zzzzz_0_xyyyy_0,  \
                             g_0_zzzzz_0_xyyyy_1,  \
                             g_0_zzzzz_0_xyyyz_0,  \
                             g_0_zzzzz_0_xyyyz_1,  \
                             g_0_zzzzz_0_xyyz_1,   \
                             g_0_zzzzz_0_xyyzz_0,  \
                             g_0_zzzzz_0_xyyzz_1,  \
                             g_0_zzzzz_0_xyzz_1,   \
                             g_0_zzzzz_0_xyzzz_0,  \
                             g_0_zzzzz_0_xyzzz_1,  \
                             g_0_zzzzz_0_xzzz_1,   \
                             g_0_zzzzz_0_xzzzz_0,  \
                             g_0_zzzzz_0_xzzzz_1,  \
                             g_0_zzzzz_0_yyyy_1,   \
                             g_0_zzzzz_0_yyyyy_0,  \
                             g_0_zzzzz_0_yyyyy_1,  \
                             g_0_zzzzz_0_yyyyz_0,  \
                             g_0_zzzzz_0_yyyyz_1,  \
                             g_0_zzzzz_0_yyyz_1,   \
                             g_0_zzzzz_0_yyyzz_0,  \
                             g_0_zzzzz_0_yyyzz_1,  \
                             g_0_zzzzz_0_yyzz_1,   \
                             g_0_zzzzz_0_yyzzz_0,  \
                             g_0_zzzzz_0_yyzzz_1,  \
                             g_0_zzzzz_0_yzzz_1,   \
                             g_0_zzzzz_0_yzzzz_0,  \
                             g_0_zzzzz_0_yzzzz_1,  \
                             g_0_zzzzz_0_zzzz_1,   \
                             g_0_zzzzz_0_zzzzz_0,  \
                             g_0_zzzzz_0_zzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xxxxx_0[i] = 5.0 * g_0_zzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxx_0[i] * pb_x + g_0_zzzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxy_0[i] = 4.0 * g_0_zzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxy_0[i] * pb_x + g_0_zzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxz_0[i] = 4.0 * g_0_zzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxz_0[i] * pb_x + g_0_zzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyy_0[i] = 3.0 * g_0_zzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyy_0[i] * pb_x + g_0_zzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyz_0[i] = 3.0 * g_0_zzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyz_0[i] * pb_x + g_0_zzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxzz_0[i] = 3.0 * g_0_zzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzz_0[i] * pb_x + g_0_zzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyy_0[i] = 2.0 * g_0_zzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyy_0[i] * pb_x + g_0_zzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyz_0[i] * pb_x + g_0_zzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyzz_0[i] = 2.0 * g_0_zzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzz_0[i] * pb_x + g_0_zzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxzzz_0[i] = 2.0 * g_0_zzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzz_0[i] * pb_x + g_0_zzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyy_0[i] = g_0_zzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyy_0[i] * pb_x + g_0_zzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyz_0[i] = g_0_zzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyz_0[i] * pb_x + g_0_zzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyzz_0[i] = g_0_zzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzz_0[i] * pb_x + g_0_zzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyzzz_0[i] = g_0_zzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzz_0[i] * pb_x + g_0_zzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xzzzz_0[i] = g_0_zzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzz_0[i] * pb_x + g_0_zzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyy_0[i] = g_0_zzzzz_0_yyyyy_0[i] * pb_x + g_0_zzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyz_0[i] = g_0_zzzzz_0_yyyyz_0[i] * pb_x + g_0_zzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyzz_0[i] = g_0_zzzzz_0_yyyzz_0[i] * pb_x + g_0_zzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyzzz_0[i] = g_0_zzzzz_0_yyzzz_0[i] * pb_x + g_0_zzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yzzzz_0[i] = g_0_zzzzz_0_yzzzz_0[i] * pb_x + g_0_zzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zzzzz_0[i] = g_0_zzzzz_0_zzzzz_0[i] * pb_x + g_0_zzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 441-462 components of targeted buffer : SISH

    auto g_0_yyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 441);

    auto g_0_yyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 442);

    auto g_0_yyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 443);

    auto g_0_yyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 444);

    auto g_0_yyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 445);

    auto g_0_yyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 446);

    auto g_0_yyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 447);

    auto g_0_yyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 448);

    auto g_0_yyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 449);

    auto g_0_yyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 450);

    auto g_0_yyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 451);

    auto g_0_yyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 452);

    auto g_0_yyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 453);

    auto g_0_yyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 454);

    auto g_0_yyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 455);

    auto g_0_yyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 456);

    auto g_0_yyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 457);

    auto g_0_yyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 458);

    auto g_0_yyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 459);

    auto g_0_yyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 460);

    auto g_0_yyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 461);

#pragma omp simd aligned(g_0_yyyy_0_xxxxx_0,       \
                             g_0_yyyy_0_xxxxx_1,   \
                             g_0_yyyy_0_xxxxy_0,   \
                             g_0_yyyy_0_xxxxy_1,   \
                             g_0_yyyy_0_xxxxz_0,   \
                             g_0_yyyy_0_xxxxz_1,   \
                             g_0_yyyy_0_xxxyy_0,   \
                             g_0_yyyy_0_xxxyy_1,   \
                             g_0_yyyy_0_xxxyz_0,   \
                             g_0_yyyy_0_xxxyz_1,   \
                             g_0_yyyy_0_xxxzz_0,   \
                             g_0_yyyy_0_xxxzz_1,   \
                             g_0_yyyy_0_xxyyy_0,   \
                             g_0_yyyy_0_xxyyy_1,   \
                             g_0_yyyy_0_xxyyz_0,   \
                             g_0_yyyy_0_xxyyz_1,   \
                             g_0_yyyy_0_xxyzz_0,   \
                             g_0_yyyy_0_xxyzz_1,   \
                             g_0_yyyy_0_xxzzz_0,   \
                             g_0_yyyy_0_xxzzz_1,   \
                             g_0_yyyy_0_xyyyy_0,   \
                             g_0_yyyy_0_xyyyy_1,   \
                             g_0_yyyy_0_xyyyz_0,   \
                             g_0_yyyy_0_xyyyz_1,   \
                             g_0_yyyy_0_xyyzz_0,   \
                             g_0_yyyy_0_xyyzz_1,   \
                             g_0_yyyy_0_xyzzz_0,   \
                             g_0_yyyy_0_xyzzz_1,   \
                             g_0_yyyy_0_xzzzz_0,   \
                             g_0_yyyy_0_xzzzz_1,   \
                             g_0_yyyy_0_yyyyy_0,   \
                             g_0_yyyy_0_yyyyy_1,   \
                             g_0_yyyy_0_yyyyz_0,   \
                             g_0_yyyy_0_yyyyz_1,   \
                             g_0_yyyy_0_yyyzz_0,   \
                             g_0_yyyy_0_yyyzz_1,   \
                             g_0_yyyy_0_yyzzz_0,   \
                             g_0_yyyy_0_yyzzz_1,   \
                             g_0_yyyy_0_yzzzz_0,   \
                             g_0_yyyy_0_yzzzz_1,   \
                             g_0_yyyy_0_zzzzz_0,   \
                             g_0_yyyy_0_zzzzz_1,   \
                             g_0_yyyyy_0_xxxx_1,   \
                             g_0_yyyyy_0_xxxxx_0,  \
                             g_0_yyyyy_0_xxxxx_1,  \
                             g_0_yyyyy_0_xxxxy_0,  \
                             g_0_yyyyy_0_xxxxy_1,  \
                             g_0_yyyyy_0_xxxxz_0,  \
                             g_0_yyyyy_0_xxxxz_1,  \
                             g_0_yyyyy_0_xxxy_1,   \
                             g_0_yyyyy_0_xxxyy_0,  \
                             g_0_yyyyy_0_xxxyy_1,  \
                             g_0_yyyyy_0_xxxyz_0,  \
                             g_0_yyyyy_0_xxxyz_1,  \
                             g_0_yyyyy_0_xxxz_1,   \
                             g_0_yyyyy_0_xxxzz_0,  \
                             g_0_yyyyy_0_xxxzz_1,  \
                             g_0_yyyyy_0_xxyy_1,   \
                             g_0_yyyyy_0_xxyyy_0,  \
                             g_0_yyyyy_0_xxyyy_1,  \
                             g_0_yyyyy_0_xxyyz_0,  \
                             g_0_yyyyy_0_xxyyz_1,  \
                             g_0_yyyyy_0_xxyz_1,   \
                             g_0_yyyyy_0_xxyzz_0,  \
                             g_0_yyyyy_0_xxyzz_1,  \
                             g_0_yyyyy_0_xxzz_1,   \
                             g_0_yyyyy_0_xxzzz_0,  \
                             g_0_yyyyy_0_xxzzz_1,  \
                             g_0_yyyyy_0_xyyy_1,   \
                             g_0_yyyyy_0_xyyyy_0,  \
                             g_0_yyyyy_0_xyyyy_1,  \
                             g_0_yyyyy_0_xyyyz_0,  \
                             g_0_yyyyy_0_xyyyz_1,  \
                             g_0_yyyyy_0_xyyz_1,   \
                             g_0_yyyyy_0_xyyzz_0,  \
                             g_0_yyyyy_0_xyyzz_1,  \
                             g_0_yyyyy_0_xyzz_1,   \
                             g_0_yyyyy_0_xyzzz_0,  \
                             g_0_yyyyy_0_xyzzz_1,  \
                             g_0_yyyyy_0_xzzz_1,   \
                             g_0_yyyyy_0_xzzzz_0,  \
                             g_0_yyyyy_0_xzzzz_1,  \
                             g_0_yyyyy_0_yyyy_1,   \
                             g_0_yyyyy_0_yyyyy_0,  \
                             g_0_yyyyy_0_yyyyy_1,  \
                             g_0_yyyyy_0_yyyyz_0,  \
                             g_0_yyyyy_0_yyyyz_1,  \
                             g_0_yyyyy_0_yyyz_1,   \
                             g_0_yyyyy_0_yyyzz_0,  \
                             g_0_yyyyy_0_yyyzz_1,  \
                             g_0_yyyyy_0_yyzz_1,   \
                             g_0_yyyyy_0_yyzzz_0,  \
                             g_0_yyyyy_0_yyzzz_1,  \
                             g_0_yyyyy_0_yzzz_1,   \
                             g_0_yyyyy_0_yzzzz_0,  \
                             g_0_yyyyy_0_yzzzz_1,  \
                             g_0_yyyyy_0_zzzz_1,   \
                             g_0_yyyyy_0_zzzzz_0,  \
                             g_0_yyyyy_0_zzzzz_1,  \
                             g_0_yyyyyy_0_xxxxx_0, \
                             g_0_yyyyyy_0_xxxxy_0, \
                             g_0_yyyyyy_0_xxxxz_0, \
                             g_0_yyyyyy_0_xxxyy_0, \
                             g_0_yyyyyy_0_xxxyz_0, \
                             g_0_yyyyyy_0_xxxzz_0, \
                             g_0_yyyyyy_0_xxyyy_0, \
                             g_0_yyyyyy_0_xxyyz_0, \
                             g_0_yyyyyy_0_xxyzz_0, \
                             g_0_yyyyyy_0_xxzzz_0, \
                             g_0_yyyyyy_0_xyyyy_0, \
                             g_0_yyyyyy_0_xyyyz_0, \
                             g_0_yyyyyy_0_xyyzz_0, \
                             g_0_yyyyyy_0_xyzzz_0, \
                             g_0_yyyyyy_0_xzzzz_0, \
                             g_0_yyyyyy_0_yyyyy_0, \
                             g_0_yyyyyy_0_yyyyz_0, \
                             g_0_yyyyyy_0_yyyzz_0, \
                             g_0_yyyyyy_0_yyzzz_0, \
                             g_0_yyyyyy_0_yzzzz_0, \
                             g_0_yyyyyy_0_zzzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xxxxx_0[i] = 5.0 * g_0_yyyy_0_xxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxx_0[i] * pb_y +
                                  g_0_yyyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxy_0[i] = 5.0 * g_0_yyyy_0_xxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxx_1[i] * fi_abcd_0 +
                                  g_0_yyyyy_0_xxxxy_0[i] * pb_y + g_0_yyyyy_0_xxxxy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxz_0[i] = 5.0 * g_0_yyyy_0_xxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxz_0[i] * pb_y +
                                  g_0_yyyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyy_0[i] = 5.0 * g_0_yyyy_0_xxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyy_0[i] * pb_y + g_0_yyyyy_0_xxxyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyz_0[i] = 5.0 * g_0_yyyy_0_xxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxz_1[i] * fi_abcd_0 +
                                  g_0_yyyyy_0_xxxyz_0[i] * pb_y + g_0_yyyyy_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxzz_0[i] = 5.0 * g_0_yyyy_0_xxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxzz_0[i] * pb_y +
                                  g_0_yyyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyy_0[i] = 5.0 * g_0_yyyy_0_xxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyy_0[i] * pb_y + g_0_yyyyy_0_xxyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyz_0[i] = 5.0 * g_0_yyyy_0_xxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyz_0[i] * pb_y + g_0_yyyyy_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyzz_0[i] = 5.0 * g_0_yyyy_0_xxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyy_0_xxyzz_0[i] * pb_y + g_0_yyyyy_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxzzz_0[i] = 5.0 * g_0_yyyy_0_xxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxzzz_0[i] * pb_y +
                                  g_0_yyyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyy_0[i] = 5.0 * g_0_yyyy_0_xyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyy_0[i] * pb_y + g_0_yyyyy_0_xyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyz_0[i] = 5.0 * g_0_yyyy_0_xyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyz_0[i] * pb_y + g_0_yyyyy_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyzz_0[i] = 5.0 * g_0_yyyy_0_xyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzz_0[i] * pb_y + g_0_yyyyy_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyzzz_0[i] = 5.0 * g_0_yyyy_0_xyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyy_0_xyzzz_0[i] * pb_y + g_0_yyyyy_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xzzzz_0[i] = 5.0 * g_0_yyyy_0_xzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzzzz_0[i] * pb_y +
                                  g_0_yyyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyy_0[i] = 5.0 * g_0_yyyy_0_yyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyy_1[i] * fti_ab_0 +
                                  5.0 * g_0_yyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyy_0[i] * pb_y + g_0_yyyyy_0_yyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyz_0[i] = 5.0 * g_0_yyyy_0_yyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyz_0[i] * pb_y + g_0_yyyyy_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyzz_0[i] = 5.0 * g_0_yyyy_0_yyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzz_0[i] * pb_y + g_0_yyyyy_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyzzz_0[i] = 5.0 * g_0_yyyy_0_yyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzz_0[i] * pb_y + g_0_yyyyy_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yzzzz_0[i] = 5.0 * g_0_yyyy_0_yzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyy_0_yzzzz_0[i] * pb_y + g_0_yyyyy_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zzzzz_0[i] = 5.0 * g_0_yyyy_0_zzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzzzz_0[i] * pb_y +
                                  g_0_yyyyy_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 462-483 components of targeted buffer : SISH

    auto g_0_yyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 462);

    auto g_0_yyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 463);

    auto g_0_yyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 464);

    auto g_0_yyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 465);

    auto g_0_yyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 466);

    auto g_0_yyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 467);

    auto g_0_yyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 468);

    auto g_0_yyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 469);

    auto g_0_yyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 470);

    auto g_0_yyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 471);

    auto g_0_yyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 472);

    auto g_0_yyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 473);

    auto g_0_yyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 474);

    auto g_0_yyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 475);

    auto g_0_yyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 476);

    auto g_0_yyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 477);

    auto g_0_yyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 478);

    auto g_0_yyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 479);

    auto g_0_yyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 480);

    auto g_0_yyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 481);

    auto g_0_yyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 482);

#pragma omp simd aligned(g_0_yyyyy_0_xxxx_1,       \
                             g_0_yyyyy_0_xxxxx_0,  \
                             g_0_yyyyy_0_xxxxx_1,  \
                             g_0_yyyyy_0_xxxxy_0,  \
                             g_0_yyyyy_0_xxxxy_1,  \
                             g_0_yyyyy_0_xxxxz_0,  \
                             g_0_yyyyy_0_xxxxz_1,  \
                             g_0_yyyyy_0_xxxy_1,   \
                             g_0_yyyyy_0_xxxyy_0,  \
                             g_0_yyyyy_0_xxxyy_1,  \
                             g_0_yyyyy_0_xxxyz_0,  \
                             g_0_yyyyy_0_xxxyz_1,  \
                             g_0_yyyyy_0_xxxz_1,   \
                             g_0_yyyyy_0_xxxzz_0,  \
                             g_0_yyyyy_0_xxxzz_1,  \
                             g_0_yyyyy_0_xxyy_1,   \
                             g_0_yyyyy_0_xxyyy_0,  \
                             g_0_yyyyy_0_xxyyy_1,  \
                             g_0_yyyyy_0_xxyyz_0,  \
                             g_0_yyyyy_0_xxyyz_1,  \
                             g_0_yyyyy_0_xxyz_1,   \
                             g_0_yyyyy_0_xxyzz_0,  \
                             g_0_yyyyy_0_xxyzz_1,  \
                             g_0_yyyyy_0_xxzz_1,   \
                             g_0_yyyyy_0_xxzzz_0,  \
                             g_0_yyyyy_0_xxzzz_1,  \
                             g_0_yyyyy_0_xyyy_1,   \
                             g_0_yyyyy_0_xyyyy_0,  \
                             g_0_yyyyy_0_xyyyy_1,  \
                             g_0_yyyyy_0_xyyyz_0,  \
                             g_0_yyyyy_0_xyyyz_1,  \
                             g_0_yyyyy_0_xyyz_1,   \
                             g_0_yyyyy_0_xyyzz_0,  \
                             g_0_yyyyy_0_xyyzz_1,  \
                             g_0_yyyyy_0_xyzz_1,   \
                             g_0_yyyyy_0_xyzzz_0,  \
                             g_0_yyyyy_0_xyzzz_1,  \
                             g_0_yyyyy_0_xzzz_1,   \
                             g_0_yyyyy_0_xzzzz_0,  \
                             g_0_yyyyy_0_xzzzz_1,  \
                             g_0_yyyyy_0_yyyy_1,   \
                             g_0_yyyyy_0_yyyyy_0,  \
                             g_0_yyyyy_0_yyyyy_1,  \
                             g_0_yyyyy_0_yyyyz_0,  \
                             g_0_yyyyy_0_yyyyz_1,  \
                             g_0_yyyyy_0_yyyz_1,   \
                             g_0_yyyyy_0_yyyzz_0,  \
                             g_0_yyyyy_0_yyyzz_1,  \
                             g_0_yyyyy_0_yyzz_1,   \
                             g_0_yyyyy_0_yyzzz_0,  \
                             g_0_yyyyy_0_yyzzz_1,  \
                             g_0_yyyyy_0_yzzz_1,   \
                             g_0_yyyyy_0_yzzzz_0,  \
                             g_0_yyyyy_0_yzzzz_1,  \
                             g_0_yyyyy_0_zzzz_1,   \
                             g_0_yyyyy_0_zzzzz_0,  \
                             g_0_yyyyy_0_zzzzz_1,  \
                             g_0_yyyyyz_0_xxxxx_0, \
                             g_0_yyyyyz_0_xxxxy_0, \
                             g_0_yyyyyz_0_xxxxz_0, \
                             g_0_yyyyyz_0_xxxyy_0, \
                             g_0_yyyyyz_0_xxxyz_0, \
                             g_0_yyyyyz_0_xxxzz_0, \
                             g_0_yyyyyz_0_xxyyy_0, \
                             g_0_yyyyyz_0_xxyyz_0, \
                             g_0_yyyyyz_0_xxyzz_0, \
                             g_0_yyyyyz_0_xxzzz_0, \
                             g_0_yyyyyz_0_xyyyy_0, \
                             g_0_yyyyyz_0_xyyyz_0, \
                             g_0_yyyyyz_0_xyyzz_0, \
                             g_0_yyyyyz_0_xyzzz_0, \
                             g_0_yyyyyz_0_xzzzz_0, \
                             g_0_yyyyyz_0_yyyyy_0, \
                             g_0_yyyyyz_0_yyyyz_0, \
                             g_0_yyyyyz_0_yyyzz_0, \
                             g_0_yyyyyz_0_yyzzz_0, \
                             g_0_yyyyyz_0_yzzzz_0, \
                             g_0_yyyyyz_0_zzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xxxxx_0[i] = g_0_yyyyy_0_xxxxx_0[i] * pb_z + g_0_yyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxy_0[i] = g_0_yyyyy_0_xxxxy_0[i] * pb_z + g_0_yyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxz_0[i] = g_0_yyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxz_0[i] * pb_z + g_0_yyyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyy_0[i] = g_0_yyyyy_0_xxxyy_0[i] * pb_z + g_0_yyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyz_0[i] = g_0_yyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyz_0[i] * pb_z + g_0_yyyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxzz_0[i] = 2.0 * g_0_yyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzz_0[i] * pb_z + g_0_yyyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyy_0[i] = g_0_yyyyy_0_xxyyy_0[i] * pb_z + g_0_yyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyz_0[i] = g_0_yyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyz_0[i] * pb_z + g_0_yyyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzz_0[i] * pb_z + g_0_yyyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxzzz_0[i] = 3.0 * g_0_yyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzz_0[i] * pb_z + g_0_yyyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyy_0[i] = g_0_yyyyy_0_xyyyy_0[i] * pb_z + g_0_yyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyz_0[i] = g_0_yyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyz_0[i] * pb_z + g_0_yyyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyzz_0[i] = 2.0 * g_0_yyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzz_0[i] * pb_z + g_0_yyyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyzzz_0[i] = 3.0 * g_0_yyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzz_0[i] * pb_z + g_0_yyyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xzzzz_0[i] = 4.0 * g_0_yyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzz_0[i] * pb_z + g_0_yyyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyy_0[i] = g_0_yyyyy_0_yyyyy_0[i] * pb_z + g_0_yyyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyz_0[i] = g_0_yyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyz_0[i] * pb_z + g_0_yyyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyzz_0[i] = 2.0 * g_0_yyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzz_0[i] * pb_z + g_0_yyyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyzzz_0[i] = 3.0 * g_0_yyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzz_0[i] * pb_z + g_0_yyyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yzzzz_0[i] = 4.0 * g_0_yyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzz_0[i] * pb_z + g_0_yyyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zzzzz_0[i] = 5.0 * g_0_yyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_zzzzz_0[i] * pb_z + g_0_yyyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 483-504 components of targeted buffer : SISH

    auto g_0_yyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 483);

    auto g_0_yyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 484);

    auto g_0_yyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 485);

    auto g_0_yyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 486);

    auto g_0_yyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 487);

    auto g_0_yyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 488);

    auto g_0_yyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 489);

    auto g_0_yyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 490);

    auto g_0_yyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 491);

    auto g_0_yyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 492);

    auto g_0_yyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 493);

    auto g_0_yyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 494);

    auto g_0_yyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 495);

    auto g_0_yyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 496);

    auto g_0_yyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 497);

    auto g_0_yyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 498);

    auto g_0_yyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 499);

    auto g_0_yyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 500);

    auto g_0_yyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 501);

    auto g_0_yyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 502);

    auto g_0_yyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 503);

#pragma omp simd aligned(g_0_yyyy_0_xxxxy_0,       \
                             g_0_yyyy_0_xxxxy_1,   \
                             g_0_yyyy_0_xxxyy_0,   \
                             g_0_yyyy_0_xxxyy_1,   \
                             g_0_yyyy_0_xxyyy_0,   \
                             g_0_yyyy_0_xxyyy_1,   \
                             g_0_yyyy_0_xyyyy_0,   \
                             g_0_yyyy_0_xyyyy_1,   \
                             g_0_yyyy_0_yyyyy_0,   \
                             g_0_yyyy_0_yyyyy_1,   \
                             g_0_yyyyz_0_xxxxy_0,  \
                             g_0_yyyyz_0_xxxxy_1,  \
                             g_0_yyyyz_0_xxxyy_0,  \
                             g_0_yyyyz_0_xxxyy_1,  \
                             g_0_yyyyz_0_xxyyy_0,  \
                             g_0_yyyyz_0_xxyyy_1,  \
                             g_0_yyyyz_0_xyyyy_0,  \
                             g_0_yyyyz_0_xyyyy_1,  \
                             g_0_yyyyz_0_yyyyy_0,  \
                             g_0_yyyyz_0_yyyyy_1,  \
                             g_0_yyyyzz_0_xxxxx_0, \
                             g_0_yyyyzz_0_xxxxy_0, \
                             g_0_yyyyzz_0_xxxxz_0, \
                             g_0_yyyyzz_0_xxxyy_0, \
                             g_0_yyyyzz_0_xxxyz_0, \
                             g_0_yyyyzz_0_xxxzz_0, \
                             g_0_yyyyzz_0_xxyyy_0, \
                             g_0_yyyyzz_0_xxyyz_0, \
                             g_0_yyyyzz_0_xxyzz_0, \
                             g_0_yyyyzz_0_xxzzz_0, \
                             g_0_yyyyzz_0_xyyyy_0, \
                             g_0_yyyyzz_0_xyyyz_0, \
                             g_0_yyyyzz_0_xyyzz_0, \
                             g_0_yyyyzz_0_xyzzz_0, \
                             g_0_yyyyzz_0_xzzzz_0, \
                             g_0_yyyyzz_0_yyyyy_0, \
                             g_0_yyyyzz_0_yyyyz_0, \
                             g_0_yyyyzz_0_yyyzz_0, \
                             g_0_yyyyzz_0_yyzzz_0, \
                             g_0_yyyyzz_0_yzzzz_0, \
                             g_0_yyyyzz_0_zzzzz_0, \
                             g_0_yyyzz_0_xxxxx_0,  \
                             g_0_yyyzz_0_xxxxx_1,  \
                             g_0_yyyzz_0_xxxxz_0,  \
                             g_0_yyyzz_0_xxxxz_1,  \
                             g_0_yyyzz_0_xxxyz_0,  \
                             g_0_yyyzz_0_xxxyz_1,  \
                             g_0_yyyzz_0_xxxz_1,   \
                             g_0_yyyzz_0_xxxzz_0,  \
                             g_0_yyyzz_0_xxxzz_1,  \
                             g_0_yyyzz_0_xxyyz_0,  \
                             g_0_yyyzz_0_xxyyz_1,  \
                             g_0_yyyzz_0_xxyz_1,   \
                             g_0_yyyzz_0_xxyzz_0,  \
                             g_0_yyyzz_0_xxyzz_1,  \
                             g_0_yyyzz_0_xxzz_1,   \
                             g_0_yyyzz_0_xxzzz_0,  \
                             g_0_yyyzz_0_xxzzz_1,  \
                             g_0_yyyzz_0_xyyyz_0,  \
                             g_0_yyyzz_0_xyyyz_1,  \
                             g_0_yyyzz_0_xyyz_1,   \
                             g_0_yyyzz_0_xyyzz_0,  \
                             g_0_yyyzz_0_xyyzz_1,  \
                             g_0_yyyzz_0_xyzz_1,   \
                             g_0_yyyzz_0_xyzzz_0,  \
                             g_0_yyyzz_0_xyzzz_1,  \
                             g_0_yyyzz_0_xzzz_1,   \
                             g_0_yyyzz_0_xzzzz_0,  \
                             g_0_yyyzz_0_xzzzz_1,  \
                             g_0_yyyzz_0_yyyyz_0,  \
                             g_0_yyyzz_0_yyyyz_1,  \
                             g_0_yyyzz_0_yyyz_1,   \
                             g_0_yyyzz_0_yyyzz_0,  \
                             g_0_yyyzz_0_yyyzz_1,  \
                             g_0_yyyzz_0_yyzz_1,   \
                             g_0_yyyzz_0_yyzzz_0,  \
                             g_0_yyyzz_0_yyzzz_1,  \
                             g_0_yyyzz_0_yzzz_1,   \
                             g_0_yyyzz_0_yzzzz_0,  \
                             g_0_yyyzz_0_yzzzz_1,  \
                             g_0_yyyzz_0_zzzz_1,   \
                             g_0_yyyzz_0_zzzzz_0,  \
                             g_0_yyyzz_0_zzzzz_1,  \
                             g_0_yyzz_0_xxxxx_0,   \
                             g_0_yyzz_0_xxxxx_1,   \
                             g_0_yyzz_0_xxxxz_0,   \
                             g_0_yyzz_0_xxxxz_1,   \
                             g_0_yyzz_0_xxxyz_0,   \
                             g_0_yyzz_0_xxxyz_1,   \
                             g_0_yyzz_0_xxxzz_0,   \
                             g_0_yyzz_0_xxxzz_1,   \
                             g_0_yyzz_0_xxyyz_0,   \
                             g_0_yyzz_0_xxyyz_1,   \
                             g_0_yyzz_0_xxyzz_0,   \
                             g_0_yyzz_0_xxyzz_1,   \
                             g_0_yyzz_0_xxzzz_0,   \
                             g_0_yyzz_0_xxzzz_1,   \
                             g_0_yyzz_0_xyyyz_0,   \
                             g_0_yyzz_0_xyyyz_1,   \
                             g_0_yyzz_0_xyyzz_0,   \
                             g_0_yyzz_0_xyyzz_1,   \
                             g_0_yyzz_0_xyzzz_0,   \
                             g_0_yyzz_0_xyzzz_1,   \
                             g_0_yyzz_0_xzzzz_0,   \
                             g_0_yyzz_0_xzzzz_1,   \
                             g_0_yyzz_0_yyyyz_0,   \
                             g_0_yyzz_0_yyyyz_1,   \
                             g_0_yyzz_0_yyyzz_0,   \
                             g_0_yyzz_0_yyyzz_1,   \
                             g_0_yyzz_0_yyzzz_0,   \
                             g_0_yyzz_0_yyzzz_1,   \
                             g_0_yyzz_0_yzzzz_0,   \
                             g_0_yyzz_0_yzzzz_1,   \
                             g_0_yyzz_0_zzzzz_0,   \
                             g_0_yyzz_0_zzzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xxxxx_0[i] = 3.0 * g_0_yyzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxx_0[i] * pb_y +
                                  g_0_yyyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxy_0[i] =
            g_0_yyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxy_0[i] * pb_z + g_0_yyyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxz_0[i] = 3.0 * g_0_yyzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxz_0[i] * pb_y +
                                  g_0_yyyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyy_0[i] =
            g_0_yyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxyy_0[i] * pb_z + g_0_yyyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxyz_0[i] = 3.0 * g_0_yyzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxz_1[i] * fi_abcd_0 +
                                  g_0_yyyzz_0_xxxyz_0[i] * pb_y + g_0_yyyzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxzz_0[i] = 3.0 * g_0_yyzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxzz_0[i] * pb_y +
                                  g_0_yyyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyy_0[i] =
            g_0_yyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxyyy_0[i] * pb_z + g_0_yyyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyz_0[i] * pb_y + g_0_yyyzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyzz_0[i] = 3.0 * g_0_yyzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxzz_1[i] * fi_abcd_0 +
                                  g_0_yyyzz_0_xxyzz_0[i] * pb_y + g_0_yyyzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxzzz_0[i] * pb_y +
                                  g_0_yyyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyy_0[i] =
            g_0_yyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xyyyy_0[i] * pb_z + g_0_yyyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xyyyz_0[i] = 3.0 * g_0_yyzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyz_0[i] * pb_y + g_0_yyyzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyzz_0[i] = 3.0 * g_0_yyzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzz_0[i] * pb_y + g_0_yyyzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyzzz_0[i] = 3.0 * g_0_yyzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyzz_0_xyzzz_0[i] * pb_y + g_0_yyyzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xzzzz_0[i] = 3.0 * g_0_yyzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzzzz_0[i] * pb_y +
                                  g_0_yyyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyy_0[i] =
            g_0_yyyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_yyyyy_0[i] * pb_z + g_0_yyyyz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yyyyz_0[i] = 3.0 * g_0_yyzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyz_0[i] * pb_y + g_0_yyyzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyzz_0[i] = 3.0 * g_0_yyzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyzz_0[i] * pb_y + g_0_yyyzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyzzz_0[i] = 3.0 * g_0_yyzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyzzz_0[i] * pb_y + g_0_yyyzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yzzzz_0[i] = 3.0 * g_0_yyzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyzz_0_yzzzz_0[i] * pb_y + g_0_yyyzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zzzzz_0[i] = 3.0 * g_0_yyzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzzzz_0[i] * pb_y +
                                  g_0_yyyzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 504-525 components of targeted buffer : SISH

    auto g_0_yyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 504);

    auto g_0_yyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 505);

    auto g_0_yyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 506);

    auto g_0_yyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 507);

    auto g_0_yyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 508);

    auto g_0_yyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 509);

    auto g_0_yyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 510);

    auto g_0_yyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 511);

    auto g_0_yyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 512);

    auto g_0_yyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 513);

    auto g_0_yyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 514);

    auto g_0_yyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 515);

    auto g_0_yyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 516);

    auto g_0_yyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 517);

    auto g_0_yyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 518);

    auto g_0_yyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 519);

    auto g_0_yyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 520);

    auto g_0_yyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 521);

    auto g_0_yyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 522);

    auto g_0_yyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 523);

    auto g_0_yyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 524);

#pragma omp simd aligned(g_0_yyyz_0_xxxxy_0,       \
                             g_0_yyyz_0_xxxxy_1,   \
                             g_0_yyyz_0_xxxyy_0,   \
                             g_0_yyyz_0_xxxyy_1,   \
                             g_0_yyyz_0_xxyyy_0,   \
                             g_0_yyyz_0_xxyyy_1,   \
                             g_0_yyyz_0_xyyyy_0,   \
                             g_0_yyyz_0_xyyyy_1,   \
                             g_0_yyyz_0_yyyyy_0,   \
                             g_0_yyyz_0_yyyyy_1,   \
                             g_0_yyyzz_0_xxxxy_0,  \
                             g_0_yyyzz_0_xxxxy_1,  \
                             g_0_yyyzz_0_xxxyy_0,  \
                             g_0_yyyzz_0_xxxyy_1,  \
                             g_0_yyyzz_0_xxyyy_0,  \
                             g_0_yyyzz_0_xxyyy_1,  \
                             g_0_yyyzz_0_xyyyy_0,  \
                             g_0_yyyzz_0_xyyyy_1,  \
                             g_0_yyyzz_0_yyyyy_0,  \
                             g_0_yyyzz_0_yyyyy_1,  \
                             g_0_yyyzzz_0_xxxxx_0, \
                             g_0_yyyzzz_0_xxxxy_0, \
                             g_0_yyyzzz_0_xxxxz_0, \
                             g_0_yyyzzz_0_xxxyy_0, \
                             g_0_yyyzzz_0_xxxyz_0, \
                             g_0_yyyzzz_0_xxxzz_0, \
                             g_0_yyyzzz_0_xxyyy_0, \
                             g_0_yyyzzz_0_xxyyz_0, \
                             g_0_yyyzzz_0_xxyzz_0, \
                             g_0_yyyzzz_0_xxzzz_0, \
                             g_0_yyyzzz_0_xyyyy_0, \
                             g_0_yyyzzz_0_xyyyz_0, \
                             g_0_yyyzzz_0_xyyzz_0, \
                             g_0_yyyzzz_0_xyzzz_0, \
                             g_0_yyyzzz_0_xzzzz_0, \
                             g_0_yyyzzz_0_yyyyy_0, \
                             g_0_yyyzzz_0_yyyyz_0, \
                             g_0_yyyzzz_0_yyyzz_0, \
                             g_0_yyyzzz_0_yyzzz_0, \
                             g_0_yyyzzz_0_yzzzz_0, \
                             g_0_yyyzzz_0_zzzzz_0, \
                             g_0_yyzzz_0_xxxxx_0,  \
                             g_0_yyzzz_0_xxxxx_1,  \
                             g_0_yyzzz_0_xxxxz_0,  \
                             g_0_yyzzz_0_xxxxz_1,  \
                             g_0_yyzzz_0_xxxyz_0,  \
                             g_0_yyzzz_0_xxxyz_1,  \
                             g_0_yyzzz_0_xxxz_1,   \
                             g_0_yyzzz_0_xxxzz_0,  \
                             g_0_yyzzz_0_xxxzz_1,  \
                             g_0_yyzzz_0_xxyyz_0,  \
                             g_0_yyzzz_0_xxyyz_1,  \
                             g_0_yyzzz_0_xxyz_1,   \
                             g_0_yyzzz_0_xxyzz_0,  \
                             g_0_yyzzz_0_xxyzz_1,  \
                             g_0_yyzzz_0_xxzz_1,   \
                             g_0_yyzzz_0_xxzzz_0,  \
                             g_0_yyzzz_0_xxzzz_1,  \
                             g_0_yyzzz_0_xyyyz_0,  \
                             g_0_yyzzz_0_xyyyz_1,  \
                             g_0_yyzzz_0_xyyz_1,   \
                             g_0_yyzzz_0_xyyzz_0,  \
                             g_0_yyzzz_0_xyyzz_1,  \
                             g_0_yyzzz_0_xyzz_1,   \
                             g_0_yyzzz_0_xyzzz_0,  \
                             g_0_yyzzz_0_xyzzz_1,  \
                             g_0_yyzzz_0_xzzz_1,   \
                             g_0_yyzzz_0_xzzzz_0,  \
                             g_0_yyzzz_0_xzzzz_1,  \
                             g_0_yyzzz_0_yyyyz_0,  \
                             g_0_yyzzz_0_yyyyz_1,  \
                             g_0_yyzzz_0_yyyz_1,   \
                             g_0_yyzzz_0_yyyzz_0,  \
                             g_0_yyzzz_0_yyyzz_1,  \
                             g_0_yyzzz_0_yyzz_1,   \
                             g_0_yyzzz_0_yyzzz_0,  \
                             g_0_yyzzz_0_yyzzz_1,  \
                             g_0_yyzzz_0_yzzz_1,   \
                             g_0_yyzzz_0_yzzzz_0,  \
                             g_0_yyzzz_0_yzzzz_1,  \
                             g_0_yyzzz_0_zzzz_1,   \
                             g_0_yyzzz_0_zzzzz_0,  \
                             g_0_yyzzz_0_zzzzz_1,  \
                             g_0_yzzz_0_xxxxx_0,   \
                             g_0_yzzz_0_xxxxx_1,   \
                             g_0_yzzz_0_xxxxz_0,   \
                             g_0_yzzz_0_xxxxz_1,   \
                             g_0_yzzz_0_xxxyz_0,   \
                             g_0_yzzz_0_xxxyz_1,   \
                             g_0_yzzz_0_xxxzz_0,   \
                             g_0_yzzz_0_xxxzz_1,   \
                             g_0_yzzz_0_xxyyz_0,   \
                             g_0_yzzz_0_xxyyz_1,   \
                             g_0_yzzz_0_xxyzz_0,   \
                             g_0_yzzz_0_xxyzz_1,   \
                             g_0_yzzz_0_xxzzz_0,   \
                             g_0_yzzz_0_xxzzz_1,   \
                             g_0_yzzz_0_xyyyz_0,   \
                             g_0_yzzz_0_xyyyz_1,   \
                             g_0_yzzz_0_xyyzz_0,   \
                             g_0_yzzz_0_xyyzz_1,   \
                             g_0_yzzz_0_xyzzz_0,   \
                             g_0_yzzz_0_xyzzz_1,   \
                             g_0_yzzz_0_xzzzz_0,   \
                             g_0_yzzz_0_xzzzz_1,   \
                             g_0_yzzz_0_yyyyz_0,   \
                             g_0_yzzz_0_yyyyz_1,   \
                             g_0_yzzz_0_yyyzz_0,   \
                             g_0_yzzz_0_yyyzz_1,   \
                             g_0_yzzz_0_yyzzz_0,   \
                             g_0_yzzz_0_yyzzz_1,   \
                             g_0_yzzz_0_yzzzz_0,   \
                             g_0_yzzz_0_yzzzz_1,   \
                             g_0_yzzz_0_zzzzz_0,   \
                             g_0_yzzz_0_zzzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xxxxx_0[i] = 2.0 * g_0_yzzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxx_0[i] * pb_y +
                                  g_0_yyzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxy_0[i] = 2.0 * g_0_yyyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxy_0[i] * pb_z +
                                  g_0_yyyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxz_0[i] = 2.0 * g_0_yzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxz_0[i] * pb_y +
                                  g_0_yyzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyy_0[i] = 2.0 * g_0_yyyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxyy_0[i] * pb_z +
                                  g_0_yyyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxyz_0[i] = 2.0 * g_0_yzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxz_1[i] * fi_abcd_0 +
                                  g_0_yyzzz_0_xxxyz_0[i] * pb_y + g_0_yyzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxzz_0[i] = 2.0 * g_0_yzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxzz_0[i] * pb_y +
                                  g_0_yyzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxyyy_0[i] * pb_z +
                                  g_0_yyyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyz_0[i] * pb_y + g_0_yyzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxzz_1[i] * fi_abcd_0 +
                                  g_0_yyzzz_0_xxyzz_0[i] * pb_y + g_0_yyzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxzzz_0[i] = 2.0 * g_0_yzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxzzz_0[i] * pb_y +
                                  g_0_yyzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyy_0[i] = 2.0 * g_0_yyyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xyyyy_0[i] * pb_z +
                                  g_0_yyyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyz_0[i] * pb_y + g_0_yyzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzz_0[i] * pb_y + g_0_yyzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzzz_1[i] * fi_abcd_0 +
                                  g_0_yyzzz_0_xyzzz_0[i] * pb_y + g_0_yyzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xzzzz_0[i] = 2.0 * g_0_yzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzzzz_0[i] * pb_y +
                                  g_0_yyzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyy_0[i] = 2.0 * g_0_yyyz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_yyyyy_0[i] * pb_z +
                                  g_0_yyyzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yyyyz_0[i] = 2.0 * g_0_yzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyz_0[i] * pb_y + g_0_yyzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyzz_0[i] = 2.0 * g_0_yzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyzz_0[i] * pb_y + g_0_yyzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyzzz_0[i] = 2.0 * g_0_yzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyzzz_0[i] * pb_y + g_0_yyzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yzzzz_0[i] = 2.0 * g_0_yzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_yyzzz_0_yzzzz_0[i] * pb_y + g_0_yyzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zzzzz_0[i] = 2.0 * g_0_yzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzzzz_0[i] * pb_y +
                                  g_0_yyzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 525-546 components of targeted buffer : SISH

    auto g_0_yyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 525);

    auto g_0_yyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 526);

    auto g_0_yyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 527);

    auto g_0_yyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 528);

    auto g_0_yyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 529);

    auto g_0_yyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 530);

    auto g_0_yyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 531);

    auto g_0_yyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 532);

    auto g_0_yyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 533);

    auto g_0_yyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 534);

    auto g_0_yyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 535);

    auto g_0_yyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 536);

    auto g_0_yyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 537);

    auto g_0_yyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 538);

    auto g_0_yyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 539);

    auto g_0_yyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 540);

    auto g_0_yyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 541);

    auto g_0_yyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 542);

    auto g_0_yyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 543);

    auto g_0_yyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 544);

    auto g_0_yyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 545);

#pragma omp simd aligned(g_0_yyzz_0_xxxxy_0,       \
                             g_0_yyzz_0_xxxxy_1,   \
                             g_0_yyzz_0_xxxyy_0,   \
                             g_0_yyzz_0_xxxyy_1,   \
                             g_0_yyzz_0_xxyyy_0,   \
                             g_0_yyzz_0_xxyyy_1,   \
                             g_0_yyzz_0_xyyyy_0,   \
                             g_0_yyzz_0_xyyyy_1,   \
                             g_0_yyzz_0_yyyyy_0,   \
                             g_0_yyzz_0_yyyyy_1,   \
                             g_0_yyzzz_0_xxxxy_0,  \
                             g_0_yyzzz_0_xxxxy_1,  \
                             g_0_yyzzz_0_xxxyy_0,  \
                             g_0_yyzzz_0_xxxyy_1,  \
                             g_0_yyzzz_0_xxyyy_0,  \
                             g_0_yyzzz_0_xxyyy_1,  \
                             g_0_yyzzz_0_xyyyy_0,  \
                             g_0_yyzzz_0_xyyyy_1,  \
                             g_0_yyzzz_0_yyyyy_0,  \
                             g_0_yyzzz_0_yyyyy_1,  \
                             g_0_yyzzzz_0_xxxxx_0, \
                             g_0_yyzzzz_0_xxxxy_0, \
                             g_0_yyzzzz_0_xxxxz_0, \
                             g_0_yyzzzz_0_xxxyy_0, \
                             g_0_yyzzzz_0_xxxyz_0, \
                             g_0_yyzzzz_0_xxxzz_0, \
                             g_0_yyzzzz_0_xxyyy_0, \
                             g_0_yyzzzz_0_xxyyz_0, \
                             g_0_yyzzzz_0_xxyzz_0, \
                             g_0_yyzzzz_0_xxzzz_0, \
                             g_0_yyzzzz_0_xyyyy_0, \
                             g_0_yyzzzz_0_xyyyz_0, \
                             g_0_yyzzzz_0_xyyzz_0, \
                             g_0_yyzzzz_0_xyzzz_0, \
                             g_0_yyzzzz_0_xzzzz_0, \
                             g_0_yyzzzz_0_yyyyy_0, \
                             g_0_yyzzzz_0_yyyyz_0, \
                             g_0_yyzzzz_0_yyyzz_0, \
                             g_0_yyzzzz_0_yyzzz_0, \
                             g_0_yyzzzz_0_yzzzz_0, \
                             g_0_yyzzzz_0_zzzzz_0, \
                             g_0_yzzzz_0_xxxxx_0,  \
                             g_0_yzzzz_0_xxxxx_1,  \
                             g_0_yzzzz_0_xxxxz_0,  \
                             g_0_yzzzz_0_xxxxz_1,  \
                             g_0_yzzzz_0_xxxyz_0,  \
                             g_0_yzzzz_0_xxxyz_1,  \
                             g_0_yzzzz_0_xxxz_1,   \
                             g_0_yzzzz_0_xxxzz_0,  \
                             g_0_yzzzz_0_xxxzz_1,  \
                             g_0_yzzzz_0_xxyyz_0,  \
                             g_0_yzzzz_0_xxyyz_1,  \
                             g_0_yzzzz_0_xxyz_1,   \
                             g_0_yzzzz_0_xxyzz_0,  \
                             g_0_yzzzz_0_xxyzz_1,  \
                             g_0_yzzzz_0_xxzz_1,   \
                             g_0_yzzzz_0_xxzzz_0,  \
                             g_0_yzzzz_0_xxzzz_1,  \
                             g_0_yzzzz_0_xyyyz_0,  \
                             g_0_yzzzz_0_xyyyz_1,  \
                             g_0_yzzzz_0_xyyz_1,   \
                             g_0_yzzzz_0_xyyzz_0,  \
                             g_0_yzzzz_0_xyyzz_1,  \
                             g_0_yzzzz_0_xyzz_1,   \
                             g_0_yzzzz_0_xyzzz_0,  \
                             g_0_yzzzz_0_xyzzz_1,  \
                             g_0_yzzzz_0_xzzz_1,   \
                             g_0_yzzzz_0_xzzzz_0,  \
                             g_0_yzzzz_0_xzzzz_1,  \
                             g_0_yzzzz_0_yyyyz_0,  \
                             g_0_yzzzz_0_yyyyz_1,  \
                             g_0_yzzzz_0_yyyz_1,   \
                             g_0_yzzzz_0_yyyzz_0,  \
                             g_0_yzzzz_0_yyyzz_1,  \
                             g_0_yzzzz_0_yyzz_1,   \
                             g_0_yzzzz_0_yyzzz_0,  \
                             g_0_yzzzz_0_yyzzz_1,  \
                             g_0_yzzzz_0_yzzz_1,   \
                             g_0_yzzzz_0_yzzzz_0,  \
                             g_0_yzzzz_0_yzzzz_1,  \
                             g_0_yzzzz_0_zzzz_1,   \
                             g_0_yzzzz_0_zzzzz_0,  \
                             g_0_yzzzz_0_zzzzz_1,  \
                             g_0_zzzz_0_xxxxx_0,   \
                             g_0_zzzz_0_xxxxx_1,   \
                             g_0_zzzz_0_xxxxz_0,   \
                             g_0_zzzz_0_xxxxz_1,   \
                             g_0_zzzz_0_xxxyz_0,   \
                             g_0_zzzz_0_xxxyz_1,   \
                             g_0_zzzz_0_xxxzz_0,   \
                             g_0_zzzz_0_xxxzz_1,   \
                             g_0_zzzz_0_xxyyz_0,   \
                             g_0_zzzz_0_xxyyz_1,   \
                             g_0_zzzz_0_xxyzz_0,   \
                             g_0_zzzz_0_xxyzz_1,   \
                             g_0_zzzz_0_xxzzz_0,   \
                             g_0_zzzz_0_xxzzz_1,   \
                             g_0_zzzz_0_xyyyz_0,   \
                             g_0_zzzz_0_xyyyz_1,   \
                             g_0_zzzz_0_xyyzz_0,   \
                             g_0_zzzz_0_xyyzz_1,   \
                             g_0_zzzz_0_xyzzz_0,   \
                             g_0_zzzz_0_xyzzz_1,   \
                             g_0_zzzz_0_xzzzz_0,   \
                             g_0_zzzz_0_xzzzz_1,   \
                             g_0_zzzz_0_yyyyz_0,   \
                             g_0_zzzz_0_yyyyz_1,   \
                             g_0_zzzz_0_yyyzz_0,   \
                             g_0_zzzz_0_yyyzz_1,   \
                             g_0_zzzz_0_yyzzz_0,   \
                             g_0_zzzz_0_yyzzz_1,   \
                             g_0_zzzz_0_yzzzz_0,   \
                             g_0_zzzz_0_yzzzz_1,   \
                             g_0_zzzz_0_zzzzz_0,   \
                             g_0_zzzz_0_zzzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xxxxx_0[i] =
            g_0_zzzz_0_xxxxx_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxx_0[i] * pb_y + g_0_yzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxy_0[i] = 3.0 * g_0_yyzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxy_0[i] * pb_z +
                                  g_0_yyzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxz_0[i] =
            g_0_zzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxz_0[i] * pb_y + g_0_yzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyy_0[i] = 3.0 * g_0_yyzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxyy_0[i] * pb_z +
                                  g_0_yyzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxyz_0[i] = g_0_zzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_xxxyz_0[i] * pb_y + g_0_yzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxzz_0[i] =
            g_0_zzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzz_0[i] * pb_y + g_0_yzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxyyy_0[i] * pb_z +
                                  g_0_yyzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxyyz_0[i] = g_0_zzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xxyz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_xxyyz_0[i] * pb_y + g_0_yzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyzz_0[i] = g_0_zzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_xxyzz_0[i] * pb_y + g_0_yzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxzzz_0[i] =
            g_0_zzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzz_0[i] * pb_y + g_0_yzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyy_0[i] = 3.0 * g_0_yyzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xyyyy_0[i] * pb_z +
                                  g_0_yyzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xyyyz_0[i] = g_0_zzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_xyyz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_xyyyz_0[i] * pb_y + g_0_yzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyzz_0[i] = g_0_zzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xyzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_xyyzz_0[i] * pb_y + g_0_yzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyzzz_0[i] = g_0_zzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_xyzzz_0[i] * pb_y + g_0_yzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xzzzz_0[i] =
            g_0_zzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzz_0[i] * pb_y + g_0_yzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyy_0[i] = 3.0 * g_0_yyzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_yyyyy_0[i] * pb_z +
                                  g_0_yyzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yyyyz_0[i] = g_0_zzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzz_0_yyyz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_yyyyz_0[i] * pb_y + g_0_yzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyzz_0[i] = g_0_zzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_yyzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_yyyzz_0[i] * pb_y + g_0_yzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyzzz_0[i] = g_0_zzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_yzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_yyzzz_0[i] * pb_y + g_0_yzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yzzzz_0[i] = g_0_zzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzz_0_yzzzz_0[i] * pb_y + g_0_yzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zzzzz_0[i] =
            g_0_zzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzz_0[i] * pb_y + g_0_yzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 546-567 components of targeted buffer : SISH

    auto g_0_yzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 546);

    auto g_0_yzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 547);

    auto g_0_yzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 548);

    auto g_0_yzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 549);

    auto g_0_yzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 550);

    auto g_0_yzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 551);

    auto g_0_yzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 552);

    auto g_0_yzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 553);

    auto g_0_yzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 554);

    auto g_0_yzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 555);

    auto g_0_yzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 556);

    auto g_0_yzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 557);

    auto g_0_yzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 558);

    auto g_0_yzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 559);

    auto g_0_yzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 560);

    auto g_0_yzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 561);

    auto g_0_yzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 562);

    auto g_0_yzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 563);

    auto g_0_yzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 564);

    auto g_0_yzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 565);

    auto g_0_yzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 566);

#pragma omp simd aligned(g_0_yzzzzz_0_xxxxx_0,     \
                             g_0_yzzzzz_0_xxxxy_0, \
                             g_0_yzzzzz_0_xxxxz_0, \
                             g_0_yzzzzz_0_xxxyy_0, \
                             g_0_yzzzzz_0_xxxyz_0, \
                             g_0_yzzzzz_0_xxxzz_0, \
                             g_0_yzzzzz_0_xxyyy_0, \
                             g_0_yzzzzz_0_xxyyz_0, \
                             g_0_yzzzzz_0_xxyzz_0, \
                             g_0_yzzzzz_0_xxzzz_0, \
                             g_0_yzzzzz_0_xyyyy_0, \
                             g_0_yzzzzz_0_xyyyz_0, \
                             g_0_yzzzzz_0_xyyzz_0, \
                             g_0_yzzzzz_0_xyzzz_0, \
                             g_0_yzzzzz_0_xzzzz_0, \
                             g_0_yzzzzz_0_yyyyy_0, \
                             g_0_yzzzzz_0_yyyyz_0, \
                             g_0_yzzzzz_0_yyyzz_0, \
                             g_0_yzzzzz_0_yyzzz_0, \
                             g_0_yzzzzz_0_yzzzz_0, \
                             g_0_yzzzzz_0_zzzzz_0, \
                             g_0_zzzzz_0_xxxx_1,   \
                             g_0_zzzzz_0_xxxxx_0,  \
                             g_0_zzzzz_0_xxxxx_1,  \
                             g_0_zzzzz_0_xxxxy_0,  \
                             g_0_zzzzz_0_xxxxy_1,  \
                             g_0_zzzzz_0_xxxxz_0,  \
                             g_0_zzzzz_0_xxxxz_1,  \
                             g_0_zzzzz_0_xxxy_1,   \
                             g_0_zzzzz_0_xxxyy_0,  \
                             g_0_zzzzz_0_xxxyy_1,  \
                             g_0_zzzzz_0_xxxyz_0,  \
                             g_0_zzzzz_0_xxxyz_1,  \
                             g_0_zzzzz_0_xxxz_1,   \
                             g_0_zzzzz_0_xxxzz_0,  \
                             g_0_zzzzz_0_xxxzz_1,  \
                             g_0_zzzzz_0_xxyy_1,   \
                             g_0_zzzzz_0_xxyyy_0,  \
                             g_0_zzzzz_0_xxyyy_1,  \
                             g_0_zzzzz_0_xxyyz_0,  \
                             g_0_zzzzz_0_xxyyz_1,  \
                             g_0_zzzzz_0_xxyz_1,   \
                             g_0_zzzzz_0_xxyzz_0,  \
                             g_0_zzzzz_0_xxyzz_1,  \
                             g_0_zzzzz_0_xxzz_1,   \
                             g_0_zzzzz_0_xxzzz_0,  \
                             g_0_zzzzz_0_xxzzz_1,  \
                             g_0_zzzzz_0_xyyy_1,   \
                             g_0_zzzzz_0_xyyyy_0,  \
                             g_0_zzzzz_0_xyyyy_1,  \
                             g_0_zzzzz_0_xyyyz_0,  \
                             g_0_zzzzz_0_xyyyz_1,  \
                             g_0_zzzzz_0_xyyz_1,   \
                             g_0_zzzzz_0_xyyzz_0,  \
                             g_0_zzzzz_0_xyyzz_1,  \
                             g_0_zzzzz_0_xyzz_1,   \
                             g_0_zzzzz_0_xyzzz_0,  \
                             g_0_zzzzz_0_xyzzz_1,  \
                             g_0_zzzzz_0_xzzz_1,   \
                             g_0_zzzzz_0_xzzzz_0,  \
                             g_0_zzzzz_0_xzzzz_1,  \
                             g_0_zzzzz_0_yyyy_1,   \
                             g_0_zzzzz_0_yyyyy_0,  \
                             g_0_zzzzz_0_yyyyy_1,  \
                             g_0_zzzzz_0_yyyyz_0,  \
                             g_0_zzzzz_0_yyyyz_1,  \
                             g_0_zzzzz_0_yyyz_1,   \
                             g_0_zzzzz_0_yyyzz_0,  \
                             g_0_zzzzz_0_yyyzz_1,  \
                             g_0_zzzzz_0_yyzz_1,   \
                             g_0_zzzzz_0_yyzzz_0,  \
                             g_0_zzzzz_0_yyzzz_1,  \
                             g_0_zzzzz_0_yzzz_1,   \
                             g_0_zzzzz_0_yzzzz_0,  \
                             g_0_zzzzz_0_yzzzz_1,  \
                             g_0_zzzzz_0_zzzz_1,   \
                             g_0_zzzzz_0_zzzzz_0,  \
                             g_0_zzzzz_0_zzzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xxxxx_0[i] = g_0_zzzzz_0_xxxxx_0[i] * pb_y + g_0_zzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxy_0[i] = g_0_zzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxy_0[i] * pb_y + g_0_zzzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxz_0[i] = g_0_zzzzz_0_xxxxz_0[i] * pb_y + g_0_zzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyy_0[i] = 2.0 * g_0_zzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyy_0[i] * pb_y + g_0_zzzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyz_0[i] = g_0_zzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyz_0[i] * pb_y + g_0_zzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxzz_0[i] = g_0_zzzzz_0_xxxzz_0[i] * pb_y + g_0_zzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyy_0[i] = 3.0 * g_0_zzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyy_0[i] * pb_y + g_0_zzzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyz_0[i] * pb_y + g_0_zzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyzz_0[i] = g_0_zzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzz_0[i] * pb_y + g_0_zzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxzzz_0[i] = g_0_zzzzz_0_xxzzz_0[i] * pb_y + g_0_zzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyy_0[i] = 4.0 * g_0_zzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyy_0[i] * pb_y + g_0_zzzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyz_0[i] = 3.0 * g_0_zzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyz_0[i] * pb_y + g_0_zzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyzz_0[i] = 2.0 * g_0_zzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzz_0[i] * pb_y + g_0_zzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyzzz_0[i] = g_0_zzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzz_0[i] * pb_y + g_0_zzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xzzzz_0[i] = g_0_zzzzz_0_xzzzz_0[i] * pb_y + g_0_zzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyy_0[i] = 5.0 * g_0_zzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyy_0[i] * pb_y + g_0_zzzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyz_0[i] = 4.0 * g_0_zzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyz_0[i] * pb_y + g_0_zzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyzz_0[i] = 3.0 * g_0_zzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzz_0[i] * pb_y + g_0_zzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyzzz_0[i] = 2.0 * g_0_zzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzz_0[i] * pb_y + g_0_zzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yzzzz_0[i] = g_0_zzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzz_0[i] * pb_y + g_0_zzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zzzzz_0[i] = g_0_zzzzz_0_zzzzz_0[i] * pb_y + g_0_zzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 567-588 components of targeted buffer : SISH

    auto g_0_zzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 567);

    auto g_0_zzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 568);

    auto g_0_zzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 569);

    auto g_0_zzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 570);

    auto g_0_zzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 571);

    auto g_0_zzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 572);

    auto g_0_zzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 573);

    auto g_0_zzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 574);

    auto g_0_zzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 575);

    auto g_0_zzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 576);

    auto g_0_zzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 577);

    auto g_0_zzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 578);

    auto g_0_zzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 579);

    auto g_0_zzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 580);

    auto g_0_zzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 581);

    auto g_0_zzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 582);

    auto g_0_zzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 583);

    auto g_0_zzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 584);

    auto g_0_zzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 585);

    auto g_0_zzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 586);

    auto g_0_zzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 587);

#pragma omp simd aligned(g_0_zzzz_0_xxxxx_0,       \
                             g_0_zzzz_0_xxxxx_1,   \
                             g_0_zzzz_0_xxxxy_0,   \
                             g_0_zzzz_0_xxxxy_1,   \
                             g_0_zzzz_0_xxxxz_0,   \
                             g_0_zzzz_0_xxxxz_1,   \
                             g_0_zzzz_0_xxxyy_0,   \
                             g_0_zzzz_0_xxxyy_1,   \
                             g_0_zzzz_0_xxxyz_0,   \
                             g_0_zzzz_0_xxxyz_1,   \
                             g_0_zzzz_0_xxxzz_0,   \
                             g_0_zzzz_0_xxxzz_1,   \
                             g_0_zzzz_0_xxyyy_0,   \
                             g_0_zzzz_0_xxyyy_1,   \
                             g_0_zzzz_0_xxyyz_0,   \
                             g_0_zzzz_0_xxyyz_1,   \
                             g_0_zzzz_0_xxyzz_0,   \
                             g_0_zzzz_0_xxyzz_1,   \
                             g_0_zzzz_0_xxzzz_0,   \
                             g_0_zzzz_0_xxzzz_1,   \
                             g_0_zzzz_0_xyyyy_0,   \
                             g_0_zzzz_0_xyyyy_1,   \
                             g_0_zzzz_0_xyyyz_0,   \
                             g_0_zzzz_0_xyyyz_1,   \
                             g_0_zzzz_0_xyyzz_0,   \
                             g_0_zzzz_0_xyyzz_1,   \
                             g_0_zzzz_0_xyzzz_0,   \
                             g_0_zzzz_0_xyzzz_1,   \
                             g_0_zzzz_0_xzzzz_0,   \
                             g_0_zzzz_0_xzzzz_1,   \
                             g_0_zzzz_0_yyyyy_0,   \
                             g_0_zzzz_0_yyyyy_1,   \
                             g_0_zzzz_0_yyyyz_0,   \
                             g_0_zzzz_0_yyyyz_1,   \
                             g_0_zzzz_0_yyyzz_0,   \
                             g_0_zzzz_0_yyyzz_1,   \
                             g_0_zzzz_0_yyzzz_0,   \
                             g_0_zzzz_0_yyzzz_1,   \
                             g_0_zzzz_0_yzzzz_0,   \
                             g_0_zzzz_0_yzzzz_1,   \
                             g_0_zzzz_0_zzzzz_0,   \
                             g_0_zzzz_0_zzzzz_1,   \
                             g_0_zzzzz_0_xxxx_1,   \
                             g_0_zzzzz_0_xxxxx_0,  \
                             g_0_zzzzz_0_xxxxx_1,  \
                             g_0_zzzzz_0_xxxxy_0,  \
                             g_0_zzzzz_0_xxxxy_1,  \
                             g_0_zzzzz_0_xxxxz_0,  \
                             g_0_zzzzz_0_xxxxz_1,  \
                             g_0_zzzzz_0_xxxy_1,   \
                             g_0_zzzzz_0_xxxyy_0,  \
                             g_0_zzzzz_0_xxxyy_1,  \
                             g_0_zzzzz_0_xxxyz_0,  \
                             g_0_zzzzz_0_xxxyz_1,  \
                             g_0_zzzzz_0_xxxz_1,   \
                             g_0_zzzzz_0_xxxzz_0,  \
                             g_0_zzzzz_0_xxxzz_1,  \
                             g_0_zzzzz_0_xxyy_1,   \
                             g_0_zzzzz_0_xxyyy_0,  \
                             g_0_zzzzz_0_xxyyy_1,  \
                             g_0_zzzzz_0_xxyyz_0,  \
                             g_0_zzzzz_0_xxyyz_1,  \
                             g_0_zzzzz_0_xxyz_1,   \
                             g_0_zzzzz_0_xxyzz_0,  \
                             g_0_zzzzz_0_xxyzz_1,  \
                             g_0_zzzzz_0_xxzz_1,   \
                             g_0_zzzzz_0_xxzzz_0,  \
                             g_0_zzzzz_0_xxzzz_1,  \
                             g_0_zzzzz_0_xyyy_1,   \
                             g_0_zzzzz_0_xyyyy_0,  \
                             g_0_zzzzz_0_xyyyy_1,  \
                             g_0_zzzzz_0_xyyyz_0,  \
                             g_0_zzzzz_0_xyyyz_1,  \
                             g_0_zzzzz_0_xyyz_1,   \
                             g_0_zzzzz_0_xyyzz_0,  \
                             g_0_zzzzz_0_xyyzz_1,  \
                             g_0_zzzzz_0_xyzz_1,   \
                             g_0_zzzzz_0_xyzzz_0,  \
                             g_0_zzzzz_0_xyzzz_1,  \
                             g_0_zzzzz_0_xzzz_1,   \
                             g_0_zzzzz_0_xzzzz_0,  \
                             g_0_zzzzz_0_xzzzz_1,  \
                             g_0_zzzzz_0_yyyy_1,   \
                             g_0_zzzzz_0_yyyyy_0,  \
                             g_0_zzzzz_0_yyyyy_1,  \
                             g_0_zzzzz_0_yyyyz_0,  \
                             g_0_zzzzz_0_yyyyz_1,  \
                             g_0_zzzzz_0_yyyz_1,   \
                             g_0_zzzzz_0_yyyzz_0,  \
                             g_0_zzzzz_0_yyyzz_1,  \
                             g_0_zzzzz_0_yyzz_1,   \
                             g_0_zzzzz_0_yyzzz_0,  \
                             g_0_zzzzz_0_yyzzz_1,  \
                             g_0_zzzzz_0_yzzz_1,   \
                             g_0_zzzzz_0_yzzzz_0,  \
                             g_0_zzzzz_0_yzzzz_1,  \
                             g_0_zzzzz_0_zzzz_1,   \
                             g_0_zzzzz_0_zzzzz_0,  \
                             g_0_zzzzz_0_zzzzz_1,  \
                             g_0_zzzzzz_0_xxxxx_0, \
                             g_0_zzzzzz_0_xxxxy_0, \
                             g_0_zzzzzz_0_xxxxz_0, \
                             g_0_zzzzzz_0_xxxyy_0, \
                             g_0_zzzzzz_0_xxxyz_0, \
                             g_0_zzzzzz_0_xxxzz_0, \
                             g_0_zzzzzz_0_xxyyy_0, \
                             g_0_zzzzzz_0_xxyyz_0, \
                             g_0_zzzzzz_0_xxyzz_0, \
                             g_0_zzzzzz_0_xxzzz_0, \
                             g_0_zzzzzz_0_xyyyy_0, \
                             g_0_zzzzzz_0_xyyyz_0, \
                             g_0_zzzzzz_0_xyyzz_0, \
                             g_0_zzzzzz_0_xyzzz_0, \
                             g_0_zzzzzz_0_xzzzz_0, \
                             g_0_zzzzzz_0_yyyyy_0, \
                             g_0_zzzzzz_0_yyyyz_0, \
                             g_0_zzzzzz_0_yyyzz_0, \
                             g_0_zzzzzz_0_yyzzz_0, \
                             g_0_zzzzzz_0_yzzzz_0, \
                             g_0_zzzzzz_0_zzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xxxxx_0[i] = 5.0 * g_0_zzzz_0_xxxxx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxx_0[i] * pb_z +
                                  g_0_zzzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxy_0[i] = 5.0 * g_0_zzzz_0_xxxxy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxy_0[i] * pb_z +
                                  g_0_zzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxz_0[i] = 5.0 * g_0_zzzz_0_xxxxz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxx_1[i] * fi_abcd_0 +
                                  g_0_zzzzz_0_xxxxz_0[i] * pb_z + g_0_zzzzz_0_xxxxz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyy_0[i] = 5.0 * g_0_zzzz_0_xxxyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxyy_0[i] * pb_z +
                                  g_0_zzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyz_0[i] = 5.0 * g_0_zzzz_0_xxxyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxy_1[i] * fi_abcd_0 +
                                  g_0_zzzzz_0_xxxyz_0[i] * pb_z + g_0_zzzzz_0_xxxyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxzz_0[i] = 5.0 * g_0_zzzz_0_xxxzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzz_0[i] * pb_z + g_0_zzzzz_0_xxxzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyy_0[i] = 5.0 * g_0_zzzz_0_xxyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxyyy_0[i] * pb_z +
                                  g_0_zzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyz_0[i] = 5.0 * g_0_zzzz_0_xxyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxyy_1[i] * fi_abcd_0 +
                                  g_0_zzzzz_0_xxyyz_0[i] * pb_z + g_0_zzzzz_0_xxyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyzz_0[i] = 5.0 * g_0_zzzz_0_xxyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzz_0[i] * pb_z + g_0_zzzzz_0_xxyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxzzz_0[i] = 5.0 * g_0_zzzz_0_xxzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzz_0[i] * pb_z + g_0_zzzzz_0_xxzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyy_0[i] = 5.0 * g_0_zzzz_0_xyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xyyyy_0[i] * pb_z +
                                  g_0_zzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyz_0[i] = 5.0 * g_0_zzzz_0_xyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xyyy_1[i] * fi_abcd_0 +
                                  g_0_zzzzz_0_xyyyz_0[i] * pb_z + g_0_zzzzz_0_xyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyzz_0[i] = 5.0 * g_0_zzzz_0_xyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzz_0[i] * pb_z + g_0_zzzzz_0_xyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyzzz_0[i] = 5.0 * g_0_zzzz_0_xyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzz_0[i] * pb_z + g_0_zzzzz_0_xyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xzzzz_0[i] = 5.0 * g_0_zzzz_0_xzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzz_0[i] * pb_z + g_0_zzzzz_0_xzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyy_0[i] = 5.0 * g_0_zzzz_0_yyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_yyyyy_0[i] * pb_z +
                                  g_0_zzzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyz_0[i] = 5.0 * g_0_zzzz_0_yyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_yyyy_1[i] * fi_abcd_0 +
                                  g_0_zzzzz_0_yyyyz_0[i] * pb_z + g_0_zzzzz_0_yyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyzz_0[i] = 5.0 * g_0_zzzz_0_yyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzz_0[i] * pb_z + g_0_zzzzz_0_yyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyzzz_0[i] = 5.0 * g_0_zzzz_0_yyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzz_0[i] * pb_z + g_0_zzzzz_0_yyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yzzzz_0[i] = 5.0 * g_0_zzzz_0_yzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzz_0[i] * pb_z + g_0_zzzzz_0_yzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zzzzz_0[i] = 5.0 * g_0_zzzz_0_zzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_zzzzz_0[i] * pb_z + g_0_zzzzz_0_zzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
