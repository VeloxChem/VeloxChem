#include "ElectronRepulsionPrimRecSKSH.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sksh(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sksh,
                                  size_t idx_eri_0_shsh,
                                  size_t idx_eri_1_shsh,
                                  size_t idx_eri_1_sisg,
                                  size_t idx_eri_0_sish,
                                  size_t idx_eri_1_sish,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
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

    auto g_0_xxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 23);

    auto g_0_xxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 26);

    auto g_0_xxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 30);

    auto g_0_xxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 35);

    auto g_0_xxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 42);

    auto g_0_xxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 43);

    auto g_0_xxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 45);

    auto g_0_xxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 48);

    auto g_0_xxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 52);

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

    auto g_0_yyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 339);

    auto g_0_yyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 342);

    auto g_0_yyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 346);

    auto g_0_yyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 351);

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

    auto g_0_yzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 401);

    auto g_0_yzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 403);

    auto g_0_yzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 404);

    auto g_0_yzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 406);

    auto g_0_yzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 407);

    auto g_0_yzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 408);

    auto g_0_yzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 410);

    auto g_0_yzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 411);

    auto g_0_yzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 412);

    auto g_0_yzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 413);

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

    auto g_0_xxxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 23);

    auto g_0_xxxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 26);

    auto g_0_xxxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 30);

    auto g_0_xxxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 35);

    auto g_0_xxxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 42);

    auto g_0_xxxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 43);

    auto g_0_xxxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 45);

    auto g_0_xxxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 48);

    auto g_0_xxxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 52);

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

    auto g_0_yyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 339);

    auto g_0_yyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 342);

    auto g_0_yyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 346);

    auto g_0_yyyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 351);

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

    auto g_0_yzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 401);

    auto g_0_yzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 403);

    auto g_0_yzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 404);

    auto g_0_yzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 406);

    auto g_0_yzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 407);

    auto g_0_yzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 408);

    auto g_0_yzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 410);

    auto g_0_yzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 411);

    auto g_0_yzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 412);

    auto g_0_yzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 413);

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

    /// Set up components of auxilary buffer : SISG

    auto g_0_xxxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg);

    auto g_0_xxxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 1);

    auto g_0_xxxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 2);

    auto g_0_xxxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 3);

    auto g_0_xxxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 4);

    auto g_0_xxxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 5);

    auto g_0_xxxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 6);

    auto g_0_xxxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 7);

    auto g_0_xxxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 8);

    auto g_0_xxxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 9);

    auto g_0_xxxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 10);

    auto g_0_xxxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 11);

    auto g_0_xxxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 12);

    auto g_0_xxxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 13);

    auto g_0_xxxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 14);

    auto g_0_xxxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 32);

    auto g_0_xxxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 34);

    auto g_0_xxxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 35);

    auto g_0_xxxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 37);

    auto g_0_xxxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 38);

    auto g_0_xxxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 39);

    auto g_0_xxxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 41);

    auto g_0_xxxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 42);

    auto g_0_xxxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 43);

    auto g_0_xxxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 44);

    auto g_0_xxxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 45);

    auto g_0_xxxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 46);

    auto g_0_xxxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 47);

    auto g_0_xxxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 48);

    auto g_0_xxxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 49);

    auto g_0_xxxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 50);

    auto g_0_xxxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 51);

    auto g_0_xxxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 52);

    auto g_0_xxxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 53);

    auto g_0_xxxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 54);

    auto g_0_xxxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 55);

    auto g_0_xxxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 56);

    auto g_0_xxxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 57);

    auto g_0_xxxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 58);

    auto g_0_xxxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 59);

    auto g_0_xxxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 75);

    auto g_0_xxxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 76);

    auto g_0_xxxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 77);

    auto g_0_xxxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 78);

    auto g_0_xxxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 79);

    auto g_0_xxxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 80);

    auto g_0_xxxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 81);

    auto g_0_xxxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 82);

    auto g_0_xxxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 83);

    auto g_0_xxxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 84);

    auto g_0_xxxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 85);

    auto g_0_xxxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 86);

    auto g_0_xxxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 87);

    auto g_0_xxxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 88);

    auto g_0_xxxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 89);

    auto g_0_xxxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 90);

    auto g_0_xxxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 91);

    auto g_0_xxxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 92);

    auto g_0_xxxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 93);

    auto g_0_xxxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 94);

    auto g_0_xxxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 95);

    auto g_0_xxxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 96);

    auto g_0_xxxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 97);

    auto g_0_xxxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 98);

    auto g_0_xxxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 99);

    auto g_0_xxxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 100);

    auto g_0_xxxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 101);

    auto g_0_xxxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 102);

    auto g_0_xxxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 103);

    auto g_0_xxxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 104);

    auto g_0_xxxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 135);

    auto g_0_xxxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 136);

    auto g_0_xxxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 137);

    auto g_0_xxxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 138);

    auto g_0_xxxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 139);

    auto g_0_xxxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 140);

    auto g_0_xxxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 141);

    auto g_0_xxxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 142);

    auto g_0_xxxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 143);

    auto g_0_xxxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 144);

    auto g_0_xxxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 145);

    auto g_0_xxxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 146);

    auto g_0_xxxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 147);

    auto g_0_xxxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 148);

    auto g_0_xxxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 149);

    auto g_0_xxyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 150);

    auto g_0_xxyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 151);

    auto g_0_xxyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 152);

    auto g_0_xxyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 153);

    auto g_0_xxyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 154);

    auto g_0_xxyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 155);

    auto g_0_xxyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 156);

    auto g_0_xxyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 157);

    auto g_0_xxyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 158);

    auto g_0_xxyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 159);

    auto g_0_xxyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 160);

    auto g_0_xxyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 161);

    auto g_0_xxyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 162);

    auto g_0_xxyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 163);

    auto g_0_xxyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 164);

    auto g_0_xxyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 184);

    auto g_0_xxyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 187);

    auto g_0_xxyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 188);

    auto g_0_xxyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 191);

    auto g_0_xxyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 192);

    auto g_0_xxyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 193);

    auto g_0_xxzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 210);

    auto g_0_xxzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 211);

    auto g_0_xxzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 212);

    auto g_0_xxzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 213);

    auto g_0_xxzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 214);

    auto g_0_xxzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 215);

    auto g_0_xxzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 216);

    auto g_0_xxzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 217);

    auto g_0_xxzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 218);

    auto g_0_xxzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 219);

    auto g_0_xxzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 220);

    auto g_0_xxzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 221);

    auto g_0_xxzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 222);

    auto g_0_xxzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 223);

    auto g_0_xxzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 224);

    auto g_0_xyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 226);

    auto g_0_xyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 228);

    auto g_0_xyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 229);

    auto g_0_xyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 231);

    auto g_0_xyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 232);

    auto g_0_xyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 233);

    auto g_0_xyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 235);

    auto g_0_xyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 236);

    auto g_0_xyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 237);

    auto g_0_xyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 238);

    auto g_0_xyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 259);

    auto g_0_xyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 262);

    auto g_0_xyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 263);

    auto g_0_xyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 266);

    auto g_0_xyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 267);

    auto g_0_xyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 268);

    auto g_0_xyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 274);

    auto g_0_xyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 277);

    auto g_0_xyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 278);

    auto g_0_xyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 281);

    auto g_0_xyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 282);

    auto g_0_xyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 283);

    auto g_0_xzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 302);

    auto g_0_xzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 304);

    auto g_0_xzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 305);

    auto g_0_xzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 307);

    auto g_0_xzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 308);

    auto g_0_xzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 309);

    auto g_0_xzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 311);

    auto g_0_xzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 312);

    auto g_0_xzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 313);

    auto g_0_xzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 314);

    auto g_0_yyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 315);

    auto g_0_yyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 316);

    auto g_0_yyyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 317);

    auto g_0_yyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 318);

    auto g_0_yyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 319);

    auto g_0_yyyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 320);

    auto g_0_yyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 321);

    auto g_0_yyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 322);

    auto g_0_yyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 323);

    auto g_0_yyyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 324);

    auto g_0_yyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 325);

    auto g_0_yyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 326);

    auto g_0_yyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 327);

    auto g_0_yyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 328);

    auto g_0_yyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 329);

    auto g_0_yyyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 332);

    auto g_0_yyyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 334);

    auto g_0_yyyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 335);

    auto g_0_yyyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 337);

    auto g_0_yyyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 338);

    auto g_0_yyyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 339);

    auto g_0_yyyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 341);

    auto g_0_yyyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 342);

    auto g_0_yyyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 343);

    auto g_0_yyyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 344);

    auto g_0_yyyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 345);

    auto g_0_yyyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 346);

    auto g_0_yyyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 347);

    auto g_0_yyyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 348);

    auto g_0_yyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 349);

    auto g_0_yyyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 350);

    auto g_0_yyyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 351);

    auto g_0_yyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 352);

    auto g_0_yyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 353);

    auto g_0_yyyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 354);

    auto g_0_yyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 355);

    auto g_0_yyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 356);

    auto g_0_yyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 357);

    auto g_0_yyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 358);

    auto g_0_yyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 359);

    auto g_0_yyyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 360);

    auto g_0_yyyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 361);

    auto g_0_yyyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 362);

    auto g_0_yyyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 363);

    auto g_0_yyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 364);

    auto g_0_yyyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 365);

    auto g_0_yyyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 366);

    auto g_0_yyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 367);

    auto g_0_yyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 368);

    auto g_0_yyyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 369);

    auto g_0_yyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 370);

    auto g_0_yyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 371);

    auto g_0_yyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 372);

    auto g_0_yyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 373);

    auto g_0_yyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 374);

    auto g_0_yyzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 375);

    auto g_0_yyzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 376);

    auto g_0_yyzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 377);

    auto g_0_yyzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 378);

    auto g_0_yyzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 379);

    auto g_0_yyzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 380);

    auto g_0_yyzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 381);

    auto g_0_yyzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 382);

    auto g_0_yyzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 383);

    auto g_0_yyzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 384);

    auto g_0_yyzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 385);

    auto g_0_yyzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 386);

    auto g_0_yyzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 387);

    auto g_0_yyzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 388);

    auto g_0_yyzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 389);

    auto g_0_yzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 391);

    auto g_0_yzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 392);

    auto g_0_yzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 393);

    auto g_0_yzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 394);

    auto g_0_yzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 395);

    auto g_0_yzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 396);

    auto g_0_yzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 397);

    auto g_0_yzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 398);

    auto g_0_yzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 399);

    auto g_0_yzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 400);

    auto g_0_yzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 401);

    auto g_0_yzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 402);

    auto g_0_yzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 403);

    auto g_0_yzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 404);

    auto g_0_zzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 405);

    auto g_0_zzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 406);

    auto g_0_zzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 407);

    auto g_0_zzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 408);

    auto g_0_zzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 409);

    auto g_0_zzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 410);

    auto g_0_zzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 411);

    auto g_0_zzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 412);

    auto g_0_zzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 413);

    auto g_0_zzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 414);

    auto g_0_zzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 415);

    auto g_0_zzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 416);

    auto g_0_zzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 417);

    auto g_0_zzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 418);

    auto g_0_zzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 419);

    /// Set up components of auxilary buffer : SISH

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

    auto g_0_xxxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 21);

    auto g_0_xxxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 22);

    auto g_0_xxxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 23);

    auto g_0_xxxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 24);

    auto g_0_xxxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 26);

    auto g_0_xxxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 27);

    auto g_0_xxxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 30);

    auto g_0_xxxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 31);

    auto g_0_xxxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 35);

    auto g_0_xxxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 36);

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

    auto g_0_xxxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 58);

    auto g_0_xxxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 59);

    auto g_0_xxxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 60);

    auto g_0_xxxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 61);

    auto g_0_xxxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 62);

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

    auto g_0_xxxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 148);

    auto g_0_xxxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 150);

    auto g_0_xxxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 153);

    auto g_0_xxxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 157);

    auto g_0_xxxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 168);

    auto g_0_xxxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 170);

    auto g_0_xxxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 173);

    auto g_0_xxxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 177);

    auto g_0_xxxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 182);

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

    auto g_0_xxyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 232);

    auto g_0_xxyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 234);

    auto g_0_xxyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 237);

    auto g_0_xxyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 241);

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

    auto g_0_xxyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 273);

    auto g_0_xxyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 275);

    auto g_0_xxyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 278);

    auto g_0_xxyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 282);

    auto g_0_xxyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sish + 287);

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

    auto g_0_xyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 315);

    auto g_0_xyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sish + 316);

    auto g_0_xyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sish + 318);

    auto g_0_xyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 319);

    auto g_0_xyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sish + 321);

    auto g_0_xyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 322);

    auto g_0_xyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 323);

    auto g_0_xyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sish + 325);

    auto g_0_xyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 326);

    auto g_0_xyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 327);

    auto g_0_xyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 328);

    auto g_0_xyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 330);

    auto g_0_xyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 331);

    auto g_0_xyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 332);

    auto g_0_xyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 333);

    auto g_0_xyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 334);

    auto g_0_xyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 335);

    auto g_0_xyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 361);

    auto g_0_xyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 364);

    auto g_0_xyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 365);

    auto g_0_xyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 368);

    auto g_0_xyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 369);

    auto g_0_xyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 370);

    auto g_0_xyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 372);

    auto g_0_xyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 373);

    auto g_0_xyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 374);

    auto g_0_xyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 375);

    auto g_0_xyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 376);

    auto g_0_xyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 377);

    auto g_0_xyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 382);

    auto g_0_xyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 385);

    auto g_0_xyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 386);

    auto g_0_xyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sish + 389);

    auto g_0_xyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sish + 390);

    auto g_0_xyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sish + 391);

    auto g_0_xyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sish + 393);

    auto g_0_xyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sish + 394);

    auto g_0_xyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sish + 395);

    auto g_0_xyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sish + 396);

    auto g_0_xyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sish + 397);

    auto g_0_xyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sish + 398);

    auto g_0_xzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sish + 420);

    auto g_0_xzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sish + 422);

    auto g_0_xzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sish + 424);

    auto g_0_xzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sish + 425);

    auto g_0_xzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sish + 427);

    auto g_0_xzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sish + 428);

    auto g_0_xzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sish + 429);

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

    /// Set up components of auxilary buffer : SISH

    auto g_0_xxxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish);

    auto g_0_xxxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 1);

    auto g_0_xxxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 2);

    auto g_0_xxxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 3);

    auto g_0_xxxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 4);

    auto g_0_xxxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 5);

    auto g_0_xxxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 6);

    auto g_0_xxxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 7);

    auto g_0_xxxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 8);

    auto g_0_xxxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 9);

    auto g_0_xxxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 10);

    auto g_0_xxxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 11);

    auto g_0_xxxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 12);

    auto g_0_xxxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 13);

    auto g_0_xxxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 14);

    auto g_0_xxxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 15);

    auto g_0_xxxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 16);

    auto g_0_xxxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 17);

    auto g_0_xxxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 18);

    auto g_0_xxxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 19);

    auto g_0_xxxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 20);

    auto g_0_xxxxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 21);

    auto g_0_xxxxxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 22);

    auto g_0_xxxxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 23);

    auto g_0_xxxxxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 24);

    auto g_0_xxxxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 26);

    auto g_0_xxxxxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 27);

    auto g_0_xxxxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 30);

    auto g_0_xxxxxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 31);

    auto g_0_xxxxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 35);

    auto g_0_xxxxxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 36);

    auto g_0_xxxxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 42);

    auto g_0_xxxxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 43);

    auto g_0_xxxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 44);

    auto g_0_xxxxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 45);

    auto g_0_xxxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 46);

    auto g_0_xxxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 47);

    auto g_0_xxxxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 48);

    auto g_0_xxxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 49);

    auto g_0_xxxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 50);

    auto g_0_xxxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 51);

    auto g_0_xxxxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 52);

    auto g_0_xxxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 53);

    auto g_0_xxxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 54);

    auto g_0_xxxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 55);

    auto g_0_xxxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 56);

    auto g_0_xxxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 58);

    auto g_0_xxxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 59);

    auto g_0_xxxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 60);

    auto g_0_xxxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 61);

    auto g_0_xxxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 62);

    auto g_0_xxxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 63);

    auto g_0_xxxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 64);

    auto g_0_xxxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 65);

    auto g_0_xxxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 66);

    auto g_0_xxxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 67);

    auto g_0_xxxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 68);

    auto g_0_xxxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 69);

    auto g_0_xxxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 70);

    auto g_0_xxxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 71);

    auto g_0_xxxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 72);

    auto g_0_xxxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 73);

    auto g_0_xxxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 74);

    auto g_0_xxxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 75);

    auto g_0_xxxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 76);

    auto g_0_xxxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 77);

    auto g_0_xxxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 78);

    auto g_0_xxxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 79);

    auto g_0_xxxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 80);

    auto g_0_xxxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 81);

    auto g_0_xxxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 82);

    auto g_0_xxxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 83);

    auto g_0_xxxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 105);

    auto g_0_xxxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 106);

    auto g_0_xxxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 107);

    auto g_0_xxxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 108);

    auto g_0_xxxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 109);

    auto g_0_xxxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 110);

    auto g_0_xxxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 111);

    auto g_0_xxxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 112);

    auto g_0_xxxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 113);

    auto g_0_xxxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 114);

    auto g_0_xxxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 115);

    auto g_0_xxxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 116);

    auto g_0_xxxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 117);

    auto g_0_xxxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 118);

    auto g_0_xxxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 119);

    auto g_0_xxxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 120);

    auto g_0_xxxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 121);

    auto g_0_xxxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 122);

    auto g_0_xxxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 123);

    auto g_0_xxxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 124);

    auto g_0_xxxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 125);

    auto g_0_xxxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 126);

    auto g_0_xxxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 127);

    auto g_0_xxxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 128);

    auto g_0_xxxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 129);

    auto g_0_xxxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 130);

    auto g_0_xxxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 131);

    auto g_0_xxxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 132);

    auto g_0_xxxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 133);

    auto g_0_xxxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 134);

    auto g_0_xxxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 135);

    auto g_0_xxxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 136);

    auto g_0_xxxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 137);

    auto g_0_xxxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 138);

    auto g_0_xxxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 139);

    auto g_0_xxxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 140);

    auto g_0_xxxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 141);

    auto g_0_xxxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 142);

    auto g_0_xxxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 143);

    auto g_0_xxxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 144);

    auto g_0_xxxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 145);

    auto g_0_xxxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 146);

    auto g_0_xxxyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 148);

    auto g_0_xxxyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 150);

    auto g_0_xxxyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 153);

    auto g_0_xxxyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 157);

    auto g_0_xxxyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 168);

    auto g_0_xxxyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 170);

    auto g_0_xxxyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 173);

    auto g_0_xxxyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 177);

    auto g_0_xxxyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 182);

    auto g_0_xxxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 189);

    auto g_0_xxxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 190);

    auto g_0_xxxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 191);

    auto g_0_xxxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 192);

    auto g_0_xxxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 193);

    auto g_0_xxxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 194);

    auto g_0_xxxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 195);

    auto g_0_xxxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 196);

    auto g_0_xxxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 197);

    auto g_0_xxxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 198);

    auto g_0_xxxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 199);

    auto g_0_xxxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 200);

    auto g_0_xxxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 201);

    auto g_0_xxxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 202);

    auto g_0_xxxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 203);

    auto g_0_xxxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 204);

    auto g_0_xxxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 205);

    auto g_0_xxxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 206);

    auto g_0_xxxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 207);

    auto g_0_xxxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 208);

    auto g_0_xxxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 209);

    auto g_0_xxyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 210);

    auto g_0_xxyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 211);

    auto g_0_xxyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 212);

    auto g_0_xxyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 213);

    auto g_0_xxyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 214);

    auto g_0_xxyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 215);

    auto g_0_xxyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 216);

    auto g_0_xxyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 217);

    auto g_0_xxyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 218);

    auto g_0_xxyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 219);

    auto g_0_xxyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 220);

    auto g_0_xxyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 221);

    auto g_0_xxyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 222);

    auto g_0_xxyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 223);

    auto g_0_xxyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 224);

    auto g_0_xxyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 225);

    auto g_0_xxyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 226);

    auto g_0_xxyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 227);

    auto g_0_xxyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 228);

    auto g_0_xxyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 229);

    auto g_0_xxyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 230);

    auto g_0_xxyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 232);

    auto g_0_xxyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 234);

    auto g_0_xxyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 237);

    auto g_0_xxyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 241);

    auto g_0_xxyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 252);

    auto g_0_xxyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 253);

    auto g_0_xxyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 254);

    auto g_0_xxyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 255);

    auto g_0_xxyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 256);

    auto g_0_xxyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 257);

    auto g_0_xxyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 258);

    auto g_0_xxyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 259);

    auto g_0_xxyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 260);

    auto g_0_xxyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 261);

    auto g_0_xxyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 262);

    auto g_0_xxyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 263);

    auto g_0_xxyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 264);

    auto g_0_xxyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 265);

    auto g_0_xxyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 266);

    auto g_0_xxyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 267);

    auto g_0_xxyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 268);

    auto g_0_xxyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 269);

    auto g_0_xxyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 270);

    auto g_0_xxyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 271);

    auto g_0_xxyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 272);

    auto g_0_xxyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 273);

    auto g_0_xxyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 275);

    auto g_0_xxyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 278);

    auto g_0_xxyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 282);

    auto g_0_xxyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 287);

    auto g_0_xxzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 294);

    auto g_0_xxzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 295);

    auto g_0_xxzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 296);

    auto g_0_xxzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 297);

    auto g_0_xxzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 298);

    auto g_0_xxzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 299);

    auto g_0_xxzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 300);

    auto g_0_xxzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 301);

    auto g_0_xxzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 302);

    auto g_0_xxzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 303);

    auto g_0_xxzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 304);

    auto g_0_xxzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 305);

    auto g_0_xxzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 306);

    auto g_0_xxzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 307);

    auto g_0_xxzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 308);

    auto g_0_xxzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 309);

    auto g_0_xxzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 310);

    auto g_0_xxzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 311);

    auto g_0_xxzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 312);

    auto g_0_xxzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 313);

    auto g_0_xxzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 314);

    auto g_0_xyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 315);

    auto g_0_xyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 316);

    auto g_0_xyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 318);

    auto g_0_xyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 319);

    auto g_0_xyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 321);

    auto g_0_xyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 322);

    auto g_0_xyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 323);

    auto g_0_xyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 325);

    auto g_0_xyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 326);

    auto g_0_xyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 327);

    auto g_0_xyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 328);

    auto g_0_xyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 330);

    auto g_0_xyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 331);

    auto g_0_xyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 332);

    auto g_0_xyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 333);

    auto g_0_xyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 334);

    auto g_0_xyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 335);

    auto g_0_xyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 361);

    auto g_0_xyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 364);

    auto g_0_xyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 365);

    auto g_0_xyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 368);

    auto g_0_xyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 369);

    auto g_0_xyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 370);

    auto g_0_xyyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 372);

    auto g_0_xyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 373);

    auto g_0_xyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 374);

    auto g_0_xyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 375);

    auto g_0_xyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 376);

    auto g_0_xyyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 377);

    auto g_0_xyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 382);

    auto g_0_xyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 385);

    auto g_0_xyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 386);

    auto g_0_xyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 389);

    auto g_0_xyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 390);

    auto g_0_xyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 391);

    auto g_0_xyyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 393);

    auto g_0_xyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 394);

    auto g_0_xyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 395);

    auto g_0_xyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 396);

    auto g_0_xyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 397);

    auto g_0_xyyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 398);

    auto g_0_xzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 420);

    auto g_0_xzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 422);

    auto g_0_xzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 424);

    auto g_0_xzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 425);

    auto g_0_xzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 427);

    auto g_0_xzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 428);

    auto g_0_xzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 429);

    auto g_0_xzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 431);

    auto g_0_xzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 432);

    auto g_0_xzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 433);

    auto g_0_xzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 434);

    auto g_0_xzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 435);

    auto g_0_xzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 436);

    auto g_0_xzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 437);

    auto g_0_xzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 438);

    auto g_0_xzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 439);

    auto g_0_xzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 440);

    auto g_0_yyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 441);

    auto g_0_yyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 442);

    auto g_0_yyyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 443);

    auto g_0_yyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 444);

    auto g_0_yyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 445);

    auto g_0_yyyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 446);

    auto g_0_yyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 447);

    auto g_0_yyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 448);

    auto g_0_yyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 449);

    auto g_0_yyyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 450);

    auto g_0_yyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 451);

    auto g_0_yyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 452);

    auto g_0_yyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 453);

    auto g_0_yyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 454);

    auto g_0_yyyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 455);

    auto g_0_yyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 456);

    auto g_0_yyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 457);

    auto g_0_yyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 458);

    auto g_0_yyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 459);

    auto g_0_yyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 460);

    auto g_0_yyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 461);

    auto g_0_yyyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 463);

    auto g_0_yyyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 464);

    auto g_0_yyyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 465);

    auto g_0_yyyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 466);

    auto g_0_yyyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 467);

    auto g_0_yyyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 468);

    auto g_0_yyyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 469);

    auto g_0_yyyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 470);

    auto g_0_yyyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 471);

    auto g_0_yyyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 472);

    auto g_0_yyyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 473);

    auto g_0_yyyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 474);

    auto g_0_yyyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 475);

    auto g_0_yyyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 476);

    auto g_0_yyyyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 477);

    auto g_0_yyyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 478);

    auto g_0_yyyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 479);

    auto g_0_yyyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 480);

    auto g_0_yyyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 481);

    auto g_0_yyyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 482);

    auto g_0_yyyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 483);

    auto g_0_yyyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 484);

    auto g_0_yyyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 485);

    auto g_0_yyyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 486);

    auto g_0_yyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 487);

    auto g_0_yyyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 488);

    auto g_0_yyyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 489);

    auto g_0_yyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 490);

    auto g_0_yyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 491);

    auto g_0_yyyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 492);

    auto g_0_yyyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 493);

    auto g_0_yyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 494);

    auto g_0_yyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 495);

    auto g_0_yyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 496);

    auto g_0_yyyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 497);

    auto g_0_yyyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 498);

    auto g_0_yyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 499);

    auto g_0_yyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 500);

    auto g_0_yyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 501);

    auto g_0_yyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 502);

    auto g_0_yyyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 503);

    auto g_0_yyyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 504);

    auto g_0_yyyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 505);

    auto g_0_yyyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 506);

    auto g_0_yyyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 507);

    auto g_0_yyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 508);

    auto g_0_yyyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 509);

    auto g_0_yyyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 510);

    auto g_0_yyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 511);

    auto g_0_yyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 512);

    auto g_0_yyyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 513);

    auto g_0_yyyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 514);

    auto g_0_yyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 515);

    auto g_0_yyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 516);

    auto g_0_yyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 517);

    auto g_0_yyyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 518);

    auto g_0_yyyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 519);

    auto g_0_yyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 520);

    auto g_0_yyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 521);

    auto g_0_yyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 522);

    auto g_0_yyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 523);

    auto g_0_yyyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 524);

    auto g_0_yyzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 525);

    auto g_0_yyzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 526);

    auto g_0_yyzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 527);

    auto g_0_yyzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 528);

    auto g_0_yyzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 529);

    auto g_0_yyzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 530);

    auto g_0_yyzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 531);

    auto g_0_yyzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 532);

    auto g_0_yyzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 533);

    auto g_0_yyzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 534);

    auto g_0_yyzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 535);

    auto g_0_yyzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 536);

    auto g_0_yyzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 537);

    auto g_0_yyzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 538);

    auto g_0_yyzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 539);

    auto g_0_yyzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 540);

    auto g_0_yyzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 541);

    auto g_0_yyzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 542);

    auto g_0_yyzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 543);

    auto g_0_yyzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 544);

    auto g_0_yyzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 545);

    auto g_0_yzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 546);

    auto g_0_yzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 547);

    auto g_0_yzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 548);

    auto g_0_yzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 549);

    auto g_0_yzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 550);

    auto g_0_yzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 551);

    auto g_0_yzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 552);

    auto g_0_yzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 553);

    auto g_0_yzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 554);

    auto g_0_yzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 555);

    auto g_0_yzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 556);

    auto g_0_yzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 557);

    auto g_0_yzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 558);

    auto g_0_yzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 559);

    auto g_0_yzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 560);

    auto g_0_yzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 561);

    auto g_0_yzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 562);

    auto g_0_yzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 563);

    auto g_0_yzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 564);

    auto g_0_yzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 565);

    auto g_0_yzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 566);

    auto g_0_zzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 567);

    auto g_0_zzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 568);

    auto g_0_zzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 569);

    auto g_0_zzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 570);

    auto g_0_zzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 571);

    auto g_0_zzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 572);

    auto g_0_zzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 573);

    auto g_0_zzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 574);

    auto g_0_zzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 575);

    auto g_0_zzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 576);

    auto g_0_zzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 577);

    auto g_0_zzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 578);

    auto g_0_zzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 579);

    auto g_0_zzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 580);

    auto g_0_zzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 581);

    auto g_0_zzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 582);

    auto g_0_zzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 583);

    auto g_0_zzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 584);

    auto g_0_zzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 585);

    auto g_0_zzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 586);

    auto g_0_zzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 587);

    /// Set up 0-21 components of targeted buffer : SKSH

    auto g_0_xxxxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh);

    auto g_0_xxxxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 1);

    auto g_0_xxxxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 2);

    auto g_0_xxxxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 3);

    auto g_0_xxxxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 4);

    auto g_0_xxxxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 5);

    auto g_0_xxxxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 6);

    auto g_0_xxxxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 7);

    auto g_0_xxxxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 8);

    auto g_0_xxxxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 9);

    auto g_0_xxxxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 10);

    auto g_0_xxxxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 11);

    auto g_0_xxxxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 12);

    auto g_0_xxxxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 13);

    auto g_0_xxxxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 14);

    auto g_0_xxxxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 15);

    auto g_0_xxxxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 16);

    auto g_0_xxxxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 17);

    auto g_0_xxxxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 18);

    auto g_0_xxxxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 19);

    auto g_0_xxxxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 20);

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxx_0, g_0_xxxxx_0_xxxxx_1, g_0_xxxxx_0_xxxxy_0, g_0_xxxxx_0_xxxxy_1, g_0_xxxxx_0_xxxxz_0, g_0_xxxxx_0_xxxxz_1, g_0_xxxxx_0_xxxyy_0, g_0_xxxxx_0_xxxyy_1, g_0_xxxxx_0_xxxyz_0, g_0_xxxxx_0_xxxyz_1, g_0_xxxxx_0_xxxzz_0, g_0_xxxxx_0_xxxzz_1, g_0_xxxxx_0_xxyyy_0, g_0_xxxxx_0_xxyyy_1, g_0_xxxxx_0_xxyyz_0, g_0_xxxxx_0_xxyyz_1, g_0_xxxxx_0_xxyzz_0, g_0_xxxxx_0_xxyzz_1, g_0_xxxxx_0_xxzzz_0, g_0_xxxxx_0_xxzzz_1, g_0_xxxxx_0_xyyyy_0, g_0_xxxxx_0_xyyyy_1, g_0_xxxxx_0_xyyyz_0, g_0_xxxxx_0_xyyyz_1, g_0_xxxxx_0_xyyzz_0, g_0_xxxxx_0_xyyzz_1, g_0_xxxxx_0_xyzzz_0, g_0_xxxxx_0_xyzzz_1, g_0_xxxxx_0_xzzzz_0, g_0_xxxxx_0_xzzzz_1, g_0_xxxxx_0_yyyyy_0, g_0_xxxxx_0_yyyyy_1, g_0_xxxxx_0_yyyyz_0, g_0_xxxxx_0_yyyyz_1, g_0_xxxxx_0_yyyzz_0, g_0_xxxxx_0_yyyzz_1, g_0_xxxxx_0_yyzzz_0, g_0_xxxxx_0_yyzzz_1, g_0_xxxxx_0_yzzzz_0, g_0_xxxxx_0_yzzzz_1, g_0_xxxxx_0_zzzzz_0, g_0_xxxxx_0_zzzzz_1, g_0_xxxxxx_0_xxxx_1, g_0_xxxxxx_0_xxxxx_0, g_0_xxxxxx_0_xxxxx_1, g_0_xxxxxx_0_xxxxy_0, g_0_xxxxxx_0_xxxxy_1, g_0_xxxxxx_0_xxxxz_0, g_0_xxxxxx_0_xxxxz_1, g_0_xxxxxx_0_xxxy_1, g_0_xxxxxx_0_xxxyy_0, g_0_xxxxxx_0_xxxyy_1, g_0_xxxxxx_0_xxxyz_0, g_0_xxxxxx_0_xxxyz_1, g_0_xxxxxx_0_xxxz_1, g_0_xxxxxx_0_xxxzz_0, g_0_xxxxxx_0_xxxzz_1, g_0_xxxxxx_0_xxyy_1, g_0_xxxxxx_0_xxyyy_0, g_0_xxxxxx_0_xxyyy_1, g_0_xxxxxx_0_xxyyz_0, g_0_xxxxxx_0_xxyyz_1, g_0_xxxxxx_0_xxyz_1, g_0_xxxxxx_0_xxyzz_0, g_0_xxxxxx_0_xxyzz_1, g_0_xxxxxx_0_xxzz_1, g_0_xxxxxx_0_xxzzz_0, g_0_xxxxxx_0_xxzzz_1, g_0_xxxxxx_0_xyyy_1, g_0_xxxxxx_0_xyyyy_0, g_0_xxxxxx_0_xyyyy_1, g_0_xxxxxx_0_xyyyz_0, g_0_xxxxxx_0_xyyyz_1, g_0_xxxxxx_0_xyyz_1, g_0_xxxxxx_0_xyyzz_0, g_0_xxxxxx_0_xyyzz_1, g_0_xxxxxx_0_xyzz_1, g_0_xxxxxx_0_xyzzz_0, g_0_xxxxxx_0_xyzzz_1, g_0_xxxxxx_0_xzzz_1, g_0_xxxxxx_0_xzzzz_0, g_0_xxxxxx_0_xzzzz_1, g_0_xxxxxx_0_yyyy_1, g_0_xxxxxx_0_yyyyy_0, g_0_xxxxxx_0_yyyyy_1, g_0_xxxxxx_0_yyyyz_0, g_0_xxxxxx_0_yyyyz_1, g_0_xxxxxx_0_yyyz_1, g_0_xxxxxx_0_yyyzz_0, g_0_xxxxxx_0_yyyzz_1, g_0_xxxxxx_0_yyzz_1, g_0_xxxxxx_0_yyzzz_0, g_0_xxxxxx_0_yyzzz_1, g_0_xxxxxx_0_yzzz_1, g_0_xxxxxx_0_yzzzz_0, g_0_xxxxxx_0_yzzzz_1, g_0_xxxxxx_0_zzzz_1, g_0_xxxxxx_0_zzzzz_0, g_0_xxxxxx_0_zzzzz_1, g_0_xxxxxxx_0_xxxxx_0, g_0_xxxxxxx_0_xxxxy_0, g_0_xxxxxxx_0_xxxxz_0, g_0_xxxxxxx_0_xxxyy_0, g_0_xxxxxxx_0_xxxyz_0, g_0_xxxxxxx_0_xxxzz_0, g_0_xxxxxxx_0_xxyyy_0, g_0_xxxxxxx_0_xxyyz_0, g_0_xxxxxxx_0_xxyzz_0, g_0_xxxxxxx_0_xxzzz_0, g_0_xxxxxxx_0_xyyyy_0, g_0_xxxxxxx_0_xyyyz_0, g_0_xxxxxxx_0_xyyzz_0, g_0_xxxxxxx_0_xyzzz_0, g_0_xxxxxxx_0_xzzzz_0, g_0_xxxxxxx_0_yyyyy_0, g_0_xxxxxxx_0_yyyyz_0, g_0_xxxxxxx_0_yyyzz_0, g_0_xxxxxxx_0_yyzzz_0, g_0_xxxxxxx_0_yzzzz_0, g_0_xxxxxxx_0_zzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_xxxxx_0[i] = 6.0 * g_0_xxxxx_0_xxxxx_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxx_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxx_0[i] * pb_x + g_0_xxxxxx_0_xxxxx_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxy_0[i] = 6.0 * g_0_xxxxx_0_xxxxy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxy_0[i] * pb_x + g_0_xxxxxx_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxz_0[i] = 6.0 * g_0_xxxxx_0_xxxxz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxz_0[i] * pb_x + g_0_xxxxxx_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyy_0[i] = 6.0 * g_0_xxxxx_0_xxxyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyy_0[i] * pb_x + g_0_xxxxxx_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyz_0[i] = 6.0 * g_0_xxxxx_0_xxxyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyz_0[i] * pb_x + g_0_xxxxxx_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxzz_0[i] = 6.0 * g_0_xxxxx_0_xxxzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxzz_0[i] * pb_x + g_0_xxxxxx_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyy_0[i] = 6.0 * g_0_xxxxx_0_xxyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyy_0[i] * pb_x + g_0_xxxxxx_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyz_0[i] = 6.0 * g_0_xxxxx_0_xxyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyz_0[i] * pb_x + g_0_xxxxxx_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyzz_0[i] = 6.0 * g_0_xxxxx_0_xxyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzz_0[i] * pb_x + g_0_xxxxxx_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxzzz_0[i] = 6.0 * g_0_xxxxx_0_xxzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzzz_0[i] * pb_x + g_0_xxxxxx_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyy_0[i] = 6.0 * g_0_xxxxx_0_xyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyy_0[i] * pb_x + g_0_xxxxxx_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyz_0[i] = 6.0 * g_0_xxxxx_0_xyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyz_0[i] * pb_x + g_0_xxxxxx_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyzz_0[i] = 6.0 * g_0_xxxxx_0_xyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzz_0[i] * pb_x + g_0_xxxxxx_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyzzz_0[i] = 6.0 * g_0_xxxxx_0_xyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzz_0[i] * pb_x + g_0_xxxxxx_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xzzzz_0[i] = 6.0 * g_0_xxxxx_0_xzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzzz_0[i] * pb_x + g_0_xxxxxx_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyy_0[i] = 6.0 * g_0_xxxxx_0_yyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyy_0[i] * pb_x + g_0_xxxxxx_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyz_0[i] = 6.0 * g_0_xxxxx_0_yyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyz_0[i] * pb_x + g_0_xxxxxx_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyzz_0[i] = 6.0 * g_0_xxxxx_0_yyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyzz_0[i] * pb_x + g_0_xxxxxx_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyzzz_0[i] = 6.0 * g_0_xxxxx_0_yyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyzzz_0[i] * pb_x + g_0_xxxxxx_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yzzzz_0[i] = 6.0 * g_0_xxxxx_0_yzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzzzz_0[i] * pb_x + g_0_xxxxxx_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_zzzzz_0[i] = 6.0 * g_0_xxxxx_0_zzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzzzz_0[i] * pb_x + g_0_xxxxxx_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : SKSH

    auto g_0_xxxxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 21);

    auto g_0_xxxxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 22);

    auto g_0_xxxxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 23);

    auto g_0_xxxxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 24);

    auto g_0_xxxxxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 25);

    auto g_0_xxxxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 26);

    auto g_0_xxxxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 27);

    auto g_0_xxxxxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 28);

    auto g_0_xxxxxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 29);

    auto g_0_xxxxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 30);

    auto g_0_xxxxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 31);

    auto g_0_xxxxxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 32);

    auto g_0_xxxxxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 33);

    auto g_0_xxxxxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 34);

    auto g_0_xxxxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 35);

    auto g_0_xxxxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 36);

    auto g_0_xxxxxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 37);

    auto g_0_xxxxxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 38);

    auto g_0_xxxxxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 39);

    auto g_0_xxxxxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 40);

    auto g_0_xxxxxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 41);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxx_1, g_0_xxxxxx_0_xxxxx_0, g_0_xxxxxx_0_xxxxx_1, g_0_xxxxxx_0_xxxxy_0, g_0_xxxxxx_0_xxxxy_1, g_0_xxxxxx_0_xxxxz_0, g_0_xxxxxx_0_xxxxz_1, g_0_xxxxxx_0_xxxy_1, g_0_xxxxxx_0_xxxyy_0, g_0_xxxxxx_0_xxxyy_1, g_0_xxxxxx_0_xxxyz_0, g_0_xxxxxx_0_xxxyz_1, g_0_xxxxxx_0_xxxz_1, g_0_xxxxxx_0_xxxzz_0, g_0_xxxxxx_0_xxxzz_1, g_0_xxxxxx_0_xxyy_1, g_0_xxxxxx_0_xxyyy_0, g_0_xxxxxx_0_xxyyy_1, g_0_xxxxxx_0_xxyyz_0, g_0_xxxxxx_0_xxyyz_1, g_0_xxxxxx_0_xxyz_1, g_0_xxxxxx_0_xxyzz_0, g_0_xxxxxx_0_xxyzz_1, g_0_xxxxxx_0_xxzz_1, g_0_xxxxxx_0_xxzzz_0, g_0_xxxxxx_0_xxzzz_1, g_0_xxxxxx_0_xyyy_1, g_0_xxxxxx_0_xyyyy_0, g_0_xxxxxx_0_xyyyy_1, g_0_xxxxxx_0_xyyyz_0, g_0_xxxxxx_0_xyyyz_1, g_0_xxxxxx_0_xyyz_1, g_0_xxxxxx_0_xyyzz_0, g_0_xxxxxx_0_xyyzz_1, g_0_xxxxxx_0_xyzz_1, g_0_xxxxxx_0_xyzzz_0, g_0_xxxxxx_0_xyzzz_1, g_0_xxxxxx_0_xzzz_1, g_0_xxxxxx_0_xzzzz_0, g_0_xxxxxx_0_xzzzz_1, g_0_xxxxxx_0_yyyy_1, g_0_xxxxxx_0_yyyyy_0, g_0_xxxxxx_0_yyyyy_1, g_0_xxxxxx_0_yyyyz_0, g_0_xxxxxx_0_yyyyz_1, g_0_xxxxxx_0_yyyz_1, g_0_xxxxxx_0_yyyzz_0, g_0_xxxxxx_0_yyyzz_1, g_0_xxxxxx_0_yyzz_1, g_0_xxxxxx_0_yyzzz_0, g_0_xxxxxx_0_yyzzz_1, g_0_xxxxxx_0_yzzz_1, g_0_xxxxxx_0_yzzzz_0, g_0_xxxxxx_0_yzzzz_1, g_0_xxxxxx_0_zzzz_1, g_0_xxxxxx_0_zzzzz_0, g_0_xxxxxx_0_zzzzz_1, g_0_xxxxxxy_0_xxxxx_0, g_0_xxxxxxy_0_xxxxy_0, g_0_xxxxxxy_0_xxxxz_0, g_0_xxxxxxy_0_xxxyy_0, g_0_xxxxxxy_0_xxxyz_0, g_0_xxxxxxy_0_xxxzz_0, g_0_xxxxxxy_0_xxyyy_0, g_0_xxxxxxy_0_xxyyz_0, g_0_xxxxxxy_0_xxyzz_0, g_0_xxxxxxy_0_xxzzz_0, g_0_xxxxxxy_0_xyyyy_0, g_0_xxxxxxy_0_xyyyz_0, g_0_xxxxxxy_0_xyyzz_0, g_0_xxxxxxy_0_xyzzz_0, g_0_xxxxxxy_0_xzzzz_0, g_0_xxxxxxy_0_yyyyy_0, g_0_xxxxxxy_0_yyyyz_0, g_0_xxxxxxy_0_yyyzz_0, g_0_xxxxxxy_0_yyzzz_0, g_0_xxxxxxy_0_yzzzz_0, g_0_xxxxxxy_0_zzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_xxxxx_0[i] = g_0_xxxxxx_0_xxxxx_0[i] * pb_y + g_0_xxxxxx_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxy_0[i] = g_0_xxxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxy_0[i] * pb_y + g_0_xxxxxx_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxz_0[i] = g_0_xxxxxx_0_xxxxz_0[i] * pb_y + g_0_xxxxxx_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyy_0[i] = 2.0 * g_0_xxxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyy_0[i] * pb_y + g_0_xxxxxx_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyz_0[i] = g_0_xxxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyz_0[i] * pb_y + g_0_xxxxxx_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxzz_0[i] = g_0_xxxxxx_0_xxxzz_0[i] * pb_y + g_0_xxxxxx_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyy_0[i] = 3.0 * g_0_xxxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyy_0[i] * pb_y + g_0_xxxxxx_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyz_0[i] = 2.0 * g_0_xxxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyz_0[i] * pb_y + g_0_xxxxxx_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyzz_0[i] = g_0_xxxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzz_0[i] * pb_y + g_0_xxxxxx_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxzzz_0[i] = g_0_xxxxxx_0_xxzzz_0[i] * pb_y + g_0_xxxxxx_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyy_0[i] = 4.0 * g_0_xxxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyy_0[i] * pb_y + g_0_xxxxxx_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyz_0[i] = 3.0 * g_0_xxxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyz_0[i] * pb_y + g_0_xxxxxx_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzz_0[i] * pb_y + g_0_xxxxxx_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyzzz_0[i] = g_0_xxxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzz_0[i] * pb_y + g_0_xxxxxx_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xzzzz_0[i] = g_0_xxxxxx_0_xzzzz_0[i] * pb_y + g_0_xxxxxx_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyy_0[i] = 5.0 * g_0_xxxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyy_0[i] * pb_y + g_0_xxxxxx_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyz_0[i] = 4.0 * g_0_xxxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyz_0[i] * pb_y + g_0_xxxxxx_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyzz_0[i] = 3.0 * g_0_xxxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyzz_0[i] * pb_y + g_0_xxxxxx_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyzzz_0[i] = 2.0 * g_0_xxxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzzz_0[i] * pb_y + g_0_xxxxxx_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yzzzz_0[i] = g_0_xxxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzzz_0[i] * pb_y + g_0_xxxxxx_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_zzzzz_0[i] = g_0_xxxxxx_0_zzzzz_0[i] * pb_y + g_0_xxxxxx_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : SKSH

    auto g_0_xxxxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 42);

    auto g_0_xxxxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 43);

    auto g_0_xxxxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 44);

    auto g_0_xxxxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 45);

    auto g_0_xxxxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 46);

    auto g_0_xxxxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 47);

    auto g_0_xxxxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 48);

    auto g_0_xxxxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 49);

    auto g_0_xxxxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 50);

    auto g_0_xxxxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 51);

    auto g_0_xxxxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 52);

    auto g_0_xxxxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 53);

    auto g_0_xxxxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 54);

    auto g_0_xxxxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 55);

    auto g_0_xxxxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 56);

    auto g_0_xxxxxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 57);

    auto g_0_xxxxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 58);

    auto g_0_xxxxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 59);

    auto g_0_xxxxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 60);

    auto g_0_xxxxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 61);

    auto g_0_xxxxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 62);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxx_1, g_0_xxxxxx_0_xxxxx_0, g_0_xxxxxx_0_xxxxx_1, g_0_xxxxxx_0_xxxxy_0, g_0_xxxxxx_0_xxxxy_1, g_0_xxxxxx_0_xxxxz_0, g_0_xxxxxx_0_xxxxz_1, g_0_xxxxxx_0_xxxy_1, g_0_xxxxxx_0_xxxyy_0, g_0_xxxxxx_0_xxxyy_1, g_0_xxxxxx_0_xxxyz_0, g_0_xxxxxx_0_xxxyz_1, g_0_xxxxxx_0_xxxz_1, g_0_xxxxxx_0_xxxzz_0, g_0_xxxxxx_0_xxxzz_1, g_0_xxxxxx_0_xxyy_1, g_0_xxxxxx_0_xxyyy_0, g_0_xxxxxx_0_xxyyy_1, g_0_xxxxxx_0_xxyyz_0, g_0_xxxxxx_0_xxyyz_1, g_0_xxxxxx_0_xxyz_1, g_0_xxxxxx_0_xxyzz_0, g_0_xxxxxx_0_xxyzz_1, g_0_xxxxxx_0_xxzz_1, g_0_xxxxxx_0_xxzzz_0, g_0_xxxxxx_0_xxzzz_1, g_0_xxxxxx_0_xyyy_1, g_0_xxxxxx_0_xyyyy_0, g_0_xxxxxx_0_xyyyy_1, g_0_xxxxxx_0_xyyyz_0, g_0_xxxxxx_0_xyyyz_1, g_0_xxxxxx_0_xyyz_1, g_0_xxxxxx_0_xyyzz_0, g_0_xxxxxx_0_xyyzz_1, g_0_xxxxxx_0_xyzz_1, g_0_xxxxxx_0_xyzzz_0, g_0_xxxxxx_0_xyzzz_1, g_0_xxxxxx_0_xzzz_1, g_0_xxxxxx_0_xzzzz_0, g_0_xxxxxx_0_xzzzz_1, g_0_xxxxxx_0_yyyy_1, g_0_xxxxxx_0_yyyyy_0, g_0_xxxxxx_0_yyyyy_1, g_0_xxxxxx_0_yyyyz_0, g_0_xxxxxx_0_yyyyz_1, g_0_xxxxxx_0_yyyz_1, g_0_xxxxxx_0_yyyzz_0, g_0_xxxxxx_0_yyyzz_1, g_0_xxxxxx_0_yyzz_1, g_0_xxxxxx_0_yyzzz_0, g_0_xxxxxx_0_yyzzz_1, g_0_xxxxxx_0_yzzz_1, g_0_xxxxxx_0_yzzzz_0, g_0_xxxxxx_0_yzzzz_1, g_0_xxxxxx_0_zzzz_1, g_0_xxxxxx_0_zzzzz_0, g_0_xxxxxx_0_zzzzz_1, g_0_xxxxxxz_0_xxxxx_0, g_0_xxxxxxz_0_xxxxy_0, g_0_xxxxxxz_0_xxxxz_0, g_0_xxxxxxz_0_xxxyy_0, g_0_xxxxxxz_0_xxxyz_0, g_0_xxxxxxz_0_xxxzz_0, g_0_xxxxxxz_0_xxyyy_0, g_0_xxxxxxz_0_xxyyz_0, g_0_xxxxxxz_0_xxyzz_0, g_0_xxxxxxz_0_xxzzz_0, g_0_xxxxxxz_0_xyyyy_0, g_0_xxxxxxz_0_xyyyz_0, g_0_xxxxxxz_0_xyyzz_0, g_0_xxxxxxz_0_xyzzz_0, g_0_xxxxxxz_0_xzzzz_0, g_0_xxxxxxz_0_yyyyy_0, g_0_xxxxxxz_0_yyyyz_0, g_0_xxxxxxz_0_yyyzz_0, g_0_xxxxxxz_0_yyzzz_0, g_0_xxxxxxz_0_yzzzz_0, g_0_xxxxxxz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_xxxxx_0[i] = g_0_xxxxxx_0_xxxxx_0[i] * pb_z + g_0_xxxxxx_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxy_0[i] = g_0_xxxxxx_0_xxxxy_0[i] * pb_z + g_0_xxxxxx_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxz_0[i] = g_0_xxxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxz_0[i] * pb_z + g_0_xxxxxx_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyy_0[i] = g_0_xxxxxx_0_xxxyy_0[i] * pb_z + g_0_xxxxxx_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyz_0[i] = g_0_xxxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyz_0[i] * pb_z + g_0_xxxxxx_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxzz_0[i] * pb_z + g_0_xxxxxx_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyy_0[i] = g_0_xxxxxx_0_xxyyy_0[i] * pb_z + g_0_xxxxxx_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyz_0[i] = g_0_xxxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyz_0[i] * pb_z + g_0_xxxxxx_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzz_0[i] * pb_z + g_0_xxxxxx_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxzzz_0[i] = 3.0 * g_0_xxxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzzz_0[i] * pb_z + g_0_xxxxxx_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyy_0[i] = g_0_xxxxxx_0_xyyyy_0[i] * pb_z + g_0_xxxxxx_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyz_0[i] = g_0_xxxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyz_0[i] * pb_z + g_0_xxxxxx_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzz_0[i] * pb_z + g_0_xxxxxx_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzz_0[i] * pb_z + g_0_xxxxxx_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xzzzz_0[i] = 4.0 * g_0_xxxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzzz_0[i] * pb_z + g_0_xxxxxx_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyy_0[i] = g_0_xxxxxx_0_yyyyy_0[i] * pb_z + g_0_xxxxxx_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyz_0[i] = g_0_xxxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyz_0[i] * pb_z + g_0_xxxxxx_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyzz_0[i] = 2.0 * g_0_xxxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyzz_0[i] * pb_z + g_0_xxxxxx_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyzzz_0[i] = 3.0 * g_0_xxxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzzz_0[i] * pb_z + g_0_xxxxxx_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yzzzz_0[i] = 4.0 * g_0_xxxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzzz_0[i] * pb_z + g_0_xxxxxx_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_zzzzz_0[i] = 5.0 * g_0_xxxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_zzzzz_0[i] * pb_z + g_0_xxxxxx_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 63-84 components of targeted buffer : SKSH

    auto g_0_xxxxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 63);

    auto g_0_xxxxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 64);

    auto g_0_xxxxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 65);

    auto g_0_xxxxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 66);

    auto g_0_xxxxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 67);

    auto g_0_xxxxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 68);

    auto g_0_xxxxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 69);

    auto g_0_xxxxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 70);

    auto g_0_xxxxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 71);

    auto g_0_xxxxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 72);

    auto g_0_xxxxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 73);

    auto g_0_xxxxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 74);

    auto g_0_xxxxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 75);

    auto g_0_xxxxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 76);

    auto g_0_xxxxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 77);

    auto g_0_xxxxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 78);

    auto g_0_xxxxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 79);

    auto g_0_xxxxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 80);

    auto g_0_xxxxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 81);

    auto g_0_xxxxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 82);

    auto g_0_xxxxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 83);

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxx_0, g_0_xxxxx_0_xxxxx_1, g_0_xxxxx_0_xxxxz_0, g_0_xxxxx_0_xxxxz_1, g_0_xxxxx_0_xxxzz_0, g_0_xxxxx_0_xxxzz_1, g_0_xxxxx_0_xxzzz_0, g_0_xxxxx_0_xxzzz_1, g_0_xxxxx_0_xzzzz_0, g_0_xxxxx_0_xzzzz_1, g_0_xxxxxy_0_xxxxx_0, g_0_xxxxxy_0_xxxxx_1, g_0_xxxxxy_0_xxxxz_0, g_0_xxxxxy_0_xxxxz_1, g_0_xxxxxy_0_xxxzz_0, g_0_xxxxxy_0_xxxzz_1, g_0_xxxxxy_0_xxzzz_0, g_0_xxxxxy_0_xxzzz_1, g_0_xxxxxy_0_xzzzz_0, g_0_xxxxxy_0_xzzzz_1, g_0_xxxxxyy_0_xxxxx_0, g_0_xxxxxyy_0_xxxxy_0, g_0_xxxxxyy_0_xxxxz_0, g_0_xxxxxyy_0_xxxyy_0, g_0_xxxxxyy_0_xxxyz_0, g_0_xxxxxyy_0_xxxzz_0, g_0_xxxxxyy_0_xxyyy_0, g_0_xxxxxyy_0_xxyyz_0, g_0_xxxxxyy_0_xxyzz_0, g_0_xxxxxyy_0_xxzzz_0, g_0_xxxxxyy_0_xyyyy_0, g_0_xxxxxyy_0_xyyyz_0, g_0_xxxxxyy_0_xyyzz_0, g_0_xxxxxyy_0_xyzzz_0, g_0_xxxxxyy_0_xzzzz_0, g_0_xxxxxyy_0_yyyyy_0, g_0_xxxxxyy_0_yyyyz_0, g_0_xxxxxyy_0_yyyzz_0, g_0_xxxxxyy_0_yyzzz_0, g_0_xxxxxyy_0_yzzzz_0, g_0_xxxxxyy_0_zzzzz_0, g_0_xxxxyy_0_xxxxy_0, g_0_xxxxyy_0_xxxxy_1, g_0_xxxxyy_0_xxxy_1, g_0_xxxxyy_0_xxxyy_0, g_0_xxxxyy_0_xxxyy_1, g_0_xxxxyy_0_xxxyz_0, g_0_xxxxyy_0_xxxyz_1, g_0_xxxxyy_0_xxyy_1, g_0_xxxxyy_0_xxyyy_0, g_0_xxxxyy_0_xxyyy_1, g_0_xxxxyy_0_xxyyz_0, g_0_xxxxyy_0_xxyyz_1, g_0_xxxxyy_0_xxyz_1, g_0_xxxxyy_0_xxyzz_0, g_0_xxxxyy_0_xxyzz_1, g_0_xxxxyy_0_xyyy_1, g_0_xxxxyy_0_xyyyy_0, g_0_xxxxyy_0_xyyyy_1, g_0_xxxxyy_0_xyyyz_0, g_0_xxxxyy_0_xyyyz_1, g_0_xxxxyy_0_xyyz_1, g_0_xxxxyy_0_xyyzz_0, g_0_xxxxyy_0_xyyzz_1, g_0_xxxxyy_0_xyzz_1, g_0_xxxxyy_0_xyzzz_0, g_0_xxxxyy_0_xyzzz_1, g_0_xxxxyy_0_yyyy_1, g_0_xxxxyy_0_yyyyy_0, g_0_xxxxyy_0_yyyyy_1, g_0_xxxxyy_0_yyyyz_0, g_0_xxxxyy_0_yyyyz_1, g_0_xxxxyy_0_yyyz_1, g_0_xxxxyy_0_yyyzz_0, g_0_xxxxyy_0_yyyzz_1, g_0_xxxxyy_0_yyzz_1, g_0_xxxxyy_0_yyzzz_0, g_0_xxxxyy_0_yyzzz_1, g_0_xxxxyy_0_yzzz_1, g_0_xxxxyy_0_yzzzz_0, g_0_xxxxyy_0_yzzzz_1, g_0_xxxxyy_0_zzzzz_0, g_0_xxxxyy_0_zzzzz_1, g_0_xxxyy_0_xxxxy_0, g_0_xxxyy_0_xxxxy_1, g_0_xxxyy_0_xxxyy_0, g_0_xxxyy_0_xxxyy_1, g_0_xxxyy_0_xxxyz_0, g_0_xxxyy_0_xxxyz_1, g_0_xxxyy_0_xxyyy_0, g_0_xxxyy_0_xxyyy_1, g_0_xxxyy_0_xxyyz_0, g_0_xxxyy_0_xxyyz_1, g_0_xxxyy_0_xxyzz_0, g_0_xxxyy_0_xxyzz_1, g_0_xxxyy_0_xyyyy_0, g_0_xxxyy_0_xyyyy_1, g_0_xxxyy_0_xyyyz_0, g_0_xxxyy_0_xyyyz_1, g_0_xxxyy_0_xyyzz_0, g_0_xxxyy_0_xyyzz_1, g_0_xxxyy_0_xyzzz_0, g_0_xxxyy_0_xyzzz_1, g_0_xxxyy_0_yyyyy_0, g_0_xxxyy_0_yyyyy_1, g_0_xxxyy_0_yyyyz_0, g_0_xxxyy_0_yyyyz_1, g_0_xxxyy_0_yyyzz_0, g_0_xxxyy_0_yyyzz_1, g_0_xxxyy_0_yyzzz_0, g_0_xxxyy_0_yyzzz_1, g_0_xxxyy_0_yzzzz_0, g_0_xxxyy_0_yzzzz_1, g_0_xxxyy_0_zzzzz_0, g_0_xxxyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_xxxxx_0[i] = g_0_xxxxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxx_0[i] * pb_y + g_0_xxxxxy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxxy_0[i] = 4.0 * g_0_xxxyy_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxy_0[i] * pb_x + g_0_xxxxyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxz_0[i] = g_0_xxxxx_0_xxxxz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxz_0[i] * pb_y + g_0_xxxxxy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxyy_0[i] = 4.0 * g_0_xxxyy_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyy_0[i] * pb_x + g_0_xxxxyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxyz_0[i] = 4.0 * g_0_xxxyy_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyz_0[i] * pb_x + g_0_xxxxyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxzz_0[i] = g_0_xxxxx_0_xxxzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxzz_0[i] * pb_y + g_0_xxxxxy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxyyy_0[i] = 4.0 * g_0_xxxyy_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyy_0[i] * pb_x + g_0_xxxxyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyyz_0[i] = 4.0 * g_0_xxxyy_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyz_0[i] * pb_x + g_0_xxxxyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyzz_0[i] = 4.0 * g_0_xxxyy_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyzz_0[i] * pb_x + g_0_xxxxyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxzzz_0[i] = g_0_xxxxx_0_xxzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxzzz_0[i] * pb_y + g_0_xxxxxy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xyyyy_0[i] = 4.0 * g_0_xxxyy_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyy_0[i] * pb_x + g_0_xxxxyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyyz_0[i] = 4.0 * g_0_xxxyy_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyz_0[i] * pb_x + g_0_xxxxyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyzz_0[i] = 4.0 * g_0_xxxyy_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyzz_0[i] * pb_x + g_0_xxxxyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyzzz_0[i] = 4.0 * g_0_xxxyy_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzzz_0[i] * pb_x + g_0_xxxxyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xzzzz_0[i] = g_0_xxxxx_0_xzzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xzzzz_0[i] * pb_y + g_0_xxxxxy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_yyyyy_0[i] = 4.0 * g_0_xxxyy_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyy_0[i] * pb_x + g_0_xxxxyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyyz_0[i] = 4.0 * g_0_xxxyy_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyz_0[i] * pb_x + g_0_xxxxyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyzz_0[i] = 4.0 * g_0_xxxyy_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyzz_0[i] * pb_x + g_0_xxxxyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyzzz_0[i] = 4.0 * g_0_xxxyy_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyzzz_0[i] * pb_x + g_0_xxxxyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yzzzz_0[i] = 4.0 * g_0_xxxyy_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzzzz_0[i] * pb_x + g_0_xxxxyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_zzzzz_0[i] = 4.0 * g_0_xxxyy_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_zzzzz_0[i] * pb_x + g_0_xxxxyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-105 components of targeted buffer : SKSH

    auto g_0_xxxxxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 84);

    auto g_0_xxxxxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 85);

    auto g_0_xxxxxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 86);

    auto g_0_xxxxxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 87);

    auto g_0_xxxxxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 88);

    auto g_0_xxxxxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 89);

    auto g_0_xxxxxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 90);

    auto g_0_xxxxxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 91);

    auto g_0_xxxxxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 92);

    auto g_0_xxxxxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 93);

    auto g_0_xxxxxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 94);

    auto g_0_xxxxxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 95);

    auto g_0_xxxxxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 96);

    auto g_0_xxxxxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 97);

    auto g_0_xxxxxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 98);

    auto g_0_xxxxxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 99);

    auto g_0_xxxxxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 100);

    auto g_0_xxxxxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 101);

    auto g_0_xxxxxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 102);

    auto g_0_xxxxxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 103);

    auto g_0_xxxxxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 104);

    #pragma omp simd aligned(g_0_xxxxxy_0_xxxxy_0, g_0_xxxxxy_0_xxxxy_1, g_0_xxxxxy_0_xxxyy_0, g_0_xxxxxy_0_xxxyy_1, g_0_xxxxxy_0_xxyyy_0, g_0_xxxxxy_0_xxyyy_1, g_0_xxxxxy_0_xyyyy_0, g_0_xxxxxy_0_xyyyy_1, g_0_xxxxxy_0_yyyyy_0, g_0_xxxxxy_0_yyyyy_1, g_0_xxxxxyz_0_xxxxx_0, g_0_xxxxxyz_0_xxxxy_0, g_0_xxxxxyz_0_xxxxz_0, g_0_xxxxxyz_0_xxxyy_0, g_0_xxxxxyz_0_xxxyz_0, g_0_xxxxxyz_0_xxxzz_0, g_0_xxxxxyz_0_xxyyy_0, g_0_xxxxxyz_0_xxyyz_0, g_0_xxxxxyz_0_xxyzz_0, g_0_xxxxxyz_0_xxzzz_0, g_0_xxxxxyz_0_xyyyy_0, g_0_xxxxxyz_0_xyyyz_0, g_0_xxxxxyz_0_xyyzz_0, g_0_xxxxxyz_0_xyzzz_0, g_0_xxxxxyz_0_xzzzz_0, g_0_xxxxxyz_0_yyyyy_0, g_0_xxxxxyz_0_yyyyz_0, g_0_xxxxxyz_0_yyyzz_0, g_0_xxxxxyz_0_yyzzz_0, g_0_xxxxxyz_0_yzzzz_0, g_0_xxxxxyz_0_zzzzz_0, g_0_xxxxxz_0_xxxxx_0, g_0_xxxxxz_0_xxxxx_1, g_0_xxxxxz_0_xxxxz_0, g_0_xxxxxz_0_xxxxz_1, g_0_xxxxxz_0_xxxyz_0, g_0_xxxxxz_0_xxxyz_1, g_0_xxxxxz_0_xxxz_1, g_0_xxxxxz_0_xxxzz_0, g_0_xxxxxz_0_xxxzz_1, g_0_xxxxxz_0_xxyyz_0, g_0_xxxxxz_0_xxyyz_1, g_0_xxxxxz_0_xxyz_1, g_0_xxxxxz_0_xxyzz_0, g_0_xxxxxz_0_xxyzz_1, g_0_xxxxxz_0_xxzz_1, g_0_xxxxxz_0_xxzzz_0, g_0_xxxxxz_0_xxzzz_1, g_0_xxxxxz_0_xyyyz_0, g_0_xxxxxz_0_xyyyz_1, g_0_xxxxxz_0_xyyz_1, g_0_xxxxxz_0_xyyzz_0, g_0_xxxxxz_0_xyyzz_1, g_0_xxxxxz_0_xyzz_1, g_0_xxxxxz_0_xyzzz_0, g_0_xxxxxz_0_xyzzz_1, g_0_xxxxxz_0_xzzz_1, g_0_xxxxxz_0_xzzzz_0, g_0_xxxxxz_0_xzzzz_1, g_0_xxxxxz_0_yyyyz_0, g_0_xxxxxz_0_yyyyz_1, g_0_xxxxxz_0_yyyz_1, g_0_xxxxxz_0_yyyzz_0, g_0_xxxxxz_0_yyyzz_1, g_0_xxxxxz_0_yyzz_1, g_0_xxxxxz_0_yyzzz_0, g_0_xxxxxz_0_yyzzz_1, g_0_xxxxxz_0_yzzz_1, g_0_xxxxxz_0_yzzzz_0, g_0_xxxxxz_0_yzzzz_1, g_0_xxxxxz_0_zzzz_1, g_0_xxxxxz_0_zzzzz_0, g_0_xxxxxz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyz_0_xxxxx_0[i] = g_0_xxxxxz_0_xxxxx_0[i] * pb_y + g_0_xxxxxz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxy_0[i] = g_0_xxxxxy_0_xxxxy_0[i] * pb_z + g_0_xxxxxy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxxz_0[i] = g_0_xxxxxz_0_xxxxz_0[i] * pb_y + g_0_xxxxxz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxyy_0[i] = g_0_xxxxxy_0_xxxyy_0[i] * pb_z + g_0_xxxxxy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxyz_0[i] = g_0_xxxxxz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxyz_0[i] * pb_y + g_0_xxxxxz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxzz_0[i] = g_0_xxxxxz_0_xxxzz_0[i] * pb_y + g_0_xxxxxz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyyy_0[i] = g_0_xxxxxy_0_xxyyy_0[i] * pb_z + g_0_xxxxxy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxyyz_0[i] = 2.0 * g_0_xxxxxz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyyz_0[i] * pb_y + g_0_xxxxxz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyzz_0[i] = g_0_xxxxxz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyzz_0[i] * pb_y + g_0_xxxxxz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxzzz_0[i] = g_0_xxxxxz_0_xxzzz_0[i] * pb_y + g_0_xxxxxz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyyy_0[i] = g_0_xxxxxy_0_xyyyy_0[i] * pb_z + g_0_xxxxxy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xyyyz_0[i] = 3.0 * g_0_xxxxxz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyyz_0[i] * pb_y + g_0_xxxxxz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyzz_0[i] = 2.0 * g_0_xxxxxz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyzz_0[i] * pb_y + g_0_xxxxxz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyzzz_0[i] = g_0_xxxxxz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyzzz_0[i] * pb_y + g_0_xxxxxz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xzzzz_0[i] = g_0_xxxxxz_0_xzzzz_0[i] * pb_y + g_0_xxxxxz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyyy_0[i] = g_0_xxxxxy_0_yyyyy_0[i] * pb_z + g_0_xxxxxy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_yyyyz_0[i] = 4.0 * g_0_xxxxxz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyyz_0[i] * pb_y + g_0_xxxxxz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyzz_0[i] = 3.0 * g_0_xxxxxz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyzz_0[i] * pb_y + g_0_xxxxxz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyzzz_0[i] = 2.0 * g_0_xxxxxz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyzzz_0[i] * pb_y + g_0_xxxxxz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yzzzz_0[i] = g_0_xxxxxz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yzzzz_0[i] * pb_y + g_0_xxxxxz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_zzzzz_0[i] = g_0_xxxxxz_0_zzzzz_0[i] * pb_y + g_0_xxxxxz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 105-126 components of targeted buffer : SKSH

    auto g_0_xxxxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 105);

    auto g_0_xxxxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 106);

    auto g_0_xxxxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 107);

    auto g_0_xxxxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 108);

    auto g_0_xxxxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 109);

    auto g_0_xxxxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 110);

    auto g_0_xxxxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 111);

    auto g_0_xxxxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 112);

    auto g_0_xxxxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 113);

    auto g_0_xxxxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 114);

    auto g_0_xxxxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 115);

    auto g_0_xxxxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 116);

    auto g_0_xxxxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 117);

    auto g_0_xxxxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 118);

    auto g_0_xxxxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 119);

    auto g_0_xxxxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 120);

    auto g_0_xxxxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 121);

    auto g_0_xxxxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 122);

    auto g_0_xxxxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 123);

    auto g_0_xxxxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 124);

    auto g_0_xxxxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 125);

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxx_0, g_0_xxxxx_0_xxxxx_1, g_0_xxxxx_0_xxxxy_0, g_0_xxxxx_0_xxxxy_1, g_0_xxxxx_0_xxxyy_0, g_0_xxxxx_0_xxxyy_1, g_0_xxxxx_0_xxyyy_0, g_0_xxxxx_0_xxyyy_1, g_0_xxxxx_0_xyyyy_0, g_0_xxxxx_0_xyyyy_1, g_0_xxxxxz_0_xxxxx_0, g_0_xxxxxz_0_xxxxx_1, g_0_xxxxxz_0_xxxxy_0, g_0_xxxxxz_0_xxxxy_1, g_0_xxxxxz_0_xxxyy_0, g_0_xxxxxz_0_xxxyy_1, g_0_xxxxxz_0_xxyyy_0, g_0_xxxxxz_0_xxyyy_1, g_0_xxxxxz_0_xyyyy_0, g_0_xxxxxz_0_xyyyy_1, g_0_xxxxxzz_0_xxxxx_0, g_0_xxxxxzz_0_xxxxy_0, g_0_xxxxxzz_0_xxxxz_0, g_0_xxxxxzz_0_xxxyy_0, g_0_xxxxxzz_0_xxxyz_0, g_0_xxxxxzz_0_xxxzz_0, g_0_xxxxxzz_0_xxyyy_0, g_0_xxxxxzz_0_xxyyz_0, g_0_xxxxxzz_0_xxyzz_0, g_0_xxxxxzz_0_xxzzz_0, g_0_xxxxxzz_0_xyyyy_0, g_0_xxxxxzz_0_xyyyz_0, g_0_xxxxxzz_0_xyyzz_0, g_0_xxxxxzz_0_xyzzz_0, g_0_xxxxxzz_0_xzzzz_0, g_0_xxxxxzz_0_yyyyy_0, g_0_xxxxxzz_0_yyyyz_0, g_0_xxxxxzz_0_yyyzz_0, g_0_xxxxxzz_0_yyzzz_0, g_0_xxxxxzz_0_yzzzz_0, g_0_xxxxxzz_0_zzzzz_0, g_0_xxxxzz_0_xxxxz_0, g_0_xxxxzz_0_xxxxz_1, g_0_xxxxzz_0_xxxyz_0, g_0_xxxxzz_0_xxxyz_1, g_0_xxxxzz_0_xxxz_1, g_0_xxxxzz_0_xxxzz_0, g_0_xxxxzz_0_xxxzz_1, g_0_xxxxzz_0_xxyyz_0, g_0_xxxxzz_0_xxyyz_1, g_0_xxxxzz_0_xxyz_1, g_0_xxxxzz_0_xxyzz_0, g_0_xxxxzz_0_xxyzz_1, g_0_xxxxzz_0_xxzz_1, g_0_xxxxzz_0_xxzzz_0, g_0_xxxxzz_0_xxzzz_1, g_0_xxxxzz_0_xyyyz_0, g_0_xxxxzz_0_xyyyz_1, g_0_xxxxzz_0_xyyz_1, g_0_xxxxzz_0_xyyzz_0, g_0_xxxxzz_0_xyyzz_1, g_0_xxxxzz_0_xyzz_1, g_0_xxxxzz_0_xyzzz_0, g_0_xxxxzz_0_xyzzz_1, g_0_xxxxzz_0_xzzz_1, g_0_xxxxzz_0_xzzzz_0, g_0_xxxxzz_0_xzzzz_1, g_0_xxxxzz_0_yyyyy_0, g_0_xxxxzz_0_yyyyy_1, g_0_xxxxzz_0_yyyyz_0, g_0_xxxxzz_0_yyyyz_1, g_0_xxxxzz_0_yyyz_1, g_0_xxxxzz_0_yyyzz_0, g_0_xxxxzz_0_yyyzz_1, g_0_xxxxzz_0_yyzz_1, g_0_xxxxzz_0_yyzzz_0, g_0_xxxxzz_0_yyzzz_1, g_0_xxxxzz_0_yzzz_1, g_0_xxxxzz_0_yzzzz_0, g_0_xxxxzz_0_yzzzz_1, g_0_xxxxzz_0_zzzz_1, g_0_xxxxzz_0_zzzzz_0, g_0_xxxxzz_0_zzzzz_1, g_0_xxxzz_0_xxxxz_0, g_0_xxxzz_0_xxxxz_1, g_0_xxxzz_0_xxxyz_0, g_0_xxxzz_0_xxxyz_1, g_0_xxxzz_0_xxxzz_0, g_0_xxxzz_0_xxxzz_1, g_0_xxxzz_0_xxyyz_0, g_0_xxxzz_0_xxyyz_1, g_0_xxxzz_0_xxyzz_0, g_0_xxxzz_0_xxyzz_1, g_0_xxxzz_0_xxzzz_0, g_0_xxxzz_0_xxzzz_1, g_0_xxxzz_0_xyyyz_0, g_0_xxxzz_0_xyyyz_1, g_0_xxxzz_0_xyyzz_0, g_0_xxxzz_0_xyyzz_1, g_0_xxxzz_0_xyzzz_0, g_0_xxxzz_0_xyzzz_1, g_0_xxxzz_0_xzzzz_0, g_0_xxxzz_0_xzzzz_1, g_0_xxxzz_0_yyyyy_0, g_0_xxxzz_0_yyyyy_1, g_0_xxxzz_0_yyyyz_0, g_0_xxxzz_0_yyyyz_1, g_0_xxxzz_0_yyyzz_0, g_0_xxxzz_0_yyyzz_1, g_0_xxxzz_0_yyzzz_0, g_0_xxxzz_0_yyzzz_1, g_0_xxxzz_0_yzzzz_0, g_0_xxxzz_0_yzzzz_1, g_0_xxxzz_0_zzzzz_0, g_0_xxxzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_xxxxx_0[i] = g_0_xxxxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxx_0[i] * pb_z + g_0_xxxxxz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxy_0[i] = g_0_xxxxx_0_xxxxy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxy_0[i] * pb_z + g_0_xxxxxz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxz_0[i] = 4.0 * g_0_xxxzz_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxz_0[i] * pb_x + g_0_xxxxzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxyy_0[i] = g_0_xxxxx_0_xxxyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxyy_0[i] * pb_z + g_0_xxxxxz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxyz_0[i] = 4.0 * g_0_xxxzz_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyz_0[i] * pb_x + g_0_xxxxzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxzz_0[i] = 4.0 * g_0_xxxzz_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxzz_0[i] * pb_x + g_0_xxxxzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyyy_0[i] = g_0_xxxxx_0_xxyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxyyy_0[i] * pb_z + g_0_xxxxxz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxyyz_0[i] = 4.0 * g_0_xxxzz_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyz_0[i] * pb_x + g_0_xxxxzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyzz_0[i] = 4.0 * g_0_xxxzz_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyzz_0[i] * pb_x + g_0_xxxxzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxzzz_0[i] = 4.0 * g_0_xxxzz_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxzzz_0[i] * pb_x + g_0_xxxxzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyyy_0[i] = g_0_xxxxx_0_xyyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xyyyy_0[i] * pb_z + g_0_xxxxxz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xyyyz_0[i] = 4.0 * g_0_xxxzz_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyz_0[i] * pb_x + g_0_xxxxzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyzz_0[i] = 4.0 * g_0_xxxzz_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyzz_0[i] * pb_x + g_0_xxxxzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyzzz_0[i] = 4.0 * g_0_xxxzz_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzzz_0[i] * pb_x + g_0_xxxxzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xzzzz_0[i] = 4.0 * g_0_xxxzz_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xzzzz_0[i] * pb_x + g_0_xxxxzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyy_0[i] = 4.0 * g_0_xxxzz_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyy_0[i] * pb_x + g_0_xxxxzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyz_0[i] = 4.0 * g_0_xxxzz_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyz_0[i] * pb_x + g_0_xxxxzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyzz_0[i] = 4.0 * g_0_xxxzz_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyzz_0[i] * pb_x + g_0_xxxxzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyzzz_0[i] = 4.0 * g_0_xxxzz_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyzzz_0[i] * pb_x + g_0_xxxxzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yzzzz_0[i] = 4.0 * g_0_xxxzz_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzzzz_0[i] * pb_x + g_0_xxxxzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_zzzzz_0[i] = 4.0 * g_0_xxxzz_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzzzz_0[i] * pb_x + g_0_xxxxzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 126-147 components of targeted buffer : SKSH

    auto g_0_xxxxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 126);

    auto g_0_xxxxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 127);

    auto g_0_xxxxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 128);

    auto g_0_xxxxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 129);

    auto g_0_xxxxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 130);

    auto g_0_xxxxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 131);

    auto g_0_xxxxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 132);

    auto g_0_xxxxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 133);

    auto g_0_xxxxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 134);

    auto g_0_xxxxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 135);

    auto g_0_xxxxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 136);

    auto g_0_xxxxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 137);

    auto g_0_xxxxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 138);

    auto g_0_xxxxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 139);

    auto g_0_xxxxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 140);

    auto g_0_xxxxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 141);

    auto g_0_xxxxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 142);

    auto g_0_xxxxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 143);

    auto g_0_xxxxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 144);

    auto g_0_xxxxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 145);

    auto g_0_xxxxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 146);

    #pragma omp simd aligned(g_0_xxxxy_0_xxxxx_0, g_0_xxxxy_0_xxxxx_1, g_0_xxxxy_0_xxxxz_0, g_0_xxxxy_0_xxxxz_1, g_0_xxxxy_0_xxxzz_0, g_0_xxxxy_0_xxxzz_1, g_0_xxxxy_0_xxzzz_0, g_0_xxxxy_0_xxzzz_1, g_0_xxxxy_0_xzzzz_0, g_0_xxxxy_0_xzzzz_1, g_0_xxxxyy_0_xxxxx_0, g_0_xxxxyy_0_xxxxx_1, g_0_xxxxyy_0_xxxxz_0, g_0_xxxxyy_0_xxxxz_1, g_0_xxxxyy_0_xxxzz_0, g_0_xxxxyy_0_xxxzz_1, g_0_xxxxyy_0_xxzzz_0, g_0_xxxxyy_0_xxzzz_1, g_0_xxxxyy_0_xzzzz_0, g_0_xxxxyy_0_xzzzz_1, g_0_xxxxyyy_0_xxxxx_0, g_0_xxxxyyy_0_xxxxy_0, g_0_xxxxyyy_0_xxxxz_0, g_0_xxxxyyy_0_xxxyy_0, g_0_xxxxyyy_0_xxxyz_0, g_0_xxxxyyy_0_xxxzz_0, g_0_xxxxyyy_0_xxyyy_0, g_0_xxxxyyy_0_xxyyz_0, g_0_xxxxyyy_0_xxyzz_0, g_0_xxxxyyy_0_xxzzz_0, g_0_xxxxyyy_0_xyyyy_0, g_0_xxxxyyy_0_xyyyz_0, g_0_xxxxyyy_0_xyyzz_0, g_0_xxxxyyy_0_xyzzz_0, g_0_xxxxyyy_0_xzzzz_0, g_0_xxxxyyy_0_yyyyy_0, g_0_xxxxyyy_0_yyyyz_0, g_0_xxxxyyy_0_yyyzz_0, g_0_xxxxyyy_0_yyzzz_0, g_0_xxxxyyy_0_yzzzz_0, g_0_xxxxyyy_0_zzzzz_0, g_0_xxxyyy_0_xxxxy_0, g_0_xxxyyy_0_xxxxy_1, g_0_xxxyyy_0_xxxy_1, g_0_xxxyyy_0_xxxyy_0, g_0_xxxyyy_0_xxxyy_1, g_0_xxxyyy_0_xxxyz_0, g_0_xxxyyy_0_xxxyz_1, g_0_xxxyyy_0_xxyy_1, g_0_xxxyyy_0_xxyyy_0, g_0_xxxyyy_0_xxyyy_1, g_0_xxxyyy_0_xxyyz_0, g_0_xxxyyy_0_xxyyz_1, g_0_xxxyyy_0_xxyz_1, g_0_xxxyyy_0_xxyzz_0, g_0_xxxyyy_0_xxyzz_1, g_0_xxxyyy_0_xyyy_1, g_0_xxxyyy_0_xyyyy_0, g_0_xxxyyy_0_xyyyy_1, g_0_xxxyyy_0_xyyyz_0, g_0_xxxyyy_0_xyyyz_1, g_0_xxxyyy_0_xyyz_1, g_0_xxxyyy_0_xyyzz_0, g_0_xxxyyy_0_xyyzz_1, g_0_xxxyyy_0_xyzz_1, g_0_xxxyyy_0_xyzzz_0, g_0_xxxyyy_0_xyzzz_1, g_0_xxxyyy_0_yyyy_1, g_0_xxxyyy_0_yyyyy_0, g_0_xxxyyy_0_yyyyy_1, g_0_xxxyyy_0_yyyyz_0, g_0_xxxyyy_0_yyyyz_1, g_0_xxxyyy_0_yyyz_1, g_0_xxxyyy_0_yyyzz_0, g_0_xxxyyy_0_yyyzz_1, g_0_xxxyyy_0_yyzz_1, g_0_xxxyyy_0_yyzzz_0, g_0_xxxyyy_0_yyzzz_1, g_0_xxxyyy_0_yzzz_1, g_0_xxxyyy_0_yzzzz_0, g_0_xxxyyy_0_yzzzz_1, g_0_xxxyyy_0_zzzzz_0, g_0_xxxyyy_0_zzzzz_1, g_0_xxyyy_0_xxxxy_0, g_0_xxyyy_0_xxxxy_1, g_0_xxyyy_0_xxxyy_0, g_0_xxyyy_0_xxxyy_1, g_0_xxyyy_0_xxxyz_0, g_0_xxyyy_0_xxxyz_1, g_0_xxyyy_0_xxyyy_0, g_0_xxyyy_0_xxyyy_1, g_0_xxyyy_0_xxyyz_0, g_0_xxyyy_0_xxyyz_1, g_0_xxyyy_0_xxyzz_0, g_0_xxyyy_0_xxyzz_1, g_0_xxyyy_0_xyyyy_0, g_0_xxyyy_0_xyyyy_1, g_0_xxyyy_0_xyyyz_0, g_0_xxyyy_0_xyyyz_1, g_0_xxyyy_0_xyyzz_0, g_0_xxyyy_0_xyyzz_1, g_0_xxyyy_0_xyzzz_0, g_0_xxyyy_0_xyzzz_1, g_0_xxyyy_0_yyyyy_0, g_0_xxyyy_0_yyyyy_1, g_0_xxyyy_0_yyyyz_0, g_0_xxyyy_0_yyyyz_1, g_0_xxyyy_0_yyyzz_0, g_0_xxyyy_0_yyyzz_1, g_0_xxyyy_0_yyzzz_0, g_0_xxyyy_0_yyzzz_1, g_0_xxyyy_0_yzzzz_0, g_0_xxyyy_0_yzzzz_1, g_0_xxyyy_0_zzzzz_0, g_0_xxyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_xxxxx_0[i] = 2.0 * g_0_xxxxy_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxxx_0[i] * pb_y + g_0_xxxxyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxxy_0[i] = 3.0 * g_0_xxyyy_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxy_0[i] * pb_x + g_0_xxxyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxz_0[i] = 2.0 * g_0_xxxxy_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxxz_0[i] * pb_y + g_0_xxxxyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxyy_0[i] = 3.0 * g_0_xxyyy_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyy_0[i] * pb_x + g_0_xxxyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxyz_0[i] = 3.0 * g_0_xxyyy_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyz_0[i] * pb_x + g_0_xxxyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxzz_0[i] = 2.0 * g_0_xxxxy_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxzz_0[i] * pb_y + g_0_xxxxyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxyyy_0[i] = 3.0 * g_0_xxyyy_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyy_0[i] * pb_x + g_0_xxxyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyyz_0[i] = 3.0 * g_0_xxyyy_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyz_0[i] * pb_x + g_0_xxxyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyzz_0[i] = 3.0 * g_0_xxyyy_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyzz_0[i] * pb_x + g_0_xxxyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxzzz_0[i] = 2.0 * g_0_xxxxy_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxzzz_0[i] * pb_y + g_0_xxxxyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xyyyy_0[i] = 3.0 * g_0_xxyyy_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyy_0[i] * pb_x + g_0_xxxyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyyz_0[i] = 3.0 * g_0_xxyyy_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyz_0[i] * pb_x + g_0_xxxyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyzz_0[i] = 3.0 * g_0_xxyyy_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyzz_0[i] * pb_x + g_0_xxxyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyzzz_0[i] = 3.0 * g_0_xxyyy_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzzz_0[i] * pb_x + g_0_xxxyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xzzzz_0[i] = 2.0 * g_0_xxxxy_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xzzzz_0[i] * pb_y + g_0_xxxxyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_yyyyy_0[i] = 3.0 * g_0_xxyyy_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyy_0[i] * pb_x + g_0_xxxyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyyz_0[i] = 3.0 * g_0_xxyyy_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyz_0[i] * pb_x + g_0_xxxyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyzz_0[i] = 3.0 * g_0_xxyyy_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyzz_0[i] * pb_x + g_0_xxxyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyzzz_0[i] = 3.0 * g_0_xxyyy_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyzzz_0[i] * pb_x + g_0_xxxyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yzzzz_0[i] = 3.0 * g_0_xxyyy_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzzzz_0[i] * pb_x + g_0_xxxyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_zzzzz_0[i] = 3.0 * g_0_xxyyy_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_zzzzz_0[i] * pb_x + g_0_xxxyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 147-168 components of targeted buffer : SKSH

    auto g_0_xxxxyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 147);

    auto g_0_xxxxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 148);

    auto g_0_xxxxyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 149);

    auto g_0_xxxxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 150);

    auto g_0_xxxxyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 151);

    auto g_0_xxxxyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 152);

    auto g_0_xxxxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 153);

    auto g_0_xxxxyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 154);

    auto g_0_xxxxyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 155);

    auto g_0_xxxxyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 156);

    auto g_0_xxxxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 157);

    auto g_0_xxxxyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 158);

    auto g_0_xxxxyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 159);

    auto g_0_xxxxyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 160);

    auto g_0_xxxxyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 161);

    auto g_0_xxxxyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 162);

    auto g_0_xxxxyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 163);

    auto g_0_xxxxyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 164);

    auto g_0_xxxxyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 165);

    auto g_0_xxxxyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 166);

    auto g_0_xxxxyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 167);

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxx_1, g_0_xxxxyy_0_xxxxx_0, g_0_xxxxyy_0_xxxxx_1, g_0_xxxxyy_0_xxxxy_0, g_0_xxxxyy_0_xxxxy_1, g_0_xxxxyy_0_xxxxz_0, g_0_xxxxyy_0_xxxxz_1, g_0_xxxxyy_0_xxxy_1, g_0_xxxxyy_0_xxxyy_0, g_0_xxxxyy_0_xxxyy_1, g_0_xxxxyy_0_xxxyz_0, g_0_xxxxyy_0_xxxyz_1, g_0_xxxxyy_0_xxxz_1, g_0_xxxxyy_0_xxxzz_0, g_0_xxxxyy_0_xxxzz_1, g_0_xxxxyy_0_xxyy_1, g_0_xxxxyy_0_xxyyy_0, g_0_xxxxyy_0_xxyyy_1, g_0_xxxxyy_0_xxyyz_0, g_0_xxxxyy_0_xxyyz_1, g_0_xxxxyy_0_xxyz_1, g_0_xxxxyy_0_xxyzz_0, g_0_xxxxyy_0_xxyzz_1, g_0_xxxxyy_0_xxzz_1, g_0_xxxxyy_0_xxzzz_0, g_0_xxxxyy_0_xxzzz_1, g_0_xxxxyy_0_xyyy_1, g_0_xxxxyy_0_xyyyy_0, g_0_xxxxyy_0_xyyyy_1, g_0_xxxxyy_0_xyyyz_0, g_0_xxxxyy_0_xyyyz_1, g_0_xxxxyy_0_xyyz_1, g_0_xxxxyy_0_xyyzz_0, g_0_xxxxyy_0_xyyzz_1, g_0_xxxxyy_0_xyzz_1, g_0_xxxxyy_0_xyzzz_0, g_0_xxxxyy_0_xyzzz_1, g_0_xxxxyy_0_xzzz_1, g_0_xxxxyy_0_xzzzz_0, g_0_xxxxyy_0_xzzzz_1, g_0_xxxxyy_0_yyyy_1, g_0_xxxxyy_0_yyyyy_0, g_0_xxxxyy_0_yyyyy_1, g_0_xxxxyy_0_yyyyz_0, g_0_xxxxyy_0_yyyyz_1, g_0_xxxxyy_0_yyyz_1, g_0_xxxxyy_0_yyyzz_0, g_0_xxxxyy_0_yyyzz_1, g_0_xxxxyy_0_yyzz_1, g_0_xxxxyy_0_yyzzz_0, g_0_xxxxyy_0_yyzzz_1, g_0_xxxxyy_0_yzzz_1, g_0_xxxxyy_0_yzzzz_0, g_0_xxxxyy_0_yzzzz_1, g_0_xxxxyy_0_zzzz_1, g_0_xxxxyy_0_zzzzz_0, g_0_xxxxyy_0_zzzzz_1, g_0_xxxxyyz_0_xxxxx_0, g_0_xxxxyyz_0_xxxxy_0, g_0_xxxxyyz_0_xxxxz_0, g_0_xxxxyyz_0_xxxyy_0, g_0_xxxxyyz_0_xxxyz_0, g_0_xxxxyyz_0_xxxzz_0, g_0_xxxxyyz_0_xxyyy_0, g_0_xxxxyyz_0_xxyyz_0, g_0_xxxxyyz_0_xxyzz_0, g_0_xxxxyyz_0_xxzzz_0, g_0_xxxxyyz_0_xyyyy_0, g_0_xxxxyyz_0_xyyyz_0, g_0_xxxxyyz_0_xyyzz_0, g_0_xxxxyyz_0_xyzzz_0, g_0_xxxxyyz_0_xzzzz_0, g_0_xxxxyyz_0_yyyyy_0, g_0_xxxxyyz_0_yyyyz_0, g_0_xxxxyyz_0_yyyzz_0, g_0_xxxxyyz_0_yyzzz_0, g_0_xxxxyyz_0_yzzzz_0, g_0_xxxxyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_xxxxx_0[i] = g_0_xxxxyy_0_xxxxx_0[i] * pb_z + g_0_xxxxyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxy_0[i] = g_0_xxxxyy_0_xxxxy_0[i] * pb_z + g_0_xxxxyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxz_0[i] = g_0_xxxxyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxz_0[i] * pb_z + g_0_xxxxyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyy_0[i] = g_0_xxxxyy_0_xxxyy_0[i] * pb_z + g_0_xxxxyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyz_0[i] = g_0_xxxxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyz_0[i] * pb_z + g_0_xxxxyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxzz_0[i] = 2.0 * g_0_xxxxyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxzz_0[i] * pb_z + g_0_xxxxyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyy_0[i] = g_0_xxxxyy_0_xxyyy_0[i] * pb_z + g_0_xxxxyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyz_0[i] = g_0_xxxxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyz_0[i] * pb_z + g_0_xxxxyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyzz_0[i] = 2.0 * g_0_xxxxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyzz_0[i] * pb_z + g_0_xxxxyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxzzz_0[i] * pb_z + g_0_xxxxyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyy_0[i] = g_0_xxxxyy_0_xyyyy_0[i] * pb_z + g_0_xxxxyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyz_0[i] = g_0_xxxxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyz_0[i] * pb_z + g_0_xxxxyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyzz_0[i] = 2.0 * g_0_xxxxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyzz_0[i] * pb_z + g_0_xxxxyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyzzz_0[i] = 3.0 * g_0_xxxxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzzz_0[i] * pb_z + g_0_xxxxyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xzzzz_0[i] = 4.0 * g_0_xxxxyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xzzzz_0[i] * pb_z + g_0_xxxxyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyy_0[i] = g_0_xxxxyy_0_yyyyy_0[i] * pb_z + g_0_xxxxyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyz_0[i] = g_0_xxxxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyyz_0[i] * pb_z + g_0_xxxxyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyzz_0[i] = 2.0 * g_0_xxxxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyzz_0[i] * pb_z + g_0_xxxxyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyzzz_0[i] = 3.0 * g_0_xxxxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyzzz_0[i] * pb_z + g_0_xxxxyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yzzzz_0[i] = 4.0 * g_0_xxxxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yzzzz_0[i] * pb_z + g_0_xxxxyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_zzzzz_0[i] = 5.0 * g_0_xxxxyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_zzzzz_0[i] * pb_z + g_0_xxxxyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 168-189 components of targeted buffer : SKSH

    auto g_0_xxxxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 168);

    auto g_0_xxxxyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 169);

    auto g_0_xxxxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 170);

    auto g_0_xxxxyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 171);

    auto g_0_xxxxyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 172);

    auto g_0_xxxxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 173);

    auto g_0_xxxxyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 174);

    auto g_0_xxxxyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 175);

    auto g_0_xxxxyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 176);

    auto g_0_xxxxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 177);

    auto g_0_xxxxyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 178);

    auto g_0_xxxxyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 179);

    auto g_0_xxxxyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 180);

    auto g_0_xxxxyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 181);

    auto g_0_xxxxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 182);

    auto g_0_xxxxyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 183);

    auto g_0_xxxxyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 184);

    auto g_0_xxxxyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 185);

    auto g_0_xxxxyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 186);

    auto g_0_xxxxyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 187);

    auto g_0_xxxxyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 188);

    #pragma omp simd aligned(g_0_xxxxyzz_0_xxxxx_0, g_0_xxxxyzz_0_xxxxy_0, g_0_xxxxyzz_0_xxxxz_0, g_0_xxxxyzz_0_xxxyy_0, g_0_xxxxyzz_0_xxxyz_0, g_0_xxxxyzz_0_xxxzz_0, g_0_xxxxyzz_0_xxyyy_0, g_0_xxxxyzz_0_xxyyz_0, g_0_xxxxyzz_0_xxyzz_0, g_0_xxxxyzz_0_xxzzz_0, g_0_xxxxyzz_0_xyyyy_0, g_0_xxxxyzz_0_xyyyz_0, g_0_xxxxyzz_0_xyyzz_0, g_0_xxxxyzz_0_xyzzz_0, g_0_xxxxyzz_0_xzzzz_0, g_0_xxxxyzz_0_yyyyy_0, g_0_xxxxyzz_0_yyyyz_0, g_0_xxxxyzz_0_yyyzz_0, g_0_xxxxyzz_0_yyzzz_0, g_0_xxxxyzz_0_yzzzz_0, g_0_xxxxyzz_0_zzzzz_0, g_0_xxxxzz_0_xxxx_1, g_0_xxxxzz_0_xxxxx_0, g_0_xxxxzz_0_xxxxx_1, g_0_xxxxzz_0_xxxxy_0, g_0_xxxxzz_0_xxxxy_1, g_0_xxxxzz_0_xxxxz_0, g_0_xxxxzz_0_xxxxz_1, g_0_xxxxzz_0_xxxy_1, g_0_xxxxzz_0_xxxyy_0, g_0_xxxxzz_0_xxxyy_1, g_0_xxxxzz_0_xxxyz_0, g_0_xxxxzz_0_xxxyz_1, g_0_xxxxzz_0_xxxz_1, g_0_xxxxzz_0_xxxzz_0, g_0_xxxxzz_0_xxxzz_1, g_0_xxxxzz_0_xxyy_1, g_0_xxxxzz_0_xxyyy_0, g_0_xxxxzz_0_xxyyy_1, g_0_xxxxzz_0_xxyyz_0, g_0_xxxxzz_0_xxyyz_1, g_0_xxxxzz_0_xxyz_1, g_0_xxxxzz_0_xxyzz_0, g_0_xxxxzz_0_xxyzz_1, g_0_xxxxzz_0_xxzz_1, g_0_xxxxzz_0_xxzzz_0, g_0_xxxxzz_0_xxzzz_1, g_0_xxxxzz_0_xyyy_1, g_0_xxxxzz_0_xyyyy_0, g_0_xxxxzz_0_xyyyy_1, g_0_xxxxzz_0_xyyyz_0, g_0_xxxxzz_0_xyyyz_1, g_0_xxxxzz_0_xyyz_1, g_0_xxxxzz_0_xyyzz_0, g_0_xxxxzz_0_xyyzz_1, g_0_xxxxzz_0_xyzz_1, g_0_xxxxzz_0_xyzzz_0, g_0_xxxxzz_0_xyzzz_1, g_0_xxxxzz_0_xzzz_1, g_0_xxxxzz_0_xzzzz_0, g_0_xxxxzz_0_xzzzz_1, g_0_xxxxzz_0_yyyy_1, g_0_xxxxzz_0_yyyyy_0, g_0_xxxxzz_0_yyyyy_1, g_0_xxxxzz_0_yyyyz_0, g_0_xxxxzz_0_yyyyz_1, g_0_xxxxzz_0_yyyz_1, g_0_xxxxzz_0_yyyzz_0, g_0_xxxxzz_0_yyyzz_1, g_0_xxxxzz_0_yyzz_1, g_0_xxxxzz_0_yyzzz_0, g_0_xxxxzz_0_yyzzz_1, g_0_xxxxzz_0_yzzz_1, g_0_xxxxzz_0_yzzzz_0, g_0_xxxxzz_0_yzzzz_1, g_0_xxxxzz_0_zzzz_1, g_0_xxxxzz_0_zzzzz_0, g_0_xxxxzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_xxxxx_0[i] = g_0_xxxxzz_0_xxxxx_0[i] * pb_y + g_0_xxxxzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxy_0[i] = g_0_xxxxzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxy_0[i] * pb_y + g_0_xxxxzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxz_0[i] = g_0_xxxxzz_0_xxxxz_0[i] * pb_y + g_0_xxxxzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyy_0[i] = 2.0 * g_0_xxxxzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyy_0[i] * pb_y + g_0_xxxxzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyz_0[i] = g_0_xxxxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyz_0[i] * pb_y + g_0_xxxxzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxzz_0[i] = g_0_xxxxzz_0_xxxzz_0[i] * pb_y + g_0_xxxxzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyy_0[i] * pb_y + g_0_xxxxzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyz_0[i] = 2.0 * g_0_xxxxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyz_0[i] * pb_y + g_0_xxxxzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyzz_0[i] = g_0_xxxxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyzz_0[i] * pb_y + g_0_xxxxzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxzzz_0[i] = g_0_xxxxzz_0_xxzzz_0[i] * pb_y + g_0_xxxxzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyy_0[i] = 4.0 * g_0_xxxxzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyy_0[i] * pb_y + g_0_xxxxzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyz_0[i] = 3.0 * g_0_xxxxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyz_0[i] * pb_y + g_0_xxxxzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyzz_0[i] = 2.0 * g_0_xxxxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyzz_0[i] * pb_y + g_0_xxxxzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyzzz_0[i] = g_0_xxxxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzzz_0[i] * pb_y + g_0_xxxxzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xzzzz_0[i] = g_0_xxxxzz_0_xzzzz_0[i] * pb_y + g_0_xxxxzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyy_0[i] = 5.0 * g_0_xxxxzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyy_0[i] * pb_y + g_0_xxxxzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyz_0[i] = 4.0 * g_0_xxxxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyz_0[i] * pb_y + g_0_xxxxzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyzz_0[i] = 3.0 * g_0_xxxxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyzz_0[i] * pb_y + g_0_xxxxzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyzzz_0[i] = 2.0 * g_0_xxxxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyzzz_0[i] * pb_y + g_0_xxxxzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yzzzz_0[i] = g_0_xxxxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yzzzz_0[i] * pb_y + g_0_xxxxzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_zzzzz_0[i] = g_0_xxxxzz_0_zzzzz_0[i] * pb_y + g_0_xxxxzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 189-210 components of targeted buffer : SKSH

    auto g_0_xxxxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 189);

    auto g_0_xxxxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 190);

    auto g_0_xxxxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 191);

    auto g_0_xxxxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 192);

    auto g_0_xxxxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 193);

    auto g_0_xxxxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 194);

    auto g_0_xxxxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 195);

    auto g_0_xxxxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 196);

    auto g_0_xxxxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 197);

    auto g_0_xxxxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 198);

    auto g_0_xxxxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 199);

    auto g_0_xxxxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 200);

    auto g_0_xxxxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 201);

    auto g_0_xxxxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 202);

    auto g_0_xxxxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 203);

    auto g_0_xxxxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 204);

    auto g_0_xxxxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 205);

    auto g_0_xxxxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 206);

    auto g_0_xxxxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 207);

    auto g_0_xxxxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 208);

    auto g_0_xxxxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 209);

    #pragma omp simd aligned(g_0_xxxxz_0_xxxxx_0, g_0_xxxxz_0_xxxxx_1, g_0_xxxxz_0_xxxxy_0, g_0_xxxxz_0_xxxxy_1, g_0_xxxxz_0_xxxyy_0, g_0_xxxxz_0_xxxyy_1, g_0_xxxxz_0_xxyyy_0, g_0_xxxxz_0_xxyyy_1, g_0_xxxxz_0_xyyyy_0, g_0_xxxxz_0_xyyyy_1, g_0_xxxxzz_0_xxxxx_0, g_0_xxxxzz_0_xxxxx_1, g_0_xxxxzz_0_xxxxy_0, g_0_xxxxzz_0_xxxxy_1, g_0_xxxxzz_0_xxxyy_0, g_0_xxxxzz_0_xxxyy_1, g_0_xxxxzz_0_xxyyy_0, g_0_xxxxzz_0_xxyyy_1, g_0_xxxxzz_0_xyyyy_0, g_0_xxxxzz_0_xyyyy_1, g_0_xxxxzzz_0_xxxxx_0, g_0_xxxxzzz_0_xxxxy_0, g_0_xxxxzzz_0_xxxxz_0, g_0_xxxxzzz_0_xxxyy_0, g_0_xxxxzzz_0_xxxyz_0, g_0_xxxxzzz_0_xxxzz_0, g_0_xxxxzzz_0_xxyyy_0, g_0_xxxxzzz_0_xxyyz_0, g_0_xxxxzzz_0_xxyzz_0, g_0_xxxxzzz_0_xxzzz_0, g_0_xxxxzzz_0_xyyyy_0, g_0_xxxxzzz_0_xyyyz_0, g_0_xxxxzzz_0_xyyzz_0, g_0_xxxxzzz_0_xyzzz_0, g_0_xxxxzzz_0_xzzzz_0, g_0_xxxxzzz_0_yyyyy_0, g_0_xxxxzzz_0_yyyyz_0, g_0_xxxxzzz_0_yyyzz_0, g_0_xxxxzzz_0_yyzzz_0, g_0_xxxxzzz_0_yzzzz_0, g_0_xxxxzzz_0_zzzzz_0, g_0_xxxzzz_0_xxxxz_0, g_0_xxxzzz_0_xxxxz_1, g_0_xxxzzz_0_xxxyz_0, g_0_xxxzzz_0_xxxyz_1, g_0_xxxzzz_0_xxxz_1, g_0_xxxzzz_0_xxxzz_0, g_0_xxxzzz_0_xxxzz_1, g_0_xxxzzz_0_xxyyz_0, g_0_xxxzzz_0_xxyyz_1, g_0_xxxzzz_0_xxyz_1, g_0_xxxzzz_0_xxyzz_0, g_0_xxxzzz_0_xxyzz_1, g_0_xxxzzz_0_xxzz_1, g_0_xxxzzz_0_xxzzz_0, g_0_xxxzzz_0_xxzzz_1, g_0_xxxzzz_0_xyyyz_0, g_0_xxxzzz_0_xyyyz_1, g_0_xxxzzz_0_xyyz_1, g_0_xxxzzz_0_xyyzz_0, g_0_xxxzzz_0_xyyzz_1, g_0_xxxzzz_0_xyzz_1, g_0_xxxzzz_0_xyzzz_0, g_0_xxxzzz_0_xyzzz_1, g_0_xxxzzz_0_xzzz_1, g_0_xxxzzz_0_xzzzz_0, g_0_xxxzzz_0_xzzzz_1, g_0_xxxzzz_0_yyyyy_0, g_0_xxxzzz_0_yyyyy_1, g_0_xxxzzz_0_yyyyz_0, g_0_xxxzzz_0_yyyyz_1, g_0_xxxzzz_0_yyyz_1, g_0_xxxzzz_0_yyyzz_0, g_0_xxxzzz_0_yyyzz_1, g_0_xxxzzz_0_yyzz_1, g_0_xxxzzz_0_yyzzz_0, g_0_xxxzzz_0_yyzzz_1, g_0_xxxzzz_0_yzzz_1, g_0_xxxzzz_0_yzzzz_0, g_0_xxxzzz_0_yzzzz_1, g_0_xxxzzz_0_zzzz_1, g_0_xxxzzz_0_zzzzz_0, g_0_xxxzzz_0_zzzzz_1, g_0_xxzzz_0_xxxxz_0, g_0_xxzzz_0_xxxxz_1, g_0_xxzzz_0_xxxyz_0, g_0_xxzzz_0_xxxyz_1, g_0_xxzzz_0_xxxzz_0, g_0_xxzzz_0_xxxzz_1, g_0_xxzzz_0_xxyyz_0, g_0_xxzzz_0_xxyyz_1, g_0_xxzzz_0_xxyzz_0, g_0_xxzzz_0_xxyzz_1, g_0_xxzzz_0_xxzzz_0, g_0_xxzzz_0_xxzzz_1, g_0_xxzzz_0_xyyyz_0, g_0_xxzzz_0_xyyyz_1, g_0_xxzzz_0_xyyzz_0, g_0_xxzzz_0_xyyzz_1, g_0_xxzzz_0_xyzzz_0, g_0_xxzzz_0_xyzzz_1, g_0_xxzzz_0_xzzzz_0, g_0_xxzzz_0_xzzzz_1, g_0_xxzzz_0_yyyyy_0, g_0_xxzzz_0_yyyyy_1, g_0_xxzzz_0_yyyyz_0, g_0_xxzzz_0_yyyyz_1, g_0_xxzzz_0_yyyzz_0, g_0_xxzzz_0_yyyzz_1, g_0_xxzzz_0_yyzzz_0, g_0_xxzzz_0_yyzzz_1, g_0_xxzzz_0_yzzzz_0, g_0_xxzzz_0_yzzzz_1, g_0_xxzzz_0_zzzzz_0, g_0_xxzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_xxxxx_0[i] = 2.0 * g_0_xxxxz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxxx_0[i] * pb_z + g_0_xxxxzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxy_0[i] = 2.0 * g_0_xxxxz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxxy_0[i] * pb_z + g_0_xxxxzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxz_0[i] = 3.0 * g_0_xxzzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxz_0[i] * pb_x + g_0_xxxzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxxz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxyy_0[i] * pb_z + g_0_xxxxzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxyz_0[i] = 3.0 * g_0_xxzzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyz_0[i] * pb_x + g_0_xxxzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxzz_0[i] = 3.0 * g_0_xxzzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxzz_0[i] * pb_x + g_0_xxxzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyyy_0[i] = 2.0 * g_0_xxxxz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxyyy_0[i] * pb_z + g_0_xxxxzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxyyz_0[i] = 3.0 * g_0_xxzzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyz_0[i] * pb_x + g_0_xxxzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyzz_0[i] = 3.0 * g_0_xxzzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyzz_0[i] * pb_x + g_0_xxxzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxzzz_0[i] = 3.0 * g_0_xxzzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxzzz_0[i] * pb_x + g_0_xxxzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyyy_0[i] = 2.0 * g_0_xxxxz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xyyyy_0[i] * pb_z + g_0_xxxxzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xyyyz_0[i] = 3.0 * g_0_xxzzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyz_0[i] * pb_x + g_0_xxxzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyzz_0[i] = 3.0 * g_0_xxzzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyzz_0[i] * pb_x + g_0_xxxzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyzzz_0[i] = 3.0 * g_0_xxzzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzzz_0[i] * pb_x + g_0_xxxzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xzzzz_0[i] = 3.0 * g_0_xxzzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xzzzz_0[i] * pb_x + g_0_xxxzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyy_0[i] = 3.0 * g_0_xxzzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyy_0[i] * pb_x + g_0_xxxzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyz_0[i] = 3.0 * g_0_xxzzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyz_0[i] * pb_x + g_0_xxxzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyzz_0[i] = 3.0 * g_0_xxzzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyzz_0[i] * pb_x + g_0_xxxzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyzzz_0[i] = 3.0 * g_0_xxzzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyzzz_0[i] * pb_x + g_0_xxxzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yzzzz_0[i] = 3.0 * g_0_xxzzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzzzz_0[i] * pb_x + g_0_xxxzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_zzzzz_0[i] = 3.0 * g_0_xxzzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzzzz_0[i] * pb_x + g_0_xxxzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 210-231 components of targeted buffer : SKSH

    auto g_0_xxxyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 210);

    auto g_0_xxxyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 211);

    auto g_0_xxxyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 212);

    auto g_0_xxxyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 213);

    auto g_0_xxxyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 214);

    auto g_0_xxxyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 215);

    auto g_0_xxxyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 216);

    auto g_0_xxxyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 217);

    auto g_0_xxxyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 218);

    auto g_0_xxxyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 219);

    auto g_0_xxxyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 220);

    auto g_0_xxxyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 221);

    auto g_0_xxxyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 222);

    auto g_0_xxxyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 223);

    auto g_0_xxxyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 224);

    auto g_0_xxxyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 225);

    auto g_0_xxxyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 226);

    auto g_0_xxxyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 227);

    auto g_0_xxxyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 228);

    auto g_0_xxxyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 229);

    auto g_0_xxxyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 230);

    #pragma omp simd aligned(g_0_xxxyy_0_xxxxx_0, g_0_xxxyy_0_xxxxx_1, g_0_xxxyy_0_xxxxz_0, g_0_xxxyy_0_xxxxz_1, g_0_xxxyy_0_xxxzz_0, g_0_xxxyy_0_xxxzz_1, g_0_xxxyy_0_xxzzz_0, g_0_xxxyy_0_xxzzz_1, g_0_xxxyy_0_xzzzz_0, g_0_xxxyy_0_xzzzz_1, g_0_xxxyyy_0_xxxxx_0, g_0_xxxyyy_0_xxxxx_1, g_0_xxxyyy_0_xxxxz_0, g_0_xxxyyy_0_xxxxz_1, g_0_xxxyyy_0_xxxzz_0, g_0_xxxyyy_0_xxxzz_1, g_0_xxxyyy_0_xxzzz_0, g_0_xxxyyy_0_xxzzz_1, g_0_xxxyyy_0_xzzzz_0, g_0_xxxyyy_0_xzzzz_1, g_0_xxxyyyy_0_xxxxx_0, g_0_xxxyyyy_0_xxxxy_0, g_0_xxxyyyy_0_xxxxz_0, g_0_xxxyyyy_0_xxxyy_0, g_0_xxxyyyy_0_xxxyz_0, g_0_xxxyyyy_0_xxxzz_0, g_0_xxxyyyy_0_xxyyy_0, g_0_xxxyyyy_0_xxyyz_0, g_0_xxxyyyy_0_xxyzz_0, g_0_xxxyyyy_0_xxzzz_0, g_0_xxxyyyy_0_xyyyy_0, g_0_xxxyyyy_0_xyyyz_0, g_0_xxxyyyy_0_xyyzz_0, g_0_xxxyyyy_0_xyzzz_0, g_0_xxxyyyy_0_xzzzz_0, g_0_xxxyyyy_0_yyyyy_0, g_0_xxxyyyy_0_yyyyz_0, g_0_xxxyyyy_0_yyyzz_0, g_0_xxxyyyy_0_yyzzz_0, g_0_xxxyyyy_0_yzzzz_0, g_0_xxxyyyy_0_zzzzz_0, g_0_xxyyyy_0_xxxxy_0, g_0_xxyyyy_0_xxxxy_1, g_0_xxyyyy_0_xxxy_1, g_0_xxyyyy_0_xxxyy_0, g_0_xxyyyy_0_xxxyy_1, g_0_xxyyyy_0_xxxyz_0, g_0_xxyyyy_0_xxxyz_1, g_0_xxyyyy_0_xxyy_1, g_0_xxyyyy_0_xxyyy_0, g_0_xxyyyy_0_xxyyy_1, g_0_xxyyyy_0_xxyyz_0, g_0_xxyyyy_0_xxyyz_1, g_0_xxyyyy_0_xxyz_1, g_0_xxyyyy_0_xxyzz_0, g_0_xxyyyy_0_xxyzz_1, g_0_xxyyyy_0_xyyy_1, g_0_xxyyyy_0_xyyyy_0, g_0_xxyyyy_0_xyyyy_1, g_0_xxyyyy_0_xyyyz_0, g_0_xxyyyy_0_xyyyz_1, g_0_xxyyyy_0_xyyz_1, g_0_xxyyyy_0_xyyzz_0, g_0_xxyyyy_0_xyyzz_1, g_0_xxyyyy_0_xyzz_1, g_0_xxyyyy_0_xyzzz_0, g_0_xxyyyy_0_xyzzz_1, g_0_xxyyyy_0_yyyy_1, g_0_xxyyyy_0_yyyyy_0, g_0_xxyyyy_0_yyyyy_1, g_0_xxyyyy_0_yyyyz_0, g_0_xxyyyy_0_yyyyz_1, g_0_xxyyyy_0_yyyz_1, g_0_xxyyyy_0_yyyzz_0, g_0_xxyyyy_0_yyyzz_1, g_0_xxyyyy_0_yyzz_1, g_0_xxyyyy_0_yyzzz_0, g_0_xxyyyy_0_yyzzz_1, g_0_xxyyyy_0_yzzz_1, g_0_xxyyyy_0_yzzzz_0, g_0_xxyyyy_0_yzzzz_1, g_0_xxyyyy_0_zzzzz_0, g_0_xxyyyy_0_zzzzz_1, g_0_xyyyy_0_xxxxy_0, g_0_xyyyy_0_xxxxy_1, g_0_xyyyy_0_xxxyy_0, g_0_xyyyy_0_xxxyy_1, g_0_xyyyy_0_xxxyz_0, g_0_xyyyy_0_xxxyz_1, g_0_xyyyy_0_xxyyy_0, g_0_xyyyy_0_xxyyy_1, g_0_xyyyy_0_xxyyz_0, g_0_xyyyy_0_xxyyz_1, g_0_xyyyy_0_xxyzz_0, g_0_xyyyy_0_xxyzz_1, g_0_xyyyy_0_xyyyy_0, g_0_xyyyy_0_xyyyy_1, g_0_xyyyy_0_xyyyz_0, g_0_xyyyy_0_xyyyz_1, g_0_xyyyy_0_xyyzz_0, g_0_xyyyy_0_xyyzz_1, g_0_xyyyy_0_xyzzz_0, g_0_xyyyy_0_xyzzz_1, g_0_xyyyy_0_yyyyy_0, g_0_xyyyy_0_yyyyy_1, g_0_xyyyy_0_yyyyz_0, g_0_xyyyy_0_yyyyz_1, g_0_xyyyy_0_yyyzz_0, g_0_xyyyy_0_yyyzz_1, g_0_xyyyy_0_yyzzz_0, g_0_xyyyy_0_yyzzz_1, g_0_xyyyy_0_yzzzz_0, g_0_xyyyy_0_yzzzz_1, g_0_xyyyy_0_zzzzz_0, g_0_xyyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_xxxxx_0[i] = 3.0 * g_0_xxxyy_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxxx_0[i] * pb_y + g_0_xxxyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxxy_0[i] = 2.0 * g_0_xyyyy_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxy_0[i] * pb_x + g_0_xxyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxz_0[i] = 3.0 * g_0_xxxyy_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxxz_0[i] * pb_y + g_0_xxxyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxyy_0[i] = 2.0 * g_0_xyyyy_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyy_0[i] * pb_x + g_0_xxyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxyz_0[i] = 2.0 * g_0_xyyyy_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyz_0[i] * pb_x + g_0_xxyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxzz_0[i] = 3.0 * g_0_xxxyy_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxzz_0[i] * pb_y + g_0_xxxyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxyyy_0[i] = 2.0 * g_0_xyyyy_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyy_0[i] * pb_x + g_0_xxyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyyz_0[i] = 2.0 * g_0_xyyyy_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyz_0[i] * pb_x + g_0_xxyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyzz_0[i] = 2.0 * g_0_xyyyy_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyzz_0[i] * pb_x + g_0_xxyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxzzz_0[i] = 3.0 * g_0_xxxyy_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxzzz_0[i] * pb_y + g_0_xxxyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xyyyy_0[i] = 2.0 * g_0_xyyyy_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyy_0[i] * pb_x + g_0_xxyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyyz_0[i] = 2.0 * g_0_xyyyy_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyz_0[i] * pb_x + g_0_xxyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyzz_0[i] = 2.0 * g_0_xyyyy_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyzz_0[i] * pb_x + g_0_xxyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyzzz_0[i] = 2.0 * g_0_xyyyy_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzzz_0[i] * pb_x + g_0_xxyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xzzzz_0[i] = 3.0 * g_0_xxxyy_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xzzzz_0[i] * pb_y + g_0_xxxyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_yyyyy_0[i] = 2.0 * g_0_xyyyy_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyy_0[i] * pb_x + g_0_xxyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyyz_0[i] = 2.0 * g_0_xyyyy_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyz_0[i] * pb_x + g_0_xxyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyzz_0[i] = 2.0 * g_0_xyyyy_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyzz_0[i] * pb_x + g_0_xxyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyzzz_0[i] = 2.0 * g_0_xyyyy_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyzzz_0[i] * pb_x + g_0_xxyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yzzzz_0[i] = 2.0 * g_0_xyyyy_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzzzz_0[i] * pb_x + g_0_xxyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_zzzzz_0[i] = 2.0 * g_0_xyyyy_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_zzzzz_0[i] * pb_x + g_0_xxyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 231-252 components of targeted buffer : SKSH

    auto g_0_xxxyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 231);

    auto g_0_xxxyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 232);

    auto g_0_xxxyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 233);

    auto g_0_xxxyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 234);

    auto g_0_xxxyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 235);

    auto g_0_xxxyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 236);

    auto g_0_xxxyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 237);

    auto g_0_xxxyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 238);

    auto g_0_xxxyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 239);

    auto g_0_xxxyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 240);

    auto g_0_xxxyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 241);

    auto g_0_xxxyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 242);

    auto g_0_xxxyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 243);

    auto g_0_xxxyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 244);

    auto g_0_xxxyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 245);

    auto g_0_xxxyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 246);

    auto g_0_xxxyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 247);

    auto g_0_xxxyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 248);

    auto g_0_xxxyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 249);

    auto g_0_xxxyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 250);

    auto g_0_xxxyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 251);

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxx_1, g_0_xxxyyy_0_xxxxx_0, g_0_xxxyyy_0_xxxxx_1, g_0_xxxyyy_0_xxxxy_0, g_0_xxxyyy_0_xxxxy_1, g_0_xxxyyy_0_xxxxz_0, g_0_xxxyyy_0_xxxxz_1, g_0_xxxyyy_0_xxxy_1, g_0_xxxyyy_0_xxxyy_0, g_0_xxxyyy_0_xxxyy_1, g_0_xxxyyy_0_xxxyz_0, g_0_xxxyyy_0_xxxyz_1, g_0_xxxyyy_0_xxxz_1, g_0_xxxyyy_0_xxxzz_0, g_0_xxxyyy_0_xxxzz_1, g_0_xxxyyy_0_xxyy_1, g_0_xxxyyy_0_xxyyy_0, g_0_xxxyyy_0_xxyyy_1, g_0_xxxyyy_0_xxyyz_0, g_0_xxxyyy_0_xxyyz_1, g_0_xxxyyy_0_xxyz_1, g_0_xxxyyy_0_xxyzz_0, g_0_xxxyyy_0_xxyzz_1, g_0_xxxyyy_0_xxzz_1, g_0_xxxyyy_0_xxzzz_0, g_0_xxxyyy_0_xxzzz_1, g_0_xxxyyy_0_xyyy_1, g_0_xxxyyy_0_xyyyy_0, g_0_xxxyyy_0_xyyyy_1, g_0_xxxyyy_0_xyyyz_0, g_0_xxxyyy_0_xyyyz_1, g_0_xxxyyy_0_xyyz_1, g_0_xxxyyy_0_xyyzz_0, g_0_xxxyyy_0_xyyzz_1, g_0_xxxyyy_0_xyzz_1, g_0_xxxyyy_0_xyzzz_0, g_0_xxxyyy_0_xyzzz_1, g_0_xxxyyy_0_xzzz_1, g_0_xxxyyy_0_xzzzz_0, g_0_xxxyyy_0_xzzzz_1, g_0_xxxyyy_0_yyyy_1, g_0_xxxyyy_0_yyyyy_0, g_0_xxxyyy_0_yyyyy_1, g_0_xxxyyy_0_yyyyz_0, g_0_xxxyyy_0_yyyyz_1, g_0_xxxyyy_0_yyyz_1, g_0_xxxyyy_0_yyyzz_0, g_0_xxxyyy_0_yyyzz_1, g_0_xxxyyy_0_yyzz_1, g_0_xxxyyy_0_yyzzz_0, g_0_xxxyyy_0_yyzzz_1, g_0_xxxyyy_0_yzzz_1, g_0_xxxyyy_0_yzzzz_0, g_0_xxxyyy_0_yzzzz_1, g_0_xxxyyy_0_zzzz_1, g_0_xxxyyy_0_zzzzz_0, g_0_xxxyyy_0_zzzzz_1, g_0_xxxyyyz_0_xxxxx_0, g_0_xxxyyyz_0_xxxxy_0, g_0_xxxyyyz_0_xxxxz_0, g_0_xxxyyyz_0_xxxyy_0, g_0_xxxyyyz_0_xxxyz_0, g_0_xxxyyyz_0_xxxzz_0, g_0_xxxyyyz_0_xxyyy_0, g_0_xxxyyyz_0_xxyyz_0, g_0_xxxyyyz_0_xxyzz_0, g_0_xxxyyyz_0_xxzzz_0, g_0_xxxyyyz_0_xyyyy_0, g_0_xxxyyyz_0_xyyyz_0, g_0_xxxyyyz_0_xyyzz_0, g_0_xxxyyyz_0_xyzzz_0, g_0_xxxyyyz_0_xzzzz_0, g_0_xxxyyyz_0_yyyyy_0, g_0_xxxyyyz_0_yyyyz_0, g_0_xxxyyyz_0_yyyzz_0, g_0_xxxyyyz_0_yyzzz_0, g_0_xxxyyyz_0_yzzzz_0, g_0_xxxyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_xxxxx_0[i] = g_0_xxxyyy_0_xxxxx_0[i] * pb_z + g_0_xxxyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxy_0[i] = g_0_xxxyyy_0_xxxxy_0[i] * pb_z + g_0_xxxyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxz_0[i] = g_0_xxxyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxz_0[i] * pb_z + g_0_xxxyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyy_0[i] = g_0_xxxyyy_0_xxxyy_0[i] * pb_z + g_0_xxxyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyz_0[i] = g_0_xxxyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyz_0[i] * pb_z + g_0_xxxyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxzz_0[i] = 2.0 * g_0_xxxyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxzz_0[i] * pb_z + g_0_xxxyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyy_0[i] = g_0_xxxyyy_0_xxyyy_0[i] * pb_z + g_0_xxxyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyz_0[i] = g_0_xxxyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyz_0[i] * pb_z + g_0_xxxyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyzz_0[i] = 2.0 * g_0_xxxyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyzz_0[i] * pb_z + g_0_xxxyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxzzz_0[i] = 3.0 * g_0_xxxyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxzzz_0[i] * pb_z + g_0_xxxyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyy_0[i] = g_0_xxxyyy_0_xyyyy_0[i] * pb_z + g_0_xxxyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyz_0[i] = g_0_xxxyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyz_0[i] * pb_z + g_0_xxxyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyzz_0[i] = 2.0 * g_0_xxxyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyzz_0[i] * pb_z + g_0_xxxyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyzzz_0[i] = 3.0 * g_0_xxxyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzzz_0[i] * pb_z + g_0_xxxyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xzzzz_0[i] * pb_z + g_0_xxxyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyy_0[i] = g_0_xxxyyy_0_yyyyy_0[i] * pb_z + g_0_xxxyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyz_0[i] = g_0_xxxyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyyz_0[i] * pb_z + g_0_xxxyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyzz_0[i] = 2.0 * g_0_xxxyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyzz_0[i] * pb_z + g_0_xxxyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyzzz_0[i] = 3.0 * g_0_xxxyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyzzz_0[i] * pb_z + g_0_xxxyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yzzzz_0[i] * pb_z + g_0_xxxyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_zzzzz_0[i] = 5.0 * g_0_xxxyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_zzzzz_0[i] * pb_z + g_0_xxxyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 252-273 components of targeted buffer : SKSH

    auto g_0_xxxyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 252);

    auto g_0_xxxyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 253);

    auto g_0_xxxyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 254);

    auto g_0_xxxyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 255);

    auto g_0_xxxyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 256);

    auto g_0_xxxyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 257);

    auto g_0_xxxyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 258);

    auto g_0_xxxyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 259);

    auto g_0_xxxyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 260);

    auto g_0_xxxyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 261);

    auto g_0_xxxyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 262);

    auto g_0_xxxyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 263);

    auto g_0_xxxyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 264);

    auto g_0_xxxyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 265);

    auto g_0_xxxyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 266);

    auto g_0_xxxyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 267);

    auto g_0_xxxyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 268);

    auto g_0_xxxyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 269);

    auto g_0_xxxyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 270);

    auto g_0_xxxyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 271);

    auto g_0_xxxyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 272);

    #pragma omp simd aligned(g_0_xxxyy_0_xxxxy_0, g_0_xxxyy_0_xxxxy_1, g_0_xxxyy_0_xxxyy_0, g_0_xxxyy_0_xxxyy_1, g_0_xxxyy_0_xxyyy_0, g_0_xxxyy_0_xxyyy_1, g_0_xxxyy_0_xyyyy_0, g_0_xxxyy_0_xyyyy_1, g_0_xxxyyz_0_xxxxy_0, g_0_xxxyyz_0_xxxxy_1, g_0_xxxyyz_0_xxxyy_0, g_0_xxxyyz_0_xxxyy_1, g_0_xxxyyz_0_xxyyy_0, g_0_xxxyyz_0_xxyyy_1, g_0_xxxyyz_0_xyyyy_0, g_0_xxxyyz_0_xyyyy_1, g_0_xxxyyzz_0_xxxxx_0, g_0_xxxyyzz_0_xxxxy_0, g_0_xxxyyzz_0_xxxxz_0, g_0_xxxyyzz_0_xxxyy_0, g_0_xxxyyzz_0_xxxyz_0, g_0_xxxyyzz_0_xxxzz_0, g_0_xxxyyzz_0_xxyyy_0, g_0_xxxyyzz_0_xxyyz_0, g_0_xxxyyzz_0_xxyzz_0, g_0_xxxyyzz_0_xxzzz_0, g_0_xxxyyzz_0_xyyyy_0, g_0_xxxyyzz_0_xyyyz_0, g_0_xxxyyzz_0_xyyzz_0, g_0_xxxyyzz_0_xyzzz_0, g_0_xxxyyzz_0_xzzzz_0, g_0_xxxyyzz_0_yyyyy_0, g_0_xxxyyzz_0_yyyyz_0, g_0_xxxyyzz_0_yyyzz_0, g_0_xxxyyzz_0_yyzzz_0, g_0_xxxyyzz_0_yzzzz_0, g_0_xxxyyzz_0_zzzzz_0, g_0_xxxyzz_0_xxxxx_0, g_0_xxxyzz_0_xxxxx_1, g_0_xxxyzz_0_xxxxz_0, g_0_xxxyzz_0_xxxxz_1, g_0_xxxyzz_0_xxxzz_0, g_0_xxxyzz_0_xxxzz_1, g_0_xxxyzz_0_xxzzz_0, g_0_xxxyzz_0_xxzzz_1, g_0_xxxyzz_0_xzzzz_0, g_0_xxxyzz_0_xzzzz_1, g_0_xxxzz_0_xxxxx_0, g_0_xxxzz_0_xxxxx_1, g_0_xxxzz_0_xxxxz_0, g_0_xxxzz_0_xxxxz_1, g_0_xxxzz_0_xxxzz_0, g_0_xxxzz_0_xxxzz_1, g_0_xxxzz_0_xxzzz_0, g_0_xxxzz_0_xxzzz_1, g_0_xxxzz_0_xzzzz_0, g_0_xxxzz_0_xzzzz_1, g_0_xxyyzz_0_xxxyz_0, g_0_xxyyzz_0_xxxyz_1, g_0_xxyyzz_0_xxyyz_0, g_0_xxyyzz_0_xxyyz_1, g_0_xxyyzz_0_xxyz_1, g_0_xxyyzz_0_xxyzz_0, g_0_xxyyzz_0_xxyzz_1, g_0_xxyyzz_0_xyyyz_0, g_0_xxyyzz_0_xyyyz_1, g_0_xxyyzz_0_xyyz_1, g_0_xxyyzz_0_xyyzz_0, g_0_xxyyzz_0_xyyzz_1, g_0_xxyyzz_0_xyzz_1, g_0_xxyyzz_0_xyzzz_0, g_0_xxyyzz_0_xyzzz_1, g_0_xxyyzz_0_yyyyy_0, g_0_xxyyzz_0_yyyyy_1, g_0_xxyyzz_0_yyyyz_0, g_0_xxyyzz_0_yyyyz_1, g_0_xxyyzz_0_yyyz_1, g_0_xxyyzz_0_yyyzz_0, g_0_xxyyzz_0_yyyzz_1, g_0_xxyyzz_0_yyzz_1, g_0_xxyyzz_0_yyzzz_0, g_0_xxyyzz_0_yyzzz_1, g_0_xxyyzz_0_yzzz_1, g_0_xxyyzz_0_yzzzz_0, g_0_xxyyzz_0_yzzzz_1, g_0_xxyyzz_0_zzzzz_0, g_0_xxyyzz_0_zzzzz_1, g_0_xyyzz_0_xxxyz_0, g_0_xyyzz_0_xxxyz_1, g_0_xyyzz_0_xxyyz_0, g_0_xyyzz_0_xxyyz_1, g_0_xyyzz_0_xxyzz_0, g_0_xyyzz_0_xxyzz_1, g_0_xyyzz_0_xyyyz_0, g_0_xyyzz_0_xyyyz_1, g_0_xyyzz_0_xyyzz_0, g_0_xyyzz_0_xyyzz_1, g_0_xyyzz_0_xyzzz_0, g_0_xyyzz_0_xyzzz_1, g_0_xyyzz_0_yyyyy_0, g_0_xyyzz_0_yyyyy_1, g_0_xyyzz_0_yyyyz_0, g_0_xyyzz_0_yyyyz_1, g_0_xyyzz_0_yyyzz_0, g_0_xyyzz_0_yyyzz_1, g_0_xyyzz_0_yyzzz_0, g_0_xyyzz_0_yyzzz_1, g_0_xyyzz_0_yzzzz_0, g_0_xyyzz_0_yzzzz_1, g_0_xyyzz_0_zzzzz_0, g_0_xyyzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzz_0_xxxxx_0[i] = g_0_xxxzz_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxx_0[i] * pb_y + g_0_xxxyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxxy_0[i] = g_0_xxxyy_0_xxxxy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxxy_0[i] * pb_z + g_0_xxxyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxxz_0[i] = g_0_xxxzz_0_xxxxz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxz_0[i] * pb_y + g_0_xxxyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxyy_0[i] = g_0_xxxyy_0_xxxyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxyy_0[i] * pb_z + g_0_xxxyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxyz_0[i] = 2.0 * g_0_xyyzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxyz_0[i] * pb_x + g_0_xxyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxzz_0[i] = g_0_xxxzz_0_xxxzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxzz_0[i] * pb_y + g_0_xxxyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxyyy_0[i] = g_0_xxxyy_0_xxyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxyyy_0[i] * pb_z + g_0_xxxyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxyyz_0[i] = 2.0 * g_0_xyyzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyyz_0[i] * pb_x + g_0_xxyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxyzz_0[i] = 2.0 * g_0_xyyzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyzz_0[i] * pb_x + g_0_xxyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxzzz_0[i] = g_0_xxxzz_0_xxzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxzzz_0[i] * pb_y + g_0_xxxyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xyyyy_0[i] = g_0_xxxyy_0_xyyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xyyyy_0[i] * pb_z + g_0_xxxyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xyyyz_0[i] = 2.0 * g_0_xyyzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyyz_0[i] * pb_x + g_0_xxyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyyzz_0[i] = 2.0 * g_0_xyyzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyzz_0[i] * pb_x + g_0_xxyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyzzz_0[i] = 2.0 * g_0_xyyzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyzzz_0[i] * pb_x + g_0_xxyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xzzzz_0[i] = g_0_xxxzz_0_xzzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xzzzz_0[i] * pb_y + g_0_xxxyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_yyyyy_0[i] = 2.0 * g_0_xyyzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyy_0[i] * pb_x + g_0_xxyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyyz_0[i] = 2.0 * g_0_xyyzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyz_0[i] * pb_x + g_0_xxyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyzz_0[i] = 2.0 * g_0_xyyzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyzz_0[i] * pb_x + g_0_xxyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyzzz_0[i] = 2.0 * g_0_xyyzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyzzz_0[i] * pb_x + g_0_xxyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yzzzz_0[i] = 2.0 * g_0_xyyzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzzzz_0[i] * pb_x + g_0_xxyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_zzzzz_0[i] = 2.0 * g_0_xyyzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_zzzzz_0[i] * pb_x + g_0_xxyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 273-294 components of targeted buffer : SKSH

    auto g_0_xxxyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 273);

    auto g_0_xxxyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 274);

    auto g_0_xxxyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 275);

    auto g_0_xxxyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 276);

    auto g_0_xxxyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 277);

    auto g_0_xxxyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 278);

    auto g_0_xxxyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 279);

    auto g_0_xxxyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 280);

    auto g_0_xxxyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 281);

    auto g_0_xxxyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 282);

    auto g_0_xxxyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 283);

    auto g_0_xxxyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 284);

    auto g_0_xxxyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 285);

    auto g_0_xxxyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 286);

    auto g_0_xxxyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 287);

    auto g_0_xxxyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 288);

    auto g_0_xxxyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 289);

    auto g_0_xxxyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 290);

    auto g_0_xxxyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 291);

    auto g_0_xxxyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 292);

    auto g_0_xxxyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 293);

    #pragma omp simd aligned(g_0_xxxyzzz_0_xxxxx_0, g_0_xxxyzzz_0_xxxxy_0, g_0_xxxyzzz_0_xxxxz_0, g_0_xxxyzzz_0_xxxyy_0, g_0_xxxyzzz_0_xxxyz_0, g_0_xxxyzzz_0_xxxzz_0, g_0_xxxyzzz_0_xxyyy_0, g_0_xxxyzzz_0_xxyyz_0, g_0_xxxyzzz_0_xxyzz_0, g_0_xxxyzzz_0_xxzzz_0, g_0_xxxyzzz_0_xyyyy_0, g_0_xxxyzzz_0_xyyyz_0, g_0_xxxyzzz_0_xyyzz_0, g_0_xxxyzzz_0_xyzzz_0, g_0_xxxyzzz_0_xzzzz_0, g_0_xxxyzzz_0_yyyyy_0, g_0_xxxyzzz_0_yyyyz_0, g_0_xxxyzzz_0_yyyzz_0, g_0_xxxyzzz_0_yyzzz_0, g_0_xxxyzzz_0_yzzzz_0, g_0_xxxyzzz_0_zzzzz_0, g_0_xxxzzz_0_xxxx_1, g_0_xxxzzz_0_xxxxx_0, g_0_xxxzzz_0_xxxxx_1, g_0_xxxzzz_0_xxxxy_0, g_0_xxxzzz_0_xxxxy_1, g_0_xxxzzz_0_xxxxz_0, g_0_xxxzzz_0_xxxxz_1, g_0_xxxzzz_0_xxxy_1, g_0_xxxzzz_0_xxxyy_0, g_0_xxxzzz_0_xxxyy_1, g_0_xxxzzz_0_xxxyz_0, g_0_xxxzzz_0_xxxyz_1, g_0_xxxzzz_0_xxxz_1, g_0_xxxzzz_0_xxxzz_0, g_0_xxxzzz_0_xxxzz_1, g_0_xxxzzz_0_xxyy_1, g_0_xxxzzz_0_xxyyy_0, g_0_xxxzzz_0_xxyyy_1, g_0_xxxzzz_0_xxyyz_0, g_0_xxxzzz_0_xxyyz_1, g_0_xxxzzz_0_xxyz_1, g_0_xxxzzz_0_xxyzz_0, g_0_xxxzzz_0_xxyzz_1, g_0_xxxzzz_0_xxzz_1, g_0_xxxzzz_0_xxzzz_0, g_0_xxxzzz_0_xxzzz_1, g_0_xxxzzz_0_xyyy_1, g_0_xxxzzz_0_xyyyy_0, g_0_xxxzzz_0_xyyyy_1, g_0_xxxzzz_0_xyyyz_0, g_0_xxxzzz_0_xyyyz_1, g_0_xxxzzz_0_xyyz_1, g_0_xxxzzz_0_xyyzz_0, g_0_xxxzzz_0_xyyzz_1, g_0_xxxzzz_0_xyzz_1, g_0_xxxzzz_0_xyzzz_0, g_0_xxxzzz_0_xyzzz_1, g_0_xxxzzz_0_xzzz_1, g_0_xxxzzz_0_xzzzz_0, g_0_xxxzzz_0_xzzzz_1, g_0_xxxzzz_0_yyyy_1, g_0_xxxzzz_0_yyyyy_0, g_0_xxxzzz_0_yyyyy_1, g_0_xxxzzz_0_yyyyz_0, g_0_xxxzzz_0_yyyyz_1, g_0_xxxzzz_0_yyyz_1, g_0_xxxzzz_0_yyyzz_0, g_0_xxxzzz_0_yyyzz_1, g_0_xxxzzz_0_yyzz_1, g_0_xxxzzz_0_yyzzz_0, g_0_xxxzzz_0_yyzzz_1, g_0_xxxzzz_0_yzzz_1, g_0_xxxzzz_0_yzzzz_0, g_0_xxxzzz_0_yzzzz_1, g_0_xxxzzz_0_zzzz_1, g_0_xxxzzz_0_zzzzz_0, g_0_xxxzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_xxxxx_0[i] = g_0_xxxzzz_0_xxxxx_0[i] * pb_y + g_0_xxxzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxy_0[i] = g_0_xxxzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxy_0[i] * pb_y + g_0_xxxzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxz_0[i] = g_0_xxxzzz_0_xxxxz_0[i] * pb_y + g_0_xxxzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyy_0[i] * pb_y + g_0_xxxzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyz_0[i] = g_0_xxxzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyz_0[i] * pb_y + g_0_xxxzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxzz_0[i] = g_0_xxxzzz_0_xxxzz_0[i] * pb_y + g_0_xxxzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyy_0[i] = 3.0 * g_0_xxxzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyy_0[i] * pb_y + g_0_xxxzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyz_0[i] = 2.0 * g_0_xxxzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyz_0[i] * pb_y + g_0_xxxzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyzz_0[i] = g_0_xxxzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyzz_0[i] * pb_y + g_0_xxxzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxzzz_0[i] = g_0_xxxzzz_0_xxzzz_0[i] * pb_y + g_0_xxxzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyy_0[i] * pb_y + g_0_xxxzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyz_0[i] = 3.0 * g_0_xxxzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyz_0[i] * pb_y + g_0_xxxzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyzz_0[i] = 2.0 * g_0_xxxzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyzz_0[i] * pb_y + g_0_xxxzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyzzz_0[i] = g_0_xxxzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzzz_0[i] * pb_y + g_0_xxxzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xzzzz_0[i] = g_0_xxxzzz_0_xzzzz_0[i] * pb_y + g_0_xxxzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyy_0[i] = 5.0 * g_0_xxxzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyy_0[i] * pb_y + g_0_xxxzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyz_0[i] = 4.0 * g_0_xxxzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyz_0[i] * pb_y + g_0_xxxzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyzz_0[i] = 3.0 * g_0_xxxzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyzz_0[i] * pb_y + g_0_xxxzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyzzz_0[i] = 2.0 * g_0_xxxzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyzzz_0[i] * pb_y + g_0_xxxzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yzzzz_0[i] = g_0_xxxzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yzzzz_0[i] * pb_y + g_0_xxxzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_zzzzz_0[i] = g_0_xxxzzz_0_zzzzz_0[i] * pb_y + g_0_xxxzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 294-315 components of targeted buffer : SKSH

    auto g_0_xxxzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 294);

    auto g_0_xxxzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 295);

    auto g_0_xxxzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 296);

    auto g_0_xxxzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 297);

    auto g_0_xxxzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 298);

    auto g_0_xxxzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 299);

    auto g_0_xxxzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 300);

    auto g_0_xxxzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 301);

    auto g_0_xxxzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 302);

    auto g_0_xxxzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 303);

    auto g_0_xxxzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 304);

    auto g_0_xxxzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 305);

    auto g_0_xxxzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 306);

    auto g_0_xxxzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 307);

    auto g_0_xxxzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 308);

    auto g_0_xxxzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 309);

    auto g_0_xxxzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 310);

    auto g_0_xxxzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 311);

    auto g_0_xxxzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 312);

    auto g_0_xxxzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 313);

    auto g_0_xxxzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 314);

    #pragma omp simd aligned(g_0_xxxzz_0_xxxxx_0, g_0_xxxzz_0_xxxxx_1, g_0_xxxzz_0_xxxxy_0, g_0_xxxzz_0_xxxxy_1, g_0_xxxzz_0_xxxyy_0, g_0_xxxzz_0_xxxyy_1, g_0_xxxzz_0_xxyyy_0, g_0_xxxzz_0_xxyyy_1, g_0_xxxzz_0_xyyyy_0, g_0_xxxzz_0_xyyyy_1, g_0_xxxzzz_0_xxxxx_0, g_0_xxxzzz_0_xxxxx_1, g_0_xxxzzz_0_xxxxy_0, g_0_xxxzzz_0_xxxxy_1, g_0_xxxzzz_0_xxxyy_0, g_0_xxxzzz_0_xxxyy_1, g_0_xxxzzz_0_xxyyy_0, g_0_xxxzzz_0_xxyyy_1, g_0_xxxzzz_0_xyyyy_0, g_0_xxxzzz_0_xyyyy_1, g_0_xxxzzzz_0_xxxxx_0, g_0_xxxzzzz_0_xxxxy_0, g_0_xxxzzzz_0_xxxxz_0, g_0_xxxzzzz_0_xxxyy_0, g_0_xxxzzzz_0_xxxyz_0, g_0_xxxzzzz_0_xxxzz_0, g_0_xxxzzzz_0_xxyyy_0, g_0_xxxzzzz_0_xxyyz_0, g_0_xxxzzzz_0_xxyzz_0, g_0_xxxzzzz_0_xxzzz_0, g_0_xxxzzzz_0_xyyyy_0, g_0_xxxzzzz_0_xyyyz_0, g_0_xxxzzzz_0_xyyzz_0, g_0_xxxzzzz_0_xyzzz_0, g_0_xxxzzzz_0_xzzzz_0, g_0_xxxzzzz_0_yyyyy_0, g_0_xxxzzzz_0_yyyyz_0, g_0_xxxzzzz_0_yyyzz_0, g_0_xxxzzzz_0_yyzzz_0, g_0_xxxzzzz_0_yzzzz_0, g_0_xxxzzzz_0_zzzzz_0, g_0_xxzzzz_0_xxxxz_0, g_0_xxzzzz_0_xxxxz_1, g_0_xxzzzz_0_xxxyz_0, g_0_xxzzzz_0_xxxyz_1, g_0_xxzzzz_0_xxxz_1, g_0_xxzzzz_0_xxxzz_0, g_0_xxzzzz_0_xxxzz_1, g_0_xxzzzz_0_xxyyz_0, g_0_xxzzzz_0_xxyyz_1, g_0_xxzzzz_0_xxyz_1, g_0_xxzzzz_0_xxyzz_0, g_0_xxzzzz_0_xxyzz_1, g_0_xxzzzz_0_xxzz_1, g_0_xxzzzz_0_xxzzz_0, g_0_xxzzzz_0_xxzzz_1, g_0_xxzzzz_0_xyyyz_0, g_0_xxzzzz_0_xyyyz_1, g_0_xxzzzz_0_xyyz_1, g_0_xxzzzz_0_xyyzz_0, g_0_xxzzzz_0_xyyzz_1, g_0_xxzzzz_0_xyzz_1, g_0_xxzzzz_0_xyzzz_0, g_0_xxzzzz_0_xyzzz_1, g_0_xxzzzz_0_xzzz_1, g_0_xxzzzz_0_xzzzz_0, g_0_xxzzzz_0_xzzzz_1, g_0_xxzzzz_0_yyyyy_0, g_0_xxzzzz_0_yyyyy_1, g_0_xxzzzz_0_yyyyz_0, g_0_xxzzzz_0_yyyyz_1, g_0_xxzzzz_0_yyyz_1, g_0_xxzzzz_0_yyyzz_0, g_0_xxzzzz_0_yyyzz_1, g_0_xxzzzz_0_yyzz_1, g_0_xxzzzz_0_yyzzz_0, g_0_xxzzzz_0_yyzzz_1, g_0_xxzzzz_0_yzzz_1, g_0_xxzzzz_0_yzzzz_0, g_0_xxzzzz_0_yzzzz_1, g_0_xxzzzz_0_zzzz_1, g_0_xxzzzz_0_zzzzz_0, g_0_xxzzzz_0_zzzzz_1, g_0_xzzzz_0_xxxxz_0, g_0_xzzzz_0_xxxxz_1, g_0_xzzzz_0_xxxyz_0, g_0_xzzzz_0_xxxyz_1, g_0_xzzzz_0_xxxzz_0, g_0_xzzzz_0_xxxzz_1, g_0_xzzzz_0_xxyyz_0, g_0_xzzzz_0_xxyyz_1, g_0_xzzzz_0_xxyzz_0, g_0_xzzzz_0_xxyzz_1, g_0_xzzzz_0_xxzzz_0, g_0_xzzzz_0_xxzzz_1, g_0_xzzzz_0_xyyyz_0, g_0_xzzzz_0_xyyyz_1, g_0_xzzzz_0_xyyzz_0, g_0_xzzzz_0_xyyzz_1, g_0_xzzzz_0_xyzzz_0, g_0_xzzzz_0_xyzzz_1, g_0_xzzzz_0_xzzzz_0, g_0_xzzzz_0_xzzzz_1, g_0_xzzzz_0_yyyyy_0, g_0_xzzzz_0_yyyyy_1, g_0_xzzzz_0_yyyyz_0, g_0_xzzzz_0_yyyyz_1, g_0_xzzzz_0_yyyzz_0, g_0_xzzzz_0_yyyzz_1, g_0_xzzzz_0_yyzzz_0, g_0_xzzzz_0_yyzzz_1, g_0_xzzzz_0_yzzzz_0, g_0_xzzzz_0_yzzzz_1, g_0_xzzzz_0_zzzzz_0, g_0_xzzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_xxxxx_0[i] = 3.0 * g_0_xxxzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxxx_0[i] * pb_z + g_0_xxxzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxy_0[i] = 3.0 * g_0_xxxzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxxy_0[i] * pb_z + g_0_xxxzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxz_0[i] = 2.0 * g_0_xzzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxz_0[i] * pb_x + g_0_xxzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxyy_0[i] = 3.0 * g_0_xxxzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxyy_0[i] * pb_z + g_0_xxxzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxyz_0[i] = 2.0 * g_0_xzzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyz_0[i] * pb_x + g_0_xxzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxzz_0[i] = 2.0 * g_0_xzzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxzz_0[i] * pb_x + g_0_xxzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxxzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxyyy_0[i] * pb_z + g_0_xxxzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxyyz_0[i] = 2.0 * g_0_xzzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyz_0[i] * pb_x + g_0_xxzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyzz_0[i] = 2.0 * g_0_xzzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyzz_0[i] * pb_x + g_0_xxzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxzzz_0[i] = 2.0 * g_0_xzzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxzzz_0[i] * pb_x + g_0_xxzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyyy_0[i] = 3.0 * g_0_xxxzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xyyyy_0[i] * pb_z + g_0_xxxzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xyyyz_0[i] = 2.0 * g_0_xzzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyz_0[i] * pb_x + g_0_xxzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyzz_0[i] = 2.0 * g_0_xzzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyzz_0[i] * pb_x + g_0_xxzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyzzz_0[i] = 2.0 * g_0_xzzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzzz_0[i] * pb_x + g_0_xxzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xzzzz_0[i] = 2.0 * g_0_xzzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xzzzz_0[i] * pb_x + g_0_xxzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyy_0[i] = 2.0 * g_0_xzzzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyy_0[i] * pb_x + g_0_xxzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyz_0[i] = 2.0 * g_0_xzzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyz_0[i] * pb_x + g_0_xxzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyzz_0[i] = 2.0 * g_0_xzzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyzz_0[i] * pb_x + g_0_xxzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyzzz_0[i] = 2.0 * g_0_xzzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyzzz_0[i] * pb_x + g_0_xxzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yzzzz_0[i] = 2.0 * g_0_xzzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzzzz_0[i] * pb_x + g_0_xxzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_zzzzz_0[i] = 2.0 * g_0_xzzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzzzz_0[i] * pb_x + g_0_xxzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-336 components of targeted buffer : SKSH

    auto g_0_xxyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 315);

    auto g_0_xxyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 316);

    auto g_0_xxyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 317);

    auto g_0_xxyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 318);

    auto g_0_xxyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 319);

    auto g_0_xxyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 320);

    auto g_0_xxyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 321);

    auto g_0_xxyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 322);

    auto g_0_xxyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 323);

    auto g_0_xxyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 324);

    auto g_0_xxyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 325);

    auto g_0_xxyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 326);

    auto g_0_xxyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 327);

    auto g_0_xxyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 328);

    auto g_0_xxyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 329);

    auto g_0_xxyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 330);

    auto g_0_xxyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 331);

    auto g_0_xxyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 332);

    auto g_0_xxyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 333);

    auto g_0_xxyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 334);

    auto g_0_xxyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 335);

    #pragma omp simd aligned(g_0_xxyyy_0_xxxxx_0, g_0_xxyyy_0_xxxxx_1, g_0_xxyyy_0_xxxxz_0, g_0_xxyyy_0_xxxxz_1, g_0_xxyyy_0_xxxzz_0, g_0_xxyyy_0_xxxzz_1, g_0_xxyyy_0_xxzzz_0, g_0_xxyyy_0_xxzzz_1, g_0_xxyyy_0_xzzzz_0, g_0_xxyyy_0_xzzzz_1, g_0_xxyyyy_0_xxxxx_0, g_0_xxyyyy_0_xxxxx_1, g_0_xxyyyy_0_xxxxz_0, g_0_xxyyyy_0_xxxxz_1, g_0_xxyyyy_0_xxxzz_0, g_0_xxyyyy_0_xxxzz_1, g_0_xxyyyy_0_xxzzz_0, g_0_xxyyyy_0_xxzzz_1, g_0_xxyyyy_0_xzzzz_0, g_0_xxyyyy_0_xzzzz_1, g_0_xxyyyyy_0_xxxxx_0, g_0_xxyyyyy_0_xxxxy_0, g_0_xxyyyyy_0_xxxxz_0, g_0_xxyyyyy_0_xxxyy_0, g_0_xxyyyyy_0_xxxyz_0, g_0_xxyyyyy_0_xxxzz_0, g_0_xxyyyyy_0_xxyyy_0, g_0_xxyyyyy_0_xxyyz_0, g_0_xxyyyyy_0_xxyzz_0, g_0_xxyyyyy_0_xxzzz_0, g_0_xxyyyyy_0_xyyyy_0, g_0_xxyyyyy_0_xyyyz_0, g_0_xxyyyyy_0_xyyzz_0, g_0_xxyyyyy_0_xyzzz_0, g_0_xxyyyyy_0_xzzzz_0, g_0_xxyyyyy_0_yyyyy_0, g_0_xxyyyyy_0_yyyyz_0, g_0_xxyyyyy_0_yyyzz_0, g_0_xxyyyyy_0_yyzzz_0, g_0_xxyyyyy_0_yzzzz_0, g_0_xxyyyyy_0_zzzzz_0, g_0_xyyyyy_0_xxxxy_0, g_0_xyyyyy_0_xxxxy_1, g_0_xyyyyy_0_xxxy_1, g_0_xyyyyy_0_xxxyy_0, g_0_xyyyyy_0_xxxyy_1, g_0_xyyyyy_0_xxxyz_0, g_0_xyyyyy_0_xxxyz_1, g_0_xyyyyy_0_xxyy_1, g_0_xyyyyy_0_xxyyy_0, g_0_xyyyyy_0_xxyyy_1, g_0_xyyyyy_0_xxyyz_0, g_0_xyyyyy_0_xxyyz_1, g_0_xyyyyy_0_xxyz_1, g_0_xyyyyy_0_xxyzz_0, g_0_xyyyyy_0_xxyzz_1, g_0_xyyyyy_0_xyyy_1, g_0_xyyyyy_0_xyyyy_0, g_0_xyyyyy_0_xyyyy_1, g_0_xyyyyy_0_xyyyz_0, g_0_xyyyyy_0_xyyyz_1, g_0_xyyyyy_0_xyyz_1, g_0_xyyyyy_0_xyyzz_0, g_0_xyyyyy_0_xyyzz_1, g_0_xyyyyy_0_xyzz_1, g_0_xyyyyy_0_xyzzz_0, g_0_xyyyyy_0_xyzzz_1, g_0_xyyyyy_0_yyyy_1, g_0_xyyyyy_0_yyyyy_0, g_0_xyyyyy_0_yyyyy_1, g_0_xyyyyy_0_yyyyz_0, g_0_xyyyyy_0_yyyyz_1, g_0_xyyyyy_0_yyyz_1, g_0_xyyyyy_0_yyyzz_0, g_0_xyyyyy_0_yyyzz_1, g_0_xyyyyy_0_yyzz_1, g_0_xyyyyy_0_yyzzz_0, g_0_xyyyyy_0_yyzzz_1, g_0_xyyyyy_0_yzzz_1, g_0_xyyyyy_0_yzzzz_0, g_0_xyyyyy_0_yzzzz_1, g_0_xyyyyy_0_zzzzz_0, g_0_xyyyyy_0_zzzzz_1, g_0_yyyyy_0_xxxxy_0, g_0_yyyyy_0_xxxxy_1, g_0_yyyyy_0_xxxyy_0, g_0_yyyyy_0_xxxyy_1, g_0_yyyyy_0_xxxyz_0, g_0_yyyyy_0_xxxyz_1, g_0_yyyyy_0_xxyyy_0, g_0_yyyyy_0_xxyyy_1, g_0_yyyyy_0_xxyyz_0, g_0_yyyyy_0_xxyyz_1, g_0_yyyyy_0_xxyzz_0, g_0_yyyyy_0_xxyzz_1, g_0_yyyyy_0_xyyyy_0, g_0_yyyyy_0_xyyyy_1, g_0_yyyyy_0_xyyyz_0, g_0_yyyyy_0_xyyyz_1, g_0_yyyyy_0_xyyzz_0, g_0_yyyyy_0_xyyzz_1, g_0_yyyyy_0_xyzzz_0, g_0_yyyyy_0_xyzzz_1, g_0_yyyyy_0_yyyyy_0, g_0_yyyyy_0_yyyyy_1, g_0_yyyyy_0_yyyyz_0, g_0_yyyyy_0_yyyyz_1, g_0_yyyyy_0_yyyzz_0, g_0_yyyyy_0_yyyzz_1, g_0_yyyyy_0_yyzzz_0, g_0_yyyyy_0_yyzzz_1, g_0_yyyyy_0_yzzzz_0, g_0_yyyyy_0_yzzzz_1, g_0_yyyyy_0_zzzzz_0, g_0_yyyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_xxxxx_0[i] = 4.0 * g_0_xxyyy_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxxx_0[i] * pb_y + g_0_xxyyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxxy_0[i] = g_0_yyyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxy_0[i] * pb_x + g_0_xyyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxz_0[i] = 4.0 * g_0_xxyyy_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxxz_0[i] * pb_y + g_0_xxyyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxyy_0[i] = g_0_yyyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyy_0[i] * pb_x + g_0_xyyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxyz_0[i] = g_0_yyyyy_0_xxxyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyz_0[i] * pb_x + g_0_xyyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxzz_0[i] = 4.0 * g_0_xxyyy_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxzz_0[i] * pb_y + g_0_xxyyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxyyy_0[i] = g_0_yyyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyy_0[i] * pb_x + g_0_xyyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyyz_0[i] = g_0_yyyyy_0_xxyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyz_0[i] * pb_x + g_0_xyyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyzz_0[i] = g_0_yyyyy_0_xxyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyzz_0[i] * pb_x + g_0_xyyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxzzz_0[i] = 4.0 * g_0_xxyyy_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxzzz_0[i] * pb_y + g_0_xxyyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xyyyy_0[i] = g_0_yyyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyyy_0[i] * pb_x + g_0_xyyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyyz_0[i] = g_0_yyyyy_0_xyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyyz_0[i] * pb_x + g_0_xyyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyzz_0[i] = g_0_yyyyy_0_xyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyzz_0[i] * pb_x + g_0_xyyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyzzz_0[i] = g_0_yyyyy_0_xyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyzzz_0[i] * pb_x + g_0_xyyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xzzzz_0[i] = 4.0 * g_0_xxyyy_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xzzzz_0[i] * pb_y + g_0_xxyyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_yyyyy_0[i] = g_0_yyyyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyy_0[i] * pb_x + g_0_xyyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyyz_0[i] = g_0_yyyyy_0_yyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyz_0[i] * pb_x + g_0_xyyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyzz_0[i] = g_0_yyyyy_0_yyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyzz_0[i] * pb_x + g_0_xyyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyzzz_0[i] = g_0_yyyyy_0_yyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzzz_0[i] * pb_x + g_0_xyyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yzzzz_0[i] = g_0_yyyyy_0_yzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzzz_0[i] * pb_x + g_0_xyyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_zzzzz_0[i] = g_0_yyyyy_0_zzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_zzzzz_0[i] * pb_x + g_0_xyyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 336-357 components of targeted buffer : SKSH

    auto g_0_xxyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 336);

    auto g_0_xxyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 337);

    auto g_0_xxyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 338);

    auto g_0_xxyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 339);

    auto g_0_xxyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 340);

    auto g_0_xxyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 341);

    auto g_0_xxyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 342);

    auto g_0_xxyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 343);

    auto g_0_xxyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 344);

    auto g_0_xxyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 345);

    auto g_0_xxyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 346);

    auto g_0_xxyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 347);

    auto g_0_xxyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 348);

    auto g_0_xxyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 349);

    auto g_0_xxyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 350);

    auto g_0_xxyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 351);

    auto g_0_xxyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 352);

    auto g_0_xxyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 353);

    auto g_0_xxyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 354);

    auto g_0_xxyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 355);

    auto g_0_xxyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 356);

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxx_1, g_0_xxyyyy_0_xxxxx_0, g_0_xxyyyy_0_xxxxx_1, g_0_xxyyyy_0_xxxxy_0, g_0_xxyyyy_0_xxxxy_1, g_0_xxyyyy_0_xxxxz_0, g_0_xxyyyy_0_xxxxz_1, g_0_xxyyyy_0_xxxy_1, g_0_xxyyyy_0_xxxyy_0, g_0_xxyyyy_0_xxxyy_1, g_0_xxyyyy_0_xxxyz_0, g_0_xxyyyy_0_xxxyz_1, g_0_xxyyyy_0_xxxz_1, g_0_xxyyyy_0_xxxzz_0, g_0_xxyyyy_0_xxxzz_1, g_0_xxyyyy_0_xxyy_1, g_0_xxyyyy_0_xxyyy_0, g_0_xxyyyy_0_xxyyy_1, g_0_xxyyyy_0_xxyyz_0, g_0_xxyyyy_0_xxyyz_1, g_0_xxyyyy_0_xxyz_1, g_0_xxyyyy_0_xxyzz_0, g_0_xxyyyy_0_xxyzz_1, g_0_xxyyyy_0_xxzz_1, g_0_xxyyyy_0_xxzzz_0, g_0_xxyyyy_0_xxzzz_1, g_0_xxyyyy_0_xyyy_1, g_0_xxyyyy_0_xyyyy_0, g_0_xxyyyy_0_xyyyy_1, g_0_xxyyyy_0_xyyyz_0, g_0_xxyyyy_0_xyyyz_1, g_0_xxyyyy_0_xyyz_1, g_0_xxyyyy_0_xyyzz_0, g_0_xxyyyy_0_xyyzz_1, g_0_xxyyyy_0_xyzz_1, g_0_xxyyyy_0_xyzzz_0, g_0_xxyyyy_0_xyzzz_1, g_0_xxyyyy_0_xzzz_1, g_0_xxyyyy_0_xzzzz_0, g_0_xxyyyy_0_xzzzz_1, g_0_xxyyyy_0_yyyy_1, g_0_xxyyyy_0_yyyyy_0, g_0_xxyyyy_0_yyyyy_1, g_0_xxyyyy_0_yyyyz_0, g_0_xxyyyy_0_yyyyz_1, g_0_xxyyyy_0_yyyz_1, g_0_xxyyyy_0_yyyzz_0, g_0_xxyyyy_0_yyyzz_1, g_0_xxyyyy_0_yyzz_1, g_0_xxyyyy_0_yyzzz_0, g_0_xxyyyy_0_yyzzz_1, g_0_xxyyyy_0_yzzz_1, g_0_xxyyyy_0_yzzzz_0, g_0_xxyyyy_0_yzzzz_1, g_0_xxyyyy_0_zzzz_1, g_0_xxyyyy_0_zzzzz_0, g_0_xxyyyy_0_zzzzz_1, g_0_xxyyyyz_0_xxxxx_0, g_0_xxyyyyz_0_xxxxy_0, g_0_xxyyyyz_0_xxxxz_0, g_0_xxyyyyz_0_xxxyy_0, g_0_xxyyyyz_0_xxxyz_0, g_0_xxyyyyz_0_xxxzz_0, g_0_xxyyyyz_0_xxyyy_0, g_0_xxyyyyz_0_xxyyz_0, g_0_xxyyyyz_0_xxyzz_0, g_0_xxyyyyz_0_xxzzz_0, g_0_xxyyyyz_0_xyyyy_0, g_0_xxyyyyz_0_xyyyz_0, g_0_xxyyyyz_0_xyyzz_0, g_0_xxyyyyz_0_xyzzz_0, g_0_xxyyyyz_0_xzzzz_0, g_0_xxyyyyz_0_yyyyy_0, g_0_xxyyyyz_0_yyyyz_0, g_0_xxyyyyz_0_yyyzz_0, g_0_xxyyyyz_0_yyzzz_0, g_0_xxyyyyz_0_yzzzz_0, g_0_xxyyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_xxxxx_0[i] = g_0_xxyyyy_0_xxxxx_0[i] * pb_z + g_0_xxyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxy_0[i] = g_0_xxyyyy_0_xxxxy_0[i] * pb_z + g_0_xxyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxz_0[i] = g_0_xxyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxz_0[i] * pb_z + g_0_xxyyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyy_0[i] = g_0_xxyyyy_0_xxxyy_0[i] * pb_z + g_0_xxyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyz_0[i] = g_0_xxyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyz_0[i] * pb_z + g_0_xxyyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxzz_0[i] = 2.0 * g_0_xxyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxzz_0[i] * pb_z + g_0_xxyyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyy_0[i] = g_0_xxyyyy_0_xxyyy_0[i] * pb_z + g_0_xxyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyz_0[i] = g_0_xxyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyz_0[i] * pb_z + g_0_xxyyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyzz_0[i] = 2.0 * g_0_xxyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyzz_0[i] * pb_z + g_0_xxyyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxzzz_0[i] * pb_z + g_0_xxyyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyy_0[i] = g_0_xxyyyy_0_xyyyy_0[i] * pb_z + g_0_xxyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyz_0[i] = g_0_xxyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyz_0[i] * pb_z + g_0_xxyyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyzz_0[i] = 2.0 * g_0_xxyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyzz_0[i] * pb_z + g_0_xxyyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzzz_0[i] * pb_z + g_0_xxyyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xzzzz_0[i] = 4.0 * g_0_xxyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xzzzz_0[i] * pb_z + g_0_xxyyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyy_0[i] = g_0_xxyyyy_0_yyyyy_0[i] * pb_z + g_0_xxyyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyz_0[i] = g_0_xxyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyyz_0[i] * pb_z + g_0_xxyyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyzz_0[i] = 2.0 * g_0_xxyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyzz_0[i] * pb_z + g_0_xxyyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyzzz_0[i] * pb_z + g_0_xxyyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yzzzz_0[i] = 4.0 * g_0_xxyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yzzzz_0[i] * pb_z + g_0_xxyyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_zzzzz_0[i] = 5.0 * g_0_xxyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_zzzzz_0[i] * pb_z + g_0_xxyyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 357-378 components of targeted buffer : SKSH

    auto g_0_xxyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 357);

    auto g_0_xxyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 358);

    auto g_0_xxyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 359);

    auto g_0_xxyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 360);

    auto g_0_xxyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 361);

    auto g_0_xxyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 362);

    auto g_0_xxyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 363);

    auto g_0_xxyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 364);

    auto g_0_xxyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 365);

    auto g_0_xxyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 366);

    auto g_0_xxyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 367);

    auto g_0_xxyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 368);

    auto g_0_xxyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 369);

    auto g_0_xxyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 370);

    auto g_0_xxyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 371);

    auto g_0_xxyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 372);

    auto g_0_xxyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 373);

    auto g_0_xxyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 374);

    auto g_0_xxyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 375);

    auto g_0_xxyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 376);

    auto g_0_xxyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 377);

    #pragma omp simd aligned(g_0_xxyyy_0_xxxxy_0, g_0_xxyyy_0_xxxxy_1, g_0_xxyyy_0_xxxyy_0, g_0_xxyyy_0_xxxyy_1, g_0_xxyyy_0_xxyyy_0, g_0_xxyyy_0_xxyyy_1, g_0_xxyyy_0_xyyyy_0, g_0_xxyyy_0_xyyyy_1, g_0_xxyyyz_0_xxxxy_0, g_0_xxyyyz_0_xxxxy_1, g_0_xxyyyz_0_xxxyy_0, g_0_xxyyyz_0_xxxyy_1, g_0_xxyyyz_0_xxyyy_0, g_0_xxyyyz_0_xxyyy_1, g_0_xxyyyz_0_xyyyy_0, g_0_xxyyyz_0_xyyyy_1, g_0_xxyyyzz_0_xxxxx_0, g_0_xxyyyzz_0_xxxxy_0, g_0_xxyyyzz_0_xxxxz_0, g_0_xxyyyzz_0_xxxyy_0, g_0_xxyyyzz_0_xxxyz_0, g_0_xxyyyzz_0_xxxzz_0, g_0_xxyyyzz_0_xxyyy_0, g_0_xxyyyzz_0_xxyyz_0, g_0_xxyyyzz_0_xxyzz_0, g_0_xxyyyzz_0_xxzzz_0, g_0_xxyyyzz_0_xyyyy_0, g_0_xxyyyzz_0_xyyyz_0, g_0_xxyyyzz_0_xyyzz_0, g_0_xxyyyzz_0_xyzzz_0, g_0_xxyyyzz_0_xzzzz_0, g_0_xxyyyzz_0_yyyyy_0, g_0_xxyyyzz_0_yyyyz_0, g_0_xxyyyzz_0_yyyzz_0, g_0_xxyyyzz_0_yyzzz_0, g_0_xxyyyzz_0_yzzzz_0, g_0_xxyyyzz_0_zzzzz_0, g_0_xxyyzz_0_xxxxx_0, g_0_xxyyzz_0_xxxxx_1, g_0_xxyyzz_0_xxxxz_0, g_0_xxyyzz_0_xxxxz_1, g_0_xxyyzz_0_xxxzz_0, g_0_xxyyzz_0_xxxzz_1, g_0_xxyyzz_0_xxzzz_0, g_0_xxyyzz_0_xxzzz_1, g_0_xxyyzz_0_xzzzz_0, g_0_xxyyzz_0_xzzzz_1, g_0_xxyzz_0_xxxxx_0, g_0_xxyzz_0_xxxxx_1, g_0_xxyzz_0_xxxxz_0, g_0_xxyzz_0_xxxxz_1, g_0_xxyzz_0_xxxzz_0, g_0_xxyzz_0_xxxzz_1, g_0_xxyzz_0_xxzzz_0, g_0_xxyzz_0_xxzzz_1, g_0_xxyzz_0_xzzzz_0, g_0_xxyzz_0_xzzzz_1, g_0_xyyyzz_0_xxxyz_0, g_0_xyyyzz_0_xxxyz_1, g_0_xyyyzz_0_xxyyz_0, g_0_xyyyzz_0_xxyyz_1, g_0_xyyyzz_0_xxyz_1, g_0_xyyyzz_0_xxyzz_0, g_0_xyyyzz_0_xxyzz_1, g_0_xyyyzz_0_xyyyz_0, g_0_xyyyzz_0_xyyyz_1, g_0_xyyyzz_0_xyyz_1, g_0_xyyyzz_0_xyyzz_0, g_0_xyyyzz_0_xyyzz_1, g_0_xyyyzz_0_xyzz_1, g_0_xyyyzz_0_xyzzz_0, g_0_xyyyzz_0_xyzzz_1, g_0_xyyyzz_0_yyyyy_0, g_0_xyyyzz_0_yyyyy_1, g_0_xyyyzz_0_yyyyz_0, g_0_xyyyzz_0_yyyyz_1, g_0_xyyyzz_0_yyyz_1, g_0_xyyyzz_0_yyyzz_0, g_0_xyyyzz_0_yyyzz_1, g_0_xyyyzz_0_yyzz_1, g_0_xyyyzz_0_yyzzz_0, g_0_xyyyzz_0_yyzzz_1, g_0_xyyyzz_0_yzzz_1, g_0_xyyyzz_0_yzzzz_0, g_0_xyyyzz_0_yzzzz_1, g_0_xyyyzz_0_zzzzz_0, g_0_xyyyzz_0_zzzzz_1, g_0_yyyzz_0_xxxyz_0, g_0_yyyzz_0_xxxyz_1, g_0_yyyzz_0_xxyyz_0, g_0_yyyzz_0_xxyyz_1, g_0_yyyzz_0_xxyzz_0, g_0_yyyzz_0_xxyzz_1, g_0_yyyzz_0_xyyyz_0, g_0_yyyzz_0_xyyyz_1, g_0_yyyzz_0_xyyzz_0, g_0_yyyzz_0_xyyzz_1, g_0_yyyzz_0_xyzzz_0, g_0_yyyzz_0_xyzzz_1, g_0_yyyzz_0_yyyyy_0, g_0_yyyzz_0_yyyyy_1, g_0_yyyzz_0_yyyyz_0, g_0_yyyzz_0_yyyyz_1, g_0_yyyzz_0_yyyzz_0, g_0_yyyzz_0_yyyzz_1, g_0_yyyzz_0_yyzzz_0, g_0_yyyzz_0_yyzzz_1, g_0_yyyzz_0_yzzzz_0, g_0_yyyzz_0_yzzzz_1, g_0_yyyzz_0_zzzzz_0, g_0_yyyzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzz_0_xxxxx_0[i] = 2.0 * g_0_xxyzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxx_0[i] * pb_y + g_0_xxyyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxxy_0[i] = g_0_xxyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxxy_0[i] * pb_z + g_0_xxyyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxxz_0[i] = 2.0 * g_0_xxyzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxz_0[i] * pb_y + g_0_xxyyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxyy_0[i] = g_0_xxyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxyy_0[i] * pb_z + g_0_xxyyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxyz_0[i] = g_0_yyyzz_0_xxxyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxyz_0[i] * pb_x + g_0_xyyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxzz_0[i] = 2.0 * g_0_xxyzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxzz_0[i] * pb_y + g_0_xxyyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxyyy_0[i] = g_0_xxyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxyyy_0[i] * pb_z + g_0_xxyyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxyyz_0[i] = g_0_yyyzz_0_xxyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyyz_0[i] * pb_x + g_0_xyyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxyzz_0[i] = g_0_yyyzz_0_xxyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyzz_0[i] * pb_x + g_0_xyyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxzzz_0[i] = 2.0 * g_0_xxyzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxzzz_0[i] * pb_y + g_0_xxyyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xyyyy_0[i] = g_0_xxyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xyyyy_0[i] * pb_z + g_0_xxyyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xyyyz_0[i] = g_0_yyyzz_0_xyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyyyz_0[i] * pb_x + g_0_xyyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyyzz_0[i] = g_0_yyyzz_0_xyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyyzz_0[i] * pb_x + g_0_xyyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyzzz_0[i] = g_0_yyyzz_0_xyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyzzz_0[i] * pb_x + g_0_xyyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xzzzz_0[i] = 2.0 * g_0_xxyzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xzzzz_0[i] * pb_y + g_0_xxyyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_yyyyy_0[i] = g_0_yyyzz_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyy_0[i] * pb_x + g_0_xyyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyyz_0[i] = g_0_yyyzz_0_yyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyz_0[i] * pb_x + g_0_xyyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyzz_0[i] = g_0_yyyzz_0_yyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyzz_0[i] * pb_x + g_0_xyyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyzzz_0[i] = g_0_yyyzz_0_yyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzzz_0[i] * pb_x + g_0_xyyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yzzzz_0[i] = g_0_yyyzz_0_yzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzzz_0[i] * pb_x + g_0_xyyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_zzzzz_0[i] = g_0_yyyzz_0_zzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_zzzzz_0[i] * pb_x + g_0_xyyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 378-399 components of targeted buffer : SKSH

    auto g_0_xxyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 378);

    auto g_0_xxyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 379);

    auto g_0_xxyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 380);

    auto g_0_xxyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 381);

    auto g_0_xxyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 382);

    auto g_0_xxyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 383);

    auto g_0_xxyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 384);

    auto g_0_xxyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 385);

    auto g_0_xxyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 386);

    auto g_0_xxyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 387);

    auto g_0_xxyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 388);

    auto g_0_xxyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 389);

    auto g_0_xxyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 390);

    auto g_0_xxyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 391);

    auto g_0_xxyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 392);

    auto g_0_xxyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 393);

    auto g_0_xxyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 394);

    auto g_0_xxyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 395);

    auto g_0_xxyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 396);

    auto g_0_xxyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 397);

    auto g_0_xxyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 398);

    #pragma omp simd aligned(g_0_xxyyz_0_xxxxy_0, g_0_xxyyz_0_xxxxy_1, g_0_xxyyz_0_xxxyy_0, g_0_xxyyz_0_xxxyy_1, g_0_xxyyz_0_xxyyy_0, g_0_xxyyz_0_xxyyy_1, g_0_xxyyz_0_xyyyy_0, g_0_xxyyz_0_xyyyy_1, g_0_xxyyzz_0_xxxxy_0, g_0_xxyyzz_0_xxxxy_1, g_0_xxyyzz_0_xxxyy_0, g_0_xxyyzz_0_xxxyy_1, g_0_xxyyzz_0_xxyyy_0, g_0_xxyyzz_0_xxyyy_1, g_0_xxyyzz_0_xyyyy_0, g_0_xxyyzz_0_xyyyy_1, g_0_xxyyzzz_0_xxxxx_0, g_0_xxyyzzz_0_xxxxy_0, g_0_xxyyzzz_0_xxxxz_0, g_0_xxyyzzz_0_xxxyy_0, g_0_xxyyzzz_0_xxxyz_0, g_0_xxyyzzz_0_xxxzz_0, g_0_xxyyzzz_0_xxyyy_0, g_0_xxyyzzz_0_xxyyz_0, g_0_xxyyzzz_0_xxyzz_0, g_0_xxyyzzz_0_xxzzz_0, g_0_xxyyzzz_0_xyyyy_0, g_0_xxyyzzz_0_xyyyz_0, g_0_xxyyzzz_0_xyyzz_0, g_0_xxyyzzz_0_xyzzz_0, g_0_xxyyzzz_0_xzzzz_0, g_0_xxyyzzz_0_yyyyy_0, g_0_xxyyzzz_0_yyyyz_0, g_0_xxyyzzz_0_yyyzz_0, g_0_xxyyzzz_0_yyzzz_0, g_0_xxyyzzz_0_yzzzz_0, g_0_xxyyzzz_0_zzzzz_0, g_0_xxyzzz_0_xxxxx_0, g_0_xxyzzz_0_xxxxx_1, g_0_xxyzzz_0_xxxxz_0, g_0_xxyzzz_0_xxxxz_1, g_0_xxyzzz_0_xxxzz_0, g_0_xxyzzz_0_xxxzz_1, g_0_xxyzzz_0_xxzzz_0, g_0_xxyzzz_0_xxzzz_1, g_0_xxyzzz_0_xzzzz_0, g_0_xxyzzz_0_xzzzz_1, g_0_xxzzz_0_xxxxx_0, g_0_xxzzz_0_xxxxx_1, g_0_xxzzz_0_xxxxz_0, g_0_xxzzz_0_xxxxz_1, g_0_xxzzz_0_xxxzz_0, g_0_xxzzz_0_xxxzz_1, g_0_xxzzz_0_xxzzz_0, g_0_xxzzz_0_xxzzz_1, g_0_xxzzz_0_xzzzz_0, g_0_xxzzz_0_xzzzz_1, g_0_xyyzzz_0_xxxyz_0, g_0_xyyzzz_0_xxxyz_1, g_0_xyyzzz_0_xxyyz_0, g_0_xyyzzz_0_xxyyz_1, g_0_xyyzzz_0_xxyz_1, g_0_xyyzzz_0_xxyzz_0, g_0_xyyzzz_0_xxyzz_1, g_0_xyyzzz_0_xyyyz_0, g_0_xyyzzz_0_xyyyz_1, g_0_xyyzzz_0_xyyz_1, g_0_xyyzzz_0_xyyzz_0, g_0_xyyzzz_0_xyyzz_1, g_0_xyyzzz_0_xyzz_1, g_0_xyyzzz_0_xyzzz_0, g_0_xyyzzz_0_xyzzz_1, g_0_xyyzzz_0_yyyyy_0, g_0_xyyzzz_0_yyyyy_1, g_0_xyyzzz_0_yyyyz_0, g_0_xyyzzz_0_yyyyz_1, g_0_xyyzzz_0_yyyz_1, g_0_xyyzzz_0_yyyzz_0, g_0_xyyzzz_0_yyyzz_1, g_0_xyyzzz_0_yyzz_1, g_0_xyyzzz_0_yyzzz_0, g_0_xyyzzz_0_yyzzz_1, g_0_xyyzzz_0_yzzz_1, g_0_xyyzzz_0_yzzzz_0, g_0_xyyzzz_0_yzzzz_1, g_0_xyyzzz_0_zzzzz_0, g_0_xyyzzz_0_zzzzz_1, g_0_yyzzz_0_xxxyz_0, g_0_yyzzz_0_xxxyz_1, g_0_yyzzz_0_xxyyz_0, g_0_yyzzz_0_xxyyz_1, g_0_yyzzz_0_xxyzz_0, g_0_yyzzz_0_xxyzz_1, g_0_yyzzz_0_xyyyz_0, g_0_yyzzz_0_xyyyz_1, g_0_yyzzz_0_xyyzz_0, g_0_yyzzz_0_xyyzz_1, g_0_yyzzz_0_xyzzz_0, g_0_yyzzz_0_xyzzz_1, g_0_yyzzz_0_yyyyy_0, g_0_yyzzz_0_yyyyy_1, g_0_yyzzz_0_yyyyz_0, g_0_yyzzz_0_yyyyz_1, g_0_yyzzz_0_yyyzz_0, g_0_yyzzz_0_yyyzz_1, g_0_yyzzz_0_yyzzz_0, g_0_yyzzz_0_yyzzz_1, g_0_yyzzz_0_yzzzz_0, g_0_yyzzz_0_yzzzz_1, g_0_yyzzz_0_zzzzz_0, g_0_yyzzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzz_0_xxxxx_0[i] = g_0_xxzzz_0_xxxxx_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxx_0[i] * pb_y + g_0_xxyzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxxy_0[i] = 2.0 * g_0_xxyyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxy_0[i] * pb_z + g_0_xxyyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxxz_0[i] = g_0_xxzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxz_0[i] * pb_y + g_0_xxyzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxyy_0[i] = 2.0 * g_0_xxyyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxyy_0[i] * pb_z + g_0_xxyyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxyz_0[i] = g_0_yyzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxyz_0[i] * pb_x + g_0_xyyzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxzz_0[i] = g_0_xxzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxzz_0[i] * pb_y + g_0_xxyzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxyyy_0[i] = 2.0 * g_0_xxyyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxyyy_0[i] * pb_z + g_0_xxyyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxyyz_0[i] = g_0_yyzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyyz_0[i] * pb_x + g_0_xyyzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxyzz_0[i] = g_0_yyzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyzz_0[i] * pb_x + g_0_xyyzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxzzz_0[i] = g_0_xxzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxzzz_0[i] * pb_y + g_0_xxyzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xyyyy_0[i] = 2.0 * g_0_xxyyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xyyyy_0[i] * pb_z + g_0_xxyyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xyyyz_0[i] = g_0_yyzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyyyz_0[i] * pb_x + g_0_xyyzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyyzz_0[i] = g_0_yyzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyyzz_0[i] * pb_x + g_0_xyyzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyzzz_0[i] = g_0_yyzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyzzz_0[i] * pb_x + g_0_xyyzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xzzzz_0[i] = g_0_xxzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xzzzz_0[i] * pb_y + g_0_xxyzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_yyyyy_0[i] = g_0_yyzzz_0_yyyyy_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyy_0[i] * pb_x + g_0_xyyzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyyz_0[i] = g_0_yyzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyz_0[i] * pb_x + g_0_xyyzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyzz_0[i] = g_0_yyzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyzz_0[i] * pb_x + g_0_xyyzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyzzz_0[i] = g_0_yyzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzzz_0[i] * pb_x + g_0_xyyzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yzzzz_0[i] = g_0_yyzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzzz_0[i] * pb_x + g_0_xyyzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_zzzzz_0[i] = g_0_yyzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_zzzzz_0[i] * pb_x + g_0_xyyzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 399-420 components of targeted buffer : SKSH

    auto g_0_xxyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 399);

    auto g_0_xxyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 400);

    auto g_0_xxyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 401);

    auto g_0_xxyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 402);

    auto g_0_xxyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 403);

    auto g_0_xxyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 404);

    auto g_0_xxyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 405);

    auto g_0_xxyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 406);

    auto g_0_xxyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 407);

    auto g_0_xxyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 408);

    auto g_0_xxyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 409);

    auto g_0_xxyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 410);

    auto g_0_xxyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 411);

    auto g_0_xxyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 412);

    auto g_0_xxyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 413);

    auto g_0_xxyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 414);

    auto g_0_xxyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 415);

    auto g_0_xxyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 416);

    auto g_0_xxyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 417);

    auto g_0_xxyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 418);

    auto g_0_xxyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 419);

    #pragma omp simd aligned(g_0_xxyzzzz_0_xxxxx_0, g_0_xxyzzzz_0_xxxxy_0, g_0_xxyzzzz_0_xxxxz_0, g_0_xxyzzzz_0_xxxyy_0, g_0_xxyzzzz_0_xxxyz_0, g_0_xxyzzzz_0_xxxzz_0, g_0_xxyzzzz_0_xxyyy_0, g_0_xxyzzzz_0_xxyyz_0, g_0_xxyzzzz_0_xxyzz_0, g_0_xxyzzzz_0_xxzzz_0, g_0_xxyzzzz_0_xyyyy_0, g_0_xxyzzzz_0_xyyyz_0, g_0_xxyzzzz_0_xyyzz_0, g_0_xxyzzzz_0_xyzzz_0, g_0_xxyzzzz_0_xzzzz_0, g_0_xxyzzzz_0_yyyyy_0, g_0_xxyzzzz_0_yyyyz_0, g_0_xxyzzzz_0_yyyzz_0, g_0_xxyzzzz_0_yyzzz_0, g_0_xxyzzzz_0_yzzzz_0, g_0_xxyzzzz_0_zzzzz_0, g_0_xxzzzz_0_xxxx_1, g_0_xxzzzz_0_xxxxx_0, g_0_xxzzzz_0_xxxxx_1, g_0_xxzzzz_0_xxxxy_0, g_0_xxzzzz_0_xxxxy_1, g_0_xxzzzz_0_xxxxz_0, g_0_xxzzzz_0_xxxxz_1, g_0_xxzzzz_0_xxxy_1, g_0_xxzzzz_0_xxxyy_0, g_0_xxzzzz_0_xxxyy_1, g_0_xxzzzz_0_xxxyz_0, g_0_xxzzzz_0_xxxyz_1, g_0_xxzzzz_0_xxxz_1, g_0_xxzzzz_0_xxxzz_0, g_0_xxzzzz_0_xxxzz_1, g_0_xxzzzz_0_xxyy_1, g_0_xxzzzz_0_xxyyy_0, g_0_xxzzzz_0_xxyyy_1, g_0_xxzzzz_0_xxyyz_0, g_0_xxzzzz_0_xxyyz_1, g_0_xxzzzz_0_xxyz_1, g_0_xxzzzz_0_xxyzz_0, g_0_xxzzzz_0_xxyzz_1, g_0_xxzzzz_0_xxzz_1, g_0_xxzzzz_0_xxzzz_0, g_0_xxzzzz_0_xxzzz_1, g_0_xxzzzz_0_xyyy_1, g_0_xxzzzz_0_xyyyy_0, g_0_xxzzzz_0_xyyyy_1, g_0_xxzzzz_0_xyyyz_0, g_0_xxzzzz_0_xyyyz_1, g_0_xxzzzz_0_xyyz_1, g_0_xxzzzz_0_xyyzz_0, g_0_xxzzzz_0_xyyzz_1, g_0_xxzzzz_0_xyzz_1, g_0_xxzzzz_0_xyzzz_0, g_0_xxzzzz_0_xyzzz_1, g_0_xxzzzz_0_xzzz_1, g_0_xxzzzz_0_xzzzz_0, g_0_xxzzzz_0_xzzzz_1, g_0_xxzzzz_0_yyyy_1, g_0_xxzzzz_0_yyyyy_0, g_0_xxzzzz_0_yyyyy_1, g_0_xxzzzz_0_yyyyz_0, g_0_xxzzzz_0_yyyyz_1, g_0_xxzzzz_0_yyyz_1, g_0_xxzzzz_0_yyyzz_0, g_0_xxzzzz_0_yyyzz_1, g_0_xxzzzz_0_yyzz_1, g_0_xxzzzz_0_yyzzz_0, g_0_xxzzzz_0_yyzzz_1, g_0_xxzzzz_0_yzzz_1, g_0_xxzzzz_0_yzzzz_0, g_0_xxzzzz_0_yzzzz_1, g_0_xxzzzz_0_zzzz_1, g_0_xxzzzz_0_zzzzz_0, g_0_xxzzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_xxxxx_0[i] = g_0_xxzzzz_0_xxxxx_0[i] * pb_y + g_0_xxzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxy_0[i] = g_0_xxzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxy_0[i] * pb_y + g_0_xxzzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxz_0[i] = g_0_xxzzzz_0_xxxxz_0[i] * pb_y + g_0_xxzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyy_0[i] = 2.0 * g_0_xxzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyy_0[i] * pb_y + g_0_xxzzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyz_0[i] = g_0_xxzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyz_0[i] * pb_y + g_0_xxzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxzz_0[i] = g_0_xxzzzz_0_xxxzz_0[i] * pb_y + g_0_xxzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyy_0[i] * pb_y + g_0_xxzzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyz_0[i] = 2.0 * g_0_xxzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyz_0[i] * pb_y + g_0_xxzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyzz_0[i] = g_0_xxzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyzz_0[i] * pb_y + g_0_xxzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxzzz_0[i] = g_0_xxzzzz_0_xxzzz_0[i] * pb_y + g_0_xxzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyy_0[i] = 4.0 * g_0_xxzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyy_0[i] * pb_y + g_0_xxzzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyz_0[i] * pb_y + g_0_xxzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyzz_0[i] = 2.0 * g_0_xxzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyzz_0[i] * pb_y + g_0_xxzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyzzz_0[i] = g_0_xxzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzzz_0[i] * pb_y + g_0_xxzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xzzzz_0[i] = g_0_xxzzzz_0_xzzzz_0[i] * pb_y + g_0_xxzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyy_0[i] = 5.0 * g_0_xxzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyy_0[i] * pb_y + g_0_xxzzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyz_0[i] = 4.0 * g_0_xxzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyz_0[i] * pb_y + g_0_xxzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyzz_0[i] * pb_y + g_0_xxzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyzzz_0[i] = 2.0 * g_0_xxzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyzzz_0[i] * pb_y + g_0_xxzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yzzzz_0[i] = g_0_xxzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yzzzz_0[i] * pb_y + g_0_xxzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_zzzzz_0[i] = g_0_xxzzzz_0_zzzzz_0[i] * pb_y + g_0_xxzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 420-441 components of targeted buffer : SKSH

    auto g_0_xxzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 420);

    auto g_0_xxzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 421);

    auto g_0_xxzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 422);

    auto g_0_xxzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 423);

    auto g_0_xxzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 424);

    auto g_0_xxzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 425);

    auto g_0_xxzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 426);

    auto g_0_xxzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 427);

    auto g_0_xxzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 428);

    auto g_0_xxzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 429);

    auto g_0_xxzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 430);

    auto g_0_xxzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 431);

    auto g_0_xxzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 432);

    auto g_0_xxzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 433);

    auto g_0_xxzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 434);

    auto g_0_xxzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 435);

    auto g_0_xxzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 436);

    auto g_0_xxzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 437);

    auto g_0_xxzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 438);

    auto g_0_xxzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 439);

    auto g_0_xxzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 440);

    #pragma omp simd aligned(g_0_xxzzz_0_xxxxx_0, g_0_xxzzz_0_xxxxx_1, g_0_xxzzz_0_xxxxy_0, g_0_xxzzz_0_xxxxy_1, g_0_xxzzz_0_xxxyy_0, g_0_xxzzz_0_xxxyy_1, g_0_xxzzz_0_xxyyy_0, g_0_xxzzz_0_xxyyy_1, g_0_xxzzz_0_xyyyy_0, g_0_xxzzz_0_xyyyy_1, g_0_xxzzzz_0_xxxxx_0, g_0_xxzzzz_0_xxxxx_1, g_0_xxzzzz_0_xxxxy_0, g_0_xxzzzz_0_xxxxy_1, g_0_xxzzzz_0_xxxyy_0, g_0_xxzzzz_0_xxxyy_1, g_0_xxzzzz_0_xxyyy_0, g_0_xxzzzz_0_xxyyy_1, g_0_xxzzzz_0_xyyyy_0, g_0_xxzzzz_0_xyyyy_1, g_0_xxzzzzz_0_xxxxx_0, g_0_xxzzzzz_0_xxxxy_0, g_0_xxzzzzz_0_xxxxz_0, g_0_xxzzzzz_0_xxxyy_0, g_0_xxzzzzz_0_xxxyz_0, g_0_xxzzzzz_0_xxxzz_0, g_0_xxzzzzz_0_xxyyy_0, g_0_xxzzzzz_0_xxyyz_0, g_0_xxzzzzz_0_xxyzz_0, g_0_xxzzzzz_0_xxzzz_0, g_0_xxzzzzz_0_xyyyy_0, g_0_xxzzzzz_0_xyyyz_0, g_0_xxzzzzz_0_xyyzz_0, g_0_xxzzzzz_0_xyzzz_0, g_0_xxzzzzz_0_xzzzz_0, g_0_xxzzzzz_0_yyyyy_0, g_0_xxzzzzz_0_yyyyz_0, g_0_xxzzzzz_0_yyyzz_0, g_0_xxzzzzz_0_yyzzz_0, g_0_xxzzzzz_0_yzzzz_0, g_0_xxzzzzz_0_zzzzz_0, g_0_xzzzzz_0_xxxxz_0, g_0_xzzzzz_0_xxxxz_1, g_0_xzzzzz_0_xxxyz_0, g_0_xzzzzz_0_xxxyz_1, g_0_xzzzzz_0_xxxz_1, g_0_xzzzzz_0_xxxzz_0, g_0_xzzzzz_0_xxxzz_1, g_0_xzzzzz_0_xxyyz_0, g_0_xzzzzz_0_xxyyz_1, g_0_xzzzzz_0_xxyz_1, g_0_xzzzzz_0_xxyzz_0, g_0_xzzzzz_0_xxyzz_1, g_0_xzzzzz_0_xxzz_1, g_0_xzzzzz_0_xxzzz_0, g_0_xzzzzz_0_xxzzz_1, g_0_xzzzzz_0_xyyyz_0, g_0_xzzzzz_0_xyyyz_1, g_0_xzzzzz_0_xyyz_1, g_0_xzzzzz_0_xyyzz_0, g_0_xzzzzz_0_xyyzz_1, g_0_xzzzzz_0_xyzz_1, g_0_xzzzzz_0_xyzzz_0, g_0_xzzzzz_0_xyzzz_1, g_0_xzzzzz_0_xzzz_1, g_0_xzzzzz_0_xzzzz_0, g_0_xzzzzz_0_xzzzz_1, g_0_xzzzzz_0_yyyyy_0, g_0_xzzzzz_0_yyyyy_1, g_0_xzzzzz_0_yyyyz_0, g_0_xzzzzz_0_yyyyz_1, g_0_xzzzzz_0_yyyz_1, g_0_xzzzzz_0_yyyzz_0, g_0_xzzzzz_0_yyyzz_1, g_0_xzzzzz_0_yyzz_1, g_0_xzzzzz_0_yyzzz_0, g_0_xzzzzz_0_yyzzz_1, g_0_xzzzzz_0_yzzz_1, g_0_xzzzzz_0_yzzzz_0, g_0_xzzzzz_0_yzzzz_1, g_0_xzzzzz_0_zzzz_1, g_0_xzzzzz_0_zzzzz_0, g_0_xzzzzz_0_zzzzz_1, g_0_zzzzz_0_xxxxz_0, g_0_zzzzz_0_xxxxz_1, g_0_zzzzz_0_xxxyz_0, g_0_zzzzz_0_xxxyz_1, g_0_zzzzz_0_xxxzz_0, g_0_zzzzz_0_xxxzz_1, g_0_zzzzz_0_xxyyz_0, g_0_zzzzz_0_xxyyz_1, g_0_zzzzz_0_xxyzz_0, g_0_zzzzz_0_xxyzz_1, g_0_zzzzz_0_xxzzz_0, g_0_zzzzz_0_xxzzz_1, g_0_zzzzz_0_xyyyz_0, g_0_zzzzz_0_xyyyz_1, g_0_zzzzz_0_xyyzz_0, g_0_zzzzz_0_xyyzz_1, g_0_zzzzz_0_xyzzz_0, g_0_zzzzz_0_xyzzz_1, g_0_zzzzz_0_xzzzz_0, g_0_zzzzz_0_xzzzz_1, g_0_zzzzz_0_yyyyy_0, g_0_zzzzz_0_yyyyy_1, g_0_zzzzz_0_yyyyz_0, g_0_zzzzz_0_yyyyz_1, g_0_zzzzz_0_yyyzz_0, g_0_zzzzz_0_yyyzz_1, g_0_zzzzz_0_yyzzz_0, g_0_zzzzz_0_yyzzz_1, g_0_zzzzz_0_yzzzz_0, g_0_zzzzz_0_yzzzz_1, g_0_zzzzz_0_zzzzz_0, g_0_zzzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_xxxxx_0[i] = 4.0 * g_0_xxzzz_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxxx_0[i] * pb_z + g_0_xxzzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxy_0[i] = 4.0 * g_0_xxzzz_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxxy_0[i] * pb_z + g_0_xxzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxz_0[i] = g_0_zzzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxz_0[i] * pb_x + g_0_xzzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxyy_0[i] = 4.0 * g_0_xxzzz_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxyy_0[i] * pb_z + g_0_xxzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxyz_0[i] = g_0_zzzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxyz_0[i] * pb_x + g_0_xzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxzz_0[i] = g_0_zzzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxzz_0[i] * pb_x + g_0_xzzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyyy_0[i] = 4.0 * g_0_xxzzz_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxyyy_0[i] * pb_z + g_0_xxzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxyyz_0[i] = g_0_zzzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyyz_0[i] * pb_x + g_0_xzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyzz_0[i] = g_0_zzzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyzz_0[i] * pb_x + g_0_xzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxzzz_0[i] = g_0_zzzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxzzz_0[i] * pb_x + g_0_xzzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyyy_0[i] = 4.0 * g_0_xxzzz_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xyyyy_0[i] * pb_z + g_0_xxzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xyyyz_0[i] = g_0_zzzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyyyz_0[i] * pb_x + g_0_xzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyzz_0[i] = g_0_zzzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyyzz_0[i] * pb_x + g_0_xzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyzzz_0[i] = g_0_zzzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyzzz_0[i] * pb_x + g_0_xzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xzzzz_0[i] = g_0_zzzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xzzzz_0[i] * pb_x + g_0_xzzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyy_0[i] = g_0_zzzzz_0_yyyyy_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyy_0[i] * pb_x + g_0_xzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyz_0[i] = g_0_zzzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyz_0[i] * pb_x + g_0_xzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyzz_0[i] = g_0_zzzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyzz_0[i] * pb_x + g_0_xzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyzzz_0[i] = g_0_zzzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzzz_0[i] * pb_x + g_0_xzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yzzzz_0[i] = g_0_zzzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzzz_0[i] * pb_x + g_0_xzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_zzzzz_0[i] = g_0_zzzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzzz_0[i] * pb_x + g_0_xzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 441-462 components of targeted buffer : SKSH

    auto g_0_xyyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 441);

    auto g_0_xyyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 442);

    auto g_0_xyyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 443);

    auto g_0_xyyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 444);

    auto g_0_xyyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 445);

    auto g_0_xyyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 446);

    auto g_0_xyyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 447);

    auto g_0_xyyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 448);

    auto g_0_xyyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 449);

    auto g_0_xyyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 450);

    auto g_0_xyyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 451);

    auto g_0_xyyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 452);

    auto g_0_xyyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 453);

    auto g_0_xyyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 454);

    auto g_0_xyyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 455);

    auto g_0_xyyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 456);

    auto g_0_xyyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 457);

    auto g_0_xyyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 458);

    auto g_0_xyyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 459);

    auto g_0_xyyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 460);

    auto g_0_xyyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 461);

    #pragma omp simd aligned(g_0_xyyyyyy_0_xxxxx_0, g_0_xyyyyyy_0_xxxxy_0, g_0_xyyyyyy_0_xxxxz_0, g_0_xyyyyyy_0_xxxyy_0, g_0_xyyyyyy_0_xxxyz_0, g_0_xyyyyyy_0_xxxzz_0, g_0_xyyyyyy_0_xxyyy_0, g_0_xyyyyyy_0_xxyyz_0, g_0_xyyyyyy_0_xxyzz_0, g_0_xyyyyyy_0_xxzzz_0, g_0_xyyyyyy_0_xyyyy_0, g_0_xyyyyyy_0_xyyyz_0, g_0_xyyyyyy_0_xyyzz_0, g_0_xyyyyyy_0_xyzzz_0, g_0_xyyyyyy_0_xzzzz_0, g_0_xyyyyyy_0_yyyyy_0, g_0_xyyyyyy_0_yyyyz_0, g_0_xyyyyyy_0_yyyzz_0, g_0_xyyyyyy_0_yyzzz_0, g_0_xyyyyyy_0_yzzzz_0, g_0_xyyyyyy_0_zzzzz_0, g_0_yyyyyy_0_xxxx_1, g_0_yyyyyy_0_xxxxx_0, g_0_yyyyyy_0_xxxxx_1, g_0_yyyyyy_0_xxxxy_0, g_0_yyyyyy_0_xxxxy_1, g_0_yyyyyy_0_xxxxz_0, g_0_yyyyyy_0_xxxxz_1, g_0_yyyyyy_0_xxxy_1, g_0_yyyyyy_0_xxxyy_0, g_0_yyyyyy_0_xxxyy_1, g_0_yyyyyy_0_xxxyz_0, g_0_yyyyyy_0_xxxyz_1, g_0_yyyyyy_0_xxxz_1, g_0_yyyyyy_0_xxxzz_0, g_0_yyyyyy_0_xxxzz_1, g_0_yyyyyy_0_xxyy_1, g_0_yyyyyy_0_xxyyy_0, g_0_yyyyyy_0_xxyyy_1, g_0_yyyyyy_0_xxyyz_0, g_0_yyyyyy_0_xxyyz_1, g_0_yyyyyy_0_xxyz_1, g_0_yyyyyy_0_xxyzz_0, g_0_yyyyyy_0_xxyzz_1, g_0_yyyyyy_0_xxzz_1, g_0_yyyyyy_0_xxzzz_0, g_0_yyyyyy_0_xxzzz_1, g_0_yyyyyy_0_xyyy_1, g_0_yyyyyy_0_xyyyy_0, g_0_yyyyyy_0_xyyyy_1, g_0_yyyyyy_0_xyyyz_0, g_0_yyyyyy_0_xyyyz_1, g_0_yyyyyy_0_xyyz_1, g_0_yyyyyy_0_xyyzz_0, g_0_yyyyyy_0_xyyzz_1, g_0_yyyyyy_0_xyzz_1, g_0_yyyyyy_0_xyzzz_0, g_0_yyyyyy_0_xyzzz_1, g_0_yyyyyy_0_xzzz_1, g_0_yyyyyy_0_xzzzz_0, g_0_yyyyyy_0_xzzzz_1, g_0_yyyyyy_0_yyyy_1, g_0_yyyyyy_0_yyyyy_0, g_0_yyyyyy_0_yyyyy_1, g_0_yyyyyy_0_yyyyz_0, g_0_yyyyyy_0_yyyyz_1, g_0_yyyyyy_0_yyyz_1, g_0_yyyyyy_0_yyyzz_0, g_0_yyyyyy_0_yyyzz_1, g_0_yyyyyy_0_yyzz_1, g_0_yyyyyy_0_yyzzz_0, g_0_yyyyyy_0_yyzzz_1, g_0_yyyyyy_0_yzzz_1, g_0_yyyyyy_0_yzzzz_0, g_0_yyyyyy_0_yzzzz_1, g_0_yyyyyy_0_zzzz_1, g_0_yyyyyy_0_zzzzz_0, g_0_yyyyyy_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_xxxxx_0[i] = 5.0 * g_0_yyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxx_0[i] * pb_x + g_0_yyyyyy_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxy_0[i] = 4.0 * g_0_yyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxy_0[i] * pb_x + g_0_yyyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxz_0[i] = 4.0 * g_0_yyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxz_0[i] * pb_x + g_0_yyyyyy_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyy_0[i] = 3.0 * g_0_yyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyy_0[i] * pb_x + g_0_yyyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyz_0[i] = 3.0 * g_0_yyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyz_0[i] * pb_x + g_0_yyyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxzz_0[i] = 3.0 * g_0_yyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxzz_0[i] * pb_x + g_0_yyyyyy_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyy_0[i] = 2.0 * g_0_yyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyy_0[i] * pb_x + g_0_yyyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyz_0[i] = 2.0 * g_0_yyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyz_0[i] * pb_x + g_0_yyyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzz_0[i] * pb_x + g_0_yyyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxzzz_0[i] = 2.0 * g_0_yyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzzz_0[i] * pb_x + g_0_yyyyyy_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyy_0[i] = g_0_yyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyy_0[i] * pb_x + g_0_yyyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyz_0[i] = g_0_yyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyz_0[i] * pb_x + g_0_yyyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyzz_0[i] = g_0_yyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzz_0[i] * pb_x + g_0_yyyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyzzz_0[i] = g_0_yyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzz_0[i] * pb_x + g_0_yyyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xzzzz_0[i] = g_0_yyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzzz_0[i] * pb_x + g_0_yyyyyy_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyy_0[i] = g_0_yyyyyy_0_yyyyy_0[i] * pb_x + g_0_yyyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyz_0[i] = g_0_yyyyyy_0_yyyyz_0[i] * pb_x + g_0_yyyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyzz_0[i] = g_0_yyyyyy_0_yyyzz_0[i] * pb_x + g_0_yyyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyzzz_0[i] = g_0_yyyyyy_0_yyzzz_0[i] * pb_x + g_0_yyyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yzzzz_0[i] = g_0_yyyyyy_0_yzzzz_0[i] * pb_x + g_0_yyyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_zzzzz_0[i] = g_0_yyyyyy_0_zzzzz_0[i] * pb_x + g_0_yyyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 462-483 components of targeted buffer : SKSH

    auto g_0_xyyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 462);

    auto g_0_xyyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 463);

    auto g_0_xyyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 464);

    auto g_0_xyyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 465);

    auto g_0_xyyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 466);

    auto g_0_xyyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 467);

    auto g_0_xyyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 468);

    auto g_0_xyyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 469);

    auto g_0_xyyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 470);

    auto g_0_xyyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 471);

    auto g_0_xyyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 472);

    auto g_0_xyyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 473);

    auto g_0_xyyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 474);

    auto g_0_xyyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 475);

    auto g_0_xyyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 476);

    auto g_0_xyyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 477);

    auto g_0_xyyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 478);

    auto g_0_xyyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 479);

    auto g_0_xyyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 480);

    auto g_0_xyyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 481);

    auto g_0_xyyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 482);

    #pragma omp simd aligned(g_0_xyyyyy_0_xxxxx_0, g_0_xyyyyy_0_xxxxx_1, g_0_xyyyyy_0_xxxxy_0, g_0_xyyyyy_0_xxxxy_1, g_0_xyyyyy_0_xxxyy_0, g_0_xyyyyy_0_xxxyy_1, g_0_xyyyyy_0_xxyyy_0, g_0_xyyyyy_0_xxyyy_1, g_0_xyyyyy_0_xyyyy_0, g_0_xyyyyy_0_xyyyy_1, g_0_xyyyyyz_0_xxxxx_0, g_0_xyyyyyz_0_xxxxy_0, g_0_xyyyyyz_0_xxxxz_0, g_0_xyyyyyz_0_xxxyy_0, g_0_xyyyyyz_0_xxxyz_0, g_0_xyyyyyz_0_xxxzz_0, g_0_xyyyyyz_0_xxyyy_0, g_0_xyyyyyz_0_xxyyz_0, g_0_xyyyyyz_0_xxyzz_0, g_0_xyyyyyz_0_xxzzz_0, g_0_xyyyyyz_0_xyyyy_0, g_0_xyyyyyz_0_xyyyz_0, g_0_xyyyyyz_0_xyyzz_0, g_0_xyyyyyz_0_xyzzz_0, g_0_xyyyyyz_0_xzzzz_0, g_0_xyyyyyz_0_yyyyy_0, g_0_xyyyyyz_0_yyyyz_0, g_0_xyyyyyz_0_yyyzz_0, g_0_xyyyyyz_0_yyzzz_0, g_0_xyyyyyz_0_yzzzz_0, g_0_xyyyyyz_0_zzzzz_0, g_0_yyyyyz_0_xxxxz_0, g_0_yyyyyz_0_xxxxz_1, g_0_yyyyyz_0_xxxyz_0, g_0_yyyyyz_0_xxxyz_1, g_0_yyyyyz_0_xxxz_1, g_0_yyyyyz_0_xxxzz_0, g_0_yyyyyz_0_xxxzz_1, g_0_yyyyyz_0_xxyyz_0, g_0_yyyyyz_0_xxyyz_1, g_0_yyyyyz_0_xxyz_1, g_0_yyyyyz_0_xxyzz_0, g_0_yyyyyz_0_xxyzz_1, g_0_yyyyyz_0_xxzz_1, g_0_yyyyyz_0_xxzzz_0, g_0_yyyyyz_0_xxzzz_1, g_0_yyyyyz_0_xyyyz_0, g_0_yyyyyz_0_xyyyz_1, g_0_yyyyyz_0_xyyz_1, g_0_yyyyyz_0_xyyzz_0, g_0_yyyyyz_0_xyyzz_1, g_0_yyyyyz_0_xyzz_1, g_0_yyyyyz_0_xyzzz_0, g_0_yyyyyz_0_xyzzz_1, g_0_yyyyyz_0_xzzz_1, g_0_yyyyyz_0_xzzzz_0, g_0_yyyyyz_0_xzzzz_1, g_0_yyyyyz_0_yyyyy_0, g_0_yyyyyz_0_yyyyy_1, g_0_yyyyyz_0_yyyyz_0, g_0_yyyyyz_0_yyyyz_1, g_0_yyyyyz_0_yyyz_1, g_0_yyyyyz_0_yyyzz_0, g_0_yyyyyz_0_yyyzz_1, g_0_yyyyyz_0_yyzz_1, g_0_yyyyyz_0_yyzzz_0, g_0_yyyyyz_0_yyzzz_1, g_0_yyyyyz_0_yzzz_1, g_0_yyyyyz_0_yzzzz_0, g_0_yyyyyz_0_yzzzz_1, g_0_yyyyyz_0_zzzz_1, g_0_yyyyyz_0_zzzzz_0, g_0_yyyyyz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyz_0_xxxxx_0[i] = g_0_xyyyyy_0_xxxxx_0[i] * pb_z + g_0_xyyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxy_0[i] = g_0_xyyyyy_0_xxxxy_0[i] * pb_z + g_0_xyyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxz_0[i] = 4.0 * g_0_yyyyyz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxz_0[i] * pb_x + g_0_yyyyyz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxyy_0[i] = g_0_xyyyyy_0_xxxyy_0[i] * pb_z + g_0_xyyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxyz_0[i] = 3.0 * g_0_yyyyyz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxyz_0[i] * pb_x + g_0_yyyyyz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxzz_0[i] = 3.0 * g_0_yyyyyz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxzz_0[i] * pb_x + g_0_yyyyyz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyyy_0[i] = g_0_xyyyyy_0_xxyyy_0[i] * pb_z + g_0_xyyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxyyz_0[i] = 2.0 * g_0_yyyyyz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyyz_0[i] * pb_x + g_0_yyyyyz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyyyz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyzz_0[i] * pb_x + g_0_yyyyyz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxzzz_0[i] = 2.0 * g_0_yyyyyz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxzzz_0[i] * pb_x + g_0_yyyyyz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyyy_0[i] = g_0_xyyyyy_0_xyyyy_0[i] * pb_z + g_0_xyyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xyyyz_0[i] = g_0_yyyyyz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyyz_0[i] * pb_x + g_0_yyyyyz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyzz_0[i] = g_0_yyyyyz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyzz_0[i] * pb_x + g_0_yyyyyz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyzzz_0[i] = g_0_yyyyyz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyzzz_0[i] * pb_x + g_0_yyyyyz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xzzzz_0[i] = g_0_yyyyyz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xzzzz_0[i] * pb_x + g_0_yyyyyz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyy_0[i] = g_0_yyyyyz_0_yyyyy_0[i] * pb_x + g_0_yyyyyz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyz_0[i] = g_0_yyyyyz_0_yyyyz_0[i] * pb_x + g_0_yyyyyz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyzz_0[i] = g_0_yyyyyz_0_yyyzz_0[i] * pb_x + g_0_yyyyyz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyzzz_0[i] = g_0_yyyyyz_0_yyzzz_0[i] * pb_x + g_0_yyyyyz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yzzzz_0[i] = g_0_yyyyyz_0_yzzzz_0[i] * pb_x + g_0_yyyyyz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_zzzzz_0[i] = g_0_yyyyyz_0_zzzzz_0[i] * pb_x + g_0_yyyyyz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 483-504 components of targeted buffer : SKSH

    auto g_0_xyyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 483);

    auto g_0_xyyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 484);

    auto g_0_xyyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 485);

    auto g_0_xyyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 486);

    auto g_0_xyyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 487);

    auto g_0_xyyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 488);

    auto g_0_xyyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 489);

    auto g_0_xyyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 490);

    auto g_0_xyyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 491);

    auto g_0_xyyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 492);

    auto g_0_xyyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 493);

    auto g_0_xyyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 494);

    auto g_0_xyyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 495);

    auto g_0_xyyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 496);

    auto g_0_xyyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 497);

    auto g_0_xyyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 498);

    auto g_0_xyyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 499);

    auto g_0_xyyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 500);

    auto g_0_xyyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 501);

    auto g_0_xyyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 502);

    auto g_0_xyyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 503);

    #pragma omp simd aligned(g_0_xyyyyzz_0_xxxxx_0, g_0_xyyyyzz_0_xxxxy_0, g_0_xyyyyzz_0_xxxxz_0, g_0_xyyyyzz_0_xxxyy_0, g_0_xyyyyzz_0_xxxyz_0, g_0_xyyyyzz_0_xxxzz_0, g_0_xyyyyzz_0_xxyyy_0, g_0_xyyyyzz_0_xxyyz_0, g_0_xyyyyzz_0_xxyzz_0, g_0_xyyyyzz_0_xxzzz_0, g_0_xyyyyzz_0_xyyyy_0, g_0_xyyyyzz_0_xyyyz_0, g_0_xyyyyzz_0_xyyzz_0, g_0_xyyyyzz_0_xyzzz_0, g_0_xyyyyzz_0_xzzzz_0, g_0_xyyyyzz_0_yyyyy_0, g_0_xyyyyzz_0_yyyyz_0, g_0_xyyyyzz_0_yyyzz_0, g_0_xyyyyzz_0_yyzzz_0, g_0_xyyyyzz_0_yzzzz_0, g_0_xyyyyzz_0_zzzzz_0, g_0_yyyyzz_0_xxxx_1, g_0_yyyyzz_0_xxxxx_0, g_0_yyyyzz_0_xxxxx_1, g_0_yyyyzz_0_xxxxy_0, g_0_yyyyzz_0_xxxxy_1, g_0_yyyyzz_0_xxxxz_0, g_0_yyyyzz_0_xxxxz_1, g_0_yyyyzz_0_xxxy_1, g_0_yyyyzz_0_xxxyy_0, g_0_yyyyzz_0_xxxyy_1, g_0_yyyyzz_0_xxxyz_0, g_0_yyyyzz_0_xxxyz_1, g_0_yyyyzz_0_xxxz_1, g_0_yyyyzz_0_xxxzz_0, g_0_yyyyzz_0_xxxzz_1, g_0_yyyyzz_0_xxyy_1, g_0_yyyyzz_0_xxyyy_0, g_0_yyyyzz_0_xxyyy_1, g_0_yyyyzz_0_xxyyz_0, g_0_yyyyzz_0_xxyyz_1, g_0_yyyyzz_0_xxyz_1, g_0_yyyyzz_0_xxyzz_0, g_0_yyyyzz_0_xxyzz_1, g_0_yyyyzz_0_xxzz_1, g_0_yyyyzz_0_xxzzz_0, g_0_yyyyzz_0_xxzzz_1, g_0_yyyyzz_0_xyyy_1, g_0_yyyyzz_0_xyyyy_0, g_0_yyyyzz_0_xyyyy_1, g_0_yyyyzz_0_xyyyz_0, g_0_yyyyzz_0_xyyyz_1, g_0_yyyyzz_0_xyyz_1, g_0_yyyyzz_0_xyyzz_0, g_0_yyyyzz_0_xyyzz_1, g_0_yyyyzz_0_xyzz_1, g_0_yyyyzz_0_xyzzz_0, g_0_yyyyzz_0_xyzzz_1, g_0_yyyyzz_0_xzzz_1, g_0_yyyyzz_0_xzzzz_0, g_0_yyyyzz_0_xzzzz_1, g_0_yyyyzz_0_yyyy_1, g_0_yyyyzz_0_yyyyy_0, g_0_yyyyzz_0_yyyyy_1, g_0_yyyyzz_0_yyyyz_0, g_0_yyyyzz_0_yyyyz_1, g_0_yyyyzz_0_yyyz_1, g_0_yyyyzz_0_yyyzz_0, g_0_yyyyzz_0_yyyzz_1, g_0_yyyyzz_0_yyzz_1, g_0_yyyyzz_0_yyzzz_0, g_0_yyyyzz_0_yyzzz_1, g_0_yyyyzz_0_yzzz_1, g_0_yyyyzz_0_yzzzz_0, g_0_yyyyzz_0_yzzzz_1, g_0_yyyyzz_0_zzzz_1, g_0_yyyyzz_0_zzzzz_0, g_0_yyyyzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_xxxxx_0[i] = 5.0 * g_0_yyyyzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxx_0[i] * pb_x + g_0_yyyyzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxy_0[i] = 4.0 * g_0_yyyyzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxy_0[i] * pb_x + g_0_yyyyzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxz_0[i] = 4.0 * g_0_yyyyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxz_0[i] * pb_x + g_0_yyyyzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyy_0[i] * pb_x + g_0_yyyyzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyz_0[i] = 3.0 * g_0_yyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyz_0[i] * pb_x + g_0_yyyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxzz_0[i] = 3.0 * g_0_yyyyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxzz_0[i] * pb_x + g_0_yyyyzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyy_0[i] = 2.0 * g_0_yyyyzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyy_0[i] * pb_x + g_0_yyyyzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyz_0[i] = 2.0 * g_0_yyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyz_0[i] * pb_x + g_0_yyyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyzz_0[i] = 2.0 * g_0_yyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyzz_0[i] * pb_x + g_0_yyyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxzzz_0[i] = 2.0 * g_0_yyyyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxzzz_0[i] * pb_x + g_0_yyyyzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyy_0[i] = g_0_yyyyzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyy_0[i] * pb_x + g_0_yyyyzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyz_0[i] = g_0_yyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyz_0[i] * pb_x + g_0_yyyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyzz_0[i] = g_0_yyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyzz_0[i] * pb_x + g_0_yyyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyzzz_0[i] = g_0_yyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzzz_0[i] * pb_x + g_0_yyyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xzzzz_0[i] = g_0_yyyyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xzzzz_0[i] * pb_x + g_0_yyyyzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyy_0[i] = g_0_yyyyzz_0_yyyyy_0[i] * pb_x + g_0_yyyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyz_0[i] = g_0_yyyyzz_0_yyyyz_0[i] * pb_x + g_0_yyyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyzz_0[i] = g_0_yyyyzz_0_yyyzz_0[i] * pb_x + g_0_yyyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyzzz_0[i] = g_0_yyyyzz_0_yyzzz_0[i] * pb_x + g_0_yyyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yzzzz_0[i] = g_0_yyyyzz_0_yzzzz_0[i] * pb_x + g_0_yyyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_zzzzz_0[i] = g_0_yyyyzz_0_zzzzz_0[i] * pb_x + g_0_yyyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 504-525 components of targeted buffer : SKSH

    auto g_0_xyyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 504);

    auto g_0_xyyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 505);

    auto g_0_xyyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 506);

    auto g_0_xyyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 507);

    auto g_0_xyyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 508);

    auto g_0_xyyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 509);

    auto g_0_xyyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 510);

    auto g_0_xyyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 511);

    auto g_0_xyyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 512);

    auto g_0_xyyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 513);

    auto g_0_xyyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 514);

    auto g_0_xyyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 515);

    auto g_0_xyyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 516);

    auto g_0_xyyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 517);

    auto g_0_xyyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 518);

    auto g_0_xyyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 519);

    auto g_0_xyyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 520);

    auto g_0_xyyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 521);

    auto g_0_xyyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 522);

    auto g_0_xyyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 523);

    auto g_0_xyyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 524);

    #pragma omp simd aligned(g_0_xyyyzzz_0_xxxxx_0, g_0_xyyyzzz_0_xxxxy_0, g_0_xyyyzzz_0_xxxxz_0, g_0_xyyyzzz_0_xxxyy_0, g_0_xyyyzzz_0_xxxyz_0, g_0_xyyyzzz_0_xxxzz_0, g_0_xyyyzzz_0_xxyyy_0, g_0_xyyyzzz_0_xxyyz_0, g_0_xyyyzzz_0_xxyzz_0, g_0_xyyyzzz_0_xxzzz_0, g_0_xyyyzzz_0_xyyyy_0, g_0_xyyyzzz_0_xyyyz_0, g_0_xyyyzzz_0_xyyzz_0, g_0_xyyyzzz_0_xyzzz_0, g_0_xyyyzzz_0_xzzzz_0, g_0_xyyyzzz_0_yyyyy_0, g_0_xyyyzzz_0_yyyyz_0, g_0_xyyyzzz_0_yyyzz_0, g_0_xyyyzzz_0_yyzzz_0, g_0_xyyyzzz_0_yzzzz_0, g_0_xyyyzzz_0_zzzzz_0, g_0_yyyzzz_0_xxxx_1, g_0_yyyzzz_0_xxxxx_0, g_0_yyyzzz_0_xxxxx_1, g_0_yyyzzz_0_xxxxy_0, g_0_yyyzzz_0_xxxxy_1, g_0_yyyzzz_0_xxxxz_0, g_0_yyyzzz_0_xxxxz_1, g_0_yyyzzz_0_xxxy_1, g_0_yyyzzz_0_xxxyy_0, g_0_yyyzzz_0_xxxyy_1, g_0_yyyzzz_0_xxxyz_0, g_0_yyyzzz_0_xxxyz_1, g_0_yyyzzz_0_xxxz_1, g_0_yyyzzz_0_xxxzz_0, g_0_yyyzzz_0_xxxzz_1, g_0_yyyzzz_0_xxyy_1, g_0_yyyzzz_0_xxyyy_0, g_0_yyyzzz_0_xxyyy_1, g_0_yyyzzz_0_xxyyz_0, g_0_yyyzzz_0_xxyyz_1, g_0_yyyzzz_0_xxyz_1, g_0_yyyzzz_0_xxyzz_0, g_0_yyyzzz_0_xxyzz_1, g_0_yyyzzz_0_xxzz_1, g_0_yyyzzz_0_xxzzz_0, g_0_yyyzzz_0_xxzzz_1, g_0_yyyzzz_0_xyyy_1, g_0_yyyzzz_0_xyyyy_0, g_0_yyyzzz_0_xyyyy_1, g_0_yyyzzz_0_xyyyz_0, g_0_yyyzzz_0_xyyyz_1, g_0_yyyzzz_0_xyyz_1, g_0_yyyzzz_0_xyyzz_0, g_0_yyyzzz_0_xyyzz_1, g_0_yyyzzz_0_xyzz_1, g_0_yyyzzz_0_xyzzz_0, g_0_yyyzzz_0_xyzzz_1, g_0_yyyzzz_0_xzzz_1, g_0_yyyzzz_0_xzzzz_0, g_0_yyyzzz_0_xzzzz_1, g_0_yyyzzz_0_yyyy_1, g_0_yyyzzz_0_yyyyy_0, g_0_yyyzzz_0_yyyyy_1, g_0_yyyzzz_0_yyyyz_0, g_0_yyyzzz_0_yyyyz_1, g_0_yyyzzz_0_yyyz_1, g_0_yyyzzz_0_yyyzz_0, g_0_yyyzzz_0_yyyzz_1, g_0_yyyzzz_0_yyzz_1, g_0_yyyzzz_0_yyzzz_0, g_0_yyyzzz_0_yyzzz_1, g_0_yyyzzz_0_yzzz_1, g_0_yyyzzz_0_yzzzz_0, g_0_yyyzzz_0_yzzzz_1, g_0_yyyzzz_0_zzzz_1, g_0_yyyzzz_0_zzzzz_0, g_0_yyyzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_xxxxx_0[i] = 5.0 * g_0_yyyzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxx_0[i] * pb_x + g_0_yyyzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxy_0[i] = 4.0 * g_0_yyyzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxy_0[i] * pb_x + g_0_yyyzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxz_0[i] = 4.0 * g_0_yyyzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxz_0[i] * pb_x + g_0_yyyzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyy_0[i] = 3.0 * g_0_yyyzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyy_0[i] * pb_x + g_0_yyyzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyz_0[i] = 3.0 * g_0_yyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyz_0[i] * pb_x + g_0_yyyzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxzz_0[i] = 3.0 * g_0_yyyzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxzz_0[i] * pb_x + g_0_yyyzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyyzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyy_0[i] * pb_x + g_0_yyyzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyz_0[i] = 2.0 * g_0_yyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyz_0[i] * pb_x + g_0_yyyzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyzz_0[i] = 2.0 * g_0_yyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyzz_0[i] * pb_x + g_0_yyyzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxzzz_0[i] = 2.0 * g_0_yyyzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxzzz_0[i] * pb_x + g_0_yyyzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyy_0[i] = g_0_yyyzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyy_0[i] * pb_x + g_0_yyyzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyz_0[i] = g_0_yyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyz_0[i] * pb_x + g_0_yyyzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyzz_0[i] = g_0_yyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyzz_0[i] * pb_x + g_0_yyyzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyzzz_0[i] = g_0_yyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzzz_0[i] * pb_x + g_0_yyyzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xzzzz_0[i] = g_0_yyyzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xzzzz_0[i] * pb_x + g_0_yyyzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyy_0[i] = g_0_yyyzzz_0_yyyyy_0[i] * pb_x + g_0_yyyzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyz_0[i] = g_0_yyyzzz_0_yyyyz_0[i] * pb_x + g_0_yyyzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyzz_0[i] = g_0_yyyzzz_0_yyyzz_0[i] * pb_x + g_0_yyyzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyzzz_0[i] = g_0_yyyzzz_0_yyzzz_0[i] * pb_x + g_0_yyyzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yzzzz_0[i] = g_0_yyyzzz_0_yzzzz_0[i] * pb_x + g_0_yyyzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_zzzzz_0[i] = g_0_yyyzzz_0_zzzzz_0[i] * pb_x + g_0_yyyzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 525-546 components of targeted buffer : SKSH

    auto g_0_xyyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 525);

    auto g_0_xyyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 526);

    auto g_0_xyyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 527);

    auto g_0_xyyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 528);

    auto g_0_xyyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 529);

    auto g_0_xyyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 530);

    auto g_0_xyyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 531);

    auto g_0_xyyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 532);

    auto g_0_xyyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 533);

    auto g_0_xyyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 534);

    auto g_0_xyyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 535);

    auto g_0_xyyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 536);

    auto g_0_xyyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 537);

    auto g_0_xyyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 538);

    auto g_0_xyyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 539);

    auto g_0_xyyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 540);

    auto g_0_xyyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 541);

    auto g_0_xyyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 542);

    auto g_0_xyyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 543);

    auto g_0_xyyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 544);

    auto g_0_xyyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 545);

    #pragma omp simd aligned(g_0_xyyzzzz_0_xxxxx_0, g_0_xyyzzzz_0_xxxxy_0, g_0_xyyzzzz_0_xxxxz_0, g_0_xyyzzzz_0_xxxyy_0, g_0_xyyzzzz_0_xxxyz_0, g_0_xyyzzzz_0_xxxzz_0, g_0_xyyzzzz_0_xxyyy_0, g_0_xyyzzzz_0_xxyyz_0, g_0_xyyzzzz_0_xxyzz_0, g_0_xyyzzzz_0_xxzzz_0, g_0_xyyzzzz_0_xyyyy_0, g_0_xyyzzzz_0_xyyyz_0, g_0_xyyzzzz_0_xyyzz_0, g_0_xyyzzzz_0_xyzzz_0, g_0_xyyzzzz_0_xzzzz_0, g_0_xyyzzzz_0_yyyyy_0, g_0_xyyzzzz_0_yyyyz_0, g_0_xyyzzzz_0_yyyzz_0, g_0_xyyzzzz_0_yyzzz_0, g_0_xyyzzzz_0_yzzzz_0, g_0_xyyzzzz_0_zzzzz_0, g_0_yyzzzz_0_xxxx_1, g_0_yyzzzz_0_xxxxx_0, g_0_yyzzzz_0_xxxxx_1, g_0_yyzzzz_0_xxxxy_0, g_0_yyzzzz_0_xxxxy_1, g_0_yyzzzz_0_xxxxz_0, g_0_yyzzzz_0_xxxxz_1, g_0_yyzzzz_0_xxxy_1, g_0_yyzzzz_0_xxxyy_0, g_0_yyzzzz_0_xxxyy_1, g_0_yyzzzz_0_xxxyz_0, g_0_yyzzzz_0_xxxyz_1, g_0_yyzzzz_0_xxxz_1, g_0_yyzzzz_0_xxxzz_0, g_0_yyzzzz_0_xxxzz_1, g_0_yyzzzz_0_xxyy_1, g_0_yyzzzz_0_xxyyy_0, g_0_yyzzzz_0_xxyyy_1, g_0_yyzzzz_0_xxyyz_0, g_0_yyzzzz_0_xxyyz_1, g_0_yyzzzz_0_xxyz_1, g_0_yyzzzz_0_xxyzz_0, g_0_yyzzzz_0_xxyzz_1, g_0_yyzzzz_0_xxzz_1, g_0_yyzzzz_0_xxzzz_0, g_0_yyzzzz_0_xxzzz_1, g_0_yyzzzz_0_xyyy_1, g_0_yyzzzz_0_xyyyy_0, g_0_yyzzzz_0_xyyyy_1, g_0_yyzzzz_0_xyyyz_0, g_0_yyzzzz_0_xyyyz_1, g_0_yyzzzz_0_xyyz_1, g_0_yyzzzz_0_xyyzz_0, g_0_yyzzzz_0_xyyzz_1, g_0_yyzzzz_0_xyzz_1, g_0_yyzzzz_0_xyzzz_0, g_0_yyzzzz_0_xyzzz_1, g_0_yyzzzz_0_xzzz_1, g_0_yyzzzz_0_xzzzz_0, g_0_yyzzzz_0_xzzzz_1, g_0_yyzzzz_0_yyyy_1, g_0_yyzzzz_0_yyyyy_0, g_0_yyzzzz_0_yyyyy_1, g_0_yyzzzz_0_yyyyz_0, g_0_yyzzzz_0_yyyyz_1, g_0_yyzzzz_0_yyyz_1, g_0_yyzzzz_0_yyyzz_0, g_0_yyzzzz_0_yyyzz_1, g_0_yyzzzz_0_yyzz_1, g_0_yyzzzz_0_yyzzz_0, g_0_yyzzzz_0_yyzzz_1, g_0_yyzzzz_0_yzzz_1, g_0_yyzzzz_0_yzzzz_0, g_0_yyzzzz_0_yzzzz_1, g_0_yyzzzz_0_zzzz_1, g_0_yyzzzz_0_zzzzz_0, g_0_yyzzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_xxxxx_0[i] = 5.0 * g_0_yyzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxx_0[i] * pb_x + g_0_yyzzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxy_0[i] = 4.0 * g_0_yyzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxy_0[i] * pb_x + g_0_yyzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxz_0[i] = 4.0 * g_0_yyzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxz_0[i] * pb_x + g_0_yyzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyy_0[i] = 3.0 * g_0_yyzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyy_0[i] * pb_x + g_0_yyzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyz_0[i] * pb_x + g_0_yyzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxzz_0[i] * pb_x + g_0_yyzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyy_0[i] = 2.0 * g_0_yyzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyy_0[i] * pb_x + g_0_yyzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyz_0[i] = 2.0 * g_0_yyzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyz_0[i] * pb_x + g_0_yyzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyzz_0[i] = 2.0 * g_0_yyzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyzz_0[i] * pb_x + g_0_yyzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxzzz_0[i] = 2.0 * g_0_yyzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxzzz_0[i] * pb_x + g_0_yyzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyy_0[i] = g_0_yyzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyy_0[i] * pb_x + g_0_yyzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyz_0[i] = g_0_yyzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyz_0[i] * pb_x + g_0_yyzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyzz_0[i] = g_0_yyzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyzz_0[i] * pb_x + g_0_yyzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyzzz_0[i] = g_0_yyzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzzz_0[i] * pb_x + g_0_yyzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xzzzz_0[i] = g_0_yyzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xzzzz_0[i] * pb_x + g_0_yyzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyy_0[i] = g_0_yyzzzz_0_yyyyy_0[i] * pb_x + g_0_yyzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyz_0[i] = g_0_yyzzzz_0_yyyyz_0[i] * pb_x + g_0_yyzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyzz_0[i] = g_0_yyzzzz_0_yyyzz_0[i] * pb_x + g_0_yyzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyzzz_0[i] = g_0_yyzzzz_0_yyzzz_0[i] * pb_x + g_0_yyzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yzzzz_0[i] = g_0_yyzzzz_0_yzzzz_0[i] * pb_x + g_0_yyzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_zzzzz_0[i] = g_0_yyzzzz_0_zzzzz_0[i] * pb_x + g_0_yyzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 546-567 components of targeted buffer : SKSH

    auto g_0_xyzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 546);

    auto g_0_xyzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 547);

    auto g_0_xyzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 548);

    auto g_0_xyzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 549);

    auto g_0_xyzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 550);

    auto g_0_xyzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 551);

    auto g_0_xyzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 552);

    auto g_0_xyzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 553);

    auto g_0_xyzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 554);

    auto g_0_xyzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 555);

    auto g_0_xyzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 556);

    auto g_0_xyzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 557);

    auto g_0_xyzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 558);

    auto g_0_xyzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 559);

    auto g_0_xyzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 560);

    auto g_0_xyzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 561);

    auto g_0_xyzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 562);

    auto g_0_xyzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 563);

    auto g_0_xyzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 564);

    auto g_0_xyzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 565);

    auto g_0_xyzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 566);

    #pragma omp simd aligned(g_0_xyzzzzz_0_xxxxx_0, g_0_xyzzzzz_0_xxxxy_0, g_0_xyzzzzz_0_xxxxz_0, g_0_xyzzzzz_0_xxxyy_0, g_0_xyzzzzz_0_xxxyz_0, g_0_xyzzzzz_0_xxxzz_0, g_0_xyzzzzz_0_xxyyy_0, g_0_xyzzzzz_0_xxyyz_0, g_0_xyzzzzz_0_xxyzz_0, g_0_xyzzzzz_0_xxzzz_0, g_0_xyzzzzz_0_xyyyy_0, g_0_xyzzzzz_0_xyyyz_0, g_0_xyzzzzz_0_xyyzz_0, g_0_xyzzzzz_0_xyzzz_0, g_0_xyzzzzz_0_xzzzz_0, g_0_xyzzzzz_0_yyyyy_0, g_0_xyzzzzz_0_yyyyz_0, g_0_xyzzzzz_0_yyyzz_0, g_0_xyzzzzz_0_yyzzz_0, g_0_xyzzzzz_0_yzzzz_0, g_0_xyzzzzz_0_zzzzz_0, g_0_xzzzzz_0_xxxxx_0, g_0_xzzzzz_0_xxxxx_1, g_0_xzzzzz_0_xxxxz_0, g_0_xzzzzz_0_xxxxz_1, g_0_xzzzzz_0_xxxzz_0, g_0_xzzzzz_0_xxxzz_1, g_0_xzzzzz_0_xxzzz_0, g_0_xzzzzz_0_xxzzz_1, g_0_xzzzzz_0_xzzzz_0, g_0_xzzzzz_0_xzzzz_1, g_0_yzzzzz_0_xxxxy_0, g_0_yzzzzz_0_xxxxy_1, g_0_yzzzzz_0_xxxy_1, g_0_yzzzzz_0_xxxyy_0, g_0_yzzzzz_0_xxxyy_1, g_0_yzzzzz_0_xxxyz_0, g_0_yzzzzz_0_xxxyz_1, g_0_yzzzzz_0_xxyy_1, g_0_yzzzzz_0_xxyyy_0, g_0_yzzzzz_0_xxyyy_1, g_0_yzzzzz_0_xxyyz_0, g_0_yzzzzz_0_xxyyz_1, g_0_yzzzzz_0_xxyz_1, g_0_yzzzzz_0_xxyzz_0, g_0_yzzzzz_0_xxyzz_1, g_0_yzzzzz_0_xyyy_1, g_0_yzzzzz_0_xyyyy_0, g_0_yzzzzz_0_xyyyy_1, g_0_yzzzzz_0_xyyyz_0, g_0_yzzzzz_0_xyyyz_1, g_0_yzzzzz_0_xyyz_1, g_0_yzzzzz_0_xyyzz_0, g_0_yzzzzz_0_xyyzz_1, g_0_yzzzzz_0_xyzz_1, g_0_yzzzzz_0_xyzzz_0, g_0_yzzzzz_0_xyzzz_1, g_0_yzzzzz_0_yyyy_1, g_0_yzzzzz_0_yyyyy_0, g_0_yzzzzz_0_yyyyy_1, g_0_yzzzzz_0_yyyyz_0, g_0_yzzzzz_0_yyyyz_1, g_0_yzzzzz_0_yyyz_1, g_0_yzzzzz_0_yyyzz_0, g_0_yzzzzz_0_yyyzz_1, g_0_yzzzzz_0_yyzz_1, g_0_yzzzzz_0_yyzzz_0, g_0_yzzzzz_0_yyzzz_1, g_0_yzzzzz_0_yzzz_1, g_0_yzzzzz_0_yzzzz_0, g_0_yzzzzz_0_yzzzz_1, g_0_yzzzzz_0_zzzzz_0, g_0_yzzzzz_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzz_0_xxxxx_0[i] = g_0_xzzzzz_0_xxxxx_0[i] * pb_y + g_0_xzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxxy_0[i] = 4.0 * g_0_yzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxy_0[i] * pb_x + g_0_yzzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxz_0[i] = g_0_xzzzzz_0_xxxxz_0[i] * pb_y + g_0_xzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxyy_0[i] = 3.0 * g_0_yzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyy_0[i] * pb_x + g_0_yzzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxyz_0[i] = 3.0 * g_0_yzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyz_0[i] * pb_x + g_0_yzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxzz_0[i] = g_0_xzzzzz_0_xxxzz_0[i] * pb_y + g_0_xzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxyyy_0[i] = 2.0 * g_0_yzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyy_0[i] * pb_x + g_0_yzzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyz_0[i] * pb_x + g_0_yzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyzz_0[i] * pb_x + g_0_yzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxzzz_0[i] = g_0_xzzzzz_0_xxzzz_0[i] * pb_y + g_0_xzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xyyyy_0[i] = g_0_yzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyy_0[i] * pb_x + g_0_yzzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyyz_0[i] = g_0_yzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyz_0[i] * pb_x + g_0_yzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyzz_0[i] = g_0_yzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyzz_0[i] * pb_x + g_0_yzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyzzz_0[i] = g_0_yzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyzzz_0[i] * pb_x + g_0_yzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xzzzz_0[i] = g_0_xzzzzz_0_xzzzz_0[i] * pb_y + g_0_xzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_yyyyy_0[i] = g_0_yzzzzz_0_yyyyy_0[i] * pb_x + g_0_yzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyyz_0[i] = g_0_yzzzzz_0_yyyyz_0[i] * pb_x + g_0_yzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyzz_0[i] = g_0_yzzzzz_0_yyyzz_0[i] * pb_x + g_0_yzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyzzz_0[i] = g_0_yzzzzz_0_yyzzz_0[i] * pb_x + g_0_yzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yzzzz_0[i] = g_0_yzzzzz_0_yzzzz_0[i] * pb_x + g_0_yzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_zzzzz_0[i] = g_0_yzzzzz_0_zzzzz_0[i] * pb_x + g_0_yzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 567-588 components of targeted buffer : SKSH

    auto g_0_xzzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 567);

    auto g_0_xzzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 568);

    auto g_0_xzzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 569);

    auto g_0_xzzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 570);

    auto g_0_xzzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 571);

    auto g_0_xzzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 572);

    auto g_0_xzzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 573);

    auto g_0_xzzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 574);

    auto g_0_xzzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 575);

    auto g_0_xzzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 576);

    auto g_0_xzzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 577);

    auto g_0_xzzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 578);

    auto g_0_xzzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 579);

    auto g_0_xzzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 580);

    auto g_0_xzzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 581);

    auto g_0_xzzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 582);

    auto g_0_xzzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 583);

    auto g_0_xzzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 584);

    auto g_0_xzzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 585);

    auto g_0_xzzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 586);

    auto g_0_xzzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 587);

    #pragma omp simd aligned(g_0_xzzzzzz_0_xxxxx_0, g_0_xzzzzzz_0_xxxxy_0, g_0_xzzzzzz_0_xxxxz_0, g_0_xzzzzzz_0_xxxyy_0, g_0_xzzzzzz_0_xxxyz_0, g_0_xzzzzzz_0_xxxzz_0, g_0_xzzzzzz_0_xxyyy_0, g_0_xzzzzzz_0_xxyyz_0, g_0_xzzzzzz_0_xxyzz_0, g_0_xzzzzzz_0_xxzzz_0, g_0_xzzzzzz_0_xyyyy_0, g_0_xzzzzzz_0_xyyyz_0, g_0_xzzzzzz_0_xyyzz_0, g_0_xzzzzzz_0_xyzzz_0, g_0_xzzzzzz_0_xzzzz_0, g_0_xzzzzzz_0_yyyyy_0, g_0_xzzzzzz_0_yyyyz_0, g_0_xzzzzzz_0_yyyzz_0, g_0_xzzzzzz_0_yyzzz_0, g_0_xzzzzzz_0_yzzzz_0, g_0_xzzzzzz_0_zzzzz_0, g_0_zzzzzz_0_xxxx_1, g_0_zzzzzz_0_xxxxx_0, g_0_zzzzzz_0_xxxxx_1, g_0_zzzzzz_0_xxxxy_0, g_0_zzzzzz_0_xxxxy_1, g_0_zzzzzz_0_xxxxz_0, g_0_zzzzzz_0_xxxxz_1, g_0_zzzzzz_0_xxxy_1, g_0_zzzzzz_0_xxxyy_0, g_0_zzzzzz_0_xxxyy_1, g_0_zzzzzz_0_xxxyz_0, g_0_zzzzzz_0_xxxyz_1, g_0_zzzzzz_0_xxxz_1, g_0_zzzzzz_0_xxxzz_0, g_0_zzzzzz_0_xxxzz_1, g_0_zzzzzz_0_xxyy_1, g_0_zzzzzz_0_xxyyy_0, g_0_zzzzzz_0_xxyyy_1, g_0_zzzzzz_0_xxyyz_0, g_0_zzzzzz_0_xxyyz_1, g_0_zzzzzz_0_xxyz_1, g_0_zzzzzz_0_xxyzz_0, g_0_zzzzzz_0_xxyzz_1, g_0_zzzzzz_0_xxzz_1, g_0_zzzzzz_0_xxzzz_0, g_0_zzzzzz_0_xxzzz_1, g_0_zzzzzz_0_xyyy_1, g_0_zzzzzz_0_xyyyy_0, g_0_zzzzzz_0_xyyyy_1, g_0_zzzzzz_0_xyyyz_0, g_0_zzzzzz_0_xyyyz_1, g_0_zzzzzz_0_xyyz_1, g_0_zzzzzz_0_xyyzz_0, g_0_zzzzzz_0_xyyzz_1, g_0_zzzzzz_0_xyzz_1, g_0_zzzzzz_0_xyzzz_0, g_0_zzzzzz_0_xyzzz_1, g_0_zzzzzz_0_xzzz_1, g_0_zzzzzz_0_xzzzz_0, g_0_zzzzzz_0_xzzzz_1, g_0_zzzzzz_0_yyyy_1, g_0_zzzzzz_0_yyyyy_0, g_0_zzzzzz_0_yyyyy_1, g_0_zzzzzz_0_yyyyz_0, g_0_zzzzzz_0_yyyyz_1, g_0_zzzzzz_0_yyyz_1, g_0_zzzzzz_0_yyyzz_0, g_0_zzzzzz_0_yyyzz_1, g_0_zzzzzz_0_yyzz_1, g_0_zzzzzz_0_yyzzz_0, g_0_zzzzzz_0_yyzzz_1, g_0_zzzzzz_0_yzzz_1, g_0_zzzzzz_0_yzzzz_0, g_0_zzzzzz_0_yzzzz_1, g_0_zzzzzz_0_zzzz_1, g_0_zzzzzz_0_zzzzz_0, g_0_zzzzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_xxxxx_0[i] = 5.0 * g_0_zzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxx_0[i] * pb_x + g_0_zzzzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxy_0[i] = 4.0 * g_0_zzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxy_0[i] * pb_x + g_0_zzzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxz_0[i] = 4.0 * g_0_zzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxz_0[i] * pb_x + g_0_zzzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyy_0[i] = 3.0 * g_0_zzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyy_0[i] * pb_x + g_0_zzzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyz_0[i] = 3.0 * g_0_zzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyz_0[i] * pb_x + g_0_zzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxzz_0[i] = 3.0 * g_0_zzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxzz_0[i] * pb_x + g_0_zzzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyy_0[i] = 2.0 * g_0_zzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyy_0[i] * pb_x + g_0_zzzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyz_0[i] * pb_x + g_0_zzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyzz_0[i] = 2.0 * g_0_zzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzz_0[i] * pb_x + g_0_zzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxzzz_0[i] = 2.0 * g_0_zzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzzz_0[i] * pb_x + g_0_zzzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyy_0[i] = g_0_zzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyy_0[i] * pb_x + g_0_zzzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyz_0[i] = g_0_zzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyz_0[i] * pb_x + g_0_zzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyzz_0[i] = g_0_zzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzz_0[i] * pb_x + g_0_zzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyzzz_0[i] = g_0_zzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzz_0[i] * pb_x + g_0_zzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xzzzz_0[i] = g_0_zzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzzz_0[i] * pb_x + g_0_zzzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyy_0[i] = g_0_zzzzzz_0_yyyyy_0[i] * pb_x + g_0_zzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyz_0[i] = g_0_zzzzzz_0_yyyyz_0[i] * pb_x + g_0_zzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyzz_0[i] = g_0_zzzzzz_0_yyyzz_0[i] * pb_x + g_0_zzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyzzz_0[i] = g_0_zzzzzz_0_yyzzz_0[i] * pb_x + g_0_zzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yzzzz_0[i] = g_0_zzzzzz_0_yzzzz_0[i] * pb_x + g_0_zzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_zzzzz_0[i] = g_0_zzzzzz_0_zzzzz_0[i] * pb_x + g_0_zzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 588-609 components of targeted buffer : SKSH

    auto g_0_yyyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 588);

    auto g_0_yyyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 589);

    auto g_0_yyyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 590);

    auto g_0_yyyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 591);

    auto g_0_yyyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 592);

    auto g_0_yyyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 593);

    auto g_0_yyyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 594);

    auto g_0_yyyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 595);

    auto g_0_yyyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 596);

    auto g_0_yyyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 597);

    auto g_0_yyyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 598);

    auto g_0_yyyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 599);

    auto g_0_yyyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 600);

    auto g_0_yyyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 601);

    auto g_0_yyyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 602);

    auto g_0_yyyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 603);

    auto g_0_yyyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 604);

    auto g_0_yyyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 605);

    auto g_0_yyyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 606);

    auto g_0_yyyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 607);

    auto g_0_yyyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 608);

    #pragma omp simd aligned(g_0_yyyyy_0_xxxxx_0, g_0_yyyyy_0_xxxxx_1, g_0_yyyyy_0_xxxxy_0, g_0_yyyyy_0_xxxxy_1, g_0_yyyyy_0_xxxxz_0, g_0_yyyyy_0_xxxxz_1, g_0_yyyyy_0_xxxyy_0, g_0_yyyyy_0_xxxyy_1, g_0_yyyyy_0_xxxyz_0, g_0_yyyyy_0_xxxyz_1, g_0_yyyyy_0_xxxzz_0, g_0_yyyyy_0_xxxzz_1, g_0_yyyyy_0_xxyyy_0, g_0_yyyyy_0_xxyyy_1, g_0_yyyyy_0_xxyyz_0, g_0_yyyyy_0_xxyyz_1, g_0_yyyyy_0_xxyzz_0, g_0_yyyyy_0_xxyzz_1, g_0_yyyyy_0_xxzzz_0, g_0_yyyyy_0_xxzzz_1, g_0_yyyyy_0_xyyyy_0, g_0_yyyyy_0_xyyyy_1, g_0_yyyyy_0_xyyyz_0, g_0_yyyyy_0_xyyyz_1, g_0_yyyyy_0_xyyzz_0, g_0_yyyyy_0_xyyzz_1, g_0_yyyyy_0_xyzzz_0, g_0_yyyyy_0_xyzzz_1, g_0_yyyyy_0_xzzzz_0, g_0_yyyyy_0_xzzzz_1, g_0_yyyyy_0_yyyyy_0, g_0_yyyyy_0_yyyyy_1, g_0_yyyyy_0_yyyyz_0, g_0_yyyyy_0_yyyyz_1, g_0_yyyyy_0_yyyzz_0, g_0_yyyyy_0_yyyzz_1, g_0_yyyyy_0_yyzzz_0, g_0_yyyyy_0_yyzzz_1, g_0_yyyyy_0_yzzzz_0, g_0_yyyyy_0_yzzzz_1, g_0_yyyyy_0_zzzzz_0, g_0_yyyyy_0_zzzzz_1, g_0_yyyyyy_0_xxxx_1, g_0_yyyyyy_0_xxxxx_0, g_0_yyyyyy_0_xxxxx_1, g_0_yyyyyy_0_xxxxy_0, g_0_yyyyyy_0_xxxxy_1, g_0_yyyyyy_0_xxxxz_0, g_0_yyyyyy_0_xxxxz_1, g_0_yyyyyy_0_xxxy_1, g_0_yyyyyy_0_xxxyy_0, g_0_yyyyyy_0_xxxyy_1, g_0_yyyyyy_0_xxxyz_0, g_0_yyyyyy_0_xxxyz_1, g_0_yyyyyy_0_xxxz_1, g_0_yyyyyy_0_xxxzz_0, g_0_yyyyyy_0_xxxzz_1, g_0_yyyyyy_0_xxyy_1, g_0_yyyyyy_0_xxyyy_0, g_0_yyyyyy_0_xxyyy_1, g_0_yyyyyy_0_xxyyz_0, g_0_yyyyyy_0_xxyyz_1, g_0_yyyyyy_0_xxyz_1, g_0_yyyyyy_0_xxyzz_0, g_0_yyyyyy_0_xxyzz_1, g_0_yyyyyy_0_xxzz_1, g_0_yyyyyy_0_xxzzz_0, g_0_yyyyyy_0_xxzzz_1, g_0_yyyyyy_0_xyyy_1, g_0_yyyyyy_0_xyyyy_0, g_0_yyyyyy_0_xyyyy_1, g_0_yyyyyy_0_xyyyz_0, g_0_yyyyyy_0_xyyyz_1, g_0_yyyyyy_0_xyyz_1, g_0_yyyyyy_0_xyyzz_0, g_0_yyyyyy_0_xyyzz_1, g_0_yyyyyy_0_xyzz_1, g_0_yyyyyy_0_xyzzz_0, g_0_yyyyyy_0_xyzzz_1, g_0_yyyyyy_0_xzzz_1, g_0_yyyyyy_0_xzzzz_0, g_0_yyyyyy_0_xzzzz_1, g_0_yyyyyy_0_yyyy_1, g_0_yyyyyy_0_yyyyy_0, g_0_yyyyyy_0_yyyyy_1, g_0_yyyyyy_0_yyyyz_0, g_0_yyyyyy_0_yyyyz_1, g_0_yyyyyy_0_yyyz_1, g_0_yyyyyy_0_yyyzz_0, g_0_yyyyyy_0_yyyzz_1, g_0_yyyyyy_0_yyzz_1, g_0_yyyyyy_0_yyzzz_0, g_0_yyyyyy_0_yyzzz_1, g_0_yyyyyy_0_yzzz_1, g_0_yyyyyy_0_yzzzz_0, g_0_yyyyyy_0_yzzzz_1, g_0_yyyyyy_0_zzzz_1, g_0_yyyyyy_0_zzzzz_0, g_0_yyyyyy_0_zzzzz_1, g_0_yyyyyyy_0_xxxxx_0, g_0_yyyyyyy_0_xxxxy_0, g_0_yyyyyyy_0_xxxxz_0, g_0_yyyyyyy_0_xxxyy_0, g_0_yyyyyyy_0_xxxyz_0, g_0_yyyyyyy_0_xxxzz_0, g_0_yyyyyyy_0_xxyyy_0, g_0_yyyyyyy_0_xxyyz_0, g_0_yyyyyyy_0_xxyzz_0, g_0_yyyyyyy_0_xxzzz_0, g_0_yyyyyyy_0_xyyyy_0, g_0_yyyyyyy_0_xyyyz_0, g_0_yyyyyyy_0_xyyzz_0, g_0_yyyyyyy_0_xyzzz_0, g_0_yyyyyyy_0_xzzzz_0, g_0_yyyyyyy_0_yyyyy_0, g_0_yyyyyyy_0_yyyyz_0, g_0_yyyyyyy_0_yyyzz_0, g_0_yyyyyyy_0_yyzzz_0, g_0_yyyyyyy_0_yzzzz_0, g_0_yyyyyyy_0_zzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_xxxxx_0[i] = 6.0 * g_0_yyyyy_0_xxxxx_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxx_0[i] * pb_y + g_0_yyyyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxy_0[i] = 6.0 * g_0_yyyyy_0_xxxxy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxy_0[i] * pb_y + g_0_yyyyyy_0_xxxxy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxz_0[i] = 6.0 * g_0_yyyyy_0_xxxxz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxz_0[i] * pb_y + g_0_yyyyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyy_0[i] = 6.0 * g_0_yyyyy_0_xxxyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyy_0[i] * pb_y + g_0_yyyyyy_0_xxxyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyz_0[i] = 6.0 * g_0_yyyyy_0_xxxyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyz_0[i] * pb_y + g_0_yyyyyy_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxzz_0[i] = 6.0 * g_0_yyyyy_0_xxxzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxzz_0[i] * pb_y + g_0_yyyyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyy_0[i] = 6.0 * g_0_yyyyy_0_xxyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyy_0[i] * pb_y + g_0_yyyyyy_0_xxyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyz_0[i] = 6.0 * g_0_yyyyy_0_xxyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyz_0[i] * pb_y + g_0_yyyyyy_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyzz_0[i] = 6.0 * g_0_yyyyy_0_xxyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzz_0[i] * pb_y + g_0_yyyyyy_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxzzz_0[i] = 6.0 * g_0_yyyyy_0_xxzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxzzz_0[i] * pb_y + g_0_yyyyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyy_0[i] = 6.0 * g_0_yyyyy_0_xyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyy_0[i] * pb_y + g_0_yyyyyy_0_xyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyz_0[i] = 6.0 * g_0_yyyyy_0_xyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyz_0[i] * pb_y + g_0_yyyyyy_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyzz_0[i] = 6.0 * g_0_yyyyy_0_xyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzz_0[i] * pb_y + g_0_yyyyyy_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyzzz_0[i] = 6.0 * g_0_yyyyy_0_xyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzz_0[i] * pb_y + g_0_yyyyyy_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xzzzz_0[i] = 6.0 * g_0_yyyyy_0_xzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzzzz_0[i] * pb_y + g_0_yyyyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyy_0[i] = 6.0 * g_0_yyyyy_0_yyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyy_0[i] * pb_y + g_0_yyyyyy_0_yyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyz_0[i] = 6.0 * g_0_yyyyy_0_yyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyz_0[i] * pb_y + g_0_yyyyyy_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyzz_0[i] = 6.0 * g_0_yyyyy_0_yyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyzz_0[i] * pb_y + g_0_yyyyyy_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyzzz_0[i] = 6.0 * g_0_yyyyy_0_yyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzzz_0[i] * pb_y + g_0_yyyyyy_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yzzzz_0[i] = 6.0 * g_0_yyyyy_0_yzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzzz_0[i] * pb_y + g_0_yyyyyy_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_zzzzz_0[i] = 6.0 * g_0_yyyyy_0_zzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzzzz_0[i] * pb_y + g_0_yyyyyy_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 609-630 components of targeted buffer : SKSH

    auto g_0_yyyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 609);

    auto g_0_yyyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 610);

    auto g_0_yyyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 611);

    auto g_0_yyyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 612);

    auto g_0_yyyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 613);

    auto g_0_yyyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 614);

    auto g_0_yyyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 615);

    auto g_0_yyyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 616);

    auto g_0_yyyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 617);

    auto g_0_yyyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 618);

    auto g_0_yyyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 619);

    auto g_0_yyyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 620);

    auto g_0_yyyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 621);

    auto g_0_yyyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 622);

    auto g_0_yyyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 623);

    auto g_0_yyyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 624);

    auto g_0_yyyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 625);

    auto g_0_yyyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 626);

    auto g_0_yyyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 627);

    auto g_0_yyyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 628);

    auto g_0_yyyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 629);

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxx_1, g_0_yyyyyy_0_xxxxx_0, g_0_yyyyyy_0_xxxxx_1, g_0_yyyyyy_0_xxxxy_0, g_0_yyyyyy_0_xxxxy_1, g_0_yyyyyy_0_xxxxz_0, g_0_yyyyyy_0_xxxxz_1, g_0_yyyyyy_0_xxxy_1, g_0_yyyyyy_0_xxxyy_0, g_0_yyyyyy_0_xxxyy_1, g_0_yyyyyy_0_xxxyz_0, g_0_yyyyyy_0_xxxyz_1, g_0_yyyyyy_0_xxxz_1, g_0_yyyyyy_0_xxxzz_0, g_0_yyyyyy_0_xxxzz_1, g_0_yyyyyy_0_xxyy_1, g_0_yyyyyy_0_xxyyy_0, g_0_yyyyyy_0_xxyyy_1, g_0_yyyyyy_0_xxyyz_0, g_0_yyyyyy_0_xxyyz_1, g_0_yyyyyy_0_xxyz_1, g_0_yyyyyy_0_xxyzz_0, g_0_yyyyyy_0_xxyzz_1, g_0_yyyyyy_0_xxzz_1, g_0_yyyyyy_0_xxzzz_0, g_0_yyyyyy_0_xxzzz_1, g_0_yyyyyy_0_xyyy_1, g_0_yyyyyy_0_xyyyy_0, g_0_yyyyyy_0_xyyyy_1, g_0_yyyyyy_0_xyyyz_0, g_0_yyyyyy_0_xyyyz_1, g_0_yyyyyy_0_xyyz_1, g_0_yyyyyy_0_xyyzz_0, g_0_yyyyyy_0_xyyzz_1, g_0_yyyyyy_0_xyzz_1, g_0_yyyyyy_0_xyzzz_0, g_0_yyyyyy_0_xyzzz_1, g_0_yyyyyy_0_xzzz_1, g_0_yyyyyy_0_xzzzz_0, g_0_yyyyyy_0_xzzzz_1, g_0_yyyyyy_0_yyyy_1, g_0_yyyyyy_0_yyyyy_0, g_0_yyyyyy_0_yyyyy_1, g_0_yyyyyy_0_yyyyz_0, g_0_yyyyyy_0_yyyyz_1, g_0_yyyyyy_0_yyyz_1, g_0_yyyyyy_0_yyyzz_0, g_0_yyyyyy_0_yyyzz_1, g_0_yyyyyy_0_yyzz_1, g_0_yyyyyy_0_yyzzz_0, g_0_yyyyyy_0_yyzzz_1, g_0_yyyyyy_0_yzzz_1, g_0_yyyyyy_0_yzzzz_0, g_0_yyyyyy_0_yzzzz_1, g_0_yyyyyy_0_zzzz_1, g_0_yyyyyy_0_zzzzz_0, g_0_yyyyyy_0_zzzzz_1, g_0_yyyyyyz_0_xxxxx_0, g_0_yyyyyyz_0_xxxxy_0, g_0_yyyyyyz_0_xxxxz_0, g_0_yyyyyyz_0_xxxyy_0, g_0_yyyyyyz_0_xxxyz_0, g_0_yyyyyyz_0_xxxzz_0, g_0_yyyyyyz_0_xxyyy_0, g_0_yyyyyyz_0_xxyyz_0, g_0_yyyyyyz_0_xxyzz_0, g_0_yyyyyyz_0_xxzzz_0, g_0_yyyyyyz_0_xyyyy_0, g_0_yyyyyyz_0_xyyyz_0, g_0_yyyyyyz_0_xyyzz_0, g_0_yyyyyyz_0_xyzzz_0, g_0_yyyyyyz_0_xzzzz_0, g_0_yyyyyyz_0_yyyyy_0, g_0_yyyyyyz_0_yyyyz_0, g_0_yyyyyyz_0_yyyzz_0, g_0_yyyyyyz_0_yyzzz_0, g_0_yyyyyyz_0_yzzzz_0, g_0_yyyyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_xxxxx_0[i] = g_0_yyyyyy_0_xxxxx_0[i] * pb_z + g_0_yyyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxy_0[i] = g_0_yyyyyy_0_xxxxy_0[i] * pb_z + g_0_yyyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxz_0[i] = g_0_yyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxz_0[i] * pb_z + g_0_yyyyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyy_0[i] = g_0_yyyyyy_0_xxxyy_0[i] * pb_z + g_0_yyyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyz_0[i] = g_0_yyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyz_0[i] * pb_z + g_0_yyyyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxzz_0[i] = 2.0 * g_0_yyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxzz_0[i] * pb_z + g_0_yyyyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyy_0[i] = g_0_yyyyyy_0_xxyyy_0[i] * pb_z + g_0_yyyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyz_0[i] = g_0_yyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyz_0[i] * pb_z + g_0_yyyyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzz_0[i] * pb_z + g_0_yyyyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzzz_0[i] * pb_z + g_0_yyyyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyy_0[i] = g_0_yyyyyy_0_xyyyy_0[i] * pb_z + g_0_yyyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyz_0[i] = g_0_yyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyz_0[i] * pb_z + g_0_yyyyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzz_0[i] * pb_z + g_0_yyyyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzz_0[i] * pb_z + g_0_yyyyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xzzzz_0[i] = 4.0 * g_0_yyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzzz_0[i] * pb_z + g_0_yyyyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyy_0[i] = g_0_yyyyyy_0_yyyyy_0[i] * pb_z + g_0_yyyyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyz_0[i] = g_0_yyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyz_0[i] * pb_z + g_0_yyyyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyzz_0[i] = 2.0 * g_0_yyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyzz_0[i] * pb_z + g_0_yyyyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyzzz_0[i] = 3.0 * g_0_yyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzzz_0[i] * pb_z + g_0_yyyyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yzzzz_0[i] = 4.0 * g_0_yyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzzz_0[i] * pb_z + g_0_yyyyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_zzzzz_0[i] = 5.0 * g_0_yyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_zzzzz_0[i] * pb_z + g_0_yyyyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 630-651 components of targeted buffer : SKSH

    auto g_0_yyyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 630);

    auto g_0_yyyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 631);

    auto g_0_yyyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 632);

    auto g_0_yyyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 633);

    auto g_0_yyyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 634);

    auto g_0_yyyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 635);

    auto g_0_yyyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 636);

    auto g_0_yyyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 637);

    auto g_0_yyyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 638);

    auto g_0_yyyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 639);

    auto g_0_yyyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 640);

    auto g_0_yyyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 641);

    auto g_0_yyyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 642);

    auto g_0_yyyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 643);

    auto g_0_yyyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 644);

    auto g_0_yyyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 645);

    auto g_0_yyyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 646);

    auto g_0_yyyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 647);

    auto g_0_yyyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 648);

    auto g_0_yyyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 649);

    auto g_0_yyyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 650);

    #pragma omp simd aligned(g_0_yyyyy_0_xxxxy_0, g_0_yyyyy_0_xxxxy_1, g_0_yyyyy_0_xxxyy_0, g_0_yyyyy_0_xxxyy_1, g_0_yyyyy_0_xxyyy_0, g_0_yyyyy_0_xxyyy_1, g_0_yyyyy_0_xyyyy_0, g_0_yyyyy_0_xyyyy_1, g_0_yyyyy_0_yyyyy_0, g_0_yyyyy_0_yyyyy_1, g_0_yyyyyz_0_xxxxy_0, g_0_yyyyyz_0_xxxxy_1, g_0_yyyyyz_0_xxxyy_0, g_0_yyyyyz_0_xxxyy_1, g_0_yyyyyz_0_xxyyy_0, g_0_yyyyyz_0_xxyyy_1, g_0_yyyyyz_0_xyyyy_0, g_0_yyyyyz_0_xyyyy_1, g_0_yyyyyz_0_yyyyy_0, g_0_yyyyyz_0_yyyyy_1, g_0_yyyyyzz_0_xxxxx_0, g_0_yyyyyzz_0_xxxxy_0, g_0_yyyyyzz_0_xxxxz_0, g_0_yyyyyzz_0_xxxyy_0, g_0_yyyyyzz_0_xxxyz_0, g_0_yyyyyzz_0_xxxzz_0, g_0_yyyyyzz_0_xxyyy_0, g_0_yyyyyzz_0_xxyyz_0, g_0_yyyyyzz_0_xxyzz_0, g_0_yyyyyzz_0_xxzzz_0, g_0_yyyyyzz_0_xyyyy_0, g_0_yyyyyzz_0_xyyyz_0, g_0_yyyyyzz_0_xyyzz_0, g_0_yyyyyzz_0_xyzzz_0, g_0_yyyyyzz_0_xzzzz_0, g_0_yyyyyzz_0_yyyyy_0, g_0_yyyyyzz_0_yyyyz_0, g_0_yyyyyzz_0_yyyzz_0, g_0_yyyyyzz_0_yyzzz_0, g_0_yyyyyzz_0_yzzzz_0, g_0_yyyyyzz_0_zzzzz_0, g_0_yyyyzz_0_xxxxx_0, g_0_yyyyzz_0_xxxxx_1, g_0_yyyyzz_0_xxxxz_0, g_0_yyyyzz_0_xxxxz_1, g_0_yyyyzz_0_xxxyz_0, g_0_yyyyzz_0_xxxyz_1, g_0_yyyyzz_0_xxxz_1, g_0_yyyyzz_0_xxxzz_0, g_0_yyyyzz_0_xxxzz_1, g_0_yyyyzz_0_xxyyz_0, g_0_yyyyzz_0_xxyyz_1, g_0_yyyyzz_0_xxyz_1, g_0_yyyyzz_0_xxyzz_0, g_0_yyyyzz_0_xxyzz_1, g_0_yyyyzz_0_xxzz_1, g_0_yyyyzz_0_xxzzz_0, g_0_yyyyzz_0_xxzzz_1, g_0_yyyyzz_0_xyyyz_0, g_0_yyyyzz_0_xyyyz_1, g_0_yyyyzz_0_xyyz_1, g_0_yyyyzz_0_xyyzz_0, g_0_yyyyzz_0_xyyzz_1, g_0_yyyyzz_0_xyzz_1, g_0_yyyyzz_0_xyzzz_0, g_0_yyyyzz_0_xyzzz_1, g_0_yyyyzz_0_xzzz_1, g_0_yyyyzz_0_xzzzz_0, g_0_yyyyzz_0_xzzzz_1, g_0_yyyyzz_0_yyyyz_0, g_0_yyyyzz_0_yyyyz_1, g_0_yyyyzz_0_yyyz_1, g_0_yyyyzz_0_yyyzz_0, g_0_yyyyzz_0_yyyzz_1, g_0_yyyyzz_0_yyzz_1, g_0_yyyyzz_0_yyzzz_0, g_0_yyyyzz_0_yyzzz_1, g_0_yyyyzz_0_yzzz_1, g_0_yyyyzz_0_yzzzz_0, g_0_yyyyzz_0_yzzzz_1, g_0_yyyyzz_0_zzzz_1, g_0_yyyyzz_0_zzzzz_0, g_0_yyyyzz_0_zzzzz_1, g_0_yyyzz_0_xxxxx_0, g_0_yyyzz_0_xxxxx_1, g_0_yyyzz_0_xxxxz_0, g_0_yyyzz_0_xxxxz_1, g_0_yyyzz_0_xxxyz_0, g_0_yyyzz_0_xxxyz_1, g_0_yyyzz_0_xxxzz_0, g_0_yyyzz_0_xxxzz_1, g_0_yyyzz_0_xxyyz_0, g_0_yyyzz_0_xxyyz_1, g_0_yyyzz_0_xxyzz_0, g_0_yyyzz_0_xxyzz_1, g_0_yyyzz_0_xxzzz_0, g_0_yyyzz_0_xxzzz_1, g_0_yyyzz_0_xyyyz_0, g_0_yyyzz_0_xyyyz_1, g_0_yyyzz_0_xyyzz_0, g_0_yyyzz_0_xyyzz_1, g_0_yyyzz_0_xyzzz_0, g_0_yyyzz_0_xyzzz_1, g_0_yyyzz_0_xzzzz_0, g_0_yyyzz_0_xzzzz_1, g_0_yyyzz_0_yyyyz_0, g_0_yyyzz_0_yyyyz_1, g_0_yyyzz_0_yyyzz_0, g_0_yyyzz_0_yyyzz_1, g_0_yyyzz_0_yyzzz_0, g_0_yyyzz_0_yyzzz_1, g_0_yyyzz_0_yzzzz_0, g_0_yyyzz_0_yzzzz_1, g_0_yyyzz_0_zzzzz_0, g_0_yyyzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_xxxxx_0[i] = 4.0 * g_0_yyyzz_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxx_0[i] * pb_y + g_0_yyyyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxy_0[i] = g_0_yyyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxxy_0[i] * pb_z + g_0_yyyyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxxz_0[i] = 4.0 * g_0_yyyzz_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxz_0[i] * pb_y + g_0_yyyyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxyy_0[i] = g_0_yyyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxyy_0[i] * pb_z + g_0_yyyyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxyz_0[i] = 4.0 * g_0_yyyzz_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyz_0[i] * pb_y + g_0_yyyyzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxzz_0[i] = 4.0 * g_0_yyyzz_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxzz_0[i] * pb_y + g_0_yyyyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyyy_0[i] = g_0_yyyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxyyy_0[i] * pb_z + g_0_yyyyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxyyz_0[i] = 4.0 * g_0_yyyzz_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyz_0[i] * pb_y + g_0_yyyyzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyzz_0[i] = 4.0 * g_0_yyyzz_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyzz_0[i] * pb_y + g_0_yyyyzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxzzz_0[i] = 4.0 * g_0_yyyzz_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxzzz_0[i] * pb_y + g_0_yyyyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyyy_0[i] = g_0_yyyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xyyyy_0[i] * pb_z + g_0_yyyyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xyyyz_0[i] = 4.0 * g_0_yyyzz_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyz_0[i] * pb_y + g_0_yyyyzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyzz_0[i] = 4.0 * g_0_yyyzz_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyzz_0[i] * pb_y + g_0_yyyyzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyzzz_0[i] = 4.0 * g_0_yyyzz_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzzz_0[i] * pb_y + g_0_yyyyzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xzzzz_0[i] = 4.0 * g_0_yyyzz_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzzzz_0[i] * pb_y + g_0_yyyyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyyy_0[i] = g_0_yyyyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_yyyyy_0[i] * pb_z + g_0_yyyyyz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_yyyyz_0[i] = 4.0 * g_0_yyyzz_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyyz_0[i] * pb_y + g_0_yyyyzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyzz_0[i] = 4.0 * g_0_yyyzz_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyzz_0[i] * pb_y + g_0_yyyyzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyzzz_0[i] = 4.0 * g_0_yyyzz_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyzzz_0[i] * pb_y + g_0_yyyyzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yzzzz_0[i] = 4.0 * g_0_yyyzz_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yzzzz_0[i] * pb_y + g_0_yyyyzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_zzzzz_0[i] = 4.0 * g_0_yyyzz_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzzzz_0[i] * pb_y + g_0_yyyyzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 651-672 components of targeted buffer : SKSH

    auto g_0_yyyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 651);

    auto g_0_yyyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 652);

    auto g_0_yyyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 653);

    auto g_0_yyyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 654);

    auto g_0_yyyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 655);

    auto g_0_yyyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 656);

    auto g_0_yyyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 657);

    auto g_0_yyyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 658);

    auto g_0_yyyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 659);

    auto g_0_yyyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 660);

    auto g_0_yyyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 661);

    auto g_0_yyyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 662);

    auto g_0_yyyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 663);

    auto g_0_yyyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 664);

    auto g_0_yyyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 665);

    auto g_0_yyyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 666);

    auto g_0_yyyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 667);

    auto g_0_yyyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 668);

    auto g_0_yyyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 669);

    auto g_0_yyyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 670);

    auto g_0_yyyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 671);

    #pragma omp simd aligned(g_0_yyyyz_0_xxxxy_0, g_0_yyyyz_0_xxxxy_1, g_0_yyyyz_0_xxxyy_0, g_0_yyyyz_0_xxxyy_1, g_0_yyyyz_0_xxyyy_0, g_0_yyyyz_0_xxyyy_1, g_0_yyyyz_0_xyyyy_0, g_0_yyyyz_0_xyyyy_1, g_0_yyyyz_0_yyyyy_0, g_0_yyyyz_0_yyyyy_1, g_0_yyyyzz_0_xxxxy_0, g_0_yyyyzz_0_xxxxy_1, g_0_yyyyzz_0_xxxyy_0, g_0_yyyyzz_0_xxxyy_1, g_0_yyyyzz_0_xxyyy_0, g_0_yyyyzz_0_xxyyy_1, g_0_yyyyzz_0_xyyyy_0, g_0_yyyyzz_0_xyyyy_1, g_0_yyyyzz_0_yyyyy_0, g_0_yyyyzz_0_yyyyy_1, g_0_yyyyzzz_0_xxxxx_0, g_0_yyyyzzz_0_xxxxy_0, g_0_yyyyzzz_0_xxxxz_0, g_0_yyyyzzz_0_xxxyy_0, g_0_yyyyzzz_0_xxxyz_0, g_0_yyyyzzz_0_xxxzz_0, g_0_yyyyzzz_0_xxyyy_0, g_0_yyyyzzz_0_xxyyz_0, g_0_yyyyzzz_0_xxyzz_0, g_0_yyyyzzz_0_xxzzz_0, g_0_yyyyzzz_0_xyyyy_0, g_0_yyyyzzz_0_xyyyz_0, g_0_yyyyzzz_0_xyyzz_0, g_0_yyyyzzz_0_xyzzz_0, g_0_yyyyzzz_0_xzzzz_0, g_0_yyyyzzz_0_yyyyy_0, g_0_yyyyzzz_0_yyyyz_0, g_0_yyyyzzz_0_yyyzz_0, g_0_yyyyzzz_0_yyzzz_0, g_0_yyyyzzz_0_yzzzz_0, g_0_yyyyzzz_0_zzzzz_0, g_0_yyyzzz_0_xxxxx_0, g_0_yyyzzz_0_xxxxx_1, g_0_yyyzzz_0_xxxxz_0, g_0_yyyzzz_0_xxxxz_1, g_0_yyyzzz_0_xxxyz_0, g_0_yyyzzz_0_xxxyz_1, g_0_yyyzzz_0_xxxz_1, g_0_yyyzzz_0_xxxzz_0, g_0_yyyzzz_0_xxxzz_1, g_0_yyyzzz_0_xxyyz_0, g_0_yyyzzz_0_xxyyz_1, g_0_yyyzzz_0_xxyz_1, g_0_yyyzzz_0_xxyzz_0, g_0_yyyzzz_0_xxyzz_1, g_0_yyyzzz_0_xxzz_1, g_0_yyyzzz_0_xxzzz_0, g_0_yyyzzz_0_xxzzz_1, g_0_yyyzzz_0_xyyyz_0, g_0_yyyzzz_0_xyyyz_1, g_0_yyyzzz_0_xyyz_1, g_0_yyyzzz_0_xyyzz_0, g_0_yyyzzz_0_xyyzz_1, g_0_yyyzzz_0_xyzz_1, g_0_yyyzzz_0_xyzzz_0, g_0_yyyzzz_0_xyzzz_1, g_0_yyyzzz_0_xzzz_1, g_0_yyyzzz_0_xzzzz_0, g_0_yyyzzz_0_xzzzz_1, g_0_yyyzzz_0_yyyyz_0, g_0_yyyzzz_0_yyyyz_1, g_0_yyyzzz_0_yyyz_1, g_0_yyyzzz_0_yyyzz_0, g_0_yyyzzz_0_yyyzz_1, g_0_yyyzzz_0_yyzz_1, g_0_yyyzzz_0_yyzzz_0, g_0_yyyzzz_0_yyzzz_1, g_0_yyyzzz_0_yzzz_1, g_0_yyyzzz_0_yzzzz_0, g_0_yyyzzz_0_yzzzz_1, g_0_yyyzzz_0_zzzz_1, g_0_yyyzzz_0_zzzzz_0, g_0_yyyzzz_0_zzzzz_1, g_0_yyzzz_0_xxxxx_0, g_0_yyzzz_0_xxxxx_1, g_0_yyzzz_0_xxxxz_0, g_0_yyzzz_0_xxxxz_1, g_0_yyzzz_0_xxxyz_0, g_0_yyzzz_0_xxxyz_1, g_0_yyzzz_0_xxxzz_0, g_0_yyzzz_0_xxxzz_1, g_0_yyzzz_0_xxyyz_0, g_0_yyzzz_0_xxyyz_1, g_0_yyzzz_0_xxyzz_0, g_0_yyzzz_0_xxyzz_1, g_0_yyzzz_0_xxzzz_0, g_0_yyzzz_0_xxzzz_1, g_0_yyzzz_0_xyyyz_0, g_0_yyzzz_0_xyyyz_1, g_0_yyzzz_0_xyyzz_0, g_0_yyzzz_0_xyyzz_1, g_0_yyzzz_0_xyzzz_0, g_0_yyzzz_0_xyzzz_1, g_0_yyzzz_0_xzzzz_0, g_0_yyzzz_0_xzzzz_1, g_0_yyzzz_0_yyyyz_0, g_0_yyzzz_0_yyyyz_1, g_0_yyzzz_0_yyyzz_0, g_0_yyzzz_0_yyyzz_1, g_0_yyzzz_0_yyzzz_0, g_0_yyzzz_0_yyzzz_1, g_0_yyzzz_0_yzzzz_0, g_0_yyzzz_0_yzzzz_1, g_0_yyzzz_0_zzzzz_0, g_0_yyzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_xxxxx_0[i] = 3.0 * g_0_yyzzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxx_0[i] * pb_y + g_0_yyyzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxy_0[i] = 2.0 * g_0_yyyyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxy_0[i] * pb_z + g_0_yyyyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxxz_0[i] = 3.0 * g_0_yyzzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxz_0[i] * pb_y + g_0_yyyzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxyy_0[i] = 2.0 * g_0_yyyyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxyy_0[i] * pb_z + g_0_yyyyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxyz_0[i] = 3.0 * g_0_yyzzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyz_0[i] * pb_y + g_0_yyyzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxzz_0[i] = 3.0 * g_0_yyzzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxzz_0[i] * pb_y + g_0_yyyzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyyyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxyyy_0[i] * pb_z + g_0_yyyyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxyyz_0[i] = 3.0 * g_0_yyzzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyz_0[i] * pb_y + g_0_yyyzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyzz_0[i] = 3.0 * g_0_yyzzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyzz_0[i] * pb_y + g_0_yyyzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxzzz_0[i] = 3.0 * g_0_yyzzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxzzz_0[i] * pb_y + g_0_yyyzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyyy_0[i] = 2.0 * g_0_yyyyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xyyyy_0[i] * pb_z + g_0_yyyyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xyyyz_0[i] = 3.0 * g_0_yyzzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyz_0[i] * pb_y + g_0_yyyzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyzz_0[i] = 3.0 * g_0_yyzzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyzz_0[i] * pb_y + g_0_yyyzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyzzz_0[i] = 3.0 * g_0_yyzzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzzz_0[i] * pb_y + g_0_yyyzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xzzzz_0[i] = 3.0 * g_0_yyzzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzzzz_0[i] * pb_y + g_0_yyyzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyyy_0[i] = 2.0 * g_0_yyyyz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_yyyyy_0[i] * pb_z + g_0_yyyyzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_yyyyz_0[i] = 3.0 * g_0_yyzzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyyz_0[i] * pb_y + g_0_yyyzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyzz_0[i] = 3.0 * g_0_yyzzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyzz_0[i] * pb_y + g_0_yyyzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyzzz_0[i] = 3.0 * g_0_yyzzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyzzz_0[i] * pb_y + g_0_yyyzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yzzzz_0[i] = 3.0 * g_0_yyzzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yzzzz_0[i] * pb_y + g_0_yyyzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_zzzzz_0[i] = 3.0 * g_0_yyzzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzzzz_0[i] * pb_y + g_0_yyyzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 672-693 components of targeted buffer : SKSH

    auto g_0_yyyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 672);

    auto g_0_yyyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 673);

    auto g_0_yyyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 674);

    auto g_0_yyyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 675);

    auto g_0_yyyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 676);

    auto g_0_yyyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 677);

    auto g_0_yyyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 678);

    auto g_0_yyyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 679);

    auto g_0_yyyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 680);

    auto g_0_yyyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 681);

    auto g_0_yyyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 682);

    auto g_0_yyyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 683);

    auto g_0_yyyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 684);

    auto g_0_yyyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 685);

    auto g_0_yyyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 686);

    auto g_0_yyyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 687);

    auto g_0_yyyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 688);

    auto g_0_yyyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 689);

    auto g_0_yyyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 690);

    auto g_0_yyyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 691);

    auto g_0_yyyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 692);

    #pragma omp simd aligned(g_0_yyyzz_0_xxxxy_0, g_0_yyyzz_0_xxxxy_1, g_0_yyyzz_0_xxxyy_0, g_0_yyyzz_0_xxxyy_1, g_0_yyyzz_0_xxyyy_0, g_0_yyyzz_0_xxyyy_1, g_0_yyyzz_0_xyyyy_0, g_0_yyyzz_0_xyyyy_1, g_0_yyyzz_0_yyyyy_0, g_0_yyyzz_0_yyyyy_1, g_0_yyyzzz_0_xxxxy_0, g_0_yyyzzz_0_xxxxy_1, g_0_yyyzzz_0_xxxyy_0, g_0_yyyzzz_0_xxxyy_1, g_0_yyyzzz_0_xxyyy_0, g_0_yyyzzz_0_xxyyy_1, g_0_yyyzzz_0_xyyyy_0, g_0_yyyzzz_0_xyyyy_1, g_0_yyyzzz_0_yyyyy_0, g_0_yyyzzz_0_yyyyy_1, g_0_yyyzzzz_0_xxxxx_0, g_0_yyyzzzz_0_xxxxy_0, g_0_yyyzzzz_0_xxxxz_0, g_0_yyyzzzz_0_xxxyy_0, g_0_yyyzzzz_0_xxxyz_0, g_0_yyyzzzz_0_xxxzz_0, g_0_yyyzzzz_0_xxyyy_0, g_0_yyyzzzz_0_xxyyz_0, g_0_yyyzzzz_0_xxyzz_0, g_0_yyyzzzz_0_xxzzz_0, g_0_yyyzzzz_0_xyyyy_0, g_0_yyyzzzz_0_xyyyz_0, g_0_yyyzzzz_0_xyyzz_0, g_0_yyyzzzz_0_xyzzz_0, g_0_yyyzzzz_0_xzzzz_0, g_0_yyyzzzz_0_yyyyy_0, g_0_yyyzzzz_0_yyyyz_0, g_0_yyyzzzz_0_yyyzz_0, g_0_yyyzzzz_0_yyzzz_0, g_0_yyyzzzz_0_yzzzz_0, g_0_yyyzzzz_0_zzzzz_0, g_0_yyzzzz_0_xxxxx_0, g_0_yyzzzz_0_xxxxx_1, g_0_yyzzzz_0_xxxxz_0, g_0_yyzzzz_0_xxxxz_1, g_0_yyzzzz_0_xxxyz_0, g_0_yyzzzz_0_xxxyz_1, g_0_yyzzzz_0_xxxz_1, g_0_yyzzzz_0_xxxzz_0, g_0_yyzzzz_0_xxxzz_1, g_0_yyzzzz_0_xxyyz_0, g_0_yyzzzz_0_xxyyz_1, g_0_yyzzzz_0_xxyz_1, g_0_yyzzzz_0_xxyzz_0, g_0_yyzzzz_0_xxyzz_1, g_0_yyzzzz_0_xxzz_1, g_0_yyzzzz_0_xxzzz_0, g_0_yyzzzz_0_xxzzz_1, g_0_yyzzzz_0_xyyyz_0, g_0_yyzzzz_0_xyyyz_1, g_0_yyzzzz_0_xyyz_1, g_0_yyzzzz_0_xyyzz_0, g_0_yyzzzz_0_xyyzz_1, g_0_yyzzzz_0_xyzz_1, g_0_yyzzzz_0_xyzzz_0, g_0_yyzzzz_0_xyzzz_1, g_0_yyzzzz_0_xzzz_1, g_0_yyzzzz_0_xzzzz_0, g_0_yyzzzz_0_xzzzz_1, g_0_yyzzzz_0_yyyyz_0, g_0_yyzzzz_0_yyyyz_1, g_0_yyzzzz_0_yyyz_1, g_0_yyzzzz_0_yyyzz_0, g_0_yyzzzz_0_yyyzz_1, g_0_yyzzzz_0_yyzz_1, g_0_yyzzzz_0_yyzzz_0, g_0_yyzzzz_0_yyzzz_1, g_0_yyzzzz_0_yzzz_1, g_0_yyzzzz_0_yzzzz_0, g_0_yyzzzz_0_yzzzz_1, g_0_yyzzzz_0_zzzz_1, g_0_yyzzzz_0_zzzzz_0, g_0_yyzzzz_0_zzzzz_1, g_0_yzzzz_0_xxxxx_0, g_0_yzzzz_0_xxxxx_1, g_0_yzzzz_0_xxxxz_0, g_0_yzzzz_0_xxxxz_1, g_0_yzzzz_0_xxxyz_0, g_0_yzzzz_0_xxxyz_1, g_0_yzzzz_0_xxxzz_0, g_0_yzzzz_0_xxxzz_1, g_0_yzzzz_0_xxyyz_0, g_0_yzzzz_0_xxyyz_1, g_0_yzzzz_0_xxyzz_0, g_0_yzzzz_0_xxyzz_1, g_0_yzzzz_0_xxzzz_0, g_0_yzzzz_0_xxzzz_1, g_0_yzzzz_0_xyyyz_0, g_0_yzzzz_0_xyyyz_1, g_0_yzzzz_0_xyyzz_0, g_0_yzzzz_0_xyyzz_1, g_0_yzzzz_0_xyzzz_0, g_0_yzzzz_0_xyzzz_1, g_0_yzzzz_0_xzzzz_0, g_0_yzzzz_0_xzzzz_1, g_0_yzzzz_0_yyyyz_0, g_0_yzzzz_0_yyyyz_1, g_0_yzzzz_0_yyyzz_0, g_0_yzzzz_0_yyyzz_1, g_0_yzzzz_0_yyzzz_0, g_0_yzzzz_0_yyzzz_1, g_0_yzzzz_0_yzzzz_0, g_0_yzzzz_0_yzzzz_1, g_0_yzzzz_0_zzzzz_0, g_0_yzzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_xxxxx_0[i] = 2.0 * g_0_yzzzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxx_0[i] * pb_y + g_0_yyzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxy_0[i] = 3.0 * g_0_yyyzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxy_0[i] * pb_z + g_0_yyyzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxxz_0[i] = 2.0 * g_0_yzzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxz_0[i] * pb_y + g_0_yyzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxyy_0[i] = 3.0 * g_0_yyyzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxyy_0[i] * pb_z + g_0_yyyzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxyz_0[i] = 2.0 * g_0_yzzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyz_0[i] * pb_y + g_0_yyzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxzz_0[i] = 2.0 * g_0_yzzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxzz_0[i] * pb_y + g_0_yyzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyyy_0[i] = 3.0 * g_0_yyyzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxyyy_0[i] * pb_z + g_0_yyyzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyz_0[i] * pb_y + g_0_yyzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyzz_0[i] * pb_y + g_0_yyzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxzzz_0[i] = 2.0 * g_0_yzzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxzzz_0[i] * pb_y + g_0_yyzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyyy_0[i] = 3.0 * g_0_yyyzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xyyyy_0[i] * pb_z + g_0_yyyzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xyyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyz_0[i] * pb_y + g_0_yyzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyzz_0[i] = 2.0 * g_0_yzzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyzz_0[i] * pb_y + g_0_yyzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyzzz_0[i] = 2.0 * g_0_yzzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzzz_0[i] * pb_y + g_0_yyzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xzzzz_0[i] = 2.0 * g_0_yzzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzzzz_0[i] * pb_y + g_0_yyzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyyy_0[i] = 3.0 * g_0_yyyzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_yyyyy_0[i] * pb_z + g_0_yyyzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_yyyyz_0[i] = 2.0 * g_0_yzzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyyz_0[i] * pb_y + g_0_yyzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyzz_0[i] = 2.0 * g_0_yzzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyzz_0[i] * pb_y + g_0_yyzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyzzz_0[i] = 2.0 * g_0_yzzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyzzz_0[i] * pb_y + g_0_yyzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yzzzz_0[i] = 2.0 * g_0_yzzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yzzzz_0[i] * pb_y + g_0_yyzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_zzzzz_0[i] = 2.0 * g_0_yzzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzzzz_0[i] * pb_y + g_0_yyzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 693-714 components of targeted buffer : SKSH

    auto g_0_yyzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 693);

    auto g_0_yyzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 694);

    auto g_0_yyzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 695);

    auto g_0_yyzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 696);

    auto g_0_yyzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 697);

    auto g_0_yyzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 698);

    auto g_0_yyzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 699);

    auto g_0_yyzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 700);

    auto g_0_yyzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 701);

    auto g_0_yyzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 702);

    auto g_0_yyzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 703);

    auto g_0_yyzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 704);

    auto g_0_yyzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 705);

    auto g_0_yyzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 706);

    auto g_0_yyzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 707);

    auto g_0_yyzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 708);

    auto g_0_yyzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 709);

    auto g_0_yyzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 710);

    auto g_0_yyzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 711);

    auto g_0_yyzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 712);

    auto g_0_yyzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 713);

    #pragma omp simd aligned(g_0_yyzzz_0_xxxxy_0, g_0_yyzzz_0_xxxxy_1, g_0_yyzzz_0_xxxyy_0, g_0_yyzzz_0_xxxyy_1, g_0_yyzzz_0_xxyyy_0, g_0_yyzzz_0_xxyyy_1, g_0_yyzzz_0_xyyyy_0, g_0_yyzzz_0_xyyyy_1, g_0_yyzzz_0_yyyyy_0, g_0_yyzzz_0_yyyyy_1, g_0_yyzzzz_0_xxxxy_0, g_0_yyzzzz_0_xxxxy_1, g_0_yyzzzz_0_xxxyy_0, g_0_yyzzzz_0_xxxyy_1, g_0_yyzzzz_0_xxyyy_0, g_0_yyzzzz_0_xxyyy_1, g_0_yyzzzz_0_xyyyy_0, g_0_yyzzzz_0_xyyyy_1, g_0_yyzzzz_0_yyyyy_0, g_0_yyzzzz_0_yyyyy_1, g_0_yyzzzzz_0_xxxxx_0, g_0_yyzzzzz_0_xxxxy_0, g_0_yyzzzzz_0_xxxxz_0, g_0_yyzzzzz_0_xxxyy_0, g_0_yyzzzzz_0_xxxyz_0, g_0_yyzzzzz_0_xxxzz_0, g_0_yyzzzzz_0_xxyyy_0, g_0_yyzzzzz_0_xxyyz_0, g_0_yyzzzzz_0_xxyzz_0, g_0_yyzzzzz_0_xxzzz_0, g_0_yyzzzzz_0_xyyyy_0, g_0_yyzzzzz_0_xyyyz_0, g_0_yyzzzzz_0_xyyzz_0, g_0_yyzzzzz_0_xyzzz_0, g_0_yyzzzzz_0_xzzzz_0, g_0_yyzzzzz_0_yyyyy_0, g_0_yyzzzzz_0_yyyyz_0, g_0_yyzzzzz_0_yyyzz_0, g_0_yyzzzzz_0_yyzzz_0, g_0_yyzzzzz_0_yzzzz_0, g_0_yyzzzzz_0_zzzzz_0, g_0_yzzzzz_0_xxxxx_0, g_0_yzzzzz_0_xxxxx_1, g_0_yzzzzz_0_xxxxz_0, g_0_yzzzzz_0_xxxxz_1, g_0_yzzzzz_0_xxxyz_0, g_0_yzzzzz_0_xxxyz_1, g_0_yzzzzz_0_xxxz_1, g_0_yzzzzz_0_xxxzz_0, g_0_yzzzzz_0_xxxzz_1, g_0_yzzzzz_0_xxyyz_0, g_0_yzzzzz_0_xxyyz_1, g_0_yzzzzz_0_xxyz_1, g_0_yzzzzz_0_xxyzz_0, g_0_yzzzzz_0_xxyzz_1, g_0_yzzzzz_0_xxzz_1, g_0_yzzzzz_0_xxzzz_0, g_0_yzzzzz_0_xxzzz_1, g_0_yzzzzz_0_xyyyz_0, g_0_yzzzzz_0_xyyyz_1, g_0_yzzzzz_0_xyyz_1, g_0_yzzzzz_0_xyyzz_0, g_0_yzzzzz_0_xyyzz_1, g_0_yzzzzz_0_xyzz_1, g_0_yzzzzz_0_xyzzz_0, g_0_yzzzzz_0_xyzzz_1, g_0_yzzzzz_0_xzzz_1, g_0_yzzzzz_0_xzzzz_0, g_0_yzzzzz_0_xzzzz_1, g_0_yzzzzz_0_yyyyz_0, g_0_yzzzzz_0_yyyyz_1, g_0_yzzzzz_0_yyyz_1, g_0_yzzzzz_0_yyyzz_0, g_0_yzzzzz_0_yyyzz_1, g_0_yzzzzz_0_yyzz_1, g_0_yzzzzz_0_yyzzz_0, g_0_yzzzzz_0_yyzzz_1, g_0_yzzzzz_0_yzzz_1, g_0_yzzzzz_0_yzzzz_0, g_0_yzzzzz_0_yzzzz_1, g_0_yzzzzz_0_zzzz_1, g_0_yzzzzz_0_zzzzz_0, g_0_yzzzzz_0_zzzzz_1, g_0_zzzzz_0_xxxxx_0, g_0_zzzzz_0_xxxxx_1, g_0_zzzzz_0_xxxxz_0, g_0_zzzzz_0_xxxxz_1, g_0_zzzzz_0_xxxyz_0, g_0_zzzzz_0_xxxyz_1, g_0_zzzzz_0_xxxzz_0, g_0_zzzzz_0_xxxzz_1, g_0_zzzzz_0_xxyyz_0, g_0_zzzzz_0_xxyyz_1, g_0_zzzzz_0_xxyzz_0, g_0_zzzzz_0_xxyzz_1, g_0_zzzzz_0_xxzzz_0, g_0_zzzzz_0_xxzzz_1, g_0_zzzzz_0_xyyyz_0, g_0_zzzzz_0_xyyyz_1, g_0_zzzzz_0_xyyzz_0, g_0_zzzzz_0_xyyzz_1, g_0_zzzzz_0_xyzzz_0, g_0_zzzzz_0_xyzzz_1, g_0_zzzzz_0_xzzzz_0, g_0_zzzzz_0_xzzzz_1, g_0_zzzzz_0_yyyyz_0, g_0_zzzzz_0_yyyyz_1, g_0_zzzzz_0_yyyzz_0, g_0_zzzzz_0_yyyzz_1, g_0_zzzzz_0_yyzzz_0, g_0_zzzzz_0_yyzzz_1, g_0_zzzzz_0_yzzzz_0, g_0_zzzzz_0_yzzzz_1, g_0_zzzzz_0_zzzzz_0, g_0_zzzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_xxxxx_0[i] = g_0_zzzzz_0_xxxxx_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxx_0[i] * pb_y + g_0_yzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxy_0[i] = 4.0 * g_0_yyzzz_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxy_0[i] * pb_z + g_0_yyzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxxz_0[i] = g_0_zzzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxz_0[i] * pb_y + g_0_yzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxyy_0[i] = 4.0 * g_0_yyzzz_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxyy_0[i] * pb_z + g_0_yyzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxyz_0[i] = g_0_zzzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyz_0[i] * pb_y + g_0_yzzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxzz_0[i] = g_0_zzzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxzz_0[i] * pb_y + g_0_yzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyyy_0[i] = 4.0 * g_0_yyzzz_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxyyy_0[i] * pb_z + g_0_yyzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxyyz_0[i] = g_0_zzzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyz_0[i] * pb_y + g_0_yzzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyzz_0[i] = g_0_zzzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyzz_0[i] * pb_y + g_0_yzzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxzzz_0[i] = g_0_zzzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzzz_0[i] * pb_y + g_0_yzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyyy_0[i] = 4.0 * g_0_yyzzz_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xyyyy_0[i] * pb_z + g_0_yyzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xyyyz_0[i] = g_0_zzzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyz_0[i] * pb_y + g_0_yzzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyzz_0[i] = g_0_zzzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyzz_0[i] * pb_y + g_0_yzzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyzzz_0[i] = g_0_zzzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyzzz_0[i] * pb_y + g_0_yzzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xzzzz_0[i] = g_0_zzzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzzz_0[i] * pb_y + g_0_yzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyyy_0[i] = 4.0 * g_0_yyzzz_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_yyyyy_0[i] * pb_z + g_0_yyzzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_yyyyz_0[i] = g_0_zzzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyyz_0[i] * pb_y + g_0_yzzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyzz_0[i] = g_0_zzzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyzz_0[i] * pb_y + g_0_yzzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyzzz_0[i] = g_0_zzzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyzzz_0[i] * pb_y + g_0_yzzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yzzzz_0[i] = g_0_zzzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yzzzz_0[i] * pb_y + g_0_yzzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_zzzzz_0[i] = g_0_zzzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzzz_0[i] * pb_y + g_0_yzzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 714-735 components of targeted buffer : SKSH

    auto g_0_yzzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 714);

    auto g_0_yzzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 715);

    auto g_0_yzzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 716);

    auto g_0_yzzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 717);

    auto g_0_yzzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 718);

    auto g_0_yzzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 719);

    auto g_0_yzzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 720);

    auto g_0_yzzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 721);

    auto g_0_yzzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 722);

    auto g_0_yzzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 723);

    auto g_0_yzzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 724);

    auto g_0_yzzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 725);

    auto g_0_yzzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 726);

    auto g_0_yzzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 727);

    auto g_0_yzzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 728);

    auto g_0_yzzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 729);

    auto g_0_yzzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 730);

    auto g_0_yzzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 731);

    auto g_0_yzzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 732);

    auto g_0_yzzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 733);

    auto g_0_yzzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 734);

    #pragma omp simd aligned(g_0_yzzzzzz_0_xxxxx_0, g_0_yzzzzzz_0_xxxxy_0, g_0_yzzzzzz_0_xxxxz_0, g_0_yzzzzzz_0_xxxyy_0, g_0_yzzzzzz_0_xxxyz_0, g_0_yzzzzzz_0_xxxzz_0, g_0_yzzzzzz_0_xxyyy_0, g_0_yzzzzzz_0_xxyyz_0, g_0_yzzzzzz_0_xxyzz_0, g_0_yzzzzzz_0_xxzzz_0, g_0_yzzzzzz_0_xyyyy_0, g_0_yzzzzzz_0_xyyyz_0, g_0_yzzzzzz_0_xyyzz_0, g_0_yzzzzzz_0_xyzzz_0, g_0_yzzzzzz_0_xzzzz_0, g_0_yzzzzzz_0_yyyyy_0, g_0_yzzzzzz_0_yyyyz_0, g_0_yzzzzzz_0_yyyzz_0, g_0_yzzzzzz_0_yyzzz_0, g_0_yzzzzzz_0_yzzzz_0, g_0_yzzzzzz_0_zzzzz_0, g_0_zzzzzz_0_xxxx_1, g_0_zzzzzz_0_xxxxx_0, g_0_zzzzzz_0_xxxxx_1, g_0_zzzzzz_0_xxxxy_0, g_0_zzzzzz_0_xxxxy_1, g_0_zzzzzz_0_xxxxz_0, g_0_zzzzzz_0_xxxxz_1, g_0_zzzzzz_0_xxxy_1, g_0_zzzzzz_0_xxxyy_0, g_0_zzzzzz_0_xxxyy_1, g_0_zzzzzz_0_xxxyz_0, g_0_zzzzzz_0_xxxyz_1, g_0_zzzzzz_0_xxxz_1, g_0_zzzzzz_0_xxxzz_0, g_0_zzzzzz_0_xxxzz_1, g_0_zzzzzz_0_xxyy_1, g_0_zzzzzz_0_xxyyy_0, g_0_zzzzzz_0_xxyyy_1, g_0_zzzzzz_0_xxyyz_0, g_0_zzzzzz_0_xxyyz_1, g_0_zzzzzz_0_xxyz_1, g_0_zzzzzz_0_xxyzz_0, g_0_zzzzzz_0_xxyzz_1, g_0_zzzzzz_0_xxzz_1, g_0_zzzzzz_0_xxzzz_0, g_0_zzzzzz_0_xxzzz_1, g_0_zzzzzz_0_xyyy_1, g_0_zzzzzz_0_xyyyy_0, g_0_zzzzzz_0_xyyyy_1, g_0_zzzzzz_0_xyyyz_0, g_0_zzzzzz_0_xyyyz_1, g_0_zzzzzz_0_xyyz_1, g_0_zzzzzz_0_xyyzz_0, g_0_zzzzzz_0_xyyzz_1, g_0_zzzzzz_0_xyzz_1, g_0_zzzzzz_0_xyzzz_0, g_0_zzzzzz_0_xyzzz_1, g_0_zzzzzz_0_xzzz_1, g_0_zzzzzz_0_xzzzz_0, g_0_zzzzzz_0_xzzzz_1, g_0_zzzzzz_0_yyyy_1, g_0_zzzzzz_0_yyyyy_0, g_0_zzzzzz_0_yyyyy_1, g_0_zzzzzz_0_yyyyz_0, g_0_zzzzzz_0_yyyyz_1, g_0_zzzzzz_0_yyyz_1, g_0_zzzzzz_0_yyyzz_0, g_0_zzzzzz_0_yyyzz_1, g_0_zzzzzz_0_yyzz_1, g_0_zzzzzz_0_yyzzz_0, g_0_zzzzzz_0_yyzzz_1, g_0_zzzzzz_0_yzzz_1, g_0_zzzzzz_0_yzzzz_0, g_0_zzzzzz_0_yzzzz_1, g_0_zzzzzz_0_zzzz_1, g_0_zzzzzz_0_zzzzz_0, g_0_zzzzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_xxxxx_0[i] = g_0_zzzzzz_0_xxxxx_0[i] * pb_y + g_0_zzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxy_0[i] = g_0_zzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxy_0[i] * pb_y + g_0_zzzzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxz_0[i] = g_0_zzzzzz_0_xxxxz_0[i] * pb_y + g_0_zzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyy_0[i] = 2.0 * g_0_zzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyy_0[i] * pb_y + g_0_zzzzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyz_0[i] = g_0_zzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyz_0[i] * pb_y + g_0_zzzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxzz_0[i] = g_0_zzzzzz_0_xxxzz_0[i] * pb_y + g_0_zzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyy_0[i] = 3.0 * g_0_zzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyy_0[i] * pb_y + g_0_zzzzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyz_0[i] * pb_y + g_0_zzzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyzz_0[i] = g_0_zzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzz_0[i] * pb_y + g_0_zzzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxzzz_0[i] = g_0_zzzzzz_0_xxzzz_0[i] * pb_y + g_0_zzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyy_0[i] = 4.0 * g_0_zzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyy_0[i] * pb_y + g_0_zzzzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyz_0[i] = 3.0 * g_0_zzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyz_0[i] * pb_y + g_0_zzzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyzz_0[i] = 2.0 * g_0_zzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzz_0[i] * pb_y + g_0_zzzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyzzz_0[i] = g_0_zzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzz_0[i] * pb_y + g_0_zzzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xzzzz_0[i] = g_0_zzzzzz_0_xzzzz_0[i] * pb_y + g_0_zzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyy_0[i] = 5.0 * g_0_zzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyy_0[i] * pb_y + g_0_zzzzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyz_0[i] = 4.0 * g_0_zzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyz_0[i] * pb_y + g_0_zzzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyzz_0[i] = 3.0 * g_0_zzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyzz_0[i] * pb_y + g_0_zzzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyzzz_0[i] = 2.0 * g_0_zzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzzz_0[i] * pb_y + g_0_zzzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yzzzz_0[i] = g_0_zzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzzz_0[i] * pb_y + g_0_zzzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_zzzzz_0[i] = g_0_zzzzzz_0_zzzzz_0[i] * pb_y + g_0_zzzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 735-756 components of targeted buffer : SKSH

    auto g_0_zzzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sksh + 735);

    auto g_0_zzzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sksh + 736);

    auto g_0_zzzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sksh + 737);

    auto g_0_zzzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sksh + 738);

    auto g_0_zzzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sksh + 739);

    auto g_0_zzzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sksh + 740);

    auto g_0_zzzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sksh + 741);

    auto g_0_zzzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sksh + 742);

    auto g_0_zzzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sksh + 743);

    auto g_0_zzzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sksh + 744);

    auto g_0_zzzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sksh + 745);

    auto g_0_zzzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sksh + 746);

    auto g_0_zzzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sksh + 747);

    auto g_0_zzzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sksh + 748);

    auto g_0_zzzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sksh + 749);

    auto g_0_zzzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sksh + 750);

    auto g_0_zzzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sksh + 751);

    auto g_0_zzzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sksh + 752);

    auto g_0_zzzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sksh + 753);

    auto g_0_zzzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sksh + 754);

    auto g_0_zzzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sksh + 755);

    #pragma omp simd aligned(g_0_zzzzz_0_xxxxx_0, g_0_zzzzz_0_xxxxx_1, g_0_zzzzz_0_xxxxy_0, g_0_zzzzz_0_xxxxy_1, g_0_zzzzz_0_xxxxz_0, g_0_zzzzz_0_xxxxz_1, g_0_zzzzz_0_xxxyy_0, g_0_zzzzz_0_xxxyy_1, g_0_zzzzz_0_xxxyz_0, g_0_zzzzz_0_xxxyz_1, g_0_zzzzz_0_xxxzz_0, g_0_zzzzz_0_xxxzz_1, g_0_zzzzz_0_xxyyy_0, g_0_zzzzz_0_xxyyy_1, g_0_zzzzz_0_xxyyz_0, g_0_zzzzz_0_xxyyz_1, g_0_zzzzz_0_xxyzz_0, g_0_zzzzz_0_xxyzz_1, g_0_zzzzz_0_xxzzz_0, g_0_zzzzz_0_xxzzz_1, g_0_zzzzz_0_xyyyy_0, g_0_zzzzz_0_xyyyy_1, g_0_zzzzz_0_xyyyz_0, g_0_zzzzz_0_xyyyz_1, g_0_zzzzz_0_xyyzz_0, g_0_zzzzz_0_xyyzz_1, g_0_zzzzz_0_xyzzz_0, g_0_zzzzz_0_xyzzz_1, g_0_zzzzz_0_xzzzz_0, g_0_zzzzz_0_xzzzz_1, g_0_zzzzz_0_yyyyy_0, g_0_zzzzz_0_yyyyy_1, g_0_zzzzz_0_yyyyz_0, g_0_zzzzz_0_yyyyz_1, g_0_zzzzz_0_yyyzz_0, g_0_zzzzz_0_yyyzz_1, g_0_zzzzz_0_yyzzz_0, g_0_zzzzz_0_yyzzz_1, g_0_zzzzz_0_yzzzz_0, g_0_zzzzz_0_yzzzz_1, g_0_zzzzz_0_zzzzz_0, g_0_zzzzz_0_zzzzz_1, g_0_zzzzzz_0_xxxx_1, g_0_zzzzzz_0_xxxxx_0, g_0_zzzzzz_0_xxxxx_1, g_0_zzzzzz_0_xxxxy_0, g_0_zzzzzz_0_xxxxy_1, g_0_zzzzzz_0_xxxxz_0, g_0_zzzzzz_0_xxxxz_1, g_0_zzzzzz_0_xxxy_1, g_0_zzzzzz_0_xxxyy_0, g_0_zzzzzz_0_xxxyy_1, g_0_zzzzzz_0_xxxyz_0, g_0_zzzzzz_0_xxxyz_1, g_0_zzzzzz_0_xxxz_1, g_0_zzzzzz_0_xxxzz_0, g_0_zzzzzz_0_xxxzz_1, g_0_zzzzzz_0_xxyy_1, g_0_zzzzzz_0_xxyyy_0, g_0_zzzzzz_0_xxyyy_1, g_0_zzzzzz_0_xxyyz_0, g_0_zzzzzz_0_xxyyz_1, g_0_zzzzzz_0_xxyz_1, g_0_zzzzzz_0_xxyzz_0, g_0_zzzzzz_0_xxyzz_1, g_0_zzzzzz_0_xxzz_1, g_0_zzzzzz_0_xxzzz_0, g_0_zzzzzz_0_xxzzz_1, g_0_zzzzzz_0_xyyy_1, g_0_zzzzzz_0_xyyyy_0, g_0_zzzzzz_0_xyyyy_1, g_0_zzzzzz_0_xyyyz_0, g_0_zzzzzz_0_xyyyz_1, g_0_zzzzzz_0_xyyz_1, g_0_zzzzzz_0_xyyzz_0, g_0_zzzzzz_0_xyyzz_1, g_0_zzzzzz_0_xyzz_1, g_0_zzzzzz_0_xyzzz_0, g_0_zzzzzz_0_xyzzz_1, g_0_zzzzzz_0_xzzz_1, g_0_zzzzzz_0_xzzzz_0, g_0_zzzzzz_0_xzzzz_1, g_0_zzzzzz_0_yyyy_1, g_0_zzzzzz_0_yyyyy_0, g_0_zzzzzz_0_yyyyy_1, g_0_zzzzzz_0_yyyyz_0, g_0_zzzzzz_0_yyyyz_1, g_0_zzzzzz_0_yyyz_1, g_0_zzzzzz_0_yyyzz_0, g_0_zzzzzz_0_yyyzz_1, g_0_zzzzzz_0_yyzz_1, g_0_zzzzzz_0_yyzzz_0, g_0_zzzzzz_0_yyzzz_1, g_0_zzzzzz_0_yzzz_1, g_0_zzzzzz_0_yzzzz_0, g_0_zzzzzz_0_yzzzz_1, g_0_zzzzzz_0_zzzz_1, g_0_zzzzzz_0_zzzzz_0, g_0_zzzzzz_0_zzzzz_1, g_0_zzzzzzz_0_xxxxx_0, g_0_zzzzzzz_0_xxxxy_0, g_0_zzzzzzz_0_xxxxz_0, g_0_zzzzzzz_0_xxxyy_0, g_0_zzzzzzz_0_xxxyz_0, g_0_zzzzzzz_0_xxxzz_0, g_0_zzzzzzz_0_xxyyy_0, g_0_zzzzzzz_0_xxyyz_0, g_0_zzzzzzz_0_xxyzz_0, g_0_zzzzzzz_0_xxzzz_0, g_0_zzzzzzz_0_xyyyy_0, g_0_zzzzzzz_0_xyyyz_0, g_0_zzzzzzz_0_xyyzz_0, g_0_zzzzzzz_0_xyzzz_0, g_0_zzzzzzz_0_xzzzz_0, g_0_zzzzzzz_0_yyyyy_0, g_0_zzzzzzz_0_yyyyz_0, g_0_zzzzzzz_0_yyyzz_0, g_0_zzzzzzz_0_yyzzz_0, g_0_zzzzzzz_0_yzzzz_0, g_0_zzzzzzz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_xxxxx_0[i] = 6.0 * g_0_zzzzz_0_xxxxx_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxx_0[i] * pb_z + g_0_zzzzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxy_0[i] = 6.0 * g_0_zzzzz_0_xxxxy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxy_0[i] * pb_z + g_0_zzzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxz_0[i] = 6.0 * g_0_zzzzz_0_xxxxz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxz_0[i] * pb_z + g_0_zzzzzz_0_xxxxz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyy_0[i] = 6.0 * g_0_zzzzz_0_xxxyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxyy_0[i] * pb_z + g_0_zzzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyz_0[i] = 6.0 * g_0_zzzzz_0_xxxyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyz_0[i] * pb_z + g_0_zzzzzz_0_xxxyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxzz_0[i] = 6.0 * g_0_zzzzz_0_xxxzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxzz_0[i] * pb_z + g_0_zzzzzz_0_xxxzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyy_0[i] = 6.0 * g_0_zzzzz_0_xxyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxyyy_0[i] * pb_z + g_0_zzzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyz_0[i] = 6.0 * g_0_zzzzz_0_xxyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyz_0[i] * pb_z + g_0_zzzzzz_0_xxyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyzz_0[i] = 6.0 * g_0_zzzzz_0_xxyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzz_0[i] * pb_z + g_0_zzzzzz_0_xxyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxzzz_0[i] = 6.0 * g_0_zzzzz_0_xxzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzzz_0[i] * pb_z + g_0_zzzzzz_0_xxzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyy_0[i] = 6.0 * g_0_zzzzz_0_xyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyyyy_0[i] * pb_z + g_0_zzzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyz_0[i] = 6.0 * g_0_zzzzz_0_xyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyz_0[i] * pb_z + g_0_zzzzzz_0_xyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyzz_0[i] = 6.0 * g_0_zzzzz_0_xyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzz_0[i] * pb_z + g_0_zzzzzz_0_xyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyzzz_0[i] = 6.0 * g_0_zzzzz_0_xyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzz_0[i] * pb_z + g_0_zzzzzz_0_xyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xzzzz_0[i] = 6.0 * g_0_zzzzz_0_xzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzzz_0[i] * pb_z + g_0_zzzzzz_0_xzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyy_0[i] = 6.0 * g_0_zzzzz_0_yyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyyyy_0[i] * pb_z + g_0_zzzzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyz_0[i] = 6.0 * g_0_zzzzz_0_yyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyz_0[i] * pb_z + g_0_zzzzzz_0_yyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyzz_0[i] = 6.0 * g_0_zzzzz_0_yyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyzz_0[i] * pb_z + g_0_zzzzzz_0_yyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyzzz_0[i] = 6.0 * g_0_zzzzz_0_yyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzzz_0[i] * pb_z + g_0_zzzzzz_0_yyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yzzzz_0[i] = 6.0 * g_0_zzzzz_0_yzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzzz_0[i] * pb_z + g_0_zzzzzz_0_yzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_zzzzz_0[i] = 6.0 * g_0_zzzzz_0_zzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_zzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_zzzzz_0[i] * pb_z + g_0_zzzzzz_0_zzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

