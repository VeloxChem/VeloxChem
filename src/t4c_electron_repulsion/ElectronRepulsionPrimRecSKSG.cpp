#include "ElectronRepulsionPrimRecSKSG.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sksg(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sksg,
                                  size_t                idx_eri_0_shsg,
                                  size_t                idx_eri_1_shsg,
                                  size_t                idx_eri_1_sisf,
                                  size_t                idx_eri_0_sisg,
                                  size_t                idx_eri_1_sisg,
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

    /// Set up components of auxilary buffer : SHSG

    auto g_0_xxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg);

    auto g_0_xxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 1);

    auto g_0_xxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 2);

    auto g_0_xxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 3);

    auto g_0_xxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 4);

    auto g_0_xxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 5);

    auto g_0_xxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 6);

    auto g_0_xxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 7);

    auto g_0_xxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 8);

    auto g_0_xxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 9);

    auto g_0_xxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 10);

    auto g_0_xxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 11);

    auto g_0_xxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 12);

    auto g_0_xxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 13);

    auto g_0_xxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 14);

    auto g_0_xxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 15);

    auto g_0_xxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 17);

    auto g_0_xxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 20);

    auto g_0_xxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 24);

    auto g_0_xxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 30);

    auto g_0_xxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 31);

    auto g_0_xxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 33);

    auto g_0_xxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 36);

    auto g_0_xxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 45);

    auto g_0_xxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 46);

    auto g_0_xxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 47);

    auto g_0_xxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 48);

    auto g_0_xxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 49);

    auto g_0_xxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 50);

    auto g_0_xxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 51);

    auto g_0_xxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 52);

    auto g_0_xxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 53);

    auto g_0_xxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 54);

    auto g_0_xxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 55);

    auto g_0_xxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 56);

    auto g_0_xxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 57);

    auto g_0_xxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 58);

    auto g_0_xxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 59);

    auto g_0_xxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 75);

    auto g_0_xxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 76);

    auto g_0_xxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 77);

    auto g_0_xxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 78);

    auto g_0_xxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 79);

    auto g_0_xxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 80);

    auto g_0_xxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 81);

    auto g_0_xxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 82);

    auto g_0_xxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 83);

    auto g_0_xxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 84);

    auto g_0_xxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 85);

    auto g_0_xxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 86);

    auto g_0_xxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 87);

    auto g_0_xxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 88);

    auto g_0_xxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 89);

    auto g_0_xxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 90);

    auto g_0_xxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 91);

    auto g_0_xxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 92);

    auto g_0_xxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 93);

    auto g_0_xxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 94);

    auto g_0_xxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 95);

    auto g_0_xxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 96);

    auto g_0_xxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 97);

    auto g_0_xxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 98);

    auto g_0_xxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 99);

    auto g_0_xxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 100);

    auto g_0_xxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 101);

    auto g_0_xxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 102);

    auto g_0_xxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 103);

    auto g_0_xxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 104);

    auto g_0_xxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 106);

    auto g_0_xxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 108);

    auto g_0_xxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 111);

    auto g_0_xxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 120);

    auto g_0_xxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 122);

    auto g_0_xxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 125);

    auto g_0_xxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 129);

    auto g_0_xxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 135);

    auto g_0_xxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 136);

    auto g_0_xxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 137);

    auto g_0_xxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 138);

    auto g_0_xxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 139);

    auto g_0_xxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 140);

    auto g_0_xxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 141);

    auto g_0_xxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 142);

    auto g_0_xxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 143);

    auto g_0_xxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 144);

    auto g_0_xxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 145);

    auto g_0_xxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 146);

    auto g_0_xxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 147);

    auto g_0_xxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 148);

    auto g_0_xxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 149);

    auto g_0_xyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 151);

    auto g_0_xyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 153);

    auto g_0_xyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 154);

    auto g_0_xyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 156);

    auto g_0_xyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 157);

    auto g_0_xyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 158);

    auto g_0_xyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 160);

    auto g_0_xyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 161);

    auto g_0_xyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 162);

    auto g_0_xyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 163);

    auto g_0_xyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 164);

    auto g_0_xyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 184);

    auto g_0_xyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 187);

    auto g_0_xyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 188);

    auto g_0_xyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 190);

    auto g_0_xyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 191);

    auto g_0_xyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 192);

    auto g_0_xyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 193);

    auto g_0_xyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 194);

    auto g_0_xzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 212);

    auto g_0_xzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 214);

    auto g_0_xzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 215);

    auto g_0_xzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 217);

    auto g_0_xzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 218);

    auto g_0_xzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 219);

    auto g_0_xzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 220);

    auto g_0_xzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 221);

    auto g_0_xzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 222);

    auto g_0_xzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 223);

    auto g_0_xzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 224);

    auto g_0_yyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 225);

    auto g_0_yyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 226);

    auto g_0_yyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 227);

    auto g_0_yyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 228);

    auto g_0_yyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 229);

    auto g_0_yyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 230);

    auto g_0_yyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 231);

    auto g_0_yyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 232);

    auto g_0_yyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 233);

    auto g_0_yyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 234);

    auto g_0_yyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 235);

    auto g_0_yyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 236);

    auto g_0_yyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 237);

    auto g_0_yyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 238);

    auto g_0_yyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 239);

    auto g_0_yyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 241);

    auto g_0_yyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 243);

    auto g_0_yyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 246);

    auto g_0_yyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 250);

    auto g_0_yyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 255);

    auto g_0_yyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 256);

    auto g_0_yyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 257);

    auto g_0_yyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 258);

    auto g_0_yyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 259);

    auto g_0_yyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 260);

    auto g_0_yyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 261);

    auto g_0_yyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 262);

    auto g_0_yyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 263);

    auto g_0_yyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 264);

    auto g_0_yyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 265);

    auto g_0_yyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 266);

    auto g_0_yyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 267);

    auto g_0_yyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 268);

    auto g_0_yyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 269);

    auto g_0_yyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 270);

    auto g_0_yyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 271);

    auto g_0_yyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 272);

    auto g_0_yyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 273);

    auto g_0_yyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 274);

    auto g_0_yyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 275);

    auto g_0_yyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 276);

    auto g_0_yyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 277);

    auto g_0_yyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 278);

    auto g_0_yyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 279);

    auto g_0_yyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 280);

    auto g_0_yyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 281);

    auto g_0_yyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 282);

    auto g_0_yyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 283);

    auto g_0_yyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 284);

    auto g_0_yzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 285);

    auto g_0_yzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 287);

    auto g_0_yzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 289);

    auto g_0_yzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 290);

    auto g_0_yzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 292);

    auto g_0_yzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 293);

    auto g_0_yzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 294);

    auto g_0_yzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 296);

    auto g_0_yzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 297);

    auto g_0_yzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 298);

    auto g_0_yzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 299);

    auto g_0_zzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 300);

    auto g_0_zzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 301);

    auto g_0_zzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 302);

    auto g_0_zzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 303);

    auto g_0_zzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 304);

    auto g_0_zzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 305);

    auto g_0_zzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 306);

    auto g_0_zzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 307);

    auto g_0_zzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 308);

    auto g_0_zzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 309);

    auto g_0_zzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 310);

    auto g_0_zzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 311);

    auto g_0_zzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 312);

    auto g_0_zzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 313);

    auto g_0_zzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 314);

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

    auto g_0_xxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 15);

    auto g_0_xxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 17);

    auto g_0_xxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 20);

    auto g_0_xxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 24);

    auto g_0_xxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 30);

    auto g_0_xxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 31);

    auto g_0_xxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 33);

    auto g_0_xxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 36);

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

    auto g_0_xxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 106);

    auto g_0_xxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 108);

    auto g_0_xxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 111);

    auto g_0_xxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 120);

    auto g_0_xxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 122);

    auto g_0_xxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 125);

    auto g_0_xxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 129);

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

    auto g_0_xyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 164);

    auto g_0_xyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 184);

    auto g_0_xyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 187);

    auto g_0_xyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 188);

    auto g_0_xyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 190);

    auto g_0_xyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 191);

    auto g_0_xyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 192);

    auto g_0_xyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 193);

    auto g_0_xyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 194);

    auto g_0_xzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 212);

    auto g_0_xzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 214);

    auto g_0_xzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 215);

    auto g_0_xzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 217);

    auto g_0_xzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 218);

    auto g_0_xzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 219);

    auto g_0_xzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 220);

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

    auto g_0_yyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 241);

    auto g_0_yyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 243);

    auto g_0_yyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 246);

    auto g_0_yyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 250);

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

    auto g_0_yzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 285);

    auto g_0_yzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 287);

    auto g_0_yzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 289);

    auto g_0_yzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 290);

    auto g_0_yzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 292);

    auto g_0_yzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 293);

    auto g_0_yzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 294);

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

    /// Set up components of auxilary buffer : SISF

    auto g_0_xxxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_sisf);

    auto g_0_xxxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 1);

    auto g_0_xxxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 2);

    auto g_0_xxxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 3);

    auto g_0_xxxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 4);

    auto g_0_xxxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 5);

    auto g_0_xxxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 6);

    auto g_0_xxxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 7);

    auto g_0_xxxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 8);

    auto g_0_xxxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 9);

    auto g_0_xxxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 22);

    auto g_0_xxxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 24);

    auto g_0_xxxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 25);

    auto g_0_xxxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 27);

    auto g_0_xxxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 28);

    auto g_0_xxxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 29);

    auto g_0_xxxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 30);

    auto g_0_xxxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 31);

    auto g_0_xxxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 32);

    auto g_0_xxxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 33);

    auto g_0_xxxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 34);

    auto g_0_xxxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 35);

    auto g_0_xxxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 36);

    auto g_0_xxxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 37);

    auto g_0_xxxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 38);

    auto g_0_xxxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 39);

    auto g_0_xxxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 50);

    auto g_0_xxxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 51);

    auto g_0_xxxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 52);

    auto g_0_xxxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 53);

    auto g_0_xxxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 54);

    auto g_0_xxxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 55);

    auto g_0_xxxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 56);

    auto g_0_xxxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 57);

    auto g_0_xxxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 58);

    auto g_0_xxxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 59);

    auto g_0_xxxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 60);

    auto g_0_xxxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 61);

    auto g_0_xxxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 62);

    auto g_0_xxxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 63);

    auto g_0_xxxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 64);

    auto g_0_xxxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 65);

    auto g_0_xxxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 66);

    auto g_0_xxxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 67);

    auto g_0_xxxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 68);

    auto g_0_xxxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 69);

    auto g_0_xxxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 90);

    auto g_0_xxxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 91);

    auto g_0_xxxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 92);

    auto g_0_xxxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 93);

    auto g_0_xxxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 94);

    auto g_0_xxxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 95);

    auto g_0_xxxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 96);

    auto g_0_xxxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 97);

    auto g_0_xxxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 98);

    auto g_0_xxxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 99);

    auto g_0_xxyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 100);

    auto g_0_xxyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 101);

    auto g_0_xxyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 102);

    auto g_0_xxyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 103);

    auto g_0_xxyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 104);

    auto g_0_xxyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 105);

    auto g_0_xxyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 106);

    auto g_0_xxyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 107);

    auto g_0_xxyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 108);

    auto g_0_xxyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 109);

    auto g_0_xxyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 124);

    auto g_0_xxyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 127);

    auto g_0_xxyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 128);

    auto g_0_xxzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 140);

    auto g_0_xxzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 141);

    auto g_0_xxzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 142);

    auto g_0_xxzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 143);

    auto g_0_xxzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 144);

    auto g_0_xxzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 145);

    auto g_0_xxzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 146);

    auto g_0_xxzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 147);

    auto g_0_xxzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 148);

    auto g_0_xxzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 149);

    auto g_0_xyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 151);

    auto g_0_xyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 153);

    auto g_0_xyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 154);

    auto g_0_xyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 156);

    auto g_0_xyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 157);

    auto g_0_xyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 158);

    auto g_0_xyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 174);

    auto g_0_xyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 177);

    auto g_0_xyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 178);

    auto g_0_xyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 184);

    auto g_0_xyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 187);

    auto g_0_xyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 188);

    auto g_0_xzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 202);

    auto g_0_xzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 204);

    auto g_0_xzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 205);

    auto g_0_xzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 207);

    auto g_0_xzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 208);

    auto g_0_xzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 209);

    auto g_0_yyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 210);

    auto g_0_yyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 211);

    auto g_0_yyyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 212);

    auto g_0_yyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 213);

    auto g_0_yyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 214);

    auto g_0_yyyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 215);

    auto g_0_yyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 216);

    auto g_0_yyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 217);

    auto g_0_yyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 218);

    auto g_0_yyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 219);

    auto g_0_yyyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 222);

    auto g_0_yyyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 224);

    auto g_0_yyyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 225);

    auto g_0_yyyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 227);

    auto g_0_yyyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 228);

    auto g_0_yyyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 229);

    auto g_0_yyyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 230);

    auto g_0_yyyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 231);

    auto g_0_yyyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 232);

    auto g_0_yyyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 233);

    auto g_0_yyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 234);

    auto g_0_yyyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 235);

    auto g_0_yyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 236);

    auto g_0_yyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 237);

    auto g_0_yyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 238);

    auto g_0_yyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 239);

    auto g_0_yyyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 240);

    auto g_0_yyyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 241);

    auto g_0_yyyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 242);

    auto g_0_yyyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 243);

    auto g_0_yyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 244);

    auto g_0_yyyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 245);

    auto g_0_yyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 246);

    auto g_0_yyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 247);

    auto g_0_yyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 248);

    auto g_0_yyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 249);

    auto g_0_yyzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 250);

    auto g_0_yyzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 251);

    auto g_0_yyzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 252);

    auto g_0_yyzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 253);

    auto g_0_yyzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 254);

    auto g_0_yyzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 255);

    auto g_0_yyzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 256);

    auto g_0_yyzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 257);

    auto g_0_yyzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 258);

    auto g_0_yyzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 259);

    auto g_0_yzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 261);

    auto g_0_yzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 262);

    auto g_0_yzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 263);

    auto g_0_yzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 264);

    auto g_0_yzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 265);

    auto g_0_yzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 266);

    auto g_0_yzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 267);

    auto g_0_yzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 268);

    auto g_0_yzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 269);

    auto g_0_zzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 270);

    auto g_0_zzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 271);

    auto g_0_zzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 272);

    auto g_0_zzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 273);

    auto g_0_zzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 274);

    auto g_0_zzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 275);

    auto g_0_zzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 276);

    auto g_0_zzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 277);

    auto g_0_zzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 278);

    auto g_0_zzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 279);

    /// Set up components of auxilary buffer : SISG

    auto g_0_xxxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg);

    auto g_0_xxxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 1);

    auto g_0_xxxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 2);

    auto g_0_xxxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 3);

    auto g_0_xxxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 4);

    auto g_0_xxxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 5);

    auto g_0_xxxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 6);

    auto g_0_xxxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 7);

    auto g_0_xxxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 8);

    auto g_0_xxxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 9);

    auto g_0_xxxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 10);

    auto g_0_xxxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 11);

    auto g_0_xxxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 12);

    auto g_0_xxxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 13);

    auto g_0_xxxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 14);

    auto g_0_xxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 15);

    auto g_0_xxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 16);

    auto g_0_xxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 17);

    auto g_0_xxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 18);

    auto g_0_xxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 20);

    auto g_0_xxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 21);

    auto g_0_xxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 24);

    auto g_0_xxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 25);

    auto g_0_xxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 30);

    auto g_0_xxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 31);

    auto g_0_xxxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 32);

    auto g_0_xxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 33);

    auto g_0_xxxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 34);

    auto g_0_xxxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 35);

    auto g_0_xxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 36);

    auto g_0_xxxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 37);

    auto g_0_xxxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 38);

    auto g_0_xxxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 39);

    auto g_0_xxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 41);

    auto g_0_xxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 42);

    auto g_0_xxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 43);

    auto g_0_xxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 44);

    auto g_0_xxxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 45);

    auto g_0_xxxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 46);

    auto g_0_xxxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 47);

    auto g_0_xxxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 48);

    auto g_0_xxxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 49);

    auto g_0_xxxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 50);

    auto g_0_xxxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 51);

    auto g_0_xxxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 52);

    auto g_0_xxxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 53);

    auto g_0_xxxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 54);

    auto g_0_xxxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 55);

    auto g_0_xxxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 56);

    auto g_0_xxxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 57);

    auto g_0_xxxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 58);

    auto g_0_xxxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 59);

    auto g_0_xxxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 75);

    auto g_0_xxxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 76);

    auto g_0_xxxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 77);

    auto g_0_xxxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 78);

    auto g_0_xxxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 79);

    auto g_0_xxxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 80);

    auto g_0_xxxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 81);

    auto g_0_xxxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 82);

    auto g_0_xxxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 83);

    auto g_0_xxxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 84);

    auto g_0_xxxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 85);

    auto g_0_xxxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 86);

    auto g_0_xxxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 87);

    auto g_0_xxxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 88);

    auto g_0_xxxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 89);

    auto g_0_xxxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 90);

    auto g_0_xxxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 91);

    auto g_0_xxxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 92);

    auto g_0_xxxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 93);

    auto g_0_xxxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 94);

    auto g_0_xxxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 95);

    auto g_0_xxxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 96);

    auto g_0_xxxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 97);

    auto g_0_xxxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 98);

    auto g_0_xxxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 99);

    auto g_0_xxxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 100);

    auto g_0_xxxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 101);

    auto g_0_xxxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 102);

    auto g_0_xxxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 103);

    auto g_0_xxxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 104);

    auto g_0_xxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 106);

    auto g_0_xxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 108);

    auto g_0_xxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 111);

    auto g_0_xxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 120);

    auto g_0_xxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 122);

    auto g_0_xxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 125);

    auto g_0_xxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 129);

    auto g_0_xxxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 135);

    auto g_0_xxxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 136);

    auto g_0_xxxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 137);

    auto g_0_xxxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 138);

    auto g_0_xxxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 139);

    auto g_0_xxxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 140);

    auto g_0_xxxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 141);

    auto g_0_xxxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 142);

    auto g_0_xxxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 143);

    auto g_0_xxxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 144);

    auto g_0_xxxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 145);

    auto g_0_xxxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 146);

    auto g_0_xxxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 147);

    auto g_0_xxxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 148);

    auto g_0_xxxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 149);

    auto g_0_xxyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 150);

    auto g_0_xxyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 151);

    auto g_0_xxyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 152);

    auto g_0_xxyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 153);

    auto g_0_xxyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 154);

    auto g_0_xxyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 155);

    auto g_0_xxyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 156);

    auto g_0_xxyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 157);

    auto g_0_xxyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 158);

    auto g_0_xxyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 159);

    auto g_0_xxyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 160);

    auto g_0_xxyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 161);

    auto g_0_xxyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 162);

    auto g_0_xxyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 163);

    auto g_0_xxyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 164);

    auto g_0_xxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 166);

    auto g_0_xxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 168);

    auto g_0_xxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 171);

    auto g_0_xxyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 180);

    auto g_0_xxyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 181);

    auto g_0_xxyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 182);

    auto g_0_xxyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 183);

    auto g_0_xxyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 184);

    auto g_0_xxyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 185);

    auto g_0_xxyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 186);

    auto g_0_xxyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 187);

    auto g_0_xxyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 188);

    auto g_0_xxyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 189);

    auto g_0_xxyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 190);

    auto g_0_xxyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 191);

    auto g_0_xxyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 192);

    auto g_0_xxyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 193);

    auto g_0_xxyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 194);

    auto g_0_xxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 195);

    auto g_0_xxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 197);

    auto g_0_xxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 200);

    auto g_0_xxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 204);

    auto g_0_xxzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 210);

    auto g_0_xxzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 211);

    auto g_0_xxzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 212);

    auto g_0_xxzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 213);

    auto g_0_xxzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 214);

    auto g_0_xxzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 215);

    auto g_0_xxzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 216);

    auto g_0_xxzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 217);

    auto g_0_xxzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 218);

    auto g_0_xxzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 219);

    auto g_0_xxzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 220);

    auto g_0_xxzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 221);

    auto g_0_xxzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 222);

    auto g_0_xxzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 223);

    auto g_0_xxzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 224);

    auto g_0_xyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 225);

    auto g_0_xyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 226);

    auto g_0_xyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 228);

    auto g_0_xyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 229);

    auto g_0_xyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 231);

    auto g_0_xyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 232);

    auto g_0_xyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 233);

    auto g_0_xyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 235);

    auto g_0_xyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 236);

    auto g_0_xyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 237);

    auto g_0_xyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 238);

    auto g_0_xyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 239);

    auto g_0_xyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 259);

    auto g_0_xyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 262);

    auto g_0_xyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 263);

    auto g_0_xyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 265);

    auto g_0_xyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 266);

    auto g_0_xyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 267);

    auto g_0_xyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 268);

    auto g_0_xyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 269);

    auto g_0_xyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 274);

    auto g_0_xyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 277);

    auto g_0_xyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 278);

    auto g_0_xyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 280);

    auto g_0_xyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 281);

    auto g_0_xyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 282);

    auto g_0_xyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 283);

    auto g_0_xyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 284);

    auto g_0_xzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 300);

    auto g_0_xzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 302);

    auto g_0_xzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 304);

    auto g_0_xzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 305);

    auto g_0_xzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 307);

    auto g_0_xzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 308);

    auto g_0_xzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 309);

    auto g_0_xzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 310);

    auto g_0_xzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 311);

    auto g_0_xzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 312);

    auto g_0_xzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 313);

    auto g_0_xzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 314);

    auto g_0_yyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 315);

    auto g_0_yyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 316);

    auto g_0_yyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 317);

    auto g_0_yyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 318);

    auto g_0_yyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 319);

    auto g_0_yyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 320);

    auto g_0_yyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 321);

    auto g_0_yyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 322);

    auto g_0_yyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 323);

    auto g_0_yyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 324);

    auto g_0_yyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 325);

    auto g_0_yyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 326);

    auto g_0_yyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 327);

    auto g_0_yyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 328);

    auto g_0_yyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 329);

    auto g_0_yyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 331);

    auto g_0_yyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 332);

    auto g_0_yyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 333);

    auto g_0_yyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 334);

    auto g_0_yyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 335);

    auto g_0_yyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 336);

    auto g_0_yyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 337);

    auto g_0_yyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 338);

    auto g_0_yyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 339);

    auto g_0_yyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 340);

    auto g_0_yyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 341);

    auto g_0_yyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 342);

    auto g_0_yyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 343);

    auto g_0_yyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 344);

    auto g_0_yyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 345);

    auto g_0_yyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 346);

    auto g_0_yyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 347);

    auto g_0_yyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 348);

    auto g_0_yyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 349);

    auto g_0_yyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 350);

    auto g_0_yyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 351);

    auto g_0_yyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 352);

    auto g_0_yyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 353);

    auto g_0_yyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 354);

    auto g_0_yyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 355);

    auto g_0_yyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 356);

    auto g_0_yyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 357);

    auto g_0_yyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 358);

    auto g_0_yyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 359);

    auto g_0_yyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 360);

    auto g_0_yyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 361);

    auto g_0_yyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 362);

    auto g_0_yyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 363);

    auto g_0_yyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 364);

    auto g_0_yyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 365);

    auto g_0_yyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 366);

    auto g_0_yyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 367);

    auto g_0_yyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 368);

    auto g_0_yyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 369);

    auto g_0_yyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 370);

    auto g_0_yyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 371);

    auto g_0_yyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 372);

    auto g_0_yyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 373);

    auto g_0_yyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 374);

    auto g_0_yyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 375);

    auto g_0_yyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 376);

    auto g_0_yyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 377);

    auto g_0_yyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 378);

    auto g_0_yyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 379);

    auto g_0_yyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 380);

    auto g_0_yyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 381);

    auto g_0_yyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 382);

    auto g_0_yyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 383);

    auto g_0_yyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 384);

    auto g_0_yyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 385);

    auto g_0_yyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 386);

    auto g_0_yyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 387);

    auto g_0_yyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 388);

    auto g_0_yyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 389);

    auto g_0_yzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 390);

    auto g_0_yzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 391);

    auto g_0_yzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 392);

    auto g_0_yzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 393);

    auto g_0_yzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 394);

    auto g_0_yzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 395);

    auto g_0_yzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 396);

    auto g_0_yzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 397);

    auto g_0_yzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 398);

    auto g_0_yzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 399);

    auto g_0_yzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 400);

    auto g_0_yzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 401);

    auto g_0_yzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 402);

    auto g_0_yzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 403);

    auto g_0_yzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 404);

    auto g_0_zzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 405);

    auto g_0_zzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 406);

    auto g_0_zzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 407);

    auto g_0_zzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 408);

    auto g_0_zzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 409);

    auto g_0_zzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 410);

    auto g_0_zzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 411);

    auto g_0_zzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 412);

    auto g_0_zzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 413);

    auto g_0_zzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 414);

    auto g_0_zzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 415);

    auto g_0_zzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 416);

    auto g_0_zzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 417);

    auto g_0_zzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 418);

    auto g_0_zzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 419);

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

    auto g_0_xxxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 15);

    auto g_0_xxxxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 16);

    auto g_0_xxxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 17);

    auto g_0_xxxxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 18);

    auto g_0_xxxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 20);

    auto g_0_xxxxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 21);

    auto g_0_xxxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 24);

    auto g_0_xxxxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 25);

    auto g_0_xxxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 30);

    auto g_0_xxxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 31);

    auto g_0_xxxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 32);

    auto g_0_xxxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 33);

    auto g_0_xxxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 34);

    auto g_0_xxxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 35);

    auto g_0_xxxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 36);

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

    auto g_0_xxxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 106);

    auto g_0_xxxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 108);

    auto g_0_xxxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 111);

    auto g_0_xxxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 120);

    auto g_0_xxxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 122);

    auto g_0_xxxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 125);

    auto g_0_xxxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 129);

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

    auto g_0_xxyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 166);

    auto g_0_xxyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 168);

    auto g_0_xxyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 171);

    auto g_0_xxyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 180);

    auto g_0_xxyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 181);

    auto g_0_xxyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 182);

    auto g_0_xxyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 183);

    auto g_0_xxyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 184);

    auto g_0_xxyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 185);

    auto g_0_xxyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 186);

    auto g_0_xxyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 187);

    auto g_0_xxyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 188);

    auto g_0_xxyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 189);

    auto g_0_xxyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 190);

    auto g_0_xxyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 191);

    auto g_0_xxyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 192);

    auto g_0_xxyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 193);

    auto g_0_xxyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 194);

    auto g_0_xxyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 195);

    auto g_0_xxyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 197);

    auto g_0_xxyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 200);

    auto g_0_xxyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 204);

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

    auto g_0_xyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 225);

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

    auto g_0_xyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 239);

    auto g_0_xyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 259);

    auto g_0_xyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 262);

    auto g_0_xyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 263);

    auto g_0_xyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 265);

    auto g_0_xyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 266);

    auto g_0_xyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 267);

    auto g_0_xyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 268);

    auto g_0_xyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 269);

    auto g_0_xyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 274);

    auto g_0_xyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 277);

    auto g_0_xyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 278);

    auto g_0_xyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 280);

    auto g_0_xyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sisg + 281);

    auto g_0_xyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sisg + 282);

    auto g_0_xyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sisg + 283);

    auto g_0_xyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sisg + 284);

    auto g_0_xzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 300);

    auto g_0_xzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 302);

    auto g_0_xzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 304);

    auto g_0_xzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 305);

    auto g_0_xzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 307);

    auto g_0_xzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 308);

    auto g_0_xzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 309);

    auto g_0_xzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 310);

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

    auto g_0_yyyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 331);

    auto g_0_yyyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 332);

    auto g_0_yyyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 333);

    auto g_0_yyyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 334);

    auto g_0_yyyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 335);

    auto g_0_yyyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 336);

    auto g_0_yyyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 337);

    auto g_0_yyyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 338);

    auto g_0_yyyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 339);

    auto g_0_yyyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 340);

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

    auto g_0_yzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 390);

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

    /// Set up 0-15 components of targeted buffer : SKSG

    auto g_0_xxxxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg);

    auto g_0_xxxxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 1);

    auto g_0_xxxxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 2);

    auto g_0_xxxxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 3);

    auto g_0_xxxxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 4);

    auto g_0_xxxxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 5);

    auto g_0_xxxxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 6);

    auto g_0_xxxxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 7);

    auto g_0_xxxxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 8);

    auto g_0_xxxxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 9);

    auto g_0_xxxxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 10);

    auto g_0_xxxxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 11);

    auto g_0_xxxxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 12);

    auto g_0_xxxxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 13);

    auto g_0_xxxxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 14);

#pragma omp simd aligned(g_0_xxxxx_0_xxxx_0,       \
                             g_0_xxxxx_0_xxxx_1,   \
                             g_0_xxxxx_0_xxxy_0,   \
                             g_0_xxxxx_0_xxxy_1,   \
                             g_0_xxxxx_0_xxxz_0,   \
                             g_0_xxxxx_0_xxxz_1,   \
                             g_0_xxxxx_0_xxyy_0,   \
                             g_0_xxxxx_0_xxyy_1,   \
                             g_0_xxxxx_0_xxyz_0,   \
                             g_0_xxxxx_0_xxyz_1,   \
                             g_0_xxxxx_0_xxzz_0,   \
                             g_0_xxxxx_0_xxzz_1,   \
                             g_0_xxxxx_0_xyyy_0,   \
                             g_0_xxxxx_0_xyyy_1,   \
                             g_0_xxxxx_0_xyyz_0,   \
                             g_0_xxxxx_0_xyyz_1,   \
                             g_0_xxxxx_0_xyzz_0,   \
                             g_0_xxxxx_0_xyzz_1,   \
                             g_0_xxxxx_0_xzzz_0,   \
                             g_0_xxxxx_0_xzzz_1,   \
                             g_0_xxxxx_0_yyyy_0,   \
                             g_0_xxxxx_0_yyyy_1,   \
                             g_0_xxxxx_0_yyyz_0,   \
                             g_0_xxxxx_0_yyyz_1,   \
                             g_0_xxxxx_0_yyzz_0,   \
                             g_0_xxxxx_0_yyzz_1,   \
                             g_0_xxxxx_0_yzzz_0,   \
                             g_0_xxxxx_0_yzzz_1,   \
                             g_0_xxxxx_0_zzzz_0,   \
                             g_0_xxxxx_0_zzzz_1,   \
                             g_0_xxxxxx_0_xxx_1,   \
                             g_0_xxxxxx_0_xxxx_0,  \
                             g_0_xxxxxx_0_xxxx_1,  \
                             g_0_xxxxxx_0_xxxy_0,  \
                             g_0_xxxxxx_0_xxxy_1,  \
                             g_0_xxxxxx_0_xxxz_0,  \
                             g_0_xxxxxx_0_xxxz_1,  \
                             g_0_xxxxxx_0_xxy_1,   \
                             g_0_xxxxxx_0_xxyy_0,  \
                             g_0_xxxxxx_0_xxyy_1,  \
                             g_0_xxxxxx_0_xxyz_0,  \
                             g_0_xxxxxx_0_xxyz_1,  \
                             g_0_xxxxxx_0_xxz_1,   \
                             g_0_xxxxxx_0_xxzz_0,  \
                             g_0_xxxxxx_0_xxzz_1,  \
                             g_0_xxxxxx_0_xyy_1,   \
                             g_0_xxxxxx_0_xyyy_0,  \
                             g_0_xxxxxx_0_xyyy_1,  \
                             g_0_xxxxxx_0_xyyz_0,  \
                             g_0_xxxxxx_0_xyyz_1,  \
                             g_0_xxxxxx_0_xyz_1,   \
                             g_0_xxxxxx_0_xyzz_0,  \
                             g_0_xxxxxx_0_xyzz_1,  \
                             g_0_xxxxxx_0_xzz_1,   \
                             g_0_xxxxxx_0_xzzz_0,  \
                             g_0_xxxxxx_0_xzzz_1,  \
                             g_0_xxxxxx_0_yyy_1,   \
                             g_0_xxxxxx_0_yyyy_0,  \
                             g_0_xxxxxx_0_yyyy_1,  \
                             g_0_xxxxxx_0_yyyz_0,  \
                             g_0_xxxxxx_0_yyyz_1,  \
                             g_0_xxxxxx_0_yyz_1,   \
                             g_0_xxxxxx_0_yyzz_0,  \
                             g_0_xxxxxx_0_yyzz_1,  \
                             g_0_xxxxxx_0_yzz_1,   \
                             g_0_xxxxxx_0_yzzz_0,  \
                             g_0_xxxxxx_0_yzzz_1,  \
                             g_0_xxxxxx_0_zzz_1,   \
                             g_0_xxxxxx_0_zzzz_0,  \
                             g_0_xxxxxx_0_zzzz_1,  \
                             g_0_xxxxxxx_0_xxxx_0, \
                             g_0_xxxxxxx_0_xxxy_0, \
                             g_0_xxxxxxx_0_xxxz_0, \
                             g_0_xxxxxxx_0_xxyy_0, \
                             g_0_xxxxxxx_0_xxyz_0, \
                             g_0_xxxxxxx_0_xxzz_0, \
                             g_0_xxxxxxx_0_xyyy_0, \
                             g_0_xxxxxxx_0_xyyz_0, \
                             g_0_xxxxxxx_0_xyzz_0, \
                             g_0_xxxxxxx_0_xzzz_0, \
                             g_0_xxxxxxx_0_yyyy_0, \
                             g_0_xxxxxxx_0_yyyz_0, \
                             g_0_xxxxxxx_0_yyzz_0, \
                             g_0_xxxxxxx_0_yzzz_0, \
                             g_0_xxxxxxx_0_zzzz_0, \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_xxxx_0[i] = 6.0 * g_0_xxxxx_0_xxxx_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxx_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxx_0[i] * pb_x + g_0_xxxxxx_0_xxxx_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxy_0[i] = 6.0 * g_0_xxxxx_0_xxxy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxy_0[i] * pb_x + g_0_xxxxxx_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxz_0[i] = 6.0 * g_0_xxxxx_0_xxxz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxz_0[i] * pb_x + g_0_xxxxxx_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyy_0[i] = 6.0 * g_0_xxxxx_0_xxyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyy_0[i] * pb_x + g_0_xxxxxx_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyz_0[i] = 6.0 * g_0_xxxxx_0_xxyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyz_0[i] * pb_x + g_0_xxxxxx_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxzz_0[i] = 6.0 * g_0_xxxxx_0_xxzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzz_0[i] * pb_x + g_0_xxxxxx_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyy_0[i] = 6.0 * g_0_xxxxx_0_xyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyy_1[i] * fi_abcd_0 +
                                  g_0_xxxxxx_0_xyyy_0[i] * pb_x + g_0_xxxxxx_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyz_0[i] = 6.0 * g_0_xxxxx_0_xyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxxxxx_0_xyyz_0[i] * pb_x + g_0_xxxxxx_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyzz_0[i] = 6.0 * g_0_xxxxx_0_xyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxxx_0_xyzz_0[i] * pb_x + g_0_xxxxxx_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xzzz_0[i] = 6.0 * g_0_xxxxx_0_xzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxxx_0_xzzz_0[i] * pb_x + g_0_xxxxxx_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyy_0[i] = 6.0 * g_0_xxxxx_0_yyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyy_0[i] * pb_x +
                                  g_0_xxxxxx_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyz_0[i] = 6.0 * g_0_xxxxx_0_yyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyz_0[i] * pb_x +
                                  g_0_xxxxxx_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyzz_0[i] = 6.0 * g_0_xxxxx_0_yyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyzz_0[i] * pb_x +
                                  g_0_xxxxxx_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yzzz_0[i] = 6.0 * g_0_xxxxx_0_yzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzzz_0[i] * pb_x +
                                  g_0_xxxxxx_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_zzzz_0[i] = 6.0 * g_0_xxxxx_0_zzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzzz_0[i] * pb_x +
                                  g_0_xxxxxx_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : SKSG

    auto g_0_xxxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 15);

    auto g_0_xxxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 16);

    auto g_0_xxxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 17);

    auto g_0_xxxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 18);

    auto g_0_xxxxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 19);

    auto g_0_xxxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 20);

    auto g_0_xxxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 21);

    auto g_0_xxxxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 22);

    auto g_0_xxxxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 23);

    auto g_0_xxxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 24);

    auto g_0_xxxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 25);

    auto g_0_xxxxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 26);

    auto g_0_xxxxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 27);

    auto g_0_xxxxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 28);

    auto g_0_xxxxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 29);

#pragma omp simd aligned(g_0_xxxxxx_0_xxx_1,       \
                             g_0_xxxxxx_0_xxxx_0,  \
                             g_0_xxxxxx_0_xxxx_1,  \
                             g_0_xxxxxx_0_xxxy_0,  \
                             g_0_xxxxxx_0_xxxy_1,  \
                             g_0_xxxxxx_0_xxxz_0,  \
                             g_0_xxxxxx_0_xxxz_1,  \
                             g_0_xxxxxx_0_xxy_1,   \
                             g_0_xxxxxx_0_xxyy_0,  \
                             g_0_xxxxxx_0_xxyy_1,  \
                             g_0_xxxxxx_0_xxyz_0,  \
                             g_0_xxxxxx_0_xxyz_1,  \
                             g_0_xxxxxx_0_xxz_1,   \
                             g_0_xxxxxx_0_xxzz_0,  \
                             g_0_xxxxxx_0_xxzz_1,  \
                             g_0_xxxxxx_0_xyy_1,   \
                             g_0_xxxxxx_0_xyyy_0,  \
                             g_0_xxxxxx_0_xyyy_1,  \
                             g_0_xxxxxx_0_xyyz_0,  \
                             g_0_xxxxxx_0_xyyz_1,  \
                             g_0_xxxxxx_0_xyz_1,   \
                             g_0_xxxxxx_0_xyzz_0,  \
                             g_0_xxxxxx_0_xyzz_1,  \
                             g_0_xxxxxx_0_xzz_1,   \
                             g_0_xxxxxx_0_xzzz_0,  \
                             g_0_xxxxxx_0_xzzz_1,  \
                             g_0_xxxxxx_0_yyy_1,   \
                             g_0_xxxxxx_0_yyyy_0,  \
                             g_0_xxxxxx_0_yyyy_1,  \
                             g_0_xxxxxx_0_yyyz_0,  \
                             g_0_xxxxxx_0_yyyz_1,  \
                             g_0_xxxxxx_0_yyz_1,   \
                             g_0_xxxxxx_0_yyzz_0,  \
                             g_0_xxxxxx_0_yyzz_1,  \
                             g_0_xxxxxx_0_yzz_1,   \
                             g_0_xxxxxx_0_yzzz_0,  \
                             g_0_xxxxxx_0_yzzz_1,  \
                             g_0_xxxxxx_0_zzz_1,   \
                             g_0_xxxxxx_0_zzzz_0,  \
                             g_0_xxxxxx_0_zzzz_1,  \
                             g_0_xxxxxxy_0_xxxx_0, \
                             g_0_xxxxxxy_0_xxxy_0, \
                             g_0_xxxxxxy_0_xxxz_0, \
                             g_0_xxxxxxy_0_xxyy_0, \
                             g_0_xxxxxxy_0_xxyz_0, \
                             g_0_xxxxxxy_0_xxzz_0, \
                             g_0_xxxxxxy_0_xyyy_0, \
                             g_0_xxxxxxy_0_xyyz_0, \
                             g_0_xxxxxxy_0_xyzz_0, \
                             g_0_xxxxxxy_0_xzzz_0, \
                             g_0_xxxxxxy_0_yyyy_0, \
                             g_0_xxxxxxy_0_yyyz_0, \
                             g_0_xxxxxxy_0_yyzz_0, \
                             g_0_xxxxxxy_0_yzzz_0, \
                             g_0_xxxxxxy_0_zzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_xxxx_0[i] = g_0_xxxxxx_0_xxxx_0[i] * pb_y + g_0_xxxxxx_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxy_0[i] = g_0_xxxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxy_0[i] * pb_y + g_0_xxxxxx_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxz_0[i] = g_0_xxxxxx_0_xxxz_0[i] * pb_y + g_0_xxxxxx_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyy_0[i] = 2.0 * g_0_xxxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyy_0[i] * pb_y + g_0_xxxxxx_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyz_0[i] = g_0_xxxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyz_0[i] * pb_y + g_0_xxxxxx_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxzz_0[i] = g_0_xxxxxx_0_xxzz_0[i] * pb_y + g_0_xxxxxx_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyy_0[i] = 3.0 * g_0_xxxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyy_0[i] * pb_y + g_0_xxxxxx_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyz_0[i] = 2.0 * g_0_xxxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyz_0[i] * pb_y + g_0_xxxxxx_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyzz_0[i] = g_0_xxxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzz_0[i] * pb_y + g_0_xxxxxx_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xzzz_0[i] = g_0_xxxxxx_0_xzzz_0[i] * pb_y + g_0_xxxxxx_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyy_0[i] = 4.0 * g_0_xxxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyy_0[i] * pb_y + g_0_xxxxxx_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyz_0[i] = 3.0 * g_0_xxxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyz_0[i] * pb_y + g_0_xxxxxx_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyzz_0[i] = 2.0 * g_0_xxxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzz_0[i] * pb_y + g_0_xxxxxx_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yzzz_0[i] = g_0_xxxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzz_0[i] * pb_y + g_0_xxxxxx_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_zzzz_0[i] = g_0_xxxxxx_0_zzzz_0[i] * pb_y + g_0_xxxxxx_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : SKSG

    auto g_0_xxxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 30);

    auto g_0_xxxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 31);

    auto g_0_xxxxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 32);

    auto g_0_xxxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 33);

    auto g_0_xxxxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 34);

    auto g_0_xxxxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 35);

    auto g_0_xxxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 36);

    auto g_0_xxxxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 37);

    auto g_0_xxxxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 38);

    auto g_0_xxxxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 39);

    auto g_0_xxxxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 40);

    auto g_0_xxxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 41);

    auto g_0_xxxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 42);

    auto g_0_xxxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 43);

    auto g_0_xxxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 44);

#pragma omp simd aligned(g_0_xxxxxx_0_xxx_1,       \
                             g_0_xxxxxx_0_xxxx_0,  \
                             g_0_xxxxxx_0_xxxx_1,  \
                             g_0_xxxxxx_0_xxxy_0,  \
                             g_0_xxxxxx_0_xxxy_1,  \
                             g_0_xxxxxx_0_xxxz_0,  \
                             g_0_xxxxxx_0_xxxz_1,  \
                             g_0_xxxxxx_0_xxy_1,   \
                             g_0_xxxxxx_0_xxyy_0,  \
                             g_0_xxxxxx_0_xxyy_1,  \
                             g_0_xxxxxx_0_xxyz_0,  \
                             g_0_xxxxxx_0_xxyz_1,  \
                             g_0_xxxxxx_0_xxz_1,   \
                             g_0_xxxxxx_0_xxzz_0,  \
                             g_0_xxxxxx_0_xxzz_1,  \
                             g_0_xxxxxx_0_xyy_1,   \
                             g_0_xxxxxx_0_xyyy_0,  \
                             g_0_xxxxxx_0_xyyy_1,  \
                             g_0_xxxxxx_0_xyyz_0,  \
                             g_0_xxxxxx_0_xyyz_1,  \
                             g_0_xxxxxx_0_xyz_1,   \
                             g_0_xxxxxx_0_xyzz_0,  \
                             g_0_xxxxxx_0_xyzz_1,  \
                             g_0_xxxxxx_0_xzz_1,   \
                             g_0_xxxxxx_0_xzzz_0,  \
                             g_0_xxxxxx_0_xzzz_1,  \
                             g_0_xxxxxx_0_yyy_1,   \
                             g_0_xxxxxx_0_yyyy_0,  \
                             g_0_xxxxxx_0_yyyy_1,  \
                             g_0_xxxxxx_0_yyyz_0,  \
                             g_0_xxxxxx_0_yyyz_1,  \
                             g_0_xxxxxx_0_yyz_1,   \
                             g_0_xxxxxx_0_yyzz_0,  \
                             g_0_xxxxxx_0_yyzz_1,  \
                             g_0_xxxxxx_0_yzz_1,   \
                             g_0_xxxxxx_0_yzzz_0,  \
                             g_0_xxxxxx_0_yzzz_1,  \
                             g_0_xxxxxx_0_zzz_1,   \
                             g_0_xxxxxx_0_zzzz_0,  \
                             g_0_xxxxxx_0_zzzz_1,  \
                             g_0_xxxxxxz_0_xxxx_0, \
                             g_0_xxxxxxz_0_xxxy_0, \
                             g_0_xxxxxxz_0_xxxz_0, \
                             g_0_xxxxxxz_0_xxyy_0, \
                             g_0_xxxxxxz_0_xxyz_0, \
                             g_0_xxxxxxz_0_xxzz_0, \
                             g_0_xxxxxxz_0_xyyy_0, \
                             g_0_xxxxxxz_0_xyyz_0, \
                             g_0_xxxxxxz_0_xyzz_0, \
                             g_0_xxxxxxz_0_xzzz_0, \
                             g_0_xxxxxxz_0_yyyy_0, \
                             g_0_xxxxxxz_0_yyyz_0, \
                             g_0_xxxxxxz_0_yyzz_0, \
                             g_0_xxxxxxz_0_yzzz_0, \
                             g_0_xxxxxxz_0_zzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_xxxx_0[i] = g_0_xxxxxx_0_xxxx_0[i] * pb_z + g_0_xxxxxx_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxy_0[i] = g_0_xxxxxx_0_xxxy_0[i] * pb_z + g_0_xxxxxx_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxz_0[i] = g_0_xxxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxz_0[i] * pb_z + g_0_xxxxxx_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyy_0[i] = g_0_xxxxxx_0_xxyy_0[i] * pb_z + g_0_xxxxxx_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyz_0[i] = g_0_xxxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyz_0[i] * pb_z + g_0_xxxxxx_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxzz_0[i] = 2.0 * g_0_xxxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzz_0[i] * pb_z + g_0_xxxxxx_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyy_0[i] = g_0_xxxxxx_0_xyyy_0[i] * pb_z + g_0_xxxxxx_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyz_0[i] = g_0_xxxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyz_0[i] * pb_z + g_0_xxxxxx_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyzz_0[i] = 2.0 * g_0_xxxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzz_0[i] * pb_z + g_0_xxxxxx_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xzzz_0[i] = 3.0 * g_0_xxxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzz_0[i] * pb_z + g_0_xxxxxx_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyy_0[i] = g_0_xxxxxx_0_yyyy_0[i] * pb_z + g_0_xxxxxx_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyz_0[i] = g_0_xxxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyz_0[i] * pb_z + g_0_xxxxxx_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyzz_0[i] = 2.0 * g_0_xxxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzz_0[i] * pb_z + g_0_xxxxxx_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yzzz_0[i] = 3.0 * g_0_xxxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzz_0[i] * pb_z + g_0_xxxxxx_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_zzzz_0[i] = 4.0 * g_0_xxxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_zzzz_0[i] * pb_z + g_0_xxxxxx_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 45-60 components of targeted buffer : SKSG

    auto g_0_xxxxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 45);

    auto g_0_xxxxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 46);

    auto g_0_xxxxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 47);

    auto g_0_xxxxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 48);

    auto g_0_xxxxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 49);

    auto g_0_xxxxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 50);

    auto g_0_xxxxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 51);

    auto g_0_xxxxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 52);

    auto g_0_xxxxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 53);

    auto g_0_xxxxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 54);

    auto g_0_xxxxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 55);

    auto g_0_xxxxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 56);

    auto g_0_xxxxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 57);

    auto g_0_xxxxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 58);

    auto g_0_xxxxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 59);

#pragma omp simd aligned(g_0_xxxxx_0_xxxx_0,       \
                             g_0_xxxxx_0_xxxx_1,   \
                             g_0_xxxxx_0_xxxz_0,   \
                             g_0_xxxxx_0_xxxz_1,   \
                             g_0_xxxxx_0_xxzz_0,   \
                             g_0_xxxxx_0_xxzz_1,   \
                             g_0_xxxxx_0_xzzz_0,   \
                             g_0_xxxxx_0_xzzz_1,   \
                             g_0_xxxxxy_0_xxxx_0,  \
                             g_0_xxxxxy_0_xxxx_1,  \
                             g_0_xxxxxy_0_xxxz_0,  \
                             g_0_xxxxxy_0_xxxz_1,  \
                             g_0_xxxxxy_0_xxzz_0,  \
                             g_0_xxxxxy_0_xxzz_1,  \
                             g_0_xxxxxy_0_xzzz_0,  \
                             g_0_xxxxxy_0_xzzz_1,  \
                             g_0_xxxxxyy_0_xxxx_0, \
                             g_0_xxxxxyy_0_xxxy_0, \
                             g_0_xxxxxyy_0_xxxz_0, \
                             g_0_xxxxxyy_0_xxyy_0, \
                             g_0_xxxxxyy_0_xxyz_0, \
                             g_0_xxxxxyy_0_xxzz_0, \
                             g_0_xxxxxyy_0_xyyy_0, \
                             g_0_xxxxxyy_0_xyyz_0, \
                             g_0_xxxxxyy_0_xyzz_0, \
                             g_0_xxxxxyy_0_xzzz_0, \
                             g_0_xxxxxyy_0_yyyy_0, \
                             g_0_xxxxxyy_0_yyyz_0, \
                             g_0_xxxxxyy_0_yyzz_0, \
                             g_0_xxxxxyy_0_yzzz_0, \
                             g_0_xxxxxyy_0_zzzz_0, \
                             g_0_xxxxyy_0_xxxy_0,  \
                             g_0_xxxxyy_0_xxxy_1,  \
                             g_0_xxxxyy_0_xxy_1,   \
                             g_0_xxxxyy_0_xxyy_0,  \
                             g_0_xxxxyy_0_xxyy_1,  \
                             g_0_xxxxyy_0_xxyz_0,  \
                             g_0_xxxxyy_0_xxyz_1,  \
                             g_0_xxxxyy_0_xyy_1,   \
                             g_0_xxxxyy_0_xyyy_0,  \
                             g_0_xxxxyy_0_xyyy_1,  \
                             g_0_xxxxyy_0_xyyz_0,  \
                             g_0_xxxxyy_0_xyyz_1,  \
                             g_0_xxxxyy_0_xyz_1,   \
                             g_0_xxxxyy_0_xyzz_0,  \
                             g_0_xxxxyy_0_xyzz_1,  \
                             g_0_xxxxyy_0_yyy_1,   \
                             g_0_xxxxyy_0_yyyy_0,  \
                             g_0_xxxxyy_0_yyyy_1,  \
                             g_0_xxxxyy_0_yyyz_0,  \
                             g_0_xxxxyy_0_yyyz_1,  \
                             g_0_xxxxyy_0_yyz_1,   \
                             g_0_xxxxyy_0_yyzz_0,  \
                             g_0_xxxxyy_0_yyzz_1,  \
                             g_0_xxxxyy_0_yzz_1,   \
                             g_0_xxxxyy_0_yzzz_0,  \
                             g_0_xxxxyy_0_yzzz_1,  \
                             g_0_xxxxyy_0_zzzz_0,  \
                             g_0_xxxxyy_0_zzzz_1,  \
                             g_0_xxxyy_0_xxxy_0,   \
                             g_0_xxxyy_0_xxxy_1,   \
                             g_0_xxxyy_0_xxyy_0,   \
                             g_0_xxxyy_0_xxyy_1,   \
                             g_0_xxxyy_0_xxyz_0,   \
                             g_0_xxxyy_0_xxyz_1,   \
                             g_0_xxxyy_0_xyyy_0,   \
                             g_0_xxxyy_0_xyyy_1,   \
                             g_0_xxxyy_0_xyyz_0,   \
                             g_0_xxxyy_0_xyyz_1,   \
                             g_0_xxxyy_0_xyzz_0,   \
                             g_0_xxxyy_0_xyzz_1,   \
                             g_0_xxxyy_0_yyyy_0,   \
                             g_0_xxxyy_0_yyyy_1,   \
                             g_0_xxxyy_0_yyyz_0,   \
                             g_0_xxxyy_0_yyyz_1,   \
                             g_0_xxxyy_0_yyzz_0,   \
                             g_0_xxxyy_0_yyzz_1,   \
                             g_0_xxxyy_0_yzzz_0,   \
                             g_0_xxxyy_0_yzzz_1,   \
                             g_0_xxxyy_0_zzzz_0,   \
                             g_0_xxxyy_0_zzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_xxxx_0[i] =
            g_0_xxxxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxx_0[i] * pb_y + g_0_xxxxxy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxy_0[i] = 4.0 * g_0_xxxyy_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxy_0[i] * pb_x + g_0_xxxxyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxz_0[i] =
            g_0_xxxxx_0_xxxz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxz_0[i] * pb_y + g_0_xxxxxy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxyy_0[i] = 4.0 * g_0_xxxyy_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyy_0[i] * pb_x + g_0_xxxxyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyz_0[i] = 4.0 * g_0_xxxyy_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyz_0[i] * pb_x + g_0_xxxxyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxzz_0[i] =
            g_0_xxxxx_0_xxzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxzz_0[i] * pb_y + g_0_xxxxxy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xyyy_0[i] = 4.0 * g_0_xxxyy_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyy_1[i] * fi_abcd_0 +
                                  g_0_xxxxyy_0_xyyy_0[i] * pb_x + g_0_xxxxyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyz_0[i] = 4.0 * g_0_xxxyy_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxxxyy_0_xyyz_0[i] * pb_x + g_0_xxxxyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyzz_0[i] = 4.0 * g_0_xxxyy_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxyy_0_xyzz_0[i] * pb_x + g_0_xxxxyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xzzz_0[i] =
            g_0_xxxxx_0_xzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xzzz_0[i] * pb_y + g_0_xxxxxy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_yyyy_0[i] = 4.0 * g_0_xxxyy_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyy_0[i] * pb_x +
                                  g_0_xxxxyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyz_0[i] = 4.0 * g_0_xxxyy_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyz_0[i] * pb_x +
                                  g_0_xxxxyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyzz_0[i] = 4.0 * g_0_xxxyy_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyzz_0[i] * pb_x +
                                  g_0_xxxxyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yzzz_0[i] = 4.0 * g_0_xxxyy_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzzz_0[i] * pb_x +
                                  g_0_xxxxyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_zzzz_0[i] = 4.0 * g_0_xxxyy_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_zzzz_0[i] * pb_x +
                                  g_0_xxxxyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 60-75 components of targeted buffer : SKSG

    auto g_0_xxxxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 60);

    auto g_0_xxxxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 61);

    auto g_0_xxxxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 62);

    auto g_0_xxxxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 63);

    auto g_0_xxxxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 64);

    auto g_0_xxxxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 65);

    auto g_0_xxxxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 66);

    auto g_0_xxxxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 67);

    auto g_0_xxxxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 68);

    auto g_0_xxxxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 69);

    auto g_0_xxxxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 70);

    auto g_0_xxxxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 71);

    auto g_0_xxxxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 72);

    auto g_0_xxxxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 73);

    auto g_0_xxxxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 74);

#pragma omp simd aligned(g_0_xxxxxy_0_xxxy_0,      \
                             g_0_xxxxxy_0_xxxy_1,  \
                             g_0_xxxxxy_0_xxyy_0,  \
                             g_0_xxxxxy_0_xxyy_1,  \
                             g_0_xxxxxy_0_xyyy_0,  \
                             g_0_xxxxxy_0_xyyy_1,  \
                             g_0_xxxxxy_0_yyyy_0,  \
                             g_0_xxxxxy_0_yyyy_1,  \
                             g_0_xxxxxyz_0_xxxx_0, \
                             g_0_xxxxxyz_0_xxxy_0, \
                             g_0_xxxxxyz_0_xxxz_0, \
                             g_0_xxxxxyz_0_xxyy_0, \
                             g_0_xxxxxyz_0_xxyz_0, \
                             g_0_xxxxxyz_0_xxzz_0, \
                             g_0_xxxxxyz_0_xyyy_0, \
                             g_0_xxxxxyz_0_xyyz_0, \
                             g_0_xxxxxyz_0_xyzz_0, \
                             g_0_xxxxxyz_0_xzzz_0, \
                             g_0_xxxxxyz_0_yyyy_0, \
                             g_0_xxxxxyz_0_yyyz_0, \
                             g_0_xxxxxyz_0_yyzz_0, \
                             g_0_xxxxxyz_0_yzzz_0, \
                             g_0_xxxxxyz_0_zzzz_0, \
                             g_0_xxxxxz_0_xxxx_0,  \
                             g_0_xxxxxz_0_xxxx_1,  \
                             g_0_xxxxxz_0_xxxz_0,  \
                             g_0_xxxxxz_0_xxxz_1,  \
                             g_0_xxxxxz_0_xxyz_0,  \
                             g_0_xxxxxz_0_xxyz_1,  \
                             g_0_xxxxxz_0_xxz_1,   \
                             g_0_xxxxxz_0_xxzz_0,  \
                             g_0_xxxxxz_0_xxzz_1,  \
                             g_0_xxxxxz_0_xyyz_0,  \
                             g_0_xxxxxz_0_xyyz_1,  \
                             g_0_xxxxxz_0_xyz_1,   \
                             g_0_xxxxxz_0_xyzz_0,  \
                             g_0_xxxxxz_0_xyzz_1,  \
                             g_0_xxxxxz_0_xzz_1,   \
                             g_0_xxxxxz_0_xzzz_0,  \
                             g_0_xxxxxz_0_xzzz_1,  \
                             g_0_xxxxxz_0_yyyz_0,  \
                             g_0_xxxxxz_0_yyyz_1,  \
                             g_0_xxxxxz_0_yyz_1,   \
                             g_0_xxxxxz_0_yyzz_0,  \
                             g_0_xxxxxz_0_yyzz_1,  \
                             g_0_xxxxxz_0_yzz_1,   \
                             g_0_xxxxxz_0_yzzz_0,  \
                             g_0_xxxxxz_0_yzzz_1,  \
                             g_0_xxxxxz_0_zzz_1,   \
                             g_0_xxxxxz_0_zzzz_0,  \
                             g_0_xxxxxz_0_zzzz_1,  \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyz_0_xxxx_0[i] = g_0_xxxxxz_0_xxxx_0[i] * pb_y + g_0_xxxxxz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxy_0[i] = g_0_xxxxxy_0_xxxy_0[i] * pb_z + g_0_xxxxxy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxz_0[i] = g_0_xxxxxz_0_xxxz_0[i] * pb_y + g_0_xxxxxz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyy_0[i] = g_0_xxxxxy_0_xxyy_0[i] * pb_z + g_0_xxxxxy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxyz_0[i] = g_0_xxxxxz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyz_0[i] * pb_y + g_0_xxxxxz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxzz_0[i] = g_0_xxxxxz_0_xxzz_0[i] * pb_y + g_0_xxxxxz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyy_0[i] = g_0_xxxxxy_0_xyyy_0[i] * pb_z + g_0_xxxxxy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xyyz_0[i] = 2.0 * g_0_xxxxxz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyz_0[i] * pb_y + g_0_xxxxxz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyzz_0[i] = g_0_xxxxxz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyzz_0[i] * pb_y + g_0_xxxxxz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xzzz_0[i] = g_0_xxxxxz_0_xzzz_0[i] * pb_y + g_0_xxxxxz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyy_0[i] = g_0_xxxxxy_0_yyyy_0[i] * pb_z + g_0_xxxxxy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_yyyz_0[i] = 3.0 * g_0_xxxxxz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyz_0[i] * pb_y + g_0_xxxxxz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyzz_0[i] = 2.0 * g_0_xxxxxz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyzz_0[i] * pb_y + g_0_xxxxxz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yzzz_0[i] = g_0_xxxxxz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yzzz_0[i] * pb_y + g_0_xxxxxz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_zzzz_0[i] = g_0_xxxxxz_0_zzzz_0[i] * pb_y + g_0_xxxxxz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 75-90 components of targeted buffer : SKSG

    auto g_0_xxxxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 75);

    auto g_0_xxxxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 76);

    auto g_0_xxxxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 77);

    auto g_0_xxxxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 78);

    auto g_0_xxxxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 79);

    auto g_0_xxxxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 80);

    auto g_0_xxxxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 81);

    auto g_0_xxxxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 82);

    auto g_0_xxxxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 83);

    auto g_0_xxxxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 84);

    auto g_0_xxxxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 85);

    auto g_0_xxxxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 86);

    auto g_0_xxxxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 87);

    auto g_0_xxxxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 88);

    auto g_0_xxxxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 89);

#pragma omp simd aligned(g_0_xxxxx_0_xxxx_0,       \
                             g_0_xxxxx_0_xxxx_1,   \
                             g_0_xxxxx_0_xxxy_0,   \
                             g_0_xxxxx_0_xxxy_1,   \
                             g_0_xxxxx_0_xxyy_0,   \
                             g_0_xxxxx_0_xxyy_1,   \
                             g_0_xxxxx_0_xyyy_0,   \
                             g_0_xxxxx_0_xyyy_1,   \
                             g_0_xxxxxz_0_xxxx_0,  \
                             g_0_xxxxxz_0_xxxx_1,  \
                             g_0_xxxxxz_0_xxxy_0,  \
                             g_0_xxxxxz_0_xxxy_1,  \
                             g_0_xxxxxz_0_xxyy_0,  \
                             g_0_xxxxxz_0_xxyy_1,  \
                             g_0_xxxxxz_0_xyyy_0,  \
                             g_0_xxxxxz_0_xyyy_1,  \
                             g_0_xxxxxzz_0_xxxx_0, \
                             g_0_xxxxxzz_0_xxxy_0, \
                             g_0_xxxxxzz_0_xxxz_0, \
                             g_0_xxxxxzz_0_xxyy_0, \
                             g_0_xxxxxzz_0_xxyz_0, \
                             g_0_xxxxxzz_0_xxzz_0, \
                             g_0_xxxxxzz_0_xyyy_0, \
                             g_0_xxxxxzz_0_xyyz_0, \
                             g_0_xxxxxzz_0_xyzz_0, \
                             g_0_xxxxxzz_0_xzzz_0, \
                             g_0_xxxxxzz_0_yyyy_0, \
                             g_0_xxxxxzz_0_yyyz_0, \
                             g_0_xxxxxzz_0_yyzz_0, \
                             g_0_xxxxxzz_0_yzzz_0, \
                             g_0_xxxxxzz_0_zzzz_0, \
                             g_0_xxxxzz_0_xxxz_0,  \
                             g_0_xxxxzz_0_xxxz_1,  \
                             g_0_xxxxzz_0_xxyz_0,  \
                             g_0_xxxxzz_0_xxyz_1,  \
                             g_0_xxxxzz_0_xxz_1,   \
                             g_0_xxxxzz_0_xxzz_0,  \
                             g_0_xxxxzz_0_xxzz_1,  \
                             g_0_xxxxzz_0_xyyz_0,  \
                             g_0_xxxxzz_0_xyyz_1,  \
                             g_0_xxxxzz_0_xyz_1,   \
                             g_0_xxxxzz_0_xyzz_0,  \
                             g_0_xxxxzz_0_xyzz_1,  \
                             g_0_xxxxzz_0_xzz_1,   \
                             g_0_xxxxzz_0_xzzz_0,  \
                             g_0_xxxxzz_0_xzzz_1,  \
                             g_0_xxxxzz_0_yyyy_0,  \
                             g_0_xxxxzz_0_yyyy_1,  \
                             g_0_xxxxzz_0_yyyz_0,  \
                             g_0_xxxxzz_0_yyyz_1,  \
                             g_0_xxxxzz_0_yyz_1,   \
                             g_0_xxxxzz_0_yyzz_0,  \
                             g_0_xxxxzz_0_yyzz_1,  \
                             g_0_xxxxzz_0_yzz_1,   \
                             g_0_xxxxzz_0_yzzz_0,  \
                             g_0_xxxxzz_0_yzzz_1,  \
                             g_0_xxxxzz_0_zzz_1,   \
                             g_0_xxxxzz_0_zzzz_0,  \
                             g_0_xxxxzz_0_zzzz_1,  \
                             g_0_xxxzz_0_xxxz_0,   \
                             g_0_xxxzz_0_xxxz_1,   \
                             g_0_xxxzz_0_xxyz_0,   \
                             g_0_xxxzz_0_xxyz_1,   \
                             g_0_xxxzz_0_xxzz_0,   \
                             g_0_xxxzz_0_xxzz_1,   \
                             g_0_xxxzz_0_xyyz_0,   \
                             g_0_xxxzz_0_xyyz_1,   \
                             g_0_xxxzz_0_xyzz_0,   \
                             g_0_xxxzz_0_xyzz_1,   \
                             g_0_xxxzz_0_xzzz_0,   \
                             g_0_xxxzz_0_xzzz_1,   \
                             g_0_xxxzz_0_yyyy_0,   \
                             g_0_xxxzz_0_yyyy_1,   \
                             g_0_xxxzz_0_yyyz_0,   \
                             g_0_xxxzz_0_yyyz_1,   \
                             g_0_xxxzz_0_yyzz_0,   \
                             g_0_xxxzz_0_yyzz_1,   \
                             g_0_xxxzz_0_yzzz_0,   \
                             g_0_xxxzz_0_yzzz_1,   \
                             g_0_xxxzz_0_zzzz_0,   \
                             g_0_xxxzz_0_zzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_xxxx_0[i] =
            g_0_xxxxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxx_0[i] * pb_z + g_0_xxxxxz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxy_0[i] =
            g_0_xxxxx_0_xxxy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxy_0[i] * pb_z + g_0_xxxxxz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxz_0[i] = 4.0 * g_0_xxxzz_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxz_0[i] * pb_x + g_0_xxxxzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyy_0[i] =
            g_0_xxxxx_0_xxyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxyy_0[i] * pb_z + g_0_xxxxxz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxyz_0[i] = 4.0 * g_0_xxxzz_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyz_0[i] * pb_x + g_0_xxxxzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxzz_0[i] = 4.0 * g_0_xxxzz_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxzz_0[i] * pb_x + g_0_xxxxzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyy_0[i] =
            g_0_xxxxx_0_xyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xyyy_0[i] * pb_z + g_0_xxxxxz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xyyz_0[i] = 4.0 * g_0_xxxzz_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxxxzz_0_xyyz_0[i] * pb_x + g_0_xxxxzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyzz_0[i] = 4.0 * g_0_xxxzz_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxzz_0_xyzz_0[i] * pb_x + g_0_xxxxzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xzzz_0[i] = 4.0 * g_0_xxxzz_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_xxxxzz_0_xzzz_0[i] * pb_x + g_0_xxxxzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyy_0[i] = 4.0 * g_0_xxxzz_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyy_0[i] * pb_x +
                                  g_0_xxxxzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyz_0[i] = 4.0 * g_0_xxxzz_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyz_0[i] * pb_x +
                                  g_0_xxxxzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyzz_0[i] = 4.0 * g_0_xxxzz_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyzz_0[i] * pb_x +
                                  g_0_xxxxzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yzzz_0[i] = 4.0 * g_0_xxxzz_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzzz_0[i] * pb_x +
                                  g_0_xxxxzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_zzzz_0[i] = 4.0 * g_0_xxxzz_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzzz_0[i] * pb_x +
                                  g_0_xxxxzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 90-105 components of targeted buffer : SKSG

    auto g_0_xxxxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 90);

    auto g_0_xxxxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 91);

    auto g_0_xxxxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 92);

    auto g_0_xxxxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 93);

    auto g_0_xxxxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 94);

    auto g_0_xxxxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 95);

    auto g_0_xxxxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 96);

    auto g_0_xxxxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 97);

    auto g_0_xxxxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 98);

    auto g_0_xxxxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 99);

    auto g_0_xxxxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 100);

    auto g_0_xxxxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 101);

    auto g_0_xxxxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 102);

    auto g_0_xxxxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 103);

    auto g_0_xxxxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 104);

#pragma omp simd aligned(g_0_xxxxy_0_xxxx_0,       \
                             g_0_xxxxy_0_xxxx_1,   \
                             g_0_xxxxy_0_xxxz_0,   \
                             g_0_xxxxy_0_xxxz_1,   \
                             g_0_xxxxy_0_xxzz_0,   \
                             g_0_xxxxy_0_xxzz_1,   \
                             g_0_xxxxy_0_xzzz_0,   \
                             g_0_xxxxy_0_xzzz_1,   \
                             g_0_xxxxyy_0_xxxx_0,  \
                             g_0_xxxxyy_0_xxxx_1,  \
                             g_0_xxxxyy_0_xxxz_0,  \
                             g_0_xxxxyy_0_xxxz_1,  \
                             g_0_xxxxyy_0_xxzz_0,  \
                             g_0_xxxxyy_0_xxzz_1,  \
                             g_0_xxxxyy_0_xzzz_0,  \
                             g_0_xxxxyy_0_xzzz_1,  \
                             g_0_xxxxyyy_0_xxxx_0, \
                             g_0_xxxxyyy_0_xxxy_0, \
                             g_0_xxxxyyy_0_xxxz_0, \
                             g_0_xxxxyyy_0_xxyy_0, \
                             g_0_xxxxyyy_0_xxyz_0, \
                             g_0_xxxxyyy_0_xxzz_0, \
                             g_0_xxxxyyy_0_xyyy_0, \
                             g_0_xxxxyyy_0_xyyz_0, \
                             g_0_xxxxyyy_0_xyzz_0, \
                             g_0_xxxxyyy_0_xzzz_0, \
                             g_0_xxxxyyy_0_yyyy_0, \
                             g_0_xxxxyyy_0_yyyz_0, \
                             g_0_xxxxyyy_0_yyzz_0, \
                             g_0_xxxxyyy_0_yzzz_0, \
                             g_0_xxxxyyy_0_zzzz_0, \
                             g_0_xxxyyy_0_xxxy_0,  \
                             g_0_xxxyyy_0_xxxy_1,  \
                             g_0_xxxyyy_0_xxy_1,   \
                             g_0_xxxyyy_0_xxyy_0,  \
                             g_0_xxxyyy_0_xxyy_1,  \
                             g_0_xxxyyy_0_xxyz_0,  \
                             g_0_xxxyyy_0_xxyz_1,  \
                             g_0_xxxyyy_0_xyy_1,   \
                             g_0_xxxyyy_0_xyyy_0,  \
                             g_0_xxxyyy_0_xyyy_1,  \
                             g_0_xxxyyy_0_xyyz_0,  \
                             g_0_xxxyyy_0_xyyz_1,  \
                             g_0_xxxyyy_0_xyz_1,   \
                             g_0_xxxyyy_0_xyzz_0,  \
                             g_0_xxxyyy_0_xyzz_1,  \
                             g_0_xxxyyy_0_yyy_1,   \
                             g_0_xxxyyy_0_yyyy_0,  \
                             g_0_xxxyyy_0_yyyy_1,  \
                             g_0_xxxyyy_0_yyyz_0,  \
                             g_0_xxxyyy_0_yyyz_1,  \
                             g_0_xxxyyy_0_yyz_1,   \
                             g_0_xxxyyy_0_yyzz_0,  \
                             g_0_xxxyyy_0_yyzz_1,  \
                             g_0_xxxyyy_0_yzz_1,   \
                             g_0_xxxyyy_0_yzzz_0,  \
                             g_0_xxxyyy_0_yzzz_1,  \
                             g_0_xxxyyy_0_zzzz_0,  \
                             g_0_xxxyyy_0_zzzz_1,  \
                             g_0_xxyyy_0_xxxy_0,   \
                             g_0_xxyyy_0_xxxy_1,   \
                             g_0_xxyyy_0_xxyy_0,   \
                             g_0_xxyyy_0_xxyy_1,   \
                             g_0_xxyyy_0_xxyz_0,   \
                             g_0_xxyyy_0_xxyz_1,   \
                             g_0_xxyyy_0_xyyy_0,   \
                             g_0_xxyyy_0_xyyy_1,   \
                             g_0_xxyyy_0_xyyz_0,   \
                             g_0_xxyyy_0_xyyz_1,   \
                             g_0_xxyyy_0_xyzz_0,   \
                             g_0_xxyyy_0_xyzz_1,   \
                             g_0_xxyyy_0_yyyy_0,   \
                             g_0_xxyyy_0_yyyy_1,   \
                             g_0_xxyyy_0_yyyz_0,   \
                             g_0_xxyyy_0_yyyz_1,   \
                             g_0_xxyyy_0_yyzz_0,   \
                             g_0_xxyyy_0_yyzz_1,   \
                             g_0_xxyyy_0_yzzz_0,   \
                             g_0_xxyyy_0_yzzz_1,   \
                             g_0_xxyyy_0_zzzz_0,   \
                             g_0_xxyyy_0_zzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_xxxx_0[i] = 2.0 * g_0_xxxxy_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxx_0[i] * pb_y +
                                  g_0_xxxxyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxy_0[i] = 3.0 * g_0_xxyyy_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxy_0[i] * pb_x + g_0_xxxyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxz_0[i] = 2.0 * g_0_xxxxy_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxz_0[i] * pb_y +
                                  g_0_xxxxyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxyy_0[i] = 3.0 * g_0_xxyyy_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyy_0[i] * pb_x + g_0_xxxyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyz_0[i] = 3.0 * g_0_xxyyy_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyz_0[i] * pb_x + g_0_xxxyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxzz_0[i] = 2.0 * g_0_xxxxy_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxzz_0[i] * pb_y +
                                  g_0_xxxxyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xyyy_0[i] = 3.0 * g_0_xxyyy_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyy_1[i] * fi_abcd_0 +
                                  g_0_xxxyyy_0_xyyy_0[i] * pb_x + g_0_xxxyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyz_0[i] = 3.0 * g_0_xxyyy_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxxyyy_0_xyyz_0[i] * pb_x + g_0_xxxyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyzz_0[i] = 3.0 * g_0_xxyyy_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxxyyy_0_xyzz_0[i] * pb_x + g_0_xxxyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xzzz_0[i] = 2.0 * g_0_xxxxy_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xzzz_0[i] * pb_y +
                                  g_0_xxxxyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_yyyy_0[i] = 3.0 * g_0_xxyyy_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyy_0[i] * pb_x +
                                  g_0_xxxyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyz_0[i] = 3.0 * g_0_xxyyy_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyz_0[i] * pb_x +
                                  g_0_xxxyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyzz_0[i] = 3.0 * g_0_xxyyy_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyzz_0[i] * pb_x +
                                  g_0_xxxyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yzzz_0[i] = 3.0 * g_0_xxyyy_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzzz_0[i] * pb_x +
                                  g_0_xxxyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_zzzz_0[i] = 3.0 * g_0_xxyyy_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_zzzz_0[i] * pb_x +
                                  g_0_xxxyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 105-120 components of targeted buffer : SKSG

    auto g_0_xxxxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 105);

    auto g_0_xxxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 106);

    auto g_0_xxxxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 107);

    auto g_0_xxxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 108);

    auto g_0_xxxxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 109);

    auto g_0_xxxxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 110);

    auto g_0_xxxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 111);

    auto g_0_xxxxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 112);

    auto g_0_xxxxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 113);

    auto g_0_xxxxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 114);

    auto g_0_xxxxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 115);

    auto g_0_xxxxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 116);

    auto g_0_xxxxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 117);

    auto g_0_xxxxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 118);

    auto g_0_xxxxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 119);

#pragma omp simd aligned(g_0_xxxxyy_0_xxx_1,       \
                             g_0_xxxxyy_0_xxxx_0,  \
                             g_0_xxxxyy_0_xxxx_1,  \
                             g_0_xxxxyy_0_xxxy_0,  \
                             g_0_xxxxyy_0_xxxy_1,  \
                             g_0_xxxxyy_0_xxxz_0,  \
                             g_0_xxxxyy_0_xxxz_1,  \
                             g_0_xxxxyy_0_xxy_1,   \
                             g_0_xxxxyy_0_xxyy_0,  \
                             g_0_xxxxyy_0_xxyy_1,  \
                             g_0_xxxxyy_0_xxyz_0,  \
                             g_0_xxxxyy_0_xxyz_1,  \
                             g_0_xxxxyy_0_xxz_1,   \
                             g_0_xxxxyy_0_xxzz_0,  \
                             g_0_xxxxyy_0_xxzz_1,  \
                             g_0_xxxxyy_0_xyy_1,   \
                             g_0_xxxxyy_0_xyyy_0,  \
                             g_0_xxxxyy_0_xyyy_1,  \
                             g_0_xxxxyy_0_xyyz_0,  \
                             g_0_xxxxyy_0_xyyz_1,  \
                             g_0_xxxxyy_0_xyz_1,   \
                             g_0_xxxxyy_0_xyzz_0,  \
                             g_0_xxxxyy_0_xyzz_1,  \
                             g_0_xxxxyy_0_xzz_1,   \
                             g_0_xxxxyy_0_xzzz_0,  \
                             g_0_xxxxyy_0_xzzz_1,  \
                             g_0_xxxxyy_0_yyy_1,   \
                             g_0_xxxxyy_0_yyyy_0,  \
                             g_0_xxxxyy_0_yyyy_1,  \
                             g_0_xxxxyy_0_yyyz_0,  \
                             g_0_xxxxyy_0_yyyz_1,  \
                             g_0_xxxxyy_0_yyz_1,   \
                             g_0_xxxxyy_0_yyzz_0,  \
                             g_0_xxxxyy_0_yyzz_1,  \
                             g_0_xxxxyy_0_yzz_1,   \
                             g_0_xxxxyy_0_yzzz_0,  \
                             g_0_xxxxyy_0_yzzz_1,  \
                             g_0_xxxxyy_0_zzz_1,   \
                             g_0_xxxxyy_0_zzzz_0,  \
                             g_0_xxxxyy_0_zzzz_1,  \
                             g_0_xxxxyyz_0_xxxx_0, \
                             g_0_xxxxyyz_0_xxxy_0, \
                             g_0_xxxxyyz_0_xxxz_0, \
                             g_0_xxxxyyz_0_xxyy_0, \
                             g_0_xxxxyyz_0_xxyz_0, \
                             g_0_xxxxyyz_0_xxzz_0, \
                             g_0_xxxxyyz_0_xyyy_0, \
                             g_0_xxxxyyz_0_xyyz_0, \
                             g_0_xxxxyyz_0_xyzz_0, \
                             g_0_xxxxyyz_0_xzzz_0, \
                             g_0_xxxxyyz_0_yyyy_0, \
                             g_0_xxxxyyz_0_yyyz_0, \
                             g_0_xxxxyyz_0_yyzz_0, \
                             g_0_xxxxyyz_0_yzzz_0, \
                             g_0_xxxxyyz_0_zzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_xxxx_0[i] = g_0_xxxxyy_0_xxxx_0[i] * pb_z + g_0_xxxxyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxy_0[i] = g_0_xxxxyy_0_xxxy_0[i] * pb_z + g_0_xxxxyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxz_0[i] = g_0_xxxxyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxz_0[i] * pb_z + g_0_xxxxyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyy_0[i] = g_0_xxxxyy_0_xxyy_0[i] * pb_z + g_0_xxxxyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyz_0[i] = g_0_xxxxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyz_0[i] * pb_z + g_0_xxxxyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxzz_0[i] = 2.0 * g_0_xxxxyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxzz_0[i] * pb_z + g_0_xxxxyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyy_0[i] = g_0_xxxxyy_0_xyyy_0[i] * pb_z + g_0_xxxxyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyz_0[i] = g_0_xxxxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyz_0[i] * pb_z + g_0_xxxxyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyzz_0[i] = 2.0 * g_0_xxxxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzz_0[i] * pb_z + g_0_xxxxyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xzzz_0[i] = 3.0 * g_0_xxxxyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xzzz_0[i] * pb_z + g_0_xxxxyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyy_0[i] = g_0_xxxxyy_0_yyyy_0[i] * pb_z + g_0_xxxxyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyz_0[i] = g_0_xxxxyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyz_0[i] * pb_z + g_0_xxxxyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyzz_0[i] = 2.0 * g_0_xxxxyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyzz_0[i] * pb_z + g_0_xxxxyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yzzz_0[i] = 3.0 * g_0_xxxxyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yzzz_0[i] * pb_z + g_0_xxxxyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_zzzz_0[i] = 4.0 * g_0_xxxxyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_zzzz_0[i] * pb_z + g_0_xxxxyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 120-135 components of targeted buffer : SKSG

    auto g_0_xxxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 120);

    auto g_0_xxxxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 121);

    auto g_0_xxxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 122);

    auto g_0_xxxxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 123);

    auto g_0_xxxxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 124);

    auto g_0_xxxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 125);

    auto g_0_xxxxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 126);

    auto g_0_xxxxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 127);

    auto g_0_xxxxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 128);

    auto g_0_xxxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 129);

    auto g_0_xxxxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 130);

    auto g_0_xxxxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 131);

    auto g_0_xxxxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 132);

    auto g_0_xxxxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 133);

    auto g_0_xxxxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 134);

#pragma omp simd aligned(g_0_xxxxyzz_0_xxxx_0,     \
                             g_0_xxxxyzz_0_xxxy_0, \
                             g_0_xxxxyzz_0_xxxz_0, \
                             g_0_xxxxyzz_0_xxyy_0, \
                             g_0_xxxxyzz_0_xxyz_0, \
                             g_0_xxxxyzz_0_xxzz_0, \
                             g_0_xxxxyzz_0_xyyy_0, \
                             g_0_xxxxyzz_0_xyyz_0, \
                             g_0_xxxxyzz_0_xyzz_0, \
                             g_0_xxxxyzz_0_xzzz_0, \
                             g_0_xxxxyzz_0_yyyy_0, \
                             g_0_xxxxyzz_0_yyyz_0, \
                             g_0_xxxxyzz_0_yyzz_0, \
                             g_0_xxxxyzz_0_yzzz_0, \
                             g_0_xxxxyzz_0_zzzz_0, \
                             g_0_xxxxzz_0_xxx_1,   \
                             g_0_xxxxzz_0_xxxx_0,  \
                             g_0_xxxxzz_0_xxxx_1,  \
                             g_0_xxxxzz_0_xxxy_0,  \
                             g_0_xxxxzz_0_xxxy_1,  \
                             g_0_xxxxzz_0_xxxz_0,  \
                             g_0_xxxxzz_0_xxxz_1,  \
                             g_0_xxxxzz_0_xxy_1,   \
                             g_0_xxxxzz_0_xxyy_0,  \
                             g_0_xxxxzz_0_xxyy_1,  \
                             g_0_xxxxzz_0_xxyz_0,  \
                             g_0_xxxxzz_0_xxyz_1,  \
                             g_0_xxxxzz_0_xxz_1,   \
                             g_0_xxxxzz_0_xxzz_0,  \
                             g_0_xxxxzz_0_xxzz_1,  \
                             g_0_xxxxzz_0_xyy_1,   \
                             g_0_xxxxzz_0_xyyy_0,  \
                             g_0_xxxxzz_0_xyyy_1,  \
                             g_0_xxxxzz_0_xyyz_0,  \
                             g_0_xxxxzz_0_xyyz_1,  \
                             g_0_xxxxzz_0_xyz_1,   \
                             g_0_xxxxzz_0_xyzz_0,  \
                             g_0_xxxxzz_0_xyzz_1,  \
                             g_0_xxxxzz_0_xzz_1,   \
                             g_0_xxxxzz_0_xzzz_0,  \
                             g_0_xxxxzz_0_xzzz_1,  \
                             g_0_xxxxzz_0_yyy_1,   \
                             g_0_xxxxzz_0_yyyy_0,  \
                             g_0_xxxxzz_0_yyyy_1,  \
                             g_0_xxxxzz_0_yyyz_0,  \
                             g_0_xxxxzz_0_yyyz_1,  \
                             g_0_xxxxzz_0_yyz_1,   \
                             g_0_xxxxzz_0_yyzz_0,  \
                             g_0_xxxxzz_0_yyzz_1,  \
                             g_0_xxxxzz_0_yzz_1,   \
                             g_0_xxxxzz_0_yzzz_0,  \
                             g_0_xxxxzz_0_yzzz_1,  \
                             g_0_xxxxzz_0_zzz_1,   \
                             g_0_xxxxzz_0_zzzz_0,  \
                             g_0_xxxxzz_0_zzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_xxxx_0[i] = g_0_xxxxzz_0_xxxx_0[i] * pb_y + g_0_xxxxzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxy_0[i] = g_0_xxxxzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxy_0[i] * pb_y + g_0_xxxxzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxz_0[i] = g_0_xxxxzz_0_xxxz_0[i] * pb_y + g_0_xxxxzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyy_0[i] = 2.0 * g_0_xxxxzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyy_0[i] * pb_y + g_0_xxxxzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyz_0[i] = g_0_xxxxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyz_0[i] * pb_y + g_0_xxxxzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxzz_0[i] = g_0_xxxxzz_0_xxzz_0[i] * pb_y + g_0_xxxxzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyy_0[i] = 3.0 * g_0_xxxxzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyy_0[i] * pb_y + g_0_xxxxzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyz_0[i] = 2.0 * g_0_xxxxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyz_0[i] * pb_y + g_0_xxxxzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyzz_0[i] = g_0_xxxxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzz_0[i] * pb_y + g_0_xxxxzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xzzz_0[i] = g_0_xxxxzz_0_xzzz_0[i] * pb_y + g_0_xxxxzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyy_0[i] = 4.0 * g_0_xxxxzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyy_0[i] * pb_y + g_0_xxxxzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyz_0[i] = 3.0 * g_0_xxxxzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyz_0[i] * pb_y + g_0_xxxxzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyzz_0[i] = 2.0 * g_0_xxxxzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyzz_0[i] * pb_y + g_0_xxxxzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yzzz_0[i] = g_0_xxxxzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yzzz_0[i] * pb_y + g_0_xxxxzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_zzzz_0[i] = g_0_xxxxzz_0_zzzz_0[i] * pb_y + g_0_xxxxzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 135-150 components of targeted buffer : SKSG

    auto g_0_xxxxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 135);

    auto g_0_xxxxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 136);

    auto g_0_xxxxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 137);

    auto g_0_xxxxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 138);

    auto g_0_xxxxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 139);

    auto g_0_xxxxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 140);

    auto g_0_xxxxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 141);

    auto g_0_xxxxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 142);

    auto g_0_xxxxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 143);

    auto g_0_xxxxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 144);

    auto g_0_xxxxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 145);

    auto g_0_xxxxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 146);

    auto g_0_xxxxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 147);

    auto g_0_xxxxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 148);

    auto g_0_xxxxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 149);

#pragma omp simd aligned(g_0_xxxxz_0_xxxx_0,       \
                             g_0_xxxxz_0_xxxx_1,   \
                             g_0_xxxxz_0_xxxy_0,   \
                             g_0_xxxxz_0_xxxy_1,   \
                             g_0_xxxxz_0_xxyy_0,   \
                             g_0_xxxxz_0_xxyy_1,   \
                             g_0_xxxxz_0_xyyy_0,   \
                             g_0_xxxxz_0_xyyy_1,   \
                             g_0_xxxxzz_0_xxxx_0,  \
                             g_0_xxxxzz_0_xxxx_1,  \
                             g_0_xxxxzz_0_xxxy_0,  \
                             g_0_xxxxzz_0_xxxy_1,  \
                             g_0_xxxxzz_0_xxyy_0,  \
                             g_0_xxxxzz_0_xxyy_1,  \
                             g_0_xxxxzz_0_xyyy_0,  \
                             g_0_xxxxzz_0_xyyy_1,  \
                             g_0_xxxxzzz_0_xxxx_0, \
                             g_0_xxxxzzz_0_xxxy_0, \
                             g_0_xxxxzzz_0_xxxz_0, \
                             g_0_xxxxzzz_0_xxyy_0, \
                             g_0_xxxxzzz_0_xxyz_0, \
                             g_0_xxxxzzz_0_xxzz_0, \
                             g_0_xxxxzzz_0_xyyy_0, \
                             g_0_xxxxzzz_0_xyyz_0, \
                             g_0_xxxxzzz_0_xyzz_0, \
                             g_0_xxxxzzz_0_xzzz_0, \
                             g_0_xxxxzzz_0_yyyy_0, \
                             g_0_xxxxzzz_0_yyyz_0, \
                             g_0_xxxxzzz_0_yyzz_0, \
                             g_0_xxxxzzz_0_yzzz_0, \
                             g_0_xxxxzzz_0_zzzz_0, \
                             g_0_xxxzzz_0_xxxz_0,  \
                             g_0_xxxzzz_0_xxxz_1,  \
                             g_0_xxxzzz_0_xxyz_0,  \
                             g_0_xxxzzz_0_xxyz_1,  \
                             g_0_xxxzzz_0_xxz_1,   \
                             g_0_xxxzzz_0_xxzz_0,  \
                             g_0_xxxzzz_0_xxzz_1,  \
                             g_0_xxxzzz_0_xyyz_0,  \
                             g_0_xxxzzz_0_xyyz_1,  \
                             g_0_xxxzzz_0_xyz_1,   \
                             g_0_xxxzzz_0_xyzz_0,  \
                             g_0_xxxzzz_0_xyzz_1,  \
                             g_0_xxxzzz_0_xzz_1,   \
                             g_0_xxxzzz_0_xzzz_0,  \
                             g_0_xxxzzz_0_xzzz_1,  \
                             g_0_xxxzzz_0_yyyy_0,  \
                             g_0_xxxzzz_0_yyyy_1,  \
                             g_0_xxxzzz_0_yyyz_0,  \
                             g_0_xxxzzz_0_yyyz_1,  \
                             g_0_xxxzzz_0_yyz_1,   \
                             g_0_xxxzzz_0_yyzz_0,  \
                             g_0_xxxzzz_0_yyzz_1,  \
                             g_0_xxxzzz_0_yzz_1,   \
                             g_0_xxxzzz_0_yzzz_0,  \
                             g_0_xxxzzz_0_yzzz_1,  \
                             g_0_xxxzzz_0_zzz_1,   \
                             g_0_xxxzzz_0_zzzz_0,  \
                             g_0_xxxzzz_0_zzzz_1,  \
                             g_0_xxzzz_0_xxxz_0,   \
                             g_0_xxzzz_0_xxxz_1,   \
                             g_0_xxzzz_0_xxyz_0,   \
                             g_0_xxzzz_0_xxyz_1,   \
                             g_0_xxzzz_0_xxzz_0,   \
                             g_0_xxzzz_0_xxzz_1,   \
                             g_0_xxzzz_0_xyyz_0,   \
                             g_0_xxzzz_0_xyyz_1,   \
                             g_0_xxzzz_0_xyzz_0,   \
                             g_0_xxzzz_0_xyzz_1,   \
                             g_0_xxzzz_0_xzzz_0,   \
                             g_0_xxzzz_0_xzzz_1,   \
                             g_0_xxzzz_0_yyyy_0,   \
                             g_0_xxzzz_0_yyyy_1,   \
                             g_0_xxzzz_0_yyyz_0,   \
                             g_0_xxzzz_0_yyyz_1,   \
                             g_0_xxzzz_0_yyzz_0,   \
                             g_0_xxzzz_0_yyzz_1,   \
                             g_0_xxzzz_0_yzzz_0,   \
                             g_0_xxzzz_0_yzzz_1,   \
                             g_0_xxzzz_0_zzzz_0,   \
                             g_0_xxzzz_0_zzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_xxxx_0[i] = 2.0 * g_0_xxxxz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxx_0[i] * pb_z +
                                  g_0_xxxxzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxy_0[i] = 2.0 * g_0_xxxxz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxy_0[i] * pb_z +
                                  g_0_xxxxzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxz_0[i] = 3.0 * g_0_xxzzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxz_0[i] * pb_x + g_0_xxxzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyy_0[i] = 2.0 * g_0_xxxxz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxyy_0[i] * pb_z +
                                  g_0_xxxxzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxyz_0[i] = 3.0 * g_0_xxzzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyz_0[i] * pb_x + g_0_xxxzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxzz_0[i] = 3.0 * g_0_xxzzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxzz_0[i] * pb_x + g_0_xxxzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyy_0[i] = 2.0 * g_0_xxxxz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xyyy_0[i] * pb_z +
                                  g_0_xxxxzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xyyz_0[i] = 3.0 * g_0_xxzzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxxzzz_0_xyyz_0[i] * pb_x + g_0_xxxzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyzz_0[i] = 3.0 * g_0_xxzzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxxzzz_0_xyzz_0[i] * pb_x + g_0_xxxzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xzzz_0[i] = 3.0 * g_0_xxzzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_xxxzzz_0_xzzz_0[i] * pb_x + g_0_xxxzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyy_0[i] = 3.0 * g_0_xxzzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyy_0[i] * pb_x +
                                  g_0_xxxzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyz_0[i] = 3.0 * g_0_xxzzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyz_0[i] * pb_x +
                                  g_0_xxxzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyzz_0[i] = 3.0 * g_0_xxzzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyzz_0[i] * pb_x +
                                  g_0_xxxzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yzzz_0[i] = 3.0 * g_0_xxzzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzzz_0[i] * pb_x +
                                  g_0_xxxzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_zzzz_0[i] = 3.0 * g_0_xxzzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzzz_0[i] * pb_x +
                                  g_0_xxxzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 150-165 components of targeted buffer : SKSG

    auto g_0_xxxyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 150);

    auto g_0_xxxyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 151);

    auto g_0_xxxyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 152);

    auto g_0_xxxyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 153);

    auto g_0_xxxyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 154);

    auto g_0_xxxyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 155);

    auto g_0_xxxyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 156);

    auto g_0_xxxyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 157);

    auto g_0_xxxyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 158);

    auto g_0_xxxyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 159);

    auto g_0_xxxyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 160);

    auto g_0_xxxyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 161);

    auto g_0_xxxyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 162);

    auto g_0_xxxyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 163);

    auto g_0_xxxyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 164);

#pragma omp simd aligned(g_0_xxxyy_0_xxxx_0,       \
                             g_0_xxxyy_0_xxxx_1,   \
                             g_0_xxxyy_0_xxxz_0,   \
                             g_0_xxxyy_0_xxxz_1,   \
                             g_0_xxxyy_0_xxzz_0,   \
                             g_0_xxxyy_0_xxzz_1,   \
                             g_0_xxxyy_0_xzzz_0,   \
                             g_0_xxxyy_0_xzzz_1,   \
                             g_0_xxxyyy_0_xxxx_0,  \
                             g_0_xxxyyy_0_xxxx_1,  \
                             g_0_xxxyyy_0_xxxz_0,  \
                             g_0_xxxyyy_0_xxxz_1,  \
                             g_0_xxxyyy_0_xxzz_0,  \
                             g_0_xxxyyy_0_xxzz_1,  \
                             g_0_xxxyyy_0_xzzz_0,  \
                             g_0_xxxyyy_0_xzzz_1,  \
                             g_0_xxxyyyy_0_xxxx_0, \
                             g_0_xxxyyyy_0_xxxy_0, \
                             g_0_xxxyyyy_0_xxxz_0, \
                             g_0_xxxyyyy_0_xxyy_0, \
                             g_0_xxxyyyy_0_xxyz_0, \
                             g_0_xxxyyyy_0_xxzz_0, \
                             g_0_xxxyyyy_0_xyyy_0, \
                             g_0_xxxyyyy_0_xyyz_0, \
                             g_0_xxxyyyy_0_xyzz_0, \
                             g_0_xxxyyyy_0_xzzz_0, \
                             g_0_xxxyyyy_0_yyyy_0, \
                             g_0_xxxyyyy_0_yyyz_0, \
                             g_0_xxxyyyy_0_yyzz_0, \
                             g_0_xxxyyyy_0_yzzz_0, \
                             g_0_xxxyyyy_0_zzzz_0, \
                             g_0_xxyyyy_0_xxxy_0,  \
                             g_0_xxyyyy_0_xxxy_1,  \
                             g_0_xxyyyy_0_xxy_1,   \
                             g_0_xxyyyy_0_xxyy_0,  \
                             g_0_xxyyyy_0_xxyy_1,  \
                             g_0_xxyyyy_0_xxyz_0,  \
                             g_0_xxyyyy_0_xxyz_1,  \
                             g_0_xxyyyy_0_xyy_1,   \
                             g_0_xxyyyy_0_xyyy_0,  \
                             g_0_xxyyyy_0_xyyy_1,  \
                             g_0_xxyyyy_0_xyyz_0,  \
                             g_0_xxyyyy_0_xyyz_1,  \
                             g_0_xxyyyy_0_xyz_1,   \
                             g_0_xxyyyy_0_xyzz_0,  \
                             g_0_xxyyyy_0_xyzz_1,  \
                             g_0_xxyyyy_0_yyy_1,   \
                             g_0_xxyyyy_0_yyyy_0,  \
                             g_0_xxyyyy_0_yyyy_1,  \
                             g_0_xxyyyy_0_yyyz_0,  \
                             g_0_xxyyyy_0_yyyz_1,  \
                             g_0_xxyyyy_0_yyz_1,   \
                             g_0_xxyyyy_0_yyzz_0,  \
                             g_0_xxyyyy_0_yyzz_1,  \
                             g_0_xxyyyy_0_yzz_1,   \
                             g_0_xxyyyy_0_yzzz_0,  \
                             g_0_xxyyyy_0_yzzz_1,  \
                             g_0_xxyyyy_0_zzzz_0,  \
                             g_0_xxyyyy_0_zzzz_1,  \
                             g_0_xyyyy_0_xxxy_0,   \
                             g_0_xyyyy_0_xxxy_1,   \
                             g_0_xyyyy_0_xxyy_0,   \
                             g_0_xyyyy_0_xxyy_1,   \
                             g_0_xyyyy_0_xxyz_0,   \
                             g_0_xyyyy_0_xxyz_1,   \
                             g_0_xyyyy_0_xyyy_0,   \
                             g_0_xyyyy_0_xyyy_1,   \
                             g_0_xyyyy_0_xyyz_0,   \
                             g_0_xyyyy_0_xyyz_1,   \
                             g_0_xyyyy_0_xyzz_0,   \
                             g_0_xyyyy_0_xyzz_1,   \
                             g_0_xyyyy_0_yyyy_0,   \
                             g_0_xyyyy_0_yyyy_1,   \
                             g_0_xyyyy_0_yyyz_0,   \
                             g_0_xyyyy_0_yyyz_1,   \
                             g_0_xyyyy_0_yyzz_0,   \
                             g_0_xyyyy_0_yyzz_1,   \
                             g_0_xyyyy_0_yzzz_0,   \
                             g_0_xyyyy_0_yzzz_1,   \
                             g_0_xyyyy_0_zzzz_0,   \
                             g_0_xyyyy_0_zzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_xxxx_0[i] = 3.0 * g_0_xxxyy_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxx_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxx_0[i] * pb_y +
                                  g_0_xxxyyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxy_0[i] = 2.0 * g_0_xyyyy_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxy_0[i] * pb_x + g_0_xxyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxz_0[i] = 3.0 * g_0_xxxyy_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxz_0[i] * pb_y +
                                  g_0_xxxyyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxyy_0[i] = 2.0 * g_0_xyyyy_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyy_0[i] * pb_x + g_0_xxyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyz_0[i] = 2.0 * g_0_xyyyy_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyz_0[i] * pb_x + g_0_xxyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxzz_0[i] = 3.0 * g_0_xxxyy_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxzz_0[i] * pb_y +
                                  g_0_xxxyyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xyyy_0[i] = 2.0 * g_0_xyyyy_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyy_1[i] * fi_abcd_0 +
                                  g_0_xxyyyy_0_xyyy_0[i] * pb_x + g_0_xxyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyz_0[i] = 2.0 * g_0_xyyyy_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxyyyy_0_xyyz_0[i] * pb_x + g_0_xxyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyzz_0[i] = 2.0 * g_0_xyyyy_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxyyyy_0_xyzz_0[i] * pb_x + g_0_xxyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xzzz_0[i] = 3.0 * g_0_xxxyy_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xzzz_0[i] * pb_y +
                                  g_0_xxxyyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_yyyy_0[i] = 2.0 * g_0_xyyyy_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyy_0[i] * pb_x +
                                  g_0_xxyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyz_0[i] = 2.0 * g_0_xyyyy_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyz_0[i] * pb_x +
                                  g_0_xxyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyzz_0[i] = 2.0 * g_0_xyyyy_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyzz_0[i] * pb_x +
                                  g_0_xxyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yzzz_0[i] = 2.0 * g_0_xyyyy_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzzz_0[i] * pb_x +
                                  g_0_xxyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_zzzz_0[i] = 2.0 * g_0_xyyyy_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_zzzz_0[i] * pb_x +
                                  g_0_xxyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 165-180 components of targeted buffer : SKSG

    auto g_0_xxxyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 165);

    auto g_0_xxxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 166);

    auto g_0_xxxyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 167);

    auto g_0_xxxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 168);

    auto g_0_xxxyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 169);

    auto g_0_xxxyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 170);

    auto g_0_xxxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 171);

    auto g_0_xxxyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 172);

    auto g_0_xxxyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 173);

    auto g_0_xxxyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 174);

    auto g_0_xxxyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 175);

    auto g_0_xxxyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 176);

    auto g_0_xxxyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 177);

    auto g_0_xxxyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 178);

    auto g_0_xxxyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 179);

#pragma omp simd aligned(g_0_xxxyyy_0_xxx_1,       \
                             g_0_xxxyyy_0_xxxx_0,  \
                             g_0_xxxyyy_0_xxxx_1,  \
                             g_0_xxxyyy_0_xxxy_0,  \
                             g_0_xxxyyy_0_xxxy_1,  \
                             g_0_xxxyyy_0_xxxz_0,  \
                             g_0_xxxyyy_0_xxxz_1,  \
                             g_0_xxxyyy_0_xxy_1,   \
                             g_0_xxxyyy_0_xxyy_0,  \
                             g_0_xxxyyy_0_xxyy_1,  \
                             g_0_xxxyyy_0_xxyz_0,  \
                             g_0_xxxyyy_0_xxyz_1,  \
                             g_0_xxxyyy_0_xxz_1,   \
                             g_0_xxxyyy_0_xxzz_0,  \
                             g_0_xxxyyy_0_xxzz_1,  \
                             g_0_xxxyyy_0_xyy_1,   \
                             g_0_xxxyyy_0_xyyy_0,  \
                             g_0_xxxyyy_0_xyyy_1,  \
                             g_0_xxxyyy_0_xyyz_0,  \
                             g_0_xxxyyy_0_xyyz_1,  \
                             g_0_xxxyyy_0_xyz_1,   \
                             g_0_xxxyyy_0_xyzz_0,  \
                             g_0_xxxyyy_0_xyzz_1,  \
                             g_0_xxxyyy_0_xzz_1,   \
                             g_0_xxxyyy_0_xzzz_0,  \
                             g_0_xxxyyy_0_xzzz_1,  \
                             g_0_xxxyyy_0_yyy_1,   \
                             g_0_xxxyyy_0_yyyy_0,  \
                             g_0_xxxyyy_0_yyyy_1,  \
                             g_0_xxxyyy_0_yyyz_0,  \
                             g_0_xxxyyy_0_yyyz_1,  \
                             g_0_xxxyyy_0_yyz_1,   \
                             g_0_xxxyyy_0_yyzz_0,  \
                             g_0_xxxyyy_0_yyzz_1,  \
                             g_0_xxxyyy_0_yzz_1,   \
                             g_0_xxxyyy_0_yzzz_0,  \
                             g_0_xxxyyy_0_yzzz_1,  \
                             g_0_xxxyyy_0_zzz_1,   \
                             g_0_xxxyyy_0_zzzz_0,  \
                             g_0_xxxyyy_0_zzzz_1,  \
                             g_0_xxxyyyz_0_xxxx_0, \
                             g_0_xxxyyyz_0_xxxy_0, \
                             g_0_xxxyyyz_0_xxxz_0, \
                             g_0_xxxyyyz_0_xxyy_0, \
                             g_0_xxxyyyz_0_xxyz_0, \
                             g_0_xxxyyyz_0_xxzz_0, \
                             g_0_xxxyyyz_0_xyyy_0, \
                             g_0_xxxyyyz_0_xyyz_0, \
                             g_0_xxxyyyz_0_xyzz_0, \
                             g_0_xxxyyyz_0_xzzz_0, \
                             g_0_xxxyyyz_0_yyyy_0, \
                             g_0_xxxyyyz_0_yyyz_0, \
                             g_0_xxxyyyz_0_yyzz_0, \
                             g_0_xxxyyyz_0_yzzz_0, \
                             g_0_xxxyyyz_0_zzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_xxxx_0[i] = g_0_xxxyyy_0_xxxx_0[i] * pb_z + g_0_xxxyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxy_0[i] = g_0_xxxyyy_0_xxxy_0[i] * pb_z + g_0_xxxyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxz_0[i] = g_0_xxxyyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxz_0[i] * pb_z + g_0_xxxyyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyy_0[i] = g_0_xxxyyy_0_xxyy_0[i] * pb_z + g_0_xxxyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyz_0[i] = g_0_xxxyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyz_0[i] * pb_z + g_0_xxxyyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxzz_0[i] = 2.0 * g_0_xxxyyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxzz_0[i] * pb_z + g_0_xxxyyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyy_0[i] = g_0_xxxyyy_0_xyyy_0[i] * pb_z + g_0_xxxyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyz_0[i] = g_0_xxxyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyz_0[i] * pb_z + g_0_xxxyyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyzz_0[i] = 2.0 * g_0_xxxyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzz_0[i] * pb_z + g_0_xxxyyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xzzz_0[i] = 3.0 * g_0_xxxyyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xzzz_0[i] * pb_z + g_0_xxxyyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyy_0[i] = g_0_xxxyyy_0_yyyy_0[i] * pb_z + g_0_xxxyyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyz_0[i] = g_0_xxxyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyz_0[i] * pb_z + g_0_xxxyyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyzz_0[i] = 2.0 * g_0_xxxyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyzz_0[i] * pb_z + g_0_xxxyyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yzzz_0[i] = 3.0 * g_0_xxxyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yzzz_0[i] * pb_z + g_0_xxxyyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_zzzz_0[i] = 4.0 * g_0_xxxyyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_zzzz_0[i] * pb_z + g_0_xxxyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 180-195 components of targeted buffer : SKSG

    auto g_0_xxxyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 180);

    auto g_0_xxxyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 181);

    auto g_0_xxxyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 182);

    auto g_0_xxxyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 183);

    auto g_0_xxxyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 184);

    auto g_0_xxxyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 185);

    auto g_0_xxxyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 186);

    auto g_0_xxxyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 187);

    auto g_0_xxxyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 188);

    auto g_0_xxxyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 189);

    auto g_0_xxxyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 190);

    auto g_0_xxxyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 191);

    auto g_0_xxxyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 192);

    auto g_0_xxxyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 193);

    auto g_0_xxxyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 194);

#pragma omp simd aligned(g_0_xxxyy_0_xxxy_0,       \
                             g_0_xxxyy_0_xxxy_1,   \
                             g_0_xxxyy_0_xxyy_0,   \
                             g_0_xxxyy_0_xxyy_1,   \
                             g_0_xxxyy_0_xyyy_0,   \
                             g_0_xxxyy_0_xyyy_1,   \
                             g_0_xxxyyz_0_xxxy_0,  \
                             g_0_xxxyyz_0_xxxy_1,  \
                             g_0_xxxyyz_0_xxyy_0,  \
                             g_0_xxxyyz_0_xxyy_1,  \
                             g_0_xxxyyz_0_xyyy_0,  \
                             g_0_xxxyyz_0_xyyy_1,  \
                             g_0_xxxyyzz_0_xxxx_0, \
                             g_0_xxxyyzz_0_xxxy_0, \
                             g_0_xxxyyzz_0_xxxz_0, \
                             g_0_xxxyyzz_0_xxyy_0, \
                             g_0_xxxyyzz_0_xxyz_0, \
                             g_0_xxxyyzz_0_xxzz_0, \
                             g_0_xxxyyzz_0_xyyy_0, \
                             g_0_xxxyyzz_0_xyyz_0, \
                             g_0_xxxyyzz_0_xyzz_0, \
                             g_0_xxxyyzz_0_xzzz_0, \
                             g_0_xxxyyzz_0_yyyy_0, \
                             g_0_xxxyyzz_0_yyyz_0, \
                             g_0_xxxyyzz_0_yyzz_0, \
                             g_0_xxxyyzz_0_yzzz_0, \
                             g_0_xxxyyzz_0_zzzz_0, \
                             g_0_xxxyzz_0_xxxx_0,  \
                             g_0_xxxyzz_0_xxxx_1,  \
                             g_0_xxxyzz_0_xxxz_0,  \
                             g_0_xxxyzz_0_xxxz_1,  \
                             g_0_xxxyzz_0_xxzz_0,  \
                             g_0_xxxyzz_0_xxzz_1,  \
                             g_0_xxxyzz_0_xzzz_0,  \
                             g_0_xxxyzz_0_xzzz_1,  \
                             g_0_xxxzz_0_xxxx_0,   \
                             g_0_xxxzz_0_xxxx_1,   \
                             g_0_xxxzz_0_xxxz_0,   \
                             g_0_xxxzz_0_xxxz_1,   \
                             g_0_xxxzz_0_xxzz_0,   \
                             g_0_xxxzz_0_xxzz_1,   \
                             g_0_xxxzz_0_xzzz_0,   \
                             g_0_xxxzz_0_xzzz_1,   \
                             g_0_xxyyzz_0_xxyz_0,  \
                             g_0_xxyyzz_0_xxyz_1,  \
                             g_0_xxyyzz_0_xyyz_0,  \
                             g_0_xxyyzz_0_xyyz_1,  \
                             g_0_xxyyzz_0_xyz_1,   \
                             g_0_xxyyzz_0_xyzz_0,  \
                             g_0_xxyyzz_0_xyzz_1,  \
                             g_0_xxyyzz_0_yyyy_0,  \
                             g_0_xxyyzz_0_yyyy_1,  \
                             g_0_xxyyzz_0_yyyz_0,  \
                             g_0_xxyyzz_0_yyyz_1,  \
                             g_0_xxyyzz_0_yyz_1,   \
                             g_0_xxyyzz_0_yyzz_0,  \
                             g_0_xxyyzz_0_yyzz_1,  \
                             g_0_xxyyzz_0_yzz_1,   \
                             g_0_xxyyzz_0_yzzz_0,  \
                             g_0_xxyyzz_0_yzzz_1,  \
                             g_0_xxyyzz_0_zzzz_0,  \
                             g_0_xxyyzz_0_zzzz_1,  \
                             g_0_xyyzz_0_xxyz_0,   \
                             g_0_xyyzz_0_xxyz_1,   \
                             g_0_xyyzz_0_xyyz_0,   \
                             g_0_xyyzz_0_xyyz_1,   \
                             g_0_xyyzz_0_xyzz_0,   \
                             g_0_xyyzz_0_xyzz_1,   \
                             g_0_xyyzz_0_yyyy_0,   \
                             g_0_xyyzz_0_yyyy_1,   \
                             g_0_xyyzz_0_yyyz_0,   \
                             g_0_xyyzz_0_yyyz_1,   \
                             g_0_xyyzz_0_yyzz_0,   \
                             g_0_xyyzz_0_yyzz_1,   \
                             g_0_xyyzz_0_yzzz_0,   \
                             g_0_xyyzz_0_yzzz_1,   \
                             g_0_xyyzz_0_zzzz_0,   \
                             g_0_xyyzz_0_zzzz_1,   \
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

        g_0_xxxyyzz_0_xxxx_0[i] =
            g_0_xxxzz_0_xxxx_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxx_0[i] * pb_y + g_0_xxxyzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxy_0[i] =
            g_0_xxxyy_0_xxxy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxy_0[i] * pb_z + g_0_xxxyyz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxz_0[i] =
            g_0_xxxzz_0_xxxz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxz_0[i] * pb_y + g_0_xxxyzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxyy_0[i] =
            g_0_xxxyy_0_xxyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxyy_0[i] * pb_z + g_0_xxxyyz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxyz_0[i] = 2.0 * g_0_xyyzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyz_0[i] * pb_x + g_0_xxyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxzz_0[i] =
            g_0_xxxzz_0_xxzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxzz_0[i] * pb_y + g_0_xxxyzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xyyy_0[i] =
            g_0_xxxyy_0_xyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xyyy_0[i] * pb_z + g_0_xxxyyz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xyyz_0[i] = 2.0 * g_0_xyyzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxyyzz_0_xyyz_0[i] * pb_x + g_0_xxyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyzz_0[i] = 2.0 * g_0_xyyzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxyyzz_0_xyzz_0[i] * pb_x + g_0_xxyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xzzz_0[i] =
            g_0_xxxzz_0_xzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xzzz_0[i] * pb_y + g_0_xxxyzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_yyyy_0[i] = 2.0 * g_0_xyyzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyy_0[i] * pb_x +
                                  g_0_xxyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyz_0[i] = 2.0 * g_0_xyyzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyz_0[i] * pb_x +
                                  g_0_xxyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyzz_0[i] = 2.0 * g_0_xyyzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyzz_0[i] * pb_x +
                                  g_0_xxyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yzzz_0[i] = 2.0 * g_0_xyyzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzzz_0[i] * pb_x +
                                  g_0_xxyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_zzzz_0[i] = 2.0 * g_0_xyyzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_zzzz_0[i] * pb_x +
                                  g_0_xxyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 195-210 components of targeted buffer : SKSG

    auto g_0_xxxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 195);

    auto g_0_xxxyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 196);

    auto g_0_xxxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 197);

    auto g_0_xxxyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 198);

    auto g_0_xxxyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 199);

    auto g_0_xxxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 200);

    auto g_0_xxxyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 201);

    auto g_0_xxxyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 202);

    auto g_0_xxxyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 203);

    auto g_0_xxxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 204);

    auto g_0_xxxyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 205);

    auto g_0_xxxyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 206);

    auto g_0_xxxyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 207);

    auto g_0_xxxyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 208);

    auto g_0_xxxyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 209);

#pragma omp simd aligned(g_0_xxxyzzz_0_xxxx_0,     \
                             g_0_xxxyzzz_0_xxxy_0, \
                             g_0_xxxyzzz_0_xxxz_0, \
                             g_0_xxxyzzz_0_xxyy_0, \
                             g_0_xxxyzzz_0_xxyz_0, \
                             g_0_xxxyzzz_0_xxzz_0, \
                             g_0_xxxyzzz_0_xyyy_0, \
                             g_0_xxxyzzz_0_xyyz_0, \
                             g_0_xxxyzzz_0_xyzz_0, \
                             g_0_xxxyzzz_0_xzzz_0, \
                             g_0_xxxyzzz_0_yyyy_0, \
                             g_0_xxxyzzz_0_yyyz_0, \
                             g_0_xxxyzzz_0_yyzz_0, \
                             g_0_xxxyzzz_0_yzzz_0, \
                             g_0_xxxyzzz_0_zzzz_0, \
                             g_0_xxxzzz_0_xxx_1,   \
                             g_0_xxxzzz_0_xxxx_0,  \
                             g_0_xxxzzz_0_xxxx_1,  \
                             g_0_xxxzzz_0_xxxy_0,  \
                             g_0_xxxzzz_0_xxxy_1,  \
                             g_0_xxxzzz_0_xxxz_0,  \
                             g_0_xxxzzz_0_xxxz_1,  \
                             g_0_xxxzzz_0_xxy_1,   \
                             g_0_xxxzzz_0_xxyy_0,  \
                             g_0_xxxzzz_0_xxyy_1,  \
                             g_0_xxxzzz_0_xxyz_0,  \
                             g_0_xxxzzz_0_xxyz_1,  \
                             g_0_xxxzzz_0_xxz_1,   \
                             g_0_xxxzzz_0_xxzz_0,  \
                             g_0_xxxzzz_0_xxzz_1,  \
                             g_0_xxxzzz_0_xyy_1,   \
                             g_0_xxxzzz_0_xyyy_0,  \
                             g_0_xxxzzz_0_xyyy_1,  \
                             g_0_xxxzzz_0_xyyz_0,  \
                             g_0_xxxzzz_0_xyyz_1,  \
                             g_0_xxxzzz_0_xyz_1,   \
                             g_0_xxxzzz_0_xyzz_0,  \
                             g_0_xxxzzz_0_xyzz_1,  \
                             g_0_xxxzzz_0_xzz_1,   \
                             g_0_xxxzzz_0_xzzz_0,  \
                             g_0_xxxzzz_0_xzzz_1,  \
                             g_0_xxxzzz_0_yyy_1,   \
                             g_0_xxxzzz_0_yyyy_0,  \
                             g_0_xxxzzz_0_yyyy_1,  \
                             g_0_xxxzzz_0_yyyz_0,  \
                             g_0_xxxzzz_0_yyyz_1,  \
                             g_0_xxxzzz_0_yyz_1,   \
                             g_0_xxxzzz_0_yyzz_0,  \
                             g_0_xxxzzz_0_yyzz_1,  \
                             g_0_xxxzzz_0_yzz_1,   \
                             g_0_xxxzzz_0_yzzz_0,  \
                             g_0_xxxzzz_0_yzzz_1,  \
                             g_0_xxxzzz_0_zzz_1,   \
                             g_0_xxxzzz_0_zzzz_0,  \
                             g_0_xxxzzz_0_zzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_xxxx_0[i] = g_0_xxxzzz_0_xxxx_0[i] * pb_y + g_0_xxxzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxy_0[i] = g_0_xxxzzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxy_0[i] * pb_y + g_0_xxxzzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxz_0[i] = g_0_xxxzzz_0_xxxz_0[i] * pb_y + g_0_xxxzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyy_0[i] = 2.0 * g_0_xxxzzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyy_0[i] * pb_y + g_0_xxxzzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyz_0[i] = g_0_xxxzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyz_0[i] * pb_y + g_0_xxxzzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxzz_0[i] = g_0_xxxzzz_0_xxzz_0[i] * pb_y + g_0_xxxzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyy_0[i] = 3.0 * g_0_xxxzzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyy_0[i] * pb_y + g_0_xxxzzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyz_0[i] = 2.0 * g_0_xxxzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyz_0[i] * pb_y + g_0_xxxzzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyzz_0[i] = g_0_xxxzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzz_0[i] * pb_y + g_0_xxxzzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xzzz_0[i] = g_0_xxxzzz_0_xzzz_0[i] * pb_y + g_0_xxxzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyy_0[i] = 4.0 * g_0_xxxzzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyy_0[i] * pb_y + g_0_xxxzzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyz_0[i] = 3.0 * g_0_xxxzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyz_0[i] * pb_y + g_0_xxxzzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyzz_0[i] = 2.0 * g_0_xxxzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyzz_0[i] * pb_y + g_0_xxxzzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yzzz_0[i] = g_0_xxxzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yzzz_0[i] * pb_y + g_0_xxxzzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_zzzz_0[i] = g_0_xxxzzz_0_zzzz_0[i] * pb_y + g_0_xxxzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 210-225 components of targeted buffer : SKSG

    auto g_0_xxxzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 210);

    auto g_0_xxxzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 211);

    auto g_0_xxxzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 212);

    auto g_0_xxxzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 213);

    auto g_0_xxxzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 214);

    auto g_0_xxxzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 215);

    auto g_0_xxxzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 216);

    auto g_0_xxxzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 217);

    auto g_0_xxxzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 218);

    auto g_0_xxxzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 219);

    auto g_0_xxxzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 220);

    auto g_0_xxxzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 221);

    auto g_0_xxxzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 222);

    auto g_0_xxxzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 223);

    auto g_0_xxxzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 224);

#pragma omp simd aligned(g_0_xxxzz_0_xxxx_0,       \
                             g_0_xxxzz_0_xxxx_1,   \
                             g_0_xxxzz_0_xxxy_0,   \
                             g_0_xxxzz_0_xxxy_1,   \
                             g_0_xxxzz_0_xxyy_0,   \
                             g_0_xxxzz_0_xxyy_1,   \
                             g_0_xxxzz_0_xyyy_0,   \
                             g_0_xxxzz_0_xyyy_1,   \
                             g_0_xxxzzz_0_xxxx_0,  \
                             g_0_xxxzzz_0_xxxx_1,  \
                             g_0_xxxzzz_0_xxxy_0,  \
                             g_0_xxxzzz_0_xxxy_1,  \
                             g_0_xxxzzz_0_xxyy_0,  \
                             g_0_xxxzzz_0_xxyy_1,  \
                             g_0_xxxzzz_0_xyyy_0,  \
                             g_0_xxxzzz_0_xyyy_1,  \
                             g_0_xxxzzzz_0_xxxx_0, \
                             g_0_xxxzzzz_0_xxxy_0, \
                             g_0_xxxzzzz_0_xxxz_0, \
                             g_0_xxxzzzz_0_xxyy_0, \
                             g_0_xxxzzzz_0_xxyz_0, \
                             g_0_xxxzzzz_0_xxzz_0, \
                             g_0_xxxzzzz_0_xyyy_0, \
                             g_0_xxxzzzz_0_xyyz_0, \
                             g_0_xxxzzzz_0_xyzz_0, \
                             g_0_xxxzzzz_0_xzzz_0, \
                             g_0_xxxzzzz_0_yyyy_0, \
                             g_0_xxxzzzz_0_yyyz_0, \
                             g_0_xxxzzzz_0_yyzz_0, \
                             g_0_xxxzzzz_0_yzzz_0, \
                             g_0_xxxzzzz_0_zzzz_0, \
                             g_0_xxzzzz_0_xxxz_0,  \
                             g_0_xxzzzz_0_xxxz_1,  \
                             g_0_xxzzzz_0_xxyz_0,  \
                             g_0_xxzzzz_0_xxyz_1,  \
                             g_0_xxzzzz_0_xxz_1,   \
                             g_0_xxzzzz_0_xxzz_0,  \
                             g_0_xxzzzz_0_xxzz_1,  \
                             g_0_xxzzzz_0_xyyz_0,  \
                             g_0_xxzzzz_0_xyyz_1,  \
                             g_0_xxzzzz_0_xyz_1,   \
                             g_0_xxzzzz_0_xyzz_0,  \
                             g_0_xxzzzz_0_xyzz_1,  \
                             g_0_xxzzzz_0_xzz_1,   \
                             g_0_xxzzzz_0_xzzz_0,  \
                             g_0_xxzzzz_0_xzzz_1,  \
                             g_0_xxzzzz_0_yyyy_0,  \
                             g_0_xxzzzz_0_yyyy_1,  \
                             g_0_xxzzzz_0_yyyz_0,  \
                             g_0_xxzzzz_0_yyyz_1,  \
                             g_0_xxzzzz_0_yyz_1,   \
                             g_0_xxzzzz_0_yyzz_0,  \
                             g_0_xxzzzz_0_yyzz_1,  \
                             g_0_xxzzzz_0_yzz_1,   \
                             g_0_xxzzzz_0_yzzz_0,  \
                             g_0_xxzzzz_0_yzzz_1,  \
                             g_0_xxzzzz_0_zzz_1,   \
                             g_0_xxzzzz_0_zzzz_0,  \
                             g_0_xxzzzz_0_zzzz_1,  \
                             g_0_xzzzz_0_xxxz_0,   \
                             g_0_xzzzz_0_xxxz_1,   \
                             g_0_xzzzz_0_xxyz_0,   \
                             g_0_xzzzz_0_xxyz_1,   \
                             g_0_xzzzz_0_xxzz_0,   \
                             g_0_xzzzz_0_xxzz_1,   \
                             g_0_xzzzz_0_xyyz_0,   \
                             g_0_xzzzz_0_xyyz_1,   \
                             g_0_xzzzz_0_xyzz_0,   \
                             g_0_xzzzz_0_xyzz_1,   \
                             g_0_xzzzz_0_xzzz_0,   \
                             g_0_xzzzz_0_xzzz_1,   \
                             g_0_xzzzz_0_yyyy_0,   \
                             g_0_xzzzz_0_yyyy_1,   \
                             g_0_xzzzz_0_yyyz_0,   \
                             g_0_xzzzz_0_yyyz_1,   \
                             g_0_xzzzz_0_yyzz_0,   \
                             g_0_xzzzz_0_yyzz_1,   \
                             g_0_xzzzz_0_yzzz_0,   \
                             g_0_xzzzz_0_yzzz_1,   \
                             g_0_xzzzz_0_zzzz_0,   \
                             g_0_xzzzz_0_zzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_xxxx_0[i] = 3.0 * g_0_xxxzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxx_0[i] * pb_z +
                                  g_0_xxxzzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxy_0[i] = 3.0 * g_0_xxxzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxy_0[i] * pb_z +
                                  g_0_xxxzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxz_0[i] = 2.0 * g_0_xzzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxz_0[i] * pb_x + g_0_xxzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyy_0[i] = 3.0 * g_0_xxxzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxyy_0[i] * pb_z +
                                  g_0_xxxzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxyz_0[i] = 2.0 * g_0_xzzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyz_0[i] * pb_x + g_0_xxzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxzz_0[i] = 2.0 * g_0_xzzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxzz_0[i] * pb_x + g_0_xxzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyy_0[i] = 3.0 * g_0_xxxzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xyyy_0[i] * pb_z +
                                  g_0_xxxzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xyyz_0[i] = 2.0 * g_0_xzzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xxzzzz_0_xyyz_0[i] * pb_x + g_0_xxzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyzz_0[i] = 2.0 * g_0_xzzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xxzzzz_0_xyzz_0[i] * pb_x + g_0_xxzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xzzz_0[i] = 2.0 * g_0_xzzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_xxzzzz_0_xzzz_0[i] * pb_x + g_0_xxzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyy_0[i] = 2.0 * g_0_xzzzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyy_0[i] * pb_x +
                                  g_0_xxzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyz_0[i] = 2.0 * g_0_xzzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyz_0[i] * pb_x +
                                  g_0_xxzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyzz_0[i] = 2.0 * g_0_xzzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyzz_0[i] * pb_x +
                                  g_0_xxzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yzzz_0[i] = 2.0 * g_0_xzzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzzz_0[i] * pb_x +
                                  g_0_xxzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_zzzz_0[i] = 2.0 * g_0_xzzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzzz_0[i] * pb_x +
                                  g_0_xxzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 225-240 components of targeted buffer : SKSG

    auto g_0_xxyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 225);

    auto g_0_xxyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 226);

    auto g_0_xxyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 227);

    auto g_0_xxyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 228);

    auto g_0_xxyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 229);

    auto g_0_xxyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 230);

    auto g_0_xxyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 231);

    auto g_0_xxyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 232);

    auto g_0_xxyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 233);

    auto g_0_xxyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 234);

    auto g_0_xxyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 235);

    auto g_0_xxyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 236);

    auto g_0_xxyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 237);

    auto g_0_xxyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 238);

    auto g_0_xxyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 239);

#pragma omp simd aligned(g_0_xxyyy_0_xxxx_0,       \
                             g_0_xxyyy_0_xxxx_1,   \
                             g_0_xxyyy_0_xxxz_0,   \
                             g_0_xxyyy_0_xxxz_1,   \
                             g_0_xxyyy_0_xxzz_0,   \
                             g_0_xxyyy_0_xxzz_1,   \
                             g_0_xxyyy_0_xzzz_0,   \
                             g_0_xxyyy_0_xzzz_1,   \
                             g_0_xxyyyy_0_xxxx_0,  \
                             g_0_xxyyyy_0_xxxx_1,  \
                             g_0_xxyyyy_0_xxxz_0,  \
                             g_0_xxyyyy_0_xxxz_1,  \
                             g_0_xxyyyy_0_xxzz_0,  \
                             g_0_xxyyyy_0_xxzz_1,  \
                             g_0_xxyyyy_0_xzzz_0,  \
                             g_0_xxyyyy_0_xzzz_1,  \
                             g_0_xxyyyyy_0_xxxx_0, \
                             g_0_xxyyyyy_0_xxxy_0, \
                             g_0_xxyyyyy_0_xxxz_0, \
                             g_0_xxyyyyy_0_xxyy_0, \
                             g_0_xxyyyyy_0_xxyz_0, \
                             g_0_xxyyyyy_0_xxzz_0, \
                             g_0_xxyyyyy_0_xyyy_0, \
                             g_0_xxyyyyy_0_xyyz_0, \
                             g_0_xxyyyyy_0_xyzz_0, \
                             g_0_xxyyyyy_0_xzzz_0, \
                             g_0_xxyyyyy_0_yyyy_0, \
                             g_0_xxyyyyy_0_yyyz_0, \
                             g_0_xxyyyyy_0_yyzz_0, \
                             g_0_xxyyyyy_0_yzzz_0, \
                             g_0_xxyyyyy_0_zzzz_0, \
                             g_0_xyyyyy_0_xxxy_0,  \
                             g_0_xyyyyy_0_xxxy_1,  \
                             g_0_xyyyyy_0_xxy_1,   \
                             g_0_xyyyyy_0_xxyy_0,  \
                             g_0_xyyyyy_0_xxyy_1,  \
                             g_0_xyyyyy_0_xxyz_0,  \
                             g_0_xyyyyy_0_xxyz_1,  \
                             g_0_xyyyyy_0_xyy_1,   \
                             g_0_xyyyyy_0_xyyy_0,  \
                             g_0_xyyyyy_0_xyyy_1,  \
                             g_0_xyyyyy_0_xyyz_0,  \
                             g_0_xyyyyy_0_xyyz_1,  \
                             g_0_xyyyyy_0_xyz_1,   \
                             g_0_xyyyyy_0_xyzz_0,  \
                             g_0_xyyyyy_0_xyzz_1,  \
                             g_0_xyyyyy_0_yyy_1,   \
                             g_0_xyyyyy_0_yyyy_0,  \
                             g_0_xyyyyy_0_yyyy_1,  \
                             g_0_xyyyyy_0_yyyz_0,  \
                             g_0_xyyyyy_0_yyyz_1,  \
                             g_0_xyyyyy_0_yyz_1,   \
                             g_0_xyyyyy_0_yyzz_0,  \
                             g_0_xyyyyy_0_yyzz_1,  \
                             g_0_xyyyyy_0_yzz_1,   \
                             g_0_xyyyyy_0_yzzz_0,  \
                             g_0_xyyyyy_0_yzzz_1,  \
                             g_0_xyyyyy_0_zzzz_0,  \
                             g_0_xyyyyy_0_zzzz_1,  \
                             g_0_yyyyy_0_xxxy_0,   \
                             g_0_yyyyy_0_xxxy_1,   \
                             g_0_yyyyy_0_xxyy_0,   \
                             g_0_yyyyy_0_xxyy_1,   \
                             g_0_yyyyy_0_xxyz_0,   \
                             g_0_yyyyy_0_xxyz_1,   \
                             g_0_yyyyy_0_xyyy_0,   \
                             g_0_yyyyy_0_xyyy_1,   \
                             g_0_yyyyy_0_xyyz_0,   \
                             g_0_yyyyy_0_xyyz_1,   \
                             g_0_yyyyy_0_xyzz_0,   \
                             g_0_yyyyy_0_xyzz_1,   \
                             g_0_yyyyy_0_yyyy_0,   \
                             g_0_yyyyy_0_yyyy_1,   \
                             g_0_yyyyy_0_yyyz_0,   \
                             g_0_yyyyy_0_yyyz_1,   \
                             g_0_yyyyy_0_yyzz_0,   \
                             g_0_yyyyy_0_yyzz_1,   \
                             g_0_yyyyy_0_yzzz_0,   \
                             g_0_yyyyy_0_yzzz_1,   \
                             g_0_yyyyy_0_zzzz_0,   \
                             g_0_yyyyy_0_zzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_xxxx_0[i] = 4.0 * g_0_xxyyy_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxx_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxx_0[i] * pb_y +
                                  g_0_xxyyyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxy_0[i] = g_0_yyyyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxy_1[i] * fi_abcd_0 +
                                  g_0_xyyyyy_0_xxxy_0[i] * pb_x + g_0_xyyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxz_0[i] = 4.0 * g_0_xxyyy_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxz_0[i] * pb_y +
                                  g_0_xxyyyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxyy_0[i] = g_0_yyyyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyy_1[i] * fi_abcd_0 +
                                  g_0_xyyyyy_0_xxyy_0[i] * pb_x + g_0_xyyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyz_0[i] = g_0_yyyyy_0_xxyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyyy_0_xxyz_0[i] * pb_x + g_0_xyyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxzz_0[i] = 4.0 * g_0_xxyyy_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxzz_0[i] * pb_y +
                                  g_0_xxyyyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xyyy_0[i] = g_0_yyyyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyy_1[i] * fi_abcd_0 +
                                  g_0_xyyyyy_0_xyyy_0[i] * pb_x + g_0_xyyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyz_0[i] = g_0_yyyyy_0_xyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyyy_0_xyyz_0[i] * pb_x + g_0_xyyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyzz_0[i] = g_0_yyyyy_0_xyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xyyyyy_0_xyzz_0[i] * pb_x + g_0_xyyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xzzz_0[i] = 4.0 * g_0_xxyyy_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xzzz_0[i] * pb_y +
                                  g_0_xxyyyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_yyyy_0[i] =
            g_0_yyyyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyy_0[i] * pb_x + g_0_xyyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyz_0[i] =
            g_0_yyyyy_0_yyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyz_0[i] * pb_x + g_0_xyyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyzz_0[i] =
            g_0_yyyyy_0_yyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzz_0[i] * pb_x + g_0_xyyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yzzz_0[i] =
            g_0_yyyyy_0_yzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzz_0[i] * pb_x + g_0_xyyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_zzzz_0[i] =
            g_0_yyyyy_0_zzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_zzzz_0[i] * pb_x + g_0_xyyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 240-255 components of targeted buffer : SKSG

    auto g_0_xxyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 240);

    auto g_0_xxyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 241);

    auto g_0_xxyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 242);

    auto g_0_xxyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 243);

    auto g_0_xxyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 244);

    auto g_0_xxyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 245);

    auto g_0_xxyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 246);

    auto g_0_xxyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 247);

    auto g_0_xxyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 248);

    auto g_0_xxyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 249);

    auto g_0_xxyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 250);

    auto g_0_xxyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 251);

    auto g_0_xxyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 252);

    auto g_0_xxyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 253);

    auto g_0_xxyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 254);

#pragma omp simd aligned(g_0_xxyyyy_0_xxx_1,       \
                             g_0_xxyyyy_0_xxxx_0,  \
                             g_0_xxyyyy_0_xxxx_1,  \
                             g_0_xxyyyy_0_xxxy_0,  \
                             g_0_xxyyyy_0_xxxy_1,  \
                             g_0_xxyyyy_0_xxxz_0,  \
                             g_0_xxyyyy_0_xxxz_1,  \
                             g_0_xxyyyy_0_xxy_1,   \
                             g_0_xxyyyy_0_xxyy_0,  \
                             g_0_xxyyyy_0_xxyy_1,  \
                             g_0_xxyyyy_0_xxyz_0,  \
                             g_0_xxyyyy_0_xxyz_1,  \
                             g_0_xxyyyy_0_xxz_1,   \
                             g_0_xxyyyy_0_xxzz_0,  \
                             g_0_xxyyyy_0_xxzz_1,  \
                             g_0_xxyyyy_0_xyy_1,   \
                             g_0_xxyyyy_0_xyyy_0,  \
                             g_0_xxyyyy_0_xyyy_1,  \
                             g_0_xxyyyy_0_xyyz_0,  \
                             g_0_xxyyyy_0_xyyz_1,  \
                             g_0_xxyyyy_0_xyz_1,   \
                             g_0_xxyyyy_0_xyzz_0,  \
                             g_0_xxyyyy_0_xyzz_1,  \
                             g_0_xxyyyy_0_xzz_1,   \
                             g_0_xxyyyy_0_xzzz_0,  \
                             g_0_xxyyyy_0_xzzz_1,  \
                             g_0_xxyyyy_0_yyy_1,   \
                             g_0_xxyyyy_0_yyyy_0,  \
                             g_0_xxyyyy_0_yyyy_1,  \
                             g_0_xxyyyy_0_yyyz_0,  \
                             g_0_xxyyyy_0_yyyz_1,  \
                             g_0_xxyyyy_0_yyz_1,   \
                             g_0_xxyyyy_0_yyzz_0,  \
                             g_0_xxyyyy_0_yyzz_1,  \
                             g_0_xxyyyy_0_yzz_1,   \
                             g_0_xxyyyy_0_yzzz_0,  \
                             g_0_xxyyyy_0_yzzz_1,  \
                             g_0_xxyyyy_0_zzz_1,   \
                             g_0_xxyyyy_0_zzzz_0,  \
                             g_0_xxyyyy_0_zzzz_1,  \
                             g_0_xxyyyyz_0_xxxx_0, \
                             g_0_xxyyyyz_0_xxxy_0, \
                             g_0_xxyyyyz_0_xxxz_0, \
                             g_0_xxyyyyz_0_xxyy_0, \
                             g_0_xxyyyyz_0_xxyz_0, \
                             g_0_xxyyyyz_0_xxzz_0, \
                             g_0_xxyyyyz_0_xyyy_0, \
                             g_0_xxyyyyz_0_xyyz_0, \
                             g_0_xxyyyyz_0_xyzz_0, \
                             g_0_xxyyyyz_0_xzzz_0, \
                             g_0_xxyyyyz_0_yyyy_0, \
                             g_0_xxyyyyz_0_yyyz_0, \
                             g_0_xxyyyyz_0_yyzz_0, \
                             g_0_xxyyyyz_0_yzzz_0, \
                             g_0_xxyyyyz_0_zzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_xxxx_0[i] = g_0_xxyyyy_0_xxxx_0[i] * pb_z + g_0_xxyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxy_0[i] = g_0_xxyyyy_0_xxxy_0[i] * pb_z + g_0_xxyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxz_0[i] = g_0_xxyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxz_0[i] * pb_z + g_0_xxyyyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyy_0[i] = g_0_xxyyyy_0_xxyy_0[i] * pb_z + g_0_xxyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyz_0[i] = g_0_xxyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyz_0[i] * pb_z + g_0_xxyyyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxzz_0[i] = 2.0 * g_0_xxyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxzz_0[i] * pb_z + g_0_xxyyyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyy_0[i] = g_0_xxyyyy_0_xyyy_0[i] * pb_z + g_0_xxyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyz_0[i] = g_0_xxyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyz_0[i] * pb_z + g_0_xxyyyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyzz_0[i] = 2.0 * g_0_xxyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzz_0[i] * pb_z + g_0_xxyyyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xzzz_0[i] = 3.0 * g_0_xxyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xzzz_0[i] * pb_z + g_0_xxyyyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyy_0[i] = g_0_xxyyyy_0_yyyy_0[i] * pb_z + g_0_xxyyyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyz_0[i] = g_0_xxyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyz_0[i] * pb_z + g_0_xxyyyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyzz_0[i] = 2.0 * g_0_xxyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyzz_0[i] * pb_z + g_0_xxyyyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yzzz_0[i] = 3.0 * g_0_xxyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yzzz_0[i] * pb_z + g_0_xxyyyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_zzzz_0[i] = 4.0 * g_0_xxyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_zzzz_0[i] * pb_z + g_0_xxyyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 255-270 components of targeted buffer : SKSG

    auto g_0_xxyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 255);

    auto g_0_xxyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 256);

    auto g_0_xxyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 257);

    auto g_0_xxyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 258);

    auto g_0_xxyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 259);

    auto g_0_xxyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 260);

    auto g_0_xxyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 261);

    auto g_0_xxyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 262);

    auto g_0_xxyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 263);

    auto g_0_xxyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 264);

    auto g_0_xxyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 265);

    auto g_0_xxyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 266);

    auto g_0_xxyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 267);

    auto g_0_xxyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 268);

    auto g_0_xxyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 269);

#pragma omp simd aligned(g_0_xxyyy_0_xxxy_0,       \
                             g_0_xxyyy_0_xxxy_1,   \
                             g_0_xxyyy_0_xxyy_0,   \
                             g_0_xxyyy_0_xxyy_1,   \
                             g_0_xxyyy_0_xyyy_0,   \
                             g_0_xxyyy_0_xyyy_1,   \
                             g_0_xxyyyz_0_xxxy_0,  \
                             g_0_xxyyyz_0_xxxy_1,  \
                             g_0_xxyyyz_0_xxyy_0,  \
                             g_0_xxyyyz_0_xxyy_1,  \
                             g_0_xxyyyz_0_xyyy_0,  \
                             g_0_xxyyyz_0_xyyy_1,  \
                             g_0_xxyyyzz_0_xxxx_0, \
                             g_0_xxyyyzz_0_xxxy_0, \
                             g_0_xxyyyzz_0_xxxz_0, \
                             g_0_xxyyyzz_0_xxyy_0, \
                             g_0_xxyyyzz_0_xxyz_0, \
                             g_0_xxyyyzz_0_xxzz_0, \
                             g_0_xxyyyzz_0_xyyy_0, \
                             g_0_xxyyyzz_0_xyyz_0, \
                             g_0_xxyyyzz_0_xyzz_0, \
                             g_0_xxyyyzz_0_xzzz_0, \
                             g_0_xxyyyzz_0_yyyy_0, \
                             g_0_xxyyyzz_0_yyyz_0, \
                             g_0_xxyyyzz_0_yyzz_0, \
                             g_0_xxyyyzz_0_yzzz_0, \
                             g_0_xxyyyzz_0_zzzz_0, \
                             g_0_xxyyzz_0_xxxx_0,  \
                             g_0_xxyyzz_0_xxxx_1,  \
                             g_0_xxyyzz_0_xxxz_0,  \
                             g_0_xxyyzz_0_xxxz_1,  \
                             g_0_xxyyzz_0_xxzz_0,  \
                             g_0_xxyyzz_0_xxzz_1,  \
                             g_0_xxyyzz_0_xzzz_0,  \
                             g_0_xxyyzz_0_xzzz_1,  \
                             g_0_xxyzz_0_xxxx_0,   \
                             g_0_xxyzz_0_xxxx_1,   \
                             g_0_xxyzz_0_xxxz_0,   \
                             g_0_xxyzz_0_xxxz_1,   \
                             g_0_xxyzz_0_xxzz_0,   \
                             g_0_xxyzz_0_xxzz_1,   \
                             g_0_xxyzz_0_xzzz_0,   \
                             g_0_xxyzz_0_xzzz_1,   \
                             g_0_xyyyzz_0_xxyz_0,  \
                             g_0_xyyyzz_0_xxyz_1,  \
                             g_0_xyyyzz_0_xyyz_0,  \
                             g_0_xyyyzz_0_xyyz_1,  \
                             g_0_xyyyzz_0_xyz_1,   \
                             g_0_xyyyzz_0_xyzz_0,  \
                             g_0_xyyyzz_0_xyzz_1,  \
                             g_0_xyyyzz_0_yyyy_0,  \
                             g_0_xyyyzz_0_yyyy_1,  \
                             g_0_xyyyzz_0_yyyz_0,  \
                             g_0_xyyyzz_0_yyyz_1,  \
                             g_0_xyyyzz_0_yyz_1,   \
                             g_0_xyyyzz_0_yyzz_0,  \
                             g_0_xyyyzz_0_yyzz_1,  \
                             g_0_xyyyzz_0_yzz_1,   \
                             g_0_xyyyzz_0_yzzz_0,  \
                             g_0_xyyyzz_0_yzzz_1,  \
                             g_0_xyyyzz_0_zzzz_0,  \
                             g_0_xyyyzz_0_zzzz_1,  \
                             g_0_yyyzz_0_xxyz_0,   \
                             g_0_yyyzz_0_xxyz_1,   \
                             g_0_yyyzz_0_xyyz_0,   \
                             g_0_yyyzz_0_xyyz_1,   \
                             g_0_yyyzz_0_xyzz_0,   \
                             g_0_yyyzz_0_xyzz_1,   \
                             g_0_yyyzz_0_yyyy_0,   \
                             g_0_yyyzz_0_yyyy_1,   \
                             g_0_yyyzz_0_yyyz_0,   \
                             g_0_yyyzz_0_yyyz_1,   \
                             g_0_yyyzz_0_yyzz_0,   \
                             g_0_yyyzz_0_yyzz_1,   \
                             g_0_yyyzz_0_yzzz_0,   \
                             g_0_yyyzz_0_yzzz_1,   \
                             g_0_yyyzz_0_zzzz_0,   \
                             g_0_yyyzz_0_zzzz_1,   \
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

        g_0_xxyyyzz_0_xxxx_0[i] = 2.0 * g_0_xxyzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxx_0[i] * pb_y +
                                  g_0_xxyyzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxy_0[i] =
            g_0_xxyyy_0_xxxy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxy_0[i] * pb_z + g_0_xxyyyz_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxz_0[i] = 2.0 * g_0_xxyzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxz_0[i] * pb_y +
                                  g_0_xxyyzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxyy_0[i] =
            g_0_xxyyy_0_xxyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxyy_0[i] * pb_z + g_0_xxyyyz_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxyz_0[i] = g_0_yyyzz_0_xxyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyzz_0_xxyz_0[i] * pb_x + g_0_xyyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxzz_0[i] = 2.0 * g_0_xxyzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxzz_0[i] * pb_y +
                                  g_0_xxyyzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xyyy_0[i] =
            g_0_xxyyy_0_xyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xyyy_0[i] * pb_z + g_0_xxyyyz_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xyyz_0[i] = g_0_yyyzz_0_xyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xyyyzz_0_xyyz_0[i] * pb_x + g_0_xyyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyzz_0[i] = g_0_yyyzz_0_xyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xyyyzz_0_xyzz_0[i] * pb_x + g_0_xyyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xzzz_0[i] = 2.0 * g_0_xxyzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xzzz_0[i] * pb_y +
                                  g_0_xxyyzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_yyyy_0[i] =
            g_0_yyyzz_0_yyyy_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyy_0[i] * pb_x + g_0_xyyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyz_0[i] =
            g_0_yyyzz_0_yyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyz_0[i] * pb_x + g_0_xyyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyzz_0[i] =
            g_0_yyyzz_0_yyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzz_0[i] * pb_x + g_0_xyyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yzzz_0[i] =
            g_0_yyyzz_0_yzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzz_0[i] * pb_x + g_0_xyyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_zzzz_0[i] =
            g_0_yyyzz_0_zzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_zzzz_0[i] * pb_x + g_0_xyyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 270-285 components of targeted buffer : SKSG

    auto g_0_xxyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 270);

    auto g_0_xxyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 271);

    auto g_0_xxyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 272);

    auto g_0_xxyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 273);

    auto g_0_xxyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 274);

    auto g_0_xxyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 275);

    auto g_0_xxyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 276);

    auto g_0_xxyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 277);

    auto g_0_xxyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 278);

    auto g_0_xxyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 279);

    auto g_0_xxyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 280);

    auto g_0_xxyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 281);

    auto g_0_xxyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 282);

    auto g_0_xxyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 283);

    auto g_0_xxyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 284);

#pragma omp simd aligned(g_0_xxyyz_0_xxxy_0,       \
                             g_0_xxyyz_0_xxxy_1,   \
                             g_0_xxyyz_0_xxyy_0,   \
                             g_0_xxyyz_0_xxyy_1,   \
                             g_0_xxyyz_0_xyyy_0,   \
                             g_0_xxyyz_0_xyyy_1,   \
                             g_0_xxyyzz_0_xxxy_0,  \
                             g_0_xxyyzz_0_xxxy_1,  \
                             g_0_xxyyzz_0_xxyy_0,  \
                             g_0_xxyyzz_0_xxyy_1,  \
                             g_0_xxyyzz_0_xyyy_0,  \
                             g_0_xxyyzz_0_xyyy_1,  \
                             g_0_xxyyzzz_0_xxxx_0, \
                             g_0_xxyyzzz_0_xxxy_0, \
                             g_0_xxyyzzz_0_xxxz_0, \
                             g_0_xxyyzzz_0_xxyy_0, \
                             g_0_xxyyzzz_0_xxyz_0, \
                             g_0_xxyyzzz_0_xxzz_0, \
                             g_0_xxyyzzz_0_xyyy_0, \
                             g_0_xxyyzzz_0_xyyz_0, \
                             g_0_xxyyzzz_0_xyzz_0, \
                             g_0_xxyyzzz_0_xzzz_0, \
                             g_0_xxyyzzz_0_yyyy_0, \
                             g_0_xxyyzzz_0_yyyz_0, \
                             g_0_xxyyzzz_0_yyzz_0, \
                             g_0_xxyyzzz_0_yzzz_0, \
                             g_0_xxyyzzz_0_zzzz_0, \
                             g_0_xxyzzz_0_xxxx_0,  \
                             g_0_xxyzzz_0_xxxx_1,  \
                             g_0_xxyzzz_0_xxxz_0,  \
                             g_0_xxyzzz_0_xxxz_1,  \
                             g_0_xxyzzz_0_xxzz_0,  \
                             g_0_xxyzzz_0_xxzz_1,  \
                             g_0_xxyzzz_0_xzzz_0,  \
                             g_0_xxyzzz_0_xzzz_1,  \
                             g_0_xxzzz_0_xxxx_0,   \
                             g_0_xxzzz_0_xxxx_1,   \
                             g_0_xxzzz_0_xxxz_0,   \
                             g_0_xxzzz_0_xxxz_1,   \
                             g_0_xxzzz_0_xxzz_0,   \
                             g_0_xxzzz_0_xxzz_1,   \
                             g_0_xxzzz_0_xzzz_0,   \
                             g_0_xxzzz_0_xzzz_1,   \
                             g_0_xyyzzz_0_xxyz_0,  \
                             g_0_xyyzzz_0_xxyz_1,  \
                             g_0_xyyzzz_0_xyyz_0,  \
                             g_0_xyyzzz_0_xyyz_1,  \
                             g_0_xyyzzz_0_xyz_1,   \
                             g_0_xyyzzz_0_xyzz_0,  \
                             g_0_xyyzzz_0_xyzz_1,  \
                             g_0_xyyzzz_0_yyyy_0,  \
                             g_0_xyyzzz_0_yyyy_1,  \
                             g_0_xyyzzz_0_yyyz_0,  \
                             g_0_xyyzzz_0_yyyz_1,  \
                             g_0_xyyzzz_0_yyz_1,   \
                             g_0_xyyzzz_0_yyzz_0,  \
                             g_0_xyyzzz_0_yyzz_1,  \
                             g_0_xyyzzz_0_yzz_1,   \
                             g_0_xyyzzz_0_yzzz_0,  \
                             g_0_xyyzzz_0_yzzz_1,  \
                             g_0_xyyzzz_0_zzzz_0,  \
                             g_0_xyyzzz_0_zzzz_1,  \
                             g_0_yyzzz_0_xxyz_0,   \
                             g_0_yyzzz_0_xxyz_1,   \
                             g_0_yyzzz_0_xyyz_0,   \
                             g_0_yyzzz_0_xyyz_1,   \
                             g_0_yyzzz_0_xyzz_0,   \
                             g_0_yyzzz_0_xyzz_1,   \
                             g_0_yyzzz_0_yyyy_0,   \
                             g_0_yyzzz_0_yyyy_1,   \
                             g_0_yyzzz_0_yyyz_0,   \
                             g_0_yyzzz_0_yyyz_1,   \
                             g_0_yyzzz_0_yyzz_0,   \
                             g_0_yyzzz_0_yyzz_1,   \
                             g_0_yyzzz_0_yzzz_0,   \
                             g_0_yyzzz_0_yzzz_1,   \
                             g_0_yyzzz_0_zzzz_0,   \
                             g_0_yyzzz_0_zzzz_1,   \
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

        g_0_xxyyzzz_0_xxxx_0[i] =
            g_0_xxzzz_0_xxxx_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxx_0[i] * pb_y + g_0_xxyzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxy_0[i] = 2.0 * g_0_xxyyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxy_0[i] * pb_z +
                                  g_0_xxyyzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxz_0[i] =
            g_0_xxzzz_0_xxxz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxz_0[i] * pb_y + g_0_xxyzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxyy_0[i] = 2.0 * g_0_xxyyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxyy_0[i] * pb_z +
                                  g_0_xxyyzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxyz_0[i] = g_0_yyzzz_0_xxyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyz_1[i] * fi_abcd_0 +
                                  g_0_xyyzzz_0_xxyz_0[i] * pb_x + g_0_xyyzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxzz_0[i] =
            g_0_xxzzz_0_xxzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxzz_0[i] * pb_y + g_0_xxyzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xyyy_0[i] = 2.0 * g_0_xxyyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xyyy_0[i] * pb_z +
                                  g_0_xxyyzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xyyz_0[i] = g_0_yyzzz_0_xyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xyyzzz_0_xyyz_0[i] * pb_x + g_0_xyyzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyzz_0[i] = g_0_yyzzz_0_xyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xyyzzz_0_xyzz_0[i] * pb_x + g_0_xyyzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xzzz_0[i] =
            g_0_xxzzz_0_xzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xzzz_0[i] * pb_y + g_0_xxyzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_yyyy_0[i] =
            g_0_yyzzz_0_yyyy_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyy_0[i] * pb_x + g_0_xyyzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyz_0[i] =
            g_0_yyzzz_0_yyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyz_0[i] * pb_x + g_0_xyyzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyzz_0[i] =
            g_0_yyzzz_0_yyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzz_0[i] * pb_x + g_0_xyyzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yzzz_0[i] =
            g_0_yyzzz_0_yzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzz_0[i] * pb_x + g_0_xyyzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_zzzz_0[i] =
            g_0_yyzzz_0_zzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_zzzz_0[i] * pb_x + g_0_xyyzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 285-300 components of targeted buffer : SKSG

    auto g_0_xxyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 285);

    auto g_0_xxyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 286);

    auto g_0_xxyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 287);

    auto g_0_xxyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 288);

    auto g_0_xxyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 289);

    auto g_0_xxyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 290);

    auto g_0_xxyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 291);

    auto g_0_xxyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 292);

    auto g_0_xxyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 293);

    auto g_0_xxyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 294);

    auto g_0_xxyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 295);

    auto g_0_xxyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 296);

    auto g_0_xxyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 297);

    auto g_0_xxyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 298);

    auto g_0_xxyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 299);

#pragma omp simd aligned(g_0_xxyzzzz_0_xxxx_0,     \
                             g_0_xxyzzzz_0_xxxy_0, \
                             g_0_xxyzzzz_0_xxxz_0, \
                             g_0_xxyzzzz_0_xxyy_0, \
                             g_0_xxyzzzz_0_xxyz_0, \
                             g_0_xxyzzzz_0_xxzz_0, \
                             g_0_xxyzzzz_0_xyyy_0, \
                             g_0_xxyzzzz_0_xyyz_0, \
                             g_0_xxyzzzz_0_xyzz_0, \
                             g_0_xxyzzzz_0_xzzz_0, \
                             g_0_xxyzzzz_0_yyyy_0, \
                             g_0_xxyzzzz_0_yyyz_0, \
                             g_0_xxyzzzz_0_yyzz_0, \
                             g_0_xxyzzzz_0_yzzz_0, \
                             g_0_xxyzzzz_0_zzzz_0, \
                             g_0_xxzzzz_0_xxx_1,   \
                             g_0_xxzzzz_0_xxxx_0,  \
                             g_0_xxzzzz_0_xxxx_1,  \
                             g_0_xxzzzz_0_xxxy_0,  \
                             g_0_xxzzzz_0_xxxy_1,  \
                             g_0_xxzzzz_0_xxxz_0,  \
                             g_0_xxzzzz_0_xxxz_1,  \
                             g_0_xxzzzz_0_xxy_1,   \
                             g_0_xxzzzz_0_xxyy_0,  \
                             g_0_xxzzzz_0_xxyy_1,  \
                             g_0_xxzzzz_0_xxyz_0,  \
                             g_0_xxzzzz_0_xxyz_1,  \
                             g_0_xxzzzz_0_xxz_1,   \
                             g_0_xxzzzz_0_xxzz_0,  \
                             g_0_xxzzzz_0_xxzz_1,  \
                             g_0_xxzzzz_0_xyy_1,   \
                             g_0_xxzzzz_0_xyyy_0,  \
                             g_0_xxzzzz_0_xyyy_1,  \
                             g_0_xxzzzz_0_xyyz_0,  \
                             g_0_xxzzzz_0_xyyz_1,  \
                             g_0_xxzzzz_0_xyz_1,   \
                             g_0_xxzzzz_0_xyzz_0,  \
                             g_0_xxzzzz_0_xyzz_1,  \
                             g_0_xxzzzz_0_xzz_1,   \
                             g_0_xxzzzz_0_xzzz_0,  \
                             g_0_xxzzzz_0_xzzz_1,  \
                             g_0_xxzzzz_0_yyy_1,   \
                             g_0_xxzzzz_0_yyyy_0,  \
                             g_0_xxzzzz_0_yyyy_1,  \
                             g_0_xxzzzz_0_yyyz_0,  \
                             g_0_xxzzzz_0_yyyz_1,  \
                             g_0_xxzzzz_0_yyz_1,   \
                             g_0_xxzzzz_0_yyzz_0,  \
                             g_0_xxzzzz_0_yyzz_1,  \
                             g_0_xxzzzz_0_yzz_1,   \
                             g_0_xxzzzz_0_yzzz_0,  \
                             g_0_xxzzzz_0_yzzz_1,  \
                             g_0_xxzzzz_0_zzz_1,   \
                             g_0_xxzzzz_0_zzzz_0,  \
                             g_0_xxzzzz_0_zzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_xxxx_0[i] = g_0_xxzzzz_0_xxxx_0[i] * pb_y + g_0_xxzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxy_0[i] = g_0_xxzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxy_0[i] * pb_y + g_0_xxzzzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxz_0[i] = g_0_xxzzzz_0_xxxz_0[i] * pb_y + g_0_xxzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyy_0[i] = 2.0 * g_0_xxzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyy_0[i] * pb_y + g_0_xxzzzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyz_0[i] = g_0_xxzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyz_0[i] * pb_y + g_0_xxzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxzz_0[i] = g_0_xxzzzz_0_xxzz_0[i] * pb_y + g_0_xxzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyy_0[i] = 3.0 * g_0_xxzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyy_0[i] * pb_y + g_0_xxzzzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyz_0[i] = 2.0 * g_0_xxzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyz_0[i] * pb_y + g_0_xxzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyzz_0[i] = g_0_xxzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzz_0[i] * pb_y + g_0_xxzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xzzz_0[i] = g_0_xxzzzz_0_xzzz_0[i] * pb_y + g_0_xxzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyy_0[i] = 4.0 * g_0_xxzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyy_0[i] * pb_y + g_0_xxzzzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyz_0[i] = 3.0 * g_0_xxzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyz_0[i] * pb_y + g_0_xxzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyzz_0[i] = 2.0 * g_0_xxzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyzz_0[i] * pb_y + g_0_xxzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yzzz_0[i] = g_0_xxzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yzzz_0[i] * pb_y + g_0_xxzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_zzzz_0[i] = g_0_xxzzzz_0_zzzz_0[i] * pb_y + g_0_xxzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 300-315 components of targeted buffer : SKSG

    auto g_0_xxzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 300);

    auto g_0_xxzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 301);

    auto g_0_xxzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 302);

    auto g_0_xxzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 303);

    auto g_0_xxzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 304);

    auto g_0_xxzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 305);

    auto g_0_xxzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 306);

    auto g_0_xxzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 307);

    auto g_0_xxzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 308);

    auto g_0_xxzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 309);

    auto g_0_xxzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 310);

    auto g_0_xxzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 311);

    auto g_0_xxzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 312);

    auto g_0_xxzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 313);

    auto g_0_xxzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 314);

#pragma omp simd aligned(g_0_xxzzz_0_xxxx_0,       \
                             g_0_xxzzz_0_xxxx_1,   \
                             g_0_xxzzz_0_xxxy_0,   \
                             g_0_xxzzz_0_xxxy_1,   \
                             g_0_xxzzz_0_xxyy_0,   \
                             g_0_xxzzz_0_xxyy_1,   \
                             g_0_xxzzz_0_xyyy_0,   \
                             g_0_xxzzz_0_xyyy_1,   \
                             g_0_xxzzzz_0_xxxx_0,  \
                             g_0_xxzzzz_0_xxxx_1,  \
                             g_0_xxzzzz_0_xxxy_0,  \
                             g_0_xxzzzz_0_xxxy_1,  \
                             g_0_xxzzzz_0_xxyy_0,  \
                             g_0_xxzzzz_0_xxyy_1,  \
                             g_0_xxzzzz_0_xyyy_0,  \
                             g_0_xxzzzz_0_xyyy_1,  \
                             g_0_xxzzzzz_0_xxxx_0, \
                             g_0_xxzzzzz_0_xxxy_0, \
                             g_0_xxzzzzz_0_xxxz_0, \
                             g_0_xxzzzzz_0_xxyy_0, \
                             g_0_xxzzzzz_0_xxyz_0, \
                             g_0_xxzzzzz_0_xxzz_0, \
                             g_0_xxzzzzz_0_xyyy_0, \
                             g_0_xxzzzzz_0_xyyz_0, \
                             g_0_xxzzzzz_0_xyzz_0, \
                             g_0_xxzzzzz_0_xzzz_0, \
                             g_0_xxzzzzz_0_yyyy_0, \
                             g_0_xxzzzzz_0_yyyz_0, \
                             g_0_xxzzzzz_0_yyzz_0, \
                             g_0_xxzzzzz_0_yzzz_0, \
                             g_0_xxzzzzz_0_zzzz_0, \
                             g_0_xzzzzz_0_xxxz_0,  \
                             g_0_xzzzzz_0_xxxz_1,  \
                             g_0_xzzzzz_0_xxyz_0,  \
                             g_0_xzzzzz_0_xxyz_1,  \
                             g_0_xzzzzz_0_xxz_1,   \
                             g_0_xzzzzz_0_xxzz_0,  \
                             g_0_xzzzzz_0_xxzz_1,  \
                             g_0_xzzzzz_0_xyyz_0,  \
                             g_0_xzzzzz_0_xyyz_1,  \
                             g_0_xzzzzz_0_xyz_1,   \
                             g_0_xzzzzz_0_xyzz_0,  \
                             g_0_xzzzzz_0_xyzz_1,  \
                             g_0_xzzzzz_0_xzz_1,   \
                             g_0_xzzzzz_0_xzzz_0,  \
                             g_0_xzzzzz_0_xzzz_1,  \
                             g_0_xzzzzz_0_yyyy_0,  \
                             g_0_xzzzzz_0_yyyy_1,  \
                             g_0_xzzzzz_0_yyyz_0,  \
                             g_0_xzzzzz_0_yyyz_1,  \
                             g_0_xzzzzz_0_yyz_1,   \
                             g_0_xzzzzz_0_yyzz_0,  \
                             g_0_xzzzzz_0_yyzz_1,  \
                             g_0_xzzzzz_0_yzz_1,   \
                             g_0_xzzzzz_0_yzzz_0,  \
                             g_0_xzzzzz_0_yzzz_1,  \
                             g_0_xzzzzz_0_zzz_1,   \
                             g_0_xzzzzz_0_zzzz_0,  \
                             g_0_xzzzzz_0_zzzz_1,  \
                             g_0_zzzzz_0_xxxz_0,   \
                             g_0_zzzzz_0_xxxz_1,   \
                             g_0_zzzzz_0_xxyz_0,   \
                             g_0_zzzzz_0_xxyz_1,   \
                             g_0_zzzzz_0_xxzz_0,   \
                             g_0_zzzzz_0_xxzz_1,   \
                             g_0_zzzzz_0_xyyz_0,   \
                             g_0_zzzzz_0_xyyz_1,   \
                             g_0_zzzzz_0_xyzz_0,   \
                             g_0_zzzzz_0_xyzz_1,   \
                             g_0_zzzzz_0_xzzz_0,   \
                             g_0_zzzzz_0_xzzz_1,   \
                             g_0_zzzzz_0_yyyy_0,   \
                             g_0_zzzzz_0_yyyy_1,   \
                             g_0_zzzzz_0_yyyz_0,   \
                             g_0_zzzzz_0_yyyz_1,   \
                             g_0_zzzzz_0_yyzz_0,   \
                             g_0_zzzzz_0_yyzz_1,   \
                             g_0_zzzzz_0_yzzz_0,   \
                             g_0_zzzzz_0_yzzz_1,   \
                             g_0_zzzzz_0_zzzz_0,   \
                             g_0_zzzzz_0_zzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_xxxx_0[i] = 4.0 * g_0_xxzzz_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxx_0[i] * pb_z +
                                  g_0_xxzzzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxy_0[i] = 4.0 * g_0_xxzzz_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxy_0[i] * pb_z +
                                  g_0_xxzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxz_0[i] = g_0_zzzzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxz_1[i] * fi_abcd_0 +
                                  g_0_xzzzzz_0_xxxz_0[i] * pb_x + g_0_xzzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyy_0[i] = 4.0 * g_0_xxzzz_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxyy_0[i] * pb_z +
                                  g_0_xxzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxyz_0[i] = g_0_zzzzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyz_1[i] * fi_abcd_0 +
                                  g_0_xzzzzz_0_xxyz_0[i] * pb_x + g_0_xzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxzz_0[i] = g_0_zzzzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzzz_0_xxzz_0[i] * pb_x + g_0_xzzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyy_0[i] = 4.0 * g_0_xxzzz_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xyyy_0[i] * pb_z +
                                  g_0_xxzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xyyz_0[i] = g_0_zzzzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_xzzzzz_0_xyyz_0[i] * pb_x + g_0_xzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyzz_0[i] = g_0_zzzzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzzz_0_xyzz_0[i] * pb_x + g_0_xzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xzzz_0[i] = g_0_zzzzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_xzzzzz_0_xzzz_0[i] * pb_x + g_0_xzzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyy_0[i] =
            g_0_zzzzz_0_yyyy_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyy_0[i] * pb_x + g_0_xzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyz_0[i] =
            g_0_zzzzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyz_0[i] * pb_x + g_0_xzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyzz_0[i] =
            g_0_zzzzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzz_0[i] * pb_x + g_0_xzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yzzz_0[i] =
            g_0_zzzzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzz_0[i] * pb_x + g_0_xzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_zzzz_0[i] =
            g_0_zzzzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzz_0[i] * pb_x + g_0_xzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 315-330 components of targeted buffer : SKSG

    auto g_0_xyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 315);

    auto g_0_xyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 316);

    auto g_0_xyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 317);

    auto g_0_xyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 318);

    auto g_0_xyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 319);

    auto g_0_xyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 320);

    auto g_0_xyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 321);

    auto g_0_xyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 322);

    auto g_0_xyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 323);

    auto g_0_xyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 324);

    auto g_0_xyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 325);

    auto g_0_xyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 326);

    auto g_0_xyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 327);

    auto g_0_xyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 328);

    auto g_0_xyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 329);

#pragma omp simd aligned(g_0_xyyyyyy_0_xxxx_0,     \
                             g_0_xyyyyyy_0_xxxy_0, \
                             g_0_xyyyyyy_0_xxxz_0, \
                             g_0_xyyyyyy_0_xxyy_0, \
                             g_0_xyyyyyy_0_xxyz_0, \
                             g_0_xyyyyyy_0_xxzz_0, \
                             g_0_xyyyyyy_0_xyyy_0, \
                             g_0_xyyyyyy_0_xyyz_0, \
                             g_0_xyyyyyy_0_xyzz_0, \
                             g_0_xyyyyyy_0_xzzz_0, \
                             g_0_xyyyyyy_0_yyyy_0, \
                             g_0_xyyyyyy_0_yyyz_0, \
                             g_0_xyyyyyy_0_yyzz_0, \
                             g_0_xyyyyyy_0_yzzz_0, \
                             g_0_xyyyyyy_0_zzzz_0, \
                             g_0_yyyyyy_0_xxx_1,   \
                             g_0_yyyyyy_0_xxxx_0,  \
                             g_0_yyyyyy_0_xxxx_1,  \
                             g_0_yyyyyy_0_xxxy_0,  \
                             g_0_yyyyyy_0_xxxy_1,  \
                             g_0_yyyyyy_0_xxxz_0,  \
                             g_0_yyyyyy_0_xxxz_1,  \
                             g_0_yyyyyy_0_xxy_1,   \
                             g_0_yyyyyy_0_xxyy_0,  \
                             g_0_yyyyyy_0_xxyy_1,  \
                             g_0_yyyyyy_0_xxyz_0,  \
                             g_0_yyyyyy_0_xxyz_1,  \
                             g_0_yyyyyy_0_xxz_1,   \
                             g_0_yyyyyy_0_xxzz_0,  \
                             g_0_yyyyyy_0_xxzz_1,  \
                             g_0_yyyyyy_0_xyy_1,   \
                             g_0_yyyyyy_0_xyyy_0,  \
                             g_0_yyyyyy_0_xyyy_1,  \
                             g_0_yyyyyy_0_xyyz_0,  \
                             g_0_yyyyyy_0_xyyz_1,  \
                             g_0_yyyyyy_0_xyz_1,   \
                             g_0_yyyyyy_0_xyzz_0,  \
                             g_0_yyyyyy_0_xyzz_1,  \
                             g_0_yyyyyy_0_xzz_1,   \
                             g_0_yyyyyy_0_xzzz_0,  \
                             g_0_yyyyyy_0_xzzz_1,  \
                             g_0_yyyyyy_0_yyy_1,   \
                             g_0_yyyyyy_0_yyyy_0,  \
                             g_0_yyyyyy_0_yyyy_1,  \
                             g_0_yyyyyy_0_yyyz_0,  \
                             g_0_yyyyyy_0_yyyz_1,  \
                             g_0_yyyyyy_0_yyz_1,   \
                             g_0_yyyyyy_0_yyzz_0,  \
                             g_0_yyyyyy_0_yyzz_1,  \
                             g_0_yyyyyy_0_yzz_1,   \
                             g_0_yyyyyy_0_yzzz_0,  \
                             g_0_yyyyyy_0_yzzz_1,  \
                             g_0_yyyyyy_0_zzz_1,   \
                             g_0_yyyyyy_0_zzzz_0,  \
                             g_0_yyyyyy_0_zzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_xxxx_0[i] = 4.0 * g_0_yyyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxx_0[i] * pb_x + g_0_yyyyyy_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxy_0[i] = 3.0 * g_0_yyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxy_0[i] * pb_x + g_0_yyyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxz_0[i] = 3.0 * g_0_yyyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxz_0[i] * pb_x + g_0_yyyyyy_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyy_0[i] = 2.0 * g_0_yyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyy_0[i] * pb_x + g_0_yyyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyz_0[i] = 2.0 * g_0_yyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyz_0[i] * pb_x + g_0_yyyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxzz_0[i] = 2.0 * g_0_yyyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzz_0[i] * pb_x + g_0_yyyyyy_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyy_0[i] = g_0_yyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyy_0[i] * pb_x + g_0_yyyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyz_0[i] = g_0_yyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyz_0[i] * pb_x + g_0_yyyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyzz_0[i] = g_0_yyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzz_0[i] * pb_x + g_0_yyyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xzzz_0[i] = g_0_yyyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzz_0[i] * pb_x + g_0_yyyyyy_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyy_0[i] = g_0_yyyyyy_0_yyyy_0[i] * pb_x + g_0_yyyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyz_0[i] = g_0_yyyyyy_0_yyyz_0[i] * pb_x + g_0_yyyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyzz_0[i] = g_0_yyyyyy_0_yyzz_0[i] * pb_x + g_0_yyyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yzzz_0[i] = g_0_yyyyyy_0_yzzz_0[i] * pb_x + g_0_yyyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_zzzz_0[i] = g_0_yyyyyy_0_zzzz_0[i] * pb_x + g_0_yyyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 330-345 components of targeted buffer : SKSG

    auto g_0_xyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 330);

    auto g_0_xyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 331);

    auto g_0_xyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 332);

    auto g_0_xyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 333);

    auto g_0_xyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 334);

    auto g_0_xyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 335);

    auto g_0_xyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 336);

    auto g_0_xyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 337);

    auto g_0_xyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 338);

    auto g_0_xyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 339);

    auto g_0_xyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 340);

    auto g_0_xyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 341);

    auto g_0_xyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 342);

    auto g_0_xyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 343);

    auto g_0_xyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 344);

#pragma omp simd aligned(g_0_xyyyyy_0_xxxx_0,      \
                             g_0_xyyyyy_0_xxxx_1,  \
                             g_0_xyyyyy_0_xxxy_0,  \
                             g_0_xyyyyy_0_xxxy_1,  \
                             g_0_xyyyyy_0_xxyy_0,  \
                             g_0_xyyyyy_0_xxyy_1,  \
                             g_0_xyyyyy_0_xyyy_0,  \
                             g_0_xyyyyy_0_xyyy_1,  \
                             g_0_xyyyyyz_0_xxxx_0, \
                             g_0_xyyyyyz_0_xxxy_0, \
                             g_0_xyyyyyz_0_xxxz_0, \
                             g_0_xyyyyyz_0_xxyy_0, \
                             g_0_xyyyyyz_0_xxyz_0, \
                             g_0_xyyyyyz_0_xxzz_0, \
                             g_0_xyyyyyz_0_xyyy_0, \
                             g_0_xyyyyyz_0_xyyz_0, \
                             g_0_xyyyyyz_0_xyzz_0, \
                             g_0_xyyyyyz_0_xzzz_0, \
                             g_0_xyyyyyz_0_yyyy_0, \
                             g_0_xyyyyyz_0_yyyz_0, \
                             g_0_xyyyyyz_0_yyzz_0, \
                             g_0_xyyyyyz_0_yzzz_0, \
                             g_0_xyyyyyz_0_zzzz_0, \
                             g_0_yyyyyz_0_xxxz_0,  \
                             g_0_yyyyyz_0_xxxz_1,  \
                             g_0_yyyyyz_0_xxyz_0,  \
                             g_0_yyyyyz_0_xxyz_1,  \
                             g_0_yyyyyz_0_xxz_1,   \
                             g_0_yyyyyz_0_xxzz_0,  \
                             g_0_yyyyyz_0_xxzz_1,  \
                             g_0_yyyyyz_0_xyyz_0,  \
                             g_0_yyyyyz_0_xyyz_1,  \
                             g_0_yyyyyz_0_xyz_1,   \
                             g_0_yyyyyz_0_xyzz_0,  \
                             g_0_yyyyyz_0_xyzz_1,  \
                             g_0_yyyyyz_0_xzz_1,   \
                             g_0_yyyyyz_0_xzzz_0,  \
                             g_0_yyyyyz_0_xzzz_1,  \
                             g_0_yyyyyz_0_yyyy_0,  \
                             g_0_yyyyyz_0_yyyy_1,  \
                             g_0_yyyyyz_0_yyyz_0,  \
                             g_0_yyyyyz_0_yyyz_1,  \
                             g_0_yyyyyz_0_yyz_1,   \
                             g_0_yyyyyz_0_yyzz_0,  \
                             g_0_yyyyyz_0_yyzz_1,  \
                             g_0_yyyyyz_0_yzz_1,   \
                             g_0_yyyyyz_0_yzzz_0,  \
                             g_0_yyyyyz_0_yzzz_1,  \
                             g_0_yyyyyz_0_zzz_1,   \
                             g_0_yyyyyz_0_zzzz_0,  \
                             g_0_yyyyyz_0_zzzz_1,  \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyz_0_xxxx_0[i] = g_0_xyyyyy_0_xxxx_0[i] * pb_z + g_0_xyyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxy_0[i] = g_0_xyyyyy_0_xxxy_0[i] * pb_z + g_0_xyyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxz_0[i] = 3.0 * g_0_yyyyyz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxz_0[i] * pb_x + g_0_yyyyyz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyy_0[i] = g_0_xyyyyy_0_xxyy_0[i] * pb_z + g_0_xyyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxyz_0[i] = 2.0 * g_0_yyyyyz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyz_0[i] * pb_x + g_0_yyyyyz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyyyz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxzz_0[i] * pb_x + g_0_yyyyyz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyy_0[i] = g_0_xyyyyy_0_xyyy_0[i] * pb_z + g_0_xyyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xyyz_0[i] = g_0_yyyyyz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyz_0[i] * pb_x + g_0_yyyyyz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyzz_0[i] = g_0_yyyyyz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyzz_0[i] * pb_x + g_0_yyyyyz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xzzz_0[i] = g_0_yyyyyz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xzzz_0[i] * pb_x + g_0_yyyyyz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyy_0[i] = g_0_yyyyyz_0_yyyy_0[i] * pb_x + g_0_yyyyyz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyz_0[i] = g_0_yyyyyz_0_yyyz_0[i] * pb_x + g_0_yyyyyz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyzz_0[i] = g_0_yyyyyz_0_yyzz_0[i] * pb_x + g_0_yyyyyz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yzzz_0[i] = g_0_yyyyyz_0_yzzz_0[i] * pb_x + g_0_yyyyyz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_zzzz_0[i] = g_0_yyyyyz_0_zzzz_0[i] * pb_x + g_0_yyyyyz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 345-360 components of targeted buffer : SKSG

    auto g_0_xyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 345);

    auto g_0_xyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 346);

    auto g_0_xyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 347);

    auto g_0_xyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 348);

    auto g_0_xyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 349);

    auto g_0_xyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 350);

    auto g_0_xyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 351);

    auto g_0_xyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 352);

    auto g_0_xyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 353);

    auto g_0_xyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 354);

    auto g_0_xyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 355);

    auto g_0_xyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 356);

    auto g_0_xyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 357);

    auto g_0_xyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 358);

    auto g_0_xyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 359);

#pragma omp simd aligned(g_0_xyyyyzz_0_xxxx_0,     \
                             g_0_xyyyyzz_0_xxxy_0, \
                             g_0_xyyyyzz_0_xxxz_0, \
                             g_0_xyyyyzz_0_xxyy_0, \
                             g_0_xyyyyzz_0_xxyz_0, \
                             g_0_xyyyyzz_0_xxzz_0, \
                             g_0_xyyyyzz_0_xyyy_0, \
                             g_0_xyyyyzz_0_xyyz_0, \
                             g_0_xyyyyzz_0_xyzz_0, \
                             g_0_xyyyyzz_0_xzzz_0, \
                             g_0_xyyyyzz_0_yyyy_0, \
                             g_0_xyyyyzz_0_yyyz_0, \
                             g_0_xyyyyzz_0_yyzz_0, \
                             g_0_xyyyyzz_0_yzzz_0, \
                             g_0_xyyyyzz_0_zzzz_0, \
                             g_0_yyyyzz_0_xxx_1,   \
                             g_0_yyyyzz_0_xxxx_0,  \
                             g_0_yyyyzz_0_xxxx_1,  \
                             g_0_yyyyzz_0_xxxy_0,  \
                             g_0_yyyyzz_0_xxxy_1,  \
                             g_0_yyyyzz_0_xxxz_0,  \
                             g_0_yyyyzz_0_xxxz_1,  \
                             g_0_yyyyzz_0_xxy_1,   \
                             g_0_yyyyzz_0_xxyy_0,  \
                             g_0_yyyyzz_0_xxyy_1,  \
                             g_0_yyyyzz_0_xxyz_0,  \
                             g_0_yyyyzz_0_xxyz_1,  \
                             g_0_yyyyzz_0_xxz_1,   \
                             g_0_yyyyzz_0_xxzz_0,  \
                             g_0_yyyyzz_0_xxzz_1,  \
                             g_0_yyyyzz_0_xyy_1,   \
                             g_0_yyyyzz_0_xyyy_0,  \
                             g_0_yyyyzz_0_xyyy_1,  \
                             g_0_yyyyzz_0_xyyz_0,  \
                             g_0_yyyyzz_0_xyyz_1,  \
                             g_0_yyyyzz_0_xyz_1,   \
                             g_0_yyyyzz_0_xyzz_0,  \
                             g_0_yyyyzz_0_xyzz_1,  \
                             g_0_yyyyzz_0_xzz_1,   \
                             g_0_yyyyzz_0_xzzz_0,  \
                             g_0_yyyyzz_0_xzzz_1,  \
                             g_0_yyyyzz_0_yyy_1,   \
                             g_0_yyyyzz_0_yyyy_0,  \
                             g_0_yyyyzz_0_yyyy_1,  \
                             g_0_yyyyzz_0_yyyz_0,  \
                             g_0_yyyyzz_0_yyyz_1,  \
                             g_0_yyyyzz_0_yyz_1,   \
                             g_0_yyyyzz_0_yyzz_0,  \
                             g_0_yyyyzz_0_yyzz_1,  \
                             g_0_yyyyzz_0_yzz_1,   \
                             g_0_yyyyzz_0_yzzz_0,  \
                             g_0_yyyyzz_0_yzzz_1,  \
                             g_0_yyyyzz_0_zzz_1,   \
                             g_0_yyyyzz_0_zzzz_0,  \
                             g_0_yyyyzz_0_zzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_xxxx_0[i] = 4.0 * g_0_yyyyzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxx_0[i] * pb_x + g_0_yyyyzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxy_0[i] = 3.0 * g_0_yyyyzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxy_0[i] * pb_x + g_0_yyyyzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxz_0[i] = 3.0 * g_0_yyyyzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxz_0[i] * pb_x + g_0_yyyyzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyy_0[i] = 2.0 * g_0_yyyyzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyy_0[i] * pb_x + g_0_yyyyzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyz_0[i] = 2.0 * g_0_yyyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyz_0[i] * pb_x + g_0_yyyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxzz_0[i] = 2.0 * g_0_yyyyzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxzz_0[i] * pb_x + g_0_yyyyzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyy_0[i] = g_0_yyyyzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyy_0[i] * pb_x + g_0_yyyyzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyz_0[i] = g_0_yyyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyz_0[i] * pb_x + g_0_yyyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyzz_0[i] = g_0_yyyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzz_0[i] * pb_x + g_0_yyyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xzzz_0[i] = g_0_yyyyzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xzzz_0[i] * pb_x + g_0_yyyyzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyy_0[i] = g_0_yyyyzz_0_yyyy_0[i] * pb_x + g_0_yyyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyz_0[i] = g_0_yyyyzz_0_yyyz_0[i] * pb_x + g_0_yyyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyzz_0[i] = g_0_yyyyzz_0_yyzz_0[i] * pb_x + g_0_yyyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yzzz_0[i] = g_0_yyyyzz_0_yzzz_0[i] * pb_x + g_0_yyyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_zzzz_0[i] = g_0_yyyyzz_0_zzzz_0[i] * pb_x + g_0_yyyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 360-375 components of targeted buffer : SKSG

    auto g_0_xyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 360);

    auto g_0_xyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 361);

    auto g_0_xyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 362);

    auto g_0_xyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 363);

    auto g_0_xyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 364);

    auto g_0_xyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 365);

    auto g_0_xyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 366);

    auto g_0_xyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 367);

    auto g_0_xyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 368);

    auto g_0_xyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 369);

    auto g_0_xyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 370);

    auto g_0_xyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 371);

    auto g_0_xyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 372);

    auto g_0_xyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 373);

    auto g_0_xyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 374);

#pragma omp simd aligned(g_0_xyyyzzz_0_xxxx_0,     \
                             g_0_xyyyzzz_0_xxxy_0, \
                             g_0_xyyyzzz_0_xxxz_0, \
                             g_0_xyyyzzz_0_xxyy_0, \
                             g_0_xyyyzzz_0_xxyz_0, \
                             g_0_xyyyzzz_0_xxzz_0, \
                             g_0_xyyyzzz_0_xyyy_0, \
                             g_0_xyyyzzz_0_xyyz_0, \
                             g_0_xyyyzzz_0_xyzz_0, \
                             g_0_xyyyzzz_0_xzzz_0, \
                             g_0_xyyyzzz_0_yyyy_0, \
                             g_0_xyyyzzz_0_yyyz_0, \
                             g_0_xyyyzzz_0_yyzz_0, \
                             g_0_xyyyzzz_0_yzzz_0, \
                             g_0_xyyyzzz_0_zzzz_0, \
                             g_0_yyyzzz_0_xxx_1,   \
                             g_0_yyyzzz_0_xxxx_0,  \
                             g_0_yyyzzz_0_xxxx_1,  \
                             g_0_yyyzzz_0_xxxy_0,  \
                             g_0_yyyzzz_0_xxxy_1,  \
                             g_0_yyyzzz_0_xxxz_0,  \
                             g_0_yyyzzz_0_xxxz_1,  \
                             g_0_yyyzzz_0_xxy_1,   \
                             g_0_yyyzzz_0_xxyy_0,  \
                             g_0_yyyzzz_0_xxyy_1,  \
                             g_0_yyyzzz_0_xxyz_0,  \
                             g_0_yyyzzz_0_xxyz_1,  \
                             g_0_yyyzzz_0_xxz_1,   \
                             g_0_yyyzzz_0_xxzz_0,  \
                             g_0_yyyzzz_0_xxzz_1,  \
                             g_0_yyyzzz_0_xyy_1,   \
                             g_0_yyyzzz_0_xyyy_0,  \
                             g_0_yyyzzz_0_xyyy_1,  \
                             g_0_yyyzzz_0_xyyz_0,  \
                             g_0_yyyzzz_0_xyyz_1,  \
                             g_0_yyyzzz_0_xyz_1,   \
                             g_0_yyyzzz_0_xyzz_0,  \
                             g_0_yyyzzz_0_xyzz_1,  \
                             g_0_yyyzzz_0_xzz_1,   \
                             g_0_yyyzzz_0_xzzz_0,  \
                             g_0_yyyzzz_0_xzzz_1,  \
                             g_0_yyyzzz_0_yyy_1,   \
                             g_0_yyyzzz_0_yyyy_0,  \
                             g_0_yyyzzz_0_yyyy_1,  \
                             g_0_yyyzzz_0_yyyz_0,  \
                             g_0_yyyzzz_0_yyyz_1,  \
                             g_0_yyyzzz_0_yyz_1,   \
                             g_0_yyyzzz_0_yyzz_0,  \
                             g_0_yyyzzz_0_yyzz_1,  \
                             g_0_yyyzzz_0_yzz_1,   \
                             g_0_yyyzzz_0_yzzz_0,  \
                             g_0_yyyzzz_0_yzzz_1,  \
                             g_0_yyyzzz_0_zzz_1,   \
                             g_0_yyyzzz_0_zzzz_0,  \
                             g_0_yyyzzz_0_zzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_xxxx_0[i] = 4.0 * g_0_yyyzzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxx_0[i] * pb_x + g_0_yyyzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxy_0[i] = 3.0 * g_0_yyyzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxy_0[i] * pb_x + g_0_yyyzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxz_0[i] = 3.0 * g_0_yyyzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxz_0[i] * pb_x + g_0_yyyzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyy_0[i] = 2.0 * g_0_yyyzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyy_0[i] * pb_x + g_0_yyyzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyz_0[i] = 2.0 * g_0_yyyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyz_0[i] * pb_x + g_0_yyyzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxzz_0[i] = 2.0 * g_0_yyyzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxzz_0[i] * pb_x + g_0_yyyzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyy_0[i] = g_0_yyyzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyy_0[i] * pb_x + g_0_yyyzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyz_0[i] = g_0_yyyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyz_0[i] * pb_x + g_0_yyyzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyzz_0[i] = g_0_yyyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzz_0[i] * pb_x + g_0_yyyzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xzzz_0[i] = g_0_yyyzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xzzz_0[i] * pb_x + g_0_yyyzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyy_0[i] = g_0_yyyzzz_0_yyyy_0[i] * pb_x + g_0_yyyzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyz_0[i] = g_0_yyyzzz_0_yyyz_0[i] * pb_x + g_0_yyyzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyzz_0[i] = g_0_yyyzzz_0_yyzz_0[i] * pb_x + g_0_yyyzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yzzz_0[i] = g_0_yyyzzz_0_yzzz_0[i] * pb_x + g_0_yyyzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_zzzz_0[i] = g_0_yyyzzz_0_zzzz_0[i] * pb_x + g_0_yyyzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 375-390 components of targeted buffer : SKSG

    auto g_0_xyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 375);

    auto g_0_xyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 376);

    auto g_0_xyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 377);

    auto g_0_xyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 378);

    auto g_0_xyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 379);

    auto g_0_xyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 380);

    auto g_0_xyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 381);

    auto g_0_xyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 382);

    auto g_0_xyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 383);

    auto g_0_xyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 384);

    auto g_0_xyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 385);

    auto g_0_xyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 386);

    auto g_0_xyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 387);

    auto g_0_xyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 388);

    auto g_0_xyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 389);

#pragma omp simd aligned(g_0_xyyzzzz_0_xxxx_0,     \
                             g_0_xyyzzzz_0_xxxy_0, \
                             g_0_xyyzzzz_0_xxxz_0, \
                             g_0_xyyzzzz_0_xxyy_0, \
                             g_0_xyyzzzz_0_xxyz_0, \
                             g_0_xyyzzzz_0_xxzz_0, \
                             g_0_xyyzzzz_0_xyyy_0, \
                             g_0_xyyzzzz_0_xyyz_0, \
                             g_0_xyyzzzz_0_xyzz_0, \
                             g_0_xyyzzzz_0_xzzz_0, \
                             g_0_xyyzzzz_0_yyyy_0, \
                             g_0_xyyzzzz_0_yyyz_0, \
                             g_0_xyyzzzz_0_yyzz_0, \
                             g_0_xyyzzzz_0_yzzz_0, \
                             g_0_xyyzzzz_0_zzzz_0, \
                             g_0_yyzzzz_0_xxx_1,   \
                             g_0_yyzzzz_0_xxxx_0,  \
                             g_0_yyzzzz_0_xxxx_1,  \
                             g_0_yyzzzz_0_xxxy_0,  \
                             g_0_yyzzzz_0_xxxy_1,  \
                             g_0_yyzzzz_0_xxxz_0,  \
                             g_0_yyzzzz_0_xxxz_1,  \
                             g_0_yyzzzz_0_xxy_1,   \
                             g_0_yyzzzz_0_xxyy_0,  \
                             g_0_yyzzzz_0_xxyy_1,  \
                             g_0_yyzzzz_0_xxyz_0,  \
                             g_0_yyzzzz_0_xxyz_1,  \
                             g_0_yyzzzz_0_xxz_1,   \
                             g_0_yyzzzz_0_xxzz_0,  \
                             g_0_yyzzzz_0_xxzz_1,  \
                             g_0_yyzzzz_0_xyy_1,   \
                             g_0_yyzzzz_0_xyyy_0,  \
                             g_0_yyzzzz_0_xyyy_1,  \
                             g_0_yyzzzz_0_xyyz_0,  \
                             g_0_yyzzzz_0_xyyz_1,  \
                             g_0_yyzzzz_0_xyz_1,   \
                             g_0_yyzzzz_0_xyzz_0,  \
                             g_0_yyzzzz_0_xyzz_1,  \
                             g_0_yyzzzz_0_xzz_1,   \
                             g_0_yyzzzz_0_xzzz_0,  \
                             g_0_yyzzzz_0_xzzz_1,  \
                             g_0_yyzzzz_0_yyy_1,   \
                             g_0_yyzzzz_0_yyyy_0,  \
                             g_0_yyzzzz_0_yyyy_1,  \
                             g_0_yyzzzz_0_yyyz_0,  \
                             g_0_yyzzzz_0_yyyz_1,  \
                             g_0_yyzzzz_0_yyz_1,   \
                             g_0_yyzzzz_0_yyzz_0,  \
                             g_0_yyzzzz_0_yyzz_1,  \
                             g_0_yyzzzz_0_yzz_1,   \
                             g_0_yyzzzz_0_yzzz_0,  \
                             g_0_yyzzzz_0_yzzz_1,  \
                             g_0_yyzzzz_0_zzz_1,   \
                             g_0_yyzzzz_0_zzzz_0,  \
                             g_0_yyzzzz_0_zzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_xxxx_0[i] = 4.0 * g_0_yyzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxx_0[i] * pb_x + g_0_yyzzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxy_0[i] = 3.0 * g_0_yyzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxy_0[i] * pb_x + g_0_yyzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxz_0[i] = 3.0 * g_0_yyzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxz_0[i] * pb_x + g_0_yyzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyy_0[i] = 2.0 * g_0_yyzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyy_0[i] * pb_x + g_0_yyzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyz_0[i] = 2.0 * g_0_yyzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyz_0[i] * pb_x + g_0_yyzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxzz_0[i] = 2.0 * g_0_yyzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxzz_0[i] * pb_x + g_0_yyzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyy_0[i] = g_0_yyzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyy_0[i] * pb_x + g_0_yyzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyz_0[i] = g_0_yyzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyz_0[i] * pb_x + g_0_yyzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyzz_0[i] = g_0_yyzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzz_0[i] * pb_x + g_0_yyzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xzzz_0[i] = g_0_yyzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xzzz_0[i] * pb_x + g_0_yyzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyy_0[i] = g_0_yyzzzz_0_yyyy_0[i] * pb_x + g_0_yyzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyz_0[i] = g_0_yyzzzz_0_yyyz_0[i] * pb_x + g_0_yyzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyzz_0[i] = g_0_yyzzzz_0_yyzz_0[i] * pb_x + g_0_yyzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yzzz_0[i] = g_0_yyzzzz_0_yzzz_0[i] * pb_x + g_0_yyzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_zzzz_0[i] = g_0_yyzzzz_0_zzzz_0[i] * pb_x + g_0_yyzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 390-405 components of targeted buffer : SKSG

    auto g_0_xyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 390);

    auto g_0_xyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 391);

    auto g_0_xyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 392);

    auto g_0_xyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 393);

    auto g_0_xyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 394);

    auto g_0_xyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 395);

    auto g_0_xyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 396);

    auto g_0_xyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 397);

    auto g_0_xyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 398);

    auto g_0_xyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 399);

    auto g_0_xyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 400);

    auto g_0_xyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 401);

    auto g_0_xyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 402);

    auto g_0_xyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 403);

    auto g_0_xyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 404);

#pragma omp simd aligned(g_0_xyzzzzz_0_xxxx_0,     \
                             g_0_xyzzzzz_0_xxxy_0, \
                             g_0_xyzzzzz_0_xxxz_0, \
                             g_0_xyzzzzz_0_xxyy_0, \
                             g_0_xyzzzzz_0_xxyz_0, \
                             g_0_xyzzzzz_0_xxzz_0, \
                             g_0_xyzzzzz_0_xyyy_0, \
                             g_0_xyzzzzz_0_xyyz_0, \
                             g_0_xyzzzzz_0_xyzz_0, \
                             g_0_xyzzzzz_0_xzzz_0, \
                             g_0_xyzzzzz_0_yyyy_0, \
                             g_0_xyzzzzz_0_yyyz_0, \
                             g_0_xyzzzzz_0_yyzz_0, \
                             g_0_xyzzzzz_0_yzzz_0, \
                             g_0_xyzzzzz_0_zzzz_0, \
                             g_0_xzzzzz_0_xxxx_0,  \
                             g_0_xzzzzz_0_xxxx_1,  \
                             g_0_xzzzzz_0_xxxz_0,  \
                             g_0_xzzzzz_0_xxxz_1,  \
                             g_0_xzzzzz_0_xxzz_0,  \
                             g_0_xzzzzz_0_xxzz_1,  \
                             g_0_xzzzzz_0_xzzz_0,  \
                             g_0_xzzzzz_0_xzzz_1,  \
                             g_0_yzzzzz_0_xxxy_0,  \
                             g_0_yzzzzz_0_xxxy_1,  \
                             g_0_yzzzzz_0_xxy_1,   \
                             g_0_yzzzzz_0_xxyy_0,  \
                             g_0_yzzzzz_0_xxyy_1,  \
                             g_0_yzzzzz_0_xxyz_0,  \
                             g_0_yzzzzz_0_xxyz_1,  \
                             g_0_yzzzzz_0_xyy_1,   \
                             g_0_yzzzzz_0_xyyy_0,  \
                             g_0_yzzzzz_0_xyyy_1,  \
                             g_0_yzzzzz_0_xyyz_0,  \
                             g_0_yzzzzz_0_xyyz_1,  \
                             g_0_yzzzzz_0_xyz_1,   \
                             g_0_yzzzzz_0_xyzz_0,  \
                             g_0_yzzzzz_0_xyzz_1,  \
                             g_0_yzzzzz_0_yyy_1,   \
                             g_0_yzzzzz_0_yyyy_0,  \
                             g_0_yzzzzz_0_yyyy_1,  \
                             g_0_yzzzzz_0_yyyz_0,  \
                             g_0_yzzzzz_0_yyyz_1,  \
                             g_0_yzzzzz_0_yyz_1,   \
                             g_0_yzzzzz_0_yyzz_0,  \
                             g_0_yzzzzz_0_yyzz_1,  \
                             g_0_yzzzzz_0_yzz_1,   \
                             g_0_yzzzzz_0_yzzz_0,  \
                             g_0_yzzzzz_0_yzzz_1,  \
                             g_0_yzzzzz_0_zzzz_0,  \
                             g_0_yzzzzz_0_zzzz_1,  \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzz_0_xxxx_0[i] = g_0_xzzzzz_0_xxxx_0[i] * pb_y + g_0_xzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxy_0[i] = 3.0 * g_0_yzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxy_0[i] * pb_x + g_0_yzzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxz_0[i] = g_0_xzzzzz_0_xxxz_0[i] * pb_y + g_0_xzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxyy_0[i] = 2.0 * g_0_yzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyy_0[i] * pb_x + g_0_yzzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyz_0[i] = 2.0 * g_0_yzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyz_0[i] * pb_x + g_0_yzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxzz_0[i] = g_0_xzzzzz_0_xxzz_0[i] * pb_y + g_0_xzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xyyy_0[i] = g_0_yzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyy_0[i] * pb_x + g_0_yzzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyz_0[i] = g_0_yzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyz_0[i] * pb_x + g_0_yzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyzz_0[i] = g_0_yzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyzz_0[i] * pb_x + g_0_yzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xzzz_0[i] = g_0_xzzzzz_0_xzzz_0[i] * pb_y + g_0_xzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_yyyy_0[i] = g_0_yzzzzz_0_yyyy_0[i] * pb_x + g_0_yzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyz_0[i] = g_0_yzzzzz_0_yyyz_0[i] * pb_x + g_0_yzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyzz_0[i] = g_0_yzzzzz_0_yyzz_0[i] * pb_x + g_0_yzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yzzz_0[i] = g_0_yzzzzz_0_yzzz_0[i] * pb_x + g_0_yzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_zzzz_0[i] = g_0_yzzzzz_0_zzzz_0[i] * pb_x + g_0_yzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 405-420 components of targeted buffer : SKSG

    auto g_0_xzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 405);

    auto g_0_xzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 406);

    auto g_0_xzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 407);

    auto g_0_xzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 408);

    auto g_0_xzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 409);

    auto g_0_xzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 410);

    auto g_0_xzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 411);

    auto g_0_xzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 412);

    auto g_0_xzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 413);

    auto g_0_xzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 414);

    auto g_0_xzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 415);

    auto g_0_xzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 416);

    auto g_0_xzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 417);

    auto g_0_xzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 418);

    auto g_0_xzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 419);

#pragma omp simd aligned(g_0_xzzzzzz_0_xxxx_0,     \
                             g_0_xzzzzzz_0_xxxy_0, \
                             g_0_xzzzzzz_0_xxxz_0, \
                             g_0_xzzzzzz_0_xxyy_0, \
                             g_0_xzzzzzz_0_xxyz_0, \
                             g_0_xzzzzzz_0_xxzz_0, \
                             g_0_xzzzzzz_0_xyyy_0, \
                             g_0_xzzzzzz_0_xyyz_0, \
                             g_0_xzzzzzz_0_xyzz_0, \
                             g_0_xzzzzzz_0_xzzz_0, \
                             g_0_xzzzzzz_0_yyyy_0, \
                             g_0_xzzzzzz_0_yyyz_0, \
                             g_0_xzzzzzz_0_yyzz_0, \
                             g_0_xzzzzzz_0_yzzz_0, \
                             g_0_xzzzzzz_0_zzzz_0, \
                             g_0_zzzzzz_0_xxx_1,   \
                             g_0_zzzzzz_0_xxxx_0,  \
                             g_0_zzzzzz_0_xxxx_1,  \
                             g_0_zzzzzz_0_xxxy_0,  \
                             g_0_zzzzzz_0_xxxy_1,  \
                             g_0_zzzzzz_0_xxxz_0,  \
                             g_0_zzzzzz_0_xxxz_1,  \
                             g_0_zzzzzz_0_xxy_1,   \
                             g_0_zzzzzz_0_xxyy_0,  \
                             g_0_zzzzzz_0_xxyy_1,  \
                             g_0_zzzzzz_0_xxyz_0,  \
                             g_0_zzzzzz_0_xxyz_1,  \
                             g_0_zzzzzz_0_xxz_1,   \
                             g_0_zzzzzz_0_xxzz_0,  \
                             g_0_zzzzzz_0_xxzz_1,  \
                             g_0_zzzzzz_0_xyy_1,   \
                             g_0_zzzzzz_0_xyyy_0,  \
                             g_0_zzzzzz_0_xyyy_1,  \
                             g_0_zzzzzz_0_xyyz_0,  \
                             g_0_zzzzzz_0_xyyz_1,  \
                             g_0_zzzzzz_0_xyz_1,   \
                             g_0_zzzzzz_0_xyzz_0,  \
                             g_0_zzzzzz_0_xyzz_1,  \
                             g_0_zzzzzz_0_xzz_1,   \
                             g_0_zzzzzz_0_xzzz_0,  \
                             g_0_zzzzzz_0_xzzz_1,  \
                             g_0_zzzzzz_0_yyy_1,   \
                             g_0_zzzzzz_0_yyyy_0,  \
                             g_0_zzzzzz_0_yyyy_1,  \
                             g_0_zzzzzz_0_yyyz_0,  \
                             g_0_zzzzzz_0_yyyz_1,  \
                             g_0_zzzzzz_0_yyz_1,   \
                             g_0_zzzzzz_0_yyzz_0,  \
                             g_0_zzzzzz_0_yyzz_1,  \
                             g_0_zzzzzz_0_yzz_1,   \
                             g_0_zzzzzz_0_yzzz_0,  \
                             g_0_zzzzzz_0_yzzz_1,  \
                             g_0_zzzzzz_0_zzz_1,   \
                             g_0_zzzzzz_0_zzzz_0,  \
                             g_0_zzzzzz_0_zzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_xxxx_0[i] = 4.0 * g_0_zzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxx_0[i] * pb_x + g_0_zzzzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxy_0[i] = 3.0 * g_0_zzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxy_0[i] * pb_x + g_0_zzzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxz_0[i] = 3.0 * g_0_zzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxz_0[i] * pb_x + g_0_zzzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyy_0[i] * pb_x + g_0_zzzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyz_0[i] = 2.0 * g_0_zzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyz_0[i] * pb_x + g_0_zzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxzz_0[i] = 2.0 * g_0_zzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzz_0[i] * pb_x + g_0_zzzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyy_0[i] = g_0_zzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyy_0[i] * pb_x + g_0_zzzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyz_0[i] = g_0_zzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyz_0[i] * pb_x + g_0_zzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyzz_0[i] = g_0_zzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzz_0[i] * pb_x + g_0_zzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xzzz_0[i] = g_0_zzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzz_0[i] * pb_x + g_0_zzzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyy_0[i] = g_0_zzzzzz_0_yyyy_0[i] * pb_x + g_0_zzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyz_0[i] = g_0_zzzzzz_0_yyyz_0[i] * pb_x + g_0_zzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyzz_0[i] = g_0_zzzzzz_0_yyzz_0[i] * pb_x + g_0_zzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yzzz_0[i] = g_0_zzzzzz_0_yzzz_0[i] * pb_x + g_0_zzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_zzzz_0[i] = g_0_zzzzzz_0_zzzz_0[i] * pb_x + g_0_zzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 420-435 components of targeted buffer : SKSG

    auto g_0_yyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 420);

    auto g_0_yyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 421);

    auto g_0_yyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 422);

    auto g_0_yyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 423);

    auto g_0_yyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 424);

    auto g_0_yyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 425);

    auto g_0_yyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 426);

    auto g_0_yyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 427);

    auto g_0_yyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 428);

    auto g_0_yyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 429);

    auto g_0_yyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 430);

    auto g_0_yyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 431);

    auto g_0_yyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 432);

    auto g_0_yyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 433);

    auto g_0_yyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 434);

#pragma omp simd aligned(g_0_yyyyy_0_xxxx_0,       \
                             g_0_yyyyy_0_xxxx_1,   \
                             g_0_yyyyy_0_xxxy_0,   \
                             g_0_yyyyy_0_xxxy_1,   \
                             g_0_yyyyy_0_xxxz_0,   \
                             g_0_yyyyy_0_xxxz_1,   \
                             g_0_yyyyy_0_xxyy_0,   \
                             g_0_yyyyy_0_xxyy_1,   \
                             g_0_yyyyy_0_xxyz_0,   \
                             g_0_yyyyy_0_xxyz_1,   \
                             g_0_yyyyy_0_xxzz_0,   \
                             g_0_yyyyy_0_xxzz_1,   \
                             g_0_yyyyy_0_xyyy_0,   \
                             g_0_yyyyy_0_xyyy_1,   \
                             g_0_yyyyy_0_xyyz_0,   \
                             g_0_yyyyy_0_xyyz_1,   \
                             g_0_yyyyy_0_xyzz_0,   \
                             g_0_yyyyy_0_xyzz_1,   \
                             g_0_yyyyy_0_xzzz_0,   \
                             g_0_yyyyy_0_xzzz_1,   \
                             g_0_yyyyy_0_yyyy_0,   \
                             g_0_yyyyy_0_yyyy_1,   \
                             g_0_yyyyy_0_yyyz_0,   \
                             g_0_yyyyy_0_yyyz_1,   \
                             g_0_yyyyy_0_yyzz_0,   \
                             g_0_yyyyy_0_yyzz_1,   \
                             g_0_yyyyy_0_yzzz_0,   \
                             g_0_yyyyy_0_yzzz_1,   \
                             g_0_yyyyy_0_zzzz_0,   \
                             g_0_yyyyy_0_zzzz_1,   \
                             g_0_yyyyyy_0_xxx_1,   \
                             g_0_yyyyyy_0_xxxx_0,  \
                             g_0_yyyyyy_0_xxxx_1,  \
                             g_0_yyyyyy_0_xxxy_0,  \
                             g_0_yyyyyy_0_xxxy_1,  \
                             g_0_yyyyyy_0_xxxz_0,  \
                             g_0_yyyyyy_0_xxxz_1,  \
                             g_0_yyyyyy_0_xxy_1,   \
                             g_0_yyyyyy_0_xxyy_0,  \
                             g_0_yyyyyy_0_xxyy_1,  \
                             g_0_yyyyyy_0_xxyz_0,  \
                             g_0_yyyyyy_0_xxyz_1,  \
                             g_0_yyyyyy_0_xxz_1,   \
                             g_0_yyyyyy_0_xxzz_0,  \
                             g_0_yyyyyy_0_xxzz_1,  \
                             g_0_yyyyyy_0_xyy_1,   \
                             g_0_yyyyyy_0_xyyy_0,  \
                             g_0_yyyyyy_0_xyyy_1,  \
                             g_0_yyyyyy_0_xyyz_0,  \
                             g_0_yyyyyy_0_xyyz_1,  \
                             g_0_yyyyyy_0_xyz_1,   \
                             g_0_yyyyyy_0_xyzz_0,  \
                             g_0_yyyyyy_0_xyzz_1,  \
                             g_0_yyyyyy_0_xzz_1,   \
                             g_0_yyyyyy_0_xzzz_0,  \
                             g_0_yyyyyy_0_xzzz_1,  \
                             g_0_yyyyyy_0_yyy_1,   \
                             g_0_yyyyyy_0_yyyy_0,  \
                             g_0_yyyyyy_0_yyyy_1,  \
                             g_0_yyyyyy_0_yyyz_0,  \
                             g_0_yyyyyy_0_yyyz_1,  \
                             g_0_yyyyyy_0_yyz_1,   \
                             g_0_yyyyyy_0_yyzz_0,  \
                             g_0_yyyyyy_0_yyzz_1,  \
                             g_0_yyyyyy_0_yzz_1,   \
                             g_0_yyyyyy_0_yzzz_0,  \
                             g_0_yyyyyy_0_yzzz_1,  \
                             g_0_yyyyyy_0_zzz_1,   \
                             g_0_yyyyyy_0_zzzz_0,  \
                             g_0_yyyyyy_0_zzzz_1,  \
                             g_0_yyyyyyy_0_xxxx_0, \
                             g_0_yyyyyyy_0_xxxy_0, \
                             g_0_yyyyyyy_0_xxxz_0, \
                             g_0_yyyyyyy_0_xxyy_0, \
                             g_0_yyyyyyy_0_xxyz_0, \
                             g_0_yyyyyyy_0_xxzz_0, \
                             g_0_yyyyyyy_0_xyyy_0, \
                             g_0_yyyyyyy_0_xyyz_0, \
                             g_0_yyyyyyy_0_xyzz_0, \
                             g_0_yyyyyyy_0_xzzz_0, \
                             g_0_yyyyyyy_0_yyyy_0, \
                             g_0_yyyyyyy_0_yyyz_0, \
                             g_0_yyyyyyy_0_yyzz_0, \
                             g_0_yyyyyyy_0_yzzz_0, \
                             g_0_yyyyyyy_0_zzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_xxxx_0[i] = 6.0 * g_0_yyyyy_0_xxxx_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxx_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxx_0[i] * pb_y +
                                  g_0_yyyyyy_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxy_0[i] = 6.0 * g_0_yyyyy_0_xxxy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxx_1[i] * fi_abcd_0 +
                                  g_0_yyyyyy_0_xxxy_0[i] * pb_y + g_0_yyyyyy_0_xxxy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxz_0[i] = 6.0 * g_0_yyyyy_0_xxxz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxz_0[i] * pb_y +
                                  g_0_yyyyyy_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyy_0[i] = 6.0 * g_0_yyyyy_0_xxyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyy_0[i] * pb_y + g_0_yyyyyy_0_xxyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyz_0[i] = 6.0 * g_0_yyyyy_0_xxyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxz_1[i] * fi_abcd_0 +
                                  g_0_yyyyyy_0_xxyz_0[i] * pb_y + g_0_yyyyyy_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxzz_0[i] = 6.0 * g_0_yyyyy_0_xxzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxzz_0[i] * pb_y +
                                  g_0_yyyyyy_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyy_0[i] = 6.0 * g_0_yyyyy_0_xyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyy_0[i] * pb_y + g_0_yyyyyy_0_xyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyz_0[i] = 6.0 * g_0_yyyyy_0_xyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyz_0[i] * pb_y + g_0_yyyyyy_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyzz_0[i] = 6.0 * g_0_yyyyy_0_xyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyyy_0_xyzz_0[i] * pb_y + g_0_yyyyyy_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xzzz_0[i] = 6.0 * g_0_yyyyy_0_xzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzzz_0[i] * pb_y +
                                  g_0_yyyyyy_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyy_0[i] = 6.0 * g_0_yyyyy_0_yyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyy_0[i] * pb_y + g_0_yyyyyy_0_yyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyz_0[i] = 6.0 * g_0_yyyyy_0_yyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyz_0[i] * pb_y + g_0_yyyyyy_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyzz_0[i] = 6.0 * g_0_yyyyy_0_yyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzz_0[i] * pb_y + g_0_yyyyyy_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yzzz_0[i] = 6.0 * g_0_yyyyy_0_yzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyyy_0_yzzz_0[i] * pb_y + g_0_yyyyyy_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_zzzz_0[i] = 6.0 * g_0_yyyyy_0_zzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzzz_0[i] * pb_y +
                                  g_0_yyyyyy_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 435-450 components of targeted buffer : SKSG

    auto g_0_yyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 435);

    auto g_0_yyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 436);

    auto g_0_yyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 437);

    auto g_0_yyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 438);

    auto g_0_yyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 439);

    auto g_0_yyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 440);

    auto g_0_yyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 441);

    auto g_0_yyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 442);

    auto g_0_yyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 443);

    auto g_0_yyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 444);

    auto g_0_yyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 445);

    auto g_0_yyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 446);

    auto g_0_yyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 447);

    auto g_0_yyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 448);

    auto g_0_yyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 449);

#pragma omp simd aligned(g_0_yyyyyy_0_xxx_1,       \
                             g_0_yyyyyy_0_xxxx_0,  \
                             g_0_yyyyyy_0_xxxx_1,  \
                             g_0_yyyyyy_0_xxxy_0,  \
                             g_0_yyyyyy_0_xxxy_1,  \
                             g_0_yyyyyy_0_xxxz_0,  \
                             g_0_yyyyyy_0_xxxz_1,  \
                             g_0_yyyyyy_0_xxy_1,   \
                             g_0_yyyyyy_0_xxyy_0,  \
                             g_0_yyyyyy_0_xxyy_1,  \
                             g_0_yyyyyy_0_xxyz_0,  \
                             g_0_yyyyyy_0_xxyz_1,  \
                             g_0_yyyyyy_0_xxz_1,   \
                             g_0_yyyyyy_0_xxzz_0,  \
                             g_0_yyyyyy_0_xxzz_1,  \
                             g_0_yyyyyy_0_xyy_1,   \
                             g_0_yyyyyy_0_xyyy_0,  \
                             g_0_yyyyyy_0_xyyy_1,  \
                             g_0_yyyyyy_0_xyyz_0,  \
                             g_0_yyyyyy_0_xyyz_1,  \
                             g_0_yyyyyy_0_xyz_1,   \
                             g_0_yyyyyy_0_xyzz_0,  \
                             g_0_yyyyyy_0_xyzz_1,  \
                             g_0_yyyyyy_0_xzz_1,   \
                             g_0_yyyyyy_0_xzzz_0,  \
                             g_0_yyyyyy_0_xzzz_1,  \
                             g_0_yyyyyy_0_yyy_1,   \
                             g_0_yyyyyy_0_yyyy_0,  \
                             g_0_yyyyyy_0_yyyy_1,  \
                             g_0_yyyyyy_0_yyyz_0,  \
                             g_0_yyyyyy_0_yyyz_1,  \
                             g_0_yyyyyy_0_yyz_1,   \
                             g_0_yyyyyy_0_yyzz_0,  \
                             g_0_yyyyyy_0_yyzz_1,  \
                             g_0_yyyyyy_0_yzz_1,   \
                             g_0_yyyyyy_0_yzzz_0,  \
                             g_0_yyyyyy_0_yzzz_1,  \
                             g_0_yyyyyy_0_zzz_1,   \
                             g_0_yyyyyy_0_zzzz_0,  \
                             g_0_yyyyyy_0_zzzz_1,  \
                             g_0_yyyyyyz_0_xxxx_0, \
                             g_0_yyyyyyz_0_xxxy_0, \
                             g_0_yyyyyyz_0_xxxz_0, \
                             g_0_yyyyyyz_0_xxyy_0, \
                             g_0_yyyyyyz_0_xxyz_0, \
                             g_0_yyyyyyz_0_xxzz_0, \
                             g_0_yyyyyyz_0_xyyy_0, \
                             g_0_yyyyyyz_0_xyyz_0, \
                             g_0_yyyyyyz_0_xyzz_0, \
                             g_0_yyyyyyz_0_xzzz_0, \
                             g_0_yyyyyyz_0_yyyy_0, \
                             g_0_yyyyyyz_0_yyyz_0, \
                             g_0_yyyyyyz_0_yyzz_0, \
                             g_0_yyyyyyz_0_yzzz_0, \
                             g_0_yyyyyyz_0_zzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_xxxx_0[i] = g_0_yyyyyy_0_xxxx_0[i] * pb_z + g_0_yyyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxy_0[i] = g_0_yyyyyy_0_xxxy_0[i] * pb_z + g_0_yyyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxz_0[i] = g_0_yyyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxz_0[i] * pb_z + g_0_yyyyyy_0_xxxz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyy_0[i] = g_0_yyyyyy_0_xxyy_0[i] * pb_z + g_0_yyyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyz_0[i] = g_0_yyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyz_0[i] * pb_z + g_0_yyyyyy_0_xxyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzz_0[i] * pb_z + g_0_yyyyyy_0_xxzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyy_0[i] = g_0_yyyyyy_0_xyyy_0[i] * pb_z + g_0_yyyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyz_0[i] = g_0_yyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyz_0[i] * pb_z + g_0_yyyyyy_0_xyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzz_0[i] * pb_z + g_0_yyyyyy_0_xyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xzzz_0[i] = 3.0 * g_0_yyyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzz_0[i] * pb_z + g_0_yyyyyy_0_xzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyy_0[i] = g_0_yyyyyy_0_yyyy_0[i] * pb_z + g_0_yyyyyy_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyz_0[i] = g_0_yyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyz_0[i] * pb_z + g_0_yyyyyy_0_yyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyzz_0[i] = 2.0 * g_0_yyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzz_0[i] * pb_z + g_0_yyyyyy_0_yyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yzzz_0[i] = 3.0 * g_0_yyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzz_0[i] * pb_z + g_0_yyyyyy_0_yzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_zzzz_0[i] = 4.0 * g_0_yyyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_zzzz_0[i] * pb_z + g_0_yyyyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 450-465 components of targeted buffer : SKSG

    auto g_0_yyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 450);

    auto g_0_yyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 451);

    auto g_0_yyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 452);

    auto g_0_yyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 453);

    auto g_0_yyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 454);

    auto g_0_yyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 455);

    auto g_0_yyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 456);

    auto g_0_yyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 457);

    auto g_0_yyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 458);

    auto g_0_yyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 459);

    auto g_0_yyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 460);

    auto g_0_yyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 461);

    auto g_0_yyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 462);

    auto g_0_yyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 463);

    auto g_0_yyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 464);

#pragma omp simd aligned(g_0_yyyyy_0_xxxy_0,       \
                             g_0_yyyyy_0_xxxy_1,   \
                             g_0_yyyyy_0_xxyy_0,   \
                             g_0_yyyyy_0_xxyy_1,   \
                             g_0_yyyyy_0_xyyy_0,   \
                             g_0_yyyyy_0_xyyy_1,   \
                             g_0_yyyyy_0_yyyy_0,   \
                             g_0_yyyyy_0_yyyy_1,   \
                             g_0_yyyyyz_0_xxxy_0,  \
                             g_0_yyyyyz_0_xxxy_1,  \
                             g_0_yyyyyz_0_xxyy_0,  \
                             g_0_yyyyyz_0_xxyy_1,  \
                             g_0_yyyyyz_0_xyyy_0,  \
                             g_0_yyyyyz_0_xyyy_1,  \
                             g_0_yyyyyz_0_yyyy_0,  \
                             g_0_yyyyyz_0_yyyy_1,  \
                             g_0_yyyyyzz_0_xxxx_0, \
                             g_0_yyyyyzz_0_xxxy_0, \
                             g_0_yyyyyzz_0_xxxz_0, \
                             g_0_yyyyyzz_0_xxyy_0, \
                             g_0_yyyyyzz_0_xxyz_0, \
                             g_0_yyyyyzz_0_xxzz_0, \
                             g_0_yyyyyzz_0_xyyy_0, \
                             g_0_yyyyyzz_0_xyyz_0, \
                             g_0_yyyyyzz_0_xyzz_0, \
                             g_0_yyyyyzz_0_xzzz_0, \
                             g_0_yyyyyzz_0_yyyy_0, \
                             g_0_yyyyyzz_0_yyyz_0, \
                             g_0_yyyyyzz_0_yyzz_0, \
                             g_0_yyyyyzz_0_yzzz_0, \
                             g_0_yyyyyzz_0_zzzz_0, \
                             g_0_yyyyzz_0_xxxx_0,  \
                             g_0_yyyyzz_0_xxxx_1,  \
                             g_0_yyyyzz_0_xxxz_0,  \
                             g_0_yyyyzz_0_xxxz_1,  \
                             g_0_yyyyzz_0_xxyz_0,  \
                             g_0_yyyyzz_0_xxyz_1,  \
                             g_0_yyyyzz_0_xxz_1,   \
                             g_0_yyyyzz_0_xxzz_0,  \
                             g_0_yyyyzz_0_xxzz_1,  \
                             g_0_yyyyzz_0_xyyz_0,  \
                             g_0_yyyyzz_0_xyyz_1,  \
                             g_0_yyyyzz_0_xyz_1,   \
                             g_0_yyyyzz_0_xyzz_0,  \
                             g_0_yyyyzz_0_xyzz_1,  \
                             g_0_yyyyzz_0_xzz_1,   \
                             g_0_yyyyzz_0_xzzz_0,  \
                             g_0_yyyyzz_0_xzzz_1,  \
                             g_0_yyyyzz_0_yyyz_0,  \
                             g_0_yyyyzz_0_yyyz_1,  \
                             g_0_yyyyzz_0_yyz_1,   \
                             g_0_yyyyzz_0_yyzz_0,  \
                             g_0_yyyyzz_0_yyzz_1,  \
                             g_0_yyyyzz_0_yzz_1,   \
                             g_0_yyyyzz_0_yzzz_0,  \
                             g_0_yyyyzz_0_yzzz_1,  \
                             g_0_yyyyzz_0_zzz_1,   \
                             g_0_yyyyzz_0_zzzz_0,  \
                             g_0_yyyyzz_0_zzzz_1,  \
                             g_0_yyyzz_0_xxxx_0,   \
                             g_0_yyyzz_0_xxxx_1,   \
                             g_0_yyyzz_0_xxxz_0,   \
                             g_0_yyyzz_0_xxxz_1,   \
                             g_0_yyyzz_0_xxyz_0,   \
                             g_0_yyyzz_0_xxyz_1,   \
                             g_0_yyyzz_0_xxzz_0,   \
                             g_0_yyyzz_0_xxzz_1,   \
                             g_0_yyyzz_0_xyyz_0,   \
                             g_0_yyyzz_0_xyyz_1,   \
                             g_0_yyyzz_0_xyzz_0,   \
                             g_0_yyyzz_0_xyzz_1,   \
                             g_0_yyyzz_0_xzzz_0,   \
                             g_0_yyyzz_0_xzzz_1,   \
                             g_0_yyyzz_0_yyyz_0,   \
                             g_0_yyyzz_0_yyyz_1,   \
                             g_0_yyyzz_0_yyzz_0,   \
                             g_0_yyyzz_0_yyzz_1,   \
                             g_0_yyyzz_0_yzzz_0,   \
                             g_0_yyyzz_0_yzzz_1,   \
                             g_0_yyyzz_0_zzzz_0,   \
                             g_0_yyyzz_0_zzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_xxxx_0[i] = 4.0 * g_0_yyyzz_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxx_0[i] * pb_y +
                                  g_0_yyyyzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxy_0[i] =
            g_0_yyyyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxy_0[i] * pb_z + g_0_yyyyyz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxz_0[i] = 4.0 * g_0_yyyzz_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxz_0[i] * pb_y +
                                  g_0_yyyyzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyy_0[i] =
            g_0_yyyyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxyy_0[i] * pb_z + g_0_yyyyyz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxyz_0[i] = 4.0 * g_0_yyyzz_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxz_1[i] * fi_abcd_0 +
                                  g_0_yyyyzz_0_xxyz_0[i] * pb_y + g_0_yyyyzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxzz_0[i] = 4.0 * g_0_yyyzz_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxzz_0[i] * pb_y +
                                  g_0_yyyyzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyy_0[i] =
            g_0_yyyyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xyyy_0[i] * pb_z + g_0_yyyyyz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xyyz_0[i] = 4.0 * g_0_yyyzz_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyz_0[i] * pb_y + g_0_yyyyzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyzz_0[i] = 4.0 * g_0_yyyzz_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyzz_0_xyzz_0[i] * pb_y + g_0_yyyyzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xzzz_0[i] = 4.0 * g_0_yyyzz_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzzz_0[i] * pb_y +
                                  g_0_yyyyzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyy_0[i] =
            g_0_yyyyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_yyyy_0[i] * pb_z + g_0_yyyyyz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_yyyz_0[i] = 4.0 * g_0_yyyzz_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyz_0[i] * pb_y + g_0_yyyyzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyzz_0[i] = 4.0 * g_0_yyyzz_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyzz_0[i] * pb_y + g_0_yyyyzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yzzz_0[i] = 4.0 * g_0_yyyzz_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_yyyyzz_0_yzzz_0[i] * pb_y + g_0_yyyyzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_zzzz_0[i] = 4.0 * g_0_yyyzz_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzzz_0[i] * pb_y +
                                  g_0_yyyyzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 465-480 components of targeted buffer : SKSG

    auto g_0_yyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 465);

    auto g_0_yyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 466);

    auto g_0_yyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 467);

    auto g_0_yyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 468);

    auto g_0_yyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 469);

    auto g_0_yyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 470);

    auto g_0_yyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 471);

    auto g_0_yyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 472);

    auto g_0_yyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 473);

    auto g_0_yyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 474);

    auto g_0_yyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 475);

    auto g_0_yyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 476);

    auto g_0_yyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 477);

    auto g_0_yyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 478);

    auto g_0_yyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 479);

#pragma omp simd aligned(g_0_yyyyz_0_xxxy_0,       \
                             g_0_yyyyz_0_xxxy_1,   \
                             g_0_yyyyz_0_xxyy_0,   \
                             g_0_yyyyz_0_xxyy_1,   \
                             g_0_yyyyz_0_xyyy_0,   \
                             g_0_yyyyz_0_xyyy_1,   \
                             g_0_yyyyz_0_yyyy_0,   \
                             g_0_yyyyz_0_yyyy_1,   \
                             g_0_yyyyzz_0_xxxy_0,  \
                             g_0_yyyyzz_0_xxxy_1,  \
                             g_0_yyyyzz_0_xxyy_0,  \
                             g_0_yyyyzz_0_xxyy_1,  \
                             g_0_yyyyzz_0_xyyy_0,  \
                             g_0_yyyyzz_0_xyyy_1,  \
                             g_0_yyyyzz_0_yyyy_0,  \
                             g_0_yyyyzz_0_yyyy_1,  \
                             g_0_yyyyzzz_0_xxxx_0, \
                             g_0_yyyyzzz_0_xxxy_0, \
                             g_0_yyyyzzz_0_xxxz_0, \
                             g_0_yyyyzzz_0_xxyy_0, \
                             g_0_yyyyzzz_0_xxyz_0, \
                             g_0_yyyyzzz_0_xxzz_0, \
                             g_0_yyyyzzz_0_xyyy_0, \
                             g_0_yyyyzzz_0_xyyz_0, \
                             g_0_yyyyzzz_0_xyzz_0, \
                             g_0_yyyyzzz_0_xzzz_0, \
                             g_0_yyyyzzz_0_yyyy_0, \
                             g_0_yyyyzzz_0_yyyz_0, \
                             g_0_yyyyzzz_0_yyzz_0, \
                             g_0_yyyyzzz_0_yzzz_0, \
                             g_0_yyyyzzz_0_zzzz_0, \
                             g_0_yyyzzz_0_xxxx_0,  \
                             g_0_yyyzzz_0_xxxx_1,  \
                             g_0_yyyzzz_0_xxxz_0,  \
                             g_0_yyyzzz_0_xxxz_1,  \
                             g_0_yyyzzz_0_xxyz_0,  \
                             g_0_yyyzzz_0_xxyz_1,  \
                             g_0_yyyzzz_0_xxz_1,   \
                             g_0_yyyzzz_0_xxzz_0,  \
                             g_0_yyyzzz_0_xxzz_1,  \
                             g_0_yyyzzz_0_xyyz_0,  \
                             g_0_yyyzzz_0_xyyz_1,  \
                             g_0_yyyzzz_0_xyz_1,   \
                             g_0_yyyzzz_0_xyzz_0,  \
                             g_0_yyyzzz_0_xyzz_1,  \
                             g_0_yyyzzz_0_xzz_1,   \
                             g_0_yyyzzz_0_xzzz_0,  \
                             g_0_yyyzzz_0_xzzz_1,  \
                             g_0_yyyzzz_0_yyyz_0,  \
                             g_0_yyyzzz_0_yyyz_1,  \
                             g_0_yyyzzz_0_yyz_1,   \
                             g_0_yyyzzz_0_yyzz_0,  \
                             g_0_yyyzzz_0_yyzz_1,  \
                             g_0_yyyzzz_0_yzz_1,   \
                             g_0_yyyzzz_0_yzzz_0,  \
                             g_0_yyyzzz_0_yzzz_1,  \
                             g_0_yyyzzz_0_zzz_1,   \
                             g_0_yyyzzz_0_zzzz_0,  \
                             g_0_yyyzzz_0_zzzz_1,  \
                             g_0_yyzzz_0_xxxx_0,   \
                             g_0_yyzzz_0_xxxx_1,   \
                             g_0_yyzzz_0_xxxz_0,   \
                             g_0_yyzzz_0_xxxz_1,   \
                             g_0_yyzzz_0_xxyz_0,   \
                             g_0_yyzzz_0_xxyz_1,   \
                             g_0_yyzzz_0_xxzz_0,   \
                             g_0_yyzzz_0_xxzz_1,   \
                             g_0_yyzzz_0_xyyz_0,   \
                             g_0_yyzzz_0_xyyz_1,   \
                             g_0_yyzzz_0_xyzz_0,   \
                             g_0_yyzzz_0_xyzz_1,   \
                             g_0_yyzzz_0_xzzz_0,   \
                             g_0_yyzzz_0_xzzz_1,   \
                             g_0_yyzzz_0_yyyz_0,   \
                             g_0_yyzzz_0_yyyz_1,   \
                             g_0_yyzzz_0_yyzz_0,   \
                             g_0_yyzzz_0_yyzz_1,   \
                             g_0_yyzzz_0_yzzz_0,   \
                             g_0_yyzzz_0_yzzz_1,   \
                             g_0_yyzzz_0_zzzz_0,   \
                             g_0_yyzzz_0_zzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_xxxx_0[i] = 3.0 * g_0_yyzzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxx_0[i] * pb_y +
                                  g_0_yyyzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxy_0[i] = 2.0 * g_0_yyyyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxy_0[i] * pb_z +
                                  g_0_yyyyzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxz_0[i] = 3.0 * g_0_yyzzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxz_0[i] * pb_y +
                                  g_0_yyyzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyy_0[i] = 2.0 * g_0_yyyyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxyy_0[i] * pb_z +
                                  g_0_yyyyzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxyz_0[i] = 3.0 * g_0_yyzzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxz_1[i] * fi_abcd_0 +
                                  g_0_yyyzzz_0_xxyz_0[i] * pb_y + g_0_yyyzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxzz_0[i] = 3.0 * g_0_yyzzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxzz_0[i] * pb_y +
                                  g_0_yyyzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyy_0[i] = 2.0 * g_0_yyyyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xyyy_0[i] * pb_z +
                                  g_0_yyyyzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xyyz_0[i] = 3.0 * g_0_yyzzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyz_0[i] * pb_y + g_0_yyyzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyzz_0[i] = 3.0 * g_0_yyzzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzz_1[i] * fi_abcd_0 +
                                  g_0_yyyzzz_0_xyzz_0[i] * pb_y + g_0_yyyzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xzzz_0[i] = 3.0 * g_0_yyzzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzzz_0[i] * pb_y +
                                  g_0_yyyzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyy_0[i] = 2.0 * g_0_yyyyz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_yyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_yyyy_0[i] * pb_z +
                                  g_0_yyyyzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_yyyz_0[i] = 3.0 * g_0_yyzzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyz_0[i] * pb_y + g_0_yyyzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyzz_0[i] = 3.0 * g_0_yyzzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyzz_0[i] * pb_y + g_0_yyyzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yzzz_0[i] = 3.0 * g_0_yyzzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_yyyzzz_0_yzzz_0[i] * pb_y + g_0_yyyzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_zzzz_0[i] = 3.0 * g_0_yyzzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzzz_0[i] * pb_y +
                                  g_0_yyyzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 480-495 components of targeted buffer : SKSG

    auto g_0_yyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 480);

    auto g_0_yyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 481);

    auto g_0_yyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 482);

    auto g_0_yyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 483);

    auto g_0_yyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 484);

    auto g_0_yyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 485);

    auto g_0_yyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 486);

    auto g_0_yyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 487);

    auto g_0_yyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 488);

    auto g_0_yyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 489);

    auto g_0_yyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 490);

    auto g_0_yyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 491);

    auto g_0_yyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 492);

    auto g_0_yyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 493);

    auto g_0_yyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 494);

#pragma omp simd aligned(g_0_yyyzz_0_xxxy_0,       \
                             g_0_yyyzz_0_xxxy_1,   \
                             g_0_yyyzz_0_xxyy_0,   \
                             g_0_yyyzz_0_xxyy_1,   \
                             g_0_yyyzz_0_xyyy_0,   \
                             g_0_yyyzz_0_xyyy_1,   \
                             g_0_yyyzz_0_yyyy_0,   \
                             g_0_yyyzz_0_yyyy_1,   \
                             g_0_yyyzzz_0_xxxy_0,  \
                             g_0_yyyzzz_0_xxxy_1,  \
                             g_0_yyyzzz_0_xxyy_0,  \
                             g_0_yyyzzz_0_xxyy_1,  \
                             g_0_yyyzzz_0_xyyy_0,  \
                             g_0_yyyzzz_0_xyyy_1,  \
                             g_0_yyyzzz_0_yyyy_0,  \
                             g_0_yyyzzz_0_yyyy_1,  \
                             g_0_yyyzzzz_0_xxxx_0, \
                             g_0_yyyzzzz_0_xxxy_0, \
                             g_0_yyyzzzz_0_xxxz_0, \
                             g_0_yyyzzzz_0_xxyy_0, \
                             g_0_yyyzzzz_0_xxyz_0, \
                             g_0_yyyzzzz_0_xxzz_0, \
                             g_0_yyyzzzz_0_xyyy_0, \
                             g_0_yyyzzzz_0_xyyz_0, \
                             g_0_yyyzzzz_0_xyzz_0, \
                             g_0_yyyzzzz_0_xzzz_0, \
                             g_0_yyyzzzz_0_yyyy_0, \
                             g_0_yyyzzzz_0_yyyz_0, \
                             g_0_yyyzzzz_0_yyzz_0, \
                             g_0_yyyzzzz_0_yzzz_0, \
                             g_0_yyyzzzz_0_zzzz_0, \
                             g_0_yyzzzz_0_xxxx_0,  \
                             g_0_yyzzzz_0_xxxx_1,  \
                             g_0_yyzzzz_0_xxxz_0,  \
                             g_0_yyzzzz_0_xxxz_1,  \
                             g_0_yyzzzz_0_xxyz_0,  \
                             g_0_yyzzzz_0_xxyz_1,  \
                             g_0_yyzzzz_0_xxz_1,   \
                             g_0_yyzzzz_0_xxzz_0,  \
                             g_0_yyzzzz_0_xxzz_1,  \
                             g_0_yyzzzz_0_xyyz_0,  \
                             g_0_yyzzzz_0_xyyz_1,  \
                             g_0_yyzzzz_0_xyz_1,   \
                             g_0_yyzzzz_0_xyzz_0,  \
                             g_0_yyzzzz_0_xyzz_1,  \
                             g_0_yyzzzz_0_xzz_1,   \
                             g_0_yyzzzz_0_xzzz_0,  \
                             g_0_yyzzzz_0_xzzz_1,  \
                             g_0_yyzzzz_0_yyyz_0,  \
                             g_0_yyzzzz_0_yyyz_1,  \
                             g_0_yyzzzz_0_yyz_1,   \
                             g_0_yyzzzz_0_yyzz_0,  \
                             g_0_yyzzzz_0_yyzz_1,  \
                             g_0_yyzzzz_0_yzz_1,   \
                             g_0_yyzzzz_0_yzzz_0,  \
                             g_0_yyzzzz_0_yzzz_1,  \
                             g_0_yyzzzz_0_zzz_1,   \
                             g_0_yyzzzz_0_zzzz_0,  \
                             g_0_yyzzzz_0_zzzz_1,  \
                             g_0_yzzzz_0_xxxx_0,   \
                             g_0_yzzzz_0_xxxx_1,   \
                             g_0_yzzzz_0_xxxz_0,   \
                             g_0_yzzzz_0_xxxz_1,   \
                             g_0_yzzzz_0_xxyz_0,   \
                             g_0_yzzzz_0_xxyz_1,   \
                             g_0_yzzzz_0_xxzz_0,   \
                             g_0_yzzzz_0_xxzz_1,   \
                             g_0_yzzzz_0_xyyz_0,   \
                             g_0_yzzzz_0_xyyz_1,   \
                             g_0_yzzzz_0_xyzz_0,   \
                             g_0_yzzzz_0_xyzz_1,   \
                             g_0_yzzzz_0_xzzz_0,   \
                             g_0_yzzzz_0_xzzz_1,   \
                             g_0_yzzzz_0_yyyz_0,   \
                             g_0_yzzzz_0_yyyz_1,   \
                             g_0_yzzzz_0_yyzz_0,   \
                             g_0_yzzzz_0_yyzz_1,   \
                             g_0_yzzzz_0_yzzz_0,   \
                             g_0_yzzzz_0_yzzz_1,   \
                             g_0_yzzzz_0_zzzz_0,   \
                             g_0_yzzzz_0_zzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_xxxx_0[i] = 2.0 * g_0_yzzzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxx_0[i] * pb_y +
                                  g_0_yyzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxy_0[i] = 3.0 * g_0_yyyzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxy_0[i] * pb_z +
                                  g_0_yyyzzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxz_0[i] = 2.0 * g_0_yzzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxz_0[i] * pb_y +
                                  g_0_yyzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyy_0[i] = 3.0 * g_0_yyyzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxyy_0[i] * pb_z +
                                  g_0_yyyzzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxyz_0[i] = 2.0 * g_0_yzzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxz_1[i] * fi_abcd_0 +
                                  g_0_yyzzzz_0_xxyz_0[i] * pb_y + g_0_yyzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxzz_0[i] = 2.0 * g_0_yzzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxzz_0[i] * pb_y +
                                  g_0_yyzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyy_0[i] = 3.0 * g_0_yyyzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xyyy_0[i] * pb_z +
                                  g_0_yyyzzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyz_0[i] * pb_y + g_0_yyzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyzz_0[i] = 2.0 * g_0_yzzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzz_1[i] * fi_abcd_0 +
                                  g_0_yyzzzz_0_xyzz_0[i] * pb_y + g_0_yyzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xzzz_0[i] = 2.0 * g_0_yzzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzzz_0[i] * pb_y +
                                  g_0_yyzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyy_0[i] = 3.0 * g_0_yyyzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_yyyy_0[i] * pb_z +
                                  g_0_yyyzzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_yyyz_0[i] = 2.0 * g_0_yzzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyz_0[i] * pb_y + g_0_yyzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyzz_0[i] = 2.0 * g_0_yzzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyzz_0[i] * pb_y + g_0_yyzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yzzz_0[i] = 2.0 * g_0_yzzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_yyzzzz_0_yzzz_0[i] * pb_y + g_0_yyzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_zzzz_0[i] = 2.0 * g_0_yzzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzzz_0[i] * pb_y +
                                  g_0_yyzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 495-510 components of targeted buffer : SKSG

    auto g_0_yyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 495);

    auto g_0_yyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 496);

    auto g_0_yyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 497);

    auto g_0_yyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 498);

    auto g_0_yyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 499);

    auto g_0_yyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 500);

    auto g_0_yyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 501);

    auto g_0_yyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 502);

    auto g_0_yyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 503);

    auto g_0_yyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 504);

    auto g_0_yyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 505);

    auto g_0_yyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 506);

    auto g_0_yyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 507);

    auto g_0_yyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 508);

    auto g_0_yyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 509);

#pragma omp simd aligned(g_0_yyzzz_0_xxxy_0,       \
                             g_0_yyzzz_0_xxxy_1,   \
                             g_0_yyzzz_0_xxyy_0,   \
                             g_0_yyzzz_0_xxyy_1,   \
                             g_0_yyzzz_0_xyyy_0,   \
                             g_0_yyzzz_0_xyyy_1,   \
                             g_0_yyzzz_0_yyyy_0,   \
                             g_0_yyzzz_0_yyyy_1,   \
                             g_0_yyzzzz_0_xxxy_0,  \
                             g_0_yyzzzz_0_xxxy_1,  \
                             g_0_yyzzzz_0_xxyy_0,  \
                             g_0_yyzzzz_0_xxyy_1,  \
                             g_0_yyzzzz_0_xyyy_0,  \
                             g_0_yyzzzz_0_xyyy_1,  \
                             g_0_yyzzzz_0_yyyy_0,  \
                             g_0_yyzzzz_0_yyyy_1,  \
                             g_0_yyzzzzz_0_xxxx_0, \
                             g_0_yyzzzzz_0_xxxy_0, \
                             g_0_yyzzzzz_0_xxxz_0, \
                             g_0_yyzzzzz_0_xxyy_0, \
                             g_0_yyzzzzz_0_xxyz_0, \
                             g_0_yyzzzzz_0_xxzz_0, \
                             g_0_yyzzzzz_0_xyyy_0, \
                             g_0_yyzzzzz_0_xyyz_0, \
                             g_0_yyzzzzz_0_xyzz_0, \
                             g_0_yyzzzzz_0_xzzz_0, \
                             g_0_yyzzzzz_0_yyyy_0, \
                             g_0_yyzzzzz_0_yyyz_0, \
                             g_0_yyzzzzz_0_yyzz_0, \
                             g_0_yyzzzzz_0_yzzz_0, \
                             g_0_yyzzzzz_0_zzzz_0, \
                             g_0_yzzzzz_0_xxxx_0,  \
                             g_0_yzzzzz_0_xxxx_1,  \
                             g_0_yzzzzz_0_xxxz_0,  \
                             g_0_yzzzzz_0_xxxz_1,  \
                             g_0_yzzzzz_0_xxyz_0,  \
                             g_0_yzzzzz_0_xxyz_1,  \
                             g_0_yzzzzz_0_xxz_1,   \
                             g_0_yzzzzz_0_xxzz_0,  \
                             g_0_yzzzzz_0_xxzz_1,  \
                             g_0_yzzzzz_0_xyyz_0,  \
                             g_0_yzzzzz_0_xyyz_1,  \
                             g_0_yzzzzz_0_xyz_1,   \
                             g_0_yzzzzz_0_xyzz_0,  \
                             g_0_yzzzzz_0_xyzz_1,  \
                             g_0_yzzzzz_0_xzz_1,   \
                             g_0_yzzzzz_0_xzzz_0,  \
                             g_0_yzzzzz_0_xzzz_1,  \
                             g_0_yzzzzz_0_yyyz_0,  \
                             g_0_yzzzzz_0_yyyz_1,  \
                             g_0_yzzzzz_0_yyz_1,   \
                             g_0_yzzzzz_0_yyzz_0,  \
                             g_0_yzzzzz_0_yyzz_1,  \
                             g_0_yzzzzz_0_yzz_1,   \
                             g_0_yzzzzz_0_yzzz_0,  \
                             g_0_yzzzzz_0_yzzz_1,  \
                             g_0_yzzzzz_0_zzz_1,   \
                             g_0_yzzzzz_0_zzzz_0,  \
                             g_0_yzzzzz_0_zzzz_1,  \
                             g_0_zzzzz_0_xxxx_0,   \
                             g_0_zzzzz_0_xxxx_1,   \
                             g_0_zzzzz_0_xxxz_0,   \
                             g_0_zzzzz_0_xxxz_1,   \
                             g_0_zzzzz_0_xxyz_0,   \
                             g_0_zzzzz_0_xxyz_1,   \
                             g_0_zzzzz_0_xxzz_0,   \
                             g_0_zzzzz_0_xxzz_1,   \
                             g_0_zzzzz_0_xyyz_0,   \
                             g_0_zzzzz_0_xyyz_1,   \
                             g_0_zzzzz_0_xyzz_0,   \
                             g_0_zzzzz_0_xyzz_1,   \
                             g_0_zzzzz_0_xzzz_0,   \
                             g_0_zzzzz_0_xzzz_1,   \
                             g_0_zzzzz_0_yyyz_0,   \
                             g_0_zzzzz_0_yyyz_1,   \
                             g_0_zzzzz_0_yyzz_0,   \
                             g_0_zzzzz_0_yyzz_1,   \
                             g_0_zzzzz_0_yzzz_0,   \
                             g_0_zzzzz_0_yzzz_1,   \
                             g_0_zzzzz_0_zzzz_0,   \
                             g_0_zzzzz_0_zzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_xxxx_0[i] =
            g_0_zzzzz_0_xxxx_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxx_0[i] * pb_y + g_0_yzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxy_0[i] = 4.0 * g_0_yyzzz_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxy_0[i] * pb_z +
                                  g_0_yyzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxz_0[i] =
            g_0_zzzzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxz_0[i] * pb_y + g_0_yzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyy_0[i] = 4.0 * g_0_yyzzz_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxyy_0[i] * pb_z +
                                  g_0_yyzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxyz_0[i] = g_0_zzzzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxz_1[i] * fi_abcd_0 +
                                  g_0_yzzzzz_0_xxyz_0[i] * pb_y + g_0_yzzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxzz_0[i] =
            g_0_zzzzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzz_0[i] * pb_y + g_0_yzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyy_0[i] = 4.0 * g_0_yyzzz_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xyyy_0[i] * pb_z +
                                  g_0_yyzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xyyz_0[i] = g_0_zzzzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xyz_1[i] * fi_abcd_0 +
                                  g_0_yzzzzz_0_xyyz_0[i] * pb_y + g_0_yzzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyzz_0[i] = g_0_zzzzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzzz_0_xyzz_0[i] * pb_y + g_0_yzzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xzzz_0[i] =
            g_0_zzzzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzz_0[i] * pb_y + g_0_yzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyy_0[i] = 4.0 * g_0_yyzzz_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_yyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_yyyy_0[i] * pb_z +
                                  g_0_yyzzzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_yyyz_0[i] = g_0_zzzzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_yyz_1[i] * fi_abcd_0 +
                                  g_0_yzzzzz_0_yyyz_0[i] * pb_y + g_0_yzzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyzz_0[i] = g_0_zzzzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_yzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzzz_0_yyzz_0[i] * pb_y + g_0_yzzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yzzz_0[i] = g_0_zzzzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzz_1[i] * fi_abcd_0 +
                                  g_0_yzzzzz_0_yzzz_0[i] * pb_y + g_0_yzzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_zzzz_0[i] =
            g_0_zzzzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzz_0[i] * pb_y + g_0_yzzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 510-525 components of targeted buffer : SKSG

    auto g_0_yzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 510);

    auto g_0_yzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 511);

    auto g_0_yzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 512);

    auto g_0_yzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 513);

    auto g_0_yzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 514);

    auto g_0_yzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 515);

    auto g_0_yzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 516);

    auto g_0_yzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 517);

    auto g_0_yzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 518);

    auto g_0_yzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 519);

    auto g_0_yzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 520);

    auto g_0_yzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 521);

    auto g_0_yzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 522);

    auto g_0_yzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 523);

    auto g_0_yzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 524);

#pragma omp simd aligned(g_0_yzzzzzz_0_xxxx_0,     \
                             g_0_yzzzzzz_0_xxxy_0, \
                             g_0_yzzzzzz_0_xxxz_0, \
                             g_0_yzzzzzz_0_xxyy_0, \
                             g_0_yzzzzzz_0_xxyz_0, \
                             g_0_yzzzzzz_0_xxzz_0, \
                             g_0_yzzzzzz_0_xyyy_0, \
                             g_0_yzzzzzz_0_xyyz_0, \
                             g_0_yzzzzzz_0_xyzz_0, \
                             g_0_yzzzzzz_0_xzzz_0, \
                             g_0_yzzzzzz_0_yyyy_0, \
                             g_0_yzzzzzz_0_yyyz_0, \
                             g_0_yzzzzzz_0_yyzz_0, \
                             g_0_yzzzzzz_0_yzzz_0, \
                             g_0_yzzzzzz_0_zzzz_0, \
                             g_0_zzzzzz_0_xxx_1,   \
                             g_0_zzzzzz_0_xxxx_0,  \
                             g_0_zzzzzz_0_xxxx_1,  \
                             g_0_zzzzzz_0_xxxy_0,  \
                             g_0_zzzzzz_0_xxxy_1,  \
                             g_0_zzzzzz_0_xxxz_0,  \
                             g_0_zzzzzz_0_xxxz_1,  \
                             g_0_zzzzzz_0_xxy_1,   \
                             g_0_zzzzzz_0_xxyy_0,  \
                             g_0_zzzzzz_0_xxyy_1,  \
                             g_0_zzzzzz_0_xxyz_0,  \
                             g_0_zzzzzz_0_xxyz_1,  \
                             g_0_zzzzzz_0_xxz_1,   \
                             g_0_zzzzzz_0_xxzz_0,  \
                             g_0_zzzzzz_0_xxzz_1,  \
                             g_0_zzzzzz_0_xyy_1,   \
                             g_0_zzzzzz_0_xyyy_0,  \
                             g_0_zzzzzz_0_xyyy_1,  \
                             g_0_zzzzzz_0_xyyz_0,  \
                             g_0_zzzzzz_0_xyyz_1,  \
                             g_0_zzzzzz_0_xyz_1,   \
                             g_0_zzzzzz_0_xyzz_0,  \
                             g_0_zzzzzz_0_xyzz_1,  \
                             g_0_zzzzzz_0_xzz_1,   \
                             g_0_zzzzzz_0_xzzz_0,  \
                             g_0_zzzzzz_0_xzzz_1,  \
                             g_0_zzzzzz_0_yyy_1,   \
                             g_0_zzzzzz_0_yyyy_0,  \
                             g_0_zzzzzz_0_yyyy_1,  \
                             g_0_zzzzzz_0_yyyz_0,  \
                             g_0_zzzzzz_0_yyyz_1,  \
                             g_0_zzzzzz_0_yyz_1,   \
                             g_0_zzzzzz_0_yyzz_0,  \
                             g_0_zzzzzz_0_yyzz_1,  \
                             g_0_zzzzzz_0_yzz_1,   \
                             g_0_zzzzzz_0_yzzz_0,  \
                             g_0_zzzzzz_0_yzzz_1,  \
                             g_0_zzzzzz_0_zzz_1,   \
                             g_0_zzzzzz_0_zzzz_0,  \
                             g_0_zzzzzz_0_zzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_xxxx_0[i] = g_0_zzzzzz_0_xxxx_0[i] * pb_y + g_0_zzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxy_0[i] = g_0_zzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxy_0[i] * pb_y + g_0_zzzzzz_0_xxxy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxz_0[i] = g_0_zzzzzz_0_xxxz_0[i] * pb_y + g_0_zzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyy_0[i] * pb_y + g_0_zzzzzz_0_xxyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyz_0[i] = g_0_zzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyz_0[i] * pb_y + g_0_zzzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxzz_0[i] = g_0_zzzzzz_0_xxzz_0[i] * pb_y + g_0_zzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyy_0[i] = 3.0 * g_0_zzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyy_0[i] * pb_y + g_0_zzzzzz_0_xyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyz_0[i] = 2.0 * g_0_zzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyz_0[i] * pb_y + g_0_zzzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyzz_0[i] = g_0_zzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzz_0[i] * pb_y + g_0_zzzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xzzz_0[i] = g_0_zzzzzz_0_xzzz_0[i] * pb_y + g_0_zzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyy_0[i] = 4.0 * g_0_zzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyy_0[i] * pb_y + g_0_zzzzzz_0_yyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyz_0[i] = 3.0 * g_0_zzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyz_0[i] * pb_y + g_0_zzzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyzz_0[i] = 2.0 * g_0_zzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzz_0[i] * pb_y + g_0_zzzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yzzz_0[i] = g_0_zzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzz_0[i] * pb_y + g_0_zzzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_zzzz_0[i] = g_0_zzzzzz_0_zzzz_0[i] * pb_y + g_0_zzzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 525-540 components of targeted buffer : SKSG

    auto g_0_zzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 525);

    auto g_0_zzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 526);

    auto g_0_zzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 527);

    auto g_0_zzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 528);

    auto g_0_zzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 529);

    auto g_0_zzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 530);

    auto g_0_zzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 531);

    auto g_0_zzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 532);

    auto g_0_zzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 533);

    auto g_0_zzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 534);

    auto g_0_zzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 535);

    auto g_0_zzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 536);

    auto g_0_zzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 537);

    auto g_0_zzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 538);

    auto g_0_zzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 539);

#pragma omp simd aligned(g_0_zzzzz_0_xxxx_0,       \
                             g_0_zzzzz_0_xxxx_1,   \
                             g_0_zzzzz_0_xxxy_0,   \
                             g_0_zzzzz_0_xxxy_1,   \
                             g_0_zzzzz_0_xxxz_0,   \
                             g_0_zzzzz_0_xxxz_1,   \
                             g_0_zzzzz_0_xxyy_0,   \
                             g_0_zzzzz_0_xxyy_1,   \
                             g_0_zzzzz_0_xxyz_0,   \
                             g_0_zzzzz_0_xxyz_1,   \
                             g_0_zzzzz_0_xxzz_0,   \
                             g_0_zzzzz_0_xxzz_1,   \
                             g_0_zzzzz_0_xyyy_0,   \
                             g_0_zzzzz_0_xyyy_1,   \
                             g_0_zzzzz_0_xyyz_0,   \
                             g_0_zzzzz_0_xyyz_1,   \
                             g_0_zzzzz_0_xyzz_0,   \
                             g_0_zzzzz_0_xyzz_1,   \
                             g_0_zzzzz_0_xzzz_0,   \
                             g_0_zzzzz_0_xzzz_1,   \
                             g_0_zzzzz_0_yyyy_0,   \
                             g_0_zzzzz_0_yyyy_1,   \
                             g_0_zzzzz_0_yyyz_0,   \
                             g_0_zzzzz_0_yyyz_1,   \
                             g_0_zzzzz_0_yyzz_0,   \
                             g_0_zzzzz_0_yyzz_1,   \
                             g_0_zzzzz_0_yzzz_0,   \
                             g_0_zzzzz_0_yzzz_1,   \
                             g_0_zzzzz_0_zzzz_0,   \
                             g_0_zzzzz_0_zzzz_1,   \
                             g_0_zzzzzz_0_xxx_1,   \
                             g_0_zzzzzz_0_xxxx_0,  \
                             g_0_zzzzzz_0_xxxx_1,  \
                             g_0_zzzzzz_0_xxxy_0,  \
                             g_0_zzzzzz_0_xxxy_1,  \
                             g_0_zzzzzz_0_xxxz_0,  \
                             g_0_zzzzzz_0_xxxz_1,  \
                             g_0_zzzzzz_0_xxy_1,   \
                             g_0_zzzzzz_0_xxyy_0,  \
                             g_0_zzzzzz_0_xxyy_1,  \
                             g_0_zzzzzz_0_xxyz_0,  \
                             g_0_zzzzzz_0_xxyz_1,  \
                             g_0_zzzzzz_0_xxz_1,   \
                             g_0_zzzzzz_0_xxzz_0,  \
                             g_0_zzzzzz_0_xxzz_1,  \
                             g_0_zzzzzz_0_xyy_1,   \
                             g_0_zzzzzz_0_xyyy_0,  \
                             g_0_zzzzzz_0_xyyy_1,  \
                             g_0_zzzzzz_0_xyyz_0,  \
                             g_0_zzzzzz_0_xyyz_1,  \
                             g_0_zzzzzz_0_xyz_1,   \
                             g_0_zzzzzz_0_xyzz_0,  \
                             g_0_zzzzzz_0_xyzz_1,  \
                             g_0_zzzzzz_0_xzz_1,   \
                             g_0_zzzzzz_0_xzzz_0,  \
                             g_0_zzzzzz_0_xzzz_1,  \
                             g_0_zzzzzz_0_yyy_1,   \
                             g_0_zzzzzz_0_yyyy_0,  \
                             g_0_zzzzzz_0_yyyy_1,  \
                             g_0_zzzzzz_0_yyyz_0,  \
                             g_0_zzzzzz_0_yyyz_1,  \
                             g_0_zzzzzz_0_yyz_1,   \
                             g_0_zzzzzz_0_yyzz_0,  \
                             g_0_zzzzzz_0_yyzz_1,  \
                             g_0_zzzzzz_0_yzz_1,   \
                             g_0_zzzzzz_0_yzzz_0,  \
                             g_0_zzzzzz_0_yzzz_1,  \
                             g_0_zzzzzz_0_zzz_1,   \
                             g_0_zzzzzz_0_zzzz_0,  \
                             g_0_zzzzzz_0_zzzz_1,  \
                             g_0_zzzzzzz_0_xxxx_0, \
                             g_0_zzzzzzz_0_xxxy_0, \
                             g_0_zzzzzzz_0_xxxz_0, \
                             g_0_zzzzzzz_0_xxyy_0, \
                             g_0_zzzzzzz_0_xxyz_0, \
                             g_0_zzzzzzz_0_xxzz_0, \
                             g_0_zzzzzzz_0_xyyy_0, \
                             g_0_zzzzzzz_0_xyyz_0, \
                             g_0_zzzzzzz_0_xyzz_0, \
                             g_0_zzzzzzz_0_xzzz_0, \
                             g_0_zzzzzzz_0_yyyy_0, \
                             g_0_zzzzzzz_0_yyyz_0, \
                             g_0_zzzzzzz_0_yyzz_0, \
                             g_0_zzzzzzz_0_yzzz_0, \
                             g_0_zzzzzzz_0_zzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_xxxx_0[i] = 6.0 * g_0_zzzzz_0_xxxx_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxx_0[i] * pb_z +
                                  g_0_zzzzzz_0_xxxx_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxy_0[i] = 6.0 * g_0_zzzzz_0_xxxy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxy_0[i] * pb_z +
                                  g_0_zzzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxz_0[i] = 6.0 * g_0_zzzzz_0_xxxz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxx_1[i] * fi_abcd_0 +
                                  g_0_zzzzzz_0_xxxz_0[i] * pb_z + g_0_zzzzzz_0_xxxz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyy_0[i] = 6.0 * g_0_zzzzz_0_xxyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxyy_0[i] * pb_z +
                                  g_0_zzzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyz_0[i] = 6.0 * g_0_zzzzz_0_xxyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxy_1[i] * fi_abcd_0 +
                                  g_0_zzzzzz_0_xxyz_0[i] * pb_z + g_0_zzzzzz_0_xxyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxzz_0[i] = 6.0 * g_0_zzzzz_0_xxzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzz_0[i] * pb_z + g_0_zzzzzz_0_xxzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyy_0[i] = 6.0 * g_0_zzzzz_0_xyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyyy_0[i] * pb_z +
                                  g_0_zzzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyz_0[i] = 6.0 * g_0_zzzzz_0_xyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyy_1[i] * fi_abcd_0 +
                                  g_0_zzzzzz_0_xyyz_0[i] * pb_z + g_0_zzzzzz_0_xyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyzz_0[i] = 6.0 * g_0_zzzzz_0_xyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzz_0[i] * pb_z + g_0_zzzzzz_0_xyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xzzz_0[i] = 6.0 * g_0_zzzzz_0_xzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzz_0[i] * pb_z + g_0_zzzzzz_0_xzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyy_0[i] = 6.0 * g_0_zzzzz_0_yyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyyy_0[i] * pb_z +
                                  g_0_zzzzzz_0_yyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyz_0[i] = 6.0 * g_0_zzzzz_0_yyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyy_1[i] * fi_abcd_0 +
                                  g_0_zzzzzz_0_yyyz_0[i] * pb_z + g_0_zzzzzz_0_yyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyzz_0[i] = 6.0 * g_0_zzzzz_0_yyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzz_0[i] * pb_z + g_0_zzzzzz_0_yyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yzzz_0[i] = 6.0 * g_0_zzzzz_0_yzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzz_0[i] * pb_z + g_0_zzzzzz_0_yzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_zzzz_0[i] = 6.0 * g_0_zzzzz_0_zzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_zzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_zzzz_0[i] * pb_z + g_0_zzzzzz_0_zzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
