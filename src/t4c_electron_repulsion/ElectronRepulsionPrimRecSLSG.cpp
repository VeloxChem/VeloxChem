#include "ElectronRepulsionPrimRecSLSG.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slsg(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_slsg,
                                  size_t idx_eri_0_sisg,
                                  size_t idx_eri_1_sisg,
                                  size_t idx_eri_1_sksf,
                                  size_t idx_eri_0_sksg,
                                  size_t idx_eri_1_sksg,
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

    auto g_0_xxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 17);

    auto g_0_xxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 20);

    auto g_0_xxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 24);

    auto g_0_xxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 30);

    auto g_0_xxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 31);

    auto g_0_xxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 33);

    auto g_0_xxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 36);

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

    auto g_0_yyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 333);

    auto g_0_yyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 336);

    auto g_0_yyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 340);

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

    auto g_0_yzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 392);

    auto g_0_yzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 394);

    auto g_0_yzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 395);

    auto g_0_yzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 397);

    auto g_0_yzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 398);

    auto g_0_yzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 399);

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

    auto g_0_xxxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 17);

    auto g_0_xxxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 20);

    auto g_0_xxxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 24);

    auto g_0_xxxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sisg + 30);

    auto g_0_xxxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sisg + 31);

    auto g_0_xxxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 33);

    auto g_0_xxxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 36);

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

    auto g_0_yyyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sisg + 333);

    auto g_0_yyyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sisg + 336);

    auto g_0_yyyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sisg + 340);

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

    auto g_0_yzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sisg + 392);

    auto g_0_yzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sisg + 394);

    auto g_0_yzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sisg + 395);

    auto g_0_yzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sisg + 397);

    auto g_0_yzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sisg + 398);

    auto g_0_yzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sisg + 399);

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

    /// Set up components of auxilary buffer : SKSF

    auto g_0_xxxxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_sksf);

    auto g_0_xxxxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 1);

    auto g_0_xxxxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 2);

    auto g_0_xxxxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 3);

    auto g_0_xxxxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 4);

    auto g_0_xxxxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 5);

    auto g_0_xxxxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 6);

    auto g_0_xxxxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 7);

    auto g_0_xxxxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 8);

    auto g_0_xxxxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 9);

    auto g_0_xxxxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 22);

    auto g_0_xxxxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 24);

    auto g_0_xxxxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 25);

    auto g_0_xxxxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 27);

    auto g_0_xxxxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 28);

    auto g_0_xxxxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 29);

    auto g_0_xxxxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 30);

    auto g_0_xxxxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 31);

    auto g_0_xxxxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 32);

    auto g_0_xxxxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 33);

    auto g_0_xxxxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 34);

    auto g_0_xxxxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 35);

    auto g_0_xxxxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 36);

    auto g_0_xxxxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 37);

    auto g_0_xxxxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 38);

    auto g_0_xxxxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 39);

    auto g_0_xxxxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 50);

    auto g_0_xxxxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 51);

    auto g_0_xxxxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 52);

    auto g_0_xxxxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 53);

    auto g_0_xxxxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 54);

    auto g_0_xxxxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 55);

    auto g_0_xxxxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 56);

    auto g_0_xxxxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 57);

    auto g_0_xxxxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 58);

    auto g_0_xxxxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 59);

    auto g_0_xxxxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 60);

    auto g_0_xxxxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 61);

    auto g_0_xxxxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 62);

    auto g_0_xxxxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 63);

    auto g_0_xxxxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 64);

    auto g_0_xxxxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 65);

    auto g_0_xxxxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 66);

    auto g_0_xxxxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 67);

    auto g_0_xxxxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 68);

    auto g_0_xxxxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 69);

    auto g_0_xxxxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 90);

    auto g_0_xxxxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 91);

    auto g_0_xxxxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 92);

    auto g_0_xxxxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 93);

    auto g_0_xxxxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 94);

    auto g_0_xxxxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 95);

    auto g_0_xxxxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 96);

    auto g_0_xxxxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 97);

    auto g_0_xxxxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 98);

    auto g_0_xxxxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 99);

    auto g_0_xxxyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 100);

    auto g_0_xxxyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 101);

    auto g_0_xxxyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 102);

    auto g_0_xxxyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 103);

    auto g_0_xxxyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 104);

    auto g_0_xxxyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 105);

    auto g_0_xxxyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 106);

    auto g_0_xxxyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 107);

    auto g_0_xxxyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 108);

    auto g_0_xxxyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 109);

    auto g_0_xxxyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 124);

    auto g_0_xxxyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 127);

    auto g_0_xxxyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 128);

    auto g_0_xxxzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 140);

    auto g_0_xxxzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 141);

    auto g_0_xxxzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 142);

    auto g_0_xxxzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 143);

    auto g_0_xxxzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 144);

    auto g_0_xxxzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 145);

    auto g_0_xxxzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 146);

    auto g_0_xxxzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 147);

    auto g_0_xxxzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 148);

    auto g_0_xxxzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 149);

    auto g_0_xxyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 150);

    auto g_0_xxyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 151);

    auto g_0_xxyyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 152);

    auto g_0_xxyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 153);

    auto g_0_xxyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 154);

    auto g_0_xxyyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 155);

    auto g_0_xxyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 156);

    auto g_0_xxyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 157);

    auto g_0_xxyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 158);

    auto g_0_xxyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 159);

    auto g_0_xxyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 174);

    auto g_0_xxyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 177);

    auto g_0_xxyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 178);

    auto g_0_xxyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 184);

    auto g_0_xxyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 187);

    auto g_0_xxyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 188);

    auto g_0_xxzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 200);

    auto g_0_xxzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 201);

    auto g_0_xxzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 202);

    auto g_0_xxzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 203);

    auto g_0_xxzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 204);

    auto g_0_xxzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 205);

    auto g_0_xxzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 206);

    auto g_0_xxzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 207);

    auto g_0_xxzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 208);

    auto g_0_xxzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 209);

    auto g_0_xyyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 211);

    auto g_0_xyyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 213);

    auto g_0_xyyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 214);

    auto g_0_xyyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 216);

    auto g_0_xyyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 217);

    auto g_0_xyyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 218);

    auto g_0_xyyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 234);

    auto g_0_xyyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 237);

    auto g_0_xyyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 238);

    auto g_0_xyyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 244);

    auto g_0_xyyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 247);

    auto g_0_xyyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 248);

    auto g_0_xyyzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 254);

    auto g_0_xyyzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 257);

    auto g_0_xyyzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 258);

    auto g_0_xzzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 272);

    auto g_0_xzzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 274);

    auto g_0_xzzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 275);

    auto g_0_xzzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 277);

    auto g_0_xzzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 278);

    auto g_0_xzzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 279);

    auto g_0_yyyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 280);

    auto g_0_yyyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 281);

    auto g_0_yyyyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 282);

    auto g_0_yyyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 283);

    auto g_0_yyyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 284);

    auto g_0_yyyyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 285);

    auto g_0_yyyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 286);

    auto g_0_yyyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 287);

    auto g_0_yyyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 288);

    auto g_0_yyyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 289);

    auto g_0_yyyyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 292);

    auto g_0_yyyyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 294);

    auto g_0_yyyyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 295);

    auto g_0_yyyyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 297);

    auto g_0_yyyyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 298);

    auto g_0_yyyyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 299);

    auto g_0_yyyyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 300);

    auto g_0_yyyyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 301);

    auto g_0_yyyyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 302);

    auto g_0_yyyyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 303);

    auto g_0_yyyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 304);

    auto g_0_yyyyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 305);

    auto g_0_yyyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 306);

    auto g_0_yyyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 307);

    auto g_0_yyyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 308);

    auto g_0_yyyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 309);

    auto g_0_yyyyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 310);

    auto g_0_yyyyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 311);

    auto g_0_yyyyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 312);

    auto g_0_yyyyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 313);

    auto g_0_yyyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 314);

    auto g_0_yyyyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 315);

    auto g_0_yyyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 316);

    auto g_0_yyyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 317);

    auto g_0_yyyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 318);

    auto g_0_yyyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 319);

    auto g_0_yyyzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 320);

    auto g_0_yyyzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 321);

    auto g_0_yyyzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 322);

    auto g_0_yyyzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 323);

    auto g_0_yyyzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 324);

    auto g_0_yyyzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 325);

    auto g_0_yyyzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 326);

    auto g_0_yyyzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 327);

    auto g_0_yyyzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 328);

    auto g_0_yyyzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 329);

    auto g_0_yyzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 330);

    auto g_0_yyzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 331);

    auto g_0_yyzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 332);

    auto g_0_yyzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 333);

    auto g_0_yyzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 334);

    auto g_0_yyzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 335);

    auto g_0_yyzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 336);

    auto g_0_yyzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 337);

    auto g_0_yyzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 338);

    auto g_0_yyzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 339);

    auto g_0_yzzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 341);

    auto g_0_yzzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 342);

    auto g_0_yzzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 343);

    auto g_0_yzzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 344);

    auto g_0_yzzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 345);

    auto g_0_yzzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 346);

    auto g_0_yzzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 347);

    auto g_0_yzzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 348);

    auto g_0_yzzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 349);

    auto g_0_zzzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sksf + 350);

    auto g_0_zzzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sksf + 351);

    auto g_0_zzzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sksf + 352);

    auto g_0_zzzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sksf + 353);

    auto g_0_zzzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sksf + 354);

    auto g_0_zzzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sksf + 355);

    auto g_0_zzzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sksf + 356);

    auto g_0_zzzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sksf + 357);

    auto g_0_zzzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sksf + 358);

    auto g_0_zzzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sksf + 359);

    /// Set up components of auxilary buffer : SKSG

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

    auto g_0_xxxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 15);

    auto g_0_xxxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 16);

    auto g_0_xxxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 17);

    auto g_0_xxxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 18);

    auto g_0_xxxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 20);

    auto g_0_xxxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 21);

    auto g_0_xxxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 24);

    auto g_0_xxxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 25);

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

    auto g_0_xxxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 41);

    auto g_0_xxxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 42);

    auto g_0_xxxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 43);

    auto g_0_xxxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 44);

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

    auto g_0_xxxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 106);

    auto g_0_xxxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 108);

    auto g_0_xxxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 111);

    auto g_0_xxxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 120);

    auto g_0_xxxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 122);

    auto g_0_xxxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 125);

    auto g_0_xxxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 129);

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

    auto g_0_xxxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 166);

    auto g_0_xxxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 168);

    auto g_0_xxxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 171);

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

    auto g_0_xxxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 195);

    auto g_0_xxxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 197);

    auto g_0_xxxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 200);

    auto g_0_xxxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 204);

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

    auto g_0_xxyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 241);

    auto g_0_xxyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 243);

    auto g_0_xxyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 246);

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

    auto g_0_xxyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 285);

    auto g_0_xxyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 287);

    auto g_0_xxyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 290);

    auto g_0_xxyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 294);

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

    auto g_0_xyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 315);

    auto g_0_xyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sksg + 316);

    auto g_0_xyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sksg + 318);

    auto g_0_xyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 319);

    auto g_0_xyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sksg + 321);

    auto g_0_xyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 322);

    auto g_0_xyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 323);

    auto g_0_xyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 325);

    auto g_0_xyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 326);

    auto g_0_xyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 327);

    auto g_0_xyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 328);

    auto g_0_xyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 329);

    auto g_0_xyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 349);

    auto g_0_xyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 352);

    auto g_0_xyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 353);

    auto g_0_xyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 355);

    auto g_0_xyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 356);

    auto g_0_xyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 357);

    auto g_0_xyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 358);

    auto g_0_xyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 359);

    auto g_0_xyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 364);

    auto g_0_xyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 367);

    auto g_0_xyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 368);

    auto g_0_xyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 370);

    auto g_0_xyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 371);

    auto g_0_xyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 372);

    auto g_0_xyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 373);

    auto g_0_xyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 374);

    auto g_0_xyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 379);

    auto g_0_xyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 382);

    auto g_0_xyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 383);

    auto g_0_xyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 385);

    auto g_0_xyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 386);

    auto g_0_xyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 387);

    auto g_0_xyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 388);

    auto g_0_xyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 389);

    auto g_0_xzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sksg + 405);

    auto g_0_xzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sksg + 407);

    auto g_0_xzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sksg + 409);

    auto g_0_xzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sksg + 410);

    auto g_0_xzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sksg + 412);

    auto g_0_xzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sksg + 413);

    auto g_0_xzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sksg + 414);

    auto g_0_xzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sksg + 415);

    auto g_0_xzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sksg + 416);

    auto g_0_xzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sksg + 417);

    auto g_0_xzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sksg + 418);

    auto g_0_xzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sksg + 419);

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

    /// Set up components of auxilary buffer : SKSG

    auto g_0_xxxxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg);

    auto g_0_xxxxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 1);

    auto g_0_xxxxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 2);

    auto g_0_xxxxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 3);

    auto g_0_xxxxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 4);

    auto g_0_xxxxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 5);

    auto g_0_xxxxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 6);

    auto g_0_xxxxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 7);

    auto g_0_xxxxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 8);

    auto g_0_xxxxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 9);

    auto g_0_xxxxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 10);

    auto g_0_xxxxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 11);

    auto g_0_xxxxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 12);

    auto g_0_xxxxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 13);

    auto g_0_xxxxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 14);

    auto g_0_xxxxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 15);

    auto g_0_xxxxxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 16);

    auto g_0_xxxxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 17);

    auto g_0_xxxxxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 18);

    auto g_0_xxxxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 20);

    auto g_0_xxxxxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 21);

    auto g_0_xxxxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 24);

    auto g_0_xxxxxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 25);

    auto g_0_xxxxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 30);

    auto g_0_xxxxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 31);

    auto g_0_xxxxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 32);

    auto g_0_xxxxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 33);

    auto g_0_xxxxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 34);

    auto g_0_xxxxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 35);

    auto g_0_xxxxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 36);

    auto g_0_xxxxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 37);

    auto g_0_xxxxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 38);

    auto g_0_xxxxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 39);

    auto g_0_xxxxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 41);

    auto g_0_xxxxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 42);

    auto g_0_xxxxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 43);

    auto g_0_xxxxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 44);

    auto g_0_xxxxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 45);

    auto g_0_xxxxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 46);

    auto g_0_xxxxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 47);

    auto g_0_xxxxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 48);

    auto g_0_xxxxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 49);

    auto g_0_xxxxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 50);

    auto g_0_xxxxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 51);

    auto g_0_xxxxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 52);

    auto g_0_xxxxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 53);

    auto g_0_xxxxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 54);

    auto g_0_xxxxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 55);

    auto g_0_xxxxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 56);

    auto g_0_xxxxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 57);

    auto g_0_xxxxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 58);

    auto g_0_xxxxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 59);

    auto g_0_xxxxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 75);

    auto g_0_xxxxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 76);

    auto g_0_xxxxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 77);

    auto g_0_xxxxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 78);

    auto g_0_xxxxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 79);

    auto g_0_xxxxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 80);

    auto g_0_xxxxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 81);

    auto g_0_xxxxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 82);

    auto g_0_xxxxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 83);

    auto g_0_xxxxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 84);

    auto g_0_xxxxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 85);

    auto g_0_xxxxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 86);

    auto g_0_xxxxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 87);

    auto g_0_xxxxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 88);

    auto g_0_xxxxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 89);

    auto g_0_xxxxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 90);

    auto g_0_xxxxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 91);

    auto g_0_xxxxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 92);

    auto g_0_xxxxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 93);

    auto g_0_xxxxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 94);

    auto g_0_xxxxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 95);

    auto g_0_xxxxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 96);

    auto g_0_xxxxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 97);

    auto g_0_xxxxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 98);

    auto g_0_xxxxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 99);

    auto g_0_xxxxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 100);

    auto g_0_xxxxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 101);

    auto g_0_xxxxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 102);

    auto g_0_xxxxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 103);

    auto g_0_xxxxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 104);

    auto g_0_xxxxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 106);

    auto g_0_xxxxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 108);

    auto g_0_xxxxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 111);

    auto g_0_xxxxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 120);

    auto g_0_xxxxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 122);

    auto g_0_xxxxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 125);

    auto g_0_xxxxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 129);

    auto g_0_xxxxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 135);

    auto g_0_xxxxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 136);

    auto g_0_xxxxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 137);

    auto g_0_xxxxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 138);

    auto g_0_xxxxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 139);

    auto g_0_xxxxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 140);

    auto g_0_xxxxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 141);

    auto g_0_xxxxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 142);

    auto g_0_xxxxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 143);

    auto g_0_xxxxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 144);

    auto g_0_xxxxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 145);

    auto g_0_xxxxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 146);

    auto g_0_xxxxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 147);

    auto g_0_xxxxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 148);

    auto g_0_xxxxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 149);

    auto g_0_xxxyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 150);

    auto g_0_xxxyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 151);

    auto g_0_xxxyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 152);

    auto g_0_xxxyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 153);

    auto g_0_xxxyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 154);

    auto g_0_xxxyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 155);

    auto g_0_xxxyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 156);

    auto g_0_xxxyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 157);

    auto g_0_xxxyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 158);

    auto g_0_xxxyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 159);

    auto g_0_xxxyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 160);

    auto g_0_xxxyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 161);

    auto g_0_xxxyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 162);

    auto g_0_xxxyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 163);

    auto g_0_xxxyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 164);

    auto g_0_xxxyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 166);

    auto g_0_xxxyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 168);

    auto g_0_xxxyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 171);

    auto g_0_xxxyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 180);

    auto g_0_xxxyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 181);

    auto g_0_xxxyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 182);

    auto g_0_xxxyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 183);

    auto g_0_xxxyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 184);

    auto g_0_xxxyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 185);

    auto g_0_xxxyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 186);

    auto g_0_xxxyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 187);

    auto g_0_xxxyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 188);

    auto g_0_xxxyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 189);

    auto g_0_xxxyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 190);

    auto g_0_xxxyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 191);

    auto g_0_xxxyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 192);

    auto g_0_xxxyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 193);

    auto g_0_xxxyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 194);

    auto g_0_xxxyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 195);

    auto g_0_xxxyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 197);

    auto g_0_xxxyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 200);

    auto g_0_xxxyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 204);

    auto g_0_xxxzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 210);

    auto g_0_xxxzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 211);

    auto g_0_xxxzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 212);

    auto g_0_xxxzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 213);

    auto g_0_xxxzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 214);

    auto g_0_xxxzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 215);

    auto g_0_xxxzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 216);

    auto g_0_xxxzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 217);

    auto g_0_xxxzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 218);

    auto g_0_xxxzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 219);

    auto g_0_xxxzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 220);

    auto g_0_xxxzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 221);

    auto g_0_xxxzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 222);

    auto g_0_xxxzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 223);

    auto g_0_xxxzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 224);

    auto g_0_xxyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 225);

    auto g_0_xxyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 226);

    auto g_0_xxyyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 227);

    auto g_0_xxyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 228);

    auto g_0_xxyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 229);

    auto g_0_xxyyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 230);

    auto g_0_xxyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 231);

    auto g_0_xxyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 232);

    auto g_0_xxyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 233);

    auto g_0_xxyyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 234);

    auto g_0_xxyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 235);

    auto g_0_xxyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 236);

    auto g_0_xxyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 237);

    auto g_0_xxyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 238);

    auto g_0_xxyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 239);

    auto g_0_xxyyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 241);

    auto g_0_xxyyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 243);

    auto g_0_xxyyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 246);

    auto g_0_xxyyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 255);

    auto g_0_xxyyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 256);

    auto g_0_xxyyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 257);

    auto g_0_xxyyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 258);

    auto g_0_xxyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 259);

    auto g_0_xxyyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 260);

    auto g_0_xxyyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 261);

    auto g_0_xxyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 262);

    auto g_0_xxyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 263);

    auto g_0_xxyyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 264);

    auto g_0_xxyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 265);

    auto g_0_xxyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 266);

    auto g_0_xxyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 267);

    auto g_0_xxyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 268);

    auto g_0_xxyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 269);

    auto g_0_xxyyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 270);

    auto g_0_xxyyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 271);

    auto g_0_xxyyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 272);

    auto g_0_xxyyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 273);

    auto g_0_xxyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 274);

    auto g_0_xxyyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 275);

    auto g_0_xxyyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 276);

    auto g_0_xxyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 277);

    auto g_0_xxyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 278);

    auto g_0_xxyyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 279);

    auto g_0_xxyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 280);

    auto g_0_xxyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 281);

    auto g_0_xxyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 282);

    auto g_0_xxyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 283);

    auto g_0_xxyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 284);

    auto g_0_xxyzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 285);

    auto g_0_xxyzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 287);

    auto g_0_xxyzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 290);

    auto g_0_xxyzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 294);

    auto g_0_xxzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 300);

    auto g_0_xxzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 301);

    auto g_0_xxzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 302);

    auto g_0_xxzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 303);

    auto g_0_xxzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 304);

    auto g_0_xxzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 305);

    auto g_0_xxzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 306);

    auto g_0_xxzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 307);

    auto g_0_xxzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 308);

    auto g_0_xxzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 309);

    auto g_0_xxzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 310);

    auto g_0_xxzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 311);

    auto g_0_xxzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 312);

    auto g_0_xxzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 313);

    auto g_0_xxzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 314);

    auto g_0_xyyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 315);

    auto g_0_xyyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 316);

    auto g_0_xyyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 318);

    auto g_0_xyyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 319);

    auto g_0_xyyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 321);

    auto g_0_xyyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 322);

    auto g_0_xyyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 323);

    auto g_0_xyyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 325);

    auto g_0_xyyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 326);

    auto g_0_xyyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 327);

    auto g_0_xyyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 328);

    auto g_0_xyyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 329);

    auto g_0_xyyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 349);

    auto g_0_xyyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 352);

    auto g_0_xyyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 353);

    auto g_0_xyyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 355);

    auto g_0_xyyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 356);

    auto g_0_xyyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 357);

    auto g_0_xyyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 358);

    auto g_0_xyyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 359);

    auto g_0_xyyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 364);

    auto g_0_xyyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 367);

    auto g_0_xyyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 368);

    auto g_0_xyyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 370);

    auto g_0_xyyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 371);

    auto g_0_xyyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 372);

    auto g_0_xyyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 373);

    auto g_0_xyyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 374);

    auto g_0_xyyzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 379);

    auto g_0_xyyzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 382);

    auto g_0_xyyzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 383);

    auto g_0_xyyzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 385);

    auto g_0_xyyzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 386);

    auto g_0_xyyzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 387);

    auto g_0_xyyzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 388);

    auto g_0_xyyzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 389);

    auto g_0_xzzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 405);

    auto g_0_xzzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 407);

    auto g_0_xzzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 409);

    auto g_0_xzzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 410);

    auto g_0_xzzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 412);

    auto g_0_xzzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 413);

    auto g_0_xzzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 414);

    auto g_0_xzzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 415);

    auto g_0_xzzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 416);

    auto g_0_xzzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 417);

    auto g_0_xzzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 418);

    auto g_0_xzzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 419);

    auto g_0_yyyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 420);

    auto g_0_yyyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 421);

    auto g_0_yyyyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 422);

    auto g_0_yyyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 423);

    auto g_0_yyyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 424);

    auto g_0_yyyyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 425);

    auto g_0_yyyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 426);

    auto g_0_yyyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 427);

    auto g_0_yyyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 428);

    auto g_0_yyyyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 429);

    auto g_0_yyyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 430);

    auto g_0_yyyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 431);

    auto g_0_yyyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 432);

    auto g_0_yyyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 433);

    auto g_0_yyyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 434);

    auto g_0_yyyyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 436);

    auto g_0_yyyyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 437);

    auto g_0_yyyyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 438);

    auto g_0_yyyyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 439);

    auto g_0_yyyyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 440);

    auto g_0_yyyyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 441);

    auto g_0_yyyyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 442);

    auto g_0_yyyyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 443);

    auto g_0_yyyyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 444);

    auto g_0_yyyyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 445);

    auto g_0_yyyyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 446);

    auto g_0_yyyyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 447);

    auto g_0_yyyyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 448);

    auto g_0_yyyyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 449);

    auto g_0_yyyyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 450);

    auto g_0_yyyyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 451);

    auto g_0_yyyyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 452);

    auto g_0_yyyyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 453);

    auto g_0_yyyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 454);

    auto g_0_yyyyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 455);

    auto g_0_yyyyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 456);

    auto g_0_yyyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 457);

    auto g_0_yyyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 458);

    auto g_0_yyyyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 459);

    auto g_0_yyyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 460);

    auto g_0_yyyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 461);

    auto g_0_yyyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 462);

    auto g_0_yyyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 463);

    auto g_0_yyyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 464);

    auto g_0_yyyyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 465);

    auto g_0_yyyyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 466);

    auto g_0_yyyyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 467);

    auto g_0_yyyyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 468);

    auto g_0_yyyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 469);

    auto g_0_yyyyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 470);

    auto g_0_yyyyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 471);

    auto g_0_yyyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 472);

    auto g_0_yyyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 473);

    auto g_0_yyyyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 474);

    auto g_0_yyyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 475);

    auto g_0_yyyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 476);

    auto g_0_yyyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 477);

    auto g_0_yyyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 478);

    auto g_0_yyyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 479);

    auto g_0_yyyzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 480);

    auto g_0_yyyzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 481);

    auto g_0_yyyzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 482);

    auto g_0_yyyzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 483);

    auto g_0_yyyzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 484);

    auto g_0_yyyzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 485);

    auto g_0_yyyzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 486);

    auto g_0_yyyzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 487);

    auto g_0_yyyzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 488);

    auto g_0_yyyzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 489);

    auto g_0_yyyzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 490);

    auto g_0_yyyzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 491);

    auto g_0_yyyzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 492);

    auto g_0_yyyzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 493);

    auto g_0_yyyzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 494);

    auto g_0_yyzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 495);

    auto g_0_yyzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 496);

    auto g_0_yyzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 497);

    auto g_0_yyzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 498);

    auto g_0_yyzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 499);

    auto g_0_yyzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 500);

    auto g_0_yyzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 501);

    auto g_0_yyzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 502);

    auto g_0_yyzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 503);

    auto g_0_yyzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 504);

    auto g_0_yyzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 505);

    auto g_0_yyzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 506);

    auto g_0_yyzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 507);

    auto g_0_yyzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 508);

    auto g_0_yyzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 509);

    auto g_0_yzzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 510);

    auto g_0_yzzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 511);

    auto g_0_yzzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 512);

    auto g_0_yzzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 513);

    auto g_0_yzzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 514);

    auto g_0_yzzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 515);

    auto g_0_yzzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 516);

    auto g_0_yzzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 517);

    auto g_0_yzzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 518);

    auto g_0_yzzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 519);

    auto g_0_yzzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 520);

    auto g_0_yzzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 521);

    auto g_0_yzzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 522);

    auto g_0_yzzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 523);

    auto g_0_yzzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 524);

    auto g_0_zzzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sksg + 525);

    auto g_0_zzzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sksg + 526);

    auto g_0_zzzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sksg + 527);

    auto g_0_zzzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sksg + 528);

    auto g_0_zzzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sksg + 529);

    auto g_0_zzzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sksg + 530);

    auto g_0_zzzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sksg + 531);

    auto g_0_zzzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sksg + 532);

    auto g_0_zzzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sksg + 533);

    auto g_0_zzzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sksg + 534);

    auto g_0_zzzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sksg + 535);

    auto g_0_zzzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sksg + 536);

    auto g_0_zzzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sksg + 537);

    auto g_0_zzzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sksg + 538);

    auto g_0_zzzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sksg + 539);

    /// Set up 0-15 components of targeted buffer : SLSG

    auto g_0_xxxxxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg);

    auto g_0_xxxxxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 1);

    auto g_0_xxxxxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 2);

    auto g_0_xxxxxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 3);

    auto g_0_xxxxxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 4);

    auto g_0_xxxxxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 5);

    auto g_0_xxxxxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 6);

    auto g_0_xxxxxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 7);

    auto g_0_xxxxxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 8);

    auto g_0_xxxxxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 9);

    auto g_0_xxxxxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 10);

    auto g_0_xxxxxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 11);

    auto g_0_xxxxxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 12);

    auto g_0_xxxxxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 13);

    auto g_0_xxxxxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 14);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxx_0, g_0_xxxxxx_0_xxxx_1, g_0_xxxxxx_0_xxxy_0, g_0_xxxxxx_0_xxxy_1, g_0_xxxxxx_0_xxxz_0, g_0_xxxxxx_0_xxxz_1, g_0_xxxxxx_0_xxyy_0, g_0_xxxxxx_0_xxyy_1, g_0_xxxxxx_0_xxyz_0, g_0_xxxxxx_0_xxyz_1, g_0_xxxxxx_0_xxzz_0, g_0_xxxxxx_0_xxzz_1, g_0_xxxxxx_0_xyyy_0, g_0_xxxxxx_0_xyyy_1, g_0_xxxxxx_0_xyyz_0, g_0_xxxxxx_0_xyyz_1, g_0_xxxxxx_0_xyzz_0, g_0_xxxxxx_0_xyzz_1, g_0_xxxxxx_0_xzzz_0, g_0_xxxxxx_0_xzzz_1, g_0_xxxxxx_0_yyyy_0, g_0_xxxxxx_0_yyyy_1, g_0_xxxxxx_0_yyyz_0, g_0_xxxxxx_0_yyyz_1, g_0_xxxxxx_0_yyzz_0, g_0_xxxxxx_0_yyzz_1, g_0_xxxxxx_0_yzzz_0, g_0_xxxxxx_0_yzzz_1, g_0_xxxxxx_0_zzzz_0, g_0_xxxxxx_0_zzzz_1, g_0_xxxxxxx_0_xxx_1, g_0_xxxxxxx_0_xxxx_0, g_0_xxxxxxx_0_xxxx_1, g_0_xxxxxxx_0_xxxy_0, g_0_xxxxxxx_0_xxxy_1, g_0_xxxxxxx_0_xxxz_0, g_0_xxxxxxx_0_xxxz_1, g_0_xxxxxxx_0_xxy_1, g_0_xxxxxxx_0_xxyy_0, g_0_xxxxxxx_0_xxyy_1, g_0_xxxxxxx_0_xxyz_0, g_0_xxxxxxx_0_xxyz_1, g_0_xxxxxxx_0_xxz_1, g_0_xxxxxxx_0_xxzz_0, g_0_xxxxxxx_0_xxzz_1, g_0_xxxxxxx_0_xyy_1, g_0_xxxxxxx_0_xyyy_0, g_0_xxxxxxx_0_xyyy_1, g_0_xxxxxxx_0_xyyz_0, g_0_xxxxxxx_0_xyyz_1, g_0_xxxxxxx_0_xyz_1, g_0_xxxxxxx_0_xyzz_0, g_0_xxxxxxx_0_xyzz_1, g_0_xxxxxxx_0_xzz_1, g_0_xxxxxxx_0_xzzz_0, g_0_xxxxxxx_0_xzzz_1, g_0_xxxxxxx_0_yyy_1, g_0_xxxxxxx_0_yyyy_0, g_0_xxxxxxx_0_yyyy_1, g_0_xxxxxxx_0_yyyz_0, g_0_xxxxxxx_0_yyyz_1, g_0_xxxxxxx_0_yyz_1, g_0_xxxxxxx_0_yyzz_0, g_0_xxxxxxx_0_yyzz_1, g_0_xxxxxxx_0_yzz_1, g_0_xxxxxxx_0_yzzz_0, g_0_xxxxxxx_0_yzzz_1, g_0_xxxxxxx_0_zzz_1, g_0_xxxxxxx_0_zzzz_0, g_0_xxxxxxx_0_zzzz_1, g_0_xxxxxxxx_0_xxxx_0, g_0_xxxxxxxx_0_xxxy_0, g_0_xxxxxxxx_0_xxxz_0, g_0_xxxxxxxx_0_xxyy_0, g_0_xxxxxxxx_0_xxyz_0, g_0_xxxxxxxx_0_xxzz_0, g_0_xxxxxxxx_0_xyyy_0, g_0_xxxxxxxx_0_xyyz_0, g_0_xxxxxxxx_0_xyzz_0, g_0_xxxxxxxx_0_xzzz_0, g_0_xxxxxxxx_0_yyyy_0, g_0_xxxxxxxx_0_yyyz_0, g_0_xxxxxxxx_0_yyzz_0, g_0_xxxxxxxx_0_yzzz_0, g_0_xxxxxxxx_0_zzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_xxxx_0[i] = 7.0 * g_0_xxxxxx_0_xxxx_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxx_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxx_0[i] * pb_x + g_0_xxxxxxx_0_xxxx_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxy_0[i] = 7.0 * g_0_xxxxxx_0_xxxy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxy_0[i] * pb_x + g_0_xxxxxxx_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxz_0[i] = 7.0 * g_0_xxxxxx_0_xxxz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxz_0[i] * pb_x + g_0_xxxxxxx_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyy_0[i] = 7.0 * g_0_xxxxxx_0_xxyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyy_0[i] * pb_x + g_0_xxxxxxx_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyz_0[i] = 7.0 * g_0_xxxxxx_0_xxyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyz_0[i] * pb_x + g_0_xxxxxxx_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxzz_0[i] = 7.0 * g_0_xxxxxx_0_xxzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzz_0[i] * pb_x + g_0_xxxxxxx_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyy_0[i] = 7.0 * g_0_xxxxxx_0_xyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyy_0[i] * pb_x + g_0_xxxxxxx_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyz_0[i] = 7.0 * g_0_xxxxxx_0_xyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyz_0[i] * pb_x + g_0_xxxxxxx_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyzz_0[i] = 7.0 * g_0_xxxxxx_0_xyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzz_0[i] * pb_x + g_0_xxxxxxx_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xzzz_0[i] = 7.0 * g_0_xxxxxx_0_xzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzz_0[i] * pb_x + g_0_xxxxxxx_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyy_0[i] = 7.0 * g_0_xxxxxx_0_yyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyy_0[i] * pb_x + g_0_xxxxxxx_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyz_0[i] = 7.0 * g_0_xxxxxx_0_yyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyz_0[i] * pb_x + g_0_xxxxxxx_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyzz_0[i] = 7.0 * g_0_xxxxxx_0_yyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyzz_0[i] * pb_x + g_0_xxxxxxx_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yzzz_0[i] = 7.0 * g_0_xxxxxx_0_yzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzzz_0[i] * pb_x + g_0_xxxxxxx_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_zzzz_0[i] = 7.0 * g_0_xxxxxx_0_zzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzzz_0[i] * pb_x + g_0_xxxxxxx_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : SLSG

    auto g_0_xxxxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 15);

    auto g_0_xxxxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 16);

    auto g_0_xxxxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 17);

    auto g_0_xxxxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 18);

    auto g_0_xxxxxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 19);

    auto g_0_xxxxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 20);

    auto g_0_xxxxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 21);

    auto g_0_xxxxxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 22);

    auto g_0_xxxxxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 23);

    auto g_0_xxxxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 24);

    auto g_0_xxxxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 25);

    auto g_0_xxxxxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 26);

    auto g_0_xxxxxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 27);

    auto g_0_xxxxxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 28);

    auto g_0_xxxxxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 29);

    #pragma omp simd aligned(g_0_xxxxxxx_0_xxx_1, g_0_xxxxxxx_0_xxxx_0, g_0_xxxxxxx_0_xxxx_1, g_0_xxxxxxx_0_xxxy_0, g_0_xxxxxxx_0_xxxy_1, g_0_xxxxxxx_0_xxxz_0, g_0_xxxxxxx_0_xxxz_1, g_0_xxxxxxx_0_xxy_1, g_0_xxxxxxx_0_xxyy_0, g_0_xxxxxxx_0_xxyy_1, g_0_xxxxxxx_0_xxyz_0, g_0_xxxxxxx_0_xxyz_1, g_0_xxxxxxx_0_xxz_1, g_0_xxxxxxx_0_xxzz_0, g_0_xxxxxxx_0_xxzz_1, g_0_xxxxxxx_0_xyy_1, g_0_xxxxxxx_0_xyyy_0, g_0_xxxxxxx_0_xyyy_1, g_0_xxxxxxx_0_xyyz_0, g_0_xxxxxxx_0_xyyz_1, g_0_xxxxxxx_0_xyz_1, g_0_xxxxxxx_0_xyzz_0, g_0_xxxxxxx_0_xyzz_1, g_0_xxxxxxx_0_xzz_1, g_0_xxxxxxx_0_xzzz_0, g_0_xxxxxxx_0_xzzz_1, g_0_xxxxxxx_0_yyy_1, g_0_xxxxxxx_0_yyyy_0, g_0_xxxxxxx_0_yyyy_1, g_0_xxxxxxx_0_yyyz_0, g_0_xxxxxxx_0_yyyz_1, g_0_xxxxxxx_0_yyz_1, g_0_xxxxxxx_0_yyzz_0, g_0_xxxxxxx_0_yyzz_1, g_0_xxxxxxx_0_yzz_1, g_0_xxxxxxx_0_yzzz_0, g_0_xxxxxxx_0_yzzz_1, g_0_xxxxxxx_0_zzz_1, g_0_xxxxxxx_0_zzzz_0, g_0_xxxxxxx_0_zzzz_1, g_0_xxxxxxxy_0_xxxx_0, g_0_xxxxxxxy_0_xxxy_0, g_0_xxxxxxxy_0_xxxz_0, g_0_xxxxxxxy_0_xxyy_0, g_0_xxxxxxxy_0_xxyz_0, g_0_xxxxxxxy_0_xxzz_0, g_0_xxxxxxxy_0_xyyy_0, g_0_xxxxxxxy_0_xyyz_0, g_0_xxxxxxxy_0_xyzz_0, g_0_xxxxxxxy_0_xzzz_0, g_0_xxxxxxxy_0_yyyy_0, g_0_xxxxxxxy_0_yyyz_0, g_0_xxxxxxxy_0_yyzz_0, g_0_xxxxxxxy_0_yzzz_0, g_0_xxxxxxxy_0_zzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_xxxx_0[i] = g_0_xxxxxxx_0_xxxx_0[i] * pb_y + g_0_xxxxxxx_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxy_0[i] = g_0_xxxxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxy_0[i] * pb_y + g_0_xxxxxxx_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxz_0[i] = g_0_xxxxxxx_0_xxxz_0[i] * pb_y + g_0_xxxxxxx_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyy_0[i] = 2.0 * g_0_xxxxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyy_0[i] * pb_y + g_0_xxxxxxx_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyz_0[i] = g_0_xxxxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyz_0[i] * pb_y + g_0_xxxxxxx_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxzz_0[i] = g_0_xxxxxxx_0_xxzz_0[i] * pb_y + g_0_xxxxxxx_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyy_0[i] = 3.0 * g_0_xxxxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyy_0[i] * pb_y + g_0_xxxxxxx_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyz_0[i] = 2.0 * g_0_xxxxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyz_0[i] * pb_y + g_0_xxxxxxx_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyzz_0[i] = g_0_xxxxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzz_0[i] * pb_y + g_0_xxxxxxx_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xzzz_0[i] = g_0_xxxxxxx_0_xzzz_0[i] * pb_y + g_0_xxxxxxx_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyy_0[i] = 4.0 * g_0_xxxxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyy_0[i] * pb_y + g_0_xxxxxxx_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyz_0[i] = 3.0 * g_0_xxxxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyz_0[i] * pb_y + g_0_xxxxxxx_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyzz_0[i] = 2.0 * g_0_xxxxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzz_0[i] * pb_y + g_0_xxxxxxx_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yzzz_0[i] = g_0_xxxxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzz_0[i] * pb_y + g_0_xxxxxxx_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_zzzz_0[i] = g_0_xxxxxxx_0_zzzz_0[i] * pb_y + g_0_xxxxxxx_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : SLSG

    auto g_0_xxxxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 30);

    auto g_0_xxxxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 31);

    auto g_0_xxxxxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 32);

    auto g_0_xxxxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 33);

    auto g_0_xxxxxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 34);

    auto g_0_xxxxxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 35);

    auto g_0_xxxxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 36);

    auto g_0_xxxxxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 37);

    auto g_0_xxxxxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 38);

    auto g_0_xxxxxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 39);

    auto g_0_xxxxxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 40);

    auto g_0_xxxxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 41);

    auto g_0_xxxxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 42);

    auto g_0_xxxxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 43);

    auto g_0_xxxxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 44);

    #pragma omp simd aligned(g_0_xxxxxxx_0_xxx_1, g_0_xxxxxxx_0_xxxx_0, g_0_xxxxxxx_0_xxxx_1, g_0_xxxxxxx_0_xxxy_0, g_0_xxxxxxx_0_xxxy_1, g_0_xxxxxxx_0_xxxz_0, g_0_xxxxxxx_0_xxxz_1, g_0_xxxxxxx_0_xxy_1, g_0_xxxxxxx_0_xxyy_0, g_0_xxxxxxx_0_xxyy_1, g_0_xxxxxxx_0_xxyz_0, g_0_xxxxxxx_0_xxyz_1, g_0_xxxxxxx_0_xxz_1, g_0_xxxxxxx_0_xxzz_0, g_0_xxxxxxx_0_xxzz_1, g_0_xxxxxxx_0_xyy_1, g_0_xxxxxxx_0_xyyy_0, g_0_xxxxxxx_0_xyyy_1, g_0_xxxxxxx_0_xyyz_0, g_0_xxxxxxx_0_xyyz_1, g_0_xxxxxxx_0_xyz_1, g_0_xxxxxxx_0_xyzz_0, g_0_xxxxxxx_0_xyzz_1, g_0_xxxxxxx_0_xzz_1, g_0_xxxxxxx_0_xzzz_0, g_0_xxxxxxx_0_xzzz_1, g_0_xxxxxxx_0_yyy_1, g_0_xxxxxxx_0_yyyy_0, g_0_xxxxxxx_0_yyyy_1, g_0_xxxxxxx_0_yyyz_0, g_0_xxxxxxx_0_yyyz_1, g_0_xxxxxxx_0_yyz_1, g_0_xxxxxxx_0_yyzz_0, g_0_xxxxxxx_0_yyzz_1, g_0_xxxxxxx_0_yzz_1, g_0_xxxxxxx_0_yzzz_0, g_0_xxxxxxx_0_yzzz_1, g_0_xxxxxxx_0_zzz_1, g_0_xxxxxxx_0_zzzz_0, g_0_xxxxxxx_0_zzzz_1, g_0_xxxxxxxz_0_xxxx_0, g_0_xxxxxxxz_0_xxxy_0, g_0_xxxxxxxz_0_xxxz_0, g_0_xxxxxxxz_0_xxyy_0, g_0_xxxxxxxz_0_xxyz_0, g_0_xxxxxxxz_0_xxzz_0, g_0_xxxxxxxz_0_xyyy_0, g_0_xxxxxxxz_0_xyyz_0, g_0_xxxxxxxz_0_xyzz_0, g_0_xxxxxxxz_0_xzzz_0, g_0_xxxxxxxz_0_yyyy_0, g_0_xxxxxxxz_0_yyyz_0, g_0_xxxxxxxz_0_yyzz_0, g_0_xxxxxxxz_0_yzzz_0, g_0_xxxxxxxz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_xxxx_0[i] = g_0_xxxxxxx_0_xxxx_0[i] * pb_z + g_0_xxxxxxx_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxy_0[i] = g_0_xxxxxxx_0_xxxy_0[i] * pb_z + g_0_xxxxxxx_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxz_0[i] = g_0_xxxxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxz_0[i] * pb_z + g_0_xxxxxxx_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyy_0[i] = g_0_xxxxxxx_0_xxyy_0[i] * pb_z + g_0_xxxxxxx_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyz_0[i] = g_0_xxxxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyz_0[i] * pb_z + g_0_xxxxxxx_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzz_0[i] * pb_z + g_0_xxxxxxx_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyy_0[i] = g_0_xxxxxxx_0_xyyy_0[i] * pb_z + g_0_xxxxxxx_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyz_0[i] = g_0_xxxxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyz_0[i] * pb_z + g_0_xxxxxxx_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzz_0[i] * pb_z + g_0_xxxxxxx_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzz_0[i] * pb_z + g_0_xxxxxxx_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyy_0[i] = g_0_xxxxxxx_0_yyyy_0[i] * pb_z + g_0_xxxxxxx_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyz_0[i] = g_0_xxxxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyz_0[i] * pb_z + g_0_xxxxxxx_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyzz_0[i] = 2.0 * g_0_xxxxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzz_0[i] * pb_z + g_0_xxxxxxx_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yzzz_0[i] = 3.0 * g_0_xxxxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzz_0[i] * pb_z + g_0_xxxxxxx_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_zzzz_0[i] = 4.0 * g_0_xxxxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_zzzz_0[i] * pb_z + g_0_xxxxxxx_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 45-60 components of targeted buffer : SLSG

    auto g_0_xxxxxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 45);

    auto g_0_xxxxxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 46);

    auto g_0_xxxxxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 47);

    auto g_0_xxxxxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 48);

    auto g_0_xxxxxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 49);

    auto g_0_xxxxxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 50);

    auto g_0_xxxxxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 51);

    auto g_0_xxxxxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 52);

    auto g_0_xxxxxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 53);

    auto g_0_xxxxxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 54);

    auto g_0_xxxxxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 55);

    auto g_0_xxxxxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 56);

    auto g_0_xxxxxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 57);

    auto g_0_xxxxxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 58);

    auto g_0_xxxxxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 59);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxx_0, g_0_xxxxxx_0_xxxx_1, g_0_xxxxxx_0_xxxz_0, g_0_xxxxxx_0_xxxz_1, g_0_xxxxxx_0_xxzz_0, g_0_xxxxxx_0_xxzz_1, g_0_xxxxxx_0_xzzz_0, g_0_xxxxxx_0_xzzz_1, g_0_xxxxxxy_0_xxxx_0, g_0_xxxxxxy_0_xxxx_1, g_0_xxxxxxy_0_xxxz_0, g_0_xxxxxxy_0_xxxz_1, g_0_xxxxxxy_0_xxzz_0, g_0_xxxxxxy_0_xxzz_1, g_0_xxxxxxy_0_xzzz_0, g_0_xxxxxxy_0_xzzz_1, g_0_xxxxxxyy_0_xxxx_0, g_0_xxxxxxyy_0_xxxy_0, g_0_xxxxxxyy_0_xxxz_0, g_0_xxxxxxyy_0_xxyy_0, g_0_xxxxxxyy_0_xxyz_0, g_0_xxxxxxyy_0_xxzz_0, g_0_xxxxxxyy_0_xyyy_0, g_0_xxxxxxyy_0_xyyz_0, g_0_xxxxxxyy_0_xyzz_0, g_0_xxxxxxyy_0_xzzz_0, g_0_xxxxxxyy_0_yyyy_0, g_0_xxxxxxyy_0_yyyz_0, g_0_xxxxxxyy_0_yyzz_0, g_0_xxxxxxyy_0_yzzz_0, g_0_xxxxxxyy_0_zzzz_0, g_0_xxxxxyy_0_xxxy_0, g_0_xxxxxyy_0_xxxy_1, g_0_xxxxxyy_0_xxy_1, g_0_xxxxxyy_0_xxyy_0, g_0_xxxxxyy_0_xxyy_1, g_0_xxxxxyy_0_xxyz_0, g_0_xxxxxyy_0_xxyz_1, g_0_xxxxxyy_0_xyy_1, g_0_xxxxxyy_0_xyyy_0, g_0_xxxxxyy_0_xyyy_1, g_0_xxxxxyy_0_xyyz_0, g_0_xxxxxyy_0_xyyz_1, g_0_xxxxxyy_0_xyz_1, g_0_xxxxxyy_0_xyzz_0, g_0_xxxxxyy_0_xyzz_1, g_0_xxxxxyy_0_yyy_1, g_0_xxxxxyy_0_yyyy_0, g_0_xxxxxyy_0_yyyy_1, g_0_xxxxxyy_0_yyyz_0, g_0_xxxxxyy_0_yyyz_1, g_0_xxxxxyy_0_yyz_1, g_0_xxxxxyy_0_yyzz_0, g_0_xxxxxyy_0_yyzz_1, g_0_xxxxxyy_0_yzz_1, g_0_xxxxxyy_0_yzzz_0, g_0_xxxxxyy_0_yzzz_1, g_0_xxxxxyy_0_zzzz_0, g_0_xxxxxyy_0_zzzz_1, g_0_xxxxyy_0_xxxy_0, g_0_xxxxyy_0_xxxy_1, g_0_xxxxyy_0_xxyy_0, g_0_xxxxyy_0_xxyy_1, g_0_xxxxyy_0_xxyz_0, g_0_xxxxyy_0_xxyz_1, g_0_xxxxyy_0_xyyy_0, g_0_xxxxyy_0_xyyy_1, g_0_xxxxyy_0_xyyz_0, g_0_xxxxyy_0_xyyz_1, g_0_xxxxyy_0_xyzz_0, g_0_xxxxyy_0_xyzz_1, g_0_xxxxyy_0_yyyy_0, g_0_xxxxyy_0_yyyy_1, g_0_xxxxyy_0_yyyz_0, g_0_xxxxyy_0_yyyz_1, g_0_xxxxyy_0_yyzz_0, g_0_xxxxyy_0_yyzz_1, g_0_xxxxyy_0_yzzz_0, g_0_xxxxyy_0_yzzz_1, g_0_xxxxyy_0_zzzz_0, g_0_xxxxyy_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_xxxx_0[i] = g_0_xxxxxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxx_0[i] * pb_y + g_0_xxxxxxy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxy_0[i] = 5.0 * g_0_xxxxyy_0_xxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxy_0[i] * pb_x + g_0_xxxxxyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxz_0[i] = g_0_xxxxxx_0_xxxz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxz_0[i] * pb_y + g_0_xxxxxxy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxyy_0[i] = 5.0 * g_0_xxxxyy_0_xxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyy_0[i] * pb_x + g_0_xxxxxyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyz_0[i] = 5.0 * g_0_xxxxyy_0_xxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyz_0[i] * pb_x + g_0_xxxxxyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxzz_0[i] = g_0_xxxxxx_0_xxzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxzz_0[i] * pb_y + g_0_xxxxxxy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xyyy_0[i] = 5.0 * g_0_xxxxyy_0_xyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyy_0[i] * pb_x + g_0_xxxxxyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyz_0[i] = 5.0 * g_0_xxxxyy_0_xyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyz_0[i] * pb_x + g_0_xxxxxyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyzz_0[i] = 5.0 * g_0_xxxxyy_0_xyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzz_0[i] * pb_x + g_0_xxxxxyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xzzz_0[i] = g_0_xxxxxx_0_xzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xzzz_0[i] * pb_y + g_0_xxxxxxy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_yyyy_0[i] = 5.0 * g_0_xxxxyy_0_yyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyy_0[i] * pb_x + g_0_xxxxxyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyz_0[i] = 5.0 * g_0_xxxxyy_0_yyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyz_0[i] * pb_x + g_0_xxxxxyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyzz_0[i] = 5.0 * g_0_xxxxyy_0_yyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyzz_0[i] * pb_x + g_0_xxxxxyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yzzz_0[i] = 5.0 * g_0_xxxxyy_0_yzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzzz_0[i] * pb_x + g_0_xxxxxyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_zzzz_0[i] = 5.0 * g_0_xxxxyy_0_zzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_zzzz_0[i] * pb_x + g_0_xxxxxyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 60-75 components of targeted buffer : SLSG

    auto g_0_xxxxxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 60);

    auto g_0_xxxxxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 61);

    auto g_0_xxxxxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 62);

    auto g_0_xxxxxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 63);

    auto g_0_xxxxxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 64);

    auto g_0_xxxxxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 65);

    auto g_0_xxxxxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 66);

    auto g_0_xxxxxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 67);

    auto g_0_xxxxxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 68);

    auto g_0_xxxxxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 69);

    auto g_0_xxxxxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 70);

    auto g_0_xxxxxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 71);

    auto g_0_xxxxxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 72);

    auto g_0_xxxxxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 73);

    auto g_0_xxxxxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 74);

    #pragma omp simd aligned(g_0_xxxxxxy_0_xxxy_0, g_0_xxxxxxy_0_xxxy_1, g_0_xxxxxxy_0_xxyy_0, g_0_xxxxxxy_0_xxyy_1, g_0_xxxxxxy_0_xyyy_0, g_0_xxxxxxy_0_xyyy_1, g_0_xxxxxxy_0_yyyy_0, g_0_xxxxxxy_0_yyyy_1, g_0_xxxxxxyz_0_xxxx_0, g_0_xxxxxxyz_0_xxxy_0, g_0_xxxxxxyz_0_xxxz_0, g_0_xxxxxxyz_0_xxyy_0, g_0_xxxxxxyz_0_xxyz_0, g_0_xxxxxxyz_0_xxzz_0, g_0_xxxxxxyz_0_xyyy_0, g_0_xxxxxxyz_0_xyyz_0, g_0_xxxxxxyz_0_xyzz_0, g_0_xxxxxxyz_0_xzzz_0, g_0_xxxxxxyz_0_yyyy_0, g_0_xxxxxxyz_0_yyyz_0, g_0_xxxxxxyz_0_yyzz_0, g_0_xxxxxxyz_0_yzzz_0, g_0_xxxxxxyz_0_zzzz_0, g_0_xxxxxxz_0_xxxx_0, g_0_xxxxxxz_0_xxxx_1, g_0_xxxxxxz_0_xxxz_0, g_0_xxxxxxz_0_xxxz_1, g_0_xxxxxxz_0_xxyz_0, g_0_xxxxxxz_0_xxyz_1, g_0_xxxxxxz_0_xxz_1, g_0_xxxxxxz_0_xxzz_0, g_0_xxxxxxz_0_xxzz_1, g_0_xxxxxxz_0_xyyz_0, g_0_xxxxxxz_0_xyyz_1, g_0_xxxxxxz_0_xyz_1, g_0_xxxxxxz_0_xyzz_0, g_0_xxxxxxz_0_xyzz_1, g_0_xxxxxxz_0_xzz_1, g_0_xxxxxxz_0_xzzz_0, g_0_xxxxxxz_0_xzzz_1, g_0_xxxxxxz_0_yyyz_0, g_0_xxxxxxz_0_yyyz_1, g_0_xxxxxxz_0_yyz_1, g_0_xxxxxxz_0_yyzz_0, g_0_xxxxxxz_0_yyzz_1, g_0_xxxxxxz_0_yzz_1, g_0_xxxxxxz_0_yzzz_0, g_0_xxxxxxz_0_yzzz_1, g_0_xxxxxxz_0_zzz_1, g_0_xxxxxxz_0_zzzz_0, g_0_xxxxxxz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyz_0_xxxx_0[i] = g_0_xxxxxxz_0_xxxx_0[i] * pb_y + g_0_xxxxxxz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxy_0[i] = g_0_xxxxxxy_0_xxxy_0[i] * pb_z + g_0_xxxxxxy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxz_0[i] = g_0_xxxxxxz_0_xxxz_0[i] * pb_y + g_0_xxxxxxz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyy_0[i] = g_0_xxxxxxy_0_xxyy_0[i] * pb_z + g_0_xxxxxxy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxyz_0[i] = g_0_xxxxxxz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyz_0[i] * pb_y + g_0_xxxxxxz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxzz_0[i] = g_0_xxxxxxz_0_xxzz_0[i] * pb_y + g_0_xxxxxxz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyy_0[i] = g_0_xxxxxxy_0_xyyy_0[i] * pb_z + g_0_xxxxxxy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xyyz_0[i] = 2.0 * g_0_xxxxxxz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyz_0[i] * pb_y + g_0_xxxxxxz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyzz_0[i] = g_0_xxxxxxz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyzz_0[i] * pb_y + g_0_xxxxxxz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xzzz_0[i] = g_0_xxxxxxz_0_xzzz_0[i] * pb_y + g_0_xxxxxxz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyy_0[i] = g_0_xxxxxxy_0_yyyy_0[i] * pb_z + g_0_xxxxxxy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_yyyz_0[i] = 3.0 * g_0_xxxxxxz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyz_0[i] * pb_y + g_0_xxxxxxz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyzz_0[i] = 2.0 * g_0_xxxxxxz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyzz_0[i] * pb_y + g_0_xxxxxxz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yzzz_0[i] = g_0_xxxxxxz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yzzz_0[i] * pb_y + g_0_xxxxxxz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_zzzz_0[i] = g_0_xxxxxxz_0_zzzz_0[i] * pb_y + g_0_xxxxxxz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 75-90 components of targeted buffer : SLSG

    auto g_0_xxxxxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 75);

    auto g_0_xxxxxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 76);

    auto g_0_xxxxxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 77);

    auto g_0_xxxxxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 78);

    auto g_0_xxxxxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 79);

    auto g_0_xxxxxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 80);

    auto g_0_xxxxxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 81);

    auto g_0_xxxxxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 82);

    auto g_0_xxxxxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 83);

    auto g_0_xxxxxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 84);

    auto g_0_xxxxxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 85);

    auto g_0_xxxxxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 86);

    auto g_0_xxxxxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 87);

    auto g_0_xxxxxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 88);

    auto g_0_xxxxxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 89);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxx_0, g_0_xxxxxx_0_xxxx_1, g_0_xxxxxx_0_xxxy_0, g_0_xxxxxx_0_xxxy_1, g_0_xxxxxx_0_xxyy_0, g_0_xxxxxx_0_xxyy_1, g_0_xxxxxx_0_xyyy_0, g_0_xxxxxx_0_xyyy_1, g_0_xxxxxxz_0_xxxx_0, g_0_xxxxxxz_0_xxxx_1, g_0_xxxxxxz_0_xxxy_0, g_0_xxxxxxz_0_xxxy_1, g_0_xxxxxxz_0_xxyy_0, g_0_xxxxxxz_0_xxyy_1, g_0_xxxxxxz_0_xyyy_0, g_0_xxxxxxz_0_xyyy_1, g_0_xxxxxxzz_0_xxxx_0, g_0_xxxxxxzz_0_xxxy_0, g_0_xxxxxxzz_0_xxxz_0, g_0_xxxxxxzz_0_xxyy_0, g_0_xxxxxxzz_0_xxyz_0, g_0_xxxxxxzz_0_xxzz_0, g_0_xxxxxxzz_0_xyyy_0, g_0_xxxxxxzz_0_xyyz_0, g_0_xxxxxxzz_0_xyzz_0, g_0_xxxxxxzz_0_xzzz_0, g_0_xxxxxxzz_0_yyyy_0, g_0_xxxxxxzz_0_yyyz_0, g_0_xxxxxxzz_0_yyzz_0, g_0_xxxxxxzz_0_yzzz_0, g_0_xxxxxxzz_0_zzzz_0, g_0_xxxxxzz_0_xxxz_0, g_0_xxxxxzz_0_xxxz_1, g_0_xxxxxzz_0_xxyz_0, g_0_xxxxxzz_0_xxyz_1, g_0_xxxxxzz_0_xxz_1, g_0_xxxxxzz_0_xxzz_0, g_0_xxxxxzz_0_xxzz_1, g_0_xxxxxzz_0_xyyz_0, g_0_xxxxxzz_0_xyyz_1, g_0_xxxxxzz_0_xyz_1, g_0_xxxxxzz_0_xyzz_0, g_0_xxxxxzz_0_xyzz_1, g_0_xxxxxzz_0_xzz_1, g_0_xxxxxzz_0_xzzz_0, g_0_xxxxxzz_0_xzzz_1, g_0_xxxxxzz_0_yyyy_0, g_0_xxxxxzz_0_yyyy_1, g_0_xxxxxzz_0_yyyz_0, g_0_xxxxxzz_0_yyyz_1, g_0_xxxxxzz_0_yyz_1, g_0_xxxxxzz_0_yyzz_0, g_0_xxxxxzz_0_yyzz_1, g_0_xxxxxzz_0_yzz_1, g_0_xxxxxzz_0_yzzz_0, g_0_xxxxxzz_0_yzzz_1, g_0_xxxxxzz_0_zzz_1, g_0_xxxxxzz_0_zzzz_0, g_0_xxxxxzz_0_zzzz_1, g_0_xxxxzz_0_xxxz_0, g_0_xxxxzz_0_xxxz_1, g_0_xxxxzz_0_xxyz_0, g_0_xxxxzz_0_xxyz_1, g_0_xxxxzz_0_xxzz_0, g_0_xxxxzz_0_xxzz_1, g_0_xxxxzz_0_xyyz_0, g_0_xxxxzz_0_xyyz_1, g_0_xxxxzz_0_xyzz_0, g_0_xxxxzz_0_xyzz_1, g_0_xxxxzz_0_xzzz_0, g_0_xxxxzz_0_xzzz_1, g_0_xxxxzz_0_yyyy_0, g_0_xxxxzz_0_yyyy_1, g_0_xxxxzz_0_yyyz_0, g_0_xxxxzz_0_yyyz_1, g_0_xxxxzz_0_yyzz_0, g_0_xxxxzz_0_yyzz_1, g_0_xxxxzz_0_yzzz_0, g_0_xxxxzz_0_yzzz_1, g_0_xxxxzz_0_zzzz_0, g_0_xxxxzz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_xxxx_0[i] = g_0_xxxxxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxx_0[i] * pb_z + g_0_xxxxxxz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxy_0[i] = g_0_xxxxxx_0_xxxy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxy_0[i] * pb_z + g_0_xxxxxxz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxz_0[i] = 5.0 * g_0_xxxxzz_0_xxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxz_0[i] * pb_x + g_0_xxxxxzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyy_0[i] = g_0_xxxxxx_0_xxyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxyy_0[i] * pb_z + g_0_xxxxxxz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxyz_0[i] = 5.0 * g_0_xxxxzz_0_xxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyz_0[i] * pb_x + g_0_xxxxxzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxzz_0[i] = 5.0 * g_0_xxxxzz_0_xxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxzz_0[i] * pb_x + g_0_xxxxxzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyy_0[i] = g_0_xxxxxx_0_xyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xyyy_0[i] * pb_z + g_0_xxxxxxz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xyyz_0[i] = 5.0 * g_0_xxxxzz_0_xyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyz_0[i] * pb_x + g_0_xxxxxzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyzz_0[i] = 5.0 * g_0_xxxxzz_0_xyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzz_0[i] * pb_x + g_0_xxxxxzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xzzz_0[i] = 5.0 * g_0_xxxxzz_0_xzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xzzz_0[i] * pb_x + g_0_xxxxxzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyy_0[i] = 5.0 * g_0_xxxxzz_0_yyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyy_0[i] * pb_x + g_0_xxxxxzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyz_0[i] = 5.0 * g_0_xxxxzz_0_yyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyz_0[i] * pb_x + g_0_xxxxxzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyzz_0[i] = 5.0 * g_0_xxxxzz_0_yyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyzz_0[i] * pb_x + g_0_xxxxxzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yzzz_0[i] = 5.0 * g_0_xxxxzz_0_yzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzzz_0[i] * pb_x + g_0_xxxxxzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_zzzz_0[i] = 5.0 * g_0_xxxxzz_0_zzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzzz_0[i] * pb_x + g_0_xxxxxzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 90-105 components of targeted buffer : SLSG

    auto g_0_xxxxxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 90);

    auto g_0_xxxxxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 91);

    auto g_0_xxxxxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 92);

    auto g_0_xxxxxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 93);

    auto g_0_xxxxxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 94);

    auto g_0_xxxxxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 95);

    auto g_0_xxxxxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 96);

    auto g_0_xxxxxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 97);

    auto g_0_xxxxxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 98);

    auto g_0_xxxxxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 99);

    auto g_0_xxxxxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 100);

    auto g_0_xxxxxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 101);

    auto g_0_xxxxxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 102);

    auto g_0_xxxxxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 103);

    auto g_0_xxxxxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 104);

    #pragma omp simd aligned(g_0_xxxxxy_0_xxxx_0, g_0_xxxxxy_0_xxxx_1, g_0_xxxxxy_0_xxxz_0, g_0_xxxxxy_0_xxxz_1, g_0_xxxxxy_0_xxzz_0, g_0_xxxxxy_0_xxzz_1, g_0_xxxxxy_0_xzzz_0, g_0_xxxxxy_0_xzzz_1, g_0_xxxxxyy_0_xxxx_0, g_0_xxxxxyy_0_xxxx_1, g_0_xxxxxyy_0_xxxz_0, g_0_xxxxxyy_0_xxxz_1, g_0_xxxxxyy_0_xxzz_0, g_0_xxxxxyy_0_xxzz_1, g_0_xxxxxyy_0_xzzz_0, g_0_xxxxxyy_0_xzzz_1, g_0_xxxxxyyy_0_xxxx_0, g_0_xxxxxyyy_0_xxxy_0, g_0_xxxxxyyy_0_xxxz_0, g_0_xxxxxyyy_0_xxyy_0, g_0_xxxxxyyy_0_xxyz_0, g_0_xxxxxyyy_0_xxzz_0, g_0_xxxxxyyy_0_xyyy_0, g_0_xxxxxyyy_0_xyyz_0, g_0_xxxxxyyy_0_xyzz_0, g_0_xxxxxyyy_0_xzzz_0, g_0_xxxxxyyy_0_yyyy_0, g_0_xxxxxyyy_0_yyyz_0, g_0_xxxxxyyy_0_yyzz_0, g_0_xxxxxyyy_0_yzzz_0, g_0_xxxxxyyy_0_zzzz_0, g_0_xxxxyyy_0_xxxy_0, g_0_xxxxyyy_0_xxxy_1, g_0_xxxxyyy_0_xxy_1, g_0_xxxxyyy_0_xxyy_0, g_0_xxxxyyy_0_xxyy_1, g_0_xxxxyyy_0_xxyz_0, g_0_xxxxyyy_0_xxyz_1, g_0_xxxxyyy_0_xyy_1, g_0_xxxxyyy_0_xyyy_0, g_0_xxxxyyy_0_xyyy_1, g_0_xxxxyyy_0_xyyz_0, g_0_xxxxyyy_0_xyyz_1, g_0_xxxxyyy_0_xyz_1, g_0_xxxxyyy_0_xyzz_0, g_0_xxxxyyy_0_xyzz_1, g_0_xxxxyyy_0_yyy_1, g_0_xxxxyyy_0_yyyy_0, g_0_xxxxyyy_0_yyyy_1, g_0_xxxxyyy_0_yyyz_0, g_0_xxxxyyy_0_yyyz_1, g_0_xxxxyyy_0_yyz_1, g_0_xxxxyyy_0_yyzz_0, g_0_xxxxyyy_0_yyzz_1, g_0_xxxxyyy_0_yzz_1, g_0_xxxxyyy_0_yzzz_0, g_0_xxxxyyy_0_yzzz_1, g_0_xxxxyyy_0_zzzz_0, g_0_xxxxyyy_0_zzzz_1, g_0_xxxyyy_0_xxxy_0, g_0_xxxyyy_0_xxxy_1, g_0_xxxyyy_0_xxyy_0, g_0_xxxyyy_0_xxyy_1, g_0_xxxyyy_0_xxyz_0, g_0_xxxyyy_0_xxyz_1, g_0_xxxyyy_0_xyyy_0, g_0_xxxyyy_0_xyyy_1, g_0_xxxyyy_0_xyyz_0, g_0_xxxyyy_0_xyyz_1, g_0_xxxyyy_0_xyzz_0, g_0_xxxyyy_0_xyzz_1, g_0_xxxyyy_0_yyyy_0, g_0_xxxyyy_0_yyyy_1, g_0_xxxyyy_0_yyyz_0, g_0_xxxyyy_0_yyyz_1, g_0_xxxyyy_0_yyzz_0, g_0_xxxyyy_0_yyzz_1, g_0_xxxyyy_0_yzzz_0, g_0_xxxyyy_0_yzzz_1, g_0_xxxyyy_0_zzzz_0, g_0_xxxyyy_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_xxxx_0[i] = 2.0 * g_0_xxxxxy_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxx_0[i] * pb_y + g_0_xxxxxyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxy_0[i] = 4.0 * g_0_xxxyyy_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxy_0[i] * pb_x + g_0_xxxxyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxz_0[i] = 2.0 * g_0_xxxxxy_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxz_0[i] * pb_y + g_0_xxxxxyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxyy_0[i] = 4.0 * g_0_xxxyyy_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyy_0[i] * pb_x + g_0_xxxxyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyz_0[i] = 4.0 * g_0_xxxyyy_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyz_0[i] * pb_x + g_0_xxxxyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxzz_0[i] = 2.0 * g_0_xxxxxy_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxzz_0[i] * pb_y + g_0_xxxxxyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xyyy_0[i] = 4.0 * g_0_xxxyyy_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyy_0[i] * pb_x + g_0_xxxxyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyz_0[i] = 4.0 * g_0_xxxyyy_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyz_0[i] * pb_x + g_0_xxxxyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyzz_0[i] = 4.0 * g_0_xxxyyy_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzz_0[i] * pb_x + g_0_xxxxyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xzzz_0[i] = 2.0 * g_0_xxxxxy_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xzzz_0[i] * pb_y + g_0_xxxxxyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_yyyy_0[i] = 4.0 * g_0_xxxyyy_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyy_0[i] * pb_x + g_0_xxxxyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyz_0[i] = 4.0 * g_0_xxxyyy_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyz_0[i] * pb_x + g_0_xxxxyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyzz_0[i] = 4.0 * g_0_xxxyyy_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyzz_0[i] * pb_x + g_0_xxxxyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yzzz_0[i] = 4.0 * g_0_xxxyyy_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzzz_0[i] * pb_x + g_0_xxxxyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_zzzz_0[i] = 4.0 * g_0_xxxyyy_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_zzzz_0[i] * pb_x + g_0_xxxxyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 105-120 components of targeted buffer : SLSG

    auto g_0_xxxxxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 105);

    auto g_0_xxxxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 106);

    auto g_0_xxxxxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 107);

    auto g_0_xxxxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 108);

    auto g_0_xxxxxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 109);

    auto g_0_xxxxxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 110);

    auto g_0_xxxxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 111);

    auto g_0_xxxxxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 112);

    auto g_0_xxxxxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 113);

    auto g_0_xxxxxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 114);

    auto g_0_xxxxxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 115);

    auto g_0_xxxxxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 116);

    auto g_0_xxxxxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 117);

    auto g_0_xxxxxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 118);

    auto g_0_xxxxxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 119);

    #pragma omp simd aligned(g_0_xxxxxyy_0_xxx_1, g_0_xxxxxyy_0_xxxx_0, g_0_xxxxxyy_0_xxxx_1, g_0_xxxxxyy_0_xxxy_0, g_0_xxxxxyy_0_xxxy_1, g_0_xxxxxyy_0_xxxz_0, g_0_xxxxxyy_0_xxxz_1, g_0_xxxxxyy_0_xxy_1, g_0_xxxxxyy_0_xxyy_0, g_0_xxxxxyy_0_xxyy_1, g_0_xxxxxyy_0_xxyz_0, g_0_xxxxxyy_0_xxyz_1, g_0_xxxxxyy_0_xxz_1, g_0_xxxxxyy_0_xxzz_0, g_0_xxxxxyy_0_xxzz_1, g_0_xxxxxyy_0_xyy_1, g_0_xxxxxyy_0_xyyy_0, g_0_xxxxxyy_0_xyyy_1, g_0_xxxxxyy_0_xyyz_0, g_0_xxxxxyy_0_xyyz_1, g_0_xxxxxyy_0_xyz_1, g_0_xxxxxyy_0_xyzz_0, g_0_xxxxxyy_0_xyzz_1, g_0_xxxxxyy_0_xzz_1, g_0_xxxxxyy_0_xzzz_0, g_0_xxxxxyy_0_xzzz_1, g_0_xxxxxyy_0_yyy_1, g_0_xxxxxyy_0_yyyy_0, g_0_xxxxxyy_0_yyyy_1, g_0_xxxxxyy_0_yyyz_0, g_0_xxxxxyy_0_yyyz_1, g_0_xxxxxyy_0_yyz_1, g_0_xxxxxyy_0_yyzz_0, g_0_xxxxxyy_0_yyzz_1, g_0_xxxxxyy_0_yzz_1, g_0_xxxxxyy_0_yzzz_0, g_0_xxxxxyy_0_yzzz_1, g_0_xxxxxyy_0_zzz_1, g_0_xxxxxyy_0_zzzz_0, g_0_xxxxxyy_0_zzzz_1, g_0_xxxxxyyz_0_xxxx_0, g_0_xxxxxyyz_0_xxxy_0, g_0_xxxxxyyz_0_xxxz_0, g_0_xxxxxyyz_0_xxyy_0, g_0_xxxxxyyz_0_xxyz_0, g_0_xxxxxyyz_0_xxzz_0, g_0_xxxxxyyz_0_xyyy_0, g_0_xxxxxyyz_0_xyyz_0, g_0_xxxxxyyz_0_xyzz_0, g_0_xxxxxyyz_0_xzzz_0, g_0_xxxxxyyz_0_yyyy_0, g_0_xxxxxyyz_0_yyyz_0, g_0_xxxxxyyz_0_yyzz_0, g_0_xxxxxyyz_0_yzzz_0, g_0_xxxxxyyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_xxxx_0[i] = g_0_xxxxxyy_0_xxxx_0[i] * pb_z + g_0_xxxxxyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxy_0[i] = g_0_xxxxxyy_0_xxxy_0[i] * pb_z + g_0_xxxxxyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxz_0[i] = g_0_xxxxxyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxz_0[i] * pb_z + g_0_xxxxxyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyy_0[i] = g_0_xxxxxyy_0_xxyy_0[i] * pb_z + g_0_xxxxxyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyz_0[i] = g_0_xxxxxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyz_0[i] * pb_z + g_0_xxxxxyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxzz_0[i] * pb_z + g_0_xxxxxyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyy_0[i] = g_0_xxxxxyy_0_xyyy_0[i] * pb_z + g_0_xxxxxyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyz_0[i] = g_0_xxxxxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyz_0[i] * pb_z + g_0_xxxxxyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzz_0[i] * pb_z + g_0_xxxxxyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xzzz_0[i] * pb_z + g_0_xxxxxyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyy_0[i] = g_0_xxxxxyy_0_yyyy_0[i] * pb_z + g_0_xxxxxyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyz_0[i] = g_0_xxxxxyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyz_0[i] * pb_z + g_0_xxxxxyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyzz_0[i] = 2.0 * g_0_xxxxxyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyzz_0[i] * pb_z + g_0_xxxxxyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yzzz_0[i] = 3.0 * g_0_xxxxxyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yzzz_0[i] * pb_z + g_0_xxxxxyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_zzzz_0[i] = 4.0 * g_0_xxxxxyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_zzzz_0[i] * pb_z + g_0_xxxxxyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 120-135 components of targeted buffer : SLSG

    auto g_0_xxxxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 120);

    auto g_0_xxxxxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 121);

    auto g_0_xxxxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 122);

    auto g_0_xxxxxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 123);

    auto g_0_xxxxxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 124);

    auto g_0_xxxxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 125);

    auto g_0_xxxxxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 126);

    auto g_0_xxxxxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 127);

    auto g_0_xxxxxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 128);

    auto g_0_xxxxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 129);

    auto g_0_xxxxxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 130);

    auto g_0_xxxxxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 131);

    auto g_0_xxxxxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 132);

    auto g_0_xxxxxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 133);

    auto g_0_xxxxxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 134);

    #pragma omp simd aligned(g_0_xxxxxyzz_0_xxxx_0, g_0_xxxxxyzz_0_xxxy_0, g_0_xxxxxyzz_0_xxxz_0, g_0_xxxxxyzz_0_xxyy_0, g_0_xxxxxyzz_0_xxyz_0, g_0_xxxxxyzz_0_xxzz_0, g_0_xxxxxyzz_0_xyyy_0, g_0_xxxxxyzz_0_xyyz_0, g_0_xxxxxyzz_0_xyzz_0, g_0_xxxxxyzz_0_xzzz_0, g_0_xxxxxyzz_0_yyyy_0, g_0_xxxxxyzz_0_yyyz_0, g_0_xxxxxyzz_0_yyzz_0, g_0_xxxxxyzz_0_yzzz_0, g_0_xxxxxyzz_0_zzzz_0, g_0_xxxxxzz_0_xxx_1, g_0_xxxxxzz_0_xxxx_0, g_0_xxxxxzz_0_xxxx_1, g_0_xxxxxzz_0_xxxy_0, g_0_xxxxxzz_0_xxxy_1, g_0_xxxxxzz_0_xxxz_0, g_0_xxxxxzz_0_xxxz_1, g_0_xxxxxzz_0_xxy_1, g_0_xxxxxzz_0_xxyy_0, g_0_xxxxxzz_0_xxyy_1, g_0_xxxxxzz_0_xxyz_0, g_0_xxxxxzz_0_xxyz_1, g_0_xxxxxzz_0_xxz_1, g_0_xxxxxzz_0_xxzz_0, g_0_xxxxxzz_0_xxzz_1, g_0_xxxxxzz_0_xyy_1, g_0_xxxxxzz_0_xyyy_0, g_0_xxxxxzz_0_xyyy_1, g_0_xxxxxzz_0_xyyz_0, g_0_xxxxxzz_0_xyyz_1, g_0_xxxxxzz_0_xyz_1, g_0_xxxxxzz_0_xyzz_0, g_0_xxxxxzz_0_xyzz_1, g_0_xxxxxzz_0_xzz_1, g_0_xxxxxzz_0_xzzz_0, g_0_xxxxxzz_0_xzzz_1, g_0_xxxxxzz_0_yyy_1, g_0_xxxxxzz_0_yyyy_0, g_0_xxxxxzz_0_yyyy_1, g_0_xxxxxzz_0_yyyz_0, g_0_xxxxxzz_0_yyyz_1, g_0_xxxxxzz_0_yyz_1, g_0_xxxxxzz_0_yyzz_0, g_0_xxxxxzz_0_yyzz_1, g_0_xxxxxzz_0_yzz_1, g_0_xxxxxzz_0_yzzz_0, g_0_xxxxxzz_0_yzzz_1, g_0_xxxxxzz_0_zzz_1, g_0_xxxxxzz_0_zzzz_0, g_0_xxxxxzz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_xxxx_0[i] = g_0_xxxxxzz_0_xxxx_0[i] * pb_y + g_0_xxxxxzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxy_0[i] = g_0_xxxxxzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxy_0[i] * pb_y + g_0_xxxxxzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxz_0[i] = g_0_xxxxxzz_0_xxxz_0[i] * pb_y + g_0_xxxxxzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyy_0[i] = 2.0 * g_0_xxxxxzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyy_0[i] * pb_y + g_0_xxxxxzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyz_0[i] = g_0_xxxxxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyz_0[i] * pb_y + g_0_xxxxxzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxzz_0[i] = g_0_xxxxxzz_0_xxzz_0[i] * pb_y + g_0_xxxxxzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyy_0[i] = 3.0 * g_0_xxxxxzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyy_0[i] * pb_y + g_0_xxxxxzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyz_0[i] = 2.0 * g_0_xxxxxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyz_0[i] * pb_y + g_0_xxxxxzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyzz_0[i] = g_0_xxxxxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzz_0[i] * pb_y + g_0_xxxxxzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xzzz_0[i] = g_0_xxxxxzz_0_xzzz_0[i] * pb_y + g_0_xxxxxzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyy_0[i] = 4.0 * g_0_xxxxxzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyy_0[i] * pb_y + g_0_xxxxxzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyz_0[i] = 3.0 * g_0_xxxxxzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyz_0[i] * pb_y + g_0_xxxxxzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyzz_0[i] = 2.0 * g_0_xxxxxzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyzz_0[i] * pb_y + g_0_xxxxxzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yzzz_0[i] = g_0_xxxxxzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yzzz_0[i] * pb_y + g_0_xxxxxzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_zzzz_0[i] = g_0_xxxxxzz_0_zzzz_0[i] * pb_y + g_0_xxxxxzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 135-150 components of targeted buffer : SLSG

    auto g_0_xxxxxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 135);

    auto g_0_xxxxxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 136);

    auto g_0_xxxxxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 137);

    auto g_0_xxxxxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 138);

    auto g_0_xxxxxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 139);

    auto g_0_xxxxxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 140);

    auto g_0_xxxxxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 141);

    auto g_0_xxxxxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 142);

    auto g_0_xxxxxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 143);

    auto g_0_xxxxxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 144);

    auto g_0_xxxxxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 145);

    auto g_0_xxxxxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 146);

    auto g_0_xxxxxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 147);

    auto g_0_xxxxxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 148);

    auto g_0_xxxxxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 149);

    #pragma omp simd aligned(g_0_xxxxxz_0_xxxx_0, g_0_xxxxxz_0_xxxx_1, g_0_xxxxxz_0_xxxy_0, g_0_xxxxxz_0_xxxy_1, g_0_xxxxxz_0_xxyy_0, g_0_xxxxxz_0_xxyy_1, g_0_xxxxxz_0_xyyy_0, g_0_xxxxxz_0_xyyy_1, g_0_xxxxxzz_0_xxxx_0, g_0_xxxxxzz_0_xxxx_1, g_0_xxxxxzz_0_xxxy_0, g_0_xxxxxzz_0_xxxy_1, g_0_xxxxxzz_0_xxyy_0, g_0_xxxxxzz_0_xxyy_1, g_0_xxxxxzz_0_xyyy_0, g_0_xxxxxzz_0_xyyy_1, g_0_xxxxxzzz_0_xxxx_0, g_0_xxxxxzzz_0_xxxy_0, g_0_xxxxxzzz_0_xxxz_0, g_0_xxxxxzzz_0_xxyy_0, g_0_xxxxxzzz_0_xxyz_0, g_0_xxxxxzzz_0_xxzz_0, g_0_xxxxxzzz_0_xyyy_0, g_0_xxxxxzzz_0_xyyz_0, g_0_xxxxxzzz_0_xyzz_0, g_0_xxxxxzzz_0_xzzz_0, g_0_xxxxxzzz_0_yyyy_0, g_0_xxxxxzzz_0_yyyz_0, g_0_xxxxxzzz_0_yyzz_0, g_0_xxxxxzzz_0_yzzz_0, g_0_xxxxxzzz_0_zzzz_0, g_0_xxxxzzz_0_xxxz_0, g_0_xxxxzzz_0_xxxz_1, g_0_xxxxzzz_0_xxyz_0, g_0_xxxxzzz_0_xxyz_1, g_0_xxxxzzz_0_xxz_1, g_0_xxxxzzz_0_xxzz_0, g_0_xxxxzzz_0_xxzz_1, g_0_xxxxzzz_0_xyyz_0, g_0_xxxxzzz_0_xyyz_1, g_0_xxxxzzz_0_xyz_1, g_0_xxxxzzz_0_xyzz_0, g_0_xxxxzzz_0_xyzz_1, g_0_xxxxzzz_0_xzz_1, g_0_xxxxzzz_0_xzzz_0, g_0_xxxxzzz_0_xzzz_1, g_0_xxxxzzz_0_yyyy_0, g_0_xxxxzzz_0_yyyy_1, g_0_xxxxzzz_0_yyyz_0, g_0_xxxxzzz_0_yyyz_1, g_0_xxxxzzz_0_yyz_1, g_0_xxxxzzz_0_yyzz_0, g_0_xxxxzzz_0_yyzz_1, g_0_xxxxzzz_0_yzz_1, g_0_xxxxzzz_0_yzzz_0, g_0_xxxxzzz_0_yzzz_1, g_0_xxxxzzz_0_zzz_1, g_0_xxxxzzz_0_zzzz_0, g_0_xxxxzzz_0_zzzz_1, g_0_xxxzzz_0_xxxz_0, g_0_xxxzzz_0_xxxz_1, g_0_xxxzzz_0_xxyz_0, g_0_xxxzzz_0_xxyz_1, g_0_xxxzzz_0_xxzz_0, g_0_xxxzzz_0_xxzz_1, g_0_xxxzzz_0_xyyz_0, g_0_xxxzzz_0_xyyz_1, g_0_xxxzzz_0_xyzz_0, g_0_xxxzzz_0_xyzz_1, g_0_xxxzzz_0_xzzz_0, g_0_xxxzzz_0_xzzz_1, g_0_xxxzzz_0_yyyy_0, g_0_xxxzzz_0_yyyy_1, g_0_xxxzzz_0_yyyz_0, g_0_xxxzzz_0_yyyz_1, g_0_xxxzzz_0_yyzz_0, g_0_xxxzzz_0_yyzz_1, g_0_xxxzzz_0_yzzz_0, g_0_xxxzzz_0_yzzz_1, g_0_xxxzzz_0_zzzz_0, g_0_xxxzzz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_xxxx_0[i] = 2.0 * g_0_xxxxxz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxx_0[i] * pb_z + g_0_xxxxxzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxy_0[i] = 2.0 * g_0_xxxxxz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxy_0[i] * pb_z + g_0_xxxxxzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxz_0[i] = 4.0 * g_0_xxxzzz_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxz_0[i] * pb_x + g_0_xxxxzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyy_0[i] = 2.0 * g_0_xxxxxz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxyy_0[i] * pb_z + g_0_xxxxxzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxyz_0[i] = 4.0 * g_0_xxxzzz_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyz_0[i] * pb_x + g_0_xxxxzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxzz_0[i] = 4.0 * g_0_xxxzzz_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxzz_0[i] * pb_x + g_0_xxxxzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyy_0[i] = 2.0 * g_0_xxxxxz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xyyy_0[i] * pb_z + g_0_xxxxxzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xyyz_0[i] = 4.0 * g_0_xxxzzz_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyz_0[i] * pb_x + g_0_xxxxzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyzz_0[i] = 4.0 * g_0_xxxzzz_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzz_0[i] * pb_x + g_0_xxxxzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xzzz_0[i] = 4.0 * g_0_xxxzzz_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xzzz_0[i] * pb_x + g_0_xxxxzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyy_0[i] = 4.0 * g_0_xxxzzz_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyy_0[i] * pb_x + g_0_xxxxzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyz_0[i] = 4.0 * g_0_xxxzzz_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyz_0[i] * pb_x + g_0_xxxxzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyzz_0[i] = 4.0 * g_0_xxxzzz_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyzz_0[i] * pb_x + g_0_xxxxzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yzzz_0[i] = 4.0 * g_0_xxxzzz_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzzz_0[i] * pb_x + g_0_xxxxzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_zzzz_0[i] = 4.0 * g_0_xxxzzz_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzzz_0[i] * pb_x + g_0_xxxxzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 150-165 components of targeted buffer : SLSG

    auto g_0_xxxxyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 150);

    auto g_0_xxxxyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 151);

    auto g_0_xxxxyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 152);

    auto g_0_xxxxyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 153);

    auto g_0_xxxxyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 154);

    auto g_0_xxxxyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 155);

    auto g_0_xxxxyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 156);

    auto g_0_xxxxyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 157);

    auto g_0_xxxxyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 158);

    auto g_0_xxxxyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 159);

    auto g_0_xxxxyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 160);

    auto g_0_xxxxyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 161);

    auto g_0_xxxxyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 162);

    auto g_0_xxxxyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 163);

    auto g_0_xxxxyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 164);

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxx_0, g_0_xxxxyy_0_xxxx_1, g_0_xxxxyy_0_xxxz_0, g_0_xxxxyy_0_xxxz_1, g_0_xxxxyy_0_xxzz_0, g_0_xxxxyy_0_xxzz_1, g_0_xxxxyy_0_xzzz_0, g_0_xxxxyy_0_xzzz_1, g_0_xxxxyyy_0_xxxx_0, g_0_xxxxyyy_0_xxxx_1, g_0_xxxxyyy_0_xxxz_0, g_0_xxxxyyy_0_xxxz_1, g_0_xxxxyyy_0_xxzz_0, g_0_xxxxyyy_0_xxzz_1, g_0_xxxxyyy_0_xzzz_0, g_0_xxxxyyy_0_xzzz_1, g_0_xxxxyyyy_0_xxxx_0, g_0_xxxxyyyy_0_xxxy_0, g_0_xxxxyyyy_0_xxxz_0, g_0_xxxxyyyy_0_xxyy_0, g_0_xxxxyyyy_0_xxyz_0, g_0_xxxxyyyy_0_xxzz_0, g_0_xxxxyyyy_0_xyyy_0, g_0_xxxxyyyy_0_xyyz_0, g_0_xxxxyyyy_0_xyzz_0, g_0_xxxxyyyy_0_xzzz_0, g_0_xxxxyyyy_0_yyyy_0, g_0_xxxxyyyy_0_yyyz_0, g_0_xxxxyyyy_0_yyzz_0, g_0_xxxxyyyy_0_yzzz_0, g_0_xxxxyyyy_0_zzzz_0, g_0_xxxyyyy_0_xxxy_0, g_0_xxxyyyy_0_xxxy_1, g_0_xxxyyyy_0_xxy_1, g_0_xxxyyyy_0_xxyy_0, g_0_xxxyyyy_0_xxyy_1, g_0_xxxyyyy_0_xxyz_0, g_0_xxxyyyy_0_xxyz_1, g_0_xxxyyyy_0_xyy_1, g_0_xxxyyyy_0_xyyy_0, g_0_xxxyyyy_0_xyyy_1, g_0_xxxyyyy_0_xyyz_0, g_0_xxxyyyy_0_xyyz_1, g_0_xxxyyyy_0_xyz_1, g_0_xxxyyyy_0_xyzz_0, g_0_xxxyyyy_0_xyzz_1, g_0_xxxyyyy_0_yyy_1, g_0_xxxyyyy_0_yyyy_0, g_0_xxxyyyy_0_yyyy_1, g_0_xxxyyyy_0_yyyz_0, g_0_xxxyyyy_0_yyyz_1, g_0_xxxyyyy_0_yyz_1, g_0_xxxyyyy_0_yyzz_0, g_0_xxxyyyy_0_yyzz_1, g_0_xxxyyyy_0_yzz_1, g_0_xxxyyyy_0_yzzz_0, g_0_xxxyyyy_0_yzzz_1, g_0_xxxyyyy_0_zzzz_0, g_0_xxxyyyy_0_zzzz_1, g_0_xxyyyy_0_xxxy_0, g_0_xxyyyy_0_xxxy_1, g_0_xxyyyy_0_xxyy_0, g_0_xxyyyy_0_xxyy_1, g_0_xxyyyy_0_xxyz_0, g_0_xxyyyy_0_xxyz_1, g_0_xxyyyy_0_xyyy_0, g_0_xxyyyy_0_xyyy_1, g_0_xxyyyy_0_xyyz_0, g_0_xxyyyy_0_xyyz_1, g_0_xxyyyy_0_xyzz_0, g_0_xxyyyy_0_xyzz_1, g_0_xxyyyy_0_yyyy_0, g_0_xxyyyy_0_yyyy_1, g_0_xxyyyy_0_yyyz_0, g_0_xxyyyy_0_yyyz_1, g_0_xxyyyy_0_yyzz_0, g_0_xxyyyy_0_yyzz_1, g_0_xxyyyy_0_yzzz_0, g_0_xxyyyy_0_yzzz_1, g_0_xxyyyy_0_zzzz_0, g_0_xxyyyy_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_xxxx_0[i] = 3.0 * g_0_xxxxyy_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxx_0[i] * pb_y + g_0_xxxxyyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxy_0[i] = 3.0 * g_0_xxyyyy_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxy_0[i] * pb_x + g_0_xxxyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxz_0[i] = 3.0 * g_0_xxxxyy_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxz_0[i] * pb_y + g_0_xxxxyyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxyy_0[i] = 3.0 * g_0_xxyyyy_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyy_0[i] * pb_x + g_0_xxxyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyz_0[i] = 3.0 * g_0_xxyyyy_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyz_0[i] * pb_x + g_0_xxxyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxzz_0[i] = 3.0 * g_0_xxxxyy_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxzz_0[i] * pb_y + g_0_xxxxyyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xyyy_0[i] = 3.0 * g_0_xxyyyy_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyy_0[i] * pb_x + g_0_xxxyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyz_0[i] = 3.0 * g_0_xxyyyy_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyz_0[i] * pb_x + g_0_xxxyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyzz_0[i] = 3.0 * g_0_xxyyyy_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzz_0[i] * pb_x + g_0_xxxyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xzzz_0[i] = 3.0 * g_0_xxxxyy_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xzzz_0[i] * pb_y + g_0_xxxxyyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_yyyy_0[i] = 3.0 * g_0_xxyyyy_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyy_0[i] * pb_x + g_0_xxxyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyz_0[i] = 3.0 * g_0_xxyyyy_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyz_0[i] * pb_x + g_0_xxxyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyzz_0[i] = 3.0 * g_0_xxyyyy_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyzz_0[i] * pb_x + g_0_xxxyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yzzz_0[i] = 3.0 * g_0_xxyyyy_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzzz_0[i] * pb_x + g_0_xxxyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_zzzz_0[i] = 3.0 * g_0_xxyyyy_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_zzzz_0[i] * pb_x + g_0_xxxyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 165-180 components of targeted buffer : SLSG

    auto g_0_xxxxyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 165);

    auto g_0_xxxxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 166);

    auto g_0_xxxxyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 167);

    auto g_0_xxxxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 168);

    auto g_0_xxxxyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 169);

    auto g_0_xxxxyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 170);

    auto g_0_xxxxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 171);

    auto g_0_xxxxyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 172);

    auto g_0_xxxxyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 173);

    auto g_0_xxxxyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 174);

    auto g_0_xxxxyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 175);

    auto g_0_xxxxyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 176);

    auto g_0_xxxxyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 177);

    auto g_0_xxxxyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 178);

    auto g_0_xxxxyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 179);

    #pragma omp simd aligned(g_0_xxxxyyy_0_xxx_1, g_0_xxxxyyy_0_xxxx_0, g_0_xxxxyyy_0_xxxx_1, g_0_xxxxyyy_0_xxxy_0, g_0_xxxxyyy_0_xxxy_1, g_0_xxxxyyy_0_xxxz_0, g_0_xxxxyyy_0_xxxz_1, g_0_xxxxyyy_0_xxy_1, g_0_xxxxyyy_0_xxyy_0, g_0_xxxxyyy_0_xxyy_1, g_0_xxxxyyy_0_xxyz_0, g_0_xxxxyyy_0_xxyz_1, g_0_xxxxyyy_0_xxz_1, g_0_xxxxyyy_0_xxzz_0, g_0_xxxxyyy_0_xxzz_1, g_0_xxxxyyy_0_xyy_1, g_0_xxxxyyy_0_xyyy_0, g_0_xxxxyyy_0_xyyy_1, g_0_xxxxyyy_0_xyyz_0, g_0_xxxxyyy_0_xyyz_1, g_0_xxxxyyy_0_xyz_1, g_0_xxxxyyy_0_xyzz_0, g_0_xxxxyyy_0_xyzz_1, g_0_xxxxyyy_0_xzz_1, g_0_xxxxyyy_0_xzzz_0, g_0_xxxxyyy_0_xzzz_1, g_0_xxxxyyy_0_yyy_1, g_0_xxxxyyy_0_yyyy_0, g_0_xxxxyyy_0_yyyy_1, g_0_xxxxyyy_0_yyyz_0, g_0_xxxxyyy_0_yyyz_1, g_0_xxxxyyy_0_yyz_1, g_0_xxxxyyy_0_yyzz_0, g_0_xxxxyyy_0_yyzz_1, g_0_xxxxyyy_0_yzz_1, g_0_xxxxyyy_0_yzzz_0, g_0_xxxxyyy_0_yzzz_1, g_0_xxxxyyy_0_zzz_1, g_0_xxxxyyy_0_zzzz_0, g_0_xxxxyyy_0_zzzz_1, g_0_xxxxyyyz_0_xxxx_0, g_0_xxxxyyyz_0_xxxy_0, g_0_xxxxyyyz_0_xxxz_0, g_0_xxxxyyyz_0_xxyy_0, g_0_xxxxyyyz_0_xxyz_0, g_0_xxxxyyyz_0_xxzz_0, g_0_xxxxyyyz_0_xyyy_0, g_0_xxxxyyyz_0_xyyz_0, g_0_xxxxyyyz_0_xyzz_0, g_0_xxxxyyyz_0_xzzz_0, g_0_xxxxyyyz_0_yyyy_0, g_0_xxxxyyyz_0_yyyz_0, g_0_xxxxyyyz_0_yyzz_0, g_0_xxxxyyyz_0_yzzz_0, g_0_xxxxyyyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_xxxx_0[i] = g_0_xxxxyyy_0_xxxx_0[i] * pb_z + g_0_xxxxyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxy_0[i] = g_0_xxxxyyy_0_xxxy_0[i] * pb_z + g_0_xxxxyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxz_0[i] = g_0_xxxxyyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxz_0[i] * pb_z + g_0_xxxxyyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyy_0[i] = g_0_xxxxyyy_0_xxyy_0[i] * pb_z + g_0_xxxxyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyz_0[i] = g_0_xxxxyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyz_0[i] * pb_z + g_0_xxxxyyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxzz_0[i] * pb_z + g_0_xxxxyyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyy_0[i] = g_0_xxxxyyy_0_xyyy_0[i] * pb_z + g_0_xxxxyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyz_0[i] = g_0_xxxxyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyz_0[i] * pb_z + g_0_xxxxyyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzz_0[i] * pb_z + g_0_xxxxyyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xzzz_0[i] * pb_z + g_0_xxxxyyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyy_0[i] = g_0_xxxxyyy_0_yyyy_0[i] * pb_z + g_0_xxxxyyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyz_0[i] = g_0_xxxxyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyz_0[i] * pb_z + g_0_xxxxyyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyzz_0[i] = 2.0 * g_0_xxxxyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyzz_0[i] * pb_z + g_0_xxxxyyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yzzz_0[i] = 3.0 * g_0_xxxxyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yzzz_0[i] * pb_z + g_0_xxxxyyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_zzzz_0[i] = 4.0 * g_0_xxxxyyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_zzzz_0[i] * pb_z + g_0_xxxxyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 180-195 components of targeted buffer : SLSG

    auto g_0_xxxxyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 180);

    auto g_0_xxxxyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 181);

    auto g_0_xxxxyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 182);

    auto g_0_xxxxyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 183);

    auto g_0_xxxxyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 184);

    auto g_0_xxxxyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 185);

    auto g_0_xxxxyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 186);

    auto g_0_xxxxyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 187);

    auto g_0_xxxxyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 188);

    auto g_0_xxxxyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 189);

    auto g_0_xxxxyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 190);

    auto g_0_xxxxyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 191);

    auto g_0_xxxxyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 192);

    auto g_0_xxxxyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 193);

    auto g_0_xxxxyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 194);

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxy_0, g_0_xxxxyy_0_xxxy_1, g_0_xxxxyy_0_xxyy_0, g_0_xxxxyy_0_xxyy_1, g_0_xxxxyy_0_xyyy_0, g_0_xxxxyy_0_xyyy_1, g_0_xxxxyyz_0_xxxy_0, g_0_xxxxyyz_0_xxxy_1, g_0_xxxxyyz_0_xxyy_0, g_0_xxxxyyz_0_xxyy_1, g_0_xxxxyyz_0_xyyy_0, g_0_xxxxyyz_0_xyyy_1, g_0_xxxxyyzz_0_xxxx_0, g_0_xxxxyyzz_0_xxxy_0, g_0_xxxxyyzz_0_xxxz_0, g_0_xxxxyyzz_0_xxyy_0, g_0_xxxxyyzz_0_xxyz_0, g_0_xxxxyyzz_0_xxzz_0, g_0_xxxxyyzz_0_xyyy_0, g_0_xxxxyyzz_0_xyyz_0, g_0_xxxxyyzz_0_xyzz_0, g_0_xxxxyyzz_0_xzzz_0, g_0_xxxxyyzz_0_yyyy_0, g_0_xxxxyyzz_0_yyyz_0, g_0_xxxxyyzz_0_yyzz_0, g_0_xxxxyyzz_0_yzzz_0, g_0_xxxxyyzz_0_zzzz_0, g_0_xxxxyzz_0_xxxx_0, g_0_xxxxyzz_0_xxxx_1, g_0_xxxxyzz_0_xxxz_0, g_0_xxxxyzz_0_xxxz_1, g_0_xxxxyzz_0_xxzz_0, g_0_xxxxyzz_0_xxzz_1, g_0_xxxxyzz_0_xzzz_0, g_0_xxxxyzz_0_xzzz_1, g_0_xxxxzz_0_xxxx_0, g_0_xxxxzz_0_xxxx_1, g_0_xxxxzz_0_xxxz_0, g_0_xxxxzz_0_xxxz_1, g_0_xxxxzz_0_xxzz_0, g_0_xxxxzz_0_xxzz_1, g_0_xxxxzz_0_xzzz_0, g_0_xxxxzz_0_xzzz_1, g_0_xxxyyzz_0_xxyz_0, g_0_xxxyyzz_0_xxyz_1, g_0_xxxyyzz_0_xyyz_0, g_0_xxxyyzz_0_xyyz_1, g_0_xxxyyzz_0_xyz_1, g_0_xxxyyzz_0_xyzz_0, g_0_xxxyyzz_0_xyzz_1, g_0_xxxyyzz_0_yyyy_0, g_0_xxxyyzz_0_yyyy_1, g_0_xxxyyzz_0_yyyz_0, g_0_xxxyyzz_0_yyyz_1, g_0_xxxyyzz_0_yyz_1, g_0_xxxyyzz_0_yyzz_0, g_0_xxxyyzz_0_yyzz_1, g_0_xxxyyzz_0_yzz_1, g_0_xxxyyzz_0_yzzz_0, g_0_xxxyyzz_0_yzzz_1, g_0_xxxyyzz_0_zzzz_0, g_0_xxxyyzz_0_zzzz_1, g_0_xxyyzz_0_xxyz_0, g_0_xxyyzz_0_xxyz_1, g_0_xxyyzz_0_xyyz_0, g_0_xxyyzz_0_xyyz_1, g_0_xxyyzz_0_xyzz_0, g_0_xxyyzz_0_xyzz_1, g_0_xxyyzz_0_yyyy_0, g_0_xxyyzz_0_yyyy_1, g_0_xxyyzz_0_yyyz_0, g_0_xxyyzz_0_yyyz_1, g_0_xxyyzz_0_yyzz_0, g_0_xxyyzz_0_yyzz_1, g_0_xxyyzz_0_yzzz_0, g_0_xxyyzz_0_yzzz_1, g_0_xxyyzz_0_zzzz_0, g_0_xxyyzz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_xxxx_0[i] = g_0_xxxxzz_0_xxxx_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxx_0[i] * pb_y + g_0_xxxxyzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxy_0[i] = g_0_xxxxyy_0_xxxy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxy_0[i] * pb_z + g_0_xxxxyyz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxz_0[i] = g_0_xxxxzz_0_xxxz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxz_0[i] * pb_y + g_0_xxxxyzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxyy_0[i] = g_0_xxxxyy_0_xxyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxyy_0[i] * pb_z + g_0_xxxxyyz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxyz_0[i] = 3.0 * g_0_xxyyzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyz_0[i] * pb_x + g_0_xxxyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxzz_0[i] = g_0_xxxxzz_0_xxzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxzz_0[i] * pb_y + g_0_xxxxyzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xyyy_0[i] = g_0_xxxxyy_0_xyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xyyy_0[i] * pb_z + g_0_xxxxyyz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xyyz_0[i] = 3.0 * g_0_xxyyzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyz_0[i] * pb_x + g_0_xxxyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyzz_0[i] = 3.0 * g_0_xxyyzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyzz_0[i] * pb_x + g_0_xxxyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xzzz_0[i] = g_0_xxxxzz_0_xzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xzzz_0[i] * pb_y + g_0_xxxxyzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_yyyy_0[i] = 3.0 * g_0_xxyyzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyy_0[i] * pb_x + g_0_xxxyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyz_0[i] = 3.0 * g_0_xxyyzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyz_0[i] * pb_x + g_0_xxxyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyzz_0[i] = 3.0 * g_0_xxyyzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyzz_0[i] * pb_x + g_0_xxxyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yzzz_0[i] = 3.0 * g_0_xxyyzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzzz_0[i] * pb_x + g_0_xxxyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_zzzz_0[i] = 3.0 * g_0_xxyyzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_zzzz_0[i] * pb_x + g_0_xxxyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 195-210 components of targeted buffer : SLSG

    auto g_0_xxxxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 195);

    auto g_0_xxxxyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 196);

    auto g_0_xxxxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 197);

    auto g_0_xxxxyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 198);

    auto g_0_xxxxyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 199);

    auto g_0_xxxxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 200);

    auto g_0_xxxxyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 201);

    auto g_0_xxxxyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 202);

    auto g_0_xxxxyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 203);

    auto g_0_xxxxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 204);

    auto g_0_xxxxyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 205);

    auto g_0_xxxxyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 206);

    auto g_0_xxxxyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 207);

    auto g_0_xxxxyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 208);

    auto g_0_xxxxyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 209);

    #pragma omp simd aligned(g_0_xxxxyzzz_0_xxxx_0, g_0_xxxxyzzz_0_xxxy_0, g_0_xxxxyzzz_0_xxxz_0, g_0_xxxxyzzz_0_xxyy_0, g_0_xxxxyzzz_0_xxyz_0, g_0_xxxxyzzz_0_xxzz_0, g_0_xxxxyzzz_0_xyyy_0, g_0_xxxxyzzz_0_xyyz_0, g_0_xxxxyzzz_0_xyzz_0, g_0_xxxxyzzz_0_xzzz_0, g_0_xxxxyzzz_0_yyyy_0, g_0_xxxxyzzz_0_yyyz_0, g_0_xxxxyzzz_0_yyzz_0, g_0_xxxxyzzz_0_yzzz_0, g_0_xxxxyzzz_0_zzzz_0, g_0_xxxxzzz_0_xxx_1, g_0_xxxxzzz_0_xxxx_0, g_0_xxxxzzz_0_xxxx_1, g_0_xxxxzzz_0_xxxy_0, g_0_xxxxzzz_0_xxxy_1, g_0_xxxxzzz_0_xxxz_0, g_0_xxxxzzz_0_xxxz_1, g_0_xxxxzzz_0_xxy_1, g_0_xxxxzzz_0_xxyy_0, g_0_xxxxzzz_0_xxyy_1, g_0_xxxxzzz_0_xxyz_0, g_0_xxxxzzz_0_xxyz_1, g_0_xxxxzzz_0_xxz_1, g_0_xxxxzzz_0_xxzz_0, g_0_xxxxzzz_0_xxzz_1, g_0_xxxxzzz_0_xyy_1, g_0_xxxxzzz_0_xyyy_0, g_0_xxxxzzz_0_xyyy_1, g_0_xxxxzzz_0_xyyz_0, g_0_xxxxzzz_0_xyyz_1, g_0_xxxxzzz_0_xyz_1, g_0_xxxxzzz_0_xyzz_0, g_0_xxxxzzz_0_xyzz_1, g_0_xxxxzzz_0_xzz_1, g_0_xxxxzzz_0_xzzz_0, g_0_xxxxzzz_0_xzzz_1, g_0_xxxxzzz_0_yyy_1, g_0_xxxxzzz_0_yyyy_0, g_0_xxxxzzz_0_yyyy_1, g_0_xxxxzzz_0_yyyz_0, g_0_xxxxzzz_0_yyyz_1, g_0_xxxxzzz_0_yyz_1, g_0_xxxxzzz_0_yyzz_0, g_0_xxxxzzz_0_yyzz_1, g_0_xxxxzzz_0_yzz_1, g_0_xxxxzzz_0_yzzz_0, g_0_xxxxzzz_0_yzzz_1, g_0_xxxxzzz_0_zzz_1, g_0_xxxxzzz_0_zzzz_0, g_0_xxxxzzz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_xxxx_0[i] = g_0_xxxxzzz_0_xxxx_0[i] * pb_y + g_0_xxxxzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxy_0[i] = g_0_xxxxzzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxy_0[i] * pb_y + g_0_xxxxzzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxz_0[i] = g_0_xxxxzzz_0_xxxz_0[i] * pb_y + g_0_xxxxzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyy_0[i] = 2.0 * g_0_xxxxzzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyy_0[i] * pb_y + g_0_xxxxzzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyz_0[i] = g_0_xxxxzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyz_0[i] * pb_y + g_0_xxxxzzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxzz_0[i] = g_0_xxxxzzz_0_xxzz_0[i] * pb_y + g_0_xxxxzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyy_0[i] = 3.0 * g_0_xxxxzzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyy_0[i] * pb_y + g_0_xxxxzzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyz_0[i] = 2.0 * g_0_xxxxzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyz_0[i] * pb_y + g_0_xxxxzzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyzz_0[i] = g_0_xxxxzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzz_0[i] * pb_y + g_0_xxxxzzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xzzz_0[i] = g_0_xxxxzzz_0_xzzz_0[i] * pb_y + g_0_xxxxzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyy_0[i] = 4.0 * g_0_xxxxzzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyy_0[i] * pb_y + g_0_xxxxzzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyz_0[i] = 3.0 * g_0_xxxxzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyz_0[i] * pb_y + g_0_xxxxzzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyzz_0[i] = 2.0 * g_0_xxxxzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyzz_0[i] * pb_y + g_0_xxxxzzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yzzz_0[i] = g_0_xxxxzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yzzz_0[i] * pb_y + g_0_xxxxzzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_zzzz_0[i] = g_0_xxxxzzz_0_zzzz_0[i] * pb_y + g_0_xxxxzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 210-225 components of targeted buffer : SLSG

    auto g_0_xxxxzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 210);

    auto g_0_xxxxzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 211);

    auto g_0_xxxxzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 212);

    auto g_0_xxxxzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 213);

    auto g_0_xxxxzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 214);

    auto g_0_xxxxzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 215);

    auto g_0_xxxxzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 216);

    auto g_0_xxxxzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 217);

    auto g_0_xxxxzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 218);

    auto g_0_xxxxzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 219);

    auto g_0_xxxxzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 220);

    auto g_0_xxxxzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 221);

    auto g_0_xxxxzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 222);

    auto g_0_xxxxzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 223);

    auto g_0_xxxxzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 224);

    #pragma omp simd aligned(g_0_xxxxzz_0_xxxx_0, g_0_xxxxzz_0_xxxx_1, g_0_xxxxzz_0_xxxy_0, g_0_xxxxzz_0_xxxy_1, g_0_xxxxzz_0_xxyy_0, g_0_xxxxzz_0_xxyy_1, g_0_xxxxzz_0_xyyy_0, g_0_xxxxzz_0_xyyy_1, g_0_xxxxzzz_0_xxxx_0, g_0_xxxxzzz_0_xxxx_1, g_0_xxxxzzz_0_xxxy_0, g_0_xxxxzzz_0_xxxy_1, g_0_xxxxzzz_0_xxyy_0, g_0_xxxxzzz_0_xxyy_1, g_0_xxxxzzz_0_xyyy_0, g_0_xxxxzzz_0_xyyy_1, g_0_xxxxzzzz_0_xxxx_0, g_0_xxxxzzzz_0_xxxy_0, g_0_xxxxzzzz_0_xxxz_0, g_0_xxxxzzzz_0_xxyy_0, g_0_xxxxzzzz_0_xxyz_0, g_0_xxxxzzzz_0_xxzz_0, g_0_xxxxzzzz_0_xyyy_0, g_0_xxxxzzzz_0_xyyz_0, g_0_xxxxzzzz_0_xyzz_0, g_0_xxxxzzzz_0_xzzz_0, g_0_xxxxzzzz_0_yyyy_0, g_0_xxxxzzzz_0_yyyz_0, g_0_xxxxzzzz_0_yyzz_0, g_0_xxxxzzzz_0_yzzz_0, g_0_xxxxzzzz_0_zzzz_0, g_0_xxxzzzz_0_xxxz_0, g_0_xxxzzzz_0_xxxz_1, g_0_xxxzzzz_0_xxyz_0, g_0_xxxzzzz_0_xxyz_1, g_0_xxxzzzz_0_xxz_1, g_0_xxxzzzz_0_xxzz_0, g_0_xxxzzzz_0_xxzz_1, g_0_xxxzzzz_0_xyyz_0, g_0_xxxzzzz_0_xyyz_1, g_0_xxxzzzz_0_xyz_1, g_0_xxxzzzz_0_xyzz_0, g_0_xxxzzzz_0_xyzz_1, g_0_xxxzzzz_0_xzz_1, g_0_xxxzzzz_0_xzzz_0, g_0_xxxzzzz_0_xzzz_1, g_0_xxxzzzz_0_yyyy_0, g_0_xxxzzzz_0_yyyy_1, g_0_xxxzzzz_0_yyyz_0, g_0_xxxzzzz_0_yyyz_1, g_0_xxxzzzz_0_yyz_1, g_0_xxxzzzz_0_yyzz_0, g_0_xxxzzzz_0_yyzz_1, g_0_xxxzzzz_0_yzz_1, g_0_xxxzzzz_0_yzzz_0, g_0_xxxzzzz_0_yzzz_1, g_0_xxxzzzz_0_zzz_1, g_0_xxxzzzz_0_zzzz_0, g_0_xxxzzzz_0_zzzz_1, g_0_xxzzzz_0_xxxz_0, g_0_xxzzzz_0_xxxz_1, g_0_xxzzzz_0_xxyz_0, g_0_xxzzzz_0_xxyz_1, g_0_xxzzzz_0_xxzz_0, g_0_xxzzzz_0_xxzz_1, g_0_xxzzzz_0_xyyz_0, g_0_xxzzzz_0_xyyz_1, g_0_xxzzzz_0_xyzz_0, g_0_xxzzzz_0_xyzz_1, g_0_xxzzzz_0_xzzz_0, g_0_xxzzzz_0_xzzz_1, g_0_xxzzzz_0_yyyy_0, g_0_xxzzzz_0_yyyy_1, g_0_xxzzzz_0_yyyz_0, g_0_xxzzzz_0_yyyz_1, g_0_xxzzzz_0_yyzz_0, g_0_xxzzzz_0_yyzz_1, g_0_xxzzzz_0_yzzz_0, g_0_xxzzzz_0_yzzz_1, g_0_xxzzzz_0_zzzz_0, g_0_xxzzzz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_xxxx_0[i] = 3.0 * g_0_xxxxzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxx_0[i] * pb_z + g_0_xxxxzzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxy_0[i] = 3.0 * g_0_xxxxzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxy_0[i] * pb_z + g_0_xxxxzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxz_0[i] = 3.0 * g_0_xxzzzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxz_0[i] * pb_x + g_0_xxxzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyy_0[i] = 3.0 * g_0_xxxxzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxyy_0[i] * pb_z + g_0_xxxxzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxyz_0[i] = 3.0 * g_0_xxzzzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyz_0[i] * pb_x + g_0_xxxzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxzz_0[i] = 3.0 * g_0_xxzzzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxzz_0[i] * pb_x + g_0_xxxzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyy_0[i] = 3.0 * g_0_xxxxzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xyyy_0[i] * pb_z + g_0_xxxxzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xyyz_0[i] = 3.0 * g_0_xxzzzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyz_0[i] * pb_x + g_0_xxxzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyzz_0[i] = 3.0 * g_0_xxzzzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzz_0[i] * pb_x + g_0_xxxzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xzzz_0[i] = 3.0 * g_0_xxzzzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xzzz_0[i] * pb_x + g_0_xxxzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyy_0[i] = 3.0 * g_0_xxzzzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyy_0[i] * pb_x + g_0_xxxzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyz_0[i] = 3.0 * g_0_xxzzzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyz_0[i] * pb_x + g_0_xxxzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyzz_0[i] * pb_x + g_0_xxxzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yzzz_0[i] = 3.0 * g_0_xxzzzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzzz_0[i] * pb_x + g_0_xxxzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_zzzz_0[i] = 3.0 * g_0_xxzzzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzzz_0[i] * pb_x + g_0_xxxzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 225-240 components of targeted buffer : SLSG

    auto g_0_xxxyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 225);

    auto g_0_xxxyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 226);

    auto g_0_xxxyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 227);

    auto g_0_xxxyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 228);

    auto g_0_xxxyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 229);

    auto g_0_xxxyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 230);

    auto g_0_xxxyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 231);

    auto g_0_xxxyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 232);

    auto g_0_xxxyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 233);

    auto g_0_xxxyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 234);

    auto g_0_xxxyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 235);

    auto g_0_xxxyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 236);

    auto g_0_xxxyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 237);

    auto g_0_xxxyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 238);

    auto g_0_xxxyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 239);

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxx_0, g_0_xxxyyy_0_xxxx_1, g_0_xxxyyy_0_xxxz_0, g_0_xxxyyy_0_xxxz_1, g_0_xxxyyy_0_xxzz_0, g_0_xxxyyy_0_xxzz_1, g_0_xxxyyy_0_xzzz_0, g_0_xxxyyy_0_xzzz_1, g_0_xxxyyyy_0_xxxx_0, g_0_xxxyyyy_0_xxxx_1, g_0_xxxyyyy_0_xxxz_0, g_0_xxxyyyy_0_xxxz_1, g_0_xxxyyyy_0_xxzz_0, g_0_xxxyyyy_0_xxzz_1, g_0_xxxyyyy_0_xzzz_0, g_0_xxxyyyy_0_xzzz_1, g_0_xxxyyyyy_0_xxxx_0, g_0_xxxyyyyy_0_xxxy_0, g_0_xxxyyyyy_0_xxxz_0, g_0_xxxyyyyy_0_xxyy_0, g_0_xxxyyyyy_0_xxyz_0, g_0_xxxyyyyy_0_xxzz_0, g_0_xxxyyyyy_0_xyyy_0, g_0_xxxyyyyy_0_xyyz_0, g_0_xxxyyyyy_0_xyzz_0, g_0_xxxyyyyy_0_xzzz_0, g_0_xxxyyyyy_0_yyyy_0, g_0_xxxyyyyy_0_yyyz_0, g_0_xxxyyyyy_0_yyzz_0, g_0_xxxyyyyy_0_yzzz_0, g_0_xxxyyyyy_0_zzzz_0, g_0_xxyyyyy_0_xxxy_0, g_0_xxyyyyy_0_xxxy_1, g_0_xxyyyyy_0_xxy_1, g_0_xxyyyyy_0_xxyy_0, g_0_xxyyyyy_0_xxyy_1, g_0_xxyyyyy_0_xxyz_0, g_0_xxyyyyy_0_xxyz_1, g_0_xxyyyyy_0_xyy_1, g_0_xxyyyyy_0_xyyy_0, g_0_xxyyyyy_0_xyyy_1, g_0_xxyyyyy_0_xyyz_0, g_0_xxyyyyy_0_xyyz_1, g_0_xxyyyyy_0_xyz_1, g_0_xxyyyyy_0_xyzz_0, g_0_xxyyyyy_0_xyzz_1, g_0_xxyyyyy_0_yyy_1, g_0_xxyyyyy_0_yyyy_0, g_0_xxyyyyy_0_yyyy_1, g_0_xxyyyyy_0_yyyz_0, g_0_xxyyyyy_0_yyyz_1, g_0_xxyyyyy_0_yyz_1, g_0_xxyyyyy_0_yyzz_0, g_0_xxyyyyy_0_yyzz_1, g_0_xxyyyyy_0_yzz_1, g_0_xxyyyyy_0_yzzz_0, g_0_xxyyyyy_0_yzzz_1, g_0_xxyyyyy_0_zzzz_0, g_0_xxyyyyy_0_zzzz_1, g_0_xyyyyy_0_xxxy_0, g_0_xyyyyy_0_xxxy_1, g_0_xyyyyy_0_xxyy_0, g_0_xyyyyy_0_xxyy_1, g_0_xyyyyy_0_xxyz_0, g_0_xyyyyy_0_xxyz_1, g_0_xyyyyy_0_xyyy_0, g_0_xyyyyy_0_xyyy_1, g_0_xyyyyy_0_xyyz_0, g_0_xyyyyy_0_xyyz_1, g_0_xyyyyy_0_xyzz_0, g_0_xyyyyy_0_xyzz_1, g_0_xyyyyy_0_yyyy_0, g_0_xyyyyy_0_yyyy_1, g_0_xyyyyy_0_yyyz_0, g_0_xyyyyy_0_yyyz_1, g_0_xyyyyy_0_yyzz_0, g_0_xyyyyy_0_yyzz_1, g_0_xyyyyy_0_yzzz_0, g_0_xyyyyy_0_yzzz_1, g_0_xyyyyy_0_zzzz_0, g_0_xyyyyy_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_xxxx_0[i] = 4.0 * g_0_xxxyyy_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxx_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxx_0[i] * pb_y + g_0_xxxyyyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxy_0[i] = 2.0 * g_0_xyyyyy_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxy_0[i] * pb_x + g_0_xxyyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxz_0[i] = 4.0 * g_0_xxxyyy_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxz_0[i] * pb_y + g_0_xxxyyyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxyy_0[i] = 2.0 * g_0_xyyyyy_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyy_0[i] * pb_x + g_0_xxyyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyz_0[i] = 2.0 * g_0_xyyyyy_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyz_0[i] * pb_x + g_0_xxyyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxzz_0[i] = 4.0 * g_0_xxxyyy_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxzz_0[i] * pb_y + g_0_xxxyyyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xyyy_0[i] = 2.0 * g_0_xyyyyy_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyy_0[i] * pb_x + g_0_xxyyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyz_0[i] = 2.0 * g_0_xyyyyy_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyz_0[i] * pb_x + g_0_xxyyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyzz_0[i] = 2.0 * g_0_xyyyyy_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzz_0[i] * pb_x + g_0_xxyyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xzzz_0[i] = 4.0 * g_0_xxxyyy_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xzzz_0[i] * pb_y + g_0_xxxyyyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_yyyy_0[i] = 2.0 * g_0_xyyyyy_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyy_0[i] * pb_x + g_0_xxyyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyz_0[i] = 2.0 * g_0_xyyyyy_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyz_0[i] * pb_x + g_0_xxyyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyzz_0[i] = 2.0 * g_0_xyyyyy_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyzz_0[i] * pb_x + g_0_xxyyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yzzz_0[i] = 2.0 * g_0_xyyyyy_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzzz_0[i] * pb_x + g_0_xxyyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_zzzz_0[i] = 2.0 * g_0_xyyyyy_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_zzzz_0[i] * pb_x + g_0_xxyyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 240-255 components of targeted buffer : SLSG

    auto g_0_xxxyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 240);

    auto g_0_xxxyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 241);

    auto g_0_xxxyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 242);

    auto g_0_xxxyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 243);

    auto g_0_xxxyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 244);

    auto g_0_xxxyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 245);

    auto g_0_xxxyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 246);

    auto g_0_xxxyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 247);

    auto g_0_xxxyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 248);

    auto g_0_xxxyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 249);

    auto g_0_xxxyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 250);

    auto g_0_xxxyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 251);

    auto g_0_xxxyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 252);

    auto g_0_xxxyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 253);

    auto g_0_xxxyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 254);

    #pragma omp simd aligned(g_0_xxxyyyy_0_xxx_1, g_0_xxxyyyy_0_xxxx_0, g_0_xxxyyyy_0_xxxx_1, g_0_xxxyyyy_0_xxxy_0, g_0_xxxyyyy_0_xxxy_1, g_0_xxxyyyy_0_xxxz_0, g_0_xxxyyyy_0_xxxz_1, g_0_xxxyyyy_0_xxy_1, g_0_xxxyyyy_0_xxyy_0, g_0_xxxyyyy_0_xxyy_1, g_0_xxxyyyy_0_xxyz_0, g_0_xxxyyyy_0_xxyz_1, g_0_xxxyyyy_0_xxz_1, g_0_xxxyyyy_0_xxzz_0, g_0_xxxyyyy_0_xxzz_1, g_0_xxxyyyy_0_xyy_1, g_0_xxxyyyy_0_xyyy_0, g_0_xxxyyyy_0_xyyy_1, g_0_xxxyyyy_0_xyyz_0, g_0_xxxyyyy_0_xyyz_1, g_0_xxxyyyy_0_xyz_1, g_0_xxxyyyy_0_xyzz_0, g_0_xxxyyyy_0_xyzz_1, g_0_xxxyyyy_0_xzz_1, g_0_xxxyyyy_0_xzzz_0, g_0_xxxyyyy_0_xzzz_1, g_0_xxxyyyy_0_yyy_1, g_0_xxxyyyy_0_yyyy_0, g_0_xxxyyyy_0_yyyy_1, g_0_xxxyyyy_0_yyyz_0, g_0_xxxyyyy_0_yyyz_1, g_0_xxxyyyy_0_yyz_1, g_0_xxxyyyy_0_yyzz_0, g_0_xxxyyyy_0_yyzz_1, g_0_xxxyyyy_0_yzz_1, g_0_xxxyyyy_0_yzzz_0, g_0_xxxyyyy_0_yzzz_1, g_0_xxxyyyy_0_zzz_1, g_0_xxxyyyy_0_zzzz_0, g_0_xxxyyyy_0_zzzz_1, g_0_xxxyyyyz_0_xxxx_0, g_0_xxxyyyyz_0_xxxy_0, g_0_xxxyyyyz_0_xxxz_0, g_0_xxxyyyyz_0_xxyy_0, g_0_xxxyyyyz_0_xxyz_0, g_0_xxxyyyyz_0_xxzz_0, g_0_xxxyyyyz_0_xyyy_0, g_0_xxxyyyyz_0_xyyz_0, g_0_xxxyyyyz_0_xyzz_0, g_0_xxxyyyyz_0_xzzz_0, g_0_xxxyyyyz_0_yyyy_0, g_0_xxxyyyyz_0_yyyz_0, g_0_xxxyyyyz_0_yyzz_0, g_0_xxxyyyyz_0_yzzz_0, g_0_xxxyyyyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_xxxx_0[i] = g_0_xxxyyyy_0_xxxx_0[i] * pb_z + g_0_xxxyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxy_0[i] = g_0_xxxyyyy_0_xxxy_0[i] * pb_z + g_0_xxxyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxz_0[i] = g_0_xxxyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxz_0[i] * pb_z + g_0_xxxyyyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyy_0[i] = g_0_xxxyyyy_0_xxyy_0[i] * pb_z + g_0_xxxyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyz_0[i] = g_0_xxxyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyz_0[i] * pb_z + g_0_xxxyyyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxzz_0[i] * pb_z + g_0_xxxyyyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyy_0[i] = g_0_xxxyyyy_0_xyyy_0[i] * pb_z + g_0_xxxyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyz_0[i] = g_0_xxxyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyz_0[i] * pb_z + g_0_xxxyyyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzz_0[i] * pb_z + g_0_xxxyyyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xzzz_0[i] * pb_z + g_0_xxxyyyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyy_0[i] = g_0_xxxyyyy_0_yyyy_0[i] * pb_z + g_0_xxxyyyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyz_0[i] = g_0_xxxyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyz_0[i] * pb_z + g_0_xxxyyyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyzz_0[i] = 2.0 * g_0_xxxyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyzz_0[i] * pb_z + g_0_xxxyyyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yzzz_0[i] = 3.0 * g_0_xxxyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yzzz_0[i] * pb_z + g_0_xxxyyyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_zzzz_0[i] = 4.0 * g_0_xxxyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_zzzz_0[i] * pb_z + g_0_xxxyyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 255-270 components of targeted buffer : SLSG

    auto g_0_xxxyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 255);

    auto g_0_xxxyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 256);

    auto g_0_xxxyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 257);

    auto g_0_xxxyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 258);

    auto g_0_xxxyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 259);

    auto g_0_xxxyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 260);

    auto g_0_xxxyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 261);

    auto g_0_xxxyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 262);

    auto g_0_xxxyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 263);

    auto g_0_xxxyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 264);

    auto g_0_xxxyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 265);

    auto g_0_xxxyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 266);

    auto g_0_xxxyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 267);

    auto g_0_xxxyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 268);

    auto g_0_xxxyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 269);

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxy_0, g_0_xxxyyy_0_xxxy_1, g_0_xxxyyy_0_xxyy_0, g_0_xxxyyy_0_xxyy_1, g_0_xxxyyy_0_xyyy_0, g_0_xxxyyy_0_xyyy_1, g_0_xxxyyyz_0_xxxy_0, g_0_xxxyyyz_0_xxxy_1, g_0_xxxyyyz_0_xxyy_0, g_0_xxxyyyz_0_xxyy_1, g_0_xxxyyyz_0_xyyy_0, g_0_xxxyyyz_0_xyyy_1, g_0_xxxyyyzz_0_xxxx_0, g_0_xxxyyyzz_0_xxxy_0, g_0_xxxyyyzz_0_xxxz_0, g_0_xxxyyyzz_0_xxyy_0, g_0_xxxyyyzz_0_xxyz_0, g_0_xxxyyyzz_0_xxzz_0, g_0_xxxyyyzz_0_xyyy_0, g_0_xxxyyyzz_0_xyyz_0, g_0_xxxyyyzz_0_xyzz_0, g_0_xxxyyyzz_0_xzzz_0, g_0_xxxyyyzz_0_yyyy_0, g_0_xxxyyyzz_0_yyyz_0, g_0_xxxyyyzz_0_yyzz_0, g_0_xxxyyyzz_0_yzzz_0, g_0_xxxyyyzz_0_zzzz_0, g_0_xxxyyzz_0_xxxx_0, g_0_xxxyyzz_0_xxxx_1, g_0_xxxyyzz_0_xxxz_0, g_0_xxxyyzz_0_xxxz_1, g_0_xxxyyzz_0_xxzz_0, g_0_xxxyyzz_0_xxzz_1, g_0_xxxyyzz_0_xzzz_0, g_0_xxxyyzz_0_xzzz_1, g_0_xxxyzz_0_xxxx_0, g_0_xxxyzz_0_xxxx_1, g_0_xxxyzz_0_xxxz_0, g_0_xxxyzz_0_xxxz_1, g_0_xxxyzz_0_xxzz_0, g_0_xxxyzz_0_xxzz_1, g_0_xxxyzz_0_xzzz_0, g_0_xxxyzz_0_xzzz_1, g_0_xxyyyzz_0_xxyz_0, g_0_xxyyyzz_0_xxyz_1, g_0_xxyyyzz_0_xyyz_0, g_0_xxyyyzz_0_xyyz_1, g_0_xxyyyzz_0_xyz_1, g_0_xxyyyzz_0_xyzz_0, g_0_xxyyyzz_0_xyzz_1, g_0_xxyyyzz_0_yyyy_0, g_0_xxyyyzz_0_yyyy_1, g_0_xxyyyzz_0_yyyz_0, g_0_xxyyyzz_0_yyyz_1, g_0_xxyyyzz_0_yyz_1, g_0_xxyyyzz_0_yyzz_0, g_0_xxyyyzz_0_yyzz_1, g_0_xxyyyzz_0_yzz_1, g_0_xxyyyzz_0_yzzz_0, g_0_xxyyyzz_0_yzzz_1, g_0_xxyyyzz_0_zzzz_0, g_0_xxyyyzz_0_zzzz_1, g_0_xyyyzz_0_xxyz_0, g_0_xyyyzz_0_xxyz_1, g_0_xyyyzz_0_xyyz_0, g_0_xyyyzz_0_xyyz_1, g_0_xyyyzz_0_xyzz_0, g_0_xyyyzz_0_xyzz_1, g_0_xyyyzz_0_yyyy_0, g_0_xyyyzz_0_yyyy_1, g_0_xyyyzz_0_yyyz_0, g_0_xyyyzz_0_yyyz_1, g_0_xyyyzz_0_yyzz_0, g_0_xyyyzz_0_yyzz_1, g_0_xyyyzz_0_yzzz_0, g_0_xyyyzz_0_yzzz_1, g_0_xyyyzz_0_zzzz_0, g_0_xyyyzz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_xxxx_0[i] = 2.0 * g_0_xxxyzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxx_0[i] * pb_y + g_0_xxxyyzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxy_0[i] = g_0_xxxyyy_0_xxxy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxy_0[i] * pb_z + g_0_xxxyyyz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxz_0[i] = 2.0 * g_0_xxxyzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxz_0[i] * pb_y + g_0_xxxyyzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxyy_0[i] = g_0_xxxyyy_0_xxyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxyy_0[i] * pb_z + g_0_xxxyyyz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxyz_0[i] = 2.0 * g_0_xyyyzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyz_0[i] * pb_x + g_0_xxyyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxzz_0[i] = 2.0 * g_0_xxxyzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxzz_0[i] * pb_y + g_0_xxxyyzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xyyy_0[i] = g_0_xxxyyy_0_xyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xyyy_0[i] * pb_z + g_0_xxxyyyz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xyyz_0[i] = 2.0 * g_0_xyyyzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyz_0[i] * pb_x + g_0_xxyyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyzz_0[i] = 2.0 * g_0_xyyyzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyzz_0[i] * pb_x + g_0_xxyyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xzzz_0[i] = 2.0 * g_0_xxxyzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xzzz_0[i] * pb_y + g_0_xxxyyzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_yyyy_0[i] = 2.0 * g_0_xyyyzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyy_0[i] * pb_x + g_0_xxyyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyz_0[i] = 2.0 * g_0_xyyyzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyz_0[i] * pb_x + g_0_xxyyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyzz_0[i] = 2.0 * g_0_xyyyzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyzz_0[i] * pb_x + g_0_xxyyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yzzz_0[i] = 2.0 * g_0_xyyyzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzzz_0[i] * pb_x + g_0_xxyyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_zzzz_0[i] = 2.0 * g_0_xyyyzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_zzzz_0[i] * pb_x + g_0_xxyyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 270-285 components of targeted buffer : SLSG

    auto g_0_xxxyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 270);

    auto g_0_xxxyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 271);

    auto g_0_xxxyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 272);

    auto g_0_xxxyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 273);

    auto g_0_xxxyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 274);

    auto g_0_xxxyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 275);

    auto g_0_xxxyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 276);

    auto g_0_xxxyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 277);

    auto g_0_xxxyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 278);

    auto g_0_xxxyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 279);

    auto g_0_xxxyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 280);

    auto g_0_xxxyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 281);

    auto g_0_xxxyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 282);

    auto g_0_xxxyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 283);

    auto g_0_xxxyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 284);

    #pragma omp simd aligned(g_0_xxxyyz_0_xxxy_0, g_0_xxxyyz_0_xxxy_1, g_0_xxxyyz_0_xxyy_0, g_0_xxxyyz_0_xxyy_1, g_0_xxxyyz_0_xyyy_0, g_0_xxxyyz_0_xyyy_1, g_0_xxxyyzz_0_xxxy_0, g_0_xxxyyzz_0_xxxy_1, g_0_xxxyyzz_0_xxyy_0, g_0_xxxyyzz_0_xxyy_1, g_0_xxxyyzz_0_xyyy_0, g_0_xxxyyzz_0_xyyy_1, g_0_xxxyyzzz_0_xxxx_0, g_0_xxxyyzzz_0_xxxy_0, g_0_xxxyyzzz_0_xxxz_0, g_0_xxxyyzzz_0_xxyy_0, g_0_xxxyyzzz_0_xxyz_0, g_0_xxxyyzzz_0_xxzz_0, g_0_xxxyyzzz_0_xyyy_0, g_0_xxxyyzzz_0_xyyz_0, g_0_xxxyyzzz_0_xyzz_0, g_0_xxxyyzzz_0_xzzz_0, g_0_xxxyyzzz_0_yyyy_0, g_0_xxxyyzzz_0_yyyz_0, g_0_xxxyyzzz_0_yyzz_0, g_0_xxxyyzzz_0_yzzz_0, g_0_xxxyyzzz_0_zzzz_0, g_0_xxxyzzz_0_xxxx_0, g_0_xxxyzzz_0_xxxx_1, g_0_xxxyzzz_0_xxxz_0, g_0_xxxyzzz_0_xxxz_1, g_0_xxxyzzz_0_xxzz_0, g_0_xxxyzzz_0_xxzz_1, g_0_xxxyzzz_0_xzzz_0, g_0_xxxyzzz_0_xzzz_1, g_0_xxxzzz_0_xxxx_0, g_0_xxxzzz_0_xxxx_1, g_0_xxxzzz_0_xxxz_0, g_0_xxxzzz_0_xxxz_1, g_0_xxxzzz_0_xxzz_0, g_0_xxxzzz_0_xxzz_1, g_0_xxxzzz_0_xzzz_0, g_0_xxxzzz_0_xzzz_1, g_0_xxyyzzz_0_xxyz_0, g_0_xxyyzzz_0_xxyz_1, g_0_xxyyzzz_0_xyyz_0, g_0_xxyyzzz_0_xyyz_1, g_0_xxyyzzz_0_xyz_1, g_0_xxyyzzz_0_xyzz_0, g_0_xxyyzzz_0_xyzz_1, g_0_xxyyzzz_0_yyyy_0, g_0_xxyyzzz_0_yyyy_1, g_0_xxyyzzz_0_yyyz_0, g_0_xxyyzzz_0_yyyz_1, g_0_xxyyzzz_0_yyz_1, g_0_xxyyzzz_0_yyzz_0, g_0_xxyyzzz_0_yyzz_1, g_0_xxyyzzz_0_yzz_1, g_0_xxyyzzz_0_yzzz_0, g_0_xxyyzzz_0_yzzz_1, g_0_xxyyzzz_0_zzzz_0, g_0_xxyyzzz_0_zzzz_1, g_0_xyyzzz_0_xxyz_0, g_0_xyyzzz_0_xxyz_1, g_0_xyyzzz_0_xyyz_0, g_0_xyyzzz_0_xyyz_1, g_0_xyyzzz_0_xyzz_0, g_0_xyyzzz_0_xyzz_1, g_0_xyyzzz_0_yyyy_0, g_0_xyyzzz_0_yyyy_1, g_0_xyyzzz_0_yyyz_0, g_0_xyyzzz_0_yyyz_1, g_0_xyyzzz_0_yyzz_0, g_0_xyyzzz_0_yyzz_1, g_0_xyyzzz_0_yzzz_0, g_0_xyyzzz_0_yzzz_1, g_0_xyyzzz_0_zzzz_0, g_0_xyyzzz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_xxxx_0[i] = g_0_xxxzzz_0_xxxx_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxx_0[i] * pb_y + g_0_xxxyzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxy_0[i] = 2.0 * g_0_xxxyyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxy_0[i] * pb_z + g_0_xxxyyzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxz_0[i] = g_0_xxxzzz_0_xxxz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxz_0[i] * pb_y + g_0_xxxyzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxyy_0[i] = 2.0 * g_0_xxxyyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxyy_0[i] * pb_z + g_0_xxxyyzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxyz_0[i] = 2.0 * g_0_xyyzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyz_0[i] * pb_x + g_0_xxyyzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxzz_0[i] = g_0_xxxzzz_0_xxzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxzz_0[i] * pb_y + g_0_xxxyzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xyyy_0[i] = 2.0 * g_0_xxxyyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xyyy_0[i] * pb_z + g_0_xxxyyzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xyyz_0[i] = 2.0 * g_0_xyyzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyz_0[i] * pb_x + g_0_xxyyzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyzz_0[i] = 2.0 * g_0_xyyzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyzz_0[i] * pb_x + g_0_xxyyzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xzzz_0[i] = g_0_xxxzzz_0_xzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xzzz_0[i] * pb_y + g_0_xxxyzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_yyyy_0[i] = 2.0 * g_0_xyyzzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyy_0[i] * pb_x + g_0_xxyyzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyz_0[i] = 2.0 * g_0_xyyzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyz_0[i] * pb_x + g_0_xxyyzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyzz_0[i] = 2.0 * g_0_xyyzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyzz_0[i] * pb_x + g_0_xxyyzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yzzz_0[i] = 2.0 * g_0_xyyzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzzz_0[i] * pb_x + g_0_xxyyzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_zzzz_0[i] = 2.0 * g_0_xyyzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_zzzz_0[i] * pb_x + g_0_xxyyzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 285-300 components of targeted buffer : SLSG

    auto g_0_xxxyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 285);

    auto g_0_xxxyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 286);

    auto g_0_xxxyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 287);

    auto g_0_xxxyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 288);

    auto g_0_xxxyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 289);

    auto g_0_xxxyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 290);

    auto g_0_xxxyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 291);

    auto g_0_xxxyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 292);

    auto g_0_xxxyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 293);

    auto g_0_xxxyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 294);

    auto g_0_xxxyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 295);

    auto g_0_xxxyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 296);

    auto g_0_xxxyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 297);

    auto g_0_xxxyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 298);

    auto g_0_xxxyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 299);

    #pragma omp simd aligned(g_0_xxxyzzzz_0_xxxx_0, g_0_xxxyzzzz_0_xxxy_0, g_0_xxxyzzzz_0_xxxz_0, g_0_xxxyzzzz_0_xxyy_0, g_0_xxxyzzzz_0_xxyz_0, g_0_xxxyzzzz_0_xxzz_0, g_0_xxxyzzzz_0_xyyy_0, g_0_xxxyzzzz_0_xyyz_0, g_0_xxxyzzzz_0_xyzz_0, g_0_xxxyzzzz_0_xzzz_0, g_0_xxxyzzzz_0_yyyy_0, g_0_xxxyzzzz_0_yyyz_0, g_0_xxxyzzzz_0_yyzz_0, g_0_xxxyzzzz_0_yzzz_0, g_0_xxxyzzzz_0_zzzz_0, g_0_xxxzzzz_0_xxx_1, g_0_xxxzzzz_0_xxxx_0, g_0_xxxzzzz_0_xxxx_1, g_0_xxxzzzz_0_xxxy_0, g_0_xxxzzzz_0_xxxy_1, g_0_xxxzzzz_0_xxxz_0, g_0_xxxzzzz_0_xxxz_1, g_0_xxxzzzz_0_xxy_1, g_0_xxxzzzz_0_xxyy_0, g_0_xxxzzzz_0_xxyy_1, g_0_xxxzzzz_0_xxyz_0, g_0_xxxzzzz_0_xxyz_1, g_0_xxxzzzz_0_xxz_1, g_0_xxxzzzz_0_xxzz_0, g_0_xxxzzzz_0_xxzz_1, g_0_xxxzzzz_0_xyy_1, g_0_xxxzzzz_0_xyyy_0, g_0_xxxzzzz_0_xyyy_1, g_0_xxxzzzz_0_xyyz_0, g_0_xxxzzzz_0_xyyz_1, g_0_xxxzzzz_0_xyz_1, g_0_xxxzzzz_0_xyzz_0, g_0_xxxzzzz_0_xyzz_1, g_0_xxxzzzz_0_xzz_1, g_0_xxxzzzz_0_xzzz_0, g_0_xxxzzzz_0_xzzz_1, g_0_xxxzzzz_0_yyy_1, g_0_xxxzzzz_0_yyyy_0, g_0_xxxzzzz_0_yyyy_1, g_0_xxxzzzz_0_yyyz_0, g_0_xxxzzzz_0_yyyz_1, g_0_xxxzzzz_0_yyz_1, g_0_xxxzzzz_0_yyzz_0, g_0_xxxzzzz_0_yyzz_1, g_0_xxxzzzz_0_yzz_1, g_0_xxxzzzz_0_yzzz_0, g_0_xxxzzzz_0_yzzz_1, g_0_xxxzzzz_0_zzz_1, g_0_xxxzzzz_0_zzzz_0, g_0_xxxzzzz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_xxxx_0[i] = g_0_xxxzzzz_0_xxxx_0[i] * pb_y + g_0_xxxzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxy_0[i] = g_0_xxxzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxy_0[i] * pb_y + g_0_xxxzzzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxz_0[i] = g_0_xxxzzzz_0_xxxz_0[i] * pb_y + g_0_xxxzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyy_0[i] = 2.0 * g_0_xxxzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyy_0[i] * pb_y + g_0_xxxzzzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyz_0[i] = g_0_xxxzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyz_0[i] * pb_y + g_0_xxxzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxzz_0[i] = g_0_xxxzzzz_0_xxzz_0[i] * pb_y + g_0_xxxzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyy_0[i] = 3.0 * g_0_xxxzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyy_0[i] * pb_y + g_0_xxxzzzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyz_0[i] = 2.0 * g_0_xxxzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyz_0[i] * pb_y + g_0_xxxzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyzz_0[i] = g_0_xxxzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzz_0[i] * pb_y + g_0_xxxzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xzzz_0[i] = g_0_xxxzzzz_0_xzzz_0[i] * pb_y + g_0_xxxzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyy_0[i] = 4.0 * g_0_xxxzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyy_0[i] * pb_y + g_0_xxxzzzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyz_0[i] = 3.0 * g_0_xxxzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyz_0[i] * pb_y + g_0_xxxzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyzz_0[i] = 2.0 * g_0_xxxzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyzz_0[i] * pb_y + g_0_xxxzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yzzz_0[i] = g_0_xxxzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yzzz_0[i] * pb_y + g_0_xxxzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_zzzz_0[i] = g_0_xxxzzzz_0_zzzz_0[i] * pb_y + g_0_xxxzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 300-315 components of targeted buffer : SLSG

    auto g_0_xxxzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 300);

    auto g_0_xxxzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 301);

    auto g_0_xxxzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 302);

    auto g_0_xxxzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 303);

    auto g_0_xxxzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 304);

    auto g_0_xxxzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 305);

    auto g_0_xxxzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 306);

    auto g_0_xxxzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 307);

    auto g_0_xxxzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 308);

    auto g_0_xxxzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 309);

    auto g_0_xxxzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 310);

    auto g_0_xxxzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 311);

    auto g_0_xxxzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 312);

    auto g_0_xxxzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 313);

    auto g_0_xxxzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 314);

    #pragma omp simd aligned(g_0_xxxzzz_0_xxxx_0, g_0_xxxzzz_0_xxxx_1, g_0_xxxzzz_0_xxxy_0, g_0_xxxzzz_0_xxxy_1, g_0_xxxzzz_0_xxyy_0, g_0_xxxzzz_0_xxyy_1, g_0_xxxzzz_0_xyyy_0, g_0_xxxzzz_0_xyyy_1, g_0_xxxzzzz_0_xxxx_0, g_0_xxxzzzz_0_xxxx_1, g_0_xxxzzzz_0_xxxy_0, g_0_xxxzzzz_0_xxxy_1, g_0_xxxzzzz_0_xxyy_0, g_0_xxxzzzz_0_xxyy_1, g_0_xxxzzzz_0_xyyy_0, g_0_xxxzzzz_0_xyyy_1, g_0_xxxzzzzz_0_xxxx_0, g_0_xxxzzzzz_0_xxxy_0, g_0_xxxzzzzz_0_xxxz_0, g_0_xxxzzzzz_0_xxyy_0, g_0_xxxzzzzz_0_xxyz_0, g_0_xxxzzzzz_0_xxzz_0, g_0_xxxzzzzz_0_xyyy_0, g_0_xxxzzzzz_0_xyyz_0, g_0_xxxzzzzz_0_xyzz_0, g_0_xxxzzzzz_0_xzzz_0, g_0_xxxzzzzz_0_yyyy_0, g_0_xxxzzzzz_0_yyyz_0, g_0_xxxzzzzz_0_yyzz_0, g_0_xxxzzzzz_0_yzzz_0, g_0_xxxzzzzz_0_zzzz_0, g_0_xxzzzzz_0_xxxz_0, g_0_xxzzzzz_0_xxxz_1, g_0_xxzzzzz_0_xxyz_0, g_0_xxzzzzz_0_xxyz_1, g_0_xxzzzzz_0_xxz_1, g_0_xxzzzzz_0_xxzz_0, g_0_xxzzzzz_0_xxzz_1, g_0_xxzzzzz_0_xyyz_0, g_0_xxzzzzz_0_xyyz_1, g_0_xxzzzzz_0_xyz_1, g_0_xxzzzzz_0_xyzz_0, g_0_xxzzzzz_0_xyzz_1, g_0_xxzzzzz_0_xzz_1, g_0_xxzzzzz_0_xzzz_0, g_0_xxzzzzz_0_xzzz_1, g_0_xxzzzzz_0_yyyy_0, g_0_xxzzzzz_0_yyyy_1, g_0_xxzzzzz_0_yyyz_0, g_0_xxzzzzz_0_yyyz_1, g_0_xxzzzzz_0_yyz_1, g_0_xxzzzzz_0_yyzz_0, g_0_xxzzzzz_0_yyzz_1, g_0_xxzzzzz_0_yzz_1, g_0_xxzzzzz_0_yzzz_0, g_0_xxzzzzz_0_yzzz_1, g_0_xxzzzzz_0_zzz_1, g_0_xxzzzzz_0_zzzz_0, g_0_xxzzzzz_0_zzzz_1, g_0_xzzzzz_0_xxxz_0, g_0_xzzzzz_0_xxxz_1, g_0_xzzzzz_0_xxyz_0, g_0_xzzzzz_0_xxyz_1, g_0_xzzzzz_0_xxzz_0, g_0_xzzzzz_0_xxzz_1, g_0_xzzzzz_0_xyyz_0, g_0_xzzzzz_0_xyyz_1, g_0_xzzzzz_0_xyzz_0, g_0_xzzzzz_0_xyzz_1, g_0_xzzzzz_0_xzzz_0, g_0_xzzzzz_0_xzzz_1, g_0_xzzzzz_0_yyyy_0, g_0_xzzzzz_0_yyyy_1, g_0_xzzzzz_0_yyyz_0, g_0_xzzzzz_0_yyyz_1, g_0_xzzzzz_0_yyzz_0, g_0_xzzzzz_0_yyzz_1, g_0_xzzzzz_0_yzzz_0, g_0_xzzzzz_0_yzzz_1, g_0_xzzzzz_0_zzzz_0, g_0_xzzzzz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_xxxx_0[i] = 4.0 * g_0_xxxzzz_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxx_0[i] * pb_z + g_0_xxxzzzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxy_0[i] = 4.0 * g_0_xxxzzz_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxy_0[i] * pb_z + g_0_xxxzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxz_0[i] = 2.0 * g_0_xzzzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxz_0[i] * pb_x + g_0_xxzzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyy_0[i] = 4.0 * g_0_xxxzzz_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxyy_0[i] * pb_z + g_0_xxxzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxyz_0[i] = 2.0 * g_0_xzzzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyz_0[i] * pb_x + g_0_xxzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxzz_0[i] = 2.0 * g_0_xzzzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxzz_0[i] * pb_x + g_0_xxzzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyy_0[i] = 4.0 * g_0_xxxzzz_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xyyy_0[i] * pb_z + g_0_xxxzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xyyz_0[i] = 2.0 * g_0_xzzzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyz_0[i] * pb_x + g_0_xxzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyzz_0[i] = 2.0 * g_0_xzzzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzz_0[i] * pb_x + g_0_xxzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xzzz_0[i] = 2.0 * g_0_xzzzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xzzz_0[i] * pb_x + g_0_xxzzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyy_0[i] = 2.0 * g_0_xzzzzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyy_0[i] * pb_x + g_0_xxzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyz_0[i] = 2.0 * g_0_xzzzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyz_0[i] * pb_x + g_0_xxzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyzz_0[i] = 2.0 * g_0_xzzzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyzz_0[i] * pb_x + g_0_xxzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yzzz_0[i] = 2.0 * g_0_xzzzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzzz_0[i] * pb_x + g_0_xxzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_zzzz_0[i] = 2.0 * g_0_xzzzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzzz_0[i] * pb_x + g_0_xxzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 315-330 components of targeted buffer : SLSG

    auto g_0_xxyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 315);

    auto g_0_xxyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 316);

    auto g_0_xxyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 317);

    auto g_0_xxyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 318);

    auto g_0_xxyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 319);

    auto g_0_xxyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 320);

    auto g_0_xxyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 321);

    auto g_0_xxyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 322);

    auto g_0_xxyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 323);

    auto g_0_xxyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 324);

    auto g_0_xxyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 325);

    auto g_0_xxyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 326);

    auto g_0_xxyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 327);

    auto g_0_xxyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 328);

    auto g_0_xxyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 329);

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxx_0, g_0_xxyyyy_0_xxxx_1, g_0_xxyyyy_0_xxxz_0, g_0_xxyyyy_0_xxxz_1, g_0_xxyyyy_0_xxzz_0, g_0_xxyyyy_0_xxzz_1, g_0_xxyyyy_0_xzzz_0, g_0_xxyyyy_0_xzzz_1, g_0_xxyyyyy_0_xxxx_0, g_0_xxyyyyy_0_xxxx_1, g_0_xxyyyyy_0_xxxz_0, g_0_xxyyyyy_0_xxxz_1, g_0_xxyyyyy_0_xxzz_0, g_0_xxyyyyy_0_xxzz_1, g_0_xxyyyyy_0_xzzz_0, g_0_xxyyyyy_0_xzzz_1, g_0_xxyyyyyy_0_xxxx_0, g_0_xxyyyyyy_0_xxxy_0, g_0_xxyyyyyy_0_xxxz_0, g_0_xxyyyyyy_0_xxyy_0, g_0_xxyyyyyy_0_xxyz_0, g_0_xxyyyyyy_0_xxzz_0, g_0_xxyyyyyy_0_xyyy_0, g_0_xxyyyyyy_0_xyyz_0, g_0_xxyyyyyy_0_xyzz_0, g_0_xxyyyyyy_0_xzzz_0, g_0_xxyyyyyy_0_yyyy_0, g_0_xxyyyyyy_0_yyyz_0, g_0_xxyyyyyy_0_yyzz_0, g_0_xxyyyyyy_0_yzzz_0, g_0_xxyyyyyy_0_zzzz_0, g_0_xyyyyyy_0_xxxy_0, g_0_xyyyyyy_0_xxxy_1, g_0_xyyyyyy_0_xxy_1, g_0_xyyyyyy_0_xxyy_0, g_0_xyyyyyy_0_xxyy_1, g_0_xyyyyyy_0_xxyz_0, g_0_xyyyyyy_0_xxyz_1, g_0_xyyyyyy_0_xyy_1, g_0_xyyyyyy_0_xyyy_0, g_0_xyyyyyy_0_xyyy_1, g_0_xyyyyyy_0_xyyz_0, g_0_xyyyyyy_0_xyyz_1, g_0_xyyyyyy_0_xyz_1, g_0_xyyyyyy_0_xyzz_0, g_0_xyyyyyy_0_xyzz_1, g_0_xyyyyyy_0_yyy_1, g_0_xyyyyyy_0_yyyy_0, g_0_xyyyyyy_0_yyyy_1, g_0_xyyyyyy_0_yyyz_0, g_0_xyyyyyy_0_yyyz_1, g_0_xyyyyyy_0_yyz_1, g_0_xyyyyyy_0_yyzz_0, g_0_xyyyyyy_0_yyzz_1, g_0_xyyyyyy_0_yzz_1, g_0_xyyyyyy_0_yzzz_0, g_0_xyyyyyy_0_yzzz_1, g_0_xyyyyyy_0_zzzz_0, g_0_xyyyyyy_0_zzzz_1, g_0_yyyyyy_0_xxxy_0, g_0_yyyyyy_0_xxxy_1, g_0_yyyyyy_0_xxyy_0, g_0_yyyyyy_0_xxyy_1, g_0_yyyyyy_0_xxyz_0, g_0_yyyyyy_0_xxyz_1, g_0_yyyyyy_0_xyyy_0, g_0_yyyyyy_0_xyyy_1, g_0_yyyyyy_0_xyyz_0, g_0_yyyyyy_0_xyyz_1, g_0_yyyyyy_0_xyzz_0, g_0_yyyyyy_0_xyzz_1, g_0_yyyyyy_0_yyyy_0, g_0_yyyyyy_0_yyyy_1, g_0_yyyyyy_0_yyyz_0, g_0_yyyyyy_0_yyyz_1, g_0_yyyyyy_0_yyzz_0, g_0_yyyyyy_0_yyzz_1, g_0_yyyyyy_0_yzzz_0, g_0_yyyyyy_0_yzzz_1, g_0_yyyyyy_0_zzzz_0, g_0_yyyyyy_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_xxxx_0[i] = 5.0 * g_0_xxyyyy_0_xxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxx_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxx_0[i] * pb_y + g_0_xxyyyyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxy_0[i] = g_0_yyyyyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxy_0[i] * pb_x + g_0_xyyyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxz_0[i] = 5.0 * g_0_xxyyyy_0_xxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxz_0[i] * pb_y + g_0_xxyyyyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxyy_0[i] = g_0_yyyyyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyy_0[i] * pb_x + g_0_xyyyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyz_0[i] = g_0_yyyyyy_0_xxyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyz_0[i] * pb_x + g_0_xyyyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxzz_0[i] = 5.0 * g_0_xxyyyy_0_xxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxzz_0[i] * pb_y + g_0_xxyyyyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xyyy_0[i] = g_0_yyyyyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyy_0[i] * pb_x + g_0_xyyyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyz_0[i] = g_0_yyyyyy_0_xyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyz_0[i] * pb_x + g_0_xyyyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyzz_0[i] = g_0_yyyyyy_0_xyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyzz_0[i] * pb_x + g_0_xyyyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xzzz_0[i] = 5.0 * g_0_xxyyyy_0_xzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xzzz_0[i] * pb_y + g_0_xxyyyyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_yyyy_0[i] = g_0_yyyyyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyy_0[i] * pb_x + g_0_xyyyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyz_0[i] = g_0_yyyyyy_0_yyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyz_0[i] * pb_x + g_0_xyyyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyzz_0[i] = g_0_yyyyyy_0_yyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzz_0[i] * pb_x + g_0_xyyyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yzzz_0[i] = g_0_yyyyyy_0_yzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzz_0[i] * pb_x + g_0_xyyyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_zzzz_0[i] = g_0_yyyyyy_0_zzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_zzzz_0[i] * pb_x + g_0_xyyyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 330-345 components of targeted buffer : SLSG

    auto g_0_xxyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 330);

    auto g_0_xxyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 331);

    auto g_0_xxyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 332);

    auto g_0_xxyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 333);

    auto g_0_xxyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 334);

    auto g_0_xxyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 335);

    auto g_0_xxyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 336);

    auto g_0_xxyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 337);

    auto g_0_xxyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 338);

    auto g_0_xxyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 339);

    auto g_0_xxyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 340);

    auto g_0_xxyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 341);

    auto g_0_xxyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 342);

    auto g_0_xxyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 343);

    auto g_0_xxyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 344);

    #pragma omp simd aligned(g_0_xxyyyyy_0_xxx_1, g_0_xxyyyyy_0_xxxx_0, g_0_xxyyyyy_0_xxxx_1, g_0_xxyyyyy_0_xxxy_0, g_0_xxyyyyy_0_xxxy_1, g_0_xxyyyyy_0_xxxz_0, g_0_xxyyyyy_0_xxxz_1, g_0_xxyyyyy_0_xxy_1, g_0_xxyyyyy_0_xxyy_0, g_0_xxyyyyy_0_xxyy_1, g_0_xxyyyyy_0_xxyz_0, g_0_xxyyyyy_0_xxyz_1, g_0_xxyyyyy_0_xxz_1, g_0_xxyyyyy_0_xxzz_0, g_0_xxyyyyy_0_xxzz_1, g_0_xxyyyyy_0_xyy_1, g_0_xxyyyyy_0_xyyy_0, g_0_xxyyyyy_0_xyyy_1, g_0_xxyyyyy_0_xyyz_0, g_0_xxyyyyy_0_xyyz_1, g_0_xxyyyyy_0_xyz_1, g_0_xxyyyyy_0_xyzz_0, g_0_xxyyyyy_0_xyzz_1, g_0_xxyyyyy_0_xzz_1, g_0_xxyyyyy_0_xzzz_0, g_0_xxyyyyy_0_xzzz_1, g_0_xxyyyyy_0_yyy_1, g_0_xxyyyyy_0_yyyy_0, g_0_xxyyyyy_0_yyyy_1, g_0_xxyyyyy_0_yyyz_0, g_0_xxyyyyy_0_yyyz_1, g_0_xxyyyyy_0_yyz_1, g_0_xxyyyyy_0_yyzz_0, g_0_xxyyyyy_0_yyzz_1, g_0_xxyyyyy_0_yzz_1, g_0_xxyyyyy_0_yzzz_0, g_0_xxyyyyy_0_yzzz_1, g_0_xxyyyyy_0_zzz_1, g_0_xxyyyyy_0_zzzz_0, g_0_xxyyyyy_0_zzzz_1, g_0_xxyyyyyz_0_xxxx_0, g_0_xxyyyyyz_0_xxxy_0, g_0_xxyyyyyz_0_xxxz_0, g_0_xxyyyyyz_0_xxyy_0, g_0_xxyyyyyz_0_xxyz_0, g_0_xxyyyyyz_0_xxzz_0, g_0_xxyyyyyz_0_xyyy_0, g_0_xxyyyyyz_0_xyyz_0, g_0_xxyyyyyz_0_xyzz_0, g_0_xxyyyyyz_0_xzzz_0, g_0_xxyyyyyz_0_yyyy_0, g_0_xxyyyyyz_0_yyyz_0, g_0_xxyyyyyz_0_yyzz_0, g_0_xxyyyyyz_0_yzzz_0, g_0_xxyyyyyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_xxxx_0[i] = g_0_xxyyyyy_0_xxxx_0[i] * pb_z + g_0_xxyyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxy_0[i] = g_0_xxyyyyy_0_xxxy_0[i] * pb_z + g_0_xxyyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxz_0[i] = g_0_xxyyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxz_0[i] * pb_z + g_0_xxyyyyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyy_0[i] = g_0_xxyyyyy_0_xxyy_0[i] * pb_z + g_0_xxyyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyz_0[i] = g_0_xxyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyz_0[i] * pb_z + g_0_xxyyyyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxzz_0[i] * pb_z + g_0_xxyyyyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyy_0[i] = g_0_xxyyyyy_0_xyyy_0[i] * pb_z + g_0_xxyyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyz_0[i] = g_0_xxyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyz_0[i] * pb_z + g_0_xxyyyyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzz_0[i] * pb_z + g_0_xxyyyyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xzzz_0[i] * pb_z + g_0_xxyyyyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyy_0[i] = g_0_xxyyyyy_0_yyyy_0[i] * pb_z + g_0_xxyyyyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyz_0[i] = g_0_xxyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyz_0[i] * pb_z + g_0_xxyyyyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyzz_0[i] = 2.0 * g_0_xxyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyzz_0[i] * pb_z + g_0_xxyyyyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yzzz_0[i] = 3.0 * g_0_xxyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yzzz_0[i] * pb_z + g_0_xxyyyyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_zzzz_0[i] = 4.0 * g_0_xxyyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_zzzz_0[i] * pb_z + g_0_xxyyyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 345-360 components of targeted buffer : SLSG

    auto g_0_xxyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 345);

    auto g_0_xxyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 346);

    auto g_0_xxyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 347);

    auto g_0_xxyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 348);

    auto g_0_xxyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 349);

    auto g_0_xxyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 350);

    auto g_0_xxyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 351);

    auto g_0_xxyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 352);

    auto g_0_xxyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 353);

    auto g_0_xxyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 354);

    auto g_0_xxyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 355);

    auto g_0_xxyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 356);

    auto g_0_xxyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 357);

    auto g_0_xxyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 358);

    auto g_0_xxyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 359);

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxy_0, g_0_xxyyyy_0_xxxy_1, g_0_xxyyyy_0_xxyy_0, g_0_xxyyyy_0_xxyy_1, g_0_xxyyyy_0_xyyy_0, g_0_xxyyyy_0_xyyy_1, g_0_xxyyyyz_0_xxxy_0, g_0_xxyyyyz_0_xxxy_1, g_0_xxyyyyz_0_xxyy_0, g_0_xxyyyyz_0_xxyy_1, g_0_xxyyyyz_0_xyyy_0, g_0_xxyyyyz_0_xyyy_1, g_0_xxyyyyzz_0_xxxx_0, g_0_xxyyyyzz_0_xxxy_0, g_0_xxyyyyzz_0_xxxz_0, g_0_xxyyyyzz_0_xxyy_0, g_0_xxyyyyzz_0_xxyz_0, g_0_xxyyyyzz_0_xxzz_0, g_0_xxyyyyzz_0_xyyy_0, g_0_xxyyyyzz_0_xyyz_0, g_0_xxyyyyzz_0_xyzz_0, g_0_xxyyyyzz_0_xzzz_0, g_0_xxyyyyzz_0_yyyy_0, g_0_xxyyyyzz_0_yyyz_0, g_0_xxyyyyzz_0_yyzz_0, g_0_xxyyyyzz_0_yzzz_0, g_0_xxyyyyzz_0_zzzz_0, g_0_xxyyyzz_0_xxxx_0, g_0_xxyyyzz_0_xxxx_1, g_0_xxyyyzz_0_xxxz_0, g_0_xxyyyzz_0_xxxz_1, g_0_xxyyyzz_0_xxzz_0, g_0_xxyyyzz_0_xxzz_1, g_0_xxyyyzz_0_xzzz_0, g_0_xxyyyzz_0_xzzz_1, g_0_xxyyzz_0_xxxx_0, g_0_xxyyzz_0_xxxx_1, g_0_xxyyzz_0_xxxz_0, g_0_xxyyzz_0_xxxz_1, g_0_xxyyzz_0_xxzz_0, g_0_xxyyzz_0_xxzz_1, g_0_xxyyzz_0_xzzz_0, g_0_xxyyzz_0_xzzz_1, g_0_xyyyyzz_0_xxyz_0, g_0_xyyyyzz_0_xxyz_1, g_0_xyyyyzz_0_xyyz_0, g_0_xyyyyzz_0_xyyz_1, g_0_xyyyyzz_0_xyz_1, g_0_xyyyyzz_0_xyzz_0, g_0_xyyyyzz_0_xyzz_1, g_0_xyyyyzz_0_yyyy_0, g_0_xyyyyzz_0_yyyy_1, g_0_xyyyyzz_0_yyyz_0, g_0_xyyyyzz_0_yyyz_1, g_0_xyyyyzz_0_yyz_1, g_0_xyyyyzz_0_yyzz_0, g_0_xyyyyzz_0_yyzz_1, g_0_xyyyyzz_0_yzz_1, g_0_xyyyyzz_0_yzzz_0, g_0_xyyyyzz_0_yzzz_1, g_0_xyyyyzz_0_zzzz_0, g_0_xyyyyzz_0_zzzz_1, g_0_yyyyzz_0_xxyz_0, g_0_yyyyzz_0_xxyz_1, g_0_yyyyzz_0_xyyz_0, g_0_yyyyzz_0_xyyz_1, g_0_yyyyzz_0_xyzz_0, g_0_yyyyzz_0_xyzz_1, g_0_yyyyzz_0_yyyy_0, g_0_yyyyzz_0_yyyy_1, g_0_yyyyzz_0_yyyz_0, g_0_yyyyzz_0_yyyz_1, g_0_yyyyzz_0_yyzz_0, g_0_yyyyzz_0_yyzz_1, g_0_yyyyzz_0_yzzz_0, g_0_yyyyzz_0_yzzz_1, g_0_yyyyzz_0_zzzz_0, g_0_yyyyzz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_xxxx_0[i] = 3.0 * g_0_xxyyzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxx_0[i] * pb_y + g_0_xxyyyzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxy_0[i] = g_0_xxyyyy_0_xxxy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxy_0[i] * pb_z + g_0_xxyyyyz_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxz_0[i] = 3.0 * g_0_xxyyzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxz_0[i] * pb_y + g_0_xxyyyzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxyy_0[i] = g_0_xxyyyy_0_xxyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxyy_0[i] * pb_z + g_0_xxyyyyz_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxyz_0[i] = g_0_yyyyzz_0_xxyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyz_0[i] * pb_x + g_0_xyyyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxzz_0[i] = 3.0 * g_0_xxyyzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxzz_0[i] * pb_y + g_0_xxyyyzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xyyy_0[i] = g_0_xxyyyy_0_xyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xyyy_0[i] * pb_z + g_0_xxyyyyz_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xyyz_0[i] = g_0_yyyyzz_0_xyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyz_0[i] * pb_x + g_0_xyyyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyzz_0[i] = g_0_yyyyzz_0_xyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyzz_0[i] * pb_x + g_0_xyyyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xzzz_0[i] = 3.0 * g_0_xxyyzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xzzz_0[i] * pb_y + g_0_xxyyyzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_yyyy_0[i] = g_0_yyyyzz_0_yyyy_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyy_0[i] * pb_x + g_0_xyyyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyz_0[i] = g_0_yyyyzz_0_yyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyz_0[i] * pb_x + g_0_xyyyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyzz_0[i] = g_0_yyyyzz_0_yyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzz_0[i] * pb_x + g_0_xyyyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yzzz_0[i] = g_0_yyyyzz_0_yzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzz_0[i] * pb_x + g_0_xyyyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_zzzz_0[i] = g_0_yyyyzz_0_zzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_zzzz_0[i] * pb_x + g_0_xyyyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 360-375 components of targeted buffer : SLSG

    auto g_0_xxyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 360);

    auto g_0_xxyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 361);

    auto g_0_xxyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 362);

    auto g_0_xxyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 363);

    auto g_0_xxyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 364);

    auto g_0_xxyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 365);

    auto g_0_xxyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 366);

    auto g_0_xxyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 367);

    auto g_0_xxyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 368);

    auto g_0_xxyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 369);

    auto g_0_xxyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 370);

    auto g_0_xxyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 371);

    auto g_0_xxyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 372);

    auto g_0_xxyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 373);

    auto g_0_xxyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 374);

    #pragma omp simd aligned(g_0_xxyyyz_0_xxxy_0, g_0_xxyyyz_0_xxxy_1, g_0_xxyyyz_0_xxyy_0, g_0_xxyyyz_0_xxyy_1, g_0_xxyyyz_0_xyyy_0, g_0_xxyyyz_0_xyyy_1, g_0_xxyyyzz_0_xxxy_0, g_0_xxyyyzz_0_xxxy_1, g_0_xxyyyzz_0_xxyy_0, g_0_xxyyyzz_0_xxyy_1, g_0_xxyyyzz_0_xyyy_0, g_0_xxyyyzz_0_xyyy_1, g_0_xxyyyzzz_0_xxxx_0, g_0_xxyyyzzz_0_xxxy_0, g_0_xxyyyzzz_0_xxxz_0, g_0_xxyyyzzz_0_xxyy_0, g_0_xxyyyzzz_0_xxyz_0, g_0_xxyyyzzz_0_xxzz_0, g_0_xxyyyzzz_0_xyyy_0, g_0_xxyyyzzz_0_xyyz_0, g_0_xxyyyzzz_0_xyzz_0, g_0_xxyyyzzz_0_xzzz_0, g_0_xxyyyzzz_0_yyyy_0, g_0_xxyyyzzz_0_yyyz_0, g_0_xxyyyzzz_0_yyzz_0, g_0_xxyyyzzz_0_yzzz_0, g_0_xxyyyzzz_0_zzzz_0, g_0_xxyyzzz_0_xxxx_0, g_0_xxyyzzz_0_xxxx_1, g_0_xxyyzzz_0_xxxz_0, g_0_xxyyzzz_0_xxxz_1, g_0_xxyyzzz_0_xxzz_0, g_0_xxyyzzz_0_xxzz_1, g_0_xxyyzzz_0_xzzz_0, g_0_xxyyzzz_0_xzzz_1, g_0_xxyzzz_0_xxxx_0, g_0_xxyzzz_0_xxxx_1, g_0_xxyzzz_0_xxxz_0, g_0_xxyzzz_0_xxxz_1, g_0_xxyzzz_0_xxzz_0, g_0_xxyzzz_0_xxzz_1, g_0_xxyzzz_0_xzzz_0, g_0_xxyzzz_0_xzzz_1, g_0_xyyyzzz_0_xxyz_0, g_0_xyyyzzz_0_xxyz_1, g_0_xyyyzzz_0_xyyz_0, g_0_xyyyzzz_0_xyyz_1, g_0_xyyyzzz_0_xyz_1, g_0_xyyyzzz_0_xyzz_0, g_0_xyyyzzz_0_xyzz_1, g_0_xyyyzzz_0_yyyy_0, g_0_xyyyzzz_0_yyyy_1, g_0_xyyyzzz_0_yyyz_0, g_0_xyyyzzz_0_yyyz_1, g_0_xyyyzzz_0_yyz_1, g_0_xyyyzzz_0_yyzz_0, g_0_xyyyzzz_0_yyzz_1, g_0_xyyyzzz_0_yzz_1, g_0_xyyyzzz_0_yzzz_0, g_0_xyyyzzz_0_yzzz_1, g_0_xyyyzzz_0_zzzz_0, g_0_xyyyzzz_0_zzzz_1, g_0_yyyzzz_0_xxyz_0, g_0_yyyzzz_0_xxyz_1, g_0_yyyzzz_0_xyyz_0, g_0_yyyzzz_0_xyyz_1, g_0_yyyzzz_0_xyzz_0, g_0_yyyzzz_0_xyzz_1, g_0_yyyzzz_0_yyyy_0, g_0_yyyzzz_0_yyyy_1, g_0_yyyzzz_0_yyyz_0, g_0_yyyzzz_0_yyyz_1, g_0_yyyzzz_0_yyzz_0, g_0_yyyzzz_0_yyzz_1, g_0_yyyzzz_0_yzzz_0, g_0_yyyzzz_0_yzzz_1, g_0_yyyzzz_0_zzzz_0, g_0_yyyzzz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_xxxx_0[i] = 2.0 * g_0_xxyzzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxx_0[i] * pb_y + g_0_xxyyzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxy_0[i] = 2.0 * g_0_xxyyyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxy_0[i] * pb_z + g_0_xxyyyzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxz_0[i] = 2.0 * g_0_xxyzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxz_0[i] * pb_y + g_0_xxyyzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxyy_0[i] = 2.0 * g_0_xxyyyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxyy_0[i] * pb_z + g_0_xxyyyzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxyz_0[i] = g_0_yyyzzz_0_xxyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyz_0[i] * pb_x + g_0_xyyyzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxzz_0[i] = 2.0 * g_0_xxyzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxzz_0[i] * pb_y + g_0_xxyyzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xyyy_0[i] = 2.0 * g_0_xxyyyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xyyy_0[i] * pb_z + g_0_xxyyyzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xyyz_0[i] = g_0_yyyzzz_0_xyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyz_0[i] * pb_x + g_0_xyyyzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyzz_0[i] = g_0_yyyzzz_0_xyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyzz_0[i] * pb_x + g_0_xyyyzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xzzz_0[i] = 2.0 * g_0_xxyzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xzzz_0[i] * pb_y + g_0_xxyyzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_yyyy_0[i] = g_0_yyyzzz_0_yyyy_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyy_0[i] * pb_x + g_0_xyyyzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyz_0[i] = g_0_yyyzzz_0_yyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyz_0[i] * pb_x + g_0_xyyyzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyzz_0[i] = g_0_yyyzzz_0_yyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzz_0[i] * pb_x + g_0_xyyyzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yzzz_0[i] = g_0_yyyzzz_0_yzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzz_0[i] * pb_x + g_0_xyyyzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_zzzz_0[i] = g_0_yyyzzz_0_zzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_zzzz_0[i] * pb_x + g_0_xyyyzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 375-390 components of targeted buffer : SLSG

    auto g_0_xxyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 375);

    auto g_0_xxyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 376);

    auto g_0_xxyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 377);

    auto g_0_xxyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 378);

    auto g_0_xxyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 379);

    auto g_0_xxyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 380);

    auto g_0_xxyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 381);

    auto g_0_xxyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 382);

    auto g_0_xxyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 383);

    auto g_0_xxyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 384);

    auto g_0_xxyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 385);

    auto g_0_xxyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 386);

    auto g_0_xxyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 387);

    auto g_0_xxyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 388);

    auto g_0_xxyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 389);

    #pragma omp simd aligned(g_0_xxyyzz_0_xxxy_0, g_0_xxyyzz_0_xxxy_1, g_0_xxyyzz_0_xxyy_0, g_0_xxyyzz_0_xxyy_1, g_0_xxyyzz_0_xyyy_0, g_0_xxyyzz_0_xyyy_1, g_0_xxyyzzz_0_xxxy_0, g_0_xxyyzzz_0_xxxy_1, g_0_xxyyzzz_0_xxyy_0, g_0_xxyyzzz_0_xxyy_1, g_0_xxyyzzz_0_xyyy_0, g_0_xxyyzzz_0_xyyy_1, g_0_xxyyzzzz_0_xxxx_0, g_0_xxyyzzzz_0_xxxy_0, g_0_xxyyzzzz_0_xxxz_0, g_0_xxyyzzzz_0_xxyy_0, g_0_xxyyzzzz_0_xxyz_0, g_0_xxyyzzzz_0_xxzz_0, g_0_xxyyzzzz_0_xyyy_0, g_0_xxyyzzzz_0_xyyz_0, g_0_xxyyzzzz_0_xyzz_0, g_0_xxyyzzzz_0_xzzz_0, g_0_xxyyzzzz_0_yyyy_0, g_0_xxyyzzzz_0_yyyz_0, g_0_xxyyzzzz_0_yyzz_0, g_0_xxyyzzzz_0_yzzz_0, g_0_xxyyzzzz_0_zzzz_0, g_0_xxyzzzz_0_xxxx_0, g_0_xxyzzzz_0_xxxx_1, g_0_xxyzzzz_0_xxxz_0, g_0_xxyzzzz_0_xxxz_1, g_0_xxyzzzz_0_xxzz_0, g_0_xxyzzzz_0_xxzz_1, g_0_xxyzzzz_0_xzzz_0, g_0_xxyzzzz_0_xzzz_1, g_0_xxzzzz_0_xxxx_0, g_0_xxzzzz_0_xxxx_1, g_0_xxzzzz_0_xxxz_0, g_0_xxzzzz_0_xxxz_1, g_0_xxzzzz_0_xxzz_0, g_0_xxzzzz_0_xxzz_1, g_0_xxzzzz_0_xzzz_0, g_0_xxzzzz_0_xzzz_1, g_0_xyyzzzz_0_xxyz_0, g_0_xyyzzzz_0_xxyz_1, g_0_xyyzzzz_0_xyyz_0, g_0_xyyzzzz_0_xyyz_1, g_0_xyyzzzz_0_xyz_1, g_0_xyyzzzz_0_xyzz_0, g_0_xyyzzzz_0_xyzz_1, g_0_xyyzzzz_0_yyyy_0, g_0_xyyzzzz_0_yyyy_1, g_0_xyyzzzz_0_yyyz_0, g_0_xyyzzzz_0_yyyz_1, g_0_xyyzzzz_0_yyz_1, g_0_xyyzzzz_0_yyzz_0, g_0_xyyzzzz_0_yyzz_1, g_0_xyyzzzz_0_yzz_1, g_0_xyyzzzz_0_yzzz_0, g_0_xyyzzzz_0_yzzz_1, g_0_xyyzzzz_0_zzzz_0, g_0_xyyzzzz_0_zzzz_1, g_0_yyzzzz_0_xxyz_0, g_0_yyzzzz_0_xxyz_1, g_0_yyzzzz_0_xyyz_0, g_0_yyzzzz_0_xyyz_1, g_0_yyzzzz_0_xyzz_0, g_0_yyzzzz_0_xyzz_1, g_0_yyzzzz_0_yyyy_0, g_0_yyzzzz_0_yyyy_1, g_0_yyzzzz_0_yyyz_0, g_0_yyzzzz_0_yyyz_1, g_0_yyzzzz_0_yyzz_0, g_0_yyzzzz_0_yyzz_1, g_0_yyzzzz_0_yzzz_0, g_0_yyzzzz_0_yzzz_1, g_0_yyzzzz_0_zzzz_0, g_0_yyzzzz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_xxxx_0[i] = g_0_xxzzzz_0_xxxx_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxx_0[i] * pb_y + g_0_xxyzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxy_0[i] = 3.0 * g_0_xxyyzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxy_0[i] * pb_z + g_0_xxyyzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxz_0[i] = g_0_xxzzzz_0_xxxz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxz_0[i] * pb_y + g_0_xxyzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxyy_0[i] = 3.0 * g_0_xxyyzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxyy_0[i] * pb_z + g_0_xxyyzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxyz_0[i] = g_0_yyzzzz_0_xxyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyz_0[i] * pb_x + g_0_xyyzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxzz_0[i] = g_0_xxzzzz_0_xxzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxzz_0[i] * pb_y + g_0_xxyzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xyyy_0[i] = 3.0 * g_0_xxyyzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xyyy_0[i] * pb_z + g_0_xxyyzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xyyz_0[i] = g_0_yyzzzz_0_xyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyz_0[i] * pb_x + g_0_xyyzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyzz_0[i] = g_0_yyzzzz_0_xyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyzz_0[i] * pb_x + g_0_xyyzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xzzz_0[i] = g_0_xxzzzz_0_xzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xzzz_0[i] * pb_y + g_0_xxyzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_yyyy_0[i] = g_0_yyzzzz_0_yyyy_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyy_0[i] * pb_x + g_0_xyyzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyz_0[i] = g_0_yyzzzz_0_yyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyz_0[i] * pb_x + g_0_xyyzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyzz_0[i] = g_0_yyzzzz_0_yyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzz_0[i] * pb_x + g_0_xyyzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yzzz_0[i] = g_0_yyzzzz_0_yzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzz_0[i] * pb_x + g_0_xyyzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_zzzz_0[i] = g_0_yyzzzz_0_zzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_zzzz_0[i] * pb_x + g_0_xyyzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 390-405 components of targeted buffer : SLSG

    auto g_0_xxyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 390);

    auto g_0_xxyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 391);

    auto g_0_xxyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 392);

    auto g_0_xxyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 393);

    auto g_0_xxyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 394);

    auto g_0_xxyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 395);

    auto g_0_xxyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 396);

    auto g_0_xxyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 397);

    auto g_0_xxyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 398);

    auto g_0_xxyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 399);

    auto g_0_xxyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 400);

    auto g_0_xxyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 401);

    auto g_0_xxyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 402);

    auto g_0_xxyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 403);

    auto g_0_xxyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 404);

    #pragma omp simd aligned(g_0_xxyzzzzz_0_xxxx_0, g_0_xxyzzzzz_0_xxxy_0, g_0_xxyzzzzz_0_xxxz_0, g_0_xxyzzzzz_0_xxyy_0, g_0_xxyzzzzz_0_xxyz_0, g_0_xxyzzzzz_0_xxzz_0, g_0_xxyzzzzz_0_xyyy_0, g_0_xxyzzzzz_0_xyyz_0, g_0_xxyzzzzz_0_xyzz_0, g_0_xxyzzzzz_0_xzzz_0, g_0_xxyzzzzz_0_yyyy_0, g_0_xxyzzzzz_0_yyyz_0, g_0_xxyzzzzz_0_yyzz_0, g_0_xxyzzzzz_0_yzzz_0, g_0_xxyzzzzz_0_zzzz_0, g_0_xxzzzzz_0_xxx_1, g_0_xxzzzzz_0_xxxx_0, g_0_xxzzzzz_0_xxxx_1, g_0_xxzzzzz_0_xxxy_0, g_0_xxzzzzz_0_xxxy_1, g_0_xxzzzzz_0_xxxz_0, g_0_xxzzzzz_0_xxxz_1, g_0_xxzzzzz_0_xxy_1, g_0_xxzzzzz_0_xxyy_0, g_0_xxzzzzz_0_xxyy_1, g_0_xxzzzzz_0_xxyz_0, g_0_xxzzzzz_0_xxyz_1, g_0_xxzzzzz_0_xxz_1, g_0_xxzzzzz_0_xxzz_0, g_0_xxzzzzz_0_xxzz_1, g_0_xxzzzzz_0_xyy_1, g_0_xxzzzzz_0_xyyy_0, g_0_xxzzzzz_0_xyyy_1, g_0_xxzzzzz_0_xyyz_0, g_0_xxzzzzz_0_xyyz_1, g_0_xxzzzzz_0_xyz_1, g_0_xxzzzzz_0_xyzz_0, g_0_xxzzzzz_0_xyzz_1, g_0_xxzzzzz_0_xzz_1, g_0_xxzzzzz_0_xzzz_0, g_0_xxzzzzz_0_xzzz_1, g_0_xxzzzzz_0_yyy_1, g_0_xxzzzzz_0_yyyy_0, g_0_xxzzzzz_0_yyyy_1, g_0_xxzzzzz_0_yyyz_0, g_0_xxzzzzz_0_yyyz_1, g_0_xxzzzzz_0_yyz_1, g_0_xxzzzzz_0_yyzz_0, g_0_xxzzzzz_0_yyzz_1, g_0_xxzzzzz_0_yzz_1, g_0_xxzzzzz_0_yzzz_0, g_0_xxzzzzz_0_yzzz_1, g_0_xxzzzzz_0_zzz_1, g_0_xxzzzzz_0_zzzz_0, g_0_xxzzzzz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_xxxx_0[i] = g_0_xxzzzzz_0_xxxx_0[i] * pb_y + g_0_xxzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxy_0[i] = g_0_xxzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxy_0[i] * pb_y + g_0_xxzzzzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxz_0[i] = g_0_xxzzzzz_0_xxxz_0[i] * pb_y + g_0_xxzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyy_0[i] = 2.0 * g_0_xxzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyy_0[i] * pb_y + g_0_xxzzzzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyz_0[i] = g_0_xxzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyz_0[i] * pb_y + g_0_xxzzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxzz_0[i] = g_0_xxzzzzz_0_xxzz_0[i] * pb_y + g_0_xxzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyy_0[i] = 3.0 * g_0_xxzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyy_0[i] * pb_y + g_0_xxzzzzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyz_0[i] = 2.0 * g_0_xxzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyz_0[i] * pb_y + g_0_xxzzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyzz_0[i] = g_0_xxzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzz_0[i] * pb_y + g_0_xxzzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xzzz_0[i] = g_0_xxzzzzz_0_xzzz_0[i] * pb_y + g_0_xxzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyy_0[i] = 4.0 * g_0_xxzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyy_0[i] * pb_y + g_0_xxzzzzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyz_0[i] = 3.0 * g_0_xxzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyz_0[i] * pb_y + g_0_xxzzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyzz_0[i] = 2.0 * g_0_xxzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyzz_0[i] * pb_y + g_0_xxzzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yzzz_0[i] = g_0_xxzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yzzz_0[i] * pb_y + g_0_xxzzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_zzzz_0[i] = g_0_xxzzzzz_0_zzzz_0[i] * pb_y + g_0_xxzzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 405-420 components of targeted buffer : SLSG

    auto g_0_xxzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 405);

    auto g_0_xxzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 406);

    auto g_0_xxzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 407);

    auto g_0_xxzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 408);

    auto g_0_xxzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 409);

    auto g_0_xxzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 410);

    auto g_0_xxzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 411);

    auto g_0_xxzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 412);

    auto g_0_xxzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 413);

    auto g_0_xxzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 414);

    auto g_0_xxzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 415);

    auto g_0_xxzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 416);

    auto g_0_xxzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 417);

    auto g_0_xxzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 418);

    auto g_0_xxzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 419);

    #pragma omp simd aligned(g_0_xxzzzz_0_xxxx_0, g_0_xxzzzz_0_xxxx_1, g_0_xxzzzz_0_xxxy_0, g_0_xxzzzz_0_xxxy_1, g_0_xxzzzz_0_xxyy_0, g_0_xxzzzz_0_xxyy_1, g_0_xxzzzz_0_xyyy_0, g_0_xxzzzz_0_xyyy_1, g_0_xxzzzzz_0_xxxx_0, g_0_xxzzzzz_0_xxxx_1, g_0_xxzzzzz_0_xxxy_0, g_0_xxzzzzz_0_xxxy_1, g_0_xxzzzzz_0_xxyy_0, g_0_xxzzzzz_0_xxyy_1, g_0_xxzzzzz_0_xyyy_0, g_0_xxzzzzz_0_xyyy_1, g_0_xxzzzzzz_0_xxxx_0, g_0_xxzzzzzz_0_xxxy_0, g_0_xxzzzzzz_0_xxxz_0, g_0_xxzzzzzz_0_xxyy_0, g_0_xxzzzzzz_0_xxyz_0, g_0_xxzzzzzz_0_xxzz_0, g_0_xxzzzzzz_0_xyyy_0, g_0_xxzzzzzz_0_xyyz_0, g_0_xxzzzzzz_0_xyzz_0, g_0_xxzzzzzz_0_xzzz_0, g_0_xxzzzzzz_0_yyyy_0, g_0_xxzzzzzz_0_yyyz_0, g_0_xxzzzzzz_0_yyzz_0, g_0_xxzzzzzz_0_yzzz_0, g_0_xxzzzzzz_0_zzzz_0, g_0_xzzzzzz_0_xxxz_0, g_0_xzzzzzz_0_xxxz_1, g_0_xzzzzzz_0_xxyz_0, g_0_xzzzzzz_0_xxyz_1, g_0_xzzzzzz_0_xxz_1, g_0_xzzzzzz_0_xxzz_0, g_0_xzzzzzz_0_xxzz_1, g_0_xzzzzzz_0_xyyz_0, g_0_xzzzzzz_0_xyyz_1, g_0_xzzzzzz_0_xyz_1, g_0_xzzzzzz_0_xyzz_0, g_0_xzzzzzz_0_xyzz_1, g_0_xzzzzzz_0_xzz_1, g_0_xzzzzzz_0_xzzz_0, g_0_xzzzzzz_0_xzzz_1, g_0_xzzzzzz_0_yyyy_0, g_0_xzzzzzz_0_yyyy_1, g_0_xzzzzzz_0_yyyz_0, g_0_xzzzzzz_0_yyyz_1, g_0_xzzzzzz_0_yyz_1, g_0_xzzzzzz_0_yyzz_0, g_0_xzzzzzz_0_yyzz_1, g_0_xzzzzzz_0_yzz_1, g_0_xzzzzzz_0_yzzz_0, g_0_xzzzzzz_0_yzzz_1, g_0_xzzzzzz_0_zzz_1, g_0_xzzzzzz_0_zzzz_0, g_0_xzzzzzz_0_zzzz_1, g_0_zzzzzz_0_xxxz_0, g_0_zzzzzz_0_xxxz_1, g_0_zzzzzz_0_xxyz_0, g_0_zzzzzz_0_xxyz_1, g_0_zzzzzz_0_xxzz_0, g_0_zzzzzz_0_xxzz_1, g_0_zzzzzz_0_xyyz_0, g_0_zzzzzz_0_xyyz_1, g_0_zzzzzz_0_xyzz_0, g_0_zzzzzz_0_xyzz_1, g_0_zzzzzz_0_xzzz_0, g_0_zzzzzz_0_xzzz_1, g_0_zzzzzz_0_yyyy_0, g_0_zzzzzz_0_yyyy_1, g_0_zzzzzz_0_yyyz_0, g_0_zzzzzz_0_yyyz_1, g_0_zzzzzz_0_yyzz_0, g_0_zzzzzz_0_yyzz_1, g_0_zzzzzz_0_yzzz_0, g_0_zzzzzz_0_yzzz_1, g_0_zzzzzz_0_zzzz_0, g_0_zzzzzz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_xxxx_0[i] = 5.0 * g_0_xxzzzz_0_xxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxx_0[i] * pb_z + g_0_xxzzzzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxy_0[i] = 5.0 * g_0_xxzzzz_0_xxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxy_0[i] * pb_z + g_0_xxzzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxz_0[i] = g_0_zzzzzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxz_0[i] * pb_x + g_0_xzzzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyy_0[i] = 5.0 * g_0_xxzzzz_0_xxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxyy_0[i] * pb_z + g_0_xxzzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxyz_0[i] = g_0_zzzzzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyz_0[i] * pb_x + g_0_xzzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxzz_0[i] = g_0_zzzzzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxzz_0[i] * pb_x + g_0_xzzzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyy_0[i] = 5.0 * g_0_xxzzzz_0_xyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xyyy_0[i] * pb_z + g_0_xxzzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xyyz_0[i] = g_0_zzzzzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyz_0[i] * pb_x + g_0_xzzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyzz_0[i] = g_0_zzzzzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyzz_0[i] * pb_x + g_0_xzzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xzzz_0[i] = g_0_zzzzzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xzzz_0[i] * pb_x + g_0_xzzzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyy_0[i] = g_0_zzzzzz_0_yyyy_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyy_0[i] * pb_x + g_0_xzzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyz_0[i] = g_0_zzzzzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyz_0[i] * pb_x + g_0_xzzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyzz_0[i] = g_0_zzzzzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzz_0[i] * pb_x + g_0_xzzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yzzz_0[i] = g_0_zzzzzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzz_0[i] * pb_x + g_0_xzzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_zzzz_0[i] = g_0_zzzzzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzz_0[i] * pb_x + g_0_xzzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 420-435 components of targeted buffer : SLSG

    auto g_0_xyyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 420);

    auto g_0_xyyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 421);

    auto g_0_xyyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 422);

    auto g_0_xyyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 423);

    auto g_0_xyyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 424);

    auto g_0_xyyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 425);

    auto g_0_xyyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 426);

    auto g_0_xyyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 427);

    auto g_0_xyyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 428);

    auto g_0_xyyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 429);

    auto g_0_xyyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 430);

    auto g_0_xyyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 431);

    auto g_0_xyyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 432);

    auto g_0_xyyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 433);

    auto g_0_xyyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 434);

    #pragma omp simd aligned(g_0_xyyyyyyy_0_xxxx_0, g_0_xyyyyyyy_0_xxxy_0, g_0_xyyyyyyy_0_xxxz_0, g_0_xyyyyyyy_0_xxyy_0, g_0_xyyyyyyy_0_xxyz_0, g_0_xyyyyyyy_0_xxzz_0, g_0_xyyyyyyy_0_xyyy_0, g_0_xyyyyyyy_0_xyyz_0, g_0_xyyyyyyy_0_xyzz_0, g_0_xyyyyyyy_0_xzzz_0, g_0_xyyyyyyy_0_yyyy_0, g_0_xyyyyyyy_0_yyyz_0, g_0_xyyyyyyy_0_yyzz_0, g_0_xyyyyyyy_0_yzzz_0, g_0_xyyyyyyy_0_zzzz_0, g_0_yyyyyyy_0_xxx_1, g_0_yyyyyyy_0_xxxx_0, g_0_yyyyyyy_0_xxxx_1, g_0_yyyyyyy_0_xxxy_0, g_0_yyyyyyy_0_xxxy_1, g_0_yyyyyyy_0_xxxz_0, g_0_yyyyyyy_0_xxxz_1, g_0_yyyyyyy_0_xxy_1, g_0_yyyyyyy_0_xxyy_0, g_0_yyyyyyy_0_xxyy_1, g_0_yyyyyyy_0_xxyz_0, g_0_yyyyyyy_0_xxyz_1, g_0_yyyyyyy_0_xxz_1, g_0_yyyyyyy_0_xxzz_0, g_0_yyyyyyy_0_xxzz_1, g_0_yyyyyyy_0_xyy_1, g_0_yyyyyyy_0_xyyy_0, g_0_yyyyyyy_0_xyyy_1, g_0_yyyyyyy_0_xyyz_0, g_0_yyyyyyy_0_xyyz_1, g_0_yyyyyyy_0_xyz_1, g_0_yyyyyyy_0_xyzz_0, g_0_yyyyyyy_0_xyzz_1, g_0_yyyyyyy_0_xzz_1, g_0_yyyyyyy_0_xzzz_0, g_0_yyyyyyy_0_xzzz_1, g_0_yyyyyyy_0_yyy_1, g_0_yyyyyyy_0_yyyy_0, g_0_yyyyyyy_0_yyyy_1, g_0_yyyyyyy_0_yyyz_0, g_0_yyyyyyy_0_yyyz_1, g_0_yyyyyyy_0_yyz_1, g_0_yyyyyyy_0_yyzz_0, g_0_yyyyyyy_0_yyzz_1, g_0_yyyyyyy_0_yzz_1, g_0_yyyyyyy_0_yzzz_0, g_0_yyyyyyy_0_yzzz_1, g_0_yyyyyyy_0_zzz_1, g_0_yyyyyyy_0_zzzz_0, g_0_yyyyyyy_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_xxxx_0[i] = 4.0 * g_0_yyyyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxx_0[i] * pb_x + g_0_yyyyyyy_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxy_0[i] = 3.0 * g_0_yyyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxy_0[i] * pb_x + g_0_yyyyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxz_0[i] = 3.0 * g_0_yyyyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxz_0[i] * pb_x + g_0_yyyyyyy_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyy_0[i] = 2.0 * g_0_yyyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyy_0[i] * pb_x + g_0_yyyyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyz_0[i] = 2.0 * g_0_yyyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyz_0[i] * pb_x + g_0_yyyyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxzz_0[i] = 2.0 * g_0_yyyyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzz_0[i] * pb_x + g_0_yyyyyyy_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyy_0[i] = g_0_yyyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyy_0[i] * pb_x + g_0_yyyyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyz_0[i] = g_0_yyyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyz_0[i] * pb_x + g_0_yyyyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyzz_0[i] = g_0_yyyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzz_0[i] * pb_x + g_0_yyyyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xzzz_0[i] = g_0_yyyyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzz_0[i] * pb_x + g_0_yyyyyyy_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyy_0[i] = g_0_yyyyyyy_0_yyyy_0[i] * pb_x + g_0_yyyyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyz_0[i] = g_0_yyyyyyy_0_yyyz_0[i] * pb_x + g_0_yyyyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyzz_0[i] = g_0_yyyyyyy_0_yyzz_0[i] * pb_x + g_0_yyyyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yzzz_0[i] = g_0_yyyyyyy_0_yzzz_0[i] * pb_x + g_0_yyyyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_zzzz_0[i] = g_0_yyyyyyy_0_zzzz_0[i] * pb_x + g_0_yyyyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 435-450 components of targeted buffer : SLSG

    auto g_0_xyyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 435);

    auto g_0_xyyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 436);

    auto g_0_xyyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 437);

    auto g_0_xyyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 438);

    auto g_0_xyyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 439);

    auto g_0_xyyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 440);

    auto g_0_xyyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 441);

    auto g_0_xyyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 442);

    auto g_0_xyyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 443);

    auto g_0_xyyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 444);

    auto g_0_xyyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 445);

    auto g_0_xyyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 446);

    auto g_0_xyyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 447);

    auto g_0_xyyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 448);

    auto g_0_xyyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 449);

    #pragma omp simd aligned(g_0_xyyyyyy_0_xxxx_0, g_0_xyyyyyy_0_xxxx_1, g_0_xyyyyyy_0_xxxy_0, g_0_xyyyyyy_0_xxxy_1, g_0_xyyyyyy_0_xxyy_0, g_0_xyyyyyy_0_xxyy_1, g_0_xyyyyyy_0_xyyy_0, g_0_xyyyyyy_0_xyyy_1, g_0_xyyyyyyz_0_xxxx_0, g_0_xyyyyyyz_0_xxxy_0, g_0_xyyyyyyz_0_xxxz_0, g_0_xyyyyyyz_0_xxyy_0, g_0_xyyyyyyz_0_xxyz_0, g_0_xyyyyyyz_0_xxzz_0, g_0_xyyyyyyz_0_xyyy_0, g_0_xyyyyyyz_0_xyyz_0, g_0_xyyyyyyz_0_xyzz_0, g_0_xyyyyyyz_0_xzzz_0, g_0_xyyyyyyz_0_yyyy_0, g_0_xyyyyyyz_0_yyyz_0, g_0_xyyyyyyz_0_yyzz_0, g_0_xyyyyyyz_0_yzzz_0, g_0_xyyyyyyz_0_zzzz_0, g_0_yyyyyyz_0_xxxz_0, g_0_yyyyyyz_0_xxxz_1, g_0_yyyyyyz_0_xxyz_0, g_0_yyyyyyz_0_xxyz_1, g_0_yyyyyyz_0_xxz_1, g_0_yyyyyyz_0_xxzz_0, g_0_yyyyyyz_0_xxzz_1, g_0_yyyyyyz_0_xyyz_0, g_0_yyyyyyz_0_xyyz_1, g_0_yyyyyyz_0_xyz_1, g_0_yyyyyyz_0_xyzz_0, g_0_yyyyyyz_0_xyzz_1, g_0_yyyyyyz_0_xzz_1, g_0_yyyyyyz_0_xzzz_0, g_0_yyyyyyz_0_xzzz_1, g_0_yyyyyyz_0_yyyy_0, g_0_yyyyyyz_0_yyyy_1, g_0_yyyyyyz_0_yyyz_0, g_0_yyyyyyz_0_yyyz_1, g_0_yyyyyyz_0_yyz_1, g_0_yyyyyyz_0_yyzz_0, g_0_yyyyyyz_0_yyzz_1, g_0_yyyyyyz_0_yzz_1, g_0_yyyyyyz_0_yzzz_0, g_0_yyyyyyz_0_yzzz_1, g_0_yyyyyyz_0_zzz_1, g_0_yyyyyyz_0_zzzz_0, g_0_yyyyyyz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyz_0_xxxx_0[i] = g_0_xyyyyyy_0_xxxx_0[i] * pb_z + g_0_xyyyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxy_0[i] = g_0_xyyyyyy_0_xxxy_0[i] * pb_z + g_0_xyyyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxz_0[i] = 3.0 * g_0_yyyyyyz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxz_0[i] * pb_x + g_0_yyyyyyz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyy_0[i] = g_0_xyyyyyy_0_xxyy_0[i] * pb_z + g_0_xyyyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxyz_0[i] = 2.0 * g_0_yyyyyyz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyz_0[i] * pb_x + g_0_yyyyyyz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyyyyz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxzz_0[i] * pb_x + g_0_yyyyyyz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyy_0[i] = g_0_xyyyyyy_0_xyyy_0[i] * pb_z + g_0_xyyyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xyyz_0[i] = g_0_yyyyyyz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyz_0[i] * pb_x + g_0_yyyyyyz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyzz_0[i] = g_0_yyyyyyz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyzz_0[i] * pb_x + g_0_yyyyyyz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xzzz_0[i] = g_0_yyyyyyz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xzzz_0[i] * pb_x + g_0_yyyyyyz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyy_0[i] = g_0_yyyyyyz_0_yyyy_0[i] * pb_x + g_0_yyyyyyz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyz_0[i] = g_0_yyyyyyz_0_yyyz_0[i] * pb_x + g_0_yyyyyyz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyzz_0[i] = g_0_yyyyyyz_0_yyzz_0[i] * pb_x + g_0_yyyyyyz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yzzz_0[i] = g_0_yyyyyyz_0_yzzz_0[i] * pb_x + g_0_yyyyyyz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_zzzz_0[i] = g_0_yyyyyyz_0_zzzz_0[i] * pb_x + g_0_yyyyyyz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 450-465 components of targeted buffer : SLSG

    auto g_0_xyyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 450);

    auto g_0_xyyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 451);

    auto g_0_xyyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 452);

    auto g_0_xyyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 453);

    auto g_0_xyyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 454);

    auto g_0_xyyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 455);

    auto g_0_xyyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 456);

    auto g_0_xyyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 457);

    auto g_0_xyyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 458);

    auto g_0_xyyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 459);

    auto g_0_xyyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 460);

    auto g_0_xyyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 461);

    auto g_0_xyyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 462);

    auto g_0_xyyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 463);

    auto g_0_xyyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 464);

    #pragma omp simd aligned(g_0_xyyyyyzz_0_xxxx_0, g_0_xyyyyyzz_0_xxxy_0, g_0_xyyyyyzz_0_xxxz_0, g_0_xyyyyyzz_0_xxyy_0, g_0_xyyyyyzz_0_xxyz_0, g_0_xyyyyyzz_0_xxzz_0, g_0_xyyyyyzz_0_xyyy_0, g_0_xyyyyyzz_0_xyyz_0, g_0_xyyyyyzz_0_xyzz_0, g_0_xyyyyyzz_0_xzzz_0, g_0_xyyyyyzz_0_yyyy_0, g_0_xyyyyyzz_0_yyyz_0, g_0_xyyyyyzz_0_yyzz_0, g_0_xyyyyyzz_0_yzzz_0, g_0_xyyyyyzz_0_zzzz_0, g_0_yyyyyzz_0_xxx_1, g_0_yyyyyzz_0_xxxx_0, g_0_yyyyyzz_0_xxxx_1, g_0_yyyyyzz_0_xxxy_0, g_0_yyyyyzz_0_xxxy_1, g_0_yyyyyzz_0_xxxz_0, g_0_yyyyyzz_0_xxxz_1, g_0_yyyyyzz_0_xxy_1, g_0_yyyyyzz_0_xxyy_0, g_0_yyyyyzz_0_xxyy_1, g_0_yyyyyzz_0_xxyz_0, g_0_yyyyyzz_0_xxyz_1, g_0_yyyyyzz_0_xxz_1, g_0_yyyyyzz_0_xxzz_0, g_0_yyyyyzz_0_xxzz_1, g_0_yyyyyzz_0_xyy_1, g_0_yyyyyzz_0_xyyy_0, g_0_yyyyyzz_0_xyyy_1, g_0_yyyyyzz_0_xyyz_0, g_0_yyyyyzz_0_xyyz_1, g_0_yyyyyzz_0_xyz_1, g_0_yyyyyzz_0_xyzz_0, g_0_yyyyyzz_0_xyzz_1, g_0_yyyyyzz_0_xzz_1, g_0_yyyyyzz_0_xzzz_0, g_0_yyyyyzz_0_xzzz_1, g_0_yyyyyzz_0_yyy_1, g_0_yyyyyzz_0_yyyy_0, g_0_yyyyyzz_0_yyyy_1, g_0_yyyyyzz_0_yyyz_0, g_0_yyyyyzz_0_yyyz_1, g_0_yyyyyzz_0_yyz_1, g_0_yyyyyzz_0_yyzz_0, g_0_yyyyyzz_0_yyzz_1, g_0_yyyyyzz_0_yzz_1, g_0_yyyyyzz_0_yzzz_0, g_0_yyyyyzz_0_yzzz_1, g_0_yyyyyzz_0_zzz_1, g_0_yyyyyzz_0_zzzz_0, g_0_yyyyyzz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_xxxx_0[i] = 4.0 * g_0_yyyyyzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxx_0[i] * pb_x + g_0_yyyyyzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxy_0[i] = 3.0 * g_0_yyyyyzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxy_0[i] * pb_x + g_0_yyyyyzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxz_0[i] = 3.0 * g_0_yyyyyzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxz_0[i] * pb_x + g_0_yyyyyzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyy_0[i] = 2.0 * g_0_yyyyyzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyy_0[i] * pb_x + g_0_yyyyyzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyz_0[i] = 2.0 * g_0_yyyyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyz_0[i] * pb_x + g_0_yyyyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxzz_0[i] = 2.0 * g_0_yyyyyzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxzz_0[i] * pb_x + g_0_yyyyyzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyy_0[i] = g_0_yyyyyzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyy_0[i] * pb_x + g_0_yyyyyzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyz_0[i] = g_0_yyyyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyz_0[i] * pb_x + g_0_yyyyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyzz_0[i] = g_0_yyyyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzz_0[i] * pb_x + g_0_yyyyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xzzz_0[i] = g_0_yyyyyzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xzzz_0[i] * pb_x + g_0_yyyyyzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyy_0[i] = g_0_yyyyyzz_0_yyyy_0[i] * pb_x + g_0_yyyyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyz_0[i] = g_0_yyyyyzz_0_yyyz_0[i] * pb_x + g_0_yyyyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyzz_0[i] = g_0_yyyyyzz_0_yyzz_0[i] * pb_x + g_0_yyyyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yzzz_0[i] = g_0_yyyyyzz_0_yzzz_0[i] * pb_x + g_0_yyyyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_zzzz_0[i] = g_0_yyyyyzz_0_zzzz_0[i] * pb_x + g_0_yyyyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 465-480 components of targeted buffer : SLSG

    auto g_0_xyyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 465);

    auto g_0_xyyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 466);

    auto g_0_xyyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 467);

    auto g_0_xyyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 468);

    auto g_0_xyyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 469);

    auto g_0_xyyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 470);

    auto g_0_xyyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 471);

    auto g_0_xyyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 472);

    auto g_0_xyyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 473);

    auto g_0_xyyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 474);

    auto g_0_xyyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 475);

    auto g_0_xyyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 476);

    auto g_0_xyyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 477);

    auto g_0_xyyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 478);

    auto g_0_xyyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 479);

    #pragma omp simd aligned(g_0_xyyyyzzz_0_xxxx_0, g_0_xyyyyzzz_0_xxxy_0, g_0_xyyyyzzz_0_xxxz_0, g_0_xyyyyzzz_0_xxyy_0, g_0_xyyyyzzz_0_xxyz_0, g_0_xyyyyzzz_0_xxzz_0, g_0_xyyyyzzz_0_xyyy_0, g_0_xyyyyzzz_0_xyyz_0, g_0_xyyyyzzz_0_xyzz_0, g_0_xyyyyzzz_0_xzzz_0, g_0_xyyyyzzz_0_yyyy_0, g_0_xyyyyzzz_0_yyyz_0, g_0_xyyyyzzz_0_yyzz_0, g_0_xyyyyzzz_0_yzzz_0, g_0_xyyyyzzz_0_zzzz_0, g_0_yyyyzzz_0_xxx_1, g_0_yyyyzzz_0_xxxx_0, g_0_yyyyzzz_0_xxxx_1, g_0_yyyyzzz_0_xxxy_0, g_0_yyyyzzz_0_xxxy_1, g_0_yyyyzzz_0_xxxz_0, g_0_yyyyzzz_0_xxxz_1, g_0_yyyyzzz_0_xxy_1, g_0_yyyyzzz_0_xxyy_0, g_0_yyyyzzz_0_xxyy_1, g_0_yyyyzzz_0_xxyz_0, g_0_yyyyzzz_0_xxyz_1, g_0_yyyyzzz_0_xxz_1, g_0_yyyyzzz_0_xxzz_0, g_0_yyyyzzz_0_xxzz_1, g_0_yyyyzzz_0_xyy_1, g_0_yyyyzzz_0_xyyy_0, g_0_yyyyzzz_0_xyyy_1, g_0_yyyyzzz_0_xyyz_0, g_0_yyyyzzz_0_xyyz_1, g_0_yyyyzzz_0_xyz_1, g_0_yyyyzzz_0_xyzz_0, g_0_yyyyzzz_0_xyzz_1, g_0_yyyyzzz_0_xzz_1, g_0_yyyyzzz_0_xzzz_0, g_0_yyyyzzz_0_xzzz_1, g_0_yyyyzzz_0_yyy_1, g_0_yyyyzzz_0_yyyy_0, g_0_yyyyzzz_0_yyyy_1, g_0_yyyyzzz_0_yyyz_0, g_0_yyyyzzz_0_yyyz_1, g_0_yyyyzzz_0_yyz_1, g_0_yyyyzzz_0_yyzz_0, g_0_yyyyzzz_0_yyzz_1, g_0_yyyyzzz_0_yzz_1, g_0_yyyyzzz_0_yzzz_0, g_0_yyyyzzz_0_yzzz_1, g_0_yyyyzzz_0_zzz_1, g_0_yyyyzzz_0_zzzz_0, g_0_yyyyzzz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_xxxx_0[i] = 4.0 * g_0_yyyyzzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxx_0[i] * pb_x + g_0_yyyyzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxy_0[i] = 3.0 * g_0_yyyyzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxy_0[i] * pb_x + g_0_yyyyzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxz_0[i] = 3.0 * g_0_yyyyzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxz_0[i] * pb_x + g_0_yyyyzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyy_0[i] = 2.0 * g_0_yyyyzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyy_0[i] * pb_x + g_0_yyyyzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyz_0[i] = 2.0 * g_0_yyyyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyz_0[i] * pb_x + g_0_yyyyzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxzz_0[i] = 2.0 * g_0_yyyyzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxzz_0[i] * pb_x + g_0_yyyyzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyy_0[i] = g_0_yyyyzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyy_0[i] * pb_x + g_0_yyyyzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyz_0[i] = g_0_yyyyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyz_0[i] * pb_x + g_0_yyyyzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyzz_0[i] = g_0_yyyyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzz_0[i] * pb_x + g_0_yyyyzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xzzz_0[i] = g_0_yyyyzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xzzz_0[i] * pb_x + g_0_yyyyzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyy_0[i] = g_0_yyyyzzz_0_yyyy_0[i] * pb_x + g_0_yyyyzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyz_0[i] = g_0_yyyyzzz_0_yyyz_0[i] * pb_x + g_0_yyyyzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyzz_0[i] = g_0_yyyyzzz_0_yyzz_0[i] * pb_x + g_0_yyyyzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yzzz_0[i] = g_0_yyyyzzz_0_yzzz_0[i] * pb_x + g_0_yyyyzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_zzzz_0[i] = g_0_yyyyzzz_0_zzzz_0[i] * pb_x + g_0_yyyyzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 480-495 components of targeted buffer : SLSG

    auto g_0_xyyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 480);

    auto g_0_xyyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 481);

    auto g_0_xyyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 482);

    auto g_0_xyyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 483);

    auto g_0_xyyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 484);

    auto g_0_xyyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 485);

    auto g_0_xyyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 486);

    auto g_0_xyyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 487);

    auto g_0_xyyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 488);

    auto g_0_xyyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 489);

    auto g_0_xyyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 490);

    auto g_0_xyyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 491);

    auto g_0_xyyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 492);

    auto g_0_xyyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 493);

    auto g_0_xyyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 494);

    #pragma omp simd aligned(g_0_xyyyzzzz_0_xxxx_0, g_0_xyyyzzzz_0_xxxy_0, g_0_xyyyzzzz_0_xxxz_0, g_0_xyyyzzzz_0_xxyy_0, g_0_xyyyzzzz_0_xxyz_0, g_0_xyyyzzzz_0_xxzz_0, g_0_xyyyzzzz_0_xyyy_0, g_0_xyyyzzzz_0_xyyz_0, g_0_xyyyzzzz_0_xyzz_0, g_0_xyyyzzzz_0_xzzz_0, g_0_xyyyzzzz_0_yyyy_0, g_0_xyyyzzzz_0_yyyz_0, g_0_xyyyzzzz_0_yyzz_0, g_0_xyyyzzzz_0_yzzz_0, g_0_xyyyzzzz_0_zzzz_0, g_0_yyyzzzz_0_xxx_1, g_0_yyyzzzz_0_xxxx_0, g_0_yyyzzzz_0_xxxx_1, g_0_yyyzzzz_0_xxxy_0, g_0_yyyzzzz_0_xxxy_1, g_0_yyyzzzz_0_xxxz_0, g_0_yyyzzzz_0_xxxz_1, g_0_yyyzzzz_0_xxy_1, g_0_yyyzzzz_0_xxyy_0, g_0_yyyzzzz_0_xxyy_1, g_0_yyyzzzz_0_xxyz_0, g_0_yyyzzzz_0_xxyz_1, g_0_yyyzzzz_0_xxz_1, g_0_yyyzzzz_0_xxzz_0, g_0_yyyzzzz_0_xxzz_1, g_0_yyyzzzz_0_xyy_1, g_0_yyyzzzz_0_xyyy_0, g_0_yyyzzzz_0_xyyy_1, g_0_yyyzzzz_0_xyyz_0, g_0_yyyzzzz_0_xyyz_1, g_0_yyyzzzz_0_xyz_1, g_0_yyyzzzz_0_xyzz_0, g_0_yyyzzzz_0_xyzz_1, g_0_yyyzzzz_0_xzz_1, g_0_yyyzzzz_0_xzzz_0, g_0_yyyzzzz_0_xzzz_1, g_0_yyyzzzz_0_yyy_1, g_0_yyyzzzz_0_yyyy_0, g_0_yyyzzzz_0_yyyy_1, g_0_yyyzzzz_0_yyyz_0, g_0_yyyzzzz_0_yyyz_1, g_0_yyyzzzz_0_yyz_1, g_0_yyyzzzz_0_yyzz_0, g_0_yyyzzzz_0_yyzz_1, g_0_yyyzzzz_0_yzz_1, g_0_yyyzzzz_0_yzzz_0, g_0_yyyzzzz_0_yzzz_1, g_0_yyyzzzz_0_zzz_1, g_0_yyyzzzz_0_zzzz_0, g_0_yyyzzzz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_xxxx_0[i] = 4.0 * g_0_yyyzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxx_0[i] * pb_x + g_0_yyyzzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxy_0[i] = 3.0 * g_0_yyyzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxy_0[i] * pb_x + g_0_yyyzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxz_0[i] = 3.0 * g_0_yyyzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxz_0[i] * pb_x + g_0_yyyzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyy_0[i] = 2.0 * g_0_yyyzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyy_0[i] * pb_x + g_0_yyyzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyz_0[i] = 2.0 * g_0_yyyzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyz_0[i] * pb_x + g_0_yyyzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxzz_0[i] = 2.0 * g_0_yyyzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxzz_0[i] * pb_x + g_0_yyyzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyy_0[i] = g_0_yyyzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyy_0[i] * pb_x + g_0_yyyzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyz_0[i] = g_0_yyyzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyz_0[i] * pb_x + g_0_yyyzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyzz_0[i] = g_0_yyyzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzz_0[i] * pb_x + g_0_yyyzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xzzz_0[i] = g_0_yyyzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xzzz_0[i] * pb_x + g_0_yyyzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyy_0[i] = g_0_yyyzzzz_0_yyyy_0[i] * pb_x + g_0_yyyzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyz_0[i] = g_0_yyyzzzz_0_yyyz_0[i] * pb_x + g_0_yyyzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyzz_0[i] = g_0_yyyzzzz_0_yyzz_0[i] * pb_x + g_0_yyyzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yzzz_0[i] = g_0_yyyzzzz_0_yzzz_0[i] * pb_x + g_0_yyyzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_zzzz_0[i] = g_0_yyyzzzz_0_zzzz_0[i] * pb_x + g_0_yyyzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 495-510 components of targeted buffer : SLSG

    auto g_0_xyyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 495);

    auto g_0_xyyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 496);

    auto g_0_xyyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 497);

    auto g_0_xyyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 498);

    auto g_0_xyyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 499);

    auto g_0_xyyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 500);

    auto g_0_xyyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 501);

    auto g_0_xyyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 502);

    auto g_0_xyyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 503);

    auto g_0_xyyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 504);

    auto g_0_xyyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 505);

    auto g_0_xyyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 506);

    auto g_0_xyyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 507);

    auto g_0_xyyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 508);

    auto g_0_xyyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 509);

    #pragma omp simd aligned(g_0_xyyzzzzz_0_xxxx_0, g_0_xyyzzzzz_0_xxxy_0, g_0_xyyzzzzz_0_xxxz_0, g_0_xyyzzzzz_0_xxyy_0, g_0_xyyzzzzz_0_xxyz_0, g_0_xyyzzzzz_0_xxzz_0, g_0_xyyzzzzz_0_xyyy_0, g_0_xyyzzzzz_0_xyyz_0, g_0_xyyzzzzz_0_xyzz_0, g_0_xyyzzzzz_0_xzzz_0, g_0_xyyzzzzz_0_yyyy_0, g_0_xyyzzzzz_0_yyyz_0, g_0_xyyzzzzz_0_yyzz_0, g_0_xyyzzzzz_0_yzzz_0, g_0_xyyzzzzz_0_zzzz_0, g_0_yyzzzzz_0_xxx_1, g_0_yyzzzzz_0_xxxx_0, g_0_yyzzzzz_0_xxxx_1, g_0_yyzzzzz_0_xxxy_0, g_0_yyzzzzz_0_xxxy_1, g_0_yyzzzzz_0_xxxz_0, g_0_yyzzzzz_0_xxxz_1, g_0_yyzzzzz_0_xxy_1, g_0_yyzzzzz_0_xxyy_0, g_0_yyzzzzz_0_xxyy_1, g_0_yyzzzzz_0_xxyz_0, g_0_yyzzzzz_0_xxyz_1, g_0_yyzzzzz_0_xxz_1, g_0_yyzzzzz_0_xxzz_0, g_0_yyzzzzz_0_xxzz_1, g_0_yyzzzzz_0_xyy_1, g_0_yyzzzzz_0_xyyy_0, g_0_yyzzzzz_0_xyyy_1, g_0_yyzzzzz_0_xyyz_0, g_0_yyzzzzz_0_xyyz_1, g_0_yyzzzzz_0_xyz_1, g_0_yyzzzzz_0_xyzz_0, g_0_yyzzzzz_0_xyzz_1, g_0_yyzzzzz_0_xzz_1, g_0_yyzzzzz_0_xzzz_0, g_0_yyzzzzz_0_xzzz_1, g_0_yyzzzzz_0_yyy_1, g_0_yyzzzzz_0_yyyy_0, g_0_yyzzzzz_0_yyyy_1, g_0_yyzzzzz_0_yyyz_0, g_0_yyzzzzz_0_yyyz_1, g_0_yyzzzzz_0_yyz_1, g_0_yyzzzzz_0_yyzz_0, g_0_yyzzzzz_0_yyzz_1, g_0_yyzzzzz_0_yzz_1, g_0_yyzzzzz_0_yzzz_0, g_0_yyzzzzz_0_yzzz_1, g_0_yyzzzzz_0_zzz_1, g_0_yyzzzzz_0_zzzz_0, g_0_yyzzzzz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_xxxx_0[i] = 4.0 * g_0_yyzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxx_0[i] * pb_x + g_0_yyzzzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxy_0[i] = 3.0 * g_0_yyzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxy_0[i] * pb_x + g_0_yyzzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxz_0[i] = 3.0 * g_0_yyzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxz_0[i] * pb_x + g_0_yyzzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyy_0[i] = 2.0 * g_0_yyzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyy_0[i] * pb_x + g_0_yyzzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyz_0[i] = 2.0 * g_0_yyzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyz_0[i] * pb_x + g_0_yyzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxzz_0[i] = 2.0 * g_0_yyzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxzz_0[i] * pb_x + g_0_yyzzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyy_0[i] = g_0_yyzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyy_0[i] * pb_x + g_0_yyzzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyz_0[i] = g_0_yyzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyz_0[i] * pb_x + g_0_yyzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyzz_0[i] = g_0_yyzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzz_0[i] * pb_x + g_0_yyzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xzzz_0[i] = g_0_yyzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xzzz_0[i] * pb_x + g_0_yyzzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyy_0[i] = g_0_yyzzzzz_0_yyyy_0[i] * pb_x + g_0_yyzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyz_0[i] = g_0_yyzzzzz_0_yyyz_0[i] * pb_x + g_0_yyzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyzz_0[i] = g_0_yyzzzzz_0_yyzz_0[i] * pb_x + g_0_yyzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yzzz_0[i] = g_0_yyzzzzz_0_yzzz_0[i] * pb_x + g_0_yyzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_zzzz_0[i] = g_0_yyzzzzz_0_zzzz_0[i] * pb_x + g_0_yyzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 510-525 components of targeted buffer : SLSG

    auto g_0_xyzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 510);

    auto g_0_xyzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 511);

    auto g_0_xyzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 512);

    auto g_0_xyzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 513);

    auto g_0_xyzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 514);

    auto g_0_xyzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 515);

    auto g_0_xyzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 516);

    auto g_0_xyzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 517);

    auto g_0_xyzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 518);

    auto g_0_xyzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 519);

    auto g_0_xyzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 520);

    auto g_0_xyzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 521);

    auto g_0_xyzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 522);

    auto g_0_xyzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 523);

    auto g_0_xyzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 524);

    #pragma omp simd aligned(g_0_xyzzzzzz_0_xxxx_0, g_0_xyzzzzzz_0_xxxy_0, g_0_xyzzzzzz_0_xxxz_0, g_0_xyzzzzzz_0_xxyy_0, g_0_xyzzzzzz_0_xxyz_0, g_0_xyzzzzzz_0_xxzz_0, g_0_xyzzzzzz_0_xyyy_0, g_0_xyzzzzzz_0_xyyz_0, g_0_xyzzzzzz_0_xyzz_0, g_0_xyzzzzzz_0_xzzz_0, g_0_xyzzzzzz_0_yyyy_0, g_0_xyzzzzzz_0_yyyz_0, g_0_xyzzzzzz_0_yyzz_0, g_0_xyzzzzzz_0_yzzz_0, g_0_xyzzzzzz_0_zzzz_0, g_0_xzzzzzz_0_xxxx_0, g_0_xzzzzzz_0_xxxx_1, g_0_xzzzzzz_0_xxxz_0, g_0_xzzzzzz_0_xxxz_1, g_0_xzzzzzz_0_xxzz_0, g_0_xzzzzzz_0_xxzz_1, g_0_xzzzzzz_0_xzzz_0, g_0_xzzzzzz_0_xzzz_1, g_0_yzzzzzz_0_xxxy_0, g_0_yzzzzzz_0_xxxy_1, g_0_yzzzzzz_0_xxy_1, g_0_yzzzzzz_0_xxyy_0, g_0_yzzzzzz_0_xxyy_1, g_0_yzzzzzz_0_xxyz_0, g_0_yzzzzzz_0_xxyz_1, g_0_yzzzzzz_0_xyy_1, g_0_yzzzzzz_0_xyyy_0, g_0_yzzzzzz_0_xyyy_1, g_0_yzzzzzz_0_xyyz_0, g_0_yzzzzzz_0_xyyz_1, g_0_yzzzzzz_0_xyz_1, g_0_yzzzzzz_0_xyzz_0, g_0_yzzzzzz_0_xyzz_1, g_0_yzzzzzz_0_yyy_1, g_0_yzzzzzz_0_yyyy_0, g_0_yzzzzzz_0_yyyy_1, g_0_yzzzzzz_0_yyyz_0, g_0_yzzzzzz_0_yyyz_1, g_0_yzzzzzz_0_yyz_1, g_0_yzzzzzz_0_yyzz_0, g_0_yzzzzzz_0_yyzz_1, g_0_yzzzzzz_0_yzz_1, g_0_yzzzzzz_0_yzzz_0, g_0_yzzzzzz_0_yzzz_1, g_0_yzzzzzz_0_zzzz_0, g_0_yzzzzzz_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzzz_0_xxxx_0[i] = g_0_xzzzzzz_0_xxxx_0[i] * pb_y + g_0_xzzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxy_0[i] = 3.0 * g_0_yzzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxy_0[i] * pb_x + g_0_yzzzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxz_0[i] = g_0_xzzzzzz_0_xxxz_0[i] * pb_y + g_0_xzzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxyy_0[i] = 2.0 * g_0_yzzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyy_0[i] * pb_x + g_0_yzzzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyz_0[i] = 2.0 * g_0_yzzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyz_0[i] * pb_x + g_0_yzzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxzz_0[i] = g_0_xzzzzzz_0_xxzz_0[i] * pb_y + g_0_xzzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xyyy_0[i] = g_0_yzzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyy_0[i] * pb_x + g_0_yzzzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyz_0[i] = g_0_yzzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyz_0[i] * pb_x + g_0_yzzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyzz_0[i] = g_0_yzzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzz_0[i] * pb_x + g_0_yzzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xzzz_0[i] = g_0_xzzzzzz_0_xzzz_0[i] * pb_y + g_0_xzzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_yyyy_0[i] = g_0_yzzzzzz_0_yyyy_0[i] * pb_x + g_0_yzzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyz_0[i] = g_0_yzzzzzz_0_yyyz_0[i] * pb_x + g_0_yzzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyzz_0[i] = g_0_yzzzzzz_0_yyzz_0[i] * pb_x + g_0_yzzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yzzz_0[i] = g_0_yzzzzzz_0_yzzz_0[i] * pb_x + g_0_yzzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_zzzz_0[i] = g_0_yzzzzzz_0_zzzz_0[i] * pb_x + g_0_yzzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 525-540 components of targeted buffer : SLSG

    auto g_0_xzzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 525);

    auto g_0_xzzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 526);

    auto g_0_xzzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 527);

    auto g_0_xzzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 528);

    auto g_0_xzzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 529);

    auto g_0_xzzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 530);

    auto g_0_xzzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 531);

    auto g_0_xzzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 532);

    auto g_0_xzzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 533);

    auto g_0_xzzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 534);

    auto g_0_xzzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 535);

    auto g_0_xzzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 536);

    auto g_0_xzzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 537);

    auto g_0_xzzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 538);

    auto g_0_xzzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 539);

    #pragma omp simd aligned(g_0_xzzzzzzz_0_xxxx_0, g_0_xzzzzzzz_0_xxxy_0, g_0_xzzzzzzz_0_xxxz_0, g_0_xzzzzzzz_0_xxyy_0, g_0_xzzzzzzz_0_xxyz_0, g_0_xzzzzzzz_0_xxzz_0, g_0_xzzzzzzz_0_xyyy_0, g_0_xzzzzzzz_0_xyyz_0, g_0_xzzzzzzz_0_xyzz_0, g_0_xzzzzzzz_0_xzzz_0, g_0_xzzzzzzz_0_yyyy_0, g_0_xzzzzzzz_0_yyyz_0, g_0_xzzzzzzz_0_yyzz_0, g_0_xzzzzzzz_0_yzzz_0, g_0_xzzzzzzz_0_zzzz_0, g_0_zzzzzzz_0_xxx_1, g_0_zzzzzzz_0_xxxx_0, g_0_zzzzzzz_0_xxxx_1, g_0_zzzzzzz_0_xxxy_0, g_0_zzzzzzz_0_xxxy_1, g_0_zzzzzzz_0_xxxz_0, g_0_zzzzzzz_0_xxxz_1, g_0_zzzzzzz_0_xxy_1, g_0_zzzzzzz_0_xxyy_0, g_0_zzzzzzz_0_xxyy_1, g_0_zzzzzzz_0_xxyz_0, g_0_zzzzzzz_0_xxyz_1, g_0_zzzzzzz_0_xxz_1, g_0_zzzzzzz_0_xxzz_0, g_0_zzzzzzz_0_xxzz_1, g_0_zzzzzzz_0_xyy_1, g_0_zzzzzzz_0_xyyy_0, g_0_zzzzzzz_0_xyyy_1, g_0_zzzzzzz_0_xyyz_0, g_0_zzzzzzz_0_xyyz_1, g_0_zzzzzzz_0_xyz_1, g_0_zzzzzzz_0_xyzz_0, g_0_zzzzzzz_0_xyzz_1, g_0_zzzzzzz_0_xzz_1, g_0_zzzzzzz_0_xzzz_0, g_0_zzzzzzz_0_xzzz_1, g_0_zzzzzzz_0_yyy_1, g_0_zzzzzzz_0_yyyy_0, g_0_zzzzzzz_0_yyyy_1, g_0_zzzzzzz_0_yyyz_0, g_0_zzzzzzz_0_yyyz_1, g_0_zzzzzzz_0_yyz_1, g_0_zzzzzzz_0_yyzz_0, g_0_zzzzzzz_0_yyzz_1, g_0_zzzzzzz_0_yzz_1, g_0_zzzzzzz_0_yzzz_0, g_0_zzzzzzz_0_yzzz_1, g_0_zzzzzzz_0_zzz_1, g_0_zzzzzzz_0_zzzz_0, g_0_zzzzzzz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_xxxx_0[i] = 4.0 * g_0_zzzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxx_0[i] * pb_x + g_0_zzzzzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxy_0[i] = 3.0 * g_0_zzzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxy_0[i] * pb_x + g_0_zzzzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxz_0[i] = 3.0 * g_0_zzzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxz_0[i] * pb_x + g_0_zzzzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyy_0[i] * pb_x + g_0_zzzzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyz_0[i] = 2.0 * g_0_zzzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyz_0[i] * pb_x + g_0_zzzzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxzz_0[i] = 2.0 * g_0_zzzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzz_0[i] * pb_x + g_0_zzzzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyy_0[i] = g_0_zzzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyy_0[i] * pb_x + g_0_zzzzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyz_0[i] = g_0_zzzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyz_0[i] * pb_x + g_0_zzzzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyzz_0[i] = g_0_zzzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzz_0[i] * pb_x + g_0_zzzzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xzzz_0[i] = g_0_zzzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzz_0[i] * pb_x + g_0_zzzzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyy_0[i] = g_0_zzzzzzz_0_yyyy_0[i] * pb_x + g_0_zzzzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyz_0[i] = g_0_zzzzzzz_0_yyyz_0[i] * pb_x + g_0_zzzzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyzz_0[i] = g_0_zzzzzzz_0_yyzz_0[i] * pb_x + g_0_zzzzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yzzz_0[i] = g_0_zzzzzzz_0_yzzz_0[i] * pb_x + g_0_zzzzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_zzzz_0[i] = g_0_zzzzzzz_0_zzzz_0[i] * pb_x + g_0_zzzzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 540-555 components of targeted buffer : SLSG

    auto g_0_yyyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 540);

    auto g_0_yyyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 541);

    auto g_0_yyyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 542);

    auto g_0_yyyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 543);

    auto g_0_yyyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 544);

    auto g_0_yyyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 545);

    auto g_0_yyyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 546);

    auto g_0_yyyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 547);

    auto g_0_yyyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 548);

    auto g_0_yyyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 549);

    auto g_0_yyyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 550);

    auto g_0_yyyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 551);

    auto g_0_yyyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 552);

    auto g_0_yyyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 553);

    auto g_0_yyyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 554);

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxx_0, g_0_yyyyyy_0_xxxx_1, g_0_yyyyyy_0_xxxy_0, g_0_yyyyyy_0_xxxy_1, g_0_yyyyyy_0_xxxz_0, g_0_yyyyyy_0_xxxz_1, g_0_yyyyyy_0_xxyy_0, g_0_yyyyyy_0_xxyy_1, g_0_yyyyyy_0_xxyz_0, g_0_yyyyyy_0_xxyz_1, g_0_yyyyyy_0_xxzz_0, g_0_yyyyyy_0_xxzz_1, g_0_yyyyyy_0_xyyy_0, g_0_yyyyyy_0_xyyy_1, g_0_yyyyyy_0_xyyz_0, g_0_yyyyyy_0_xyyz_1, g_0_yyyyyy_0_xyzz_0, g_0_yyyyyy_0_xyzz_1, g_0_yyyyyy_0_xzzz_0, g_0_yyyyyy_0_xzzz_1, g_0_yyyyyy_0_yyyy_0, g_0_yyyyyy_0_yyyy_1, g_0_yyyyyy_0_yyyz_0, g_0_yyyyyy_0_yyyz_1, g_0_yyyyyy_0_yyzz_0, g_0_yyyyyy_0_yyzz_1, g_0_yyyyyy_0_yzzz_0, g_0_yyyyyy_0_yzzz_1, g_0_yyyyyy_0_zzzz_0, g_0_yyyyyy_0_zzzz_1, g_0_yyyyyyy_0_xxx_1, g_0_yyyyyyy_0_xxxx_0, g_0_yyyyyyy_0_xxxx_1, g_0_yyyyyyy_0_xxxy_0, g_0_yyyyyyy_0_xxxy_1, g_0_yyyyyyy_0_xxxz_0, g_0_yyyyyyy_0_xxxz_1, g_0_yyyyyyy_0_xxy_1, g_0_yyyyyyy_0_xxyy_0, g_0_yyyyyyy_0_xxyy_1, g_0_yyyyyyy_0_xxyz_0, g_0_yyyyyyy_0_xxyz_1, g_0_yyyyyyy_0_xxz_1, g_0_yyyyyyy_0_xxzz_0, g_0_yyyyyyy_0_xxzz_1, g_0_yyyyyyy_0_xyy_1, g_0_yyyyyyy_0_xyyy_0, g_0_yyyyyyy_0_xyyy_1, g_0_yyyyyyy_0_xyyz_0, g_0_yyyyyyy_0_xyyz_1, g_0_yyyyyyy_0_xyz_1, g_0_yyyyyyy_0_xyzz_0, g_0_yyyyyyy_0_xyzz_1, g_0_yyyyyyy_0_xzz_1, g_0_yyyyyyy_0_xzzz_0, g_0_yyyyyyy_0_xzzz_1, g_0_yyyyyyy_0_yyy_1, g_0_yyyyyyy_0_yyyy_0, g_0_yyyyyyy_0_yyyy_1, g_0_yyyyyyy_0_yyyz_0, g_0_yyyyyyy_0_yyyz_1, g_0_yyyyyyy_0_yyz_1, g_0_yyyyyyy_0_yyzz_0, g_0_yyyyyyy_0_yyzz_1, g_0_yyyyyyy_0_yzz_1, g_0_yyyyyyy_0_yzzz_0, g_0_yyyyyyy_0_yzzz_1, g_0_yyyyyyy_0_zzz_1, g_0_yyyyyyy_0_zzzz_0, g_0_yyyyyyy_0_zzzz_1, g_0_yyyyyyyy_0_xxxx_0, g_0_yyyyyyyy_0_xxxy_0, g_0_yyyyyyyy_0_xxxz_0, g_0_yyyyyyyy_0_xxyy_0, g_0_yyyyyyyy_0_xxyz_0, g_0_yyyyyyyy_0_xxzz_0, g_0_yyyyyyyy_0_xyyy_0, g_0_yyyyyyyy_0_xyyz_0, g_0_yyyyyyyy_0_xyzz_0, g_0_yyyyyyyy_0_xzzz_0, g_0_yyyyyyyy_0_yyyy_0, g_0_yyyyyyyy_0_yyyz_0, g_0_yyyyyyyy_0_yyzz_0, g_0_yyyyyyyy_0_yzzz_0, g_0_yyyyyyyy_0_zzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_xxxx_0[i] = 7.0 * g_0_yyyyyy_0_xxxx_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxx_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxx_0[i] * pb_y + g_0_yyyyyyy_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxy_0[i] = 7.0 * g_0_yyyyyy_0_xxxy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxy_0[i] * pb_y + g_0_yyyyyyy_0_xxxy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxz_0[i] = 7.0 * g_0_yyyyyy_0_xxxz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxz_0[i] * pb_y + g_0_yyyyyyy_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyy_0[i] = 7.0 * g_0_yyyyyy_0_xxyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyy_0[i] * pb_y + g_0_yyyyyyy_0_xxyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyz_0[i] = 7.0 * g_0_yyyyyy_0_xxyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyz_0[i] * pb_y + g_0_yyyyyyy_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxzz_0[i] = 7.0 * g_0_yyyyyy_0_xxzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxzz_0[i] * pb_y + g_0_yyyyyyy_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyy_0[i] = 7.0 * g_0_yyyyyy_0_xyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyy_0[i] * pb_y + g_0_yyyyyyy_0_xyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyz_0[i] = 7.0 * g_0_yyyyyy_0_xyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyz_0[i] * pb_y + g_0_yyyyyyy_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyzz_0[i] = 7.0 * g_0_yyyyyy_0_xyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzz_0[i] * pb_y + g_0_yyyyyyy_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xzzz_0[i] = 7.0 * g_0_yyyyyy_0_xzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzzz_0[i] * pb_y + g_0_yyyyyyy_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyy_0[i] = 7.0 * g_0_yyyyyy_0_yyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyy_0[i] * pb_y + g_0_yyyyyyy_0_yyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyz_0[i] = 7.0 * g_0_yyyyyy_0_yyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyz_0[i] * pb_y + g_0_yyyyyyy_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyzz_0[i] = 7.0 * g_0_yyyyyy_0_yyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzz_0[i] * pb_y + g_0_yyyyyyy_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yzzz_0[i] = 7.0 * g_0_yyyyyy_0_yzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzz_0[i] * pb_y + g_0_yyyyyyy_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_zzzz_0[i] = 7.0 * g_0_yyyyyy_0_zzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_zzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzzz_0[i] * pb_y + g_0_yyyyyyy_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 555-570 components of targeted buffer : SLSG

    auto g_0_yyyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 555);

    auto g_0_yyyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 556);

    auto g_0_yyyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 557);

    auto g_0_yyyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 558);

    auto g_0_yyyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 559);

    auto g_0_yyyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 560);

    auto g_0_yyyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 561);

    auto g_0_yyyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 562);

    auto g_0_yyyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 563);

    auto g_0_yyyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 564);

    auto g_0_yyyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 565);

    auto g_0_yyyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 566);

    auto g_0_yyyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 567);

    auto g_0_yyyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 568);

    auto g_0_yyyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 569);

    #pragma omp simd aligned(g_0_yyyyyyy_0_xxx_1, g_0_yyyyyyy_0_xxxx_0, g_0_yyyyyyy_0_xxxx_1, g_0_yyyyyyy_0_xxxy_0, g_0_yyyyyyy_0_xxxy_1, g_0_yyyyyyy_0_xxxz_0, g_0_yyyyyyy_0_xxxz_1, g_0_yyyyyyy_0_xxy_1, g_0_yyyyyyy_0_xxyy_0, g_0_yyyyyyy_0_xxyy_1, g_0_yyyyyyy_0_xxyz_0, g_0_yyyyyyy_0_xxyz_1, g_0_yyyyyyy_0_xxz_1, g_0_yyyyyyy_0_xxzz_0, g_0_yyyyyyy_0_xxzz_1, g_0_yyyyyyy_0_xyy_1, g_0_yyyyyyy_0_xyyy_0, g_0_yyyyyyy_0_xyyy_1, g_0_yyyyyyy_0_xyyz_0, g_0_yyyyyyy_0_xyyz_1, g_0_yyyyyyy_0_xyz_1, g_0_yyyyyyy_0_xyzz_0, g_0_yyyyyyy_0_xyzz_1, g_0_yyyyyyy_0_xzz_1, g_0_yyyyyyy_0_xzzz_0, g_0_yyyyyyy_0_xzzz_1, g_0_yyyyyyy_0_yyy_1, g_0_yyyyyyy_0_yyyy_0, g_0_yyyyyyy_0_yyyy_1, g_0_yyyyyyy_0_yyyz_0, g_0_yyyyyyy_0_yyyz_1, g_0_yyyyyyy_0_yyz_1, g_0_yyyyyyy_0_yyzz_0, g_0_yyyyyyy_0_yyzz_1, g_0_yyyyyyy_0_yzz_1, g_0_yyyyyyy_0_yzzz_0, g_0_yyyyyyy_0_yzzz_1, g_0_yyyyyyy_0_zzz_1, g_0_yyyyyyy_0_zzzz_0, g_0_yyyyyyy_0_zzzz_1, g_0_yyyyyyyz_0_xxxx_0, g_0_yyyyyyyz_0_xxxy_0, g_0_yyyyyyyz_0_xxxz_0, g_0_yyyyyyyz_0_xxyy_0, g_0_yyyyyyyz_0_xxyz_0, g_0_yyyyyyyz_0_xxzz_0, g_0_yyyyyyyz_0_xyyy_0, g_0_yyyyyyyz_0_xyyz_0, g_0_yyyyyyyz_0_xyzz_0, g_0_yyyyyyyz_0_xzzz_0, g_0_yyyyyyyz_0_yyyy_0, g_0_yyyyyyyz_0_yyyz_0, g_0_yyyyyyyz_0_yyzz_0, g_0_yyyyyyyz_0_yzzz_0, g_0_yyyyyyyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_xxxx_0[i] = g_0_yyyyyyy_0_xxxx_0[i] * pb_z + g_0_yyyyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxy_0[i] = g_0_yyyyyyy_0_xxxy_0[i] * pb_z + g_0_yyyyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxz_0[i] = g_0_yyyyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxz_0[i] * pb_z + g_0_yyyyyyy_0_xxxz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyy_0[i] = g_0_yyyyyyy_0_xxyy_0[i] * pb_z + g_0_yyyyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyz_0[i] = g_0_yyyyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyz_0[i] * pb_z + g_0_yyyyyyy_0_xxyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzz_0[i] * pb_z + g_0_yyyyyyy_0_xxzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyy_0[i] = g_0_yyyyyyy_0_xyyy_0[i] * pb_z + g_0_yyyyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyz_0[i] = g_0_yyyyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyz_0[i] * pb_z + g_0_yyyyyyy_0_xyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzz_0[i] * pb_z + g_0_yyyyyyy_0_xyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzz_0[i] * pb_z + g_0_yyyyyyy_0_xzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyy_0[i] = g_0_yyyyyyy_0_yyyy_0[i] * pb_z + g_0_yyyyyyy_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyz_0[i] = g_0_yyyyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyz_0[i] * pb_z + g_0_yyyyyyy_0_yyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyzz_0[i] = 2.0 * g_0_yyyyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzz_0[i] * pb_z + g_0_yyyyyyy_0_yyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yzzz_0[i] = 3.0 * g_0_yyyyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzz_0[i] * pb_z + g_0_yyyyyyy_0_yzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_zzzz_0[i] = 4.0 * g_0_yyyyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_zzzz_0[i] * pb_z + g_0_yyyyyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 570-585 components of targeted buffer : SLSG

    auto g_0_yyyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 570);

    auto g_0_yyyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 571);

    auto g_0_yyyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 572);

    auto g_0_yyyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 573);

    auto g_0_yyyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 574);

    auto g_0_yyyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 575);

    auto g_0_yyyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 576);

    auto g_0_yyyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 577);

    auto g_0_yyyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 578);

    auto g_0_yyyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 579);

    auto g_0_yyyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 580);

    auto g_0_yyyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 581);

    auto g_0_yyyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 582);

    auto g_0_yyyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 583);

    auto g_0_yyyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 584);

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxy_0, g_0_yyyyyy_0_xxxy_1, g_0_yyyyyy_0_xxyy_0, g_0_yyyyyy_0_xxyy_1, g_0_yyyyyy_0_xyyy_0, g_0_yyyyyy_0_xyyy_1, g_0_yyyyyy_0_yyyy_0, g_0_yyyyyy_0_yyyy_1, g_0_yyyyyyz_0_xxxy_0, g_0_yyyyyyz_0_xxxy_1, g_0_yyyyyyz_0_xxyy_0, g_0_yyyyyyz_0_xxyy_1, g_0_yyyyyyz_0_xyyy_0, g_0_yyyyyyz_0_xyyy_1, g_0_yyyyyyz_0_yyyy_0, g_0_yyyyyyz_0_yyyy_1, g_0_yyyyyyzz_0_xxxx_0, g_0_yyyyyyzz_0_xxxy_0, g_0_yyyyyyzz_0_xxxz_0, g_0_yyyyyyzz_0_xxyy_0, g_0_yyyyyyzz_0_xxyz_0, g_0_yyyyyyzz_0_xxzz_0, g_0_yyyyyyzz_0_xyyy_0, g_0_yyyyyyzz_0_xyyz_0, g_0_yyyyyyzz_0_xyzz_0, g_0_yyyyyyzz_0_xzzz_0, g_0_yyyyyyzz_0_yyyy_0, g_0_yyyyyyzz_0_yyyz_0, g_0_yyyyyyzz_0_yyzz_0, g_0_yyyyyyzz_0_yzzz_0, g_0_yyyyyyzz_0_zzzz_0, g_0_yyyyyzz_0_xxxx_0, g_0_yyyyyzz_0_xxxx_1, g_0_yyyyyzz_0_xxxz_0, g_0_yyyyyzz_0_xxxz_1, g_0_yyyyyzz_0_xxyz_0, g_0_yyyyyzz_0_xxyz_1, g_0_yyyyyzz_0_xxz_1, g_0_yyyyyzz_0_xxzz_0, g_0_yyyyyzz_0_xxzz_1, g_0_yyyyyzz_0_xyyz_0, g_0_yyyyyzz_0_xyyz_1, g_0_yyyyyzz_0_xyz_1, g_0_yyyyyzz_0_xyzz_0, g_0_yyyyyzz_0_xyzz_1, g_0_yyyyyzz_0_xzz_1, g_0_yyyyyzz_0_xzzz_0, g_0_yyyyyzz_0_xzzz_1, g_0_yyyyyzz_0_yyyz_0, g_0_yyyyyzz_0_yyyz_1, g_0_yyyyyzz_0_yyz_1, g_0_yyyyyzz_0_yyzz_0, g_0_yyyyyzz_0_yyzz_1, g_0_yyyyyzz_0_yzz_1, g_0_yyyyyzz_0_yzzz_0, g_0_yyyyyzz_0_yzzz_1, g_0_yyyyyzz_0_zzz_1, g_0_yyyyyzz_0_zzzz_0, g_0_yyyyyzz_0_zzzz_1, g_0_yyyyzz_0_xxxx_0, g_0_yyyyzz_0_xxxx_1, g_0_yyyyzz_0_xxxz_0, g_0_yyyyzz_0_xxxz_1, g_0_yyyyzz_0_xxyz_0, g_0_yyyyzz_0_xxyz_1, g_0_yyyyzz_0_xxzz_0, g_0_yyyyzz_0_xxzz_1, g_0_yyyyzz_0_xyyz_0, g_0_yyyyzz_0_xyyz_1, g_0_yyyyzz_0_xyzz_0, g_0_yyyyzz_0_xyzz_1, g_0_yyyyzz_0_xzzz_0, g_0_yyyyzz_0_xzzz_1, g_0_yyyyzz_0_yyyz_0, g_0_yyyyzz_0_yyyz_1, g_0_yyyyzz_0_yyzz_0, g_0_yyyyzz_0_yyzz_1, g_0_yyyyzz_0_yzzz_0, g_0_yyyyzz_0_yzzz_1, g_0_yyyyzz_0_zzzz_0, g_0_yyyyzz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_xxxx_0[i] = 5.0 * g_0_yyyyzz_0_xxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxx_0[i] * pb_y + g_0_yyyyyzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxy_0[i] = g_0_yyyyyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxy_0[i] * pb_z + g_0_yyyyyyz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxz_0[i] = 5.0 * g_0_yyyyzz_0_xxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxz_0[i] * pb_y + g_0_yyyyyzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyy_0[i] = g_0_yyyyyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxyy_0[i] * pb_z + g_0_yyyyyyz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxyz_0[i] = 5.0 * g_0_yyyyzz_0_xxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyz_0[i] * pb_y + g_0_yyyyyzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxzz_0[i] = 5.0 * g_0_yyyyzz_0_xxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxzz_0[i] * pb_y + g_0_yyyyyzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyy_0[i] = g_0_yyyyyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xyyy_0[i] * pb_z + g_0_yyyyyyz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xyyz_0[i] = 5.0 * g_0_yyyyzz_0_xyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyz_0[i] * pb_y + g_0_yyyyyzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyzz_0[i] = 5.0 * g_0_yyyyzz_0_xyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzz_0[i] * pb_y + g_0_yyyyyzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xzzz_0[i] = 5.0 * g_0_yyyyzz_0_xzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzzz_0[i] * pb_y + g_0_yyyyyzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyy_0[i] = g_0_yyyyyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_yyyy_0[i] * pb_z + g_0_yyyyyyz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_yyyz_0[i] = 5.0 * g_0_yyyyzz_0_yyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyz_0[i] * pb_y + g_0_yyyyyzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyzz_0[i] = 5.0 * g_0_yyyyzz_0_yyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyzz_0[i] * pb_y + g_0_yyyyyzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yzzz_0[i] = 5.0 * g_0_yyyyzz_0_yzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yzzz_0[i] * pb_y + g_0_yyyyyzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_zzzz_0[i] = 5.0 * g_0_yyyyzz_0_zzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzzz_0[i] * pb_y + g_0_yyyyyzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 585-600 components of targeted buffer : SLSG

    auto g_0_yyyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 585);

    auto g_0_yyyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 586);

    auto g_0_yyyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 587);

    auto g_0_yyyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 588);

    auto g_0_yyyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 589);

    auto g_0_yyyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 590);

    auto g_0_yyyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 591);

    auto g_0_yyyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 592);

    auto g_0_yyyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 593);

    auto g_0_yyyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 594);

    auto g_0_yyyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 595);

    auto g_0_yyyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 596);

    auto g_0_yyyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 597);

    auto g_0_yyyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 598);

    auto g_0_yyyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 599);

    #pragma omp simd aligned(g_0_yyyyyz_0_xxxy_0, g_0_yyyyyz_0_xxxy_1, g_0_yyyyyz_0_xxyy_0, g_0_yyyyyz_0_xxyy_1, g_0_yyyyyz_0_xyyy_0, g_0_yyyyyz_0_xyyy_1, g_0_yyyyyz_0_yyyy_0, g_0_yyyyyz_0_yyyy_1, g_0_yyyyyzz_0_xxxy_0, g_0_yyyyyzz_0_xxxy_1, g_0_yyyyyzz_0_xxyy_0, g_0_yyyyyzz_0_xxyy_1, g_0_yyyyyzz_0_xyyy_0, g_0_yyyyyzz_0_xyyy_1, g_0_yyyyyzz_0_yyyy_0, g_0_yyyyyzz_0_yyyy_1, g_0_yyyyyzzz_0_xxxx_0, g_0_yyyyyzzz_0_xxxy_0, g_0_yyyyyzzz_0_xxxz_0, g_0_yyyyyzzz_0_xxyy_0, g_0_yyyyyzzz_0_xxyz_0, g_0_yyyyyzzz_0_xxzz_0, g_0_yyyyyzzz_0_xyyy_0, g_0_yyyyyzzz_0_xyyz_0, g_0_yyyyyzzz_0_xyzz_0, g_0_yyyyyzzz_0_xzzz_0, g_0_yyyyyzzz_0_yyyy_0, g_0_yyyyyzzz_0_yyyz_0, g_0_yyyyyzzz_0_yyzz_0, g_0_yyyyyzzz_0_yzzz_0, g_0_yyyyyzzz_0_zzzz_0, g_0_yyyyzzz_0_xxxx_0, g_0_yyyyzzz_0_xxxx_1, g_0_yyyyzzz_0_xxxz_0, g_0_yyyyzzz_0_xxxz_1, g_0_yyyyzzz_0_xxyz_0, g_0_yyyyzzz_0_xxyz_1, g_0_yyyyzzz_0_xxz_1, g_0_yyyyzzz_0_xxzz_0, g_0_yyyyzzz_0_xxzz_1, g_0_yyyyzzz_0_xyyz_0, g_0_yyyyzzz_0_xyyz_1, g_0_yyyyzzz_0_xyz_1, g_0_yyyyzzz_0_xyzz_0, g_0_yyyyzzz_0_xyzz_1, g_0_yyyyzzz_0_xzz_1, g_0_yyyyzzz_0_xzzz_0, g_0_yyyyzzz_0_xzzz_1, g_0_yyyyzzz_0_yyyz_0, g_0_yyyyzzz_0_yyyz_1, g_0_yyyyzzz_0_yyz_1, g_0_yyyyzzz_0_yyzz_0, g_0_yyyyzzz_0_yyzz_1, g_0_yyyyzzz_0_yzz_1, g_0_yyyyzzz_0_yzzz_0, g_0_yyyyzzz_0_yzzz_1, g_0_yyyyzzz_0_zzz_1, g_0_yyyyzzz_0_zzzz_0, g_0_yyyyzzz_0_zzzz_1, g_0_yyyzzz_0_xxxx_0, g_0_yyyzzz_0_xxxx_1, g_0_yyyzzz_0_xxxz_0, g_0_yyyzzz_0_xxxz_1, g_0_yyyzzz_0_xxyz_0, g_0_yyyzzz_0_xxyz_1, g_0_yyyzzz_0_xxzz_0, g_0_yyyzzz_0_xxzz_1, g_0_yyyzzz_0_xyyz_0, g_0_yyyzzz_0_xyyz_1, g_0_yyyzzz_0_xyzz_0, g_0_yyyzzz_0_xyzz_1, g_0_yyyzzz_0_xzzz_0, g_0_yyyzzz_0_xzzz_1, g_0_yyyzzz_0_yyyz_0, g_0_yyyzzz_0_yyyz_1, g_0_yyyzzz_0_yyzz_0, g_0_yyyzzz_0_yyzz_1, g_0_yyyzzz_0_yzzz_0, g_0_yyyzzz_0_yzzz_1, g_0_yyyzzz_0_zzzz_0, g_0_yyyzzz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_xxxx_0[i] = 4.0 * g_0_yyyzzz_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxx_0[i] * pb_y + g_0_yyyyzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxy_0[i] = 2.0 * g_0_yyyyyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxy_0[i] * pb_z + g_0_yyyyyzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxz_0[i] = 4.0 * g_0_yyyzzz_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxz_0[i] * pb_y + g_0_yyyyzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyy_0[i] = 2.0 * g_0_yyyyyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxyy_0[i] * pb_z + g_0_yyyyyzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxyz_0[i] = 4.0 * g_0_yyyzzz_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyz_0[i] * pb_y + g_0_yyyyzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxzz_0[i] = 4.0 * g_0_yyyzzz_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxzz_0[i] * pb_y + g_0_yyyyzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyy_0[i] = 2.0 * g_0_yyyyyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xyyy_0[i] * pb_z + g_0_yyyyyzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xyyz_0[i] = 4.0 * g_0_yyyzzz_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyz_0[i] * pb_y + g_0_yyyyzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyzz_0[i] = 4.0 * g_0_yyyzzz_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzz_0[i] * pb_y + g_0_yyyyzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xzzz_0[i] = 4.0 * g_0_yyyzzz_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzzz_0[i] * pb_y + g_0_yyyyzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyy_0[i] = 2.0 * g_0_yyyyyz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_yyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_yyyy_0[i] * pb_z + g_0_yyyyyzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_yyyz_0[i] = 4.0 * g_0_yyyzzz_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyz_0[i] * pb_y + g_0_yyyyzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyzz_0[i] = 4.0 * g_0_yyyzzz_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyzz_0[i] * pb_y + g_0_yyyyzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yzzz_0[i] = 4.0 * g_0_yyyzzz_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yzzz_0[i] * pb_y + g_0_yyyyzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_zzzz_0[i] = 4.0 * g_0_yyyzzz_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzzz_0[i] * pb_y + g_0_yyyyzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 600-615 components of targeted buffer : SLSG

    auto g_0_yyyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 600);

    auto g_0_yyyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 601);

    auto g_0_yyyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 602);

    auto g_0_yyyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 603);

    auto g_0_yyyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 604);

    auto g_0_yyyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 605);

    auto g_0_yyyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 606);

    auto g_0_yyyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 607);

    auto g_0_yyyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 608);

    auto g_0_yyyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 609);

    auto g_0_yyyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 610);

    auto g_0_yyyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 611);

    auto g_0_yyyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 612);

    auto g_0_yyyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 613);

    auto g_0_yyyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 614);

    #pragma omp simd aligned(g_0_yyyyzz_0_xxxy_0, g_0_yyyyzz_0_xxxy_1, g_0_yyyyzz_0_xxyy_0, g_0_yyyyzz_0_xxyy_1, g_0_yyyyzz_0_xyyy_0, g_0_yyyyzz_0_xyyy_1, g_0_yyyyzz_0_yyyy_0, g_0_yyyyzz_0_yyyy_1, g_0_yyyyzzz_0_xxxy_0, g_0_yyyyzzz_0_xxxy_1, g_0_yyyyzzz_0_xxyy_0, g_0_yyyyzzz_0_xxyy_1, g_0_yyyyzzz_0_xyyy_0, g_0_yyyyzzz_0_xyyy_1, g_0_yyyyzzz_0_yyyy_0, g_0_yyyyzzz_0_yyyy_1, g_0_yyyyzzzz_0_xxxx_0, g_0_yyyyzzzz_0_xxxy_0, g_0_yyyyzzzz_0_xxxz_0, g_0_yyyyzzzz_0_xxyy_0, g_0_yyyyzzzz_0_xxyz_0, g_0_yyyyzzzz_0_xxzz_0, g_0_yyyyzzzz_0_xyyy_0, g_0_yyyyzzzz_0_xyyz_0, g_0_yyyyzzzz_0_xyzz_0, g_0_yyyyzzzz_0_xzzz_0, g_0_yyyyzzzz_0_yyyy_0, g_0_yyyyzzzz_0_yyyz_0, g_0_yyyyzzzz_0_yyzz_0, g_0_yyyyzzzz_0_yzzz_0, g_0_yyyyzzzz_0_zzzz_0, g_0_yyyzzzz_0_xxxx_0, g_0_yyyzzzz_0_xxxx_1, g_0_yyyzzzz_0_xxxz_0, g_0_yyyzzzz_0_xxxz_1, g_0_yyyzzzz_0_xxyz_0, g_0_yyyzzzz_0_xxyz_1, g_0_yyyzzzz_0_xxz_1, g_0_yyyzzzz_0_xxzz_0, g_0_yyyzzzz_0_xxzz_1, g_0_yyyzzzz_0_xyyz_0, g_0_yyyzzzz_0_xyyz_1, g_0_yyyzzzz_0_xyz_1, g_0_yyyzzzz_0_xyzz_0, g_0_yyyzzzz_0_xyzz_1, g_0_yyyzzzz_0_xzz_1, g_0_yyyzzzz_0_xzzz_0, g_0_yyyzzzz_0_xzzz_1, g_0_yyyzzzz_0_yyyz_0, g_0_yyyzzzz_0_yyyz_1, g_0_yyyzzzz_0_yyz_1, g_0_yyyzzzz_0_yyzz_0, g_0_yyyzzzz_0_yyzz_1, g_0_yyyzzzz_0_yzz_1, g_0_yyyzzzz_0_yzzz_0, g_0_yyyzzzz_0_yzzz_1, g_0_yyyzzzz_0_zzz_1, g_0_yyyzzzz_0_zzzz_0, g_0_yyyzzzz_0_zzzz_1, g_0_yyzzzz_0_xxxx_0, g_0_yyzzzz_0_xxxx_1, g_0_yyzzzz_0_xxxz_0, g_0_yyzzzz_0_xxxz_1, g_0_yyzzzz_0_xxyz_0, g_0_yyzzzz_0_xxyz_1, g_0_yyzzzz_0_xxzz_0, g_0_yyzzzz_0_xxzz_1, g_0_yyzzzz_0_xyyz_0, g_0_yyzzzz_0_xyyz_1, g_0_yyzzzz_0_xyzz_0, g_0_yyzzzz_0_xyzz_1, g_0_yyzzzz_0_xzzz_0, g_0_yyzzzz_0_xzzz_1, g_0_yyzzzz_0_yyyz_0, g_0_yyzzzz_0_yyyz_1, g_0_yyzzzz_0_yyzz_0, g_0_yyzzzz_0_yyzz_1, g_0_yyzzzz_0_yzzz_0, g_0_yyzzzz_0_yzzz_1, g_0_yyzzzz_0_zzzz_0, g_0_yyzzzz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_xxxx_0[i] = 3.0 * g_0_yyzzzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxx_0[i] * pb_y + g_0_yyyzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxy_0[i] = 3.0 * g_0_yyyyzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxy_0[i] * pb_z + g_0_yyyyzzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxz_0[i] = 3.0 * g_0_yyzzzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxz_0[i] * pb_y + g_0_yyyzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxyy_0[i] * pb_z + g_0_yyyyzzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyz_0[i] * pb_y + g_0_yyyzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxzz_0[i] * pb_y + g_0_yyyzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyy_0[i] = 3.0 * g_0_yyyyzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xyyy_0[i] * pb_z + g_0_yyyyzzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xyyz_0[i] = 3.0 * g_0_yyzzzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyz_0[i] * pb_y + g_0_yyyzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyzz_0[i] = 3.0 * g_0_yyzzzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzz_0[i] * pb_y + g_0_yyyzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xzzz_0[i] = 3.0 * g_0_yyzzzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzzz_0[i] * pb_y + g_0_yyyzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyy_0[i] = 3.0 * g_0_yyyyzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_yyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_yyyy_0[i] * pb_z + g_0_yyyyzzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_yyyz_0[i] = 3.0 * g_0_yyzzzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyz_0[i] * pb_y + g_0_yyyzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyzz_0[i] = 3.0 * g_0_yyzzzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyzz_0[i] * pb_y + g_0_yyyzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yzzz_0[i] = 3.0 * g_0_yyzzzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yzzz_0[i] * pb_y + g_0_yyyzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_zzzz_0[i] = 3.0 * g_0_yyzzzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzzz_0[i] * pb_y + g_0_yyyzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 615-630 components of targeted buffer : SLSG

    auto g_0_yyyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 615);

    auto g_0_yyyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 616);

    auto g_0_yyyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 617);

    auto g_0_yyyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 618);

    auto g_0_yyyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 619);

    auto g_0_yyyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 620);

    auto g_0_yyyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 621);

    auto g_0_yyyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 622);

    auto g_0_yyyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 623);

    auto g_0_yyyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 624);

    auto g_0_yyyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 625);

    auto g_0_yyyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 626);

    auto g_0_yyyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 627);

    auto g_0_yyyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 628);

    auto g_0_yyyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 629);

    #pragma omp simd aligned(g_0_yyyzzz_0_xxxy_0, g_0_yyyzzz_0_xxxy_1, g_0_yyyzzz_0_xxyy_0, g_0_yyyzzz_0_xxyy_1, g_0_yyyzzz_0_xyyy_0, g_0_yyyzzz_0_xyyy_1, g_0_yyyzzz_0_yyyy_0, g_0_yyyzzz_0_yyyy_1, g_0_yyyzzzz_0_xxxy_0, g_0_yyyzzzz_0_xxxy_1, g_0_yyyzzzz_0_xxyy_0, g_0_yyyzzzz_0_xxyy_1, g_0_yyyzzzz_0_xyyy_0, g_0_yyyzzzz_0_xyyy_1, g_0_yyyzzzz_0_yyyy_0, g_0_yyyzzzz_0_yyyy_1, g_0_yyyzzzzz_0_xxxx_0, g_0_yyyzzzzz_0_xxxy_0, g_0_yyyzzzzz_0_xxxz_0, g_0_yyyzzzzz_0_xxyy_0, g_0_yyyzzzzz_0_xxyz_0, g_0_yyyzzzzz_0_xxzz_0, g_0_yyyzzzzz_0_xyyy_0, g_0_yyyzzzzz_0_xyyz_0, g_0_yyyzzzzz_0_xyzz_0, g_0_yyyzzzzz_0_xzzz_0, g_0_yyyzzzzz_0_yyyy_0, g_0_yyyzzzzz_0_yyyz_0, g_0_yyyzzzzz_0_yyzz_0, g_0_yyyzzzzz_0_yzzz_0, g_0_yyyzzzzz_0_zzzz_0, g_0_yyzzzzz_0_xxxx_0, g_0_yyzzzzz_0_xxxx_1, g_0_yyzzzzz_0_xxxz_0, g_0_yyzzzzz_0_xxxz_1, g_0_yyzzzzz_0_xxyz_0, g_0_yyzzzzz_0_xxyz_1, g_0_yyzzzzz_0_xxz_1, g_0_yyzzzzz_0_xxzz_0, g_0_yyzzzzz_0_xxzz_1, g_0_yyzzzzz_0_xyyz_0, g_0_yyzzzzz_0_xyyz_1, g_0_yyzzzzz_0_xyz_1, g_0_yyzzzzz_0_xyzz_0, g_0_yyzzzzz_0_xyzz_1, g_0_yyzzzzz_0_xzz_1, g_0_yyzzzzz_0_xzzz_0, g_0_yyzzzzz_0_xzzz_1, g_0_yyzzzzz_0_yyyz_0, g_0_yyzzzzz_0_yyyz_1, g_0_yyzzzzz_0_yyz_1, g_0_yyzzzzz_0_yyzz_0, g_0_yyzzzzz_0_yyzz_1, g_0_yyzzzzz_0_yzz_1, g_0_yyzzzzz_0_yzzz_0, g_0_yyzzzzz_0_yzzz_1, g_0_yyzzzzz_0_zzz_1, g_0_yyzzzzz_0_zzzz_0, g_0_yyzzzzz_0_zzzz_1, g_0_yzzzzz_0_xxxx_0, g_0_yzzzzz_0_xxxx_1, g_0_yzzzzz_0_xxxz_0, g_0_yzzzzz_0_xxxz_1, g_0_yzzzzz_0_xxyz_0, g_0_yzzzzz_0_xxyz_1, g_0_yzzzzz_0_xxzz_0, g_0_yzzzzz_0_xxzz_1, g_0_yzzzzz_0_xyyz_0, g_0_yzzzzz_0_xyyz_1, g_0_yzzzzz_0_xyzz_0, g_0_yzzzzz_0_xyzz_1, g_0_yzzzzz_0_xzzz_0, g_0_yzzzzz_0_xzzz_1, g_0_yzzzzz_0_yyyz_0, g_0_yzzzzz_0_yyyz_1, g_0_yzzzzz_0_yyzz_0, g_0_yzzzzz_0_yyzz_1, g_0_yzzzzz_0_yzzz_0, g_0_yzzzzz_0_yzzz_1, g_0_yzzzzz_0_zzzz_0, g_0_yzzzzz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_xxxx_0[i] = 2.0 * g_0_yzzzzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxx_0[i] * pb_y + g_0_yyzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxy_0[i] = 4.0 * g_0_yyyzzz_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxy_0[i] * pb_z + g_0_yyyzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxz_0[i] = 2.0 * g_0_yzzzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxz_0[i] * pb_y + g_0_yyzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyy_0[i] = 4.0 * g_0_yyyzzz_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxyy_0[i] * pb_z + g_0_yyyzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxyz_0[i] = 2.0 * g_0_yzzzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyz_0[i] * pb_y + g_0_yyzzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxzz_0[i] = 2.0 * g_0_yzzzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxzz_0[i] * pb_y + g_0_yyzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyy_0[i] = 4.0 * g_0_yyyzzz_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xyyy_0[i] * pb_z + g_0_yyyzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyz_0[i] * pb_y + g_0_yyzzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzz_0[i] * pb_y + g_0_yyzzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xzzz_0[i] = 2.0 * g_0_yzzzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzzz_0[i] * pb_y + g_0_yyzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyy_0[i] = 4.0 * g_0_yyyzzz_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_yyyy_0[i] * pb_z + g_0_yyyzzzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_yyyz_0[i] = 2.0 * g_0_yzzzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyz_0[i] * pb_y + g_0_yyzzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyzz_0[i] = 2.0 * g_0_yzzzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyzz_0[i] * pb_y + g_0_yyzzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yzzz_0[i] = 2.0 * g_0_yzzzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yzzz_0[i] * pb_y + g_0_yyzzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_zzzz_0[i] = 2.0 * g_0_yzzzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzzz_0[i] * pb_y + g_0_yyzzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 630-645 components of targeted buffer : SLSG

    auto g_0_yyzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 630);

    auto g_0_yyzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 631);

    auto g_0_yyzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 632);

    auto g_0_yyzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 633);

    auto g_0_yyzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 634);

    auto g_0_yyzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 635);

    auto g_0_yyzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 636);

    auto g_0_yyzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 637);

    auto g_0_yyzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 638);

    auto g_0_yyzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 639);

    auto g_0_yyzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 640);

    auto g_0_yyzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 641);

    auto g_0_yyzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 642);

    auto g_0_yyzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 643);

    auto g_0_yyzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 644);

    #pragma omp simd aligned(g_0_yyzzzz_0_xxxy_0, g_0_yyzzzz_0_xxxy_1, g_0_yyzzzz_0_xxyy_0, g_0_yyzzzz_0_xxyy_1, g_0_yyzzzz_0_xyyy_0, g_0_yyzzzz_0_xyyy_1, g_0_yyzzzz_0_yyyy_0, g_0_yyzzzz_0_yyyy_1, g_0_yyzzzzz_0_xxxy_0, g_0_yyzzzzz_0_xxxy_1, g_0_yyzzzzz_0_xxyy_0, g_0_yyzzzzz_0_xxyy_1, g_0_yyzzzzz_0_xyyy_0, g_0_yyzzzzz_0_xyyy_1, g_0_yyzzzzz_0_yyyy_0, g_0_yyzzzzz_0_yyyy_1, g_0_yyzzzzzz_0_xxxx_0, g_0_yyzzzzzz_0_xxxy_0, g_0_yyzzzzzz_0_xxxz_0, g_0_yyzzzzzz_0_xxyy_0, g_0_yyzzzzzz_0_xxyz_0, g_0_yyzzzzzz_0_xxzz_0, g_0_yyzzzzzz_0_xyyy_0, g_0_yyzzzzzz_0_xyyz_0, g_0_yyzzzzzz_0_xyzz_0, g_0_yyzzzzzz_0_xzzz_0, g_0_yyzzzzzz_0_yyyy_0, g_0_yyzzzzzz_0_yyyz_0, g_0_yyzzzzzz_0_yyzz_0, g_0_yyzzzzzz_0_yzzz_0, g_0_yyzzzzzz_0_zzzz_0, g_0_yzzzzzz_0_xxxx_0, g_0_yzzzzzz_0_xxxx_1, g_0_yzzzzzz_0_xxxz_0, g_0_yzzzzzz_0_xxxz_1, g_0_yzzzzzz_0_xxyz_0, g_0_yzzzzzz_0_xxyz_1, g_0_yzzzzzz_0_xxz_1, g_0_yzzzzzz_0_xxzz_0, g_0_yzzzzzz_0_xxzz_1, g_0_yzzzzzz_0_xyyz_0, g_0_yzzzzzz_0_xyyz_1, g_0_yzzzzzz_0_xyz_1, g_0_yzzzzzz_0_xyzz_0, g_0_yzzzzzz_0_xyzz_1, g_0_yzzzzzz_0_xzz_1, g_0_yzzzzzz_0_xzzz_0, g_0_yzzzzzz_0_xzzz_1, g_0_yzzzzzz_0_yyyz_0, g_0_yzzzzzz_0_yyyz_1, g_0_yzzzzzz_0_yyz_1, g_0_yzzzzzz_0_yyzz_0, g_0_yzzzzzz_0_yyzz_1, g_0_yzzzzzz_0_yzz_1, g_0_yzzzzzz_0_yzzz_0, g_0_yzzzzzz_0_yzzz_1, g_0_yzzzzzz_0_zzz_1, g_0_yzzzzzz_0_zzzz_0, g_0_yzzzzzz_0_zzzz_1, g_0_zzzzzz_0_xxxx_0, g_0_zzzzzz_0_xxxx_1, g_0_zzzzzz_0_xxxz_0, g_0_zzzzzz_0_xxxz_1, g_0_zzzzzz_0_xxyz_0, g_0_zzzzzz_0_xxyz_1, g_0_zzzzzz_0_xxzz_0, g_0_zzzzzz_0_xxzz_1, g_0_zzzzzz_0_xyyz_0, g_0_zzzzzz_0_xyyz_1, g_0_zzzzzz_0_xyzz_0, g_0_zzzzzz_0_xyzz_1, g_0_zzzzzz_0_xzzz_0, g_0_zzzzzz_0_xzzz_1, g_0_zzzzzz_0_yyyz_0, g_0_zzzzzz_0_yyyz_1, g_0_zzzzzz_0_yyzz_0, g_0_zzzzzz_0_yyzz_1, g_0_zzzzzz_0_yzzz_0, g_0_zzzzzz_0_yzzz_1, g_0_zzzzzz_0_zzzz_0, g_0_zzzzzz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_xxxx_0[i] = g_0_zzzzzz_0_xxxx_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxx_0[i] * pb_y + g_0_yzzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxy_0[i] = 5.0 * g_0_yyzzzz_0_xxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxy_0[i] * pb_z + g_0_yyzzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxz_0[i] = g_0_zzzzzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxz_0[i] * pb_y + g_0_yzzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyy_0[i] = 5.0 * g_0_yyzzzz_0_xxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxyy_0[i] * pb_z + g_0_yyzzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxyz_0[i] = g_0_zzzzzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyz_0[i] * pb_y + g_0_yzzzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxzz_0[i] = g_0_zzzzzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzz_0[i] * pb_y + g_0_yzzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyy_0[i] = 5.0 * g_0_yyzzzz_0_xyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xyyy_0[i] * pb_z + g_0_yyzzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xyyz_0[i] = g_0_zzzzzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyz_0[i] * pb_y + g_0_yzzzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyzz_0[i] = g_0_zzzzzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzz_0[i] * pb_y + g_0_yzzzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xzzz_0[i] = g_0_zzzzzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzz_0[i] * pb_y + g_0_yzzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyy_0[i] = 5.0 * g_0_yyzzzz_0_yyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_yyyy_0[i] * pb_z + g_0_yyzzzzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_yyyz_0[i] = g_0_zzzzzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyz_0[i] * pb_y + g_0_yzzzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyzz_0[i] = g_0_zzzzzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyzz_0[i] * pb_y + g_0_yzzzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yzzz_0[i] = g_0_zzzzzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yzzz_0[i] * pb_y + g_0_yzzzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_zzzz_0[i] = g_0_zzzzzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzz_0[i] * pb_y + g_0_yzzzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 645-660 components of targeted buffer : SLSG

    auto g_0_yzzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 645);

    auto g_0_yzzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 646);

    auto g_0_yzzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 647);

    auto g_0_yzzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 648);

    auto g_0_yzzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 649);

    auto g_0_yzzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 650);

    auto g_0_yzzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 651);

    auto g_0_yzzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 652);

    auto g_0_yzzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 653);

    auto g_0_yzzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 654);

    auto g_0_yzzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 655);

    auto g_0_yzzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 656);

    auto g_0_yzzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 657);

    auto g_0_yzzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 658);

    auto g_0_yzzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 659);

    #pragma omp simd aligned(g_0_yzzzzzzz_0_xxxx_0, g_0_yzzzzzzz_0_xxxy_0, g_0_yzzzzzzz_0_xxxz_0, g_0_yzzzzzzz_0_xxyy_0, g_0_yzzzzzzz_0_xxyz_0, g_0_yzzzzzzz_0_xxzz_0, g_0_yzzzzzzz_0_xyyy_0, g_0_yzzzzzzz_0_xyyz_0, g_0_yzzzzzzz_0_xyzz_0, g_0_yzzzzzzz_0_xzzz_0, g_0_yzzzzzzz_0_yyyy_0, g_0_yzzzzzzz_0_yyyz_0, g_0_yzzzzzzz_0_yyzz_0, g_0_yzzzzzzz_0_yzzz_0, g_0_yzzzzzzz_0_zzzz_0, g_0_zzzzzzz_0_xxx_1, g_0_zzzzzzz_0_xxxx_0, g_0_zzzzzzz_0_xxxx_1, g_0_zzzzzzz_0_xxxy_0, g_0_zzzzzzz_0_xxxy_1, g_0_zzzzzzz_0_xxxz_0, g_0_zzzzzzz_0_xxxz_1, g_0_zzzzzzz_0_xxy_1, g_0_zzzzzzz_0_xxyy_0, g_0_zzzzzzz_0_xxyy_1, g_0_zzzzzzz_0_xxyz_0, g_0_zzzzzzz_0_xxyz_1, g_0_zzzzzzz_0_xxz_1, g_0_zzzzzzz_0_xxzz_0, g_0_zzzzzzz_0_xxzz_1, g_0_zzzzzzz_0_xyy_1, g_0_zzzzzzz_0_xyyy_0, g_0_zzzzzzz_0_xyyy_1, g_0_zzzzzzz_0_xyyz_0, g_0_zzzzzzz_0_xyyz_1, g_0_zzzzzzz_0_xyz_1, g_0_zzzzzzz_0_xyzz_0, g_0_zzzzzzz_0_xyzz_1, g_0_zzzzzzz_0_xzz_1, g_0_zzzzzzz_0_xzzz_0, g_0_zzzzzzz_0_xzzz_1, g_0_zzzzzzz_0_yyy_1, g_0_zzzzzzz_0_yyyy_0, g_0_zzzzzzz_0_yyyy_1, g_0_zzzzzzz_0_yyyz_0, g_0_zzzzzzz_0_yyyz_1, g_0_zzzzzzz_0_yyz_1, g_0_zzzzzzz_0_yyzz_0, g_0_zzzzzzz_0_yyzz_1, g_0_zzzzzzz_0_yzz_1, g_0_zzzzzzz_0_yzzz_0, g_0_zzzzzzz_0_yzzz_1, g_0_zzzzzzz_0_zzz_1, g_0_zzzzzzz_0_zzzz_0, g_0_zzzzzzz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_xxxx_0[i] = g_0_zzzzzzz_0_xxxx_0[i] * pb_y + g_0_zzzzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxy_0[i] = g_0_zzzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxy_0[i] * pb_y + g_0_zzzzzzz_0_xxxy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxz_0[i] = g_0_zzzzzzz_0_xxxz_0[i] * pb_y + g_0_zzzzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyy_0[i] * pb_y + g_0_zzzzzzz_0_xxyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyz_0[i] = g_0_zzzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyz_0[i] * pb_y + g_0_zzzzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxzz_0[i] = g_0_zzzzzzz_0_xxzz_0[i] * pb_y + g_0_zzzzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyy_0[i] = 3.0 * g_0_zzzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyy_0[i] * pb_y + g_0_zzzzzzz_0_xyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyz_0[i] = 2.0 * g_0_zzzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyz_0[i] * pb_y + g_0_zzzzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyzz_0[i] = g_0_zzzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzz_0[i] * pb_y + g_0_zzzzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xzzz_0[i] = g_0_zzzzzzz_0_xzzz_0[i] * pb_y + g_0_zzzzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyy_0[i] = 4.0 * g_0_zzzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyy_0[i] * pb_y + g_0_zzzzzzz_0_yyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyz_0[i] = 3.0 * g_0_zzzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyz_0[i] * pb_y + g_0_zzzzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyzz_0[i] = 2.0 * g_0_zzzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzz_0[i] * pb_y + g_0_zzzzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yzzz_0[i] = g_0_zzzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzz_0[i] * pb_y + g_0_zzzzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_zzzz_0[i] = g_0_zzzzzzz_0_zzzz_0[i] * pb_y + g_0_zzzzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 660-675 components of targeted buffer : SLSG

    auto g_0_zzzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_slsg + 660);

    auto g_0_zzzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_slsg + 661);

    auto g_0_zzzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_slsg + 662);

    auto g_0_zzzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_slsg + 663);

    auto g_0_zzzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_slsg + 664);

    auto g_0_zzzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_slsg + 665);

    auto g_0_zzzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_slsg + 666);

    auto g_0_zzzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_slsg + 667);

    auto g_0_zzzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_slsg + 668);

    auto g_0_zzzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_slsg + 669);

    auto g_0_zzzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_slsg + 670);

    auto g_0_zzzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_slsg + 671);

    auto g_0_zzzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_slsg + 672);

    auto g_0_zzzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_slsg + 673);

    auto g_0_zzzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_slsg + 674);

    #pragma omp simd aligned(g_0_zzzzzz_0_xxxx_0, g_0_zzzzzz_0_xxxx_1, g_0_zzzzzz_0_xxxy_0, g_0_zzzzzz_0_xxxy_1, g_0_zzzzzz_0_xxxz_0, g_0_zzzzzz_0_xxxz_1, g_0_zzzzzz_0_xxyy_0, g_0_zzzzzz_0_xxyy_1, g_0_zzzzzz_0_xxyz_0, g_0_zzzzzz_0_xxyz_1, g_0_zzzzzz_0_xxzz_0, g_0_zzzzzz_0_xxzz_1, g_0_zzzzzz_0_xyyy_0, g_0_zzzzzz_0_xyyy_1, g_0_zzzzzz_0_xyyz_0, g_0_zzzzzz_0_xyyz_1, g_0_zzzzzz_0_xyzz_0, g_0_zzzzzz_0_xyzz_1, g_0_zzzzzz_0_xzzz_0, g_0_zzzzzz_0_xzzz_1, g_0_zzzzzz_0_yyyy_0, g_0_zzzzzz_0_yyyy_1, g_0_zzzzzz_0_yyyz_0, g_0_zzzzzz_0_yyyz_1, g_0_zzzzzz_0_yyzz_0, g_0_zzzzzz_0_yyzz_1, g_0_zzzzzz_0_yzzz_0, g_0_zzzzzz_0_yzzz_1, g_0_zzzzzz_0_zzzz_0, g_0_zzzzzz_0_zzzz_1, g_0_zzzzzzz_0_xxx_1, g_0_zzzzzzz_0_xxxx_0, g_0_zzzzzzz_0_xxxx_1, g_0_zzzzzzz_0_xxxy_0, g_0_zzzzzzz_0_xxxy_1, g_0_zzzzzzz_0_xxxz_0, g_0_zzzzzzz_0_xxxz_1, g_0_zzzzzzz_0_xxy_1, g_0_zzzzzzz_0_xxyy_0, g_0_zzzzzzz_0_xxyy_1, g_0_zzzzzzz_0_xxyz_0, g_0_zzzzzzz_0_xxyz_1, g_0_zzzzzzz_0_xxz_1, g_0_zzzzzzz_0_xxzz_0, g_0_zzzzzzz_0_xxzz_1, g_0_zzzzzzz_0_xyy_1, g_0_zzzzzzz_0_xyyy_0, g_0_zzzzzzz_0_xyyy_1, g_0_zzzzzzz_0_xyyz_0, g_0_zzzzzzz_0_xyyz_1, g_0_zzzzzzz_0_xyz_1, g_0_zzzzzzz_0_xyzz_0, g_0_zzzzzzz_0_xyzz_1, g_0_zzzzzzz_0_xzz_1, g_0_zzzzzzz_0_xzzz_0, g_0_zzzzzzz_0_xzzz_1, g_0_zzzzzzz_0_yyy_1, g_0_zzzzzzz_0_yyyy_0, g_0_zzzzzzz_0_yyyy_1, g_0_zzzzzzz_0_yyyz_0, g_0_zzzzzzz_0_yyyz_1, g_0_zzzzzzz_0_yyz_1, g_0_zzzzzzz_0_yyzz_0, g_0_zzzzzzz_0_yyzz_1, g_0_zzzzzzz_0_yzz_1, g_0_zzzzzzz_0_yzzz_0, g_0_zzzzzzz_0_yzzz_1, g_0_zzzzzzz_0_zzz_1, g_0_zzzzzzz_0_zzzz_0, g_0_zzzzzzz_0_zzzz_1, g_0_zzzzzzzz_0_xxxx_0, g_0_zzzzzzzz_0_xxxy_0, g_0_zzzzzzzz_0_xxxz_0, g_0_zzzzzzzz_0_xxyy_0, g_0_zzzzzzzz_0_xxyz_0, g_0_zzzzzzzz_0_xxzz_0, g_0_zzzzzzzz_0_xyyy_0, g_0_zzzzzzzz_0_xyyz_0, g_0_zzzzzzzz_0_xyzz_0, g_0_zzzzzzzz_0_xzzz_0, g_0_zzzzzzzz_0_yyyy_0, g_0_zzzzzzzz_0_yyyz_0, g_0_zzzzzzzz_0_yyzz_0, g_0_zzzzzzzz_0_yzzz_0, g_0_zzzzzzzz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_xxxx_0[i] = 7.0 * g_0_zzzzzz_0_xxxx_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxx_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxx_0[i] * pb_z + g_0_zzzzzzz_0_xxxx_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxy_0[i] = 7.0 * g_0_zzzzzz_0_xxxy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxy_0[i] * pb_z + g_0_zzzzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxz_0[i] = 7.0 * g_0_zzzzzz_0_xxxz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxz_0[i] * pb_z + g_0_zzzzzzz_0_xxxz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyy_0[i] = 7.0 * g_0_zzzzzz_0_xxyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxyy_0[i] * pb_z + g_0_zzzzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyz_0[i] = 7.0 * g_0_zzzzzz_0_xxyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyz_0[i] * pb_z + g_0_zzzzzzz_0_xxyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxzz_0[i] = 7.0 * g_0_zzzzzz_0_xxzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzz_0[i] * pb_z + g_0_zzzzzzz_0_xxzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyy_0[i] = 7.0 * g_0_zzzzzz_0_xyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyyy_0[i] * pb_z + g_0_zzzzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyz_0[i] = 7.0 * g_0_zzzzzz_0_xyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyz_0[i] * pb_z + g_0_zzzzzzz_0_xyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyzz_0[i] = 7.0 * g_0_zzzzzz_0_xyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzz_0[i] * pb_z + g_0_zzzzzzz_0_xyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xzzz_0[i] = 7.0 * g_0_zzzzzz_0_xzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzz_0[i] * pb_z + g_0_zzzzzzz_0_xzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyy_0[i] = 7.0 * g_0_zzzzzz_0_yyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyyy_0[i] * pb_z + g_0_zzzzzzz_0_yyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyz_0[i] = 7.0 * g_0_zzzzzz_0_yyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyz_0[i] * pb_z + g_0_zzzzzzz_0_yyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyzz_0[i] = 7.0 * g_0_zzzzzz_0_yyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzz_0[i] * pb_z + g_0_zzzzzzz_0_yyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yzzz_0[i] = 7.0 * g_0_zzzzzz_0_yzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzz_0[i] * pb_z + g_0_zzzzzzz_0_yzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_zzzz_0[i] = 7.0 * g_0_zzzzzz_0_zzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_zzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_zzzz_0[i] * pb_z + g_0_zzzzzzz_0_zzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

