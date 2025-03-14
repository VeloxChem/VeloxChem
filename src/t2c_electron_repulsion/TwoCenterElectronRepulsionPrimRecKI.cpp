#include "TwoCenterElectronRepulsionPrimRecKI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ki(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ki,
                                const size_t idx_eri_0_hi,
                                const size_t idx_eri_1_hi,
                                const size_t idx_eri_1_ih,
                                const size_t idx_eri_1_ii,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : HI

    auto g_xxxxx_xxxxxx_0 = pbuffer.data(idx_eri_0_hi);

    auto g_xxxxx_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 1);

    auto g_xxxxx_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 2);

    auto g_xxxxx_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 3);

    auto g_xxxxx_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 4);

    auto g_xxxxx_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 5);

    auto g_xxxxx_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 6);

    auto g_xxxxx_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 7);

    auto g_xxxxx_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 8);

    auto g_xxxxx_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 9);

    auto g_xxxxx_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 10);

    auto g_xxxxx_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 11);

    auto g_xxxxx_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 12);

    auto g_xxxxx_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 13);

    auto g_xxxxx_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 14);

    auto g_xxxxx_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 15);

    auto g_xxxxx_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 16);

    auto g_xxxxx_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 17);

    auto g_xxxxx_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 18);

    auto g_xxxxx_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 19);

    auto g_xxxxx_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 20);

    auto g_xxxxx_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 21);

    auto g_xxxxx_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 22);

    auto g_xxxxx_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 23);

    auto g_xxxxx_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 24);

    auto g_xxxxx_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 25);

    auto g_xxxxx_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 26);

    auto g_xxxxx_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 27);

    auto g_xxxxy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 28);

    auto g_xxxxy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 30);

    auto g_xxxxy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 33);

    auto g_xxxxy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 37);

    auto g_xxxxy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 42);

    auto g_xxxxy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 48);

    auto g_xxxxz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 56);

    auto g_xxxxz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 57);

    auto g_xxxxz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 59);

    auto g_xxxxz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 62);

    auto g_xxxxz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 66);

    auto g_xxxxz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 71);

    auto g_xxxyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 84);

    auto g_xxxyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 85);

    auto g_xxxyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 86);

    auto g_xxxyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 87);

    auto g_xxxyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 88);

    auto g_xxxyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 89);

    auto g_xxxyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 90);

    auto g_xxxyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 91);

    auto g_xxxyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 92);

    auto g_xxxyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 93);

    auto g_xxxyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 94);

    auto g_xxxyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 95);

    auto g_xxxyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 96);

    auto g_xxxyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 97);

    auto g_xxxyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 98);

    auto g_xxxyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 99);

    auto g_xxxyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 100);

    auto g_xxxyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 101);

    auto g_xxxyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 102);

    auto g_xxxyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 103);

    auto g_xxxyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 104);

    auto g_xxxyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 105);

    auto g_xxxyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 106);

    auto g_xxxyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 107);

    auto g_xxxyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 108);

    auto g_xxxyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 109);

    auto g_xxxyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 110);

    auto g_xxxyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 111);

    auto g_xxxzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 140);

    auto g_xxxzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 141);

    auto g_xxxzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 142);

    auto g_xxxzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 143);

    auto g_xxxzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 144);

    auto g_xxxzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 145);

    auto g_xxxzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 146);

    auto g_xxxzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 147);

    auto g_xxxzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 148);

    auto g_xxxzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 149);

    auto g_xxxzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 150);

    auto g_xxxzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 151);

    auto g_xxxzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 152);

    auto g_xxxzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 153);

    auto g_xxxzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 154);

    auto g_xxxzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 155);

    auto g_xxxzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 156);

    auto g_xxxzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 157);

    auto g_xxxzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 158);

    auto g_xxxzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 159);

    auto g_xxxzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 160);

    auto g_xxxzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 161);

    auto g_xxxzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 162);

    auto g_xxxzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 163);

    auto g_xxxzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 164);

    auto g_xxxzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 165);

    auto g_xxxzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 166);

    auto g_xxxzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 167);

    auto g_xxyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 168);

    auto g_xxyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 169);

    auto g_xxyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 170);

    auto g_xxyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 171);

    auto g_xxyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 172);

    auto g_xxyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 173);

    auto g_xxyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 174);

    auto g_xxyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 175);

    auto g_xxyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 176);

    auto g_xxyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 177);

    auto g_xxyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 178);

    auto g_xxyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 179);

    auto g_xxyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 180);

    auto g_xxyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 181);

    auto g_xxyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 182);

    auto g_xxyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 183);

    auto g_xxyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 184);

    auto g_xxyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 185);

    auto g_xxyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 186);

    auto g_xxyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 187);

    auto g_xxyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 188);

    auto g_xxyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 189);

    auto g_xxyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 190);

    auto g_xxyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 191);

    auto g_xxyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 192);

    auto g_xxyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 193);

    auto g_xxyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 194);

    auto g_xxyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 195);

    auto g_xxyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 197);

    auto g_xxyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 199);

    auto g_xxyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 202);

    auto g_xxyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 206);

    auto g_xxyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 211);

    auto g_xxyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 224);

    auto g_xxyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 226);

    auto g_xxyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 229);

    auto g_xxyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 233);

    auto g_xxyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 238);

    auto g_xxyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 244);

    auto g_xxzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 252);

    auto g_xxzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 253);

    auto g_xxzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 254);

    auto g_xxzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 255);

    auto g_xxzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 256);

    auto g_xxzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 257);

    auto g_xxzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 258);

    auto g_xxzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 259);

    auto g_xxzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 260);

    auto g_xxzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 261);

    auto g_xxzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 262);

    auto g_xxzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 263);

    auto g_xxzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 264);

    auto g_xxzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 265);

    auto g_xxzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 266);

    auto g_xxzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 267);

    auto g_xxzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 268);

    auto g_xxzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 269);

    auto g_xxzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 270);

    auto g_xxzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 271);

    auto g_xxzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 272);

    auto g_xxzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 273);

    auto g_xxzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 274);

    auto g_xxzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 275);

    auto g_xxzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 276);

    auto g_xxzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 277);

    auto g_xxzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 278);

    auto g_xxzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 279);

    auto g_xyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 281);

    auto g_xyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 283);

    auto g_xyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 284);

    auto g_xyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 286);

    auto g_xyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 287);

    auto g_xyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 288);

    auto g_xyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 290);

    auto g_xyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 291);

    auto g_xyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 292);

    auto g_xyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 293);

    auto g_xyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 295);

    auto g_xyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 296);

    auto g_xyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 297);

    auto g_xyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 298);

    auto g_xyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 299);

    auto g_xyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 301);

    auto g_xyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 302);

    auto g_xyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 303);

    auto g_xyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 304);

    auto g_xyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 305);

    auto g_xyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 306);

    auto g_xyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 307);

    auto g_xyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 340);

    auto g_xyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 343);

    auto g_xyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 344);

    auto g_xyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 347);

    auto g_xyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 348);

    auto g_xyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 349);

    auto g_xyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 352);

    auto g_xyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 353);

    auto g_xyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 354);

    auto g_xyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 355);

    auto g_xyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 357);

    auto g_xyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 358);

    auto g_xyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 359);

    auto g_xyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 360);

    auto g_xyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 361);

    auto g_xyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 362);

    auto g_xyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 363);

    auto g_xzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 394);

    auto g_xzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 396);

    auto g_xzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 397);

    auto g_xzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 399);

    auto g_xzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 400);

    auto g_xzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 401);

    auto g_xzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 403);

    auto g_xzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 404);

    auto g_xzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 405);

    auto g_xzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 406);

    auto g_xzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 408);

    auto g_xzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 409);

    auto g_xzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 410);

    auto g_xzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 411);

    auto g_xzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 412);

    auto g_xzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 413);

    auto g_xzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 414);

    auto g_xzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 415);

    auto g_xzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 416);

    auto g_xzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 417);

    auto g_xzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 418);

    auto g_xzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 419);

    auto g_yyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 420);

    auto g_yyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 421);

    auto g_yyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 422);

    auto g_yyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 423);

    auto g_yyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 424);

    auto g_yyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 425);

    auto g_yyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 426);

    auto g_yyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 427);

    auto g_yyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 428);

    auto g_yyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 429);

    auto g_yyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 430);

    auto g_yyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 431);

    auto g_yyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 432);

    auto g_yyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 433);

    auto g_yyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 434);

    auto g_yyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 435);

    auto g_yyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 436);

    auto g_yyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 437);

    auto g_yyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 438);

    auto g_yyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 439);

    auto g_yyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 440);

    auto g_yyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 441);

    auto g_yyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 442);

    auto g_yyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 443);

    auto g_yyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 444);

    auto g_yyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 445);

    auto g_yyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 446);

    auto g_yyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 447);

    auto g_yyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 449);

    auto g_yyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 451);

    auto g_yyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 454);

    auto g_yyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 458);

    auto g_yyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 463);

    auto g_yyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 469);

    auto g_yyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 476);

    auto g_yyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 477);

    auto g_yyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 478);

    auto g_yyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 479);

    auto g_yyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 480);

    auto g_yyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 481);

    auto g_yyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 482);

    auto g_yyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 483);

    auto g_yyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 484);

    auto g_yyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 485);

    auto g_yyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 486);

    auto g_yyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 487);

    auto g_yyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 488);

    auto g_yyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 489);

    auto g_yyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 490);

    auto g_yyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 491);

    auto g_yyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 492);

    auto g_yyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 493);

    auto g_yyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 494);

    auto g_yyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 495);

    auto g_yyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 496);

    auto g_yyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 497);

    auto g_yyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 498);

    auto g_yyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 499);

    auto g_yyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 500);

    auto g_yyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 501);

    auto g_yyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 502);

    auto g_yyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 503);

    auto g_yyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 504);

    auto g_yyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 505);

    auto g_yyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 506);

    auto g_yyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 507);

    auto g_yyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 508);

    auto g_yyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 509);

    auto g_yyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 510);

    auto g_yyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 511);

    auto g_yyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 512);

    auto g_yyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 513);

    auto g_yyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 514);

    auto g_yyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 515);

    auto g_yyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 516);

    auto g_yyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 517);

    auto g_yyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 518);

    auto g_yyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 519);

    auto g_yyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 520);

    auto g_yyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 521);

    auto g_yyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 522);

    auto g_yyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 523);

    auto g_yyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 524);

    auto g_yyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 525);

    auto g_yyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 526);

    auto g_yyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 527);

    auto g_yyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 528);

    auto g_yyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 529);

    auto g_yyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 530);

    auto g_yyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 531);

    auto g_yzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 532);

    auto g_yzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 534);

    auto g_yzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 536);

    auto g_yzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 537);

    auto g_yzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 539);

    auto g_yzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 540);

    auto g_yzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 541);

    auto g_yzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 543);

    auto g_yzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 544);

    auto g_yzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 545);

    auto g_yzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 546);

    auto g_yzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 548);

    auto g_yzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 549);

    auto g_yzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 550);

    auto g_yzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 551);

    auto g_yzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 552);

    auto g_yzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 554);

    auto g_yzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 555);

    auto g_yzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 556);

    auto g_yzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 557);

    auto g_yzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 558);

    auto g_yzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 559);

    auto g_zzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 560);

    auto g_zzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 561);

    auto g_zzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 562);

    auto g_zzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 563);

    auto g_zzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 564);

    auto g_zzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 565);

    auto g_zzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 566);

    auto g_zzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 567);

    auto g_zzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 568);

    auto g_zzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 569);

    auto g_zzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 570);

    auto g_zzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 571);

    auto g_zzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 572);

    auto g_zzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 573);

    auto g_zzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 574);

    auto g_zzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 575);

    auto g_zzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 576);

    auto g_zzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 577);

    auto g_zzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 578);

    auto g_zzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 579);

    auto g_zzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 580);

    auto g_zzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 581);

    auto g_zzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 582);

    auto g_zzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 583);

    auto g_zzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 584);

    auto g_zzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 585);

    auto g_zzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 586);

    auto g_zzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 587);

    // Set up components of auxiliary buffer : HI

    auto g_xxxxx_xxxxxx_1 = pbuffer.data(idx_eri_1_hi);

    auto g_xxxxx_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 1);

    auto g_xxxxx_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 2);

    auto g_xxxxx_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 3);

    auto g_xxxxx_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 4);

    auto g_xxxxx_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 5);

    auto g_xxxxx_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 6);

    auto g_xxxxx_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 7);

    auto g_xxxxx_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 8);

    auto g_xxxxx_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 9);

    auto g_xxxxx_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 10);

    auto g_xxxxx_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 11);

    auto g_xxxxx_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 12);

    auto g_xxxxx_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 13);

    auto g_xxxxx_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 14);

    auto g_xxxxx_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 15);

    auto g_xxxxx_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 16);

    auto g_xxxxx_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 17);

    auto g_xxxxx_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 18);

    auto g_xxxxx_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 19);

    auto g_xxxxx_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 20);

    auto g_xxxxx_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 21);

    auto g_xxxxx_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 22);

    auto g_xxxxx_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 23);

    auto g_xxxxx_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 24);

    auto g_xxxxx_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 25);

    auto g_xxxxx_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 26);

    auto g_xxxxx_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 27);

    auto g_xxxxy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 28);

    auto g_xxxxy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 30);

    auto g_xxxxy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 33);

    auto g_xxxxy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 37);

    auto g_xxxxy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 42);

    auto g_xxxxy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 48);

    auto g_xxxxz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 56);

    auto g_xxxxz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 57);

    auto g_xxxxz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 59);

    auto g_xxxxz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 62);

    auto g_xxxxz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 66);

    auto g_xxxxz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 71);

    auto g_xxxyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 84);

    auto g_xxxyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 85);

    auto g_xxxyy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 86);

    auto g_xxxyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 87);

    auto g_xxxyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 88);

    auto g_xxxyy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 89);

    auto g_xxxyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 90);

    auto g_xxxyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 91);

    auto g_xxxyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 92);

    auto g_xxxyy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 93);

    auto g_xxxyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 94);

    auto g_xxxyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 95);

    auto g_xxxyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 96);

    auto g_xxxyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 97);

    auto g_xxxyy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 98);

    auto g_xxxyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 99);

    auto g_xxxyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 100);

    auto g_xxxyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 101);

    auto g_xxxyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 102);

    auto g_xxxyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 103);

    auto g_xxxyy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 104);

    auto g_xxxyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 105);

    auto g_xxxyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 106);

    auto g_xxxyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 107);

    auto g_xxxyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 108);

    auto g_xxxyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 109);

    auto g_xxxyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 110);

    auto g_xxxyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 111);

    auto g_xxxzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 140);

    auto g_xxxzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 141);

    auto g_xxxzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 142);

    auto g_xxxzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 143);

    auto g_xxxzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 144);

    auto g_xxxzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 145);

    auto g_xxxzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 146);

    auto g_xxxzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 147);

    auto g_xxxzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 148);

    auto g_xxxzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 149);

    auto g_xxxzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 150);

    auto g_xxxzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 151);

    auto g_xxxzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 152);

    auto g_xxxzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 153);

    auto g_xxxzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 154);

    auto g_xxxzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 155);

    auto g_xxxzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 156);

    auto g_xxxzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 157);

    auto g_xxxzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 158);

    auto g_xxxzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 159);

    auto g_xxxzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 160);

    auto g_xxxzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 161);

    auto g_xxxzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 162);

    auto g_xxxzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 163);

    auto g_xxxzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 164);

    auto g_xxxzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 165);

    auto g_xxxzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 166);

    auto g_xxxzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 167);

    auto g_xxyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 168);

    auto g_xxyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 169);

    auto g_xxyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 170);

    auto g_xxyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 171);

    auto g_xxyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 172);

    auto g_xxyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 173);

    auto g_xxyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 174);

    auto g_xxyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 175);

    auto g_xxyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 176);

    auto g_xxyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 177);

    auto g_xxyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 178);

    auto g_xxyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 179);

    auto g_xxyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 180);

    auto g_xxyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 181);

    auto g_xxyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 182);

    auto g_xxyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 183);

    auto g_xxyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 184);

    auto g_xxyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 185);

    auto g_xxyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 186);

    auto g_xxyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 187);

    auto g_xxyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 188);

    auto g_xxyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 189);

    auto g_xxyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 190);

    auto g_xxyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 191);

    auto g_xxyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 192);

    auto g_xxyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 193);

    auto g_xxyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 194);

    auto g_xxyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 195);

    auto g_xxyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 197);

    auto g_xxyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 199);

    auto g_xxyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 202);

    auto g_xxyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 206);

    auto g_xxyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 211);

    auto g_xxyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 224);

    auto g_xxyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 226);

    auto g_xxyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 229);

    auto g_xxyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 233);

    auto g_xxyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 238);

    auto g_xxyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 244);

    auto g_xxzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 252);

    auto g_xxzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 253);

    auto g_xxzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 254);

    auto g_xxzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 255);

    auto g_xxzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 256);

    auto g_xxzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 257);

    auto g_xxzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 258);

    auto g_xxzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 259);

    auto g_xxzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 260);

    auto g_xxzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 261);

    auto g_xxzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 262);

    auto g_xxzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 263);

    auto g_xxzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 264);

    auto g_xxzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 265);

    auto g_xxzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 266);

    auto g_xxzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 267);

    auto g_xxzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 268);

    auto g_xxzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 269);

    auto g_xxzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 270);

    auto g_xxzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 271);

    auto g_xxzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 272);

    auto g_xxzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 273);

    auto g_xxzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 274);

    auto g_xxzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 275);

    auto g_xxzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 276);

    auto g_xxzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 277);

    auto g_xxzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 278);

    auto g_xxzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 279);

    auto g_xyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 281);

    auto g_xyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 283);

    auto g_xyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 284);

    auto g_xyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 286);

    auto g_xyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 287);

    auto g_xyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 288);

    auto g_xyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 290);

    auto g_xyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 291);

    auto g_xyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 292);

    auto g_xyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 293);

    auto g_xyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 295);

    auto g_xyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 296);

    auto g_xyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 297);

    auto g_xyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 298);

    auto g_xyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 299);

    auto g_xyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 301);

    auto g_xyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 302);

    auto g_xyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 303);

    auto g_xyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 304);

    auto g_xyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 305);

    auto g_xyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 306);

    auto g_xyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 307);

    auto g_xyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 340);

    auto g_xyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 343);

    auto g_xyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 344);

    auto g_xyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 347);

    auto g_xyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 348);

    auto g_xyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 349);

    auto g_xyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 352);

    auto g_xyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 353);

    auto g_xyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 354);

    auto g_xyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 355);

    auto g_xyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 357);

    auto g_xyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 358);

    auto g_xyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 359);

    auto g_xyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 360);

    auto g_xyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 361);

    auto g_xyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 362);

    auto g_xyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 363);

    auto g_xzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 394);

    auto g_xzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 396);

    auto g_xzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 397);

    auto g_xzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 399);

    auto g_xzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 400);

    auto g_xzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 401);

    auto g_xzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 403);

    auto g_xzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 404);

    auto g_xzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 405);

    auto g_xzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 406);

    auto g_xzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 408);

    auto g_xzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 409);

    auto g_xzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 410);

    auto g_xzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 411);

    auto g_xzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 412);

    auto g_xzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 413);

    auto g_xzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 414);

    auto g_xzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 415);

    auto g_xzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 416);

    auto g_xzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 417);

    auto g_xzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 418);

    auto g_xzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 419);

    auto g_yyyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 420);

    auto g_yyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 421);

    auto g_yyyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 422);

    auto g_yyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 423);

    auto g_yyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 424);

    auto g_yyyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 425);

    auto g_yyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 426);

    auto g_yyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 427);

    auto g_yyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 428);

    auto g_yyyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 429);

    auto g_yyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 430);

    auto g_yyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 431);

    auto g_yyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 432);

    auto g_yyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 433);

    auto g_yyyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 434);

    auto g_yyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 435);

    auto g_yyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 436);

    auto g_yyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 437);

    auto g_yyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 438);

    auto g_yyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 439);

    auto g_yyyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 440);

    auto g_yyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 441);

    auto g_yyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 442);

    auto g_yyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 443);

    auto g_yyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 444);

    auto g_yyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 445);

    auto g_yyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 446);

    auto g_yyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 447);

    auto g_yyyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 449);

    auto g_yyyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 451);

    auto g_yyyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 454);

    auto g_yyyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 458);

    auto g_yyyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 463);

    auto g_yyyyz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 469);

    auto g_yyyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 476);

    auto g_yyyzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 477);

    auto g_yyyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 478);

    auto g_yyyzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 479);

    auto g_yyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 480);

    auto g_yyyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 481);

    auto g_yyyzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 482);

    auto g_yyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 483);

    auto g_yyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 484);

    auto g_yyyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 485);

    auto g_yyyzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 486);

    auto g_yyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 487);

    auto g_yyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 488);

    auto g_yyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 489);

    auto g_yyyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 490);

    auto g_yyyzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 491);

    auto g_yyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 492);

    auto g_yyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 493);

    auto g_yyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 494);

    auto g_yyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 495);

    auto g_yyyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 496);

    auto g_yyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 497);

    auto g_yyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 498);

    auto g_yyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 499);

    auto g_yyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 500);

    auto g_yyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 501);

    auto g_yyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 502);

    auto g_yyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 503);

    auto g_yyzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 504);

    auto g_yyzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 505);

    auto g_yyzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 506);

    auto g_yyzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 507);

    auto g_yyzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 508);

    auto g_yyzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 509);

    auto g_yyzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 510);

    auto g_yyzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 511);

    auto g_yyzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 512);

    auto g_yyzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 513);

    auto g_yyzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 514);

    auto g_yyzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 515);

    auto g_yyzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 516);

    auto g_yyzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 517);

    auto g_yyzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 518);

    auto g_yyzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 519);

    auto g_yyzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 520);

    auto g_yyzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 521);

    auto g_yyzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 522);

    auto g_yyzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 523);

    auto g_yyzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 524);

    auto g_yyzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 525);

    auto g_yyzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 526);

    auto g_yyzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 527);

    auto g_yyzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 528);

    auto g_yyzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 529);

    auto g_yyzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 530);

    auto g_yyzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 531);

    auto g_yzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 532);

    auto g_yzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 534);

    auto g_yzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 536);

    auto g_yzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 537);

    auto g_yzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 539);

    auto g_yzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 540);

    auto g_yzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 541);

    auto g_yzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 543);

    auto g_yzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 544);

    auto g_yzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 545);

    auto g_yzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 546);

    auto g_yzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 548);

    auto g_yzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 549);

    auto g_yzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 550);

    auto g_yzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 551);

    auto g_yzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 552);

    auto g_yzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 554);

    auto g_yzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 555);

    auto g_yzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 556);

    auto g_yzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 557);

    auto g_yzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 558);

    auto g_yzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 559);

    auto g_zzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 560);

    auto g_zzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 561);

    auto g_zzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 562);

    auto g_zzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 563);

    auto g_zzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 564);

    auto g_zzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 565);

    auto g_zzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 566);

    auto g_zzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 567);

    auto g_zzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 568);

    auto g_zzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 569);

    auto g_zzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 570);

    auto g_zzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 571);

    auto g_zzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 572);

    auto g_zzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 573);

    auto g_zzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 574);

    auto g_zzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 575);

    auto g_zzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 576);

    auto g_zzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 577);

    auto g_zzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 578);

    auto g_zzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 579);

    auto g_zzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 580);

    auto g_zzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 581);

    auto g_zzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 582);

    auto g_zzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 583);

    auto g_zzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 584);

    auto g_zzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 585);

    auto g_zzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 586);

    auto g_zzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 587);

    // Set up components of auxiliary buffer : IH

    auto g_xxxxxx_xxxxx_1 = pbuffer.data(idx_eri_1_ih);

    auto g_xxxxxx_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 1);

    auto g_xxxxxx_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 2);

    auto g_xxxxxx_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 3);

    auto g_xxxxxx_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 4);

    auto g_xxxxxx_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 5);

    auto g_xxxxxx_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 6);

    auto g_xxxxxx_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 7);

    auto g_xxxxxx_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 8);

    auto g_xxxxxx_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 9);

    auto g_xxxxxx_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 10);

    auto g_xxxxxx_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 11);

    auto g_xxxxxx_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 12);

    auto g_xxxxxx_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 13);

    auto g_xxxxxx_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 14);

    auto g_xxxxxx_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 15);

    auto g_xxxxxx_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 16);

    auto g_xxxxxx_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 17);

    auto g_xxxxxx_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 18);

    auto g_xxxxxx_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 19);

    auto g_xxxxxx_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 20);

    auto g_xxxxxz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 44);

    auto g_xxxxxz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 46);

    auto g_xxxxxz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 47);

    auto g_xxxxxz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 49);

    auto g_xxxxxz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 50);

    auto g_xxxxxz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 51);

    auto g_xxxxxz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 53);

    auto g_xxxxxz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 54);

    auto g_xxxxxz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 55);

    auto g_xxxxxz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 56);

    auto g_xxxxxz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 58);

    auto g_xxxxxz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 59);

    auto g_xxxxxz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 60);

    auto g_xxxxxz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 61);

    auto g_xxxxxz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 62);

    auto g_xxxxyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 63);

    auto g_xxxxyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 64);

    auto g_xxxxyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 65);

    auto g_xxxxyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 66);

    auto g_xxxxyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 67);

    auto g_xxxxyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 68);

    auto g_xxxxyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 69);

    auto g_xxxxyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 70);

    auto g_xxxxyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 71);

    auto g_xxxxyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 72);

    auto g_xxxxyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 73);

    auto g_xxxxyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 74);

    auto g_xxxxyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 75);

    auto g_xxxxyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 76);

    auto g_xxxxyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 77);

    auto g_xxxxyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 78);

    auto g_xxxxyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 79);

    auto g_xxxxyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 80);

    auto g_xxxxyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 81);

    auto g_xxxxyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 82);

    auto g_xxxxyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 83);

    auto g_xxxxzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 105);

    auto g_xxxxzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 106);

    auto g_xxxxzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 107);

    auto g_xxxxzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 108);

    auto g_xxxxzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 109);

    auto g_xxxxzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 110);

    auto g_xxxxzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 111);

    auto g_xxxxzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 112);

    auto g_xxxxzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 113);

    auto g_xxxxzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 114);

    auto g_xxxxzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 115);

    auto g_xxxxzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 116);

    auto g_xxxxzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 117);

    auto g_xxxxzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 118);

    auto g_xxxxzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 119);

    auto g_xxxxzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 120);

    auto g_xxxxzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 121);

    auto g_xxxxzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 122);

    auto g_xxxxzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 123);

    auto g_xxxxzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 124);

    auto g_xxxxzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 125);

    auto g_xxxyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 126);

    auto g_xxxyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 127);

    auto g_xxxyyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 128);

    auto g_xxxyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 129);

    auto g_xxxyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 130);

    auto g_xxxyyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 131);

    auto g_xxxyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 132);

    auto g_xxxyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 133);

    auto g_xxxyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 134);

    auto g_xxxyyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 135);

    auto g_xxxyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 136);

    auto g_xxxyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 137);

    auto g_xxxyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 138);

    auto g_xxxyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 139);

    auto g_xxxyyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 140);

    auto g_xxxyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 141);

    auto g_xxxyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 142);

    auto g_xxxyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 143);

    auto g_xxxyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 144);

    auto g_xxxyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 145);

    auto g_xxxyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 146);

    auto g_xxxzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 189);

    auto g_xxxzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 190);

    auto g_xxxzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 191);

    auto g_xxxzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 192);

    auto g_xxxzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 193);

    auto g_xxxzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 194);

    auto g_xxxzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 195);

    auto g_xxxzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 196);

    auto g_xxxzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 197);

    auto g_xxxzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 198);

    auto g_xxxzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 199);

    auto g_xxxzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 200);

    auto g_xxxzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 201);

    auto g_xxxzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 202);

    auto g_xxxzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 203);

    auto g_xxxzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 204);

    auto g_xxxzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 205);

    auto g_xxxzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 206);

    auto g_xxxzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 207);

    auto g_xxxzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 208);

    auto g_xxxzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 209);

    auto g_xxyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 210);

    auto g_xxyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 211);

    auto g_xxyyyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 212);

    auto g_xxyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 213);

    auto g_xxyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 214);

    auto g_xxyyyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 215);

    auto g_xxyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 216);

    auto g_xxyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 217);

    auto g_xxyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 218);

    auto g_xxyyyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 219);

    auto g_xxyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 220);

    auto g_xxyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 221);

    auto g_xxyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 222);

    auto g_xxyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 223);

    auto g_xxyyyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 224);

    auto g_xxyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 225);

    auto g_xxyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 226);

    auto g_xxyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 227);

    auto g_xxyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 228);

    auto g_xxyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 229);

    auto g_xxyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 230);

    auto g_xxyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 256);

    auto g_xxyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 259);

    auto g_xxyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 260);

    auto g_xxyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 263);

    auto g_xxyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 264);

    auto g_xxyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 265);

    auto g_xxyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 268);

    auto g_xxyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 269);

    auto g_xxyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 270);

    auto g_xxyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 271);

    auto g_xxzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 294);

    auto g_xxzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 295);

    auto g_xxzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 296);

    auto g_xxzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 297);

    auto g_xxzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 298);

    auto g_xxzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 299);

    auto g_xxzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 300);

    auto g_xxzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 301);

    auto g_xxzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 302);

    auto g_xxzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 303);

    auto g_xxzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 304);

    auto g_xxzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 305);

    auto g_xxzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 306);

    auto g_xxzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 307);

    auto g_xxzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 308);

    auto g_xxzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 309);

    auto g_xxzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 310);

    auto g_xxzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 311);

    auto g_xxzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 312);

    auto g_xxzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 313);

    auto g_xxzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 314);

    auto g_xyyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 316);

    auto g_xyyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 318);

    auto g_xyyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 319);

    auto g_xyyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 321);

    auto g_xyyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 322);

    auto g_xyyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 323);

    auto g_xyyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 325);

    auto g_xyyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 326);

    auto g_xyyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 327);

    auto g_xyyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 328);

    auto g_xyyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 330);

    auto g_xyyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 331);

    auto g_xyyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 332);

    auto g_xyyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 333);

    auto g_xyyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 334);

    auto g_xyyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 361);

    auto g_xyyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 364);

    auto g_xyyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 365);

    auto g_xyyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 368);

    auto g_xyyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 369);

    auto g_xyyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 370);

    auto g_xyyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 373);

    auto g_xyyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 374);

    auto g_xyyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 375);

    auto g_xyyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 376);

    auto g_xyyzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 382);

    auto g_xyyzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 385);

    auto g_xyyzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 386);

    auto g_xyyzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 389);

    auto g_xyyzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 390);

    auto g_xyyzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 391);

    auto g_xyyzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 394);

    auto g_xyyzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 395);

    auto g_xyyzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 396);

    auto g_xyyzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 397);

    auto g_xzzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 422);

    auto g_xzzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 424);

    auto g_xzzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 425);

    auto g_xzzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 427);

    auto g_xzzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 428);

    auto g_xzzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 429);

    auto g_xzzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 431);

    auto g_xzzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 432);

    auto g_xzzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 433);

    auto g_xzzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 434);

    auto g_xzzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 436);

    auto g_xzzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 437);

    auto g_xzzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 438);

    auto g_xzzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 439);

    auto g_xzzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 440);

    auto g_yyyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 441);

    auto g_yyyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 442);

    auto g_yyyyyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 443);

    auto g_yyyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 444);

    auto g_yyyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 445);

    auto g_yyyyyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 446);

    auto g_yyyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 447);

    auto g_yyyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 448);

    auto g_yyyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 449);

    auto g_yyyyyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 450);

    auto g_yyyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 451);

    auto g_yyyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 452);

    auto g_yyyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 453);

    auto g_yyyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 454);

    auto g_yyyyyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 455);

    auto g_yyyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 456);

    auto g_yyyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 457);

    auto g_yyyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 458);

    auto g_yyyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 459);

    auto g_yyyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 460);

    auto g_yyyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 461);

    auto g_yyyyyz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 464);

    auto g_yyyyyz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 466);

    auto g_yyyyyz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 467);

    auto g_yyyyyz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 469);

    auto g_yyyyyz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 470);

    auto g_yyyyyz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 471);

    auto g_yyyyyz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 473);

    auto g_yyyyyz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 474);

    auto g_yyyyyz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 475);

    auto g_yyyyyz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 476);

    auto g_yyyyyz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 478);

    auto g_yyyyyz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 479);

    auto g_yyyyyz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 480);

    auto g_yyyyyz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 481);

    auto g_yyyyyz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 482);

    auto g_yyyyzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 483);

    auto g_yyyyzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 484);

    auto g_yyyyzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 485);

    auto g_yyyyzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 486);

    auto g_yyyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 487);

    auto g_yyyyzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 488);

    auto g_yyyyzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 489);

    auto g_yyyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 490);

    auto g_yyyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 491);

    auto g_yyyyzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 492);

    auto g_yyyyzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 493);

    auto g_yyyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 494);

    auto g_yyyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 495);

    auto g_yyyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 496);

    auto g_yyyyzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 497);

    auto g_yyyyzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 498);

    auto g_yyyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 499);

    auto g_yyyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 500);

    auto g_yyyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 501);

    auto g_yyyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 502);

    auto g_yyyyzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 503);

    auto g_yyyzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 504);

    auto g_yyyzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 505);

    auto g_yyyzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 506);

    auto g_yyyzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 507);

    auto g_yyyzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 508);

    auto g_yyyzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 509);

    auto g_yyyzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 510);

    auto g_yyyzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 511);

    auto g_yyyzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 512);

    auto g_yyyzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 513);

    auto g_yyyzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 514);

    auto g_yyyzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 515);

    auto g_yyyzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 516);

    auto g_yyyzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 517);

    auto g_yyyzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 518);

    auto g_yyyzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 519);

    auto g_yyyzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 520);

    auto g_yyyzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 521);

    auto g_yyyzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 522);

    auto g_yyyzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 523);

    auto g_yyyzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 524);

    auto g_yyzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 525);

    auto g_yyzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 526);

    auto g_yyzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 527);

    auto g_yyzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 528);

    auto g_yyzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 529);

    auto g_yyzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 530);

    auto g_yyzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 531);

    auto g_yyzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 532);

    auto g_yyzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 533);

    auto g_yyzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 534);

    auto g_yyzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 535);

    auto g_yyzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 536);

    auto g_yyzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 537);

    auto g_yyzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 538);

    auto g_yyzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 539);

    auto g_yyzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 540);

    auto g_yyzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 541);

    auto g_yyzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 542);

    auto g_yyzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 543);

    auto g_yyzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 544);

    auto g_yyzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 545);

    auto g_yzzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 547);

    auto g_yzzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 548);

    auto g_yzzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 549);

    auto g_yzzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 550);

    auto g_yzzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 551);

    auto g_yzzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 552);

    auto g_yzzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 553);

    auto g_yzzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 554);

    auto g_yzzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 555);

    auto g_yzzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 556);

    auto g_yzzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 557);

    auto g_yzzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 558);

    auto g_yzzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 559);

    auto g_yzzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 560);

    auto g_yzzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 561);

    auto g_yzzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 562);

    auto g_yzzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 563);

    auto g_yzzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 564);

    auto g_yzzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 565);

    auto g_yzzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 566);

    auto g_zzzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 567);

    auto g_zzzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 568);

    auto g_zzzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 569);

    auto g_zzzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 570);

    auto g_zzzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 571);

    auto g_zzzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 572);

    auto g_zzzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 573);

    auto g_zzzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 574);

    auto g_zzzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 575);

    auto g_zzzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 576);

    auto g_zzzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 577);

    auto g_zzzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 578);

    auto g_zzzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 579);

    auto g_zzzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 580);

    auto g_zzzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 581);

    auto g_zzzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 582);

    auto g_zzzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 583);

    auto g_zzzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 584);

    auto g_zzzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 585);

    auto g_zzzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 586);

    auto g_zzzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 587);

    // Set up components of auxiliary buffer : II

    auto g_xxxxxx_xxxxxx_1 = pbuffer.data(idx_eri_1_ii);

    auto g_xxxxxx_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 1);

    auto g_xxxxxx_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 2);

    auto g_xxxxxx_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 3);

    auto g_xxxxxx_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 4);

    auto g_xxxxxx_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 5);

    auto g_xxxxxx_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 6);

    auto g_xxxxxx_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 7);

    auto g_xxxxxx_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 8);

    auto g_xxxxxx_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 9);

    auto g_xxxxxx_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 10);

    auto g_xxxxxx_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 11);

    auto g_xxxxxx_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 12);

    auto g_xxxxxx_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 13);

    auto g_xxxxxx_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 14);

    auto g_xxxxxx_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 15);

    auto g_xxxxxx_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 16);

    auto g_xxxxxx_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 17);

    auto g_xxxxxx_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 18);

    auto g_xxxxxx_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 19);

    auto g_xxxxxx_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 20);

    auto g_xxxxxx_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 21);

    auto g_xxxxxx_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 22);

    auto g_xxxxxx_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 23);

    auto g_xxxxxx_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 24);

    auto g_xxxxxx_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 25);

    auto g_xxxxxx_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 26);

    auto g_xxxxxx_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 27);

    auto g_xxxxxy_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 28);

    auto g_xxxxxy_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 29);

    auto g_xxxxxy_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 30);

    auto g_xxxxxy_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 31);

    auto g_xxxxxy_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 33);

    auto g_xxxxxy_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 34);

    auto g_xxxxxy_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 37);

    auto g_xxxxxy_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 38);

    auto g_xxxxxy_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 42);

    auto g_xxxxxy_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 43);

    auto g_xxxxxy_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 48);

    auto g_xxxxxy_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 49);

    auto g_xxxxxz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 56);

    auto g_xxxxxz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 57);

    auto g_xxxxxz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 58);

    auto g_xxxxxz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 59);

    auto g_xxxxxz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 60);

    auto g_xxxxxz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 61);

    auto g_xxxxxz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 62);

    auto g_xxxxxz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 63);

    auto g_xxxxxz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 64);

    auto g_xxxxxz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 65);

    auto g_xxxxxz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 66);

    auto g_xxxxxz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 67);

    auto g_xxxxxz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 68);

    auto g_xxxxxz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 69);

    auto g_xxxxxz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 70);

    auto g_xxxxxz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 71);

    auto g_xxxxxz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 72);

    auto g_xxxxxz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 73);

    auto g_xxxxxz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 74);

    auto g_xxxxxz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 75);

    auto g_xxxxxz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 76);

    auto g_xxxxxz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 78);

    auto g_xxxxxz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 79);

    auto g_xxxxxz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 80);

    auto g_xxxxxz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 81);

    auto g_xxxxxz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 82);

    auto g_xxxxxz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 83);

    auto g_xxxxyy_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 84);

    auto g_xxxxyy_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 85);

    auto g_xxxxyy_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 86);

    auto g_xxxxyy_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 87);

    auto g_xxxxyy_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 88);

    auto g_xxxxyy_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 89);

    auto g_xxxxyy_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 90);

    auto g_xxxxyy_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 91);

    auto g_xxxxyy_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 92);

    auto g_xxxxyy_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 93);

    auto g_xxxxyy_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 94);

    auto g_xxxxyy_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 95);

    auto g_xxxxyy_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 96);

    auto g_xxxxyy_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 97);

    auto g_xxxxyy_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 98);

    auto g_xxxxyy_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 99);

    auto g_xxxxyy_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 100);

    auto g_xxxxyy_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 101);

    auto g_xxxxyy_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 102);

    auto g_xxxxyy_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 103);

    auto g_xxxxyy_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 104);

    auto g_xxxxyy_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 105);

    auto g_xxxxyy_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 106);

    auto g_xxxxyy_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 107);

    auto g_xxxxyy_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 108);

    auto g_xxxxyy_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 109);

    auto g_xxxxyy_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 110);

    auto g_xxxxyy_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 111);

    auto g_xxxxzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 140);

    auto g_xxxxzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 141);

    auto g_xxxxzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 142);

    auto g_xxxxzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 143);

    auto g_xxxxzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 144);

    auto g_xxxxzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 145);

    auto g_xxxxzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 146);

    auto g_xxxxzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 147);

    auto g_xxxxzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 148);

    auto g_xxxxzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 149);

    auto g_xxxxzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 150);

    auto g_xxxxzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 151);

    auto g_xxxxzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 152);

    auto g_xxxxzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 153);

    auto g_xxxxzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 154);

    auto g_xxxxzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 155);

    auto g_xxxxzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 156);

    auto g_xxxxzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 157);

    auto g_xxxxzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 158);

    auto g_xxxxzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 159);

    auto g_xxxxzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 160);

    auto g_xxxxzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 161);

    auto g_xxxxzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 162);

    auto g_xxxxzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 163);

    auto g_xxxxzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 164);

    auto g_xxxxzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 165);

    auto g_xxxxzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 166);

    auto g_xxxxzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 167);

    auto g_xxxyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 168);

    auto g_xxxyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 169);

    auto g_xxxyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 170);

    auto g_xxxyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 171);

    auto g_xxxyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 172);

    auto g_xxxyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 173);

    auto g_xxxyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 174);

    auto g_xxxyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 175);

    auto g_xxxyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 176);

    auto g_xxxyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 177);

    auto g_xxxyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 178);

    auto g_xxxyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 179);

    auto g_xxxyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 180);

    auto g_xxxyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 181);

    auto g_xxxyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 182);

    auto g_xxxyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 183);

    auto g_xxxyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 184);

    auto g_xxxyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 185);

    auto g_xxxyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 186);

    auto g_xxxyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 187);

    auto g_xxxyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 188);

    auto g_xxxyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 189);

    auto g_xxxyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 190);

    auto g_xxxyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 191);

    auto g_xxxyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 192);

    auto g_xxxyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 193);

    auto g_xxxyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 194);

    auto g_xxxyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 195);

    auto g_xxxyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 197);

    auto g_xxxyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 199);

    auto g_xxxyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 202);

    auto g_xxxyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 206);

    auto g_xxxyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 211);

    auto g_xxxyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 224);

    auto g_xxxyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 226);

    auto g_xxxyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 229);

    auto g_xxxyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 233);

    auto g_xxxyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 238);

    auto g_xxxyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 244);

    auto g_xxxzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 252);

    auto g_xxxzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 253);

    auto g_xxxzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 254);

    auto g_xxxzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 255);

    auto g_xxxzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 256);

    auto g_xxxzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 257);

    auto g_xxxzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 258);

    auto g_xxxzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 259);

    auto g_xxxzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 260);

    auto g_xxxzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 261);

    auto g_xxxzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 262);

    auto g_xxxzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 263);

    auto g_xxxzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 264);

    auto g_xxxzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 265);

    auto g_xxxzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 266);

    auto g_xxxzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 267);

    auto g_xxxzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 268);

    auto g_xxxzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 269);

    auto g_xxxzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 270);

    auto g_xxxzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 271);

    auto g_xxxzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 272);

    auto g_xxxzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 273);

    auto g_xxxzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 274);

    auto g_xxxzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 275);

    auto g_xxxzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 276);

    auto g_xxxzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 277);

    auto g_xxxzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 278);

    auto g_xxxzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 279);

    auto g_xxyyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 280);

    auto g_xxyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 281);

    auto g_xxyyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 282);

    auto g_xxyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 283);

    auto g_xxyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 284);

    auto g_xxyyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 285);

    auto g_xxyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 286);

    auto g_xxyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 287);

    auto g_xxyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 288);

    auto g_xxyyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 289);

    auto g_xxyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 290);

    auto g_xxyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 291);

    auto g_xxyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 292);

    auto g_xxyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 293);

    auto g_xxyyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 294);

    auto g_xxyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 295);

    auto g_xxyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 296);

    auto g_xxyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 297);

    auto g_xxyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 298);

    auto g_xxyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 299);

    auto g_xxyyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 300);

    auto g_xxyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 301);

    auto g_xxyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 302);

    auto g_xxyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 303);

    auto g_xxyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 304);

    auto g_xxyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 305);

    auto g_xxyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 306);

    auto g_xxyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 307);

    auto g_xxyyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 309);

    auto g_xxyyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 311);

    auto g_xxyyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 314);

    auto g_xxyyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 318);

    auto g_xxyyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 323);

    auto g_xxyyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 336);

    auto g_xxyyzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 337);

    auto g_xxyyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 338);

    auto g_xxyyzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 339);

    auto g_xxyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 340);

    auto g_xxyyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 341);

    auto g_xxyyzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 342);

    auto g_xxyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 343);

    auto g_xxyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 344);

    auto g_xxyyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 345);

    auto g_xxyyzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 346);

    auto g_xxyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 347);

    auto g_xxyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 348);

    auto g_xxyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 349);

    auto g_xxyyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 350);

    auto g_xxyyzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 351);

    auto g_xxyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 352);

    auto g_xxyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 353);

    auto g_xxyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 354);

    auto g_xxyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 355);

    auto g_xxyyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 356);

    auto g_xxyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 357);

    auto g_xxyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 358);

    auto g_xxyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 359);

    auto g_xxyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 360);

    auto g_xxyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 361);

    auto g_xxyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 362);

    auto g_xxyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 363);

    auto g_xxyzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 364);

    auto g_xxyzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 366);

    auto g_xxyzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 369);

    auto g_xxyzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 373);

    auto g_xxyzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 378);

    auto g_xxyzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 384);

    auto g_xxzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 392);

    auto g_xxzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 393);

    auto g_xxzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 394);

    auto g_xxzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 395);

    auto g_xxzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 396);

    auto g_xxzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 397);

    auto g_xxzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 398);

    auto g_xxzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 399);

    auto g_xxzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 400);

    auto g_xxzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 401);

    auto g_xxzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 402);

    auto g_xxzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 403);

    auto g_xxzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 404);

    auto g_xxzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 405);

    auto g_xxzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 406);

    auto g_xxzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 407);

    auto g_xxzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 408);

    auto g_xxzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 409);

    auto g_xxzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 410);

    auto g_xxzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 411);

    auto g_xxzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 412);

    auto g_xxzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 413);

    auto g_xxzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 414);

    auto g_xxzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 415);

    auto g_xxzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 416);

    auto g_xxzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 417);

    auto g_xxzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 418);

    auto g_xxzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 419);

    auto g_xyyyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 420);

    auto g_xyyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 421);

    auto g_xyyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 423);

    auto g_xyyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 424);

    auto g_xyyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 426);

    auto g_xyyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 427);

    auto g_xyyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 428);

    auto g_xyyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 430);

    auto g_xyyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 431);

    auto g_xyyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 432);

    auto g_xyyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 433);

    auto g_xyyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 435);

    auto g_xyyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 436);

    auto g_xyyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 437);

    auto g_xyyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 438);

    auto g_xyyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 439);

    auto g_xyyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 441);

    auto g_xyyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 442);

    auto g_xyyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 443);

    auto g_xyyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 444);

    auto g_xyyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 445);

    auto g_xyyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 446);

    auto g_xyyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 447);

    auto g_xyyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 480);

    auto g_xyyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 483);

    auto g_xyyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 484);

    auto g_xyyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 487);

    auto g_xyyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 488);

    auto g_xyyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 489);

    auto g_xyyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 492);

    auto g_xyyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 493);

    auto g_xyyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 494);

    auto g_xyyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 495);

    auto g_xyyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 497);

    auto g_xyyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 498);

    auto g_xyyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 499);

    auto g_xyyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 500);

    auto g_xyyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 501);

    auto g_xyyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 502);

    auto g_xyyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 503);

    auto g_xyyzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 508);

    auto g_xyyzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 511);

    auto g_xyyzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 512);

    auto g_xyyzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 515);

    auto g_xyyzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 516);

    auto g_xyyzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 517);

    auto g_xyyzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 520);

    auto g_xyyzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 521);

    auto g_xyyzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 522);

    auto g_xyyzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 523);

    auto g_xyyzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 525);

    auto g_xyyzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 526);

    auto g_xyyzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 527);

    auto g_xyyzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 528);

    auto g_xyyzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 529);

    auto g_xyyzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 530);

    auto g_xyyzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 531);

    auto g_xzzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 560);

    auto g_xzzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 562);

    auto g_xzzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 564);

    auto g_xzzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 565);

    auto g_xzzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 567);

    auto g_xzzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 568);

    auto g_xzzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 569);

    auto g_xzzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 571);

    auto g_xzzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 572);

    auto g_xzzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 573);

    auto g_xzzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 574);

    auto g_xzzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 576);

    auto g_xzzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 577);

    auto g_xzzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 578);

    auto g_xzzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 579);

    auto g_xzzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 580);

    auto g_xzzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 581);

    auto g_xzzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 582);

    auto g_xzzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 583);

    auto g_xzzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 584);

    auto g_xzzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 585);

    auto g_xzzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 586);

    auto g_xzzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 587);

    auto g_yyyyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 588);

    auto g_yyyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 589);

    auto g_yyyyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 590);

    auto g_yyyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 591);

    auto g_yyyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 592);

    auto g_yyyyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 593);

    auto g_yyyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 594);

    auto g_yyyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 595);

    auto g_yyyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 596);

    auto g_yyyyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 597);

    auto g_yyyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 598);

    auto g_yyyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 599);

    auto g_yyyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 600);

    auto g_yyyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 601);

    auto g_yyyyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 602);

    auto g_yyyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 603);

    auto g_yyyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 604);

    auto g_yyyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 605);

    auto g_yyyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 606);

    auto g_yyyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 607);

    auto g_yyyyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 608);

    auto g_yyyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 609);

    auto g_yyyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 610);

    auto g_yyyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 611);

    auto g_yyyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 612);

    auto g_yyyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 613);

    auto g_yyyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 614);

    auto g_yyyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 615);

    auto g_yyyyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 617);

    auto g_yyyyyz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 618);

    auto g_yyyyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 619);

    auto g_yyyyyz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 620);

    auto g_yyyyyz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 621);

    auto g_yyyyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 622);

    auto g_yyyyyz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 623);

    auto g_yyyyyz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 624);

    auto g_yyyyyz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 625);

    auto g_yyyyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 626);

    auto g_yyyyyz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 627);

    auto g_yyyyyz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 628);

    auto g_yyyyyz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 629);

    auto g_yyyyyz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 630);

    auto g_yyyyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 631);

    auto g_yyyyyz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 632);

    auto g_yyyyyz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 633);

    auto g_yyyyyz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 634);

    auto g_yyyyyz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 635);

    auto g_yyyyyz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 636);

    auto g_yyyyyz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 637);

    auto g_yyyyyz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 638);

    auto g_yyyyyz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 639);

    auto g_yyyyyz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 640);

    auto g_yyyyyz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 641);

    auto g_yyyyyz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 642);

    auto g_yyyyyz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 643);

    auto g_yyyyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 644);

    auto g_yyyyzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 645);

    auto g_yyyyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 646);

    auto g_yyyyzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 647);

    auto g_yyyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 648);

    auto g_yyyyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 649);

    auto g_yyyyzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 650);

    auto g_yyyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 651);

    auto g_yyyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 652);

    auto g_yyyyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 653);

    auto g_yyyyzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 654);

    auto g_yyyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 655);

    auto g_yyyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 656);

    auto g_yyyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 657);

    auto g_yyyyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 658);

    auto g_yyyyzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 659);

    auto g_yyyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 660);

    auto g_yyyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 661);

    auto g_yyyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 662);

    auto g_yyyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 663);

    auto g_yyyyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 664);

    auto g_yyyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 665);

    auto g_yyyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 666);

    auto g_yyyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 667);

    auto g_yyyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 668);

    auto g_yyyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 669);

    auto g_yyyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 670);

    auto g_yyyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 671);

    auto g_yyyzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 672);

    auto g_yyyzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 673);

    auto g_yyyzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 674);

    auto g_yyyzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 675);

    auto g_yyyzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 676);

    auto g_yyyzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 677);

    auto g_yyyzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 678);

    auto g_yyyzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 679);

    auto g_yyyzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 680);

    auto g_yyyzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 681);

    auto g_yyyzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 682);

    auto g_yyyzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 683);

    auto g_yyyzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 684);

    auto g_yyyzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 685);

    auto g_yyyzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 686);

    auto g_yyyzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 687);

    auto g_yyyzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 688);

    auto g_yyyzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 689);

    auto g_yyyzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 690);

    auto g_yyyzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 691);

    auto g_yyyzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 692);

    auto g_yyyzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 693);

    auto g_yyyzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 694);

    auto g_yyyzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 695);

    auto g_yyyzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 696);

    auto g_yyyzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 697);

    auto g_yyyzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 698);

    auto g_yyyzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 699);

    auto g_yyzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 700);

    auto g_yyzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 701);

    auto g_yyzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 702);

    auto g_yyzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 703);

    auto g_yyzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 704);

    auto g_yyzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 705);

    auto g_yyzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 706);

    auto g_yyzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 707);

    auto g_yyzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 708);

    auto g_yyzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 709);

    auto g_yyzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 710);

    auto g_yyzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 711);

    auto g_yyzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 712);

    auto g_yyzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 713);

    auto g_yyzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 714);

    auto g_yyzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 715);

    auto g_yyzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 716);

    auto g_yyzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 717);

    auto g_yyzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 718);

    auto g_yyzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 719);

    auto g_yyzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 720);

    auto g_yyzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 721);

    auto g_yyzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 722);

    auto g_yyzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 723);

    auto g_yyzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 724);

    auto g_yyzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 725);

    auto g_yyzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 726);

    auto g_yyzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 727);

    auto g_yzzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 728);

    auto g_yzzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 729);

    auto g_yzzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 730);

    auto g_yzzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 731);

    auto g_yzzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 732);

    auto g_yzzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 733);

    auto g_yzzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 734);

    auto g_yzzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 735);

    auto g_yzzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 736);

    auto g_yzzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 737);

    auto g_yzzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 738);

    auto g_yzzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 739);

    auto g_yzzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 740);

    auto g_yzzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 741);

    auto g_yzzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 742);

    auto g_yzzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 743);

    auto g_yzzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 744);

    auto g_yzzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 745);

    auto g_yzzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 746);

    auto g_yzzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 747);

    auto g_yzzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 748);

    auto g_yzzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 749);

    auto g_yzzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 750);

    auto g_yzzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 751);

    auto g_yzzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 752);

    auto g_yzzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 753);

    auto g_yzzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 754);

    auto g_yzzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 755);

    auto g_zzzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_ii + 756);

    auto g_zzzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_ii + 757);

    auto g_zzzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_ii + 758);

    auto g_zzzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_ii + 759);

    auto g_zzzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_ii + 760);

    auto g_zzzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_ii + 761);

    auto g_zzzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_ii + 762);

    auto g_zzzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_ii + 763);

    auto g_zzzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_ii + 764);

    auto g_zzzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_ii + 765);

    auto g_zzzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_ii + 766);

    auto g_zzzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_ii + 767);

    auto g_zzzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_ii + 768);

    auto g_zzzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_ii + 769);

    auto g_zzzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_ii + 770);

    auto g_zzzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_ii + 771);

    auto g_zzzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_ii + 772);

    auto g_zzzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_ii + 773);

    auto g_zzzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_ii + 774);

    auto g_zzzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_ii + 775);

    auto g_zzzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_ii + 776);

    auto g_zzzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_ii + 777);

    auto g_zzzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_ii + 778);

    auto g_zzzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_ii + 779);

    auto g_zzzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_ii + 780);

    auto g_zzzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_ii + 781);

    auto g_zzzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_ii + 782);

    auto g_zzzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_ii + 783);

    // Set up 0-28 components of targeted buffer : KI

    auto g_xxxxxxx_xxxxxx_0 = pbuffer.data(idx_eri_0_ki);

    auto g_xxxxxxx_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 1);

    auto g_xxxxxxx_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 2);

    auto g_xxxxxxx_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 3);

    auto g_xxxxxxx_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 4);

    auto g_xxxxxxx_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 5);

    auto g_xxxxxxx_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 6);

    auto g_xxxxxxx_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 7);

    auto g_xxxxxxx_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 8);

    auto g_xxxxxxx_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 9);

    auto g_xxxxxxx_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 10);

    auto g_xxxxxxx_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 11);

    auto g_xxxxxxx_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 12);

    auto g_xxxxxxx_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 13);

    auto g_xxxxxxx_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 14);

    auto g_xxxxxxx_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 15);

    auto g_xxxxxxx_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 16);

    auto g_xxxxxxx_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 17);

    auto g_xxxxxxx_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 18);

    auto g_xxxxxxx_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 19);

    auto g_xxxxxxx_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 20);

    auto g_xxxxxxx_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 21);

    auto g_xxxxxxx_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 22);

    auto g_xxxxxxx_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 23);

    auto g_xxxxxxx_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 24);

    auto g_xxxxxxx_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 25);

    auto g_xxxxxxx_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 26);

    auto g_xxxxxxx_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 27);

    #pragma omp simd aligned(g_xxxxx_xxxxxx_0, g_xxxxx_xxxxxx_1, g_xxxxx_xxxxxy_0, g_xxxxx_xxxxxy_1, g_xxxxx_xxxxxz_0, g_xxxxx_xxxxxz_1, g_xxxxx_xxxxyy_0, g_xxxxx_xxxxyy_1, g_xxxxx_xxxxyz_0, g_xxxxx_xxxxyz_1, g_xxxxx_xxxxzz_0, g_xxxxx_xxxxzz_1, g_xxxxx_xxxyyy_0, g_xxxxx_xxxyyy_1, g_xxxxx_xxxyyz_0, g_xxxxx_xxxyyz_1, g_xxxxx_xxxyzz_0, g_xxxxx_xxxyzz_1, g_xxxxx_xxxzzz_0, g_xxxxx_xxxzzz_1, g_xxxxx_xxyyyy_0, g_xxxxx_xxyyyy_1, g_xxxxx_xxyyyz_0, g_xxxxx_xxyyyz_1, g_xxxxx_xxyyzz_0, g_xxxxx_xxyyzz_1, g_xxxxx_xxyzzz_0, g_xxxxx_xxyzzz_1, g_xxxxx_xxzzzz_0, g_xxxxx_xxzzzz_1, g_xxxxx_xyyyyy_0, g_xxxxx_xyyyyy_1, g_xxxxx_xyyyyz_0, g_xxxxx_xyyyyz_1, g_xxxxx_xyyyzz_0, g_xxxxx_xyyyzz_1, g_xxxxx_xyyzzz_0, g_xxxxx_xyyzzz_1, g_xxxxx_xyzzzz_0, g_xxxxx_xyzzzz_1, g_xxxxx_xzzzzz_0, g_xxxxx_xzzzzz_1, g_xxxxx_yyyyyy_0, g_xxxxx_yyyyyy_1, g_xxxxx_yyyyyz_0, g_xxxxx_yyyyyz_1, g_xxxxx_yyyyzz_0, g_xxxxx_yyyyzz_1, g_xxxxx_yyyzzz_0, g_xxxxx_yyyzzz_1, g_xxxxx_yyzzzz_0, g_xxxxx_yyzzzz_1, g_xxxxx_yzzzzz_0, g_xxxxx_yzzzzz_1, g_xxxxx_zzzzzz_0, g_xxxxx_zzzzzz_1, g_xxxxxx_xxxxx_1, g_xxxxxx_xxxxxx_1, g_xxxxxx_xxxxxy_1, g_xxxxxx_xxxxxz_1, g_xxxxxx_xxxxy_1, g_xxxxxx_xxxxyy_1, g_xxxxxx_xxxxyz_1, g_xxxxxx_xxxxz_1, g_xxxxxx_xxxxzz_1, g_xxxxxx_xxxyy_1, g_xxxxxx_xxxyyy_1, g_xxxxxx_xxxyyz_1, g_xxxxxx_xxxyz_1, g_xxxxxx_xxxyzz_1, g_xxxxxx_xxxzz_1, g_xxxxxx_xxxzzz_1, g_xxxxxx_xxyyy_1, g_xxxxxx_xxyyyy_1, g_xxxxxx_xxyyyz_1, g_xxxxxx_xxyyz_1, g_xxxxxx_xxyyzz_1, g_xxxxxx_xxyzz_1, g_xxxxxx_xxyzzz_1, g_xxxxxx_xxzzz_1, g_xxxxxx_xxzzzz_1, g_xxxxxx_xyyyy_1, g_xxxxxx_xyyyyy_1, g_xxxxxx_xyyyyz_1, g_xxxxxx_xyyyz_1, g_xxxxxx_xyyyzz_1, g_xxxxxx_xyyzz_1, g_xxxxxx_xyyzzz_1, g_xxxxxx_xyzzz_1, g_xxxxxx_xyzzzz_1, g_xxxxxx_xzzzz_1, g_xxxxxx_xzzzzz_1, g_xxxxxx_yyyyy_1, g_xxxxxx_yyyyyy_1, g_xxxxxx_yyyyyz_1, g_xxxxxx_yyyyz_1, g_xxxxxx_yyyyzz_1, g_xxxxxx_yyyzz_1, g_xxxxxx_yyyzzz_1, g_xxxxxx_yyzzz_1, g_xxxxxx_yyzzzz_1, g_xxxxxx_yzzzz_1, g_xxxxxx_yzzzzz_1, g_xxxxxx_zzzzz_1, g_xxxxxx_zzzzzz_1, g_xxxxxxx_xxxxxx_0, g_xxxxxxx_xxxxxy_0, g_xxxxxxx_xxxxxz_0, g_xxxxxxx_xxxxyy_0, g_xxxxxxx_xxxxyz_0, g_xxxxxxx_xxxxzz_0, g_xxxxxxx_xxxyyy_0, g_xxxxxxx_xxxyyz_0, g_xxxxxxx_xxxyzz_0, g_xxxxxxx_xxxzzz_0, g_xxxxxxx_xxyyyy_0, g_xxxxxxx_xxyyyz_0, g_xxxxxxx_xxyyzz_0, g_xxxxxxx_xxyzzz_0, g_xxxxxxx_xxzzzz_0, g_xxxxxxx_xyyyyy_0, g_xxxxxxx_xyyyyz_0, g_xxxxxxx_xyyyzz_0, g_xxxxxxx_xyyzzz_0, g_xxxxxxx_xyzzzz_0, g_xxxxxxx_xzzzzz_0, g_xxxxxxx_yyyyyy_0, g_xxxxxxx_yyyyyz_0, g_xxxxxxx_yyyyzz_0, g_xxxxxxx_yyyzzz_0, g_xxxxxxx_yyzzzz_0, g_xxxxxxx_yzzzzz_0, g_xxxxxxx_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxxx_xxxxxx_0[i] = 6.0 * g_xxxxx_xxxxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxxxxx_xxxxx_1[i] * fe_0 + g_xxxxxx_xxxxxx_1[i] * pa_x[i];

        g_xxxxxxx_xxxxxy_0[i] = 6.0 * g_xxxxx_xxxxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxxxx_xxxxy_1[i] * fe_0 + g_xxxxxx_xxxxxy_1[i] * pa_x[i];

        g_xxxxxxx_xxxxxz_0[i] = 6.0 * g_xxxxx_xxxxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_xxxxz_1[i] * fe_0 + g_xxxxxx_xxxxxz_1[i] * pa_x[i];

        g_xxxxxxx_xxxxyy_0[i] = 6.0 * g_xxxxx_xxxxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxxxx_xxxyy_1[i] * fe_0 + g_xxxxxx_xxxxyy_1[i] * pa_x[i];

        g_xxxxxxx_xxxxyz_0[i] = 6.0 * g_xxxxx_xxxxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_xxxyz_1[i] * fe_0 + g_xxxxxx_xxxxyz_1[i] * pa_x[i];

        g_xxxxxxx_xxxxzz_0[i] = 6.0 * g_xxxxx_xxxxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_xxxzz_1[i] * fe_0 + g_xxxxxx_xxxxzz_1[i] * pa_x[i];

        g_xxxxxxx_xxxyyy_0[i] = 6.0 * g_xxxxx_xxxyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxyyy_1[i] * fe_0 + g_xxxxxx_xxxyyy_1[i] * pa_x[i];

        g_xxxxxxx_xxxyyz_0[i] = 6.0 * g_xxxxx_xxxyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxyyz_1[i] * fe_0 + g_xxxxxx_xxxyyz_1[i] * pa_x[i];

        g_xxxxxxx_xxxyzz_0[i] = 6.0 * g_xxxxx_xxxyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxyzz_1[i] * fe_0 + g_xxxxxx_xxxyzz_1[i] * pa_x[i];

        g_xxxxxxx_xxxzzz_0[i] = 6.0 * g_xxxxx_xxxzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxzzz_1[i] * fe_0 + g_xxxxxx_xxxzzz_1[i] * pa_x[i];

        g_xxxxxxx_xxyyyy_0[i] = 6.0 * g_xxxxx_xxyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyyyy_1[i] * fe_0 + g_xxxxxx_xxyyyy_1[i] * pa_x[i];

        g_xxxxxxx_xxyyyz_0[i] = 6.0 * g_xxxxx_xxyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyyyz_1[i] * fe_0 + g_xxxxxx_xxyyyz_1[i] * pa_x[i];

        g_xxxxxxx_xxyyzz_0[i] = 6.0 * g_xxxxx_xxyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyyzz_1[i] * fe_0 + g_xxxxxx_xxyyzz_1[i] * pa_x[i];

        g_xxxxxxx_xxyzzz_0[i] = 6.0 * g_xxxxx_xxyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyzzz_1[i] * fe_0 + g_xxxxxx_xxyzzz_1[i] * pa_x[i];

        g_xxxxxxx_xxzzzz_0[i] = 6.0 * g_xxxxx_xxzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xzzzz_1[i] * fe_0 + g_xxxxxx_xxzzzz_1[i] * pa_x[i];

        g_xxxxxxx_xyyyyy_0[i] = 6.0 * g_xxxxx_xyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyyyy_1[i] * fz_be_0 + g_xxxxxx_yyyyy_1[i] * fe_0 + g_xxxxxx_xyyyyy_1[i] * pa_x[i];

        g_xxxxxxx_xyyyyz_0[i] = 6.0 * g_xxxxx_xyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyyyz_1[i] * fz_be_0 + g_xxxxxx_yyyyz_1[i] * fe_0 + g_xxxxxx_xyyyyz_1[i] * pa_x[i];

        g_xxxxxxx_xyyyzz_0[i] = 6.0 * g_xxxxx_xyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyyzz_1[i] * fz_be_0 + g_xxxxxx_yyyzz_1[i] * fe_0 + g_xxxxxx_xyyyzz_1[i] * pa_x[i];

        g_xxxxxxx_xyyzzz_0[i] = 6.0 * g_xxxxx_xyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyzzz_1[i] * fz_be_0 + g_xxxxxx_yyzzz_1[i] * fe_0 + g_xxxxxx_xyyzzz_1[i] * pa_x[i];

        g_xxxxxxx_xyzzzz_0[i] = 6.0 * g_xxxxx_xyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyzzzz_1[i] * fz_be_0 + g_xxxxxx_yzzzz_1[i] * fe_0 + g_xxxxxx_xyzzzz_1[i] * pa_x[i];

        g_xxxxxxx_xzzzzz_0[i] = 6.0 * g_xxxxx_xzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xzzzzz_1[i] * fz_be_0 + g_xxxxxx_zzzzz_1[i] * fe_0 + g_xxxxxx_xzzzzz_1[i] * pa_x[i];

        g_xxxxxxx_yyyyyy_0[i] = 6.0 * g_xxxxx_yyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyyyy_1[i] * fz_be_0 + g_xxxxxx_yyyyyy_1[i] * pa_x[i];

        g_xxxxxxx_yyyyyz_0[i] = 6.0 * g_xxxxx_yyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyyyz_1[i] * fz_be_0 + g_xxxxxx_yyyyyz_1[i] * pa_x[i];

        g_xxxxxxx_yyyyzz_0[i] = 6.0 * g_xxxxx_yyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyyzz_1[i] * fz_be_0 + g_xxxxxx_yyyyzz_1[i] * pa_x[i];

        g_xxxxxxx_yyyzzz_0[i] = 6.0 * g_xxxxx_yyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyzzz_1[i] * fz_be_0 + g_xxxxxx_yyyzzz_1[i] * pa_x[i];

        g_xxxxxxx_yyzzzz_0[i] = 6.0 * g_xxxxx_yyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyzzzz_1[i] * fz_be_0 + g_xxxxxx_yyzzzz_1[i] * pa_x[i];

        g_xxxxxxx_yzzzzz_0[i] = 6.0 * g_xxxxx_yzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yzzzzz_1[i] * fz_be_0 + g_xxxxxx_yzzzzz_1[i] * pa_x[i];

        g_xxxxxxx_zzzzzz_0[i] = 6.0 * g_xxxxx_zzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_zzzzzz_1[i] * fz_be_0 + g_xxxxxx_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : KI

    auto g_xxxxxxy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 28);

    auto g_xxxxxxy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 29);

    auto g_xxxxxxy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 30);

    auto g_xxxxxxy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 31);

    auto g_xxxxxxy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 32);

    auto g_xxxxxxy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 33);

    auto g_xxxxxxy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 34);

    auto g_xxxxxxy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 35);

    auto g_xxxxxxy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 36);

    auto g_xxxxxxy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 37);

    auto g_xxxxxxy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 38);

    auto g_xxxxxxy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 39);

    auto g_xxxxxxy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 40);

    auto g_xxxxxxy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 41);

    auto g_xxxxxxy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 42);

    auto g_xxxxxxy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 43);

    auto g_xxxxxxy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 44);

    auto g_xxxxxxy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 45);

    auto g_xxxxxxy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 46);

    auto g_xxxxxxy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 47);

    auto g_xxxxxxy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 48);

    auto g_xxxxxxy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 49);

    auto g_xxxxxxy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 50);

    auto g_xxxxxxy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 51);

    auto g_xxxxxxy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 52);

    auto g_xxxxxxy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 53);

    auto g_xxxxxxy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 54);

    auto g_xxxxxxy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 55);

    #pragma omp simd aligned(g_xxxxxx_xxxxx_1, g_xxxxxx_xxxxxx_1, g_xxxxxx_xxxxxy_1, g_xxxxxx_xxxxxz_1, g_xxxxxx_xxxxy_1, g_xxxxxx_xxxxyy_1, g_xxxxxx_xxxxyz_1, g_xxxxxx_xxxxz_1, g_xxxxxx_xxxxzz_1, g_xxxxxx_xxxyy_1, g_xxxxxx_xxxyyy_1, g_xxxxxx_xxxyyz_1, g_xxxxxx_xxxyz_1, g_xxxxxx_xxxyzz_1, g_xxxxxx_xxxzz_1, g_xxxxxx_xxxzzz_1, g_xxxxxx_xxyyy_1, g_xxxxxx_xxyyyy_1, g_xxxxxx_xxyyyz_1, g_xxxxxx_xxyyz_1, g_xxxxxx_xxyyzz_1, g_xxxxxx_xxyzz_1, g_xxxxxx_xxyzzz_1, g_xxxxxx_xxzzz_1, g_xxxxxx_xxzzzz_1, g_xxxxxx_xyyyy_1, g_xxxxxx_xyyyyy_1, g_xxxxxx_xyyyyz_1, g_xxxxxx_xyyyz_1, g_xxxxxx_xyyyzz_1, g_xxxxxx_xyyzz_1, g_xxxxxx_xyyzzz_1, g_xxxxxx_xyzzz_1, g_xxxxxx_xyzzzz_1, g_xxxxxx_xzzzz_1, g_xxxxxx_xzzzzz_1, g_xxxxxx_yyyyy_1, g_xxxxxx_yyyyyy_1, g_xxxxxx_yyyyyz_1, g_xxxxxx_yyyyz_1, g_xxxxxx_yyyyzz_1, g_xxxxxx_yyyzz_1, g_xxxxxx_yyyzzz_1, g_xxxxxx_yyzzz_1, g_xxxxxx_yyzzzz_1, g_xxxxxx_yzzzz_1, g_xxxxxx_yzzzzz_1, g_xxxxxx_zzzzz_1, g_xxxxxx_zzzzzz_1, g_xxxxxxy_xxxxxx_0, g_xxxxxxy_xxxxxy_0, g_xxxxxxy_xxxxxz_0, g_xxxxxxy_xxxxyy_0, g_xxxxxxy_xxxxyz_0, g_xxxxxxy_xxxxzz_0, g_xxxxxxy_xxxyyy_0, g_xxxxxxy_xxxyyz_0, g_xxxxxxy_xxxyzz_0, g_xxxxxxy_xxxzzz_0, g_xxxxxxy_xxyyyy_0, g_xxxxxxy_xxyyyz_0, g_xxxxxxy_xxyyzz_0, g_xxxxxxy_xxyzzz_0, g_xxxxxxy_xxzzzz_0, g_xxxxxxy_xyyyyy_0, g_xxxxxxy_xyyyyz_0, g_xxxxxxy_xyyyzz_0, g_xxxxxxy_xyyzzz_0, g_xxxxxxy_xyzzzz_0, g_xxxxxxy_xzzzzz_0, g_xxxxxxy_yyyyyy_0, g_xxxxxxy_yyyyyz_0, g_xxxxxxy_yyyyzz_0, g_xxxxxxy_yyyzzz_0, g_xxxxxxy_yyzzzz_0, g_xxxxxxy_yzzzzz_0, g_xxxxxxy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxy_xxxxxx_0[i] = g_xxxxxx_xxxxxx_1[i] * pa_y[i];

        g_xxxxxxy_xxxxxy_0[i] = g_xxxxxx_xxxxx_1[i] * fe_0 + g_xxxxxx_xxxxxy_1[i] * pa_y[i];

        g_xxxxxxy_xxxxxz_0[i] = g_xxxxxx_xxxxxz_1[i] * pa_y[i];

        g_xxxxxxy_xxxxyy_0[i] = 2.0 * g_xxxxxx_xxxxy_1[i] * fe_0 + g_xxxxxx_xxxxyy_1[i] * pa_y[i];

        g_xxxxxxy_xxxxyz_0[i] = g_xxxxxx_xxxxz_1[i] * fe_0 + g_xxxxxx_xxxxyz_1[i] * pa_y[i];

        g_xxxxxxy_xxxxzz_0[i] = g_xxxxxx_xxxxzz_1[i] * pa_y[i];

        g_xxxxxxy_xxxyyy_0[i] = 3.0 * g_xxxxxx_xxxyy_1[i] * fe_0 + g_xxxxxx_xxxyyy_1[i] * pa_y[i];

        g_xxxxxxy_xxxyyz_0[i] = 2.0 * g_xxxxxx_xxxyz_1[i] * fe_0 + g_xxxxxx_xxxyyz_1[i] * pa_y[i];

        g_xxxxxxy_xxxyzz_0[i] = g_xxxxxx_xxxzz_1[i] * fe_0 + g_xxxxxx_xxxyzz_1[i] * pa_y[i];

        g_xxxxxxy_xxxzzz_0[i] = g_xxxxxx_xxxzzz_1[i] * pa_y[i];

        g_xxxxxxy_xxyyyy_0[i] = 4.0 * g_xxxxxx_xxyyy_1[i] * fe_0 + g_xxxxxx_xxyyyy_1[i] * pa_y[i];

        g_xxxxxxy_xxyyyz_0[i] = 3.0 * g_xxxxxx_xxyyz_1[i] * fe_0 + g_xxxxxx_xxyyyz_1[i] * pa_y[i];

        g_xxxxxxy_xxyyzz_0[i] = 2.0 * g_xxxxxx_xxyzz_1[i] * fe_0 + g_xxxxxx_xxyyzz_1[i] * pa_y[i];

        g_xxxxxxy_xxyzzz_0[i] = g_xxxxxx_xxzzz_1[i] * fe_0 + g_xxxxxx_xxyzzz_1[i] * pa_y[i];

        g_xxxxxxy_xxzzzz_0[i] = g_xxxxxx_xxzzzz_1[i] * pa_y[i];

        g_xxxxxxy_xyyyyy_0[i] = 5.0 * g_xxxxxx_xyyyy_1[i] * fe_0 + g_xxxxxx_xyyyyy_1[i] * pa_y[i];

        g_xxxxxxy_xyyyyz_0[i] = 4.0 * g_xxxxxx_xyyyz_1[i] * fe_0 + g_xxxxxx_xyyyyz_1[i] * pa_y[i];

        g_xxxxxxy_xyyyzz_0[i] = 3.0 * g_xxxxxx_xyyzz_1[i] * fe_0 + g_xxxxxx_xyyyzz_1[i] * pa_y[i];

        g_xxxxxxy_xyyzzz_0[i] = 2.0 * g_xxxxxx_xyzzz_1[i] * fe_0 + g_xxxxxx_xyyzzz_1[i] * pa_y[i];

        g_xxxxxxy_xyzzzz_0[i] = g_xxxxxx_xzzzz_1[i] * fe_0 + g_xxxxxx_xyzzzz_1[i] * pa_y[i];

        g_xxxxxxy_xzzzzz_0[i] = g_xxxxxx_xzzzzz_1[i] * pa_y[i];

        g_xxxxxxy_yyyyyy_0[i] = 6.0 * g_xxxxxx_yyyyy_1[i] * fe_0 + g_xxxxxx_yyyyyy_1[i] * pa_y[i];

        g_xxxxxxy_yyyyyz_0[i] = 5.0 * g_xxxxxx_yyyyz_1[i] * fe_0 + g_xxxxxx_yyyyyz_1[i] * pa_y[i];

        g_xxxxxxy_yyyyzz_0[i] = 4.0 * g_xxxxxx_yyyzz_1[i] * fe_0 + g_xxxxxx_yyyyzz_1[i] * pa_y[i];

        g_xxxxxxy_yyyzzz_0[i] = 3.0 * g_xxxxxx_yyzzz_1[i] * fe_0 + g_xxxxxx_yyyzzz_1[i] * pa_y[i];

        g_xxxxxxy_yyzzzz_0[i] = 2.0 * g_xxxxxx_yzzzz_1[i] * fe_0 + g_xxxxxx_yyzzzz_1[i] * pa_y[i];

        g_xxxxxxy_yzzzzz_0[i] = g_xxxxxx_zzzzz_1[i] * fe_0 + g_xxxxxx_yzzzzz_1[i] * pa_y[i];

        g_xxxxxxy_zzzzzz_0[i] = g_xxxxxx_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : KI

    auto g_xxxxxxz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 56);

    auto g_xxxxxxz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 57);

    auto g_xxxxxxz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 58);

    auto g_xxxxxxz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 59);

    auto g_xxxxxxz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 60);

    auto g_xxxxxxz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 61);

    auto g_xxxxxxz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 62);

    auto g_xxxxxxz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 63);

    auto g_xxxxxxz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 64);

    auto g_xxxxxxz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 65);

    auto g_xxxxxxz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 66);

    auto g_xxxxxxz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 67);

    auto g_xxxxxxz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 68);

    auto g_xxxxxxz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 69);

    auto g_xxxxxxz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 70);

    auto g_xxxxxxz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 71);

    auto g_xxxxxxz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 72);

    auto g_xxxxxxz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 73);

    auto g_xxxxxxz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 74);

    auto g_xxxxxxz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 75);

    auto g_xxxxxxz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 76);

    auto g_xxxxxxz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 77);

    auto g_xxxxxxz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 78);

    auto g_xxxxxxz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 79);

    auto g_xxxxxxz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 80);

    auto g_xxxxxxz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 81);

    auto g_xxxxxxz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 82);

    auto g_xxxxxxz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 83);

    #pragma omp simd aligned(g_xxxxxx_xxxxx_1, g_xxxxxx_xxxxxx_1, g_xxxxxx_xxxxxy_1, g_xxxxxx_xxxxxz_1, g_xxxxxx_xxxxy_1, g_xxxxxx_xxxxyy_1, g_xxxxxx_xxxxyz_1, g_xxxxxx_xxxxz_1, g_xxxxxx_xxxxzz_1, g_xxxxxx_xxxyy_1, g_xxxxxx_xxxyyy_1, g_xxxxxx_xxxyyz_1, g_xxxxxx_xxxyz_1, g_xxxxxx_xxxyzz_1, g_xxxxxx_xxxzz_1, g_xxxxxx_xxxzzz_1, g_xxxxxx_xxyyy_1, g_xxxxxx_xxyyyy_1, g_xxxxxx_xxyyyz_1, g_xxxxxx_xxyyz_1, g_xxxxxx_xxyyzz_1, g_xxxxxx_xxyzz_1, g_xxxxxx_xxyzzz_1, g_xxxxxx_xxzzz_1, g_xxxxxx_xxzzzz_1, g_xxxxxx_xyyyy_1, g_xxxxxx_xyyyyy_1, g_xxxxxx_xyyyyz_1, g_xxxxxx_xyyyz_1, g_xxxxxx_xyyyzz_1, g_xxxxxx_xyyzz_1, g_xxxxxx_xyyzzz_1, g_xxxxxx_xyzzz_1, g_xxxxxx_xyzzzz_1, g_xxxxxx_xzzzz_1, g_xxxxxx_xzzzzz_1, g_xxxxxx_yyyyy_1, g_xxxxxx_yyyyyy_1, g_xxxxxx_yyyyyz_1, g_xxxxxx_yyyyz_1, g_xxxxxx_yyyyzz_1, g_xxxxxx_yyyzz_1, g_xxxxxx_yyyzzz_1, g_xxxxxx_yyzzz_1, g_xxxxxx_yyzzzz_1, g_xxxxxx_yzzzz_1, g_xxxxxx_yzzzzz_1, g_xxxxxx_zzzzz_1, g_xxxxxx_zzzzzz_1, g_xxxxxxz_xxxxxx_0, g_xxxxxxz_xxxxxy_0, g_xxxxxxz_xxxxxz_0, g_xxxxxxz_xxxxyy_0, g_xxxxxxz_xxxxyz_0, g_xxxxxxz_xxxxzz_0, g_xxxxxxz_xxxyyy_0, g_xxxxxxz_xxxyyz_0, g_xxxxxxz_xxxyzz_0, g_xxxxxxz_xxxzzz_0, g_xxxxxxz_xxyyyy_0, g_xxxxxxz_xxyyyz_0, g_xxxxxxz_xxyyzz_0, g_xxxxxxz_xxyzzz_0, g_xxxxxxz_xxzzzz_0, g_xxxxxxz_xyyyyy_0, g_xxxxxxz_xyyyyz_0, g_xxxxxxz_xyyyzz_0, g_xxxxxxz_xyyzzz_0, g_xxxxxxz_xyzzzz_0, g_xxxxxxz_xzzzzz_0, g_xxxxxxz_yyyyyy_0, g_xxxxxxz_yyyyyz_0, g_xxxxxxz_yyyyzz_0, g_xxxxxxz_yyyzzz_0, g_xxxxxxz_yyzzzz_0, g_xxxxxxz_yzzzzz_0, g_xxxxxxz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxz_xxxxxx_0[i] = g_xxxxxx_xxxxxx_1[i] * pa_z[i];

        g_xxxxxxz_xxxxxy_0[i] = g_xxxxxx_xxxxxy_1[i] * pa_z[i];

        g_xxxxxxz_xxxxxz_0[i] = g_xxxxxx_xxxxx_1[i] * fe_0 + g_xxxxxx_xxxxxz_1[i] * pa_z[i];

        g_xxxxxxz_xxxxyy_0[i] = g_xxxxxx_xxxxyy_1[i] * pa_z[i];

        g_xxxxxxz_xxxxyz_0[i] = g_xxxxxx_xxxxy_1[i] * fe_0 + g_xxxxxx_xxxxyz_1[i] * pa_z[i];

        g_xxxxxxz_xxxxzz_0[i] = 2.0 * g_xxxxxx_xxxxz_1[i] * fe_0 + g_xxxxxx_xxxxzz_1[i] * pa_z[i];

        g_xxxxxxz_xxxyyy_0[i] = g_xxxxxx_xxxyyy_1[i] * pa_z[i];

        g_xxxxxxz_xxxyyz_0[i] = g_xxxxxx_xxxyy_1[i] * fe_0 + g_xxxxxx_xxxyyz_1[i] * pa_z[i];

        g_xxxxxxz_xxxyzz_0[i] = 2.0 * g_xxxxxx_xxxyz_1[i] * fe_0 + g_xxxxxx_xxxyzz_1[i] * pa_z[i];

        g_xxxxxxz_xxxzzz_0[i] = 3.0 * g_xxxxxx_xxxzz_1[i] * fe_0 + g_xxxxxx_xxxzzz_1[i] * pa_z[i];

        g_xxxxxxz_xxyyyy_0[i] = g_xxxxxx_xxyyyy_1[i] * pa_z[i];

        g_xxxxxxz_xxyyyz_0[i] = g_xxxxxx_xxyyy_1[i] * fe_0 + g_xxxxxx_xxyyyz_1[i] * pa_z[i];

        g_xxxxxxz_xxyyzz_0[i] = 2.0 * g_xxxxxx_xxyyz_1[i] * fe_0 + g_xxxxxx_xxyyzz_1[i] * pa_z[i];

        g_xxxxxxz_xxyzzz_0[i] = 3.0 * g_xxxxxx_xxyzz_1[i] * fe_0 + g_xxxxxx_xxyzzz_1[i] * pa_z[i];

        g_xxxxxxz_xxzzzz_0[i] = 4.0 * g_xxxxxx_xxzzz_1[i] * fe_0 + g_xxxxxx_xxzzzz_1[i] * pa_z[i];

        g_xxxxxxz_xyyyyy_0[i] = g_xxxxxx_xyyyyy_1[i] * pa_z[i];

        g_xxxxxxz_xyyyyz_0[i] = g_xxxxxx_xyyyy_1[i] * fe_0 + g_xxxxxx_xyyyyz_1[i] * pa_z[i];

        g_xxxxxxz_xyyyzz_0[i] = 2.0 * g_xxxxxx_xyyyz_1[i] * fe_0 + g_xxxxxx_xyyyzz_1[i] * pa_z[i];

        g_xxxxxxz_xyyzzz_0[i] = 3.0 * g_xxxxxx_xyyzz_1[i] * fe_0 + g_xxxxxx_xyyzzz_1[i] * pa_z[i];

        g_xxxxxxz_xyzzzz_0[i] = 4.0 * g_xxxxxx_xyzzz_1[i] * fe_0 + g_xxxxxx_xyzzzz_1[i] * pa_z[i];

        g_xxxxxxz_xzzzzz_0[i] = 5.0 * g_xxxxxx_xzzzz_1[i] * fe_0 + g_xxxxxx_xzzzzz_1[i] * pa_z[i];

        g_xxxxxxz_yyyyyy_0[i] = g_xxxxxx_yyyyyy_1[i] * pa_z[i];

        g_xxxxxxz_yyyyyz_0[i] = g_xxxxxx_yyyyy_1[i] * fe_0 + g_xxxxxx_yyyyyz_1[i] * pa_z[i];

        g_xxxxxxz_yyyyzz_0[i] = 2.0 * g_xxxxxx_yyyyz_1[i] * fe_0 + g_xxxxxx_yyyyzz_1[i] * pa_z[i];

        g_xxxxxxz_yyyzzz_0[i] = 3.0 * g_xxxxxx_yyyzz_1[i] * fe_0 + g_xxxxxx_yyyzzz_1[i] * pa_z[i];

        g_xxxxxxz_yyzzzz_0[i] = 4.0 * g_xxxxxx_yyzzz_1[i] * fe_0 + g_xxxxxx_yyzzzz_1[i] * pa_z[i];

        g_xxxxxxz_yzzzzz_0[i] = 5.0 * g_xxxxxx_yzzzz_1[i] * fe_0 + g_xxxxxx_yzzzzz_1[i] * pa_z[i];

        g_xxxxxxz_zzzzzz_0[i] = 6.0 * g_xxxxxx_zzzzz_1[i] * fe_0 + g_xxxxxx_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : KI

    auto g_xxxxxyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 84);

    auto g_xxxxxyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 85);

    auto g_xxxxxyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 86);

    auto g_xxxxxyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 87);

    auto g_xxxxxyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 88);

    auto g_xxxxxyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 89);

    auto g_xxxxxyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 90);

    auto g_xxxxxyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 91);

    auto g_xxxxxyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 92);

    auto g_xxxxxyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 93);

    auto g_xxxxxyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 94);

    auto g_xxxxxyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 95);

    auto g_xxxxxyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 96);

    auto g_xxxxxyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 97);

    auto g_xxxxxyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 98);

    auto g_xxxxxyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 99);

    auto g_xxxxxyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 100);

    auto g_xxxxxyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 101);

    auto g_xxxxxyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 102);

    auto g_xxxxxyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 103);

    auto g_xxxxxyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 104);

    auto g_xxxxxyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 105);

    auto g_xxxxxyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 106);

    auto g_xxxxxyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 107);

    auto g_xxxxxyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 108);

    auto g_xxxxxyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 109);

    auto g_xxxxxyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 110);

    auto g_xxxxxyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 111);

    #pragma omp simd aligned(g_xxxxx_xxxxxx_0, g_xxxxx_xxxxxx_1, g_xxxxx_xxxxxz_0, g_xxxxx_xxxxxz_1, g_xxxxx_xxxxzz_0, g_xxxxx_xxxxzz_1, g_xxxxx_xxxzzz_0, g_xxxxx_xxxzzz_1, g_xxxxx_xxzzzz_0, g_xxxxx_xxzzzz_1, g_xxxxx_xzzzzz_0, g_xxxxx_xzzzzz_1, g_xxxxxy_xxxxxx_1, g_xxxxxy_xxxxxz_1, g_xxxxxy_xxxxzz_1, g_xxxxxy_xxxzzz_1, g_xxxxxy_xxzzzz_1, g_xxxxxy_xzzzzz_1, g_xxxxxyy_xxxxxx_0, g_xxxxxyy_xxxxxy_0, g_xxxxxyy_xxxxxz_0, g_xxxxxyy_xxxxyy_0, g_xxxxxyy_xxxxyz_0, g_xxxxxyy_xxxxzz_0, g_xxxxxyy_xxxyyy_0, g_xxxxxyy_xxxyyz_0, g_xxxxxyy_xxxyzz_0, g_xxxxxyy_xxxzzz_0, g_xxxxxyy_xxyyyy_0, g_xxxxxyy_xxyyyz_0, g_xxxxxyy_xxyyzz_0, g_xxxxxyy_xxyzzz_0, g_xxxxxyy_xxzzzz_0, g_xxxxxyy_xyyyyy_0, g_xxxxxyy_xyyyyz_0, g_xxxxxyy_xyyyzz_0, g_xxxxxyy_xyyzzz_0, g_xxxxxyy_xyzzzz_0, g_xxxxxyy_xzzzzz_0, g_xxxxxyy_yyyyyy_0, g_xxxxxyy_yyyyyz_0, g_xxxxxyy_yyyyzz_0, g_xxxxxyy_yyyzzz_0, g_xxxxxyy_yyzzzz_0, g_xxxxxyy_yzzzzz_0, g_xxxxxyy_zzzzzz_0, g_xxxxyy_xxxxxy_1, g_xxxxyy_xxxxy_1, g_xxxxyy_xxxxyy_1, g_xxxxyy_xxxxyz_1, g_xxxxyy_xxxyy_1, g_xxxxyy_xxxyyy_1, g_xxxxyy_xxxyyz_1, g_xxxxyy_xxxyz_1, g_xxxxyy_xxxyzz_1, g_xxxxyy_xxyyy_1, g_xxxxyy_xxyyyy_1, g_xxxxyy_xxyyyz_1, g_xxxxyy_xxyyz_1, g_xxxxyy_xxyyzz_1, g_xxxxyy_xxyzz_1, g_xxxxyy_xxyzzz_1, g_xxxxyy_xyyyy_1, g_xxxxyy_xyyyyy_1, g_xxxxyy_xyyyyz_1, g_xxxxyy_xyyyz_1, g_xxxxyy_xyyyzz_1, g_xxxxyy_xyyzz_1, g_xxxxyy_xyyzzz_1, g_xxxxyy_xyzzz_1, g_xxxxyy_xyzzzz_1, g_xxxxyy_yyyyy_1, g_xxxxyy_yyyyyy_1, g_xxxxyy_yyyyyz_1, g_xxxxyy_yyyyz_1, g_xxxxyy_yyyyzz_1, g_xxxxyy_yyyzz_1, g_xxxxyy_yyyzzz_1, g_xxxxyy_yyzzz_1, g_xxxxyy_yyzzzz_1, g_xxxxyy_yzzzz_1, g_xxxxyy_yzzzzz_1, g_xxxxyy_zzzzzz_1, g_xxxyy_xxxxxy_0, g_xxxyy_xxxxxy_1, g_xxxyy_xxxxyy_0, g_xxxyy_xxxxyy_1, g_xxxyy_xxxxyz_0, g_xxxyy_xxxxyz_1, g_xxxyy_xxxyyy_0, g_xxxyy_xxxyyy_1, g_xxxyy_xxxyyz_0, g_xxxyy_xxxyyz_1, g_xxxyy_xxxyzz_0, g_xxxyy_xxxyzz_1, g_xxxyy_xxyyyy_0, g_xxxyy_xxyyyy_1, g_xxxyy_xxyyyz_0, g_xxxyy_xxyyyz_1, g_xxxyy_xxyyzz_0, g_xxxyy_xxyyzz_1, g_xxxyy_xxyzzz_0, g_xxxyy_xxyzzz_1, g_xxxyy_xyyyyy_0, g_xxxyy_xyyyyy_1, g_xxxyy_xyyyyz_0, g_xxxyy_xyyyyz_1, g_xxxyy_xyyyzz_0, g_xxxyy_xyyyzz_1, g_xxxyy_xyyzzz_0, g_xxxyy_xyyzzz_1, g_xxxyy_xyzzzz_0, g_xxxyy_xyzzzz_1, g_xxxyy_yyyyyy_0, g_xxxyy_yyyyyy_1, g_xxxyy_yyyyyz_0, g_xxxyy_yyyyyz_1, g_xxxyy_yyyyzz_0, g_xxxyy_yyyyzz_1, g_xxxyy_yyyzzz_0, g_xxxyy_yyyzzz_1, g_xxxyy_yyzzzz_0, g_xxxyy_yyzzzz_1, g_xxxyy_yzzzzz_0, g_xxxyy_yzzzzz_1, g_xxxyy_zzzzzz_0, g_xxxyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxyy_xxxxxx_0[i] = g_xxxxx_xxxxxx_0[i] * fbe_0 - g_xxxxx_xxxxxx_1[i] * fz_be_0 + g_xxxxxy_xxxxxx_1[i] * pa_y[i];

        g_xxxxxyy_xxxxxy_0[i] = 4.0 * g_xxxyy_xxxxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxxyy_xxxxy_1[i] * fe_0 + g_xxxxyy_xxxxxy_1[i] * pa_x[i];

        g_xxxxxyy_xxxxxz_0[i] = g_xxxxx_xxxxxz_0[i] * fbe_0 - g_xxxxx_xxxxxz_1[i] * fz_be_0 + g_xxxxxy_xxxxxz_1[i] * pa_y[i];

        g_xxxxxyy_xxxxyy_0[i] = 4.0 * g_xxxyy_xxxxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxxyy_xxxyy_1[i] * fe_0 + g_xxxxyy_xxxxyy_1[i] * pa_x[i];

        g_xxxxxyy_xxxxyz_0[i] = 4.0 * g_xxxyy_xxxxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_xxxyz_1[i] * fe_0 + g_xxxxyy_xxxxyz_1[i] * pa_x[i];

        g_xxxxxyy_xxxxzz_0[i] = g_xxxxx_xxxxzz_0[i] * fbe_0 - g_xxxxx_xxxxzz_1[i] * fz_be_0 + g_xxxxxy_xxxxzz_1[i] * pa_y[i];

        g_xxxxxyy_xxxyyy_0[i] = 4.0 * g_xxxyy_xxxyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_xxyyy_1[i] * fe_0 + g_xxxxyy_xxxyyy_1[i] * pa_x[i];

        g_xxxxxyy_xxxyyz_0[i] = 4.0 * g_xxxyy_xxxyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_xxyyz_1[i] * fe_0 + g_xxxxyy_xxxyyz_1[i] * pa_x[i];

        g_xxxxxyy_xxxyzz_0[i] = 4.0 * g_xxxyy_xxxyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_xxyzz_1[i] * fe_0 + g_xxxxyy_xxxyzz_1[i] * pa_x[i];

        g_xxxxxyy_xxxzzz_0[i] = g_xxxxx_xxxzzz_0[i] * fbe_0 - g_xxxxx_xxxzzz_1[i] * fz_be_0 + g_xxxxxy_xxxzzz_1[i] * pa_y[i];

        g_xxxxxyy_xxyyyy_0[i] = 4.0 * g_xxxyy_xxyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyyyy_1[i] * fe_0 + g_xxxxyy_xxyyyy_1[i] * pa_x[i];

        g_xxxxxyy_xxyyyz_0[i] = 4.0 * g_xxxyy_xxyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyyyz_1[i] * fe_0 + g_xxxxyy_xxyyyz_1[i] * pa_x[i];

        g_xxxxxyy_xxyyzz_0[i] = 4.0 * g_xxxyy_xxyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyyzz_1[i] * fe_0 + g_xxxxyy_xxyyzz_1[i] * pa_x[i];

        g_xxxxxyy_xxyzzz_0[i] = 4.0 * g_xxxyy_xxyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyzzz_1[i] * fe_0 + g_xxxxyy_xxyzzz_1[i] * pa_x[i];

        g_xxxxxyy_xxzzzz_0[i] = g_xxxxx_xxzzzz_0[i] * fbe_0 - g_xxxxx_xxzzzz_1[i] * fz_be_0 + g_xxxxxy_xxzzzz_1[i] * pa_y[i];

        g_xxxxxyy_xyyyyy_0[i] = 4.0 * g_xxxyy_xyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyyyy_1[i] * fz_be_0 + g_xxxxyy_yyyyy_1[i] * fe_0 + g_xxxxyy_xyyyyy_1[i] * pa_x[i];

        g_xxxxxyy_xyyyyz_0[i] = 4.0 * g_xxxyy_xyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyyyz_1[i] * fz_be_0 + g_xxxxyy_yyyyz_1[i] * fe_0 + g_xxxxyy_xyyyyz_1[i] * pa_x[i];

        g_xxxxxyy_xyyyzz_0[i] = 4.0 * g_xxxyy_xyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyyzz_1[i] * fz_be_0 + g_xxxxyy_yyyzz_1[i] * fe_0 + g_xxxxyy_xyyyzz_1[i] * pa_x[i];

        g_xxxxxyy_xyyzzz_0[i] = 4.0 * g_xxxyy_xyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyzzz_1[i] * fz_be_0 + g_xxxxyy_yyzzz_1[i] * fe_0 + g_xxxxyy_xyyzzz_1[i] * pa_x[i];

        g_xxxxxyy_xyzzzz_0[i] = 4.0 * g_xxxyy_xyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyzzzz_1[i] * fz_be_0 + g_xxxxyy_yzzzz_1[i] * fe_0 + g_xxxxyy_xyzzzz_1[i] * pa_x[i];

        g_xxxxxyy_xzzzzz_0[i] = g_xxxxx_xzzzzz_0[i] * fbe_0 - g_xxxxx_xzzzzz_1[i] * fz_be_0 + g_xxxxxy_xzzzzz_1[i] * pa_y[i];

        g_xxxxxyy_yyyyyy_0[i] = 4.0 * g_xxxyy_yyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyyyy_1[i] * fz_be_0 + g_xxxxyy_yyyyyy_1[i] * pa_x[i];

        g_xxxxxyy_yyyyyz_0[i] = 4.0 * g_xxxyy_yyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyyyz_1[i] * fz_be_0 + g_xxxxyy_yyyyyz_1[i] * pa_x[i];

        g_xxxxxyy_yyyyzz_0[i] = 4.0 * g_xxxyy_yyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyyzz_1[i] * fz_be_0 + g_xxxxyy_yyyyzz_1[i] * pa_x[i];

        g_xxxxxyy_yyyzzz_0[i] = 4.0 * g_xxxyy_yyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyzzz_1[i] * fz_be_0 + g_xxxxyy_yyyzzz_1[i] * pa_x[i];

        g_xxxxxyy_yyzzzz_0[i] = 4.0 * g_xxxyy_yyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyzzzz_1[i] * fz_be_0 + g_xxxxyy_yyzzzz_1[i] * pa_x[i];

        g_xxxxxyy_yzzzzz_0[i] = 4.0 * g_xxxyy_yzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yzzzzz_1[i] * fz_be_0 + g_xxxxyy_yzzzzz_1[i] * pa_x[i];

        g_xxxxxyy_zzzzzz_0[i] = 4.0 * g_xxxyy_zzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_zzzzzz_1[i] * fz_be_0 + g_xxxxyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : KI

    auto g_xxxxxyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 112);

    auto g_xxxxxyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 113);

    auto g_xxxxxyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 114);

    auto g_xxxxxyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 115);

    auto g_xxxxxyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 116);

    auto g_xxxxxyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 117);

    auto g_xxxxxyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 118);

    auto g_xxxxxyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 119);

    auto g_xxxxxyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 120);

    auto g_xxxxxyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 121);

    auto g_xxxxxyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 122);

    auto g_xxxxxyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 123);

    auto g_xxxxxyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 124);

    auto g_xxxxxyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 125);

    auto g_xxxxxyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 126);

    auto g_xxxxxyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 127);

    auto g_xxxxxyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 128);

    auto g_xxxxxyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 129);

    auto g_xxxxxyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 130);

    auto g_xxxxxyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 131);

    auto g_xxxxxyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 132);

    auto g_xxxxxyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 133);

    auto g_xxxxxyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 134);

    auto g_xxxxxyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 135);

    auto g_xxxxxyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 136);

    auto g_xxxxxyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 137);

    auto g_xxxxxyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 138);

    auto g_xxxxxyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 139);

    #pragma omp simd aligned(g_xxxxxy_xxxxxy_1, g_xxxxxy_xxxxyy_1, g_xxxxxy_xxxyyy_1, g_xxxxxy_xxyyyy_1, g_xxxxxy_xyyyyy_1, g_xxxxxy_yyyyyy_1, g_xxxxxyz_xxxxxx_0, g_xxxxxyz_xxxxxy_0, g_xxxxxyz_xxxxxz_0, g_xxxxxyz_xxxxyy_0, g_xxxxxyz_xxxxyz_0, g_xxxxxyz_xxxxzz_0, g_xxxxxyz_xxxyyy_0, g_xxxxxyz_xxxyyz_0, g_xxxxxyz_xxxyzz_0, g_xxxxxyz_xxxzzz_0, g_xxxxxyz_xxyyyy_0, g_xxxxxyz_xxyyyz_0, g_xxxxxyz_xxyyzz_0, g_xxxxxyz_xxyzzz_0, g_xxxxxyz_xxzzzz_0, g_xxxxxyz_xyyyyy_0, g_xxxxxyz_xyyyyz_0, g_xxxxxyz_xyyyzz_0, g_xxxxxyz_xyyzzz_0, g_xxxxxyz_xyzzzz_0, g_xxxxxyz_xzzzzz_0, g_xxxxxyz_yyyyyy_0, g_xxxxxyz_yyyyyz_0, g_xxxxxyz_yyyyzz_0, g_xxxxxyz_yyyzzz_0, g_xxxxxyz_yyzzzz_0, g_xxxxxyz_yzzzzz_0, g_xxxxxyz_zzzzzz_0, g_xxxxxz_xxxxxx_1, g_xxxxxz_xxxxxz_1, g_xxxxxz_xxxxyz_1, g_xxxxxz_xxxxz_1, g_xxxxxz_xxxxzz_1, g_xxxxxz_xxxyyz_1, g_xxxxxz_xxxyz_1, g_xxxxxz_xxxyzz_1, g_xxxxxz_xxxzz_1, g_xxxxxz_xxxzzz_1, g_xxxxxz_xxyyyz_1, g_xxxxxz_xxyyz_1, g_xxxxxz_xxyyzz_1, g_xxxxxz_xxyzz_1, g_xxxxxz_xxyzzz_1, g_xxxxxz_xxzzz_1, g_xxxxxz_xxzzzz_1, g_xxxxxz_xyyyyz_1, g_xxxxxz_xyyyz_1, g_xxxxxz_xyyyzz_1, g_xxxxxz_xyyzz_1, g_xxxxxz_xyyzzz_1, g_xxxxxz_xyzzz_1, g_xxxxxz_xyzzzz_1, g_xxxxxz_xzzzz_1, g_xxxxxz_xzzzzz_1, g_xxxxxz_yyyyyz_1, g_xxxxxz_yyyyz_1, g_xxxxxz_yyyyzz_1, g_xxxxxz_yyyzz_1, g_xxxxxz_yyyzzz_1, g_xxxxxz_yyzzz_1, g_xxxxxz_yyzzzz_1, g_xxxxxz_yzzzz_1, g_xxxxxz_yzzzzz_1, g_xxxxxz_zzzzz_1, g_xxxxxz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxyz_xxxxxx_0[i] = g_xxxxxz_xxxxxx_1[i] * pa_y[i];

        g_xxxxxyz_xxxxxy_0[i] = g_xxxxxy_xxxxxy_1[i] * pa_z[i];

        g_xxxxxyz_xxxxxz_0[i] = g_xxxxxz_xxxxxz_1[i] * pa_y[i];

        g_xxxxxyz_xxxxyy_0[i] = g_xxxxxy_xxxxyy_1[i] * pa_z[i];

        g_xxxxxyz_xxxxyz_0[i] = g_xxxxxz_xxxxz_1[i] * fe_0 + g_xxxxxz_xxxxyz_1[i] * pa_y[i];

        g_xxxxxyz_xxxxzz_0[i] = g_xxxxxz_xxxxzz_1[i] * pa_y[i];

        g_xxxxxyz_xxxyyy_0[i] = g_xxxxxy_xxxyyy_1[i] * pa_z[i];

        g_xxxxxyz_xxxyyz_0[i] = 2.0 * g_xxxxxz_xxxyz_1[i] * fe_0 + g_xxxxxz_xxxyyz_1[i] * pa_y[i];

        g_xxxxxyz_xxxyzz_0[i] = g_xxxxxz_xxxzz_1[i] * fe_0 + g_xxxxxz_xxxyzz_1[i] * pa_y[i];

        g_xxxxxyz_xxxzzz_0[i] = g_xxxxxz_xxxzzz_1[i] * pa_y[i];

        g_xxxxxyz_xxyyyy_0[i] = g_xxxxxy_xxyyyy_1[i] * pa_z[i];

        g_xxxxxyz_xxyyyz_0[i] = 3.0 * g_xxxxxz_xxyyz_1[i] * fe_0 + g_xxxxxz_xxyyyz_1[i] * pa_y[i];

        g_xxxxxyz_xxyyzz_0[i] = 2.0 * g_xxxxxz_xxyzz_1[i] * fe_0 + g_xxxxxz_xxyyzz_1[i] * pa_y[i];

        g_xxxxxyz_xxyzzz_0[i] = g_xxxxxz_xxzzz_1[i] * fe_0 + g_xxxxxz_xxyzzz_1[i] * pa_y[i];

        g_xxxxxyz_xxzzzz_0[i] = g_xxxxxz_xxzzzz_1[i] * pa_y[i];

        g_xxxxxyz_xyyyyy_0[i] = g_xxxxxy_xyyyyy_1[i] * pa_z[i];

        g_xxxxxyz_xyyyyz_0[i] = 4.0 * g_xxxxxz_xyyyz_1[i] * fe_0 + g_xxxxxz_xyyyyz_1[i] * pa_y[i];

        g_xxxxxyz_xyyyzz_0[i] = 3.0 * g_xxxxxz_xyyzz_1[i] * fe_0 + g_xxxxxz_xyyyzz_1[i] * pa_y[i];

        g_xxxxxyz_xyyzzz_0[i] = 2.0 * g_xxxxxz_xyzzz_1[i] * fe_0 + g_xxxxxz_xyyzzz_1[i] * pa_y[i];

        g_xxxxxyz_xyzzzz_0[i] = g_xxxxxz_xzzzz_1[i] * fe_0 + g_xxxxxz_xyzzzz_1[i] * pa_y[i];

        g_xxxxxyz_xzzzzz_0[i] = g_xxxxxz_xzzzzz_1[i] * pa_y[i];

        g_xxxxxyz_yyyyyy_0[i] = g_xxxxxy_yyyyyy_1[i] * pa_z[i];

        g_xxxxxyz_yyyyyz_0[i] = 5.0 * g_xxxxxz_yyyyz_1[i] * fe_0 + g_xxxxxz_yyyyyz_1[i] * pa_y[i];

        g_xxxxxyz_yyyyzz_0[i] = 4.0 * g_xxxxxz_yyyzz_1[i] * fe_0 + g_xxxxxz_yyyyzz_1[i] * pa_y[i];

        g_xxxxxyz_yyyzzz_0[i] = 3.0 * g_xxxxxz_yyzzz_1[i] * fe_0 + g_xxxxxz_yyyzzz_1[i] * pa_y[i];

        g_xxxxxyz_yyzzzz_0[i] = 2.0 * g_xxxxxz_yzzzz_1[i] * fe_0 + g_xxxxxz_yyzzzz_1[i] * pa_y[i];

        g_xxxxxyz_yzzzzz_0[i] = g_xxxxxz_zzzzz_1[i] * fe_0 + g_xxxxxz_yzzzzz_1[i] * pa_y[i];

        g_xxxxxyz_zzzzzz_0[i] = g_xxxxxz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : KI

    auto g_xxxxxzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 140);

    auto g_xxxxxzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 141);

    auto g_xxxxxzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 142);

    auto g_xxxxxzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 143);

    auto g_xxxxxzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 144);

    auto g_xxxxxzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 145);

    auto g_xxxxxzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 146);

    auto g_xxxxxzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 147);

    auto g_xxxxxzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 148);

    auto g_xxxxxzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 149);

    auto g_xxxxxzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 150);

    auto g_xxxxxzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 151);

    auto g_xxxxxzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 152);

    auto g_xxxxxzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 153);

    auto g_xxxxxzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 154);

    auto g_xxxxxzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 155);

    auto g_xxxxxzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 156);

    auto g_xxxxxzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 157);

    auto g_xxxxxzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 158);

    auto g_xxxxxzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 159);

    auto g_xxxxxzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 160);

    auto g_xxxxxzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 161);

    auto g_xxxxxzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 162);

    auto g_xxxxxzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 163);

    auto g_xxxxxzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 164);

    auto g_xxxxxzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 165);

    auto g_xxxxxzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 166);

    auto g_xxxxxzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 167);

    #pragma omp simd aligned(g_xxxxx_xxxxxx_0, g_xxxxx_xxxxxx_1, g_xxxxx_xxxxxy_0, g_xxxxx_xxxxxy_1, g_xxxxx_xxxxyy_0, g_xxxxx_xxxxyy_1, g_xxxxx_xxxyyy_0, g_xxxxx_xxxyyy_1, g_xxxxx_xxyyyy_0, g_xxxxx_xxyyyy_1, g_xxxxx_xyyyyy_0, g_xxxxx_xyyyyy_1, g_xxxxxz_xxxxxx_1, g_xxxxxz_xxxxxy_1, g_xxxxxz_xxxxyy_1, g_xxxxxz_xxxyyy_1, g_xxxxxz_xxyyyy_1, g_xxxxxz_xyyyyy_1, g_xxxxxzz_xxxxxx_0, g_xxxxxzz_xxxxxy_0, g_xxxxxzz_xxxxxz_0, g_xxxxxzz_xxxxyy_0, g_xxxxxzz_xxxxyz_0, g_xxxxxzz_xxxxzz_0, g_xxxxxzz_xxxyyy_0, g_xxxxxzz_xxxyyz_0, g_xxxxxzz_xxxyzz_0, g_xxxxxzz_xxxzzz_0, g_xxxxxzz_xxyyyy_0, g_xxxxxzz_xxyyyz_0, g_xxxxxzz_xxyyzz_0, g_xxxxxzz_xxyzzz_0, g_xxxxxzz_xxzzzz_0, g_xxxxxzz_xyyyyy_0, g_xxxxxzz_xyyyyz_0, g_xxxxxzz_xyyyzz_0, g_xxxxxzz_xyyzzz_0, g_xxxxxzz_xyzzzz_0, g_xxxxxzz_xzzzzz_0, g_xxxxxzz_yyyyyy_0, g_xxxxxzz_yyyyyz_0, g_xxxxxzz_yyyyzz_0, g_xxxxxzz_yyyzzz_0, g_xxxxxzz_yyzzzz_0, g_xxxxxzz_yzzzzz_0, g_xxxxxzz_zzzzzz_0, g_xxxxzz_xxxxxz_1, g_xxxxzz_xxxxyz_1, g_xxxxzz_xxxxz_1, g_xxxxzz_xxxxzz_1, g_xxxxzz_xxxyyz_1, g_xxxxzz_xxxyz_1, g_xxxxzz_xxxyzz_1, g_xxxxzz_xxxzz_1, g_xxxxzz_xxxzzz_1, g_xxxxzz_xxyyyz_1, g_xxxxzz_xxyyz_1, g_xxxxzz_xxyyzz_1, g_xxxxzz_xxyzz_1, g_xxxxzz_xxyzzz_1, g_xxxxzz_xxzzz_1, g_xxxxzz_xxzzzz_1, g_xxxxzz_xyyyyz_1, g_xxxxzz_xyyyz_1, g_xxxxzz_xyyyzz_1, g_xxxxzz_xyyzz_1, g_xxxxzz_xyyzzz_1, g_xxxxzz_xyzzz_1, g_xxxxzz_xyzzzz_1, g_xxxxzz_xzzzz_1, g_xxxxzz_xzzzzz_1, g_xxxxzz_yyyyyy_1, g_xxxxzz_yyyyyz_1, g_xxxxzz_yyyyz_1, g_xxxxzz_yyyyzz_1, g_xxxxzz_yyyzz_1, g_xxxxzz_yyyzzz_1, g_xxxxzz_yyzzz_1, g_xxxxzz_yyzzzz_1, g_xxxxzz_yzzzz_1, g_xxxxzz_yzzzzz_1, g_xxxxzz_zzzzz_1, g_xxxxzz_zzzzzz_1, g_xxxzz_xxxxxz_0, g_xxxzz_xxxxxz_1, g_xxxzz_xxxxyz_0, g_xxxzz_xxxxyz_1, g_xxxzz_xxxxzz_0, g_xxxzz_xxxxzz_1, g_xxxzz_xxxyyz_0, g_xxxzz_xxxyyz_1, g_xxxzz_xxxyzz_0, g_xxxzz_xxxyzz_1, g_xxxzz_xxxzzz_0, g_xxxzz_xxxzzz_1, g_xxxzz_xxyyyz_0, g_xxxzz_xxyyyz_1, g_xxxzz_xxyyzz_0, g_xxxzz_xxyyzz_1, g_xxxzz_xxyzzz_0, g_xxxzz_xxyzzz_1, g_xxxzz_xxzzzz_0, g_xxxzz_xxzzzz_1, g_xxxzz_xyyyyz_0, g_xxxzz_xyyyyz_1, g_xxxzz_xyyyzz_0, g_xxxzz_xyyyzz_1, g_xxxzz_xyyzzz_0, g_xxxzz_xyyzzz_1, g_xxxzz_xyzzzz_0, g_xxxzz_xyzzzz_1, g_xxxzz_xzzzzz_0, g_xxxzz_xzzzzz_1, g_xxxzz_yyyyyy_0, g_xxxzz_yyyyyy_1, g_xxxzz_yyyyyz_0, g_xxxzz_yyyyyz_1, g_xxxzz_yyyyzz_0, g_xxxzz_yyyyzz_1, g_xxxzz_yyyzzz_0, g_xxxzz_yyyzzz_1, g_xxxzz_yyzzzz_0, g_xxxzz_yyzzzz_1, g_xxxzz_yzzzzz_0, g_xxxzz_yzzzzz_1, g_xxxzz_zzzzzz_0, g_xxxzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxzz_xxxxxx_0[i] = g_xxxxx_xxxxxx_0[i] * fbe_0 - g_xxxxx_xxxxxx_1[i] * fz_be_0 + g_xxxxxz_xxxxxx_1[i] * pa_z[i];

        g_xxxxxzz_xxxxxy_0[i] = g_xxxxx_xxxxxy_0[i] * fbe_0 - g_xxxxx_xxxxxy_1[i] * fz_be_0 + g_xxxxxz_xxxxxy_1[i] * pa_z[i];

        g_xxxxxzz_xxxxxz_0[i] = 4.0 * g_xxxzz_xxxxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_xxxxz_1[i] * fe_0 + g_xxxxzz_xxxxxz_1[i] * pa_x[i];

        g_xxxxxzz_xxxxyy_0[i] = g_xxxxx_xxxxyy_0[i] * fbe_0 - g_xxxxx_xxxxyy_1[i] * fz_be_0 + g_xxxxxz_xxxxyy_1[i] * pa_z[i];

        g_xxxxxzz_xxxxyz_0[i] = 4.0 * g_xxxzz_xxxxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_xxxyz_1[i] * fe_0 + g_xxxxzz_xxxxyz_1[i] * pa_x[i];

        g_xxxxxzz_xxxxzz_0[i] = 4.0 * g_xxxzz_xxxxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_xxxzz_1[i] * fe_0 + g_xxxxzz_xxxxzz_1[i] * pa_x[i];

        g_xxxxxzz_xxxyyy_0[i] = g_xxxxx_xxxyyy_0[i] * fbe_0 - g_xxxxx_xxxyyy_1[i] * fz_be_0 + g_xxxxxz_xxxyyy_1[i] * pa_z[i];

        g_xxxxxzz_xxxyyz_0[i] = 4.0 * g_xxxzz_xxxyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_xxyyz_1[i] * fe_0 + g_xxxxzz_xxxyyz_1[i] * pa_x[i];

        g_xxxxxzz_xxxyzz_0[i] = 4.0 * g_xxxzz_xxxyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_xxyzz_1[i] * fe_0 + g_xxxxzz_xxxyzz_1[i] * pa_x[i];

        g_xxxxxzz_xxxzzz_0[i] = 4.0 * g_xxxzz_xxxzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_xxzzz_1[i] * fe_0 + g_xxxxzz_xxxzzz_1[i] * pa_x[i];

        g_xxxxxzz_xxyyyy_0[i] = g_xxxxx_xxyyyy_0[i] * fbe_0 - g_xxxxx_xxyyyy_1[i] * fz_be_0 + g_xxxxxz_xxyyyy_1[i] * pa_z[i];

        g_xxxxxzz_xxyyyz_0[i] = 4.0 * g_xxxzz_xxyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xyyyz_1[i] * fe_0 + g_xxxxzz_xxyyyz_1[i] * pa_x[i];

        g_xxxxxzz_xxyyzz_0[i] = 4.0 * g_xxxzz_xxyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xyyzz_1[i] * fe_0 + g_xxxxzz_xxyyzz_1[i] * pa_x[i];

        g_xxxxxzz_xxyzzz_0[i] = 4.0 * g_xxxzz_xxyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xyzzz_1[i] * fe_0 + g_xxxxzz_xxyzzz_1[i] * pa_x[i];

        g_xxxxxzz_xxzzzz_0[i] = 4.0 * g_xxxzz_xxzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xzzzz_1[i] * fe_0 + g_xxxxzz_xxzzzz_1[i] * pa_x[i];

        g_xxxxxzz_xyyyyy_0[i] = g_xxxxx_xyyyyy_0[i] * fbe_0 - g_xxxxx_xyyyyy_1[i] * fz_be_0 + g_xxxxxz_xyyyyy_1[i] * pa_z[i];

        g_xxxxxzz_xyyyyz_0[i] = 4.0 * g_xxxzz_xyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyyyyz_1[i] * fz_be_0 + g_xxxxzz_yyyyz_1[i] * fe_0 + g_xxxxzz_xyyyyz_1[i] * pa_x[i];

        g_xxxxxzz_xyyyzz_0[i] = 4.0 * g_xxxzz_xyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyyyzz_1[i] * fz_be_0 + g_xxxxzz_yyyzz_1[i] * fe_0 + g_xxxxzz_xyyyzz_1[i] * pa_x[i];

        g_xxxxxzz_xyyzzz_0[i] = 4.0 * g_xxxzz_xyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyyzzz_1[i] * fz_be_0 + g_xxxxzz_yyzzz_1[i] * fe_0 + g_xxxxzz_xyyzzz_1[i] * pa_x[i];

        g_xxxxxzz_xyzzzz_0[i] = 4.0 * g_xxxzz_xyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyzzzz_1[i] * fz_be_0 + g_xxxxzz_yzzzz_1[i] * fe_0 + g_xxxxzz_xyzzzz_1[i] * pa_x[i];

        g_xxxxxzz_xzzzzz_0[i] = 4.0 * g_xxxzz_xzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xzzzzz_1[i] * fz_be_0 + g_xxxxzz_zzzzz_1[i] * fe_0 + g_xxxxzz_xzzzzz_1[i] * pa_x[i];

        g_xxxxxzz_yyyyyy_0[i] = 4.0 * g_xxxzz_yyyyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyyyy_1[i] * fz_be_0 + g_xxxxzz_yyyyyy_1[i] * pa_x[i];

        g_xxxxxzz_yyyyyz_0[i] = 4.0 * g_xxxzz_yyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyyyz_1[i] * fz_be_0 + g_xxxxzz_yyyyyz_1[i] * pa_x[i];

        g_xxxxxzz_yyyyzz_0[i] = 4.0 * g_xxxzz_yyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyyzz_1[i] * fz_be_0 + g_xxxxzz_yyyyzz_1[i] * pa_x[i];

        g_xxxxxzz_yyyzzz_0[i] = 4.0 * g_xxxzz_yyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyzzz_1[i] * fz_be_0 + g_xxxxzz_yyyzzz_1[i] * pa_x[i];

        g_xxxxxzz_yyzzzz_0[i] = 4.0 * g_xxxzz_yyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyzzzz_1[i] * fz_be_0 + g_xxxxzz_yyzzzz_1[i] * pa_x[i];

        g_xxxxxzz_yzzzzz_0[i] = 4.0 * g_xxxzz_yzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yzzzzz_1[i] * fz_be_0 + g_xxxxzz_yzzzzz_1[i] * pa_x[i];

        g_xxxxxzz_zzzzzz_0[i] = 4.0 * g_xxxzz_zzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_zzzzzz_1[i] * fz_be_0 + g_xxxxzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : KI

    auto g_xxxxyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 168);

    auto g_xxxxyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 169);

    auto g_xxxxyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 170);

    auto g_xxxxyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 171);

    auto g_xxxxyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 172);

    auto g_xxxxyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 173);

    auto g_xxxxyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 174);

    auto g_xxxxyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 175);

    auto g_xxxxyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 176);

    auto g_xxxxyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 177);

    auto g_xxxxyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 178);

    auto g_xxxxyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 179);

    auto g_xxxxyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 180);

    auto g_xxxxyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 181);

    auto g_xxxxyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 182);

    auto g_xxxxyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 183);

    auto g_xxxxyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 184);

    auto g_xxxxyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 185);

    auto g_xxxxyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 186);

    auto g_xxxxyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 187);

    auto g_xxxxyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 188);

    auto g_xxxxyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 189);

    auto g_xxxxyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 190);

    auto g_xxxxyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 191);

    auto g_xxxxyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 192);

    auto g_xxxxyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 193);

    auto g_xxxxyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 194);

    auto g_xxxxyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 195);

    #pragma omp simd aligned(g_xxxxy_xxxxxx_0, g_xxxxy_xxxxxx_1, g_xxxxy_xxxxxz_0, g_xxxxy_xxxxxz_1, g_xxxxy_xxxxzz_0, g_xxxxy_xxxxzz_1, g_xxxxy_xxxzzz_0, g_xxxxy_xxxzzz_1, g_xxxxy_xxzzzz_0, g_xxxxy_xxzzzz_1, g_xxxxy_xzzzzz_0, g_xxxxy_xzzzzz_1, g_xxxxyy_xxxxxx_1, g_xxxxyy_xxxxxz_1, g_xxxxyy_xxxxzz_1, g_xxxxyy_xxxzzz_1, g_xxxxyy_xxzzzz_1, g_xxxxyy_xzzzzz_1, g_xxxxyyy_xxxxxx_0, g_xxxxyyy_xxxxxy_0, g_xxxxyyy_xxxxxz_0, g_xxxxyyy_xxxxyy_0, g_xxxxyyy_xxxxyz_0, g_xxxxyyy_xxxxzz_0, g_xxxxyyy_xxxyyy_0, g_xxxxyyy_xxxyyz_0, g_xxxxyyy_xxxyzz_0, g_xxxxyyy_xxxzzz_0, g_xxxxyyy_xxyyyy_0, g_xxxxyyy_xxyyyz_0, g_xxxxyyy_xxyyzz_0, g_xxxxyyy_xxyzzz_0, g_xxxxyyy_xxzzzz_0, g_xxxxyyy_xyyyyy_0, g_xxxxyyy_xyyyyz_0, g_xxxxyyy_xyyyzz_0, g_xxxxyyy_xyyzzz_0, g_xxxxyyy_xyzzzz_0, g_xxxxyyy_xzzzzz_0, g_xxxxyyy_yyyyyy_0, g_xxxxyyy_yyyyyz_0, g_xxxxyyy_yyyyzz_0, g_xxxxyyy_yyyzzz_0, g_xxxxyyy_yyzzzz_0, g_xxxxyyy_yzzzzz_0, g_xxxxyyy_zzzzzz_0, g_xxxyyy_xxxxxy_1, g_xxxyyy_xxxxy_1, g_xxxyyy_xxxxyy_1, g_xxxyyy_xxxxyz_1, g_xxxyyy_xxxyy_1, g_xxxyyy_xxxyyy_1, g_xxxyyy_xxxyyz_1, g_xxxyyy_xxxyz_1, g_xxxyyy_xxxyzz_1, g_xxxyyy_xxyyy_1, g_xxxyyy_xxyyyy_1, g_xxxyyy_xxyyyz_1, g_xxxyyy_xxyyz_1, g_xxxyyy_xxyyzz_1, g_xxxyyy_xxyzz_1, g_xxxyyy_xxyzzz_1, g_xxxyyy_xyyyy_1, g_xxxyyy_xyyyyy_1, g_xxxyyy_xyyyyz_1, g_xxxyyy_xyyyz_1, g_xxxyyy_xyyyzz_1, g_xxxyyy_xyyzz_1, g_xxxyyy_xyyzzz_1, g_xxxyyy_xyzzz_1, g_xxxyyy_xyzzzz_1, g_xxxyyy_yyyyy_1, g_xxxyyy_yyyyyy_1, g_xxxyyy_yyyyyz_1, g_xxxyyy_yyyyz_1, g_xxxyyy_yyyyzz_1, g_xxxyyy_yyyzz_1, g_xxxyyy_yyyzzz_1, g_xxxyyy_yyzzz_1, g_xxxyyy_yyzzzz_1, g_xxxyyy_yzzzz_1, g_xxxyyy_yzzzzz_1, g_xxxyyy_zzzzzz_1, g_xxyyy_xxxxxy_0, g_xxyyy_xxxxxy_1, g_xxyyy_xxxxyy_0, g_xxyyy_xxxxyy_1, g_xxyyy_xxxxyz_0, g_xxyyy_xxxxyz_1, g_xxyyy_xxxyyy_0, g_xxyyy_xxxyyy_1, g_xxyyy_xxxyyz_0, g_xxyyy_xxxyyz_1, g_xxyyy_xxxyzz_0, g_xxyyy_xxxyzz_1, g_xxyyy_xxyyyy_0, g_xxyyy_xxyyyy_1, g_xxyyy_xxyyyz_0, g_xxyyy_xxyyyz_1, g_xxyyy_xxyyzz_0, g_xxyyy_xxyyzz_1, g_xxyyy_xxyzzz_0, g_xxyyy_xxyzzz_1, g_xxyyy_xyyyyy_0, g_xxyyy_xyyyyy_1, g_xxyyy_xyyyyz_0, g_xxyyy_xyyyyz_1, g_xxyyy_xyyyzz_0, g_xxyyy_xyyyzz_1, g_xxyyy_xyyzzz_0, g_xxyyy_xyyzzz_1, g_xxyyy_xyzzzz_0, g_xxyyy_xyzzzz_1, g_xxyyy_yyyyyy_0, g_xxyyy_yyyyyy_1, g_xxyyy_yyyyyz_0, g_xxyyy_yyyyyz_1, g_xxyyy_yyyyzz_0, g_xxyyy_yyyyzz_1, g_xxyyy_yyyzzz_0, g_xxyyy_yyyzzz_1, g_xxyyy_yyzzzz_0, g_xxyyy_yyzzzz_1, g_xxyyy_yzzzzz_0, g_xxyyy_yzzzzz_1, g_xxyyy_zzzzzz_0, g_xxyyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyyy_xxxxxx_0[i] = 2.0 * g_xxxxy_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxxxx_1[i] * fz_be_0 + g_xxxxyy_xxxxxx_1[i] * pa_y[i];

        g_xxxxyyy_xxxxxy_0[i] = 3.0 * g_xxyyy_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxyyy_xxxxy_1[i] * fe_0 + g_xxxyyy_xxxxxy_1[i] * pa_x[i];

        g_xxxxyyy_xxxxxz_0[i] = 2.0 * g_xxxxy_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxxxz_1[i] * fz_be_0 + g_xxxxyy_xxxxxz_1[i] * pa_y[i];

        g_xxxxyyy_xxxxyy_0[i] = 3.0 * g_xxyyy_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxyyy_xxxyy_1[i] * fe_0 + g_xxxyyy_xxxxyy_1[i] * pa_x[i];

        g_xxxxyyy_xxxxyz_0[i] = 3.0 * g_xxyyy_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_xxxyz_1[i] * fe_0 + g_xxxyyy_xxxxyz_1[i] * pa_x[i];

        g_xxxxyyy_xxxxzz_0[i] = 2.0 * g_xxxxy_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxxzz_1[i] * fz_be_0 + g_xxxxyy_xxxxzz_1[i] * pa_y[i];

        g_xxxxyyy_xxxyyy_0[i] = 3.0 * g_xxyyy_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_xxyyy_1[i] * fe_0 + g_xxxyyy_xxxyyy_1[i] * pa_x[i];

        g_xxxxyyy_xxxyyz_0[i] = 3.0 * g_xxyyy_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_xxyyz_1[i] * fe_0 + g_xxxyyy_xxxyyz_1[i] * pa_x[i];

        g_xxxxyyy_xxxyzz_0[i] = 3.0 * g_xxyyy_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_xxyzz_1[i] * fe_0 + g_xxxyyy_xxxyzz_1[i] * pa_x[i];

        g_xxxxyyy_xxxzzz_0[i] = 2.0 * g_xxxxy_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxzzz_1[i] * fz_be_0 + g_xxxxyy_xxxzzz_1[i] * pa_y[i];

        g_xxxxyyy_xxyyyy_0[i] = 3.0 * g_xxyyy_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyyyy_1[i] * fe_0 + g_xxxyyy_xxyyyy_1[i] * pa_x[i];

        g_xxxxyyy_xxyyyz_0[i] = 3.0 * g_xxyyy_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyyyz_1[i] * fe_0 + g_xxxyyy_xxyyyz_1[i] * pa_x[i];

        g_xxxxyyy_xxyyzz_0[i] = 3.0 * g_xxyyy_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyyzz_1[i] * fe_0 + g_xxxyyy_xxyyzz_1[i] * pa_x[i];

        g_xxxxyyy_xxyzzz_0[i] = 3.0 * g_xxyyy_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyzzz_1[i] * fe_0 + g_xxxyyy_xxyzzz_1[i] * pa_x[i];

        g_xxxxyyy_xxzzzz_0[i] = 2.0 * g_xxxxy_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxzzzz_1[i] * fz_be_0 + g_xxxxyy_xxzzzz_1[i] * pa_y[i];

        g_xxxxyyy_xyyyyy_0[i] = 3.0 * g_xxyyy_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyyyy_1[i] * fz_be_0 + g_xxxyyy_yyyyy_1[i] * fe_0 + g_xxxyyy_xyyyyy_1[i] * pa_x[i];

        g_xxxxyyy_xyyyyz_0[i] = 3.0 * g_xxyyy_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyyyz_1[i] * fz_be_0 + g_xxxyyy_yyyyz_1[i] * fe_0 + g_xxxyyy_xyyyyz_1[i] * pa_x[i];

        g_xxxxyyy_xyyyzz_0[i] = 3.0 * g_xxyyy_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyyzz_1[i] * fz_be_0 + g_xxxyyy_yyyzz_1[i] * fe_0 + g_xxxyyy_xyyyzz_1[i] * pa_x[i];

        g_xxxxyyy_xyyzzz_0[i] = 3.0 * g_xxyyy_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyzzz_1[i] * fz_be_0 + g_xxxyyy_yyzzz_1[i] * fe_0 + g_xxxyyy_xyyzzz_1[i] * pa_x[i];

        g_xxxxyyy_xyzzzz_0[i] = 3.0 * g_xxyyy_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyzzzz_1[i] * fz_be_0 + g_xxxyyy_yzzzz_1[i] * fe_0 + g_xxxyyy_xyzzzz_1[i] * pa_x[i];

        g_xxxxyyy_xzzzzz_0[i] = 2.0 * g_xxxxy_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xzzzzz_1[i] * fz_be_0 + g_xxxxyy_xzzzzz_1[i] * pa_y[i];

        g_xxxxyyy_yyyyyy_0[i] = 3.0 * g_xxyyy_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyyyy_1[i] * fz_be_0 + g_xxxyyy_yyyyyy_1[i] * pa_x[i];

        g_xxxxyyy_yyyyyz_0[i] = 3.0 * g_xxyyy_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyyyz_1[i] * fz_be_0 + g_xxxyyy_yyyyyz_1[i] * pa_x[i];

        g_xxxxyyy_yyyyzz_0[i] = 3.0 * g_xxyyy_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyyzz_1[i] * fz_be_0 + g_xxxyyy_yyyyzz_1[i] * pa_x[i];

        g_xxxxyyy_yyyzzz_0[i] = 3.0 * g_xxyyy_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyzzz_1[i] * fz_be_0 + g_xxxyyy_yyyzzz_1[i] * pa_x[i];

        g_xxxxyyy_yyzzzz_0[i] = 3.0 * g_xxyyy_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyzzzz_1[i] * fz_be_0 + g_xxxyyy_yyzzzz_1[i] * pa_x[i];

        g_xxxxyyy_yzzzzz_0[i] = 3.0 * g_xxyyy_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yzzzzz_1[i] * fz_be_0 + g_xxxyyy_yzzzzz_1[i] * pa_x[i];

        g_xxxxyyy_zzzzzz_0[i] = 3.0 * g_xxyyy_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_zzzzzz_1[i] * fz_be_0 + g_xxxyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : KI

    auto g_xxxxyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 196);

    auto g_xxxxyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 197);

    auto g_xxxxyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 198);

    auto g_xxxxyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 199);

    auto g_xxxxyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 200);

    auto g_xxxxyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 201);

    auto g_xxxxyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 202);

    auto g_xxxxyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 203);

    auto g_xxxxyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 204);

    auto g_xxxxyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 205);

    auto g_xxxxyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 206);

    auto g_xxxxyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 207);

    auto g_xxxxyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 208);

    auto g_xxxxyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 209);

    auto g_xxxxyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 210);

    auto g_xxxxyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 211);

    auto g_xxxxyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 212);

    auto g_xxxxyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 213);

    auto g_xxxxyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 214);

    auto g_xxxxyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 215);

    auto g_xxxxyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 216);

    auto g_xxxxyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 217);

    auto g_xxxxyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 218);

    auto g_xxxxyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 219);

    auto g_xxxxyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 220);

    auto g_xxxxyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 221);

    auto g_xxxxyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 222);

    auto g_xxxxyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 223);

    #pragma omp simd aligned(g_xxxxyy_xxxxx_1, g_xxxxyy_xxxxxx_1, g_xxxxyy_xxxxxy_1, g_xxxxyy_xxxxxz_1, g_xxxxyy_xxxxy_1, g_xxxxyy_xxxxyy_1, g_xxxxyy_xxxxyz_1, g_xxxxyy_xxxxz_1, g_xxxxyy_xxxxzz_1, g_xxxxyy_xxxyy_1, g_xxxxyy_xxxyyy_1, g_xxxxyy_xxxyyz_1, g_xxxxyy_xxxyz_1, g_xxxxyy_xxxyzz_1, g_xxxxyy_xxxzz_1, g_xxxxyy_xxxzzz_1, g_xxxxyy_xxyyy_1, g_xxxxyy_xxyyyy_1, g_xxxxyy_xxyyyz_1, g_xxxxyy_xxyyz_1, g_xxxxyy_xxyyzz_1, g_xxxxyy_xxyzz_1, g_xxxxyy_xxyzzz_1, g_xxxxyy_xxzzz_1, g_xxxxyy_xxzzzz_1, g_xxxxyy_xyyyy_1, g_xxxxyy_xyyyyy_1, g_xxxxyy_xyyyyz_1, g_xxxxyy_xyyyz_1, g_xxxxyy_xyyyzz_1, g_xxxxyy_xyyzz_1, g_xxxxyy_xyyzzz_1, g_xxxxyy_xyzzz_1, g_xxxxyy_xyzzzz_1, g_xxxxyy_xzzzz_1, g_xxxxyy_xzzzzz_1, g_xxxxyy_yyyyy_1, g_xxxxyy_yyyyyy_1, g_xxxxyy_yyyyyz_1, g_xxxxyy_yyyyz_1, g_xxxxyy_yyyyzz_1, g_xxxxyy_yyyzz_1, g_xxxxyy_yyyzzz_1, g_xxxxyy_yyzzz_1, g_xxxxyy_yyzzzz_1, g_xxxxyy_yzzzz_1, g_xxxxyy_yzzzzz_1, g_xxxxyy_zzzzz_1, g_xxxxyy_zzzzzz_1, g_xxxxyyz_xxxxxx_0, g_xxxxyyz_xxxxxy_0, g_xxxxyyz_xxxxxz_0, g_xxxxyyz_xxxxyy_0, g_xxxxyyz_xxxxyz_0, g_xxxxyyz_xxxxzz_0, g_xxxxyyz_xxxyyy_0, g_xxxxyyz_xxxyyz_0, g_xxxxyyz_xxxyzz_0, g_xxxxyyz_xxxzzz_0, g_xxxxyyz_xxyyyy_0, g_xxxxyyz_xxyyyz_0, g_xxxxyyz_xxyyzz_0, g_xxxxyyz_xxyzzz_0, g_xxxxyyz_xxzzzz_0, g_xxxxyyz_xyyyyy_0, g_xxxxyyz_xyyyyz_0, g_xxxxyyz_xyyyzz_0, g_xxxxyyz_xyyzzz_0, g_xxxxyyz_xyzzzz_0, g_xxxxyyz_xzzzzz_0, g_xxxxyyz_yyyyyy_0, g_xxxxyyz_yyyyyz_0, g_xxxxyyz_yyyyzz_0, g_xxxxyyz_yyyzzz_0, g_xxxxyyz_yyzzzz_0, g_xxxxyyz_yzzzzz_0, g_xxxxyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyyz_xxxxxx_0[i] = g_xxxxyy_xxxxxx_1[i] * pa_z[i];

        g_xxxxyyz_xxxxxy_0[i] = g_xxxxyy_xxxxxy_1[i] * pa_z[i];

        g_xxxxyyz_xxxxxz_0[i] = g_xxxxyy_xxxxx_1[i] * fe_0 + g_xxxxyy_xxxxxz_1[i] * pa_z[i];

        g_xxxxyyz_xxxxyy_0[i] = g_xxxxyy_xxxxyy_1[i] * pa_z[i];

        g_xxxxyyz_xxxxyz_0[i] = g_xxxxyy_xxxxy_1[i] * fe_0 + g_xxxxyy_xxxxyz_1[i] * pa_z[i];

        g_xxxxyyz_xxxxzz_0[i] = 2.0 * g_xxxxyy_xxxxz_1[i] * fe_0 + g_xxxxyy_xxxxzz_1[i] * pa_z[i];

        g_xxxxyyz_xxxyyy_0[i] = g_xxxxyy_xxxyyy_1[i] * pa_z[i];

        g_xxxxyyz_xxxyyz_0[i] = g_xxxxyy_xxxyy_1[i] * fe_0 + g_xxxxyy_xxxyyz_1[i] * pa_z[i];

        g_xxxxyyz_xxxyzz_0[i] = 2.0 * g_xxxxyy_xxxyz_1[i] * fe_0 + g_xxxxyy_xxxyzz_1[i] * pa_z[i];

        g_xxxxyyz_xxxzzz_0[i] = 3.0 * g_xxxxyy_xxxzz_1[i] * fe_0 + g_xxxxyy_xxxzzz_1[i] * pa_z[i];

        g_xxxxyyz_xxyyyy_0[i] = g_xxxxyy_xxyyyy_1[i] * pa_z[i];

        g_xxxxyyz_xxyyyz_0[i] = g_xxxxyy_xxyyy_1[i] * fe_0 + g_xxxxyy_xxyyyz_1[i] * pa_z[i];

        g_xxxxyyz_xxyyzz_0[i] = 2.0 * g_xxxxyy_xxyyz_1[i] * fe_0 + g_xxxxyy_xxyyzz_1[i] * pa_z[i];

        g_xxxxyyz_xxyzzz_0[i] = 3.0 * g_xxxxyy_xxyzz_1[i] * fe_0 + g_xxxxyy_xxyzzz_1[i] * pa_z[i];

        g_xxxxyyz_xxzzzz_0[i] = 4.0 * g_xxxxyy_xxzzz_1[i] * fe_0 + g_xxxxyy_xxzzzz_1[i] * pa_z[i];

        g_xxxxyyz_xyyyyy_0[i] = g_xxxxyy_xyyyyy_1[i] * pa_z[i];

        g_xxxxyyz_xyyyyz_0[i] = g_xxxxyy_xyyyy_1[i] * fe_0 + g_xxxxyy_xyyyyz_1[i] * pa_z[i];

        g_xxxxyyz_xyyyzz_0[i] = 2.0 * g_xxxxyy_xyyyz_1[i] * fe_0 + g_xxxxyy_xyyyzz_1[i] * pa_z[i];

        g_xxxxyyz_xyyzzz_0[i] = 3.0 * g_xxxxyy_xyyzz_1[i] * fe_0 + g_xxxxyy_xyyzzz_1[i] * pa_z[i];

        g_xxxxyyz_xyzzzz_0[i] = 4.0 * g_xxxxyy_xyzzz_1[i] * fe_0 + g_xxxxyy_xyzzzz_1[i] * pa_z[i];

        g_xxxxyyz_xzzzzz_0[i] = 5.0 * g_xxxxyy_xzzzz_1[i] * fe_0 + g_xxxxyy_xzzzzz_1[i] * pa_z[i];

        g_xxxxyyz_yyyyyy_0[i] = g_xxxxyy_yyyyyy_1[i] * pa_z[i];

        g_xxxxyyz_yyyyyz_0[i] = g_xxxxyy_yyyyy_1[i] * fe_0 + g_xxxxyy_yyyyyz_1[i] * pa_z[i];

        g_xxxxyyz_yyyyzz_0[i] = 2.0 * g_xxxxyy_yyyyz_1[i] * fe_0 + g_xxxxyy_yyyyzz_1[i] * pa_z[i];

        g_xxxxyyz_yyyzzz_0[i] = 3.0 * g_xxxxyy_yyyzz_1[i] * fe_0 + g_xxxxyy_yyyzzz_1[i] * pa_z[i];

        g_xxxxyyz_yyzzzz_0[i] = 4.0 * g_xxxxyy_yyzzz_1[i] * fe_0 + g_xxxxyy_yyzzzz_1[i] * pa_z[i];

        g_xxxxyyz_yzzzzz_0[i] = 5.0 * g_xxxxyy_yzzzz_1[i] * fe_0 + g_xxxxyy_yzzzzz_1[i] * pa_z[i];

        g_xxxxyyz_zzzzzz_0[i] = 6.0 * g_xxxxyy_zzzzz_1[i] * fe_0 + g_xxxxyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 224-252 components of targeted buffer : KI

    auto g_xxxxyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 224);

    auto g_xxxxyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 225);

    auto g_xxxxyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 226);

    auto g_xxxxyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 227);

    auto g_xxxxyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 228);

    auto g_xxxxyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 229);

    auto g_xxxxyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 230);

    auto g_xxxxyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 231);

    auto g_xxxxyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 232);

    auto g_xxxxyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 233);

    auto g_xxxxyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 234);

    auto g_xxxxyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 235);

    auto g_xxxxyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 236);

    auto g_xxxxyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 237);

    auto g_xxxxyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 238);

    auto g_xxxxyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 239);

    auto g_xxxxyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 240);

    auto g_xxxxyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 241);

    auto g_xxxxyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 242);

    auto g_xxxxyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 243);

    auto g_xxxxyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 244);

    auto g_xxxxyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 245);

    auto g_xxxxyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 246);

    auto g_xxxxyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 247);

    auto g_xxxxyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 248);

    auto g_xxxxyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 249);

    auto g_xxxxyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 250);

    auto g_xxxxyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 251);

    #pragma omp simd aligned(g_xxxxyzz_xxxxxx_0, g_xxxxyzz_xxxxxy_0, g_xxxxyzz_xxxxxz_0, g_xxxxyzz_xxxxyy_0, g_xxxxyzz_xxxxyz_0, g_xxxxyzz_xxxxzz_0, g_xxxxyzz_xxxyyy_0, g_xxxxyzz_xxxyyz_0, g_xxxxyzz_xxxyzz_0, g_xxxxyzz_xxxzzz_0, g_xxxxyzz_xxyyyy_0, g_xxxxyzz_xxyyyz_0, g_xxxxyzz_xxyyzz_0, g_xxxxyzz_xxyzzz_0, g_xxxxyzz_xxzzzz_0, g_xxxxyzz_xyyyyy_0, g_xxxxyzz_xyyyyz_0, g_xxxxyzz_xyyyzz_0, g_xxxxyzz_xyyzzz_0, g_xxxxyzz_xyzzzz_0, g_xxxxyzz_xzzzzz_0, g_xxxxyzz_yyyyyy_0, g_xxxxyzz_yyyyyz_0, g_xxxxyzz_yyyyzz_0, g_xxxxyzz_yyyzzz_0, g_xxxxyzz_yyzzzz_0, g_xxxxyzz_yzzzzz_0, g_xxxxyzz_zzzzzz_0, g_xxxxzz_xxxxx_1, g_xxxxzz_xxxxxx_1, g_xxxxzz_xxxxxy_1, g_xxxxzz_xxxxxz_1, g_xxxxzz_xxxxy_1, g_xxxxzz_xxxxyy_1, g_xxxxzz_xxxxyz_1, g_xxxxzz_xxxxz_1, g_xxxxzz_xxxxzz_1, g_xxxxzz_xxxyy_1, g_xxxxzz_xxxyyy_1, g_xxxxzz_xxxyyz_1, g_xxxxzz_xxxyz_1, g_xxxxzz_xxxyzz_1, g_xxxxzz_xxxzz_1, g_xxxxzz_xxxzzz_1, g_xxxxzz_xxyyy_1, g_xxxxzz_xxyyyy_1, g_xxxxzz_xxyyyz_1, g_xxxxzz_xxyyz_1, g_xxxxzz_xxyyzz_1, g_xxxxzz_xxyzz_1, g_xxxxzz_xxyzzz_1, g_xxxxzz_xxzzz_1, g_xxxxzz_xxzzzz_1, g_xxxxzz_xyyyy_1, g_xxxxzz_xyyyyy_1, g_xxxxzz_xyyyyz_1, g_xxxxzz_xyyyz_1, g_xxxxzz_xyyyzz_1, g_xxxxzz_xyyzz_1, g_xxxxzz_xyyzzz_1, g_xxxxzz_xyzzz_1, g_xxxxzz_xyzzzz_1, g_xxxxzz_xzzzz_1, g_xxxxzz_xzzzzz_1, g_xxxxzz_yyyyy_1, g_xxxxzz_yyyyyy_1, g_xxxxzz_yyyyyz_1, g_xxxxzz_yyyyz_1, g_xxxxzz_yyyyzz_1, g_xxxxzz_yyyzz_1, g_xxxxzz_yyyzzz_1, g_xxxxzz_yyzzz_1, g_xxxxzz_yyzzzz_1, g_xxxxzz_yzzzz_1, g_xxxxzz_yzzzzz_1, g_xxxxzz_zzzzz_1, g_xxxxzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyzz_xxxxxx_0[i] = g_xxxxzz_xxxxxx_1[i] * pa_y[i];

        g_xxxxyzz_xxxxxy_0[i] = g_xxxxzz_xxxxx_1[i] * fe_0 + g_xxxxzz_xxxxxy_1[i] * pa_y[i];

        g_xxxxyzz_xxxxxz_0[i] = g_xxxxzz_xxxxxz_1[i] * pa_y[i];

        g_xxxxyzz_xxxxyy_0[i] = 2.0 * g_xxxxzz_xxxxy_1[i] * fe_0 + g_xxxxzz_xxxxyy_1[i] * pa_y[i];

        g_xxxxyzz_xxxxyz_0[i] = g_xxxxzz_xxxxz_1[i] * fe_0 + g_xxxxzz_xxxxyz_1[i] * pa_y[i];

        g_xxxxyzz_xxxxzz_0[i] = g_xxxxzz_xxxxzz_1[i] * pa_y[i];

        g_xxxxyzz_xxxyyy_0[i] = 3.0 * g_xxxxzz_xxxyy_1[i] * fe_0 + g_xxxxzz_xxxyyy_1[i] * pa_y[i];

        g_xxxxyzz_xxxyyz_0[i] = 2.0 * g_xxxxzz_xxxyz_1[i] * fe_0 + g_xxxxzz_xxxyyz_1[i] * pa_y[i];

        g_xxxxyzz_xxxyzz_0[i] = g_xxxxzz_xxxzz_1[i] * fe_0 + g_xxxxzz_xxxyzz_1[i] * pa_y[i];

        g_xxxxyzz_xxxzzz_0[i] = g_xxxxzz_xxxzzz_1[i] * pa_y[i];

        g_xxxxyzz_xxyyyy_0[i] = 4.0 * g_xxxxzz_xxyyy_1[i] * fe_0 + g_xxxxzz_xxyyyy_1[i] * pa_y[i];

        g_xxxxyzz_xxyyyz_0[i] = 3.0 * g_xxxxzz_xxyyz_1[i] * fe_0 + g_xxxxzz_xxyyyz_1[i] * pa_y[i];

        g_xxxxyzz_xxyyzz_0[i] = 2.0 * g_xxxxzz_xxyzz_1[i] * fe_0 + g_xxxxzz_xxyyzz_1[i] * pa_y[i];

        g_xxxxyzz_xxyzzz_0[i] = g_xxxxzz_xxzzz_1[i] * fe_0 + g_xxxxzz_xxyzzz_1[i] * pa_y[i];

        g_xxxxyzz_xxzzzz_0[i] = g_xxxxzz_xxzzzz_1[i] * pa_y[i];

        g_xxxxyzz_xyyyyy_0[i] = 5.0 * g_xxxxzz_xyyyy_1[i] * fe_0 + g_xxxxzz_xyyyyy_1[i] * pa_y[i];

        g_xxxxyzz_xyyyyz_0[i] = 4.0 * g_xxxxzz_xyyyz_1[i] * fe_0 + g_xxxxzz_xyyyyz_1[i] * pa_y[i];

        g_xxxxyzz_xyyyzz_0[i] = 3.0 * g_xxxxzz_xyyzz_1[i] * fe_0 + g_xxxxzz_xyyyzz_1[i] * pa_y[i];

        g_xxxxyzz_xyyzzz_0[i] = 2.0 * g_xxxxzz_xyzzz_1[i] * fe_0 + g_xxxxzz_xyyzzz_1[i] * pa_y[i];

        g_xxxxyzz_xyzzzz_0[i] = g_xxxxzz_xzzzz_1[i] * fe_0 + g_xxxxzz_xyzzzz_1[i] * pa_y[i];

        g_xxxxyzz_xzzzzz_0[i] = g_xxxxzz_xzzzzz_1[i] * pa_y[i];

        g_xxxxyzz_yyyyyy_0[i] = 6.0 * g_xxxxzz_yyyyy_1[i] * fe_0 + g_xxxxzz_yyyyyy_1[i] * pa_y[i];

        g_xxxxyzz_yyyyyz_0[i] = 5.0 * g_xxxxzz_yyyyz_1[i] * fe_0 + g_xxxxzz_yyyyyz_1[i] * pa_y[i];

        g_xxxxyzz_yyyyzz_0[i] = 4.0 * g_xxxxzz_yyyzz_1[i] * fe_0 + g_xxxxzz_yyyyzz_1[i] * pa_y[i];

        g_xxxxyzz_yyyzzz_0[i] = 3.0 * g_xxxxzz_yyzzz_1[i] * fe_0 + g_xxxxzz_yyyzzz_1[i] * pa_y[i];

        g_xxxxyzz_yyzzzz_0[i] = 2.0 * g_xxxxzz_yzzzz_1[i] * fe_0 + g_xxxxzz_yyzzzz_1[i] * pa_y[i];

        g_xxxxyzz_yzzzzz_0[i] = g_xxxxzz_zzzzz_1[i] * fe_0 + g_xxxxzz_yzzzzz_1[i] * pa_y[i];

        g_xxxxyzz_zzzzzz_0[i] = g_xxxxzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : KI

    auto g_xxxxzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 252);

    auto g_xxxxzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 253);

    auto g_xxxxzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 254);

    auto g_xxxxzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 255);

    auto g_xxxxzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 256);

    auto g_xxxxzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 257);

    auto g_xxxxzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 258);

    auto g_xxxxzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 259);

    auto g_xxxxzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 260);

    auto g_xxxxzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 261);

    auto g_xxxxzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 262);

    auto g_xxxxzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 263);

    auto g_xxxxzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 264);

    auto g_xxxxzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 265);

    auto g_xxxxzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 266);

    auto g_xxxxzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 267);

    auto g_xxxxzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 268);

    auto g_xxxxzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 269);

    auto g_xxxxzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 270);

    auto g_xxxxzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 271);

    auto g_xxxxzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 272);

    auto g_xxxxzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 273);

    auto g_xxxxzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 274);

    auto g_xxxxzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 275);

    auto g_xxxxzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 276);

    auto g_xxxxzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 277);

    auto g_xxxxzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 278);

    auto g_xxxxzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 279);

    #pragma omp simd aligned(g_xxxxz_xxxxxx_0, g_xxxxz_xxxxxx_1, g_xxxxz_xxxxxy_0, g_xxxxz_xxxxxy_1, g_xxxxz_xxxxyy_0, g_xxxxz_xxxxyy_1, g_xxxxz_xxxyyy_0, g_xxxxz_xxxyyy_1, g_xxxxz_xxyyyy_0, g_xxxxz_xxyyyy_1, g_xxxxz_xyyyyy_0, g_xxxxz_xyyyyy_1, g_xxxxzz_xxxxxx_1, g_xxxxzz_xxxxxy_1, g_xxxxzz_xxxxyy_1, g_xxxxzz_xxxyyy_1, g_xxxxzz_xxyyyy_1, g_xxxxzz_xyyyyy_1, g_xxxxzzz_xxxxxx_0, g_xxxxzzz_xxxxxy_0, g_xxxxzzz_xxxxxz_0, g_xxxxzzz_xxxxyy_0, g_xxxxzzz_xxxxyz_0, g_xxxxzzz_xxxxzz_0, g_xxxxzzz_xxxyyy_0, g_xxxxzzz_xxxyyz_0, g_xxxxzzz_xxxyzz_0, g_xxxxzzz_xxxzzz_0, g_xxxxzzz_xxyyyy_0, g_xxxxzzz_xxyyyz_0, g_xxxxzzz_xxyyzz_0, g_xxxxzzz_xxyzzz_0, g_xxxxzzz_xxzzzz_0, g_xxxxzzz_xyyyyy_0, g_xxxxzzz_xyyyyz_0, g_xxxxzzz_xyyyzz_0, g_xxxxzzz_xyyzzz_0, g_xxxxzzz_xyzzzz_0, g_xxxxzzz_xzzzzz_0, g_xxxxzzz_yyyyyy_0, g_xxxxzzz_yyyyyz_0, g_xxxxzzz_yyyyzz_0, g_xxxxzzz_yyyzzz_0, g_xxxxzzz_yyzzzz_0, g_xxxxzzz_yzzzzz_0, g_xxxxzzz_zzzzzz_0, g_xxxzzz_xxxxxz_1, g_xxxzzz_xxxxyz_1, g_xxxzzz_xxxxz_1, g_xxxzzz_xxxxzz_1, g_xxxzzz_xxxyyz_1, g_xxxzzz_xxxyz_1, g_xxxzzz_xxxyzz_1, g_xxxzzz_xxxzz_1, g_xxxzzz_xxxzzz_1, g_xxxzzz_xxyyyz_1, g_xxxzzz_xxyyz_1, g_xxxzzz_xxyyzz_1, g_xxxzzz_xxyzz_1, g_xxxzzz_xxyzzz_1, g_xxxzzz_xxzzz_1, g_xxxzzz_xxzzzz_1, g_xxxzzz_xyyyyz_1, g_xxxzzz_xyyyz_1, g_xxxzzz_xyyyzz_1, g_xxxzzz_xyyzz_1, g_xxxzzz_xyyzzz_1, g_xxxzzz_xyzzz_1, g_xxxzzz_xyzzzz_1, g_xxxzzz_xzzzz_1, g_xxxzzz_xzzzzz_1, g_xxxzzz_yyyyyy_1, g_xxxzzz_yyyyyz_1, g_xxxzzz_yyyyz_1, g_xxxzzz_yyyyzz_1, g_xxxzzz_yyyzz_1, g_xxxzzz_yyyzzz_1, g_xxxzzz_yyzzz_1, g_xxxzzz_yyzzzz_1, g_xxxzzz_yzzzz_1, g_xxxzzz_yzzzzz_1, g_xxxzzz_zzzzz_1, g_xxxzzz_zzzzzz_1, g_xxzzz_xxxxxz_0, g_xxzzz_xxxxxz_1, g_xxzzz_xxxxyz_0, g_xxzzz_xxxxyz_1, g_xxzzz_xxxxzz_0, g_xxzzz_xxxxzz_1, g_xxzzz_xxxyyz_0, g_xxzzz_xxxyyz_1, g_xxzzz_xxxyzz_0, g_xxzzz_xxxyzz_1, g_xxzzz_xxxzzz_0, g_xxzzz_xxxzzz_1, g_xxzzz_xxyyyz_0, g_xxzzz_xxyyyz_1, g_xxzzz_xxyyzz_0, g_xxzzz_xxyyzz_1, g_xxzzz_xxyzzz_0, g_xxzzz_xxyzzz_1, g_xxzzz_xxzzzz_0, g_xxzzz_xxzzzz_1, g_xxzzz_xyyyyz_0, g_xxzzz_xyyyyz_1, g_xxzzz_xyyyzz_0, g_xxzzz_xyyyzz_1, g_xxzzz_xyyzzz_0, g_xxzzz_xyyzzz_1, g_xxzzz_xyzzzz_0, g_xxzzz_xyzzzz_1, g_xxzzz_xzzzzz_0, g_xxzzz_xzzzzz_1, g_xxzzz_yyyyyy_0, g_xxzzz_yyyyyy_1, g_xxzzz_yyyyyz_0, g_xxzzz_yyyyyz_1, g_xxzzz_yyyyzz_0, g_xxzzz_yyyyzz_1, g_xxzzz_yyyzzz_0, g_xxzzz_yyyzzz_1, g_xxzzz_yyzzzz_0, g_xxzzz_yyzzzz_1, g_xxzzz_yzzzzz_0, g_xxzzz_yzzzzz_1, g_xxzzz_zzzzzz_0, g_xxzzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzzz_xxxxxx_0[i] = 2.0 * g_xxxxz_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxxxx_1[i] * fz_be_0 + g_xxxxzz_xxxxxx_1[i] * pa_z[i];

        g_xxxxzzz_xxxxxy_0[i] = 2.0 * g_xxxxz_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxxxy_1[i] * fz_be_0 + g_xxxxzz_xxxxxy_1[i] * pa_z[i];

        g_xxxxzzz_xxxxxz_0[i] = 3.0 * g_xxzzz_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_xxxxz_1[i] * fe_0 + g_xxxzzz_xxxxxz_1[i] * pa_x[i];

        g_xxxxzzz_xxxxyy_0[i] = 2.0 * g_xxxxz_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxxyy_1[i] * fz_be_0 + g_xxxxzz_xxxxyy_1[i] * pa_z[i];

        g_xxxxzzz_xxxxyz_0[i] = 3.0 * g_xxzzz_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_xxxyz_1[i] * fe_0 + g_xxxzzz_xxxxyz_1[i] * pa_x[i];

        g_xxxxzzz_xxxxzz_0[i] = 3.0 * g_xxzzz_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_xxxzz_1[i] * fe_0 + g_xxxzzz_xxxxzz_1[i] * pa_x[i];

        g_xxxxzzz_xxxyyy_0[i] = 2.0 * g_xxxxz_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxyyy_1[i] * fz_be_0 + g_xxxxzz_xxxyyy_1[i] * pa_z[i];

        g_xxxxzzz_xxxyyz_0[i] = 3.0 * g_xxzzz_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_xxyyz_1[i] * fe_0 + g_xxxzzz_xxxyyz_1[i] * pa_x[i];

        g_xxxxzzz_xxxyzz_0[i] = 3.0 * g_xxzzz_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_xxyzz_1[i] * fe_0 + g_xxxzzz_xxxyzz_1[i] * pa_x[i];

        g_xxxxzzz_xxxzzz_0[i] = 3.0 * g_xxzzz_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_xxzzz_1[i] * fe_0 + g_xxxzzz_xxxzzz_1[i] * pa_x[i];

        g_xxxxzzz_xxyyyy_0[i] = 2.0 * g_xxxxz_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxyyyy_1[i] * fz_be_0 + g_xxxxzz_xxyyyy_1[i] * pa_z[i];

        g_xxxxzzz_xxyyyz_0[i] = 3.0 * g_xxzzz_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xyyyz_1[i] * fe_0 + g_xxxzzz_xxyyyz_1[i] * pa_x[i];

        g_xxxxzzz_xxyyzz_0[i] = 3.0 * g_xxzzz_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xyyzz_1[i] * fe_0 + g_xxxzzz_xxyyzz_1[i] * pa_x[i];

        g_xxxxzzz_xxyzzz_0[i] = 3.0 * g_xxzzz_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xyzzz_1[i] * fe_0 + g_xxxzzz_xxyzzz_1[i] * pa_x[i];

        g_xxxxzzz_xxzzzz_0[i] = 3.0 * g_xxzzz_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xzzzz_1[i] * fe_0 + g_xxxzzz_xxzzzz_1[i] * pa_x[i];

        g_xxxxzzz_xyyyyy_0[i] = 2.0 * g_xxxxz_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xyyyyy_1[i] * fz_be_0 + g_xxxxzz_xyyyyy_1[i] * pa_z[i];

        g_xxxxzzz_xyyyyz_0[i] = 3.0 * g_xxzzz_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyyyyz_1[i] * fz_be_0 + g_xxxzzz_yyyyz_1[i] * fe_0 + g_xxxzzz_xyyyyz_1[i] * pa_x[i];

        g_xxxxzzz_xyyyzz_0[i] = 3.0 * g_xxzzz_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyyyzz_1[i] * fz_be_0 + g_xxxzzz_yyyzz_1[i] * fe_0 + g_xxxzzz_xyyyzz_1[i] * pa_x[i];

        g_xxxxzzz_xyyzzz_0[i] = 3.0 * g_xxzzz_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyyzzz_1[i] * fz_be_0 + g_xxxzzz_yyzzz_1[i] * fe_0 + g_xxxzzz_xyyzzz_1[i] * pa_x[i];

        g_xxxxzzz_xyzzzz_0[i] = 3.0 * g_xxzzz_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyzzzz_1[i] * fz_be_0 + g_xxxzzz_yzzzz_1[i] * fe_0 + g_xxxzzz_xyzzzz_1[i] * pa_x[i];

        g_xxxxzzz_xzzzzz_0[i] = 3.0 * g_xxzzz_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xzzzzz_1[i] * fz_be_0 + g_xxxzzz_zzzzz_1[i] * fe_0 + g_xxxzzz_xzzzzz_1[i] * pa_x[i];

        g_xxxxzzz_yyyyyy_0[i] = 3.0 * g_xxzzz_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyyyy_1[i] * fz_be_0 + g_xxxzzz_yyyyyy_1[i] * pa_x[i];

        g_xxxxzzz_yyyyyz_0[i] = 3.0 * g_xxzzz_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyyyz_1[i] * fz_be_0 + g_xxxzzz_yyyyyz_1[i] * pa_x[i];

        g_xxxxzzz_yyyyzz_0[i] = 3.0 * g_xxzzz_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyyzz_1[i] * fz_be_0 + g_xxxzzz_yyyyzz_1[i] * pa_x[i];

        g_xxxxzzz_yyyzzz_0[i] = 3.0 * g_xxzzz_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyzzz_1[i] * fz_be_0 + g_xxxzzz_yyyzzz_1[i] * pa_x[i];

        g_xxxxzzz_yyzzzz_0[i] = 3.0 * g_xxzzz_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyzzzz_1[i] * fz_be_0 + g_xxxzzz_yyzzzz_1[i] * pa_x[i];

        g_xxxxzzz_yzzzzz_0[i] = 3.0 * g_xxzzz_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yzzzzz_1[i] * fz_be_0 + g_xxxzzz_yzzzzz_1[i] * pa_x[i];

        g_xxxxzzz_zzzzzz_0[i] = 3.0 * g_xxzzz_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_zzzzzz_1[i] * fz_be_0 + g_xxxzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : KI

    auto g_xxxyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 280);

    auto g_xxxyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 281);

    auto g_xxxyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 282);

    auto g_xxxyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 283);

    auto g_xxxyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 284);

    auto g_xxxyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 285);

    auto g_xxxyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 286);

    auto g_xxxyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 287);

    auto g_xxxyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 288);

    auto g_xxxyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 289);

    auto g_xxxyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 290);

    auto g_xxxyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 291);

    auto g_xxxyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 292);

    auto g_xxxyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 293);

    auto g_xxxyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 294);

    auto g_xxxyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 295);

    auto g_xxxyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 296);

    auto g_xxxyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 297);

    auto g_xxxyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 298);

    auto g_xxxyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 299);

    auto g_xxxyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 300);

    auto g_xxxyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 301);

    auto g_xxxyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 302);

    auto g_xxxyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 303);

    auto g_xxxyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 304);

    auto g_xxxyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 305);

    auto g_xxxyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 306);

    auto g_xxxyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 307);

    #pragma omp simd aligned(g_xxxyy_xxxxxx_0, g_xxxyy_xxxxxx_1, g_xxxyy_xxxxxz_0, g_xxxyy_xxxxxz_1, g_xxxyy_xxxxzz_0, g_xxxyy_xxxxzz_1, g_xxxyy_xxxzzz_0, g_xxxyy_xxxzzz_1, g_xxxyy_xxzzzz_0, g_xxxyy_xxzzzz_1, g_xxxyy_xzzzzz_0, g_xxxyy_xzzzzz_1, g_xxxyyy_xxxxxx_1, g_xxxyyy_xxxxxz_1, g_xxxyyy_xxxxzz_1, g_xxxyyy_xxxzzz_1, g_xxxyyy_xxzzzz_1, g_xxxyyy_xzzzzz_1, g_xxxyyyy_xxxxxx_0, g_xxxyyyy_xxxxxy_0, g_xxxyyyy_xxxxxz_0, g_xxxyyyy_xxxxyy_0, g_xxxyyyy_xxxxyz_0, g_xxxyyyy_xxxxzz_0, g_xxxyyyy_xxxyyy_0, g_xxxyyyy_xxxyyz_0, g_xxxyyyy_xxxyzz_0, g_xxxyyyy_xxxzzz_0, g_xxxyyyy_xxyyyy_0, g_xxxyyyy_xxyyyz_0, g_xxxyyyy_xxyyzz_0, g_xxxyyyy_xxyzzz_0, g_xxxyyyy_xxzzzz_0, g_xxxyyyy_xyyyyy_0, g_xxxyyyy_xyyyyz_0, g_xxxyyyy_xyyyzz_0, g_xxxyyyy_xyyzzz_0, g_xxxyyyy_xyzzzz_0, g_xxxyyyy_xzzzzz_0, g_xxxyyyy_yyyyyy_0, g_xxxyyyy_yyyyyz_0, g_xxxyyyy_yyyyzz_0, g_xxxyyyy_yyyzzz_0, g_xxxyyyy_yyzzzz_0, g_xxxyyyy_yzzzzz_0, g_xxxyyyy_zzzzzz_0, g_xxyyyy_xxxxxy_1, g_xxyyyy_xxxxy_1, g_xxyyyy_xxxxyy_1, g_xxyyyy_xxxxyz_1, g_xxyyyy_xxxyy_1, g_xxyyyy_xxxyyy_1, g_xxyyyy_xxxyyz_1, g_xxyyyy_xxxyz_1, g_xxyyyy_xxxyzz_1, g_xxyyyy_xxyyy_1, g_xxyyyy_xxyyyy_1, g_xxyyyy_xxyyyz_1, g_xxyyyy_xxyyz_1, g_xxyyyy_xxyyzz_1, g_xxyyyy_xxyzz_1, g_xxyyyy_xxyzzz_1, g_xxyyyy_xyyyy_1, g_xxyyyy_xyyyyy_1, g_xxyyyy_xyyyyz_1, g_xxyyyy_xyyyz_1, g_xxyyyy_xyyyzz_1, g_xxyyyy_xyyzz_1, g_xxyyyy_xyyzzz_1, g_xxyyyy_xyzzz_1, g_xxyyyy_xyzzzz_1, g_xxyyyy_yyyyy_1, g_xxyyyy_yyyyyy_1, g_xxyyyy_yyyyyz_1, g_xxyyyy_yyyyz_1, g_xxyyyy_yyyyzz_1, g_xxyyyy_yyyzz_1, g_xxyyyy_yyyzzz_1, g_xxyyyy_yyzzz_1, g_xxyyyy_yyzzzz_1, g_xxyyyy_yzzzz_1, g_xxyyyy_yzzzzz_1, g_xxyyyy_zzzzzz_1, g_xyyyy_xxxxxy_0, g_xyyyy_xxxxxy_1, g_xyyyy_xxxxyy_0, g_xyyyy_xxxxyy_1, g_xyyyy_xxxxyz_0, g_xyyyy_xxxxyz_1, g_xyyyy_xxxyyy_0, g_xyyyy_xxxyyy_1, g_xyyyy_xxxyyz_0, g_xyyyy_xxxyyz_1, g_xyyyy_xxxyzz_0, g_xyyyy_xxxyzz_1, g_xyyyy_xxyyyy_0, g_xyyyy_xxyyyy_1, g_xyyyy_xxyyyz_0, g_xyyyy_xxyyyz_1, g_xyyyy_xxyyzz_0, g_xyyyy_xxyyzz_1, g_xyyyy_xxyzzz_0, g_xyyyy_xxyzzz_1, g_xyyyy_xyyyyy_0, g_xyyyy_xyyyyy_1, g_xyyyy_xyyyyz_0, g_xyyyy_xyyyyz_1, g_xyyyy_xyyyzz_0, g_xyyyy_xyyyzz_1, g_xyyyy_xyyzzz_0, g_xyyyy_xyyzzz_1, g_xyyyy_xyzzzz_0, g_xyyyy_xyzzzz_1, g_xyyyy_yyyyyy_0, g_xyyyy_yyyyyy_1, g_xyyyy_yyyyyz_0, g_xyyyy_yyyyyz_1, g_xyyyy_yyyyzz_0, g_xyyyy_yyyyzz_1, g_xyyyy_yyyzzz_0, g_xyyyy_yyyzzz_1, g_xyyyy_yyzzzz_0, g_xyyyy_yyzzzz_1, g_xyyyy_yzzzzz_0, g_xyyyy_yzzzzz_1, g_xyyyy_zzzzzz_0, g_xyyyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyyy_xxxxxx_0[i] = 3.0 * g_xxxyy_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxxxx_1[i] * fz_be_0 + g_xxxyyy_xxxxxx_1[i] * pa_y[i];

        g_xxxyyyy_xxxxxy_0[i] = 2.0 * g_xyyyy_xxxxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxyyyy_xxxxy_1[i] * fe_0 + g_xxyyyy_xxxxxy_1[i] * pa_x[i];

        g_xxxyyyy_xxxxxz_0[i] = 3.0 * g_xxxyy_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxxxz_1[i] * fz_be_0 + g_xxxyyy_xxxxxz_1[i] * pa_y[i];

        g_xxxyyyy_xxxxyy_0[i] = 2.0 * g_xyyyy_xxxxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxyyyy_xxxyy_1[i] * fe_0 + g_xxyyyy_xxxxyy_1[i] * pa_x[i];

        g_xxxyyyy_xxxxyz_0[i] = 2.0 * g_xyyyy_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_xxxyz_1[i] * fe_0 + g_xxyyyy_xxxxyz_1[i] * pa_x[i];

        g_xxxyyyy_xxxxzz_0[i] = 3.0 * g_xxxyy_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxxzz_1[i] * fz_be_0 + g_xxxyyy_xxxxzz_1[i] * pa_y[i];

        g_xxxyyyy_xxxyyy_0[i] = 2.0 * g_xyyyy_xxxyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_xxyyy_1[i] * fe_0 + g_xxyyyy_xxxyyy_1[i] * pa_x[i];

        g_xxxyyyy_xxxyyz_0[i] = 2.0 * g_xyyyy_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_xxyyz_1[i] * fe_0 + g_xxyyyy_xxxyyz_1[i] * pa_x[i];

        g_xxxyyyy_xxxyzz_0[i] = 2.0 * g_xyyyy_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_xxyzz_1[i] * fe_0 + g_xxyyyy_xxxyzz_1[i] * pa_x[i];

        g_xxxyyyy_xxxzzz_0[i] = 3.0 * g_xxxyy_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxzzz_1[i] * fz_be_0 + g_xxxyyy_xxxzzz_1[i] * pa_y[i];

        g_xxxyyyy_xxyyyy_0[i] = 2.0 * g_xyyyy_xxyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyyyy_1[i] * fe_0 + g_xxyyyy_xxyyyy_1[i] * pa_x[i];

        g_xxxyyyy_xxyyyz_0[i] = 2.0 * g_xyyyy_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyyyz_1[i] * fe_0 + g_xxyyyy_xxyyyz_1[i] * pa_x[i];

        g_xxxyyyy_xxyyzz_0[i] = 2.0 * g_xyyyy_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyyzz_1[i] * fe_0 + g_xxyyyy_xxyyzz_1[i] * pa_x[i];

        g_xxxyyyy_xxyzzz_0[i] = 2.0 * g_xyyyy_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyzzz_1[i] * fe_0 + g_xxyyyy_xxyzzz_1[i] * pa_x[i];

        g_xxxyyyy_xxzzzz_0[i] = 3.0 * g_xxxyy_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxzzzz_1[i] * fz_be_0 + g_xxxyyy_xxzzzz_1[i] * pa_y[i];

        g_xxxyyyy_xyyyyy_0[i] = 2.0 * g_xyyyy_xyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyyyy_1[i] * fz_be_0 + g_xxyyyy_yyyyy_1[i] * fe_0 + g_xxyyyy_xyyyyy_1[i] * pa_x[i];

        g_xxxyyyy_xyyyyz_0[i] = 2.0 * g_xyyyy_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyyyz_1[i] * fz_be_0 + g_xxyyyy_yyyyz_1[i] * fe_0 + g_xxyyyy_xyyyyz_1[i] * pa_x[i];

        g_xxxyyyy_xyyyzz_0[i] = 2.0 * g_xyyyy_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyyzz_1[i] * fz_be_0 + g_xxyyyy_yyyzz_1[i] * fe_0 + g_xxyyyy_xyyyzz_1[i] * pa_x[i];

        g_xxxyyyy_xyyzzz_0[i] = 2.0 * g_xyyyy_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyzzz_1[i] * fz_be_0 + g_xxyyyy_yyzzz_1[i] * fe_0 + g_xxyyyy_xyyzzz_1[i] * pa_x[i];

        g_xxxyyyy_xyzzzz_0[i] = 2.0 * g_xyyyy_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyzzzz_1[i] * fz_be_0 + g_xxyyyy_yzzzz_1[i] * fe_0 + g_xxyyyy_xyzzzz_1[i] * pa_x[i];

        g_xxxyyyy_xzzzzz_0[i] = 3.0 * g_xxxyy_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xzzzzz_1[i] * fz_be_0 + g_xxxyyy_xzzzzz_1[i] * pa_y[i];

        g_xxxyyyy_yyyyyy_0[i] = 2.0 * g_xyyyy_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyyyy_1[i] * fz_be_0 + g_xxyyyy_yyyyyy_1[i] * pa_x[i];

        g_xxxyyyy_yyyyyz_0[i] = 2.0 * g_xyyyy_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyyyz_1[i] * fz_be_0 + g_xxyyyy_yyyyyz_1[i] * pa_x[i];

        g_xxxyyyy_yyyyzz_0[i] = 2.0 * g_xyyyy_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyyzz_1[i] * fz_be_0 + g_xxyyyy_yyyyzz_1[i] * pa_x[i];

        g_xxxyyyy_yyyzzz_0[i] = 2.0 * g_xyyyy_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyzzz_1[i] * fz_be_0 + g_xxyyyy_yyyzzz_1[i] * pa_x[i];

        g_xxxyyyy_yyzzzz_0[i] = 2.0 * g_xyyyy_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyzzzz_1[i] * fz_be_0 + g_xxyyyy_yyzzzz_1[i] * pa_x[i];

        g_xxxyyyy_yzzzzz_0[i] = 2.0 * g_xyyyy_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yzzzzz_1[i] * fz_be_0 + g_xxyyyy_yzzzzz_1[i] * pa_x[i];

        g_xxxyyyy_zzzzzz_0[i] = 2.0 * g_xyyyy_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_zzzzzz_1[i] * fz_be_0 + g_xxyyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : KI

    auto g_xxxyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 308);

    auto g_xxxyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 309);

    auto g_xxxyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 310);

    auto g_xxxyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 311);

    auto g_xxxyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 312);

    auto g_xxxyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 313);

    auto g_xxxyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 314);

    auto g_xxxyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 315);

    auto g_xxxyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 316);

    auto g_xxxyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 317);

    auto g_xxxyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 318);

    auto g_xxxyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 319);

    auto g_xxxyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 320);

    auto g_xxxyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 321);

    auto g_xxxyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 322);

    auto g_xxxyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 323);

    auto g_xxxyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 324);

    auto g_xxxyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 325);

    auto g_xxxyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 326);

    auto g_xxxyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 327);

    auto g_xxxyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 328);

    auto g_xxxyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 329);

    auto g_xxxyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 330);

    auto g_xxxyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 331);

    auto g_xxxyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 332);

    auto g_xxxyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 333);

    auto g_xxxyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 334);

    auto g_xxxyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 335);

    #pragma omp simd aligned(g_xxxyyy_xxxxx_1, g_xxxyyy_xxxxxx_1, g_xxxyyy_xxxxxy_1, g_xxxyyy_xxxxxz_1, g_xxxyyy_xxxxy_1, g_xxxyyy_xxxxyy_1, g_xxxyyy_xxxxyz_1, g_xxxyyy_xxxxz_1, g_xxxyyy_xxxxzz_1, g_xxxyyy_xxxyy_1, g_xxxyyy_xxxyyy_1, g_xxxyyy_xxxyyz_1, g_xxxyyy_xxxyz_1, g_xxxyyy_xxxyzz_1, g_xxxyyy_xxxzz_1, g_xxxyyy_xxxzzz_1, g_xxxyyy_xxyyy_1, g_xxxyyy_xxyyyy_1, g_xxxyyy_xxyyyz_1, g_xxxyyy_xxyyz_1, g_xxxyyy_xxyyzz_1, g_xxxyyy_xxyzz_1, g_xxxyyy_xxyzzz_1, g_xxxyyy_xxzzz_1, g_xxxyyy_xxzzzz_1, g_xxxyyy_xyyyy_1, g_xxxyyy_xyyyyy_1, g_xxxyyy_xyyyyz_1, g_xxxyyy_xyyyz_1, g_xxxyyy_xyyyzz_1, g_xxxyyy_xyyzz_1, g_xxxyyy_xyyzzz_1, g_xxxyyy_xyzzz_1, g_xxxyyy_xyzzzz_1, g_xxxyyy_xzzzz_1, g_xxxyyy_xzzzzz_1, g_xxxyyy_yyyyy_1, g_xxxyyy_yyyyyy_1, g_xxxyyy_yyyyyz_1, g_xxxyyy_yyyyz_1, g_xxxyyy_yyyyzz_1, g_xxxyyy_yyyzz_1, g_xxxyyy_yyyzzz_1, g_xxxyyy_yyzzz_1, g_xxxyyy_yyzzzz_1, g_xxxyyy_yzzzz_1, g_xxxyyy_yzzzzz_1, g_xxxyyy_zzzzz_1, g_xxxyyy_zzzzzz_1, g_xxxyyyz_xxxxxx_0, g_xxxyyyz_xxxxxy_0, g_xxxyyyz_xxxxxz_0, g_xxxyyyz_xxxxyy_0, g_xxxyyyz_xxxxyz_0, g_xxxyyyz_xxxxzz_0, g_xxxyyyz_xxxyyy_0, g_xxxyyyz_xxxyyz_0, g_xxxyyyz_xxxyzz_0, g_xxxyyyz_xxxzzz_0, g_xxxyyyz_xxyyyy_0, g_xxxyyyz_xxyyyz_0, g_xxxyyyz_xxyyzz_0, g_xxxyyyz_xxyzzz_0, g_xxxyyyz_xxzzzz_0, g_xxxyyyz_xyyyyy_0, g_xxxyyyz_xyyyyz_0, g_xxxyyyz_xyyyzz_0, g_xxxyyyz_xyyzzz_0, g_xxxyyyz_xyzzzz_0, g_xxxyyyz_xzzzzz_0, g_xxxyyyz_yyyyyy_0, g_xxxyyyz_yyyyyz_0, g_xxxyyyz_yyyyzz_0, g_xxxyyyz_yyyzzz_0, g_xxxyyyz_yyzzzz_0, g_xxxyyyz_yzzzzz_0, g_xxxyyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyyz_xxxxxx_0[i] = g_xxxyyy_xxxxxx_1[i] * pa_z[i];

        g_xxxyyyz_xxxxxy_0[i] = g_xxxyyy_xxxxxy_1[i] * pa_z[i];

        g_xxxyyyz_xxxxxz_0[i] = g_xxxyyy_xxxxx_1[i] * fe_0 + g_xxxyyy_xxxxxz_1[i] * pa_z[i];

        g_xxxyyyz_xxxxyy_0[i] = g_xxxyyy_xxxxyy_1[i] * pa_z[i];

        g_xxxyyyz_xxxxyz_0[i] = g_xxxyyy_xxxxy_1[i] * fe_0 + g_xxxyyy_xxxxyz_1[i] * pa_z[i];

        g_xxxyyyz_xxxxzz_0[i] = 2.0 * g_xxxyyy_xxxxz_1[i] * fe_0 + g_xxxyyy_xxxxzz_1[i] * pa_z[i];

        g_xxxyyyz_xxxyyy_0[i] = g_xxxyyy_xxxyyy_1[i] * pa_z[i];

        g_xxxyyyz_xxxyyz_0[i] = g_xxxyyy_xxxyy_1[i] * fe_0 + g_xxxyyy_xxxyyz_1[i] * pa_z[i];

        g_xxxyyyz_xxxyzz_0[i] = 2.0 * g_xxxyyy_xxxyz_1[i] * fe_0 + g_xxxyyy_xxxyzz_1[i] * pa_z[i];

        g_xxxyyyz_xxxzzz_0[i] = 3.0 * g_xxxyyy_xxxzz_1[i] * fe_0 + g_xxxyyy_xxxzzz_1[i] * pa_z[i];

        g_xxxyyyz_xxyyyy_0[i] = g_xxxyyy_xxyyyy_1[i] * pa_z[i];

        g_xxxyyyz_xxyyyz_0[i] = g_xxxyyy_xxyyy_1[i] * fe_0 + g_xxxyyy_xxyyyz_1[i] * pa_z[i];

        g_xxxyyyz_xxyyzz_0[i] = 2.0 * g_xxxyyy_xxyyz_1[i] * fe_0 + g_xxxyyy_xxyyzz_1[i] * pa_z[i];

        g_xxxyyyz_xxyzzz_0[i] = 3.0 * g_xxxyyy_xxyzz_1[i] * fe_0 + g_xxxyyy_xxyzzz_1[i] * pa_z[i];

        g_xxxyyyz_xxzzzz_0[i] = 4.0 * g_xxxyyy_xxzzz_1[i] * fe_0 + g_xxxyyy_xxzzzz_1[i] * pa_z[i];

        g_xxxyyyz_xyyyyy_0[i] = g_xxxyyy_xyyyyy_1[i] * pa_z[i];

        g_xxxyyyz_xyyyyz_0[i] = g_xxxyyy_xyyyy_1[i] * fe_0 + g_xxxyyy_xyyyyz_1[i] * pa_z[i];

        g_xxxyyyz_xyyyzz_0[i] = 2.0 * g_xxxyyy_xyyyz_1[i] * fe_0 + g_xxxyyy_xyyyzz_1[i] * pa_z[i];

        g_xxxyyyz_xyyzzz_0[i] = 3.0 * g_xxxyyy_xyyzz_1[i] * fe_0 + g_xxxyyy_xyyzzz_1[i] * pa_z[i];

        g_xxxyyyz_xyzzzz_0[i] = 4.0 * g_xxxyyy_xyzzz_1[i] * fe_0 + g_xxxyyy_xyzzzz_1[i] * pa_z[i];

        g_xxxyyyz_xzzzzz_0[i] = 5.0 * g_xxxyyy_xzzzz_1[i] * fe_0 + g_xxxyyy_xzzzzz_1[i] * pa_z[i];

        g_xxxyyyz_yyyyyy_0[i] = g_xxxyyy_yyyyyy_1[i] * pa_z[i];

        g_xxxyyyz_yyyyyz_0[i] = g_xxxyyy_yyyyy_1[i] * fe_0 + g_xxxyyy_yyyyyz_1[i] * pa_z[i];

        g_xxxyyyz_yyyyzz_0[i] = 2.0 * g_xxxyyy_yyyyz_1[i] * fe_0 + g_xxxyyy_yyyyzz_1[i] * pa_z[i];

        g_xxxyyyz_yyyzzz_0[i] = 3.0 * g_xxxyyy_yyyzz_1[i] * fe_0 + g_xxxyyy_yyyzzz_1[i] * pa_z[i];

        g_xxxyyyz_yyzzzz_0[i] = 4.0 * g_xxxyyy_yyzzz_1[i] * fe_0 + g_xxxyyy_yyzzzz_1[i] * pa_z[i];

        g_xxxyyyz_yzzzzz_0[i] = 5.0 * g_xxxyyy_yzzzz_1[i] * fe_0 + g_xxxyyy_yzzzzz_1[i] * pa_z[i];

        g_xxxyyyz_zzzzzz_0[i] = 6.0 * g_xxxyyy_zzzzz_1[i] * fe_0 + g_xxxyyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 336-364 components of targeted buffer : KI

    auto g_xxxyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 336);

    auto g_xxxyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 337);

    auto g_xxxyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 338);

    auto g_xxxyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 339);

    auto g_xxxyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 340);

    auto g_xxxyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 341);

    auto g_xxxyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 342);

    auto g_xxxyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 343);

    auto g_xxxyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 344);

    auto g_xxxyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 345);

    auto g_xxxyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 346);

    auto g_xxxyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 347);

    auto g_xxxyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 348);

    auto g_xxxyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 349);

    auto g_xxxyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 350);

    auto g_xxxyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 351);

    auto g_xxxyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 352);

    auto g_xxxyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 353);

    auto g_xxxyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 354);

    auto g_xxxyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 355);

    auto g_xxxyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 356);

    auto g_xxxyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 357);

    auto g_xxxyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 358);

    auto g_xxxyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 359);

    auto g_xxxyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 360);

    auto g_xxxyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 361);

    auto g_xxxyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 362);

    auto g_xxxyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 363);

    #pragma omp simd aligned(g_xxxyy_xxxxxy_0, g_xxxyy_xxxxxy_1, g_xxxyy_xxxxyy_0, g_xxxyy_xxxxyy_1, g_xxxyy_xxxyyy_0, g_xxxyy_xxxyyy_1, g_xxxyy_xxyyyy_0, g_xxxyy_xxyyyy_1, g_xxxyy_xyyyyy_0, g_xxxyy_xyyyyy_1, g_xxxyyz_xxxxxy_1, g_xxxyyz_xxxxyy_1, g_xxxyyz_xxxyyy_1, g_xxxyyz_xxyyyy_1, g_xxxyyz_xyyyyy_1, g_xxxyyzz_xxxxxx_0, g_xxxyyzz_xxxxxy_0, g_xxxyyzz_xxxxxz_0, g_xxxyyzz_xxxxyy_0, g_xxxyyzz_xxxxyz_0, g_xxxyyzz_xxxxzz_0, g_xxxyyzz_xxxyyy_0, g_xxxyyzz_xxxyyz_0, g_xxxyyzz_xxxyzz_0, g_xxxyyzz_xxxzzz_0, g_xxxyyzz_xxyyyy_0, g_xxxyyzz_xxyyyz_0, g_xxxyyzz_xxyyzz_0, g_xxxyyzz_xxyzzz_0, g_xxxyyzz_xxzzzz_0, g_xxxyyzz_xyyyyy_0, g_xxxyyzz_xyyyyz_0, g_xxxyyzz_xyyyzz_0, g_xxxyyzz_xyyzzz_0, g_xxxyyzz_xyzzzz_0, g_xxxyyzz_xzzzzz_0, g_xxxyyzz_yyyyyy_0, g_xxxyyzz_yyyyyz_0, g_xxxyyzz_yyyyzz_0, g_xxxyyzz_yyyzzz_0, g_xxxyyzz_yyzzzz_0, g_xxxyyzz_yzzzzz_0, g_xxxyyzz_zzzzzz_0, g_xxxyzz_xxxxxx_1, g_xxxyzz_xxxxxz_1, g_xxxyzz_xxxxzz_1, g_xxxyzz_xxxzzz_1, g_xxxyzz_xxzzzz_1, g_xxxyzz_xzzzzz_1, g_xxxzz_xxxxxx_0, g_xxxzz_xxxxxx_1, g_xxxzz_xxxxxz_0, g_xxxzz_xxxxxz_1, g_xxxzz_xxxxzz_0, g_xxxzz_xxxxzz_1, g_xxxzz_xxxzzz_0, g_xxxzz_xxxzzz_1, g_xxxzz_xxzzzz_0, g_xxxzz_xxzzzz_1, g_xxxzz_xzzzzz_0, g_xxxzz_xzzzzz_1, g_xxyyzz_xxxxyz_1, g_xxyyzz_xxxyyz_1, g_xxyyzz_xxxyz_1, g_xxyyzz_xxxyzz_1, g_xxyyzz_xxyyyz_1, g_xxyyzz_xxyyz_1, g_xxyyzz_xxyyzz_1, g_xxyyzz_xxyzz_1, g_xxyyzz_xxyzzz_1, g_xxyyzz_xyyyyz_1, g_xxyyzz_xyyyz_1, g_xxyyzz_xyyyzz_1, g_xxyyzz_xyyzz_1, g_xxyyzz_xyyzzz_1, g_xxyyzz_xyzzz_1, g_xxyyzz_xyzzzz_1, g_xxyyzz_yyyyyy_1, g_xxyyzz_yyyyyz_1, g_xxyyzz_yyyyz_1, g_xxyyzz_yyyyzz_1, g_xxyyzz_yyyzz_1, g_xxyyzz_yyyzzz_1, g_xxyyzz_yyzzz_1, g_xxyyzz_yyzzzz_1, g_xxyyzz_yzzzz_1, g_xxyyzz_yzzzzz_1, g_xxyyzz_zzzzzz_1, g_xyyzz_xxxxyz_0, g_xyyzz_xxxxyz_1, g_xyyzz_xxxyyz_0, g_xyyzz_xxxyyz_1, g_xyyzz_xxxyzz_0, g_xyyzz_xxxyzz_1, g_xyyzz_xxyyyz_0, g_xyyzz_xxyyyz_1, g_xyyzz_xxyyzz_0, g_xyyzz_xxyyzz_1, g_xyyzz_xxyzzz_0, g_xyyzz_xxyzzz_1, g_xyyzz_xyyyyz_0, g_xyyzz_xyyyyz_1, g_xyyzz_xyyyzz_0, g_xyyzz_xyyyzz_1, g_xyyzz_xyyzzz_0, g_xyyzz_xyyzzz_1, g_xyyzz_xyzzzz_0, g_xyyzz_xyzzzz_1, g_xyyzz_yyyyyy_0, g_xyyzz_yyyyyy_1, g_xyyzz_yyyyyz_0, g_xyyzz_yyyyyz_1, g_xyyzz_yyyyzz_0, g_xyyzz_yyyyzz_1, g_xyyzz_yyyzzz_0, g_xyyzz_yyyzzz_1, g_xyyzz_yyzzzz_0, g_xyyzz_yyzzzz_1, g_xyyzz_yzzzzz_0, g_xyyzz_yzzzzz_1, g_xyyzz_zzzzzz_0, g_xyyzz_zzzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyzz_xxxxxx_0[i] = g_xxxzz_xxxxxx_0[i] * fbe_0 - g_xxxzz_xxxxxx_1[i] * fz_be_0 + g_xxxyzz_xxxxxx_1[i] * pa_y[i];

        g_xxxyyzz_xxxxxy_0[i] = g_xxxyy_xxxxxy_0[i] * fbe_0 - g_xxxyy_xxxxxy_1[i] * fz_be_0 + g_xxxyyz_xxxxxy_1[i] * pa_z[i];

        g_xxxyyzz_xxxxxz_0[i] = g_xxxzz_xxxxxz_0[i] * fbe_0 - g_xxxzz_xxxxxz_1[i] * fz_be_0 + g_xxxyzz_xxxxxz_1[i] * pa_y[i];

        g_xxxyyzz_xxxxyy_0[i] = g_xxxyy_xxxxyy_0[i] * fbe_0 - g_xxxyy_xxxxyy_1[i] * fz_be_0 + g_xxxyyz_xxxxyy_1[i] * pa_z[i];

        g_xxxyyzz_xxxxyz_0[i] = 2.0 * g_xyyzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_xxxyz_1[i] * fe_0 + g_xxyyzz_xxxxyz_1[i] * pa_x[i];

        g_xxxyyzz_xxxxzz_0[i] = g_xxxzz_xxxxzz_0[i] * fbe_0 - g_xxxzz_xxxxzz_1[i] * fz_be_0 + g_xxxyzz_xxxxzz_1[i] * pa_y[i];

        g_xxxyyzz_xxxyyy_0[i] = g_xxxyy_xxxyyy_0[i] * fbe_0 - g_xxxyy_xxxyyy_1[i] * fz_be_0 + g_xxxyyz_xxxyyy_1[i] * pa_z[i];

        g_xxxyyzz_xxxyyz_0[i] = 2.0 * g_xyyzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_xxyyz_1[i] * fe_0 + g_xxyyzz_xxxyyz_1[i] * pa_x[i];

        g_xxxyyzz_xxxyzz_0[i] = 2.0 * g_xyyzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_xxyzz_1[i] * fe_0 + g_xxyyzz_xxxyzz_1[i] * pa_x[i];

        g_xxxyyzz_xxxzzz_0[i] = g_xxxzz_xxxzzz_0[i] * fbe_0 - g_xxxzz_xxxzzz_1[i] * fz_be_0 + g_xxxyzz_xxxzzz_1[i] * pa_y[i];

        g_xxxyyzz_xxyyyy_0[i] = g_xxxyy_xxyyyy_0[i] * fbe_0 - g_xxxyy_xxyyyy_1[i] * fz_be_0 + g_xxxyyz_xxyyyy_1[i] * pa_z[i];

        g_xxxyyzz_xxyyyz_0[i] = 2.0 * g_xyyzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_xyyyz_1[i] * fe_0 + g_xxyyzz_xxyyyz_1[i] * pa_x[i];

        g_xxxyyzz_xxyyzz_0[i] = 2.0 * g_xyyzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_xyyzz_1[i] * fe_0 + g_xxyyzz_xxyyzz_1[i] * pa_x[i];

        g_xxxyyzz_xxyzzz_0[i] = 2.0 * g_xyyzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_xyzzz_1[i] * fe_0 + g_xxyyzz_xxyzzz_1[i] * pa_x[i];

        g_xxxyyzz_xxzzzz_0[i] = g_xxxzz_xxzzzz_0[i] * fbe_0 - g_xxxzz_xxzzzz_1[i] * fz_be_0 + g_xxxyzz_xxzzzz_1[i] * pa_y[i];

        g_xxxyyzz_xyyyyy_0[i] = g_xxxyy_xyyyyy_0[i] * fbe_0 - g_xxxyy_xyyyyy_1[i] * fz_be_0 + g_xxxyyz_xyyyyy_1[i] * pa_z[i];

        g_xxxyyzz_xyyyyz_0[i] = 2.0 * g_xyyzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyyyyz_1[i] * fz_be_0 + g_xxyyzz_yyyyz_1[i] * fe_0 + g_xxyyzz_xyyyyz_1[i] * pa_x[i];

        g_xxxyyzz_xyyyzz_0[i] = 2.0 * g_xyyzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyyyzz_1[i] * fz_be_0 + g_xxyyzz_yyyzz_1[i] * fe_0 + g_xxyyzz_xyyyzz_1[i] * pa_x[i];

        g_xxxyyzz_xyyzzz_0[i] = 2.0 * g_xyyzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyyzzz_1[i] * fz_be_0 + g_xxyyzz_yyzzz_1[i] * fe_0 + g_xxyyzz_xyyzzz_1[i] * pa_x[i];

        g_xxxyyzz_xyzzzz_0[i] = 2.0 * g_xyyzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyzzzz_1[i] * fz_be_0 + g_xxyyzz_yzzzz_1[i] * fe_0 + g_xxyyzz_xyzzzz_1[i] * pa_x[i];

        g_xxxyyzz_xzzzzz_0[i] = g_xxxzz_xzzzzz_0[i] * fbe_0 - g_xxxzz_xzzzzz_1[i] * fz_be_0 + g_xxxyzz_xzzzzz_1[i] * pa_y[i];

        g_xxxyyzz_yyyyyy_0[i] = 2.0 * g_xyyzz_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyyyy_1[i] * fz_be_0 + g_xxyyzz_yyyyyy_1[i] * pa_x[i];

        g_xxxyyzz_yyyyyz_0[i] = 2.0 * g_xyyzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyyyz_1[i] * fz_be_0 + g_xxyyzz_yyyyyz_1[i] * pa_x[i];

        g_xxxyyzz_yyyyzz_0[i] = 2.0 * g_xyyzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyyzz_1[i] * fz_be_0 + g_xxyyzz_yyyyzz_1[i] * pa_x[i];

        g_xxxyyzz_yyyzzz_0[i] = 2.0 * g_xyyzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyzzz_1[i] * fz_be_0 + g_xxyyzz_yyyzzz_1[i] * pa_x[i];

        g_xxxyyzz_yyzzzz_0[i] = 2.0 * g_xyyzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyzzzz_1[i] * fz_be_0 + g_xxyyzz_yyzzzz_1[i] * pa_x[i];

        g_xxxyyzz_yzzzzz_0[i] = 2.0 * g_xyyzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yzzzzz_1[i] * fz_be_0 + g_xxyyzz_yzzzzz_1[i] * pa_x[i];

        g_xxxyyzz_zzzzzz_0[i] = 2.0 * g_xyyzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_zzzzzz_1[i] * fz_be_0 + g_xxyyzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : KI

    auto g_xxxyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 364);

    auto g_xxxyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 365);

    auto g_xxxyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 366);

    auto g_xxxyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 367);

    auto g_xxxyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 368);

    auto g_xxxyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 369);

    auto g_xxxyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 370);

    auto g_xxxyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 371);

    auto g_xxxyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 372);

    auto g_xxxyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 373);

    auto g_xxxyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 374);

    auto g_xxxyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 375);

    auto g_xxxyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 376);

    auto g_xxxyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 377);

    auto g_xxxyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 378);

    auto g_xxxyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 379);

    auto g_xxxyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 380);

    auto g_xxxyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 381);

    auto g_xxxyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 382);

    auto g_xxxyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 383);

    auto g_xxxyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 384);

    auto g_xxxyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 385);

    auto g_xxxyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 386);

    auto g_xxxyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 387);

    auto g_xxxyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 388);

    auto g_xxxyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 389);

    auto g_xxxyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 390);

    auto g_xxxyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 391);

    #pragma omp simd aligned(g_xxxyzzz_xxxxxx_0, g_xxxyzzz_xxxxxy_0, g_xxxyzzz_xxxxxz_0, g_xxxyzzz_xxxxyy_0, g_xxxyzzz_xxxxyz_0, g_xxxyzzz_xxxxzz_0, g_xxxyzzz_xxxyyy_0, g_xxxyzzz_xxxyyz_0, g_xxxyzzz_xxxyzz_0, g_xxxyzzz_xxxzzz_0, g_xxxyzzz_xxyyyy_0, g_xxxyzzz_xxyyyz_0, g_xxxyzzz_xxyyzz_0, g_xxxyzzz_xxyzzz_0, g_xxxyzzz_xxzzzz_0, g_xxxyzzz_xyyyyy_0, g_xxxyzzz_xyyyyz_0, g_xxxyzzz_xyyyzz_0, g_xxxyzzz_xyyzzz_0, g_xxxyzzz_xyzzzz_0, g_xxxyzzz_xzzzzz_0, g_xxxyzzz_yyyyyy_0, g_xxxyzzz_yyyyyz_0, g_xxxyzzz_yyyyzz_0, g_xxxyzzz_yyyzzz_0, g_xxxyzzz_yyzzzz_0, g_xxxyzzz_yzzzzz_0, g_xxxyzzz_zzzzzz_0, g_xxxzzz_xxxxx_1, g_xxxzzz_xxxxxx_1, g_xxxzzz_xxxxxy_1, g_xxxzzz_xxxxxz_1, g_xxxzzz_xxxxy_1, g_xxxzzz_xxxxyy_1, g_xxxzzz_xxxxyz_1, g_xxxzzz_xxxxz_1, g_xxxzzz_xxxxzz_1, g_xxxzzz_xxxyy_1, g_xxxzzz_xxxyyy_1, g_xxxzzz_xxxyyz_1, g_xxxzzz_xxxyz_1, g_xxxzzz_xxxyzz_1, g_xxxzzz_xxxzz_1, g_xxxzzz_xxxzzz_1, g_xxxzzz_xxyyy_1, g_xxxzzz_xxyyyy_1, g_xxxzzz_xxyyyz_1, g_xxxzzz_xxyyz_1, g_xxxzzz_xxyyzz_1, g_xxxzzz_xxyzz_1, g_xxxzzz_xxyzzz_1, g_xxxzzz_xxzzz_1, g_xxxzzz_xxzzzz_1, g_xxxzzz_xyyyy_1, g_xxxzzz_xyyyyy_1, g_xxxzzz_xyyyyz_1, g_xxxzzz_xyyyz_1, g_xxxzzz_xyyyzz_1, g_xxxzzz_xyyzz_1, g_xxxzzz_xyyzzz_1, g_xxxzzz_xyzzz_1, g_xxxzzz_xyzzzz_1, g_xxxzzz_xzzzz_1, g_xxxzzz_xzzzzz_1, g_xxxzzz_yyyyy_1, g_xxxzzz_yyyyyy_1, g_xxxzzz_yyyyyz_1, g_xxxzzz_yyyyz_1, g_xxxzzz_yyyyzz_1, g_xxxzzz_yyyzz_1, g_xxxzzz_yyyzzz_1, g_xxxzzz_yyzzz_1, g_xxxzzz_yyzzzz_1, g_xxxzzz_yzzzz_1, g_xxxzzz_yzzzzz_1, g_xxxzzz_zzzzz_1, g_xxxzzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzzz_xxxxxx_0[i] = g_xxxzzz_xxxxxx_1[i] * pa_y[i];

        g_xxxyzzz_xxxxxy_0[i] = g_xxxzzz_xxxxx_1[i] * fe_0 + g_xxxzzz_xxxxxy_1[i] * pa_y[i];

        g_xxxyzzz_xxxxxz_0[i] = g_xxxzzz_xxxxxz_1[i] * pa_y[i];

        g_xxxyzzz_xxxxyy_0[i] = 2.0 * g_xxxzzz_xxxxy_1[i] * fe_0 + g_xxxzzz_xxxxyy_1[i] * pa_y[i];

        g_xxxyzzz_xxxxyz_0[i] = g_xxxzzz_xxxxz_1[i] * fe_0 + g_xxxzzz_xxxxyz_1[i] * pa_y[i];

        g_xxxyzzz_xxxxzz_0[i] = g_xxxzzz_xxxxzz_1[i] * pa_y[i];

        g_xxxyzzz_xxxyyy_0[i] = 3.0 * g_xxxzzz_xxxyy_1[i] * fe_0 + g_xxxzzz_xxxyyy_1[i] * pa_y[i];

        g_xxxyzzz_xxxyyz_0[i] = 2.0 * g_xxxzzz_xxxyz_1[i] * fe_0 + g_xxxzzz_xxxyyz_1[i] * pa_y[i];

        g_xxxyzzz_xxxyzz_0[i] = g_xxxzzz_xxxzz_1[i] * fe_0 + g_xxxzzz_xxxyzz_1[i] * pa_y[i];

        g_xxxyzzz_xxxzzz_0[i] = g_xxxzzz_xxxzzz_1[i] * pa_y[i];

        g_xxxyzzz_xxyyyy_0[i] = 4.0 * g_xxxzzz_xxyyy_1[i] * fe_0 + g_xxxzzz_xxyyyy_1[i] * pa_y[i];

        g_xxxyzzz_xxyyyz_0[i] = 3.0 * g_xxxzzz_xxyyz_1[i] * fe_0 + g_xxxzzz_xxyyyz_1[i] * pa_y[i];

        g_xxxyzzz_xxyyzz_0[i] = 2.0 * g_xxxzzz_xxyzz_1[i] * fe_0 + g_xxxzzz_xxyyzz_1[i] * pa_y[i];

        g_xxxyzzz_xxyzzz_0[i] = g_xxxzzz_xxzzz_1[i] * fe_0 + g_xxxzzz_xxyzzz_1[i] * pa_y[i];

        g_xxxyzzz_xxzzzz_0[i] = g_xxxzzz_xxzzzz_1[i] * pa_y[i];

        g_xxxyzzz_xyyyyy_0[i] = 5.0 * g_xxxzzz_xyyyy_1[i] * fe_0 + g_xxxzzz_xyyyyy_1[i] * pa_y[i];

        g_xxxyzzz_xyyyyz_0[i] = 4.0 * g_xxxzzz_xyyyz_1[i] * fe_0 + g_xxxzzz_xyyyyz_1[i] * pa_y[i];

        g_xxxyzzz_xyyyzz_0[i] = 3.0 * g_xxxzzz_xyyzz_1[i] * fe_0 + g_xxxzzz_xyyyzz_1[i] * pa_y[i];

        g_xxxyzzz_xyyzzz_0[i] = 2.0 * g_xxxzzz_xyzzz_1[i] * fe_0 + g_xxxzzz_xyyzzz_1[i] * pa_y[i];

        g_xxxyzzz_xyzzzz_0[i] = g_xxxzzz_xzzzz_1[i] * fe_0 + g_xxxzzz_xyzzzz_1[i] * pa_y[i];

        g_xxxyzzz_xzzzzz_0[i] = g_xxxzzz_xzzzzz_1[i] * pa_y[i];

        g_xxxyzzz_yyyyyy_0[i] = 6.0 * g_xxxzzz_yyyyy_1[i] * fe_0 + g_xxxzzz_yyyyyy_1[i] * pa_y[i];

        g_xxxyzzz_yyyyyz_0[i] = 5.0 * g_xxxzzz_yyyyz_1[i] * fe_0 + g_xxxzzz_yyyyyz_1[i] * pa_y[i];

        g_xxxyzzz_yyyyzz_0[i] = 4.0 * g_xxxzzz_yyyzz_1[i] * fe_0 + g_xxxzzz_yyyyzz_1[i] * pa_y[i];

        g_xxxyzzz_yyyzzz_0[i] = 3.0 * g_xxxzzz_yyzzz_1[i] * fe_0 + g_xxxzzz_yyyzzz_1[i] * pa_y[i];

        g_xxxyzzz_yyzzzz_0[i] = 2.0 * g_xxxzzz_yzzzz_1[i] * fe_0 + g_xxxzzz_yyzzzz_1[i] * pa_y[i];

        g_xxxyzzz_yzzzzz_0[i] = g_xxxzzz_zzzzz_1[i] * fe_0 + g_xxxzzz_yzzzzz_1[i] * pa_y[i];

        g_xxxyzzz_zzzzzz_0[i] = g_xxxzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : KI

    auto g_xxxzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 392);

    auto g_xxxzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 393);

    auto g_xxxzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 394);

    auto g_xxxzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 395);

    auto g_xxxzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 396);

    auto g_xxxzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 397);

    auto g_xxxzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 398);

    auto g_xxxzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 399);

    auto g_xxxzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 400);

    auto g_xxxzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 401);

    auto g_xxxzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 402);

    auto g_xxxzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 403);

    auto g_xxxzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 404);

    auto g_xxxzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 405);

    auto g_xxxzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 406);

    auto g_xxxzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 407);

    auto g_xxxzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 408);

    auto g_xxxzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 409);

    auto g_xxxzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 410);

    auto g_xxxzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 411);

    auto g_xxxzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 412);

    auto g_xxxzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 413);

    auto g_xxxzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 414);

    auto g_xxxzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 415);

    auto g_xxxzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 416);

    auto g_xxxzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 417);

    auto g_xxxzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 418);

    auto g_xxxzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 419);

    #pragma omp simd aligned(g_xxxzz_xxxxxx_0, g_xxxzz_xxxxxx_1, g_xxxzz_xxxxxy_0, g_xxxzz_xxxxxy_1, g_xxxzz_xxxxyy_0, g_xxxzz_xxxxyy_1, g_xxxzz_xxxyyy_0, g_xxxzz_xxxyyy_1, g_xxxzz_xxyyyy_0, g_xxxzz_xxyyyy_1, g_xxxzz_xyyyyy_0, g_xxxzz_xyyyyy_1, g_xxxzzz_xxxxxx_1, g_xxxzzz_xxxxxy_1, g_xxxzzz_xxxxyy_1, g_xxxzzz_xxxyyy_1, g_xxxzzz_xxyyyy_1, g_xxxzzz_xyyyyy_1, g_xxxzzzz_xxxxxx_0, g_xxxzzzz_xxxxxy_0, g_xxxzzzz_xxxxxz_0, g_xxxzzzz_xxxxyy_0, g_xxxzzzz_xxxxyz_0, g_xxxzzzz_xxxxzz_0, g_xxxzzzz_xxxyyy_0, g_xxxzzzz_xxxyyz_0, g_xxxzzzz_xxxyzz_0, g_xxxzzzz_xxxzzz_0, g_xxxzzzz_xxyyyy_0, g_xxxzzzz_xxyyyz_0, g_xxxzzzz_xxyyzz_0, g_xxxzzzz_xxyzzz_0, g_xxxzzzz_xxzzzz_0, g_xxxzzzz_xyyyyy_0, g_xxxzzzz_xyyyyz_0, g_xxxzzzz_xyyyzz_0, g_xxxzzzz_xyyzzz_0, g_xxxzzzz_xyzzzz_0, g_xxxzzzz_xzzzzz_0, g_xxxzzzz_yyyyyy_0, g_xxxzzzz_yyyyyz_0, g_xxxzzzz_yyyyzz_0, g_xxxzzzz_yyyzzz_0, g_xxxzzzz_yyzzzz_0, g_xxxzzzz_yzzzzz_0, g_xxxzzzz_zzzzzz_0, g_xxzzzz_xxxxxz_1, g_xxzzzz_xxxxyz_1, g_xxzzzz_xxxxz_1, g_xxzzzz_xxxxzz_1, g_xxzzzz_xxxyyz_1, g_xxzzzz_xxxyz_1, g_xxzzzz_xxxyzz_1, g_xxzzzz_xxxzz_1, g_xxzzzz_xxxzzz_1, g_xxzzzz_xxyyyz_1, g_xxzzzz_xxyyz_1, g_xxzzzz_xxyyzz_1, g_xxzzzz_xxyzz_1, g_xxzzzz_xxyzzz_1, g_xxzzzz_xxzzz_1, g_xxzzzz_xxzzzz_1, g_xxzzzz_xyyyyz_1, g_xxzzzz_xyyyz_1, g_xxzzzz_xyyyzz_1, g_xxzzzz_xyyzz_1, g_xxzzzz_xyyzzz_1, g_xxzzzz_xyzzz_1, g_xxzzzz_xyzzzz_1, g_xxzzzz_xzzzz_1, g_xxzzzz_xzzzzz_1, g_xxzzzz_yyyyyy_1, g_xxzzzz_yyyyyz_1, g_xxzzzz_yyyyz_1, g_xxzzzz_yyyyzz_1, g_xxzzzz_yyyzz_1, g_xxzzzz_yyyzzz_1, g_xxzzzz_yyzzz_1, g_xxzzzz_yyzzzz_1, g_xxzzzz_yzzzz_1, g_xxzzzz_yzzzzz_1, g_xxzzzz_zzzzz_1, g_xxzzzz_zzzzzz_1, g_xzzzz_xxxxxz_0, g_xzzzz_xxxxxz_1, g_xzzzz_xxxxyz_0, g_xzzzz_xxxxyz_1, g_xzzzz_xxxxzz_0, g_xzzzz_xxxxzz_1, g_xzzzz_xxxyyz_0, g_xzzzz_xxxyyz_1, g_xzzzz_xxxyzz_0, g_xzzzz_xxxyzz_1, g_xzzzz_xxxzzz_0, g_xzzzz_xxxzzz_1, g_xzzzz_xxyyyz_0, g_xzzzz_xxyyyz_1, g_xzzzz_xxyyzz_0, g_xzzzz_xxyyzz_1, g_xzzzz_xxyzzz_0, g_xzzzz_xxyzzz_1, g_xzzzz_xxzzzz_0, g_xzzzz_xxzzzz_1, g_xzzzz_xyyyyz_0, g_xzzzz_xyyyyz_1, g_xzzzz_xyyyzz_0, g_xzzzz_xyyyzz_1, g_xzzzz_xyyzzz_0, g_xzzzz_xyyzzz_1, g_xzzzz_xyzzzz_0, g_xzzzz_xyzzzz_1, g_xzzzz_xzzzzz_0, g_xzzzz_xzzzzz_1, g_xzzzz_yyyyyy_0, g_xzzzz_yyyyyy_1, g_xzzzz_yyyyyz_0, g_xzzzz_yyyyyz_1, g_xzzzz_yyyyzz_0, g_xzzzz_yyyyzz_1, g_xzzzz_yyyzzz_0, g_xzzzz_yyyzzz_1, g_xzzzz_yyzzzz_0, g_xzzzz_yyzzzz_1, g_xzzzz_yzzzzz_0, g_xzzzz_yzzzzz_1, g_xzzzz_zzzzzz_0, g_xzzzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzzz_xxxxxx_0[i] = 3.0 * g_xxxzz_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxxxx_1[i] * fz_be_0 + g_xxxzzz_xxxxxx_1[i] * pa_z[i];

        g_xxxzzzz_xxxxxy_0[i] = 3.0 * g_xxxzz_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxxxy_1[i] * fz_be_0 + g_xxxzzz_xxxxxy_1[i] * pa_z[i];

        g_xxxzzzz_xxxxxz_0[i] = 2.0 * g_xzzzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_xxxxz_1[i] * fe_0 + g_xxzzzz_xxxxxz_1[i] * pa_x[i];

        g_xxxzzzz_xxxxyy_0[i] = 3.0 * g_xxxzz_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxxyy_1[i] * fz_be_0 + g_xxxzzz_xxxxyy_1[i] * pa_z[i];

        g_xxxzzzz_xxxxyz_0[i] = 2.0 * g_xzzzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_xxxyz_1[i] * fe_0 + g_xxzzzz_xxxxyz_1[i] * pa_x[i];

        g_xxxzzzz_xxxxzz_0[i] = 2.0 * g_xzzzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_xxxzz_1[i] * fe_0 + g_xxzzzz_xxxxzz_1[i] * pa_x[i];

        g_xxxzzzz_xxxyyy_0[i] = 3.0 * g_xxxzz_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxyyy_1[i] * fz_be_0 + g_xxxzzz_xxxyyy_1[i] * pa_z[i];

        g_xxxzzzz_xxxyyz_0[i] = 2.0 * g_xzzzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_xxyyz_1[i] * fe_0 + g_xxzzzz_xxxyyz_1[i] * pa_x[i];

        g_xxxzzzz_xxxyzz_0[i] = 2.0 * g_xzzzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_xxyzz_1[i] * fe_0 + g_xxzzzz_xxxyzz_1[i] * pa_x[i];

        g_xxxzzzz_xxxzzz_0[i] = 2.0 * g_xzzzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_xxzzz_1[i] * fe_0 + g_xxzzzz_xxxzzz_1[i] * pa_x[i];

        g_xxxzzzz_xxyyyy_0[i] = 3.0 * g_xxxzz_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxyyyy_1[i] * fz_be_0 + g_xxxzzz_xxyyyy_1[i] * pa_z[i];

        g_xxxzzzz_xxyyyz_0[i] = 2.0 * g_xzzzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xyyyz_1[i] * fe_0 + g_xxzzzz_xxyyyz_1[i] * pa_x[i];

        g_xxxzzzz_xxyyzz_0[i] = 2.0 * g_xzzzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xyyzz_1[i] * fe_0 + g_xxzzzz_xxyyzz_1[i] * pa_x[i];

        g_xxxzzzz_xxyzzz_0[i] = 2.0 * g_xzzzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xyzzz_1[i] * fe_0 + g_xxzzzz_xxyzzz_1[i] * pa_x[i];

        g_xxxzzzz_xxzzzz_0[i] = 2.0 * g_xzzzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xzzzz_1[i] * fe_0 + g_xxzzzz_xxzzzz_1[i] * pa_x[i];

        g_xxxzzzz_xyyyyy_0[i] = 3.0 * g_xxxzz_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xyyyyy_1[i] * fz_be_0 + g_xxxzzz_xyyyyy_1[i] * pa_z[i];

        g_xxxzzzz_xyyyyz_0[i] = 2.0 * g_xzzzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyyyyz_1[i] * fz_be_0 + g_xxzzzz_yyyyz_1[i] * fe_0 + g_xxzzzz_xyyyyz_1[i] * pa_x[i];

        g_xxxzzzz_xyyyzz_0[i] = 2.0 * g_xzzzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyyyzz_1[i] * fz_be_0 + g_xxzzzz_yyyzz_1[i] * fe_0 + g_xxzzzz_xyyyzz_1[i] * pa_x[i];

        g_xxxzzzz_xyyzzz_0[i] = 2.0 * g_xzzzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyyzzz_1[i] * fz_be_0 + g_xxzzzz_yyzzz_1[i] * fe_0 + g_xxzzzz_xyyzzz_1[i] * pa_x[i];

        g_xxxzzzz_xyzzzz_0[i] = 2.0 * g_xzzzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyzzzz_1[i] * fz_be_0 + g_xxzzzz_yzzzz_1[i] * fe_0 + g_xxzzzz_xyzzzz_1[i] * pa_x[i];

        g_xxxzzzz_xzzzzz_0[i] = 2.0 * g_xzzzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xzzzzz_1[i] * fz_be_0 + g_xxzzzz_zzzzz_1[i] * fe_0 + g_xxzzzz_xzzzzz_1[i] * pa_x[i];

        g_xxxzzzz_yyyyyy_0[i] = 2.0 * g_xzzzz_yyyyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyyyy_1[i] * fz_be_0 + g_xxzzzz_yyyyyy_1[i] * pa_x[i];

        g_xxxzzzz_yyyyyz_0[i] = 2.0 * g_xzzzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyyyz_1[i] * fz_be_0 + g_xxzzzz_yyyyyz_1[i] * pa_x[i];

        g_xxxzzzz_yyyyzz_0[i] = 2.0 * g_xzzzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyyzz_1[i] * fz_be_0 + g_xxzzzz_yyyyzz_1[i] * pa_x[i];

        g_xxxzzzz_yyyzzz_0[i] = 2.0 * g_xzzzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyzzz_1[i] * fz_be_0 + g_xxzzzz_yyyzzz_1[i] * pa_x[i];

        g_xxxzzzz_yyzzzz_0[i] = 2.0 * g_xzzzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyzzzz_1[i] * fz_be_0 + g_xxzzzz_yyzzzz_1[i] * pa_x[i];

        g_xxxzzzz_yzzzzz_0[i] = 2.0 * g_xzzzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yzzzzz_1[i] * fz_be_0 + g_xxzzzz_yzzzzz_1[i] * pa_x[i];

        g_xxxzzzz_zzzzzz_0[i] = 2.0 * g_xzzzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_zzzzzz_1[i] * fz_be_0 + g_xxzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : KI

    auto g_xxyyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 420);

    auto g_xxyyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 421);

    auto g_xxyyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 422);

    auto g_xxyyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 423);

    auto g_xxyyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 424);

    auto g_xxyyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 425);

    auto g_xxyyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 426);

    auto g_xxyyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 427);

    auto g_xxyyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 428);

    auto g_xxyyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 429);

    auto g_xxyyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 430);

    auto g_xxyyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 431);

    auto g_xxyyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 432);

    auto g_xxyyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 433);

    auto g_xxyyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 434);

    auto g_xxyyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 435);

    auto g_xxyyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 436);

    auto g_xxyyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 437);

    auto g_xxyyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 438);

    auto g_xxyyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 439);

    auto g_xxyyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 440);

    auto g_xxyyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 441);

    auto g_xxyyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 442);

    auto g_xxyyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 443);

    auto g_xxyyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 444);

    auto g_xxyyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 445);

    auto g_xxyyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 446);

    auto g_xxyyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 447);

    #pragma omp simd aligned(g_xxyyy_xxxxxx_0, g_xxyyy_xxxxxx_1, g_xxyyy_xxxxxz_0, g_xxyyy_xxxxxz_1, g_xxyyy_xxxxzz_0, g_xxyyy_xxxxzz_1, g_xxyyy_xxxzzz_0, g_xxyyy_xxxzzz_1, g_xxyyy_xxzzzz_0, g_xxyyy_xxzzzz_1, g_xxyyy_xzzzzz_0, g_xxyyy_xzzzzz_1, g_xxyyyy_xxxxxx_1, g_xxyyyy_xxxxxz_1, g_xxyyyy_xxxxzz_1, g_xxyyyy_xxxzzz_1, g_xxyyyy_xxzzzz_1, g_xxyyyy_xzzzzz_1, g_xxyyyyy_xxxxxx_0, g_xxyyyyy_xxxxxy_0, g_xxyyyyy_xxxxxz_0, g_xxyyyyy_xxxxyy_0, g_xxyyyyy_xxxxyz_0, g_xxyyyyy_xxxxzz_0, g_xxyyyyy_xxxyyy_0, g_xxyyyyy_xxxyyz_0, g_xxyyyyy_xxxyzz_0, g_xxyyyyy_xxxzzz_0, g_xxyyyyy_xxyyyy_0, g_xxyyyyy_xxyyyz_0, g_xxyyyyy_xxyyzz_0, g_xxyyyyy_xxyzzz_0, g_xxyyyyy_xxzzzz_0, g_xxyyyyy_xyyyyy_0, g_xxyyyyy_xyyyyz_0, g_xxyyyyy_xyyyzz_0, g_xxyyyyy_xyyzzz_0, g_xxyyyyy_xyzzzz_0, g_xxyyyyy_xzzzzz_0, g_xxyyyyy_yyyyyy_0, g_xxyyyyy_yyyyyz_0, g_xxyyyyy_yyyyzz_0, g_xxyyyyy_yyyzzz_0, g_xxyyyyy_yyzzzz_0, g_xxyyyyy_yzzzzz_0, g_xxyyyyy_zzzzzz_0, g_xyyyyy_xxxxxy_1, g_xyyyyy_xxxxy_1, g_xyyyyy_xxxxyy_1, g_xyyyyy_xxxxyz_1, g_xyyyyy_xxxyy_1, g_xyyyyy_xxxyyy_1, g_xyyyyy_xxxyyz_1, g_xyyyyy_xxxyz_1, g_xyyyyy_xxxyzz_1, g_xyyyyy_xxyyy_1, g_xyyyyy_xxyyyy_1, g_xyyyyy_xxyyyz_1, g_xyyyyy_xxyyz_1, g_xyyyyy_xxyyzz_1, g_xyyyyy_xxyzz_1, g_xyyyyy_xxyzzz_1, g_xyyyyy_xyyyy_1, g_xyyyyy_xyyyyy_1, g_xyyyyy_xyyyyz_1, g_xyyyyy_xyyyz_1, g_xyyyyy_xyyyzz_1, g_xyyyyy_xyyzz_1, g_xyyyyy_xyyzzz_1, g_xyyyyy_xyzzz_1, g_xyyyyy_xyzzzz_1, g_xyyyyy_yyyyy_1, g_xyyyyy_yyyyyy_1, g_xyyyyy_yyyyyz_1, g_xyyyyy_yyyyz_1, g_xyyyyy_yyyyzz_1, g_xyyyyy_yyyzz_1, g_xyyyyy_yyyzzz_1, g_xyyyyy_yyzzz_1, g_xyyyyy_yyzzzz_1, g_xyyyyy_yzzzz_1, g_xyyyyy_yzzzzz_1, g_xyyyyy_zzzzzz_1, g_yyyyy_xxxxxy_0, g_yyyyy_xxxxxy_1, g_yyyyy_xxxxyy_0, g_yyyyy_xxxxyy_1, g_yyyyy_xxxxyz_0, g_yyyyy_xxxxyz_1, g_yyyyy_xxxyyy_0, g_yyyyy_xxxyyy_1, g_yyyyy_xxxyyz_0, g_yyyyy_xxxyyz_1, g_yyyyy_xxxyzz_0, g_yyyyy_xxxyzz_1, g_yyyyy_xxyyyy_0, g_yyyyy_xxyyyy_1, g_yyyyy_xxyyyz_0, g_yyyyy_xxyyyz_1, g_yyyyy_xxyyzz_0, g_yyyyy_xxyyzz_1, g_yyyyy_xxyzzz_0, g_yyyyy_xxyzzz_1, g_yyyyy_xyyyyy_0, g_yyyyy_xyyyyy_1, g_yyyyy_xyyyyz_0, g_yyyyy_xyyyyz_1, g_yyyyy_xyyyzz_0, g_yyyyy_xyyyzz_1, g_yyyyy_xyyzzz_0, g_yyyyy_xyyzzz_1, g_yyyyy_xyzzzz_0, g_yyyyy_xyzzzz_1, g_yyyyy_yyyyyy_0, g_yyyyy_yyyyyy_1, g_yyyyy_yyyyyz_0, g_yyyyy_yyyyyz_1, g_yyyyy_yyyyzz_0, g_yyyyy_yyyyzz_1, g_yyyyy_yyyzzz_0, g_yyyyy_yyyzzz_1, g_yyyyy_yyzzzz_0, g_yyyyy_yyzzzz_1, g_yyyyy_yzzzzz_0, g_yyyyy_yzzzzz_1, g_yyyyy_zzzzzz_0, g_yyyyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyyy_xxxxxx_0[i] = 4.0 * g_xxyyy_xxxxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxxxx_1[i] * fz_be_0 + g_xxyyyy_xxxxxx_1[i] * pa_y[i];

        g_xxyyyyy_xxxxxy_0[i] = g_yyyyy_xxxxxy_0[i] * fbe_0 - g_yyyyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyyyyy_xxxxy_1[i] * fe_0 + g_xyyyyy_xxxxxy_1[i] * pa_x[i];

        g_xxyyyyy_xxxxxz_0[i] = 4.0 * g_xxyyy_xxxxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxxxz_1[i] * fz_be_0 + g_xxyyyy_xxxxxz_1[i] * pa_y[i];

        g_xxyyyyy_xxxxyy_0[i] = g_yyyyy_xxxxyy_0[i] * fbe_0 - g_yyyyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyyyyy_xxxyy_1[i] * fe_0 + g_xyyyyy_xxxxyy_1[i] * pa_x[i];

        g_xxyyyyy_xxxxyz_0[i] = g_yyyyy_xxxxyz_0[i] * fbe_0 - g_yyyyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_xxxyz_1[i] * fe_0 + g_xyyyyy_xxxxyz_1[i] * pa_x[i];

        g_xxyyyyy_xxxxzz_0[i] = 4.0 * g_xxyyy_xxxxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxxzz_1[i] * fz_be_0 + g_xxyyyy_xxxxzz_1[i] * pa_y[i];

        g_xxyyyyy_xxxyyy_0[i] = g_yyyyy_xxxyyy_0[i] * fbe_0 - g_yyyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_xxyyy_1[i] * fe_0 + g_xyyyyy_xxxyyy_1[i] * pa_x[i];

        g_xxyyyyy_xxxyyz_0[i] = g_yyyyy_xxxyyz_0[i] * fbe_0 - g_yyyyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_xxyyz_1[i] * fe_0 + g_xyyyyy_xxxyyz_1[i] * pa_x[i];

        g_xxyyyyy_xxxyzz_0[i] = g_yyyyy_xxxyzz_0[i] * fbe_0 - g_yyyyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_xxyzz_1[i] * fe_0 + g_xyyyyy_xxxyzz_1[i] * pa_x[i];

        g_xxyyyyy_xxxzzz_0[i] = 4.0 * g_xxyyy_xxxzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxzzz_1[i] * fz_be_0 + g_xxyyyy_xxxzzz_1[i] * pa_y[i];

        g_xxyyyyy_xxyyyy_0[i] = g_yyyyy_xxyyyy_0[i] * fbe_0 - g_yyyyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyyyy_1[i] * fe_0 + g_xyyyyy_xxyyyy_1[i] * pa_x[i];

        g_xxyyyyy_xxyyyz_0[i] = g_yyyyy_xxyyyz_0[i] * fbe_0 - g_yyyyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyyyz_1[i] * fe_0 + g_xyyyyy_xxyyyz_1[i] * pa_x[i];

        g_xxyyyyy_xxyyzz_0[i] = g_yyyyy_xxyyzz_0[i] * fbe_0 - g_yyyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyyzz_1[i] * fe_0 + g_xyyyyy_xxyyzz_1[i] * pa_x[i];

        g_xxyyyyy_xxyzzz_0[i] = g_yyyyy_xxyzzz_0[i] * fbe_0 - g_yyyyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyzzz_1[i] * fe_0 + g_xyyyyy_xxyzzz_1[i] * pa_x[i];

        g_xxyyyyy_xxzzzz_0[i] = 4.0 * g_xxyyy_xxzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxzzzz_1[i] * fz_be_0 + g_xxyyyy_xxzzzz_1[i] * pa_y[i];

        g_xxyyyyy_xyyyyy_0[i] = g_yyyyy_xyyyyy_0[i] * fbe_0 - g_yyyyy_xyyyyy_1[i] * fz_be_0 + g_xyyyyy_yyyyy_1[i] * fe_0 + g_xyyyyy_xyyyyy_1[i] * pa_x[i];

        g_xxyyyyy_xyyyyz_0[i] = g_yyyyy_xyyyyz_0[i] * fbe_0 - g_yyyyy_xyyyyz_1[i] * fz_be_0 + g_xyyyyy_yyyyz_1[i] * fe_0 + g_xyyyyy_xyyyyz_1[i] * pa_x[i];

        g_xxyyyyy_xyyyzz_0[i] = g_yyyyy_xyyyzz_0[i] * fbe_0 - g_yyyyy_xyyyzz_1[i] * fz_be_0 + g_xyyyyy_yyyzz_1[i] * fe_0 + g_xyyyyy_xyyyzz_1[i] * pa_x[i];

        g_xxyyyyy_xyyzzz_0[i] = g_yyyyy_xyyzzz_0[i] * fbe_0 - g_yyyyy_xyyzzz_1[i] * fz_be_0 + g_xyyyyy_yyzzz_1[i] * fe_0 + g_xyyyyy_xyyzzz_1[i] * pa_x[i];

        g_xxyyyyy_xyzzzz_0[i] = g_yyyyy_xyzzzz_0[i] * fbe_0 - g_yyyyy_xyzzzz_1[i] * fz_be_0 + g_xyyyyy_yzzzz_1[i] * fe_0 + g_xyyyyy_xyzzzz_1[i] * pa_x[i];

        g_xxyyyyy_xzzzzz_0[i] = 4.0 * g_xxyyy_xzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xzzzzz_1[i] * fz_be_0 + g_xxyyyy_xzzzzz_1[i] * pa_y[i];

        g_xxyyyyy_yyyyyy_0[i] = g_yyyyy_yyyyyy_0[i] * fbe_0 - g_yyyyy_yyyyyy_1[i] * fz_be_0 + g_xyyyyy_yyyyyy_1[i] * pa_x[i];

        g_xxyyyyy_yyyyyz_0[i] = g_yyyyy_yyyyyz_0[i] * fbe_0 - g_yyyyy_yyyyyz_1[i] * fz_be_0 + g_xyyyyy_yyyyyz_1[i] * pa_x[i];

        g_xxyyyyy_yyyyzz_0[i] = g_yyyyy_yyyyzz_0[i] * fbe_0 - g_yyyyy_yyyyzz_1[i] * fz_be_0 + g_xyyyyy_yyyyzz_1[i] * pa_x[i];

        g_xxyyyyy_yyyzzz_0[i] = g_yyyyy_yyyzzz_0[i] * fbe_0 - g_yyyyy_yyyzzz_1[i] * fz_be_0 + g_xyyyyy_yyyzzz_1[i] * pa_x[i];

        g_xxyyyyy_yyzzzz_0[i] = g_yyyyy_yyzzzz_0[i] * fbe_0 - g_yyyyy_yyzzzz_1[i] * fz_be_0 + g_xyyyyy_yyzzzz_1[i] * pa_x[i];

        g_xxyyyyy_yzzzzz_0[i] = g_yyyyy_yzzzzz_0[i] * fbe_0 - g_yyyyy_yzzzzz_1[i] * fz_be_0 + g_xyyyyy_yzzzzz_1[i] * pa_x[i];

        g_xxyyyyy_zzzzzz_0[i] = g_yyyyy_zzzzzz_0[i] * fbe_0 - g_yyyyy_zzzzzz_1[i] * fz_be_0 + g_xyyyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 448-476 components of targeted buffer : KI

    auto g_xxyyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 448);

    auto g_xxyyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 449);

    auto g_xxyyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 450);

    auto g_xxyyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 451);

    auto g_xxyyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 452);

    auto g_xxyyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 453);

    auto g_xxyyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 454);

    auto g_xxyyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 455);

    auto g_xxyyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 456);

    auto g_xxyyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 457);

    auto g_xxyyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 458);

    auto g_xxyyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 459);

    auto g_xxyyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 460);

    auto g_xxyyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 461);

    auto g_xxyyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 462);

    auto g_xxyyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 463);

    auto g_xxyyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 464);

    auto g_xxyyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 465);

    auto g_xxyyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 466);

    auto g_xxyyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 467);

    auto g_xxyyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 468);

    auto g_xxyyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 469);

    auto g_xxyyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 470);

    auto g_xxyyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 471);

    auto g_xxyyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 472);

    auto g_xxyyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 473);

    auto g_xxyyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 474);

    auto g_xxyyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 475);

    #pragma omp simd aligned(g_xxyyyy_xxxxx_1, g_xxyyyy_xxxxxx_1, g_xxyyyy_xxxxxy_1, g_xxyyyy_xxxxxz_1, g_xxyyyy_xxxxy_1, g_xxyyyy_xxxxyy_1, g_xxyyyy_xxxxyz_1, g_xxyyyy_xxxxz_1, g_xxyyyy_xxxxzz_1, g_xxyyyy_xxxyy_1, g_xxyyyy_xxxyyy_1, g_xxyyyy_xxxyyz_1, g_xxyyyy_xxxyz_1, g_xxyyyy_xxxyzz_1, g_xxyyyy_xxxzz_1, g_xxyyyy_xxxzzz_1, g_xxyyyy_xxyyy_1, g_xxyyyy_xxyyyy_1, g_xxyyyy_xxyyyz_1, g_xxyyyy_xxyyz_1, g_xxyyyy_xxyyzz_1, g_xxyyyy_xxyzz_1, g_xxyyyy_xxyzzz_1, g_xxyyyy_xxzzz_1, g_xxyyyy_xxzzzz_1, g_xxyyyy_xyyyy_1, g_xxyyyy_xyyyyy_1, g_xxyyyy_xyyyyz_1, g_xxyyyy_xyyyz_1, g_xxyyyy_xyyyzz_1, g_xxyyyy_xyyzz_1, g_xxyyyy_xyyzzz_1, g_xxyyyy_xyzzz_1, g_xxyyyy_xyzzzz_1, g_xxyyyy_xzzzz_1, g_xxyyyy_xzzzzz_1, g_xxyyyy_yyyyy_1, g_xxyyyy_yyyyyy_1, g_xxyyyy_yyyyyz_1, g_xxyyyy_yyyyz_1, g_xxyyyy_yyyyzz_1, g_xxyyyy_yyyzz_1, g_xxyyyy_yyyzzz_1, g_xxyyyy_yyzzz_1, g_xxyyyy_yyzzzz_1, g_xxyyyy_yzzzz_1, g_xxyyyy_yzzzzz_1, g_xxyyyy_zzzzz_1, g_xxyyyy_zzzzzz_1, g_xxyyyyz_xxxxxx_0, g_xxyyyyz_xxxxxy_0, g_xxyyyyz_xxxxxz_0, g_xxyyyyz_xxxxyy_0, g_xxyyyyz_xxxxyz_0, g_xxyyyyz_xxxxzz_0, g_xxyyyyz_xxxyyy_0, g_xxyyyyz_xxxyyz_0, g_xxyyyyz_xxxyzz_0, g_xxyyyyz_xxxzzz_0, g_xxyyyyz_xxyyyy_0, g_xxyyyyz_xxyyyz_0, g_xxyyyyz_xxyyzz_0, g_xxyyyyz_xxyzzz_0, g_xxyyyyz_xxzzzz_0, g_xxyyyyz_xyyyyy_0, g_xxyyyyz_xyyyyz_0, g_xxyyyyz_xyyyzz_0, g_xxyyyyz_xyyzzz_0, g_xxyyyyz_xyzzzz_0, g_xxyyyyz_xzzzzz_0, g_xxyyyyz_yyyyyy_0, g_xxyyyyz_yyyyyz_0, g_xxyyyyz_yyyyzz_0, g_xxyyyyz_yyyzzz_0, g_xxyyyyz_yyzzzz_0, g_xxyyyyz_yzzzzz_0, g_xxyyyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyyz_xxxxxx_0[i] = g_xxyyyy_xxxxxx_1[i] * pa_z[i];

        g_xxyyyyz_xxxxxy_0[i] = g_xxyyyy_xxxxxy_1[i] * pa_z[i];

        g_xxyyyyz_xxxxxz_0[i] = g_xxyyyy_xxxxx_1[i] * fe_0 + g_xxyyyy_xxxxxz_1[i] * pa_z[i];

        g_xxyyyyz_xxxxyy_0[i] = g_xxyyyy_xxxxyy_1[i] * pa_z[i];

        g_xxyyyyz_xxxxyz_0[i] = g_xxyyyy_xxxxy_1[i] * fe_0 + g_xxyyyy_xxxxyz_1[i] * pa_z[i];

        g_xxyyyyz_xxxxzz_0[i] = 2.0 * g_xxyyyy_xxxxz_1[i] * fe_0 + g_xxyyyy_xxxxzz_1[i] * pa_z[i];

        g_xxyyyyz_xxxyyy_0[i] = g_xxyyyy_xxxyyy_1[i] * pa_z[i];

        g_xxyyyyz_xxxyyz_0[i] = g_xxyyyy_xxxyy_1[i] * fe_0 + g_xxyyyy_xxxyyz_1[i] * pa_z[i];

        g_xxyyyyz_xxxyzz_0[i] = 2.0 * g_xxyyyy_xxxyz_1[i] * fe_0 + g_xxyyyy_xxxyzz_1[i] * pa_z[i];

        g_xxyyyyz_xxxzzz_0[i] = 3.0 * g_xxyyyy_xxxzz_1[i] * fe_0 + g_xxyyyy_xxxzzz_1[i] * pa_z[i];

        g_xxyyyyz_xxyyyy_0[i] = g_xxyyyy_xxyyyy_1[i] * pa_z[i];

        g_xxyyyyz_xxyyyz_0[i] = g_xxyyyy_xxyyy_1[i] * fe_0 + g_xxyyyy_xxyyyz_1[i] * pa_z[i];

        g_xxyyyyz_xxyyzz_0[i] = 2.0 * g_xxyyyy_xxyyz_1[i] * fe_0 + g_xxyyyy_xxyyzz_1[i] * pa_z[i];

        g_xxyyyyz_xxyzzz_0[i] = 3.0 * g_xxyyyy_xxyzz_1[i] * fe_0 + g_xxyyyy_xxyzzz_1[i] * pa_z[i];

        g_xxyyyyz_xxzzzz_0[i] = 4.0 * g_xxyyyy_xxzzz_1[i] * fe_0 + g_xxyyyy_xxzzzz_1[i] * pa_z[i];

        g_xxyyyyz_xyyyyy_0[i] = g_xxyyyy_xyyyyy_1[i] * pa_z[i];

        g_xxyyyyz_xyyyyz_0[i] = g_xxyyyy_xyyyy_1[i] * fe_0 + g_xxyyyy_xyyyyz_1[i] * pa_z[i];

        g_xxyyyyz_xyyyzz_0[i] = 2.0 * g_xxyyyy_xyyyz_1[i] * fe_0 + g_xxyyyy_xyyyzz_1[i] * pa_z[i];

        g_xxyyyyz_xyyzzz_0[i] = 3.0 * g_xxyyyy_xyyzz_1[i] * fe_0 + g_xxyyyy_xyyzzz_1[i] * pa_z[i];

        g_xxyyyyz_xyzzzz_0[i] = 4.0 * g_xxyyyy_xyzzz_1[i] * fe_0 + g_xxyyyy_xyzzzz_1[i] * pa_z[i];

        g_xxyyyyz_xzzzzz_0[i] = 5.0 * g_xxyyyy_xzzzz_1[i] * fe_0 + g_xxyyyy_xzzzzz_1[i] * pa_z[i];

        g_xxyyyyz_yyyyyy_0[i] = g_xxyyyy_yyyyyy_1[i] * pa_z[i];

        g_xxyyyyz_yyyyyz_0[i] = g_xxyyyy_yyyyy_1[i] * fe_0 + g_xxyyyy_yyyyyz_1[i] * pa_z[i];

        g_xxyyyyz_yyyyzz_0[i] = 2.0 * g_xxyyyy_yyyyz_1[i] * fe_0 + g_xxyyyy_yyyyzz_1[i] * pa_z[i];

        g_xxyyyyz_yyyzzz_0[i] = 3.0 * g_xxyyyy_yyyzz_1[i] * fe_0 + g_xxyyyy_yyyzzz_1[i] * pa_z[i];

        g_xxyyyyz_yyzzzz_0[i] = 4.0 * g_xxyyyy_yyzzz_1[i] * fe_0 + g_xxyyyy_yyzzzz_1[i] * pa_z[i];

        g_xxyyyyz_yzzzzz_0[i] = 5.0 * g_xxyyyy_yzzzz_1[i] * fe_0 + g_xxyyyy_yzzzzz_1[i] * pa_z[i];

        g_xxyyyyz_zzzzzz_0[i] = 6.0 * g_xxyyyy_zzzzz_1[i] * fe_0 + g_xxyyyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 476-504 components of targeted buffer : KI

    auto g_xxyyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 476);

    auto g_xxyyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 477);

    auto g_xxyyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 478);

    auto g_xxyyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 479);

    auto g_xxyyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 480);

    auto g_xxyyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 481);

    auto g_xxyyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 482);

    auto g_xxyyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 483);

    auto g_xxyyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 484);

    auto g_xxyyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 485);

    auto g_xxyyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 486);

    auto g_xxyyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 487);

    auto g_xxyyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 488);

    auto g_xxyyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 489);

    auto g_xxyyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 490);

    auto g_xxyyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 491);

    auto g_xxyyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 492);

    auto g_xxyyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 493);

    auto g_xxyyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 494);

    auto g_xxyyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 495);

    auto g_xxyyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 496);

    auto g_xxyyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 497);

    auto g_xxyyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 498);

    auto g_xxyyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 499);

    auto g_xxyyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 500);

    auto g_xxyyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 501);

    auto g_xxyyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 502);

    auto g_xxyyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 503);

    #pragma omp simd aligned(g_xxyyy_xxxxxy_0, g_xxyyy_xxxxxy_1, g_xxyyy_xxxxyy_0, g_xxyyy_xxxxyy_1, g_xxyyy_xxxyyy_0, g_xxyyy_xxxyyy_1, g_xxyyy_xxyyyy_0, g_xxyyy_xxyyyy_1, g_xxyyy_xyyyyy_0, g_xxyyy_xyyyyy_1, g_xxyyyz_xxxxxy_1, g_xxyyyz_xxxxyy_1, g_xxyyyz_xxxyyy_1, g_xxyyyz_xxyyyy_1, g_xxyyyz_xyyyyy_1, g_xxyyyzz_xxxxxx_0, g_xxyyyzz_xxxxxy_0, g_xxyyyzz_xxxxxz_0, g_xxyyyzz_xxxxyy_0, g_xxyyyzz_xxxxyz_0, g_xxyyyzz_xxxxzz_0, g_xxyyyzz_xxxyyy_0, g_xxyyyzz_xxxyyz_0, g_xxyyyzz_xxxyzz_0, g_xxyyyzz_xxxzzz_0, g_xxyyyzz_xxyyyy_0, g_xxyyyzz_xxyyyz_0, g_xxyyyzz_xxyyzz_0, g_xxyyyzz_xxyzzz_0, g_xxyyyzz_xxzzzz_0, g_xxyyyzz_xyyyyy_0, g_xxyyyzz_xyyyyz_0, g_xxyyyzz_xyyyzz_0, g_xxyyyzz_xyyzzz_0, g_xxyyyzz_xyzzzz_0, g_xxyyyzz_xzzzzz_0, g_xxyyyzz_yyyyyy_0, g_xxyyyzz_yyyyyz_0, g_xxyyyzz_yyyyzz_0, g_xxyyyzz_yyyzzz_0, g_xxyyyzz_yyzzzz_0, g_xxyyyzz_yzzzzz_0, g_xxyyyzz_zzzzzz_0, g_xxyyzz_xxxxxx_1, g_xxyyzz_xxxxxz_1, g_xxyyzz_xxxxzz_1, g_xxyyzz_xxxzzz_1, g_xxyyzz_xxzzzz_1, g_xxyyzz_xzzzzz_1, g_xxyzz_xxxxxx_0, g_xxyzz_xxxxxx_1, g_xxyzz_xxxxxz_0, g_xxyzz_xxxxxz_1, g_xxyzz_xxxxzz_0, g_xxyzz_xxxxzz_1, g_xxyzz_xxxzzz_0, g_xxyzz_xxxzzz_1, g_xxyzz_xxzzzz_0, g_xxyzz_xxzzzz_1, g_xxyzz_xzzzzz_0, g_xxyzz_xzzzzz_1, g_xyyyzz_xxxxyz_1, g_xyyyzz_xxxyyz_1, g_xyyyzz_xxxyz_1, g_xyyyzz_xxxyzz_1, g_xyyyzz_xxyyyz_1, g_xyyyzz_xxyyz_1, g_xyyyzz_xxyyzz_1, g_xyyyzz_xxyzz_1, g_xyyyzz_xxyzzz_1, g_xyyyzz_xyyyyz_1, g_xyyyzz_xyyyz_1, g_xyyyzz_xyyyzz_1, g_xyyyzz_xyyzz_1, g_xyyyzz_xyyzzz_1, g_xyyyzz_xyzzz_1, g_xyyyzz_xyzzzz_1, g_xyyyzz_yyyyyy_1, g_xyyyzz_yyyyyz_1, g_xyyyzz_yyyyz_1, g_xyyyzz_yyyyzz_1, g_xyyyzz_yyyzz_1, g_xyyyzz_yyyzzz_1, g_xyyyzz_yyzzz_1, g_xyyyzz_yyzzzz_1, g_xyyyzz_yzzzz_1, g_xyyyzz_yzzzzz_1, g_xyyyzz_zzzzzz_1, g_yyyzz_xxxxyz_0, g_yyyzz_xxxxyz_1, g_yyyzz_xxxyyz_0, g_yyyzz_xxxyyz_1, g_yyyzz_xxxyzz_0, g_yyyzz_xxxyzz_1, g_yyyzz_xxyyyz_0, g_yyyzz_xxyyyz_1, g_yyyzz_xxyyzz_0, g_yyyzz_xxyyzz_1, g_yyyzz_xxyzzz_0, g_yyyzz_xxyzzz_1, g_yyyzz_xyyyyz_0, g_yyyzz_xyyyyz_1, g_yyyzz_xyyyzz_0, g_yyyzz_xyyyzz_1, g_yyyzz_xyyzzz_0, g_yyyzz_xyyzzz_1, g_yyyzz_xyzzzz_0, g_yyyzz_xyzzzz_1, g_yyyzz_yyyyyy_0, g_yyyzz_yyyyyy_1, g_yyyzz_yyyyyz_0, g_yyyzz_yyyyyz_1, g_yyyzz_yyyyzz_0, g_yyyzz_yyyyzz_1, g_yyyzz_yyyzzz_0, g_yyyzz_yyyzzz_1, g_yyyzz_yyzzzz_0, g_yyyzz_yyzzzz_1, g_yyyzz_yzzzzz_0, g_yyyzz_yzzzzz_1, g_yyyzz_zzzzzz_0, g_yyyzz_zzzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyzz_xxxxxx_0[i] = 2.0 * g_xxyzz_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxxxx_1[i] * fz_be_0 + g_xxyyzz_xxxxxx_1[i] * pa_y[i];

        g_xxyyyzz_xxxxxy_0[i] = g_xxyyy_xxxxxy_0[i] * fbe_0 - g_xxyyy_xxxxxy_1[i] * fz_be_0 + g_xxyyyz_xxxxxy_1[i] * pa_z[i];

        g_xxyyyzz_xxxxxz_0[i] = 2.0 * g_xxyzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxxxz_1[i] * fz_be_0 + g_xxyyzz_xxxxxz_1[i] * pa_y[i];

        g_xxyyyzz_xxxxyy_0[i] = g_xxyyy_xxxxyy_0[i] * fbe_0 - g_xxyyy_xxxxyy_1[i] * fz_be_0 + g_xxyyyz_xxxxyy_1[i] * pa_z[i];

        g_xxyyyzz_xxxxyz_0[i] = g_yyyzz_xxxxyz_0[i] * fbe_0 - g_yyyzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_xxxyz_1[i] * fe_0 + g_xyyyzz_xxxxyz_1[i] * pa_x[i];

        g_xxyyyzz_xxxxzz_0[i] = 2.0 * g_xxyzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxxzz_1[i] * fz_be_0 + g_xxyyzz_xxxxzz_1[i] * pa_y[i];

        g_xxyyyzz_xxxyyy_0[i] = g_xxyyy_xxxyyy_0[i] * fbe_0 - g_xxyyy_xxxyyy_1[i] * fz_be_0 + g_xxyyyz_xxxyyy_1[i] * pa_z[i];

        g_xxyyyzz_xxxyyz_0[i] = g_yyyzz_xxxyyz_0[i] * fbe_0 - g_yyyzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_xxyyz_1[i] * fe_0 + g_xyyyzz_xxxyyz_1[i] * pa_x[i];

        g_xxyyyzz_xxxyzz_0[i] = g_yyyzz_xxxyzz_0[i] * fbe_0 - g_yyyzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_xxyzz_1[i] * fe_0 + g_xyyyzz_xxxyzz_1[i] * pa_x[i];

        g_xxyyyzz_xxxzzz_0[i] = 2.0 * g_xxyzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxzzz_1[i] * fz_be_0 + g_xxyyzz_xxxzzz_1[i] * pa_y[i];

        g_xxyyyzz_xxyyyy_0[i] = g_xxyyy_xxyyyy_0[i] * fbe_0 - g_xxyyy_xxyyyy_1[i] * fz_be_0 + g_xxyyyz_xxyyyy_1[i] * pa_z[i];

        g_xxyyyzz_xxyyyz_0[i] = g_yyyzz_xxyyyz_0[i] * fbe_0 - g_yyyzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_xyyyz_1[i] * fe_0 + g_xyyyzz_xxyyyz_1[i] * pa_x[i];

        g_xxyyyzz_xxyyzz_0[i] = g_yyyzz_xxyyzz_0[i] * fbe_0 - g_yyyzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_xyyzz_1[i] * fe_0 + g_xyyyzz_xxyyzz_1[i] * pa_x[i];

        g_xxyyyzz_xxyzzz_0[i] = g_yyyzz_xxyzzz_0[i] * fbe_0 - g_yyyzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_xyzzz_1[i] * fe_0 + g_xyyyzz_xxyzzz_1[i] * pa_x[i];

        g_xxyyyzz_xxzzzz_0[i] = 2.0 * g_xxyzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxzzzz_1[i] * fz_be_0 + g_xxyyzz_xxzzzz_1[i] * pa_y[i];

        g_xxyyyzz_xyyyyy_0[i] = g_xxyyy_xyyyyy_0[i] * fbe_0 - g_xxyyy_xyyyyy_1[i] * fz_be_0 + g_xxyyyz_xyyyyy_1[i] * pa_z[i];

        g_xxyyyzz_xyyyyz_0[i] = g_yyyzz_xyyyyz_0[i] * fbe_0 - g_yyyzz_xyyyyz_1[i] * fz_be_0 + g_xyyyzz_yyyyz_1[i] * fe_0 + g_xyyyzz_xyyyyz_1[i] * pa_x[i];

        g_xxyyyzz_xyyyzz_0[i] = g_yyyzz_xyyyzz_0[i] * fbe_0 - g_yyyzz_xyyyzz_1[i] * fz_be_0 + g_xyyyzz_yyyzz_1[i] * fe_0 + g_xyyyzz_xyyyzz_1[i] * pa_x[i];

        g_xxyyyzz_xyyzzz_0[i] = g_yyyzz_xyyzzz_0[i] * fbe_0 - g_yyyzz_xyyzzz_1[i] * fz_be_0 + g_xyyyzz_yyzzz_1[i] * fe_0 + g_xyyyzz_xyyzzz_1[i] * pa_x[i];

        g_xxyyyzz_xyzzzz_0[i] = g_yyyzz_xyzzzz_0[i] * fbe_0 - g_yyyzz_xyzzzz_1[i] * fz_be_0 + g_xyyyzz_yzzzz_1[i] * fe_0 + g_xyyyzz_xyzzzz_1[i] * pa_x[i];

        g_xxyyyzz_xzzzzz_0[i] = 2.0 * g_xxyzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xzzzzz_1[i] * fz_be_0 + g_xxyyzz_xzzzzz_1[i] * pa_y[i];

        g_xxyyyzz_yyyyyy_0[i] = g_yyyzz_yyyyyy_0[i] * fbe_0 - g_yyyzz_yyyyyy_1[i] * fz_be_0 + g_xyyyzz_yyyyyy_1[i] * pa_x[i];

        g_xxyyyzz_yyyyyz_0[i] = g_yyyzz_yyyyyz_0[i] * fbe_0 - g_yyyzz_yyyyyz_1[i] * fz_be_0 + g_xyyyzz_yyyyyz_1[i] * pa_x[i];

        g_xxyyyzz_yyyyzz_0[i] = g_yyyzz_yyyyzz_0[i] * fbe_0 - g_yyyzz_yyyyzz_1[i] * fz_be_0 + g_xyyyzz_yyyyzz_1[i] * pa_x[i];

        g_xxyyyzz_yyyzzz_0[i] = g_yyyzz_yyyzzz_0[i] * fbe_0 - g_yyyzz_yyyzzz_1[i] * fz_be_0 + g_xyyyzz_yyyzzz_1[i] * pa_x[i];

        g_xxyyyzz_yyzzzz_0[i] = g_yyyzz_yyzzzz_0[i] * fbe_0 - g_yyyzz_yyzzzz_1[i] * fz_be_0 + g_xyyyzz_yyzzzz_1[i] * pa_x[i];

        g_xxyyyzz_yzzzzz_0[i] = g_yyyzz_yzzzzz_0[i] * fbe_0 - g_yyyzz_yzzzzz_1[i] * fz_be_0 + g_xyyyzz_yzzzzz_1[i] * pa_x[i];

        g_xxyyyzz_zzzzzz_0[i] = g_yyyzz_zzzzzz_0[i] * fbe_0 - g_yyyzz_zzzzzz_1[i] * fz_be_0 + g_xyyyzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 504-532 components of targeted buffer : KI

    auto g_xxyyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 504);

    auto g_xxyyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 505);

    auto g_xxyyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 506);

    auto g_xxyyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 507);

    auto g_xxyyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 508);

    auto g_xxyyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 509);

    auto g_xxyyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 510);

    auto g_xxyyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 511);

    auto g_xxyyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 512);

    auto g_xxyyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 513);

    auto g_xxyyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 514);

    auto g_xxyyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 515);

    auto g_xxyyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 516);

    auto g_xxyyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 517);

    auto g_xxyyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 518);

    auto g_xxyyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 519);

    auto g_xxyyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 520);

    auto g_xxyyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 521);

    auto g_xxyyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 522);

    auto g_xxyyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 523);

    auto g_xxyyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 524);

    auto g_xxyyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 525);

    auto g_xxyyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 526);

    auto g_xxyyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 527);

    auto g_xxyyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 528);

    auto g_xxyyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 529);

    auto g_xxyyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 530);

    auto g_xxyyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 531);

    #pragma omp simd aligned(g_xxyyz_xxxxxy_0, g_xxyyz_xxxxxy_1, g_xxyyz_xxxxyy_0, g_xxyyz_xxxxyy_1, g_xxyyz_xxxyyy_0, g_xxyyz_xxxyyy_1, g_xxyyz_xxyyyy_0, g_xxyyz_xxyyyy_1, g_xxyyz_xyyyyy_0, g_xxyyz_xyyyyy_1, g_xxyyzz_xxxxxy_1, g_xxyyzz_xxxxyy_1, g_xxyyzz_xxxyyy_1, g_xxyyzz_xxyyyy_1, g_xxyyzz_xyyyyy_1, g_xxyyzzz_xxxxxx_0, g_xxyyzzz_xxxxxy_0, g_xxyyzzz_xxxxxz_0, g_xxyyzzz_xxxxyy_0, g_xxyyzzz_xxxxyz_0, g_xxyyzzz_xxxxzz_0, g_xxyyzzz_xxxyyy_0, g_xxyyzzz_xxxyyz_0, g_xxyyzzz_xxxyzz_0, g_xxyyzzz_xxxzzz_0, g_xxyyzzz_xxyyyy_0, g_xxyyzzz_xxyyyz_0, g_xxyyzzz_xxyyzz_0, g_xxyyzzz_xxyzzz_0, g_xxyyzzz_xxzzzz_0, g_xxyyzzz_xyyyyy_0, g_xxyyzzz_xyyyyz_0, g_xxyyzzz_xyyyzz_0, g_xxyyzzz_xyyzzz_0, g_xxyyzzz_xyzzzz_0, g_xxyyzzz_xzzzzz_0, g_xxyyzzz_yyyyyy_0, g_xxyyzzz_yyyyyz_0, g_xxyyzzz_yyyyzz_0, g_xxyyzzz_yyyzzz_0, g_xxyyzzz_yyzzzz_0, g_xxyyzzz_yzzzzz_0, g_xxyyzzz_zzzzzz_0, g_xxyzzz_xxxxxx_1, g_xxyzzz_xxxxxz_1, g_xxyzzz_xxxxzz_1, g_xxyzzz_xxxzzz_1, g_xxyzzz_xxzzzz_1, g_xxyzzz_xzzzzz_1, g_xxzzz_xxxxxx_0, g_xxzzz_xxxxxx_1, g_xxzzz_xxxxxz_0, g_xxzzz_xxxxxz_1, g_xxzzz_xxxxzz_0, g_xxzzz_xxxxzz_1, g_xxzzz_xxxzzz_0, g_xxzzz_xxxzzz_1, g_xxzzz_xxzzzz_0, g_xxzzz_xxzzzz_1, g_xxzzz_xzzzzz_0, g_xxzzz_xzzzzz_1, g_xyyzzz_xxxxyz_1, g_xyyzzz_xxxyyz_1, g_xyyzzz_xxxyz_1, g_xyyzzz_xxxyzz_1, g_xyyzzz_xxyyyz_1, g_xyyzzz_xxyyz_1, g_xyyzzz_xxyyzz_1, g_xyyzzz_xxyzz_1, g_xyyzzz_xxyzzz_1, g_xyyzzz_xyyyyz_1, g_xyyzzz_xyyyz_1, g_xyyzzz_xyyyzz_1, g_xyyzzz_xyyzz_1, g_xyyzzz_xyyzzz_1, g_xyyzzz_xyzzz_1, g_xyyzzz_xyzzzz_1, g_xyyzzz_yyyyyy_1, g_xyyzzz_yyyyyz_1, g_xyyzzz_yyyyz_1, g_xyyzzz_yyyyzz_1, g_xyyzzz_yyyzz_1, g_xyyzzz_yyyzzz_1, g_xyyzzz_yyzzz_1, g_xyyzzz_yyzzzz_1, g_xyyzzz_yzzzz_1, g_xyyzzz_yzzzzz_1, g_xyyzzz_zzzzzz_1, g_yyzzz_xxxxyz_0, g_yyzzz_xxxxyz_1, g_yyzzz_xxxyyz_0, g_yyzzz_xxxyyz_1, g_yyzzz_xxxyzz_0, g_yyzzz_xxxyzz_1, g_yyzzz_xxyyyz_0, g_yyzzz_xxyyyz_1, g_yyzzz_xxyyzz_0, g_yyzzz_xxyyzz_1, g_yyzzz_xxyzzz_0, g_yyzzz_xxyzzz_1, g_yyzzz_xyyyyz_0, g_yyzzz_xyyyyz_1, g_yyzzz_xyyyzz_0, g_yyzzz_xyyyzz_1, g_yyzzz_xyyzzz_0, g_yyzzz_xyyzzz_1, g_yyzzz_xyzzzz_0, g_yyzzz_xyzzzz_1, g_yyzzz_yyyyyy_0, g_yyzzz_yyyyyy_1, g_yyzzz_yyyyyz_0, g_yyzzz_yyyyyz_1, g_yyzzz_yyyyzz_0, g_yyzzz_yyyyzz_1, g_yyzzz_yyyzzz_0, g_yyzzz_yyyzzz_1, g_yyzzz_yyzzzz_0, g_yyzzz_yyzzzz_1, g_yyzzz_yzzzzz_0, g_yyzzz_yzzzzz_1, g_yyzzz_zzzzzz_0, g_yyzzz_zzzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzzz_xxxxxx_0[i] = g_xxzzz_xxxxxx_0[i] * fbe_0 - g_xxzzz_xxxxxx_1[i] * fz_be_0 + g_xxyzzz_xxxxxx_1[i] * pa_y[i];

        g_xxyyzzz_xxxxxy_0[i] = 2.0 * g_xxyyz_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxxxxy_1[i] * fz_be_0 + g_xxyyzz_xxxxxy_1[i] * pa_z[i];

        g_xxyyzzz_xxxxxz_0[i] = g_xxzzz_xxxxxz_0[i] * fbe_0 - g_xxzzz_xxxxxz_1[i] * fz_be_0 + g_xxyzzz_xxxxxz_1[i] * pa_y[i];

        g_xxyyzzz_xxxxyy_0[i] = 2.0 * g_xxyyz_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxxxyy_1[i] * fz_be_0 + g_xxyyzz_xxxxyy_1[i] * pa_z[i];

        g_xxyyzzz_xxxxyz_0[i] = g_yyzzz_xxxxyz_0[i] * fbe_0 - g_yyzzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_xxxyz_1[i] * fe_0 + g_xyyzzz_xxxxyz_1[i] * pa_x[i];

        g_xxyyzzz_xxxxzz_0[i] = g_xxzzz_xxxxzz_0[i] * fbe_0 - g_xxzzz_xxxxzz_1[i] * fz_be_0 + g_xxyzzz_xxxxzz_1[i] * pa_y[i];

        g_xxyyzzz_xxxyyy_0[i] = 2.0 * g_xxyyz_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxxyyy_1[i] * fz_be_0 + g_xxyyzz_xxxyyy_1[i] * pa_z[i];

        g_xxyyzzz_xxxyyz_0[i] = g_yyzzz_xxxyyz_0[i] * fbe_0 - g_yyzzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_xxyyz_1[i] * fe_0 + g_xyyzzz_xxxyyz_1[i] * pa_x[i];

        g_xxyyzzz_xxxyzz_0[i] = g_yyzzz_xxxyzz_0[i] * fbe_0 - g_yyzzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_xxyzz_1[i] * fe_0 + g_xyyzzz_xxxyzz_1[i] * pa_x[i];

        g_xxyyzzz_xxxzzz_0[i] = g_xxzzz_xxxzzz_0[i] * fbe_0 - g_xxzzz_xxxzzz_1[i] * fz_be_0 + g_xxyzzz_xxxzzz_1[i] * pa_y[i];

        g_xxyyzzz_xxyyyy_0[i] = 2.0 * g_xxyyz_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxyyyy_1[i] * fz_be_0 + g_xxyyzz_xxyyyy_1[i] * pa_z[i];

        g_xxyyzzz_xxyyyz_0[i] = g_yyzzz_xxyyyz_0[i] * fbe_0 - g_yyzzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_xyyyz_1[i] * fe_0 + g_xyyzzz_xxyyyz_1[i] * pa_x[i];

        g_xxyyzzz_xxyyzz_0[i] = g_yyzzz_xxyyzz_0[i] * fbe_0 - g_yyzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_xyyzz_1[i] * fe_0 + g_xyyzzz_xxyyzz_1[i] * pa_x[i];

        g_xxyyzzz_xxyzzz_0[i] = g_yyzzz_xxyzzz_0[i] * fbe_0 - g_yyzzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_xyzzz_1[i] * fe_0 + g_xyyzzz_xxyzzz_1[i] * pa_x[i];

        g_xxyyzzz_xxzzzz_0[i] = g_xxzzz_xxzzzz_0[i] * fbe_0 - g_xxzzz_xxzzzz_1[i] * fz_be_0 + g_xxyzzz_xxzzzz_1[i] * pa_y[i];

        g_xxyyzzz_xyyyyy_0[i] = 2.0 * g_xxyyz_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xyyyyy_1[i] * fz_be_0 + g_xxyyzz_xyyyyy_1[i] * pa_z[i];

        g_xxyyzzz_xyyyyz_0[i] = g_yyzzz_xyyyyz_0[i] * fbe_0 - g_yyzzz_xyyyyz_1[i] * fz_be_0 + g_xyyzzz_yyyyz_1[i] * fe_0 + g_xyyzzz_xyyyyz_1[i] * pa_x[i];

        g_xxyyzzz_xyyyzz_0[i] = g_yyzzz_xyyyzz_0[i] * fbe_0 - g_yyzzz_xyyyzz_1[i] * fz_be_0 + g_xyyzzz_yyyzz_1[i] * fe_0 + g_xyyzzz_xyyyzz_1[i] * pa_x[i];

        g_xxyyzzz_xyyzzz_0[i] = g_yyzzz_xyyzzz_0[i] * fbe_0 - g_yyzzz_xyyzzz_1[i] * fz_be_0 + g_xyyzzz_yyzzz_1[i] * fe_0 + g_xyyzzz_xyyzzz_1[i] * pa_x[i];

        g_xxyyzzz_xyzzzz_0[i] = g_yyzzz_xyzzzz_0[i] * fbe_0 - g_yyzzz_xyzzzz_1[i] * fz_be_0 + g_xyyzzz_yzzzz_1[i] * fe_0 + g_xyyzzz_xyzzzz_1[i] * pa_x[i];

        g_xxyyzzz_xzzzzz_0[i] = g_xxzzz_xzzzzz_0[i] * fbe_0 - g_xxzzz_xzzzzz_1[i] * fz_be_0 + g_xxyzzz_xzzzzz_1[i] * pa_y[i];

        g_xxyyzzz_yyyyyy_0[i] = g_yyzzz_yyyyyy_0[i] * fbe_0 - g_yyzzz_yyyyyy_1[i] * fz_be_0 + g_xyyzzz_yyyyyy_1[i] * pa_x[i];

        g_xxyyzzz_yyyyyz_0[i] = g_yyzzz_yyyyyz_0[i] * fbe_0 - g_yyzzz_yyyyyz_1[i] * fz_be_0 + g_xyyzzz_yyyyyz_1[i] * pa_x[i];

        g_xxyyzzz_yyyyzz_0[i] = g_yyzzz_yyyyzz_0[i] * fbe_0 - g_yyzzz_yyyyzz_1[i] * fz_be_0 + g_xyyzzz_yyyyzz_1[i] * pa_x[i];

        g_xxyyzzz_yyyzzz_0[i] = g_yyzzz_yyyzzz_0[i] * fbe_0 - g_yyzzz_yyyzzz_1[i] * fz_be_0 + g_xyyzzz_yyyzzz_1[i] * pa_x[i];

        g_xxyyzzz_yyzzzz_0[i] = g_yyzzz_yyzzzz_0[i] * fbe_0 - g_yyzzz_yyzzzz_1[i] * fz_be_0 + g_xyyzzz_yyzzzz_1[i] * pa_x[i];

        g_xxyyzzz_yzzzzz_0[i] = g_yyzzz_yzzzzz_0[i] * fbe_0 - g_yyzzz_yzzzzz_1[i] * fz_be_0 + g_xyyzzz_yzzzzz_1[i] * pa_x[i];

        g_xxyyzzz_zzzzzz_0[i] = g_yyzzz_zzzzzz_0[i] * fbe_0 - g_yyzzz_zzzzzz_1[i] * fz_be_0 + g_xyyzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 532-560 components of targeted buffer : KI

    auto g_xxyzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 532);

    auto g_xxyzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 533);

    auto g_xxyzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 534);

    auto g_xxyzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 535);

    auto g_xxyzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 536);

    auto g_xxyzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 537);

    auto g_xxyzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 538);

    auto g_xxyzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 539);

    auto g_xxyzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 540);

    auto g_xxyzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 541);

    auto g_xxyzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 542);

    auto g_xxyzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 543);

    auto g_xxyzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 544);

    auto g_xxyzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 545);

    auto g_xxyzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 546);

    auto g_xxyzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 547);

    auto g_xxyzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 548);

    auto g_xxyzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 549);

    auto g_xxyzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 550);

    auto g_xxyzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 551);

    auto g_xxyzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 552);

    auto g_xxyzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 553);

    auto g_xxyzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 554);

    auto g_xxyzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 555);

    auto g_xxyzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 556);

    auto g_xxyzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 557);

    auto g_xxyzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 558);

    auto g_xxyzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 559);

    #pragma omp simd aligned(g_xxyzzzz_xxxxxx_0, g_xxyzzzz_xxxxxy_0, g_xxyzzzz_xxxxxz_0, g_xxyzzzz_xxxxyy_0, g_xxyzzzz_xxxxyz_0, g_xxyzzzz_xxxxzz_0, g_xxyzzzz_xxxyyy_0, g_xxyzzzz_xxxyyz_0, g_xxyzzzz_xxxyzz_0, g_xxyzzzz_xxxzzz_0, g_xxyzzzz_xxyyyy_0, g_xxyzzzz_xxyyyz_0, g_xxyzzzz_xxyyzz_0, g_xxyzzzz_xxyzzz_0, g_xxyzzzz_xxzzzz_0, g_xxyzzzz_xyyyyy_0, g_xxyzzzz_xyyyyz_0, g_xxyzzzz_xyyyzz_0, g_xxyzzzz_xyyzzz_0, g_xxyzzzz_xyzzzz_0, g_xxyzzzz_xzzzzz_0, g_xxyzzzz_yyyyyy_0, g_xxyzzzz_yyyyyz_0, g_xxyzzzz_yyyyzz_0, g_xxyzzzz_yyyzzz_0, g_xxyzzzz_yyzzzz_0, g_xxyzzzz_yzzzzz_0, g_xxyzzzz_zzzzzz_0, g_xxzzzz_xxxxx_1, g_xxzzzz_xxxxxx_1, g_xxzzzz_xxxxxy_1, g_xxzzzz_xxxxxz_1, g_xxzzzz_xxxxy_1, g_xxzzzz_xxxxyy_1, g_xxzzzz_xxxxyz_1, g_xxzzzz_xxxxz_1, g_xxzzzz_xxxxzz_1, g_xxzzzz_xxxyy_1, g_xxzzzz_xxxyyy_1, g_xxzzzz_xxxyyz_1, g_xxzzzz_xxxyz_1, g_xxzzzz_xxxyzz_1, g_xxzzzz_xxxzz_1, g_xxzzzz_xxxzzz_1, g_xxzzzz_xxyyy_1, g_xxzzzz_xxyyyy_1, g_xxzzzz_xxyyyz_1, g_xxzzzz_xxyyz_1, g_xxzzzz_xxyyzz_1, g_xxzzzz_xxyzz_1, g_xxzzzz_xxyzzz_1, g_xxzzzz_xxzzz_1, g_xxzzzz_xxzzzz_1, g_xxzzzz_xyyyy_1, g_xxzzzz_xyyyyy_1, g_xxzzzz_xyyyyz_1, g_xxzzzz_xyyyz_1, g_xxzzzz_xyyyzz_1, g_xxzzzz_xyyzz_1, g_xxzzzz_xyyzzz_1, g_xxzzzz_xyzzz_1, g_xxzzzz_xyzzzz_1, g_xxzzzz_xzzzz_1, g_xxzzzz_xzzzzz_1, g_xxzzzz_yyyyy_1, g_xxzzzz_yyyyyy_1, g_xxzzzz_yyyyyz_1, g_xxzzzz_yyyyz_1, g_xxzzzz_yyyyzz_1, g_xxzzzz_yyyzz_1, g_xxzzzz_yyyzzz_1, g_xxzzzz_yyzzz_1, g_xxzzzz_yyzzzz_1, g_xxzzzz_yzzzz_1, g_xxzzzz_yzzzzz_1, g_xxzzzz_zzzzz_1, g_xxzzzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzzz_xxxxxx_0[i] = g_xxzzzz_xxxxxx_1[i] * pa_y[i];

        g_xxyzzzz_xxxxxy_0[i] = g_xxzzzz_xxxxx_1[i] * fe_0 + g_xxzzzz_xxxxxy_1[i] * pa_y[i];

        g_xxyzzzz_xxxxxz_0[i] = g_xxzzzz_xxxxxz_1[i] * pa_y[i];

        g_xxyzzzz_xxxxyy_0[i] = 2.0 * g_xxzzzz_xxxxy_1[i] * fe_0 + g_xxzzzz_xxxxyy_1[i] * pa_y[i];

        g_xxyzzzz_xxxxyz_0[i] = g_xxzzzz_xxxxz_1[i] * fe_0 + g_xxzzzz_xxxxyz_1[i] * pa_y[i];

        g_xxyzzzz_xxxxzz_0[i] = g_xxzzzz_xxxxzz_1[i] * pa_y[i];

        g_xxyzzzz_xxxyyy_0[i] = 3.0 * g_xxzzzz_xxxyy_1[i] * fe_0 + g_xxzzzz_xxxyyy_1[i] * pa_y[i];

        g_xxyzzzz_xxxyyz_0[i] = 2.0 * g_xxzzzz_xxxyz_1[i] * fe_0 + g_xxzzzz_xxxyyz_1[i] * pa_y[i];

        g_xxyzzzz_xxxyzz_0[i] = g_xxzzzz_xxxzz_1[i] * fe_0 + g_xxzzzz_xxxyzz_1[i] * pa_y[i];

        g_xxyzzzz_xxxzzz_0[i] = g_xxzzzz_xxxzzz_1[i] * pa_y[i];

        g_xxyzzzz_xxyyyy_0[i] = 4.0 * g_xxzzzz_xxyyy_1[i] * fe_0 + g_xxzzzz_xxyyyy_1[i] * pa_y[i];

        g_xxyzzzz_xxyyyz_0[i] = 3.0 * g_xxzzzz_xxyyz_1[i] * fe_0 + g_xxzzzz_xxyyyz_1[i] * pa_y[i];

        g_xxyzzzz_xxyyzz_0[i] = 2.0 * g_xxzzzz_xxyzz_1[i] * fe_0 + g_xxzzzz_xxyyzz_1[i] * pa_y[i];

        g_xxyzzzz_xxyzzz_0[i] = g_xxzzzz_xxzzz_1[i] * fe_0 + g_xxzzzz_xxyzzz_1[i] * pa_y[i];

        g_xxyzzzz_xxzzzz_0[i] = g_xxzzzz_xxzzzz_1[i] * pa_y[i];

        g_xxyzzzz_xyyyyy_0[i] = 5.0 * g_xxzzzz_xyyyy_1[i] * fe_0 + g_xxzzzz_xyyyyy_1[i] * pa_y[i];

        g_xxyzzzz_xyyyyz_0[i] = 4.0 * g_xxzzzz_xyyyz_1[i] * fe_0 + g_xxzzzz_xyyyyz_1[i] * pa_y[i];

        g_xxyzzzz_xyyyzz_0[i] = 3.0 * g_xxzzzz_xyyzz_1[i] * fe_0 + g_xxzzzz_xyyyzz_1[i] * pa_y[i];

        g_xxyzzzz_xyyzzz_0[i] = 2.0 * g_xxzzzz_xyzzz_1[i] * fe_0 + g_xxzzzz_xyyzzz_1[i] * pa_y[i];

        g_xxyzzzz_xyzzzz_0[i] = g_xxzzzz_xzzzz_1[i] * fe_0 + g_xxzzzz_xyzzzz_1[i] * pa_y[i];

        g_xxyzzzz_xzzzzz_0[i] = g_xxzzzz_xzzzzz_1[i] * pa_y[i];

        g_xxyzzzz_yyyyyy_0[i] = 6.0 * g_xxzzzz_yyyyy_1[i] * fe_0 + g_xxzzzz_yyyyyy_1[i] * pa_y[i];

        g_xxyzzzz_yyyyyz_0[i] = 5.0 * g_xxzzzz_yyyyz_1[i] * fe_0 + g_xxzzzz_yyyyyz_1[i] * pa_y[i];

        g_xxyzzzz_yyyyzz_0[i] = 4.0 * g_xxzzzz_yyyzz_1[i] * fe_0 + g_xxzzzz_yyyyzz_1[i] * pa_y[i];

        g_xxyzzzz_yyyzzz_0[i] = 3.0 * g_xxzzzz_yyzzz_1[i] * fe_0 + g_xxzzzz_yyyzzz_1[i] * pa_y[i];

        g_xxyzzzz_yyzzzz_0[i] = 2.0 * g_xxzzzz_yzzzz_1[i] * fe_0 + g_xxzzzz_yyzzzz_1[i] * pa_y[i];

        g_xxyzzzz_yzzzzz_0[i] = g_xxzzzz_zzzzz_1[i] * fe_0 + g_xxzzzz_yzzzzz_1[i] * pa_y[i];

        g_xxyzzzz_zzzzzz_0[i] = g_xxzzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 560-588 components of targeted buffer : KI

    auto g_xxzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 560);

    auto g_xxzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 561);

    auto g_xxzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 562);

    auto g_xxzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 563);

    auto g_xxzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 564);

    auto g_xxzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 565);

    auto g_xxzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 566);

    auto g_xxzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 567);

    auto g_xxzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 568);

    auto g_xxzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 569);

    auto g_xxzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 570);

    auto g_xxzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 571);

    auto g_xxzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 572);

    auto g_xxzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 573);

    auto g_xxzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 574);

    auto g_xxzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 575);

    auto g_xxzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 576);

    auto g_xxzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 577);

    auto g_xxzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 578);

    auto g_xxzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 579);

    auto g_xxzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 580);

    auto g_xxzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 581);

    auto g_xxzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 582);

    auto g_xxzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 583);

    auto g_xxzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 584);

    auto g_xxzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 585);

    auto g_xxzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 586);

    auto g_xxzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 587);

    #pragma omp simd aligned(g_xxzzz_xxxxxx_0, g_xxzzz_xxxxxx_1, g_xxzzz_xxxxxy_0, g_xxzzz_xxxxxy_1, g_xxzzz_xxxxyy_0, g_xxzzz_xxxxyy_1, g_xxzzz_xxxyyy_0, g_xxzzz_xxxyyy_1, g_xxzzz_xxyyyy_0, g_xxzzz_xxyyyy_1, g_xxzzz_xyyyyy_0, g_xxzzz_xyyyyy_1, g_xxzzzz_xxxxxx_1, g_xxzzzz_xxxxxy_1, g_xxzzzz_xxxxyy_1, g_xxzzzz_xxxyyy_1, g_xxzzzz_xxyyyy_1, g_xxzzzz_xyyyyy_1, g_xxzzzzz_xxxxxx_0, g_xxzzzzz_xxxxxy_0, g_xxzzzzz_xxxxxz_0, g_xxzzzzz_xxxxyy_0, g_xxzzzzz_xxxxyz_0, g_xxzzzzz_xxxxzz_0, g_xxzzzzz_xxxyyy_0, g_xxzzzzz_xxxyyz_0, g_xxzzzzz_xxxyzz_0, g_xxzzzzz_xxxzzz_0, g_xxzzzzz_xxyyyy_0, g_xxzzzzz_xxyyyz_0, g_xxzzzzz_xxyyzz_0, g_xxzzzzz_xxyzzz_0, g_xxzzzzz_xxzzzz_0, g_xxzzzzz_xyyyyy_0, g_xxzzzzz_xyyyyz_0, g_xxzzzzz_xyyyzz_0, g_xxzzzzz_xyyzzz_0, g_xxzzzzz_xyzzzz_0, g_xxzzzzz_xzzzzz_0, g_xxzzzzz_yyyyyy_0, g_xxzzzzz_yyyyyz_0, g_xxzzzzz_yyyyzz_0, g_xxzzzzz_yyyzzz_0, g_xxzzzzz_yyzzzz_0, g_xxzzzzz_yzzzzz_0, g_xxzzzzz_zzzzzz_0, g_xzzzzz_xxxxxz_1, g_xzzzzz_xxxxyz_1, g_xzzzzz_xxxxz_1, g_xzzzzz_xxxxzz_1, g_xzzzzz_xxxyyz_1, g_xzzzzz_xxxyz_1, g_xzzzzz_xxxyzz_1, g_xzzzzz_xxxzz_1, g_xzzzzz_xxxzzz_1, g_xzzzzz_xxyyyz_1, g_xzzzzz_xxyyz_1, g_xzzzzz_xxyyzz_1, g_xzzzzz_xxyzz_1, g_xzzzzz_xxyzzz_1, g_xzzzzz_xxzzz_1, g_xzzzzz_xxzzzz_1, g_xzzzzz_xyyyyz_1, g_xzzzzz_xyyyz_1, g_xzzzzz_xyyyzz_1, g_xzzzzz_xyyzz_1, g_xzzzzz_xyyzzz_1, g_xzzzzz_xyzzz_1, g_xzzzzz_xyzzzz_1, g_xzzzzz_xzzzz_1, g_xzzzzz_xzzzzz_1, g_xzzzzz_yyyyyy_1, g_xzzzzz_yyyyyz_1, g_xzzzzz_yyyyz_1, g_xzzzzz_yyyyzz_1, g_xzzzzz_yyyzz_1, g_xzzzzz_yyyzzz_1, g_xzzzzz_yyzzz_1, g_xzzzzz_yyzzzz_1, g_xzzzzz_yzzzz_1, g_xzzzzz_yzzzzz_1, g_xzzzzz_zzzzz_1, g_xzzzzz_zzzzzz_1, g_zzzzz_xxxxxz_0, g_zzzzz_xxxxxz_1, g_zzzzz_xxxxyz_0, g_zzzzz_xxxxyz_1, g_zzzzz_xxxxzz_0, g_zzzzz_xxxxzz_1, g_zzzzz_xxxyyz_0, g_zzzzz_xxxyyz_1, g_zzzzz_xxxyzz_0, g_zzzzz_xxxyzz_1, g_zzzzz_xxxzzz_0, g_zzzzz_xxxzzz_1, g_zzzzz_xxyyyz_0, g_zzzzz_xxyyyz_1, g_zzzzz_xxyyzz_0, g_zzzzz_xxyyzz_1, g_zzzzz_xxyzzz_0, g_zzzzz_xxyzzz_1, g_zzzzz_xxzzzz_0, g_zzzzz_xxzzzz_1, g_zzzzz_xyyyyz_0, g_zzzzz_xyyyyz_1, g_zzzzz_xyyyzz_0, g_zzzzz_xyyyzz_1, g_zzzzz_xyyzzz_0, g_zzzzz_xyyzzz_1, g_zzzzz_xyzzzz_0, g_zzzzz_xyzzzz_1, g_zzzzz_xzzzzz_0, g_zzzzz_xzzzzz_1, g_zzzzz_yyyyyy_0, g_zzzzz_yyyyyy_1, g_zzzzz_yyyyyz_0, g_zzzzz_yyyyyz_1, g_zzzzz_yyyyzz_0, g_zzzzz_yyyyzz_1, g_zzzzz_yyyzzz_0, g_zzzzz_yyyzzz_1, g_zzzzz_yyzzzz_0, g_zzzzz_yyzzzz_1, g_zzzzz_yzzzzz_0, g_zzzzz_yzzzzz_1, g_zzzzz_zzzzzz_0, g_zzzzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzzz_xxxxxx_0[i] = 4.0 * g_xxzzz_xxxxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxxxx_1[i] * fz_be_0 + g_xxzzzz_xxxxxx_1[i] * pa_z[i];

        g_xxzzzzz_xxxxxy_0[i] = 4.0 * g_xxzzz_xxxxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxxxy_1[i] * fz_be_0 + g_xxzzzz_xxxxxy_1[i] * pa_z[i];

        g_xxzzzzz_xxxxxz_0[i] = g_zzzzz_xxxxxz_0[i] * fbe_0 - g_zzzzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_xxxxz_1[i] * fe_0 + g_xzzzzz_xxxxxz_1[i] * pa_x[i];

        g_xxzzzzz_xxxxyy_0[i] = 4.0 * g_xxzzz_xxxxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxxyy_1[i] * fz_be_0 + g_xxzzzz_xxxxyy_1[i] * pa_z[i];

        g_xxzzzzz_xxxxyz_0[i] = g_zzzzz_xxxxyz_0[i] * fbe_0 - g_zzzzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_xxxyz_1[i] * fe_0 + g_xzzzzz_xxxxyz_1[i] * pa_x[i];

        g_xxzzzzz_xxxxzz_0[i] = g_zzzzz_xxxxzz_0[i] * fbe_0 - g_zzzzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_xxxzz_1[i] * fe_0 + g_xzzzzz_xxxxzz_1[i] * pa_x[i];

        g_xxzzzzz_xxxyyy_0[i] = 4.0 * g_xxzzz_xxxyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxyyy_1[i] * fz_be_0 + g_xxzzzz_xxxyyy_1[i] * pa_z[i];

        g_xxzzzzz_xxxyyz_0[i] = g_zzzzz_xxxyyz_0[i] * fbe_0 - g_zzzzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_xxyyz_1[i] * fe_0 + g_xzzzzz_xxxyyz_1[i] * pa_x[i];

        g_xxzzzzz_xxxyzz_0[i] = g_zzzzz_xxxyzz_0[i] * fbe_0 - g_zzzzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_xxyzz_1[i] * fe_0 + g_xzzzzz_xxxyzz_1[i] * pa_x[i];

        g_xxzzzzz_xxxzzz_0[i] = g_zzzzz_xxxzzz_0[i] * fbe_0 - g_zzzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_xxzzz_1[i] * fe_0 + g_xzzzzz_xxxzzz_1[i] * pa_x[i];

        g_xxzzzzz_xxyyyy_0[i] = 4.0 * g_xxzzz_xxyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxyyyy_1[i] * fz_be_0 + g_xxzzzz_xxyyyy_1[i] * pa_z[i];

        g_xxzzzzz_xxyyyz_0[i] = g_zzzzz_xxyyyz_0[i] * fbe_0 - g_zzzzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xyyyz_1[i] * fe_0 + g_xzzzzz_xxyyyz_1[i] * pa_x[i];

        g_xxzzzzz_xxyyzz_0[i] = g_zzzzz_xxyyzz_0[i] * fbe_0 - g_zzzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xyyzz_1[i] * fe_0 + g_xzzzzz_xxyyzz_1[i] * pa_x[i];

        g_xxzzzzz_xxyzzz_0[i] = g_zzzzz_xxyzzz_0[i] * fbe_0 - g_zzzzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xyzzz_1[i] * fe_0 + g_xzzzzz_xxyzzz_1[i] * pa_x[i];

        g_xxzzzzz_xxzzzz_0[i] = g_zzzzz_xxzzzz_0[i] * fbe_0 - g_zzzzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xzzzz_1[i] * fe_0 + g_xzzzzz_xxzzzz_1[i] * pa_x[i];

        g_xxzzzzz_xyyyyy_0[i] = 4.0 * g_xxzzz_xyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xyyyyy_1[i] * fz_be_0 + g_xxzzzz_xyyyyy_1[i] * pa_z[i];

        g_xxzzzzz_xyyyyz_0[i] = g_zzzzz_xyyyyz_0[i] * fbe_0 - g_zzzzz_xyyyyz_1[i] * fz_be_0 + g_xzzzzz_yyyyz_1[i] * fe_0 + g_xzzzzz_xyyyyz_1[i] * pa_x[i];

        g_xxzzzzz_xyyyzz_0[i] = g_zzzzz_xyyyzz_0[i] * fbe_0 - g_zzzzz_xyyyzz_1[i] * fz_be_0 + g_xzzzzz_yyyzz_1[i] * fe_0 + g_xzzzzz_xyyyzz_1[i] * pa_x[i];

        g_xxzzzzz_xyyzzz_0[i] = g_zzzzz_xyyzzz_0[i] * fbe_0 - g_zzzzz_xyyzzz_1[i] * fz_be_0 + g_xzzzzz_yyzzz_1[i] * fe_0 + g_xzzzzz_xyyzzz_1[i] * pa_x[i];

        g_xxzzzzz_xyzzzz_0[i] = g_zzzzz_xyzzzz_0[i] * fbe_0 - g_zzzzz_xyzzzz_1[i] * fz_be_0 + g_xzzzzz_yzzzz_1[i] * fe_0 + g_xzzzzz_xyzzzz_1[i] * pa_x[i];

        g_xxzzzzz_xzzzzz_0[i] = g_zzzzz_xzzzzz_0[i] * fbe_0 - g_zzzzz_xzzzzz_1[i] * fz_be_0 + g_xzzzzz_zzzzz_1[i] * fe_0 + g_xzzzzz_xzzzzz_1[i] * pa_x[i];

        g_xxzzzzz_yyyyyy_0[i] = g_zzzzz_yyyyyy_0[i] * fbe_0 - g_zzzzz_yyyyyy_1[i] * fz_be_0 + g_xzzzzz_yyyyyy_1[i] * pa_x[i];

        g_xxzzzzz_yyyyyz_0[i] = g_zzzzz_yyyyyz_0[i] * fbe_0 - g_zzzzz_yyyyyz_1[i] * fz_be_0 + g_xzzzzz_yyyyyz_1[i] * pa_x[i];

        g_xxzzzzz_yyyyzz_0[i] = g_zzzzz_yyyyzz_0[i] * fbe_0 - g_zzzzz_yyyyzz_1[i] * fz_be_0 + g_xzzzzz_yyyyzz_1[i] * pa_x[i];

        g_xxzzzzz_yyyzzz_0[i] = g_zzzzz_yyyzzz_0[i] * fbe_0 - g_zzzzz_yyyzzz_1[i] * fz_be_0 + g_xzzzzz_yyyzzz_1[i] * pa_x[i];

        g_xxzzzzz_yyzzzz_0[i] = g_zzzzz_yyzzzz_0[i] * fbe_0 - g_zzzzz_yyzzzz_1[i] * fz_be_0 + g_xzzzzz_yyzzzz_1[i] * pa_x[i];

        g_xxzzzzz_yzzzzz_0[i] = g_zzzzz_yzzzzz_0[i] * fbe_0 - g_zzzzz_yzzzzz_1[i] * fz_be_0 + g_xzzzzz_yzzzzz_1[i] * pa_x[i];

        g_xxzzzzz_zzzzzz_0[i] = g_zzzzz_zzzzzz_0[i] * fbe_0 - g_zzzzz_zzzzzz_1[i] * fz_be_0 + g_xzzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 588-616 components of targeted buffer : KI

    auto g_xyyyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 588);

    auto g_xyyyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 589);

    auto g_xyyyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 590);

    auto g_xyyyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 591);

    auto g_xyyyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 592);

    auto g_xyyyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 593);

    auto g_xyyyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 594);

    auto g_xyyyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 595);

    auto g_xyyyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 596);

    auto g_xyyyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 597);

    auto g_xyyyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 598);

    auto g_xyyyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 599);

    auto g_xyyyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 600);

    auto g_xyyyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 601);

    auto g_xyyyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 602);

    auto g_xyyyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 603);

    auto g_xyyyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 604);

    auto g_xyyyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 605);

    auto g_xyyyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 606);

    auto g_xyyyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 607);

    auto g_xyyyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 608);

    auto g_xyyyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 609);

    auto g_xyyyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 610);

    auto g_xyyyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 611);

    auto g_xyyyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 612);

    auto g_xyyyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 613);

    auto g_xyyyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 614);

    auto g_xyyyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 615);

    #pragma omp simd aligned(g_xyyyyyy_xxxxxx_0, g_xyyyyyy_xxxxxy_0, g_xyyyyyy_xxxxxz_0, g_xyyyyyy_xxxxyy_0, g_xyyyyyy_xxxxyz_0, g_xyyyyyy_xxxxzz_0, g_xyyyyyy_xxxyyy_0, g_xyyyyyy_xxxyyz_0, g_xyyyyyy_xxxyzz_0, g_xyyyyyy_xxxzzz_0, g_xyyyyyy_xxyyyy_0, g_xyyyyyy_xxyyyz_0, g_xyyyyyy_xxyyzz_0, g_xyyyyyy_xxyzzz_0, g_xyyyyyy_xxzzzz_0, g_xyyyyyy_xyyyyy_0, g_xyyyyyy_xyyyyz_0, g_xyyyyyy_xyyyzz_0, g_xyyyyyy_xyyzzz_0, g_xyyyyyy_xyzzzz_0, g_xyyyyyy_xzzzzz_0, g_xyyyyyy_yyyyyy_0, g_xyyyyyy_yyyyyz_0, g_xyyyyyy_yyyyzz_0, g_xyyyyyy_yyyzzz_0, g_xyyyyyy_yyzzzz_0, g_xyyyyyy_yzzzzz_0, g_xyyyyyy_zzzzzz_0, g_yyyyyy_xxxxx_1, g_yyyyyy_xxxxxx_1, g_yyyyyy_xxxxxy_1, g_yyyyyy_xxxxxz_1, g_yyyyyy_xxxxy_1, g_yyyyyy_xxxxyy_1, g_yyyyyy_xxxxyz_1, g_yyyyyy_xxxxz_1, g_yyyyyy_xxxxzz_1, g_yyyyyy_xxxyy_1, g_yyyyyy_xxxyyy_1, g_yyyyyy_xxxyyz_1, g_yyyyyy_xxxyz_1, g_yyyyyy_xxxyzz_1, g_yyyyyy_xxxzz_1, g_yyyyyy_xxxzzz_1, g_yyyyyy_xxyyy_1, g_yyyyyy_xxyyyy_1, g_yyyyyy_xxyyyz_1, g_yyyyyy_xxyyz_1, g_yyyyyy_xxyyzz_1, g_yyyyyy_xxyzz_1, g_yyyyyy_xxyzzz_1, g_yyyyyy_xxzzz_1, g_yyyyyy_xxzzzz_1, g_yyyyyy_xyyyy_1, g_yyyyyy_xyyyyy_1, g_yyyyyy_xyyyyz_1, g_yyyyyy_xyyyz_1, g_yyyyyy_xyyyzz_1, g_yyyyyy_xyyzz_1, g_yyyyyy_xyyzzz_1, g_yyyyyy_xyzzz_1, g_yyyyyy_xyzzzz_1, g_yyyyyy_xzzzz_1, g_yyyyyy_xzzzzz_1, g_yyyyyy_yyyyy_1, g_yyyyyy_yyyyyy_1, g_yyyyyy_yyyyyz_1, g_yyyyyy_yyyyz_1, g_yyyyyy_yyyyzz_1, g_yyyyyy_yyyzz_1, g_yyyyyy_yyyzzz_1, g_yyyyyy_yyzzz_1, g_yyyyyy_yyzzzz_1, g_yyyyyy_yzzzz_1, g_yyyyyy_yzzzzz_1, g_yyyyyy_zzzzz_1, g_yyyyyy_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyy_xxxxxx_0[i] = 6.0 * g_yyyyyy_xxxxx_1[i] * fe_0 + g_yyyyyy_xxxxxx_1[i] * pa_x[i];

        g_xyyyyyy_xxxxxy_0[i] = 5.0 * g_yyyyyy_xxxxy_1[i] * fe_0 + g_yyyyyy_xxxxxy_1[i] * pa_x[i];

        g_xyyyyyy_xxxxxz_0[i] = 5.0 * g_yyyyyy_xxxxz_1[i] * fe_0 + g_yyyyyy_xxxxxz_1[i] * pa_x[i];

        g_xyyyyyy_xxxxyy_0[i] = 4.0 * g_yyyyyy_xxxyy_1[i] * fe_0 + g_yyyyyy_xxxxyy_1[i] * pa_x[i];

        g_xyyyyyy_xxxxyz_0[i] = 4.0 * g_yyyyyy_xxxyz_1[i] * fe_0 + g_yyyyyy_xxxxyz_1[i] * pa_x[i];

        g_xyyyyyy_xxxxzz_0[i] = 4.0 * g_yyyyyy_xxxzz_1[i] * fe_0 + g_yyyyyy_xxxxzz_1[i] * pa_x[i];

        g_xyyyyyy_xxxyyy_0[i] = 3.0 * g_yyyyyy_xxyyy_1[i] * fe_0 + g_yyyyyy_xxxyyy_1[i] * pa_x[i];

        g_xyyyyyy_xxxyyz_0[i] = 3.0 * g_yyyyyy_xxyyz_1[i] * fe_0 + g_yyyyyy_xxxyyz_1[i] * pa_x[i];

        g_xyyyyyy_xxxyzz_0[i] = 3.0 * g_yyyyyy_xxyzz_1[i] * fe_0 + g_yyyyyy_xxxyzz_1[i] * pa_x[i];

        g_xyyyyyy_xxxzzz_0[i] = 3.0 * g_yyyyyy_xxzzz_1[i] * fe_0 + g_yyyyyy_xxxzzz_1[i] * pa_x[i];

        g_xyyyyyy_xxyyyy_0[i] = 2.0 * g_yyyyyy_xyyyy_1[i] * fe_0 + g_yyyyyy_xxyyyy_1[i] * pa_x[i];

        g_xyyyyyy_xxyyyz_0[i] = 2.0 * g_yyyyyy_xyyyz_1[i] * fe_0 + g_yyyyyy_xxyyyz_1[i] * pa_x[i];

        g_xyyyyyy_xxyyzz_0[i] = 2.0 * g_yyyyyy_xyyzz_1[i] * fe_0 + g_yyyyyy_xxyyzz_1[i] * pa_x[i];

        g_xyyyyyy_xxyzzz_0[i] = 2.0 * g_yyyyyy_xyzzz_1[i] * fe_0 + g_yyyyyy_xxyzzz_1[i] * pa_x[i];

        g_xyyyyyy_xxzzzz_0[i] = 2.0 * g_yyyyyy_xzzzz_1[i] * fe_0 + g_yyyyyy_xxzzzz_1[i] * pa_x[i];

        g_xyyyyyy_xyyyyy_0[i] = g_yyyyyy_yyyyy_1[i] * fe_0 + g_yyyyyy_xyyyyy_1[i] * pa_x[i];

        g_xyyyyyy_xyyyyz_0[i] = g_yyyyyy_yyyyz_1[i] * fe_0 + g_yyyyyy_xyyyyz_1[i] * pa_x[i];

        g_xyyyyyy_xyyyzz_0[i] = g_yyyyyy_yyyzz_1[i] * fe_0 + g_yyyyyy_xyyyzz_1[i] * pa_x[i];

        g_xyyyyyy_xyyzzz_0[i] = g_yyyyyy_yyzzz_1[i] * fe_0 + g_yyyyyy_xyyzzz_1[i] * pa_x[i];

        g_xyyyyyy_xyzzzz_0[i] = g_yyyyyy_yzzzz_1[i] * fe_0 + g_yyyyyy_xyzzzz_1[i] * pa_x[i];

        g_xyyyyyy_xzzzzz_0[i] = g_yyyyyy_zzzzz_1[i] * fe_0 + g_yyyyyy_xzzzzz_1[i] * pa_x[i];

        g_xyyyyyy_yyyyyy_0[i] = g_yyyyyy_yyyyyy_1[i] * pa_x[i];

        g_xyyyyyy_yyyyyz_0[i] = g_yyyyyy_yyyyyz_1[i] * pa_x[i];

        g_xyyyyyy_yyyyzz_0[i] = g_yyyyyy_yyyyzz_1[i] * pa_x[i];

        g_xyyyyyy_yyyzzz_0[i] = g_yyyyyy_yyyzzz_1[i] * pa_x[i];

        g_xyyyyyy_yyzzzz_0[i] = g_yyyyyy_yyzzzz_1[i] * pa_x[i];

        g_xyyyyyy_yzzzzz_0[i] = g_yyyyyy_yzzzzz_1[i] * pa_x[i];

        g_xyyyyyy_zzzzzz_0[i] = g_yyyyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 616-644 components of targeted buffer : KI

    auto g_xyyyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 616);

    auto g_xyyyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 617);

    auto g_xyyyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 618);

    auto g_xyyyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 619);

    auto g_xyyyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 620);

    auto g_xyyyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 621);

    auto g_xyyyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 622);

    auto g_xyyyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 623);

    auto g_xyyyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 624);

    auto g_xyyyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 625);

    auto g_xyyyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 626);

    auto g_xyyyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 627);

    auto g_xyyyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 628);

    auto g_xyyyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 629);

    auto g_xyyyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 630);

    auto g_xyyyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 631);

    auto g_xyyyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 632);

    auto g_xyyyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 633);

    auto g_xyyyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 634);

    auto g_xyyyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 635);

    auto g_xyyyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 636);

    auto g_xyyyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 637);

    auto g_xyyyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 638);

    auto g_xyyyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 639);

    auto g_xyyyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 640);

    auto g_xyyyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 641);

    auto g_xyyyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 642);

    auto g_xyyyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 643);

    #pragma omp simd aligned(g_xyyyyy_xxxxxx_1, g_xyyyyy_xxxxxy_1, g_xyyyyy_xxxxyy_1, g_xyyyyy_xxxyyy_1, g_xyyyyy_xxyyyy_1, g_xyyyyy_xyyyyy_1, g_xyyyyyz_xxxxxx_0, g_xyyyyyz_xxxxxy_0, g_xyyyyyz_xxxxxz_0, g_xyyyyyz_xxxxyy_0, g_xyyyyyz_xxxxyz_0, g_xyyyyyz_xxxxzz_0, g_xyyyyyz_xxxyyy_0, g_xyyyyyz_xxxyyz_0, g_xyyyyyz_xxxyzz_0, g_xyyyyyz_xxxzzz_0, g_xyyyyyz_xxyyyy_0, g_xyyyyyz_xxyyyz_0, g_xyyyyyz_xxyyzz_0, g_xyyyyyz_xxyzzz_0, g_xyyyyyz_xxzzzz_0, g_xyyyyyz_xyyyyy_0, g_xyyyyyz_xyyyyz_0, g_xyyyyyz_xyyyzz_0, g_xyyyyyz_xyyzzz_0, g_xyyyyyz_xyzzzz_0, g_xyyyyyz_xzzzzz_0, g_xyyyyyz_yyyyyy_0, g_xyyyyyz_yyyyyz_0, g_xyyyyyz_yyyyzz_0, g_xyyyyyz_yyyzzz_0, g_xyyyyyz_yyzzzz_0, g_xyyyyyz_yzzzzz_0, g_xyyyyyz_zzzzzz_0, g_yyyyyz_xxxxxz_1, g_yyyyyz_xxxxyz_1, g_yyyyyz_xxxxz_1, g_yyyyyz_xxxxzz_1, g_yyyyyz_xxxyyz_1, g_yyyyyz_xxxyz_1, g_yyyyyz_xxxyzz_1, g_yyyyyz_xxxzz_1, g_yyyyyz_xxxzzz_1, g_yyyyyz_xxyyyz_1, g_yyyyyz_xxyyz_1, g_yyyyyz_xxyyzz_1, g_yyyyyz_xxyzz_1, g_yyyyyz_xxyzzz_1, g_yyyyyz_xxzzz_1, g_yyyyyz_xxzzzz_1, g_yyyyyz_xyyyyz_1, g_yyyyyz_xyyyz_1, g_yyyyyz_xyyyzz_1, g_yyyyyz_xyyzz_1, g_yyyyyz_xyyzzz_1, g_yyyyyz_xyzzz_1, g_yyyyyz_xyzzzz_1, g_yyyyyz_xzzzz_1, g_yyyyyz_xzzzzz_1, g_yyyyyz_yyyyyy_1, g_yyyyyz_yyyyyz_1, g_yyyyyz_yyyyz_1, g_yyyyyz_yyyyzz_1, g_yyyyyz_yyyzz_1, g_yyyyyz_yyyzzz_1, g_yyyyyz_yyzzz_1, g_yyyyyz_yyzzzz_1, g_yyyyyz_yzzzz_1, g_yyyyyz_yzzzzz_1, g_yyyyyz_zzzzz_1, g_yyyyyz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyz_xxxxxx_0[i] = g_xyyyyy_xxxxxx_1[i] * pa_z[i];

        g_xyyyyyz_xxxxxy_0[i] = g_xyyyyy_xxxxxy_1[i] * pa_z[i];

        g_xyyyyyz_xxxxxz_0[i] = 5.0 * g_yyyyyz_xxxxz_1[i] * fe_0 + g_yyyyyz_xxxxxz_1[i] * pa_x[i];

        g_xyyyyyz_xxxxyy_0[i] = g_xyyyyy_xxxxyy_1[i] * pa_z[i];

        g_xyyyyyz_xxxxyz_0[i] = 4.0 * g_yyyyyz_xxxyz_1[i] * fe_0 + g_yyyyyz_xxxxyz_1[i] * pa_x[i];

        g_xyyyyyz_xxxxzz_0[i] = 4.0 * g_yyyyyz_xxxzz_1[i] * fe_0 + g_yyyyyz_xxxxzz_1[i] * pa_x[i];

        g_xyyyyyz_xxxyyy_0[i] = g_xyyyyy_xxxyyy_1[i] * pa_z[i];

        g_xyyyyyz_xxxyyz_0[i] = 3.0 * g_yyyyyz_xxyyz_1[i] * fe_0 + g_yyyyyz_xxxyyz_1[i] * pa_x[i];

        g_xyyyyyz_xxxyzz_0[i] = 3.0 * g_yyyyyz_xxyzz_1[i] * fe_0 + g_yyyyyz_xxxyzz_1[i] * pa_x[i];

        g_xyyyyyz_xxxzzz_0[i] = 3.0 * g_yyyyyz_xxzzz_1[i] * fe_0 + g_yyyyyz_xxxzzz_1[i] * pa_x[i];

        g_xyyyyyz_xxyyyy_0[i] = g_xyyyyy_xxyyyy_1[i] * pa_z[i];

        g_xyyyyyz_xxyyyz_0[i] = 2.0 * g_yyyyyz_xyyyz_1[i] * fe_0 + g_yyyyyz_xxyyyz_1[i] * pa_x[i];

        g_xyyyyyz_xxyyzz_0[i] = 2.0 * g_yyyyyz_xyyzz_1[i] * fe_0 + g_yyyyyz_xxyyzz_1[i] * pa_x[i];

        g_xyyyyyz_xxyzzz_0[i] = 2.0 * g_yyyyyz_xyzzz_1[i] * fe_0 + g_yyyyyz_xxyzzz_1[i] * pa_x[i];

        g_xyyyyyz_xxzzzz_0[i] = 2.0 * g_yyyyyz_xzzzz_1[i] * fe_0 + g_yyyyyz_xxzzzz_1[i] * pa_x[i];

        g_xyyyyyz_xyyyyy_0[i] = g_xyyyyy_xyyyyy_1[i] * pa_z[i];

        g_xyyyyyz_xyyyyz_0[i] = g_yyyyyz_yyyyz_1[i] * fe_0 + g_yyyyyz_xyyyyz_1[i] * pa_x[i];

        g_xyyyyyz_xyyyzz_0[i] = g_yyyyyz_yyyzz_1[i] * fe_0 + g_yyyyyz_xyyyzz_1[i] * pa_x[i];

        g_xyyyyyz_xyyzzz_0[i] = g_yyyyyz_yyzzz_1[i] * fe_0 + g_yyyyyz_xyyzzz_1[i] * pa_x[i];

        g_xyyyyyz_xyzzzz_0[i] = g_yyyyyz_yzzzz_1[i] * fe_0 + g_yyyyyz_xyzzzz_1[i] * pa_x[i];

        g_xyyyyyz_xzzzzz_0[i] = g_yyyyyz_zzzzz_1[i] * fe_0 + g_yyyyyz_xzzzzz_1[i] * pa_x[i];

        g_xyyyyyz_yyyyyy_0[i] = g_yyyyyz_yyyyyy_1[i] * pa_x[i];

        g_xyyyyyz_yyyyyz_0[i] = g_yyyyyz_yyyyyz_1[i] * pa_x[i];

        g_xyyyyyz_yyyyzz_0[i] = g_yyyyyz_yyyyzz_1[i] * pa_x[i];

        g_xyyyyyz_yyyzzz_0[i] = g_yyyyyz_yyyzzz_1[i] * pa_x[i];

        g_xyyyyyz_yyzzzz_0[i] = g_yyyyyz_yyzzzz_1[i] * pa_x[i];

        g_xyyyyyz_yzzzzz_0[i] = g_yyyyyz_yzzzzz_1[i] * pa_x[i];

        g_xyyyyyz_zzzzzz_0[i] = g_yyyyyz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 644-672 components of targeted buffer : KI

    auto g_xyyyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 644);

    auto g_xyyyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 645);

    auto g_xyyyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 646);

    auto g_xyyyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 647);

    auto g_xyyyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 648);

    auto g_xyyyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 649);

    auto g_xyyyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 650);

    auto g_xyyyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 651);

    auto g_xyyyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 652);

    auto g_xyyyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 653);

    auto g_xyyyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 654);

    auto g_xyyyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 655);

    auto g_xyyyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 656);

    auto g_xyyyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 657);

    auto g_xyyyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 658);

    auto g_xyyyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 659);

    auto g_xyyyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 660);

    auto g_xyyyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 661);

    auto g_xyyyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 662);

    auto g_xyyyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 663);

    auto g_xyyyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 664);

    auto g_xyyyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 665);

    auto g_xyyyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 666);

    auto g_xyyyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 667);

    auto g_xyyyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 668);

    auto g_xyyyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 669);

    auto g_xyyyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 670);

    auto g_xyyyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 671);

    #pragma omp simd aligned(g_xyyyyzz_xxxxxx_0, g_xyyyyzz_xxxxxy_0, g_xyyyyzz_xxxxxz_0, g_xyyyyzz_xxxxyy_0, g_xyyyyzz_xxxxyz_0, g_xyyyyzz_xxxxzz_0, g_xyyyyzz_xxxyyy_0, g_xyyyyzz_xxxyyz_0, g_xyyyyzz_xxxyzz_0, g_xyyyyzz_xxxzzz_0, g_xyyyyzz_xxyyyy_0, g_xyyyyzz_xxyyyz_0, g_xyyyyzz_xxyyzz_0, g_xyyyyzz_xxyzzz_0, g_xyyyyzz_xxzzzz_0, g_xyyyyzz_xyyyyy_0, g_xyyyyzz_xyyyyz_0, g_xyyyyzz_xyyyzz_0, g_xyyyyzz_xyyzzz_0, g_xyyyyzz_xyzzzz_0, g_xyyyyzz_xzzzzz_0, g_xyyyyzz_yyyyyy_0, g_xyyyyzz_yyyyyz_0, g_xyyyyzz_yyyyzz_0, g_xyyyyzz_yyyzzz_0, g_xyyyyzz_yyzzzz_0, g_xyyyyzz_yzzzzz_0, g_xyyyyzz_zzzzzz_0, g_yyyyzz_xxxxx_1, g_yyyyzz_xxxxxx_1, g_yyyyzz_xxxxxy_1, g_yyyyzz_xxxxxz_1, g_yyyyzz_xxxxy_1, g_yyyyzz_xxxxyy_1, g_yyyyzz_xxxxyz_1, g_yyyyzz_xxxxz_1, g_yyyyzz_xxxxzz_1, g_yyyyzz_xxxyy_1, g_yyyyzz_xxxyyy_1, g_yyyyzz_xxxyyz_1, g_yyyyzz_xxxyz_1, g_yyyyzz_xxxyzz_1, g_yyyyzz_xxxzz_1, g_yyyyzz_xxxzzz_1, g_yyyyzz_xxyyy_1, g_yyyyzz_xxyyyy_1, g_yyyyzz_xxyyyz_1, g_yyyyzz_xxyyz_1, g_yyyyzz_xxyyzz_1, g_yyyyzz_xxyzz_1, g_yyyyzz_xxyzzz_1, g_yyyyzz_xxzzz_1, g_yyyyzz_xxzzzz_1, g_yyyyzz_xyyyy_1, g_yyyyzz_xyyyyy_1, g_yyyyzz_xyyyyz_1, g_yyyyzz_xyyyz_1, g_yyyyzz_xyyyzz_1, g_yyyyzz_xyyzz_1, g_yyyyzz_xyyzzz_1, g_yyyyzz_xyzzz_1, g_yyyyzz_xyzzzz_1, g_yyyyzz_xzzzz_1, g_yyyyzz_xzzzzz_1, g_yyyyzz_yyyyy_1, g_yyyyzz_yyyyyy_1, g_yyyyzz_yyyyyz_1, g_yyyyzz_yyyyz_1, g_yyyyzz_yyyyzz_1, g_yyyyzz_yyyzz_1, g_yyyyzz_yyyzzz_1, g_yyyyzz_yyzzz_1, g_yyyyzz_yyzzzz_1, g_yyyyzz_yzzzz_1, g_yyyyzz_yzzzzz_1, g_yyyyzz_zzzzz_1, g_yyyyzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyzz_xxxxxx_0[i] = 6.0 * g_yyyyzz_xxxxx_1[i] * fe_0 + g_yyyyzz_xxxxxx_1[i] * pa_x[i];

        g_xyyyyzz_xxxxxy_0[i] = 5.0 * g_yyyyzz_xxxxy_1[i] * fe_0 + g_yyyyzz_xxxxxy_1[i] * pa_x[i];

        g_xyyyyzz_xxxxxz_0[i] = 5.0 * g_yyyyzz_xxxxz_1[i] * fe_0 + g_yyyyzz_xxxxxz_1[i] * pa_x[i];

        g_xyyyyzz_xxxxyy_0[i] = 4.0 * g_yyyyzz_xxxyy_1[i] * fe_0 + g_yyyyzz_xxxxyy_1[i] * pa_x[i];

        g_xyyyyzz_xxxxyz_0[i] = 4.0 * g_yyyyzz_xxxyz_1[i] * fe_0 + g_yyyyzz_xxxxyz_1[i] * pa_x[i];

        g_xyyyyzz_xxxxzz_0[i] = 4.0 * g_yyyyzz_xxxzz_1[i] * fe_0 + g_yyyyzz_xxxxzz_1[i] * pa_x[i];

        g_xyyyyzz_xxxyyy_0[i] = 3.0 * g_yyyyzz_xxyyy_1[i] * fe_0 + g_yyyyzz_xxxyyy_1[i] * pa_x[i];

        g_xyyyyzz_xxxyyz_0[i] = 3.0 * g_yyyyzz_xxyyz_1[i] * fe_0 + g_yyyyzz_xxxyyz_1[i] * pa_x[i];

        g_xyyyyzz_xxxyzz_0[i] = 3.0 * g_yyyyzz_xxyzz_1[i] * fe_0 + g_yyyyzz_xxxyzz_1[i] * pa_x[i];

        g_xyyyyzz_xxxzzz_0[i] = 3.0 * g_yyyyzz_xxzzz_1[i] * fe_0 + g_yyyyzz_xxxzzz_1[i] * pa_x[i];

        g_xyyyyzz_xxyyyy_0[i] = 2.0 * g_yyyyzz_xyyyy_1[i] * fe_0 + g_yyyyzz_xxyyyy_1[i] * pa_x[i];

        g_xyyyyzz_xxyyyz_0[i] = 2.0 * g_yyyyzz_xyyyz_1[i] * fe_0 + g_yyyyzz_xxyyyz_1[i] * pa_x[i];

        g_xyyyyzz_xxyyzz_0[i] = 2.0 * g_yyyyzz_xyyzz_1[i] * fe_0 + g_yyyyzz_xxyyzz_1[i] * pa_x[i];

        g_xyyyyzz_xxyzzz_0[i] = 2.0 * g_yyyyzz_xyzzz_1[i] * fe_0 + g_yyyyzz_xxyzzz_1[i] * pa_x[i];

        g_xyyyyzz_xxzzzz_0[i] = 2.0 * g_yyyyzz_xzzzz_1[i] * fe_0 + g_yyyyzz_xxzzzz_1[i] * pa_x[i];

        g_xyyyyzz_xyyyyy_0[i] = g_yyyyzz_yyyyy_1[i] * fe_0 + g_yyyyzz_xyyyyy_1[i] * pa_x[i];

        g_xyyyyzz_xyyyyz_0[i] = g_yyyyzz_yyyyz_1[i] * fe_0 + g_yyyyzz_xyyyyz_1[i] * pa_x[i];

        g_xyyyyzz_xyyyzz_0[i] = g_yyyyzz_yyyzz_1[i] * fe_0 + g_yyyyzz_xyyyzz_1[i] * pa_x[i];

        g_xyyyyzz_xyyzzz_0[i] = g_yyyyzz_yyzzz_1[i] * fe_0 + g_yyyyzz_xyyzzz_1[i] * pa_x[i];

        g_xyyyyzz_xyzzzz_0[i] = g_yyyyzz_yzzzz_1[i] * fe_0 + g_yyyyzz_xyzzzz_1[i] * pa_x[i];

        g_xyyyyzz_xzzzzz_0[i] = g_yyyyzz_zzzzz_1[i] * fe_0 + g_yyyyzz_xzzzzz_1[i] * pa_x[i];

        g_xyyyyzz_yyyyyy_0[i] = g_yyyyzz_yyyyyy_1[i] * pa_x[i];

        g_xyyyyzz_yyyyyz_0[i] = g_yyyyzz_yyyyyz_1[i] * pa_x[i];

        g_xyyyyzz_yyyyzz_0[i] = g_yyyyzz_yyyyzz_1[i] * pa_x[i];

        g_xyyyyzz_yyyzzz_0[i] = g_yyyyzz_yyyzzz_1[i] * pa_x[i];

        g_xyyyyzz_yyzzzz_0[i] = g_yyyyzz_yyzzzz_1[i] * pa_x[i];

        g_xyyyyzz_yzzzzz_0[i] = g_yyyyzz_yzzzzz_1[i] * pa_x[i];

        g_xyyyyzz_zzzzzz_0[i] = g_yyyyzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 672-700 components of targeted buffer : KI

    auto g_xyyyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 672);

    auto g_xyyyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 673);

    auto g_xyyyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 674);

    auto g_xyyyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 675);

    auto g_xyyyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 676);

    auto g_xyyyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 677);

    auto g_xyyyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 678);

    auto g_xyyyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 679);

    auto g_xyyyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 680);

    auto g_xyyyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 681);

    auto g_xyyyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 682);

    auto g_xyyyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 683);

    auto g_xyyyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 684);

    auto g_xyyyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 685);

    auto g_xyyyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 686);

    auto g_xyyyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 687);

    auto g_xyyyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 688);

    auto g_xyyyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 689);

    auto g_xyyyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 690);

    auto g_xyyyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 691);

    auto g_xyyyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 692);

    auto g_xyyyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 693);

    auto g_xyyyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 694);

    auto g_xyyyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 695);

    auto g_xyyyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 696);

    auto g_xyyyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 697);

    auto g_xyyyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 698);

    auto g_xyyyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 699);

    #pragma omp simd aligned(g_xyyyzzz_xxxxxx_0, g_xyyyzzz_xxxxxy_0, g_xyyyzzz_xxxxxz_0, g_xyyyzzz_xxxxyy_0, g_xyyyzzz_xxxxyz_0, g_xyyyzzz_xxxxzz_0, g_xyyyzzz_xxxyyy_0, g_xyyyzzz_xxxyyz_0, g_xyyyzzz_xxxyzz_0, g_xyyyzzz_xxxzzz_0, g_xyyyzzz_xxyyyy_0, g_xyyyzzz_xxyyyz_0, g_xyyyzzz_xxyyzz_0, g_xyyyzzz_xxyzzz_0, g_xyyyzzz_xxzzzz_0, g_xyyyzzz_xyyyyy_0, g_xyyyzzz_xyyyyz_0, g_xyyyzzz_xyyyzz_0, g_xyyyzzz_xyyzzz_0, g_xyyyzzz_xyzzzz_0, g_xyyyzzz_xzzzzz_0, g_xyyyzzz_yyyyyy_0, g_xyyyzzz_yyyyyz_0, g_xyyyzzz_yyyyzz_0, g_xyyyzzz_yyyzzz_0, g_xyyyzzz_yyzzzz_0, g_xyyyzzz_yzzzzz_0, g_xyyyzzz_zzzzzz_0, g_yyyzzz_xxxxx_1, g_yyyzzz_xxxxxx_1, g_yyyzzz_xxxxxy_1, g_yyyzzz_xxxxxz_1, g_yyyzzz_xxxxy_1, g_yyyzzz_xxxxyy_1, g_yyyzzz_xxxxyz_1, g_yyyzzz_xxxxz_1, g_yyyzzz_xxxxzz_1, g_yyyzzz_xxxyy_1, g_yyyzzz_xxxyyy_1, g_yyyzzz_xxxyyz_1, g_yyyzzz_xxxyz_1, g_yyyzzz_xxxyzz_1, g_yyyzzz_xxxzz_1, g_yyyzzz_xxxzzz_1, g_yyyzzz_xxyyy_1, g_yyyzzz_xxyyyy_1, g_yyyzzz_xxyyyz_1, g_yyyzzz_xxyyz_1, g_yyyzzz_xxyyzz_1, g_yyyzzz_xxyzz_1, g_yyyzzz_xxyzzz_1, g_yyyzzz_xxzzz_1, g_yyyzzz_xxzzzz_1, g_yyyzzz_xyyyy_1, g_yyyzzz_xyyyyy_1, g_yyyzzz_xyyyyz_1, g_yyyzzz_xyyyz_1, g_yyyzzz_xyyyzz_1, g_yyyzzz_xyyzz_1, g_yyyzzz_xyyzzz_1, g_yyyzzz_xyzzz_1, g_yyyzzz_xyzzzz_1, g_yyyzzz_xzzzz_1, g_yyyzzz_xzzzzz_1, g_yyyzzz_yyyyy_1, g_yyyzzz_yyyyyy_1, g_yyyzzz_yyyyyz_1, g_yyyzzz_yyyyz_1, g_yyyzzz_yyyyzz_1, g_yyyzzz_yyyzz_1, g_yyyzzz_yyyzzz_1, g_yyyzzz_yyzzz_1, g_yyyzzz_yyzzzz_1, g_yyyzzz_yzzzz_1, g_yyyzzz_yzzzzz_1, g_yyyzzz_zzzzz_1, g_yyyzzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzzz_xxxxxx_0[i] = 6.0 * g_yyyzzz_xxxxx_1[i] * fe_0 + g_yyyzzz_xxxxxx_1[i] * pa_x[i];

        g_xyyyzzz_xxxxxy_0[i] = 5.0 * g_yyyzzz_xxxxy_1[i] * fe_0 + g_yyyzzz_xxxxxy_1[i] * pa_x[i];

        g_xyyyzzz_xxxxxz_0[i] = 5.0 * g_yyyzzz_xxxxz_1[i] * fe_0 + g_yyyzzz_xxxxxz_1[i] * pa_x[i];

        g_xyyyzzz_xxxxyy_0[i] = 4.0 * g_yyyzzz_xxxyy_1[i] * fe_0 + g_yyyzzz_xxxxyy_1[i] * pa_x[i];

        g_xyyyzzz_xxxxyz_0[i] = 4.0 * g_yyyzzz_xxxyz_1[i] * fe_0 + g_yyyzzz_xxxxyz_1[i] * pa_x[i];

        g_xyyyzzz_xxxxzz_0[i] = 4.0 * g_yyyzzz_xxxzz_1[i] * fe_0 + g_yyyzzz_xxxxzz_1[i] * pa_x[i];

        g_xyyyzzz_xxxyyy_0[i] = 3.0 * g_yyyzzz_xxyyy_1[i] * fe_0 + g_yyyzzz_xxxyyy_1[i] * pa_x[i];

        g_xyyyzzz_xxxyyz_0[i] = 3.0 * g_yyyzzz_xxyyz_1[i] * fe_0 + g_yyyzzz_xxxyyz_1[i] * pa_x[i];

        g_xyyyzzz_xxxyzz_0[i] = 3.0 * g_yyyzzz_xxyzz_1[i] * fe_0 + g_yyyzzz_xxxyzz_1[i] * pa_x[i];

        g_xyyyzzz_xxxzzz_0[i] = 3.0 * g_yyyzzz_xxzzz_1[i] * fe_0 + g_yyyzzz_xxxzzz_1[i] * pa_x[i];

        g_xyyyzzz_xxyyyy_0[i] = 2.0 * g_yyyzzz_xyyyy_1[i] * fe_0 + g_yyyzzz_xxyyyy_1[i] * pa_x[i];

        g_xyyyzzz_xxyyyz_0[i] = 2.0 * g_yyyzzz_xyyyz_1[i] * fe_0 + g_yyyzzz_xxyyyz_1[i] * pa_x[i];

        g_xyyyzzz_xxyyzz_0[i] = 2.0 * g_yyyzzz_xyyzz_1[i] * fe_0 + g_yyyzzz_xxyyzz_1[i] * pa_x[i];

        g_xyyyzzz_xxyzzz_0[i] = 2.0 * g_yyyzzz_xyzzz_1[i] * fe_0 + g_yyyzzz_xxyzzz_1[i] * pa_x[i];

        g_xyyyzzz_xxzzzz_0[i] = 2.0 * g_yyyzzz_xzzzz_1[i] * fe_0 + g_yyyzzz_xxzzzz_1[i] * pa_x[i];

        g_xyyyzzz_xyyyyy_0[i] = g_yyyzzz_yyyyy_1[i] * fe_0 + g_yyyzzz_xyyyyy_1[i] * pa_x[i];

        g_xyyyzzz_xyyyyz_0[i] = g_yyyzzz_yyyyz_1[i] * fe_0 + g_yyyzzz_xyyyyz_1[i] * pa_x[i];

        g_xyyyzzz_xyyyzz_0[i] = g_yyyzzz_yyyzz_1[i] * fe_0 + g_yyyzzz_xyyyzz_1[i] * pa_x[i];

        g_xyyyzzz_xyyzzz_0[i] = g_yyyzzz_yyzzz_1[i] * fe_0 + g_yyyzzz_xyyzzz_1[i] * pa_x[i];

        g_xyyyzzz_xyzzzz_0[i] = g_yyyzzz_yzzzz_1[i] * fe_0 + g_yyyzzz_xyzzzz_1[i] * pa_x[i];

        g_xyyyzzz_xzzzzz_0[i] = g_yyyzzz_zzzzz_1[i] * fe_0 + g_yyyzzz_xzzzzz_1[i] * pa_x[i];

        g_xyyyzzz_yyyyyy_0[i] = g_yyyzzz_yyyyyy_1[i] * pa_x[i];

        g_xyyyzzz_yyyyyz_0[i] = g_yyyzzz_yyyyyz_1[i] * pa_x[i];

        g_xyyyzzz_yyyyzz_0[i] = g_yyyzzz_yyyyzz_1[i] * pa_x[i];

        g_xyyyzzz_yyyzzz_0[i] = g_yyyzzz_yyyzzz_1[i] * pa_x[i];

        g_xyyyzzz_yyzzzz_0[i] = g_yyyzzz_yyzzzz_1[i] * pa_x[i];

        g_xyyyzzz_yzzzzz_0[i] = g_yyyzzz_yzzzzz_1[i] * pa_x[i];

        g_xyyyzzz_zzzzzz_0[i] = g_yyyzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 700-728 components of targeted buffer : KI

    auto g_xyyzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 700);

    auto g_xyyzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 701);

    auto g_xyyzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 702);

    auto g_xyyzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 703);

    auto g_xyyzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 704);

    auto g_xyyzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 705);

    auto g_xyyzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 706);

    auto g_xyyzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 707);

    auto g_xyyzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 708);

    auto g_xyyzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 709);

    auto g_xyyzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 710);

    auto g_xyyzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 711);

    auto g_xyyzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 712);

    auto g_xyyzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 713);

    auto g_xyyzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 714);

    auto g_xyyzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 715);

    auto g_xyyzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 716);

    auto g_xyyzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 717);

    auto g_xyyzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 718);

    auto g_xyyzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 719);

    auto g_xyyzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 720);

    auto g_xyyzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 721);

    auto g_xyyzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 722);

    auto g_xyyzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 723);

    auto g_xyyzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 724);

    auto g_xyyzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 725);

    auto g_xyyzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 726);

    auto g_xyyzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 727);

    #pragma omp simd aligned(g_xyyzzzz_xxxxxx_0, g_xyyzzzz_xxxxxy_0, g_xyyzzzz_xxxxxz_0, g_xyyzzzz_xxxxyy_0, g_xyyzzzz_xxxxyz_0, g_xyyzzzz_xxxxzz_0, g_xyyzzzz_xxxyyy_0, g_xyyzzzz_xxxyyz_0, g_xyyzzzz_xxxyzz_0, g_xyyzzzz_xxxzzz_0, g_xyyzzzz_xxyyyy_0, g_xyyzzzz_xxyyyz_0, g_xyyzzzz_xxyyzz_0, g_xyyzzzz_xxyzzz_0, g_xyyzzzz_xxzzzz_0, g_xyyzzzz_xyyyyy_0, g_xyyzzzz_xyyyyz_0, g_xyyzzzz_xyyyzz_0, g_xyyzzzz_xyyzzz_0, g_xyyzzzz_xyzzzz_0, g_xyyzzzz_xzzzzz_0, g_xyyzzzz_yyyyyy_0, g_xyyzzzz_yyyyyz_0, g_xyyzzzz_yyyyzz_0, g_xyyzzzz_yyyzzz_0, g_xyyzzzz_yyzzzz_0, g_xyyzzzz_yzzzzz_0, g_xyyzzzz_zzzzzz_0, g_yyzzzz_xxxxx_1, g_yyzzzz_xxxxxx_1, g_yyzzzz_xxxxxy_1, g_yyzzzz_xxxxxz_1, g_yyzzzz_xxxxy_1, g_yyzzzz_xxxxyy_1, g_yyzzzz_xxxxyz_1, g_yyzzzz_xxxxz_1, g_yyzzzz_xxxxzz_1, g_yyzzzz_xxxyy_1, g_yyzzzz_xxxyyy_1, g_yyzzzz_xxxyyz_1, g_yyzzzz_xxxyz_1, g_yyzzzz_xxxyzz_1, g_yyzzzz_xxxzz_1, g_yyzzzz_xxxzzz_1, g_yyzzzz_xxyyy_1, g_yyzzzz_xxyyyy_1, g_yyzzzz_xxyyyz_1, g_yyzzzz_xxyyz_1, g_yyzzzz_xxyyzz_1, g_yyzzzz_xxyzz_1, g_yyzzzz_xxyzzz_1, g_yyzzzz_xxzzz_1, g_yyzzzz_xxzzzz_1, g_yyzzzz_xyyyy_1, g_yyzzzz_xyyyyy_1, g_yyzzzz_xyyyyz_1, g_yyzzzz_xyyyz_1, g_yyzzzz_xyyyzz_1, g_yyzzzz_xyyzz_1, g_yyzzzz_xyyzzz_1, g_yyzzzz_xyzzz_1, g_yyzzzz_xyzzzz_1, g_yyzzzz_xzzzz_1, g_yyzzzz_xzzzzz_1, g_yyzzzz_yyyyy_1, g_yyzzzz_yyyyyy_1, g_yyzzzz_yyyyyz_1, g_yyzzzz_yyyyz_1, g_yyzzzz_yyyyzz_1, g_yyzzzz_yyyzz_1, g_yyzzzz_yyyzzz_1, g_yyzzzz_yyzzz_1, g_yyzzzz_yyzzzz_1, g_yyzzzz_yzzzz_1, g_yyzzzz_yzzzzz_1, g_yyzzzz_zzzzz_1, g_yyzzzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzzz_xxxxxx_0[i] = 6.0 * g_yyzzzz_xxxxx_1[i] * fe_0 + g_yyzzzz_xxxxxx_1[i] * pa_x[i];

        g_xyyzzzz_xxxxxy_0[i] = 5.0 * g_yyzzzz_xxxxy_1[i] * fe_0 + g_yyzzzz_xxxxxy_1[i] * pa_x[i];

        g_xyyzzzz_xxxxxz_0[i] = 5.0 * g_yyzzzz_xxxxz_1[i] * fe_0 + g_yyzzzz_xxxxxz_1[i] * pa_x[i];

        g_xyyzzzz_xxxxyy_0[i] = 4.0 * g_yyzzzz_xxxyy_1[i] * fe_0 + g_yyzzzz_xxxxyy_1[i] * pa_x[i];

        g_xyyzzzz_xxxxyz_0[i] = 4.0 * g_yyzzzz_xxxyz_1[i] * fe_0 + g_yyzzzz_xxxxyz_1[i] * pa_x[i];

        g_xyyzzzz_xxxxzz_0[i] = 4.0 * g_yyzzzz_xxxzz_1[i] * fe_0 + g_yyzzzz_xxxxzz_1[i] * pa_x[i];

        g_xyyzzzz_xxxyyy_0[i] = 3.0 * g_yyzzzz_xxyyy_1[i] * fe_0 + g_yyzzzz_xxxyyy_1[i] * pa_x[i];

        g_xyyzzzz_xxxyyz_0[i] = 3.0 * g_yyzzzz_xxyyz_1[i] * fe_0 + g_yyzzzz_xxxyyz_1[i] * pa_x[i];

        g_xyyzzzz_xxxyzz_0[i] = 3.0 * g_yyzzzz_xxyzz_1[i] * fe_0 + g_yyzzzz_xxxyzz_1[i] * pa_x[i];

        g_xyyzzzz_xxxzzz_0[i] = 3.0 * g_yyzzzz_xxzzz_1[i] * fe_0 + g_yyzzzz_xxxzzz_1[i] * pa_x[i];

        g_xyyzzzz_xxyyyy_0[i] = 2.0 * g_yyzzzz_xyyyy_1[i] * fe_0 + g_yyzzzz_xxyyyy_1[i] * pa_x[i];

        g_xyyzzzz_xxyyyz_0[i] = 2.0 * g_yyzzzz_xyyyz_1[i] * fe_0 + g_yyzzzz_xxyyyz_1[i] * pa_x[i];

        g_xyyzzzz_xxyyzz_0[i] = 2.0 * g_yyzzzz_xyyzz_1[i] * fe_0 + g_yyzzzz_xxyyzz_1[i] * pa_x[i];

        g_xyyzzzz_xxyzzz_0[i] = 2.0 * g_yyzzzz_xyzzz_1[i] * fe_0 + g_yyzzzz_xxyzzz_1[i] * pa_x[i];

        g_xyyzzzz_xxzzzz_0[i] = 2.0 * g_yyzzzz_xzzzz_1[i] * fe_0 + g_yyzzzz_xxzzzz_1[i] * pa_x[i];

        g_xyyzzzz_xyyyyy_0[i] = g_yyzzzz_yyyyy_1[i] * fe_0 + g_yyzzzz_xyyyyy_1[i] * pa_x[i];

        g_xyyzzzz_xyyyyz_0[i] = g_yyzzzz_yyyyz_1[i] * fe_0 + g_yyzzzz_xyyyyz_1[i] * pa_x[i];

        g_xyyzzzz_xyyyzz_0[i] = g_yyzzzz_yyyzz_1[i] * fe_0 + g_yyzzzz_xyyyzz_1[i] * pa_x[i];

        g_xyyzzzz_xyyzzz_0[i] = g_yyzzzz_yyzzz_1[i] * fe_0 + g_yyzzzz_xyyzzz_1[i] * pa_x[i];

        g_xyyzzzz_xyzzzz_0[i] = g_yyzzzz_yzzzz_1[i] * fe_0 + g_yyzzzz_xyzzzz_1[i] * pa_x[i];

        g_xyyzzzz_xzzzzz_0[i] = g_yyzzzz_zzzzz_1[i] * fe_0 + g_yyzzzz_xzzzzz_1[i] * pa_x[i];

        g_xyyzzzz_yyyyyy_0[i] = g_yyzzzz_yyyyyy_1[i] * pa_x[i];

        g_xyyzzzz_yyyyyz_0[i] = g_yyzzzz_yyyyyz_1[i] * pa_x[i];

        g_xyyzzzz_yyyyzz_0[i] = g_yyzzzz_yyyyzz_1[i] * pa_x[i];

        g_xyyzzzz_yyyzzz_0[i] = g_yyzzzz_yyyzzz_1[i] * pa_x[i];

        g_xyyzzzz_yyzzzz_0[i] = g_yyzzzz_yyzzzz_1[i] * pa_x[i];

        g_xyyzzzz_yzzzzz_0[i] = g_yyzzzz_yzzzzz_1[i] * pa_x[i];

        g_xyyzzzz_zzzzzz_0[i] = g_yyzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 728-756 components of targeted buffer : KI

    auto g_xyzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 728);

    auto g_xyzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 729);

    auto g_xyzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 730);

    auto g_xyzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 731);

    auto g_xyzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 732);

    auto g_xyzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 733);

    auto g_xyzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 734);

    auto g_xyzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 735);

    auto g_xyzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 736);

    auto g_xyzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 737);

    auto g_xyzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 738);

    auto g_xyzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 739);

    auto g_xyzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 740);

    auto g_xyzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 741);

    auto g_xyzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 742);

    auto g_xyzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 743);

    auto g_xyzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 744);

    auto g_xyzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 745);

    auto g_xyzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 746);

    auto g_xyzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 747);

    auto g_xyzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 748);

    auto g_xyzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 749);

    auto g_xyzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 750);

    auto g_xyzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 751);

    auto g_xyzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 752);

    auto g_xyzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 753);

    auto g_xyzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 754);

    auto g_xyzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 755);

    #pragma omp simd aligned(g_xyzzzzz_xxxxxx_0, g_xyzzzzz_xxxxxy_0, g_xyzzzzz_xxxxxz_0, g_xyzzzzz_xxxxyy_0, g_xyzzzzz_xxxxyz_0, g_xyzzzzz_xxxxzz_0, g_xyzzzzz_xxxyyy_0, g_xyzzzzz_xxxyyz_0, g_xyzzzzz_xxxyzz_0, g_xyzzzzz_xxxzzz_0, g_xyzzzzz_xxyyyy_0, g_xyzzzzz_xxyyyz_0, g_xyzzzzz_xxyyzz_0, g_xyzzzzz_xxyzzz_0, g_xyzzzzz_xxzzzz_0, g_xyzzzzz_xyyyyy_0, g_xyzzzzz_xyyyyz_0, g_xyzzzzz_xyyyzz_0, g_xyzzzzz_xyyzzz_0, g_xyzzzzz_xyzzzz_0, g_xyzzzzz_xzzzzz_0, g_xyzzzzz_yyyyyy_0, g_xyzzzzz_yyyyyz_0, g_xyzzzzz_yyyyzz_0, g_xyzzzzz_yyyzzz_0, g_xyzzzzz_yyzzzz_0, g_xyzzzzz_yzzzzz_0, g_xyzzzzz_zzzzzz_0, g_xzzzzz_xxxxxx_1, g_xzzzzz_xxxxxz_1, g_xzzzzz_xxxxzz_1, g_xzzzzz_xxxzzz_1, g_xzzzzz_xxzzzz_1, g_xzzzzz_xzzzzz_1, g_yzzzzz_xxxxxy_1, g_yzzzzz_xxxxy_1, g_yzzzzz_xxxxyy_1, g_yzzzzz_xxxxyz_1, g_yzzzzz_xxxyy_1, g_yzzzzz_xxxyyy_1, g_yzzzzz_xxxyyz_1, g_yzzzzz_xxxyz_1, g_yzzzzz_xxxyzz_1, g_yzzzzz_xxyyy_1, g_yzzzzz_xxyyyy_1, g_yzzzzz_xxyyyz_1, g_yzzzzz_xxyyz_1, g_yzzzzz_xxyyzz_1, g_yzzzzz_xxyzz_1, g_yzzzzz_xxyzzz_1, g_yzzzzz_xyyyy_1, g_yzzzzz_xyyyyy_1, g_yzzzzz_xyyyyz_1, g_yzzzzz_xyyyz_1, g_yzzzzz_xyyyzz_1, g_yzzzzz_xyyzz_1, g_yzzzzz_xyyzzz_1, g_yzzzzz_xyzzz_1, g_yzzzzz_xyzzzz_1, g_yzzzzz_yyyyy_1, g_yzzzzz_yyyyyy_1, g_yzzzzz_yyyyyz_1, g_yzzzzz_yyyyz_1, g_yzzzzz_yyyyzz_1, g_yzzzzz_yyyzz_1, g_yzzzzz_yyyzzz_1, g_yzzzzz_yyzzz_1, g_yzzzzz_yyzzzz_1, g_yzzzzz_yzzzz_1, g_yzzzzz_yzzzzz_1, g_yzzzzz_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzzz_xxxxxx_0[i] = g_xzzzzz_xxxxxx_1[i] * pa_y[i];

        g_xyzzzzz_xxxxxy_0[i] = 5.0 * g_yzzzzz_xxxxy_1[i] * fe_0 + g_yzzzzz_xxxxxy_1[i] * pa_x[i];

        g_xyzzzzz_xxxxxz_0[i] = g_xzzzzz_xxxxxz_1[i] * pa_y[i];

        g_xyzzzzz_xxxxyy_0[i] = 4.0 * g_yzzzzz_xxxyy_1[i] * fe_0 + g_yzzzzz_xxxxyy_1[i] * pa_x[i];

        g_xyzzzzz_xxxxyz_0[i] = 4.0 * g_yzzzzz_xxxyz_1[i] * fe_0 + g_yzzzzz_xxxxyz_1[i] * pa_x[i];

        g_xyzzzzz_xxxxzz_0[i] = g_xzzzzz_xxxxzz_1[i] * pa_y[i];

        g_xyzzzzz_xxxyyy_0[i] = 3.0 * g_yzzzzz_xxyyy_1[i] * fe_0 + g_yzzzzz_xxxyyy_1[i] * pa_x[i];

        g_xyzzzzz_xxxyyz_0[i] = 3.0 * g_yzzzzz_xxyyz_1[i] * fe_0 + g_yzzzzz_xxxyyz_1[i] * pa_x[i];

        g_xyzzzzz_xxxyzz_0[i] = 3.0 * g_yzzzzz_xxyzz_1[i] * fe_0 + g_yzzzzz_xxxyzz_1[i] * pa_x[i];

        g_xyzzzzz_xxxzzz_0[i] = g_xzzzzz_xxxzzz_1[i] * pa_y[i];

        g_xyzzzzz_xxyyyy_0[i] = 2.0 * g_yzzzzz_xyyyy_1[i] * fe_0 + g_yzzzzz_xxyyyy_1[i] * pa_x[i];

        g_xyzzzzz_xxyyyz_0[i] = 2.0 * g_yzzzzz_xyyyz_1[i] * fe_0 + g_yzzzzz_xxyyyz_1[i] * pa_x[i];

        g_xyzzzzz_xxyyzz_0[i] = 2.0 * g_yzzzzz_xyyzz_1[i] * fe_0 + g_yzzzzz_xxyyzz_1[i] * pa_x[i];

        g_xyzzzzz_xxyzzz_0[i] = 2.0 * g_yzzzzz_xyzzz_1[i] * fe_0 + g_yzzzzz_xxyzzz_1[i] * pa_x[i];

        g_xyzzzzz_xxzzzz_0[i] = g_xzzzzz_xxzzzz_1[i] * pa_y[i];

        g_xyzzzzz_xyyyyy_0[i] = g_yzzzzz_yyyyy_1[i] * fe_0 + g_yzzzzz_xyyyyy_1[i] * pa_x[i];

        g_xyzzzzz_xyyyyz_0[i] = g_yzzzzz_yyyyz_1[i] * fe_0 + g_yzzzzz_xyyyyz_1[i] * pa_x[i];

        g_xyzzzzz_xyyyzz_0[i] = g_yzzzzz_yyyzz_1[i] * fe_0 + g_yzzzzz_xyyyzz_1[i] * pa_x[i];

        g_xyzzzzz_xyyzzz_0[i] = g_yzzzzz_yyzzz_1[i] * fe_0 + g_yzzzzz_xyyzzz_1[i] * pa_x[i];

        g_xyzzzzz_xyzzzz_0[i] = g_yzzzzz_yzzzz_1[i] * fe_0 + g_yzzzzz_xyzzzz_1[i] * pa_x[i];

        g_xyzzzzz_xzzzzz_0[i] = g_xzzzzz_xzzzzz_1[i] * pa_y[i];

        g_xyzzzzz_yyyyyy_0[i] = g_yzzzzz_yyyyyy_1[i] * pa_x[i];

        g_xyzzzzz_yyyyyz_0[i] = g_yzzzzz_yyyyyz_1[i] * pa_x[i];

        g_xyzzzzz_yyyyzz_0[i] = g_yzzzzz_yyyyzz_1[i] * pa_x[i];

        g_xyzzzzz_yyyzzz_0[i] = g_yzzzzz_yyyzzz_1[i] * pa_x[i];

        g_xyzzzzz_yyzzzz_0[i] = g_yzzzzz_yyzzzz_1[i] * pa_x[i];

        g_xyzzzzz_yzzzzz_0[i] = g_yzzzzz_yzzzzz_1[i] * pa_x[i];

        g_xyzzzzz_zzzzzz_0[i] = g_yzzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 756-784 components of targeted buffer : KI

    auto g_xzzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 756);

    auto g_xzzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 757);

    auto g_xzzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 758);

    auto g_xzzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 759);

    auto g_xzzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 760);

    auto g_xzzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 761);

    auto g_xzzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 762);

    auto g_xzzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 763);

    auto g_xzzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 764);

    auto g_xzzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 765);

    auto g_xzzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 766);

    auto g_xzzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 767);

    auto g_xzzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 768);

    auto g_xzzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 769);

    auto g_xzzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 770);

    auto g_xzzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 771);

    auto g_xzzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 772);

    auto g_xzzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 773);

    auto g_xzzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 774);

    auto g_xzzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 775);

    auto g_xzzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 776);

    auto g_xzzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 777);

    auto g_xzzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 778);

    auto g_xzzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 779);

    auto g_xzzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 780);

    auto g_xzzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 781);

    auto g_xzzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 782);

    auto g_xzzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 783);

    #pragma omp simd aligned(g_xzzzzzz_xxxxxx_0, g_xzzzzzz_xxxxxy_0, g_xzzzzzz_xxxxxz_0, g_xzzzzzz_xxxxyy_0, g_xzzzzzz_xxxxyz_0, g_xzzzzzz_xxxxzz_0, g_xzzzzzz_xxxyyy_0, g_xzzzzzz_xxxyyz_0, g_xzzzzzz_xxxyzz_0, g_xzzzzzz_xxxzzz_0, g_xzzzzzz_xxyyyy_0, g_xzzzzzz_xxyyyz_0, g_xzzzzzz_xxyyzz_0, g_xzzzzzz_xxyzzz_0, g_xzzzzzz_xxzzzz_0, g_xzzzzzz_xyyyyy_0, g_xzzzzzz_xyyyyz_0, g_xzzzzzz_xyyyzz_0, g_xzzzzzz_xyyzzz_0, g_xzzzzzz_xyzzzz_0, g_xzzzzzz_xzzzzz_0, g_xzzzzzz_yyyyyy_0, g_xzzzzzz_yyyyyz_0, g_xzzzzzz_yyyyzz_0, g_xzzzzzz_yyyzzz_0, g_xzzzzzz_yyzzzz_0, g_xzzzzzz_yzzzzz_0, g_xzzzzzz_zzzzzz_0, g_zzzzzz_xxxxx_1, g_zzzzzz_xxxxxx_1, g_zzzzzz_xxxxxy_1, g_zzzzzz_xxxxxz_1, g_zzzzzz_xxxxy_1, g_zzzzzz_xxxxyy_1, g_zzzzzz_xxxxyz_1, g_zzzzzz_xxxxz_1, g_zzzzzz_xxxxzz_1, g_zzzzzz_xxxyy_1, g_zzzzzz_xxxyyy_1, g_zzzzzz_xxxyyz_1, g_zzzzzz_xxxyz_1, g_zzzzzz_xxxyzz_1, g_zzzzzz_xxxzz_1, g_zzzzzz_xxxzzz_1, g_zzzzzz_xxyyy_1, g_zzzzzz_xxyyyy_1, g_zzzzzz_xxyyyz_1, g_zzzzzz_xxyyz_1, g_zzzzzz_xxyyzz_1, g_zzzzzz_xxyzz_1, g_zzzzzz_xxyzzz_1, g_zzzzzz_xxzzz_1, g_zzzzzz_xxzzzz_1, g_zzzzzz_xyyyy_1, g_zzzzzz_xyyyyy_1, g_zzzzzz_xyyyyz_1, g_zzzzzz_xyyyz_1, g_zzzzzz_xyyyzz_1, g_zzzzzz_xyyzz_1, g_zzzzzz_xyyzzz_1, g_zzzzzz_xyzzz_1, g_zzzzzz_xyzzzz_1, g_zzzzzz_xzzzz_1, g_zzzzzz_xzzzzz_1, g_zzzzzz_yyyyy_1, g_zzzzzz_yyyyyy_1, g_zzzzzz_yyyyyz_1, g_zzzzzz_yyyyz_1, g_zzzzzz_yyyyzz_1, g_zzzzzz_yyyzz_1, g_zzzzzz_yyyzzz_1, g_zzzzzz_yyzzz_1, g_zzzzzz_yyzzzz_1, g_zzzzzz_yzzzz_1, g_zzzzzz_yzzzzz_1, g_zzzzzz_zzzzz_1, g_zzzzzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzzz_xxxxxx_0[i] = 6.0 * g_zzzzzz_xxxxx_1[i] * fe_0 + g_zzzzzz_xxxxxx_1[i] * pa_x[i];

        g_xzzzzzz_xxxxxy_0[i] = 5.0 * g_zzzzzz_xxxxy_1[i] * fe_0 + g_zzzzzz_xxxxxy_1[i] * pa_x[i];

        g_xzzzzzz_xxxxxz_0[i] = 5.0 * g_zzzzzz_xxxxz_1[i] * fe_0 + g_zzzzzz_xxxxxz_1[i] * pa_x[i];

        g_xzzzzzz_xxxxyy_0[i] = 4.0 * g_zzzzzz_xxxyy_1[i] * fe_0 + g_zzzzzz_xxxxyy_1[i] * pa_x[i];

        g_xzzzzzz_xxxxyz_0[i] = 4.0 * g_zzzzzz_xxxyz_1[i] * fe_0 + g_zzzzzz_xxxxyz_1[i] * pa_x[i];

        g_xzzzzzz_xxxxzz_0[i] = 4.0 * g_zzzzzz_xxxzz_1[i] * fe_0 + g_zzzzzz_xxxxzz_1[i] * pa_x[i];

        g_xzzzzzz_xxxyyy_0[i] = 3.0 * g_zzzzzz_xxyyy_1[i] * fe_0 + g_zzzzzz_xxxyyy_1[i] * pa_x[i];

        g_xzzzzzz_xxxyyz_0[i] = 3.0 * g_zzzzzz_xxyyz_1[i] * fe_0 + g_zzzzzz_xxxyyz_1[i] * pa_x[i];

        g_xzzzzzz_xxxyzz_0[i] = 3.0 * g_zzzzzz_xxyzz_1[i] * fe_0 + g_zzzzzz_xxxyzz_1[i] * pa_x[i];

        g_xzzzzzz_xxxzzz_0[i] = 3.0 * g_zzzzzz_xxzzz_1[i] * fe_0 + g_zzzzzz_xxxzzz_1[i] * pa_x[i];

        g_xzzzzzz_xxyyyy_0[i] = 2.0 * g_zzzzzz_xyyyy_1[i] * fe_0 + g_zzzzzz_xxyyyy_1[i] * pa_x[i];

        g_xzzzzzz_xxyyyz_0[i] = 2.0 * g_zzzzzz_xyyyz_1[i] * fe_0 + g_zzzzzz_xxyyyz_1[i] * pa_x[i];

        g_xzzzzzz_xxyyzz_0[i] = 2.0 * g_zzzzzz_xyyzz_1[i] * fe_0 + g_zzzzzz_xxyyzz_1[i] * pa_x[i];

        g_xzzzzzz_xxyzzz_0[i] = 2.0 * g_zzzzzz_xyzzz_1[i] * fe_0 + g_zzzzzz_xxyzzz_1[i] * pa_x[i];

        g_xzzzzzz_xxzzzz_0[i] = 2.0 * g_zzzzzz_xzzzz_1[i] * fe_0 + g_zzzzzz_xxzzzz_1[i] * pa_x[i];

        g_xzzzzzz_xyyyyy_0[i] = g_zzzzzz_yyyyy_1[i] * fe_0 + g_zzzzzz_xyyyyy_1[i] * pa_x[i];

        g_xzzzzzz_xyyyyz_0[i] = g_zzzzzz_yyyyz_1[i] * fe_0 + g_zzzzzz_xyyyyz_1[i] * pa_x[i];

        g_xzzzzzz_xyyyzz_0[i] = g_zzzzzz_yyyzz_1[i] * fe_0 + g_zzzzzz_xyyyzz_1[i] * pa_x[i];

        g_xzzzzzz_xyyzzz_0[i] = g_zzzzzz_yyzzz_1[i] * fe_0 + g_zzzzzz_xyyzzz_1[i] * pa_x[i];

        g_xzzzzzz_xyzzzz_0[i] = g_zzzzzz_yzzzz_1[i] * fe_0 + g_zzzzzz_xyzzzz_1[i] * pa_x[i];

        g_xzzzzzz_xzzzzz_0[i] = g_zzzzzz_zzzzz_1[i] * fe_0 + g_zzzzzz_xzzzzz_1[i] * pa_x[i];

        g_xzzzzzz_yyyyyy_0[i] = g_zzzzzz_yyyyyy_1[i] * pa_x[i];

        g_xzzzzzz_yyyyyz_0[i] = g_zzzzzz_yyyyyz_1[i] * pa_x[i];

        g_xzzzzzz_yyyyzz_0[i] = g_zzzzzz_yyyyzz_1[i] * pa_x[i];

        g_xzzzzzz_yyyzzz_0[i] = g_zzzzzz_yyyzzz_1[i] * pa_x[i];

        g_xzzzzzz_yyzzzz_0[i] = g_zzzzzz_yyzzzz_1[i] * pa_x[i];

        g_xzzzzzz_yzzzzz_0[i] = g_zzzzzz_yzzzzz_1[i] * pa_x[i];

        g_xzzzzzz_zzzzzz_0[i] = g_zzzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 784-812 components of targeted buffer : KI

    auto g_yyyyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 784);

    auto g_yyyyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 785);

    auto g_yyyyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 786);

    auto g_yyyyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 787);

    auto g_yyyyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 788);

    auto g_yyyyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 789);

    auto g_yyyyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 790);

    auto g_yyyyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 791);

    auto g_yyyyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 792);

    auto g_yyyyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 793);

    auto g_yyyyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 794);

    auto g_yyyyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 795);

    auto g_yyyyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 796);

    auto g_yyyyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 797);

    auto g_yyyyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 798);

    auto g_yyyyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 799);

    auto g_yyyyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 800);

    auto g_yyyyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 801);

    auto g_yyyyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 802);

    auto g_yyyyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 803);

    auto g_yyyyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 804);

    auto g_yyyyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 805);

    auto g_yyyyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 806);

    auto g_yyyyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 807);

    auto g_yyyyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 808);

    auto g_yyyyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 809);

    auto g_yyyyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 810);

    auto g_yyyyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 811);

    #pragma omp simd aligned(g_yyyyy_xxxxxx_0, g_yyyyy_xxxxxx_1, g_yyyyy_xxxxxy_0, g_yyyyy_xxxxxy_1, g_yyyyy_xxxxxz_0, g_yyyyy_xxxxxz_1, g_yyyyy_xxxxyy_0, g_yyyyy_xxxxyy_1, g_yyyyy_xxxxyz_0, g_yyyyy_xxxxyz_1, g_yyyyy_xxxxzz_0, g_yyyyy_xxxxzz_1, g_yyyyy_xxxyyy_0, g_yyyyy_xxxyyy_1, g_yyyyy_xxxyyz_0, g_yyyyy_xxxyyz_1, g_yyyyy_xxxyzz_0, g_yyyyy_xxxyzz_1, g_yyyyy_xxxzzz_0, g_yyyyy_xxxzzz_1, g_yyyyy_xxyyyy_0, g_yyyyy_xxyyyy_1, g_yyyyy_xxyyyz_0, g_yyyyy_xxyyyz_1, g_yyyyy_xxyyzz_0, g_yyyyy_xxyyzz_1, g_yyyyy_xxyzzz_0, g_yyyyy_xxyzzz_1, g_yyyyy_xxzzzz_0, g_yyyyy_xxzzzz_1, g_yyyyy_xyyyyy_0, g_yyyyy_xyyyyy_1, g_yyyyy_xyyyyz_0, g_yyyyy_xyyyyz_1, g_yyyyy_xyyyzz_0, g_yyyyy_xyyyzz_1, g_yyyyy_xyyzzz_0, g_yyyyy_xyyzzz_1, g_yyyyy_xyzzzz_0, g_yyyyy_xyzzzz_1, g_yyyyy_xzzzzz_0, g_yyyyy_xzzzzz_1, g_yyyyy_yyyyyy_0, g_yyyyy_yyyyyy_1, g_yyyyy_yyyyyz_0, g_yyyyy_yyyyyz_1, g_yyyyy_yyyyzz_0, g_yyyyy_yyyyzz_1, g_yyyyy_yyyzzz_0, g_yyyyy_yyyzzz_1, g_yyyyy_yyzzzz_0, g_yyyyy_yyzzzz_1, g_yyyyy_yzzzzz_0, g_yyyyy_yzzzzz_1, g_yyyyy_zzzzzz_0, g_yyyyy_zzzzzz_1, g_yyyyyy_xxxxx_1, g_yyyyyy_xxxxxx_1, g_yyyyyy_xxxxxy_1, g_yyyyyy_xxxxxz_1, g_yyyyyy_xxxxy_1, g_yyyyyy_xxxxyy_1, g_yyyyyy_xxxxyz_1, g_yyyyyy_xxxxz_1, g_yyyyyy_xxxxzz_1, g_yyyyyy_xxxyy_1, g_yyyyyy_xxxyyy_1, g_yyyyyy_xxxyyz_1, g_yyyyyy_xxxyz_1, g_yyyyyy_xxxyzz_1, g_yyyyyy_xxxzz_1, g_yyyyyy_xxxzzz_1, g_yyyyyy_xxyyy_1, g_yyyyyy_xxyyyy_1, g_yyyyyy_xxyyyz_1, g_yyyyyy_xxyyz_1, g_yyyyyy_xxyyzz_1, g_yyyyyy_xxyzz_1, g_yyyyyy_xxyzzz_1, g_yyyyyy_xxzzz_1, g_yyyyyy_xxzzzz_1, g_yyyyyy_xyyyy_1, g_yyyyyy_xyyyyy_1, g_yyyyyy_xyyyyz_1, g_yyyyyy_xyyyz_1, g_yyyyyy_xyyyzz_1, g_yyyyyy_xyyzz_1, g_yyyyyy_xyyzzz_1, g_yyyyyy_xyzzz_1, g_yyyyyy_xyzzzz_1, g_yyyyyy_xzzzz_1, g_yyyyyy_xzzzzz_1, g_yyyyyy_yyyyy_1, g_yyyyyy_yyyyyy_1, g_yyyyyy_yyyyyz_1, g_yyyyyy_yyyyz_1, g_yyyyyy_yyyyzz_1, g_yyyyyy_yyyzz_1, g_yyyyyy_yyyzzz_1, g_yyyyyy_yyzzz_1, g_yyyyyy_yyzzzz_1, g_yyyyyy_yzzzz_1, g_yyyyyy_yzzzzz_1, g_yyyyyy_zzzzz_1, g_yyyyyy_zzzzzz_1, g_yyyyyyy_xxxxxx_0, g_yyyyyyy_xxxxxy_0, g_yyyyyyy_xxxxxz_0, g_yyyyyyy_xxxxyy_0, g_yyyyyyy_xxxxyz_0, g_yyyyyyy_xxxxzz_0, g_yyyyyyy_xxxyyy_0, g_yyyyyyy_xxxyyz_0, g_yyyyyyy_xxxyzz_0, g_yyyyyyy_xxxzzz_0, g_yyyyyyy_xxyyyy_0, g_yyyyyyy_xxyyyz_0, g_yyyyyyy_xxyyzz_0, g_yyyyyyy_xxyzzz_0, g_yyyyyyy_xxzzzz_0, g_yyyyyyy_xyyyyy_0, g_yyyyyyy_xyyyyz_0, g_yyyyyyy_xyyyzz_0, g_yyyyyyy_xyyzzz_0, g_yyyyyyy_xyzzzz_0, g_yyyyyyy_xzzzzz_0, g_yyyyyyy_yyyyyy_0, g_yyyyyyy_yyyyyz_0, g_yyyyyyy_yyyyzz_0, g_yyyyyyy_yyyzzz_0, g_yyyyyyy_yyzzzz_0, g_yyyyyyy_yzzzzz_0, g_yyyyyyy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyyy_xxxxxx_0[i] = 6.0 * g_yyyyy_xxxxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxxx_1[i] * fz_be_0 + g_yyyyyy_xxxxxx_1[i] * pa_y[i];

        g_yyyyyyy_xxxxxy_0[i] = 6.0 * g_yyyyy_xxxxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxxy_1[i] * fz_be_0 + g_yyyyyy_xxxxx_1[i] * fe_0 + g_yyyyyy_xxxxxy_1[i] * pa_y[i];

        g_yyyyyyy_xxxxxz_0[i] = 6.0 * g_yyyyy_xxxxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxxz_1[i] * fz_be_0 + g_yyyyyy_xxxxxz_1[i] * pa_y[i];

        g_yyyyyyy_xxxxyy_0[i] = 6.0 * g_yyyyy_xxxxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xxxxy_1[i] * fe_0 + g_yyyyyy_xxxxyy_1[i] * pa_y[i];

        g_yyyyyyy_xxxxyz_0[i] = 6.0 * g_yyyyy_xxxxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxyz_1[i] * fz_be_0 + g_yyyyyy_xxxxz_1[i] * fe_0 + g_yyyyyy_xxxxyz_1[i] * pa_y[i];

        g_yyyyyyy_xxxxzz_0[i] = 6.0 * g_yyyyy_xxxxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxzz_1[i] * fz_be_0 + g_yyyyyy_xxxxzz_1[i] * pa_y[i];

        g_yyyyyyy_xxxyyy_0[i] = 6.0 * g_yyyyy_xxxyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_xxxyy_1[i] * fe_0 + g_yyyyyy_xxxyyy_1[i] * pa_y[i];

        g_yyyyyyy_xxxyyz_0[i] = 6.0 * g_yyyyy_xxxyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xxxyz_1[i] * fe_0 + g_yyyyyy_xxxyyz_1[i] * pa_y[i];

        g_yyyyyyy_xxxyzz_0[i] = 6.0 * g_yyyyy_xxxyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxyzz_1[i] * fz_be_0 + g_yyyyyy_xxxzz_1[i] * fe_0 + g_yyyyyy_xxxyzz_1[i] * pa_y[i];

        g_yyyyyyy_xxxzzz_0[i] = 6.0 * g_yyyyy_xxxzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxzzz_1[i] * fz_be_0 + g_yyyyyy_xxxzzz_1[i] * pa_y[i];

        g_yyyyyyy_xxyyyy_0[i] = 6.0 * g_yyyyy_xxyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_xxyyy_1[i] * fe_0 + g_yyyyyy_xxyyyy_1[i] * pa_y[i];

        g_yyyyyyy_xxyyyz_0[i] = 6.0 * g_yyyyy_xxyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_xxyyz_1[i] * fe_0 + g_yyyyyy_xxyyyz_1[i] * pa_y[i];

        g_yyyyyyy_xxyyzz_0[i] = 6.0 * g_yyyyy_xxyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xxyzz_1[i] * fe_0 + g_yyyyyy_xxyyzz_1[i] * pa_y[i];

        g_yyyyyyy_xxyzzz_0[i] = 6.0 * g_yyyyy_xxyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyzzz_1[i] * fz_be_0 + g_yyyyyy_xxzzz_1[i] * fe_0 + g_yyyyyy_xxyzzz_1[i] * pa_y[i];

        g_yyyyyyy_xxzzzz_0[i] = 6.0 * g_yyyyy_xxzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxzzzz_1[i] * fz_be_0 + g_yyyyyy_xxzzzz_1[i] * pa_y[i];

        g_yyyyyyy_xyyyyy_0[i] = 6.0 * g_yyyyy_xyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyyy_xyyyy_1[i] * fe_0 + g_yyyyyy_xyyyyy_1[i] * pa_y[i];

        g_yyyyyyy_xyyyyz_0[i] = 6.0 * g_yyyyy_xyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_xyyyz_1[i] * fe_0 + g_yyyyyy_xyyyyz_1[i] * pa_y[i];

        g_yyyyyyy_xyyyzz_0[i] = 6.0 * g_yyyyy_xyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_xyyzz_1[i] * fe_0 + g_yyyyyy_xyyyzz_1[i] * pa_y[i];

        g_yyyyyyy_xyyzzz_0[i] = 6.0 * g_yyyyy_xyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xyzzz_1[i] * fe_0 + g_yyyyyy_xyyzzz_1[i] * pa_y[i];

        g_yyyyyyy_xyzzzz_0[i] = 6.0 * g_yyyyy_xyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyzzzz_1[i] * fz_be_0 + g_yyyyyy_xzzzz_1[i] * fe_0 + g_yyyyyy_xyzzzz_1[i] * pa_y[i];

        g_yyyyyyy_xzzzzz_0[i] = 6.0 * g_yyyyy_xzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xzzzzz_1[i] * fz_be_0 + g_yyyyyy_xzzzzz_1[i] * pa_y[i];

        g_yyyyyyy_yyyyyy_0[i] = 6.0 * g_yyyyy_yyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyyy_yyyyy_1[i] * fe_0 + g_yyyyyy_yyyyyy_1[i] * pa_y[i];

        g_yyyyyyy_yyyyyz_0[i] = 6.0 * g_yyyyy_yyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_yyyyz_1[i] * fe_0 + g_yyyyyy_yyyyyz_1[i] * pa_y[i];

        g_yyyyyyy_yyyyzz_0[i] = 6.0 * g_yyyyy_yyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_yyyzz_1[i] * fe_0 + g_yyyyyy_yyyyzz_1[i] * pa_y[i];

        g_yyyyyyy_yyyzzz_0[i] = 6.0 * g_yyyyy_yyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_yyzzz_1[i] * fe_0 + g_yyyyyy_yyyzzz_1[i] * pa_y[i];

        g_yyyyyyy_yyzzzz_0[i] = 6.0 * g_yyyyy_yyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_yzzzz_1[i] * fe_0 + g_yyyyyy_yyzzzz_1[i] * pa_y[i];

        g_yyyyyyy_yzzzzz_0[i] = 6.0 * g_yyyyy_yzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yzzzzz_1[i] * fz_be_0 + g_yyyyyy_zzzzz_1[i] * fe_0 + g_yyyyyy_yzzzzz_1[i] * pa_y[i];

        g_yyyyyyy_zzzzzz_0[i] = 6.0 * g_yyyyy_zzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_zzzzzz_1[i] * fz_be_0 + g_yyyyyy_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 812-840 components of targeted buffer : KI

    auto g_yyyyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 812);

    auto g_yyyyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 813);

    auto g_yyyyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 814);

    auto g_yyyyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 815);

    auto g_yyyyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 816);

    auto g_yyyyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 817);

    auto g_yyyyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 818);

    auto g_yyyyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 819);

    auto g_yyyyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 820);

    auto g_yyyyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 821);

    auto g_yyyyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 822);

    auto g_yyyyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 823);

    auto g_yyyyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 824);

    auto g_yyyyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 825);

    auto g_yyyyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 826);

    auto g_yyyyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 827);

    auto g_yyyyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 828);

    auto g_yyyyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 829);

    auto g_yyyyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 830);

    auto g_yyyyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 831);

    auto g_yyyyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 832);

    auto g_yyyyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 833);

    auto g_yyyyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 834);

    auto g_yyyyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 835);

    auto g_yyyyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 836);

    auto g_yyyyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 837);

    auto g_yyyyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 838);

    auto g_yyyyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 839);

    #pragma omp simd aligned(g_yyyyyy_xxxxx_1, g_yyyyyy_xxxxxx_1, g_yyyyyy_xxxxxy_1, g_yyyyyy_xxxxxz_1, g_yyyyyy_xxxxy_1, g_yyyyyy_xxxxyy_1, g_yyyyyy_xxxxyz_1, g_yyyyyy_xxxxz_1, g_yyyyyy_xxxxzz_1, g_yyyyyy_xxxyy_1, g_yyyyyy_xxxyyy_1, g_yyyyyy_xxxyyz_1, g_yyyyyy_xxxyz_1, g_yyyyyy_xxxyzz_1, g_yyyyyy_xxxzz_1, g_yyyyyy_xxxzzz_1, g_yyyyyy_xxyyy_1, g_yyyyyy_xxyyyy_1, g_yyyyyy_xxyyyz_1, g_yyyyyy_xxyyz_1, g_yyyyyy_xxyyzz_1, g_yyyyyy_xxyzz_1, g_yyyyyy_xxyzzz_1, g_yyyyyy_xxzzz_1, g_yyyyyy_xxzzzz_1, g_yyyyyy_xyyyy_1, g_yyyyyy_xyyyyy_1, g_yyyyyy_xyyyyz_1, g_yyyyyy_xyyyz_1, g_yyyyyy_xyyyzz_1, g_yyyyyy_xyyzz_1, g_yyyyyy_xyyzzz_1, g_yyyyyy_xyzzz_1, g_yyyyyy_xyzzzz_1, g_yyyyyy_xzzzz_1, g_yyyyyy_xzzzzz_1, g_yyyyyy_yyyyy_1, g_yyyyyy_yyyyyy_1, g_yyyyyy_yyyyyz_1, g_yyyyyy_yyyyz_1, g_yyyyyy_yyyyzz_1, g_yyyyyy_yyyzz_1, g_yyyyyy_yyyzzz_1, g_yyyyyy_yyzzz_1, g_yyyyyy_yyzzzz_1, g_yyyyyy_yzzzz_1, g_yyyyyy_yzzzzz_1, g_yyyyyy_zzzzz_1, g_yyyyyy_zzzzzz_1, g_yyyyyyz_xxxxxx_0, g_yyyyyyz_xxxxxy_0, g_yyyyyyz_xxxxxz_0, g_yyyyyyz_xxxxyy_0, g_yyyyyyz_xxxxyz_0, g_yyyyyyz_xxxxzz_0, g_yyyyyyz_xxxyyy_0, g_yyyyyyz_xxxyyz_0, g_yyyyyyz_xxxyzz_0, g_yyyyyyz_xxxzzz_0, g_yyyyyyz_xxyyyy_0, g_yyyyyyz_xxyyyz_0, g_yyyyyyz_xxyyzz_0, g_yyyyyyz_xxyzzz_0, g_yyyyyyz_xxzzzz_0, g_yyyyyyz_xyyyyy_0, g_yyyyyyz_xyyyyz_0, g_yyyyyyz_xyyyzz_0, g_yyyyyyz_xyyzzz_0, g_yyyyyyz_xyzzzz_0, g_yyyyyyz_xzzzzz_0, g_yyyyyyz_yyyyyy_0, g_yyyyyyz_yyyyyz_0, g_yyyyyyz_yyyyzz_0, g_yyyyyyz_yyyzzz_0, g_yyyyyyz_yyzzzz_0, g_yyyyyyz_yzzzzz_0, g_yyyyyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyyz_xxxxxx_0[i] = g_yyyyyy_xxxxxx_1[i] * pa_z[i];

        g_yyyyyyz_xxxxxy_0[i] = g_yyyyyy_xxxxxy_1[i] * pa_z[i];

        g_yyyyyyz_xxxxxz_0[i] = g_yyyyyy_xxxxx_1[i] * fe_0 + g_yyyyyy_xxxxxz_1[i] * pa_z[i];

        g_yyyyyyz_xxxxyy_0[i] = g_yyyyyy_xxxxyy_1[i] * pa_z[i];

        g_yyyyyyz_xxxxyz_0[i] = g_yyyyyy_xxxxy_1[i] * fe_0 + g_yyyyyy_xxxxyz_1[i] * pa_z[i];

        g_yyyyyyz_xxxxzz_0[i] = 2.0 * g_yyyyyy_xxxxz_1[i] * fe_0 + g_yyyyyy_xxxxzz_1[i] * pa_z[i];

        g_yyyyyyz_xxxyyy_0[i] = g_yyyyyy_xxxyyy_1[i] * pa_z[i];

        g_yyyyyyz_xxxyyz_0[i] = g_yyyyyy_xxxyy_1[i] * fe_0 + g_yyyyyy_xxxyyz_1[i] * pa_z[i];

        g_yyyyyyz_xxxyzz_0[i] = 2.0 * g_yyyyyy_xxxyz_1[i] * fe_0 + g_yyyyyy_xxxyzz_1[i] * pa_z[i];

        g_yyyyyyz_xxxzzz_0[i] = 3.0 * g_yyyyyy_xxxzz_1[i] * fe_0 + g_yyyyyy_xxxzzz_1[i] * pa_z[i];

        g_yyyyyyz_xxyyyy_0[i] = g_yyyyyy_xxyyyy_1[i] * pa_z[i];

        g_yyyyyyz_xxyyyz_0[i] = g_yyyyyy_xxyyy_1[i] * fe_0 + g_yyyyyy_xxyyyz_1[i] * pa_z[i];

        g_yyyyyyz_xxyyzz_0[i] = 2.0 * g_yyyyyy_xxyyz_1[i] * fe_0 + g_yyyyyy_xxyyzz_1[i] * pa_z[i];

        g_yyyyyyz_xxyzzz_0[i] = 3.0 * g_yyyyyy_xxyzz_1[i] * fe_0 + g_yyyyyy_xxyzzz_1[i] * pa_z[i];

        g_yyyyyyz_xxzzzz_0[i] = 4.0 * g_yyyyyy_xxzzz_1[i] * fe_0 + g_yyyyyy_xxzzzz_1[i] * pa_z[i];

        g_yyyyyyz_xyyyyy_0[i] = g_yyyyyy_xyyyyy_1[i] * pa_z[i];

        g_yyyyyyz_xyyyyz_0[i] = g_yyyyyy_xyyyy_1[i] * fe_0 + g_yyyyyy_xyyyyz_1[i] * pa_z[i];

        g_yyyyyyz_xyyyzz_0[i] = 2.0 * g_yyyyyy_xyyyz_1[i] * fe_0 + g_yyyyyy_xyyyzz_1[i] * pa_z[i];

        g_yyyyyyz_xyyzzz_0[i] = 3.0 * g_yyyyyy_xyyzz_1[i] * fe_0 + g_yyyyyy_xyyzzz_1[i] * pa_z[i];

        g_yyyyyyz_xyzzzz_0[i] = 4.0 * g_yyyyyy_xyzzz_1[i] * fe_0 + g_yyyyyy_xyzzzz_1[i] * pa_z[i];

        g_yyyyyyz_xzzzzz_0[i] = 5.0 * g_yyyyyy_xzzzz_1[i] * fe_0 + g_yyyyyy_xzzzzz_1[i] * pa_z[i];

        g_yyyyyyz_yyyyyy_0[i] = g_yyyyyy_yyyyyy_1[i] * pa_z[i];

        g_yyyyyyz_yyyyyz_0[i] = g_yyyyyy_yyyyy_1[i] * fe_0 + g_yyyyyy_yyyyyz_1[i] * pa_z[i];

        g_yyyyyyz_yyyyzz_0[i] = 2.0 * g_yyyyyy_yyyyz_1[i] * fe_0 + g_yyyyyy_yyyyzz_1[i] * pa_z[i];

        g_yyyyyyz_yyyzzz_0[i] = 3.0 * g_yyyyyy_yyyzz_1[i] * fe_0 + g_yyyyyy_yyyzzz_1[i] * pa_z[i];

        g_yyyyyyz_yyzzzz_0[i] = 4.0 * g_yyyyyy_yyzzz_1[i] * fe_0 + g_yyyyyy_yyzzzz_1[i] * pa_z[i];

        g_yyyyyyz_yzzzzz_0[i] = 5.0 * g_yyyyyy_yzzzz_1[i] * fe_0 + g_yyyyyy_yzzzzz_1[i] * pa_z[i];

        g_yyyyyyz_zzzzzz_0[i] = 6.0 * g_yyyyyy_zzzzz_1[i] * fe_0 + g_yyyyyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 840-868 components of targeted buffer : KI

    auto g_yyyyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 840);

    auto g_yyyyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 841);

    auto g_yyyyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 842);

    auto g_yyyyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 843);

    auto g_yyyyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 844);

    auto g_yyyyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 845);

    auto g_yyyyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 846);

    auto g_yyyyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 847);

    auto g_yyyyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 848);

    auto g_yyyyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 849);

    auto g_yyyyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 850);

    auto g_yyyyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 851);

    auto g_yyyyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 852);

    auto g_yyyyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 853);

    auto g_yyyyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 854);

    auto g_yyyyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 855);

    auto g_yyyyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 856);

    auto g_yyyyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 857);

    auto g_yyyyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 858);

    auto g_yyyyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 859);

    auto g_yyyyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 860);

    auto g_yyyyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 861);

    auto g_yyyyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 862);

    auto g_yyyyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 863);

    auto g_yyyyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 864);

    auto g_yyyyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 865);

    auto g_yyyyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 866);

    auto g_yyyyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 867);

    #pragma omp simd aligned(g_yyyyy_xxxxxy_0, g_yyyyy_xxxxxy_1, g_yyyyy_xxxxyy_0, g_yyyyy_xxxxyy_1, g_yyyyy_xxxyyy_0, g_yyyyy_xxxyyy_1, g_yyyyy_xxyyyy_0, g_yyyyy_xxyyyy_1, g_yyyyy_xyyyyy_0, g_yyyyy_xyyyyy_1, g_yyyyy_yyyyyy_0, g_yyyyy_yyyyyy_1, g_yyyyyz_xxxxxy_1, g_yyyyyz_xxxxyy_1, g_yyyyyz_xxxyyy_1, g_yyyyyz_xxyyyy_1, g_yyyyyz_xyyyyy_1, g_yyyyyz_yyyyyy_1, g_yyyyyzz_xxxxxx_0, g_yyyyyzz_xxxxxy_0, g_yyyyyzz_xxxxxz_0, g_yyyyyzz_xxxxyy_0, g_yyyyyzz_xxxxyz_0, g_yyyyyzz_xxxxzz_0, g_yyyyyzz_xxxyyy_0, g_yyyyyzz_xxxyyz_0, g_yyyyyzz_xxxyzz_0, g_yyyyyzz_xxxzzz_0, g_yyyyyzz_xxyyyy_0, g_yyyyyzz_xxyyyz_0, g_yyyyyzz_xxyyzz_0, g_yyyyyzz_xxyzzz_0, g_yyyyyzz_xxzzzz_0, g_yyyyyzz_xyyyyy_0, g_yyyyyzz_xyyyyz_0, g_yyyyyzz_xyyyzz_0, g_yyyyyzz_xyyzzz_0, g_yyyyyzz_xyzzzz_0, g_yyyyyzz_xzzzzz_0, g_yyyyyzz_yyyyyy_0, g_yyyyyzz_yyyyyz_0, g_yyyyyzz_yyyyzz_0, g_yyyyyzz_yyyzzz_0, g_yyyyyzz_yyzzzz_0, g_yyyyyzz_yzzzzz_0, g_yyyyyzz_zzzzzz_0, g_yyyyzz_xxxxxx_1, g_yyyyzz_xxxxxz_1, g_yyyyzz_xxxxyz_1, g_yyyyzz_xxxxz_1, g_yyyyzz_xxxxzz_1, g_yyyyzz_xxxyyz_1, g_yyyyzz_xxxyz_1, g_yyyyzz_xxxyzz_1, g_yyyyzz_xxxzz_1, g_yyyyzz_xxxzzz_1, g_yyyyzz_xxyyyz_1, g_yyyyzz_xxyyz_1, g_yyyyzz_xxyyzz_1, g_yyyyzz_xxyzz_1, g_yyyyzz_xxyzzz_1, g_yyyyzz_xxzzz_1, g_yyyyzz_xxzzzz_1, g_yyyyzz_xyyyyz_1, g_yyyyzz_xyyyz_1, g_yyyyzz_xyyyzz_1, g_yyyyzz_xyyzz_1, g_yyyyzz_xyyzzz_1, g_yyyyzz_xyzzz_1, g_yyyyzz_xyzzzz_1, g_yyyyzz_xzzzz_1, g_yyyyzz_xzzzzz_1, g_yyyyzz_yyyyyz_1, g_yyyyzz_yyyyz_1, g_yyyyzz_yyyyzz_1, g_yyyyzz_yyyzz_1, g_yyyyzz_yyyzzz_1, g_yyyyzz_yyzzz_1, g_yyyyzz_yyzzzz_1, g_yyyyzz_yzzzz_1, g_yyyyzz_yzzzzz_1, g_yyyyzz_zzzzz_1, g_yyyyzz_zzzzzz_1, g_yyyzz_xxxxxx_0, g_yyyzz_xxxxxx_1, g_yyyzz_xxxxxz_0, g_yyyzz_xxxxxz_1, g_yyyzz_xxxxyz_0, g_yyyzz_xxxxyz_1, g_yyyzz_xxxxzz_0, g_yyyzz_xxxxzz_1, g_yyyzz_xxxyyz_0, g_yyyzz_xxxyyz_1, g_yyyzz_xxxyzz_0, g_yyyzz_xxxyzz_1, g_yyyzz_xxxzzz_0, g_yyyzz_xxxzzz_1, g_yyyzz_xxyyyz_0, g_yyyzz_xxyyyz_1, g_yyyzz_xxyyzz_0, g_yyyzz_xxyyzz_1, g_yyyzz_xxyzzz_0, g_yyyzz_xxyzzz_1, g_yyyzz_xxzzzz_0, g_yyyzz_xxzzzz_1, g_yyyzz_xyyyyz_0, g_yyyzz_xyyyyz_1, g_yyyzz_xyyyzz_0, g_yyyzz_xyyyzz_1, g_yyyzz_xyyzzz_0, g_yyyzz_xyyzzz_1, g_yyyzz_xyzzzz_0, g_yyyzz_xyzzzz_1, g_yyyzz_xzzzzz_0, g_yyyzz_xzzzzz_1, g_yyyzz_yyyyyz_0, g_yyyzz_yyyyyz_1, g_yyyzz_yyyyzz_0, g_yyyzz_yyyyzz_1, g_yyyzz_yyyzzz_0, g_yyyzz_yyyzzz_1, g_yyyzz_yyzzzz_0, g_yyyzz_yyzzzz_1, g_yyyzz_yzzzzz_0, g_yyyzz_yzzzzz_1, g_yyyzz_zzzzzz_0, g_yyyzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyzz_xxxxxx_0[i] = 4.0 * g_yyyzz_xxxxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxxxx_1[i] * fz_be_0 + g_yyyyzz_xxxxxx_1[i] * pa_y[i];

        g_yyyyyzz_xxxxxy_0[i] = g_yyyyy_xxxxxy_0[i] * fbe_0 - g_yyyyy_xxxxxy_1[i] * fz_be_0 + g_yyyyyz_xxxxxy_1[i] * pa_z[i];

        g_yyyyyzz_xxxxxz_0[i] = 4.0 * g_yyyzz_xxxxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxxxz_1[i] * fz_be_0 + g_yyyyzz_xxxxxz_1[i] * pa_y[i];

        g_yyyyyzz_xxxxyy_0[i] = g_yyyyy_xxxxyy_0[i] * fbe_0 - g_yyyyy_xxxxyy_1[i] * fz_be_0 + g_yyyyyz_xxxxyy_1[i] * pa_z[i];

        g_yyyyyzz_xxxxyz_0[i] = 4.0 * g_yyyzz_xxxxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxxyz_1[i] * fz_be_0 + g_yyyyzz_xxxxz_1[i] * fe_0 + g_yyyyzz_xxxxyz_1[i] * pa_y[i];

        g_yyyyyzz_xxxxzz_0[i] = 4.0 * g_yyyzz_xxxxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxxzz_1[i] * fz_be_0 + g_yyyyzz_xxxxzz_1[i] * pa_y[i];

        g_yyyyyzz_xxxyyy_0[i] = g_yyyyy_xxxyyy_0[i] * fbe_0 - g_yyyyy_xxxyyy_1[i] * fz_be_0 + g_yyyyyz_xxxyyy_1[i] * pa_z[i];

        g_yyyyyzz_xxxyyz_0[i] = 4.0 * g_yyyzz_xxxyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_xxxyz_1[i] * fe_0 + g_yyyyzz_xxxyyz_1[i] * pa_y[i];

        g_yyyyyzz_xxxyzz_0[i] = 4.0 * g_yyyzz_xxxyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxyzz_1[i] * fz_be_0 + g_yyyyzz_xxxzz_1[i] * fe_0 + g_yyyyzz_xxxyzz_1[i] * pa_y[i];

        g_yyyyyzz_xxxzzz_0[i] = 4.0 * g_yyyzz_xxxzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxzzz_1[i] * fz_be_0 + g_yyyyzz_xxxzzz_1[i] * pa_y[i];

        g_yyyyyzz_xxyyyy_0[i] = g_yyyyy_xxyyyy_0[i] * fbe_0 - g_yyyyy_xxyyyy_1[i] * fz_be_0 + g_yyyyyz_xxyyyy_1[i] * pa_z[i];

        g_yyyyyzz_xxyyyz_0[i] = 4.0 * g_yyyzz_xxyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_xxyyz_1[i] * fe_0 + g_yyyyzz_xxyyyz_1[i] * pa_y[i];

        g_yyyyyzz_xxyyzz_0[i] = 4.0 * g_yyyzz_xxyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_xxyzz_1[i] * fe_0 + g_yyyyzz_xxyyzz_1[i] * pa_y[i];

        g_yyyyyzz_xxyzzz_0[i] = 4.0 * g_yyyzz_xxyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxyzzz_1[i] * fz_be_0 + g_yyyyzz_xxzzz_1[i] * fe_0 + g_yyyyzz_xxyzzz_1[i] * pa_y[i];

        g_yyyyyzz_xxzzzz_0[i] = 4.0 * g_yyyzz_xxzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxzzzz_1[i] * fz_be_0 + g_yyyyzz_xxzzzz_1[i] * pa_y[i];

        g_yyyyyzz_xyyyyy_0[i] = g_yyyyy_xyyyyy_0[i] * fbe_0 - g_yyyyy_xyyyyy_1[i] * fz_be_0 + g_yyyyyz_xyyyyy_1[i] * pa_z[i];

        g_yyyyyzz_xyyyyz_0[i] = 4.0 * g_yyyzz_xyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_xyyyz_1[i] * fe_0 + g_yyyyzz_xyyyyz_1[i] * pa_y[i];

        g_yyyyyzz_xyyyzz_0[i] = 4.0 * g_yyyzz_xyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_xyyzz_1[i] * fe_0 + g_yyyyzz_xyyyzz_1[i] * pa_y[i];

        g_yyyyyzz_xyyzzz_0[i] = 4.0 * g_yyyzz_xyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_xyzzz_1[i] * fe_0 + g_yyyyzz_xyyzzz_1[i] * pa_y[i];

        g_yyyyyzz_xyzzzz_0[i] = 4.0 * g_yyyzz_xyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyzzzz_1[i] * fz_be_0 + g_yyyyzz_xzzzz_1[i] * fe_0 + g_yyyyzz_xyzzzz_1[i] * pa_y[i];

        g_yyyyyzz_xzzzzz_0[i] = 4.0 * g_yyyzz_xzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xzzzzz_1[i] * fz_be_0 + g_yyyyzz_xzzzzz_1[i] * pa_y[i];

        g_yyyyyzz_yyyyyy_0[i] = g_yyyyy_yyyyyy_0[i] * fbe_0 - g_yyyyy_yyyyyy_1[i] * fz_be_0 + g_yyyyyz_yyyyyy_1[i] * pa_z[i];

        g_yyyyyzz_yyyyyz_0[i] = 4.0 * g_yyyzz_yyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_yyyyz_1[i] * fe_0 + g_yyyyzz_yyyyyz_1[i] * pa_y[i];

        g_yyyyyzz_yyyyzz_0[i] = 4.0 * g_yyyzz_yyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_yyyzz_1[i] * fe_0 + g_yyyyzz_yyyyzz_1[i] * pa_y[i];

        g_yyyyyzz_yyyzzz_0[i] = 4.0 * g_yyyzz_yyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_yyzzz_1[i] * fe_0 + g_yyyyzz_yyyzzz_1[i] * pa_y[i];

        g_yyyyyzz_yyzzzz_0[i] = 4.0 * g_yyyzz_yyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_yzzzz_1[i] * fe_0 + g_yyyyzz_yyzzzz_1[i] * pa_y[i];

        g_yyyyyzz_yzzzzz_0[i] = 4.0 * g_yyyzz_yzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yzzzzz_1[i] * fz_be_0 + g_yyyyzz_zzzzz_1[i] * fe_0 + g_yyyyzz_yzzzzz_1[i] * pa_y[i];

        g_yyyyyzz_zzzzzz_0[i] = 4.0 * g_yyyzz_zzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_zzzzzz_1[i] * fz_be_0 + g_yyyyzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 868-896 components of targeted buffer : KI

    auto g_yyyyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 868);

    auto g_yyyyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 869);

    auto g_yyyyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 870);

    auto g_yyyyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 871);

    auto g_yyyyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 872);

    auto g_yyyyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 873);

    auto g_yyyyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 874);

    auto g_yyyyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 875);

    auto g_yyyyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 876);

    auto g_yyyyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 877);

    auto g_yyyyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 878);

    auto g_yyyyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 879);

    auto g_yyyyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 880);

    auto g_yyyyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 881);

    auto g_yyyyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 882);

    auto g_yyyyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 883);

    auto g_yyyyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 884);

    auto g_yyyyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 885);

    auto g_yyyyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 886);

    auto g_yyyyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 887);

    auto g_yyyyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 888);

    auto g_yyyyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 889);

    auto g_yyyyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 890);

    auto g_yyyyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 891);

    auto g_yyyyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 892);

    auto g_yyyyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 893);

    auto g_yyyyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 894);

    auto g_yyyyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 895);

    #pragma omp simd aligned(g_yyyyz_xxxxxy_0, g_yyyyz_xxxxxy_1, g_yyyyz_xxxxyy_0, g_yyyyz_xxxxyy_1, g_yyyyz_xxxyyy_0, g_yyyyz_xxxyyy_1, g_yyyyz_xxyyyy_0, g_yyyyz_xxyyyy_1, g_yyyyz_xyyyyy_0, g_yyyyz_xyyyyy_1, g_yyyyz_yyyyyy_0, g_yyyyz_yyyyyy_1, g_yyyyzz_xxxxxy_1, g_yyyyzz_xxxxyy_1, g_yyyyzz_xxxyyy_1, g_yyyyzz_xxyyyy_1, g_yyyyzz_xyyyyy_1, g_yyyyzz_yyyyyy_1, g_yyyyzzz_xxxxxx_0, g_yyyyzzz_xxxxxy_0, g_yyyyzzz_xxxxxz_0, g_yyyyzzz_xxxxyy_0, g_yyyyzzz_xxxxyz_0, g_yyyyzzz_xxxxzz_0, g_yyyyzzz_xxxyyy_0, g_yyyyzzz_xxxyyz_0, g_yyyyzzz_xxxyzz_0, g_yyyyzzz_xxxzzz_0, g_yyyyzzz_xxyyyy_0, g_yyyyzzz_xxyyyz_0, g_yyyyzzz_xxyyzz_0, g_yyyyzzz_xxyzzz_0, g_yyyyzzz_xxzzzz_0, g_yyyyzzz_xyyyyy_0, g_yyyyzzz_xyyyyz_0, g_yyyyzzz_xyyyzz_0, g_yyyyzzz_xyyzzz_0, g_yyyyzzz_xyzzzz_0, g_yyyyzzz_xzzzzz_0, g_yyyyzzz_yyyyyy_0, g_yyyyzzz_yyyyyz_0, g_yyyyzzz_yyyyzz_0, g_yyyyzzz_yyyzzz_0, g_yyyyzzz_yyzzzz_0, g_yyyyzzz_yzzzzz_0, g_yyyyzzz_zzzzzz_0, g_yyyzzz_xxxxxx_1, g_yyyzzz_xxxxxz_1, g_yyyzzz_xxxxyz_1, g_yyyzzz_xxxxz_1, g_yyyzzz_xxxxzz_1, g_yyyzzz_xxxyyz_1, g_yyyzzz_xxxyz_1, g_yyyzzz_xxxyzz_1, g_yyyzzz_xxxzz_1, g_yyyzzz_xxxzzz_1, g_yyyzzz_xxyyyz_1, g_yyyzzz_xxyyz_1, g_yyyzzz_xxyyzz_1, g_yyyzzz_xxyzz_1, g_yyyzzz_xxyzzz_1, g_yyyzzz_xxzzz_1, g_yyyzzz_xxzzzz_1, g_yyyzzz_xyyyyz_1, g_yyyzzz_xyyyz_1, g_yyyzzz_xyyyzz_1, g_yyyzzz_xyyzz_1, g_yyyzzz_xyyzzz_1, g_yyyzzz_xyzzz_1, g_yyyzzz_xyzzzz_1, g_yyyzzz_xzzzz_1, g_yyyzzz_xzzzzz_1, g_yyyzzz_yyyyyz_1, g_yyyzzz_yyyyz_1, g_yyyzzz_yyyyzz_1, g_yyyzzz_yyyzz_1, g_yyyzzz_yyyzzz_1, g_yyyzzz_yyzzz_1, g_yyyzzz_yyzzzz_1, g_yyyzzz_yzzzz_1, g_yyyzzz_yzzzzz_1, g_yyyzzz_zzzzz_1, g_yyyzzz_zzzzzz_1, g_yyzzz_xxxxxx_0, g_yyzzz_xxxxxx_1, g_yyzzz_xxxxxz_0, g_yyzzz_xxxxxz_1, g_yyzzz_xxxxyz_0, g_yyzzz_xxxxyz_1, g_yyzzz_xxxxzz_0, g_yyzzz_xxxxzz_1, g_yyzzz_xxxyyz_0, g_yyzzz_xxxyyz_1, g_yyzzz_xxxyzz_0, g_yyzzz_xxxyzz_1, g_yyzzz_xxxzzz_0, g_yyzzz_xxxzzz_1, g_yyzzz_xxyyyz_0, g_yyzzz_xxyyyz_1, g_yyzzz_xxyyzz_0, g_yyzzz_xxyyzz_1, g_yyzzz_xxyzzz_0, g_yyzzz_xxyzzz_1, g_yyzzz_xxzzzz_0, g_yyzzz_xxzzzz_1, g_yyzzz_xyyyyz_0, g_yyzzz_xyyyyz_1, g_yyzzz_xyyyzz_0, g_yyzzz_xyyyzz_1, g_yyzzz_xyyzzz_0, g_yyzzz_xyyzzz_1, g_yyzzz_xyzzzz_0, g_yyzzz_xyzzzz_1, g_yyzzz_xzzzzz_0, g_yyzzz_xzzzzz_1, g_yyzzz_yyyyyz_0, g_yyzzz_yyyyyz_1, g_yyzzz_yyyyzz_0, g_yyzzz_yyyyzz_1, g_yyzzz_yyyzzz_0, g_yyzzz_yyyzzz_1, g_yyzzz_yyzzzz_0, g_yyzzz_yyzzzz_1, g_yyzzz_yzzzzz_0, g_yyzzz_yzzzzz_1, g_yyzzz_zzzzzz_0, g_yyzzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzzz_xxxxxx_0[i] = 3.0 * g_yyzzz_xxxxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxxxx_1[i] * fz_be_0 + g_yyyzzz_xxxxxx_1[i] * pa_y[i];

        g_yyyyzzz_xxxxxy_0[i] = 2.0 * g_yyyyz_xxxxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxxxxy_1[i] * fz_be_0 + g_yyyyzz_xxxxxy_1[i] * pa_z[i];

        g_yyyyzzz_xxxxxz_0[i] = 3.0 * g_yyzzz_xxxxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxxxz_1[i] * fz_be_0 + g_yyyzzz_xxxxxz_1[i] * pa_y[i];

        g_yyyyzzz_xxxxyy_0[i] = 2.0 * g_yyyyz_xxxxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxxxyy_1[i] * fz_be_0 + g_yyyyzz_xxxxyy_1[i] * pa_z[i];

        g_yyyyzzz_xxxxyz_0[i] = 3.0 * g_yyzzz_xxxxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxxyz_1[i] * fz_be_0 + g_yyyzzz_xxxxz_1[i] * fe_0 + g_yyyzzz_xxxxyz_1[i] * pa_y[i];

        g_yyyyzzz_xxxxzz_0[i] = 3.0 * g_yyzzz_xxxxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxxzz_1[i] * fz_be_0 + g_yyyzzz_xxxxzz_1[i] * pa_y[i];

        g_yyyyzzz_xxxyyy_0[i] = 2.0 * g_yyyyz_xxxyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxxyyy_1[i] * fz_be_0 + g_yyyyzz_xxxyyy_1[i] * pa_z[i];

        g_yyyyzzz_xxxyyz_0[i] = 3.0 * g_yyzzz_xxxyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_xxxyz_1[i] * fe_0 + g_yyyzzz_xxxyyz_1[i] * pa_y[i];

        g_yyyyzzz_xxxyzz_0[i] = 3.0 * g_yyzzz_xxxyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxyzz_1[i] * fz_be_0 + g_yyyzzz_xxxzz_1[i] * fe_0 + g_yyyzzz_xxxyzz_1[i] * pa_y[i];

        g_yyyyzzz_xxxzzz_0[i] = 3.0 * g_yyzzz_xxxzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxzzz_1[i] * fz_be_0 + g_yyyzzz_xxxzzz_1[i] * pa_y[i];

        g_yyyyzzz_xxyyyy_0[i] = 2.0 * g_yyyyz_xxyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxyyyy_1[i] * fz_be_0 + g_yyyyzz_xxyyyy_1[i] * pa_z[i];

        g_yyyyzzz_xxyyyz_0[i] = 3.0 * g_yyzzz_xxyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_xxyyz_1[i] * fe_0 + g_yyyzzz_xxyyyz_1[i] * pa_y[i];

        g_yyyyzzz_xxyyzz_0[i] = 3.0 * g_yyzzz_xxyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_xxyzz_1[i] * fe_0 + g_yyyzzz_xxyyzz_1[i] * pa_y[i];

        g_yyyyzzz_xxyzzz_0[i] = 3.0 * g_yyzzz_xxyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxyzzz_1[i] * fz_be_0 + g_yyyzzz_xxzzz_1[i] * fe_0 + g_yyyzzz_xxyzzz_1[i] * pa_y[i];

        g_yyyyzzz_xxzzzz_0[i] = 3.0 * g_yyzzz_xxzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxzzzz_1[i] * fz_be_0 + g_yyyzzz_xxzzzz_1[i] * pa_y[i];

        g_yyyyzzz_xyyyyy_0[i] = 2.0 * g_yyyyz_xyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xyyyyy_1[i] * fz_be_0 + g_yyyyzz_xyyyyy_1[i] * pa_z[i];

        g_yyyyzzz_xyyyyz_0[i] = 3.0 * g_yyzzz_xyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_xyyyz_1[i] * fe_0 + g_yyyzzz_xyyyyz_1[i] * pa_y[i];

        g_yyyyzzz_xyyyzz_0[i] = 3.0 * g_yyzzz_xyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_xyyzz_1[i] * fe_0 + g_yyyzzz_xyyyzz_1[i] * pa_y[i];

        g_yyyyzzz_xyyzzz_0[i] = 3.0 * g_yyzzz_xyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_xyzzz_1[i] * fe_0 + g_yyyzzz_xyyzzz_1[i] * pa_y[i];

        g_yyyyzzz_xyzzzz_0[i] = 3.0 * g_yyzzz_xyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyzzzz_1[i] * fz_be_0 + g_yyyzzz_xzzzz_1[i] * fe_0 + g_yyyzzz_xyzzzz_1[i] * pa_y[i];

        g_yyyyzzz_xzzzzz_0[i] = 3.0 * g_yyzzz_xzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xzzzzz_1[i] * fz_be_0 + g_yyyzzz_xzzzzz_1[i] * pa_y[i];

        g_yyyyzzz_yyyyyy_0[i] = 2.0 * g_yyyyz_yyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_yyyyyy_1[i] * fz_be_0 + g_yyyyzz_yyyyyy_1[i] * pa_z[i];

        g_yyyyzzz_yyyyyz_0[i] = 3.0 * g_yyzzz_yyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_yyyyz_1[i] * fe_0 + g_yyyzzz_yyyyyz_1[i] * pa_y[i];

        g_yyyyzzz_yyyyzz_0[i] = 3.0 * g_yyzzz_yyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_yyyzz_1[i] * fe_0 + g_yyyzzz_yyyyzz_1[i] * pa_y[i];

        g_yyyyzzz_yyyzzz_0[i] = 3.0 * g_yyzzz_yyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_yyzzz_1[i] * fe_0 + g_yyyzzz_yyyzzz_1[i] * pa_y[i];

        g_yyyyzzz_yyzzzz_0[i] = 3.0 * g_yyzzz_yyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_yzzzz_1[i] * fe_0 + g_yyyzzz_yyzzzz_1[i] * pa_y[i];

        g_yyyyzzz_yzzzzz_0[i] = 3.0 * g_yyzzz_yzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yzzzzz_1[i] * fz_be_0 + g_yyyzzz_zzzzz_1[i] * fe_0 + g_yyyzzz_yzzzzz_1[i] * pa_y[i];

        g_yyyyzzz_zzzzzz_0[i] = 3.0 * g_yyzzz_zzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_zzzzzz_1[i] * fz_be_0 + g_yyyzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 896-924 components of targeted buffer : KI

    auto g_yyyzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 896);

    auto g_yyyzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 897);

    auto g_yyyzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 898);

    auto g_yyyzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 899);

    auto g_yyyzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 900);

    auto g_yyyzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 901);

    auto g_yyyzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 902);

    auto g_yyyzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 903);

    auto g_yyyzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 904);

    auto g_yyyzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 905);

    auto g_yyyzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 906);

    auto g_yyyzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 907);

    auto g_yyyzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 908);

    auto g_yyyzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 909);

    auto g_yyyzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 910);

    auto g_yyyzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 911);

    auto g_yyyzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 912);

    auto g_yyyzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 913);

    auto g_yyyzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 914);

    auto g_yyyzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 915);

    auto g_yyyzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 916);

    auto g_yyyzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 917);

    auto g_yyyzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 918);

    auto g_yyyzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 919);

    auto g_yyyzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 920);

    auto g_yyyzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 921);

    auto g_yyyzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 922);

    auto g_yyyzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 923);

    #pragma omp simd aligned(g_yyyzz_xxxxxy_0, g_yyyzz_xxxxxy_1, g_yyyzz_xxxxyy_0, g_yyyzz_xxxxyy_1, g_yyyzz_xxxyyy_0, g_yyyzz_xxxyyy_1, g_yyyzz_xxyyyy_0, g_yyyzz_xxyyyy_1, g_yyyzz_xyyyyy_0, g_yyyzz_xyyyyy_1, g_yyyzz_yyyyyy_0, g_yyyzz_yyyyyy_1, g_yyyzzz_xxxxxy_1, g_yyyzzz_xxxxyy_1, g_yyyzzz_xxxyyy_1, g_yyyzzz_xxyyyy_1, g_yyyzzz_xyyyyy_1, g_yyyzzz_yyyyyy_1, g_yyyzzzz_xxxxxx_0, g_yyyzzzz_xxxxxy_0, g_yyyzzzz_xxxxxz_0, g_yyyzzzz_xxxxyy_0, g_yyyzzzz_xxxxyz_0, g_yyyzzzz_xxxxzz_0, g_yyyzzzz_xxxyyy_0, g_yyyzzzz_xxxyyz_0, g_yyyzzzz_xxxyzz_0, g_yyyzzzz_xxxzzz_0, g_yyyzzzz_xxyyyy_0, g_yyyzzzz_xxyyyz_0, g_yyyzzzz_xxyyzz_0, g_yyyzzzz_xxyzzz_0, g_yyyzzzz_xxzzzz_0, g_yyyzzzz_xyyyyy_0, g_yyyzzzz_xyyyyz_0, g_yyyzzzz_xyyyzz_0, g_yyyzzzz_xyyzzz_0, g_yyyzzzz_xyzzzz_0, g_yyyzzzz_xzzzzz_0, g_yyyzzzz_yyyyyy_0, g_yyyzzzz_yyyyyz_0, g_yyyzzzz_yyyyzz_0, g_yyyzzzz_yyyzzz_0, g_yyyzzzz_yyzzzz_0, g_yyyzzzz_yzzzzz_0, g_yyyzzzz_zzzzzz_0, g_yyzzzz_xxxxxx_1, g_yyzzzz_xxxxxz_1, g_yyzzzz_xxxxyz_1, g_yyzzzz_xxxxz_1, g_yyzzzz_xxxxzz_1, g_yyzzzz_xxxyyz_1, g_yyzzzz_xxxyz_1, g_yyzzzz_xxxyzz_1, g_yyzzzz_xxxzz_1, g_yyzzzz_xxxzzz_1, g_yyzzzz_xxyyyz_1, g_yyzzzz_xxyyz_1, g_yyzzzz_xxyyzz_1, g_yyzzzz_xxyzz_1, g_yyzzzz_xxyzzz_1, g_yyzzzz_xxzzz_1, g_yyzzzz_xxzzzz_1, g_yyzzzz_xyyyyz_1, g_yyzzzz_xyyyz_1, g_yyzzzz_xyyyzz_1, g_yyzzzz_xyyzz_1, g_yyzzzz_xyyzzz_1, g_yyzzzz_xyzzz_1, g_yyzzzz_xyzzzz_1, g_yyzzzz_xzzzz_1, g_yyzzzz_xzzzzz_1, g_yyzzzz_yyyyyz_1, g_yyzzzz_yyyyz_1, g_yyzzzz_yyyyzz_1, g_yyzzzz_yyyzz_1, g_yyzzzz_yyyzzz_1, g_yyzzzz_yyzzz_1, g_yyzzzz_yyzzzz_1, g_yyzzzz_yzzzz_1, g_yyzzzz_yzzzzz_1, g_yyzzzz_zzzzz_1, g_yyzzzz_zzzzzz_1, g_yzzzz_xxxxxx_0, g_yzzzz_xxxxxx_1, g_yzzzz_xxxxxz_0, g_yzzzz_xxxxxz_1, g_yzzzz_xxxxyz_0, g_yzzzz_xxxxyz_1, g_yzzzz_xxxxzz_0, g_yzzzz_xxxxzz_1, g_yzzzz_xxxyyz_0, g_yzzzz_xxxyyz_1, g_yzzzz_xxxyzz_0, g_yzzzz_xxxyzz_1, g_yzzzz_xxxzzz_0, g_yzzzz_xxxzzz_1, g_yzzzz_xxyyyz_0, g_yzzzz_xxyyyz_1, g_yzzzz_xxyyzz_0, g_yzzzz_xxyyzz_1, g_yzzzz_xxyzzz_0, g_yzzzz_xxyzzz_1, g_yzzzz_xxzzzz_0, g_yzzzz_xxzzzz_1, g_yzzzz_xyyyyz_0, g_yzzzz_xyyyyz_1, g_yzzzz_xyyyzz_0, g_yzzzz_xyyyzz_1, g_yzzzz_xyyzzz_0, g_yzzzz_xyyzzz_1, g_yzzzz_xyzzzz_0, g_yzzzz_xyzzzz_1, g_yzzzz_xzzzzz_0, g_yzzzz_xzzzzz_1, g_yzzzz_yyyyyz_0, g_yzzzz_yyyyyz_1, g_yzzzz_yyyyzz_0, g_yzzzz_yyyyzz_1, g_yzzzz_yyyzzz_0, g_yzzzz_yyyzzz_1, g_yzzzz_yyzzzz_0, g_yzzzz_yyzzzz_1, g_yzzzz_yzzzzz_0, g_yzzzz_yzzzzz_1, g_yzzzz_zzzzzz_0, g_yzzzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzzz_xxxxxx_0[i] = 2.0 * g_yzzzz_xxxxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxxxx_1[i] * fz_be_0 + g_yyzzzz_xxxxxx_1[i] * pa_y[i];

        g_yyyzzzz_xxxxxy_0[i] = 3.0 * g_yyyzz_xxxxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxxxxy_1[i] * fz_be_0 + g_yyyzzz_xxxxxy_1[i] * pa_z[i];

        g_yyyzzzz_xxxxxz_0[i] = 2.0 * g_yzzzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxxxz_1[i] * fz_be_0 + g_yyzzzz_xxxxxz_1[i] * pa_y[i];

        g_yyyzzzz_xxxxyy_0[i] = 3.0 * g_yyyzz_xxxxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxxxyy_1[i] * fz_be_0 + g_yyyzzz_xxxxyy_1[i] * pa_z[i];

        g_yyyzzzz_xxxxyz_0[i] = 2.0 * g_yzzzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxxyz_1[i] * fz_be_0 + g_yyzzzz_xxxxz_1[i] * fe_0 + g_yyzzzz_xxxxyz_1[i] * pa_y[i];

        g_yyyzzzz_xxxxzz_0[i] = 2.0 * g_yzzzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxxzz_1[i] * fz_be_0 + g_yyzzzz_xxxxzz_1[i] * pa_y[i];

        g_yyyzzzz_xxxyyy_0[i] = 3.0 * g_yyyzz_xxxyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxxyyy_1[i] * fz_be_0 + g_yyyzzz_xxxyyy_1[i] * pa_z[i];

        g_yyyzzzz_xxxyyz_0[i] = 2.0 * g_yzzzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_xxxyz_1[i] * fe_0 + g_yyzzzz_xxxyyz_1[i] * pa_y[i];

        g_yyyzzzz_xxxyzz_0[i] = 2.0 * g_yzzzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxyzz_1[i] * fz_be_0 + g_yyzzzz_xxxzz_1[i] * fe_0 + g_yyzzzz_xxxyzz_1[i] * pa_y[i];

        g_yyyzzzz_xxxzzz_0[i] = 2.0 * g_yzzzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxzzz_1[i] * fz_be_0 + g_yyzzzz_xxxzzz_1[i] * pa_y[i];

        g_yyyzzzz_xxyyyy_0[i] = 3.0 * g_yyyzz_xxyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxyyyy_1[i] * fz_be_0 + g_yyyzzz_xxyyyy_1[i] * pa_z[i];

        g_yyyzzzz_xxyyyz_0[i] = 2.0 * g_yzzzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_xxyyz_1[i] * fe_0 + g_yyzzzz_xxyyyz_1[i] * pa_y[i];

        g_yyyzzzz_xxyyzz_0[i] = 2.0 * g_yzzzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_xxyzz_1[i] * fe_0 + g_yyzzzz_xxyyzz_1[i] * pa_y[i];

        g_yyyzzzz_xxyzzz_0[i] = 2.0 * g_yzzzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxyzzz_1[i] * fz_be_0 + g_yyzzzz_xxzzz_1[i] * fe_0 + g_yyzzzz_xxyzzz_1[i] * pa_y[i];

        g_yyyzzzz_xxzzzz_0[i] = 2.0 * g_yzzzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxzzzz_1[i] * fz_be_0 + g_yyzzzz_xxzzzz_1[i] * pa_y[i];

        g_yyyzzzz_xyyyyy_0[i] = 3.0 * g_yyyzz_xyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xyyyyy_1[i] * fz_be_0 + g_yyyzzz_xyyyyy_1[i] * pa_z[i];

        g_yyyzzzz_xyyyyz_0[i] = 2.0 * g_yzzzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_xyyyz_1[i] * fe_0 + g_yyzzzz_xyyyyz_1[i] * pa_y[i];

        g_yyyzzzz_xyyyzz_0[i] = 2.0 * g_yzzzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_xyyzz_1[i] * fe_0 + g_yyzzzz_xyyyzz_1[i] * pa_y[i];

        g_yyyzzzz_xyyzzz_0[i] = 2.0 * g_yzzzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_xyzzz_1[i] * fe_0 + g_yyzzzz_xyyzzz_1[i] * pa_y[i];

        g_yyyzzzz_xyzzzz_0[i] = 2.0 * g_yzzzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyzzzz_1[i] * fz_be_0 + g_yyzzzz_xzzzz_1[i] * fe_0 + g_yyzzzz_xyzzzz_1[i] * pa_y[i];

        g_yyyzzzz_xzzzzz_0[i] = 2.0 * g_yzzzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xzzzzz_1[i] * fz_be_0 + g_yyzzzz_xzzzzz_1[i] * pa_y[i];

        g_yyyzzzz_yyyyyy_0[i] = 3.0 * g_yyyzz_yyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_yyyyyy_1[i] * fz_be_0 + g_yyyzzz_yyyyyy_1[i] * pa_z[i];

        g_yyyzzzz_yyyyyz_0[i] = 2.0 * g_yzzzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_yyyyz_1[i] * fe_0 + g_yyzzzz_yyyyyz_1[i] * pa_y[i];

        g_yyyzzzz_yyyyzz_0[i] = 2.0 * g_yzzzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_yyyzz_1[i] * fe_0 + g_yyzzzz_yyyyzz_1[i] * pa_y[i];

        g_yyyzzzz_yyyzzz_0[i] = 2.0 * g_yzzzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_yyzzz_1[i] * fe_0 + g_yyzzzz_yyyzzz_1[i] * pa_y[i];

        g_yyyzzzz_yyzzzz_0[i] = 2.0 * g_yzzzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_yzzzz_1[i] * fe_0 + g_yyzzzz_yyzzzz_1[i] * pa_y[i];

        g_yyyzzzz_yzzzzz_0[i] = 2.0 * g_yzzzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yzzzzz_1[i] * fz_be_0 + g_yyzzzz_zzzzz_1[i] * fe_0 + g_yyzzzz_yzzzzz_1[i] * pa_y[i];

        g_yyyzzzz_zzzzzz_0[i] = 2.0 * g_yzzzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_zzzzzz_1[i] * fz_be_0 + g_yyzzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 924-952 components of targeted buffer : KI

    auto g_yyzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 924);

    auto g_yyzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 925);

    auto g_yyzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 926);

    auto g_yyzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 927);

    auto g_yyzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 928);

    auto g_yyzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 929);

    auto g_yyzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 930);

    auto g_yyzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 931);

    auto g_yyzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 932);

    auto g_yyzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 933);

    auto g_yyzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 934);

    auto g_yyzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 935);

    auto g_yyzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 936);

    auto g_yyzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 937);

    auto g_yyzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 938);

    auto g_yyzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 939);

    auto g_yyzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 940);

    auto g_yyzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 941);

    auto g_yyzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 942);

    auto g_yyzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 943);

    auto g_yyzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 944);

    auto g_yyzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 945);

    auto g_yyzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 946);

    auto g_yyzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 947);

    auto g_yyzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 948);

    auto g_yyzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 949);

    auto g_yyzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 950);

    auto g_yyzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 951);

    #pragma omp simd aligned(g_yyzzz_xxxxxy_0, g_yyzzz_xxxxxy_1, g_yyzzz_xxxxyy_0, g_yyzzz_xxxxyy_1, g_yyzzz_xxxyyy_0, g_yyzzz_xxxyyy_1, g_yyzzz_xxyyyy_0, g_yyzzz_xxyyyy_1, g_yyzzz_xyyyyy_0, g_yyzzz_xyyyyy_1, g_yyzzz_yyyyyy_0, g_yyzzz_yyyyyy_1, g_yyzzzz_xxxxxy_1, g_yyzzzz_xxxxyy_1, g_yyzzzz_xxxyyy_1, g_yyzzzz_xxyyyy_1, g_yyzzzz_xyyyyy_1, g_yyzzzz_yyyyyy_1, g_yyzzzzz_xxxxxx_0, g_yyzzzzz_xxxxxy_0, g_yyzzzzz_xxxxxz_0, g_yyzzzzz_xxxxyy_0, g_yyzzzzz_xxxxyz_0, g_yyzzzzz_xxxxzz_0, g_yyzzzzz_xxxyyy_0, g_yyzzzzz_xxxyyz_0, g_yyzzzzz_xxxyzz_0, g_yyzzzzz_xxxzzz_0, g_yyzzzzz_xxyyyy_0, g_yyzzzzz_xxyyyz_0, g_yyzzzzz_xxyyzz_0, g_yyzzzzz_xxyzzz_0, g_yyzzzzz_xxzzzz_0, g_yyzzzzz_xyyyyy_0, g_yyzzzzz_xyyyyz_0, g_yyzzzzz_xyyyzz_0, g_yyzzzzz_xyyzzz_0, g_yyzzzzz_xyzzzz_0, g_yyzzzzz_xzzzzz_0, g_yyzzzzz_yyyyyy_0, g_yyzzzzz_yyyyyz_0, g_yyzzzzz_yyyyzz_0, g_yyzzzzz_yyyzzz_0, g_yyzzzzz_yyzzzz_0, g_yyzzzzz_yzzzzz_0, g_yyzzzzz_zzzzzz_0, g_yzzzzz_xxxxxx_1, g_yzzzzz_xxxxxz_1, g_yzzzzz_xxxxyz_1, g_yzzzzz_xxxxz_1, g_yzzzzz_xxxxzz_1, g_yzzzzz_xxxyyz_1, g_yzzzzz_xxxyz_1, g_yzzzzz_xxxyzz_1, g_yzzzzz_xxxzz_1, g_yzzzzz_xxxzzz_1, g_yzzzzz_xxyyyz_1, g_yzzzzz_xxyyz_1, g_yzzzzz_xxyyzz_1, g_yzzzzz_xxyzz_1, g_yzzzzz_xxyzzz_1, g_yzzzzz_xxzzz_1, g_yzzzzz_xxzzzz_1, g_yzzzzz_xyyyyz_1, g_yzzzzz_xyyyz_1, g_yzzzzz_xyyyzz_1, g_yzzzzz_xyyzz_1, g_yzzzzz_xyyzzz_1, g_yzzzzz_xyzzz_1, g_yzzzzz_xyzzzz_1, g_yzzzzz_xzzzz_1, g_yzzzzz_xzzzzz_1, g_yzzzzz_yyyyyz_1, g_yzzzzz_yyyyz_1, g_yzzzzz_yyyyzz_1, g_yzzzzz_yyyzz_1, g_yzzzzz_yyyzzz_1, g_yzzzzz_yyzzz_1, g_yzzzzz_yyzzzz_1, g_yzzzzz_yzzzz_1, g_yzzzzz_yzzzzz_1, g_yzzzzz_zzzzz_1, g_yzzzzz_zzzzzz_1, g_zzzzz_xxxxxx_0, g_zzzzz_xxxxxx_1, g_zzzzz_xxxxxz_0, g_zzzzz_xxxxxz_1, g_zzzzz_xxxxyz_0, g_zzzzz_xxxxyz_1, g_zzzzz_xxxxzz_0, g_zzzzz_xxxxzz_1, g_zzzzz_xxxyyz_0, g_zzzzz_xxxyyz_1, g_zzzzz_xxxyzz_0, g_zzzzz_xxxyzz_1, g_zzzzz_xxxzzz_0, g_zzzzz_xxxzzz_1, g_zzzzz_xxyyyz_0, g_zzzzz_xxyyyz_1, g_zzzzz_xxyyzz_0, g_zzzzz_xxyyzz_1, g_zzzzz_xxyzzz_0, g_zzzzz_xxyzzz_1, g_zzzzz_xxzzzz_0, g_zzzzz_xxzzzz_1, g_zzzzz_xyyyyz_0, g_zzzzz_xyyyyz_1, g_zzzzz_xyyyzz_0, g_zzzzz_xyyyzz_1, g_zzzzz_xyyzzz_0, g_zzzzz_xyyzzz_1, g_zzzzz_xyzzzz_0, g_zzzzz_xyzzzz_1, g_zzzzz_xzzzzz_0, g_zzzzz_xzzzzz_1, g_zzzzz_yyyyyz_0, g_zzzzz_yyyyyz_1, g_zzzzz_yyyyzz_0, g_zzzzz_yyyyzz_1, g_zzzzz_yyyzzz_0, g_zzzzz_yyyzzz_1, g_zzzzz_yyzzzz_0, g_zzzzz_yyzzzz_1, g_zzzzz_yzzzzz_0, g_zzzzz_yzzzzz_1, g_zzzzz_zzzzzz_0, g_zzzzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzzz_xxxxxx_0[i] = g_zzzzz_xxxxxx_0[i] * fbe_0 - g_zzzzz_xxxxxx_1[i] * fz_be_0 + g_yzzzzz_xxxxxx_1[i] * pa_y[i];

        g_yyzzzzz_xxxxxy_0[i] = 4.0 * g_yyzzz_xxxxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxxxxy_1[i] * fz_be_0 + g_yyzzzz_xxxxxy_1[i] * pa_z[i];

        g_yyzzzzz_xxxxxz_0[i] = g_zzzzz_xxxxxz_0[i] * fbe_0 - g_zzzzz_xxxxxz_1[i] * fz_be_0 + g_yzzzzz_xxxxxz_1[i] * pa_y[i];

        g_yyzzzzz_xxxxyy_0[i] = 4.0 * g_yyzzz_xxxxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxxxyy_1[i] * fz_be_0 + g_yyzzzz_xxxxyy_1[i] * pa_z[i];

        g_yyzzzzz_xxxxyz_0[i] = g_zzzzz_xxxxyz_0[i] * fbe_0 - g_zzzzz_xxxxyz_1[i] * fz_be_0 + g_yzzzzz_xxxxz_1[i] * fe_0 + g_yzzzzz_xxxxyz_1[i] * pa_y[i];

        g_yyzzzzz_xxxxzz_0[i] = g_zzzzz_xxxxzz_0[i] * fbe_0 - g_zzzzz_xxxxzz_1[i] * fz_be_0 + g_yzzzzz_xxxxzz_1[i] * pa_y[i];

        g_yyzzzzz_xxxyyy_0[i] = 4.0 * g_yyzzz_xxxyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxxyyy_1[i] * fz_be_0 + g_yyzzzz_xxxyyy_1[i] * pa_z[i];

        g_yyzzzzz_xxxyyz_0[i] = g_zzzzz_xxxyyz_0[i] * fbe_0 - g_zzzzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_xxxyz_1[i] * fe_0 + g_yzzzzz_xxxyyz_1[i] * pa_y[i];

        g_yyzzzzz_xxxyzz_0[i] = g_zzzzz_xxxyzz_0[i] * fbe_0 - g_zzzzz_xxxyzz_1[i] * fz_be_0 + g_yzzzzz_xxxzz_1[i] * fe_0 + g_yzzzzz_xxxyzz_1[i] * pa_y[i];

        g_yyzzzzz_xxxzzz_0[i] = g_zzzzz_xxxzzz_0[i] * fbe_0 - g_zzzzz_xxxzzz_1[i] * fz_be_0 + g_yzzzzz_xxxzzz_1[i] * pa_y[i];

        g_yyzzzzz_xxyyyy_0[i] = 4.0 * g_yyzzz_xxyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxyyyy_1[i] * fz_be_0 + g_yyzzzz_xxyyyy_1[i] * pa_z[i];

        g_yyzzzzz_xxyyyz_0[i] = g_zzzzz_xxyyyz_0[i] * fbe_0 - g_zzzzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_xxyyz_1[i] * fe_0 + g_yzzzzz_xxyyyz_1[i] * pa_y[i];

        g_yyzzzzz_xxyyzz_0[i] = g_zzzzz_xxyyzz_0[i] * fbe_0 - g_zzzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_xxyzz_1[i] * fe_0 + g_yzzzzz_xxyyzz_1[i] * pa_y[i];

        g_yyzzzzz_xxyzzz_0[i] = g_zzzzz_xxyzzz_0[i] * fbe_0 - g_zzzzz_xxyzzz_1[i] * fz_be_0 + g_yzzzzz_xxzzz_1[i] * fe_0 + g_yzzzzz_xxyzzz_1[i] * pa_y[i];

        g_yyzzzzz_xxzzzz_0[i] = g_zzzzz_xxzzzz_0[i] * fbe_0 - g_zzzzz_xxzzzz_1[i] * fz_be_0 + g_yzzzzz_xxzzzz_1[i] * pa_y[i];

        g_yyzzzzz_xyyyyy_0[i] = 4.0 * g_yyzzz_xyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xyyyyy_1[i] * fz_be_0 + g_yyzzzz_xyyyyy_1[i] * pa_z[i];

        g_yyzzzzz_xyyyyz_0[i] = g_zzzzz_xyyyyz_0[i] * fbe_0 - g_zzzzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_xyyyz_1[i] * fe_0 + g_yzzzzz_xyyyyz_1[i] * pa_y[i];

        g_yyzzzzz_xyyyzz_0[i] = g_zzzzz_xyyyzz_0[i] * fbe_0 - g_zzzzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_xyyzz_1[i] * fe_0 + g_yzzzzz_xyyyzz_1[i] * pa_y[i];

        g_yyzzzzz_xyyzzz_0[i] = g_zzzzz_xyyzzz_0[i] * fbe_0 - g_zzzzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_xyzzz_1[i] * fe_0 + g_yzzzzz_xyyzzz_1[i] * pa_y[i];

        g_yyzzzzz_xyzzzz_0[i] = g_zzzzz_xyzzzz_0[i] * fbe_0 - g_zzzzz_xyzzzz_1[i] * fz_be_0 + g_yzzzzz_xzzzz_1[i] * fe_0 + g_yzzzzz_xyzzzz_1[i] * pa_y[i];

        g_yyzzzzz_xzzzzz_0[i] = g_zzzzz_xzzzzz_0[i] * fbe_0 - g_zzzzz_xzzzzz_1[i] * fz_be_0 + g_yzzzzz_xzzzzz_1[i] * pa_y[i];

        g_yyzzzzz_yyyyyy_0[i] = 4.0 * g_yyzzz_yyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_yyyyyy_1[i] * fz_be_0 + g_yyzzzz_yyyyyy_1[i] * pa_z[i];

        g_yyzzzzz_yyyyyz_0[i] = g_zzzzz_yyyyyz_0[i] * fbe_0 - g_zzzzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_yyyyz_1[i] * fe_0 + g_yzzzzz_yyyyyz_1[i] * pa_y[i];

        g_yyzzzzz_yyyyzz_0[i] = g_zzzzz_yyyyzz_0[i] * fbe_0 - g_zzzzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_yyyzz_1[i] * fe_0 + g_yzzzzz_yyyyzz_1[i] * pa_y[i];

        g_yyzzzzz_yyyzzz_0[i] = g_zzzzz_yyyzzz_0[i] * fbe_0 - g_zzzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_yyzzz_1[i] * fe_0 + g_yzzzzz_yyyzzz_1[i] * pa_y[i];

        g_yyzzzzz_yyzzzz_0[i] = g_zzzzz_yyzzzz_0[i] * fbe_0 - g_zzzzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_yzzzz_1[i] * fe_0 + g_yzzzzz_yyzzzz_1[i] * pa_y[i];

        g_yyzzzzz_yzzzzz_0[i] = g_zzzzz_yzzzzz_0[i] * fbe_0 - g_zzzzz_yzzzzz_1[i] * fz_be_0 + g_yzzzzz_zzzzz_1[i] * fe_0 + g_yzzzzz_yzzzzz_1[i] * pa_y[i];

        g_yyzzzzz_zzzzzz_0[i] = g_zzzzz_zzzzzz_0[i] * fbe_0 - g_zzzzz_zzzzzz_1[i] * fz_be_0 + g_yzzzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 952-980 components of targeted buffer : KI

    auto g_yzzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 952);

    auto g_yzzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 953);

    auto g_yzzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 954);

    auto g_yzzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 955);

    auto g_yzzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 956);

    auto g_yzzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 957);

    auto g_yzzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 958);

    auto g_yzzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 959);

    auto g_yzzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 960);

    auto g_yzzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 961);

    auto g_yzzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 962);

    auto g_yzzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 963);

    auto g_yzzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 964);

    auto g_yzzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 965);

    auto g_yzzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 966);

    auto g_yzzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 967);

    auto g_yzzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 968);

    auto g_yzzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 969);

    auto g_yzzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 970);

    auto g_yzzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 971);

    auto g_yzzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 972);

    auto g_yzzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 973);

    auto g_yzzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 974);

    auto g_yzzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 975);

    auto g_yzzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 976);

    auto g_yzzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 977);

    auto g_yzzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 978);

    auto g_yzzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 979);

    #pragma omp simd aligned(g_yzzzzzz_xxxxxx_0, g_yzzzzzz_xxxxxy_0, g_yzzzzzz_xxxxxz_0, g_yzzzzzz_xxxxyy_0, g_yzzzzzz_xxxxyz_0, g_yzzzzzz_xxxxzz_0, g_yzzzzzz_xxxyyy_0, g_yzzzzzz_xxxyyz_0, g_yzzzzzz_xxxyzz_0, g_yzzzzzz_xxxzzz_0, g_yzzzzzz_xxyyyy_0, g_yzzzzzz_xxyyyz_0, g_yzzzzzz_xxyyzz_0, g_yzzzzzz_xxyzzz_0, g_yzzzzzz_xxzzzz_0, g_yzzzzzz_xyyyyy_0, g_yzzzzzz_xyyyyz_0, g_yzzzzzz_xyyyzz_0, g_yzzzzzz_xyyzzz_0, g_yzzzzzz_xyzzzz_0, g_yzzzzzz_xzzzzz_0, g_yzzzzzz_yyyyyy_0, g_yzzzzzz_yyyyyz_0, g_yzzzzzz_yyyyzz_0, g_yzzzzzz_yyyzzz_0, g_yzzzzzz_yyzzzz_0, g_yzzzzzz_yzzzzz_0, g_yzzzzzz_zzzzzz_0, g_zzzzzz_xxxxx_1, g_zzzzzz_xxxxxx_1, g_zzzzzz_xxxxxy_1, g_zzzzzz_xxxxxz_1, g_zzzzzz_xxxxy_1, g_zzzzzz_xxxxyy_1, g_zzzzzz_xxxxyz_1, g_zzzzzz_xxxxz_1, g_zzzzzz_xxxxzz_1, g_zzzzzz_xxxyy_1, g_zzzzzz_xxxyyy_1, g_zzzzzz_xxxyyz_1, g_zzzzzz_xxxyz_1, g_zzzzzz_xxxyzz_1, g_zzzzzz_xxxzz_1, g_zzzzzz_xxxzzz_1, g_zzzzzz_xxyyy_1, g_zzzzzz_xxyyyy_1, g_zzzzzz_xxyyyz_1, g_zzzzzz_xxyyz_1, g_zzzzzz_xxyyzz_1, g_zzzzzz_xxyzz_1, g_zzzzzz_xxyzzz_1, g_zzzzzz_xxzzz_1, g_zzzzzz_xxzzzz_1, g_zzzzzz_xyyyy_1, g_zzzzzz_xyyyyy_1, g_zzzzzz_xyyyyz_1, g_zzzzzz_xyyyz_1, g_zzzzzz_xyyyzz_1, g_zzzzzz_xyyzz_1, g_zzzzzz_xyyzzz_1, g_zzzzzz_xyzzz_1, g_zzzzzz_xyzzzz_1, g_zzzzzz_xzzzz_1, g_zzzzzz_xzzzzz_1, g_zzzzzz_yyyyy_1, g_zzzzzz_yyyyyy_1, g_zzzzzz_yyyyyz_1, g_zzzzzz_yyyyz_1, g_zzzzzz_yyyyzz_1, g_zzzzzz_yyyzz_1, g_zzzzzz_yyyzzz_1, g_zzzzzz_yyzzz_1, g_zzzzzz_yyzzzz_1, g_zzzzzz_yzzzz_1, g_zzzzzz_yzzzzz_1, g_zzzzzz_zzzzz_1, g_zzzzzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzzz_xxxxxx_0[i] = g_zzzzzz_xxxxxx_1[i] * pa_y[i];

        g_yzzzzzz_xxxxxy_0[i] = g_zzzzzz_xxxxx_1[i] * fe_0 + g_zzzzzz_xxxxxy_1[i] * pa_y[i];

        g_yzzzzzz_xxxxxz_0[i] = g_zzzzzz_xxxxxz_1[i] * pa_y[i];

        g_yzzzzzz_xxxxyy_0[i] = 2.0 * g_zzzzzz_xxxxy_1[i] * fe_0 + g_zzzzzz_xxxxyy_1[i] * pa_y[i];

        g_yzzzzzz_xxxxyz_0[i] = g_zzzzzz_xxxxz_1[i] * fe_0 + g_zzzzzz_xxxxyz_1[i] * pa_y[i];

        g_yzzzzzz_xxxxzz_0[i] = g_zzzzzz_xxxxzz_1[i] * pa_y[i];

        g_yzzzzzz_xxxyyy_0[i] = 3.0 * g_zzzzzz_xxxyy_1[i] * fe_0 + g_zzzzzz_xxxyyy_1[i] * pa_y[i];

        g_yzzzzzz_xxxyyz_0[i] = 2.0 * g_zzzzzz_xxxyz_1[i] * fe_0 + g_zzzzzz_xxxyyz_1[i] * pa_y[i];

        g_yzzzzzz_xxxyzz_0[i] = g_zzzzzz_xxxzz_1[i] * fe_0 + g_zzzzzz_xxxyzz_1[i] * pa_y[i];

        g_yzzzzzz_xxxzzz_0[i] = g_zzzzzz_xxxzzz_1[i] * pa_y[i];

        g_yzzzzzz_xxyyyy_0[i] = 4.0 * g_zzzzzz_xxyyy_1[i] * fe_0 + g_zzzzzz_xxyyyy_1[i] * pa_y[i];

        g_yzzzzzz_xxyyyz_0[i] = 3.0 * g_zzzzzz_xxyyz_1[i] * fe_0 + g_zzzzzz_xxyyyz_1[i] * pa_y[i];

        g_yzzzzzz_xxyyzz_0[i] = 2.0 * g_zzzzzz_xxyzz_1[i] * fe_0 + g_zzzzzz_xxyyzz_1[i] * pa_y[i];

        g_yzzzzzz_xxyzzz_0[i] = g_zzzzzz_xxzzz_1[i] * fe_0 + g_zzzzzz_xxyzzz_1[i] * pa_y[i];

        g_yzzzzzz_xxzzzz_0[i] = g_zzzzzz_xxzzzz_1[i] * pa_y[i];

        g_yzzzzzz_xyyyyy_0[i] = 5.0 * g_zzzzzz_xyyyy_1[i] * fe_0 + g_zzzzzz_xyyyyy_1[i] * pa_y[i];

        g_yzzzzzz_xyyyyz_0[i] = 4.0 * g_zzzzzz_xyyyz_1[i] * fe_0 + g_zzzzzz_xyyyyz_1[i] * pa_y[i];

        g_yzzzzzz_xyyyzz_0[i] = 3.0 * g_zzzzzz_xyyzz_1[i] * fe_0 + g_zzzzzz_xyyyzz_1[i] * pa_y[i];

        g_yzzzzzz_xyyzzz_0[i] = 2.0 * g_zzzzzz_xyzzz_1[i] * fe_0 + g_zzzzzz_xyyzzz_1[i] * pa_y[i];

        g_yzzzzzz_xyzzzz_0[i] = g_zzzzzz_xzzzz_1[i] * fe_0 + g_zzzzzz_xyzzzz_1[i] * pa_y[i];

        g_yzzzzzz_xzzzzz_0[i] = g_zzzzzz_xzzzzz_1[i] * pa_y[i];

        g_yzzzzzz_yyyyyy_0[i] = 6.0 * g_zzzzzz_yyyyy_1[i] * fe_0 + g_zzzzzz_yyyyyy_1[i] * pa_y[i];

        g_yzzzzzz_yyyyyz_0[i] = 5.0 * g_zzzzzz_yyyyz_1[i] * fe_0 + g_zzzzzz_yyyyyz_1[i] * pa_y[i];

        g_yzzzzzz_yyyyzz_0[i] = 4.0 * g_zzzzzz_yyyzz_1[i] * fe_0 + g_zzzzzz_yyyyzz_1[i] * pa_y[i];

        g_yzzzzzz_yyyzzz_0[i] = 3.0 * g_zzzzzz_yyzzz_1[i] * fe_0 + g_zzzzzz_yyyzzz_1[i] * pa_y[i];

        g_yzzzzzz_yyzzzz_0[i] = 2.0 * g_zzzzzz_yzzzz_1[i] * fe_0 + g_zzzzzz_yyzzzz_1[i] * pa_y[i];

        g_yzzzzzz_yzzzzz_0[i] = g_zzzzzz_zzzzz_1[i] * fe_0 + g_zzzzzz_yzzzzz_1[i] * pa_y[i];

        g_yzzzzzz_zzzzzz_0[i] = g_zzzzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 980-1008 components of targeted buffer : KI

    auto g_zzzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ki + 980);

    auto g_zzzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ki + 981);

    auto g_zzzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ki + 982);

    auto g_zzzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ki + 983);

    auto g_zzzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ki + 984);

    auto g_zzzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ki + 985);

    auto g_zzzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ki + 986);

    auto g_zzzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ki + 987);

    auto g_zzzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ki + 988);

    auto g_zzzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ki + 989);

    auto g_zzzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ki + 990);

    auto g_zzzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ki + 991);

    auto g_zzzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ki + 992);

    auto g_zzzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ki + 993);

    auto g_zzzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ki + 994);

    auto g_zzzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ki + 995);

    auto g_zzzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ki + 996);

    auto g_zzzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ki + 997);

    auto g_zzzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ki + 998);

    auto g_zzzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ki + 999);

    auto g_zzzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ki + 1000);

    auto g_zzzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ki + 1001);

    auto g_zzzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ki + 1002);

    auto g_zzzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ki + 1003);

    auto g_zzzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ki + 1004);

    auto g_zzzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ki + 1005);

    auto g_zzzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ki + 1006);

    auto g_zzzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ki + 1007);

    #pragma omp simd aligned(g_zzzzz_xxxxxx_0, g_zzzzz_xxxxxx_1, g_zzzzz_xxxxxy_0, g_zzzzz_xxxxxy_1, g_zzzzz_xxxxxz_0, g_zzzzz_xxxxxz_1, g_zzzzz_xxxxyy_0, g_zzzzz_xxxxyy_1, g_zzzzz_xxxxyz_0, g_zzzzz_xxxxyz_1, g_zzzzz_xxxxzz_0, g_zzzzz_xxxxzz_1, g_zzzzz_xxxyyy_0, g_zzzzz_xxxyyy_1, g_zzzzz_xxxyyz_0, g_zzzzz_xxxyyz_1, g_zzzzz_xxxyzz_0, g_zzzzz_xxxyzz_1, g_zzzzz_xxxzzz_0, g_zzzzz_xxxzzz_1, g_zzzzz_xxyyyy_0, g_zzzzz_xxyyyy_1, g_zzzzz_xxyyyz_0, g_zzzzz_xxyyyz_1, g_zzzzz_xxyyzz_0, g_zzzzz_xxyyzz_1, g_zzzzz_xxyzzz_0, g_zzzzz_xxyzzz_1, g_zzzzz_xxzzzz_0, g_zzzzz_xxzzzz_1, g_zzzzz_xyyyyy_0, g_zzzzz_xyyyyy_1, g_zzzzz_xyyyyz_0, g_zzzzz_xyyyyz_1, g_zzzzz_xyyyzz_0, g_zzzzz_xyyyzz_1, g_zzzzz_xyyzzz_0, g_zzzzz_xyyzzz_1, g_zzzzz_xyzzzz_0, g_zzzzz_xyzzzz_1, g_zzzzz_xzzzzz_0, g_zzzzz_xzzzzz_1, g_zzzzz_yyyyyy_0, g_zzzzz_yyyyyy_1, g_zzzzz_yyyyyz_0, g_zzzzz_yyyyyz_1, g_zzzzz_yyyyzz_0, g_zzzzz_yyyyzz_1, g_zzzzz_yyyzzz_0, g_zzzzz_yyyzzz_1, g_zzzzz_yyzzzz_0, g_zzzzz_yyzzzz_1, g_zzzzz_yzzzzz_0, g_zzzzz_yzzzzz_1, g_zzzzz_zzzzzz_0, g_zzzzz_zzzzzz_1, g_zzzzzz_xxxxx_1, g_zzzzzz_xxxxxx_1, g_zzzzzz_xxxxxy_1, g_zzzzzz_xxxxxz_1, g_zzzzzz_xxxxy_1, g_zzzzzz_xxxxyy_1, g_zzzzzz_xxxxyz_1, g_zzzzzz_xxxxz_1, g_zzzzzz_xxxxzz_1, g_zzzzzz_xxxyy_1, g_zzzzzz_xxxyyy_1, g_zzzzzz_xxxyyz_1, g_zzzzzz_xxxyz_1, g_zzzzzz_xxxyzz_1, g_zzzzzz_xxxzz_1, g_zzzzzz_xxxzzz_1, g_zzzzzz_xxyyy_1, g_zzzzzz_xxyyyy_1, g_zzzzzz_xxyyyz_1, g_zzzzzz_xxyyz_1, g_zzzzzz_xxyyzz_1, g_zzzzzz_xxyzz_1, g_zzzzzz_xxyzzz_1, g_zzzzzz_xxzzz_1, g_zzzzzz_xxzzzz_1, g_zzzzzz_xyyyy_1, g_zzzzzz_xyyyyy_1, g_zzzzzz_xyyyyz_1, g_zzzzzz_xyyyz_1, g_zzzzzz_xyyyzz_1, g_zzzzzz_xyyzz_1, g_zzzzzz_xyyzzz_1, g_zzzzzz_xyzzz_1, g_zzzzzz_xyzzzz_1, g_zzzzzz_xzzzz_1, g_zzzzzz_xzzzzz_1, g_zzzzzz_yyyyy_1, g_zzzzzz_yyyyyy_1, g_zzzzzz_yyyyyz_1, g_zzzzzz_yyyyz_1, g_zzzzzz_yyyyzz_1, g_zzzzzz_yyyzz_1, g_zzzzzz_yyyzzz_1, g_zzzzzz_yyzzz_1, g_zzzzzz_yyzzzz_1, g_zzzzzz_yzzzz_1, g_zzzzzz_yzzzzz_1, g_zzzzzz_zzzzz_1, g_zzzzzz_zzzzzz_1, g_zzzzzzz_xxxxxx_0, g_zzzzzzz_xxxxxy_0, g_zzzzzzz_xxxxxz_0, g_zzzzzzz_xxxxyy_0, g_zzzzzzz_xxxxyz_0, g_zzzzzzz_xxxxzz_0, g_zzzzzzz_xxxyyy_0, g_zzzzzzz_xxxyyz_0, g_zzzzzzz_xxxyzz_0, g_zzzzzzz_xxxzzz_0, g_zzzzzzz_xxyyyy_0, g_zzzzzzz_xxyyyz_0, g_zzzzzzz_xxyyzz_0, g_zzzzzzz_xxyzzz_0, g_zzzzzzz_xxzzzz_0, g_zzzzzzz_xyyyyy_0, g_zzzzzzz_xyyyyz_0, g_zzzzzzz_xyyyzz_0, g_zzzzzzz_xyyzzz_0, g_zzzzzzz_xyzzzz_0, g_zzzzzzz_xzzzzz_0, g_zzzzzzz_yyyyyy_0, g_zzzzzzz_yyyyyz_0, g_zzzzzzz_yyyyzz_0, g_zzzzzzz_yyyzzz_0, g_zzzzzzz_yyzzzz_0, g_zzzzzzz_yzzzzz_0, g_zzzzzzz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzzz_xxxxxx_0[i] = 6.0 * g_zzzzz_xxxxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxxx_1[i] * fz_be_0 + g_zzzzzz_xxxxxx_1[i] * pa_z[i];

        g_zzzzzzz_xxxxxy_0[i] = 6.0 * g_zzzzz_xxxxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxxy_1[i] * fz_be_0 + g_zzzzzz_xxxxxy_1[i] * pa_z[i];

        g_zzzzzzz_xxxxxz_0[i] = 6.0 * g_zzzzz_xxxxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxxz_1[i] * fz_be_0 + g_zzzzzz_xxxxx_1[i] * fe_0 + g_zzzzzz_xxxxxz_1[i] * pa_z[i];

        g_zzzzzzz_xxxxyy_0[i] = 6.0 * g_zzzzz_xxxxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxyy_1[i] * fz_be_0 + g_zzzzzz_xxxxyy_1[i] * pa_z[i];

        g_zzzzzzz_xxxxyz_0[i] = 6.0 * g_zzzzz_xxxxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxyz_1[i] * fz_be_0 + g_zzzzzz_xxxxy_1[i] * fe_0 + g_zzzzzz_xxxxyz_1[i] * pa_z[i];

        g_zzzzzzz_xxxxzz_0[i] = 6.0 * g_zzzzz_xxxxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xxxxz_1[i] * fe_0 + g_zzzzzz_xxxxzz_1[i] * pa_z[i];

        g_zzzzzzz_xxxyyy_0[i] = 6.0 * g_zzzzz_xxxyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxyyy_1[i] * fz_be_0 + g_zzzzzz_xxxyyy_1[i] * pa_z[i];

        g_zzzzzzz_xxxyyz_0[i] = 6.0 * g_zzzzz_xxxyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxyyz_1[i] * fz_be_0 + g_zzzzzz_xxxyy_1[i] * fe_0 + g_zzzzzz_xxxyyz_1[i] * pa_z[i];

        g_zzzzzzz_xxxyzz_0[i] = 6.0 * g_zzzzz_xxxyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xxxyz_1[i] * fe_0 + g_zzzzzz_xxxyzz_1[i] * pa_z[i];

        g_zzzzzzz_xxxzzz_0[i] = 6.0 * g_zzzzz_xxxzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_xxxzz_1[i] * fe_0 + g_zzzzzz_xxxzzz_1[i] * pa_z[i];

        g_zzzzzzz_xxyyyy_0[i] = 6.0 * g_zzzzz_xxyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyyyy_1[i] * fz_be_0 + g_zzzzzz_xxyyyy_1[i] * pa_z[i];

        g_zzzzzzz_xxyyyz_0[i] = 6.0 * g_zzzzz_xxyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyyyz_1[i] * fz_be_0 + g_zzzzzz_xxyyy_1[i] * fe_0 + g_zzzzzz_xxyyyz_1[i] * pa_z[i];

        g_zzzzzzz_xxyyzz_0[i] = 6.0 * g_zzzzz_xxyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xxyyz_1[i] * fe_0 + g_zzzzzz_xxyyzz_1[i] * pa_z[i];

        g_zzzzzzz_xxyzzz_0[i] = 6.0 * g_zzzzz_xxyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_xxyzz_1[i] * fe_0 + g_zzzzzz_xxyzzz_1[i] * pa_z[i];

        g_zzzzzzz_xxzzzz_0[i] = 6.0 * g_zzzzz_xxzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_xxzzz_1[i] * fe_0 + g_zzzzzz_xxzzzz_1[i] * pa_z[i];

        g_zzzzzzz_xyyyyy_0[i] = 6.0 * g_zzzzz_xyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyyyy_1[i] * fz_be_0 + g_zzzzzz_xyyyyy_1[i] * pa_z[i];

        g_zzzzzzz_xyyyyz_0[i] = 6.0 * g_zzzzz_xyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyyyz_1[i] * fz_be_0 + g_zzzzzz_xyyyy_1[i] * fe_0 + g_zzzzzz_xyyyyz_1[i] * pa_z[i];

        g_zzzzzzz_xyyyzz_0[i] = 6.0 * g_zzzzz_xyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xyyyz_1[i] * fe_0 + g_zzzzzz_xyyyzz_1[i] * pa_z[i];

        g_zzzzzzz_xyyzzz_0[i] = 6.0 * g_zzzzz_xyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_xyyzz_1[i] * fe_0 + g_zzzzzz_xyyzzz_1[i] * pa_z[i];

        g_zzzzzzz_xyzzzz_0[i] = 6.0 * g_zzzzz_xyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_xyzzz_1[i] * fe_0 + g_zzzzzz_xyzzzz_1[i] * pa_z[i];

        g_zzzzzzz_xzzzzz_0[i] = 6.0 * g_zzzzz_xzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_xzzzz_1[i] * fe_0 + g_zzzzzz_xzzzzz_1[i] * pa_z[i];

        g_zzzzzzz_yyyyyy_0[i] = 6.0 * g_zzzzz_yyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyyyy_1[i] * fz_be_0 + g_zzzzzz_yyyyyy_1[i] * pa_z[i];

        g_zzzzzzz_yyyyyz_0[i] = 6.0 * g_zzzzz_yyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyyyz_1[i] * fz_be_0 + g_zzzzzz_yyyyy_1[i] * fe_0 + g_zzzzzz_yyyyyz_1[i] * pa_z[i];

        g_zzzzzzz_yyyyzz_0[i] = 6.0 * g_zzzzz_yyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_yyyyz_1[i] * fe_0 + g_zzzzzz_yyyyzz_1[i] * pa_z[i];

        g_zzzzzzz_yyyzzz_0[i] = 6.0 * g_zzzzz_yyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_yyyzz_1[i] * fe_0 + g_zzzzzz_yyyzzz_1[i] * pa_z[i];

        g_zzzzzzz_yyzzzz_0[i] = 6.0 * g_zzzzz_yyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_yyzzz_1[i] * fe_0 + g_zzzzzz_yyzzzz_1[i] * pa_z[i];

        g_zzzzzzz_yzzzzz_0[i] = 6.0 * g_zzzzz_yzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_yzzzz_1[i] * fe_0 + g_zzzzzz_yzzzzz_1[i] * pa_z[i];

        g_zzzzzzz_zzzzzz_0[i] = 6.0 * g_zzzzz_zzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_zzzzz_1[i] * fe_0 + g_zzzzzz_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

